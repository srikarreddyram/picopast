"""Open sources that need no login: NOAA OISST v2.1, the Dipole Mode Index and ONI."""

from __future__ import annotations

import struct
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

from picopast.config import Config

_DAYS_PER_REQUEST = 31
_PAUSE_S = 3.0


class RateLimited(RuntimeError):
    """The server answered HTTP 429; retrying at once only makes it worse."""


_SEASON_CENTRE_MONTH = {
    "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6,
    "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
}


def _dap_array(url: str, expr: str, dtype: str) -> np.ndarray:
    """Read one array through OPeNDAP's binary (DAP2 .dods) response.

    Plain HTTPS instead of the netCDF library's DAP client, which fails against
    this server on macOS wheels ("NetCDF: I/O failure"). The payload after the
    "Data:" marker is the element count (twice, big-endian int32) followed by
    big-endian values.
    """
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            r = requests.get(f"{url}.dods?{expr}", timeout=(15, 90))
            if r.status_code == 429:
                raise RateLimited(
                    "NOAA PSL is rate-limiting OPeNDAP requests from this address. Wait an hour "
                    "or use `picopast fetch-open --oisst-method file`."
                )
            r.raise_for_status()
            _, marker, data = r.content.partition(b"\nData:\n")
            if not marker:
                raise RuntimeError(f"unexpected OPeNDAP response for {expr}")
            n = struct.unpack(">i", data[:4])[0]
            values = np.frombuffer(data, dtype=dtype, count=n, offset=8)
            time.sleep(_PAUSE_S)  # stay well under the server's request limit
            return values.astype(values.dtype.newbyteorder("="))
        except RateLimited:
            raise
        except (requests.RequestException, RuntimeError, ValueError) as e:
            last_error = e  # truncated or dropped responses are common on this server
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"OPeNDAP read failed after retries: {expr}") from last_error


def _oisst_opendap(cfg: Config, year: int) -> xr.DataArray:
    """Read only the box through OPeNDAP, about a month per request."""
    box = cfg.box
    url = cfg["open_sources"]["oisst_opendap"].format(year=year)
    lat = _dap_array(url, "lat", ">f4")
    lon = _dap_array(url, "lon", ">f4")  # 0-360; the box (50-100E) needs no wrapping
    days = _dap_array(url, "time", ">f8")  # days since 1800-01-01
    iy = np.flatnonzero((lat >= box.lat_min) & (lat <= box.lat_max))
    ix = np.flatnonzero((lon >= box.lon_min) & (lon <= box.lon_max))
    parts = []
    for t0 in range(0, days.size, _DAYS_PER_REQUEST):
        t1 = min(t0 + _DAYS_PER_REQUEST, days.size) - 1
        expr = f"sst.sst[{t0}:{t1}][{iy[0]}:{iy[-1]}][{ix[0]}:{ix[-1]}]"
        parts.append(_dap_array(url, expr, ">f4").reshape(t1 - t0 + 1, iy.size, ix.size))
    sst = np.concatenate(parts)
    sst = np.where(np.abs(sst) > 100.0, np.nan, sst).astype(np.float32)  # fill value is -9.97e36
    return xr.DataArray(
        sst,
        dims=("time", "lat", "lon"),
        coords={
            "time": pd.Timestamp("1800-01-01") + pd.to_timedelta(days, unit="D"),
            "lat": lat[iy].astype(np.float64),
            "lon": lon[ix].astype(np.float64),
        },
        attrs={"source_url": url},
    )


def _download_resumable(url: str, local: Path, attempts: int = 6) -> None:
    """Stream ``url`` to ``local``, resuming with an HTTP Range request after a dropped connection.

    NOAA PSL often cuts long transfers partway through a 0.5 GB file. If the
    server ignores the Range header (answers 200, not 206), the file restarts.
    """
    local.unlink(missing_ok=True)
    last_error: Exception | None = None
    for attempt in range(attempts):
        have = local.stat().st_size if local.exists() else 0
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with requests.get(url, stream=True, timeout=(15, 120), headers=headers) as r:
                if r.status_code == 416:  # requested range starts at the end: already complete
                    return
                r.raise_for_status()
                mode = "ab" if have and r.status_code == 206 else "wb"
                with open(local, mode) as f:
                    for chunk in r.iter_content(chunk_size=2**22):
                        f.write(chunk)
            return
        except requests.RequestException as e:
            last_error = e
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"download failed after {attempts} attempts: {url}") from last_error


def _oisst_file(cfg: Config, year: int) -> xr.DataArray:
    """Download the whole global yearly file (about 0.5 GB), cut the box, delete the file."""
    box = cfg.box
    url = cfg["open_sources"]["oisst_file"].format(year=year)
    cache = cfg.path("data") / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    local = cache / f"oisst_global_{year}.nc"
    try:
        _download_resumable(url, local)
        with xr.open_dataset(local) as ds:
            sub = ds["sst"].sel(
                lat=slice(box.lat_min, box.lat_max), lon=slice(box.lon_min, box.lon_max)
            ).load()
    finally:
        local.unlink(missing_ok=True)
    return xr.DataArray(
        sub.values.astype(np.float32),
        dims=("time", "lat", "lon"),
        coords={"time": sub["time"].values, "lat": sub["lat"].values.astype(np.float64),
                "lon": sub["lon"].values.astype(np.float64)},
        attrs={"source_url": url},
    )


def fetch_oisst_year(
    cfg: Config, year: int, refresh: bool = False, method: str = "opendap"
) -> Path:
    """One year of daily OISST v2.1 (0.25 degree) for the study box.

    OISST is the SST predictor for 2012-2026 and the basis of marine heatwave
    detection (1991-2020 baseline). ``method="opendap"`` transfers only the box
    (about 30 MB a year) but NOAA limits how many requests it accepts;
    ``method="file"`` downloads each global yearly file (about 0.5 GB) and
    keeps only the box. The current year is always re-fetched because NOAA
    appends to it.
    """
    out_dir = cfg.path("raw") / "oisst"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"oisst_v2.1_box_{year}.nc"
    if out.exists() and not refresh and year < date.today().year:
        return out
    da = {"opendap": _oisst_opendap, "file": _oisst_file}[method](cfg, year)
    da.attrs.update(units="degC", long_name="Daily sea surface temperature",
                    source="NOAA OISST v2.1 (high resolution)")
    tmp = out.with_suffix(".part")
    da.to_dataset(name="sst").to_netcdf(tmp, encoding={"sst": {"zlib": True, "complevel": 4}})
    tmp.rename(out)
    return out


def parse_dmi(text: str) -> pd.DataFrame:
    """NOAA PSL monthly table (year + 12 values, -9999 missing) -> tidy frame."""
    rows = []
    for line in text.splitlines()[1:]:
        parts = line.split()
        if len(parts) != 13:
            break  # footer starts
        year = int(parts[0])
        for month, value in enumerate(parts[1:], start=1):
            v = float(value)
            rows.append((pd.Timestamp(year, month, 1), np.nan if v <= -999 else v))
    return pd.DataFrame(rows, columns=["date", "dmi"]).dropna().reset_index(drop=True)


def parse_oni(text: str) -> pd.DataFrame:
    """NOAA CPC ONI table (3-month seasons) -> tidy frame dated at the centre month."""
    rows = []
    for line in text.splitlines()[1:]:
        parts = line.split()
        if len(parts) != 4:
            continue
        season, year, _total, anom = parts
        rows.append((pd.Timestamp(int(year), _SEASON_CENTRE_MONTH[season], 1), float(anom)))
    return pd.DataFrame(rows, columns=["date", "oni"])


def fetch_indices(cfg: Config) -> dict[str, Path]:
    """Download DMI (IOD) and ONI (ENSO) and store raw text plus tidy CSV."""
    out_dir = cfg.path("raw") / "indices"
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, parser in (("dmi", parse_dmi), ("oni", parse_oni)):
        r = requests.get(cfg["open_sources"][name], timeout=60)
        r.raise_for_status()
        (out_dir / f"{name}_raw.txt").write_text(r.text)
        paths[name] = out_dir / f"{name}.csv"
        parser(r.text).to_csv(paths[name], index=False)
    return paths
