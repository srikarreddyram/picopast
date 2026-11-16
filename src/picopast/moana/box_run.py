"""Run MOANA over the downloaded study-box PACE files and store one map per file.

Inputs are the box subsets written by ``picopast fetch``:
``data/raw/pace_rrs/PACE_OCI.<date>.L3m.<period>.AOP.*.box.nc`` and the
daily CMC SST analyses in ``data/raw/cmc_sst/``. A daily scene is paired
with the CMC file of the same date (it agrees better with NASA's own L4
output than the previous day's file). An 8-day composite is paired with the
mean of the daily CMC files in its period, and at least half of them must be
present.

Output: ``data/interim/moana/<period>/moana_<start>.nc`` with pro, syn and
peuk in cells/mL, MOANA flags, and a length-1 time axis at the period start,
so a run of files opens as one cube with ``xr.open_mfdataset``.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.moana.grid_run import moana_on_grid

_DATES = re.compile(r"PACE_OCI\.(\d{8})(?:_(\d{8}))?\.L3m\.(DAY|8D)\.")


def scene_period(path: Path) -> tuple[pd.Timestamp, pd.Timestamp, str]:
    """(first day, last day, period) from a PACE L3M file name."""
    m = _DATES.search(path.name)
    if not m:
        raise ValueError(f"not a PACE L3M file name: {path.name}")
    start = pd.Timestamp(m.group(1))
    end = pd.Timestamp(m.group(2)) if m.group(2) else start
    return start, end, m.group(3)


def _cmc_files(cfg: Config) -> dict[pd.Timestamp, Path]:
    out = {}
    for p in (cfg.path("raw") / "cmc_sst").glob("*CMC*.box.nc"):
        out[pd.Timestamp(p.name[:8])] = p
    return out


def sst_for_period(cfg: Config, start: pd.Timestamp, end: pd.Timestamp,
                   cmc: dict[pd.Timestamp, Path] | None = None) -> xr.DataArray | None:
    """CMC SST in deg C for the period: the day itself, or the mean over an 8-day period."""
    cmc = _cmc_files(cfg) if cmc is None else cmc
    days = pd.date_range(start, end)
    paths = [cmc[d] for d in days if d in cmc]
    if len(paths) < max(1, len(days) / 2):
        return None
    fields = []
    for p in paths:
        with xr.open_dataset(p) as ds:
            fields.append(ds["analysed_sst"].isel(time=0).load() - 273.15)
    return xr.concat(fields, dim="day").mean("day") if len(fields) > 1 else fields[0]


def run_file(cfg: Config, rrs_path: Path, pc_mapping: str = "operational",
             refresh: bool = False, cmc: dict | None = None) -> Path | None:
    """MOANA for one PACE box file. Returns the output path, or None if SST is missing."""
    start, end, period = scene_period(rrs_path)
    out_dir = cfg.path("interim") / "moana" / period
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"moana_{start:%Y%m%d}.nc"
    if out.exists() and not refresh:
        return out
    sst = sst_for_period(cfg, start, end, cmc)
    if sst is None:
        return None
    with xr.open_dataset(rrs_path) as rrs_ds:
        ds = moana_on_grid(rrs_ds["Rrs"], sst, pc_mapping=pc_mapping)
    ds = ds.expand_dims(time=[np.datetime64(start, "ns")])
    ds.attrs.update(period=period, period_start=f"{start:%Y-%m-%d}", period_end=f"{end:%Y-%m-%d}",
                    rrs_file=rrs_path.name, sst="GHRSST CMC 0.1 deg (same day; mean for 8D)")
    encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    tmp = out.with_suffix(".part")
    ds.to_netcdf(tmp, encoding=encoding)
    tmp.rename(out)
    return out


def run(cfg: Config, start: str, end: str, period: str = "DAY", pc_mapping: str = "operational",
        refresh: bool = False) -> tuple[list[Path], list[str]]:
    """MOANA for every downloaded box file of ``period`` starting within [start, end].

    Returns (written paths, names of scenes skipped for missing SST).
    """
    lo, hi = pd.Timestamp(start), pd.Timestamp(end)
    cmc = _cmc_files(cfg)
    written, skipped = [], []
    for p in sorted((cfg.path("raw") / "pace_rrs").glob(f"PACE_OCI.*.L3m.{period}.AOP.*.box.nc")):
        s, _, _ = scene_period(p)
        if not lo <= s <= hi:
            continue
        out = run_file(cfg, p, pc_mapping=pc_mapping, refresh=refresh, cmc=cmc)
        (written.append(out) if out else skipped.append(p.name))
    return written, skipped
