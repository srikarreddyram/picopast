"""HPLC pigments from SeaBASS: the ship side of the pigment-proxy check (pre-registration, check A).

Divinyl chlorophyll a (``dv_chl_a``) is made only by Prochlorococcus, so it
tracks where Prochlorococcus is more or less abundant. It is not a cell count:
cells per unit pigment change with light and nutrients. It can therefore test
ranks (Spearman), never absolute errors.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from picopast.config import Box
from picopast.insitu.seabass import read_sb

PIGMENTS = ("dv_chl_a", "zea", "tchl")


#: Field-name spellings seen in SeaBASS HPLC files (matched case-insensitively) -> our column.
_ALIASES = {"dv_chl_a": "dv_chl_a", "zea": "zea", "tchl": "tchl", "tot_chl_a": "tchl"}


def _col(df: pd.DataFrame, name: str, fallback: float = np.nan) -> pd.Series:
    return pd.to_numeric(df[name], errors="coerce") if name in df else pd.Series(fallback, index=df.index)


def _header_float(header: dict[str, str], key: str) -> float:
    try:
        return float(header.get(key, "").split("[")[0])
    except ValueError:
        return np.nan


def load_hplc(paths: list[Path]) -> pd.DataFrame:
    """All samples from SeaBASS HPLC files, one row per sample, pigments in mg m-3.

    Handles the three layouts met so far: a ``date`` column; ``year``/``month``/``day``
    columns; and older one-station files whose position, date, time and depth are only
    in the header. Field names are matched case-insensitively.
    """
    frames = []
    for path in paths:
        header, df = read_sb(path)
        df.columns = [c.lower() for c in df.columns]
        n = len(df)
        if "date" in df:
            day = df["date"].astype("Int64").astype(str)
        elif {"year", "month", "day"} <= set(df.columns):
            day = (df["year"].astype("Int64").astype(str) + df["month"].astype("Int64").astype(str).str.zfill(2)
                   + df["day"].astype("Int64").astype(str).str.zfill(2))
        else:
            day = pd.Series([header.get("start_date", "")] * n, index=df.index)
        clock = df["time"].astype(str) if "time" in df else pd.Series(
            [header.get("start_time", "00:00:00").split("[")[0]] * n, index=df.index)
        time = pd.to_datetime(day + " " + clock, format="%Y%m%d %H:%M:%S", errors="coerce")

        lat_h = (_header_float(header, "north_latitude") if _header_float(header, "north_latitude")
                 == _header_float(header, "south_latitude") else np.nan)
        lon_h = (_header_float(header, "east_longitude") if _header_float(header, "east_longitude")
                 == _header_float(header, "west_longitude") else np.nan)
        own = Path(path).name + "#" + pd.Series(range(n), index=df.index).astype(str)
        gsfc = df["hplc_gsfc_id"].astype(str).str.strip() if "hplc_gsfc_id" in df else own
        sample_id = gsfc.where(~gsfc.isin(["nan", "", "-9999", "<NA>"]), own)  # a missing id must not merge rows
        pig = {}
        for field, ours in _ALIASES.items():
            if field in df and ours not in pig:
                pig[ours] = pd.to_numeric(df[field], errors="coerce")
        frames.append(pd.DataFrame({
            "source": "seabass",
            "dataset_id": Path(path).name,
            "cruise": header.get("cruise", Path(path).stem),
            "station": (df["station"] if "station" in df else pd.Series(header.get("station", ""), index=df.index)
                        ).astype(str),
            "sample_id": sample_id,
            "time_utc": time,
            "lat": _col(df, "lat", lat_h),
            "lon": _col(df, "lon", lon_h),
            "depth_m": _col(df, "depth", _header_float(header, "measurement_depth")),
            **{p: pig.get(p, pd.Series(np.nan, index=df.index)) for p in PIGMENTS},
        }))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def station_days(samples: pd.DataFrame, box: Box, start: str, max_depth_m: float = 10.0,
                 resolution_deg: float = 1 / 24) -> pd.DataFrame:
    """Surface samples in the box from ``start`` on, merged to one row per station-day.

    A station-day is one UTC date and one pixel of the 4.6 km L3M grid: that is
    the unit the satellite sees, so co-located samples (replicates, or the same
    station in two file releases) are one observation. Pigments are the median
    of the merged samples. Samples sharing an ``hplc_gsfc_id`` are counted once.
    """
    s = samples.drop_duplicates("sample_id") if "sample_id" in samples else samples
    s = s[
        (s["depth_m"] <= max_depth_m)
        & (s["time_utc"] >= pd.Timestamp(start))
        & s["lat"].between(box.lat_min, box.lat_max)
        & s["lon"].between(box.lon_min, box.lon_max)
        & s["dv_chl_a"].notna()
    ].copy()
    s["date"] = s["time_utc"].dt.normalize()
    s["iy"] = np.floor((90.0 - s["lat"]) / resolution_deg).astype(int)
    s["ix"] = np.floor((s["lon"] + 180.0) / resolution_deg).astype(int)
    agg = {p: "median" for p in PIGMENTS if p in s}
    agg.update(lat="mean", lon="mean", time_utc="min", cruise="first", depth_m="min",
               station=lambda x: ",".join(sorted(set(x))), sample_id="size")
    out = s.groupby(["date", "iy", "ix"], as_index=False).agg(agg)
    out = out.rename(columns={"sample_id": "n_samples"})
    out["station_day"] = out["date"].dt.strftime("%Y%m%d") + "_" + out["iy"].astype(str) + "_" + out["ix"].astype(str)
    return out.sort_values("time_utc").reset_index(drop=True)
