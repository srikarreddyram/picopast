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


def load_hplc(paths: list[Path]) -> pd.DataFrame:
    """All samples from SeaBASS HPLC files, one row per sample, pigments in mg m-3."""
    frames = []
    for path in paths:
        header, df = read_sb(path)
        time = pd.to_datetime(
            df["date"].astype("Int64").astype(str) + " " + df.get("time", "00:00:00").astype(str),
            format="%Y%m%d %H:%M:%S", errors="coerce",
        )
        frames.append(pd.DataFrame({
            "source": "seabass",
            "dataset_id": Path(path).name,
            "cruise": header.get("cruise", Path(path).stem),
            "station": df.get("station", pd.Series(index=df.index, dtype=object)).astype(str),
            "sample_id": df.get("hplc_gsfc_id", pd.Series(index=df.index, dtype=object)).astype(str),
            "time_utc": time,
            "lat": pd.to_numeric(df["lat"], errors="coerce"),
            "lon": pd.to_numeric(df["lon"], errors="coerce"),
            "depth_m": pd.to_numeric(df["depth"], errors="coerce"),
            **{p: pd.to_numeric(df[p], errors="coerce") if p in df else np.nan for p in PIGMENTS},
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
