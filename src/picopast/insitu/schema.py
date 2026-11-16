"""One table layout for every ship count, whatever archive it came from."""

from __future__ import annotations

import pandas as pd

COLUMNS = [
    "source",        # archive: bcodmo, pangaea, seabass, published_table
    "dataset_id",    # archive identifier or DOI
    "cruise",        # cruise or expedition label; the bootstrap resamples these
    "station",
    "time_utc",      # datetime64; date-only records sit at 00:00
    "lat",
    "lon",
    "depth_m",       # NaN when the archive gives none
    "pro_cells_ml",
    "syn_cells_ml",
    "peuk_cells_ml",
    "method",        # e.g. flow cytometry
    "url",           # where this exact record can be downloaded again
]
GROUP_COLUMNS = {"pro": "pro_cells_ml", "syn": "syn_cells_ml", "peuk": "peuk_cells_ml"}


def standardise(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` with exactly the standard columns, typed, sorted and de-duplicated."""
    out = df.reindex(columns=COLUMNS).copy()
    if not pd.api.types.is_datetime64_any_dtype(out["time_utc"]):
        # ISO8601 accepts date-only and date-time strings in one column.
        out["time_utc"] = pd.to_datetime(out["time_utc"], format="ISO8601", utc=True, errors="coerce")
    if out["time_utc"].dt.tz is not None:
        out["time_utc"] = out["time_utc"].dt.tz_convert("UTC").dt.tz_localize(None)
    for col in ("lat", "lon", "depth_m", *GROUP_COLUMNS.values()):
        out[col] = pd.to_numeric(out[col], errors="coerce")
    for col in GROUP_COLUMNS.values():
        out.loc[out[col] < 0, col] = float("nan")  # archive fill values
    has_count = out[list(GROUP_COLUMNS.values())].notna().any(axis=1)
    out = out[has_count & out["lat"].notna() & out["lon"].notna()]
    out = out.drop_duplicates().sort_values(["time_utc", "lat", "lon", "depth_m"])
    return out.reset_index(drop=True)


def era(time_utc: pd.Series, pace_start: str = "2024-03-05") -> pd.Series:
    """Which research question a record can serve.

    pace_era -> RQ1 trust gate; viirs_era (2012 to PACE launch) -> RQ4 check of
    the reconstruction; pre_viirs -> climatology and niche-model baseline only.
    """
    t = pd.to_datetime(time_utc)
    out = pd.Series("pre_viirs", index=t.index, dtype="object")
    out[t >= pd.Timestamp("2012-01-02")] = "viirs_era"
    out[t >= pd.Timestamp(pace_start)] = "pace_era"
    out[t.isna()] = "undated"
    return out
