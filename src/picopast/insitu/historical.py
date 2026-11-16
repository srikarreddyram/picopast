"""The 1991-2005 ship counts for the climatological pattern check (pre-registration, check C).

BCO-DMO's global compilation already contains the 1995 JGOFS Arabian Sea
counts that PANGAEA holds station by station (same values; positions and
depths rounded differently). Counting both would double the sample and
invent cruises, so a BCO-DMO record is dropped when PANGAEA has a record
from the same date within ``tol_deg`` in latitude and longitude and
``tol_depth_m`` in depth. The PANGAEA copy is kept because it carries the
real cruise label, which is the bootstrap unit.
"""

from __future__ import annotations

import pandas as pd

from picopast.insitu.schema import standardise


def drop_duplicates_of(bcodmo: pd.DataFrame, pangaea: pd.DataFrame, tol_deg: float = 0.05,
                       tol_depth_m: float = 0.5) -> pd.DataFrame:
    """BCO-DMO records with no PANGAEA twin (same date, within the position and depth tolerances)."""
    b = bcodmo.reset_index(drop=True).assign(_row=lambda d: d.index, date=lambda d: d["time_utc"].dt.normalize())
    p = pangaea.assign(date=pangaea["time_utc"].dt.normalize())[["date", "lat", "lon", "depth_m"]]
    pairs = b.merge(p, on="date", suffixes=("", "_p"))
    twin = (
        ((pairs["lat"] - pairs["lat_p"]).abs() < tol_deg)
        & ((pairs["lon"] - pairs["lon_p"]).abs() < tol_deg)
        & ((pairs["depth_m"] - pairs["depth_m_p"]).abs() < tol_depth_m)
    )
    dup_rows = set(pairs.loc[twin, "_row"])
    return b[~b["_row"].isin(dup_rows)].drop(columns=["_row", "date"])


def historical_counts(bcodmo: pd.DataFrame, pangaea: pd.DataFrame) -> pd.DataFrame:
    """Both archives in the common layout, de-duplicated, PANGAEA preferred."""
    b, p = standardise(bcodmo), standardise(pangaea)
    return standardise(pd.concat([p, drop_duplicates_of(b, p)], ignore_index=True))
