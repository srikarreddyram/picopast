"""The 1991-2005 ship counts for the climatological pattern check (pre-registration, check C).

BCO-DMO's global compilation already contains the 1995 JGOFS Arabian Sea
counts that PANGAEA holds station by station (same values; positions
rounded differently, and depths sometimes very differently: PANGAEA gives
TT054 depths in whole metres, 0-3 m, where BCO-DMO has 2.5-2.8 m). Counting
both would double the sample and invent cruises, so the match is made per
station, not per sample: a BCO-DMO record is dropped when PANGAEA has any
record from the same date within ``tol_deg`` in latitude and longitude.
PANGAEA's profile is kept because it carries the real cruise label, which
is the bootstrap unit.
"""

from __future__ import annotations

import pandas as pd

from picopast.insitu.schema import standardise


def drop_duplicates_of(bcodmo: pd.DataFrame, pangaea: pd.DataFrame, tol_deg: float = 0.05) -> pd.DataFrame:
    """BCO-DMO records from stations PANGAEA does not hold (same date, within ``tol_deg``)."""
    b = bcodmo.reset_index(drop=True).assign(_row=lambda d: d.index, date=lambda d: d["time_utc"].dt.normalize())
    p = pangaea.assign(date=pangaea["time_utc"].dt.normalize())[["date", "lat", "lon"]].drop_duplicates()
    pairs = b.merge(p, on="date", suffixes=("", "_p"))
    twin = ((pairs["lat"] - pairs["lat_p"]).abs() < tol_deg) & ((pairs["lon"] - pairs["lon_p"]).abs() < tol_deg)
    dup_rows = set(pairs.loc[twin, "_row"])
    return b[~b["_row"].isin(dup_rows)].drop(columns=["_row", "date"])


def historical_counts(bcodmo: pd.DataFrame, pangaea: pd.DataFrame) -> pd.DataFrame:
    """Both archives in the common layout, de-duplicated, PANGAEA preferred."""
    b, p = standardise(bcodmo), standardise(pangaea)
    return standardise(pd.concat([p, drop_duplicates_of(b, p)], ignore_index=True))
