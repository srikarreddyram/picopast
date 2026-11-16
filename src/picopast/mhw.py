"""Marine heatwaves by the Hobday et al. (2016) definition (PRD section 6 step 8).

A marine heatwave is a period of at least 5 days with SST above the seasonally
varying 90th percentile, where gaps of 2 days or fewer between such periods
are bridged. The climatology and threshold come from a fixed baseline
(1991-2020, OISST): for each day of year, all baseline values within an 11-day
window centred on that day are pooled, their mean and 90th percentile taken,
and both are smoothed with a 31-day running mean.

All functions work on arrays with time as the first axis, so one call covers
the whole grid.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def day_of_year_366(time) -> np.ndarray:
    """Day of year on a leap-year calendar (1-366), so 1 March is day 61 every year."""
    t = pd.DatetimeIndex(pd.to_datetime(np.asarray(time)))
    doy = t.day_of_year.to_numpy().copy()
    doy[(~t.is_leap_year) & (t.month > 2)] += 1
    return doy


def _circular_running_mean(x: np.ndarray, width: int) -> np.ndarray:
    half = width // 2
    padded = np.concatenate([x[-half:], x, x[:half]], axis=0)
    csum = np.cumsum(np.insert(padded, 0, 0.0, axis=0), axis=0)
    return (csum[width:] - csum[:-width]) / width


def climatology(
    sst: np.ndarray,
    time,
    baseline: tuple[str, str] = ("1991-01-01", "2020-12-31"),
    window_half_days: int = 5,
    smooth_days: int = 31,
    percentile: float = 90.0,
) -> tuple[np.ndarray, np.ndarray]:
    """(seasonal mean, threshold), each of shape (366, ...) indexed by day_of_year_366 - 1.

    29 February is filled by averaging 28 February and 1 March before smoothing.
    """
    t = pd.DatetimeIndex(pd.to_datetime(np.asarray(time)))
    in_base = np.asarray((t >= baseline[0]) & (t <= baseline[1]))
    if in_base.sum() < 365 * 10:
        raise ValueError("baseline holds fewer than 10 years of daily data")
    base = np.asarray(sst)[in_base]
    doy = day_of_year_366(t[in_base])
    shape = (366,) + base.shape[1:]
    seas, thresh = np.full(shape, np.nan), np.full(shape, np.nan)
    offsets = np.arange(-window_half_days, window_half_days + 1)
    for d in range(1, 367):
        if d == 60:
            continue  # 29 February, filled below
        wanted = (d + offsets - 1) % 366 + 1
        sample = base[np.isin(doy, wanted)]
        with np.errstate(all="ignore"):
            seas[d - 1] = np.nanmean(sample, axis=0)
            thresh[d - 1] = np.nanpercentile(sample, percentile, axis=0)
    for arr in (seas, thresh):
        arr[59] = 0.5 * (arr[58] + arr[60])
    return _circular_running_mean(seas, smooth_days), _circular_running_mean(thresh, smooth_days)


def _run_length(mask: np.ndarray) -> np.ndarray:
    """For every True element, the length of the run of Trues it belongs to (along axis 0)."""
    fwd = np.zeros(mask.shape, dtype=np.int32)
    bwd = np.zeros(mask.shape, dtype=np.int32)
    fwd[0] = mask[0]
    for i in range(1, mask.shape[0]):
        fwd[i] = (fwd[i - 1] + 1) * mask[i]
    bwd[-1] = mask[-1]
    for i in range(mask.shape[0] - 2, -1, -1):
        bwd[i] = (bwd[i + 1] + 1) * mask[i]
    return (fwd + bwd - 1) * mask


def heatwave_days(
    sst: np.ndarray, time, thresh: np.ndarray, min_days: int = 5, max_gap_days: int = 2
) -> np.ndarray:
    """Boolean array, True on days that belong to a marine heatwave.

    Bridged gap days count as heatwave days, as in Hobday et al. (2016).
    Missing SST never exceeds the threshold.
    """
    doy = day_of_year_366(time)
    with np.errstate(invalid="ignore"):
        exceed = np.asarray(sst) > thresh[doy - 1]
    event = exceed & (_run_length(exceed) >= min_days)
    gap = ~event
    interior = ~np.logical_and.accumulate(gap, axis=0)
    interior &= ~np.logical_and.accumulate(gap[::-1], axis=0)[::-1]
    bridged = gap & interior & (_run_length(gap) <= max_gap_days)
    return event | bridged


def events(sst_1d, time, seas_1d, thresh_1d, min_days: int = 5, max_gap_days: int = 2) -> pd.DataFrame:
    """Event table for one pixel or one area-mean series: dates, duration and intensity."""
    t = pd.DatetimeIndex(pd.to_datetime(np.asarray(time)))
    sst_1d = np.asarray(sst_1d, dtype=float)
    hot = heatwave_days(sst_1d, t, np.asarray(thresh_1d), min_days, max_gap_days)
    anomaly = sst_1d - np.asarray(seas_1d)[day_of_year_366(t) - 1]
    edges = np.flatnonzero(np.diff(np.concatenate([[False], hot, [False]]).astype(int)))
    rows = []
    for start, stop in zip(edges[::2], edges[1::2]):
        a = anomaly[start:stop]
        rows.append(
            {
                "start": t[start], "end": t[stop - 1], "duration_days": int(stop - start),
                "intensity_max": float(np.nanmax(a)), "intensity_mean": float(np.nanmean(a)),
                "intensity_cumulative": float(np.nansum(a)),
            }
        )
    return pd.DataFrame(
        rows, columns=["start", "end", "duration_days", "intensity_max", "intensity_mean",
                       "intensity_cumulative"]
    )
