"""Match ship stations to satellite maps (PRD section 6 step 3).

For each station: take the satellite map of the same day; if it has too few
valid pixels around the station, try the day before, then the day after. The
matched value is the median of the 3x3 pixel block centred on the station.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def match(
    stations: pd.DataFrame,
    field: xr.DataArray,
    window_px: int = 3,
    min_valid_px: int = 5,
    day_offsets: tuple[int, ...] = (0, -1, 1),
) -> pd.DataFrame:
    """Satellite value for each station.

    Parameters
    ----------
    stations : frame with ``time_utc``, ``lat``, ``lon`` (other columns are kept).
    field : daily maps with dims (time, lat, lon); any latitude order.

    Returns
    -------
    ``stations`` plus ``sat_value`` (median of valid pixels, NaN if no match),
    ``sat_n_valid`` and ``sat_day_offset`` (0, -1, +1, or NaN if no match).
    """
    field = field.transpose("time", "lat", "lon")
    lat, lon = field["lat"].values, field["lon"].values
    dlat, dlon = abs(lat[1] - lat[0]), abs(lon[1] - lon[0])
    days = pd.DatetimeIndex(field["time"].values).normalize()
    day_index = {d: i for i, d in enumerate(days)}
    half = window_px // 2
    data = field.values  # read once; one PACE-era box of one variable fits in memory

    value = np.full(len(stations), np.nan)
    n_valid = np.zeros(len(stations), dtype=int)
    offset_used = np.full(len(stations), np.nan)
    for k, (t, y, x) in enumerate(zip(stations["time_utc"], stations["lat"], stations["lon"])):
        if pd.isna(t) or not np.isfinite(y) or not np.isfinite(x):
            continue
        iy, ix = int(np.argmin(np.abs(lat - y))), int(np.argmin(np.abs(lon - x)))
        # The station must lie inside its nearest pixel, and the block inside the map.
        if abs(lat[iy] - y) > dlat / 2 or abs(lon[ix] - x) > dlon / 2:
            continue
        if iy < half or ix < half or iy >= lat.size - half or ix >= lon.size - half:
            continue
        day = pd.Timestamp(t).normalize()
        for off in day_offsets:
            it = day_index.get(day + pd.Timedelta(days=off))
            if it is None:
                continue
            block = data[it, iy - half : iy + half + 1, ix - half : ix + half + 1]
            good = np.isfinite(block)
            if good.sum() >= min_valid_px:
                value[k], n_valid[k], offset_used[k] = np.median(block[good]), good.sum(), off
                break
    return stations.assign(sat_value=value, sat_n_valid=n_valid, sat_day_offset=offset_used)
