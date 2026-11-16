"""Run MOANA over mapped PACE Rrs, and check our code against NASA's Atlantic L4 files."""

from __future__ import annotations

import numpy as np
import xarray as xr

from picopast.grid import nearest_to_grid
from picopast.moana.algorithm import GROUPS, NASA_L4_NAMES, log10_abundance, run_moana


def moana_on_grid(
    rrs: xr.DataArray,
    sst_celsius: xr.DataArray | None,
    pc_mapping: str = "operational",
    rows_per_chunk: int = 64,
) -> xr.Dataset:
    """MOANA on a mapped PACE scene.

    ``rrs`` has dims (lat, lon, wavelength). ``sst_celsius`` has dims (lat, lon)
    on any regular grid; it is looked up at the Rrs pixel centres (nearest
    neighbour), which is how a 0.1-degree SST analysis meets a 4.6 km grid.
    """
    rrs = rrs.transpose("lat", "lon", "wavelength")
    lat, lon = rrs["lat"].values, rrs["lon"].values
    sst = None
    if sst_celsius is not None:
        sst = nearest_to_grid(sst_celsius, lat, lon).transpose("lat", "lon").values

    shape = (lat.size, lon.size)
    out = {g: np.full(shape, np.nan, dtype=np.float32) for g in GROUPS}
    flags = np.zeros(shape, dtype=np.uint8)
    wavelength = rrs["wavelength"].values
    for start in range(0, lat.size, rows_per_chunk):
        rows = slice(start, start + rows_per_chunk)
        res = run_moana(
            wavelength,
            rrs.isel(lat=rows).values,
            None if sst is None else sst[rows],
            pc_mapping=pc_mapping,
        )
        for g in GROUPS:
            out[g][rows] = res[g]
        flags[rows] = res["flags"]

    ds = xr.Dataset(
        {g: (("lat", "lon"), out[g], {"units": "cells mL-1"}) for g in GROUPS},
        coords={"lat": lat, "lon": lon},
    )
    ds["flags"] = (("lat", "lon"), flags)
    ds.attrs.update(algorithm="MOANA (PicoPast implementation)", pc_mapping=pc_mapping)
    return ds


def compare_with_nasa_l4(ours: xr.Dataset, nasa: xr.Dataset) -> dict[str, dict[str, float]]:
    """Per-group agreement between our MOANA and NASA's L4 file on the same grid and date.

    NASA clamps negative Prochlorococcus to zero and truncates all three outputs
    to integers, so we apply the same post-processing to ours before comparing.
    Agreement is reported as the absolute difference in log10(cells/mL) over
    pixels where both are positive.
    """
    report = {}
    for g in GROUPS:
        a = ours[g].values.astype(np.float64)
        b = nasa[NASA_L4_NAMES[g]].values.astype(np.float64)
        a = np.where(np.isfinite(a), np.trunc(np.maximum(a, 0.0)), np.nan)
        both = np.isfinite(a) & np.isfinite(b)
        la, lb = log10_abundance(a), log10_abundance(b)
        pos = np.isfinite(la) & np.isfinite(lb)
        d = np.abs(la[pos] - lb[pos])
        report[g] = {
            "n_both_valid": int(both.sum()),
            "n_compared": int(pos.sum()),
            "n_only_ours": int((np.isfinite(a) & ~np.isfinite(b)).sum()),
            "n_only_nasa": int((~np.isfinite(a) & np.isfinite(b)).sum()),
            "median_abs_dlog10": float(np.median(d)) if d.size else float("nan"),
            "p95_abs_dlog10": float(np.percentile(d, 95)) if d.size else float("nan"),
            "share_within_1pct": float(np.mean(d < np.log10(1.01))) if d.size else float("nan"),
        }
    return report
