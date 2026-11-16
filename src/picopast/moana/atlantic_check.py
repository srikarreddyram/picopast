"""Check our MOANA code against NASA's own Atlantic L4 files (PRD section 6 step 2).

For one date: read PACE L3M Rrs and CMC SST over a patch of the Atlantic, run
our MOANA, and compare pixel by pixel with NASA's L4 MOANA file for the same
day and grid. Needs the Earthdata login.

Not yet run against live data (no login when it was written); the first run
is the test of the file-layout assumptions noted inline.
"""

from __future__ import annotations

import xarray as xr

from picopast.config import Box, Config
from picopast.grid import subset_box
from picopast.ingest import earthdata
from picopast.moana.grid_run import compare_with_nasa_l4, moana_on_grid

#: A subtropical North Atlantic patch well inside NASA's MOANA domain (70S-70N, 85W-25E).
ATLANTIC_PATCH = Box(lat_min=20.0, lat_max=30.0, lon_min=-60.0, lon_max=-40.0)


def run(cfg: Config, day: str, patch: Box = ATLANTIC_PATCH, pc_mapping: str = "operational") -> dict:
    """Agreement report (see ``compare_with_nasa_l4``) for one date, e.g. "2025-01-15"."""
    kw = {"region": patch, "tag": "atl"}
    rrs_path = earthdata.fetch(cfg, "pace_rrs", day, day, "DAY", **kw)[0]
    l4_path = earthdata.fetch(cfg, "pace_moana_l4", day, day, "DAY", **kw)[0]
    sst_path = earthdata.fetch(cfg, "cmc_sst", day, day, **kw)[0]

    with xr.open_dataset(rrs_path) as rrs_ds, xr.open_dataset(l4_path) as nasa, xr.open_dataset(
        sst_path
    ) as sst_ds:
        # Assumed layouts: AOP has Rrs(lat, lon, wavelength); CMC has
        # analysed_sst(time, lat, lon) in kelvin.
        sst_c = sst_ds["analysed_sst"].isel(time=0) - 273.15
        ours = moana_on_grid(rrs_ds["Rrs"], sst_c, pc_mapping=pc_mapping)
        nasa = subset_box(nasa, patch)
        report = compare_with_nasa_l4(ours, nasa)
    report["_inputs"] = {"rrs": rrs_path.name, "l4": l4_path.name, "sst": sst_path.name,
                         "pc_mapping": pc_mapping}
    return report
