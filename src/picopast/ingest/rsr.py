"""VIIRS relative spectral response (RSR) tables for the ocean-colour bands M1-M5.

The PRD's link (oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables) now
redirects to a general page. The band-averaged RSRs come instead from NOAA's
calibration centre (NCC), which publishes them openly:

- SNPP (flight model F1): Northrop Grumman October 2011 "Fused" band-averaged
  RSR, the set used operationally in IDPS.
- NOAA-20 (J1): VIIRS Data Analysis Working Group (DAWG) at-launch public
  release V2.1 (November 2016), "Fused" band-averaged V2F files.

The archives are downloaded and cached under data/raw/rsr, not shipped in the
repository. Columns are labelled with NASA's nominal band centres, the names its
L3m Rrs files use.
"""

from __future__ import annotations

import io
import zipfile

import numpy as np
import pandas as pd

from picopast.config import Config
from picopast.virtual_viirs import VIIRS_BANDS

SOURCES = {
    "snpp": {
        "url": "https://ncc.nesdis.noaa.gov/documents/documentation/NG_VIIRS_NPP_RSR_filtered_Oct2011_BA.zip",
        "member": "NG_VIIRS_NPP_RSR_filtered_Oct2011_BA/NPP_VIIRS_NG_RSR_M{n}_Fused_Oct2011f_BA.dat",
        "release": "Northrop Grumman Oct 2011 fused band-averaged RSR (SNPP F1), via NOAA NCC",
    },
    "noaa20": {
        "url": "https://ncc.nesdis.noaa.gov/NOAA-20/docs/J1_VIIRS_RSR_DAWG_At-Launch_Public_Release_V2.1_Nov2016.zip",
        "member": "J1_VIIRS_BA_RSR_V2F/J1_VIIRS_RSR_M{n}_BA_Fused_V2F.txt",
        "release": "DAWG J1 VIIRS at-launch public release V2.1, fused band-averaged V2F, via NOAA NCC",
    },
}
GRID = np.round(np.arange(340.0, 760.0001, 0.1), 1)  # nm


def parse_rsr_text(text: str) -> tuple[np.ndarray, np.ndarray]:
    """(wavelength nm, response) from one band file: '%' comment lines, last two columns numeric."""
    wl, resp = [], []
    for line in text.splitlines():
        parts = line.split()
        if not parts or line.lstrip().startswith(("%", "#")):
            continue
        try:
            w, r = float(parts[-2]), float(parts[-1])
        except (ValueError, IndexError):
            continue
        wl.append(w)
        resp.append(r)
    wl, resp = np.asarray(wl), np.asarray(resp)
    order = np.argsort(wl)
    return wl[order], np.clip(resp[order], 0.0, None)


def table_from_zip(blob: bytes, sensor: str) -> pd.DataFrame:
    """M1-M5 responses from a release archive on one 0.1 nm grid, each band peak-normalised to 1."""
    spec, cols = SOURCES[sensor], {}
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for n, label in zip(range(1, 6), VIIRS_BANDS[sensor]):
            wl, resp = parse_rsr_text(z.read(spec["member"].format(n=n)).decode("latin-1"))
            col = np.interp(GRID, wl, resp, left=0.0, right=0.0)
            cols[label] = col / col.max()
    return pd.DataFrame(cols, index=pd.Index(GRID, name="wavelength_nm"))


def centre_wavelength(rsr: pd.DataFrame) -> pd.Series:
    """Response-weighted mean wavelength of each band [nm]."""
    w = rsr.index.to_numpy(dtype=float)
    return pd.Series({c: float((w * rsr[c]).sum() / rsr[c].sum()) for c in rsr.columns})


def load_rsr(cfg: Config, sensor: str) -> pd.DataFrame:
    """The sensor's M1-M5 RSR table, downloading NOAA's release once."""
    out = cfg.path("raw") / "rsr" / f"{sensor}_M1-M5.csv"
    if not out.exists():
        import requests

        r = requests.get(SOURCES[sensor]["url"], timeout=120, headers={"User-Agent": "picopast research client"})
        r.raise_for_status()
        table = table_from_zip(r.content, sensor)
        out.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out)
        (out.parent / f"{sensor}_SOURCE.txt").write_text(f"{SOURCES[sensor]['release']}\n{SOURCES[sensor]['url']}\n")
    t = pd.read_csv(out, index_col=0)
    t.columns = [int(c) for c in t.columns]
    return t
