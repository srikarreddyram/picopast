"""Virtual VIIRS vs real VIIRS on same-day pixels (PRD section 6 step 4).

For every day with PACE and VIIRS daily files, PACE spectra are band-averaged
with each sensor's M1-M5 responses and compared with that sensor's own L3m
Rrs on the same pixels (offshore, both valid). The residual is documented, not
corrected: it mixes the band-averaging shortcut, calibration and
atmospheric-correction differences, and the 0-1 h gap between overpasses.
One day is held in memory at a time (a PACE box day is ~0.5 GB as float32).
"""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.ingest.rsr import load_rsr
from picopast.reconstruct import PREFIX
from picopast.virtual_viirs import VIIRS_BANDS, band_average, residual_summary

_DAY = re.compile(r"PACE_OCI\.(\d{8})\.L3m\.DAY\.")
NAMES = {"snpp": "SNPP", "noaa20": "NOAA-20"}


def _viirs_day(cfg: Config, sensor: str, stamp: str) -> np.ndarray | None:
    bands = []
    for b in VIIRS_BANDS[sensor]:
        paths = list((cfg.path("raw") / f"viirs_{sensor}_rrs").glob(
            f"{PREFIX[sensor]}.{stamp}.L3m.DAY.RRS.Rrs_{b}.4km.box.nc"))
        if not paths:
            return None
        with xr.open_dataset(paths[0]) as ds:
            bands.append(ds[f"Rrs_{b}"].values.astype(np.float32))
    return np.stack(bands, axis=-1)


def same_day_pairs(cfg: Config, sensor: str, per_day: int = 50_000, seed: int = 0) -> pd.DataFrame:
    """Sampled (virtual, real) pairs per band for every day with both PACE and VIIRS daily files."""
    from picopast.fallback_run import _offshore

    rsr = load_rsr(cfg, sensor)
    rng, rows, offshore = np.random.default_rng(seed), [], None
    for path in sorted((cfg.path("raw") / "pace_rrs").glob("PACE_OCI.*.L3m.DAY.*.box.nc")):
        stamp = _DAY.search(path.name).group(1)
        real = _viirs_day(cfg, sensor, stamp)
        if real is None:
            continue
        with xr.open_dataset(path) as ds:
            if offshore is None:
                offshore = _offshore(cfg, ds)
            virtual = band_average(ds["wavelength"].values, ds["Rrs"].values, rsr).astype(np.float32)
        ok = offshore & np.isfinite(virtual).all(-1) & np.isfinite(real).all(-1)
        idx = np.flatnonzero(ok.ravel())
        if idx.size > per_day:
            idx = rng.choice(idx, per_day, replace=False)
        v, r = virtual.reshape(-1, 5)[idx], real.reshape(-1, 5)[idx]
        day = pd.DataFrame({f"v{b}": v[:, j] for j, b in enumerate(VIIRS_BANDS[sensor])}
                           | {f"r{b}": r[:, j] for j, b in enumerate(VIIRS_BANDS[sensor])})
        day["date"] = pd.Timestamp(stamp)
        rows.append(day)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def report(cfg: Config) -> Path:
    L = ["# Virtual VIIRS vs real VIIRS on the same day", "",
         f"Generated {date.today().isoformat()} by `picopast virtual-viirs-check`. PACE OCI daily Rrs, band-averaged "
         "with each sensor's M1–M5 response (NOAA NCC releases, `picopast.ingest.rsr`), against that sensor's own "
         "daily L3m Rrs on the same 4.6 km pixels (offshore, both valid; up to 50,000 sampled pixels a day). The "
         "residual mixes the simple band-averaging, calibration and atmospheric-correction differences between "
         "missions, and the 0–1 h gap between overpasses. It is documented here, not corrected.", ""]
    for sensor in ("snpp", "noaa20"):
        pairs = same_day_pairs(cfg, sensor)
        if pairs.empty:
            L += [f"## {NAMES[sensor]}", "", "No days with both PACE and VIIRS daily files.", ""]
            continue
        out = cfg.path("interim") / "virtual_viirs" / f"pairs_{sensor}.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        pairs.to_parquet(out, index=False)
        days = pairs["date"].dt.date
        L += [f"## {NAMES[sensor]}", "",
              f"{days.nunique()} days ({days.min()} to {days.max()}), {len(pairs):,} pixel pairs.", "",
              "| band (nm) | median real Rrs (sr⁻¹) | bias (virtual − real) | bias / median | RMSD | median ratio | r |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
        for b in VIIRS_BANDS[sensor]:
            s = residual_summary(pairs[f"v{b}"], pairs[f"r{b}"])
            med = float(np.median(pairs[f"r{b}"]))
            L.append(f"| {b} | {med:.5f} | {s['bias']:+.5f} | {s['bias'] / med:+.1%} | {s['rmsd']:.5f} | "
                     f"{s['median_ratio']:.3f} | {s['r']:.3f} |")
        ratio = pairs[[f"v{b}" for b in VIIRS_BANDS[sensor][1:4]]].to_numpy()
        real = pairs[[f"r{b}" for b in VIIRS_BANDS[sensor][1:4]]].to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            v_ratio = np.log10(ratio[:, 0] / ratio[:, 2])
            r_ratio = np.log10(real[:, 0] / real[:, 2])
        ok = np.isfinite(v_ratio) & np.isfinite(r_ratio)
        L += ["", f"Blue-green ratio log10(Rrs{VIIRS_BANDS[sensor][1]}/Rrs{VIIRS_BANDS[sensor][3]}), the core of "
              f"chlorophyll and of the models' strongest inputs: bias {np.mean(v_ratio[ok] - r_ratio[ok]):+.3f}, "
              f"r {np.corrcoef(v_ratio[ok], r_ratio[ok])[0, 1]:.3f}.", ""]
    out = cfg.path("reports") / "virtual_viirs.md"
    out.write_text("\n".join(L))
    return out
