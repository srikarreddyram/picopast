"""Virtual VIIRS (PRD section 6 step 4): PACE spectra seen through VIIRS band responses.

``Rrs_band = integral(Rrs * RSR) / integral(RSR)`` for each VIIRS band, using
the relative spectral response (RSR) of SNPP or NOAA-20.

This simple RSR weighting is not exactly how NASA forms band values (the
operational processing weights radiance and solar irradiance separately), and
PACE and VIIRS differ in calibration and overpass time. The PRD therefore
requires the residual against real same-day VIIRS pixels to be measured and
documented; ``residual_summary`` does that.

RSR tables: picopast.ingest.rsr.load_rsr (NOAA NCC releases, cached under data/raw/rsr).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: Nominal centres of the VIIRS visible ocean-colour bands M1-M5 [nm].
VIIRS_BANDS = {"snpp": (410, 443, 486, 551, 671), "noaa20": (411, 445, 489, 556, 667)}


def band_average(wavelength, rrs, rsr: pd.DataFrame) -> np.ndarray:
    """Band-average spectra of any leading shape.

    Parameters
    ----------
    wavelength : (w,) nm, ascending.
    rrs : (..., w) reflectance; a pixel with any NaN inside a band's response gives NaN there.
    rsr : frame indexed by wavelength [nm], one column per band, values >= 0.

    Returns
    -------
    (..., n_bands) array in the column order of ``rsr``.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    rrs = np.asarray(rrs, dtype=float)
    rsr_wave = rsr.index.to_numpy(dtype=float)
    out = np.empty(rrs.shape[:-1] + (rsr.shape[1],))
    for j, band in enumerate(rsr.columns):
        # Response sampled at the PACE wavelengths; zero outside the tabulated range.
        resp = np.interp(wavelength, rsr_wave, rsr[band].to_numpy(dtype=float), left=0.0, right=0.0)
        use = resp > 1e-4 * resp.max()
        if use.sum() < 2:
            raise ValueError(f"band {band}: the spectra do not sample its response")
        weight = resp[use] * np.gradient(wavelength)[use]
        out[..., j] = (rrs[..., use] * weight).sum(axis=-1) / weight.sum()
    return out


def residual_summary(virtual, real) -> dict[str, float]:
    """Virtual-minus-real VIIRS on same-day pixels, for one band."""
    v, r = np.asarray(virtual, dtype=float).ravel(), np.asarray(real, dtype=float).ravel()
    ok = np.isfinite(v) & np.isfinite(r)
    d = v[ok] - r[ok]
    if not d.size:
        return {"n": 0}
    return {
        "n": int(d.size),
        "bias": float(d.mean()),
        "rmsd": float(np.sqrt(np.mean(d**2))),
        "median_ratio": float(np.median(v[ok][r[ok] > 0] / r[ok][r[ok] > 0])),
        "r": float(np.corrcoef(v[ok], r[ok])[0, 1])
        if d.size > 2 and np.ptp(v[ok]) > 0 and np.ptp(r[ok]) > 0
        else float("nan"),
    }
