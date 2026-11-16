"""Our own MOANA run (PRD section 6 step 2).

MOANA (Lange et al. 2020; NASA ATBD v1.2) maps hyperspectral Rrs, plus SST for
*Prochlorococcus* only, to cell abundances of three picophytoplankton groups:

1. interpolate Rrs linearly onto 414, 416, ..., 660 nm (124 bands);
2. standardise each spectrum across wavelength (mean 0, sample sd 1);
3. project onto NASA's PCA loadings: ``U_i = sum_lambda Rrs'(lambda) V[lambda, i]``;
4. evaluate one multilinear regression per group.

Things the PRD summary glosses over, taken from NASA's operational tables
(``resources/moana``, checksum-pinned) and the ATBD equations:

* *Synechococcus* and picoeukaryotes are regressed in log10(cells/mL);
  **Prochlorococcus is regressed in linear cells/mL** with a log10(SST) term.
  Its raw output can therefore be negative. We keep the raw value and raise
  ``FLAG_NEGATIVE_PRO`` instead of clamping to zero as NASA's code does.
* The ATBD equations and the operational coefficient file disagree on the
  principal component of two terms (Pro: U17 vs U7; Syn: U13 vs U16). The
  EPFT-UP project showed the shipping PACE product follows the operational
  file, so ``pc_mapping="operational"`` is the default; ``"atbd"`` is kept so
  the difference can be reported.

MOANA output is a model estimate, never truth.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import numpy as np

GROUPS = ("pro", "syn", "peuk")
GROUP_NAMES = {
    "pro": "Prochlorococcus",
    "syn": "Synechococcus",
    "peuk": "picoeukaryotes",
}
#: Variable names in NASA's L4 MOANA files, for the Atlantic code check.
NASA_L4_NAMES = {"pro": "prococcus_moana", "syn": "syncoccus_moana", "peuk": "picoeuk_moana"}

FLAG_NEGATIVE_PRO = 1  # Pro regression went negative (value kept, not clamped)
FLAG_MISSING_RRS = 2  # a band needed for 414-660 nm was missing; no retrieval
FLAG_BAD_SST = 4  # SST missing or <= 0 degC; Pro not retrieved
FLAG_FLAT_SPECTRUM = 8  # zero spectral variance; cannot standardise

# SHA-256 of the two NASA OCSSW files (tag T2023.31, bundle "common"),
# downloaded from oceandata.sci.gsfc.nasa.gov/manifest on 2026-11-10.
_LUT_SHA256 = {
    "picophyt.json": "cd6b5c1cca879bc97acf6c5a3d4f47bf658d5a839b99ff4b1f534eac881756d4",
    "pca_picophyto.h5": "478ace8798ce525bc6a09567b3773cadd532056c13d20ce6daf0a529f46deb57",
}


@dataclass(frozen=True)
class MoanaLuts:
    wavelength: np.ndarray  # (124,) nm
    loadings: np.ndarray  # (124, npc) PCA loadings V
    npc: int
    intercept: dict[str, float]
    sst_coef: float  # Pro only, multiplies log10(SST in degC)
    weights: dict[str, np.ndarray]  # per group, (npc,) weight on U_1..U_npc


def _floats(text: str) -> np.ndarray:
    return np.array([float(x) for x in text.split(",")], dtype=np.float64)


@lru_cache(maxsize=1)
def load_luts() -> MoanaLuts:
    """Read NASA's MOANA tables shipped with the package, verifying checksums."""
    import h5py

    base = resources.files("picopast.resources.moana")
    for name, expected in _LUT_SHA256.items():
        got = hashlib.sha256((base / name).read_bytes()).hexdigest()
        if got != expected:
            raise RuntimeError(f"{name} does not match the pinned NASA checksum")

    coef = json.loads((base / "picophyt.json").read_text())
    npc = int(coef["npc"])
    pro, syn, peuk = (_floats(coef[k]) for k in ("pro_coef", "syn_coef", "apeuk_coef"))

    with resources.as_file(base / "pca_picophyto.h5") as path, h5py.File(path, "r") as f:
        wavelength = f["wavelength"][:].astype(np.float64)
        loadings = f["component"][:, :npc].astype(np.float64)

    # Slot layout from OCSSW get_Cpicophyt.c: Pro = [intercept, log10 SST, U1..U17];
    # Syn and picoeukaryotes = [intercept, U1..U17].
    return MoanaLuts(
        wavelength=wavelength,
        loadings=loadings,
        npc=npc,
        intercept={"pro": pro[0], "syn": syn[0], "peuk": peuk[0]},
        sst_coef=pro[1],
        weights={"pro": pro[2 : 2 + npc], "syn": syn[1 : 1 + npc], "peuk": peuk[1 : 1 + npc]},
    )


def _weights(luts: MoanaLuts, pc_mapping: str) -> dict[str, np.ndarray]:
    if pc_mapping == "operational":
        return luts.weights
    if pc_mapping == "atbd":
        w = {g: v.copy() for g, v in luts.weights.items()}
        # ATBD v1.2 puts the last Pro term on U17 (not U7) and the last Syn term on U13 (not U16).
        for group, src, dst in (("pro", 7, 17), ("syn", 16, 13)):
            w[group][dst - 1], w[group][src - 1] = w[group][src - 1], 0.0
        return w
    raise ValueError("pc_mapping must be 'operational' or 'atbd'")


def interp_matrix(wave_in, wave_out) -> np.ndarray:
    """Linear-interpolation weights W so that ``rrs_out = rrs_in @ W.T``.

    One matrix for the whole grid replaces a per-pixel ``np.interp`` call; the
    two are identical for spectra with no missing bands.
    """
    wave_in = np.asarray(wave_in, dtype=np.float64)
    wave_out = np.asarray(wave_out, dtype=np.float64)
    if np.any(np.diff(wave_in) <= 0):
        raise ValueError("input wavelengths must be strictly ascending")
    if wave_in[0] > wave_out[0] or wave_in[-1] < wave_out[-1]:
        raise ValueError(
            f"input wavelengths {wave_in[0]:g}-{wave_in[-1]:g} nm do not cover "
            f"{wave_out[0]:g}-{wave_out[-1]:g} nm; MOANA is not extrapolated"
        )
    hi = np.clip(np.searchsorted(wave_in, wave_out, side="left"), 1, wave_in.size - 1)
    lo = hi - 1
    frac = (wave_out - wave_in[lo]) / (wave_in[hi] - wave_in[lo])
    w = np.zeros((wave_out.size, wave_in.size))
    rows = np.arange(wave_out.size)
    w[rows, lo] = 1.0 - frac
    w[rows, hi] += frac
    return w


def run_moana(wavelength, rrs, sst=None, pc_mapping: str = "operational") -> dict[str, np.ndarray]:
    """MOANA on spectra of any leading shape.

    Parameters
    ----------
    wavelength : (w,) nm, ascending, covering 414-660 nm.
    rrs : (..., w) remote-sensing reflectance [sr^-1]; NaN where missing.
    sst : (...) sea-surface temperature in **degrees Celsius**, or None.
        Needed for Prochlorococcus only.
    pc_mapping : "operational" (as NASA ships it) or "atbd".

    Returns
    -------
    dict with ``pro``, ``syn``, ``peuk`` in cells/mL (shape ``rrs.shape[:-1]``;
    ``pro`` may be negative, see module docstring), ``flags`` (uint8 bitmask)
    and ``scores`` (..., npc).
    """
    luts = load_luts()
    weights = _weights(luts, pc_mapping)
    rrs = np.asarray(rrs, dtype=np.float64)
    lead = rrs.shape[:-1]
    flat = rrs.reshape(-1, rrs.shape[-1])
    flags = np.zeros(flat.shape[0], dtype=np.uint8)

    w = interp_matrix(wavelength, luts.wavelength)
    used = w.any(axis=0)  # input bands that contribute to 414-660 nm
    ok = np.isfinite(flat[:, used]).all(axis=1)
    flags[~ok] |= FLAG_MISSING_RRS

    scores = np.full((flat.shape[0], luts.npc), np.nan)
    if ok.any():
        spec = flat[ok][:, used] @ w[:, used].T
        mean = spec.mean(axis=1, keepdims=True)
        sd = spec.std(axis=1, ddof=1, keepdims=True)  # sample sd, as in OCSSW
        flat_spec = sd[:, 0] == 0
        with np.errstate(invalid="ignore", divide="ignore"):
            std = (spec - mean) / sd
        s = std @ luts.loadings
        s[flat_spec] = np.nan
        scores[ok] = s
        idx = np.flatnonzero(ok)
        flags[idx[flat_spec]] |= FLAG_FLAT_SPECTRUM

    syn = 10.0 ** (luts.intercept["syn"] + scores @ weights["syn"])
    peuk = 10.0 ** (luts.intercept["peuk"] + scores @ weights["peuk"])

    if sst is None:
        log_sst = np.full(flat.shape[0], np.nan)
        flags |= FLAG_BAD_SST
    else:
        sst = np.broadcast_to(np.asarray(sst, dtype=np.float64), lead).reshape(-1).copy()
        finite = sst[np.isfinite(sst)]
        if finite.size and np.median(finite) > 60.0:
            raise ValueError("SST looks like Kelvin; MOANA needs degrees Celsius")
        bad = ~np.isfinite(sst) | (sst <= 0)
        flags[bad] |= FLAG_BAD_SST
        sst[bad] = np.nan
        log_sst = np.log10(sst)
    pro = luts.intercept["pro"] + luts.sst_coef * log_sst + scores @ weights["pro"]
    flags[np.isfinite(pro) & (pro < 0)] |= FLAG_NEGATIVE_PRO

    return {
        "pro": pro.reshape(lead),
        "syn": syn.reshape(lead),
        "peuk": peuk.reshape(lead),
        "flags": flags.reshape(lead),
        "scores": scores.reshape(*lead, luts.npc),
    }


def log10_abundance(cells_ml) -> np.ndarray:
    """log10(cells/mL), NaN where the abundance is missing or not positive."""
    x = np.asarray(cells_ml, dtype=np.float64)
    out = np.full(x.shape, np.nan)
    pos = np.isfinite(x) & (x > 0)
    out[pos] = np.log10(x[pos])
    return out
