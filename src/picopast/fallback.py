"""Trust-gate fallback checks for when cell counts are missing (pre-registration, checks A-C).

These checks compare ranks only. They can support a group as a *relative
index*; they can never give it a usable cells/mL label. Thresholds come from
``trust_gate.pigment_proxy`` and ``trust_gate.climatology_pattern`` in
``config/picopast.yaml``, fixed on 2026-11-16 before any comparison.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def _rho(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan")
    return float(spearmanr(x, y).statistic)


def spearman_bootstrap(x, y, cluster, n_boot: int = 2000, ci: float = 0.95, seed: int = 0) -> dict:
    """Spearman rho with a percentile interval from resampling whole clusters.

    ``cluster`` is the resampling unit: cruise for check C, station-day for
    check A (two cruises are too few to resample).
    """
    x, y, cluster = np.asarray(x, float), np.asarray(y, float), np.asarray(cluster)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, cluster = x[ok], y[ok], cluster[ok]
    labels = np.unique(cluster)
    members = [np.flatnonzero(cluster == c) for c in labels]
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot if labels.size > 1 else 0):
        idx = np.concatenate([members[i] for i in rng.integers(0, labels.size, labels.size)])
        draws.append(_rho(x[idx], y[idx]))
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    lo, hi = np.percentile(d, [(1 - ci) / 2 * 100, (1 + ci) / 2 * 100]) if d.size else (np.nan, np.nan)
    return {"rho": _rho(x, y), "lo": float(lo), "hi": float(hi), "n": int(x.size),
            "n_clusters": int(labels.size)}


def partial_spearman(x, y, z) -> float:
    """Rank correlation of x and y with z held fixed (partial correlation of the ranks)."""
    rxy, rxz, ryz = _rho(x, y), _rho(x, z), _rho(y, z)
    return float((rxy - rxz * ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2)))


def paired_spearman_difference(a, b, y, cluster, n_boot: int = 2000, ci: float = 0.95, seed: int = 0) -> dict:
    """rho(a, y) - rho(b, y) on the same points, with an interval from resampling whole clusters.

    Also gives the partial rank correlation of a and y with b held fixed: what ``a``
    tells about ``y`` beyond ``b``.
    """
    a, b, y, cluster = np.asarray(a, float), np.asarray(b, float), np.asarray(y, float), np.asarray(cluster)
    ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(y)
    a, b, y, cluster = a[ok], b[ok], y[ok], cluster[ok]
    labels = np.unique(cluster)
    members = [np.flatnonzero(cluster == c) for c in labels]
    rng = np.random.default_rng(seed)
    diff, part = [], []
    for _ in range(n_boot if labels.size > 1 else 0):
        idx = np.concatenate([members[i] for i in rng.integers(0, labels.size, labels.size)])
        diff.append(_rho(a[idx], y[idx]) - _rho(b[idx], y[idx]))
        part.append(partial_spearman(a[idx], y[idx], b[idx]))
    q = [(1 - ci) / 2 * 100, (1 + ci) / 2 * 100]

    def interval(d):
        d = np.asarray(d, float)
        d = d[np.isfinite(d)]
        return np.percentile(d, q) if d.size else (np.nan, np.nan)

    (dlo, dhi), (plo, phi) = interval(diff), interval(part)
    return {"rho_a": _rho(a, y), "rho_b": _rho(b, y), "diff": _rho(a, y) - _rho(b, y), "diff_lo": float(dlo),
            "diff_hi": float(dhi), "partial": partial_spearman(a, y, b), "partial_lo": float(plo),
            "partial_hi": float(phi), "n": int(y.size), "n_clusters": int(labels.size)}


def pigment_proxy_verdict(stats: dict, rule: dict) -> tuple[str, list[str]]:
    """Check A: MOANA Prochlorococcus vs HPLC divinyl chlorophyll a."""
    reasons = []
    if stats["n"] < rule["min_station_days"]:
        reasons.append(f"{stats['n']} matched station-days < {rule['min_station_days']}")
    if not stats["rho"] >= rule["min_spearman"]:
        reasons.append(f"rho {stats['rho']:.2f} < {rule['min_spearman']}")
    if not stats["lo"] > rule["min_spearman_ci_low"]:
        reasons.append(f"interval lower end {stats['lo']:.2f} <= {rule['min_spearman_ci_low']}")
    if reasons:
        return "not supported", reasons
    return "relative index supported (provisional, pigment proxy)", []


def climatology_verdict(stats: dict, rule: dict) -> tuple[str, list[str]]:
    """Check C: same-month PACE-era MOANA climatology vs 1991-2005 ship counts."""
    if stats["n"] < rule["min_stations"] or stats["n_clusters"] < rule["min_cruises"]:
        return "inconclusive", [f"{stats['n']} stations from {stats['n_clusters']} cruises "
                                f"(need {rule['min_stations']} from {rule['min_cruises']})"]
    if stats["lo"] > rule["consistent_ci_low"]:
        return "pattern consistent", []
    if stats["hi"] < rule["inconsistent_ci_high"]:
        return "pattern not consistent", [f"interval upper end {stats['hi']:.2f} < {rule['inconsistent_ci_high']}"]
    return "inconclusive", [f"interval [{stats['lo']:.2f}, {stats['hi']:.2f}] spans the thresholds"]
