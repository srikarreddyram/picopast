"""Validation metrics (PRD sections 3, 6.3 and 7). All abundance errors are in log10 space."""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def _pair(est_log10, obs_log10):
    est = np.asarray(est_log10, dtype=float)
    obs = np.asarray(obs_log10, dtype=float)
    ok = np.isfinite(est) & np.isfinite(obs)
    return est[ok], obs[ok]


def log_rmse(est_log10, obs_log10) -> float:
    e, o = _pair(est_log10, obs_log10)
    return float(np.sqrt(np.mean((e - o) ** 2))) if e.size else float("nan")


def log_bias(est_log10, obs_log10) -> float:
    """Mean of (estimate - observation) in log10 units; 10**bias is the typical ratio."""
    e, o = _pair(est_log10, obs_log10)
    return float(np.mean(e - o)) if e.size else float("nan")


def spearman(est_log10, obs_log10) -> float:
    e, o = _pair(est_log10, obs_log10)
    if e.size < 3 or np.ptp(e) == 0 or np.ptp(o) == 0:
        return float("nan")
    return float(spearmanr(e, o).statistic)


def within_factor(est_log10, obs_log10, factor: float = 2.0) -> float:
    """Share of pairs whose ratio lies within [1/factor, factor]."""
    e, o = _pair(est_log10, obs_log10)
    return float(np.mean(np.abs(e - o) <= np.log10(factor))) if e.size else float("nan")


def summary(est_log10, obs_log10) -> dict[str, float]:
    e, o = _pair(est_log10, obs_log10)
    return {
        "n": int(e.size),
        "log_rmse": log_rmse(e, o),
        "log_bias": log_bias(e, o),
        "spearman": spearman(e, o),
        "within_factor_2": within_factor(e, o, 2.0),
    }


def cruise_bootstrap(
    est_log10, obs_log10, cruise, n_boot: int = 2000, ci: float = 0.95, seed: int = 0
) -> dict[str, dict[str, float]]:
    """Point estimates with intervals from a bootstrap that resamples whole cruises.

    Stations on one cruise share water masses, instruments and analysts, so they
    are not independent; resampling cruises keeps that dependence intact.
    """
    est = np.asarray(est_log10, dtype=float)
    obs = np.asarray(obs_log10, dtype=float)
    cruise = np.asarray(cruise)
    ok = np.isfinite(est) & np.isfinite(obs)
    est, obs, cruise = est[ok], obs[ok], cruise[ok]
    point = summary(est, obs)
    labels = np.unique(cruise)
    members = [np.flatnonzero(cruise == c) for c in labels]
    rng = np.random.default_rng(seed)
    draws = {k: [] for k in point if k != "n"}
    for _ in range(n_boot):
        idx = np.concatenate([members[i] for i in rng.integers(0, len(labels), len(labels))])
        s = summary(est[idx], obs[idx])
        for k in draws:
            draws[k].append(s[k])
    lo_q, hi_q = (1 - ci) / 2 * 100, (1 + ci) / 2 * 100
    out = {"n": {"value": point["n"]}, "n_cruises": {"value": int(len(labels))}}
    for k, values in draws.items():
        v = np.asarray(values, dtype=float)
        v = v[np.isfinite(v)]
        lo, hi = (np.percentile(v, [lo_q, hi_q]) if v.size else (np.nan, np.nan))
        out[k] = {"value": point[k], "lo": float(lo), "hi": float(hi)}
    return out


def skill(est_log10, obs_log10, reference_log10) -> float:
    """Skill = 1 - RMSE_model / RMSE_reference (PRD section 7); > 0 beats the reference."""
    est, obs, ref = (np.asarray(x, dtype=float) for x in (est_log10, obs_log10, reference_log10))
    ok = np.isfinite(est) & np.isfinite(obs) & np.isfinite(ref)
    if not ok.any():
        return float("nan")
    return float(1.0 - log_rmse(est[ok], obs[ok]) / log_rmse(ref[ok], obs[ok]))


def interval_coverage(lower, upper, obs) -> float:
    """Share of observations inside [lower, upper]; a 90% interval should give 0.85-0.95."""
    lower, upper, obs = (np.asarray(x, dtype=float) for x in (lower, upper, obs))
    ok = np.isfinite(lower) & np.isfinite(upper) & np.isfinite(obs)
    return float(np.mean((obs[ok] >= lower[ok]) & (obs[ok] <= upper[ok]))) if ok.any() else float("nan")
