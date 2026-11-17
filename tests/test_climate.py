"""Climate statistics. Made-up numbers only: these check arithmetic, not results."""

import numpy as np
import pytest

from picopast.climate import corr_effective, fdr_bh, mann_kendall


def test_mann_kendall_and_sen():
    up = mann_kendall(np.arange(10.0) * 0.5)
    assert up["S"] == 45 and up["p"] < 0.001 and up["sen"] == pytest.approx(0.5)
    flat = mann_kendall([1.0, 3.0, 2.0, 3.0, 1.0, 2.0])
    assert flat["p"] > 0.5
    assert np.isnan(mann_kendall([1.0, 2.0])["p"])


def test_fdr_bh_matches_hand_calculation():
    p = np.array([0.01, 0.04, 0.03, 0.20, np.nan])
    # sorted 0.01, 0.03, 0.04, 0.20 -> 0.04, 0.0533, 0.0533, 0.20
    assert fdr_bh(p)[:4] == pytest.approx([0.04, 0.0533333, 0.0533333, 0.20], rel=1e-4)
    assert np.isnan(fdr_bh(p)[4])


def test_effective_sample_size_shrinks_for_smooth_series():
    rng = np.random.default_rng(0)
    white = corr_effective(rng.normal(size=200), rng.normal(size=200))
    assert white["n_eff"] > 150
    smooth = np.convolve(rng.normal(size=220), np.ones(20) / 20, mode="valid")[:200]
    other = np.convolve(rng.normal(size=220), np.ones(20) / 20, mode="valid")[:200]
    red = corr_effective(smooth, other)
    assert red["n_eff"] < 40 and red["n"] == 200
