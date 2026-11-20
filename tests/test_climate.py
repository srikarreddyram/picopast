"""Climate statistics. Made-up numbers only: these check arithmetic, not results."""

import numpy as np
import pandas as pd
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


# ---------------------------------------------------------------- pixel-level maps

from picopast.climate_maps import all_8d_starts, coarsen, composite, heatwave_classes, index_classes  # noqa: E402


def test_coarsen_needs_enough_pixels():
    v = np.arange(144, dtype=float).reshape(12, 12)
    ok = np.ones_like(v, bool)
    ok[:6, :6] = False
    ok[0, 0] = True  # one pixel only: below the minimum
    out = coarsen(v, ok, cell=6, min_pixels=9)
    assert out.shape == (2, 2) and np.isnan(out[0, 0])
    assert out[1, 1] == pytest.approx(v[6:, 6:].mean())


def test_composite_recovers_a_planted_difference():
    rng = np.random.default_rng(1)
    T, years = 200, np.repeat(np.arange(2010, 2020), 20)
    anom = rng.normal(0, 0.3, (T, 3, 4))
    a = np.zeros(T, bool)
    a[::4] = True
    b = np.zeros(T, bool)
    b[2::4] = True
    anom[a, 0, 0] += 1.0  # only cell (0, 0) responds
    c = composite(anom, years, a, b, n_boot=300)
    assert c["diff"][0, 0] == pytest.approx(1.0, abs=0.2)
    assert c["lo"][0, 0] > 0
    assert abs(c["diff"][2, 3]) < 0.25 and c["lo"][2, 3] < 0 < c["hi"][2, 3]
    few = composite(anom, years, a & (years == 2010), b, n_boot=50)  # one year only
    assert np.isnan(few["diff"]).all()


def test_index_lead_and_heatwave_classes():
    idx = pd.Series([1.0, 0.0, -1.0], index=pd.to_datetime(["2015-01-01", "2015-02-01", "2015-03-01"]))
    time = pd.to_datetime(["2015-03-06", "2015-05-01"])
    pos, neg = index_classes(idx, time, 0.5, lead=2)
    assert list(pos) == [True, False] and list(neg) == [False, True]  # March <- January, May <- March
    frac = np.array([[0.6], [0.05], [0.05], [0.3]]).reshape(4, 1, 1)
    t = pd.to_datetime(["2015-06-02", "2016-06-10", "2016-07-01", "2017-06-18"])
    hot, normal = heatwave_classes(frac, t)
    assert hot[:, 0, 0].tolist() == [True, False, False, False]
    assert normal[:, 0, 0].tolist() == [False, True, False, False]  # July has no heatwave here


def test_8d_starts():
    s = all_8d_starts(2015, 2016)
    assert len(s) == 92 and s[45] == pd.Timestamp("2015-12-27") and s[46] == pd.Timestamp("2016-01-01")
