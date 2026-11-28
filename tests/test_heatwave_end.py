"""Storms-end-heatwaves arithmetic on made-up arrays: these check the bookkeeping, not results."""

import numpy as np
import pandas as pd
import pytest

from picopast import heatwave_end as he


def _cube(T: int = 300, ny: int = 8, nx: int = 8) -> he.Cube:
    c = object.__new__(he.Cube)
    c.lat = 10.0 + 0.25 * np.arange(ny)
    c.lon = 80.0 + 0.25 * np.arange(nx)
    c.ocean = np.ones((ny, nx), bool)
    c.time = pd.date_range("1982-01-01", periods=T, freq="D")
    c.T, c.nx = T, nx
    c.flat = np.zeros((T, ny * nx), np.uint8)
    c.ssta = None
    c.la2, c.lo2 = np.meshgrid(c.lat, c.lon, indexing="ij")
    return c


def test_mantel_haenszel_single_stratum():
    n1, a1, n0, a0 = np.array([100.0]), np.array([50.0]), np.array([1000.0]), np.array([100.0])
    r = he._mh(n0, a0, n1, a1, a0 + a1)
    assert r["hr"] == pytest.approx(5.0)                  # 0.5 per day against 0.1
    assert r["caused"] == pytest.approx(40.0)              # 50 seen, 10 expected at the unexposed rate
    assert r["share"] == pytest.approx(40 / 150)
    assert r["share_exposed"] == pytest.approx(40 / 50)


def test_endings_and_exposure_are_counted_once():
    c = _cube()
    cell = 3 * c.nx + 3
    c.flat[10:20, cell] = 1                                # one heatwave: days 10-19, ends after day 19
    group = np.zeros(c.flat.shape[1], int)
    totals = he.daily_totals(c, group, 1)
    assert totals[:, 0, 0].sum() == 10 and totals[:, 0, 1].sum() == 1
    assert totals[19, 0, 1] == 1
    fixes = pd.DataFrame({"time": [c.time[20] + pd.Timedelta(hours=6), c.time[20] + pd.Timedelta(hours=12)],
                          "lat": [c.lat[3]] * 2, "lon": [c.lon[3]] * 2, "wind": [70.0, 40.0]})
    t, cells, w = he.exposed_transitions(c, fixes, "wind", 50.0, (-1, 1))
    assert set(t) == {19, 20, 21}                           # a fix on day 20 exposes t = 19, 20, 21
    assert np.all(w == 70.0)                               # the strongest wind is kept, each pair once
    assert len(set(zip(t, cells))) == t.size
    C = he.count_array(c, group, 1, totals, t, cells, w)
    S = C.sum(axis=(0, 1))[0]
    assert S[0] == 10 and S[1] == 1 and S[2] == 1 and S[3] == 1 and S[5] == 1


def test_cell_stats_survival_return_and_occupancy():
    c = _cube()
    k0 = 50
    cells = np.array([0, 1, 2, 3])
    c.flat[k0 - 3:k0 + 8, 0] = 1                           # still hot on D0+5
    c.flat[k0 - 3:k0 + 2, 1] = 1                           # ended by D0+5, back on D0+40
    c.flat[k0 + 40:k0 + 50, 1] = 1
    c.flat[k0 - 3:k0 + 2, 2] = 1                           # ended, never back: capped at 240
    st = he.cell_stats(c.flat, k0, cells)
    assert st["n_hot"] == 3 and st["curve"][4] == pytest.approx(1 / 3)
    assert st["n_ended"] == 2 and st["ret"] == pytest.approx((40 + 240) / 2)
    assert st["occ_not|month 1"] == 0.0
    assert he.weighted(0.75, st["occ_hot|month 1"], st["occ_not|month 1"]) == pytest.approx(
        0.75 * st["occ_hot|month 1"])
    assert np.isnan(he.weighted(0.5, np.nan, 0.1)) and he.weighted(0.0, np.nan, 0.1) == pytest.approx(0.1)


def test_storm_bootstrap_mean_and_slope():
    rng = np.random.default_rng(1)
    sids = np.repeat(np.arange(40), 5)
    x = rng.normal(size=sids.size)
    y = -0.5 * x + rng.normal(scale=0.1, size=sids.size)
    e = he.boot_mean(y, sids, 500)
    assert e["lo"] <= e["effect"] <= e["hi"] and e["n_storms"] == 40
    s = he.boot_slope(y, np.column_stack([np.ones_like(x), x]), sids, 500)
    assert s["hi"] < 0 and s["slope"] == pytest.approx(-0.5, abs=0.05)
    assert he.s2_verdict({"lo": 2.0, "hi": 8.0}).endswith("(by less than 10 days)")
    assert he.s2_verdict({"lo": -3.0, "hi": 4.0}) == "no meaningful difference"
