"""Storm-table helpers. Made-up numbers only: these check arithmetic and selection rules, not results."""

import numpy as np
import pandas as pd
import pytest

from picopast.storms import _period_values, disc, haversine_km, pre_post, season


def test_haversine_one_degree_of_latitude():
    assert haversine_km(0, 80, 1, 80) == pytest.approx(111.19, abs=0.05)
    assert haversine_km(10, 80, 10, 80) == 0.0


def test_disc_keeps_cells_inside_the_radius_and_the_valid_mask():
    lat = np.arange(0.125, 25, 0.25)
    lon = np.arange(50.125, 100, 0.25)
    iy, ix = disc(lat, lon, 15.0, 85.0, 150.0)
    d = haversine_km(15.0, 85.0, lat[iy], lon[ix])
    assert iy.size > 50 and (d <= 150.0).all()
    valid = np.ones((lat.size, lon.size), bool)
    valid[:, lon > 85.0] = False  # land east of 85E
    iy2, ix2 = disc(lat, lon, 15.0, 85.0, 150.0, valid)
    assert (lon[ix2] < 85.0).all() and iy2.size < iy.size
    assert disc(lat, lon, 15.0, 120.0, 150.0)[0].size == 0  # far outside the grid


def test_season_classes():
    assert list(season([2, 3, 6, 7, 9, 12])) == ["other", "pre-monsoon", "pre-monsoon", "other", "post-monsoon",
                                                 "post-monsoon"]


def test_pre_post_periods_skip_the_one_holding_the_storm():
    starts = pd.DatetimeIndex(["2020-05-08", "2020-05-16", "2020-05-24", "2020-06-01"])
    ends = starts + pd.Timedelta(days=7)
    # storm on 20 May: the 16-23 May period holds it, so pre = 8-15 May, post = 24 May and 1 June
    assert pre_post(starts, ends, pd.Timestamp("2020-05-20")) == (0, 2, 3)
    # storm on 17 May: 8-15 May ends on D0-2, so it is still a clean pre period
    assert pre_post(starts, ends, pd.Timestamp("2020-05-17")) == (0, 2, 3)
    # nothing before, one period after
    assert pre_post(starts, ends, pd.Timestamp("2020-05-28")) == (1, 3, -1)


def test_period_values_need_enough_valid_pixels():
    field = np.full((2, 5), np.nan)
    field[0, :] = [1.0, 2.0, 3.0, np.nan, np.nan]
    idx = {(7, 150.0): np.array([0, 1, 2, 3, 4]), (8, 150.0): np.array([5, 6, 7, 8, 9])}
    needs = {0: [(7, 150.0, "pre"), (8, 150.0, "pre")]}
    out = _period_values([(0, field)], needs, idx, min_share=0.5)
    assert out[(7, 150.0, "pre")] == pytest.approx(2.0)  # 3 of 5 valid
    assert np.isnan(out[(8, 150.0, "pre")])  # 0 of 5 valid


# ---------------------------------------------------------------- statistics (made-up data)

from picopast.storms_analysis import effect, km_median, storm_draws, verdict_bloom  # noqa: E402


def _fake_table(effect_size, n_storms=40, seed=3):
    rng = np.random.default_rng(seed)
    rows = []
    for s in range(n_storms):
        hot = s % 2  # half the storms over heatwave water
        storm_shift = rng.normal(0, 0.05)
        for _ in range(8):
            wind = rng.uniform(35, 120)
            rows.append({"sid": f"S{s}", "heatwave": float(hot), "wind": wind, "speed_kmh": rng.uniform(5, 25),
                         "lat": rng.uniform(5, 20), "chl_pre": rng.normal(0, 0.1), "ssha_cm": rng.normal(0, 5),
                         "dist2land_km": rng.uniform(50, 600), "basin": "Bay of Bengal" if s % 3 else "Arabian Sea",
                         "season": "pre-monsoon" if s % 4 else "post-monsoon",
                         "bloom": 0.002 * wind + effect_size * hot + storm_shift + rng.normal(0, 0.03)})
    return pd.DataFrame(rows)


def test_effect_recovers_a_planted_heatwave_effect():
    e = effect(_fake_table(-0.10), "bloom", n_boot=300)
    assert e["effect"] == pytest.approx(-0.10, abs=0.03)
    assert e["hi"] < 0 and e["storms_heatwave"] == 20 and e["storms_normal"] == 20
    assert verdict_bloom(e).startswith("weaker bloom")
    null = effect(_fake_table(0.0, seed=4), "bloom", n_boot=300)
    assert null["lo"] < 0 < null["hi"] and verdict_bloom(null) == "no detectable difference"


def test_too_few_storms_is_not_a_verdict():
    e = effect(_fake_table(-0.10, n_storms=20), "bloom", n_boot=100)
    assert e["storms_heatwave"] == 10 and verdict_bloom(e) == "not enough data"


def test_storm_draws_keep_storms_whole():
    sids = np.array(["a", "a", "b", "c", "c", "c"])
    for idx in storm_draws(sids, n_boot=20, seed=1):
        picked = sids[idx]
        for s in set(picked):
            assert (picked == s).sum() % (sids == s).sum() == 0


def test_km_median_with_censoring():
    t = np.array([2, 4, 6, 8, 60, 60])
    ev = np.array([1, 1, 1, 1, 0, 0], bool)
    assert km_median(t, ev) == 6.0  # S: 5/6, 4/6, 3/6 -> reaches 0.5 at day 6
    assert np.isnan(km_median(np.array([60, 60]), np.array([0, 0], bool)))


def test_mean_ci_and_group_verdicts():
    from picopast.storms_analysis import mean_ci, verdict_group

    sids = np.repeat(np.arange(30), 4)
    rng = np.random.default_rng(5)
    up = mean_ci(0.2 + rng.normal(0, 0.05, sids.size), sids, n_boot=200)
    assert up["lo"] > 0 and up["n"] == 120 and up["storms"] == 30
    assert np.isnan(mean_ci([1.0, np.nan], [1, 2])["mean"])
    assert verdict_group({"storms_heatwave": 20, "storms_normal": 20, "lo": -0.3, "hi": -0.1}).startswith("weaker")
    assert verdict_group({"storms_heatwave": 9, "storms_normal": 20, "lo": -0.3, "hi": -0.1}) == "not enough data"
