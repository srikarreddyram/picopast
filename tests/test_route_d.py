"""Route D building blocks. Made-up numbers only."""

import numpy as np
import pytest

from picopast.route_d import (
    basin,
    effective_n,
    importance_weights,
    in_indian_ocean,
    season,
    window_median,
)


def test_indian_ocean_exclusion_and_basins():
    lat = np.array([10.0, -40.0, 35.0, 20.0, 40.0, 0.0])
    lon = np.array([65.0, 80.0, 65.0, -40.0, 15.0, -150.0])
    assert in_indian_ocean(lat, lon).tolist() == [True, True, False, False, False, False]
    assert basin(lat, lon).tolist() == ["indian", "indian", "other", "atlantic", "other", "other"]


def test_season_flips_south_of_the_equator():
    s_n, c_n = season([1], [30.0])
    s_s, c_s = season([7], [-30.0])
    assert s_n[0] == pytest.approx(s_s[0], abs=0.02) and c_n[0] == pytest.approx(c_s[0], abs=0.02)


def test_window_median_wraps_longitude_and_skips_gaps():
    glat = np.array([1.0, 0.0, -1.0])
    glon = np.arange(0.5, 360, 1.0)
    arr = np.full((3, 360), np.nan)
    arr[1, 359] = 4.0
    arr[1, 0] = 2.0
    out = window_median(arr, glat, glon, [0.0, 0.0], [0.4, 180.0])
    assert out[0] == pytest.approx(3.0)  # 359.5 and 0.5 are neighbours across the date line
    assert np.isnan(out[1])


def test_importance_weights_favour_target_like_samples():
    rng = np.random.default_rng(0)
    train = np.r_[rng.normal(0, 1, (500, 2)), rng.normal(4, 1, (500, 2))]
    target = rng.normal(4, 1, (500, 2))
    w = importance_weights(train, target)
    assert w.mean() == pytest.approx(1.0)
    assert w[500:].mean() > 5 * w[:500].mean()
    assert 100 < effective_n(w) < 1000


def test_mission_climatology_choice():
    from picopast.ingest.route_d_inputs import mission_climatologies

    names = ["AQUA_MODIS.20030101_20250131.L3m.MC.CHL.chlor_a.9km.nc",
             "AQUA_MODIS.20030101_20260131.L3m.MC.CHL.chlor_a.9km.nc",
             "AQUA_MODIS.20100101_20260131.L3m.MC.CHL.chlor_a.9km.nc",  # not mission-long
             "MODISA_L3m_CHL_AQUA_MODIS.20020701_20250731.L3m.MC.CHL.chlor_a.9km.nc_2022.0"]
    got = mission_climatologies(names)
    assert got[1] == names[1] and got[7] == names[3] and len(got) == 2
