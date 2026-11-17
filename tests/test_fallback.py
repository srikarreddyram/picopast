"""Fallback checks, pigment station-days and the MOANA box runner's file handling.

Made-up numbers only: these check arithmetic and thresholds, not results.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from picopast.config import Box
from picopast.fallback import climatology_verdict, pigment_proxy_verdict, spearman_bootstrap
from picopast.insitu.pigments import station_days
from picopast.moana.box_run import scene_period

RULE_A = {"min_station_days": 15, "min_spearman": 0.40, "min_spearman_ci_low": 0.0}
RULE_C = {"min_stations": 30, "min_cruises": 3, "consistent_ci_low": 0.0, "inconsistent_ci_high": 0.30}


def test_spearman_bootstrap_monotone_and_clusters():
    x = np.arange(40.0)
    s = spearman_bootstrap(x, x**2, cluster=np.repeat(np.arange(8), 5), n_boot=200)
    assert s["rho"] == pytest.approx(1.0)
    assert s["lo"] == pytest.approx(1.0) and s["n"] == 40 and s["n_clusters"] == 8
    one = spearman_bootstrap(x, -x, cluster=np.zeros(40), n_boot=200)
    assert one["rho"] == pytest.approx(-1.0) and np.isnan(one["lo"])  # one cluster: no interval


def test_pigment_proxy_verdict():
    good = {"rho": 0.6, "lo": 0.2, "hi": 0.8, "n": 20, "n_clusters": 20}
    assert pigment_proxy_verdict(good, RULE_A)[0].startswith("relative index supported")
    assert pigment_proxy_verdict({**good, "n": 14}, RULE_A)[0] == "not supported"
    assert pigment_proxy_verdict({**good, "rho": 0.39}, RULE_A)[0] == "not supported"
    assert pigment_proxy_verdict({**good, "lo": 0.0}, RULE_A)[0] == "not supported"  # strictly above 0
    assert pigment_proxy_verdict({**good, "rho": float("nan")}, RULE_A)[0] == "not supported"


def test_climatology_verdict():
    base = {"rho": 0.3, "lo": 0.1, "hi": 0.5, "n": 40, "n_clusters": 5}
    assert climatology_verdict(base, RULE_C)[0] == "pattern consistent"
    assert climatology_verdict({**base, "lo": -0.4, "hi": 0.2}, RULE_C)[0] == "pattern not consistent"
    assert climatology_verdict({**base, "lo": -0.1, "hi": 0.6}, RULE_C)[0] == "inconclusive"
    assert climatology_verdict({**base, "n_clusters": 2}, RULE_C)[0] == "inconclusive"


def test_station_days_merges_replicates_and_filters():
    t = pd.Timestamp("2024-05-01 06:00")
    s = pd.DataFrame({
        "sample_id": ["a", "a", "b", "c", "d", "e", "f"],
        "cruise": "X", "station": ["1", "1", "1b", "2", "3", "4", "5"],
        "time_utc": [t, t, t, t, t, pd.Timestamp("2024-01-01"), t],
        "lat": [10.01, 10.01, 10.012, 12.01, 30.0, 10.01, 12.01],  # 10.0 is a pixel edge
        "lon": [85.01, 85.01, 85.012, 86.01, 85.0, 85.01, 86.01],
        "depth_m": [0, 0, 5, 0, 0, 0, 50],
        "dv_chl_a": [0.02, 0.02, 0.04, 0.10, 0.5, 0.5, 0.9],
        "zea": 0.05, "tchl": 0.2,
    })
    box = Box(lat_min=0, lat_max=25, lon_min=50, lon_max=100)
    out = station_days(s, box, start="2024-03-05")
    # "a" duplicate dropped; a+b share a pixel and day -> one row (median 0.03);
    # 30N is outside the box, January is before PACE, 50 m is not surface.
    assert len(out) == 2
    merged = out[out["n_samples"] == 2].iloc[0]
    assert merged["dv_chl_a"] == pytest.approx(0.03)
    assert set(out["dv_chl_a"].round(3)) == {0.03, 0.10}


def test_scene_period():
    d = scene_period(Path("PACE_OCI.20250115.L3m.DAY.AOP.V3_2.4km.box.nc"))
    assert d == (pd.Timestamp("2025-01-15"), pd.Timestamp("2025-01-15"), "DAY")
    e = scene_period(Path("PACE_OCI.20240305_20240312.L3m.8D.AOP.V3_2.4km.box.nc"))
    assert e == (pd.Timestamp("2024-03-05"), pd.Timestamp("2024-03-12"), "8D")


def test_historical_counts_drops_bcodmo_twins():
    from picopast.insitu.historical import historical_counts

    t = "1995-01-10T05:00:00"
    pangaea = pd.DataFrame({"source": "pangaea", "cruise": "TT043", "time_utc": [t, t],
                            "lat": [21.1826, 21.1826], "lon": [63.5491, 63.5491], "depth_m": [2.87, 40.0],
                            "syn_cells_ml": [45188.0, 1000.0]})
    bcodmo = pd.DataFrame({"source": "bcodmo", "cruise": "bcodmo_1995-01",
                           "time_utc": ["1995-01-10", "1995-01-10", "2001-05-01"],
                           "lat": [21.183, 21.183, 15.0], "lon": [63.549, 63.549, 65.0],
                           "depth_m": [2.9, 20.0, 5.0], "syn_cells_ml": [45188.0, 9.0, 7.0]})
    out = historical_counts(bcodmo, pangaea)
    # PANGAEA holds the 1995-01-10 station, so both BCO-DMO samples from it go, whatever their
    # depth (PANGAEA's TT054 depths are whole metres, BCO-DMO's are not); the 2001 one stays.
    assert len(out) == 3
    assert out.loc[out["source"] == "bcodmo", "depth_m"].tolist() == [5.0]


def test_water_mask_uses_land_flag_not_sst_gaps(tmp_path):
    """CMC fills SST over land, so the coast must come from the mask's land flag."""
    import xarray as xr

    from picopast.fallback_run import water_mask

    lat = np.arange(10.05, 9.0, -0.1)
    lon = np.arange(80.05, 81.0, 0.1)
    flags = np.ones((1, lat.size, lon.size), dtype=np.float32)
    flags[:, :, :3] = 2.0  # land on the western edge
    raw = tmp_path / "raw" / "cmc_sst"
    raw.mkdir(parents=True)
    xr.Dataset({"analysed_sst": (("time", "lat", "lon"), np.full(flags.shape, 300.0, np.float32)),
                "mask": (("time", "lat", "lon"), flags)},
               coords={"time": [np.datetime64("2025-01-01T12")], "lat": lat, "lon": lon}
               ).to_netcdf(raw / "20250101120000-CMC-test.box.nc")

    class Cfg:
        def path(self, name):
            return tmp_path / name

    grid_lat, grid_lon = np.array([9.5]), np.array([80.1, 80.5, 80.9])
    np.testing.assert_array_equal(water_mask(Cfg(), grid_lat, grid_lon), [[False, True, True]])


def test_decision_follows_preregistered_rule():
    from picopast.fallback_run import decision

    v = {("A", "pro"): "relative index supported (provisional, pigment proxy)",
         ("C", "pro"): "pattern consistent", ("C", "syn"): "inconclusive", ("C", "peuk"): "pattern consistent"}
    d = decision(v)
    assert d["pro"][0] and not d["syn"][0] and d["peuk"][0]
    # Prochlorococcus needs A as well as C.
    assert not decision({**v, ("A", "pro"): "not supported"})["pro"][0]
    assert not decision({k: x for k, x in v.items() if k != ("A", "pro")})["pro"][0]


def test_paired_difference_and_partial():
    from picopast.fallback import paired_spearman_difference, partial_spearman

    rng = np.random.default_rng(1)
    z = rng.normal(size=300)                 # shared driver ("chlorophyll")
    a = z + rng.normal(scale=0.3, size=300)  # "MOANA": z plus its own signal
    extra = rng.normal(size=300)
    y = z + 0.8 * extra                       # counts depend on z and on something else
    a_extra = a + 0.8 * extra
    d = paired_spearman_difference(a_extra, z, y, np.repeat(np.arange(30), 10), n_boot=300)
    assert d["diff"] > 0 and d["diff_lo"] > 0 and d["partial"] > 0.3
    assert abs(partial_spearman(a, y, z)) < 0.15  # a adds nothing beyond z here
