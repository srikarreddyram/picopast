"""Checks for config, grid, metrics, matchup, trust gate, splits, virtual VIIRS and parsers."""

import io

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from picopast import metrics, trust_gate
from picopast.config import Box, load_config
from picopast.grid import block_id, box_axes, coast_mask, global_axes, nearest_to_grid, subset_box
from picopast.ingest.open_data import parse_dmi, parse_oni
from picopast.insitu.schema import era, standardise
from picopast.insitu.seabass import _in_box, cell_count_fields, read_sb
from picopast.matchup import match
from picopast.splits import spatial_block_folds, temporal_holdout
from picopast.virtual_viirs import band_average, residual_summary

CFG = load_config()


# --- config and grid -------------------------------------------------------------------------

def test_config_box_regions_seasons():
    assert CFG.box == Box(0.0, 25.0, 50.0, 100.0)
    assert set(CFG.regions(phase=1)) == {"bay_of_bengal", "arabian_sea_central_east"}
    assert len(CFG.regions()) == 4
    assert CFG.season_of_month(7) == "sw_monsoon" and CFG.season_of_month(1) == "ne_monsoon"
    for name, region in CFG.regions().items():
        assert CFG.box.contains(region.lat_min, region.lon_min), name
        assert CFG.box.contains(region.lat_max, region.lon_max), name


def test_grid_is_the_native_4km_grid():
    res = CFG["grid"]["resolution_deg"]
    lat, lon = global_axes(res)
    assert (lat.size, lon.size) == (4320, 8640)
    blat, blon = box_axes(CFG.box, res)
    assert (blat.size, blon.size) == (600, 1200)
    assert blat[0] > blat[-1]  # north to south, like the L3M files


def test_subset_box_handles_both_latitude_orders():
    lat = np.arange(-9.5, 40, 1.0)
    da = xr.DataArray(np.zeros((lat.size, 80)), dims=("lat", "lon"),
                      coords={"lat": lat, "lon": np.arange(40.5, 120, 1.0)})
    a = subset_box(da, CFG.box)
    b = subset_box(da.isel(lat=slice(None, None, -1)), CFG.box)
    assert a.shape == b.shape == (25, 50)


def test_nearest_to_grid_and_blocks():
    coarse = xr.DataArray([[1.0, 2.0], [3.0, 4.0]], dims=("lat", "lon"),
                          coords={"lat": [10.0, 11.0], "lon": [60.0, 61.0]})
    out = nearest_to_grid(coarse, np.array([10.1, 10.9]), np.array([60.2, 60.6]))
    assert out.values.tolist() == [[1.0, 2.0], [3.0, 4.0]]
    ids = block_id([10.1, 11.9, 12.1], [60.1, 61.9, 60.1], 2.0)
    assert ids[0] == ids[1] != ids[2]


def test_coast_mask_removes_a_25_km_band():
    res = CFG["grid"]["resolution_deg"]
    lat = 15.0 - np.arange(60) * res
    ocean = np.ones((60, 60), dtype=bool)
    ocean[:, :10] = False  # land on the western side
    keep = coast_mask(ocean, lat, res, mask_km=25.0)
    dx_km = np.radians(res) * 6371.0088 * np.cos(np.radians(lat.mean()))
    first_kept = int(np.argmax(keep[30]))
    assert not keep[:, :10].any()
    assert (first_kept - 9) * dx_km >= 25.0 > (first_kept - 10) * dx_km


# --- metrics and trust gate ------------------------------------------------------------------

def test_metrics_known_values():
    obs = np.array([3.0, 4.0, 5.0, 6.0, np.nan])
    est = obs + np.array([0.1, -0.1, 0.5, 0.2, 0.0])
    s = metrics.summary(est, obs)
    assert s["n"] == 4
    assert s["log_bias"] == pytest.approx(0.175)
    assert s["log_rmse"] == pytest.approx(np.sqrt((0.01 + 0.01 + 0.25 + 0.04) / 4))
    assert s["within_factor_2"] == pytest.approx(0.75)  # 0.5 > log10(2)
    assert s["spearman"] == pytest.approx(1.0)
    assert metrics.skill(est, obs, np.full(5, 4.5)) > 0
    assert metrics.interval_coverage(est - 0.3, est + 0.3, obs) == pytest.approx(0.75)


def _synthetic_matchups(noise, bias=0.0, n_cruises=6, per_cruise=12, seed=1):
    rng = np.random.default_rng(seed)
    obs = rng.uniform(3, 5.5, n_cruises * per_cruise)
    cruise = np.repeat(np.arange(n_cruises), per_cruise)
    est = obs + bias + rng.normal(0, noise, obs.size)
    return est, obs, cruise


def test_trust_gate_verdicts():
    good = trust_gate.evaluate(CFG, *_synthetic_matchups(0.15), n_boot=300)
    assert good["verdict"] == "usable_label"
    lo, hi = good["stats"]["log_rmse"]["lo"], good["stats"]["log_rmse"]["hi"]
    assert lo <= good["stats"]["log_rmse"]["value"] <= hi

    biased = trust_gate.evaluate(CFG, *_synthetic_matchups(0.15, bias=0.9), n_boot=300)
    assert biased["verdict"] == "relative_index_only"

    rng = np.random.default_rng(0)
    est, obs, cruise = _synthetic_matchups(0.15)
    useless = trust_gate.evaluate(CFG, rng.permutation(est), obs, cruise, n_boot=300)
    assert useless["verdict"] == "excluded"

    few = trust_gate.evaluate(CFG, est[:20], obs[:20], cruise[:20], n_boot=100)
    assert few["verdict"] == "insufficient_matchups"


def test_cruise_bootstrap_is_wider_than_a_station_bootstrap():
    # One cruise carries a large offset: resampling cruises must show that uncertainty.
    est, obs, cruise = _synthetic_matchups(0.05)
    est = est + np.where(cruise == 0, 0.8, 0.0)
    by_cruise = metrics.cruise_bootstrap(est, obs, cruise, n_boot=500)
    by_station = metrics.cruise_bootstrap(est, obs, np.arange(est.size), n_boot=500)
    width = lambda s: s["log_bias"]["hi"] - s["log_bias"]["lo"]  # noqa: E731
    assert width(by_cruise) > 1.5 * width(by_station)


# --- matchup ---------------------------------------------------------------------------------

def test_matchup_same_day_then_fallback_and_block_median():
    lat = np.array([12.0, 11.0, 10.0, 9.0, 8.0])
    lon = np.array([85.0, 86.0, 87.0, 88.0, 89.0])
    time = pd.date_range("2024-05-01", periods=3)
    data = np.full((3, 5, 5), np.nan)
    data[1] = np.arange(25.0).reshape(5, 5)  # 2 May: full coverage
    data[0, 1:4, 1:4] = 100.0  # 1 May: only the block around (10N, 87E)
    field = xr.DataArray(data, dims=("time", "lat", "lon"),
                         coords={"time": time, "lat": lat, "lon": lon})
    stations = pd.DataFrame(
        {
            "time_utc": pd.to_datetime(["2024-05-02 06:00", "2024-05-01 23:00", "2024-05-03 10:00",
                                        "2024-05-02 06:00", "2024-06-01 00:00"]),
            "lat": [10.1, 10.0, 10.0, 12.0, 10.0],
            "lon": [87.2, 87.0, 87.0, 87.0, 87.0],
        }
    )
    out = match(stations, field)
    assert out["sat_value"][0] == 12.0 and out["sat_day_offset"][0] == 0  # median of the 3x3 block
    assert out["sat_value"][1] == 100.0 and out["sat_day_offset"][1] == 0
    assert out["sat_value"][2] == 12.0 and out["sat_day_offset"][2] == -1  # fell back one day
    assert np.isnan(out["sat_value"][3])  # block would leave the map
    assert np.isnan(out["sat_value"][4])  # no map within +-1 day


def test_matchup_needs_enough_valid_pixels():
    field = xr.DataArray(
        np.full((1, 3, 3), np.nan), dims=("time", "lat", "lon"),
        coords={"time": [pd.Timestamp("2024-05-01")], "lat": [11.0, 10.0, 9.0],
                "lon": [86.0, 87.0, 88.0]},
    )
    field.values[0, 0, :] = 1.0
    field.values[0, 1, 0] = 1.0  # 4 valid pixels: below the minimum of 5
    stations = pd.DataFrame({"time_utc": [pd.Timestamp("2024-05-01")], "lat": [10.0], "lon": [87.0]})
    assert np.isnan(match(stations, field)["sat_value"][0])
    field.values[0, 1, 1] = 1.0
    assert match(stations, field)["sat_value"][0] == 1.0


# --- splits ----------------------------------------------------------------------------------

def test_spatial_blocks_never_straddle_folds():
    rng = np.random.default_rng(2)
    lat, lon = rng.uniform(5, 22, 4000), rng.uniform(80, 95, 4000)
    blocks = block_id(lat, lon, 2.0)
    seen = np.zeros(lat.size, dtype=int)
    for train, test in spatial_block_folds(lat, lon, n_folds=5):
        assert not set(blocks[train]) & set(blocks[test])
        seen[test] += 1
    assert (seen == 1).all()


def test_temporal_holdout_takes_the_latest_six_months():
    t = pd.date_range("2024-03-05", "2026-10-31", freq="D")
    train, test = temporal_holdout(t, months=6)
    assert t[test].min() > t[train].max()
    assert t[test].min() == pd.Timestamp("2026-05-01")
    assert train.size + test.size == t.size


# --- virtual VIIRS ---------------------------------------------------------------------------

def test_band_average_recovers_flat_and_linear_spectra():
    wave = np.arange(400.0, 700.0, 2.5)
    rsr_wave = np.arange(400.0, 700.0, 1.0)
    rsr = pd.DataFrame(
        {c: np.exp(-0.5 * ((rsr_wave - c) / 8.0) ** 2) for c in (443, 551)}, index=rsr_wave
    )
    flat = np.full((2, 3, wave.size), 0.004)
    np.testing.assert_allclose(band_average(wave, flat, rsr), 0.004)
    linear = 1e-5 * wave  # a symmetric response returns the value at the band centre
    np.testing.assert_allclose(band_average(wave, linear, rsr), [1e-5 * 443, 1e-5 * 551], rtol=1e-3)
    gap = flat.copy()
    gap[0, 0, np.argmin(np.abs(wave - 445))] = np.nan
    out = band_average(wave, gap, rsr)
    assert np.isnan(out[0, 0, 0]) and np.isfinite(out[0, 0, 1])
    assert residual_summary(out[..., 1], flat[..., 0])["bias"] == pytest.approx(0.0, abs=1e-12)


# --- parsers and schema ----------------------------------------------------------------------

def test_parse_indices():
    dmi = parse_dmi(" 2025 2026\n2025 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2\n"
                    "2026 0.5 -9999.000 " + "-9999.000 " * 10 + "\n-9999\nDMI HadISST1.1\n")
    assert len(dmi) == 13 and dmi["date"].iloc[-1] == pd.Timestamp("2026-01-01")
    oni = parse_oni(" SEAS  YR   TOTAL   ANOM\n  DJF 1950  25.01  -1.32\n  NDJ 1950  25.36  -1.20\n")
    assert oni["date"].tolist() == [pd.Timestamp("1950-01-01"), pd.Timestamp("1950-12-01")]
    assert oni["oni"].tolist() == [-1.32, -1.20]


SB_SAMPLE = """/begin_header
/investigators=Test_Person
/experiment=DEMO
/cruise=demo_cruise
/missing=-9999
/delimiter=comma
/fields=date,time,lat,lon,depth,abun_prochlorococcus,abun_synechococcus,abun_picoeukaryotes
/units=yyyymmdd,hh:mm:ss,degrees,degrees,m,cell/ml,cell/ml,cell/ml
/end_header
20240429,06:10:00,14.1165,85.0092,5,150000,42000,2100
20240430,05:55:00,16.7473,87.0128,5,-9999,38000,1800
"""


def test_read_seabass_file():
    header, df = read_sb(io.StringIO(SB_SAMPLE))
    assert header["cruise"] == "demo_cruise" and df.shape == (2, 8)
    assert np.isnan(df["abun_prochlorococcus"][1])
    assert cell_count_fields(header) == {
        "pro_cells_ml": "abun_prochlorococcus",
        "syn_cells_ml": "abun_synechococcus",
        "peuk_cells_ml": "abun_picoeukaryotes",
    }


def test_cmr_footprint_test():
    assert _in_box("14.1165 85.0092", CFG.box) is True
    assert _in_box("2.0 82.6 17.8 89.5", CFG.box) is True
    assert _in_box("40.0 -70.0", CFG.box) is False
    assert _in_box("", CFG.box) is None


def test_schema_standardise_and_era():
    raw = pd.DataFrame(
        {
            "source": "x", "dataset_id": "d", "cruise": "c",
            "time_utc": ["1995-03-20", "2015-06-01", "2024-04-29", "2024-04-29", None],
            "lat": [17.3, 10.0, 14.1, 14.1, 5.0], "lon": [67.9, 88.0, 85.0, 85.0, 60.0],
            "depth_m": [5, 5, 5, 5, 5],
            "pro_cells_ml": [1e5, -999, 2e5, 2e5, None], "syn_cells_ml": [1e4, None, 3e4, 3e4, None],
        }
    )
    df = standardise(raw)
    assert len(df) == 2  # all-missing rows and the duplicate are dropped
    assert era(df["time_utc"]).tolist() == ["pre_viirs", "pace_era"]
