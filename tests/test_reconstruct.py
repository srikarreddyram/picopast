"""Reconstruction helpers. Made-up numbers only: these check arithmetic, not results."""

import numpy as np
import pandas as pd
import pytest

from picopast.reconstruct import SeasonalClimatology, features, to_target, viirs_periods


def test_to_target_units():
    assert to_target("pro", [250_000.0, -50_000.0]) == pytest.approx([2.5, -0.5])  # negatives kept
    out = to_target("syn", [1000.0, 0.0, -5.0])
    assert out[0] == pytest.approx(3.0) and np.isnan(out[1]) and np.isnan(out[2])


def test_features_ratios_need_positive_reflectance():
    v = {f"rrs{b}": np.array([0.004, 0.004], np.float32) for b in (410, 443, 486, 551, 671)}
    v["rrs551"] = np.array([0.002, 0.002], np.float32)
    v["rrs671"] = np.array([0.0002, -0.0001], np.float32)  # red band can be negative
    v["par"] = np.array([40.0, 40.0], np.float32)
    f = features(v, np.array([28.0, 28.0], np.float32), np.array([10.0, 11.0]), np.array([80.0, 81.0]),
                 pd.Timestamp("2025-01-05"))
    assert f["r443_551"].tolist() == pytest.approx([np.log10(2)] * 2)
    assert f["r671_551"].iloc[0] == pytest.approx(-1.0) and np.isnan(f["r671_551"].iloc[1])
    assert f["mbr"].tolist() == pytest.approx([np.log10(2)] * 2)
    assert f["doy_cos"].iloc[0] == pytest.approx(np.cos(2 * np.pi * 4.5 / 365.25), rel=1e-5)


def test_seasonal_climatology_falls_back_to_month_mean():
    train = pd.DataFrame({"month": [1, 1, 1, 2], "lat": [1.0, 1.0, 11.0, 1.0], "lon": [51.0, 51.0, 61.0, 51.0]})
    clim = SeasonalClimatology(cell_deg=2.5).fit(train, np.array([1.0, 3.0, 8.0, 5.0]))
    test = pd.DataFrame({"month": [1, 1, 2, 7], "lat": [1.5, 21.0, 11.0, 1.0], "lon": [51.5, 90.0, 61.0, 51.0]})
    # same cell and month; empty cell -> January mean; empty cell -> February mean; unseen month -> overall
    assert clim.predict(test) == pytest.approx([2.0, 4.0, 5.0, 4.25])


def test_viirs_periods_need_all_bands(tmp_path):
    raw = tmp_path / "raw"
    for d in ("viirs_snpp_rrs", "viirs_snpp_par"):
        (raw / d).mkdir(parents=True)
    for stamp, bands in (("20240101_20240108", (410, 443, 486, 551, 671)), ("20240109_20240116", (410, 443))):
        (raw / "viirs_snpp_par" / f"SNPP_VIIRS.{stamp}.L3m.8D.PAR.par.4km.box.nc").touch()
        for b in bands:
            (raw / "viirs_snpp_rrs" / f"SNPP_VIIRS.{stamp}.L3m.8D.RRS.Rrs_{b}.4km.box.nc").touch()

    class Cfg:
        def path(self, name):
            return tmp_path / name

    p = viirs_periods(Cfg())
    assert p["start"].tolist() == [pd.Timestamp("2024-01-01")] and p["end"].tolist() == [pd.Timestamp("2024-01-08")]
