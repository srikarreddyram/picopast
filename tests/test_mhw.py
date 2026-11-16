"""Marine heatwave detection: run logic on hand-built series, climatology on a known cycle."""

import numpy as np
import pandas as pd

from picopast import mhw


def test_leap_calendar_day_of_year():
    doy = mhw.day_of_year_366(["2023-02-28", "2023-03-01", "2024-02-29", "2024-03-01", "2023-12-31"])
    assert doy.tolist() == [59, 61, 60, 61, 366]


def test_event_rules_min_duration_and_gap_bridging():
    time = pd.date_range("2021-01-01", periods=60)
    thresh = np.full(366, 28.0)
    sst = np.full(60, 27.0)
    sst[2:6] = 29.0     # 4 days: too short
    sst[10:15] = 29.0   # 5 days: an event
    sst[17:22] = 29.0   # 5 days after a 2-day gap: bridged into one 12-day event
    sst[30:36] = 29.0   # 6 days
    sst[39:44] = 29.0   # 3-day gap: stays a separate event
    sst[57:60] = 29.0   # 3 days at the end: too short
    hot = mhw.heatwave_days(sst, time, thresh)
    expected = np.zeros(60, dtype=bool)
    expected[10:22] = expected[30:36] = expected[39:44] = True
    assert hot.tolist() == expected.tolist()

    table = mhw.events(sst, time, np.full(366, 27.0), thresh)
    assert table["duration_days"].tolist() == [12, 6, 5]
    assert table["start"][0] == pd.Timestamp("2021-01-11")
    assert table["intensity_max"][1] == 2.0
    assert table["intensity_cumulative"][2] == 10.0


def test_detection_is_vectorised_and_ignores_missing_values():
    time = pd.date_range("2021-01-01", periods=30)
    sst = np.full((30, 2, 2), 27.0)
    sst[5:12, 0, 0] = 30.0
    sst[5:12, 1, 1] = 30.0
    sst[8, 1, 1] = np.nan  # a missing day splits the run into 3 + 3 days: no event
    hot = mhw.heatwave_days(sst, time, np.full((366, 2, 2), 28.0))
    assert hot[:, 0, 0].sum() == 7
    assert hot[:, 1, 1].sum() == 0 and hot[:, 0, 1].sum() == 0


def test_climatology_recovers_a_known_seasonal_cycle():
    time = pd.date_range("1991-01-01", "2020-12-31")
    doy = mhw.day_of_year_366(time)
    rng = np.random.default_rng(0)
    cycle = 28.0 + 2.0 * np.sin(2 * np.pi * (doy - 100) / 366)
    sst = (cycle + rng.normal(0, 0.5, time.size))[:, None]
    seas, thresh = mhw.climatology(sst, time)
    assert seas.shape == thresh.shape == (366, 1) and np.isfinite(seas).all()
    truth = 28.0 + 2.0 * np.sin(2 * np.pi * (np.arange(1, 367) - 100) / 366)
    assert np.abs(seas[:, 0] - truth).max() < 0.1
    # The 90th percentile of N(0, 0.5) noise sits 1.28 * 0.5 = 0.64 above the mean.
    assert np.abs((thresh - seas)[:, 0] - 0.64).max() < 0.15
    # About 10% of baseline days exceed the threshold, by construction.
    exceed = sst[:, 0] > thresh[doy - 1, 0]
    assert 0.07 < exceed.mean() < 0.13
