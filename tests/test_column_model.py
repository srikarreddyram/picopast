"""Column model physics and biology on made-up profiles: these check behaviour, not results."""

import numpy as np
import pytest

from picopast.column_model import (
    Biology,
    Column,
    _static,
    biology_step,
    physics_step,
    storm_stress,
    vortex_wind,
)


def _column(warm_cap: float = 0.0, cap_m: float = 20.0, lat: float = 15.0) -> Column:
    z = np.arange(0.5, 200, 1.0)
    T = 29.0 - 10.0 * np.clip((z - 30) / 120, 0, 1)        # 30 m mixed layer over a thermocline
    T = T + warm_cap * (z < cap_m)                          # a heatwave-like warm cap
    S = np.full(z.size, 35.0)
    N = 0.05 + 20.0 * np.clip((z - 50) / 100, 0, 1)        # nitracline below 50 m
    P = 0.1 * np.exp(-((z - 60) / 15) ** 2) + 0.01
    return Column(z, T, S, N, P, np.full(z.size, 0.05), np.full(z.size, 0.05), lat)


def _blow(col: Column, tau: float, hours: int = 24, dt: float = 900.0) -> Column:
    for _ in range(int(hours * 3600 / dt)):
        physics_step(col, tau, 0.0, dt)
    return col


def test_static_instability_is_removed():
    col = _column()
    col.T[:10] -= 3.0                       # cold (dense) water on top
    _static(col)
    rho = col.density()
    assert np.all(np.diff(rho[:60]) >= -1e-9)


def test_wind_deepens_and_cools_more_when_stronger():
    weak, strong = _blow(_column(), 0.3), _blow(_column(), 3.0)
    assert strong.mld() > weak.mld() > 25
    assert strong.T[0] < weak.T[0] < 29.0
    assert strong.N[0] > weak.N[0]          # deeper mixing reaches the nitracline


def test_a_warm_cap_resists_mixing():
    plain, capped = _blow(_column(), 1.0, hours=12), _blow(_column(warm_cap=1.5), 1.0, hours=12)
    assert capped.mld() <= plain.mld()


def test_biology_conserves_nitrogen_and_needs_light():
    bio = Biology(w_sink=0.0)
    col = _column()
    total0 = (col.N + col.P + col.Z + col.D).sum()
    for _ in range(48):
        biology_step(col, 100.0, 3600.0, bio)
    assert (col.N + col.P + col.Z + col.D).sum() == pytest.approx(total0, rel=1e-6)
    lit, dark = _column(), _column()
    lit.N[:] = dark.N[:] = 5.0
    for _ in range(48):
        biology_step(lit, 100.0, 3600.0, bio)
        biology_step(dark, 0.0, 3600.0, bio)
    assert lit.P[:20].sum() > dark.P[:20].sum()


def test_vortex_peaks_at_the_radius_of_maximum_wind():
    r = np.array([5.0, 20.0, 40.0, 80.0, 300.0])
    w = vortex_wind(r, 50.0, 40.0)
    assert w.argmax() == 2 and w[-1] < w[2]
    t = np.arange(0, 48.0)
    lat = np.linspace(10, 14, t.size)
    taux, tauy, speed = storm_stress(t, lat, np.full(t.size, 88.0), np.full(t.size, 90.0), np.full(t.size, 40.0),
                                     12.0, 88.0)
    assert speed.max() > 30 and np.hypot(taux, tauy).max() > 2.0
