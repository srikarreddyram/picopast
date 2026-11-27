"""A 1-D ocean column for virtual cyclone experiments: PWP mixing with NPZD biology.

Physics: Price, Weller & Pinkel (1986). Each step, the wind stress accelerates the mixed layer, currents
rotate inertially, the column is made statically stable, the mixed layer deepens until its bulk Richardson
number reaches 0.65, and adjacent layers are mixed until the gradient Richardson number reaches 0.25.
Every tracer (temperature, salinity, nutrients, plankton) is mixed by the same operations.

Biology: a standard nitrogen-based NPZD model (nutrient, phytoplankton, zooplankton, detritus) with
temperature-dependent growth (Eppley 1972), Michaelis-Menten nutrient limitation, a saturating light curve,
Holling III grazing, linear mortality, remineralisation and detritus sinking. Light decays with depth at
Kd = 0.04 + 0.03 Chl. Chlorophyll = P x 1.59 mg Chl / mmol N (C:Chl 50, Redfield C:N 6.625).

Wind: a modified Rankine vortex from a best track (maximum wind, radius of maximum wind), with a 20 degree
inflow angle and half the translation velocity added on the moving storm's right (left in the south).
Stress: rho_air * Cd * U^2, Cd from Large & Pond (1981), capped at 2.4e-3 at high winds (Powell et al. 2003).

Units: depth m (positive down), time s, temperature degC, salinity psu, tracers mmol N m-3, light W m-2 PAR.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

G = 9.81
RHO0 = 1025.0
RHO_AIR = 1.22
ALPHA = 2.4e-4      # thermal expansion, 1/degC (warm tropical water)
BETA = 7.6e-4       # haline contraction, 1/psu
RB_CRIT = 0.65
RG_CRIT = 0.25
CHL_PER_N = 1.59    # mg Chl per mmol N
DAY = 86400.0


@dataclass
class Biology:
    mu_max_25: float = 1.0    # 1/day at 25 C
    q10_rate: float = 1.066   # Eppley: growth x 1.066 per degC
    k_n: float = 0.5          # half-saturation for nitrate, mmol N m-3
    i_k: float = 30.0         # light saturation, W m-2 PAR
    g_max: float = 1.0        # 1/day
    k_p: float = 1.0          # grazing half-saturation, mmol N m-3
    assim: float = 0.7        # zooplankton assimilation efficiency
    m_p: float = 0.05         # phytoplankton mortality, 1/day
    m_z: float = 0.05         # zooplankton mortality, 1/day
    remin: float = 0.05       # detritus remineralisation, 1/day
    w_sink: float = 5.0       # detritus sinking, m/day
    kd_water: float = 0.04    # 1/m
    kd_chl: float = 0.03      # m2 / mg Chl


@dataclass
class Column:
    z: np.ndarray                     # layer centres, m (positive down), uniform spacing
    T: np.ndarray
    S: np.ndarray
    N: np.ndarray
    P: np.ndarray
    Z: np.ndarray
    D: np.ndarray
    lat: float
    u: np.ndarray = field(default=None)
    v: np.ndarray = field(default=None)

    def __post_init__(self):
        n = self.z.size
        self.u = np.zeros(n) if self.u is None else self.u
        self.v = np.zeros(n) if self.v is None else self.v
        self.dz = float(self.z[1] - self.z[0])

    def copy(self) -> "Column":
        return Column(self.z.copy(), self.T.copy(), self.S.copy(), self.N.copy(), self.P.copy(), self.Z.copy(),
                      self.D.copy(), self.lat, self.u.copy(), self.v.copy())

    @property
    def tracers(self):
        return (self.T, self.S, self.N, self.P, self.Z, self.D, self.u, self.v)

    def density(self) -> np.ndarray:
        return RHO0 * (1 - ALPHA * (self.T - 25.0) + BETA * (self.S - 35.0))

    def chl(self) -> np.ndarray:
        return self.P * CHL_PER_N

    def mld(self, dT: float = 0.2) -> float:
        """Depth where temperature first falls dT below the top layer."""
        below = np.flatnonzero(self.T < self.T[0] - dT)
        return float(self.z[below[0]]) if below.size else float(self.z[-1])

    def surface_chl(self) -> float:
        """Mean chlorophyll over the first optical depth (what a satellite sees)."""
        kd = 0.04 + 0.03 * self.chl()
        od = np.cumsum(kd) * self.dz
        top = od <= 1.0
        top[0] = True
        return float(self.chl()[top].mean())


# ---------------------------------------------------------------- PWP mixing

def _mix(col: Column, j: int) -> None:
    """Homogenise layers 0..j (inclusive) for every tracer."""
    for a in col.tracers:
        a[: j + 1] = a[: j + 1].mean()


def _static(col: Column) -> int:
    """Mix down until the column below the mixed layer is denser; return the mixed-layer bottom index."""
    rho = col.density()
    j = 0
    while j + 1 < rho.size and rho[j + 1] <= rho[: j + 1].mean() + 1e-9:
        j += 1
    if j > 0:
        _mix(col, j)
    return j


def _bulk(col: Column, j: int) -> int:
    """Deepen the mixed layer while its bulk Richardson number is below 0.65."""
    while j + 1 < col.z.size:
        rho = col.density()
        h = (j + 1) * col.dz
        drho = rho[j + 1] - rho[j]
        dv2 = (col.u[j] - col.u[j + 1]) ** 2 + (col.v[j] - col.v[j + 1]) ** 2
        if dv2 < 1e-12 or G * drho * h / (RHO0 * dv2) >= RB_CRIT:
            break
        j += 1
        _mix(col, j)
    return j


def _gradient(col: Column, j: int, max_iter: int = 200) -> None:
    """Partially mix adjacent layers below the mixed layer while the gradient Richardson number < 0.25."""
    dz = col.dz
    for _ in range(max_iter):
        rho = col.density()
        dv2 = np.diff(col.u) ** 2 + np.diff(col.v) ** 2
        rg = np.where(dv2 > 1e-12, G * np.diff(rho) * dz / (RHO0 * np.maximum(dv2, 1e-12)), np.inf)
        rg[:j] = np.inf
        k = int(np.argmin(rg))
        if rg[k] >= RG_CRIT:
            return
        # PWP's partial mixing toward Rg = 0.25. Capped at 1 (the pair fully mixed): below a density inversion
        # Rg is negative and the uncapped fraction overshoots, which ran away in ~4% of real-storm runs.
        frac = min(1.0, 0.5 * (1 - rg[k] / RG_CRIT))
        for a in col.tracers:
            d = a[k] - a[k + 1]
            a[k] -= frac * d / 2
            a[k + 1] += frac * d / 2


def physics_step(col: Column, taux: float, tauy: float, dt: float) -> None:
    """One PWP step: wind stress into the mixed layer, inertial rotation, then the three mixing rules."""
    f = 2 * 7.2921e-5 * np.sin(np.radians(col.lat))
    j = _static(col)
    h = (j + 1) * col.dz
    col.u[: j + 1] += taux / (RHO0 * h) * dt
    col.v[: j + 1] += tauy / (RHO0 * h) * dt
    ang = -f * dt  # inertial rotation of the whole velocity profile
    c, s = np.cos(ang), np.sin(ang)
    u = col.u * c - col.v * s
    col.v = col.u * s + col.v * c
    col.u = u
    j = _static(col)
    j = _bulk(col, j)
    _gradient(col, j)


# ---------------------------------------------------------------- NPZD biology

def biology_step(col: Column, par0: float, dt: float, bio: Biology | None = None) -> None:
    """One explicit step of the NPZD model (rates per day, dt in seconds)."""
    bio = bio or Biology()
    days = dt / DAY
    chl = col.chl()
    kd = bio.kd_water + bio.kd_chl * chl
    light = par0 * np.exp(-np.cumsum(kd) * col.dz + 0.5 * kd * col.dz)
    mu = bio.mu_max_25 * bio.q10_rate ** (col.T - 25.0)
    N, P, Z, D = np.maximum(col.N, 0), np.maximum(col.P, 0), np.maximum(col.Z, 0), np.maximum(col.D, 0)
    growth = mu * N / (bio.k_n + N) * light / (bio.i_k + light) * P
    graze = bio.g_max * P ** 2 / (bio.k_p ** 2 + P ** 2) * Z
    dN = -growth + bio.remin * D
    dP = growth - graze - bio.m_p * P
    dZ = bio.assim * graze - bio.m_z * Z
    dD = (1 - bio.assim) * graze + bio.m_p * P + bio.m_z * Z - bio.remin * D
    col.N, col.P, col.Z, col.D = N + dN * days, P + dP * days, Z + dZ * days, D + dD * days
    # detritus sinking (upwind), out of the bottom
    flux = bio.w_sink * days / col.dz * col.D
    col.D -= flux
    col.D[1:] += flux[:-1]
    np.maximum(col.N, 0, out=col.N)


def step(col: Column, taux: float, tauy: float, par0: float, dt: float, bio: Biology | None = None) -> None:
    physics_step(col, taux, tauy, dt)
    biology_step(col, par0, dt, bio)


# ---------------------------------------------------------------- wind from a best track

def drag(u10: np.ndarray, cap: float = 2.4e-3) -> np.ndarray:
    """Large & Pond (1981) drag coefficient, capped (2.4e-3 by default) for hurricane winds (Powell et al. 2003)."""
    u10 = np.asarray(u10, float)
    cd = np.where(u10 < 11, 1.2e-3, (0.49 + 0.065 * u10) * 1e-3)
    return np.minimum(cd, cap)


def vortex_wind(dist_km: np.ndarray, vmax_ms: np.ndarray, rmw_km: np.ndarray, x: float = 0.6) -> np.ndarray:
    """Modified Rankine vortex: linear inside the radius of maximum wind, decaying as (Rm/r)^x outside."""
    r = np.maximum(np.asarray(dist_km, float), 1e-3)
    rm = np.asarray(rmw_km, float)
    return np.where(r < rm, vmax_ms * r / rm, vmax_ms * (rm / r) ** x)


def storm_stress(t_hours: np.ndarray, storm_lat: np.ndarray, storm_lon: np.ndarray, vmax_kt: np.ndarray,
                 rmw_km: np.ndarray, lat0: float, lon0: float, south: bool = False,
                 inflow_deg: float = 20.0, x: float = 0.6, cd_cap: float = 2.4e-3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Wind stress (taux, tauy, N m-2) and wind speed at (lat0, lon0) from an hourly track."""
    km_per_deg = 111.2
    dy = (lat0 - storm_lat) * km_per_deg
    dx = (lon0 - storm_lon) * km_per_deg * np.cos(np.radians(lat0))
    r = np.hypot(dx, dy)
    speed = vortex_wind(r, np.asarray(vmax_kt) * 0.5144, rmw_km, x)
    # tangential direction: counterclockwise in the north, clockwise in the south; plus inflow toward the centre
    sign = -1.0 if south else 1.0
    th = np.arctan2(dy, dx) + sign * (np.pi / 2 + np.radians(inflow_deg))
    u, v = speed * np.cos(th), speed * np.sin(th)
    # half the storm's translation velocity, added on the whole field (asymmetry)
    vy = np.gradient(storm_lat * km_per_deg * 1000, t_hours * 3600)
    vx = np.gradient(storm_lon * km_per_deg * 1000 * np.cos(np.radians(lat0)), t_hours * 3600)
    fade = np.exp(-r / 500.0)
    u, v = u + 0.5 * vx * fade, v + 0.5 * vy * fade
    w = np.hypot(u, v)
    cd = drag(w, cd_cap)
    return RHO_AIR * cd * w * u, RHO_AIR * cd * w * v, w


# ---------------------------------------------------------------- a whole experiment

def run(col: Column, taux: np.ndarray, tauy: np.ndarray, par0: float, dt: float,
        bio: Biology | None = None) -> dict[str, np.ndarray]:
    """Step through a forcing series; return daily-resolved surface temperature, surface chl and MLD."""
    sst, chl, mld, n0, ichl = [], [], [], [], []
    for tx, ty in zip(taux, tauy):
        step(col, float(tx), float(ty), par0, dt, bio)
        sst.append(col.T[0])
        chl.append(col.surface_chl())
        mld.append(col.mld())
        n0.append(col.N[0])
        ichl.append(col.chl().sum() * col.dz)
    return {"sst": np.asarray(sst), "chl": np.asarray(chl), "mld": np.asarray(mld), "n_surface": np.asarray(n0),
            "int_chl": np.asarray(ichl)}


def spin_up(col: Column, par0: float, days: float = 30.0, wind_ms: float = 5.0, dt: float = 3600.0,
            bio: Biology | None = None) -> Column:
    """Settle the biology under a steady light wind (no storm) so a deep chlorophyll maximum forms."""
    tau = RHO_AIR * drag(wind_ms) * wind_ms ** 2
    n = int(days * DAY / dt)
    for i in range(n):
        # a slowly rotating light wind, so the mixed layer does not accumulate momentum
        ang = 2 * np.pi * i * dt / DAY
        step(col, tau * np.cos(ang), tau * np.sin(ang), par0, dt, bio)
    col.u[:] = 0.0
    col.v[:] = 0.0
    return col
