"""Storms end heatwaves: the whole record, 1982-2026.

Rules: ``docs/cyclone_heatwave_preregistration.md``, section "Storms end heatwaves" (added 28 November 2026 and
committed before anything here was computed).

- W1 (primary): the daily rate at which heatwave cells end, storm-exposed against not, Mantel-Haenszel over
  basin x calendar month x 5-degree latitude band, with calendar years resampled.
- S1-S4 (Holm): the 1982-2001 survival replication, return time, next-month suppression, dose-response.
- D1-D6: descriptive.

Heatwave flags come from ``picopast.mhw`` over 1982-2026 with the fixed 1991-2020 baseline. Each box keeps one
uint8 flag array and one int16 SST-anomaly array (hundredths of a degree) under ``interim/heatwave_end``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from picopast import storms
from picopast.config import Config, load_config

START = pd.Timestamp("1982-01-01")
EARLY_LAST_YEAR = 2001
R_KM = 150.0
WINDOW = (-1, 1)                    # a fix on day d exposes the transition t -> t+1 when d - t lies in this range
STRONG_KT = 64.0
BAND_DEG = 5.0
CURVE_DAYS = 30
SURVIVAL_DAY = 5
RETURN_CAP = 240
OCC_WINDOWS = {"month 1": (6, 35), "months 2-3": (36, 95), "months 4-8": (96, 245)}
SESOI_DAYS = 10.0
N_BOOT = 5000
N_BOOT_YEARS = 2000
SEED = 0
SSTA_FILL = -32768
BASIN_CODE = {"Arabian Sea": 0, "Bay of Bengal": 1, "South Indian Ocean": 2}
SOUTH_CONFIG = "config/picopast_sio.yaml"
VARIANTS = {                         # W1 exposure variants: (wind column, radius km, window)
    "primary": ("wind", R_KM, WINDOW),
    "75 km": ("wind", 75.0, WINDOW),
    "300 km": ("wind", 300.0, WINDOW),
    "window t to t+1": ("wind", R_KM, (0, 1)),
    "window t-3 to t+1": ("wind", R_KM, (-3, 1)),
    "IMD wind": ("wind_imd", R_KM, WINDOW),   # north only
}


# ---------------------------------------------------------------- the daily heatwave cube, 1982-2026

def cube_dir(cfg: Config) -> Path:
    return cfg.path("interim") / "heatwave_end"


def build_cube(cfg: Config, rows: int = 10) -> Path:
    """Heatwave flags and SST anomaly for every day from 1 January 1982, from the box's OISST files."""
    import xarray as xr

    from picopast import mhw

    d = cube_dir(cfg)
    if (d / "grid.npz").exists():
        return d
    d.mkdir(parents=True, exist_ok=True)
    paths = sorted((cfg.path("raw") / "oisst").glob("oisst_v2.1_box_*.nc"))
    src = xr.open_mfdataset(paths, combine="by_coords")["sst"]
    raw_time = src["time"].values
    time = pd.DatetimeIndex(raw_time).normalize()
    if time[0] != START or not time.equals(pd.date_range(time[0], time[-1], freq="D")):
        raise ValueError("the OISST box files must run daily, without gaps, from 1 January 1982")
    ny, nx = src.sizes["lat"], src.sizes["lon"]
    hot = np.lib.format.open_memmap(d / "hot.npy", mode="w+", dtype=np.uint8, shape=(time.size, ny, nx))
    ssta = np.lib.format.open_memmap(d / "ssta.npy", mode="w+", dtype=np.int16, shape=(time.size, ny, nx))
    valid = np.zeros((ny, nx))
    doy = mhw.day_of_year_366(raw_time)
    for i in range(0, ny, rows):
        arr = src.isel(lat=slice(i, i + rows)).values.astype(np.float32)
        seas, thresh = mhw.climatology(arr, raw_time)
        hot[:, i:i + rows] = mhw.heatwave_days(arr, raw_time, thresh)
        a = arr - seas[doy - 1]
        ssta[:, i:i + rows] = np.where(np.isfinite(a), np.round(a * 100), SSTA_FILL).astype(np.int16)
        valid[i:i + rows] = np.isfinite(arr).mean(axis=0)
        del arr, seas, thresh, a
    hot.flush()
    ssta.flush()
    del hot, ssta
    np.savez(d / "grid.npz", lat=src["lat"].values, lon=src["lon"].values,
             time=time.values.astype("datetime64[D]"), ocean=valid > 0.5)
    return d


class Cube:
    """The box's daily heatwave flags (and SST anomaly) in memory, with the grid."""

    def __init__(self, cfg: Config, with_ssta: bool = True):
        d = cube_dir(cfg)
        g = np.load(d / "grid.npz")
        self.lat, self.lon, self.ocean = g["lat"], g["lon"], g["ocean"]
        self.time = pd.DatetimeIndex(g["time"].astype("datetime64[ns]"))
        self.T, self.nx = self.time.size, self.lon.size
        self.flat = np.load(d / "hot.npy").reshape(self.T, -1)
        self.ssta = np.load(d / "ssta.npy").reshape(self.T, -1) if with_ssta else None
        self.la2, self.lo2 = np.meshgrid(self.lat, self.lon, indexing="ij")

    def k(self, day) -> int:
        return int((pd.Timestamp(day).normalize() - self.time[0]).days)

    def cells(self, lat0: float, lon0: float, r_km: float) -> np.ndarray:
        iy, ix = storms.disc(self.lat, self.lon, lat0, lon0, r_km, self.ocean)
        return iy * self.nx + ix


def basin_cells(cfg: Config, cube: Cube) -> np.ndarray:
    """Basin code per cell (flattened), -1 on land, split as in the storm table."""
    b = (np.full(cube.la2.shape, BASIN_CODE["South Indian Ocean"]) if storms.is_south(cfg)
         else np.where(cube.lo2 < storms.SPLIT_LON, BASIN_CODE["Arabian Sea"], BASIN_CODE["Bay of Bengal"]))
    return np.where(cube.ocean, b, -1).ravel()


def lat_band(lat) -> np.ndarray:
    return np.floor(np.asarray(lat, float) / BAND_DEG).astype(int)


def strata_groups(cfg: Config, cube: Cube) -> tuple[np.ndarray, list[str]]:
    """Group id per cell (basin x 5-degree latitude band; -1 on land) and the group labels."""
    basin = basin_cells(cfg, cube)
    band = lat_band(cube.la2).ravel()
    code = np.where(basin >= 0, basin * 100 + band + 50, -1)
    uniq = np.unique(code[code >= 0])
    gid = np.full(code.size, -1)
    gid[code >= 0] = np.searchsorted(uniq, code[code >= 0])
    names = {v: k for k, v in BASIN_CODE.items()}
    labels = [f"{names[u // 100]} {(u % 100 - 50) * BAND_DEG:+.0f}°" for u in uniq]
    return gid, labels


# ---------------------------------------------------------------- W1: daily ending rate under storms

def daily_totals(cube: Cube, group: np.ndarray, n_groups: int, chunk: int = 365) -> np.ndarray:
    """(T-1, groups, 2): hot cell-days and endings (hot on t, not hot on t+1) per day and group."""
    idx = np.flatnonzero(group >= 0)
    G = np.zeros((idx.size, n_groups), np.float32)
    G[np.arange(idx.size), group[idx]] = 1
    out = np.zeros((cube.T - 1, n_groups, 2))
    for k in range(0, cube.T - 1, chunk):
        stop = min(k + chunk, cube.T - 1)
        h = cube.flat[k:stop + 1][:, idx].astype(np.float32)
        out[k:stop, :, 0] = h[:-1] @ G
        out[k:stop, :, 1] = (h[:-1] * (1 - h[1:])) @ G
    return out


def exposed_transitions(cube: Cube, fixes: pd.DataFrame, wind_col: str = "wind", r_km: float = R_KM,
                        window: tuple[int, int] = WINDOW) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unique storm-exposed transitions (t, cell) with the strongest wind among the fixes that expose them."""
    f = fixes[fixes[wind_col] >= storms.MIN_WIND_KT]
    days = ((f["time"].dt.normalize() - cube.time[0]).dt.days).to_numpy()
    ncell = cube.flat.shape[1]
    lo, hi = window
    keys, winds = [], []
    for la, lon, d, w in zip(f["lat"].to_numpy(), f["lon"].to_numpy(), days, f[wind_col].to_numpy()):
        ts = np.arange(d - hi, d - lo + 1)
        ts = ts[(ts >= 0) & (ts <= cube.T - 2)]
        if not ts.size:
            continue
        cells = cube.cells(la, lon, r_km)
        if not cells.size:
            continue
        keys.append((ts[:, None] * ncell + cells[None, :]).ravel())
        winds.append(np.full(ts.size * cells.size, w, np.float32))
    if not keys:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0, np.float32)
    k, w = np.concatenate(keys), np.concatenate(winds)
    order = np.argsort(k, kind="stable")
    k, w = k[order], w[order]
    first = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    return k[first] // ncell, k[first] % ncell, np.maximum.reduceat(w, first)


def count_array(cube: Cube, group: np.ndarray, n_groups: int, totals: np.ndarray, t: np.ndarray,
                cell: np.ndarray, wmax: np.ndarray) -> np.ndarray:
    """(years, 12, groups, 6): hot cell-days, endings, exposed hot cell-days, exposed endings, and the same two
    for exposures by a fix of 64 kt or more."""
    years = cube.time.year.to_numpy()[:-1] - START.year
    months = cube.time.month.to_numpy()[:-1] - 1
    C = np.zeros((years.max() + 1, 12, n_groups, 6))
    np.add.at(C[..., 0], (years, months), totals[:, :, 0])
    np.add.at(C[..., 1], (years, months), totals[:, :, 1])
    h0, h1 = cube.flat[t, cell], cube.flat[t + 1, cell]
    n1 = (h0 == 1).astype(float)
    a1 = ((h0 == 1) & (h1 == 0)).astype(float)
    strong = wmax >= STRONG_KT
    yi, mi, gi = years[t], months[t], group[cell]
    for col, v in ((2, n1), (3, a1), (4, n1 * strong), (5, a1 * strong)):
        np.add.at(C[..., col], (yi, mi, gi), v)
    return C


def _mh(n0, a0, n1, a1, a_all) -> dict:
    """Mantel-Haenszel rate ratio over the last axis (strata), storm-caused endings and attributable shares."""
    with np.errstate(divide="ignore", invalid="ignore"):
        N = n0 + n1
        ok = (n0 > 0) & (n1 > 0)
        num = np.where(ok, a1 * n0 / N, 0.0).sum(-1)
        den = np.where(ok, a0 * n1 / N, 0.0).sum(-1)
        caused = np.where(n0 > 0, a1 - n1 * a0 / n0, 0.0).sum(-1)
        return {"hr": num / den, "caused": caused, "share": caused / a_all.sum(-1),
                "share_exposed": caused / a1.sum(-1),
                "rate_exposed": a1.sum(-1) / n1.sum(-1), "rate_unexposed": a0.sum(-1) / n0.sum(-1),
                "exposed_share_of_days": n1.sum(-1) / (n0 + n1).sum(-1)}


def mh_all(S: np.ndarray) -> dict:
    """S: (..., strata, 6). The rate ratio for all exposures and per wind class, each against unexposed days."""
    n_all, a_all, n1, a1, n1s, a1s = (S[..., i] for i in range(6))
    n0, a0 = n_all - n1, a_all - a1
    return {"all": _mh(n0, a0, n1, a1, a_all), "64 kt and above": _mh(n0, a0, n1s, a1s, a_all),
            "34-63 kt": _mh(n0, a0, n1 - n1s, a1 - a1s, a_all)}


def year_weights(n_years: int, n_boot: int = N_BOOT_YEARS, seed: int = SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, n_years, (n_boot, n_years))
    return np.stack([np.bincount(r, minlength=n_years) for r in draws]).astype(float)


def w1_summary(C: np.ndarray, W: np.ndarray) -> dict:
    """Point values and year-resampled 95% intervals for a (years, ..., 6) count array."""
    flat = C.reshape(C.shape[0], -1, 6)
    point = mh_all(flat.sum(0))
    boot = mh_all(np.einsum("by,ysc->bsc", W, flat))
    out = {}
    for cls, vals in point.items():
        out[cls] = {}
        for k, v in vals.items():
            d = boot[cls][k]
            d = d[np.isfinite(d)]
            lo, hi = (np.percentile(d, [2.5, 97.5]) if d.size else (np.nan, np.nan))
            out[cls][k] = {"value": float(v), "lo": float(lo), "hi": float(hi)}
    out["endings"] = float(flat[..., 1].sum())
    out["hot_cell_days"] = float(flat[..., 0].sum())
    out["exposed_endings"] = float(flat[..., 3].sum())
    return out


def w1_verdict(hr: dict) -> str:
    if hr["lo"] > 1:
        return "storms make heatwaves end more often"
    if hr["hi"] < 1:
        return "storms make heatwaves end less often"
    return "no detectable difference"


# ---------------------------------------------------------------- track points: survival, return, occupancy

def cell_stats(flat: np.ndarray, k0: int, cells: np.ndarray) -> dict:
    """For a set of cells around day k0 (D0): survival curve of the D0-3 heatwave cells, their return time once
    ended by D0+5, and heatwave occupancy in each window for D0-3 heatwave and non-heatwave cells."""
    T = flat.shape[0]
    h0 = flat[k0 - storms.LAG, cells] == 1
    hot_c, cold_c = cells[h0], cells[~h0]
    out = {"n": int(cells.size), "n_hot": int(hot_c.size), "curve": np.full(CURVE_DAYS, np.nan),
           "ret": np.nan, "n_ended": 0}
    last = min(CURVE_DAYS, T - 1 - k0)
    if hot_c.size and last > 0:
        out["curve"][:last] = flat[k0 + 1:k0 + 1 + last][:, hot_c].mean(axis=1)
    if hot_c.size and k0 + RETURN_CAP < T:
        ended = hot_c[flat[k0 + SURVIVAL_DAY, hot_c] == 0]
        out["n_ended"] = int(ended.size)
        if ended.size:
            block = flat[k0 + SURVIVAL_DAY + 1:k0 + RETURN_CAP + 1][:, ended]
            first = block.argmax(axis=0) + SURVIVAL_DAY + 1
            out["ret"] = float(np.where(block.any(axis=0), first, RETURN_CAP).mean())
    for name, (a, b) in OCC_WINDOWS.items():
        hv = nv = np.nan
        if k0 + b < T:
            blk = flat[k0 + a:k0 + b + 1]
            hv = float(blk[:, hot_c].mean()) if hot_c.size else np.nan
            nv = float(blk[:, cold_c].mean()) if cold_c.size else np.nan
        out[f"occ_hot|{name}"], out[f"occ_not|{name}"] = hv, nv
    return out


def weighted(E: float, hot_v: float, not_v: float) -> float:
    """E * hot_v + (1 - E) * not_v, needing only the groups that carry weight."""
    if not np.isfinite(E) or (E > 0 and not np.isfinite(hot_v)) or (E < 1 and not np.isfinite(not_v)):
        return np.nan
    return (E * hot_v if E > 0 else 0.0) + ((1 - E) * not_v if E < 1 else 0.0)


def far_mask(cube: Cube, fixes: pd.DataFrame, d0: pd.Timestamp) -> np.ndarray:
    """Cells more than 750 km from every fix (any strength) between D0-10 and D0+10 (flattened)."""
    near = fixes[(fixes["time"] >= d0 - pd.Timedelta(days=10)) & (fixes["time"] <= d0 + pd.Timedelta(days=10))]
    far = np.ones(cube.la2.shape, bool)
    for la, lon in zip(near["lat"].to_numpy(), near["lon"].to_numpy()):
        iy, ix = storms.disc(cube.lat, cube.lon, la, lon, storms.CONTROL_KM)
        far[iy, ix] = False
    return far.ravel()


def rewarm(ssta: np.ndarray, k0: int, cells: np.ndarray) -> tuple[float, bool]:
    """The main study's re-warming days (disc SSTA back within 0.2 C of its D0-3 value; 60-day cap)."""
    T = ssta.shape[0]

    def mean(k):
        v = ssta[k, cells]
        v = v[v != SSTA_FILL]
        return float(v.mean()) / 100 if v.size else np.nan

    base = mean(k0 - storms.LAG)
    s = np.array([mean(k0 + j) if k0 + j < T else np.nan for j in range(storms.REWARM_MAX_DAYS + 1)])
    if not (np.isfinite(base) and np.isfinite(s[:6]).any()):
        return np.nan, False
    j_min = int(np.nanargmin(s[:6]))
    back = np.flatnonzero(s[j_min:] >= base - storms.REWARM_TOL)
    if back.size:
        return float(j_min + back[0]), True
    return float(min(storms.REWARM_MAX_DAYS, np.isfinite(s).sum() - 1)), False


def point_table(cfg: Config, cube: Cube, fixes: pd.DataFrame) -> pd.DataFrame:
    """One row per JTWC >= 34 kt track point, 1982 to the end of the record, with storm and control measures."""
    pts = storms.track_points(cfg, start=str(START.date()))
    pts = pts[pts["primary"]].reset_index(drop=True)
    basin = basin_cells(cfg, cube)
    band = lat_band(cube.la2).ravel()
    rows = []
    for (d0, bname), grp in pts.groupby(["d0", "basin"], sort=True):
        k0 = cube.k(d0)
        if k0 - storms.LAG < 0 or k0 + SURVIVAL_DAY >= cube.T:
            continue
        ctrl = np.flatnonzero(far_mask(cube, fixes, d0) & (basin == BASIN_CODE[bname]))
        cs = cell_stats(cube.flat, k0, ctrl)
        band_stats: dict = {}
        basin_share = float(cube.flat[k0 - storms.LAG, basin == BASIN_CODE[bname]].mean())
        for p in grp.itertuples():
            cells = cube.cells(p.lat, p.lon, R_KM)
            if not cells.size:
                continue
            pb = int(lat_band(p.lat))
            if pb not in band_stats:
                band_stats[pb] = cell_stats(cube.flat, k0, ctrl[band[ctrl] == pb])
            bs = band_stats[pb]
            st = cell_stats(cube.flat, k0, cells)
            E = st["n_hot"] / st["n"]
            v = cube.ssta[k0 - storms.LAG, cells]
            v = v[v != SSTA_FILL]
            rw_days, rw_event = rewarm(cube.ssta, k0, cells)
            row = {"sid": p.sid, "name": p.name, "time": p.time, "d0": d0, "year": d0.year, "lat": p.lat,
                   "lon": p.lon, "wind": p.wind, "speed_kmh": p.speed_kmh, "basin": bname, "season": p.season,
                   "south": int(bname == "South Indian Ocean"), "E": E, "basin_share": basin_share,
                   "ssta_pre": float(v.mean()) / 100 if v.size else np.nan, "n_cells": st["n"],
                   "rewarm_days": rw_days, "rewarm_event": rw_event,
                   "ret_storm": st["ret"], "ret_ctrl": cs["ret"], "ret_ctrl_band": bs["ret"],
                   "n_ended": st["n_ended"]}
            for j in range(CURVE_DAYS):
                row[f"cs{j + 1}"], row[f"cc{j + 1}"] = st["curve"][j], cs["curve"][j]
            row["surv_storm"], row["surv_ctrl"] = st["curve"][SURVIVAL_DAY - 1], cs["curve"][SURVIVAL_DAY - 1]
            row["surv_ctrl_band"] = bs["curve"][SURVIVAL_DAY - 1]
            for name in OCC_WINDOWS:
                row[f"occ_storm|{name}"] = weighted(E, st[f"occ_hot|{name}"], st[f"occ_not|{name}"])
                row[f"occ_ctrl|{name}"] = weighted(E, cs[f"occ_hot|{name}"], cs[f"occ_not|{name}"])
                row[f"occ_ctrl_band|{name}"] = weighted(E, bs[f"occ_hot|{name}"], bs[f"occ_not|{name}"])
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- checks against what the project already used

def check_flags(cfg: Config, cube: Cube) -> dict:
    """Rebuilt heatwave flags against the storm-table caches, on the days they share."""
    import xarray as xr

    out = {}
    for name in ("oisst_daily_mhw.nc", "oisst_daily_mhw_20011101_20120331.nc"):
        f = cfg.path("interim") / "storms" / name
        if not f.exists():
            continue
        with xr.open_dataset(f) as ds:
            same_grid = (np.allclose(ds["lat"].values, cube.lat) and np.allclose(ds["lon"].values, cube.lon))
            k = (pd.DatetimeIndex(ds["time"].values).normalize() - cube.time[0]).days.to_numpy()
            keep = (k >= 0) & (k < cube.T)
            old = ds["hot"].values[keep].reshape(int(keep.sum()), -1)
        new = cube.flat[k[keep]]
        out[name] = {"same_grid": bool(same_grid), "days": int(keep.sum()), "cell_days": int(old.size),
                     "mismatch": int((old != new).sum())}
    return out


def check_survival(cfg: Config, pts: pd.DataFrame) -> dict:
    """Rebuilt storm and control survival against the stored storm tables, point by point."""
    names = (["storm_table_all.parquet"] if storms.is_south(cfg)
             else ["storm_table.parquet", "storm_table_2002.parquet"])
    old = pd.concat([pd.read_parquet(cfg.path("interim") / "storms" / n,
                                     columns=["sid", "time", "survival_storm", "survival_control"])
                     for n in names if (cfg.path("interim") / "storms" / n).exists()], ignore_index=True)
    m = pts.merge(old, on=["sid", "time"], how="inner")
    out = {"matched_points": int(len(m))}
    for new, ref in (("surv_storm", "survival_storm"), ("surv_ctrl", "survival_control")):
        a, b = m[new].to_numpy(float), m[ref].to_numpy(float)
        same = (np.isnan(a) & np.isnan(b)) | (np.abs(a - b) < 1e-9)
        out[ref] = {"compared": int(len(m)), "mismatch": int((~same).sum())}
    return out


# ---------------------------------------------------------------- the point-level statistics

def pvalue(draws: np.ndarray) -> float:
    d = np.asarray(draws)
    d = d[np.isfinite(d)]
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return float(min(1.0, max(p, 1.0 / d.size)))


def holm(pvals: dict[str, float]) -> dict[str, float]:
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, running = len(items), {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        out[k] = running
    return out


def _storm_counts(sids: np.ndarray, n_boot: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Storm labels per row, and (n_boot, storms) draw counts (whole storms with replacement)."""
    labels, inv = np.unique(sids, return_inverse=True)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, labels.size, (n_boot, labels.size))
    W = np.zeros((n_boot, labels.size))
    np.add.at(W, (np.repeat(np.arange(n_boot), labels.size), draws.ravel()), 1)
    return inv, W


def boot_mean(v, sids, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    v, sids = np.asarray(v, float), np.asarray(sids)
    keep = np.isfinite(v)
    v, sids = v[keep], sids[keep]
    if not v.size:
        return {"effect": np.nan, "lo": np.nan, "hi": np.nan, "p": np.nan, "n_points": 0, "n_storms": 0}
    inv, W = _storm_counts(sids, n_boot, seed)
    sums = np.bincount(inv, weights=v)
    cnt = np.bincount(inv).astype(float)
    d = (W @ sums) / (W @ cnt)
    return {"effect": float(v.mean()), "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5)),
            "p": pvalue(d), "n_points": int(v.size), "n_storms": int(np.unique(sids).size)}


def boot_slope(y, X, sids, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """OLS of y on X (with an intercept as X[:, 0]); the coefficient of X[:, 1], storm-resampled."""
    y, X, sids = np.asarray(y, float), np.asarray(X, float), np.asarray(sids)
    keep = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X, sids = y[keep], X[keep], sids[keep]
    inv, W = _storm_counts(sids, n_boot, seed)
    n_s, k = W.shape[1], X.shape[1]
    XtX = np.zeros((n_s, k, k))
    Xty = np.zeros((n_s, k))
    np.add.at(XtX, inv, X[:, :, None] * X[:, None, :])
    np.add.at(Xty, inv, X * y[:, None])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    A, b = np.einsum("bs,sij->bij", W, XtX), np.einsum("bs,si->bi", W, Xty)
    d = np.array([np.linalg.lstsq(A[i], b[i], rcond=None)[0][1] for i in range(n_boot)])
    return {"slope": float(beta[1]), "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5)),
            "p": pvalue(d), "n_points": int(y.size), "n_storms": int(np.unique(sids).size)}


def s1_verdict(e: dict) -> str:
    return "replicated" if e["hi"] < 0 else "not replicated"


def s2_verdict(e: dict, sesoi: float = SESOI_DAYS) -> str:
    if e["lo"] > 0:
        return "storms keep heatwaves away longer" + (f" (by less than {sesoi:.0f} days)" if e["hi"] < sesoi else "")
    if e["hi"] < 0:
        return "heatwaves ended by storms come back sooner"
    if -sesoi < e["lo"] and e["hi"] < sesoi:
        return "no meaningful difference"
    return "inconclusive"


def s3_verdict(e: dict) -> str:
    if e["hi"] < 0:
        return "storms hold back heatwaves for the next month"
    if e["lo"] > 0:
        return "storms are followed by more heatwave days"
    return "no detectable difference"


def s4_verdict(e: dict) -> str:
    if e["hi"] < 0:
        return "stronger, slower storms end more of the heatwave"
    if e["lo"] > 0:
        return "the opposite: stronger, slower storms end less of it"
    return "no detectable dependence on storm power"


def mixing_power(d: pd.DataFrame) -> pd.Series:
    return np.log(d["wind"] ** 2 / d["speed_kmh"].where(d["speed_kmh"] > 0))


def point_tests(P: pd.DataFrame, n_boot: int = N_BOOT) -> dict:
    hw = P[P["E"] >= storms.HOT]
    res: dict = {"counts": {
        "points": int(len(P)), "storms": int(P["sid"].nunique()),
        "heatwave_points": int(len(hw)), "heatwave_storms": int(hw["sid"].nunique()),
        "heatwave_points_by_era": {"1982-2001": int((hw["year"] <= EARLY_LAST_YEAR).sum()),
                                   "2002-2026": int((hw["year"] > EARLY_LAST_YEAR).sum())}}}
    surv = hw["surv_storm"] - hw["surv_ctrl"]
    early = hw["year"] <= EARLY_LAST_YEAR
    s1 = boot_mean(surv[early], hw.loc[early, "sid"], n_boot)
    s1 |= {"storm_share": float(hw.loc[early, "surv_storm"].mean()),
           "control_share": float(hw.loc[early & surv.notna(), "surv_ctrl"].mean())}
    s2 = boot_mean(hw["ret_storm"] - hw["ret_ctrl"], hw["sid"], n_boot)
    ok2 = (hw["ret_storm"] - hw["ret_ctrl"]).notna()
    s2 |= {"storm_days": float(hw.loc[ok2, "ret_storm"].mean()), "control_days": float(hw.loc[ok2, "ret_ctrl"].mean())}
    occ = P["occ_storm|month 1"] - P["occ_ctrl|month 1"]
    s3 = boot_mean(occ, P["sid"], n_boot)
    ctrl_mean = float(P.loc[occ.notna(), "occ_ctrl|month 1"].mean())
    s3 |= {"storm_occupancy": float(P.loc[occ.notna(), "occ_storm|month 1"].mean()), "control_occupancy": ctrl_mean,
           "percent_of_control": float(100 * s3["effect"] / ctrl_mean)}
    m = mixing_power(hw)
    used = surv.notna() & np.isfinite(m)
    mz = (m - m[used].mean()) / m[used].std()
    X = np.column_stack([np.ones(len(hw)), mz, hw["south"]])
    s4 = boot_slope(surv, X, hw["sid"], n_boot)
    adj = holm({"S1": s1["p"], "S2": s2["p"], "S3": s3["p"], "S4": s4["p"]})
    for key, e, fn in (("S1", s1, s1_verdict), ("S2", s2, s2_verdict), ("S3", s3, s3_verdict), ("S4", s4, s4_verdict)):
        e["p_holm"] = adj[key]
        e["verdict"] = fn(e)
        e["holm_significant"] = bool(adj[key] < 0.05)
    res |= {"S1": s1, "S2": s2, "S3": s3, "S4": s4}
    res["survival_2002_2026"] = boot_mean(surv[~early], hw.loc[~early, "sid"], n_boot)
    res["survival_all"] = boot_mean(surv, hw["sid"], n_boot) | {
        "storm_share": float(hw.loc[surv.notna(), "surv_storm"].mean()),
        "control_share": float(hw.loc[surv.notna(), "surv_ctrl"].mean())}
    return res


def descriptive(P: pd.DataFrame, n_boot: int = N_BOOT) -> dict:
    from scipy.stats import theilslopes

    from picopast.storms_analysis import km_median

    hw = P[P["E"] >= storms.HOT]
    surv = hw["surv_storm"] - hw["surv_ctrl"]
    out: dict = {}
    ok = surv.notna()
    out["D1_curve"] = {"storm": [float(hw.loc[ok, f"cs{j}"].mean()) for j in range(1, CURVE_DAYS + 1)],
                       "control": [float(hw.loc[ok, f"cc{j}"].mean()) for j in range(1, CURVE_DAYS + 1)],
                       "difference": [boot_mean(hw[f"cs{j}"] - hw[f"cc{j}"], hw["sid"], n_boot) for j in (1, 5, 10, 20, 30)]}
    years = np.arange(START.year, int(P["year"].max()) + 1)
    enc = {}
    for reg, sub in (("north", P[P["south"] == 0]), ("south", P[P["south"] == 1])):
        h = sub[sub["E"] >= storms.HOT]
        n_st = h.groupby("year")["sid"].nunique().reindex(years, fill_value=0)
        n_pt = h.groupby("year").size().reindex(years, fill_value=0)
        share = (sub.assign(h=sub["E"] >= storms.HOT).groupby("year")["h"].mean()).reindex(years)
        full = years[years <= 2025]           # 2026 is incomplete
        ts = theilslopes(n_st.loc[full].to_numpy(float), full)
        sh = share.loc[full].dropna()
        tp = theilslopes(sh.to_numpy(float), sh.index.to_numpy(float))
        enc[reg] = {"storms_per_year": {int(k): int(v) for k, v in n_st.items()},
                    "points_per_year": {int(k): int(v) for k, v in n_pt.items()},
                    "share_of_points": {int(k): float(v) for k, v in share.items()},
                    "storms_per_decade": {"slope": 10 * ts.slope, "lo": 10 * ts.low_slope, "hi": 10 * ts.high_slope},
                    "share_per_decade": {"slope": 10 * tp.slope, "lo": 10 * tp.low_slope, "hi": 10 * tp.high_slope},
                    "by_decade": {int(k): int(v) for k, v in h.groupby(h["year"] // 10 * 10)["sid"].nunique().items()}}
        enc[reg]["chance"] = boot_mean(sub["E"] - sub["basin_share"], sub["sid"], n_boot) | {
            "mean_E": float(sub["E"].mean()), "mean_basin_share": float(sub["basin_share"].mean())}
    out["D2_encounters"] = enc
    groups = {"basin": hw["basin"], "season": hw["season"],
              "era": np.where(hw["year"] <= EARLY_LAST_YEAR, "1982-2001", "2002-2026"),
              "wind class": np.where(hw["wind"] >= STRONG_KT, "64 kt and above", "34-63 kt")}
    t = hw["ssta_pre"]
    groups["pre-storm SST anomaly"] = pd.qcut(t, 3, labels=["lowest third", "middle third", "highest third"]).astype(str)
    out["D3_modifiers"] = {name: {str(level): boot_mean(surv[g == level], hw.loc[g == level, "sid"], n_boot)
                                  for level in sorted(pd.unique(g)) if str(level) != "nan"}
                           for name, g in ((k, np.asarray(v)) for k, v in groups.items())}
    out["D3_modifiers"]["pre-storm SST anomaly cut points"] = [float(x) for x in t.quantile([1 / 3, 2 / 3])]
    rec = {}
    for reg, sub in (("pooled", P), ("north", P[P["south"] == 0]), ("south", P[P["south"] == 1])):
        rec[reg] = {}
        for cls, s in (("heatwave", sub[sub["E"] >= storms.HOT]), ("normal", sub[sub["E"] < storms.NORMAL])):
            s = s.dropna(subset=["rewarm_days"])
            rec[reg][cls] = {"median_days": km_median(s["rewarm_days"].to_numpy(), s["rewarm_event"].to_numpy(bool)),
                             "points": int(len(s)), "storms": int(s["sid"].nunique())}
    out["D4_recovery"] = rec
    out["D5_longer"] = {}
    for name in ("months 2-3", "months 4-8"):
        occ = P[f"occ_storm|{name}"] - P[f"occ_ctrl|{name}"]
        e = boot_mean(occ, P["sid"], n_boot)
        cm = float(P.loc[occ.notna(), f"occ_ctrl|{name}"].mean())
        out["D5_longer"][name] = e | {"control_occupancy": cm, "percent_of_control": float(100 * e["effect"] / cm)}
    out["D5_longer"]["heatwave points, month 1"] = boot_mean(
        hw["occ_storm|month 1"] - hw["occ_ctrl|month 1"], hw["sid"], n_boot)
    nm = P[P["E"] < storms.NORMAL]
    out["D5_longer"]["normal points, month 1"] = boot_mean(
        nm["occ_storm|month 1"] - nm["occ_ctrl|month 1"], nm["sid"], n_boot)
    early = hw["year"] <= EARLY_LAST_YEAR
    out["D6_band_control"] = {
        "S1": boot_mean((hw["surv_storm"] - hw["surv_ctrl_band"])[early], hw.loc[early, "sid"], n_boot),
        "S2": boot_mean(hw["ret_storm"] - hw["ret_ctrl_band"], hw["sid"], n_boot),
        "S3": boot_mean(P["occ_storm|month 1"] - P["occ_ctrl_band|month 1"], P["sid"], n_boot),
        "survival 1982-2026": boot_mean(hw["surv_storm"] - hw["surv_ctrl_band"], hw["sid"], n_boot)}
    return out


# ---------------------------------------------------------------- running it

def region_inputs(cfg: Config, name: str, refresh: bool = False) -> tuple[dict, pd.DataFrame, dict]:
    """Count arrays for every W1 variant, the point table and the checks for one box (cached)."""
    d = cube_dir(cfg)
    counts_path, pts_path, checks_path = d / "counts.npz", d / "points.parquet", d / "checks.json"
    if not refresh and counts_path.exists() and pts_path.exists() and checks_path.exists():
        z = np.load(counts_path)
        return {k: z[k] for k in z.files}, pd.read_parquet(pts_path), json.loads(checks_path.read_text())
    build_cube(cfg)
    cube = Cube(cfg)
    fixes = storms.all_fixes(cfg)
    checks = {"flags": check_flags(cfg, cube)}
    group, labels = strata_groups(cfg, cube)
    totals = daily_totals(cube, group, len(labels))
    counts = {}
    for vname, (col, r, win) in VARIANTS.items():
        if col == "wind_imd" and storms.is_south(cfg):
            continue
        t, cell, w = exposed_transitions(cube, fixes, col, r, win)
        counts[vname] = count_array(cube, group, len(labels), totals, t, cell, w)
        print(f"{name}: W1 variant {vname}: {t.size:,} exposed transitions", flush=True)
    np.savez(counts_path, **counts)
    pts = point_table(cfg, cube, fixes)
    pts.to_parquet(pts_path)
    checks["survival"] = check_survival(cfg, pts)
    checks["strata"] = labels
    checks["record"] = [str(cube.time[0].date()), str(cube.time[-1].date())]
    checks_path.write_text(json.dumps(checks, indent=1))
    print(f"{name}: {len(pts):,} points; checks {checks}", flush=True)
    return counts, pts, checks


def analyse(n_boot: int = N_BOOT, n_boot_years: int = N_BOOT_YEARS, refresh: bool = False) -> dict:
    cfg_n, cfg_s = load_config(), load_config(SOUTH_CONFIG)
    Cn, Pn, chk_n = region_inputs(cfg_n, "north", refresh)
    Cs, Ps, chk_s = region_inputs(cfg_s, "south", refresh)
    ny = Cn["primary"].shape[0]
    if Cs["primary"].shape[0] != ny:
        raise ValueError("the two boxes cover different years")
    W = year_weights(ny, n_boot_years)
    years = np.arange(START.year, START.year + ny)
    res: dict = {"checks": {"north": chk_n, "south": chk_s}, "years": [int(years[0]), int(years[-1])]}

    def three(cn, cs):
        return {"pooled": w1_summary(np.concatenate([cn, cs], axis=2), W), "north": w1_summary(cn, W),
                "south": w1_summary(cs, W)}

    w1 = {"primary": three(Cn["primary"], Cs["primary"])}
    for reg in ("pooled", "north", "south"):
        w1["primary"][reg]["verdict"] = w1_verdict(w1["primary"][reg]["all"]["hr"])
    for v in VARIANTS:
        if v == "primary":
            continue
        if v == "IMD wind":
            w1[v] = {"north": w1_summary(Cn[v], W)}
        else:
            w1[v] = three(Cn[v], Cs[v])
    for era, sel in (("1982-2001", years <= EARLY_LAST_YEAR), ("2002-2026", years > EARLY_LAST_YEAR)):
        We = year_weights(int(sel.sum()), n_boot_years)
        cn, cs = Cn["primary"][sel], Cs["primary"][sel]
        w1[f"era {era}"] = {"pooled": w1_summary(np.concatenate([cn, cs], axis=2), We),
                            "north": w1_summary(cn, We), "south": w1_summary(cs, We)}
    res["W1"] = w1
    P = pd.concat([Pn, Ps], ignore_index=True)
    res["points"] = point_tests(P, n_boot)
    res["descriptive"] = descriptive(P, n_boot)
    return res


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


def main(refresh: bool = False, n_boot: int = N_BOOT, n_boot_years: int = N_BOOT_YEARS) -> Path:
    from picopast import heatwave_end_report

    cfg_n = load_config()
    res = analyse(n_boot, n_boot_years, refresh)
    out = cube_dir(cfg_n) / "results.json"
    out.write_text(json.dumps(res, indent=1, default=_json_default))
    return heatwave_end_report.report(cfg_n, res)
