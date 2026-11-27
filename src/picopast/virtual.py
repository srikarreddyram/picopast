"""Virtual experiments: each storm run with and without its heatwave in a 1-D column model.

Rules: ``docs/cyclone_heatwave_preregistration.md``, section "Virtual experiments" (binding from 2026-11-27,
committed before any real storm was simulated). Model: ``picopast.column_model``.
"""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast import column_model as cm
from picopast.config import Config, load_config

DZ = 2.0
ZMAX = 200.0
PAR0 = 100.0
SPIN_DAYS = 30.0
STORM_HALF_DAYS = 2.0
POST_DAYS = 17.0
DT_STORM = 900.0
DT_CALM = 3600.0
BOX = (-40.0, 30.0, 25.0, 125.0)  # lat_min, lat_max, lon_min, lon_max of the profile cache
VARIANTS = {"base": {}, "cap20": {"cap_m": 20.0}, "par60": {"par0": 60.0}, "cd20": {"cd_cap": 2.0e-3},
            "x05": {"x": 0.5}}


# ---------------------------------------------------------------- climatological profiles

def woa_cache(cfg: Config) -> Path:
    """Monthly 0-200 m temperature, salinity and nitrate over the Indian Ocean, cached as one npz."""
    out = cfg.path("interim") / "virtual" / "woa23_indian_ocean.npz"
    if out.exists():
        return out
    raw = cfg.path("raw")
    T, S, N = [], [], []
    for m in range(1, 13):
        with xr.open_dataset(raw / "woa23_ts" / f"woa23_decav91C0_t{m:02d}_01.nc", decode_times=False) as t, \
             xr.open_dataset(raw / "woa23_ts" / f"woa23_decav91C0_s{m:02d}_01.nc", decode_times=False) as s, \
             xr.open_dataset(raw / "woa23_nitrate" / f"woa23_all_n{m:02d}_01.nc", decode_times=False) as n:
            sel = dict(lat=slice(BOX[0], BOX[1]), lon=slice(BOX[2], BOX[3]), depth=slice(0, ZMAX))
            T.append(t["t_an"].isel(time=0).sel(**sel).values)
            S.append(s["s_an"].isel(time=0).sel(**sel).values)
            N.append(n["n_an"].isel(time=0).sel(**sel).interp(depth=t["depth"].sel(depth=slice(0, ZMAX))).values)
            lat, lon = t["lat"].sel(lat=slice(BOX[0], BOX[1])).values, t["lon"].sel(lon=slice(BOX[2], BOX[3])).values
            depth = t["depth"].sel(depth=slice(0, ZMAX)).values
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, T=np.stack(T), S=np.stack(S), N=np.stack(N), lat=lat, lon=lon, depth=depth)
    return out


@dataclass
class Profiles:
    T: np.ndarray   # (month, depth, lat, lon)
    S: np.ndarray
    N: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    depth: np.ndarray

    @classmethod
    def load(cls, path: Path) -> "Profiles":
        z = np.load(path)
        return cls(z["T"], z["S"], z["N"], z["lat"], z["lon"], z["depth"])

    def at(self, month: int, lat: float, lon: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """(z, T, S, N) on the 2 m grid at the nearest cell with a surface value (searched outwards)."""
        m = month - 1
        iy, ix = int(np.abs(self.lat - lat).argmin()), int(np.abs(self.lon - lon).argmin())
        for rad in range(0, 4):
            ys = range(max(iy - rad, 0), min(iy + rad + 1, self.lat.size))
            xs = range(max(ix - rad, 0), min(ix + rad + 1, self.lon.size))
            cands = [(abs(y - iy) + abs(x - ix), y, x) for y in ys for x in xs
                     if np.isfinite(self.T[m, 0, y, x]) and np.isfinite(self.N[m, 0, y, x])]
            if cands:
                _, y, x = min(cands)
                break
        else:
            raise ValueError("no ocean profile near the point")
        z = np.arange(DZ / 2, ZMAX, DZ)

        def prof(a):
            v = a[m, :, y, x]
            ok = np.isfinite(v)
            return np.interp(z, self.depth[ok], v[ok])

        return z, prof(self.T), prof(self.S), np.maximum(prof(self.N), 0.0)


def shift_surface(z: np.ndarray, T: np.ndarray, dT: float, cap_m: float | None = None) -> np.ndarray:
    """Add dT over the surface layer (to the first 0.5 C drop, or to cap_m), tapering to 0 over the next 10 m."""
    below = np.flatnonzero(T < T[0] - 0.5)
    h = cap_m if cap_m is not None else (float(z[below[0]]) if below.size else float(z[-1]))
    w = np.clip(1 - (z - h) / 10.0, 0, 1)
    return T + dT * w


# ---------------------------------------------------------------- one run

def track_hourly(fixes: pd.DataFrame, sid: str, d0: pd.Timestamp) -> dict | None:
    t0 = d0 - pd.Timedelta(days=STORM_HALF_DAYS)
    t1 = d0 + pd.Timedelta(days=STORM_HALF_DAYS)
    f = fixes[(fixes["sid"] == sid)].sort_values("time")
    if f.empty:
        return None
    hours = pd.date_range(t0, t1, freq="1h")
    ft = (f["time"] - t0).dt.total_seconds().to_numpy() / 3600
    th = (hours - t0).total_seconds().to_numpy() / 3600
    wind = f["wind"].fillna(f["wind_imd"] * 1.11).to_numpy(float)
    rmw = f["rmw_km"].fillna(40.0).to_numpy(float)
    inside = (th >= ft.min()) & (th <= ft.max())
    out = {"t": th, "lat": np.interp(th, ft, f["lat"].to_numpy(float)),
           "lon": np.interp(th, ft, f["lon"].to_numpy(float)),
           "vmax": np.where(inside, np.interp(th, ft, np.nan_to_num(wind, nan=0.0)), 0.0),
           "rmw": np.interp(th, ft, rmw)}
    return out


def simulate(task: dict) -> dict:
    """One run: build the column, spin it up, apply the storm and 17 calm days; return the outputs."""
    prof = Profiles.load(Path(task["woa"]))
    z, T, S, N = prof.at(task["month"], task["lat"], task["lon"])
    T = shift_surface(z, T, task["dT"], task.get("cap_m"))
    P = np.where(z < 100, 0.1, 0.01)
    col = cm.Column(z, T, S.copy(), N.copy(), P, np.full(z.size, 0.05), np.full(z.size, 0.05), task["lat"])
    par0 = task.get("par0", PAR0)
    cm.spin_up(col, par0, SPIN_DAYS, 5.0, DT_CALM)
    tr = task["track"]
    taux_h, tauy_h, _ = cm.storm_stress(tr["t"], tr["lat"], tr["lon"], tr["vmax"], tr["rmw"], task["lat"],
                                        task["lon"], south=task["lat"] < 0, x=task.get("x", 0.6),
                                        cd_cap=task.get("cd_cap", 2.4e-3))
    per_hour = int(3600 / DT_STORM)
    taux, tauy = np.repeat(taux_h, per_hour), np.repeat(tauy_h, per_hour)
    sst_ref, chl_ref = col.T[0], col.surface_chl()
    int_ref = float(col.chl().sum() * col.dz)
    n_ref = col.N[0]
    storm = cm.run(col, taux, tauy, par0, DT_STORM)
    t_storm = (np.arange(taux.size) + 1) * DT_STORM / cm.DAY - STORM_HALF_DAYS   # days relative to D0
    tau = cm.RHO_AIR * cm.drag(5.0) * 25.0
    n_post = int(POST_DAYS * cm.DAY / DT_CALM)
    ang = 2 * np.pi * np.arange(n_post) * DT_CALM / cm.DAY
    post = cm.run(col, tau * np.cos(ang), tau * np.sin(ang), par0, DT_CALM)
    t_post = STORM_HALF_DAYS + (np.arange(n_post) + 1) * DT_CALM / cm.DAY
    t_all = np.concatenate([t_storm, t_post])
    cat = {k: np.concatenate([storm[k], post[k]]) for k in storm}
    w_cool = (t_all >= 0) & (t_all <= 5)
    w_bloom = (t_all >= 1) & (t_all <= 16)
    k16 = int(np.abs(t_all - 16.0).argmin())
    return {"key": task["key"], "version": task["version"],
            "cooling": float(cat["sst"][w_cool].min() - sst_ref),
            "bloom": float(np.log10(cat["chl"][w_bloom].mean()) - np.log10(chl_ref)),
            "int_chl_change": float(cat["int_chl"][k16] - int_ref),           # mg Chl m-2, D0+16 minus D0-2
            "mld_max": float(cat["mld"][t_all <= STORM_HALF_DAYS].max()),
            "n_surface_rise": float(cat["n_surface"][w_cool].max() - n_ref),  # max over D0..D0+5 minus D0-2
            "chl_ref": float(chl_ref), "int_chl_ref": int_ref}


# ---------------------------------------------------------------- the experiment

def select_points(d: pd.DataFrame) -> pd.DataFrame:
    """Per storm and class, the point with the highest JTWC wind (ties: the earliest)."""
    d = d[d["heatwave"].isin([0.0, 1.0])].sort_values(["sid", "heatwave", "wind", "time"],
                                                      ascending=[True, True, False, True])
    return d.groupby(["sid", "heatwave"], as_index=False).first()


def tasks_for(d: pd.DataFrame, fixes: pd.DataFrame, woa: Path, variant: str) -> list[dict]:
    opts = VARIANTS[variant]
    med = d[d["heatwave"] == 1].groupby("south")["ssta_pre"].median().to_dict()
    tasks = []
    prof = Profiles.load(woa)
    for p in d.itertuples():
        tr = track_hourly(fixes, p.sid, p.d0)
        if tr is None or not np.isfinite(p.sst_pre) or not np.isfinite(p.ssta_pre):
            continue
        try:
            _, T, _, _ = prof.at(p.d0.month, p.lat, p.lon)
        except ValueError:
            continue
        d_obs = p.sst_pre - T[0]
        if p.heatwave == 1:
            versions = {"actual": d_obs, "counterfactual": d_obs - p.ssta_pre}          # heatwave removed
        else:
            versions = {"actual": d_obs, "counterfactual": d_obs + med.get(p.south, 0.0)}  # heatwave added
        key = f"{p.sid}|{int(p.heatwave)}"
        for v, dT in versions.items():
            tasks.append({"key": key, "version": v, "woa": str(woa), "month": int(p.d0.month), "lat": float(p.lat),
                          "lon": float(p.lon), "dT": float(dT), "track": tr} | opts)
    return tasks


def experiment(cfg_n: Config, cfg_s: Config | None, variant: str = "base", workers: int = 6) -> pd.DataFrame:
    from picopast import final_test, storms

    out = cfg_n.path("interim") / "virtual" / f"runs_{variant}.parquet"
    if out.exists():
        return pd.read_parquet(out)
    d = select_points(final_test.pooled_table(cfg_n, cfg_s))
    fixes = storms.all_fixes(cfg_n, ("NI", "SI"))
    tasks = tasks_for(d, fixes, woa_cache(cfg_n), variant)
    with ProcessPoolExecutor(workers) as ex:
        runs = pd.DataFrame(list(ex.map(simulate, tasks, chunksize=4)))
    meta = d.assign(key=d["sid"] + "|" + d["heatwave"].astype(int).astype(str))
    runs = runs.merge(meta[["key", "sid", "heatwave", "south", "bloom", "cooling"]].rename(
        columns={"bloom": "obs_bloom", "cooling": "obs_cooling"}), on="key", how="left")
    out.parent.mkdir(parents=True, exist_ok=True)
    runs.to_parquet(out, index=False)
    return runs


def paired(runs: pd.DataFrame, cls: float, n_boot: int = 5000) -> dict:
    """Mean (with heatwave - without heatwave) per output, storm-resampled 95% intervals."""
    from picopast.storms_analysis import storm_draws

    r = runs[runs["heatwave"] == cls]
    w = r.pivot_table(index=["key", "sid"], columns="version",
                      values=["cooling", "bloom", "int_chl_change", "mld_max", "n_surface_rise"]).dropna()
    sign = 1.0 if cls == 1 else -1.0   # heatwave points: actual has it; normal points: counterfactual has it
    res = {}
    sids = w.index.get_level_values("sid").to_numpy()
    for v in ("bloom", "int_chl_change", "cooling", "mld_max", "n_surface_rise"):
        diff = sign * (w[(v, "actual")] - w[(v, "counterfactual")]).to_numpy()
        draws = [diff[i].mean() for i in storm_draws(sids, n_boot)]
        lo, hi = np.percentile(draws, [2.5, 97.5])
        res[v] = {"effect": float(diff.mean()), "lo": float(lo), "hi": float(hi), "n": int(diff.size)}
    for k in ("bloom", "int_chl_change"):
        b = res[k]
        b["verdict"] = ("in the model, a heatwave weakens it" if b["hi"] < 0 else
                        "in the model, a heatwave strengthens it" if b["lo"] > 0 else "no detectable model effect")
    return res


def validation(runs: pd.DataFrame) -> dict:
    from scipy.stats import spearmanr

    a = runs[runs["version"] == "actual"]
    out = {}
    for m, o in (("cooling", "obs_cooling"), ("bloom", "obs_bloom")):
        ok = a[[m, o]].dropna()
        out[m] = {"rho": float(spearmanr(ok[m], ok[o])[0]) if len(ok) > 5 else np.nan, "n": int(len(ok))}
    return out


def report(cfg_n: Config, cfg_s: Config | None, workers: int = 6) -> Path:
    res = {}
    for variant in VARIANTS:
        runs = experiment(cfg_n, cfg_s, variant, workers)
        res[variant] = {"heatwave_points": paired(runs, 1.0), "normal_points": paired(runs, 0.0)}
        if variant == "base":
            res["validation"] = validation(runs)
    (cfg_n.path("interim") / "virtual" / "results.json").write_text(json.dumps(res, indent=1, default=float))
    b = res["base"]["heatwave_points"]
    L = ["# Virtual experiments: each storm run with and without its heatwave", "",
         f"Generated {date.today().isoformat()} by `picopast virtual`. Rules: the \"Virtual experiments\" section of "
         "`docs/cyclone_heatwave_preregistration.md`, committed before any real storm was simulated. Model: 1-D PWP "
         "mixing with NPZD biology (`picopast.column_model`). Differences are (with heatwave − without heatwave); "
         "intervals resample storms.", "",
         f"**Heatwave points ({b['bloom']['n']} storms), heatwave removed. Surface bloom: {b['bloom']['verdict']}; "
         f"depth-integrated chlorophyll: {b['int_chl_change']['verdict']}.**", "",
         "| output | heatwave points: effect [95%] | normal points, heatwave added: effect [95%] |",
         "| --- | --- | --- |"]
    nrm = res["base"]["normal_points"]
    labels = {"bloom": "bloom (log10 surface chl)", "int_chl_change": "depth-integrated chl change (mg m⁻²)",
              "cooling": "cooling (°C)", "mld_max": "deepest mixed layer (m)",
              "n_surface_rise": "surface nitrate rise (mmol m⁻³)"}
    for k, lab in labels.items():
        h, n = b[k], nrm[k]
        L.append(f"| {lab} | {h['effect']:+.3f} [{h['lo']:+.3f}, {h['hi']:+.3f}] | "
                 f"{n['effect']:+.3f} [{n['lo']:+.3f}, {n['hi']:+.3f}] |")
    v = res["validation"]
    L += ["", f"Validation (actual runs vs observations, Spearman): cooling ρ {v['cooling']['rho']:+.2f} "
          f"(n {v['cooling']['n']}), bloom ρ {v['bloom']['rho']:+.2f} (n {v['bloom']['n']}).", "",
          "| sensitivity | bloom effect, heatwave points [95%] |", "| --- | --- |"]
    for variant in VARIANTS:
        e = res[variant]["heatwave_points"]["bloom"]
        L.append(f"| {variant} | {e['effect']:+.3f} [{e['lo']:+.3f}, {e['hi']:+.3f}] |")
    L += ["", "Limits: one column, no advection or eddies; the heatwave enters as a warmer surface layer; fixed "
          "literature biology.", ""]
    out = cfg_n.path("reports") / "virtual_experiments.md"
    out.write_text("\n".join(L))
    return out


def main(workers: int = 6) -> Path:
    return report(load_config(), load_config("config/picopast_sio.yaml"), workers)
