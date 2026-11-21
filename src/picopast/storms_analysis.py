"""When cyclones hit heatwaves: the pre-registered tests and the report.

Rules: ``docs/cyclone_heatwave_preregistration.md``. Every interval resamples whole storms
(2,000 draws, seed 0), so points from one storm never count as independent.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from picopast.config import Config
from picopast.storms import HOT, NORMAL, R_KM

N_BOOT = 2000
SEED = 0
MIN_STORMS = 15
COVS = ["wind", "speed_kmh", "lat", "chl_pre", "ssha_cm", "dist2land_km"]
HW_INPUTS = ["ssta_pre", "E", "excess_pre"]
ML_INPUTS = ["wind", "speed_kmh", "lat", "bay", "doy_sin", "doy_cos", "dist2land_km", "sst_pre", "ssta_pre", "E",
             "excess_pre", "chl_pre", "ssha_cm"]
LGBM = dict(n_estimators=200, num_leaves=15, learning_rate=0.05, min_child_samples=20, colsample_bytree=0.8,
            random_state=0, verbose=-1)


# ---------------------------------------------------------------- regression with storm resampling

def design(d: pd.DataFrame, exposure: str, covs: list[str]) -> tuple[np.ndarray, list[str]]:
    """Intercept, exposure, covariates, Bay of Bengal and season dummies (pre-monsoon is the reference)."""
    cols = {"intercept": np.ones(len(d)), exposure: d[exposure].to_numpy(float)}
    for c in covs:
        cols[c] = d[c].to_numpy(float)
    bay = (d["basin"] == "Bay of Bengal").to_numpy(float)
    if 0 < bay.sum() < len(d):
        cols["bay"] = bay
    for s in ("post-monsoon", "other"):
        x = (d["season"] == s).to_numpy(float)
        if 0 < x.sum() < len(d):
            cols[s] = x
    return np.column_stack(list(cols.values())), list(cols)


def ols(X: np.ndarray, y: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    if w is not None:
        s = np.sqrt(w)
        X, y = X * s[:, None], y * s
    return np.linalg.lstsq(X, y, rcond=None)[0]


def storm_draws(storms: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED):
    """Yield row indices of storm-resampled data sets (whole storms drawn with replacement)."""
    labels, inverse = np.unique(storms, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(labels.size)]
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        yield np.concatenate([members[i] for i in rng.integers(0, labels.size, labels.size)])


def effect(df: pd.DataFrame, response: str, exposure: str = "heatwave", covs=COVS, weights: str | None = None,
           n_boot: int = N_BOOT) -> dict:
    """The exposure coefficient with a storm-resampled 95% interval, plus the counts the decision needs."""
    need = [response, exposure] + list(covs)
    d = df.dropna(subset=need).reset_index(drop=True)
    out = {"response": response, "exposure": exposure, "n_points": int(len(d)), "n_storms": int(d["sid"].nunique())}
    if exposure == "heatwave":
        out["storms_heatwave"] = int(d.loc[d[exposure] == 1, "sid"].nunique())
        out["storms_normal"] = int(d.loc[d[exposure] == 0, "sid"].nunique())
        out["points_heatwave"] = int((d[exposure] == 1).sum())
        out["points_normal"] = int((d[exposure] == 0).sum())
    if len(d) < 20 or d[exposure].nunique() < 2:
        return out | {"effect": np.nan, "lo": np.nan, "hi": np.nan}
    X, names = design(d, exposure, list(covs))
    y = d[response].to_numpy(float)
    w = None
    if weights == "storm_equal":
        w = (1.0 / d.groupby("sid")["sid"].transform("size")).to_numpy(float)
    j = names.index(exposure)
    est = float(ols(X, y, w)[j])
    draws = []
    sids = d["sid"].to_numpy()
    for idx in storm_draws(sids, n_boot):
        Xb, yb = X[idx], y[idx]
        if np.ptp(Xb[:, j]) == 0:
            continue
        draws.append(ols(Xb, yb, None if w is None else w[idx])[j])
    lo, hi = np.percentile(draws, [2.5, 97.5]) if draws else (np.nan, np.nan)
    return out | {"effect": est, "lo": float(lo), "hi": float(hi)}


def raw_difference(df: pd.DataFrame, response: str, n_boot: int = N_BOOT) -> dict:
    d = df.dropna(subset=[response, "heatwave"]).reset_index(drop=True)
    y, h = d[response].to_numpy(float), d["heatwave"].to_numpy(float)

    def diff(i):
        a, b = y[i][h[i] == 1], y[i][h[i] == 0]
        return a.mean() - b.mean() if a.size and b.size else np.nan

    est = diff(np.arange(len(d)))
    draws = [diff(i) for i in storm_draws(d["sid"].to_numpy(), n_boot)]
    draws = np.asarray(draws)[np.isfinite(draws)]
    return {"diff": float(est), "lo": float(np.percentile(draws, 2.5)), "hi": float(np.percentile(draws, 97.5)),
            "mean_heatwave": float(y[h == 1].mean()), "mean_normal": float(y[h == 0].mean())}


def verdict_bloom(e: dict) -> str:
    if e.get("storms_heatwave", 0) < MIN_STORMS or e.get("storms_normal", 0) < MIN_STORMS:
        return "not enough data"
    if e["hi"] < 0:
        return "weaker bloom over heatwaves (hypothesis A)"
    if e["lo"] > 0:
        return "stronger bloom over heatwaves (hypothesis B)"
    return "no detectable difference"


def verdict_cooling(e: dict) -> str:
    if e.get("storms_heatwave", 0) < MIN_STORMS or e.get("storms_normal", 0) < MIN_STORMS:
        return "not enough data"
    if e["hi"] < 0:
        return "stronger cooling over heatwaves"
    if e["lo"] > 0:
        return "weaker cooling over heatwaves"
    return "no detectable difference"


def strata(df: pd.DataFrame, response: str) -> pd.DataFrame:
    d = df.dropna(subset=[response, "heatwave"])
    rows = []
    for (basin, seas), g in d.groupby(["basin", "season"]):
        a, b = g.loc[g["heatwave"] == 1, response], g.loc[g["heatwave"] == 0, response]
        rows.append({"basin": basin, "season": seas, "points_heatwave": len(a), "points_normal": len(b),
                     "storms_heatwave": g.loc[g["heatwave"] == 1, "sid"].nunique(),
                     "storms_normal": g.loc[g["heatwave"] == 0, "sid"].nunique(),
                     "mean_heatwave": a.mean() if len(a) else np.nan, "mean_normal": b.mean() if len(b) else np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- RQ1 extras

def survival_effect(df: pd.DataFrame, n_boot: int = N_BOOT) -> dict:
    d = df[(df["heatwave"] == 1)].dropna(subset=["survival_storm", "survival_control"]).reset_index(drop=True)
    v = (d["survival_storm"] - d["survival_control"]).to_numpy(float)
    if not v.size:
        return {"n_points": 0, "n_storms": 0, "effect": np.nan, "lo": np.nan, "hi": np.nan}
    draws = [v[i].mean() for i in storm_draws(d["sid"].to_numpy(), n_boot)]
    return {"n_points": int(v.size), "n_storms": int(d["sid"].nunique()), "effect": float(v.mean()),
            "lo": float(np.percentile(draws, 2.5)), "hi": float(np.percentile(draws, 97.5)),
            "storm_share": float(d["survival_storm"].mean()), "control_share": float(d["survival_control"].mean())}


def km_median(t: np.ndarray, event: np.ndarray) -> float:
    """Kaplan-Meier median time (NaN if the curve never reaches 0.5)."""
    t, event = np.asarray(t, float), np.asarray(event, bool)
    s = 1.0
    for u in np.unique(t[event]):
        at_risk = (t >= u).sum()
        s *= 1 - ((t == u) & event).sum() / at_risk
        if s <= 0.5:
            return float(u)
    return np.nan


# ---------------------------------------------------------------- RQ4 prediction

def _ridge_fit_predict(Xtr, ytr, Xte, alpha=1.0):
    med = np.nanmedian(Xtr, axis=0)
    Xtr, Xte = np.where(np.isnan(Xtr), med, Xtr), np.where(np.isnan(Xte), med, Xte)
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1
    A, B = (Xtr - mu) / sd, (Xte - mu) / sd
    ym = ytr.mean()
    coef = np.linalg.solve(A.T @ A + alpha * np.eye(A.shape[1]), A.T @ (ytr - ym))
    return B @ coef + ym


def loso(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Out-of-storm predictions: baseline, ridge, LightGBM with and without the heatwave inputs."""
    from lightgbm import LGBMRegressor

    d = df.dropna(subset=[target]).reset_index(drop=True).copy()
    d["bay"] = (d["basin"] == "Bay of Bengal").astype(float)
    no_hw = [c for c in ML_INPUTS if c not in HW_INPUTS]
    pred = {k: np.full(len(d), np.nan) for k in ("baseline", "ridge", "lgbm", "lgbm_no_heatwave")}
    for s in d["sid"].unique():
        te = (d["sid"] == s).to_numpy()
        tr = ~te
        ytr = d.loc[tr, target].to_numpy(float)
        means = d[tr].groupby(["basin", "season"])[target].mean()
        pred["baseline"][te] = [means.get((b, se), ytr.mean()) for b, se in zip(d.loc[te, "basin"], d.loc[te, "season"])]
        pred["ridge"][te] = _ridge_fit_predict(d.loc[tr, ML_INPUTS].to_numpy(float), ytr,
                                               d.loc[te, ML_INPUTS].to_numpy(float))
        for key, cols in (("lgbm", ML_INPUTS), ("lgbm_no_heatwave", no_hw)):
            m = LGBMRegressor(**LGBM).fit(d.loc[tr, cols], ytr)
            pred[key][te] = m.predict(d.loc[te, cols])
    out = d[["point", "sid", target]].copy()
    for k, v in pred.items():
        out[k] = v
    return out


def skill_table(p: pd.DataFrame, target: str, n_boot: int = N_BOOT) -> dict:
    from scipy.stats import spearmanr

    y = p[target].to_numpy(float)

    def skill(idx, key):
        base = np.mean((y[idx] - p["baseline"].to_numpy()[idx]) ** 2)
        return 1 - np.mean((y[idx] - p[key].to_numpy()[idx]) ** 2) / base

    allidx = np.arange(len(p))
    res = {k: {"skill": float(skill(allidx, k)), "rho": float(spearmanr(p[k], y)[0])}
           for k in ("baseline", "ridge", "lgbm", "lgbm_no_heatwave")}
    draws = [skill(i, "lgbm") - skill(i, "lgbm_no_heatwave") for i in storm_draws(p["sid"].to_numpy(), n_boot)]
    gain = skill(allidx, "lgbm") - skill(allidx, "lgbm_no_heatwave")
    res["heatwave_gain"] = {"gain": float(gain), "lo": float(np.percentile(draws, 2.5)),
                            "hi": float(np.percentile(draws, 97.5))}
    res["n_points"], res["n_storms"] = int(len(p)), int(p["sid"].nunique())
    return res


def shap_importance(df: pd.DataFrame, target: str) -> pd.Series:
    from lightgbm import LGBMRegressor

    d = df.dropna(subset=[target]).copy()
    d["bay"] = (d["basin"] == "Bay of Bengal").astype(float)
    m = LGBMRegressor(**LGBM).fit(d[ML_INPUTS], d[target])
    contrib = m.predict(d[ML_INPUTS], pred_contrib=True)[:, :-1]
    return pd.Series(np.abs(contrib).mean(0), index=ML_INPUTS).sort_values(ascending=False)


# ---------------------------------------------------------------- everything

def run(cfg: Config, table: pd.DataFrame) -> dict:
    t = table.copy()
    prim = t[t["primary"]].copy()
    res: dict = {"counts": {
        "points": int(len(prim)), "storms": int(prim["sid"].nunique()),
        "points_heatwave": int((prim["heatwave"] == 1).sum()), "points_normal": int((prim["heatwave"] == 0).sum()),
        "points_between": int(prim["heatwave"].isna().sum()),
        "storms_heatwave": int(prim.loc[prim["heatwave"] == 1, "sid"].nunique()),
        "storms_normal": int(prim.loc[prim["heatwave"] == 0, "sid"].nunique()),
        "bloom_valid": int(prim["bloom"].notna().sum()), "cooling_valid": int(prim["cooling"].notna().sum()),
    }}
    # RQ2 primary and sensitivities
    e = effect(prim, "bloom")
    res["rq2"] = e | {"verdict": verdict_bloom(e)}
    res["rq2_raw"] = raw_difference(prim, "bloom")
    sens = {"storm-equal weights": effect(prim, "bloom", weights="storm_equal"),
            "continuous E (0 to 1), all points": effect(prim, "bloom", exposure="E"),
            "larger post period": effect(prim, "bloom_max")}
    for r in (75, 300):
        s = prim.copy()
        s["heatwave"] = np.where(s[f"E_r{r}"] >= HOT, 1.0, np.where(s[f"E_r{r}"] < NORMAL, 0.0, np.nan))
        s["chl_pre"] = s[f"chl_pre_r{r}"]
        sens[f"{r} km disc"] = effect(s, f"bloom_r{r}")
    s = prim.copy()
    s["heatwave"] = np.where(s["E_l5"] >= HOT, 1.0, np.where(s["E_l5"] < NORMAL, 0.0, np.nan))
    sens["exposure on D0 - 5"] = effect(s, "bloom")
    s = t[t["imd"]].copy()
    s["wind"] = s["wind_imd"]
    sens["IMD winds"] = effect(s, "bloom")
    s = prim[prim["d0"] >= pd.Timestamp("2015-05-03")]
    sens["with SMAP salinity (from May 2015)"] = effect(s, "bloom", covs=COVS + ["sss"])
    res["rq2_sensitivity"] = sens
    res["rq2_strata"] = strata(prim, "bloom").to_dict("records")
    # RQ1
    c = effect(prim, "cooling")
    res["rq1_cooling"] = c | {"verdict": verdict_cooling(c)}
    res["rq1_cooling_raw"] = raw_difference(prim, "cooling")
    res["rq1_survival"] = survival_effect(prim)
    km = {}
    for label, cls in (("heatwave", 1), ("normal", 0)):
        g = prim[(prim["heatwave"] == cls) & prim["rewarm_days"].notna()]
        km[label] = {"n": int(len(g)), "median_days": km_median(g["rewarm_days"], g["rewarm_event"].astype(bool)),
                     "share_rewarmed": float(g["rewarm_event"].astype(bool).mean()) if len(g) else np.nan}
    res["rq1_rewarm"] = km
    res["rq1_excess_cooling_posthoc"] = excess_cooling(cfg, prim)
    # RQ3
    pace = prim.dropna(subset=["d_peuk", "d_pro"], how="all")
    rows = []
    for sid, g in pace.groupby("sid"):
        rows.append({"storm": g["name"].iloc[0], "sid": sid, "date": str(g["d0"].min().date()), "points": len(g),
                     "d_peuk": g["d_peuk"].median(), "d_syn": g["d_syn"].median(), "d_pro": g["d_pro"].median(),
                     "bloom": g["bloom"].median(), "heatwave_share": g["E"].mean()})
    res["rq3"] = rows
    res["rq3_extended"] = rq3_extended(prim)
    # RQ4
    res["rq4"] = {}
    for target in ("cooling", "bloom"):
        p = loso(prim, target)
        res["rq4"][target] = skill_table(p, target) | {"shap": shap_importance(prim, target).round(4).to_dict()}
    return res


GROUPS3 = {"peuk": "picoeukaryotes", "syn": "Synechococcus", "pro": "Prochlorococcus"}


def mean_ci(values, sids, n_boot: int = N_BOOT) -> dict:
    """Mean of per-point values with a storm-resampled 95% interval."""
    v, s = np.asarray(values, float), np.asarray(sids)
    ok = np.isfinite(v)
    v, s = v[ok], s[ok]
    if v.size < 3:
        return {"mean": np.nan, "lo": np.nan, "hi": np.nan, "n": int(v.size), "storms": int(np.unique(s).size)}
    draws = [v[i].mean() for i in storm_draws(s, n_boot)]
    return {"mean": float(v.mean()), "lo": float(np.percentile(draws, 2.5)), "hi": float(np.percentile(draws, 97.5)),
            "n": int(v.size), "storms": int(np.unique(s).size)}


def verdict_group(e: dict) -> str:
    if e.get("storms_heatwave", 0) < MIN_STORMS or e.get("storms_normal", 0) < MIN_STORMS:
        return "not enough data"
    if e["hi"] < 0:
        return "weaker response over heatwaves"
    if e["lo"] > 0:
        return "stronger response over heatwaves"
    return "no detectable difference"


def rq3_extended(prim: pd.DataFrame) -> dict:
    """RQ3 extended (pre-registered 22 November 2026): group changes after all storms from PicoPast's records."""
    out = {}
    for variant, suffix in (("primary", "_nosst"), ("secondary", "")):
        groups = {}
        for g in GROUPS3:
            col, pre = f"rec_{g}{suffix}", f"recpre_{g}{suffix}"
            if col not in prim:
                continue
            d = prim.dropna(subset=[col])
            allm = mean_ci(d[col], d["sid"])
            allm["verdict"] = ("rises after storms" if allm["lo"] > 0 else "falls after storms" if allm["hi"] < 0
                               else "no detectable change")
            e = effect(prim, col, covs=[pre if c == "chl_pre" else c for c in COVS])
            groups[g] = {"all": allm, "heatwave": mean_ci(d.loc[d["heatwave"] == 1, col], d.loc[d["heatwave"] == 1, "sid"]),
                         "normal": mean_ci(d.loc[d["heatwave"] == 0, col], d.loc[d["heatwave"] == 0, "sid"]),
                         "effect": e | {"verdict": verdict_group(e)}}
        shifts = {}
        for a in ("peuk", "syn"):
            ca, cp = f"rec_{a}{suffix}", f"rec_pro{suffix}"
            if ca in prim and cp in prim:
                d = prim.dropna(subset=[ca, cp])
                shifts[f"{a}-pro"] = mean_ci(d[ca] - d[cp], d["sid"])
        out[variant] = {"groups": groups, "shifts": shifts}
    if {"d_peuk", "rec_peuk_nosst"} <= set(prim.columns):
        both = prim.dropna(subset=["d_peuk", "rec_peuk_nosst"])
        per = both.groupby("sid").agg(name=("name", "first"), moana=("d_peuk", "median"),
                                      picopast=("rec_peuk_nosst", "median"))
        out["cross_check"] = {"storms": int(len(per)),
                              "same_sign": int((np.sign(per["moana"]) == np.sign(per["picopast"])).sum()),
                              "rows": per.reset_index().round(3).to_dict("records")}
    return out


def _fmt(e: dict, k="effect", digits=3) -> str:
    return f"{e[k]:+.{digits}f} [{e['lo']:+.{digits}f}, {e['hi']:+.{digits}f}]"


def _days(v) -> str:
    return "over 60 days" if not np.isfinite(v) else f"{v:.0f} days"


def excess_cooling(cfg: Config, prim: pd.DataFrame, n_boot: int = N_BOOT) -> dict:
    """POST HOC (not pre-registered): storm cooling beyond what far-away water of the same heatwave status did.

    Per cell, min(D0..D0+5) minus mean(D0-7..D0-3) of OISST; storm = mean over the 150 km disc's ocean cells;
    control = mean over cells in the same basin, > 750 km from every fix within +-10 days, with the same heatwave
    status on D0-3. Difference of the excess (heatwave minus normal points) with a storm-resampled interval.
    """
    from picopast import storms

    oc = storms.daily_ocean(cfg)
    time = pd.DatetimeIndex(oc["time"].values)
    day = {d: i for i, d in enumerate(time)}
    lat, lon = oc["lat"].values, oc["lon"].values
    sst = oc["sst"].values.astype(np.float32)
    hot = oc["hot"].values
    is_ocean = np.isfinite(sst).mean(0) > 0.5
    la2, lo2 = np.meshgrid(lat, lon, indexing="ij")
    fixes = storms.all_fixes(cfg)
    cache, rows = {}, []
    for r in prim.dropna(subset=["heatwave"]).itertuples():
        k0, ka = day.get(r.d0), day.get(r.d0 - pd.Timedelta(days=storms.LAG))
        if k0 is None or ka is None or k0 < 7 or k0 + 6 > time.size:
            continue
        key = (r.d0, r.basin)
        if key not in cache:
            near = fixes[(fixes["time"] >= r.d0 - pd.Timedelta(days=10)) & (fixes["time"] <= r.d0 + pd.Timedelta(days=10))]
            far = is_ocean & ((lo2 < storms.SPLIT_LON) if r.basin == "Arabian Sea" else (lo2 >= storms.SPLIT_LON))
            for f in near.itertuples():
                far &= storms.haversine_km(f.lat, f.lon, la2, lo2) > storms.CONTROL_KM
            pre = np.nanmean(sst[k0 - 7:k0 - 2], axis=0)
            cache[key] = (np.nanmin(sst[k0:k0 + 6], axis=0) - pre, far)
        change, far = cache[key]
        iy, ix = storms.disc(lat, lon, r.lat, r.lon, R_KM, is_ocean)
        ctrl = change[far & (hot[ka] == int(r.heatwave))]
        if iy.size and ctrl.size:
            rows.append({"sid": r.sid, "heatwave": r.heatwave, "storm": float(np.nanmean(change[iy, ix])),
                         "control": float(np.nanmean(ctrl))})
    d = pd.DataFrame(rows).dropna()
    d["excess"] = d["storm"] - d["control"]
    y, h = d["excess"].to_numpy(), d["heatwave"].to_numpy()
    est = y[h == 1].mean() - y[h == 0].mean()
    draws = [y[i][h[i] == 1].mean() - y[i][h[i] == 0].mean() for i in storm_draws(d["sid"].to_numpy(), n_boot)]
    m = d.groupby("heatwave")[["storm", "control", "excess"]].mean()
    return {"n_points": int(len(d)), "diff": float(est), "lo": float(np.percentile(draws, 2.5)),
            "hi": float(np.percentile(draws, 97.5)),
            "storm_heatwave": float(m.loc[1.0, "storm"]), "control_heatwave": float(m.loc[1.0, "control"]),
            "excess_heatwave": float(m.loc[1.0, "excess"]), "storm_normal": float(m.loc[0.0, "storm"]),
            "control_normal": float(m.loc[0.0, "control"]), "excess_normal": float(m.loc[0.0, "excess"])}


def report(cfg: Config, refresh: bool = False) -> Path:
    from picopast import storms

    table = storms.build(cfg, refresh=refresh)
    res = run(cfg, table)
    work = cfg.path("interim") / "storms"
    (work / "results.json").write_text(json.dumps(res, indent=1, default=str))
    figs = cfg.path("reports") / "figures"
    plot_points(cfg, table, figs / "cyclone_heatwave_points.png")
    plot_responses(table, figs / "cyclone_heatwave_responses.png")
    n = res["counts"]
    L = ["# When cyclones hit heatwaves: results", "",
         f"Generated {date.today().isoformat()} by `picopast cyclone-heatwave`. Rules: "
         "`docs/cyclone_heatwave_preregistration.md`, committed before any response was computed. Intervals resample "
         "whole storms (2,000 draws).", "",
         f"**Data:** {n['points']:,} track points from {n['storms']} storms (JTWC ≥ 34 kt, ≥ 50 km from land, "
         f"2012–2026). Heatwave points (E ≥ {HOT}): {n['points_heatwave']} from {n['storms_heatwave']} storms; "
         f"normal points (E < {NORMAL}): {n['points_normal']} from {n['storms_normal']} storms; "
         f"{n['points_between']} in between. Bloom measurable at {n['bloom_valid']} points (clouds remove the rest); "
         f"cooling at {n['cooling_valid']}.", "",
         "![points](figures/cyclone_heatwave_points.png)", "",
         "## RQ2: does a heatwave change the bloom? (primary)", ""]
    e = res["rq2"]
    raw = res["rq2_raw"]
    L += [f"**Verdict: {e['verdict']}.** Heatwave effect on the bloom, adjusted for wind, speed, latitude, basin, "
          f"season, pre-storm chlorophyll, sea level and distance to land: **{_fmt(e)}** log10 units "
          f"({e['n_points']} points, {e['storms_heatwave']} storms with heatwave points, "
          f"{e['storms_normal']} with normal points).", "",
          f"Unadjusted: mean bloom {raw['mean_heatwave']:+.3f} over heatwave water and {raw['mean_normal']:+.3f} over "
          f"normal water; difference {_fmt(raw, 'diff')}. (+0.1 log10 = 26% more chlorophyll.)", "",
          "![responses](figures/cyclone_heatwave_responses.png)", "",
          "| sensitivity check | heatwave effect [95%] | points | storms |", "| --- | --- | --- | --- |"]
    for label, s in res["rq2_sensitivity"].items():
        val = _fmt(s) if np.isfinite(s.get("effect", np.nan)) else "n/a"
        L.append(f"| {label} | {val} | {s['n_points']} | {s['n_storms']} |")
    L += ["", "By basin and season (unadjusted means):", "",
          "| basin | season | heatwave points (storms) | normal points (storms) | bloom, heatwave | bloom, normal |",
          "| --- | --- | --- | --- | --- | --- |"]
    for r in res["rq2_strata"]:
        L.append(f"| {r['basin']} | {r['season']} | {r['points_heatwave']} ({r['storms_heatwave']}) | "
                 f"{r['points_normal']} ({r['storms_normal']}) | {r['mean_heatwave']:+.3f} | {r['mean_normal']:+.3f} |")
    c, cr, sv, km = res["rq1_cooling"], res["rq1_cooling_raw"], res["rq1_survival"], res["rq1_rewarm"]
    L += ["", "## RQ1: cooling, heatwave ending and re-warming", "",
          f"- **Cooling** (°C, minimum of D0 to D0+5 minus the D0−7 to D0−3 mean): heatwave effect "
          f"{_fmt(c)} °C, *{c['verdict']}*; unadjusted means {cr['mean_heatwave']:+.2f} (heatwave) and "
          f"{cr['mean_normal']:+.2f} (normal).",
          f"- **Heatwave survival** to D0+5: {sv['storm_share']:.0%} of the heatwave cells under storms stayed in a "
          f"heatwave, against {sv['control_share']:.0%} of far-away heatwave cells on the same days; difference "
          f"{_fmt(sv)} ({sv['n_points']} points, {sv['n_storms']} storms).",
          f"- **Re-warming** (descriptive): back to the pre-storm anomaly within 60 days at "
          f"{km['heatwave']['share_rewarmed']:.0%} of heatwave points (median {_days(km['heatwave']['median_days'])}) "
          f"and {km['normal']['share_rewarmed']:.0%} of normal points (median {_days(km['normal']['median_days'])}). "
          "Over heatwave water the pre-storm anomaly is the heatwave itself, so most of those heatwaves did not return.",
          ""]
    x = res["rq1_excess_cooling_posthoc"]
    L += ["**Post hoc check, not pre-registered** (added after the first results; it does not change the verdict "
          "above). Heatwave water cools on its own as heatwaves decay, storm or no storm. The same cooling measure, "
          "per cell, for far-away water with the same heatwave status on the same days (same basin, more than 750 km "
          "from any storm):", "",
          "| | under the storm | far-away water, same days | storm's own extra cooling |", "| --- | --- | --- | --- |",
          f"| heatwave points | {x['storm_heatwave']:+.2f} °C | {x['control_heatwave']:+.2f} °C | "
          f"{x['excess_heatwave']:+.2f} °C |",
          f"| normal points | {x['storm_normal']:+.2f} °C | {x['control_normal']:+.2f} °C | {x['excess_normal']:+.2f} °C |",
          "",
          f"The storm's own extra cooling is still stronger over heatwave water by **{_fmt(x, 'diff', 2)} °C** "
          f"({x['n_points']} points): about half the unadjusted gap is the heatwave decaying anyway.", "",
          ]
    L += rq3_section(res["rq3_extended"], cfg)
    L += ["### Cross-check: PACE MOANA, storms since 2024 (descriptive)", "",
          "Median change in log10 cells over each storm's points (post minus pre, 150 km). *Synechococcus* and "
          "*Prochlorococcus* from MOANA failed check C in this box: exploratory.", "",
          "| storm | first day | points | picoeukaryotes | Synechococcus | Prochlorococcus | bloom (chl) | heatwave share |",
          "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in res["rq3"]:
        L.append(f"| {r['storm']} | {r['date']} | {r['points']} | {r['d_peuk']:+.2f} | {r['d_syn']:+.2f} | "
                 f"{r['d_pro']:+.2f} | {r['bloom']:+.2f} | {r['heatwave_share']:.0%} |")
    L += ["", "## RQ4: can the response be predicted? (one storm left out at a time)", "",
          "| target | model | skill vs basin × season mean | Spearman ρ |", "| --- | --- | --- | --- |"]
    for target, r in res["rq4"].items():
        for k in ("ridge", "lgbm", "lgbm_no_heatwave"):
            L.append(f"| {target} | {k} | {r[k]['skill']:+.2f} | {r[k]['rho']:+.2f} |")
    for target, r in res["rq4"].items():
        g = r["heatwave_gain"]
        top = ", ".join(list(r["shap"])[:4])
        L.append(f"\n- **{target}:** skill gained from the heatwave inputs {g['gain']:+.3f} "
                 f"[{g['lo']:+.3f}, {g['hi']:+.3f}] ({r['n_points']} points, {r['n_storms']} storms); "
                 f"most influential inputs (SHAP): {top}.")
    L += ["", "## Limits", "",
          "- OISST has no microwave input and under-reads cold wakes beneath cloud (same for both classes).",
          "- 8-day chlorophyll composites smear the bloom's timing; cloud removes points (counts above).",
          "- No mixed-layer depth after 2017; sea level and salinity stand in.",
          "- Heatwave and normal points differ in season and basin; the model and the strata address this, not fully.",
          ""]
    out = cfg.path("reports") / "cyclone_heatwave.md"
    out.write_text("\n".join(L))
    return out


# ---------------------------------------------------------------- figures

HEAT, NORM, MID = "#eb6834", "#2a78d6", "#b9b8b0"
INK, MUTED, SURFACE = "#2b2b29", "#6b6a63", "#fcfcfb"


def _axes_style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def plot_points(cfg: Config, table: pd.DataFrame, out: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from picopast.route_d import box_fields

    bf = box_fields(cfg)
    lat, lon = bf["ds"]["lat"].values, bf["ds"]["lon"].values
    land = ~bf["offshore"]
    p = table[table["primary"]]
    fig, ax = plt.subplots(figsize=(9, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.pcolormesh(lon, lat, np.where(land, 1.0, np.nan), cmap="Greys", vmin=0, vmax=3, shading="nearest")
    for label, cls, col, z in (("between", None, MID, 1), ("normal water (E < 0.1)", 0, NORM, 2),
                               ("heatwave water (E ≥ 0.5)", 1, HEAT, 3)):
        s = p[p["heatwave"].isna()] if cls is None else p[p["heatwave"] == cls]
        ax.scatter(s["lon"], s["lat"], s=10, c=col, lw=0, alpha=0.85, zorder=z, label=f"{label}: {len(s)} points")
    ax.set_xlim(50, 100)
    ax.set_ylim(0, 25)
    ax.set_aspect("equal")
    _axes_style(ax)
    ax.legend(loc="lower left", fontsize=8, frameon=False, labelcolor=INK)
    ax.set_title(f"Cyclone track points, 2012–2026 ({p['sid'].nunique()} storms): heatwave exposure 3 days before",
                 loc="left", fontsize=10, color=INK)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out


def rq3_section(x: dict, cfg: Config) -> list[str]:
    """Report lines for RQ3 extended, plus its figure."""
    L = ["## RQ3: which phytoplankton groups respond? (all storms, PicoPast's records)", "",
         "Change in log10 cells/mL after the storm (post minus pre, anomalies from each record's own 8-day "
         "climatology, 150 km). Primary: the versions built without sea temperature, so storm cooling cannot move "
         "them through the model's temperature input. These are model-implied group responses, not counts, and not "
         "independent of the chlorophyll response; *Prochlorococcus* is exploratory (failed check E).", "",
         "![groups](figures/cyclone_heatwave_groups.png)", "",
         "| group | change after storms [95%] | over normal water | over heatwave water | heatwave effect, adjusted [95%] |",
         "| --- | --- | --- | --- | --- |"]
    for g, r in x["primary"]["groups"].items():
        a, n, h, e = r["all"], r["normal"], r["heatwave"], r["effect"]
        L.append(f"| {GROUPS3[g]} | {a['mean']:+.3f} [{a['lo']:+.3f}, {a['hi']:+.3f}], *{a['verdict']}* "
                 f"({a['n']} points, {a['storms']} storms) | {n['mean']:+.3f} [{n['lo']:+.3f}, {n['hi']:+.3f}] | "
                 f"{h['mean']:+.3f} [{h['lo']:+.3f}, {h['hi']:+.3f}] | {_fmt(e)}, *{e['verdict']}* |")
    sh = x["primary"]["shifts"]
    L += ["", "**Community shift** (per point, then averaged):"]
    for k, label in (("peuk-pro", "picoeukaryotes minus *Prochlorococcus*"),
                     ("syn-pro", "*Synechococcus* minus *Prochlorococcus*")):
        if k in sh:
            s = sh[k]
            L.append(f"- {label}: {s['mean']:+.3f} [{s['lo']:+.3f}, {s['hi']:+.3f}] ({s['n']} points, {s['storms']} storms)")
    L += ["", "Secondary, the main records (they include sea temperature):", "",
          "| group | change after storms [95%] | heatwave effect, adjusted [95%] |", "| --- | --- | --- |"]
    for g, r in x["secondary"]["groups"].items():
        a, e = r["all"], r["effect"]
        L.append(f"| {GROUPS3[g]} | {a['mean']:+.3f} [{a['lo']:+.3f}, {a['hi']:+.3f}], *{a['verdict']}* | "
                 f"{_fmt(e)}, *{e['verdict']}* |")
    if "cross_check" in x:
        c = x["cross_check"]
        L += ["", f"Against PACE MOANA, the PicoPast picoeukaryote change has the same sign for {c['same_sign']} of "
              f"{c['storms']} storms.", ""]
    plot_groups(x["primary"]["groups"], cfg.path("reports") / "figures" / "cyclone_heatwave_groups.png")
    return L


def plot_groups(groups: dict, out: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 3.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for i, (g, r) in enumerate(groups.items()):
        for dx, key, col in ((-0.12, "normal", NORM), (0.12, "heatwave", HEAT)):
            m = r[key]
            ax.errorbar(i + dx, m["mean"], yerr=[[m["mean"] - m["lo"]], [m["hi"] - m["mean"]]], fmt="o", ms=6,
                        color=col, ecolor=col, elinewidth=1.5, capsize=0)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xticks(range(len(groups)), [GROUPS3[g] for g in groups])
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_ylabel("change after the storm\n(log10 cells/mL)", fontsize=8, color=MUTED)
    ax.plot([], [], "o", color=NORM, label="normal water")
    ax.plot([], [], "o", color=HEAT, label="heatwave water")
    ax.legend(loc="upper right", fontsize=8, frameon=False, labelcolor=INK)
    ax.set_title("Phytoplankton groups after cyclones, 2012–2026 (PicoPast records without SST; 95% intervals)",
                 loc="left", fontsize=10, color=INK)
    _axes_style(ax)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out


def plot_responses(table: pd.DataFrame, out: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = table[table["primary"]]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), facecolor=SURFACE)
    for ax, col, title, unit in ((axes[0], "cooling", "Cooling under the storm", "°C"),
                                 (axes[1], "bloom", "Chlorophyll change after the storm", "log10 units")):
        ax.set_facecolor(SURFACE)
        data = [p.loc[(p["heatwave"] == 0) & p[col].notna(), col], p.loc[(p["heatwave"] == 1) & p[col].notna(), col]]
        parts = ax.violinplot(data, positions=[0, 1], showmedians=True, widths=0.8)
        for body, c in zip(parts["bodies"], (NORM, HEAT)):
            body.set_facecolor(c)
            body.set_alpha(0.35)
            body.set_edgecolor(c)
        for k in ("cbars", "cmins", "cmaxes", "cmedians"):
            parts[k].set_color(MUTED)
            parts[k].set_linewidth(1)
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_xticks([0, 1], [f"normal water\n(n {len(data[0])})", f"heatwave water\n(n {len(data[1])})"])
        ax.set_ylabel(unit, fontsize=8, color=MUTED)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        _axes_style(ax)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out
