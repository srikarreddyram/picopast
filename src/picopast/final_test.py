"""Final confirmatory test over the whole Indian Ocean.

Rules: ``docs/cyclone_heatwave_preregistration.md``, section "Final confirmatory test" (binding from
2026-11-27, committed before any South Indian Ocean value or any test here was computed).

T1 (decisive): the pooled heatwave effect on the bloom, 99% storm-resampled interval, smallest meaningful
effect +-0.05 log10. T2 storm mixing power, T3 timing, T4 storm-minus-no-storm control (95%, Holm).
T5 power simulation, descriptive, computed first.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from picopast.config import Config, load_config
from picopast.storms_analysis import storm_draws

N_BOOT = 5000
SEED = 0
SESOI = 0.05
LEVEL_T1 = 0.99
LEVEL_SECONDARY = 0.95
MIN_STORMS = 15
COVS = ["wind", "speed_kmh", "abs_lat", "chl_pre", "ssha_cm", "dist2land_km", "era_modis"]
SOUTH_CONFIG = "config/picopast_sio.yaml"


# ---------------------------------------------------------------- data

def pooled_table(cfg_n: Config, cfg_s: Config | None) -> pd.DataFrame:
    """Primary points of the north (main 2012-, replication 2002-2011) and, if built, the south (both eras)."""
    parts = []
    main = pd.read_parquet(cfg_n.path("interim") / "storms" / "storm_table.parquet").assign(era="viirs")
    rep = pd.read_parquet(cfg_n.path("interim") / "storms" / "storm_table_2002.parquet").assign(era="modis")
    parts += [main, rep]
    if cfg_s is not None:
        path = cfg_s.path("interim") / "storms" / "storm_table_all.parquet"
        if path.exists():
            parts.append(pd.read_parquet(path))
    cols = set.intersection(*(set(p.columns) for p in parts))
    d = pd.concat([p[sorted(cols)] for p in parts], ignore_index=True)
    d = d[d["primary"]].copy()
    d["south"] = (d["basin"] == "South Indian Ocean").astype(float)
    d["era_modis"] = (d["era"] == "modis").astype(float)
    d["abs_lat"] = d["lat"].abs()
    with np.errstate(divide="ignore", invalid="ignore"):
        m = np.log(d["wind"] ** 2 / d["speed_kmh"].where(d["speed_kmh"] > 0))
        m3 = np.log(d["wind"] ** 3 / d["speed_kmh"].where(d["speed_kmh"] > 0))
    d["M"] = (m - m.mean()) / m.std()
    d["M3"] = (m3 - m3.mean()) / m3.std()
    d["hw_x_M"] = d["heatwave"] * d["M"]
    d["hw_x_M3"] = d["heatwave"] * d["M3"]
    d["hw_x_south"] = d["heatwave"] * d["south"]
    d["early"] = d["chl_post1"] - d["chl_pre"]       # T3: first post period minus pre
    d["late"] = d["chl_post2"] - d["chl_post1"]      # T3: second post period minus the first
    if "bloom_control" in d:
        d["excess"] = d["bloom"] - d["bloom_control"]  # T4
    return d.reset_index(drop=True)


# ---------------------------------------------------------------- the model

def design(d: pd.DataFrame, terms: list[str], covs: list[str]) -> tuple[np.ndarray, list[str]]:
    """Intercept, the terms of interest, covariates, a Bay of Bengal indicator and region-season indicators
    (north pre-monsoon is the reference; the south's three seasons span the South Indian Ocean indicator)."""
    cols = {"intercept": np.ones(len(d))}
    for c in terms + covs:
        cols[c] = d[c].to_numpy(float)
    bay = (d["basin"] == "Bay of Bengal").to_numpy(float)
    if 0 < bay.sum() < len(d):
        cols["bay"] = bay
    for s in sorted(set(d["season"]) - {"pre-monsoon"}):
        x = (d["season"] == s).to_numpy(float)
        if 0 < x.sum() < len(d):
            cols[f"season: {s}"] = x
    names = list(cols)
    X = np.column_stack([cols[n] for n in names])
    # drop exact duplicates of earlier columns (e.g. era_modis when only one era is present)
    keep = [j for j in range(X.shape[1]) if j == 0 or np.ptp(X[:, j]) > 0]
    return X[:, keep], [names[j] for j in keep]


def fit(d: pd.DataFrame, response: str, terms: list[str], covs=COVS, n_boot: int = N_BOOT, seed: int = SEED,
        weights: np.ndarray | None = None) -> dict:
    """OLS estimates of `terms` with storm-resampled draws."""
    need = [response, "heatwave"] + terms + list(covs)
    d = d.dropna(subset=list(dict.fromkeys(need))).reset_index(drop=True)
    out = {"response": response, "n_points": int(len(d)), "n_storms": int(d["sid"].nunique()),
           "storms_heatwave": int(d.loc[d["heatwave"] == 1, "sid"].nunique()),
           "storms_normal": int(d.loc[d["heatwave"] == 0, "sid"].nunique())}
    if len(d) < 30 or d["heatwave"].nunique() < 2:
        return out | {"est": {}, "draws": {}}
    X, names = design(d, terms, list(covs))
    y = d[response].to_numpy(float)
    pos = [names.index(t) for t in terms]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    draws = []
    for idx in storm_draws(d["sid"].to_numpy(), n_boot, seed):
        b = np.linalg.lstsq(X[idx], y[idx], rcond=None)[0]
        draws.append(b[pos])
    draws = np.asarray(draws)
    return out | {"est": {t: float(beta[j]) for t, j in zip(terms, pos)},
                  "draws": {t: draws[:, i] for i, t in enumerate(terms)}}


def interval(draws: np.ndarray, level: float) -> tuple[float, float]:
    a = (1 - level) / 2 * 100
    lo, hi = np.percentile(draws, [a, 100 - a])
    return float(lo), float(hi)


def pvalue(draws: np.ndarray) -> float:
    """Two-sided bootstrap p-value for 'the coefficient is zero' (never below 1 / number of draws)."""
    d = np.asarray(draws)
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return float(min(1.0, max(p, 1.0 / d.size)))


def holm(pvals: dict[str, float]) -> dict[str, float]:
    """Holm-adjusted p-values."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, running = len(items), {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        out[k] = running
    return out


def t1_verdict(lo: float, hi: float, sesoi: float = SESOI) -> str:
    inside = -sesoi < lo and hi < sesoi
    if hi < 0:
        return "a small weakening, less than a third of the typical bloom" if inside else "weaker bloom over heatwaves"
    if lo > 0:
        return "a small strengthening, less than a third of the typical bloom" if inside else "stronger bloom over heatwaves"
    if inside:
        return "no meaningful effect"
    return "inconclusive"


def summary(f: dict, term: str, level: float) -> dict:
    if term not in f.get("est", {}):
        return {k: v for k, v in f.items() if k not in ("est", "draws")} | {"effect": np.nan}
    lo, hi = interval(f["draws"][term], level)
    return ({k: v for k, v in f.items() if k not in ("est", "draws")}
            | {"effect": f["est"][term], "lo": lo, "hi": hi, "level": level, "p": pvalue(f["draws"][term])})


# ---------------------------------------------------------------- T5: power simulation

def power_simulation(d_north: pd.DataFrame, n_hw_storms: int, n_norm_storms: int, effects=(0.0, -SESOI),
                     reps: int = 500, n_boot: int = 300, seed: int = SEED) -> dict:
    """Chance of each T1 outcome if the full sample behaved like the north, at planted effects.

    North storms with any heatwave point and those with only normal points are resampled (whole storms) to the
    full sample's counts; the observed heatwave effect is removed and the planted one added; T1 is refitted.
    """
    d = d_north.dropna(subset=["bloom", "heatwave"] + COVS).copy()
    f0 = fit(d, "bloom", ["heatwave"], n_boot=1)
    d["bloom_null"] = d["bloom"] - f0["est"]["heatwave"] * d["heatwave"]
    hw_sids = d.loc[d["heatwave"] == 1, "sid"].unique()
    norm_sids = np.setdiff1d(d["sid"].unique(), hw_sids)
    groups = {s: g for s, g in d.groupby("sid")}
    rng = np.random.default_rng(seed)
    out = {}
    for eff in effects:
        tally: dict[str, int] = {}
        for _ in range(reps):
            pick = list(rng.choice(hw_sids, n_hw_storms)) + list(rng.choice(norm_sids, n_norm_storms))
            sim = pd.concat([groups[s].assign(sid=f"{s}#{i}") for i, s in enumerate(pick)], ignore_index=True)
            sim["bloom"] = sim["bloom_null"] + eff * sim["heatwave"]
            f = fit(sim, "bloom", ["heatwave"], n_boot=n_boot, seed=int(rng.integers(1 << 31)))
            if "heatwave" not in f["est"]:
                continue
            v = t1_verdict(*interval(f["draws"]["heatwave"], LEVEL_T1))
            tally[v] = tally.get(v, 0) + 1
        n = sum(tally.values())
        out[f"{eff:+.2f}"] = {k: v / n for k, v in sorted(tally.items(), key=lambda kv: -kv[1])}
    return {"heatwave_storms": n_hw_storms, "normal_storms": n_norm_storms, "reps": reps, "outcomes": out}


# ---------------------------------------------------------------- everything

def run(cfg_n: Config, cfg_s: Config | None, n_boot: int = N_BOOT, power_reps: int = 500) -> dict:
    d = pooled_table(cfg_n, cfg_s)
    north = d[d["south"] == 0]
    res: dict = {"counts": {
        "points": int(len(d)), "storms": int(d["sid"].nunique()),
        "storms_north": int(north["sid"].nunique()), "storms_south": int(d.loc[d["south"] == 1, "sid"].nunique()),
        "heatwave_storms": int(d.loc[d["heatwave"] == 1, "sid"].nunique()),
        "normal_storms": int(d.loc[d["heatwave"] == 0, "sid"].nunique()),
        "heatwave_points": int((d["heatwave"] == 1).sum()), "normal_points": int((d["heatwave"] == 0).sum()),
        "heatwave_share_by_region_era": d.groupby(["basin", "era"])["heatwave"].apply(
            lambda s: float((s == 1).mean())).round(3).to_dict()}}
    res["counts"]["heatwave_share_by_region_era"] = {f"{k[0]} / {k[1]}": v for k, v in
                                                    res["counts"]["heatwave_share_by_region_era"].items()}
    # T5 first (exposure counts only)
    res["t5_power"] = power_simulation(north, res["counts"]["heatwave_storms"], res["counts"]["normal_storms"],
                                       reps=power_reps)
    # T1
    f1 = fit(d, "bloom", ["heatwave"], n_boot=n_boot)
    t1 = summary(f1, "heatwave", LEVEL_T1)
    t1["verdict"] = t1_verdict(t1["lo"], t1["hi"]) if np.isfinite(t1["effect"]) else "not enough data"
    if t1.get("storms_heatwave", 0) < MIN_STORMS or t1.get("storms_normal", 0) < MIN_STORMS:
        t1["verdict"] = "not enough data"
    t1["pct"] = {k: float(100 * (10 ** t1[k] - 1)) for k in ("effect", "lo", "hi")}
    res["t1"] = t1
    if d["south"].nunique() > 1:
        fh = fit(d, "bloom", ["heatwave", "hw_x_south"], n_boot=n_boot)
        het = summary(fh, "hw_x_south", LEVEL_SECONDARY)
        res["t1_heterogeneity"] = het
        res["t1_regions"] = {}
        for name, sub in (("North Indian Ocean", d[d["south"] == 0]), ("South Indian Ocean", d[d["south"] == 1])):
            s = summary(fit(sub, "bloom", ["heatwave"], n_boot=n_boot), "heatwave", LEVEL_T1)
            if np.isfinite(s.get("effect", np.nan)):
                s["verdict"] = t1_verdict(s["lo"], s["hi"])
            res["t1_regions"][name] = s
        if het.get("hi", 0) < 0 or het.get("lo", 0) > 0:
            res["t1"]["replaced_by_regions"] = True
    # T2-T4
    f2 = fit(d, "bloom", ["heatwave", "M", "hw_x_M"], n_boot=n_boot)
    t2 = summary(f2, "hw_x_M", LEVEL_SECONDARY)
    if f2.get("est"):
        qs = np.nanpercentile(d["M"], [25, 75])
        t2["effect_at_M"] = {}
        for label, q in zip(("weak or fast storms (M 25th pct)", "strong or slow storms (M 75th pct)"), qs):
            dr = f2["draws"]["heatwave"] + q * f2["draws"]["hw_x_M"]
            lo, hi = interval(dr, LEVEL_SECONDARY)
            t2["effect_at_M"][label] = {"effect": float(f2["est"]["heatwave"] + q * f2["est"]["hw_x_M"]),
                                        "lo": lo, "hi": hi}
    t2_m3 = summary(fit(d, "bloom", ["heatwave", "M3", "hw_x_M3"], n_boot=n_boot), "hw_x_M3", LEVEL_SECONDARY)
    t3a = summary(fit(d, "early", ["heatwave"], n_boot=n_boot), "heatwave", LEVEL_SECONDARY)
    t3b = summary(fit(d, "late", ["heatwave"], n_boot=n_boot), "heatwave", LEVEL_SECONDARY)
    t4 = (summary(fit(d, "excess", ["heatwave"], n_boot=n_boot), "heatwave", LEVEL_SECONDARY)
          if "excess" in d else {"effect": np.nan})
    tests = {"T2 heatwave x mixing power": t2, "T3 first period (mostly mixed-up chlorophyll)": t3a,
             "T3 second period (later growth)": t3b, "T4 storm minus no-storm control": t4}
    adj = holm({k: v["p"] for k, v in tests.items() if np.isfinite(v.get("effect", np.nan))})
    for k, v in tests.items():
        v["p_holm"] = adj.get(k, np.nan)
        v["clear"] = bool(np.isfinite(v["p_holm"]) and v["p_holm"] < 0.05)
    t2["verdict"] = ("confirmed: heatwaves weaken blooms mainly under weak or fast storms"
                     if t2["clear"] and t2["effect"] > 0 else
                     "the opposite: the heatwave effect grows with storm power" if t2["clear"]
                     else "no detectable dependence on storm power")
    for k in ("T3 first period (mostly mixed-up chlorophyll)", "T3 second period (later growth)",
              "T4 storm minus no-storm control"):
        v = tests[k]
        v["verdict"] = ("weaker over heatwaves" if v["clear"] and v["effect"] < 0 else
                        "stronger over heatwaves" if v["clear"] else "no detectable difference")
    res["secondary"] = tests
    res["t2_sensitivity_M3"] = t2_m3
    # RQ1 for the whole sample and the south
    rq1 = {}
    for name, sub in (("pooled", d), ("South Indian Ocean", d[d["south"] == 1])):
        if len(sub):
            rq1[name] = summary(fit(sub, "cooling", ["heatwave"], n_boot=min(n_boot, 2000)), "heatwave", 0.95)
    res["rq1_cooling"] = rq1
    return res


def _f(s: dict, k="effect", digits=3) -> str:
    if not np.isfinite(s.get(k, np.nan)):
        return "n/a"
    return f"{s[k]:+.{digits}f} [{s['lo']:+.{digits}f}, {s['hi']:+.{digits}f}]"


def report(cfg_n: Config, cfg_s: Config | None, n_boot: int = N_BOOT, power_reps: int = 500) -> Path:
    res = run(cfg_n, cfg_s, n_boot, power_reps)
    work = cfg_n.path("interim") / "storms"
    (work / "final_results.json").write_text(json.dumps(res, indent=1, default=float))
    n, t1 = res["counts"], res["t1"]
    L = ["# Final confirmatory test: heatwaves, cyclones and the phytoplankton bloom, whole Indian Ocean", "",
         f"Generated {date.today().isoformat()} by `picopast final-test`. Rules: the \"Final confirmatory test\" "
         "section of `docs/cyclone_heatwave_preregistration.md`, committed before any South Indian Ocean value or "
         "any test here was computed. Intervals resample whole storms (5,000 draws).", "",
         f"**Sample:** {n['points']:,} track points from {n['storms']} storms ({n['storms_north']} North, "
         f"{n['storms_south']} South Indian Ocean), 2002–2026. Heatwave points: {n['heatwave_points']} from "
         f"{n['heatwave_storms']} storms; normal points: {n['normal_points']} from {n['normal_storms']} storms.", "",
         "## T1: does a heatwave change the bloom? (decisive)", "",
         f"**Verdict: {t1['verdict']}.** Heatwave effect on the bloom {_f(t1)} log10 (99% interval), i.e. "
         f"{t1['pct']['effect']:+.0f}% [{t1['pct']['lo']:+.0f}%, {t1['pct']['hi']:+.0f}%] chlorophyll; "
         f"{t1['n_points']:,} points, {t1['storms_heatwave']} heatwave storms, {t1['storms_normal']} normal storms. "
         f"The smallest meaningful effect was fixed in advance at ±{SESOI} log10 (a third of the typical bloom).", ""]
    if "t1_heterogeneity" in res:
        h = res["t1_heterogeneity"]
        L += [f"Heterogeneity guard: heatwave × South Indian Ocean interaction {_f(h)} (95%)"
              + (" — **excludes zero, so the regional verdicts below replace the pooled one**." if
                 t1.get("replaced_by_regions") else " — includes zero, so the pooled verdict stands."), "",
              "| region | heatwave effect, 99% | verdict | points | heatwave / normal storms |",
              "| --- | --- | --- | --- | --- |"]
        for name, s in res["t1_regions"].items():
            L.append(f"| {name} | {_f(s)} | {s.get('verdict', 'n/a')} | {s.get('n_points', 0)} | "
                     f"{s.get('storms_heatwave', 0)} / {s.get('storms_normal', 0)} |")
        L.append("")
    p = res["t5_power"]
    L += ["### T5: power simulation (computed before T1, from exposure counts only)", "",
          f"Resampling North Indian Ocean storms to {p['heatwave_storms']} heatwave and {p['normal_storms']} normal "
          f"storms ({p['reps']} runs):", "", "| planted effect | chance of each T1 outcome |", "| --- | --- |"]
    for eff, o in p["outcomes"].items():
        L.append(f"| {eff} log10 | " + "; ".join(f"{k} {v:.0%}" for k, v in o.items()) + " |")
    L += ["", "## T2–T4: mechanism, timing and control (95%, Holm across the four)", "",
          "| test | heatwave effect or interaction [95%] | Holm p | verdict |", "| --- | --- | --- | --- |"]
    for k, v in res["secondary"].items():
        L.append(f"| {k} | {_f(v)} | {v.get('p_holm', np.nan):.3f} | {v.get('verdict', 'n/a')} |")
    t2 = res["secondary"]["T2 heatwave x mixing power"]
    if "effect_at_M" in t2:
        L += ["", "Heatwave effect on the bloom by storm mixing power (T2 model, 95%):", ""]
        for label, e in t2["effect_at_M"].items():
            L.append(f"- {label}: {e['effect']:+.3f} [{e['lo']:+.3f}, {e['hi']:+.3f}]")
    L += ["", f"Sensitivity, M3 = log(wind³ / speed): interaction {_f(res['t2_sensitivity_M3'])}.", "",
          "## RQ1: cooling over heatwave water (95%)", ""]
    for name, s in res["rq1_cooling"].items():
        L.append(f"- {name}: {_f(s, digits=2)} °C ({s.get('n_points', 0)} points)")
    L += ["", "Heatwave share of track points by region and era: "
          + "; ".join(f"{k} {v:.0%}" for k, v in n["heatwave_share_by_region_era"].items()) + ".", ""]
    out = cfg_n.path("reports") / "cyclone_heatwave_final.md"
    out.write_text("\n".join(L))
    return out


def main(build: bool = False, n_boot: int = N_BOOT, power_reps: int = 500) -> Path:
    from picopast import storms

    cfg_n = load_config()
    cfg_s = load_config(SOUTH_CONFIG)
    if build:
        storms.build_region(cfg_s)
        storms.add_controls(cfg_s, cfg_s.path("interim") / "storms" / "storm_table_all.parquet",
                            {"modis": ("modis_aqua_chl", "2001-11-01", "2012-03-31"),
                             "viirs": ("viirs_snpp_chl", None, None)})
        storms.add_controls(cfg_n, cfg_n.path("interim") / "storms" / "storm_table.parquet",
                            {None: ("viirs_snpp_chl", None, None)})
        storms.add_controls(cfg_n, cfg_n.path("interim") / "storms" / "storm_table_2002.parquet",
                            {None: ("modis_aqua_chl", "2001-11-01", "2012-03-31")})
    return report(cfg_n, cfg_s, n_boot, power_reps)
