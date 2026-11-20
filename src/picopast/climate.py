"""RQ5: the back-extended record against monsoon, IOD, ENSO and marine heatwaves (PRD section 6 step 8).

Associations, not causes. Inputs are the regional means of the SNPP
reconstruction on reliable pixels (every model input inside the range seen in
training), as a relative index in SD units. Each group's tier is carried into
every table: only picoeukaryotes passed the trust gate.

- Anomaly: 8-day index minus its 2012-2026 mean for the same 8-day slot of the
  year; monthly anomaly = mean of the month's 8-day anomalies.
- Trend: Mann-Kendall test and Sen's slope on annual means of the core-season
  anomalies (SW monsoon left out: cloud), complete years only. Descriptive.
- Lagged correlation: monthly anomaly vs DMI and ONI, index leading by 0-6
  months; Pearson r with an effective sample size for autocorrelation
  (Bretherton et al. 1999); p-values FDR-corrected (Benjamini-Hochberg) over
  all lag tests together.
- IOD 2019: September-December 2019 anomaly against the same months of
  IOD-neutral years (|DMI| < 0.4 for the September-November mean).
- Marine heatwaves: Hobday definition per OISST pixel (1991-2020 baseline,
  picopast.mhw). An 8-day period is a heatwave period for a region when at
  least half of its ocean pixel-days are in a heatwave; it is compared with
  periods of the same calendar months with fewer than 10%. Interval from
  resampling whole years.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

from picopast.config import Config

GROUPS = ("peuk", "syn", "pro")
NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes"}
REGIONS = {"study_box": "Study box", "bay_of_bengal": "Bay of Bengal",
           "arabian_sea_central_east": "Arabian Sea C/E"}


# ---------------------------------------------------------------- statistics


def mann_kendall(x) -> dict:
    """Two-sided Mann-Kendall test (tie-corrected variance) and Sen's slope per step."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 4:
        return {"n": n, "S": np.nan, "z": np.nan, "p": np.nan, "sen": np.nan}
    s = sum(np.sign(x[j] - x[i]) for i in range(n - 1) for j in range(i + 1, n))
    _, counts = np.unique(x, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(counts * (counts - 1) * (2 * counts + 5))) / 18
    z = (s - np.sign(s)) / np.sqrt(var) if s != 0 else 0.0
    slopes = [(x[j] - x[i]) / (j - i) for i in range(n - 1) for j in range(i + 1, n)]
    return {"n": n, "S": float(s), "z": float(z), "p": float(2 * stats.norm.sf(abs(z))),
            "sen": float(np.median(slopes))}


def fdr_bh(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN kept)."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order] * q.size / np.arange(1, q.size + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty_like(q)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def _lag1(x) -> float:
    x = np.asarray(x, float)
    x = x - x.mean()
    return float(np.sum(x[1:] * x[:-1]) / np.sum(x * x))


def corr_effective(x, y) -> dict:
    """Pearson r with p from an effective sample size n (1 - r1 r2) / (1 + r1 r2) (Bretherton et al. 1999)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 8:
        return {"r": np.nan, "n": n, "n_eff": np.nan, "p": np.nan}
    r = float(np.corrcoef(x, y)[0, 1])
    a, b = _lag1(x), _lag1(y)
    n_eff = float(np.clip(n * (1 - a * b) / (1 + a * b), 3, n))
    t = r * np.sqrt((n_eff - 2) / max(1 - r * r, 1e-12))
    return {"r": r, "n": n, "n_eff": n_eff, "p": float(2 * stats.t.sf(abs(t), n_eff - 2))}


# ---------------------------------------------------------------- inputs


def anomalies(cfg: Config, sensor: str = "snpp") -> pd.DataFrame:
    """8-day anomalies of each group's reliable-pixel index, per region (columns: time, region, groups)."""
    meta = json.loads((cfg.path("interim") / "reconstruct" / f"final_{sensor}_operational.json").read_text())
    s = pd.read_csv(cfg.path("products") / "reconstruction" / f"{sensor}_8D_regional_means.csv", parse_dates=["time"])
    s = s[s["region"].isin(REGIONS)].copy()
    s["slot"] = (s["time"].dt.day_of_year - 1) // 8
    out = s[["time", "region", "slot"]].copy()
    for g in GROUPS:
        z = (s[f"{g}_reliable"] - meta[g]["index_mean"]) / meta[g]["index_sd"]
        out[g] = z - z.groupby([s["region"], s["slot"]]).transform("mean")
    return out.drop(columns="slot")


def monthly(anom: pd.DataFrame) -> pd.DataFrame:
    """Monthly mean anomaly per region (index: month start)."""
    return (anom.set_index("time").groupby("region")[list(GROUPS)].resample("MS").mean().reset_index())


def load_indices(cfg: Config) -> pd.DataFrame:
    raw = cfg.path("raw") / "indices"
    dmi = pd.read_csv(raw / "dmi.csv", parse_dates=["date"]).set_index("date")["dmi"]
    oni = pd.read_csv(raw / "oni.csv", parse_dates=["date"]).set_index("date")["oni"]
    return pd.concat({"DMI": dmi, "ONI": oni}, axis=1)


def heatwave_share(cfg: Config, periods: pd.DatetimeIndex) -> pd.DataFrame:
    """Share of each region's ocean pixel-days in a marine heatwave, per 8-day period (OISST, Hobday)."""
    from picopast import mhw

    paths = sorted((cfg.path("raw") / "oisst").glob("oisst_v2.1_box_*.nc"))
    sst = xr.open_mfdataset(paths, combine="by_coords")["sst"]
    rows = {}
    for region in REGIONS:
        if region == "study_box":
            sub = sst
        else:
            r = cfg["regions"][region]
            sub = sst.sel(lat=slice(r["lat_min"], r["lat_max"]), lon=slice(r["lon_min"], r["lon_max"]))
        time = pd.DatetimeIndex(sub["time"].values)
        n_hot = np.zeros(time.size)
        n_ocean = np.zeros(time.size)
        for i in range(0, sub.sizes["lat"], 20):  # strips of 20 rows keep memory near 1 GB (16 GB machine)
            arr = sub.isel(lat=slice(i, i + 20)).values.astype(np.float32)
            _, thresh = mhw.climatology(arr, time)
            hot = mhw.heatwave_days(arr, time, thresh)
            ocean = np.isfinite(arr)
            n_hot += (hot & ocean).sum(axis=(1, 2))
            n_ocean += ocean.sum(axis=(1, 2))
            del arr, hot, ocean
        daily = pd.Series(n_hot / np.maximum(n_ocean, 1), index=time)
        rows[region] = [daily[p:p + pd.Timedelta(days=7)].mean() for p in periods]
    return pd.DataFrame(rows, index=periods)


def sst_free_maps(cfg: Config, sensor: str = "snpp", group: str = "peuk"):
    """Yield (period start, map in SD units) from a Route A model trained without SST, reliable pixels only.

    SST is a model input, so a heatwave or ENSO association could partly be put
    there by the model. This record sees only what VIIRS sees (bands, ratios, PAR,
    position, season), on the same pixels as the main record.
    """
    from picopast import reconstruct as R

    table = pd.read_parquet(cfg.path("interim") / "reconstruct" / f"training_{sensor}_operational.parquet")
    feats = [f for f in R.FEATURES["A"] if f != "sst"]
    d = table[np.isfinite(table[group])]
    model = R._lgbm().fit(d[feats], d[group].values)
    ranges = {f: (np.nanpercentile(d[f], 0.5), np.nanpercentile(d[f], 99.5))
              for f in feats if f not in ("lat", "lon", "doy_sin", "doy_cos")}
    mean, sd = float(d[group].mean()), float(d[group].std())
    lat, lon, offshore = R._grid_and_offshore(cfg)
    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    oisst = R.Oisst(cfg, lat, lon)
    for p in R.viirs_periods(cfg, sensor).itertuples():
        v = R.load_viirs(cfg, sensor, p.start)
        ok = R._inputs_valid(v, oisst.mean(p.start, p.end)) & offshore  # same pixels as the main record
        idx = np.flatnonzero(ok.ravel())
        X = R.features({k: a.ravel()[idx] for k, a in v.items()}, np.zeros(idx.size, np.float32),
                       lat2.ravel()[idx], lon2.ravel()[idx], p.start + (p.end - p.start) / 2)
        reliable = np.ones(idx.size, bool)
        for f, (a, b) in ranges.items():
            reliable &= (X[f].values >= a) & (X[f].values <= b)
        z = np.full(ok.size, np.nan, np.float32)
        z[idx[reliable]] = (model.predict(X.loc[reliable, feats]) - mean) / sd
        yield p.start, z.reshape(ok.shape), lat, lon


def sst_free_series(cfg: Config, sensor: str = "snpp", group: str = "peuk") -> pd.DataFrame:
    """Regional means of the SST-free record (see sst_free_maps). Cached."""
    out = cfg.path("interim") / "climate" / f"{group}_no_sst_{sensor}_regional.csv"
    if out.exists():
        return pd.read_csv(out, parse_dates=["time"])
    boxes = {"study_box": None} | {k: cfg["regions"][k] for k in REGIONS if k != "study_box"}
    rows, lat2 = [], None
    for start, z, lat, lon in sst_free_maps(cfg, sensor, group):
        if lat2 is None:
            lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
        for name, r in boxes.items():
            sub = z if r is None else z[(lat2 >= r["lat_min"]) & (lat2 <= r["lat_max"])
                                        & (lon2 >= r["lon_min"]) & (lon2 <= r["lon_max"])]
            rows.append({"time": start, "region": name,
                         group: float(np.nanmean(sub)) if np.isfinite(sub).any() else np.nan})
    df = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    return df


def _anomaly_of(series: pd.DataFrame, group: str) -> pd.DataFrame:
    s = series.copy()
    slot = (s["time"].dt.day_of_year - 1) // 8
    s[group] = s[group] - s[group].groupby([s["region"], slot]).transform("mean")
    return s


# ---------------------------------------------------------------- analyses


def trends(cfg: Config, anom: pd.DataFrame) -> pd.DataFrame:
    sw = cfg["seasons"]["sw_monsoon"]
    core = anom[~anom["time"].dt.month.isin(sw)].assign(year=lambda d: d["time"].dt.year)
    counts = core.groupby(["region", "year"]).size()
    full = counts[counts >= 0.9 * counts.max()].reset_index()[["region", "year"]]
    core = core.merge(full, on=["region", "year"])
    rows = []
    for region, d in core.groupby("region"):
        annual = d.groupby("year")[list(GROUPS)].mean()
        for g in GROUPS:
            mk = mann_kendall(annual[g].values)
            rows.append({"region": region, "group": g, "years": f"{annual.index.min()}–{annual.index.max()}",
                         **mk, "sen_per_decade": 10 * mk["sen"]})
    out = pd.DataFrame(rows)
    out["p_fdr"] = fdr_bh(out["p"])
    return out


def lag_correlations(mon: pd.DataFrame, idx: pd.DataFrame, max_lag: int = 6) -> pd.DataFrame:
    rows = []
    for region, d in mon.groupby("region"):
        d = d.set_index("time")[list(GROUPS)]
        for name in ("DMI", "ONI"):
            for lag in range(max_lag + 1):
                lead = idx[name].shift(lag, freq="MS")  # index value from `lag` months earlier
                both = d.join(lead.rename("index"), how="inner")
                for g in GROUPS:
                    rows.append({"region": region, "index": name, "lag_months": lag, "group": g,
                                 **corr_effective(both["index"], both[g])})
    out = pd.DataFrame(rows)
    out["p_fdr"] = fdr_bh(out["p"])
    return out


def iod_2019(mon: pd.DataFrame, idx: pd.DataFrame) -> pd.DataFrame:
    son = idx["DMI"][idx.index.month.isin([9, 10, 11])].groupby(lambda t: t.year).mean()
    neutral = [y for y, v in son.items() if abs(v) < 0.4 and y != 2019]
    rows = []
    for region, d in mon.groupby("region"):
        d = d.set_index("time")[list(GROUPS)]
        sond = d[d.index.month.isin([9, 10, 11, 12])].groupby(lambda t: t.year).mean()
        base = sond.loc[[y for y in neutral if y in sond.index]]
        for g in GROUPS:
            if 2019 not in sond.index or len(base) < 3:
                continue
            v, b = sond.loc[2019, g], base[g]
            rows.append({"region": region, "group": g, "anomaly_2019": v, "neutral_mean": b.mean(),
                         "neutral_sd": b.std(), "z": (v - b.mean()) / b.std(), "n_neutral_years": len(b),
                         "neutral_years": ", ".join(str(y) for y in b.index)})
    return pd.DataFrame(rows)


def heatwave_composite(cfg: Config, anom: pd.DataFrame, share: pd.DataFrame, n_boot: int = 2000,
                       seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for region, d in anom.groupby("region"):
        d = d.set_index("time").join(share[region].rename("hw"))
        hot = d[d["hw"] >= 0.5]
        normal = d[(d["hw"] < 0.1) & d.index.month.isin(hot.index.month.unique())]
        for g in GROUPS:
            h, n = hot[g].dropna(), normal[g].dropna()
            if len(h) < 3 or len(n) < 3:
                rows.append({"region": region, "group": g, "hw_periods": len(h), "normal_periods": len(n)})
                continue
            years = np.unique(np.concatenate([h.index.year, n.index.year]))
            diffs = []
            for _ in range(n_boot):
                pick = rng.choice(years, years.size)
                hh = np.concatenate([h[h.index.year == y].values for y in pick])
                nn = np.concatenate([n[n.index.year == y].values for y in pick])
                if hh.size and nn.size:
                    diffs.append(hh.mean() - nn.mean())
            lo, hi = np.percentile(diffs, [2.5, 97.5])
            rows.append({"region": region, "group": g, "hw_periods": len(h), "hw_years": len(np.unique(h.index.year)),
                         "normal_periods": len(n), "difference": h.mean() - n.mean(), "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- report


def _drift_text(cfg: Config) -> str:
    """SNPP minus NOAA-20 Theil-Sen slope per decade, per group (from the RQ4 comparison)."""
    path = cfg.path("products") / "reconstruction" / "snpp_vs_noaa20_8D.csv"
    if not path.exists():
        return "an unknown amount (run sensor-comparison)"
    from picopast.sensor_comparison import _drift

    pp = pd.read_csv(path, parse_dates=["time"])
    work = cfg.path("interim") / "reconstruct"
    sds = {g: json.loads((work / "final_snpp_operational.json").read_text())[g]["index_sd"] for g in GROUPS}
    return ", ".join(f"{_drift(pp, g, sds[g])['slope']:+.2f} ({NAMES[g]})" for g in GROUPS)


def report(cfg: Config, sensor: str = "snpp") -> Path:
    tiers = cfg["reconstruction"]["tiers"]
    anom = anomalies(cfg, sensor)
    mon = monthly(anom)
    idx = load_indices(cfg)
    periods = pd.DatetimeIndex(sorted(anom["time"].unique()))
    share_path = cfg.path("interim") / "climate" / "heatwave_share_8D.csv"
    if share_path.exists():
        share = pd.read_csv(share_path, index_col=0, parse_dates=True)
    else:
        share_path.parent.mkdir(parents=True, exist_ok=True)
        share = heatwave_share(cfg, periods)
        share.to_csv(share_path)
    tr, lc, iod, hw = trends(cfg, anom), lag_correlations(mon, idx), iod_2019(mon, idx), heatwave_composite(
        cfg, anom, share)
    out_dir = cfg.path("interim") / "climate"
    for name, df in (("trends", tr), ("lag_correlations", lc), ("iod_2019", iod), ("heatwave_composite", hw)):
        df.to_csv(out_dir / f"{name}_{sensor}.csv", index=False)

    tier = {g: "relative index" if tiers[g].startswith("relative") else "EXPLORATORY" for g in GROUPS}
    L = ["# Climate associations (RQ5): first pass", "",
         f"Generated {date.today().isoformat()} by `picopast climate`, from the VIIRS {sensor.upper()} "
         "back-extension, reliable pixels only (every model input inside the range seen in training). Values are "
         "anomalies of the relative index in SD units. Associations, not causes. Only picoeukaryotes passed the "
         "trust gate and agree between the two VIIRS sensors; Synechococcus and Prochlorococcus rows are "
         "exploratory and shown for completeness.", "",
         "## Trends (core seasons, annual means)", "",
         "Mann–Kendall test and Sen's slope on annual means of the anomalies outside the SW monsoon. Descriptive: "
         "a short record and a sensor that may drift. For comparison, SNPP drifts against NOAA-20 by "
         + _drift_text(cfg) + " per decade in the study box (`sensor_comparison.md`): a trend no larger than "
         "that drift is not evidence of change in the ocean.", "",
         "| region | group | tier | years | Sen's slope per decade | p | p (FDR) |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in tr.itertuples():
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {tier[r.group]} | {r.years} | {r.sen_per_decade:+.2f} | "
                 f"{r.p:.3f} | {r.p_fdr:.3f} |")
    L += ["", "## Lagged correlation with DMI (IOD) and ONI (ENSO)", "",
          "Monthly anomalies against the index 0–6 months earlier. Shown: the strongest lag per region, group and "
          "index. p uses an effective sample size for autocorrelation; FDR over all "
          f"{len(lc)} lag tests.", "",
          "| region | group | tier | index | best lag (months) | r | n (effective) | p | p (FDR) |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    best = lc.loc[lc.groupby(["region", "group", "index"])["r"].apply(lambda s: s.abs().idxmax())]
    for r in best.sort_values(["group", "region", "index"]).itertuples():
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {tier[r.group]} | {r.index} | {r.lag_months} | "
                 f"{r.r:+.2f} | {r.n} ({r.n_eff:.0f}) | {r.p:.3f} | {r.p_fdr:.3f} |")
    L += ["", "## The 2019 positive IOD", "",
          "September–December 2019 anomaly against the same months of IOD-neutral years (|DMI| < 0.4, "
          "September–November mean). z = (2019 − neutral mean) / neutral SD.", "",
          "| region | group | tier | 2019 | neutral mean ± SD | z | neutral years |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in iod.itertuples():
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {tier[r.group]} | {r.anomaly_2019:+.2f} | "
                 f"{r.neutral_mean:+.2f} ± {r.neutral_sd:.2f} | {r.z:+.1f} | {r.neutral_years} |")
    L += ["", "## Marine heatwaves", "",
          "Heatwave periods: at least half of the region's ocean pixel-days in a heatwave (Hobday definition, "
          "OISST, 1991–2020 baseline). Normal periods: the same calendar months with under 10%. Difference = "
          "heatwave mean − normal mean; 95% interval from resampling whole years.", "",
          "| region | group | tier | heatwave periods (years) | normal periods | difference | 95% interval |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in hw.itertuples():
        if pd.isna(getattr(r, "difference", np.nan)):
            L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {tier[r.group]} | {r.hw_periods} | "
                     f"{r.normal_periods} | too few | — |")
            continue
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {tier[r.group]} | {r.hw_periods} ({r.hw_years}) | "
                 f"{r.normal_periods} | {r.difference:+.2f} | [{r.lo:+.2f}, {r.hi:+.2f}] |")
    no_sst = cfg.path("interim") / "climate" / f"peuk_no_sst_{sensor}_regional.csv"
    if no_sst.exists():
        a2 = _anomaly_of(pd.read_csv(no_sst, parse_dates=["time"]), "peuk")
        a2 = a2.assign(syn=np.nan, pro=np.nan)
        hw2 = heatwave_composite(cfg, a2, share)
        lc2 = lag_correlations(monthly(a2), idx)
        L += ["", "### Sensitivity: picoeukaryotes from a model without SST", "",
              "SST is a model input, so part of a heatwave or ENSO association could be put there by the model. "
              "This record comes from the same model trained without SST (VIIRS bands, ratios, PAR, position, "
              "season only).", "",
              "| region | heatwave difference [95%] | best ONI lag: r (p) | best DMI lag: r (p) |",
              "| --- | --- | --- | --- |"]
        for region in REGIONS:
            h = hw2[(hw2["region"] == region) & (hw2["group"] == "peuk")].iloc[0]
            cells = []
            for name in ("ONI", "DMI"):
                q = lc2[(lc2["region"] == region) & (lc2["group"] == "peuk") & (lc2["index"] == name)]
                b = q.loc[q["r"].abs().idxmax()]
                cells.append(f"lag {b.lag_months}: {b.r:+.2f} ({b.p:.3f})")
            L.append(f"| {REGIONS[region]} | {h.difference:+.2f} [{h.lo:+.2f}, {h.hi:+.2f}] | {cells[0]} | {cells[1]} |")
    hw_years = share[share >= 0.5].count()
    L += ["", "Heatwave periods per region: " + ", ".join(f"{REGIONS[k]} {v}" for k, v in hw_years.items())
          + f" of {len(share)} 8-day periods.", "",
          "## Not yet done", "",
          "- Mixed-layer depth as a covariate (needs the Copernicus Marine account).",
          "- The same analysis on the NOAA-20 record for 2018–2026 (picoeukaryotes), as a sensor check.",
          "- Seasonal (monsoon-phase) composites beyond the climatology.", "",
          "Pixel-level maps of the ENSO, IOD and heatwave composites, with the SST-free check and a NOAA-20 "
          "comparison: `climate_maps.md`.", ""]
    out = cfg.path("reports") / "climate.md"
    out.write_text("\n".join(L))
    return out
