"""Route D: a globally ship-trained, Indian-Ocean-adjusted model (docs/route_d_preregistration.md).

Training: MAREDAT picophytoplankton counts (Buitenhuis et al. 2012) at <= 10 m,
1980 or later, with every Indian Ocean sample removed. Inputs a satellite era
can supply: SST (OISST, the sample's month), log10 chlorophyll and PAR
(MODIS-Aqua monthly climatologies), surface nitrate (WOA23 monthly
climatology) and season (hemisphere-aligned). No position.

- D1: LightGBM, fixed settings, unweighted.
- D2 (primary): D1 with importance weights that make the training samples
  resemble the box's PACE-era conditions.
- D3 (secondary): D2 plus the Indian Ocean stations of the other cruises,
  leave one cruise out.

Each model builds a same-calendar-month field over the box from PACE-era
inputs (OISST, VIIRS SNPP chlorophyll and PAR, WOA nitrate) and goes through
check C unchanged (fallback_run.check_c_fields).
"""

from __future__ import annotations

import json
import re
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr

from picopast.config import Config

GROUPS = ("peuk", "syn", "pro")
NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes"}
MAREDAT_COLS = {"pro": "PromL", "syn": "SynmL", "peuk": "PEukmL"}
FEATURES = ["sst", "logchl", "par", "no3", "season_sin", "season_cos"]
LGBM = {"n_estimators": 300, "num_leaves": 31, "learning_rate": 0.05, "min_child_samples": 50, "verbose": -1,
        "n_jobs": 6, "random_state": 0}
D3_WEIGHT = 5.0


# ---------------------------------------------------------------- small pieces


def in_indian_ocean(lat, lon) -> np.ndarray:
    """The pre-registered exclusion: longitude 20-120°E and latitude <= 31°N (longitudes in -180..180)."""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    return (lon >= 20) & (lon <= 120) & (lat <= 31)


def basin(lat, lon) -> np.ndarray:
    """Rough basin label for the leave-one-basin-out diagnostic: indian / atlantic / other."""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    med = (lon > -6) & (lon < 37) & (lat > 30) & (lat < 46)
    atlantic = (((lon >= -70) & (lon <= 20)) | ((lon >= -100) & (lon < -70) & (lat >= 9))) & ~med
    return np.where(in_indian_ocean(lat, lon), "indian", np.where(atlantic, "atlantic", "other"))


def season(month, lat) -> tuple[np.ndarray, np.ndarray]:
    """Sine and cosine of the mid-month day of year, shifted half a year south of the equator."""
    doy = (np.asarray(month, float) - 0.5) * 365.25 / 12 + np.where(np.asarray(lat, float) < 0, 365.25 / 2, 0.0)
    ang = 2 * np.pi * doy / 365.25
    return np.sin(ang).astype(np.float32), np.cos(ang).astype(np.float32)


def window_median(arr: np.ndarray, glat: np.ndarray, glon: np.ndarray, plat, plon, half: int = 1) -> np.ndarray:
    """Median of finite values in a (2 half + 1)^2 window around the nearest cell; longitudes wrap.

    arr: (lat, lon) on a regular grid; glon and plon in the same convention.
    """
    plat, plon = np.asarray(plat, float), np.asarray(plon, float)
    iy = np.abs(glat[None, :] - plat[:, None]).argmin(1)
    dlon = glon[1] - glon[0]
    ix = np.round((plon - glon[0]) / dlon).astype(int) % glon.size
    off = np.arange(-half, half + 1)
    yy = np.clip(iy[:, None, None] + off[None, :, None], 0, glat.size - 1)
    xx = (ix[:, None, None] + off[None, None, :]) % glon.size
    vals = arr[yy, xx].reshape(len(plat), -1).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN windows (land) give NaN
        return np.nanmedian(vals, axis=1)


def importance_weights(X_train: np.ndarray, X_target: np.ndarray, clip_pct: float = 99.0) -> np.ndarray:
    """p(target|x) / p(train|x) from a quadratic logistic classifier, clipped and normalised to mean 1."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import PolynomialFeatures, StandardScaler

    X = np.vstack([X_train, X_target])
    y = np.r_[np.zeros(len(X_train)), np.ones(len(X_target))]
    clf = make_pipeline(StandardScaler(), PolynomialFeatures(2), LogisticRegression(max_iter=3000, C=1.0))
    clf.fit(X, y)
    p = np.clip(clf.predict_proba(X_train)[:, 1], 1e-6, 1 - 1e-6)
    w = p / (1 - p) * len(X_train) / len(X_target)
    w = np.minimum(w, np.percentile(w, clip_pct))
    return w / w.mean()


def effective_n(w: np.ndarray) -> float:
    return float(w.sum() ** 2 / (w ** 2).sum())


# ---------------------------------------------------------------- training data


def load_maredat(cfg: Config) -> pd.DataFrame:
    """Surface (<= 10 m) MAREDAT stations, 1980 or later; one row per station (its shallowest sample)."""
    d = pd.read_excel(cfg.path("raw") / "maredat_pico" / "extracted" / "picophyto111130.xlsx")
    d = d.rename(columns={"Lat": "lat", "Long": "lon", "Depth": "depth", "Investigators": "investigator"})
    d = d[(d["depth"] <= 10) & (d["year"] >= 1980) & d["month"].between(1, 12)
          & d["lat"].between(-90, 90) & d["lon"].between(-180, 180)]
    d = d.sort_values("depth").groupby([d["lat"].round(3), d["lon"].round(3), "year", "month", "day"],
                                       as_index=False).first()
    out = d[["lat", "lon", "year", "month", "investigator"]].copy()
    for g, col in MAREDAT_COLS.items():
        out[g] = np.where(d[col] > 0, d[col], np.nan)
    out["basin"] = basin(out["lat"], out["lon"])
    return out.reset_index(drop=True)


class GlobalInputs:
    """Route D inputs at any ocean point and month (training samples and Indian Ocean stations)."""

    def __init__(self, cfg: Config):
        raw = cfg.path("raw")
        with xr.open_dataset(raw / "oisst_v2_1deg" / "sst.mnmean.nc") as ds:
            self.sst = ds["sst"].load()
        with xr.open_dataset(raw / "oisst_v2_1deg" / "lsmask.nc") as m:
            ocean = m["mask"].squeeze().values > 0
        self.sst = self.sst.where(ocean)
        self.sst_time = pd.DatetimeIndex(self.sst["time"].values)
        rx = re.compile(r"AQUA_MODIS\.\d{4}(\d{2})01_\d{8}\.L3m\.MC\.(CHL|PAR)\.")
        self.modis = {}
        for p in sorted((raw / "modis_clim").glob("AQUA_MODIS.*.L3m.MC.*.9km.nc")):
            k = rx.search(p.name)
            if k:
                self.modis[(k.group(2), int(k.group(1)))] = p
        self.woa = {m: raw / "woa23_nitrate" / f"woa23_all_n{m:02d}_01.nc" for m in range(1, 13)}

    def at(self, lat, lon, year, month) -> pd.DataFrame:
        lat, lon = np.asarray(lat, float), np.asarray(lon, float)
        year, month = np.asarray(year, int), np.asarray(month, int)
        n = lat.size
        out = {f: np.full(n, np.nan, np.float32) for f in FEATURES}
        # SST: the sample's own month, 1° OISST (longitudes 0-360), nearest ocean cell within one cell.
        t = pd.to_datetime(pd.DataFrame({"year": year, "month": month, "day": np.ones(n, int)}))
        ti = self.sst_time.get_indexer(t)
        glat, glon = self.sst["lat"].values, self.sst["lon"].values
        arr = self.sst.values
        for k in np.unique(ti[ti >= 0]):
            sel = ti == k
            out["sst"][sel] = window_median(arr[k], glat, glon, lat[sel], lon[sel] % 360, half=1)
        # Before OISST (MAREDAT's year 1980, one contributor, appears to mean "year unknown"): the
        # 1982-2011 mean of the same calendar month (change log, 2026-11-20).
        for m in np.unique(month[ti < 0]):
            sel = (ti < 0) & (month == m)
            base = (self.sst_time.month == m) & (self.sst_time.year >= 1982) & (self.sst_time.year <= 2011)
            clim = np.nanmean(arr[base], axis=0)
            out["sst"][sel] = window_median(clim, glat, glon, lat[sel], lon[sel] % 360, half=1)
        for m in np.unique(month):
            sel = month == m
            for prod, name, log in (("CHL", "logchl", True), ("PAR", "par", False)):
                path = self.modis.get((prod, int(m)))
                if path is None:
                    continue
                with xr.open_dataset(path) as ds:
                    var = "chlor_a" if prod == "CHL" else "par"
                    a = ds[var].values.astype(np.float32)
                    v = window_median(a, ds["lat"].values, ds["lon"].values, lat[sel], lon[sel], half=1)
                out[name][sel] = np.log10(np.where(v > 0, v, np.nan)) if log else v
            with xr.open_dataset(self.woa[int(m)], decode_times=False) as ds:
                a = ds["n_an"].isel(time=0, depth=0).values.astype(np.float32)
                out["no3"][sel] = window_median(a, ds["lat"].values, ds["lon"].values, lat[sel], lon[sel], half=1)
        out["season_sin"], out["season_cos"] = season(month, lat)
        return pd.DataFrame(out)


def training_table(cfg: Config, refresh: bool = False) -> pd.DataFrame:
    out = cfg.path("interim") / "route_d" / "training.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    st = load_maredat(cfg)
    X = GlobalInputs(cfg).at(st["lat"], st["lon"], st["year"], st["month"])
    t = pd.concat([st, X], axis=1)
    out.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(out, index=False)
    return t


# ---------------------------------------------------------------- box fields (PACE era, by calendar month)


def box_fields(cfg: Config) -> dict:
    """{feature: (month, lat, lon) float32} on the analysis grid, PACE-era means; plus lat, lon, offshore."""
    from picopast.fallback_run import viirs_monthly_climatology
    from picopast.reconstruct import _grid_and_offshore

    cache = cfg.path("interim") / "route_d" / "box_fields.nc"
    if cache.exists():
        ds = xr.open_dataset(cache).load()
        return {"ds": ds, "offshore": ds["offshore"].values.astype(bool)}
    lat, lon, offshore = _grid_and_offshore(cfg)
    like = xr.DataArray(np.zeros((lat.size, lon.size)), dims=("lat", "lon"), coords={"lat": lat, "lon": lon})
    chl = viirs_monthly_climatology(cfg, like, offshore, "chl", log=True)
    par = viirs_monthly_climatology(cfg, like, offshore, "par", log=False)
    months = chl["month"].values
    pace = cfg["periods"]["pace"]
    years = range(pd.Timestamp(pace["start"]).year, pd.Timestamp(pace["end"]).year + 1)
    sst = xr.concat([xr.open_dataset(cfg.path("raw") / "oisst" / f"oisst_v2.1_box_{y}.nc")["sst"] for y in years],
                    "time").sel(time=slice(pace["start"], pace["end"]))
    sst_m = sst.groupby("time.month").mean("time").sel(month=months)
    olat, olon = sst_m["lat"].values, sst_m["lon"].values
    iy = np.clip(np.round((lat - olat[0]) / (olat[1] - olat[0])).astype(int), 0, olat.size - 1)
    ix = np.clip(np.round((lon - olon[0]) / (olon[1] - olon[0])).astype(int), 0, olon.size - 1)
    sst_g = sst_m.values[:, iy[:, None], ix[None, :]]
    no3 = []
    for m in months:
        with xr.open_dataset(cfg.path("raw") / "woa23_nitrate" / f"woa23_all_n{int(m):02d}_01.nc",
                             decode_times=False) as ds:
            a = ds["n_an"].isel(time=0, depth=0)
            a = a.sel(lat=slice(lat.min() - 2, lat.max() + 2), lon=slice(lon.min() - 2, lon.max() + 2)).load()
        # Nearest 1° value, then fill coastal gaps from the 3x3 neighbourhood.
        ga, go = a["lat"].values, a["lon"].values
        lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
        v = window_median(a.values, ga, go, lat2.ravel(), lon2.ravel(), half=1).reshape(lat2.shape)
        no3.append(v)
    s_sin, s_cos = season(months, np.full(months.size, 10.0))
    ds = xr.Dataset(
        {"sst": (("month", "lat", "lon"), sst_g.astype(np.float32)),
         "logchl": (("month", "lat", "lon"), chl.values.astype(np.float32)),
         "par": (("month", "lat", "lon"), par.values.astype(np.float32)),
         "no3": (("month", "lat", "lon"), np.stack(no3).astype(np.float32)),
         "season_sin": (("month",), s_sin), "season_cos": (("month",), s_cos),
         "offshore": (("lat", "lon"), offshore.astype(np.uint8))},
        coords={"month": months, "lat": lat, "lon": lon})
    cache.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(cache)
    return {"ds": ds, "offshore": offshore}


def _box_rows(fields: dict, months=None):
    """Feature rows for every offshore pixel-month with all inputs finite: (frame, (month idx, iy, ix))."""
    ds, offshore = fields["ds"], fields["offshore"]
    mi = np.arange(ds.sizes["month"]) if months is None else np.flatnonzero(np.isin(ds["month"].values, months))
    cols = {}
    grid = {f: ds[f].values[mi] for f in ("sst", "logchl", "par", "no3")}
    ok = offshore[None] & np.all([np.isfinite(v) for v in grid.values()], axis=0)
    k, y, x = np.nonzero(ok)
    for f, v in grid.items():
        cols[f] = v[k, y, x]
    cols["season_sin"] = ds["season_sin"].values[mi][k]
    cols["season_cos"] = ds["season_cos"].values[mi][k]
    return pd.DataFrame(cols)[FEATURES], (mi[k], y, x)


def predict_field(model, fields: dict, months=None) -> xr.DataArray:
    """A model's same-calendar-month field over the box (log10 cells/mL)."""
    ds = fields["ds"]
    X, (k, y, x) = _box_rows(fields, months)
    out = np.full((ds.sizes["month"], ds.sizes["lat"], ds.sizes["lon"]), np.nan, np.float32)
    for i in range(0, len(X), 2_000_000):
        out[k[i:i + 2_000_000], y[i:i + 2_000_000], x[i:i + 2_000_000]] = model.predict(X.iloc[i:i + 2_000_000])
    da = xr.DataArray(out, dims=("month", "lat", "lon"),
                      coords={"month": ds["month"].values, "lat": ds["lat"].values, "lon": ds["lon"].values})
    return da if months is None else da.sel(month=[m for m in ds["month"].values if m in set(months)])


# ---------------------------------------------------------------- models


def _fit(X: pd.DataFrame, y: np.ndarray, w: np.ndarray | None = None):
    import lightgbm as lgb

    return lgb.LGBMRegressor(**LGBM).fit(X, y, sample_weight=w)


def indian_stations(cfg: Config, gi: GlobalInputs) -> pd.DataFrame:
    """The check C stations with Route D inputs at their own month (for D3's fine-tuning)."""
    from picopast.fallback_run import historical_surface_stations

    st = historical_surface_stations(cfg)
    X = gi.at(st["lat"], st["lon"], st["time_utc"].dt.year, st["time_utc"].dt.month)
    return pd.concat([st.reset_index(drop=True), X], axis=1)


def evaluate(cfg: Config) -> dict:
    from picopast.fallback import paired_spearman_difference, spearman_bootstrap
    from picopast.fallback_run import check_c, check_c_fields

    work = cfg.path("interim") / "route_d"
    t = training_table(cfg)
    train = t[(t["basin"] != "indian") & t[FEATURES].notna().all(axis=1)].reset_index(drop=True)
    fields = box_fields(cfg)
    Xb, _ = _box_rows(fields)
    rng = np.random.default_rng(0)
    box_sample = Xb.iloc[rng.choice(len(Xb), min(50_000, len(Xb)), replace=False)].to_numpy()
    res = {"n_maredat_surface": int(len(t)), "n_indian_excluded": int((t["basin"] == "indian").sum()),
           "n_train_complete": int(len(train)), "groups": {}}
    moana = check_c(cfg, cfg["reconstruction"]["mapping"])
    chl_field = fields["ds"]["logchl"]
    gi = None
    for g in GROUPS:
        d = train[train[g].notna()].reset_index(drop=True)
        y = np.log10(d[g].to_numpy())
        w = importance_weights(d[FEATURES].to_numpy(), box_sample)
        r = {"n_train": int(len(d)), "ess_weighted": effective_n(w),
             "by_basin": d["basin"].value_counts().to_dict(),
             "weight_by_basin": pd.Series(w).groupby(d["basin"].values).sum().div(w.sum()).round(3).to_dict()}
        # Leave-one-basin-out (diagnostic, D1).
        lobo = {}
        for held in ("atlantic", "other"):
            tr, te = d["basin"] != held, d["basin"] == held
            if tr.sum() > 100 and te.sum() > 100:
                p = _fit(d.loc[tr, FEATURES], y[tr]).predict(d.loc[te, FEATURES])
                lobo[held] = float(spearmanr(p, y[te]).statistic)
        r["lobo"] = lobo
        models = {"D1": _fit(d[FEATURES], y), "D2": _fit(d[FEATURES], y, w)}
        r["importance"] = {k: dict(zip(FEATURES, (m.feature_importances_ / m.feature_importances_.sum()).round(3)
                                       .tolist())) for k, m in models.items()}
        checks = {}
        for name, m in models.items():
            field = predict_field(m, fields)
            field.to_netcdf(work / f"field_{name}_{g}.nc")
            checks[name] = check_c_fields(cfg, {g: field}, fields["offshore"])["groups"][g]
        # D3: D2 plus the other cruises' Indian stations, leave one cruise out (secondary).
        if gi is None:
            gi = GlobalInputs(cfg)
            ind = indian_stations(cfg, gi)
        col = f"{g}_cells_ml"
        di = ind[(ind[col] > 0) & ind[FEATURES].notna().all(axis=1)]
        parts = []
        for cruise in sorted(di["cruise"].unique()):
            other = di[di["cruise"] != cruise]
            X3 = pd.concat([d[FEATURES], other[FEATURES]], ignore_index=True)
            y3 = np.r_[y, np.log10(other[col].to_numpy())]
            w3 = np.r_[w, np.full(len(other), D3_WEIGHT)]
            months = sorted(di.loc[di["cruise"] == cruise, "time_utc"].dt.month.unique())
            f3 = predict_field(_fit(X3, y3, w3), fields, months)
            mm = check_c_fields(cfg, {g: f3}, fields["offshore"])["groups"][g]["matchups"]
            parts.append(mm[mm["cruise"] == cruise])
        m3 = pd.concat(parts, ignore_index=True)
        stats3 = spearman_bootstrap(m3["sat_value"], m3[col], m3["cruise"], n_boot=cfg["trust_gate"]["n_bootstrap"])
        r["D3"] = {"stats": stats3}
        # Paired comparisons on stations matched by both.
        mo = moana["groups"][g]["matchups"][["station_key", "sat_value"]].rename(columns={"sat_value": "moana"})
        chl = check_c_fields(cfg, {g: chl_field}, fields["offshore"])["groups"][g]["matchups"][
            ["station_key", "sat_value"]].rename(columns={"sat_value": "chl"})
        r["moana_stats"] = moana["groups"][g]["stats"]
        r["moana_verdict"] = moana["groups"][g]["verdict"]
        for name, c in checks.items():
            m = c["matchups"].merge(mo, on="station_key").merge(chl, on="station_key")
            vs_moana = paired_spearman_difference(m["sat_value"].values, m["moana"].values, m[col].values,
                                                  m["cruise"].values, n_boot=cfg["trust_gate"]["n_bootstrap"])
            sign = -1.0 if spearman_bootstrap(m["chl"], m[col], m["cruise"], n_boot=10)["rho"] < 0 else 1.0
            vs_chl = paired_spearman_difference(m["sat_value"].values, sign * m["chl"].values, m[col].values,
                                                m["cruise"].values, n_boot=cfg["trust_gate"]["n_bootstrap"])
            r[name] = {"stats": c["stats"], "verdict": c["verdict"], "vs_moana": vs_moana, "vs_chl": vs_chl,
                       "chl_sign": sign}
        res["groups"][g] = r
    clean = json.loads(json.dumps(res, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    (work / "results.json").write_text(json.dumps(clean, indent=1))
    return clean


def fragility(cfg: Config, refresh: bool = False) -> pd.DataFrame:
    """D2 refitted without each heavily weighted contributor, and without the Pacific warm pool (descriptive).

    Added after the pre-registered test, because a handful of samples carry most of D2's weight; it does not
    change the decision.
    """
    from picopast.fallback_run import check_c_fields

    out = cfg.path("interim") / "route_d" / "fragility.csv"
    if out.exists() and not refresh:
        return pd.read_csv(out)
    t = training_table(cfg)
    train = t[(t["basin"] != "indian") & t[FEATURES].notna().all(axis=1)].reset_index(drop=True)
    fields = box_fields(cfg)
    Xb, _ = _box_rows(fields)
    rng = np.random.default_rng(0)
    box_sample = Xb.iloc[rng.choice(len(Xb), min(50_000, len(Xb)), replace=False)].to_numpy()
    drops = {
        "(none)": lambda d: np.zeros(len(d), bool),
        "Blanchot (all entries)": lambda d: d["investigator"].fillna("").str.startswith("Blanchot").to_numpy(),
        "Landry": lambda d: (d["investigator"] == "Landry").to_numpy(),
        "Buck": lambda d: (d["investigator"] == "Buck").to_numpy(),
        "Shimada/Maruyama": lambda d: (d["investigator"] == "Shimada/Maruyama").to_numpy(),
        "Pacific warm pool (east of 120°E or west of 150°W, within 20° of the equator)":
            lambda d: (((d["lon"] > 120) | (d["lon"] < -150)) & (d["lat"].abs() < 20)).to_numpy(),
    }
    rows = []
    for g in GROUPS:
        d0 = train[train[g].notna()].reset_index(drop=True)
        for name, rule in drops.items():
            d = d0[~rule(d0)].reset_index(drop=True)
            w = importance_weights(d[FEATURES].to_numpy(), box_sample)
            m = _fit(d[FEATURES], np.log10(d[g].to_numpy()), w)
            r = check_c_fields(cfg, {g: predict_field(m, fields)}, fields["offshore"])["groups"][g]
            rows.append({"group": g, "left_out": name, "n_train": len(d), "ess": effective_n(w),
                         "rho": r["stats"]["rho"], "lo": r["stats"]["lo"], "hi": r["stats"]["hi"],
                         "verdict": r["verdict"]})
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    return df


# ---------------------------------------------------------------- report


def decision(res: dict) -> dict[str, str]:
    """The pre-registered decision, per group."""
    out = {}
    for g in GROUPS:
        r = res["groups"][g]
        d2 = r["D2"]
        if g == "peuk":
            better = d2["vs_moana"]["diff_lo"] > 0
            out[g] = ("Route D replaces MOANA as the base" if better else
                      "MOANA stays the base (D2 is not clearly better on the paired test)")
        else:
            out[g] = ("reconstruct from D2 as a relative index (Route D)" if d2["verdict"] == "pattern consistent"
                      else "stays exploratory; Route D is a negative result for this group")
    return out


def _ci(s: dict) -> str:
    return f"{s['rho']:.2f} [{s['lo']:.2f}, {s['hi']:.2f}]"


def report(cfg: Config, refresh: bool = False) -> Path:
    path = cfg.path("interim") / "route_d" / "results.json"
    res = evaluate(cfg) if refresh or not path.exists() else json.loads(path.read_text())
    dec = decision(res)
    L = ["# Route D: a globally trained, Indian-Ocean-adjusted model", "",
         f"Generated {date.today().isoformat()} by `picopast route-d`. Rules: `docs/route_d_preregistration.md` "
         "(committed before this test). Trained on MAREDAT ship counts from the rest of the world ocean "
         f"({res['n_train_complete']:,} surface stations with every input; {res['n_indian_excluded']} Indian Ocean "
         "stations excluded), inputs SST, chlorophyll, PAR, surface nitrate and season, no position. Tested "
         "exactly like MOANA in check C: the model's PACE-era same-month field against the 1991–2005 Indian "
         "Ocean stations, Spearman ρ with 95% intervals from resampling cruises.", "",
         "## Check C, side by side", "",
         "| group | MOANA | D1 (global) | D2 (adjusted, primary) | D3 (fine-tuned, leave one cruise out) | "
         "D2 verdict |", "| --- | --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        r = res["groups"][g]
        L.append(f"| {NAMES[g]} | {_ci(r['moana_stats'])} ({r['moana_verdict']}) | {_ci(r['D1']['stats'])} | "
                 f"{_ci(r['D2']['stats'])} | {_ci(r['D3']['stats'])} | {r['D2']['verdict']} |")
    L += ["", "## Paired tests (stations matched by every field)", "",
          "| group | model | stations (cruises) | ρ(model) − ρ(MOANA) | ρ(model) − ρ(chlorophyll) |",
          "| --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        r = res["groups"][g]
        for name in ("D1", "D2"):
            a, b = r[name]["vs_moana"], r[name]["vs_chl"]
            chl = "chlorophyll" if r[name]["chl_sign"] > 0 else "low chlorophyll"
            L.append(f"| {NAMES[g]} | {name} | {a['n']} ({a['n_clusters']}) | {a['diff']:+.2f} "
                     f"[{a['diff_lo']:+.2f}, {a['diff_hi']:+.2f}] | {b['diff']:+.2f} [{b['diff_lo']:+.2f}, "
                     f"{b['diff_hi']:+.2f}] vs {chl} |")
    L += ["", "## Decision (pre-registered rule)", ""]
    L += [f"- **{NAMES[g]}:** {dec[g]}." for g in GROUPS]
    L += ["", "## Diagnostics", "",
          "| group | training stations (by basin) | effective n after weighting | weight share by basin | "
          "leave-one-basin-out ρ | input importance (D2) |", "| --- | --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        r = res["groups"][g]
        imp = ", ".join(f"{k} {v:.2f}" for k, v in sorted(r["importance"]["D2"].items(), key=lambda kv: -kv[1]))
        lobo = ", ".join(f"held-out {k} {v:.2f}" for k, v in r["lobo"].items())
        L.append(f"| {NAMES[g]} | {r['n_train']:,} ({', '.join(f'{k} {v}' for k, v in r['by_basin'].items())}) | "
                 f"{r['ess_weighted']:.0f} | {', '.join(f'{k} {v:.0%}' for k, v in r['weight_by_basin'].items())} | "
                 f"{lobo} | {imp} |")
    frag = fragility(cfg)
    L += ["", "### Does the pass rest on a few samples? (added after the test; descriptive)", "",
          "After weighting, a few dozen samples carry most of D2's weight, mostly from the western equatorial "
          "Pacific warm pool (SST about 29 °C, low chlorophyll and nitrate: the closest analogue to this box). "
          "D2 refitted without each heavily weighted contributor, and without the warm pool altogether:", "",
          "| left out of training | " + " | ".join(NAMES[g] for g in GROUPS) + " |",
          "| --- | " + " | ".join("---" for _ in GROUPS) + " |"]
    for name, part in frag.groupby("left_out", sort=False):
        cells = []
        for g in GROUPS:
            q = part[part["group"] == g].iloc[0]
            mark = "✓" if q["verdict"] == "pattern consistent" else "✗"
            cells.append(f"{q['rho']:.2f} [{q['lo']:.2f}, {q['hi']:.2f}] {mark}")
        L.append(f"| {name} | " + " | ".join(cells) + " |")
    L += ["", "✓ = still *pattern consistent* under the check C rule."]
    L += ["", "Leave-one-basin-out: D1 trained without one basin and scored on it (Spearman ρ on the samples), "
          "a check that conditions, not positions, carry the signal between oceans.", "",
          "## Limits", "",
          "- Training counts are 1987–2011 and chlorophyll and PAR are climatologies: the model learns typical "
          "conditions, not same-day ones. Counts are 1991–2005; the fields are 2024–2026.",
          "- Most Atlantic samples are from the cold North Atlantic; the weights help only where conditions "
          "overlap (effective n above).",
          "- Nitrate is a fixed climatology. SST is an input, so any climate association of a Route D record "
          "must be repeated without SST.", ""]
    out = cfg.path("reports") / "route_d.md"
    out.write_text("\n".join(L))
    return out
