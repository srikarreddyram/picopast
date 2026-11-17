"""Back-extension of the MOANA 8-day maps to 2012 from VIIRS (PRD section 6, steps 5-7).

The models learn "MOANA's estimate from what VIIRS sees" in the PACE era
(2024-2026) and are then applied to every VIIRS 8-day period from 2012. They
reproduce MOANA, so they can be no more right than MOANA is: whether MOANA is
right in this ocean is the trust gate's question, answered separately.

Built here:
- Route A (direct): VIIRS Rrs bands and ratios + OISST SST + VIIRS PAR + position
  + season -> group abundance; LightGBM, with 5% and 95% quantile models.
- Route C (niche, Flombaum-style): SST + PAR only.
- Baselines: seasonal climatology (training mean by calendar month in 2.5 degree
  cells) and ridge regression.
Not built yet: Route B (needs the VIIRS response tables, PRD item 10), the MLP,
mixed-layer depth (needs a Copernicus Marine account), the NOAA-20 copy.

Every group is reconstructed, but each output carries its trust-gate tier
(config ``reconstruction.tiers``): picoeukaryotes passed the gate and are a
relative index; Synechococcus and Prochlorococcus did not pass and are
exploratory only (pre-registration change log, 16 November 2026).

Targets: log10 cells/mL for Synechococcus and picoeukaryotes; Prochlorococcus in
units of 1e5 cells/mL, linear, because MOANA's Prochlorococcus regression is
linear. Training pixels follow the trust gate's pixel rule (MOANA flags 0, so
no negative Prochlorococcus) and must lie inside NASA's declared valid range
for its MOANA product for every group (config ``reconstruction.valid_range``):
outside it, the spectrum has pushed MOANA's regressions far beyond the Atlantic
data they were fitted to (Synechococcus down to 1e-38 cells/mL). Memory: one
8-day period (600 x 1200 pixels) at a time; the training table is a per-period
sample (~10k pixels each).
"""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config

GROUPS = ("pro", "syn", "peuk")
BANDS = (410, 443, 486, 551, 671)  # feature names follow SNPP's band centres (nm)
# Each sensor's own bands, in the same order. NOAA-20's sit a few nm away from SNPP's;
# each sensor gets its own models (PRD: the two reconstructions are built independently).
SENSOR_BANDS = {"snpp": BANDS, "noaa20": (411, 445, 489, 556, 667)}
PREFIX = {"snpp": "SNPP_VIIRS", "noaa20": "JPSS1_VIIRS"}
FEATURES = {
    "A": [f"rrs{b}" for b in BANDS] + ["r443_551", "r486_551", "r410_443", "r671_551", "mbr",
                                        "par", "sst", "lat", "lon", "doy_sin", "doy_cos"],
    "C": ["sst", "par"],
}
RIDGE_FEATURES = FEATURES["A"]
UNITS = {"pro": "1e5 cells/mL (linear)",
         "syn": "log10 cells/mL", "peuk": "log10 cells/mL"}
QUANTILES = (0.05, 0.95)
_PERIOD = re.compile(r"\.(\d{8})_(\d{8})\.L3m\.8D\.")



def usable_pixels(flags: np.ndarray, cells: dict[str, np.ndarray], valid_range: dict) -> np.ndarray:
    """Trust-gate pixel rule (MOANA flags 0) and every group inside NASA's declared valid range."""
    ok = flags == 0
    for g, (lo, hi) in valid_range.items():
        with np.errstate(invalid="ignore"):
            ok &= (cells[g] >= lo) & (cells[g] <= hi)
    return ok


def to_target(group: str, cells_ml):
    """Model units: 1e5 cells/mL for Prochlorococcus, log10 cells/mL otherwise (non-positive -> NaN)."""
    x = np.asarray(cells_ml, dtype=np.float64)
    if group == "pro":
        return x / 1e5
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(x > 0, np.log10(x), np.nan)


# ---------------------------------------------------------------- inputs


def viirs_periods(cfg: Config, sensor: str = "snpp") -> pd.DataFrame:
    """8-day periods with all five Rrs bands and PAR on disk: columns start, end."""
    raw, prefix = cfg.path("raw"), PREFIX[sensor]
    have: dict[str, set] = {}
    for path in (raw / f"viirs_{sensor}_par").glob(f"{prefix}.*.L3m.8D.PAR.par.4km.box.nc"):
        have.setdefault(_PERIOD.search(path.name).group(1), set()).add("par")
    for path in (raw / f"viirs_{sensor}_rrs").glob(f"{prefix}.*.L3m.8D.RRS.Rrs_*.4km.box.nc"):
        have.setdefault(_PERIOD.search(path.name).group(1), set()).add(path.name.split("Rrs_")[1][:3])
    need = {"par"} | {str(b) for b in SENSOR_BANDS[sensor]}
    rows = []
    for path in (raw / f"viirs_{sensor}_par").glob(f"{prefix}.*.L3m.8D.PAR.par.4km.box.nc"):
        m = _PERIOD.search(path.name)
        if need <= have.get(m.group(1), set()):
            rows.append((pd.Timestamp(m.group(1)), pd.Timestamp(m.group(2))))
    return pd.DataFrame(rows, columns=["start", "end"]).sort_values("start", ignore_index=True)


def load_viirs(cfg: Config, sensor: str, start: pd.Timestamp) -> dict[str, np.ndarray]:
    """The five Rrs bands and PAR of one 8-day period on the analysis grid (float32)."""
    raw, prefix, stamp = cfg.path("raw"), PREFIX[sensor], start.strftime("%Y%m%d")
    out = {}
    for name, b in zip(BANDS, SENSOR_BANDS[sensor]):
        (path,) = (raw / f"viirs_{sensor}_rrs").glob(f"{prefix}.{stamp}_*.L3m.8D.RRS.Rrs_{b}.4km.box.nc")
        with xr.open_dataset(path) as ds:
            out[f"rrs{name}"] = ds[f"Rrs_{b}"].values.astype(np.float32)
    (path,) = (raw / f"viirs_{sensor}_par").glob(f"{prefix}.{stamp}_*.L3m.8D.PAR.par.4km.box.nc")
    with xr.open_dataset(path) as ds:
        out["par"] = ds["par"].values.astype(np.float32)
    return out


class Oisst:
    """OISST v2.1 period means (degC) on the analysis grid, nearest 0.25 degree cell."""

    def __init__(self, cfg: Config, lat: np.ndarray, lon: np.ndarray):
        self.dir, self.lat, self.lon, self._years = cfg.path("raw") / "oisst", lat, lon, {}
        self._iy = self._ix = None

    def _year(self, year: int) -> xr.DataArray:
        if year not in self._years:
            if len(self._years) > 2:
                self._years.pop(min(self._years))
            self._years[year] = xr.open_dataset(self.dir / f"oisst_v2.1_box_{year}.nc")["sst"].load()
        return self._years[year]

    def mean(self, start: pd.Timestamp, end: pd.Timestamp) -> np.ndarray:
        parts = [self._year(y).sel(time=slice(start, end + pd.Timedelta(hours=23)))
                 for y in range(start.year, end.year + 1)]
        da = xr.concat(parts, "time") if len(parts) > 1 else parts[0]
        field = da.mean("time").values
        if self._iy is None:
            olat, olon = da["lat"].values, da["lon"].values
            self._iy = np.clip(np.round((self.lat - olat[0]) / (olat[1] - olat[0])).astype(int), 0, olat.size - 1)
            self._ix = np.clip(np.round((self.lon - olon[0]) / (olon[1] - olon[0])).astype(int), 0, olon.size - 1)
        return field[self._iy[:, None], self._ix[None, :]].astype(np.float32)


def _log_ratio(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((a > 0) & (b > 0), np.log10(a / b), np.nan).astype(np.float32)


def features(v: dict[str, np.ndarray], sst: np.ndarray, lat: np.ndarray, lon: np.ndarray,
             mid: pd.Timestamp) -> pd.DataFrame:
    """Feature table for the given pixels (1-D arrays of equal length)."""
    doy = 2 * np.pi * (mid.day_of_year - 0.5) / 365.25
    f = {f"rrs{b}": v[f"rrs{b}"] for b in BANDS}
    f.update(r443_551=_log_ratio(v["rrs443"], v["rrs551"]), r486_551=_log_ratio(v["rrs486"], v["rrs551"]),
             r410_443=_log_ratio(v["rrs410"], v["rrs443"]), r671_551=_log_ratio(v["rrs671"], v["rrs551"]),
             mbr=_log_ratio(np.maximum(v["rrs443"], v["rrs486"]), v["rrs551"]),  # OC3-style band ratio
             par=v["par"], sst=sst, lat=lat.astype(np.float32), lon=lon.astype(np.float32),
             doy_sin=np.full(lat.shape, np.sin(doy), np.float32), doy_cos=np.full(lat.shape, np.cos(doy), np.float32))
    return pd.DataFrame(f)


def _inputs_valid(v: dict[str, np.ndarray], sst: np.ndarray) -> np.ndarray:
    ok = np.isfinite(sst) & np.isfinite(v["par"])
    for b in BANDS:
        ok &= np.isfinite(v[f"rrs{b}"])
    return ok


def _grid_and_offshore(cfg: Config):
    from picopast.fallback_run import _offshore

    first = sorted((cfg.path("interim") / "moana" / "operational" / "8D").glob("moana_*.nc"))[0]
    with xr.open_dataset(first) as ds:
        lat, lon = ds["lat"].values, ds["lon"].values
        offshore = _offshore(cfg, ds)
    return lat, lon, offshore


# ---------------------------------------------------------------- training table


def training_table(cfg: Config, sensor: str = "snpp", mapping: str = "operational",
                   per_period: int = 10_000, seed: int = 0, refresh: bool = False) -> pd.DataFrame:
    """PACE-era pairs of VIIRS inputs and MOANA targets, sampled per 8-day period."""
    out = cfg.path("interim") / "reconstruct" / f"training_{sensor}_{mapping}.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    lat, lon, offshore = _grid_and_offshore(cfg)
    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    periods = viirs_periods(cfg, sensor).set_index("start")
    oisst, rng, parts = Oisst(cfg, lat, lon), np.random.default_rng(seed), []
    screen = {"retrieved": 0, "usable": 0, "negative_pro": 0, "outside_nasa_range": 0}
    for path in sorted((cfg.path("interim") / "moana" / mapping / "8D").glob("moana_*.nc")):
        with xr.open_dataset(path) as ds:
            start, end = pd.Timestamp(ds.attrs["period_start"]), pd.Timestamp(ds.attrs["period_end"])
            if start not in periods.index:
                continue
            cells = {g: ds[g].values[0] for g in GROUPS}
            usable = usable_pixels(ds["flags"].values[0], cells, cfg["reconstruction"]["valid_range"])
            targets = {g: np.where(usable, to_target(g, cells[g]), np.nan) for g in GROUPS}
        v = load_viirs(cfg, sensor, start)
        sst = oisst.mean(start, end)
        candidate = _inputs_valid(v, sst) & offshore & np.isfinite(cells["peuk"])  # MOANA ran here
        ok = candidate & usable
        screen["retrieved"] += int(candidate.sum())
        screen["usable"] += int(ok.sum())
        screen["negative_pro"] += int((candidate & (cells["pro"] < 0)).sum())
        screen["outside_nasa_range"] += int((candidate & ~usable).sum())
        idx = np.flatnonzero(ok.ravel())
        if idx.size > per_period:
            idx = np.sort(rng.choice(idx, per_period, replace=False))
        if not idx.size:
            continue
        mid = start + (end - start) / 2
        df = features({k: a.ravel()[idx] for k, a in v.items()}, sst.ravel()[idx],
                      lat2.ravel()[idx], lon2.ravel()[idx], mid)
        for g in GROUPS:
            df[g] = targets[g].ravel()[idx].astype(np.float32)
        df["time"], df["month"] = start, mid.month
        parts.append(df)
    table = pd.concat(parts, ignore_index=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(out, index=False)
    import json

    out.with_suffix(".screen.json").write_text(json.dumps(screen, indent=2))
    return table


# ---------------------------------------------------------------- models


class SeasonalClimatology:
    """Training mean by calendar month in ``cell_deg`` cells; the month's overall mean where a cell is empty."""

    def __init__(self, cell_deg: float = 2.5):
        self.cell_deg = cell_deg

    def _cell(self, df):
        return (np.floor(df["lat"].values / self.cell_deg) * 1000 + np.floor(df["lon"].values / self.cell_deg)).astype(int)

    def fit(self, df: pd.DataFrame, y: np.ndarray):
        d = pd.DataFrame({"m": df["month"].values, "c": self._cell(df), "y": y})
        self.by_cell = d.groupby(["m", "c"])["y"].mean()
        self.by_month = d.groupby("m")["y"].mean()
        self.overall = float(np.mean(y))
        return self

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        key = pd.MultiIndex.from_arrays([df["month"].values, self._cell(df)])
        out = self.by_cell.reindex(key).values
        fill = self.by_month.reindex(df["month"].values).fillna(self.overall).values
        return np.where(np.isfinite(out), out, fill)


def _lgbm(alpha: float | None = None, n_estimators: int = 400):
    from lightgbm import LGBMRegressor

    kw = dict(n_estimators=n_estimators, learning_rate=0.05, num_leaves=63, min_child_samples=200,
              subsample=0.8, subsample_freq=1, colsample_bytree=0.8, n_jobs=6, verbose=-1, random_state=0)
    return LGBMRegressor(objective="quantile", alpha=alpha, **kw) if alpha else LGBMRegressor(**kw)


def _ridge():
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), Ridge(alpha=1.0))


def _scores(y, pred, clim) -> dict:
    from scipy.stats import spearmanr

    rmse = float(np.sqrt(np.mean((pred - y) ** 2)))
    rmse_c = float(np.sqrt(np.mean((clim - y) ** 2)))
    return {"rmse": rmse, "skill": 1 - rmse / rmse_c, "spearman": float(spearmanr(pred, y)[0])}


def evaluate(table: pd.DataFrame, group: str, n_folds: int = 5, block_deg: float = 2.0,
             holdout_months: int = 6) -> pd.DataFrame:
    """Out-of-sample scores per model: pooled spatial-block folds, and the temporal hold-out.

    Everything (climatology, scaler, models) is fitted inside each training part.
    """
    from picopast.splits import spatial_block_folds, temporal_holdout

    d = table[np.isfinite(table[group])].reset_index(drop=True)
    y = d[group].values
    splits = {"spatial": list(spatial_block_folds(d["lat"].values, d["lon"].values, n_folds, block_deg)),
              "temporal": [temporal_holdout(d["time"].values, holdout_months)]}
    rows = []
    for kind, folds in splits.items():
        pred = {k: np.full(y.size, np.nan) for k in ("clim", "ridge", "C", "A", "A_lo", "A_hi")}
        for tr, te in folds:
            dtr, dte = d.iloc[tr], d.iloc[te]
            pred["clim"][te] = SeasonalClimatology().fit(dtr, y[tr]).predict(dte)
            pred["ridge"][te] = _ridge().fit(dtr[RIDGE_FEATURES], y[tr]).predict(dte[RIDGE_FEATURES])
            pred["C"][te] = _lgbm().fit(dtr[FEATURES["C"]], y[tr]).predict(dte[FEATURES["C"]])
            pred["A"][te] = _lgbm().fit(dtr[FEATURES["A"]], y[tr]).predict(dte[FEATURES["A"]])
            for q, k in zip(QUANTILES, ("A_lo", "A_hi")):
                pred[k][te] = _lgbm(q, 200).fit(dtr[FEATURES["A"]], y[tr]).predict(dte[FEATURES["A"]])
        m = np.isfinite(pred["clim"])
        for model in ("clim", "ridge", "C", "A"):
            s = _scores(y[m], pred[model][m], pred["clim"][m])
            cov = float(np.mean((y[m] >= pred["A_lo"][m]) & (y[m] <= pred["A_hi"][m]))) if model == "A" else np.nan
            rows.append({"group": group, "split": kind, "model": model, "n": int(m.sum()), **s,
                         "coverage_90": cov})
    return pd.DataFrame(rows)


def fit_final(table: pd.DataFrame, group: str, route: str) -> dict:
    """Median and 5%/95% models on all PACE-era data, with the feature ranges seen in training."""
    d = table[np.isfinite(table[group])]
    X, y = d[FEATURES[route]], d[group].values
    model = _lgbm().fit(X, y)
    lo, hi = (_lgbm(q, 200).fit(X, y) for q in QUANTILES)
    checked = [f for f in FEATURES[route] if f not in ("lat", "lon", "doy_sin", "doy_cos")]
    ranges = {f: (float(np.nanpercentile(d[f], 0.5)), float(np.nanpercentile(d[f], 99.5))) for f in checked}
    gain = pd.Series(model.booster_.feature_importance("gain"), index=FEATURES[route])
    return {"route": route, "median": model, "lo": lo, "hi": hi, "ranges": ranges,
            "importance": (gain / gain.sum()).sort_values(ascending=False),
            "index_mean": float(np.mean(y)), "index_sd": float(np.std(y))}


# ---------------------------------------------------------------- reconstruction


def reconstruct(cfg: Config, finals: dict, sensor: str = "snpp", tiers: dict | None = None,
                refresh: bool = False) -> list[Path]:
    """Apply the final models to every VIIRS 8-day period on disk; one NetCDF per period."""
    lat, lon, offshore = _grid_and_offshore(cfg)
    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    out_dir = cfg.path("products") / "reconstruction" / sensor / "8D"
    out_dir.mkdir(parents=True, exist_ok=True)
    oisst, written = Oisst(cfg, lat, lon), []
    for p in viirs_periods(cfg, sensor).itertuples():
        out = out_dir / f"recon_{p.start:%Y%m%d}.nc"
        if out.exists() and not refresh:
            written.append(out)
            continue
        v = load_viirs(cfg, sensor, p.start)
        sst = oisst.mean(p.start, p.end)
        ok = _inputs_valid(v, sst) & offshore
        idx = np.flatnonzero(ok.ravel())
        mid = p.start + (p.end - p.start) / 2
        X = features({k: a.ravel()[idx] for k, a in v.items()}, sst.ravel()[idx],
                     lat2.ravel()[idx], lon2.ravel()[idx], mid)
        data = {}
        for g, fin in finals.items():
            Xg = X[FEATURES[fin["route"]]]
            for name, model in (("", fin["median"]), ("_lo", fin["lo"]), ("_hi", fin["hi"])):
                grid = np.full(ok.size, np.nan, np.float32)
                if idx.size:
                    grid[idx] = model.predict(Xg)
                data[g + name] = (("lat", "lon"), grid.reshape(ok.shape))
            n_out = np.zeros(ok.size, np.uint8)
            for f, (a, b) in fin["ranges"].items():
                vals = X[f].values
                n_out[idx] += ((vals < a) | (vals > b)).astype(np.uint8)
            data[f"{g}_out_of_range"] = (("lat", "lon"), n_out.reshape(ok.shape))
        ds = xr.Dataset(data, coords={"lat": lat, "lon": lon})
        ds = ds.expand_dims(time=[p.start])
        for g, fin in finals.items():
            ds[g].attrs.update(units=UNITS[g], route=fin["route"], tier=(tiers or {}).get(g, ""),
                               relative_index=f"(value - {fin['index_mean']:.4f}) / {fin['index_sd']:.4f}",
                               interval="5% and 95% quantile models (_lo, _hi)")
            ds[f"{g}_out_of_range"].attrs["meaning"] = "number of inputs outside the 0.5-99.5% range seen in training"
        ds.attrs.update(period_start=f"{p.start:%Y-%m-%d}", period_end=f"{p.end:%Y-%m-%d}", sensor=sensor,
                        created=date.today().isoformat(),
                        source="PicoPast back-extension of MOANA (operational mapping) from VIIRS; see "
                               "docs/reports/reconstruction.md")
        enc = {k: {"zlib": True, "complevel": 4} for k in ds.data_vars}
        part = out.with_name(out.name + ".part")
        ds.to_netcdf(part, encoding=enc)
        part.rename(out)
        written.append(out)
    return written


def regional_series(cfg: Config, sensor: str = "snpp") -> pd.DataFrame:
    """Per-period regional means of the reconstruction (and valid-pixel share), for plots and RQ5."""
    paths = sorted((cfg.path("products") / "reconstruction" / sensor / "8D").glob("recon_*.nc"))
    rows = []
    for path in paths:
        with xr.open_dataset(path) as ds:
            for name, r in [("study_box", None)] + list(cfg["regions"].items()):
                sub = ds if r is None else ds.sel(lat=slice(r["lat_max"], r["lat_min"]),
                                                  lon=slice(r["lon_min"], r["lon_max"]))
                row = {"time": pd.Timestamp(ds.attrs["period_start"]), "region": name,
                       "valid_share": float(np.isfinite(sub["peuk"]).mean())}
                for g in GROUPS:
                    if g in sub:
                        vals = sub[g].values
                        row[g] = float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan
                        row[f"{g}_n"] = int(np.isfinite(vals).sum())
                        row[f"{g}_n_out"] = int((sub[f"{g}_out_of_range"].values[np.isfinite(vals)] > 0).sum())
                rows.append(row)
    return pd.DataFrame(rows)
