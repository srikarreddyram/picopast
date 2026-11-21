"""Route D back-extension: Synechococcus and Prochlorococcus from D2, every 8-day period 2012-2026.

The pre-registered decision (docs/route_d_preregistration.md): both groups
passed check C with D2, so both are reconstructed from D2. Synechococcus is a relative index;
Prochlorococcus failed check E (HPLC pigments, docs/reports/check_e.md) and is exploratory.
Per 8-day period the inputs are OISST (period mean), VIIRS chlorophyll and PAR
(8-day composites), WOA23 nitrate (month of the period's middle day) and
season. A pixel is *reliable* when every input lies inside the 1-99%
importance-weighted range of the training samples (the conditions the adjusted
model effectively learned from).

D2 learned from climatological chlorophyll and PAR; here it sees 8-day values,
which vary more. The reliable flag keeps the record to conditions inside the
training range. An SST-free D2 is run alongside, for the climate analyses.

Outputs: one NetCDF per period under data/products/reconstruction/<sensor>_route_d/8D,
0.25° cells for the climate maps, and regional means per period.
"""

from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.route_d import (
    FEATURES,
    _box_rows,
    _fit,
    box_fields,
    effective_n,
    importance_weights,
    season,
    training_table,
)

GROUPS = ("syn", "pro")
NAMES = {"syn": "Synechococcus", "pro": "Prochlorococcus"}
VARIANTS = {"": FEATURES, "_nosst": [f for f in FEATURES if f != "sst"]}
RANGE_INPUTS = ("sst", "logchl", "par", "no3")
PREFIX = {"snpp": "SNPP_VIIRS", "noaa20": "JPSS1_VIIRS"}
_PERIOD = re.compile(r"\.(\d{8})_(\d{8})\.L3m\.8D\.")


def weighted_quantile(x, w, q) -> np.ndarray:
    """Weighted quantiles with mid-point cumulative weights; zero-weight samples are ignored."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    keep = w > 0
    x, w = x[keep], w[keep]
    order = np.argsort(x)
    x, w = x[order], w[order]
    cw = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(np.asarray(q), cw, x)


def final_models(cfg: Config) -> tuple[dict, dict, dict]:
    """D2 (and D2 without SST) on all non-Indian training stations; weighted input ranges; summary."""
    t = training_table(cfg)
    train = t[(t["basin"] != "indian") & t[FEATURES].notna().all(axis=1)].reset_index(drop=True)
    fields = box_fields(cfg)
    Xb, _ = _box_rows(fields)
    rng = np.random.default_rng(0)  # same box sample as the test
    box_sample = Xb.iloc[rng.choice(len(Xb), min(50_000, len(Xb)), replace=False)].to_numpy()
    models, ranges, info = {}, {}, {}
    for g in GROUPS:
        d = train[train[g].notna()].reset_index(drop=True)
        w = importance_weights(d[FEATURES].to_numpy(), box_sample)
        y = np.log10(d[g].to_numpy())
        for suffix, feats in VARIANTS.items():
            models[g + suffix] = (_fit(d[feats], y, w), feats)
        ranges[g] = {f: weighted_quantile(d[f], w, [0.01, 0.99]).tolist() for f in RANGE_INPUTS}
        info[g] = {"n_train": int(len(d)), "ess": effective_n(w)}
    return models, ranges, info


def chl_par_periods(cfg: Config, sensor: str) -> pd.DataFrame:
    """8-day periods with both chlorophyll and PAR on disk: start, end, chl path, par path."""
    raw, prefix = cfg.path("raw"), PREFIX[sensor]
    chl = {_PERIOD.search(p.name).group(1): p
           for p in (raw / f"viirs_{sensor}_chl").glob(f"{prefix}.*.L3m.8D.CHL.chlor_a.4km.box.nc")}
    rows = []
    for p in (raw / f"viirs_{sensor}_par").glob(f"{prefix}.*.L3m.8D.PAR.par.4km.box.nc"):
        k = _PERIOD.search(p.name)
        if k.group(1) in chl:
            rows.append((pd.Timestamp(k.group(1)), pd.Timestamp(k.group(2)), chl[k.group(1)], p))
    return pd.DataFrame(rows, columns=["start", "end", "chl", "par"]).sort_values("start", ignore_index=True)


def _region_masks(cfg: Config, lat2: np.ndarray, lon2: np.ndarray) -> dict[str, np.ndarray]:
    boxes = {"study_box": None} | dict(cfg["regions"].items())
    return {k: np.ones_like(lat2, bool) if r is None else
            (lat2 >= r["lat_min"]) & (lat2 <= r["lat_max"]) & (lon2 >= r["lon_min"]) & (lon2 <= r["lon_max"])
            for k, r in boxes.items()}


def all_valid(cfg: Config, sensor: str = "snpp") -> tuple[pd.DataFrame, xr.Dataset]:
    """Regional means and 0.25° cells from *every* valid pixel, ignoring the reliable flag.

    A sensitivity check: the flag drops pixels whose inputs leave the training range, and in this box that is
    mostly the warmest water (Prochlorococcus training SST ends near 29.9 °C), so composites over hot periods
    compare different pixels. Outside the range LightGBM holds its edge value, so this version is conservative.
    Read from the saved 8-day files; cached next to the reliable-pixel versions.
    """
    from picopast.climate_maps import CELL, coarsen

    work = cfg.path("interim") / "route_d"
    csv, nc = work / f"regional_all_{sensor}.csv", cfg.path("interim") / "climate" / f"cells_route_d_all_{sensor}.nc"
    if csv.exists() and nc.exists():
        return pd.read_csv(csv, parse_dates=["time"]), xr.open_dataset(nc).load()
    files = sorted((cfg.path("products") / "reconstruction" / f"{sensor}_route_d" / "8D").glob("recon_d_*.nc"))
    keys = [g + s for g in GROUPS for s in VARIANTS]
    rows, cells, times, masks = [], {k: [] for k in keys}, [], None
    for f in files:
        with xr.open_dataset(f) as ds:
            if masks is None:
                lat, lon = ds["lat"].values, ds["lon"].values
                masks = _region_masks(cfg, *np.meshgrid(lat, lon, indexing="ij"))
            t = pd.Timestamp(ds["time"].values[0])
            for k in keys:
                v = ds[k].values[0]
                good = np.isfinite(v)
                cells[k].append(coarsen(v, good))
                for name, m in masks.items():
                    sel = good & m
                    rows.append({"time": t, "region": name, "series": k,
                                 "value": float(v[sel].mean()) if sel.any() else np.nan, "n": int(sel.sum())})
        times.append(t)
    clat = lat.reshape(-1, CELL).mean(1)
    clon = lon.reshape(-1, CELL).mean(1)
    order = np.argsort(clat)
    cube = xr.Dataset({k: (("time", "lat", "lon"), np.stack(v)[:, order]) for k, v in cells.items()},
                      coords={"time": times, "lat": np.round(clat[order], 4), "lon": np.round(clon, 4)})
    cube.to_netcdf(nc)
    table = pd.DataFrame(rows)
    table.to_csv(csv, index=False)
    return table, cube


def reconstruct(cfg: Config, sensor: str = "snpp", refresh: bool = False) -> dict:
    """Write the per-period maps, the 0.25° cell cube and the regional series. Returns the metadata."""
    from picopast.climate_maps import CELL, coarsen
    from picopast.reconstruct import Oisst

    work = cfg.path("interim") / "route_d"
    meta_path = work / f"recon_{sensor}_meta.json"
    if meta_path.exists() and not refresh:
        return json.loads(meta_path.read_text())
    out_dir = cfg.path("products") / "reconstruction" / f"{sensor}_route_d" / "8D"
    out_dir.mkdir(parents=True, exist_ok=True)
    models, ranges, info = final_models(cfg)
    fields = box_fields(cfg)
    ds_f, offshore = fields["ds"], fields["offshore"]
    lat, lon = ds_f["lat"].values, ds_f["lon"].values
    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    no3_by_month = {int(m): ds_f["no3"].sel(month=m).values for m in ds_f["month"].values}
    oisst = Oisst(cfg, lat, lon)
    masks = _region_masks(cfg, lat2, lon2)
    keys = list(models)
    cells, times, rows = {k: [] for k in keys}, [], []
    total = {k: 0.0 for k in keys}
    total2 = {k: 0.0 for k in keys}
    count = {k: 0 for k in keys}
    for p in chl_par_periods(cfg, sensor).itertuples():
        with xr.open_dataset(p.chl) as ds:
            chl = ds["chlor_a"].values.astype(np.float32)
        with xr.open_dataset(p.par) as ds:
            par = ds["par"].values.astype(np.float32)
        with np.errstate(divide="ignore", invalid="ignore"):
            logchl = np.where(chl > 0, np.log10(chl), np.nan).astype(np.float32)
        sst = oisst.mean(p.start, p.end)
        month = (p.start + (p.end - p.start) / 2).month
        no3 = no3_by_month[month]
        grid = {"sst": sst, "logchl": logchl, "par": par, "no3": no3}
        ok = offshore & np.all([np.isfinite(v) for v in grid.values()], axis=0)
        idx = np.flatnonzero(ok.ravel())
        s_sin, s_cos = season([month], [10.0])
        X = pd.DataFrame({f: grid[f].ravel()[idx] for f in RANGE_INPUTS})
        X["season_sin"], X["season_cos"] = s_sin[0], s_cos[0]
        data = {}
        for g in GROUPS:
            rel = np.ones(idx.size, bool)
            for f, (a, b) in ranges[g].items():
                rel &= (X[f].values >= a) & (X[f].values <= b)
            r_map = np.zeros(ok.size, np.uint8)
            r_map[idx[rel]] = 1
            data[f"{g}_reliable"] = (("lat", "lon"), r_map.reshape(ok.shape))
            for suffix in VARIANTS:
                model, feats = models[g + suffix]
                v = np.full(ok.size, np.nan, np.float32)
                if idx.size:
                    v[idx] = model.predict(X[feats])
                v = v.reshape(ok.shape)
                data[g + suffix] = (("lat", "lon"), v)
                good = (r_map.reshape(ok.shape) == 1) & np.isfinite(v)
                total[g + suffix] += float(v[good].sum())
                total2[g + suffix] += float((v[good].astype(np.float64) ** 2).sum())
                count[g + suffix] += int(good.sum())
                cells[g + suffix].append(coarsen(v, good))
                for name, m in masks.items():
                    sel = good & m
                    rows.append({"time": p.start, "region": name, "series": g + suffix,
                                 "value": float(v[sel].mean()) if sel.any() else np.nan,
                                 "n": int(sel.sum()), "n_valid": int((np.isfinite(v) & m).sum())})
        times.append(p.start)
        out = xr.Dataset({k: (("time",) + dims, arr[None]) for k, (dims, arr) in data.items()},
                         coords={"time": [p.start], "lat": lat, "lon": lon})
        for g in GROUPS:
            out[g].attrs = {"units": "log10 cells/mL (Route D, D2); use as a relative index",
                            "long_name": f"{NAMES[g]}, Route D"}
            out[f"{g}_reliable"].attrs["meaning"] = "1 = every input inside the 1-99% weighted training range"
        out.attrs = {"period_start": str(p.start.date()), "period_end": str(p.end.date()), "sensor": sensor,
                     "created": date.today().isoformat(),
                     "source": "PicoPast Route D (globally ship-trained, Indian-Ocean-adjusted); "
                               "see docs/reports/route_d.md"}
        enc = {k: {"zlib": True, "complevel": 4} for k in out.data_vars}
        out.to_netcdf(out_dir / f"recon_d_{p.start:%Y%m%d}.nc", encoding=enc)
    clat = lat.reshape(-1, CELL).mean(1)
    clon = lon.reshape(-1, CELL).mean(1)
    order = np.argsort(clat)
    cube = xr.Dataset({k: (("time", "lat", "lon"), np.stack(v)[:, order]) for k, v in cells.items()},
                      coords={"time": times, "lat": np.round(clat[order], 4), "lon": np.round(clon, 4)})
    cube.to_netcdf(cfg.path("interim") / "climate" / f"cells_route_d_{sensor}.nc")
    pd.DataFrame(rows).to_csv(work / f"regional_{sensor}.csv", index=False)
    meta = {"sensor": sensor, "n_periods": len(times), "ranges": ranges, "train": info,
            "index": {k: {"mean": total[k] / count[k], "sd": float(np.sqrt(total2[k] / count[k]
                                                                       - (total[k] / count[k]) ** 2))}
                      for k in keys if count[k]}}
    meta_path.write_text(json.dumps(meta, indent=1))
    return meta


# ---------------------------------------------------------------- analysis and report

#: Tier of each Route D record. Prochlorococcus passed check C (1991-2005 counts) but failed the
#: pre-registered check E (2007-2025 HPLC divinyl chlorophyll a, docs/reports/check_e.md), so the
#: binding rule makes it exploratory.
TIERS = {"syn": "relative index (Route D): passed check C",
         "pro": "exploratory (Route D): passed check C, failed check E"}
TIER_SHORT = {"syn": "relative index (Route D: passed check C)",
              "pro": "EXPLORATORY (Route D: passed check C, failed check E)"}

REGION_LABELS = {"study_box": "Study box", "bay_of_bengal": "Bay of Bengal",
                 "arabian_sea_central_east": "Arabian Sea C/E"}
COLORS = {"study_box": "#2a78d6", "bay_of_bengal": "#eb6834", "arabian_sea_central_east": "#1baf7a"}


def index_series(cfg: Config, sensor: str = "snpp", pixels: str = "reliable") -> pd.DataFrame:
    """Regional means as a relative index (SD units of the whole record): time, region, series columns.

    pixels="all" uses every valid pixel (see all_valid), on the same scale as the reliable-pixel index.
    """
    meta = json.loads((cfg.path("interim") / "route_d" / f"recon_{sensor}_meta.json").read_text())
    if pixels == "all":
        r = all_valid(cfg, sensor)[0]
    else:
        r = pd.read_csv(cfg.path("interim") / "route_d" / f"regional_{sensor}.csv", parse_dates=["time"])
    r["z"] = [(v - meta["index"][s]["mean"]) / meta["index"][s]["sd"] for v, s in zip(r["value"], r["series"])]
    return r.pivot_table(index=["time", "region"], columns="series", values="z").reset_index()


def _anomaly(z: pd.DataFrame, cols) -> pd.DataFrame:
    s = z[z["region"].isin(REGION_LABELS)].copy()
    slot = (s["time"].dt.day_of_year - 1) // 8
    out = s[["time", "region"]].copy()
    for c in cols:
        out[c] = s[c] - s[c].groupby([s["region"], slot]).transform("mean")
    return out


def climate_tables(cfg: Config, z: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Trends, lagged ENSO/IOD correlations and heatwave composites, for D2 and for D2 without SST."""
    from picopast import climate

    idx = climate.load_indices(cfg)
    share = pd.read_csv(cfg.path("interim") / "climate" / "heatwave_share_8D.csv", index_col=0, parse_dates=True)
    out = {}
    for label, cols in (("main", {"syn": "syn", "pro": "pro"}), ("no SST", {"syn": "syn_nosst", "pro": "pro_nosst"})):
        a = _anomaly(z, cols.values()).rename(columns={v: k for k, v in cols.items()}).assign(peuk=np.nan)
        out[f"trends {label}"] = climate.trends(cfg, a)
        out[f"lags {label}"] = climate.lag_correlations(climate.monthly(a), idx)
        out[f"heatwave {label}"] = climate.heatwave_composite(cfg, a, share)
    return out


def map_composites(cfg: Config, sensor: str = "snpp", pixels: str = "reliable") -> dict:
    """Per-cell ENSO, IOD and heatwave composites of the Route D records (SD units), as in climate_maps."""
    from picopast import climate_maps as cm
    from picopast.climate import load_indices

    meta = json.loads((cfg.path("interim") / "route_d" / f"recon_{sensor}_meta.json").read_text())
    if pixels == "all":
        cube = all_valid(cfg, sensor)[1]
    else:
        cube = xr.open_dataset(cfg.path("interim") / "climate" / f"cells_route_d_{sensor}.nc").load()
    time = pd.DatetimeIndex(cube["time"].values)
    hw = cm.heatwave_fraction(cfg, time)
    idx = load_indices(cfg)
    classes = {e: (cm.heatwave_classes(hw.values, time) if name is None
                   else cm.index_classes(idx[name].dropna(), time, thr)) for e, (_, name, thr) in cm.EVENTS.items()}
    res = {"lat": cube["lat"].values, "lon": cube["lon"].values, "time": time,
           "ever_valid": np.isfinite(cube["syn"].values).any(axis=0)}
    years = time.year.to_numpy()
    for k in cube.data_vars:
        zc = (cube[k].values - meta["index"][k]["mean"]) / meta["index"][k]["sd"]
        anom = cm.slot_anomaly(zc, time)
        for e, (a, b) in classes.items():
            res[(e, k)] = cm.composite(anom, years, a, b)
    return res


def plot_series(z: pd.DataFrame, out: Path, title: str) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, surface = "#2b2b29", "#6b6a63", "#fcfcfb"
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, facecolor=surface)
    for ax, g in zip(axes, GROUPS):
        ax.set_facecolor(surface)
        ends = []
        for region, label in REGION_LABELS.items():
            m = z[z["region"] == region].set_index("time")[g].resample("MS").mean()
            ax.plot(m.index, m.values, color=COLORS[region], lw=1.6)
            last = m.dropna()
            ends.append([last.values[-1], last.index[-1], label])
        lo, hi = ax.get_ylim()
        ends.sort(key=lambda e: e[0])
        for i in range(1, len(ends)):
            ends[i][0] = max(ends[i][0], ends[i - 1][0] + 0.08 * (hi - lo))
        for y, x, label in ends:
            ax.annotate(label, (x, y), xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                        color=ink, annotation_clip=False)
        ax.set_title(f"{NAMES[g]} — {TIER_SHORT[g]}", loc="left", fontsize=10, color=ink)
        ax.axhline(0, color=muted, lw=0.6)
        ax.grid(axis="y", color="#e4e3dc", lw=0.6)
        ax.set_ylabel("index (SD)", fontsize=8, color=muted)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(muted)
        ax.tick_params(colors=muted, labelsize=8)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=ink)
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    fig.savefig(out, dpi=150, facecolor=surface)
    plt.close(fig)
    return out


def _monthly(z: pd.DataFrame, region: str, col: str) -> pd.Series:
    return z[z["region"] == region].set_index("time")[col].resample("MS").mean()


def report(cfg: Config, refresh: bool = False) -> Path:
    from scipy.stats import theilslopes

    from picopast import climate
    from picopast import climate_maps as cm
    from picopast.climate import REGIONS

    meta = reconstruct(cfg, "snpp", refresh)
    have_n20 = not chl_par_periods(cfg, "noaa20").empty
    if have_n20:
        reconstruct(cfg, "noaa20", refresh)
    z = index_series(cfg, "snpp")
    res = json.loads((cfg.path("interim") / "route_d" / "results.json").read_text())
    figs = cfg.path("reports") / "figures"
    plot_series(z, figs / "route_d_series.png",
                "Cyanobacteria 2012–2026 from Route D (VIIRS SNPP), monthly regional means, reliable pixels")
    regional = pd.read_csv(cfg.path("interim") / "route_d" / "regional_snpp.csv", parse_dates=["time"])
    L = ["# Route D back-extension: Synechococcus and Prochlorococcus, 2012–2026", "",
         f"Generated {date.today().isoformat()} by `picopast route-d-reconstruct`. The pre-registered decision "
         "(`docs/route_d_preregistration.md`): both groups passed check C with D2 (`route_d.md`), so both are "
         "reconstructed from D2. *Synechococcus* is a **relative index**. *Prochlorococcus* then failed the "
         "pre-registered check E against 2007–2025 HPLC divinyl chlorophyll a (`check_e.md`) and is "
         "**exploratory**. Per 8-day period: OISST, VIIRS chlorophyll and PAR, WOA23 "
         "nitrate for the month, season. Reliable pixels: every input inside the 1–99% importance-weighted range "
         "of the training samples. Index = SD units of the whole record.", "",
         "| group | tier | check C ρ, D2 | check C ρ, MOANA | check E (record vs HPLC pigment) | "
         "training stations (effective n) | 8-day periods | reliable share of valid pixels |",
         "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    e_path = cfg.path("interim") / "check_e" / "results.json"
    e_res = json.loads(e_path.read_text()) if e_path.exists() else {}
    for g in GROUPS:
        r = res["groups"][g]
        q = regional[(regional["region"] == "study_box") & (regional["series"] == g)]
        e = e_res.get(g, {}).get("E2", {})
        e_txt = (f"{e['record']['rho']:+.2f} [{e['record']['lo']:+.2f}, {e['record']['hi']:+.2f}], "
                 f"*{e['verdict']}*" if "verdict" in e else "no verdict (zeaxanthin is shared; no new counts)")
        L.append(f"| {NAMES[g]} | {TIERS[g].split(' (')[0]} | {r['D2']['stats']['rho']:.2f} | "
                 f"{r['moana_stats']['rho']:.2f} | {e_txt} | "
                 f"{meta['train'][g]['n_train']:,} ({meta['train'][g]['ess']:.0f}) | {meta['n_periods']} | "
                 f"{q['n'].sum() / q['n_valid'].sum():.0%} |")
    L += ["", "![series](figures/route_d_series.png)", "",
          "## Against the earlier exploratory records (MOANA learned from VIIRS, Route A)", "",
          "Correlation of monthly anomalies (seasonal cycle removed) between the Route D record and the earlier "
          "exploratory one, same sensor, same regions. Low values mean the two disagree on *when* the groups "
          "were high or low, not only on the level.", "",
          "| group | " + " | ".join(REGION_LABELS.values()) + " |", "| --- | --- | --- | --- |"]
    old = pd.read_csv(cfg.path("products") / "reconstruction" / "snpp_8D_regional_means.csv", parse_dates=["time"])
    for g in GROUPS:
        cells = []
        for region in REGION_LABELS:
            a = _monthly(z, region, g)
            b = old[old["region"] == region].set_index("time")[f"{g}_reliable"].resample("MS").mean()
            both = pd.concat([a, b], axis=1).dropna()
            both = both - both.groupby(both.index.month).transform("mean")
            cells.append(f"{both.corr().iloc[0, 1]:.2f}")
        L.append(f"| {NAMES[g]} | " + " | ".join(cells) + " |")
    if have_n20:
        z2 = index_series(cfg, "noaa20")
        L += ["", "## SNPP vs NOAA-20 (same model, the other satellite's chlorophyll and PAR)", "",
              "Monthly study-box index from each satellite over the overlap. Gap = mean SNPP − NOAA-20; drift = "
              "Theil–Sen slope of that difference per decade (95% interval). For comparison, the MOANA-based "
              "exploratory records drifted +0.13 (*Synechococcus*) and +0.25 (*Prochlorococcus*) SD per decade. "
              "Both copies share OISST and the nitrate climatology, so part of the agreement is built in; the "
              "models without SST see only what differs between the satellites (chlorophyll and PAR) and are the "
              "fairer test.", "",
              "| group | model | months | r | gap (SD) | drift per decade [95%] |",
              "| --- | --- | --- | --- | --- | --- |"]
        for g in GROUPS:
            for suffix, label in (("", "main"), ("_nosst", "without SST")):
                both = pd.concat([_monthly(z, "study_box", g + suffix), _monthly(z2, "study_box", g + suffix)],
                                 axis=1).dropna()
                diff = both.iloc[:, 0] - both.iloc[:, 1]
                yrs = (diff.index - diff.index[0]).days / 365.25
                s, _, lo, hi = theilslopes(diff.values, yrs)
                L.append(f"| {NAMES[g]} | {label} | {len(both)} | {both.corr().iloc[0, 1]:.2f} | {diff.mean():+.2f} | "
                         f"{10 * s:+.2f} [{10 * lo:+.2f}, {10 * hi:+.2f}] |")
    t = climate_tables(cfg, z)
    for name in ("trends main", "trends no SST", "lags main", "lags no SST", "heatwave main", "heatwave no SST"):
        t[name].to_csv(cfg.path("interim") / "route_d" / f"climate_{name.replace(' ', '_')}.csv", index=False)
    L += ["", "## Climate associations (RQ5 for the cyanobacteria)", "",
          "Same methods as `climate.md`. SST is a Route D input and heatwaves and ENSO are SST signals, so every "
          "result is shown next to the same model trained without SST.", "",
          "### Trends (core seasons, annual means)", "",
          "| region | group | years | Sen's slope per decade | p (FDR) | without SST: slope, p (FDR) |",
          "| --- | --- | --- | --- | --- | --- |"]
    tm, tn = t["trends main"], t["trends no SST"]
    for r in tm[tm["group"].isin(GROUPS)].itertuples():
        q = tn[(tn["region"] == r.region) & (tn["group"] == r.group)].iloc[0]
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {r.years} | {r.sen_per_decade:+.2f} | {r.p_fdr:.3f} | "
                 f"{q.sen_per_decade:+.2f}, {q.p_fdr:.3f} |")
    L += ["", "### Strongest lagged correlation with ONI and DMI", "",
          "| region | group | index | lag (months) | r | p (FDR) | without SST: r at the same lag |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    lm, ln = t["lags main"], t["lags no SST"]
    lm = lm[lm["group"].isin(GROUPS)]
    best = lm.loc[lm.groupby(["region", "group", "index"])["r"].apply(lambda s: s.abs().idxmax())]
    for r in best.sort_values(["group", "region", "index"]).itertuples():
        q = ln[(ln["region"] == r.region) & (ln["group"] == r.group) & (ln["index"] == r.index)
               & (ln["lag_months"] == r.lag_months)].iloc[0]
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {r.index} | {r.lag_months} | {r.r:+.2f} | "
                 f"{r.p_fdr:.3f} | {q.r:+.2f} |")
    share = pd.read_csv(cfg.path("interim") / "climate" / "heatwave_share_8D.csv", index_col=0, parse_dates=True)
    z_all = index_series(cfg, "snpp", pixels="all")
    a_all = _anomaly(z_all, GROUPS).assign(peuk=np.nan)
    h_all = climate.heatwave_composite(cfg, a_all, share)
    h_all.to_csv(cfg.path("interim") / "route_d" / "climate_heatwave_all_pixels.csv", index=False)
    rel = regional[regional["series"].isin(GROUPS)].merge(share.stack().rename("hw").reset_index()
                                                          .set_axis(["time", "region", "hw"], axis=1))
    rel["share"] = rel["n"] / rel["n_valid"]
    L += ["", "### Marine heatwaves", "",
          "Reliable pixels thin out during heatwaves because the warmest water leaves the training range "
          "(*Prochlorococcus* training SST ends at "
          f"{meta['ranges']['pro']['sst'][1]:.1f} °C, *Synechococcus* at {meta['ranges']['syn']['sst'][1]:.1f} °C). "
          "A composite over reliable pixels then compares cooler water during heatwaves with all water otherwise. "
          "The last column repeats it with every valid pixel (the model holds its edge value beyond the range, so "
          "this is conservative). Read a heatwave result only where both agree in sign.", "",
          "| region | group | reliable pixels, heatwave / normal | heatwave − normal [95%] | without SST | "
          "all valid pixels |", "| --- | --- | --- | --- | --- | --- |"]
    hm, hn = t["heatwave main"], t["heatwave no SST"]
    for r in hm[hm["group"].isin(GROUPS)].itertuples():
        q = hn[(hn["region"] == r.region) & (hn["group"] == r.group)].iloc[0]
        w = h_all[(h_all["region"] == r.region) & (h_all["group"] == r.group)].iloc[0]
        if pd.isna(getattr(r, "difference", np.nan)):
            continue
        s = rel[(rel["region"] == r.region) & (rel["series"] == r.group)]
        hot = s[s["hw"] >= 0.5]
        normal = s[(s["hw"] < 0.1) & s["time"].dt.month.isin(hot["time"].dt.month.unique())]
        L.append(f"| {REGIONS[r.region]} | {NAMES[r.group]} | {hot['share'].mean():.0%} / {normal['share'].mean():.0%} | "
                 f"{r.difference:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}] | {q.difference:+.2f} [{q.lo:+.2f}, {q.hi:+.2f}] | "
                 f"{w.difference:+.2f} [{w.lo:+.2f}, {w.hi:+.2f}] |")
    maps = map_composites(cfg, "snpp")
    maps_all = map_composites(cfg, "snpp", pixels="all")
    tiers = TIERS
    cm.plot(maps, ("syn", "syn_nosst", "pro", "pro_nosst"), tiers, figs / "route_d_climate_maps.png",
            "Cyanobacteria (Route D) during climate events, VIIRS SNPP 2012–2026, 0.25° cells (dots: 95% interval "
            "excludes zero); each group with and without SST")
    L += ["", "### Maps per 0.25° cell", "", "Same composites as `climate_maps.md`.", "",
          "![maps](figures/route_d_climate_maps.png)", "",
          "| group | event | mean difference (SD) | clear lower / clear higher | without SST: mean, pattern r | "
          "all valid pixels: mean, clear lower / higher, pattern r |",
          "| --- | --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        for e, (label, _, _) in cm.EVENTS.items():
            a, b, c = maps[(e, g)], maps[(e, g + "_nosst")], maps_all[(e, g)]
            ok, okc = np.isfinite(a["diff"]), np.isfinite(c["diff"])
            L.append(f"| {NAMES[g]} | {label} | {np.nanmean(a['diff']):+.2f} | {(a['hi'][ok] < 0).mean():.0%} / "
                     f"{(a['lo'][ok] > 0).mean():.0%} | {np.nanmean(b['diff']):+.2f}, "
                     f"{cm._pattern_r(a['diff'], b['diff']):.2f} | {np.nanmean(c['diff']):+.2f}, "
                     f"{(c['hi'][okc] < 0).mean():.0%} / {(c['lo'][okc] > 0).mean():.0%}, "
                     f"{cm._pattern_r(a['diff'], c['diff']):.2f} |")
    L += ["", "## Limits", "",
          "- D2 learned from climatological chlorophyll and PAR and 1987–2011 counts; here it reads 8-day values. "
          "The reliable flag keeps it inside the weighted training range, but 8-day extremes are rarer there.",
          "- Its check C skill is about that of chlorophyll alone (`route_d.md`): the record says where and when "
          "each group is high or low in the way the global data relate them to chlorophyll, SST, light and "
          "nitrate, not more.",
          "- Nitrate is a fixed climatology; year-to-year change enters through SST, chlorophyll and PAR only.", ""]
    out = cfg.path("reports") / "route_d_reconstruction.md"
    out.write_text("\n".join(L))
    return out
