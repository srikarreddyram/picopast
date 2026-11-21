"""Pixel-level climate composites: where the record moves with ENSO, the IOD and marine heatwaves (RQ5 follow-up).

The regional analysis (climate.py) averages whole regions. Here the same
question is asked per 0.25° cell, the OISST grid, which the 1/24°
reconstruction grid divides exactly (6 x 6 pixels per cell):

- Cell value: mean of the reliable pixels (every model input inside the range
  seen in training) as a relative index in SD units; a cell needs at least
  MIN_PIXELS reliable pixels in the period, else it is missing.
- Anomaly: cell value minus that cell's mean for the same 8-day slot of the
  year over the whole record.
- ENSO / IOD: periods whose month had ONI (or DMI) two months earlier at or
  above +0.5 (+0.4) against periods at or below -0.5 (-0.4). The two-month lead
  is the strongest lag in the regional analysis for both indices.
- Marine heatwave: periods when at least half of the cell's own days were in a
  heatwave (Hobday, OISST, 1991-2020 baseline) against periods of the same
  calendar months with under 10%, both in that cell.
- Interval: 95% from resampling whole years (1000 draws), so the 8-day periods
  within a year are not counted as independent. A cell is "clear" when the
  interval excludes zero. Neighbouring cells are correlated, so the share of
  clear cells is descriptive (about 5% would be expected by chance alone).

Associations, not causes; SST is a model input (see the SST-free sensitivity in
climate.md for the regional version).
"""

from __future__ import annotations

import json
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config

GROUPS = ("peuk", "syn", "pro")
NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes",
         "peuk_nosst": "picoeukaryotes, no SST", "syn_nosst": "Synechococcus, no SST",
         "pro_nosst": "Prochlorococcus, no SST"}
CELL = 6  # 1/24° pixels per 0.25° cell
MIN_PIXELS = 9  # of 36
LEAD_MONTHS = 2
EVENTS = {
    "enso": ("El Niño − La Niña", "ONI", 0.5),
    "iod": ("positive − negative IOD", "DMI", 0.4),
    "heatwave": ("marine heatwave − normal", None, None),
}


# ---------------------------------------------------------------- building blocks


def coarsen(values: np.ndarray, ok: np.ndarray, cell: int = CELL, min_pixels: int = MIN_PIXELS) -> np.ndarray:
    """Mean of the `ok` pixels in each cell x cell block; NaN where fewer than `min_pixels` are ok."""
    ny, nx = values.shape[0] // cell, values.shape[1] // cell
    v = np.where(ok, values, 0.0).reshape(ny, cell, nx, cell)
    n = ok.reshape(ny, cell, nx, cell).sum(axis=(1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = v.sum(axis=(1, 3)) / n
    return np.where(n >= min_pixels, out, np.nan).astype(np.float32)


def slot_anomaly(cube: np.ndarray, time: pd.DatetimeIndex) -> np.ndarray:
    """Subtract each cell's mean for the same 8-day slot of the year (time is the first axis)."""
    slot = ((time.day_of_year - 1) // 8).to_numpy()
    out = np.empty_like(cube)
    for s in np.unique(slot):
        k = slot == s
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN cells stay NaN
            out[k] = cube[k] - np.nanmean(cube[k], axis=0)
    return out


def composite(anom: np.ndarray, years: np.ndarray, in_a: np.ndarray, in_b: np.ndarray, n_boot: int = 1000,
              seed: int = 0, min_periods: int = 5, min_years: int = 2) -> dict[str, np.ndarray]:
    """Per-cell mean(a) - mean(b) with a 95% interval from resampling whole years.

    in_a / in_b are boolean masks of the same shape as `anom` or of its time axis only.
    Cells with fewer than `min_periods` valid periods or `min_years` years in either class are NaN.
    """
    shape = anom.shape[1:]
    T = anom.shape[0]
    flat = anom.reshape(T, -1)
    ok = np.isfinite(flat)

    def mask(m):
        m = np.asarray(m, bool)
        return (np.broadcast_to(m[:, None], flat.shape) if m.ndim == 1 else m.reshape(T, -1)) & ok

    a, b = mask(in_a), mask(in_b)
    uy = np.unique(years)
    by_year = np.stack([years == y for y in uy])  # (n_years, T)
    sa, sb = by_year @ np.where(a, flat, 0.0), by_year @ np.where(b, flat, 0.0)
    ca, cb = by_year @ a.astype(float), by_year @ b.astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        diff = sa.sum(0) / ca.sum(0) - sb.sum(0) / cb.sum(0)
        enough = ((ca.sum(0) >= min_periods) & (cb.sum(0) >= min_periods)
                  & ((ca > 0).sum(0) >= min_years) & ((cb > 0).sum(0) >= min_years))
        rng = np.random.default_rng(seed)
        draws = np.empty((n_boot, flat.shape[1]), np.float32)
        for i in range(0, n_boot, 200):  # 200 draws at a time keeps memory small
            pick = rng.integers(0, uy.size, (min(200, n_boot - i), uy.size))
            w = np.stack([np.bincount(p, minlength=uy.size) for p in pick]).astype(float)
            draws[i:i + w.shape[0]] = (w @ sa) / (w @ ca) - (w @ sb) / (w @ cb)
    lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    diff, lo, hi = (np.where(enough, x, np.nan).reshape(shape) for x in (diff, lo, hi))
    return {"diff": diff, "lo": lo, "hi": hi, "n_a": ca.sum(0).reshape(shape), "n_b": cb.sum(0).reshape(shape)}


def index_classes(idx: pd.Series, time: pd.DatetimeIndex, threshold: float, lead: int = LEAD_MONTHS):
    """(positive, negative) masks: the index `lead` months before each period's mid-month."""
    mid = time + pd.Timedelta(days=3.5)
    months = (mid.to_period("M") - lead).to_timestamp()
    v = idx.reindex(months).to_numpy()
    return v >= threshold, v <= -threshold


def heatwave_classes(frac: np.ndarray, time: pd.DatetimeIndex):
    """(heatwave, normal) masks per period and cell; normal = same calendar months as the cell's heatwaves, < 10%."""
    hot = frac >= 0.5
    month = time.month.to_numpy()
    has_hot = np.stack([hot[month == m].any(axis=0) for m in range(1, 13)])  # (12, Y, X)
    normal = (frac < 0.1) & has_hot[month - 1]
    return hot, normal


# ---------------------------------------------------------------- data


def all_8d_starts(first: int, last: int) -> pd.DatetimeIndex:
    """Start dates of NASA's 8-day composites (day of year 1, 9, ..., 361) for each year."""
    return pd.DatetimeIndex([pd.Timestamp(y, 1, 1) + pd.Timedelta(days=d) for y in range(first, last + 1)
                             for d in range(0, 366, 8) if d < 365 + (y % 4 == 0)])


def _periods(cfg: Config, sensor: str) -> list[Path]:
    return sorted((cfg.path("products") / "reconstruction" / sensor / "8D").glob("recon_*.nc"))


def cell_cube(cfg: Config, sensor: str = "snpp") -> xr.Dataset:
    """(time, lat, lon) 0.25° cell values of each group, SD units, reliable pixels only. Cached."""
    out = cfg.path("interim") / "climate" / f"cells_{sensor}.nc"
    if out.exists():
        return xr.open_dataset(out).load()
    meta = json.loads((cfg.path("interim") / "reconstruct" / f"final_{sensor}_operational.json").read_text())
    times, data, lat, lon = [], {g: [] for g in GROUPS}, None, None
    for path in _periods(cfg, sensor):
        with xr.open_dataset(path) as ds:
            if lat is None:
                lat = ds["lat"].values.reshape(-1, CELL).mean(1)
                lon = ds["lon"].values.reshape(-1, CELL).mean(1)
            times.append(pd.Timestamp(ds["time"].values[0]))
            for g in GROUPS:
                v = ds[g].values[0]
                ok = np.isfinite(v) & (ds[f"{g}_out_of_range"].values[0] == 0)
                data[g].append(coarsen((v - meta[g]["index_mean"]) / meta[g]["index_sd"], ok))
    order = np.argsort(lat)  # ascending, like OISST
    lat = np.round(lat[order], 4)
    lon = np.round(lon, 4)
    cube = xr.Dataset({g: (("time", "lat", "lon"), np.stack(data[g])[:, order]) for g in GROUPS},
                      coords={"time": times, "lat": lat, "lon": lon})
    out.parent.mkdir(parents=True, exist_ok=True)
    cube.to_netcdf(out)
    return cube


def sst_free_cube(cfg: Config, sensor: str = "snpp", group: str = "peuk") -> xr.DataArray:
    """0.25° cells of the SST-free record (climate.sst_free_maps), SD units. Cached."""
    from picopast.climate import sst_free_maps

    out = cfg.path("interim") / "climate" / f"cells_{group}_no_sst_{sensor}.nc"
    if out.exists():
        return xr.open_dataarray(out).load()
    times, maps, lat, lon = [], [], None, None
    for start, z, la, lo in sst_free_maps(cfg, sensor, group):
        if lat is None:
            lat, lon = la.reshape(-1, CELL).mean(1), lo.reshape(-1, CELL).mean(1)
        times.append(start)
        maps.append(coarsen(z, np.isfinite(z)))
    order = np.argsort(lat)
    da = xr.DataArray(np.stack(maps)[:, order], dims=("time", "lat", "lon"), name=group,
                      coords={"time": times, "lat": np.round(lat[order], 4), "lon": np.round(lon, 4)})
    out.parent.mkdir(parents=True, exist_ok=True)
    da.to_netcdf(out)
    return da


def heatwave_fraction(cfg: Config, periods: pd.DatetimeIndex) -> xr.DataArray:
    """Share of each OISST cell's days in a marine heatwave, per 8-day period starting at `periods`. Cached."""
    from picopast import mhw

    out = cfg.path("interim") / "climate" / "heatwave_cells_8D.nc"
    if out.exists():
        return xr.open_dataarray(out).load().sel(time=periods)
    wanted, periods = periods, all_8d_starts(2012, 2026)
    paths = sorted((cfg.path("raw") / "oisst").glob("oisst_v2.1_box_*.nc"))
    sst = xr.open_mfdataset(paths, combine="by_coords")["sst"]
    time = pd.DatetimeIndex(sst["time"].values)
    start = time.searchsorted(periods)
    stop = np.minimum(periods + pd.Timedelta(days=8), pd.to_datetime(periods.year.astype(str) + "-12-31")
                      + pd.Timedelta(days=1))  # the last period of a year ends on 31 December
    end = time.searchsorted(stop)
    frac = np.full((periods.size, sst.sizes["lat"], sst.sizes["lon"]), np.nan, np.float32)
    for i in range(0, sst.sizes["lat"], 20):  # strips of 20 rows (16 GB machine)
        arr = sst.isel(lat=slice(i, i + 20)).values.astype(np.float32)
        _, thresh = mhw.climatology(arr, time)
        hot = mhw.heatwave_days(arr, time, thresh).astype(np.int32)
        csum = np.concatenate([np.zeros((1,) + hot.shape[1:], np.int32), np.cumsum(hot, axis=0)])
        days = (end - start)[:, None, None]
        f = (csum[end] - csum[start]) / np.maximum(days, 1)
        ocean = np.isfinite(arr).all(axis=0)
        frac[:, i:i + 20] = np.where(ocean & (days > 0), f, np.nan)  # no OISST yet -> NaN
        del arr, hot, csum
    da = xr.DataArray(frac, dims=("time", "lat", "lon"),
                      coords={"time": periods, "lat": sst["lat"].values, "lon": sst["lon"].values},
                      name="heatwave_fraction")
    out.parent.mkdir(parents=True, exist_ok=True)
    da.to_netcdf(out)
    return da.sel(time=wanted)


def all_composites(cfg: Config, sensor: str = "snpp", years: tuple[int, int] | None = None,
                   sst_free: bool = False) -> dict:
    """{(event, group): composite dict} plus the cube's coordinates, for one sensor (optionally a year range)."""
    from picopast.climate import load_indices

    cube = cell_cube(cfg, sensor)
    if years is not None:
        cube = cube.sel(time=slice(f"{years[0]}-01-01", f"{years[1]}-12-31"))
    time = pd.DatetimeIndex(cube["time"].values)
    hw = heatwave_fraction(cfg, time)
    assert np.allclose(hw["lat"].values, cube["lat"].values) and np.allclose(hw["lon"].values, cube["lon"].values)
    idx = load_indices(cfg)
    yrs = time.year.to_numpy()
    classes = {}
    for event, (_, name, thr) in EVENTS.items():
        classes[event] = (heatwave_classes(hw.values, time) if name is None
                          else index_classes(idx[name].dropna(), time, thr))
    res = {"lat": cube["lat"].values, "lon": cube["lon"].values, "time": time}
    arrays = {g: cube[g].values for g in GROUPS}
    if sst_free:
        nosst = sst_free_cube(cfg, sensor).sel(time=time)
        assert np.allclose(nosst["lat"].values, cube["lat"].values)
        arrays["peuk_nosst"] = nosst.values
    for g, values in arrays.items():
        anom = slot_anomaly(values, time)
        for event, (a, b) in classes.items():
            res[(event, g)] = composite(anom, yrs, a, b)
            if a.ndim == 1:
                res[(event, g)]["years_a"] = sorted(set(yrs[a]))
                res[(event, g)]["years_b"] = sorted(set(yrs[b]))
    res["ever_valid"] = np.isfinite(cube["peuk"].values).any(axis=0)
    return res


# ---------------------------------------------------------------- figure and report

DIVERGING = ["#173f73", "#2a78d6", "#f0efec", "#e34948", "#8c1f1f"]  # blue <-> red, gray midpoint


def plot(res: dict, groups: tuple[str, ...], tiers: dict, out: Path, title: str) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    ink, muted, surface, land = "#2b2b29", "#6b6a63", "#fcfcfb", "#d6d4cc"
    cmap = LinearSegmentedColormap.from_list("div", DIVERGING)
    cmap.set_bad(land)
    vals = np.concatenate([np.abs(res[(e, g)]["diff"]).ravel() for e in EVENTS for g in groups])
    vmax = float(np.ceil(np.nanpercentile(vals, 98) * 10) / 10)
    lat, lon = res["lat"], res["lon"]
    fig, axes = plt.subplots(len(groups), 3, figsize=(15, 3.3 * len(groups) + 0.9), facecolor=surface,
                             squeeze=False, layout="constrained")
    for r, g in enumerate(groups):
        kind = tiers[g.split("_")[0]].split(":")[0].strip()
        for c, (event, (label, _, _)) in enumerate(EVENTS.items()):
            ax = axes[r, c]
            comp = res[(event, g)]
            d = comp["diff"]
            shown = np.where(res["ever_valid"], d, np.nan)
            im = ax.pcolormesh(lon, lat, shown, cmap=cmap, vmin=-vmax, vmax=vmax, shading="nearest")
            clear = np.isfinite(d) & ((comp["lo"] > 0) | (comp["hi"] < 0))
            yy, xx = np.nonzero(clear[::2, ::2])
            ax.scatter(lon[::2][xx], lat[::2][yy], s=1.2, c=ink, lw=0, alpha=0.55)
            ax.set_aspect("equal")
            ax.set_facecolor(land)
            share = clear.sum() / max(np.isfinite(d).sum(), 1)
            name = NAMES[g][0].upper() + NAMES[g][1:]
            tag = "EXPLORATORY" if kind == "exploratory" else kind
            ax.set_title(f"{name} — {label}\n{tag} · clear cells (dots): {share:.0%}", loc="left", fontsize=9,
                         color=ink, fontweight="bold" if kind == "exploratory" else "normal")
            ax.tick_params(colors=muted, labelsize=7)
            for s in ax.spines.values():
                s.set_color(muted)
                s.set_linewidth(0.6)
    cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.04, aspect=50, pad=0.02)
    cb.set_label("difference in anomaly (SD of the index)   ← lower during the event · higher during the event →",
                 fontsize=8, color=ink)
    cb.ax.tick_params(colors=muted, labelsize=7)
    cb.outline.set_edgecolor(muted)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=ink)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=surface)
    plt.close(fig)
    return out


def _region_masks(cfg: Config, lat: np.ndarray, lon: np.ndarray) -> dict[str, np.ndarray]:
    from picopast.climate import REGIONS

    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    out = {}
    for k, label in REGIONS.items():
        if k == "study_box":
            continue
        r = cfg["regions"][k]
        out[label] = (lat2 >= r["lat_min"]) & (lat2 <= r["lat_max"]) & (lon2 >= r["lon_min"]) & (lon2 <= r["lon_max"])
    lat_band = {"north of 15°N": lat2 >= 15, "5–15°N": (lat2 >= 5) & (lat2 < 15), "south of 5°N": lat2 < 5}
    west = lon2 < 78
    for k, m in lat_band.items():
        out[f"Arabian Sea side, {k}"] = m & west
        out[f"Bay of Bengal side, {k}"] = m & ~west
    return out


def _pattern_r(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 10 else np.nan


def report(cfg: Config) -> Path:
    tiers = cfg["reconstruction"]["tiers"]
    figs = cfg.path("reports") / "figures"
    snpp = all_composites(cfg, "snpp", sst_free=True)
    plot(snpp, ("peuk", "peuk_nosst"), tiers, figs / "climate_maps_peuk.png",
         "Picoeukaryote anomaly during climate events, VIIRS SNPP 2012–2026, 0.25° cells (dots: 95% interval "
         "excludes zero). Top: main record. Bottom: the same model without SST.")
    plot(snpp, ("syn", "pro"), tiers, figs / "climate_maps_exploratory.png",
         "EXPLORATORY: cyanobacteria during climate events (did not pass the trust gate), VIIRS SNPP 2012–2026")
    have_n20 = bool(_periods(cfg, "noaa20"))
    if have_n20:
        n20 = all_composites(cfg, "noaa20")
        t0 = pd.DatetimeIndex(n20["time"])
        overlap = (int(t0.year.min()), int(t0.year.max()))
        snpp_same = all_composites(cfg, "snpp", overlap)

    L = ["# Climate maps (RQ5, per 0.25° cell)", "",
         f"Generated {date.today().isoformat()} by `picopast climate-maps`. Where in the North Indian Ocean the "
         "back-extended record moves with ENSO, the IOD and marine heatwaves. Same idea as the regional analysis "
         "(`climate.md`), cell by cell on the OISST grid. Associations, not causes; only picoeukaryotes passed the "
         "trust gate.", "",
         "- Cell value: mean of the reliable 1/24° pixels (all model inputs inside the training range) in each 0.25° "
         f"cell, as a relative index in SD units; at least {MIN_PIXELS} of 36 pixels needed. Anomaly: minus the "
         "cell's mean for the same 8-day slot of the year.",
         f"- ENSO: ONI {LEAD_MONTHS} months earlier ≥ +0.5 vs ≤ −0.5. IOD: DMI {LEAD_MONTHS} months earlier ≥ +0.4 "
         "vs ≤ −0.4 (the strongest lag regionally, for both). ENSO and IOD events overlap (2015, 2019, 2023), so "
         "the two maps are not independent.",
         "- Heatwave: periods with at least half of the cell's days in a marine heatwave (Hobday, OISST, 1991–2020 "
         "baseline) vs periods of the same calendar months under 10%, in that same cell.",
         "- 95% interval from resampling whole years (1000 draws). \"Clear\" = interval excludes zero. Neighbouring "
         "cells are correlated, so the share of clear cells is descriptive; about 5% would be expected by chance.",
         "", "![picoeukaryotes](figures/climate_maps_peuk.png)", "",
         "## How much of the ocean moves, and which way", "",
         "| group | tier | event | cells with a value | mean difference (SD) | lower / higher (share of cells) | "
         "clear lower / clear higher |", "| --- | --- | --- | --- | --- | --- | --- |"]
    tier = {g: "relative index" if tiers[g].startswith("relative") else "EXPLORATORY" for g in GROUPS}
    for g in GROUPS:
        for event, (label, _, _) in EVENTS.items():
            c = snpp[(event, g)]
            d = c["diff"]
            ok = np.isfinite(d)
            n = ok.sum()
            lower, higher = (d[ok] < 0).mean(), (d[ok] > 0).mean()
            cl, ch = (c["hi"][ok] < 0).mean(), (c["lo"][ok] > 0).mean()
            L.append(f"| {NAMES[g]} | {tier[g]} | {label} | {n:,} | {np.nanmean(d):+.2f} | {lower:.0%} / "
                     f"{higher:.0%} | {cl:.0%} / {ch:.0%} |")
    e = snpp[("enso", "peuk")]
    i = snpp[("iod", "peuk")]
    L += ["", f"Years in each class (any period): El Niño {', '.join(map(str, e['years_a']))}; La Niña "
          f"{', '.join(map(str, e['years_b']))}; positive IOD {', '.join(map(str, i['years_a']))}; negative IOD "
          f"{', '.join(map(str, i['years_b']))}.", "",
          "## Picoeukaryotes by area", "",
          "Mean difference over the cells of each area (SD of the index), with the share of its cells that are clear "
          "in the same direction as the mean. Arabian Sea / Bay of Bengal sides split at 78°E.", "",
          "| area | El Niño − La Niña | positive − negative IOD | heatwave − normal |", "| --- | --- | --- | --- |"]
    for area, m in _region_masks(cfg, snpp["lat"], snpp["lon"]).items():
        cells = []
        for event in EVENTS:
            c = snpp[(event, "peuk")]
            d = np.where(m, c["diff"], np.nan)
            ok = np.isfinite(d)
            if ok.sum() < 10:
                cells.append("—")
                continue
            mean = np.nanmean(d)
            same = (c["hi"][ok] < 0).mean() if mean < 0 else (c["lo"][ok] > 0).mean()
            cells.append(f"{mean:+.2f} ({same:.0%} clear)")
        L.append(f"| {area} | " + " | ".join(cells) + " |")
    if have_n20:
        L += ["", f"## Sensor check: NOAA-20 vs SNPP, {overlap[0]}–{overlap[1]}", "",
              "The same composites from the NOAA-20 record and from SNPP restricted to the same years. Pattern r = "
              "correlation between the two maps over cells where both have a value: high means the two satellites "
              "put the change in the same places.", "",
              "| group | tier | event | pattern r | mean difference, SNPP | mean difference, NOAA-20 |",
              "| --- | --- | --- | --- | --- | --- |"]
        for g in GROUPS:
            for event, (label, _, _) in EVENTS.items():
                a, b = snpp_same[(event, g)]["diff"], n20[(event, g)]["diff"]
                L.append(f"| {NAMES[g]} | {tier[g]} | {label} | {_pattern_r(a, b):.2f} | {np.nanmean(a):+.2f} | "
                         f"{np.nanmean(b):+.2f} |")
    L += ["", "## Is the picoeukaryote signal put there by SST? (model without SST)", "",
          "SST is one of the model's inputs, and heatwaves (and much of ENSO's and the IOD's footprint here) are "
          "SST signals, so the model could produce a drop simply by reacting to warm water. The second row of the "
          "figure repeats every composite with the same model trained without SST (VIIRS bands, ratios, PAR, "
          "position and season only), on the same pixels. What survives is carried by the ocean's colour.", "",
          "| event | mean difference, main | mean difference, without SST | survives | clear lower, main | "
          "clear lower, without SST | pattern r |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for event, (label, _, _) in EVENTS.items():
        a, b = snpp[(event, "peuk")], snpp[(event, "peuk_nosst")]
        ma, mb = np.nanmean(a["diff"]), np.nanmean(b["diff"])
        ca = (a["hi"][np.isfinite(a["diff"])] < 0).mean()
        cb = (b["hi"][np.isfinite(b["diff"])] < 0).mean()
        L.append(f"| {label} | {ma:+.2f} | {mb:+.2f} | {mb / ma:.0%} | {ca:.0%} | {cb:.0%} | "
                 f"{_pattern_r(a['diff'], b['diff']):.2f} |")
    L += ["", "\"Survives\" = mean difference without SST / with SST, over the whole map. Pattern r compares where "
          "the two put the change.", "",
          "## Cyanobacteria (exploratory)", "",
          "These are the MOANA-based records learned from VIIRS (Route A), shown for completeness: they did not "
          "pass the trust gate and the two satellites disagree on their long-term change "
          "(`sensor_comparison.md`). The cyanobacteria records that did pass are from Route D, with their own "
          "maps in `route_d_reconstruction.md`.", "",
          "![exploratory](figures/climate_maps_exploratory.png)", ""]
    out = cfg.path("reports") / "climate_maps.md"
    out.write_text("\n".join(L))
    return out
