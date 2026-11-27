"""When cyclones hit heatwaves: the storm table.

Rules: ``docs/cyclone_heatwave_preregistration.md`` (binding from 2026-11-22, committed
before any response below was computed).

One row per 6-hourly IBTrACS fix inside the box, over the sea, at storm strength. For each
fix: heatwave exposure on D0 - 3 and the ocean state before the storm, then the responses:
cooling, heatwave survival, re-warming, the chlorophyll bloom and, in the PACE era, the
picophytoplankton groups. Everything is a disc mean around the fix (150 km; 75 and 300 km
as sensitivities).
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config

R_KM = 150.0
RADII = (75.0, 150.0, 300.0)
LAG = 3
LAGS = (3, 5)
HOT, NORMAL = 0.5, 0.1
MIN_WIND_KT = 34.0
MIN_LAND_KM = 50.0
SPLIT_LON = 77.5
CONTROL_KM = 750.0
MIN_VALID_SHARE = 0.2
REWARM_TOL = 0.2
REWARM_MAX_DAYS = 60
EARTH_KM = 6371.0
CUBE_START = pd.Timestamp("2011-11-01")
_CHL_PERIOD = re.compile(r"\.(\d{8})_(\d{8})\.L3m\.8D\.")


# ---------------------------------------------------------------- geometry and calendar

def haversine_km(lat1, lon1, lat2, lon2) -> np.ndarray:
    """Great-circle distance in km (inputs in degrees, broadcast)."""
    p = np.pi / 180
    la1, lo1, la2, lo2 = (np.asarray(x, float) * p for x in (lat1, lon1, lat2, lon2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * EARTH_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def season_south(month) -> np.ndarray:
    """South Indian Ocean cyclone seasons: early (November-January), late (February-April), other."""
    m = np.asarray(month)
    return np.where((m >= 11) | (m == 1), "SIO early", np.where((m >= 2) & (m <= 4), "SIO late", "SIO other"))


def is_south(cfg: Config) -> bool:
    return cfg.box.lat_max <= 0


def season(month) -> np.ndarray:
    """Pre-monsoon (March-June), post-monsoon (September-December), other."""
    m = np.asarray(month)
    return np.where((m >= 3) & (m <= 6), "pre-monsoon", np.where(m >= 9, "post-monsoon", "other"))


def disc(glat, glon, lat0: float, lon0: float, r_km: float, valid=None) -> tuple[np.ndarray, np.ndarray]:
    """(iy, ix) of the cells of a 1-D lat/lon grid whose centres lie within r_km of (lat0, lon0)."""
    glat, glon = np.asarray(glat, float), np.asarray(glon, float)
    step = max(abs(glat[1] - glat[0]) if glat.size > 1 else 0.0, abs(glon[1] - glon[0]) if glon.size > 1 else 0.0)
    dlat = r_km / 111.0 + step
    dlon = r_km / (111.0 * max(np.cos(np.radians(lat0)), 0.2)) + step
    iy = np.flatnonzero(np.abs(glat - lat0) <= dlat)
    ix = np.flatnonzero(np.abs(glon - lon0) <= dlon)
    if not iy.size or not ix.size:
        return iy[:0], ix[:0]
    la, lo = np.meshgrid(glat[iy], glon[ix], indexing="ij")
    keep = haversine_km(lat0, lon0, la, lo) <= r_km
    if valid is not None:
        keep &= np.asarray(valid, bool)[np.ix_(iy, ix)]
    yy, xx = np.nonzero(keep)
    return iy[yy], ix[xx]


def pre_post(starts: pd.DatetimeIndex, ends: pd.DatetimeIndex, d0: pd.Timestamp) -> tuple[int, int, int]:
    """Indices of the pre period (latest ending on or before D0 - 2) and the first two post periods
    (starting on or after D0 + 1); -1 where none exists."""
    pre = np.flatnonzero(ends <= d0 - pd.Timedelta(days=2))
    post = np.flatnonzero(starts >= d0 + pd.Timedelta(days=1))
    return (int(pre[-1]) if pre.size else -1, int(post[0]) if post.size else -1,
            int(post[1]) if post.size > 1 else -1)


# ---------------------------------------------------------------- tracks

def all_fixes(cfg: Config, lists: tuple[str, ...] | None = None) -> pd.DataFrame:
    """Every 6-hourly IBTrACS fix, any strength, with translation speed, from the config's basin list
    (``storms.ibtracs``: NI for the north box, SI for the south box) or from the given lists."""
    lists = lists or (cfg.get("storms", {}).get("ibtracs", "NI"),)
    cols = ["SID", "NAME", "ISO_TIME", "LAT", "LON", "WMO_WIND", "USA_WIND", "DIST2LAND"]
    t = pd.concat([pd.read_csv(cfg.path("raw") / "ibtracs" / f"ibtracs.{b}.list.csv", skiprows=[1], usecols=cols,
                               low_memory=False, na_values=[" ", ""]) for b in lists], ignore_index=True)
    t["time"] = pd.to_datetime(t["ISO_TIME"], errors="coerce")
    for c in ("LAT", "LON", "WMO_WIND", "USA_WIND", "DIST2LAND"):
        t[c] = pd.to_numeric(t[c], errors="coerce")
    t = t[t["time"].notna() & t["LAT"].notna() & t["LON"].notna() & (t["time"].dt.minute == 0)
          & (t["time"].dt.hour % 6 == 0)]
    t = t.sort_values(["SID", "time"]).drop_duplicates(["SID", "time"]).reset_index(drop=True)
    g = t.groupby("SID")
    a_lat, a_lon, a_t = (g[c].shift(1).fillna(t[c]) for c in ("LAT", "LON", "time"))
    b_lat, b_lon, b_t = (g[c].shift(-1).fillna(t[c]) for c in ("LAT", "LON", "time"))
    hours = (pd.to_datetime(b_t) - pd.to_datetime(a_t)).dt.total_seconds() / 3600
    dist = haversine_km(a_lat, a_lon, b_lat, b_lon)
    t["speed_kmh"] = np.where(hours > 0, dist / hours.where(hours > 0), np.nan)
    return t.rename(columns={"SID": "sid", "NAME": "name", "LAT": "lat", "LON": "lon", "USA_WIND": "wind",
                             "WMO_WIND": "wind_imd", "DIST2LAND": "dist2land_km"})


def track_points(cfg: Config, start: str = "2012-01-01", end: str | None = None) -> pd.DataFrame:
    """The analysis points: from `start` (to `end`), inside the box, >= 50 km from land, >= 34 kt by JTWC
    (primary) or IMD. The main study starts in 2012; the replication covers 4 July 2002 to 31 December 2011."""
    t = all_fixes(cfg)
    b = cfg.box
    in_time = (t["time"] >= pd.Timestamp(start)) & (t["time"].dt.normalize() <= pd.Timestamp(end or "2100-01-01"))
    keep = (in_time & t["lat"].between(b.lat_min, b.lat_max)
            & t["lon"].between(b.lon_min, b.lon_max) & (t["dist2land_km"] >= MIN_LAND_KM)
            & ((t["wind"] >= MIN_WIND_KT) | (t["wind_imd"] >= MIN_WIND_KT)))
    p = t[keep].copy()
    p["primary"] = p["wind"] >= MIN_WIND_KT
    p["imd"] = p["wind_imd"] >= MIN_WIND_KT
    p["d0"] = p["time"].dt.normalize()
    if is_south(cfg):
        p["basin"], p["season"] = "South Indian Ocean", season_south(p["time"].dt.month)
    else:
        p["basin"] = np.where(p["lon"] < SPLIT_LON, "Arabian Sea", "Bay of Bengal")
        p["season"] = season(p["time"].dt.month)
    doy = p["time"].dt.dayofyear.to_numpy()
    p["doy_sin"], p["doy_cos"] = np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25)
    p = p.reset_index(drop=True)
    p["point"] = np.arange(len(p))
    return p[["point", "sid", "name", "time", "d0", "lat", "lon", "wind", "wind_imd", "primary", "imd",
              "speed_kmh", "dist2land_km", "basin", "season", "doy_sin", "doy_cos"]]


# ---------------------------------------------------------------- daily OISST heatwave cube

def daily_ocean(cfg: Config, start: pd.Timestamp = CUBE_START, end: pd.Timestamp | None = None) -> xr.Dataset:
    """Daily OISST from `start` (default November 2011) to `end`: sst, ssta (vs 1991-2020), heatwave excess and
    heatwave flag. Heatwaves always come from the full record and the fixed 1991-2020 baseline. Cached."""
    from picopast import mhw

    name = ("oisst_daily_mhw.nc" if start == CUBE_START and end is None
            else f"oisst_daily_mhw_{pd.Timestamp(start):%Y%m%d}_{pd.Timestamp(end):%Y%m%d}.nc")
    out = cfg.path("interim") / "storms" / name
    if out.exists():
        return xr.open_dataset(out)
    paths = sorted((cfg.path("raw") / "oisst").glob("oisst_v2.1_box_*.nc"))
    src = xr.open_mfdataset(paths, combine="by_coords")["sst"]
    time = pd.DatetimeIndex(src["time"].values)
    keep = np.asarray((time >= pd.Timestamp(start)) & (time <= pd.Timestamp(end or "2100-01-01")))
    tk = time[keep]
    doy = mhw.day_of_year_366(tk)
    ny, nx = src.sizes["lat"], src.sizes["lon"]
    shape = (tk.size, ny, nx)
    sst, ssta, exc = (np.full(shape, np.nan, np.float32) for _ in range(3))
    hot = np.zeros(shape, np.uint8)
    for i in range(0, ny, 20):  # strips of 20 rows (16 GB machine)
        arr = src.isel(lat=slice(i, i + 20)).values.astype(np.float32)
        seas, thresh = mhw.climatology(arr, time)
        h = mhw.heatwave_days(arr, time, thresh)[keep]
        a = arr[keep]
        sst[:, i:i + 20] = a
        ssta[:, i:i + 20] = a - seas[doy - 1]
        exc[:, i:i + 20] = np.where(h, a - thresh[doy - 1], 0.0)
        hot[:, i:i + 20] = h
        del arr, h, a
    ds = xr.Dataset({"sst": (("time", "lat", "lon"), sst), "ssta": (("time", "lat", "lon"), ssta),
                     "excess": (("time", "lat", "lon"), exc), "hot": (("time", "lat", "lon"), hot)},
                    coords={"time": tk, "lat": src["lat"].values, "lon": src["lon"].values})
    ds.attrs = {"source": "NOAA OISST v2.1; heatwaves by picopast.mhw (Hobday et al. 2016, 1991-2020 baseline)"}
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    for v in ("sst", "ssta", "excess"):
        enc[v].update(dtype="int16", scale_factor=0.01, _FillValue=-32768)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out, encoding=enc)
    return xr.open_dataset(out)


def _control_share(hot: np.ndarray, k_from: int, k_to: int, basin_mask: np.ndarray, far: np.ndarray) -> float:
    """Share of heatwave cells on day k_from (inside basin_mask and far from storms) still hot on day k_to."""
    was = (hot[k_from] == 1) & basin_mask & far
    return float((hot[k_to][was] == 1).mean()) if was.any() else np.nan


def ocean_state(cfg: Config, pts: pd.DataFrame, fixes: pd.DataFrame, ocean: xr.Dataset | None = None) -> pd.DataFrame:
    """Exposure, pre-storm state, cooling, heatwave survival and re-warming for every point."""
    ocean = daily_ocean(cfg) if ocean is None else ocean
    time = pd.DatetimeIndex(ocean["time"].values)
    day = {d: i for i, d in enumerate(time)}
    lat, lon = ocean["lat"].values, ocean["lon"].values
    sst = ocean["sst"].values.astype(np.float32)
    ssta = ocean["ssta"].values.astype(np.float32)
    exc = ocean["excess"].values.astype(np.float32)
    hot = ocean["hot"].values
    is_ocean = np.isfinite(sst).mean(axis=0) > 0.5
    la2, lo2 = np.meshgrid(lat, lon, indexing="ij")
    basins = {"Arabian Sea": is_ocean & (lo2 < SPLIT_LON), "Bay of Bengal": is_ocean & (lo2 >= SPLIT_LON),
              "South Indian Ocean": is_ocean}
    control_cache: dict = {}

    def disc_mean(arr, k, cells):
        v = arr[k][cells]
        return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan

    rows = []
    for p in pts.itertuples():
        row = {"point": p.point}
        d0 = p.d0
        k0 = day.get(d0)
        for r in RADII:
            cells = disc(lat, lon, p.lat, p.lon, r, is_ocean)
            tag = "" if r == R_KM else f"_r{int(r)}"
            for lag in LAGS:
                ltag = tag + ("" if lag == LAG else f"_l{lag}")
                k = day.get(d0 - pd.Timedelta(days=lag))
                if k is None or not cells[0].size:
                    row[f"E{ltag}"] = np.nan
                    continue
                h = hot[k][cells]
                row[f"E{ltag}"] = float(h.mean())
                if r == R_KM and lag == LAG:
                    row["sst_pre"] = disc_mean(sst, k, cells)
                    row["ssta_pre"] = disc_mean(ssta, k, cells)
                    row["excess_pre"] = disc_mean(exc, k, cells)
                    row["n_cells"] = int(cells[0].size)
            # cooling: min over D0..D0+5 minus mean over D0-7..D0-3
            ks_pre = [day.get(d0 - pd.Timedelta(days=j)) for j in range(7, 2, -1)]
            ks_post = [day.get(d0 + pd.Timedelta(days=j)) for j in range(0, 6)]
            if k0 is None or any(k is None for k in ks_pre + ks_post) or not cells[0].size:
                row[f"cooling{tag}"] = np.nan
            else:
                pre = np.nanmean([disc_mean(sst, k, cells) for k in ks_pre])
                post = [disc_mean(sst, k, cells) for k in ks_post]
                row[f"cooling{tag}"] = float(np.nanmin(post) - pre) if np.isfinite(post).any() else np.nan
            if r != R_KM:
                continue
            # heatwave survival, D0-3 -> D0+5, storm disc vs far-away control cells in the same basin
            ka, kb = day.get(d0 - pd.Timedelta(days=LAG)), day.get(d0 + pd.Timedelta(days=5))
            row["survival_storm"] = row["survival_control"] = np.nan
            if ka is not None and kb is not None and cells[0].size:
                was = hot[ka][cells] == 1
                if was.any():
                    row["survival_storm"] = float((hot[kb][cells][was] == 1).mean())
                key = (d0, p.basin)
                if key not in control_cache:
                    near = fixes[(fixes["time"] >= d0 - pd.Timedelta(days=10)) & (fixes["time"] <= d0 + pd.Timedelta(days=10))]
                    far = np.ones_like(is_ocean)
                    for f in near.itertuples():
                        far &= haversine_km(f.lat, f.lon, la2, lo2) > CONTROL_KM
                    control_cache[key] = _control_share(hot, ka, kb, basins[p.basin], far)
                row["survival_control"] = control_cache[key]
            # re-warming: days from D0 until the disc SSTA is back within 0.2 C of its D0-3 value
            row["rewarm_days"], row["rewarm_event"] = np.nan, False
            if ka is not None and k0 is not None and cells[0].size:
                base = disc_mean(ssta, ka, cells)
                series = [disc_mean(ssta, k0 + j, cells) if k0 + j < time.size else np.nan
                          for j in range(REWARM_MAX_DAYS + 1)]
                s = np.asarray(series)
                if np.isfinite(base) and np.isfinite(s[:6]).any():
                    j_min = int(np.nanargmin(s[:6]))
                    back = np.flatnonzero(s[j_min:] >= base - REWARM_TOL)
                    observed = np.isfinite(s).sum()
                    if back.size:
                        row["rewarm_days"], row["rewarm_event"] = float(j_min + back[0]), True
                    else:
                        row["rewarm_days"] = float(min(REWARM_MAX_DAYS, observed - 1))
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- chlorophyll bloom and PACE groups

CHL_PREFIX = {"viirs_snpp_chl": "SNPP_VIIRS", "modis_aqua_chl": "AQUA_MODIS"}


def chl_periods(cfg: Config, product: str = "viirs_snpp_chl") -> pd.DataFrame:
    rows = []
    pattern = f"{CHL_PREFIX[product]}.*.L3m.8D.CHL.chlor_a.4km.box.nc"
    for f in sorted((cfg.path("raw") / product).glob(pattern)):
        m = _CHL_PERIOD.search(f.name)
        rows.append((pd.Timestamp(m.group(1)), pd.Timestamp(m.group(2)), f))
    return pd.DataFrame(rows, columns=["start", "end", "path"]).sort_values("start", ignore_index=True)


def slot(start: pd.Timestamp) -> int:
    return (pd.Timestamp(start).dayofyear - 1) // 8


def chl_climatology(cfg: Config, periods: pd.DataFrame, offshore: np.ndarray,
                    cache_name: str = "chl_slot_climatology.nc") -> np.ndarray:
    """Per-pixel mean log10 chlorophyll per 8-day slot (46, lat, lon), all years on disk. Cached."""
    out = cfg.path("interim") / "storms" / cache_name
    if out.exists():
        return xr.open_dataarray(out).values
    total = np.zeros((46,) + offshore.shape, np.float64)
    count = np.zeros((46,) + offshore.shape, np.int32)
    for p in periods.itertuples():
        with xr.open_dataset(p.path) as ds:
            x = ds["chlor_a"].values.astype(np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where((x > 0) & offshore, np.log10(x), np.nan)
        s = slot(p.start)
        ok = np.isfinite(v)
        total[s][ok] += v[ok]
        count[s] += ok
    with np.errstate(invalid="ignore"):
        clim = np.where(count > 0, total / np.maximum(count, 1), np.nan).astype(np.float32)
    out.parent.mkdir(parents=True, exist_ok=True)
    xr.DataArray(clim, dims=("slot", "lat", "lon"), name="log10_chl_slot_mean").to_netcdf(out)
    return clim


def _disc_flat(lat, lon, pts: pd.DataFrame, offshore: np.ndarray, radii) -> dict:
    """{(point, r): flat indices of offshore pixels in the disc} on the 4 km grid."""
    nx = lon.size
    out = {}
    for p in pts.itertuples():
        for r in radii:
            iy, ix = disc(lat, lon, p.lat, p.lon, r, offshore)
            out[(p.point, r)] = (iy * nx + ix).astype(np.int64)
    return out


def _period_values(field_by_period, needs, idx, min_share=MIN_VALID_SHARE) -> dict:
    """{(point, r, role): disc mean} from one field per period; a mean needs >= min_share valid pixels."""
    out = {}
    for k, field in field_by_period:
        flat = field.ravel()
        for point, r, role in needs.get(k, ()):
            ii = idx[(point, r)]
            if not ii.size:
                continue
            v = flat[ii]
            ok = np.isfinite(v)
            out[(point, r, role)] = float(v[ok].mean()) if ok.mean() >= min_share else np.nan
    return out


def analysis_grid(cfg: Config, product: str = "viirs_snpp_chl") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(lat, lon, offshore) of the 4 km chlorophyll grid.

    North box: PicoPast's offshore mask (unchanged from the main study). South box (final test): water is any
    pixel with valid chlorophyll in at least one 8-day period of that sensor's files, then the 25 km coast rule.
    """
    if not is_south(cfg):
        from picopast.route_d import box_fields

        bf = box_fields(cfg)
        return bf["ds"]["lat"].values, bf["ds"]["lon"].values, bf["offshore"]
    from picopast.grid import coast_mask

    cache = cfg.path("interim") / "storms" / f"offshore_{product}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["lat"], z["lon"], z["offshore"]
    periods = chl_periods(cfg, product)
    water = None
    for f in periods["path"]:
        with xr.open_dataset(f) as ds:
            ok = np.isfinite(ds["chlor_a"].values)
            if water is None:
                water, lat, lon = ok, ds["lat"].values, ds["lon"].values
            else:
                water |= ok
    offshore = coast_mask(water, lat, cfg["grid"]["resolution_deg"], cfg["grid"]["coast_mask_km"])
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, lat=lat, lon=lon, offshore=offshore)
    return lat, lon, offshore


def bloom(cfg: Config, pts: pd.DataFrame, product: str = "viirs_snpp_chl") -> pd.DataFrame:
    """Pre and post chlorophyll anomalies (log10) per point and radius, and the bloom B."""
    lat, lon, offshore = analysis_grid(cfg, product)
    periods = chl_periods(cfg, product)
    clim = chl_climatology(cfg, periods, offshore, "chl_slot_climatology.nc" if product == "viirs_snpp_chl"
                           else f"chl_slot_climatology_{product}.nc")
    starts, ends = pd.DatetimeIndex(periods["start"]), pd.DatetimeIndex(periods["end"])
    idx = _disc_flat(lat, lon, pts, offshore, RADII)
    needs: dict = {}
    chosen = {}
    for p in pts.itertuples():
        k_pre, k1, k2 = pre_post(starts, ends, p.d0)
        chosen[p.point] = (k_pre, k1, k2)
        for role, k in (("pre", k_pre), ("post1", k1), ("post2", k2)):
            if k >= 0:
                for r in RADII:
                    needs.setdefault(k, []).append((p.point, r, role))

    def fields():
        for k in sorted(needs):
            row = periods.iloc[k]
            with xr.open_dataset(row.path) as ds:
                x = ds["chlor_a"].values.astype(np.float32)
            with np.errstate(divide="ignore", invalid="ignore"):
                yield k, np.where((x > 0) & offshore, np.log10(x) - clim[slot(row.start)], np.nan)

    vals = _period_values(fields(), needs, idx)
    rows = []
    for p in pts.itertuples():
        k_pre, k1, k2 = chosen[p.point]
        row = {"point": p.point, "chl_pre_start": starts[k_pre] if k_pre >= 0 else pd.NaT,
               "chl_post_start": starts[k1] if k1 >= 0 else pd.NaT}
        for r in RADII:
            tag = "" if r == R_KM else f"_r{int(r)}"
            pre = vals.get((p.point, r, "pre"), np.nan)
            post = [vals.get((p.point, r, role), np.nan) for role in ("post1", "post2")]
            row[f"chl_pre{tag}"] = pre
            row[f"chl_post1{tag}"], row[f"chl_post2{tag}"] = post
            good = [v for v in post if np.isfinite(v)]
            row[f"bloom{tag}"] = float(np.mean(good) - pre) if good and np.isfinite(pre) else np.nan
            if r == R_KM:
                row["bloom_max"] = float(max(good) - pre) if good and np.isfinite(pre) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def pace_groups(cfg: Config, pts: pd.DataFrame) -> pd.DataFrame:
    """Change in MOANA log10 cells (post mean minus pre) per group for points in the PACE era (150 km)."""
    from picopast.route_d import box_fields

    files = sorted((cfg.path("interim") / "moana" / "operational" / "8D").glob("moana_*.nc"))
    if not files:
        return pd.DataFrame({"point": pts["point"]})
    meta = []
    for f in files:
        with xr.open_dataset(f) as ds:
            meta.append((pd.Timestamp(ds.attrs["period_start"]), pd.Timestamp(ds.attrs["period_end"]), f))
    periods = pd.DataFrame(meta, columns=["start", "end", "path"]).sort_values("start", ignore_index=True)
    starts, ends = pd.DatetimeIndex(periods["start"]), pd.DatetimeIndex(periods["end"])
    bf = box_fields(cfg)
    lat, lon, offshore = bf["ds"]["lat"].values, bf["ds"]["lon"].values, bf["offshore"]
    in_era = pts[(pts["d0"] >= starts[0] + pd.Timedelta(days=10)) & (pts["d0"] <= ends[-1] - pd.Timedelta(days=1))]
    idx = _disc_flat(lat, lon, in_era, offshore, (R_KM,))
    groups = ("peuk", "syn", "pro")
    needs: dict = {}
    for p in in_era.itertuples():
        for role, k in zip(("pre", "post1", "post2"), pre_post(starts, ends, p.d0)):
            if k >= 0:
                needs.setdefault(k, []).append((p.point, R_KM, role))
    out = {g: {} for g in groups}
    for g in groups:
        def fields(g=g):
            for k in sorted(needs):
                with xr.open_dataset(periods.iloc[k].path) as ds:
                    x = ds[g].isel(time=0).values.astype(np.float32) if "time" in ds[g].dims else ds[g].values
                    flags = ds["flags"].isel(time=0).values if "time" in ds["flags"].dims else ds["flags"].values
                with np.errstate(divide="ignore", invalid="ignore"):
                    yield k, np.where((x > 0) & (flags == 0) & offshore, np.log10(x), np.nan)
        out[g] = _period_values(fields(), needs, idx)
    rows = []
    for p in pts.itertuples():
        row = {"point": p.point}
        for g in groups:
            pre = out[g].get((p.point, R_KM, "pre"), np.nan)
            post = [out[g].get((p.point, R_KM, role), np.nan) for role in ("post1", "post2")]
            good = [v for v in post if np.isfinite(v)]
            row[f"d_{g}"] = float(np.mean(good) - pre) if good and np.isfinite(pre) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- RQ3 extended: PicoPast's own group records

#: record name -> (kind, variable). "A" = Route A (MOANA learned from VIIRS), "D" = Route D, "A0" = Route A
#: retrained without SST (computed here). The *_nosst records are primary (pre-registration, RQ3 extended).
GROUP_RECORDS = {"peuk_nosst": ("A0", "peuk"), "syn_nosst": ("D", "syn_nosst"), "pro_nosst": ("D", "pro_nosst"),
                 "peuk": ("A", "peuk"), "syn": ("D", "syn"), "pro": ("D", "pro")}


def peuk_nosst_maps(cfg: Config):
    """Yield (period start, log10 picoeukaryote map) from Route A retrained without SST, all valid pixels.

    Same training table, features minus SST and LightGBM settings as ``climate.sst_free_maps``; maps are cached
    as float16 per period, so a second call only reads them.
    """
    from picopast import reconstruct as R

    cache = cfg.path("interim") / "storms" / "peuk_nosst_8D"
    periods = R.viirs_periods(cfg, "snpp")
    files = [cache / f"{p:%Y%m%d}.npy" for p in periods["start"]]
    if cache.exists() and all(f.exists() for f in files):
        for p, f in zip(periods["start"], files):
            yield p, np.load(f).astype(np.float32)
        return
    cache.mkdir(parents=True, exist_ok=True)
    table = pd.read_parquet(cfg.path("interim") / "reconstruct" / "training_snpp_operational.parquet")
    feats = [f for f in R.FEATURES["A"] if f != "sst"]
    d = table[np.isfinite(table["peuk"])]
    model = R._lgbm().fit(d[feats], d["peuk"].values)
    lat, lon, offshore = R._grid_and_offshore(cfg)
    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    oisst = R.Oisst(cfg, lat, lon)
    for p, f in zip(periods.itertuples(), files):
        v = R.load_viirs(cfg, "snpp", p.start)
        ok = R._inputs_valid(v, oisst.mean(p.start, p.end)) & offshore  # the same pixels as the main record
        idx = np.flatnonzero(ok.ravel())
        out = np.full(ok.size, np.nan, np.float32)
        if idx.size:
            X = R.features({k: a.ravel()[idx] for k, a in v.items()}, np.zeros(idx.size, np.float32),
                           lat2.ravel()[idx], lon2.ravel()[idx], p.start + (p.end - p.start) / 2)
            out[idx] = model.predict(X[feats])
        out = out.reshape(ok.shape)
        np.save(f, out.astype(np.float16))
        yield p.start, out


def _stored_record(cfg: Config, kind: str, var: str, offshore: np.ndarray):
    """Loader start -> log10 map (NaN off the offshore mask) for a stored Route A or Route D record."""
    base = cfg.path("products") / "reconstruction"

    def load(start: pd.Timestamp):
        f = (base / "snpp" / "8D" / f"recon_{start:%Y%m%d}.nc" if kind == "A"
             else base / "snpp_route_d" / "8D" / f"recon_d_{start:%Y%m%d}.nc")
        if not f.exists():
            return None
        with xr.open_dataset(f) as ds:
            x = ds[var].isel(time=0).values.astype(np.float32) if "time" in ds[var].dims else ds[var].values
        return np.where(offshore, x, np.nan)

    return load


def group_changes(cfg: Config, pts: pd.DataFrame) -> pd.DataFrame:
    """Per point and record: pre-period anomaly and the change after the storm (log10), 150 km disc.

    Anomalies are taken against each record's own per-pixel 8-day slot climatology (2012-2026 as on disk).
    """
    from picopast.route_d import box_fields

    bf = box_fields(cfg)
    lat, lon, offshore = bf["ds"]["lat"].values, bf["ds"]["lon"].values, bf["offshore"]
    periods = chl_periods(cfg)
    starts, ends = pd.DatetimeIndex(periods["start"]), pd.DatetimeIndex(periods["end"])
    idx = _disc_flat(lat, lon, pts, offshore, (R_KM,))
    needs, chosen = {}, {}
    for p in pts.itertuples():
        chosen[p.point] = pre_post(starts, ends, p.d0)
        for role, k in zip(("pre", "post1", "post2"), chosen[p.point]):
            if k >= 0:
                needs.setdefault(k, []).append((p.point, R_KM, role))
    out = pd.DataFrame({"point": pts["point"].to_numpy()})
    cache = cfg.path("interim") / "storms" / "peuk_nosst_8D"
    for name, (kind, var) in GROUP_RECORDS.items():
        total = np.zeros((46,) + offshore.shape, np.float32)
        count = np.zeros((46,) + offshore.shape, np.int16)
        if kind == "A0":
            source = ((s, m) for s, m in peuk_nosst_maps(cfg))
        else:
            load = _stored_record(cfg, kind, var, offshore)
            source = ((s, load(s)) for s in starts)
        kept = set(needs)
        for s, m in source:  # pass 1: climatology
            if m is None:
                continue
            ok = np.isfinite(m)
            total[slot(s)][ok] += m[ok]
            count[slot(s)] += ok
        with np.errstate(invalid="ignore", divide="ignore"):
            clim = np.where(count > 0, total / np.maximum(count, 1), np.nan).astype(np.float32)
        del total, count

        def fields(kind=kind, var=var, clim=clim, kept=kept):  # pass 2: anomalies of the needed periods only
            load = None if kind == "A0" else _stored_record(cfg, kind, var, offshore)
            for k in sorted(kept):
                s = starts[k]
                if kind == "A0":
                    f = cache / f"{s:%Y%m%d}.npy"
                    m = np.load(f).astype(np.float32) if f.exists() else None
                else:
                    m = load(s)
                if m is not None:
                    yield k, m - clim[slot(s)]

        vals = _period_values(fields(), needs, idx)
        pre, d = [], []
        for p in pts.itertuples():
            a = vals.get((p.point, R_KM, "pre"), np.nan)
            post = [vals.get((p.point, R_KM, r), np.nan) for r in ("post1", "post2")]
            good = [v for v in post if np.isfinite(v)]
            pre.append(a)
            d.append(float(np.mean(good) - a) if good and np.isfinite(a) else np.nan)
        out[f"recpre_{name}"] = pre
        out[f"rec_{name}"] = d
    return out


def add_group_changes(cfg: Config) -> pd.DataFrame:
    """Add the RQ3-extended columns to the cached storm table (recomputing only them)."""
    path = cfg.path("interim") / "storms" / "storm_table.parquet"
    table = pd.read_parquet(path)
    table = table[[c for c in table.columns if not c.startswith(("rec_", "recpre_"))]]
    table = table.merge(group_changes(cfg, table[["point", "lat", "lon", "d0"]]), on="point", how="left")
    table.to_parquet(path, index=False)
    return table


# ---------------------------------------------------------------- covariates

def _ssh_files(cfg: Config) -> pd.Series:
    files = sorted((cfg.path("raw") / "nasa_ssh").glob("NASA-SSH_alt_ref_simple_grid_v1_1_*.nc"))
    return pd.Series(files, index=pd.DatetimeIndex([pd.Timestamp(f.stem.rsplit("_", 1)[1]) for f in files])).sort_index()


def sea_level(cfg: Config, pts: pd.DataFrame) -> pd.DataFrame:
    """Disc-mean NASA-SSH sea level anomaly (cm) from the grid nearest in time on or before D0 - 3."""
    files = _ssh_files(cfg)
    rows, cache = [], {}
    for p in pts.itertuples():
        target = p.d0 - pd.Timedelta(days=LAG)
        before = files[files.index <= target]
        row = {"point": p.point, "ssha_cm": np.nan, "ssh_date": pd.NaT}
        if len(before) and (target - before.index[-1]).days <= 10:
            date = before.index[-1]
            if date not in cache:
                with xr.open_dataset(before.iloc[-1]) as ds:
                    v = ds["ssha"].squeeze()
                    cache[date] = (v.values.astype(np.float32), ds["latitude"].values, ds["longitude"].values,
                                   str(v.attrs.get("units", "m")))
            arr, la, lo, units = cache[date]
            iy, ix = disc(la, lo, p.lat, p.lon % 360, R_KM)
            vals = arr[iy, ix]
            scale = 100.0 if units.strip().lower() in ("m", "meter", "meters", "metre", "metres") else 1.0
            if np.isfinite(vals).any():
                row["ssha_cm"], row["ssh_date"] = float(np.nanmean(vals) * scale), date
        rows.append(row)
    return pd.DataFrame(rows)


def salinity(cfg: Config, pts: pd.DataFrame) -> pd.DataFrame:
    """Disc-mean SMAP salinity on D0 - 3 (secondary covariate); NaN where the file is not on disk."""
    folder = cfg.path("raw") / "smap_sss"
    rows, cache = [], {}
    for p in pts.itertuples():
        target = p.d0 - pd.Timedelta(days=LAG)
        f = folder / f"SMAP_L3_SSS_{target:%Y%m%d}_8DAYS_V5.0.nc"
        row = {"point": p.point, "sss": np.nan}
        if f.exists():
            if f not in cache:
                with xr.open_dataset(f) as ds:
                    cache[f] = (ds["smap_sss"].squeeze().values.astype(np.float32), ds["latitude"].values,
                                ds["longitude"].values)
            arr, la, lo = cache[f]
            iy, ix = disc(la, lo, p.lat, p.lon, R_KM)
            vals = arr[iy, ix]
            if np.isfinite(vals).any():
                row["sss"] = float(np.nanmean(vals))
        rows.append(row)
    return pd.DataFrame(rows)


def smap_dates_needed(pts: pd.DataFrame) -> list[pd.Timestamp]:
    """The SMAP files the salinity covariate needs (D0 - 3 of every point from 30 April 2015)."""
    d = (pts["d0"] - pd.Timedelta(days=LAG)).drop_duplicates()
    return sorted(d[d >= pd.Timestamp("2015-04-30")])


# ---------------------------------------------------------------- the table

def build(cfg: Config, refresh: bool = False) -> pd.DataFrame:
    out = cfg.path("interim") / "storms" / "storm_table.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    pts = track_points(cfg)
    fixes = all_fixes(cfg)
    table = pts
    for part in (ocean_state(cfg, pts, fixes), bloom(cfg, pts), pace_groups(cfg, pts), sea_level(cfg, pts),
                 salinity(cfg, pts), group_changes(cfg, pts)):
        table = table.merge(part, on="point", how="left")
    table["heatwave"] = np.where(table["E"] >= HOT, 1.0, np.where(table["E"] < NORMAL, 0.0, np.nan))
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(out, index=False)
    return table


# ---------------------------------------------------------------- the 2002-2011 replication

REPLICATION = ("2002-07-04", "2011-12-31")


def build_replication(cfg: Config, refresh: bool = False) -> pd.DataFrame:
    """The 2002-2011 storm table with MODIS-Aqua chlorophyll (pre-registration: replication section)."""
    out = cfg.path("interim") / "storms" / "storm_table_2002.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    pts = track_points(cfg, *REPLICATION)
    fixes = all_fixes(cfg)
    ocean = daily_ocean(cfg, pd.Timestamp("2001-11-01"), pd.Timestamp("2012-03-31"))
    table = pts
    for part in (ocean_state(cfg, pts, fixes, ocean), bloom(cfg, pts, "modis_aqua_chl"), sea_level(cfg, pts)):
        table = table.merge(part, on="point", how="left")
    table["heatwave"] = np.where(table["E"] >= HOT, 1.0, np.where(table["E"] < NORMAL, 0.0, np.nan))
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(out, index=False)
    return table


# ---------------------------------------------------------------- final test: any region, both eras

ERAS = {"modis": ("2002-07-04", "2011-12-31", "modis_aqua_chl", ("2001-11-01", "2012-03-31")),
        "viirs": ("2012-01-01", None, "viirs_snpp_chl", None)}


def build_region(cfg: Config, refresh: bool = False) -> pd.DataFrame:
    """Both eras for one region's config (used for the South Indian Ocean in the final test)."""
    out = cfg.path("interim") / "storms" / "storm_table_all.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    fixes = all_fixes(cfg)
    parts = []
    for era, (start, end, product, cube) in ERAS.items():
        pts = track_points(cfg, start, end)
        ocean = daily_ocean(cfg) if cube is None else daily_ocean(cfg, pd.Timestamp(cube[0]), pd.Timestamp(cube[1]))
        t = pts
        for part in (ocean_state(cfg, pts, fixes, ocean), bloom(cfg, pts, product), sea_level(cfg, pts)):
            t = t.merge(part, on="point", how="left")
        t["era"] = era
        parts.append(t)
        del ocean
    table = pd.concat(parts, ignore_index=True)
    table["point"] = np.arange(len(table))  # unique across eras
    table["heatwave"] = np.where(table["E"] >= HOT, 1.0, np.where(table["E"] < NORMAL, 0.0, np.nan))
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(out, index=False)
    return table


NOSTORM_KM = 500.0  # T4: a control year has no storm fix within this distance of the disc centre


def nostorm_controls(cfg: Config, pts: pd.DataFrame, product: str, ocean: xr.Dataset) -> pd.DataFrame:
    """T4 control bloom per point: the same disc, the same calendar 8-day periods, in other years of the same
    sensor era with no storm fix (any strength, both IBTrACS lists) within 500 km between the pre period's
    start and the second post period's end, and the same heatwave class on the same calendar day (D0 - 3)."""
    lat, lon, offshore = analysis_grid(cfg, product)
    periods = chl_periods(cfg, product)
    clim = chl_climatology(cfg, periods, offshore, "chl_slot_climatology.nc" if product == "viirs_snpp_chl"
                           else f"chl_slot_climatology_{product}.nc")
    starts, ends = pd.DatetimeIndex(periods["start"]), pd.DatetimeIndex(periods["end"])
    by_key = {(s.year, slot(s)): k for k, s in enumerate(starts)}
    years = sorted({s.year for s in starts})
    fixes = all_fixes(cfg, ("NI", "SI"))
    fx_t, fx_la, fx_lo = fixes["time"].to_numpy(), fixes["lat"].to_numpy(), fixes["lon"].to_numpy()
    otime = pd.DatetimeIndex(ocean["time"].values)
    oday = {d: i for i, d in enumerate(otime)}
    hot = ocean["hot"].values
    olat, olon = ocean["lat"].values, ocean["lon"].values
    is_ocean = np.isfinite(ocean["sst"].isel(time=slice(0, 30)).values).mean(axis=0) > 0.5
    idx = _disc_flat(lat, lon, pts, offshore, (R_KM,))
    needs, plan = {}, {}
    for p in pts.itertuples():
        if not np.isfinite(getattr(p, "heatwave", np.nan)):
            continue
        roles = pre_post(starts, ends, p.d0)
        if min(roles) < 0:
            continue
        cells = disc(olat, olon, p.lat, p.lon, R_KM, is_ocean)
        for y in years:
            dy = y - p.d0.year
            if dy == 0:
                continue
            ks = [by_key.get((starts[k].year + dy, slot(starts[k]))) for k in roles]
            if any(k is None for k in ks):
                continue
            t0, t1 = starts[ks[0]], ends[ks[2]] + pd.Timedelta(days=1)
            near = (fx_t >= np.datetime64(t0)) & (fx_t < np.datetime64(t1))
            if near.any() and (haversine_km(p.lat, p.lon, fx_la[near], fx_lo[near]) <= NOSTORM_KM).any():
                continue
            ke = oday.get((p.d0 - pd.Timedelta(days=LAG)) + pd.DateOffset(years=dy))
            if ke is None or not cells[0].size:
                continue
            e = float(hot[ke][cells].mean())
            cls = 1.0 if e >= HOT else 0.0 if e < NORMAL else np.nan
            if cls != p.heatwave:
                continue
            plan.setdefault(p.point, []).append(dy)
            for role, k in zip(("pre", "post1", "post2"), ks):
                needs.setdefault(k, []).append((p.point, R_KM, f"{role}@{dy}"))

    def fields():
        for k in sorted(needs):
            row = periods.iloc[k]
            with xr.open_dataset(row.path) as ds:
                x = ds["chlor_a"].values.astype(np.float32)
            with np.errstate(divide="ignore", invalid="ignore"):
                yield k, np.where((x > 0) & offshore, np.log10(x) - clim[slot(row.start)], np.nan)

    vals = _period_values(fields(), needs, idx)
    rows = []
    for p in pts.itertuples():
        bs = []
        for dy in plan.get(p.point, []):
            pre = vals.get((p.point, R_KM, f"pre@{dy}"), np.nan)
            post = [vals.get((p.point, R_KM, f"{r}@{dy}"), np.nan) for r in ("post1", "post2")]
            good = [v for v in post if np.isfinite(v)]
            if good and np.isfinite(pre):
                bs.append(np.mean(good) - pre)
        rows.append({"point": p.point, "bloom_control": float(np.mean(bs)) if bs else np.nan,
                     "n_control": len(bs)})
    return pd.DataFrame(rows)


def add_controls(cfg: Config, table_path, eras: dict) -> pd.DataFrame:
    """Add T4's bloom_control to a cached storm table, era by era ({era: (product, cube start, cube end)})."""
    table = pd.read_parquet(table_path)
    table = table[[c for c in table.columns if c not in ("bloom_control", "n_control")]]
    parts = []
    for era, (product, c0, c1) in eras.items():
        sub = table if era is None else table[table["era"] == era]
        if sub.empty:
            continue
        ocean = daily_ocean(cfg) if c0 is None else daily_ocean(cfg, pd.Timestamp(c0), pd.Timestamp(c1))
        parts.append(nostorm_controls(cfg, sub, product, ocean))
        del ocean
    table = table.merge(pd.concat(parts, ignore_index=True), on="point", how="left")
    table.to_parquet(table_path, index=False)
    return table
