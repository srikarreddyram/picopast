"""Run the pre-registered fallback checks A (pigment proxy) and B (phycoerythrin, exploratory).

Inputs: SeaBASS HPLC and flow-through (CLS) files in ``data/insitu/seabass/``
and daily MOANA maps from ``picopast moana-run``. Output: a Markdown report
in ``docs/reports/`` and the matchup tables in ``data/interim/``.

The rules (matchup, masks, metrics, verdicts) are those of
``docs/trust_gate_preregistration.md``; nothing here may change them.
"""

from __future__ import annotations

import glob
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.fallback import pigment_proxy_verdict, spearman_bootstrap
from picopast.grid import coast_mask, nearest_to_grid
from picopast.insitu.pigments import load_hplc, station_days
from picopast.insitu.seabass import read_sb
from picopast.matchup import match

SENSITIVITY = [(1, 1), (3, 5), (5, 13)]  # (window_px, min_valid_px): pixel, rule, wider


def _moana_days(cfg: Config, mapping: str, days: pd.DatetimeIndex) -> xr.Dataset | None:
    paths = []
    for d in days:
        p = cfg.path("interim") / "moana" / mapping / "DAY" / f"moana_{d:%Y%m%d}.nc"
        if p.exists():
            paths.append(p)
    return xr.open_mfdataset(paths, combine="by_coords").load() if paths else None


def _offshore(cfg: Config, maps: xr.Dataset) -> np.ndarray:
    """Pixels at least the coast-mask distance from land, land taken from CMC's 0.1 degree mask."""
    cmc = sorted((cfg.path("raw") / "cmc_sst").glob("*CMC*.box.nc"))
    with xr.open_dataset(cmc[0]) as ds:
        water = np.isfinite(nearest_to_grid(ds["analysed_sst"].isel(time=0), maps["lat"].values,
                                            maps["lon"].values).transpose("lat", "lon").values)
    return coast_mask(water, maps["lat"].values, cfg["grid"]["resolution_deg"], cfg["grid"]["coast_mask_km"])


def _valid(maps: xr.Dataset, group: str, offshore: np.ndarray) -> xr.DataArray:
    """Pre-registered pixel rule: MOANA flags zero (a negative Prochlorococcus is flagged), offshore."""
    return maps[group].where((maps["flags"] == 0) & offshore)


def _station_points(sd: pd.DataFrame, offshore: np.ndarray, maps: xr.Dataset) -> pd.DataFrame:
    lat, lon = maps["lat"].values, maps["lon"].values
    iy = np.abs(lat[None, :] - sd["lat"].values[:, None]).argmin(1)
    ix = np.abs(lon[None, :] - sd["lon"].values[:, None]).argmin(1)
    return sd.assign(offshore=offshore[iy, ix])


def check_a(cfg: Config, mapping: str = "operational") -> dict:
    """Check A: MOANA Prochlorococcus vs HPLC DVChl a at PACE-era surface station-days."""
    m_cfg, rule = cfg["matchup"], cfg["trust_gate"]["pigment_proxy"]
    samples = load_hplc(sorted(Path(p) for p in glob.glob(str(cfg.path("insitu") / "seabass" / "*HPLC*.sb"))))
    sd = station_days(samples, cfg.box, start=cfg["periods"]["pace"]["start"],
                      max_depth_m=m_cfg["max_depth_m"], resolution_deg=cfg["grid"]["resolution_deg"])
    days = pd.DatetimeIndex(sorted(set(sd["date"]) | set(sd["date"] - pd.Timedelta(days=1))
                                   | set(sd["date"] + pd.Timedelta(days=1))))
    maps = _moana_days(cfg, mapping, days)
    if maps is None:
        raise FileNotFoundError(f"no {mapping} MOANA maps for the station days; run picopast moana-run")
    offshore = _offshore(cfg, maps)
    sd = _station_points(sd, offshore, maps)
    stations = sd[sd["offshore"]]
    pro = _valid(maps, "pro", offshore)

    res = {"mapping": mapping, "n_station_days": int(len(sd)), "n_coastal_excluded": int((~sd["offshore"]).sum()),
           "maps_available": int(maps.sizes["time"]), "days_needed": int(len(days))}
    rows = []
    for window, min_valid in SENSITIVITY:
        m = match(stations, pro, window, min_valid, tuple(m_cfg["day_offsets"]))
        s = spearman_bootstrap(m["sat_value"], m["dv_chl_a"], m["station_day"],
                               n_boot=cfg["trust_gate"]["n_bootstrap"])
        rows.append({"window_px": window, "min_valid_px": min_valid, **s})
        if (window, min_valid) == (m_cfg["window_px"], m_cfg["min_valid_px"]):
            res["stats"] = s
            res["verdict"], res["reasons"] = pigment_proxy_verdict(s, rule)
            res["matchups"] = m
            res["day_offsets"] = m["sat_day_offset"].value_counts(dropna=False).to_dict()
            z = m[["sat_value", "zea"]].dropna()
            res["zea_rho_descriptive"] = spearman_bootstrap(z["sat_value"], z["zea"], m.loc[z.index, "station_day"],
                                                            n_boot=cfg["trust_gate"]["n_bootstrap"])["rho"]
            m_un = m[m["sat_value"].isna()]
            res["unmatched_by_cruise"] = m_un.groupby("cruise").size().to_dict()
    res["sensitivity"] = pd.DataFrame(rows)
    retrieved = maps["pro"].where(offshore)
    n_ret = int(np.isfinite(retrieved).sum())
    res["negative_pro_share"] = float((retrieved < 0).sum() / n_ret) if n_ret else float("nan")
    return res


def check_b(cfg: Config, mapping: str = "operational") -> dict:
    """Check B (exploratory, no verdict): MOANA Synechococcus vs flow-through phycoerythrin Fpe565."""
    m_cfg, res_deg = cfg["matchup"], cfg["grid"]["resolution_deg"]
    frames = []
    for p in sorted(glob.glob(str(cfg.path("insitu") / "seabass" / "*_CLS.sb"))):
        header, df = read_sb(p)
        df = df.assign(cruise=header.get("cruise", Path(p).stem),
                       time_utc=pd.to_datetime(df["date"].astype("Int64").astype(str) + " " + df["time"].astype(str),
                                               format="%Y%m%d %H:%M:%S", errors="coerce"))
        frames.append(df[["cruise", "time_utc", "lat", "lon", "Fpe565"]])
    ft = pd.concat(frames, ignore_index=True).dropna()
    ft = ft[(ft["time_utc"] >= pd.Timestamp(cfg["periods"]["pace"]["start"]))
            & ft["lat"].between(cfg.box.lat_min, cfg.box.lat_max)
            & ft["lon"].between(cfg.box.lon_min, cfg.box.lon_max)]
    ft = ft.assign(date=ft["time_utc"].dt.normalize(),
                   iy=np.floor((90 - ft["lat"]) / res_deg).astype(int),
                   ix=np.floor((ft["lon"] + 180) / res_deg).astype(int))
    pts = ft.groupby(["date", "iy", "ix"], as_index=False).agg(
        Fpe565=("Fpe565", "median"), lat=("lat", "mean"), lon=("lon", "mean"),
        time_utc=("time_utc", "min"), cruise=("cruise", "first"))
    pts["station_day"] = pts["date"].dt.strftime("%Y%m%d") + "_" + pts["iy"].astype(str) + "_" + pts["ix"].astype(str)
    days = pd.DatetimeIndex(sorted(set(pts["date"]) | set(pts["date"] - pd.Timedelta(days=1))
                                   | set(pts["date"] + pd.Timedelta(days=1))))
    maps = _moana_days(cfg, mapping, days)
    if maps is None:
        raise FileNotFoundError(f"no {mapping} MOANA maps for the flow-through days")
    offshore = _offshore(cfg, maps)
    pts = _station_points(pts, offshore, maps)
    m = match(pts[pts["offshore"]], _valid(maps, "syn", offshore), m_cfg["window_px"],
              m_cfg["min_valid_px"], tuple(m_cfg["day_offsets"]))
    s = spearman_bootstrap(m["sat_value"], m["Fpe565"], m["station_day"], n_boot=cfg["trust_gate"]["n_bootstrap"])
    return {"mapping": mapping, "n_pixel_days": int(len(pts)), "stats": s, "matchups": m}


def _fmt(s: dict) -> str:
    return f"ρ = {s['rho']:.2f} [{s['lo']:.2f}, {s['hi']:.2f}], n = {s['n']}"


def report(cfg: Config, mappings=("operational", "atbd")) -> Path:
    """Run A and B under each mapping that has maps, and write the report."""
    lines = ["# Trust-gate fallback checks A and B", "",
             f"Generated {date.today().isoformat()} by `picopast fallback-checks`, under the rule fixed in "
             "`docs/trust_gate_preregistration.md` (16 November 2026). Check A can support Prochlorococcus "
             "only as a relative index, never as a cells/mL label; check B is exploratory and gives no verdict.", ""]
    out_dir = cfg.path("interim")
    for mapping in mappings:
        try:
            a = check_a(cfg, mapping)
        except FileNotFoundError as e:
            lines += [f"## {mapping} mapping", "", f"Not run: {e}", ""]
            continue
        a["matchups"].to_csv(out_dir / f"matchups_check_a_{mapping}.csv", index=False)
        lines += [f"## {mapping} mapping", "", "### A. Prochlorococcus vs HPLC divinyl chlorophyll a", "",
                  f"**Verdict: {a['verdict']}**" + (f" ({'; '.join(a['reasons'])})" if a["reasons"] else ""), "",
                  f"- Station-days: {a['n_station_days']} ({a['n_coastal_excluded']} within "
                  f"{cfg['grid']['coast_mask_km']:.0f} km of land, excluded); MOANA maps found for "
                  f"{a['maps_available']} of {a['days_needed']} days needed.",
                  f"- Spearman, MOANA Prochlorococcus vs DVChl a: {_fmt(a['stats'])} station-days matched "
                  "(95% interval, station-day bootstrap).",
                  f"- Day used: {a['day_offsets']} (0 = same day; NaN = unmatched). Unmatched by cruise: "
                  f"{a['unmatched_by_cruise']}.",
                  f"- Zeaxanthin (descriptive only): ρ = {a['zea_rho_descriptive']:.2f}.",
                  f"- Share of negative Prochlorococcus retrievals among offshore retrieved pixels on these days: "
                  f"{a['negative_pro_share']:.1%}.", "", "Sensitivity to the block rule:", "",
                  "| window | min valid | ρ | interval | n |", "| --- | --- | --- | --- | --- |"]
        for r in a["sensitivity"].itertuples():
            lines.append(f"| {r.window_px}×{r.window_px} | {r.min_valid_px} | {r.rho:.2f} | [{r.lo:.2f}, {r.hi:.2f}] | {r.n} |")
        try:
            b = check_b(cfg, mapping)
            b["matchups"].to_csv(out_dir / f"matchups_check_b_{mapping}.csv", index=False)
            lines += ["", "### B. Synechococcus vs flow-through phycoerythrin (Fpe565), exploratory", "",
                      f"- Pixel-days with flow-through data: {b['n_pixel_days']}; {_fmt(b['stats'])} matched.", ""]
        except FileNotFoundError as e:
            lines += ["", f"### B. Not run: {e}", ""]
    out = cfg.path("reports") / "trust_gate_fallback_ab.md"
    out.write_text("\n".join(lines) + "\n")
    return out
