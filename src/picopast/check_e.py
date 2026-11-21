"""Check E: new, independent HPLC stations scored against Route D.

Rules: ``docs/check_e_preregistration.md`` (binding, committed before this ran).

* E1, pattern: the D2 same-calendar-month field vs every new station, matched as in check C.
* E2, record: stations dated 2012 or later vs the SNPP Route D record for the station's own
  8-day period (all valid pixels primary; reliable pixels only as a sensitivity).

Divinyl chlorophyll a scores *Prochlorococcus* (verdict). Zeaxanthin is shared by both
cyanobacteria, so the *Synechococcus* rows are descriptive only.
"""

from __future__ import annotations

import glob
import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.fallback import climatology_verdict, paired_spearman_difference, spearman_bootstrap
from picopast.fallback_run import _chl_monthly_climatology, _monthly_climatology, _station_points
from picopast.insitu.pigments import load_hplc, station_days
from picopast.matchup import match

#: (group, pigment, gives a verdict, chlorophyll sign fixed in check C)
TESTS = (("pro", "dv_chl_a", True, -1.0), ("syn", "zea", False, 1.0))
NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus"}
PIGMENT_NAMES = {"dv_chl_a": "divinyl chlorophyll a", "zea": "zeaxanthin"}
RECORD_START = pd.Timestamp("2012-01-01")


def stations(cfg: Config) -> pd.DataFrame:
    """Surface station-days with HPLC pigments in the box, all years (new files plus the earlier PACE-era ones)."""
    ins = cfg.path("insitu")
    paths = sorted({Path(p) for p in glob.glob(str(ins / "seabass_e" / "*.sb"))}
                   | {Path(p) for p in glob.glob(str(ins / "seabass" / "*HPLC*.sb"))})
    return station_days(load_hplc(paths), cfg.box, start="1990-01-01", max_depth_m=cfg["matchup"]["max_depth_m"])


def period_start(t: pd.Timestamp) -> pd.Timestamp:
    """Start of the 8-day period (restarting each 1 January, as the L3 8-day products do) holding ``t``."""
    t = pd.Timestamp(t).normalize()
    return pd.Timestamp(t.year, 1, 1) + pd.Timedelta(days=8 * ((t.dayofyear - 1) // 8))


def verdict(stats: dict, rule: dict) -> str:
    if stats["n"] < rule["min_stations"] or stats["n_clusters"] < rule["min_cruises"]:
        return "not enough data"
    return climatology_verdict(stats, rule)[0]


def _same_month(s: pd.DataFrame, clim: xr.DataArray, m_cfg: dict) -> np.ndarray:
    """check C's lookup: same calendar month, 3x3 median, no day offsets."""
    field = clim.rename(month="time").assign_coords(
        time=pd.to_datetime([f"2000-{int(m):02d}-01" for m in clim["month"].values]))
    t = s.assign(time_utc=pd.to_datetime("2000-" + s["time_utc"].dt.strftime("%m") + "-01"))
    return match(t, field, m_cfg["window_px"], m_cfg["min_valid_px"], (0,))["sat_value"].to_numpy()


def _same_period(s: pd.DataFrame, field: xr.DataArray, m_cfg: dict) -> np.ndarray:
    t = s.assign(time_utc=s["period"])
    return match(t, field, m_cfg["window_px"], m_cfg["min_valid_px"], (0,))["sat_value"].to_numpy()


def _stack(starts, load) -> xr.DataArray:
    maps, times, coords = [], [], None
    for p in starts:
        a = load(p)
        if a is None:
            continue
        maps.append(a.values.astype(np.float32))
        times.append(p)
        coords = a
    return xr.DataArray(np.stack(maps), dims=("time", "lat", "lon"),
                        coords={"time": times, "lat": coords["lat"].values, "lon": coords["lon"].values})


def record_fields(cfg: Config, starts, offshore: np.ndarray) -> dict[str, xr.DataArray]:
    """Per-period fields at the needed 8-day starts: Route D (all / reliable), Route A, VIIRS log10 chl."""
    from picopast.route_d_reconstruct import chl_par_periods

    prod = cfg.path("products") / "reconstruction"
    chl_paths = dict(zip(chl_par_periods(cfg, "snpp")["start"], chl_par_periods(cfg, "snpp")["chl"]))

    def route_d(var, reliable=False):
        def load(p):
            f = prod / "snpp_route_d" / "8D" / f"recon_d_{p:%Y%m%d}.nc"
            if not f.exists():
                return None
            with xr.open_dataset(f) as ds:
                a = ds[var].isel(time=0).load()
                return a.where(ds[f"{var}_reliable"].isel(time=0) == 1) if reliable else a
        return load

    def route_a(var):
        def load(p):
            f = prod / "snpp" / "8D" / f"recon_{p:%Y%m%d}.nc"
            if not f.exists():
                return None
            with xr.open_dataset(f) as ds:
                return ds[var].isel(time=0).load()
        return load

    def chl(p):
        f = chl_paths.get(p)
        if f is None:
            return None
        with xr.open_dataset(f) as ds:
            x = ds["chlor_a"].values.astype(np.float64)
            like = ds["chlor_a"]
            with np.errstate(divide="ignore", invalid="ignore"):
                v = np.where((x > 0) & offshore, np.log10(x), np.nan)
            return xr.DataArray(v, dims=("lat", "lon"), coords={"lat": like["lat"].values, "lon": like["lon"].values})

    out = {}
    for g in ("pro", "syn"):
        out[f"D_{g}"] = _stack(starts, route_d(g))
        out[f"Drel_{g}"] = _stack(starts, route_d(g, reliable=True))
        out[f"A_{g}"] = _stack(starts, route_a(g))
    out["chl"] = _stack(starts, chl)
    return out


def _stat(s: pd.DataFrame, col: str, pig: str, mask, n_boot: int) -> dict:
    m = np.ones(len(s), bool) if mask is None else mask
    return spearman_bootstrap(s[col].to_numpy()[m], s[pig].to_numpy()[m], s["cruise"].to_numpy()[m], n_boot=n_boot)


def _paired(s: pd.DataFrame, a: str, b: str, pig: str, mask, n_boot: int) -> dict:
    m = np.ones(len(s), bool) if mask is None else mask
    return paired_spearman_difference(s[a].to_numpy()[m], s[b].to_numpy()[m], s[pig].to_numpy()[m],
                                      s["cruise"].to_numpy()[m], n_boot=n_boot)


def run(cfg: Config) -> dict:
    from picopast.route_d import box_fields

    m_cfg, rule = cfg["matchup"], cfg["trust_gate"]["climatology_pattern"]
    n_boot = cfg["trust_gate"]["n_bootstrap"]
    work = cfg.path("interim") / "check_e"
    work.mkdir(parents=True, exist_ok=True)
    d2 = {g: xr.open_dataarray(cfg.path("interim") / "route_d" / f"field_D2_{g}.nc").load() for g in ("pro", "syn")}
    offshore = box_fields(cfg)["offshore"]
    sd = _station_points(stations(cfg), offshore, d2["pro"])
    res = {"n_station_days": int(len(sd)), "n_coastal_excluded": int((~sd["offshore"]).sum())}
    sd = sd[sd["offshore"]].copy()
    sd["period"] = sd["time_utc"].map(period_start)
    moana = _monthly_climatology(cfg, "operational", offshore)[0]
    chl_clim = _chl_monthly_climatology(cfg, d2["pro"], offshore)
    starts = sorted(sd.loc[sd["time_utc"] >= RECORD_START, "period"].unique())
    rec = record_fields(cfg, starts, offshore)

    for g, pig, gives_verdict, sign in TESTS:
        s = sd[sd[pig] > 0].copy()
        s["clim_D2"] = _same_month(s, d2[g], m_cfg)
        s["clim_chl"] = sign * _same_month(s, chl_clim, m_cfg)
        s["clim_moana"] = _same_month(s, moana[g], m_cfg)
        e2 = s["time_utc"] >= RECORD_START
        r = s[e2].copy()
        r["rec_D"] = _same_period(r, rec[f"D_{g}"], m_cfg)
        r["rec_Drel"] = _same_period(r, rec[f"Drel_{g}"], m_cfg)
        r["rec_A"] = _same_period(r, rec[f"A_{g}"], m_cfg)
        r["rec_chl"] = sign * _same_period(r, rec["chl"], m_cfg)
        s = s.join(r[["rec_D", "rec_Drel", "rec_A", "rec_chl"]])
        s.assign(group=g, pigment=pig).to_csv(work / f"matchups_{g}.csv", index=False)
        e2m = e2.to_numpy()

        def st(col, mask=None, s=s, pig=pig):
            return _stat(s, col, pig, mask, n_boot)

        def paired(a, b, mask=None, s=s, pig=pig):
            return _paired(s, a, b, pig, mask, n_boot)

        out = {
            "E1": {"D2": st("clim_D2"), "chl": st("clim_chl"), "moana": st("clim_moana"),
                   "D2_vs_chl": paired("clim_D2", "clim_chl"), "D2_vs_moana": paired("clim_D2", "clim_moana")},
            "E2": {"record": st("rec_D", e2m), "record_reliable": st("rec_Drel", e2m), "chl": st("rec_chl", e2m),
                   "route_a": st("rec_A", e2m), "record_vs_chl": paired("rec_D", "rec_chl", e2m),
                   "record_vs_route_a": paired("rec_D", "rec_A", e2m),
                   "record_given_clim": paired("rec_D", "clim_D2", e2m)},
            "n_stations": int(len(s)), "n_e2_stations": int(e2m.sum()),
            "cruises": s.groupby("cruise").size().to_dict(),
        }
        if gives_verdict:
            out["E1"]["verdict"] = verdict(out["E1"]["D2"], rule)
            out["E2"]["verdict"] = verdict(out["E2"]["record"], rule)
            primary = "E2" if out["E2"]["verdict"] != "not enough data" else "E1"
            out["primary"] = primary
            out["decision"] = {"pattern consistent": "supported: the Route D record stays a relative index",
                               "pattern not consistent": "downgraded to exploratory",
                               }.get(out[primary]["verdict"], "no change (reported as such)")
        res[g] = out
    (work / "results.json").write_text(json.dumps(res, indent=1, default=float))
    return res


def _f(s: dict) -> str:
    return f"{s['rho']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] (n {s['n']}, {s['n_clusters']} cruises)"


def _p(s: dict) -> str:
    return f"{s['diff']:+.2f} [{s['diff_lo']:+.2f}, {s['diff_hi']:+.2f}]"


def _partial(s: dict) -> str:
    return f"{s['partial']:+.2f} [{s['partial_lo']:+.2f}, {s['partial_hi']:+.2f}]"


def report(cfg: Config) -> Path:
    res = run(cfg)
    L = ["# Check E: new, independent HPLC stations vs Route D", "",
         f"Generated {date.today().isoformat()} by `picopast check-e`. Rules: `docs/check_e_preregistration.md`, "
         "committed before this ran. Divinyl chlorophyll a (made only by *Prochlorococcus*) scores "
         "*Prochlorococcus*; zeaxanthin, made by both cyanobacteria, is shown for *Synechococcus* as a "
         "description only. Pigments are not cell counts: these are rank tests.", "",
         f"Station-days (surface, one UTC day × one 4.6 km pixel): {res['n_station_days']}, of which "
         f"{res['n_coastal_excluded']} coastal (excluded by check C's offshore rule).", ""]
    for g, pig, gives_verdict, _ in TESTS:
        r = res[g]
        label = f"{NAMES[g]} vs {PIGMENT_NAMES[pig]}" + ("" if gives_verdict else " (descriptive only)")
        L += [f"## {label}", "",
              f"{r['n_stations']} offshore station-days, {r['n_e2_stations']} dated 2012 or later. Per cruise: "
              + ", ".join(f"{k} {v}" for k, v in sorted(r["cruises"].items())) + ".", "",
              "| test | model | Spearman ρ [95%, cruises resampled] |", "| --- | --- | --- |"]
        L += [f"| E1 pattern (same month, all dates) | Route D (D2) | {_f(r['E1']['D2'])} |",
              f"| | chlorophyll alone | {_f(r['E1']['chl'])} |",
              f"| | MOANA | {_f(r['E1']['moana'])} |",
              f"| E2 record (own 8-day period, 2012 on) | **Route D record, all valid pixels** | "
              f"{_f(r['E2']['record'])} |",
              f"| | Route D record, reliable pixels only | {_f(r['E2']['record_reliable'])} |",
              f"| | VIIRS 8-day chlorophyll alone | {_f(r['E2']['chl'])} |",
              f"| | old MOANA-learned record (Route A) | {_f(r['E2']['route_a'])} |", ""]
        L += ["| paired comparison, same stations | difference in ρ [95%] | partial ρ of the first given the second |",
              "| --- | --- | --- |",
              f"| E1: D2 − chlorophyll | {_p(r['E1']['D2_vs_chl'])} | {_partial(r['E1']['D2_vs_chl'])} |",
              f"| E1: D2 − MOANA | {_p(r['E1']['D2_vs_moana'])} | {_partial(r['E1']['D2_vs_moana'])} |",
              f"| E2: record − chlorophyll | {_p(r['E2']['record_vs_chl'])} | {_partial(r['E2']['record_vs_chl'])} |",
              f"| E2: record − Route A | {_p(r['E2']['record_vs_route_a'])} | "
              f"{_partial(r['E2']['record_vs_route_a'])} |",
              f"| E2: record − its own same-month climatology | {_p(r['E2']['record_given_clim'])} | "
              f"{_partial(r['E2']['record_given_clim'])} |", ""]
        if gives_verdict:
            L += [f"**Verdicts:** E1 *{r['E1']['verdict']}*; E2 *{r['E2']['verdict']}*. Primary test "
                  f"({r['primary']}) → **{r['decision']}**.", ""]
    L += ["## Limits", "",
          "- Pigment per cell changes with light and nutrients, so a pigment rank is not a cell-count rank.",
          "- Stations come from a handful of cruises; intervals resample whole cruises and are wide.",
          "- Same-period matching loses stations under cloud; there is no gap filling.",
          "- No new *Synechococcus* or *Prochlorococcus* cell counts exist in the open archives searched "
          "(PANGAEA, BCO-DMO, gVOD, SeaBASS), so *Synechococcus* cannot be re-tested here.", ""]
    out = cfg.path("reports") / "check_e.md"
    out.write_text("\n".join(L))
    return out
