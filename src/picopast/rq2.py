"""RQ2: how much group information survives in VIIRS bands? (PRD section 3 and section 6 step 4)

Identical models predict our MOANA output from two views of the same PACE
spectrum: (a) the full spectrum (172 bands, 346-719 nm) and (b) that spectrum
band-averaged through SNPP's M1-M5 responses (virtual VIIRS). Ancillary inputs
are the same in both (CMC SST, position, season). Because both views come from
one spectrum, atmosphere and calibration are identical: the skill gap is the
information the five broad bands lose. Real-VIIRS differences are separate
(docs/reports/virtual_viirs.md).

PACE-era 8-day maps, offshore pixels that pass the reconstruction's screen,
sampled per period; scored like the reconstruction (spatial blocks, temporal
hold-out, skill against a seasonal climatology). One period in memory at a time.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.reconstruct import (
    GROUPS,
    SeasonalClimatology,
    _lgbm,
    _ridge,
    _scores,
    to_target,
    usable_pixels,
)

NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes"}
ANCILLARY = ["sst", "lat", "lon", "doy_sin", "doy_cos"]


def _ratio(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((a > 0) & (b > 0), np.log10(a / b), np.nan).astype(np.float32)


def table(cfg: Config, per_period: int = 5000, seed: int = 0, refresh: bool = False) -> pd.DataFrame:
    """Sampled pixels with the full spectrum, the virtual-VIIRS view, ancillaries and MOANA targets."""
    from picopast.fallback_run import _offshore
    from picopast.grid import nearest_to_grid
    from picopast.ingest.rsr import load_rsr
    from picopast.moana.box_run import _cmc_files, scene_period, sst_for_period
    from picopast.virtual_viirs import VIIRS_BANDS, band_average

    out = cfg.path("interim") / "rq2" / "table.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    rsr, bands = load_rsr(cfg, "snpp"), VIIRS_BANDS["snpp"]
    cmc, rng, parts, offshore = _cmc_files(cfg), np.random.default_rng(seed), [], None
    moana_dir = cfg.path("interim") / "moana" / cfg["reconstruction"]["mapping"] / "8D"
    for path in sorted((cfg.path("raw") / "pace_rrs").glob("PACE_OCI.*.L3m.8D.*.box.nc")):
        start, end, _ = scene_period(path)
        mpath = moana_dir / f"moana_{start:%Y%m%d}.nc"
        if not mpath.exists():
            continue
        with xr.open_dataset(mpath) as m:
            cells = {g: m[g].values[0] for g in GROUPS}
            usable = usable_pixels(m["flags"].values[0], cells, cfg["reconstruction"]["valid_range"])
        sst_c = sst_for_period(cfg, start, end, cmc)
        with xr.open_dataset(path) as ds:
            if offshore is None:
                offshore = _offshore(cfg, ds)
                lat, lon = ds["lat"].values, ds["lon"].values
                lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
                wl = ds["wavelength"].values
            rrs = ds["Rrs"].values  # (lat, lon, 172) float32, ~0.5 GB
        sst = nearest_to_grid(sst_c, lat, lon).transpose("lat", "lon").values.astype(np.float32)
        ok = usable & offshore & np.isfinite(rrs).all(-1) & np.isfinite(sst)
        idx = np.flatnonzero(ok.ravel())
        if idx.size > per_period:
            idx = np.sort(rng.choice(idx, per_period, replace=False))
        if not idx.size:
            continue
        spec = rrs.reshape(-1, wl.size)[idx]
        vv = band_average(wl, spec, rsr).astype(np.float32)
        mid = start + (end - start) / 2
        doy = 2 * np.pi * (mid.day_of_year - 0.5) / 365.25
        d = {f"p{w:.0f}_{i}": spec[:, i] for i, w in enumerate(wl)}
        d |= {f"v{b}": vv[:, j] for j, b in enumerate(bands)}
        d |= {"v_r443_551": _ratio(vv[:, 1], vv[:, 3]), "v_r486_551": _ratio(vv[:, 2], vv[:, 3]),
              "v_r410_443": _ratio(vv[:, 0], vv[:, 1]), "v_r671_551": _ratio(vv[:, 4], vv[:, 3]),
              "v_mbr": _ratio(np.maximum(vv[:, 1], vv[:, 2]), vv[:, 3])}
        d |= {"sst": sst.ravel()[idx], "lat": lat2.ravel()[idx].astype(np.float32),
              "lon": lon2.ravel()[idx].astype(np.float32), "doy_sin": np.full(idx.size, np.sin(doy), np.float32),
              "doy_cos": np.full(idx.size, np.cos(doy), np.float32)}
        # Real VIIRS SNPP bands of the same period on the same pixels, for Route B trained on real matchups.
        # Extra columns only: the sample is drawn before, so the RQ2 rows do not depend on them.
        try:
            from picopast.reconstruct import BANDS, load_viirs

            real = load_viirs(cfg, "snpp", start)
            d |= {f"r{b}": real[f"rrs{b}"].ravel()[idx] for b in BANDS}
        except ValueError:  # no SNPP 8-day file for this period (archive gap)
            d |= {f"r{b}": np.full(idx.size, np.nan, np.float32) for b in VIIRS_BANDS["snpp"]}
        df = pd.DataFrame(d)
        for g in GROUPS:
            df[g] = to_target(g, cells[g].ravel()[idx]).astype(np.float32)
        df["time"], df["month"] = start, mid.month
        parts.append(df)
        del rrs
    t = pd.concat(parts, ignore_index=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(out, index=False)
    return t


def views(t: pd.DataFrame) -> dict[str, list[str]]:
    full = [c for c in t.columns if c.startswith("p") and "_" in c]
    virtual = [c for c in t.columns if c.startswith("v")]
    return {"full PACE spectrum": full + ANCILLARY, "virtual VIIRS (5 bands)": virtual + ANCILLARY}


def evaluate(cfg: Config, t: pd.DataFrame) -> pd.DataFrame:
    from picopast.splits import spatial_block_folds, temporal_holdout

    rows = []
    for g in GROUPS:
        d = t[np.isfinite(t[g])].reset_index(drop=True)
        y = d[g].values
        splits = {"spatial": list(spatial_block_folds(d["lat"].values, d["lon"].values, 5,
                                                      cfg["grid"]["block_deg"])),
                  "temporal": [temporal_holdout(d["time"].values, cfg["periods"]["temporal_holdout_months"])]}
        for kind, folds in splits.items():
            clim = np.full(y.size, np.nan)
            preds = {(v, m): np.full(y.size, np.nan) for v in views(t) for m in ("ridge", "LightGBM")}
            for tr, te in folds:
                clim[te] = SeasonalClimatology().fit(d.iloc[tr], y[tr]).predict(d.iloc[te])
                for v, cols in views(t).items():
                    X = d[cols]
                    preds[(v, "ridge")][te] = _ridge().fit(X.iloc[tr], y[tr]).predict(X.iloc[te])
                    preds[(v, "LightGBM")][te] = _lgbm().fit(X.iloc[tr], y[tr]).predict(X.iloc[te])
            m = np.isfinite(clim)
            for (v, model), p in preds.items():
                rows.append({"group": g, "split": kind, "inputs": v, "model": model, "n": int(m.sum()),
                             **_scores(y[m], p[m], clim[m])})
    return pd.DataFrame(rows)


def report(cfg: Config, refresh: bool = False) -> Path:
    t = table(cfg)
    scores_path = cfg.path("interim") / "rq2" / "scores.csv"
    if refresh or not scores_path.exists():
        evaluate(cfg, t).to_csv(scores_path, index=False)  # ~80 min, ~6 GB on the 16 GB machine
    scores = pd.read_csv(scores_path)
    L = ["# RQ2: how much group information survives in VIIRS bands?", "",
         f"Generated {date.today().isoformat()} by `picopast rq2`. Identical models predict our MOANA output "
         "(operational mapping) from two views of the same PACE spectrum: all 172 bands, or the spectrum "
         "band-averaged through SNPP's M1–M5 responses (virtual VIIRS, with the same band ratios the reconstruction "
         "uses). Both views get the same SST, position and season inputs. Same atmosphere and calibration in both, "
         "so the gap is the information the broad bands lose.", "",
         f"- {t['time'].nunique()} PACE-era 8-day periods, {len(t):,} sampled offshore pixels (the reconstruction's "
         "screen: MOANA flags 0 and inside NASA's valid range).",
         "- Skill = 1 − RMSE / RMSE of a seasonal climatology fitted in the same training part; spatial = 2° blocks "
         "in 5 folds, pooled; temporal = the last 6 months held out.", "",
         "| group | split | model | skill, full PACE | skill, virtual VIIRS | lost | ρ, full | ρ, virtual |",
         "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        for split in ("spatial", "temporal"):
            for model in ("LightGBM", "ridge"):
                s = scores[(scores["group"] == g) & (scores["split"] == split) & (scores["model"] == model)]
                f = s[s["inputs"] == "full PACE spectrum"].iloc[0]
                v = s[s["inputs"] == "virtual VIIRS (5 bands)"].iloc[0]
                L.append(f"| {NAMES[g]} | {split} | {model} | {f.skill:.2f} | {v.skill:.2f} | "
                         f"{f.skill - v.skill:+.2f} | {f.spearman:.2f} | {v.spearman:.2f} |")
    L += ["", "\"Lost\" is the skill the five bands give up against the full spectrum (larger = more group "
          "information lost). The PRD's hypothesis: the biggest loss is for Synechococcus.", ""]
    real = cfg.path("interim") / "reconstruct" / "validation_snpp_operational.csv"
    if real.exists():
        r = pd.read_csv(real)
        L += ["## Where the skill goes, PACE to real VIIRS (LightGBM)", "",
              "Three steps on the same PACE-era 8-day periods and the same screen: the full spectrum, the five "
              "bands simulated from that spectrum, and real VIIRS SNPP (the reconstruction's Route A, which also "
              "has PAR). The first gap is the bands; the second is everything that differs between two sensors "
              "looking at the ocean: overpass time, which days are cloud-free within the 8-day composite, "
              "calibration and atmospheric correction.", "",
              "| group | split | full PACE | virtual VIIRS | real VIIRS | lost to bands | lost to sensor mismatch |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
        for g in GROUPS:
            for split in ("spatial", "temporal"):
                s = scores[(scores["group"] == g) & (scores["split"] == split) & (scores["model"] == "LightGBM")]
                full = s[s["inputs"] == "full PACE spectrum"]["skill"].iloc[0]
                virt = s[s["inputs"] == "virtual VIIRS (5 bands)"]["skill"].iloc[0]
                rv = r[(r["group"] == g) & (r["split"] == split) & (r["model"] == "A")]["skill"].iloc[0]
                L.append(f"| {NAMES[g]} | {split} | {full:.2f} | {virt:.2f} | {rv:.2f} | {full - virt:.2f} | "
                         f"{virt - rv:.2f} |")
        L += ["", "Inputs differ slightly between steps (CMC SST here, OISST and PAR in the reconstruction) and the "
              "pixel samples are drawn separately, so the split is approximate.", ""]
    out = cfg.path("reports") / "rq2_information.md"
    out.write_text("\n".join(L))
    return out
