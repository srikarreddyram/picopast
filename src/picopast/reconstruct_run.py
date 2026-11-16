"""Run the back-extension end to end and write ``docs/reports/reconstruction.md``.

Stages: training table -> validation (spatial blocks, temporal hold-out) ->
final models -> reconstruction of every VIIRS 8-day period -> regional series
-> report. Each stage caches its output under ``data/interim/reconstruct`` or
``data/products``, so a rerun only does what is new.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from picopast.config import Config
from picopast.reconstruct import (
    FEATURES,
    GROUPS,
    UNITS,
    evaluate,
    fit_final,
    reconstruct,
    regional_series,
    training_table,
    viirs_periods,
)

NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes"}


def _season(cfg: Config, month: int) -> str:
    return next(s for s, months in cfg["seasons"].items() if month in months)


def choose_route(scores: pd.DataFrame) -> tuple[str, bool]:
    """Route A or C, whichever has the higher pooled spatial-block skill; and whether it beats climatology."""
    s = scores[(scores["split"] == "spatial") & scores["model"].isin(["A", "C"])].set_index("model")["skill"]
    route = str(s.idxmax())
    return route, bool(s[route] > 0)


def run(cfg: Config, stages=("table", "validate", "apply", "report"), refresh: bool = False) -> Path:
    rc = cfg["reconstruction"]
    sensor, mapping, tiers = rc["sensor"], rc["mapping"], rc["tiers"]
    work = cfg.path("interim") / "reconstruct"
    work.mkdir(parents=True, exist_ok=True)

    table = training_table(cfg, sensor, mapping, rc["per_period_sample"], refresh=refresh and "table" in stages)

    scores_path = work / f"validation_{sensor}_{mapping}.csv"
    if "validate" in stages and (refresh or not scores_path.exists()):
        scores = pd.concat([evaluate(table, g, block_deg=cfg["grid"]["block_deg"],
                                     holdout_months=cfg["periods"]["temporal_holdout_months"]) for g in GROUPS])
        scores.to_csv(scores_path, index=False)
    scores = pd.read_csv(scores_path)

    finals, routes = {}, {}
    for g in GROUPS:
        route, beats = choose_route(scores[scores["group"] == g])
        routes[g] = (route, beats)
        finals[g] = fit_final(table, g, route)
    if "apply" in stages:
        reconstruct(cfg, finals, sensor, tiers, refresh=refresh)
    series = regional_series(cfg, sensor)
    series.to_csv(cfg.path("products") / "reconstruction" / f"{sensor}_8D_regional_means.csv", index=False)
    meta = {g: {"route": routes[g][0], "beats_climatology": routes[g][1], "index_mean": finals[g]["index_mean"],
                "index_sd": finals[g]["index_sd"], "importance": finals[g]["importance"].round(4).to_dict()}
            for g in GROUPS}
    (work / f"final_{sensor}_{mapping}.json").write_text(json.dumps(meta, indent=2))
    return write_report(cfg, table, scores, finals, routes, series, sensor, mapping, tiers)


REGIONS = {"study_box": ("Study box", "#2a78d6"), "bay_of_bengal": ("Bay of Bengal", "#eb6834"),
           "arabian_sea_central_east": ("Arabian Sea C/E", "#1baf7a")}  # validated categorical slots 1-3


def plot_series(series: pd.DataFrame, finals: dict, tiers: dict, train: tuple, out: Path) -> Path:
    """Monthly regional means as a relative index, one panel per group; tier in each title."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, surface = "#2b2b29", "#6b6a63", "#fcfcfb"
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True, facecolor=surface)
    for ax, g in zip(axes, ("peuk", "syn", "pro")):
        ax.set_facecolor(surface)
        ax.axvspan(*train, color="#e4e3dc", lw=0)
        for region, (label, color) in REGIONS.items():
            d = series[series["region"] == region].set_index("time")[g]
            m = ((d - finals[g]["index_mean"]) / finals[g]["index_sd"]).resample("MS").mean()
            ax.plot(m.index, m.values, color=color, lw=1.6, label=label)
            last = m.dropna()
            if len(last):
                ax.annotate(label, (last.index[-1], last.values[-1]), xytext=(6, 0), textcoords="offset points",
                            va="center", fontsize=8, color=ink)
        tier = tiers[g]
        flag = "" if tier.startswith("relative index") else "EXPLORATORY — "
        ax.set_title(f"{flag}{NAMES[g]}: {tier}", loc="left", fontsize=10, color=ink,
                     fontweight="bold" if flag else "normal")
        ax.axhline(0, color=muted, lw=0.6)
        ax.grid(axis="y", color="#e4e3dc", lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(muted)
        ax.tick_params(colors=muted, labelsize=8)
        ax.set_ylabel("relative index (SD units)", fontsize=8, color=muted)
    axes[0].text(train[0], axes[0].get_ylim()[1], " training (PACE era)", va="top", fontsize=8, color=muted)
    axes[0].legend(loc="upper left", bbox_to_anchor=(0, 1.32), ncol=3, frameon=False, fontsize=8, labelcolor=ink)
    fig.suptitle("Back-extension from VIIRS SNPP, monthly means of 8-day maps", x=0.01, ha="left", fontsize=11,
                 color=ink, y=0.995)
    fig.tight_layout(rect=(0, 0, 0.93, 0.97))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=surface)
    plt.close(fig)
    return out


def _out_of_range_by_year(cfg: Config, sensor: str) -> pd.DataFrame:
    import xarray as xr

    rows = []
    for path in sorted((cfg.path("products") / "reconstruction" / sensor / "8D").glob("recon_*.nc")):
        with xr.open_dataset(path) as ds:
            row = {"year": pd.Timestamp(ds.attrs["period_start"]).year}
            for g in GROUPS:
                valid = np.isfinite(ds[g].values)
                n_out = ds[f"{g}_out_of_range"].values
                row[g] = (n_out[valid] > 0).sum()
                row[f"{g}_n"] = valid.sum()
            rows.append(row)
    d = pd.DataFrame(rows).groupby("year").sum()
    return pd.DataFrame({g: d[g] / d[f"{g}_n"] for g in GROUPS})


def write_report(cfg, table, scores, finals, routes, series, sensor, mapping, tiers) -> Path:
    periods = viirs_periods(cfg, sensor)
    first, last = table["time"].min(), table["time"].max()
    L = ["# Back-extension to 2012: first attempt", "",
         f"Generated {date.today().isoformat()} by `picopast reconstruct`. Models learn MOANA ({mapping} mapping) "
         f"from VIIRS {sensor.upper()} 8-day inputs in the PACE era and are applied to every VIIRS 8-day period on "
         "disk. They reproduce MOANA, so they can be no more right than MOANA is.", "",
         "## What each group's record is", "",
         "| group | tier | route used | beats climatology on held-out blocks |", "| --- | --- | --- | --- |"]
    for g in GROUPS:
        route, beats = routes[g]
        L.append(f"| {NAMES[g]} | {tiers[g]} | {route} | {'yes' if beats else '**no**'} |")
    L += ["", "Only a group that passed the trust gate is used for headline results. The exploratory records are "
          "an attempt, shown with their bottlenecks below; they are not validated estimates.", "",
          "## Training data", "",
          f"- PACE-era 8-day periods with both MOANA and VIIRS: {table['time'].nunique()} "
          f"({first:%d %b %Y} to {last:%d %b %Y}); up to {cfg['reconstruction']['per_period_sample']:,} offshore "
          "pixels sampled per period.",
          "- Pixels per group: " + ", ".join(f"{NAMES[g]} {int(np.isfinite(table[g]).sum()):,}" for g in GROUPS) + ".",
          f"- Units: Prochlorococcus {UNITS['pro']}; the others {UNITS['syn']}.",
          f"- Route A inputs: {', '.join(FEATURES['A'])}. Route C inputs: {', '.join(FEATURES['C'])}.", "",
          "## Validation (out of sample)", "",
          "Skill = 1 − RMSE_model / RMSE_climatology (> 0 beats the seasonal climatology). Spatial: 2° blocks in 5 "
          "folds, predictions pooled. Temporal: the last 6 months of the PACE era held out. Coverage: share of "
          "held-out values inside Route A's 5–95% interval (target 85–95%).", "",
          "| group | split | climatology RMSE | ridge skill | Route C skill | Route A skill | Route A ρ | Route A 90% coverage |",
          "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        for split in ("spatial", "temporal"):
            s = scores[(scores["group"] == g) & (scores["split"] == split)].set_index("model")
            L.append(f"| {NAMES[g]} | {split} | {s.loc['clim', 'rmse']:.3f} | {s.loc['ridge', 'skill']:.2f} | "
                     f"{s.loc['C', 'skill']:.2f} | {s.loc['A', 'skill']:.2f} | {s.loc['A', 'spearman']:.2f} | "
                     f"{s.loc['A', 'coverage_90']:.0%} |")
    L += ["", "Inputs that matter most (share of Route A gain, or Route C if used):", ""]
    for g in GROUPS:
        imp = finals[g]["importance"].head(5)
        L.append(f"- {NAMES[g]}: " + ", ".join(f"{k} {v:.0%}" for k, v in imp.items()))

    L += ["", "## The reconstruction", ""]
    if len(series):
        yrs = series["time"].dt.year
        box = series[series["region"] == "study_box"].assign(season=lambda d: d["time"].dt.month.map(
            lambda m: _season(cfg, m)))
        L += [f"- {len(periods)} VIIRS 8-day periods on disk, {yrs.min()}–{yrs.max()}; maps in "
              f"`data/products/reconstruction/{sensor}/8D/` (median, 5% and 95% per group, plus an out-of-range count).",
              "- Share of study-box pixels with a value, by season: " + ", ".join(
                  f"{s} {v:.0%}" for s, v in box.groupby("season")["valid_share"].mean().items()) + "."]
        fig = plot_series(series, finals, tiers, (first, last + pd.Timedelta(days=7)),
                          cfg.path("reports") / "figures" / "reconstruction_series.png")
        L += ["", f"![Regional relative index, 2012–2026](figures/{fig.name})", "",
              "Each line is a regional monthly mean of the 8-day medians, as a relative index "
              "((value − PACE-era mean) / PACE-era SD). The shaded span is the training period. The table view is "
              f"`data/products/reconstruction/{sensor}_8D_regional_means.csv`. Cloud makes June–September means "
              "rest on few pixels."]
        oor = _out_of_range_by_year(cfg, sensor)
        L += ["", "Share of reconstructed pixels with at least one input outside the PACE-era training range:", "",
              "| year | " + " | ".join(NAMES[g] for g in GROUPS) + " |", "| --- |" + " --- |" * len(GROUPS)]
        for y, r in oor.iterrows():
            L.append(f"| {y} | " + " | ".join(f"{r[g]:.1%}" for g in GROUPS) + " |")
    else:
        L.append("Not run yet: no reconstructed periods on disk.")

    neg = float((table["pro"] < 0).mean())
    sst_share = float(finals["pro"]["importance"].get("sst", np.nan))
    sw = table[table["month"].isin(cfg["seasons"]["sw_monsoon"])]
    L += ["", "## Bottlenecks", "",
          "All groups:", "",
          "- **The labels are MOANA estimates, not cell counts.** The models learn MOANA, errors included. No "
          "PACE-era cell counts exist in the box (trust gate, `docs/reports/trust_gate_fallback.md`).",
          f"- **Short training record.** {table['time'].nunique()} periods over about 2.6 years. Extending to 2012 "
          "assumes the VIIRS-to-MOANA relation holds in other years; the out-of-range table tracks input shift, "
          "not a change in the relation.",
          f"- **Monsoon cloud.** June–September periods give {len(sw) / len(table):.0%} of training pixels; the "
          "SW monsoon is thinly trained and thinly reconstructed (reported separately, PRD section 2.1).",
          "- **Missing inputs.** No mixed-layer depth yet (needs a Copernicus Marine account). Route B and the "
          "virtual-VIIRS step wait for the VIIRS response tables.",
          "- **Sensor drift.** SNPP over 2012–2026 is not yet checked against NOAA-20 (RQ4).", "",
          "Prochlorococcus (exploratory):", "",
          "- Check A *not supported* (12 matched station-days, one cruise); check C *inconclusive*.",
          f"- MOANA's Prochlorococcus is a linear fit fitted in the Atlantic: {neg:.1%} of training targets are "
          "negative, and positive values sit well above NASA's Atlantic output (open item 15).",
          f"- MOANA computes Prochlorococcus with an SST term; SST carries {sst_share:.0%} of the model's gain, so "
          "the record leans on OISST.", "",
          "Synechococcus (exploratory):", "",
          "- Check C *inconclusive*; check B found no relation with flow-through phycoerythrin.",
          "- How much Synechococcus information survives in the broad VIIRS bands (RQ2) is untested until virtual "
          "VIIRS is built.", "",
          "Picoeukaryotes (relative index):", "",
          "- Passed check C against 1991–2005 counts, a pattern check only: the record is a relative index, "
          "never cells/mL. MOANA picoeukaryotes follow chlorophyll; a chlorophyll-only comparison is still due.", ""]
    out = cfg.path("reports") / "reconstruction.md"
    out.write_text("\n".join(L) + "\n")
    return out
