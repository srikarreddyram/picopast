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
    table_path = work / f"training_{sensor}_{mapping}.parquet"
    if stages == ("table",):
        return table_path

    scores_path = work / f"validation_{sensor}_{mapping}.csv"
    if "validate" in stages and (refresh or not scores_path.exists()):
        scores = pd.concat([evaluate(table, g, block_deg=cfg["grid"]["block_deg"],
                                     holdout_months=cfg["periods"]["temporal_holdout_months"]) for g in GROUPS])
        scores.to_csv(scores_path, index=False)
    if not scores_path.exists():
        raise FileNotFoundError(f"{scores_path} missing; run picopast reconstruct --stage validate")
    scores = pd.read_csv(scores_path)
    if stages == ("validate",):
        return scores_path

    finals, routes = {}, {}
    for g in GROUPS:
        route, beats = choose_route(scores[scores["group"] == g])
        routes[g] = (route, beats)
        finals[g] = fit_final(table, g, route)
    if "apply" in stages:
        reconstruct(cfg, finals, sensor, tiers, refresh=refresh)
    products = cfg.path("products") / "reconstruction"
    products.mkdir(parents=True, exist_ok=True)
    series_path = products / f"{sensor}_8D_regional_means.csv"
    if "apply" in stages or refresh or not series_path.exists():
        regional_series(cfg, sensor).to_csv(series_path, index=False)
    series = pd.read_csv(series_path, parse_dates=["time"])
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
        ends = []
        for region, (label, color) in REGIONS.items():
            d = series[series["region"] == region].set_index("time")[g]
            m = ((d - finals[g]["index_mean"]) / finals[g]["index_sd"]).resample("MS").mean()
            ax.plot(m.index, m.values, color=color, lw=1.6, label=label)
            last = m.dropna()
            if len(last):
                ends.append([last.values[-1], last.index[-1], label])
        lo, hi = ax.get_ylim()
        gap = 0.08 * (hi - lo)
        ends.sort(key=lambda e: e[0])
        for i in range(1, len(ends)):  # push labels apart so they never overlap
            ends[i][0] = max(ends[i][0], ends[i - 1][0] + gap)
        for y, x, label in ends:
            ax.annotate(label, (x, y), xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                        color=ink, annotation_clip=False)
        kind, _, rest = tiers[g].partition(":")
        status = rest.split("(")[0].strip()
        exploratory = kind.strip() == "exploratory"
        ax.set_title(f"{NAMES[g][0].upper() + NAMES[g][1:]} — {'EXPLORATORY' if exploratory else kind.strip()} "
                     f"({status})", loc="left", fontsize=10, color=ink, fontweight="bold" if exploratory else "normal")
        ax.axhline(0, color=muted, lw=0.6)
        ax.grid(axis="y", color="#e4e3dc", lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(muted)
        ax.tick_params(colors=muted, labelsize=8)
        ax.set_ylabel("relative index (SD units)", fontsize=8, color=muted)
    axes[0].text(train[0], axes[0].get_ylim()[1], " training (PACE era)", va="top", fontsize=8, color=muted)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.01, 0.965), ncol=3, frameon=False, fontsize=8,
               labelcolor=ink)
    fig.suptitle("Back-extension from VIIRS SNPP, monthly means of 8-day maps", x=0.01, ha="left", fontsize=11,
                 color=ink, y=0.995)
    fig.tight_layout(rect=(0, 0, 0.9, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=surface)
    plt.close(fig)
    return out


def _out_of_range_by_year(series: pd.DataFrame) -> pd.DataFrame:
    """Share of reconstructed study-box pixels with at least one input outside the training range."""
    box = series[series["region"] == "study_box"].assign(year=lambda d: d["time"].dt.year).groupby("year").sum(
        numeric_only=True)
    return pd.DataFrame({g: box[f"{g}_n_out"] / box[f"{g}_n"] for g in GROUPS})


def _shift_text(box_series: pd.DataFrame, finals: dict) -> str:
    """Early-to-late change of the study-box index per group: a trend or sensor drift, not yet separable."""
    if not len(box_series):
        return ""
    d = box_series.set_index("time")
    shift = {g: ((d[g] - finals[g]["index_mean"]) / finals[g]["index_sd"]) for g in GROUPS}
    parts = [f"{NAMES[g]} {z['2023':'2025'].mean() - z['2012':'2014'].mean():+.2f}" for g, z in shift.items()]
    return ("Study-box index, 2023–25 mean minus 2012–14 mean (SD units): " + ", ".join(parts) + ". Such a "
            "shift could be a real change (warming favours cyanobacteria over "
            "picoeukaryotes) or slow drift in SNPP's bands; until the NOAA-20 copy is compared, it is not "
            "reported as a trend.")


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
        oor = _out_of_range_by_year(series)
        same = all(finals[g]["route"] == finals["peuk"]["route"] for g in GROUPS)
        cols = ["peuk"] if same else list(GROUPS)
        pace = oor.loc[oor.index >= first.year, cols[0]].mean()
        L += ["", "Share of reconstructed pixels with at least one input outside the 0.5–99.5% range seen in "
              "training" + (" (all groups use the same inputs, so one column):" if same else ":"), "",
              "| year | " + " | ".join("share" if same else NAMES[g] for g in cols) + " |",
              "| --- |" + " --- |" * len(cols)]
        for y, r in oor.iterrows():
            L.append(f"| {y} | " + " | ".join(f"{r[g]:.1%}" for g in cols) + " |")
        L += ["", f"About {pace:.0%} is the floor even inside the training years: with {len(finals['peuk']['ranges'])} "
              "inputs each allowed a 1% tail, some in-range pixels are always flagged. Years well above that floor "
              "are where the inputs drift from what the models saw."]
    else:
        L.append("Not run yet: no reconstructed periods on disk.")

    screen_path = cfg.path("interim") / "reconstruct" / f"training_{sensor}_{mapping}.screen.json"
    screen = json.loads(screen_path.read_text()) if screen_path.exists() else None
    sst_share = float(finals["pro"]["importance"].get("sst", np.nan))
    sw = table[table["month"].isin(cfg["seasons"]["sw_monsoon"])]
    box_series = series[series["region"] == "study_box"] if len(series) else series
    sw_cover = (box_series.loc[box_series["time"].dt.month.isin(cfg["seasons"]["sw_monsoon"]), "valid_share"].mean()
                if len(box_series) else np.nan)
    other_cover = (box_series.loc[~box_series["time"].dt.month.isin(cfg["seasons"]["sw_monsoon"]), "valid_share"]
                   .mean() if len(box_series) else np.nan)
    L += ["", "## Bottlenecks", "",
          "All groups:", "",
          "- **The labels are MOANA estimates, not cell counts.** The models learn MOANA, errors included. No "
          "PACE-era cell counts exist in the box (trust gate, `docs/reports/trust_gate_fallback.md`).",
          f"- **Short training record.** {table['time'].nunique()} periods over about 2.6 years. Extending to 2012 "
          "assumes the VIIRS-to-MOANA relation holds in other years; the out-of-range table tracks input shift, "
          "not a change in the relation.",
          f"- **Monsoon cloud.** The per-period sample gives June–September {len(sw) / len(table):.0%} of training "
          "pixels, but only from the cloud-free parts of each scene; in the reconstruction "
          f"{sw_cover:.0%} of study-box pixels have a value in the SW monsoon against {other_cover:.0%} in other "
          "seasons. SW monsoon results are reported separately (PRD section 2.1).",
          "- **Missing inputs.** No mixed-layer depth yet (needs a Copernicus Marine account). Route B and the "
          "virtual-VIIRS step wait for the VIIRS response tables.",
          "- **Sensor drift.** SNPP over 2012–2026 is not yet checked against NOAA-20 (RQ4). " + _shift_text(
              box_series, finals), "",
          "Prochlorococcus (exploratory):", "",
          "- Check A *not supported* (12 matched station-days, one cruise); check C *inconclusive*.",
          ("- MOANA's Prochlorococcus is a linear fit made in the Atlantic. Of the PACE-era pixels where MOANA "
           f"ran, {screen['outside_nasa_range'] / screen['retrieved']:.1%} were left out of training: "
           f"{screen['negative_pro'] / screen['retrieved']:.1%} for a negative Prochlorococcus, the rest above "
           "NASA's valid maximum for some group or otherwise flagged. Positive values sit well above NASA's "
           "Atlantic output (open item 15)." if screen else
           "- MOANA's Prochlorococcus is a linear fit made in the Atlantic; it goes negative in this ocean "
           "(open item 15)."),
          f"- MOANA computes Prochlorococcus with an SST term, yet SST carries only {sst_share:.0%} of the model's "
          "gain; the VIIRS band ratios carry most of it." if sst_share < 0.3 else
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
