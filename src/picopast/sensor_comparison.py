"""RQ4: do the SNPP and NOAA-20 back-extensions agree? (PRD section 7)

The two records are built independently (own bands, own models, same MOANA
targets). Where both have a value on the same 8-day period and pixel:

- agreement: the NOAA-20 median lies inside SNPP's 5-95% interval (SNPP is the
  main record and its interval is its stated uncertainty); the looser "the two
  intervals overlap" is reported too. PRD pass condition: >= 80% of open-ocean
  pixels agree.
- drift: the study-box mean of SNPP minus NOAA-20 per period, in units of the
  SNPP index's SD; its Theil-Sen slope per decade. A record that drifts against
  the other shows up as a slope; a real ocean change appears in both and not in
  the difference.

One period is read at a time (memory).
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from picopast.config import Config
from picopast.reconstruct import GROUPS

NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes"}


def compare(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-period agreement counts and study-box mean differences, per group."""
    base = cfg.path("products") / "reconstruction"
    snpp = {p.name: p for p in (base / "snpp" / "8D").glob("recon_*.nc")}
    noaa = {p.name: p for p in (base / "noaa20" / "8D").glob("recon_*.nc")}
    rows = []
    for name in sorted(set(snpp) & set(noaa)):
        with xr.open_dataset(snpp[name]) as a, xr.open_dataset(noaa[name]) as b:
            row = {"time": pd.Timestamp(a.attrs["period_start"])}
            for g in GROUPS:
                sa, lo, hi = (a[g + k].values for k in ("", "_lo", "_hi"))
                sb, blo, bhi = (b[g + k].values for k in ("", "_lo", "_hi"))
                both = np.isfinite(sa) & np.isfinite(sb)
                row[f"{g}_n"] = int(both.sum())
                row[f"{g}_in_snpp_interval"] = int((both & (sb >= lo) & (sb <= hi)).sum())
                row[f"{g}_overlap"] = int((both & (blo <= hi) & (bhi >= lo)).sum())
                row[f"{g}_diff"] = float(np.mean(sa[both] - sb[both])) if both.any() else np.nan
            rows.append(row)
    per_period = pd.DataFrame(rows)
    out = []
    for g in GROUPS:
        n = per_period[f"{g}_n"].sum()
        out.append({"group": g, "pixels": int(n), "periods": int((per_period[f"{g}_n"] > 0).sum()),
                    "in_snpp_interval": per_period[f"{g}_in_snpp_interval"].sum() / n,
                    "intervals_overlap": per_period[f"{g}_overlap"].sum() / n})
    return per_period, pd.DataFrame(out)


def _drift(per_period: pd.DataFrame, g: str, sd: float) -> dict:
    from scipy.stats import theilslopes

    d = per_period.set_index("time")[f"{g}_diff"].dropna() / sd
    monthly = d.resample("MS").mean().dropna()
    years = (monthly.index - monthly.index[0]).days / 365.25
    slope, _, lo, hi = theilslopes(monthly.values, years)
    return {"mean": float(d.mean()), "slope": 10 * slope, "slope_lo": 10 * lo, "slope_hi": 10 * hi,
            "monthly": monthly}


def plot(per_period: pd.DataFrame, sds: dict, out: Path, train: tuple | None = None) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, surface = "#2b2b29", "#6b6a63", "#fcfcfb"
    fig, axes = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True, facecolor=surface)
    for ax, g in zip(axes, ("peuk", "syn", "pro")):
        m = _drift(per_period, g, sds[g])["monthly"]
        ax.set_facecolor(surface)
        if train:
            ax.axvspan(*train, color="#e4e3dc", lw=0)
        ax.plot(m.index, m.values, color="#2a78d6", lw=1.6)
        ax.axhline(0, color=muted, lw=0.6)
        ax.set_title(f"{NAMES[g][0].upper() + NAMES[g][1:]}: SNPP minus NOAA-20, study-box mean", loc="left",
                     fontsize=10, color=ink)
        ax.grid(axis="y", color="#e4e3dc", lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(muted)
        ax.tick_params(colors=muted, labelsize=8)
        ax.set_ylabel("SD units", fontsize=8, color=muted)
    if train:
        axes[0].text(train[0], axes[0].get_ylim()[1], " training (PACE era)", va="top", fontsize=8, color=muted)
    fig.suptitle("Difference between the two independent back-extensions, monthly", x=0.01, ha="left",
                 fontsize=11, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=surface)
    plt.close(fig)
    return out


def report(cfg: Config) -> Path:
    per_period, summary = compare(cfg)
    work = cfg.path("interim") / "reconstruct"
    meta = {s: json.loads((work / f"final_{s}_operational.json").read_text()) for s in ("snpp", "noaa20")}
    sds = {g: meta["snpp"][g]["index_sd"] for g in GROUPS}
    tiers = cfg["reconstruction"]["tiers"]
    t0, t1 = per_period["time"].min(), per_period["time"].max()
    pace = cfg["periods"]["pace"]
    fig = plot(per_period, sds, cfg.path("reports") / "figures" / "sensor_comparison.png",
               (pd.Timestamp(pace["start"]), pd.Timestamp(pace["end"])))
    L = ["# SNPP vs NOAA-20 back-extensions (RQ4)", "",
         f"Generated {date.today().isoformat()} by `picopast sensor-comparison`. The two records are built "
         "independently (each sensor's own bands and models, trained on the same PACE-era MOANA). Compared on "
         f"{len(per_period)} common 8-day periods, {t0:%b %Y} to {t1:%b %Y}, pixel by pixel where both have a value.",
         "", "## Agreement", "",
         "PRD pass condition: at least 80% of open-ocean pixels agree within the stated uncertainty. Here: the "
         "NOAA-20 median inside SNPP's 5–95% interval. The looser test (the two intervals overlap) is shown too.", "",
         "| group | tier | pixel-periods | NOAA-20 inside SNPP interval | intervals overlap | PRD (≥80%) | "
         "typical gap between the records |",
         "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in summary.itertuples():
        gap = np.nanmedian(np.abs(per_period[f"{r.group}_diff"])) / sds[r.group]
        L.append(f"| {NAMES[r.group]} | {tiers[r.group].split(':')[0]} | {r.pixels:,} | {r.in_snpp_interval:.0%} | "
                 f"{r.intervals_overlap:.0%} | {'pass' if r.in_snpp_interval >= 0.8 else '**fail**'} | "
                 f"{gap:.2f} SD |")
    L += ["", "The test is easy to pass: each interval carries the model's error against MOANA, which is much "
          "larger than the gap between two models fed by similar sensors. The last column (median over periods of "
          "|SNPP − NOAA-20| in the study-box mean, in SD units of the index) is the sharper measure."]
    L += ["", "## Drift between the records", "",
          "Study-box mean of SNPP minus NOAA-20 per period, in SD units of each group's index. A slope means one "
          "record drifts against the other; a change in the ocean would appear in both records and cancel here.",
          "", "| group | mean difference | Theil–Sen slope per decade (95%) |", "| --- | --- | --- |"]
    for g in GROUPS:
        d = _drift(per_period, g, sds[g])
        L.append(f"| {NAMES[g]} | {d['mean']:+.2f} | {d['slope']:+.2f} [{d['slope_lo']:+.2f}, {d['slope_hi']:+.2f}] |")
    ser = {k: pd.read_csv(cfg.path("products") / "reconstruction" / f"{k}_8D_regional_means.csv", parse_dates=["time"])
           for k in ("snpp", "noaa20")}
    L += ["", "## Change since 2018 in each record", "",
          "Study-box index, 2023–2025 mean minus 2018–2020 mean, each record in its own SD units. A change seen by "
          "both sensors is a candidate ocean signal (still to be tested in RQ5); a change seen by one only is "
          "treated as that sensor's drift.", "",
          "| group | SNPP | NOAA-20 | seen by both |", "| --- | --- | --- | --- |"]
    for g in GROUPS:
        ch = {}
        for k in ("snpp", "noaa20"):
            d = ser[k][ser[k]["region"] == "study_box"].set_index("time")[g]
            z = (d - meta[k][g]["index_mean"]) / meta[k][g]["index_sd"]
            ch[k] = z["2023":"2025"].mean() - z["2018":"2020"].mean()
        both = np.sign(ch["snpp"]) == np.sign(ch["noaa20"]) and min(abs(ch["snpp"]), abs(ch["noaa20"])) >= 0.05
        L.append(f"| {NAMES[g]} | {ch['snpp']:+.2f} | {ch['noaa20']:+.2f} | {'yes' if both else 'no'} |")
    L += ["", "\"Seen by both\" means the same sign and at least 0.05 SD in each record.", "",
          f"![SNPP minus NOAA-20]({'figures/' + fig.name})", "",
          "The early-to-late shift in the SNPP record (`docs/reports/reconstruction.md`, Sensor drift) is read "
          "against these tables: a slope near zero and a change seen by both mean SNPP drift does not explain it.",
          ""]
    per_period.to_csv(cfg.path("products") / "reconstruction" / "snpp_vs_noaa20_8D.csv", index=False)
    out = cfg.path("reports") / "sensor_comparison.md"
    out.write_text("\n".join(L))
    return out
