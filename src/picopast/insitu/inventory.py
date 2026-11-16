"""The week-1 ship-count inventory (PRD section 5.3: "count what is actually downloadable").

Open archives only; nobody is contacted. Writes CSV tables under
``data/insitu/`` and a dated summary to ``docs/reports/data_inventory.md``.
"""

from __future__ import annotations

from datetime import date

import pandas as pd

from picopast.config import Config
from picopast.insitu import bcodmo, pangaea, seabass
from picopast.insitu.schema import GROUP_COLUMNS, era, standardise


def _count_table(df: pd.DataFrame, cfg: Config, surface_m: float) -> pd.DataFrame:
    """Records per era, with surface-only and Phase-1 region counts per group."""
    df = df.assign(era=era(df["time_utc"], cfg["periods"]["pace"]["start"]))
    surface = df["depth_m"] <= surface_m
    rows = []
    for name, part in df.groupby("era"):
        row = {"era": name, "records": len(part), f"surface (<= {surface_m:g} m)": int(surface[part.index].sum())}
        for group, col in GROUP_COLUMNS.items():
            row[f"surface {group}"] = int((surface[part.index] & part[col].notna()).sum())
        for region, box in cfg.regions(phase=1).items():
            row[f"surface in {region}"] = int(
                (surface[part.index] & box.contains(part["lat"], part["lon"])).sum()
            )
        years = part["time_utc"].dt.year.dropna()
        row["years"] = f"{int(years.min())}-{int(years.max())}" if len(years) else "-"
        rows.append(row)
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df.empty:
        return "_none_\n"
    df = df.head(max_rows).astype(str)
    head = "| " + " | ".join(df.columns) + " |\n| " + " | ".join("---" for _ in df.columns) + " |\n"
    return head + "\n".join("| " + " | ".join(row) + " |" for row in df.values) + "\n"


def run(cfg: Config, skip_seabass: bool = False) -> str:
    box = cfg.box
    out_dir = cfg.path("insitu")
    out_dir.mkdir(parents=True, exist_ok=True)
    surface_m = float(cfg["matchup"]["max_depth_m"])
    pace_start = cfg["periods"]["pace"]["start"]
    md = [
        "# Ship-count inventory (open archives only)",
        "",
        f"Generated {date.today().isoformat()} by `picopast inventory`. Box: "
        f"{box.lat_min:g}-{box.lat_max:g}N, {box.lon_min:g}-{box.lon_max:g}E. "
        f"Eras: pre_viirs (before 2012), viirs_era (2012 to {pace_start}, serves RQ4), "
        f"pace_era (from {pace_start}, serves the RQ1 trust gate).",
        "",
    ]

    # BCO-DMO: records can be counted exactly because the table is downloadable.
    counts = standardise(bcodmo.fetch_global_compilation(box))
    counts.to_csv(out_dir / "bcodmo_global_compilation_box.csv", index=False)
    md += [
        "## BCO-DMO global compilation (dataset 793451)",
        "",
        f"{len(counts)} records with at least one group count, downloaded to "
        "`data/insitu/bcodmo_global_compilation_box.csv`.",
        "",
        _md(_count_table(counts, cfg, surface_m)),
        "Other BCO-DMO datasets that mention the groups and touch the box:",
        "",
        _md(bcodmo.search(box)),
    ]

    # PANGAEA: one row per dataset (usually one station profile).
    pg = pangaea.inventory(box, cfg.path("data") / "cache" / "pangaea")
    pg.to_csv(out_dir / "pangaea_datasets_box.csv", index=False)
    md += ["## PANGAEA", ""]
    if pg.empty:
        md += ["No datasets found.", ""]
    else:
        pg_era = era(pd.to_datetime(pg["time_start"], errors="coerce", utc=True).dt.tz_localize(None),
                     pace_start)
        summary = (
            pg.assign(era=pg_era)
            .groupby("era")
            .agg(datasets=("doi", "size"), with_pro=("has_pro", "sum"), with_syn=("has_syn", "sum"),
                 with_peuk=("has_peuk", "sum"), first=("time_start", "min"), last=("time_start", "max"))
            .reset_index()
        )
        summary[["first", "last"]] = summary[["first", "last"]].apply(lambda s: s.str[:10])
        md += [
            f"{len(pg)} datasets in the box report at least one group "
            "(`data/insitu/pangaea_datasets_box.csv`). A dataset is usually one station profile.",
            "",
            _md(summary),
        ]
        recent = pg[pg_era != "pre_viirs"][["doi", "title", "time_start", "lat", "lon"]]
        md += ["Datasets from 2012 onward:", "", _md(recent)]

    # SeaBASS: CMR lists files but not their fields, so this is a shortlist, not a count.
    if not skip_seabass:
        cols, files = seabass.candidates(box)
        cols.to_csv(out_dir / "seabass_collections.csv", index=False)
        files.to_csv(out_dir / "seabass_candidate_files.csv", index=False)
        regional = cols[cols["regional"]].sort_values("time_start")
        other = cols[~cols["regional"] & (cols["candidate_files"] > 0)].sort_values("time_start")
        failed = cols[cols["error"] != ""]
        md += [
            "## SeaBASS (via NASA CMR)",
            "",
            "CMR lists every SeaBASS file with its name, date and footprint, but not the fields "
            "inside it, and SeaBASS collection footprints are all global. So this section is a "
            "shortlist: **nothing here is a confirmed cell count until the file header is read, "
            "which needs the Earthdata login.**",
            "",
            f"Experiments whose title or summary names the region ({len(regional)}). Every one "
            "should be opened after login:",
            "",
            _md(regional[["collection", "title", "time_start", "time_end", "files_listed",
                          "files_in_box", "files_no_footprint", "candidate_files"]]),
            f"Other experiments with files named like flow cytometry (`*fcm*`, `*cyto*`) whose "
            f"footprint lies in the box ({len(other)}):",
            "",
            _md(other[["collection", "title", "time_start", "time_end", "candidate_files"]]),
            f"{len(files)} candidate files in total (`data/insitu/seabass_candidate_files.csv`).",
            "",
        ]
        if not files.empty:
            by_col = (
                files.groupby("collection")
                .agg(files=("file", "size"), first=("time_start", "min"), last=("time_start", "max"),
                     example=("file", "first"))
                .reset_index()
                .sort_values("last", ascending=False)
            )
            md += [_md(by_col)]
        if len(failed):
            md += [f"{len(failed)} collections could not be listed (see the `error` column).", ""]

    report = cfg.path("reports") / "data_inventory.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(md))
    return str(report)
