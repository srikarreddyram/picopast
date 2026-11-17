"""Command line for the pipeline: ``picopast <command>``. Run ``picopast --help``."""

from __future__ import annotations

import argparse
import json
from datetime import date

from picopast.config import load_config


def _login(args):
    from picopast.ingest import earthdata

    earthdata.login(interactive=True)
    print("Earthdata login works. Credentials are stored in ~/.netrc.")


def _check_access(args):
    """What exists for one sample month, per product, and whether we can download it."""
    from picopast.ingest import cmr, earthdata

    cfg = load_config(args.config)
    start, end = f"{args.month}-01", f"{args.month}-28"
    lines = [
        "# Data access check",
        "",
        f"Generated {date.today().isoformat()} by `picopast check-access`. Sample window: "
        f"{start} to {end}. Counts come from NASA CMR (no login needed to search).",
        "",
        "| product | collection | version | daily files | 8-day files | example | MB each |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for key, spec in cfg["products"].items():
        info = cmr.collections(short_name=spec["short_name"])
        version = info[0].get("version_id", "?") if info else "NOT FOUND"
        if spec.get("name_pattern", "*") == "*":
            daily, eight = earthdata.find(cfg, key, start, end), []
        else:
            daily = earthdata.find(cfg, key, start, end, "DAY")
            eight = earthdata.find(cfg, key, start, end, "8D")
        example = daily[0] if daily else None
        lines.append(
            f"| {key} | {spec['short_name']} | {version} | {len(daily)} | {len(eight)} | "
            f"{example.name if example else '-'} | {example.size_mb:.0f} |"
            if example
            else f"| {key} | {spec['short_name']} | {version} | 0 | {len(eight)} | - | - |"
        )
    logged_in = earthdata.is_logged_in()
    lines += ["", f"Earthdata login available: **{'yes' if logged_in else 'no'}**"
              + ("" if logged_in else " (run `picopast login`; downloads are blocked until then)"), ""]
    out = cfg.path("reports") / "data_access.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines))
    print("\n".join(lines))


def _inventory(args):
    from picopast.insitu import inventory

    print("wrote", inventory.run(load_config(args.config), skip_seabass=args.skip_seabass))


def _fetch_open(args):
    from picopast.ingest import open_data

    cfg = load_config(args.config)
    print("indices:", {k: str(v) for k, v in open_data.fetch_indices(cfg).items()})
    for year in range(args.start_year, args.end_year + 1):
        print("oisst:", open_data.fetch_oisst_year(cfg, year, method=args.oisst_method))


def _fetch(args):
    from picopast.ingest import earthdata

    cfg = load_config(args.config)
    paths = earthdata.fetch(cfg, args.product, args.start, args.end, args.period, limit=args.limit)
    print(f"{len(paths)} files in {paths[0].parent if paths else cfg.path('raw') / args.product}")


def _moana_check(args):
    from picopast.moana import atlantic_check

    report = atlantic_check.run(load_config(args.config), args.date, pc_mapping=args.pc_mapping)
    print(json.dumps(report, indent=2))


def _moana_run(args):
    from picopast.moana import box_run

    cfg = load_config(args.config)
    written, skipped = box_run.run(cfg, args.start, args.end, args.period, args.pc_mapping, args.refresh)
    print(f"{len(written)} MOANA maps in {cfg.path('interim') / 'moana' / args.pc_mapping / args.period}")
    if skipped:
        print(f"{len(skipped)} scenes skipped (no CMC SST yet): {skipped[0]} ...")


def _fallback_checks(args):
    from picopast import fallback_run

    out = fallback_run.report(load_config(args.config))
    print(out.read_text())


def _reconstruct(args):
    from picopast import reconstruct_run

    stages = ("table", "validate", "apply", "report") if args.stage == "all" else (args.stage,)
    out = reconstruct_run.run(load_config(args.config), stages=stages, refresh=args.refresh, sensor=args.sensor)
    print(out.read_text() if out.suffix == ".md" else f"wrote {out}")


def _chl_comparison(args):
    from picopast import fallback_run

    print(fallback_run.chlorophyll_comparison(load_config(args.config)).read_text())


def _sensor_comparison(args):
    from picopast import sensor_comparison

    print(sensor_comparison.report(load_config(args.config)).read_text())


def _climate(args):
    from picopast import climate

    print(climate.report(load_config(args.config)).read_text())


def _virtual_viirs_check(args):
    from picopast import virtual_viirs_run

    print(virtual_viirs_run.report(load_config(args.config)).read_text())


def _rq2(args):
    from picopast import rq2

    print(rq2.report(load_config(args.config)).read_text())


def main(argv=None):
    p = argparse.ArgumentParser(prog="picopast", description=__doc__)
    p.add_argument("--config", default=None, help="path to picopast.yaml (default: config/)")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("login", help="store an Earthdata login in ~/.netrc").set_defaults(func=_login)

    c = sub.add_parser("check-access", help="list what NASA holds for a sample month")
    c.add_argument("--month", default="2025-01", help="YYYY-MM")
    c.set_defaults(func=_check_access)

    c = sub.add_parser("inventory", help="count ship measurements in open archives")
    c.add_argument("--skip-seabass", action="store_true", help="skip the slow SeaBASS listing")
    c.set_defaults(func=_inventory)

    c = sub.add_parser("fetch-open", help="download DMI, ONI and OISST (no login)")
    c.add_argument("--start-year", type=int, default=1991)
    c.add_argument("--end-year", type=int, default=date.today().year)
    c.add_argument("--oisst-method", default="opendap", choices=["opendap", "file"],
                   help="opendap: box only, rate-limited; file: 0.5 GB per year, then cut")
    c.set_defaults(func=_fetch_open)

    c = sub.add_parser("fetch", help="download a NASA product cut to the study box (needs login)")
    c.add_argument("product")
    c.add_argument("start", help="YYYY-MM-DD")
    c.add_argument("end", help="YYYY-MM-DD")
    c.add_argument("--period", default="DAY", choices=["DAY", "8D"])
    c.add_argument("--limit", type=int, default=None, help="at most this many files (for testing)")
    c.set_defaults(func=_fetch)

    c = sub.add_parser("moana-check", help="compare our MOANA with NASA's Atlantic L4 file")
    c.add_argument("date", help="YYYY-MM-DD")
    c.add_argument("--pc-mapping", default="operational", choices=["operational", "atbd"])
    c.set_defaults(func=_moana_check)

    c = sub.add_parser("moana-run", help="run MOANA on downloaded PACE box files")
    c.add_argument("start", help="YYYY-MM-DD")
    c.add_argument("end", help="YYYY-MM-DD")
    c.add_argument("--period", default="DAY", choices=["DAY", "8D"])
    c.add_argument("--pc-mapping", default="operational", choices=["operational", "atbd"])
    c.add_argument("--refresh", action="store_true", help="recompute existing maps")
    c.set_defaults(func=_moana_run)

    sub.add_parser(
        "fallback-checks", help="pre-registered fallback checks A, B and C"
    ).set_defaults(func=_fallback_checks)

    sub.add_parser(
        "chl-comparison", help="check C with a chlorophyll climatology in place of MOANA (descriptive)"
    ).set_defaults(func=_chl_comparison)

    sub.add_parser(
        "sensor-comparison", help="RQ4: SNPP vs NOAA-20 back-extensions (run reconstruct for both first)"
    ).set_defaults(func=_sensor_comparison)

    sub.add_parser(
        "virtual-viirs-check", help="PACE spectra through VIIRS band responses vs real same-day VIIRS"
    ).set_defaults(func=_virtual_viirs_check)

    sub.add_parser("rq2", help="RQ2: full PACE spectra vs virtual VIIRS bands, same models").set_defaults(func=_rq2)

    sub.add_parser("climate", help="RQ5: trends, IOD/ENSO lags, heatwave composites").set_defaults(func=_climate)

    c = sub.add_parser("reconstruct", help="back-extend MOANA to 2012 from VIIRS (all groups, tiered)")
    c.add_argument("--stage", default="all", choices=["all", "table", "validate", "apply", "report"])
    c.add_argument("--sensor", default=None, choices=["snpp", "noaa20"], help="default: config reconstruction.sensor")
    c.add_argument("--refresh", action="store_true", help="redo cached stages")
    c.set_defaults(func=_reconstruct)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
