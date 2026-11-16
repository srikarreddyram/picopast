"""Download NASA products cut to the study box (PRD section 6 step 1: "subset at download").

Needs a free Earthdata login. Run ``picopast login`` once; it stores the
credentials in ``~/.netrc`` and nothing in this repository.

A global PACE L3M AOP daily file is about 0.8 GB, so we never keep whole files:
each granule is opened over HTTPS, only the box is read, and a small NetCDF
subset is written to ``data/raw/<product>/``. Small files (VIIRS single-band
maps) are downloaded whole to a temp file first, which is faster than many
range requests.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import xarray as xr

from picopast.config import Box, Config
from picopast.grid import subset_box
from picopast.ingest import cmr

REMOTE_READ_ABOVE_MB = 100.0


def login(interactive: bool = True):
    """Log in to Earthdata. Uses ~/.netrc if present, otherwise prompts and saves it."""
    import earthaccess

    try:
        auth = earthaccess.login(strategy="netrc")
        if auth.authenticated:
            return auth
    except Exception:
        pass
    if not interactive:
        raise RuntimeError("No Earthdata credentials found. Run `picopast login` first.")
    return earthaccess.login(strategy="interactive", persist=True)


def is_logged_in() -> bool:
    try:
        login(interactive=False)
        return True
    except Exception:
        return False


def find(cfg: Config, product: str, start: str, end: str, period: str = "DAY") -> list[cmr.Granule]:
    """Granules of a configured product. ``period`` is DAY or 8D."""
    spec = cfg["products"][product]
    pattern = spec.get("name_pattern", "*").format(period=period)
    return cmr.granules(
        spec["short_name"],
        f"{start}T00:00:00Z",
        f"{end}T23:59:59Z",
        name_pattern=None if pattern == "*" else pattern,
    )


def _write_subset(ds: xr.Dataset, region: Box, variables: list[str] | None, out: Path) -> None:
    if variables:
        ds = ds[[v for v in variables if v in ds]]
    sub = subset_box(ds, region).load()
    encoding = {v: {"zlib": True, "complevel": 4} for v in sub.data_vars}
    tmp = out.with_suffix(".part")
    sub.to_netcdf(tmp, encoding=encoding)
    tmp.rename(out)


def fetch(
    cfg: Config,
    product: str,
    start: str,
    end: str,
    period: str = "DAY",
    box: bool = True,
    limit: int | None = None,
    region: Box | None = None,
    tag: str = "box",
) -> list[Path]:
    """Fetch a product for a date range as box subsets. Existing files are skipped.

    ``region`` replaces the study box (used for the Atlantic MOANA code check);
    give it a ``tag`` so its files do not collide with the study-box subsets.
    """
    import earthaccess

    login(interactive=False)
    region = region or cfg.box
    spec = cfg["products"][product]
    out_dir = cfg.path("raw") / product
    out_dir.mkdir(parents=True, exist_ok=True)
    granules = [g for g in find(cfg, product, start, end, period) if g.url]
    if limit:
        granules = granules[:limit]
    fs = earthaccess.get_fsspec_https_session()
    written = []
    for g in granules:
        stem = g.name[:-3] if g.name.endswith(".nc") else g.name
        out = out_dir / (f"{stem}.{tag}.nc" if box else f"{stem}.nc")
        if out.exists():
            written.append(out)
            continue
        if not box:
            fs.get(g.url, str(out))
        elif g.size_mb > REMOTE_READ_ABOVE_MB:
            with fs.open(g.url, block_size=4 * 2**20) as f, xr.open_dataset(
                f, engine="h5netcdf"
            ) as ds:
                _write_subset(ds, region, spec.get("variables"), out)
        else:
            with tempfile.TemporaryDirectory() as tmp:
                local = Path(tmp) / "granule.nc"
                fs.get(g.url, str(local))
                with xr.open_dataset(local) as ds:
                    _write_subset(ds, region, spec.get("variables"), out)
        written.append(out)
    return written
