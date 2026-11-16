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
import zlib
from pathlib import Path

import numpy as np
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


_PACKING_KEYS = ("dtype", "scale_factor", "add_offset", "_FillValue")


def _write_netcdf(sub: xr.Dataset, out: Path) -> None:
    """Write with compression while keeping NASA's packing (int16 + scale_factor).

    ``to_netcdf(encoding=...)`` replaces a variable's encoding, so the packing
    keys are carried over explicitly; without them Rrs is stored as float32 at
    about twice the size, with no gain in precision.
    """
    encoding = {}
    for v in sub.data_vars:
        enc = {k: sub[v].encoding[k] for k in _PACKING_KEYS if k in sub[v].encoding}
        encoding[v] = {**enc, "zlib": True, "complevel": 4}
    tmp = out.with_suffix(".part")
    sub.to_netcdf(tmp, encoding=encoding)
    tmp.rename(out)


def _write_subset(ds: xr.Dataset, region: Box, variables: list[str] | None, out: Path) -> None:
    if variables:
        ds = ds[[v for v in variables if v in ds]]
    _write_netcdf(subset_box(ds, region).load(), out)


def _attr(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray) and v.size == 1:
        return v.reshape(())[()]  # numpy scalar: keeps float32 scale_factor as float32
    return v


def _unshuffle(raw: bytes, dtype: np.dtype) -> np.ndarray:
    """Undo HDF5's byte-shuffle filter: byte k of every element is stored together."""
    b = np.frombuffer(raw, dtype=np.uint8)
    return b.reshape(dtype.itemsize, -1).T.copy().view(dtype).ravel()


def _merge_ranges(spans: list[tuple[int, int]], max_gap: int) -> list[tuple[int, int]]:
    """Merge sorted (start, end) byte spans whose gap is at most ``max_gap``."""
    merged: list[list[int]] = []
    for s, e in sorted(spans):
        if merged and s - merged[-1][1] <= max_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def _u(b: bytes, pos: int, n: int) -> int:
    return int.from_bytes(b[pos: pos + n], "little")


def _header_messages(b: bytes, version: int, first: bool, flags: int = 0):
    """(type, body) of each message in one HDF5 object-header block (v1 or v2 header).

    ``flags`` is the v2 header's flag byte (from the first block); continuation
    blocks need it to know whether messages carry a creation-order field.
    """
    if version == 2:
        if first:
            pos = 6 + (16 if flags & 0x20 else 0) + (4 if flags & 0x10 else 0)
            size_len = 1 << (flags & 0x03)
            end = pos + size_len + _u(b, pos, size_len)
            pos += size_len
        else:  # "OCHK" continuation block, checksum at the end
            pos, end = 4, len(b) - 4
        extra = 2 if flags & 0x04 else 0
        while pos + 4 + extra <= end:
            mtype, size = b[pos], _u(b, pos + 1, 2)
            pos += 4 + extra
            yield mtype, b[pos: pos + size]
            pos += size
    else:
        pos, end = (16, 16 + _u(b, 8, 4)) if first else (0, len(b))
        while pos + 8 <= end:
            mtype, size = _u(b, pos, 2), _u(b, pos + 2, 2)
            pos += 8
            yield mtype, b[pos: pos + size]
            pos += size


def _chunk_btree_root(fs, url: str, header_addr: int) -> tuple[int, int]:
    """Root address of a chunked dataset's v1 B-tree and its key rank (rank + 1).

    Read from the Data Layout message (type 8, version 3, class 2 = chunked) in
    the dataset's object header, following continuation blocks (type 16).
    """
    b = fs.cat_file(url, start=header_addr, end=header_addr + 2**16)
    version = 2 if b[:4] == b"OHDR" else 1
    flags = b[5] if version == 2 else 0
    blocks, first = [b], True
    while blocks:
        block = blocks.pop(0)
        for mtype, body in _header_messages(block, version, first, flags):
            if mtype == 0x08:
                if body[0] != 3 or body[1] != 2:
                    raise NotImplementedError(f"data layout v{body[0]} class {body[1]}")
                return _u(body, 3, 8), body[2]
            if mtype == 0x10:
                addr, length = _u(body, 0, 8), _u(body, 8, 8)
                blocks.append(fs.cat_file(url, start=addr, end=addr + length))
        first = False
    raise NotImplementedError("no data layout message in object header")


def _btree_chunk_addresses(fs, url: str, root: int, key_rank: int, wanted: list[tuple]) -> dict:
    """Byte address, size and filter mask of each wanted chunk, from an HDF5 v1 chunk B-tree.

    The tree is walked one level at a time; every node of a level that can hold
    a wanted chunk is fetched in one concurrent batch, so the cost is a few
    round trips instead of hundreds of sequential reads. Chunks missing from
    the tree were never written (all fill value).
    """
    rank = key_rank - 1
    key_size = 8 + 8 * key_rank
    node_guess = 24 + 64 * (key_size + 8) + key_size  # default K=32: 64 entries
    found: dict[tuple, tuple[int, int, int]] = {}
    level = [(root, sorted(wanted))]
    while level:
        blobs = fs.cat_ranges([url] * len(level), [a for a, _ in level],
                              [a + node_guess for a, _ in level])
        nxt = []
        for (addr, want), b in zip(level, blobs):
            if isinstance(b, Exception):
                raise b
            if b[:4] != b"TREE" or b[4] != 1:
                raise NotImplementedError("not a v1 raw-data chunk B-tree node")
            node_level, n = b[5], _u(b, 6, 2)
            need = 24 + n * (key_size + 8) + key_size
            if need > len(b):
                b = fs.cat_file(url, start=addr, end=addr + need)
            keys, kids, pos = [], [], 24
            for i in range(n + 1):
                offset = tuple(_u(b, pos + 8 + 8 * j, 8) for j in range(rank))
                keys.append((offset, _u(b, pos, 4), _u(b, pos + 4, 4)))
                pos += key_size
                if i < n:
                    kids.append(_u(b, pos, 8))
                    pos += 8
            if node_level == 0:
                wanted_here = set(want)
                for (offset, size, mask), child in zip(keys[:n], kids):
                    if offset in wanted_here:
                        found[offset] = (child, size, mask)
            else:
                for i, child in enumerate(kids):
                    lo, hi = keys[i][0], keys[i + 1][0]
                    sub = [c for c in want if lo <= c <= hi]  # inclusive: leaves match exactly
                    if sub:
                        nxt.append((child, sub))
        level = nxt
    return found


def read_box_chunked(fs, url: str, variables: list[str], region: Box, max_gap: int = 2**16) -> xr.Dataset:
    """Read only the HDF5 chunks of ``variables`` that touch ``region``, many at once.

    A PACE L3M AOP file stores Rrs in (16 lat, 1024 lon, 8 band) gzip chunks,
    and our box touches about 1,700 of them (~20 MB). Reading through
    ``xarray`` + a block-cached file object pulls one 4 MB block per chunk
    touched, which is gigabytes per day. Here the chunk byte offsets are looked
    up in the file's index, the byte ranges fetched concurrently, then inflated
    and un-shuffled locally. Values are kept as stored (int16) and decoded by
    the scale_factor/add_offset attributes, exactly as NASA packs them.
    """
    import h5py

    plans, coords, gattrs = {}, {}, {}
    with fs.open(url, block_size=2**20) as f, h5py.File(f, "r") as h:
        lat = h["lat"][:]
        lon = h["lon"][:]
        iy = np.flatnonzero((lat >= region.lat_min) & (lat <= region.lat_max))
        ix = np.flatnonzero((lon >= region.lon_min) & (lon <= region.lon_max))
        y0, y1, x0, x1 = int(iy[0]), int(iy[-1]) + 1, int(ix[0]), int(ix[-1]) + 1
        coords["lat"] = lat[y0:y1]
        coords["lon"] = lon[x0:x1]
        gattrs = {k: _attr(v) for k, v in h.attrs.items() if not k.startswith("_")}
        for var in variables:
            d = h[var]
            if d.compression != "gzip" or d.fletcher32 or d.chunks is None:
                raise NotImplementedError(f"{var}: unsupported storage layout")
            extra = d.shape[2:]
            dims = ["lat", "lon"]
            for axis in range(2, d.ndim):
                name = d.dims[axis][0].name.lstrip("/")
                dims.append(name)
                coords[name] = h[name][:]
            ranges = [range(y0 // d.chunks[0] * d.chunks[0], y1, d.chunks[0]),
                      range(x0 // d.chunks[1] * d.chunks[1], x1, d.chunks[1])]
            ranges += [range(0, n, c) for n, c in zip(extra, d.chunks[2:])]
            origins = [tuple(int(o) for o in origin) for origin in
                       np.stack(np.meshgrid(*ranges, indexing="ij"), -1).reshape(-1, d.ndim)]
            # Header address from the root group's link. h5o.get_info would also
            # report the index size, which walks the whole chunk B-tree remotely.
            header_addr = h.id.links.get_info(var.encode()).u
            root, key_rank = _chunk_btree_root(fs, url, header_addr)
            index = _btree_chunk_addresses(fs, url, root, key_rank, origins)
            chunks = []
            for origin in origins:
                addr, size, mask = index.get(origin, (None, 0, 0))
                if mask:
                    raise NotImplementedError(f"{var}: chunk with skipped filters")
                chunks.append((origin, addr, size))
            plans[var] = dict(
                chunks=chunks, chunk_shape=d.chunks, dtype=d.dtype, shuffle=d.shuffle,
                fill=d.fillvalue, dims=dims, shape=(y1 - y0, x1 - x0, *extra),
                attrs={k: _attr(v) for k, v in d.attrs.items()
                       if k not in ("DIMENSION_LIST", "_Netcdf4Dimid", "_Netcdf4Coordinates")},
            )

    spans = sorted({(off, off + size) for p in plans.values() for _, off, size in p["chunks"] if off is not None})
    merged = _merge_ranges(spans, max_gap)
    blobs = fs.cat_ranges([url] * len(merged), [s for s, _ in merged], [e for _, e in merged])
    for b in blobs:
        if isinstance(b, Exception):
            raise b
    starts = np.array([s for s, _ in merged])

    def chunk_bytes(off: int, size: int) -> bytes:
        i = int(np.searchsorted(starts, off, side="right")) - 1
        rel = off - merged[i][0]
        return blobs[i][rel: rel + size]

    data_vars = {}
    origin0 = (y0, x0)
    for var, p in plans.items():
        out = np.full(p["shape"], p["fill"], dtype=p["dtype"])
        cs = p["chunk_shape"]
        for origin, off, size in p["chunks"]:
            if off is None:
                continue  # never written: the whole chunk is fill value
            raw = zlib.decompress(chunk_bytes(off, size))
            flat = _unshuffle(raw, p["dtype"]) if p["shuffle"] else np.frombuffer(raw, p["dtype"])
            block = flat.reshape(cs)
            dst, src = [], []
            for axis, o in enumerate(origin):
                lo = origin0[axis] if axis < 2 else 0
                a = max(o, lo)
                b = min(o + cs[axis], lo + p["shape"][axis])
                dst.append(slice(a - lo, b - lo))
                src.append(slice(a - o, b - o))
            out[tuple(dst)] = block[tuple(src)]
        data_vars[var] = (p["dims"], out, p["attrs"])
    ds = xr.Dataset(data_vars, coords={k: v for k, v in coords.items()}, attrs=gattrs)
    return xr.decode_cf(ds)


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
        elif g.size_mb > REMOTE_READ_ABOVE_MB and spec.get("variables"):
            _write_netcdf(read_box_chunked(fs, g.url, spec["variables"], region), out)
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
