"""The chunk-level PACE reader, checked against h5py on local HDF5 files.

The arrays are made up: they only check that byte offsets, decompression,
un-shuffling, box cutting and unpacking reproduce what h5py reads.
"""

import fsspec
import h5py
import numpy as np
import pytest

from picopast.config import Box
from picopast.ingest.earthdata import (
    _btree_chunk_addresses,
    _chunk_btree_root,
    _merge_ranges,
    _unshuffle,
    read_box_chunked,
)

SCALE, OFFSET, FILL = np.float32(2e-6), np.float32(0.05), np.int16(-32767)


def _make_l3m_like(path, libver):
    """A small file shaped like a PACE L3M AOP file: Rrs(lat, lon, wavelength) int16, gzip + shuffle."""
    lat = (10.0 - (np.arange(40) + 0.5) * 0.5).astype(np.float32)  # north to south, like L3M
    lon = (40.0 + (np.arange(80) + 0.5) * 0.5).astype(np.float32)
    wave = np.linspace(400, 700, 9).astype(np.float32)
    rng = np.random.default_rng(0)
    raw = rng.integers(-30000, 25000, size=(40, 80, 9), dtype=np.int16)
    raw[rng.random(raw.shape) < 0.1] = FILL
    with h5py.File(path, "w", libver=libver) as h:
        for name, values in (("lat", lat), ("lon", lon), ("wavelength", wave)):
            h[name] = values
            h[name].make_scale(name)
        d = h.create_dataset("Rrs", shape=raw.shape, dtype=np.int16, chunks=(4, 8, 3),
                             compression="gzip", shuffle=True, fillvalue=FILL)
        d[:24] = raw[:24]  # rows 24+ are never written: those chunks stay unallocated
        d.attrs["scale_factor"] = SCALE
        d.attrs["add_offset"] = OFFSET
        d.attrs["_FillValue"] = FILL
        for axis, name in enumerate(("lat", "lon", "wavelength")):
            d.dims[axis].attach_scale(h[name])
    raw[24:] = FILL
    return lat, lon, wave, raw


@pytest.mark.parametrize("libver", ["earliest", "v108"])  # v1 and v2 object headers
def test_btree_walk_matches_h5py(tmp_path, libver):
    path = str(tmp_path / "f.h5")
    _make_l3m_like(path, libver)
    fs = fsspec.filesystem("file")
    with h5py.File(path) as h:
        d = h["Rrs"]
        root, key_rank = _chunk_btree_root(fs, path, h.id.links.get_info(b"Rrs").u)
        assert key_rank == 4
        wanted = [(y, x, z) for y in range(0, 40, 4) for x in range(0, 80, 8) for z in range(0, 9, 3)]
        found = _btree_chunk_addresses(fs, path, root, key_rank, wanted)
        for c in wanted:
            info = d.id.get_chunk_info_by_coord(c)
            if info.byte_offset is None:
                assert c not in found
            else:
                assert found[c] == (info.byte_offset, info.size, info.filter_mask)
    assert len(found) == 6 * 10 * 3  # rows 0-23 written: 6 chunk rows


@pytest.mark.parametrize("libver", ["earliest", "v108"])
def test_read_box_chunked_matches_direct_read(tmp_path, libver):
    path = str(tmp_path / "f.h5")
    lat, lon, wave, raw = _make_l3m_like(path, libver)
    box = Box(lat_min=-3.0, lat_max=7.2, lon_min=43.1, lon_max=61.0)  # crosses chunk edges and unwritten rows
    ds = read_box_chunked(fsspec.filesystem("file"), path, ["Rrs"], box)
    iy = np.flatnonzero((lat >= box.lat_min) & (lat <= box.lat_max))
    ix = np.flatnonzero((lon >= box.lon_min) & (lon <= box.lon_max))
    sub = raw[iy[0]: iy[-1] + 1, ix[0]: ix[-1] + 1]
    expected = np.where(sub == FILL, np.nan, sub.astype(np.float32) * SCALE + OFFSET)
    np.testing.assert_array_equal(ds["lat"].values, lat[iy])
    np.testing.assert_array_equal(ds["lon"].values, lon[ix])
    np.testing.assert_array_equal(ds["wavelength"].values, wave)
    assert ds["Rrs"].dims == ("lat", "lon", "wavelength")
    np.testing.assert_array_equal(np.isnan(ds["Rrs"].values), np.isnan(expected))
    np.testing.assert_allclose(ds["Rrs"].values, expected, rtol=1e-6, equal_nan=True)


def test_unshuffle_inverts_byte_shuffle():
    values = np.arange(-500, 500, 7, dtype=np.int16)
    shuffled = values.view(np.uint8).reshape(-1, 2).T.tobytes()
    np.testing.assert_array_equal(_unshuffle(shuffled, np.dtype(np.int16)), values)


def test_merge_ranges():
    spans = [(100, 200), (0, 50), (205, 300), (1000, 1100)]
    assert _merge_ranges(spans, max_gap=10) == [(0, 50), (100, 300), (1000, 1100)]
    assert _merge_ranges(spans, max_gap=0) == [(0, 50), (100, 200), (205, 300), (1000, 1100)]


def test_rsr_parsing_both_release_formats():
    import io
    import zipfile

    from picopast.ingest.rsr import SOURCES, centre_wavelength, parse_rsr_text, table_from_zip

    ng = "    440.000    0.0000\n    445.000    1.0000\n    450.000    0.0000\n"
    dawg = "%DAWG header\n%Band wvl_nm RSR\n   M2  440.0  0.1E-09\n   M2  445.0  0.5\n   M2  450.0  0.0\n"
    w, r = parse_rsr_text(dawg)
    assert w.tolist() == [440.0, 445.0, 450.0] and r[1] == 0.5
    for sensor, text in (("snpp", ng), ("noaa20", dawg)):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            for n in range(1, 6):
                z.writestr(SOURCES[sensor]["member"].format(n=n), text)
        t = table_from_zip(buf.getvalue(), sensor)
        assert t.shape[1] == 5 and t.max().tolist() == [1.0] * 5  # peak-normalised
        assert centre_wavelength(t).iloc[0] == pytest.approx(445.0, abs=0.05)
