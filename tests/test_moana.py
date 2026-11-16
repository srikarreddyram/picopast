"""MOANA code checks.

The reference below is written straight from the ATBD v1.2 equations and
Table 1 (coefficients typed by hand, per-spectrum ``np.interp``), so it shares
no code and no coefficient parsing with ``picopast.moana``. The spectra are
smooth test curves, not ocean data: they exercise the arithmetic only. The
scientific check is the Atlantic comparison against NASA's own L4 files.
"""

import numpy as np
import pytest
import xarray as xr

from picopast.moana import (
    FLAG_BAD_SST,
    FLAG_MISSING_RRS,
    FLAG_NEGATIVE_PRO,
    compare_with_nasa_l4,
    load_luts,
    moana_on_grid,
    run_moana,
)

# ATBD v1.2 Table 1 (s3 and s11 are printed there to four decimals; the fifth
# is from NASA's file). Principal-component indices follow the operational file
# (Pro: last term on U7; Syn: last term on U16); the ATBD text says U17 and U13.
PEUK = (3.31202, {2: 0.21332, 3: 0.09385, 4: 0.72960, 5: 0.97040, 10: 0.67274, 13: -2.25057,
                  15: -3.45864})
SYN = (-17.52889, {1: 1.95949, 2: 0.33859, 3: 0.49127, 5: 1.04400, 8: -1.06453, 9: -1.55528,
                   10: -2.25345, 11: -3.41758, 12: -1.67013, 16: -3.19241})
PRO = (-4562144.0, 770448.0, {1: 343590.0, 2: 54975.0, 6: 290751.0, 7: -1041404.0})


def test_spectra(n=60, seed=7):
    rng = np.random.default_rng(seed)
    wave = np.arange(400.0, 702.5, 2.5)  # PACE-like sampling
    blue = rng.uniform(0.002, 0.012, (n, 1)) * np.exp(-(wave - 400) / rng.uniform(60, 160, (n, 1)))
    green = rng.uniform(0, 0.004, (n, 1)) * np.exp(
        -0.5 * ((wave - rng.uniform(540, 570, (n, 1))) / 35) ** 2
    )
    rrs = blue + green + rng.normal(0, 2e-5, (n, wave.size))
    return wave, rrs, rng.uniform(18, 31, n)


test_spectra.__test__ = False


def reference(wave, rrs, sst, pro_last=7, syn_last=16):
    luts_wave = np.arange(414.0, 661.0, 2.0)
    v = load_luts().loadings
    out = {"pro": [], "syn": [], "peuk": []}
    for spectrum, t in zip(rrs, sst):
        r = np.interp(luts_wave, wave, spectrum)
        r = (r - r.mean()) / r.std(ddof=1)
        u = {i + 1: float(r @ v[:, i]) for i in range(v.shape[1])}
        out["peuk"].append(10 ** (PEUK[0] + sum(c * u[i] for i, c in PEUK[1].items())))
        syn = {(syn_last if i == 16 else i): c for i, c in SYN[1].items()}
        out["syn"].append(10 ** (SYN[0] + sum(c * u[i] for i, c in syn.items())))
        pro = {(pro_last if i == 7 else i): c for i, c in PRO[2].items()}
        out["pro"].append(PRO[0] + PRO[1] * np.log10(t) + sum(c * u[i] for i, c in pro.items()))
    return {k: np.array(x) for k, x in out.items()}


def test_luts_shape_and_orthonormal_loadings():
    luts = load_luts()
    assert luts.npc == 17
    assert luts.wavelength[0] == 414 and luts.wavelength[-1] == 660 and luts.wavelength.size == 124
    gram = luts.loadings.T @ luts.loadings
    assert np.allclose(gram, np.eye(17), atol=1e-5)


@pytest.mark.parametrize("mapping,pro_last,syn_last", [("operational", 7, 16), ("atbd", 17, 13)])
def test_matches_reference_from_atbd_equations(mapping, pro_last, syn_last):
    wave, rrs, sst = test_spectra()
    ours = run_moana(wave, rrs, sst, pc_mapping=mapping)
    ref = reference(wave, rrs, sst, pro_last, syn_last)
    for g in ("pro", "syn", "peuk"):
        np.testing.assert_allclose(ours[g], ref[g], rtol=1e-9, atol=1e-6)


def test_mappings_differ_only_for_pro_and_syn():
    wave, rrs, sst = test_spectra()
    a, b = run_moana(wave, rrs, sst), run_moana(wave, rrs, sst, pc_mapping="atbd")
    np.testing.assert_array_equal(a["peuk"], b["peuk"])
    assert not np.allclose(a["syn"], b["syn"])
    assert not np.allclose(a["pro"], b["pro"])


def test_standardisation_makes_output_scale_invariant():
    wave, rrs, sst = test_spectra(10)
    a, b = run_moana(wave, rrs, sst), run_moana(wave, 3.7 * rrs + 0.001, sst)
    for g in ("pro", "syn", "peuk"):
        np.testing.assert_allclose(a[g], b[g], rtol=1e-8, atol=1e-4)


def test_flags_and_missing_data():
    wave, rrs, sst = test_spectra(6)
    rrs = rrs.copy()
    rrs[0, np.argmin(np.abs(wave - 500))] = np.nan  # inside 414-660 nm
    rrs[1, 0] = np.nan  # 400 nm: outside the fit range, must not matter
    sst = sst.copy()
    sst[2] = np.nan
    out = run_moana(wave, rrs, sst)
    assert out["flags"][0] & FLAG_MISSING_RRS and np.isnan(out["syn"][0])
    assert not out["flags"][1] & FLAG_MISSING_RRS and np.isfinite(out["syn"][1])
    assert out["flags"][2] & FLAG_BAD_SST and np.isnan(out["pro"][2])
    assert np.isfinite(out["syn"][2])  # only Pro needs SST
    negative = out["pro"] < 0
    assert np.array_equal(negative, (out["flags"] & FLAG_NEGATIVE_PRO) > 0)


def test_rejects_kelvin_and_short_wavelength_range():
    wave, rrs, sst = test_spectra(4)
    with pytest.raises(ValueError, match="Kelvin"):
        run_moana(wave, rrs, sst + 273.15)
    keep = wave >= 420
    with pytest.raises(ValueError, match="do not cover"):
        run_moana(wave[keep], rrs[:, keep], sst)


def test_grid_run_and_nasa_comparison():
    wave, rrs, sst = test_spectra(24, seed=3)
    lat = np.array([20.0, 19.0, 18.0, 17.0])
    lon = np.arange(60.0, 66.0)
    rrs_da = xr.DataArray(
        rrs.reshape(4, 6, -1), dims=("lat", "lon", "wavelength"),
        coords={"lat": lat, "lon": lon, "wavelength": wave},
    )
    # SST on a coarser, south-to-north grid to exercise the nearest-neighbour lookup.
    sst_da = xr.DataArray(
        np.full((3, 4), 27.0), dims=("lat", "lon"),
        coords={"lat": [16.5, 18.5, 20.5], "lon": [59.0, 61.0, 63.0, 65.0]},
    )
    ds = moana_on_grid(rrs_da, sst_da, rows_per_chunk=3)
    flat = run_moana(wave, rrs, 27.0)
    np.testing.assert_allclose(ds["syn"].values.ravel(), flat["syn"], rtol=1e-6)
    np.testing.assert_allclose(ds["pro"].values.ravel(), flat["pro"], rtol=1e-6, atol=1)

    # A stand-in for a NASA L4 file: our own output with NASA's clamp-and-truncate.
    nasa = xr.Dataset(
        {
            "prococcus_moana": ds["pro"].clip(min=0).pipe(np.trunc),
            "syncoccus_moana": ds["syn"].pipe(np.trunc),
            "picoeuk_moana": ds["peuk"].pipe(np.trunc),
        }
    )
    report = compare_with_nasa_l4(ds, nasa)
    for g in ("syn", "peuk"):
        assert report[g]["median_abs_dlog10"] == 0.0
