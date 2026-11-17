"""Route B (two-step): rebuild a PACE-like spectrum from the VIIRS bands, then run MOANA on it.

PRD section 6 step 5. A small MLP maps the five virtual-VIIRS bands (PACE
spectra seen through SNPP's responses) to the leading principal components of
the PACE spectrum over 400-680 nm; the rebuilt spectrum goes through our MOANA
with SST. Everything fitted (scalers, PCA, MLP, climatology) is fitted inside
the training part of each split.

Scored two ways against our MOANA on the PACE spectrum:
- virtual VIIRS (the RQ2 table): comparable with RQ2's direct models;
- real VIIRS SNPP pixels (the reconstruction's training table), with the MLP
  still trained on virtual bands only: comparable with Route A.
Both tables share one block-to-fold assignment and the same temporal cut.
MOANA on a rebuilt spectrum can run far outside the training range, so
predictions are clipped to the range of the training targets (as any fitted
model's output is bounded); the share clipped is reported.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from picopast.config import Config
from picopast.reconstruct import GROUPS, SeasonalClimatology, _scores, to_target

NAMES = {"pro": "Prochlorococcus", "syn": "Synechococcus", "peuk": "picoeukaryotes"}
N_PC = 12


def _wavelengths(cols: list[str]) -> np.ndarray:
    return np.array([float(c[1:].split("_")[0]) for c in cols])


def fold_ids(lat_a, lon_a, lat_b, lon_b, n_folds: int = 5, block_deg: float = 2.0, seed: int = 0):
    """One block-to-fold assignment shared by two tables."""
    from picopast.grid import block_id

    ba, bb = block_id(lat_a, lon_a, block_deg), block_id(lat_b, lon_b, block_deg)
    blocks = np.unique(np.concatenate([ba, bb]))
    fold_of = dict(zip(blocks, np.random.default_rng(seed).permutation(blocks.size) % n_folds))
    return np.array([fold_of[b] for b in ba]), np.array([fold_of[b] for b in bb])


class SpectrumRebuilder:
    """Five bands (+ log ratios) -> standardised PCs of the 400-680 nm spectrum -> spectrum."""

    def fit(self, bands: np.ndarray, spectra: np.ndarray, seed: int = 0):
        from sklearn.decomposition import PCA
        from sklearn.neural_network import MLPRegressor
        from sklearn.preprocessing import StandardScaler

        self.xs = StandardScaler().fit(self._x(bands))
        self.pca = PCA(N_PC, random_state=seed).fit(spectra)
        scores = self.pca.transform(spectra)
        self.ys = StandardScaler().fit(scores)
        self.mlp = MLPRegressor(hidden_layer_sizes=(64, 64), early_stopping=True, max_iter=200, random_state=seed,
                                batch_size=2048, learning_rate_init=1e-3)
        self.mlp.fit(self.xs.transform(self._x(bands)), self.ys.transform(scores))
        return self

    @staticmethod
    def _x(bands):
        b = np.clip(bands, 1e-5, None)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.log10(b[:, [0, 1, 2, 4]] / b[:, [1, 3, 3, 3]])
        return np.column_stack([bands * 100, ratios])

    def predict(self, bands: np.ndarray) -> np.ndarray:
        z = self.ys.inverse_transform(self.mlp.predict(self.xs.transform(self._x(bands))))
        return self.pca.inverse_transform(z)


def _moana_targets(wl, spectra, sst, pc_mapping) -> dict[str, np.ndarray]:
    from picopast.moana.algorithm import run_moana

    out = run_moana(wl, spectra, sst, pc_mapping)
    return {g: to_target(g, out[g]) for g in GROUPS}


def evaluate(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    from picopast import rq2
    from picopast.reconstruct import BANDS, training_table

    mapping = cfg["reconstruction"]["mapping"]
    v = rq2.table(cfg)
    r = training_table(cfg, "snpp", mapping, cfg["reconstruction"]["per_period_sample"])
    pcols = [c for c in v.columns if c.startswith("p") and "_" in c]
    wl_all = _wavelengths(pcols)
    keep = (wl_all >= 400) & (wl_all <= 680)
    pcols, wl = [c for c, k in zip(pcols, keep) if k], wl_all[keep]
    vbands = [f"v{b}" for b in BANDS]
    rbands = [f"rrs{b}" for b in BANDS]
    fv, fr = fold_ids(v["lat"].values, v["lon"].values, r["lat"].values, r["lon"].values,
                      block_deg=cfg["grid"]["block_deg"])
    cut = v["time"].max() - pd.DateOffset(months=cfg["periods"]["temporal_holdout_months"])
    splits = [("spatial", fv != k, fv == k, fr != k, fr == k) for k in range(5)]
    splits.append(("temporal", (v["time"] <= cut).values, (v["time"] > cut).values,
                   (r["time"] <= cut).values, (r["time"] > cut).values))
    wheres = ("virtual", "real", "real-trained")
    preds = {(kind, where, g): [] for kind in ("spatial", "temporal") for where in wheres for g in GROUPS}
    vreal = [f"r{b}" for b in BANDS]
    has_real = np.isfinite(v[vreal].to_numpy()).all(axis=1)
    spec_err = []
    for kind, vtr, vte, rtr, rte in splits:
        rb = SpectrumRebuilder().fit(v.loc[vtr, vbands].to_numpy(np.float64), v.loc[vtr, pcols].to_numpy(np.float64))
        rebuilt_v = rb.predict(v.loc[vte, vbands].to_numpy(np.float64))
        truth = v.loc[vte, pcols].to_numpy(np.float64)
        spec_err.append(np.median(np.abs(rebuilt_v - truth) / np.abs(truth), axis=0))
        rebuilt_r = rb.predict(r.loc[rte, rbands].to_numpy(np.float64))
        # The PRD's Route B: the same rebuild, trained on real VIIRS bands paired with the PACE spectrum.
        tr_real = vtr & has_real
        rb_real = SpectrumRebuilder().fit(v.loc[tr_real, vreal].to_numpy(np.float64),
                                          v.loc[tr_real, pcols].to_numpy(np.float64))
        rebuilt_rr = rb_real.predict(r.loc[rte, rbands].to_numpy(np.float64))
        est = {"virtual": _moana_targets(wl, rebuilt_v, v.loc[vte, "sst"].to_numpy(), mapping),
               "real": _moana_targets(wl, rebuilt_r, r.loc[rte, "sst"].to_numpy(), mapping),
               "real-trained": _moana_targets(wl, rebuilt_rr, r.loc[rte, "sst"].to_numpy(), mapping)}
        tables = {"virtual": (v, vtr, vte), "real": (r, rtr, rte), "real-trained": (r, rtr, rte)}
        for where, (tab, tr, te) in tables.items():
            for g in GROUPS:
                y_tr, y_te = tab.loc[tr, g].to_numpy(), tab.loc[te, g].to_numpy()
                ok_tr, ok_te = np.isfinite(y_tr), np.isfinite(y_te)
                clim = SeasonalClimatology().fit(tab.loc[tr][ok_tr], y_tr[ok_tr]).predict(tab.loc[te][ok_te])
                p = est[where][g][ok_te]
                lo, hi = np.nanmin(y_tr[ok_tr]), np.nanmax(y_tr[ok_tr])
                clipped = ~np.isfinite(p) | (p < lo) | (p > hi)
                p = np.clip(np.where(np.isfinite(p), p, np.nanmedian(y_tr[ok_tr])), lo, hi)
                preds[(kind, where, g)].append((y_te[ok_te], p, clim, clipped))
    rows = []
    for (kind, where, g), parts in preds.items():
        y, p, c, clipped = (np.concatenate([q[i] for q in parts]) for i in range(4))
        rows.append({"group": g, "split": kind, "inputs": where, "n": int(y.size), "clipped": float(clipped.mean()),
                     **_scores(y, p, c)})
    spectrum = pd.DataFrame({"wavelength": wl, "median_rel_error": np.median(np.vstack(spec_err), axis=0)})
    return pd.DataFrame(rows), spectrum


def report(cfg: Config, refresh: bool = False) -> Path:
    work = cfg.path("interim") / "route_b"
    work.mkdir(parents=True, exist_ok=True)
    if refresh or not (work / "scores.csv").exists():
        scores, spectrum = evaluate(cfg)
        scores.to_csv(work / "scores.csv", index=False)
        spectrum.to_csv(work / "spectrum_error.csv", index=False)
    scores, spectrum = pd.read_csv(work / "scores.csv"), pd.read_csv(work / "spectrum_error.csv")
    rq2 = pd.read_csv(cfg.path("interim") / "rq2" / "scores.csv")
    ra = pd.read_csv(cfg.path("interim") / "reconstruct" / "validation_snpp_operational.csv")
    L = ["# Route B: rebuild the spectrum, then run MOANA", "",
         f"Generated {date.today().isoformat()} by `picopast route-b`. An MLP (two layers of 64) learns the leading "
         f"{N_PC} principal components of the PACE spectrum (400–680 nm) from the five virtual-VIIRS bands; MOANA "
         "then runs on the rebuilt spectrum with SST. The MLP never sees real VIIRS. Scored against our MOANA on "
         "the PACE spectrum, with the same spatial blocks and temporal cut as Routes A and C; predictions clipped "
         "to the training range (share shown).", "",
         "Rebuilt spectrum, median relative error on held-out pixels: " + ", ".join(
             f"{w:.0f} nm {e:.1%}" for w, e in spectrum.iloc[::10][["wavelength", "median_rel_error"]].values) + ".",
         "", "Two versions on real VIIRS: B trained on virtual bands only (never sees real VIIRS), and B trained on "
         "real VIIRS SNPP bands paired with the PACE spectrum on the same pixels and 8-day periods (the PRD's "
         "intent).", "",
         "| group | split | B, virtual VIIRS | direct model, virtual VIIRS (RQ2) | B (virtual-trained), real VIIRS | "
         "B (real-trained), real VIIRS | Route A, real VIIRS | Route C, real VIIRS |",
         "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for g in GROUPS:
        for split in ("spatial", "temporal"):
            s = scores[(scores["group"] == g) & (scores["split"] == split)].set_index("inputs")
            d = rq2[(rq2["group"] == g) & (rq2["split"] == split) & (rq2["model"] == "LightGBM")
                    & (rq2["inputs"] == "virtual VIIRS (5 bands)")]["skill"].iloc[0]
            a = ra[(ra["group"] == g) & (ra["split"] == split)].set_index("model")["skill"]
            L.append(f"| {NAMES[g]} | {split} | {s.loc['virtual', 'skill']:.2f} | {d:.2f} | "
                     f"{s.loc['real', 'skill']:.2f} | {s.loc['real-trained', 'skill']:.2f} | {a['A']:.2f} | "
                     f"{a['C']:.2f} |")
    clip = scores.groupby("inputs")["clipped"].mean()
    L += ["", "Share of predictions clipped to the training range: " + ", ".join(
        f"{k} {v:.0%}" for k, v in clip.items()) + "."]
    L += ["", "Skill = 1 − RMSE / RMSE of a seasonal climatology fitted in the same training part (> 0 beats it).", ""]
    out = cfg.path("reports") / "route_b.md"
    out.write_text("\n".join(L))
    return out
