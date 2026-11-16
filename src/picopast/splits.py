"""Leakage-safe splits (PRD section 7).

Anything fitted to data (scalers, PCA bases, climatologies, models) must be
fitted on the training part of each split only. These functions give the
indices; the modelling code is responsible for fitting inside them.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import pandas as pd

from picopast.grid import block_id


def spatial_block_folds(
    lat, lon, n_folds: int = 5, block_deg: float = 2.0, seed: int = 0
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield (train_idx, test_idx); every 2x2-degree block falls entirely in one fold.

    Blocks are dealt to folds largest-first so folds end up with similar sample
    counts even when cloud leaves blocks very unequal.
    """
    blocks = block_id(lat, lon, block_deg)
    labels, counts = np.unique(blocks, return_counts=True)
    if labels.size < n_folds:
        raise ValueError(f"only {labels.size} blocks for {n_folds} folds")
    rng = np.random.default_rng(seed)
    order = np.lexsort((rng.random(labels.size), -counts))
    fold_of, load = {}, np.zeros(n_folds)
    for i in order:
        f = int(np.argmin(load))
        fold_of[labels[i]] = f
        load[f] += counts[i]
    fold = np.array([fold_of[b] for b in blocks])
    for f in range(n_folds):
        yield np.flatnonzero(fold != f), np.flatnonzero(fold == f)


def temporal_holdout(time, months: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """(train_idx, test_idx) with the latest ``months`` of the record held out."""
    t = pd.DatetimeIndex(pd.to_datetime(np.asarray(time)))
    cut = t.max() - pd.DateOffset(months=months)
    test = np.asarray(t > cut)
    return np.flatnonzero(~test), np.flatnonzero(test)
