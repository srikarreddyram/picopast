# PicoPast

Satellite-based evaluation and back-extension of PACE picophytoplankton
abundance estimates for the North Indian Ocean (0–25°N, 50–100°E), 2012–2026.

Major project, B.Tech Computer Engineering, SRM Institute of Science and
Technology. Team: Aryan Sinha, Ramachandra TejSrikar Reddy. Guide: Dr. M.
Suganiya. The full plan is in [PicoPast_PRD.md](PicoPast_PRD.md).

We run NASA's MOANA algorithm ourselves over the North Indian Ocean, test it
against ship counts, and (if it passes) extend the estimates back to 2012 with
VIIRS. MOANA output is treated as a model estimate, not truth, and climate
results are reported as associations, not causes.

## Setup

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m pytest          # 29 tests, a few seconds, no network
```

NASA downloads need a free Earthdata account (<https://urs.earthdata.nasa.gov>).
Log in once; the credentials go to `~/.netrc`, never into this repository:

```bash
.venv/bin/picopast login
```

## Commands

| Command | What it does | Login |
| --- | --- | --- |
| `picopast check-access` | Lists what NASA holds per product for a sample month → `docs/reports/data_access.md` | no |
| `picopast inventory` | Counts ship measurements in BCO-DMO, PANGAEA and SeaBASS → `docs/reports/data_inventory.md` | no |
| `picopast fetch-open` | Downloads DMI, ONI and OISST v2.1 for the box → `data/raw/` | no |
| `picopast fetch <product> <start> <end>` | Downloads a NASA product cut to the box → `data/raw/<product>/` | yes |
| `picopast moana-check <date>` | Compares our MOANA with NASA's Atlantic L4 file for one day | yes |

Products are named in [config/picopast.yaml](config/picopast.yaml), which is
the single place that defines the box, regions, seasons, periods, matchup
settings and trust-gate thresholds.

## Layout

```
config/picopast.yaml      the one study config
src/picopast/
  moana/                  our MOANA run, NASA's tables, the Atlantic check
  ingest/                 NASA catalogue search, box-subset downloads, OISST, DMI, ONI
  insitu/                 ship-count archives, common table layout, inventory
  matchup.py              3x3 median, same day then ±1 day
  metrics.py              log-RMSE, bias, Spearman, factor-of-2, cruise bootstrap, skill
  trust_gate.py           RQ1 verdict per group
  splits.py               spatial-block folds, temporal hold-out
  virtual_viirs.py        PACE spectra through VIIRS band responses
docs/
  trust_gate_preregistration.md   the RQ1 rule, to be signed before week 3
  moana_notes.md                  what the algorithm really does, and what is verified
  open_items.md                   decisions and blockers
  reports/                        generated: data access, ship-count inventory
data/                     not in git; rebuilt from the config
```

## Status (10 November 2026, week 1)

Done: repository and environment; our MOANA implementation, checked against an
independent reference written from the ATBD equations; data-access check;
ship-count inventory of the open archives; DMI, ONI and OISST downloads;
matchup, metrics, trust-gate, split and virtual-VIIRS code with tests.

Not done yet: anything that needs the Earthdata login (PACE and VIIRS
downloads, the Atlantic MOANA check, opening SeaBASS files); the Zarr cube;
models. See [docs/open_items.md](docs/open_items.md).

## Rules we hold ourselves to

- Scalers, PCA bases and climatologies are fitted inside training folds only.
- No synthetic inputs; observed or reanalysis data only. (Tests use made-up
  arrays to check arithmetic; no result is ever drawn from them.)
- The trust-gate rule is fixed before any matchup is computed.
- Novelty is worded as "limited evidence in the literature we reviewed".
- Ship data come from open archives and published tables; nobody is contacted.
