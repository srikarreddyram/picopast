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

**New direction (22 November 2026): how marine heatwaves change the phytoplankton
response to cyclones.** The project now asks whether a marine heatwave changes how
phytoplankton respond to tropical cyclones in the North Indian Ocean (83 storms,
2012–2026), with NASA satellite data and PicoPast's own picophytoplankton records
(the work below), which answer RQ3 for every storm. Literature check: [docs/literature_survey.md](docs/literature_survey.md);
rules fixed before testing: [docs/cyclone_heatwave_preregistration.md](docs/cyclone_heatwave_preregistration.md);
first results: [docs/reports/cyclone_heatwave.md](docs/reports/cyclone_heatwave.md).

First pass: storms over heatwave water cool the sea more (-0.57 °C
[-0.86, -0.25]; about half of that remains after allowing for heatwaves
decaying anyway) and end most of the heatwaves they cross. Whether the heatwave changes
the plankton bloom afterwards is not settled: the pre-registered test gives
-0.058 [-0.125, +0.004] (log10 chlorophyll), leaning towards a
weaker bloom with the interval just including zero.
After storms, picoeukaryotes (+0.15 log10) and *Synechococcus*
(+0.14) rise and *Prochlorococcus* (-0.16) falls in PicoPast's
records built without sea temperature: a shift towards larger cells, with intervals
clear of zero. A heatwave does not detectably change that shift.
A pre-registered replication on 50 earlier storms (2002–2011, MODIS-Aqua) cannot settle
the bloom question either: heatwaves were rare then (6% of track points against
30% after 2012). Pooled over 2002–2026 the bloom is about 11% weaker over
heatwave water (95% interval: 23% weaker to 1% stronger), so no detectable
difference; a large effect either way, and any real strengthening, are ruled out.

**Definitive test (27 November 2026): no meaningful effect.**
- **Design:** the pre-registered final test adds the South Indian Ocean (360 more
  storms) and is judged on a 99% interval.
- **Result:** across 448 storms with a usable bloom, the heatwave effect on the bloom
  is −0.008 [−0.030, +0.014] log10 (−2% chlorophyll, between −7% and +3%). The whole
  interval lies inside the ±0.05 threshold fixed beforehand.
- **Caveats:**
  - The North Indian Ocean on its own remains inconclusive (27 heatwave storms).
  - Measured against the same water in years without a storm, the storm's own
    contribution is about 8% smaller over heatwaves (2–12%, Holm p 0.016).
- Report: [docs/reports/cyclone_heatwave_final.md](docs/reports/cyclone_heatwave_final.md).

**Why (virtual experiments, 28 November 2026).**
- **Method:** a 1-D ocean column model (PWP mixing plus NPZD biology) runs 484 storms
  (610 storm points) twice: over the real pre-storm column, and with the heatwave removed
  (or, for normal storms, added).
- **Result:** in the model, the heatwave weakens the bloom by about 7% (−0.031
  [−0.042, −0.021] log10). The warm surface layer makes the water more stable: the
  storm mixes 4 m less deep and brings up less nitrate.
- **Agreement:** the size matches the observed storm-minus-no-storm comparison
  (about 8%).
- **Robustness:** every sensitivity variant keeps the sign.
- Report: [docs/reports/virtual_experiments.md](docs/reports/virtual_experiments.md).

**Storms end heatwaves (1982–2026, 28 November 2026).**
- **Method:** pre-registered; 45 years of daily heatwave maps and 915 storms over the
  North and South Indian Ocean.
- **Results:**
  - A heatwave patch under a storm is 3.5 times as likely to end that day (3.1–3.8),
    and 4.6 times for storms of 64 kt and above.
  - Five days after a storm, 21% of its heatwave water is still in a heatwave, against
    45% far away. This holds in 1982–2001 alone (15% against 34%).
  - Stronger, slower storms end more of the heatwave.
  - Storms cut heatwave days in the following month by 30%, but not beyond that.
    Storm-ended heatwaves do not stay away longer than others.
- **In context:** storms are rare, so they cause only 0.4% of all heatwave endings.
- Report: [docs/reports/heatwave_end.md](docs/reports/heatwave_end.md).

## Setup

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m pytest          # 92 tests, no network
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
| `picopast cyclone-heatwave [--refresh]` | Builds the storm table and runs the pre-registered cyclone–heatwave tests → `docs/reports/cyclone_heatwave.md` | no (data on disk) |
| `picopast final-test [--build]` | Final confirmatory test over the whole Indian Ocean (T1–T5) → `docs/reports/cyclone_heatwave_final.md` | no (data on disk) |
| `picopast heatwave-end [--refresh]` | Storms end heatwaves, 1982–2026: ending rate, return time, suppression (pre-registered) → `docs/reports/heatwave_end.md` | no (data on disk) |
| `picopast virtual [--workers 6]` | Column-model virtual experiments, each storm with and without its heatwave → `docs/reports/virtual_experiments.md` | no (data on disk) |

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

## Status (16 November 2026, week 2)

Done: repository and environment; our MOANA implementation, checked against an
independent reference written from the ATBD equations and against NASA's own
Atlantic L4 output (within 1% median for all three groups); data-access check;
ship-count inventory of the open archives; DMI, ONI and the complete OISST
record (1991 to November 2026); a chunk-level PACE reader that fetches only the
study box; matchup, metrics, trust-gate, split and virtual-VIIRS code with
tests.

Trust gate (binding rule, run 16 November): there are no PACE-era cell counts
in the box, so the pre-registered fallback checks decided which groups go
forward. **Picoeukaryotes pass** (same-month pattern vs 1991–2005 counts,
ρ = 0.51 [0.24, 0.64]) and are reconstructed as a relative index.
Synechococcus and Prochlorococcus are inconclusive with MOANA; they pass with Route D (below).
Details in [docs/reports/trust_gate_fallback.md](docs/reports/trust_gate_fallback.md).

Back-extension, first attempt (17 November): `picopast reconstruct` learns
MOANA from VIIRS SNPP 8-day inputs in the PACE era and rebuilds all three
groups for 2012–2026, each tagged with its trust-gate tier (Synechococcus and
Prochlorococcus are exploratory here; use their Route D record, below). Route A beats a seasonal climatology on
held-out blocks for every group (skill 0.39 picoeukaryotes, 0.20
Synechococcus, 0.18 Prochlorococcus); the SST + PAR route does not. Details
and bottlenecks in [docs/reports/reconstruction.md](docs/reports/reconstruction.md).
An independent NOAA-20 copy (RQ4) agrees with SNPP for picoeukaryotes to
within ±0.07 SD over 2018–2026; the two cyanobacteria records drift apart
away from the training years
([docs/reports/sensor_comparison.md](docs/reports/sensor_comparison.md)). A
chlorophyll climatology matches the 1991–2005 counts about as well as MOANA
picoeukaryotes and better than MOANA's cyanobacteria
([docs/reports/chlorophyll_comparison.md](docs/reports/chlorophyll_comparison.md)).

All data are downloaded: PACE (daily and 8-day, to NASA's last file, 30 September
2026), VIIRS SNPP and NOAA-20, CMC SST and OISST.

Climate associations (RQ5, first pass): the picoeukaryote index has no
2012–2025 trend, is lower after El Niño and positive IOD (r ≈ −0.55 to −0.6,
same on NOAA-20) and lower during marine heatwaves (−0.25 to −0.36 SD); most
of each association survives a model without SST
([docs/reports/climate.md](docs/reports/climate.md)). Mapped per 0.25° cell, the
El Niño and IOD drops sit in the Arabian Sea while most of the Bay of Bengal does
not move; heatwave drops cover almost the whole region; both survive without SST
and agree on NOAA-20 ([docs/reports/climate_maps.md](docs/reports/climate_maps.md)).

RQ2: the five VIIRS bands cost every group about 0.15 of skill against the full
PACE spectrum; the larger loss, for the cyanobacteria, comes from the change of
sensor itself ([docs/reports/rq2_information.md](docs/reports/rq2_information.md),
[docs/reports/virtual_viirs.md](docs/reports/virtual_viirs.md)).

Route D (cyanobacteria): MOANA, fitted to Atlantic cruises, does not rank the
Indian Ocean stations for *Synechococcus* or *Prochlorococcus*. A model trained on
ship counts from the rest of the world ocean (MAREDAT, Indian Ocean held out) and
weighted towards this box's conditions does: ρ 0.50 and 0.45, both passing the
same pre-registered test, and clearly better than MOANA; it is no better than
chlorophyll alone ([docs/reports/route_d.md](docs/reports/route_d.md)).
Both groups are now rebuilt from it for 2012–2026 as a relative index
(`picopast route-d-reconstruct`). The two satellites agree, and the old records' rises
over 2012–2025 were sensor drift. *Prochlorococcus* rises during marine heatwaves, with
or without SST as an input. *Synechococcus* falls, but mostly through the SST input
([docs/reports/route_d_reconstruction.md](docs/reports/route_d_reconstruction.md)).

Check E (pre-registered): 149 new HPLC pigment stations from 12 cruises (2007–2025). Route D's
*Prochlorococcus* record does not rank divinyl chlorophyll a (ρ −0.07 [−0.30, +0.16]), so
*Prochlorococcus* is **exploratory** again; the pigment follows total chlorophyll, opposite to the
cell counts, so the two tests disagree and only new counts can settle it
([docs/reports/check_e.md](docs/reports/check_e.md)). *Synechococcus* stays a relative index.

**Everything side by side** (which record each group uses, and where performance
improved): [docs/reports/comparison.md](docs/reports/comparison.md).

RQ3: on real VIIRS the direct route (A) beats the two-step route (B: rebuild
the PACE spectrum, then MOANA) for every group, and the SST + PAR route (C)
does no better than climatology ([docs/reports/route_b.md](docs/reports/route_b.md)).

Not done yet: mixed-layer depth (Copernicus account), the MLP version of
Route A. See
[docs/open_items.md](docs/open_items.md).

## Rules we hold ourselves to

- Scalers, PCA bases and climatologies are fitted inside training folds only.
- No synthetic inputs; observed or reanalysis data only. (Tests use made-up
  arrays to check arithmetic; no result is ever drawn from them.)
- The trust-gate rule is fixed before any matchup is computed.
- Novelty is worded as "limited evidence in the literature we reviewed".
- Ship data come from open archives and published tables; nobody is contacted.
