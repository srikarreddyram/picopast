# PicoPast: Product Requirements Document and Data Access Guide (v3.1)

*PicoPast: Satellite-Based Evaluation and Back-Extension of PACE Picophytoplankton Abundance Estimates for the North Indian Ocean Using Explainable Machine Learning*

Major Project, B.Tech Computer Engineering · SRM Institute of Science and Technology
Team: Aryan Sinha, Ramachandra TejSrikar Reddy · Batch B38 · Panel 3
Guide: Dr. M. Suganiya, Assistant Professor, Dept. of Computing Technologies
Version 3.1 · 10 November 2026 · supersedes v3 · **Target: all project work wrapped by 31 December 2026**

---

## 1. One-paragraph summary

NASA's PACE satellite (data since March 2024) can estimate the abundance of three picophytoplankton groups (*Prochlorococcus*, *Synechococcus*, picoeukaryotes) with the MOANA algorithm. MOANA was built on Atlantic cruise data, and NASA distributes its group product only for the Atlantic (70°S–70°N, 85°W–25°E). PicoPast will (1) run MOANA ourselves over the North Indian Ocean and test it against ship measurements (the trust gate), (2) measure how much group information survives in the coarser VIIRS bands, (3) compare three routes for extending the estimates back to 2012, (4) validate the extended record with leakage-safe splits, a second satellite and older ship counts, and (5) use the 2012–2026 record to study associations with the monsoon, the Indian Ocean Dipole (IOD), ENSO and marine heatwaves. We report associations, not causes, and we treat MOANA output as a model estimate, not truth.

## 2. What changed

### 2.1 From v2 to v3

| Area | v2 | v3 | Why |
| --- | --- | --- | --- |
| MOANA labels | Download NASA's L4 MOANA product as training labels | **We run MOANA ourselves** on PACE L3M Rrs plus SST, and check our implementation against NASA's Atlantic files on the same dates | NASA's product covers the Atlantic only; nothing exists for the North Indian Ocean |
| MOANA SST input | "Same SST family MOANA uses" | log10 SST is used **for *Prochlorococcus* only**; GHRSST CMC 0.1° starts in 2016, so **NOAA OISST v2.1** is the SST predictor for the whole 2012–2026 record | Corrects the earlier description and keeps one consistent SST source |
| Ship data | Email authors, INCOIS, NCPOR | **Open archives and published tables only** (SeaBASS, NASA PACE validation data, BCO-DMO, PANGAEA, tables and supplements in published papers). No author or institutional requests | Team decision |
| Baselines | Climatology and persistence | **Climatology and ridge regression** | Persistence has no meaning for a reconstruction |
| Grid | 4 km | **4.6 km** (native L3M grid), 8-day composites for the long record, daily for training and matchups | Matches the native mapped products |
| Research questions | RQ4 was SNPP vs NOAA-20 only | **RQ4 = SNPP vs NOAA-20 agreement and 2012–2023 ship counts** | Ship counts from years PACE never saw are the strongest test of the reconstruction |
| Regions | One box | **Phase 1:** Bay of Bengal (pilot) + central/eastern Arabian Sea. **Phase 2 (conditional):** western Arabian Sea upwelling + equatorial band | Keeps scope finishable |
| Seasons | All | Core: north-east monsoon and both inter-monsoons. **South-west monsoon reported separately** (heavy cloud) with valid-pixel counts | Optical sensors cannot see through monsoon cloud |
| Climate methods | Generic | Hobday marine heatwave definition (OISST, 1991–2020 baseline), composites, lagged correlation, Mann–Kendall/Sen's slope (descriptive), FDR correction | Specific and defensible |

### 2.2 From v3 to v3.1

| Area | v3 | v3.1 | Why |
| --- | --- | --- | --- |
| Finish date | Paper by 24 Apr 2027, work spread to April | **All project work (pipeline, results, paper draft, repo, dataset) wrapped by 31 Dec 2026** | Team decision; work every day |
| Scope | Full PRD scope | **Unchanged** (Phase 2 and stretch items remain optional) | Scope is kept; only the schedule is compressed |
| Role of Reviews I–III | Show progress | **Present finished work** and use Jan to Apr for revision and submission | Work is done by December |

## 3. Research questions

| ID | Question | Evidence and pass condition |
| --- | --- | --- |
| RQ1 | How accurate is MOANA in the North Indian Ocean, per group? **(trust gate)** | Log-RMSE, bias, Spearman ρ, % within a factor of 2 vs ship counts, with cruise-level bootstrap intervals. Verdict per group: *usable label / relative index only / excluded* |
| RQ2 | How much group information survives in VIIRS bands? | Skill of identical models on virtual-VIIRS vs full PACE spectra |
| RQ3 | Which extension route is best: direct (A), two-step (B) or niche (C)? | Identical spatial-block and temporal hold-out scores; all three reported even if A loses |
| RQ4 | Is the back-extended record consistent with independent evidence? | SNPP vs NOAA-20 agreement (target ≥80% of open-ocean pixels within stated uncertainty) and comparison with 2012–2023 ship counts |
| RQ5 | How do the groups relate to monsoon phases, IOD, ENSO and marine heatwaves? | Composites, lagged correlations and trends on reliable pixels only; associations, not causation |

Hypotheses (to be tested, not assumed): picoeukaryotes transfer best, cyanobacteria are shakier (RQ1); the biggest information loss is for *Synechococcus* (RQ2); direct learning matches or beats the other routes (RQ3).

## 4. Research gaps addressed (limited evidence in the literature we reviewed)

| Gap | Statement | Answered by |
| --- | --- | --- |
| G1 | MOANA has not been evaluated outside the Atlantic | RQ1 |
| G2 | The group information content of VIIRS bands, and the best way to extend back in time, are untested | RQ2, RQ3 |
| G3 | Long records use random splits and no independent-sensor validation | RQ3, RQ4 |
| G4 | No multi-year group-level satellite record exists for the North Indian Ocean | RQ5 |

Always word novelty carefully: "limited evidence in the literature we reviewed", never "first ever". A logged Scopus / Web of Science / Shodhganga search is task T01.

## 5. Data sources

**Step zero:** create a free NASA Earthdata account at <https://urs.earthdata.nasa.gov> (one login for all NASA datasets; the Python library `earthaccess` uses it) and a free Copernicus Marine account.

### 5.1 Core satellite data (NASA)

| What | Where | Notes |
| --- | --- | --- |
| **PACE OCI L3M Rrs** (hyperspectral reflectance, ~5 nm steps) | Earthdata Search: <https://search.earthdata.nasa.gov> (search "PACE OCI L3M RRS"); OB.DAAC: <https://oceancolor.gsfc.nasa.gov/l3/> | **Main input to our own MOANA run** and to virtual VIIRS. Files are large; subset to the study box |
| **PACE OCI L4 MOANA (Atlantic-only)** | <https://www.earthdata.nasa.gov/data/catalog/ob-cloud-pace-oci-l4m-moana-3.2> · <https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-4%20Mapped/PACE-OCI/> | **Not a training label for us.** Used only to check our MOANA code against NASA's Atlantic files on the same dates. Confirm the current version number when downloading |
| **VIIRS SNPP L3 mapped** (Rrs bands, CHL, Kd490, PAR) | <https://oceancolor.gsfc.nasa.gov/l3/> · Earthdata Search | Main reconstruction sensor, 2012 to Oct 2026 |
| **VIIRS NOAA-20 L3 mapped** | Same as above | Second sensor, 2018 onward, for RQ4 |
| **VIIRS spectral response functions** | <https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables> | Used to build virtual VIIRS; take both SNPP and NOAA-20 tables |
| **MOANA algorithm document (ATBD)** | <https://oceandata.sci.gsfc.nasa.gov/files/atbd/atbd-obdaac-phytoplankton-groups-MOANA.pdf> | Read first; defines the exact inputs and coefficients |
| **PACE data overview** | <https://pace.oceansciences.org/data.htm> | Product list and tutorials |

### 5.2 Physical and climate predictors

| What | Where | Notes |
| --- | --- | --- |
| **SST, NOAA OISST v2.1** | NOAA PSL / NCEI | **SST predictor for 2012–2026** and the basis for marine heatwave detection (1991–2020 baseline) |
| **SST, GHRSST CMC 0.1° L4** | NASA PO.DAAC: <https://podaac.jpl.nasa.gov> (search "CMC0.1deg") | MOANA's SST input for the PACE era (starts 2016). Confirm the version covering 2024–2026 |
| **Mixed-layer depth** | Copernicus Marine, GLORYS12 reanalysis (GLOBAL_MULTIYEAR_PHY_001_030): <https://data.marine.copernicus.eu> | Stratification proxy; reanalysis, not observation |
| **Dipole Mode Index (IOD)** | <https://psl.noaa.gov/gcos_wgsp/Timeseries/DMI/> | Monthly |
| **Oceanic Niño Index (ENSO)** | <https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt> | Monthly |

### 5.3 Ship picoplankton counts (open sources only; no one is contacted)

| What | Where | Use |
| --- | --- | --- |
| **SeaBASS** (NASA in-situ archive) | <https://seabass.gsfc.nasa.gov> (search flow cytometry / picoplankton, 0–25°N, 50–100°E) | Main source. 2024+ counts feed the trust gate; 2012–2023 counts validate the reconstruction |
| **NASA PACE validation, Northern Indian Ocean** | <https://www.earthdata.nasa.gov/data/catalog/ob-daac-pvst-northern-indian-ocean-0> | Use whatever is publicly released; check whether it includes cell counts |
| **BCO-DMO** | <https://www.bco-dmo.org> | Public cruise datasets, mainly older counts |
| **PANGAEA** | <https://www.pangaea.de> | Backup repository for flow-cytometry cell counts |
| **Published cruise studies** (tables and supplements already public) | Nazirahmed et al., JGR Biogeosciences 2026: <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2025JG009023> · Bay of Bengal study, J. Earth Syst. Sci. 2026: <https://link.springer.com/article/10.1007/s12040-026-02759-7> · Strauss 2023 · Wei 2019 | Station values read from public tables and supplements only |
| **MAREDAT compilation** | <https://essd.copernicus.org/articles/4/37/2012/> · <https://doi.pangaea.de/10.1594/PANGAEA.777385> | Mostly pre-2012; climatology and niche-model baseline only |
| **Updated global compilation** | Visintini et al., L&O Letters 2021: <https://doaj.org/article/812606d399224cb4bade3768754cc69b> | Check the data availability statement for post-2012 Indian Ocean points |
| **US JGOFS Arabian Sea 1995** | <https://www.bco-dmo.org/dataset/2524/description> | Historical baseline only (pre-VIIRS) |

**First task is a data inventory:** count what is actually downloadable for the box before any modelling. That decides how strong RQ1 and RQ4 can be.

### 5.4 Tools and starter code

| What | Where |
| --- | --- |
| `earthaccess` (NASA downloads in Python) | <https://earthaccess.readthedocs.io> |
| NASA OB.DAAC tutorial notebooks | <https://oceancolor.gsfc.nasa.gov> → Resources → tutorials |
| EPFT-UP (independent MOANA reimplementation and held-out tests) | <https://github.com/ocean-colour/EPFT-UP> |

## 6. Method and pipeline

1. **Ingest and preprocess.** One config (box 0–25°N, 50–100°E, dates, products) → download → subset → regrid to the 4.6 km grid → quality flags (cloud, sun glint) → mask ~25 km of coast → analysis cube (Zarr).
2. **Run MOANA ourselves.** Standardise Rrs 414–660 nm → PCA → multilinear regression per group → log10 cells/mL; log10 SST added for *Prochlorococcus* only. Verify against NASA's Atlantic files on the same dates before trusting it.
3. **Trust gate (RQ1).** Match MOANA to 2024+ ship counts (same day, ±1 day fallback, median of a 3×3 block). Metrics: log-RMSE, bias, Spearman ρ, share within a factor of 2, with a bootstrap that resamples whole cruises. Decision rule written **before** results are seen. Fallback: if a group has fewer than ~30 matchups, compare spatial and seasonal patterns with published cruise maps instead of absolute errors.
4. **Virtual VIIRS (RQ2).** Weight PACE spectra by SNPP and NOAA-20 band responses; compare with real VIIRS on same-day pixels and document the residual. Train identical models on full PACE vs virtual VIIRS; the skill drop is the cost of broad bands.
5. **Routes on identical splits (RQ3).**
   - **A (direct):** VIIRS bands and ratios + OISST + PAR + mixed-layer depth + lat/lon + day of year → log10 abundance. Models: XGBoost, LightGBM, small MLP. Train on virtual VIIRS, fine-tune on real VIIRS matchups.
   - **B (two-step):** small MLP rebuilds a PACE-like spectrum from the VIIRS bands, then MOANA runs on it.
   - **C (niche):** SST + PAR (+ MLD) only, Flombaum-style.
   - **Baselines:** climatology and ridge regression.
6. **Uncertainty and explanation.** Quantile regression (5/50/95%) for trees; ensemble of 5 for the MLP; out-of-range flags where inputs fall outside the PACE-era range; SHAP for interpretation (explains the model, not the ocean).
7. **Reconstruct 2012–2026** with the winning route as 8-day maps with intervals; SNPP and NOAA-20 reconstructions built independently.
8. **Climate analysis (RQ5).** Marine heatwaves by the Hobday definition (at least 5 days above the 90th percentile, OISST 1991–2020); composites (heatwave periods, 2019 positive IOD vs matched normal periods); lagged correlation with DMI and ONI; Mann–Kendall and Sen's slope (descriptive only); FDR correction.

Models are deliberately classical: ~2.5 years of PACE labels is too few for deep spatial models. U-Net and deep ensembles are stretch goals.

## 7. Validation and success criteria

| Test | Pass condition |
| --- | --- |
| Spatial-block CV (2°×2° blocks, size checked with a variogram) | Positive skill over climatology for every usable group |
| Temporal hold-out (latest 6 months of PACE) | Positive skill over climatology |
| Leakage control | Scalers, PCA bases and climatologies fitted inside training folds only |
| Skill definition | Skill = 1 − RMSE_model / RMSE_climatology |
| Interval calibration | 90% intervals cover 85–95% of held-out truth |
| SNPP vs NOAA-20 (2018–2026) | ≥80% of open-ocean pixels agree within stated uncertainty |
| Ship counts (2012–2023) | Reported per group, region and season, whatever the result |
| Route comparison | A vs B vs C reported in full even if A loses |

A failed trust gate or a null climate result is still a reportable finding.

## 8. Scope and non-goals

- **In scope:** surface open-ocean waters, 0–25°N, 50–100°E, Phase 1 regions first, 2012–Oct 2026 (SNPP end; NOAA-20 continues after).
- **Phase 2 (only if ahead of schedule):** western Arabian Sea upwelling zone and equatorial band.
- **Non-goals:** coastal and river-influenced waters, subsurface plankton, fisheries forecasts, causal claims, operational products.
- **Stretch:** U-Net on spatial patches; MODIS-Aqua extension to 2002 (orbit drift makes this uncertain).

## 9. Timeline: wrap by 31 December 2026 (work every day)

| Week | Dates | Work | Gate |
| --- | --- | --- | --- |
| 1 | 10–16 Nov | Earthdata and Copernicus accounts; **open-archive ship-count inventory**; test data access; logged Scopus/WoS/Shodhganga search (T01); read the MOANA ATBD and EPFT-UP; one-month test cube; draft Introduction and Related Work | Ship-count availability known; monsoon coverage acceptable |
| 2 | 17–23 Nov | Full downloads (PACE Rrs, VIIRS SNPP and NOAA-20, OISST, CMC SST, GLORYS12, DMI, ONI); preprocessing, coast mask, 4.6 km grid; Zarr cube; draft Data and Methods | Clean analysis cube |
| 3 | 24–30 Nov | Own MOANA implementation and Atlantic check; Bay of Bengal matchups; **trust gate** | **RQ1 verdict per group, fixed before any modelling** |
| 4 | 1–7 Dec | Virtual VIIRS and residual check; RQ2 information-content test; baselines (climatology, ridge); Route C; Route A (trees) | Route A beats climatology on held-out blocks |
| 5 | 8–14 Dec | Route A (MLP, fine-tuning); Route B; quantile intervals; SHAP; Arabian Sea expansion | All three routes scored on identical splits |
| 6 | 15–21 Dec | Full validation suite (spatial blocks, temporal hold-out, coverage); 2012–2026 reconstruction; SNPP vs NOAA-20; 2012–2023 ship-count check; out-of-range flags | RQ3 and RQ4 complete |
| 7 | 22–31 Dec | Climate analysis (heatwaves, IOD, ENSO); finish paper draft; clean GitHub repo; dataset release (Zenodo DOI); basic web map | **All project work wrapped by 31 Dec** |

Writing runs alongside the work: each results section is drafted as soon as its RQ finishes.

**After December (SRMIST calendar):**

| Dates | Event | Use |
| --- | --- | --- |
| 6–13 Jan 2027 | Title-change window | Confirm final title with guide, now that results are known |
| 23–30 Jan | Review I | Present finished Phase 1 pilot and trust-gate results |
| 2–9 Mar | Review II | Present information-content and route-comparison results |
| 20–24 Apr | Review III | Present validated record and climate analysis (TRL 4) |
| By 24 Apr | **Paper submitted** | January to April used for revision and submission |
| 28 Apr–6 May | Final report | — |
| 7–15 May | Viva | — |

IGARSS (deadline 10 Jan) is not targeted. Target journals from the department list: *Remote Sensing* (MDPI), *Remote Sensing Applications: Society and Environment*, *Earth Science Informatics*. Confirm whether *Frontiers in Marine Science* is on the list before choosing it.

## 10. Risks

| Risk | Likelihood | Mitigation |
| --- | --- | --- |
| Too few ship counts from open sources | Medium–High | Data inventory in week 1; pattern-comparison fallback; state limits openly; no author requests |
| MOANA fails the gate (cyanobacteria especially) | High | Publish the evaluation itself; use as relative index or exclude; picoeukaryotes as headline |
| Our MOANA code differs from NASA's | Medium | Check against NASA Atlantic files on the same dates |
| Seven-week schedule slips | Medium–High | Gates in the timeline; Phase 2 and stretch dropped first; Route B kept simple; writing alongside results |
| Large downloads and compute | Medium | Subset at download; Zarr chunks; trees run on a laptop; start downloads in week 1 to 2 |
| Similar paper appears (a NASA validation team collects Arabian Sea data) | Medium | Alerts; Bay of Bengal emphasis, three-route benchmark and dual-sensor test still differ |
| Monsoon cloud gaps | High | 8-day composites; report valid-pixel counts; no silent gap-filling; SW monsoon reported separately |
| Short PACE record misses extremes (e.g. 2019 IOD) | High | Out-of-range flags; skill reported separately for in-range and out-of-range inputs |
| Overclaiming causation | Medium | Wording rules: associations only; matched-period comparisons |

## 11. Lessons from the minor project (BlueEco) built into the design

| BlueEco limitation (from its own limitations section) | PicoPast response |
| --- | --- |
| Forecaster had no held-out test | Spatial-block CV and temporal hold-out |
| Scaler fitted on the full series | Scalers, PCA and climatologies fitted inside training folds only |
| Synthetic CO₂ input | No synthetic inputs; observed or reanalysis data only |
| Hand-picked thresholds | Pre-registered trust-gate rule |
| No uncertainty | Calibrated 90% intervals, coverage checked |

We reuse our experience with NASA ocean-colour data and gradient-boosted models, not results. Never describe the minor paper as "published" unless it is.

## 12. Deliverables

- Research paper draft complete by 31 Dec 2026; submitted by 24 Apr 2027
- Review Report, literature survey, Gantt chart (done for Zeroth Review; Gantt to be updated to this timeline)
- PicoPast v1.0 dataset (8-day maps, uncertainty and flags) with a Zenodo DOI
- Reproducible GitHub repository
- Interactive year-slider web map (basic version if time is short)

## 13. Open items to confirm

- BlueEco paper venue and submission status
- Course registration code (21CSP401L or 21CSP402L)
- Gantt task owners (Srikar / Aryan split)
- Guide's recommendation on treating this as an extension of the minor project
- Whether to keep Phase 2 regions as "conditional" (recommended) or drop them
- Whether *Frontiers in Marine Science* is on the department journal list

## References

1. [SA-CGAN: From Aqua MODIS to PACE OCI](https://spj.science.org/doi/10.34133/remotesensing.1036), Journal of Remote Sensing, 2026
2. [Global PFTs using PACE OCI and deep ensemble learning](https://www.sciencedirect.com/science/article/pii/S156984322600141X), 2026
3. [MOANA ATBD](https://oceandata.sci.gsfc.nasa.gov/files/atbd/atbd-obdaac-phytoplankton-groups-MOANA.pdf), NASA
4. [EPFT-UP](https://github.com/ocean-colour/EPFT-UP), GitHub, 2026
5. [Copernicus long-term PFT observations](https://sp.copernicus.org/articles/6-osr9/7/2025/), 2025
6. [NASA Arabian Sea validation data](https://www.earthdata.nasa.gov/data/catalog/ob-daac-pvst-northern-indian-ocean-0), 2025
7. [Flombaum et al., PNAS 2013](https://www.pnas.org/doi/10.1073/pnas.1307701110)
8. [Nazirahmed et al., JGR Biogeosciences 2026](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2025JG009023)
9. [Bay of Bengal picophytoplankton, J. Earth Syst. Sci. 2026](https://link.springer.com/article/10.1007/s12040-026-02759-7)
10. [MAREDAT picophytoplankton database, ESSD 2012](https://essd.copernicus.org/articles/4/37/2012/)
11. [Visintini et al., L&O Letters 2021](https://doaj.org/article/812606d399224cb4bade3768754cc69b)
12. [US JGOFS Arabian Sea picoplankton, BCO-DMO](https://www.bco-dmo.org/dataset/2524/description)
13. [OBPG relative spectral response tables](https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables)
14. [PACE data overview](https://pace.oceansciences.org/data.htm)
15. Lange et al., 2020, *Optics Express* (MOANA algorithm)
16. Hobday et al., 2016, *Progress in Oceanography* (marine heatwave definition)