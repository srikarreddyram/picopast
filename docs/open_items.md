# Open items

Updated 16 November 2026 (week 2), evening. Ordered by how soon each one blocks work.

## Done since week 1 (16 November)

- **Earthdata login** stored in `~/.netrc` on the main machine.
- **Atlantic MOANA check passed** (`picopast moana-check 2025-01-15`, operational
  mapping). Same 32,550 valid pixels as NASA's L4 file; median difference
  0.8% (Prochlorococcus), 0.4% (Synechococcus), 0.3% (picoeukaryotes); 95% of
  pixels within 2.7%, 3.6% and 2.3%.
- **OISST complete**, 1991-01-01 to 2026-11-14, no missing days. The file
  download now resumes after a dropped connection (HTTP Range).
- **Fast PACE box reader.** A PACE AOP file is ~840 MB; our box is ~20 MB of
  it. The reader walks the HDF5 chunk index in parallel and fetches only the
  box's chunks: ~10x faster than reading through xarray, values identical,
  and files stored as NASA's int16 (half the size). Tests in
  `tests/test_ingest.py`.
- **Downloads queued** (box subsets): January 2025 for PACE, CMC SST and all
  VIIRS products; then PACE 8-day and CMC SST and PACE daily for the whole PACE
  era.

## Blocking now

1. **Copernicus Marine account** for mixed-layer depth (week 2).
2. **SeaBASS files: first pass done (16 November).** The three regional
   collections that overlap the PACE era (PVST_NORTHERN_INDIAN_OCEAN,
   EKAMSAT_Pilot_ASTRAL, OCEANX_2025) hold **no flow cytometry or cell
   counts**. The two "candidate" files (`*_CLS.sb`) are flow-through
   fluorometer records, though their phycoerythrin channels (Fpe565/578/590)
   track Synechococcus-type pigment. They **do** hold HPLC pigments, including
   divinyl chlorophyll a (unique to Prochlorococcus) and zeaxanthin
   (cyanobacteria). Surface (<= 10 m) samples inside the box:
   - PACE era: 95 samples, roughly 27 station-days on 16 days. EKAMSAT Bay of
     Bengal, 28 April to 11 May 2024 (two releases, mostly distinct sample IDs),
     plus OceanX 9-10 January 2025 near the equator.
   - VIIRS era: 43 samples, EKAMSAT pilot, Arabian Sea, 11-24 June 2023.
   Files are in `data/insitu/seabass/`. ASIRI (2013) and the older collections
   are not yet opened.

## Decisions (made 16 November by Ramachandra; Aryan and the guide did not review)

3. **Trust-gate rule: fixed.** `docs/trust_gate_preregistration.md` is
   binding. Thresholds were accepted unchanged, and fallback checks were
   added: a DVChl a pigment proxy for Prochlorococcus, phycoerythrin for
   Synechococcus (exploratory only), and a same-month climatological pattern
   check against the 1991–2005 counts.
4. **Region boxes: fixed** as drafted. Bay of Bengal 5–22°N, 80–95°E;
   central/eastern Arabian Sea 8–22°N, 62–75°E.
5. **Ship-count shortfall.** See the inventory report. If the 2024+ count
   stays below about 30 matchups per group after SeaBASS files and published
   tables are read, the PRD's fallback applies (pattern comparison instead of
   absolute errors) and RQ1 and RQ4 must be worded accordingly. Decide early
   whether that changes the paper's headline.
   **Decided 16 November:** the open archives give no PACE-era cell counts in
   the box, so the fallback applies to all three groups. The paper frames
   MOANA as a **relative index checked against patterns**, not as validated
   abundances. The DVChl a pigment proxy (Prochlorococcus) and the
   climatological pattern check are written into the pre-registration. If
   the published tables (item 6) later yield 2024+ cell counts, they are
   added under the main rule, and the change log records it.
6. **Published tables.** The two 2026 papers named in the PRD could not be
   read from here (publisher pages refuse automated access or need a login).
   Someone with library access needs to check each for station-level tables
   and note the cruise dates. Record every table in a registry file with its
   source, page and who typed it.

## Corrections to the PRD found while building

7. **Prochlorococcus is regressed in linear cells/mL, not log10** (PRD section
   6 step 2 says log10 for all groups). Negative retrievals are possible and
   are flagged. Details in `docs/moana_notes.md`.
8. **The ATBD equations and NASA's operational coefficients disagree** for two
   terms. We default to the operational table and can run both.
9. **PACE Rrs collection name.** In version 3.2 the hyperspectral Rrs is in
   `PACE_OCI_L3M_AOP`; there is no `PACE_OCI_L3M_RRS` collection in CMR.
10. **The RSR table link in the PRD is dead; replaced (17 November).**
    `oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables` redirects to a
    general NASA page. The band-averaged M1–M5 responses now come from NOAA's
    calibration centre (`picopast.ingest.rsr`, cached in `data/raw/rsr`):
    SNPP from Northrop Grumman's October 2011 fused release (used in IDPS),
    NOAA-20 from the DAWG at-launch public release V2.1. Response-weighted
    centres: SNPP 412.6/444.4/487.2/550.8/670.3 nm, NOAA-20
    411.4/445.0/488.9/556.6/667.3 nm, matching NASA's band labels (SNPP M1 is
    labelled 410 but centres near 412.6); widths 15–20 nm.
11. **DMI lags.** NOAA PSL's HadISST-based Dipole Mode Index currently ends in
    May 2026. The last months of the record will have no DMI unless another
    source is chosen.
12. **Download volume.** A global PACE AOP daily 4 km file is 0.5–0.8 GB, so
    the PACE era is about 970 files and several hundred GB if kept whole. The
    pipeline reads only the box from each file over HTTPS; how fast that is
    will be known after the first logged-in run.

13. **OISST server limits.** NOAA PSL's OPeNDAP service started answering
    HTTP 429 (too many requests) after two years of box-only reads. The
    download now stops cleanly on 429, and `picopast fetch-open --oisst-method
    file` fetches each global yearly file (about 0.5 GB), keeps the box and
    deletes the rest. Both methods give identical grids.

14. **Download environment.** The main machine sleeps on battery and freezes
    transfers, so long downloads run under `caffeinate -i`. The line gives
    about 3 MB/s in total, so time, not disk, limits the download plan: full
    daily PACE for the box is ~18 GB but takes about a day (~85 s per file).

15. **MOANA Prochlorococcus outside the Atlantic (first look, January 2025,
    satellite only).** Under the operational mapping, 12.6% of retrieved
    pixel-days in the box give a negative Prochlorococcus; under ATBD, 11.3%.
    The negatives cluster near land: 57% of retrievals in water within 25 km
    of land, 16% at 25–200 km, 6% beyond 200 km. Positive values are higher
    than NASA's Atlantic output and nearly flat across regions (median
    400,000–490,000 cells/mL; NASA Atlantic 2025-01-15: 274,000). Synechococcus
    and picoeukaryotes show the expected regional contrasts, for example highest
    Synechococcus in the western Arabian Sea upwelling. The linear
    Prochlorococcus regression, fitted to Atlantic data, appears to extrapolate
    here. This is consistent with NASA producing MOANA L4 only for the
    Atlantic, and it is a result to report, not a bug. Our code matches NASA's
    within 1% on the same inputs.
16. **Land mask.** CMC fills SST over land, so the coast mask uses CMC's
    `mask` flag (1 = water, 2 = land), not SST gaps. `cmc_sst` downloads now
    keep `mask`.

17. **Checks A and B: results (16 November, run once under the fixed rule).**
    See `docs/reports/trust_gate_fallback.md`.
    - **A (Prochlorococcus vs HPLC DVChl a): not supported** under both
      mappings. Only 12 of 27 station-days matched; the other 15 had fewer
      than 5 cloud-free pixels on all three allowed days (checked one by one).
      ρ = −0.05 (operational) and −0.34 (ATBD), with intervals spanning about
      −0.7 to +0.6. Read with care: this is a low-power test, not evidence
      that MOANA is wrong. DVChl a ranged only 0.023–0.054 mg m⁻³ on the
      cruise, and the matches come from one cruise in one season. Two
      station-days (8–9 May, 2.8°N 83.3°E) share a pixel and a satellite
      value. By the pre-registered rule, Prochlorococcus is **not
      reconstructed**.
    - **B (Synechococcus vs phycoerythrin, exploratory): no relationship**,
      ρ ≈ 0 over 692 pixel-days, under both mappings.
    - **C (same-month PACE-era climatology vs 1991–2005 surface counts),
      operational mapping:** picoeukaryotes **pattern consistent**, ρ = 0.51
      [0.24, 0.64], 181 stations from 13 cruises; leaving out any one cruise
      gives 0.44–0.55. Synechococcus **inconclusive**, ρ = 0.18 [−0.09, 0.38],
      197 stations. Prochlorococcus **inconclusive**, ρ = 0.09 [−0.21, 0.36],
      179 stations. Unmatched stations are all July–September (monsoon cloud
      leaves the climatology too sparse there). Under ATBD (sensitivity only)
      Prochlorococcus is consistent (0.30 [0.06, 0.47]); Synechococcus and
      picoeukaryotes are as above.
    - **First run of C counted stations twice.** The BCO-DMO compilation repeats
      PANGAEA's 1995 JGOFS stations, and the per-sample de-duplication missed
      most of them because PANGAEA gives TT053/TT054 depths in whole metres.
      That run (385 stations, 4 invented cruise labels) gave Synechococcus and
      Prochlorococcus higher ρ (operational 0.18 and 0.12; ATBD both
      "consistent"). De-duplication is now per station (1,433 records removed),
      C was rerun once, and the change log records it. The picoeukaryote result
      did not change.
    - **Decision (operational, by the pre-registered rule):** only
      **picoeukaryotes** go to the reconstruction, as a standardised relative
      index (never cells/mL). Synechococcus and Prochlorococcus are evaluated and
      reported but not reconstructed.
    - **Caveat for the paper:** MOANA picoeukaryotes track chlorophyll, so the
      pattern match may partly be a chlorophyll pattern. The week-4 baselines
      should include a chlorophyll-only climatology in the same check C, to show
      how much MOANA adds.

18. **Back-extension, first attempt (17 November).** `picopast reconstruct`
    trains on the PACE-era overlap (VIIRS SNPP 8-day + OISST vs our MOANA,
    115 periods, 1.15 million sampled pixels) and back-extends **all three
    groups** to 2012. Each output carries its trust-gate tier: picoeukaryotes
    a relative index; Synechococcus and Prochlorococcus *exploratory* (team
    request after the gate; logged in the pre-registration change log).
    Report: `docs/reports/reconstruction.md`.
    - **Training screen.** Pixels follow the trust gate's rule (MOANA flags 0)
      and must sit inside NASA's declared valid range for its MOANA product
      (Pro 0–600,000; Syn 0–300,000; peuk 0–40,000 cells/mL). Without it, 7.8%
      of PACE-era pixels (5.4% with a negative Prochlorococcus) would enter
      training as MOANA extrapolations (Synechococcus down to
      1e-38 cells/mL, Prochlorococcus down to −7.5 million), and they inflated
      the scores: Route A spatial skill was 0.39 (Pro) and 0.38 (Syn) before
      the screen and 0.18 and 0.20 after. Picoeukaryotes barely changed (0.40).
    - **Validation, Route A vs seasonal climatology (screened):** skill 0.40
      spatial / 0.39 temporal for picoeukaryotes (ρ 0.91 / 0.83), 0.20 / 0.11
      for Synechococcus, 0.18 / 0.18 for Prochlorococcus; 90% intervals cover
      87–89%. Route C (SST + PAR) never beats climatology.
    - **Early-to-late shift, not yet a trend.** Study-box index, 2023–25 minus
      2012–14: picoeukaryotes −0.14 SD, Synechococcus +0.12, Prochlorococcus
      +0.17, same sign in every region. Item 21: only the picoeukaryote
      decline is also seen by NOAA-20.
21. **SNPP vs NOAA-20 (RQ4, 17 November).** NOAA-20 back-extended
    independently (own bands, own models; 2018–2026, 402 periods; validation
    as SNPP's: picoeukaryote skill 0.40 / 0.41). Report:
    `docs/reports/sensor_comparison.md`.
    - PRD agreement test (NOAA-20 median inside SNPP's 90% interval) passes
      for all groups (98–100%), but it is lenient: the intervals carry each
      model's error against MOANA. Typical gap between the two study-box
      means: 0.02 SD (picoeukaryotes), 0.04 (Synechococcus), 0.10 (Pro).
    - **Picoeukaryotes are sensor-independent:** the records stay within
      ±0.07 SD of each other over 2018–2026, and both show a decline from
      2018–20 to 2023–25 (SNPP −0.11, NOAA-20 −0.07 SD). A candidate signal
      for the RQ5 trend tests, still a relative index.
    - **Synechococcus and Prochlorococcus are not:** the two records agree in
      the training years (by construction) and drift apart going back, with
      a winter swing of up to 0.3 SD. Their older values depend on which
      sensor is used, and the SNPP-only rise in both is not seen in NOAA-20.
      Another reason they stay exploratory.
30. **Check E: Prochlorococcus fails the new pigment test (21 November;
    `docs/reports/check_e.md`, rules in `docs/check_e_preregistration.md`, committed
    before the test).**
    - **New data:** no new cell counts exist in PANGAEA, BCO-DMO, gVOD or SeaBASS.
      SeaBASS holds HPLC pigments: 149 surface station-days in the box from 12
      cruises (CLIVAR 2007, Tara 2010, GO-SHIP I07N 2018, EKAMSAT 2023 and 2024,
      GO-BGC 2023, Bio-GO-SHIP I09N 2025, SO-PACE TN444 2025, OceanX 2025).
    - **Result:** Route D's *Prochlorococcus* record vs divinyl chlorophyll a,
      same 8-day period, 2012 on: ρ −0.07 [−0.30, +0.16] (69 stations, 6
      cruises), *not consistent*; the same-month pattern −0.02. Chlorophyll alone
      and MOANA fail it too. **Binding decision: *Prochlorococcus* is exploratory.**
    - **Why (post hoc):** the pigment follows total chlorophyll (ρ 0.76); the
      1991–2005 cell counts fall where chlorophyll rises. Cells and pigment rank
      stations in opposite ways here, so the two pre-registered tests disagree.
    - **To settle it:** *Prochlorococcus* cell counts from 2012 on (NIO Goa and
      the EKAMSAT/BoBBLE groups are the likeliest holders). *Synechococcus* is
      untested by check E (zeaxanthin is shared; no new counts) and stays a
      relative index on check C.
    - Also fixed: the HPLC reader now reads three SeaBASS layouts, and rows
      without a sample ID are no longer merged (PACE-era count unchanged, 27).
29. **Route D record for the cyanobacteria, 2012–2026 (21 November;
    `docs/reports/route_d_reconstruction.md`, overview in `docs/reports/comparison.md`).**
    `picopast route-d-reconstruct`: D2 rebuilt *Synechococcus* and
    *Prochlorococcus* for all 674 eight-day periods (SNPP) and 402 (NOAA-20).
    The pre-registered decision; tier "relative index, Route D".
    - **Sensors agree:** r 0.97–1.00, drift −0.03 and +0.05 SD per decade, down
      from +0.13 and +0.25 in the old MOANA-based records. Same result with the
      no-SST models, which see only what differs between the satellites.
    - **The old records' significant 2012–2025 rises are gone:** they were the
      SNPP drift. The new and old records disagree on timing (monthly anomaly
      r −0.48 to +0.15).
    - **Synechococcus:** lower after El Niño and positive IOD (r −0.56 / −0.60)
      and during heatwaves (−0.56 on the maps), but most of that is the SST
      input. Without SST: r −0.41 / −0.34, heatwave −0.25. The Bay of Bengal
      trend (−0.24 per decade) disappears without SST.
    - **Prochlorococcus:** higher during heatwaves (+0.26 on the maps, 68% of
      cells clear; study box +0.15 [+0.06, +0.20]). This holds without SST,
      so it comes through low chlorophyll and high light. Weak ENSO and IOD
      links.
    - **Found and fixed:** Route D's reliable pixels thin out during heatwaves.
      *Prochlorococcus* training SST ends at 29.9 °C, so the warmest water drops
      out, and this had hidden the *Prochlorococcus* rise in the regional
      average. The report now gives all-valid-pixel results next to the
      reliable-pixel ones. The picoeukaryote record is not affected (its
      in-range share is ≈ 90% in both heatwave and normal periods).
    - **Limits:** 57% of *Prochlorococcus* pixels and 74% of *Synechococcus*
      pixels are inside the training range.
28. **Route D passes for the cyanobacteria (20 November; `docs/reports/route_d.md`,
    rules in `docs/route_d_preregistration.md`, committed before the test).**
    A LightGBM model trained on MAREDAT ship counts from the rest of the
    world ocean (5,445 surface stations, Indian Ocean excluded) on SST,
    chlorophyll, PAR, nitrate and season, with importance weights towards
    the box's conditions (D2), scored exactly like MOANA in check C:
    *Synechococcus* 0.50 [0.31, 0.64], *Prochlorococcus* 0.45 [0.25, 0.63],
    both *pattern consistent*; paired against MOANA +0.32 and +0.35 with
    intervals above zero. Picoeukaryotes 0.52, no better than MOANA (stays).
    Caveats: D2 is no better than chlorophyll alone (paired +0.14 and +0.01,
    intervals include zero), so the gain is in getting each group's map,
    not in new information; a few dozen warm-pool Pacific samples carry
    most of the weight, but leaving out any contributor or the whole warm
    pool keeps every group passing (lowest: *Synechococcus* 0.32 [0.13, 0.51]).
    Next: reconstruct both groups from D2 for 2012–2026 (pre-registered
    decision), with the SST-free repeat for any climate result.
27. **RQ5 per 0.25° cell (20 November; `docs/reports/climate_maps.md`).**
    `picopast climate-maps`: the reconstruction's reliable pixels averaged
    into OISST cells, composites per cell with year-resampled intervals.
    Picoeukaryotes: El Niño and positive IOD lower the index in the
    Arabian Sea (5–15°N: −0.32 and −0.44 SD, 83–93% of cells clear) while
    most of the Bay of Bengal does not move and its north rises slightly
    (+0.11): the regional average hid a west–east split. Heatwaves lower it
    almost everywhere (96% of cells clear; Bay of Bengal −0.25, Arabian Sea
    −0.38). Checks: a model without SST keeps 71–74% of each signal with
    the same pattern (r 0.98–0.99); NOAA-20 gives the same maps for
    2018–2026 (pattern r 0.96–0.99). Cyanobacteria maps are exploratory:
    their heatwave maps split north (higher) / south (lower).
26. **RQ3 complete: Route B (17 November; `docs/reports/route_b.md`).** An MLP
    rebuilds the PACE spectrum (400–680 nm, 12 principal components) from
    the five bands; MOANA runs on it. The rebuild is accurate (median error
    under 1% at 400–550 nm, 5–10% at 620–660 nm). Trained on virtual bands
    only, it fails on real VIIRS (Synechococcus skill −0.59 to −0.91): MOANA
    amplifies the small virtual-vs-real band differences. Trained on real
    VIIRS bands paired with PACE spectra (the PRD's intent) it beats
    climatology for every group (picoeukaryotes 0.28 / 0.32,
    Synechococcus 0.06 / 0.10, Prochlorococcus 0.06 / 0.11; spatial /
    temporal) but stays below Route A (0.39, 0.20 / 0.11, 0.18). Ranking on
    real VIIRS: **A > B > C**, C no better than climatology. The PRD's
    hypothesis (direct learning matches or beats the others) holds.
25. **RQ2 (17 November; `docs/reports/rq2_information.md`).** Identical models
    reproduce MOANA from the full PACE spectrum or from the same spectrum
    through SNPP's five bands (~590,000 pixels, 118 periods). LightGBM skill
    against a seasonal climatology: full 0.71–0.77, virtual VIIRS 0.56–0.64;
    the bands cost 0.13–0.18 for every group. In rank terms Synechococcus
    loses most (temporal ρ 0.96 → 0.83, others 0.98 → 0.93 and 0.96 → 0.92):
    the PRD's hypothesis holds for ranks, not for skill.
    **Most of the PACE-to-VIIRS loss is not the bands but the sensor
    mismatch** (overpass time, cloud sampling within the 8-day composite,
    calibration, atmospheric correction): real VIIRS skill is 0.39 for
    picoeukaryotes (mismatch costs 0.17) but 0.11–0.20 for the cyanobacteria
    (mismatch costs 0.39–0.49). Their MOANA signal lives in spectral detail
    that does not survive a change of sensor, consistent with items 17, 20
    and 21. Run: 80 min, 5.7 GB peak memory.
24. **Virtual VIIRS vs real VIIRS (17 November, January 2025, 29 days, ~1.3
    million same-day pixel pairs per sensor; `docs/reports/virtual_viirs.md`).**
    Blue-green bands within ±3% (r 0.87–0.97); the blue/green ratio behind
    chlorophyll: bias −0.006 (SNPP) and +0.026 (NOAA-20), r 0.93–0.94. Two
    outliers, documented, not corrected: the red band reads +16% to +34% in
    relative terms on tiny values (~0.0001 sr⁻¹), and **NOAA-20's 411 nm band
    reads about 17% below PACE and SNPP** (virtual − real +25%), an offset in
    NOAA-20's own L3m product. Each sensor has its own models, so a constant
    offset is learned; a drifting one would not be.
23. **RQ5, first pass (17 November).** `picopast climate`, SNPP record,
    reliable pixels only (`docs/reports/climate.md`). Picoeukaryotes:
    - **No trend** 2012–2025 (study box −0.02 SD/decade, p = 1.0).
    - **ENSO and IOD:** lower after El Niño and positive IOD (study box r =
      −0.59 with ONI at lag 0, −0.55 with DMI at lag 2; FDR p < 0.01;
      effective n 40–57). NOAA-20 alone over 2018–2026 gives the same
      (−0.53, −0.62). ONI and DMI co-vary, so the two cannot be separated.
    - **Marine heatwaves:** 0.25–0.36 SD lower during heatwave periods in
      every region, intervals clear of zero.
    - **The 2019 positive IOD on its own** is unremarkable (z −1.0 to +0.2).
    - **SST is a model input, so the associations were rechecked with a model
      trained without SST:** heatwave −0.18 to −0.25 SD, ONI r −0.52, DMI r
      −0.47 in the study box. About three quarters of each association
      survives on VIIRS optics alone.
    - Synechococcus and Prochlorococcus show significant upward trends, but no
      larger than SNPP's drift against NOAA-20 (+0.13 and +0.25 SD/decade):
      not evidence of change in the ocean.
    - Still to do: mixed-layer depth as a covariate (Copernicus account),
      pixel-level composite maps, monsoon-phase composites.
22. **PACE download complete (17 November, 12:52).** 922 daily and 118
    8-day box files: everything NASA's archive holds for 5 March 2024 to
    30 September 2026. NASA has no PACE L3M data yet for October 2026, so the
    PACE era ends on 30 September in practice (config still says 31 October).
    The 18 other missing days are archive gaps, the longest 25–29 May 2024.
20. **Does MOANA add to chlorophyll? (17 November, descriptive, after the
    gate.)** Check C rerun with a PACE-era VIIRS chlorophyll climatology in
    place of MOANA, same stations (`docs/reports/chlorophyll_comparison.md`).
    Picoeukaryotes: MOANA ρ 0.52, chlorophyll 0.47, difference +0.05
    [−0.08, +0.14]; partial ρ of MOANA given chlorophyll 0.25 [−0.01, 0.45], a
    hint of added information, not established. Synechococcus: chlorophyll
    does better (0.35 vs 0.19, difference −0.17 [−0.49, −0.03]).
    Prochlorococcus: low chlorophyll does far better (|ρ| 0.43 vs 0.07,
    difference −0.35 [−0.73, −0.11]); MOANA adds nothing (partial −0.01). For
    the paper: at the level of 1991–2005 patterns, MOANA's picoeukaryote index
    is about as good as chlorophyll and the cyanobacteria are worse. Say so.
19. **SNPP 8-day record in 2026 stops at 28 August** in NASA's archive (plus
    one period, 22–29 September). The training overlap therefore ends there;
    NOAA-20 has 34 periods in 2026.

## Not started (by design)

- Zarr cube builder: waits for the first real PACE and VIIRS files so it is
  written against actual file layouts.
- Routes A, B, C, baselines, intervals, SHAP (weeks 4–5).
- Reconstruction, dual-sensor test, climate statistics beyond heatwave
  detection (weeks 6–7).
- Scopus / Web of Science / Shodhganga search log (T01): needs institutional
  access.
