# When cyclones hit heatwaves — rules fixed before testing

Status: **binding from 2026-11-22**, committed before any storm's cooling,
heatwave ending, bloom or plankton response is computed. Changes go in the
change log at the end, with the reason, and every earlier run stays reported.
Proposal: the "Project Proposal — When Cyclones Hit Heatwaves" doc; literature
check: `docs/literature_survey.md`.

## What has been seen already (disclosed)

- Storm counts from IBTrACS (22 November 2026): 81 storms reached 34 kt inside
  the box since 2012 with IMD winds where present and JTWC winds otherwise; 84
  with JTWC winds alone; 61 with IMD winds alone.
- A storm-level exposure count with a generous rule (any 0.25° cell on the
  track with heatwave conditions on at least half the days of the 8-day period
  before arrival): 42 storms over heatwave water, 32 over normal water.
- PicoPast's earlier climate work found lower picoeukaryote and *Synechococcus*
  indices during marine heatwaves (`climate_maps.md`,
  `route_d_reconstruction.md`). These are not cyclone responses.
- No cooling, heatwave ending, chlorophyll or PACE response has been computed
  for any storm in this project.

## Data

| Role | Product | Notes |
| --- | --- | --- |
| Tracks | IBTrACS v04r01, North Indian list, downloaded 22 November 2026 | 6-hourly fixes |
| Sea surface temperature, heatwaves | NOAA OISST v2.1, 0.25° daily, box files 1991–2026 | Heatwaves by `picopast.mhw`: Hobday et al. (2016), 90th percentile of 1991–2020, at least 5 days, gaps of up to 2 days bridged |
| Chlorophyll | VIIRS SNPP L3m 8-day chlor_a, 4 km, 2012–2026 (as on disk) | Pixel 8-day climatology from all years on disk |
| Sea level anomaly | NASA-SSH Simple Grid v1.1, 0.5°, about weekly (`ssha`) | Eddy covariate |
| Sea surface salinity | SMAP JPL CAP v5, 8-day running mean, 0.25° | Secondary covariate, storms from May 2015 |
| Phytoplankton groups | PACE OCI MOANA, 8-day, 4 km (operational mapping, as on disk) | RQ3 only |
| Coast | PicoPast offshore mask (pixels at least 25 km from land) | Chlorophyll and PACE pixels |

## Storms and track points (the unit of analysis)

- **Points:** IBTrACS fixes at 00, 06, 12 and 18 UTC with a valid position,
  from 1 January 2012, inside 0–25°N, 50–100°E.
- **Strength:** JTWC 1-minute wind (`USA_WIND`) of 34 kt or more. IMD 3-minute
  wind (`WMO_WIND`) of 34 kt or more is a sensitivity check.
- **Over the sea:** `DIST2LAND` of 50 km or more.
- **Day 0 (D0):** the fix's UTC date.
- **Disc:** every grid cell or pixel whose centre lies within **150 km**
  (great circle) of the fix. Sensitivity: 75 km and 300 km. Discs are cut at
  the box edge.
- **Basin:** Arabian Sea if the fix is west of 77.5°E, Bay of Bengal otherwise.
- **Season:** pre-monsoon (March–June), post-monsoon (September–December),
  other (January, February, July, August), by D0.
- **Translation speed (km/h):** great-circle distance between the storm's
  previous and next fixes divided by the time between them (one-sided at the
  ends).
- **Uncertainty:** always from resampling whole storms (2,000 draws, seed 0), so
  points from one storm never count as independent.

## Heatwave exposure (before the storm)

- **E:** the share of the disc's OISST ocean cells in a heatwave on **D0 − 3**.
  Three days keeps the storm's own cooling out. Sensitivity: D0 − 5.
- **Classes:**
  - *Heatwave point:* E ≥ 0.5.
  - *Normal point:* E < 0.1.
  - Points in between are left out of class comparisons and kept in the
    continuous analyses.
- **Also recorded on D0 − 3, as disc means:**
  - sea surface temperature;
  - its anomaly from the 1991–2020 daily climatology (SSTA);
  - heatwave excess (temperature minus the 90th-percentile threshold, over
    heatwave cells, 0 elsewhere).

## Responses (computed only after this file is committed)

**Cooling (RQ1).** The minimum disc-mean OISST over D0 to D0 + 5, minus the
mean disc-mean OISST over D0 − 7 to D0 − 3, in °C. Negative values are cooling.

**Heatwave survival (RQ1).**
- *Storm:* for heatwave points, the share of the disc's cells in a heatwave on
  D0 − 3 that are still in one on D0 + 5.
- *Control:* on the same days, the same share for every heatwave cell in the
  same basin lying more than 750 km from every fix of every storm active
  between D0 − 10 and D0 + 10.
- *Effect:* the mean over heatwave points of (storm share − control share).

**Re-warming (RQ1, descriptive).** Days from D0 until the disc-mean SSTA climbs
back to within 0.2 °C of its D0 − 3 value, searched from the day of minimum
temperature and censored at 60 days. Compared by Kaplan–Meier medians; no
verdict.

**Bloom (RQ2).**
- *Pixel anomaly:* log10 chlorophyll minus that pixel's mean for the same 8-day
  slot across 2012–2026. Offshore pixels only.
- *Periods:*
  - pre = the latest 8-day period ending on or before D0 − 2;
  - post = the first two 8-day periods starting on or after D0 + 1.
- *Disc value per period:* the mean anomaly over valid pixels. A period counts
  only if at least 20% of the disc's offshore pixels are valid.
- *Bloom B:* the mean of the counted post periods minus the pre period, in
  log10 units (+0.1 = 26% more chlorophyll). It needs a counted pre period and
  at least one counted post period.
- *Sensitivity:* the larger post period instead of the mean.

**Plankton groups (RQ3, PACE storms only).** The same pre and post periods,
applied to MOANA picoeukaryotes, *Synechococcus* and *Prochlorococcus*:
log10 cells, flags = 0, offshore pixels, disc means. *Synechococcus* and
*Prochlorococcus* from MOANA failed PicoPast's check C in this box, so they
are shown as exploratory.

**Covariates.**
- Wind (JTWC, kt) and translation speed.
- Latitude.
- Basin and season.
- `DIST2LAND`.
- The pre-period chlorophyll anomaly.
- Sea level: the disc-mean `ssha` of the NASA-SSH grid nearest in time on or
  before D0 − 3, in cm.
- Secondary: SMAP salinity, the disc mean of the file for D0 − 3 (fresher water
  marks the Bay of Bengal's barrier layer).

## Tests and decisions

**RQ2, the primary test.**
- *Model:* ordinary least squares of B on the heatwave class (heatwave = 1,
  normal = 0), wind, translation speed, latitude, basin, season, pre-period
  anomaly, sea level and `DIST2LAND`.
- *Effect:* the heatwave coefficient, with a 95% interval from resampling
  storms.
- *Minimum data:* at least 15 storms contributing heatwave points and 15
  contributing normal points, all with a valid B. With fewer, the result is
  *not enough data*, reported descriptively.
- *Decision:*
  - interval entirely below 0 → *weaker bloom over heatwaves* (hypothesis A);
  - entirely above 0 → *stronger bloom over heatwaves* (hypothesis B);
  - otherwise → *no detectable difference*, reported with the interval as the
    range of effects the data allow.
- *Also reported:*
  - the raw difference;
  - basin × season strata;
  - storm-equal weights (each point weighted 1 / that storm's points);
  - the continuous E in place of the class, on all points;
  - the 75 and 300 km discs, D0 − 5, IMD winds and the larger post period;
  - the SMAP model on storms from May 2015, with salinity added.

**RQ1, key secondary tests,** reported with 95% intervals:
- the cooling model, the RQ2 model with cooling as the response;
- the heatwave-survival effect, expected below 0 if storms end heatwaves;
- re-warming, descriptive only.

**RQ3, descriptive.**
- For each PACE-era storm: the median change per group, and whether
  picoeukaryotes rose relative to *Prochlorococcus*.
- No verdict.

**RQ4, prediction.**
- *Targets:* cooling and B at each point.
- *Inputs:*
  - wind, translation speed, latitude, basin, season (sine and cosine of the
    day of year);
  - `DIST2LAND`;
  - the pre-storm state: temperature, SSTA, E, heatwave excess, pre-period
    anomaly, sea level.
- *Models:*
  - a baseline (mean response by basin × season, from the training storms);
  - ridge regression on standardised inputs, alpha 1;
  - LightGBM with fixed settings: 200 trees, 15 leaves, learning rate 0.05,
    `min_child_samples` 20, `feature_fraction` 0.8, seed 0.
- *Evaluation:* leave one storm out. Skill = 1 − mean squared error / that of
  the baseline. Spearman ρ is also reported.
- *Heatwave value:* LightGBM skill with all inputs minus skill without SSTA, E
  and heatwave excess, with a 95% interval from resampling storms. The verdict
  "heatwave information adds skill" requires an interval above 0.
- *Which inputs matter:* SHAP values (LightGBM contributions) from the
  all-storm model, descriptive.
- No settings are tuned on these storms.

## RQ3 extended to all storms with PicoPast's own records (added 22 November 2026)

Added after the first results, whose RQ3 covered only the six PACE-era storms with
data, and before any PicoPast record value was computed at a storm point. Project
retitled at the same time: *How marine heatwaves change the phytoplankton response to
cyclones in the North Indian Ocean*.

**Records.** All are VIIRS SNPP, 8-day, 4 km, log10 cells/mL, using all valid offshore
pixels.

- *Primary: versions built without sea temperature,* so storm cooling cannot move them
  through the model's temperature input.
  - Picoeukaryotes: Route A retrained without SST. Same training table, features minus
    SST and LightGBM settings as `climate.sst_free_maps`, but all valid pixels instead
    of reliable-only, in log10 units.
  - *Synechococcus* and *Prochlorococcus*: Route D `syn_nosst` and `pro_nosst`.
- *Secondary:* the main records, Route A `peuk` and Route D `syn` and `pro`, which
  include SST.

**Response ΔG per group.**
- The same points, 150 km disc, pre and post periods and 20% validity rule as the
  bloom.
- Pixel anomalies come from each record's own 8-day slot climatology over 2012–2026,
  as on disk.
- ΔG = the mean of the counted post periods minus the pre period.

**Tests.** All are secondary, with storm-resampled 95% intervals. *Prochlorococcus*
results are exploratory, because its Route D record failed check E.

1. **Change after storms:** the mean ΔG per group over all points.
   - interval above 0 → "rises after storms";
   - below 0 → "falls";
   - otherwise → "no detectable change".
2. **Community shift:** the mean, per point, of Δpicoeukaryotes − Δ*Prochlorococcus* and
   of Δ*Synechococcus* − Δ*Prochlorococcus*, where both are valid. An interval above 0
   means the larger cells gain relative to *Prochlorococcus*.
3. **Heatwave effect per group:** the RQ2 model with ΔG as the response, and the
   group's own pre-period anomaly in place of the chlorophyll one. Decisions are worded
   as for RQ2 ("weaker" or "stronger response over heatwaves", "no detectable
   difference"), with at least 15 storms per class.
4. **Cross-check (descriptive):** for PACE-era storms with both, whether the PicoPast
   picoeukaryote change has the same sign as PACE MOANA's.

**Limit.** These records are model estimates built from chlorophyll, light and ocean
colour, plus SST in the secondary versions. Their changes are model-implied group
responses, not counts, and are not independent of the RQ2 chlorophyll response.

## Replication on 2002–2011 storms with MODIS-Aqua (added 22 November 2026)

**Why and when.** Decided after the first results. For 2012–2026 the bloom effect was
−0.058 [−0.125, +0.004], leaning weaker with the interval just reaching zero. Written
before any 2002–2011 response was computed. Adding storms after seeing a near-miss can
bias a study towards finding an effect, so the 2002–2011 storms are also judged on
their own (test 3).

**Storms and points.**
- IBTrACS fixes with D0 from 4 July 2002 (the start of MODIS-Aqua) to 31 December 2011.
- Every rule of the main study applies: JTWC ≥ 34 kt, ≥ 50 km from land, the box,
  6-hourly fixes, 150 km discs, exposure on D0 − 3 and the same classes.
- Points from 2012 on stay in the main study.

**Data.**
- *Chlorophyll:* MODIS-Aqua L3m 8-day chlor_a, 4 km (`MODISA_L3m_CHL`), cut to the
  box, files from 4 July 2002 to 31 January 2012.
  - The pixel 8-day slot climatology comes from those MODIS-Aqua files.
  - The bloom is computed exactly as for VIIRS: pre and post periods, the 20% rule.
- *Sea temperature and heatwaves:* the same OISST and `picopast.mhw` definition, with
  the 1991–2020 baseline; a daily cube from 1 November 2001 to 31 March 2012.
- *Sea level:* NASA-SSH Simple Grid v1.1 from October 2001.
- *Not available for these years:* salinity (SMAP starts in 2015), PACE, and PicoPast's
  group records (which start in 2012).

**Tests.**
1. **Replication (2002–2011 points only).** The RQ2 model, with storm-resampled 95%
   intervals and at least 15 storms per class (else *not enough data*).
   - interval excluding zero, with the same sign as 2012–2026 (negative) →
     *replicated*;
   - same sign but interval including zero → *same direction, not significant*;
   - otherwise → *opposite direction*.
2. **Pooled (2002–2026).** The RQ2 model on all points, adding an era indicator
   (MODIS 2002–2011 = 1) that absorbs any constant offset between the sensors.
3. **Updated RQ2 verdict.** *Weaker* (or *stronger*) *bloom over heatwaves* only if
   both hold:
   - the pooled interval excludes zero;
   - the 2002–2011 estimate alone has the same sign.

   Otherwise *no detectable difference* stands.
4. **RQ1 for 2002–2011 (secondary).** The cooling model and heatwave survival,
   reported alone and pooled with the era indicator.
5. **Pooled sensitivity checks.** Storm-equal weights, continuous E, and 75 km and
   300 km discs.

**Limits.**
- MODIS-Aqua and VIIRS chlorophyll differ in absolute level. The bloom is a change
  within one sensor's anomalies, and the era indicator absorbs a constant offset.
- 2002–2011 had fewer heatwaves (the warming trend), so its heatwave class may be
  small.

## Final confirmatory test: the whole Indian Ocean (added 27 November 2026)

**When.** Written after the main study and the 2002–2011 replication. The pooled North
Indian Ocean bloom effect was −0.052 [−0.114, +0.005]. This was written before any
South Indian Ocean value was analysed, and before any test in this section was computed
on any data. It is the last test of RQ2: its verdict is final, and no later data change
it. Methods and sources: `docs/methods_research.md`.

### Sample

- **North Indian Ocean:** the existing points, 4 July 2002 to 2026 (the main and
  replication tables).
- **South Indian Ocean.**
  - *Box and tracks:* 0–35°S, 30–120°E; IBTrACS v04r01 South Indian list (downloaded
    22 November 2026).
  - *Dates:* D0 from 4 July 2002 to the last date with data.
  - *Rules:* the same as the north: 6-hourly fixes, JTWC ≥ 34 kt, ≥ 50 km from land,
    150 km discs, exposure on D0 − 3, classes E ≥ 0.5 and E < 0.1.
  - *Data:*
    - OISST v2.1 for the box, 1991–2026, with heatwaves by the same `picopast.mhw`
      definition and 1991–2020 baseline;
    - VIIRS SNPP 8-day chlorophyll (2012 to September 2026) and MODIS-Aqua 8-day
      chlorophyll (July 2002 to January 2012), each with its own pixel 8-day slot
      climatology;
    - NASA-SSH, as in the north.
  - *Offshore pixels:* water is any pixel with valid chlorophyll in at least one 8-day
    period of that sensor's files; then the same 25 km coast rule.
  - *Seasons:* early (November–January), late (February–April), other (May–October).
- **Pooled model.** The RQ2 covariates, with these changes:
  - latitude enters as *absolute* latitude (distance from the equator; unchanged for
    the north);
  - an era indicator (MODIS);
  - indicators for the Bay of Bengal and the South Indian Ocean (the Arabian Sea is
    the reference);
  - an indicator for each region-season (north pre-monsoon is the reference).
- **Intervals:** resample whole storms, 5,000 draws, seed 0.

### T1, the primary and decisive test: does a heatwave change the bloom?

The pooled heatwave effect on the bloom, judged on its **99% interval**. This is the
third look at RQ2, so α ≈ 0.05/3, rounded to 0.01. The smallest effect that matters is
**±0.05 log10**: one third of the mean 2012–2026 north bloom (+0.15 log10), about a 12%
change in chlorophyll.

| 99% interval | Verdict |
| --- | --- |
| entirely below 0, and inside ±0.05 | *a small weakening, less than a third of the typical bloom* |
| entirely below 0 (otherwise) | *weaker bloom over heatwaves* |
| entirely above 0, and inside ±0.05 | *a small strengthening* |
| entirely above 0 (otherwise) | *stronger bloom over heatwaves* |
| entirely inside ±0.05, spanning 0 | *no meaningful effect* |
| anything else | *inconclusive* |

**Heterogeneity guard.** If the 95% interval of the heatwave × South Indian Ocean
interaction excludes zero, the pooled verdict is replaced by separate north and south
verdicts, by the same rule.

### Secondary tests T2–T4

These are decided at 95%, with a Holm correction across the three. p-values come from
the bootstrap draws.

**T2, mechanism: storm mixing power.**
- *Index:* M = log(wind² / translation speed), standardised over the pooled points. It
  is a proxy for the wind impulse a storm gives each patch of ocean.
- *Model:* T1's model plus M and heatwave × M.
- *Hypothesis:* heatwaves weaken blooms mainly under weak or fast storms, so the
  interaction is above 0.
- *Verdicts:*
  - interaction above 0 → *confirmed*;
  - below 0 → *the opposite: the heatwave effect grows with storm power*;
  - otherwise → *no detectable dependence on storm power*.
- *Also reported:* the heatwave effect at the 25th and 75th percentiles of M; the
  sensitivity M3 = log(wind³ / speed).

**T3, timing: mixed-up chlorophyll or new growth.** The heatwave effect is estimated
separately on two parts of the response:
- *the first post period* minus pre: the immediate signal, mostly chlorophyll mixed up
  from below;
- *the second post period* minus the first: later growth.

Each component has its own verdict.

**T4, storm-minus-no-storm control (difference-in-differences).**
- *Control years:* every other year of the same sensor era in which the same disc
  qualifies, judged by two conditions:
  - no IBTrACS fix of any strength, from either list, passed within 500 km of the disc
    centre between the pre period's start and the second post period's end;
  - the disc had the same heatwave class on the same calendar day (D0 − 3).
- *Control bloom:* for those years, the bloom computed on the same calendar 8-day
  periods. The control value is their mean.
- *Test:* T1's model on B − control.
- *Why:* this removes what heatwave and normal water do anyway over those weeks, such
  as heatwave decay and seasonal rebound.

### T5, power simulation (descriptive; does not change any rule)

Computed before T1, from the exposure counts of the full sample (exposure is not a
response):
- *Method:* resample the North Indian Ocean storms to the full sample's heatwave and
  normal storm counts, plant effects of 0 and −0.05, refit T1, and repeat 500 times.
- *Output:* the chance of each T1 outcome.

### Also reported

RQ1 (cooling, heatwave survival) for the South Indian Ocean alone and pooled, under the
original rules.

## Virtual experiments: the heatwave's effect with the storm held fixed (added 27 November 2026)

**What.** A 1-D column model (`picopast.column_model`: PWP mixing plus NPZD biology)
runs each selected storm twice:
- over the real pre-storm column;
- over the same column with the heatwave changed.

The storm, the place and the date are identical, so the paired difference is the
heatwave's effect *in the model world*. This answers the mechanism question that
observations cannot isolate. It does not replace T1, which stays the observational
verdict.

**When.** Written before any real storm was simulated; the model was checked only on
made-up profiles (`tests/test_column_model.py`).

**Points.**
- *Sample:* the final-test sample (north and south, 2002–2026).
- *Selection:* for each storm and each class (heatwave E ≥ 0.5, normal E < 0.1), the one
  point with the highest JTWC wind (ties: the earliest).
- *Rule:* heatwave points are the primary experiment; normal points are secondary.

**Columns.**
- *Profiles:* World Ocean Atlas 2023 (1991–2020 climate normal, 1°) temperature,
  salinity and nitrate for the month of D0, at the nearest ocean cell, on a 2 m grid to
  200 m.
- *Actual column:* the climatological surface layer (down to where temperature first
  drops 0.5 °C) is shifted by the observed pre-storm OISST minus the atlas surface value,
  tapering to zero over the next 10 m.
- *Counterfactuals:*
  - heatwave points (primary): the same shift minus the observed pre-storm anomaly
    (`ssta_pre`), so the column has climatological temperature (heatwave removed);
  - normal points (secondary): plus the median `ssta_pre` of the heatwave points in the
    same region (heatwave added).
- *Biology:*
  - starting values: nitrate from the atlas; phytoplankton 0.1 mmol N m⁻³ in the top
    100 m and 0.01 below; zooplankton and detritus 0.05;
  - spin-up: 30 days under a light 5 m s⁻¹ rotating wind, with surface PAR fixed at
    100 W m⁻²;
  - each version spins up separately, so the heatwave also acts on growth through
    temperature.

**Storm forcing.**
- *Wind:* the storm's own IBTrACS fixes (both lists), interpolated hourly from D0 − 2 to
  D0 + 2. A modified Rankine vortex from the JTWC wind and radius of maximum wind
  (40 km where missing), with 20° inflow and half the translation velocity added.
  Clockwise rotation in the south.
- *Drag and steps:* Large & Pond drag capped at 2.4 × 10⁻³; 900 s steps during the
  storm, then 17 days of light wind at 3,600 s steps.

**Outputs per run.**
- *cooling:* minimum surface temperature over D0 to D0 + 5, minus the value at D0 − 2;
- *bloom:* log10 of the mean surface chlorophyll (first optical depth) over D0 + 1 to
  D0 + 16, minus log10 of the value at D0 − 2;
- *depth-integrated chlorophyll change* to 200 m (near zero means redistribution, not
  growth);
- *maximum mixed-layer depth;*
- *surface nitrate increase.*

**Tests.** For each output, the mean paired difference (heatwave − no heatwave), with a
95% storm-resampled interval (5,000 draws).
- *Verdict:* for the bloom, and for depth-integrated chlorophyll:
  - interval excluding zero → "in the model, a heatwave weakens (or strengthens) the
    bloom";
  - otherwise → "no detectable model effect".
- *Validation (descriptive):* Spearman ρ of model against observed cooling, and model
  against observed bloom, over the actual runs.
- *Sensitivity checks:*
  - anomaly confined to the top 20 m;
  - surface PAR 60 W m⁻²;
  - drag cap 2.0 × 10⁻³;
  - vortex decay exponent 0.5.

**Limits.** One column has no horizontal advection or eddies. The heatwave enters as a
warmer surface layer, not as a separately observed shallower mixed layer. The biology
uses fixed literature parameters, not ones fitted to this region.

## Known limits, stated with every result

- OISST has no microwave input and under-reads cold wakes beneath cloud. The
  bias is the same for both classes.
- 8-day composites smear the bloom's timing, and cloud removes some storms.
  Both counts are reported.
- Mixed-layer depth after 2017 is not available here. Sea level and salinity
  stand in for it.
- Heatwave and normal points differ in season and basin. The model and the
  strata address this but cannot remove every difference.

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-11-22 | First version, binding | Written before any storm response was computed. |
| 2026-11-22 | Added a post-hoc, non-binding check of cooling against far-away water with the same heatwave status on the same days (per-cell cooling, same basin, more than 750 km from any storm) | Added after the first results: heatwave water cools on its own as heatwaves decay, so part of the pre-registered cooling difference may not be the storm's. Reported beside the pre-registered result, which it does not replace. |
| 2026-11-22 | Added RQ3 extended to all storms with PicoPast's own records (section above); project retitled to put phytoplankton first | Asked by the project lead after the first results; written before any PicoPast record value was computed at a storm point. |
| 2026-11-22 | Added the 2002–2011 MODIS-Aqua replication (section above) | Decided after the first results to sharpen RQ2; written before any 2002–2011 response was computed, with a both-must-agree rule against optional stopping. |
| 2026-11-27 | Added the final confirmatory test (section above): South Indian Ocean storms, a three-outcome rule with a smallest meaningful effect and a 99% interval, and secondary tests T2–T4 with a power simulation | Asked for a definitive answer after the replication; written before any South Indian Ocean value or new test was computed. Methods from `docs/methods_research.md`. |
| 2026-11-27 | Added the virtual experiments (section above) | The simulation part of the methods plan; written before any real storm was simulated. |
| 2026-11-27 | Bug fix in the column model, then every virtual-experiment variant rerun from scratch. The gradient-Richardson partial mix is now capped at full mixing of the two layers. The report also counts runs with non-finite outputs. | *Found:* after the first run (base and three sensitivity variants). About 4% of runs had non-physical values (NaN, a surface warming of 580 °C, nitrate of 10⁶⁸). *Cause:* below a density inversion the Richardson number is negative and the uncapped mixing fraction overshoots. *What was seen:* the base paired summary was looked at before the fix. *Why that does not bias the result:* the fix follows from the numerical failure alone; it changes no rule, input or output definition; it adds a unit test. The first run's files are kept outside the repository (`data/interim/virtual/bugged_20261127/`) and are not reported as results. |
