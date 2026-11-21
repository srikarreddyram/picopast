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
