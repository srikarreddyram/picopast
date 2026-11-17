# Back-extension to 2012: first attempt

Generated 2026-11-17 by `picopast reconstruct`. Models learn MOANA (operational mapping) from VIIRS SNPP 8-day inputs in the PACE era and are applied to every VIIRS 8-day period on disk. They reproduce MOANA, so they can be no more right than MOANA is.

## What each group's record is

| group | tier | route used | beats climatology on held-out blocks |
| --- | --- | --- | --- |
| Prochlorococcus | exploratory: did not pass the trust gate (check A not supported; check C inconclusive) | A | yes |
| Synechococcus | exploratory: did not pass the trust gate (check C inconclusive) | A | yes |
| picoeukaryotes | relative index: passed the trust gate (check C pattern consistent) | A | yes |

Only a group that passed the trust gate is used for headline results. The exploratory records are an attempt, shown with their bottlenecks below; they are not validated estimates.

## Training data

- PACE-era 8-day periods with both MOANA and VIIRS: 115 (05 Mar 2024 to 22 Sep 2026); up to 10,000 offshore pixels sampled per period.
- Pixels per group: Prochlorococcus 1,150,000, Synechococcus 1,150,000, picoeukaryotes 1,150,000.
- Units: Prochlorococcus 1e5 cells/mL (linear); the others log10 cells/mL.
- Route A inputs: rrs410, rrs443, rrs486, rrs551, rrs671, r443_551, r486_551, r410_443, r671_551, mbr, par, sst, lat, lon, doy_sin, doy_cos. Route C inputs: sst, par.

## Validation (out of sample)

Skill = 1 − RMSE_model / RMSE_climatology (> 0 beats the seasonal climatology). Spatial: 2° blocks in 5 folds, predictions pooled. Temporal: the last 6 months of the PACE era held out. Coverage: share of held-out values inside Route A's 5–95% interval (target 85–95%).

| group | split | climatology RMSE | ridge skill | Route C skill | Route A skill | Route A ρ | Route A 90% coverage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Prochlorococcus | spatial | 0.760 | 0.01 | -0.10 | 0.18 | 0.74 | 89% |
| Prochlorococcus | temporal | 0.753 | 0.07 | -0.04 | 0.18 | 0.70 | 87% |
| Synechococcus | spatial | 0.459 | 0.02 | -0.11 | 0.20 | 0.67 | 89% |
| Synechococcus | temporal | 0.406 | -0.02 | -0.01 | 0.11 | 0.46 | 89% |
| picoeukaryotes | spatial | 0.215 | 0.27 | -0.09 | 0.39 | 0.91 | 89% |
| picoeukaryotes | temporal | 0.223 | 0.34 | 0.01 | 0.39 | 0.83 | 89% |

Inputs that matter most (share of Route A gain, or Route C if used):

- Prochlorococcus: r410_443 19%, r443_551 15%, mbr 12%, lon 9%, par 7%
- Synechococcus: r443_551 21%, r410_443 17%, lon 13%, sst 13%, r486_551 8%
- picoeukaryotes: r486_551 60%, mbr 20%, r443_551 9%, sst 5%, lon 1%

## The reconstruction

- 674 VIIRS 8-day periods on disk, 2012–2026; maps in `data/products/reconstruction/snpp/8D/` (median, 5% and 95% per group, plus an out-of-range count).
- Share of study-box pixels with a value, by season: fall_intermonsoon 52%, ne_monsoon 59%, spring_intermonsoon 51%, sw_monsoon 29%.

![Regional relative index, 2012–2026](figures/reconstruction_series.png)

Each line is a regional monthly mean of the 8-day medians, as a relative index ((value − PACE-era mean) / PACE-era SD). The shaded span is the training period. The table view is `data/products/reconstruction/snpp_8D_regional_means.csv`. Cloud makes June–September means rest on few pixels.

Share of reconstructed pixels with at least one input outside the 0.5–99.5% range seen in training (all groups use the same inputs, so one column):

| year | share |
| --- | --- |
| 2012 | 11.4% |
| 2013 | 11.1% |
| 2014 | 9.9% |
| 2015 | 8.9% |
| 2016 | 9.8% |
| 2017 | 9.6% |
| 2018 | 9.0% |
| 2019 | 9.2% |
| 2020 | 9.9% |
| 2021 | 8.9% |
| 2022 | 9.4% |
| 2023 | 9.1% |
| 2024 | 9.2% |
| 2025 | 9.8% |
| 2026 | 9.5% |

About 9% is the floor even inside the training years: with 12 inputs each allowed a 1% tail, some in-range pixels are always flagged. Years well above that floor are where the inputs drift from what the models saw.

## Bottlenecks

All groups:

- **The labels are MOANA estimates, not cell counts.** The models learn MOANA, errors included. No PACE-era cell counts exist in the box (trust gate, `docs/reports/trust_gate_fallback.md`).
- **Short training record.** 115 periods over about 2.6 years. Extending to 2012 assumes the VIIRS-to-MOANA relation holds in other years; the out-of-range table tracks input shift, not a change in the relation.
- **Monsoon cloud.** The per-period sample gives June–September 37% of training pixels, but only from the cloud-free parts of each scene; in the reconstruction 29% of study-box pixels have a value in the SW monsoon against 54% in other seasons. SW monsoon results are reported separately (PRD section 2.1).
- **Missing inputs.** No mixed-layer depth yet (needs a Copernicus Marine account). Route B and the virtual-VIIRS step wait for the VIIRS response tables.
- **Sensor drift.** Study-box index, 2023–2025 mean minus 2012–2014 mean (SD units): Prochlorococcus +0.17, Synechococcus +0.12, picoeukaryotes -0.14. Whether such a shift is in the ocean or in one sensor is judged against the other VIIRS record in `docs/reports/sensor_comparison.md` (RQ4); nothing is called a trend before the RQ5 trend tests.

Prochlorococcus (exploratory):

- Check A *not supported* (12 matched station-days, one cruise); check C *inconclusive*.
- MOANA's Prochlorococcus is a linear fit made in the Atlantic. Of the PACE-era pixels where MOANA ran, 7.8% were left out of training: 5.4% for a negative Prochlorococcus, the rest above NASA's valid maximum for some group or otherwise flagged. Positive values sit well above NASA's Atlantic output (open item 15).
- MOANA computes Prochlorococcus with an SST term, yet SST carries only 6% of the model's gain; the VIIRS band ratios carry most of it.

Synechococcus (exploratory):

- Check C *inconclusive*; check B found no relation with flow-through phycoerythrin.
- How much Synechococcus information survives in the broad VIIRS bands (RQ2) is untested until virtual VIIRS is built.

Picoeukaryotes (relative index):

- Passed check C against 1991–2005 counts, a pattern check only: the record is a relative index, never cells/mL. MOANA picoeukaryotes follow chlorophyll; a chlorophyll-only comparison is still due.

