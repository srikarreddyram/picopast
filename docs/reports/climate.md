# Climate associations (RQ5): first pass

Generated 2026-11-17 by `picopast climate`, from the VIIRS SNPP back-extension, reliable pixels only (every model input inside the range seen in training). Values are anomalies of the relative index in SD units. Associations, not causes. Only picoeukaryotes passed the trust gate and agree between the two VIIRS sensors; Synechococcus and Prochlorococcus rows are exploratory and shown for completeness.

## Trends (core seasons, annual means)

Mann–Kendall test and Sen's slope on annual means of the anomalies outside the SW monsoon. Descriptive: a short record and a sensor that may drift. For comparison, SNPP drifts against NOAA-20 by -0.04 (picoeukaryotes), +0.13 (Synechococcus), +0.25 (Prochlorococcus) per decade in the study box (`sensor_comparison.md`): a trend no larger than that drift is not evidence of change in the ocean.

| region | group | tier | years | Sen's slope per decade | p | p (FDR) |
| --- | --- | --- | --- | --- | --- | --- |
| Arabian Sea C/E | picoeukaryotes | relative index | 2012–2025 | +0.05 | 0.743 | 0.835 |
| Arabian Sea C/E | Synechococcus | EXPLORATORY | 2012–2025 | +0.16 | 0.001 | 0.002 |
| Arabian Sea C/E | Prochlorococcus | EXPLORATORY | 2012–2025 | +0.24 | 0.000 | 0.000 |
| Bay of Bengal | picoeukaryotes | relative index | 2012–2025 | -0.10 | 0.080 | 0.103 |
| Bay of Bengal | Synechococcus | EXPLORATORY | 2012–2025 | +0.12 | 0.000 | 0.000 |
| Bay of Bengal | Prochlorococcus | EXPLORATORY | 2012–2025 | +0.17 | 0.001 | 0.002 |
| Study box | picoeukaryotes | relative index | 2012–2025 | -0.02 | 1.000 | 1.000 |
| Study box | Synechococcus | EXPLORATORY | 2012–2025 | +0.09 | 0.001 | 0.002 |
| Study box | Prochlorococcus | EXPLORATORY | 2012–2025 | +0.15 | 0.000 | 0.000 |

## Lagged correlation with DMI (IOD) and ONI (ENSO)

Monthly anomalies against the index 0–6 months earlier. Shown: the strongest lag per region, group and index. p uses an effective sample size for autocorrelation; FDR over all 126 lag tests.

| region | group | tier | index | best lag (months) | r | n (effective) | p | p (FDR) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Arabian Sea C/E | picoeukaryotes | relative index | DMI | 2 | -0.48 | 175 (57) | 0.000 | 0.004 |
| Arabian Sea C/E | picoeukaryotes | relative index | ONI | 0 | -0.54 | 177 (43) | 0.000 | 0.004 |
| Bay of Bengal | picoeukaryotes | relative index | DMI | 5 | -0.08 | 177 (62) | 0.512 | 0.679 |
| Bay of Bengal | picoeukaryotes | relative index | ONI | 6 | -0.24 | 177 (47) | 0.108 | 0.356 |
| Study box | picoeukaryotes | relative index | DMI | 2 | -0.55 | 175 (51) | 0.000 | 0.003 |
| Study box | picoeukaryotes | relative index | ONI | 0 | -0.59 | 177 (40) | 0.000 | 0.004 |
| Arabian Sea C/E | Prochlorococcus | EXPLORATORY | DMI | 6 | +0.17 | 177 (70) | 0.159 | 0.408 |
| Arabian Sea C/E | Prochlorococcus | EXPLORATORY | ONI | 6 | +0.13 | 177 (56) | 0.326 | 0.579 |
| Bay of Bengal | Prochlorococcus | EXPLORATORY | DMI | 3 | +0.35 | 176 (52) | 0.010 | 0.068 |
| Bay of Bengal | Prochlorococcus | EXPLORATORY | ONI | 0 | +0.07 | 177 (39) | 0.663 | 0.787 |
| Study box | Prochlorococcus | EXPLORATORY | DMI | 0 | +0.10 | 173 (54) | 0.474 | 0.664 |
| Study box | Prochlorococcus | EXPLORATORY | ONI | 1 | -0.15 | 177 (43) | 0.337 | 0.579 |
| Arabian Sea C/E | Synechococcus | EXPLORATORY | DMI | 1 | +0.13 | 174 (80) | 0.251 | 0.531 |
| Arabian Sea C/E | Synechococcus | EXPLORATORY | ONI | 6 | +0.19 | 177 (66) | 0.118 | 0.356 |
| Bay of Bengal | Synechococcus | EXPLORATORY | DMI | 1 | +0.32 | 174 (57) | 0.016 | 0.099 |
| Bay of Bengal | Synechococcus | EXPLORATORY | ONI | 0 | +0.23 | 177 (46) | 0.127 | 0.356 |
| Study box | Synechococcus | EXPLORATORY | DMI | 3 | -0.20 | 176 (56) | 0.139 | 0.382 |
| Study box | Synechococcus | EXPLORATORY | ONI | 3 | -0.25 | 177 (44) | 0.109 | 0.356 |

## The 2019 positive IOD

September–December 2019 anomaly against the same months of IOD-neutral years (|DMI| < 0.4, September–November mean). z = (2019 − neutral mean) / neutral SD.

| region | group | tier | 2019 | neutral mean ± SD | z | neutral years |
| --- | --- | --- | --- | --- | --- | --- |
| Arabian Sea C/E | picoeukaryotes | relative index | -0.02 | +0.02 ± 0.23 | -0.1 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Arabian Sea C/E | Synechococcus | EXPLORATORY | +0.09 | +0.01 ± 0.09 | +0.9 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Arabian Sea C/E | Prochlorococcus | EXPLORATORY | +0.14 | -0.03 ± 0.16 | +1.1 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Bay of Bengal | picoeukaryotes | relative index | +0.01 | -0.01 ± 0.14 | +0.2 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Bay of Bengal | Synechococcus | EXPLORATORY | +0.01 | -0.02 ± 0.05 | +0.5 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Bay of Bengal | Prochlorococcus | EXPLORATORY | +0.06 | -0.02 ± 0.11 | +0.8 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Study box | picoeukaryotes | relative index | -0.12 | +0.01 ± 0.13 | -1.0 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Study box | Synechococcus | EXPLORATORY | -0.06 | +0.00 ± 0.05 | -1.1 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |
| Study box | Prochlorococcus | EXPLORATORY | +0.03 | -0.02 ± 0.10 | +0.6 | 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2024, 2025 |

## Marine heatwaves

Heatwave periods: at least half of the region's ocean pixel-days in a heatwave (Hobday definition, OISST, 1991–2020 baseline). Normal periods: the same calendar months with under 10%. Difference = heatwave mean − normal mean; 95% interval from resampling whole years.

| region | group | tier | heatwave periods (years) | normal periods | difference | 95% interval |
| --- | --- | --- | --- | --- | --- | --- |
| Arabian Sea C/E | picoeukaryotes | relative index | 92 (11) | 369 | -0.36 | [-0.43, -0.25] |
| Arabian Sea C/E | Synechococcus | EXPLORATORY | 92 (11) | 369 | +0.10 | [+0.04, +0.16] |
| Arabian Sea C/E | Prochlorococcus | EXPLORATORY | 92 (11) | 369 | +0.15 | [+0.08, +0.24] |
| Bay of Bengal | picoeukaryotes | relative index | 58 (10) | 336 | -0.27 | [-0.35, -0.21] |
| Bay of Bengal | Synechococcus | EXPLORATORY | 58 (10) | 336 | +0.02 | [-0.02, +0.07] |
| Bay of Bengal | Prochlorococcus | EXPLORATORY | 58 (10) | 336 | +0.12 | [+0.04, +0.18] |
| Study box | picoeukaryotes | relative index | 42 (6) | 356 | -0.25 | [-0.33, -0.15] |
| Study box | Synechococcus | EXPLORATORY | 42 (6) | 356 | -0.03 | [-0.08, +0.00] |
| Study box | Prochlorococcus | EXPLORATORY | 42 (6) | 356 | +0.06 | [+0.01, +0.15] |

### Sensitivity: picoeukaryotes from a model without SST

SST is a model input, so part of a heatwave or ENSO association could be put there by the model. This record comes from the same model trained without SST (VIIRS bands, ratios, PAR, position, season only).

| region | heatwave difference [95%] | best ONI lag: r (p) | best DMI lag: r (p) |
| --- | --- | --- | --- |
| Study box | -0.18 [-0.26, -0.05] | lag 0: -0.52 (0.000) | lag 2: -0.47 (0.000) |
| Bay of Bengal | -0.19 [-0.27, -0.13] | lag 5: -0.17 (0.228) | lag 0: +0.05 (0.676) |
| Arabian Sea C/E | -0.25 [-0.32, -0.14] | lag 0: -0.49 (0.001) | lag 2: -0.41 (0.001) |

Heatwave periods per region: Study box 42, Bay of Bengal 58, Arabian Sea C/E 93 of 674 8-day periods.

## Not yet done

- Mixed-layer depth as a covariate (needs the Copernicus Marine account).
- The same analysis on the NOAA-20 record for 2018–2026 (picoeukaryotes), as a sensor check.
- Pixel-level composites and maps; seasonal (monsoon-phase) composites beyond the climatology.
