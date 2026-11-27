# Final confirmatory test: heatwaves, cyclones and the phytoplankton bloom, whole Indian Ocean

Generated 2026-11-27 by `picopast final-test`. Rules: the "Final confirmatory test" section of `docs/cyclone_heatwave_preregistration.md`, committed before any South Indian Ocean value or any test here was computed. Intervals resample whole storms (5,000 draws).

**Sample:** 7,840 track points from 493 storms (133 North, 360 South Indian Ocean), 2002–2026. Heatwave points: 1066 from 169 storms; normal points: 5674 from 441 storms.

## T1: does a heatwave change the bloom? (decisive)

**Verdict: no meaningful effect.** Heatwave effect on the bloom -0.008 [-0.030, +0.014] log10 (99% interval), i.e. -2% [-7%, +3%] chlorophyll; 5,974 points, 148 heatwave storms, 415 normal storms. The smallest meaningful effect was fixed in advance at ±0.05 log10 (a third of the typical bloom).

Heterogeneity guard: heatwave × South Indian Ocean interaction +0.027 [-0.026, +0.086] (95%) — includes zero, so the pooled verdict stands.

| region | heatwave effect, 99% | verdict | points | heatwave / normal storms |
| --- | --- | --- | --- | --- |
| North Indian Ocean | -0.052 [-0.129, +0.024] | inconclusive | 752 | 27 / 88 |
| South Indian Ocean | -0.005 [-0.026, +0.017] | no meaningful effect | 5222 | 121 / 327 |

### T5: power simulation (computed before T1, from exposure counts only)

Resampling North Indian Ocean storms to 169 heatwave and 441 normal storms (500 runs):

| planted effect | chance of each T1 outcome |
| --- | --- |
| +0.00 log10 | no meaningful effect 94%; inconclusive 5%; weaker bloom over heatwaves 1%; stronger bloom over heatwaves 1% |
| -0.05 log10 | weaker bloom over heatwaves 98%; inconclusive 2%; a small weakening, less than a third of the typical bloom 0%; no meaningful effect 0% |

## T2–T4: mechanism, timing and control (95%, Holm across the four)

The pre-registration names three secondary tests; T3 has two components with their own verdicts, so the Holm family here has four members (the stricter reading).

| test | heatwave effect or interaction [95%] | Holm p | verdict | points | heatwave / normal storms |
| --- | --- | --- | --- | --- | --- |
| T2 heatwave x mixing power | -0.012 [-0.026, +0.003] | 0.335 | no detectable dependence on storm power | 5974 | 148 / 415 |
| T3 first period (mostly mixed-up chlorophyll) | -0.012 [-0.033, +0.009] | 0.570 | no detectable difference | 5837 | 143 / 408 |
| T3 second period (later growth) | +0.004 [-0.010, +0.019] | 0.605 | no detectable difference | 5677 | 138 / 402 |
| T4 storm minus no-storm control | -0.035 [-0.057, -0.011] | 0.016 | weaker over heatwaves | 5404 | 99 / 409 |

Heatwave effect on the bloom by storm mixing power (T2 model, 95%):

- weak or fast storms (M 25th pct): -0.001 [-0.016, +0.015]
- strong or slow storms (M 75th pct): -0.018 [-0.040, +0.005]

Sensitivity, M3 = log(wind³ / speed): interaction -0.012 [-0.026, +0.002].

## RQ1: cooling over heatwave water (95%)

- pooled: -0.41 [-0.52, -0.29] °C (6173 points)
- South Indian Ocean: -0.35 [-0.48, -0.22] °C (5241 points)

## RQ1: does the heatwave survive the storm? (95%)

At heatwave points: the share of the disc's heatwave cells on D0 − 3 still in a heatwave on D0 + 5, under the storm vs water in the same region more than 750 km from any storm.

- pooled: 22% under the storm vs 47% far away; difference -0.25 [-0.30, -0.19] (1060 points, 168 storms)
- South Indian Ocean: 28% under the storm vs 52% far away; difference -0.25 [-0.31, -0.18] (771 points, 123 storms)

Heatwave share of track points by region and era: Arabian Sea / modis 10%; Arabian Sea / viirs 30%; Bay of Bengal / modis 4%; Bay of Bengal / viirs 29%; South Indian Ocean / modis 6%; South Indian Ocean / viirs 15%.
