# Check E: new, independent HPLC stations vs Route D

Generated 2026-11-21 by `picopast check-e`. Rules: `docs/check_e_preregistration.md`, committed before this ran. Divinyl chlorophyll a (made only by *Prochlorococcus*) scores *Prochlorococcus*; zeaxanthin, made by both cyanobacteria, is shown for *Synechococcus* as a description only. Pigments are not cell counts: these are rank tests.

Station-days (surface, one UTC day × one 4.6 km pixel): 149, of which 0 coastal (excluded by check C's offshore rule).

## Prochlorococcus vs divinyl chlorophyll a

149 offshore station-days, 127 dated 2012 or later. Per cruise: EKAMSAT-2024-Bay-of-Bengal 25, EKAMSAT_PILOT_2023 21, I07N_2018 16, I09N 39, OXR20250108 2, RR2307_2023 1, TN444 23, abudhabi-masqat 1, i8si9n 14, male-stbrandon 1, masqat-mumbay 3, mumbay-male 3.

| test | model | Spearman ρ [95%, cruises resampled] |
| --- | --- | --- |
| E1 pattern (same month, all dates) | Route D (D2) | -0.02 [-0.18, +0.16] (n 147, 12 cruises) |
| | chlorophyll alone | -0.09 [-0.31, +0.15] (n 147, 12 cruises) |
| | MOANA | -0.04 [-0.24, +0.25] (n 147, 12 cruises) |
| E2 record (own 8-day period, 2012 on) | **Route D record, all valid pixels** | -0.07 [-0.30, +0.16] (n 69, 6 cruises) |
| | Route D record, reliable pixels only | -0.49 [-0.80, -0.25] (n 18, 4 cruises) |
| | VIIRS 8-day chlorophyll alone | -0.13 [-0.45, +0.00] (n 69, 6 cruises) |
| | old MOANA-learned record (Route A) | +0.13 [-0.12, +0.39] (n 69, 6 cruises) |

| paired comparison, same stations | difference in ρ [95%] | partial ρ of the first given the second |
| --- | --- | --- |
| E1: D2 − chlorophyll | +0.07 [-0.07, +0.20] | +0.03 [-0.08, +0.20] |
| E1: D2 − MOANA | -0.00 [-0.38, +0.30] | -0.04 [-0.21, +0.14] |
| E2: record − chlorophyll | +0.06 [-0.11, +0.37] | -0.02 [-0.19, +0.35] |
| E2: record − Route A | -0.20 [-0.54, +0.04] | -0.07 [-0.28, +0.15] |
| E2: record − its own same-month climatology | +0.01 [-0.22, +0.13] | -0.03 [-0.30, +0.12] |

**Verdicts:** E1 *pattern not consistent*; E2 *pattern not consistent*. Primary test (E2) → **downgraded to exploratory**.

## Synechococcus vs zeaxanthin (descriptive only)

149 offshore station-days, 127 dated 2012 or later. Per cruise: EKAMSAT-2024-Bay-of-Bengal 25, EKAMSAT_PILOT_2023 21, I07N_2018 16, I09N 39, OXR20250108 2, RR2307_2023 1, TN444 23, abudhabi-masqat 1, i8si9n 14, male-stbrandon 1, masqat-mumbay 3, mumbay-male 3.

| test | model | Spearman ρ [95%, cruises resampled] |
| --- | --- | --- |
| E1 pattern (same month, all dates) | Route D (D2) | -0.22 [-0.47, +0.12] (n 147, 12 cruises) |
| | chlorophyll alone | -0.05 [-0.46, +0.29] (n 147, 12 cruises) |
| | MOANA | +0.13 [-0.17, +0.39] (n 147, 12 cruises) |
| E2 record (own 8-day period, 2012 on) | **Route D record, all valid pixels** | -0.09 [-0.47, +0.33] (n 69, 6 cruises) |
| | Route D record, reliable pixels only | +0.05 [-0.29, +0.81] (n 34, 5 cruises) |
| | VIIRS 8-day chlorophyll alone | +0.09 [-0.11, +0.52] (n 69, 6 cruises) |
| | old MOANA-learned record (Route A) | +0.13 [-0.09, +0.58] (n 69, 6 cruises) |

| paired comparison, same stations | difference in ρ [95%] | partial ρ of the first given the second |
| --- | --- | --- |
| E1: D2 − chlorophyll | -0.17 [-0.40, +0.13] | -0.24 [-0.41, +0.12] |
| E1: D2 − MOANA | -0.36 [-0.64, +0.07] | -0.25 [-0.48, +0.09] |
| E2: record − chlorophyll | -0.18 [-0.42, -0.01] | -0.14 [-0.47, +0.23] |
| E2: record − Route A | -0.22 [-0.74, +0.19] | -0.13 [-0.48, +0.25] |
| E2: record − its own same-month climatology | +0.07 [-0.34, +0.60] | +0.02 [-0.49, +0.49] |

## Why the pigment test fails (post hoc, not part of the verdict)

These were looked at after the verdict above and do not change it.

- Divinyl chlorophyll a at these stations spans only 0.020–0.116 mg m⁻³ and follows **total** chlorophyll a measured on the same samples (Spearman ρ = 0.76, n 149; *Prochlorococcus*'s share of total chlorophyll a, median 0.41). It is tracking overall phytoplankton biomass.
- The 1991–2005 cell counts in check C did the opposite: *Prochlorococcus* cells were most numerous where chlorophyll was lowest (which is why chlorophyll's sign is fixed negative). Cells in richer water hold more pigment each, so in this box the pigment and the cell count rank stations in nearly opposite directions. No satellite field ranks the pigment: Route D, MOANA and chlorophyll alone all fail it.
- Within single cruises the record–pigment correlation has no consistent sign: EKAMSAT-2024-Bay-of-Bengal +0.24 (n 20); I07N_2018 +0.10 (n 5); I09N -0.16 (n 32); TN444 +0.06 (n 9).
- So check E most likely measures the limits of divinyl chlorophyll a as a *cell-count* proxy here, as much as the model. It is still binding: the pre-registered rule downgrades *Prochlorococcus* to exploratory. Only new cell counts can settle which test is right.

## Limits

- Pigment per cell changes with light and nutrients, so a pigment rank is not a cell-count rank.
- Stations come from a handful of cruises; intervals resample whole cruises and are wide.
- Same-period matching loses stations under cloud; there is no gap filling.
- No new *Synechococcus* or *Prochlorococcus* cell counts exist in the open archives searched (PANGAEA, BCO-DMO, gVOD, SeaBASS), so *Synechococcus* cannot be re-tested here.
