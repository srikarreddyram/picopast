# Trust-gate fallback checks A, B and C

Generated 2026-11-16 by `picopast fallback-checks`, under the rule fixed in `docs/trust_gate_preregistration.md` (16 November 2026). Check A can support Prochlorococcus only as a relative index, never as a cells/mL label; check B is exploratory and gives no verdict; check C tests patterns against 1991–2005 counts. The operational coefficient mapping is the pre-registered one; ATBD results are reported as a sensitivity check and do not change the decision.

## Decision (operational mapping)

A group is reconstructed, as a standardised relative index and never in cells/mL, only if check C is *pattern consistent* and, for Prochlorococcus, check A is not *not supported*.

| group | reconstructed | because |
| --- | --- | --- |
| Prochlorococcus | no | C: inconclusive; A: not supported |
| Synechococcus | no | C: inconclusive |
| picoeukaryotes | yes, relative index | C: pattern consistent |

## operational mapping

### A. Prochlorococcus vs HPLC divinyl chlorophyll a

**Verdict: not supported** (12 matched station-days < 15; rho -0.05 < 0.4; interval lower end -0.68 <= 0.0)

- Station-days: 27 (0 within 25 km of land, excluded); MOANA maps found for 20 of 20 days needed.
- Spearman, MOANA Prochlorococcus vs DVChl a: ρ = -0.05 [-0.68, 0.62], n = 12 station-days matched (95% interval, station-day bootstrap).
- Day used: {nan: 15, -1.0: 5, 0.0: 4, 1.0: 3} (0 = same day; NaN = unmatched). Unmatched by cruise: {'EKAMSAT-2024-Bay-of-Bengal': 13, 'OXR20250108': 2}.
- Zeaxanthin (descriptive only): ρ = 0.23.
- Share of negative Prochlorococcus retrievals among offshore retrieved pixels on these days: 5.4%.

Sensitivity to the block rule:

| window | min valid | ρ | interval | n |
| --- | --- | --- | --- | --- |
| 1×1 | 1 | -0.11 | [-0.69, 0.57] | 12 |
| 3×3 | 5 | -0.05 | [-0.68, 0.62] | 12 |
| 5×5 | 13 | -0.16 | [-0.70, 0.50] | 12 |

### B. Synechococcus vs flow-through phycoerythrin (Fpe565), exploratory

- Pixel-days with flow-through data: 1320; ρ = -0.01 [-0.10, 0.07], n = 692 matched.


### C. Same-month climatology vs 1991–2005 surface counts

Stations: 257 (7 coastal, excluded). Months covered by the PACE-era climatology: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]. The counts are 20+ years older than the satellite record; this tests where and when values are high or low, not today's values.

| group | verdict | ρ | interval | stations | cruises | zeros dropped | unmatched | ρ, leaving out one cruise |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pro | inconclusive (interval [-0.21, 0.36] spans the thresholds) | 0.09 | [-0.21, 0.36] | 179 | 14 | 14 | 36 | -0.03 to 0.18 |
| syn | inconclusive (interval [-0.09, 0.38] spans the thresholds) | 0.18 | [-0.09, 0.38] | 197 | 16 | 2 | 40 | 0.06 to 0.24 |
| peuk | pattern consistent | 0.51 | [0.24, 0.64] | 181 | 13 | 0 | 39 | 0.44 to 0.55 |

## atbd mapping

### A. Prochlorococcus vs HPLC divinyl chlorophyll a

**Verdict: not supported** (12 matched station-days < 15; rho -0.34 < 0.4; interval lower end -0.90 <= 0.0)

- Station-days: 27 (0 within 25 km of land, excluded); MOANA maps found for 20 of 20 days needed.
- Spearman, MOANA Prochlorococcus vs DVChl a: ρ = -0.34 [-0.90, 0.42], n = 12 station-days matched (95% interval, station-day bootstrap).
- Day used: {nan: 15, -1.0: 5, 0.0: 4, 1.0: 3} (0 = same day; NaN = unmatched). Unmatched by cruise: {'EKAMSAT-2024-Bay-of-Bengal': 13, 'OXR20250108': 2}.
- Zeaxanthin (descriptive only): ρ = -0.48.
- Share of negative Prochlorococcus retrievals among offshore retrieved pixels on these days: 4.9%.

Sensitivity to the block rule:

| window | min valid | ρ | interval | n |
| --- | --- | --- | --- | --- |
| 1×1 | 1 | -0.27 | [-0.86, 0.47] | 12 |
| 3×3 | 5 | -0.34 | [-0.90, 0.42] | 12 |
| 5×5 | 13 | -0.40 | [-0.91, 0.34] | 12 |

### B. Synechococcus vs flow-through phycoerythrin (Fpe565), exploratory

- Pixel-days with flow-through data: 1320; ρ = -0.03 [-0.11, 0.06], n = 692 matched.


### C. Same-month climatology vs 1991–2005 surface counts

Stations: 257 (7 coastal, excluded). Months covered by the PACE-era climatology: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]. The counts are 20+ years older than the satellite record; this tests where and when values are high or low, not today's values.

| group | verdict | ρ | interval | stations | cruises | zeros dropped | unmatched | ρ, leaving out one cruise |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pro | pattern consistent | 0.30 | [0.06, 0.47] | 179 | 14 | 14 | 36 | 0.21 to 0.36 |
| syn | inconclusive (interval [-0.02, 0.43] spans the thresholds) | 0.22 | [-0.02, 0.43] | 197 | 16 | 2 | 40 | 0.14 to 0.29 |
| peuk | pattern consistent | 0.51 | [0.25, 0.65] | 181 | 13 | 0 | 39 | 0.44 to 0.55 |

