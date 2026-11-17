# Does MOANA add to chlorophyll? Check C with a chlorophyll climatology

Generated 2026-11-17 by `picopast chl-comparison` (operational mapping). Descriptive, added 17 November 2026 after the trust-gate verdicts; it does not change them. Same stations, same-month matching and same 3×3 rule as check C, with a PACE-era VIIRS SNPP chlorophyll climatology (mean log10 chlorophyll by calendar month) in place of MOANA. Intervals: 95%, resampling whole cruises; only stations matched by both are used.

| group | stations | cruises | ρ MOANA | ρ chlorophyll | MOANA − chlorophyll | partial ρ, MOANA given chlorophyll |
| --- | --- | --- | --- | --- | --- | --- |
| Prochlorococcus | 175 | 14 | 0.07 | -0.43 | +0.50 [+0.12, +0.93] | -0.01 [-0.25, 0.17] |
| Synechococcus | 193 | 16 | 0.19 | 0.35 | -0.17 [-0.49, -0.03] | 0.10 [-0.19, 0.37] |
| picoeukaryotes | 177 | 13 | 0.52 | 0.47 | +0.05 [-0.08, +0.14] | 0.25 [-0.01, 0.45] |

- Prochlorococcus: counts fall as chlorophyll rises (ρ -0.43), so the fair baseline is *low* chlorophyll. Against it, MOANA − baseline = -0.35 [-0.73, -0.11]: low chlorophyll alone ranks the stations better. The partial ρ does not depend on the baseline's sign.

Reading: a difference interval above zero means MOANA ranks the stations better than chlorophyll does; a partial ρ interval above zero means MOANA carries information about the counts that chlorophyll does not. Counts are 1991–2005; the climatologies are 2024–2026.
