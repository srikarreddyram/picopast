# Route D: a globally trained, Indian-Ocean-adjusted model

Generated 2026-11-20 by `picopast route-d`. Rules: `docs/route_d_preregistration.md` (committed before this test). Trained on MAREDAT ship counts from the rest of the world ocean (5,445 surface stations with every input; 230 Indian Ocean stations excluded), inputs SST, chlorophyll, PAR, surface nitrate and season, no position. Tested exactly like MOANA in check C: the model's PACE-era same-month field against the 1991–2005 Indian Ocean stations, Spearman ρ with 95% intervals from resampling cruises.

## Check C, side by side

| group | MOANA | D1 (global) | D2 (adjusted, primary) | D3 (fine-tuned, leave one cruise out) | D2 verdict |
| --- | --- | --- | --- | --- | --- |
| picoeukaryotes | 0.51 [0.24, 0.64] (pattern consistent) | 0.19 [-0.09, 0.38] | 0.52 [0.24, 0.66] | 0.48 [0.08, 0.65] | pattern consistent |
| Synechococcus | 0.18 [-0.09, 0.38] (inconclusive) | 0.32 [0.08, 0.58] | 0.50 [0.31, 0.64] | 0.51 [0.33, 0.63] | pattern consistent |
| Prochlorococcus | 0.09 [-0.21, 0.36] (inconclusive) | 0.37 [0.22, 0.69] | 0.45 [0.25, 0.63] | 0.41 [0.13, 0.71] | pattern consistent |

## Paired tests (stations matched by every field)

| group | model | stations (cruises) | ρ(model) − ρ(MOANA) | ρ(model) − ρ(chlorophyll) |
| --- | --- | --- | --- | --- |
| picoeukaryotes | D1 | 177 (13) | -0.29 [-0.58, -0.02] | -0.29 [-0.55, -0.03] vs chlorophyll |
| picoeukaryotes | D2 | 177 (13) | -0.01 [-0.15, +0.17] | +0.04 [-0.10, +0.17] vs chlorophyll |
| Synechococcus | D1 | 193 (16) | +0.13 [-0.12, +0.47] | -0.04 [-0.31, +0.26] vs chlorophyll |
| Synechococcus | D2 | 193 (16) | +0.32 [+0.18, +0.51] | +0.14 [-0.05, +0.26] vs chlorophyll |
| Prochlorococcus | D1 | 175 (14) | +0.29 [+0.05, +0.69] | -0.07 [-0.16, +0.06] vs low chlorophyll |
| Prochlorococcus | D2 | 175 (14) | +0.35 [+0.13, +0.64] | +0.01 [-0.14, +0.06] vs low chlorophyll |

## Decision (pre-registered rule)

- **picoeukaryotes:** MOANA stays the base (D2 is not clearly better on the paired test).
- **Synechococcus:** reconstruct from D2 as a relative index (Route D).
- **Prochlorococcus:** reconstruct from D2 as a relative index (Route D).

## Diagnostics

| group | training stations (by basin) | effective n after weighting | weight share by basin | leave-one-basin-out ρ | input importance (D2) |
| --- | --- | --- | --- | --- | --- |
| picoeukaryotes | 5,212 (other 2646, atlantic 2566) | 93 | atlantic 7%, other 93% | held-out atlantic -0.02, held-out other 0.38 | logchl 0.28, par 0.27, sst 0.21, no3 0.17, season_sin 0.04, season_cos 0.04 |
| Synechococcus | 4,965 (atlantic 2587, other 2378) | 93 | atlantic 8%, other 92% | held-out atlantic 0.44, held-out other 0.33 | par 0.24, logchl 0.24, sst 0.23, no3 0.19, season_cos 0.05, season_sin 0.04 |
| Prochlorococcus | 1,840 (other 1074, atlantic 766) | 41 | atlantic 12%, other 88% | held-out atlantic 0.41, held-out other 0.23 | sst 0.25, logchl 0.25, par 0.21, no3 0.20, season_sin 0.05, season_cos 0.04 |

### Does the pass rest on a few samples? (added after the test; descriptive)

After weighting, a few dozen samples carry most of D2's weight, mostly from the western equatorial Pacific warm pool (SST about 29 °C, low chlorophyll and nitrate: the closest analogue to this box). D2 refitted without each heavily weighted contributor, and without the warm pool altogether:

| left out of training | picoeukaryotes | Synechococcus | Prochlorococcus |
| --- | --- | --- | --- |
| (none) | 0.52 [0.24, 0.66] ✓ | 0.50 [0.31, 0.64] ✓ | 0.45 [0.25, 0.63] ✓ |
| Blanchot (all entries) | 0.61 [0.24, 0.71] ✓ | 0.32 [0.13, 0.51] ✓ | 0.40 [0.21, 0.63] ✓ |
| Landry | 0.58 [0.29, 0.70] ✓ | 0.39 [0.16, 0.64] ✓ | 0.47 [0.29, 0.67] ✓ |
| Buck | 0.54 [0.28, 0.68] ✓ | 0.43 [0.18, 0.62] ✓ | 0.39 [0.17, 0.64] ✓ |
| Shimada/Maruyama | 0.47 [0.22, 0.67] ✓ | 0.40 [0.16, 0.58] ✓ | 0.49 [0.23, 0.69] ✓ |
| Pacific warm pool (east of 120°E or west of 150°W, within 20° of the equator) | 0.50 [0.14, 0.68] ✓ | 0.42 [0.20, 0.65] ✓ | 0.43 [0.25, 0.66] ✓ |

✓ = still *pattern consistent* under the check C rule.

Leave-one-basin-out: D1 trained without one basin and scored on it (Spearman ρ on the samples), a check that conditions, not positions, carry the signal between oceans.

## Limits

- Training counts are 1987–2011 and chlorophyll and PAR are climatologies: the model learns typical conditions, not same-day ones. Counts are 1991–2005; the fields are 2024–2026.
- Most Atlantic samples are from the cold North Atlantic; the weights help only where conditions overlap (effective n above).
- Nitrate is a fixed climatology. SST is an input, so any climate association of a Route D record must be repeated without SST.
