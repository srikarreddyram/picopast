# RQ2: how much group information survives in VIIRS bands?

Generated 2026-11-17 by `picopast rq2`. Identical models predict our MOANA output (operational mapping) from two views of the same PACE spectrum: all 172 bands, or the spectrum band-averaged through SNPP's M1–M5 responses (virtual VIIRS, with the same band ratios the reconstruction uses). Both views get the same SST, position and season inputs. Same atmosphere and calibration in both, so the gap is the information the broad bands lose.

- 118 PACE-era 8-day periods, 590,000 sampled offshore pixels (the reconstruction's screen: MOANA flags 0 and inside NASA's valid range).
- Skill = 1 − RMSE / RMSE of a seasonal climatology fitted in the same training part; spatial = 2° blocks in 5 folds, pooled; temporal = the last 6 months held out.

| group | split | model | skill, full PACE | skill, virtual VIIRS | lost | ρ, full | ρ, virtual |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Prochlorococcus | spatial | LightGBM | 0.72 | 0.59 | +0.13 | 0.97 | 0.93 |
| Prochlorococcus | spatial | ridge | 0.37 | 0.37 | +0.00 | 0.91 | 0.86 |
| Prochlorococcus | temporal | LightGBM | 0.71 | 0.57 | +0.14 | 0.96 | 0.92 |
| Prochlorococcus | temporal | ridge | 0.45 | 0.43 | +0.02 | 0.92 | 0.87 |
| Synechococcus | spatial | LightGBM | 0.77 | 0.64 | +0.13 | 0.97 | 0.90 |
| Synechococcus | spatial | ridge | 0.33 | 0.34 | -0.01 | 0.87 | 0.77 |
| Synechococcus | temporal | LightGBM | 0.77 | 0.61 | +0.17 | 0.96 | 0.83 |
| Synechococcus | temporal | ridge | 0.37 | 0.39 | -0.02 | 0.85 | 0.68 |
| picoeukaryotes | spatial | LightGBM | 0.72 | 0.56 | +0.16 | 0.99 | 0.96 |
| picoeukaryotes | spatial | ridge | 0.64 | 0.42 | +0.21 | 0.98 | 0.94 |
| picoeukaryotes | temporal | LightGBM | 0.75 | 0.56 | +0.18 | 0.98 | 0.93 |
| picoeukaryotes | temporal | ridge | 0.68 | 0.49 | +0.19 | 0.98 | 0.91 |

"Lost" is the skill the five bands give up against the full spectrum (larger = more group information lost). The PRD's hypothesis: the biggest loss is for Synechococcus.

## Where the skill goes, PACE to real VIIRS (LightGBM)

Three steps on the same PACE-era 8-day periods and the same screen: the full spectrum, the five bands simulated from that spectrum, and real VIIRS SNPP (the reconstruction's Route A, which also has PAR). The first gap is the bands; the second is everything that differs between two sensors looking at the ocean: overpass time, which days are cloud-free within the 8-day composite, calibration and atmospheric correction.

| group | split | full PACE | virtual VIIRS | real VIIRS | lost to bands | lost to sensor mismatch |
| --- | --- | --- | --- | --- | --- | --- |
| Prochlorococcus | spatial | 0.72 | 0.59 | 0.18 | 0.13 | 0.41 |
| Prochlorococcus | temporal | 0.71 | 0.57 | 0.18 | 0.14 | 0.39 |
| Synechococcus | spatial | 0.77 | 0.64 | 0.20 | 0.13 | 0.44 |
| Synechococcus | temporal | 0.77 | 0.61 | 0.11 | 0.17 | 0.49 |
| picoeukaryotes | spatial | 0.72 | 0.56 | 0.39 | 0.16 | 0.17 |
| picoeukaryotes | temporal | 0.75 | 0.56 | 0.39 | 0.18 | 0.17 |

Inputs differ slightly between steps (CMC SST here, OISST and PAR in the reconstruction) and the pixel samples are drawn separately, so the split is approximate.
