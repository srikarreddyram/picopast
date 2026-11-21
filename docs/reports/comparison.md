# Every route, every group, side by side

Written 2026-11-21 from the reports linked in each section; no new computation.
The question: after Routes A–D, which record should each group use, and where did
performance actually improve?

## Two different yardsticks

1. **Ship counts (check C).** The only test against real cells: Spearman ρ between
   a model's same-month field and 1991–2005 surface counts in the box, 95% interval
   from resampling cruises, *pattern consistent* if the lower end is above 0
   (`trust_gate_fallback.md`, `route_d.md`).
2. **Reproducing PACE MOANA from VIIRS (held-out skill).** How well five VIIRS bands
   rebuild NASA's MOANA map (skill = 1 − RMSE / RMSE of a seasonal climatology,
   held-out 2° blocks and the last 6 months). This scores *copying MOANA*, not
   truth, and applies only to the routes that copy MOANA (A, B, C).

A route can score well on (2) and still inherit MOANA's failure on (1). That is
what happened to the cyanobacteria.

## 1. Against ship counts (check C)

ρ [95% interval]; stations / cruises as in check C. **Bold** = the record used.

| model | picoeukaryotes | Synechococcus | Prochlorococcus |
| --- | --- | --- | --- |
| chlorophyll alone (sign fixed in advance) | 0.48 | 0.36 | 0.44 (lower chlorophyll → more) |
| MOANA (Atlantic-fitted) | **0.51 [0.24, 0.64]** | 0.18 [−0.09, 0.38] | 0.09 [−0.21, 0.36] |
| D1, global ships, unweighted | 0.19 [−0.09, 0.38] | 0.32 [0.08, 0.58] | 0.37 [0.22, 0.69] |
| D2, global ships, weighted to the box | 0.52 [0.24, 0.66] | **0.50 [0.31, 0.64]** | **0.45 [0.25, 0.63]** |
| D3, D2 + other Indian cruises (leave one out) | 0.48 [0.08, 0.65] | 0.51 [0.33, 0.63] | 0.41 [0.13, 0.71] |
| verdict for the record used | pattern consistent | pattern consistent | pattern consistent |

Paired differences on the same stations (cruise-resampled):

| | picoeukaryotes | Synechococcus | Prochlorococcus |
| --- | --- | --- | --- |
| D2 − MOANA | −0.01 [−0.15, +0.17] | **+0.32 [+0.18, +0.51]** | **+0.35 [+0.13, +0.64]** |
| D2 − chlorophyll | +0.04 [−0.10, +0.17] | +0.14 [−0.05, +0.26] | +0.01 [−0.14, +0.06] |
| D2 given chlorophyll (partial ρ) | 0.23 [−0.07, 0.49] | **0.38 [0.12, 0.48]** | 0.13 [−0.10, 0.26] |

- **Improved:** *Synechococcus* and *Prochlorococcus* went from inconclusive to
  passing, each clearly better than MOANA.
- **Not improved:** picoeukaryotes (D2 ties MOANA, so MOANA stays, as
  pre-registered).
- **Beyond chlorophyll:** only *Synechococcus* D2 carries information that
  chlorophyll alone does not (partial ρ interval above 0). *Prochlorococcus* D2
  ranks stations about as well as "lower chlorophyll, more *Prochlorococcus*".
  Its added value is a mapped, seasonally resolved record with the right
  direction, not extra skill.
- **The weighting is what makes it work:** D1 fails for picoeukaryotes and is
  weaker for both cyanobacteria. Robustness: leaving out any single contributor,
  or the whole Pacific warm pool, keeps all three D2 groups passing (lowest:
  *Synechococcus* 0.32 [0.13, 0.51]).

## 2. Reproducing MOANA from VIIRS (only for MOANA-copying routes)

Skill on held-out blocks, real VIIRS SNPP (spatial / temporal).

| route | picoeukaryotes | Synechococcus | Prochlorococcus |
| --- | --- | --- | --- |
| A: direct, LightGBM on VIIRS bands + SST + PAR | **0.39 / 0.39** | 0.20 / 0.11 | 0.18 / 0.18 |
| B: rebuild the PACE spectrum, then MOANA | 0.28 / 0.32 | 0.06 / 0.10 | 0.06 / 0.11 |
| C: SST + PAR only | −0.09 / 0.01 | −0.11 / −0.01 | −0.10 / −0.04 |
| upper bound: full PACE spectrum | 0.72 / 0.75 | 0.77 / 0.77 | 0.72 / 0.71 |

Route A wins for every group. Of the gap to the PACE upper bound:

- about 0.15 comes from having five bands instead of the full spectrum;
- the larger part (0.39–0.49 for the cyanobacteria) comes from the change of
  sensor itself.

Sources: `reconstruction.md`, `route_b.md`, `rq2_information.md`.

## 3. The 2012–2026 records

| | picoeukaryotes (Route A, MOANA) | Synechococcus, old (Route A) | Synechococcus, Route D | Prochlorococcus, old (Route A) | Prochlorococcus, Route D |
| --- | --- | --- | --- | --- | --- |
| tier | relative index | exploratory | **relative index** | exploratory | **relative index** |
| SNPP − NOAA-20 drift, SD per decade | −0.04 [−0.06, −0.02] | +0.13 [+0.09, +0.17] | −0.03 [−0.05, +0.00] | +0.25 [+0.19, +0.30] | +0.05 [+0.01, +0.09] |
| same, model without SST | — | — | −0.01 [−0.05, +0.01] | — | +0.05 [+0.01, +0.09] |
| pixels inside the training range | ≈ 90% | ≈ 90% | 74% | ≈ 90% | 57% |

- **Sensor drift:** five times smaller for *Prochlorococcus* and gone for
  *Synechococcus*. The old records' "SNPP-only rise" was a sensor artefact.
- **The old and new records disagree on *when* the cyanobacteria were high or
  low** (monthly anomaly r between −0.48 and +0.15). They are not two versions of
  one signal.
- **Pixel coverage:** Route D's training data cover less of the box's
  conditions. *Prochlorococcus* training SST ends at 29.9 °C, so the warmest
  water is outside it.

Sources: `sensor_comparison.md`, `route_d_reconstruction.md`.

## 4. Climate associations, by the record each group now uses

SD of each index; "no SST" = the same model trained without SST. "Survives"
compares the two.

| group | El Niño / positive IOD | marine heatwaves | trend 2012–2025 | how robust |
| --- | --- | --- | --- | --- |
| picoeukaryotes | lower: r −0.59 (ONI), −0.55 (DMI), study box. On the maps the drop is in the Arabian Sea; most of the Bay of Bengal does not move | lower almost everywhere: −0.34 mean, 96% of cells clear | none | 71–74% survives without SST; NOAA-20 maps agree (r 0.96–0.99) |
| Synechococcus (Route D) | lower: r −0.56 (ONI), −0.60 (DMI), study box | lower: −0.56 mean, 99% of cells clear | Bay of Bengal −0.24 per decade (FDR p 0.006) | **mostly SST:** without SST r −0.41 / −0.34, heatwave map −0.25 (pattern r 0.58), trend gone (−0.10, p 0.48) |
| Prochlorococcus (Route D) | weak: higher in the Arabian Sea after positive IOD (+0.18 mean, 22% of cells clear) | **higher**: +0.26 mean, 68% of cells clear; study box +0.15 [+0.06, +0.20] with all valid pixels | none | survives without SST (maps +0.24, pattern r 0.89): it comes through low chlorophyll and high light, not the SST input |

What changed against the old exploratory cyanobacteria records:

- **Trends.** The old records' significant 2012–2025 rises (+0.09 to +0.24 per
  decade) are gone. They matched the SNPP-only sensor drift.
- **Synechococcus during heatwaves.** The sign flipped, from +0.10 (old,
  Arabian Sea) to −0.42 / −0.52 (Route D). Most of the new value is the SST
  input itself, so the safe claim is weak: without SST it is −0.10
  [−0.17, −0.03] in the Arabian Sea and −0.17 [−0.34, −0.08] in the Bay of
  Bengal.
- **Prochlorococcus during heatwaves.** The direction (higher) is the same in
  both records and holds without SST: the cyanobacteria finding with the best
  support.

Found while checking this:

- **Route D heatwave results.** Route D's reliable pixels thin out during
  heatwaves (Prochlorococcus 48% vs 60%), which had hidden the
  *Prochlorococcus* heatwave rise in the regional average. The record now
  reports reliable-pixel and all-valid-pixel results side by side
  (`route_d_reconstruction.md`).
- **Picoeukaryote heatwave result.** It is not affected: Route A's in-range
  share is the same during heatwaves (≈ 90%).

## Which record each group uses

| group | record | why |
| --- | --- | --- |
| picoeukaryotes | MOANA learned from VIIRS (Route A), `data/products/reconstruction/snpp_8D/` | passed check C; Route D ties it, the pre-registered rule keeps MOANA |
| Synechococcus | Route D (D2), `data/products/reconstruction/snpp_route_d/` | passed check C, beats MOANA; climate signals mostly via SST, report the no-SST version alongside |
| Prochlorococcus | Route D (D2), same | passed check C, beats MOANA; no better than chlorophyll alone; 57% of pixels inside the training range |

## Still weak, in order of how much it limits the conclusions

1. **No same-era counts.**
   - There are no ship counts in the box after 2005.
   - Every check is against 1991–2005 climatological patterns. None tests
     year-to-year change.
2. **Prochlorococcus training is thin.**
   - Effective sample size is 41, mostly from the Pacific warm pool.
   - Training SST stops at 29.9 °C.
3. **The Synechococcus climate links run through SST.**
   - The Route D model learned that SST matters globally, and ENSO, the IOD
     and heatwaves are all SST signals.
   - So the no-SST version is the defensible one.
4. **No mixed-layer depth yet** (needs the Copernicus account). It is the
   obvious missing control on all three groups, and the likeliest route to
   real gains for *Prochlorococcus*.
5. **Route A as an MLP** is untested. Any gain would be in copying MOANA, which
   helps picoeukaryotes only.
