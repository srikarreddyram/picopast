# Route D: a globally trained, Indian-Ocean-adjusted model — rules fixed before testing

Status: **binding from 2026-11-20**, written and committed before any Route D
prediction is compared with Indian Ocean counts. Changes go in the change log
at the end, with the reason, and every earlier run stays reported.

## Why

MOANA was fitted to Atlantic cruises. In this box its *Synechococcus* and
*Prochlorococcus* failed check C (ρ 0.18 and 0.09, inconclusive), and its
*Prochlorococcus* extrapolates to a nearly flat 400,000–490,000 cells/mL
(`docs/open_items.md`, items 15 and 17). Route D asks whether a model trained
on ship counts from the rest of the world ocean, with inputs a satellite era
can supply, and adjusted towards Indian Ocean conditions, ranks the Indian
Ocean stations better than MOANA does.

## Data

- **Training:** MAREDAT picophytoplankton compilation (Buitenhuis et al.
  2012, ESSD 4:37; PANGAEA doi:10.1594/PANGAEA.777385, CC-BY 3.0), samples at
  ≤ 10 m depth, years 1980 or later, counts > 0, per group.
  **Every sample in the Indian Ocean is excluded** (longitude 20–120°E and
  latitude ≤ 31°N), so no test cruise can leak into training.
- **Test:** the check C stations exactly as used for MOANA (BCO-DMO and
  PANGAEA, 1991–2005, ≤ 10 m, in the box, de-duplicated; same offshore rule,
  same 3×3 matching, same month).

## Inputs (all available globally for 1987–2011 and over the box today)

| input | training (at each sample) | test and reconstruction (box) |
| --- | --- | --- |
| SST | NOAA OISST v2 1° monthly, the sample's own month | OISST v2.1 0.25°, monthly mean of the same years as the field being built |
| log10 chlorophyll | MODIS-Aqua monthly climatology (2002–2025), 9 km | VIIRS SNPP chlorophyll |
| PAR | MODIS-Aqua monthly climatology, 9 km | VIIRS SNPP PAR |
| surface nitrate | World Ocean Atlas 2023, monthly 1° climatology, 0 m | the same |
| month | as seasonal sine and cosine, sign-flipped south of the equator | the same |

No latitude or longitude: the model must transfer between basins on
conditions, not positions.

## Models

- **D1:** LightGBM on log10 cells/mL, fixed settings (300 trees, 31 leaves,
  learning rate 0.05, min_child_samples 50), unweighted.
- **D2 (primary):** D1 with importance weights that make the training
  samples resemble the box. A logistic classifier on the inputs separates
  training samples from box pixels (PACE-era monthly fields, offshore);
  weight = p(box | x) / p(train | x), clipped at its 99th percentile and
  normalised to mean 1.
- **D3 (secondary):** D2 fine-tuned with Indian Ocean stations, **leave one
  cruise out**: for each test cruise, D2 plus the other cruises' stations
  (weight 5 each) predicts that cruise. Reported, never used for a verdict on
  its own.

Settings are not tuned on the test stations. A leave-one-basin-out check on
the training data (train without the Atlantic, test on it; and the reverse)
is reported to show whether transfer between basins works at all.

## Test

Each model builds a same-calendar-month field over the box from PACE-era
(2024–2026) inputs, as MOANA's check C climatology is built, and goes
through check C unchanged: Spearman ρ across stations, 95% interval from
resampling cruises, the same verdict rule (*pattern consistent* if the
lower end of the interval is above 0 with ≥ 30 stations from ≥ 3 cruises).
Also reported: the paired difference ρ(D) − ρ(MOANA) on stations matched by
both, with a cruise-resampled interval, and the same against chlorophyll
alone.

## Decision

- **Synechococcus, Prochlorococcus** (MOANA inconclusive): if **D2** is
  *pattern consistent*, the group is reconstructed from D2 as a relative
  index (tier "relative index, Route D"), 2012–2026, from VIIRS SNPP
  chlorophyll and PAR, OISST and WOA nitrate per 8-day period. Otherwise it
  stays exploratory, and Route D is reported as a negative result.
- **Picoeukaryotes** (MOANA passed): MOANA stays the base unless the paired
  difference interval for D2 is entirely above zero.
- Because SST is a D input, any climate association of a Route D record is
  repeated with D2 trained without SST, as for the main record.

## Known limits, stated with every result

- Training counts are 1987–2011 and chlorophyll/PAR are climatologies; the
  model learns typical conditions, not same-day ones.
- Most Atlantic training samples are from the cold North Atlantic; D2's
  weights address this only where the two overlap.
- Nitrate is a fixed climatology, so year-to-year changes reach the model
  only through SST, chlorophyll and PAR.

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-11-20 | First version, binding | Written before any Route D prediction was compared with Indian Ocean counts. |
