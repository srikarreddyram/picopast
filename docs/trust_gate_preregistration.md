# Trust gate pre-registration (RQ1)

**Status: BINDING from 16 November 2026.** Drafted 10 November and fixed
16 November 2026, before any satellite-to-ship matchup was computed. From now
on this file and the `trust_gate` and `matchup` blocks of
`config/picopast.yaml` are frozen; any later change is recorded in the change
log below with a reason.

The thresholds were first drafted as proposals; none comes from the PRD. They
were accepted unchanged on 16 November. The pigment-proxy and climatological
pattern checks below were added that day, after the SeaBASS search found no
PACE-era cell counts in the box but did find HPLC pigments (see
`docs/open_items.md`, items 2 and 5). They were added before any of those
samples was compared with satellite data.

## Question

How accurate is our MOANA run in the North Indian Ocean, for each of
*Prochlorococcus*, *Synechococcus* and picoeukaryotes?

## Data

- **Satellite:** our own MOANA run on PACE OCI L3M Rrs (daily, native 4.6 km
  grid) with GHRSST CMC 0.1° SST, operational coefficient mapping.
- **Ship:** surface cell counts from open archives and published tables only,
  dated 5 March 2024 or later, inside 0–25°N, 50–100°E.
- **Surface** means sample depth ≤ 10 m. Where several depths
  ≤ 10 m exist at one station, the shallowest is used.
- **Excluded before matching:** stations within 25 km of the coast; counts
  that are zero, negative or flagged by the archive.

## Matchup

1. Same calendar day (UTC) as the ship sample. If that fails, the day before,
   then the day after.
2. Value = median of the 3×3 pixel block centred on the station.
3. A block counts only with at least 5 of 9 valid pixels.
4. A pixel is invalid if its MOANA flags are non-zero. For *Prochlorococcus*,
   a negative retrieval is invalid, and the share of negative retrievals is
   reported alongside the verdict.
5. Stations with no valid block in the three days are reported as unmatched,
   with counts per region and season.

## Metrics

All in log10(cells/mL), estimate minus observation: log-RMSE, mean bias,
Spearman ρ, and the share of matchups within a factor of 2. Intervals are 95%
percentile intervals from 2,000 bootstrap resamples of **whole cruises**.

## Decision rule (per group)

| Verdict | Condition |
| --- | --- |
| **Insufficient matchups** | Fewer than 30 matchups. No absolute-error verdict. Fallback: compare spatial and seasonal patterns with published cruise maps, and say so. |
| **Usable label** | All of: log-RMSE ≤ 0.50; \|bias\| ≤ 0.30; lower end of the ρ interval ≥ 0.30; at least 50% within a factor of 2. |
| **Relative index only** | Not a usable label, but ρ ≥ 0.40 and the lower end of its interval is above 0. The group is modelled as a standardised index, not as cells/mL. |
| **Excluded** | Neither of the above. The group is dropped from the reconstruction and the evaluation itself is reported. |

A group needs matchups from at least 2 cruises for a cruise bootstrap to mean
anything; with one cruise the verdict is reported as provisional.

The rule is implemented in `picopast.trust_gate.verdict` and reads its
thresholds from `config/picopast.yaml`.

## When cell counts are missing (added 16 November, before any comparison)

The open archives hold no PACE-era cell counts in the box, so every group is
expected to fall under **Insufficient matchups**. The fallback is fixed as
follows. None of these checks can give a **Usable label**; they can only
support or fail to support **Relative index only**.

### A. Pigment proxy for *Prochlorococcus* (PACE era)

- **Ship data:** HPLC divinyl chlorophyll a (DVChl a, mg m⁻³), which only
  *Prochlorococcus* makes, from SeaBASS (EKAMSAT 2024 Bay of Bengal, OceanX
  January 2025), surface (≤ 10 m), inside the box, dated 5 March 2024 or later.
  Samples with the same `hplc_gsfc_id` are counted once. A **station-day** is
  one UTC date and one pixel of the 4.6 km L3M grid; where it has several
  surface samples (replicates, or one station in two file releases), their
  median is used.
- **Matchup:** as above (same day, then ±1 day; 3×3 median; ≥ 5 of 9 valid
  pixels; coast mask).
- **Metric:** Spearman ρ between MOANA *Prochlorococcus* (cells/mL) and DVChl
  a. There are only two cruises, so the 95% interval comes from 2,000
  bootstrap resamples of **station-days**, not cruises. No log-RMSE or bias is
  computed: the units differ and cells per unit DVChl a varies with light.
- **Verdict:** *relative index supported (provisional, pigment proxy)* if
  ρ ≥ 0.40 and the lower end of the interval is above 0, with at least 15
  matched station-days. Otherwise *not supported*. Always reported as
  provisional: one main cruise, one season (inter-monsoon, April–May).
- Zeaxanthin is reported alongside as a descriptive check only. It is shared
  by *Prochlorococcus* and *Synechococcus*, so it is not used for a verdict.

### B. Phycoerythrin fluorescence for *Synechococcus* (exploratory)

The flow-through fluorometer files (`*_CLS.sb`, Fpe565/578/590 in RFU) track
phycoerythrin, which *Synechococcus* carries. They are uncalibrated and the
pigment is not unique to it, so this comparison is **reported but gives no
verdict**: Spearman ρ between MOANA *Synechococcus* and Fpe565, flow-through
points binned to the 4.6 km grid per day, same matchup rules.

### C. Climatological pattern check (all three groups)

The historical ship counts (BCO-DMO compilation and PANGAEA, 1991–2005,
surface ≤ 10 m, inside the box) are compared with the PACE-era MOANA
climatology for the **same calendar month** at the same pixel (3×3 median of
the monthly-mean MOANA). Metric: Spearman ρ across stations, with a
bootstrap over **cruises**. This tests whether MOANA puts high and low values
in the right places and seasons, not whether today's values match the 1990s.

- **Verdict per group:** *pattern consistent* if the lower end of the ρ
  interval is above 0 with at least 30 stations from at least 3 cruises;
  *pattern not consistent* if the upper end is below 0.30; otherwise
  *inconclusive*.
- The gap of 20+ years is stated with every result.

### How the fallback feeds the reconstruction

A group is reconstructed as a **standardised relative index** (never in
cells/mL) if check C is *pattern consistent* and, for *Prochlorococcus*, check
A is not *not supported*. Otherwise the group is excluded from the
reconstruction, and the evaluation itself is reported.

## What is reported whatever the outcome

Number of ship records, matched and unmatched; metrics with intervals per
group, and per region and season where the sample allows; results under both
coefficient mappings (operational and ATBD); the day-offset distribution; a
sensitivity table for the 3×3 window and the 5-pixel minimum.

## Sign-off

| Name | Role | Date | Agreed |
| --- | --- | --- | --- |
| Ramachandra TejSrikar Reddy | Team | 2026-11-16 | Yes. Fixed the rule and delegated its drafting |
| Aryan Sinha | Team | — | Did not review (left the decision to Ramachandra) |
| Dr. M. Suganiya | Guide | — | Did not review |

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-11-10 | First draft | — |
| 2026-11-16 | Thresholds accepted unchanged; fallback checks A–C added; status set to binding | SeaBASS search found no PACE-era cell counts in the box but did find HPLC pigments and phycoerythrin fluorescence. Added before any comparison was made. |
| 2026-11-16 | Station-day defined as one UTC date in one 4.6 km grid pixel | The two EKAMSAT 2024 releases name stations differently. Defined before any comparison. |
| 2026-11-16 | Data correction, not a rule change: BCO-DMO records are de-duplicated against PANGAEA per station (same date, within 0.05°) instead of per sample (also within 0.5 m depth) | Found while checking the first run of check C. PANGAEA gives TT053/TT054 depths in whole metres (0–3 m) where BCO-DMO has 2.5–2.8 m, so 1,433 duplicate 1995 records survived and appeared as extra stations under four invented cruise labels. Check C was rerun once. Both runs are recorded in `docs/open_items.md` (item 17); the decision uses the corrected run. Checks A and B do not use these data and are unchanged. |
