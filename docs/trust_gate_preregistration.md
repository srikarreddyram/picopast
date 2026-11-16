# Trust gate pre-registration (RQ1)

**Status: DRAFT, not yet binding.** Drafted 10 November 2026, before any
satellite-to-ship matchup has been computed. It becomes binding when both team
members and the guide have read it and the sign-off block at the bottom is
filled in. That must happen before the week 3 matchups are run. After sign-off
this file and the `trust_gate` and `matchup` blocks of `config/picopast.yaml`
are frozen; any later change is recorded in the change log below with a reason.

Every number marked *proposed* is a starting suggestion for the team to accept
or change now. None of them comes from the PRD.

## Question

How accurate is our MOANA run in the North Indian Ocean, for each of
*Prochlorococcus*, *Synechococcus* and picoeukaryotes?

## Data

- **Satellite:** our own MOANA run on PACE OCI L3M Rrs (daily, native 4.6 km
  grid) with GHRSST CMC 0.1° SST, operational coefficient mapping.
- **Ship:** surface cell counts from open archives and published tables only,
  dated 5 March 2024 or later, inside 0–25°N, 50–100°E.
- **Surface** means sample depth ≤ 10 m (*proposed*). Where several depths
  ≤ 10 m exist at one station, the shallowest is used.
- **Excluded before matching:** stations within 25 km of the coast; counts
  that are zero, negative or flagged by the archive.

## Matchup

1. Same calendar day (UTC) as the ship sample. If that fails, the day before,
   then the day after.
2. Value = median of the 3×3 pixel block centred on the station.
3. A block counts only with at least 5 of 9 valid pixels (*proposed*).
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
| **Usable label** | All of: log-RMSE ≤ 0.50; \|bias\| ≤ 0.30; lower end of the ρ interval ≥ 0.30; at least 50% within a factor of 2 (*all proposed*). |
| **Relative index only** | Not a usable label, but ρ ≥ 0.40 and the lower end of its interval is above 0 (*proposed*). The group is modelled as a standardised index, not as cells/mL. |
| **Excluded** | Neither of the above. The group is dropped from the reconstruction and the evaluation itself is reported. |

A group needs matchups from at least 2 cruises for a cruise bootstrap to mean
anything; with one cruise the verdict is reported as provisional (*proposed*).

The rule is implemented in `picopast.trust_gate.verdict` and reads its
thresholds from `config/picopast.yaml`.

## What is reported whatever the outcome

Number of ship records, matched and unmatched; metrics with intervals per
group, and per region and season where the sample allows; results under both
coefficient mappings (operational and ATBD); the day-offset distribution; a
sensitivity table for the 3×3 window and the 5-pixel minimum.

## Sign-off

| Name | Role | Date | Agreed |
| --- | --- | --- | --- |
| Aryan Sinha | Team | | |
| Ramachandra TejSrikar Reddy | Team | | |
| Dr. M. Suganiya | Guide | | |

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-11-10 | First draft | — |
