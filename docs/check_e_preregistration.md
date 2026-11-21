# Check E: new, independent stations for the Route D record — rules fixed before testing

Status: **binding from 2026-11-21**, committed before any value from the new
stations is compared with Route D, MOANA or chlorophyll. Changes go in the
change log at the end, with the reason; every earlier run stays reported.

## Why

Check C tested Route D against 1991–2005 ship counts with a same-calendar-month
field. That tests the *pattern* (where and in which season each group is high),
not the 2012–2026 *record* itself. Check E adds stations that check C did not
use, and, for stations dated 2012 or later, tests the record's own 8-day values.

## What was already seen (disclosed)

- Searched before this file: SeaBASS (all OB_DAAC collections, files named for
  HPLC or flow cytometry, in the box, all years), PANGAEA, BCO-DMO ERDDAP and
  gVOD (PANGAEA.915758). Only metadata and column names were read: gVOD holds
  no picophytoplankton counts in the box; Tarran (2004, PANGAEA 1358xx) reports
  only total cyanobacteria; BCO-DMO holds no in-box dataset beyond the
  compilation check C already uses.
- The 27 PACE-era HPLC station-days from EKAMSAT 2024 and OceanX 2025 were
  used in check A against MOANA (`trust_gate_fallback.md`); their divinyl
  chlorophyll a ranged only 0.023–0.054 mg m⁻³. They have never been compared
  with Route D, and they are included here.

## Data

- **E-counts:** flow-cytometry counts of *Prochlorococcus* or *Synechococcus*
  in the box that are not check C stations (same UTC date and within 0.05° of
  a check C station = duplicate), any year, from any source found later (each
  logged here with its source before use).
- **E-pigments:** HPLC from SeaBASS (and any other HPLC source logged here),
  in the box, any year. Divinyl chlorophyll a (DVChl a, made only by
  *Prochlorococcus*) for *Prochlorococcus*; zeaxanthin for *Synechococcus*.
  Zeaxanthin is also made by *Prochlorococcus*, so the *Synechococcus*
  pigment test is **descriptive only and never gives a verdict**.
- For both: shallowest sample ≤ 10 m; one row per station-day (one UTC date ×
  one 4.6 km pixel; median of merged samples, as `insitu/pigments.py`);
  pigment or count > 0; same offshore rule as check C.

## Tests

**E1, pattern (all dates).** The Route D D2 same-calendar-month field from
PACE-era inputs, exactly as in check C (3×3 matching, same month). Also scored:
the chlorophyll-alone field (sign fixed as in check C: *Prochlorococcus*
negative, *Synechococcus* positive) and MOANA's field.

**E2, record (stations dated 2012-01-01 or later).** The SNPP Route D record
for the station's 8-day period: median of the valid pixels in the 3×3 window
around the station. **All valid pixels are used** (primary), because the
reliable flag drops the warmest water and the test should cover the conditions
the record is used in; reliable pixels only is reported as a sensitivity.
Same period only, no borrowing from neighbouring periods. Secondary
comparisons on the same stations: VIIRS 8-day chlorophyll (sign as above), the
earlier MOANA-learned Route A record (exploratory) for the same period, and
the partial ρ of the record given the same-month D2 field (does the 8-day
record carry information beyond its own climatology?).

**Statistics and verdict (both tests),** unchanged from check C: Spearman ρ;
95% interval from resampling cruises (2,000 draws); *pattern consistent* if the
lower end is above 0 with ≥ 30 stations from ≥ 3 cruises; *not consistent* if
the upper end is below 0.30; otherwise *inconclusive*; fewer stations or
cruises = *not enough data*, reported descriptively without a verdict. Paired
differences (record − chlorophyll, record − Route A) with cruise-resampled
intervals.

## Decision

- ***Prochlorococcus*:** E2 with DVChl a is the primary test of the record (E1
  if E2 has not enough data). *Pattern consistent* supports the Route D
  record as a relative index; *not consistent* downgrades it to exploratory;
  inconclusive or not enough data changes nothing and is reported as such.
- ***Synechococcus*:** only E-counts can change its tier, by the same rule.
  The zeaxanthin results are reported, never used for a verdict.
- **Picoeukaryotes:** if E-counts hold any, the MOANA-based record is scored
  the same way and reported; no change of tier without a pre-registered rule.

## Known limits, stated with every result

- Pigment per cell changes with light and nutrients (cells near the surface
  carry less DVChl a), so pigment tests are about ranks, not amounts.
- Stations cluster in a few cruises; the interval resamples cruises, so a
  test with few cruises will be wide.
- Same-period matching loses stations under cloud; no gap filling.

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-11-21 | First version, binding | Written before any new station was compared with a model. |
| 2026-11-21 | Sources logged (E-pigments): SeaBASS HPLC from CLIVAR I8S/I9N 2007, Tara Oceans 2010 (Abu Dhabi–Muscat–Mumbai–Malé legs), GO-SHIP I07N 2018, EKAMSAT pilot 2023, GO-BGC RR2307 2023, Bio-GO-SHIP I09N 2025, PVST SO-PACE TN444 2025, OceanX 2025, plus the earlier EKAMSAT 2024 and OceanX files. 149 surface station-days with DVChl a in the box: 127 dated 2012 or later from 7 cruises, 22 earlier. E-counts: none found (PANGAEA, BCO-DMO, gVOD, SeaBASS hold no new in-box *Prochlorococcus* or *Synechococcus* counts). | Full SeaBASS scan finished; only positions, dates and depths were tabulated. |
| 2026-11-21 | HPLC reader extended to three SeaBASS layouts (a `date` column; `year`/`month`/`day` columns; older one-station files with position, date and depth in the header and `/end_header@`); field names matched case-insensitively; rows without an `hplc_gsfc_id` are no longer merged. The PACE-era station-day count is unchanged (27). | Needed to read the new files. While checking the CLIVAR layout one sample's pigment values were displayed; nothing was compared with any model. |
