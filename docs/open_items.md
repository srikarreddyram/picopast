# Open items

Updated 10 November 2026 (week 1). Ordered by how soon each one blocks work.

## Blocking now

1. **Earthdata login.** Run `.venv/bin/picopast login` once on this machine.
   Until then: no PACE or VIIRS downloads, no Atlantic MOANA check, no SeaBASS
   files. Also create the Copernicus Marine account (mixed-layer depth, week 2).
2. **First commands after login**, in this order:
   - `picopast moana-check 2025-01-15` — our MOANA vs NASA's Atlantic file.
     This is the first live run of the download code, so expect to fix
     file-layout details.
   - `picopast fetch pace_rrs 2025-01-01 2025-01-31` and the VIIRS products for
     the same month — the one-month test cube and the monsoon-coverage check.
   - Open every regional SeaBASS collection listed in
     `docs/reports/data_inventory.md` and read the file headers for cell counts.

## Decisions for the team and guide (before week 3)

3. **Sign the trust-gate rule** in `docs/trust_gate_preregistration.md`. The
   thresholds in it are proposals. It must be fixed before any matchup is run.
4. **Region boxes.** The PRD names the regions but gives no coordinates. The
   boxes in `config/picopast.yaml` are proposals: Bay of Bengal 5–22°N,
   80–95°E; central/eastern Arabian Sea 8–22°N, 62–75°E.
5. **Ship-count shortfall.** See the inventory report. If the 2024+ count
   stays below about 30 matchups per group after SeaBASS files and published
   tables are read, the PRD's fallback applies (pattern comparison instead of
   absolute errors) and RQ1 and RQ4 must be worded accordingly. Decide early
   whether that changes the paper's headline.
6. **Published tables.** The two 2026 papers named in the PRD could not be
   read from here (publisher pages refuse automated access or need a login).
   Someone with library access needs to check each for station-level tables
   and note the cruise dates. Record every table in a registry file with its
   source, page and who typed it.

## Corrections to the PRD found while building

7. **Prochlorococcus is regressed in linear cells/mL, not log10** (PRD section
   6 step 2 says log10 for all groups). Negative retrievals are possible and
   are flagged. Details in `docs/moana_notes.md`.
8. **The ATBD equations and NASA's operational coefficients disagree** for two
   terms. We default to the operational table and can run both.
9. **PACE Rrs collection name.** In version 3.2 the hyperspectral Rrs is in
   `PACE_OCI_L3M_AOP`; there is no `PACE_OCI_L3M_RRS` collection in CMR.
10. **The RSR table link in the PRD is dead.**
    `oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables` now redirects to a
    general NASA page. The VIIRS SNPP and NOAA-20 response tables need a new
    source before week 4 (virtual VIIRS). `picopast.virtual_viirs` takes any
    table of wavelength by band.
11. **DMI lags.** NOAA PSL's HadISST-based Dipole Mode Index currently ends in
    May 2026. The last months of the record will have no DMI unless another
    source is chosen.
12. **Download volume.** A global PACE AOP daily 4 km file is 0.5–0.8 GB, so
    the PACE era is about 970 files and several hundred GB if kept whole. The
    pipeline reads only the box from each file over HTTPS; how fast that is
    will be known after the first logged-in run.

13. **OISST server limits.** NOAA PSL's OPeNDAP service started answering
    HTTP 429 (too many requests) after two years of box-only reads. The
    download now stops cleanly on 429, and `picopast fetch-open --oisst-method
    file` fetches each global yearly file (about 0.5 GB), keeps the box and
    deletes the rest. Both methods give identical grids.

## Not started (by design)

- Zarr cube builder: waits for the first real PACE and VIIRS files so it is
  written against actual file layouts.
- Routes A, B, C, baselines, intervals, SHAP (weeks 4–5).
- Reconstruction, dual-sensor test, climate statistics beyond heatwave
  detection (weeks 6–7).
- Scopus / Web of Science / Shodhganga search log (T01): needs institutional
  access.
