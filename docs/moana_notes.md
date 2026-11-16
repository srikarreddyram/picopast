# MOANA implementation notes

What we learned from the ATBD (v1.2, 16 Dec 2024), NASA's operational tables
and the EPFT-UP project while writing `picopast.moana`. Each point matters for
how results are described in the paper.

## The algorithm as NASA runs it

1. Interpolate Rrs to 414–660 nm in 2 nm steps (124 bands). The ATBD abstract
   says 395–705 nm and its introduction says 400–600 nm; the shipped PCA table
   is 414–660 nm, which is also what the ATBD's method section states.
2. Standardise each spectrum across wavelength (mean 0, sample standard
   deviation 1). The scores are a plain dot product with the PCA loadings; no
   training-set mean is subtracted.
3. Use the first 17 principal components (45 are stored).
4. Regressions:
   - picoeukaryotes: log10(cells/mL) on U2, U3, U4, U5, U10, U13, U15;
   - *Synechococcus*: log10(cells/mL) on U1, U2, U3, U5, U8–U12 and one more;
   - *Prochlorococcus*: **linear cells/mL** on log10(SST), U1, U2, U6 and one more.

## Points that differ from the PRD summary

- **Prochlorococcus is not regressed in log10 space.** Its regression is
  linear in cells/mL, so the raw output can be negative. NASA clamps negatives
  to zero; we keep the value and set a flag, and treat negatives as invalid in
  matchups. The share of negative retrievals in the North Indian Ocean is
  itself a result worth reporting.
- **The ATBD and the operational table disagree** on the principal component
  of the last term for two groups: *Prochlorococcus* U17 (ATBD) vs U7 (table);
  *Synechococcus* U13 (ATBD) vs U16 (table). Coefficient values agree. EPFT-UP
  reports that the shipping PACE product follows the table. We default to the
  table (`pc_mapping="operational"`) and can run the ATBD version for
  comparison. Which one the authors intended is unknown to us.
- **SST is in degrees Celsius**, from GHRSST CMC 0.1° in NASA's processing.
  CMC files are in kelvin, so the pipeline converts; `run_moana` refuses input
  that looks like kelvin.

## Verification status

| Check | Status |
| --- | --- |
| NASA tables match published checksums | Done (pinned in code) |
| Our code vs an independent reference written from the ATBD equations | Done, agreement to 1e-9 relative (`tests/test_moana.py`) |
| Our code vs NASA's Atlantic L4 files, same dates | **Not done: needs the Earthdata login.** `picopast moana-check 2025-01-15` |
| Our code vs EPFT-UP's numbers | Not done. Running their code was not possible in this session; the Atlantic check supersedes it |

The Atlantic check is the one that counts. Small differences from NASA's L4
product are possible even with correct code, because we do not know at which
processing stage NASA applies MOANA (the ATBD mentions both Level-2 and
Level-3 Rrs as input), and running on mapped Rrs is not identical to mapping a
retrieval. EPFT-UP reports a median difference of about 0.006 in log10 for
*Synechococcus* on one daily file pair, so we expect close agreement. The size
of the difference is what `moana-check` measures and what the paper must state.

## Sources

- MOANA ATBD v1.2: <https://oceandata.sci.gsfc.nasa.gov/files/atbd/atbd-obdaac-phytoplankton-groups-MOANA.pdf>
- Lange et al. 2020, *Optics Express* 28(18), doi:10.1364/OE.398127
- EPFT-UP: <https://github.com/ocean-colour/EPFT-UP> (not peer reviewed)
- NASA OCSSW tables: `https://oceandata.sci.gsfc.nasa.gov/manifest/tags/T2023.31/common/`
