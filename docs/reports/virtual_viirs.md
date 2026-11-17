# Virtual VIIRS vs real VIIRS on the same day

Generated 2026-11-17 by `picopast virtual-viirs-check`. PACE OCI daily Rrs, band-averaged with each sensor's M1–M5 response (NOAA NCC releases, `picopast.ingest.rsr`), against that sensor's own daily L3m Rrs on the same 4.6 km pixels (offshore, both valid; up to 50,000 sampled pixels a day). The residual mixes the simple band-averaging, calibration and atmospheric-correction differences between missions, and the 0–1 h gap between overpasses. It is documented here, not corrected.

## SNPP

29 days (2025-01-01 to 2025-01-31), 1,316,982 pixel pairs.

| band (nm) | median real Rrs (sr⁻¹) | bias (virtual − real) | bias / median | RMSD | median ratio | r |
| --- | --- | --- | --- | --- | --- | --- |
| 410 | 0.00406 | +0.00010 | +2.6% | 0.00120 | 1.018 | 0.868 |
| 443 | 0.00377 | +0.00001 | +0.4% | 0.00092 | 1.000 | 0.878 |
| 486 | 0.00367 | -0.00003 | -0.9% | 0.00070 | 0.993 | 0.918 |
| 551 | 0.00175 | +0.00005 | +2.6% | 0.00049 | 1.021 | 0.965 |
| 671 | 0.00027 | +0.00009 | +33.6% | 0.00024 | 1.300 | 0.976 |

Blue-green ratio log10(Rrs443/Rrs551), the core of chlorophyll and of the models' strongest inputs: bias -0.006, r 0.928.

## NOAA-20

29 days (2025-01-01 to 2025-01-31), 1,325,579 pixel pairs.

| band (nm) | median real Rrs (sr⁻¹) | bias (virtual − real) | bias / median | RMSD | median ratio | r |
| --- | --- | --- | --- | --- | --- | --- |
| 411 | 0.00338 | +0.00084 | +25.0% | 0.00152 | 1.208 | 0.855 |
| 445 | 0.00371 | +0.00011 | +2.9% | 0.00097 | 1.017 | 0.870 |
| 489 | 0.00366 | -0.00004 | -1.1% | 0.00074 | 0.982 | 0.916 |
| 556 | 0.00171 | -0.00002 | -1.4% | 0.00052 | 0.973 | 0.967 |
| 667 | 0.00030 | +0.00005 | +16.0% | 0.00025 | 1.112 | 0.975 |

Blue-green ratio log10(Rrs445/Rrs556), the core of chlorophyll and of the models' strongest inputs: bias +0.026, r 0.938.
