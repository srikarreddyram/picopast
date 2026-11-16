# Ship-count inventory (open archives only)

Generated 2026-11-10 by `picopast inventory`. Box: 0-25N, 50-100E. Eras: pre_viirs (before 2012), viirs_era (2012 to 2024-03-05, serves RQ4), pace_era (from 2024-03-05, serves the RQ1 trust gate).

## BCO-DMO global compilation (dataset 793451)

2792 records with at least one group count, downloaded to `data/insitu/bcodmo_global_compilation_box.csv`.

| era | records | surface (<= 10 m) | surface pro | surface syn | surface peuk | surface in bay_of_bengal | surface in arabian_sea_central_east | years |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pre_viirs | 2785 | 347 | 309 | 277 | 78 | 8 | 157 | 1993-2005 |
| undated | 7 | 4 | 0 | 0 | 4 | 0 | 0 | - |

Other BCO-DMO datasets that mention the groups and touch the box:

| dataset_id | title |
| --- | --- |
| bcodmo_dataset_526852 | [JeDI] - Jellyfish Database Initiative: Global records on gelatinous zooplankton for the past 200 years, collected from global sources and literature (Trophic BATS project) (Plankton Community Composition and Trophic Interactions as Modifiers of Carbon Export in the Sargasso Sea 	) |
| bcodmo_dataset_783463 | [Picoeukaryotic phytoplankton observations] - Picoeukaryotic phytoplankton observations from available public repositories and primary sources from from 1988-2007 (Convergence: RAISE: Linking the adaptive dynamics of plankton with emergent global ocean biogeochemistry) |
| bcodmo_dataset_793451 | [Global Prochlorococcus Synechococcus and picoeukaryotic phytoplankton] - Global observations Prochlorococcus, Synechococcus, and picoeukaryotic phytoplankton with ancillary environmental data from 1987 to 2008 (Convergence: RAISE: Linking the adaptive dynamics of plankton with emergent global ocean biogeochemistry) |
| bcodmo_dataset_793847 | [Pro Syn Peuk CC Surface] - Prochlorococcus, Synechococcus, and picoeukaryotic phytoplankton for four future climatic scenarios and five Earth System Models (cell/ml) in a global 1x1 grid for the ocean surface (50m) (Convergence: RAISE: Linking the adaptive dynamics of plankton with emergent global ocean biogeochemistry) |
| bcodmo_dataset_806510 | [Pacific Nitrite Oxidoreductase] - Nitrite Oxidoreductase targeted metaproteomics from R/V Kilo Moana cruise KM1128 and R/V Falkor cruise FK160115 in the Central Pacific Ocean in 2011 and 2016 (Connecting Trace Elements and Metalloenzymes Across Marine Biogeochemical Gradients) |
| bcodmo_dataset_924554_v1 | [Plankton size spectra] - Plankton size spectra compiled from projects CCE LTER, HOT, CRD FluZIE, BLOOFINZ-GoM, and SalpPOOP in the California Current Ecosystem, North Pacific subtropical gyre, Costa Rica Dome, Gulf of Mexico, and Southern Ocean subtropical front from 2004-2018 (Collaborative Research: Quantifying trophic roles and food web ecology of salp blooms of the Chatham Rise) |

## PANGAEA

184 datasets in the box report at least one group (`data/insitu/pangaea_datasets_box.csv`). A dataset is usually one station profile.

| era | datasets | with_pro | with_syn | with_peuk | first | last |
| --- | --- | --- | --- | --- | --- | --- |
| pre_viirs | 183 | 183 | 182 | 182 | 1991-05-14 | 2002-07-08 |
| viirs_era | 1 | 1 | 0 | 0 | 2016-04-30 | 2016-04-30 |

Datasets from 2012 onward:

| doi | title | time_start | lat | lon |
| --- | --- | --- | --- | --- |
| 10.1594/PANGAEA.982433 | A validation data set of phytoplankton pigment concentrations and phytoplankton groups measured on water samples collected from various expeditions | 2016-04-30T10:59:00 | 5.981769999999997 | 8.047416500000011 |

## SeaBASS (via NASA CMR)

CMR lists every SeaBASS file with its name, date and footprint, but not the fields inside it, and SeaBASS collection footprints are all global. So this section is a shortlist: **nothing here is a confirmed cell count until the file header is read, which needs the Earthdata login.**

Experiments whose title or summary names the region (14). Every one should be opened after login:

| collection | title | time_start | time_end | files_listed | files_in_box | files_no_footprint | candidate_files |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BOA | Measurements used to develop the Bio-Optical Algorithm (BOA) | 1991-05-17 |  | 3 | 1 | 2 | 0 |
| JGOFS_BOFS | Joint Global Ocean Flux Study (JGOFS) Arabian Sea measurements - Biogeochemical Ocean Flux Study (BOFS) | 1991-06-19 |  | 40 | 0 | 0 | 0 |
| JGOFS_Arabian_Sea | Joint Global Ocean Flux Study (JGOFS) Arabian Sea measurements | 1994-03-15 |  | 1053 | 1043 | 0 | 0 |
| INDOEX | India Ocean Experiment | 1999-01-16 |  | 721 | 272 | 67 | 0 |
| THAILAND | Optical measurements in the Gulf of Thailand and the Andaman Sea in 2003 | 2003-10-09 |  | 1 | 0 | 1 | 0 |
| EGEE3 | Gulf of Guinea climate and ocean circulation study (EGEE) project - Gulf of Guinea off of west-central Africa | 2006-05-29 |  | 66 | 0 | 1 | 0 |
| AMMA | African Monsoon Multidisciplinary Analyses (AMMA) | 2006-06-25 |  | 140 | 0 | 0 | 0 |
| EGEE5 | Gulf of Guinea climate and ocean circulation study (EGEES) | 2007-06-06 |  | 13 | 0 | 1 | 0 |
| OCEANSAT-2 | Monsoonal Arabian Sea measurements to validate OceanSat-2 | 2009-02-10 |  | 22 | 20 | 2 | 0 |
| ArabianSea_2011 | Measurements made in the monsoonal Arabian Sea in 2011 | 2011-03-07 |  | 134 | 134 | 0 | 0 |
| ASIRI | Air-Sea Interaction Research Initiative (ASIRI), Bay of Bengal | 2013-11-29 |  | 88 | 88 | 0 | 0 |
| EKAMSAT_Pilot_ASTRAL | Enhancing Knowledge of the Arabian Sea Marine environment through Science and Advanced Training (EKAMSAT) | 2023-06-10 |  | 100 | 98 | 2 | 0 |
| PVST_NORTHERN_INDIAN_OCEAN | Arabian Sea bio-optical and biogeochemical data for ocean color validation | 2024-01-01 |  | 590 | 471 | 119 | 1 |
| OCEANX_2025 | OCEANX_2025 | 2025-01-01 | 2025-12-31 | 1738 | 47 | 736 | 1 |

Other experiments with files named like flow cytometry (`*fcm*`, `*cyto*`) whose footprint lies in the box (0):

_none_

2 candidate files in total (`data/insitu/seabass_candidate_files.csv`).

| collection | files | first | last | example |
| --- | --- | --- | --- | --- |
| OCEANX_2025 | 1 | 2025-01-09 | 2025-01-09 | 8193875b42_OCEANX_2025_OXR20250108_CLS.sb |
| PVST_NORTHERN_INDIAN_OCEAN | 1 | 2024-04-28 | 2024-04-28 | 5ea44233a1_EKAMSAT_PVST_NORTHERN_INDIAN_OCEAN_CLS.sb |

3 collections could not be listed (see the `error` column).
