# Methods research: how to get a definitive answer on heatwaves, cyclones and phytoplankton

22 November 2026. The question this note answers: the main test (does a marine heatwave
change the phytoplankton bloom after a cyclone?) ended *no detectable difference*, with
every estimate leaning towards a weaker bloom. Which techniques, data and simulations
could settle it? Every source below came from searches run today. For most, only the
abstract or an excerpt was read; the reference list says where.

## Why the answer is not yet definitive

1. **Too few storms crossed heatwave water.** 46 in all (39 since 2012, 7 before), and
   cloud left a valid bloom at only about 30 of them.
2. **The satellite sees only the surface.** A surface "bloom" after a storm is often
   chlorophyll mixed up from the deep chlorophyll maximum, not new growth (Typhoon Trami
   study, Biogeosciences 2021). Heatwaves change that deep layer.
3. **The average may hide the answer.** Pre-existing stratification blocks the bloom
   unless the storm is strong or slow enough to break through. Cyclone Phailin was
   intense but gave a small bloom over strong stratification; the weaker Thane gave a
   large one (North Indian Ocean review, arXiv 2012.04384). A heatwave may matter for
   weak or fast storms and not for strong or slow ones. Averaging over both dilutes it.
4. **Measurement noise.** 8-day composites smear the timing, clouds remove storms, and
   OISST (infrared plus in situ) under-reads cold wakes. Infrared-mostly products miss
   about 0.4 °C of hurricane cooling (Remote Sensing 2019).

## The upgrades, by priority

### Tier 1: sharper tests on data we have or are downloading (days)

| # | Technique | What it answers | Evidence it works | Cost |
| --- | --- | --- | --- | --- |
| 1 | **Mechanism test: heatwave × storm mixing power** | Does a heatwave weaken the bloom only when the storm is too weak or fast to break the warm layer? A pre-registered interaction with a mixing index (wind² / translation speed, an energy-input proxy) | Stratification and translation speed control blooms in this basin (Kuttippurath 2021; review 2012.04384) | Existing tables |
| 2 | **More storms: the South Indian Ocean** (0–35°S, 30–120°E) | The same test on 360 more storms, about 4× the sample, same satellites and rules | IBTrACS SI list; sea-temperature download already running | About 1 day of downloads, then the existing pipeline |
| 3 | **Entrainment or new growth: timing** | Mixed-up chlorophyll appears at once (first 8 days); growth peaks days to weeks later (second period). Test the heatwave effect on each separately | Trami study; BGC-Argo cases in the Bay of Bengal show an early mixing pulse, then a later bloom | Existing tables |
| 4 | **Storm-minus-no-storm controls (difference-in-differences)** | Compare each storm disc with the same water, same heatwave status, same days, but no storm. This removes heatwave decay and seasonal rebound, as the post-hoc cooling check did | Background-variability composites (Southampton thesis 2020); coastal vs offshore contrasts | Existing data |
| 5 | **A decision rule that can end definitively both ways** | Fix in advance what counts as "too small to matter", for example a third of the typical bloom. The result is then a confirmed weakening, a confirmed strengthening, or a confirmed "no meaningful effect" (an equivalence test). A stricter interval (99%) allows for this being the third look | Standard practice (equivalence testing; corrections for repeated looks) | None |
| 6 | **Power simulation** | How many heatwave storms are needed for rule 5 to decide? Monte Carlo on our own bootstrap variances | Our variance estimates | Minutes |

A rough power calculation from today's numbers:
- *Now:* the pooled standard error is about 0.030 log10, with about 30 heatwave storms
  that have a valid bloom.
- *To detect* a −0.05 effect at 99% with 80% power, it must fall to about 0.015. That
  takes roughly 4 times the heatwave storms, about 120–130.
- *The South Indian Ocean* could supply about 70. If heatwaves are as common there after
  2012 as in the north (30% of points), that gets close.

### Tier 2: simulations (1–2 weeks)

**7. Virtual experiments with a 1-D ocean column model.** This is the "simulation" core.
- *The model:* Price–Weller–Pinkel (PWP), the bulk mixed-layer model used to explain a
  typhoon bloom (Lin et al. 2003, GRL), extended with a simple nutrient–phytoplankton
  (NPZD) biology. PWP-based 1-D physical–biological models already exist (Marine Ecology
  Progress Series 2011).
- *Code:* Python ports exist (ParagGawai/PWP_Model; RomiP/pyPWP, based on Earle Wilson's
  code). The original Fortran is on James Price's WHOI page. GOTM with FABM is the
  heavier alternative, already used in 1-D to train biogeochemistry emulators.
- *Forcing:* each real storm's winds, rebuilt from IBTrACS with the Holland parametric
  profile (maximum wind, radius of maximum wind, shape B). It can be calibrated with
  best-track data (Nederhoff et al. 2019, NHESS); its limits are discussed in
  Willoughby & Rahn (2004).
- *Initial profiles:*
  - temperature and salinity: Argo (the 2026 paired Argo dataset below, or Argo GDAC)
    or the Copernicus reanalysis;
  - nitrate: World Ocean Atlas depth profiles;
  - a deep chlorophyll maximum, set from BGC-Argo climatology.
- *The experiment:* run every storm twice.
  - once over its real pre-storm profile;
  - once with the heatwave's warm anomaly removed, all else equal.

  The paired difference is the heatwave's causal effect in the model world, on mixing
  depth, nutrient supply, cooling and surface chlorophyll. That cannot be measured from
  observations alone. A Hurricane Michael study used this idea, comparing idealized
  mixed versus stratified profiles in a 1-D model (Nature Communications 2020).
- *Validation:* compare the model's cooling and mixed-layer deepening with OISST,
  microwave SST and Argo before and after each storm.
- *Ensemble:* repeat over uncertain settings (drag coefficient, nutrient profile, light)
  to give a range, not a single run.
- *Limit:* one column has no horizontal advection or eddies. Eddies matter here (Hudhud
  was modelled in 3-D with ROMS-NPZD). The 1-D runs isolate the vertical mechanism, and
  the satellite statistics give the real-world total.

**8. Microwave sea temperature for cooling.** REMSS MW_IR_OI and MW_OI v5.0 (NASA
PO.DAAC) see through cloud. They fix the cooling under-read and give the PWP runs a
cleaner target.

**9. Cloud gap-filling for daily chlorophyll.**
- *Methods:* DINEOF fills cloudy pixels from the main patterns of the time series.
  ExDINEOF filled daily chlorophyll during Cyclone Hudhud in the Bay of Bengal better
  than older methods (Wang et al. 2023, Frontiers). DINEOF was also used for Arabian
  Sea cyclone blooms (ISPRS 2020). DINCAE is a neural-network alternative; a 2026 GMD
  comparison found each has strengths, with code on Zenodo.
- *What it gives:* more storms with a bloom, and daily timing, which feeds test 3.

### Tier 3: below the surface (needs the Copernicus login for 10)

**10. Copernicus reanalyses as a 3-D laboratory.**
- *Products:* GLORYS (physics: mixed-layer depth, temperature and salinity profiles)
  and the global biogeochemistry hindcast GLOBAL_MULTIYEAR_BGC_001_029 (PISCES; nitrate
  and chlorophyll at depth, 1/4°, from 1993). BIORYS4, which assimilates ocean colour,
  is announced.
- *What it answers:* for every storm since 1993, whether heatwave water had a deeper
  nitracline and deep chlorophyll maximum, and whether chlorophyll summed over depth
  rises, or only the surface value (entrainment).

**11. BGC-Argo floats under storms.** Floats measured chlorophyll and nitrate through
Hudhud, Vardah, Leher and Madi in the Bay of Bengal, and Kyarr and Maha in the Arabian
Sea. They are too few for statistics, but they are direct case studies of what happens
under the satellite's view. Heatwaves are known to reshape phytoplankton's vertical
structure: a global BGC-Argo study (Communications Earth & Environment 2025); spring
blooms suppressed in the Mediterranean (GRL 2024).

**12. The global paired Argo dataset (2000–2024).** Pre- and post-storm temperature and
salinity profiles for IBTrACS storms, about 3.16 million pairs (ESSD preprint 2026,
figshare). It gives the real pre-storm stratification under heatwave and normal water,
for the mechanism and for the PWP starting profiles.

### Tier 4: exploratory machine learning

**13. Causal forests** (Wager & Athey 2018) estimate how the heatwave effect varies with
storm and ocean conditions, and which conditions matter most. They have been used for
weather-shock and aerosol–cloud questions. They are exploratory only and come after the
pre-registered tests; small samples make them unstable (compare the Bay of Bengal
emulator "too rare to learn", arXiv 2026).

## Recommended order

1. **Pre-register tests 1, 3, 4 and 5 now**, with the power simulation (6) attached, and
   finish the South Indian Ocean downloads (2). One final confirmatory test then
   combines everything. No later data changes its verdict.
2. **Build the PWP + NPZD virtual experiments (7)**, validated against microwave SST (8)
   and Argo (12). Pre-register them as model experiments before running.
3. **Copernicus (10)** as soon as the login exists; BGC-Argo case studies (11) alongside.
4. **Gap-filling (9) and causal forests (13)** if time allows.

## What this costs and what it buys

| | Days | Decides |
| --- | --- | --- |
| Tier 1 | 2–4 (mostly download time) | The average effect, *and* whether it depends on mixing power: a definitive answer either way, or a precisely bounded one |
| Tier 2 | 7–10 | Why: the heatwave effect isolated in a controlled model, with the storm held fixed |
| Tier 3 | 3–5 once logged in | Whether the surface bloom is real growth or mixed-up chlorophyll, and how a heatwave changes the deep layer |

## Sources (as found on 22 November 2026; most read at abstract or excerpt level)

Simulation and models
- Lin et al. 2003, GRL — [mixed-layer (PWP) model explains a typhoon bloom](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003GL017141)
- PWP-based 1-D physical–biological model, Marine Ecology Progress Series 427 (2011) — [PDF](https://digitalcommons.library.umaine.edu/cgi/viewcontent.cgi?article=1040&context=sms_facpub)
- Compounding impact of severe weather events fuels marine heatwave in the coastal ocean, Nature Communications 2020 (Hurricane Michael; idealized 1-D model) — [article](https://www.nature.com/articles/s41467-020-18339-2)
- Rainfall-enhanced blooming in typhoon wakes — [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4992858/)
- Cyclone Hudhud–eddy bloom with a coupled ROMS-NPZD model — [article](https://www.sciencedirect.com/science/article/abs/pii/S0079661121001178)
- Resolved tropical cyclones trigger CO2 uptake and blooms in an Earth system model (ICON-HAMOCC, 2025) — [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12718355/)
- GOTM–FABM–ERSEM 1-D runs as emulator training — [arXiv 2606.27168](https://arxiv.org/pdf/2606.27168)
- PWP code: [ParagGawai/PWP_Model](https://github.com/ParagGawai/PWP_Model); James Price's [upper-ocean projects and Fortran PWP](https://www2.whoi.edu/staff/jprice/upper-ocean-projects)
- Holland profile calibrated with best-track data, Nederhoff et al. 2019 — [NHESS](https://nhess.copernicus.org/articles/19/2359/2019/); its limits, Willoughby & Rahn 2004 — [PDF](https://faculty.fiu.edu/~willough/PUBS/HEW&RAHN_04.pdf)

Observations and data
- Global paired Argo dataset for tropical cyclones (2000–2024), ESSD preprint 2026 — [preprint](https://essd.copernicus.org/preprints/essd-2026-520/)
- Spatio-temporal methods for subsurface responses to tropical cyclones — [arXiv 2012.15130](https://arxiv.org/pdf/2012.15130)
- REMSS microwave SST on NASA PO.DAAC — [MW_IR_OI v5.0](https://podaac.jpl.nasa.gov/dataset/MW_IR_OI-REMSS-L4-GLOB-v5.0), [MW_OI v5.0](https://podaac.jpl.nasa.gov/dataset/MW_OI-REMSS-L4-GLOB-v5.0)
- Infrared-mostly SST under-reads hurricane cooling by about 0.4 °C, Remote Sensing 2019 — [article](https://www.mdpi.com/2072-4292/11/20/2368)
- Copernicus biogeochemistry hindcast GLOBAL_MULTIYEAR_BGC_001_029 — [product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description); BIORYS4 — [EGU 2025](https://meetingorganizer.copernicus.org/EGU25/EGU25-6514.html)
- BGC-Argo under Bay of Bengal and Arabian Sea cyclones — [UCAR record](https://impacts.ucar.edu/en/publications/quantifying-tropical-cyclones-effect-on-the-biogeochemical-proces/), [Remote Sensing 14(3) 529](https://www.mdpi.com/2072-4292/14/3/529), [Cyclone-enhanced chlorophyll in the Bay of Bengal (record)](https://katalog.lib.cas.cz/KNAV/EdsRecord/edsair,edsair.doi...........166bb04e7e4e0ce76f58bd5a89a1de36)
- North Indian Ocean SST response, 96 storms — [record](https://agris.fao.org/search/en/records/67597c22c7a957febdf8881c)

Gap-filling
- ExDINEOF for Cyclone Hudhud, Frontiers in Marine Science 2023 — [article](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1230116/full)
- DINEOF for Arabian Sea cyclone blooms, ISPRS 2020 — [PDF](https://isprs-archives.copernicus.org/articles/XLIII-B3-2020/1391/2020/isprs-archives-XLIII-B3-2020-1391-2020.pdf)
- DINCAE vs DINEOF for ocean colour, GMD 2026 — [article](https://gmd.copernicus.org/articles/19/1619/2026/), [code](https://zenodo.org/records/14905369)

Mechanism
- Typhoon Trami: surface chlorophyll gains mostly entrainment of the deep maximum, Biogeosciences 2021 — [article](https://bg.copernicus.org/articles/18/849/2021/)
- Global impact of tropical cyclones on primary production (about 1%) — [PDF](https://hal.sorbonne-universite.fr/hal-01331852/document)
- Mediterranean cyclones uplift the deep chlorophyll maximum, 2025 preprint — [EGUsphere](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-3260/)
- Marine heatwaves are shaping the vertical structure of phytoplankton in the global ocean, Communications Earth & Environment 2025 — [article](https://www.nature.com/articles/s43247-025-02718-y)
- Marine heatwaves inhibit the spring bloom in the North-Western Mediterranean, GRL 2024 — [PDF](https://iris.cnr.it/retrieve/5221ea40-f00f-4250-b567-b0f895ed764b/Li_et_al_GRL_2024.pdf)
- Review of ocean–atmosphere interactions during North Indian Ocean cyclones (Thane vs Phailin) — [arXiv 2012.04384](https://arxiv.org/pdf/2012.04384)
- Kuttippurath et al. 2021 (translation speed, eddies) — [article](https://www.nature.com/articles/s41612-021-00166-x)
- Storm-centred composites with background variability, Southampton PhD thesis 2020 — [PDF](https://eprints.soton.ac.uk/437859/1/Jaume_Catany_Rafael_PhD_Thesis_Jan_2020.pdf)

Statistics
- Causal forests, Wager & Athey 2018 — [arXiv 1510.04342](https://arxiv.org/pdf/1510.04342); aerosol–cloud application — [arXiv 2110.15084](https://www.arxiv.org/pdf/2110.15084v3)
