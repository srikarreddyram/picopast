# Literature survey: has the new direction already been done?

22 November 2026. The question: before PicoPast expands beyond picophytoplankton,
which of the candidate directions are already published, and where is there room
for something new? Novelty is worded as "limited evidence in the literature we
reviewed", never "first ever".

**How this was done.** About 25 web searches (22 November 2026), aimed at
2022–2026 papers, preprints and conference abstracts. For the three papers closest
to the proposed gap, the abstracts were read from Crossref and NCBI. Full texts
behind paywalls were not read. Every reference below comes from those searches;
none is cited from memory. Limits and a ready-made Scopus query are at the end.

## Verdict at a glance

| Candidate direction | Status in the literature | Verdict |
| --- | --- | --- |
| Marine heatwaves make cyclones intensify faster | Global study of 1,600 landfalling storms (Science Advances, April 2026); global ML model (Earth's Future, 2024); North Indian Ocean case studies of Amphan, Tauktae and Biparjoy | **Taken** |
| Forecasting North Indian Ocean heatwaves with ML | U-Net forecasts up to 10 weeks ahead for the North Indian Ocean and Arabian Sea (Environmental Data Science, 2026); global deep-learning forecasts (arXiv, 2024) | **Taken** |
| What drives North Indian Ocean heatwaves | Arabian Sea 2010 case (Ocean Science, 2022); Bay of Bengal 1982–2021 (Scientific Reports, 2023); basin-wide modes (Frontiers in Climate, 2026); heat-budget work (EGU 2025) | **Taken** |
| Is Arabian Sea productivity rising or falling? | Long-running debate (Goes; Roxy 2016; Modi & Roxy 2023; others) | **Crowded**: one more voice |
| Cyclone-driven blooms in the North Indian Ocean | All cyclones 1997–2019 (npj Climate and Atmospheric Science, 2021); full review (JGR Oceans, 2025) | **Well studied** |
| Heatwaves lower chlorophyll in the North Indian Ocean | 1982–2021, pre-monsoon (Marine Pollution Bulletin, 2023) | **Studied** |
| Cyclones end heatwaves | Western Pacific (Science Advances, 2025); in this basin only the Amphan case | **Partly open** for this basin |
| ML prediction of cyclone cooling | Northwest Pacific (GRL 2023; IEEE TGRS 2024) | **Open** for this basin |
| ML prediction of post-cyclone *blooms* | Only a 2016 conference abstract (Taiwan Strait) | **Open** |
| **A cyclone crossing a heatwave: what happens to the plankton bloom?** | Each half is published (heatwaves suppress chlorophyll; cyclones cause blooms). No study found that tests the combination across many storms, in any basin. The closest is one Mediterranean storm (Medicane Daniel, 2026). | **Open: the gap** |
| PACE (hyperspectral) view of cyclones or plankton in this region | PACE after a hurricane: one Gulf of Mexico study (Milton, 2024). Nothing published with PACE in the North Indian Ocean. | **Open** |

## 1. Cyclones and marine heatwaves

**Heatwaves fuel cyclones (taken).**
- **Global studies.** Radfar et al. (2026, Science Advances) analysed about 1,600 landfalling storms since 1981. Storms that rapidly intensified over marine heatwaves caused about 60% more billion-dollar disasters, and they name the North Indian Ocean among the regions with rising trends. Radfar et al. (2024, Earth's Future) built an ensemble ML model that predicts heatwave-driven rapid intensification worldwide. Choi et al. (2024, Communications Earth & Environment) report much stronger storms over heatwaves.
- **North Indian Ocean case studies.** Amphan 2020 (Rathore et al. 2022, Frontiers in Climate): a heatwave with anomalies above 2.5 °C helped it reach super-cyclone strength in about 24 hours, and the storm then ended the heatwave. Tauktae 2021 (Theoretical and Applied Climatology, 2023). Biparjoy 2023 (Gupta & Deogharia 2025, Natural Hazards; a Climate Dynamics 2024 study on "spiciness"; INCOIS scientists via press reports).

**Cyclones end heatwaves (open for this basin).**
- **Pacific.** Pun, Lin & Wu (2025, Science Advances): Typhoon Bavi cooled the sea by up to 7 °C, ended a heatwave and held back new ones for up to eight months.
- **Subsurface warming.** A 2026 EMS abstract (Lian & Wang) reports that storms can also *create* subsurface heatwaves.
- **This basin.** Only the Amphan case, plus a remark that pre-monsoon storms often end Arabian Sea heatwaves (Chatterjee et al. 2022). A review of ocean–cyclone interactions in the basin explicitly lists "a statistical assessment of the recovery time of the SSTs and the ocean heat content" as still to be done (arXiv 2012.04384).

**Cyclone blooms in the North Indian Ocean (well studied).**
- **Basin-wide.** Kuttippurath et al. (2021) analysed all cyclones from 1997 to 2019. Blooms average about 1.65 mg m⁻³. They are larger under slower storms and when storms cross cold-core eddies. They are modulated by ENSO and the IOD. The abstract does not test heatwaves or pre-storm warm anomalies, and uses no predictive model.
- **Review.** Bhattacharya, Ghosal & Chakraborty (2025, JGR Oceans) review the whole biogeochemical response to cyclones. Their abstract does not mention heatwaves.
- **Stratification matters, in single cases.** Phailin (2013) produced a small bloom under strong stratification, while Thane (2011), a weaker storm, produced a large one.

**Heatwaves and chlorophyll (studied, without cyclones).**
- **This basin.** Krishnapriya et al. (2023, Marine Pollution Bulletin): pre-monsoon heatwaves lower chlorophyll and net primary production, by 10% in the Arabian Sea and 2% in the Bay of Bengal. Cyclones are not considered.
- **Global.** Le Grix et al. (2021, Biogeosciences) on compound high-temperature, low-chlorophyll extremes; a 2026 Remote Sensing paper on compound heatwave–bloom events; a 2025 study using BGC-Argo floats (Communications Earth & Environment) showing heatwaves reshape the vertical structure of phytoplankton.

**The combination (not found).**
- **Searches run.** Heatwave preconditioning followed by a cyclone and a bloom was searched for in the North Indian Ocean, the South China Sea and western Pacific, the Gulf of Mexico and Atlantic, and globally.
- **What turned up.** Single-storm studies (Idalia, Ida, Milton, Rey) and global bloom statistics (a 2016 study of about 1,000 cyclones; Da, Foltz & Balaguru 2021). None tests heatwave status as a control.
- **Nearest match.** A 2026 Biogeosciences study of Medicane Daniel, a single Mediterranean storm.

**ML for cyclone responses.**
- **Cooling.** Cui and colleagues predicted cold wakes in the northwest Pacific with a random forest (GRL 2023) and with XGBoost plus SHAP explanations (IEEE TGRS 2024). Their predictors include mixed-layer depth but not salinity.
- **Salinity.** Eusebi et al. (2025, Environmental Research Letters) use satellite salinity in ML for rapid intensification.
- **Blooms.** For post-cyclone blooms only a 2016 AGU abstract was found (Taiwan Strait).
- **A warning for this basin.** A Bay of Bengal neural ocean emulator got *worse* when given cyclone tracks, because the cyclones were "too rare to learn" (Islam, arXiv 2609.04635, 2026). With few storms, small and simple models are the honest choice.

## 2. Heatwave forecasting and drivers (taken)

- **Forecasting.** Howard et al. (2026, Environmental Data Science) give skilful subseasonal forecasts of Indian Ocean heatwaves up to 10 weeks ahead with a U-Net, focused on the North Indian Ocean and Arabian Sea. Ning et al. (2024, arXiv) report global forecasts up to six months ahead.
- **Drivers.** Chatterjee et al. (2022, Ocean Science; reduced latent heat loss and a shallow mixed layer); Bay of Bengal 1982–2021 (Scientific Reports, 2023; more than 70% of heatwave days follow El Niño); Joseph et al. (2026, Frontiers in Climate; ENSO and intraseasonal modes); a heat-budget study (EGU 2025).

## 3. Arabian Sea productivity (crowded)

- **Falling.** Roxy et al. (2016, GRL) report declining productivity under Indian Ocean warming. Modi & Roxy (2023, book chapter) expect the decline to continue.
- **Rising, or not clear.** Goes and colleagues argue monsoon upwelling has strengthened (NASA Ocean Color Research Team presentations, 2007 and 2012). Other groups question both.
- **Records.** Satellite records are still short compared with monsoon variability.

A new long record would add one more voice to a live debate rather than open new ground.

## 4. PACE in the North Indian Ocean (open)

- **What's published.** The only PACE storm study found is Hurricane Milton (2024), which tracked a Karenia brevis red tide and shifts in the phytoplankton community. No PACE study of plankton communities in this region was found.
- **What exists for this region.** NASA's PVST project collected Arabian Sea validation data for PACE in 2024; PicoPast's check E already uses those HPLC samples. A 2026 Remote Sensing paper maps Arabian Sea phytoplankton size classes, but from VIIRS, not PACE.
- **Aerosols.** PACE aerosol products (HARP2, SPEXone, OCI) were released in 2025. No Arabian Sea study was found.

## Proposed problem statement

> **When cyclones hit heatwaves: how marine heatwaves change the physical and
> biological response of the North Indian Ocean to tropical cyclones, 2012–2026,
> from NASA satellites.**

| | Question | Why it is open |
| --- | --- | --- |
| RQ1 | Do cyclones that cross heatwave water cool the surface differently, how often do they end the heatwave, and how long until the ocean re-warms? | The basin-wide statistics are listed as missing by a regional review; the systematic study exists only for the Pacific. |
| RQ2 (core) | Is the post-cyclone bloom weaker, stronger or later when the storm crosses a heatwave, after allowing for storm strength, speed, eddies (sea level) and the salinity barrier layer (SMAP)? | No study found in any basin. |
| RQ3 | Which phytoplankton groups respond after storms (PACE, 2024–2026; picophytoplankton included as one part)? | Only one PACE storm study, in the Gulf of Mexico. |
| RQ4 | Can a small, explainable model predict each storm's cooling and bloom from conditions before it arrives, and does heatwave status add skill? | ML for blooms was not found; for cooling, only in the Pacific. |

**Competing hypotheses for RQ2, both physically reasonable:**
- **Weaker bloom:** heatwave water is stratified and nutrient-poor, so mixing has further to reach.
- **Stronger bloom:** heatwave mixed layers are shallow (Chatterjee et al. 2022), so mixing breaks them more easily.

Either result is publishable, and the test is decided before the data are compared.

**Data (all NASA except the storm tracks):**
- VIIRS chlorophyll, already processed for 2012–2026
- PACE OCI
- SMAP salinity
- NASA-SSH sea level
- CCMP winds
- GPM IMERG rainfall
- MERRA-2 heat fluxes
- OISST (NOAA), already processed; MUR as the NASA alternative
- Storm tracks from IBTrACS (NOAA, free)

**Feasibility (quick count, 22 November; definitions to be fixed in a pre-registration):**
- **Storms:** 81 storms reached cyclonic-storm strength (≥ 34 kt) inside 0–25°N, 50–100°E between 2012 and November 2026: 28 in the Arabian Sea and 53 in the Bay of Bengal, 4 to 9 a year.
- **Heatwave or not:** In the 8-day period before arrival, 42 crossed a 0.25° cell with heatwave conditions on at least half its days, and 32 crossed only normal water (below 10%). This criterion is generous, and the final definition will lower the first number.
- **PACE era:** 13 storms since March 2024: Remal, Asna, Dana, Fengal, Shakhti, Montha, Senyar, Ditwah, Arnab, Murasu, Akvan and two unnamed.

**Risks:**
- Clouds hide the ocean under storms. We use post-storm windows of 8-day composites, the standard approach.
- Heatwave storms cluster in certain seasons. We compare like with like by season and basin.
- With 81 storms, ML must stay small and be tested on storms it never saw.
- Others are moving fast in this area (two of the closest papers are from 2026). See below.

## How we protect priority

A dated, public pre-registration in the repository (as for Route D and check E)
records the questions and methods before any result exists. Posting a short
methods preprint (for example on ESS Open Archive) early would make the claim
citable. Re-run this search before submission.

## Limits of this survey and a query for the library

- **Search engines.** The search used general web search, not Scopus or Web of Science. Paywalled full texts were not read; abstracts were read for the three nearest papers.
- **Not covered.** Theses, Indian conference proceedings and non-English papers. Very recent preprints may be missed.
- **Before final sign-off,** the library or guide could run in Scopus:
  `TITLE-ABS-KEY(("marine heatwave*" OR "warm anomal*") AND (cyclone* OR typhoon* OR hurricane*) AND (chlorophyll OR phytoplankton OR bloom* OR "primary production"))`

## References (as found, with links)

Cyclones and heatwaves
- Radfar et al. 2026, Science Advances — [Synergistic impact of marine heat waves and rapid intensification exacerbates tropical cyclone destructive power worldwide](https://www.science.org/doi/10.1126/sciadv.adu1733)
- Radfar et al. 2024, Earth's Future — [Global Predictability of Marine Heatwave Induced Rapid Intensification of Tropical Cyclones](https://par.nsf.gov/biblio/10640505)
- Choi et al. 2024, Communications Earth & Environment — [article](https://www.nature.com/articles/s43247-024-01239-4)
- Rathore et al. 2022, Frontiers in Climate — [Interactions Between a Marine Heatwave and Tropical Cyclone Amphan in the Bay of Bengal in 2020](https://www.frontiersin.org/journals/climate/articles/10.3389/fclim.2022.861477/full)
- Tauktae ocean conditions, Theoretical and Applied Climatology 2023 — [article](https://link.springer.com/article/10.1007/s00704-023-04480-7)
- Gupta & Deogharia 2025, Natural Hazards (Biparjoy) — [record](https://ideas.repec.org/a/spr/nathaz/v121y2025i15d10.1007_s11069-025-07513-8.html)
- Biparjoy and ocean spiciness, Climate Dynamics 2024 — [article](https://link.springer.com/content/pdf/10.1007/s00382-024-07109-1.pdf)
- INCOIS on heatwaves and Arabian Sea cyclones (press) — [PreventionWeb](https://www.preventionweb.net/news/scientists-link-marine-heatwaves-rising-cyclone-risks)
- Pun, Lin & Wu 2025, Science Advances — Suppression of marine heatwave activity by tropical cyclone-induced upper ocean cooling ([group page](https://poseidon.ihs.ncu.edu.tw/marine-heatwave.html))
- Lian & Wang 2026, EMS abstract — [cyclones and subsurface heatwaves](https://meetingorganizer.copernicus.org/EMS2026/EMS2026-61.html)
- Review of ocean–atmosphere interactions during North Indian Ocean cyclones — [arXiv 2012.04384](https://arxiv.org/pdf/2012.04384)
- Bay of Bengal cyclone cold wakes, 1982–2020 — [Environmental Science: Atmospheres 2022](https://pubs.rsc.org/en/content/articlehtml/2022/ea/d1ea00066g)

Cyclone blooms and heatwave biology
- Kuttippurath et al. 2021, npj Climate and Atmospheric Science — [Tropical storms trigger phytoplankton blooms in the deserts of north Indian Ocean](https://www.nature.com/articles/s41612-021-00166-x)
- Bhattacharya, Ghosal & Chakraborty 2025, JGR Oceans — [Biogeochemical Response to Tropical Cyclones in the Northern Indian Ocean: Understanding and Gray Areas](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JC022540)
- Krishnapriya et al. 2023, Marine Pollution Bulletin — [Marine heatwaves during the pre-monsoon season and their impact on Chlorophyll-a in the north Indian Ocean in 1982–2021](https://pubmed.ncbi.nlm.nih.gov/37988881)
- Le Grix et al. 2021, Biogeosciences — [Compound high-temperature and low-chlorophyll extremes](https://bg.copernicus.org/articles/18/2119/2021/)
- Compound heatwave–bloom events, Remote Sensing 2026 — [PDF](https://mdpi-res.com/d_attachment/remotesensing/remotesensing-18-01322/article_deploy/remotesensing-18-01322.pdf)
- Heatwaves and the vertical structure of phytoplankton, Communications Earth & Environment 2025 — [article](https://www.nature.com/articles/s43247-025-02718-y)
- Medicane Daniel, Biogeosciences 2026 — [article](https://bg.copernicus.org/articles/23/4271/2026/)
- Da, Foltz & Balaguru 2021, GRL — [Observed Global Increases in Tropical Cyclone Induced Ocean Cooling and Primary Production](https://ui.adsabs.harvard.edu/abs/2021GeoRL..4892574D)
- Global impact of tropical cyclones on primary production (2016, ~1,000 cyclones) — [PDF](https://hal.sorbonne-universite.fr/hal-01331852/document)
- Hurricane Idalia biogeochemistry (Saildrone, BGC-Argo) — [summary](https://www.saildrone.com/eu-en/scientific-papers/subsurface-biogeochemical-response-to-hurricane-idalia-within-a-cyclonic-eddy-and-river-plume-stratified-environment)

ML
- Cui and colleagues 2023, GRL — [random-forest cold-wake prediction](https://doaj.org/article/66b4675e33574aea912ab4516d5656e2); 2024 IEEE TGRS XGBoost follow-up — [record](https://openresearch.surrey.ac.uk/esploro/outputs/journalArticle/Predicting-Tropical-Cyclone-Induced-Sea-Surface-Temperature/99929042702346)
- Eusebi et al. 2025, Environmental Research Letters — [salinity in ML for rapid intensification](https://clima.caltech.edu/wp-content/uploads/2025/10/Eusebi_2025_Environ._Res._Lett._20_034010.pdf)
- Islam 2026, arXiv — [Too Rare to Learn: Prescribed Cyclone Tracks Degrade a Bay of Bengal Ocean Emulator](https://arxiv.org/pdf/2609.04635)
- Typhoon chlorophyll forecasting, AGU Ocean Sciences 2016 abstract — [abstract](https://agu.confex.com/agu/os16/preliminaryview.cgi/Paper86825.html)

Heatwave forecasting and drivers
- Howard et al. 2026, Environmental Data Science — [Skillful subseasonal Indian Ocean marine heatwave forecasts using a neural network](https://www.cambridge.org/core/journals/environmental-data-science/article/skillful-subseasonal-indian-ocean-marine-heatwave-forecasts-using-a-neural-network/6DBF5C0FAC53BF861EFF50EB8F741475)
- Ning et al. 2024, arXiv — [Advancing Marine Heatwave Forecasts: An Integrated Deep Learning Approach](https://arxiv.org/abs/2412.04475)
- Chatterjee et al. 2022, Ocean Science — [Arabian Sea 2010 heatwave](https://os.copernicus.org/articles/18/639/2022/)
- Bay of Bengal heatwaves 1982–2021, Scientific Reports 2023 — [article](https://link.springer.com/article/10.1038/s41598-023-39884-y)
- Joseph et al. 2026, Frontiers in Climate (preprint) — [North Indian Ocean heatwave modes](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-3886/egusphere-2025-3886.pdf)
- Heat-budget drivers, EGU 2025 — [abstract](https://meetingorganizer.copernicus.org/EGU25/EGU25-15053.html)

Arabian Sea productivity
- Roxy et al. 2016, GRL — [A reduction in marine primary productivity driven by rapid warming over the tropical Indian Ocean](https://typeset.io/pdf/a-reduction-in-marine-primary-productivity-driven-by-rapid-2hyr78yr71.pdf)
- Modi & Roxy 2023 — [book chapter](https://www.rocksea.org/bin/research/modi_roxy_indian_ocean_primary_productivity_springer_2023.pdf)
- Goes, NASA Ocean Color Research Team 2012 — [presentation](https://oceancolor.gsfc.nasa.gov/files/meetings/ocrt/201204/goes_as_presentations/9924-goes-aslo-arabian_sea/GOES-NASA-OCRT-2012-26-April2012pdf.pdf)

PACE and salinity
- Hurricane Milton and PACE (USGS 2026) — [Rapid intensification of Karenia brevis red tides after Hurricane Milton revealed by PACE](https://www.usgs.gov/publications/rapid-intensification-karenia-brevis-red-tides-after-hurricane-milton-revealed-pace)
- NASA PVST Northern Indian Ocean (PACE validation) — [SeaBASS](https://seabass.gsfc.nasa.gov/experiment/PVST_NORTHERN_INDIAN_OCEAN)
- Arabian Sea size classes from VIIRS, Remote Sensing 2026 — [PDF](https://mdpi-res.com/d_attachment/remotesensing/remotesensing-18-01451/article_deploy/remotesensing-18-01451.pdf)
- PACE aerosol products released (NASA Earthdata) — [notice](https://earthdata.nasa.gov/data/alerts-outages/pace-harp2-spexone-oci-products-released)
- SMAP and SMOS salinity under Bay of Bengal cyclones, Remote Sensing 2022 — [article](https://mdpi.com/2072-4292/14/15/3733)
