# Route B: rebuild the spectrum, then run MOANA

Generated 2026-11-17 by `picopast route-b`. An MLP (two layers of 64) learns the leading 12 principal components of the PACE spectrum (400–680 nm) from the five virtual-VIIRS bands; MOANA then runs on the rebuilt spectrum with SST. The MLP never sees real VIIRS. Scored against our MOANA on the PACE spectrum, with the same spatial blocks and temporal cut as Routes A and C; predictions clipped to the training range (share shown).

Rebuilt spectrum, median relative error on held-out pixels: 400 nm 0.8%, 425 nm 0.5%, 450 nm 0.5%, 475 nm 0.5%, 500 nm 0.5%, 525 nm 0.7%, 550 nm 0.6%, 575 nm 2.3%, 623 nm 5.0%, 643 nm 5.7%, 656 nm 10.3%, 668 nm 1.7%.

Two versions on real VIIRS: B trained on virtual bands only (never sees real VIIRS), and B trained on real VIIRS SNPP bands paired with the PACE spectrum on the same pixels and 8-day periods (the PRD's intent).

| group | split | B, virtual VIIRS | direct model, virtual VIIRS (RQ2) | B (virtual-trained), real VIIRS | B (real-trained), real VIIRS | Route A, real VIIRS | Route C, real VIIRS |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Prochlorococcus | spatial | 0.48 | 0.59 | -0.16 | 0.06 | 0.18 | -0.10 |
| Prochlorococcus | temporal | 0.51 | 0.57 | -0.17 | 0.11 | 0.18 | -0.04 |
| Synechococcus | spatial | 0.55 | 0.64 | -0.59 | 0.06 | 0.20 | -0.11 |
| Synechococcus | temporal | 0.56 | 0.61 | -0.91 | 0.10 | 0.11 | -0.01 |
| picoeukaryotes | spatial | 0.46 | 0.56 | 0.17 | 0.28 | 0.39 | -0.09 |
| picoeukaryotes | temporal | 0.51 | 0.56 | 0.24 | 0.32 | 0.39 | 0.01 |

Share of predictions clipped to the training range: real 1%, real-trained 0%, virtual 0%.

Skill = 1 − RMSE / RMSE of a seasonal climatology fitted in the same training part (> 0 beats it).
