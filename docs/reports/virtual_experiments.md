# Virtual experiments: each storm run with and without its heatwave

Generated 2026-11-28 by `picopast virtual`. Rules: the "Virtual experiments" section of `docs/cyclone_heatwave_preregistration.md`, committed before any real storm was simulated. Model: 1-D PWP mixing with NPZD biology (`picopast.column_model`). Differences are (with heatwave − without heatwave); intervals resample storms.

**Heatwave points (169 storms), heatwave removed. Surface bloom: in the model, a heatwave weakens it; depth-integrated chlorophyll: in the model, a heatwave weakens it.**

| output | heatwave points: effect [95%] | normal points, heatwave added: effect [95%] |
| --- | --- | --- |
| bloom (log10 surface chl) | -0.031 [-0.042, -0.021] | -0.030 [-0.037, -0.023] |
| depth-integrated chl change (mg m⁻²) | -0.997 [-1.865, -0.142] | -1.407 [-2.179, -0.647] |
| cooling (°C) | -0.265 [-0.318, -0.217] | -0.266 [-0.297, -0.237] |
| deepest mixed layer (m) | -4.130 [-4.805, -3.479] | -4.032 [-4.549, -3.537] |
| surface nitrate rise (mmol m⁻³) | -0.063 [-0.096, -0.034] | -0.104 [-0.133, -0.077] |

Validation (actual runs vs observations, Spearman): cooling ρ +0.51 (n 609), bloom ρ +0.29 (n 567).

| sensitivity | bloom effect, heatwave points [95%] |
| --- | --- |
| base | -0.031 [-0.042, -0.021] |
| cap20 | -0.014 [-0.026, -0.003] |
| par60 | -0.029 [-0.041, -0.018] |
| cd20 | -0.034 [-0.044, -0.024] |
| x05 | -0.035 [-0.045, -0.025] |

Runs with a non-finite output: base 0 of 1220; cap20 0 of 1220; par60 0 of 1220; cd20 0 of 1220; x05 0 of 1220.

Limits: one column, no advection or eddies; the heatwave enters as a warmer surface layer; fixed literature biology.
