# MOANA tables (NASA OCSSW)

These two files hold every tunable constant in NASA's MOANA algorithm. Neither
is printed in Lange et al. (2020) or in the ATBD, so the algorithm cannot be
run without them.

| File | Content |
| --- | --- |
| `pca_picophyto.h5` | `component` (124 x 45 PCA loadings), `wavelength` (414-660 nm, 2 nm steps) |
| `picophyt.json` | `npc` = 17 and the zero-padded regression coefficients per group |

Source: `https://oceandata.sci.gsfc.nasa.gov/manifest/tags/T2023.31/common/<file>`
(NASA Ocean Color Science Software, bundle `common`), downloaded 2026-11-10.

```
478ace8798ce525bc6a09567b3773cadd532056c13d20ce6daf0a529f46deb57  pca_picophyto.h5
cd6b5c1cca879bc97acf6c5a3d4f47bf658d5a839b99ff4b1f534eac881756d4  picophyt.json
```

`picopast.moana.load_luts()` refuses to run if either checksum changes. The
same checksums are published by the EPFT-UP project
(<https://github.com/ocean-colour/EPFT-UP>), whose documentation of the slot
layout and of the ATBD-vs-operational coefficient mismatch we rely on.
