"""Route D inputs: OISST v2 1° monthly (NOAA PSL), WOA23 monthly surface nitrate (NOAA NCEI) and MODIS-Aqua
monthly climatologies of chlorophyll and PAR (NASA OB.DAAC, Earthdata login from ~/.netrc).

Files land in their own directories under data/raw and are read only as data.
"""

from __future__ import annotations

import re
from pathlib import Path

import requests

from picopast.config import Config

UA = {"User-Agent": "picopast research client"}
OISST = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2/{name}"
WOA = "https://www.ncei.noaa.gov/thredds-ocean/fileServer/woa23/DATA/nitrate/netcdf/all/1.00/woa23_all_n{m:02d}_01.nc"
MODIS_MC = re.compile(r"AQUA_MODIS\.(\d{4})(\d{2})01_(\d{8})\.L3m\.MC\.")


def _get(url: str, out: Path) -> None:
    if out.exists():
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    with requests.get(url, stream=True, timeout=300, headers=UA) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    tmp.rename(out)


def mission_climatologies(names: list[str]) -> dict[int, str]:
    """The mission-long monthly climatology per calendar month (starting July 2002 - June 2003), latest end."""
    best: dict[int, tuple[str, str]] = {}
    for name in names:
        k = MODIS_MC.search(name)
        if not k:
            continue
        year, month = k.group(1), int(k.group(2))
        if not ((year == "2002" and month >= 7) or (year == "2003" and month < 7)):
            continue
        if month not in best or k.group(3) > best[month][0]:
            best[month] = (k.group(3), name)
    return {m: v[1] for m, v in best.items()}


def fetch(cfg: Config) -> None:
    raw = cfg.path("raw")
    for name in ("sst.mnmean.nc", "lsmask.nc"):
        _get(OISST.format(name=name), raw / "oisst_v2_1deg" / name)
    for m in range(1, 13):
        _get(WOA.format(m=m), raw / "woa23_nitrate" / f"woa23_all_n{m:02d}_01.nc")
    import earthaccess

    earthaccess.login(strategy="netrc")
    for short in ("MODISA_L3m_CHL", "MODISA_L3m_PAR"):
        grans = earthaccess.search_data(short_name=short, granule_name="*.L3m.MC.*9km*", count=500)
        by_name = {g["umm"]["GranuleUR"]: g for g in grans}
        keep = mission_climatologies(list(by_name))
        todo = [by_name[n] for n in keep.values()
                if not list((raw / "modis_clim").glob(f"*{MODIS_MC.search(n).group(0)}*"))]
        if todo:
            earthaccess.download(todo, str(raw / "modis_clim"))
