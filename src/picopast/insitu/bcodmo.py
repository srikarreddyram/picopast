"""BCO-DMO via its public ERDDAP server (no login)."""

from __future__ import annotations

import io

import pandas as pd
import requests

from picopast.config import Box

ERDDAP = "https://erddap.bco-dmo.org/erddap"
#: Global Pro/Syn/picoeukaryote compilation (Flombaum et al.), the main BCO-DMO source.
GLOBAL_COMPILATION = "bcodmo_dataset_793451"
SEARCH_TERMS = ["prochlorococcus", "synechococcus", "picoeukaryotes", "picoplankton", "flow cytometry"]


def search(box: Box, terms: list[str] = SEARCH_TERMS) -> pd.DataFrame:
    """Datasets whose footprint meets the box and whose metadata mention any term."""
    found: dict[str, str] = {}
    for term in terms:
        r = requests.get(
            f"{ERDDAP}/search/advanced.json",
            params={
                "searchFor": term, "page": 1, "itemsPerPage": 1000,
                "minLat": box.lat_min, "maxLat": box.lat_max,
                "minLon": box.lon_min, "maxLon": box.lon_max,
            },
            timeout=120,
        )
        if r.status_code == 404:  # ERDDAP answers 404 for "no matches"
            continue
        r.raise_for_status()
        table = r.json()["table"]
        i_id, i_title = table["columnNames"].index("Dataset ID"), table["columnNames"].index("Title")
        for row in table["rows"]:
            found[row[i_id]] = row[i_title]
    return pd.DataFrame(sorted(found.items()), columns=["dataset_id", "title"])


def fetch_global_compilation(box: Box) -> pd.DataFrame:
    """Rows of the global compilation inside the box, in the standard layout."""
    cols = "Year,Day,latitude,longitude,depth,Prochlorococcus,Synechococcus,Pico_eukaryotes"
    query = (
        f"{cols}&latitude>={box.lat_min}&latitude<={box.lat_max}"
        f"&longitude>={box.lon_min}&longitude<={box.lon_max}"
    )
    url = f"{ERDDAP}/tabledap/{GLOBAL_COMPILATION}.csv?{query}"
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    raw = pd.read_csv(io.StringIO(r.text), skiprows=[1])  # row 1 holds units
    time = pd.to_datetime(raw["Year"].astype(int).astype(str), format="%Y") + pd.to_timedelta(
        raw["Day"] - 1, unit="D"
    )
    return pd.DataFrame(
        {
            "source": "bcodmo",
            "dataset_id": GLOBAL_COMPILATION,
            # The compilation carries no cruise label. Year-month is a stand-in
            # until records are traced back to their original cruises.
            "cruise": "bcodmo_" + time.dt.strftime("%Y-%m"),
            "station": None,
            "time_utc": time,
            "lat": raw["latitude"],
            "lon": raw["longitude"],
            "depth_m": raw["depth"],
            "pro_cells_ml": raw["Prochlorococcus"],
            "syn_cells_ml": raw["Synechococcus"],
            "peuk_cells_ml": raw["Pico_eukaryotes"],
            "method": "compilation (mostly flow cytometry)",
            "url": f"{ERDDAP}/tabledap/{GLOBAL_COMPILATION}.html",
        }
    )
