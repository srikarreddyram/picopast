"""PANGAEA: search by box and read each dataset's public metadata (no login)."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

from picopast.config import Box

SEARCH = "https://www.pangaea.de/advanced/search.php"
QUERIES = ["Prochlorococcus", "Synechococcus", "picoeukaryotes", "picophytoplankton flow cytometry"]
_GROUP_PATTERNS = {
    "has_pro": r"prochlorococc",
    "has_syn": r"synechococc",
    "has_peuk": r"picoeuk|pico-euk|ultra eukaryot|eukaryot.*pico|picophytoeuk",
}


def search(box: Box, queries: list[str] = QUERIES, page: int = 500) -> dict[str, float]:
    """DOI -> best relevance score for datasets inside the box matching any query."""
    found: dict[str, float] = {}
    for q in queries:
        offset = 0
        while True:
            r = requests.get(
                SEARCH,
                params={
                    "q": q, "count": page, "offset": offset,
                    "minlat": box.lat_min, "maxlat": box.lat_max,
                    "minlon": box.lon_min, "maxlon": box.lon_max,
                },
                timeout=120,
            )
            r.raise_for_status()
            body = r.json()
            for hit in body["results"]:
                doi = hit["URI"].removeprefix("doi:")
                found[doi] = max(found.get(doi, 0.0), float(hit.get("score", 0.0)))
            offset += page
            if offset >= body["totalCount"] or not body["results"]:
                break
    return found


def metadata(doi: str, cache_dir: Path, pause_s: float = 0.25) -> dict:
    """JSON-LD metadata of one dataset, cached on disk."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / (doi.replace("/", "_") + ".json")
    if cached.exists():
        return json.loads(cached.read_text())
    r = requests.get(f"https://doi.pangaea.de/{doi}", params={"format": "metadata_jsonld"},
                     timeout=120)
    r.raise_for_status()
    cached.write_text(r.text)
    time.sleep(pause_s)  # be polite to a free service
    return r.json()


def _coverage(meta: dict) -> tuple[str, str]:
    t = meta.get("temporalCoverage") or ""
    start, _, end = t.partition("/")
    return start, end or start


def _centre(meta: dict) -> tuple[float, float]:
    geo = (meta.get("spatialCoverage") or {}).get("geo") or {}
    if "latitude" in geo:
        return float(geo["latitude"]), float(geo["longitude"])
    if "box" in geo:  # "south west north east"
        s, w, n, e = (float(x) for x in geo["box"].split())
        return (s + n) / 2, (w + e) / 2
    return float("nan"), float("nan")


def inventory(box: Box, cache_dir: Path) -> pd.DataFrame:
    """One row per PANGAEA dataset in the box that reports at least one of the three groups."""
    rows = []
    for doi, score in sorted(search(box).items()):
        try:
            meta = metadata(doi, cache_dir)
        except (requests.RequestException, json.JSONDecodeError):
            rows.append({"doi": doi, "error": "metadata unavailable"})
            continue
        names = " | ".join(
            str(v.get("name", "")) for v in meta.get("variableMeasured", []) if isinstance(v, dict)
        )
        start, end = _coverage(meta)
        lat, lon = _centre(meta)
        row = {
            "doi": doi, "title": meta.get("name", ""), "time_start": start, "time_end": end,
            "lat": lat, "lon": lon, "size": str(meta.get("size", {}).get("value", "")),
            "score": score, "url": f"https://doi.pangaea.de/{doi}",
        }
        for key, pattern in _GROUP_PATTERNS.items():
            row[key] = bool(re.search(pattern, names, re.I))
        rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    keep = df[["has_pro", "has_syn", "has_peuk"]].fillna(False).any(axis=1)
    return df[keep].sort_values("time_start").reset_index(drop=True)


#: Count columns in PANGAEA station tables -> our group columns.
_COUNT_COLUMNS = {
    "pro_cells_ml": re.compile(r"^Prochlorococcus \[#/ml\]"),
    "syn_cells_ml": re.compile(r"^Synechococcus \[#/ml\]"),
    "peuk_cells_ml": re.compile(r"^(Phytopl(ankton)? ultra euk|Picoeukaryotes?|Phytopl pico euk)\w* \[#/ml\]", re.I),
}


def _header_value(header: str, key: str) -> str | None:
    m = re.search(rf"{key}:\s*([^*\n]+)", header)
    return m.group(1).strip() if m else None


def read_station_table(text: str, doi: str) -> pd.DataFrame:
    """One PANGAEA station dataset (``?format=textfile``) in the common ship-count layout.

    Single-station datasets give position, time and campaign in the header's
    Event line and the counts in a tab-separated table below ``*/``.
    """
    from picopast.insitu.schema import standardise

    header, _, table = text.partition("*/\n")
    event = next((ln for ln in header.splitlines() if ln.startswith("Event(s):")), "")
    df = pd.read_csv(pd.io.common.StringIO(table), sep="\t")
    out = pd.DataFrame({
        "source": "pangaea",
        "dataset_id": doi,
        "cruise": _header_value(event, "CAMPAIGN"),
        "station": event.split("\t", 1)[-1].split(" ", 1)[0] if event else None,
        "time_utc": _header_value(event, "DATE/TIME"),
        "lat": float(_header_value(event, "LATITUDE") or "nan"),
        "lon": float(_header_value(event, "LONGITUDE") or "nan"),
        "depth_m": df["Depth water [m]"] if "Depth water [m]" in df else float("nan"),
        "method": "flow cytometry",
        "url": f"https://doi.pangaea.de/{doi}",
    }, index=df.index)
    for column, pattern in _COUNT_COLUMNS.items():
        match = next((c for c in df.columns if pattern.search(c)), None)
        out[column] = df[match] if match else float("nan")
    return standardise(out)


def fetch_station_counts(dois: list[str], cache_dir: Path, pause_s: float = 0.25) -> pd.DataFrame:
    """Counts from single-station PANGAEA datasets, text files cached on disk."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for doi in dois:
        cached = cache_dir / (doi.replace("/", "_") + ".txt")
        if not cached.exists():
            r = requests.get(f"https://doi.pangaea.de/{doi}", params={"format": "textfile"}, timeout=120)
            r.raise_for_status()
            cached.write_text(r.text)
            time.sleep(pause_s)  # be polite to a free service
        text = cached.read_text()
        if "Event(s):" not in text or "\t" not in text.split("*/\n", 1)[-1][:500]:
            continue  # not a single-station table (e.g. a global compilation)
        frames.append(read_station_table(text, doi))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
