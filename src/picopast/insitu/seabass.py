"""SeaBASS: a reader for .sb files and a candidate search through NASA's CMR.

Listing SeaBASS files needs no login; downloading them does. CMR knows each
file's name, time and footprint but not which fields it holds, so the search
below can only shortlist files by name. A candidate becomes a confirmed source
only after its header is read (``read_sb``) and found to hold cell counts.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import pandas as pd
import requests

from picopast.config import Box
from picopast.ingest import cmr

#: File-name hints for flow cytometry / cell abundance. Deliberately broad.
NAME_HINT = re.compile(r"fcm|flow_?cyt|cyto|pico|abun|cell|phyto|_CLS", re.I)
#: CMR allows at most five leading-wildcard conditions per query, and each
#: ``*x*`` pattern counts twice, so only two patterns can be sent at once.
SERVER_PATTERNS = ["*fcm*", "*cyto*"]
#: Words that mark a SeaBASS experiment as regional to the study box. SeaBASS
#: collection footprints in CMR are all global, so the text is the only clue.
REGION_WORDS = re.compile(
    r"indian ocean|arabian|bengal|ekamsat|asiri|monsoon|\bindia\b|\boman\b|andaman|sri lanka|"
    r"maldives|lakshadweep",
    re.I,
)
_FIELD_PATTERNS = {
    "pro_cells_ml": r"prochlorococc|(^|_)pro(_|$)",
    "syn_cells_ml": r"synechococc|(^|_)syn(_|$)",
    "peuk_cells_ml": r"picoeuk|pico_euk|(^|_)peuk(_|$)",
}


def read_sb(source: str | Path | io.TextIOBase) -> tuple[dict[str, str], pd.DataFrame]:
    """Parse a SeaBASS file into (header dict, data frame with missing values as NaN)."""
    if isinstance(source, (str, Path)):
        text = Path(source).read_text(errors="replace")
    else:
        text = source.read()
    header: dict[str, str] = {}
    lines = text.splitlines()
    for i, line in enumerate(lines):
        line = line.strip()
        if line.lower() == "/end_header":
            body = lines[i + 1 :]
            break
        if line.startswith("/") and "=" in line:
            key, _, value = line[1:].partition("=")
            header[key.strip().lower()] = value.strip()
    else:
        raise ValueError("no /end_header line: not a SeaBASS file")
    fields = [f.strip() for f in header["fields"].split(",")]
    sep = {"comma": ",", "tab": "\t", "space": r"\s+"}[header.get("delimiter", "comma").lower()]
    df = pd.read_csv(
        io.StringIO("\n".join(body)), sep=sep, names=fields, engine="python", comment="!"
    )
    missing = header.get("missing")
    if missing is not None:
        df = df.replace(float(missing), float("nan"))
    return header, df


def cell_count_fields(header: dict[str, str]) -> dict[str, str]:
    """Map our group columns to the file's field names, where the file has them."""
    fields = [f.strip() for f in header.get("fields", "").split(",")]
    found = {}
    for column, pattern in _FIELD_PATTERNS.items():
        for field in fields:
            if re.search(pattern, field, re.I):
                found[column] = field
                break
    return found


def _in_box(spatial: str, box: Box) -> bool | None:
    """True/False if the CMR footprint meets the box; None if the granule has no footprint."""
    if not spatial:
        return None
    for part in spatial.split("; "):
        v = [float(x) for x in part.split()]
        if len(v) == 2:  # point: lat lon
            if box.contains(v[0], v[1]):
                return True
        elif len(v) == 4:  # box: south west north east
            if v[0] <= box.lat_max and v[2] >= box.lat_min and v[1] <= box.lon_max and v[3] >= box.lon_min:
                return True
    return False


def candidates(box: Box) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(collections worth opening, candidate files).

    Two passes, because CMR cannot search SeaBASS by region or by measured field:

    * **Regional experiments** (title or summary names the region): every file
      is listed; a file is a candidate if its name hints at cell counts and its
      footprint is inside the box or missing.
    * **All other experiments** (global programmes, other regions): only files
      named ``*fcm*`` or ``*cyto*`` are listed, and only those whose footprint
      is inside the box are kept.

    A file with no footprint in CMR is marked ``in_box = unknown``; its
    positions can only be checked after download.
    """
    col_rows, file_rows = [], []
    for c in cmr.collections(provider="OB_DAAC"):
        text = " ".join(str(c.get(k, "")) for k in ("short_name", "title", "summary"))
        regional = bool(REGION_WORDS.search(text))
        error = ""
        try:
            files = cmr.granules(c["short_name"], name_pattern=None if regional else SERVER_PATTERNS)
        except requests.RequestException as e:
            files, error = [], str(e)[:120]
        n_in = n_unknown = kept = 0
        for g in files:
            inside = _in_box(g.spatial, box)
            n_in += inside is True
            n_unknown += inside is None
            hint = bool(NAME_HINT.search(g.name))
            if (hint and inside is not False) if regional else (inside is True):
                kept += 1
                file_rows.append(
                    {
                        "collection": c["short_name"], "file": g.name, "time_start": g.start[:10],
                        "in_box": "yes" if inside else "unknown", "regional": regional,
                        "footprint": g.spatial, "url": g.url,
                    }
                )
        if regional or kept or error:
            col_rows.append(
                {
                    "collection": c["short_name"], "title": c.get("title", ""),
                    "regional": regional, "time_start": c.get("time_start", "")[:10],
                    "time_end": c.get("time_end", "")[:10], "files_listed": len(files),
                    "files_in_box": n_in, "files_no_footprint": n_unknown,
                    "candidate_files": kept, "error": error,
                }
            )
    file_cols = ["collection", "file", "time_start", "in_box", "regional", "footprint", "url"]
    return pd.DataFrame(col_rows), pd.DataFrame(file_rows, columns=file_cols)
