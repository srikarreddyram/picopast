"""Search NASA's Common Metadata Repository (CMR). No login needed for searching."""

from __future__ import annotations

from dataclasses import dataclass

import requests

CMR = "https://cmr.earthdata.nasa.gov/search"
_TIMEOUT = 90
_PAGE = 2000


@dataclass(frozen=True)
class Granule:
    name: str
    collection: str
    start: str  # ISO time, may be ""
    end: str
    size_mb: float
    url: str  # https data link, "" if none
    spatial: str  # CMR points/boxes as text, "" when the granule has no footprint


def _get(path: str, params: list[tuple[str, str]]) -> requests.Response:
    r = requests.get(f"{CMR}/{path}", params=params, timeout=_TIMEOUT)
    r.raise_for_status()
    return r


def collections(
    short_name: str | None = None,
    provider: str | None = None,
    bbox: str | None = None,
    keyword: str | None = None,
) -> list[dict]:
    """Collection records (short_name, version_id, time_start, time_end, title, id)."""
    base = [("page_size", str(_PAGE))]
    for key, value in (
        ("short_name", short_name),
        ("provider", provider),
        ("bounding_box", bbox),
        ("keyword", keyword),
    ):
        if value:
            base.append((key, value))
    out, page = [], 1
    while True:
        entries = _get("collections.json", base + [("page_num", str(page))]).json()["feed"]["entry"]
        out += entries
        if len(entries) < _PAGE:
            return out
        page += 1


def granules(
    short_name: str,
    start: str | None = None,
    end: str | None = None,
    name_pattern: str | list[str] | None = None,
    bbox: str | None = None,
) -> list[Granule]:
    """All granules of a collection matching a time range, file-name pattern and/or box.

    ``name_pattern`` uses ``*`` and ``?`` wildcards and is matched by CMR against
    the file name, which is how daily files are told apart from 8-day composites
    and from the climatologies that also overlap any date range.
    """
    base = [("short_name", short_name), ("page_size", str(_PAGE)), ("sort_key", "start_date")]
    if start or end:
        base.append(("temporal", f"{start or ''},{end or ''}"))
    if bbox:
        base.append(("bounding_box", bbox))
    if name_pattern:
        patterns = [name_pattern] if isinstance(name_pattern, str) else list(name_pattern)
        base += [("readable_granule_name[]", p) for p in patterns]  # several patterns are OR-ed
        base.append(("options[readable_granule_name][pattern]", "true"))
    out, page = [], 1
    while True:
        entries = _get("granules.json", base + [("page_num", str(page))]).json()["feed"]["entry"]
        for e in entries:
            links = [
                link["href"]
                for link in e.get("links", [])
                if link.get("href", "").startswith("https")
                and link.get("rel", "").endswith("/data#")
                and not link.get("inherited")
            ]
            out.append(
                Granule(
                    name=e.get("producer_granule_id") or e["title"],
                    collection=short_name,
                    start=e.get("time_start", ""),
                    end=e.get("time_end", ""),
                    size_mb=float(e.get("granule_size") or 0.0),
                    url=links[0] if links else "",
                    spatial="; ".join(e.get("points", []) + e.get("boxes", [])),
                )
            )
        if len(entries) < _PAGE:
            return out
        page += 1


def hits(short_name: str, start: str | None = None, end: str | None = None,
         name_pattern: str | None = None) -> int:
    """Number of matching granules, without listing them."""
    params = [("short_name", short_name), ("page_size", "0")]
    if start or end:
        params.append(("temporal", f"{start or ''},{end or ''}"))
    if name_pattern:
        params += [
            ("readable_granule_name", name_pattern),
            ("options[readable_granule_name][pattern]", "true"),
        ]
    return int(_get("granules.json", params).headers.get("CMR-Hits", 0))
