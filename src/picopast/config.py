"""Load the single study config (``config/picopast.yaml``).

Everything that defines the study (box, regions, seasons, products, periods)
is read from that file so a run is reproducible from one place.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "picopast.yaml"


@dataclass(frozen=True)
class Box:
    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float

    def contains(self, lat, lon):
        """Boolean (array-friendly) test for points inside the box."""
        return (
            (lat >= self.lat_min)
            & (lat <= self.lat_max)
            & (lon >= self.lon_min)
            & (lon <= self.lon_max)
        )

    @property
    def cmr_bbox(self) -> str:
        """west,south,east,north as CMR and ERDDAP expect it."""
        return f"{self.lon_min},{self.lat_min},{self.lon_max},{self.lat_max}"


class Config(dict):
    """A dict with a few typed helpers."""

    @property
    def box(self) -> Box:
        return Box(**self["study_box"])

    def region(self, name: str) -> Box:
        r = self["regions"][name]
        return Box(r["lat_min"], r["lat_max"], r["lon_min"], r["lon_max"])

    def regions(self, phase: int | None = None) -> dict[str, Box]:
        return {
            name: self.region(name)
            for name, r in self["regions"].items()
            if phase is None or r["phase"] <= phase
        }

    def path(self, key: str) -> Path:
        p = Path(self["paths"][key])
        return p if p.is_absolute() else PROJECT_ROOT / p

    def season_of_month(self, month: int) -> str:
        for season, months in self["seasons"].items():
            if month in months:
                return season
        raise ValueError(f"month {month} is not assigned to a season")


@lru_cache(maxsize=4)
def _load(path: str) -> Config:
    with open(path) as f:
        raw: dict[str, Any] = yaml.safe_load(f)
    months = sorted(m for ms in raw["seasons"].values() for m in ms)
    if months != list(range(1, 13)):
        raise ValueError("seasons must cover each month exactly once")
    return Config(raw)


def load_config(path: str | Path | None = None) -> Config:
    return _load(str(path or DEFAULT_CONFIG))
