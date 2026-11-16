"""The analysis grid: the native L3M "4km" grid (1/24 degree, ~4.6 km) cut to the box.

We never interpolate ocean-colour products onto a new grid. NASA L3M 4km files
share one global 8640 x 4320 grid, so the cube is a plain index subset of it.
Coarser fields (OISST 0.25 deg, CMC 0.1 deg, GLORYS 1/12 deg) are brought onto
it by nearest-neighbour lookup, which keeps their native values visible.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from picopast.config import Box

EARTH_RADIUS_KM = 6371.0088


def global_axes(resolution_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Pixel-centre latitudes (north to south) and longitudes of the global L3M grid."""
    n_lat = int(round(180.0 / resolution_deg))
    n_lon = int(round(360.0 / resolution_deg))
    lat = 90.0 - (np.arange(n_lat) + 0.5) * resolution_deg
    lon = -180.0 + (np.arange(n_lon) + 0.5) * resolution_deg
    return lat, lon


def box_axes(box: Box, resolution_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Pixel centres of the global grid that fall inside ``box``."""
    lat, lon = global_axes(resolution_deg)
    return (
        lat[(lat >= box.lat_min) & (lat <= box.lat_max)],
        lon[(lon >= box.lon_min) & (lon <= box.lon_max)],
    )


def subset_box(ds: xr.Dataset | xr.DataArray, box: Box):
    """Cut a lat/lon dataset to the box, whatever the latitude order."""
    lat = ds["lat"].values
    lat_slice = (
        slice(box.lat_max, box.lat_min) if lat[0] > lat[-1] else slice(box.lat_min, box.lat_max)
    )
    return ds.sel(lat=lat_slice, lon=slice(box.lon_min, box.lon_max))


def nearest_to_grid(da: xr.DataArray, lat: np.ndarray, lon: np.ndarray) -> xr.DataArray:
    """Nearest-neighbour lookup of a coarser field at the analysis grid's pixel centres."""
    return da.sel(
        lat=xr.DataArray(lat, dims="lat"), lon=xr.DataArray(lon, dims="lon"), method="nearest"
    ).assign_coords(lat=lat, lon=lon)


def block_id(lat, lon, block_deg: float = 2.0) -> np.ndarray:
    """Integer id of the ``block_deg`` x ``block_deg`` block holding each point.

    Used as the group label for spatial-block cross-validation: every pixel in
    a block goes to the same fold, so neighbours never straddle train and test.
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    row = np.floor((lat + 90.0) / block_deg).astype(np.int64)
    col = np.floor((lon + 180.0) / block_deg).astype(np.int64)
    n_col = int(np.ceil(360.0 / block_deg))
    return row * n_col + col


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (array-friendly)."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dlmb = np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def coast_mask(ocean: np.ndarray, lat: np.ndarray, resolution_deg: float, mask_km: float = 25.0):
    """True where a pixel is ocean and at least ``mask_km`` from the nearest land pixel.

    ``ocean`` is a 2-D (lat, lon) boolean array on the analysis grid, True where
    the pixel is water. The distance is a Euclidean distance transform with the
    east-west pixel size scaled by cos(latitude); over a 25-degree tropical box
    that is accurate to well under one pixel.
    """
    from scipy.ndimage import distance_transform_edt

    ocean = np.asarray(ocean, dtype=bool)
    dy_km = np.radians(resolution_deg) * EARTH_RADIUS_KM
    dx_km = dy_km * float(np.cos(np.radians(np.mean(lat))))
    dist_km = distance_transform_edt(ocean, sampling=(dy_km, dx_km))
    return ocean & (dist_km >= mask_km)
