"""
Récupère les prévisions vent 10m depuis Open-Meteo pour différents modèles
et expose des fonctions identiques à vent_grib_deg / vent_grib_nm de vents.py.

Usage :
    from inputs.vents_openmeteo import get_V_deg, get_V_nm, get_meta
    meta  = get_meta('ecmwf_ifs025')
    V     = get_V_deg('ecmwf_ifs025')
    V_nm  = get_V_nm('ecmwf_ifs025')
"""

import time
import numpy as np
import requests
from datetime import datetime, timezone
from scipy.interpolate import RegularGridInterpolator
from shapely import contains_xy
from inputs.vents import land_geom

# ── Modèles disponibles ──────────────────────────────────────────────────────

MODELS = {
    "ecmwf_ifs025": {"label": "ECMWF IFS 0.25°", "days": 10},
    "gfs_seamless":  {"label": "NOAA GFS 0.25°",  "days": 16},
}

# ── Paramètres communs ────────────────────────────────────────────────────────

URL        = "https://api.open-meteo.com/v1/forecast"
BBOX       = (-80.0, 20.0, 25.0, 75.0)   # lon0, lat0, lon1, lat1
RESOLUTION = 5.0
MAX_LOC    = 10
CACHE_AGE  = 6                             # heures
HEADERS    = {"User-Agent": "Mozilla/5.0 (compatible; sailing-router/1.0)"}

# ── Cache interne (clé = model_id) ───────────────────────────────────────────

_cache: dict = {}


# ── Téléchargement ────────────────────────────────────────────────────────────

def _fetch_batch(lats_b, lons_b, model_id: str, days: int, max_retries=6):
    params = {
        "latitude":        ",".join(f"{la:.4f}" for la in lats_b),
        "longitude":       ",".join(f"{lo:.4f}" for lo in lons_b),
        "hourly":          "wind_speed_10m,wind_direction_10m",
        "models":          model_id,
        "forecast_days":   days,
        "wind_speed_unit": "ms",
        "timezone":        "UTC",
    }
    for attempt in range(max_retries):
        if attempt > 0:
            time.sleep(2 ** attempt)
        r = requests.get(URL, params=params, headers=HEADERS, timeout=60)
        if r.status_code == 429 and attempt < max_retries - 1:
            continue
        r.raise_for_status()
        data = r.json()
        return data if isinstance(data, list) else [data]


def _safe(lst):
    return np.array([x if x is not None else 0.0 for x in lst], dtype=float)


def _download(model_id: str):
    """Télécharge la grille du modèle et retourne (meta, interp_u, interp_v)."""
    cfg = MODELS[model_id]
    days = cfg["days"]

    lon0, lat0, lon1, lat1 = BBOX
    lons = np.arange(np.ceil(lon0 / RESOLUTION) * RESOLUTION,
                     lon1 + RESOLUTION / 2, RESOLUTION)
    lats = np.arange(np.ceil(lat0 / RESOLUTION) * RESOLUTION,
                     lat1 + RESOLUTION / 2, RESOLUTION)

    LON, LAT = np.meshgrid(lons, lats)
    flat_lons = LON.ravel()
    flat_lats = LAT.ravel()
    n_pts = len(flat_lons)
    n_lon, n_lat = len(lons), len(lats)

    times_h = None
    times_iso = None
    u_flat = v_flat = None

    for start in range(0, n_pts, MAX_LOC):
        time.sleep(0.6)
        locs = _fetch_batch(flat_lats[start:start + MAX_LOC],
                            flat_lons[start:start + MAX_LOC],
                            model_id, days)

        if times_h is None:
            times_iso = locs[0]["hourly"]["time"]
            ref = datetime.fromisoformat(times_iso[0])
            times_h = np.array(
                [(datetime.fromisoformat(t) - ref).total_seconds() / 3600
                 for t in times_iso], dtype=float
            )
            u_flat = np.zeros((n_pts, len(times_h)))
            v_flat = np.zeros((n_pts, len(times_h)))

        for k, loc in enumerate(locs):
            spd  = _safe(loc["hourly"]["wind_speed_10m"])
            dir_ = _safe(loc["hourly"]["wind_direction_10m"])
            rad  = np.radians(dir_)
            u_flat[start + k] = -spd * np.sin(rad)
            v_flat[start + k] = -spd * np.cos(rad)

    n_t = len(times_h)
    u_grid = u_flat.reshape(n_lat, n_lon, n_t).transpose(1, 0, 2)
    v_grid = v_flat.reshape(n_lat, n_lon, n_t).transpose(1, 0, 2)

    kw = dict(method="linear", bounds_error=False, fill_value=None)
    interp_u = RegularGridInterpolator((lons, lats, times_h), u_grid, **kw)
    interp_v = RegularGridInterpolator((lons, lats, times_h), v_grid, **kw)

    valid_dts = [datetime.fromisoformat(t) for t in times_iso]
    meta = {
        "valid_times": [t.strftime("%Y-%m-%dT%H:%M:%S+00:00") for t in valid_dts],
        "times":       times_h.tolist(),
        "bbox":        [float(lons[0]), float(lats[0]), float(lons[-1]), float(lats[-1])],
        "model":       cfg["label"],
        "days":        days,
    }
    return meta, interp_u, interp_v


# ── Interface publique ────────────────────────────────────────────────────────

def _get_cached(model_id: str):
    now = datetime.now(timezone.utc)
    if model_id in _cache:
        meta, iu, iv, fetched = _cache[model_id]
        if (now - fetched).total_seconds() < CACHE_AGE * 3600:
            return meta, iu, iv
    meta, iu, iv = _download(model_id)
    _cache[model_id] = (meta, iu, iv, now)
    return meta, iu, iv


def get_meta(model_id: str):
    meta, _, _ = _get_cached(model_id)
    return meta


def _uv(interp_u, interp_v, lons_arr, lats_arr, t):
    pts  = np.column_stack([lons_arr, lats_arr, np.full(len(lons_arr), t)])
    u    = np.nan_to_num(interp_u(pts))
    v    = np.nan_to_num(interp_v(pts))
    land = contains_xy(land_geom, lons_arr, lats_arr)
    u[land] = v[land] = 0.0
    return u, v


def _to_dir_force(u, v):
    force     = np.sqrt(u**2 + v**2) * 1.94384
    direction = (np.degrees(np.arctan2(u, v)) + 180) % 360
    return direction, force


def get_V_deg(model_id: str):
    _, iu, iv = _get_cached(model_id)

    def V(p, t):
        p     = np.asarray(p, dtype=float)
        batch = p.ndim == 2
        xs    = p[:, 0] if batch else p[0:1]
        ys    = p[:, 1] if batch else p[1:2]
        u, v  = _uv(iu, iv, xs, ys, t)
        d, f  = _to_dir_force(u, v)
        res   = np.column_stack([d, f])
        return res if batch else res[0]

    return V


def get_V_nm(model_id: str):
    _, iu, iv = _get_cached(model_id)

    def V(p, t):
        p     = np.asarray(p, dtype=float)
        batch = p.ndim == 2
        xs    = p[:, 0] if batch else p[0:1]
        ys    = p[:, 1] if batch else p[1:2]
        lons_arr = xs / (60.0 * 0.7)
        lats_arr = ys / 60.0
        u, v  = _uv(iu, iv, lons_arr, lats_arr, t)
        d, f  = _to_dir_force(u, v)
        res   = np.column_stack([d, f])
        return res if batch else res[0]

    return V
