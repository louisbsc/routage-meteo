"""
Vent 10m IFS HRES 0.25° téléchargé directement depuis l'API Open Data d'ECMWF
(data.ecmwf.int), sans passer par Open-Meteo.

Chaque échéance de prévision est publiée en un fichier GRIB2 global (~140 Mo,
tous paramètres). On évite de le télécharger en entier : un fichier .index
(JSON) liste l'offset/longueur en octets de chaque paramètre dans le GRIB2,
ce qui permet de ne récupérer que les messages 10u/10v (~1,7 Mo/échéance) via
des requêtes HTTP Range.

Expose la même interface que vents_gfs : get_meta / get_V_deg / get_V_nm.

Usage :
    from inputs.vents_ecmwf import get_meta, get_V_deg, get_V_nm
"""

import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import eccodes
import numpy as np
import requests
from scipy.interpolate import RegularGridInterpolator
from shapely import contains_xy

from inputs.vents import land_geom

BASE_URL  = "https://data.ecmwf.int/forecasts"
RESOL     = "0p25"
STREAM    = "oper"
BBOX      = (-80.0, 20.0, 25.0, 75.0)   # lon0, lat0, lon1, lat1 (zone Atlantique/Europe)
HEADERS   = {"User-Agent": "Mozilla/5.0 (compatible; sailing-router/1.0)"}

PUBLISH_DELAY_H = 7      # marge avant qu'un run soit considéré comme publié
CACHE_AGE       = 6 * 3600
FAIL_COOLDOWN   = 120    # secondes avant de retenter après un échec
N_WORKERS       = 8      # téléchargements d'échéances en parallèle

_cache: dict = {}        # "ifs025" → (meta, interp_u, interp_v, fetched_at)
_fail_cache: dict = {}   # "ifs025" → (timestamp, message)


def _forecast_steps():
    """Échéances publiées par le stream oper IFS HRES : 3h jusqu'à 144h, puis 6h jusqu'à 360h (15j)."""
    return list(range(0, 145, 3)) + list(range(150, 361, 6))


def _latest_run():
    """Dernier run (00/06/12/18 UTC) dont la publication est probablement terminée."""
    candidate = datetime.now(timezone.utc) - timedelta(hours=PUBLISH_DELAY_H)
    hour = (candidate.hour // 6) * 6
    return candidate.replace(hour=hour, minute=0, second=0, microsecond=0)


def _url_base(run, step):
    d, h = run.strftime("%Y%m%d"), run.strftime("%H")
    return f"{BASE_URL}/{d}/{h}z/ifs/{RESOL}/{STREAM}/{d}{h}0000-{step}h-{STREAM}-fc"


def _get(session, url, headers=None, max_retries=4, timeout=60):
    for attempt in range(max_retries):
        if attempt:
            time.sleep(2 ** attempt)
        r = session.get(url, headers={**HEADERS, **(headers or {})}, timeout=timeout)
        if r.status_code == 429 and attempt < max_retries - 1:
            continue
        r.raise_for_status()
        return r


def _grid_bbox_indices():
    """Indices (dans la grille globale 1440×721) couvrant BBOX, lat/lon croissants."""
    lon0, lat0, lon1, lat1 = BBOX
    lats_full = 90.0 - np.arange(721) * 0.25
    lons_full = -180.0 + np.arange(1440) * 0.25
    lat_idx = np.where((lats_full >= lat0) & (lats_full <= lat1))[0][::-1]  # → lat croissante
    lon_idx = np.where((lons_full >= lon0) & (lons_full <= lon1))[0]
    return lat_idx, lon_idx, lats_full[lat_idx], lons_full[lon_idx]


def _fetch_step_uv(session, run, step, lat_idx, lon_idx):
    """Télécharge les messages 10u/10v d'une échéance (Range HTTP) et les recadre sur BBOX."""
    base = _url_base(run, step)
    idx_text = _get(session, base + ".index").text
    entries = {}
    for line in idx_text.splitlines():
        rec = json.loads(line)
        entries[rec["param"]] = rec

    out = []
    for param in ("10u", "10v"):
        entry = entries.get(param)
        if entry is None:
            raise RuntimeError(f"paramètre {param} absent de l'index ({base})")
        off, ln = entry["_offset"], entry["_length"]
        r = _get(session, base + ".grib2", headers={"Range": f"bytes={off}-{off + ln - 1}"})
        gid = eccodes.codes_new_from_message(r.content)
        try:
            Ni = eccodes.codes_get(gid, "Ni")
            Nj = eccodes.codes_get(gid, "Nj")
            vals = np.array(eccodes.codes_get_values(gid), dtype=float).reshape(Nj, Ni)
        finally:
            eccodes.codes_release(gid)
        out.append(vals[lat_idx, :][:, lon_idx])
    return step, out[0], out[1]


def _download():
    run   = _latest_run()
    steps = _forecast_steps()
    lat_idx, lon_idx, lats, lons = _grid_bbox_indices()

    u_by_step = {}
    v_by_step = {}
    with requests.Session() as session, ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = [pool.submit(_fetch_step_uv, session, run, s, lat_idx, lon_idx) for s in steps]
        for fut in futures:
            step, u, v = fut.result()
            u_by_step[step] = u
            v_by_step[step] = v

    n_lon, n_lat, n_t = len(lons), len(lats), len(steps)
    u_grid = np.empty((n_lon, n_lat, n_t))
    v_grid = np.empty((n_lon, n_lat, n_t))
    for k, step in enumerate(steps):
        u_grid[:, :, k] = u_by_step[step].T
        v_grid[:, :, k] = v_by_step[step].T

    times_h   = np.array(steps, dtype=float)
    valid_dts = [run + timedelta(hours=int(s)) for s in steps]
    meta = {
        "valid_times": [dt.strftime("%Y-%m-%dT%H:%M:%S+00:00") for dt in valid_dts],
        "times":       times_h.tolist(),
        "bbox":        [float(lons[0]), float(lats[0]), float(lons[-1]), float(lats[-1])],
        "model":       "ECMWF IFS HRES 0.25°",
        "days":        int(steps[-1] // 24),
    }

    kw = dict(method="linear", bounds_error=False, fill_value=None)
    interp_u = RegularGridInterpolator((lons, lats, times_h), u_grid, **kw)
    interp_v = RegularGridInterpolator((lons, lats, times_h), v_grid, **kw)
    return meta, interp_u, interp_v


def _get_cached():
    now    = datetime.now(timezone.utc)
    cached = _cache.get("ifs025")
    if cached and (now - cached[3]).total_seconds() < CACHE_AGE:
        return cached[0], cached[1], cached[2]

    fail = _fail_cache.get("ifs025")
    if fail and (now - fail[0]).total_seconds() < FAIL_COOLDOWN:
        raise RuntimeError(fail[1])

    try:
        meta, iu, iv = _download()
    except Exception as e:
        _fail_cache["ifs025"] = (now, str(e))
        raise
    _fail_cache.pop("ifs025", None)
    _cache["ifs025"] = (meta, iu, iv, now)
    return meta, iu, iv


def get_meta():
    meta, _, _ = _get_cached()
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


def get_V_deg():
    _, iu, iv = _get_cached()

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


def get_V_nm():
    _, iu, iv = _get_cached()

    def V(p, t):
        p        = np.asarray(p, dtype=float)
        batch    = p.ndim == 2
        xs       = p[:, 0] if batch else p[0:1]
        ys       = p[:, 1] if batch else p[1:2]
        lons_arr = xs / (60.0 * 0.7)
        lats_arr = ys / 60.0
        u, v     = _uv(iu, iv, lons_arr, lats_arr, t)
        d, f     = _to_dir_force(u, v)
        res      = np.column_stack([d, f])
        return res if batch else res[0]

    return V
