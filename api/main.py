from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pathlib import Path
from pydantic import BaseModel
from typing import List, Optional
from datetime import datetime, timedelta, timezone
import sys
import math
import time
import json
import queue as _queue
import threading
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
from inputs.vents import table, vent_grib_nm, vent_grib_deg as load_vent_deg, vent_uniforme, land_geom, land_geoms  # noqa: E402
from inputs.wind_fetcher import ecmwf as ecmwf_model, gfs as gfs_model  # noqa: E402
from shapely import contains_xy

# Modèles à téléchargement direct (ECMWF Open Data / NOAA NOMADS).
# Les anciennes clés "openmeteo_*" sont conservées pour la compatibilité frontend.
DIRECT_WIND_MODELS = {
    "ecmwf":           ecmwf_model,
    "gfs":             gfs_model,
    "openmeteo_ecmwf": ecmwf_model,
    "openmeteo_gfs":   gfs_model,
}
from inputs.polaires import polaire as load_polaire, polaire_uniforme   # noqa: E402
from core.isochrone import routage_def    # noqa: E402
from inputs.courants import (  # noqa: E402
    table as current_table,
    courant_grib_deg as load_courant_deg,
    courant_grib_nm  as load_courant_nm,
    courant_uniforme,
)

from apscheduler.schedulers.background import BackgroundScheduler  # noqa: E402

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

_scheduler = BackgroundScheduler(timezone="UTC")


@app.on_event("startup")
def _on_startup():
    def _boot():
        # Téléchargements initiaux en parallèle (réduit le temps de démarrage)
        threads = [
            threading.Thread(target=ecmwf_model.ensure_fresh, daemon=True),
            threading.Thread(target=gfs_model.ensure_fresh, daemon=True),
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        ecmwf_model.register_jobs(_scheduler)
        gfs_model.register_jobs(_scheduler)
        _scheduler.start()

    threading.Thread(target=_boot, daemon=True, name="wind-boot").start()


@app.on_event("shutdown")
def _on_shutdown():
    if _scheduler.running:
        _scheduler.shutdown(wait=False)

GRIB_DIR         = ROOT / "data" / "grib_vent"
GRIB_COURANT_DIR = ROOT / "data" / "grib_courant"
POLAIRE_DIR      = ROOT / "data" / "polaires"

_df_cache:   dict = {}
_v_cache:    dict = {}
_vd_cache:   dict = {}
_grid_cache: dict = {}   # (filename, stride) → (la, lo)
_sea_mask_cache: dict = {}  # (lat0r,lat1r,lon0r,lon1r,stepr) → (lo_f, la_f, sea)
_p_cache:    dict = {}

_cd_cache:    dict = {}  # current: filename → DataFrame
_cv_cache:    dict = {}  # current: grib_path → courant_grib_deg (display)
_cv_nm_cache: dict = {}  # current: grib_path → courant_grib_nm  (routing)
_cgrid_cache: dict = {}  # current: (filename, stride) → (la, lo)


# ── helpers ────────────────────────────────────────────────────────────────

def _load(filename: str):
    """DataFrame GRIB (longitude déjà normalisée -180/180 par table())."""
    if filename not in _df_cache:
        _df_cache[filename] = table(str(GRIB_DIR / filename))
    return _df_cache[filename]


def _get_V_deg(grib_path: str):
    """Fonction vent interpolée en degrés (cache), utilisée pour l'interpolation temporelle."""
    if grib_path not in _vd_cache:
        _vd_cache[grib_path] = load_vent_deg(grib_path)
    return _vd_cache[grib_path]


def _get_grid(filename: str, stride: int):
    """Grille de points (la, lo) en convention -180/180."""
    key = (filename, stride)
    if key not in _grid_cache:
        df = _load(filename)
        lats = sorted(df["latitude"].unique())
        lons = sorted(df["longitude"].unique())
        la, lo = np.meshgrid(
            np.array(lats[::stride], dtype=float),
            np.array(lons[::stride], dtype=float),
            indexing='ij',
        )
        _grid_cache[key] = (la.ravel(), lo.ravel())
    return _grid_cache[key]


def _get_V(grib_path: str):
    """Fonction vent interpolée (cache)."""
    if grib_path not in _v_cache:
        _v_cache[grib_path] = vent_grib_nm(grib_path)
    return _v_cache[grib_path]


def _get_sea_grid(lat0: float, lat1: float, lon0: float, lon1: float, step: float):
    """Grille mer en cache : évite contains_xy répété pour le même viewport+step."""
    key = (round(lat0, 2), round(lat1, 2), round(lon0, 2), round(lon1, 2), round(step, 3))
    if key in _sea_mask_cache:
        return _sea_mask_cache[key]
    lat_start = math.ceil(lat0 / step) * step
    lon_start = math.ceil(lon0 / step) * step
    lats = np.arange(lat_start, lat1 + step / 2, step)
    lons = np.arange(lon_start, lon1 + step / 2, step)
    if lats.size == 0 or lons.size == 0:
        return None
    lo, la = np.meshgrid(lons, lats)
    la_f = la.ravel().astype(float)
    lo_f = lo.ravel().astype(float)
    sea  = ~contains_xy(land_geom, lo_f, la_f)
    _sea_mask_cache[key] = (lo_f, la_f, sea)
    if len(_sea_mask_cache) > 80:
        del _sea_mask_cache[next(iter(_sea_mask_cache))]
    return lo_f, la_f, sea


def _get_P(pol_path: str):
    """Fonction polaire interpolée (cache)."""
    if pol_path not in _p_cache:
        _p_cache[pol_path] = load_polaire(pol_path)
    return _p_cache[pol_path]


# ── endpoints vent ──────────────────────────────────────────────────────────

def _om_call(fn, *args):
    """Appelle une fonction de téléchargement vent et convertit les échecs réseau en 503 explicite."""
    try:
        return fn(*args)
    except Exception as e:
        raise HTTPException(503, f"Source de vent indisponible : {e}")


@app.get("/files")
def list_files():
    return sorted(f.name for f in GRIB_DIR.glob("*.grb2")) + list(DIRECT_WIND_MODELS)


@app.get("/wind/{filename}/meta")
def get_meta(filename: str):
    if filename in DIRECT_WIND_MODELS:
        return _om_call(DIRECT_WIND_MODELS[filename].get_meta)
    if not (GRIB_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    df = _load(filename)
    step_df = (
        df[["temps", "valid_time"]]
        .drop_duplicates("temps")
        .sort_values("temps")
    )
    times = step_df["temps"].tolist()

    def _iso(ts):
        t = pd.Timestamp(ts)
        return (t if t.tzinfo else t.tz_localize("UTC")).isoformat()

    valid_times = [_iso(vt) for vt in step_df["valid_time"]]
    return {
        "times": times,
        "valid_times": valid_times,
        "bbox": [
            float(df["longitude"].min()),
            float(df["latitude"].min()),
            float(df["longitude"].max()),
            float(df["latitude"].max()),
        ],
    }


@app.get("/wind/{filename}/step/{time_idx}")
def get_step(filename: str, time_idx: int, stride: int = 2):
    if not (GRIB_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    df = _load(filename)
    times = sorted(df["temps"].unique())
    if not (0 <= time_idx < len(times)):
        raise HTTPException(400, "Invalid time index")

    t = times[time_idx]
    df_t = df[df["temps"] == t].dropna(subset=["u10", "v10"])

    lats = sorted(df_t["latitude"].unique())
    lons = sorted(df_t["longitude"].unique())
    df_t = df_t[
        df_t["latitude"].isin(set(lats[::stride])) &
        df_t["longitude"].isin(set(lons[::stride]))
    ]

    cols = df_t[["latitude", "longitude", "force", "direction"]].copy()
    cols.columns = ["lat", "lon", "speed", "dir"]
    cols = cols.round({"lat": 4, "lon": 4, "speed": 2, "dir": 1})
    return {"time_h": float(t), "data": cols.to_dict(orient="records")}


@app.get("/wind/{filename}/interpolated")
def get_interpolated(filename: str, t: float, stride: int = 2):
    """Vent interpolé linéairement à l'heure t (peut être entre deux pas GRIB)."""
    if not (GRIB_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    V = _get_V_deg(str(GRIB_DIR / filename))
    la, lo = _get_grid(filename, stride)
    wind = V(np.column_stack([lo, la]), t)   # [[dir, force], …]
    df_out = pd.DataFrame({
        "lat":   np.round(la, 4),
        "lon":   np.round(lo, 4),
        "dir":   np.round(wind[:, 0], 1),
        "speed": np.round(wind[:, 1], 2),
    })
    return {"time_h": float(t), "data": df_out.to_dict(orient="records")}


# ── endpoint grille vent uniforme ───────────────────────────────────────────

@app.get("/wind/uniform/grid")
def uniform_grid(
    lat0: float, lat1: float, lon0: float, lon1: float,
    step: float = 1.0,
    direction: float = 270,
    force: float = 15,
):
    lat_start = math.ceil(lat0 / step) * step
    lon_start = math.ceil(lon0 / step) * step
    lats = np.arange(lat_start, lat1 + step / 2, step)
    lons = np.arange(lon_start, lon1 + step / 2, step)
    if lats.size == 0 or lons.size == 0:
        return {"data": []}
    lo, la = np.meshgrid(lons, lats)
    la_f = la.ravel().astype(float)
    lo_f = lo.ravel().astype(float)
    sea = ~contains_xy(land_geom, lo_f, la_f)
    return {"data": [
        {"lat": round(float(la_f[i]), 4), "lon": round(float(lo_f[i]), 4),
         "dir": direction, "speed": force}
        for i in np.where(sea)[0]
    ]}


# ── endpoint grille vent GRIB viewport (après uniform/grid pour éviter conflit) ──

@app.get("/wind/{filename}/grid")
def get_wind_grid(filename: str, t: float, lat0: float, lat1: float, lon0: float, lon1: float, step: float = 1.0):
    if filename in DIRECT_WIND_MODELS:
        model = DIRECT_WIND_MODELS[filename]
        V     = _om_call(model.get_V_deg)
        mbbox = _om_call(model.get_meta)["bbox"]
        clat0 = max(lat0, mbbox[1]);  clat1 = min(lat1, mbbox[3])
        clon0 = max(lon0, mbbox[0]);  clon1 = min(lon1, mbbox[2])
        grid = _get_sea_grid(clat0, clat1, clon0, clon1, step)
        if grid is None:
            return {"lat": [], "lon": [], "dir": [], "speed": []}
        lo_f, la_f, sea = grid
        wind = V(np.column_stack([lo_f[sea], la_f[sea]]), t)
        return {
            "time_h": float(t),
            "lat":   np.round(la_f[sea], 4).tolist(),
            "lon":   np.round(lo_f[sea], 4).tolist(),
            "dir":   np.round(wind[:, 0], 1).tolist(),
            "speed": np.round(wind[:, 1], 2).tolist(),
        }
    if not (GRIB_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    V = _get_V_deg(str(GRIB_DIR / filename))
    grid = _get_sea_grid(lat0, lat1, lon0, lon1, step)
    if grid is None:
        return {"lat": [], "lon": [], "dir": [], "speed": []}
    lo_f, la_f, sea = grid
    df = _load(filename)
    in_bbox = (
        (lo_f >= float(df["longitude"].min())) & (lo_f <= float(df["longitude"].max())) &
        (la_f >= float(df["latitude"].min()))  & (la_f <= float(df["latitude"].max()))
    )
    keep = in_bbox & sea
    if not keep.any():
        return {"lat": [], "lon": [], "dir": [], "speed": []}
    wind = V(np.column_stack([lo_f[keep], la_f[keep]]), t)
    return {
        "time_h": float(t),
        "lat":   np.round(la_f[keep], 4).tolist(),
        "lon":   np.round(lo_f[keep], 4).tolist(),
        "dir":   np.round(wind[:, 0], 1).tolist(),
        "speed": np.round(wind[:, 1], 2).tolist(),
    }


# ── helpers courant ────────────────────────────────────────────────────────

def _load_current(filename: str):
    if filename not in _cd_cache:
        _cd_cache[filename] = current_table(str(GRIB_COURANT_DIR / filename))
    return _cd_cache[filename]


def _get_courant(grib_path: str):
    if grib_path not in _cv_cache:
        _cv_cache[grib_path] = load_courant_deg(grib_path)
    return _cv_cache[grib_path]


def _get_courant_nm(grib_path: str):
    """Courant interpolé en coordonnées NM (pour le routage)."""
    if grib_path not in _cv_nm_cache:
        _cv_nm_cache[grib_path] = load_courant_nm(grib_path)
    return _cv_nm_cache[grib_path]


def _get_current_grid(filename: str, stride: int):
    key = (filename, stride)
    if key not in _cgrid_cache:
        df = _load_current(filename)
        lats = sorted(df["lat"].unique())
        lons = sorted(df["lon"].unique())
        la, lo = np.meshgrid(
            np.array(lats[::stride], dtype=float),
            np.array(lons[::stride], dtype=float),
            indexing='ij',
        )
        _cgrid_cache[key] = (la.ravel(), lo.ravel())
    return _cgrid_cache[key]


# ── endpoints courant ───────────────────────────────────────────────────────

@app.get("/current-files")
def list_current_files():
    return sorted(
        f.name for f in GRIB_COURANT_DIR.iterdir()
        if f.suffix in ('.grb', '.grb2', '.grib', '.grib2')
    )


@app.get("/current/{filename}/meta")
def get_current_meta(filename: str):
    if not (GRIB_COURANT_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    df = _load_current(filename)

    # Datetime de référence : date + heure du run
    date_int = int(df["date"].iloc[0])
    time_int = int(df["run_time"].iloc[0]) if "run_time" in df.columns else 0
    year   = date_int // 10000
    month  = (date_int % 10000) // 100
    day    = date_int % 100
    hour   = time_int // 100
    minute = time_int % 100
    ref_dt = datetime(year, month, day, hour, minute, tzinfo=timezone.utc)

    times = sorted(float(t) for t in df["step_h"].unique())

    def _iso(step_h):
        dt = ref_dt + timedelta(hours=step_h)
        return dt.isoformat()

    return {
        "times":       times,
        "valid_times": [_iso(t) for t in times],
        "bbox": [
            float(df["lon"].min()),
            float(df["lat"].min()),
            float(df["lon"].max()),
            float(df["lat"].max()),
        ],
    }


@app.get("/current/{filename}/interpolated")
def get_current_interpolated(filename: str, t: float, stride: int = 2):
    if not (GRIB_COURANT_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    C = _get_courant(str(GRIB_COURANT_DIR / filename))
    la, lo = _get_current_grid(filename, stride)
    current = C(np.column_stack([lo, la]), t)  # [[dir, force], ...]
    df_out = pd.DataFrame({
        "lat":   np.round(la, 4),
        "lon":   np.round(lo, 4),
        "dir":   np.round(current[:, 0], 1),
        "speed": np.round(current[:, 1], 3),
    })
    df_out = df_out[df_out["speed"] > 0.01]
    return {"time_h": float(t), "data": df_out.to_dict(orient="records")}


# ── endpoint grille courant GRIB viewport ─────────────────────────────────

@app.get("/current/{filename}/grid")
def get_current_grid_view(filename: str, t: float, lat0: float, lat1: float, lon0: float, lon1: float, step: float = 1.0):
    if not (GRIB_COURANT_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    C = _get_courant(str(GRIB_COURANT_DIR / filename))
    grid = _get_sea_grid(lat0, lat1, lon0, lon1, step)
    if grid is None:
        return {"lat": [], "lon": [], "dir": [], "speed": []}
    lo_f, la_f, sea = grid
    df = _load_current(filename)
    in_bbox = (
        (lo_f >= float(df["lon"].min())) & (lo_f <= float(df["lon"].max())) &
        (la_f >= float(df["lat"].min())) & (la_f <= float(df["lat"].max()))
    )
    keep = in_bbox & sea
    if not keep.any():
        return {"lat": [], "lon": [], "dir": [], "speed": []}
    current = C(np.column_stack([lo_f[keep], la_f[keep]]), t)
    fast    = current[:, 1] > 0.01
    return {
        "time_h": float(t),
        "lat":   np.round(la_f[keep][fast], 4).tolist(),
        "lon":   np.round(lo_f[keep][fast], 4).tolist(),
        "dir":   np.round(current[fast, 0], 1).tolist(),
        "speed": np.round(current[fast, 1], 3).tolist(),
    }


# ── cartographie terrestre ──────────────────────────────────────────────────

@app.get("/land/geojson")
def get_land_geojson(lat0: float, lat1: float, lon0: float, lon1: float, resolution: str = '10m'):
    from shapely.geometry import box, mapping
    geom = land_geoms.get(resolution, land_geom)
    clip = box(lon0, lat0, lon1, lat1)
    clipped = geom.intersection(clip)
    if clipped.is_empty:
        return {"type": "FeatureCollection", "features": []}
    geoms = list(clipped.geoms) if hasattr(clipped, 'geoms') else [clipped]
    features = [{"type": "Feature", "geometry": mapping(g), "properties": {}} for g in geoms]
    return {"type": "FeatureCollection", "features": features}


# ── endpoints routage ───────────────────────────────────────────────────────

@app.get("/polaires")
def list_polaires():
    return sorted(f.name for f in POLAIRE_DIR.glob("*.csv"))


@app.get("/polaires/{filename}/stats")
def get_polaire_stats(filename: str):
    pol_path = POLAIRE_DIR / filename
    if not pol_path.exists():
        raise HTTPException(404, "Polaire introuvable")
    df = pd.read_csv(pol_path)
    speeds = df.iloc[:, 1:].astype(float).values
    nonzero = speeds[speeds > 0]
    return {"v_mean": round(float(nonzero.mean()), 2)}


class WindUniform(BaseModel):
    direction: float   # degrés (convention météo : d'où vient le vent)
    force:     float   # nœuds

class RoutingRequest(BaseModel):
    grib_file:          Optional[str] = None
    polaire_file:       Optional[str] = None
    motor_speed:        Optional[float] = None
    p_dep:   List[float]   # [lat, lon] en degrés -180/180
    p_arr:   List[float]
    t:       float = 0.0   # heure de départ (offset GRIB)
    dt:      float = 1.0
    n:       int   = 100
    ang_deg:  float = 90.0
    dang_deg: float = 0.3
    polar_pct: float = 100.0
    wind_uniform:    Optional[WindUniform] = None
    grib_courant_file: Optional[str] = None
    courant_uniform:   Optional[WindUniform] = None  # même structure direction/force
    seuil_nm:          Optional[float] = None   # None → 50 % de la distance départ-arrivée
    facteur_raf:       float = 2.0


def _build_polar(req: "RoutingRequest"):
    if req.motor_speed is not None:
        return polaire_uniforme(req.motor_speed)
    if req.polaire_file:
        pol_path = POLAIRE_DIR / req.polaire_file
        if not pol_path.exists():
            raise HTTPException(404, "Polaire introuvable")
        return _scale_polar(_get_P(str(pol_path)), req.polar_pct)
    raise HTTPException(400, "polaire_file ou motor_speed requis")


def _build_courant(req: "RoutingRequest"):
    """Construit la fonction courant C(p, t) pour le routage, ou None si absent."""
    if req.courant_uniform:
        return courant_uniforme(req.courant_uniform.direction, req.courant_uniform.force)
    if req.grib_courant_file:
        path = GRIB_COURANT_DIR / req.grib_courant_file
        if not path.exists():
            raise HTTPException(404, "GRIB courant introuvable")
        return _get_courant_nm(str(path))
    return None


def _scale_polar(P, pct: float):
    """Retourne une fonction polaire dont les vitesses sont multipliées par pct/100."""
    if pct == 100.0:
        return P
    factor = pct / 100.0
    base = P
    def scaled(ang, v, _f=factor, _b=base):
        return _b(ang, v) * _f
    scaled.v_max = base.v_max * factor
    return scaled


def _build_isochrones(L, grib_file=None, break_inactive=False):
    """Construit les isochrones filtrées (terre + bbox GRIB) en segments contigus."""
    lon_min, lat_min, lon_max, lat_max = -180.0, -90.0, 180.0, 90.0
    if grib_file:
        df = _load(grib_file)
        lon_min = float(df["longitude"].min())
        lat_min = float(df["latitude"].min())
        lon_max = float(df["longitude"].max())
        lat_max = float(df["latitude"].max())
    isochrones = []
    for idx in np.unique(L[:, 2].astype(int)):
        pts  = L[L[:, 2].astype(int) == idx]
        lats = pts[:, 0]
        lons = pts[:, 1]
        in_bbox = (lons >= lon_min) & (lons <= lon_max) & \
                  (lats >= lat_min) & (lats <= lat_max)
        on_land = np.zeros(len(pts), dtype=bool)
        if in_bbox.any():
            on_land[in_bbox] = contains_xy(land_geom, lons[in_bbox], lats[in_bbox])
        keep = in_bbox & ~on_land
        if break_inactive:
            keep &= (pts[:, 3] >= 0)
        # Découpe en segments contigus : un point supprimé brise le tracé
        segment = []
        for i in range(len(pts)):
            if keep[i]:
                segment.append([round(float(lons[i]), 4), round(float(lats[i]), 4)])
            else:
                if len(segment) >= 2:
                    isochrones.append(segment)
                segment = []
        if len(segment) >= 2:
            isochrones.append(segment)
    return isochrones


@app.post("/routing")
def run_routing(req: RoutingRequest):
    if req.wind_uniform:
        V = vent_uniforme(req.wind_uniform.direction, req.wind_uniform.force)
    elif req.grib_file in DIRECT_WIND_MODELS:
        V = _om_call(DIRECT_WIND_MODELS[req.grib_file].get_V_nm)
    else:
        if not req.grib_file:
            raise HTTPException(400, "grib_file requis si wind_uniform absent")
        grib_path = GRIB_DIR / req.grib_file
        if not grib_path.exists():
            raise HTTPException(404, "GRIB introuvable")
        V = _get_V(str(grib_path))
    P = _build_polar(req)

    p_dep = [req.p_dep[0], req.p_dep[1]]
    p_arr = [req.p_arr[0], req.p_arr[1]]

    if contains_xy(land_geom, req.p_dep[1], req.p_dep[0]) or \
       contains_xy(land_geom, req.p_arr[1], req.p_arr[0]):
        raise HTTPException(400, "Point à terre")

    C = _build_courant(req)

    t0 = time.perf_counter()
    if req.seuil_nm is None:
        x_dep = p_dep[1] * 60 * 0.7;  y_dep = p_dep[0] * 60
        x_arr = p_arr[1] * 60 * 0.7;  y_arr = p_arr[0] * 60
        seuil = 0.5 * math.sqrt((x_arr - x_dep) ** 2 + (y_arr - y_dep) ** 2)
    else:
        seuil = req.seuil_nm
    lat, lon, time_list, L = routage_def(
        p_dep, p_arr, req.t,
        dt=req.dt, n=req.n, V=V, P=P,
        ang=math.radians(req.ang_deg),
        dang=math.radians(req.dang_deg),
        seuil=seuil,
        facteur_raf=req.facteur_raf,
        C=C,
    )
    calc_time_s = round(time.perf_counter() - t0, 2)

    route       = [[float(lo), float(la)] for lo, la in zip(lon.tolist(), lat.tolist())]
    route_times = time_list[:len(lat)].tolist()
    duration    = float(time_list[-1] - req.t)

    total_min = int(round(duration * 60))
    days      = total_min // (24 * 60)
    hours     = (total_min % (24 * 60)) // 60
    minutes   = total_min % 60

    iso_grib = None if req.grib_file in DIRECT_WIND_MODELS else req.grib_file
    isochrones = _build_isochrones(L, grib_file=iso_grib, break_inactive=False)

    return {
        "route":        route,
        "time_list":    route_times,
        "isochrones":   isochrones,
        "duration_h":   duration,
        "days":         days,
        "hours":        hours,
        "minutes":      minutes,
        "calc_time_s":  calc_time_s,
    }


@app.post("/routing/stream")
def run_routing_stream(req: RoutingRequest):
    if req.wind_uniform:
        V = vent_uniforme(req.wind_uniform.direction, req.wind_uniform.force)
    elif req.grib_file in DIRECT_WIND_MODELS:
        V = _om_call(DIRECT_WIND_MODELS[req.grib_file].get_V_nm)
    else:
        if not req.grib_file:
            raise HTTPException(400, "grib_file requis si wind_uniform absent")
        grib_path = GRIB_DIR / req.grib_file
        if not grib_path.exists():
            raise HTTPException(404, "GRIB introuvable")
        V = _get_V(str(grib_path))
    P = _build_polar(req)

    p_dep = [req.p_dep[0], req.p_dep[1]]
    p_arr = [req.p_arr[0], req.p_arr[1]]

    if contains_xy(land_geom, req.p_dep[1], req.p_dep[0]) or \
       contains_xy(land_geom, req.p_arr[1], req.p_arr[0]):
        raise HTTPException(400, "Point à terre")

    C = _build_courant(req)

    q = _queue.Queue()

    def _progress(pct):
        q.put_nowait({"type": "progress", "pct": pct})

    def _worker():
        try:
            t0 = time.perf_counter()
            if req.seuil_nm is None:
                x_dep = p_dep[1] * 60 * 0.7;  y_dep = p_dep[0] * 60
                x_arr = p_arr[1] * 60 * 0.7;  y_arr = p_arr[0] * 60
                seuil = 0.5 * math.sqrt((x_arr - x_dep) ** 2 + (y_arr - y_dep) ** 2)
            else:
                seuil = req.seuil_nm
            lat, lon, time_list, L = routage_def(
                p_dep, p_arr, req.t,
                dt=req.dt, n=req.n, V=V, P=P,
                ang=math.radians(req.ang_deg),
                dang=math.radians(req.dang_deg),
                seuil=seuil,
                facteur_raf=req.facteur_raf,
                C=C,
                progress_cb=_progress,
            )
            calc_time_s = round(time.perf_counter() - t0, 2)
            route       = [[float(lo), float(la)] for lo, la in zip(lon.tolist(), lat.tolist())]
            route_times = time_list[:len(lat)].tolist()
            duration    = float(time_list[-1] - req.t)
            total_min   = int(round(duration * 60))
            iso_grib    = None if req.grib_file in DIRECT_WIND_MODELS else req.grib_file
            isochrones  = _build_isochrones(L, grib_file=iso_grib, break_inactive=False)
            q.put_nowait({"type": "result",
                "route": route, "time_list": route_times, "isochrones": isochrones,
                "duration_h": duration,
                "days":    total_min // (24 * 60),
                "hours":   (total_min % (24 * 60)) // 60,
                "minutes": total_min % 60,
                "calc_time_s": calc_time_s,
            })
        except Exception as e:
            q.put_nowait({"type": "error", "detail": str(e)})
        finally:
            q.put_nowait(None)

    threading.Thread(target=_worker, daemon=True).start()

    def _stream():
        while True:
            msg = q.get()
            if msg is None:
                break
            yield json.dumps(msg) + '\n'

    return StreamingResponse(_stream(), media_type="application/x-ndjson")
