from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pathlib import Path
from pydantic import BaseModel
from typing import List, Optional
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
from inputs.vents import table, vent_grib_nm, vent_grib_deg as load_vent_deg, vent_uniforme, land_geom  # noqa: E402
from shapely import contains_xy
from inputs.polaires import polaire as load_polaire   # noqa: E402
from core.isochrone import routage                    # noqa: E402
from inputs.courants import table as current_table, courant_grib_deg as load_courant_deg  # noqa: E402

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

GRIB_DIR         = ROOT / "data" / "grib_vent"
GRIB_COURANT_DIR = ROOT / "data" / "grib_courant"
POLAIRE_DIR      = ROOT / "data" / "polaires"

_df_cache:   dict = {}
_v_cache:    dict = {}
_vd_cache:   dict = {}
_grid_cache: dict = {}   # (filename, stride) → (la, lo)
_p_cache:    dict = {}

_cd_cache:    dict = {}  # current: filename → DataFrame
_cv_cache:    dict = {}  # current: grib_path → courant function
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


def _get_P(pol_path: str):
    """Fonction polaire interpolée (cache)."""
    if pol_path not in _p_cache:
        _p_cache[pol_path] = load_polaire(pol_path)
    return _p_cache[pol_path]


# ── endpoints vent ──────────────────────────────────────────────────────────

@app.get("/files")
def list_files():
    return sorted(f.name for f in GRIB_DIR.glob("*.grb2"))


@app.get("/wind/{filename}/meta")
def get_meta(filename: str):
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


# ── helpers courant ────────────────────────────────────────────────────────

def _load_current(filename: str):
    if filename not in _cd_cache:
        _cd_cache[filename] = current_table(str(GRIB_COURANT_DIR / filename))
    return _cd_cache[filename]


def _get_courant(grib_path: str):
    if grib_path not in _cv_cache:
        _cv_cache[grib_path] = load_courant_deg(grib_path)
    return _cv_cache[grib_path]


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
    times = sorted(float(t) for t in df["step_h"].unique())
    return {
        "times": times,
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
    grib_file:    Optional[str] = None
    polaire_file: str
    p_dep:   List[float]   # [lat, lon] en degrés -180/180
    p_arr:   List[float]
    t:       float = 0.0   # heure de départ (offset GRIB)
    dt:      float = 1.0
    n:       int   = 100
    ang_deg:  float = 90.0
    dang_deg: float = 0.3
    polar_pct: float = 100.0  # pourcentage de performance polaire (100 = nominal)
    wind_uniform: Optional[WindUniform] = None


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


def _build_isochrones(L, grib_file=None):
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
    pol_path = POLAIRE_DIR / req.polaire_file
    if not pol_path.exists():
        raise HTTPException(404, "Polaire introuvable")

    if req.wind_uniform:
        V = vent_uniforme(req.wind_uniform.direction, req.wind_uniform.force)
    else:
        if not req.grib_file:
            raise HTTPException(400, "grib_file requis si wind_uniform absent")
        grib_path = GRIB_DIR / req.grib_file
        if not grib_path.exists():
            raise HTTPException(404, "GRIB introuvable")
        V = _get_V(str(grib_path))
    P = _scale_polar(_get_P(str(pol_path)), req.polar_pct)

    p_dep = [req.p_dep[0], req.p_dep[1]]
    p_arr = [req.p_arr[0], req.p_arr[1]]

    if contains_xy(land_geom, req.p_dep[1], req.p_dep[0]) or \
       contains_xy(land_geom, req.p_arr[1], req.p_arr[0]):
        raise HTTPException(400, "Point à terre")

    t0 = time.perf_counter()
    lat, lon, time_list, L = routage(
        p_dep, p_arr, req.t,
        dt=req.dt, n=req.n, V=V, P=P,
        ang=math.radians(req.ang_deg),
        dang=math.radians(req.dang_deg),
    )
    calc_time_s = round(time.perf_counter() - t0, 2)

    route       = [[float(lo), float(la)] for lo, la in zip(lon.tolist(), lat.tolist())]
    route_times = time_list[:len(lat)].tolist()
    duration    = float(time_list[-1] - req.t)

    total_min = int(round(duration * 60))
    days      = total_min // (24 * 60)
    hours     = (total_min % (24 * 60)) // 60
    minutes   = total_min % 60

    isochrones = _build_isochrones(L, grib_file=req.grib_file)

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
    pol_path = POLAIRE_DIR / req.polaire_file
    if not pol_path.exists():
        raise HTTPException(404, "Polaire introuvable")

    if req.wind_uniform:
        V = vent_uniforme(req.wind_uniform.direction, req.wind_uniform.force)
    else:
        if not req.grib_file:
            raise HTTPException(400, "grib_file requis si wind_uniform absent")
        grib_path = GRIB_DIR / req.grib_file
        if not grib_path.exists():
            raise HTTPException(404, "GRIB introuvable")
        V = _get_V(str(grib_path))
    P = _scale_polar(_get_P(str(pol_path)), req.polar_pct)

    p_dep = [req.p_dep[0], req.p_dep[1]]
    p_arr = [req.p_arr[0], req.p_arr[1]]

    if contains_xy(land_geom, req.p_dep[1], req.p_dep[0]) or \
       contains_xy(land_geom, req.p_arr[1], req.p_arr[0]):
        raise HTTPException(400, "Point à terre")

    q = _queue.Queue()

    def _progress(pct):
        q.put_nowait({"type": "progress", "pct": pct})

    def _worker():
        try:
            t0 = time.perf_counter()
            lat, lon, time_list, L = routage(
                p_dep, p_arr, req.t,
                dt=req.dt, n=req.n, V=V, P=P,
                ang=math.radians(req.ang_deg),
                dang=math.radians(req.dang_deg),
                        progress_cb=_progress,
            )
            calc_time_s = round(time.perf_counter() - t0, 2)
            route       = [[float(lo), float(la)] for lo, la in zip(lon.tolist(), lat.tolist())]
            route_times = time_list[:len(lat)].tolist()
            duration    = float(time_list[-1] - req.t)
            total_min   = int(round(duration * 60))
            isochrones  = _build_isochrones(L, grib_file=req.grib_file)
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
