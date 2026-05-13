from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from pydantic import BaseModel
from typing import List
import sys
import math
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
from inputs.vents import table, vent_grib_nm, vent_grib_deg as load_vent_deg  # noqa: E402
from inputs.polaires import polaire as load_polaire   # noqa: E402
from core.isochrone import routage                    # noqa: E402

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

GRIB_DIR    = ROOT / "data" / "grib_vent"
POLAIRE_DIR = ROOT / "data" / "polaires"

_df_cache:   dict = {}
_uses_360:   dict = {}   # filename → bool (convention lon GRIB)
_v_cache:    dict = {}
_vd_cache:   dict = {}   # vent_grib_deg cache
_grid_cache: dict = {}   # (filename, stride) → (la, lo_norm, lo_native)
_p_cache:    dict = {}


# ── helpers ────────────────────────────────────────────────────────────────

def _load(filename: str):
    """DataFrame GRIB avec longitude normalisée -180/180 (pour affichage)."""
    if filename not in _df_cache:
        df = table(str(GRIB_DIR / filename))
        _uses_360[filename] = bool((df["longitude"] > 180).any())
        df["longitude"] = df["longitude"].apply(lambda x: x - 360 if x > 180 else x)
        _df_cache[filename] = df
    return _df_cache[filename]


def _get_V_deg(grib_path: str):
    """Fonction vent interpolée en degrés (cache), utilisée pour l'interpolation temporelle."""
    if grib_path not in _vd_cache:
        _vd_cache[grib_path] = load_vent_deg(grib_path)
    return _vd_cache[grib_path]


def _get_grid(filename: str, stride: int):
    """Grille de points (la, lo_norm, lo_native) pour l'interpolation (cache)."""
    key = (filename, stride)
    if key not in _grid_cache:
        df = _load(filename)
        lats = sorted(df["latitude"].unique())
        lons = sorted(df["longitude"].unique())   # -180/180
        la, lo = np.meshgrid(
            np.array(lats[::stride], dtype=float),
            np.array(lons[::stride], dtype=float),
            indexing='ij',
        )
        la_flat = la.ravel()
        lo_flat = lo.ravel()                              # -180/180 (réponse JSON)
        # Convention native du GRIB pour interroger l'interpolateur
        lo_native = lo_flat % 360 if _uses_360.get(filename, True) else lo_flat
        _grid_cache[key] = (la_flat, lo_flat, lo_native)
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
    la, lo_norm, lo_native = _get_grid(filename, stride)
    wind = V(np.column_stack([lo_native, la]), t)   # [[dir, force], …]
    df_out = pd.DataFrame({
        "lat":   np.round(la, 4),
        "lon":   np.round(lo_norm, 4),
        "dir":   np.round(wind[:, 0], 1),
        "speed": np.round(wind[:, 1], 2),
    })
    return {"time_h": float(t), "data": df_out.to_dict(orient="records")}


# ── endpoints routage ───────────────────────────────────────────────────────

@app.get("/polaires")
def list_polaires():
    return sorted(f.name for f in POLAIRE_DIR.glob("*.csv"))


class RoutingRequest(BaseModel):
    grib_file:    str
    polaire_file: str
    p_dep:   List[float]   # [lat, lon] en degrés -180/180
    p_arr:   List[float]
    t:       float = 0.0   # heure de départ (offset GRIB)
    dt:      float = 1.0
    n:       int   = 100
    e_arr:   float = 20.0
    r:       float = 2.0
    ang_deg:  float = 90.0
    dang_deg: float = 0.3
    delta:   float = 2.0


@app.post("/routing")
def run_routing(req: RoutingRequest):
    grib_path = GRIB_DIR / req.grib_file
    pol_path  = POLAIRE_DIR / req.polaire_file

    if not grib_path.exists():
        raise HTTPException(404, "GRIB introuvable")
    if not pol_path.exists():
        raise HTTPException(404, "Polaire introuvable")

    V = _get_V(str(grib_path))
    P = _get_P(str(pol_path))

    # routage() attend la longitude en convention 0-360
    p_dep = [req.p_dep[0], req.p_dep[1] % 360]
    p_arr = [req.p_arr[0], req.p_arr[1] % 360]

    lat, lon, time_list, L = routage(
        p_dep, p_arr, req.t,
        dt=req.dt, n=req.n, V=V, P=P,
        e_arr=req.e_arr, r=req.r,
        ang=math.radians(req.ang_deg),
        dang=math.radians(req.dang_deg),
        delta=req.delta,
    )

    lon_norm    = [x - 360 if x > 180 else x for x in lon.tolist()]
    route       = [[lo, la] for lo, la in zip(lon_norm, lat.tolist())]
    route_times = time_list[:len(lat)].tolist()
    duration    = float(time_list[-1] - req.t)

    # Isochrones : L[:, 0]=lat, L[:, 1]=lon_0_360, L[:, 2]=index_iso
    # Les points de chaque isochrone sont déjà ordonnés par l'enveloppe
    isochrones = []
    for idx in np.unique(L[:, 2].astype(int)):
        pts = L[L[:, 2].astype(int) == idx]
        lo_norm_iso = [x - 360 if x > 180 else x for x in pts[:, 1]]
        isochrones.append([
            [round(float(lo), 4), round(float(la), 4)]
            for lo, la in zip(lo_norm_iso, pts[:, 0])
        ])

    return {
        "route":       route,
        "time_list":   route_times,
        "isochrones":  isochrones,
        "duration_h":  duration,
        "days":        int(duration // 24),
        "hours":       int(duration % 24),
    }
