"""
Vent 10m GFS 0.25° téléchargé directement depuis NOMADS (NOAA), sans passer
par Open-Meteo.

Comme pour vents_ecmwf.py, chaque échéance est publiée en GRIB2 global
(~500 Mo, tous paramètres). NOMADS expose un fichier .idx (texte, format
"n:offset:date:param:level:plage:") qui donne l'offset en octets de chaque
message ; on ne télécharge donc que UGRD/VGRD à 10m (~1,1 Mo/échéance) via
des requêtes HTTP Range.

Particularité GFS : le pas de temps natif est horaire de 0 à 120h, puis
3-horaire de 123h à 384h (16 jours) — contrairement à l'IFS (3h/6h).

Expose la même interface que vents_ecmwf : get_meta / get_V_deg / get_V_nm.
"""

import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone

import eccodes
import numpy as np
import requests
from scipy.interpolate import RegularGridInterpolator
from shapely import contains_xy

from inputs.vents import land_geom

BASE_URL = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
BBOX     = (-80.0, 20.0, 25.0, 75.0)   # lon0, lat0, lon1, lat1 (même zone que vents_ecmwf)
HEADERS  = {"User-Agent": "Mozilla/5.0 (compatible; sailing-router/1.0)"}

PUBLISH_DELAY_H  = 5     # marge avant qu'un run soit considéré comme publié
MAX_RUN_FALLBACK = 3     # nb de runs précédents tentés si le plus récent n'est pas dispo
CACHE_AGE        = 6 * 3600
FAIL_COOLDOWN    = 120
BLOCKED_COOLDOWN = 600   # NOMADS bloque plus longtemps qu'une simple erreur transitoire
N_WORKERS        = 3     # NOMADS est un serveur public partagé : on reste très modéré
MIN_REQUEST_GAP  = 0.4   # secondes minimum entre deux requêtes (toutes confondues)


class _NomadsBlocked(RuntimeError):
    """NOMADS répond par sa page « Over Rate Limit » (déguisée en 200/302)."""


def _forecast_steps():
    """Échéances publiées par GFS 0.25° : horaire 0→120h, puis 3h jusqu'à 384h (16j)."""
    return list(range(0, 121, 1)) + list(range(123, 385, 3))


_cache: dict = {}        # "gfs025" → (meta, interp_u, interp_v, fetched_at)
_fail_cache: dict = {}   # "gfs025" → (timestamp, message, cooldown_s)


class _RateLimiter:
    """Espace les requêtes HTTP dans le temps, tous threads confondus."""

    def __init__(self, min_interval):
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._last = 0.0

    def wait(self):
        with self._lock:
            now = time.monotonic()
            delay = self._last + self.min_interval - now
            if delay > 0:
                time.sleep(delay)
            self._last = time.monotonic()


_rate_limiter = _RateLimiter(MIN_REQUEST_GAP)


def _run_candidates():
    """Runs (00/06/12/18 UTC) du plus récent probable au plus ancien, par ordre de tentative."""
    now = datetime.now(timezone.utc) - timedelta(hours=PUBLISH_DELAY_H)
    hour = (now.hour // 6) * 6
    latest = now.replace(hour=hour, minute=0, second=0, microsecond=0)
    return [latest - timedelta(hours=6 * k) for k in range(MAX_RUN_FALLBACK + 1)]


def _file_base(run, step):
    d, h = run.strftime("%Y%m%d"), run.strftime("%H")
    return f"{BASE_URL}/gfs.{d}/{h}/atmos/gfs.t{h}z.pgrb2.0p25.f{step:03d}"


def _get(session, url, headers=None, max_retries=4, timeout=60):
    for attempt in range(max_retries):
        if attempt:
            time.sleep(2 ** attempt)
        _rate_limiter.wait()
        r = session.get(url, headers={**HEADERS, **(headers or {})}, timeout=timeout)
        # NOMADS renvoie sa page « Over Rate Limit » en HTML avec un statut 200/302
        # (redirection suivie automatiquement par requests) : on la détecte explicitement
        # plutôt que de la confondre avec une réponse .idx/.grib2 valide.
        if r.history or "text/html" in r.headers.get("Content-Type", ""):
            raise _NomadsBlocked(
                "NOMADS a limité le débit (Over Rate Limit) — réessayer dans quelques minutes"
            )
        if r.status_code in (429, 503) and attempt < max_retries - 1:
            continue
        r.raise_for_status()
        return r


def _pick_run(session):
    """Premier run de _run_candidates() dont le f000 répond (200)."""
    for run in _run_candidates():
        url = _file_base(run, 0) + ".idx"
        try:
            _get(session, url, max_retries=1)
            return run
        except _NomadsBlocked:
            raise
        except requests.RequestException:
            continue
    raise RuntimeError("Aucun run GFS récent n'est accessible sur NOMADS")


def _lon_indices_0_360(lon0, lon1, n=1440, res=0.25):
    """Indices de la grille (convention 0..360°) couvrant [lon0, lon1] signé (-180..180),
    avec gestion du passage par 0° (notre BBOX va de -80° à 25°)."""
    raws = np.arange(n) * res
    lon0_raw, lon1_raw = lon0 % 360, lon1 % 360
    if lon0_raw <= lon1_raw:
        idx = np.where((raws >= lon0_raw) & (raws <= lon1_raw))[0]
        signed = np.where(raws[idx] > 180, raws[idx] - 360, raws[idx])
    else:
        idx1 = np.where(raws >= lon0_raw)[0]
        idx2 = np.where(raws <= lon1_raw)[0]
        idx = np.concatenate([idx1, idx2])
        signed = np.concatenate([raws[idx1] - 360, raws[idx2]])
    return idx, signed


def _grid_bbox_indices():
    lon0, lat0, lon1, lat1 = BBOX
    lats_full = 90.0 - np.arange(721) * 0.25
    lat_idx = np.where((lats_full >= lat0) & (lats_full <= lat1))[0][::-1]  # → lat croissante
    lon_idx, lons = _lon_indices_0_360(lon0, lon1)
    return lat_idx, lon_idx, lats_full[lat_idx], lons


def _parse_idx(text):
    entries = []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split(":")
        entries.append({"offset": int(parts[1]), "param": parts[3], "level": parts[4]})
    return entries


def _decode_grib(message, lat_idx, lon_idx):
    gid = eccodes.codes_new_from_message(message)
    try:
        Ni = eccodes.codes_get(gid, "Ni")
        Nj = eccodes.codes_get(gid, "Nj")
        vals = np.array(eccodes.codes_get_values(gid), dtype=float).reshape(Nj, Ni)
    finally:
        eccodes.codes_release(gid)
    return vals[lat_idx, :][:, lon_idx]


def _fetch_step_uv_once(session, run, step, lat_idx, lon_idx):
    """Une seule requête Range couvrant UGRD+VGRD (adjacents dans l'index), pour limiter
    le nombre de requêtes envoyées à NOMADS (serveur public, sensible au débit)."""
    base = _file_base(run, step)
    entries = _parse_idx(_get(session, base + ".idx").text)

    positions = {}
    for i, e in enumerate(entries):
        if e["param"] in ("UGRD", "VGRD") and e["level"] == "10 m above ground":
            positions[e["param"]] = i
    if "UGRD" not in positions or "VGRD" not in positions:
        raise RuntimeError(f"UGRD/VGRD absents de l'index ({base})")

    i_lo, i_hi = min(positions.values()), max(positions.values())
    off = entries[i_lo]["offset"]
    end = entries[i_hi + 1]["offset"] - 1 if i_hi + 1 < len(entries) else off + 4_000_000
    r = _get(session, base, headers={"Range": f"bytes={off}-{end}"})
    data = r.content

    out = {}
    for param, i in positions.items():
        rel_off = entries[i]["offset"] - off
        rel_end = (entries[i + 1]["offset"] - off) if i + 1 < len(entries) else len(data)
        out[param] = _decode_grib(data[rel_off:rel_end], lat_idx, lon_idx)
    return step, out["UGRD"], out["VGRD"]


def _fetch_step_uv(session, run, step, lat_idx, lon_idx, max_retries=4):
    """NOMADS répond parfois (rarement) de façon tronquée/incomplète : on retente l'échéance entière.
    Un blocage de débit (_NomadsBlocked) n'est en revanche jamais retenté ici — inutile de
    continuer à insister tant qu'on est bloqué, voir _download()."""
    for attempt in range(max_retries):
        try:
            return _fetch_step_uv_once(session, run, step, lat_idx, lon_idx)
        except _NomadsBlocked:
            raise
        except Exception:
            if attempt == max_retries - 1:
                raise
            time.sleep(1 + attempt)


def _download():
    steps = _forecast_steps()
    lat_idx, lon_idx, lats, lons = _grid_bbox_indices()

    with requests.Session() as session:
        run = _pick_run(session)
        u_by_step, v_by_step = {}, {}
        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(_fetch_step_uv, session, run, s, lat_idx, lon_idx): s for s in steps}
            try:
                for fut in as_completed(futures):
                    step, u, v = fut.result()
                    u_by_step[step] = u
                    v_by_step[step] = v
            except _NomadsBlocked:
                pool.shutdown(wait=False, cancel_futures=True)
                raise

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
        "model":       "NOAA GFS 0.25°",
        "days":        int(steps[-1] // 24),
    }

    kw = dict(method="linear", bounds_error=False, fill_value=None)
    interp_u = RegularGridInterpolator((lons, lats, times_h), u_grid, **kw)
    interp_v = RegularGridInterpolator((lons, lats, times_h), v_grid, **kw)
    return meta, interp_u, interp_v


def _get_cached():
    now    = datetime.now(timezone.utc)
    cached = _cache.get("gfs025")
    if cached and (now - cached[3]).total_seconds() < CACHE_AGE:
        return cached[0], cached[1], cached[2]

    fail = _fail_cache.get("gfs025")
    if fail and (now - fail[0]).total_seconds() < fail[2]:
        raise RuntimeError(fail[1])

    try:
        meta, iu, iv = _download()
    except _NomadsBlocked as e:
        _fail_cache["gfs025"] = (now, str(e), BLOCKED_COOLDOWN)
        raise
    except Exception as e:
        _fail_cache["gfs025"] = (now, str(e), FAIL_COOLDOWN)
        raise
    _fail_cache.pop("gfs025", None)
    _cache["gfs025"] = (meta, iu, iv, now)
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
