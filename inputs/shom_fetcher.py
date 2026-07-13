"""
Téléchargement des prévisions de courant SHOM via le service ncWMS2 public.

Source   : https://services.data.shom.fr/ncwms2/wms
Données  : HYCOM de surface — HYDRODYN-SURF_HYCOM3D-SURF_R1000_MANGASC
Zone     : Manche + Golfe de Gascogne  (-15 °W → 3 °E  /  43 °N → 51 °N)
Résolution temporelle : 1 h  (5 jours de prévision)
Résolution spatiale   : ~1 km (grille 180×80 px ≈ 0,1 °)
Accès : public, sans clé API.

Interface publique (identique à WindModel dans wind_fetcher.py) :
    from inputs.shom_fetcher import shom_mangasc
    shom_mangasc.get_meta()    → dict (times, valid_times, bbox, model, days)
    shom_mangasc.get_V_deg()   → C(pts_deg, t_h)   courant en nœuds
    shom_mangasc.get_V_nm()    → C(pts_nm,  t_h)   courant en nœuds (coords NM)
"""

import io
import logging
import re
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import requests
from PIL import Image
from scipy.interpolate import RegularGridInterpolator
from shapely import contains_xy

from inputs.vents import land_geom

logger = logging.getLogger(__name__)

NCWMS_URL = "https://services.data.shom.fr/ncwms2/wms"
CACHE_DIR  = Path(__file__).parent.parent / "data" / "cache"
HEADERS    = {"User-Agent": "Mozilla/5.0 (compatible; sailing-router/1.0)"}
_EPOCH     = datetime(1970, 1, 1, tzinfo=timezone.utc)


def _parse_iso_duration(s: str) -> timedelta:
    """Parse durée ISO 8601 (ex. PT1H, PT6H, P1D)."""
    s = s.lstrip("P")
    date_part, time_part = (s.split("T") + [""])[:2]
    result = timedelta()
    for m in re.finditer(r"(\d+)D", date_part):
        result += timedelta(days=int(m.group(1)))
    for m in re.finditer(r"(\d+)H", time_part):
        result += timedelta(hours=int(m.group(1)))
    for m in re.finditer(r"(\d+)M", time_part):
        result += timedelta(minutes=int(m.group(1)))
    return result


def _strip_ns(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


class SHOMCurrentModel:
    """
    Courant de surface SHOM HYCOM — zone Manche/Atlantique.

    Télécharge u, v depuis ncWMS2 (PNG mode=32bit) pour les prochaines MAX_HOURS.
    Cache les données dans data/cache/shom_mangasc_*.npz.
    Rafraîchissement planifié quotidien à 12h30 UTC (SHOM publie vers 10h UTC).
    """

    LAYER_PREFIX = "HYDRODYN-SURF_HYCOM3D-SURF_R1000_MANGASC"
    BBOX         = (-15.0, 43.01, 2.97, 51.0)   # lon0, lat0, lon1, lat1
    GRID_W       = 180    # pixels longitude  (≈ 0,1 ° / px)
    GRID_H       = 80     # pixels latitude   (≈ 0,1 ° / px)
    SCALE_MIN    = -3.0   # m/s  (composantes u et v)
    SCALE_MAX    =  3.0
    MAX_HOURS    = 120    # horizon de téléchargement
    N_WORKERS    = 4      # requêtes parallèles (ne pas surcharger le serveur SHOM)

    def __init__(self):
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        for f in CACHE_DIR.glob("_tmp_shom_*.npz"):
            f.unlink(missing_ok=True)
        self._mem_path  = None
        self._mem_data  = None
        self._load_lock = threading.Lock()
        self._dl_lock   = threading.Lock()

    # ── Découverte des pas de temps disponibles ──────────────────────────────

    def _get_available_times(self, session: requests.Session) -> list:
        """Interroge GetCapabilities et retourne la liste des datetime disponibles."""
        r = session.get(NCWMS_URL, headers=HEADERS, params={
            "SERVICE": "WMS", "VERSION": "1.3.0", "REQUEST": "GetCapabilities",
        }, timeout=30)
        r.raise_for_status()

        root = ET.fromstring(r.content)
        time_str = None

        # Cherche la Dimension time de notre couche (ou la première PT1H)
        target = f"{self.LAYER_PREFIX}/u"
        found_target = False
        for el in root.iter():
            tag = _strip_ns(el.tag)
            if tag == "Name" and el.text and el.text.strip() == target:
                found_target = True
            if found_target and tag == "Dimension" and el.get("name", "").lower() == "time":
                time_str = (el.text or "").strip()
                break

        if not time_str:
            for el in root.iter():
                tag = _strip_ns(el.tag)
                if (tag == "Dimension"
                        and el.get("name", "").lower() == "time"
                        and el.text and "PT1H" in el.text):
                    time_str = el.text.strip()
                    break

        if not time_str:
            raise RuntimeError("Impossible de trouver la dimension temporelle SHOM ncWMS2")

        parts   = time_str.split("/")
        t_start = datetime.fromisoformat(parts[0].replace(".000Z", "+00:00"))
        t_end   = datetime.fromisoformat(parts[1].replace(".000Z", "+00:00"))
        step    = _parse_iso_duration(parts[2])

        now     = datetime.now(timezone.utc).replace(minute=0, second=0, microsecond=0)
        t_first = max(t_start, now)
        t_first = min(t_first, t_end - step)

        times, t = [], t_first
        while t <= t_end and len(times) < self.MAX_HOURS:
            times.append(t)
            t += step
        return times

    # ── Téléchargement d'un tableau de valeurs ───────────────────────────────

    def _get_map(self, session: requests.Session, var: str, time_dt: datetime) -> np.ndarray:
        """Retourne la grille numpy (H, W) pour la composante var au temps time_dt."""
        lon0, lat0, lon1, lat1 = self.BBOX
        time_iso = time_dt.strftime("%Y-%m-%dT%H:%M:%S.000Z")
        params = {
            "SERVICE": "WMS", "VERSION": "1.3.0", "REQUEST": "GetMap",
            "LAYERS":  f"{self.LAYER_PREFIX}/{var}",
            "STYLES":  "",
            "CRS":     "EPSG:4326",
            # WMS 1.3.0 + EPSG:4326 → ordre des axes : latitude, longitude
            "BBOX":    f"{lat0},{lon0},{lat1},{lon1}",
            "WIDTH":   self.GRID_W,
            "HEIGHT":  self.GRID_H,
            "TIME":    time_iso,
            "FORMAT":  "image/png;mode=32bit",
            "COLORSCALERANGE": f"{self.SCALE_MIN},{self.SCALE_MAX}",
            "TRANSPARENT": "true",
        }
        for attempt in range(3):
            try:
                r = session.get(NCWMS_URL, headers=HEADERS, params=params, timeout=60)
                r.raise_for_status()
                break
            except Exception as e:
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)

        img = Image.open(io.BytesIO(r.content)).convert("RGBA")
        arr = np.array(img, dtype=np.float64)           # (H, W, 4)
        R, G, B, A = arr[..., 0], arr[..., 1], arr[..., 2], arr[..., 3]
        # Décodage mode=32bit : value = min + (R*65536 + G*256 + B) / 16777215 * (max - min)
        values = (self.SCALE_MIN
                  + (R * 65536 + G * 256 + B) / 16777215.0
                  * (self.SCALE_MAX - self.SCALE_MIN))
        values[A == 0] = 0.0    # terre / hors domaine → courant nul
        return values            # (H, W), ligne 0 = nord

    def _fetch_step(self, session: requests.Session, time_dt: datetime):
        u = self._get_map(session, "u", time_dt)
        v = self._get_map(session, "v", time_dt)
        return time_dt, u, v

    # ── Téléchargement complet et sauvegarde .npz ────────────────────────────

    def _npz_name(self, t0: datetime) -> str:
        return f"shom_mangasc_{t0.strftime('%Y%m%d_%Hz')}.npz"

    def _find_latest_npz(self):
        files = sorted(CACHE_DIR.glob("shom_mangasc_*.npz"), reverse=True)
        return files[0] if files else None

    def _cleanup_old(self):
        for f in sorted(CACHE_DIR.glob("shom_mangasc_*.npz"), reverse=True)[2:]:
            f.unlink(missing_ok=True)

    def _download_and_save(self):
        if not self._dl_lock.acquire(blocking=False):
            logger.info("[shom] téléchargement déjà en cours, ignoré")
            return
        try:
            with requests.Session() as session:
                times = self._get_available_times(session)
            if not times:
                raise RuntimeError("Aucun pas de temps SHOM disponible")
            logger.info("[shom] %d pas de temps : %s → %s",
                        len(times), times[0].isoformat()[:16], times[-1].isoformat()[:16])

            results: dict = {}
            with requests.Session() as session:
                with ThreadPoolExecutor(max_workers=self.N_WORKERS) as pool:
                    futures = {pool.submit(self._fetch_step, session, t): t for t in times}
                    for fut in as_completed(futures):
                        t_dt, u_img, v_img = fut.result()
                        results[t_dt] = (u_img, v_img)
                        logger.debug("[shom] ✓ %s", t_dt.isoformat()[:16])

            lon0, lat0, lon1, lat1 = self.BBOX
            lons   = np.linspace(lon0, lon1, self.GRID_W)   # croissant ouest→est
            lats   = np.linspace(lat0, lat1, self.GRID_H)   # croissant sud→nord
            n_t    = len(times)
            u_grid = np.zeros((self.GRID_W, self.GRID_H, n_t), dtype=np.float32)
            v_grid = np.zeros((self.GRID_W, self.GRID_H, n_t), dtype=np.float32)

            for k, t in enumerate(times):
                u_img, v_img = results[t]
                # L'image WMS a ligne 0 = nord → retourner verticalement puis transposer
                u_grid[:, :, k] = u_img[::-1, :].T.astype(np.float32)
                v_grid[:, :, k] = v_img[::-1, :].T.astype(np.float32)

            t0_dt     = times[0]
            times_h   = np.arange(n_t, dtype=np.float64)   # offset horaire depuis t0
            valid_times = np.array([t.strftime("%Y-%m-%dT%H:%M:%S+00:00") for t in times])

            path = CACHE_DIR / self._npz_name(t0_dt)
            tmp  = path.parent / ("_tmp_" + path.name)
            np.savez_compressed(
                tmp,
                u_grid=u_grid, v_grid=v_grid,
                lons=lons, lats=lats,
                times_h=times_h,
                run_iso=np.array([t0_dt.isoformat()]),
                valid_times=valid_times,
                model_name=np.array(["SHOM HYCOM Manche/Atlantique"]),
            )
            tmp.rename(path)
            logger.info("[shom] sauvegardé → %s", path.name)
            self._cleanup_old()
        finally:
            self._dl_lock.release()

    def _try_fetch_with_retries(self):
        for attempt in range(5):
            try:
                self._download_and_save()
                return
            except Exception as e:
                logger.warning("[shom] tentative %d/5 échouée : %s", attempt + 1, e)
                if attempt < 4:
                    time.sleep(60 * (attempt + 1))
        logger.error("[shom] toutes les tentatives ont échoué")

    # ── Cache mémoire ────────────────────────────────────────────────────────

    def _load_npz(self, path: Path):
        d        = np.load(path, allow_pickle=False)
        lons     = d["lons"].astype(np.float64)
        lats     = d["lats"].astype(np.float64)
        times_h  = d["times_h"].astype(np.float64)
        u_grid   = d["u_grid"].astype(np.float32)
        v_grid   = d["v_grid"].astype(np.float32)
        kw = dict(method="linear", bounds_error=False, fill_value=0.0)
        interp_u = RegularGridInterpolator((lons, lats, times_h), u_grid, **kw)
        interp_v = RegularGridInterpolator((lons, lats, times_h), v_grid, **kw)
        meta = {
            "valid_times": d["valid_times"].tolist(),
            "times":       times_h.tolist(),
            "bbox":        [float(lons[0]), float(lats[0]), float(lons[-1]), float(lats[-1])],
            "model":       str(d["model_name"][0]),
            "days":        int(times_h[-1] // 24),
            "run_time":    str(d["run_iso"][0]),
        }
        return meta, interp_u, interp_v

    def _get_current(self):
        latest = self._find_latest_npz()
        if latest is None:
            raise RuntimeError("Aucune donnée SHOM disponible — téléchargement en cours")
        if self._mem_path == latest:
            return self._mem_data
        with self._load_lock:
            if self._mem_path == latest:
                return self._mem_data
            data = self._load_npz(latest)
            self._mem_data = data
            self._mem_path = latest
            return data

    # ── Démarrage et planification ───────────────────────────────────────────

    def ensure_fresh(self):
        if self._find_latest_npz() is None:
            logger.info("[shom] aucun cache disque, téléchargement initial…")
            self._try_fetch_with_retries()
        else:
            logger.info("[shom] cache existant trouvé : %s", self._find_latest_npz().name)

    def register_jobs(self, scheduler):
        from apscheduler.triggers.cron import CronTrigger
        scheduler.add_job(
            self._try_fetch_with_retries,
            CronTrigger(hour="12", minute="30", timezone="UTC"),
            id="shom_mangasc_fetch",
            replace_existing=True,
        )
        logger.info("[shom] planifié à 12h30 UTC")

    # ── Interface publique ───────────────────────────────────────────────────

    def get_meta(self) -> dict:
        meta, _, _ = self._get_current()
        return meta

    def _uv(self, interp_u, interp_v, lons_arr, lats_arr, t):
        pts  = np.column_stack([lons_arr, lats_arr, np.full(len(lons_arr), t)])
        u    = np.nan_to_num(interp_u(pts)).astype(float)
        v    = np.nan_to_num(interp_v(pts)).astype(float)
        land = contains_xy(land_geom, lons_arr, lats_arr)
        u[land] = v[land] = 0.0
        return u, v

    @staticmethod
    def _to_dir_force(u, v):
        force     = np.sqrt(u ** 2 + v ** 2) * 1.94384   # m/s → nœuds
        direction = (np.degrees(np.arctan2(u, v)) + 180) % 360
        return direction, force

    def get_V_deg(self):
        """Renvoie C(pts, t) — pts en degrés (lon, lat), t en heures depuis valid_times[0]."""
        _, iu, iv = self._get_current()

        def C(p, t):
            p     = np.asarray(p, dtype=float)
            batch = p.ndim == 2
            xs    = p[:, 0] if batch else p[0:1]
            ys    = p[:, 1] if batch else p[1:2]
            u, v  = self._uv(iu, iv, xs, ys, t)
            d, f  = self._to_dir_force(u, v)
            res   = np.column_stack([d, f])
            return res if batch else res[0]

        return C

    def get_V_nm(self):
        """Renvoie C(pts, t) — pts en milles nautiques (x, y), t en heures depuis valid_times[0]."""
        _, iu, iv = self._get_current()

        def C(p, t):
            p        = np.asarray(p, dtype=float)
            batch    = p.ndim == 2
            xs       = p[:, 0] if batch else p[0:1]
            ys       = p[:, 1] if batch else p[1:2]
            lons_arr = xs / (60.0 * 0.7)
            lats_arr = ys / 60.0
            u, v     = self._uv(iu, iv, lons_arr, lats_arr, t)
            d, f     = self._to_dir_force(u, v)
            res      = np.column_stack([d, f])
            return res if batch else res[0]

        return C


# Instance exposée
shom_mangasc = SHOMCurrentModel()
