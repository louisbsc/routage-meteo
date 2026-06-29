"""
Téléchargement planifié des prévisions ECMWF IFS HRES et GFS 0.25°.

Architecture :
- Un worker tourne en arrière-plan (APScheduler, dans le même processus FastAPI).
- Il se réveille aux heures de publication connues et télécharge le dernier run.
- Les données sont sauvegardées sur disque au format .npz (u/v grid + méta).
- L'API lit toujours depuis le disque — elle ne déclenche jamais de téléchargement.
- Le run précédent est conservé comme fallback pendant le téléchargement du suivant.

Remplace vents_ecmwf.py et vents_gfs.py. Expose deux instances :
    from inputs.wind_fetcher import ecmwf, gfs
    ecmwf.get_meta()   ecmwf.get_V_deg()   ecmwf.get_V_nm()
    gfs.get_meta()     gfs.get_V_deg()     gfs.get_V_nm()
"""

import json
import logging
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

import eccodes
import numpy as np
import requests
from scipy.interpolate import RegularGridInterpolator
from shapely import contains_xy

from inputs.vents import land_geom

logger = logging.getLogger(__name__)

CACHE_DIR = Path(__file__).parent.parent / "data" / "cache"
BBOX = (-180.0, -80.0, 179.75, 80.0)  # lon0, lat0, lon1, lat1 (couverture mondiale)
HEADERS = {"User-Agent": "Mozilla/5.0 (compatible; sailing-router/1.0)"}

_RETRY_INTERVAL_S = 120    # secondes entre deux tentatives
_RETRY_MAX = 15            # 15 × 2 min = 30 min max par cycle


# ── Classe de base ──────────────────────────────────────────────────────────

class WindModel:
    """
    Logique commune aux deux modèles : sauvegarde/chargement .npz,
    cache mémoire, orchestration du téléchargement, planning APScheduler,
    interface publique get_meta / get_V_deg / get_V_nm.

    Les sous-classes implémentent :
        _get_run(session)         → datetime du run à télécharger
        _forecast_steps(run)      → liste des échéances (entiers, heures)
        _fetch_step_uv(...)       → (step, u_array, v_array) pour une échéance
        _grid_bbox_indices()      → (lat_idx, lon_idx, lats, lons)  [optionnel]
    """

    def __init__(self, key, model_name, n_workers, schedule_utc_hours):
        self.key = key
        self.model_name = model_name
        self.n_workers = n_workers
        self.schedule_utc_hours = schedule_utc_hours

        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        for f in CACHE_DIR.glob("_tmp_*.npz"):   # fichiers temporaires interrompus
            f.unlink(missing_ok=True)

        self._mem_path = None       # chemin du .npz actuellement en mémoire
        self._mem_data = None       # (meta, interp_u, interp_v)
        self._load_lock = threading.Lock()
        self._dl_lock = threading.Lock()   # empêche deux téléchargements simultanés

    # ── À surcharger ────────────────────────────────────────────────────────

    def _get_run(self, session):
        raise NotImplementedError

    def _forecast_steps(self, run):
        raise NotImplementedError

    def _fetch_step_uv(self, session, run, step, lat_idx, lon_idx):
        raise NotImplementedError

    # ── Grille commune (ECMWF) — GFSModel la surcharge ─────────────────────

    def _grid_bbox_indices(self):
        lon0, lat0, lon1, lat1 = BBOX
        lats_full = 90.0 - np.arange(721) * 0.25
        lons_full = -180.0 + np.arange(1440) * 0.25
        lat_idx = np.where((lats_full >= lat0) & (lats_full <= lat1))[0][::-1]
        lon_idx = np.where((lons_full >= lon0) & (lons_full <= lon1))[0]
        # Sous-échantillonner à 0.5° pour limiter la taille des fichiers .npz mondiaux
        lat_idx = lat_idx[::2]
        lon_idx = lon_idx[::2]
        return lat_idx, lon_idx, lats_full[lat_idx], lons_full[lon_idx]

    # ── Orchestration du téléchargement ─────────────────────────────────────

    def _run_download(self):
        """Télécharge toutes les échéances en parallèle. Retourne les grilles brutes."""
        lat_idx, lon_idx, lats, lons = self._grid_bbox_indices()
        with requests.Session() as session:
            run = self._get_run(session)
            steps = self._forecast_steps(run)
            u_by_step, v_by_step = {}, {}
            with ThreadPoolExecutor(max_workers=self.n_workers) as pool:
                futures = [
                    pool.submit(self._fetch_step_uv, session, run, s, lat_idx, lon_idx)
                    for s in steps
                ]
                for fut in futures:
                    step, u, v = fut.result()
                    u_by_step[step] = u
                    v_by_step[step] = v
        return run, steps, lats, lons, u_by_step, v_by_step

    def _download_and_save(self):
        """Télécharge le run courant et le sauvegarde en .npz. Thread-safe."""
        if not self._dl_lock.acquire(blocking=False):
            logger.info("[%s] téléchargement déjà en cours, ignoré", self.key)
            return

        try:
            run, steps, lats, lons, u_by_step, v_by_step = self._run_download()

            n_lon, n_lat, n_t = len(lons), len(lats), len(steps)
            u_grid = np.empty((n_lon, n_lat, n_t))
            v_grid = np.empty((n_lon, n_lat, n_t))
            for k, step in enumerate(steps):
                u_grid[:, :, k] = u_by_step[step].T
                v_grid[:, :, k] = v_by_step[step].T

            times_h = np.array(steps, dtype=np.float64)
            self._save_npz(CACHE_DIR / self._npz_name(run),
                           u_grid, v_grid, lons, lats, times_h, run, steps)
            self._cleanup_old_files()
        finally:
            self._dl_lock.release()

    def _try_fetch_with_retries(self):
        """Tentative de téléchargement avec retries. Appelé par le scheduler."""
        for attempt in range(_RETRY_MAX):
            try:
                self._download_and_save()
                return
            except Exception as e:
                logger.warning("[%s] tentative %d/%d échouée : %s",
                               self.key, attempt + 1, _RETRY_MAX, e)
                if attempt < _RETRY_MAX - 1:
                    time.sleep(_RETRY_INTERVAL_S)
        logger.error("[%s] toutes les tentatives ont échoué, données existantes conservées",
                     self.key)

    # ── Fichiers .npz ────────────────────────────────────────────────────────

    def _npz_name(self, run):
        return f"{self.key}_{run.strftime('%Y%m%d_%Hz')}.npz"

    def _find_latest_npz(self):
        files = sorted(CACHE_DIR.glob(f"{self.key}_*.npz"), reverse=True)
        return files[0] if files else None

    def _cleanup_old_files(self):
        files = sorted(CACHE_DIR.glob(f"{self.key}_*.npz"), reverse=True)
        for old in files[2:]:
            old.unlink(missing_ok=True)

    def _save_npz(self, path, u_grid, v_grid, lons, lats, times_h, run, steps):
        valid_dts = [run + timedelta(hours=int(s)) for s in steps]
        valid_times = np.array(
            [dt.strftime("%Y-%m-%dT%H:%M:%S+00:00") for dt in valid_dts]
        )
        # Préfixe _tmp_ pour : (1) garder l'extension .npz (numpy n'en rajoute pas),
        # (2) exclure du glob "ecmwf_*.npz" pendant l'écriture.
        tmp = path.parent / ("_tmp_" + path.name)
        np.savez_compressed(
            tmp,
            u_grid=u_grid.astype(np.float32),
            v_grid=v_grid.astype(np.float32),
            lons=lons,
            lats=lats,
            times_h=times_h,
            run_iso=np.array([run.isoformat()]),
            valid_times=valid_times,
            model_name=np.array([self.model_name]),
        )
        tmp.rename(path)   # rename atomique : l'API ne voit jamais un fichier partiel
        logger.info("[%s] run %s sauvegardé → %s", self.key, run.isoformat(), path.name)

    def _load_npz(self, path):
        d = np.load(path, allow_pickle=False)
        lons = d["lons"].astype(np.float64)
        lats = d["lats"].astype(np.float64)
        times_h = d["times_h"].astype(np.float64)
        u_grid = d["u_grid"]  # float32 : précision suffisante, réduit la RAM de 2×
        v_grid = d["v_grid"]

        kw = dict(method="linear", bounds_error=False, fill_value=None)
        interp_u = RegularGridInterpolator((lons, lats, times_h), u_grid, **kw)
        interp_v = RegularGridInterpolator((lons, lats, times_h), v_grid, **kw)

        meta = {
            "valid_times": d["valid_times"].tolist(),
            "times":       times_h.tolist(),
            "bbox":        [float(lons[0]), float(lats[0]), float(lons[-1]), float(lats[-1])],
            "model":       str(d["model_name"][0]),
            "days":        int(times_h[-1] // 24),
            "run_time":    str(d["run_iso"][0]),   # timestamp du run pour affichage
        }
        return meta, interp_u, interp_v

    # ── Cache mémoire ────────────────────────────────────────────────────────

    def _get_current(self):
        """
        Retourne (meta, interp_u, interp_v) depuis le cache mémoire si le .npz
        n'a pas changé, sinon recharge depuis le disque.
        """
        latest = self._find_latest_npz()
        if latest is None:
            raise RuntimeError(
                f"Aucune donnée {self.key.upper()} disponible — téléchargement initial en cours"
            )

        if self._mem_path == latest:
            return self._mem_data   # chemin rapide, pas de lock nécessaire

        with self._load_lock:
            if self._mem_path == latest:   # double-check après acquisition du lock
                return self._mem_data
            meta, iu, iv = self._load_npz(latest)
            self._mem_data = (meta, iu, iv)
            self._mem_path = latest
            return meta, iu, iv

    # ── Démarrage ────────────────────────────────────────────────────────────

    def ensure_fresh(self):
        """
        Appelé au démarrage du serveur.
        Si aucun .npz n'existe, lance un téléchargement immédiat.
        Si le BBOX du fichier existant diffère du BBOX configuré (ex. extension mondiale),
        supprime l'ancien fichier et re-télécharge automatiquement.
        """
        latest = self._find_latest_npz()
        if latest is None:
            logger.info("[%s] aucun cache disque, téléchargement initial...", self.key)
            self._try_fetch_with_retries()
            return
        # Vérifie que le BBOX stocké correspond au BBOX configuré
        try:
            d = np.load(latest, allow_pickle=False)
            _, _, lats_exp, lons_exp = self._grid_bbox_indices()
            if (abs(float(d["lons"][0]) - lons_exp[0]) > 2 or
                    abs(float(d["lons"][-1]) - lons_exp[-1]) > 2 or
                    abs(float(d["lats"][0]) - lats_exp[0]) > 2 or
                    abs(float(d["lats"][-1]) - lats_exp[-1]) > 2):
                logger.info("[%s] BBOX changé, suppression du cache et re-téléchargement…", self.key)
                latest.unlink(missing_ok=True)
                self._try_fetch_with_retries()
                return
        except Exception:
            self._try_fetch_with_retries()
            return
        logger.info("[%s] cache disque trouvé : %s", self.key, latest.name)

    # ── Planning APScheduler ─────────────────────────────────────────────────

    def register_jobs(self, scheduler):
        from apscheduler.triggers.cron import CronTrigger
        hours = ",".join(str(h) for h in self.schedule_utc_hours)
        scheduler.add_job(
            self._try_fetch_with_retries,
            CronTrigger(hour=hours, minute=5, timezone="UTC"),
            id=f"{self.key}_fetch",
            replace_existing=True,
        )
        logger.info("[%s] planifié à %sh05 UTC", self.key, hours)

    # ── Interface publique ───────────────────────────────────────────────────

    def get_meta(self):
        meta, _, _ = self._get_current()
        return meta

    def _uv(self, interp_u, interp_v, lons_arr, lats_arr, t):
        pts = np.column_stack([lons_arr, lats_arr, np.full(len(lons_arr), t)])
        u = np.nan_to_num(interp_u(pts))
        v = np.nan_to_num(interp_v(pts))
        land = contains_xy(land_geom, lons_arr, lats_arr)
        u[land] = v[land] = 0.0
        return u, v

    @staticmethod
    def _to_dir_force(u, v):
        force = np.sqrt(u**2 + v**2) * 1.94384
        direction = (np.degrees(np.arctan2(u, v)) + 180) % 360
        return direction, force

    def get_V_deg(self):
        _, iu, iv = self._get_current()

        def V(p, t):
            p = np.asarray(p, dtype=float)
            batch = p.ndim == 2
            xs = p[:, 0] if batch else p[0:1]
            ys = p[:, 1] if batch else p[1:2]
            u, v = self._uv(iu, iv, xs, ys, t)
            d, f = self._to_dir_force(u, v)
            res = np.column_stack([d, f])
            return res if batch else res[0]

        return V

    def get_V_nm(self):
        _, iu, iv = self._get_current()

        def V(p, t):
            p = np.asarray(p, dtype=float)
            batch = p.ndim == 2
            xs = p[:, 0] if batch else p[0:1]
            ys = p[:, 1] if batch else p[1:2]
            lons_arr = xs / (60.0 * 0.7)
            lats_arr = ys / 60.0
            u, v = self._uv(iu, iv, lons_arr, lats_arr, t)
            d, f = self._to_dir_force(u, v)
            res = np.column_stack([d, f])
            return res if batch else res[0]

        return V


# ── Modèle ECMWF ────────────────────────────────────────────────────────────

class ECMWFModel(WindModel):

    _BASE_URL       = "https://data.ecmwf.int/forecasts"
    _RESOL          = "0p25"
    _STREAM         = "oper"
    _PUBLISH_DELAY_H = 7

    def __init__(self):
        super().__init__(
            key="ecmwf",
            model_name="ECMWF IFS HRES 0.25°",
            n_workers=8,
            schedule_utc_hours=[1, 7, 13, 19],
        )

    def _latest_run(self):
        candidate = datetime.now(timezone.utc) - timedelta(hours=self._PUBLISH_DELAY_H)
        hour = (candidate.hour // 6) * 6
        return candidate.replace(hour=hour, minute=0, second=0, microsecond=0)

    def _get_run(self, session):
        return self._latest_run()

    def _forecast_steps(self, run):
        base = list(range(0, 145, 3))
        if run.hour in (0, 12):
            base += list(range(150, 361, 6))
        return base

    def _url_base(self, run, step):
        d, h = run.strftime("%Y%m%d"), run.strftime("%H")
        return (
            f"{self._BASE_URL}/{d}/{h}z/ifs/{self._RESOL}/{self._STREAM}"
            f"/{d}{h}0000-{step}h-{self._STREAM}-fc"
        )

    @staticmethod
    def _http_get(session, url, headers=None, max_retries=4, timeout=60):
        for attempt in range(max_retries):
            if attempt:
                time.sleep(2 ** attempt)
            r = session.get(url, headers={**HEADERS, **(headers or {})}, timeout=timeout)
            if r.status_code == 429 and attempt < max_retries - 1:
                continue
            r.raise_for_status()
            return r

    def _fetch_step_uv(self, session, run, step, lat_idx, lon_idx):
        base = self._url_base(run, step)
        entries = {}
        for line in self._http_get(session, base + ".index").text.splitlines():
            rec = json.loads(line)
            entries[rec["param"]] = rec

        out = []
        for param in ("10u", "10v"):
            entry = entries.get(param)
            if entry is None:
                raise RuntimeError(f"paramètre {param} absent de l'index ({base})")
            off, ln = entry["_offset"], entry["_length"]
            r = self._http_get(
                session, base + ".grib2",
                headers={"Range": f"bytes={off}-{off + ln - 1}"}
            )
            gid = eccodes.codes_new_from_message(r.content)
            try:
                Ni = eccodes.codes_get(gid, "Ni")
                Nj = eccodes.codes_get(gid, "Nj")
                vals = np.array(eccodes.codes_get_values(gid), dtype=float).reshape(Nj, Ni)
            finally:
                eccodes.codes_release(gid)
            out.append(vals[lat_idx, :][:, lon_idx])
        return step, out[0], out[1]


# ── Modèle GFS ───────────────────────────────────────────────────────────────

class _NomadsBlocked(RuntimeError):
    """NOMADS répond par sa page HTML « Over Rate Limit » (déguisée en 200/302)."""


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


class GFSModel(WindModel):

    _BASE_URL        = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
    _PUBLISH_DELAY_H = 5
    _MAX_RUN_FALLBACK = 3
    _MIN_REQUEST_GAP = 0.4
    _BLOCKED_WAIT_S  = 600   # attente spéciale si NOMADS nous bloque

    def __init__(self):
        super().__init__(
            key="gfs",
            model_name="NOAA GFS 0.25°",
            n_workers=3,
            schedule_utc_hours=[5, 11, 17, 23],
        )
        self._rate_limiter = _RateLimiter(self._MIN_REQUEST_GAP)

    # ── Grille GFS (convention 0..360°) ─────────────────────────────────────

    def _grid_bbox_indices(self):
        lon0, lat0, lon1, lat1 = BBOX
        lats_full = 90.0 - np.arange(721) * 0.25
        lat_idx = np.where((lats_full >= lat0) & (lats_full <= lat1))[0][::-1]
        lon_idx, lons = self._lon_indices_0_360(lon0, lon1)
        # Sous-échantillonner à 0.5° pour limiter la taille des fichiers .npz mondiaux
        lat_idx = lat_idx[::2]
        lon_idx = lon_idx[::2]
        lons    = lons[::2]
        return lat_idx, lon_idx, lats_full[lat_idx], lons

    @staticmethod
    def _lon_indices_0_360(lon0, lon1, n=1440, res=0.25):
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

    # ── Choix du run ────────────────────────────────────────────────────────

    def _run_candidates(self):
        now = datetime.now(timezone.utc) - timedelta(hours=self._PUBLISH_DELAY_H)
        hour = (now.hour // 6) * 6
        latest = now.replace(hour=hour, minute=0, second=0, microsecond=0)
        return [latest - timedelta(hours=6 * k) for k in range(self._MAX_RUN_FALLBACK + 1)]

    def _get_run(self, session):
        for run in self._run_candidates():
            url = self._file_base(run, 0) + ".idx"
            try:
                self._http_get(session, url, max_retries=1)
                return run
            except _NomadsBlocked:
                raise
            except requests.RequestException:
                continue
        raise RuntimeError("Aucun run GFS récent n'est accessible sur NOMADS")

    def _forecast_steps(self, run):
        return list(range(0, 121, 1)) + list(range(123, 385, 3))

    # ── Requêtes HTTP ────────────────────────────────────────────────────────

    def _file_base(self, run, step):
        d, h = run.strftime("%Y%m%d"), run.strftime("%H")
        return f"{self._BASE_URL}/gfs.{d}/{h}/atmos/gfs.t{h}z.pgrb2.0p25.f{step:03d}"

    def _http_get(self, session, url, headers=None, max_retries=4, timeout=60):
        for attempt in range(max_retries):
            if attempt:
                time.sleep(2 ** attempt)
            self._rate_limiter.wait()
            r = session.get(url, headers={**HEADERS, **(headers or {})}, timeout=timeout)
            if r.history or "text/html" in r.headers.get("Content-Type", ""):
                raise _NomadsBlocked(
                    "NOMADS a limité le débit (Over Rate Limit) — réessayer plus tard"
                )
            if r.status_code in (429, 503) and attempt < max_retries - 1:
                continue
            r.raise_for_status()
            return r

    # ── Décodage GRIB ────────────────────────────────────────────────────────

    @staticmethod
    def _parse_idx(text):
        entries = []
        for line in text.splitlines():
            if not line.strip():
                continue
            parts = line.split(":")
            entries.append({"offset": int(parts[1]), "param": parts[3], "level": parts[4]})
        return entries

    @staticmethod
    def _decode_grib(message, lat_idx, lon_idx):
        gid = eccodes.codes_new_from_message(message)
        try:
            Ni = eccodes.codes_get(gid, "Ni")
            Nj = eccodes.codes_get(gid, "Nj")
            vals = np.array(eccodes.codes_get_values(gid), dtype=float).reshape(Nj, Ni)
        finally:
            eccodes.codes_release(gid)
        return vals[lat_idx, :][:, lon_idx]

    def _fetch_step_uv_once(self, session, run, step, lat_idx, lon_idx):
        base = self._file_base(run, step)
        entries = self._parse_idx(self._http_get(session, base + ".idx").text)

        positions = {}
        for i, e in enumerate(entries):
            if e["param"] in ("UGRD", "VGRD") and e["level"] == "10 m above ground":
                positions[e["param"]] = i
        if "UGRD" not in positions or "VGRD" not in positions:
            raise RuntimeError(f"UGRD/VGRD absents de l'index ({base})")

        i_lo, i_hi = min(positions.values()), max(positions.values())
        off = entries[i_lo]["offset"]
        end = entries[i_hi + 1]["offset"] - 1 if i_hi + 1 < len(entries) else off + 4_000_000
        r = self._http_get(session, base, headers={"Range": f"bytes={off}-{end}"})
        data = r.content

        out = {}
        for param, i in positions.items():
            rel_off = entries[i]["offset"] - off
            rel_end = (entries[i + 1]["offset"] - off) if i + 1 < len(entries) else len(data)
            out[param] = self._decode_grib(data[rel_off:rel_end], lat_idx, lon_idx)
        return step, out["UGRD"], out["VGRD"]

    def _fetch_step_uv(self, session, run, step, lat_idx, lon_idx, max_retries=4):
        for attempt in range(max_retries):
            try:
                return self._fetch_step_uv_once(session, run, step, lat_idx, lon_idx)
            except _NomadsBlocked:
                raise
            except Exception:
                if attempt == max_retries - 1:
                    raise
                time.sleep(1 + attempt)

    # ── Orchestration spécifique GFS (annulation si bloqué) ─────────────────

    def _run_download(self):
        lat_idx, lon_idx, lats, lons = self._grid_bbox_indices()
        with requests.Session() as session:
            run = self._get_run(session)
            steps = self._forecast_steps(run)
            u_by_step, v_by_step = {}, {}
            with ThreadPoolExecutor(max_workers=self.n_workers) as pool:
                futures = {
                    pool.submit(self._fetch_step_uv, session, run, s, lat_idx, lon_idx): s
                    for s in steps
                }
                try:
                    for fut in as_completed(futures):
                        step, u, v = fut.result()
                        u_by_step[step] = u
                        v_by_step[step] = v
                except _NomadsBlocked:
                    pool.shutdown(wait=False, cancel_futures=True)
                    raise
        return run, steps, lats, lons, u_by_step, v_by_step

    def _try_fetch_with_retries(self):
        for attempt in range(_RETRY_MAX):
            try:
                self._download_and_save()
                return
            except _NomadsBlocked as e:
                logger.warning("[gfs] NOMADS bloqué : %s — attente %ds", e, self._BLOCKED_WAIT_S)
                time.sleep(self._BLOCKED_WAIT_S)
            except Exception as e:
                logger.warning("[gfs] tentative %d/%d échouée : %s", attempt + 1, _RETRY_MAX, e)
                if attempt < _RETRY_MAX - 1:
                    time.sleep(_RETRY_INTERVAL_S)
        logger.error("[gfs] toutes les tentatives ont échoué, données existantes conservées")


# ── Instances exposées ───────────────────────────────────────────────────────

ecmwf = ECMWFModel()
gfs   = GFSModel()
