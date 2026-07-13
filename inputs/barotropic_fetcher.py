"""
Téléchargement des prévisions de courant Barotropic (service commercial).

Source   : https://barotropic.fr (compte payant requis)
Données  : modèle hydrodynamique haute résolution, composantes u/v du courant de surface
Accès    : identifiants dans .env (BAROTROPIC_USER / BAROTROPIC_PASSWORD),
           login WooCommerce puis récupération du lien de téléchargement du produit
           depuis la page "mon-compte/téléchargements" (le lien exact change à
           chaque renouvellement de commande, on le relit donc à chaque fetch).

Note technique : le lien de téléchargement redirige vers barotropic.org, dont le
certificat TLS ne couvre que barotropic.fr (même serveur, même IP). On réécrit
donc l'hôte de la redirection vers barotropic.fr pour conserver une vérification
TLS valide plutôt que de la désactiver.

Interface publique (identique à SHOMCurrentModel dans shom_fetcher.py) :
    from inputs.barotropic_fetcher import barotropic_manche_120h
    barotropic_manche_120h.get_meta()    → dict (times, valid_times, bbox, model, days)
    barotropic_manche_120h.get_V_deg()   → C(pts_deg, t_h)   courant en nœuds
    barotropic_manche_120h.get_V_nm()    → C(pts_nm,  t_h)   courant en nœuds (coords NM)
"""

from __future__ import annotations

import bz2
import logging
import os
import re
import threading
import time
import unicodedata
from datetime import datetime, timedelta, timezone
from pathlib import Path

import eccodes
import numpy as np
import requests
from dotenv import load_dotenv
from scipy.interpolate import RegularGridInterpolator
from shapely import contains_xy

from inputs.vents import land_geom

load_dotenv()

logger = logging.getLogger(__name__)

ACCOUNT_URL   = "https://barotropic.fr/index.php/mon-compte/"
DOWNLOADS_URL = "https://barotropic.fr/index.php/mon-compte/downloads/"
CACHE_DIR     = Path(__file__).parent.parent / "data" / "cache"
HEADERS       = {"User-Agent": "Mozilla/5.0 (compatible; sailing-router/1.0)"}


class BarotropicCurrentModel:
    """
    Courant de surface Barotropic — un produit (couverture géographique) donné.

    Se connecte au compte Barotropic, récupère le lien de téléchargement courant
    du produit PRODUCT_NAME, télécharge le GRIB2 compressé (.bz2), le décode et
    construit un interpolateur (lon, lat, t_h). Sous-échantillonne la grille
    (native ~700 m) d'un facteur STRIDE pour rester utilisable en mémoire.
    Cache le résultat dans data/cache/barotropic_<slug>_*.npz.
    Rafraîchissement planifié quotidien (Barotropic publie vers 10h UTC).
    """

    PRODUCT_NAME = "Manche Atlantique 120h"
    STRIDE       = 8      # sous-échantillonnage de la grille native (~0.05-0.07°)
    MAX_HOURS    = 130    # marge au-dessus de l'horizon max des produits (120h)

    def __init__(self, product_name: str | None = None, stride: int | None = None):
        if product_name:
            self.PRODUCT_NAME = product_name
        if stride:
            self.STRIDE = stride
        ascii_name = unicodedata.normalize("NFKD", self.PRODUCT_NAME).encode("ascii", "ignore").decode()
        self._slug = re.sub(r"[^a-z0-9]+", "_", ascii_name.lower()).strip("_")

        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        for f in CACHE_DIR.glob(f"_tmp_barotropic_{self._slug}_*"):
            f.unlink(missing_ok=True)
        self._mem_path  = None
        self._mem_data  = None
        self._load_lock = threading.Lock()
        self._dl_lock   = threading.Lock()

    # ── Authentification ──────────────────────────────────────────────────────

    def _login(self, session: requests.Session):
        user = os.environ.get("BAROTROPIC_USER")
        pwd  = os.environ.get("BAROTROPIC_PASSWORD")
        if not user or not pwd:
            raise RuntimeError(
                "BAROTROPIC_USER / BAROTROPIC_PASSWORD manquants (fichier .env)"
            )
        r = session.get(ACCOUNT_URL, headers=HEADERS, timeout=30)
        r.raise_for_status()
        nonce_m   = re.search(r'name="woocommerce-login-nonce" value="([^"]+)"', r.text)
        referer_m = re.search(r'name="_wp_http_referer" value="([^"]+)"', r.text)
        if not nonce_m:
            raise RuntimeError("Formulaire de connexion Barotropic introuvable")

        data = {
            "username": user,
            "password": pwd,
            "woocommerce-login-nonce": nonce_m.group(1),
            "_wp_http_referer": referer_m.group(1) if referer_m else "/index.php/mon-compte/",
            "login": "Se connecter",
        }
        r2 = session.post(ACCOUNT_URL, headers=HEADERS, data=data, timeout=30)
        r2.raise_for_status()
        if "woocommerce-form-login" in r2.text:
            raise RuntimeError("Échec de connexion Barotropic (identifiants invalides ?)")

    # ── Découverte du lien de téléchargement courant ─────────────────────────

    def _get_download_link(self, session: requests.Session) -> str:
        r = session.get(DOWNLOADS_URL, headers=HEADERS, timeout=30)
        r.raise_for_status()
        rows = re.findall(r"<tr>\s*<td class=\"download-product\".*?</tr>", r.text, re.S)
        for row in rows:
            name_m = re.search(r'download-product"[^>]*>\s*<a[^>]*>([^<]+)</a>', row)
            link_m = re.search(r'href="([^"]*download_file=[^"]*)"', row)
            if name_m and link_m and name_m.group(1).strip() == self.PRODUCT_NAME:
                return link_m.group(1).replace("&#038;", "&")
        raise RuntimeError(
            f"Produit Barotropic « {self.PRODUCT_NAME} » introuvable dans mon-compte/downloads/ "
            "(abonnement absent ou expiré ?)"
        )

    # ── Téléchargement + décompression du GRIB ───────────────────────────────

    def _download_grib_bytes(self, session: requests.Session) -> bytes:
        link = self._get_download_link(session)
        r = session.get(link, headers=HEADERS, timeout=30, allow_redirects=False)
        if r.status_code in (301, 302, 303, 307, 308):
            location = r.headers["Location"].replace("https://barotropic.org", "https://barotropic.fr")
            r = session.get(location, headers=HEADERS, timeout=180)
        r.raise_for_status()
        if r.headers.get("content-type", "").startswith("text/html"):
            raise RuntimeError("Réponse HTML inattendue (session expirée ou lien invalide)")
        return bz2.decompress(r.content)

    # ── Parsing GRIB2 → grilles u/v ──────────────────────────────────────────

    def _parse_grib(self, raw: bytes):
        """
        Chaque message du GRIB Barotropic représente un instant unique via
        (dataDate, dataTime) avec forecastTime=0 (pas de vrai "step" de prévision).
        parameterNumber : 2 = u (courant Est), 3 = v (courant Nord).
        """
        tmp_path = CACHE_DIR / f"_tmp_barotropic_{self._slug}_{os.getpid()}.grb2"
        tmp_path.write_bytes(raw)
        try:
            by_time: dict = {}   # datetime -> {"u": arr, "v": arr}
            grid_meta = None
            with open(tmp_path, "rb") as f:
                while True:
                    msg = eccodes.codes_grib_new_from_file(f)
                    if msg is None:
                        break
                    try:
                        param = eccodes.codes_get(msg, "parameterNumber")
                        if param not in (2, 3):
                            continue
                        # validityDate/Time (plutôt que dataDate/dataTime) : certains produits
                        # Barotropic encodent l'heure via un vrai "step" de prévision (Hycom),
                        # d'autres via dataTime avec step=0 (Manche/Atlantique, Finistère HR).
                        vdate = eccodes.codes_get(msg, "validityDate")
                        vtime = eccodes.codes_get(msg, "validityTime")
                        dt = datetime.strptime(f"{vdate:08d}{vtime:04d}", "%Y%m%d%H%M").replace(
                            tzinfo=timezone.utc
                        )
                        if grid_meta is None:
                            ni  = eccodes.codes_get(msg, "Ni")
                            nj  = eccodes.codes_get(msg, "Nj")
                            lat0 = eccodes.codes_get(msg, "latitudeOfFirstGridPointInDegrees")
                            lon0 = eccodes.codes_get(msg, "longitudeOfFirstGridPointInDegrees")
                            dlat = eccodes.codes_get(msg, "jDirectionIncrementInDegrees")
                            dlon = eccodes.codes_get(msg, "iDirectionIncrementInDegrees")
                            grid_meta = (ni, nj, lat0, lon0, dlat, dlon)

                        values = eccodes.codes_get_values(msg)
                        values = np.where(values == 9999, 0.0, values).astype(np.float32)
                        by_time.setdefault(dt, {})["u" if param == 2 else "v"] = values
                    finally:
                        eccodes.codes_release(msg)

            if grid_meta is None or not by_time:
                raise RuntimeError("Aucune donnée de courant (paramètres 2/3) dans le GRIB Barotropic")

            ni, nj, lat0, lon0, dlat, dlon = grid_meta
            stride = self.STRIDE
            lons_full = lon0 + np.arange(ni) * dlon
            lons_full = np.where(lons_full > 180, lons_full - 360, lons_full)
            lats_full = lat0 + np.arange(nj) * dlat
            lons = lons_full[::stride]
            lats = lats_full[::stride]

            times = sorted(t for t, d in by_time.items() if "u" in d and "v" in d)
            t0 = times[0]
            # Plafonne par durée écoulée réelle (et non par nb de messages : certains
            # produits ont un pas de 15 min, d'autres 1h).
            times = [t for t in times if (t - t0).total_seconds() / 3600.0 <= self.MAX_HOURS]

            u_grid = np.zeros((len(lons), len(lats), len(times)), dtype=np.float32)
            v_grid = np.zeros((len(lons), len(lats), len(times)), dtype=np.float32)
            for k, t in enumerate(times):
                u2d = by_time[t]["u"].reshape(nj, ni)[::stride, ::stride]
                v2d = by_time[t]["v"].reshape(nj, ni)[::stride, ::stride]
                u_grid[:, :, k] = u2d.T
                v_grid[:, :, k] = v2d.T

            times_h = np.array([(t - t0).total_seconds() / 3600.0 for t in times])
            valid_times = np.array([t.strftime("%Y-%m-%dT%H:%M:%S+00:00") for t in times])
            return lons, lats, times_h, u_grid, v_grid, t0, valid_times
        finally:
            tmp_path.unlink(missing_ok=True)

    # ── Téléchargement complet et sauvegarde .npz ────────────────────────────

    def _npz_name(self, t0: datetime) -> str:
        return f"barotropic_{self._slug}_{t0.strftime('%Y%m%d_%Hz')}.npz"

    def _find_latest_npz(self):
        files = sorted(CACHE_DIR.glob(f"barotropic_{self._slug}_*.npz"), reverse=True)
        return files[0] if files else None

    def _cleanup_old(self):
        for f in sorted(CACHE_DIR.glob(f"barotropic_{self._slug}_*.npz"), reverse=True)[2:]:
            f.unlink(missing_ok=True)

    def _download_and_save(self):
        if not self._dl_lock.acquire(blocking=False):
            logger.info("[barotropic:%s] téléchargement déjà en cours, ignoré", self._slug)
            return
        try:
            with requests.Session() as session:
                self._login(session)
                raw = self._download_grib_bytes(session)
            logger.info("[barotropic:%s] GRIB téléchargé (%d octets décompressés)", self._slug, len(raw))

            lons, lats, times_h, u_grid, v_grid, t0, valid_times = self._parse_grib(raw)
            logger.info(
                "[barotropic:%s] %d pas de temps : %s → %s",
                self._slug, len(times_h), valid_times[0], valid_times[-1],
            )

            path = CACHE_DIR / self._npz_name(t0)
            tmp  = path.parent / ("_tmp_" + path.name)
            np.savez_compressed(
                tmp,
                u_grid=u_grid, v_grid=v_grid,
                lons=lons, lats=lats,
                times_h=times_h,
                run_iso=np.array([t0.isoformat()]),
                valid_times=valid_times,
                model_name=np.array([f"Barotropic {self.PRODUCT_NAME}"]),
            )
            tmp.rename(path)
            logger.info("[barotropic:%s] sauvegardé → %s", self._slug, path.name)
            self._cleanup_old()
        finally:
            self._dl_lock.release()

    def _try_fetch_with_retries(self):
        for attempt in range(5):
            try:
                self._download_and_save()
                return
            except Exception as e:
                logger.warning("[barotropic:%s] tentative %d/5 échouée : %s", self._slug, attempt + 1, e)
                if attempt < 4:
                    time.sleep(60 * (attempt + 1))
        logger.error("[barotropic:%s] toutes les tentatives ont échoué", self._slug)

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
            raise RuntimeError(f"Aucune donnée Barotropic disponible ({self._slug}) — téléchargement en cours")
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
            logger.info("[barotropic:%s] aucun cache disque, téléchargement initial…", self._slug)
            self._try_fetch_with_retries()
        else:
            logger.info("[barotropic:%s] cache existant trouvé : %s", self._slug, self._find_latest_npz().name)

    def register_jobs(self, scheduler, minute_offset: int = 0):
        """minute_offset : décale l'heure de rafraîchissement (les produits partagent
        un même compte Barotropic, on évite des logins simultanés à 11h00 UTC pile)."""
        from apscheduler.triggers.cron import CronTrigger
        hour, minute = divmod(11 * 60 + minute_offset, 60)
        scheduler.add_job(
            self._try_fetch_with_retries,
            CronTrigger(hour=str(hour), minute=str(minute), timezone="UTC"),
            id=f"barotropic_{self._slug}_fetch",
            replace_existing=True,
        )
        logger.info("[barotropic:%s] planifié à %02dh%02d UTC", self._slug, hour, minute)

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


# Instances exposées — une par produit disponible sur le compte Barotropic.
# STRIDE ajusté par produit pour viser une résolution effective ~0.03-0.08°
# tout en gardant une empreinte mémoire raisonnable (grille native jusqu'à ~1400×1400
# et jusqu'à ~290 pas de temps pour le Finistère haute résolution en 15 min).
barotropic_manche_120h      = BarotropicCurrentModel("Manche Atlantique 120h", stride=8)
barotropic_manche_72h       = BarotropicCurrentModel("Manche Atlantique 72h", stride=8)
barotropic_manche_48h       = BarotropicCurrentModel("Manche Atlantique 48h", stride=8)
barotropic_atlantique_ne    = BarotropicCurrentModel("Atlantique Nord Est 120h", stride=6)
barotropic_finistere_hr     = BarotropicCurrentModel("Finistère haute résolution 72h", stride=10)
barotropic_hycom_manche     = BarotropicCurrentModel("Hycom Manche Atlantique 120h", stride=3)

BAROTROPIC_MODELS = {
    "barotropic_manche_120h":   barotropic_manche_120h,
    "barotropic_manche_72h":    barotropic_manche_72h,
    "barotropic_manche_48h":    barotropic_manche_48h,
    "barotropic_atlantique_ne": barotropic_atlantique_ne,
    "barotropic_finistere_hr":  barotropic_finistere_hr,
    "barotropic_hycom_manche":  barotropic_hycom_manche,
}
