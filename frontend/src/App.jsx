import { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import WindMap from './WindMap.jsx';

const API = 'http://localhost:8000';
const INIT_VIEW = { longitude: -5, latitude: 47, zoom: 4, pitch: 0, bearing: 0 };

const DAYS   = ['Dim', 'Lun', 'Mar', 'Mer', 'Jeu', 'Ven', 'Sam'];
const MONTHS = ['jan', 'fév', 'mar', 'avr', 'mai', 'juin', 'juil', 'aoû', 'sep', 'oct', 'nov', 'déc'];

function fmtDatetime(iso) {
  const d = new Date(iso);
  return `${DAYS[d.getUTCDay()]} ${String(d.getUTCDate()).padStart(2, '0')} ${MONTHS[d.getUTCMonth()]}  ${String(d.getUTCHours()).padStart(2, '0')}h UTC`;
}

function fmtShort(ts_h) {
  const d = new Date(ts_h * 3600000);
  return `${String(d.getUTCDate()).padStart(2,'0')}/${String(d.getUTCMonth()+1).padStart(2,'0')} ${String(d.getUTCHours()).padStart(2,'0')}h`;
}

function fmtCoord(pt) {
  if (!pt) return null;
  const fmt = (dd, pos, neg) => {
    const a = Math.abs(dd);
    const deg = Math.floor(a);
    const min = ((a % 1) * 60).toFixed(2);
    return `${deg}°${min}'${dd >= 0 ? pos : neg}`;
  };
  return `${fmt(pt[0], 'N', 'S')} ${fmt(pt[1], 'E', 'W')}`;
}


const FONT = "'SF Mono', 'Consolas', monospace";

// Modèles de courant à téléchargement direct (SHOM + produits Barotropic).
const CURRENT_MODELS = [
  { key: 'shom',        label: 'SHOM HYCOM',       file: 'shom_mangasc',              zone: 'Manche / Atlantique NE',        short: 'SHOM HYCOM' },
  { key: 'bt_120h',     label: 'Barotropic 120h',  file: 'barotropic_manche_120h',    zone: 'Manche / Atlantique',           short: 'Barotropic' },
  { key: 'bt_ne',       label: 'Barotropic NE',    file: 'barotropic_atlantique_ne',  zone: 'Atlantique Nord-Est',           short: 'Barotropic' },
  { key: 'bt_finistere',label: 'Barotropic Finistère HR', file: 'barotropic_finistere_hr', zone: 'Finistère (haute résolution)', short: 'Barotropic' },
  { key: 'bt_hycom',    label: 'Barotropic Hycom', file: 'barotropic_hycom_manche',   zone: 'Manche / Atlantique (Hycom)',   short: 'Barotropic' },
];

function computeBoatPosition(routeResult, routeDepAbsH, currentTimeH) {
  const tl = routeResult?.time_list;
  const rt = routeResult?.route;
  if (!tl || !rt || tl.length < 2 || routeDepAbsH === null) return null;
  const base = routeDepAbsH - tl[0];
  const tlAbs = tl.map(t => base + t);
  if (currentTimeH <= tlAbs[0]) return rt[0];
  if (currentTimeH >= tlAbs[tlAbs.length - 1]) return rt[rt.length - 1];
  for (let i = 0; i < tlAbs.length - 1; i++) {
    if (currentTimeH >= tlAbs[i] && currentTimeH < tlAbs[i + 1]) {
      const f = (currentTimeH - tlAbs[i]) / (tlAbs[i + 1] - tlAbs[i]);
      return [rt[i][0] + f * (rt[i + 1][0] - rt[i][0]),
              rt[i][1] + f * (rt[i + 1][1] - rt[i][1])];
    }
  }
  return null;
}

function parseWindResp(resp) {
  if (Array.isArray(resp.lat)) {
    return resp.lat.map((lat, i) => ({ lat, lon: resp.lon[i], dir: resp.dir[i], speed: resp.speed[i] }));
  }
  return resp.data || [];
}

const SAVED_ROUTE_COLORS = [
  [79,  195, 247, 230],
  [255, 183,  77, 230],
  [174, 213, 129, 230],
  [186, 104, 200, 230],
  [255, 112,  67, 230],
];


const card = {
  background: 'rgba(8,13,30,0.92)', backdropFilter: 'blur(12px)',
  border: '1px solid rgba(100,160,255,0.15)',
  borderRadius: 12, padding: '14px 16px',
  color: '#c8d8ff', boxShadow: '0 8px 32px rgba(0,0,30,0.6)',
  fontFamily: FONT, minWidth: 268,
};

const inputStyle = {
  width: '100%', padding: '7px 10px', borderRadius: 7, fontSize: 12,
  background: 'rgba(30,40,80,0.8)', color: '#c8d8ff',
  border: '1px solid rgba(100,160,255,0.2)', outline: 'none',
};

const labelStyle = { fontSize: 11, opacity: 0.65, marginBottom: 5, display: 'block' };

// ── Petit composant point (départ/arrivée) ──────────────────────────────────
function PointRow({ label, point, mode, activeMode, onToggle, accentColor, onManualSet }) {
  const active = activeMode === mode;
  const [latDeg, setLatDeg] = useState('');
  const [latMin, setLatMin] = useState('');
  const [latHem, setLatHem] = useState('N');
  const [lonDeg, setLonDeg] = useState('');
  const [lonMin, setLonMin] = useState('');
  const [lonHem, setLonHem] = useState('E');
  const anyFocused = useRef(false);

  useEffect(() => {
    if (anyFocused.current) return;
    if (!point) {
      setLatDeg(''); setLatMin(''); setLatHem('N');
      setLonDeg(''); setLonMin(''); setLonHem('E');
      return;
    }
    const aLat = Math.abs(point[0]), aLon = Math.abs(point[1]);
    setLatDeg(String(Math.floor(aLat)));
    setLatMin(((aLat % 1) * 60).toFixed(4));
    setLatHem(point[0] >= 0 ? 'N' : 'S');
    setLonDeg(String(Math.floor(aLon)));
    setLonMin(((aLon % 1) * 60).toFixed(4));
    setLonHem(point[1] >= 0 ? 'E' : 'W');
  }, [point]);

  const trySet = (ld, lm, lh, od, om, oh) => {
    const ldN = parseInt(ld, 10);
    const lmN = parseFloat(lm);
    const odN = parseInt(od, 10);
    const omN = parseFloat(om);
    if ([ldN, lmN, odN, omN].some(isNaN)) return;
    if (ldN < 0 || ldN > 90  || lmN < 0 || lmN >= 60) return;
    if (odN < 0 || odN > 180 || omN < 0 || omN >= 60) return;
    onManualSet([(ldN + lmN / 60) * (lh === 'N' ? 1 : -1),
                 (odN + omN / 60) * (oh === 'E' ? 1 : -1)]);
  };

  const fo = {
    onFocus: () => { anyFocused.current = true; },
    onBlur:  () => { anyFocused.current = false; },
  };
  const is = { ...inputStyle, padding: '3px 4px', fontSize: 13, color: point ? accentColor : 'rgba(200,216,255,0.5)', MozAppearance: 'textfield', WebkitAppearance: 'none' };
  const hemStyle = { padding: '3px 5px', borderRadius: 4, fontSize: 10, cursor: 'pointer', fontFamily: FONT, background: 'rgba(30,40,80,0.8)', color: accentColor, border: `1px solid rgba(100,160,255,0.2)` };

  return (
    <div style={{ marginBottom: 10 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 5 }}>
        <span style={{ fontSize: 11, opacity: 0.65 }}>{label}</span>
        <button onClick={onToggle} style={{
          fontSize: 10, padding: '3px 9px', borderRadius: 5, cursor: 'pointer',
          background: active ? accentColor : 'rgba(30,40,80,0.8)',
          color: active ? '#080d1a' : '#c8d8ff',
          border: `1px solid ${active ? accentColor : 'rgba(100,160,255,0.2)'}`,
          fontFamily: FONT, transition: 'all 0.15s',
        }}>
          {active ? '⊕ cliquer…' : 'Pointer'}
        </button>
      </div>
      <div style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
        {/* Lat */}
        <input type="number" min={0} max={90}  step={1} value={latDeg} placeholder="—"
          {...fo} style={{ ...is, width: 46 }}
          onChange={e => { setLatDeg(e.target.value); trySet(e.target.value, latMin, latHem, lonDeg, lonMin, lonHem); }} />
        <span style={{ fontSize: 11, opacity: 0.35 }}>°</span>
        <input type="number" min={0} max={59.9999} step={0.0001} value={latMin} placeholder="—"
          {...fo} style={{ ...is, width: 46 }}
          onChange={e => { setLatMin(e.target.value); trySet(latDeg, e.target.value, latHem, lonDeg, lonMin, lonHem); }} />
        <span style={{ fontSize: 11, opacity: 0.35 }}>'</span>
        <button style={hemStyle} onClick={() => { const h = latHem === 'N' ? 'S' : 'N'; setLatHem(h); trySet(latDeg, latMin, h, lonDeg, lonMin, lonHem); }}>{latHem}</button>
        <span style={{ opacity: 0.2, fontSize: 10 }}>|</span>
        {/* Lon */}
        <input type="number" min={0} max={180} step={1} value={lonDeg} placeholder="—"
          {...fo} style={{ ...is, width: 46 }}
          onChange={e => { setLonDeg(e.target.value); trySet(latDeg, latMin, latHem, e.target.value, lonMin, lonHem); }} />
        <span style={{ fontSize: 11, opacity: 0.35 }}>°</span>
        <input type="number" min={0} max={59.9999} step={0.0001} value={lonMin} placeholder="—"
          {...fo} style={{ ...is, width: 46 }}
          onChange={e => { setLonMin(e.target.value); trySet(latDeg, latMin, latHem, lonDeg, e.target.value, lonHem); }} />
        <span style={{ fontSize: 11, opacity: 0.35 }}>'</span>
        <button style={hemStyle} onClick={() => { const h = lonHem === 'E' ? 'W' : 'E'; setLonHem(h); trySet(latDeg, latMin, latHem, lonDeg, lonMin, h); }}>{lonHem}</button>
      </div>
    </div>
  );
}

// ── App principale ──────────────────────────────────────────────────────────
export default function App() {
  // Wind
  const [files, setFiles]           = useState([]);
  const [file, setFile]             = useState('');
  const [meta, setMeta]             = useState(null);
  const [currentTimeH, setCurrentTimeH] = useState(0);  // heure absolue (offset GRIB)
  const [showGrib, setShowGrib]         = useState(true);
  const [showParticles, setShowParticles] = useState(true);
  const [windData, setWindData]     = useState([]);
  const [windLoading, setWindLoading] = useState(false);
  const [metaLoading, setMetaLoading] = useState(false);
  const [windError, setWindError]     = useState(null);
  const [viewState, setViewState]   = useState(INIT_VIEW);

  // Flèches courant : ~160 colonnes visibles
  const autoStep = useMemo(() => {
    const lonSpan = (360 / Math.pow(2, viewState.zoom)) * (window.innerWidth / 256);
    return Math.max(0.05, lonSpan / 160);
  }, [viewState.zoom]);

  // Raster vent : résolution proche du natif modèle (0.25°), ~320 colonnes max
  const windStep = useMemo(() => {
    const lonSpan = (360 / Math.pow(2, viewState.zoom)) * (window.innerWidth / 256);
    return Math.max(0.25, lonSpan / 320);
  }, [viewState.zoom]);

  const viewport = useMemo(() => {
    const { longitude, latitude, zoom } = viewState;
    const lonSpan = (360 / Math.pow(2, zoom)) * (window.innerWidth  / 256) * 1.4;
    const latSpan = (360 / Math.pow(2, zoom)) * (window.innerHeight / 256) * 1.4;
    return {
      lon0: Math.max(-180, longitude - lonSpan / 2),
      lon1: Math.min( 180, longitude + lonSpan / 2),
      lat0: Math.max( -85, latitude  - latSpan / 2),
      lat1: Math.min(  85, latitude  + latSpan / 2),
    };
  }, [viewState]);

  const [windMode, setWindMode]       = useState('ecmwf');  // 'grib' | 'ecmwf' | 'uniform'
  const [uniformWind, setUniformWind] = useState({ direction: 270, force: 15 });
  const [uniformWindData, setUniformWindData] = useState([]);

  // Current
  const [currentFiles, setCurrentFiles]         = useState([]);
  const [currentFile, setCurrentFile]           = useState('');
  const [currentMeta, setCurrentMeta]           = useState(null);
  const [currentData, setCurrentData]           = useState([]);
  const [currentLoading, setCurrentLoading]     = useState(false);
  const [currentMetaLoading, setCurrentMetaLoading] = useState(false);
  const [currentError, setCurrentError]         = useState(null);
  const [showCurrent, setShowCurrent]           = useState(false);
  const [currentMode, setCurrentMode]           = useState('models');  // 'models' | 'grib' | 'uniform'
  const [uniformCurrent, setUniformCurrent]     = useState({ direction: 180, force: 0.5 });
  const [uniformCurrentData, setUniformCurrentData]     = useState([]);

  // Routing
  const [propulsionMode, setPropulsionMode] = useState('voile'); // 'voile' | 'moteur'
  const [motorSpeed, setMotorSpeed]         = useState(7);
  const [polaires, setPolaires]     = useState([]);
  const [polaire, setPolaire]       = useState('');
  const [depPoint, setDepPoint]     = useState(null);  // [lat, lon]
  const [arrPoint, setArrPoint]     = useState(null);
  const [clickMode, setClickMode]   = useState(null);  // 'dep' | 'arr' | null
  const [routeResult, setRouteResult] = useState(null);
  const [routeDepAbsH, setRouteDepAbsH] = useState(null); // heure absolue de départ du routage
  const [routing, setRouting]       = useState(false);
  const [routeError, setRouteError] = useState(null);
  const [routingProgress, setRoutingProgress] = useState(0);
  const [depTimeIdx, setDepTimeIdx] = useState(0);    // index dans meta.times
  const [showIsochrones, setShowIsochrones] = useState(true);
  const [showCurrentRoute, setShowCurrentRoute] = useState(true);
  const [advOpen, setAdvOpen]       = useState(false);
  const [params, setParams]         = useState({
    dt: 1, n: 100, ang_deg: 90, dang_deg: 0.3, seuil_nm: 20, facteur_raf: 2,
  });
  const [polarPct, setPolarPct]     = useState(100);

  // Cartographie terrestre
  const [landData, setLandData] = useState(null);

  // Caches viewport : évite les re-fetch pour le même timestep au même zoom
  const windCache       = useRef(new Map());
  const windPrefetching = useRef(new Set());
  const curCache        = useRef(new Map());
  const curPrefetching  = useRef(new Set());

  // Saved routes
  const [savedRoutes, setSavedRoutes] = useState([]);
  const [routeCounter, setRouteCounter] = useState(1);
  const [activeRouteIds, setActiveRouteIds] = useState(new Set());
  const [selectedRouteId, setSelectedRouteId] = useState(null);
  const [savedIsoVisible, setSavedIsoVisible] = useState({});

  // Vide les caches quand la source change (viewport différent ou nouveau fichier)
  useEffect(() => { windCache.current.clear(); windPrefetching.current.clear(); }, [file]);
  useEffect(() => { curCache.current.clear();  curPrefetching.current.clear();  }, [currentFile]);

  // ── Init ──────────────────────────────────────────────────────────────
  useEffect(() => {
    fetch(`${API}/files`).then(r => r.json()).then(setFiles).catch(() => {});
    fetch(`${API}/current-files`).then(r => r.json()).then(setCurrentFiles).catch(() => {});
    fetch(`${API}/polaires`).then(r => r.json()).then(setPolaires).catch(() => {});
  }, []);

  // ── Terre Natural Earth ───────────────────────────────────────────────
  // Résolution avec hystérésis : seuils d'entrée ≠ seuils de sortie pour
  // éviter les basculements répétés au voisinage d'un seuil.
  //   110m  →(zoom≥3)→  50m  →(zoom≥5.5)→  10m
  //   110m  ←(zoom<2)←  50m  ←(zoom<4.5)←  10m
  const [landResolution, setLandResolution] = useState(() => {
    const z = INIT_VIEW.zoom;
    return z < 3 ? '110m' : z < 5.5 ? '50m' : '10m';
  });

  useEffect(() => {
    const z = viewState.zoom;
    setLandResolution(res => {
      if (res === '110m' && z >= 3)   return '50m';
      if (res === '50m'  && z <  2)   return '110m';
      if (res === '50m'  && z >= 5.5) return '10m';
      if (res === '10m'  && z <  4.5) return '50m';
      return res;
    });
  }, [viewState.zoom]);

  useEffect(() => {
    const { lat0, lat1, lon0, lon1 } = viewport;
    const timer = setTimeout(() => {
      fetch(`${API}/land/geojson?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&resolution=${landResolution}`)
        .then(r => r.json())
        .then(setLandData)
        .catch(() => {});
    }, 150);
    return () => clearTimeout(timer);
  }, [viewport, landResolution]);

  useEffect(() => {
    if (!file) { setMeta(null); setWindData([]); setWindError(null); return; }
    // Si un modèle était déjà chargé, on change juste de source : la vue carte
    // et l'heure du curseur restent inchangées (seul le changement de fichier initial les initialise).
    const hadMeta = meta !== null;
    setMeta(null); setWindData([]); setWindError(null); setMetaLoading(true);
    fetch(`${API}/wind/${encodeURIComponent(file)}/meta`)
      .then(async r => {
        if (!r.ok) throw new Error((await r.json().catch(() => ({}))).detail ?? r.statusText);
        return r.json();
      })
      .then(m => {
        setMeta(m);
        // currentTimeH = heures absolues depuis l'époque Unix
        const refH = new Date(m.valid_times[0]).getTime() / 3600000;
        if (!hadMeta) {
          setDepTimeIdx(0);
          setCurrentTimeH(refH + (m.times[0] ?? 0));
          // Les modèles ECMWF/GFS ont un bbox global (centre ~0,0, golfe de Guinée) :
          // leur centrage se fait via centerOnEuropeOnce(), pas sur ce bbox.
          if (!file.startsWith('openmeteo_')) {
            const [lon0, lat0, lon1, lat1] = m.bbox;
            setViewState(v => ({ ...v, longitude: (lon0 + lon1) / 2, latitude: (lat0 + lat1) / 2, zoom: 4 }));
          }
        } else {
          // Recale depTimeIdx sur l'heure du curseur déjà en place ; si cette heure
          // sort de la fenêtre du nouveau modèle, recale aussi le curseur lui-même
          // sur la date la plus proche disponible (au lieu de laisser la carte vide).
          setCurrentTimeH(t => {
            let best = 0, bestDiff = Infinity;
            m.times.forEach((tm, i) => {
              const d = Math.abs((refH + tm) - t);
              if (d < bestDiff) { bestDiff = d; best = i; }
            });
            setDepTimeIdx(best);
            const tMin = refH + m.times[0];
            const tMax = refH + m.times[m.times.length - 1];
            return (t < tMin || t > tMax) ? refH + m.times[best] : t;
          });
        }
      })
      .catch(e => setWindError(e.message))
      .finally(() => setMetaLoading(false));
  }, [file]);

  // Heures absolues (depuis époque Unix) des références vent et courant
  const windRefH = useMemo(
    () => meta ? new Date(meta.valid_times[0]).getTime() / 3600000 : null,
    [meta],
  );
  const curRefH = useMemo(
    () => currentMeta ? new Date(currentMeta.valid_times[0]).getTime() / 3600000 : null,
    [currentMeta],
  );

  const { sliderMin, sliderMax } = useMemo(() => {
    const candidates = [
      ...(meta        && windRefH !== null ? [windRefH + meta.times[0],        windRefH + meta.times[meta.times.length - 1]]               : []),
      ...(currentMeta && curRefH  !== null ? [curRefH  + currentMeta.times[0], curRefH  + currentMeta.times[currentMeta.times.length - 1]] : []),
    ];
    if (!candidates.length) return { sliderMin: 0, sliderMax: 0 };
    return { sliderMin: Math.min(...candidates), sliderMax: Math.max(...candidates) };
  }, [meta, windRefH, currentMeta, curRefH]);

  // Index GRIB le plus proche de l'heure courante du slider
  const nearestGribIdx = useMemo(() => {
    if (!meta || windRefH === null) return 0;
    const windT = currentTimeH - windRefH;
    let best = 0, bestDiff = Infinity;
    meta.times.forEach((t, i) => {
      const d = Math.abs(t - windT);
      if (d < bestDiff) { bestDiff = d; best = i; }
    });
    return best;
  }, [currentTimeH, windRefH, meta]);

  useEffect(() => {
    if (!showGrib || !file || !meta || windRefH === null) { setWindLoading(false); return; }
    const tWind = currentTimeH - windRefH;
    if (tWind < meta.times[0] || tWind > meta.times[meta.times.length - 1]) {
      setWindData([]); setWindLoading(false); return;
    }
    const { lat0, lat1, lon0, lon1 } = viewport;
    const cKey = `${file}|${Math.round(tWind)}|${lat0.toFixed(1)}|${lat1.toFixed(1)}|${lon0.toFixed(1)}|${lon1.toFixed(1)}|${windStep.toFixed(2)}`;

    if (windCache.current.has(cKey)) {
      setWindData(windCache.current.get(cKey));
      return;
    }

    setWindLoading(true);
    const ctrl = new AbortController();
    const timer = setTimeout(() => {
      fetch(`${API}/wind/${encodeURIComponent(file)}/grid?t=${tWind}&lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${windStep}`, { signal: ctrl.signal })
        .then(r => r.json())
        .then(resp => {
          const data = parseWindResp(resp);
          setWindData(data);
          setWindLoading(false);
          windCache.current.set(cKey, data);
          if (windCache.current.size > 40) windCache.current.delete(windCache.current.keys().next().value);

          // Précharge les pas de temps adjacents en arrière-plan
          setTimeout(() => {
            if (!showGrib) return;
            [nearestGribIdx - 1, nearestGribIdx + 1, nearestGribIdx - 2, nearestGribIdx + 2].forEach(idx => {
              if (idx < 0 || idx >= meta.times.length) return;
              const tPre = meta.times[idx];
              const pKey = `${file}|${idx}|${lat0.toFixed(1)}|${lat1.toFixed(1)}|${lon0.toFixed(1)}|${lon1.toFixed(1)}|${windStep.toFixed(2)}`;
              if (windCache.current.has(pKey) || windPrefetching.current.has(pKey)) return;
              windPrefetching.current.add(pKey);
              fetch(`${API}/wind/${encodeURIComponent(file)}/grid?t=${tPre}&lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${windStep}`)
                .then(r => r.json())
                .then(r2 => {
                  windPrefetching.current.delete(pKey);
                  const d = parseWindResp(r2);
                  windCache.current.set(pKey, d);
                  if (windCache.current.size > 40) windCache.current.delete(windCache.current.keys().next().value);
                }).catch(() => windPrefetching.current.delete(pKey));
            });
          }, 200);
        })
        .catch(() => setWindLoading(false));
    }, 80);
    return () => { clearTimeout(timer); ctrl.abort(); };
  }, [showGrib, file, meta, currentTimeH, windRefH, viewport, windStep, nearestGribIdx]);

  // ── Courant : meta ────────────────────────────────────────────────────
  useEffect(() => {
    if (!currentFile) { setCurrentMeta(null); setCurrentData([]); setCurrentError(null); return; }
    const hadCurrentMeta = currentMeta !== null;
    setCurrentMeta(null); setCurrentData([]); setCurrentError(null); setCurrentMetaLoading(true);
    fetch(`${API}/current/${encodeURIComponent(currentFile)}/meta`)
      .then(async r => {
        if (!r.ok) throw new Error((await r.json().catch(() => ({}))).detail ?? r.statusText);
        return r.json();
      })
      .then(m => {
        setCurrentMeta(m);
        const refH = new Date(m.valid_times[0]).getTime() / 3600000;
        if (!hadCurrentMeta && !meta) {
          // Tout premier chargement, sans vent déjà affiché : initialise au début du modèle
          setCurrentTimeH(refH + (m.times[0] ?? 0));
        } else {
          // Si l'heure du curseur sort de la fenêtre du nouveau modèle, recale sur
          // la date la plus proche disponible (au lieu de laisser la carte vide).
          setCurrentTimeH(t => {
            const tMin = refH + m.times[0];
            const tMax = refH + m.times[m.times.length - 1];
            if (t >= tMin && t <= tMax) return t;
            let best = 0, bestDiff = Infinity;
            m.times.forEach((tm, i) => {
              const d = Math.abs((refH + tm) - t);
              if (d < bestDiff) { bestDiff = d; best = i; }
            });
            return refH + m.times[best];
          });
        }
      })
      .catch(e => setCurrentError(e.message))
      .finally(() => setCurrentMetaLoading(false));
  }, [currentFile]);  // eslint-disable-line react-hooks/exhaustive-deps

  // ── Courant : données interpolées ─────────────────────────────────────
  useEffect(() => {
    if (!showCurrent || !currentFile || !currentMeta || curRefH === null) { setCurrentLoading(false); return; }
    const tCur = currentTimeH - curRefH;
    const tMin = currentMeta.times[0];
    const tMax = currentMeta.times[currentMeta.times.length - 1];
    if (tCur < tMin || tCur > tMax) {
      setCurrentData([]); setCurrentLoading(false); return;
    }
    const { lat0, lat1, lon0, lon1 } = viewport;
    // Index le plus proche dans meta courant (pour clé cache stable)
    let bestCurIdx = 0, bestCurDiff = Infinity;
    currentMeta.times.forEach((tm, i) => { const d = Math.abs(tm - tCur); if (d < bestCurDiff) { bestCurDiff = d; bestCurIdx = i; } });
    const cKey = `${currentFile}|${bestCurIdx}|${lat0.toFixed(1)}|${lat1.toFixed(1)}|${lon0.toFixed(1)}|${lon1.toFixed(1)}|${autoStep.toFixed(2)}`;

    if (curCache.current.has(cKey)) {
      setCurrentData(curCache.current.get(cKey));
      return;
    }

    setCurrentLoading(true);
    const ctrl = new AbortController();
    const timer = setTimeout(() => {
      fetch(`${API}/current/${encodeURIComponent(currentFile)}/grid?t=${tCur}&lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}`, { signal: ctrl.signal })
        .then(r => r.json())
        .then(resp => {
          const data = parseWindResp(resp);
          setCurrentData(data);
          setCurrentLoading(false);
          curCache.current.set(cKey, data);
          if (curCache.current.size > 40) curCache.current.delete(curCache.current.keys().next().value);

          setTimeout(() => {
            if (!showCurrent) return;
            [bestCurIdx - 1, bestCurIdx + 1].forEach(idx => {
              if (idx < 0 || idx >= currentMeta.times.length) return;
              const tPre = currentMeta.times[idx];
              const pKey = `${currentFile}|${idx}|${lat0.toFixed(1)}|${lat1.toFixed(1)}|${lon0.toFixed(1)}|${lon1.toFixed(1)}|${autoStep.toFixed(2)}`;
              if (curCache.current.has(pKey) || curPrefetching.current.has(pKey)) return;
              curPrefetching.current.add(pKey);
              fetch(`${API}/current/${encodeURIComponent(currentFile)}/grid?t=${tPre}&lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}`)
                .then(r => r.json())
                .then(r2 => {
                  curPrefetching.current.delete(pKey);
                  const d = parseWindResp(r2);
                  curCache.current.set(pKey, d);
                  if (curCache.current.size > 40) curCache.current.delete(curCache.current.keys().next().value);
                }).catch(() => curPrefetching.current.delete(pKey));
            });
          }, 200);
        })
        .catch(() => setCurrentLoading(false));
    }, 80);
    return () => { clearTimeout(timer); ctrl.abort(); };
  }, [showCurrent, currentFile, currentMeta, currentTimeH, curRefH, viewport, autoStep]);

  // ── Reset routeResult si départ/arrivée changent ─────────────────────
  useEffect(() => { setRouteResult(null); }, [depPoint, arrPoint]);

  // ── Pré-sélection dt ──────────────────────────────────────────────────
  useEffect(() => {
    if (!depPoint || !arrPoint) return;
    if (propulsionMode === 'moteur') {
      const [lat1, lon1] = depPoint;
      const [lat2, lon2] = arrPoint;
      const R = 3440.065;
      const φ1 = lat1 * Math.PI / 180, φ2 = lat2 * Math.PI / 180;
      const Δφ = (lat2 - lat1) * Math.PI / 180;
      const Δλ = (lon2 - lon1) * Math.PI / 180;
      const a = Math.sin(Δφ / 2) ** 2 + Math.cos(φ1) * Math.cos(φ2) * Math.sin(Δλ / 2) ** 2;
      const dist_nm = 2 * R * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
      const dt = dist_nm / (30 * motorSpeed);
      const seuil_nm = Math.max(5, Math.round(dist_nm / 10));
      setParams(p => ({ ...p, dt, seuil_nm }));
      return;
    }
    if (!polaire) return;
    const [lat1, lon1] = depPoint;
    const [lat2, lon2] = arrPoint;
    const R = 3440.065; // rayon terrestre en milles nautiques
    const φ1 = lat1 * Math.PI / 180, φ2 = lat2 * Math.PI / 180;
    const Δφ = (lat2 - lat1) * Math.PI / 180;
    const Δλ = (lon2 - lon1) * Math.PI / 180;
    const a = Math.sin(Δφ / 2) ** 2 + Math.cos(φ1) * Math.cos(φ2) * Math.sin(Δλ / 2) ** 2;
    const dist_nm = 2 * R * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    const seuil_nm = Math.max(5, Math.round(dist_nm / 10));
    setParams(p => ({ ...p, seuil_nm }));
    fetch(`${API}/polaires/${encodeURIComponent(polaire)}/stats`)
      .then(r => r.json())
      .then(({ v_mean }) => {
        if (!v_mean) return;
        const dt = dist_nm / (30 * v_mean);
        setParams(p => ({ ...p, dt }));
      })
      .catch(() => {});
  }, [depPoint, arrPoint, polaire, propulsionMode, motorSpeed]);

  // ── Navigation clavier du slider de temps ────────────────────────────
  useEffect(() => {
    const handleKey = (e) => {
      const tag = document.activeElement?.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
      if (e.key === 'ArrowUp'   || e.key === 'ArrowDown')  e.preventDefault();
      if (e.key === 'ArrowRight') { e.preventDefault(); setCurrentTimeH(t => Math.min(sliderMax, t + 1)); }
      if (e.key === 'ArrowLeft')  { e.preventDefault(); setCurrentTimeH(t => Math.max(sliderMin, t - 1)); }
    };
    window.addEventListener('keydown', handleKey);
    return () => window.removeEventListener('keydown', handleKey);
  }, [sliderMin, sliderMax]);

  // ── Handlers ──────────────────────────────────────────────────────────
  const toggleSavedRoute = (id) => {
    setActiveRouteIds(prev => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id); else next.add(id);
      return next;
    });
  };

  const deleteSavedRoute = (id) => {
    setSavedRoutes(r => r.filter(s => s.id !== id));
    setActiveRouteIds(prev => { const next = new Set(prev); next.delete(id); return next; });
    setSavedIsoVisible(p => { const next = { ...p }; delete next[id]; return next; });
    if (selectedRouteId === id) setSelectedRouteId(null);
  };

  const handleMapClick = useCallback(([lat, lon]) => {
    if (clickMode === 'dep') { setDepPoint([lat, lon]); setClickMode(null); }
    if (clickMode === 'arr') { setArrPoint([lat, lon]); setClickMode(null); }
  }, [clickMode]);

  const toggleClick = (mode) =>
    setClickMode(m => m === mode ? null : mode);

  // Ne recentre sur l'Europe qu'à la toute première sélection d'un modèle météo
  // (vent ou courant) — les sélections suivantes laissent la vue où l'utilisateur l'a mise.
  const firstModelSelectRef = useRef(true);
  const centerOnEuropeOnce = () => {
    if (!firstModelSelectRef.current) return;
    firstModelSelectRef.current = false;
    setViewState(v => ({ ...v, longitude: INIT_VIEW.longitude, latitude: INIT_VIEW.latitude, zoom: INIT_VIEW.zoom }));
  };

  const selectWind = () => {
    if (showGrib) { setShowGrib(false); setShowParticles(false); }
    else { setShowGrib(true); setShowParticles(true); setShowCurrent(false); }
  };
  const selectCurrent = () => {
    if (showCurrent) { setShowCurrent(false); }
    else { setShowCurrent(true); setShowGrib(false); setShowParticles(false); }
  };
  // Sélection seule (jamais de désélection) : utilisé pour les clics sur les
  // champs internes des cartes, qui ne doivent pouvoir qu'activer la carte.
  const ensureWindSelected = () => { if (!showGrib) selectWind(); };
  const ensureCurrentSelected = () => { if (!showCurrent) selectCurrent(); };

  const runRouting = async () => {
    if (!depPoint || !arrPoint) return;
    if (propulsionMode === 'voile' && !polaire) return;
    if (windMode === 'grib' && !file) return;
    if (windMode === 'ecmwf' && !meta) return;
    const savedDepPoint = depPoint;
    const savedArrPoint = arrPoint;
    setRouting(true); setRouteResult(null); setRouteError(null); setRoutingProgress(0); setShowCurrentRoute(true);
    const depAbsH = (windMode === 'grib' || windMode === 'ecmwf') && windRefH !== null
      ? windRefH + (meta?.times?.[depTimeIdx] ?? 0)
      : currentMode !== 'uniform' && curRefH !== null
        ? curRefH + (currentMeta?.times?.[depTimeIdx] ?? 0)
        : currentTimeH;
    setRouteDepAbsH(depAbsH);
    try {
      const body = {
        ...(propulsionMode === 'voile'
          ? { polaire_file: polaire, polar_pct: polarPct }
          : { motor_speed: motorSpeed }),
        p_dep: depPoint, p_arr: arrPoint,
        t: (windMode === 'grib' || windMode === 'ecmwf')
          ? (meta?.times?.[depTimeIdx] ?? 0)
          : (currentMode !== 'uniform' ? (currentMeta?.times?.[depTimeIdx] ?? 0) : 0),
        ...params,
        ...(windMode === 'uniform'
          ? { wind_uniform: { direction: uniformWind.direction, force: uniformWind.force } }
          : { grib_file: file }),
        ...(currentMode === 'uniform'
          ? { courant_uniform: { direction: uniformCurrent.direction, force: uniformCurrent.force } }
          : currentFile
          ? { grib_courant_file: currentFile }
          : {}),
      };
      const res = await fetch(`${API}/routing/stream`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(body),
      });
      if (!res.ok) throw new Error((await res.json()).detail ?? res.statusText);

      const reader = res.body.getReader();
      const dec = new TextDecoder();
      let buf = '';
      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        buf += dec.decode(value, { stream: true });
        const lines = buf.split('\n');
        buf = lines.pop();
        for (const line of lines) {
          if (!line.trim()) continue;
          const msg = JSON.parse(line);
          if (msg.type === 'progress') {
            setRoutingProgress(msg.pct);
          } else if (msg.type === 'result') {
            const { type, ...data } = msg;
            setRouteResult(data);
            const id = Date.now();
            setSavedRoutes(r => [...r, {
              id,
              name: `Route ${routeCounter}`,
              routeResult: data,
              depPoint: savedDepPoint,
              arrPoint: savedArrPoint,
              routeDepAbsH: depAbsH,
              color: SAVED_ROUTE_COLORS[r.length % SAVED_ROUTE_COLORS.length],
            }]);
            setSavedIsoVisible(p => ({ ...p, [id]: true }));
            setRouteCounter(c => c + 1);
            setActiveRouteIds(prev => new Set([...prev, id]));
            setSelectedRouteId(id);
            setShowCurrentRoute(false);
          } else if (msg.type === 'error') {
            throw new Error(msg.detail);
          }
        }
      }
    } catch (e) {
      setRouteError(e.message);
    } finally {
      setRouting(false);
    }
  };

  // ── Grille vent uniforme ──────────────────────────────────────────────
  useEffect(() => {
    if (windMode !== 'uniform') { setUniformWindData([]); return; }
    if (!showGrib) return;
    const { lat0, lat1, lon0, lon1 } = viewport;
    const url = `${API}/wind/uniform/grid?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${windStep}&direction=${uniformWind.direction}&force=${uniformWind.force}`;
    const timer = setTimeout(() => {
      fetch(url).then(r => r.json()).then(({ data }) => setUniformWindData(data)).catch(() => {});
    }, 120);
    return () => clearTimeout(timer);
  }, [windMode, showGrib, viewport, uniformWind, windStep]);

  // ── Grille courant uniforme ───────────────────────────────────────────
  useEffect(() => {
    if (currentMode !== 'uniform') { setUniformCurrentData([]); return; }
    if (!showCurrent) return;
    const { lat0, lat1, lon0, lon1 } = viewport;
    const url = `${API}/wind/uniform/grid?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}&direction=${uniformCurrent.direction}&force=${uniformCurrent.force}`;
    const timer = setTimeout(() => {
      fetch(url).then(r => r.json()).then(({ data }) => setUniformCurrentData(data)).catch(() => {});
    }, 120);
    return () => clearTimeout(timer);
  }, [currentMode, showCurrent, viewport, uniformCurrent, autoStep]);

  // ── Computed ──────────────────────────────────────────────────────────
  // currentTimeH est en heures absolues depuis l'époque Unix
  const timeLabel = useMemo(() => {
    if (!meta && !currentMeta) return '—';
    return fmtDatetime(new Date(currentTimeH * 3600000).toISOString());
  }, [meta, currentMeta, currentTimeH]);

  const extraRoutes = useMemo(
    () => savedRoutes.filter(s => activeRouteIds.has(s.id)),
    [savedRoutes, activeRouteIds],
  );

  const focusedSaved = useMemo(
    () => selectedRouteId ? savedRoutes.find(s => s.id === selectedRouteId) ?? null : null,
    [selectedRouteId, savedRoutes],
  );
  const focusedResult  = focusedSaved?.routeResult  ?? routeResult;
  const focusedDepAbsH = focusedSaved?.routeDepAbsH ?? routeDepAbsH;

  const isoChecked = selectedRouteId
    ? (savedIsoVisible[selectedRouteId] ?? true)
    : showIsochrones;
  const setIsoChecked = (v) => {
    if (selectedRouteId) setSavedIsoVisible(p => ({ ...p, [selectedRouteId]: v }));
    else setShowIsochrones(v);
  };

  // ── Position du bateau sur la route focalisée (pour le HUD) ──────────
  const boatPosition = useMemo(() => {
    if (!focusedResult) return null;
    return computeBoatPosition(focusedResult, focusedDepAbsH, currentTimeH);
  }, [focusedResult, focusedDepAbsH, currentTimeH]);

  // ── Infos bateau à l'instant courant ──────────────────────────────────
  const boatInfo = useMemo(() => {
    const tl = focusedResult?.time_list;
    const rt = focusedResult?.route;
    if (!tl || !rt || tl.length < 2 || !boatPosition || focusedDepAbsH === null) return null;
    const base = focusedDepAbsH - tl[0];
    const tlAbs = tl.map(t => base + t);

    let i = 0;
    for (; i < tlAbs.length - 1; i++) {
      if (currentTimeH >= tlAbs[i] && currentTimeH < tlAbs[i + 1]) break;
    }
    if (i >= tlAbs.length - 1) return null;

    const [lon1, lat1] = rt[i];
    const [lon2, lat2] = rt[i + 1];
    const dt = tlAbs[i + 1] - tlAbs[i];

    // Cap (Nord = 0°, Est = 90°)
    const midLatRad = ((lat1 + lat2) / 2) * Math.PI / 180;
    const dx = (lon2 - lon1) * Math.cos(midLatRad);
    const dy = lat2 - lat1;
    const heading = (Math.atan2(dx, dy) * 180 / Math.PI + 360) % 360;

    // Vitesse sur l'eau (nœuds)
    const speed = Math.sqrt((dy * 60) ** 2 + (dx * 60) ** 2) / dt;
    if (!isFinite(speed)) return null;

    // Vent au point le plus proche du bateau
    const [bLon, bLat] = boatPosition;
    let windSpeed = null, windDir = null;
    if (windMode === 'uniform') {
      windSpeed = uniformWind.force;
      windDir   = uniformWind.direction;
    } else if (windData.length) {
      let minD = Infinity;
      for (const w of windData) {
        const d = (w.lat - bLat) ** 2 + (w.lon - bLon) ** 2;
        if (d < minD) { minD = d; windSpeed = w.speed; windDir = w.dir; }
      }
    }

    // Angle au vent (0° = face au vent, 180° = vent arrière)
    const angVent = windDir != null
      ? Math.abs((((windDir - heading + 180) % 360) + 360) % 360 - 180)
      : null;

    return { speed, windSpeed, angVent };
  }, [focusedResult, boatPosition, currentTimeH, windMode, uniformWind, windData, focusedDepAbsH]);

  const activeIsochrones = useMemo(() => {
    const all = [];
    if (showCurrentRoute && showIsochrones && routeResult?.isochrones)
      all.push(...routeResult.isochrones);
    extraRoutes.forEach(s => {
      if ((savedIsoVisible[s.id] ?? true) && s.routeResult?.isochrones)
        all.push(...s.routeResult.isochrones);
    });
    return all;
  }, [showCurrentRoute, showIsochrones, routeResult, extraRoutes, savedIsoVisible]);

  const boats = useMemo(() => {
    const result = [];
    if (showCurrentRoute && routeResult) {
      const pos = computeBoatPosition(routeResult, routeDepAbsH, currentTimeH);
      if (pos) result.push({ pos, color: [255, 255, 255, 255], outline: [80, 160, 255, 255] });
    }
    extraRoutes.forEach(saved => {
      const pos = computeBoatPosition(saved.routeResult, saved.routeDepAbsH, currentTimeH);
      if (pos) result.push({ pos, color: [...saved.color.slice(0, 3), 255], outline: [255, 255, 255, 200] });
    });
    return result;
  }, [showCurrentRoute, routeResult, routeDepAbsH, extraRoutes, currentTimeH]);

  const windT = windRefH !== null ? (currentTimeH - windRefH) : 0;
  const timeOffset = meta
    ? (() => {
        const label = windMode === 'ecmwf' ? 'MODEL' : 'GRIB';
        const base  = `T+${windT.toFixed(1)}h  (${label}: T+${meta.times[nearestGribIdx]}h)`;
        if (windMode === 'ecmwf' && meta.run_time) {
          const d   = new Date(meta.run_time);
          const run = `${String(d.getUTCDate()).padStart(2,'0')}/${String(d.getUTCMonth()+1).padStart(2,'0')} ${String(d.getUTCHours()).padStart(2,'0')}h`;
          return `${base}  •  run ${run} UTC`;
        }
        return base;
      })()
    : '';
  const canRoute   = (windMode === 'uniform' || !!file) && !!depPoint && !!arrPoint &&
    (propulsionMode === 'moteur' || !!polaire) && !routing;

  return (
    <div style={{ width: '100vw', height: '100vh', position: 'relative', background: '#080d1a', fontFamily: FONT }}>
      <WindMap
        data={windMode === 'uniform' ? uniformWindData : windData}
        viewState={viewState}
        onViewStateChange={setViewState}
        depPoint={depPoint}
        arrPoint={arrPoint}
        route={showCurrentRoute ? (routeResult?.route ?? null) : null}
        boats={boats}
        isochrones={activeIsochrones}
        showIsochrones={showIsochrones}
        showGrib={showGrib}
        currentData={currentMode === 'uniform' ? uniformCurrentData : currentData}
        showCurrent={showCurrent}
        clickMode={clickMode}
        onMapClick={handleMapClick}
        extraRoutes={extraRoutes}
        landData={landData}
        showParticles={showParticles}
      />

      {/* ══ HUD bateau ══════════════════════════════════════════════════════ */}
      {focusedResult && boatInfo && (
        <div style={{
          position: 'absolute', top: 16, left: '50%', transform: 'translateX(-50%)',
          zIndex: 10, pointerEvents: 'none',
          background: 'rgba(8,13,30,0.88)', backdropFilter: 'blur(12px)',
          border: '1px solid rgba(100,160,255,0.15)', borderRadius: 12,
          padding: '10px 28px', boxShadow: '0 8px 32px rgba(0,0,30,0.6)',
          fontFamily: FONT, display: 'flex', gap: 36, alignItems: 'center',
          whiteSpace: 'nowrap',
        }}>
          {[
            ['Vit. bateau',   `${boatInfo.speed.toFixed(1)} kt`],
            ['Vit. vent',     boatInfo.windSpeed != null ? `${boatInfo.windSpeed.toFixed(1)} kt` : '—'],
            ['Angle au vent', boatInfo.angVent   != null ? `${Math.round(boatInfo.angVent)}°`    : '—'],
          ].map(([lbl, val]) => (
            <div key={lbl} style={{ textAlign: 'center' }}>
              <div style={{ fontSize: 10, opacity: 0.5, marginBottom: 2, color: '#c8d8ff' }}>{lbl}</div>
              <div style={{ fontSize: 18, color: '#4fc3f7', fontWeight: 'bold', fontVariantNumeric: 'tabular-nums' }}>
                {val}
              </div>
            </div>
          ))}
        </div>
      )}

      {/* ══ Sidebar gauche ══════════════════════════════════════════════════ */}
      <div style={{
        position: 'absolute', top: 0, left: 16, zIndex: 10,
        display: 'flex', flexDirection: 'column', gap: 10,
        width: 300,
        maxHeight: `calc(100vh - 16px - ${(meta || currentMeta) ? 40 : 0}px)`,
        overflowY: 'auto', scrollbarWidth: 'none',
      }}>

        {/* ── Carte vent ───────────────────────────────────────────── */}
        <div onClick={selectWind} style={{ ...card, opacity: showGrib ? 1 : 0.55, cursor: 'pointer', transition: 'opacity 0.15s' }}>
          <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase', marginBottom: 12 }}>
            Vent
          </div>

          <div
            onClickCapture={e => { if (!showGrib) { e.stopPropagation(); ensureWindSelected(); } }}
            onClick={e => e.stopPropagation()}
          >
          {/* Toggle MODELS / GRIB / Uniforme */}
          <div style={{ display: 'flex', gap: 6, marginBottom: 12 }}>
            {[['ecmwf', 'MODELS'], ['grib', 'GRIB'], ['uniform', 'Uniforme']].map(([mode, label]) => (
              <button key={mode} onClick={() => {
                setWindMode(mode);
                if (mode === 'ecmwf') setFile('');
              }}
                style={{
                  flex: 1, padding: '5px 0', borderRadius: 6, fontSize: 11,
                  fontFamily: FONT, cursor: 'pointer',
                  background: windMode === mode ? '#60a5fa' : 'rgba(30,40,80,0.8)',
                  color: windMode === mode ? '#080d1a' : '#c8d8ff',
                  border: `1px solid ${windMode === mode ? '#60a5fa' : 'rgba(100,160,255,0.2)'}`,
                  transition: 'all 0.15s',
                }}>
                {label}
              </button>
            ))}
          </div>

          {windMode === 'ecmwf' ? (
            <>
              {/* Sélecteur de modèle */}
              {[
                ['ecmwf', 'ECMWF IFS', 'openmeteo_ecmwf'],
                ['gfs',   'GFS 0.25°', 'openmeteo_gfs'],
              ].map(([key, label, virtualFile]) => (
                <button key={key} onClick={() => {
                  setFile(virtualFile);
                  centerOnEuropeOnce();
                }}
                  style={{
                    display: 'block', width: '100%', marginBottom: 5,
                    padding: '4px 0', borderRadius: 5, fontSize: 10,
                    fontFamily: FONT, cursor: 'pointer',
                    background: file === virtualFile ? '#93c5fd' : 'rgba(30,40,80,0.8)',
                    color: file === virtualFile ? '#080d1a' : '#c8d8ff',
                    border: `1px solid ${file === virtualFile ? '#93c5fd' : 'rgba(100,160,255,0.2)'}`,
                    transition: 'all 0.15s',
                  }}>
                  {label}
                </button>
              ))}
              {/* Statut */}
              <div style={{
                marginBottom: 10, fontSize: 11, opacity: 0.7, textAlign: 'center', lineHeight: 1.5,
                color: windError ? '#ff8080' : undefined,
              }}>
                {metaLoading || windLoading
                  ? `⏳ Téléchargement ${file === 'openmeteo_gfs' ? 'GFS' : 'ECMWF'}… (~20s)`
                  : windError
                  ? `⚠ ${windError}`
                  : meta
                  ? (<>
                      {`✓ ${meta.model} — pas ${meta.times.length > 1 ? meta.times[1] - meta.times[0] : '?'}h (${meta.days}j)`}
                      {meta.run_time && (() => {
                        const d = new Date(meta.run_time);
                        const run = `${String(d.getUTCDate()).padStart(2,'0')}/${String(d.getUTCMonth()+1).padStart(2,'0')} ${String(d.getUTCHours()).padStart(2,'0')}h UTC`;
                        return <span style={{ display: 'block', opacity: 0.6, fontSize: 10 }}>Run : {run}</span>;
                      })()}
                    </>)
                  : null}
              </div>
            </>
          ) : windMode === 'grib' ? (
            <>
              <label style={labelStyle}>Fichier</label>
              <select value={file} onChange={e => setFile(e.target.value)}
                style={{ ...inputStyle, marginBottom: 12 }}>
                <option value="">— Choisir un fichier —</option>
                {files.filter(f => !f.startsWith('openmeteo_')).map(f => <option key={f} value={f}>{f}</option>)}
              </select>
              {windLoading && <div style={{ marginTop: 8, fontSize: 10, opacity: 0.4, textAlign: 'center' }}>Chargement…</div>}
            </>
          ) : (
            <>
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '6px 10px', marginBottom: 10 }}>
                <label style={{ fontSize: 11, opacity: 0.7, display: 'flex', flexDirection: 'column', gap: 4 }}>
                  Direction (°)
                  <input type="number" min={0} max={359} step={5} value={uniformWind.direction}
                    onChange={e => setUniformWind(w => ({ ...w, direction: +e.target.value }))}
                    style={{ ...inputStyle, padding: '5px 8px', fontSize: 12 }} />
                </label>
                <label style={{ fontSize: 11, opacity: 0.7, display: 'flex', flexDirection: 'column', gap: 4 }}>
                  Force (nœuds)
                  <input type="number" min={0} max={60} step={1} value={uniformWind.force}
                    onChange={e => setUniformWind(w => ({ ...w, force: +e.target.value }))}
                    style={{ ...inputStyle, padding: '5px 8px', fontSize: 12 }} />
                </label>
              </div>
            </>
          )}
          </div>
        </div>

        {/* ── Carte courant ────────────────────────────────────────── */}
        <div onClick={selectCurrent} style={{ ...card, opacity: showCurrent ? 1 : 0.55, cursor: 'pointer', transition: 'opacity 0.15s' }}>
          <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase', marginBottom: 12 }}>
            Courant
          </div>

          <div
            onClickCapture={e => { if (!showCurrent) { e.stopPropagation(); ensureCurrentSelected(); } }}
            onClick={e => e.stopPropagation()}
          >
          {/* Toggle MODELS / GRIB / Uniforme */}
          <div style={{ display: 'flex', gap: 6, marginBottom: 12 }}>
            {[['models', 'MODELS'], ['grib', 'GRIB'], ['uniform', 'Uniforme']].map(([mode, label]) => (
              <button key={mode} onClick={() => {
                setCurrentMode(mode);
                if (mode === 'models') setCurrentFile('');
                if (mode === 'grib') setCurrentFile('');
              }}
                style={{
                  flex: 1, padding: '5px 0', borderRadius: 6, fontSize: 11,
                  fontFamily: FONT, cursor: 'pointer',
                  background: currentMode === mode ? '#60a5fa' : 'rgba(30,40,80,0.8)',
                  color: currentMode === mode ? '#080d1a' : '#c8d8ff',
                  border: `1px solid ${currentMode === mode ? '#60a5fa' : 'rgba(100,160,255,0.2)'}`,
                  transition: 'all 0.15s',
                }}>
                {label}
              </button>
            ))}
          </div>

          {currentMode === 'models' ? (
            <>
              {/* Sélecteur de modèle */}
              <div style={{ maxHeight: 168, overflowY: 'auto', marginBottom: 8 }}>
                {CURRENT_MODELS.map(({ key, label, file: virtualFile }) => (
                  <button key={key} onClick={() => { setCurrentFile(virtualFile); centerOnEuropeOnce(); }}
                    style={{
                      display: 'block', width: '100%', marginBottom: 5,
                      padding: '4px 0', borderRadius: 5, fontSize: 10,
                      fontFamily: FONT, cursor: 'pointer',
                      background: currentFile === virtualFile ? '#93c5fd' : 'rgba(30,40,80,0.8)',
                      color: currentFile === virtualFile ? '#080d1a' : '#c8d8ff',
                      border: `1px solid ${currentFile === virtualFile ? '#93c5fd' : 'rgba(100,160,255,0.2)'}`,
                      transition: 'all 0.15s',
                    }}>
                    {label}
                  </button>
                ))}
              </div>
              {/* Zone fixe du modèle direct */}
              {(() => {
                const selectedModel = CURRENT_MODELS.find(m => m.file === currentFile);
                if (!selectedModel) return null;
                return (
                  <>
                    <div style={{ marginBottom: 8 }}>
                      <div style={{
                        fontSize: 10, padding: '4px 8px', borderRadius: 5,
                        background: currentMeta ? 'rgba(20,85,164,0.25)' : 'rgba(30,40,80,0.8)',
                        border: `1px solid ${currentMeta ? 'rgba(100,160,255,0.35)' : 'rgba(100,160,255,0.15)'}`,
                        cursor: 'default',
                        color: currentMeta ? '#93c5fd' : '#c8d8ff',
                      }}>
                        {selectedModel.zone}
                      </div>
                    </div>
                    {/* Statut */}
                    <div style={{
                      fontSize: 11, opacity: 0.7, textAlign: 'center', lineHeight: 1.5,
                      color: currentError ? '#ff8080' : undefined,
                    }}>
                      {currentMetaLoading || currentLoading
                        ? `⏳ Téléchargement ${selectedModel.short}…`
                        : currentError
                        ? `⚠ ${currentError}`
                        : currentMeta
                        ? (<>
                            {(() => {
                              const stepH = currentMeta.times.length > 1 ? currentMeta.times[1] - currentMeta.times[0] : null;
                              const stepLabel = stepH === null ? '?' : stepH < 1 ? `${Math.round(stepH * 60)}min` : `${stepH}h`;
                              return `✓ ${currentMeta.model} — ${stepLabel} (${currentMeta.days}j)`;
                            })()}
                            {currentMeta.run_time && (() => {
                              const d = new Date(currentMeta.run_time);
                              const run = `${String(d.getUTCDate()).padStart(2,'0')}/${String(d.getUTCMonth()+1).padStart(2,'0')} ${String(d.getUTCHours()).padStart(2,'0')}h UTC`;
                              return <span style={{ display: 'block', opacity: 0.6, fontSize: 10 }}>Depuis : {run}</span>;
                            })()}
                          </>)
                        : null}
                    </div>
                  </>
                );
              })()}
            </>
          ) : currentMode === 'grib' ? (
            <>
              <label style={labelStyle}>Fichier</label>
              <select value={currentFile} onChange={e => setCurrentFile(e.target.value)}
                style={{ ...inputStyle, marginBottom: 12 }}>
                <option value="">— Choisir un fichier —</option>
                {currentFiles.filter(f => !f.startsWith('shom_') && !f.startsWith('barotropic_')).map(f => <option key={f} value={f}>{f}</option>)}
              </select>
              {currentLoading && <div style={{ marginTop: 8, fontSize: 10, opacity: 0.4, textAlign: 'center' }}>Chargement…</div>}
            </>
          ) : (
            <>
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '6px 10px', marginBottom: 10 }}>
                <label style={{ fontSize: 11, opacity: 0.7, display: 'flex', flexDirection: 'column', gap: 4 }}>
                  Direction (°)
                  <input type="number" min={0} max={359} step={5} value={uniformCurrent.direction}
                    onChange={e => setUniformCurrent(c => ({ ...c, direction: +e.target.value }))}
                    style={{ ...inputStyle, padding: '5px 8px', fontSize: 12 }} />
                </label>
                <label style={{ fontSize: 11, opacity: 0.7, display: 'flex', flexDirection: 'column', gap: 4 }}>
                  Force (nœuds)
                  <input type="number" min={0} max={10} step={0.1} value={uniformCurrent.force}
                    onChange={e => setUniformCurrent(c => ({ ...c, force: +e.target.value }))}
                    style={{ ...inputStyle, padding: '5px 8px', fontSize: 12 }} />
                </label>
              </div>
            </>
          )}
          </div>
        </div>

        {/* ── Carte polaire ─────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase', marginBottom: 12 }}>
            Polaire
</div>

          {/* Toggle Voile / Moteur */}
          <div style={{ display: 'flex', gap: 6, marginBottom: 12 }}>
            {['voile', 'moteur'].map(mode => (
              <button key={mode} onClick={() => setPropulsionMode(mode)}
                style={{
                  flex: 1, padding: '5px 0', borderRadius: 6, fontSize: 11,
                  fontFamily: FONT, cursor: 'pointer',
                  background: propulsionMode === mode ? '#60a5fa' : 'rgba(30,40,80,0.8)',
                  color: propulsionMode === mode ? '#080d1a' : '#c8d8ff',
                  border: `1px solid ${propulsionMode === mode ? '#60a5fa' : 'rgba(100,160,255,0.2)'}`,
                  transition: 'all 0.15s',
                }}>
                {mode === 'voile' ? 'Voile' : 'Moteur'}
              </button>
            ))}
          </div>

          {propulsionMode === 'voile' ? (
            <>
              <label style={labelStyle}>Fichier polaire</label>
              <select value={polaire} onChange={e => setPolaire(e.target.value)}
                style={{ ...inputStyle, marginBottom: 12 }}>
                <option value="">— Choisir la polaire —</option>
                {polaires.map(p => <option key={p} value={p}>{p.replace('.csv', '')}</option>)}
              </select>
              <label style={labelStyle}>
                Performance polaire&nbsp;
                <strong style={{ color: polarPct < 100 ? '#ff8080' : polarPct > 100 ? '#3ddc84' : '#60a5fa' }}>
                  {polarPct}%
                </strong>
              </label>
              <input type="range" min={50} max={150} step={5} value={polarPct}
                onChange={e => setPolarPct(+e.target.value)}
                style={{ width: '100%', accentColor: '#60a5fa', cursor: 'pointer' }} />
            </>
          ) : (
            <>
              <label style={labelStyle}>
                Vitesse moteur&nbsp;
                <strong style={{ color: '#60a5fa' }}>{motorSpeed} kt</strong>
              </label>
              <input type="range" min={1} max={30} step={0.5} value={motorSpeed}
                onChange={e => setMotorSpeed(+e.target.value)}
                style={{ width: '100%', accentColor: '#60a5fa', cursor: 'pointer' }} />
            </>
          )}
        </div>

        {/* ── Carte routage ────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
            <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase' }}>Routage</div>
            {(routeResult || focusedSaved) && (
              <label style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 11, opacity: 0.7, cursor: 'pointer' }}>
                <input type="checkbox" checked={isoChecked} onChange={e => setIsoChecked(e.target.checked)}
                  style={{ accentColor: '#4fc3f7', cursor: 'pointer' }} />
                Isochrones
              </label>
            )}
          </div>

          <PointRow
            label="Départ" mode="dep" activeMode={clickMode}
            point={depPoint} accentColor="#3ddc84"
            onToggle={() => toggleClick('dep')}
            onManualSet={setDepPoint}
          />
          <PointRow
            label="Arrivée" mode="arr" activeMode={clickMode}
            point={arrPoint} accentColor="#ff6b6b"
            onToggle={() => toggleClick('arr')}
            onManualSet={setArrPoint}
          />

          {/* Date de départ */}
          {(() => {
            const useCurrentGrib = windMode === 'uniform' && currentMode === 'grib';
            const depMeta   = useCurrentGrib ? currentMeta : meta;
            const depRefH   = useCurrentGrib ? curRefH : windRefH;
            return (
              <>
                <label style={labelStyle}>Date de départ</label>
                <select
                  value={depTimeIdx}
                  onChange={e => {
                    const idx = +e.target.value;
                    setDepTimeIdx(idx);
                    if (depRefH !== null && depMeta) setCurrentTimeH(depRefH + (depMeta.times[idx] ?? 0));
                  }}
                  style={{ ...inputStyle, marginBottom: 10 }}
                  disabled={!depMeta}
                >
                  {depMeta?.valid_times?.map((vt, i) => (
                    <option key={i} value={i}>{fmtDatetime(vt)}</option>
                  )) ?? <option>— Choisir date de départ —</option>}
                </select>
              </>
            );
          })()}

          {/* Paramètres avancés */}
          <button onClick={() => setAdvOpen(v => !v)}
            style={{
              background: 'none', border: 'none', cursor: 'pointer',
              color: 'rgba(200,216,255,0.45)', fontSize: 11, fontFamily: FONT,
              padding: 0, marginBottom: advOpen ? 10 : 4,
              display: 'flex', alignItems: 'center', gap: 5,
            }}>
            {advOpen ? '▾' : '▸'} Paramètres avancés
          </button>

          {advOpen && (
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '6px 10px', marginBottom: 10 }}>
              {[
                ['dt (h)',       'dt',       0.25, 6,   0.25],
                ['n caps',       'n',        20,   360, 10  ],
                ['ang init (°)', 'ang_deg',  10,   180, 5   ],
                ['dang/pas (°)', 'dang_deg', 0,    2,   0.05],
                ['seuil (NM)',   'seuil_nm',   5, 5000, 10  ],
                ['facteur raf',  'facteur_raf', 1,  20,  1  ],
              ].map(([lbl, key, min, max, step]) => (
                <label key={key} style={{ fontSize: 10, opacity: 0.7, display: 'flex', flexDirection: 'column', gap: 2 }}>
                  {lbl}
                  <input type="number" min={min} max={max} step={step} value={params[key]}
                    onChange={e => setParams(p => ({ ...p, [key]: +e.target.value }))}
                    style={{ ...inputStyle, padding: '4px 6px', fontSize: 11 }} />
                </label>
              ))}
            </div>
          )}

          {/* Bouton lancer */}
          <div style={{ display: 'flex', gap: 6, marginTop: 2, marginBottom: routing ? 6 : (routeResult || routeError) ? 10 : 0 }}>
            <button onClick={runRouting} disabled={!canRoute || routing}
              style={{
                flex: 1, padding: '9px 0', borderRadius: 8,
                fontSize: 12, fontWeight: 'bold', fontFamily: FONT,
                cursor: (canRoute && !routing) ? 'pointer' : 'not-allowed',
                background: (canRoute && !routing)
                  ? 'linear-gradient(135deg, #1455a4, #1e90d8)'
                  : 'rgba(30,40,80,0.4)',
                color: (canRoute && !routing) ? 'white' : 'rgba(200,216,255,0.25)',
                border: 'none',
              }}>
              {routing ? 'Calcul…' : 'Lancer'}
            </button>
          </div>

          {/* Barre de progression */}
          {routing && (
            <div style={{ marginBottom: 10 }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: 10, marginBottom: 4 }}>
                <span style={{ opacity: 0.45 }}>Progression isochrones</span>
                <span style={{ color: '#4fc3f7', fontWeight: 'bold' }}>{routingProgress}%</span>
              </div>
              <div style={{ background: 'rgba(100,160,255,0.12)', borderRadius: 4, height: 4, overflow: 'hidden' }}>
                <div style={{
                  height: '100%', borderRadius: 4,
                  background: 'linear-gradient(90deg, #1455a4, #4fc3f7)',
                  width: `${routingProgress}%`,
                  transition: routingProgress > 0 ? 'width 0.5s ease' : 'none',
                }} />
              </div>
            </div>
          )}

          {/* Résultat */}
          {routeResult && !routing && (
            <>
              <div style={{
                background: 'rgba(20,85,164,0.2)', borderRadius: 8,
                padding: '10px 12px', border: '1px solid rgba(100,160,255,0.2)',
              }}>
                <div style={{ fontSize: 10, opacity: 0.5, marginBottom: 4 }}>Durée estimée</div>
                <div style={{ fontSize: 22, color: '#4fc3f7', fontWeight: 'bold' }}>
                  {routeResult.days > 0 && <>{routeResult.days}<span style={{ fontSize: 13, opacity: 0.7 }}>j </span></>}
                  {routeResult.hours}<span style={{ fontSize: 13, opacity: 0.7 }}>h </span>
                  {String(routeResult.minutes).padStart(2, '0')}<span style={{ fontSize: 13, opacity: 0.7 }}>min</span>
                </div>
                <div style={{ fontSize: 10, opacity: 0.4, marginTop: 6 }}>
                  Calcul : {routeResult.calc_time_s}s
                </div>
              </div>
            </>
          )}

          {routeError && (
            <div style={{ fontSize: 11, color: '#ff8080', marginTop: 4 }}>
              Erreur : {routeError}
            </div>
          )}
        </div>
      </div>

{/* ══ Barre de temps ══════════════════════════════════════════════════ */}
      {(meta || currentMeta) && (() => {
        const sliderRange = sliderMax - sliderMin;

        // Bandes GRIB pour l'indicateur (GRIB uniquement, pas uniforme)
        const windBand = (meta && windRefH !== null && windMode === 'grib') ? {
          t0: windRefH + meta.times[0],
          t1: windRefH + meta.times[meta.times.length - 1],
          color: '#e8f2ff', opacity: 0.95, key: 'vent',
        } : null;
        const curBand = (currentMeta && curRefH !== null && currentMode !== 'uniform') ? {
          t0: curRefH + currentMeta.times[0],
          t1: curRefH + currentMeta.times[currentMeta.times.length - 1],
          color: '#7fa8c0', opacity: 0.70, key: 'courant',
        } : null;

        const pct = t => ((t - sliderMin) / sliderRange * 100);

        return (
          <div style={{
            position: 'absolute', bottom: 0, left: 0, right: 0, zIndex: 10,
            background: 'rgba(8,13,30,0.92)', backdropFilter: 'blur(12px)',
            borderTop: '1px solid rgba(100,160,255,0.15)',
            padding: '8px 24px 8px', fontFamily: FONT,
          }}>
            {/* En-tête */}
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', marginBottom: 6 }}>
              <strong style={{ color: '#4fc3f7', fontSize: 13 }}>{timeLabel}</strong>
            </div>

            {/* ── Barre unique : piste colorée + thumb blanc + input invisible ── */}
            <div style={{ position: 'relative', height: 20, marginTop: 2 }}>
              {/* Piste visuelle avec bandes GRIB */}
              <div style={{
                position: 'absolute', top: '50%', transform: 'translateY(-50%)',
                left: 0, right: 0, height: 6,
                background: 'rgba(100,160,255,0.12)', borderRadius: 3, overflow: 'hidden',
              }}>
                {sliderRange > 0 && [windBand, curBand].filter(Boolean).map(b => (
                  <div key={b.key} style={{
                    position: 'absolute',
                    left:  `${pct(b.t0)}%`,
                    width: `${pct(b.t1) - pct(b.t0)}%`,
                    top: 0, height: '100%',
                    background: b.color, opacity: b.opacity, borderRadius: 3,
                  }} />
                ))}
              </div>
              {/* Thumb : cercle blanc */}
              <div style={{
                position: 'absolute', top: '50%',
                left: `${sliderRange > 0 ? pct(currentTimeH) : 0}%`,
                transform: 'translate(-50%, -50%)',
                width: 12, height: 12,
                background: 'white', borderRadius: '50%',
                boxShadow: '0 1px 4px rgba(0,0,0,0.5)',
                pointerEvents: 'none', zIndex: 2,
              }} />
              {/* Input range invisible par-dessus pour l'interaction */}
              <input
                type="range"
                min={sliderMin} max={sliderMax}
                step={1}
                value={currentTimeH}
                onChange={e => setCurrentTimeH(+e.target.value)}
                style={{
                  position: 'absolute', top: 0, left: 0,
                  width: '100%', height: '100%',
                  opacity: 0, cursor: 'pointer',
                  margin: 0, padding: 0,
                  WebkitAppearance: 'none', appearance: 'none',
                }}
              />
            </div>

          </div>
        );
      })()}

      {/* ══ Panneau routes sauvegardées (droite) ══════════════════════════ */}
      {savedRoutes.length > 0 && (
        <div style={{
          position: 'absolute', top: 0, right: 16, zIndex: 10,
          width: 260,
          maxHeight: `calc(100vh - 16px - ${(meta || currentMeta) ? 56 : 0}px)`,
          overflowY: 'auto', scrollbarWidth: 'none',
        }}>
          <div style={{ ...card, minWidth: 0 }}>
            <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase', marginBottom: 12 }}>
              Routes sauvegardées
            </div>
            {savedRoutes.map((saved, idx) => {
              const active   = activeRouteIds.has(saved.id);
              const selected = selectedRouteId === saved.id;
              const [r, g, b] = saved.color;
              const accent = `rgb(${r},${g},${b})`;
              return (
                <div key={saved.id}
                  onClick={() => setSelectedRouteId(id => id === saved.id ? null : saved.id)}
                  style={{
                    background: selected ? `rgba(${r},${g},${b},0.18)` : active ? `rgba(${r},${g},${b},0.07)` : 'rgba(20,30,60,0.6)',
                    borderRadius: 8,
                    padding: '9px 10px',
                    marginBottom: idx < savedRoutes.length - 1 ? 8 : 0,
                    border: `${selected ? 2 : 1}px solid ${selected ? accent : active ? `rgba(${r},${g},${b},0.55)` : 'rgba(100,160,255,0.15)'}`,
                    transition: 'border 0.2s, background 0.2s',
                    cursor: 'pointer',
                  }}>
                  <div style={{ display: 'flex', alignItems: 'center', gap: 6, marginBottom: 5 }}>
                    <div style={{
                      width: 8, height: 8, borderRadius: '50%', flexShrink: 0,
                      background: accent, opacity: active ? 1 : 0.35,
                      transition: 'opacity 0.2s',
                    }} />
                    <div style={{ fontSize: 12, color: selected ? accent : active ? accent : '#7ec8f7', fontWeight: 'bold' }}>
                      {saved.name}
                    </div>
                  </div>
                  <div style={{ fontSize: 10, opacity: 0.55, marginBottom: 1 }}>
                    Dép. {fmtCoord(saved.depPoint)}
                  </div>
                  <div style={{ fontSize: 10, opacity: 0.55, marginBottom: 6 }}>
                    Arr. {fmtCoord(saved.arrPoint)}
                  </div>
                  <div style={{ fontSize: 13, color: '#c8d8ff', fontWeight: 'bold', marginBottom: 8, fontVariantNumeric: 'tabular-nums' }}>
                    {saved.routeResult.days > 0 && (
                      <>{saved.routeResult.days}<span style={{ fontSize: 10, opacity: 0.6 }}>j </span></>
                    )}
                    {saved.routeResult.hours}<span style={{ fontSize: 10, opacity: 0.6 }}>h </span>
                    {String(saved.routeResult.minutes).padStart(2, '0')}<span style={{ fontSize: 10, opacity: 0.6 }}>min</span>
                  </div>
                  <div style={{ display: 'flex', gap: 6 }}>
                    <button onClick={() => toggleSavedRoute(saved.id)} style={{
                      flex: 1, padding: '5px 0', borderRadius: 5, fontSize: 10,
                      fontFamily: FONT, cursor: 'pointer',
                      background: active ? `rgba(${r},${g},${b},0.25)` : 'rgba(20,85,164,0.45)',
                      color: active ? accent : '#7ec8f7',
                      border: `1px solid ${active ? `rgba(${r},${g},${b},0.5)` : 'rgba(100,160,255,0.3)'}`,
                      transition: 'all 0.15s',
                    }}>
                      {active ? 'Masquer' : 'Afficher'}
                    </button>
                    <button onClick={() => deleteSavedRoute(saved.id)} style={{
                      flex: 1, padding: '5px 0', borderRadius: 5, fontSize: 10,
                      fontFamily: FONT, cursor: 'pointer',
                      background: 'rgba(180,30,30,0.25)',
                      color: '#ff8080',
                      border: '1px solid rgba(255,100,100,0.3)',
                      transition: 'all 0.15s',
                    }}>
                      Supprimer
                    </button>
                  </div>
                </div>
              );
            })}
          </div>
        </div>
      )}

    </div>
  );
}
