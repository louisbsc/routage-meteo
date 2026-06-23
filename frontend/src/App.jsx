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

function computeBoatPosition(routeResult, routeDepAbsH, currentTimeH) {
  const tl = routeResult?.time_list;
  const rt = routeResult?.route;
  if (!tl || !rt || tl.length < 2 || routeDepAbsH === null) return null;
  const base = routeDepAbsH - tl[0];
  const tlAbs = tl.map(t => base + t);
  if (currentTimeH <= tlAbs[0]) return rt[0];
  if (currentTimeH >= tlAbs[tlAbs.length - 1]) return null;
  for (let i = 0; i < tlAbs.length - 1; i++) {
    if (currentTimeH >= tlAbs[i] && currentTimeH < tlAbs[i + 1]) {
      const f = (currentTimeH - tlAbs[i]) / (tlAbs[i + 1] - tlAbs[i]);
      return [rt[i][0] + f * (rt[i + 1][0] - rt[i][0]),
              rt[i][1] + f * (rt[i + 1][1] - rt[i][1])];
    }
  }
  return null;
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
  const [showGrib, setShowGrib]     = useState(true);
  const [windData, setWindData]     = useState([]);
  const [windLoading, setWindLoading] = useState(false);
  const [viewState, setViewState]   = useState(INIT_VIEW);

  // ~15 colonnes de flèches visibles quelle que soit le zoom
  const autoStep = useMemo(() => {
    const lonSpan = (360 / Math.pow(2, viewState.zoom)) * (window.innerWidth / 256);
    return Math.max(0.05, lonSpan / 160);
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

  const [windMode, setWindMode]       = useState('grib');  // 'grib' | 'uniform'
  const [uniformWind, setUniformWind] = useState({ direction: 270, force: 15 });
  const [uniformWindData, setUniformWindData] = useState([]);

  // Current
  const [currentFiles, setCurrentFiles]         = useState([]);
  const [currentFile, setCurrentFile]           = useState('');
  const [currentMeta, setCurrentMeta]           = useState(null);
  const [currentData, setCurrentData]           = useState([]);
  const [currentLoading, setCurrentLoading]     = useState(false);
  const [showCurrent, setShowCurrent]           = useState(true);
  const [currentMode, setCurrentMode]           = useState('grib');  // 'grib' | 'uniform'
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
    dt: 1, n: 100, ang_deg: 90, dang_deg: 0.3, seuil_nm: 20, facteur_raf: 5,
  });
  const [polarPct, setPolarPct]     = useState(100);

  // Cartographie terrestre
  const [landData, setLandData] = useState(null);

  // Saved routes
  const [savedRoutes, setSavedRoutes] = useState([]);
  const [routeCounter, setRouteCounter] = useState(1);
  const [activeRouteIds, setActiveRouteIds] = useState(new Set());
  const [selectedRouteId, setSelectedRouteId] = useState(null);
  const [savedIsoVisible, setSavedIsoVisible] = useState({});

  // ── Init ──────────────────────────────────────────────────────────────
  useEffect(() => {
    fetch(`${API}/files`).then(r => r.json()).then(setFiles).catch(() => {});
    fetch(`${API}/current-files`).then(r => r.json()).then(setCurrentFiles).catch(() => {});
    fetch(`${API}/polaires`).then(r => r.json()).then(setPolaires).catch(() => {});
  }, []);

  // ── Terre Natural Earth ───────────────────────────────────────────────
  useEffect(() => {
    const { lat0, lat1, lon0, lon1 } = viewport;
    const timer = setTimeout(() => {
      fetch(`${API}/land/geojson?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}`)
        .then(r => r.json())
        .then(setLandData)
        .catch(() => {});
    }, 150);
    return () => clearTimeout(timer);
  }, [viewport]);

  useEffect(() => {
    if (!file) { setMeta(null); setWindData([]); return; }
    setMeta(null); setWindData([]); setRouteResult(null);
    fetch(`${API}/wind/${encodeURIComponent(file)}/meta`)
      .then(r => r.json())
      .then(m => {
        setMeta(m);
        setDepTimeIdx(0);
        // currentTimeH = heures absolues depuis l'époque Unix
        const refH = new Date(m.valid_times[0]).getTime() / 3600000;
        setCurrentTimeH(refH + (m.times[0] ?? 0));
        const [lon0, lat0, lon1, lat1] = m.bbox;
        setViewState(v => ({ ...v, longitude: (lon0 + lon1) / 2, latitude: (lat0 + lat1) / 2, zoom: 4 }));
      }).catch(() => {});
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
    if (!file || !meta || windRefH === null) return;
    const tWind = currentTimeH - windRefH;
    // Masquer si hors de la plage du GRIB vent
    if (tWind < meta.times[0] || tWind > meta.times[meta.times.length - 1]) {
      setWindData([]);
      setWindLoading(false);
      return;
    }
    setWindLoading(true);
    const { lat0, lat1, lon0, lon1 } = viewport;
    const timer = setTimeout(() => {
      fetch(`${API}/wind/${encodeURIComponent(file)}/grid?t=${tWind}&lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}`)
        .then(r => r.json())
        .then(({ data }) => { setWindData(data); setWindLoading(false); })
        .catch(() => setWindLoading(false));
    }, 80);
    return () => clearTimeout(timer);
  }, [file, meta, currentTimeH, windRefH, viewport, autoStep]);

  // ── Courant : meta ────────────────────────────────────────────────────
  useEffect(() => {
    if (!currentFile) { setCurrentMeta(null); setCurrentData([]); return; }
    setCurrentMeta(null); setCurrentData([]);
    fetch(`${API}/current/${encodeURIComponent(currentFile)}/meta`)
      .then(r => r.json())
      .then(m => {
        setCurrentMeta(m);
        // N'initialise le slider que si aucun vent n'est chargé
        if (!meta) {
          const refH = new Date(m.valid_times[0]).getTime() / 3600000;
          setCurrentTimeH(refH + (m.times[0] ?? 0));
        }
      })
      .catch(() => {});
  }, [currentFile]);  // eslint-disable-line react-hooks/exhaustive-deps

  // ── Courant : données interpolées ─────────────────────────────────────
  useEffect(() => {
    if (!currentFile || !currentMeta || curRefH === null) return;
    const tCur = currentTimeH - curRefH;
    const tMin = currentMeta.times[0];
    const tMax = currentMeta.times[currentMeta.times.length - 1];
    // Masquer si hors de la plage du GRIB courant
    if (tCur < tMin || tCur > tMax) {
      setCurrentData([]);
      setCurrentLoading(false);
      return;
    }
    setCurrentLoading(true);
    const { lat0, lat1, lon0, lon1 } = viewport;
    const timer = setTimeout(() => {
      fetch(`${API}/current/${encodeURIComponent(currentFile)}/grid?t=${tCur}&lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}`)
        .then(r => r.json())
        .then(({ data }) => { setCurrentData(data); setCurrentLoading(false); })
        .catch(() => setCurrentLoading(false));
    }, 80);
    return () => clearTimeout(timer);
  }, [currentFile, currentMeta, currentTimeH, curRefH, viewport, autoStep]);

  // ── Reset routeResult si paramètres de départ changent ───────────────
  useEffect(() => { setRouteResult(null); }, [depPoint, arrPoint, depTimeIdx]);

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
      const dt_raw = dist_nm / (50 * motorSpeed);
      const dt = Math.min(6, Math.max(0.25, Math.round(dt_raw / 0.25) * 0.25));
      const seuil_nm = Math.max(5, Math.round(dist_nm / 20));
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
    const seuil_nm = Math.max(5, Math.round(dist_nm / 20));
    setParams(p => ({ ...p, seuil_nm }));
    fetch(`${API}/polaires/${encodeURIComponent(polaire)}/stats`)
      .then(r => r.json())
      .then(({ v_mean }) => {
        if (!v_mean) return;
        const dt_raw = dist_nm / (50 * v_mean);
        const dt = Math.min(6, Math.max(0.25, Math.round(dt_raw / 0.25) * 0.25));
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

  const runRouting = async () => {
    if (!depPoint || !arrPoint) return;
    if (propulsionMode === 'voile' && !polaire) return;
    if (windMode === 'grib' && !file) return;
    const savedDepPoint = depPoint;
    const savedArrPoint = arrPoint;
    setRouting(true); setRouteResult(null); setRouteError(null); setRoutingProgress(0); setShowCurrentRoute(true);
    const depAbsH = windMode === 'grib' && windRefH !== null
      ? windRefH + (meta?.times?.[depTimeIdx] ?? 0)
      : currentMode === 'grib' && curRefH !== null
        ? curRefH + (currentMeta?.times?.[depTimeIdx] ?? 0)
        : currentTimeH;
    setRouteDepAbsH(depAbsH);
    try {
      const body = {
        ...(propulsionMode === 'voile'
          ? { polaire_file: polaire, polar_pct: polarPct }
          : { motor_speed: motorSpeed }),
        p_dep: depPoint, p_arr: arrPoint,
        t: windMode === 'grib'
          ? (meta?.times?.[depTimeIdx] ?? 0)
          : (currentMode === 'grib' ? (currentMeta?.times?.[depTimeIdx] ?? 0) : 0),
        ...params,
        ...(routeResult ? { route: routeResult.route } : {}),
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
    const { lat0, lat1, lon0, lon1 } = viewport;
    const url = `${API}/wind/uniform/grid?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}&direction=${uniformWind.direction}&force=${uniformWind.force}`;
    const timer = setTimeout(() => {
      fetch(url).then(r => r.json()).then(({ data }) => setUniformWindData(data)).catch(() => {});
    }, 120);
    return () => clearTimeout(timer);
  }, [windMode, viewport, uniformWind, autoStep]);

  // ── Grille courant uniforme ───────────────────────────────────────────
  useEffect(() => {
    if (currentMode !== 'uniform') { setUniformCurrentData([]); return; }
    const { lat0, lat1, lon0, lon1 } = viewport;
    const url = `${API}/wind/uniform/grid?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${autoStep}&direction=${uniformCurrent.direction}&force=${uniformCurrent.force}`;
    const timer = setTimeout(() => {
      fetch(url).then(r => r.json()).then(({ data }) => setUniformCurrentData(data)).catch(() => {});
    }, 120);
    return () => clearTimeout(timer);
  }, [currentMode, viewport, uniformCurrent, autoStep]);

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
    ? `T+${windT.toFixed(1)}h  (GRIB: T+${meta.times[nearestGribIdx]}h)`
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
        <div style={card}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
            <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase' }}>Vent</div>
            <label style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 11, opacity: 0.7, cursor: 'pointer' }}>
              <input type="checkbox" checked={showGrib} onChange={e => setShowGrib(e.target.checked)}
                style={{ accentColor: '#60a5fa', cursor: 'pointer' }} />
              Afficher
            </label>
          </div>

          {/* Toggle GRIB / Uniforme */}
          <div style={{ display: 'flex', gap: 6, marginBottom: 12 }}>
            {['grib', 'uniform'].map(mode => (
              <button key={mode} onClick={() => setWindMode(mode)}
                style={{
                  flex: 1, padding: '5px 0', borderRadius: 6, fontSize: 11,
                  fontFamily: FONT, cursor: 'pointer',
                  background: windMode === mode ? '#60a5fa' : 'rgba(30,40,80,0.8)',
                  color: windMode === mode ? '#080d1a' : '#c8d8ff',
                  border: `1px solid ${windMode === mode ? '#60a5fa' : 'rgba(100,160,255,0.2)'}`,
                  transition: 'all 0.15s',
                }}>
                {mode === 'grib' ? 'GRIB' : 'Uniforme'}
              </button>
            ))}
          </div>

          {windMode === 'grib' ? (
            <>
              <label style={labelStyle}>Fichier</label>
              <select value={file} onChange={e => setFile(e.target.value)}
                style={{ ...inputStyle, marginBottom: 12 }}>
                <option value="">— Choisir un fichier —</option>
                {files.map(f => <option key={f} value={f}>{f}</option>)}
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

        {/* ── Carte courant ────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
            <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase' }}>Courant</div>
            <label style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 11, opacity: 0.7, cursor: 'pointer' }}>
              <input type="checkbox" checked={showCurrent} onChange={e => setShowCurrent(e.target.checked)}
                style={{ accentColor: '#60a5fa', cursor: 'pointer' }} />
              Afficher
            </label>
          </div>

          {/* Toggle GRIB / Uniforme */}
          <div style={{ display: 'flex', gap: 6, marginBottom: 12 }}>
            {['grib', 'uniform'].map(mode => (
              <button key={mode} onClick={() => setCurrentMode(mode)}
                style={{
                  flex: 1, padding: '5px 0', borderRadius: 6, fontSize: 11,
                  fontFamily: FONT, cursor: 'pointer',
                  background: currentMode === mode ? '#60a5fa' : 'rgba(30,40,80,0.8)',
                  color: currentMode === mode ? '#080d1a' : '#c8d8ff',
                  border: `1px solid ${currentMode === mode ? '#60a5fa' : 'rgba(100,160,255,0.2)'}`,
                  transition: 'all 0.15s',
                }}>
                {mode === 'grib' ? 'GRIB' : 'Uniforme'}
              </button>
            ))}
          </div>

          {currentMode === 'grib' ? (
            <>
              <label style={labelStyle}>Fichier</label>
              <select value={currentFile} onChange={e => setCurrentFile(e.target.value)}
                style={{ ...inputStyle, marginBottom: 12 }}>
                <option value="">— Choisir un fichier —</option>
                {currentFiles.map(f => <option key={f} value={f}>{f}</option>)}
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
                ['seuil (NM)',   'seuil_nm',   5,  200, 5   ],
                ['facteur raf',  'facteur_raf', 1, 20,  1   ],
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

          {/* Bouton lancer / raffiner */}
          <button onClick={runRouting} disabled={!canRoute}
            style={{
              width: '100%', padding: '9px 0', borderRadius: 8,
              fontSize: 13, fontWeight: 'bold', fontFamily: FONT,
              cursor: canRoute ? 'pointer' : 'not-allowed',
              background: canRoute
                ? (routeResult && !routing
                    ? 'linear-gradient(135deg, #0d6b3a, #1daa5e)'
                    : 'linear-gradient(135deg, #1455a4, #1e90d8)')
                : 'rgba(30,40,80,0.4)',
              color: canRoute ? 'white' : 'rgba(200,216,255,0.25)',
              border: 'none',
              marginTop: 2,
              marginBottom: routing ? 6 : (routeResult || routeError) ? 10 : 0,
            }}>
            {routing
              ? 'Calcul en cours…'
              : routeResult
                ? 'Raffiner la route'
                : 'Lancer le routage'}
          </button>

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
        const curBand = (currentMeta && curRefH !== null && currentMode === 'grib') ? {
          t0: curRefH + currentMeta.times[0],
          t1: curRefH + currentMeta.times[currentMeta.times.length - 1],
          color: '#7fa8c0', opacity: 0.70, key: 'courant',
        } : null;

        const pct = t => ((t - sliderMin) / sliderRange * 100);

        // Ticks (start + end de chaque bande), dédupliqués si confondus
        const ticks = [
          ...(windBand ? [
            { t: windBand.t0, color: windBand.color },
            { t: windBand.t1, color: windBand.color },
          ] : []),
          ...(curBand ? [
            { t: curBand.t0, color: curBand.color },
            { t: curBand.t1, color: curBand.color },
          ] : []),
        ].filter((tk, i, arr) =>
          // Dédoublonner les ticks trop proches (< 0.5 % du slider)
          arr.findIndex(o => Math.abs(pct(o.t) - pct(tk.t)) < 0.5) === i
        );

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
              <span style={{ opacity: 0.4, fontSize: 11 }}>{timeOffset}</span>
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
