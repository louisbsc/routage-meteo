import { useState, useEffect, useCallback, useMemo } from 'react';
import WindMap from './WindMap.jsx';

const API = 'http://localhost:8000';
const INIT_VIEW = { longitude: -5, latitude: 47, zoom: 4, pitch: 0, bearing: 0 };

const DAYS   = ['Dim', 'Lun', 'Mar', 'Mer', 'Jeu', 'Ven', 'Sam'];
const MONTHS = ['jan', 'fév', 'mar', 'avr', 'mai', 'juin', 'juil', 'aoû', 'sep', 'oct', 'nov', 'déc'];

function fmtDatetime(iso) {
  const d = new Date(iso);
  return `${DAYS[d.getUTCDay()]} ${String(d.getUTCDate()).padStart(2, '0')} ${MONTHS[d.getUTCMonth()]}  ${String(d.getUTCHours()).padStart(2, '0')}h UTC`;
}

function fmtCoord(pt) {
  if (!pt) return null;
  const [lat, lon] = pt;
  return `${Math.abs(lat).toFixed(2)}°${lat >= 0 ? 'N' : 'S'}  ${Math.abs(lon).toFixed(2)}°${lon >= 0 ? 'E' : 'W'}`;
}


const FONT = "'SF Mono', 'Consolas', monospace";

const STEP_OPTIONS = [
  { label: '1min',  h: 1 / 60 },
  { label: '10min', h: 10 / 60 },
  { label: '30min', h: 0.5 },
  { label: '1h',    h: 1 },
  { label: '2h',    h: 2 },
  { label: '3h',    h: 3 },
  { label: '5h',    h: 5 },
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
function PointRow({ label, point, mode, activeMode, onToggle, accentColor }) {
  const active = activeMode === mode;
  return (
    <div style={{ marginBottom: 10 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 3 }}>
        <span style={{ fontSize: 11, opacity: 0.65 }}>{label}</span>
        <button
          onClick={onToggle}
          style={{
            fontSize: 10, padding: '3px 9px', borderRadius: 5, cursor: 'pointer',
            background: active ? accentColor : 'rgba(30,40,80,0.8)',
            color: active ? '#080d1a' : '#c8d8ff',
            border: `1px solid ${active ? accentColor : 'rgba(100,160,255,0.2)'}`,
            fontFamily: FONT, transition: 'all 0.15s',
          }}
        >
          {active ? '⊕ cliquer sur la carte…' : 'Pointer'}
        </button>
      </div>
      <div style={{
        fontSize: 12, fontVariantNumeric: 'tabular-nums', letterSpacing: '0.03em',
        color: point ? accentColor : 'rgba(200,216,255,0.25)',
      }}>
        {fmtCoord(point) ?? '—'}
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
  const [stride, setStride]         = useState(2);
  const [showGrib, setShowGrib]     = useState(true);
  const [stepIdx, setStepIdx] = useState(3);  // index dans STEP_OPTIONS (défaut 1h)
  const [windData, setWindData]     = useState([]);
  const [windLoading, setWindLoading] = useState(false);
  const [viewState, setViewState]   = useState(INIT_VIEW);
  const [windMode, setWindMode]       = useState('grib');  // 'grib' | 'uniform'
  const [uniformWind, setUniformWind] = useState({ direction: 270, force: 15 });
  const [uniformStride, setUniformStride] = useState(3);  // 1=dense … 5=sparse
  const [uniformWindData, setUniformWindData] = useState([]);

  // Current
  const [currentFiles, setCurrentFiles]         = useState([]);
  const [currentFile, setCurrentFile]           = useState('');
  const [currentMeta, setCurrentMeta]           = useState(null);
  const [currentData, setCurrentData]           = useState([]);
  const [currentLoading, setCurrentLoading]     = useState(false);
  const [currentStride, setCurrentStride]       = useState(2);
  const [showCurrent, setShowCurrent]           = useState(true);
  const [currentMode, setCurrentMode]           = useState('grib');  // 'grib' | 'uniform'
  const [uniformCurrent, setUniformCurrent]     = useState({ direction: 180, force: 0.5 });
  const [uniformCurrentStride, setUniformCurrentStride] = useState(3);
  const [uniformCurrentData, setUniformCurrentData]     = useState([]);

  // Routing
  const [polaires, setPolaires]     = useState([]);
  const [polaire, setPolaire]       = useState('');
  const [depPoint, setDepPoint]     = useState(null);  // [lat, lon]
  const [arrPoint, setArrPoint]     = useState(null);
  const [clickMode, setClickMode]   = useState(null);  // 'dep' | 'arr' | null
  const [routeResult, setRouteResult] = useState(null);
  const [routing, setRouting]       = useState(false);
  const [routeError, setRouteError] = useState(null);
  const [routingProgress, setRoutingProgress] = useState(0);
  const [depTimeIdx, setDepTimeIdx] = useState(0);    // index dans meta.times
  const [showIsochrones, setShowIsochrones] = useState(true);
  const [advOpen, setAdvOpen]       = useState(false);
  const [params, setParams]         = useState({
    dt: 1, n: 100, ang_deg: 90, dang_deg: 0.3,
  });
  const [polarPct, setPolarPct]     = useState(100);

  // ── Init ──────────────────────────────────────────────────────────────
  useEffect(() => {
    fetch(`${API}/files`).then(r => r.json()).then(setFiles).catch(() => {});
    fetch(`${API}/current-files`).then(r => r.json()).then(setCurrentFiles).catch(() => {});
    fetch(`${API}/polaires`).then(r => r.json()).then(list => {
      setPolaires(list);
      if (list.length > 0) setPolaire(list[0]);
    }).catch(() => {});
  }, []);

  useEffect(() => {
    if (!file) return;
    setMeta(null); setWindData([]); setRouteResult(null);
    fetch(`${API}/wind/${encodeURIComponent(file)}/meta`)
      .then(r => r.json())
      .then(m => {
        setMeta(m); setCurrentTimeH(m.times[0] ?? 0); setDepTimeIdx(0);
        const [lon0, lat0, lon1, lat1] = m.bbox;
        setViewState(v => ({ ...v, longitude: (lon0 + lon1) / 2, latitude: (lat0 + lat1) / 2, zoom: 4 }));
      }).catch(() => {});
  }, [file]);

  // Index GRIB le plus proche de l'heure courante du slider
  const nearestGribIdx = useMemo(() => {
    if (!meta) return 0;
    let best = 0, bestDiff = Infinity;
    meta.times.forEach((t, i) => {
      const d = Math.abs(t - currentTimeH);
      if (d < bestDiff) { bestDiff = d; best = i; }
    });
    return best;
  }, [currentTimeH, meta]);

  useEffect(() => {
    if (!file || !meta) return;
    setWindLoading(true);
    const timer = setTimeout(() => {
      fetch(`${API}/wind/${encodeURIComponent(file)}/interpolated?t=${currentTimeH}&stride=${stride}`)
        .then(r => r.json())
        .then(({ data }) => { setWindData(data); setWindLoading(false); })
        .catch(() => setWindLoading(false));
    }, 80);
    return () => clearTimeout(timer);
  }, [file, meta, currentTimeH, stride]);

  // ── Courant : meta ────────────────────────────────────────────────────
  useEffect(() => {
    if (!currentFile) { setCurrentMeta(null); setCurrentData([]); return; }
    setCurrentMeta(null); setCurrentData([]);
    fetch(`${API}/current/${encodeURIComponent(currentFile)}/meta`)
      .then(r => r.json())
      .then(m => {
        setCurrentMeta(m);
        if (!meta) setCurrentTimeH(m.times[0] ?? 0);
      })
      .catch(() => {});
  }, [currentFile]);  // eslint-disable-line react-hooks/exhaustive-deps

  // ── Courant : données interpolées ─────────────────────────────────────
  useEffect(() => {
    if (!currentFile || !currentMeta) return;
    setCurrentLoading(true);
    const tMin = currentMeta.times[0];
    const tMax = currentMeta.times[currentMeta.times.length - 1];
    const tClamped = Math.max(tMin, Math.min(tMax, currentTimeH));
    const timer = setTimeout(() => {
      fetch(`${API}/current/${encodeURIComponent(currentFile)}/interpolated?t=${tClamped}&stride=${currentStride}`)
        .then(r => r.json())
        .then(({ data }) => { setCurrentData(data); setCurrentLoading(false); })
        .catch(() => setCurrentLoading(false));
    }, 80);
    return () => clearTimeout(timer);
  }, [currentFile, currentMeta, currentTimeH, currentStride]);

  // ── Pré-sélection dt ──────────────────────────────────────────────────
  useEffect(() => {
    if (!depPoint || !arrPoint || !polaire) return;
    const [lat1, lon1] = depPoint;
    const [lat2, lon2] = arrPoint;
    const R = 3440.065; // rayon terrestre en milles nautiques
    const φ1 = lat1 * Math.PI / 180, φ2 = lat2 * Math.PI / 180;
    const Δφ = (lat2 - lat1) * Math.PI / 180;
    const Δλ = (lon2 - lon1) * Math.PI / 180;
    const a = Math.sin(Δφ / 2) ** 2 + Math.cos(φ1) * Math.cos(φ2) * Math.sin(Δλ / 2) ** 2;
    const dist_nm = 2 * R * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    fetch(`${API}/polaires/${encodeURIComponent(polaire)}/stats`)
      .then(r => r.json())
      .then(({ v_mean }) => {
        if (!v_mean) return;
        const dt_raw = dist_nm / (50 * v_mean);
        const dt = Math.min(6, Math.max(0.25, Math.round(dt_raw / 0.25) * 0.25));
        setParams(p => ({ ...p, dt }));
      })
      .catch(() => {});
  }, [depPoint, arrPoint, polaire]);

  // ── Handlers ──────────────────────────────────────────────────────────
  const handleMapClick = useCallback(([lat, lon]) => {
    if (clickMode === 'dep') { setDepPoint([lat, lon]); setClickMode(null); }
    if (clickMode === 'arr') { setArrPoint([lat, lon]); setClickMode(null); }
  }, [clickMode]);

  const toggleClick = (mode) =>
    setClickMode(m => m === mode ? null : mode);

  const runRouting = async () => {
    if (!depPoint || !arrPoint || !polaire) return;
    if (windMode === 'grib' && !file) return;
    setRouting(true); setRouteResult(null); setRouteError(null); setRoutingProgress(0);
    try {
      const body = {
        polaire_file: polaire,
        p_dep: depPoint, p_arr: arrPoint,
        t: windMode === 'grib' ? (meta?.times?.[depTimeIdx] ?? 0) : 0,
        ...params,
        polar_pct: polarPct,
        ...(windMode === 'uniform'
          ? { wind_uniform: { direction: uniformWind.direction, force: uniformWind.force } }
          : { grib_file: file }),
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

  // ── Grille vent uniforme (fetch backend, ancrée en geo, sans terre) ──────
  const UNIFORM_STEPS = [4, 2, 1, 0.5, 0.25];  // degrés selon stride 1→5

  useEffect(() => {
    if (windMode !== 'uniform') { setUniformWindData([]); return; }
    const { longitude, latitude, zoom } = viewState;
    const lonSpan = (360 / Math.pow(2, zoom)) * (window.innerWidth  / 256) * 1.3;
    const latSpan = (360 / Math.pow(2, zoom)) * (window.innerHeight / 256) * 1.3;
    const lon0 = Math.max(-180, longitude - lonSpan / 2);
    const lon1 = Math.min(180,  longitude + lonSpan / 2);
    const lat0 = Math.max(-85,  latitude  - latSpan / 2);
    const lat1 = Math.min(85,   latitude  + latSpan / 2);
    const step = UNIFORM_STEPS[uniformStride - 1];
    const url = `${API}/wind/uniform/grid?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${step}&direction=${uniformWind.direction}&force=${uniformWind.force}`;
    const timer = setTimeout(() => {
      fetch(url).then(r => r.json()).then(({ data }) => setUniformWindData(data)).catch(() => {});
    }, 120);
    return () => clearTimeout(timer);
  }, [windMode, viewState, uniformWind, uniformStride]);

  // ── Grille courant uniforme ───────────────────────────────────────────
  useEffect(() => {
    if (currentMode !== 'uniform') { setUniformCurrentData([]); return; }
    const { longitude, latitude, zoom } = viewState;
    const lonSpan = (360 / Math.pow(2, zoom)) * (window.innerWidth  / 256) * 1.3;
    const latSpan = (360 / Math.pow(2, zoom)) * (window.innerHeight / 256) * 1.3;
    const lon0 = Math.max(-180, longitude - lonSpan / 2);
    const lon1 = Math.min(180,  longitude + lonSpan / 2);
    const lat0 = Math.max(-85,  latitude  - latSpan / 2);
    const lat1 = Math.min(85,   latitude  + latSpan / 2);
    const step = UNIFORM_STEPS[uniformCurrentStride - 1];
    const url = `${API}/wind/uniform/grid?lat0=${lat0}&lat1=${lat1}&lon0=${lon0}&lon1=${lon1}&step=${step}&direction=${uniformCurrent.direction}&force=${uniformCurrent.force}`;
    const timer = setTimeout(() => {
      fetch(url).then(r => r.json()).then(({ data }) => setUniformCurrentData(data)).catch(() => {});
    }, 120);
    return () => clearTimeout(timer);
  }, [currentMode, viewState, uniformCurrent, uniformCurrentStride]);

  // ── Position du bateau sur la route ───────────────────────────────────
  const boatPosition = useMemo(() => {
    const tl = routeResult?.time_list;
    const rt = routeResult?.route;
    if (!tl || !rt || tl.length < 2) return null;
    if (currentTimeH <= tl[0]) return rt[0];
    if (currentTimeH >= tl[tl.length - 1]) return null;
    for (let i = 0; i < tl.length - 1; i++) {
      if (currentTimeH >= tl[i] && currentTimeH < tl[i + 1]) {
        const f = (currentTimeH - tl[i]) / (tl[i + 1] - tl[i]);
        return [rt[i][0] + f * (rt[i + 1][0] - rt[i][0]),
                rt[i][1] + f * (rt[i + 1][1] - rt[i][1])];
      }
    }
    return null;
  }, [routeResult, currentTimeH]);

  // ── Infos bateau à l'instant courant ──────────────────────────────────
  const boatInfo = useMemo(() => {
    const tl = routeResult?.time_list;
    const rt = routeResult?.route;
    if (!tl || !rt || tl.length < 2 || !boatPosition) return null;

    let i = 0;
    for (; i < tl.length - 1; i++) {
      if (currentTimeH >= tl[i] && currentTimeH < tl[i + 1]) break;
    }
    if (i >= tl.length - 1) return null;

    const [lon1, lat1] = rt[i];
    const [lon2, lat2] = rt[i + 1];
    const dt = tl[i + 1] - tl[i];

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
  }, [routeResult, boatPosition, currentTimeH, windMode, uniformWind, windData]);

  // ── Computed ──────────────────────────────────────────────────────────
  const timeLabel = useMemo(() => {
    if (meta) {
      const ref = new Date(meta.valid_times[0]);
      const delta = currentTimeH - meta.times[0];
      ref.setUTCMinutes(ref.getUTCMinutes() + Math.round(delta * 60));
      return fmtDatetime(ref.toISOString());
    }
    if (currentMeta) return `T+${currentTimeH}h`;
    return '—';
  }, [meta, currentMeta, currentTimeH]);

  const timeOffset = meta
    ? `T+${currentTimeH}h  (GRIB: T+${meta.times[nearestGribIdx]}h)`
    : '';
  const canRoute   = (windMode === 'uniform' || !!file) && !!depPoint && !!arrPoint && !!polaire && !routing;

  return (
    <div style={{ width: '100vw', height: '100vh', position: 'relative', background: '#080d1a', fontFamily: FONT }}>
      <WindMap
        data={windMode === 'uniform' ? uniformWindData : windData}
        viewState={viewState}
        onViewStateChange={setViewState}
        depPoint={depPoint}
        arrPoint={arrPoint}
        route={routeResult?.route ?? null}
        boatPosition={boatPosition}
        isochrones={routeResult?.isochrones ?? []}
        showIsochrones={showIsochrones}
        showGrib={showGrib}
        currentData={currentMode === 'uniform' ? uniformCurrentData : currentData}
        showCurrent={showCurrent}
        clickMode={clickMode}
        onMapClick={handleMapClick}
      />

      {/* ══ HUD bateau ══════════════════════════════════════════════════════ */}
      {routeResult && boatInfo && (
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
        position: 'absolute', top: 16, left: 16, zIndex: 10,
        display: 'flex', flexDirection: 'column', gap: 10,
        maxHeight: 'calc(100vh - 80px)', overflowY: 'auto', scrollbarWidth: 'none',
      }}>

        {/* ── Carte vent ───────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
            <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase' }}>Vent</div>
            <label style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 11, opacity: 0.7, cursor: 'pointer' }}>
              <input type="checkbox" checked={showGrib} onChange={e => setShowGrib(e.target.checked)}
                style={{ accentColor: '#4fc3f7', cursor: 'pointer' }} />
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
                  background: windMode === mode ? '#4fc3f7' : 'rgba(30,40,80,0.8)',
                  color: windMode === mode ? '#080d1a' : '#c8d8ff',
                  border: `1px solid ${windMode === mode ? '#4fc3f7' : 'rgba(100,160,255,0.2)'}`,
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
              {meta && (
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0 14px' }}>
                  <div>
                    <label style={labelStyle}>
                      Densité <strong style={{ color: '#4fc3f7' }}>1/{stride}</strong>
                      <span style={{ opacity: 0.35, marginLeft: 6 }}>({windData.length.toLocaleString()})</span>
                    </label>
                    <input type="range" min={1} max={6} value={stride}
                      onChange={e => setStride(+e.target.value)}
                      style={{ width: '100%', accentColor: '#4fc3f7', cursor: 'pointer' }} />
                  </div>
                  <div>
                    <label style={labelStyle}>
                      Pas <strong style={{ color: '#4fc3f7' }}>{STEP_OPTIONS[stepIdx].label}</strong>
                    </label>
                    <input type="range" min={0} max={STEP_OPTIONS.length - 1} value={stepIdx}
                      onChange={e => { setStepIdx(+e.target.value); setCurrentTimeH(meta.times[0] ?? 0); }}
                      style={{ width: '100%', accentColor: '#4fc3f7', cursor: 'pointer' }} />
                  </div>
                </div>
              )}
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
              <label style={labelStyle}>
                Densité <strong style={{ color: '#4fc3f7' }}>{['très sparse', 'sparse', 'normale', 'dense', 'très dense'][uniformStride - 1]}</strong>
                <span style={{ opacity: 0.35, marginLeft: 6 }}>({uniformWindData.length} pts)</span>
              </label>
              <input type="range" min={1} max={5} value={uniformStride}
                onChange={e => setUniformStride(+e.target.value)}
                style={{ width: '100%', accentColor: '#4fc3f7', cursor: 'pointer' }} />
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
              {currentMeta && (
                <div>
                  <label style={labelStyle}>
                    Densité <strong style={{ color: '#60a5fa' }}>1/{currentStride}</strong>
                    <span style={{ opacity: 0.35, marginLeft: 6 }}>({currentData.length.toLocaleString()})</span>
                  </label>
                  <input type="range" min={1} max={6} value={currentStride}
                    onChange={e => setCurrentStride(+e.target.value)}
                    style={{ width: '100%', accentColor: '#60a5fa', cursor: 'pointer' }} />
                </div>
              )}
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
              <label style={labelStyle}>
                Densité <strong style={{ color: '#60a5fa' }}>{['très sparse', 'sparse', 'normale', 'dense', 'très dense'][uniformCurrentStride - 1]}</strong>
                <span style={{ opacity: 0.35, marginLeft: 6 }}>({uniformCurrentData.length} pts)</span>
              </label>
              <input type="range" min={1} max={5} value={uniformCurrentStride}
                onChange={e => setUniformCurrentStride(+e.target.value)}
                style={{ width: '100%', accentColor: '#60a5fa', cursor: 'pointer' }} />
            </>
          )}
        </div>

        {/* ── Carte polaire ─────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase', marginBottom: 12 }}>
            Polaire
          </div>
          <label style={labelStyle}>Fichier polaire</label>
          <select value={polaire} onChange={e => setPolaire(e.target.value)}
            style={{ ...inputStyle, marginBottom: 12 }}>
            {polaires.map(p => <option key={p} value={p}>{p.replace('.csv', '')}</option>)}
          </select>
          <label style={labelStyle}>
            Performance polaire&nbsp;
            <strong style={{ color: polarPct < 100 ? '#ff8080' : polarPct > 100 ? '#3ddc84' : '#4fc3f7' }}>
              {polarPct}%
            </strong>
          </label>
          <input type="range" min={50} max={150} step={5} value={polarPct}
            onChange={e => setPolarPct(+e.target.value)}
            style={{ width: '100%', accentColor: '#4fc3f7', cursor: 'pointer' }} />
        </div>

        {/* ── Carte routage ────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
            <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase' }}>Routage</div>
            {routeResult && (
              <label style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 11, opacity: 0.7, cursor: 'pointer' }}>
                <input type="checkbox" checked={showIsochrones} onChange={e => setShowIsochrones(e.target.checked)}
                  style={{ accentColor: '#4fc3f7', cursor: 'pointer' }} />
                Isochrones
              </label>
            )}
          </div>

          <PointRow
            label="Départ" mode="dep" activeMode={clickMode}
            point={depPoint} accentColor="#3ddc84"
            onToggle={() => toggleClick('dep')}
          />
          <PointRow
            label="Arrivée" mode="arr" activeMode={clickMode}
            point={arrPoint} accentColor="#ff6b6b"
            onToggle={() => toggleClick('arr')}
          />

          {/* Date de départ (GRIB uniquement) */}
          {windMode === 'grib' && (
            <>
              <label style={labelStyle}>Date de départ</label>
              <select
                value={depTimeIdx}
                onChange={e => setDepTimeIdx(+e.target.value)}
                style={{ ...inputStyle, marginBottom: 10 }}
                disabled={!meta}
              >
                {meta?.valid_times?.map((vt, i) => (
                  <option key={i} value={i}>{fmtDatetime(vt)}</option>
                )) ?? <option>— choisir un fichier GRIB —</option>}
              </select>
            </>
          )}

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
          <button onClick={runRouting} disabled={!canRoute}
            style={{
              width: '100%', padding: '9px 0', borderRadius: 8,
              fontSize: 13, fontWeight: 'bold', fontFamily: FONT,
              cursor: canRoute ? 'pointer' : 'not-allowed',
              background: canRoute
                ? 'linear-gradient(135deg, #1455a4, #1e90d8)'
                : 'rgba(30,40,80,0.4)',
              color: canRoute ? 'white' : 'rgba(200,216,255,0.25)',
              border: 'none',
              marginTop: 2,
              marginBottom: routing ? 6 : (routeResult || routeError) ? 10 : 0,
            }}>
            {routing ? 'Calcul en cours…' : 'Lancer le routage'}
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
            <div style={{
              background: 'rgba(20,85,164,0.2)', borderRadius: 8,
              padding: '10px 12px', border: '1px solid rgba(100,160,255,0.2)',
            }}>
              <div style={{ fontSize: 10, opacity: 0.5, marginBottom: 4 }}>Durée estimée</div>
              <div style={{ fontSize: 22, color: '#4fc3f7', fontWeight: 'bold' }}>
                {routeResult.days > 0 && <>{routeResult.days}<span style={{ fontSize: 13, opacity: 0.7 }}>j </span></>}
                {routeResult.hours}<span style={{ fontSize: 13, opacity: 0.7 }}>h</span>
                {params.dt < 1 && routeResult.minutes > 0 && (
                  <>{String(routeResult.minutes).padStart(2, '0')}<span style={{ fontSize: 13, opacity: 0.7 }}>min</span></>
                )}
              </div>
              <div style={{ fontSize: 10, opacity: 0.4, marginTop: 6 }}>
                Calcul : {routeResult.calc_time_s}s
              </div>
            </div>
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
        const effectiveMeta = meta || currentMeta;
        const tMin = effectiveMeta.times[0];
        const tMax = effectiveMeta.times[effectiveMeta.times.length - 1];
        return (
          <div style={{
            position: 'absolute', bottom: 0, left: 0, right: 0, zIndex: 10,
            background: 'rgba(8,13,30,0.92)', backdropFilter: 'blur(12px)',
            borderTop: '1px solid rgba(100,160,255,0.15)',
            padding: '10px 24px 14px', fontFamily: FONT,
          }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', marginBottom: 6 }}>
              <strong style={{ color: '#4fc3f7', fontSize: 13 }}>{timeLabel}</strong>
              <span style={{ opacity: 0.4, fontSize: 11 }}>{timeOffset}</span>
            </div>
            <input
              type="range"
              min={tMin}
              max={tMax}
              step={meta ? STEP_OPTIONS[stepIdx].h : 1}
              value={currentTimeH}
              onChange={e => setCurrentTimeH(+e.target.value)}
              style={{ width: '100%', accentColor: '#4fc3f7', cursor: 'pointer', display: 'block' }}
            />
          </div>
        );
      })()}

      <div style={{
        position: 'absolute', bottom: 60, right: 16, zIndex: 10,
        background: 'rgba(8,13,30,0.7)', borderRadius: 8, padding: '6px 12px',
        color: 'rgba(150,170,220,0.55)', fontSize: 10, fontFamily: FONT,
      }}>
        Hover pour les valeurs
      </div>
    </div>
  );
}
