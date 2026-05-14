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

  // Routing
  const [polaires, setPolaires]     = useState([]);
  const [polaire, setPolaire]       = useState('');
  const [depPoint, setDepPoint]     = useState(null);  // [lat, lon]
  const [arrPoint, setArrPoint]     = useState(null);
  const [clickMode, setClickMode]   = useState(null);  // 'dep' | 'arr' | null
  const [routeResult, setRouteResult] = useState(null);
  const [routing, setRouting]       = useState(false);
  const [routeError, setRouteError] = useState(null);
  const [depTimeIdx, setDepTimeIdx] = useState(0);    // index dans meta.times
  const [showIsochrones, setShowIsochrones] = useState(true);
  const [advOpen, setAdvOpen]       = useState(false);
  const [params, setParams]         = useState({
    dt: 1, n: 100, ang_deg: 90, dang_deg: 0.3, delta: 2,
  });

  // ── Init ──────────────────────────────────────────────────────────────
  useEffect(() => {
    fetch(`${API}/files`).then(r => r.json()).then(setFiles).catch(() => {});
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
    setRouting(true); setRouteResult(null); setRouteError(null);
    try {
      const body = {
        polaire_file: polaire,
        p_dep: depPoint, p_arr: arrPoint,
        t: windMode === 'grib' ? (meta?.times?.[depTimeIdx] ?? 0) : 0,
        ...params,
        ...(windMode === 'uniform'
          ? { wind_uniform: { direction: uniformWind.direction, force: uniformWind.force } }
          : { grib_file: file }),
      };
      const res = await fetch(`${API}/routing`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(body),
      });
      if (!res.ok) throw new Error((await res.json()).detail ?? res.statusText);
      setRouteResult(await res.json());
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

  // ── Computed ──────────────────────────────────────────────────────────
  const timeLabel = useMemo(() => {
    if (!meta) return '—';
    // Datetime de référence (premier pas GRIB) + décalage courant
    const ref = new Date(meta.valid_times[0]);
    const delta = currentTimeH - meta.times[0];
    ref.setUTCMinutes(ref.getUTCMinutes() + Math.round(delta * 60));
    return fmtDatetime(ref.toISOString());
  }, [meta, currentTimeH]);

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
        clickMode={clickMode}
        onMapClick={handleMapClick}
      />

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

        {/* ── Carte polaire ─────────────────────────────────────────── */}
        <div style={card}>
          <div style={{ fontSize: 11, letterSpacing: 2, opacity: 0.4, textTransform: 'uppercase', marginBottom: 12 }}>
            Polaire
          </div>
          <label style={labelStyle}>Fichier polaire</label>
          <select value={polaire} onChange={e => setPolaire(e.target.value)} style={inputStyle}>
            {polaires.map(p => <option key={p} value={p}>{p.replace('.csv', '')}</option>)}
          </select>
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
                ['delta',        'delta',    1,    10,  0.5 ],
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
              marginBottom: (routeResult || routing || routeError) ? 10 : 0,
            }}>
            {routing ? 'Calcul en cours…' : 'Lancer le routage'}
          </button>

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
      {meta && (
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
            min={meta.times[0]}
            max={meta.times[meta.times.length - 1]}
            step={STEP_OPTIONS[stepIdx].h}
            value={currentTimeH}
            onChange={e => setCurrentTimeH(+e.target.value)}
            style={{ width: '100%', accentColor: '#4fc3f7', cursor: 'pointer', display: 'block' }}
          />
        </div>
      )}

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
