import { useState, useEffect } from 'react';
import WindMap from './WindMap.jsx';

const API = 'http://localhost:8000';

const INIT_VIEW = { longitude: -5, latitude: 47, zoom: 4, pitch: 0, bearing: 0 };

const DAYS  = ['Dim', 'Lun', 'Mar', 'Mer', 'Jeu', 'Ven', 'Sam'];
const MONTHS = ['jan', 'fév', 'mar', 'avr', 'mai', 'juin', 'juil', 'aoû', 'sep', 'oct', 'nov', 'déc'];

function fmtDatetime(iso) {
  const d = new Date(iso);
  const day   = DAYS[d.getUTCDay()];
  const date  = String(d.getUTCDate()).padStart(2, '0');
  const month = MONTHS[d.getUTCMonth()];
  const hour  = String(d.getUTCHours()).padStart(2, '0');
  return `${day} ${date} ${month}  ${hour}h UTC`;
}

const LEGEND = [
  ['≤ 8 kt',  '#50a0ff'],
  ['16 kt',   '#32dc5a'],
  ['25 kt',   '#ffdc00'],
  ['35 kt',   '#ff8200'],
  ['≥ 50 kt', '#dc1e1e'],
];

/* ── Styles inline ── */
const S = {
  root: {
    width: '100vw', height: '100vh', position: 'relative', background: '#080d1a',
    fontFamily: "'SF Mono', 'Consolas', monospace",
  },
  panel: {
    position: 'absolute', top: 16, left: 16, zIndex: 10,
    background: 'rgba(8,13,30,0.92)', backdropFilter: 'blur(12px)',
    border: '1px solid rgba(100,160,255,0.15)',
    borderRadius: 12, padding: '16px 18px', minWidth: 270,
    color: '#c8d8ff', boxShadow: '0 8px 32px rgba(0,0,30,0.6)',
  },
  title: { fontSize: 11, letterSpacing: 2, opacity: 0.5, marginBottom: 14, textTransform: 'uppercase' },
  label: { fontSize: 11, opacity: 0.65, marginBottom: 5, display: 'block' },
  select: {
    width: '100%', padding: '7px 10px', borderRadius: 7, fontSize: 12,
    background: 'rgba(30,40,80,0.8)', color: '#c8d8ff',
    border: '1px solid rgba(100,160,255,0.2)', marginBottom: 16, outline: 'none',
  },
  range: { width: '100%', accentColor: '#4fc3f7', cursor: 'pointer' },
  sectionGap: { marginTop: 14 },
  spinner: { marginTop: 12, fontSize: 11, opacity: 0.45, textAlign: 'center' },
  timeBar: {
    position: 'absolute', bottom: 0, left: 0, right: 0, zIndex: 10,
    background: 'rgba(8,13,30,0.92)', backdropFilter: 'blur(12px)',
    borderTop: '1px solid rgba(100,160,255,0.15)',
    padding: '10px 24px 14px',
  },
  timeBarTop: {
    display: 'flex', justifyContent: 'space-between', alignItems: 'baseline',
    marginBottom: 6,
  },
  legend: {
    position: 'absolute', bottom: 72, left: 16, zIndex: 10,
    background: 'rgba(8,13,30,0.88)', backdropFilter: 'blur(8px)',
    border: '1px solid rgba(100,160,255,0.12)',
    borderRadius: 10, padding: '10px 14px',
    color: '#c8d8ff', fontSize: 11,
  },
  legendRow: { display: 'flex', alignItems: 'center', gap: 8, marginBottom: 5 },
  dot: { width: 12, height: 12, borderRadius: 3, flexShrink: 0 },
  badge: {
    position: 'absolute', bottom: 72, right: 16, zIndex: 10,
    background: 'rgba(8,13,30,0.7)', borderRadius: 8, padding: '6px 12px',
    color: 'rgba(150,170,220,0.6)', fontSize: 10,
  },
};

export default function App() {
  const [files, setFiles]           = useState([]);
  const [file, setFile]             = useState('');
  const [meta, setMeta]             = useState(null);
  const [timeIdx, setTimeIdx]       = useState(0);
  const [stride, setStride]         = useState(2);
  const [windData, setWindData]     = useState([]);
  const [loading, setLoading]       = useState(false);
  const [viewState, setViewState]   = useState(INIT_VIEW);

  // Load file list once
  useEffect(() => {
    fetch(`${API}/files`).then(r => r.json()).then(setFiles).catch(() => {});
  }, []);

  // Load metadata when file changes
  useEffect(() => {
    if (!file) return;
    setMeta(null);
    setWindData([]);
    fetch(`${API}/wind/${encodeURIComponent(file)}/meta`)
      .then(r => r.json())
      .then(m => {
        setMeta(m);
        setTimeIdx(0);
        const [lon0, lat0, lon1, lat1] = m.bbox;
        setViewState(v => ({
          ...v,
          longitude: (lon0 + lon1) / 2,
          latitude:  (lat0 + lat1) / 2,
          zoom: 4,
        }));
      })
      .catch(() => {});
  }, [file]);

  // Load wind data when time step or stride changes
  useEffect(() => {
    if (!file || !meta) return;
    setLoading(true);
    fetch(`${API}/wind/${encodeURIComponent(file)}/step/${timeIdx}?stride=${stride}`)
      .then(r => r.json())
      .then(({ data }) => { setWindData(data); setLoading(false); })
      .catch(() => setLoading(false));
  }, [file, meta, timeIdx, stride]);

  const timeLabel = meta?.valid_times?.[timeIdx]
    ? fmtDatetime(meta.valid_times[timeIdx])
    : '—';
  const timeOffset = meta?.times?.[timeIdx] !== undefined
    ? `T+${meta.times[timeIdx]}h`
    : '';

  const pointCount = windData.length;

  return (
    <div style={S.root}>
      <WindMap
        data={windData}
        viewState={viewState}
        onViewStateChange={setViewState}
      />

      {/* ── Controls ── */}
      <div style={S.panel}>
        <div style={S.title}>Wind GRIB Viewer</div>

        <label style={S.label}>Fichier GRIB</label>
        <select value={file} onChange={e => setFile(e.target.value)} style={S.select}>
          <option value="">— Choisir un fichier —</option>
          {files.map(f => <option key={f} value={f}>{f}</option>)}
        </select>

        {meta && (
          <>
            <div style={S.sectionGap}>
              <label style={S.label}>
                Densité&nbsp;
                <strong style={{ color: '#4fc3f7' }}>1 / {stride}</strong>
                <span style={{ opacity: 0.4, marginLeft: 6 }}>
                  ({pointCount.toLocaleString()} pts)
                </span>
              </label>
              <input
                type="range" min={1} max={6} value={stride}
                onChange={e => setStride(+e.target.value)}
                style={S.range}
              />
            </div>
          </>
        )}

        {loading && <div style={S.spinner}>Chargement…</div>}
      </div>

      {/* ── Legend ── */}
      <div style={S.legend}>
        {LEGEND.map(([label, color]) => (
          <div key={label} style={S.legendRow}>
            <div style={{ ...S.dot, background: color }} />
            {label}
          </div>
        ))}
      </div>

      {meta && (
        <div style={S.timeBar}>
          <div style={S.timeBarTop}>
            <strong style={{ color: '#4fc3f7', fontSize: 13 }}>{timeLabel}</strong>
            <span style={{ opacity: 0.4, fontSize: 11 }}>
              {timeOffset} · {timeIdx + 1}/{meta.times.length}
            </span>
          </div>
          <input
            type="range" min={0} max={meta.times.length - 1}
            value={timeIdx} onChange={e => setTimeIdx(+e.target.value)}
            style={{ ...S.range, display: 'block' }}
          />
        </div>
      )}

      <div style={S.badge}>Hover pour les valeurs</div>
    </div>
  );
}
