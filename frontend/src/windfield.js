// ─────────────────────────────────────────────────────────────────────────────
//  windfield.js — palette, projection Web-Mercator (compatible deck.gl) et
//  champ de vent interpolé, partagés par le fond coloré et l'animation Windy.
// ─────────────────────────────────────────────────────────────────────────────

// Palette de vitesse de vent (nœuds) — même esprit que Windy.
const STOPS = [
  [0,  [ 70, 140, 235]],
  [8,  [ 60, 200, 110]],
  [16, [235, 215,  40]],
  [25, [240, 130,  25]],
  [35, [225,  45,  45]],
  [50, [150,  20, 150]],
  [70, [120,  20,  90]],
];

function lerpColor(speed) {
  const s = Math.max(0, speed);
  for (let i = 1; i < STOPS.length; i++) {
    const [s0, c0] = STOPS[i - 1];
    const [s1, c1] = STOPS[i];
    if (s <= s1) {
      const f = (s - s0) / (s1 - s0);
      return [
        Math.round(c0[0] + f * (c1[0] - c0[0])),
        Math.round(c0[1] + f * (c1[1] - c0[1])),
        Math.round(c0[2] + f * (c1[2] - c0[2])),
      ];
    }
  }
  return STOPS[STOPS.length - 1][1];
}

export function speedColorArr(speed, alpha = 255) {
  const [r, g, b] = lerpColor(speed);
  return [r, g, b, alpha];
}

export function speedColorCss(speed, alpha = 0.9) {
  const [r, g, b] = lerpColor(speed);
  return `rgba(${r},${g},${b},${alpha})`;
}

// ── Projection Web-Mercator identique à la MapView de deck.gl (TILE = 512) ────
const TILE = 512;
const D2R = Math.PI / 180;

function worldX(lng, worldSize) {
  return (lng + 180) / 360 * worldSize;
}
function worldY(lat, worldSize) {
  const phi = lat * D2R;
  return (Math.PI - Math.log(Math.tan(Math.PI / 4 + phi / 2))) / (2 * Math.PI) * worldSize;
}

export function makeProjector(viewState, width, height) {
  const worldSize = TILE * Math.pow(2, viewState.zoom);
  const cx = worldX(viewState.longitude, worldSize);
  const cy = worldY(viewState.latitude, worldSize);
  const ox = width / 2 - cx;
  const oy = height / 2 - cy;
  return {
    worldSize,
    project(lng, lat) {
      return [worldX(lng, worldSize) + ox, worldY(lat, worldSize) + oy];
    },
    unproject(px, py) {
      const x = px - ox;
      const y = py - oy;
      const lng = x / worldSize * 360 - 180;
      const n = Math.PI - 2 * Math.PI * (y / worldSize);
      const lat = Math.atan(Math.sinh(n)) / D2R;
      return [lng, lat];
    },
  };
}

// ── Champ de vent régulier (lat/lon) construit depuis le tableau de points ───
//    Chaque point : { lon, lat, speed (nœuds), dir (°, direction d'où vient le
//    vent) }. Les cases manquantes (terre) sont marquées absentes.
export function buildField(data) {
  if (!data || data.length < 4) return null;

  const lonSet = Array.from(new Set(data.map(d => d.lon))).sort((a, b) => a - b);
  const latSet = Array.from(new Set(data.map(d => d.lat))).sort((a, b) => a - b);
  const nLon = lonSet.length, nLat = latSet.length;
  if (nLon < 2 || nLat < 2) return null;

  const lon0 = lonSet[0], lat0 = latSet[0];
  const dLon = (lonSet[nLon - 1] - lon0) / (nLon - 1);
  const dLat = (latSet[nLat - 1] - lat0) / (nLat - 1);

  const N = nLon * nLat;
  const U = new Float32Array(N);   // composante est (nœuds)
  const V = new Float32Array(N);   // composante nord (nœuds)
  const SP = new Float32Array(N);  // vitesse (nœuds)
  const M = new Uint8Array(N);     // 1 = donnée présente

  for (const d of data) {
    const i = Math.round((d.lon - lon0) / dLon);
    const j = Math.round((d.lat - lat0) / dLat);
    if (i < 0 || i >= nLon || j < 0 || j >= nLat) continue;
    const idx = j * nLon + i;
    // dir = direction d'OÙ vient le vent → il souffle VERS dir + 180
    const bearing = (d.dir + 180) * D2R;
    U[idx] = d.speed * Math.sin(bearing);
    V[idx] = d.speed * Math.cos(bearing);
    SP[idx] = d.speed;
    M[idx] = 1;
  }

  return { lon0, lat0, dLon, dLat, nLon, nLat, U, V, SP, M };
}

// Échantillonnage bilinéaire (en ignorant les coins absents près des côtes).
export function sample(field, lng, lat) {
  const { lon0, lat0, dLon, dLat, nLon, nLat, U, V, SP, M } = field;
  const fx = (lng - lon0) / dLon;
  const fy = (lat - lat0) / dLat;
  if (fx < -0.5 || fy < -0.5 || fx > nLon - 0.5 || fy > nLat - 0.5) return null;

  let i0 = Math.floor(fx), j0 = Math.floor(fy);
  i0 = Math.max(0, Math.min(nLon - 2, i0));
  j0 = Math.max(0, Math.min(nLat - 2, j0));
  const tx = fx - i0, ty = fy - j0;

  const corners = [
    [i0,     j0,     (1 - tx) * (1 - ty)],
    [i0 + 1, j0,     tx * (1 - ty)],
    [i0,     j0 + 1, (1 - tx) * ty],
    [i0 + 1, j0 + 1, tx * ty],
  ];

  let u = 0, v = 0, sp = 0, w = 0;
  for (const [i, j, wgt] of corners) {
    const idx = j * nLon + i;
    if (!M[idx]) continue;
    u += U[idx] * wgt;
    v += V[idx] * wgt;
    sp += SP[idx] * wgt;
    w += wgt;
  }
  if (w < 1e-4) return null;
  return { u: u / w, v: v / w, speed: sp / w };
}

// ── Fond coloré (raster) : une image lissée positionnée sous la couche terre ─
// Les cellules sans donnée (terre côté GRIB) sont remplies par BFS depuis la
// mer la plus proche, ce qui évite les trous et donne une frontière nette une
// fois que la couche GeoJSON de la terre recouvre le raster.
export function buildSpeedRaster(data, alpha = 185) {
  const field = buildField(data);
  if (!field) return null;
  const { lon0, lat0, dLon, dLat, nLon, nLat, SP, M } = field;
  const N = nLon * nLat;

  // BFS : propage les vitesses mer vers les cellules manquantes (terre GRIB)
  const SP_filled = new Float32Array(SP);
  const M_filled  = new Uint8Array(M);
  const queue = [];
  for (let k = 0; k < N; k++) { if (M[k]) queue.push(k); }
  let qi = 0;
  while (qi < queue.length) {
    const idx = queue[qi++];
    const j = Math.floor(idx / nLon);
    const i = idx % nLon;
    const nb = [
      i > 0       ? idx - 1    : -1,
      i < nLon-1  ? idx + 1    : -1,
      j > 0       ? idx - nLon : -1,
      j < nLat-1  ? idx + nLon : -1,
    ];
    for (const n of nb) {
      if (n >= 0 && !M_filled[n]) {
        SP_filled[n] = SP_filled[idx];
        M_filled[n]  = 1;
        queue.push(n);
      }
    }
  }

  // Upsampling 4× par interpolation bilinéaire : gradient lisse sans extra-requête API
  const SCALE = Math.min(8, Math.max(1, Math.floor(2048 / Math.max(nLon, nLat))));
  const W = nLon * SCALE;
  const H = nLat * SCALE;

  const cv  = document.createElement('canvas');
  cv.width  = W;
  cv.height = H;
  const ctx = cv.getContext('2d');
  const img = ctx.createImageData(W, H);

  // BitmapLayer interpole les UV linéairement en espace écran = en Y Mercator.
  // Le canvas doit donc échantillonner les latitudes via l'inverse Mercator pour
  // que les données apparaissent à la bonne position géographique.
  const lat_max = lat0 + (nLat - 1) * dLat;
  const D2R_half = Math.PI / 360;
  const mercYmax = Math.log(Math.tan(Math.PI / 4 + lat_max * D2R_half));
  const mercYmin = Math.log(Math.tan(Math.PI / 4 + lat0   * D2R_half));
  const mercSpan = mercYmax - mercYmin;

  for (let r = 0; r < H; r++) {   // r=0 = haut du canvas = Nord
    // Latitude géographique correspondant à cette ligne canvas (via Mercator inverse)
    const t   = r / (H - 1);
    const mercY   = mercYmax - t * mercSpan;
    const lat_r   = (2 * Math.atan(Math.exp(mercY)) - Math.PI / 2) * 180 / Math.PI;
    const fy  = (lat_r - lat0) / dLat;
    if (fy < 0 || fy > nLat - 1) continue;

    const j0 = Math.max(0, Math.min(nLat - 2, Math.floor(fy)));
    const ty = fy - j0;

    for (let pi = 0; pi < W; pi++) {
      const fx = (pi / (W - 1)) * (nLon - 1);
      const i0 = Math.max(0, Math.min(nLon - 2, Math.floor(fx)));
      const tx = fx - i0;

      // Interpolation bilinéaire sur SP_filled (toutes cellules remplies après BFS)
      const sp =
        SP_filled[ j0      * nLon +  i0     ] * (1 - tx) * (1 - ty) +
        SP_filled[ j0      * nLon + (i0 + 1)] *      tx  * (1 - ty) +
        SP_filled[(j0 + 1) * nLon +  i0     ] * (1 - tx) *      ty  +
        SP_filled[(j0 + 1) * nLon + (i0 + 1)] *      tx  *      ty;

      const dst = (r * W + pi) * 4;
      if (M_filled[j0 * nLon + i0]) {
        const [rc, g, b] = lerpColor(sp);
        img.data[dst]     = rc;
        img.data[dst + 1] = g;
        img.data[dst + 2] = b;
        img.data[dst + 3] = alpha;
      }
    }
  }
  ctx.putImageData(img, 0, 0);

  const bounds = [
    lon0 - dLon / 2,
    lat0 - dLat / 2,
    lon0 + (nLon - 1) * dLon + dLon / 2,
    lat0 + (nLat - 1) * dLat + dLat / 2,
  ];
  return { image: cv, bounds };
}
