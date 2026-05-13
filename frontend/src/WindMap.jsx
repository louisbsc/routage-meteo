import DeckGL from '@deck.gl/react';
import { IconLayer } from '@deck.gl/layers';
import Map from 'react-map-gl/maplibre';
import { useMemo } from 'react';

const MAP_STYLE = 'https://basemaps.cartocdn.com/gl/dark-matter-nolabels-gl-style/style.json';

// Wind speed (knots) → RGBA color
const STOPS = [
  [0,  [80,  160, 255, 210]],
  [8,  [50,  220, 90,  210]],
  [16, [255, 220, 0,   220]],
  [25, [255, 130, 0,   230]],
  [35, [220, 30,  30,  240]],
  [50, [140, 0,   0,   255]],
];

function speedColor(speed) {
  const s = Math.max(0, speed);
  for (let i = 1; i < STOPS.length; i++) {
    const [s0, c0] = STOPS[i - 1];
    const [s1, c1] = STOPS[i];
    if (s <= s1) {
      const f = (s - s0) / (s1 - s0);
      return c0.map((v, j) => Math.round(v + f * (c1[j] - v)));
    }
  }
  return STOPS[STOPS.length - 1][1];
}

// Arrow icon canvas: white upward-pointing arrow (used as mask → tinted by getColor)
function buildIconAtlas() {
  const S = 32;
  const canvas = document.createElement('canvas');
  canvas.width = S;
  canvas.height = S;
  const ctx = canvas.getContext('2d');
  const cx = S / 2;
  ctx.fillStyle = 'white';
  // Arrowhead
  ctx.beginPath();
  ctx.moveTo(cx, 2);
  ctx.lineTo(cx + 10, 15);
  ctx.lineTo(cx + 4,  13);
  ctx.lineTo(cx + 4,  30);
  ctx.lineTo(cx - 4,  30);
  ctx.lineTo(cx - 4,  13);
  ctx.lineTo(cx - 10, 15);
  ctx.closePath();
  ctx.fill();
  return canvas;
}

const iconAtlas = buildIconAtlas();
const iconMapping = { arrow: { x: 0, y: 0, width: 32, height: 32, mask: true } };

export default function WindMap({ data, viewState, onViewStateChange }) {
  const layers = useMemo(() => [
    new IconLayer({
      id: 'wind-arrows',
      data,
      iconAtlas,
      iconMapping,
      getIcon: () => 'arrow',
      getPosition: d => [d.lon, d.lat, 0],
      // Pixels: base 18 + slight scaling with speed
      getSize: d => Math.min(38, 18 + d.speed * 0.35),
      // dir = FROM direction (météo, CW depuis Nord) → arrow points TO
      getAngle: d => -(d.dir + 180),
      getColor: d => speedColor(d.speed),
      pickable: true,
      billboard: true,
      updateTriggers: { getColor: data, getAngle: data, getPosition: data },
    }),
  ], [data]);

  return (
    <DeckGL
      viewState={viewState}
      controller={{ dragPan: true, scrollZoom: true, doubleClickZoom: true }}
      onViewStateChange={e => onViewStateChange(e.viewState)}
      layers={layers}
      getTooltip={({ object }) =>
        object && {
          html: `<b>${object.speed.toFixed(1)} nœuds</b><br/>Direction : ${object.dir.toFixed(0)}°`,
          style: {
            background: 'rgba(10,12,28,0.92)',
            color: '#e0eaff',
            borderRadius: '6px',
            fontSize: '12px',
            padding: '6px 10px',
            fontFamily: 'monospace',
          },
        }
      }
    >
      <Map mapStyle={MAP_STYLE} reuseMaps />
    </DeckGL>
  );
}
