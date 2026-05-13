import DeckGL from '@deck.gl/react';
import { IconLayer, PathLayer, ScatterplotLayer, TextLayer } from '@deck.gl/layers';
import Map from 'react-map-gl/maplibre';
import { useMemo } from 'react';

const MAP_STYLE = 'https://basemaps.cartocdn.com/gl/dark-matter-nolabels-gl-style/style.json';

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

function buildIconAtlas() {
  const S = 32;
  const canvas = document.createElement('canvas');
  canvas.width = S; canvas.height = S;
  const ctx = canvas.getContext('2d');
  const cx = S / 2;
  ctx.fillStyle = 'white';
  ctx.beginPath();
  ctx.moveTo(cx, 2);
  ctx.lineTo(cx + 10, 15); ctx.lineTo(cx + 4, 13);
  ctx.lineTo(cx + 4,  30); ctx.lineTo(cx - 4, 30);
  ctx.lineTo(cx - 4,  13); ctx.lineTo(cx - 10, 15);
  ctx.closePath();
  ctx.fill();
  return canvas;
}

const iconAtlas   = buildIconAtlas();
const iconMapping = { arrow: { x: 0, y: 0, width: 32, height: 32, mask: true } };

export default function WindMap({
  data, viewState, onViewStateChange,
  depPoint, arrPoint, route, boatPosition,
  isochrones, showIsochrones, showGrib,
  clickMode, onMapClick,
}) {
  const layers = useMemo(() => {
    const result = [];

    // ── Vent ──────────────────────────────────────────────────────────────
    if (showGrib) result.push(new IconLayer({
      id: 'wind-arrows',
      data,
      iconAtlas, iconMapping,
      getIcon: () => 'arrow',
      getPosition: d => [d.lon, d.lat, 0],
      getSize:  d => Math.min(38, 18 + d.speed * 0.35),
      getAngle: d => -(d.dir + 180),
      getColor: d => speedColor(d.speed),
      pickable: true,
      billboard: true,
      updateTriggers: { getColor: data, getAngle: data, getPosition: data },
    }));

    // ── Isochrones ────────────────────────────────────────────────────────
    if (showIsochrones && isochrones?.length) {
      result.push(new PathLayer({
        id: 'isochrones',
        data: isochrones.map((path, i) => ({ path, i })),
        getPath:  d => d.path,
        getColor: [255, 255, 255, 55],
        getWidth: 1,
        widthUnits: 'pixels',
        widthMinPixels: 0.5,
      }));
    }

    // ── Route ─────────────────────────────────────────────────────────────
    if (route?.length) {
      result.push(new PathLayer({
        id: 'route-path',
        data: [{ path: route }],
        getPath:  d => d.path,
        getColor: [255, 255, 255, 230],
        getWidth: 3,
        widthUnits: 'pixels',
        widthMinPixels: 2,
        jointRounded: true,
        capRounded: true,
      }));
    }

    // ── Marqueurs départ / arrivée ────────────────────────────────────────
    const markers = [
      depPoint && { id: 'dep', pos: [depPoint[1], depPoint[0]], fill: [60, 220, 110, 255], label: 'DEP' },
      arrPoint && { id: 'arr', pos: [arrPoint[1], arrPoint[0]], fill: [230, 70,  70,  255], label: 'ARR' },
    ].filter(Boolean);

    if (markers.length) {
      result.push(
        new ScatterplotLayer({
          id: 'markers',
          data: markers,
          getPosition: d => d.pos,
          getRadius: 9,
          radiusUnits: 'pixels',
          getFillColor: d => d.fill,
          stroked: true,
          getLineColor: [255, 255, 255, 200],
          lineWidthMinPixels: 2,
          pickable: false,
        }),
        new TextLayer({
          id: 'marker-labels',
          data: markers,
          getPosition: d => d.pos,
          getText: d => d.label,
          getSize: 11,
          getColor: [255, 255, 255, 210],
          getPixelOffset: [0, -20],
          fontFamily: 'monospace',
          fontWeight: 'bold',
          background: true,
          getBackgroundColor: [10, 12, 28, 180],
          backgroundPadding: [3, 2],
        }),
      );
    }

    // ── Position du bateau ────────────────────────────────────────────────
    if (boatPosition) {
      result.push(new ScatterplotLayer({
        id: 'boat',
        data: [{ pos: boatPosition }],
        getPosition: d => d.pos,
        getRadius: 8,
        radiusUnits: 'pixels',
        getFillColor: [255, 255, 255, 255],
        stroked: true,
        getLineColor: [80, 160, 255, 255],
        lineWidthMinPixels: 2.5,
        pickable: false,
      }));
    }

    return result;
  }, [data, route, isochrones, showIsochrones, showGrib, depPoint, arrPoint, boatPosition]);

  return (
    <DeckGL
      viewState={viewState}
      controller={{ dragPan: true, scrollZoom: true, doubleClickZoom: true }}
      onViewStateChange={e => onViewStateChange(e.viewState)}
      layers={layers}
      onClick={e => {
        if (onMapClick && e.coordinate) {
          onMapClick([e.coordinate[1], e.coordinate[0]]); // [lat, lon]
        }
      }}
      getCursor={({ isDragging }) =>
        clickMode ? 'crosshair' : (isDragging ? 'grabbing' : 'grab')
      }
      getTooltip={({ object }) =>
        object && {
          html: `<b>${object.speed.toFixed(1)} nœuds</b><br/>Direction : ${object.dir.toFixed(0)}°`,
          style: {
            background: 'rgba(10,12,28,0.92)', color: '#e0eaff',
            borderRadius: '6px', fontSize: '12px',
            padding: '6px 10px', fontFamily: 'monospace',
          },
        }
      }
    >
      <Map mapStyle={MAP_STYLE} reuseMaps />
    </DeckGL>
  );
}
