import DeckGL from '@deck.gl/react';
import { BitmapLayer, GeoJsonLayer, IconLayer, PathLayer, ScatterplotLayer, TextLayer } from '@deck.gl/layers';
import { useMemo } from 'react';
import { buildSpeedRaster } from './windfield';

const CURRENT_STOPS = [
  [0,   [190, 225, 255, 170]],
  [0.3, [120, 185, 255, 195]],
  [0.8, [50,  130, 240, 215]],
  [1.5, [15,  70,  200, 230]],
  [2.5, [0,   25,  155, 245]],
  [4.0, [0,   5,   100, 255]],
];

function currentSpeedColor(speed) {
  const s = Math.max(0, speed);
  for (let i = 1; i < CURRENT_STOPS.length; i++) {
    const [s0, c0] = CURRENT_STOPS[i - 1];
    const [s1, c1] = CURRENT_STOPS[i];
    if (s <= s1) {
      const f = (s - s0) / (s1 - s0);
      return c0.map((v, j) => Math.round(v + f * (c1[j] - v)));
    }
  }
  return CURRENT_STOPS[CURRENT_STOPS.length - 1][1];
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
  depPoint, arrPoint, route, boats = [],
  isochrones, showIsochrones, showGrib,
  currentData, showCurrent,
  clickMode, onMapClick,
  extraRoutes = [],
  landData = null,
}) {
  const windRaster = useMemo(() => buildSpeedRaster(data), [data]);

  const layers = useMemo(() => {
    const result = [];

    // ── Gradient vent (sous la terre pour frontière nette) ────────────────
    if (showGrib && windRaster) result.push(new BitmapLayer({
      id: 'wind-raster',
      image: windRaster.image,
      bounds: windRaster.bounds,
    }));

    // ── Terre (Natural Earth) ─────────────────────────────────────────────
    if (landData) result.push(new GeoJsonLayer({
      id: 'land',
      data: landData,
      filled: true,
      getFillColor: [18, 18, 18, 255],
      stroked: false,
    }));

    // ── Courant ───────────────────────────────────────────────────────────
    if (showCurrent && currentData?.length) result.push(new IconLayer({
      id: 'current-arrows',
      data: currentData.filter(d => d.speed > 0.05),
      iconAtlas, iconMapping,
      getIcon: () => 'arrow',
      getPosition: d => [d.lon, d.lat, 0],
      getSize:  d => Math.min(36, 20 + d.speed * 5),
      getAngle: d => -(d.dir + 180),
      getColor: d => currentSpeedColor(d.speed),
      pickable: true,
      billboard: true,
      updateTriggers: { getColor: currentData, getAngle: currentData, getPosition: currentData },
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

    // ── Routes sauvegardées actives ───────────────────────────────────────
    extraRoutes.forEach(saved => {
      if (saved.routeResult?.route?.length) {
        result.push(new PathLayer({
          id: `saved-route-${saved.id}`,
          data: [{ path: saved.routeResult.route }],
          getPath:  d => d.path,
          getColor: saved.color,
          getWidth: 3,
          widthUnits: 'pixels',
          widthMinPixels: 2,
          jointRounded: true,
          capRounded: true,
        }));
      }
    });

    // ── Route courante ────────────────────────────────────────────────────
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

    // ── Bateaux ───────────────────────────────────────────────────────────
    if (boats.length) {
      result.push(new ScatterplotLayer({
        id: 'boats',
        data: boats,
        getPosition: d => d.pos,
        getRadius: 8,
        radiusUnits: 'pixels',
        getFillColor: d => d.color,
        stroked: true,
        getLineColor: d => d.outline,
        lineWidthMinPixels: 2.5,
        pickable: false,
        updateTriggers: { getFillColor: boats, getLineColor: boats },
      }));
    }

    return result;
  }, [windRaster, route, isochrones, showIsochrones, showGrib, currentData, showCurrent, depPoint, arrPoint, boats, extraRoutes, landData]);

  return (
    <DeckGL
      viewState={viewState}
      controller={{ dragPan: true, scrollZoom: true, doubleClickZoom: true, keyboard: false }}
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
      getTooltip={({ object, layer }) =>
        object && {
          html: layer?.id === 'current-arrows'
            ? `<b>${object.speed.toFixed(2)} nœuds</b><br/>Courant : ${object.dir.toFixed(0)}°`
            : `<b>${object.speed.toFixed(1)} nœuds</b><br/>Vent : ${object.dir.toFixed(0)}°`,
          style: {
            background: 'rgba(10,12,28,0.92)', color: '#e0eaff',
            borderRadius: '6px', fontSize: '12px',
            padding: '6px 10px', fontFamily: 'monospace',
          },
        }
      }
      style={{ background: '#3a3a3a' }}
    />

  );
}
