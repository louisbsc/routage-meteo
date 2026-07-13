import DeckGL from '@deck.gl/react';
import { BitmapLayer, GeoJsonLayer, PathLayer, ScatterplotLayer, TextLayer } from '@deck.gl/layers';
import { useMemo } from 'react';
import { buildSpeedRaster } from './windfield';
import WindParticles from './WindParticles';

// Fond coloré + particules du courant : dégradé bleu clair → bleu foncé
// (intensité croissante), distinct du dégradé arc-en-ciel du vent.
const CURRENT_RASTER_STOPS = [
  [0,   [190, 225, 255]],
  [0.3, [120, 185, 255]],
  [0.8, [50,  130, 240]],
  [1.5, [15,  70,  200]],
  [2.5, [0,   25,  155]],
  [4.0, [0,   5,   100]],
];

function currentRasterColor(speed) {
  const s = Math.max(0, speed);
  for (let i = 1; i < CURRENT_RASTER_STOPS.length; i++) {
    const [s0, c0] = CURRENT_RASTER_STOPS[i - 1];
    const [s1, c1] = CURRENT_RASTER_STOPS[i];
    if (s <= s1) {
      const f = (s - s0) / (s1 - s0);
      return c0.map((v, j) => Math.round(v + f * (c1[j] - v)));
    }
  }
  return CURRENT_RASTER_STOPS[CURRENT_RASTER_STOPS.length - 1][1];
}

const FILL = { position: 'absolute', inset: 0, width: '100%', height: '100%' };

export default function WindMap({
  data, viewState, onViewStateChange,
  depPoint, arrPoint, route, boats = [],
  isochrones, showIsochrones, showGrib,
  currentData, showCurrent,
  clickMode, onMapClick,
  extraRoutes = [],
  landData = null,
  showParticles = false,
}) {
  const windRaster    = useMemo(() => buildSpeedRaster(data), [data]);
  const currentRaster = useMemo(() => buildSpeedRaster(currentData, 165, currentRasterColor), [currentData]);

  // Couche du bas : fond couleur vent + fond couleur courant
  const rasterLayer = useMemo(() => {
    const layers = [];
    if (showGrib && windRaster) layers.push(new BitmapLayer({
      id: 'wind-raster',
      image: windRaster.image,
      bounds: windRaster.bounds,
    }));
    if (showCurrent && currentRaster) layers.push(new BitmapLayer({
      id: 'current-raster',
      image: currentRaster.image,
      bounds: currentRaster.bounds,
    }));
    return layers;
  }, [showGrib, windRaster, showCurrent, currentRaster]);

  // Couches du dessus : terre, courant, isochrones, routes, marqueurs, bateaux
  const topLayers = useMemo(() => {
    const result = [];

    if (landData) result.push(new GeoJsonLayer({
      id: 'land',
      data: landData,
      filled: true,
      getFillColor: [18, 18, 18, 255],
      stroked: false,
    }));

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
  }, [landData, route, isochrones, showIsochrones, depPoint, arrPoint, boats, extraRoutes]);

  return (
    <div style={{ position: 'absolute', inset: 0 }}>
      {/* Fond + raster vent */}
      <DeckGL
        viewState={viewState}
        controller={false}
        layers={rasterLayer}
        style={{ ...FILL, background: '#3a3a3a' }}
      />

      {/* Particules : entre le raster et la terre/routes/marqueurs */}
      {showParticles && <WindParticles data={data} viewState={viewState} />}
      {showCurrent && <WindParticles data={currentData} viewState={viewState} speedScale={6} fixedSpeed={0} minSpeed={0.15} densityScale={2.5} />}

      {/* Terre, isochrones, routes, marqueurs (capte les événements) */}
      <DeckGL
        viewState={viewState}
        controller={{ dragPan: true, scrollZoom: true, doubleClickZoom: true, keyboard: false }}
        onViewStateChange={e => onViewStateChange(e.viewState)}
        layers={topLayers}
        onClick={e => {
          if (onMapClick && e.coordinate) {
            onMapClick([e.coordinate[1], e.coordinate[0]]);
          }
        }}
        getCursor={({ isDragging }) =>
          clickMode ? 'crosshair' : (isDragging ? 'grabbing' : 'grab')
        }
        style={{ ...FILL, background: 'transparent' }}
      />
    </div>
  );
}
