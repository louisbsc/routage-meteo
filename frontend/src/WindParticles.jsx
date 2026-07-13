// ─────────────────────────────────────────────────────────────────────────────
//  WindParticles.jsx — animation de particules « façon Windy ».
//  Canvas transparent superposé à la carte deck.gl : des milliers de particules
//  suivent le champ de vent, laissant des traînées colorées par la vitesse.
// ─────────────────────────────────────────────────────────────────────────────
import { useEffect, useRef } from 'react';
import { buildField, sample, makeProjector, speedColorCss } from './windfield';

// Réglages animation
const PARTICLE_DENSITY = 1 / 2000;  // particules par pixel²
const MAX_PARTICLES    = 5000;
const MAX_AGE          = 220;       // images avant ré-apparition
const TIME_STEP        = 0.025;     // « heures » de vent simulées par image (à zoom de référence)
const REF_WORLD_SIZE   = 512 * Math.pow(2, 5); // zoom 5 = référence vitesse
const MAX_STEP_PX      = 4;         // déplacement max par image (anti-saut)
const TRAIL_LENGTH     = 25;        // positions mémorisées par particule
const LINE_WIDTH       = 1.15;


export default function WindParticles({ data, viewState, speedScale = 1, fixedSpeed, minSpeed = 0, densityScale = 1 }) {
  const canvasRef  = useRef(null);
  const fieldRef   = useRef(null);
  const vsRef      = useRef(viewState);
  const projRef    = useRef(null);    // projecteur courant (mis à jour chaque frame)
  const tooltipRef = useRef(null);   // div tooltip (mise à jour imperative, sans re-render)

  useEffect(() => { fieldRef.current = buildField(data); }, [data]);
  useEffect(() => { vsRef.current = viewState; }, [viewState]);

  // ── Animation loop ────────────────────────────────────────────────────────
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const dpr = Math.min(2, window.devicePixelRatio || 1);

    let width = 0, height = 0;
    const resize = () => {
      const r = canvas.getBoundingClientRect();
      width = r.width; height = r.height;
      canvas.width = Math.round(width * dpr);
      canvas.height = Math.round(height * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    };
    resize();
    const ro = new ResizeObserver(resize);
    ro.observe(canvas);

    const particles = [];
    const reseed = (p, proj) => {
      const px = Math.random() * width;
      const py = Math.random() * height;
      const [lng, lat] = proj.unproject(px, py);
      p.lng = lng; p.lat = lat;
      p.age = Math.floor(Math.random() * MAX_AGE);
      p.trail = [];   // efface la traînée à la mort
    };

    let raf = 0;
    const loop = () => {
      const field = fieldRef.current;
      const proj  = makeProjector(vsRef.current, width, height);
      projRef.current = proj;   // expose le projecteur au tooltip

      const target = Math.min(MAX_PARTICLES * densityScale, Math.floor(width * height * PARTICLE_DENSITY * densityScale));
      while (particles.length < target) { const p = {}; reseed(p, proj); particles.push(p); }
      if (particles.length > target) particles.length = target;

      // Efface le canvas entier chaque frame : pas de traînée orpheline
      ctx.clearRect(0, 0, width, height);

      if (!field) {
        raf = requestAnimationFrame(loop);
        return;
      }

      ctx.lineWidth = LINE_WIDTH;
      ctx.lineCap   = 'round';

      for (const p of particles) {
        if (p.age >= MAX_AGE) { reseed(p, proj); }
        const w = sample(field, p.lng, p.lat);
        if (!w) { reseed(p, proj); continue; }

        const latR     = p.lat * Math.PI / 180;
        const geoScale = REF_WORLD_SIZE / proj.worldSize;
        const step     = TIME_STEP * speedScale;
        const nLat     = p.lat + w.v / 60 * step * geoScale;
        const nLng     = p.lng + w.u / 60 * step * geoScale / Math.max(0.2, Math.cos(latR));

        const [x0, y0] = proj.project(p.lng, p.lat);
        let   [x1, y1] = proj.project(nLng, nLat);

        const dx = x1 - x0, dy = y1 - y0;
        const d  = Math.hypot(dx, dy);
        if (d > MAX_STEP_PX) { const k = MAX_STEP_PX / d; x1 = x0 + dx * k; y1 = y0 + dy * k; }

        if (x1 < -30 || x1 > width + 30 || y1 < -30 || y1 > height + 30) { reseed(p, proj); continue; }

        // Enregistre la position géographique dans la traînée (pas les coords écran)
        p.trail.push({ lng: p.lng, lat: p.lat, speed: w.speed });
        if (p.trail.length > TRAIL_LENGTH) p.trail.shift();

        // Dessine la traînée : re-projection à chaque frame → stable au zoom/pan
        // (segments sous minSpeed omis : pas de particules visibles en zone de courant trop faible)
        const n = p.trail.length;
        for (let i = 0; i < n - 1; i++) {
          if (p.trail[i].speed < minSpeed) continue;
          const [tx0, ty0] = proj.project(p.trail[i].lng, p.trail[i].lat);
          const [tx1, ty1] = proj.project(p.trail[i + 1].lng, p.trail[i + 1].lat);
          ctx.strokeStyle = speedColorCss(fixedSpeed ?? p.trail[i].speed, (i + 1) / n * 0.88);
          ctx.beginPath();
          ctx.moveTo(tx0, ty0);
          ctx.lineTo(tx1, ty1);
          ctx.stroke();
        }
        // Segment de tête (pleine opacité)
        if (w.speed >= minSpeed) {
          ctx.strokeStyle = speedColorCss(fixedSpeed ?? w.speed, 0.92);
          ctx.beginPath();
          ctx.moveTo(x0, y0);
          ctx.lineTo(x1, y1);
          ctx.stroke();
        }

        p.lng = nLng; p.lat = nLat; p.age++;
      }

      raf = requestAnimationFrame(loop);
    };
    raf = requestAnimationFrame(loop);

    return () => { cancelAnimationFrame(raf); ro.disconnect(); };
  }, []);

  // ── Tooltip vent au survol ─────────────────────────────────────────────────
  useEffect(() => {
    const onMove = (e) => {
      const canvas  = canvasRef.current;
      const tooltip = tooltipRef.current;
      const field   = fieldRef.current;
      const proj    = projRef.current;
      if (!canvas || !tooltip || !field || !proj) return;

      const rect = canvas.getBoundingClientRect();
      const px = e.clientX - rect.left;
      const py = e.clientY - rect.top;

      if (px < 0 || px > rect.width || py < 0 || py > rect.height) {
        tooltip.style.display = 'none';
        return;
      }

      const [lng, lat] = proj.unproject(px, py);
      const w = sample(field, lng, lat);

      if (!w) {
        tooltip.style.display = 'none';
        return;
      }

      // direction d'où vient le vent (convention météo)
      const dir = Math.round((Math.atan2(w.u, w.v) * 180 / Math.PI + 180 + 360) % 360);

      // Position du tooltip : à droite du curseur, avec clamp pour rester à l'écran
      let tx = e.clientX + 18;
      let ty = e.clientY - 30;
      const W = window.innerWidth, H = window.innerHeight;
      if (tx + 145 > W) tx = e.clientX - 155;
      if (ty + 70  > H) ty = H - 78;
      if (ty < 8)        ty = 8;

      tooltip.style.display = 'block';
      tooltip.style.left    = tx + 'px';
      tooltip.style.top     = ty + 'px';
      tooltip.innerHTML = `
        <div style="font-size:11px;color:#c8d8ff;opacity:0.85;font-variant-numeric:tabular-nums">
          ${w.speed.toFixed(1)} kt &thinsp; ${dir}°
        </div>`;
    };

    const onLeave = () => {
      if (tooltipRef.current) tooltipRef.current.style.display = 'none';
    };

    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseleave', onLeave);
    return () => {
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('mouseleave', onLeave);
    };
  }, []);

  return (
    <>
      <canvas
        ref={canvasRef}
        style={{
          position: 'absolute', inset: 0, width: '100%', height: '100%',
          pointerEvents: 'none',
        }}
      />
      <div
        ref={tooltipRef}
        style={{
          display: 'none',
          position: 'fixed', zIndex: 20,
          background: 'rgba(4,6,16,0.97)', backdropFilter: 'blur(10px)',
          border: '1px solid rgba(100,160,255,0.10)',
          borderRadius: 6, padding: '4px 8px',
          fontFamily: "'SF Mono','Consolas',monospace",
          pointerEvents: 'none',
          whiteSpace: 'nowrap',
          boxShadow: '0 4px 20px rgba(0,0,30,0.5)',
        }}
      />
    </>
  );
}
