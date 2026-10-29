/* Shared runtime for the cryoEM explainer frames.
   Each page defines TARGETS (one object of numeric params per deck step) and
   draw(P, t). The deck posts {deck: 'cryo', step} on every click; params ease
   toward the step's targets so stepping back and forth always lands cleanly.
   Colours follow the Manim clips (helpers/style.py). */
const C = { bg: '#0E1116', panel: '#1A1F27', grid: '#3A424D', obs: '#4FD1E8', truth: '#5BD67A',
  est: '#F5A142', alt: '#B48CF2', red: '#F2645A', text: '#ECEDEF', muted: '#9AA3AE', label: '#FFD166' };
const NS = 'http://www.w3.org/2000/svg';
function S(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs || {}) {
    if (k === 'text') e.textContent = attrs[k];
    else if (k === 'html') e.innerHTML = attrs[k];
    else e.setAttribute(k, attrs[k]);
  }
  if (parent) parent.appendChild(e);
  return e;
}
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const lerp = (a, b, u) => a + (b - a) * u;

const Cryo = {
  step: 0, P: {}, targets: [], tau: 0.45, static: false, paused: false, clock: 0, draw: null, onStep: null, t0: performance.now(),
  init(targets, draw, opts) {
    this.targets = targets; this.draw = draw; Object.assign(this, opts || {});
    this.P = Object.assign({}, targets[0]);
    let last = performance.now();
    const loop = now => {
      const dt = this.paused ? 0 : Math.min(0.1, (now - last) / 1000); last = now; this.clock += dt;
      const T = this.targets[Math.min(this.step, this.targets.length - 1)];
      const u = this.static ? 1 : this.paused ? 0 : 1 - Math.exp(-dt / this.tau);
      for (const k in T) {
        if (typeof T[k] !== 'number') { this.P[k] = T[k]; continue; }
        const c = this.P[k] == null ? T[k] : this.P[k];
        const r = this.rates && this.rates[k];   // constant-speed params (units per second)
        if (r && !this.static) this.P[k] = c + clamp(T[k] - c, -r * dt, r * dt);
        else this.P[k] = Math.abs(T[k] - c) < 1e-4 ? T[k] : lerp(c, T[k], u);
      }
      this.draw(this.P, this.clock, dt);
      requestAnimationFrame(loop);
    };
    requestAnimationFrame(loop);
  }
};

addEventListener('message', e => {
  const m = e.data || {};
  if (m.deck && m.deck !== 'cryo') return;
  if ('pause' in m) { Cryo.paused = !!m.pause; return; }   // deck's Pause/Play button
  if (m.static) Cryo.static = true;
  const s = m.step || 0;
  if (s !== Cryo.step) { const old = Cryo.step; Cryo.step = s; if (Cryo.onStep) Cryo.onStep(s, old); }
});
/* keep the deck's keys working after a click inside the frame */
addEventListener('keydown', e => {
  if (["ArrowRight", "ArrowLeft", "ArrowUp", "ArrowDown", "PageDown", "PageUp", " ", "Enter", "Escape", "Backspace", "f", "p"].includes(e.key) || /^[0-9]$/.test(e.key)) {
    e.preventDefault(); parent.postMessage({ deckKey: e.key }, '*');
  }
});
/* a click in the frame advances the deck like a click on the slide */
addEventListener('click', e => { if (!e.target.closest('input, button, a')) parent.postMessage({ deckKey: e.clientX < innerWidth * 0.15 ? 'ArrowLeft' : 'ArrowRight' }, '*'); });
/* render checks: page.html?step=N&static */
{ const q = new URLSearchParams(location.search); if (q.has('step')) { Cryo.step = +q.get('step'); Cryo.static = q.has('static'); } }
