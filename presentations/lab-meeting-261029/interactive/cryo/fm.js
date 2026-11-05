/* Shared image pipeline of the forward-model frames (forward.html, fmterms.html): a Gaussian-blob
   molecule, rotated by R_i and projected, then CTF, in-plane shift t_i and Gaussian noise, all on a
   64 px box, plus the display colour maps. Needs cryo.js (clamp). */
/* ---------------- the molecule: Gaussian blobs (pixel units of a 64 px box) ---------------- */
const BL = [];
for (let n = 0; n < 11; n++) { const a = (210 - n * 24) * Math.PI / 180; BL.push([13 * Math.cos(a), 13 * Math.sin(a) - 2, 4 * Math.sin(a * 1.3), 2.1, 1]); }
BL.push([-6, 9, 7, 4.3, .8], [8.5, 6.5, -6, 3.7, .8]);
function rotM(ax, th) { const [x, y, z] = ax, c = Math.cos(th), s = Math.sin(th), t = 1 - c;
  return [[t*x*x + c, t*x*y - s*z, t*x*z + s*y], [t*x*y + s*z, t*y*y + c, t*y*z - s*x], [t*x*z - s*y, t*y*z + s*x, t*z*z + c]]; }
const AX = (() => { const v = [0.35, 1, 0.25], n = Math.hypot(...v); return v.map(a => a / n); })(), TH = 72 * Math.PI / 180;
const mul = (M, p) => [0, 1, 2].map(r => M[r][0] * p[0] + M[r][1] * p[1] + M[r][2] * p[2]);
/* oblique camera: screen right, screen down = beam z, depth toward viewer */
const AZ = 0.6, EL = 0.45;
function cam(p) { const w = Math.sin(AZ) * p[0] + Math.cos(AZ) * p[1];
  return [Math.cos(AZ) * p[0] - Math.sin(AZ) * p[1], Math.cos(EL) * p[2] + Math.sin(EL) * w, Math.cos(EL) * w - Math.sin(EL) * p[2]]; }

/* ---------------- the images: rotate, project, CTF, shift, noise (computed once) ---------------- */
const N = 64, RI = rotM(AX, TH), TX = 7, TY = -5;
function fft(re, im, inv) { const n = re.length;
  for (let i = 1, j = 0; i < n; i++) { let b = n >> 1; for (; j & b; b >>= 1) j ^= b; j ^= b;
    if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; } }
  for (let len = 2; len <= n; len <<= 1) { const a = (inv ? 2 : -2) * Math.PI / len, wr = Math.cos(a), wi = Math.sin(a);
    for (let i = 0; i < n; i += len) { let cr = 1, ci = 0;
      for (let k = 0; k < len / 2; k++) { const ur = re[i + k], ui = im[i + k], vr = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci, vi = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
        re[i + k] = ur + vr; im[i + k] = ui + vi; re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi; [cr, ci] = [cr * wr - ci * wi, cr * wi + ci * wr]; } } }
  if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; } }
function fft2(re, im, inv) { const r = new Float64Array(N), i = new Float64Array(N);
  for (let y = 0; y < N; y++) { for (let x = 0; x < N; x++) { r[x] = re[y * N + x]; i[x] = im[y * N + x]; } fft(r, i, inv); for (let x = 0; x < N; x++) { re[y * N + x] = r[x]; im[y * N + x] = i[x]; } }
  for (let x = 0; x < N; x++) { for (let y = 0; y < N; y++) { r[y] = re[y * N + x]; i[y] = im[y * N + x]; } fft(r, i, inv); for (let y = 0; y < N; y++) { re[y * N + x] = r[y]; im[y * N + x] = i[y]; } } }
const fq = u => (u < N / 2 ? u : u - N) / N;
const proj = new Float64Array(N * N);
BL.forEach(b => { const [px, py] = mul(RI, b), s = b[3], A = b[4] * s;
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) { const X = x - N / 2 - px, Y = (N / 2 - y) - py; proj[y * N + x] += A * Math.exp(-(X * X + Y * Y) / (2 * s * s)); } });
const Fr = Float64Array.from(proj), Fi = new Float64Array(N * N); fft2(Fr, Fi, false);
const ctf = new Float64Array(N * N), ramp = new Float64Array(N * N);
for (let v = 0; v < N; v++) for (let u = 0; u < N; u++) { const q2 = fq(u) ** 2 + fq(v) ** 2, chi = Math.PI * 42 * q2 - Math.PI * 260 * q2 * q2;
  ctf[v * N + u] = -(0.995 * Math.sin(chi) + 0.1 * Math.cos(chi)); ramp[v * N + u] = -2 * Math.PI * (fq(u) * TX + fq(v) * TY); }
const Cr = Fr.map((a, k) => a * ctf[k]), Ci = Fi.map((a, k) => a * ctf[k]);
const Sr = Cr.map((a, k) => a * Math.cos(ramp[k]) - Ci[k] * Math.sin(ramp[k])), Si = Cr.map((a, k) => a * Math.sin(ramp[k]) + Ci[k] * Math.cos(ramp[k]));
const real = (re, im) => { const r = Float64Array.from(re), i = Float64Array.from(im); fft2(r, i, true); return r; };
const img3 = real(Cr, Ci), img4 = real(Sr, Si);
const sd = a => { const m = a.reduce((s, x) => s + x, 0) / a.length; return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / a.length); };
let seed = 7; const rnd = () => { seed = (seed + 0x6D2B79F5) | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t ^= t + Math.imul(t ^ t >>> 7, 61 | t); return ((t ^ t >>> 14) >>> 0) / 4294967296; };
const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
const s4 = sd(img4), img5 = img4.map(a => a + 0.9 * s4 * gauss());
const Yr = Float64Array.from(img5), Yi = new Float64Array(N * N); fft2(Yr, Yi, false);
const mag = (re, im) => re.map((a, k) => Math.hypot(a, im[k]));
const grey = v => { const g = Math.round(255 * clamp(v, 0, 1)); return [g, g, g]; };
function realMap(a, sym) { const m = sym ? 0 : 0, s = sd(a), mx = Math.max(...a);
  return v => grey(sym ? 0.5 + v / (6 * s) : v / mx); }
const shift = (a, k) => { const y = Math.floor(k / N), x = k % N; return a[((y + N / 2) % N) * N + (x + N / 2) % N]; };   // fftshift for display
function fourMap(a) { const L = a.map(v => Math.log(v + 1e-3)); const lo = Math.min(...L), hi = Math.max(...L);
  return (_, k) => { const v = (shift(L, k) - lo) / (hi - lo); return [Math.round(255 * clamp(v * 1.05, 0, 1) ** 1.2), Math.round(255 * clamp(v, 0, 1) ** 1.1), Math.round(255 * clamp(v * 1.1, 0, 1))]; }; }
