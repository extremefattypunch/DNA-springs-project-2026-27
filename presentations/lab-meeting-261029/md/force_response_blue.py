#!/usr/bin/env python3
"""Axial force clamp on the unclicked 2x Tet2-Et protein only (the blue series of
fig1_force_response), drawn for the RMSF slide. Writes md/sfgfp/fig1_clamp_only.png
and its mark positions into md/anns_sfgfp.json.
"""
import csv, json, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
REPO = ROOT.parent.parent
rows = [r for r in csv.DictReader(open(REPO / 'MD simulations/sfGFP/figures/fig1_force_response.csv')) if r['series'] == 'S6_clamp']
F = np.array([float(r['applied_pN']) for r in rows]); X = np.array([float(r['extension_A']) for r in rows]); E = np.array([float(r['sem_A']) for r in rows])
w = 1 / E ** 2
A = np.vstack([F, np.ones_like(F)]).T
cov = np.linalg.inv(A.T @ (A * w[:, None])); slope, icpt = cov @ (A.T @ (w * X)); se = np.sqrt(cov[0, 0])
print(f'weighted slope {slope:.4f} +- {se:.4f} A/pN; span {X.max() - X.min():.3f} A')

fig, ax = plt.subplots(figsize=(12.4, 6.6), dpi=150)
c = '#2c79d8'
ax.axhspan(X[0] - 0.5, X[0] + 0.5, color=c, alpha=.07, lw=0)
ax.text(19.8, X[0] + 0.46, '±0.5 Å around the unloaded distance', ha='right', va='top', fontsize=12, color='#555')
ax.errorbar(F, X, yerr=E, fmt='o-', color=c, ms=10, lw=2.2, capsize=5, mec='white', mew=1.2, zorder=5)
ax.set_xlim(-0.8, 21); ax.set_ylim(X[0] - 1.0, X[0] + 1.0)
ax.set_xlabel('force applied between the two attachment sites (pN)', fontsize=14)
ax.set_ylabel('Cβ–Cβ distance, Asp133–Asn149 (Å)', fontsize=14)
ax.set_title('Axial force clamp on the unclicked 2× Tet2-Et protein (n = 3 per force)', fontsize=15, loc='left')
ax.tick_params(labelsize=12); ax.grid(alpha=.25)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.canvas.draw()
W, H = fig.canvas.get_width_height()
marks = {'cl_0': (F[0], X[0]), 'cl_20': (F[-1], X[-1]), 'cl_4': (F[2], X[2]), 'cl_band': (10, X[0] + 0.5)}
anns = {k: [round(float(ax.transData.transform(v)[0]) / W, 3), round(1 - float(ax.transData.transform(v)[1]) / H, 3)] for k, v in marks.items()}
out = ROOT / 'md/sfgfp/fig1_clamp_only.png'; fig.savefig(out); print('wrote', out, anns)
a = json.loads((ROOT / 'md/anns_sfgfp.json').read_text()); a['fig1_clamp_only.png'] = anns
(ROOT / 'md/anns_sfgfp.json').write_text(json.dumps(a))
