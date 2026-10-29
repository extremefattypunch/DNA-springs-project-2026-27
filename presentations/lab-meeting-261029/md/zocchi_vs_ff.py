#!/usr/bin/env python3
"""Zocchi's two force branches against the force field's bare-duplex measurements.

Writes md/sfgfp/zocchi_vs_ff.png and the mark positions (image fractions) into
md/anns_sfgfp.json under that file name. Uses the repo's own spring model.
"""
import csv, json, math, pathlib, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
REPO = ROOT.parent.parent
sys.path.insert(0, str(REPO / 'MD simulations/sfGFP/analysis'))
from spring_model import Spring  # noqa: E402

C = {27: '#2c79d8', 40: '#ec6934'}
rows = list(csv.DictReader(open(REPO / 'MD simulations/sfGFP/figures/fig6_force_extension.csv')))
fig, ax = plt.subplots(figsize=(13, 6.2), dpi=150)
marks = {}
for n in (27, 40):
    s = Spring(n, tau_c=34.0)
    L = s.L
    xs = [L * (0.30 + i * 0.0068) for i in range(int((0.995 - 0.30) / 0.0068))]
    bent = [(x * 10, s.unkinked_force(x)) for x in xs]
    kinked = [(x * 10, s.force(x)) for x in xs if x < s.x_c]
    euler = math.pi ** 2 * s.B / L ** 2
    ax.plot(*zip(*bent), color=C[n], lw=2.6, label=f'{n} bp · bent (unkinked) branch, eq. 2.203')
    ax.plot(*zip(*kinked), color=C[n], lw=2.6, ls=(0, (2, 2)), label=f'{n} bp · kinked branch, Zocchi τc = 34 pN nm')
    ax.axhline(euler, color=C[n], lw=1.2, ls=':', alpha=.9)
    ax.text(30, euler + 0.6, f'Euler load π²B/L² = {euler:.1f} pN', color=C[n], fontsize=11)
    ax.axvline(s.x_c * 10, color=C[n], lw=1, alpha=.35)
    ax.text(s.x_c * 10 + 0.6, 33, f'x_c, {n} bp', color=C[n], fontsize=10, rotation=90, va='top')
    pts = [(float(r['mean_r_A']), float(r['force_pN'])) for r in rows if r['duplex'] == f'{n} bp']
    ax.scatter(*zip(*pts), s=34, color=C[n], edgecolor='white', lw=.8, zorder=5, label=f'{n} bp · measured, force field')
    mx = sum(p[0] for p in pts) / len(pts); my = sum(p[1] for p in pts) / len(pts)
    marks[f'zf_meas{n}'] = (mx, my); marks[f'zf_kink{n}'] = (kinked[len(kinked) // 2])
    marks[f'zf_bent{n}'] = bent[int(len(bent) * .55)]; marks[f'zf_xc{n}'] = (s.x_c * 10, 22)
ax.axvspan(62.9, 72.8, color='#888', alpha=.13, lw=0)
ax.text(67.8, -1.6, 'sfGFP chimeras\n63–73 Å', ha='center', va='top', fontsize=10, color='#555')
marks['zf_chimera'] = (67.8, 4)
ax.set_xlim(25, 135); ax.set_ylim(-5, 44)
ax.set_xlabel('end-to-end distance of the duplex, x (Å)', fontsize=13)
ax.set_ylabel('force the duplex pushes with, f (pN)', fontsize=13)
ax.set_title('Zocchi’s two branches against the force field', fontsize=15, loc='left')
ax.grid(alpha=.25); ax.legend(fontsize=9.5, ncol=2, loc='upper right', frameon=False)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.canvas.draw()
W, H = fig.canvas.get_width_height()
anns = {k: [round(float(ax.transData.transform(v)[0]) / W, 3), round(1 - float(ax.transData.transform(v)[1]) / H, 3)] for k, v in marks.items()}
out = ROOT / 'md/sfgfp/zocchi_vs_ff.png'
fig.savefig(out); print('wrote', out, anns)
a = json.loads((ROOT / 'md/anns_sfgfp.json').read_text()); a['zocchi_vs_ff.png'] = anns
(ROOT / 'md/anns_sfgfp.json').write_text(json.dumps(a))
for n in (27, 40):
    s = Spring(n, tau_c=34.0)
    print(n, 'L', round(s.L, 2), 'gamma', round(s.gamma, 3), 'x_c A', round(s.x_c * 10, 1), '10B/L2', round(10 * s.B / s.L ** 2, 2),
          'pi2B/L2', round(math.pi ** 2 * s.B / s.L ** 2, 2), 'curv tau/B deg/bp', round(math.degrees(s.tau_c / s.B * 0.33), 2))
