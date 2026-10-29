#!/usr/bin/env python3
"""Per-step axis bend in the MD chimeras (fig4_dna_bend.csv) against what a kinked
duplex would look like: Zocchi's kinked rod at the same end-to-end distance, two straight
arms with the whole bend at one step. Writes md/sfgfp/fig4_dna_bend_kink.png and marks.
"""
import csv, json, math, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
REPO = ROOT.parent.parent
rows = list(csv.DictReader(open(REPO / 'MD simulations/sfGFP/figures/fig4_dna_bend.csv')))
SER = [('S3_spring27', '27 bp spring (MD, mean of 6)', '#1aa37a'), ('S4_spring40', '40 bp spring (MD, mean of 6)', '#e89c0c'),
       ('S5_spring40nick', '40 bp nicked (MD, mean of 6)', '#e8739d')]
def kink_bend(n_bp, x_nm):                       # Zocchi's kinked rod: bend at the apex between two straight arms
    half = 0.33 * n_bp / 2; h = math.sqrt(half ** 2 - (x_nm / 2) ** 2); return math.degrees(2 * math.atan2(h, x_nm / 2))
K = [(27, 6.8, 14, '#1aa37a'), (40, 6.6, 20, '#e89c0c')]   # n_bp, chimera end-to-end (nm), kink position (bp), colour
fig, ax = plt.subplots(figsize=(13, 6.8), dpi=150)
marks = {}
for key, lab, col in SER:
    pts = [(int(r['bp']), float(r[key])) for r in rows if r[key]]
    ax.plot(*zip(*pts), color=col, lw=2.6, label=lab, zorder=4)
    marks['dk_md_' + key] = max(pts, key=lambda p: p[1])
for n, x, pos, col in K:
    b = kink_bend(n, x); first = 6 if n == 27 else 7; last = 22 if n == 27 else 34
    xs = list(range(first, last + 1)); ys = [b if s == pos else 0 for s in xs]
    ax.step(xs, ys, where='mid', color=col, lw=2, ls=(0, (5, 3)), alpha=.9, zorder=3,
            label=f'{n} bp kinked (Zocchi’s rod at x = {x * 10:.0f} Å): {b:.0f}° at one step, 0° elsewhere')
    ax.text(pos + .7, b + 1, f'{b:.0f}°', color=col, fontsize=13, fontweight='bold', va='top')
    marks[f'dk_kink{n}'] = (pos, b)
ax.axhspan(0, 10, color='#888', alpha=.08, lw=0)
ax.text(34.6, 11, 'the whole MD range (0–10°)', ha='right', fontsize=11, color='#555')
marks['dk_band'] = (27, 6)
ax.set_xlim(5, 35); ax.set_ylim(0, 165)
ax.set_xlabel('base pair along the duplex', fontsize=14); ax.set_ylabel('axis bend at each base-pair step (°)', fontsize=14)
ax.set_title('MD duplexes bend smoothly: compare with a kinked duplex at the same end-to-end distance', fontsize=15, loc='left')
ax.legend(fontsize=11, frameon=False, loc='upper left', bbox_to_anchor=(0.0, 0.97)); ax.tick_params(labelsize=12); ax.grid(alpha=.25)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.canvas.draw(); W, H = fig.canvas.get_width_height()
anns = {k: [round(float(ax.transData.transform(v)[0]) / W, 3), round(1 - float(ax.transData.transform(v)[1]) / H, 3)] for k, v in marks.items()}
fig.savefig(ROOT / 'md/sfgfp/fig4_dna_bend_kink.png')
a = json.loads((ROOT / 'md/anns_sfgfp.json').read_text()); a['fig4_dna_bend_kink.png'] = anns
(ROOT / 'md/anns_sfgfp.json').write_text(json.dumps(a))
print({n: round(kink_bend(n, x), 1) for n, x, *_ in K}, anns)
