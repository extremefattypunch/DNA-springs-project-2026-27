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
import sys
sys.path.insert(0, str(REPO / 'MD simulations/sfGFP/analysis'))
from spring_model import Spring  # noqa: E402
def kink_bend(n_bp, x_nm):
    """Bend at the kink, Zocchi Molecular Machines Fig. 2.31: internal angle pi - 2*theta, so bend = 2*theta,
    theta = arccos(x/2R) - Delta (eq. 2.177), Delta = arccos(x_c/2R) (eq. 2.178), R and x_c from eqs 2.197, 2.194."""
    sp = Spring(n_bp, tau_c=34.0)
    return math.degrees(2 * (math.acos(x_nm / (2 * sp.R)) - math.acos(sp.x_c / (2 * sp.R))))
K = [(27, 6.8, 14, '#1aa37a'), (40, 6.6, 20, '#e89c0c')]   # n_bp, chimera end-to-end (nm), kink position (bp), colour
def draw(kinks):
    fig, ax = plt.subplots(figsize=(13, 6.8), dpi=150)
    marks = {}
    for key, lab, col in SER:
        pts = [(int(r['bp']), float(r[key])) for r in rows if r[key]]
        ax.plot(*zip(*pts), color=col, lw=2.6, label=lab, zorder=4)
        marks['dk_md_' + key] = max(pts, key=lambda p: p[1])
    if kinks:
        for n, x, pos, col in K:
            b = kink_bend(n, x); first = 6 if n == 27 else 7; last = 22 if n == 27 else 34
            xs = list(range(first, last + 1)); ys = [b if q == pos else 0 for q in xs]
            ax.step(xs, ys, where='mid', color=col, lw=2, ls=(0, (5, 3)), alpha=.9, zorder=3,
                    label=f'{n} bp kinked (Zocchi, x = {x * 10:.0f} Å): {b:.0f}° at one step, 0° elsewhere')
            ax.text(pos + .7, b + 1, f'{b:.0f}°', color=col, fontsize=13, fontweight='bold', va='top')
            marks[f'dk_kink{n}'] = (pos, b)
        ax.axhspan(0, 10, color='#888', alpha=.08, lw=0)
        ax.text(34.6, 11, 'the whole MD range (0–10°)', ha='right', fontsize=11, color='#555')
        marks['dk_band'] = (27, 6)
        ax.set_ylim(0, 100)
        ax.set_title('Same MD profiles, now with a kinked duplex at the same end-to-end distance', fontsize=15, loc='left')
        ax.legend(fontsize=10.5, frameon=False, loc='upper left', bbox_to_anchor=(0.0, 0.97))
    else:
        ax.set_ylim(0, 10)
        ax.set_title('Axis bend per base-pair step in the MD chimeras (mean of 6 runs each)', fontsize=15, loc='left')
        ax.legend(fontsize=11, frameon=False, loc='upper right')
    ax.set_xlim(5, 35)
    ax.set_xlabel('base pair along the duplex', fontsize=14); ax.set_ylabel('axis bend at each base-pair step (°)', fontsize=14)
    ax.tick_params(labelsize=12); ax.grid(alpha=.25)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.canvas.draw(); W, H = fig.canvas.get_width_height()
    pre = 'k_' if kinks else 'm_'
    anns = {pre + k: [round(float(ax.transData.transform(v)[0]) / W, 3), round(1 - float(ax.transData.transform(v)[1]) / H, 3)] for k, v in marks.items()}
    name = 'fig4_dna_bend_kink.png' if kinks else 'fig4_dna_bend_md.png'
    fig.savefig(ROOT / 'md/sfgfp' / name); plt.close(fig)
    a = json.loads((ROOT / 'md/anns_sfgfp.json').read_text()); a[name] = anns
    (ROOT / 'md/anns_sfgfp.json').write_text(json.dumps(a))


draw(False); draw(True)
print({n: round(kink_bend(n, x), 1) for n, x, *_ in K})
