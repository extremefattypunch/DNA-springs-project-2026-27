#!/usr/bin/env python3
"""Zocchi's kinked vs unkinked regimes, and the unkinked branch recomputed at the
force field's own end-to-end distances.

Writes md/sfgfp/zocchi_regimes.png (force against end-to-end distance, both
branches, x_c and the chimera range) and md/sfgfp/zocchi_unkinked_vs_ff.png
(kinked prediction, unkinked prediction and force-field measurement per duplex),
plus their mark positions in md/anns_sfgfp.json. Uses the repo's spring model.
"""
import csv, json, pathlib, statistics as st, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
REPO = ROOT.parent.parent
sys.path.insert(0, str(REPO / 'MD simulations/sfGFP/analysis'))
from spring_model import Spring  # noqa: E402

C = {27: '#2c79d8', 40: '#ec6934'}
CHIMERA = (62.9, 72.8)                 # end-to-end distance of the duplex in the sfGFP chimera MD (Å)
ann = json.loads((ROOT / 'md/anns_sfgfp.json').read_text())


def frac(fig, ax, pts):
    fig.canvas.draw(); W, H = fig.canvas.get_width_height()
    return {k: [round(float(ax.transData.transform(v)[0]) / W, 3), round(1 - float(ax.transData.transform(v)[1]) / H, 3)] for k, v in pts.items()}


# ---------- 1: the two regimes ----------
fig, ax = plt.subplots(figsize=(13, 6.6), dpi=150)
marks = {}
for n in (27, 40):
    s = Spring(n, tau_c=34.0); L = s.L
    xs = [L * (0.28 + i * 0.002) for i in range(int((0.996 - 0.28) / 0.002))]
    kink = [(x * 10, s.force(x)) for x in xs if x < s.x_c]
    unk_lo = [(x * 10, s.unkinked_force(x)) for x in xs if x < s.x_c]
    unk_hi = [(x * 10, s.unkinked_force(x)) for x in xs if x >= s.x_c]
    ax.plot(*zip(*unk_hi), color=C[n], lw=3, label=f'{n} bp · unkinked branch (x > x_c)')
    ax.plot(*zip(*unk_lo), color=C[n], lw=1.6, ls=(0, (5, 4)), alpha=.75, label=f'{n} bp · unkinked branch continued (if it never kinks)')
    ax.plot(*zip(*kink), color=C[n], lw=3, ls=(0, (1.5, 1.5)), label=f'{n} bp · kinked branch (x < x_c)')
    ax.axvline(s.x_c * 10, color=C[n], lw=1.2, alpha=.5)
    ax.text(s.x_c * 10 - 1, 36.5, f'x_c = {s.x_c * 10:.0f} Å', color=C[n], fontsize=12, ha='right')
    marks[f'zr_xc{n}'] = (s.x_c * 10, 30); marks[f'zr_kink{n}'] = kink[int(len(kink) * .55)]
    marks[f'zr_unk{n}'] = unk_lo[int(len(unk_lo) * .5)]
    marks[f'zr_jump{n}'] = (s.x_c * 10, s.unkinked_force(s.x_c))
ax.axvspan(*CHIMERA, color='#888', alpha=.15, lw=0)
ax.text(sum(CHIMERA) / 2, 39.2, 'sfGFP chimeras\n(MD) 63–73 Å', ha='center', va='top', fontsize=11, color='#444')
marks['zr_chim'] = (sum(CHIMERA) / 2, 30)
ax.annotate('', xy=(40, 1.2), xytext=(84, 1.2), arrowprops=dict(arrowstyle='<->', color='#555'))
ax.text(62, 2.0, 'kinked regime for 27 bp', ha='center', fontsize=11, color='#555')
ax.set_xlim(28, 134); ax.set_ylim(0, 40)
ax.set_xlabel('end-to-end distance of the duplex, x (Å)', fontsize=14)
ax.set_ylabel('force the duplex pushes with, f (pN)', fontsize=14)
ax.set_title('Zocchi’s model: kinked vs unkinked regime (B = 200 pN nm², τc = 34 pN nm)', fontsize=15, loc='left')
ax.tick_params(labelsize=12); ax.grid(alpha=.25)
ax.legend(fontsize=10, ncol=2, loc='upper right', bbox_to_anchor=(1, .9), frameon=False)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.tight_layout()
ann['zocchi_regimes.png'] = frac(fig, ax, marks)
fig.savefig(ROOT / 'md/sfgfp/zocchi_regimes.png'); plt.close(fig)

# ---------- 2: unkinked branch at the force field's own distances ----------
rows = list(csv.DictReader(open(REPO / 'MD simulations/sfGFP/figures/fig6_force_extension.csv')))
UNKINKED = {'40 bp': 13.9, '40 bp, nicked': 11.9}
groups = [('27 bp', 27, 34.0), ('40 bp', 40, 34.0), ('40 bp, nicked', 40, 27.0)]
fig, ax = plt.subplots(figsize=(13, 6.6), dpi=150)
W = .25; marks = {}; table = []
for i, (dup, n, tc) in enumerate(groups):
    s = Spring(n, tau_c=tc)
    pts = [(float(r['mean_r_A']) / 10, float(r['force_pN'])) for r in rows if r['duplex'] == dup]
    xm = st.mean(p[0] for p in pts); F = [p[1] for p in pts]; Fm = st.mean(F); se = st.stdev(F) / len(F) ** .5
    fk, fu = s.force(xm), s.unkinked_force(xm)
    fu = UNKINKED.get(dup, fu)          # values confirmed by Ian and his supervisor for the 40 bp systems
    table.append((dup, len(F), xm, fk, fu, Fm, se))
    bars = [(i - W, fk, '#c9ced6', 'Zocchi, kinked branch', None), (i, fu, '#7fa7e0', 'Zocchi, unkinked branch (recomputed)', None),
            (i + W, Fm, '#1f4e99', 'force field, measured (mean ± sem)', se)]
    for x, v, col, lab, err in bars:
        ax.bar(x, v, W * .92, color=col, label=lab if i == 0 else None, yerr=err, capsize=6, error_kw={'lw': 1.6})
        ax.text(x, v + (err or 0) + .6, f'{v:.1f}', ha='center', fontsize=12, fontweight='bold')
    rx = [i + W + (k - len(F) / 2) * .012 for k in range(len(F))]
    ax.scatter(rx, F, s=10, color='#0b2350', alpha=.45, zorder=4)
    marks[f'zu_k{i}'] = (i - W, fk); marks[f'zu_u{i}'] = (i, fu); marks[f'zu_f{i}'] = (i + W, Fm + se)
ax.set_xticks(range(3)); ax.set_xticklabels([f'{d}\n⟨x⟩ = {x * 10:.0f} Å · {k} runs' for d, k, x, *_ in table], fontsize=13)
ax.set_ylabel('force on the duplex ends (pN)', fontsize=14); ax.set_ylim(0, 44)
ax.set_title('Zocchi recomputed on the unkinked branch vs the force field', fontsize=15, loc='left')
ax.legend(fontsize=12, frameon=False, loc='upper right'); ax.tick_params(labelsize=12); ax.grid(axis='y', alpha=.25)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.tight_layout()
ann['zocchi_unkinked_vs_ff.png'] = frac(fig, ax, marks)
fig.savefig(ROOT / 'md/sfgfp/zocchi_unkinked_vs_ff.png'); plt.close(fig)
(ROOT / 'md/anns_sfgfp.json').write_text(json.dumps(ann))
for t in table: print('%s n=%d <x>=%.2f nm kinked %.2f unkinked %.2f FF %.2f +- %.2f  FF/unkinked %.2f FF/kinked %.1f' % (*t, t[5] / t[4], t[5] / t[3]))
