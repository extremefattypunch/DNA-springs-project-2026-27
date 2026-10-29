#!/usr/bin/env python3
"""Static excitation + emission spectra of the 1-tet concentration series (2026-08-25),
for the band-ratio slide. Data: the PAYLOAD embedded in interactive/dilution.html
(copied from docs/fluorometry/dilution-260825.html). Means of the three scans per
concentration, linear y, no normalisation. Excitation is S1 / R1 (lamp-corrected),
emission is S1. Writes charts/band_ratio_spectra.png and its mark positions to
md/anns_sfgfp.json under that name.
"""
import json, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
s = (ROOT / 'interactive/dilution.html').read_text()
i = s.index('PAYLOAD'); i += s[i:].index('{')
P, _ = json.JSONDecoder().raw_decode(s[i:])
T = [t for t in P['traces'] if t['group'] != 'calibration']
GROUPS = [g for g in P['groups'] if g != 'calibration']
cols = plt.cm.Blues(np.linspace(.38, .98, len(GROUPS)))

fig, ax = plt.subplots(figsize=(13.5, 7.4), dpi=150)
ax2 = ax.twinx()
ax2.axvspan(398, 500, color='#dbe8fb', alpha=.55, lw=0, zorder=0)
ax2.axvspan(500, 602, color='#dff3e4', alpha=.6, lw=0, zorder=0)
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
marks, rows = {}, []
for g, c in zip(GROUPS, cols):
    ex = np.mean([np.array(t['s1'], float) / np.array(t['r1'], float) for t in T if t['group'] == g and t['mode'] == 'excitation'], 0)
    em = np.mean([np.array(t['s1'], float) for t in T if t['group'] == g and t['mode'] == 'emission'], 0)
    xe = 400 + np.arange(len(ex)); xm = 500 + np.arange(len(em))
    ax.plot(xe, ex, color=c, lw=2.6)
    ax2.plot(xm, em, color=c, lw=2.6, ls=(0, (6, 2.5)))
    k = int(ex.argmax()); pk, i400 = ex[k], ex[0]
    ax.plot([xe[k]], [pk], 'o', ms=8, mfc='white', mec=c, mew=2.2, zorder=6)
    ax.plot([400], [i400], 's', ms=8, mfc='white', mec=c, mew=2.2, zorder=6)
    j = int(em.argmax()); ax2.text(xm[j] + 2.5, em[j], g, color=c, fontsize=11, fontweight='bold', va='bottom')
    ax.text(xe[k] - 1.5, pk, g, color=c, fontsize=11, fontweight='bold', ha='right', va='bottom')
    rows.append((g, i400, pk, xe[k], i400 / pk, xm[int(em.argmax())]))
    gk = g.replace(' µM', '')
    marks[f'br_pk{gk}'] = (ax, xe[k], pk); marks[f'br_400_{gk}'] = (ax, 400, i400); marks[f'br_em{gk}'] = (ax2, xm[int(em.argmax())], em.max())
ax.set_xlim(395, 605); ax.set_yscale('log'); ax2.set_yscale('log'); ax.set_ylim(6e4, 1.5e9); ax2.set_ylim(3e3, 1.5e9 / 15)
ax.set_xlabel('wavelength (nm): excitation wavelength swept (left band), emission wavelength swept (right band)', fontsize=13)
ax.set_ylabel('excitation intensity, S1 ÷ R1 (CPS / µA)', fontsize=13, color='#1f4e99')
ax2.set_ylabel('emission intensity, S1 (CPS)', fontsize=13, color='#1d6b3a')
ax.tick_params(labelsize=11); ax2.tick_params(labelsize=11)
ax.text(402, 1.3e9, 'EXCITATION spectra (solid)\ndetect at 510 nm, sweep excitation 400 → 500 nm', fontsize=12, fontweight='bold', color='#1f4e99', va='top')
ax.text(503, 1.3e9, 'EMISSION spectra (dashed)\nexcite at 485 nm, sweep emission 500 → 600 nm', fontsize=12, fontweight='bold', color='#1d6b3a', va='top')
tab = 'band ratio = I(400) / I(peak)\n' + '\n'.join(f'{g:>6}:  {r:.3f}' for g, _, _, _, r, _ in rows)
ax.text(556, 3.2e8, tab, fontsize=11.5, family='monospace', va='top', bbox=dict(fc='white', ec='#9aa3ae', boxstyle='round,pad=.5'))
ax.plot([], [], 's', ms=8, mfc='white', mec='#333', label='I(400): excitation at 400 nm')
ax.plot([], [], 'o', ms=8, mfc='white', mec='#333', label='I(peak): excitation maximum')
ax.legend(loc='upper left', bbox_to_anchor=(0.33, 0.86), fontsize=11, frameon=True)
ax.set_title('1-tet sfGFP concentration series, 2026-08-25 · mean of 3 scans per concentration · log y-axis, not normalised', fontsize=14, loc='left')
for a in (ax, ax2):
    for sp in ('top',): a.spines[sp].set_visible(False)
ax.grid(alpha=.25)
# inset: peak intensity vs concentration, log-log, against strict proportionality
conc = np.array([float(g.split()[0]) for g, *_ in rows]); pks = np.array([pk for _, _, pk, *_ in rows])
ins = ax.inset_axes([0.09, 0.56, 0.22, 0.27])
ins.loglog(conc, pks, 'o-', color='#1f4e99', lw=2, ms=6, label='measured I(peak)')
ins.loglog(conc, pks[0] * conc / conc[0], '--', color='#e5484d', lw=1.6, label='if fluorescence ∝ concentration')
ins.set_xticks(conc); ins.set_xticklabels([f'{c:g}' for c in conc], fontsize=9); ins.minorticks_off()
ins.tick_params(labelsize=9); ins.set_xlabel('concentration (µM)', fontsize=9); ins.set_ylabel('I(peak)', fontsize=9)
ins.legend(fontsize=8, frameon=False, loc='upper left'); ins.set_title('16× protein → only 11× signal', fontsize=10, fontweight='bold')
fig.tight_layout(); fig.canvas.draw()
W, H = fig.canvas.get_width_height()
anns = {k: [round(float(a.transData.transform((x, y))[0]) / W, 3), round(1 - float(a.transData.transform((x, y))[1]) / H, 3)] for k, (a, x, y) in marks.items()}
out = ROOT / 'charts/band_ratio_spectra.png'; fig.savefig(out)
A = json.loads((ROOT / 'md/anns_sfgfp.json').read_text()); A['band_ratio_spectra.png'] = anns
(ROOT / 'md/anns_sfgfp.json').write_text(json.dumps(A))
for r in rows: print('%s I400 %.3g Ipeak %.3g at %d ratio %.4f em max %d' % r)
