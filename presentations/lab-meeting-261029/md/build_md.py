#!/usr/bin/env python3
"""Write the MD analysis slides (explain phase + answers) into deck.json.

Each figure is shown full-slide first with its axes and terms explained and its
points of interest marked, then the slide shrinks to the figure | answers view.
Mark positions are fractions of the image, in md/anns_*.json.
"""
import json, pathlib
from PIL import Image
ROOT = pathlib.Path(__file__).resolve().parent.parent
ANN = {**json.loads((ROOT / 'md/anns_sfgfp.json').read_text()), **json.loads((ROOT / 'md/anns_pp2a.json').read_text())}

def fig(d, name, cap):
    w, h = Image.open(ROOT / 'md' / d / name).size
    return {'img': f'md/{d}/{name}', 'w': w, 'h': h, 'cap': cap, 'anns': ANN.get(name, {})}

def M(a, t='', dx=None, dy=None):
    m = {'a': a, 't': t}
    if dx is not None: m['dx'] = dx
    if dy is not None: m['dy'] = dy
    return m

RAMA = ('<dfn>Ramachandran plot</dfn>: every dot is one simulation frame, placed by the two backbone dihedral angles of a single residue. '
        '<dfn>φ (phi)</dfn>, x-axis: rotation about the N–Cα bond. <dfn>ψ (psi)</dfn>, y-axis: rotation about the Cα–C bond. Both run −180° to +180°. '
        'The boxes are the allowed regions: β / PPII (extended strand), right-handed α (helix) and left-handed α (rare; mostly Gly, Asp, Asn). '
        'Large symbols are each system’s mean. If a modification strains the backbone, its cloud leaves the WT basin.')
RMSF = ('<dfn>RMSF</dfn> (root-mean-square fluctuation): how far each residue’s Cα atom wanders from its own average position over the trajectory, '
        'after every frame is aligned on the rigid core. Low = rigid, high = mobile. ')

S, P = 'sfgfp', 'pp2a'
slides = {
'12b': {'title': 'sfGFP: Ramachandran plots at the modified sites', 'q': 'Does adding the tet groups change the protein’s local conformation?',
  'items': [fig(S, 'ramachandran.png', 'Backbone φ/ψ at the two modified sites')],
  'explain': [
    {'item': 0, 'text': RAMA, 'marks': [M('rama_lha133', 'left-handed α box', 40, -70)]},
    {'item': 0, 'text': 'Left panel: <b>site 133</b> (Asp133). Right panel: <b>site 149</b> (Asn149). Colours: <b>blue</b> wild type, <b>orange</b> 2× Tet2-Et unclicked, <b>green</b> clicked adduct. The two sites sit on adjacent β-strands with His148, which bonds the chromophore, between them.',
     'marks': [M('rama_wt133', 'WT', -30, -60), M('rama_tet133', 'Tet, unclicked', 30, 14), M('rama_clk133', 'clicked', 30, -60), M('rama_wt149', 'Asn149 cluster', 30, -50)]}],
  'points': [
    {'text': '<b>Asp133</b> sits in the rare left-handed α basin: WT +57, +16', 'hl': 0, 'marks': [M('rama_wt133', 'WT', -30, -60), M('rama_lha133', '')]},
    {'text': 'Unclicked 2× Tet2-Et pushes it to +74, +2, about <b>17.5° out</b> of the basin', 'hl': 0, 'marks': [M('rama_tet133', 'Tet', 30, 14)]},
    {'text': 'After the click it returns to +60, +18, <b>back beside WT</b>', 'hl': 0, 'marks': [M('rama_clk133', 'clicked', 30, -60)]},
    {'text': '<b>Asn149</b> stays in β: −127, +137 (WT), −125, +124 (tet), −114, +133 (clicked)', 'hl': 0, 'marks': [M('rama_wt149', 'Asn149', 30, -50)]},
    {'text': '<b>Answer</b>: no significant change; only a local, reversible shift at the one positive-φ site'}]},
'12c': {'title': 'sfGFP: backbone rigidity (RMSF)', 'q': 'Does the spring deform the scaffold? Does the fold stay intact?',
  'items': [fig(S, 'fig2_rmsf.png', 'Cα RMSF, aligned on the barrel core'), fig(S, 'fig1_force_response.png', 'Axial force clamp between the two Cβ atoms')],
  'explain': [
    {'item': 0, 'text': RMSF + 'x-axis: residue number (2B3P numbering); y-axis: Cα RMSF (Å). Dashed lines mark the chromophore (CRO 66) and the attachment sites Asp133, His148, Asn149. One colour per system: WT, 2× Tet2-Et and the 27 / 40 / 40-nicked bp springs (n = replicates).',
     'marks': [M('rmsf_barrel', 'β-strand: rigid', 30, -60), M('rmsf_loop', 'loop: mobile', 30, 10), M('rmsf_cro', 'chromophore', 30, -50), M('rmsf_sites', 'attachment sites', 30, -50)]},
    {'item': 1, 'text': '<dfn>Axial force clamp</dfn>: a constant force applied directly between the two Cβ atoms, a controlled stand-in for a spring. x-axis: applied force (pN); y-axis: Cβ–Cβ distance Asp133–Asn149 (Å). <b>Blue</b>: clamp on the unclicked protein; <b>green</b>: clamp on the clicked tether; <b>orange</b>: the DNA-spring chimeras at their modelled force. <dfn>Compliance</dfn> is the slope (Å per pN); the <dfn>stiffness κ</dfn> is its inverse.',
     'marks': [M('fr_clamp20', '20 pN clamp', -30, -60), M('fr_dip4', '4 pN: noise', 30, 10), M('fr_27', '27 bp spring', 30, -50), M('fr_text', 'compliance fit', -30, -60)]}],
  'points': [
    {'text': 'The β-barrel stays <b>uniformly rigid</b> (≈ 0.4 Å on the strands) in every system', 'hl': 0, 'marks': [M('rmsf_barrel', 'rigid strand', 30, -60)]},
    {'text': 'RMSF traces <b>overlie each other</b>; only loops and tethers move', 'hl': 0, 'marks': [M('rmsf_loop', 'loop', 30, 10)]},
    {'text': 'Clamp: 0–20 pN moves the sites by hundredths of an Å; <b>κ > 107 kBT nm⁻²</b>', 'hl': 1, 'marks': [M('fr_clamp20', '20 pN', -30, -60), M('fr_text', '')]},
    {'text': 'On the clicked topology <b>κ ≈ 350 kBT nm⁻²</b>, so a 4 pN spring moves the sites only 0.03 Å', 'hl': 1},
    {'text': '<b>Answer</b>: the scaffold is not deformed; structural integrity persists'}]},
'12d': {'title': 'sfGFP: forces applied, force field vs model', 'q': 'How much force reaches the protein, and do the models agree with the force field?',
  'items': [fig(S, 'fig6_force_extension.png', 'Force measured on the bare duplex, restraint windows'), fig(S, 'fig7_headline.png', 'Site displacement against modelled force')],
  'explain': [
    {'item': 0, 'text': 'A bare duplex, no protein, held at a series of end-to-end separations by a harmonic restraint between its two linker nitrogens. x-axis: the separation the duplex actually adopted (Å); y-axis: the force the duplex exerts on the restraint, ⟨k(r − x₀)⟩ (pN). Points are runs, lines are means. Dashed: the <dfn>Euler buckling load</dfn> π²B/L², the force at which a straight elastic rod buckles. Grey band: what <dfn>Zocchi’s kinked-rod model</dfn> predicts for the same duplexes.',
     'marks': [M('fe_27', '27 bp: 28 pN', 30, -60), M('fe_euler27', 'Euler load, 27 bp', 30, 14), M('fe_40', '40 bp: 15.8 pN', 30, 10), M('fe_zocchi', 'Zocchi model', 30, -55)]},
    {'item': 1, 'text': 'x-axis: the force the spring model says each DNA delivers (pN). y-axis: the shift of the attachment sites against each construct’s own clicked control (Å); below zero = pulled together. <b>Blue</b> sfGFP, <b>orange</b> ubiquitin (preprint). Ringed points carry no spring at all (slack strand or arms that cannot close). A real force effect would climb towards the upper right.',
     'marks': [M('hl_zero', 'no change', 30, 14), M('hl_ringed', 'no spring at all', 30, -55), M('hl_27two', '27 bp two-arm', -30, 14), M('hl_stats', 'r = +0.08', -30, -60)]}],
  'points': [
    {'text': 'Measured directly: <b>27 bp 28.0 pN</b>, 40 bp 15.8 pN, 40 bp nicked 13.2 pN', 'hl': 0, 'marks': [M('fe_27', '28 pN', 30, -60), M('fe_40', '15.8 pN', 30, 10)]},
    {'text': 'That is ≈ the <b>Euler buckling load</b> (27 bp: 26.8 pN); up to 35 pN in the largest box', 'hl': 0, 'marks': [M('fe_euler27', 'Euler', 30, 14)]},
    {'text': 'Zocchi’s kinked branch predicts only <b>2.4–6.6 pN</b>: the force field gives ≈ 6× more', 'hl': 0, 'marks': [M('fe_zocchi', 'model', 30, -55)]},
    {'text': 'Every construct pulls its sites slightly <b>together</b>; displacement vs force r = +0.08, p = 0.81', 'hl': 1, 'marks': [M('hl_stats', 'r = +0.08', -30, -60), M('hl_zero', '')]},
    {'text': 'The force-free ssDNA loop compresses the sites by 0.32 Å, as much as any duplex', 'hl': 1, 'marks': [M('hl_ringed', 'no spring', 30, -55)]},
    {'text': '<b>Answer</b>: tens of pN on paper, but the stiff barrel does not move; the model underestimates the force field ≈ 6×'}]},
'12e': {'title': 'sfGFP: DNA dynamics', 'q': 'Is the DNA kinked or coiled? Does it interact with the protein beyond pushing?',
  'items': [fig(S, 'fig4_dna_bend.png', 'Axis bend per base-pair step along each duplex'), fig(S, 'fig5_kink_threshold.png', 'Base-pair opening against imposed curvature')],
  'explain': [
    {'item': 0, 'text': 'x-axis: base pair along the duplex in the chimera; y-axis: <dfn>axis bend</dfn> at that base-pair step (°), the angle between the local helix axes either side of it. Line: mean across replicates; band: ± sd. A <dfn>kink</dfn> is a sharp local break of the duplex, usually with a base pair opened: it would show as one narrow spike. Smooth bending spreads the angle evenly.',
     'marks': [M('db_27peak', 'highest step, 27 bp', 30, -60), M('db_40', '40 bp: even bending', 30, 14)]},
    {'item': 1, 'text': 'Bare 27 and 40 bp duplexes held at set curvatures. x-axis: curvature (° per base-pair step); y-axis: fraction of runs in which an <b>interior</b> base pair opens. Ticks are single runs (top = opened, bottom = stayed paired). <b>Blue band</b>: where our chimera sits; <b>grey band</b>: the type II kink range Curuksu et al. found by umbrella sampling.',
     'marks': [M('kt_chimera', 'our chimera 4.8–5.2°/bp', 30, -60), M('kt_curuksu', 'Curuksu 8.6–9.3°/bp', -30, -60), M('kt_first', 'first opening 5.5°/bp', 30, 14), M('kt_last', 'always opens', -30, 14)]}],
  'points': [
    {'text': 'The duplex bends <b>smoothly</b> at 4.8–5.2°/bp, nucleosome-like, with no local spike', 'hl': 0, 'marks': [M('db_40', 'even', 30, 14)]},
    {'text': 'Only 0.03–0.83 base pairs open at any instant: the duplex stays <b>intact</b>', 'hl': 0},
    {'text': 'Pairs start opening at 5.5°/bp and always open above <b>8.7°/bp</b>', 'hl': 1, 'marks': [M('kt_first', '5.5', 30, 14), M('kt_last', '8.7+', -30, 14)]},
    {'text': 'So the chimera is <b>not bent enough to kink</b> on this force field', 'hl': 1, 'marks': [M('kt_chimera', 'chimera', 30, -60)]},
    {'text': 'The DNA touches the protein in 72–98 % of frames (50 nt ssDNA loop 93 %): <b>local contact</b>, not action at a distance'},
    {'text': '<b>Answer</b>: no kink, no coiling of the duplex; but the DNA lies against the protein'}]},
'31c': {'title': 'PP2A: Ramachandran plots at the modified sites', 'q': 'Does adding the tet groups change the protein’s local conformation?',
  'items': [fig(P, 'pp2a_ramachandran.png', 'Backbone φ/ψ at the two modified sites')],
  'explain': [
    {'item': 0, 'text': RAMA + ' The star is the AlphaFold starting model.', 'marks': [M('pr_lha7', 'left-handed α box', 40, -70), M('pr_af7', 'AlphaFold start', 30, 14)]},
    {'item': 0, 'text': 'Left: <b>Asp7</b> (construct D6X), on the flexible N-terminal tail. Right: <b>Asn553</b> (N552X), in the first turn of HEAT 14’s helix. Colours: <b>blue</b> apo Aα, <b>orange</b> 2× Tet2-Et unclicked, <b>dark green</b> clicked no DNA, <b>teal</b> 50 nt ssDNA loop, <b>amber</b> 50 bp dsDNA loop.',
     'marks': [M('pr_apo7', 'apo', -30, -60), M('pr_tet7', 'Tet, unclicked', 30, 14), M('pr_clk7', 'clicked', -30, -60), M('pr_ss7', 'ssDNA', 30, 14), M('pr_ds7', 'dsDNA', -30, 14), M('pr_553', 'Asn553', 30, -50)]}],
  'points': [
    {'text': '<b>Asp7</b> starts in left-handed α (+65, +57), on the floppy tail (pLDDT 52)', 'hl': 0, 'marks': [M('pr_af7', 'start', 30, 14), M('pr_lha7', '')]},
    {'text': 'Unclicked 2× Tet2-Et flips it to <b>right-handed α</b>: −53, −46', 'hl': 0, 'marks': [M('pr_tet7', 'Tet', 30, 14)]},
    {'text': 'Clicked (+35, +53) and the ssDNA loop (+49, +33) put it <b>back in left-handed α</b>', 'hl': 0, 'marks': [M('pr_clk7', 'clicked', -30, -60), M('pr_ss7', 'ssDNA', 30, 14)]},
    {'text': 'With the dsDNA loop it <b>spreads widely</b>: −60, +38, sd 78° / 73°', 'hl': 0, 'marks': [M('pr_ds7', 'dsDNA', -30, 14)]},
    {'text': '<b>Asn553</b> stays helical (−52 to −60, −29 to −46) except with dsDNA (−69, +2, sd 40° / 79°)', 'hl': 0, 'marks': [M('pr_553', 'Asn553', 30, -50)]},
    {'text': '<b>Answer</b>: once clicked the tet groups barely matter; only the dsDNA loop disturbs the sites'}]},
'31d': {'title': 'PP2A: backbone rigidity (RMSF)', 'q': 'Does the spring deform the scaffold? Does the fold stay intact?',
  'items': [fig(P, 'pp2a_fig2_rmsf.png', 'Cα RMSF aligned on the core, HEAT repeats shaded'), fig(P, 'pp2a_fig5_heat_bend.png', 'Bend and twist per HEAT repeat; horseshoe aperture')],
  'explain': [
    {'item': 0, 'text': RMSF + 'x-axis: residue (UniProt numbering), HEAT repeats 1–15 shaded and numbered along the top; y-axis: Cα RMSF (Å), aligned on the core (residues 14–589). Dashed lines: the attachment sites Asp7 and Asn553. Apo, unclicked, clicked and 50 nt ssDNA are drawn; the 50 bp dsDNA run is off this scale.',
     'marks': [M('prm_tail', 'N-terminal tail', 30, -50), M('prm_core', 'rigid core', 30, -60), M('prm_asp7', 'Asp7', 30, 14), M('prm_asn553', 'Asn553', -30, -60), M('prm_h12', 'HEAT 12', 30, 10)]},
    {'item': 1, 'text': 'Left: <dfn>bend</dfn> at each HEAT repeat, the angle between the solenoid axis entering and leaving it (°). Middle: <dfn>superhelical twist</dfn> per step (°). Right: the horseshoe <dfn>aperture</dfn>, the HEAT1–HEAT15 centroid distance over time (Å); dotted line: 62.7 Å in the B56γ holoenzyme (Day et al.). Same colours; amber is the 50 bp dsDNA loop.',
     'marks': [M('phb_h12', 'HEAT 12', 30, 14), M('phb_ap_ds', 'dsDNA run opens', -30, 14), M('phb_holo', 'holoenzyme 62.7 Å', -30, 14)]}],
  'points': [
    {'text': 'Core RMSF <b>2.1–2.6 Å</b> for apo, tet, clicked and ssDNA; only the N-tail (residues 1–13) moves 10+ Å', 'hl': 0, 'marks': [M('prm_core', 'core', 30, -60), M('prm_tail', 'tail', 30, -50)]},
    {'text': '<b>HEAT 12</b> is the most variable repeat, the hinge region the literature reports (HEAT 11–13)', 'hl': 1, 'marks': [M('phb_h12', 'HEAT 12', 30, 14)]},
    {'text': 'With ssDNA the aperture changes by only <b>+0.8 ± 2.8 Å</b> (p = 0.77)', 'hl': 1},
    {'text': 'The <b>50 bp dsDNA</b> run comes apart: core RMSD 48 Å, helix fraction 0.26 (0.77 elsewhere), aperture 127 ± 43 Å', 'hl': 1, 'marks': [M('phb_ap_ds', 'dsDNA', -30, 14)]},
    {'text': 'Along the spring coordinate κ ≈ <b>1–2 kBT nm⁻²</b>, about 100× softer than sfGFP’s barrel'},
    {'text': '<b>Answer</b>: a single strand leaves the horseshoe intact; in this one 500 ns run the dsDNA loop unfolds it'}]},
'31e': {'title': 'PP2A: forces applied, force field vs model', 'q': 'How much force reaches the protein, and do the models agree with the force field?',
  'items': [fig(P, 'pp2a_fig7_headline.png', 'Cβ7–Cβ553 against the clicked zero'), fig(P, 'pp2a_fig4_force_extension.png', 'Bare-strand restraint windows against the models')],
  'explain': [
    {'item': 0, 'text': '<dfn>Cβ7–Cβ553</dfn>: the distance between the β-carbons of the two attachment residues, the coordinate the spring acts along (Å). Symbol: each system’s mean; thin line: its full range; bars: block SEM. <dfn>Δ</dfn> is the shift from the clicked protein without DNA (the right zero), with Welch’s p on ten block means.',
     'marks': [M('ph_clk', 'the zero: clicked, no DNA', 30, -60), M('ph_apo', 'apo', 30, 14), M('ph_ss', 'ssDNA', 30, 14), M('ph_ds', 'dsDNA', -30, -60)]},
    {'item': 1, 'text': 'x-axis: the spring’s end-to-end distance (Å); y-axis: force on the anchors (pN), + pushes them apart, − pulls them together. Lines are models: <dfn>Marko–Siggia</dfn> (worm-like chain: a single strand is an entropic spring that can only pull) and <dfn>Zocchi f(x)</dfn> (a bent duplex rod that pushes, and softens once it kinks). Points: forces measured on the bare strands held at 60 / 80 / 100 Å. Shading: where each chimera actually sits.',
     'marks': [M('pfe_ss50', 'measured ss50', 30, 14), M('pfe_ds50', 'measured ds50', 30, -60), M('pfe_ms50', 'Marko–Siggia (pull)', 30, 14), M('pfe_zocchi50', 'Zocchi duplex (push)', 30, -55)]}],
  'points': [
    {'text': 'Against the clicked zero, ssDNA: <b>−3.1 ± 3.5 Å</b> (p = 0.39), not resolved', 'hl': 0, 'marks': [M('ph_ss', 'ssDNA', 30, 14), M('ph_clk', '')]},
    {'text': 'dsDNA: <b>+34.4 ± 8.8 Å</b> (p = 0.003): the sites are pushed apart, as the horseshoe unfolds', 'hl': 0, 'marks': [M('ph_ds', 'dsDNA', -30, -60)]},
    {'text': 'The strand sits at 95 Å, its own rms size: Marko–Siggia gives a <b>1.74 pN pull</b>', 'hl': 1, 'marks': [M('pfe_ms50', 'pull', 30, 14)]},
    {'text': 'Bare ss50 agrees with Marko–Siggia within error (−1.23 ± 1.70 pN at 100 Å)', 'hl': 1, 'marks': [M('pfe_ss50', 'ss50', 30, 14)]},
    {'text': 'Bare ds50 pushes <b>+10 to +13 pN</b> at 60–100 Å, about 13× Zocchi’s 0.9 pN', 'hl': 1, 'marks': [M('pfe_ds50', 'ds50', 30, -60), M('pfe_zocchi50', 'model', 30, -55)]},
    {'text': '<b>Answer</b>: the duplex pushes far harder than modelled; a single strand pulls weakly and moves nothing measurable'}]},
'31f': {'title': 'PP2A: DNA dynamics', 'q': 'Is the DNA kinked or coiled? Does it interact with the protein beyond pulling?',
  'items': [fig(P, 'pp2a_fig6_dna_kink.png', 'Bend, sharp turns, torsions and stacking along each strand'), fig(P, 'pp2a_fig8_strand.png', 'Strand size and modelled force over time')],
  'explain': [
    {'item': 0, 'text': 'Top row: 50 nt ssDNA; bottom row: the modified strand of the 50 bp dsDNA. Left: <dfn>bend</dfn> of the C4′ backbone trace at each nucleotide (°); pink bars (right axis): fraction of frames with a sharp > 90° turn, i.e. a <dfn>kink</dfn>. Middle: <dfn>BII fraction</dfn> (backbone torsion ε − ζ > 0, a sign of local strain) and how often a base is <dfn>unstacked</dfn> from its neighbour. Right: mean curvature against end-to-end distance; dotted lines: kink thresholds 2.6 (Zocchi), 5.5 (onset), 8.7 °/bp (certain).',
     'marks': [M('pk_ss_hot', 'turns at nt 3', 30, 14), M('pk_ss_stack', 'BII / unstacking', 30, -55), M('pk_ds_open', 'sharpest dsDNA turn', 30, -55), M('pk_ss_curv', 'curvature vs span', -30, 14)]},
    {'item': 1, 'text': 'Left: each strand’s end-to-end distance over time (Å, solid) and radius of gyration (dotted). Grey line: the rms size of a free 50-mer (94 Å); a single strand at that size is neither slack nor taut. Right: the force the models assign frame by frame (pN): Marko–Siggia pull for ssDNA, Zocchi push for dsDNA (drawn negative).',
     'marks': [M('ps_rms', 'free 50-mer rms 94 Å', 30, 14), M('ps_ss', 'ssDNA', 30, -55), M('ps_ds', 'dsDNA', 30, 14), M('ps_ssforce', 'ssDNA pull ≈ 1.7 pN', -30, -60)]}],
  'points': [
    {'text': 'ssDNA end-to-end <b>95 ± 9 Å</b>, right at its free rms size: neither slack nor taut', 'hl': 1, 'marks': [M('ps_ss', 'ssDNA', 30, -55), M('ps_rms', '94 Å', 30, 14)]},
    {'text': 'Sharp turns in 18 % of frames, hottest at nt 3 by the tether: <b>no interior kink</b>', 'hl': 0, 'marks': [M('pk_ss_hot', 'nt 3', 30, 14)]},
    {'text': '78 % of adjacent bases stay stacked; 15 % of the backbone in BII', 'hl': 0, 'marks': [M('pk_ss_stack', '', 30, -55)]},
    {'text': '<b>3 of 50 nt</b> touch the protein at any time, mostly the N-tail (Met1–Asp5, Tet sites, Gln557)'},
    {'text': 'dsDNA: 8.4°/bp with ≈ 36 base pairs open on average: the duplex <b>melts</b> rather than kinking cleanly', 'hl': 0, 'marks': [M('pk_ds_open', 'dsDNA', 30, -55)]},
    {'text': '<b>Answer</b>: the strand is a relaxed coil that touches the protein; the dsDNA loop comes apart'}]},
}
deck = json.loads((ROOT / 'deck.json').read_text())
for sl in deck['slides']:
    pg = str(sl.get('page'))
    if pg in slides:
        sec = sl.get('section', 'md')
        sl.clear(); sl.update(kind='media', page=pg, section=sec, figW=1180, exH=590, **slides[pg])
    if pg == '31a':
        for p in sl['points']:
            p['text'] = p['text'].replace('425 ns over 5 systems', '500 ns over 5 systems')
(ROOT / 'deck.json').write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print('updated', len(slides), 'slides')
