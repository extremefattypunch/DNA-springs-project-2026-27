#!/usr/bin/env python3
"""Write the MD analysis slides (explain phase + answers) into deck.json.

Each figure is shown full-slide first with its axes and terms explained and its
points of interest marked, then the slide shrinks to the figure | answers view.
Mark positions are fractions of the image, in md/anns_*.json.
"""
import html, json, pathlib
from PIL import Image
ROOT = pathlib.Path(__file__).resolve().parent.parent
ANN = {**json.loads((ROOT / 'md/anns_sfgfp.json').read_text()), **json.loads((ROOT / 'md/anns_pp2a.json').read_text())}

def fig(d, name, cap):
    w, h = Image.open(ROOT / 'md' / d / name).size
    return {'img': f'md/{d}/{name}', 'w': w, 'h': h, 'cap': cap, 'anns': ANN.get(name, {})}

def tex(src, display=False):
    """LaTeX, typeset in the deck by KaTeX (MathML output)."""
    return f'<span class="tex{" d" if display else ""}">{html.escape(src, quote=False)}</span>'

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
'12a': {'title': 'sfGFP MD: the protein and the simulated systems', 'layout': 'row', 'rowH': 430, 'exH': 560,
  'items': [fig(S, 'structure_27bp.png', '<b>Structure</b>: what each part is'),
    {'video': 'md/sfgfp/S3_spring27_overview.mp4', 'poster': 'md/sfgfp/S3_spring27_overview_still.png', 'cap': '<b>Overview</b>: barrel, two Tet2-Et/sTCO tethers, 27 bp duplex'},
    {'video': 'md/sfgfp/S3_spring27_strain.mp4', 'poster': 'md/sfgfp/S3_spring27_strain_still.png', 'cap': '<b>Strain</b>: putty tube scaled by per-residue RMSF'},
    {'video': 'md/sfgfp/S3_spring27_chromophore.mp4', 'poster': 'md/sfgfp/S3_spring27_chromophore_still.png', 'cap': '<b>Chromophore</b>: the three H-bonds that stop it twisting'}],
  'explain': [
    {'item': 0, 'text': '<b>sfGFP architecture</b> (238 residues). <dfn>β-strands</dfn>: extended segments hydrogen-bonded side by side; eleven of them form the <dfn>β-barrel</dfn> wall (≈ 42 Å tall, 24 Å across), the rigid part. '
      'A central <dfn>α-helix</dfn> threads the barrel axis and carries the <dfn>chromophore</dfn> (made from Thr65–Tyr66–Gly67), sealed off from solvent; short helices also sit on the caps. '
      '<dfn>Loops</dfn>: the short, irregular segments that join one strand to the next at the two ends (caps) of the barrel. Nothing holds them in the sheet, so they are the most mobile parts (e.g. residues 189–194 between strands 9 and 10).',
     'marks': [M('st_strand', 'β-strand', -40, -90), M('st_loop', 'loop', 40, 0), M('st_helix', 'short α-helix (cap)', -30, 34), M('st_cro', 'chromophore (inside)', 40, 36)]},
    {'item': 0, 'text': '<b>The spring.</b> Our tethers (orange, Tet2-Et clicked to sTCO) attach at <b>Asp133</b> and <b>Asn149</b> (red), 31.3 Å apart across the barrel, with His148, which bonds the chromophore, between them. The DNA spring (blue) spans the two tethers and pushes or pulls on the sites.',
     'marks': [M('st_site1', 'Asp133', -40, -60), M('st_site2', 'Asn149', 40, 30), M('st_tether', 'tether', -40, -30), M('st_dna', 'DNA spring', 40, -70)]}],
  'points': [
    {'text': '<b>Sites</b>: Asp133 and Asn149, 31.3 Å apart, with His148 between them bonded to the chromophore', 'hl': 0},
    {'text': '<b>Systems compared</b>: WT sfGFP · tethers unclicked vs clicked (no DNA, the zero) · DNA springs of <b>various lengths</b> · dsDNA vs ssDNA · intact vs nicked duplex · a direct force clamp'},
    {'text': '<b>Animations</b> show the 27 bp spring chimera (50 ns, water stripped, superposed on the barrel core)', 'hl': [1, 2, 3]},
    {'text': '<b>Scale</b>: 4,050 ns of production, 81 replicates across 17 systems'}]},
'12c': {'title': 'sfGFP: backbone rigidity (RMSF)', 'q': 'Is the backbone rigid, with and without force?',
  'items': [fig(S, 'fig2_rmsf.png', 'Cα RMSF, aligned on the barrel core'), fig(S, 'fig1_clamp_only.png', 'Axial force clamp on the unclicked protein')],
  'explain': [
    {'item': 0, 'text': RMSF + 'x-axis: residue number (2B3P numbering); y-axis: Cα RMSF (Å). Dashed lines mark the chromophore (CRO 66) and the attachment sites Asp133, His148, Asn149. One colour per system. '
      'The low, flat stretches (≈ 0.4 Å) are the <dfn>β-strands</dfn> of the barrel wall. The peaks are <dfn>loops</dfn>, the short connectors between strands at the barrel caps: the tallest is residues <b>189–194</b> (between strands 9 and 10, ≈ 1.5 Å), then 155–158 and 211–212.',
     'marks': [M('rmsf_barrel', 'β-strand: rigid', 30, -60), M('rmsf_loop', 'loop 189–194: mobile', 30, 10), M('rmsf_cro', 'chromophore', 30, -50), M('rmsf_sites', 'attachment sites', 30, -50)]},
    {'item': 1, 'text': '<dfn>Axial force clamp</dfn>: a constant force F applied directly between the Cβ atoms of Asp133 and Asn149 of the unclicked protein, pulling them apart, with no DNA: a controlled stand-in for a spring. '
      'x-axis: F (pN); y-axis: Cβ–Cβ distance (Å), mean ± sem of 3 replicates. Shaded: ±0.5 Å around the unloaded (0 pN) distance.',
     'marks': [M('cl_0', '0 pN: 31.61 Å', 30, -60), M('cl_20', '20 pN: 31.75 Å', -30, -60), M('cl_band', '±0.5 Å', 30, -50)]}],
  'points': [
    {'text': 'The β-barrel stays <b>uniformly rigid</b> (≈ 0.4 Å on the strands) in every system, including WT with no force at all', 'hl': 0, 'marks': [M('rmsf_barrel', 'rigid strand', 30, -60)]},
    {'text': 'RMSF traces <b>overlie each other</b>; only the loops and the tethers move', 'hl': 0, 'marks': [M('rmsf_loop', 'loop 189–194', 30, 10)]},
    {'text': 'With <b>no force</b>, the site distance is 31.61 ± 0.07 Å: the unloaded barrel barely fluctuates', 'hl': 1, 'marks': [M('cl_0', '0 pN', 30, -60)]},
    {'text': '0 → 20 pN shifts it by only <b>0.14 Å</b>; the 4 pN dip is replicate noise. Compliance 0.0105 ± 0.0057 Å/pN (includes zero), so <b>κ > 107 k<sub>B</sub>T nm⁻²</b>', 'hl': 1, 'marks': [M('cl_20', '20 pN', -30, -60), M('cl_4', 'noise', 30, 14)]},
    {'text': '<b>Answer</b>: the backbone is already super rigid without any force, and force barely changes it'}]},
'12d': {'title': 'sfGFP: forces applied, force field vs model', 'q': 'How hard does a bent duplex push in the force field?',
  'items': [fig(S, 'fig6_force_extension.png', 'Force measured on the bare duplex, one point per run')],
  'explain': [
    {'item': 0, 'text': '<b>The experiment.</b> A bare <dfn>duplex</dfn> (the double-stranded DNA alone: 27 bp, 40 bp or 40 bp nicked; no protein) is held by a <dfn>harmonic restraint</dfn> on its <dfn>end-to-end separation</dfn> r, the distance between the two linker nitrogens that would attach to the protein:'
      + tex(r'U(r) = \tfrac{1}{2}\,k\,(r - x_0)^2, \qquad F = k\,(\langle r\rangle - x_0)', True)
      + tex('r') + ': the end-to-end separation in each frame; ' + tex(r'x_0') + ': the separation the restraint is centred on (one per run); ' + tex('k') + ': restraint stiffness, 200 (stiff) or 33 (soft) kJ mol⁻¹ nm⁻². '
      'The duplex wants to straighten, so it sits at ' + tex(r'\langle r\rangle > x_0') + '; the restraint then holds it with exactly the force ' + tex('F') + ' that the duplex pushes with.'},
    {'item': 0, 'text': '<b>Reading the plot.</b> x-axis: ' + tex(r'\langle r\rangle') + ', where the duplex actually sat (Å); y-axis: ' + tex('F') + ' (pN). Every point is one run. <b>Colour</b> = duplex (blue 27 bp, red 40 bp, green 40 bp nicked); <b>shape</b> = restraint (circle stiff, square soft). '
      'Solid lines: mean force per duplex. Dashed: the <dfn>Euler buckling load</dfn> ' + tex(r'F_E = \pi^2 B / L^2') + ' (B = 200 pN nm², L = contour length), the push of a smoothly bent elastic rod. Grey band: what <dfn>Zocchi’s kinked branch</dfn> predicts for the same duplexes.',
     'marks': [M('fe_27', '27 bp: 28 pN', 30, -60), M('fe_euler27', 'Euler load, 27 bp', 30, 14), M('fe_40', '40 bp: 15.8 pN', 30, 10), M('fe_zocchi', 'Zocchi kinked branch', 30, -55)]}],
  'points': [
    {'text': 'Measured: <b>27 bp 28.0 ± 1.2 pN</b> (25 runs), 40 bp 15.8 ± 1.6 pN (10 runs), 40 bp nicked 13.2 pN', 'hl': 0, 'marks': [M('fe_27', '28 pN', 30, -60), M('fe_40', '15.8 pN', 30, 10)]},
    {'text': 'Roughly <b>flat</b> across ⟨r⟩: the duplex pushes with about the same force however far its ends are pushed together, as a buckled rod does', 'hl': 0},
    {'text': 'Close to the <b>Euler load</b> (dashed: 26.8 pN for 27 bp, 11.9 pN for 40 bp)', 'hl': 0, 'marks': [M('fe_euler27', 'Euler', 30, 14)]},
    {'text': 'Zocchi’s <b>kinked</b> branch (grey) predicts only 2.4–6.6 pN: 5–6× less', 'hl': 0, 'marks': [M('fe_zocchi', 'kinked', 30, -55)]},
    {'text': '<b>What we learn</b>: the force-field duplex behaves as a smoothly bent, unkinked rod pushing tens of pN. The next two slides show why Zocchi’s model expects less'}]},
'12e': {'title': 'sfGFP: DNA dynamics', 'q': 'Is the DNA kinked or coiled? Does it interact with the protein beyond pushing?',
  'items': [fig(S, 'fig4_dna_bend.png', 'Axis bend per base-pair step along each duplex')],
  'explain': [
    {'item': 0, 'text': 'x-axis: base pair along the duplex in the chimera; y-axis: <dfn>axis bend</dfn> at that base-pair step (°), the angle between the local helix axes either side of it. Line: mean across replicates; band: ± sd. A <dfn>kink</dfn> is a sharp local break of the duplex: it would show as one narrow spike. Smooth bending spreads the angle evenly.',
     'marks': [M('db_27peak', 'highest step, 27 bp', 30, -60), M('db_40', '40 bp: even bending', 30, 14)]}],
  'points': [
    {'text': 'The duplex bends <b>smoothly</b> at 4.8–5.2°/bp, nucleosome-like, with no local spike', 'hl': 0, 'marks': [M('db_40', 'even', 30, 14)]},
    {'text': 'The highest single step stays a broad bump, not a kink', 'hl': 0, 'marks': [M('db_27peak', 'broad', 30, -60)]},
    {'text': 'The DNA touches the protein in 72–98 % of frames (ssDNA loop 93 %): <b>local contact</b>, not action at a distance'},
    {'text': '<b>Answer</b>: no kink, no coiling of the duplex; but the DNA lies against the protein'}]},
'12f': {'title': 'sfGFP: why Zocchi’s model puts our springs in the kinked regime', 'items': [], 'exH': 590,
  'q': 'Zocchi’s equations (Molecular Machines §2.15) evaluated with our MD numbers',
  'points': [
    {'text': '<b>Inputs.</b> Zocchi’s rod parameters, and what our MD measured: in the sfGFP chimeras the duplex ends sit <b>x = 63–73 Å</b> apart, bent <b>4.8–5.2° per bp</b>'
      + tex(r'L = 0.33\,\mathrm{nm}\times N_{\mathrm{bp}}, \qquad B = k_BT\,\ell_p \approx 200\ \mathrm{pN\,nm^2}, \qquad \tau_c = 34\ (\text{intact}),\ 27\ (\text{nicked})\ \mathrm{pN\,nm}', True)},
    {'text': '<b>Step 1: the kink parameter</b> (eq. 2.194)'
      + tex(r'\gamma = \frac{L\,\tau_c}{2B} = \frac{8.91 \times 34}{2 \times 200} = 0.76\ \ (27\ \text{bp}), \qquad \frac{13.2 \times 34}{2 \times 200} = 1.12\ \ (40\ \text{bp})', True)},
    {'text': '<b>Step 2: the critical distance</b> (eq. 2.188): below x<sub>c</sub> the bending torque would exceed τ<sub>c</sub>, so the rod kinks instead'
      + tex(r'x_c = L\Big(1-\frac{\gamma^2}{15}\Big) = 8.91\,(1 - 0.038) = 8.57\ \text{nm} = 86\ \text{Å}\ \ (27\ \text{bp}), \qquad 121\ \text{Å}\ \ (40\ \text{bp})', True)},
    {'text': '<b>Step 3: compare with the MD.</b>'
      + tex(r'x_{\mathrm{MD}} = 63\text{–}73\ \text{Å} \;<\; x_c = 86\ \text{Å}\ (27\ \text{bp}),\ 121\ \text{Å}\ (40\ \text{bp}) \quad\Longrightarrow\quad \text{kinked regime}', True)},
    {'text': '<b>Same answer from curvature</b>: a rod kinks once its bending torque Bκ reaches τ<sub>c</sub>'
      + tex(r'\kappa_c = \frac{\tau_c}{B} = \frac{34}{200}\ \text{nm}^{-1} \;\Rightarrow\; \kappa_c \times 0.33\ \tfrac{\text{nm}}{\text{bp}} = 3.2^\circ/\text{bp}\ \ (2.6^\circ/\text{bp nicked}) \;<\; 4.8\text{–}5.2^\circ/\text{bp in the MD}', True)},
    {'text': '<b>So the model predicts the kinked-branch force</b> (eq. 2.196), only a few pN:'
      + tex(r'f_{\mathrm{kinked}}(x) = \frac{\tau_c}{2R\sqrt{1-(x/2R)^2}},\quad R = \frac{L}{2}\Big(1-\frac{\gamma^2}{90}\Big) \;\Rightarrow\; 5.5\text{–}6.8\ \text{pN}\ (27\ \text{bp}),\ \ 3.0\text{–}3.2\ \text{pN}\ (40\ \text{bp})', True)
      + 'Yet the simulated duplexes <b>never kinked</b> (smooth bending, DNA dynamics slide): next slide'}]},
'12f2': {'title': 'sfGFP: kinked vs unkinked regime, and the force field', 'q': 'What changes between the regimes, and which one does the force field follow?',
  'items': [fig(S, 'zocchi_regimes.png', 'Zocchi’s force against end-to-end distance, both branches'),
            fig(S, 'zocchi_unkinked_vs_ff.png', 'Zocchi recomputed at the force field’s own distances')],
  'explain': [
    {'item': 0, 'text': 'x-axis: end-to-end distance x of the duplex; y-axis: the force f it pushes its ends apart with. <b>Solid</b>: the unkinked branch where Zocchi applies it (x &gt; x<sub>c</sub>), '
      + tex(r'f = 10B/L^2 - k_BT/(L-x)') + '. <b>Dotted</b>: the kinked branch (x &lt; x<sub>c</sub>), ' + tex(r'f = \tau_c / \big(2R\sqrt{1-(x/2R)^2}\big)')
      + '. <b>Dashed</b>: the unkinked branch continued below x<sub>c</sub>, i.e. the push of a rod that never kinks. Grey band: where our chimeras hold the duplex.',
     'marks': [M('zr_xc27', 'x_c, 27 bp', 30, -60), M('zr_kink27', 'kinked: 5–7 pN', -30, 20), M('zr_unk27', 'unkinked: ≈ 24 pN', -30, -60), M('zr_chim', 'our chimeras', 30, 14)]},
    {'item': 1, 'text': 'Bare duplexes in the force field (no protein), held at ⟨x⟩ = 49 Å (27 bp) and 62 Å (40 bp), both below x<sub>c</sub>. For each duplex: <b>grey</b> = Zocchi’s kinked prediction at that ⟨x⟩; '
      '<b>light blue</b> = Zocchi recomputed on the unkinked branch at the same ⟨x⟩; <b>dark blue</b> = what the force field measured, ' + tex(r'\langle k(r - x_0)\rangle') + ' (mean ± sem; dots are single runs).',
     'marks': [M('zu_k0', 'kinked 4.6', -30, -60), M('zu_u0', 'unkinked 24.1', -30, -70), M('zu_f0', 'measured 28.0', 30, -60)]}],
  'points': [
    {'text': '<b>Kinked regime</b> (x &lt; x<sub>c</sub>): the kink takes up the bending, so the push collapses to a few pN', 'hl': 0, 'marks': [M('zr_kink27', '5–7 pN', 30, 14), M('zr_kink40', '≈ 3 pN', 30, 14)]},
    {'text': '<b>Unkinked regime</b>: the push stays near the Euler load π²B/L² ≈ 10B/L² (24.9 pN for 27 bp, 11.3 pN for 40 bp), almost independent of x', 'hl': 0, 'marks': [M('zr_unk27', '≈ 24 pN', -30, -60), M('zr_unk40', '≈ 11 pN', 30, -50)]},
    {'text': 'Our chimeras (grey band) sit well inside the kinked regime of the model: <b>5×</b> less force than if the rod stays unkinked', 'hl': 0, 'marks': [M('zr_chim', 'chimeras', 30, 14), M('zr_xc27', 'x_c', 30, -60)]},
    {'text': 'Recomputed on the <b>unkinked</b> branch at the force field’s own ⟨x⟩: <b>24.1 vs 28.0 pN</b> (27 bp), 10.9 vs 15.8 (40 bp), 10.9 vs 13.2 (nicked): within 1.2–1.5×', 'hl': 1,
     'marks': [M('zu_u0', '24.1', -30, -70), M('zu_f0', '28.0', 30, -60), M('zu_u1', '10.9', -30, -60), M('zu_f1', '15.8', 30, -60)]},
    {'text': 'The kinked branch is <b>5–6× too low</b> (4.6, 3.0, 2.3 pN)', 'hl': 1, 'marks': [M('zu_k0', '4.6', -30, -60), M('zu_k1', '3.0', -30, -60), M('zu_k2', '2.3', -30, -60)]},
    {'text': '<b>Answer</b>: the force field does not contradict Zocchi’s spring. Its duplex never kinks, so it follows his <b>unkinked</b> branch; the only disagreement is whether DNA kinks at 4.8–5.2°/bp'}]},
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
        sl.clear(); sl.update({'kind': 'media', 'page': pg, 'section': sec, 'figW': 1180, 'exH': 590, **slides[pg]})
    if pg == '31a':
        for p in sl['points']:
            p['text'] = p['text'].replace('425 ns over 5 systems', '500 ns over 5 systems')
# slides added after the original deck: insert after their predecessor
for pg, after in (('12f2', '12f'),):
    if not any(str(x.get('page')) == pg for x in deck['slides']):
        i = next(i for i, x in enumerate(deck['slides']) if str(x.get('page')) == after)
        deck['slides'].insert(i + 1, {'kind': 'media', 'page': pg, 'section': 'md', 'figW': 1180, 'exH': 590, **slides[pg]})
(ROOT / 'deck.json').write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print('updated', len(slides), 'slides')
