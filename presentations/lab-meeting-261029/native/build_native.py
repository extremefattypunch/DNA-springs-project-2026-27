#!/usr/bin/env python3
"""Rebuild the original Google-Slides pages as native deck slides.

Replaces every `kind: "image"` slide (a render of the original PDF page) with a
deck-styled slide carrying the same content, so the whole deck shares one look,
and inserts the "Why DNA springs?" intro before the SAM slide. Each top-level
bullet group appears on its own click. Re-running is idempotent.
"""
import html, json, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent


def tex(src, display=False):
    return f'<span class="tex{" d" if display else ""}">{html.escape(src, quote=False)}</span>'


def done(t): return f'<span class="t-done">[DONE]: {t}</span>'
def pend(t): return f'<span class="t-pend">[PENDING]: {t}</span>'
def todo(t): return f'<span class="t-todo">[To-Do]: {t}</span>'


def grp(head, *items):
    sub = ''.join(f'<li>{i}</li>' for i in items)
    return f'<span class="lead">{head}</span>' + (f'<ul class="sub">{sub}</ul>' if items else '')


AGO = {'text': 'Project status: 4 months ago', 'color': 'var(--sec-status)'}
NOW = {'text': 'Project status: current', 'color': '#00838f'}


def text(page, title, section, groups, q=None, badge=None, aside=None, items=None, figW=None):
    s = {'kind': 'media', 'page': page, 'section': section, 'title': title, 'items': items or [], 'big': True,
         'points': [{'text': g} for g in groups]}
    if q: s['q'] = q
    if badge: s['badge'] = badge
    if aside: s['aside'] = aside
    if figW: s['figW'] = figW
    return s


DAY = 'Day et al., 2025, <i>Structure</i> 33:1688'
DAY_TXT = 'Day et al., 2025, Structure 33:1688'
DAYQ = 'Method and results from ' + DAY_TXT + ' (Aα vs Aβ); our version compares Aα vs Aα + DNA spring'
def anim(page, title, src, q, cap, pts):
    return {'kind': 'media', 'page': page, 'section': 'assay', 'title': title, 'figW': 1180, 'msgChannel': 'cryo', 'q': q,
            'items': [{'iframe': f'interactive/{src}', 'title': title, 'cap': cap}],
            'points': [{'text': t, 'msg': {'step': i + 1}} for i, t in enumerate(pts)]}


def pic(name, cap):
    from PIL import Image
    w, h = Image.open(ROOT / 'native' / name).size
    return {'img': f'native/{name}', 'w': w, 'h': h, 'cap': cap}


A23 = {'kind': 'media', 'page': 23, 'section': 'assay', 'title': '[PENDING] PP2A assays: following Day et al. 2025', 'figW': 1180, 'msgChannel': 'cryo',
  'q': 'Template: Day et al., 2025, Structure 33:1688, which compared the two A-subunit isoforms',
  'items': [{'iframe': 'interactive/pp2a.html', 'title': 'Aα vs Aβ vs Aα + spring', 'cap': 'Schematic after Day et al. Fig. 1D and 2 (not to scale)'}],
  'points': [{'text': t, 'msg': {'step': st}} for t, st in [
    ('<b>Day et al.</b>: Aα vs Aβ holoenzymes. Aβ’s extra 12 N-terminal residues push the B56γ helix away; its scaffold is more concave', 4),
    ('<b>Ours</b>: Aα vs Aα <b>with a DNA spring attached</b>: a perturbed A subunit, the same question', 5),
    ('<b>Three readouts</b>, as in Day et al.: activity (DiFMUP), B56 binding (BLI), structure (cryoEM). Activity alone is not enough: Day et al. found the same activity but up to 4.7× weaker binding', 6)]]}

PP_FN = anim('22a', 'PP2A: what it does in the cell', 'pp2afn.html', 'Protein phosphatase 2A: a serine/threonine phosphatase, the cell’s main “off switch”',
    'Schematic of kinase/phosphatase signalling', [
    'Signalling proteins are switched by <b>phosphate groups</b> on serine/threonine residues',
    'A <b>kinase</b> adds the phosphate (from ATP): <b>ON</b>',
    '<b>PP2A</b> removes it (dephosphorylation): <b>OFF</b>',
    'PP2A and PP1 together carry out ≈ <b>90 %</b> of the cell’s Ser/Thr dephosphorylation',
    'Targets include AKT, ERK, c-Myc, p53, tau and mitotic CDK substrates: growth, proliferation, cell cycle, neurons',
    'A <b>tumour suppressor</b>: inactivated in many cancers (A-subunit hot spots P179, R183) and by viral oncoproteins (SV40 small t); hyperphosphorylated tau in Alzheimer’s'])

PP_ABC = anim('22b', 'PP2A architecture: A, C and B subunits', 'pp2a.html', 'PP2A works as a three-subunit holoenzyme: scaffold A, catalytic C, regulatory B',
    'Schematic after Day et al. Fig. 1D and 2 (not to scale)', [
    '<b>A, scaffold</b> (PPP2R1A/B, ≈ 65 kDa): a horseshoe of 15 HEAT repeats with no activity of its own; it holds the other two in place',
    '<b>C, catalytic</b> (PPP2CA/B): the phosphatase itself; binds HEAT 11–15. A + C = the <b>AC core dimer</b>, which forms first. Its C-terminal tail (methylated at Leu309) helps B bind',
    '<b>B, regulatory</b>: binds the AC core on HEAT 1–10 to complete the <b>ABC holoenzyme</b>. 4 families (B55, B56, PR72, striatin), > 40 isoforms; B decides <b>which substrate</b> and <b>where</b> in the cell',
    '<b>Is assembly needed?</b> Free C can dephosphorylate small artificial substrates (e.g. DiFMUP), but in cells it is kept scarce; physiological targeting needs the holoenzyme, with B docking the substrate (e.g. LxxIxE motifs on B56) next to C'])

PP_HEAT = anim('22c', 'HEAT repeats, and our Tet sites on Aα', 'heat.html', 'Why the A subunit is a flexible scaffold, and where the DNA spring attaches',
    'Schematic of the HEAT-repeat fold (not to scale)', [
    'A <b>HEAT repeat</b>: ≈ 39 residues forming two antiparallel α-helices joined by a turn (named after Huntingtin, EF3, PP2A A, TOR)',
    'Repeats stack side by side through hydrophobic contacts into one continuous, spring-like <b>solenoid</b>',
    'Aα: 15 repeats curve into a <b>horseshoe</b>. B binds HEAT 1–10, C binds HEAT 11–15; a hinge near HEAT 11–12 lets the arc flex, so it can fit > 40 different B subunits',
    '<b>Our construct</b>: Aα with Tet groups at <b>D6X</b> (N-terminus) and <b>N552X</b> (HEAT 14), so the DNA spring spans the horseshoe'])

PP_ISO = {'kind': 'media', 'page': '22d', 'section': 'assay', 'title': 'Isoforms: Aα vs Aβ, and B56α vs B56γ', 'figW': 760,
  'q': 'Which subunits we use, and what Day et al. compared',
  'items': [pic('pp2a_tet_sites.png', '<b>Our system</b>: Aα holoenzyme model (white A, tan B, lilac C) with the Tet sites D6X and N552X (green)'),
            pic('day_fig1.png', 'Day et al. Fig. 1C–D: Aβ’s N-terminus is mobile in MD and blocks the B56γ C-terminal helix')],
  'points': [{'text': t, 'hl': h} for t, h in [
    ('<b>Aα</b> (PPP2R1A): the dominant A isoform, high in all tissues; <b>the one we use</b>. <b>Aβ</b> (PPP2R1B): 87 % identical, lower and more variable expression, with a 12-residue N-terminal extension', 0),
    ('Functional differences: SV40 small t binds only Aα; Aβ-specific dephosphorylation of RalA; Aβ monomer is less stable and aggregates (Day et al.)', 1),
    ('<b>B56</b> family (PPP2R5): five isoforms α, β, γ, δ, ε, sharing a central core but differing at their N- and C-termini, which set localisation and partners', 1),
    ('<b>B56α</b>: mostly cytoplasmic; targets e.g. c-Myc (Ser62). <b>B56γ</b>: mostly nuclear; targets e.g. p53 (Thr55), a tumour suppressor. Its C-terminal helix contacts the A N-terminus', 1),
    ('<b>Day et al.</b>: Aβ binds B56α 4.7× and B56γ 1.5× more weakly than Aα. We compare Aα ± spring with the same B56 subunits', 1)]]}

A24 = {'kind': 'media', 'page': 24, 'section': 'assay', 'title': '[Plans] Holoenzyme formation analysis (cryoEM)', 'figW': 1180, 'exH': 600,
  'q': 'Day et al. solved both holoenzymes by cryoEM; we would do the same for Aα ± spring',
  'items': [pic('day_fig2.png', 'Day et al. Fig. 2: cryoEM maps, Aα (red) and Aβ (grey) holoenzymes, B56γ blue, C orange'),
            pic('day_fig3.png', 'Day et al. Fig. 3A–B: A subunits superposed, and per-HEAT-repeat shift')],
  'explain': [
    {'item': 0, 'text': '<b>What cryoEM gives</b>: a 3D density map of the whole holoenzyme (A, B56γ, C), into which an atomic model is built. Day et al. reached <b>3.4 Å</b> (Aα) and <b>3.2 Å</b> (Aβ), maps sharpened with EMReady. '
       '<b>Right (6 Å low-pass)</b>: blurring the map reveals the thin B56γ C-terminal helix (blue, down to the A N-terminus) in Aα, but not in Aβ.'},
    {'item': 1, 'text': '<b>Measuring the deformation</b>: superpose the two A subunits on their C-binding end (HEAT 12–15), then compare the centre of mass of each HEAT repeat. '
       '<b>Left</b>: the horseshoes overlaid; Aβ spans 59.6 Å from HEAT 1 to 15 vs 62.7 Å for Aα (more concave). <b>Right</b>: the shift per repeat, largest at HEAT 1–8 (up to ≈ 3 Å), near zero at the aligned end.'}],
  'points': [
    {'text': '<b>Structural deformation</b>: align on HEAT 12–15 and measure HEAT-repeat centre-of-mass shifts, WT vs spring attached', 'hl': 1},
    {'text': '<b>Day et al.</b>: Aβ vs Aα shifts up to ≈ 3 Å; span 59.6 vs 62.7 Å', 'hl': 1},
    {'text': '<b>Ascertain spring attachment</b>: resolution level, neural-network particle picking, refinement to an atomic model', 'hl': 0},
    {'text': '<b>Day et al.</b>: 3D classification resolved the B56γ C-terminal helix (Aα only) and the C-tail', 'hl': 0}]}

A29 = anim(29, '[Plans] B56-binding assay: biolayer interferometry (BLI)', 'bli.html', 'Method of Day et al. 2025: does the A subunit change how tightly B56 binds the AC core?',
    'Schematic of BLI; curves computed from Day et al.’s fitted rate constants (Fig. 5B, B56α); timings illustrative', [
    '<b>Biosensor</b>: an optical-fibre tip coated with B56 (AR2G tip, amine coupling). Bound mass shifts the reflected light by Δλ',
    '<b>Baseline</b> in buffer: nothing binds',
    '<b>Association</b>: dip into the AC core dimer (31–500 nM); AC binds B56 and Δλ rises, at a rate set by k<sub>on</sub>',
    '<b>Dissociation</b>: back into buffer; AC falls off and Δλ decays, at a rate set by k<sub>off</sub>',
    '<b>Fit</b> (1:1 binding):' + tex(r'AC + B \;\xrightleftharpoons[k_{\mathrm{off}}]{\;k_{\mathrm{on}}\;}\; ABC, \qquad K_D = \frac{k_{\mathrm{off}}}{k_{\mathrm{on}}}', True) + 'smaller K<sub>D</sub> = tighter binding'])
A29b = {'kind': 'media', 'page': '29b', 'section': 'assay', 'title': 'BLI results (Day et al.): Aβ binds B56 more weakly', 'figW': 900,
  'q': 'Day et al. Fig. 5 (n = 2): AαChimera and AβNT are their swap constructs',
  'items': [pic('day_fig5.png', 'Day et al. Fig. 5A: K_D of AC core dimers for B56α and B56γ')],
  'points': [{'text': t} for t in [
    '<table class="tb"><tr><th></th><th>AαC K<sub>D</sub></th><th>AβC K<sub>D</sub></th><th>Aβ vs Aα</th></tr><tr><td>B56α</td><td>25.7 nM</td><td>119.9 nM</td><td>≈ 4.7× weaker</td></tr><tr><td>B56γ</td><td>9.2 nM</td><td>14.0 nM</td><td>≈ 1.5× weaker</td></tr></table>',
    'Moving the Aβ N-terminus onto Aα (AαChimera) barely changes binding; removing it from Aβ (AβNT) only partly restores B56α binding',
    '<b>Ours</b>: spring-attached Aα–C vs WT Aα–C against the same immobilised B56']]}

A30 = anim(30, '[Plans] Phosphatase activity assay (DiFMUP)', 'difmup.html', 'Method of Day et al. 2025: is the assembled enzyme still active?',
    'Schematic of the assay; curves computed from Day et al.’s fitted K_M and V_max (Fig. S4D)', [
    '<b>Enzyme</b>: assembled holoenzyme, 20 nM AC core with excess (200 nM) B56γ',
    '<b>Substrate DiFMUP</b> is dark; PP2A removes its phosphate → <b>DiFMU</b>, which fluoresces blue',
    'The <b>initial slope</b> of fluorescence vs time is the velocity v at that substrate concentration',
    '<b>Repeat</b> at 2–128 µM (triplicate) and fit Michaelis–Menten:' + tex(r'v = \frac{V_{\max}\,[S]}{K_M + [S]}', True),
    '<b>Day et al.</b>: Aα and Aβ give virtually the same curve. <b>Ours</b>: the same assay with and without the DNA spring'])
A30b = {'kind': 'media', 'page': '30b', 'section': 'assay', 'title': 'DiFMUP results (Day et al.): same activity for Aα and Aβ', 'figW': 760,
  'q': 'Day et al. Fig. S4D, mean ± SD of 3 replicates',
  'items': [{'img': 'native/difmup_kinetics.png', 'w': 537, 'h': 427, 'cap': 'Day et al. Fig. S4D: velocity vs [DiFMUP]'}],
  'points': [{'text': t} for t in [
    '<table class="tb"><tr><th></th><th>Aα–B56γ–C</th><th>Aβ–B56γ–C</th></tr><tr><td>K<sub>M</sub></td><td>25.6 ± 3.7 µM</td><td>34.9 ± 3.1 µM</td></tr><tr><td>V<sub>max</sub></td><td>60.1 ± 3.5 /s</td><td>56.8 ± 0.7 /s</td></tr></table>',
    'V<sub>max</sub> within 6 %: once assembled, both holoenzymes are <b>equally active</b>',
    'Yet B56 binding differs up to 4.7×: activity alone would have missed the difference, hence all three assays']]}

SLIDES = {
1: {'kind': 'cover', 'page': 1, 'title': 'Lab Meeting 261029', 'subtitle': 'Ian Poon',
    'meta': 'DNA-spring protein chimeras: sfGFP and PP2A Aα'},
'1b': {'kind': 'media', 'page': '1b', 'section': 'background', 'title': 'Why attach DNA springs to proteins?', 'figW': 1180,
       'msgChannel': 'cryo', 'q': 'Force and atomic structure, on the same molecules, at the same time',
       'items': [{'iframe': 'interactive/springs.html', 'title': 'Force methods vs structural methods vs DNA-spring chimeras',
                  'cap': 'Schematic: what each approach can measure'}],
       'points': [{'text': t, 'msg': {'step': i + 1}} for i, t in enumerate([
           '<b>Single-molecule force methods</b> (optical and magnetic tweezers, AFM) apply and measure force, but read out only one end-to-end distance, not the structure',
           '<b>Structural methods</b> (cryoEM, X-ray) give atomic structures, but of proteins at rest, with no force applied',
           '<b>Our approach</b>: attach a DNA spring to the protein (sfGFP, PP2A) across two chosen residues, so the spring itself applies the force',
           'Every molecule carries its own spring, so a whole sample is under load and can be frozen for cryoEM: <b>force and atomic structure simultaneously</b>, which the traditional single-molecule methods cannot do',
           'For each protein we follow the same plan, <b>SAM</b>: Synthesis, Assay, MD analysis (next slide)'])]},
2: text(2, 'For each protein studied, we seek “SAM”', 'background', [
    grp('<b>S</b>ynthesis goals', 'First obtain evidence, by gel, that it can be synthesized', 'Then obtain a large purified batch for analysis'),
    grp('<b>A</b>ssay goals', 'Functional assays of the protein under mechanical tension from the attached DNA springs', 'Also structural imaging assays (cryoEM etc.)'),
    grp('<b>M</b>D analysis goals', 'A suitability test: can the addition of DNA springs meaningfully perturb our protein’s structure?',
        'A computational/theoretical complement to the experimental results (e.g. theoretical force modelling of DNA-spring dynamics)')]),
3: text(3, 'sfGFP initial progress and open problems', 'status', [
    grp('<b>Synthesis</b>', 'SeArre and Micheal: purified and obtained WT and tet-mutant sfGFP',
        'Micheal Schnebly: first obtained evidence of circular product on tet-mutant sfGFP on gel just before I arrived',
        todo('obtain a batch of purified sfGFP chimera samples for further analysis')),
    grp('<b>Assay</b>', 'SeArre Abebe: initialized fluorometer. No protein sample readings taken yet',
        todo('perform fluorometry readings once purified sfGFP chimera samples are obtained')),
    grp('<b>MD analysis</b>', 'Darartu Gamada: initial talks with Doeke in 2020 for setting up MD on a ubiquitin with cysteine attachments. No other result documented in Notebook',
        todo('set up MD system for tet-TCO/BCN protein chimeras'))], badge=AGO),
4: text(4, 'How do we identify circular products?', 'synthesis', [
    grp('By weight?', 'Insufficient: linear and circular products cannot be distinguished'),
    grp('Solution:', 'Comparison of 1 tet vs 2 tet'),
    grp('Hybridization: theory', 'Linear products run differently from circular products')]),
11: text(11, 'What is MD analysis? How?', 'md', [
    grp('<b>Ramachandran plots</b>', 'Does the addition of tet groups significantly affect the structural conformation of our protein?'),
    grp('<b>Backbone rigidity: RMSF</b>', 'To what extent does the spring deform the entire scaffold? Does the overall structural integrity of the protein still persist?'),
    grp('<b>Forces applied: direct force field + modelled forces</b>', 'How much force can we actually transmit to the protein via the DNA springs?',
        'Do theoretical models of force agree with the forces measured in the MD force field?'),
    grp('<b>DNA dynamics</b>', 'Is the DNA kinked, or coiled up (under torsion), which would then affect DNA-spring force application?',
        'Is the DNA otherwise interacting with the protein, other than applying direct mechanical force?')],
    aside={'w': 560, 'html': '<b>Recap: MD goals</b><ul><li>Force modelling (to make the system reproducible and configurable to different proteins)</li>'
                             '<li>Suitability of a candidate protein for study under a DNA spring</li></ul>'}),
13: text(13, 'Ian’s sfGFP progress summary', 'status', [
    grp('<b>Synthesis</b>', done('obtain a batch of purified sfGFP chimera samples'), 'High yield without need for purification: high throughput and cheap'),
    grp('<b>Assay</b>: fluorometry', done('perform fluorometry readings once purified sfGFP chimera samples are obtained'), 'Reproduced tet quenching results',
        'Weak evidence of excitation/emission spectra change due to force applied'),
    grp('<b>MD analysis</b>', done('set up MD system for tet-TCO/BCN sfGFP protein chimeras'), 'Rigid backbone, minimal deformation under force',
        'Chromophore environment largely unchanged under force, consistent with fluorometry results')], badge=NOW),
14: text(14, 'PP2A initial progress and open problems', 'status', [
    grp('<b>Synthesis</b>', 'SeArre and Micheal: purified and obtained WT and tet-mutant pp2a', 'Micheal Schnebly: could not obtain evidence of circular product on gel yet',
        todo('obtain circular product pp2a')),
    grp('<b>Assay</b>', 'Not attempted before in lab', todo('plan out pp2a assays')),
    grp('<b>MD analysis</b>', 'Not attempted before in lab, although there is literature on MD of pp2a A-subunit isoforms and mutants; no direct mechanical tension or spring systems',
        todo('set up MD system for tet-TCO/BCN protein chimeras'))], badge=AGO),
23: A23,
'22a': PP_FN,
'22b': PP_ABC,
'22c': PP_HEAT,
'22d': PP_ISO,
24: A24,
29: A29,
'29b': A29b,
30: A30,
'30b': A30b,
32: text(32, 'Ian’s pp2a progress summary', 'status', [
    grp('<b>Synthesis</b>', done('synthesized the lab’s first pp2a chimera circular product (gel verified)'), pend('obtain a purified sample batch for further analysis')),
    grp('<b>Assay</b>', done('planned out 3 distinct assay routes for pp2a'), 'Ian: set up all necessary downstream analysis tools for cryoEM',
        'Micheal: loaded the sfGFP samples I made on the cryoEM grid, pending imaging session',
        pend('substrate binding assay once purified sample batch obtained'), pend('holoenzyme formation analysis once purified sample batch imaged with cryoEM'),
        pend('phosphatase activity assay once purified sample batch obtained')),
    grp('<b>MD analysis</b>', pend('set up MD system for tet-TCO/BCN pp2a protein chimeras'),
        '<span class="t-pend">Previous MD of pp2a in the literature shows a much more flexible structure, so it will likely respond to spring force application</span>')], badge=NOW),
33: text(33, 'Setting up cryoEM downstream analysis (new to Hekstra lab)', 'background', [
    grp('Ideally I would have explained the theory <b>alongside a live CryoSPARC demo on imaged sfGFP chimeras</b>, had we had real screening time. But the grid-screening session planned for 2 days ago was postponed: the cryoEM microscope broke down yet again 😢'),
    grp('So instead: the <b>fundamental concepts that govern a cryoEM analysis workflow</b>, which I picked up while setting up and learning CryoSPARC, new to the lab'),
    grp('A good starting point if we want to build our <b>own open-source or custom software</b>, as with careless and the other Hekstra-lab X-ray crystallography tools, since <b>CryoSPARC is paid software</b>'),
    grp('Next: from electron optics (defocus, contrast, CTF) to reconstruction and refinement')])
}

deck = json.loads((ROOT / 'deck.json').read_text())
out = []
for s in deck['slides']:
    pg = s.get('page')
    if pg in ('1b', '29b', '30b', '22a', '22b', '22c', '22d'):
        continue
    if pg in SLIDES:
        s = SLIDES[pg]
    if pg == 23:
        out += [SLIDES[x] for x in ('22a', '22b', '22c', '22d')]
    out.append(s)
    for parent, child in ((1, '1b'), (29, '29b'), (30, '30b')):
        if pg == parent:
            out.append(SLIDES[child])
deck['slides'] = out
(ROOT / 'deck.json').write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print(len(out), 'slides;', sum(1 for s in out if s['kind'] == 'image'), 'image slides left')
