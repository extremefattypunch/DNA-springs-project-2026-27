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
23: text(23, '[PENDING] PP2A assays', 'assay', [
    grp('Phosphatase activity assay'), grp('Substrate binding assay'), grp('Holoenzyme formation analysis')],
    aside={'w': 640, 'html': '<b>Why not just a direct phosphatase assay?</b><ul><li>An artificial substrate is not equivalent to the in-vivo substrate</li>'
                             '<li>It assumes proper assembly of the C catalytic unit with the A subunit, which is not always the case when the A subunit is perturbed</li></ul>'}),
24: text(24, '[Plans] Holoenzyme formation analysis', 'assay', [
    grp('<b>Structural deformation</b>', 'Align the relatively conserved C subunit, then measure centre-of-mass distances between HEAT repeats on the A subunit, WT vs spring attached'),
    grp('<b>Ascertain spring attachment</b>', 'Run downstream analysis tools: determine the resolution level, use neural-network particle picking, refine and resolve to obtain the atomic model')],
    q='cryoEM analysis'),
29: text(29, '[Plans] Substrate binding assay', 'assay', [
    grp('Assemble the spring-attached A subunit with the catalytic subunit C. Attach the B regulatory subunit to the biosensor AR2G'),
    grp('Binding and the dissociation constant:' + tex(r'AC + B \;\xrightleftharpoons[k_{\mathrm{off}}]{\;k_{\mathrm{on}}\;}\; ABC, \qquad K_D = \frac{k_{\mathrm{off}}}{k_{\mathrm{on}}}', True)),
    grp('<table class="tb"><tr><th>Regulatory subunit</th><th>AαC K<sub>D</sub></th><th>AβC K<sub>D</sub></th><th>Aβ relative to Aα</th></tr>'
        '<tr><td>B56α</td><td>25.7 nM</td><td>119.9 nM</td><td>≈ 4.7× weaker</td></tr><tr><td>B56γ</td><td>9.2 nM</td><td>14.0 nM</td><td>≈ 1.5× weaker</td></tr></table>')],
    q='BLI (biolayer interferometry) assay: measure the dissociation constant'),
30: text(30, '[Plans] Phosphatase activity assay', 'assay', [
    grp('Fitted with the standard Michaelis–Menten kinetics equation:' + tex(r'v = \frac{V_{\max}\,[S]}{K_M + [S]}', True)),
    grp('<table class="tb"><tr><th>Measurement</th><th>Aα–B56γ–C</th><th>Aβ–B56γ–C</th></tr>'
        '<tr><td>K<sub>M</sub></td><td>25.6 ± 3.7 µM</td><td>34.9 ± 3.1 µM</td></tr><tr><td>V<sub>max</sub></td><td>60.1 ± 3.5 fluorescence/s</td><td>56.8 ± 0.7 fluorescence/s</td></tr></table>')],
    q='Fluorogenic phosphatase substrate DiFMUP', figW=760,
    items=[{'img': 'native/difmup_kinetics.png', 'w': 537, 'h': 427, 'cap': 'Velocity vs [substrate], replicates'}]),
32: text(32, 'Ian’s pp2a progress summary', 'status', [
    grp('<b>Synthesis</b>', done('synthesized the lab’s first pp2a chimera circular product (gel verified)'), pend('obtain a purified sample batch for further analysis')),
    grp('<b>Assay</b>', done('planned out 3 distinct assay routes for pp2a'), 'Ian: set up all necessary downstream analysis tools for cryoEM',
        'Micheal: loaded the sfGFP samples I made on the cryoEM grid, pending imaging session',
        pend('substrate binding assay once purified sample batch obtained'), pend('holoenzyme formation analysis once purified sample batch imaged with cryoEM'),
        pend('phosphatase activity assay once purified sample batch obtained')),
    grp('<b>MD analysis</b>', pend('set up MD system for tet-TCO/BCN pp2a protein chimeras'),
        '<span class="t-pend">Previous MD of pp2a in the literature shows a much more flexible structure, so it will likely respond to spring force application</span>')], badge=NOW),
33: text(33, 'Technical tutorials: brief walkthrough', 'background', [
    grp('New experimental techniques to the lab that I troubleshot independently over the last 4 months. Possibly useful for the lab:',
        'Setting up the MD system for DNA-spring chimeras', 'Fluorometry data collection', 'Performing downstream cryoEM analysis with CryoSPARC'),
    grp('Additional resources (to be documented via GitHub):', 'Claude skills for optimal FAIRSHARE usage and CryoSPARC usage', 'Set-up instructions for CryoSPARC')]),
}

deck = json.loads((ROOT / 'deck.json').read_text())
out = []
for s in deck['slides']:
    pg = s.get('page')
    if pg == '1b':
        continue
    if pg in SLIDES:
        s = SLIDES[pg]
    out.append(s)
    if pg == 1:
        out.append(SLIDES['1b'])
deck['slides'] = out
(ROOT / 'deck.json').write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print(len(out), 'slides;', sum(1 for s in out if s['kind'] == 'image'), 'image slides left')
