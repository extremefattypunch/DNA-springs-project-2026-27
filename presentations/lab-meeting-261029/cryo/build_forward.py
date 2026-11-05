#!/usr/bin/env python3
"""Write interactive/cryo/forward.html: the forward-model explainer (slide c7b).

The equations (41) and (42) of the imaging-theory chapter are typeset to MathML with the
deck's own KaTeX (tools/tex.cjs). Every term is wrapped in \\textcolor{#0000NN}, and the
marker colour is swapped for CSS classes, so the frame can colour, glow and measure each
term as the deck steps through it. Any $...$ in the card, title and label text below is
typeset the same way.
"""
import json, pathlib, re, subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent
MARK = {}


def c(cls, src):
    """Wrap src in a marker colour that becomes class="cls" in the MathML."""
    if cls not in MARK:
        MARK[cls] = '#%06x' % (len(MARK) + 1)
    return r'\textcolor{%s}{%s}' % (MARK[cls], src)


q = r'\mathbf q'
EQ41 = (c('g41a cP', '(' + c('cP', 'P_{' + c('cR', 'R_i') + '}') + c('cV', 'V') + ')(' + c('cX', 'x,y') + ')')
        + ' = ' + c('g41b cP', r'\int ' + c('cV', r'V\big(' + c('cR', 'R_i^{-1}') + '(' + c('cX', 'x,y,z') + r')\big)')
                    + r'\,\mathrm d' + c('cX', 'z')))
EQ42 = (c('g42Y cY', 'Y_i(' + c('cQ', q) + ')') + ' = '
        + c('g42C cC', 'C_i(' + c('cQ', q) + ')') + r'\,'
        + c('g42V cVt', c('uV', r'\tilde V') + r'\Big(' + c('cR2 uR', r'R_i^{\mathsf T}') + c('cQ', r'(q_x,q_y,0)^{\mathsf T}') + r'\Big)') + r'\,'
        + c('g42E cT', r'e^{-2\pi i\,' + c('cQ', q) + r'\cdot' + c('uT', r'\mathbf t_i') + '}')
        + ' + ' + c('g42N cN', 'N_i(' + c('cQ', q) + ')'))

# card per deck step: heading, colour key, rows (term TeX, colour key, definition with $..$), or free text
CARD = [
    dict(h='The forward model', key='', text='How does the microscope turn the 3D molecule into each 2D particle image we record? '
         'Eq. (41) is one <b>projection</b> in real space; eq. (42) is the full model in <b>Fourier space</b>, with the microscope (CTF), '
         'the particle’s shift and the noise. Each click defines the next symbols, in the colour used in the equations and the panels.'),
    dict(h='1 · The molecule', key='V', rows=[
        (r'V(\mathbf r)', 'V', 'the 3D <b>specimen potential</b>: the electrostatic (Coulomb) potential that scatters the electrons, at every point of the molecule. '
         'It is the 3D <b>map</b> we want'),
        (r'\mathbf r = (x,y,z)', 'V', 'a position in 3D')]),
    dict(h='2 · Its orientation in the ice', key='R', rows=[
        (r'i', '', '<b>particle index</b>: one particle image each, $i = 1,\\dots,N$ (typically $10^5$–$10^6$ particles)'),
        (r'R_i', 'R', '3 × 3 <b>rotation matrix</b>: the orientation particle $i$ happens to have in the ice'),
        (r'R_i^{-1}', 'R', 'the inverse rotation: takes a point $(x,y,z)$ of the microscope frame back into the molecule’s own frame, where $V$ is defined')]),
    dict(h='3 · Projection along the beam: eq. (41)', key='P', rows=[
        (r'(x,y),\ z', 'X', 'detector (image) coordinates, and the <b>beam direction</b> (optic axis)'),
        (r'\int \cdots\,\mathrm d z', 'P', 'adds up $V$ along each straight electron path through the particle'),
        (r'(P_{R_i}V)(x,y)', 'P', 'the <b>projection</b>: the 2D image of the rotated particle, one value per detector pixel')]),
    dict(h='4 · The same image in Fourier space: eq. (42)', key='Q', rows=[
        (r'\mathbf q = (q_x,q_y)', 'Q', '2D <b>spatial frequency</b> in the image (Å⁻¹): detail of size $d$ sits at $|\\mathbf q| = 1/d$'),
        (r'\tilde V', 'V', 'the 3D <b>Fourier transform</b> of $V$'),
        (r'(q_x,q_y,0)^{\mathsf T}', 'Q', 'the image’s frequencies: the plane $q_z = 0$, perpendicular to the beam ($^{\\mathsf T}$: transpose, a column vector)'),
        (r'R_i^{\mathsf T} = R_i^{-1}', 'R', 'tilts that plane into the molecule’s frame: $\\tilde V$ on it is a <b>central slice</b> (two slides on)')]),
    dict(h='5 · The microscope: contrast transfer function', key='C', rows=[
        (r'C_i(\mathbf q)', 'C', 'the <b>CTF</b> of image $i$: $\\sin\\chi(q)$, with $\\chi(q) = \\pi\\lambda\\Delta f_i\\,q^2 - \\tfrac12\\pi C_s\\lambda^3 q^4$ (previous slides). '
         'It <b>multiplies</b> each frequency: at its zeros the information is lost, in negative lobes the contrast reverses. In real space: fringes around the particle')]),
    dict(h='6 · Where the particle sits in its box', key='T', rows=[
        (r'\mathbf t_i = (t_x,t_y)', 'T', 'in-plane <b>translation</b>: particle $i$ is not exactly centred in its box'),
        (r'e^{-2\pi i\,\mathbf q\cdot\mathbf t_i}', 'T', 'a shift by $\\mathbf t_i$ in real space is a <b>phase ramp</b> in Fourier space; the magnitude is unchanged. '
         '$\\mathbf q\\cdot\\mathbf t_i = q_x t_x + q_y t_y$'),
        (r'i\ \text{in the exponent}', '', 'the imaginary unit $\\sqrt{-1}$, not the particle index')]),
    dict(h='7 · Noise, and what we record', key='N', rows=[
        (r'N_i(\mathbf q)', 'N', '<b>noise</b>, added on top: mostly shot noise from the low electron dose, plus the ice (real images are far noisier than drawn)'),
        (r'Y_i(\mathbf q)', 'Y', 'the 2D Fourier transform of the <b>recorded particle image</b> $i$')]),
    dict(h='Eq. (42) in one sentence', key='', text='It models the <b>projection of the 3D volume</b> $V$ (the scattering potential) onto every <b>2D image</b> we record: '
         'rotate ($R_i$), project, filter by the CTF ($C_i$), shift ($\\mathbf t_i$) and add noise ($N_i$).'),
    dict(h='What we know', key='K', rows=[
        (r'Y_i(\mathbf q)\ \checkmark', 'Y', 'the Fourier transform of each experimentally recorded <b>raw 2D particle image</b>'),
        (r'C_i(\mathbf q)\ \checkmark', 'C', 'the CTF, <b>calculated</b> at every frequency $q$ from the $\\chi(q)$ formula we just derived, with the defocus $\\Delta f_i$ fitted to each micrograph')]),
    dict(h='The goal: what we do not know', key='U', rows=[
        (r'(R_i,\ \mathbf t_i)\ ?', 'R', 'the <b>pose</b> of every particle: its orientation and its shift'),
        (r'V\ ?', 'V', 'the 3D <b>electron scattering potential</b>: the map into which we fit <b>atomic models</b> to get the 3D structure')]),
]
TOP = [r'the molecule $V(\mathbf r)$', r'projection $P_{R_i}V$', 'with the CTF', r'shifted by $\mathbf t_i$', 'recorded image']
BOT = [r'$\tilde V$ and the slice plane', r'central slice of $\tilde V$', r'$\times\,C_i$: rings and zeros', r'phase ramp, same $|\cdot|$', r'$Y_i(\mathbf q)$']
OPS = [r'$R_i$, then $\int \mathrm d z$', r'$\times\,C_i(\mathbf q)$', r'$\times\,e^{-2\pi i\,\mathbf q\cdot\mathbf t_i}$', r'$+\,N_i(\mathbf q)$']
LAB = {  # brace labels under the equation terms: (step shown, text by step)
    'g41a': 'projection image', 'g41b': 'V summed along each beam path',
    'g42Y': 'recorded image', 'g42C': 'CTF', 'g42V': 'central slice of Ṽ', 'g42E': 'shift by tᵢ', 'g42N': 'noise'}

# ---- typeset everything in one KaTeX call ----
jobs, slots = [], []
def need(src, display=False):
    slots.append(len(jobs)); jobs.append([src, display]); return len(jobs) - 1
EQ = [need(EQ41, True), need(EQ42, True)]
INL = {}
def inline(txt):
    def sub(m):
        k = m.group(1)
        if k not in INL: INL[k] = need(k)
        return '{{%d}}' % INL[k]
    return re.sub(r'\$(.+?)\$', sub, txt)
for cd in CARD:
    if 'text' in cd: cd['text'] = inline(cd['text'])
    cd['rows'] = [(need(t), k, inline(d)) for t, k, d in cd.get('rows', [])]
TOP = [inline(t) for t in TOP]; BOT = [inline(t) for t in BOT]; OPS = [inline(t) for t in OPS]
out = json.loads(subprocess.run(['node', str(ROOT / 'tools/tex.cjs')], input=json.dumps(jobs), capture_output=True, text=True, check=True).stdout)
inv = {v: k for k, v in MARK.items()}
def classes(ml):
    ml = re.sub(r'<annotation encoding="application/x-tex">.*?</annotation>', '', ml, flags=re.S)
    return re.sub(r'mathcolor="(#[0-9a-f]{6})"', lambda m: 'class="%s"' % inv[m.group(1)], ml)
ML = [classes(m) for m in out]
fill = lambda s: re.sub(r'\{\{(\d+)\}\}', lambda m: ML[int(m.group(1))], s)
for cd in CARD:
    if 'text' in cd: cd['text'] = fill(cd['text'])
    cd['rows'] = [(ML[t], k, fill(d)) for t, k, d in cd['rows']]
DATA = dict(eq41=ML[EQ[0]], eq42=ML[EQ[1]], card=CARD, top=[fill(t) for t in TOP], bot=[fill(t) for t in BOT],
            ops=[fill(t) for t in OPS], lab=LAB)

src = (ROOT / 'cryo/forward_template.html').read_text()
(ROOT / 'interactive/cryo/forward.html').write_text(src.replace('/*DATA*/null', json.dumps(DATA, ensure_ascii=False)))
print('forward.html written:', len(MARK), 'term classes,', len(jobs), 'formulas')
