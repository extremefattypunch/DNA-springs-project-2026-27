#!/usr/bin/env python3
"""Write the cryoEM downstream-analysis section into deck.json.

Replaces the four static cryoEM slides (PDF pages 25-28) and appends the
animated section at the end of the deck. Re-running is idempotent: every
slide whose page starts with "c" is rebuilt.

Sources: "cryoEM/Lab Meeting Ian 261029(cryoEM downstream) - Google Slides.pdf"
(content) and Ian's Manim clips in cryoEM/manim_cryoem_math/output/ (copied
to cryo/). Iframe slides use interactive/cryo/*.html, driven by {deck: 'cryo',
step}. Video points carry seek [start, end] in seconds of the clip (beat
boundaries from each scene's docstring and contact sheets).
"""
import html, json, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
DROP = {25, 26, 27, 28}   # What is cryoEM downstream analysis? / Steps 1-3


def tex(src, display=False):
    return f'<span class="tex{" d" if display else ""}">{html.escape(src, quote=False)}</span>'


def eq(html):
    return f'<span class="eq blk">{html}</span>'


def frame(page, title, src, q, cap, pts):
    return {"kind": "media", "page": page, "section": "cryoem", "title": title, "figW": 1180,
            "msgChannel": "cryo", "q": q,
            "items": [{"iframe": f"interactive/cryo/{src}", "title": title, "cap": cap}],
            "points": [{"text": t, "msg": {"step": i + 1}} for i, t in enumerate(pts)]}


def clip(page, title, name, q, cap, pts, overlay=None):
    s = clip0(page, title, name, q, cap, pts)
    if overlay:
        s['items'].append(dict(overlay, over=0))
        for pt in s['points']:
            if pt['text'].startswith('ALT:'):
                pt['text'] = pt['text'][4:]; pt['alt'] = True
    return s


def clip0(page, title, name, q, cap, pts):
    return {"kind": "media", "page": page, "section": "cryoem", "title": title, "figW": 1200, "q": q,
            "items": [{"video": f"cryo/{name}.mp4", "poster": f"cryo/{name}_poster.jpg", "cap": cap, "dark": True}],
            "points": [dict(text=t, **({"seek": s} if s else {})) for t, s in pts]}


MANIM = "Ian’s Manim clip {n}: mathematical illustration in dimensionless units (each click plays the next beat)"
DRAWN = "Animation drawn for this deck from the content of the cryoEM slides"

def gl(rows):
    return '<table class="gl">' + ''.join(f'<tr><td>{a}</td><td>{b}</td></tr>' for a, b in rows) + '</table>'


def eqslide(page, title, q, pts, terms):
    return {"kind": "media", "page": page, "section": "cryoem", "items": [], "final": True, "big": True, "title": title, "q": q,
            "aside": {"w": 560, "html": "<b>Terms</b>" + gl(terms)},
            "points": [{"text": t} for t in pts]}


T = lambda s: tex(s)
FM = lambda k, s: f'<span class="fm{k}">{tex(s)}</span>'   # forward-model term, coloured as in forward.html
# WPOA phasor diagram for slide c5 (static, drawn to scale with phi = 0.2 rad)
PHASOR = '<svg viewBox="0 0 540 300" width="540" height="300" style="display:block;margin:14px auto 0;overflow:visible" font-family="STIX Two Text, Georgia, serif" font-size="22"> <defs> <marker id="wpA" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10z" style="fill:var(--amber)"/></marker> <marker id="wpV" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10z" style="fill:var(--violet)"/></marker> <marker id="wpT" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10z" style="fill:var(--teal)"/></marker> </defs> <g style="stroke:var(--ink-3);fill:none" stroke-width="1.4"><path d="M26 236H506"/><path d="M46 266V116"/></g> <text x="490" y="264" style="fill:var(--ink-3)" font-style="italic">Re</text> <text x="38" y="108" style="fill:var(--ink-3)" font-style="italic">Im</text> <path d="M443.1 283.9A400 400 0 0 0 428.1 117.8" style="stroke:var(--ink-3);fill:none" stroke-width="1.6" stroke-dasharray="6 5"/> <text x="358.1" y="107.8" style="fill:var(--ink-3)" font-size="17" font-family="Archivo, sans-serif">unit circle, |z| = 1</text> <path d="M116 236A70 70 0 0 0 114.6 222.3" style="stroke:var(--teal);fill:none" stroke-width="1.6"/> <text x="124" y="232" style="fill:var(--teal)" font-style="italic" font-size="20">φ</text> <path d="M432 236V222H446" style="stroke:var(--ink-2);fill:none" stroke-width="1.4"/> <line x1="46" y1="236" x2="443" y2="236" style="stroke:var(--amber)" stroke-width="4" marker-end="url(#wpA)"/> <line x1="446" y1="236" x2="446" y2="159.0" style="stroke:var(--violet)" stroke-width="4" marker-end="url(#wpV)"/> <line x1="46" y1="236" x2="443.0" y2="156.6" style="stroke:var(--teal)" stroke-width="4" marker-end="url(#wpT)"/> <text x="246.0" y="266" text-anchor="middle" style="fill:var(--amber)"><tspan font-weight="700">1</tspan><tspan font-family="Archivo, sans-serif" font-size="17"> unscattered reference</tspan></text> <text x="458" y="178" style="fill:var(--violet)"><tspan font-style="italic">iφ</tspan><tspan x="458" dy="20" font-family="Archivo, sans-serif" font-size="16">scattered</tspan><tspan x="458" dy="18" font-family="Archivo, sans-serif" font-size="16">(+90°)</tspan></text> <text x="236.0" y="182.0" text-anchor="middle" transform="rotate(-11.31 236.0 182.0)" style="fill:var(--teal)"><tspan font-style="italic">1 + iφ ≈ e</tspan><tspan font-style="italic" dy="-9" font-size="15">iφ</tspan><tspan dy="9" font-family="Archivo, sans-serif" font-size="17"> exit wave</tspan></text> </svg>'
D = lambda s: tex(s, True)
RECON = [
    {"kind": "media", "page": "c11a", "section": "cryoem", "title": "The forward model in Fourier space: CTF, shift and noise", "figW": 1120,
     "msgChannel": "cryo", "q": "",
     "items": [{"iframe": "interactive/cryo/fmterms.html", "title": "The CTF, shift and noise terms of the forward model",
                "cap": DRAWN + "; images computed live"}],
     "points": [
         {"text": "Recall the forward model (42): every term is a function of the spatial frequency " + FM("Q", r"\mathbf q"), "msg": {"step": 1}},
         {"text": "The <b>CTF</b> " + FM("C", r"C_i(\mathbf q)") + " is a simple <b>product</b> in Fourier space", "msg": {"step": 2}},
         {"banner": True, "text": '<span class="ttl">Why Fourier space pays off</span>So now that we have expressed our 2D projection in <b>Fourier space</b>, '
                                  'we can directly <b>multiply</b> it by the <b>CTF</b> C<sub>i</sub>(q), the term that carries the defocus <b>phase shift</b> χ(q).'},
         {"text": "But how exactly do the <b>shift</b> " + FM("T", r"\mathbf t_i") + " and the <b>noise</b> " + FM("N", r"N_i(\mathbf q)") + " make sense here?", "msg": {"step": 3}},
         {"text": "Two <b>Fourier properties</b>: a <b>shift</b> becomes a phase factor, a <b>convolution</b> becomes a product", "msg": {"step": 4}},
         {"text": FM("T", r"e^{-2\pi i\,\mathbf q\cdot\mathbf t_i}") + " is just the in-plane <b>shift</b> " + FM("T", r"\mathbf t_i") + " of the particle in real space, written in Fourier space", "msg": {"step": 5}},
         {"text": "So the <b>pose</b> " + FM("R", r"(R_i,\mathbf t_i)") + " says <b>where</b> exactly each particle is (" + FM("T", r"\mathbf t_i") + ") and how it is <b>oriented</b> relative to the beam (" + FM("R", "R_i") + ")", "msg": {"step": 6}},
         {"text": '<span class="headline"><b>Noise</b> ' + FM("N", r"N_i(\mathbf q)") + " is taken to be <b>Gaussian</b>, as assumed throughout the literature</span>", "msg": {"step": 7}}]},
    dict(eqslide("c11", "Reconstruction, known poses (1/2): a linear Gaussian model",
            "If every pose (Rᵢ, tᵢ) were known, how would we recover the map V?", [
        "With every pose <b>known</b>, the forward model in <b>real space</b> (the voxels of V stacked into one vector <b>v</b>) is <b>linear</b>: "
        + T(r"\mathbf y_i = A_i\,\mathbf v + \boldsymbol\varepsilon_i") + ", where A<sub>i</sub> (rotate, project, shift, CTF-filter) is a fixed matrix",
        "<b>The residual</b> " + T(r"\mathbf y_i - A_i\mathbf v") + " is the recorded image minus its <b>prediction</b> " + T(r"A_i\mathbf v") + ". If <b>v</b> is the true map, the residual is exactly the noise:"
        + D(r"\mathbf y_i - A_i\mathbf v = \boldsymbol\varepsilon_i, \qquad \boldsymbol\varepsilon_i \sim \mathcal N(\mathbf 0, \Sigma_i)"),
        "<b>One pixel</b>: the noise " + T(r"\varepsilon \sim \mathcal N(0,\sigma^2)") + " has the standard Gaussian density. A measurement " + T(r"y = \mu + \varepsilon") + " is then as probable as its noise, so evaluate that density at the residual " + T(r"y - \mu") + ":"
        + D(r"p(\varepsilon) = \frac{1}{\sqrt{2\pi\sigma^2}}\,e^{-\varepsilon^2/2\sigma^2} \quad\Longrightarrow\quad p(y\mid\mu) = \frac{1}{\sqrt{2\pi\sigma^2}}\,e^{-(y-\mu)^2/2\sigma^2}"),
        "<b>m pixels</b>, independent noise " + T(r"\Sigma = \mathrm{diag}(\sigma_1^2,\dots,\sigma_m^2)") + ": multiply the m one-pixel densities. The normalisers multiply to " + T(r"(2\pi)^{m/2}|\Sigma|^{1/2}") + " and the exponents add:"
        + D(r"p(\mathbf y\mid\boldsymbol\mu) = \prod_{j=1}^{m}\frac{e^{-(y_j-\mu_j)^2/2\sigma_j^2}}{\sqrt{2\pi\sigma_j^2}} = \frac{1}{(2\pi)^{m/2}|\Sigma|^{1/2}}\exp\Big[-\tfrac12\sum_{j=1}^{m}\frac{(y_j-\mu_j)^2}{\sigma_j^2}\Big]")
        + "and the sum is " + T(r"(\mathbf y-\boldsymbol\mu)^{\mathsf T}\Sigma^{-1}(\mathbf y-\boldsymbol\mu)") + ", which also covers correlated pixels (any Σ)",
        "<b>Set μ = A<sub>i</sub>v</b>, the prediction, to get the likelihood of image i given the map:"
        + D(r"\boxed{\;p(\mathbf y_i\mid\mathbf v) = \frac{1}{(2\pi)^{m/2}|\Sigma_i|^{1/2}}\exp\Big[-\tfrac12\,(\mathbf y_i-A_i\mathbf v)^{\mathsf T}\,\Sigma_i^{-1}\,(\mathbf y_i-A_i\mathbf v)\Big]\;}")
        + "The exponent is the <b>noise-weighted squared residual</b>, and <b>v</b> enters only through it: the closer the prediction is to the image, the more probable the image",
    ], [
        (T(r"\mathbf y_i,\ Y_i"), "particle image i in real space (m pixels) / Fourier space"),
        (T(r"V,\ \mathbf v"), "3D Coulomb potential (the map); <b>v</b> = its n voxels as one vector"),
        (T(r"A_i"), "m × n matrix: rotate, project, shift and CTF-filter v"),
        (T(r"A_i\mathbf v"), "the <b>prediction</b>: image i without noise"),
        (T(r"\mathbf y_i - A_i\mathbf v"), "the <b>residual</b>: what the prediction does not explain"),
        (T(r"\boldsymbol\varepsilon_i,\ N_i"), "noise in real / Fourier space"),
        (T(r"\mathcal N(\boldsymbol\mu,\Sigma)"), "Gaussian with mean <b>μ</b> and covariance Σ"),
        (T(r"\mu,\ \sigma^2"), "mean and variance of one pixel’s Gaussian"),
        (T(r"m,\ j"), "number of pixels in an image; pixel index"),
        (T(r"\Sigma_i,\ |\Sigma_i|,\ \Sigma_i^{-1}"), "noise covariance of image i; its determinant; its inverse"),
        (T(r"p(\mathbf y_i\mid\mathbf v)"), "likelihood: probability of image i if the map were <b>v</b>"),
    ]), big=False),
    eqslide("c11b", "Reconstruction, known poses (2/2): the normal equations",
            "Maximum likelihood: which map v makes all N images most probable?", [
        "The N images are independent given <b>v</b>, so the joint likelihood is a product:"
        + D(r"p(\mathbf y_1,\dots,\mathbf y_N\mid\mathbf v) = \prod_{i=1}^{N} p(\mathbf y_i\mid\mathbf v)"),
        "A product of thousands of tiny probabilities is numerically unstable; the log turns it into a sum, "
        + T(r"\log\prod_i p_i = \sum_i \log\, p_i") + ". Keeping only the terms that depend on <b>v</b>:"
        + D(r"-\log\, p(\mathbf y_1,\dots,\mathbf y_N\mid\mathbf v) = L(\mathbf v) + \text{const}, \qquad L(\mathbf v) = \frac12\sum_{i=1}^{N}(\mathbf y_i-A_i\mathbf v)^{\mathsf T}\,\Sigma_i^{-1}\,(\mathbf y_i-A_i\mathbf v)")
        + "Maximizing the likelihood = minimizing L, a noise-weighted least-squares misfit",
        "Take the gradient with respect to <b>v</b>:"
        + D(r"\nabla_{\mathbf v} L = \sum_{i=1}^{N} A_i^{\mathsf T}\,\Sigma_i^{-1}\,(A_i\mathbf v-\mathbf y_i)"),
        "Setting " + T(r"\nabla_{\mathbf v} L = \mathbf 0") + " gives the <b>normal equations</b>, which are solved for the estimate " + T(r"\hat{\mathbf v}") + ":"
        + D(r"\boxed{\;\Big(\sum_{i=1}^{N} A_i^{\mathsf T}\,\Sigma_i^{-1}A_i\Big)\,\hat{\mathbf v} \;=\; \sum_{i=1}^{N} A_i^{\mathsf T}\,\Sigma_i^{-1}\,\mathbf y_i\;}")
        + "Right: every image <b>back-projected</b> into 3D, weighted by its noise. Left: how often each part of <b>v</b> was measured. Solving divides one by the other",
    ], [
        (T(r"N"), "number of particle images"),
        (T(r"p(\mathbf y_i\mid\mathbf v)"), "likelihood: probability of image i if the map were <b>v</b>"),
        (T(r"L(\mathbf v)"), "negative log-likelihood (the misfit to minimize)"),
        (T(r"\Sigma_i^{-1}"), "precision matrix: down-weights noisy pixels"),
        (T(r"A_i^{\mathsf T}"), "back-projection: smears a 2D image back into 3D"),
        (T(r"\nabla_{\mathbf v}"), "gradient with respect to the voxels of <b>v</b>"),
        (T(r"\hat{\mathbf v}"), "maximum-likelihood estimate of the map"),
        (T(r"\text{const}"), "terms independent of <b>v</b> (e.g. the normalizer of each Gaussian)"),
    ]),
    eqslide("c12", "Reconstruction, unknown poses: alternate poses and map",
            "Neither the poses nor the map is known: alternate between them until both stop changing", [
        "Now the map <b>v</b> and every pose " + T(r"\phi_i = (R_i,\mathbf t_i)") + " are unknown. The same misfit, with " + T(r"A_i = A(\phi_i)") + ", is minimized over <b>both</b>:"
        + D(r"\min_{\mathbf v,\ \phi_1,\dots,\phi_N}\; L(\mathbf v,\phi_1,\dots,\phi_N) = \frac12\sum_{i=1}^{N}\big\|\mathbf y_i - A(\phi_i)\,\mathbf v\big\|^2_{\Sigma_i^{-1}}")
        + "A(φ) is nonlinear in the pose, so there is no closed-form solution",
        "Start from an <b>initial volume</b> " + T(r"\mathbf v^{(0)}") + " (ab initio, or a low-pass filtered reference), then repeat for k = 0, 1, 2, …",
        "<b>Pose step</b>: hold the map, give each image its best pose"
        + D(r"\phi_i^{(k)} = \arg\min_{\phi\,\in\,SO(3)\times\mathbb R^2}\ \big\|\mathbf y_i - A(\phi)\,\mathbf v^{(k)}\big\|^2_{\Sigma_i^{-1}}"),
        "<b>Map step</b>: hold the poses and solve the known-pose problem (normal equations)"
        + D(r"\Big(\sum_{i} A(\phi_i^{(k)})^{\mathsf T}\Sigma_i^{-1}A(\phi_i^{(k)})\Big)\,\mathbf v^{(k+1)} = \sum_{i} A(\phi_i^{(k)})^{\mathsf T}\Sigma_i^{-1}\,\mathbf y_i"),
        "Stop when L no longer decreases. Neither step can raise L, so the loop converges, but to a <b>local</b> optimum that depends on " + T(r"\mathbf v^{(0)}")
        + ". RELION and CryoSPARC soften the pose step: they weight every pose by its probability (expectation–maximization) instead of keeping only the best one",
    ], [
        (T(r"\phi_i = (R_i,\mathbf t_i)"), "pose of particle i, now an unknown"),
        (T(r"A(\phi)"), "the projection matrix for pose φ (same as A<sub>i</sub>)"),
        (T(r"SO(3)\times\mathbb R^2"), "all 3D rotations × all in-plane shifts"),
        (T(r"\|\mathbf r\|^2_{\Sigma^{-1}}"), T(r"\mathbf r^{\mathsf T}\Sigma^{-1}\mathbf r") + ": noise-weighted squared misfit"),
        (T(r"\mathbf v^{(0)}"), "initial volume"),
        (T(r"k"), "iteration number"),
        (T(r"\mathbf v^{(k)},\ \phi_i^{(k)}"), "map and poses after iteration k"),
        (T(r"\arg\min"), "the value that minimizes the expression"),
    ]),
]


SLIDES = [
    dict(frame("c1", "cryoEM downstream analysis: how it differs from X-ray crystallography", "pipeline.html", "",
          DRAWN, [
              '<table class="tb cmp"><tr><th></th><th>X-ray crystallography</th><th>cryoEM</th></tr>'
              '<tr><td>Radiation</td><td>Diffracted X-ray photons</td><td>Scattered and transmitted electron waves</td></tr>'
              '<tr><td>3D reconstruction</td><td>Recover missing Fourier phases → <b>electron-density map</b></td><td>Recover particle orientations → <b>electrostatic-potential map</b></td></tr>'
              '<tr><td>Refinement</td><td><b>Reciprocal space</b>: adjust atoms to match the measured diffraction intensities</td><td><b>Real space</b>: adjust atoms to fit the reconstructed 3D map</td></tr></table>'
              '<p class="after">Both fit an <b>atomic model</b> into their reconstructed 3D maps.</p>']), final=True),
    frame("c2", "What does it mean to defocus a lens?", "lens.html",
          "Thin lens: 1/f = 1/s + 1/s′, with the object and the sensor fixed", DRAWN, [
              "<b>Stronger lens</b> (shorter f): the rays cross before the sensor and spread into a disc",
              "<b>Weaker lens</b> (longer f): the rays would meet beyond the sensor; again a disc. In the TEM this is <b>underfocus</b>",
              "The objective lens current sets f, so turning the focus knob moves the focus by Δf. cryoEM images are taken deliberately ≈ 0.5–3 µm out of focus"]),
    frame("c3", "Effect 1: the protein shifts the phase AND absorbs amplitude (intensity)", "wave.html",
          "What happens to the electron beam as it passes through protein?", DRAWN + " (phase shift exaggerated)", [
              "An <b>amplitude object</b> absorbs electrons: the wave leaves weaker, and the image is darker there",
              "A <b>phase object</b> keeps the amplitude but delays the wave by φ (amplitude decrease is minor; the <b>phase shift</b> dominates)",
              "A protein in ice is almost a pure phase object, and the detector records |ψ|² = |e<sup>iφ</sup>|² = 1: <b>in focus it is nearly invisible</b>"]),
    frame("c4", "Effect 2: the protein scatters the wave (diffraction)", "scatter.html",
          "Where does each spatial frequency go?", DRAWN, [
              "A repeat of spacing d scatters the wave by α ≈ λ/d (Bragg: 2d sin θ = λ, with α = 2θ)",
              "Finer detail scatters further: the direction encodes the <b>spatial frequency q = 1/d</b>",
              "The objective lens sends each direction to one point of its back focal plane: fine detail lands far from the centre",
              "At 300 kV (λ = 0.0197 Å) even 3 Å detail scatters by only ≈ 0.4°"]),
    {"kind": "media", "page": "c5", "section": "cryoem", "items": [], "final": True, "big": True,
     "title": "Key assumption: the weak phase object approximation",
     "q": "Implication: the scattered wave is always 90° out of phase with the unscattered wave",
     "aside": {"w": 600, "html": "<b>Also implies</b> the amplitude of the exit wave after scattering at the protein remains at unity, so there is <b>no contrast</b> between the unscattered reference and the scattered wave: the detector only picks up <b>intensity</b>, not the phase difference!"
               + tex(r"I(x,y) = \big|e^{i\phi(x,y)}\big|^2 = 1", True)
               + PHASOR + "<div style=\"text-align:center;font-size:16px;color:var(--ink-3);margin-top:2px\">Drawn to scale for φ = 0.2 rad</div>"
               + "<div style=\"margin-top:14px\">Intensity of the <b>exit wave</b> (what the detector records):</div>"
               + tex(r"I_{\mathrm{exit}} = |1 + i\phi|^2 = (1 + i\phi)(1 - i\phi) = 1 + \phi^2 = 1.04", True)
               + "<div style=\"margin-top:8px;padding:10px 14px;border-left:5px solid var(--amber);border-radius:6px;background:color-mix(in srgb, var(--amber) 13%, transparent);color:var(--ink)\">"
                 "<b>I<sub>exit</sub> ≈ 1 = I<sub>unscattered</sub></b>: the exit wave is as bright as the unscattered reference, so the protein gives <b>no contrast</b></div>"},
     "points": [{"text": t} for t in [
         "Incident plane wave:" + tex(r"\Psi_{\mathrm{in}}(\mathbf r, z) = e^{ikz}, \qquad k = \frac{2\pi}{\lambda}", True),
         "Exit wave after the protein:" + tex(r"\Psi_{\mathrm{exit}}(x,y) = e^{i\phi(x,y)}", True)
         + "Under the WPOA, " + tex(r"|\phi| \ll 1") + ", so a Taylor expansion gives" + tex(r"e^{i\phi} \approx 1 + i\phi", True),
         "Therefore," + tex(r"\Psi_{\mathrm{exit}}(x,y) \approx 1 + i\phi(x,y)", True) + "where <b>1</b> is the <b>unscattered reference</b> and <b>iφ(x, y)</b> is the <b>scattered wave</b>",
         "Because " + tex(r"i = e^{i\pi/2}") + ", this perturbation is in <b>quadrature</b> (90° out of phase) with the reference, when φ is real"]] + [{"banner": True, "text": '<span class="ttl">Under this approximation</span><b>1.</b> The scattered wave always differs in phase from the unscattered wave by <b>π/2</b>.<br><br><b>2.</b> Their sum, the exit wave, is <b>virtually indistinguishable in intensity</b> from a wave that never passed through a protein, so we <b>cannot see the protein</b> yet.'}]},
    dict(frame("c6", "Key principle: defocus turns phase shifts into contrast", "contrast.html",
          "How do we make the scattered wave visible?", DRAWN + " (φ exaggerated)", [
              "Extra phase χ = π/2: the scattered wave becomes <b>antiparallel</b> to the reference, so the particle is <b>dark</b>",
              "χ = 3π/2: <b>parallel</b>, so the particle is <b>bright</b> (contrast reversed)",
              "Contrast ≈ −2φ sin χ: <b>sin χ is the contrast transfer function</b>, largest at χ = π/2 + kπ",
              "Defocus supplies χ(q) = πλΔf q² − ½πC<sub>s</sub>λ³q⁴: a different extra phase for every spatial frequency q"]), intro='<b>Key idea</b>: adjusting the <b>defocus</b> changes the <b>relative phase shift</b> between the scattered and unscattered waves.<br><br>That is what creates the <b>image contrast</b> we need to see our protein.'),
    eqslide("c6b", "Where the defocus phase comes from: the Fresnel propagator",
            "Why does defocus give each spatial frequency q its own extra phase?", [
        "A <b>3D plane wave</b> of wavelength λ travels along its <b>wavevector</b> " + T(r"\mathbf k") + ". Specimen detail of <b>spatial frequency</b> "
        + T(r"\mathbf q") + " scatters it off the optic axis z, giving it a transverse part " + T(r"\mathbf k_\perp")
        + D(r"\psi(\mathbf r) = e^{\,i\mathbf k\cdot\mathbf r} = e^{\,i(k_x x + k_y y + k_z z)}, \qquad \mathbf k_\perp = (k_x, k_y) = 2\pi\mathbf q"),
        "<b>Elastic</b> scattering keeps the energy, so " + T(r"|\mathbf k| = k_0") + " is unchanged. Solve for " + T(r"k_z") + " and expand for small angles (" + T(r"\lambda q \ll 1") + "):"
        + D(r"k_z^2 + |\mathbf k_\perp|^2 = k_0^2, \qquad k_0 = \frac{2\pi}{\lambda}")
        + D(r"k_z = \sqrt{k_0^2 - (2\pi q)^2} = k_0\sqrt{1 - \lambda^2 q^2} \;\simeq\; k_0 - \tfrac12 k_0\lambda^2 q^2 = k_0 - \pi\lambda q^2"),
        "Over a distance z the component gains " + T(r"e^{\,i k_z z}") + ". Dropping the carrier phase " + T(r"e^{\,i k_0 z}") + " shared by every component:"
        + D(r"\boxed{H_z(q) = e^{-i\pi\lambda z q^2}}"),
        "<b>Defocus</b> records the wave Δf from the in-focus plane (underfocus: z = −Δf), so every q gets its own extra phase, growing as q²:"
        + D(r"\chi_{\text{def}}(q) = \pi\lambda\,\Delta f\, q^2")
        + "Spherical aberration adds " + T(r"-\tfrac12\pi C_s\lambda^3 q^4") + ", giving the χ(q) of the CTF (next slide)",
    ], [
        (T(r"\psi(\mathbf r)"), "electron wavefunction (complex amplitude) at position " + T(r"\mathbf r = (x, y, z)")),
        (T(r"z;\ x, y"), "optic axis (beam direction); transverse image-plane coordinates"),
        (T(r"\mathbf k = (k_x, k_y, k_z)"), "wavevector: points along the direction of travel"),
        (T(r"\mathbf k_\perp,\ k_z"), "transverse (tilt off the axis) and axial components of " + T(r"\mathbf k")),
        (T(r"k_0 = 2\pi/\lambda"), "wavenumber: length of " + T(r"\mathbf k") + ", fixed by the beam energy"),
        (T(r"\lambda"), "electron wavelength (0.0197 Å at 300 kV, so λq ≪ 1 for q ≲ 0.5 Å⁻¹)"),
        (T(r"\mathbf q = (q_x, q_y),\ q"), "2D spatial frequency of a specimen detail (Å⁻¹) and its magnitude"),
        (T(r"H_z(q)"), "Fresnel propagator: phase picked up over a distance z"),
        (T(r"\Delta f"), "defocus (positive = underfocus); sign conventions differ between texts"),
        (T(r"\chi_{\text{def}}(q)"), "extra phase from defocus"),
        (T(r"C_s"), "spherical-aberration coefficient: adds " + T(r"-\tfrac12\pi C_s\lambda^3 q^4")),
    ]),
    frame("c7", "Each defocus value gives contrast to different spatial frequencies", "ctf.html",
          "CTF(q) = sin χ(q), computed for 300 kV and Cₛ = 2.7 mm", DRAWN + "; CTF computed live", [
              "<b>Zeros</b>: the image carries no information at those q. <b>Negative lobes</b>: contrast reversed",
              "More defocus: stronger low-resolution contrast and faster oscillation",
              "Collecting over a <b>range of Δf</b> puts the zeros at different q, so every q is measured somewhere",
              "An envelope also damps high q: high resolution needs many averaged particles"]),
    dict(frame("c7b", "The forward model: from the 3D molecule to each recorded image", "forward.html", "",
          DRAWN + "; projection, CTF, shift and noise computed live", [
              FM("V", r"V(\mathbf r)") + ": the molecule’s 3D <b>scattering potential</b>, at each position " + FM("V", r"\mathbf r"),
              "Particle " + T("i") + " lies in the ice with its own <b>orientation</b> " + FM("R", "R_i"),
              "<b>Projection</b> (41): the beam adds up " + FM("V", "V") + " along " + FM("P", "z") + ", giving the 2D image " + FM("P", r"P_{R_i}V"),
              "<b>Fourier space</b> (42): that image is a <b>central slice</b> " + FM("V", r"\tilde V(R_i^{\mathsf T}\mathbf q)") + " of the 3D transform",
              "The <b>CTF</b> " + FM("C", r"C_i(\mathbf q)") + " multiplies every spatial frequency " + FM("Q", r"\mathbf q"),
              "A <b>shift</b> " + FM("T", r"\mathbf t_i") + " in the box is a <b>phase ramp</b> " + FM("T", r"e^{-2\pi i\,\mathbf q\cdot\mathbf t_i}"),
              "<b>Noise</b> " + FM("N", r"N_i(\mathbf q)") + " on top gives the recorded image " + FM("Y", r"Y_i(\mathbf q)"),
              "So (42) models the <b>projection of the 3D volume</b> (scattering potential) onto the <b>2D images</b> we record",
              "<b>Known</b>: " + FM("Y", "Y_i") + ", the Fourier transform of each raw 2D particle image, and " + FM("C", "C_i") + ", the CTF calculated at each q from χ(q)",
              '<span class="headline"><b>Goal</b>: find the <b>pose</b> ' + FM("R", r"(R_i,\mathbf t_i)") + " of every particle and the <b>volume</b> " + FM("V", "V") + ", the electron scattering potential into which we fit <b>atomic models</b></span>"]),
         figW=1120),
    clip("c8", "Projection and orientation: what one particle image records", "01_forward_model",
         "What does the beam record, and why does each image need its own orientation R?", MANIM.format(n="01"), [
             ("<b>Projection</b>: the beam integrates the 3D potential along z, so one view cannot fix depth", [0, 6.3]),
             ("<b>A beam from another direction</b>: the detector frame rotates to x′, y′, z′. The rotation <b>R</b> turns the particle so we always integrate along the beam (z′), whatever its direction", [6.3, 14.3]),
             ("ALT:In a real sample, particles sit in the ice in <b>every orientation</b>, so the beam effectively strikes each one from a different direction: each image needs its own R", None)], overlay={'video': 'cryo/cryosparc_orientations.mp4', 'poster': 'cryo/cryosparc_orientations_poster.jpg', 'dark': True, 'cap': 'Animation from the CryoSPARC guide (guide.cryosparc.com): rays striking objects on their way to the detector'}),
    {"kind": "media", "page": "c10", "section": "cryoem", "items": [], "title": "Central slice theorem",
     "q": "Each 2D Fourier image samples a central plane through the 3D Fourier volume", "points": [{"text": t} for t in [
         "<b>Projection.</b> The image integrates the scattering potential V along the beam direction z. R is a rotation matrix, so changing R changes the direction the beam passes through the molecule"
         + tex(r"(P_R V)(x,y) = \int_{-\infty}^{\infty} V\big(R^{\mathsf T}(x,y,z)^{\mathsf T}\big)\,dz", True),
         "<b>The 3D Fourier transform</b> of the potential:"
         + tex(r"\tilde V(k_x,k_y,k_z) = \iiint V(x,y,z)\,e^{-2\pi i\,(k_x x + k_y y + k_z z)}\,dx\,dy\,dz", True),
         "<b>The 2D Fourier transform of a projection</b> (R = I, beam along z):"
         + tex(r"\mathcal F_2\{PV\}(q_x,q_y) = \iint \Big[\int V(x,y,z)\,dz\Big]\, e^{-2\pi i\,(q_x x + q_y y)}\,dx\,dy", True),
         "The exponential does not depend on z, so it moves inside the z-integral: this is the 3D transform with k<sub>z</sub> = 0"
         + tex(r"\mathcal F_2\{PV\}(q_x,q_y) = \iiint V(x,y,z)\,e^{-2\pi i\,(q_x x + q_y y + 0\cdot z)}\,dx\,dy\,dz = \tilde V(q_x,q_y,0)", True),
         "So the 2D transform of a projection is the <b>plane through the origin</b> of the 3D transform, perpendicular to the beam. For any orientation R the plane simply rotates:"
         + tex(r"\mathcal F_2\{P_R V\}(\mathbf q) = \tilde V\big(R^{\mathsf T}(q_x,q_y,0)^{\mathsf T}\big)", True),
         "<b>Consequence</b>: every particle image fills one central plane of V̂. Images from many orientations fill the 3D Fourier volume, and an inverse 3D transform gives back V (next slide); orientations clustered around one axis leave a missing cone unmeasured"]]},
] + RECON + [
    frame("c13", "Reconstruction, unknown poses: the loop in action", "align.html",
          "Can the poses and the map be found together from noisy images alone?", DRAWN + "; 100 simulated images, alignment computed live", [
              "<b>2D toy</b>: every image " + T(r"\mathbf y_i") + " is the same object rotated by an unknown angle " + T(r"\phi_i") + ", plus heavy noise. Start from a blurred blob " + T(r"\mathbf v^{(0)}"),
              "<b>Pose step</b>: compare each image with " + T(r"\mathbf v^{(0)}") + " rotated through every angle (5° steps); the best match is its pose. The images turn to their estimated poses",
              "<b>Map step</b>: rotate the images back and average them → " + T(r"\mathbf v^{(1)}") + ". Already far clearer than any single image",
              "<b>Repeat</b>: a better map gives better poses, and better poses a better map. The mean pose error falls from 22° to ≈ 2° within two iterations",
              "<b>Compare with the truth</b>: the same object, up to an overall rotation that the images alone can never fix. A poor " + T(r"\mathbf v^{(0)}") + " can trap the loop in a local optimum"]),
    {"kind": "media", "page": "c14", "section": "cryoem", "title": "In CryoSPARC: from particle images to a 3D map",
     "layout": "row", "aspect": 1.86, "rowH": 540,
     "q": "", "introTtl": "In practice",
     "intro": "Altogether, in <b>CryoSPARC</b> it looks something like this: a past project of mine, recovering <b>3D maps</b> from particle images of "
              "<b>actin filaments</b>, essentially a proxy for <b>force-transducing DNA strands</b>.",
     "items": [{"img": "cryo/cryosparc_j1_import.png", "w": 1903, "h": 1021, "cap": "J1, Import Particle Stack: 5,000 particle images, with no orientations yet"},
               {"img": "cryo/cryosparc_j2_abinitio.png", "w": 1961, "h": 1029, "cap": "J2, Ab-initio Reconstruction (1 class, C1): the 3D map taking shape"}],
     "points": [
         {"text": "<b>J1</b> imports 5,000 particles (160 × 160 px at 5.76 Å/px). Its warnings say the file carries <b>no orientations</b>: the poses (R<sub>i</sub>, <b>t</b><sub>i</sub>) are exactly what must be found. It also imports each particle’s <b>CTF parameters</b> (defocus, amplitude contrast 0.1), so J2 can apply C<sub>i</sub>(q) during reconstruction", "hl": 0},
         {"text": "<b>J2</b>, ab-initio reconstruction, alternates pose and map steps from a random start (the loop of the previous slides). By iteration 900 a filament density emerges at ≈ 17 Å", "hl": 1},
         {"text": "Run on the <b>Hekstra lab cluster</b> from my home directory; the raw data are on <b>netscratch</b>"}]},
]

for _s in SLIDES:   # slide c1 is static: open the flow charts on their last state (brackets + roadmap)
    if _s["page"] == "c1":
        _s["points"][-1]["msg"] = {"step": 4}
        _s["points"].append({"banner": True, "text": '<span class="ttl">Key differences</span>'   # Ian, 2026-11-05
            "Both X-ray crystallography and cryo-EM measure <b>intensities</b> rather than phases directly. However, unlike X-ray crystallography, "
            "where Fourier <b>phase information is lost</b> from diffraction intensity measurements, cryo-EM exploits <b>electron-wave interference</b> "
            "to encode phase information into detectable <b>image intensity contrast</b>.<br><br>"
            "Furthermore, whereas proteins in X-ray crystals form an <b>ordered lattice</b>, individual proteins in cryo-EM are <b>randomly oriented</b> "
            "in vitreous ice. Consequently, cryo-EM requires computational estimation of <b>particle orientations</b> to reconstruct a 3D structural map, "
            "whereas X-ray crystallography primarily requires recovering <b>missing Fourier phases</b>."})
    if _s["page"] == "c5":    # opening banner (Ian, 2026-11-05): the two questions the WPOA answers
        _s["introTtl"] = "Two questions"
        _s["intro"] = ("We saw that when the beam passes through the protein its <b>phase shifts</b>, yet the protein is then "
                       "<b>not visible</b> on a detector that measures <b>intensity</b>.<br><br>"
                       "<b>(1)</b> Why is that so?<br>"
                       "<b>(2)</b> How much is the phase of each beam that passed through the protein (<b>scattered</b>) shifted, relative to the <b>unscattered</b> beam?<br><br>"
                       "It turns out the analysis simplifies greatly once we know that this phase shift is <b>SMALL</b>.")
    if _s["page"] == "c8":    # opening banner (Ian, 2026-11-05)
        _s["intro"] = ("So what is the <b>projection</b> " + tex(r"P_{R}V") + ", and why do we include the <b>orientation</b> "
                       + tex("R") + " in our model?")
    if _s["page"] == "c11":   # opening banner (Ian, 2026-11-05): the real-space form of the forward model, in terms of its Fourier-space pieces
        _s["introTtl"] = "Reconstruction in practice"
        _s["intro"] = ("So how do we recover the 3D volume V (i.e. <b>reconstruction</b>) in practice? First, let us assume we know one of the two unknowns: "
                       "the <b>pose</b> (R<sub>i</sub>, <b>t</b><sub>i</sub>).<br><br>"
                       "Expressing the forward equation in <b>real space</b> (inverse 2D transform " + tex(r"\mathcal F_2^{-1}") + ", with the voxels of V stacked into one vector <b>v</b>), we have:"
                       + tex(r"\underbrace{\mathbf y_i}_{\mathcal F_2^{-1}[Y_i]} \;=\; \underbrace{A_i\,\mathbf v}_{\mathcal F_2^{-1}\left[C_i\,\tilde V(R_i^{\mathsf T}[q_x,q_y,0]^{\mathsf T})\,e^{-2\pi i\,\mathbf q\cdot\mathbf t_i}\right]} \;+\; \underbrace{\boldsymbol\varepsilon_i}_{\mathcal F_2^{-1}[N_i]}", True)
                       + "With the pose and CTF known, A<sub>i</sub> is a fixed matrix: the problem is <b>linear</b> in <b>v</b>.")
    if _s["page"] == "c6":    # Ian, 2026-11-05: why only parallel / antiparallel scattering gives first-order contrast
        _row = lambda lab, src, col: (f'<tr><td style="padding:3px 14px 3px 0;white-space:nowrap;color:{col}"><b style="color:{col}">{lab}</b></td>'
                                      f'<td style="padding:3px 0">{tex(src)}</td></tr>')
        _s["points"].append({"msg": {"step": 4}, "text":
            '<div style="padding:12px 16px;border-left:5px solid var(--amber);border-radius:8px;background:color-mix(in srgb, var(--amber) 10%, transparent);color:var(--ink)">'
            "<b>Exit-wave intensity</b>, scattered wave along the reference or at 90°:"
            '<table style="border-collapse:collapse;margin:8px 0 8px;font-size:.95em">'
            + _row("antiparallel", r"(1-\phi)^2 = 1 - 2\phi + \phi^2 \approx 1 - 2\phi", "var(--sec-status)")
            + _row("parallel", r"(1+\phi)^2 = 1 + 2\phi + \phi^2 \approx 1 + 2\phi", "var(--sec-synthesis)")
            + _row("perpendicular", r"|1+i\phi|^2 = 1 + \phi^2 \approx 1", "var(--ink-3)")
            + "</table>Under the WPOA " + tex(r"\phi \ll 1") + ", so to first order " + tex(r"\phi^2 \approx 0") + " while " + tex(r"2\phi") + " is not: "
            "only the parallel and antiparallel cases move the intensity away from unity, giving <b>contrast</b>. At 90° (in focus) the change is second order: <b>invisible</b></div>"})
    if _s["page"] == "c7b":   # opening banner (Ian, 2026-11-05)
        _s["introTtl"] = "Next: the forward model"
        _s["intro"] = ("With an understanding of the basics of <b>image formation</b> in cryoEM, we now move on to the key theoretical model "
                       "that governs <b>3D structure determination</b> from the <b>2D images</b> we obtain from the microscope.")
    if _s["page"] == "c2":    # opening banner (Ian, 2026-11-05): what the cryoEM theory section is for
        _s["introTtl"] = "From 2D images to a 3D map"
        _s["intro"] = ("I will now try to explain the <b>CryoSPARC downstream analysis</b> in a hopefully <b>intuitive</b> way, as that is the software "
                       "I have mainly been working with: essentially, <b>how we get from 2D images to a 3D density map</b>.<br><br>"
                       "Crystallographers can hopefully draw parallels with <b>reconstruction in X-ray crystallography</b>, as some concepts are common "
                       "to the two workflows: <b>diffraction</b>, and its connection to <b>Fourier space</b>.")
    if _s["page"] == "c10":   # opening banner: why the projection is now written in Fourier space
        _s["introTtl"] = "Why Fourier space?"
        _s["intro"] = ("Here we want to express our projected 2D image in <b>reciprocal (Fourier) space</b>.<br><br>"
                       "This is because, as we saw, <b>CTF(q)</b> adds a different <b>relative phase shift</b> to the components at each <b>spatial frequency q</b>.")

deck = json.loads((ROOT / "deck.json").read_text())
keep = [s for s in deck["slides"] if s.get("page") not in DROP and not str(s.get("page", "")).startswith("c")]
deck["slides"] = keep + SLIDES
(ROOT / "deck.json").write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print(f"{len(keep)} kept + {len(SLIDES)} cryoEM = {len(deck['slides'])} slides")
