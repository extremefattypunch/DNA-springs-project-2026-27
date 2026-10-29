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
D = lambda s: tex(s, True)
RECON = [
    eqslide("c11", "Reconstruction, known poses (1/2): a linear Gaussian model",
            "If every pose (Rᵢ, tᵢ) were known, how would we recover the map V?", [
        "Start from the forward model, with each pose now <b>known</b>:"
        + D(r"Y_i(\mathbf q) = C_i(\mathbf q)\,\tilde V\big(R_i^{\mathsf T}[q_x,q_y,0]^{\mathsf T}\big)\,e^{-2\pi i\,\mathbf q\cdot\mathbf t_i} + N_i(\mathbf q)"),
        "Two Fourier rules make every step <b>linear</b> in V: a shift becomes a phase factor, and the CTF blur (a convolution) becomes a product"
        + D(r"\mathcal F[g(\mathbf x-\mathbf t)](\mathbf q) = e^{-2\pi i\,\mathbf q\cdot\mathbf t}\,\tilde g(\mathbf q), \qquad \mathcal F[h*g](\mathbf q) = \tilde h(\mathbf q)\,\tilde g(\mathbf q)"),
        "Write it in <b>real space</b> (inverse 2D transform), with the voxels of V stacked into one vector <b>v</b>:"
        + D(r"\underbrace{\mathbf y_i}_{\mathcal F_2^{-1}[Y_i]} \;=\; \underbrace{A_i\,\mathbf v}_{\mathcal F_2^{-1}\left[C_i\,\tilde V(R_i^{\mathsf T}[q_x,q_y,0]^{\mathsf T})\,e^{-2\pi i\,\mathbf q\cdot\mathbf t_i}\right]} \;+\; \underbrace{\boldsymbol\varepsilon_i}_{\mathcal F_2^{-1}[N_i]}")
        + "With the pose and CTF known, A<sub>i</sub> is a fixed matrix: the problem is <b>linear</b> in <b>v</b>",
        "The noise is Gaussian, " + T(r"\boldsymbol\varepsilon_i \sim \mathcal N(\mathbf 0, \Sigma_i)") + ", so each image is Gaussian about its prediction A<sub>i</sub><b>v</b> (the multi-pixel form of "
        + T(r"p(y\mid\mu,\sigma^2) = e^{-(y-\mu)^2/2\sigma^2}/\sqrt{2\pi\sigma^2}") + "):"
        + D(r"p(\mathbf y_i\mid\mathbf v) = \frac{1}{(2\pi)^{m/2}\,|\Sigma_i|^{1/2}}\exp\!\Big[-\tfrac12\,(\mathbf y_i-A_i\mathbf v)^{\mathsf T}\,\Sigma_i^{-1}\,(\mathbf y_i-A_i\mathbf v)\Big]"),
    ], [
        (T(r"Y_i(\mathbf q),\ \mathbf y_i"), "particle image i in Fourier space / real space (m pixels)"),
        (T(r"\mathbf q = (q_x,q_y)"), "2D spatial frequency in the image"),
        (T(r"V,\ \mathbf v"), "3D Coulomb potential (the map); <b>v</b> = its n voxels as one vector"),
        (T(r"\tilde V"), "3D Fourier transform of V"),
        (T(r"R_i,\ \mathbf t_i"), "pose of particle i: rotation in SO(3) and in-plane shift"),
        (T(r"C_i(\mathbf q)"), "CTF of image i (from its fitted defocus)"),
        (T(r"N_i,\ \boldsymbol\varepsilon_i"), "noise in Fourier / real space"),
        (T(r"\mathcal F,\ \mathcal F_2^{-1}"), "Fourier transform; inverse 2D transform"),
        (T(r"A_i"), "m × n matrix: rotate, project, shift and CTF-filter v"),
        (T(r"\Sigma_i,\ |\Sigma_i|"), "noise covariance of image i; its determinant"),
    ]),
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
               + tex(r"I(x,y) = \big|e^{i\phi(x,y)}\big|^2 = 1", True)},
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
    frame("c7", "Each defocus value gives contrast to different spatial frequencies", "ctf.html",
          "CTF(q) = sin χ(q), computed for 300 kV and Cₛ = 2.7 mm", DRAWN + "; CTF computed live", [
              "<b>Zeros</b>: the image carries no information at those q. <b>Negative lobes</b>: contrast reversed",
              "More defocus: stronger low-resolution contrast and faster oscillation",
              "Collecting over a <b>range of Δf</b> puts the zeros at different q, so every q is measured somewhere",
              "An envelope also damps high q: high resolution needs many averaged particles"]),
    clip("c8", "Forward model: what one particle image contains", "01_forward_model",
         "How do projection, the CTF and noise alter an observation?", MANIM.format(n="01"), [
             ("<b>Projection</b>: the beam integrates the 3D potential along z, so one view cannot fix depth", [0, 6.3]),
             ("<b>A beam from another direction</b>: the detector frame rotates to x′, y′, z′. The rotation <b>R</b> turns the particle so we always integrate along the beam (z′), whatever its direction", [6.3, 14.3]),
             ("ALT:In a real sample, particles sit in the ice in <b>every orientation</b>, so the beam effectively strikes each one from a different direction: each image needs its own R", None),
             ("<b>CTF</b>: every spatial frequency is scaled, removed (at a zero) or sign-reversed (negative lobe)", [14.3, 24]),
             ("<b>Noise</b>: independent and comparable to the signal here; far larger than it in real micrographs", [24, 31.5]),
             ("Together, in Fourier space:" + eq(
                 '<span class="kn"><i>Y</i><sub>i</sub>(<b>q</b>)</span> = <span class="kn"><i>C</i><sub>i</sub>(<b>q</b>)</span> '
                 '<span class="un"><i>Ṽ</i></span>(<span class="un"><i>R</i><sub>i</sub></span>ᵀ[<i>q<sub>x</sub></i>, <i>q<sub>y</sub></i>, 0]ᵀ) '
                 '<i>e</i><sup>−2π<i>i</i> <b>q</b>·<span class="un"><b>t</b><sub>i</sub></span></sup> + <span class="nz"><i>N</i><sub>i</sub>(<b>q</b>)</span>'), [31.5, 37.8]),
             ("<b>Takeaway</b>: projection loses depth, the CTF filters frequencies and noise swamps each image, so the 3D source must be recovered from many images", [37.8, 42]),
             ("<b class=\"kn\">Known</b>: the image <i>Y</i><sub>i</sub> and its CTF <i>C</i><sub>i</sub> (from the fitted Δf). "
              "<b class=\"un\">Unknown</b>: each particle’s pose (<i>R</i><sub>i</sub>, <b>t</b><sub>i</sub>) and the 3D potential <i>V</i>", None),
             ("<b>Ultimate goal</b>: the 3D potential V in real space, to fit the atomic model into", None)], overlay={'video': 'cryo/cryosparc_orientations.mp4', 'poster': 'cryo/cryosparc_orientations_poster.jpg', 'dark': True, 'cap': 'Animation from the CryoSPARC guide (guide.cryosparc.com): rays striking objects on their way to the detector'}),
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
]

for _s in SLIDES:   # slide c1 is static: open the flow charts on their last state (brackets + roadmap)
    if _s["page"] == "c1":
        _s["points"][-1]["msg"] = {"step": 4}

deck = json.loads((ROOT / "deck.json").read_text())
keep = [s for s in deck["slides"] if s.get("page") not in DROP and not str(s.get("page", "")).startswith("c")]
deck["slides"] = keep + SLIDES
(ROOT / "deck.json").write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print(f"{len(keep)} kept + {len(SLIDES)} cryoEM = {len(deck['slides'])} slides")
