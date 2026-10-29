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


def clip(page, title, name, q, cap, pts):
    return {"kind": "media", "page": page, "section": "cryoem", "title": title, "figW": 1200, "q": q,
            "items": [{"video": f"cryo/{name}.mp4", "poster": f"cryo/{name}_poster.jpg", "cap": cap, "dark": True}],
            "points": [dict(text=t, **({"seek": s} if s else {})) for t, s in pts]}


MANIM = "Ian’s Manim clip {n}: mathematical illustration in dimensionless units (each click plays the next beat)"
DRAWN = "Animation drawn for this deck from the content of the cryoEM slides"

SLIDES = [
    frame("c1", "cryoEM downstream analysis: how it differs from X-ray crystallography", "pipeline.html",
          "What is cryoEM downstream analysis?",
          DRAWN + "; diffraction image and micrographs from the original slide (not PP2A data)", [
              "<b>X-ray</b>: the detector records only spot intensities, so the phases must be recovered (<b>phasing</b>) before any map exists",
              "<b>cryoEM</b>: every particle image is a real-space projection with its phases intact, but the particle’s <b>pose is unknown</b> and its signal is buried in noise",
              "<b>CryoSPARC, TOPAZ</b>: movies → 3D map. <b>RELION, EMReady, COOT</b>: map → atomic model",
              "The next slides explain how the 3D map is made: defocus and the CTF, the forward model, the Fourier slice theorem, and reconstruction with known and unknown poses"]),
    frame("c2", "What does it mean to defocus a lens?", "lens.html",
          "Thin lens: 1/f = 1/s + 1/s′, with the object and the sensor fixed", DRAWN, [
              "<b>Stronger lens</b> (shorter f): the rays cross before the sensor and spread into a disc",
              "<b>Weaker lens</b> (longer f): the rays would meet beyond the sensor; again a disc. In the TEM this is <b>underfocus</b>",
              "The objective lens current sets f, so turning the focus knob moves the focus by Δf. cryoEM images are taken deliberately ≈ 0.5–3 µm out of focus"]),
    frame("c3", "Effect 1: the protein shifts the phase of the electron wave", "wave.html",
          "What happens to the electron beam as it passes through protein?", DRAWN + " (phase shift exaggerated)", [
              "An <b>amplitude object</b> absorbs electrons: the wave leaves weaker, and the image is darker there",
              "A <b>phase object</b> keeps the amplitude but delays the wave by φ (amplitude decrease is minor; the <b>phase shift</b> dominates)",
              "A protein in ice is almost a pure phase object, and the detector records |ψ|² = |e<sup>iφ</sup>|² = 1: <b>in focus it is nearly invisible</b>"]),
    frame("c4", "Effect 2: the protein scatters the wave (diffraction)", "scatter.html",
          "Where does each spatial frequency go?", DRAWN, [
              "A repeat of spacing d scatters the wave by α ≈ λ/d (Bragg: 2d sin θ = λ, with α = 2θ)",
              "Finer detail scatters further: the direction encodes the <b>spatial frequency q = 1/d</b>",
              "The objective lens sends each direction to one point of its back focal plane: a <b>Fourier transform</b> of the object",
              "At 300 kV (λ = 0.0197 Å) even 3 Å detail scatters by only ≈ 0.4°"]),
    frame("c5", "Key assumption: the weak phase object approximation", "wpoa.html",
          "Why is an in-focus protein invisible?", DRAWN + " (φ exaggerated)", [
              "Exit wave Ψ<sub>exit</sub>(x, y) = e<sup>iφ(x,y)</sup>: the phasor turns by φ but keeps length 1",
              "For |φ| ≪ 1: e<sup>iφ</sup> ≈ 1 + iφ = <b>unscattered reference</b> + <b>scattered wave</b>, 90° apart (i = e<sup>iπ/2</sup>)",
              "Intensity |1 + iφ|² = 1 + φ² ≈ 1: <b>no first-order contrast</b>"]),
    frame("c6", "Key principle: defocus turns phase shifts into contrast", "contrast.html",
          "How do we make the scattered wave visible?", DRAWN + " (φ exaggerated)", [
              "Extra phase χ = π/2: the scattered wave becomes <b>antiparallel</b> to the reference, so the particle is <b>dark</b>",
              "χ = 3π/2: <b>parallel</b>, so the particle is <b>bright</b> (contrast reversed)",
              "Contrast ≈ −2φ sin χ: <b>sin χ is the contrast transfer function</b>, largest at χ = π/2 + kπ",
              "Defocus supplies χ(q) = πλΔf q² − ½πC<sub>s</sub>λ³q⁴: a different extra phase for every spatial frequency q"]),
    frame("c7", "Each defocus value gives contrast to different spatial frequencies", "ctf.html",
          "CTF(q) = sin χ(q), computed for 300 kV and Cₛ = 2.7 mm", DRAWN + "; CTF computed live", [
              "<b>Zeros</b>: the image carries no information at those q. <b>Negative lobes</b>: contrast reversed",
              "More defocus: stronger low-resolution contrast and faster oscillation. Thon rings reveal each micrograph’s Δf (CTF estimation)",
              "Collecting over a <b>range of Δf</b> puts the zeros at different q, so every q is measured somewhere",
              "An envelope also damps high q: high resolution needs many averaged particles"]),
    clip("c8", "Forward model: what one particle image contains", "01_forward_model",
         "How do projection, the CTF and noise alter an observation?", MANIM.format(n="01"), [
             ("<b>Projection</b>: the beam integrates the 3D potential along z, so one view cannot fix depth", [0, 6.3]),
             ("<b>CTF</b>: every spatial frequency is scaled, removed (at a zero) or sign-reversed (negative lobe)", [6.3, 16]),
             ("<b>Noise</b>: independent and comparable to the signal here; far larger than it in real micrographs", [16, 23.5]),
             ("Together, in Fourier space:" + eq(
                 '<span class="kn"><i>Y</i><sub>i</sub>(<b>q</b>)</span> = <span class="kn"><i>C</i><sub>i</sub>(<b>q</b>)</span> '
                 '<span class="un"><i>Ṽ</i></span>(<span class="un"><i>R</i><sub>i</sub></span>ᵀ[<i>q<sub>x</sub></i>, <i>q<sub>y</sub></i>, 0]ᵀ) '
                 '<i>e</i><sup>−2π<i>i</i> <b>q</b>·<span class="un"><b>t</b><sub>i</sub></span></sup> + <span class="nz"><i>N</i><sub>i</sub>(<b>q</b>)</span>'), [23.5, 34]),
             ("<b class=\"kn\">Known</b>: the image <i>Y</i><sub>i</sub> and its CTF <i>C</i><sub>i</sub> (from the fitted Δf). "
              "<b class=\"un\">Unknown</b>: each particle’s pose (<i>R</i><sub>i</sub>, <b>t</b><sub>i</sub>) and the 3D potential <i>V</i>", None),
             ("Two Fourier rules keep it simple: a <b>shift</b> becomes a phase ramp e<sup>−2πi q·t</sup>, and a <b>blur</b> (convolution) becomes a product, so the CTF just multiplies", None),
             ("<b>Ultimate goal</b>: the 3D potential V in real space, to fit the atomic model into", None)]),
    {"kind": "media", "page": "c10", "section": "cryoem", "items": [], "title": "Central slice theorem",
     "q": "Each 2D Fourier image samples a central plane through the 3D Fourier volume", "points": [{"text": t} for t in [
         "<b>Projection.</b> The image integrates the scattering potential V along the beam direction z. R is a rotation matrix, so changing R changes the direction the beam passes through the molecule"
         + tex(r"(P_R V)(x,y) = \int_{-\infty}^{\infty} V\big(R^{\mathsf T}(x,y,z)^{\mathsf T}\big)\,dz", True),
         "<b>The 3D Fourier transform</b> of the potential:"
         + tex(r"\hat V(k_x,k_y,k_z) = \iiint V(x,y,z)\,e^{-2\pi i\,(k_x x + k_y y + k_z z)}\,dx\,dy\,dz", True),
         "<b>The 2D Fourier transform of a projection</b> (R = I, beam along z):"
         + tex(r"\mathcal F_2\{PV\}(q_x,q_y) = \iint \Big[\int V(x,y,z)\,dz\Big]\, e^{-2\pi i\,(q_x x + q_y y)}\,dx\,dy", True),
         "The exponential does not depend on z, so it moves inside the z-integral: this is the 3D transform with k<sub>z</sub> = 0"
         + tex(r"\mathcal F_2\{PV\}(q_x,q_y) = \iiint V(x,y,z)\,e^{-2\pi i\,(q_x x + q_y y + 0\cdot z)}\,dx\,dy\,dz = \hat V(q_x,q_y,0)", True),
         "So the 2D transform of a projection is the <b>plane through the origin</b> of the 3D transform, perpendicular to the beam. For any orientation R the plane simply rotates:"
         + tex(r"\mathcal F_2\{P_R V\}(\mathbf q) = \hat V\big(R^{\mathsf T}(q_x,q_y,0)^{\mathsf T}\big)", True),
         "<b>Consequence</b>: every particle image fills one central plane of V̂. Images from many orientations fill the 3D Fourier volume, and an inverse 3D transform gives back V (next slide); orientations clustered around one axis leave a missing cone unmeasured"]]},
    frame("c11", "Reconstruction scheme 1: you know the poses", "recon.html",
          "Least squares: L(v) = ½ Σᵢ (yᵢ − Aᵢv)ᵀ Σᵢ⁻¹ (yᵢ − Aᵢv); ∇L = 0 gives the normal equations",
          "Live 2D reconstruction (64 × 64, nearest-neighbour gridding, noisy images) computed in the browser", [
              "One image = one central line (a plane in 3D). Backprojecting it gives only streaks",
              "More views fan out from the origin: low frequencies are measured by every image, high ones by few",
              "The <b>backprojection</b> (right-hand side) over-counts the centre, so it is blurred; dividing by the <b>count</b> (left-hand matrix) sharpens it",
              "Noise averages down where a frequency is measured often; the edges stay noisy (hence weighting and regularisation)",
              "<b>Preferred orientation</b> leaves a missing wedge (a cone in 3D) that no amount of averaging fills"]),
    clip("c12", "Reconstruction scheme 2: you don’t know any poses", "04_iterative_refinement",
         "What does iterative reconstruction optimise?", MANIM.format(n="04"), [
             ("Compare each image with the projection predicted from the current map; the residual measures the mismatch", [0, 9]),
             ("<b>Pose step</b>: hold the map and search each particle’s orientation and shift to shrink its residual", [9, 19]),
             ("<b>Map step</b>: hold the poses and update the map (scheme 1). Alternate until the loss stops falling; one particle here stays in a local minimum", [19, 31]),
             ("<b>Ab initio</b> (no starting map): CryoSPARC combines branch-and-bound pose search with stochastic gradient descent on random mini-batches ℬ:" +
              eq('<i>V</i><sup>(<i>t</i>+1)</sup> = <i>V</i><sup>(<i>t</i>)</sup> − η<sub><i>t</i></sub> (1/|ℬ|) Σ<sub><i>i</i>∈ℬ</sub> ∇<sub><i>V</i></sub> ℓ<sub><i>i</i></sub>(<i>V</i><sup>(<i>t</i>)</sup>)') +
              "Cheap, noisy steps; no guarantee of the global optimum", None),
             ("<b>Regularisation</b>: a smaller data misfit is not a better map; without a penalty the fit chases the noise", [31, 39]),
             ("<b>Risk</b>: aligning noisy images to a reference that already contains the same noise reinforces it. The standard defence is to refine two independent half-sets and compare them", [39, 44.9])]),
]

deck = json.loads((ROOT / "deck.json").read_text())
keep = [s for s in deck["slides"] if s.get("page") not in DROP and not str(s.get("page", "")).startswith("c")]
deck["slides"] = keep + SLIDES
(ROOT / "deck.json").write_text(json.dumps(deck, indent=1, ensure_ascii=False))
print(f"{len(keep)} kept + {len(SLIDES)} cryoEM = {len(deck['slides'])} slides")
