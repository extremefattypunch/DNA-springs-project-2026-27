# Authoring guide for the six clips

The spec is `../../manim-gpt-animation-plan.md` (repo root). Its storyboards,
timings, equations, guardrails and acceptance tests are requirements. This file
covers the project's conventions and the helper API. Read it before writing a
scene.

## Run / inspect

```bash
cd /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/manim_cryoem_math
source env/activate_manim.sh            # Manim 0.22 + TinyTeX; media + cache on netscratch
manim -ql scenes/s02_alignment_averaging.py Scene02AlignmentAveraging      # 480p15 preview
manim -qm -s scenes/s02_alignment_averaging.py Scene02AlignmentAveraging   # 720p last frame
# videos: $MANIM_CRYOEM_SCRATCH/media/videos/<module>/<480p15|720p30|1080p30>/<Scene>.mp4
# stills: $MANIM_CRYOEM_SCRATCH/media/images/<module>/<Scene>_ManimCE_v0.22.0.png
# contact sheet (one frame every 2 s, 4x4 grid) to look at with the Read tool:
ffmpeg -loglevel error -y -i VIDEO -vf "fps=1/2,scale=480:-1,tile=4x4" -frames:v 1 sheet.png
# a single frame at time T (seconds):
ffmpeg -loglevel error -y -ss T -i VIDEO -frames:v 1 frame.png
ffprobe -v error -show_entries format=duration -of csv=p=0 VIDEO     # clip length
```

Put contact sheets and frames in your own scratch folder:
`$MANIM_CRYOEM_SCRATCH/review/<scene>/`. Never write into the repo except your
scene file. Final 1080p renders are done later by a Slurm job, not by you.

## Scene file skeleton

```python
"""Clip NN — <title>. Storyboard beats -> code sections (see spec §3)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import styling as st, demos as dm, projections as pj, fourier as fo, optics as op
from helpers.cache import cached

class SceneNN...(Scene):           # ThreeDScene for clip 03 (and only where needed)
    def construct(self):
        data = cached("alignment", dm.alignment)     # precomputed arrays; never compute per frame
        ...
```

- **Class names and modules:** use the names in `main.py`'s `SCENES`.
- **Clip-specific data:** if a clip needs extra arrays, write a module-level
  function in your scene file and wrap it with `cached("<scene>_<name>", fn, **params)`.
  Do **not** edit `helpers/`, because the other clip authors share it. If a
  helper is wrong or missing something, work around it locally and report it.
- **Determinism:** use only seeded generators (`np.random.default_rng(2027 + k)`),
  no `random` and no time-dependent behaviour.

## Visual conventions (shared across all six clips)

- **Frame:** 14.2 × 8 scene units; 1 unit ≈ 135 px at 1080p. Background `st.BG` comes from manim.cfg.
- **Colours** (`helpers/styling.py`). Never repurpose them:

  | Constant | Colour | Meaning |
  |---|---|---|
  | `st.OBS` | cyan | observation / measured image |
  | `st.TRUTH` | green | ground truth / signal |
  | `st.EST` | orange | estimate (pose, reconstruction) |
  | `st.ALT` | purple | alternate conformer / appendage state B |
  | `st.RESID` | red | residual / noise |
  | `st.TEXT`, `st.MUTED` | white, grey | text |
  | `st.BADGE` | yellow | the SIMULATION tag and "Takeaway" heading |

- **Text:** `st.T(text, size)` (Nimbus Sans) and `st.M(tex, size=...)` (MathTex).
  - Body text ≥ 28; math ≥ 36; nothing smaller than `st.SMALL_SIZE` = 26 (badges and footnotes only).
  - Keep 0.45 units of outer margin. Nothing may touch or cross the frame edge.
- **Opening:** `st.title_card(title, question)` at the top-left, with the clip's question from the spec.
- **SIMULATION badge:** `st.simulation_badge()` sits in the upper-right corner for as long as
  any synthetic image or map is on screen (practically: the whole clip). In ThreeDScene, add it
  (and all text) with `self.add_fixed_in_frame_mobjects(...)`.
- **Phantom label:** the first time the phantom appears in a clip, show `st.phantom_disclaimer()`
  ("synthetic PP2A-inspired shape; not a measured structure").
- **Shots:** one primary idea per shot. FadeOut old annotations before the next concept. No flashing,
  no fast spins, no strobing noise; noise changes only on intentional animations.
- **Ending:** a 3–5 s static freeze on `st.takeaway_card(math, text, xray)`, with one mathematical
  takeaway and an optional short X-ray analogy, while the key visual stays visible if space allows.
  The last frame must be clean and readable on its own.
- **Images:** `st.array_to_image_mobject(arr_xy, height=..., stops=st.GRAY | st.DIVERGING | st.MAGMA,
  symmetric=True for signed data, nearest=True for noisy pixel images)`.
  - Arrays are (x, y)-indexed. The helper flips them so +y is up.
  - Swap images with `FadeTransform` or by cross-fading two ImageMobjects; never rebuild an image in an updater.
- **3D density:** `st.density_image(dm.cached_render(rot, tilt, psi, state="A", wA=None), color, height,
  label_colors={15: st.ALT})`.
  - Shaded isodensity of the smooth phantom; blob index 15 is the appendage, 0–14 the scaffold domains.
  - `wA` in [0, 1] renders the mixture wA·ρ_A + (1−wA)·ρ_B.
  - An uncached render takes about 7 s, so prefer a handful of fixed views.
- **3D geometry (ThreeDScene):** `drawing.phantom_3d(state, appendage_color=st.ALT)` (translucent
  ellipsoids), `Arrow3D`, `ThreeDAxes`, `Surface`/`Polygon` planes. Keep meshes small.

## Mathematical conventions (verified by tests — keep the on-screen math consistent)

- **Object rotation:** ρ_R(r) = ρ(Rᵀr). The beam runs along z.
- **Projection:** (P_R ρ)(x, y) = ∫ ρ(Rᵀ(x, y, z)) dz.
- **Fourier slice:** FT₂[P_R ρ](k_x, k_y) = ρ̂(Rᵀ(k_x, k_y, 0)ᵀ). `tests/test_math_helpers.py` shows that Rᵀ is right and R is wrong.
- **Orientation:** R ∈ SO(3) has 3 DOF. ZYZ Euler angles: R = R_z(ψ) R_y(tilt) R_z(rot). (rot, tilt) set the
  viewing direction v = Rᵀe_z (2 DOF); ψ is the in-plane rotation (the 3rd DOF).
- **Imaging model (simplified weak-phase, parallel projection):**
  Ĩ_i(k) = CTF_i(k) · P̂_{R_i}ρ(k) · e^{−2πi k·t_i} + ε̃_i(k).
  - CTF(k) = −sin χ(k) with χ = πλΔf k² − (π/2)C_s λ³ k⁴ + arcsin(A).
  - Parameters: 300 kV (λ = 0.0197 Å), C_s = 2.7 mm, A = 0.07.
  - The CTF oscillates with zeros and sign reversals. It is not a blur.
- **Noise:** i.i.d. zero-mean Gaussian. This is an approximation; real detector and shot noise is dose-dependent.
- **SNR convention:** signal RMS (std over the box) / noise σ.
- **Shift operator:** T_t via the phase ramp e^{−2πi k·t}.
- **FSC:** standard per-shell formula, half-map threshold 0.143. Map-model FSC uses 0.5.

## Precomputed data (`helpers/demos.py`, use via `cached(name, dm.<name>)`)

All arrays are (x, y)-indexed numpy arrays. Units are Å; images are line integrals.

| Clip | Function and contents |
|---|---|
| 01 | `image_formation()` (128 px, 1.5 Å/px) — see below |
| 02 | `alignment()` (96 px, 1.6 Å/px, one view, 256 particles, raw SNR 0.3) — see below |
| 03 | `reconstruction()` (64³, 2 Å/voxel, poses KNOWN) — see below |
| 04 | `refinement()` (48³, 2.5 Å/voxel, 120 particles) — see below |
| 05 | `classification()` (64 px, 2 Å/px, 400 particles of one view, 55% state A) — see below |
| 06 | `fsc_resolution()` (128³, 1 Å/voxel) — see below |

**01 `image_formation()`**
- `clean`
- `defoci` [0.6, 1.0, 1.6, 2.4] µm; `k` (1/Å), `ctf_curves[4, 600]`, `ctf_images[4]`
- `ctf2d[4]` (fftshifted 2D CTF), `first_zero[4]` (1/Å)
- `snrs` [4, 1, .5, .25, .1, .05], `noise_sigma[6]`
- `noisy[6]`: the 1.0 µm CTF image plus one FIXED white-noise pattern scaled to each SNR
- `signal_rms`

**02 `alignment()`**
- `signal` S, `sigma`
- `phi`, `t` (true), `phi_est`, `t_est` (estimates with 1.5° / 0.4 Å jitter)
- `images[12]` (raw windows), `aligned_examples[12]`
- `ns` [1, 4, 16, 64, 256], `unaligned_mean[5]`, `aligned_mean[5]`
- `snr_measured[5]`, `snr_theory[5]` (0.30 → 4.83, ∝√N), `sigma_mean_measured`, `sigma_mean_theory`
- Rotation uses an exact Fourier-shear method, so noise stays i.i.d.

**03 `reconstruction()`**
- `truth`
- `show_R[4]` with `show_proj[4]` and `show_ft[4]` (log|FT| of each projection)
- `slice_ft2` vs `slice_central`: the numerical Fourier-slice check (128² log-magnitudes)
- `rotations[300]`
- `counts` [2, 6, 20, 80, 300] → `recons[5]`, `recon_corr` [.54 .86 .96 .998 .999], `coverage_counts`
- Biased set: `biased_rotations[300]` (tilt ≤ 12°) → `biased_recon`, `biased_corr` 0.50, `coverage_biased` vs `coverage_uniform`
- Side-on projections of every map: `side_truth`, `side_recons[5]`, `side_biased`

**04 `refinement()`**
- Toy projection matching. Each iteration is a pose search (rotation + shift by cross-correlation) followed by a
  regularised direct-Fourier volume update (Gaussian smoothness prior).
- History arrays of length 17:
  - `label` (0 = start, 1 = after pose step, 2 = after volume step)
  - `loss_data`, `loss_reg`, `loss_total` (monotone 1.27e6 → 7.74e5), `loss_noise_floor` (the data term can't go below this)
  - `pose_err` (mean angular error, 20° → 7.6°), `map_corr` vs truth (0.72 → 0.93)
- Tracked particle 0: `track_image`, `track_clean`, `track_proj[17]`, `track_resid[17]`
- `truth`, `vol_init`, `vol_final`, `reg_strength` (illustrative spatial regularisation map, 3D)

**05 `classification()`**
- `states` (0 = A, 1 = B)
- `pcs` (PCA of low-pass-filtered images) and `pcs_raw` (PCA of raw pixels: more overlap)
- `resp` (P(class 1 | particle)) and `resp_snaps[9, 400]` (EM iterations)
- `gmm_mu`, `gmm_cov`, `gmm_pi`; `ambiguous` (13 indices with 0.2 < resp < 0.8); `accuracy` 0.96
- `class_wA` [0.93, 0.003] (the fraction of state A in each class)
- `consensus_wA` 0.555, `contrast_consensus` / `contrast_classes` / `contrast_pure`
- `examples[10]` with `example_states`; `base_A`, `base_B`
- Class and consensus maps are rendered with `dm.cached_render(..., wA=w)`.

**06 `fsc_resolution()`**
- `s` (1/Å, band centres), `fsc_half` (masked half-map FSC, crosses 0.143 at 3.20 Å by construction: noise gain solved for it)
- `fsc_model` (map-model FSC, crosses 0.5 at 3.59 Å); `half_crossing_s`, `model_crossing_s`
- `fsc_unmasked` with `s_unmasked` (for reference)
- `local_res_A[16]` per blob (core 2.8–3.1, tips ~3.9, appendage 4.1 Å); `local_curves`, `blob_centers`, `flex_blur_A`
- `atoms` (toy model points), `atom_owner`
- `sec_A`, `sec_B`, `sec_truth` (central sections), `ft_A`, `ft_B`

Labels to keep honest:
- Everything is synthetic.
- The 3.2 Å crossing is engineered to echo the paper and must be marked "illustrative".
- The paper's own numbers (Table S1: global 3.4 / 3.2 Å, local 3.04–5.87 / 2.82–5.22 Å, model 3.8 / 3.6 Å) may appear only
  on a separate text-only "Paper reports (Day et al., 2025)" card, never on a synthetic curve.
