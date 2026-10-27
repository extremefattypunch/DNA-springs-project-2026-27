# Authoring guide for the six clips

The spec is `../../manim-gpt-animation-plan.md` (repo root). Its restrictions (§0),
storyboards (§3), Manim instructions (§4) and Definition of Done (§7) are
requirements. This file gives the project conventions and the helper API.
Read both before writing a scene.

## Run / inspect

```bash
cd /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/manim_cryoem_math
source env/activate_manim.sh            # Manim 0.22 + TinyTeX; media + cache on netscratch
M=$MANIM_CRYOEM_SCRATCH/media_<your scene module>      # your OWN media dir (avoids Tex-cache races)
manim render -ql --media_dir $M scenes/s02_alignment_averaging.py Scene02AlignmentAveraging   # 480p15
manim render -qm -s --media_dir $M scenes/s02_alignment_averaging.py Scene02AlignmentAveraging # 720p last frame
# video: $M/videos/<module>/480p15/<Scene>.mp4     still: $M/images/<module>/<Scene>_ManimCE_v0.22.0.png
ffmpeg -loglevel error -y -i VIDEO -vf "fps=1/2,scale=480:-1,tile=4x4" -frames:v 1 sheet.png   # contact sheet
ffmpeg -loglevel error -y -ss T -i VIDEO -frames:v 1 frame.png                                 # frame at T s
ffprobe -v error -show_entries format=duration -of csv=p=0 VIDEO                               # length
python -m pytest -m "not render" -q      # math + scope tests (must stay green)
```

Always run manim from the project root (so `manim.cfg` applies: 1920×1080, 30 fps,
background `#0E1116`). Put sheets and frames in `$MANIM_CRYOEM_SCRATCH/review/<module>/`.
Never write into the repo except your own scene file. Final 1080p renders run later
on Slurm.

## Scene file skeleton

```python
"""Clip NN — <title>. Storyboard beats -> code sections (spec §3)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne, projections as pj
from helpers import point_constellation as pc, noise_and_ctf as nc, fourier as fo
from helpers.cache import cached


class SceneNN...(Scene):        # ThreeDScene only where 3D geometry is needed (clips 01, 03, 04?)
    def construct(self):
        d = cached("alignment", ne.alignment)    # precomputed arrays; never compute per frame
        self.add(st.illustration_label())        # on screen for the whole clip
        ...
```

- Class and module names: exactly as in `main.py` `SCENES`.
- Clip-specific numbers: write a module-level function in your scene file and wrap it
  with `cached("<module>_<name>", fn, **params)`. Do **not** edit `helpers/` (shared by
  all authors). If a helper is wrong or missing something, work around it locally and
  report it.
- Determinism: only seeded generators (`nc.rng(k)` = `default_rng(2027 + k)`); no
  `random`, no time-dependent behaviour.
- Every number shown on screen (SNR, loss, error, FSC, d, probabilities) must come
  from the data dicts or from a local computation on them — no invented values.

## Hard content rules (spec §0; enforced partly by tests/test_visual_constraints.py)

- The source is **only dots** (`Dot3D` in 3D, `Dot` in 2D). Never draw a surface, blob
  outline, mesh, sphere, bonds or lines connecting dots. Gaussian kernels exist only
  inside numerical images/maps (tiny heatmaps) and curves.
- No biological/sample/dataset vocabulary anywhere, including comments: no
  protein, molecule, atom, domain, chain, structure names, papers, software names,
  "micrograph", "Å"/angstrom. Units are dimensionless ("units", "1/units", "a.u.").
  Say "observation" or "image", not "particle".
- Each clip shows `st.illustration_label()` (upper-right, whole clip). Mark simplified
  equations with `st.simplified_note(...)`. Never present toy curves as experimental results.
- Allowed objects: dots, arrows, axes, lines, planes (`Polygon`/`Square`), circles/arcs,
  plots, MathTex/Text, and small heatmaps from `st.heatmap(array)`.

## Visual conventions

- Frame 14.2 × 8 units (1 unit ≈ 135 px at 1080p). Keep 0.45 units of outer margin
  (`st.MARGIN`); nothing touches the frame edge.
- **Colours** (`helpers/style.py`), never repurposed:

  | Constant | Colour | Meaning |
  |---|---|---|
  | `st.OBS` | cyan | observed data |
  | `st.TRUTH` | green | underlying signal / ground truth (the true dots) |
  | `st.EST` | orange | estimates and alignment / orientation-update arrows |
  | `st.ALT` | purple | alternate mathematical state (the moved dot 5 in clip 05) |
  | `st.RESID` | red | noise and residuals |
  | `st.GRID`, `st.MUTED` | greys | axes, grids, secondary text |
  | `st.LABEL` | yellow | the MATHEMATICAL ILLUSTRATION tag and "Takeaway" heading |

  Clip 06's two half-maps are cyan (A) and orange (B), as the spec storyboard says.
- **Text:** `st.T(text, size)` (Nimbus Sans) and `st.M(tex, size=...)`. Body ≥ 28,
  math ≥ 36, nothing below `st.SMALL_SIZE` = 26 except the tag.
- **Opening:** `st.title_card(title, question)` top-left (title auto-shrinks to stay clear
  of the tag; break long questions with "\n").
- **Takeaway:** end with a 3–5 s static freeze containing `st.takeaway_card(lines, math=...)`;
  pass short lines (≤ ~60 characters). The card raises if lines are too wide.
- One idea per shot, at most one complex equation on screen at a time; equations sit
  next to the quantity they describe. FadeOut old annotations before the next concept.
  No flashing, no fast spins, no strobing noise; stable camera (move arrows/dots, not
  the camera) unless a slow camera move is essential.
- Target durations (spec §1): 30 / 40 / 70 / 40 / 45 / 45 s, ±15 %. Clarity beats timing.

## Lessons from the previous build of this project (Manim 0.22)

- `ImageMobject` cannot go in a `VGroup`; use `Group(...)` for anything mixing images.
- Swapping two labels/equations at the same place: use
  `Succession(FadeOut(a, shift=0.1*UP), FadeIn(b, shift=0.1*UP))`, not simultaneous fades
  (they ghost on top of each other).
- Arrows/lines drawn over noisy heatmaps need a dark outline:
  `mob.set_stroke(st.BG, width=8, background=True)`.
- ThreeDScene: text, plots and heatmaps go in `self.add_fixed_in_frame_mobjects(...)`.
  Moving `camera.frame_center` also moves fixed-in-frame overlays, so offset 3D content by
  shifting the geometry itself. After `set_camera_orientation`, call
  `self.renderer.camera.reset_rotation_matrix()` if you need the rotation matrix before
  the first frame.
- Labels over scatter clouds: `add_background_rectangle(color=st.BG, opacity=1, buff=0.05)`
  and `set_z_index(2)`.
- Check text widths against the frame (`mob.get_right()[0] <= 7.11 - 0.45`), especially
  next to the tag (its left edge is about x = 3.0).

## Helper API

### `point_constellation as pc`
`pc.POINTS` (5×3), `pc.points("A"|"B")` (B moves dot 5 = index `pc.MOVED` = 4 to
`pc.POINT_B5`), `pc.AMPLITUDES`, `pc.SIGMAS` (kernel weights, numerical use only),
`pc.PATTERN_2D` (4×2: x,y of dots 1–4; clip 02's 2D pattern), `pc.field(state, n, extent)`
(3D array [x,y,z]), `pc.grid_coords(n, extent)`, `pc.spacing(n, extent)`,
`pc.IMAGE_N`=64, `pc.IMAGE_EXTENT`=2.0 (images span [-2,2)), `pc.FIELD_N`=64,
`pc.FIELD_EXTENT`=3.0.

### `projections as pj`
`pj.rotation_from_euler(alpha, theta, psi)` (ZYZ: R = Rz(ψ)Ry(θ)Rz(α), degrees),
`pj.viewing_direction(R)` = Rᵀe_z, `pj.project_points(R, X, t)` → (m,2) = (R x)_xy + t,
`pj.transform_2d(P, phi, t)`, `pj.rot2(phi)`, `pj.dots_image(P, amplitudes, sigmas, n, extent)`,
`pj.project_image(R, state, t, n, extent)`, `pj.random_rotations`, `pj.biased_rotations`,
`pj.rotate_image(img, phi)` (exact, noise-preserving), `pj.shift_image(img, t)`.
Convention: projection along +z of the detector frame; source rotated by R
(`x → R x`); Fourier slice FT2[P_Rρ](k) = ρ̂(Rᵀ(kx,ky,0)) (verified: complex
correlation 0.9999 with Rᵀ vs ≤ 0.41 with R).

### `noise_and_ctf as nc`
`nc.ctf_1d(k)` = −sin(πDk² − (π/2)Ck⁴ + w), D=0.12, C=6e-4, w=asin 0.07 (dimensionless,
illustrative; zeros at k ≈ 2.88, 4.15, 5.16, 6.04 …), `nc.zeros(kmax)`, `nc.ctf_2d`,
`nc.apply_ctf(img)`, `nc.apply_ctf_1d(sig, dx)`, `nc.rng(k)`, `nc.add_noise`,
`nc.snr_amplitude(signal, sigma)` = std(signal)/σ (the amplitude SNR used everywhere).

### `fourier as fo`
`fo.fft2c/ifft2c/fft3c/ifft3c`, `fo.log_magnitude(F)` (display only; the theorem is about
complex values), `fo.central_slice(F3, R)`, `fo.plane_normal(R)`,
`fo.reconstruct_direct_fourier`, `fo.sampling_coverage`, `fo.ball_fraction`,
`fo.fsc(A, B, spacing)` → (s, fsc, counts), `fo.crossing(s, curve, thr)`, `fo.soft_mask`.

### `style as st`
Colours above; `st.T`, `st.M`, `st.title_card`, `st.illustration_label`,
`st.simplified_note`, `st.caption`, `st.takeaway_card(lines, math)`, `st.legend`,
`st.heatmap(arr, height, symmetric=False, stops=st.GRAY|st.DIVERGING)` (nearest-neighbour
pixels; arrays are [x, y]-indexed, x right, y up), `st.framed(img)`, `st.fit_width`.

### `numeric_examples as ne` — data per clip (load with `cached(name, fn)`)

- **`cached("forward_model", ne.forward_model)`** (clip 01): `x` (256 samples on [-2,2)),
  `s` (1D profile of the projected dots), `s_ctf`, `eps`, `I` = s_ctf + eps, `sigma_1d`,
  `k`, `ctf` (CTF on k∈[0,8]), `ctf_zeros`, `k_bands` = [lobe-1 centre (neg), first zero,
  lobe-2 centre (pos)], `ctf_at_bands` (≈ −0.75, 0, +0.99), `bands_before`/`bands_after`
  (3×256 band components before/after the CTF; very different amplitudes — normalise each
  band for display and say so), `img`, `img_ctf`, `noise_frames` (6 frames: CTF image +
  α·noise, α = `noise_alphas`), `sigma_2d`, `ctf2d` (centred 2D CTF), `points3d`, `proj2d`, `R`.
- **`cached("alignment", ne.alignment)`** (clip 02, 64 px, extent 2.4): `S` (clean 4-dot
  image), `sigma`, `snr1`=0.5, `phi` (deg), `t` (n_obs×2), `obs12`, `clean12`, `aligned12`
  (first 12 observations; aligned with the TRUE poses), `ns`=[1,4,16,64,256],
  `unaligned_mean`, `aligned_mean`, `residual` (= aligned_mean − S), `snr_measured`
  (≈ [0.50, 1.03, 2.08, 4.15, 8.18]), `snr_theory` = 0.5√N, `sd_measured`, `sd_theory`,
  `corr_unaligned`, `corr_aligned` (corr with S: 0.44 vs 0.99 at N=256), `pattern`, `extent`.
- **`cached("orientation_fourier", ne.orientation_fourier)`** (clip 03): `X`, `views`
  (4 Euler triples), `Rs`, `proj` (4×5×2 projected dots), `imgs`, `ft_mag`, `ft_re`, `ft_im`,
  `vdirs`, `corr_RT` / `corr_R` (numerical slice check for views 2–4), `cov_uniform` (0.92)
  / `cov_biased` (0.26) ball fractions sampled by 60 planes, `normals_uniform` /
  `normals_biased` (12 plane normals each, for drawing), `max_theta_biased`=15°, `ks`
  = [1,2,3,5,10,20] views, `est_uniform` / `est_biased` (6×5×3 least-squares dot estimates
  from k known views with projected-coordinate noise `meas_noise`), `err_uniform`
  (0.61 → 0.011) / `err_biased` (0.57 → 0.044) mean dot errors.
- **`cached("refinement", ne.refinement)`** (clip 04; 24 observations, 40 px): `imgs`,
  `clean`, `sigma`, `lam`, `R_true`, `t_true`, `X_true`, `hist_X` (13×5×3, true frame),
  `loss`, `loss_data`, `loss_reg`, `step` (0 init, 1 after pose step, 2 after source step;
  6 iterations), `pose_err` (13×24 deg; mean 16 → 1.7), `shift_err`, `track` (index 0),
  `track_R`, `track_t`, `track_pred`, `track_res` (residual images), `track_res_norm`,
  `noise_floor` (= expected data term at the true solution), and the 1D regularization
  example `reg_x`, `reg_f`, `reg_xs`, `reg_y`, `reg_fit_ls` (unregularised, oscillating),
  `reg_fit_rg` (smoothness-penalised), `reg_data_res_ls` < `reg_data_res_rg`,
  `reg_true_err_ls` > `reg_true_err_rg`, `reg_penalty_ls/rg`, `reg_lam_smooth`.
  The displayed loss is exactly Σ‖I_i − P[X]‖² + λΣ|x_j|² (CTF = 1 in this toy).
- **`cached("classification", ne.classification)`** (clip 05; 240 observations, one fixed
  view): `cleanA`, `cleanB`, `imgs12`, `labels` (0=A,1=B), `sigma`, `mixed`, `mean_true_A/B`,
  `scores` (240×2 PC1/PC2), `explained`, `pc_images`, `pA` = p(c=A|I_i), `resp_hist`
  (EM iterations × 240), `em_means`, `em_covs`, `em_pis`, `class_mean` (2×48×48,
  responsibility-weighted), `assign`, `ambiguous` (indices with 0.25<pA<0.75), `accuracy`
  (≈0.95), `frac_A`, `site_A`/`site_B` (2D position of dot 5 in each state), `peaks`
  (rows cleanA, cleanB, mixed, class0, class1; cols value at site A, site B), `points_A`,
  `points_B`, `proj_A`, `proj_B`.
- **`cached("fsc_resolution", ne.fsc_resolution)`** (clip 06; 64³ fields, extent 3): `s`
  (1/units), `fsc_half`, `fsc_region1`, `fsc_region2`, `fsc_model`, `fsc_unmasked`, `counts`,
  `s_half`, `d_half` (≈0.34 units), `d_region1` (≈0.33), `d_region2` (≈0.65), `s_model`,
  `d_model` (map vs ideal point-model field at 0.5; ≈ d_half, as expected for a perfect
  model), `X`, `est_A`/`est_B` (24×5×3 independently jittered position estimates per
  half; region 2 = dots 2,5 more scattered), `jitter`, `region1`=[0,2,3], `region2`=[1,4],
  `noise_sd`, `nyquist`, `sec_A/B/model` (central z sections), `proj_A/B/model`
  (projections), `ft2_A/B` (log-magnitude of projection FTs).
