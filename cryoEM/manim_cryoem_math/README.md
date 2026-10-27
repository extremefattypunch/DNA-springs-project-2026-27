# Mathematics of single-particle cryo-EM: six Manim illustrations

Six silent 1920×1080, 30 fps clips built only from abstract dots, arrows, planes,
plots, small numerical matrices and equations. Every clip carries a
**MATHEMATICAL ILLUSTRATION** tag for its full length and ends on a static takeaway
card. The spec is `manim-gpt-animation-plan.md` at the repo root.

| # | File | Question | Target | Rendered |
|---|---|---|---:|---:|
| 01 | `output/01_forward_model.mp4` | How do projection, the CTF and noise alter an observation? | 30 s | 34.0 s |
| 02 | `output/02_alignment_averaging.mp4` | Why align before averaging, and how fast does the noise fall? | 40 s | 44.2 s |
| 03 | `output/03_orientation_fourier_slice.mp4` | How do unknown 3D orientations and Fourier slices define reconstruction? | 70 s | 76.1 s |
| 04 | `output/04_iterative_refinement.mp4` | What does iterative reconstruction optimize? | 40 s | 44.9 s |
| 05 | `output/05_classification_pca.mp4` | How do heterogeneity and classification interact? | 45 s | 49.5 s |
| 06 | `output/06_fsc_resolution.mp4` | How are resolution and reconstruction reliability assessed? | 45 s | 46.1 s |

The total is 4:55 (294.7 s). Each clip is within ±15 % of its target, because the spec
puts clarity ahead of exact timing. `output/thumbnails/` holds the last frame of each
clip (the takeaway card). `output/all_scenes.mp4` is the optional concatenated master,
and `output/render_report.txt` is the ffprobe report (h264, yuv420p, no audio stream).

## Setup (FASRC Cannon)

The environment lives on holylabs. Renders, the Tex cache and the precomputed
numerical data go to netscratch.

```bash
cd cryoEM/manim_cryoem_math
source env/activate_manim.sh     # micromamba env envs/manim (Manim CE 0.22, Pillow 11) + TinyTeX
                                 # sets MANIM_CRYOEM_SCRATCH / MANIM_CRYOEM_CACHE on netscratch
```

On another machine, install `requirements.txt` into Python 3.12, plus a LaTeX that
provides `standalone`, `amsmath` and `dvisvgm`, plus ffmpeg. Then set
`MANIM_CRYOEM_CACHE` to any writable directory.

## Render

Run from this directory, so that `manim.cfg` (background colour, 1080p30) is read.

```bash
# one clip
manim render -qh --fps 30 scenes/s03_orientation_fourier_slice.py Scene03OrientationFourierSlice
# preview quality
manim render -ql scenes/s01_forward_model.py Scene01ForwardModel
# all six, in order (python main.py --list prints module, class, output name)
python main.py --list | while read -r MOD CLS OUT; do
    manim render -qh --fps 30 "scenes/$MOD.py" "$CLS" < /dev/null
done
```

On the cluster the final renders run as one CPU job on the free `test` partition.
The six clips take about 15–25 minutes. Afterwards, `postprocess.sh` writes the
thumbnails, the master video and the report.

```bash
sbatch cryoEM/manim_cryoem_math/slurm/render_final.sbatch      # -> output/*.mp4
bash   cryoEM/manim_cryoem_math/slurm/postprocess.sh            # thumbnails, all_scenes.mp4, report
```

## Tests

```bash
python -m pytest -m "not render" -q   # numerical truth + scope checks (spec §5 items 1-9), ~15 s
python -m pytest -m render -q         # smoke render: one still per clip at low quality
```

- **`tests/test_math.py`** covers spec §5 items 1–8:
  - SO(3) properties, and that a tilt is not an in-plane spin.
  - Projected dots against direct coordinates.
  - The analytic image against a numerical line integral.
  - The Fourier slice check uses **Rᵀ**: correlation 0.9999 against ≤ 0.41 with R.
  - Phases carry the shift.
  - Direct Fourier inversion with known poses, and the coverage loss under biased views.
  - CTF zeros and sign changes, and that a contrast-reversed band is not a blur.
  - σ/√N scaling, and the signal loss from misalignment.
  - Displayed refinement losses are recomputed from the stored predictions.
  - The mixed mean is ambiguous while the class means are not.
  - FSC is 1 for identical maps, about 0 for independent noise, and empty shells are safe.
- **`tests/test_visual_constraints.py`** covers spec item 9:
  - It scans every scene and helper for biological, sample or dataset vocabulary and for Å.
  - It checks for any Manim surface, solid or SVG/image-file class.
  - It checks that every clip uses the illustration tag, a title card and a takeaway card.

## Conventions

- **Source.** The source is five dots x₁…x₅ in arbitrary dimensionless units (`helpers/point_constellation.py`). State B moves only dot 5. Clip 02 uses the x,y coordinates of dots 1–4 as its 2D pattern. All lengths and frequencies are dimensionless, and nothing is in Å.
- **Projection.** The source is rotated, x → Rx, and integrated along z: p_j = (R x_j)_{x,y} + t. The viewing direction is v = Rᵀe_z. Euler angles are ZYZ, R = R_z(ψ)R_y(θ)R_z(α). The first angle is called α so it does not clash with the 2D in-plane angle φ in clip 02.
- **Fourier slice.** FT₂[P_R ρ](k) = ρ̂(Rᵀ(k_x, k_y, 0)). Spectra are complex. Magnitude panels are display-only, and each is normalised to its own maximum.
- **Signal-to-noise.** SNR is the amplitude ratio std(signal)/σ.
- **Images.** Images are tiny Gaussian-dot matrices, from exact line integrals of isotropic Gaussian kernels. They are shown with nearest-neighbour pixels.
- **Colours** are reserved across clips: cyan observed, green truth, orange estimate/arrows, purple alternate state, red noise/residual, grey axes. Yellow is used only for the tag and the "Takeaway" heading. Clip 06 draws half-map B in orange; that clip notes the exception on screen.

## Deliberate simplifications

- **CTF.** The CTF is illustrative: CTF(k) = −sin(πDk² − (π/2)Ck⁴ + asin 0.07), with D = 0.12 and C = 6×10⁻⁴. It is dimensionless and has no envelope or astigmatism. Its zeros fall at k ≈ 2.88, 4.15, 5.16 and 6.04.
- **Noise.** Noise is i.i.d. Gaussian, seeded (`noise_and_ctf.rng`, base seed 2027). Every run renders identical frames.
- **Clip 02.** The model is a 2D toy with known poses, no CTF and one fixed view. Alignment is shown as the true inverse transform, and the noise fall-off follows σ/√N.
- **Clip 03.**
  - The reconstruction uses direct Fourier insertion and a least-squares inversion over nested view sets. Coverage is the fraction of the frequency ball that the slices sample.
  - The biased set keeps v within 15° of k₃.
  - The numerical check uses views 2–4. View 1 is the identity, where R = Rᵀ, so it cannot tell the two apart.
- **Clip 04.**
  - Refinement alternates L-BFGS over poses and dot positions. It uses 24 observations, CTF = 1, λ = 2 and a quadratic penalty.
  - Pose errors are measured after Kabsch superposition, which removes the global rotation gauge.
  - The loss is a toy objective, not any software package's loss. One observation stays in a local minimum (about 22°), and the clip shows it.
- **Clip 05.**
  - A two-state mixture is fitted at amplitude SNR 0.15, with 240 observations of one view.
  - PCA is done on the image vectors, and a 2-component Gaussian mixture is fitted by EM on the first two scores.
  - The class means are posterior-weighted soft means.
  - PC2 is essentially a noise direction, and the clip says so.
- **Clip 06.**
  - Each half-map is the average of jittered-position fields plus independent noise.
  - The two regions differ only in positional jitter, 0.03 against 0.18. Local resolution is shown as a masked FSC per region.
  - Map against model uses the ideal point model with the illustrative 0.5 threshold. Half-maps use 0.143.
  - In this toy, d_model ≈ d_half is a coincidence of the setup, not a general rule.

## Layout

```text
main.py                 scene registry (module, class, output name); `--list`
manim.cfg               1920x1080, 30 fps, background #0E1116, media_dir on netscratch
scenes/s0N_*.py         one Scene class per clip
helpers/                point_constellation, projections, noise_and_ctf, fourier,
                        numeric_examples (cached per-clip data), style, cache
tests/                  test_math, test_visual_constraints, test_smoke_render
slurm/                  render_final.sbatch, postprocess.sh
AUTHORING.md            authoring guide: palette, layout rules, lessons learned, helper API
output/                 final MP4s, thumbnails/, all_scenes.mp4, render_report.txt
```
