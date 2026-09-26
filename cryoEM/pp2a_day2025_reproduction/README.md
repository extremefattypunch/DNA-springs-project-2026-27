# Reproducing Day et al. 2025 (PP2A Aα vs Aβ holoenzymes) — cryo-EM and MD

Paper: Day, Huang, Leonard, O'Connor, Narla & Taylor, *Structure* 33:1688–1699 (2025), doi 10.1016/j.str.2025.06.013.
Plan: [`../gpt-plans/pp2a-day2025-reproduction.md`](../gpt-plans/pp2a-day2025-reproduction.md).
**`report.md` is the findings document** (answers questions A–H of the plan). `results/*.csv` hold every
number quoted there. `DATA.md` says where the bulk data is and when netscratch will purge it. This file is
orientation and the runbook only.

## Short answer

The deposited maps and models reproduce the paper's map/model-level claims: Table S1 validation to the last
digit, the concave-Aβ geometry (HEAT1–15 62.7→63.3 Å vs 59.6→60.4 Å; salt bridge 13.2→12.7 Å vs 12.4→12.4 Å),
and the Fig 2B contrast (B56γ C-terminal helix density at 6 Å: ratio 0.60 in Aα vs 0.30 in Aβ). Nothing at the
particle level (2D classes, viewing directions, Fig 4/S7 3D classification, Fig S6 counts) can be reproduced
because no movies or particle stacks were deposited. The MD (Fig 1, Videos S1/S2) runs in OpenMM on free
`gpu_test` MIG slices and is reported as it finishes — see `report.md` §A.

## What is where

```
scripts/      00-04 inputs, models, tleap systems, OpenMM engine, MD submission
              10-15 MD analysis and rendering (RMSD, RMSF, Fig 1C, Videos S1/S2)
              20-28 cryo-EM: CryoSPARC local-res/FSC, RELION post-processing, EMReady, Phenix, geometry,
                    ChimeraX renders, matplotlib plots, montages
              40    software ledger  ->  software_versions.json
slurm/        build/install jobs (RELION+GROMACS, ChimeraX, EMReady), CPU wrapper, MD chain
env/          activate_pp2a.sh  - the only environment file; sets paths for every tool
figures/      figN_*.png (paper figure IDs), videoS1_aa.mp4 / videoS2_ab.mp4
results/      tableS1_comparison.csv, em_metrics.csv, geometry_metrics.csv, fig2b_metrics.csv, md_*.csv
```

## Running it

```bash
source env/activate_pp2a.sh                      # paths for dnaspring python, Phenix, RELION, GROMACS, ChimeraX, EMReady
# software (once; each job is idempotent)
sbatch slurm/pp2a_build_software.sbatch          # RELION 3.1 + GROMACS 2023.3 -> holylabs software/
sbatch slurm/pp2a_install_chimerax.sbatch        # ChimeraX 1.12 in a Stream 9 container + rpm sysroot (OSMesa offscreen)
sbatch slurm/pp2a_install_emready.sbatch         # EMReady v1.1 venv
# inputs (EMDB maps, PDB models, AlphaFold DB models, UniProt)
$DNASPRING_PY scripts/00_fetch_inputs.py
# cryo-EM, deposited data only (CryoSPARC needs the master up: `launchcryosparc` on the VDI, GUI open)
PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/21_cryosparc_localres.py       # Import -> Local Resolution -> Local Filter -> Validation (FSC), lane Cannon-test
for i in Aa Ab; do sbatch -J pp2a-relion-$i -c 4 --mem=16G slurm/pp2a_cpu.sbatch bash scripts/22_relion_postprocess.sh $i; done   # ~55 min each (locres)
for i in Aa Ab; do sbatch -J pp2a-masksweep-$i -c 2 slurm/pp2a_cpu.sbatch bash scripts/22b_relion_mask_sweep.sh $i; done
for i in Aa Ab; do sbatch -J pp2a-phenix-$i -c 2 slurm/pp2a_cpu.sbatch bash scripts/24_phenix_validate.sh $i; done              # mtriage + molprobity, ~3 min
for i in Aa Ab; do sbatch -J pp2a-emready-$i -p gpu_test --gres=gpu:1 -c 4 --mem=24G slurm/pp2a_cpu.sbatch bash scripts/23_emready.sh $i; done
$DNASPRING_PY scripts/25_geometry.py            # Fig 3A-E numbers        -> results/geometry_metrics.csv
$DNASPRING_PY scripts/25b_fig2b_density.py      # Fig 2B density ratios   -> results/fig2b_metrics.csv
$DNASPRING_PY scripts/26_chimerax_figures.py    # all ChimeraX panels (offscreen; ~3-6 min per panel, run one at a time)
$DNASPRING_PY scripts/27_em_figures.py          # FSC plot, local-res histograms, Fig 3B, Table S1 comparison
$DNASPRING_PY scripts/28_montage.py             # paper-style multi-panel figures
# MD (Tier 1: one 100 ns run per isoform; ~2-3.5 days each on one MIG slice, run back-to-back because
# gpu_test allows two jobs per user and the OOD desktop holds one)
$DNASPRING_PY scripts/01_prepare_models.py && $DNASPRING_PY scripts/02_build_systems.py
bash scripts/04_submit_md.sh                    # 1 ns benchmark, then the self-resubmitting 100 ns chain, per system
$DNASPRING_PY scripts/10_md_rmsd_rmsf.py && $DNASPRING_PY scripts/13_md_figures.py   # when progress.json says 100 ns
$DNASPRING_PY scripts/14_md_fig1c_chimerax.py && $DNASPRING_PY scripts/15_md_movies.py
$DNASPRING_PY scripts/40_software_versions.py
```

## Gotchas met on the way (so nobody meets them twice)

* **CryoSPARC v5 half-map inputs are one group with two slots.** Import Volumes makes one slot per job, so the
  second half map (and the local-res map for Local Filtering) must be attached with
  `cs.api.jobs.connect_result(P, J, "volume", 0, "map_half_B", ...)`; a second `connect` fails with
  "Cannot connect ... more than 1 time(s)". Projects are created with `cs.api.projects.create(title=, parent_dir=)`.
* **RELION treats `.map` as SPIDER** — symlink the EMDB maps as `.mrc`. `--o postprocess` must be `--o ./postprocess`
  or RELION tries to write `postprocess/logfile.pdf`. RELION `--locres` prefixes outputs with the rootname
  (`locres_locres.mrc`). The gold-standard resolution depends on the mask: contour-matched tight mask 3.93/3.55 Å,
  loose 15 Å-lowpass mask 3.54/3.33 Å (plateau), deposited CryoSPARC curve 3.42/3.22 Å.
* **ChimeraX 1.12 has no RHEL 8 build.** The el9 rpm runs on Rocky 8 inside a CentOS Stream 9 Singularity image
  (Rocky 9 lacks the OpenSSL 3.3 it links), with X11/GL/OSMesa rpms unpacked into a sysroot — no root needed.
  Two offscreen renders in parallel hung; run panels sequentially. `distance` creates a pseudobond model, so
  open files with explicit `id N`. `supersample 2` is enough and half the time of 3.
* **gpu_test QOS = 2 jobs per user**, and the Open OnDemand desktop is one of them. `test` allows 5.
* **EMReady v1.1 needs `np.asfarray` → `np.asarray`** with NumPy 2, and its `EMReady.sh` hard-codes the author's
  conda paths (patched by the install job).
* **AlphaFold DB is now v6**; ask the API for the current file name rather than hard-coding `_v4`.
* The Aβ AF model's 25-residue tail sticks out, so the 16 Å-padded box has 226k particles vs 163k for Aα.
