# Plan: reproduce the key cryo-EM and MD results of Day et al. 2025 (PP2A Aα vs Aβ)

Paper: Day A, Huang W, Leonard D, O'Connor CM, Narla G, Taylor DJ. *Structure* 33:1688–1699 (2 Oct 2025), doi 10.1016/j.str.2025.06.013. Source read in full: `reading list/pp2a MD and cryo/mmc4.pdf` (main text, STAR methods, Figs S1–S10, Table S1).

## Context

Ian wants the paper's **key** computational results (figures, plots, videos) reproduced on the Hekstra lab's FASRC Cannon allocation with his CryoSPARC v5 instance, following the `hekstra-cluster` and `cryosparc` skills. Decisions made with Ian (2026-09-19, 2026-09-26):

1. **No raw cryo-EM data is public** (EMDB has only primary map, two half-maps and an FSC XML per entry; no masks; **no EMPIAR entry**; no MD trajectories or code). Ian: **analyse the deposited maps and models directly; no synthetic-particle test of the CryoSPARC pipeline.** Particle-level results (2D classes, viewing distributions, Fig 4/S7 3D classification, Fig S6 counts) are documented as not reproducible.
2. **Software policy (Ian): install missing open-source tools on holylabs; do not depend on SBGrid.** Phenix 2.2 (release 6143) comes from the labmate's install (`source /n/lab_storage/hekstra_lab/people/aldama/software/phenix2/phenix_env.sh`, tested by Ian; the script only sets `PHENIX*` and `PATH`). Coot is launched through Phenix (https://phenix-online.org/documentation/coot.html) for interactive inspection. Installed by this plan on `/n/holylabs/hekstra_lab/Everyone/ianpoon/` (~1.9 TiB free as of 2026-09-26): **RELION 3.1**, **GROMACS 2023.3**, **ChimeraX**, **EMReady**. **NAMD 3.0.1 and VMD 1.9.4** are free for academic use but their downloads require a registration login; if Ian downloads the tarballs into `inputs/software/`, they are installed and used (paper-faithful); otherwise MD runs in the already-installed **OpenMM 8.6** (same ff19SB/OPC Hamiltonian) and VMD's roles go to ChimeraX/MDAnalysis.
3. **MD Tier 1 only** (one 100 ns run per isoform = Fig 1A–C + Videos S1/S2; Fig S3 replicates deferred). Starting models from **AlphaFold Server (AF3) supplied by Ian**, AlphaFold DB (AF2) fallback.
4. **Compute**: everything on the free `gpu_test`/`test` partitions (CryoSPARC lane `Cannon-test`); the paid `gpu` partition and lane `Cannon` are never used.

## 1. Key results and reproduction route

| ID | Paper item | Route | Deliverable |
|---|---|---|---|
| K1 | Fig 1A per-residue RMSF Aα vs Aβ | MD → `gmx rmsf` (GROMACS 2023.3, as the paper) | `fig1a_rmsf.png` |
| K2 | Fig 1B RMSD vs time: N-terminus (1–13 / 1–25) and excluding N-terminus (14–589 / 26–601) | VMD RMSD Trajectory Tool if VMD supplied, else MDAnalysis with the same convention | `fig1b_rmsd_nterm.png`, `fig1b_rmsd_core.png` |
| K3 | Fig 1C N-terminus overlay every 20 ns, red (0 ns) → blue (100 ns) | ChimeraX offscreen | `fig1c_nterm_overlay_aa.png`, `_ab.png` |
| K4 | Videos S1/S2 smoothed trajectory movies, N-terminal residues blue | VMD (smoothing 10) if supplied, else ChimeraX `coordset` on a 10-frame running-mean trajectory; ffmpeg | `videoS1_aa.mp4`, `videoS2_ab.mp4` |
| K7 | Fig 2A EMReady-sharpened maps with models; Fig 2B 6 Å lowpass showing the B56γ C-terminal helix in Aα only | EMReady + RELION image handler + ChimeraX + python ratio | `fig2a_maps.png`, `fig2b_lowpass6_b56_helix.png` |
| K8 | Fig 3A–C: HEAT1–15 CoM distance 62.7 Å (Aα) vs 59.6 Å (Aβ); per-HEAT CoM shift (max HEAT 3–6, ~2–3 Å); salt-bridge Cα distance 13.2 Å (R113–E191) vs 12.4 Å (R125–E203) after alignment on HEAT 12–15 | python geometry + ChimeraX | `fig3a_heat_com.png`, `fig3b_com_shift.png`, `fig3c_salt_bridge.png` |
| K9 | Fig 3D,E interface residues conserved | python contacts + ChimeraX | `fig3de_interfaces.png`, contact table |
| K11 | Fig S5A/E local-resolution maps; S5C/G half-map FSC and map-model FSC | RELION localres/postprocess + CryoSPARC Local Resolution + Phenix mtriage | `figS5_localres.png`, `figS5_fsc.png` |
| K13 | Table S1 map resolution, model resolution (FSC 0.5), MolProbity, clashscore, Ramachandran, rotamers, bond RMSDs | RELION postprocess, Phenix mtriage + molprobity | rows in `metrics.csv` |

Not reproducible (report §H): Fig S5B/F 2D classes, S5D/H viewing distributions, Fig 4/S7 3D classification and C-tail occupancy counts, Fig S6 particle counts, TOPAZ/blob picking. Out of scope: Figs 5–7, S4, S8–S10 (wet lab); Fig S3 (Tier 2 MD); model rebuilding (deposited 9MIP/9MF5 reused).

## 2. Paper parameters the implementation follows

**MD**: AF3 full-length Aα (589 aa, P30153) and Aβ (601 aa, P30154; N-terminal extension MAGASELGTGPG). ff19SB + OPC via tleap; rectangular box, 16 Å padding (measured 18–19.5 Å); Na⁺ neutralise (23 Aα / 25 Aβ) + NaCl to 150 mM. NAMD 3.0.1: minimise 10,000 steps with heavy atoms restrained; heat 0→300 K, +10 K per 10,000 steps with temperature reassignment; restraints scaled 1.0→0 in 0.2 decrements every 50,000 steps; production 100 ns NPT, Langevin 300 K damping 5 ps⁻¹, Langevin piston 1 atm period 100 fs decay 50 fs. RMSD vs frame 1 (VMD RMSD Trajectory Tool), RMSF (GROMACS 2023.3), movies VMD smoothing window 10; Fig 1C overlay at 20 ns increments.

**Cryo-EM**: 1.056 Å/px (Aα) / 1.058 Å/px (Aβ), box 256; 3.4 / 3.2 Å (FSC 0.143); local-res range 5.87–3.04 / 5.22–2.82 Å; model resolution 3.8 / 3.6 Å (FSC 0.5); MolProbity 1.25 / 1.15, clashscore 4.88 / 3.59, Ramachandran favoured 98.08 / 98.22 %, bond RMSD 0.010 / 0.009 Å, angle 1.291 / 1.241°. Sharpening EMReady; local resolution RELION 3.1; Fig 2B = 6 Å lowpass of the raw map. Models: 9MIP (A 7–589, B56γ 21–449, Cα 2–295+304–309), 9MF5 (A 25–601, B56γ 38–403, Cα 3–295+302–309); start model 3FGA. Recommended contours (EMDB): 0.0777 (Aα), 0.154 (Aβ).

## 3. Layout and conventions (as in the repo's three prior cryo-EM tests)

Repo (code, figures, small results; home is snapshotted):
```
cryoEM/gpt-plans/pp2a-day2025-reproduction.md      # this plan, copied in at step 0
cryoEM/pp2a_day2025_reproduction/
  README.md  report.md  DATA.md  metrics.csv  software_versions.json
  inputs/af3/                 # Ian drops AlphaFold Server CIFs here (P30153, P30154); AFDB fallback otherwise
  inputs/software/            # optional: NAMD_3.0.1_Linux-x86_64-multicore-CUDA.tar.gz, vmd-1.9.4*.tar.gz (Ian downloads)
  env/activate_pp2a.sh        # DNASPRING_PY, MPLBACKEND=Agg, Phenix env, PATH for holylabs installs (see §5A)
  namd/                       # NAMD .conf templates (used only if NAMD is supplied)
  slurm/pp2a_md_chain.sbatch  pp2a_build_software.sbatch  pp2a_cpu.sbatch  pp2a_emready.sbatch
  scripts/                    # numbered, §6
  figures/
```
Data (bulk, regenerable, **90-day purge**): `/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025/{inputs,md/{systems,runs},em/{emdb,models,relion,phenix,emready,chimerax},cs_state}`; CryoSPARC project `/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-pp2a-day2025`; logs `/n/netscratch/hekstra_lab/Lab/ian_poon/logs/`. Software on holylabs under `/n/holylabs/hekstra_lab/Everyone/ianpoon/software/` (§5A). Nothing written to `/n/lab_storage` (Phenix is only read from there).

Reuse: `cryoEM/minimal_force_cryoem_test/scripts/{cs_connect.py,mrcio.py}`; `run_job` pattern (create → queue on `Cannon-test` → poll, `cs_state.json` resume) from `cryoEM/minimal_mechanics_test_3/scripts/04_run_compressed_3dva.py`; `MD simulations/sfGFP/runs/openmm_engine.py` (OpenMM fallback engine: already has `k_res` restraint ramp, Monte Carlo barostat, checkpoint resume, `--no-hmr`); `MD simulations/sfGFP/slurm/md_chain.sbatch` (gpu_test conventions: `-c 1`, `--mem=8G`, 12 h, restart-safe); `MD simulations/sfGFP/analysis/vizstyle.py` (plot styling). Never source `MD simulations/sfGFP/env/activate.sh` (points at full lab_storage). Envs are activated by full prefix (`micromamba activate -p …/envs/dnaspring`).

## 4. MD track (Tier 1)

**Inputs** (`00_fetch_inputs.py`): use `inputs/af3/*.cif` if present, else download `https://alphafold.ebi.ac.uk/files/AF-P30153-F1-model_v4.pdb` and `AF-P30154-F1-model_v4.pdb`; verify sequences (589/601 aa; Aβ begins MAGASELGTGPG); record which model was used.

**Build** (`02_build_systems.py`, tleap in dnaspring): strip H and hetero atoms; rotate principal axes onto x/y/z (mdtraj) so `solvateBox … OPCBOX 16.0` pads the true extent; `source leaprc.protein.ff19SB; source leaprc.water.opc`; `addIons Na+ 0` then `addIonsRand Na+ N Cl- N`, N = round(0.15 × N_water / 55.5); record padding per axis; expect ~165k particles (OPC has 4 sites), assert < 200k; write prmtop/inpcrd/PDB and `residue_map.json` (assert tleap numbering equals UniProt numbering).

**Engine A: NAMD 3.0.1 CUDA (if Ian supplies the tarball)** (`03_write_namd_configs.py` fills `namd/*.conf`): `amber yes`, `parmfile`, `ambercoor`, `readexclusions yes`, `exclude scaled1-4`, `1-4scaling 0.833333`, `cutoff 10`, `switching off`, `pairlistdist 12`, `PME yes`, `PMEGridSpacing 1.0`, `rigidBonds all`, `timestep 2.0`, `fullElectFrequency 2`, `wrapAll on`. Stages: `min.conf` (`constraints on`, heavy-atom restraint PDB k = 10 kcal/mol/Å², `minimize 10000`); `heat.conf` (`temperature 0`, `reassignFreq 10000`, `reassignIncr 10`, `reassignHold 300`, `run 300000`); `release.conf` (Langevin 300 K damping 5; `langevinPiston on`, target 1.01325 bar, period 100, decay 50; `constraintScaling` 0.8/0.6/0.4/0.2/0.0 each 50,000 steps); `prod.conf` (100 ns, `dcdfreq 50000` → 1,000 frames, `restartfreq 500000`). Restart-safe: the sbatch detects the latest `.restart.*`, sets `bincoordinates/binvelocities/extendedSystem/firsttimestep`, runs the remaining steps; DCD segments concatenated at analysis. Run as `namd3 +p1 +devices 0 <conf>`.

**Engine B: OpenMM 8.6 (default if NAMD is not supplied when systems are ready)** (`03_openmm_pp2a_engine.py`, fork of `openmm_engine.py`; deviations recorded in a README table): 2 fs, `--no-hmr`; minimise 10,000 steps **with** `k_res` 10 kcal/mol/Å² on protein heavy atoms, then reset the restraint reference to minimised coordinates; heat NVT 30 × 10,000 steps (`setTemperature(10·s K)`); release NPT with `MonteCarloBarostat(1 bar, 300 K, 25)` as the piston substitute, `k_res` × 0.8/0.6/0.4/0.2/0.0 each 50,000 steps; production `LangevinMiddleIntegrator` 300 K friction 5 ps⁻¹, PME 10 Å, 100 ns, DCD/state every 100 ps, checkpoint every 1 ns; seeds 2001/2002.

**Slurm** (`slurm/pp2a_md_chain.sbatch`, `04_submit_md.sh`): `#!/bin/bash -l`, `-p gpu_test --gres=gpu:1 -c 1 --mem=8G -t 12:00:00`, logs on netscratch, one run per MIG slice (no packing). **First a 1 ns benchmark per system**; from measured ns/day the chain self-submits with `--dependency=afterany` and `--time` = min(12 h, 1.15 × remaining/speed) so the last link does not trip the Job Defense Shield. Expected 30–70 ns/day → 100 ns in 1.5–3.5 days per isoform, ~50–80 MIG-hours each, zero billing. `seff`/`jobstats` per job into `DATA.md`; `cluster-health --quick` before starting. Afterwards a protein-only trajectory (`.xtc` via MDAnalysis, ~150 MB) is written for analysis; full DCDs stay on netscratch.

**Analysis**:
- RMSD (`10_md_rmsd.py` or `10_md_rmsd_vmd.tcl` with `vmd -dispdev text` if VMD supplied): selections N-term `resid 1:13` / `1:25`, core `14:589` / `26:601`, backbone; RMSD vs frame 0 with each selection's own fit (RMSD Trajectory Tool convention; the paper's Fig 1B and its 10–15 vs 20–25 Å whole-protein numbers); core-fit RMSD of the N-terminus as a supplementary panel.
- RMSF (`11_md_rmsf_gromacs.sh`): `gmx rmsf -s ref.pdb -f protein.xtc -res -fit` → per-residue Cα RMSF (full 0–100 ns as in the paper; 10–100 ns also reported).
- `12_md_collect.py` → `metrics.csv`; `13_md_figures.py` (imports `vizstyle`) → Fig 1A/1B analogues with the paper's axis ranges (RMSF 0–25 Å, RMSD 0–30 Å, 0–100 ns).
- Fig 1C (`14_md_fig1c_chimerax.py`, `chimerax --nogui --offscreen`): frames 0/20/40/60/80/100 ns aligned on core Cα (`matchmaker`/`align`), core cartoon from frame 0 (grey Aβ, red Aα), N-terminus per frame coloured red→blue, colour key 0–50–100 ns, `save … width 2400 height 1800 supersample 3`.
- Videos (`15_md_movies.py`): 10-frame running mean of coordinates (VMD smoothing 10 equivalent) written as a 600-frame trajectory; VMD `render Tachyon` per frame if supplied, else ChimeraX `movie record … coordset … movie encode` offscreen; N-terminal residues blue, time stamp; ffmpeg 30 fps ≈ 20 s, `-crf 20 -pix_fmt yuv420p`. Runs on `test` (CPU).

## 5. Cryo-EM track (deposited maps and models)

**Fetch** (`00_fetch_inputs.py`, aria2c): EMD-48300/48224 primary maps, half-maps and FSC XML from `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-4xxxx/{map,other,fsc}/` (~380 MB); 9MIP, 9MF5, 3FGA mmCIF/PDB from `https://files.rcsb.org/download/`. md5 and sizes into `DATA.md`.

**Raw-map proxy**: the unsharpened map is not deposited; `raw = (half1 + half2)/2` (`relion_image_handler --add`, cross-checked in numpy via `mrcio.py`) stands in for it. Recorded as a deviation.

**Step A. Software installs on holylabs** (`slurm/pp2a_build_software.sbatch` on `test`, 8 cores, 4 h, `#!/bin/bash -l`, `module load gcc/12.2.0-fasrc01 cmake/3.28.3-fasrc01`; builds in `/scratch/ianpoon_$SLURM_JOB_ID`; every source tree and tarball deleted after install; `cluster-health --quick` first):
- **RELION 3.1** (GPL): `git clone -b ver3.1 https://github.com/3dem/relion.git`; `cmake -DCMAKE_INSTALL_PREFIX=…/software/relion-3.1 -DGUI=OFF -DCUDA=OFF -DFORCE_OWN_FFTW=ON ..`; `make -j8 install` (~0.5 GB). Fallback: gcc 10.2 module, else `ver4.0`.
- **GROMACS 2023.3** (LGPL): tarball from ftp.gromacs.org; `cmake -DCMAKE_INSTALL_PREFIX=…/software/gromacs-2023.3 -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF -DGMX_MPI=OFF -DREGRESSIONTEST_DOWNLOAD=OFF ..`; `make -j8 install` (~0.3 GB). Only `gmx rmsf`/`gmx trjconv` are needed.
- **ChimeraX** (free for academic use, direct download): Linux generic tarball from https://www.cgl.ucsf.edu/chimerax/download.html → `…/software/chimerax-<ver>/`; verify `chimerax --nogui --offscreen --cmd "open 9mip; save /tmp/t.png; exit"` (offscreen rendering uses the bundled OSMesa). ISOLDE is not needed.
- **EMReady** (GPL, GitHub huanglab): `python -m venv …/software/emready-venv` (python from dnaspring), `pip install torch --index-url https://download.pytorch.org/whl/cu121` + EMReady requirements, clone repo and model weights (~5 GB total). Run on `gpu_test`.
- **NAMD 3.0.1 / VMD 1.9.4 (optional)**: if `inputs/software/NAMD_*.tar.gz` / `vmd-*.tar.gz` exist (Ian's registration downloads), extract to `…/software/namd-3.0.1` and install VMD to `…/software/vmd-1.9.4` (`./configure LINUXAMD64 OPENGL … ; cd src && make install` with text-mode use only). Otherwise skip and log the fallback.
- `env/activate_pp2a.sh` exports PATH entries for each installed tool; `40_software_versions.py` records `relion_postprocess --version`, `gmx --version`, `chimerax --version`, EMReady commit, `namd3`/`vmd` if present, `phenix.version`, CryoSPARC version, dnaspring package versions.

**Step B. CryoSPARC** (Ian starts it with `launchcryosparc` on the VDI so the GUI is open; `export CS_MASTER_TIME=1-00:00` first). `20_cs_connect_and_register.py` asserts the v5 job-type strings from `cs.job_register`; `21_cryosparc_localres.py` creates project `CS-pp2a-day2025` under the required parent and, per isoform: Import Volumes ×3 (primary, half A, half B) → **Local Resolution Estimation** (half-maps, FSC 0.5) → **Local Filtering** → **Validation (FSC)**. Lane `Cannon-test`; < 1 MIG-hour total.

**Step C. RELION post-processing** (`22_relion_postprocess.sh`, `pp2a_cpu.sbatch` on `test`): per isoform, mask `relion_mask_create --i raw.mrc --ini_threshold <contour scaled to the unsharpened map> --extend_inimask 3 --width_soft_edge 6 --angpix 1.056`; gold-standard FSC `relion_postprocess --i half1 --i2 half2 --mask mask.mrc --angpix 1.056 --auto_bfac --o postprocess` → resolution at 0.143 (target 3.4 / 3.2 Å); local resolution `relion_postprocess … --locres --locres_sampling 15 --o locres` → `relion_locres.mrc`, `locres_filtered.mrc` (Fig S5A/E; target range 5.87–3.04 / 5.22–2.82 Å); Fig 2B input `relion_image_handler --i raw.mrc --lowpass 6 --o raw_lp6.mrc`.

**Step D. EMReady** (`23_emready.sh`, `pp2a_emready.sbatch` on `gpu_test`): `EMReady.sh raw.mrc emready.mrc` per isoform (default settings, output resampled to the original grid); compare with the deposited primary map (real-space CC, FSC) to test whether the deposited map is the EMReady output; used for Fig 2A beside the deposited map.

**Step E. Phenix** (`24_phenix_validate.sh`, `test` partition, `source …/phenix_env.sh`): `phenix.mtriage map.mrc half1 half2 model.cif resolution=3.4` → half-map FSC, map-model FSC at 0.5 (Table S1 model resolution 3.8 / 3.6 Å), d99; `phenix.molprobity model.cif` → MolProbity score, clashscore, Ramachandran, rotamers, bond/angle RMSDs. Optional single-macrocycle `phenix.real_space_refine` sanity run only if the statistics disagree with Table S1 beyond rounding. Coot via Phenix in the VDI for interactive inspection; no rebuilding planned.

**Step F. Geometry** (`25_geometry.py`, dnaspring python with mdtraj/gemmi):
- HEAT boundaries from DSSP on 9MIP chain A (helices ≥ 7 residues, expect 30, paired consecutively), cross-checked against the nominal grid HEAT k = 8+39(k−1) … 8+39k−1 (anchors R113 ∈ HEAT3, E191 ∈ HEAT5, W257 ∈ HEAT7, K416/R418 ∈ HEAT11, D572/D574 ∈ HEAT15); Aβ = Aα + 12. Both definitions reported; must agree < 0.5 Å.
- HEAT1–HEAT15 Cα-CoM distance (intramolecular; targets 62.7 / 59.6 Å); per-repeat CoM shift after Kabsch on HEAT 12–15 Cα (target max at HEAT 3–6, 2–3 Å); salt bridge Cα–Cα R113–E191 vs R125–E203 (13.2 / 12.4 Å) plus minimum side-chain N–O distance; interface residues of chain A within 4 Å of B56γ and of Cα in each model, checked against the paper's 11 pairs.
- **Fig 2B quantification**: candidate C-helix = 9MIP chain B residues 404–449 (9MIP models B56γ to 449, 9MF5 stops at 403; verify by inspection that it is the helix reaching the A N-terminus, else take the helix from 3FGA/7K36 superposed on the B56γ core). Ratio r = ⟨ρ⟩(helix mask, 3 Å dilation) / ⟨ρ⟩(B56γ core 38–403) in `raw_lp6.mrc`; for Aβ, transfer the helix atoms by Kabsch fit of 9MIP B 38–403 onto 9MF5 B (assert RMSD < 1.5 Å). Pre-registered expectation: r_Aα ≥ 0.6, r_Aβ ≤ 0.25.

**Step G. Figures** (`26_chimerax_figures.py` run via `chimerax --nogui --offscreen --script`; `27_em_figures.py` with matplotlib + `vizstyle`): Fig 2A (deposited and EMReady-resharpened maps at recommended contour with models; Aα red / Aβ grey / B56γ blue / Cα orange as in the paper); Fig 2B (6 Å lowpass raw maps at matched level, helix region indicated); Fig 3A (models superposed on HEAT 12–15 with HEAT labels and CoM distance bars); Fig 3B bar chart; Fig 3C salt-bridge residues with mesh and distance labels; Fig 3D/E interface residues with mesh; Fig S5A/E map coloured by RELION local resolution (`color sample`, palette range as the paper); Fig S5C/G FSC curves (RELION, CryoSPARC, mtriage, deposited XML on one axis).

## 6. Scripts (execution order) and deliverables

`00_fetch_inputs.py` · `01_prepare_models.py` · `02_build_systems.py` · `03_write_namd_configs.py` / `03_openmm_pp2a_engine.py` · `04_submit_md.sh` · `10_md_rmsd.py` (or `.tcl`) · `11_md_rmsf_gromacs.sh` · `12_md_collect.py` · `13_md_figures.py` · `14_md_fig1c_chimerax.py` · `15_md_movies.py` · `20_cs_connect_and_register.py` · `21_cryosparc_localres.py` · `22_relion_postprocess.sh` · `23_emready.sh` · `24_phenix_validate.sh` · `25_geometry.py` · `26_chimerax_figures.py` · `27_em_figures.py` · `40_software_versions.py` · `41_readme_and_mirror.py` (generates README/DATA, mirrors figures, CSVs, scripts, configs only). Shared: `cs_connect.py`, `mrcio.py`, `vizstyle.py`.

`metrics.csv` schema: `track,system,quantity,selection_or_method,value,unit,paper_value,notes` (one row per number quoted in `report.md`).

`report.md` questions: **A** Do N-terminal RMSD/RMSF separate Aβ from Aα as in Fig 1, with Aβ whole-protein RMSD reaching 20–25 Å? **B** Which engine and starting models were used, and how does a single replicate limit A? **C** Are the deposited map resolutions (RELION, CryoSPARC, mtriage) and Table S1 model statistics reproduced? **D** Is the B56γ C-helix density present at 6 Å in Aα and absent in Aβ, quantitatively? **E** Are the HEAT CoM distances, per-repeat shifts and salt-bridge distances reproduced within 0.5 Å? **F** Are the 11 interface residues in contact in both models? **G** Does EMReady on the half-map average reproduce the deposited sharpened map? **H** What could not be reproduced and why (no EMPIAR data: 2D classes, viewing distributions, 3D classification and C-tail counts, Fig S6; any engine/tool substitutions; Tier 2 MD deferred).

Pre-registered success criteria: A: Aβ N-term RMSF mean > 1.5 × Aα and Aβ whole-protein RMSD plateau exceeds Aα by ≥ 5 Å. C: half-map resolution within ±0.15 Å of 3.4/3.2; MolProbity score within ±0.1 of Table S1. D: r_Aα ≥ 0.6, r_Aβ ≤ 0.25. E: |Δ| ≤ 0.5 Å on the three headline distances. G: real-space CC(EMReady, deposited) > 0.95. Negative or partial outcomes are reported as such.

Do-not-run list: anything on lane `Cannon` / partition `gpu`; synthetic particles or CryoSPARC particle-level jobs; Tier 2 MD; model rebuilding beyond the optional sanity macrocycle; sourcing `sfGFP/env/activate.sh`; writing to `/n/lab_storage`; `module load cuda` in any CryoSPARC lane; leaving source trees or tarballs on holylabs after installs; installing anything outside `/n/holylabs/hekstra_lab/Everyone/ianpoon/software/`.

## 7. Execution order

1. Step 0: copy this plan to `cryoEM/gpt-plans/pp2a-day2025-reproduction.md`; create repo and netscratch dirs; copy the 26 MB cryosparc-tools/gemmi `pylibs` to `…/software/cs_pylibs`; `cluster-health --quick`. Ian (optional, before step 3): AlphaFold Server CIFs into `inputs/af3/`; NAMD/VMD tarballs into `inputs/software/`.
2. Submit the software build job (§5A) and, in parallel, `00`–`02` (downloads, model prep, tleap systems).
3. `03`–`04`: choose engine (NAMD if supplied, else OpenMM); 1 ns benchmark → two 100 ns chains on `gpu_test` (1.5–3.5 days unattended).
4. While MD runs: Ian runs `launchcryosparc`; `20`–`21` (CryoSPARC), `22` (RELION), `23` (EMReady), `24` (Phenix), `25` (geometry), `26`–`27` (Figs 2, 3, S5).
5. When MD finishes: `10`–`15` (Fig 1A–C, videos).
6. `40`, `41`, `report.md`; `stopcryosparc`.

## 8. Verification

- Installs: `relion_postprocess --version`, `gmx --version`, `chimerax --version` + offscreen PNG test, EMReady help text; no leftover source/tarballs; `quota /n/holylabs/hekstra_lab` before/after recorded in `DATA.md`.
- MD: chain reaches 100 ns per isoform; temperature/pressure/volume traces stable after release; protein-only trajectory has 1,000 frames; `seff` shows one busy core, `jobstats` a busy GPU; Fig 1A/1B compared side by side with the paper; videos play (≈20 s, 30 fps).
- RELION/Phenix/CryoSPARC: resolutions and validation numbers vs Table S1; CryoSPARC jobs "completed" in the GUI and recorded in `cs_state.json`; RELION and CryoSPARC local-res agree within ~0.3 Å.
- Geometry: HEAT CoM and salt-bridge numbers vs 62.7/59.6 and 13.2/12.4 Å; contact table; Fig 2B ratio table.
- Housekeeping: `cluster-health --quick` before and after; `DATA.md` lists every netscratch path with the purge date; repo mirror < 50 MB; `software_versions.json` present.

## 9. Risks

1. RELION 3.1 or GROMACS fail to compile with gcc 12 → gcc 10.2 module; RELION 4.0 fallback (same CLI).
2. ChimeraX offscreen rendering fails on the node (OSMesa/GLIBC) → run it inside Ian's VDI session with a display, or fall back to PyMOL (in dnaspring) for that figure.
3. EMReady torch/CUDA mismatch on MIG → CPU mode (slow but 256³ maps are small) or drop question G.
4. NAMD/VMD tarballs not supplied → OpenMM engine, MDAnalysis RMSD, ChimeraX movies (documented deviations).
5. CryoSPARC v5 job-type names → asserted from `job_register`; create one by hand in the GUI if needed.
6. MIG throughput lower than estimated → the 1 ns benchmark sets the chain; 4 fs HMR (OpenMM) fallback only if < 25 ns/day.
7. Box inflation from an extended AF tail → principal-axis orientation and the < 200k-particle assert; longer chain rather than reduced padding.
8. Fig 2B helix region assumption (9MIP 404–449) → verified by inspection; 3FGA/7K36 fallback.
9. Labmate's Phenix install moves → SBGrid ticket as the backup for Phenix.
10. Netscratch purge → only DCDs/maps live there; mirrored artefacts small and regenerable; `DATA.md` carries the purge date.
