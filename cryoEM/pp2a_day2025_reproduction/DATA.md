# Where the data lives

Everything bulky is on **netscratch**, which is purged 90 days after a file's last modification
(files written 2026-09-26 disappear around **2026-12-25** unless regenerated). Nothing here is
irreplaceable: every file below is produced by a numbered script in `scripts/` from public
downloads. The repo keeps only code, small CSV/JSON results and figures. Storage rules follow
`/n/home05/ianpoon/bin/hekstra_labCLUSTER_guide.md` (the `hekstra-cluster` skill).

| Path (netscratch, `$PP2A_DATA` = `/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025`) | Contents | Made by |
|---|---|---|
| `inputs/manifest.json`, `inputs/af/` | download manifest with md5s; AlphaFold DB v6 models AF-P30153-F1, AF-P30154-F1 (pdb+cif) and UniProt FASTAs | `00_fetch_inputs.py` |
| `em/emdb/{Aa,Ab}/` | EMD-48300 / EMD-48224 primary map, half maps 1+2, deposited FSC XML, EMDB header XML (gz + unpacked, ~350 MB) | `00_fetch_inputs.py` |
| `em/models/` | 9MIP, 9MF5, 3FGA (cif + pdb) | `00_fetch_inputs.py` |
| `em/phenix/{Aa,Ab}/` | `phenix.mtriage` logs/XML FSCs and `phenix.molprobity` output for the deposited models | `24_phenix_validate.sh` |
| `em/relion/{Aa,Ab}/` | `raw.mrc` (half-map average), `mask.mrc`, `postprocess.*` (gold-standard FSC, B-factor-sharpened map), `locres_locres.mrc` (= `relion_locres.mrc`), `locres_locres_filtered.mrc`, `raw_lp6.mrc`, `sharp_lp6.mrc`, `mask_sweep/` | `22_relion_postprocess.sh`, `22b_relion_mask_sweep.sh` |
| `em/cryosparc/` | copies of the CryoSPARC job outputs: local-resolution maps, locally filtered maps, Validation (FSC) curves and masks | `21_cryosparc_localres.py` |
| `em/emready/{Aa,Ab}/` | EMReady v1.1 re-sharpening of `raw.mrc` and `compare.json` vs the deposited map | `23_emready.sh` |
| `em/geometry/` | `geometry.json`, `fig2b_density.json` | `25_geometry.py`, `25b_fig2b_density.py` |
| `em/chimerax/` | ChimeraX command files and logs for every rendered panel | `26_chimerax_figures.py` |
| `cs_state/cs_state.json` | CryoSPARC project/workspace/job UIDs (P5 / W1 / J1–J14) | `21_cryosparc_localres.py` |
| `md/systems/{Aa,Ab}/` | `protein.pdb`, `protein_oriented.pdb`, `system.prmtop/inpcrd/pdb`, tleap logs; `build.json`, `prepare_models.json` | `01_prepare_models.py`, `02_build_systems.py` |
| `md/runs/<sys>_bench1ns/`, `md/runs/<sys>_rep1/` | OpenMM stage states, `traj.dcd` (100 ps frames), `04_production.log`, `state.chk`, `progress.json`, `gpu_usage.csv` | `slurm/pp2a_md_chain.sbatch` |

CryoSPARC project directory: `/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-pp2a-day2025`
(project **P5**, workspace W1 "deposited maps"; the instance's database is in `~/cryosparc_db`).

Slurm logs: `/n/netscratch/hekstra_lab/Lab/ian_poon/logs/pp2a-*`.

## Software installed for this reproduction (holylabs, persistent, no backup)

`/n/holylabs/hekstra_lab/Everyone/ianpoon/software/`

| Directory | What | Size | Installed by |
|---|---|---|---|
| `relion-3.1/` | RELION 3.1.4 (CPU, no GUI, own FFTW; gcc 12 + OpenMPI 4.1.5 runtime) | 113 MB | `slurm/pp2a_build_software.sbatch` |
| `gromacs-2023.3/` | GROMACS 2023.3 (CPU, AVX2_256, own FFTW) | 34 MB | `slurm/pp2a_build_software.sbatch` |
| `chimerax/` + `containers/centos-stream9.sif` | ChimeraX 1.12 el9 rpm unpacked + Stream 9 image + `sysroot/` of runtime rpms (OSMesa etc.) | 2.0 GB + 63 MB | `slurm/pp2a_install_chimerax.sbatch` |
| `emready/` | EMReady v1.1 + venv (torch 2.14 cu130) | 7.2 GB | `slurm/pp2a_install_emready.sbatch` |
| `cs_pylibs/` | cryosparc-tools 5.0.3 + gemmi for CryoSPARC's worker python | 26 MB | copied from the earlier tests |

Phenix 2.2 is used read-only from a labmate's install: `/n/lab_storage/hekstra_lab/people/aldama/software/phenix2/`.
Nothing was written to `/n/lab_storage`.
