# sfGFP + 50 bp DNA spring chimera — cryo-EM of the 12-11-2026 screening session

CryoSPARC v5 single-particle processing of 9,019 Falcon/EPU movies (0.71 Å/px, 50 fractions, 39.7 e/Å²) of
sfGFP clicked to a closed 50 bp DNA loop, following the Day et al. 2025 (Structure) Fig S6 workflow, with the
same outputs (Fig S5–S7, Table S1) as `../pp2a_day2025_reproduction/`. Plan: `~/.claude/plans/i-want-you-to-valiant-ullman.md`.
`report.md` will hold the findings; `DATA.md` says where the bulk data is and when it purges.

```
env/activate_sfgfp.sh     paths (sources the PP2A env for dnaspring python, RELION, ChimeraX, EMReady, Phenix, CS python)
scripts/00_inventory.py   facility listing -> results/raw_inventory.csv
scripts/download_login.sh rsync from the facility (LOGIN NODE ONLY; compute nodes have no route)
model/                    L6_loop50 atomic model (fitting/validation only, never a processing reference)
```
