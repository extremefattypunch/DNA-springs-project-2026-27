#!/bin/bash
# Plan section 4 analysis - per-residue RMSF with GROMACS 2023.3 (the paper's RMSF tool), Fig 1A.
# Uses the protein-only xtc + reference pdb written by 10_md_rmsd_rmsf.py.  `-fit` superposes each frame on the
# reference before computing fluctuations (all Calpha), `-res` averages per residue.  Output: results/md_rmsf_gmx_<sys>.csv
set -euo pipefail
: "${PP2A_DATA:?source env/activate_pp2a.sh}"; command -v gmx >/dev/null || { echo "gmx not on PATH"; exit 1; }
for s in "${@:-Aa Ab}"; do
  RUN=$PP2A_DATA/md/runs/${s}_rep1; [ -f "$RUN/protein.xtc" ] || { echo "$s: run 10_md_rmsd_rmsf.py first"; continue; }
  cd "$RUN"
  # index group 3 = C-alpha in gmx's default groups for a protein-only structure
  echo "C-alpha" | gmx rmsf -s protein_ref.pdb -f protein.xtc -o rmsf_gmx.xvg -res -fit yes 2> rmsf_gmx.log >/dev/null
  grep -vE "^[@#]" rmsf_gmx.xvg | awk 'BEGIN{print "resid,rmsf_A"}{printf "%d,%.3f\n",$1,$2*10}' > "$PP2A_REPO/results/md_rmsf_gmx_${s}.csv"
  echo "$s: $(wc -l < "$PP2A_REPO/results/md_rmsf_gmx_${s}.csv") residues -> results/md_rmsf_gmx_${s}.csv"
done
