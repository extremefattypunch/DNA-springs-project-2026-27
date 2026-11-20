#!/bin/bash -l
#SBATCH -J figS7
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
# Fig S7 analogue: prep display maps -> ChimeraX renders -> montage.  STAGES="prep render montage" (default all)
set -e
REPO=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening
PY="/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python -I"
source $REPO/env/activate_sfgfp.sh
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT; export TMPDIR="$TMP"
export OMP_NUM_THREADS=4
cd $REPO
for st in ${STAGES:-prep render montage}; do
  echo "== $st $(date)"
  case $st in
    prep)    $PY scripts/50_figS7_class3d.py prep ;;
    render)  $CHIMERAX --nogui --offscreen scripts/50_figS7_class3d.cxc > /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/figS7/chimerax_log.txt 2>&1 ;;
    montage) $PY scripts/50_figS7_class3d.py montage ;;
  esac
done
echo "== done $(date)"
