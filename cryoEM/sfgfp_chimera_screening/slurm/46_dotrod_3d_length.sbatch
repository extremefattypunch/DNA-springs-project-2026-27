#!/bin/bash -l
#SBATCH -J dotrod3d
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
# W4 3D attached-DNA length: analysis -> ChimeraX renders + sfGFP fit -> finalize.
#   STAGES="analyse poses render fit finalize" (default all) sbatch slurm/46_dotrod_3d_length.sbatch
set -e
REPO=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening
PY="/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python -I"
source $REPO/env/activate_sfgfp.sh
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT; export TMPDIR="$TMP"
export OMP_NUM_THREADS=4
cd $REPO
for st in ${STAGES:-analyse poses render fit finalize}; do
  echo "== $st $(date)"
  case $st in
    analyse)  $PY scripts/46_dotrod_3d_length.py analyse ;;
    poses)    $PY scripts/46_dotrod_3d_length.py poses ;;
    render)   $CHIMERAX --nogui --offscreen --silent scripts/46_dotrod_render_fit.cxc ;;
    fit)      $CHIMERAX --nogui --offscreen scripts/46_dotrod_fit_sfgfp.cxc > $SFGFP_DATA/dotrod_3d_length/chimerax_fit_log.txt 2>&1 ;;
    finalize) $PY scripts/46_dotrod_3d_length.py finalize ;;
  esac
done
echo "== done $(date)"
