#!/bin/bash -l
# Run finalize.sh (collect -> findings -> figures -> animations -> wiki page) as a batch job.
# The chimera trajectory is 259k atoms; at stride 2 a single module holds ~8 GB of
# coordinates, and PyMOL renders 80 frames x 3 views, so this does not belong on a login
# node.  `shared` CPU, no GPU: billing weight 1.0 per core-hour.
#
#   sbatch slurm/finalize.sbatch                 # everything
#   ONLY="P0_apo" ANIMATE="P0_apo" sbatch slurm/finalize.sbatch
#
#SBATCH -J pp2amd-finalize
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 6:00:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
cd "/n/home05/ianpoon/DNA-springs-project-2026-27/MD simulations/pp2a"
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; export TMPDIR="$TMP"; trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=4
bash finalize.sh
