#!/bin/bash -l
# Generic CPU wrapper for the cryo-EM analysis steps (RELION post-processing, Phenix
# validation, geometry, figure rendering).  Usage:
#   sbatch [-J name] [-c N] [--mem=..] [-t ..] slurm/pp2a_cpu.sbatch <command...>
# `test` partition: free (all billing weights 0), 12 h cap, PreemptMode=OFF.
# Defaults are sized for relion_postprocess --locres on a 256^3 map (single-threaded
# per map, ~2 GB): request 1 core unless the command is known to thread.
#SBATCH -J pp2a-cpu
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-02:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
source /etc/profile.d/lmod.sh
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/env/activate_pp2a.sh
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT; export TMPDIR="$TMP"
cd "$PP2A_REPO"
echo "[$(date)] $*"
"$@"
echo "[$(date)] done"
