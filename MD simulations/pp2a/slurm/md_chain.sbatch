#!/bin/bash -l
# One system, one replicate, three stages (minimise -> equilibrate -> produce), restart-safe.
# The sfGFP chain (../sfGFP/slurm/md_chain.sbatch) with three changes for this project:
#   * run directories are ABSOLUTE on netscratch ($PP2A_MD_SCRATCH/runs/<TAG>/rep<N>);
#     the home checkout has no data/ symlink and must not grow trajectories;
#   * logs go to netscratch, as the cluster guide requires;
#   * partition, time and constraint come from the submitter, so the same script runs a
#     free gpu_test benchmark, a half-price gpu_requeue production chunk, or (explicitly
#     allowed) a full-price gpu job.
#
# Cost, so nobody is surprised: Gres/gpu billing weights are gpu_test 0, gpu_requeue
# 418.25, gpu 836.5, gpu_h200 2648.8 per GPU-hour.  The engine checkpoints every 200 ps
# and resumes from state.chk, so preemption on gpu_requeue loses at most 200 ps.
#
#   SYS=P3_ss50 NS=100 REP=1 sbatch -p gpu_requeue -t 18:00:00 slurm/md_chain.sbatch
#   SYS=D_ss50 NS=20 RESTRAINT_X0=80 RESTRAINT_K=33 sbatch -p gpu_test --constraint=a100-mig -t 8:00:00 ...
#
#SBATCH --job-name=pp2amd
#SBATCH --output=/n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH --error=/n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
#SBATCH --partition=gpu_test
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
set -euo pipefail
source "/n/home05/ianpoon/DNA-springs-project-2026-27/MD simulations/pp2a/env/activate.sh"
cd "$PP2A_MD_REPO"

: "${SYS:?set SYS, e.g. SYS=P3_ss50}"
REP="${REP:-1}"
NS="${NS:-100}"
SYSDIR="build/systems/$SYS"
TAG="$SYS${CLAMP_PN:+_${CLAMP_PN}pN}${RESTRAINT_X0:+_x${RESTRAINT_X0}A}${RESTRAINT_K:+_k${RESTRAINT_K}}${RUN_SUFFIX:-}"
case "$TAG" in
    *_x_*|*_xA_*|*_x_A*|*_k|*_xA_k|*pN_*pN*)
        echo "== malformed TAG '$TAG': a condition variable expanded empty." >&2; exit 2 ;;
esac
[[ "$REP" =~ ^[0-9]+$ ]] || { echo "== REP must be a number, got '$REP'" >&2; exit 2; }
RUNDIR="$PP2A_MD_SCRATCH/runs/$TAG/rep$REP"
mkdir -p "$RUNDIR"
exec 9>"$RUNDIR/.lock"; flock -n 9 || { echo "== another job already owns $RUNDIR; exiting"; exit 0; }
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; export TMPDIR="$TMP"
trap 'rm -rf "$TMP"' EXIT

echo "== $(date -Is) $(hostname) part=${SLURM_JOB_PARTITION:-none} job=${SLURM_JOB_ID:-none}"
nvidia-smi --query-gpu=name,memory.total,mig.mode.current --format=csv,noheader | head -1
echo "== SYS=$SYS TAG=$TAG REP=$REP NS=$NS  restart_count=${SLURM_RESTART_COUNT:-0}"

case "${SLURM_JOB_PARTITION:-}" in
    gpu_test|"") ;;
    gpu_requeue|serial_requeue)
        echo "== billed partition ${SLURM_JOB_PARTITION} (Gres/gpu=418.25), preemptible; checkpoint-resume makes that safe." ;;
    gpu)
        if [[ "${ALLOW_EXPENSIVE_GPU:-0}" != "1" ]]; then
            echo "== REFUSING ${SLURM_JOB_PARTITION} (836.5 per GPU-hour, twice gpu_requeue). Resubmit with ALLOW_EXPENSIVE_GPU=1 if it is really wanted." >&2
            exit 2
        fi
        echo "== WARNING: gpu partition at 836.5 per GPU-hour, explicitly allowed." ;;
    *)
        echo "== REFUSING ${SLURM_JOB_PARTITION}: not a partition this project uses." >&2; exit 2 ;;
esac

USAGE="$RUNDIR/resource_usage.csv"
[[ -s "$USAGE" ]] || echo "iso,gpu_util_pct,gpu_mem_used_MiB,gpu_mem_total_MiB,host_rss_MiB" > "$USAGE"
(
  while sleep 60; do
    g=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    r=$(ps -o rss= -u "$USER" 2>/dev/null | awk '{s+=$1} END {printf "%d", s/1024}')
    [[ -n "$g" ]] && echo "$(date -Is),$g,$r" >> "$USAGE"
  done
) & SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true; rm -rf "$TMP"' EXIT

CLAMP_ARGS=()
if [[ -n "${RESTRAINT_X0:-}" ]]; then
    CLAMP_ATOMS="$(cat "$SYSDIR/clamp_atoms.txt")"
    CLAMP_ARGS=(--restraint-x0 "$RESTRAINT_X0" --clamp-atoms "$CLAMP_ATOMS")
    [[ -n "${RESTRAINT_K:-}" ]] && CLAMP_ARGS+=(--restraint-k "$RESTRAINT_K")
    [[ -n "${RESTRAINT_RAMP_PS:-}" ]] && CLAMP_ARGS+=(--restraint-ramp-ps "$RESTRAINT_RAMP_PS")
    echo "== harmonic window x0 = ${RESTRAINT_X0} A, k = ${RESTRAINT_K:-default} between atoms ${CLAMP_ATOMS}"
elif [[ -n "${CLAMP_PN:-}" ]]; then
    CLAMP_ATOMS="$(cat "$SYSDIR/clamp_atoms.txt")"
    CLAMP_ARGS=(--clamp-pN "$CLAMP_PN" --clamp-atoms "$CLAMP_ATOMS")
    echo "== force clamp ${CLAMP_PN} pN between atoms ${CLAMP_ATOMS}"
fi

COMMON=(--prmtop "$SYSDIR/system.prmtop" --inpcrd "$SYSDIR/system.inpcrd" --platform CUDA)
if [[ ! -f "$RUNDIR/01_minimise/minimised.xml" ]]; then
    "$DNASPRING_PY" runs/openmm_engine.py minimise "${COMMON[@]}" --outdir "$RUNDIR/01_minimise"
fi
if [[ ! -f "$RUNDIR/02_equilibrate/equilibrated.xml" ]]; then
    "$DNASPRING_PY" runs/openmm_engine.py equilibrate "${COMMON[@]}" \
        --outdir "$RUNDIR/02_equilibrate" --seed "$((1000 + REP))"
fi
"$DNASPRING_PY" runs/openmm_engine.py produce "${COMMON[@]}" \
    --outdir "$RUNDIR/03_production" --ns "$NS" --seed "$((2000 + REP))" "${CLAMP_ARGS[@]}"

if [[ -s "$USAGE" ]]; then
    awk -F, 'NR>1 {n++; u+=$2; if($2>mu)mu=$2; if($3>mg)mg=$3; if($5>mr)mr=$5}
             END {if(n) printf "== used: GPU %.0f%% mean / %.0f%% peak, GPU mem %.1f GB peak of %.1f, host RSS %.1f GB peak (%d samples)\n", u/n, mu, mg/1024, $4/1024, mr/1024, n}' "$USAGE"
fi
echo "== $(date -Is) done"
