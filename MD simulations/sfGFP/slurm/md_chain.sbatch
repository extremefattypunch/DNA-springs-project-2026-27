#!/bin/bash
# One system, one replicate: minimise -> equilibrate -> produce, restart-safe.
#
# Resource notes (why these numbers):
#  * --gres names the A100-80GB GRES type explicitly. gpu_requeue and
#    serial_requeue also hand out a100_1g.5gb / a100_3g.20gb MIG slices and
#    16 GB V100s; any of them would *run* these systems, but the throughput
#    spread across them is several-fold and would make replicates
#    incomparable in wall-clock planning.
#  * --cpus-per-task=4 is plenty: OpenMM's CUDA platform uses one worker thread
#    plus PME bookkeeping. Asking for more only delays scheduling.
#  * 3-day cap is the partition maximum on gpu/gpu_requeue. The lab's older
#    launch_production.sh asks for 7 days, which no GPU partition now allows.
#  * --requeue with append-mode logs: production checkpoints every 10 report
#    intervals and resumes from state.chk, so preemption costs at most that much.
#SBATCH --job-name=sfgfp-md
#SBATCH --output=slurm/logs/md_%x_%A_%a.out
#SBATCH --error=slurm/logs/md_%x_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p slurm/logs
source ./env/activate.sh

: "${SYS:?set SYS, e.g. SYS=S0_wt}"
REP="${REP:-${SLURM_ARRAY_TASK_ID:-1}}"
NS="${NS:-50}"
SYSDIR="build/systems/$SYS"
RUNDIR="data/runs/$SYS/rep$REP"
mkdir -p "$RUNDIR"

echo "== $(date -Is) $(hostname) part=${SLURM_JOB_PARTITION:-none} job=${SLURM_JOB_ID:-none}"
nvidia-smi --query-gpu=name,memory.total,mig.mode.current --format=csv,noheader | head -1
echo "== SYS=$SYS REP=$REP NS=$NS  restart_count=${SLURM_RESTART_COUNT:-0}"

CLAMP_ARGS=()
if [[ -n "${CLAMP_PN:-}" ]]; then
    # Force-clamp ladder: the two attachment-site atom indices are written by the
    # system builder so the physics and the topology cannot drift apart.
    CLAMP_ATOMS="$(cat "$SYSDIR/clamp_atoms.txt")"
    CLAMP_ARGS=(--clamp-pN "$CLAMP_PN" --clamp-atoms "$CLAMP_ATOMS")
    echo "== force clamp ${CLAMP_PN} pN between atoms ${CLAMP_ATOMS}"
fi

COMMON=(--prmtop "$SYSDIR/system.prmtop" --inpcrd "$SYSDIR/system.inpcrd"
        --platform CUDA)

# Each stage is skipped if its output already exists, so a requeued job walks
# straight to the unfinished stage instead of redoing the cheap ones.
if [[ ! -f "$RUNDIR/01_minimise/minimised.xml" ]]; then
    "$DNASPRING_PY" runs/openmm_engine.py minimise "${COMMON[@]}" \
        --outdir "$RUNDIR/01_minimise"
fi
if [[ ! -f "$RUNDIR/02_equilibrate/equilibrated.xml" ]]; then
    "$DNASPRING_PY" runs/openmm_engine.py equilibrate "${COMMON[@]}" \
        --outdir "$RUNDIR/02_equilibrate" --seed "$((1000 + REP))"
fi
"$DNASPRING_PY" runs/openmm_engine.py produce "${COMMON[@]}" \
    --outdir "$RUNDIR/03_production" --ns "$NS" --seed "$((2000 + REP))" \
    "${CLAMP_ARGS[@]}"

echo "== $(date -Is) done"
