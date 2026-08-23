#!/bin/bash
# One system, one replicate: minimise -> equilibrate -> produce, restart-safe.
#
# Resource notes (why these numbers):
#  * PARTITION: gpu_test by default, and this is not a detail.  `scontrol show
#    partition` gives TRESBillingWeights Gres/gpu=836.5 on `gpu` and Gres/gpu=0 on
#    `gpu_test` -- the billed partition charges the whole lab's fairshare and the test
#    one is free.  230 GPU-hours of this project's exploratory sweeps went onto `gpu`
#    in a single day and exhausted the lab's allocation, which is a cost borne by
#    everyone else in it.  Development, parameter scans and anything you might rerun
#    belong on gpu_test.  Move to `gpu` only for a final production set you are
#    confident in, and check `sshare -U -u $USER` before you do.
#  * The trade is a 12 h wall limit and a100_3g.20gb MIG slices rather than whole
#    80 GB cards -- roughly a third of the throughput.  50 ns of a 120k-atom system
#    still fits; the restart machinery below covers anything that does not.
#  * --gres is therefore generic.  It used to name the A100-80GB type explicitly to keep
#    replicate throughput comparable, which is a real concern on a mixed pool but not
#    one worth an entire lab's fairshare.
#  * --cpus-per-task=4 is plenty: OpenMM's CUDA platform uses one worker thread
#    plus PME bookkeeping. Asking for more only delays scheduling.
#  * 3-day cap is the partition maximum on gpu/gpu_requeue. The lab's older
#    launch_production.sh asks for 7 days, which no GPU partition now allows.
#  * --requeue with append-mode logs: production checkpoints every 10 report
#    intervals and resumes from state.chk, so preemption costs at most that much.
#SBATCH --job-name=sfgfp-md
#SBATCH --output=slurm/logs/md_%x_%A_%a.out
#SBATCH --error=slurm/logs/md_%x_%A_%a.err
#SBATCH --partition=gpu_test
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
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
# The force-clamp ladder runs six forces off ONE topology, so the force has to be in
# the run directory name.  Without it every clamp point writes to data/runs/S6_clamp
# and six different physical conditions overwrite each other's trajectory.
# The same reasoning applies to the force-extension windows: several separations run
# off ONE topology, so x0 has to be in the directory name too.
TAG="$SYS${CLAMP_PN:+_${CLAMP_PN}pN}${RESTRAINT_X0:+_x${RESTRAINT_X0}A}${RESTRAINT_K:+_k${RESTRAINT_K}}"
RUNDIR="data/runs/$TAG/rep$REP"
mkdir -p "$RUNDIR"

echo "== $(date -Is) $(hostname) part=${SLURM_JOB_PARTITION:-none} job=${SLURM_JOB_ID:-none}"
nvidia-smi --query-gpu=name,memory.total,mig.mode.current --format=csv,noheader | head -1
echo "== SYS=$SYS TAG=$TAG REP=$REP NS=$NS  restart_count=${SLURM_RESTART_COUNT:-0}"

CLAMP_ARGS=()
if [[ -n "${RESTRAINT_X0:-}" ]]; then
    # Force-extension window: hold the two ends of the spring at a fixed separation and
    # let the mean restraint force report what the force field says the spring pushes
    # with.  Same anchor atoms as the clamp ladder, written by the system builder.
    CLAMP_ATOMS="$(cat "$SYSDIR/clamp_atoms.txt")"
    CLAMP_ARGS=(--restraint-x0 "$RESTRAINT_X0" --clamp-atoms "$CLAMP_ATOMS")
    if [[ -n "${RESTRAINT_K:-}" ]]; then
        CLAMP_ARGS+=(--restraint-k "$RESTRAINT_K")
    fi
    # The deepest compression windows have to be walked in from a built geometry 40 A
    # away, so they get a longer ramp than the 2 ns default; going faster risks tearing a
    # base pair on the way in, and a pair torn by the approach looks exactly like the
    # kink these runs exist to look for.
    if [[ -n "${RESTRAINT_RAMP_PS:-}" ]]; then
        CLAMP_ARGS+=(--restraint-ramp-ps "$RESTRAINT_RAMP_PS")
    fi
    echo "== harmonic window x0 = ${RESTRAINT_X0} A, k = ${RESTRAINT_K:-default} \
between atoms ${CLAMP_ATOMS}"
elif [[ -n "${CLAMP_PN:-}" ]]; then
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
