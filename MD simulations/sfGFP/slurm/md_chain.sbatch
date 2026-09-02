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
#  * WHEN A BILLED PARTITION IS UNAVOIDABLE, IT IS `gpu_requeue`, NOT `gpu`.  Measured
#    from `scontrol show partition`, Gres/gpu weights are: gpu_test 0, gpu_requeue
#    418.25, gpu 836.5, seas_gpu 1727.8, gpu_h200 2648.8.  gpu_requeue is exactly half
#    of `gpu` and its pool includes H100 and H200 cards, so `gpu_h200` charges 6.3x for
#    hardware gpu_requeue hands out at 418.25.  The premium buys only immunity from
#    preemption, and the produce step checkpoints every 10 report intervals and resumes
#    by appending -- so a preemption costs a few hundred picoseconds of wall clock and
#    no work at all.  There is no good reason for this project to use `gpu`.
#  * --cpus-per-task=1.  OpenMM's CUDA platform is single-threaded on the host, so one
#    run needs exactly one core.  FAS RC's Job Defense Shield flagged 157 of my jobs for
#    this: every one showed CPU utilisation of almost exactly 100%/cores -- 24.9% on four
#    cores -- which is the signature of a single-threaded process, and it wasted about
#    1,032 CPU-hours.  A PACKED job is the one exception and wants one core per concurrent
#    run; see md_pack.sbatch.
#  * --constraint=a100.  Billing is FLAT per GPU: gpu_requeue charges 418.25 whether it
#    hands you an A100-SXM4-80GB (2039 GB/s) or an A40 (696 GB/s).  MD throughput is set
#    by memory bandwidth, so taking whatever is free costs up to 2.9x more per nanosecond
#    for the same money.  Asking for an A100 is the cheapest thing in this file.  Relax it
#    if the queue is long -- but relax it knowingly.
#  * --mem=8G, not 32G.  Honestly, this barely matters:
#    sacct decomposes a job of the old shape as 836.5 (GPU) + 3.6 (4 cores) +
#    1.92 (32 GB) = 842, so the GPU is 99.3 % of the bill and right-sizing host
#    resources saves 0.5 %.  It is done because it is free, because it lets more jobs
#    share a node, and because 32G was a number I guessed rather than measured.  The
#    sampler below now records peak RSS and GPU utilisation so the next person sizes
#    this from evidence.  Note that asking for less GPU *memory* saves nothing
#    whatsoever: Gres/gpu=1 is charged for a whole device however little of it is used,
#    and this cluster publishes no per-GPU-type weight.
#  * --time is 12 h to match gpu_test's cap. It is worth setting accurately, but not for
#    the reason it seems: billing accrues on ELAPSED time, so an over-long limit costs
#    nothing directly (my jobs requested 2,880 GPU-hours, used 230.5, and were billed
#    842 x 230.5).  What an accurate limit buys is backfill -- the scheduler can slot a
#    3 h job into a gap it cannot slot a 3-day job into, so the job starts sooner.
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
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --constraint=a100
#
# GOTCHA -- the a100 feature name differs by partition.  gpu_test's nodes advertise
# `a100-mig` (they are a100_3g.20gb slices); only the billed partitions advertise `a100`
# for full cards.  `--constraint=a100` on gpu_test fails at submit time with "Requested
# node configuration is not available", which reads like a quota problem and is not.
# Override it on the command line:  sbatch --partition=gpu_test --constraint=a100-mig ...
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
# A tag with an empty condition in it is a submission bug, not a run.  On 2026-08-23 an
# invocation with RESTRAINT_X0, RESTRAINT_K and REP all empty produced
# data/runs/B_27bp_xA_k/rep -- 15 MB of minimisation under a name that matches the
# force-extension glob and means nothing.  Catch it at submission instead.
case "$TAG" in
    *_x_*|*_xA_*|*_x_A*|*_k|*_xA_k|*pN_*pN*)
        echo "== malformed TAG '$TAG': a condition variable expanded empty." >&2
        echo "   RESTRAINT_X0='${RESTRAINT_X0:-}' RESTRAINT_K='${RESTRAINT_K:-}'" >&2
        echo "   CLAMP_PN='${CLAMP_PN:-}' -- set them or leave them unset, not empty." >&2
        exit 2 ;;
esac
if ! [[ "$REP" =~ ^[0-9]+$ ]]; then
    echo "== REP must be a number, got '$REP'" >&2; exit 2
fi
RUNDIR="data/runs/$TAG/rep$REP"
mkdir -p "$RUNDIR"

echo "== $(date -Is) $(hostname) part=${SLURM_JOB_PARTITION:-none} job=${SLURM_JOB_ID:-none}"
nvidia-smi --query-gpu=name,memory.total,mig.mode.current --format=csv,noheader | head -1
echo "== SYS=$SYS TAG=$TAG REP=$REP NS=$NS  restart_count=${SLURM_RESTART_COUNT:-0}"

# Refuse to spend the lab's fairshare by accident.  The cost of a GPU job is invisible
# at submission -- squeue looks identical whether the partition weight is 0 or 2648.8 --
# which is how 120 jobs came to cost 194,084 billing-hours before anyone noticed.  An
# expensive partition now has to be asked for in as many words.
case "${SLURM_JOB_PARTITION:-}" in
    gpu_test|"") ;;
    gpu_requeue|serial_requeue)
        echo "== billed partition ${SLURM_JOB_PARTITION} (Gres/gpu=418.25), preemptible;"
        echo "   checkpoint-resume makes that safe. python slurm/cost.py prices it." ;;
    *)
        if [[ "${ALLOW_EXPENSIVE_GPU:-0}" != "1" ]]; then
            echo "== REFUSING to run on ${SLURM_JOB_PARTITION}." >&2
            echo "   gpu=836.5, seas_gpu=1727.8, gpu_h200=2648.8 per GPU-hour, against" >&2
            echo "   418.25 on gpu_requeue for the same or better hardware and 0 on" >&2
            echo "   gpu_test. Price it with 'python slurm/cost.py' first; if this really" >&2
            echo "   is the right partition, resubmit with ALLOW_EXPENSIVE_GPU=1." >&2
            exit 2
        fi
        echo "== WARNING: ${SLURM_JOB_PARTITION} at $( \
            scontrol show partition "$SLURM_JOB_PARTITION" | tr ' ' '\n' \
            | grep -o 'Gres/gpu=[0-9.]*' | head -1) per GPU-hour, explicitly allowed." ;;
esac

# Sample what this job actually uses, so the next request is measured rather than
# guessed.  32G and 4 cores were guesses; nothing in sacct recorded MaxRSS for these
# steps, so there was no evidence to correct them with.  Cheap: one line per minute.
USAGE="$RUNDIR/resource_usage.csv"
if [[ ! -s "$USAGE" ]]; then
    echo "iso,gpu_util_pct,gpu_mem_used_MiB,gpu_mem_total_MiB,host_rss_MiB" > "$USAGE"
fi
(
  while sleep 60; do
    g=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total \
        --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    r=$(ps -o rss= -u "$USER" 2>/dev/null | awk '{s+=$1} END {printf "%d", s/1024}')
    [[ -n "$g" ]] && echo "$(date -Is),$g,$r" >> "$USAGE"
  done
) & SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true' EXIT

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

# Summarise what was actually used, so it lands in the job log where it will be seen.
if [[ -s "$USAGE" ]]; then
    awk -F, 'NR>1 {n++; u+=$2; if($2>mu)mu=$2; if($3>mg)mg=$3; if($5>mr)mr=$5}
             END {if(n) printf "== used: GPU %.0f%% mean / %.0f%% peak, "\
                  "GPU mem %.1f GB peak of %.1f, host RSS %.1f GB peak "\
                  "(%d samples)\n", u/n, mu, mg/1024, $4/1024, mr/1024, n}' "$USAGE"
    echo "   -> if GPU mean is well under 100%, several replicates would share one card;"
    echo "      that divides the 418.25/hour between them, which is the only lever that"
    echo "      touches the 99.3% of the bill the GPU accounts for."
fi
echo "== $(date -Is) done"
