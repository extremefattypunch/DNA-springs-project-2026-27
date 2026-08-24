#!/bin/bash
# Several replicates of one system, sharing a single GPU, in one job.
#
# Why this exists
# ---------------
# `Gres/gpu=1` is billed for a whole device however little of it you use, and the GPU is
# 99.3 % of a job's bill (836.5 of 842 on `gpu`).  So the one change that meaningfully
# reduces fairshare, short of not running at all, is to get more science out of each
# GPU-second -- and our jobs do not come close to filling an A100.  U0_apo is 19,497
# atoms and runs at 1,240 ns/day; S0_wt is 44,016 and runs at 751.  A card with 80 GB and
# 108 SMs is not the bottleneck for a system that size; kernel launch latency is.
#
# Running N replicates concurrently on one card therefore costs 418.25/hour once instead
# of N times.  Whether that is a real saving depends entirely on how throughput scales:
#
#   * if each of 3 replicates keeps 60 % of its solo speed, aggregate is 1.8x and the
#     cost per nanosecond falls by 1.8x;
#   * if each drops to exactly 1/3, elapsed time triples and the cost is IDENTICAL --
#     nothing gained, nothing lost;
#   * it cannot be worse than break-even on cost, only on wall clock.
#
# THIS HAS NOT YET BEEN MEASURED ON THIS PROJECT.  Do not assume the factor.  Run this
# script once on gpu_test (free) with PACK=3 and compare the reported aggregate against
# the solo ns/day in the run logs; the summary at the end prints exactly that comparison.
# Only then use it for billed work, with the measured factor in hand.
#
# Host resources scale with the pack: one CPU worker per replicate plus one for I/O, and
# host RAM per replicate.  GPU memory is the constraint to watch -- 3 x ~2 GB against a
# 20 GB MIG slice on gpu_test is comfortable, but a 211k-atom system is not 2 GB.
#
#   SYS=S0_wt PACK=3 NS=50 sbatch slurm/md_pack.sbatch
#
#SBATCH --job-name=sfgfp-pack
#SBATCH --output=slurm/logs/pack_%x_%j.out
#SBATCH --error=slurm/logs/pack_%x_%j.err
#SBATCH --partition=gpu_test
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p slurm/logs
source ./env/activate.sh

: "${SYS:?set SYS, e.g. SYS=S0_wt}"
PACK="${PACK:-3}"
NS="${NS:-50}"
FIRST_REP="${FIRST_REP:-1}"
SYSDIR="build/systems/$SYS"

case "${SLURM_JOB_PARTITION:-}" in
    gpu_test|gpu_requeue|serial_requeue|"") ;;
    *)  if [[ "${ALLOW_EXPENSIVE_GPU:-0}" != "1" ]]; then
            echo "== REFUSING ${SLURM_JOB_PARTITION}: price it with slurm/cost.py; use" >&2
            echo "   gpu_test (0) or gpu_requeue (418.25), or set" >&2
            echo "   ALLOW_EXPENSIVE_GPU=1 deliberately." >&2
            exit 2
        fi ;;
esac

echo "== $(date -Is) $(hostname) part=${SLURM_JOB_PARTITION:-none} job=${SLURM_JOB_ID:-none}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1
echo "== SYS=$SYS PACK=$PACK reps ${FIRST_REP}..$((FIRST_REP + PACK - 1)) NS=$NS"
echo "== all $PACK replicates share ONE GPU, so this job bills one GPU, not $PACK."

# CUDA MPS lets concurrent contexts interleave kernels instead of time-slicing, which is
# where most of the packing gain comes from for small systems.  Started if available;
# harmless to skip, the pack just runs a little slower.
export CUDA_MPS_PIPE_DIRECTORY="${TMPDIR:-/tmp}/mps-$SLURM_JOB_ID"
export CUDA_MPS_LOG_DIRECTORY="$CUDA_MPS_PIPE_DIRECTORY/log"
mkdir -p "$CUDA_MPS_LOG_DIRECTORY"
if command -v nvidia-cuda-mps-control >/dev/null 2>&1; then
    nvidia-cuda-mps-control -d >/dev/null 2>&1 && echo "== CUDA MPS started" \
        || echo "== CUDA MPS unavailable (fine; contexts will time-slice)"
fi
stop_mps() { echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true; }

pids=(); reps=()
for i in $(seq 0 $((PACK - 1))); do
    REP=$((FIRST_REP + i))
    RUNDIR="data/runs/$SYS/rep$REP"
    mkdir -p "$RUNDIR"
    reps+=("$REP")
    (
        set -euo pipefail
        COMMON=(--sysdir "$SYSDIR" --outdir "$RUNDIR")
        if [[ ! -f "$RUNDIR/minimised.xml" ]]; then
            "$DNASPRING_PY" runs/openmm_engine.py minimise "${COMMON[@]}"
        fi
        if [[ ! -f "$RUNDIR/equilibrated.xml" ]]; then
            "$DNASPRING_PY" runs/openmm_engine.py equilibrate "${COMMON[@]}"
        fi
        "$DNASPRING_PY" runs/openmm_engine.py produce "${COMMON[@]}" --ns "$NS"
    ) > "data/runs/$SYS/rep$REP/pack.log" 2>&1 &
    pids+=($!)
    echo "== launched rep$REP as pid ${pids[-1]}"
    sleep 5                      # stagger context creation; simultaneous init can OOM
done

# Sample the shared card, which is the whole point: if utilisation is still well short of
# 100 % with the pack running, the pack is too small.
USAGE="data/runs/$SYS/pack_usage_${SLURM_JOB_ID:-local}.csv"
echo "iso,gpu_util_pct,gpu_mem_used_MiB,gpu_mem_total_MiB" > "$USAGE"
( while sleep 60; do
    g=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total \
        --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    [[ -n "$g" ]] && echo "$(date -Is),$g" >> "$USAGE"
  done ) & SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true; stop_mps' EXIT

fail=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        echo "== rep${reps[$i]} finished"
    else
        echo "== rep${reps[$i]} FAILED (see data/runs/$SYS/rep${reps[$i]}/pack.log)" >&2
        fail=1
    fi
done
kill "$SAMPLER" 2>/dev/null || true

echo
echo "== packing efficiency, measured =="
"$DNASPRING_PY" - "$SYS" "$PACK" "$USAGE" <<'PYEOF'
import glob, os, sys
sysname, pack, usage = sys.argv[1], int(sys.argv[2]), sys.argv[3]

def speed(log):
    try:
        last = None
        for line in open(log):
            if line[:1].isdigit():
                last = line
        return float(last.split()[7]) if last else None
    except (OSError, IndexError, ValueError, TypeError):
        return None

sp = [s for s in (speed(p) for p in
                  sorted(glob.glob(f"data/runs/{sysname}/rep*/03_production/production.log")))
      if s]
if sp:
    agg = sum(sp[:pack])
    print(f"  per-replicate ns/day in this pack: "
          + ", ".join(f"{s:.0f}" for s in sp[:pack]))
    print(f"  aggregate {agg:.0f} ns/day from ONE billed GPU")
    print(f"  Compare this against the same system's solo ns/day. If aggregate exceeds")
    print(f"  solo, the ratio is exactly the factor by which packing cuts cost per ns.")
if os.path.exists(usage):
    rows = [l.split(",") for l in open(usage).read().splitlines()[1:] if l]
    if rows:
        u = [float(r[1]) for r in rows]
        m = max(float(r[2]) for r in rows) / 1024
        tot = max(float(r[3]) for r in rows) / 1024
        print(f"  GPU utilisation {sum(u)/len(u):.0f}% mean, {max(u):.0f}% peak; "
              f"GPU memory {m:.1f} of {tot:.1f} GB")
        if sum(u)/len(u) < 85:
            print(f"  Still short of saturation -- a larger PACK would divide the same "
                  f"GPU bill further.")
PYEOF
echo "== $(date -Is) done (fail=$fail)"
exit "$fail"
