#!/bin/bash
# Several independent runs sharing ONE GPU, in one job.
#
# Why
# ---
# `Gres/gpu=1` is billed for a whole device however little of it you use, and the GPU is
# 99.3 % of a job's bill (836.5 of 842).  So the only change that meaningfully cuts
# fairshare, short of running less science, is to get more out of each GPU-second -- and
# our systems do not come close to filling an A100.  U0_apo is 19,497 atoms at
# 1,240 ns/day; S0_wt is 44,016 at 751.  For systems that size the card is idle waiting
# on kernel launches, not saturated.
#
# Running N runs concurrently on one card bills one GPU instead of N.  The worst case is
# break-even on cost (if each run drops to exactly 1/N speed, elapsed time grows by N and
# the bill is unchanged); it can only lose wall-clock, never fairshare.  The job prints
# the measured GPU utilisation and per-run throughput at the end so the factor is known
# rather than assumed.
#
# JOBS is the work list: one spec per line, whitespace-separated KEY=VALUE pairs.
# Recognised keys mirror md_chain.sbatch exactly, so a packed run is byte-identical to a
# chained one and the analysis finds it in the same place:
#
#   SYS       (required)  system name under build/systems/
#   REP       (required)  replicate number; sets the seeds, as in md_chain
#   NS                    nanoseconds of production (default $NS, itself default 30)
#   X0                    harmonic window target separation, A  -> RESTRAINT_X0
#   K                     harmonic window stiffness, kJ/mol/nm2 -> RESTRAINT_K
#   RAMP_PS               restraint approach ramp, ps           -> RESTRAINT_RAMP_PS
#   CLAMP_PN              constant-force clamp, pN
#
# Example -- four force-extension windows and two replicates each, one GPU:
#
#   JOBS=$'SYS=B_27bp REP=1 X0=50 K=200\nSYS=B_27bp REP=2 X0=50 K=200\n...' \
#     sbatch slurm/md_pack.sbatch
#
# Host resources scale with the pack: one CPU worker per run plus one for I/O, and host
# RAM per run.  GPU memory is the real ceiling -- check it against the card before
# choosing N, especially for the 211k-atom systems.
#
#SBATCH --job-name=sfgfp-pack
#SBATCH --output=slurm/logs/pack_%x_%j.out
#SBATCH --error=slurm/logs/pack_%x_%j.err
#SBATCH --partition=gpu_test
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p slurm/logs
source ./env/activate.sh

: "${JOBS:?set JOBS to a newline-separated list of 'SYS=.. REP=.. [X0=..] ..' specs}"
NS_DEFAULT="${NS:-30}"

case "${SLURM_JOB_PARTITION:-}" in
    gpu_test|gpu_requeue|serial_requeue|"") ;;
    *)  if [[ "${ALLOW_EXPENSIVE_GPU:-0}" != "1" ]]; then
            echo "== REFUSING ${SLURM_JOB_PARTITION}: gpu=836.5, seas_gpu=1727.8," >&2
            echo "   gpu_h200=2648.8 per GPU-hour against 418.25 on gpu_requeue for the" >&2
            echo "   same or better hardware, and 0 on gpu_test. Price it with" >&2
            echo "   'python slurm/cost.py'; set ALLOW_EXPENSIVE_GPU=1 to override." >&2
            exit 2
        fi ;;
esac

echo "== $(date -Is) $(hostname) part=${SLURM_JOB_PARTITION:-none} job=${SLURM_JOB_ID:-none}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1
NJOBS=$(printf '%s\n' "$JOBS" | grep -c 'SYS=' || true)
echo "== $NJOBS runs sharing one GPU; this job bills ONE device, not $NJOBS."

# CUDA MPS lets concurrent contexts interleave kernels rather than time-slice, which is
# where most of the packing gain comes from for systems too small to fill the card.
export CUDA_MPS_PIPE_DIRECTORY="${TMPDIR:-/tmp}/mps-${SLURM_JOB_ID:-$$}"
export CUDA_MPS_LOG_DIRECTORY="$CUDA_MPS_PIPE_DIRECTORY/log"
mkdir -p "$CUDA_MPS_LOG_DIRECTORY"
MPS=no
if command -v nvidia-cuda-mps-control >/dev/null 2>&1; then
    if nvidia-cuda-mps-control -d >/dev/null 2>&1; then MPS=yes; fi
fi
echo "== CUDA MPS: $MPS"
stop_mps() { [[ "$MPS" == yes ]] && (echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true); }

pids=(); labels=()
while IFS= read -r spec; do
    [[ "$spec" =~ SYS= ]] || continue
    SYS=""; REP=""; X0=""; K=""; RAMP_PS=""; CLAMP_PN=""; RNS="$NS_DEFAULT"
    for kv in $spec; do
        case "$kv" in
            SYS=*) SYS="${kv#SYS=}" ;;
            REP=*) REP="${kv#REP=}" ;;
            NS=*)  RNS="${kv#NS=}" ;;
            X0=*)  X0="${kv#X0=}" ;;
            K=*)   K="${kv#K=}" ;;
            RAMP_PS=*) RAMP_PS="${kv#RAMP_PS=}" ;;
            CLAMP_PN=*) CLAMP_PN="${kv#CLAMP_PN=}" ;;
            *) echo "== unknown key in spec: $kv" >&2; exit 2 ;;
        esac
    done
    [[ -n "$SYS" && -n "$REP" ]] || { echo "== spec needs SYS and REP: $spec" >&2; exit 2; }
    SYSDIR="build/systems/$SYS"
    [[ -f "$SYSDIR/system.prmtop" ]] || { echo "== no topology for $SYS" >&2; exit 2; }

    # Identical TAG construction to md_chain.sbatch.  Several windows and clamp forces
    # run off ONE topology, so the condition has to be in the directory name or they
    # overwrite each other's trajectory.
    TAG="$SYS${CLAMP_PN:+_${CLAMP_PN}pN}${X0:+_x${X0}A}${K:+_k${K}}"
    case "$TAG" in
        *_x_*|*_xA_*|*_k) echo "== malformed TAG '$TAG' from spec: $spec" >&2; exit 2 ;;
    esac
    [[ "$REP" =~ ^[0-9]+$ ]] || { echo "== REP not a number in: $spec" >&2; exit 2; }
    RUNDIR="data/runs/$TAG/rep$REP"
    mkdir -p "$RUNDIR"
    labels+=("$TAG/rep$REP")
    echo "== queued $TAG/rep$REP  ns=$RNS"

    (
        set -euo pipefail
        CLAMP_ARGS=()
        if [[ -n "$X0" ]]; then
            CA="$(cat "$SYSDIR/clamp_atoms.txt")"
            CLAMP_ARGS=(--restraint-x0 "$X0" --clamp-atoms "$CA")
            [[ -n "$K" ]] && CLAMP_ARGS+=(--restraint-k "$K")
            [[ -n "$RAMP_PS" ]] && CLAMP_ARGS+=(--restraint-ramp-ps "$RAMP_PS")
        elif [[ -n "$CLAMP_PN" ]]; then
            CA="$(cat "$SYSDIR/clamp_atoms.txt")"
            CLAMP_ARGS=(--clamp-pN "$CLAMP_PN" --clamp-atoms "$CA")
        fi
        COMMON=(--prmtop "$SYSDIR/system.prmtop" --inpcrd "$SYSDIR/system.inpcrd"
                --platform CUDA)
        if [[ ! -f "$RUNDIR/01_minimise/minimised.xml" ]]; then
            "$DNASPRING_PY" runs/openmm_engine.py minimise "${COMMON[@]}" \
                --outdir "$RUNDIR/01_minimise"
        fi
        if [[ ! -f "$RUNDIR/02_equilibrate/equilibrated.xml" ]]; then
            "$DNASPRING_PY" runs/openmm_engine.py equilibrate "${COMMON[@]}" \
                --outdir "$RUNDIR/02_equilibrate" --seed "$((1000 + REP))"
        fi
        "$DNASPRING_PY" runs/openmm_engine.py produce "${COMMON[@]}" \
            --outdir "$RUNDIR/03_production" --ns "$RNS" --seed "$((2000 + REP))" \
            "${CLAMP_ARGS[@]}"
    ) > "$RUNDIR/pack.log" 2>&1 &
    pids+=($!)
    sleep 8          # stagger context creation: simultaneous CUDA init can exhaust GPU RAM
done <<< "$JOBS"

# A note on utilisation: on a MIG-enabled device `utilization.gpu` returns [N/A] --
# nvidia-smi cannot report SM utilisation per MIG instance -- so on gpu_test this column
# is empty and the throughput comparison below is the only measure available.  On
# gpu_requeue's whole cards it is reported normally.  Per-process GPU memory comes from
# --query-compute-apps, which works in both cases and is what shows how little of the
# card one run needs.
USAGE="slurm/logs/pack_usage_${SLURM_JOB_ID:-local}.csv"
echo "iso,gpu_util_pct,gpu_mem_used_MiB,gpu_mem_total_MiB,n_procs,mem_per_proc_MiB" > "$USAGE"
( while sleep 60; do
    g=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total \
        --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    apps=$(nvidia-smi --query-compute-apps=used_memory --format=csv,noheader,nounits \
        2>/dev/null | tr -d ' ' | grep -c . || true)
    per=$(nvidia-smi --query-compute-apps=used_memory --format=csv,noheader,nounits \
        2>/dev/null | tr -d ' ' | awk '{s+=$1; n++} END {if(n) printf "%.0f", s/n}')
    [[ -n "$g" ]] && echo "$(date -Is),$g,${apps:-0},${per:-0}" >> "$USAGE"
  done ) & SAMPLER=$!
trap 'kill "$SAMPLER" 2>/dev/null || true; stop_mps' EXIT

fail=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        echo "== ${labels[$i]} finished"
    else
        echo "== ${labels[$i]} FAILED -- see data/runs/${labels[$i]}/pack.log" >&2
        fail=1
    fi
done
kill "$SAMPLER" 2>/dev/null || true

echo
echo "== what this job actually got out of its one GPU =="
"$DNASPRING_PY" - "$USAGE" "${labels[@]}" <<'PYEOF'
import os, sys
usage, labels = sys.argv[1], sys.argv[2:]

def speed(tag):
    p = f"data/runs/{tag}/03_production/production.log"
    try:
        last = None
        for line in open(p):
            if line[:1].isdigit():
                last = line
        return float(last.split()[7]) if last else None
    except (OSError, IndexError, ValueError, TypeError):
        return None

sp = [(t, speed(t)) for t in labels]
got = [s for _, s in sp if s]
for t, s in sp:
    print(f"  {t:34s} {s:7.0f} ns/day" if s else f"  {t:34s}       — no speed logged")
if got:
    print(f"\n  aggregate {sum(got):.0f} ns/day out of ONE billed GPU "
          f"({len(got)} runs, mean {sum(got)/len(got):.0f} each)")
if os.path.exists(usage):
    rows = [l.split(",") for l in open(usage).read().splitlines()[1:] if l.strip()]
    if rows:
        u = [float(r[1]) for r in rows if r[1] not in ("[N/A]", "N/A", "")]
        mem = max(float(r[2]) for r in rows) / 1024
        tot = max(float(r[3]) for r in rows) / 1024
        print(f"  GPU memory {mem:.2f} of {tot:.1f} GB used by the whole pack")
        per = [float(r[5]) for r in rows if len(r) > 5 and r[5] not in ("", "0")]
        np = [int(r[4]) for r in rows if len(r) > 4 and r[4] not in ("", "0")]
        if per and np:
            print(f"  {max(np)} processes on the card, ~{sum(per)/len(per):.0f} MiB each "
                  f"-- one run alone would use {100*(sum(per)/len(per))/(tot*1024):.1f}% "
                  f"of its memory")
        if u:
            print(f"  GPU {sum(u)/len(u):.0f}% mean / {max(u):.0f}% peak utilisation "
                  f"({len(u)} samples)")
            if sum(u) / len(u) < 80:
                print("  Under 80% mean: the card still has headroom, so a larger pack "
                      "would divide the same GPU bill further.")
        else:
            print("  Utilisation not reportable on this device (MIG instances return "
                  "[N/A]); judge packing by the aggregate throughput above instead.")
PYEOF
echo "== $(date -Is) done (fail=$fail)"
exit "$fail"
