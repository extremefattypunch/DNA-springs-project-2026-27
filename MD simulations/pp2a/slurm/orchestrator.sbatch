#!/bin/bash -l
# Walks slurm/jobs.txt and keeps the runs moving without anyone watching.
#
# Why it exists: gpu_test's QOS allows 2 jobs per user and the Open OnDemand desktop
# holds one, so nothing can be queued behind a running gpu_test job and --dependency is
# rejected.  This 1-core job on the free `test` partition polls every 2 min and submits
# the next unfinished line of jobs.txt whenever its partition has room:
#   * gpu_test lines: one at a time (the desktop takes the other slot);
#   * gpu_requeue / gpu lines: up to MAX_PAID at once, and only after the free 1 ns
#     benchmark of the same system exists, so --time is sized from a measured ns/day
#     (bench speed x PAID_FACTOR for a full card vs the MIG slice, then the run's own
#     production.log once it has one) x 1.15.  Over-asking trips the Job Defense
#     Shield's "elapsed < 50% of limit"; under-asking just costs one more chunk, since
#     the engine resumes from state.chk.
# A line is finished when its progress.json reports ns_done >= NS.  `test` caps jobs at
# 12 h, so the orchestrator resubmits itself 15 min before its own limit.
#
#SBATCH -J pp2amd-orch
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 12:00:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -uo pipefail
source "/n/home05/ianpoon/DNA-springs-project-2026-27/MD simulations/pp2a/env/activate.sh"
cd "$PP2A_MD_REPO"
START=$(date +%s); LIMIT=$((11*3600+45*60))
MAX_PAID=${MAX_PAID:-3}; PAID_FACTOR=${PAID_FACTOR:-2.0}; JOBS=${JOBS:-slurm/jobs.txt}

py() { "$DNASPRING_PY" -c "$1"; }
ns_done() { py "import json;print(json.load(open('$1')).get('ns_done',0))" 2>/dev/null || echo 0; }
log_speed() {  # ns/day from the last line of a production.log (column 8), or empty
  local s; s=$(awk -F'\t' 'END{print $8}' "$1" 2>/dev/null)
  [[ "$s" =~ ^[0-9.]+$ ]] && [[ "$s" != "0" ]] && echo "$s"
}
hours_for() { py "
rem=float('$1'); sp=float('$2'); cap=float('$3')
h=min(cap,max(0.5,1.15*rem/sp*24+0.4))
print(f'{int(h)}:{int((h-int(h))*60):02d}:00')"; }
is_done() { py "import sys;sys.exit(0 if float('$1')>=float('$2')-1e-6 else 1)"; }

while true; do
  if [ $(( $(date +%s) - START )) -gt $LIMIT ]; then
    echo "$(date -Is) resubmitting orchestrator"
    until jid=$(sbatch --parsable "$PP2A_MD_REPO/slurm/orchestrator.sbatch"); do
      echo "$(date -Is) orchestrator resubmit FAILED; retrying in 60 s"; sleep 60
    done
    echo "   orchestrator job $jid"; exit 0
  fi
  # One squeue per pass, and the pass is skipped if it fails: a slurmctld timeout that
  # returned nothing would otherwise read as "nothing queued" and resubmit every running job.
  if ! jobs_now=$(squeue -u "$USER" -h -o "%j %P"); then
    echo "$(date -Is) squeue failed; skipping this pass"; sleep 120; continue
  fi
  queued_test=$(awk '$1 ~ /^pp2amd/ && $2 ~ /(^|,)gpu_test(,|$)/' <<<"$jobs_now" | wc -l)
  queued_paid=$(awk '$1 ~ /^pp2amd/ && $2 ~ /(^|,)gpu(_requeue)?(,|$)/' <<<"$jobs_now" | wc -l)
  pending=0; submitted=0; need_publish=0; deferred=0
  PUBLISHED="$PP2A_MD_SCRATCH/published.txt"; touch "$PUBLISHED"
  while read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    SYS=; NS=; REP=1; PART=gpu_test; X0=; K=; RAMP_PS=; HOURS=; SUFFIX=
    for kv in $line; do eval "$kv"; done
    TAG="$SYS${X0:+_x${X0}A}${K:+_k${K}}${SUFFIX}"
    RUN="$PP2A_MD_SCRATCH/runs/$TAG/rep$REP"
    done_ns=$(ns_done "$RUN/03_production/progress.json")
    if is_done "$done_ns" "$NS"; then
      # a run that just finished: analyse everything and push the page (not for benchmarks)
      if [ -z "$SUFFIX" ] && ! grep -qx "$TAG" "$PUBLISHED"; then echo "$TAG" >> "$PUBLISHED"; need_publish=1; fi
      continue
    fi
    pending=$((pending+1))
    name="pp2amd-$TAG"
    awk '{print $1}' <<<"$jobs_now" | grep -qx "$name" && continue
    bench_log="$PP2A_MD_SCRATCH/runs/${SYS}_bench/rep1/03_production/production.log"
    bench_prog="$PP2A_MD_SCRATCH/runs/${SYS}_bench/rep1/03_production/progress.json"
    case $PART in
      gpu_test) [ "$queued_test" -ge 1 ] && continue; cons="--constraint=a100-mig"; cap=12 ;;
      gpu_requeue|gpu)
        [ "$queued_paid" -ge "$MAX_PAID" ] && continue
        if ! is_done "$(ns_done "$bench_prog")" 1; then continue; fi   # wait for the free benchmark
        cons="--constraint=h200|h100|a100"; cap=72 ;;   # no MIG slices at full price
      *) echo "bad PART $PART"; continue ;;
    esac
    speed=$(log_speed "$RUN/03_production/production.log")
    if [ -z "$speed" ]; then
      b=$(log_speed "$bench_log")
      case $PART in gpu_test) speed=${b:-40} ;; *) speed=$(py "print(${b:-40}*$PAID_FACTOR)") ;; esac
    fi
    rem=$(py "print(float('$NS')-float('$done_ns'))")
    h=${HOURS:+$HOURS:00:00}; h=${h:-$(hours_for "$rem" "$speed" "$cap")}
    env_args="SYS=$SYS,NS=$NS,REP=$REP,RUN_SUFFIX=$SUFFIX"
    [ -n "$X0" ] && env_args="$env_args,RESTRAINT_X0=$X0,RESTRAINT_K=$K,RESTRAINT_RAMP_PS=$RAMP_PS"
    [ "$PART" = gpu ] && env_args="$env_args,ALLOW_EXPENSIVE_GPU=1"
    echo "$(date -Is) submit $name on $PART: $rem ns left at $speed ns/day -> --time=$h"
    jid=$(sbatch --parsable -J "$name" -p "$PART" $cons -t "$h" --export=ALL,$env_args "$PP2A_MD_REPO/slurm/md_chain.sbatch") \
      && echo "   job $jid" || echo "   sbatch FAILED"
    submitted=$((submitted+1))
    case $PART in gpu_test) queued_test=$((queued_test+1));; *) queued_paid=$((queued_paid+1));; esac
  done < "$JOBS"
  if [ "$need_publish" -eq 1 ]; then
    if awk '{print $1}' <<<"$jobs_now" | grep -q "^pp2amd-finalize"; then
      echo "$(date -Is) a finalize job is already queued; deferring the publish to the next pass"
      sed -i '$ d' "$PUBLISHED"; deferred=1
    else
      echo "$(date -Is) submit finalize + publish"
      PUBLISH=1 sbatch "$PP2A_MD_REPO/slurm/finalize.sbatch" || { sed -i '$ d' "$PUBLISHED"; deferred=1; }
    fi
  fi
  if [ "$pending" -eq 0 ] && [ "$deferred" -eq 0 ]; then echo "$(date -Is) all runs in $JOBS complete"; exit 0; fi
  [ "$submitted" -gt 0 ] && sleep 180 || sleep 120
done
