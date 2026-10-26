#!/bin/bash -l
# Plan section 4 - MD orchestrator.  gpu_test's QOS allows 2 submitted jobs per user and the Open OnDemand
# desktop holds one, so nothing can be queued behind a running GPU job (a --dependency submission is rejected
# with QOSMaxSubmitJobPerUserLimit).  This 1-core `test` job polls every 2 min and, whenever no pp2a-md job of
# ours is in the gpu_test queue, submits the next step:
#     Aa bench (1 ns) -> Aa 100 ns chunks until progress.json says 100 ns -> Ab bench -> Ab chunks.
# Each chunk is slurm/pp2a_md_chain.sbatch with --time sized from the measured ns/day (x1.15, <= 12 h).
# `test` caps jobs at 12 h, so the orchestrator resubmits itself 15 min before its own limit.  Free on both partitions.
#SBATCH -J pp2a-md-orchestrator
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 12:00:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -uo pipefail
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/env/activate_pp2a.sh
cd "$PP2A_REPO"
START=$(date +%s); LIMIT=$((11*3600+45*60)); NS_TOTAL=${NS_TOTAL:-100}
progress() { "$DNASPRING_PY" -c "import json,sys;d=json.load(open('$1'));print(d.get('ns_done',0),d.get('ns_per_day_this_job',0))" 2>/dev/null || echo "0 0"; }
hours_for() { "$DNASPRING_PY" -c "rem=float('$1'); sp=float('$2') or 20.0; h=min(12.0,max(0.5,1.15*rem/sp*24)); print(f'{int(h)}:{int((h-int(h))*60):02d}:00')"; }
while true; do
  if [ $(( $(date +%s) - START )) -gt $LIMIT ]; then
    echo "$(date -Is) resubmitting orchestrator"; sbatch "$PP2A_REPO/slurm/pp2a_md_orchestrator.sbatch"; exit 0
  fi
  if squeue -u "$USER" -h -p gpu_test -o "%j" | grep -q "^pp2a-md"; then sleep 120; continue; fi
  next=""
  for s in Aa Ab; do
    bench=$PP2A_DATA/md/runs/${s}_bench1ns/progress.json; run=$PP2A_DATA/md/runs/${s}_rep1/progress.json
    if [ ! -f "$bench" ] || [ "$(progress "$bench" | awk '{print ($1>=1)}')" != "1" ]; then next="bench $s"; break; fi
    read -r done speed <<<"$(progress "$run")"
    if "$DNASPRING_PY" -c "import sys;sys.exit(0 if float('$done')<$NS_TOTAL-1e-6 else 1)"; then
      [ "$speed" = "0" ] && read -r _ speed <<<"$(progress "$bench")"
      next="chain $s $("$DNASPRING_PY" -c "print($NS_TOTAL-float('$done'))") $speed"; break
    fi
  done
  if [ -z "$next" ]; then echo "$(date -Is) all MD done"; exit 0; fi
  set -- $next
  case $1 in
    bench) echo "$(date -Is) submit benchmark $2"; SYS=$2 NS=1 sbatch -J pp2a-md-bench-$2 -t 2:00:00 "$PP2A_REPO/slurm/pp2a_md_chain.sbatch" ;;
    chain) h=$(hours_for "$3" "$4"); echo "$(date -Is) submit chain $2: $3 ns left at $4 ns/day -> --time=$h"
           SYS=$2 NS=$NS_TOTAL sbatch -J pp2a-md-$2 -t "$h" "$PP2A_REPO/slurm/pp2a_md_chain.sbatch" ;;
  esac
  sleep 180
done
