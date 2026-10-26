#!/bin/bash
# Queue + per-run progress for the PP2A spring runs.
source "$(dirname "$0")/../env/activate.sh"
echo "== queue"
squeue -u "$USER" -o "%.10i %.28j %.12P %.9T %.11M %.10l %R" | grep -E "JOBID|pp2amd" || true
echo "== runs in $PP2A_MD_SCRATCH/runs"
printf "  %-28s %8s %8s %6s %10s\n" run ns_done ns_total pct ns_per_day
for p in "$PP2A_MD_SCRATCH"/runs/*/rep*/03_production/progress.json; do
  [ -f "$p" ] || continue
  run=$(basename "$(dirname "$(dirname "$(dirname "$p")")")")/$(basename "$(dirname "$(dirname "$p")")")
  speed=$(awk -F'\t' 'END{print $8}' "$(dirname "$p")/production.log" 2>/dev/null)
  "$DNASPRING_PY" -c "
import json; d=json.load(open('$p')); tot=d['steps_total']*d['dt_fs']/1e6
print(f\"  {'$run':<28} {d['ns_done']:8.1f} {tot:8.1f} {100*d['ns_done']/tot:5.0f}% {'$speed':>10}\")"
done
echo "== fairshare"; sshare -A hekstra_lab -u "$USER" 2>/dev/null | tail -2
