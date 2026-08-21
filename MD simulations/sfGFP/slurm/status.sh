#!/bin/bash
# One-screen progress report: queue state plus how far each production run has got.
set -euo pipefail
cd "$(dirname "$0")/.."
echo "=== queue ==="
squeue -u "$USER" -o "%.10i %.14j %.10P %.2t %.11M %.11l %R" || true
echo
echo "=== production progress ==="
printf "%-22s %-6s %10s %10s %8s\n" system rep ns_done ns_total pct
shopt -s nullglob
for p in data/runs/*/rep*/03_production/progress.json; do
    sys=$(basename "$(dirname "$(dirname "$(dirname "$p")")")")
    rep=$(basename "$(dirname "$(dirname "$p")")")
    read -r d t < <(python3 -c "
import json,sys
j=json.load(open(sys.argv[1]))
print(j['ns_done'], round(j['steps_total']*j['dt_fs']/1e6,3))" "$p")
    pct=$(python3 -c "print(f'{100*$d/$t:.1f}' if $t else '0')")
    printf "%-22s %-6s %10s %10s %7s%%\n" "$sys" "$rep" "$d" "$t" "$pct"
done
