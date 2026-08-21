#!/bin/bash
# Re-run every analysis over whatever has finished, redraw the figures and animations,
# and rebuild the report. Safe to run repeatedly while production is still going.
set -euo pipefail
cd "$(dirname "$0")"
source ./env/activate.sh
export PATH="$DNASPRING_ENV/bin:$PATH"
P="$DNASPRING_PY"

STRIDE="${STRIDE:-5}"
$P analysis/collect.py --stride "$STRIDE"
$P analysis/figures.py

# Animate the longest chimera replicate available, plus a control for comparison.
for sys in S3_spring27 S0_wt; do
  traj=$(ls -S data/runs/$sys/rep*/03_production/traj.dcd 2>/dev/null | head -1) || true
  [ -n "${traj:-}" ] || continue
  $P analysis/animate.py --traj "$traj" --top "build/systems/$sys/system.prmtop" \
     --name "$sys" --frames "${FRAMES:-80}"
done

$P report/build_report.py
echo
echo "report: $(pwd)/report/sfgfp_dna_spring_report.html"
