#!/bin/bash
# Wait for the queue to drain, then run every analysis over the finished trajectories.
# Deliberately stops short of publishing: the numbers get read before they go out.
set -euo pipefail
cd "$(dirname "$0")"
source ./env/activate.sh
P="$DNASPRING_PY"

while squeue -u "$USER" -h -o "%j" | grep -qvE '^\.fasrcood'; do
  sleep 300
done
echo "== queue drained at $(date -Is)"

$P analysis/collect.py --stride 5 --skip-ns 5 --min-ns 20 > analysis/_collect.log 2>&1
$P analysis/findings.py                                    >> analysis/_collect.log 2>&1
$P analysis/force_extension.py --skip-ns 10 --min-ns 20    >> analysis/_collect.log 2>&1
$P analysis/figures.py                                     >> analysis/_collect.log 2>&1
echo "== analysis complete at $(date -Is)"
