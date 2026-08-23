#!/bin/bash
# Measure f(x) for the bare springs directly from the force field (Doeke's point 6).
#
# The report currently quotes Zocchi's analytic f(x) evaluated at the end-to-end
# distance the chimera settles at.  That is a model, and he is right to ask which it is.
# These runs measure the same quantity with no model in the way: hold the two ends of a
# solvated duplex -- no protein -- at a fixed separation x with a harmonic restraint,
# and the mean restraint force is what the spring pushes back with, f = k(<r> - x0).
#
# The low-x windows do double duty for point 8.  They compress the duplex far past
# where it sits in the chimera, so if a kink is going to nucleate anywhere in this force
# field, it is there.  Windows are spaced ~8 A apart, which is well outside the
# sqrt(kT/k) ~ 1.1 A spread of a single window.
set -euo pipefail
cd "$(dirname "$0")/.."

NS=${NS:-30}; REPS=${REPS:-2}; DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

# x0 ladders chosen from each construct's own contour: from strongly buckled up to
# nearly straight.  27 bp spans 91 A end to end when straight, 40 bp spans 134 A.
declare -A WINDOWS=(
  [D_27bp]="50 58 66 74 82 88"
  [D_40bp]="50 66 82 98 114 128"
  [D_40bp_nick]="50 66 82 98 114 128"
)

n=0
for sys in D_27bp D_40bp D_40bp_nick; do
  for x0 in ${WINDOWS[$sys]}; do
    name="${sys}_x${x0}"
    if [[ $DRY -eq 1 ]]; then
      echo "sbatch --array=1-$REPS -J $name SYS=$sys RESTRAINT_X0=$x0 NS=$NS"
    else
      sbatch --array="1-$REPS" -J "$name" \
        --export=ALL,SYS="$sys",RESTRAINT_X0="$x0",NS="$NS" slurm/md_chain.sbatch
    fi
    n=$((n+1))
  done
done
echo "$n windows x $REPS replicates x $NS ns"
