#!/bin/bash
# Plan section 4 - submit the Tier 1 MD: for each system a 1 ns benchmark, then the 100 ns chain.
#
# gpu_test's QOS allows 2 jobs per user and Ian's Open OnDemand desktop holds one of them, so only ONE
# GPU job runs at a time: Aa and Ab are queued strictly one after another (afterany dependencies), and
# each 100 ns chain resubmits itself every 12 h (slurm/pp2a_md_chain.sbatch).  Expected wall time at
# 30-60 ns/day per MIG slice: 2-3.5 days per system -> 4-7 days for both.  Zero billing.
# Usage:  bash scripts/04_submit_md.sh            (waits for any running gpu_test job of ours first)
set -euo pipefail
source "$(dirname "$0")/../env/activate_pp2a.sh"
cd "$PP2A_REPO"
for s in Aa Ab; do [ -f "$PP2A_DATA/md/systems/$s/system.prmtop" ] || { echo "run 02_build_systems.py first"; exit 1; }; done
prev=""
for s in Aa Ab; do
  dep=${prev:+--dependency=afterany:$prev}
  b=$(SYS=$s NS=1 sbatch --parsable $dep -t 1:30:00 -J pp2a-md-bench-$s slurm/pp2a_md_chain.sbatch)
  echo "benchmark $s: job $b ${dep:+(after $prev)}"
  c=$(SYS=$s NS=100 sbatch --parsable --dependency=afterany:$b -J pp2a-md-$s slurm/pp2a_md_chain.sbatch)
  echo "100 ns chain $s: job $c (after $b; resubmits itself until 100 ns)"
  prev=$c
done
squeue -u "$USER" -o "%i %P %j %T %M %E" | grep -E "pp2a-md|JOBID"
