#!/bin/bash -l
# Plan section 4 - one PP2A A-subunit MD run (Aa or Ab), restart-safe across 12 h gpu_test jobs.
#
#  * gpu_test only (Gres/gpu billing weight 0; the `gpu` partition costs 836.5 per GPU-hour).  MIG A100
#    3g.20gb slices: a 165k-particle OPC system uses ~3 GB and one host core, so -c 1 / --mem=8G.
#  * The engine resumes from state.chk + progress.json.  At exit, if production is short of NS, this
#    script resubmits itself with --dependency=afterany and a --time sized from the measured ns/day
#    (x1.15, capped at 12 h) so the last link does not trip the Job Defense Shield's "<50 % of the
#    walltime used" rule.  A 1 ns benchmark (NS=1) is the first submission of every system.
# Usage:  SYS=Aa NS=100 sbatch slurm/pp2a_md_chain.sbatch      (SYS in {Aa, Ab}; SEED default 2001)
#SBATCH -J pp2a-md
#SBATCH -p gpu_test
#SBATCH --gres=gpu:1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH --constraint=a100-mig
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/env/activate_pp2a.sh
: "${SYS:?SYS=Aa or Ab}"; NS="${NS:-100}"; SEED="${SEED:-2001}"; HMR="${HMR:-}"
SYSDIR=$PP2A_DATA/md/systems/$SYS
RUN=$PP2A_DATA/md/runs/${SYS}_rep1${HMR:+_hmr}
[ "$NS" = "1" ] && RUN=$PP2A_DATA/md/runs/${SYS}_bench1ns
mkdir -p "$RUN"
case "${SLURM_JOB_PARTITION:-}" in gpu_test) ;; *) echo "refusing to run on ${SLURM_JOB_PARTITION} (billed)"; exit 2 ;; esac
echo "== $(date -Is) $(hostname) job=$SLURM_JOB_ID SYS=$SYS NS=$NS RUN=$RUN"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1
export TMPDIR=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT

# GPU-utilisation sampler for the report (MIG shows N/A for util; memory still tells us it is in use)
( while true; do echo "$(date +%s),$(nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits | head -1)" >> "$RUN/gpu_usage.csv"; sleep 120; done ) &
SAMPLER=$!

set +e
"$DNASPRING_PY" "$PP2A_REPO/scripts/03_openmm_pp2a_engine.py" --sysdir "$SYSDIR" --outdir "$RUN" --ns "$NS" --seed "$SEED" ${HMR:+--hmr}
rc=$?
set -e
kill $SAMPLER 2>/dev/null || true

# self-resubmit if production is not finished (engine exits 0 after each checkpointed chunk only when done)
if [ -f "$RUN/progress.json" ]; then
  read -r done total speed <<<"$("$DNASPRING_PY" -c "import json;d=json.load(open('$RUN/progress.json'));print(d.get('ns_done',0),d.get('ns_total',$NS),d.get('ns_per_day_this_job',0))")"
  echo "== progress: $done / $total ns at $speed ns/day (rc=$rc)"
  if "$DNASPRING_PY" -c "import sys; sys.exit(0 if float('$done') < float('$total') - 1e-6 else 1)"; then
    hours=$("$DNASPRING_PY" -c "
rem=float('$total')-float('$done'); sp=float('$speed') or 20.0
h=min(12.0, max(0.5, 1.15*rem/sp*24)); print(f'{int(h)}:{int((h-int(h))*60):02d}:00')")
    echo "== resubmitting for the remaining $("$DNASPRING_PY" -c "print(round(float('$total')-float('$done'),2))") ns with --time=$hours"
    SYS=$SYS NS=$NS SEED=$SEED HMR=$HMR sbatch --dependency=afterany:$SLURM_JOB_ID -t "$hours" \
        "$PP2A_REPO/slurm/pp2a_md_chain.sbatch"
  fi
fi
exit $rc
