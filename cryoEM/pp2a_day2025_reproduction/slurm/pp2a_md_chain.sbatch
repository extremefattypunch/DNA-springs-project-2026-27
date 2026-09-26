#!/bin/bash -l
# Plan section 4 - one PP2A A-subunit MD run (Aa or Ab), restart-safe across 12 h gpu_test jobs.
#
#  * gpu_test only (Gres/gpu billing weight 0; the `gpu` partition costs 836.5 per GPU-hour).  MIG A100
#    3g.20gb slices: a 165k-particle OPC system uses ~3 GB and one host core, so -c 1 / --mem=8G.
#  * The engine resumes from state.chk + progress.json.  slurm/pp2a_md_orchestrator.sbatch submits the
#    next chunk when the GPU slot is free, with --time sized from the measured ns/day (x1.15, <= 12 h) so the
#    last link does not trip the Job Defense Shield's "<50 % of the walltime used" rule.
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

# Progress summary; the next chunk is submitted by slurm/pp2a_md_orchestrator.sbatch (gpu_test's 2-job QOS
# rejects a --dependency resubmission while this job is still running).
if [ -f "$RUN/progress.json" ]; then
  "$DNASPRING_PY" -c "import json;d=json.load(open('$RUN/progress.json'));print('== progress:',d.get('ns_done'),'/',d.get('ns_total'),'ns at',d.get('ns_per_day_this_job'),'ns/day')"
fi
exit $rc
