#!/bin/bash -l
# Runs one CryoSPARC stage driver (scripts/<stage>.py) off the desktop, so it outlives the OOD session.
# Drivers only poll CryoSPARC (1 core, <1 GB); they re-attach to running jobs via cs_state.json when restarted.
#   sbatch -J sfgfp-drv-11 slurm/driver.sbatch scripts/11_preprocess.py
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 0-12:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
cd /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening
source env/activate_sfgfp.sh >/dev/null 2>&1
PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY -u "$@" 2>&1 | grep --line-buffered -v "^Success" | tee -a "$SFGFP_LOGS/sfgfp-$(basename "$1" .py)-driver.log"
