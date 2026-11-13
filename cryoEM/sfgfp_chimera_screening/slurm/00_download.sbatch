#!/bin/bash -l
# Holds an ssh session to a login node open while scripts/download_login.sh runs THERE.
# Why: compute nodes have no route to the facility server (10.245.3.142:873, probe 60306535), and
# login nodes kill a user's processes (nohup/setsid/tmux included) when the session closes.
# rsync is resumable, so if this job ends early just resubmit it.
#SBATCH -J sfgfp-download
#SBATCH -p test
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 0-04:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
S=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening/scripts/download_login.sh
ssh -o BatchMode=yes -o ServerAliveInterval=60 -o ServerAliveCountMax=10 holylogin05.rc.fas.harvard.edu "bash $S"
echo "ssh exit $? $(date -Is)"
