#!/bin/bash -l
#SBATCH -J census
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -e
REPO=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening
export OMP_NUM_THREADS=4
cd $REPO
/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python -I scripts/51_census_bootstrap.py
