#!/bin/bash
# Environment for the PP2A Aalpha DNA-spring MD.  `source` this file.
#
# Deliberately does NOT source ../sfGFP/env/activate.sh: that exports MAMBA paths on
# /n/lab_storage, which is full.  Everything here is home (code), holylabs (software) or
# netscratch (data), per the hekstra-cluster rules.
export PP2A_MD_REPO="/n/home05/ianpoon/DNA-springs-project-2026-27/MD simulations/pp2a"
export SFGFP_REPO="/n/home05/ianpoon/DNA-springs-project-2026-27/MD simulations/sfGFP"
export PP2A_MD_SCRATCH=/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a-md
export PP2A_DATA=/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025     # the AF model + apo run
export PP2A_LOGS=/n/netscratch/hekstra_lab/Lab/ian_poon/logs
export DNASPRING_ENV=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring
export DNASPRING_PY="$DNASPRING_ENV/bin/python"
export AMBERHOME="$DNASPRING_ENV"
export PATH="$DNASPRING_ENV/bin:$PATH"     # tleap, antechamber, pymol, ffmpeg, cpptraj
export MPLBACKEND=Agg
