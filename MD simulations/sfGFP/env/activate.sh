#!/bin/bash
# Source this instead of relying on ~/.bashrc: the mamba block there still points at
# /n/hekstra_lab/... which was retired when the share moved to /n/lab_storage/hekstra_lab.
export MAMBA_EXE=/n/lab_storage/hekstra_lab/people/ian_poon/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/n/lab_storage/hekstra_lab/people/ian_poon/micromamba
export DNASPRING_ENV=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring
export DNASPRING_QM_ENV=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring-qm
export DNASPRING_PY="$DNASPRING_ENV/bin/python"
export DNASPRING_SCRATCH=/n/netscratch/hekstra_lab/Lab/ian_poon/sfGFP-md
export MPLBACKEND=Agg
