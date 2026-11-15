#!/bin/bash
# Fig S5 ChimeraX renders / fits, one file per job (renders run one at a time).
#   sbatch slurm/52_figS5_chimerax.sbatch scripts/51_figS5_render_W3_duplex.cxc
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 1:00:00
#SBATCH -J figS5_cx
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/slurm_%x_%j.log
set -euo pipefail
f=${1:?cxc or py file}
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening/env/activate_sfgfp.sh
cd /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening
echo "host $(hostname)  file $f  $(date)"
$CHIMERAX --nogui --offscreen --silent "$f"
echo "done $(date)"
