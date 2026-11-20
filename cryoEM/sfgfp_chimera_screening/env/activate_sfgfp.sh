# Source me: paths for the sfGFP chimera screening cryo-EM project.
export SFGFP_REPO=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/sfgfp_chimera_screening
export SFGFP_DATA=/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112
export SFGFP_RAW=$SFGFP_DATA/raw
export SFGFP_LOGS=/n/netscratch/hekstra_lab/Lab/ian_poon/logs
export SFGFP_RSYNC_URL=rsync://em00000013@10.245.3.142/em00000013
export SFGFP_RSYNC_PASS=$HOME/.config/rsync/em00000013.pass
export CS_PROJECT_TITLE=sfgfp-chimera-screening
export CS_PROJECT_PARENT=/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects
# Reuse the PP2A reproduction's tool paths (dnaspring python, Phenix, RELION, GROMACS, ChimeraX, EMReady, CS python).
source /n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction/env/activate_pp2a.sh
export TOPAZ_CS=/n/holylabs/hekstra_lab/Everyone/ianpoon/software/bin/topaz-cs
