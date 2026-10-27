# Source me:  source cryoEM/manim_cryoem_math/env/activate_manim.sh
# Manim CE env + TinyTeX live on holylabs (software tier); caches on netscratch.
HOLY=/n/holylabs/hekstra_lab/Everyone/ianpoon
export MAMBA_ROOT_PREFIX=$HOLY/micromamba
eval "$($HOLY/bin/micromamba shell hook --shell bash)"
micromamba activate $HOLY/envs/manim
export PATH=$HOLY/software/tinytex/bin/x86_64-linux:$PATH
export MANIM_CRYOEM_SCRATCH=/n/netscratch/hekstra_lab/Lab/ian_poon/manim_cryoem_math
export MANIM_CRYOEM_CACHE=$MANIM_CRYOEM_SCRATCH/cache
mkdir -p "$MANIM_CRYOEM_SCRATCH/media" "$MANIM_CRYOEM_CACHE"
# keep numpy single-threaded so parallel renders don't oversubscribe the cores
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
