#!/bin/bash -l
# Final 1080p30 renders of the six clips: one array task per clip on the free `test`
# partition (CPU-only; Manim's Cairo renderer is single-threaded, +1 core for x264).
#   sbatch cryoEM/manim_cryoem_math/slurm/render_final.sbatch          # all six
#   sbatch --array=2 cryoEM/manim_cryoem_math/slurm/render_final.sbatch  # just clip 03
# -t is set from measured -ql timings before submitting (Job Defense Shield: elapsed >= 50% of limit)
#SBATCH -J manim-final
#SBATCH -p test
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 0-01:30
#SBATCH --array=0-5
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%A_%a.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%A_%a.err
set -euo pipefail
source /etc/profile.d/lmod.sh
ROOT=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/manim_cryoem_math
cd "$ROOT"
set +u; source env/activate_manim.sh; set -u
TMP=/scratch/${USER}_${SLURM_JOB_ID}
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
export TMPDIR="$TMP"

read -r MOD CLS OUT < <(python main.py --list | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
echo "[$(date)] rendering $CLS ($MOD) -> output/$OUT on $(hostname)"
# one media_dir per task: concurrent renders sharing Tex/ and texts/ caches can race
MEDIA="$MANIM_CRYOEM_SCRATCH/media_final/$MOD"
manim render -qh --progress_bar none --media_dir "$MEDIA" "scenes/$MOD.py" "$CLS"
SRC="$MEDIA/videos/$MOD/1080p30/$CLS.mp4"
cp "$SRC" "output/$OUT"
ffprobe -v error -select_streams v:0 -show_entries stream=codec_name,width,height,r_frame_rate,pix_fmt \
    -show_entries format=duration -of compact "output/$OUT"
echo "[$(date)] done"
