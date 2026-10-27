#!/bin/bash -l
# Final 1080p30 renders of the six clips, one after another in a single job on the free
# `test` partition (CPU-only; Manim's Cairo renderer is single-threaded, +1 core for x264).
# Sequential rather than an array: the short clips would finish far inside a shared
# array -t and be flagged by the Job Defense Shield (elapsed < 50% of the limit).
#   sbatch cryoEM/manim_cryoem_math/slurm/render_final.sbatch                 # all six
#   sbatch --export=ALL,ONLY=2 -t 0-00:30 cryoEM/manim_cryoem_math/slurm/render_final.sbatch  # clip 03 only
# -t from measured timings: 720p30 renders at ~2.6x real time, 1080p ~2x that; the six
# clips are ~300 s of video, so ~25-40 min in total.
#SBATCH -J manim-final
#SBATCH -p test
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 0-00:55
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
source /etc/profile.d/lmod.sh
ROOT=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/manim_cryoem_math
cd "$ROOT"
set +u; source env/activate_manim.sh; set -u
TMP=/scratch/${USER}_${SLURM_JOB_ID}
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
export TMPDIR="$TMP"

python main.py --list | { [ -n "${ONLY:-}" ] && sed -n "$((ONLY + 1))p" || cat; } > "$TMP/scenes.txt"
while read -r MOD CLS OUT; do
    echo "[$(date)] rendering $CLS ($MOD) -> output/$OUT on $(hostname)"
    # one media_dir per clip keeps Tex/ and texts/ caches separate from the review renders
    MEDIA="$MANIM_CRYOEM_SCRATCH/media_final/$MOD"
    SECONDS=0
    manim render -qh --progress_bar none --media_dir "$MEDIA" "scenes/$MOD.py" "$CLS" < /dev/null
    cp "$MEDIA/videos/$MOD/1080p30/$CLS.mp4" "output/$OUT"
    echo "[$(date)] $OUT rendered in ${SECONDS}s"
    ffprobe -v error -select_streams v:0 -show_entries stream=codec_name,width,height,r_frame_rate,pix_fmt \
        -show_entries format=duration -of compact "output/$OUT"
done < "$TMP/scenes.txt"
echo "[$(date)] done"
