#!/bin/bash
# After render_final.sbatch: last-frame thumbnails, concatenated master video and
# an ffprobe report. Light work; run inline:  bash cryoEM/manim_cryoem_math/slurm/postprocess.sh
set -euo pipefail
ROOT=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/manim_cryoem_math
cd "$ROOT"
set +u; source env/activate_manim.sh; set -u
mkdir -p output/thumbnails
LIST=$(mktemp)
trap 'rm -f "$LIST"' EXIT
python main.py --list | while read -r MOD CLS OUT; do
    f="output/$OUT"
    [ -f "$f" ] || { echo "missing $f" >&2; exit 1; }
    # last decodable frame = the takeaway freeze
    ffmpeg -loglevel error -y -sseof -0.5 -i "$f" -update 1 -frames:v 1 "output/thumbnails/${OUT%.mp4}.png"
    echo "file '$ROOT/$f'" >> "$LIST"
done
ffmpeg -loglevel error -y -f concat -safe 0 -i "$LIST" -c copy output/all_scenes.mp4
{
  printf "%-36s %6s %10s %8s %7s %9s %s\n" file width height fps codec pix_fmt duration_s
  for f in output/*.mp4; do
    ffprobe -v error -select_streams v:0 \
      -show_entries stream=width,height,r_frame_rate,codec_name,pix_fmt -show_entries format=duration \
      -of default=nw=1:nk=1 "$f" | paste -sd' ' | awk -v f="$(basename "$f")" '{printf "%-36s %6s %10s %8s %7s %9s %s\n", f, $2, $3, $5, $1, $4, $6}'
    if ffprobe -v error -select_streams a -show_entries stream=index -of csv=p=0 "$f" | grep -q .; then
      echo "  WARNING: $f has an audio stream"
    fi
  done
} | tee output/render_report.txt
