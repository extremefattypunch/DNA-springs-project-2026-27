#!/bin/bash
# Re-run every analysis over whatever has finished, redraw the figures and animations,
# and rebuild the report. Safe to run repeatedly while production is still going.
set -euo pipefail
cd "$(dirname "$0")"
source ./env/activate.sh
export PATH="$DNASPRING_ENV/bin:$PATH"
P="$DNASPRING_PY"

STRIDE="${STRIDE:-5}"
$P analysis/collect.py --stride "$STRIDE"
$P analysis/figures.py

# Animate the longest chimera replicate available, plus a control for comparison.
for sys in S3_spring27 S0_wt; do
  traj=$(ls -S data/runs/$sys/rep*/03_production/traj.dcd 2>/dev/null | head -1) || true
  [ -n "${traj:-}" ] || continue
  $P analysis/animate.py --traj "$traj" --top "build/systems/$sys/system.prmtop" \
     --name "$sys" --frames "${FRAMES:-80}"
done

# Web-sized GIFs for the report.  The artifact viewer renders <img> reliably and
# silently drops base64 <video>, so the animations ship as GIF; the full-quality mp4s
# stay in figures/ for slides.  dither=none matters: PyMOL renders are mostly flat
# background, and dithering it costs megabytes for nothing.
mkdir -p figures/web
gifify () {  # view width colors fps crop
  local src="figures/${ANIM:-S3_spring27}_$1.mp4"
  [ -f "$src" ] || return 0
  ffmpeg -y -loglevel error -i "$src" \
    -vf "fps=$4,crop=iw*$5:ih*$5,scale=$2:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=$3[p];[b][p]paletteuse=dither=none" \
    -loop 0 "figures/web/${ANIM:-S3_spring27}_$1.gif"
  printf "  %-12s %6d kB\n" "$1" "$(( $(stat -c%s "figures/web/${ANIM:-S3_spring27}_$1.gif")/1024 ))"
}
gifify overview    500 64 10 0.94
gifify chromophore 400 32  6 0.78
gifify strain      460 48  8 0.94

$P analysis/findings.py
$P report/build_report.py
echo
echo "report: $(pwd)/report/sfgfp_dna_spring_report.html"
