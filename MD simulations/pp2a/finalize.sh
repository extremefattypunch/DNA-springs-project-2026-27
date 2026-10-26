#!/bin/bash
# Everything after the trajectories: analyse, aggregate, draw, animate, write the wiki page,
# and publish the figures to docs/ (GitHub Pages serves docs/figures, which is where the
# wiki page's images point -- the repository is private, so raw.githubusercontent links
# do not render).  Re-runnable; each step overwrites its own outputs.
set -euo pipefail
cd "$(dirname "$0")"
source ./env/activate.sh
P="$DNASPRING_PY"
ROOT="$(cd ../.. && pwd)"
STRIDE=${STRIDE:-2}
FRAMES=${FRAMES:-80}

echo "##### collect"
$P analysis/collect.py --stride "$STRIDE" ${ONLY:+--only $ONLY}
echo "##### force windows"
$P analysis/force_windows.py || echo "  (no finished windows yet)"
echo "##### findings"
$P analysis/findings.py
echo "##### figures"
$P analysis/figures.py

echo "##### animations"
for sys in ${ANIMATE:-P0_apo P3_ss50 P2_clicked}; do
  if [ "$sys" = P0_apo ]; then
    traj="$PP2A_DATA/md/runs/Aa_rep1/traj.dcd"; top="$PP2A_DATA/md/systems/Aa/system.prmtop"; rm_args=(); views="overview strain"
  else
    traj="$PP2A_MD_SCRATCH/runs/$sys/rep1/03_production/traj.dcd"; top="build/systems/$sys/system.prmtop"
    rm_args=(--resmap "build/systems/$sys/residue_map.json"); views="overview sites strain"
    [ "$sys" = P2_clicked ] && views="overview sites"
  fi
  [ -f "$traj" ] || { echo "  $sys: no trajectory"; continue; }
  $P analysis/animate.py --traj "$traj" --top "$top" "${rm_args[@]}" --name "$sys" --frames "$FRAMES" --views $views
done
# small GIFs for the wiki page (the mp4s stay in figures/ for anyone who wants them)
mkdir -p figures/web
for f in figures/*_overview.gif figures/*_sites.gif figures/*_strain.gif; do
  [ -f "$f" ] || continue
  b=$(basename "$f" .gif)
  ffmpeg -y -loglevel error -i "${f%.gif}.mp4" \
     -vf "fps=10,scale=520:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=64[p];[b][p]paletteuse=dither=none" \
     -loop 0 "figures/web/pp2a_$b.gif"
done

echo "##### publish figures to docs/"
mkdir -p "$ROOT/docs/figures"
cp figures/pp2a_*.png "$ROOT/docs/figures/" 2>/dev/null || true
for f in figures/*_still.png; do [ -f "$f" ] && cp "$f" "$ROOT/docs/figures/pp2a_$(basename "$f")"; done
cp figures/web/pp2a_*.gif "$ROOT/docs/figures/" 2>/dev/null || true
touch "$ROOT/docs/.nojekyll"

echo "##### wiki page  (after the figures are in docs/, which is where it looks for them)"
$P report/build_wiki_pp2a.py
echo "done -- commit docs/figures and docs/wiki, then push the wiki page to the wiki repository"
