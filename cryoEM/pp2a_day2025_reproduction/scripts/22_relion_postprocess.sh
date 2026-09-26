#!/bin/bash
# Plan section 5C - RELION 3.1 post-processing of the deposited half maps.
#
# Produces, per isoform, in $PP2A_DATA/em/relion/<iso>/:
#   raw.mrc            (half1+half2)/2  - stand-in for the undeposited unsharpened map
#   mask.mrc           soft solvent mask from raw.mrc (extend 3 px, 6 px cosine edge)
#   postprocess.star   gold-standard FSC (unmasked / masked / phase-randomisation corrected)
#   postprocess.mrc    B-factor-sharpened map (independent comparison to the deposited EMReady map)
#   relion_locres.mrc, locres_filtered.mrc   local resolution as in the paper (RELION 3.1)  [Fig S5A/E]
#   raw_lp6.mrc        6 A low-pass of the raw map                                       [Fig 2B]
# Usage: bash scripts/22_relion_postprocess.sh Aa|Ab
set -euo pipefail
iso=${1:?Aa or Ab}
case $iso in
  Aa) emd=48300; apix=1.056; contour=0.0777 ;;   # EMDB recommended contour of the deposited (sharpened) map
  Ab) emd=48224; apix=1.058; contour=0.154 ;;
  *) echo "iso must be Aa or Ab"; exit 1 ;;
esac
: "${PP2A_DATA:?source env/activate_pp2a.sh first}"
command -v relion_postprocess >/dev/null || { echo "relion not on PATH (build job finished?)"; exit 1; }
E=$PP2A_DATA/em/emdb/$iso; OUT=$PP2A_DATA/em/relion/$iso; mkdir -p "$OUT"; cd "$OUT"
# RELION decides the format from the extension and treats ".map" as SPIDER -> use .mrc symlinks
ln -sf "$E/emd_${emd}_half_map_1.map" half1.mrc; ln -sf "$E/emd_${emd}_half_map_2.map" half2.mrc; ln -sf "$E/emd_${emd}.map" sharp.mrc
H1=half1.mrc; H2=half2.mrc
relion_postprocess --version 2>&1 | head -1 || true

# 1. raw map = average of the half maps
[ -f raw.mrc ] || relion_image_handler --i "$H1" --add "$H2" --multiply_constant 0.5 --o raw.mrc --angpix $apix
# 2. mask.  The deposited contour belongs to the sharpened map; for the unsharpened average
#    use a threshold at the same enclosed volume, found by matching the voxel count.
if [ ! -f mask.mrc ]; then
  thr=$("$DNASPRING_PY" - sharp.mrc raw.mrc $contour <<'EOF'
import sys, numpy as np, mrcfile
sharp = mrcfile.open(sys.argv[1]).data; raw = mrcfile.open(sys.argv[2]).data; c = float(sys.argv[3])
n = int((sharp > c).sum())
print(f"{np.sort(raw.ravel())[::-1][n]:.6g}")
EOF
)
  echo "raw-map threshold enclosing the same volume as contour $contour: $thr"
  relion_mask_create --i raw.mrc --o mask.mrc --ini_threshold $thr --extend_inimask 3 --width_soft_edge 6 --angpix $apix --j ${OMP_NUM_THREADS:-1}
fi
# 3. gold-standard FSC + B-factor sharpening
relion_postprocess --i "$H1" --i2 "$H2" --mask mask.mrc --angpix $apix --auto_bfac --autob_lowres 10 --o ./postprocess > postprocess.log 2>&1 || grep -q "FINAL RESOLUTION" postprocess.log   # rootname needs "./": RELION writes <dir>/logfile.pdf
grep -E "FINAL RESOLUTION|B-factor|FINAL_RESOLUTION|bfactor" postprocess.log | head -5
# 4. local resolution (paper: RELION 3.1)
relion_postprocess --i "$H1" --i2 "$H2" --mask mask.mrc --angpix $apix --locres --locres_sampling 15 --locres_maskrad 25 --o ./locres > locres.log 2>&1 || ls relion_locres.mrc
ls -la relion_locres.mrc locres_filtered.mrc 2>/dev/null || ls
# 5. Fig 2B input
[ -f raw_lp6.mrc ] || relion_image_handler --i raw.mrc --lowpass 6 --angpix $apix --o raw_lp6.mrc
# also 6 A lowpass of the deposited sharpened map for comparison
[ -f sharp_lp6.mrc ] || relion_image_handler --i sharp.mrc --lowpass 6 --angpix $apix --o sharp_lp6.mrc
echo "done: $OUT"
