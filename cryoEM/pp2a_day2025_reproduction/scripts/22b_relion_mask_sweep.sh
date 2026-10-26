#!/bin/bash
# Plan section 5C addendum - how much does the RELION gold-standard resolution depend on the mask?
#
# The deposited FSC curve (CryoSPARC, "corrected" with CryoSPARC's auto mask) reads 3.4 / 3.2 A.
# RELION postprocess with our first mask (raw-map threshold matched to the EMDB contour volume,
# extend 3 px, 6 px soft edge) gave 3.98 / 3.56 A and a masked FSC that dips at low resolution -
# the signature of a mask that cuts through the molecule.  Sweep: mask from a 15 A low-passed raw
# map at thresholds of 0.5x..0.1x the contour-matched threshold, extend 3, soft edge 6/10, record
# FINAL RESOLUTION for each.  Output: $PP2A_DATA/em/relion/<iso>/mask_sweep/sweep.csv
set -euo pipefail
iso=${1:?Aa or Ab}
case $iso in Aa) apix=1.056; contour=0.0777 ;; Ab) apix=1.058; contour=0.154 ;; *) exit 1 ;; esac
: "${PP2A_DATA:?}"
R=$PP2A_DATA/em/relion/$iso; OUT=$R/mask_sweep; mkdir -p "$OUT"; cd "$OUT"
thr0=$("$DNASPRING_PY" - "$R/sharp.mrc" "$R/raw.mrc" $contour <<'EOF'
import sys, numpy as np, mrcfile
sharp = mrcfile.open(sys.argv[1]).data; raw = mrcfile.open(sys.argv[2]).data; c = float(sys.argv[3])
print(f"{np.sort(raw.ravel())[::-1][int((sharp > c).sum())]:.6g}")
EOF
)
[ -f raw_lp15.mrc ] || relion_image_handler --i "$R/raw.mrc" --lowpass 15 --angpix $apix --o raw_lp15.mrc
echo "mask,threshold_frac,threshold,extend_px,soft_px,mask_voxels,final_resolution_A,bfactor" > sweep.csv
i=0
for frac in 0.5 0.3 0.2 0.1; do
  for soft in 6 10; do
    i=$((i+1)); thr=$("$DNASPRING_PY" -c "print($thr0*$frac)")
    tag=m${i}_f${frac}_s${soft}
    relion_mask_create --i raw_lp15.mrc --o mask_$tag.mrc --ini_threshold $thr --extend_inimask 3 --width_soft_edge $soft --angpix $apix > /dev/null 2>&1
    nvox=$("$DNASPRING_PY" -c "import mrcfile,sys; print(int((mrcfile.open('mask_$tag.mrc').data>0.5).sum()))")
    relion_postprocess --i "$R/half1.mrc" --i2 "$R/half2.mrc" --mask mask_$tag.mrc --angpix $apix --auto_bfac --autob_lowres 10 --o ./pp_$tag > pp_$tag.log 2>&1 || true
    res=$(grep -E "FINAL RESOLUTION" pp_$tag.log | awk '{print $NF}' | tail -1)
    bf=$(grep -iE "B-factor" pp_$tag.log | grep -oE "\-?[0-9]+\.[0-9]+" | tail -1)
    echo "$tag,$frac,$thr,3,$soft,$nvox,$res,$bf" | tee -a sweep.csv
    rm -f pp_${tag}_masked.mrc pp_$tag.mrc   # keep star/xml/log only
  done
done
echo "deposited (CryoSPARC corrected) 0.143 resolution for reference: Aa 3.4, Ab 3.2"
