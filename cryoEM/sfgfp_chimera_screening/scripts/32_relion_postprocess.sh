#!/bin/bash
# RELION 3.1 post-processing of a CryoSPARC refinement's half maps (generalised from
# ../pp2a_day2025_reproduction/scripts/22_relion_postprocess.sh).
#
#   bash scripts/32_relion_postprocess.sh TAG HALF_A HALF_B APIX CS_MASK
# Writes $SFGFP_DATA/relion/TAG/:
#   raw.mrc                    (half A + half B) / 2
#   pp_csmask/                 FSC with CryoSPARC's validation mask: unmasked / masked / phase-randomisation corrected
#   mask_loose.mrc, pp_loose/  same with a loose RELION mask (raw map 15 A low-pass, threshold at the 98.5th
#                              percentile, extend 6 px, 8 px soft edge) to expose mask-driven FSC inflation
#   relion_locres.mrc          RELION local resolution (as in Day et al. 2025), loose mask
set -euo pipefail
tag=${1:?TAG}; HA=${2:?half A}; HB=${3:?half B}; apix=${4:?apix}; CSMASK=${5:?cryosparc mask}
: "${SFGFP_DATA:?source env/activate_sfgfp.sh first}"
command -v relion_postprocess >/dev/null || { echo "relion not on PATH"; exit 1; }
OUT=$SFGFP_DATA/relion/$tag; mkdir -p "$OUT"; cd "$OUT"
ln -sf "$HA" half1.mrc; ln -sf "$HB" half2.mrc; ln -sf "$CSMASK" csmask.mrc
[ -f raw.mrc ] || relion_image_handler --i half1.mrc --add half2.mrc --multiply_constant 0.5 --o raw.mrc --angpix "$apix"
mkdir -p pp_csmask pp_loose
relion_postprocess --i half1.mrc --i2 half2.mrc --mask csmask.mrc --angpix "$apix" --auto_bfac --autob_lowres 10 \
  --o ./pp_csmask/postprocess > pp_csmask/postprocess.log 2>&1 || grep -q "FINAL RESOLUTION" pp_csmask/postprocess.log
if [ ! -f mask_loose.mrc ]; then
  relion_image_handler --i raw.mrc --lowpass 15 --angpix "$apix" --o raw_lp15.mrc
  thr=$("$DNASPRING_PY" -I -c "import mrcfile,numpy as np; v=mrcfile.read('raw_lp15.mrc'); print(f'{np.percentile(v,98.5):.6g}')")
  relion_mask_create --i raw_lp15.mrc --o mask_loose.mrc --ini_threshold "$thr" --extend_inimask 6 --width_soft_edge 8 --angpix "$apix"
fi
relion_postprocess --i half1.mrc --i2 half2.mrc --mask mask_loose.mrc --angpix "$apix" --auto_bfac --autob_lowres 10 \
  --o ./pp_loose/postprocess > pp_loose/postprocess.log 2>&1 || grep -q "FINAL RESOLUTION" pp_loose/postprocess.log
relion_postprocess --i half1.mrc --i2 half2.mrc --mask mask_loose.mrc --angpix "$apix" --locres --locres_sampling 25 \
  --locres_maskrad 30 --o ./locres > locres.log 2>&1 || true
[ -f locres_locres.mrc ] && ln -sf locres_locres.mrc relion_locres.mrc
for d in pp_csmask pp_loose; do echo "== $d"; grep -E "FINAL RESOLUTION|apply B-factor" $d/postprocess.log | head -3; done
echo "done: $OUT"
