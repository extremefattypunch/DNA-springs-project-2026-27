#!/bin/bash
# Plan section 5D - EMReady v1.1 (the paper's sharpening tool) on the half-map average.
#
# Question G of the report: does EMReady applied to (half1+half2)/2 reproduce the deposited
# (EMReady-sharpened) primary map?  Output per isoform in $PP2A_DATA/em/emready/<iso>/:
#   emready_raw.mrc      EMReady output on the raw-map proxy (1.0 A grid model, default stride 12)
#   compare.json         real-space CC and FSC(0.5) between emready_raw.mrc and the deposited map
# Usage (GPU, MIG slice is plenty for a 256^3 box):
#   sbatch -p gpu_test --gres=gpu:1 -c 4 --mem=24G -t 0-01:00 slurm/pp2a_cpu.sbatch bash scripts/23_emready.sh Aa
set -euo pipefail
iso=${1:?Aa or Ab}
case $iso in Aa) emd=48300; apix=1.056 ;; Ab) emd=48224; apix=1.058 ;; *) exit 1 ;; esac
: "${PP2A_DATA:?}" "${EMREADY_HOME:?EMReady not installed}"
E=$PP2A_DATA/em/emdb/$iso; R=$PP2A_DATA/em/relion/$iso; OUT=$PP2A_DATA/em/emready/$iso; mkdir -p "$OUT"; cd "$OUT"
RAW=$R/raw.mrc
if [ ! -f "$RAW" ]; then   # RELION not run yet: make the half-map average with numpy
  RAW=$OUT/raw.mrc
  "$DNASPRING_PY" - "$E/emd_${emd}_half_map_1.map" "$E/emd_${emd}_half_map_2.map" "$RAW" <<'EOF'
import sys, numpy as np, mrcfile
a = mrcfile.open(sys.argv[1]); b = mrcfile.open(sys.argv[2])
with mrcfile.new(sys.argv[3], overwrite=True) as o:
    o.set_data(((a.data.astype(np.float32) + b.data.astype(np.float32)) * 0.5)); o.voxel_size = a.voxel_size; o.header.origin = a.header.origin
EOF
fi
GPU=${CUDA_VISIBLE_DEVICES:-0}; [ -n "$GPU" ] && gopt="-g 0" || gopt="--use_cpu"
if [ ! -f emready_raw.mrc ]; then
  bash "$EMREADY_HOME/EMReady.sh" "$RAW" emready_raw.mrc $gopt -s 12 -b 10 2>&1 | tail -20
fi
# compare with the deposited sharpened map: resample EMReady output (1.0 A grid) back onto the deposited grid
"$DNASPRING_PY" - "$E/emd_${emd}.map" emready_raw.mrc compare.json <<'EOF'
import sys, json, numpy as np, mrcfile
from scipy.ndimage import zoom
dep = mrcfile.open(sys.argv[1]); em = mrcfile.open(sys.argv[2])
d = dep.data.astype(np.float32); e = em.data.astype(np.float32)
vd = float(dep.voxel_size.x); ve = float(em.voxel_size.x)
# EMReady writes a 1 A grid covering the same physical box (origin preserved); rescale to the deposited grid
if e.shape != d.shape:
    e = zoom(e, [ve / vd] * 3, order=1)
    # crop / pad to identical shape
    out = np.zeros_like(d); s = [min(a, b) for a, b in zip(d.shape, e.shape)]
    out[:s[0], :s[1], :s[2]] = e[:s[0], :s[1], :s[2]]; e = out
mask = d > np.percentile(d, 99)      # region of the deposited density
cc_all = float(np.corrcoef(d.ravel(), e.ravel())[0, 1]); cc_mask = float(np.corrcoef(d[mask], e[mask])[0, 1])
# FSC between the two maps
F1, F2 = np.fft.fftn(d), np.fft.fftn(e)
n = d.shape[0]; fr = np.fft.fftfreq(n)
kx, ky, kz = np.meshgrid(fr, fr, fr, indexing="ij"); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
bins = np.minimum((k * n).astype(int), n // 2)
num = np.bincount(bins.ravel(), (F1 * np.conj(F2)).real.ravel()); den = np.sqrt(np.bincount(bins.ravel(), np.abs(F1).ravel() ** 2) * np.bincount(bins.ravel(), np.abs(F2).ravel() ** 2))
fsc = num / np.maximum(den, 1e-12); shells = np.arange(len(fsc)) / (n * vd)
res05 = float(1 / shells[np.argmax(fsc[1:] < 0.5) + 1]) if (fsc[1:] < 0.5).any() else float("inf")
json.dump({"cc_all": cc_all, "cc_top1pct": cc_mask, "fsc05_between_maps_A": res05, "shape_dep": d.shape, "apix_dep": vd, "apix_emready": ve}, open(sys.argv[3], "w"), indent=1)
print(json.load(open(sys.argv[3])))
EOF
