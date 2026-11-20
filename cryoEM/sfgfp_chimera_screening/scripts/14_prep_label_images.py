"""Make the labelling images for selections/label_set.csv (CPU only).

For each micrograph: Fourier-crop 2048 -> 1024 (2.84 A/px), Topaz pretrained denoise (CPU, `topaz-cs denoise`),
and a plain 15 A Gaussian low-pass view. Writes 8-bit JPEGs (1024x1024; label coordinates in this frame x2 =
CryoSPARC micrograph pixels at 1.42 A/px) to $SFGFP_DATA/label_images/{mic_id}_{dn,lp}.jpg and a manifest.json.
  sbatch -p test -c 8 --mem=16G -t 1:00:00 --wrap "$DNASPRING_PY -I scripts/14_prep_label_images.py"
"""
import csv
import json
import os
import subprocess
import tempfile

import mrcfile
import numpy as np
from PIL import Image

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")
OUT = os.path.join(DATA, "label_images")
TOPAZ = "/n/holylabs/hekstra_lab/Everyone/ianpoon/software/bin/topaz-cs"
APIX_IN, BIN = 1.42, 2


def fcrop(a, n):
    F = np.fft.fftshift(np.fft.fft2(a))
    c = a.shape[0] // 2
    F = F[c - n // 2: c + n // 2, c - n // 2: c + n // 2]
    return np.real(np.fft.ifft2(np.fft.ifftshift(F))).astype(np.float32) * (n / a.shape[0]) ** 2


def lowpass(a, apix, res):
    ky = np.fft.fftfreq(a.shape[0], apix)[:, None]
    kx = np.fft.rfftfreq(a.shape[1], apix)[None, :]
    F = np.fft.rfft2(a) * np.exp(-(ky ** 2 + kx ** 2) * res ** 2)   # Gaussian, 1/e at 1/res
    return np.fft.irfft2(F, s=a.shape)


def to_jpg(a, path):
    lo, hi = np.percentile(a, [0.5, 99.5])
    u8 = (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
    Image.fromarray(u8).save(path, quality=90)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(REPO, "selections", "label_set.csv"))))
    tmp = tempfile.mkdtemp(dir=os.environ.get("TMPDIR", "/tmp"))
    inputs = []
    for r in rows:
        a = mrcfile.read(r["mrc"]).astype(np.float32)
        b = fcrop(a, a.shape[0] // BIN)
        p = os.path.join(tmp, f"{r['mic_id']}.mrc")
        with mrcfile.new(p, overwrite=True) as m:
            m.set_data(b)
            m.voxel_size = APIX_IN * BIN
        inputs.append(p)
        to_jpg(lowpass(b, APIX_IN * BIN, 15.0), os.path.join(OUT, f"{r['mic_id']}_lp.jpg"))
    dn = os.path.join(tmp, "dn")
    os.makedirs(dn, exist_ok=True)
    subprocess.run([TOPAZ, "denoise", "--device", "-1", "-o", dn, *inputs], check=True)
    for r in rows:
        d = mrcfile.read(os.path.join(dn, f"{r['mic_id']}.mrc")).astype(np.float32)
        to_jpg(d, os.path.join(OUT, f"{r['mic_id']}_dn.jpg"))
    man = [{"mic_id": r["mic_id"], "split": r["split"], "square": r["square"],
            "target_defocus_um": float(r["target_defocus_um"]), "stem": r["stem"],
            "scale_to_cryosparc_px": BIN, "apix_display": APIX_IN * BIN} for r in rows]
    json.dump(man, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    print("wrote", len(rows), "x2 images to", OUT)


if __name__ == "__main__":
    main()
