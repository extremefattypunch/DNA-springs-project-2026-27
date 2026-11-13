"""Choose 40 motion-corrected micrographs for Ian's hand labelling (ground truth for the census, Topaz, ring rate).

Stratified over grid squares and EPU target defocus (AppliedDefocus from each movie's xml; CTF for the full set
does not exist yet). 15 are flagged as a frozen held-out test set (never used for training).
Writes selections/label_set.csv (mic_id, square, target_defocus_um, split, mrc path).
  $DNASPRING_PY -I scripts/13_select_label_set.py
"""
import csv
import glob
import os
import random
import re

PROJ = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-cs-sfgfp-chimera-screening"
RAW = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/raw/Images-Disc1"
REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
N, N_TEST, SEED = 40, 15, 20261113

stem2sq = {}
for t in glob.glob(os.path.join(RAW, "GridSquare_*", "Data", "*_Fractions.tiff")):
    stem = os.path.basename(t)[: -len("_Fractions.tiff")]
    stem2sq[stem] = (t.split("/")[-3], t[: -len("_Fractions.tiff")] + ".xml")

rows = []
for job in ("J10", "J12", "J3"):
    for m in glob.glob(os.path.join(PROJ, job, "motioncorrected", "*_patch_aligned_doseweighted.mrc")):
        stem = re.sub(r"^\d+_", "", os.path.basename(m)).replace("_Fractions_patch_aligned_doseweighted.mrc", "")
        if stem not in stem2sq:
            continue
        sq, xml = stem2sq[stem]
        s = open(xml).read()
        mm = re.search(r"AppliedDefocus</a:Key><a:Value[^>]*>([-\d.E]+)<", s)
        df = abs(float(mm.group(1))) * 1e6 if mm else float("nan")
        rows.append({"stem": stem, "square": sq, "target_defocus_um": round(df, 2), "mrc": m, "job": job})

rng = random.Random(SEED)
by_sq = {}
for r in rows:
    by_sq.setdefault(r["square"], []).append(r)
squares = sorted(by_sq)
picked = []
# round-robin over squares; inside a square alternate defocus targets
while len(picked) < N and any(by_sq.values()):
    for sq in squares:
        pool = by_sq[sq]
        if not pool or len(picked) >= N:
            continue
        dfs = sorted({r["target_defocus_um"] for r in pool})
        want = dfs[len([p for p in picked if p["square"] == sq]) % len(dfs)]
        cand = [r for r in pool if r["target_defocus_um"] == want]
        r = rng.choice(cand)
        pool.remove(r)
        picked.append(r)
rng.shuffle(picked)
for i, r in enumerate(picked):
    r["mic_id"] = f"m{i:02d}"
    r["split"] = "test" if i < N_TEST else "train"
os.makedirs(os.path.join(REPO, "selections"), exist_ok=True)
with open(os.path.join(REPO, "selections", "label_set.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, ["mic_id", "split", "square", "target_defocus_um", "job", "stem", "mrc"])
    w.writeheader()
    for r in sorted(picked, key=lambda r: r["mic_id"]):
        w.writerow({k: r[k] for k in w.fieldnames})
print(f"{len(rows)} candidate micrographs over {len(squares)} squares; picked {len(picked)} "
      f"from {len({r['square'] for r in picked})} squares; defocus targets {sorted({r['target_defocus_um'] for r in picked})}")
