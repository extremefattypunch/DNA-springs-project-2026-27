#!/usr/bin/env python
"""Register variants of the fitted 32 bp duplex (scripts/53_figS5_fit_W3_dna.py), used to test whether the W3
map-model FSC depends on the helical register (groove phase / polarity), which the fit does not determine.
Variants (no refitting): spin 90 and 180 deg about the fitted helix axis through the model centre, and a polarity
flip (180 deg about an axis perpendicular to the helix axis). Run: $DNASPRING_PY -I scripts/53b_figS5_register_variants.py"""
import json
import numpy as np
from scipy.spatial.transform import Rotation
S = "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W3_duplex/fit"
lines = open(f"{S}/dna_32bp_fit.pdb").read().splitlines()
atom = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
X = np.array([[float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54])] for i in atom])
P = np.array([X[k] for k, i in enumerate(atom) if lines[i][12:16].strip() == "P"])
c = X.mean(0)
w, V = np.linalg.eigh(np.cov((P - P.mean(0)).T))
ax = V[:, 2]
perp = np.cross(ax, [0, 0, 1.0]); perp /= np.linalg.norm(perp)
var = {"spin090": Rotation.from_rotvec(ax * np.pi / 2), "spin180": Rotation.from_rotvec(ax * np.pi),
       "flip": Rotation.from_rotvec(perp * np.pi)}
for name, R in var.items():
    Y = R.apply(X - c) + c
    out = list(lines)
    for k, i in enumerate(atom):
        l = lines[i]
        out[i] = l[:30] + f"{Y[k,0]:8.3f}{Y[k,1]:8.3f}{Y[k,2]:8.3f}" + l[54:]
    open(f"{S}/dna_32bp_fit_{name}.pdb", "w").write("\n".join(out) + "\n")
json.dump({"axis": ax.tolist(), "centre": c.tolist(), "variants": list(var)}, open(f"{S}/register_variants.json", "w"), indent=1)
print("axis", ax, "centre", c)
