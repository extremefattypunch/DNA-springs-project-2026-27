"""Split test for the dot-on-rod NU map (P6 J56): reconstruct with the EXISTING NU poses from
(a) particles of the A1 dot-at-rod-end 2D classes only, (b) particles of the haze classes 105/12 only.
Both skeptics of the 3D length analysis (results/dotrod_3d_length_verified.json) named this as the decisive check:
45 % of the NU particles came from haze classes whose poses re-project poorly, and every map shows a weak kinked
'far arm' of unknown origin. If the far arm vanishes from (a) it came from haze particles; if it persists in (a) it
belongs to the dot-on-rod population. homo_reconstruct = no further alignment, same mask as the NU refinement.
  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/36_split_reconstruct.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

SEL = json.load(open(os.path.join(R.REPO, "results", "A1_length_vs_hinge.json")))["selection"]
DOT = [s["cls"] for s in SEL if s.get("group") == "DOT"]
HAZE = [105, 12]


def main():
    cs, st = R.session()
    nu = st["jobs"]["W4_nu_dotrod"]["uid"]
    ds = cs.find_job(st["project"], nu).load_output("particles", slots=["alignments2D"])
    print("NU particles", len(ds), "class field present:", "alignments2D/class" in ds.fields(), flush=True)
    for name, cls in (("dot", DOT), ("haze", HAZE)):
        sub = R.save_subset(cs, st, f"W4_nu_{name}_only", "W4", [("W4_nu_dotrod", "particles")],
                            lambda d, c=cls: np.isin(np.asarray(d["alignments2D/class"]), c), "particle", "particles",
                            f"W4 NU particles from A1 '{name}' 2D classes only (NU poses kept)")
        R.run(cs, st, f"W4_recon_{name}", "W4", "homo_reconstruct",
              connections={"particles": (sub.uid, "particles"), "mask": (nu, "mask")},
              params=dict(refine_symmetry="C1"),
              title=f"W4 reconstruct-only with NU poses: {name} particles", walltime="01:00:00", num_gpu=1)
    print("split reconstruct complete", flush=True)


if __name__ == "__main__":
    main()
