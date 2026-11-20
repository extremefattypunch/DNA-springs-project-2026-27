"""W3 positive control: free 50 bp DNA duplex through the Day et al. Fig S6 core (free gpu_test lane).

Particles come from round-1 2D of the rod stream (P6 J27, 256 -> 128 px, 2.84 A/px, 364 A box) and its blind
consensus labels (results/class2d_labels_W2_2D_rods.csv, 3-reviewer majority):
  target  : rod_no_dot classes (free straight DNA)                  -> homo_abinit K=1
  decoys  : noise_junk, ice_contaminant, dot_pair classes (paper: "decoy" ab-initio classes from junk 2D classes)
  hetero  : target + 3 decoys, on the pooled target+decoy particles, forced hard classification
  NU      : particles of the duplex class, C1, default dynamic mask
  checks  : validation (gold-standard FSC with a fresh mask), orientation_diagnostics (cFAR/SCF), local_resolution
Expected: a ~20 A-wide, ~170 A rod; major/minor grooves (34 A pitch) only if the map reaches ~10 A or better.
Strong preferred (side-view) orientation is expected for a rod: judge resolution together with cFAR.
  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/31_duplex_control.py   (or via slurm/driver.sbatch)
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

LABELS = os.path.join(R.REPO, "results", "class2d_labels_W2_2D_rods.csv")


def classes(*labels):
    return [int(r["class"]) for r in csv.DictReader(open(LABELS)) if r["consensus_label"] in labels]


def by_class(cls_list):
    s = set(cls_list)
    return lambda ds: np.isin(np.asarray(ds["alignments2D/class"]), list(s))


def one_volume(job):
    job.refresh()
    outs = list(job.outputs.keys())
    for name in ("volume_class_0", "volume", "volumes_all_classes"):
        if name in outs:
            return name
    raise SystemExit(f"no volume output on {job.uid}: {outs}")


def main():
    cs, st = R.session()
    sets = {"duplex": classes("rod_no_dot"), "junk": classes("noise_junk"), "ice": classes("ice_contaminant"),
            "pairs": classes("dot_pair")}
    print({k: len(v) for k, v in sets.items()}, flush=True)
    subs = {k: R.save_subset(cs, st, f"W3_{k}_particles", "W3", [("W2_2D_rods", "particles")], by_class(v),
                             "particle", "particles", f"W3 rod-stream particles in '{k}' 2D classes",
                             extra={"classes": v}) for k, v in sets.items()}
    pooled = R.save_subset(cs, st, "W3_pooled_particles", "W3", [("W2_2D_rods", "particles")],
                           by_class(sum(sets.values(), [])), "particle", "particles",
                           "W3 pooled duplex + decoy-class particles for hetero refinement")
    ab = {}
    for k in sets:
        ab[k] = R.run(cs, st, f"W3_abinit_{k}", "W3", "homo_abinit",
                      connections={"particles": (subs[k].uid, "particles")},
                      params=dict(abinit_K=1, abinit_max_res=12.0, abinit_init_res=35.0,
                                  abinit_num_particles=60000 if k == "duplex" else 20000),
                      title=f"W3 ab initio K=1: {k}" + ("" if k == "duplex" else " (decoy)"),
                      walltime="04:00:00", num_gpu=1, block=False)
    for k, j in ab.items():
        ab[k] = R.wait(cs, st["project"], j.uid, f"W3_abinit_{k}", st)
    vols = [(ab[k].uid, one_volume(ab[k])) for k in ("duplex", "junk", "ice", "pairs")]
    het = R.run(cs, st, "W3_hetero", "W3", "hetero_refine",
                connections={"particles": (pooled.uid, "particles"), "volume": vols},
                params=dict(multirefine_force_hard_class=1, multirefine_res_init=20),
                title="W3 hetero refine: duplex + 3 decoys (junk, ice, dot pairs)", walltime="06:00:00", num_gpu=1)
    het.refresh()
    print("hetero outputs:", list(het.outputs.keys()), flush=True)
    pcls = [o for o in het.outputs.keys() if o.startswith("particles_class_")]
    vcls = [o for o in het.outputs.keys() if o.startswith("volume_class_")]
    counts = {o: het.outputs[o].num_items for o in pcls}
    st["jobs"]["W3_hetero"]["class_counts"] = counts
    R.save_state(st)
    print("hetero class sizes:", counts, flush=True)
    # class 0 was seeded by the duplex ab initio; check it is still the largest non-decoy-like class by eye in the GUI
    nu = R.run(cs, st, "W3_nu_duplex", "W3", "nonuniform_refine_new",
               connections={"particles": (het.uid, "particles_class_0"), "volume": (het.uid, "volume_class_0")},
               params=dict(refine_symmetry="C1", refine_res_init=20),
               title="W3 NU refine: free 50 bp duplex (control)", walltime="06:00:00", num_gpu=1)
    R.run(cs, st, "W3_validation", "W3", "validation",
          connections={"volume": (nu.uid, "volume")}, params=dict(validate_generate_new_mask=1),
          title="W3 validation: gold-standard FSC, fresh mask", walltime="02:00:00", num_gpu=1)
    R.run(cs, st, "W3_orientdiag", "W3", "orientation_diagnostics",
          connections={"volume": (nu.uid, "volume"), "particles": (nu.uid, "particles"), "mask": (nu.uid, "mask")},
          title="W3 orientation diagnostics (cFAR, SCF)", walltime="02:00:00", num_gpu=1)
    R.run(cs, st, "W3_locres", "W3", "local_resolution",
          connections={"volume": (nu.uid, "volume"), "mask": (nu.uid, "mask")},
          title="W3 local resolution", walltime="02:00:00", num_gpu=1)
    print("W3 duplex control complete", flush=True)


if __name__ == "__main__":
    main()
