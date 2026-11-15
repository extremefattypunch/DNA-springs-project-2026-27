"""W4: 3D of the dot-on-rod particles (Fig S6 core), to measure attached-DNA length and geometry in 3D.

Why: A1 (report §13a) showed DNA on sfGFP projects 16-33 A shorter than free rods, but a 2D projection cannot tell a
genuinely shorter/incomplete duplex (hypothesis 1, yielded loops) from a full duplex tilted out of plane (2/3).
Particles: A1 window-free 2D (P6 J45; 454 A box, 2.84 A/px) classes, grouped by the A1 analysis
(results/A1_length_vs_hinge.json "selection"):
  target  DOT  : 15 dot-at-rod-end classes (76,607)                       -> ab initio K=1
  competitor FREE: 16 dot-free straight-rod classes (11,995)              -> ab initio K=1 (a real class, not a decoy)
  decoys       : 'haze' classes 105 + 12 (lone dot, rod not averaged) and the large ice/disc classes
Then: hetero refinement (all four volumes, hard classification) on the pooled particles -> NU refine of the
dot-on-rod class -> validation / orientation diagnostics -> class_3D in PCA mode (10 classes, 3 components, as in
Day et al. Fig S7, at 12 A because the consensus is ~8-12 A) to split short vs long attached DNA in 3D.
  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/35_dot_rod_3d.py   (or via slurm/driver.sbatch)
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

SEL = json.load(open(os.path.join(R.REPO, "results", "A1_length_vs_hinge.json")))["selection"]
DOT = [s["cls"] for s in SEL if s.get("group") == "DOT"]
FREE = [s["cls"] for s in SEL if s.get("group") == "FREE"]
HAZE = [105, 12]
ICE = [58, 56, 67, 117, 23, 25, 37, 88]


def by_class(cls_list):
    s = list(set(cls_list))
    return lambda ds: np.isin(np.asarray(ds["alignments2D/class"]), s)


def vol_out(job):
    job.refresh()
    for name in ("volume_class_0", "volume", "volumes_all_classes"):
        if name in job.outputs.keys():
            return name
    raise SystemExit(f"no volume output on {job.uid}")


def main():
    cs, st = R.session()
    sets = {"dot": DOT, "free": FREE, "haze": HAZE, "ice": ICE}
    print({k: len(v) for k, v in sets.items()}, flush=True)
    subs = {k: R.save_subset(cs, st, f"W4_{k}_particles", "W4", [("A1_2D", "particles")], by_class(v),
                             "particle", "particles", f"W4 A1-2D particles in '{k}' classes", extra={"classes": v})
            for k, v in sets.items()}
    pooled = R.save_subset(cs, st, "W4_pooled_particles", "W4", [("A1_2D", "particles")],
                           by_class(DOT + FREE + HAZE + ICE), "particle", "particles",
                           "W4 pooled dot-on-rod + free rods + decoy particles")
    nmax = {"dot": 60000, "free": 12000, "haze": 20000, "ice": 9000}
    ab = {k: R.run(cs, st, f"W4_abinit_{k}", "W4", "homo_abinit",
                   connections={"particles": (subs[k].uid, "particles")},
                   params=dict(abinit_K=1, abinit_max_res=12.0, abinit_init_res=35.0, abinit_num_particles=nmax[k]),
                   title=f"W4 ab initio K=1: {k}" + (" (decoy)" if k in ("haze", "ice") else ""),
                   walltime="04:00:00", num_gpu=1, block=False) for k in sets}
    for k, j in ab.items():
        ab[k] = R.wait(cs, st["project"], j.uid, f"W4_abinit_{k}", st)
    vols = [(ab[k].uid, vol_out(ab[k])) for k in ("dot", "free", "haze", "ice")]
    het = R.run(cs, st, "W4_hetero", "W4", "hetero_refine",
                connections={"particles": (pooled.uid, "particles"), "volume": vols},
                params=dict(multirefine_force_hard_class=1, multirefine_res_init=20),
                title="W4 hetero refine: dot-on-rod + free rod + 2 decoys", walltime="08:00:00", num_gpu=1)
    het.refresh()
    counts = {o: het.outputs[o].num_items for o in het.outputs.keys() if o.startswith("particles_class_")}
    st["jobs"]["W4_hetero"]["class_counts"] = counts
    R.save_state(st)
    print("hetero class sizes:", counts, flush=True)
    nu = R.run(cs, st, "W4_nu_dotrod", "W4", "nonuniform_refine_new",
               connections={"particles": (het.uid, "particles_class_0"), "volume": (het.uid, "volume_class_0")},
               params=dict(refine_symmetry="C1", refine_res_init=20),
               title="W4 NU refine: sfGFP + attached DNA (dot-on-rod)", walltime="06:00:00", num_gpu=1)
    R.run(cs, st, "W4_validation", "W4", "validation", connections={"volume": (nu.uid, "volume")},
          params=dict(validate_generate_new_mask=1), title="W4 validation (GSFSC, fresh mask)",
          walltime="02:00:00", num_gpu=1)
    R.run(cs, st, "W4_orientdiag", "W4", "orientation_diagnostics",
          connections={"volume": (nu.uid, "volume"), "particles": (nu.uid, "particles"), "mask": (nu.uid, "mask")},
          title="W4 orientation diagnostics (cFAR, SCF)", walltime="02:00:00", num_gpu=1)
    R.run(cs, st, "W4_class3D_pca", "W4", "class_3D",
          connections={"particles": (nu.uid, "particles"), "mask": (nu.uid, "mask")},
          params=dict(class3D_N_K=10, class3D_init_mode="PCA", class3D_PCA_num_components=3,
                      class3D_target_res=12, class3D_force_hard_class=1),
          title="W4 class_3D PCA: 10 classes, 3 components, 12 A (Fig S7-style)", walltime="08:00:00", num_gpu=1)
    print("W4 dot-on-rod 3D complete", flush=True)


if __name__ == "__main__":
    main()
