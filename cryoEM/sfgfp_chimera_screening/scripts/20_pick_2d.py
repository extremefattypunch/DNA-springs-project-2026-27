"""W2 (Fig S6 'picking + 2D', adapted to the mixture): curate -> two pick streams -> power filter -> extract -> 2D.

Streams (judge-panel plan, 2026-11-13; see report.md):
  anchors : blob 30-50 A circles, min_distance 1.0  -> every sfGFP-sized dot (one per chimera; census denominator)
  rods    : blob ellipses 20-180 A, min_distance 4  -> DNA-sized elongated objects (free duplex, tails, loops)
Each stream: External power filter (selections/pick_power_filter.json, per 0.25 um defocus bin) -> extract ->
2D classification (rods: K 200, 12 A, recentre, window 200 A; anchors: random 25% of curated micrographs,
K 200, 12 A, no recentre, window 220 A). Exposure curation thresholds: selections/S0_exposures.json.

  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/20_pick_2d.py         (or via slurm/driver.sbatch)
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

SEL = os.path.join(R.REPO, "selections")


def curate_mask(th):
    def f(ds):
        df = (ds["ctf/df1_A"] + ds["ctf/df2_A"]) / 2
        ast = np.abs(ds["ctf/df1_A"] - ds["ctf/df2_A"])
        ok = (ds["ctf/ctf_fit_to_A"] <= th["ctf_fit_to_A_max"]) & (ast <= th["astig_A_max"])
        ok &= (df >= th["defocus_A_min"]) & (df <= th["defocus_A_max"])
        ok &= ds["ctf_stats/ice_thickness_rel"] <= th["ice_thickness_rel_max"]
        return ok
    return f


def power_mask(exps_ds, keep_frac):
    """Keep the top `keep_frac` of picks by local power WITHIN each 0.25 um defocus bin, computed on this stream's own
    picks. (The pilot's absolute thresholds, selections/pick_power_filter.json, were calibrated on 80-140 A circular
    blobs and do not transfer: the 20-180 A ellipse picker's power runs ~2.5x lower, so they kept only 1.3% of picks.)"""
    uid2df = dict(zip(exps_ds["uid"].tolist(), ((exps_ds["ctf/df1_A"] + exps_ds["ctf/df2_A"]) / 2e4).tolist()))

    def f(ds):
        df = np.array([uid2df.get(u, np.nan) for u in ds["location/micrograph_uid"].tolist()])
        pw = np.asarray(ds["pick_stats/power"])
        keep = np.zeros(len(pw), bool)
        b = np.floor(df / 0.25)
        for v in np.unique(b[np.isfinite(b)]):
            s = b == v
            keep[s] = pw[s] > np.quantile(pw[s], 1 - keep_frac)
        return keep
    return f


def main():
    cs, st = R.session()
    th = json.load(open(os.path.join(SEL, "S0_exposures.json")))
    pf = json.load(open(os.path.join(SEL, "pick_power_filter_v2.json")))
    cur = R.save_subset(cs, st, "W2_curated", "W1", [("W1_ctf_0", "exposures"), ("W1_ctf_1", "exposures")],
                        curate_mask(th), "exposure", "exposures",
                        "W1 curated exposures (selections/S0_exposures.json)", extra={"thresholds": th})
    exps = cur.load_output("exposures")
    print(f"curated: {len(exps)} micrographs", flush=True)

    streams = {
        "rods": dict(pick=dict(diameter=20, diameter_max=180, use_circle=0, use_ellipse=1, min_distance=4.0,
                               lowpass_res=15, max_num_hits=400),
                     box=(256, 128), d2=dict(class2D_K=200, class2D_max_res=12, class2D_window_inner_A=200,
                                             class2D_recenter=1, class2D_force_max=0, class2D_min_dist_A=80,
                                             class2D_num_full_iter_batch=30,
                                             class2D_num_full_iter_batchsize_per_class=200)),
        "anchors": dict(pick=dict(diameter=30, diameter_max=50, use_circle=1, use_ellipse=0, min_distance=1.0,
                                  lowpass_res=15, max_num_hits=800),
                        box=(192, 96), d2=dict(class2D_K=200, class2D_max_res=12, class2D_window_inner_A=220,
                                               class2D_recenter=0, class2D_force_max=0,
                                               class2D_num_full_iter_batch=30,
                                               class2D_num_full_iter_batchsize_per_class=200)),
    }
    picks = {}
    for name, s in streams.items():
        picks[name] = R.run(cs, st, f"W2_pick_{name}", "W2", "blob_picker_gpu",
                            connections={"micrographs": (cur.uid, "exposures")}, params=s["pick"],
                            title=f"W2 blob pick: {name}", walltime="03:00:00", num_gpu=1)
    for name, s in streams.items():
        filt = R.save_subset(cs, st, f"W2_power_{name}", "W2", [(f"W2_pick_{name}", "particles")],
                             power_mask(exps, pf["keep_frac"][name]), "particle", "particles",
                             f"W2 {name} picks, top {pf['keep_frac'][name]:.0%} local power per defocus bin",
                             extra={"keep_frac": pf["keep_frac"][name]})
        if name == "anchors":   # census subset: random 25% of curated micrographs, fixed seed
            rng = np.random.default_rng(20261113)
            mics = np.array(sorted(set(exps["uid"].tolist())))
            keep = set(rng.choice(mics, size=len(mics) // 4, replace=False).tolist())
            filt = R.save_subset(cs, st, "W2_anchors_25pct", "W2", [("W2_power_anchors", "particles")],
                                 lambda ds: np.isin(ds["location/micrograph_uid"], list(keep)),
                                 "particle", "particles", "W2 anchors on a random 25% of curated micrographs",
                                 extra={"seed": 20261113, "micrographs": len(keep)})
        ext = R.run(cs, st, f"W2_extract_{name}", "W2", "extract_micrographs_multi",
                    connections={"micrographs": (cur.uid, "exposures"), "particles": (filt.uid, "particles")},
                    params=dict(box_size_pix=s["box"][0], bin_size_pix=s["box"][1], output_f16=1, compute_num_gpus=1),
                    title=f"W2 extract {name} {s['box'][0]} -> {s['box'][1]} px", walltime="06:00:00", num_gpu=1)
        st["jobs"][f"W2_extract_{name}"]["n_particles"] = ext.outputs["particles"].num_items
        R.save_state(st)
    jobs = []
    for name, s in streams.items():
        jobs.append(R.run(cs, st, f"W2_2D_{name}", "W2", "class_2D_new",
                          connections={"particles": (st["jobs"][f"W2_extract_{name}"]["uid"], "particles")},
                          params={**s["d2"], "compute_num_gpus": 1},
                          title=f"W2 2D round 1: {name} (K 200, 12 A)", walltime="11:00:00", num_gpu=1, block=False))
    for name, j in zip(streams, jobs):
        R.wait(cs, st["project"], j.uid, f"W2_2D_{name}", st)
    print("W2 round-1 2D complete", flush=True)


if __name__ == "__main__":
    main()
