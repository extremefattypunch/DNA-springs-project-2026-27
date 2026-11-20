"""A1 (gel-vs-EM panel, top discriminator): window-free length of DNA attached to sfGFP vs free DNA.

Round-1 rod 2D (P6 J27) used a 200 A circular window, which clips anything >100 A from the alignment centre; in
dot classes that centre sits on the dot, so attached DNA longer than ~85 A was truncated. Here each rod particle is
re-centred so a full 50 bp duplex (170 A) plus its dot fits inside a much larger window:
  dot_long_rod / dot_short_rod : new centre = dot position - 90 A along the +axis (toward the DNA)
  rod_no_dot                   : new centre = rod midpoint
Geometry from the co-localisation analysis ($SFGFP_DATA/rod_end_coloc/rods_analysed.npz: X, Y rod centre, ux, uy
+axis toward the dot, h half-length, all in micrograph px at 1.42 A; pidx = row of the J27 particle output; the dot
sits at centre + (h + 17 A) * u, script 17). Particles whose new box would cross the micrograph edge are dropped.
Then: extract 320 px @1.42 A -> F-crop 160 (2.84 A/px, 454 A box) with recenter_using_shifts OFF, and 2D with a
400 A window, no recentring, 150 classes, 12 A. Interpretation (results/gel_em_reconciliation.json): attached DNA
as long as free rods (~140-150 A, 4-5 helical repeats) -> linear adducts; a distinct shorter mode -> yielded loops.
  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/34_A1_window_free_rods.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

APIX_MIC = 1.42
BOX_PX = 320           # extraction box at 1.42 A/px (454 A)
GEO = os.path.join(R.DATA, "rod_end_coloc", "rods_analysed.npz")


def main():
    cs, st = R.session()
    key = "A1_recentred"
    if st["jobs"].get(key, {}).get("status") != "completed":
        g = np.load(GEO)
        src = cs.find_job(st["project"], st["jobs"]["W2_2D_rods"]["uid"]).load_output("particles")
        n_src = len(src)
        X, Y, ux, uy, h, lab, pidx = g["X"], g["Y"], g["ux"], g["uy"], g["h"], g["lab"], g["pidx"]
        dotcls = np.isin(lab, ["dot_long_rod", "dot_short_rod"])
        dot_x = X + ux * (h + 17 / APIX_MIC)
        dot_y = Y + uy * (h + 17 / APIX_MIC)
        cx = np.where(dotcls, dot_x - ux * (90 / APIX_MIC), X)
        cy = np.where(dotcls, dot_y - uy * (90 / APIX_MIC), Y)
        W = np.asarray(src["location/micrograph_shape"])[pidx, 1].astype(float)
        H = np.asarray(src["location/micrograph_shape"])[pidx, 0].astype(float)
        half = BOX_PX / 2
        ok = (cx > half) & (cx < W - half) & (cy > half) & (cy < H - half)
        print(f"{ok.sum()} of {len(ok)} rods keep a full {BOX_PX}px box inside the micrograph", flush=True)
        sub = src.take(pidx[ok]) if hasattr(src, "take") else src.mask(np.isin(np.arange(n_src), pidx[ok]))
        # Dataset.take keeps the order of pidx[ok]; guard against a mask-based fallback reordering rows
        if not hasattr(src, "take"):
            order = np.argsort(pidx[ok])
            cx, cy, W, H, lab = cx[ok][order], cy[ok][order], W[ok][order], H[ok][order], lab[ok][order]
        else:
            cx, cy, W, H, lab = cx[ok], cy[ok], W[ok], H[ok], lab[ok]
        sub["location/center_x_frac"] = (cx / W).astype(np.float32)
        sub["location/center_y_frac"] = (cy / H).astype(np.float32)
        juid = cs.find_project(st["project"]).save_external_result(
            st["workspaces"]["W4"], sub, "particle", name="particles",
            title="A1 rod particles re-centred for a window-free 2D (dot -90 A along axis / rod midpoint)")
        st["jobs"][key] = {"uid": juid, "type": "external", "ws": "W4", "status": "completed",
                           "rows_out": int(len(sub)), "by_label": {l: int((lab == l).sum()) for l in np.unique(lab)}}
        R.save_state(st)
        print(f"[{key}] {len(sub)} particles -> External {juid}", flush=True)
    rec = st["jobs"][key]["uid"]
    ext = R.run(cs, st, "A1_extract", "W4", "extract_micrographs_multi",
                connections={"micrographs": (st["jobs"]["W2_curated"]["uid"], "exposures"),
                             "particles": (rec, "particles")},
                params=dict(box_size_pix=BOX_PX, bin_size_pix=160, output_f16=1, recenter_using_shifts=0,
                            compute_num_gpus=1),
                title="A1 extract 320 px @1.42 -> 160 (454 A box), no shift recentring", walltime="04:00:00", num_gpu=1)
    R.run(cs, st, "A1_2D", "W4", "class_2D_new",
          connections={"particles": (ext.uid, "particles")},
          params=dict(class2D_K=150, class2D_max_res=12, class2D_window_inner_A=400, class2D_recenter=0,
                      class2D_force_max=0, class2D_num_full_iter_batch=30,
                      class2D_num_full_iter_batchsize_per_class=200, compute_num_gpus=1),
          title="A1 window-free 2D: 150 classes, 400 A window, no recentring", walltime="12:00:00", num_gpu=1)
    print("A1 complete", flush=True)


if __name__ == "__main__":
    main()
