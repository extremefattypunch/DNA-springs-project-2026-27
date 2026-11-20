"""P0 pilot (free lane): one grid square through import -> patch motion -> patch CTF -> blob pick ->
extract -> 2D, to measure per-movie cost on a MIG slice and check Gate G1 before any paid GPU time.

Acquisition (G0, results/g0_acquisition.json): Glacios 200 kV, Falcon 4i (energy filter), 4096^2 uint8 TIFF,
50 fractions, 0.71 A/px, 39.72 e/A^2 (Ian's vacuum calibration; 28.2 e/A^2 detected through ice), raw counts
(not gain-corrected) -> EPU EER gain reference applied with NO flip/rotation (flat-field correlation test).

  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/10_pilot.py [GridSquare_xxx]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

RAW = os.environ.get("SFGFP_RAW", os.path.join(R.DATA, "raw"))
SQUARE = sys.argv[1] if len(sys.argv) > 1 else "GridSquare_12363290"
GAIN = os.path.join(RAW, "20260624_135522_EER_GainReference.gain")
ACQ = dict(psize_A=0.71, accel_kv=200.0, cs_mm=2.7, total_dose_e_per_A2=39.72)
N_PILOT = 300


def main():
    cs, st = R.session()
    imp = R.run(cs, st, "P0_import", "P0", "import_movies",
                params=dict(blob_paths=os.path.join(RAW, "Images-Disc1", SQUARE, "Data", "*_Fractions.tiff"),
                            gainref_path=GAIN, gainref_flip_x=0, gainref_flip_y=0, gainref_rotate_num=0,
                            parse_xml_files=1, skip_header_check=0, **ACQ,
                            # EPU xml is <stem>.xml, the movie <stem>_Fractions.tiff: trim 15 chars to match
                            xml_paths=os.path.join(RAW, "Images-Disc1", SQUARE, "Data", "FoilHole_*.xml"),
                            mov_cut_suffix_xml=15, xml_cut_suffix_xml=4),
                title=f"P0 import {SQUARE}", walltime="01:00:00")
    mot = R.run(cs, st, "P0_motion", "P0", "patch_motion_correction_multi",
                connections={"movies": (imp.uid, "imported_movies")},
                params=dict(random_num=N_PILOT, output_fcrop_factor=0.5, output_f16=1, compute_num_gpus=1),
                title=f"P0 patch motion ({N_PILOT} movies, F-crop 1/2 -> 1.42 A/px)", walltime="04:00:00", num_gpu=1)
    ctf = R.run(cs, st, "P0_ctf", "P0", "patch_ctf_estimation_multi",
                connections={"exposures": (mot.uid, "micrographs")},
                params=dict(compute_num_gpus=1), title="P0 patch CTF", walltime="02:00:00", num_gpu=1)
    pick = R.run(cs, st, "P0_blob", "P0", "blob_picker_gpu",
                 connections={"micrographs": (ctf.uid, "exposures")},
                 params=dict(diameter=80, diameter_max=140, use_circle=1, use_ellipse=1),
                 title="P0 blob pick 80-140 A (circle + ellipse)", walltime="02:00:00", num_gpu=1)
    ext = R.run(cs, st, "P0_extract", "P0", "extract_micrographs_multi",
                connections={"micrographs": (pick.uid, "micrographs"), "particles": (pick.uid, "particles")},
                params=dict(box_size_pix=256, bin_size_pix=128, output_f16=1, compute_num_gpus=1),
                title="P0 extract 256 @1.42 A -> F-crop 128 (2.84 A/px)", walltime="02:00:00", num_gpu=1)
    R.run(cs, st, "P0_2D", "P0", "class_2D_new",
          connections={"particles": (ext.uid, "particles")},
          params=dict(class2D_K=100, class2D_max_res=12, class2D_window_inner_A=200, class2D_force_max=0,
                      class2D_num_full_iter_batch=40, class2D_num_full_iter_batchsize_per_class=200,
                      compute_num_gpus=1),
          title="P0 2D (100 classes, 12 A, marginalised, batch 200/class)", walltime="04:00:00", num_gpu=1)
    print("pilot complete", {k: v.get("run_s") for k, v in st["jobs"].items() if k.startswith("P0_")})


if __name__ == "__main__":
    main()
