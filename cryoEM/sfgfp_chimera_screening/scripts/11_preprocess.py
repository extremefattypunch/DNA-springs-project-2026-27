"""W1 preprocessing of the whole session: import all 9,019 movies -> split in two -> patch motion x2 and
patch CTF x2, run concurrently on the two free gpu_test MIG slots.

Lane choice (deviation from the plan's "paid motion"): the pilot measured 3.7 s/movie for patch motion on a MIG
3g.20gb slice (P6 J3: 300 movies in 18.3 min), i.e. ~9.2 GPU-h for the session -> ~4.6 h per half on the two free
slots, which is no slower than queueing for a paid A100 with fairshare at 0.24. The paid budget is kept for
2D/3D.

Same acquisition parameters and gain handling as scripts/10_pilot.py (results/g0_acquisition.json).
  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/11_preprocess.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

RAW = os.environ.get("SFGFP_RAW", os.path.join(R.DATA, "raw"))
GAIN = os.path.join(RAW, "20260624_135522_EER_GainReference.gain")
ACQ = dict(psize_A=0.71, accel_kv=200.0, cs_mm=2.7, total_dose_e_per_A2=39.72)
EXPECTED = 9019


def main():
    cs, st = R.session()
    imp = R.run(cs, st, "W1_import", "W1", "import_movies",
                params=dict(blob_paths=os.path.join(RAW, "Images-Disc1", "GridSquare_*", "Data", "*_Fractions.tiff"),
                            gainref_path=GAIN, gainref_flip_x=0, gainref_flip_y=0, gainref_rotate_num=0,
                            parse_xml_files=1, skip_header_check=0, **ACQ,
                            xml_paths=os.path.join(RAW, "Images-Disc1", "GridSquare_*", "Data", "FoilHole_*.xml"),
                            mov_cut_suffix_xml=15, xml_cut_suffix_xml=4),
                title="W1 import all movies (gain, EPU xml beam shifts)", walltime="02:00:00")
    n = len(imp.load_output("imported_movies"))
    print(f"imported {n} movies (expected {EXPECTED})", flush=True)
    if n != EXPECTED:
        raise SystemExit("movie count does not match the facility inventory - check the download first")

    split = R.run(cs, st, "W1_split", "W1", "exposure_sets",
                  connections={"exposures_A": (imp.uid, "imported_movies")},
                  params=dict(set_operation="split", set_split_num=2),
                  title="W1 split movies into 2 halves (one per free MIG slot)", walltime="00:30:00")
    split.refresh()
    halves = sorted(h for h in split.outputs.keys() if "split" in h)[:2]
    print("split outputs:", halves, flush=True)

    mots = []
    for i, h in enumerate(halves):
        mots.append(R.run(cs, st, f"W1_motion_{i}", "W1", "patch_motion_correction_multi",
                          connections={"movies": (split.uid, h)},
                          params=dict(output_fcrop_factor=0.5, output_f16=1, compute_num_gpus=1),
                          title=f"W1 patch motion half {i} (F-crop 1/2 -> 1.42 A/px)",
                          walltime="09:00:00", num_gpu=1, block=False))
    for i, j in enumerate(mots):
        R.wait(cs, st["project"], j.uid, f"W1_motion_{i}", st)

    ctfs = []
    for i, j in enumerate(mots):
        ctfs.append(R.run(cs, st, f"W1_ctf_{i}", "W1", "patch_ctf_estimation_multi",
                          connections={"exposures": (j.uid, "micrographs")},
                          params=dict(compute_num_gpus=1), title=f"W1 patch CTF half {i}",
                          walltime="06:00:00", num_gpu=1, block=False))
    for i, j in enumerate(ctfs):
        R.wait(cs, st["project"], j.uid, f"W1_ctf_{i}", st)
    print("W1 motion + CTF complete", flush=True)


if __name__ == "__main__":
    main()
