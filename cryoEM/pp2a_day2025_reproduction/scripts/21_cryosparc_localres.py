#!/usr/bin/env python
"""Plan section 5B - CryoSPARC on the deposited half maps (watched live in the web GUI).

Per isoform, in project CS-pp2a-day2025 / workspace "deposited maps":
  Import 3D Volumes (half A, half B, and the deposited sharpened map)
    -> Local Resolution Estimation  (FSC 0.5, adaptive window; GPU on lane Cannon-test)   [Fig S5A/E cross-check]
    -> Local Filtering              (locally filtered map for rendering)
    -> Validation (FSC)             (gold-standard FSC with an auto mask)                  [Fig S5C/G cross-check]
Job-type strings and parameter names were read from cs.job_register on this v5.0.7 master
(import_volumes, local_resolution, local_filter, validation).

Runs with CryoSPARC's worker python + cryosparc-tools:
  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/21_cryosparc_localres.py
State (project/workspace/job UIDs) is kept in $PP2A_DATA/cs_state/cs_state.json so re-runs
resume instead of duplicating jobs.  All jobs are created on the free lane; the master must
already be up (launchcryosparc) and the GUI open.
"""
import json
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cs_connect import connect  # noqa: E402

DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
PROJECT_PARENT = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects"
PROJECT_TITLE = "CS-pp2a-day2025"
WORKSPACE_TITLE = "deposited maps"
LANE = "Cannon-test"            # gpu_test MIG slice: free.  Never 'Cannon' (836.5 billing / GPU-h).
STATE = os.path.join(DATA, "cs_state", "cs_state.json")
ISO = {"Aa": {"emd": "48300", "psize": 1.056}, "Ab": {"emd": "48224", "psize": 1.058}}


def load_state():
    return json.load(open(STATE)) if os.path.exists(STATE) else {"jobs": {}, "timings": {}}


def save_state(s):
    os.makedirs(os.path.dirname(STATE), exist_ok=True)
    json.dump(s, open(STATE, "w"), indent=2)


def wait(job, key, state, t0):
    last = None
    while True:
        job.refresh()
        st = job.status
        if st != last:
            print(f"[{key}] {job.uid} {st} (+{time.time() - t0:.0f}s)", flush=True)
            last = st
        if st in ("completed", "failed", "killed"):
            break
        time.sleep(15)
    state["timings"][key] = round(time.time() - t0, 1)
    if st != "completed":
        try:
            print("\n".join(str(job.log()).splitlines()[-30:]))
        except Exception:
            pass
        save_state(state)
        raise SystemExit(f"{key} ({job.uid}) ended {st}")
    state["jobs"][key] = job.uid
    save_state(state)
    return job


def run_job(cs, state, key, puid, wuid, jtype, connections, params, title):
    if key in state["jobs"]:
        print(f"[{key}] already done: {state['jobs'][key]}")
        return cs.find_job(puid, state["jobs"][key])
    job = cs.create_job(puid, wuid, jtype, connections=connections, params=params, title=title)
    print(f"[{key}] created {job.uid} ({jtype}) -> queue on {LANE}", flush=True)
    t0 = time.time()
    job.queue(lane=LANE)
    return wait(job, key, state, t0)


def main():
    cs, c = connect()
    state = load_state()
    # project / workspace
    if "project" not in state:
        os.makedirs(PROJECT_PARENT, exist_ok=True)
        proj = cs.create_project(PROJECT_TITLE, PROJECT_PARENT,
                                 desc="Day et al. 2025 PP2A Aa/Ab-B56g-Ca: analysis of deposited EMDB maps (no particles).")
        state["project"] = proj.uid
        save_state(state)
    puid = state["project"]
    proj = cs.find_project(puid)
    if "workspace" not in state:
        ws = cs.create_workspace(puid, WORKSPACE_TITLE, desc="Local resolution / FSC on EMD-48300 and EMD-48224 half maps")
        state["workspace"] = ws.uid
        save_state(state)
    wuid = state["workspace"]
    print(f"project {puid} ({proj.dir()}), workspace {wuid}")

    for iso, meta in ISO.items():
        emd, psize = meta["emd"], meta["psize"]
        src = os.path.join(DATA, "em", "emdb", iso)
        # import: the 'volume_out_name' tells CryoSPARC which slot the blob fills
        imports = {}
        for slot, fname in (("map_half_A", f"emd_{emd}_half_map_1.map"),
                            ("map_half_B", f"emd_{emd}_half_map_2.map"),
                            ("map", f"emd_{emd}.map")):
            key = f"{iso}_import_{slot}"
            job = run_job(cs, state, key, puid, wuid, "import_volumes", connections={},
                          params={"volume_blob_path": os.path.join(src, fname), "volume_psize": psize,
                                  "volume_out_name": slot},
                          title=f"{iso} import EMD-{emd} {slot}")
            imports[slot] = job
        # local resolution from the two half maps
        conn_halves = {"volume": [(imports["map_half_A"].uid, "imported_volume_1"),
                                  (imports["map_half_B"].uid, "imported_volume_1")]}
        locres = run_job(cs, state, f"{iso}_locres", puid, wuid, "local_resolution",
                         connections=conn_halves,
                         params={"locres_fsc_thresh": 0.5, "locres_use_mask": 0, "locres_compute_facility": "GPU"},
                         title=f"{iso} local resolution (EMD-{emd} half maps, FSC 0.5)")
        # locally filtered map: half maps + local-res map
        run_job(cs, state, f"{iso}_locfilt", puid, wuid, "local_filter",
                connections={"volume": [(imports["map_half_A"].uid, "imported_volume_1"),
                                        (imports["map_half_B"].uid, "imported_volume_1"),
                                        (locres.uid, "volume")]},
                params={"localfilter_use_mask": 0},
                title=f"{iso} local filtering")
        # gold-standard FSC with an auto-generated mask
        run_job(cs, state, f"{iso}_fsc", puid, wuid, "validation",
                connections=conn_halves,
                params={"validate_generate_new_mask": 1, "validate_compute_facility": "GPU"},
                title=f"{iso} validation (FSC) EMD-{emd}")

    # copy the small outputs (local-res map, FSC txt/plots) next to the other EM results
    out = os.path.join(DATA, "em", "cryosparc")
    os.makedirs(out, exist_ok=True)
    for key, uid in state["jobs"].items():
        if key.endswith(("_locres", "_locfilt", "_fsc")):
            jdir = os.path.join(proj.dir(), uid)
            for f in os.listdir(jdir):
                if f.endswith((".mrc", ".txt", ".png", ".json")) and os.path.getsize(os.path.join(jdir, f)) < 300e6:
                    shutil.copy(os.path.join(jdir, f), os.path.join(out, f"{key}_{f}"))
    print("done; outputs mirrored to", out)
    print(json.dumps(state, indent=1))


if __name__ == "__main__":
    sys.exit(main())
