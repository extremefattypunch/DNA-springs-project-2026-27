"""Shared CryoSPARC driver for the sfGFP chimera screening project.

Extends the PP2A reproduction's run_job (pp2a_day2025_reproduction/scripts/21_cryosparc_localres.py) with:
  * every job type and param key validated against the live v5 job register before create_job
    (a misspelt key would otherwise silently fall back to the default)
  * per-job lane + walltime (lanes take a `walltime` cluster var since 2026-11-13)
  * a paid-GPU ledger: jobs on lane `Cannon` (A100, 836.5 fairshare units / GPU-h) are refused once the
    projected total would pass PAID_GPU_H_CAP (Ian's ceiling: ~40 GPU-h, then ask)
  * resume: the job UID is written to state the moment it is created, so a crashed/restarted driver
    re-attaches to the running job instead of creating a duplicate
  * refusal to queue a job whose walltime would outlive the CryoSPARC master's Slurm job

Run under CryoSPARC's python:  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/<stage>.py
"""
import datetime as dt
import fcntl
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cs_connect import connect  # noqa: E402

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")
STATE = os.path.join(DATA, "cs_state", "cs_state.json")
STATE_MIRROR = os.path.join(REPO, "results", "cs_state.json")
PROJECT_TITLE = "CS-sfgfp-chimera-screening"
PROJECT_PARENT = "/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects"
PROJECT_DIR = os.path.join(PROJECT_PARENT, "CS-cs-sfgfp-chimera-screening")  # CryoSPARC prefixed "CS-"
PAID_LANE, FREE_LANE = "Cannon", "Cannon-test"
PAID_GPU_H_CAP = float(os.environ.get("SFGFP_PAID_GPU_H_CAP", "40"))
BILLING_PER_GPU_H = 836.5
WORKSPACES = {
    "W1": "W1 preprocess (motion, CTF, curation)",
    "W2": "W2 picking + 2D classification",
    "W3": "W3 ab initio + decoys + hetero/NU cycles",
    "W4": "W4 3D classification (PCA) + local refinement",
    "W5": "W5 Topaz re-pick",
    "W6": "W6 final refinement + validation",
    "W7": "W7 heterogeneity (Fig S7, 3DVA)",
    "P0": "P0 pilot (first grid squares, free lane)",
}


# ---------------------------------------------------------------- state --
def load_state():
    if os.path.exists(STATE):
        with open(STATE) as fh:
            return json.load(fh)
    return {"project": None, "workspaces": {}, "jobs": {}, "paid_gpu_h": 0.0}


def save_state(state):
    """Merge-save: several drivers (e.g. pilot + preprocess) may run at once, so re-read the file under a lock,
    overlay this process's job records (newer 'finished'/status wins) and write atomically."""
    os.makedirs(os.path.dirname(STATE), exist_ok=True)
    with open(STATE + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        disk = load_state()
        for k, v in disk.get("jobs", {}).items():
            mine = state["jobs"].get(k)
            if mine is None or (mine.get("uid") == v.get("uid") and v.get("status") == "completed"
                                and mine.get("status") != "completed"):
                state["jobs"][k] = v
        for k, v in disk.get("workspaces", {}).items():
            state["workspaces"].setdefault(k, v)
        state["project"] = state.get("project") or disk.get("project")
        state["paid_gpu_h"] = round(sum(j.get("paid_gpu_h", 0.0) for j in state["jobs"].values()), 3)
        tmp = STATE + f".tmp{os.getpid()}"
        with open(tmp, "w") as fh:
            json.dump(state, fh, indent=1, default=str)
        os.replace(tmp, STATE)
    os.makedirs(os.path.dirname(STATE_MIRROR), exist_ok=True)
    with open(STATE_MIRROR, "w") as fh:
        json.dump(state, fh, indent=1, default=str)


# --------------------------------------------------------------- setup --
def setup(cs, state):
    if not state.get("project"):
        os.makedirs(PROJECT_PARENT, exist_ok=True)
        p = cs.api.projects.create(
            title=PROJECT_TITLE, parent_dir=PROJECT_PARENT,
            description="sfGFP + 50 bp DNA spring chimera, 12-11-2026 Glacios/Falcon 4i screening session. "
                        "Driven by cryoEM/sfgfp_chimera_screening/scripts (Day et al. 2025 Fig S6 workflow).")
        state["project"] = p.uid
        save_state(state)
        print(f"created project {p.uid} at {p.project_dir}", flush=True)
    puid = state["project"]
    for key, title in WORKSPACES.items():
        if key not in state["workspaces"]:
            w = cs.create_workspace(puid, title)
            state["workspaces"][key] = w.uid
            save_state(state)
            print(f"created workspace {key} -> {w.uid}", flush=True)
    return puid


# ---------------------------------------------------------- validation --
_REG = {}


def register(cs):
    if not _REG:
        for s in cs.job_register.specs:
            d = s.model_dump()
            _REG[d["type"]] = d
    return _REG


def validate(cs, jtype, params, connections):
    reg = register(cs)
    if jtype not in reg:
        raise SystemExit(f"unknown job type {jtype!r} on this master")
    spec = reg[jtype]
    bad = [k for k in params if k not in (spec.get("params") or {})]
    if bad:
        raise SystemExit(f"{jtype}: unknown param(s) {bad}; valid: {sorted(spec['params'])}")
    badin = [k for k in connections if k not in (spec.get("inputs") or {})]
    if badin:
        raise SystemExit(f"{jtype}: unknown input(s) {badin}; valid: {sorted(spec['inputs'])}")


# ------------------------------------------------------------ guards --
def master_end():
    out = subprocess.run(["squeue", "-u", os.environ.get("USER", "ianpoon"), "-h", "-n", "CryoSPARC", "-o", "%e"],
                         capture_output=True, text=True).stdout.split()
    try:
        return dt.datetime.fromisoformat(out[0])
    except (IndexError, ValueError):
        return None


def _hours(walltime):
    d, _, hms = walltime.rpartition("-")
    h, m, s = (list(map(int, hms.split(":"))) + [0, 0])[:3]
    return (int(d) if d else 0) * 24 + h + m / 60 + s / 3600


def guard(state, lane, walltime, est_gpu_h, key):
    end = master_end()
    if end is not None:
        left = (end - dt.datetime.now()).total_seconds() / 3600
        if _hours(walltime) > left:
            raise SystemExit(f"[{key}] walltime {walltime} outlives the CryoSPARC master ({left:.1f} h left); "
                             "relaunch the master (with a DB backup) first")
    if lane == PAID_LANE:
        proj = state.get("paid_gpu_h", 0.0) + est_gpu_h
        print(f"[{key}] PAID lane: est {est_gpu_h:.1f} GPU-h = {est_gpu_h * BILLING_PER_GPU_H:,.0f} units; "
              f"ledger {state.get('paid_gpu_h', 0.0):.1f} -> {proj:.1f} / cap {PAID_GPU_H_CAP:.0f} GPU-h", flush=True)
        if proj > PAID_GPU_H_CAP:
            raise SystemExit(f"[{key}] would exceed the paid cap ({proj:.1f} > {PAID_GPU_H_CAP} GPU-h) - ask Ian")


GPU_TEST_QOS_MAX = 2   # gpu_test QOS: 2 jobs per user, submitted OR running; a 3rd is rejected at sbatch


def wait_gpu_test_slot(key, poll=30):
    """Block until fewer than GPU_TEST_QOS_MAX of Ian's jobs (any project) sit in gpu_test."""
    told = False
    while True:
        out = subprocess.run(["squeue", "-u", os.environ.get("USER", "ianpoon"), "-h", "-p", "gpu_test", "-o", "%i"],
                             capture_output=True, text=True).stdout.split()
        if len(out) < GPU_TEST_QOS_MAX:
            return
        if not told:
            print(f"[{key}] waiting for a gpu_test slot ({len(out)}/{GPU_TEST_QOS_MAX} in use)", flush=True)
            told = True
        time.sleep(poll)


# ------------------------------------------------------------- running --
def _log_params(key, jtype, lane, walltime, connections, params):
    d = os.path.join(REPO, "results", "job_params")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{key}.json"), "w") as fh:
        json.dump({"type": jtype, "lane": lane, "walltime": walltime,
                   "connections": connections, "params": params}, fh, indent=1, default=str)


def wait(cs, puid, juid, key, state, poll=20):
    job = cs.find_job(puid, juid)
    t0, last = time.time(), None
    while True:
        job.refresh()
        st = job.status
        if st != last:
            print(f"[{key}] {puid} {juid} {st} (+{time.time() - t0:.0f}s)", flush=True)
            last = st
        if st in ("completed", "failed", "killed"):
            break
        time.sleep(poll)
    rec = state["jobs"][key]
    rec["status"] = st
    rec["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    try:
        m = job.model
        if getattr(m, "started_at", None) and getattr(m, "completed_at", None):
            rec["run_s"] = (m.completed_at - m.started_at).total_seconds()
    except Exception:
        pass
    if rec.get("lane") == PAID_LANE and rec.get("num_gpu"):
        used = rec.get("run_s", time.time() - t0) / 3600 * rec["num_gpu"]
        state["paid_gpu_h"] = round(state.get("paid_gpu_h", 0.0) + used, 3)
        rec["paid_gpu_h"] = round(used, 3)
    save_state(state)
    if st != "completed":
        logf = os.path.join(PROJECT_DIR, juid, "job.log")
        if os.path.exists(logf):
            with open(logf, errors="replace") as fh:
                print("".join(fh.readlines()[-40:]), flush=True)
        raise SystemExit(f"[{key}] {juid} ended {st}")
    return job


def run(cs, state, key, ws, jtype, connections=None, params=None, title="", lane=FREE_LANE,
        walltime="08:00:00", est_gpu_h=0.0, num_gpu=None, extra_results=(), queue=True, block=True):
    """Create (or re-attach to) one job, queue it on `lane`, and optionally wait for it.

    connections: {input_name: (job_uid, output_name)} or a list of such tuples (multi-connect).
    extra_results: (input_name, idx, result_name, src_job, src_output, src_result) slot overrides.
    """
    connections, params = connections or {}, params or {}
    puid, wuid = state["project"], state["workspaces"][ws]
    rec = state["jobs"].get(key)
    if rec and rec.get("status") == "completed":
        print(f"[{key}] already completed: {rec['uid']}", flush=True)
        return cs.find_job(puid, rec["uid"])
    if rec and rec.get("uid"):
        job = cs.find_job(puid, rec["uid"])
        job.refresh()
        if job.status in ("queued", "launched", "started", "running", "waiting"):
            print(f"[{key}] re-attaching to {rec['uid']} ({job.status})", flush=True)
            return wait(cs, puid, rec["uid"], key, state) if block else job
        if job.status == "completed":
            return wait(cs, puid, rec["uid"], key, state)
        if job.status == "building" and queue:
            # created by a driver that died before queueing it: queue that job rather than orphaning it
            print(f"[{key}] {rec['uid']} was created but never queued; queueing it now", flush=True)
            guard(state, lane, walltime, est_gpu_h, key)
            if lane == FREE_LANE and (num_gpu or 0) > 0:
                wait_gpu_test_slot(key)
            job.queue(lane=lane, cluster_vars={"walltime": walltime})
            rec["status"] = "queued"
            save_state(state)
            return wait(cs, puid, rec["uid"], key, state) if block else job
        print(f"[{key}] previous attempt {rec['uid']} is {job.status}; creating a fresh job", flush=True)
    validate(cs, jtype, params, connections)
    guard(state, lane, walltime, est_gpu_h, key)
    conn_args = {}
    multi = []
    for inp, src in connections.items():
        if isinstance(src, list):
            conn_args[inp] = src[0]
            multi += [(inp, s) for s in src[1:]]
        else:
            conn_args[inp] = src
    job = cs.create_job(puid, wuid, jtype, connections=conn_args, params=params, title=title or key)
    for inp, (sj, so) in multi:
        job.connect(inp, sj, so)
    for inp, idx, res, sj, so, sr in extra_results:
        cs.api.jobs.connect_result(puid, job.uid, inp, idx, res,
                                   source_output_name=so, source_result_name=sr, source_job_uid=sj)
    state["jobs"][key] = {"uid": job.uid, "type": jtype, "ws": ws, "lane": lane, "walltime": walltime,
                          "num_gpu": num_gpu, "status": "created", "title": title or key,
                          "parents": sorted({s[0] for s in conn_args.values()} | {s[1][0] for s in multi}),
                          "created": dt.datetime.now().isoformat(timespec="seconds")}
    save_state(state)
    _log_params(key, jtype, lane, walltime, connections, params)
    if not queue:
        print(f"[{key}] created {puid} {job.uid} ({jtype}) - not queued (interactive/manual)", flush=True)
        return job
    if lane == FREE_LANE and (num_gpu or 0) > 0:
        wait_gpu_test_slot(key)
    job.queue(lane=lane, cluster_vars={"walltime": walltime})
    if lane == FREE_LANE and (num_gpu or 0) > 0:
        # CryoSPARC calls sbatch a few seconds after queue(); wait until it has, so the next slot check sees it
        for _ in range(60):
            job.refresh()
            if job.status != "queued":
                break
            time.sleep(5)
    state["jobs"][key]["status"] = "queued"
    save_state(state)
    print(f"[{key}] created + queued {puid} {wuid} {job.uid} ({jtype}) on {lane}, walltime {walltime}", flush=True)
    return wait(cs, puid, job.uid, key, state) if block else job


def session():
    cs, cfg = connect()
    state = load_state()
    setup(cs, state)
    return cs, state
