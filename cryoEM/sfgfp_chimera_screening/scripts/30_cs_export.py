"""Export CryoSPARC outputs for figure scripts (which then run under $DNASPRING_PY without the master).

  PYTHONPATH=$CS_PYLIBS $CS_WORKER_PY scripts/30_cs_export.py KEY[:output] [KEY[:output] ...]

KEY is a cs_state.json job key (e.g. P0_ctf, P0_2D). For each requested output the whole dataset is written to
$SFGFP_DATA/exports/<KEY>__<output>.npz (one array per field; non-numeric fields as str) and its field list to
exports/field_manifest.json. For class_2D jobs the class-average stack (.mrc) is downloaded too; for refinement
jobs the FSC txt files and maps are copied by filename pattern. Nothing here guesses field names.
"""
import json
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs_runner as R  # noqa: E402

EXP = os.path.join(R.DATA, "exports")
COPY_PATTERNS = ("fsc", "map_sharp", "map_half", "volume_map", "mask_fsc", "map_locres", ".mrc", "_class_averages")


def dump(job, out, key):
    ds = job.load_output(out)
    arrays = {}
    for f in ds.fields():
        col = ds[f]
        a = np.asarray(col)
        if a.dtype == object:
            a = a.astype(str)
        arrays[f.replace("/", "__")] = a
    path = os.path.join(EXP, f"{key}__{out}.npz")
    np.savez_compressed(path, **arrays)
    return path, len(ds), list(ds.fields())


def main(specs):
    cs, st = R.session()
    puid = st["project"]
    os.makedirs(EXP, exist_ok=True)
    manpath = os.path.join(EXP, "field_manifest.json")
    manifest = json.load(open(manpath)) if os.path.exists(manpath) else {}
    for spec in specs:
        key, _, out = spec.partition(":")
        rec = st["jobs"][key]
        job = cs.find_job(puid, rec["uid"])
        outs = [out] if out else list(job.outputs.keys())
        for o in outs:
            try:
                path, n, fields = dump(job, o, key)
            except Exception as e:  # some outputs (e.g. volumes-only) have no rows
                print(f"[{key}:{o}] skipped: {e}")
                continue
            manifest[f"{key}:{o}"] = {"uid": rec["uid"], "rows": n, "fields": fields, "npz": path}
            print(f"[{key}:{o}] {n} rows, {len(fields)} fields -> {path}")
        # copy small result files (FSC curves, maps, class averages) from the job directory
        jdir = os.path.join(R.PROJECT_DIR, rec["uid"])
        if os.path.isdir(jdir):
            dest = os.path.join(EXP, key)
            os.makedirs(dest, exist_ok=True)
            for fn in os.listdir(jdir):
                if any(p in fn for p in COPY_PATTERNS) and os.path.getsize(os.path.join(jdir, fn)) < 2e9:
                    shutil.copy2(os.path.join(jdir, fn), os.path.join(dest, fn))
            print(f"[{key}] copied {len(os.listdir(dest))} result files -> {dest}")
    json.dump(manifest, open(manpath, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
