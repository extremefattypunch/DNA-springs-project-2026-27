#!/usr/bin/env python
"""Record every tool version used, into software_versions.json (repo) - the reproducibility ledger."""
import json
import os
import re
import shlex
import subprocess
import sys

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
SW = os.environ.get("PP2A_SW", "/n/holylabs/hekstra_lab/Everyone/ianpoon/software")


def run(cmd, pattern=None, n=1):
    try:
        out = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, timeout=300)
        text = (out.stdout + out.stderr).strip()
        if pattern:
            m = re.search(pattern, text)
            return m.group(n) if m else text.splitlines()[0][:120] if text else None
        return text.splitlines()[0][:160] if text else None
    except Exception as e:  # noqa: BLE001
        return f"unavailable ({e.__class__.__name__})"


def main():
    py = os.environ.get("DNASPRING_PY", sys.executable)
    v = {"host_os": run("cat /etc/redhat-release"),
         "cryosparc_master": run("cryosparcm version"),
         "relion": run("relion_postprocess --version", r"RELION version:\s*(\S+)"),
         "gromacs": run("gmx --version", r"GROMACS version:\s*(\S+)"),
         "phenix": run("phenix.version", r"Version:\s*(\S+)"),
         "phenix_release": run("phenix.version", r"Release tag:\s*(\S+)"),
         "chimerax": run(shlex.split(os.environ.get("CHIMERAX", "false")) + ["--nogui", "--version"], r"ChimeraX version:\s*([^\n]+)") if os.environ.get("CHIMERAX") else None,
         "chimerax_container": os.path.basename(os.environ.get("CHIMERAX", "").split(".sif")[0] + ".sif") if os.environ.get("CHIMERAX") else None,
         "emready": "v1.1 (He, Li & Huang 2023; huanglab tarball EMReady_v1.1.tgz) + torch " + (run(f"{SW}/emready/venv/bin/python -c 'import torch;print(torch.__version__)'") or "?"),
         "namd": (run("namd3 +version") if subprocess.run("which namd3", shell=True, capture_output=True).returncode == 0
                  else "not installed (registration-gated download; OpenMM 8.6 used instead)"),
         "vmd": "not installed (registration download; ChimeraX/MDAnalysis used)",
         "python_dnaspring": run(f"{py} --version"),
         "cryosparc_tools": run(f"PYTHONPATH={os.environ.get('CS_PYLIBS', '')} {os.environ.get('CS_WORKER_PY', py)} -c 'import importlib.metadata as m;print(m.version(\"cryosparc-tools\"))'") if os.environ.get("CS_PYLIBS") else None,
         "alphafold_models": "AlphaFold DB v6 (AF2) AF-P30153-F1 / AF-P30154-F1 (AlphaFold Server AF3 not supplied)",
         }
    pkgs = ["openmm", "parmed", "mdtraj", "MDAnalysis", "mrcfile", "numpy", "scipy", "matplotlib", "biopython", "pymol"]
    v["dnaspring_packages"] = json.loads(run(f"""{py} -c "import importlib.metadata as m,json;print(json.dumps({{p:(lambda p:(m.version(p) if True else None))(p) for p in {pkgs!r} if True}}))" """) or "{}") if False else {}
    for p in pkgs:
        v["dnaspring_packages"][p] = run(f"{py} -c 'import importlib.metadata as m;print(m.version(\"{p}\"))'")
    v["ambertools"] = run(f"ls {os.path.dirname(os.path.dirname(py))}/conda-meta | grep -i '^ambertools'", r"ambertools-([\d.]+)")
    v["git_commit"] = run("git -C /n/home05/ianpoon/DNA-springs-project-2026-27 rev-parse --short HEAD")
    json.dump(v, open(os.path.join(REPO, "software_versions.json"), "w"), indent=2)
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    sys.exit(main())
