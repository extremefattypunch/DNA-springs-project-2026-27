#!/bin/bash -l
# Plan section 5A - install EMReady v1.1 (He, Li & Huang 2023; GPL-3; the version the paper cites)
# into holylabs software/emready with its own venv.
#
# EMReady v1.1 ships as a 35 MB tarball (code + two trained models) and asks for PyTorch,
# einops, timm, mrcfile, tqdm, numpy plus an f2py-compiled Fortran interpolation module.
# The 2021-era pins (torch 1.8, python 3.8) are not honoured: a current torch cu12 wheel runs
# the same pure-PyTorch model; the smoke test below checks the import chain.
# CPU job: pip downloads (~2.5 GB torch wheel) and a Fortran compile only.  The venv lives on
# holylabs (software tier); the tarball is deleted after extraction (hekstra-cluster rule 3).
#SBATCH -J pp2a-install-emready
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-02:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.2.0-fasrc01          # gfortran for f2py
SW=/n/holylabs/hekstra_lab/Everyone/ianpoon/software
DEST=$SW/emready
PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
export TMPDIR=$TMP PIP_CACHE_DIR=$TMP/pipcache

if [ ! -f "$DEST/pred.py" ]; then
  cd "$TMP"; curl -sSL -o EMReady_v1.1.tgz http://huanglab.phys.hust.edu.cn/EMReady/EMReady_v1.1.tgz
  mkdir -p "$SW"; tar -xzf EMReady_v1.1.tgz && rm -f EMReady_v1.1.tgz
  rm -rf "$DEST"; mv EMReady "$DEST"
fi
if [ ! -x "$DEST/venv/bin/python" ]; then
  "$PY" -m venv "$DEST/venv"
  "$DEST/venv/bin/pip" install -q --upgrade pip
  "$DEST/venv/bin/pip" install -q torch --index-url https://download.pytorch.org/whl/cu124
  "$DEST/venv/bin/pip" install -q numpy einops timm mrcfile tqdm meson ninja
fi
cd "$DEST"
if ! ls interp3d*.so >/dev/null 2>&1; then
  "$DEST/venv/bin/python" -m numpy.f2py -c ./interp3d.f90 -m interp3d
fi
# Point EMReady.sh at this install (it expects a conda activator; we substitute the venv).
sed -i "s|^EMReady_home=.*|EMReady_home=\"$DEST\"|" EMReady.sh
sed -i "s|^activate=.*|activate=\"$DEST/venv/bin/activate\"|" EMReady.sh
sed -i "s|^EMReady_env=.*|EMReady_env=\"\"|" EMReady.sh
grep -nE "^(EMReady_home|activate|EMReady_env|source)" EMReady.sh | head
# smoke test: imports + model load on CPU
"$DEST/venv/bin/python" - <<'EOF'
import torch, einops, timm, mrcfile, numpy, interp3d
print("torch", torch.__version__, "cuda build", torch.version.cuda, "| einops", einops.__version__, "| timm", timm.__version__)
import importlib.util, sys
sys.path.insert(0, ".")
import scunet  # noqa: F401
sd = torch.load("model_state_dicts/model_grid_size_1.0.pth", map_location="cpu")
print("state dict tensors:", len(sd) if hasattr(sd, "__len__") else type(sd))
EOF
find "$DEST" -name "*.tgz" -delete
du -sh "$DEST"; echo "holylabs:"; quota /n/holylabs/hekstra_lab | tail -1
