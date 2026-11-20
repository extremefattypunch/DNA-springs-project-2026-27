#!/bin/bash -l
# Install Topaz (topaz-em 0.3.x, Bepler et al. 2019) into a venv on holylabs for CryoSPARC's
# Topaz Train / Topaz Extract wrappers, plus software/bin/topaz-cs, a launcher that strips the
# CryoSPARC worker's PYTHONPATH/PYTHONHOME/LD_LIBRARY_PATH so Topaz runs on its own torch.
# CPU-only job (pip downloads); the GPU smoke test is slurm/check_topaz_gpu.sbatch.
#SBATCH -J sfgfp-install-topaz
#SBATCH -p test
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-01:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
SW=/n/holylabs/hekstra_lab/Everyone/ianpoon/software
DEST=$SW/topaz-0.3
PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring/bin/python
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
export TMPDIR=$TMP PIP_CACHE_DIR=$TMP/pipcache
quota /n/holylabs/hekstra_lab | tail -1
if [ ! -x "$DEST/bin/topaz" ]; then
  "$PY" -m venv "$DEST"
  "$DEST/bin/pip" install -q --upgrade pip
  "$DEST/bin/pip" install -q torch torchvision --index-url https://download.pytorch.org/whl/cu124
  "$DEST/bin/pip" install -q "topaz-em>=0.3,<0.4"
fi
mkdir -p "$SW/bin"
cat > "$SW/bin/topaz-cs" <<'WRAP'
#!/bin/bash
# Topaz launcher for CryoSPARC (exec_path). Runs the holylabs venv's topaz with a clean Python env.
exec env -u PYTHONPATH -u PYTHONHOME -u LD_LIBRARY_PATH PYTHONNOUSERSITE=1 \
  /n/holylabs/hekstra_lab/Everyone/ianpoon/software/topaz-0.3/bin/topaz "$@"
WRAP
chmod 755 "$SW/bin/topaz-cs"
"$SW/bin/topaz-cs" --version
"$DEST/bin/python" -c "import torch, topaz; print('torch', torch.__version__, 'cuda build', torch.version.cuda, 'topaz', topaz.__version__)"
du -sh "$DEST"; quota /n/holylabs/hekstra_lab | tail -1
