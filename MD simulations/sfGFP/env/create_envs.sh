#!/bin/bash
# Build both conda envs. Sequential on purpose: micromamba shares one package
# cache, and two concurrent solves race on it.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/activate.sh"
MM="$MAMBA_EXE"

echo "=== $(date -Is) creating $DNASPRING_ENV ==="
"$MM" create -y -p "$DNASPRING_ENV" -f "$HERE/environment.yml"

echo "=== $(date -Is) creating $DNASPRING_QM_ENV ==="
"$MM" create -y -p "$DNASPRING_QM_ENV" -f "$HERE/environment-qm.yml" || {
    echo "WARN: QM env solve failed; RESP/torsion-scan step will need a fallback." >&2
}

echo "=== $(date -Is) done ==="
"$DNASPRING_PY" -c "import openmm; print('openmm', openmm.__version__)"
