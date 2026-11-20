#!/bin/bash -l
# Build the sfGFP + 50 bp closed-loop DNA chimera "L6_loop50" (fitting/validation model only).
# Minimal copy of the L3/L4 invocation in "MD simulations/sfGFP/build/build_loop_systems.sh":
# same sites (133/149), same 66 A span, same tethers -- only n_bp = 50.  Every output goes to
# cryoEM/sfgfp_chimera_screening/model/L6_loop50; nothing under "MD simulations/" is written.
#   usage: build_L6_loop50.sh [dna|assemble|system|all]   (default: dna + assemble)
set -euo pipefail
REPO=/n/home05/ianpoon/DNA-springs-project-2026-27
SF="$REPO/MD simulations/sfGFP"
OUT=$REPO/cryoEM/sfgfp_chimera_screening/model/L6_loop50
export DNASPRING_ENV=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring
export DNASPRING_PY=$DNASPRING_ENV/bin/python
export PATH="$DNASPRING_ENV/bin:$PATH" AMBERHOME="$DNASPRING_ENV" MPLBACKEND=Agg
export PYTHONDONTWRITEBYTECODE=1        # don't touch __pycache__ under MD simulations/
P=$DNASPRING_PY
B="$SF/build"; PARAMS="$B/02_params"
SPAN=66; NBP=50; SITES="133 149"
stage=${1:-default}
mkdir -p "$OUT"

if [[ $stage == dna || $stage == default || $stage == all ]]; then
  $P "$B/03_dna/build_dna.py" --n-bp $NBP --span-A $SPAN --outdir "$OUT/dna"
fi
if [[ $stage == assemble || $stage == default || $stage == all ]]; then
  $P "$B/04_assemble/assemble_chimera.py" --protein "$B/01_protein/sfgfp_prepped.pdb" \
     --dna "$OUT/dna/spring_${NBP}bp.pdb" --out "$OUT/chimera_loop50.pdb"
fi
if [[ $stage == system || $stage == all ]]; then
  bonds=$($P -c "import json;print(' '.join(json.load(open('$OUT/chimera_loop50.json'))['tleap_bonds']))")
  $P "$B/04_assemble/build_system.py" --name L6_loop50 --solute "$OUT/chimera_loop50.pdb" \
     --outdir "$OUT/tleap" --mg-M 0.005 \
     --frcmod "$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" "$PARAMS/DNP/DNP.frcmod" "$PARAMS/junction.frcmod" \
     --prep "$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin" "$PARAMS/DNP/DNP.prepin" \
     --bond $bonds --expect-residue CRO TDP DNL DNP --expect-loop TDP \
     --clamp-sites $SITES
fi
