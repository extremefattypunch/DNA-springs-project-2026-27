#!/bin/bash
# Plan section 5E - Phenix validation of the deposited models against the deposited maps.
#
# Reproduces the Table S1 refinement/validation rows without re-refining anything:
#   phenix.mtriage    -> half-map FSC (0.143), map-model FSC (0.5; "Model resolution"), d99, d_model
#   phenix.molprobity -> MolProbity score, clashscore, Ramachandran, rotamers, bond/angle RMSDs
# Phenix 2.2 (release 6143) comes from the labmate's install sourced in env/activate_pp2a.sh.
# Run per isoform:  bash scripts/24_phenix_validate.sh Aa   |   Ab
set -euo pipefail
iso=${1:?Aa or Ab}
case $iso in
  Aa) emd=48300; pdb=9MIP; res=3.4 ;;
  Ab) emd=48224; pdb=9MF5; res=3.2 ;;
  *) echo "iso must be Aa or Ab"; exit 1 ;;
esac
: "${PP2A_DATA:?source env/activate_pp2a.sh first}"
command -v phenix.mtriage >/dev/null || { echo "phenix not on PATH"; exit 1; }
E=$PP2A_DATA/em/emdb/$iso; M=$PP2A_DATA/em/models; OUT=$PP2A_DATA/em/phenix/$iso
mkdir -p "$OUT"; cd "$OUT"
phenix.version | head -3

echo "== mtriage ($iso, $pdb vs EMD-$emd) =="
phenix.mtriage "$E/emd_${emd}.map" "$E/emd_${emd}_half_map_1.map" "$E/emd_${emd}_half_map_2.map" \
  "$M/$pdb.cif" resolution=$res output.prefix=mtriage_$iso > mtriage_$iso.log 2>&1 || { tail -30 mtriage_$iso.log; exit 1; }
grep -E "Resolution|d99|d_model|d_fsc|FSC|Map-model|half" mtriage_$iso.log | head -40

echo "== molprobity ($pdb) =="
phenix.molprobity "$M/$pdb.cif" > molprobity_$iso.log 2>&1 || { tail -30 molprobity_$iso.log; exit 1; }
grep -E -A14 "Summary" molprobity.out | head -30 || grep -E "Ramachandran|Rotamer|Clashscore|MolProbity|RMS\(bonds\)|RMS\(angles\)" molprobity_$iso.log
echo "outputs in $OUT"
