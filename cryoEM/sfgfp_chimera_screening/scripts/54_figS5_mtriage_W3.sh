#!/bin/bash
# Fig S5C map-model FSC for the W3 DNA-control map: phenix.mtriage of the rigidly fitted straight 32 bp B-DNA
# (scripts/53_figS5_fit_W3_dna.py) and of its register variants (scripts/53b_figS5_register_variants.py: spin 90,
# spin 180, polarity flip) against the RELION half-map average and the J40 half maps. The variants test whether the
# map-model FSC depends on the helical register, which the fit leaves undetermined.
# Phenix is the labmate's install, used read-only. Outputs go to netscratch (never /n/lab_storage).
#   sbatch -p test -c 4 --mem 16G -t 1:00:00 --wrap "bash scripts/54_figS5_mtriage_W3.sh"
set -eo pipefail
source /n/lab_storage/hekstra_lab/people/aldama/software/phenix2/phenix_env.sh
set -u
R=/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/relion/W3_duplex
F=/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/figures_scratch/W3_duplex/fit
H1=$(readlink -f $R/half1.mrc); H2=$(readlink -f $R/half2.mrc)
# resolution = RELION loose-mask FSC 0.143 (8.66 A, relion/W3_duplex/pp_loose/postprocess.star)
for v in best spin090 spin180 flip; do
  pdb=$F/dna_32bp_fit.pdb; [ "$v" = best ] || pdb=$F/dna_32bp_fit_$v.pdb
  OUT=$F/mtriage_$v; mkdir -p "$OUT"; cd "$OUT"
  phenix.mtriage "$R/raw.mrc" "$H1" "$H2" "$pdb" resolution=8.66 > mtriage.log 2>&1 || { echo "== $v FAILED"; tail -40 mtriage.log; exit 1; }
  echo "== $v"; grep -E "FSC\(map,model map\)=0\.5|FSC\(half map 1,2\)=0\.143|CC_mask|d_model " mtriage.log | head -8
done
