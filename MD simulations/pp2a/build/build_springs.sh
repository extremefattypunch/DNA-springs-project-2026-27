#!/bin/bash
# Build the ss and ds chimeras and the bare springs for one DNA length:  build_springs.sh N [SEQUENCE]
#
# The experiment Zocchi's model is about is the pair: the same strand clicked at both ends,
# unhybridised (an entropic coil that at most pulls) and hybridised (a leaf spring that
# pushes).  The sfGFP/ubiquitin work ran that pair at 50 nt; here it is run at 50 and 80 so
# the length dependence of every metric can be read, as the sfGFP page did for 27/40 bp.
#   P3_ss50 / P4_ds50     50 nt: the centre of Zocchi's verified 60-mer
#   P5_ss80 / P6_ds80     80 nt: that 60-mer padded by 10 nt each side (chosen for the
#                         shortest self-complementary word; a reuse plus padding, not a
#                         designed oligo)
#   D_ss<N> / D_ds<N>     the bare strand / duplex for the restraint-window force measurement
# Same 90 A anchor span throughout, so length and hybridisation are the only variables.
set -euo pipefail
cd "$(dirname "$0")/.."
source ./env/activate.sh
N=${1:?length in nt}; SEQ=${2:-}
case $N in 50) SS=P3_ss50; DS=P4_ds50 ;; 80) SS=P5_ss80; DS=P6_ds80 ;; *) SS=P_ss$N; DS=P_ds$N ;; esac
P="$DNASPRING_PY"; SF="$SFGFP_REPO/build"; PARAMS=build/02_params
SITES="7 553"; SPAN=${SPAN:-90}
frc=("$PARAMS/TET/TET.frcmod" "$PARAMS/TDP/TDP.frcmod" "$PARAMS/DNL/DNL.frcmod" \
     "$PARAMS/DNP/DNP.frcmod" "$PARAMS/DNH/DNH.frcmod" "$PARAMS/junction.frcmod")
prep=("$PARAMS/TET/TET.prepin" "$PARAMS/TDP/TDP.prepin" "$PARAMS/DNL/DNL.prepin" \
      "$PARAMS/DNP/DNP.prepin" "$PARAMS/DNH/DNH.prepin")
bonds_of () { $P -c "import json;print(' '.join(json.load(open('$1'))['tleap_bonds']))"; }
seqarg=(); [ -n "$SEQ" ] && seqarg=(--sequence "$SEQ")

echo "##### springs, $N nt/bp, bent to $SPAN A"
$P "$SF/03_dna/build_dna.py" --n-bp $N --span-A $SPAN --single-strand "${seqarg[@]}" --outdir build/03_dna
$P "$SF/03_dna/build_dna.py" --n-bp $N --span-A $SPAN "${seqarg[@]}" --outdir build/03_dna

build_chimera () {   # name, spring pdb
  local name=$1 spring=$2 pdb=build/04_assemble/pp2a_${1#P?_}.pdb
  $P "$SF/04_assemble/assemble_chimera.py" --protein build/01_protein/pp2a_prepped.pdb \
     --dna "$spring" --sites $SITES --out "$pdb"
  $P "$SF/04_assemble/build_system.py" --name "$name" --solute "$pdb" --outdir "build/systems/$name" \
     --mg-M 0.005 --frcmod "${frc[@]}" --prep "${prep[@]}" --bond $(bonds_of "${pdb%.pdb}.json") \
     --expect-residue TDP DNL DNP --expect-loop TDP --clamp-sites $SITES
}
echo "##### $SS"; build_chimera $SS build/03_dna/spring_${N}nt.pdb
echo "##### $DS"; build_chimera $DS build/03_dna/spring_${N}bp.pdb
echo "##### bare springs"
$P build/03_dna/bare_strand.py --spring build/03_dna/spring_${N}nt.pdb --name D_ss$N --pad "${BARE_PAD:-24}"
$P build/03_dna/bare_strand.py --spring build/03_dna/spring_${N}bp.pdb --name D_ds$N --pad "${BARE_PAD:-24}"
for d in build/systems/$SS build/systems/$DS build/systems/D_ss$N build/systems/D_ds$N; do
  $P -c "
import json
r = json.load(open('$d/build_report.json')); l = r.get('covalent_loop', {})
tag = ('loop ' + ('CLOSED' if l.get('closed') else 'OPEN') + f\" (wanted {l.get('expected')})\") if l else 'no loop test'
print(f\"  {r['system']:<11} {r['atoms']:>7} atoms  box {r['box'][0]:7.2f} A  charge {r['net_charge']:+.4f}  {tag:<32} problems {len(r['problems'])}\")"
done
