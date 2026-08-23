#!/bin/bash
# Submit the Tier-A schedule: 3 replicates of every system, plus the force-clamp
# ladder. Pass --tier B to submit the long runs instead; pass --dry-run to print.
#
# The clamp ladder brackets the spring forces the Zocchi model predicts for this
# geometry (~2.5 pN at 40 bp, ~5-7 pN at 27 bp) and extends well past them, so the
# chimera runs can be read against a curve rather than a single point. 0 pN is the
# matched control: same topology, same clamp code path, no force.
set -euo pipefail
cd "$(dirname "$0")/.."

TIER=A; DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tier) TIER="$2"; shift 2;;
    --dry-run) DRY=1; shift;;
    --legacy) CHIMERAS+=("${LEGACY[@]}"); shift;;
    --no-ubq) UBQ=(); shift;;
    *) echo "unknown arg $1" >&2; exit 1;;
  esac
done
case "$TIER" in
  A) NS=50;   REPS=3;;
  B) NS=300;  REPS=3;;
  *) echo "tier must be A or B" >&2; exit 1;;
esac

# The systems of record.  S3-S5 are the superseded two-arm builds and are NOT resubmitted
# by default: they are already run, and their value now is as the non-circular comparison,
# not as the construct.  Pass --legacy to include them.
CHIMERAS=(L3_loop27 L4_loop40 L5_loop40nick C1_ss50 C2_ssarms)
LEGACY=(S3_spring27 S4_spring40 S5_spring40nick)
UBQ=(U0_apo U1_clicked U2_ss50 U3_ds50 U4_ssarms)
# S2_clicked is the zero-force control every sfGFP chimera is measured against, so it is
# a control and not an afterthought.
CONTROLS=(S0_wt S1_tet S2_clicked)
# The clamp ladder runs on S1's topology -- the same two Tet2-Et residues the spring
# pulls on -- so its force-response curve and the chimeras are directly comparable.
# 0 pN is the matched control: identical topology and code path, no force.
CLAMPS=(0 2 4 7 12 20)

submit () {
  local name="$1"; shift
  if [[ $DRY -eq 1 ]]; then echo "sbatch --array=1-$REPS -J $name $*"; return; fi
  sbatch --array="1-$REPS" -J "$name" "$@" slurm/md_chain.sbatch
}

for s in "${CONTROLS[@]}" "${CHIMERAS[@]}" "${UBQ[@]}"; do
  submit "$s" --export=ALL,SYS="$s",NS="$NS"
done
for f in "${CLAMPS[@]}"; do
  submit "clamp${f}pN" --export=ALL,SYS=S6_clamp,NS="$NS",CLAMP_PN="$f"
done
echo "tier $TIER: ${#CONTROLS[@]} controls + ${#CHIMERAS[@]} sfGFP chimeras + \
${#UBQ[@]} ubiquitin + ${#CLAMPS[@]} clamp points, $REPS reps each, $NS ns per run"
echo "  (the force-extension windows are submitted separately: slurm/submit_force_extension.sh)"
