#!/usr/bin/env python3
"""Geometry of candidate DNA-spring attachment sites on sfGFP.

Why this exists
---------------
A spring can only be built between two sites if the geometry allows it: the pair must be
far enough apart to span with a duplex of a length we can make, both residues must be
solvent-exposed enough to click a linker onto, and -- for the purpose of moving the
chromophore's spectrum -- they should sit on the beta-strands that carry the chromophore's
own hydrogen bonds.

Everything here is measured from the structure rather than asserted, because the candidate
list came from a language model whose two supporting citations both turned out to be the
wrong paper for the claim, and the residue numbering in this project has already caused
one documented off-by-one.

    python analysis/attachment_candidates.py
"""
from __future__ import annotations

import json
from pathlib import Path

import mdtraj as md
import numpy as np

HERE = Path(__file__).resolve().parent.parent
PDB = HERE / "build" / "01_protein" / "2b3p.pdb"

# Candidate pairs, in 2B3P numbering, which is canonical avGFP numbering (verified against
# the structure: SER147 HIS148 ASN149 / SER202 THR203 GLN204 all appear with those
# identities).  The construct's own numbering runs +1, so anything ordered must be
# incremented -- that offset already produced one wrong pair of PyMOL selections.
PAIRS = [
    ("N149-S202", 149, 202, "flanks H148 and T203 on both supporting strands"),
    ("S147-Q204", 147, 204, "brackets H148/T203 one residue further out; roGFP1 pair"),
    ("N149-Q204", 149, 204, "mixed: keeps the existing 149 handle, reaches past T203"),
    ("S147-S202", 147, 202, "the other mixed diagonal, for completeness"),
    ("D133-N149", 133, 149, "the pair this project has actually built"),
]
# The chromophore-adjacent residues whose H-bonds the spring is meant to perturb.
CHROMO = 66                      # CRO 66: Thr65-Tyr66-Gly67 collapsed
NEIGHBOURS = {96: "Arg96", 148: "His148", 203: "Thr203", 205: "Ser205", 222: "Glu222"}


def cb_or_ca(top, resseq):
    """Index of Cbeta, falling back to Calpha for glycine."""
    for r in top.residues:
        if r.resSeq == resseq and r.chain.index == 0:
            for name in ("CB", "CA"):
                a = next((x.index for x in r.atoms if x.name == name), None)
                if a is not None:
                    return a, r.name, name
    raise KeyError(f"residue {resseq} not found")


def main():
    t = md.load(str(PDB))
    prot = t.atom_slice(t.top.select("protein"))
    top = t.top

    # ---- secondary structure, so we can say which strand each site is on ----
    dssp = md.compute_dssp(prot, simplified=False)[0]
    ss = {}
    for r, code in zip(prot.top.residues, dssp):
        ss[r.resSeq] = code

    # ---- solvent accessibility: can a linker actually be attached? ----
    sasa = md.shrake_rupley(prot, mode="residue")[0] * 100.0   # nm^2 -> A^2
    sasa_by_res = {r.resSeq: float(s) for r, s in zip(prot.top.residues, sasa)}

    print(f"  structure {PDB.name}, {top.n_residues} residues\n")
    print("  === candidate residues ===")
    print(f"  {'res':>10} {'DSSP':>5} {'SASA (A^2)':>11} {'to CRO66 (A)':>13}")
    cro = [a.index for r in top.residues if r.resSeq == CHROMO for a in r.atoms]
    interest = sorted({rs for _, a, b, _ in PAIRS for rs in (a, b)} | set(NEIGHBOURS))
    rows = {}
    for rs in interest:
        try:
            idx, rname, atom = cb_or_ca(top, rs)
        except KeyError:
            continue
        d = float(np.min(np.linalg.norm(
            t.xyz[0, cro, :] - t.xyz[0, idx, :], axis=1)) * 10)
        code = ss.get(rs, "-")
        rows[rs] = {"name": rname, "atom": atom, "dssp": code,
                    "sasa_A2": round(sasa_by_res.get(rs, float("nan")), 1),
                    "min_dist_to_chromophore_A": round(d, 2)}
        tag = "  <- chromophore partner" if rs in NEIGHBOURS else ""
        print(f"  {rname + str(rs):>10} {code:>5} {sasa_by_res.get(rs, 0):11.1f} "
              f"{d:13.2f}{tag}")

    print("\n  === candidate pairs ===")
    print(f"  {'pair':>12} {'residues':>16} {'Cb-Cb (A)':>10} {'both exposed?':>14}  note")
    out = {}
    for label, a, b, note in PAIRS:
        ia, na, _ = cb_or_ca(top, a)
        ib, nb, _ = cb_or_ca(top, b)
        d = float(np.linalg.norm(t.xyz[0, ia] - t.xyz[0, ib]) * 10)
        exposed = (sasa_by_res.get(a, 0) > 20 and sasa_by_res.get(b, 0) > 20)
        out[label] = {"resseq": [a, b], "residues": f"{na}{a}-{nb}{b}",
                      "cb_cb_A": round(d, 2),
                      "sasa_A2": [rows[a]["sasa_A2"], rows[b]["sasa_A2"]],
                      "dssp": [rows[a]["dssp"], rows[b]["dssp"]],
                      "both_exposed": bool(exposed), "note": note}
        print(f"  {label:>12} {na + str(a) + '-' + nb + str(b):>16} {d:10.2f} "
              f"{('yes' if exposed else 'NO'):>14}  {note}")

    ref = out["D133-N149"]["cb_cb_A"]
    print(f"\n  The built pair spans {ref:.1f} A. A 27 bp duplex was chosen for that span,")
    print("  so any new pair needs its own spring length -- the ratio below is how far off")
    print("  the existing 27 bp geometry each candidate is:")
    for label, e in out.items():
        if label == "D133-N149":
            continue
        print(f"    {label:>12} {e['cb_cb_A']:6.1f} A  = {e['cb_cb_A'] / ref:.2f} x the "
              f"built span")

    (HERE / "analysis" / "attachment_candidates.json").write_text(
        json.dumps({"structure": PDB.name, "residues": rows, "pairs": out,
                    "chromophore_resseq": CHROMO,
                    "numbering": "2B3P == canonical avGFP; construct numbering is +1"},
                   indent=2))
    print(f"\n  wrote analysis/attachment_candidates.json")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
def strands(prot, dssp):
    """Contiguous DSSP 'E' runs -> numbered beta-strands, with their residue extents."""
    out, cur = [], []
    for r, code in zip(prot.top.residues, dssp):
        if code == "E":
            cur.append(r.resSeq)
        else:
            if len(cur) >= 3:
                out.append((min(cur), max(cur)))
            cur = []
    if len(cur) >= 3:
        out.append((min(cur), max(cur)))
    return out


def design_search(t, prot, dssp, sasa_by_res, lo=24.0, hi=42.0, min_sasa=30.0):
    """Which exposed pairs on the two chromophore-supporting strands can a spring span?

    The roGFP positions are 4.6-10.6 A apart because a disulfide is a 2 A bond.  A DNA
    spring cannot work there: each tether -- the Tet2-Et side chain, the sTCO adduct and
    the DNA's own 5'/3' arm -- is itself 10-15 A, and the built geometry needs 62-76 A
    between the duplex ends to sit on a 31 A Cb-Cb separation.

    So the question is not "attach at 147 and 204" but "which anchors, far enough apart to
    build, would PRY APART the strands that carry His148 and Thr203 when pushed".  The
    figure of merit is how much of the anchor-separation axis lies along the H148->T203
    direction: a spring pushing along that axis separates the two chromophore H-bond
    donors, while one pushing across it merely shears or twists the barrel.
    """
    top = t.top
    sl = strands(prot, dssp)
    s148 = next((s for s in sl if s[0] <= 148 <= s[1]), None)
    s203 = next((s for s in sl if s[0] <= 203 <= s[1]), None)
    if not (s148 and s203):
        print("  could not assign strands for 148/203"); return {}
    print(f"\n  === the two chromophore-supporting strands ===")
    print(f"  His148 lies on the strand spanning residues {s148[0]}-{s148[1]}")
    print(f"  Thr203 lies on the strand spanning residues {s203[0]}-{s203[1]}")

    def cb(rs):
        for r in top.residues:
            if r.resSeq == rs and r.chain.index == 0:
                for nm in ("CB", "CA"):
                    a = next((x.index for x in r.atoms if x.name == nm), None)
                    if a is not None:
                        return a
        return None

    i148, i203 = cb(148), cb(203)
    axis = t.xyz[0, i203] - t.xyz[0, i148]
    axis /= np.linalg.norm(axis)
    d148_203 = float(np.linalg.norm(t.xyz[0, i203] - t.xyz[0, i148]) * 10)
    print(f"  His148 Cb to Thr203 Cb: {d148_203:.2f} A -- this is the distance a spring")
    print("  would have to increase to loosen the chromophore's two strand H-bonds.\n")

    cands = []
    A = [r for r in range(s148[0], s148[1] + 1) if sasa_by_res.get(r, 0) >= min_sasa]
    B = [r for r in range(s203[0], s203[1] + 1) if sasa_by_res.get(r, 0) >= min_sasa]
    for a in A:
        for b in B:
            ia, ib = cb(a), cb(b)
            if ia is None or ib is None:
                continue
            v = t.xyz[0, ib] - t.xyz[0, ia]
            d = float(np.linalg.norm(v) * 10)
            if not (lo <= d <= hi):
                continue
            proj = abs(float(np.dot(v / np.linalg.norm(v), axis)))
            cands.append((proj, d, a, b))
    cands.sort(reverse=True)
    print(f"  exposed pairs on those two strands with Cb-Cb in {lo:.0f}-{hi:.0f} A "
          f"(buildable range): {len(cands)}")
    if not cands:
        print("  NONE. The two strands are too short and too close for a spring to span")
        print("  them at all -- which is the central geometric objection to the proposal.")
        return {}
    print(f"\n  {'pair':>14} {'Cb-Cb':>7} {'alignment':>10}  interpretation")
    best = []
    for proj, d, a, b in cands[:10]:
        rn = {r.resSeq: r.name for r in top.residues if r.chain.index == 0}
        lab = f"{rn.get(a,'?')}{a}-{rn.get(b,'?')}{b}"
        interp = ("pushes along the H148-T203 axis" if proj > 0.7 else
                  "mostly shears the strands" if proj > 0.4 else
                  "pushes across the axis: twist, not separation")
        print(f"  {lab:>14} {d:6.1f} A {proj:10.2f}  {interp}")
        best.append({"pair": lab, "resseq": [a, b], "cb_cb_A": round(d, 2),
                     "alignment_with_H148_T203_axis": round(proj, 3),
                     "interpretation": interp})
    return {"strand_H148": list(s148), "strand_T203": list(s203),
            "H148_T203_cb_cb_A": round(d148_203, 2), "buildable": best}


if __name__ == "__main__":
    _t = md.load(str(PDB))
    _prot = _t.atom_slice(_t.top.select("protein"))
    _dssp = md.compute_dssp(_prot, simplified=False)[0]
    _sasa = md.shrake_rupley(_prot, mode="residue")[0] * 100.0
    _sb = {r.resSeq: float(s) for r, s in zip(_prot.top.residues, _sasa)}
    _d = design_search(_t, _prot, _dssp, _sb)
    if _d:
        p = HERE / "analysis" / "attachment_candidates.json"
        j = json.loads(p.read_text()); j["design_search"] = _d
        p.write_text(json.dumps(j, indent=2))
        print(f"\n  appended design_search to {p.name}")
