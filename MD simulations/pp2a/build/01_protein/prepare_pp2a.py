#!/usr/bin/env python3
"""PP2A A-alpha (PR65, UniProt P30153) -> the D6X/N552X construct, ready for tleap.

The starting model is the AlphaFold DB v6 model of the full 589-residue scaffold subunit,
already stripped to heavy atoms by the cryo-EM reproduction
(cryoEM/pp2a_day2025_reproduction/scripts/01_prepare_models.py).  It is the model the
finished 100 ns apo run used, so the chimera and its baseline start from the same
coordinates.  The monomer is simulated alone: no B or C subunit.

Numbering -- read this before anything else
-------------------------------------------
The construct's sites are labelled D6X and N552X (pdb files/pp2aTETattachments.png).
In UniProt/PDB numbering residue 6 is Gly and 552 is Asp; residues 7 and 553 are Asp
and Asn.  The labels count from Ala2, i.e. construct n = UniProt n+1, the same
off-by-one the sfGFP build met (construct D134/N150 = 2B3P Asp133/Asn149).  So

    D6X  = Asp7     N552X = Asn553     (file numbering, used everywhere downstream)

Confirmed with Ian on 2026-10-26.  The Yilmaz/Bahar end-to-end metric "N29-F577" is the
same Met1-less convention (P30153 30 = Asn, 578 = Phe) and is recorded here as the
literature pair so the two can be compared without re-deriving the offset.

What this script does
---------------------
* keeps the full chain 1-589 with free termini, as the apo run did (no caps: nothing in
  the repo says the expressed construct is tagged, and Met1 is far from both sites);
* assigns HID/HIE per histidine from the local geometry (the ubiquitin rule), rather
  than tleap's all-HIE default the apo run used -- recorded per residue;
* measures the attachment geometry (CB-CB, SASA, pLDDT, native phi/psi) and the HEAT
  grid, and writes attachment_sites.json in the sfGFP schema;
* writes pp2a_prepped.pdb: chain A, heavy atoms only, UniProt numbering.  The two
  sites stay ASP/ASN here; the assembler replaces them with TDP/TET itself.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SFGFP = HERE.parents[2] / "sfGFP"
sys.path.insert(0, str(SFGFP / "build" / "01_protein"))
from prepare_ubiquitin import fmt, his_tautomer, read_pdb  # noqa: E402

UNIPROT = "P30153"
N_RES = 589
SITES = (7, 553)                       # file numbering
CONSTRUCT_LABEL = {7: "D6X", 553: "N552X"}
LITERATURE_PAIR = (30, 578)            # Yilmaz et al. 2025 "N29-F577", Met1-less
HEAT_START, HEAT_LEN, N_HEAT = 8, 39, 15   # nominal grid, cryoEM/.../25_geometry.py
WATER = {"HOH", "WAT"}


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2      # b0 points back, per the standard formula
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def sasa_per_residue(pdb):
    import mdtraj as md
    t = md.load(str(pdb))
    s = md.shrake_rupley(t, mode="residue")[0] * 100.0  # nm^2 -> A^2
    return {r.resSeq: float(s[r.index]) for r in t.topology.residues}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default=os.path.join(
        os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025"),
        "md", "systems", "Aa", "protein.pdb"))
    ap.add_argument("--out", default=str(HERE / "pp2a_prepped.pdb"))
    args = ap.parse_args()

    atoms = read_pdb(args.pdb)
    prot = [a for a in atoms if a["resn"] not in WATER and a["elem"] != "H"]
    res = {}
    for a in prot:
        res.setdefault(a["resi"], {})[a["name"]] = a
    nums = sorted(res)
    if nums != list(range(1, N_RES + 1)):
        sys.exit(f"expected residues 1-{N_RES}, got {nums[0]}-{nums[-1]} ({len(nums)})")
    seq = {n: next(iter(res[n].values()))["resn"] for n in nums}
    for n, want in ((6, "GLY"), (7, "ASP"), (552, "ASP"), (553, "ASN")):
        if seq[n] != want:
            sys.exit(f"residue {n} is {seq[n]}, expected {want}: not {UNIPROT}?")
    print(f"{UNIPROT}: residues 1-{N_RES}, {len(prot)} heavy atoms; "
          f"res 6/7 = {seq[6]}/{seq[7]}, 552/553 = {seq[552]}/{seq[553]}")

    renamed, hist = {}, {}
    for n in nums:
        if seq[n] != "HIS":
            continue
        others = [a for k, rr in res.items() if k != n for a in rr.values()]
        taut, why = his_tautomer(res[n], others)
        renamed[n] = taut
        hist[str(n)] = {"name": taut, "reason": why}
        print(f"  His{n} -> {taut}: {why}")

    sasa = sasa_per_residue(args.pdb)
    # pLDDT lives in the B-factor column; read_pdb does not keep it, so re-read it here.
    bfac = {}
    for ln in Path(args.pdb).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            bfac.setdefault(int(ln[22:26]), float(ln[60:66]))

    def site(n):
        r = res[n]
        anchor = "CB" if "CB" in r else "CA"
        d = {"resi": n, "resname": seq[n], "construct_label": CONSTRUCT_LABEL.get(n),
             "anchor_atom": anchor, "sasa_A2": round(sasa[n], 1),
             "plddt": bfac.get(n), "xyz": [round(float(v), 3) for v in r[anchor]["xyz"]]}
        if n - 1 in res and n + 1 in res:
            d["phi_deg"] = round(dihedral(res[n - 1]["C"]["xyz"], r["N"]["xyz"],
                                          r["CA"]["xyz"], r["C"]["xyz"]), 1)
            d["psi_deg"] = round(dihedral(r["N"]["xyz"], r["CA"]["xyz"], r["C"]["xyz"],
                                          res[n + 1]["N"]["xyz"]), 1)
        return d

    a, b = SITES
    sa, sb = site(a), site(b)
    d_cb = float(np.linalg.norm(np.array(sa["xyz"]) - np.array(sb["xyz"])))
    d_ca = float(np.linalg.norm(res[a]["CA"]["xyz"] - res[b]["CA"]["xyz"]))
    la, lb = LITERATURE_PAIR
    d_lit = float(np.linalg.norm(res[la]["CA"]["xyz"] - res[lb]["CA"]["xyz"]))
    print(f"  sites Asp{a}/Asn{b}: CB-CB {d_cb:.2f} A, CA-CA {d_ca:.2f} A; "
          f"SASA {sa['sasa_A2']}/{sb['sasa_A2']} A^2; pLDDT {sa['plddt']}/{sb['plddt']}; "
          f"phi/psi {sa.get('phi_deg')}/{sa.get('psi_deg')} and "
          f"{sb.get('phi_deg')}/{sb.get('psi_deg')}")
    print(f"  literature pair CA{la}-CA{lb} (Yilmaz 'N29-F577'): {d_lit:.2f} A")

    heat = {k: [HEAT_START + HEAT_LEN * (k - 1), min(HEAT_START + HEAT_LEN * k - 1, N_RES)]
            for k in range(1, N_HEAT + 1)}
    aromatic = []
    for n in SITES:
        c = np.array(sa["xyz"] if n == a else sb["xyz"])
        for m, rr in res.items():
            if rr and next(iter(rr.values()))["resn"] in ("PHE", "TYR", "TRP", "HIS") and m != n:
                dmin = min(float(np.linalg.norm(at["xyz"] - c)) for at in rr.values())
                if dmin < 8.0:
                    aromatic.append({"site": n, "resn": seq[m], "resi": m,
                                     "min_dist_to_anchor_A": round(dmin, 2)})

    report = {
        "source_pdb": args.pdb,
        "source_model": "AlphaFold DB v6, AF-P30153-F1 (the apo run's start model)",
        "uniprot": UNIPROT,
        "numbering_note": (
            "File numbering is UniProt P30153 (Met1 = 1).  The construct labels D6X/N552X "
            "count from Ala2 (construct n = UniProt n+1): D6X = Asp7, N552X = Asn553.  "
            "UniProt 6 is Gly (no CB) and 552 is Asp, so the literal reading is impossible.  "
            "Confirmed with Ian 2026-10-26.  Yilmaz et al.'s 'N29-F577' is the same "
            "convention = Asn30-Phe578 here."),
        "chosen_pair": f"{seq[a]}{a}-{seq[b]}{b}",
        "construct_labels": {str(a): CONSTRUCT_LABEL[a], str(b): CONSTRUCT_LABEL[b]},
        "termini": "free 1-589, charged termini, no caps (as the apo run Aa_rep1)",
        "pairs": {f"{seq[a]}{a}-{seq[b]}{b}": {
            "site_a": sa, "site_b": sb,
            "d_anchor_anchor_A": round(d_cb, 2), "d_CA_CA_A": round(d_ca, 2)}},
        "literature_pair": {"resi": list(LITERATURE_PAIR),
                            "label": "Yilmaz/Bahar 2025 'N29-F577' (Met1-less)",
                            "d_CA_CA_A": round(d_lit, 2),
                            "apo_MD_reference_A": "66.5 +/- 9.7 (CHARMM36, 3 x 704 ns)"},
        "heat_repeats_nominal": heat,
        "aromatic_within_8A_of_anchor": aromatic,
        "histidine_tautomers": hist,
        "plddt": {"mean_all": round(float(np.mean(list(bfac.values()))), 1),
                  "res_1_13": [bfac[n] for n in range(1, 14)],
                  "at_sites": {str(a): bfac.get(a), str(b): bfac.get(b)}},
    }

    lines, serial = [], 1
    for n in nums:
        rn = renamed.get(n, seq[n])
        for nm, at in res[n].items():
            at = dict(at, chain="A")
            lines.append(fmt(serial, at, rn))
            serial += 1
    lines.append(f"TER   {serial:5d}")
    Path(args.out).write_text("\n".join(lines) + "\nEND\n")
    Path(HERE / "attachment_sites.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {args.out} ({serial - 1} atoms) and attachment_sites.json")


if __name__ == "__main__":
    main()
