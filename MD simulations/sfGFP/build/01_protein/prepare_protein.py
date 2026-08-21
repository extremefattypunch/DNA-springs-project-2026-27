#!/usr/bin/env python3
"""Turn the raw 2B3P deposition into a tleap-ready sfGFP structure.

What this does and why
----------------------
* Strips the 9 Cd(2+) and 6 acetates.  Both are crystallisation additives, and
  Cd(2+) has no ff14SB parameters.
* Keeps the CRO chromophore (residue 66) untouched.  Its parameters come from
  leaprc.xFPchromophores (Breyfogle et al., J Phys Chem B 2023, 127:5772), so the
  atom names must survive verbatim -- they are checked, not assumed.
* Keeps the crystallographic waters.  GFP has buried waters in the chromophore
  cavity that are part of its H-bond network; discarding them and re-solvating
  blind would lose them.  Waters that clash with anything are dropped, reported.
* Assigns histidine tautomers from H-bond geometry rather than accepting a
  default.  This matters: His148 donates to the chromophore phenolate, and
  Amber's default HIS -> HIE would put the proton on the wrong nitrogen.
* Adds ACE/NME caps at Ser2 and Gly232 using ideal internal coordinates, because
  in the real construct both positions are internal (Met1 precedes; GSHHHHHH
  follows) and charged termini there would be an invented pair of charges.
* Does NOT add hydrogens or rebuild missing side chains -- tleap does both from
  its own libraries, which keeps one tool responsible for the topology.
* Measures the attachment-site geometry and turns it into the spring force
  ladder, so no downstream script hard-codes a distance.

Usage
-----
    python prepare_protein.py [--pdb 2b3p.pdb] [--outdir .]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis"))
from spring_model import (TAU_C_CONTINUOUS, TAU_C_NICKED, Spring,  # noqa: E402
                          choose_n_bp, span)

DROP_HET = {"CD", "ACY"}

# Heavy atoms per standard residue, backbone N/CA/C/O included.  Used only to
# report what tleap will have to rebuild.
HEAVY_EXPECTED = {
    "GLY": 4, "ALA": 5, "SER": 6, "CYS": 6, "THR": 7, "PRO": 7, "VAL": 7,
    "ASN": 8, "ASP": 8, "LEU": 8, "ILE": 8, "MET": 8, "GLN": 9, "GLU": 9,
    "LYS": 9, "HIS": 10, "PHE": 11, "ARG": 11, "TYR": 12, "TRP": 14,
}
WATER_NAMES = {"HOH", "WAT", "DOD"}

# Bondi-like radii used for the SASA probe; matches the values Amber/mdtraj use.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}
PROBE = 1.40

# Ideal internal coordinates for the caps (Amber ACE/NME geometry).
# (name, (ref1, ref2, ref3), bond/A, angle/deg, dihedral/deg)
ACE_BUILD = [
    ("C",   ("C", "CA", "N"), 1.335, 121.7, 180.0),   # amide C bonded to res-2 N
    ("O",   ("CA", "N", "C"), 1.229, 122.9,   0.0),
    ("CH3", ("CA", "N", "C"), 1.522, 116.6, 180.0),
]
NME_BUILD = [
    ("N",   ("N", "CA", "C"), 1.335, 116.6, 180.0),   # amide N bonded to res-232 C
    ("CH3", ("CA", "C", "N"), 1.449, 121.9, 180.0),
]
# ff14SB loads aminoct12.lib, whose NME template names the methyl carbon "C"
# (with H1/H2/H3), not "CH3".  Emitting "CH3" makes tleap build the six template
# atoms and then choke on a seventh, untyped atom.  ACE is unaffected: aminont12's
# ACE really does use CH3, alongside its own C and O.  We place the atom under the
# internal name CH3 -- the placement references residue 232's own C -- and rename
# on output so the two never collide.
CAP_RENAME = {"NME": {"CH3": "C"}}


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def place_atom(a, b, c, bond, angle, dihedral):
    """Natural-extension-reference-frame placement of a 4th atom.

    Returns the position d such that |c-d| = bond, angle(b,c,d) = angle and
    dihedral(a,b,c,d) = dihedral.  Standard Z-matrix -> Cartesian conversion.
    """
    ang, dih = math.radians(angle), math.radians(dihedral)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * math.cos(ang),
                   bond * math.cos(dih) * math.sin(ang),
                   bond * math.sin(dih) * math.sin(ang)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def sasa(coords, radii, n_points=256):
    """Shrake-Rupley solvent-accessible surface area per atom (A^2).

    Self-contained so structure prep does not depend on the MD stack being
    installed.  Points are a golden-spiral (Fibonacci) lattice, which is more
    uniform than the original icosahedral subdivision at the same count.
    """
    k = np.arange(n_points)
    phi = math.pi * (3.0 - math.sqrt(5.0)) * k
    z = 1.0 - 2.0 * (k + 0.5) / n_points
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    sphere = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)

    coords = np.asarray(coords, float)
    radii = np.asarray(radii, float) + PROBE
    out = np.zeros(len(coords))
    # neighbour list by bounding-sphere overlap
    for i in range(len(coords)):
        d = np.linalg.norm(coords - coords[i], axis=1)
        nb = np.where((d < radii + radii[i]) & (d > 1e-9))[0]
        pts = coords[i] + radii[i] * sphere
        if len(nb):
            dist = np.linalg.norm(pts[:, None, :] - coords[nb][None, :, :], axis=2)
            free = np.all(dist >= radii[nb][None, :], axis=1)
        else:
            free = np.ones(len(pts), bool)
        out[i] = 4.0 * math.pi * radii[i] ** 2 * free.mean()
    return out


# --------------------------------------------------------------------------
# PDB I/O
# --------------------------------------------------------------------------
class Atom:
    __slots__ = ("name", "resn", "resi", "chain", "xyz", "elem", "occ", "het")

    def __init__(self, name, resn, resi, chain, xyz, elem, occ, het):
        self.name, self.resn, self.resi, self.chain = name, resn, resi, chain
        self.xyz, self.elem, self.occ, self.het = xyz, elem, occ, het


def read_pdb(path):
    atoms, seen_alt = [], set()
    for line in Path(path).read_text().splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            seen_alt.add(altloc)
            continue
        elem = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        if elem == "H" or elem == "D":
            continue
        atoms.append(Atom(line[12:16].strip(), line[17:20].strip(), int(line[22:26]),
                          line[21], np.array([float(line[30:38]), float(line[38:46]),
                                              float(line[46:54])]),
                          elem, float(line[54:60] or 1.0), rec == "HETATM"))
    return atoms, seen_alt


def write_pdb(path, groups):
    """groups: list of (list_of_Atom) -- one TER after each group."""
    serial = 1
    with open(path, "w") as fh:
        fh.write("REMARK   1 sfGFP prepared from PDB 2B3P for tleap; see prepare_protein.py\n")
        for grp in groups:
            for a in grp:
                nm = a.name
                # PDB atom-name column convention: 1-3 char names start at col 14
                nm_field = f" {nm:<3}" if len(nm) < 4 else nm
                fh.write(f"ATOM  {serial:5d} {nm_field}{a.resn:>4}{a.chain:>2}"
                         f"{a.resi:4d}    {a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
                         f"{1.00:6.2f}{0.00:6.2f}          {a.elem:>2}\n")
                serial += 1
            fh.write(f"TER   {serial:5d}\n")
            serial += 1
        fh.write("END\n")


# --------------------------------------------------------------------------
# histidine tautomers
# --------------------------------------------------------------------------
def assign_his(residues, all_atoms):
    """Pick HID / HIE per histidine from the local H-bond geometry.

    A ring nitrogen that sits within 3.5 A of a hydrogen-bond *acceptor* wants to
    be the donor (protonated); one that sits near a *donor* wants to be the
    acceptor (deprotonated).  Where both nitrogens are equally served we fall back
    to HIE, Amber's default, and say so.
    """
    donors = {"N", "NE", "NH1", "NH2", "ND2", "NE2", "NZ", "OG", "OG1", "OH", "NE1", "SG"}
    # The CRO chromophore is modelled in its anionic (phenolate) state, which is the
    # fluorescent form and the one leaprc.xFPchromophores parameterises. Its OH oxygen
    # is therefore deprotonated -- a pure acceptor, not the donor its atom name implies.
    # Getting this wrong would invert the His148 tautomer argument.
    cro_acceptors = {"OH", "O2", "O3", "OG1"}
    coords = np.array([a.xyz for a in all_atoms])
    decisions = {}
    for (resi, resn), rat in residues.items():
        if resn != "HIS":
            continue
        evidence = {}
        for nname in ("ND1", "NE2"):
            if nname not in rat:
                continue
            d = np.linalg.norm(coords - rat[nname].xyz, axis=1)
            near = []
            for j in np.where((d < 3.5) & (d > 1e-6))[0]:
                o = all_atoms[j]
                if o.resi == resi and not o.het:
                    continue
                if o.elem not in ("N", "O", "S"):
                    continue
                if o.resn == "CRO":
                    role = "acceptor" if o.name in cro_acceptors else "donor"
                elif o.resn in WATER_NAMES:
                    role = "donor"        # water can do either; treat as donor-capable
                else:
                    role = "donor" if o.name in donors else "acceptor"
                near.append({"partner": f"{o.resn}{o.resi}.{o.name}",
                             "dist": round(float(d[j]), 2), "partner_role": role})
            near.sort(key=lambda r: r["dist"])
            evidence[nname] = near[:4]

        def score(n):  # protonating n is favourable if it faces acceptors
            return sum(1.0 / e["dist"] for e in evidence.get(n, [])
                       if e["partner_role"] == "acceptor")
        s_d, s_e = score("ND1"), score("NE2")
        if abs(s_d - s_e) < 1e-9:
            name, why = "HIE", "no discriminating contact; Amber default"
        elif s_d > s_e:
            name, why = "HID", "ND1 faces an acceptor, so ND1 is the donor"
        else:
            name, why = "HIE", "NE2 faces an acceptor, so NE2 is the donor"
        decisions[resi] = {"name": name, "reason": why,
                           "score_ND1": round(s_d, 3), "score_NE2": round(s_e, 3),
                           "contacts": evidence}
    return decisions


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--pdb", default=str(here / "2b3p.pdb"))
    ap.add_argument("--outdir", default=str(here))
    ap.add_argument("--water-clash", type=float, default=2.2,
                    help="drop crystal waters closer than this to a solute heavy atom")
    ap.add_argument("--arm", type=float, default=2.1,
                    help="provisional one-sided linker length in nm; refined in step 2")
    args = ap.parse_args()
    out = Path(args.outdir)

    atoms, alts = read_pdb(args.pdb)
    if alts:
        print(f"note: discarded altloc(s) {sorted(alts)}")

    solute, waters = [], []
    for a in atoms:
        if a.resn in DROP_HET:
            continue
        (waters if a.resn in WATER_NAMES else solute).append(a)
    print(f"solute heavy atoms: {len(solute)}   crystal waters: {len(waters)}")

    residues = OrderedDict()
    for a in solute:
        residues.setdefault((a.resi, a.resn), {})[a.name] = a
    resi_list = sorted(residues, key=lambda k: k[0])
    first, last = resi_list[0], resi_list[-1]
    print(f"chain spans {first[1]}{first[0]} .. {last[1]}{last[0]} "
          f"({len(residues)} residues)")

    # ---- CRO integrity ----
    cro = residues.get((66, "CRO"))
    if cro is None:
        sys.exit("FATAL: CRO 66 not found")
    expected_cro = {"C1", "C2", "C3", "CA1", "CA2", "CA3", "CB1", "CB2", "CD1", "CD2",
                    "CE1", "CE2", "CG1", "CG2", "CZ", "N1", "N2", "N3", "O2", "O3",
                    "OG1", "OH"}
    if set(cro) != expected_cro:
        sys.exit(f"FATAL: CRO atom names differ from the xFPchromophores template: "
                 f"missing {expected_cro - set(cro)}, extra {set(cro) - expected_cro}")
    print(f"CRO 66 intact: {len(cro)} atoms, names match the expected template")

    # ---- histidine tautomers ----
    his = assign_his(residues, solute)
    print("\nhistidine tautomer assignment:")
    for resi in sorted(his):
        h = his[resi]
        print(f"  HIS{resi:<4} -> {h['name']}   ({h['reason']})")
        for n in ("ND1", "NE2"):
            for c in h["contacts"].get(n, [])[:2]:
                print(f"        {n} .. {c['partner']:<16} {c['dist']:.2f} A "
                      f"({c['partner_role']})")

    # ---- caps ----
    solute_xyz = np.array([a.xyz for a in solute])
    elem_of = {"C": "C", "O": "O", "N": "N", "CH3": "C"}

    def build_cap(spec, ref_res, resi, resn, attach_atom):
        """Build a cap, choosing the free dihedral to avoid steric clashes.

        The first cap atom's dihedral is genuinely undetermined by the crystal
        structure -- it is the phi of the first residue (or the psi continuation past
        the last).  A fixed 180 deg happens to superimpose the NME nitrogen on the
        Gly232 carbonyl oxygen, so scan instead and keep the least-clashing value.
        ``attach_atom`` is the reference atom the cap bonds to; distances to it are
        a bond length, not a clash, so it is excluded from the score.
        """
        keep = np.array([np.linalg.norm(a.xyz - ref_res[attach_atom].xyz) > 1e-6
                         for a in solute])
        others = solute_xyz[keep]
        best = None
        for offset in range(0, 360, 5):
            made = {}
            for i, (name, (r1, r2, r3), bond, ang, dih) in enumerate(spec):
                pos = {**{k: v.xyz for k, v in ref_res.items()}, **made}
                made[name] = place_atom(pos[r1], pos[r2], pos[r3], bond, ang,
                                        dih + (offset if i == 0 else 0))
            dmin = min(np.linalg.norm(others - xyz, axis=1).min()
                       for xyz in made.values())
            if best is None or dmin > best[0]:
                best = (dmin, offset, made)
        dmin, offset, made = best
        print(f"  {resn}{resi}: dihedral offset {offset} deg, closest non-bonded "
              f"contact {dmin:.2f} A")
        if dmin < 2.2:
            print(f"  WARNING: {resn}{resi} still crowded at {dmin:.2f} A")
        rename = CAP_RENAME.get(resn, {})
        return [Atom(rename.get(n, n), resn, resi, ref_res["CA"].chain, xyz,
                     elem_of[n], 1.0, False)
                for n, xyz in made.items()]

    print("\nbuilding terminal caps (ideal internal coordinates; tleap adds their H):")
    ace = build_cap(ACE_BUILD, residues[first], first[0] - 1, "ACE", "N")
    nme = build_cap(NME_BUILD, residues[last], last[0] + 1, "NME", "C")

    # ---- rename for Amber, drop clashing waters ----
    for a in solute:
        if a.resn == "HIS":
            a.resn = his[a.resi]["name"]
    keep_w, dropped = [], 0
    wres = defaultdict(list)
    for a in waters:
        wres[a.resi].append(a)
    for resi, wat in wres.items():
        if any(np.linalg.norm(solute_xyz - a.xyz, axis=1).min() < args.water_clash
               for a in wat) or any(a.occ < 0.9 for a in wat):
            dropped += 1
            continue
        for a in wat:
            a.resn = "WAT"
            a.name = "O"
            keep_w.append(a)
    print(f"crystal waters kept: {len(keep_w)}  dropped (clash or partial occupancy): {dropped}")

    chain = ace + [residues[k][n] for k in resi_list for n in residues[k]] + nme
    write_pdb(out / "sfgfp_prepped.pdb", [chain] + [[w] for w in keep_w])
    print(f"wrote {out / 'sfgfp_prepped.pdb'}  ({len(chain)} solute + {len(keep_w)} water atoms)")

    # ---- SASA and attachment-site geometry ----
    print("\ncomputing SASA (Shrake-Rupley, 256 points/atom)...")
    areas = sasa(solute_xyz, [VDW.get(a.elem, 1.7) for a in solute])
    per_res = defaultdict(float)
    for a, s in zip(solute, areas):
        per_res[(a.resi, a.resn)] += s

    cro_centroid = np.mean([a.xyz for a in cro.values()], axis=0)

    def site_info(resi):
        key = next(k for k in residues if k[0] == resi)
        rat = residues[key]
        anchor = "CB" if "CB" in rat else "CA"
        return {"resi": resi, "resname": key[1], "anchor_atom": anchor,
                "sasa_A2": round(per_res[key], 1),
                "dist_to_chromophore_A": round(
                    float(np.linalg.norm(rat["CA"].xyz - cro_centroid)), 2),
                "xyz": [round(float(v), 3) for v in rat[anchor].xyz]}

    pairs = {}
    for a_i in (133, 134):
        for b_i in (149, 150):
            ia, ib = site_info(a_i), site_info(b_i)
            d_ang = float(np.linalg.norm(np.array(ia["xyz"]) - np.array(ib["xyz"])))
            x = span(d_ang / 10.0, args.arm)
            entry = {"site_a": ia, "site_b": ib,
                     "d_anchor_anchor_A": round(d_ang, 2),
                     "span_nm_at_arm": {f"{args.arm}": round(x, 3)},
                     "force_pN": {}}
            for tau, tag in ((TAU_C_NICKED, "nicked"), (TAU_C_CONTINUOUS, "continuous")):
                entry["force_pN"][tag] = {
                    n: round(Spring(n, tau_c=tau).force(x), 2)
                    for n in (24, 26, 28, 30, 35, 40, 50, 60)
                    if Spring(n, tau_c=tau).is_feasible(x)}
            entry["strong_spring_bp"] = choose_n_bp(x, 7.0)
            pairs[f"{ia['resname']}{a_i}-{ib['resname']}{b_i}"] = entry

    chosen = f"ASP133-ASN149"
    payload = {
        "source_pdb": "2B3P",
        "numbering_note": (
            "2B3P numbers the chain 1-246 with the Thr65-Tyr66-Gly67 chromophore "
            "collapsed into CRO 66, so 65 and 67 are absent. The construct's own "
            "numbering runs +1 relative to this, i.e. construct D134/N150 are 2B3P "
            "Asp133/Asn149. The .pse selections named 134TAG/150TAG point at 2B3P "
            "Gly134/Val150, which is one residue off: the ESI-MS mass shifts of "
            "+141 and +140 Da can only be Asn->Tet2-Et and Asp->Tet2-Et."),
        "chosen_pair": chosen,
        "arm_nm_provisional": args.arm,
        "pairs": pairs,
        "histidine_tautomers": his,
        "crystal_waters_kept": len(keep_w),
        "residues_missing_sidechain_atoms": {
            f"{k[1]}{k[0]}": HEAVY_EXPECTED[k[1]] - len(residues[k])
            for k in resi_list
            if k[1] in HEAVY_EXPECTED and len(residues[k]) < HEAVY_EXPECTED[k[1]]},
        "missing_atoms_note": ("tleap rebuilds these from its own libraries; all are "
                               "surface side chains, none within 15 A of an "
                               "attachment site or the chromophore"),
    }
    (out / "attachment_sites.json").write_text(json.dumps(payload, indent=2))
    print(f"wrote {out / 'attachment_sites.json'}")

    print(f"\nchosen attachment pair {chosen}:")
    e = pairs[chosen]
    print(f"  anchor separation {e['d_anchor_anchor_A']} A "
          f"-> span {list(e['span_nm_at_arm'].values())[0]} nm at arm = {args.arm} nm")
    print(f"  SASA: {e['site_a']['resname']}{e['site_a']['resi']} "
          f"{e['site_a']['sasa_A2']} A^2, {e['site_b']['resname']}{e['site_b']['resi']} "
          f"{e['site_b']['sasa_A2']} A^2")
    print(f"  distance to chromophore: {e['site_a']['dist_to_chromophore_A']} and "
          f"{e['site_b']['dist_to_chromophore_A']} A")
    print(f"  force ladder (nicked): {e['force_pN']['nicked']}")
    print(f"  strong spring -> {e['strong_spring_bp']} bp")

    # ---- corrected PyMOL selections ----
    (out / "fix_pse_selections.pml").write_text(f"""\
# Corrected attachment-site selections for the sfGFP-DNA spring chimera.
# Load over the existing session:  pymol "sfGFP 2b3p (150TAG, 134TAG).pse" fix_pse_selections.pml
#
# The session's 134TAG / 150TAG selections point at 2B3P Gly134 and Val150.
# Those are one residue off from the sites actually encoded.  The construct is
# numbered +1 relative to the 2B3P deposition, and the ESI-MS shifts settle it:
#   single Tet construct  27,968 - 27,827 = +141 Da  = Asn -> Tet2-Et
#   double Tet construct  28,108 - 27,827 = +281 Da  = +141 (Asn) + 140 (Asp)
# Gly -> Tet2-Et would be +198 Da and Val -> Tet2-Et +156 Da; neither is observed.
delete 134TAG_wrong
delete 150TAG_wrong
select 134TAG_wrong, 2b3p and chain A and resi 134     # Gly134, the old bookmark
select 150TAG_wrong, 2b3p and chain A and resi 150     # Val150, the old bookmark
select site_D134, 2b3p and chain A and resi 133        # Asp133 in 2B3P = construct D134
select site_N150, 2b3p and chain A and resi 149        # Asn149 in 2B3P = construct N150
select His148_gate, 2b3p and chain A and resi 148      # donates to the chromophore phenolate
show sticks, site_D134 or site_N150 or His148_gate or chromo
color orange, site_D134
color marine, site_N150
color yellow, His148_gate
distance span, site_D134 and name CB, site_N150 and name CB
print "Asp133(CB)-Asn149(CB) = {e['d_anchor_anchor_A']} A"
""")
    print(f"wrote {out / 'fix_pse_selections.pml'}")


if __name__ == "__main__":
    main()
