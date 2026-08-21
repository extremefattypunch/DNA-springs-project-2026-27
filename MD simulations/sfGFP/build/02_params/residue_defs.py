#!/usr/bin/env python3
"""Atom-level definitions of the three non-standard residues in the chimera.

This file is the single source of truth for the chemistry.  Everything else --
3D building, RESP charge fitting, prepgen mainchain files, tleap bond commands,
analysis atom selections -- reads these tables, so the chemistry is stated once.

The three residues
------------------
``TET``  Tet2-Et, the genetically encoded tetrazine amino acid, oxidized
         (reactive) form.  Used for the unclicked control.
``TDP``  the click product: Tet2-Et + sTCO after inverse-electron-demand
         Diels-Alder and loss of N2, carrying the carbamate arm, and ending at
         the carbamate carbonyl carbon.
``DNL``  the 5'-amino-modifier-C6 linker on the DNA side, from the carbamate
         nitrogen to the bridging oxygen that bonds the first nucleotide's
         phosphorus.

Chemistry and its sources
-------------------------
Tet2-Et = 4-(6-ethyl-1,2,4,5-tetrazin-3-yl)-L-phenylalanine.  A phenylalanine
scaffold with a para 1,2,4,5-tetrazine bearing an ethyl at the distal carbon --
not a lysine carbamate, and no PEG in the amino acid itself.  Residue mass
255.27 Da, confirmed twice in Eddins et al., Bio-protocol 14(16):e5048 (2024):
the 100 mM stock recipe implies 309 g/mol for the HCl salt (p.8), and the
sfGFP ESI-MS shift is +141 Da for Asn -> Tet2-Et (p.21).  Encoded at TAG by
Mj-TyrRS(E7) from plasmid pAJE3-E7.  Built in the *oxidized* aromatic tetrazine
state, which is what exists after aerobic purification (p.20).

sTCO = strained trans-cyclooctene, i.e. the cyclopropane-fused
trans-bicyclo[6.1.0]non-4-ene, ~100x faster than plain TCO (Eddins p.19).
Numbering: C1..C8 form the eight-membered ring, C9 bridges C1-C8 as the
cyclopropane, and the alkene sits at C4=C5.  C1-C8 are bonded directly (the
"0" bridge of bicyclo[6.1.0]), so C1/C8/C9 are the cyclopropane.

The IEDDA forms bonds at the two tetrazine carbons and the retro-Diels-Alder
then expels one of the two diaza bridges as N2, leaving a
**4,5-dihydropyridazine**: a conjugated C=N-N=C azadiene with no N-H, fused
across its C4-C5 to the former alkene carbons.  Aromatisation to pyridazine is
blocked because those two carbons are sp3 ring-fusion centres of the
bicyclononane.  This is the structure drawn in Eddins Fig. 1B (p.3), the
graphical overview (p.2), and slide 21 of the lab meeting deck.

The DNA-side arm is a design decision, not literature: no source in the reading
list says how sTCO reaches DNA (the only characterised reagent is sTCO-PEG5000,
a gel-shift diagnostic).  We use the standard route from sTCO's own precursor,
(1R,8S,9S)-bicyclo[6.1.0]non-4-en-9-yl)methanol, activated as an NHS carbonate
and reacted with a 5'-amino-modifier-C6 oligo, giving a **carbamate**:

    ...C9(H)-CH2-O-C(=O)-NH-(CH2)6-O-PO2(-)-O5'-DNA
                        |
                        cut here: TDP | DNL

Both residue cuts sit at amide-like bonds (the carbamate C-N, and the
carbamate's own ester oxygen stays with TDP), never through the conjugated ring
system.  DNL's bridging oxygen takes the structural role of the preceding
nucleotide's O3', so nucleotide 1 uses the *internal* DA/DT/DG/DC template and
keeps its -1 phosphate: no phosphate charge surgery is needed, and both custom
residues are neutral.

Stereochemistry
---------------
sTCO is (1R,8S,9S).  The cycloaddition can occur on either tetrazine face, so
two diastereomers form; the literature does not resolve them.  We build one and
say so.
"""
from __future__ import annotations

# Element of each atom is inferred from the leading character of its name,
# except where the name would lie: listed explicitly where needed.
ELEMENT_OVERRIDE: dict[str, str] = {}

# ---------------------------------------------------------------------------
# TET -- Tet2-Et, oxidized tetrazine
# ---------------------------------------------------------------------------
# bonds are (a, b, order); order 1 = single, 2 = double, 'ar' = aromatic
TET = {
    "name": "TET",
    "description": "Tet2-Et: 4-(6-ethyl-1,2,4,5-tetrazin-3-yl)-L-phenylalanine, oxidized",
    "smiles": r"N[C@@H](Cc1ccc(cc1)-c1nnc(CC)nn1)C(=O)O",
    "net_charge": 0,
    "expected_formula": "C13H13N5O",     # in-chain residue: no OXT, one backbone H
    "expected_residue_mass": 255.27,
    "head": "N", "tail": "C", "main_chain": ["N", "CA", "C"],
    "heavy": [
        # backbone
        "N", "CA", "C", "O",
        # CH2 and para-phenylene
        "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
        # 1,2,4,5-tetrazine: ring order NT1-NT2-CT3-NT4-NT5-CT6-NT1
        "CT3", "NT2", "NT1", "CT6", "NT5", "NT4",
        # ethyl on CT6
        "CX1", "CX2",
    ],
    "hydrogens": {
        "N": ["H"], "CA": ["HA"], "CB": ["HB2", "HB3"],
        "CD1": ["HD1"], "CD2": ["HD2"], "CE1": ["HE1"], "CE2": ["HE2"],
        "CX1": ["HX1", "HX2"], "CX2": ["HX3", "HX4", "HX5"],
    },
    "bonds": [
        ("N", "CA", 1), ("CA", "C", 1), ("C", "O", 2), ("CA", "CB", 1),
        ("CB", "CG", 1),
        ("CG", "CD1", "ar"), ("CD1", "CE1", "ar"), ("CE1", "CZ", "ar"),
        ("CZ", "CE2", "ar"), ("CE2", "CD2", "ar"), ("CD2", "CG", "ar"),
        ("CZ", "CT3", 1),
        ("CT3", "NT2", "ar"), ("NT2", "NT1", "ar"), ("NT1", "CT6", "ar"),
        ("CT6", "NT5", "ar"), ("NT5", "NT4", "ar"), ("NT4", "CT3", "ar"),
        ("CT6", "CX1", 1), ("CX1", "CX2", 1),
    ],
}

# ---------------------------------------------------------------------------
# TDP -- the click adduct, protein side, ending at the carbamate carbonyl
# ---------------------------------------------------------------------------
TDP = {
    "name": "TDP",
    "description": ("Tet2-Et + sTCO click adduct: 4,5-dihydropyridazine fused to "
                    "trans-bicyclo[6.1.0]nonane, ending at the carbamate C=O"),
    # R on the carbamate nitrogen is DNL; the capped model compound methylates it.
    "smiles_capped": (r"N[C@@H](Cc1ccc(cc1)C2=NN=C(CC)[C@@H]3CC[C@H]4[C@@H]"
                      r"(COC(=O)NC)[C@H]4CC[C@H]23)C(=O)O"),
    "net_charge": 0,
    "head": "N", "tail": "C", "main_chain": ["N", "CA", "C"],
    # the extra covalent bond tleap must make, beyond head/tail:
    "extra_bond": {"this_atom": "CN", "partner_residue": "DNL", "partner_atom": "N"},
    "heavy": [
        "N", "CA", "C", "O",
        "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
        # 4,5-dihydropyridazine: CP3(=NP2)-NP1(=CP6); CQ4/CQ5 sp3 ring fusion
        "CP3", "NP2", "NP1", "CP6", "CQ5", "CQ4",
        # ethyl on CP6
        "CX1", "CX2",
        # bicyclo[6.1.0]nonane remnant: eight-ring CO1-CO2-CO3-CQ4-CQ5-CO6-CO7-CO8-CO1
        # with cyclopropane CO1-CO8-CO9
        "CO1", "CO2", "CO3", "CO6", "CO7", "CO8", "CO9",
        # carbamate arm
        "CM", "OM", "CN", "ON",
    ],
    "hydrogens": {
        "N": ["H"], "CA": ["HA"], "CB": ["HB2", "HB3"],
        "CD1": ["HD1"], "CD2": ["HD2"], "CE1": ["HE1"], "CE2": ["HE2"],
        "CX1": ["HX1", "HX2"], "CX2": ["HX3", "HX4", "HX5"],
        "CQ4": ["HQ4"], "CQ5": ["HQ5"],
        "CO1": ["HO1"], "CO2": ["HO2A", "HO2B"], "CO3": ["HO3A", "HO3B"],
        "CO6": ["HO6A", "HO6B"], "CO7": ["HO7A", "HO7B"],
        "CO8": ["HO8"], "CO9": ["HO9"],
        "CM": ["HM1", "HM2"],
    },
    "bonds": [
        ("N", "CA", 1), ("CA", "C", 1), ("C", "O", 2), ("CA", "CB", 1),
        ("CB", "CG", 1),
        ("CG", "CD1", "ar"), ("CD1", "CE1", "ar"), ("CE1", "CZ", "ar"),
        ("CZ", "CE2", "ar"), ("CE2", "CD2", "ar"), ("CD2", "CG", "ar"),
        ("CZ", "CP3", 1),
        # the azadiene: C=N-N=C, no N-H  (4,5-dihydropyridazine)
        ("CP3", "NP2", 2), ("NP2", "NP1", 1), ("NP1", "CP6", 2),
        ("CP6", "CQ5", 1), ("CQ5", "CQ4", 1), ("CQ4", "CP3", 1),
        ("CP6", "CX1", 1), ("CX1", "CX2", 1),
        # eight-membered carbocycle
        ("CQ4", "CO3", 1), ("CO3", "CO2", 1), ("CO2", "CO1", 1),
        ("CQ5", "CO6", 1), ("CO6", "CO7", 1), ("CO7", "CO8", 1),
        ("CO1", "CO8", 1),
        # cyclopropane
        ("CO1", "CO9", 1), ("CO8", "CO9", 1),
        # carbamate arm off the cyclopropane methine
        ("CO9", "CM", 1), ("CM", "OM", 1), ("OM", "CN", 1), ("CN", "ON", 2),
    ],
}

# ---------------------------------------------------------------------------
# DNL -- 5'-amino-modifier-C6 linker, DNA side
# ---------------------------------------------------------------------------
DNL = {
    "name": "DNL",
    "description": "5'-amino-modifier-C6 arm: carbamate N through hexyl to the "
                   "oxygen that bridges the first nucleotide's phosphorus",
    "smiles_capped": r"CC(=O)NCCCCCCO[P](=O)([O-])OC",
    # -0.3079, not 0.  Amber's DNA libraries split one unit of charge between the two
    # chain ends: every 3'-terminal template (DA3/DG3/DC3/DT3) carries -0.6921 and
    # every 5'-terminal one -0.3079, so that a complete strand comes out integral.
    # In this chimera the 5' nucleotide uses the *internal* template -- it has to, to
    # keep the phosphate the tether bonds to -- and DNL stands in its place.  DNL must
    # therefore carry the 5'-cap's share, or the strand is left 0.3079 e short and the
    # whole system has a fractional net charge that no ion count can neutralise.
    # -0.3079 + (-0.6921) = -1 exactly, as DA5 + DA3 does.
    "net_charge": -0.3079,
    "head": "N", "tail": "OL",
    "extra_bond": {"this_atom": "OL", "partner_residue": "DNA_5prime",
                   "partner_atom": "P"},
    # The hexyl carbons are C1..C6, not CL1..CL6.  antechamber (and PyMOL, and VMD)
    # infer an atom's element from its name, and a name beginning "CL" reads as
    # chlorine: an earlier version of this table produced six chlorines carrying two
    # hydrogens each, with GAFF type "cl", and parmchk2 duly flagged 14 invented
    # parameters for bonds that do not exist.  Never name a carbon CL-anything.
    "heavy": ["N", "C1", "C2", "C3", "C4", "C5", "C6", "OL"],
    "hydrogens": {
        "N": ["HN"],
        "C1": ["H11", "H12"], "C2": ["H21", "H22"], "C3": ["H31", "H32"],
        "C4": ["H41", "H42"], "C5": ["H51", "H52"], "C6": ["H61", "H62"],
    },
    "bonds": [
        ("N", "C1", 1), ("C1", "C2", 1), ("C2", "C3", 1), ("C3", "C4", 1),
        ("C4", "C5", 1), ("C5", "C6", 1), ("C6", "OL", 1),
    ],
}

RESIDUES = {r["name"]: r for r in (TET, TDP, DNL)}

# Standard atomic masses, enough for the elements present.
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974}
# Expected heavy-atom valences.  Aromatic bonds count as 1.5 for the check.
VALENCE = {"C": 4, "N": 3, "O": 2, "P": 5}


def element_of(name: str) -> str:
    return ELEMENT_OVERRIDE.get(name, name[0])


def formula(res: dict) -> str:
    counts: dict[str, int] = {}
    for a in res["heavy"]:
        counts[element_of(a)] = counts.get(element_of(a), 0) + 1
    nh = sum(len(v) for v in res["hydrogens"].values())
    counts["H"] = counts.get("H", 0) + nh
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}"
                   for e in ("C", "H", "N", "O", "P") if counts.get(e))


def mass(res: dict) -> float:
    m = sum(MASS[element_of(a)] for a in res["heavy"])
    m += sum(MASS["H"] * len(v) for v in res["hydrogens"].values())
    return m


def validate(res: dict) -> list[str]:
    """Check the tables are internally consistent: names unique, bonds refer to
    declared atoms, and every heavy atom's valence is satisfied once its declared
    hydrogens and its inter-residue connections are counted."""
    errs = []
    heavy = res["heavy"]
    if len(set(heavy)) != len(heavy):
        errs.append("duplicate heavy-atom names")
    hnames = [h for v in res["hydrogens"].values() for h in v]
    if len(set(hnames)) != len(hnames):
        dupes = {h for h in hnames if hnames.count(h) > 1}
        errs.append(f"duplicate hydrogen names: {sorted(dupes)}")
    if set(res["hydrogens"]) - set(heavy):
        errs.append(f"hydrogens on undeclared atoms: {set(res['hydrogens']) - set(heavy)}")

    order = {1: 1.0, 2: 2.0, 3: 3.0, "ar": 1.5}
    val: dict[str, float] = {a: 0.0 for a in heavy}
    for a, b, o in res["bonds"]:
        for x in (a, b):
            if x not in val:
                errs.append(f"bond refers to undeclared atom {x}")
        if a in val and b in val:
            val[a] += order[o]
            val[b] += order[o]
    # hydrogens
    for a, hs in res["hydrogens"].items():
        val[a] += len(hs)
    # inter-residue connections: head takes one bond from the preceding residue,
    # tail one from the following, and any declared extra bond one more.
    external = {res["head"]: 1, res.get("tail", ""): 1}
    # An extra_bond is only *extra* if it is not already the head or tail connection.
    # For DNL the tail bond and the bond to the nucleotide phosphorus are the same
    # bond, so counting both would give the bridging oxygen a spurious third bond.
    if "extra_bond" in res:
        ea = res["extra_bond"]["this_atom"]
        if ea not in (res["head"], res.get("tail")):
            external[ea] = external.get(ea, 0) + 1
    for a, n in external.items():
        if a in val:
            val[a] += n

    for a in heavy:
        want = VALENCE[element_of(a)]
        got = val[a]
        # aromatic rings give x.5 sums; accept within 0.6 of the target
        if abs(got - want) > 0.6:
            errs.append(f"{a} ({element_of(a)}) valence {got} != {want}")
    return errs


if __name__ == "__main__":
    import sys
    bad = 0
    for name, res in RESIDUES.items():
        errs = validate(res)
        f, m = formula(res), mass(res)
        print(f"{name}: {f}  mass {m:.2f} Da  "
              f"{len(res['heavy'])} heavy + "
              f"{sum(len(v) for v in res['hydrogens'].values())} H")
        print(f"     {res['description']}")
        if "expected_formula" in res:
            ok = f == res["expected_formula"]
            print(f"     formula check vs literature {res['expected_formula']}: "
                  f"{'OK' if ok else 'MISMATCH'}")
            bad += not ok
        if "expected_residue_mass" in res:
            d = abs(m - res["expected_residue_mass"])
            print(f"     mass check vs literature {res['expected_residue_mass']} Da: "
                  f"{'OK' if d < 0.1 else 'MISMATCH'} (diff {d:.2f})")
            bad += d >= 0.1
        if errs:
            bad += len(errs)
            for e in errs:
                print(f"     ERROR: {e}")
        else:
            print("     valence/connectivity: OK")
        print()

    # Mass bookkeeping across the click reaction: the IEDDA is an addition followed by
    # retro-Diels-Alder loss of N2, so
    #     TDP + DNL  ==  TET + reagent_arm - N2
    # where reagent_arm is everything the sTCO carbamate contributes.  The carbamate
    # nitrogen sits in DNL, so the TDP-side arm is C11H15O2 and the DNL side C6H13NO.
    m_n2 = 2 * MASS["N"]
    arm_tdp = mass(TDP) - (mass(TET) - m_n2)
    want_tdp = 11 * MASS["C"] + 15 * MASS["H"] + 2 * MASS["O"]      # C11H15O2
    want_dnl = 6 * MASS["C"] + 13 * MASS["H"] + MASS["N"] + MASS["O"]  # C6H13NO
    print(f"TDP - (TET - N2) = {arm_tdp:.2f} Da   expected C11H15O2 = {want_tdp:.2f}  "
          f"{'OK' if abs(arm_tdp - want_tdp) < 0.05 else 'MISMATCH'}")
    print(f"DNL              = {mass(DNL):.2f} Da   expected C6H13NO  = {want_dnl:.2f}  "
          f"{'OK' if abs(mass(DNL) - want_dnl) < 0.05 else 'MISMATCH'}")
    print(f"whole tether TDP+DNL = {mass(TDP) + mass(DNL):.2f} Da "
          f"(protein CB through to the bridging O on the DNA phosphate)")
    bad += abs(arm_tdp - want_tdp) >= 0.05
    bad += abs(mass(DNL) - want_dnl) >= 0.05
    print(f"\n{'FAILED: ' + str(bad) + ' problem(s)' if bad else 'all definitions consistent'}")
    sys.exit(1 if bad else 0)
