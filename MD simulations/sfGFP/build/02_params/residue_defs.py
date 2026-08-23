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
keeps its -1 phosphate.  The charges are not all neutral: DNL carries -0.3079 and
DNP -0.6921, the two shares Amber splits between a strand's 5' and 3' ends.

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

# ---------------------------------------------------------------------------
# DNH -- the DNA-side arm, capped as a free alcohol instead of joined to DNA
# ---------------------------------------------------------------------------
# This exists to make one comparison clean.  The force-clamp ladder runs on TET
# (unclicked tetrazine), so it cannot serve as the zero-force control for a chimera
# built from TDP: the two differ by the whole sTCO adduct, and that difference turns
# out to move the attachment sites about as far as the spring does.  DNH gives a
# system with the identical tether chemistry -- TDP plus the full carbamate arm --
# and no DNA and no load, which is what isolates the spring's force from the
# tether's bulk.
DNH = {
    "name": "DNH",
    "description": "5'-amino-C6 arm terminated as a free alcohol: the clicked, "
                   "unloaded control for the chimeras",
    "smiles_capped": r"CC(=O)NCCCCCCO",
    "net_charge": 0,
    "head": "N", "tail": None,
    "heavy": ["N", "C1", "C2", "C3", "C4", "C5", "C6", "OL"],
    "hydrogens": {
        "N": ["HN"], "OL": ["HOL"],
        "C1": ["H11", "H12"], "C2": ["H21", "H22"], "C3": ["H31", "H32"],
        "C4": ["H41", "H42"], "C5": ["H51", "H52"], "C6": ["H61", "H62"],
    },
    "bonds": [
        ("N", "C1", 1), ("C1", "C2", 1), ("C2", "C3", 1), ("C3", "C4", 1),
        ("C4", "C5", 1), ("C5", "C6", 1), ("C6", "OL", 1),
    ],
}

# ---------------------------------------------------------------------------
# DNP -- the 3'-amino-modifier arm, for the corrected same-strand attachment
# ---------------------------------------------------------------------------
# The experiment modifies ONE strand at both ends and leaves the complement
# untouched; that is what closes the covalent loop and what stops one chimera
# base-pairing with another to form the dimers Zocchi showed are favourable.  Our
# first build put a linker on the 5' end of each of the two strands, which loads the
# duplex mechanically but leaves the "loop" closed only by base pairing -- every
# tether bond is a bridge in the bond graph, not part of a cycle.
#
# A 3' modifier is not a mirror image of DNL: it brings its own phosphate, because
# it is installed as a phosphoramidite on the 3'-OH.  So unlike DNL it is charged,
# and its head is the phosphorus rather than the nitrogen.
#
#     nt(N)-O3' -- P(=O)(O-) -- O -- (CH2)7 -- NH -- C(=O)... [TDP side]
#                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#                  this residue
#
# Net charge -0.6921, and DNL keeps its -0.3079.  The modified strand carries 28
# phosphates (nucleotide 1's, the 26 internal ones, and this residue's), so its
# residues must sum to -28:  -0.3079 + 27*(-1) + q = -28  =>  q = -0.6921.  That is
# exactly the share Amber's 3'-terminal DA3/DG3/DC3/DT3 templates carry, which is
# fitting -- this residue is what replaces them.
#
# Named DNP, not DN3: three places in the analysis do name.rstrip("35") to strip
# Amber's terminal-variant suffixes, so "DN3" would become "DN" and fall out of the
# nucleic-acid lists, and build_system's compatible() treats X3/X5 as a legal tleap
# rename.  A trailing digit here is a silent-failure waiting to happen.
DNP = {
    "name": "DNP",
    "description": "3'-amino-modifier C7 arm: its own phosphate on the strand's 3' "
                   "oxygen, through a heptyl spacer to the carbamate nitrogen",
    "smiles_capped": r"COP(=O)([O-])OCCCCCCCNC(C)=O",
    "net_charge": -0.6921,
    "head": "P", "tail": None,
    "formal_charges": {"OP2": -1},
    "extra_bond": {"this_atom": "N", "partner_residue": "TDP",
                   "partner_atom": "CN"},
    # main_chain is given explicitly because prepgen is only handed MAIN_CHAIN lines
    # when a residue has a tail, and this one does not.  Without them prepgen guesses
    # the tree, and for a chain hanging off a four-substituent phosphorus that guess
    # is not safe -- the same condition that once made DNL's bridging oxygen a child
    # of its nitrogen.
    "main_chain": ["OL", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "N"],
    "heavy": ["P", "OP1", "OP2", "OL",
              "C1", "C2", "C3", "C4", "C5", "C6", "C7", "N"],
    "hydrogens": {
        "N": ["HN"],
        "C1": ["H11", "H12"], "C2": ["H21", "H22"], "C3": ["H31", "H32"],
        "C4": ["H41", "H42"], "C5": ["H51", "H52"], "C6": ["H61", "H62"],
        "C7": ["H71", "H72"],
    },
    "bonds": [
        ("P", "OP1", 2), ("P", "OP2", 1), ("P", "OL", 1), ("OL", "C1", 1),
        ("C1", "C2", 1), ("C2", "C3", 1), ("C3", "C4", 1), ("C4", "C5", 1),
        ("C5", "C6", 1), ("C6", "C7", 1), ("C7", "N", 1),
    ],
}

# ---------------------------------------------------------------------------
# The ubiquitin chimera: a different protein, a different chemistry, the same DNA
# ---------------------------------------------------------------------------
# Boral, Schnebly et al. attach DNA to ubiquitin with DBCO-maleimide (Vector Labs
# CCT-A108, CAS 1395786-30-7, C25H21N3O4, MW 427.75), which forms two bonds: a
# thioether from a surface cysteine's S-gamma across the maleimide, and a triazole from
# the strained cyclooctyne to an azide on the oligo.  The whole covalent path from the
# protein to the first nucleotide is
#
#   Cys C-beta - S - [succinimide] - CH2CH2 - NH - CO - CH2CH2 - CO - N[DBCO]
#     ... DBCO cyclooctyne fused to [1,2,3-triazole] ... - CH2CH2 - CO - NH - [DNL/DNP]
#
# and it is split into two residues at the first amide, CMS + CDB, for a practical
# reason: the whole arm is 45 heavy atoms, which is a large and slow AM1-BCC problem,
# whereas 16 + 29 is two comfortable ones.  The split is free because it falls on an
# amide -- a c-ns bond, exactly the bond type junction.frcmod already carries for the
# carbamate in the sfGFP tether, so no new cross terms are needed.
#
# It also means DNL and DNP are reused unchanged.  Their nitrogens are amidated by
# whatever carbonyl precedes them; in sfGFP that is TDP's carbamate and here it is the
# NHS-ester-derived amide, and the bond is the same in both.
#
# Four structural details are CHOSEN rather than established, and each is flagged in
# the report as an approximation:
#
#   1. The spacer lengths inside DBCO-maleimide.  The vendor gives a formula and a "C6
#      spacer" but not a structure sheet.  Only one arrangement of the known pieces --
#      azadibenzocyclooctyne, a 4-oxobutanoyl, an amide, an ethylene, a maleimide --
#      reproduces C25H21N3O4 exactly, and that is what is built here:
#        C15H10N (azocinyl) + C4H4O2 + NH + C2H4 + C4H2NO2 = C25H21N3O4.
#      Getting a CH2 wrong here would change the arm's reach by ~1.2 A on a linker that
#      is already ~20 A of freely rotating chain, so the mechanical consequence is small,
#      but it is an assumption and not a measurement.
#   2. Triazole regiochemistry.  Strain-promoted cycloaddition on DBCO gives both
#      orientations and they are not separated; one is built.
#   3. The stereocentre created at the thioether carbon.  Racemic in the experiment; one
#      epimer is built.
#   4. The NHS-ester azide's own spacer, taken as 3-azidopropanoyl.
#
# None of the four is load-bearing for the question the ubiquitin runs are meant to
# answer, which is whether local protein-DNA contacts or spring action dominates.
CMS = {
    "name": "CMS",
    "description": "cysteine with its S-gamma added across a maleimide: the "
                   "succinimidyl thioether, through an ethylene to the amide nitrogen "
                   "that continues into the DBCO half",
    "smiles_capped": r"CC(=O)NC(CSC1CC(=O)N(CCNC(C)=O)C1=O)C(=O)NC",
    "net_charge": 0.0,
    "head": "N", "tail": "C", "aa": True,
    "extra_bond": {"this_atom": "N7", "partner_residue": "CDB",
                   "partner_atom": "CW"},
    "heavy": ["N", "CA", "C", "O", "CB", "SG",
              "NS", "CS1", "OS1", "CS2", "CS3", "CS4", "OS4", "CV5", "CV6", "N7"],
    "hydrogens": {
        "N": ["H"], "CA": ["HA"], "CB": ["HB2", "HB3"],
        "CS2": ["HS2"], "CS3": ["HS31", "HS32"],
        "CV5": ["HV51", "HV52"], "CV6": ["HV61", "HV62"], "N7": ["HN7"],
    },
    "bonds": [
        ("N", "CA", 1), ("CA", "C", 1), ("C", "O", 2), ("CA", "CB", 1),
        ("CB", "SG", 1), ("SG", "CS2", 1),
        # the succinimide ring, saturated by the thiol addition
        ("NS", "CS1", 1), ("CS1", "OS1", 2), ("CS1", "CS2", 1), ("CS2", "CS3", 1),
        ("CS3", "CS4", 1), ("CS4", "OS4", 2), ("CS4", "NS", 1),
        # out of the ring nitrogen to the amide that joins CDB
        ("NS", "CV5", 1), ("CV5", "CV6", 1), ("CV6", "N7", 1),
    ],
}

CDB = {
    "name": "CDB",
    "description": "the DBCO half: an amide and a 4-oxobutanoyl into the "
                   "azadibenzocyclooctyne, whose strained alkyne has become a "
                   "1,2,3-triazole with the oligo's azide, then out to the carbonyl "
                   "that amidates DNL or DNP",
    "smiles_capped": r"CNC(=O)CCC(=O)N1Cc2ccccc2c3nnn(CCC(=O)NC)c3c4ccccc41",
    "net_charge": 0.0,
    "head": "CW", "tail": None,
    "extra_bond": {"this_atom": "CG", "partner_residue": "DNL",
                   "partner_atom": "N"},
    # Stated explicitly because CDB has no tail: prepgen is handed no MAIN_CHAIN in that
    # case and would guess the tree through a fused polycyclic cage.
    # Amber's MAIN_CHAIN lists the atoms BETWEEN head and tail, so the head CW is not
    # named here (compare DNP, whose list starts after its head phosphorus).  The path
    # runs through the fused cage by way of one benzo ring and the triazole.
    "main_chain": ["CV1", "CV2", "CT", "NZ", "CP1", "CP2", "CJ2", "NJ1",
                   "CY1", "CY2", "CG"],
    "heavy": ["CW", "OW", "CV1", "CV2", "CT", "OT", "NZ", "CK",
              "CP1", "CP2", "CP3", "CP4", "CP5", "CP6",
              "CQ1", "CQ2", "CQ3", "CQ4", "CQ5", "CQ6",
              "CJ1", "CJ2", "NJ1", "NJ2", "NJ3",
              "CY1", "CY2", "CG", "OG"],
    "hydrogens": {
        "CV1": ["HV11", "HV12"], "CV2": ["HV21", "HV22"], "CK": ["HK1", "HK2"],
        "CP3": ["HP3"], "CP4": ["HP4"], "CP5": ["HP5"], "CP6": ["HP6"],
        "CQ3": ["HQ3"], "CQ4": ["HQ4"], "CQ5": ["HQ5"], "CQ6": ["HQ6"],
        "CY1": ["HY11", "HY12"], "CY2": ["HY21", "HY22"],
    },
    "bonds": [
        ("CW", "OW", 2), ("CW", "CV1", 1), ("CV1", "CV2", 1),
        ("CV2", "CT", 1), ("CT", "OT", 2), ("CT", "NZ", 1),
        # the eight-membered azocine: N and a methylene on one side, the former alkyne
        # (now two triazole carbons) on the other, with a benzo ring bridging each pair
        ("NZ", "CP1", 1), ("NZ", "CK", 1), ("CK", "CQ1", 1),
        ("CP1", "CP2", "ar"), ("CP2", "CP3", "ar"), ("CP3", "CP4", "ar"),
        ("CP4", "CP5", "ar"), ("CP5", "CP6", "ar"), ("CP6", "CP1", "ar"),
        ("CQ1", "CQ2", "ar"), ("CQ2", "CQ3", "ar"), ("CQ3", "CQ4", "ar"),
        ("CQ4", "CQ5", "ar"), ("CQ5", "CQ6", "ar"), ("CQ6", "CQ1", "ar"),
        ("CP2", "CJ2", 1), ("CQ2", "CJ1", 1),
        # the triazole formed by the cycloaddition; NJ1 is the nitrogen the azide
        # brought in, so it is the one carrying the chain back to the DNA
        ("NJ1", "NJ2", 1), ("NJ2", "NJ3", 2), ("NJ3", "CJ1", 1),
        ("CJ1", "CJ2", 2), ("CJ2", "NJ1", 1),
        ("NJ1", "CY1", 1), ("CY1", "CY2", 1), ("CY2", "CG", 1), ("CG", "OG", 2),
    ],
}

RESIDUES = {r["name"]: r for r in (TET, TDP, DNL, DNH, DNP, CMS, CDB)}

# Standard atomic masses, enough for the elements present.
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, "S": 32.06}
# Expected heavy-atom valences.  Aromatic bonds count as 1.5 for the check.
VALENCE = {"C": 4, "N": 3, "O": 2, "P": 5, "S": 2}


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
    external = {res["head"]: 1}
    if res.get("tail"):
        external[res["tail"]] = 1
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

    # A formally charged atom has a different valence: a carboxylate/phosphate O(-)
    # takes one bond, not two, and an ammonium N(+) takes four, not three.  Without
    # this the anionic phosphate oxygen of DNP fails the check by a whole bond, which
    # the +/-0.6 aromatic tolerance cannot absorb.
    formal = res.get("formal_charges", {})
    for a in heavy:
        want = VALENCE[element_of(a)] + formal.get(a, 0)
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

    # The ubiquitin crosslinker is a purchased reagent with a published formula, so the
    # split into CMS + CDB can be checked against it rather than merely asserted.  Drop
    # the atoms that do NOT come from DBCO-maleimide -- the cysteine's own backbone and
    # side chain, the three triazole nitrogens the azide brought in, and the arm running
    # from the triazole out to the DNA -- and what remains must be exactly C25H21N3O4.
    # One hydrogen moves in the reaction: the maleimide's C=C picks up the thiol's H.
    not_crosslinker = {"CMS": {"N", "CA", "C", "O", "CB", "SG"},
                       "CDB": {"NJ1", "NJ2", "NJ3", "CY1", "CY2", "CG", "OG"}}
    tally = {}
    for nm in ("CMS", "CDB"):
        r, drop = RESIDUES[nm], not_crosslinker[nm]
        for a in r["heavy"]:
            if a not in drop:
                tally[element_of(a)] = tally.get(element_of(a), 0) + 1
        for parent, hs in r["hydrogens"].items():
            if parent not in drop:
                tally["H"] = tally.get("H", 0) + len(hs) - (
                    1 if (nm == "CMS" and parent == "CS3") else 0)
    got = "".join(f"{e}{tally[e]}" for e in ("C", "H", "N", "O") if e in tally)
    ok = got == "C25H21N3O4"
    print(f"CMS+CDB minus protein, azide and DNA arm = {got}   "
          f"expected DBCO-maleimide C25H21N3O4   {'OK' if ok else 'MISMATCH'}")
    bad += (not ok)

    print(f"\n{'FAILED: ' + str(bad) + ' problem(s)' if bad else 'all definitions consistent'}")
    sys.exit(1 if bad else 0)
