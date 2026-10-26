#!/usr/bin/env python
"""Plan section 4 - starting structures for the MD (Tier 1: full-length AlphaFold models).

Takes the AlphaFold models fetched by 00_fetch_inputs.py (AF3 CIFs from Ian if present, else the
AlphaFold DB models), checks the sequence against UniProt (Aalpha P30153 = 589 aa; Abeta P30154
= 601 aa with the 12-residue N-terminal extension MAGASELGTGPG), strips hydrogens and
non-protein atoms, and writes clean single-chain PDBs for tleap to $PP2A_DATA/md/systems/<sys>/protein.pdb.
Also records per-residue pLDDT of the N-terminus (1-13 / 1-25) - the region whose dynamics the paper
studies - in prepare_models.json.
"""
import json
import os
import sys

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio import SeqIO

DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
AF = os.path.join(DATA, "inputs", "af")
SYS = os.path.join(DATA, "md", "systems")
SYSTEMS = {"Aa": {"acc": "P30153", "length": 589, "nterm": 13},
           "Ab": {"acc": "P30154", "length": 601, "nterm": 25, "extension": "MAGASELGTGPG"}}


class ProteinHeavy(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "

    def accept_atom(self, atom):
        return atom.element != "H" and not atom.get_id().startswith("H")


def find_model(acc):
    cands = sorted(f for f in os.listdir(AF) if acc in f and f.endswith((".cif", ".pdb")))
    af3 = [f for f in cands if "AF-" not in f]          # user-supplied AlphaFold Server files
    afdb = [f for f in cands if f.startswith("AF-")]
    if af3:
        return os.path.join(AF, af3[0]), "AlphaFold Server (AF3)"
    pdbs = [f for f in afdb if f.endswith(".pdb")]
    return os.path.join(AF, (pdbs or afdb)[0]), "AlphaFold DB (AF2)"


def load(path):
    parser = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return parser.get_structure("m", path)


def main():
    report = {}
    for name, meta in SYSTEMS.items():
        path, source = find_model(meta["acc"])
        st = load(path)
        chain = next(st[0].get_chains())
        residues = [r for r in chain if r.id[0] == " "]
        seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in residues)
        uni = str(next(SeqIO.parse(os.path.join(AF, f"{meta['acc']}.fasta"), "fasta")).seq)
        assert len(uni) == meta["length"], (name, len(uni))
        assert seq == uni, f"{name}: model sequence differs from UniProt"
        if "extension" in meta:
            assert seq.startswith(meta["extension"]), seq[:12]
        assert residues[0].id[1] == 1, "numbering must start at Met1"
        # pLDDT lives in the B-factor column of AlphaFold files
        plddt = [float(next(iter(r)).get_bfactor()) for r in residues]
        nt = meta["nterm"]
        out = os.path.join(SYS, name)
        os.makedirs(out, exist_ok=True)
        io = PDBIO()
        io.set_structure(st)
        io.save(os.path.join(out, "protein.pdb"), ProteinHeavy())
        report[name] = {"source_file": os.path.basename(path), "source": source, "n_res": len(residues),
                        "plddt_mean_all": sum(plddt) / len(plddt),
                        f"plddt_mean_nterm_1_{nt}": sum(plddt[:nt]) / nt,
                        "plddt_min_nterm": min(plddt[:nt]), "sequence_ok": True}
        print(f"{name}: {source} {os.path.basename(path)}  {len(residues)} res  "
              f"pLDDT all {report[name]['plddt_mean_all']:.1f}, N-term 1-{nt} {report[name][f'plddt_mean_nterm_1_{nt}']:.1f}")
    json.dump(report, open(os.path.join(SYS, "prepare_models.json"), "w"), indent=2)


if __name__ == "__main__":
    sys.exit(main())
