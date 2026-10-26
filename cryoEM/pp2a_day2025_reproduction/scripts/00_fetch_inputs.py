#!/usr/bin/env python
"""Plan section 4/5 - fetch every public input once, with checksums.

Cryo-EM: EMD-48300 (Aalpha-B56g-Ca, 3.4 A) and EMD-48224 (Abeta-B56g-Ca, 3.2 A):
primary map, both half maps, deposited FSC XML.  Models 9MIP, 9MF5, start model 3FGA.
MD: AlphaFold Server (AF3) CIFs if Ian placed them in inputs/af3/, else AlphaFold DB
(AF2, v4) models for P30153 (Aalpha) and P30154 (Abeta).  UniProt sequences for checks.

Everything lands on netscratch ($PP2A_DATA/inputs, em/emdb, em/models) - bulk and
regenerable.  Only the small AF models are also copied into the repo inputs/ dir.
Idempotent: existing files with a recorded md5 are skipped.
"""
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
DATA = os.environ.get("PP2A_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025")
EMDB = os.path.join(DATA, "em", "emdb")
MODELS = os.path.join(DATA, "em", "models")
AF = os.path.join(DATA, "inputs", "af")
MANIFEST = os.path.join(DATA, "inputs", "manifest.json")

ENTRIES = {"Aa": "48300", "Ab": "48224"}
PDBS = {"Aa": "9MIP", "Ab": "9MF5", "start": "3FGA"}
UNIPROT = {"Aa": "P30153", "Ab": "P30154"}
EMDB_FTP = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{id}/"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, manifest):
    if os.path.exists(dest) and dest in manifest:
        print(f"  have {os.path.basename(dest)}")
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    print(f"  get  {url}")
    if shutil.which("aria2c"):
        subprocess.run(["aria2c", "-q", "-x4", "-s4", "--allow-overwrite=true",
                        "-d", os.path.dirname(dest), "-o", os.path.basename(dest), url], check=True)
    else:
        urllib.request.urlretrieve(url, dest)
    manifest[dest] = {"url": url, "md5": md5(dest), "bytes": os.path.getsize(dest)}


def gunzip(src):
    dst = src[:-3]
    if not os.path.exists(dst):
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 22)
    return dst


def main():
    manifest = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else {}
    for iso, eid in ENTRIES.items():
        print(f"EMD-{eid} ({iso})")
        base = EMDB_FTP.format(id=eid)
        for rel in (f"map/emd_{eid}.map.gz", f"other/emd_{eid}_half_map_1.map.gz",
                    f"other/emd_{eid}_half_map_2.map.gz", f"fsc/emd_{eid}_fsc.xml"):
            dest = os.path.join(EMDB, iso, os.path.basename(rel))
            fetch(base + rel, dest, manifest)
            if dest.endswith(".gz"):
                gunzip(dest)
        # EMDB header XML (pixel size, contour level) for the record
        fetch(base + f"header/emd-{eid}-v30.xml", os.path.join(EMDB, iso, f"emd-{eid}-v30.xml"), manifest)
    print("models")
    for tag, pid in PDBS.items():
        for ext in ("cif", "pdb"):
            fetch(f"https://files.rcsb.org/download/{pid}.{ext}", os.path.join(MODELS, f"{pid}.{ext}"), manifest)
    print("sequences")
    for iso, acc in UNIPROT.items():
        fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", os.path.join(AF, f"{acc}.fasta"), manifest)
    print("AlphaFold models")
    af3 = sorted(f for f in os.listdir(os.path.join(REPO, "inputs", "af3")) if f.endswith((".cif", ".pdb")))
    source = "AlphaFold Server (AF3), supplied by Ian" if af3 else "AlphaFold DB (AF2)"
    if af3:
        for f in af3:
            shutil.copy(os.path.join(REPO, "inputs", "af3", f), os.path.join(AF, f))
            print(f"  AF3 model {f}")
    else:
        # AlphaFold DB re-versions its files (v4 -> v6 in 2025); ask the API for the current one.
        for iso, acc in UNIPROT.items():
            with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}") as r:
                entry = json.load(r)[0]
            for key in ("pdbUrl", "cifUrl"):
                url = entry[key]
                fetch(url, os.path.join(AF, os.path.basename(url)), manifest)
            manifest[f"_afdb_{acc}"] = {"latestVersion": entry.get("latestVersion"),
                                        "modelCreatedDate": entry.get("modelCreatedDate")}
        source += f" (AFDB version {entry.get('latestVersion')})"
    manifest["_af_source"] = source
    os.makedirs(os.path.dirname(MANIFEST), exist_ok=True)
    json.dump(manifest, open(MANIFEST, "w"), indent=2)
    total = sum(v["bytes"] for k, v in manifest.items() if not k.startswith("_"))
    print(f"done: {len(manifest) - 1} files, {total / 1e6:.0f} MB, AF source = {source}")


if __name__ == "__main__":
    sys.exit(main())
