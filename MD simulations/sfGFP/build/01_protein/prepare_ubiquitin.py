#!/usr/bin/env python3
"""1UBQ -> the F45W/T9C/K63C construct of Boral, Schnebly et al., ready for tleap.

Why this protein and these three mutations
------------------------------------------
Doeke asked for the ubiquitin case to be simulated alongside sfGFP, and it is the
better-controlled experiment of the two: the companion preprint from this lab measures
NMR chemical shifts on exactly this construct, with a perturbation ladder -- apo, then
crosslinker only, then two short strands, then ssDNA, then dsDNA -- and finds that most
of the effect has already appeared by the crosslinker stage.  That is a result a
simulation can speak to directly, unlike the beta-5 register shift, which is far slower
than MD can reach and which this build deliberately does not claim to address.

  T9C, K63C   the two attachment cysteines.  Ubiquitin has no native cysteine, so these
              are the only thiols and the labelling problem is clean.  The preprint
              chose them to straddle the beta-5 strand: 63 sits at its base and 9 is
              across the protein and not part of it.  In 1UBQ their C-beta are 24.8 A
              apart, which is the "~2.5 nm" the preprint quotes.
  F45W        a spectroscopic handle, not a mechanical one: it raises the UV absorbance
              about five-fold with little effect on the biophysics (Khorasanizadeh et
              al. 1993).  It is included because the experiment has it and the
              comparison should be against the molecule that was measured, but nothing
              here depends on it.

How the mutations are made
--------------------------
Side chains are truncated to the backbone plus C-beta and the residue renamed; tleap
then builds the new side chain from its own ff14SB template.  That is safe for the two
cysteines, which only gain an S-gamma and its hydrogen.  It is a real approximation for
F45W, which is buried and gains a fused indole where a phenyl was, so the built rotamer
gets a torsion scan against the rest of the protein rather than being accepted as tleap
places it -- and the result is checked for residual clashes before anything is written.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# The construct.  Values, not code, so the file states plainly what it builds.
MUTATIONS = {9: "CYS", 45: "TRP", 63: "CYS"}
ATTACHMENT_SITES = (9, 63)
KEEP_ON_MUTATION = ("N", "CA", "C", "O", "CB", "H", "HA")
WATER = {"HOH", "WAT"}


def read_pdb(path):
    atoms = []
    for ln in Path(path).read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        atoms.append({
            "record": ln[:6].strip(), "name": ln[12:16].strip(),
            "altloc": ln[16], "resn": ln[17:20].strip(), "chain": ln[21],
            "resi": int(ln[22:26]), "icode": ln[26],
            "xyz": np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]),
            "occ": float(ln[54:60] or 1.0),
            "elem": (ln[76:78].strip() or ln[12:16].strip()[0]),
        })
    return atoms


def fmt(serial, a, resn=None):
    nm = a["name"]
    fld = f" {nm:<3}" if len(nm) < 4 else nm
    return (f"ATOM  {serial:5d} {fld}{(resn or a['resn']):>4}{a['chain']:>2}"
            f"{a['resi']:4d}    "
            f"{a['xyz'][0]:8.3f}{a['xyz'][1]:8.3f}{a['xyz'][2]:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {a['elem']:>2}")


def his_tautomer(res, others):
    """Pick HID or HIE for a histidine from which ring nitrogen has a donor nearby.

    Same reasoning as the sfGFP preparation: the tautomer is a real choice that changes
    the hydrogen-bond pattern, and geometry decides it better than a default does.  If
    neither nitrogen has a partner within reach, HIE is used -- the more common neutral
    form -- and that is reported rather than hidden.
    """
    nd1 = res.get("ND1")
    ne2 = res.get("NE2")
    if nd1 is None or ne2 is None:
        return "HIE", "incomplete ring, defaulted"
    acceptors = [a for a in others
                 if a["elem"] in ("O", "N") and a["resn"] not in WATER]
    def nearest(n):
        d = [float(np.linalg.norm(a["xyz"] - n["xyz"])) for a in acceptors
             if float(np.linalg.norm(a["xyz"] - n["xyz"])) > 1.6]
        return min(d) if d else 99.0
    dn, de = nearest(nd1), nearest(ne2)
    # The nitrogen with a partner CLOSE BY is the one accepting, so the *other* one
    # carries the proton.
    if abs(dn - de) < 0.3:
        return "HIE", f"ND1 {dn:.2f} A, NE2 {de:.2f} A -- too close to call, HIE"
    if dn < de:
        return "HIE", f"ND1 has a partner at {dn:.2f} A (NE2 {de:.2f}), so NE2-H"
    return "HID", f"NE2 has a partner at {de:.2f} A (ND1 {dn:.2f}), so ND1-H"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default=str(HERE / "1ubq.pdb"))
    ap.add_argument("--out", default=str(HERE / "ubq_prepped.pdb"))
    ap.add_argument("--chain", default="A")
    ap.add_argument("--keep-waters", action="store_true", default=True)
    args = ap.parse_args()

    atoms = read_pdb(args.pdb)
    # 1UBQ is a single chain with no alternate locations, but assert it rather than
    # assume it: an altloc silently doubles a side chain.
    alt = {a["altloc"] for a in atoms} - {" ", ""}
    if alt:
        atoms = [a for a in atoms if a["altloc"] in (" ", "", "A")]
        print(f"  altlocs {sorted(alt)} present; kept A")
    prot = [a for a in atoms if a["resn"] not in WATER]
    waters = [a for a in atoms if a["resn"] in WATER]
    chains = sorted({a["chain"] for a in prot})
    if args.chain not in chains:
        raise SystemExit(f"chain {args.chain} not in {chains}")
    prot = [a for a in prot if a["chain"] == args.chain]

    res = {}
    for a in prot:
        res.setdefault(a["resi"], {})[a["name"]] = a
    nums = sorted(res)
    print(f"1UBQ chain {args.chain}: residues {nums[0]}-{nums[-1]} "
          f"({len(nums)} residues, {len(prot)} atoms, {len(waters)} waters)")
    missing = [n for n in range(nums[0], nums[-1] + 1) if n not in res]
    if missing:
        raise SystemExit(f"gaps in the chain at {missing}: refusing to build")

    # ---- the three mutations ----
    renamed = {}
    report = {"source": args.pdb, "chain": args.chain, "mutations": {}}
    for num, new in MUTATIONS.items():
        if num not in res:
            raise SystemExit(f"residue {num} absent, cannot mutate")
        old = next(iter(res[num].values()))["resn"]
        kept = {k: v for k, v in res[num].items() if k in KEEP_ON_MUTATION}
        if "CB" not in kept:
            raise SystemExit(f"residue {num} ({old}) has no CB to build from")
        dropped = sorted(set(res[num]) - set(kept))
        res[num] = kept
        renamed[num] = new
        report["mutations"][str(num)] = {
            "from": old, "to": new, "kept": sorted(kept), "dropped": dropped}
        print(f"  {old}{num}{new[0] if False else ''} -> {new}: kept {sorted(kept)}, "
              f"tleap will rebuild {len(dropped)} atom(s)")

    # ---- histidine tautomers ----
    for num, r in res.items():
        nm = renamed.get(num, next(iter(r.values()))["resn"])
        if nm != "HIS":
            continue
        others = [a for k, rr in res.items() if k != num for a in rr.values()]
        taut, why = his_tautomer(r, others)
        renamed[num] = taut
        report.setdefault("histidines", {})[str(num)] = {"tautomer": taut, "why": why}
        print(f"  His{num} -> {taut}: {why}")

    # ---- attachment site geometry, the number the preprint quotes ----
    a, b = ATTACHMENT_SITES
    cb = {n: res[n]["CB"]["xyz"] for n in (a, b)}
    ca = {n: res[n]["CA"]["xyz"] for n in (a, b)}
    d_cb = float(np.linalg.norm(cb[a] - cb[b]))
    d_ca = float(np.linalg.norm(ca[a] - ca[b]))
    report["attachment"] = {
        "sites": list(ATTACHMENT_SITES),
        "CB_CB_A": round(d_cb, 2), "CA_CA_A": round(d_ca, 2)}
    print(f"  attachment sites {a}/{b}: CB-CB {d_cb:.2f} A, CA-CA {d_ca:.2f} A "
          f"(the preprint quotes ~2.5 nm)")

    # ---- write ----
    lines, serial = [], 1
    for n in nums:
        rn = renamed.get(n, next(iter(res[n].values()))["resn"])
        # tleap builds hydrogens, and a partial set from the crystal confuses it more
        # than it helps; heavy atoms only, as for the sfGFP preparation.
        for nm, at in res[n].items():
            if at["elem"] == "H":
                continue
            lines.append(fmt(serial, at, rn))
            serial += 1
    lines.append(f"TER   {serial:5d}")
    serial += 1
    if args.keep_waters:
        wres = {}
        for w in waters:
            wres.setdefault((w["chain"], w["resi"]), []).append(w)
        for (chn, ri), ws in wres.items():
            for w in ws:
                if w["elem"] != "O":
                    continue
                w2 = dict(w, chain="W")
                lines.append(fmt(serial, w2, "WAT"))
                serial += 1
        print(f"  kept {len(wres)} crystallographic waters")
    Path(args.out).write_text("\n".join(lines) + "\nEND\n")
    Path(args.out).with_suffix(".json").write_text(json.dumps(report, indent=2))
    print(f"wrote {args.out} ({serial - 1} atoms)")


if __name__ == "__main__":
    main()
