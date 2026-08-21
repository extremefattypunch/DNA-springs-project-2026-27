#!/usr/bin/env python3
"""Chromophore geometry and H-bond network -- the link from mechanics to fluorescence.

Why these observables
---------------------
GFP fluoresces because the phenolate and imidazolinone halves of the CRO chromophore
are held coplanar and rigid inside the barrel.  Twisting about the methine bridge
opens a non-radiative decay channel, so the two bridge dihedrals

    tau  =  CD2-CG2-CB2-CA2   (phenol ring vs. the bridge)
    phi  =  CG2-CB2-CA2-N2    (bridge vs. the imidazolinone)

are the structural proxies for quantum yield: the wider their distributions, the more
non-radiative decay.  Holding the chromophore planar is the job of a specific H-bond
network -- His148 and Thr203 to the phenolate oxygen, Arg96 and Gln94 to the
imidazolinone carbonyl, Glu222 and Ser205 behind it -- so this module measures both
the geometry and the network that enforces it, and counts the waters that get in when
the network loosens.

Atom names follow the CRO template in xFPchromophores.lib: CA1/CB1/OG1/CG1 come from
Thr65, CA2/CB2/CG2/CD*/CE*/CZ/OH from Tyr66, CA3/C3/O3/N3 from Gly67, and
N2/C2/O2/N1/C1 form the imidazolinone.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

# Hydrogen bonds as (donor heavy atom, acceptor).  The donor's hydrogens are found
# from the topology rather than named here, and occupancy is judged on the
# H...acceptor distance and the D-H...A angle -- not on the heavy-atom distance.
#
# That distinction is not pedantic.  His148 in this system sits at 3.2-4.9 A from the
# phenolate oxygen by heavy-atom distance, which a 3.5 A cutoff calls broken 90% of
# the time, while its HD1 proton is 2.3-2.5 A from that oxygen with a near-linear
# geometry -- a perfectly good, if long and mobile, hydrogen bond.  His148 is the
# solvent-exposed gate at the phenolate end of the cavity and it breathes; measuring
# it by the heavy-atom distance alone would report a broken network that is not broken.
HBOND_MAX_H_A = 2.5      # A, hydrogen to acceptor
HBOND_MIN_ANGLE = 120.0  # deg, donor-H-acceptor

HBONDS = [
    ("His148 ND1-H -> chromophore phenolate", ("HIS", 148, "ND1"), ("CRO", None, "OH")),
    ("Thr203 OG1-H -> chromophore phenolate", ("THR", 203, "OG1"), ("CRO", None, "OH")),
    ("Arg96 NH2-H -> imidazolinone O2", ("ARG", 96, "NH2"), ("CRO", None, "O2")),
    ("Gln94 NE2-H -> imidazolinone O2", ("GLN", 94, "NE2"), ("CRO", None, "O2")),
    ("Ser205 OG-H -> Glu222 OE1", ("SER", 205, "OG"), ("GLU", 222, "OE1")),
    ("chromophore OG1-H -> Glu222 OE2", ("CRO", None, "OG1"), ("GLU", 222, "OE2")),
    ("His148 NE2 <- Arg168 backbone NH", ("ARG", 168, "N"), ("HIS", 148, "NE2")),
]
DIHEDRALS = {
    "tau_phenol_bridge": ("CD2", "CG2", "CB2", "CA2"),
    "phi_bridge_imidazolinone": ("CG2", "CB2", "CA2", "N2"),
}


# mdtraj normalises the Amber histidine tautomer names to HIS, so a name check has
# to treat the whole family as one.
_ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HIS": "HIS"}


def load_residue_map(path):
    """index -> (orig_resnum, orig_resname), written by build_system.py."""
    rows = json.loads(Path(path).read_text())
    return {r["index"]: (r["orig_resnum"], r["orig_resname"]) for r in rows}


def resolve(top, resmap, resname, resnum, atom):
    """Find an atom by (residue name, crystallographic residue number, atom name).

    Uses the explicit residue map rather than offset arithmetic on the prmtop's own
    sequential numbering.  Both the number and the name are checked, so a shifted
    map cannot silently return the wrong residue.
    """
    for idx, (num, name) in resmap.items():
        if resnum is not None and num != resnum:
            continue
        if _ALIASES.get(name, name) != _ALIASES.get(resname, resname):
            continue
        for a in top.residue(idx).atoms:
            if a.name == atom:
                return a.index
        raise KeyError(f"{resname}{resnum} has no atom {atom}")
    raise KeyError(f"{resname}{resnum}.{atom} not found in the residue map")


def circular_stats(deg):
    """Circular mean and sd (degrees) plus the planarity twist.

    A dihedral sitting near +/-180 wraps, which makes an arithmetic mean and sd
    meaningless -- the phenol dihedral in this chromophore reported a 175 deg sd
    that way.  ``twist`` is min(|theta|, 180-|theta|): zero when the two halves are
    coplanar in either the cis or trans sense, which is the quantity that actually
    reports on non-radiative decay, and it is immune to both the wraparound and to
    the two-fold symmetry of the phenol ring.
    """
    r = np.radians(deg)
    c, s = np.cos(r).mean(), np.sin(r).mean()
    R = np.hypot(c, s)
    mean = np.degrees(np.arctan2(s, c))
    sd = np.degrees(np.sqrt(-2.0 * np.log(max(R, 1e-12))))
    a = np.abs(((deg + 180.0) % 360.0) - 180.0)
    twist = np.minimum(a, 180.0 - a)
    return mean, sd, twist


def drop_equilibration(t, skip):
    """Drop the first ``skip`` frames (already strided) as settling time.

    The box is still contracting and the solvent-exposed side chains are still
    finding their positions for the first few ns of production; averaging that into
    an equilibrium observable biases it and inflates its spread.  Never drops so much
    that fewer than 10 frames remain.
    """
    if skip and t.n_frames - skip >= 10:
        return t[skip:]
    return t


def analyse(traj_path, top_path, resmap_path, out_prefix, stride=1, skip=0):
    t = md.load(str(traj_path), top=str(top_path), stride=stride)
    t = drop_equilibration(t, skip)
    top = t.topology
    resmap = load_residue_map(resmap_path)
    print(f"loaded {t.n_frames} frames, {t.n_atoms} atoms")

    cro = next(r for r in top.residues if r.name == "CRO")
    cro_idx = {a.name: a.index for a in cro.atoms}

    rows = {}
    # --- bridge dihedrals ---
    dihedral_stats = {}
    for label, names in DIHEDRALS.items():
        quad = np.array([[cro_idx[n] for n in names]])
        ang = np.degrees(md.compute_dihedrals(t, quad)[:, 0])
        mean, sd, twist = circular_stats(ang)
        rows[label] = ang
        rows[f"twist_{label}"] = twist
        dihedral_stats[label] = {"circular_mean_deg": round(float(mean), 3),
                                 "circular_sd_deg": round(float(sd), 3),
                                 "twist_mean_deg": round(float(twist.mean()), 3),
                                 "twist_sd_deg": round(float(twist.std()), 3),
                                 "twist_p95_deg": round(float(np.percentile(twist, 95)), 3)}
        print(f"  {label:<26} circ mean {mean:+7.2f} deg  circ sd {sd:5.2f}  "
              f"twist {twist.mean():5.2f} +/- {twist.std():4.2f} (p95 {np.percentile(twist, 95):5.2f})")

    # --- H-bond geometry and occupancy ---
    bonded_h = {}
    for bond in top.bonds:
        for x, y in ((bond[0], bond[1]), (bond[1], bond[0])):
            if y.element.symbol == "H":
                bonded_h.setdefault(x.index, []).append(y.index)

    occupancy = {}
    for label, (dr, dn, da), (ar, an, aa) in HBONDS:
        try:
            don = resolve(top, resmap, dr, dn, da)
            acc = resolve(top, resmap, ar, an, aa)
        except (KeyError, StopIteration) as e:
            print(f"  {label:<44} SKIPPED ({e})")
            continue
        hs = bonded_h.get(don, [])
        d_heavy = md.compute_distances(t, np.array([[don, acc]]))[:, 0] * 10.0
        if not hs:
            print(f"  {label:<44} {d_heavy.mean():5.2f} A heavy-atom only "
                  f"(donor carries no hydrogen)")
            rows[f"d_{label}"] = d_heavy
            occupancy[label] = {"heavy_mean_A": round(float(d_heavy.mean()), 3),
                                "note": "no donor hydrogen; heavy-atom distance only"}
            continue
        # pick, per frame, whichever donor hydrogen is closest to the acceptor
        dh = md.compute_distances(t, np.array([[h, acc] for h in hs])) * 10.0
        which = dh.argmin(axis=1)
        d_ha = dh[np.arange(t.n_frames), which]
        ang = np.degrees(md.compute_angles(
            t, np.array([[don, h, acc] for h in hs]))[np.arange(t.n_frames), which])
        occ = float(((d_ha < HBOND_MAX_H_A) & (ang > HBOND_MIN_ANGLE)).mean())
        rows[f"d_{label}"] = d_heavy
        rows[f"dHA_{label}"] = d_ha
        rows[f"ang_{label}"] = ang
        occupancy[label] = {
            "heavy_mean_A": round(float(d_heavy.mean()), 3),
            "heavy_sd_A": round(float(d_heavy.std()), 3),
            "H_acceptor_mean_A": round(float(d_ha.mean()), 3),
            "H_acceptor_sd_A": round(float(d_ha.std()), 3),
            "DHA_angle_mean_deg": round(float(ang.mean()), 2),
            "occupancy": round(occ, 4),
            "criterion": f"H...A < {HBOND_MAX_H_A} A and D-H...A > {HBOND_MIN_ANGLE} deg",
        }
        print(f"  {label:<44} heavy {d_heavy.mean():5.2f}  H...A {d_ha.mean():5.2f} "
              f"+/- {d_ha.std():4.2f} A  angle {ang.mean():5.1f}  occ {100 * occ:5.1f}%")

    # --- waters near the chromophore ---
    # compute_neighbors, not a pair matrix: waters x chromophore atoms is ~300k pairs,
    # and materialising that for every frame needs hundreds of MB and got the
    # aggregation killed by the OOM reaper.  The neighbour search is what this is for.
    water_o = top.select("water and name O")
    if len(water_o):
        near = md.compute_neighbors(t, 0.5, np.array(list(cro_idx.values())),
                                    haystack_indices=water_o, periodic=True)
        n_near = np.array([len(x) for x in near], float)
        rows["waters_within_5A"] = n_near
        print(f"  waters within 5 A of the chromophore: "
              f"{n_near.mean():.2f} +/- {n_near.std():.2f}")

    df = pd.DataFrame(rows)
    df.insert(0, "time_ps", t.time if t.time is not None else np.arange(t.n_frames))
    csv = Path(f"{out_prefix}_chromophore.csv")
    df.to_csv(csv, index=False)
    summary = {
        "n_frames": int(t.n_frames),
        "dihedrals": dihedral_stats,
        "hbonds": occupancy,
        "waters_within_5A_mean": (round(float(np.mean(rows["waters_within_5A"])), 3)
                                  if "waters_within_5A" in rows else None),
    }
    Path(f"{out_prefix}_chromophore.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {csv} and {out_prefix}_chromophore.json")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--resmap", required=True,
                    help="residue_map.json written next to the prmtop")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--skip", type=int, default=0,
                    help="frames to drop from the start, after striding")
    a = ap.parse_args()
    analyse(a.traj, a.top, a.resmap, a.out_prefix, a.stride)
