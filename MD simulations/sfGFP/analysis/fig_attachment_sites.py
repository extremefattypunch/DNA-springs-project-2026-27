#!/usr/bin/env python3
"""Where the chromophore's hydrogen bonds are, and why the roGFP positions cannot take
a DNA spring.

The argument this figure has to carry
-------------------------------------
A memo proposed moving our spring anchors to the positions roGFP uses for its redox
disulfide -- Ser147/Gln204 and Asn149/Ser202 -- because those strands carry the
chromophore's own hydrogen bonds and constraining them is known to swing the 400/490 nm
excitation equilibrium.  The reasoning is sound and the target is right; the geometry is
not.  A disulfide is a 2 A covalent bond, so roGFP's pairs sit 4.6-5.9 A apart.  A DNA
spring needs ~31 A between anchors, because each tether -- the Tet2-Et side chain, the
sTCO adduct, the DNA's own terminal arm -- is 10-15 A on its own.

A reader should be able to see that in one glance, so the two panels are drawn to the
same scale with the distances measured on the structure rather than annotated by hand.
"""
from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
PDB = HERE.parent / "build" / "01_protein" / "2b3p.pdb"
OUT = HERE.parent / "figures"

# 2B3P numbering == canonical avGFP numbering (verified against the structure).
CRO = 66
NETWORK = {96: "Arg96", 148: "His148", 203: "Thr203", 205: "Ser205", 222: "Glu222"}
ROGFP = [(147, 204, "roGFP1 disulfide"), (149, 202, "roGFP2 disulfide")]
BUILT = (133, 149, "our spring")


def render():
    import pymol
    pymol.finish_launching(["pymol", "-qck"])
    from pymol import cmd

    cmd.load(str(PDB), "p")
    cmd.remove("solvent or hetatm and not resn CRO")
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)   # blend with the figure surface, no seam
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("cartoon_transparency", 0.62)
    cmd.set("stick_radius", 0.14)
    cmd.set("sphere_scale", 0.42)
    cmd.set("dash_width", 2.6)
    cmd.set("dash_gap", 0.32)
    cmd.set("label_size", 17)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_position", (0, 0, 2.2))
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.2)

    cmd.show("cartoon", "polymer")
    cmd.color("grey80", "polymer")
    cmd.show("sticks", f"resn CRO")
    cmd.color("palegreen", "resn CRO")
    cmd.util.cnc("resn CRO")

    # ---------- panel A: the network the spring is meant to reach ----------
    # Label offsets are per-residue and hand-tuned: PyMOL places every label at the same
    # vector from its atom, which in a packed pocket drops the text on top of the sticks
    # it is naming. These push each one into the nearest clear space.
    OFFSET = {96: (-7.0, -1.5, 3.0), 148: (-1.0, -6.0, 3.0), 203: (5.5, -2.0, 3.0),
              205: (5.0, 3.5, 3.0), 222: (1.0, 7.0, 3.0)}
    for rs, nm in NETWORK.items():
        cmd.show("sticks", f"resi {rs} and not (name C+N+O)")
        cmd.color("skyblue", f"resi {rs} and elem C")
        cmd.set("label_position", OFFSET.get(rs, (0, 0, 2.2)), f"resi {rs}")
        cmd.label(f"resi {rs} and name CB", f'"{nm}"')
    for rs, atom in ((148, "ND1"), (203, "OG1"), (96, "NH2"), (205, "OG")):
        cmd.distance(f"hb{rs}", f"resi {rs} and name {atom}",
                     f"resn CRO and name OH+O2", 3.6, mode=0)
        cmd.color("grey30", f"hb{rs}")
    cmd.distance("hb222", "resi 205 and name OG", "resi 222 and name OE1", 3.6, mode=0)
    cmd.color("grey30", "hb222")
    cmd.hide("labels", "hb*")
    cmd.orient("resn CRO or resi 148+203+205+222+96")
    cmd.zoom("resn CRO or resi 148+203+205+222+96", 2.6)
    cmd.png(str(OUT / "_att_panelA.png"), width=1500, height=1250, dpi=200, ray=1)

    # ---------- panel B: the candidate anchors, to scale ----------
    cmd.delete("hb*")
    cmd.hide("labels")
    cmd.hide("sticks")
    cmd.show("sticks", "resn CRO")
    cmd.color("palegreen", "resn CRO")
    cmd.util.cnc("resn CRO")
    for a, b, lab in ROGFP:
        for rs in (a, b):
            cmd.show("spheres", f"resi {rs} and name CB")
            cmd.color("tv_orange", f"resi {rs} and name CB")
        cmd.distance(f"d{a}_{b}", f"resi {a} and name CB", f"resi {b} and name CB")
        cmd.color("tv_orange", f"d{a}_{b}")
    a, b, _ = BUILT
    for rs in (a, b):
        cmd.show("spheres", f"resi {rs} and name CB")
        cmd.color("marine", f"resi {rs} and name CB")
    cmd.distance(f"d{a}_{b}", f"resi {a} and name CB", f"resi {b} and name CB")
    cmd.color("marine", f"d{a}_{b}")
    # Deliberately unlabelled. Five residue labels and three distance readouts on one
    # barrel collided into illegibility; the panel only has to show that two pairs are
    # nearly touching and one spans the protein, which colour and the caption do.
    cmd.hide("labels")
    cmd.orient("polymer")
    cmd.turn("y", 25)
    cmd.zoom("polymer", 1.5)
    cmd.png(str(OUT / "_att_panelB.png"), width=1500, height=1400, dpi=200, ray=1)
    print("  rendered both panels")


def compose():
    """Two panels, one caption, distances stated in the text not guessed from the render."""
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import vizstyle as VS

    j = json.loads((HERE / "attachment_candidates.json").read_text())
    P = j["pairs"]
    fig = plt.figure(figsize=(12.6, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.04,
                          left=0.01, right=0.99, top=0.90, bottom=0.17)

    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(mpimg.imread(OUT / "_att_panelA.png")); ax.axis("off")
    ax.set_title("the hydrogen bonds a spring would have to disturb",
                 fontsize=11.5, color=VS.INK, pad=8)
    ax.text(0.5, -0.02, "His148 and Thr203 donate to the chromophore phenolate;\n"
            "Ser205–Glu222 completes the proton wire that sets its charge state",
            transform=ax.transAxes, ha="center", va="top", fontsize=9.4,
            color=VS.INK_2, linespacing=1.5)

    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(mpimg.imread(OUT / "_att_panelB.png")); ax.axis("off")
    ax.set_title("why the roGFP positions cannot take a DNA spring",
                 fontsize=11.5, color=VS.INK, pad=8)
    ax.text(0.5, -0.02,
            f"orange, roGFP's disulfide pairs — S147–Q204 "
            f"{P['S147-Q204']['cb_cb_A']:.1f} Å and N149–S202 "
            f"{P['N149-S202']['cb_cb_A']:.1f} Å apart\n"
            f"blue, the pair we build — D133–N149, "
            f"{P['D133-N149']['cb_cb_A']:.1f} Å, and its tethers alone are 10–15 Å each",
            transform=ax.transAxes, ha="center", va="top", fontsize=9.4,
            color=VS.INK_2, linespacing=1.5)

    fig.patch.set_facecolor(VS.SURFACE)
    fig.savefig(OUT / "fig8_attachment_sites.png", dpi=170,
                bbox_inches="tight", facecolor=VS.SURFACE)
    plt.close(fig)
    for f in ("_att_panelA.png", "_att_panelB.png"):
        (OUT / f).unlink(missing_ok=True)
    print(f"  wrote figures/fig8_attachment_sites.png")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    render()
    compose()
