#!/usr/bin/env python
"""Assemble the ChimeraX single-view renders into paper-style multi-panel figures.

  fig2a_maps.png        Aalpha (top) / Abeta (bottom) deposited sharpened maps, two views each        (Fig 2A)
  fig2b_lowpass6.png    6 A low-pass raw maps, Aalpha | Abeta, B56g C-helix region marked             (Fig 2B)
  fig3c_salt_bridge.png Aalpha | Abeta salt-bridge close-ups                                          (Fig 3C)
  fig3de_interfaces.png A-B56g (top) and A-Calpha (bottom) interfaces, Aalpha | Abeta                 (Fig 3D/E)
  figS5_localres.png    locally filtered maps coloured by local resolution, two views per isoform    (Fig S5A/E)
Panels that have not been rendered yet are skipped.  Run: $DNASPRING_PY scripts/28_montage.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402

REPO = os.environ.get("PP2A_REPO", "/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction")
FIG = os.path.join(REPO, "figures")

LAYOUTS = {
    "fig2a_maps.png": ("Fig 2A analogue: deposited (EMReady-sharpened) maps at the EMDB contour, coloured by subunit "
                       "(A red/grey, B56γ blue, Cα orange)",
                       [["fig2a_maps_Aa.png", "fig2a_maps_Aa_view2.png"], ["fig2a_maps_Ab.png", "fig2a_maps_Ab_view2.png"]],
                       [["Aα-B56γ-Cα (EMD-48300)", "rotated 90°"], ["Aβ-B56γ-Cα (EMD-48224)", "rotated 90°"]]),
    "fig2b_lowpass6.png": ("Fig 2B analogue: half-map averages low-passed to 6 Å at matched σ; blue ribbon = B56γ C-terminal "
                           "helix (9MIP B402–436), carried into the Aβ frame by superposing the B56γ cores",
                           [["fig2b_lowpass6_Aa.png", "fig2b_lowpass6_Ab.png"]], [["Aα", "Aβ"]]),
    "fig3c_salt_bridge.png": ("Fig 3C analogue: HEAT3–HEAT5 salt bridge with deposited map mesh; Cα–Cα distance labelled",
                              [["fig3c_salt_bridge_Aa.png", "fig3c_salt_bridge_Ab.png"]], [["Aα R113–E191", "Aβ R125–E203"]]),
    "fig3de_interfaces.png": ("Fig 3D/E analogue: A-subunit residues at the B56γ (top) and Cα (bottom) interfaces",
                              [["fig3d_ifaceB56_Aa.png", "fig3d_ifaceB56_Ab.png"], ["fig3e_ifaceC_Aa.png", "fig3e_ifaceC_Ab.png"]],
                              [["Aα – B56γ", "Aβ – B56γ"], ["Aα – Cα", "Aβ – Cα"]]),
    "figS5_localres.png": ("Fig S5A/E analogue: locally filtered maps coloured by RELION 3.1 local resolution "
                           "(red 3.0 Å … blue 5.0 Å)",
                           [["figS5_localres_Aa.png", "figS5_localres_Aa_view2.png"], ["figS5_localres_Ab.png", "figS5_localres_Ab_view2.png"]],
                           [["Aα", "Aα rotated 90°"], ["Aβ", "Aβ rotated 90°"]]),
}


def main():
    made = []
    for out, (title, grid, labels) in LAYOUTS.items():
        files = [[os.path.join(FIG, f) for f in row] for row in grid]
        if not all(os.path.exists(f) for row in files for f in row):
            print(f"skip {out}: missing panels")
            continue
        nr, nc = len(grid), len(grid[0])
        fig, axes = plt.subplots(nr, nc, figsize=(5.2 * nc, 4.0 * nr))
        axes = [[axes]] if nr == nc == 1 else ([axes] if nr == 1 else ([[a] for a in axes] if nc == 1 else axes))
        for r in range(nr):
            for c in range(nc):
                ax = axes[r][c]
                ax.imshow(mpimg.imread(files[r][c]))
                ax.set_title(labels[r][c], fontsize=11)
                ax.axis("off")
        fig.suptitle(title, fontsize=11, wrap=True)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, out), dpi=150)
        plt.close(fig)
        made.append(out)
    print("made:", made)


if __name__ == "__main__":
    sys.exit(main())
