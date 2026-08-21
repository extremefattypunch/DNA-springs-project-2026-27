"""Shared plotting style: the validated categorical palette and recessive chrome.

Palette taken unchanged from the data-viz reference instance (adopted July 2026),
used in its documented fixed slot order and within its documented slot limits:
line and bar forms use the adjacent pairlist, where the eight-slot order is
validated (worst adjacent CVD dE 9.1, normal-vision 19.6, OKLab x100); the one
scatter-form figure keeps to the first three slots, which are the ones validated
on the all-pairs list (CVD dE 9.2, normal-vision 24.0).

The validator script is Node-based and this cluster has no node, so the palette is
used exactly as documented rather than re-stepped -- re-stepping is what would have
needed re-validation.  Three light-mode slots (magenta, yellow, aqua) sit below 3:1
on a light surface, so the relief rule applies: every series is directly labelled and
every figure ships the CSV it was drawn from.
"""
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
          "#008300", "#4a3aa7", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8983"
GRID = "#e4e3dd"

# Fixed system -> slot assignment, so a figure that drops a system never repaints
# the survivors.  Colour follows the entity, not its rank.
SYSTEM_COLOR = {
    "S0_wt": SERIES[0], "S1_tet": SERIES[1], "S2_clicked": SERIES[5],
    "S3_spring27": SERIES[2], "S4_spring40": SERIES[3],
    "S5_spring40nick": SERIES[4], "S6_clamp": SERIES[6],
}
SYSTEM_LABEL = {
    "S0_wt": "WT sfGFP", "S1_tet": "2× Tet2-Et", "S2_clicked": "clicked, no DNA",
    "S3_spring27": "27 bp spring", "S4_spring40": "40 bp spring",
    "S5_spring40nick": "40 bp nicked", "S6_clamp": "force clamp",
}


def apply(mpl):
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "savefig.dpi": 200,
        "font.size": 9, "axes.labelsize": 9.5, "axes.titlesize": 10.5,
        "axes.titleweight": "medium", "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.edgecolor": GRID, "axes.linewidth": 0.8,
        "axes.labelcolor": INK_2, "text.color": INK,
        "xtick.color": INK_MUTED, "ytick.color": INK_MUTED,
        "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "grid.color": GRID, "grid.linewidth": 0.7,
        "legend.frameon": False, "legend.fontsize": 8.5,
        "lines.linewidth": 2.0, "lines.markersize": 5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.constrained_layout.use": True,
    })
