"""Shared matplotlib style (validated categorical palette, thin marks, recessive grid).

Colours are the first three slots of a palette validated for colour-vision
deficiency (all pairs). Aqua is below 3:1 contrast on the light surface, so every
series also has its own marker shape and a direct or legend label.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e6e5e1"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def apply():
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
            "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
            "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
            "axes.labelsize": 10, "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
            "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
            "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
            "lines.linewidth": 1.6, "lines.markersize": 7, "legend.frameon": False,
            "legend.fontsize": 9, "font.size": 10,
        }
    )
