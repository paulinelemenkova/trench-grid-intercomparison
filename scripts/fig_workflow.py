#!/usr/bin/env python3
"""fig_workflow.py -- processing-workflow schematic for the Pacific-trench grid
intercomparison. Renders fig_workflow.pdf (vector) and fig_workflow.png.
Pure matplotlib; no data inputs. Visual language matched to fig_grid_lineage
(Nimbus Sans, Okabe-Ito accents, 1.0 pt outlines, 0.8 pt arrows).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

for _p in ("/root/.fonts/NimbusSans-Regular.otf", "/root/.fonts/NimbusSans-Bold.otf",
           os.path.expanduser("~/.fonts/NimbusSans-Regular.otf"),
           os.path.expanduser("~/.fonts/NimbusSans-Bold.otf")):
    if os.path.exists(_p):
        try: fm.fontManager.addfont(_p)
        except Exception: pass
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans", "Arial", "Liberation Sans",
                                   "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"    # Okabe-Ito, CVD-safe
SFILL = {BLUE: "#DCEAF5", ORANGE: "#FBECD2", GREEN: "#D6F0E6"}
INK, GREY = "#1a1a1a", "#6f6f6f"
PROC_F, PROC_E = "#EDEDED", "#3a3a3a"                    # process (grey) boxes
STEP_F, STEP_E = "#DCE4EC", "#34586B"                    # analysis (blue-grey) boxes
BOX_LW, ARR_LW = 1.0, 0.8                                # requested weights

def box(ax, x, y, w, h, edge, fill, radius=0.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}", linewidth=BOX_LW,
        edgecolor=edge, facecolor=fill, zorder=2))

def txt(ax, x, y, s, fs=9.3, col=INK, weight="normal", style="normal"):
    ax.text(x, y, s, ha="center", va="center", fontsize=fs, color=col,
            weight=weight, style=style, zorder=4, linespacing=1.32)

def arrow(ax, x0, y0, x1, y1, color=GREY, ms=10):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
        mutation_scale=ms, linewidth=ARR_LW, color=color, zorder=1,
        shrinkA=0, shrinkB=0))

fig, ax = plt.subplots(figsize=(8.8, 9.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

# provenance/tooling note (top-right, italic)
txt(ax, 95, 98.5, "Python + GMT / pyGMT", fs=9, col=GREY, style="italic")

# ---- source boxes ----
src = [("GEBCO", BLUE, 6), ("SRTM15+", ORANGE, 38), ("GMRT", GREEN, 70)]
for name, c, x0 in src:
    box(ax, x0, 87, 24, 9, c, SFILL[c])
    txt(ax, x0 + 12, 91.5, name, fs=13, col=c, weight="normal")

# ---- process + analysis boxes ----
box(ax, 16, 71, 68, 11, PROC_E, PROC_F)
txt(ax, 50, 76.5,
    "Harmonise \u2014 clip to [120\u2013290\u00b0E, 60\u00b0S\u201360\u00b0N]\n"
    "and resample onto a common 15\u2033 mesh", fs=9.3)

box(ax, 5, 55, 44, 11, STEP_E, STEP_F)
txt(ax, 27, 60.5, "Pairwise differencing and\ninter-grid dispersion $\\sigma_{\\mathrm{G}}$")

box(ax, 51, 55, 44, 11, STEP_E, STEP_F)
txt(ax, 73, 60.5, "Sample the common\ntrench-normal transects")

box(ax, 5, 40, 44, 11, STEP_E, STEP_F)
txt(ax, 27, 45.5, "Difference and dispersion\nmaps (GMT / pyGMT)")

# widened morphometrics box, text stretched over three lines
box(ax, 51, 40, 44, 11, STEP_E, STEP_F)
txt(ax, 71.5, 44.5,
    "Morphometric metrics: axial depth, width,\n"
    "wall gradient, cross-sectional area, curvature,\n"
    "and grid-induced uncertainty $u_{m}$", fs=9.0)

box(ax, 20, 22, 60, 11, PROC_E, PROC_F)
txt(ax, 50, 27.5,
    "Difference statistics by setting and\nregression on spatial controls")

box(ax, 6, 6, 88, 10, GREEN, SFILL[GREEN])
txt(ax, 50, 11,
    "Difference maps   \u00b7   Metric-spread figures   \u00b7   Uncertainty tables",
    fs=11, col=GREEN, weight="bold")

# ---- flows (all 0.8 pt) ----
# sources -> harmonise (accent colours, converging on the top edge)
arrow(ax, 18, 87, 38, 82.4, BLUE)
arrow(ax, 50, 87, 50, 82.4, ORANGE)
arrow(ax, 82, 87, 62, 82.4, GREEN)
# harmonise -> two branches
arrow(ax, 42, 71, 27, 66.4)
arrow(ax, 60, 71, 73, 66.4)
# branch A -> branch B
arrow(ax, 27, 55, 27, 51.4)
arrow(ax, 73, 55, 73, 51.4)
# branch B -> statistics (converging)
arrow(ax, 27, 40, 40, 33.4)
arrow(ax, 73, 40, 60, 33.4)
# statistics -> outputs
arrow(ax, 50, 22, 50, 16.4)

fig.savefig("fig_workflow.pdf", bbox_inches="tight", pad_inches=0.05)
fig.savefig("fig_workflow.png", dpi=300, bbox_inches="tight", pad_inches=0.05)
print("rendered")
