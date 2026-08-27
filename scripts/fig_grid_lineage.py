#!/usr/bin/env python3
"""fig_grid_lineage.py -- construction and provenance schematic for the three
bathymetric grids (GEBCO, SRTM15+, GMRT). Renders fig_grid_lineage.pdf (vector)
and fig_grid_lineage.png. Pure matplotlib; no data inputs.

Four columns: grid identity/steward | data inputs (dashed source chips) |
construction (dashed method sub-block) | released-product spec. Dashed grey
links and blue-outlined GEBCO chips mark SHARED inputs, so the figure also shows
the three grids are only partly independent (GEBCO & SRTM15+ share altimetric
gravity from Sandwell & Smith; GMRT uses GEBCO for gap-fill).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# --- Arial-metric font: Nimbus Sans (URW) if present, else Liberation/Arial ---
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

BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"   # Okabe-Ito, CVD-safe
FILL = {BLUE: "#DCEAF5", ORANGE: "#FBECD2", GREEN: "#D6F0E6"}
INK, GREY = "#1a1a1a", "#6f6f6f"
DASH = (0, (3, 2))
SOLID_LW, DASH_LW = 1.0, 0.9                     # requested line weights

def box(ax, x, y, w, h, edge, fill, radius=0.5, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        linewidth=(DASH_LW if dashed else SOLID_LW),
        edgecolor=edge, facecolor=fill,
        linestyle=(DASH if dashed else "solid"), zorder=2))

def txt(ax, x, y, s, fs=7.3, col=INK, ha="center", va="center", weight="normal",
        rot=0):
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=col, weight=weight,
            zorder=4, linespacing=1.22, rotation=rot)

def arrow(ax, x0, y0, x1, y1, color=GREY, dashed=False, ms=10, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
        mutation_scale=ms, linewidth=(DASH_LW if dashed else SOLID_LW),
        color=color, linestyle=(DASH if dashed else "solid"),
        zorder=1, shrinkA=0, shrinkB=0))

def chip(ax, x, y, w, h, s, edge, fs=5.9):
    box(ax, x, y, w, h, edge, "#FFFFFF", radius=0.3, dashed=True)
    txt(ax, x + w/2, y + h/2, s, fs=fs, col=INK)

fig, ax = plt.subplots(figsize=(9.2, 5.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

IDX, IDW = 1.5, 12.0
INX, INW = 18.0, 21.0
COX, COW = 42.5, 21.0
OUX, OUW = 67.0, 27.0

for cx, t in [(IDX+IDW/2, "Grid / source"), (INX+INW/2, "Data inputs"),
              (COX+COW/2, "Construction & gridding"), (OUX+OUW/2, "Released product")]:
    txt(ax, cx, 96, t, fs=10, weight="bold")

lanes = [
    ("GEBCO", BLUE, "Seabed 2030", "global \u00b7 IHO\u2013IOC",
     "Compiled soundings +\npredicted depth in gaps",
     [("altimetric\ngravity (gaps)", BLUE), ("ship, multibeam,\nregional grids", BLUE)],
     "Regional grids merged\nonto one global mesh",
     "per-source blending;\nsource-type (TID) grid",
     ["GEBCO_2026", "15\u2033 node (~450 m)", "global coverage",
      "~20\u201325% sounding-based"]),
    ("SRTM15+", ORANGE, "Scripps \u00b7 UCSD", "global",
     "Altimetric gravity +\nsparse ship soundings",
     [("altimetric\ngravity", ORANGE), ("sparse ship\nsoundings", ORANGE)],
     "Gravity-to-depth\nprediction, tied to soundings",
     "biharmonic spline;\ngravity\u2013depth transfer",
     ["SRTM15+ V2", "15\u2033 node (~450 m)", "global coverage",
      "altimetry-led in gaps"]),
    ("GMRT", GREEN, "LDEO", "swath-limited",
     "Multibeam swaths +\nGEBCO to fill gaps",
     [("multibeam\nswath", GREEN), ("GEBCO\n(gap-fill)", BLUE)],
     "Swaths mosaicked at\nnative source resolution",
     "native-resolution\ntiling; no re-gridding",
     ["GMRT v4", "\u2264100 m where surveyed", "swath-limited", "GEBCO elsewhere"]),
]

Y = {"GEBCO": 80, "SRTM15+": 54, "GMRT": 28}
grav = {}

for (name, c, steward, cover, insum, chips, consum, method, spec), key in zip(
        lanes, ["GEBCO", "SRTM15+", "GMRT"]):
    cy = Y[key]; mcy = cy + 5.75
    box(ax, IDX, cy - 10, IDW, 20, c, FILL[c])
    txt(ax, IDX + IDW/2, cy + 4.6, name, fs=8.0, col=c)
    txt(ax, IDX + IDW/2, cy - 1.0, steward, fs=6.0, col=GREY)
    txt(ax, IDX + IDW/2, cy - 5.4, cover, fs=5.9, col=GREY)
    box(ax, INX, cy + 1.5, INW, 8.5, c, "#FFFFFF")
    txt(ax, INX + INW/2, mcy, insum, fs=7.1)
    cw = (INW - 1.4) / 2
    for i, (ch, ce) in enumerate(chips):
        cx0 = INX + i * (cw + 1.4)
        chip(ax, cx0, cy - 8.4, cw, 6.6, ch, ce)
        if i == 0 and key in ("GEBCO", "SRTM15+"):
            grav[key] = cy - 8.4 + 3.3
    box(ax, COX, cy + 1.5, COW, 8.5, c, "#FFFFFF")
    txt(ax, COX + COW/2, mcy, consum, fs=7.1)
    chip(ax, COX + 1.2, cy - 8.4, COW - 2.4, 6.6, method, c)
    box(ax, OUX, cy - 9.5, OUW, 19, c, FILL[c])
    for j, line in enumerate(spec):
        txt(ax, OUX + OUW/2, cy + 6.0 - j * 4.3, line,
            fs=(8.0 if j == 0 else 6.8), col=(c if j == 0 else INK))
    arrow(ax, IDX + IDW, cy, INX, mcy, c)
    arrow(ax, INX + INW, mcy, COX, mcy, c)
    arrow(ax, COX + COW, mcy, OUX, cy, c)

# shared altimetric gravity -- dashed link in the clear gutter, compact label
xl = INX - 3.2
arrow(ax, INX, grav["GEBCO"], xl, grav["GEBCO"], GREY, dashed=True, style="-")
arrow(ax, INX, grav["SRTM15+"], xl, grav["SRTM15+"], GREY, dashed=True, style="-")
arrow(ax, xl, grav["GEBCO"], xl, grav["SRTM15+"], GREY, dashed=True, ms=8,
      style="<|-|>")
txt(ax, xl + 1.5, grav["GEBCO"] - 6.5,
    "shared altimetric\ngravity", fs=5.6, col=GREY, rot=90)

# legend + note, pulled up close to the bottom row
ly = 11.5
arrow(ax, 31, ly, 35, ly, INK)
txt(ax, 36, ly, "primary build flow", fs=6.7, ha="left")
arrow(ax, 57, ly, 61, ly, GREY, dashed=True)
txt(ax, 62, ly, "shared / auxiliary input", fs=6.7, col=GREY, ha="left")
txt(ax, 50, 6.5,
    "Dashed links and blue GEBCO chips mark shared inputs \u2014 the three grids are not fully independent.",
    fs=6.5, col=GREY)

fig.savefig("fig_grid_lineage.pdf", bbox_inches="tight", pad_inches=0.05)
fig.savefig("fig_grid_lineage.png", dpi=300, bbox_inches="tight", pad_inches=0.05)
print("rendered")
