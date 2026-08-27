#!/usr/bin/env python3
"""
fig_diff_histograms.py  --  Distributions of pairwise grid differences by setting.
================================================================================
Figure for "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to \\label{fig:histograms}.

The maps and profiles show WHERE and along-section the grids diverge; this figure
shows the DISTRIBUTION of the disagreement and how its spread changes with position
on the trench profile. Nodes are partitioned into three morphological settings,
operationalised by water-depth band (a robust, axis-trace-free proxy):

    outer rise / abyssal  |  trench slope / wall  |  trench axis

and in each band the three pairwise differences (GEBCO-SRTM15+, GEBCO-GMRT,
SRTM15+-GMRT) are shown as normalised step histograms. The narrowing/broadening of
the distributions from the abyssal band to the axis is the message.

GRID SOURCES (as elsewhere)
---------------------------
  * GEBCO   -> @earth_gebco_15s   (GEBCO 2026, GMT server)
  * SRTM15+ -> @earth_relief_15s  (IGPP SRTM15+ V2.7, Tozer et al. 2019, GMT server)
  * GMRT    -> local tile GMRT_<trench>.nc (GridServer layer=topo); else the labelled
              @earth_synbath_15s stand-in. GMRT is low-pass filtered to the common
              15" resolution before differencing (same anti-alias as the other figs).

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Local use:  pip install pygmt   (needs GMT >= 6);   python fig_diff_histograms.py
"""

import os
import sys
import numpy as np
import pygmt

# --- representative-trench windows [W, E, S, N] -----------------------------
TRENCHES = {
    "Mariana":        [141, 148,   9, 16],
    "KurilKamchatka": [150, 161,  43, 52],
    "IzuBonin":       [139, 146,  26, 34],
    "Tonga":          [181, 189, -26, -15],
    "Kermadec":       [178, 186, -36, -26],
    "PeruChile":      [281, 289, -34, -20],
}
TRENCH = "Mariana"
REGION = TRENCHES[TRENCH]

# --- grid sources -----------------------------------------------------------
GRIDS = {"GEBCO": "@earth_gebco_15s", "SRTM15+": "@earth_relief_15s", "GMRT": None}
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"
PAIRS = [("GEBCO", "SRTM15+"), ("GEBCO", "GMRT"), ("SRTM15+", "GMRT")]

# --- morphological settings, operationalised by |depth| band (metres) -------
BANDS = [
    ("Outer rise / abyssal", 3000, 6000),
    ("Trench slope / wall",  6000, 8500),
    ("Trench axis",          8500, 12000),
]

# --- output / layout / style ------------------------------------------------
STEM = "fig_diff_histograms"
SPACING = "15s"
ANTIALIAS_KM = 0.9
CAP = 250                               # x-axis half-range for the difference (m)
BW = 10                                 # histogram bin width (m)
PW, PH = 5.5, 5.0                       # panel width/height (cm)
PAIR_COLORS = {("GEBCO", "SRTM15+"): "204/121/167",
               ("GEBCO", "GMRT"):    "0/158/115",
               ("SRTM15+", "GMRT"):  "230/159/0"}

SEARCH_DIRS = ["/sandbox/input", ".", os.path.dirname(os.path.abspath(__file__))]
OUT_DIR = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."


def find(names):
    for n in names:
        for d in SEARCH_DIRS:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def resolve_gmrt():
    g = find(GMRT_FILE)
    if g:
        return g, False
    print(f"WARNING: no GMRT tile found ({', '.join(GMRT_FILE)}); using "
          f"{GMRT_STANDIN} stand-in.", file=sys.stderr)
    return GMRT_STANDIN, True


def harmonised(src, antialias=False):
    cut = pygmt.grdcut(grid=src, region=REGION)
    if antialias:
        cut = pygmt.grdfilter(grid=cut, filter=f"g{ANTIALIAS_KM}", distance="1")
    return pygmt.grdsample(grid=cut, region=REGION, spacing=SPACING,
                           registration="pixel").values


def stairs(edges, freq):
    """Step-outline coordinates from bin edges and per-bin heights."""
    return np.repeat(edges, 2)[1:-1], np.repeat(freq, 2)


def main():
    gmrt_src, standin = resolve_gmrt()
    srcs = dict(GRIDS); srcs["GMRT"] = gmrt_src
    Z = {name: harmonised(src, antialias=(name == "GMRT"))
         for name, src in srcs.items()}

    diffs = {p: (Z[p[0]] - Z[p[1]]).ravel() for p in PAIRS}
    zmean = np.abs((Z["GEBCO"] + Z["SRTM15+"] + Z["GMRT"]) / 3.0).ravel()
    edges = np.arange(-CAP, CAP + BW, BW)

    print(f"[{TRENCH}] node counts per band and pair spread (std, m):")
    fig = pygmt.Figure()
    with pygmt.config(FONT_TITLE="11p,Helvetica,black", FONT_TAG="12p,Helvetica-Bold",
                      FONT_ANNOT_PRIMARY="8p,Helvetica", FONT_LABEL="9p,Helvetica",
                      MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.8p,black",
                      MAP_GRID_PEN_PRIMARY="thinnest,gray85"):
        with fig.subplot(nrows=1, ncols=3, subsize=(f"{PW}c", f"{PH}c"),
                         margins=["0.9c", "0.8c"]):
            for j, (label, lo, hi) in enumerate(BANDS):
                band = (zmean >= lo) & (zmean < hi)
                nfrac = 100.0 * band.sum() / band.size
                hists = {}
                ymax = 1.0
                for p in PAIRS:
                    v = diffs[p][band]
                    v = v[np.isfinite(v)]
                    counts, _ = np.histogram(v, bins=edges)
                    freq = 100.0 * counts / max(counts.sum(), 1)
                    hists[p] = freq
                    ymax = max(ymax, freq.max())
                    print(f"  {label:<22} {p[0]}-{p[1]:<8} n={v.size:>7}  "
                          f"std={v.std():6.1f}")
                region = [-CAP, CAP, 0, ymax * 1.18]
                proj = f"X{PW}c/{PH}c"
                with fig.set_panel(panel=j):
                    fig.basemap(region=region, projection=proj,
                                frame=["xafg+lGrid difference (m)",
                                       "yafg+lFrequency (%)", f"WSne+t{label}"])
                    fig.plot(x=[0, 0], y=[0, region[3]], pen="0.5p,gray,--",
                             region=region, projection=proj)
                    for p in PAIRS:
                        xs, ys = stairs(edges, hists[p])
                        fig.plot(x=xs, y=ys, pen=f"0.9p,{PAIR_COLORS[p]}",
                                 region=region, projection=proj)
                    # tag outside above top-left; band share annotated inside
                    fig.text(position="TL", justify="BL", offset="0c/0.12c",
                             no_clip=True, font="12p,Helvetica-Bold,black",
                             text=f"{'abc'[j]})", region=region, projection=proj)
                    fig.text(position="TR", justify="TR", offset="-0.15c/-0.2c",
                             font="7p,Helvetica,gray30", region=region, projection=proj,
                             text=f"|z| {lo//1000}-{hi//1000} km, {nfrac:.0f}% of nodes")

        fig_w = 3 * PW + 1.8
        credit = ("Normalised step histograms of pairwise differences by depth-band "
                  "setting" + (".  GMRT = SYNBATH stand-in" if standin else "")
                  + ".  Grids: GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4.  "
                  + "Rendered in pyGMT.  Source: authors.")
        # one shared horizontal legend below the panels, then the credit line
        fig.shift_origin(yshift="-1.8c")
        foot = dict(region=[0, fig_w, 0, 1], projection=f"X{fig_w}c/1c")
        legitems = [("GEBCO-SRTM15+", PAIR_COLORS[("GEBCO", "SRTM15+")], 4.30),
                    ("GEBCO-GMRT",    PAIR_COLORS[("GEBCO", "GMRT")],    8.05),
                    ("SRTM15+-GMRT",  PAIR_COLORS[("SRTM15+", "GMRT")],  11.30)]
        for lbl, col, xs in legitems:
            fig.plot(x=[xs, xs + 0.7], y=[0.62, 0.62], pen=f"0.9p,{col}", **foot)
            fig.text(x=xs + 0.85, y=0.62, text=lbl, justify="ML", no_clip=True,
                     font="8p,Helvetica,black", **foot)
        fig.text(x=fig_w / 2.0, y=0.10, text=credit, justify="TC", no_clip=True,
                 font="6p,Helvetica,gray40", **foot)

    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")
    png = os.path.join(OUT_DIR, f"{STEM}.png")
    fig.savefig(pdf); fig.savefig(png, dpi=300)
    print("wrote:", pdf, png)


if __name__ == "__main__":
    main()
