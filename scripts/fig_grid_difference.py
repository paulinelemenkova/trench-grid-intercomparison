#!/usr/bin/env python3
"""
fig_grid_difference.py  --  Pairwise differences between the three grids.
================================================================================
Companion to fig_grid_triptych.py for "Sensitivity of Deep-Sea Trench Morphometry
to the Choice of Global Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to
\\label{fig:diffmaps}.

Where the matched-scale triptych makes the grids look near-identical, the tens-of-
metres disagreements are invisible against the ~11 km depth range. This figure makes
them visible: three panels show the pairwise differences

    a) GEBCO - SRTM15+     b) GEBCO - GMRT     c) SRTM15+ - GMRT

on ONE symmetric diverging colour scale saturated at +/-DIFF_CAP metres. Blue = the
first grid is deeper than the second; red = the first grid is shallower; white =
agreement. Differences beyond +/-DIFF_CAP (steep flanks, ridge crests, the trench
wall) saturate to the end colours, flagged by the colour-bar triangles.

GRID SOURCES (identical to the triptych)
----------------------------------------
  * GEBCO   -> @earth_gebco_15s   (GEBCO 2026 compilation, GMT server)
  * SRTM15+ -> @earth_relief_15s  (IGPP SRTM15+ V2.7, Tozer et al. 2019, GMT server)
  * GMRT    -> local tile GMRT_<trench>.nc (GridServer layer=topo). If absent, the
              script falls back to @earth_synbath_15s as a labelled stand-in.
All three are cut to the common window and resampled onto a shared pixel-registered
15 arc-second mesh so the subtraction is node-for-node exact.

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Local use:  pip install pygmt   (needs GMT >= 6);   python fig_grid_difference.py
"""

import os
import sys
import numpy as np
import pygmt

# --- representative-trench windows [W, E, S, N] (lon 0-360) ------------------
TRENCHES = {
    "Mariana":        [141, 148,   9, 16],   # Challenger Deep (default)
    "KurilKamchatka": [150, 161,  43, 52],
    "IzuBonin":       [139, 146,  26, 34],
    "Tonga":          [181, 189, -26, -15],
    "Kermadec":       [178, 186, -36, -26],
    "PeruChile":      [281, 289, -34, -20],
}
TRENCH = "Mariana"                      # <- switch trench here
REGION = TRENCHES[TRENCH]               # or set a custom [W, E, S, N] box

# --- grid sources -----------------------------------------------------------
GRIDS = {
    "GEBCO":   "@earth_gebco_15s",      # GEBCO 2026
    "SRTM15+": "@earth_relief_15s",     # IGPP SRTM15+ V2.7 (Tozer et al. 2019)
    "GMRT":    None,                    # resolved below to a local tile or stand-in
}
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"     # only if no GMRT tile is present

# the three pairwise differences drawn, as (A, B) -> panel shows A - B
PAIRS = [("GEBCO", "SRTM15+"), ("GEBCO", "GMRT"), ("SRTM15+", "GMRT")]

# --- output / layout --------------------------------------------------------
STEM = "fig_grid_difference"            # writes STEM.pdf (vector) and STEM.png
PANEL_W = 5.8                           # panel width (cm); 3 panels ~ text width
SPACING = "15s"                         # shared node spacing for the subtraction
ANTIALIAS_KM = 0.9                      # Gaussian low-pass full-width for GMRT (~2x 15s)
DIFF_CAP = 150                          # diverging colour scale saturates at +/- this (m)
CMAP_DIFF = "vik"                       # perceptual diverging map (Crameri), hinged at 0
CREDIT_DY = "-2.9c"                     # origin shift for the credit line (more negative = lower)

# --- input-file discovery (sandbox /sandbox/input first, then cwd) -----------
SEARCH_DIRS = ["/sandbox/input", ".", os.path.dirname(os.path.abspath(__file__))]
OUT_DIR = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."


def find(names):
    """First existing path among `names` across SEARCH_DIRS, else None."""
    for n in names:
        for d in SEARCH_DIRS:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def resolve_gmrt():
    """Local GMRT tile if present, else the labelled stand-in."""
    gmrt = find(GMRT_FILE)
    if gmrt:
        return gmrt, False
    print(f"WARNING: no GMRT tile found ({', '.join(GMRT_FILE)}); using "
          f"{GMRT_STANDIN} as a STAND-IN. Drop in the real GMRT tile before "
          f"submission.", file=sys.stderr)
    return GMRT_STANDIN, True


def harmonised(src, antialias=False):
    """Cut to the window and resample onto the shared pixel-registered 15s mesh so
    every grid has identical node coordinates (exact node-for-node subtraction).
    For a grid finer than the target (GMRT, ~100 m), a Gaussian low-pass to the
    target resolution is applied first so decimation does not alias the multibeam
    swaths into the difference."""
    cut = pygmt.grdcut(grid=src, region=REGION)
    if antialias:
        cut = pygmt.grdfilter(grid=cut, filter=f"g{ANTIALIAS_KM}", distance="2")
    return pygmt.grdsample(grid=cut, region=REGION, spacing=SPACING,
                           registration="pixel")


def stats(diff):
    """mean, RMS, 95th-percentile |diff| and max |diff| of a difference grid (m)."""
    v = np.asarray(diff.values, dtype=float)
    v = v[np.isfinite(v)]
    a = np.abs(v)
    return (float(v.mean()), float(np.sqrt((v ** 2).mean())),
            float(np.percentile(a, 95)), float(a.max()))


def main():
    gmrt_src, standin = resolve_gmrt()
    srcs = dict(GRIDS)
    srcs["GMRT"] = gmrt_src

    # load + harmonise each grid once (GMRT is low-pass filtered before decimation)
    H = {name: harmonised(src, antialias=(name == "GMRT"))
         for name, src in srcs.items()}

    # pairwise differences and their statistics
    diffs, tags = {}, ["a", "b", "c"]
    print(f"[{TRENCH}] region={REGION}  pairwise differences (m):")
    print(f"  {'pair':<18}{'mean':>8}{'rms':>8}{'p95|d|':>9}{'max|d|':>9}")
    for a, b in PAIRS:
        da, db = H[a], H[b]
        d = da.copy(data=da.values - db.values)   # positional subtract (grids aligned)
        diffs[(a, b)] = d
        m, r, p95, mx = stats(d)
        print(f"  {a+' - '+b:<18}{m:8.1f}{r:8.1f}{p95:9.1f}{mx:9.1f}")

    fig = pygmt.Figure()
    with pygmt.config(
        FONT_TITLE="11p,Helvetica,black",            # plain (non-bold) panel titles
        FONT_TAG="12p,Helvetica-Bold,black",         # a)/b)/c) tags stay bold
        FONT_ANNOT_PRIMARY="7p,Helvetica,black",
        FONT_LABEL="8p,Helvetica,black",
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="0.8p,black",
        MAP_TICK_PEN_PRIMARY="thinner,black",
        MAP_GRID_PEN_PRIMARY="thinnest,gray70",      # faint grid, visible on white
        MAP_TITLE_OFFSET="4p",
    ):
        # one symmetric diverging CPT, hinged at 0, shared by all three panels
        pygmt.makecpt(cmap=CMAP_DIFF, series=[-DIFF_CAP, DIFF_CAP], continuous=True)

        proj = "M?"
        with fig.subplot(
            nrows=1, ncols=3,
            figsize=(f"{3 * PANEL_W + 1.0}c", f"{PANEL_W + 1.4}c"),
            frame="WSne",
            margins="0.25c",
            autolabel="a)",
            sharey="l",
        ):
            for j, (a, b) in enumerate(PAIRS):
                with fig.set_panel(panel=j):
                    fig.grdimage(grid=diffs[(a, b)], region=REGION, projection=proj,
                                 cmap=True, frame=["afg", f"+t{a} - {b}"])
                    fig.coast(region=REGION, projection=proj,
                              shorelines="0.4p,gray30", resolution="l")
                    if j == 0:
                        fig.basemap(region=REGION, projection=proj,
                                    map_scale="jBL+w100k+o0.35c/0.35c+u")

        # shared diverging colour bar with end triangles (values exceed +/-DIFF_CAP)
        fig.colorbar(
            cmap=True,
            position="JBC+w11c/0.35c+h+o0c/0.55c+e",
            frame=[f"xa{DIFF_CAP // 2}f{DIFF_CAP // 6}+lGrid difference", "y+lm"],
        )

        credit = ("Pairwise differences (A - B), matched +/-"
                  + f"{DIFF_CAP}"
                  + " m diverging scale.  Blue: A deeper; red: A shallower.  "
                  + "Grids: GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4"
                  + (" [GMRT = SYNBATH stand-in]" if standin else "")
                  + ".  Rendered in pyGMT.  Source: authors.")
        fig_w = 3 * PANEL_W + 1.0
        fig.shift_origin(yshift=CREDIT_DY)
        fig.text(x=fig_w / 2.0, y=1.0, region=[0, fig_w, 0, 1],
                 projection=f"X{fig_w}c/1c", justify="TC", no_clip=True,
                 font="6p,Helvetica,gray40", text=credit)

    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")
    png = os.path.join(OUT_DIR, f"{STEM}.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    print("wrote:", pdf, png)


if __name__ == "__main__":
    main()
