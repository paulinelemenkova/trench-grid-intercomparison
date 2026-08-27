#!/usr/bin/env python3
"""
fig_grid_triptych.py  --  The three grids over one representative trench.
================================================================================
Figure for "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to \\label{fig:triptych}.

Three side-by-side shaded-relief panels -- GEBCO, SRTM15+ and GMRT -- of the SAME
trench window at ONE matched colour scale and identical hillshade, so that where
the grids diverge (the deep, sparsely surveyed axis) is directly comparable by eye.
A single shared depth colour bar and a km scale bar serve all three panels.

GRID SOURCES
------------
Two of the three grids the paper compares are served, as independent global
compilations, by the GMT remote data server and are used here directly:
  * GEBCO   -> @earth_gebco_15s   (GEBCO 2026 compilation)
  * SRTM15+ -> @earth_relief_15s  (IGPP SRTM15+ V2.7, Tozer et al. 2019)
GMRT (Ryan et al. 2009) is NOT on the GMT server. Drop the project's own GMRT tile
for this window in beside the script (see GMRT_FILE candidates) and it is used
automatically. If no GMRT tile is found the panel falls back to @earth_synbath_15s
purely so the figure renders standalone; the panel is then labelled "(stand-in)"
and a warning is printed. Replace it with the real GMRT tile before submission.

An optional GMRT multibeam-only tile (GridServer layer=topo-mask, saved as
GMRT_<trench>_mask.nc) adds a small coverage box in the margin BELOW the panels
(bottom-right, clear of the maps): two flat tones -- blue where GMRT has genuine
multibeam data and grey where GMRT inherits GEBCO fill -- so the reader can see how
much of panel c) is real GMRT vs. gap-fill. If the mask tile is absent the box is
silently skipped and the rest of the figure is unchanged.

REPRESENTATIVE TRENCH
---------------------
Default window is the southern Mariana / Challenger Deep -- the deepest, most
divergence-prone setting. Switch trench with one edit: set TRENCH to any key of
TRENCHES below (Mariana, KurilKamchatka, IzuBonin, Tonga, Kermadec, PeruChile), or
set REGION directly to a custom [W, E, S, N] box.

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Runtime `gmt` (GMT 6.6.0, reaches the GMT data server) or `python`; PyGMT is not
preinstalled, so the render call runs `pip install pygmt --break-system-packages`
first. Any local GMRT tile passed through input_files arrives in /sandbox/input/;
the figure is written to /sandbox/output/. The script searches a few candidate
dirs so it runs unchanged in the sandbox and on your machine.

Local use:  pip install pygmt   (needs GMT >= 6);   python fig_grid_triptych.py
"""

import os
import sys
import numpy as np
import xarray as xr
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

# --- grid sources: name -> remote @name or local file -----------------------
# GEBCO and SRTM15+ are the real grids (GMT server); GMRT is a drop-in file.
GRIDS = {
    "GEBCO":   "@earth_gebco_15s",      # GEBCO 2026
    "SRTM15+": "@earth_relief_15s",     # IGPP SRTM15+ V2.7 (Tozer et al. 2019)
    "GMRT":    None,                    # resolved below to a local tile or stand-in
}
# candidate filenames for the project's own GMRT tile for this window
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"     # only if no GMRT tile is present
# optional GMRT multibeam-only tile (GridServer layer=topo-mask) for the coverage
# inset: NaN where GMRT falls back to GEBCO fill, real relief where multibeam exists.
GMRT_MASK_FILE = [f"GMRT_{TRENCH}_mask.nc", "GMRT_mask.nc", "gmrt_mask.grd"]

# --- output / layout --------------------------------------------------------
STEM = "fig_grid_triptych"              # writes STEM.pdf (vector) and STEM.png
PANEL_W = 5.8                           # panel width (cm); 3 panels ~ text width
GAP = 0.45                              # horizontal gap between panels (cm)
CMAP = "turbo"                          # perceptual relief ramp across the depth range
AZIMUTH = 315                           # hillshade illumination azimuth (all panels)
COVER_W = 2.6                           # coverage-box width (cm), placed below the maps
DPI = 300                               # raster preview resolution

# --- input-file discovery (sandbox /sandbox/input first, then cwd) -----------
SEARCH_DIRS = ["/sandbox/input", ".", os.path.dirname(os.path.abspath(__file__))]
OUT_DIR = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."


def find(names):
    """Return the first existing path among `names` across SEARCH_DIRS, else None."""
    for n in names:
        for d in SEARCH_DIRS:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def resolve_grids():
    """Fill the GMRT slot with a local tile if present, else the labelled stand-in."""
    grids, standin = dict(GRIDS), False
    gmrt = find(GMRT_FILE)
    if gmrt:
        grids["GMRT"] = gmrt
    else:
        grids["GMRT"] = GMRT_STANDIN
        standin = True
        print(f"WARNING: no GMRT tile found ({', '.join(GMRT_FILE)}); "
              f"panel c uses {GMRT_STANDIN} as a STAND-IN. "
              f"Drop in the real GMRT tile before submission.", file=sys.stderr)
    return grids, standin


def matched_series(cut):
    """One matched CPT range for all panels: span the ensemble min/max, nice 500 m,
    kept across 0 so the `geo` master hinges land/sea at sea level."""
    zmin = min(float(g.min()) for g in cut.values())
    zmax = max(float(g.max()) for g in cut.values())
    lo = 500.0 * (zmin // 500.0)                 # floor to 500 m
    hi = max(500.0, 500.0 * -(-zmax // 500.0))   # ceil to 500 m, and span 0
    return [lo, hi]


def main():
    import math
    grids, standin = resolve_grids()

    # cut every grid to the common window ONCE (remote @grids download here)
    cut = {name: pygmt.grdcut(grid=src, region=REGION) for name, src in grids.items()}
    series = matched_series(cut)
    print(f"[{TRENCH}] region={REGION}  matched depth scale {series} m")

    # optional GMRT multibeam-only tile -> small coverage inset below the panels
    mask_path = find(GMRT_MASK_FILE)
    mask_da = pygmt.grdcut(grid=mask_path, region=REGION) if mask_path else None
    if mask_da is None:
        print("note: no GMRT topo-mask tile found; coverage inset skipped. "
              "Download GridServer layer=topo-mask as "
              f"GMRT_{TRENCH}_mask.nc to enable it.", file=sys.stderr)

    # Mercator height/width ratio for the window -> map height and coverage-box height
    W, E, S, N = REGION
    mercY = lambda lat: math.log(math.tan(math.radians(45.0 + lat / 2.0)))
    ratio = (mercY(N) - mercY(S)) / math.radians(E - W)
    step = PANEL_W + GAP
    proj = f"M{PANEL_W}c"
    fig_w = 2 * step + PANEL_W
    tags = ["a)", "b)", "c)"]

    fig = pygmt.Figure()
    with pygmt.config(
        FONT_TITLE="11p,Helvetica,black",           # plain (non-bold) panel titles
        FONT_ANNOT_PRIMARY="7p,Helvetica,black",
        FONT_LABEL="8p,Helvetica,black",
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="0.8p,black",
        MAP_TICK_PEN_PRIMARY="thinner,black",
        MAP_GRID_PEN_PRIMARY="thinnest,white@40",   # thin white graticule
        MAP_TITLE_OFFSET="4p",
    ):
        # one matched CPT, built once, shared by all three panels
        pygmt.makecpt(cmap=CMAP, series=series, continuous=True)

        # --- three fixed-width panels, placed by hand (no subplot: cannot squash) ---
        for j, (name, g) in enumerate(cut.items()):
            if j > 0:
                fig.shift_origin(xshift=f"{step}c")
            yside = "W" if j == 0 else "w"          # y annotations on the first panel only
            title = name + (" (stand-in)" if (name == "GMRT" and standin) else "")
            fig.grdimage(grid=g, region=REGION, projection=proj, cmap=True,
                         shading=f"+a{AZIMUTH}+nt1",
                         frame=["xafg", "yafg", f"{yside}Sne+t{title}"])
            fig.coast(region=REGION, projection=proj,
                      shorelines="0.4p,gray15", resolution="h")
            if j == 0:
                fig.basemap(region=REGION, projection=proj,
                            map_scale="jBL+w100k+o0.4c/0.4c")
            # a)/b)/c) tag INSIDE the panel, top-left, subtle white halo for legibility
            fig.text(region=REGION, projection=proj, position="TL", justify="TL",
                     offset="0.15c/-0.15c", no_clip=True,
                     font="12p,Helvetica-Bold,black", fill="white@30",
                     clearance="2p/2p", text=tags[j])

        # back to the first panel's origin for the shared furniture below
        fig.shift_origin(xshift=f"{-2 * step}c")

        # shared depth colour bar, centred under all three panels
        pygmt.makecpt(cmap=CMAP, series=series, continuous=True)
        fig.colorbar(cmap=True,
                     position=f"x{fig_w / 2.0}c/-1.2c+w11c/0.35c+h+jTC",
                     frame=["xa2000f1000+lDepth", "y+lm"])

        credit = ("Grids: GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4"
                  + (" [GMRT = SYNBATH stand-in]" if standin else "")
                  + ".  Matched colour scale and hillshade.  "
                  + "Rendered in pyGMT.  Source: authors.")

        # coverage box (below-right) + credit (below-left), on one row in the margin
        if mask_da is not None:
            cover = xr.where(mask_da.notnull(), 1.0, 0.0)
            cov_cpt = os.path.join(OUT_DIR, "coverage.cpt")
            pygmt.makecpt(cmap="gray85,dodgerblue3", series=[0, 2, 1],
                          categorical=True, output=cov_cpt)
            cw = COVER_W
            ch = cw * ratio
            rx = fig_w - cw
            fig.shift_origin(xshift=f"{rx}c", yshift=f"{-(2.4 + ch)}c")
            fig.grdimage(grid=cover, region=REGION, projection=f"M{cw}c", cmap=cov_cpt)
            fig.basemap(region=REGION, projection=f"M{cw}c", frame="lrtb")
            fig.text(region=REGION, projection=f"M{cw}c", position="BC", justify="TC",
                     offset="0c/-0.1c", no_clip=True, font="5p,Helvetica,black",
                     text="GMRT multibeam coverage")
            fig.shift_origin(xshift=f"{-rx}c")       # same row, back to the left
        else:
            fig.shift_origin(yshift="-2.6c")

        # credit line, left-justified so it clears the coverage box on the right
        fig.text(x=0.1, y=0.25, region=[0, fig_w, 0, 1], projection=f"X{fig_w}c/1c",
                 justify="LB", no_clip=True, font="6p,Helvetica,gray40", text=credit)

    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")
    png = os.path.join(OUT_DIR, f"{STEM}.png")
    fig.savefig(pdf)                    # vector, for LaTeX/print
    fig.savefig(png, dpi=DPI)           # raster preview
    print("wrote:", pdf, png)


if __name__ == "__main__":
    main()
