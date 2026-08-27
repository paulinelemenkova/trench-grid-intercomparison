#!/usr/bin/env python3
"""
fig_profiles.py  --  Trench-normal and along-axis depth profiles from the 3 grids.
================================================================================
Figure for "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to \\label{fig:profiles}.

Where the maps show WHERE the grids disagree, the profiles show BY HOW MUCH along a
section. Four panels: three trench-normal (cross-axis) sections at different points
along the Mariana arc, and one along-axis profile following the trench. In every
panel the three grids (GEBCO, SRTM15+, GMRT) are overlaid and the between-grid
spread (min-max envelope across the three) is shaded, so the reader sees the grids
track one another over the flanks and separate at the axis and on steep walls.

GRID SOURCES (as elsewhere)
---------------------------
  * GEBCO   -> @earth_gebco_15s   (GEBCO 2026, GMT server)
  * SRTM15+ -> @earth_relief_15s  (IGPP SRTM15+ V2.7, Tozer et al. 2019, GMT server)
  * GMRT    -> local tile GMRT_<trench>.nc (GridServer layer=topo); else the labelled
              @earth_synbath_15s stand-in. GMRT is low-pass filtered to the common
              15" resolution before sampling (same anti-alias step as the difference
              figure) so its finer swaths do not alias the profiles.

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Local use:  pip install pygmt   (needs GMT >= 6);   python fig_profiles.py
"""

import os
import sys
import numpy as np
import pandas as pd
import pygmt

# --- representative-trench windows [W, E, S, N] (lon 0-360) ------------------
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

# --- cross-axis (trench-normal) transects: (name, lon1, lat1, lon2, lat2) ----
# and one along-axis polyline following the trench. Edit per trench.
CROSS = {
    "Mariana": [
        ("Challenger Deep", 142.60, 10.50, 142.60, 12.30),
        ("Central arc",     143.70, 12.90, 145.30, 11.50),
        ("Northern limb",   146.00, 14.50, 148.00, 14.50),
    ],
}
ALONG = {
    "Mariana": [(141.30, 11.00), (142.60, 11.37), (144.50, 12.20),
                (146.30, 13.60), (147.40, 15.40)],
}

# --- grid sources -----------------------------------------------------------
GRIDS = {"GEBCO": "@earth_gebco_15s", "SRTM15+": "@earth_relief_15s", "GMRT": None}
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"

# --- output / layout / style ------------------------------------------------
STEM = "fig_profiles"
SPACING = "15s"                         # common node spacing
ANTIALIAS_KM = 0.9                      # Gaussian low-pass full-width for GMRT
STEP_KM = 1.0                           # profile sampling step
PW, PH = 7.6, 4.7                       # panel width/height (cm)
COLORS = {"GEBCO": "0/119/187", "SRTM15+": "238/119/51", "GMRT": "0/153/136"}
BAND = "gray80"                         # between-grid spread fill

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
                           registration="pixel")


def haversine_cum(lon, lat):
    """Cumulative great-circle distance (km) along a lon/lat path."""
    R = 6371.0
    d = np.zeros(len(lon))
    for i in range(1, len(lon)):
        dlon = np.radians(lon[i] - lon[i - 1])
        dlat = np.radians(lat[i] - lat[i - 1])
        la = np.radians(0.5 * (lat[i] + lat[i - 1]))
        d[i] = d[i - 1] + R * np.hypot(dlat, dlon * np.cos(la))
    return d


def straight_track(a, b, step_km=STEP_KM):
    """Points along a straight lon/lat segment at ~step_km spacing."""
    R = 6371.0
    dlon = np.radians(b[0] - a[0]); dlat = np.radians(b[1] - a[1])
    la = np.radians(0.5 * (a[1] + b[1]))
    length = R * np.hypot(dlat, dlon * np.cos(la))
    n = max(2, int(length / step_km))
    return np.linspace(a[0], b[0], n), np.linspace(a[1], b[1], n)


def polyline_track(points, step_km=STEP_KM):
    """Points along a multi-segment lon/lat polyline."""
    lons, lats = [], []
    for a, b in zip(points[:-1], points[1:]):
        lo, la = straight_track(a, b, step_km)
        if lons:                                   # drop duplicated segment joins
            lo, la = lo[1:], la[1:]
        lons.append(lo); lats.append(la)
    return np.concatenate(lons), np.concatenate(lats)


def sample(grid, lon, lat):
    pts = pd.DataFrame({"x": lon, "y": lat})
    out = pygmt.grdtrack(points=pts, grid=grid, newcolname="z")
    return out["z"].to_numpy()


def draw_panel(fig, tag, title, lon, lat, H):
    d = haversine_cum(lon, lat)
    z = {name: sample(H[name], lon, lat) for name in GRIDS}
    stack = np.vstack([z[n] for n in GRIDS])
    zmin, zmax = np.nanmin(stack, axis=0), np.nanmax(stack, axis=0)
    lo = float(np.nanmin(zmin)); hi = float(np.nanmax(zmax))
    pad = 0.05 * (hi - lo + 1)
    region = [0, float(d[-1]), lo - pad, hi + pad]
    proj = f"X{PW}c/{PH}c"
    fig.basemap(region=region, projection=proj,
                frame=["xafg+lDistance (km)", "yafg+lDepth (m)", f"WSne+t{title}"])
    # between-grid spread envelope
    xx = np.concatenate([d, d[::-1]])
    yy = np.concatenate([zmax, zmin[::-1]])
    fig.plot(x=xx, y=yy, fill=BAND, region=region, projection=proj)
    for name in GRIDS:
        fig.plot(x=d, y=z[name], pen=f"0.9p,{COLORS[name]}", label=name,
                 region=region, projection=proj)
    # panel tag, placed OUTSIDE just above the top-left corner
    fig.text(position="TL", justify="BL", offset="0c/0.12c", no_clip=True,
             font="12p,Helvetica-Bold,black", text=tag,
             region=region, projection=proj)


def main():
    gmrt_src, standin = resolve_gmrt()
    srcs = dict(GRIDS); srcs["GMRT"] = gmrt_src
    H = {name: harmonised(src, antialias=(name == "GMRT"))
         for name, src in srcs.items()}

    cross = CROSS[TRENCH]
    fig = pygmt.Figure()
    with pygmt.config(FONT_TITLE="11p,Helvetica,black", FONT_TAG="12p,Helvetica-Bold",
                      FONT_ANNOT_PRIMARY="8p,Helvetica", FONT_LABEL="9p,Helvetica",
                      MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.8p,black",
                      MAP_GRID_PEN_PRIMARY="thinnest,gray85"):
        with fig.subplot(nrows=2, ncols=2, subsize=(f"{PW}c", f"{PH}c"),
                         margins=["1.1c", "1.2c"]):
            # per-panel legend anchor, kept clear of the data
            LEG = ["jBR", "jBL", "jBL", "jBR"]      # a, b, c, d
            # a-c: cross-axis sections
            for j, (name, x1, y1, x2, y2) in enumerate(cross):
                lon, lat = straight_track((x1, y1), (x2, y2))
                with fig.set_panel(panel=j):
                    draw_panel(fig, f"{'abc'[j]})", f"Cross-axis: {name}", lon, lat, H)
                    fig.legend(position=f"{LEG[j]}+o0.15c",
                               box="+gwhite@20+p0.4p,gray")
            # d: along-axis profile
            lon, lat = polyline_track(ALONG[TRENCH])
            with fig.set_panel(panel=3):
                draw_panel(fig, "d)", "Along-axis (S -> N)", lon, lat, H)
                fig.legend(position=f"{LEG[3]}+o0.15c",
                           box="+gwhite@20+p0.4p,gray")

        # credit
        fig_w = 2 * PW + 2.0
        credit = ("Trench sections from GEBCO, SRTM15+ and GMRT; shaded band = "
                  "between-grid spread (min-max)"
                  + (".  GMRT = SYNBATH stand-in" if standin else "")
                  + ".  Grids: GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4.  "
                  + "Rendered in pyGMT.  Source: authors.")
        fig.shift_origin(yshift="-2.3c")
        fig.text(x=fig_w / 2.0, y=1.0, region=[0, fig_w, 0, 1],
                 projection=f"X{fig_w}c/1c", justify="TC", no_clip=True,
                 font="6p,Helvetica,gray40", text=credit)

    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")
    png = os.path.join(OUT_DIR, f"{STEM}.png")
    fig.savefig(pdf); fig.savefig(png, dpi=300)
    print("wrote:", pdf, png)


if __name__ == "__main__":
    main()
