#!/usr/bin/env python3
"""
fig_controls.py  --  What controls where the grids disagree.
================================================================================
Figure for "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to \\label{fig:controls}.

The inter-grid dispersion sigma_G (Eq. of the difference model; std of the three grid
depths at a node) is regressed against three candidate controls to see which governs
where the grids diverge:
  a) sounding-track density  -- proxied by local GMRT multibeam coverage (%)
  b) sediment thickness      -- from a local GlobSed tile if present, else seafloor
                                age @earth_age is used as a clearly-labelled proxy
  c) water depth             -- ensemble-mean depth (km)
Each panel shows the binned-median trend of sigma_G against the control with its
inter-quartile band, plus a faint subsampled scatter.

GRID SOURCES: GEBCO @earth_gebco_15s, SRTM15+ @earth_relief_15s, GMRT local tile
(anti-aliased) or @earth_synbath_15s stand-in; GMRT mask GMRT_<trench>_mask.nc
(GridServer layer=topo-mask) for coverage; optional GlobSed_<trench>.nc for sediment.

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Local use:  pip install pygmt   (needs GMT >= 6);   python fig_controls.py
"""

import os
import sys
import io
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

# --- grid + ancillary sources -----------------------------------------------
GRIDS = {"GEBCO": "@earth_gebco_15s", "SRTM15+": "@earth_relief_15s", "GMRT": None}
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"
MASK_FILE = [f"GMRT_{TRENCH}_mask.nc", "GMRT_mask.nc", "gmrt_mask.grd"]
SED_FILE = ["/Volumes/TOSHIBA/DATA/GlobSed-v3.nc", f"GlobSed_{TRENCH}.nc",
            "GlobSed.nc", "sediment.nc"]
SED_PROXY = "@earth_age_02m"            # labelled proxy if no GlobSed tile present

# --- analysis / style -------------------------------------------------------
STEM = "fig_controls"
ANALYSIS_INC = "3m"                     # working node spacing for the regression (~5 km)
ANTIALIAS_KM = 0.9
COV_FILT_KM = 40.0                      # window for local multibeam-coverage density
NBINS = 14
# fixed x-axis ranges for panels where a data-driven 1--99 pct clips the story
# (depth would otherwise cut the >9 km trench axis where sigma_G peaks)
XRANGE = {"Multibeam coverage (%)": (0.0, 100.0), "Water depth (km)": (2.0, 11.0)}
PW, PH = 5.6, 5.6                       # panel size (cm)
SCATTER_N = 3000                        # faint scatter points per panel
BAND = "200/210/225"

SEARCH_DIRS = ["/sandbox/input", ".", os.path.dirname(os.path.abspath(__file__))]
OUT_DIR = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."


def find(names):
    for n in names:
        for d in SEARCH_DIRS:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


def prep(src, antialias=False):
    cut = pygmt.grdcut(grid=src, region=REGION)
    if antialias:
        cut = pygmt.grdfilter(grid=cut, filter=f"g{ANTIALIAS_KM}", distance="1")
    return pygmt.grdsample(grid=cut, region=REGION, spacing=ANALYSIS_INC,
                           registration="pixel")


def binned(ctrl, sig, nbins=NBINS, lo=None, hi=None):
    good = np.isfinite(ctrl) & np.isfinite(sig)
    c, s = ctrl[good], sig[good]
    if lo is None:
        lo = np.percentile(c, 1)
    if hi is None:
        hi = np.percentile(c, 99)
    edges = np.linspace(lo, hi, nbins + 1)
    cen, med, q1, q3 = [], [], [], []
    for i in range(nbins):
        m = (c >= edges[i]) & (c < edges[i + 1])
        if m.sum() > 30:
            cen.append(0.5 * (edges[i] + edges[i + 1]))
            med.append(np.median(s[m])); q1.append(np.percentile(s[m], 25))
            q3.append(np.percentile(s[m], 75))
    return map(np.array, (cen, med, q1, q3))


def main():
    gmrt = find(GMRT_FILE)
    standin = gmrt is None
    gmrt_src = gmrt if gmrt else GMRT_STANDIN
    if standin:
        print(f"WARNING: no GMRT tile; using {GMRT_STANDIN} stand-in.", file=sys.stderr)

    Z = {"GEBCO": prep(GRIDS["GEBCO"]), "SRTM15+": prep(GRIDS["SRTM15+"]),
         "GMRT": prep(gmrt_src, antialias=True)}
    g, s, r = (Z["GEBCO"].values, Z["SRTM15+"].values, Z["GMRT"].values)
    mean = (g + s + r) / 3.0
    sigma = np.sqrt(((g - mean) ** 2 + (s - mean) ** 2 + (r - mean) ** 2) / 3.0)
    depth = np.abs(mean) / 1000.0                          # km

    # sounding-density proxy: local GMRT multibeam coverage (%)
    mpath = find(MASK_FILE)
    if mpath:
        mask = pygmt.grdsample(grid=pygmt.grdcut(grid=mpath, region=REGION),
                               region=REGION, spacing=ANALYSIS_INC, registration="pixel")
        cov = mask.copy()                      # copy keeps the grid geometry
        cov.values = np.where(np.isfinite(mask.values), 1.0, 0.0)
        # Cartesian filter (width in degrees) sidesteps the geographic-distance
        # check, which some local GMRT-mask netcdfs trip because they don't
        # advertise lon/lat units; at this latitude the difference is <3%.
        filt_deg = COV_FILT_KM / 111.0
        dens = pygmt.grdfilter(grid=cov, filter=f"g{filt_deg:.3f}",
                               distance="0").values * 100
    else:
        print("note: no GMRT mask tile; coverage-density panel skipped.", file=sys.stderr)
        dens = None

    # sediment thickness (GlobSed) or labelled proxy (seafloor age)
    spath = find(SED_FILE)
    if spath:
        sed = prep(spath).values; sed_label = "Sediment thickness (m)"
    else:
        print(f"note: no GlobSed tile; using {SED_PROXY} as a labelled proxy.",
              file=sys.stderr)
        sed = prep(SED_PROXY).values; sed_label = "Seafloor age (Myr, sediment proxy)"

    panels = []
    if dens is not None:
        panels.append(("Multibeam coverage (%)", dens))
    panels.append((sed_label, sed))
    panels.append(("Water depth (km)", depth))

    sig_flat = sigma.ravel()
    rng = np.random.default_rng(0)
    # precompute the binned trend per panel; the shared y-axis is scaled to the
    # trend medians (not the heavy dispersion tail) so the trends stay legible
    panel_stats = []
    maxmed = 1.0
    for label, ctrl in panels:
        cf = ctrl.ravel()
        good = np.isfinite(cf) & np.isfinite(sig_flat)
        cg, sg = cf[good], sig_flat[good]
        fixed = XRANGE.get(label)
        xlo, xhi = fixed if fixed else np.percentile(cg, [1, 99])
        cen, med, q1, q3 = binned(cg, sg, lo=xlo, hi=xhi)
        if med.size:
            maxmed = max(maxmed, float(np.nanmax(med)))
        panel_stats.append((label, cg, sg, xlo, xhi, cen, med, q1, q3))
    ymax = maxmed * 2.6

    fig = pygmt.Figure()
    tags = "abc"
    with pygmt.config(FONT_ANNOT_PRIMARY="8p,Helvetica", FONT_LABEL="9p,Helvetica",
                      FONT_TITLE="10p,Helvetica", MAP_FRAME_TYPE="plain",
                      MAP_FRAME_PEN="0.8p,black", MAP_GRID_PEN_PRIMARY="thinnest,gray88"):
        with fig.subplot(nrows=1, ncols=len(panels), subsize=(f"{PW}c", f"{PH}c"),
                         margins=["0.35c", "0.8c"]):
            for j, (label, cg, sg, xlo, xhi, cen, med, q1, q3) in enumerate(panel_stats):
                region = [xlo, xhi, 0, ymax]
                proj = f"X{PW}c/{PH}c"
                if j == 0:
                    frame = [f"xafg+l{label}",
                             "yafg+lInter-grid dispersion @~s@~@-G@- (m)", "WSne"]
                else:
                    frame = [f"xafg+l{label}", "yafg", "wSne"]  # shared y: no annot
                with fig.set_panel(panel=j):
                    fig.basemap(region=region, projection=proj, frame=frame)
                    k = rng.choice(cg.size, size=min(SCATTER_N, cg.size), replace=False)
                    fig.plot(x=cg[k], y=sg[k], style="c0.03c", fill="gray70",
                             transparency=70, region=region, projection=proj)
                    if cen.size:
                        fig.plot(x=np.r_[cen, cen[::-1]], y=np.r_[q3, q1[::-1]],
                                 fill=BAND, transparency=35, region=region, projection=proj)
                        fig.plot(x=cen, y=med, pen="0.8p,firebrick", region=region,
                                 projection=proj)
                    fig.text(position="TL", justify="BL", offset="0c/0.12c", no_clip=True,
                             font="11p,Helvetica,black", text=f"{tags[j]})",
                             region=region, projection=proj)

        fig_w = len(panels) * PW + 1.6
        credit = ("Inter-grid dispersion vs. candidate controls."
                  + (".  GMRT = SYNBATH stand-in" if standin else "")
                  + ("  Sediment panel uses seafloor age as a proxy." if not spath else "")
                  + "  Grids: GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4.  "
                  + "Rendered in pyGMT.  Source: authors.")
        # one-row legend below the x-labels, then the credit under it
        fig.shift_origin(yshift="-2.0c")
        foot = dict(region=[0, fig_w, 0, 1], projection=f"X{fig_w}c/1c")
        yr = 0.60
        fig.plot(x=[3.70, 4.30], y=[yr, yr], pen="0.8p,firebrick", **foot)
        fig.text(x=4.45, y=yr, text="Binned-median trend", justify="ML", no_clip=True,
                 font="8p,Helvetica,black", **foot)
        fig.plot(x=[8.55], y=[yr], style="c0.13c", fill="gray70", **foot)
        fig.text(x=8.75, y=yr, text="Grid nodes", justify="ML", no_clip=True,
                 font="8p,Helvetica,black", **foot)
        fig.plot(x=[11.55], y=[yr], style="s0.28c", fill=BAND, pen="0.3p,gray", **foot)
        fig.text(x=11.78, y=yr, text="Interquartile range", justify="ML", no_clip=True,
                 font="8p,Helvetica,black", **foot)
        fig.text(x=fig_w / 2.0, y=0.12, text=credit, justify="TC", no_clip=True,
                 font="6p,Helvetica,gray40", **foot)

    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")
    png = os.path.join(OUT_DIR, f"{STEM}.png")
    fig.savefig(pdf); fig.savefig(png, dpi=300)
    print("wrote:", pdf, png)


if __name__ == "__main__":
    main()
