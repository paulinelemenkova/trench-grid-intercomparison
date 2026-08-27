#!/usr/bin/env python3
"""
fig_metric_spread.py  --  Grid-induced spread of the morphometric metric suite.
================================================================================
Figure for "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to \\label{fig:metricspread}.

The maps/profiles/histograms characterise depth disagreement; this figure carries it
into the FIVE morphometric metrics the paper reports and ranks them by how strongly
each depends on the source grid. An ensemble of trench-normal transects is laid along
the trench; on each transect every metric is computed on GEBCO, SRTM15+ and GMRT, and
the grid-induced dispersion u_m = std across the three grids (Eq. of the metric
section) is expressed relative to the ensemble-mean metric. The distribution of that
relative spread over all transects is shown as a box plot, one box per metric.

Metrics per transect profile (s in km, z in m):
  * axial depth   = min(z)
  * trench width  = horizontal extent with z <= axial + 4000 m           [km]
  * wall gradient = max |dz/ds|                                          [m/km]
  * cross-sec area= integral of (axial+4000 - z) over that extent        [m*km]
  * profile curv. = max |d2z/ds2|                                        [m/km^2]

GRID SOURCES (as elsewhere): GEBCO @earth_gebco_15s, SRTM15+ @earth_relief_15s,
GMRT local tile (anti-aliased to 15" before sampling) or @earth_synbath_15s stand-in.

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Local use:  pip install pygmt   (needs GMT >= 6);   python fig_metric_spread.py
"""

import os
import sys
import numpy as np
import pandas as pd
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

# axis polyline the trench-normal transects are hung from (lon, lat)
AXIS = {
    "Mariana": [(141.30, 11.00), (142.60, 11.37), (144.50, 12.20),
                (146.30, 13.60), (147.40, 15.40)],
}

# --- grid sources -----------------------------------------------------------
GRIDS = {"GEBCO": "@earth_gebco_15s", "SRTM15+": "@earth_relief_15s", "GMRT": None}
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"

# --- sampling / metric parameters -------------------------------------------
SPACING = "15s"
ANTIALIAS_KM = 0.9
AXIS_STEP = 20.0                        # spacing of transects along the axis (km)
HALF_LEN = 110.0                        # half-length of each trench-normal transect (km)
PROF_STEP = 2.0                         # along-transect sampling step (km)
AREA_REF = -6000.0                      # fixed reference depth for width/area (m)
TRENCH_MIN = -7000.0                    # keep transect only if all grids reach below this (m)

# --- output / style ---------------------------------------------------------
STEM = "fig_metric_spread"
METRICS = ["axial", "width", "wall", "area", "curv"]
MLABEL = {"axial": "Axial depth", "width": "Trench width", "wall": "Wall gradient",
          "area": "Cross-sectional area", "curv": "Profile curvature"}
PW, PH = 12.0, 8.0                      # panel size (cm)
SET3 = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3"]  # per-box colours

SEARCH_DIRS = ["/sandbox/input", ".", os.path.dirname(os.path.abspath(__file__))]
OUT_DIR = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."
R_EARTH = 6371.0


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


def axis_points():
    """Densify the axis polyline and return (lon, lat, perpendicular-azimuth) at
    ~AXIS_STEP spacing."""
    lons, lats = [], []
    for (a, b) in zip(AXIS[TRENCH][:-1], AXIS[TRENCH][1:]):
        dlat = np.radians(b[1] - a[1]); dlon = np.radians(b[0] - a[0])
        la = np.radians(0.5 * (a[1] + b[1]))
        seg = R_EARTH * np.hypot(dlat, dlon * np.cos(la))
        n = max(2, int(seg / 5.0))
        lons.append(np.linspace(a[0], b[0], n)); lats.append(np.linspace(a[1], b[1], n))
    lon = np.concatenate(lons); lat = np.concatenate(lats)
    # cumulative distance to pick ~AXIS_STEP samples
    d = np.zeros(len(lon))
    for i in range(1, len(lon)):
        dlon = np.radians(lon[i] - lon[i - 1]); dlat = np.radians(lat[i] - lat[i - 1])
        la = np.radians(0.5 * (lat[i] + lat[i - 1]))
        d[i] = d[i - 1] + R_EARTH * np.hypot(dlat, dlon * np.cos(la))
    picks = np.arange(d[1], d[-1] - 1, AXIS_STEP)
    idx = np.searchsorted(d, picks)
    out = []
    for i in idx:
        if i <= 0 or i >= len(lon) - 1:
            continue
        dlon = lon[i + 1] - lon[i - 1]; dlat = lat[i + 1] - lat[i - 1]
        az = np.degrees(np.arctan2(dlon * np.cos(np.radians(lat[i])), dlat))
        out.append((lon[i], lat[i], az + 90.0))     # perpendicular
    return out


def transect(lon0, lat0, az):
    """Straight trench-normal transect through (lon0,lat0) at azimuth az."""
    d = HALF_LEN
    dlat = np.degrees(d / R_EARTH) * np.cos(np.radians(az))
    dlon = np.degrees(d / R_EARTH) * np.sin(np.radians(az)) / np.cos(np.radians(lat0))
    a = (lon0 - dlon, lat0 - dlat); b = (lon0 + dlon, lat0 + dlat)
    n = max(5, int(2 * HALF_LEN / PROF_STEP))
    lon = np.linspace(a[0], b[0], n); lat = np.linspace(a[1], b[1], n)
    s = np.linspace(-HALF_LEN, HALF_LEN, n)
    return lon, lat, s


def metrics(s, z):
    z = np.asarray(z, float)
    m = np.isfinite(z)
    if m.sum() < 10:
        return None
    s = np.asarray(s, float)[m]; z = z[m]
    axial = float(np.min(z))
    slope = np.gradient(z, s)
    wall = float(np.nanmax(np.abs(slope)))
    curv = float(np.nanmax(np.abs(np.gradient(slope, s))))
    # continuous width and cross-sectional area below AREA_REF (interpolated
    # crossings, so the metrics are not integer-quantised)
    width = 0.0; area = 0.0
    for i in range(len(z) - 1):
        z1, z2 = z[i], z[i + 1]
        ds = abs(s[i + 1] - s[i])
        b1, b2 = z1 <= AREA_REF, z2 <= AREA_REF
        if b1 and b2:
            width += ds
            area += 0.5 * ((AREA_REF - z1) + (AREA_REF - z2)) * ds
        elif b1 or b2:
            tc = (AREA_REF - z1) / (z2 - z1)
            if b1:
                width += ds * tc;      area += 0.5 * (AREA_REF - z1) * ds * tc
            else:
                width += ds * (1 - tc); area += 0.5 * (AREA_REF - z2) * ds * (1 - tc)
    return dict(axial=axial, width=width, wall=wall, area=area, curv=curv)


def rel_spread(vals):
    a = np.array(vals, float)
    mean = a.mean()
    return 100.0 * a.std() / abs(mean) if mean != 0 else np.nan


def main():
    gmrt_src, standin = resolve_gmrt()
    srcs = dict(GRIDS); srcs["GMRT"] = gmrt_src
    H = {name: harmonised(src, antialias=(name == "GMRT"))
         for name, src in srcs.items()}

    rels = {m: [] for m in METRICS}
    pts = axis_points()
    for lon0, lat0, az in pts:
        lon, lat, s = transect(lon0, lat0, az)
        tp = pd.DataFrame({"x": lon, "y": lat, "s": s})   # carry s through grdtrack
        mets = {}
        ok = True
        for name in GRIDS:
            out = pygmt.grdtrack(points=tp, grid=H[name], newcolname="z")
            mm = metrics(out["s"].to_numpy(), out["z"].to_numpy())
            if mm is None:
                ok = False; break
            mets[name] = mm
        if not ok:
            continue
        if any(mets[g]["axial"] > TRENCH_MIN for g in GRIDS):
            continue                                    # not a real trench section
        for m in METRICS:
            rels[m].append(rel_spread([mets[g][m] for g in GRIDS]))

    # order metrics by median relative spread (ranking, most sensitive first)
    stats = {}
    for m in METRICS:
        v = np.array(rels[m], float); v = v[np.isfinite(v)]
        stats[m] = np.percentile(v, [5, 25, 50, 75, 95]) if v.size else np.zeros(5)
    order = sorted(METRICS, key=lambda m: stats[m][2], reverse=True)
    print(f"[{TRENCH}] {len(rels[METRICS[0]])} transects; median grid-induced "
          f"spread u_m/|mean| (%):")
    for m in order:
        print(f"  {MLABEL[m]:<22} median={stats[m][2]:5.1f}  IQR="
              f"[{stats[m][1]:.1f},{stats[m][3]:.1f}]")

    ymin = max(0.02, min(stats[m][0] for m in METRICS if stats[m][0] > 0) / 1.5)
    ymax = max(stats[m][4] for m in METRICS) * 1.6
    region = [0.5, 5.5, ymin, ymax]

    fig = pygmt.Figure()
    with pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica", FONT_LABEL="10p,Helvetica",
                      MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.9p,black",
                      MAP_GRID_PEN_PRIMARY="thinnest,gray85"):
        fig.basemap(region=region, projection=f"X{PW}c/{PH}cl",
                    frame=["WSne", "yafg+lGrid-induced spread  @~s@~@-m@-/|mean| (%)"])
        for j, m in enumerate(order):
            q05, q25, q50, q75, q95 = stats[m]
            x = j + 1
            fig.plot(x=[x - 0.28, x + 0.28, x + 0.28, x - 0.28, x - 0.28],
                     y=[q25, q25, q75, q75, q25], pen="0.9p,black",
                     fill=SET3[j % len(SET3)])
            fig.plot(x=[x - 0.28, x + 0.28], y=[q50, q50], pen="1.6p,black")
            fig.plot(x=[x, x], y=[q05, q25], pen="0.8p,black")
            fig.plot(x=[x, x], y=[q75, q95], pen="0.8p,black")
            for cap in (q05, q95):
                fig.plot(x=[x - 0.1, x + 0.1], y=[cap, cap], pen="0.8p,black")
            fig.text(x=x, y=ymin, text=MLABEL[m], justify="RM", angle=35,
                     offset="0c/-0.15c", no_clip=True, font="9p,Helvetica,black")

        fig_w = PW + 1.6
        credit = ("Box = IQR, whiskers 5--95th percentile of the grid-induced spread "
                  "over " + f"{len(rels['axial'])} trench-normal transects"
                  + (".  GMRT = SYNBATH stand-in" if standin else "")
                  + ".  GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4.  "
                  + "Rendered in pyGMT.  Source: authors.")
        fig.shift_origin(yshift="-3.0c")
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
