#!/usr/bin/env python3
"""
fig_envelope.py  --  Effect of grid choice on a downstream morphometric relationship.
================================================================================
Figure for "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid" (Lemenkova & Piskarev). Corresponds to \\label{fig:envelope}.

A representative published-style trench-morphometric scaling -- cross-sectional area
versus axial depth -- is re-derived independently on GEBCO, SRTM15+ and GMRT over the
same ensemble of trench-normal transects. The three least-squares fits and the
envelope between them show how far the choice of grid alone would move a downstream
result: if the envelope is narrow relative to the scatter of the relationship, the
grid choice does not change the conclusion; if it is wide, it does.

Per transect the metrics are, as in the metric-spread figure:
  * axial depth   = min(z)                                            [-> |z| in km]
  * cross-sec area= integral of (REF - z) over z <= REF (REF = -6000) [-> km^2]

GRID SOURCES (as elsewhere): GEBCO @earth_gebco_15s, SRTM15+ @earth_relief_15s,
GMRT local tile (anti-aliased) or @earth_synbath_15s stand-in.

RENDER (skill rule S -- maps render in the sandbox, never locally)
------------------------------------------------------------------
Local use:  pip install pygmt   (needs GMT >= 6);   python fig_envelope.py
"""

import os
import sys
import io
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
AXIS = {"Mariana": [(141.30, 11.00), (142.60, 11.37), (144.50, 12.20),
                    (146.30, 13.60), (147.40, 15.40)]}

# --- grid sources -----------------------------------------------------------
GRIDS = ["GEBCO", "SRTM15+", "GMRT"]
GRID_SRC = {"GEBCO": "@earth_gebco_15s", "SRTM15+": "@earth_relief_15s", "GMRT": None}
GMRT_FILE = [f"GMRT_{TRENCH}.nc", "GMRT.nc", "GMRTv4.nc", "gmrt.grd"]
GMRT_STANDIN = "@earth_synbath_15s"
GRID_COLOR = {"GEBCO": "#0072B2", "SRTM15+": "#E69F00", "GMRT": "#009E73"}

# --- sampling / metric parameters (match fig_metric_spread) -----------------
SPACING = "15s"
ANTIALIAS_KM = 0.9
AXIS_STEP = 20.0
HALF_LEN = 110.0
PROF_STEP = 2.0
AREA_REF = -6000.0
TRENCH_MIN = -7000.0

STEM = "fig_envelope"
PW, PH = 11.0, 8.5
R_EARTH = 6371.0
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
    print(f"WARNING: no GMRT tile; using {GMRT_STANDIN} stand-in.", file=sys.stderr)
    return GMRT_STANDIN, True


def prep(src, antialias=False):
    cut = pygmt.grdcut(grid=src, region=REGION)
    if antialias:
        cut = pygmt.grdfilter(grid=cut, filter=f"g{ANTIALIAS_KM}", distance="1")
    return pygmt.grdsample(grid=cut, region=REGION, spacing=SPACING,
                           registration="pixel")


def axis_points():
    lons, lats = [], []
    for a, b in zip(AXIS[TRENCH][:-1], AXIS[TRENCH][1:]):
        dlat = np.radians(b[1] - a[1]); dlon = np.radians(b[0] - a[0])
        la = np.radians(0.5 * (a[1] + b[1]))
        seg = R_EARTH * np.hypot(dlat, dlon * np.cos(la))
        n = max(2, int(seg / 5.0))
        lons.append(np.linspace(a[0], b[0], n)); lats.append(np.linspace(a[1], b[1], n))
    lon = np.concatenate(lons); lat = np.concatenate(lats)
    d = np.zeros(len(lon))
    for i in range(1, len(lon)):
        dlon = np.radians(lon[i] - lon[i - 1]); dlat = np.radians(lat[i] - lat[i - 1])
        la = np.radians(0.5 * (lat[i] + lat[i - 1]))
        d[i] = d[i - 1] + R_EARTH * np.hypot(dlat, dlon * np.cos(la))
    out = []
    for i in np.searchsorted(d, np.arange(d[1], d[-1] - 1, AXIS_STEP)):
        if 0 < i < len(lon) - 1:
            dlon = lon[i + 1] - lon[i - 1]; dlat = lat[i + 1] - lat[i - 1]
            az = np.degrees(np.arctan2(dlon * np.cos(np.radians(lat[i])), dlat))
            out.append((lon[i], lat[i], az + 90.0))
    return out


def transect(lon0, lat0, az):
    dlat = np.degrees(HALF_LEN / R_EARTH) * np.cos(np.radians(az))
    dlon = np.degrees(HALF_LEN / R_EARTH) * np.sin(np.radians(az)) / np.cos(np.radians(lat0))
    n = max(5, int(2 * HALF_LEN / PROF_STEP))
    return (np.linspace(lon0 - dlon, lon0 + dlon, n),
            np.linspace(lat0 - dlat, lat0 + dlat, n),
            np.linspace(-HALF_LEN, HALF_LEN, n))


def metrics(s, z):
    z = np.asarray(z, float)
    m = np.isfinite(z)
    if m.sum() < 10:
        return None
    s = np.asarray(s, float)[m]; z = z[m]
    axial = float(np.min(z))
    area = 0.0
    for i in range(len(z) - 1):
        z1, z2 = z[i], z[i + 1]; ds = abs(s[i + 1] - s[i])
        b1, b2 = z1 <= AREA_REF, z2 <= AREA_REF
        if b1 and b2:
            area += 0.5 * ((AREA_REF - z1) + (AREA_REF - z2)) * ds
        elif b1 or b2:
            tc = (AREA_REF - z1) / (z2 - z1)
            area += (0.5 * (AREA_REF - z1) * ds * tc if b1
                     else 0.5 * (AREA_REF - z2) * ds * (1 - tc))
    return axial, area / 1000.0            # |axial| handled later; area -> km^2


def main():
    gmrt_src, standin = resolve_gmrt()
    src = dict(GRID_SRC); src["GMRT"] = gmrt_src
    H = {g: prep(src[g], antialias=(g == "GMRT")) for g in GRIDS}

    D = {g: [] for g in GRIDS}             # axial depth (km, positive)
    A = {g: [] for g in GRIDS}             # cross-sectional area (km^2)
    for lon0, lat0, az in axis_points():
        lon, lat, s = transect(lon0, lat0, az)
        tp = pd.DataFrame({"x": lon, "y": lat, "s": s})   # carry s through grdtrack
        vals = {}
        ok = True
        for g in GRIDS:
            out = pygmt.grdtrack(points=tp, grid=H[g], newcolname="z")
            mm = metrics(out["s"].to_numpy(), out["z"].to_numpy())
            if mm is None or mm[0] > TRENCH_MIN:
                ok = False; break
            vals[g] = mm
        if not ok:
            continue
        for g in GRIDS:
            D[g].append(-vals[g][0] / 1000.0); A[g].append(vals[g][1])

    fits = {}
    for g in GRIDS:
        x = np.array(D[g]); y = np.array(A[g])
        b, a = np.polyfit(x, y, 1)
        fits[g] = (a, b, x, y)
    slopes = np.array([fits[g][1] for g in GRIDS])
    slopevar = 100.0 * (slopes.max() - slopes.min()) / abs(slopes.mean())

    allx = np.concatenate([D[g] for g in GRIDS])
    ally = np.concatenate([A[g] for g in GRIDS])
    xlo, xhi = allx.min(), allx.max()
    xs = np.linspace(xlo, xhi, 60)
    lines = {g: fits[g][0] + fits[g][1] * xs for g in GRIDS}
    env_lo = np.min([lines[g] for g in GRIDS], axis=0)
    env_hi = np.max([lines[g] for g in GRIDS], axis=0)
    # pooled fit + residual scatter, and envelope width, for the caption
    pb, pa = np.polyfit(allx, ally, 1)
    resid = ally - (pa + pb * allx)
    scatter = float(np.std(resid))
    env_width = float(np.mean(env_hi - env_lo))
    print(f"[{TRENCH}] {len(D['GEBCO'])} transects")
    for g in GRIDS:
        print(f"  {g:<8} area = {fits[g][0]:7.1f} + {fits[g][1]:6.2f} * depth")
    print(f"  slope spread across grids = {slopevar:.1f}%")
    print(f"  mean envelope width = {env_width:.2f} km^2 ;  fit scatter (1sigma) "
          f"= {scatter:.2f} km^2  (ratio {env_width/scatter:.2f})")

    ypad = 0.08 * (ally.max() - ally.min())
    region = [xlo - 0.1, xhi + 0.1, ally.min() - ypad, ally.max() + ypad]

    fig = pygmt.Figure()
    with pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica", FONT_LABEL="10p,Helvetica",
                      MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.9p,black",
                      MAP_GRID_PEN_PRIMARY="thinnest,gray88"):
        fig.basemap(region=region, projection=f"X{PW}c/{PH}c",
                    frame=["xafg+lAxial depth (km)",
                           "yafg+lCross-sectional area (km@+2@+)", "WSne"])
        # grid-choice envelope
        fig.plot(x=np.r_[xs, xs[::-1]], y=np.r_[env_hi, env_lo[::-1]],
                 fill="gray75", transparency=55)
        # per-grid scatter (slight transparency + thin grey outline)
        for g in GRIDS:
            fig.plot(x=fits[g][2], y=fits[g][3], style="c0.10c", fill=GRID_COLOR[g],
                     pen="0.25p,gray40", transparency=25)
        # fit lines, GEBCO drawn last so its (near-coincident) line stays visible
        for g in reversed(GRIDS):
            fig.plot(x=xs, y=lines[g], pen=f"0.8p,{GRID_COLOR[g]}")
        # legend
        spec = io.StringIO()
        for g in GRIDS:
            spec.write(f"S 0.3c - 0.7c - 1.6p,{GRID_COLOR[g]} 0.85c {g}\n")
        spec.write("S 0.3c s 0.28c gray75 0.3p,gray 0.85c Grid-choice envelope\n")
        spec.seek(0)
        fig.legend(spec=spec, position="jTL+o0.2c/0.2c+w4.0c",
                   box="+gwhite@10+p0.4p,gray")
        fig.text(position="BR", justify="BR", offset="-0.25c/0.3c", no_clip=False,
                 font="9p,Helvetica,black",
                 text=f"slope spread across grids: {slopevar:.0f}%")

        fig_w = PW + 1.7
        line1 = ("Cross-sectional area vs. axial depth re-fitted per grid; "
                 "band = envelope between the three fits.")
        line2 = (("GMRT = SYNBATH stand-in.  " if standin else "")
                 + "GEBCO 2026, SRTM15+ V2.7 (IGPP), GMRT v4.  "
                 + "Rendered in pyGMT.  Source: authors.")
        fig.shift_origin(yshift="-1.7c")
        foot = dict(region=[0, fig_w, 0, 1], projection=f"X{fig_w}c/1c")
        fig.text(x=fig_w / 2.0, y=0.58, text=line1, justify="TC", no_clip=True,
                 font="6p,Helvetica,gray40", **foot)
        fig.text(x=fig_w / 2.0, y=0.18, text=line2, justify="TC", no_clip=True,
                 font="6p,Helvetica,gray40", **foot)

    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, f"{STEM}.pdf")
    png = os.path.join(OUT_DIR, f"{STEM}.png")
    fig.savefig(pdf); fig.savefig(png, dpi=300)
    print("wrote:", pdf, png)


if __name__ == "__main__":
    main()
