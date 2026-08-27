#!/usr/bin/env python3
r"""
fig_studyarea.py  --  Study-area / locator map for the Pacific subduction system.
================================================================================
Figure 1 of "Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global
Bathymetric Grid".

Mollweide (Pacific-centred, oval) shaded-relief basemap overlaid with:
  * the 20 subduction trench axes from the project's own digitised set
    (axes_full/, one lon/lat polyline per trench), drawn with the GMT front
    symbol whose teeth point toward the over-riding plate (per-trench side);
  * optional plate boundaries (Bird 2003, PB2002) and arc volcanoes if present;
  * white graticule, depth colour bar, a compact scale bar, a relief-filled
    locator globe, and a provenance credit.

RELIEF: uses the local GEBCO_2026.nc if found (resampled to a display grid),
otherwise the remote GEBCO relief from the GMT server. NOTE: a full 15-arc-second
global GEBCO cannot be held in memory; if GEBCO_2026.nc is the global tile, either
pre-coarsen it once (gmt grdsample GEBCO_2026.nc -R120/290/-60/60 -I04m -Grelief.nc)
or let the remote-GEBCO fallback run.

The trench-normal transects are NOT drawn here: the analysis transects live only in
the per-window figures; a global "transect set" on the locator would be indicative
only, so it is omitted.

RENDER (skill rule S -- maps render in the sandbox, never locally)
Local use:  pip install pygmt  (needs GMT >= 6);  python fig_studyarea.py
"""

import os
import glob
import numpy as np
import pygmt

# --- region / projection / output -------------------------------------------
REGION = [120, 290, -60, 60]              # Pacific trench belt (0-360 lon)
PROJECTION = "W205/24c"                    # Mollweide, Pacific-centred oval
TITLE = "Pacific subduction system: study area"
STEM = "fig_studyarea"
RELIEF_INC = "04m"                         # display resolution of the basemap relief

# --- relief: local GEBCO_2026 first, then remote GEBCO ----------------------
RELIEF_LOCAL = ["/Volumes/TOSHIBA/DATA/GEBCO_2026.nc", "GEBCO_2026.nc"]

# --- input discovery --------------------------------------------------------
SEARCH_DIRS = ["/sandbox/input", ".", os.path.dirname(os.path.abspath(__file__))]

# --- teeth (subduction polarity) side per trench ----------------------------
# side is relative to the digitisation direction of each axis file; teeth point
# toward the over-riding plate. Flip +l/+r here if any trench reads the wrong way.
SIDE = {"Aleutian": "+l", "Kuril-Kamchatka": "+r", "Japan": "+r", "Izu-Bonin": "+r",
        "Mariana": "+r", "Yap": "+r", "Palau": "+r", "Ryukyu": "+r", "Manila": "+l",
        "Philippine": "+r", "New Britain": "+l", "San Cristobal": "+l", "Vityaz": "+l",
        "New Hebrides": "+l", "Tonga": "+r", "Kermadec": "+r", "Hikurangi": "+r",
        "Puysegur": "+l", "Middle America": "+l", "Peru-Chile": "+l"}

# --- labels for all 20 trenches (lon, lat, justify) -------------------------
LABELS = [(188, 58, "Aleutian", "CM"), (160, 48, "Kuril-Kamchatka", "LM"),
          (146, 37, "Japan", "LM"), (145, 30, "Izu-Bonin", "LM"),
          (150, 17, "Mariana", "LM"), (141, 9, "Yap", "LM"),
          (133, 6.5, "Palau", "RM"), (133, 27, "Ryukyu", "LM"),
          (122, 19, "Manila", "LM"), (129, 10, "Philippine", "LM"),
          (150, -4.5, "New Britain", "CB"), (157, -13, "San Cristobal", "LM"),
          (168, -8, "Vityaz", "CB"), (174, -18, "New Hebrides", "LM"),
          (190, -18, "Tonga", "LM"), (187, -33, "Kermadec", "LM"),
          (174, -41, "Hikurangi", "RM"), (161, -49, "Puysegur", "RM"),
          (258, 22, "Middle America", "CB"), (285, -31, "Peru-Chile", "RM")]


def find(name):
    for d in SEARCH_DIRS:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


def find_axes_dir():
    for d in SEARCH_DIRS:
        for nm in ("axes_full", "axes"):
            cand = os.path.join(d, nm)
            if os.path.isdir(cand):
                return cand
    return None


# --- plate boundaries: accept normalised .txt or raw PB2002 .gmt ------------
def normalize_pb2002(src, dst):
    keep, seg = [], False
    with open(src, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip() or "end of line segment" in line:
                continue
            head = line.lstrip()[:1]
            if head in "+-0123456789." and "," in line:
                lon, _, lat = line.partition(",")
                try:
                    keep.append(f"{float(lon):.5f} {float(lat):.5f}"); seg = True
                except ValueError:
                    pass
            else:
                keep.append("> " + line); seg = True
    if seg:
        with open(dst, "w") as fh:
            fh.write("\n".join(keep) + "\n")
    return dst


def resolve_plate_boundaries():
    p = find("PB2002.txt")
    if p:
        return p
    raw = find("PB2002_boundaries.gmt")
    if raw:
        return normalize_pb2002(raw, "PB2002_norm.txt")
    return None


def basemap_relief():
    """Local GEBCO_2026 resampled to a display grid, else remote GEBCO."""
    for p in RELIEF_LOCAL:
        if os.path.isfile(p):
            try:
                out = "relief_display.nc"
                pygmt.grdsample(grid=p, region=REGION, spacing=RELIEF_INC, outgrid=out)
                print(f"relief: local GEBCO_2026 ({p}) -> {out} at {RELIEF_INC}")
                return out
            except Exception as e:                       # convention / memory issue
                print(f"note: could not use {p} ({e}); using remote GEBCO instead.")
                break
    print("relief: remote GEBCO from the GMT server")
    return pygmt.datasets.load_earth_relief(resolution=RELIEF_INC, region=REGION,
                                            data_source="gebco")


# ============================================================================
plates = resolve_plate_boundaries()
volc = find("volcanoes.gmt")
axes_dir = find_axes_dir()
grid = basemap_relief()

fig = pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain", MAP_GRID_PEN_PRIMARY="0.35p,white@30",
             FONT_TITLE="16p,Helvetica-Bold", FONT_ANNOT_PRIMARY="9p",
             MAP_TITLE_OFFSET="10p")

# 1. shaded-relief bathymetry with a white graticule, coordinates on all sides
pygmt.makecpt(cmap="seafloor", series=[-9000, 0])
fig.grdimage(grid=grid, projection=PROJECTION, region=REGION, cmap=True,
             shading="+a315+nt0.6", frame=["afg", f"WESN+t{TITLE}"])

# 2. land + coastlines over the ocean relief
fig.coast(land="gray70", shorelines="thinnest,gray30", resolution="l", area_thresh=5000)

# 3. plate boundaries (optional; Bird 2003, PB2002)
if plates:
    fig.plot(data=plates, pen="0.9p,gray15")

# 4. subduction trench axes with the front (teeth) symbol, per-trench polarity
if axes_dir:
    files = sorted(glob.glob(os.path.join(axes_dir, "[0-9]*_*.txt")),
                   key=lambda p: int(os.path.basename(p).split("_")[0]))
    for path in files:
        name = os.path.basename(path)[:-4].split("_", 1)[1].replace("_", " ")
        side = SIDE.get(name, "+l")
        fig.plot(data=path, pen="1.3p,firebrick",
                 style=f"f0.5c/0.13c{side}+t", fill="firebrick")
else:
    print("WARNING: no axes_full/ directory found next to the script; "
          "trench axes not drawn.")

# 5. ocean label
fig.text(x=204, y=7, text="P A C I F I C   O C E A N", no_clip=True,
         font="22p,Helvetica,white@25")

# 5. arc volcanoes (optional context)
if volc:
    fig.plot(data=volc, style="t0.12c", fill="orange", pen="0.2p,black")

# 6. trench name labels (faint white pillow)
for lon, lat, name, just in LABELS:
    fig.text(x=lon, y=lat, text=name, font="8p,Helvetica-Bold,black",
             fill="white@25", clearance="0.06c/0.06c", justify=just, no_clip=True)

# 7. depth colour bar
fig.colorbar(cmap=True, position="JBC+w11c/0.35c+h+o0/1.2c",
             frame=["xa2000f1000+lDepth", "y+lm"])

# 8. compact scale bar (2000 km, small font)
with pygmt.config(FONT_ANNOT_PRIMARY="8p", FONT_LABEL="8p"):
    fig.basemap(map_scale="g222/-48+w2000k+c0")

# 9. locator globe inset WITH relief (upper-right corner)
with fig.inset(position="jTR+w4.2c+o0.3c", box="+p0.8p,black+gwhite"):
    fig.grdimage(grid="@earth_relief_30m", region="g", projection="G205/10/4.2c",
                 cmap=True, shading="+a315+nt0.5")
    fig.coast(region="g", projection="G205/10/4.2c", land="gray70",
              shorelines="thinnest,gray40", area_thresh=10000)
    top = np.arange(120, 291, 4); rig = np.arange(-60, 61, 4)
    bot = np.arange(290, 119, -4); lef = np.arange(60, -61, -4)
    bx = np.concatenate([top, np.full(rig.size, 290), bot, np.full(lef.size, 120)])
    by = np.concatenate([np.full(top.size, -60), rig, np.full(bot.size, 60), lef])
    fig.plot(x=bx, y=by, pen="1.2p,red", region="g", projection="G205/10/4.2c")

# 10. credit / provenance (below the colour bar)
credit = ("Relief: GEBCO 2026." + (" Plate boundaries: Bird (2003)." if plates else "")
          + " Trench axes: project data (20 trenches). Projection: Mollweide centred "
            "on 205\u00b0E. Rendered in GMT/PyGMT. Source: authors.")
fig.text(position="BC", justify="TC", offset="0/-2.9c", no_clip=True,
         font="7p,Helvetica,gray30", text=credit)

fig.savefig(f"{STEM}.pdf")
fig.savefig(f"{STEM}.png", dpi=300)
print("wrote", f"{STEM}.pdf", f"{STEM}.png")
