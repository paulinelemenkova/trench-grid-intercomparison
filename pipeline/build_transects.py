#!/usr/bin/env python3
"""
build_transects.py  --  pipeline stage 2: common trench-normal transect set.
================================================================================
Manuscript Listing 5 (lst:transects). Densify the trench axis at a fixed
along-axis spacing, emit a trench-normal transect of fixed half-length at every
node, and sample each harmonised grid along every transect -- one depth profile
per (grid, transect). The stacked profiles are written to profiles.csv, the
input to the difference statistics (pipeline/diffstats.py) and to the
morphometric metrics.

The trench axis is read from ../data/trench_axes (one lon/lat polyline per
trench); by default the Mariana axis (5_Mariana.txt) is used, matching the
results window. The harmonised grids are those written by prepare_grids.py.
The self-contained figure scripts in ../scripts embed this same construction.

Requires: pygmt (GMT >= 6), numpy, pandas; {GEBCO,SRTM15+,GMRT}_15s.nc from
prepare_grids.py.
"""
import os
import numpy as np
import pandas as pd
import pygmt

AXES = os.path.join(os.path.dirname(__file__), "..", "data", "trench_axes",
                    "5_Mariana.txt")     # lon/lat polyline for the Mariana axis
STEP_KM = 25.0                            # along-axis spacing between transects
HALF_KM = 120.0                           # half-length of each trench-normal line
GRIDS = {"GEBCO": "GEBCO_15s.nc", "SRTM15+": "SRTM15p_15s.nc", "GMRT": "GMRT_15s.nc"}


def normals(axes=AXES, step=STEP_KM, half=HALF_KM):
    """Densify the axis to `step` km, emit a perpendicular of length 2*half at each."""
    dense = pygmt.project(data=axes, unit=True, generate=str(step))
    hl = half / 111.195                                   # km -> degrees
    segs = []
    for i, (lon, lat, az) in enumerate(zip(dense.r, dense.s, dense.az)):
        n = np.radians(az + 90.0)                         # normal direction
        dlon = hl * np.sin(n) / np.cos(np.radians(lat))
        dlat = hl * np.cos(n)
        segs.append((i, lon - dlon, lat - dlat, lon + dlon, lat + dlat))
    return pd.DataFrame(segs, columns=["id", "lon1", "lat1", "lon2", "lat2"])


def sample(grid, seg):
    prof = pygmt.project(center=[seg.lon1, seg.lat1],
                         endpoint=[seg.lon2, seg.lat2], generate="1")   # 1 km step
    return pygmt.grdtrack(points=prof, grid=grid, newcolname="z")


if __name__ == "__main__":
    rows = []
    for name, gpath in GRIDS.items():
        for _, seg in normals().iterrows():
            s = sample(gpath, seg)
            s["grid"], s["transect"] = name, int(seg.id)
            rows.append(s)
    pd.concat(rows, ignore_index=True).to_csv("profiles.csv", index=False)
    print("wrote profiles.csv")
