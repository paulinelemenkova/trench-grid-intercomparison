#!/usr/bin/env python3
"""
prepare_grids.py  --  pipeline stage 1: harmonise the three grids.
================================================================================
Manuscript Listing 4 (lst:prep). Clip GEBCO, SRTM15+ and GMRT to the common
region and resample them onto a shared 15 arc-second, pixel-registered mesh,
writing one harmonised NetCDF per grid ({name}_15s.nc). GMRT, delivered finer
than the target mesh, is low-pass filtered (Gaussian, ~2x target width) BEFORE
decimation so its multibeam swaths are not aliased into the difference fields.

The paper shows this stage at Pacific-belt scale ([120, 290, -60, 60]); here it
is parametrised to the Mariana / Challenger Deep analysis window used for the
results, so it runs against the GMRT tile from pipeline/download_gmrt.sh. The
self-contained figure scripts in ../scripts embed this same recipe, so this
stage is optional -- it exists to make the documented pipeline reproducible end
to end and to emit shareable harmonised grids for downstream stages.

Requires: pygmt (GMT >= 6); GMRT_Mariana.nc in the working directory.
"""
import pygmt

REGION = [141, 148, 9, 16]        # Mariana analysis window [W, E, S, N]
SPACING = "15s"                    # common node spacing (15 arc-seconds)

SOURCES = {
    "GEBCO":   "@earth_gebco_15s",     # GEBCO 2026 via the GMT remote server
    "SRTM15+": "@earth_relief_15s",    # IGPP SRTM15+ V2.7 (Tozer et al. 2019)
    "GMRT":    "GMRT_Mariana.nc",      # local tile (pipeline/download_gmrt.sh)
}
FINE = {"GMRT"}                        # grids finer than the target mesh


def prepare(name, path, region=REGION, spacing=SPACING, out=None):
    clipped = pygmt.grdcut(grid=path, region=region)
    if name in FINE:                                        # anti-alias first
        clipped = pygmt.grdfilter(grid=clipped, filter="g0.9", distance="1")
    harmon = pygmt.grdsample(grid=clipped, region=region, spacing=spacing,
                             registration="pixel")
    if out:
        harmon.to_netcdf(out)
    return harmon


if __name__ == "__main__":
    harmonised = {name: prepare(name, p, out=f"{name}_15s.nc")
                  for name, p in SOURCES.items()}
    print("wrote:", ", ".join(f"{n.replace('+', 'p')}_15s.nc" for n in harmonised))
