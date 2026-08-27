#!/usr/bin/env bash
# ============================================================================
# download_gmrt.sh -- fetch the GMRT v4 tile for the representative window.
# Manuscript Listing 1 (lst:gmrtdl).
#
# GEBCO 2026 and SRTM15+ V2.7 are pulled automatically from the GMT remote data
# server by the figure scripts (@earth_gebco_15s / @earth_relief_15s), so only
# GMRT -- which is not served that way -- has to be downloaded here. The topo
# layer is the complete grid; the topo-mask layer marks genuine multibeam cover
# (used by the coverage inset in fig_grid_triptych.py and by fig_controls.py).
# ============================================================================
set -euo pipefail

W=141; E=148; S=9; N=16                       # Mariana / Challenger Deep window
BASE="https://www.gmrt.org/services/GridServer"

curl -o GMRT_Mariana.nc \
  "${BASE}?west=${W}&east=${E}&south=${S}&north=${N}&layer=topo&format=netcdf&resolution=high"

curl -o GMRT_Mariana_mask.nc \
  "${BASE}?west=${W}&east=${E}&south=${S}&north=${N}&layer=topo-mask&format=netcdf&resolution=high"

echo "wrote GMRT_Mariana.nc and GMRT_Mariana_mask.nc"
