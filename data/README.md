# Data

## Included here

`trench_axes/` — digitised axis polylines (lon lat, one file per Pacific
subduction trench, `N_Name.txt`) used to place the study-area trench traces
(Fig. 1) and the trench-normal transects. `_ranges.txt` lists the per-trench
longitude/latitude extents. These are the authors' own curated data and are
released with the code.

## Not redistributed here (open access — fetch as noted)

These grids are large and openly hosted, so they are not bundled. The scripts
either download them automatically or expect a local tile.

| Dataset | How to obtain | Used by |
|---|---|---|
| **GEBCO 2026** (15″) | auto: GMT remote server `@earth_gebco_15s` | all analysis scripts |
| **SRTM15+ V2.7** (15″) | auto: GMT remote server `@earth_relief_15s` | all analysis scripts |
| **GMRT v4** tile + multibeam mask | run `../pipeline/download_gmrt.sh` → `GMRT_Mariana.nc`, `GMRT_Mariana_mask.nc` | all analysis scripts; mask used by `fig_grid_triptych.py`, `fig_controls.py` |
| **GlobSed v3** (sediment thickness) | download from the GlobSed data repository (Straume et al., 2019) and set its path in `fig_controls.py` / `controls_regression.py` | `fig_controls.py`, `controls_regression.py` |

The GEBCO 2026 relief used as the Fig. 1 basemap can be the GMT-server grid or a
local `GEBCO_2026.nc`; see the header of `../scripts/fig_studyarea.py`.

Place the downloaded `GMRT_Mariana*.nc` next to the scripts you run (or in the
directory the script searches — see each script header).
