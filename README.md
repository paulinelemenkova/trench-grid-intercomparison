# Trench morphometry vs. bathymetric-grid choice — reproducibility pipeline

Code and curated data for:

> **Sensitivity of Deep-Sea Trench Morphometry to the Choice of Global Bathymetric Grid**
> Polina Lemenkova & Alexey L. Piskarev (under review, 2026).

The study quantifies how far the three open global bathymetric grids — **GEBCO
2026**, **SRTM15+ V2.7** and **GMRT v4** — disagree across Pacific subduction
trenches, and propagates that disagreement into the standard morphometric
metrics. Every figure and statistic in the paper is produced by the scripts
here from open inputs.

## Layout

```
scripts/     one self-contained generator per manuscript figure, + the regression
pipeline/    the modular pipeline stages shown in the Appendix listings
tools/        PNG size utility for the repository/preview figures
data/         curated trench axes (+ notes on the open grids not redistributed)
```

## Figure → script

| Fig. | Script | Content |
|------|--------|---------|
| 1  | `scripts/fig_studyarea.py`        | Pacific study area, barbed trench axes (needs `data/trench_axes/`) |
| 2  | `scripts/fig_grid_lineage.py`     | grid construction/provenance schematic (matplotlib) |
| 3  | `scripts/fig_grid_triptych.py`    | GEBCO / SRTM15+ / GMRT over the Mariana window + coverage inset |
| 4  | `scripts/fig_workflow.py`         | processing-workflow schematic (matplotlib) |
| 5  | `scripts/fig_grid_difference.py`  | pairwise difference maps (±150 m, vik) |
| 6  | `scripts/fig_profiles.py`         | cross- and along-axis depth profiles |
| 7  | `scripts/fig_diff_histograms.py`  | difference distributions by morphological setting |
| 8  | `scripts/fig_metric_spread.py`    | grid-induced spread of each metric |
| 9  | `scripts/fig_controls.py`         | dispersion vs. coverage / sediment / depth |
| 10 | `scripts/fig_envelope.py`         | area–depth scaling re-fitted per grid |
| §4.5 | `scripts/controls_regression.py` | standardised regression of σ_G on the controls |

## Manuscript listing → code

| Listing | What it shows | Where it lives |
|---|---|---|
| 1 (`lst:gmrtdl`) | GMRT tile download | `pipeline/download_gmrt.sh` |
| 2 (`lst:diff`) | harmonisation + inter-grid dispersion | embedded in every `scripts/fig_*.py`; standalone in `pipeline/prepare_grids.py` |
| 3 (`lst:morphometry`) | metric extraction + propagation | embedded in `fig_metric_spread.py`, `fig_envelope.py` |
| 4 (`lst:prep`) | harmonise → `{name}_15s.nc` | `pipeline/prepare_grids.py` |
| 5 (`lst:transects`) | transects → `profiles.csv` | `pipeline/build_transects.py` |
| 6 (`lst:stats`) | difference stats + regression | `pipeline/diffstats.py`; definitive regression in `scripts/controls_regression.py` |
| 7 (`lst:plots`) | difference / spread rendering | the real figures: `fig_grid_difference.py`, `fig_metric_spread.py` |

The `scripts/fig_*.py` files are the **turnkey** path — each is self-contained
(harmonise → sample → measure → render) and reproduces one figure from open
inputs. The `pipeline/` scripts are the same stages factored out to match the
paper's appendix and to emit shareable intermediates (`{name}_15s.nc`,
`profiles.csv`, `diff_stats.csv`).

## Getting the data

GEBCO 2026 and SRTM15+ are pulled automatically from the GMT remote server
(`@earth_gebco_15s`, `@earth_relief_15s`). GMRT is not, so run:

```bash
cd pipeline && ./download_gmrt.sh      # writes GMRT_Mariana.nc + _mask.nc
```

GlobSed v3 (sediment thickness, for Fig. 9 and the regression) must be
downloaded separately and its path set in `fig_controls.py` /
`controls_regression.py`. See `data/README.md`.

## Environment & running

```bash
# GMT >= 6.5 must be installed and on PATH (https://www.generic-mapping-tools.org)
python -m pip install -r requirements.txt
python scripts/fig_grid_triptych.py        # each script writes STEM.pdf + STEM.png
```

Each figure script writes a **vector PDF** (`fig_*.pdf`, the submission format)
and a 300-dpi **PNG** preview. `tools/compress_figures.py` shrinks oversized
preview PNGs below a size cap (default 1.5 MB) without banding — it keeps full
colour and only reduces resolution, dropping to a 256-colour palette solely as a
last resort. Vector PDFs are resolution-independent and are what the journal
receives.

## Citation

See `CITATION.cff`. The archived release DOI (Zenodo) will be added on
acceptance.

## License

Code is released under the MIT License (`LICENSE`). The curated trench-axis data
in `data/trench_axes/` may be reused with attribution to the authors.
