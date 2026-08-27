#!/usr/bin/env python3
"""
diffstats.py  --  pipeline stage 3: difference statistics and regression.
================================================================================
Manuscript Listing 6 (lst:stats). From the per-transect profiles (profiles.csv,
written by build_transects.py) compute, for each grid pair, the mean, standard
deviation and RMS of the depth difference, and the per-node inter-grid
dispersion sigma_G; write diff_stats.csv. Then regress the grid-induced metric
spread u_m on the candidate spatial controls and write the coefficients.

The regression here follows the paper's Listing 6 (statsmodels OLS on
metric_spread.csv). The definitive, dependency-light version actually used for
the Section 4.5 numbers -- standardised coefficients, Pearson r and R^2 on the
real GlobSed layer -- is ../scripts/controls_regression.py (numpy only); prefer
that one to reproduce the reported values.

Requires: numpy, pandas, statsmodels; profiles.csv (from build_transects.py) and
metric_spread.csv (per-transect u_m and covariates).
"""
import numpy as np
import pandas as pd

# --- pairwise difference statistics and inter-grid dispersion ---------------
prof = pd.read_csv("profiles.csv")          # depth along every transect, per grid
wide = prof.pivot_table("z", ["transect", "dist"], "grid")

pairs = [("GEBCO", "SRTM15+"), ("GEBCO", "GMRT"), ("SRTM15+", "GMRT")]
diffstats = pd.DataFrame({
    f"{a}-{b}": {"mean": (wide[a] - wide[b]).mean(),
                 "std":  (wide[a] - wide[b]).std(ddof=0),
                 "rms":  np.sqrt(((wide[a] - wide[b]) ** 2).mean())}
    for a, b in pairs}).T
sigma_G = wide.std(axis=1, ddof=0)          # inter-grid dispersion sigma_G(x)
diffstats.to_csv("diff_stats.csv")
print("wrote diff_stats.csv")

# --- regress grid-induced metric spread u_m on candidate controls -----------
# (see ../scripts/controls_regression.py for the definitive numpy version)
try:
    import statsmodels.api as sm
    M = pd.read_csv("metric_spread.csv")    # per-transect u_m and covariates
    X = sm.add_constant(M[["depth", "wall_gradient", "sounding_density"]])
    model = sm.OLS(M["u_m"], X).fit()
    print(model.summary())
    model.params.to_csv("regression_coefficients.csv")
    print("wrote regression_coefficients.csv")
except FileNotFoundError:
    print("metric_spread.csv not found; skipping regression "
          "(use scripts/controls_regression.py for the Section 4.5 numbers).")
