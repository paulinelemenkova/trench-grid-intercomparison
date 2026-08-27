#!/usr/bin/env python3
"""
controls_regression.py -- 3-covariate regression of the inter-grid dispersion.

Recomputes sigma_G and the three controls exactly as fig_controls.py does (so the
numbers match the figure), then prints the Pearson r of each control, the
standardized regression coefficients, the combined R^2 and the residual, using the
REAL GlobSed sediment layer. Run it in the same folder as fig_controls.py:

    python controls_regression.py

"""
import numpy as np
import pygmt

REGION = [141, 148, 9, 16]
INC = "3m"                                    # matches ANALYSIS_INC in fig_controls.py
GMRT = "GMRT_Mariana.nc"
MASK = "GMRT_Mariana_mask.nc"
SED = "/Volumes/TOSHIBA/DATA/GlobSed-v3.nc"


def prep(src, antialias=False):
    g = pygmt.grdcut(grid=src, region=REGION)
    if antialias:
        g = pygmt.grdfilter(grid=g, filter="g0.9", distance="1")
    return pygmt.grdsample(grid=g, region=REGION, spacing=INC, registration="pixel")


# --- fields (identical recipe to the figure) --------------------------------
G = prep("@earth_gebco_15s").values
S = prep("@earth_relief_15s").values
R = prep(GMRT, antialias=True).values
mean = (G + S + R) / 3.0
sigma = np.sqrt(((G - mean) ** 2 + (S - mean) ** 2 + (R - mean) ** 2) / 3.0)
depth = np.abs(mean) / 1000.0                                     # km

mask = prep(MASK)
cov = mask.copy()
cov.values = np.where(np.isfinite(mask.values), 1.0, 0.0)
dens = pygmt.grdfilter(grid=cov, filter="g0.360", distance="0").values * 100  # % coverage
sed = prep(SED).values                                            # sediment thickness (m)

# --- assemble and clean -----------------------------------------------------
y = sigma.ravel()
X = np.column_stack([dens.ravel(), sed.ravel(), depth.ravel()])
names = ["coverage", "sediment", "depth"]
good = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
y, X = y[good], X[good]
print(f"n = {y.size} nodes\n")

# --- univariate Pearson r ---------------------------------------------------
for n, col in zip(names, X.T):
    print(f"  r({n:9s}) = {np.corrcoef(col, y)[0, 1]:+.3f}")

# --- standardized multiple regression --------------------------------------
Xs = (X - X.mean(0)) / X.std(0)
ys = (y - y.mean()) / y.std()
beta, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
R2 = 1.0 - np.sum((ys - Xs @ beta) ** 2) / np.sum(ys ** 2)
print()
for n, b in zip(names, beta):
    print(f"  std_beta({n:9s}) = {b:+.3f}")
print(f"\n  R^2 (3 covariates) = {R2:.3f}  ({100 * R2:.1f}%)   "
      f"residual = {100 * (1 - R2):.1f}%")
print(f"  strongest predictor: {names[int(np.argmax(np.abs(beta)))]}")
