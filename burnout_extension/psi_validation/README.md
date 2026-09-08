# Numerical validation for psi = 0.5

This directory reproduces the boundary-layer-independent comparisons corresponding
to Figures 12, 13, and 15 for the extended incidence model in
`../burnout_psi_theory_inverseflow.tex`.

Run from the repository root:

```r
Rscript psi_validation/tests/test_psi.R
Rscript psi_validation/run_scan.R
Rscript psi_validation/plot_figures.R
```

The full 624-point grid uses the exact CTMC for K from 1,000 through 30,000,
with the same rho, theta, and R0 grid as the original Figures 12/13. Peaks and
troughs are classified using the state-dependent growth rate
`R0*x*(x+y)^(-psi)-1`.
