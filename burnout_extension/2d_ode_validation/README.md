# Full 2D ODE validation for the finite-prevalence trough theory

This directory validates `../burnout_finite_prevalence_trough_theory.tex` for
`psi = 0.5`.  It directly integrates the full nonlinear two-dimensional ODE
from `(x(0), y(0)) = (1 - 1/K, 1/K)`, detects the first infective peak, physical
trough, and recovery peak as continuous events, and evaluates both:

- the closed physical-trough saddle formula
  `B_tr = K*y_t*sqrt(alpha_t/(2*pi))`;
- the direct finite-horizon Kendall integral on the same deterministic path.

The CTMC estimates are reused from
`../psi_validation/data/psi05_stochastic_results.csv`. No stochastic simulation
is rerun by this workflow.

Run from the repository root:

```r
Rscript 2d_ode_validation/tests/test_2d_ode.R
Rscript 2d_ode_validation/build_curves.R
Rscript 2d_ode_validation/summarize_results.R
Rscript 2d_ode_validation/plot_figures.R
Rscript 2d_ode_validation/plot_comb.R
```

The validation points are the same 624-point grid as `psi_validation`: four
values of `rho`, three values of `theta`, four population sizes, and thirteen
values of `R0`. Smooth analytical curves additionally use the same dense
101-point log-spaced `R0 - 1` grid as the earlier nonlinear-tail figures.
Figures 12 and 13 compare the reused CTMC estimates with both
semi-analytical predictions. Figure 15 checks the closed trough reduction
against direct Kendall quadrature on the 48-point theory grid. The condensed
figure uses colour for `rho`, line type/point shape for `theta`, heavy curves for
direct Kendall quadrature, and no closed-trough overlay.

The corresponding `psi = 1` workflow is provided by the scripts whose names
end in `_psi1.R`. Its CTMC results are newly simulated and checkpointed in
`data/psi1_scan_checkpoint.rds`; they are not reused from the earlier scan.
