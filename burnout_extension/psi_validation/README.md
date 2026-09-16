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

## Nonlinear-tail theory update

For `../burnout_psi_nonlinear_tail_theory.tex`, reuse the same stochastic CSV
and rebuild only the analytical curves and figures:

```r
Rscript psi_validation/build_nonlinear_tail_curves.R
Rscript psi_validation/plot_nonlinear_tail_figures.R
Rscript psi_validation/summarize_nonlinear_tail.R
```

The nonlinear continuation retains `(1+y/x)^(-psi)` and selects the geometric
handoff when it is reached before the safe pre-trough cutoff. Otherwise it uses
the deepest available point with `g < 0.98`. Curves are left undefined when the
second-order regular-outer initialization is not on the descending branch.

The rebuilt curve cache also contains the finite-prevalence Kendall prediction
from the updated theory.  The deterministic system is continued in time to the
first recovery trough, located by a continuous `g = 1` event, and records
`x_t`, `y_t`, `log_y_t`, `s_t = y_t/x_t`, and the analytic crossing speed
`alpha_t`. The log prevalence is retained so exponentially deep troughs remain
numerically meaningful even when `y_t` itself underflows. The
new `trough` column is computed from
`B_tr = K*y_t*sqrt(alpha_t/(2*pi))`.  Figures compare this refinement against
the same CTMC cache; stochastic trajectories do not need to be regenerated.

For a compact conditional overview similar to the original `fig13_comb.pdf`:

```r
Rscript psi_validation/plot_fig13_comb_nonlinear.R
```
