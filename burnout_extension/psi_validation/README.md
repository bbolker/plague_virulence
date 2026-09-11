# Numerical validation for psi = 0.25, 0.5, 0.75, 1

This directory reproduces the boundary-layer-independent comparisons corresponding
to Figures 12, 13, and 15 for the extended incidence model in
`../burnout_psi_theory_inverseflow.tex`.

Run from the repository root:

```r
Rscript psi_validation/tests/test_psi.R
Rscript psi_validation/run_scan.R
Rscript psi_validation/plot_figures.R
```

`run_scan.R` loops over `psi = 0.25, 0.5, 0.75, 1`, running the same 624-point
grid for each and writing independent, resumable checkpoint/output files per
value (`psi_validation/data/psi<label>_scan_checkpoint.rds` and
`psi_validation/data/psi<label>_stochastic_results.csv`, with `<label>` the
decimal point stripped from the psi value, e.g. `psi025`, `psi05`, `psi075`,
`psi1`). The downstream analysis scripts (`plot_figures.R`,
`summarize_results.R`, and the nonlinear-tail scripts below) currently only
consume the `psi05` outputs.

By default a checkpoint (results CSV rewrite + resumable state save) is
written every 1,000 batches; pass `--checkpoint-every=N` to change that, e.g.:

```r
Rscript psi_validation/run_scan.R --checkpoint-every=50
```

The full 624-point grid uses the exact CTMC for K from 1,000 through 30,000,
with the same rho, theta, and R0 grid as the original Figures 12/13. Peaks and
troughs are classified using the state-dependent growth rate
`R0*x*(x+y)^(-psi)-1`.

## Timing and OpenMP threads

`R/ctmc_psi.cpp` parallelizes each batch of attempts internally with OpenMP
(`#pragma omp parallel for`); `run_scan.R` itself does not call
`omp_set_num_threads()`, so the thread count follows whatever `OMP_NUM_THREADS`
is set to in the environment the script is launched from (or the OpenMP
runtime default if unset, which is not reliable on machines where `nproc`
under-reports the available cores). Set it explicitly, e.g.:

```sh
OMP_NUM_THREADS=28 Rscript psi_validation/run_scan.R
```

With `OMP_NUM_THREADS=28`, all four psi scans (0.25, 0.5, 0.75, 1) together
took about 7.5 minutes wall-clock time on a 32-core machine.

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
