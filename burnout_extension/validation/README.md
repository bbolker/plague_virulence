# Numerical validation (R)

This directory contains the reproducible numerical validation of
`../burnout_theta_theory_unified_en.tex`. All production calculations are in R.

## Current production workflow

The current error-map workflow is the dedicated `(K, R0-1)` scan. Run from the
repository root:

```r
Rscript validation/tests/test_theory.R
Rscript validation/run_KR_scan.R
Rscript -e "source('validation/R/plotting.R'); make_current_plots()"
```

This produces `data/K_R0_scan.csv` and the current multi-page PDFs in
`figures/`. Each PDF has five pages for fixed rho values 0.02 through 0.10.
Only eight paper-facing figure sets are retained: first- and second-order entry
errors, first- and second-order action errors, first- and second-order exact-Kendall B errors,
Laplace-only B error, and the selected stochastic P1 comparison.

The broader exact-CTMC validation scan is checkpointed and resumable:

```r
Rscript validation/run_stochastic_scan.R
Rscript validation/plot_stochastic_scan.R
```

It retains aggregate stochastic counts in `data/stochastic_scan_results.csv`,
the resumable state in `data/stochastic_scan_checkpoint.rds`, dense analytical
curves in `data/stochastic_scan_analytic.csv`, and a small set of complete event
histories in `data/stochastic_scan_diagnostic_trajectories.rds`. Figures 9 and
10 show unconditional persistence and persistence conditional on escaping early
fizzle, respectively.

The scan contains 624 parameter combinations with
`rho = {0.01, 0.02, 0.05, 0.10}`, `theta = {0, 0.5, 1}`, four population sizes,
and 13 values of `R0 - 1`.

In Figures 9 and 10, solid analytical curves use the standard matching height
`sqrt(y*/K)`. Dashed curves show the alternative `y*` matching height over its
entire admissible-root domain, allowing the two choices to be compared directly.

## Supporting validation checks

These remain active because they support conclusions in the scientific report:

```r
Rscript validation/run_dense.R          # broad small-rho validity grid
Rscript validation/run_ode_sensitivity.R
Rscript validation/run_stochastic.R
```

Requirements: R 4.5 or later, `deSolve`, `data.table`, and `ggplot2`.
All random experiments use recorded deterministic seeds.

Key files are `R/theory.R`, `R/ode_reference.R`, `R/kendall.R`, `R/scan.R`,
and `R/stochastic.R`. Generated CSVs are under `data/`, figures under
`figures/`, and the scientific interpretation is in `validation_report.md`.
