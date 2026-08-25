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

Figure 11 scans the extended model at `rho = 0.01`, with
`theta = {0, 0.5, 1}`, common population sizes
`K = {10^6, 10^7, 10^8, 10^9}`, and 14 common values of `R0 - 1` from
0.03 to 5:

```r
Rscript validation/run_fig11_scan.R
Rscript validation/plot_fig11_scan.R
```

The stochastic results and resumable state are stored separately as
`data/fig11_scan_results.csv` and
`data/fig11_scan_checkpoint.rds`; dense curves are cached in
`data/fig11_scan_analytic.csv`. Figure 11 contains six pages: unconditional and
conditional-on-not-fizzling probabilities for each value of `theta`. The scan uses adaptive
tau-leaping with tolerance 0.01, 3,000 attempts per point, and 10,000 where the
initial Wilson interval is wider than 0.05.

## Supporting validation checks

These remain active because they support conclusions in the scientific report:

```r
Rscript validation/run_dense.R          # broad small-rho validity grid
Rscript validation/run_ode_sensitivity.R
Rscript validation/run_stochastic.R
```

Requirements: R 4.5 or later, `deSolve`, `data.table`, `ggplot2`, and
`adaptivetau`.
All random experiments use recorded deterministic seeds.

Key files are `R/theory.R`, `R/ode_reference.R`, `R/kendall.R`, `R/scan.R`,
and `R/stochastic.R`. Generated CSVs are under `data/` and figures under
`figures/`.
