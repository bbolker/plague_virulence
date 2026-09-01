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
Within each error map, facet rows compare `sqrt(y*/K)`, the `2/3 compromise`,
the `3/4 compromise`, and `y*` using the same scales and missing-value treatment.
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

The boundary-layer-independent (BI) leading-order prediction is evaluated on
those same cached points with:

```r
Rscript validation/plot_BI_validation.R
```

This does not rerun stochastic simulations. It writes separate unconditional
and conditional curve PDFs to
`figures/fig12_BI_unconditional_stochastic_validation.pdf` and
`figures/fig13_BI_conditional_stochastic_validation.pdf`.
BI has no matching-height argument and remains defined wherever the
finite-boundary curves stop.

The script also writes two conditional display alternatives without replacing
the linear figure. `fig12_BI_stochastic_validation_cond_logit.pdf` uses
probability-labelled logit coordinates with explicit clipping at 0.001 and
0.999; `fig12_BI_stochastic_validation_cond_free.pdf` uses a separate linear
y range in each K panel. Absolute probability error remains the numerical
accuracy metric; these transformations affect visualization only.

The scan contains 624 parameter combinations with
`rho = {0.01, 0.02, 0.05, 0.10}`, `theta = {0, 0.5, 1}`, four population sizes,
and 13 values of `R0 - 1`.

In Figures 9 and 10, analytical curves compare `sqrt(y*/K)`, `y*`, and the
intermediate **2/3 compromise** `K^(-1/3) * y*^(2/3)` and **3/4 compromise**
`K^(-1/4) * y*^(3/4)` over each choice's entire admissible-root domain.

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

For legibility, Figures 9, 10, and 11 are directories under `figures/`, each
containing one PDF per matching height. They use identical axes, analytical
methods, stochastic points, and uncertainty intervals; no combined overlaid
version is produced.

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

## First-trough boundary-entry diagnostic

`R/trough_diagnostics.R` measures the first post-epidemic deterministic trough
at the first upward crossing of the infective nullcline `x = 1/R0`. The crossing
time is root-refined by the ODE solver rather than read from the output grid.
It compares the existing four matching heights with two configurable diagnostic
inverse-log candidates, `c_log/log(K)` and `c_log*y_star/log(K)`. These candidates
do not replace the matching height used by the burnout calculation.

```r
Rscript validation/run_trough_diagnostics.R          # representative smoke grid
Rscript validation/run_trough_diagnostics.R --full   # existing K/R0 grid
Rscript validation/plot_trough_diagnostics.R validation/data/first_trough_diagnostics.csv
```

The CSV records trough depth, nullcline residual, boundary ratio, actual downward
entry coordinate, and first/second-order phase-plane trough approximations. The
literal `c_log/log(K)` option need not lie below `y_star`; it is retained as a
clearly labeled interpretation to expose, rather than hide, that ambiguity.
