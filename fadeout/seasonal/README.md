# Seasonal recurrent-fade-out diagnostics

This workflow diagnoses whether established patches continue to fade out under
seasonal forcing. It uses the existing sinusoidally forced, stochastic
single-strain metapopulation model and keeps strain 2 inactive.

## Operational definition

For the corresponding non-seasonal deterministic single-patch model,

$$
\frac{dS}{dt}=rS\left(1-\frac{S}{K}\right)-\frac{\beta SI}{K},
\qquad
\frac{dI}{dt}=\frac{\beta SI}{K}-\gamma I,
$$

the damped intrinsic frequency and period are

$$
\omega_0=
\sqrt{\frac{\gamma r(R_0-1)}{R_0}-\frac{r^2}{4R_0^2}},
\qquad
T_0=\frac{2\pi}{\omega_0}.
$$

The implementation verifies this frequency against the imaginary part of the
Jacobian eigenvalues. By default, an observed infection episode ending no
later than $1.5T_0$ is an early extinction, while an observed extinction after
$1.5T_0$ is a fade-out. An episode still active at the simulation endpoint is
right-censored, not an observed fade-out. If it has already exceeded
$1.5T_0$, it is included among established episodes.

The multiplier 1.5 is an intentionally approximate operational definition for
detecting post-burnout fade-out. It does not assume that a seasonal trajectory
has clearly identifiable epidemic waves.

## Running the workflow

Validate the example grid without simulating:

```bash
Rscript fadeout/seasonal/run_seasonal_fadeout.R --dry-run
```

Run the small example grid:

```bash
Rscript fadeout/seasonal/run_seasonal_fadeout.R \
  --grid fadeout/seasonal/seasonal_fadeout_example_grid.csv \
  --run-id example
```

Each CSV row is one parameter combination and can specify `R0`, `gamma`, `r`,
`K`, `alpha`, `seasonal_amp`, `season_period`, `peak_day`, `n_patch`, `dt`,
`t_max`, `n_reps`, `base_seed`, and `threshold_multiplier`. Replicate seeds are
explicit and saved. Joint grids over amplitude, `K`, and `alpha` are supported;
the example is deliberately small and is not a parameter scan.

Outputs are written to
`fadeout/output/seasonal_fadeout/<run-id>/`. The `data/` directory contains
episode-level records, annual occupancy, annual fade-out counts,
replicate-level summaries, parameter/seed records, and global-extinction
records. The `figures/` directory contains:

- annual mean occupancy with an exploratory linear fit;
- literal annual maximum-minus-minimum occupancy amplitude;
- a black-and-white patch occupancy raster for every replicate;
- Kaplan-Meier survival of all infection episodes;
- conditional survival after episodes cross the establishment threshold.

All simulation years are retained. The occupancy trend uses complete years
only. If global extinction occurs, the fit stops in the year containing global
extinction so that subsequent zero occupancy is not treated as an ordinary
trend.

## Seasonal vs constant comparison

```bash
Rscript fadeout/seasonal/occupancy_seasonal_vs_constant_100y.R
```

Two 100-year realizations differing only in `seasonal_amp` (0.40 vs 0). The
cosine forcing averages to zero over a year, so both runs share the same mean
transmission and the comparison isolates seasonality rather than a change in
mean. Output is thinned to every 10th day, which makes the one-step
extinction/recolonization counts meaningless; those are not computed here.
Figures and summaries go to `fadeout/output/seasonal_vs_constant_100y/`.

## Pairwise invasion

```bash
Rscript fadeout/seasonal/pairwise_invasion_100y.R [--pilot]
Rscript fadeout/seasonal/pairwise_invasion_100y_plots.R [--pilot]
```

(or `make pip_plots` from the repository root). At K = 5000, seasonal
amplitude 0.2 and alpha = 5e-5/day, where the batch runs show intermediate
occupancy, each resident R0 in {1.2, 1.5, 2, 2.5, 3, 3.5} runs alone for 100
years. Then 10 infecteds of each invader R0 go into one randomly chosen patch
the resident occupies, and both strains are followed for 50 more years (100
replicates per pair, 15 workers by default via `N_WORKERS`). Outputs go to
`fadeout/output/pairwise_invasion_100y/`: per-run and per-pair outcome
summaries, a thinned occupancy series, PIP heatmaps, invader survival curves
and occupancy time series.

`pairwise_invasion_100y.R --add-missing` runs only the grid cells missing from
the existing output (e.g. after adding R0 values) and merges them;
`--susc-eq` starts each patch's S at Uniform(0.7, 1.3) x (S* + I*), capped at
K, and writes to a separate `_susc_eq` directory.

The K × alpha × amplitude factorial version runs on Nibi as a job array
(`odin/sharcnet/pairwise_invasion_factorial_run_array.R`; see "pairwise_invasion_factorial"
in `odin/sharcnet/README.md`). Both versions share
`pairwise_invasion_functions.R`; `pairwise_invasion_factorial_plots.R` draws the
factorial figures from the combined cluster output.

## Initial susceptible density and persistence

```bash
Rscript fadeout/seasonal/s0_persistence_check.R
Rscript fadeout/seasonal/s0_persistence_check_plots.R
```

Single-strain runs in the same regime, with every patch starting at
S = f (S* + I*) (capped at K) and I = I*, for f = 0.5, ..., 1.5 and the eight
pairwise-invasion R0 values (20 runs each, 100 years, stopped at global
extinction). Outputs go to `fadeout/output/s0_persistence_check/`.

## Validation

Run focused synthetic checks:

```bash
Rscript fadeout/seasonal/validate_seasonal_fadeout.R
```

Run the synthetic checks plus the 200-patch, 10-year reference seasonal
scenario used for the original occupancy raster:

```bash
Rscript fadeout/seasonal/validate_seasonal_fadeout.R --smoke
```

The synthetic checks cover early extinction, later fade-out, established
right-censoring, two episodes separated by zero infection, known annual
occupancy summaries, and agreement between the analytic and Jacobian periods.
