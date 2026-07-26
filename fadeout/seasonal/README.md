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
