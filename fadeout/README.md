# Fadeout and patch-occupancy analyses

This directory contains exploratory scripts for local infection fadeout,
recolonization, and metapopulation patch occupancy. The current main workflow is
the non-seasonal, single-strain occupancy exploration and its analytical
comparison. Seasonal and two-strain scripts are retained as related exploratory
work.

Unless noted otherwise, occupancy is defined patch by patch as
`1` when the number of infected hosts is greater than zero (`I > 0`) and `0`
when `I = 0`.

The simulations use functions from the [`plagueMetapop`](../plagueMetapop/)
package. The main non-seasonal stochastic model is
[`euler_odin_def.R`](../plagueMetapop/inst/odin/euler_odin_def.R).

## Current patch-occupancy workflow

Run these commands from the repository root:

```bash
Rscript fadeout/run_occupancy_exploration.R
Rscript fadeout/compare_patch_occupancy_approximation.R
```

The first command runs a small one-factor-at-a-time collection of representative
trajectories. It is not a high-dimensional parameter sweep. The second command
reads those saved trajectories and compares each one with the analytical
post-burnout approximation; it does not rerun the metapopulation simulations.

| File | Purpose |
|---|---|
| [`run_occupancy_exploration.R`](run_occupancy_exploration.R) | Defines the baseline parameter set and the separate `R0`, `K`, `r`, and `alpha` series, runs one stochastic trajectory per parameter value, applies the display-only absorbing boundary at full occupancy, and saves comparison plots and data under `output/patch_occupancy/`. |
| [`occupancy_exploration_functions.R`](occupancy_exploration_functions.R) | Reusable helpers for explicit fadeout initial conditions, running the odin model, calculating patch occupancy and curve summaries, applying the full-occupancy absorbing rule, and making the grouped comparison plots. The current baseline uses virgin-soil initialization: 10 infected hosts and `K - 10` susceptible hosts in every patch. |
| [`patch_occupancy_dynamics.md`](patch_occupancy_dynamics.md) | Short derivation of the logistic post-burnout approximation, with growth rate `lambda = alpha * I_star * P1`. |
| [`compare_patch_occupancy_approximation.R`](compare_patch_occupancy_approximation.R) | Compares all saved one-factor-at-a-time trajectories with the analytical approximation. It estimates `P1` from the existing single-patch extinction data with the GAM used in the project notes, obtains `I_star` from the implemented deterministic model, aligns time zero to the simulated occupancy minimum, and writes PDFs plus pointwise diagnostics to `output/patch_occupancy_approximation/`. |

### Current occupancy outputs

All current fadeout results are kept below [`output/`](output/). New plots for
this work should be saved as PDF.

| File or directory | Contents |
|---|---|
| [`output/patch_occupancy/data/occupancy_results.rds`](output/patch_occupancy/data/occupancy_results.rds) | Complete R object containing parameters, raw summaries, displayed trajectories, and per-run summaries for the occupancy exploration. |
| [`output/patch_occupancy/data/occupancy_curves_absorbing.csv`](output/patch_occupancy/data/occupancy_curves_absorbing.csv) | Long-format trajectories after applying the display-only rule that occupancy remains at 1 after recovery reaches full occupancy. |
| [`output/patch_occupancy/data/curve_summaries.csv`](output/patch_occupancy/data/curve_summaries.csv) | One-row summaries of the trough, trough time, recovery, and stopping information for every parameter value. |
| [`output/patch_occupancy/figures/occupancy_compare_R0.pdf`](output/patch_occupancy/figures/occupancy_compare_R0.pdf) | All `R0` trajectories on one panel, with the other parameters fixed at their baseline values. |
| [`output/patch_occupancy/figures/occupancy_compare_K.pdf`](output/patch_occupancy/figures/occupancy_compare_K.pdf) | All `K` trajectories on one panel. |
| [`output/patch_occupancy/figures/occupancy_compare_r.pdf`](output/patch_occupancy/figures/occupancy_compare_r.pdf) | All host-growth-rate (`r`) trajectories on one panel. |
| [`output/patch_occupancy/figures/occupancy_compare_alpha.pdf`](output/patch_occupancy/figures/occupancy_compare_alpha.pdf) | All between-patch transmission (`alpha`) trajectories on one panel. |
| `output/patch_occupancy/figures/occupancy_compare_{R0,K,r,alpha}.png` | Older PNG copies of the four comparison figures. The PDF files are the preferred outputs. |
| [`output/patch_occupancy_approximation/patch_occupancy_approximation_compare_R0.pdf`](output/patch_occupancy_approximation/patch_occupancy_approximation_compare_R0.pdf) | Simulated and analytical post-burnout curves for the `R0` series. |
| [`output/patch_occupancy_approximation/patch_occupancy_approximation_compare_K.pdf`](output/patch_occupancy_approximation/patch_occupancy_approximation_compare_K.pdf) | Simulated and analytical post-burnout curves for the `K` series. |
| [`output/patch_occupancy_approximation/patch_occupancy_approximation_compare_r.pdf`](output/patch_occupancy_approximation/patch_occupancy_approximation_compare_r.pdf) | Simulated and analytical post-burnout curves for the `r` series. The analytical formula has no direct `r` term, although `r` can affect the simulation and assumptions behind the approximation. |
| [`output/patch_occupancy_approximation/patch_occupancy_approximation_compare_alpha.pdf`](output/patch_occupancy_approximation/patch_occupancy_approximation_compare_alpha.pdf) | Simulated and analytical post-burnout curves for the `alpha` series. |
| [`output/patch_occupancy_approximation/patch_occupancy_approximation_diagnostics.csv`](output/patch_occupancy_approximation/patch_occupancy_approximation_diagnostics.csv) | Parameters, estimated `P1`, equilibrium `I_star`, alignment time, RMSE, maximum absolute deviation, and correlation for every comparison. |

The single-patch data used to estimate `P1` are in
[`odin/sharcnet/outputs/euler_onepatch_onestrain_extinct_logistic_continuous.rds`](../odin/sharcnet/outputs/euler_onepatch_onestrain_extinct_logistic_continuous.rds).
The estimate inherits that experiment's Poisson introduction and 200-generation
extinction horizon, so it is a preliminary proxy for successful persistent
recolonization rather than an exact match to every receiving-patch state.

## Seasonal single-strain workflow

| File | Purpose |
|---|---|
| [`seasonal_model_metapop.R`](seasonal_model_metapop.R) | Stand-alone odin model definition for a stochastic Euler metapopulation with synchronous sinusoidal forcing of transmission. Its between-patch force of infection has the same mean-field normalization as the non-seasonal model. |
| [`seasonal_single_strain_metapop.R`](seasonal_single_strain_metapop.R) | Runs a ten-year, single-strain seasonal metapopulation simulation from explicit initial conditions; calculates occupancy, extinction/recolonization, seasonal `R0`, and patch population summaries; and writes the results listed below. |

Files currently in [`output/single_strain/`](output/single_strain/) are generated
by `seasonal_single_strain_metapop.R`:

| File | Contents |
|---|---|
| [`trajectory.rds`](output/single_strain/trajectory.rds) | Converted full stochastic trajectory. |
| [`patch_state.csv`](output/single_strain/patch_state.csv) | Long-format susceptible, infected, population, and occupancy states by patch and time. |
| [`metapop_summary.csv`](output/single_strain/metapop_summary.csv) | Time series of metapopulation totals, occupancy, turnover events, seasonal transmission, and seasonal `R0`. |
| [`patch_summary.csv`](output/single_strain/patch_summary.csv) | Patch-level extinction, recolonization, and occupancy summaries. |
| [`overall_summary.csv`](output/single_strain/overall_summary.csv) | Overall parameter and trajectory summary. |
| [`occupancy_raster.pdf`](output/single_strain/occupancy_raster.pdf) | Patch-by-time occupancy raster. |
| [`occupied_patches.pdf`](output/single_strain/occupied_patches.pdf) | Number of infected patches through time. |
| [`occupancy_seasonal_R0.pdf`](output/single_strain/occupancy_seasonal_R0.pdf) | Occupancy overlaid with seasonally varying `R0`. |
| [`extinction_recolonization.pdf`](output/single_strain/extinction_recolonization.pdf) | Local extinction and recolonization counts through time. |
| [`mean_patch_population.pdf`](output/single_strain/mean_patch_population.pdf) | Mean host population per patch through time. |

## Other exploratory scripts

| File | Purpose and output location |
|---|---|
| [`single_strain_metapop.R`](single_strain_metapop.R) | Earlier non-seasonal single-trajectory analysis with explicit susceptible and infected initial values, occupancy, and local turnover. It currently writes to `odin/outputs/onestrain_rescue/`, not to `fadeout/output/`. Its initialization and occupancy calculations were the template for the current occupancy workflow. |
| [`two_strain_metapop.R`](two_strain_metapop.R) | Two-strain rescue/invasion example: strain 1 begins near its endemic equilibrium and strain 2 is introduced later into one patch. It analyzes four patch states, strain-specific occupancy, prevalence, extinction, and recolonization, and writes to `odin/outputs/twostrain_rescue/`. |
| [`stochastic_SIR.R`](stochastic_SIR.R) | Single-patch demographic SIR prototype comparing a Gillespie stochastic trajectory with the corresponding deterministic solution and endemic equilibrium. It displays plots interactively and does not save files. |
| [`stochastic_seasonal.R`](stochastic_seasonal.R) | Single-patch extension using piecewise-constant high- and low-risk seasons; compares stochastic Gillespie and deterministic seasonal dynamics. It displays plots interactively and does not save files. This is separate from the sinusoidally forced metapopulation model above. |

## Maintenance convention

When a fadeout script, derived data file, or figure is added, removed, renamed, or
changes purpose, update this README in the same change. Keep generated results
under `fadeout/output/`, prefer PDF for figures, and avoid overwriting existing
simulation data unless regeneration is explicitly intended.
