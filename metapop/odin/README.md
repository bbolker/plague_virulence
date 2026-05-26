## odin (and other) subdirectory

This directory contains run and plot scripts for plague virulence simulations.
Core simulation infrastructure (model definitions, platform backends, dispatcher)
lives in the [`plagueMetapop`](../plagueMetapop/) R package; all scripts load it
via `library(plagueMetapop)`.

### plagueMetapop package ([`../plagueMetapop/`](../plagueMetapop/))

| Component | Description |
|-----------|-------------|
| `R/discrete_odin.R` | `compile_odin`, `make/run/conv_simulator_odin`, chunked early-stopping, `stop_either/both_extinct` |
| `R/discrete_macpan2.R` | Equivalent make/run/conv functions for the macpan2 platform |
| `R/discrete_pureR.R` | Pure-R reference implementation |
| `R/discrete_run.R` | Top-level `discrete_run()` dispatcher and `sumfun_discrete()` summary function |
| `inst/odin/discrete_odin_def.R` | Stochastic discrete-time two-strain n-patch SIR model |
| `inst/odin/euler_odin_def.R` | Stochastic continuous-time version with per-strain recovery and `dt` |
| `inst/odin/euler_det_odin_def.R` | Deterministic version of `euler_odin_def.R` |
| `inst/tinytest/test_platforms.R` | Automated tests: correctness checks for each platform, parallel execution, summary statistics |

### Testing and validation

| File | Description |
|------|-------------|
| `compare_platforms.R` | Compares odin, macpan2, and pureR outputs for consistency; also times simulator creation and trajectory generation |

Run package tests with `tinytest::test_package("plagueMetapop")`.

### HPC job arrays (Alliance Canada / SLURM) ([`sharcnet/`](sharcnet/))

Each pair of scripts runs the corresponding grid on Compute Canada using SLURM job arrays.
Submit from the `metapop/odin/` directory after `mkdir -p sharcnet/logs sharcnet/outputs`.

| File | Description |
|------|-------------|
| `sharcnet/submit_euler_extinct.sh` | 280-task array (1 per grid point) for the euler one-strain extinction grid; 4 CPUs/task, 30 min |
| `sharcnet/euler_onepatch_onestrain_extinct_run_array.R` | Array R script: reads `SLURM_ARRAY_TASK_ID`, runs one (R0, K) point with `nsim=100` |
| `sharcnet/euler_onepatch_onestrain_extinct_combine.R` | Combines per-task `.rds` files into `euler_onepatch_onestrain_extinct.rds` |
| `sharcnet/submit_twostrain_extinct.sh` | 1000-task array (batched) for the two-strain extinction grid (25 921 points → ~26 rows/job); 1 CPU/task, 1 h |
| `sharcnet/discrete_onepatch_twostrain_extinct_run_array.R` | Array R script: reads task ID and `--stepR0`/`--nsim`/`--njobs` options, runs one batch of grid rows |
| `sharcnet/discrete_onepatch_twostrain_extinct_combine.R` | Combines per-task `.rds` files into `discrete_onepatch_twostrain_extinct.rds` |

### Run and plot scripts (current)

| File | Description |
|------|-------------|
| `discrete_onepatch_onestrain_example.R` | Single-patch, one-strain example run using the discrete stochastic model |
| `discrete_onepatch_onestrain_extinct_run.R` | Grid run over R0 × K to compute extinction probability and time for one strain (discrete model) |
| `discrete_onepatch_onestrain_extinct_plot.R` | Raster plots of extinction probability and mean extinction time from the above grid |
| `discrete_onepatch_twostrain_extinct_run.R` | Grid run over R0₁ × R0₂ for two-strain invasion/coexistence analysis (discrete model) |
| `discrete_onepatch_twostrain_extinct_plot.R` | Plots of two-strain extinction results |
| `euler_onepatch_onestrain_example.R` | Single-patch example using the continuous-time (`euler_odin_def.R`) model; checks S monotonicity |
| `euler_onepatch_onestrain_extinct_run.R` | Grid run over R0 × K using the continuous-time stochastic model |

### Run and plot scripts (older/exploratory)

| File | Description |
|------|-------------|
| `odin_twostrain0.R` | Earlier single-patch two-strain odin model definition; superseded by `discrete_odin_def.R` |
| `discrete_onepatch_odin.R` | Early prototype: single-patch discrete-time SIR; superseded |
| `odin_twostrain_run0.R` | Early two-strain run script; superseded by `plagueMetapop` framework |
| `odin_twostrain_run.R` | Intermediate two-strain run script; superseded by `plagueMetapop` framework |
| `odin_invasion_run.R` | Invasion experiment: resident run to quasi-equilibrium then invader seeded |
| `discrete_onespecies_sim.R` | Single-species (no strains) simulation; exploratory |
