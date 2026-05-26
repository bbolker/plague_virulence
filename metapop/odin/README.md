## odin (and other) subdirectory

This directory implements plague virulence simulations using [odin](https://mrc-ide.github.io/odin/) (and other platforms). (The directory should probably be renamed during cleanup.)

### Model definition files (odin DSL)

| File | Description |
|------|-------------|
| `discrete_odin_def.R` | Stochastic discrete-time two-strain n-patch SIR model (mass-action, binomial/Poisson draws, logistic vital dynamics); `I2_ini[]` allows strain-2 initial condition to be set on chunked restarts |
| `euler_odin_def.R` | Stochastic continuous-time version: adds per-strain recovery rate `gamma[]` and time step `dt`; transitions use exponential-hazard probabilities scaled by `dt`; also carries `I2_ini[]` |
| `euler_det_odin_def.R` | Deterministic version of `euler_odin_def.R`: replaces `rpois`/`rbinom` draws with their expectations; state variables are real-valued; also carries `I2_ini[]` |
| `odin_twostrain0.R` | Earlier single-patch two-strain odin model definition (continuous-time SIR with explicit `mu`); superseded by `discrete_odin_def.R` |
| `discrete_onepatch_odin.R` | Early prototype: single-patch discrete-time SIR with logistic growth; superseded |

### Platform infrastructure

| File | Description |
|------|-------------|
| `discrete_odin.R` | `compile_odin`, `make_simulator_odin` (stores generator and init args as attributes for chunked restarts), `run_simulator_odin` (chunked early-stopping via per-chunk reinitialisation), `conv_odin`, and stopping-condition helpers (`stop_either_extinct`, `stop_both_extinct`) |
| `discrete_macpan2.R` | Equivalent make/run/conv functions for the macpan2 platform |
| `discrete_pureR.R` | Pure-R reference implementation of the same model (no compiled backend) |
| `discrete_run.R` | Top-level `discrete_run()` dispatcher: compiles the chosen platform, parallelises over `nsim` simulations via furrr/future, and provides `sumfun_discrete()` for summarising results; default `stop_cond = stop_both_extinct` halts odin runs as soon as both strains are globally extinct |

### Testing and validation

| File | Description |
|------|-------------|
| `compare_platforms.R` | Compares odin, macpan2, and pureR outputs for consistency; also times simulator creation and trajectory generation |
| `tinytest/tests.R` | Automated tests via tinytest: correctness checks for each platform, parallel execution, and summary statistics |

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
| `odin_twostrain_run0.R` | Early two-strain run script paired with `odin_twostrain0.R`; superseded by `discrete_run.R` framework |
| `odin_twostrain_run.R` | Intermediate two-strain run script; superseded by `discrete_run.R` framework |
| `odin_invasion_run.R` | Invasion experiment: resident run to quasi-equilibrium then invader seeded (uses `odin_twostrain0.R`) |
| `discrete_onespecies_sim.R` | Single-species (no strains) simulation; exploratory |
| `tmp.R` | Scratch/exploratory code; not part of the main workflow |
