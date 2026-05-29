## odin (and other) subdirectory

This directory contains run and plot scripts for plague virulence simulations.
Core simulation infrastructure (model definitions, platform backends, dispatcher)
lives in the [`plagueMetapop`](../plagueMetapop/) R package; all scripts load it
via `library(plagueMetapop)`.

There is a sharcnet subdirectory that BB is leaning on right now.

## 2026 May 27 (Wed)

the current simulations are at the scale of disease generations (or smaller). We have given up for now on the annual (or pseudo-annual) time scale because of conceptual problems in linking scales (particularly rat and flea movement).

Using an Euler step with “some” hazard calculations. Started with Reed-Frost (meaning fixed-length generations; everyone recovers or dies at the end of each time step). Things that start with discrete (rather than Euler) are left over from that.

Reed-Frost leads to _lots_ of burnout: JD suggested possibly leaning in to this, since it might lead to more coexistence tradeoffs. Bellman-Harris and Crump–Mode–Jagers (CMJ) are more general words for “branching process” without the exponential assumption. Back now to geometric distribution (~5-10 Euler steps per disease generation.

An interesting observation: replacing “leaky bucket” (standard replenishment) dynamics with logistic (survivor-based replenishment) is expected to favor relatively more burnout at higher R values. This is apparently reflected in [a burnout plot you can make](euler_onepatch_onestrain_extinct.png). 

Also true that ε is much larger here (Ben says ~ 0.02), but the logisticity seems like a larger problem. Todd asks whether we can examine what the rat troughs look like across these simulations.

Ben wants to do pairwise invasion plots, panelled by α (patch-linkage) and K (rat carrying capacity per patch).

Right now, we're doing mass-action (in JD terms B(N) = βN/K). Should consider other functional forms (one example would be β(N/K)̂^φ).

### plagueMetapop package ([`../plagueMetapop/`](../plagueMetapop/))

| Component | Description |
|-----------|-------------|
| `R/discrete_odin.R` | `compile_odin`, `make/run/conv_simulator_odin`, chunked early-stopping, `stop_either/both_extinct` |
| `R/discrete_macpan2.R` | Equivalent make/run/conv functions for the macpan2 platform |
| `R/discrete_pureR.R` | Pure-R reference implementation |
| `R/discrete_run.R` | Top-level `discrete_run()` dispatcher and `sumfun_discrete()` summary function |
| `inst/odin/discrete_odin_def.R` | Stochastic discrete-time two-strain n-patch SIR model |
| `inst/odin/euler_odin_def.R` | Stochastic continuous-time version with per-strain recovery, `dt`, `logistic_growth` (1=logistic, 0=linear restoring force), and `reedfrost` (1=100% removal per step) |
| `inst/odin/euler_det_odin_def.R` | Deterministic version of `euler_odin_def.R` |
| `inst/tinytest/test_platforms.R` | Tests: platform correctness, parallel execution, summary statistics |
| `inst/tinytest/test_stop_conditions.R` | Tests: `stop_either_extinct()` factory and pre-seeding behaviour |
| `inst/tinytest/test_discrete-vs-euler.R` | Tests: euler with `reedfrost=1, gamma=1, dt=1` matches discrete model |

Run package tests with `tinytest::test_package("plagueMetapop")`.

### Testing and validation

| File | Description |
|------|-------------|
| `compare_platforms.R` | Compares odin, macpan2, and pureR outputs for consistency; times simulator creation and trajectory generation |

### HPC job arrays (Alliance Canada / SLURM) ([`sharcnet/`](sharcnet/))

See [`sharcnet/README.md`](sharcnet/README.md) for full details on all jobs,
grids, parameters, and submission workflow.

Submit all jobs from the `sharcnet/` directory after `mkdir -p logs outputs`.

### Run and plot scripts (current)

`discrete` is the Reed-Frost model, `euler` is the continuous (ish) model (time scale = disease generation time; `gamma=1` wlog; multiple steps per disease generation). * indicates there's a corresponding batch run in the `sharcnet` subdir

| File | Description |
|------|-------------|
| `discrete_onepatch_onestrain_example.R` | Single-patch, one-strain example run using the discrete stochastic model |
| `discrete_onepatch_onestrain_extinct_run.R` | Grid run over R0 × K to compute extinction probability and time (discrete model) |
| `discrete_onepatch_onestrain_extinct_plot.R` | Raster plots of extinction probability and mean extinction time |
| `discrete_onepatch_twostrain_extinct_run.R` | Grid run over R0₁ × R0₂ for two-strain invasion/coexistence analysis |
| `discrete_onepatch_twostrain_extinct_plot.R` | Plots of two-strain extinction results |
| `euler_onepatch_onestrain_example.R` | Single-patch example using the continuous-time (`euler_odin_def.R`) model |
| `euler_onepatch_onestrain_extinct_run.R` | Grid run over R0 × K; `--lineargrowth` (linear demography), `--reedfrost` (Reed-Frost dynamics), `--mini` flags * |
| `euler_onepatch_onestrain_extinct_plot.R` | Raster plots of extinction results from the euler one-strain grid |
| `euler_onestrain_run.R` | Multi-patch grid run over R0 × K × alpha; supports `--mini` flag via optparse  * |
| `odin_twostrain_run.R` | Two-strain run script; grid over R01 x R02 x K x alpha *: `meta_euler_twostrain`|

### Run and plot scripts (older/exploratory)

| File | Description |
|------|-------------|
| `odin_twostrain0.R` | Earlier single-patch two-strain odin model definition; superseded by `discrete_odin_def.R` |
| `discrete_onepatch_odin.R` | Early prototype: single-patch discrete-time SIR; superseded |
| `odin_twostrain_run0.R` | Early two-strain run script; superseded by `plagueMetapop` framework |
| `odin_invasion_run.R` | Invasion experiment: resident run to quasi-equilibrium then invader seeded |
| `discrete_onespecies_sim.R` | Single-species (no strains) simulation; exploratory |
