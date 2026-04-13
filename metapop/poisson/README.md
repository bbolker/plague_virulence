## Overview

This is a stochastic metapopulation model to ask when **lower** $R_0$ can lead to **more spread and longer persistence** of plague across patches.

The key mechanism is: **between‑patch transmission depends on healthy hosts moving while carrying infected fleas**. Strong within‑patch epidemics (caused by higher $R_0$) kill more hosts, leaving fewer healthy movers, which can reduce export and reseeding.

# simulation outputs

Across replicate simulations, we summarize:

-   **Extinction probability within 1000 steps** (years):\
    fraction of replicates that go extinct by time `t <= 1000`.
-   **Mean extinction time** among extinct runs:\
    average time-to-extinction conditional on extinction.
-   **Infection level among non-extinct runs** (quasi-equilibrium), including
    -   average number of infected patches,
    -   total infections ,
    -   total host population .

# Result

There are still parameter regimes where **lower R0 can favor epidemic persistence** (e.g., lower extinction rate and longer extinction time). However, this effect is no longer universal across all parameter combinations

# Files and links

## Derivations

-   Single-strain stochastic model derivation:
    -   [`single_strain.pdf`](./single_strain.pdf)
-   Deterministic / invasion-condition notes (simplified analytical insight):
    -   [`deterministic_invasion_condition.pdf`](./deterministic_invasion_condition.pdf)

## Simulation scripts

-   Core simulation functions (Poisson introductions):
    -   [`simulation_funs.R`](./simulation_funs.R)
-   Parameter sweep / batch simulation driver:
    -   [`sim.R`](./sim.R)
-   Plotting/report script for batch outputs:
    -   [`plot.R`](./plot.R)

## Outputs

-   Batch summary table:
    -   [`outputs/sim.rds`](./outputs/sim.rds)
-   Raw replicate results (per parameter combination):
    -   [`outputs/raw/`](./outputs/raw/)
-   Plots:
    -   [`outputs/plot.pdf`](./outputs/plot.pdf)
