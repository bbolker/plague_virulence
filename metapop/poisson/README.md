# Overview

This subdirectory contains an updated stochastic metapopulation implementation that modifies two key assumptions used in earlier versions of the model:

1.  **Initial number of infected hosts per newly infected patch**
    -   **Old version:** every newly infected patch starts with exactly **1** infected host.
    -   **New version:** the number of initially infected hosts is **Poisson-distributed**
2.  **Burnout / establishment mechanism**
    -   **Old version:** used a burnout probability (assume some infected rats live until next epidemic season)
    -   **New version:** assume all infected hosts die and infected fleas introduce the next wave

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
