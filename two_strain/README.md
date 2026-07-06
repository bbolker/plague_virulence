# Overview

This subdirectory contains a **two-strain (resident + invader) stochastic metapopulation** implementation, extending the **Poisson-introduction** logic used in `metapop/poisson/` to a **competitive two-strain** setting.

We let the **resident strain** run for **100 years** to generate a quasi-stationary spatial epidemic background, then introduce an **invader strain** (with slightly lower $R_0$) and continue the simulation for an additional **400 years**.

If both strains establish in a patch in the same year, the model assigns an effective reproduction number using an average (`R0_mean <- (R01 + R02)/2`), then **partitions total epidemic mortality** between strains in proportion to the Poisson seed counts (allocation based on `a/(a+b)` vs `b/(a+b)`).

# Measurement

Across replicate simulations, we summarize per strain:

-   **Extinction probability** (fraction of replicates that go extinct) in 500 years
-   **Mean extinction time** among extinct runs
-   **Infection level among non-extinct runs** (quasi-equilibrium), including
    -   average number of infected patches
    -   total infections (mortality totals)
    -   total host population (reported for convenience for both strains)

To indicate the invasion capability of the invader and construct PIP, we use the **persistence probability**(fraction of replicates that don't go extinct in 500 years) if there is a considerable number of simulations that don’t go extinct (\> 10%), and otherwise use **mean extinction time**.

# Results

-   An example where invader strain establishes successfully:
    -   [`two_strain_test.pdf`](./two_strain_test.pdf)
-   Simulations over larger parameter sets (Invader R0 = Resident R0 − 0.2):
    -   [`outputs_2strain/plot_2strai_points.pdf`](./outputs_2strain/plot_2strai_points.pdf)
-   PIP: [`outputs_pip/pip_hybrid.pdf`](./outputs_pip/pip_hybrid.pdf)

# Files and links

## Derivations

-   Two-strain notes/derivations:
    -   [`two_strain.pdf`](./two_strain.pdf)
    -   LaTeX source:
        -   [`two_strain.tex`](./two_strain.tex)

## Simulation scripts

-   Core simulation function:
    -   [`sim_fun.R`](./sim_fun.R)
    
-   Batch simulation:
    -   [`sim_2strain.R`](./sim_2strain.R)
-   Plotting scripts:
    -   [`plot_2strain.R`](./plot_2strain.R)

-   PIP simulation:
    -   [`sim_pip.R`](./sim_pip.R)
-   PIP plotting:
    -   [`plot_pip.R`](./plot_pip.R)
-   Exploring two-strain ODE outcomes:
    - run ODEs, [`ode_test.R`](./ode_test.R)
    - compare approximations, [`ode_test_proc.R`](./ode_test_proc.R)
