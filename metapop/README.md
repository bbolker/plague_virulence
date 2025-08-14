# overview of plague metapopulation model

## three dimensional model

We have a discrete-time model keeping track of three state variables: fraction of infected patches $p$, average population density of infected patches $N_I$ and average population density of susceptible patches $N_S$.

There are several events in each iteration that can change the state of the system: movement of hosts and infection of susceptible patches; death of infected hosts and burnout of infected patches; and logistic growth of hosts in both susceptible and infected patches.

Detailed description of the model in the documents: [pdf](./3d.pdf), [latex source (zip file)](./3d.zip) and [simulation code](./simulation_fun.R).

Simulation indicates that population density of susceptible patches doesn't change a lot over time, so we simplify the model by assuming that $N_S$ is constant. This leads to a two-dimensional model with state variables $p$ and $N_I$.

## two dimensional model -- fixed population density of susceptible patches

In the two-dimensional model, we have only two state variables: fraction of infected patches $p$ and average population density of infected patches $N$.

The equations are much less complicated than the three-dimensional one, and we get everything analytically in this version. We have now got the disease free equilibrium (DFE), the endemic equilibrium (EE) and the threshold condition for the stability of DFE.

Detailed results in [pdf](./fixedns.pdf), [source (zip file)](./fixedns.zip) and [simulation](./fixed_ns.R).

# Update

Updates since last time:

-   replaced original infection probability with hazard model
-   discussed order of events (non-cyclic permutation makes slight difference)
-   find a scenario where reducing death probability $D$ makes originally stable DFE unstable
-   get an estimate of R0 (not very confident)

Some challenges now:

-   stability of EE
-   smaller $R_0$ decreases persistence probability in our current model (that seems not to be what we expected)
-   estimation of other parameters ($B$,$c$,$r$ seems related to length of a time step, $S$ is related to patch size), and the "definition" of $\varepsilon$ (should that include fleas?)
