# overview of plague metapopulation model

## three dimensional model

We have a discrete-time model keeping track of three state variables: fraction of infected patches $p$, average population density of infected patches $N_I$ and average population density of susceptible patches $N_S$.

There are several events in each iteration that can change the state of the system: movement of hosts and infection of susceptible patches; death of infected hosts and burnout of infected patches; and logistic growth of hosts in both susceptible and infected patches.

Detailed description of the model in [3d.pdf](./3d.pdf) and latex source [3d_latex_source](./3d_latex_source). Simulation code in [simulation_fun.R](./simulation_fun.R).

Simulation indicates that population density of susceptible patches doesn't change a lot over time, so we simplify the model by assuming that $N_S$ is constant. This leads to a two-dimensional model with state variables $p$ and $N_I$.

## two dimensional model -- fixed population density of susceptible patches

In the two-dimensional model, we have only two state variables: fraction of infected patches $p$ and average population density of infected patches $N$. 

The equations are much less complicated than the three-dimensional one, and we get everything analytically in this version. We have now got the disease free equilibrium (DFE), the endemic equilibrium (EE) and the threshold condition for the stability of DFE.

Detailed results in [fixed_ns](./fixed_ns.pdf), source file [fixedns.tex](./fixedns.tex) and simulation code in [fiexed_ns.R](./fixed_ns.R).
