# overview of plague metapopulation model

## three dimensional model

We have a discrete-time model keeping track of three state variables: fraction of infected patches $p$, average population density of infected patches $N_I$ and average population density of susceptible patches $N_S$. JD: This means we are neglecting variation in population density within each group; this is an assumption we need to return to.

There are several events in each iteration that can change the state of the system: movement of hosts and infection of susceptible patches; death of infected hosts and burnout of infected patches; and logistic growth of hosts in both susceptible and infected patches.

JD: I am also a bit worried about the implications of doing everything in turns in our discrete-time model. In particular, this implies that a within-patch epidemic that burns out at the most natural time cannot possibly transmit, which seems unlikely. The simplest alternative (although it may be too extreme) would be to make transmission for next season proportional to the size of the epidemic in a patch. This undercuts our mechanism substantially (because transmission would always increase with R0) – but it's an open question whether it completely eliminates the mechanism (because persistence in the patch could still decrease with R0). Instead of our current landscape reproductive number L0 = (1+B)p_1, we would have something like p_1 + p_nf\*Z, where p_nf=1-1/R0 is the “non-fizzle” probability and Z is the final size. We could start by asking whether this is ever non-monotonic, which I think is not obvious. We could also think about assumptions somewhere between the current (no within-season transmission) and this simple alternative (perfect proportionality).

Detailed description of the model in the documents: [pdf](./3d.pdf), [latex source (zip file)](./3d.zip) and [simulation code](./simulation_fun.R).

Simulation indicates that population density of susceptible patches doesn't change a lot over time, so we simplify the model by assuming that $N_S$ is constant. This leads to a two-dimensional model with state variables $p$ and $N_I$.

## two dimensional model -- fixed population density of susceptible patches

In the two-dimensional model, we have only two state variables: fraction of infected patches $p$ and average population density of infected patches $N$.

The equations are much less complicated than the three-dimensional one, and we get everything analytically in this version. We have now got the disease free equilibrium (DFE), the endemic equilibrium (EE) and the threshold condition for the stability of DFE.

Detailed results in [pdf](./fixedns.pdf), [source (zip file)](./fixedns.zip) and [simulation](./fixed_ns.R).

# Update

Updates since last time:

-   find a scenario where reducing death probability $D$ makes originally stable DFE unstable
-   incorporated change of $N$ in [P1-R0 figure](./P1atDFE.R) (looks similar to figure 5)

Some challenges now:

-   stability of EE
-   Parameters are not likely to be in the region where decreasing $R_0$ increases persistence
-   estimation of other parameters ($B$,$c$,$r$ seems related to length of a time step, $S$ is related to patch size), and the "definition" of $\varepsilon$ (should that include fleas?)

Current ideas:

-   try to find a scenario where decreasing $D$ destabilizes DFE with parameters in reasonable range
-   try some different models that enables infected patches infect other patches before they burn out. (For example, make transmission for next season proportional to the size of the epidemic in a patch, or burning-out patches always survive current iteration)
