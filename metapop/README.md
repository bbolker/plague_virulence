# overview of plague metapopulation model

## three dimensional model

We have a discrete-time model keeping track of three state variables: fraction of infected patches $p$, average population density of infected patches $N_I$ and average population density of susceptible patches $N_S$. JD: This means we are neglecting variation in population density within each group; this is an assumption we need to return to.

There are several events in each iteration that can change the state of the system: movement of hosts and infection of susceptible patches; death of infected hosts and burnout of infected patches; and logistic growth of hosts in both susceptible and infected patches.

JD: I am also a bit worried about the implications of doing everything in turns in our discrete-time model. In particular, this implies that a within-patch epidemic that burns out at the most natural time cannot possibly transmit, which seems unlikely. The simplest alternative (although it may be too extreme) would be to make transmission for next season proportional to the size of the epidemic in a patch. This undercuts our mechanism substantially (because transmission would always increase with R0) – but it's an open question whether it completely eliminates the mechanism (because persistence in the patch could still decrease with R0). Instead of our current landscape reproductive number L0 = (1+B)p_1, we would have something like p_1 + p_nf\*Z, where p_nf=1-1/R0 is the “non-fizzle” probability and Z is the final size. We could start by asking whether this is ever non-monotonic, which I think is not obvious. We could also think about assumptions somewhere between the current (no within-season transmission) and this simple alternative (perfect proportionality).

Detailed description of the model in the [documents](./3d)

Simulation indicates that population density of susceptible patches doesn't change a lot over time, so we simplify the model by assuming that $N_S$ is constant. This leads to a two-dimensional model with state variables $p$ and $N_I$.

## two dimensional model -- fixed population density of susceptible patches

In the two-dimensional model, we have only two state variables: fraction of infected patches $p$ and average population density of infected patches $N$.

The equations are much less complicated than the three-dimensional one, and we get everything analytically in this version. We have now got the disease free equilibrium (DFE), the endemic equilibrium (EE) and the threshold condition for the stability of DFE.

[Detailed results](./fixedns)

## something about parameters

The reference in file [parameters.md](./parameters.md) gives some estimation of bubonic R0 (around 1.5). The paper estimated R0 by two methods:

-   using epidemic doubling time and the infectious period
-   using SIR model assuming human-flea-human

They didn't mention rats in their models or data, so it might not be appropriate to say R0 is the same value for rats.

If we assume human cases are caused by spillover from rats and are proportional to rat cases, we can use the same formula $R_0=1+(D/T_d)ln2$ and say that epidemic doubling-time $T_d=3.714$ is the same for humans and rats, but the infectious period $D$ should be the infectious periods of rats instead of humans (for example, 8 days), then we get $R_0=2.5$ (instead of 1.5)

This is close to the value $R_0=2.75$ we get by estimating from simplified ODEs considering spread between rats. I think this might be a more reasonable estimation although there is no reference giving this number...

# Update

Updates since last time:

-   find a scenario where reducing death probability $D$ makes originally stable DFE unstable
-   incorporated change of $N$ in [P1-R0 figure](./fixedns/P1atDFE.R) (looks similar to figure 5)
-   the scenario where decreasing $D$ destabilizes DFE with reasonable parameters It is a very narrow range for $R_0$ and $\varepsilon$ shown in [figure](./metapop/fixedns/matter-range.png) with other parameters fixed $c=0.5,S=300,B=10,r=0.5$, which also seems unlikely...Maybe test other parameters to see if this range can become wider?

Some challenges now:

-   stability of EE
-   estimation of other parameters, such as $B$,$c$,$r$,$S$ (Too many parameter to visualize...)

Current ideas:

-   our original idea about decreasing $R_0$ might still work?
-   try some different models that enables infected patches to infect other patches before they burn out. (For example, make transmission for next season proportional to the size of the epidemic in a patch, or burning-out patches always survive current iteration)
