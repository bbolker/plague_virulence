# overview of plague metapopulation model

## three dimensional model

We have a discrete-time model keeping track of three state variables: fraction of infected patches $p$, average population density of infected patches $N_I$ and average population density of susceptible patches $N_S$.

There are several events in each iteration that can change the state of the system: movement of hosts and infection of susceptible patches; death of infected hosts and burnout of infected patches; and growth of hosts in both susceptible and infected patches.

Detailed description of the model in the [documents](./3d)

Simulation indicates that population density of susceptible patches doesn't change a lot over time, so we simplify the model by assuming that $N_S$ is constant. This leads to a two-dimensional model with state variables $p$ and $N_I$.

## two dimensional model -- fixed population density of susceptible patches

In the two-dimensional model, we have only two state variables: fraction of infected patches $p$ and average population density of infected patches $N$.

The equations are much less complicated than the three-dimensional one, and we get everything analytically in this version. We have now got the disease free equilibrium (DFE), the endemic equilibrium (EE) and the threshold condition for the stability of DFE.

[Detailed results](./fixedns)

## model 3 -- transmission proportional to size of within-patch epidemic

In previous models, burned-out patches are set to be susceptible in the next season and have no chance to spread disease. But it is very likely that a patch with a large epidemic can infect other patches before it burns out.

We can model this by assuming that the transmission for next season is proportional to the size of the epidemic in a patch. Instead of our current landscape reproductive number L0 = (1+B)p_1, we would have something like p_1 + p_nf\*Z, where p_nf=1-1/R0 is the “non-fizzle” probability and Z is the final size. We'll need to see whether this works out with reasonable parameters.

## model 4 -- burning-out patches always survive current iteration

## something about parameters

The reference in file [parameters.md](./parameters/parameters.md) gives some estimation of bubonic R0 (around 1.5). The paper estimated R0 by two methods:

-   using epidemic doubling time and the infectious period
-   using SIR model assuming human-flea-human

They didn't mention rats in their models or data, so it might not be appropriate to say R0 is the same value for rats.

If we assume human cases are caused by spillover from rats and are proportional to rat cases, we can use the same formula $R_0=1+(D/T_d)ln2$ and say that epidemic doubling-time $T_d=3.714$ is the same for humans and rats, but the infectious period $D$ should be the infectious periods of rats instead of humans (for example, 8 days), then we get $R_0=2.5$ (instead of 1.5)

This is close to the value $R_0=2.75$ we get by estimating from simplified ODEs considering spread between rats. I think this might be a more reasonable estimation although there is no reference giving this number...

# Update

Updates since 15 Aug 2025:

-   some ideas about estimating $R_0$ and $epsilon$
-   incorporated change of $N$ in [P1-R0 figure](./fixedns/P1atDFE.R) (looks similar to figure 5)
-   tried to find the scenario where decreasing $D$ destabilizes DFE with reasonable parameters. It is a very narrow range for $R_0$ and $\varepsilon$ shown in figures created by [script](./fixedns/more_matter-range.R) with other parameters fixed $c=0.5,S=300,B=10,r=0.5$, which also seems unlikely...

Updates since 7 Aug 2025:

-   Replace original infection probability $Bp(1-p)$ with hazard model $(1-e^{-Bp})(1-p)$
-   discussed the order of events
-   Next step: we should try to find a scenario where Pla-depletion strain (decreased R0 or D) is favored. The problem is whether reasonable parameters are in the range that our model can work.

What we've get before 7 Aug 2025 

-   Two-dimensional model and parameter range
-   condition for persistence: $L_0>1$ where $L_0=(1+B)P_1|_{N=\frac{cS}{1-(1+r)(1-zD)(1-c)}}$

Some challenges now:

-   stability of EE
-   estimation of other parameters, such as $B$,$c$,$r$,$S$ (Too many parameter to visualize...)

Current ideas:

-   our original idea about decreasing $R_0$ might still work?
-   Re-run the code for different c and r.
-   try some different models that enables infected patches to infect other patches before they burn out.

# Points and assumptions we might need to return to

-   We are dealing with "average" population density in patches. This means that we are neglecting variation in population density within each group.(This seems fine in two-dimensional model but really strange in three-dimensional one)
-   We are assuming that the population density of susceptible patches is constant.
-   We are taking use of burnout probability from a SIR model, while our model of bubonic plague is slightly different (spread by fleas/different birth-and-death process) so some parameters might have a different meaning
