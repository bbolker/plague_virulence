# overview of plague metapopulation model

## three dimensional model

We have a discrete-time model keeping track of three state variables: fraction of infected patches $p$, average population density of infected patches $N_I$ and average population density of susceptible patches $N_S$.

There are several events in each iteration that can change the state of the system: movement of hosts and infection of susceptible patches; death of infected hosts and burnout of infected patches; and growth of hosts in both susceptible and infected patches.

Detailed description of the model in the [documents](./3d)

Simulation indicates that population density of susceptible patches doesn't change a lot over time, so we simplify the model by assuming that $N_S$ is constant. This leads to a two-dimensional model with state variables $p$ and $N_I$.

## two dimensional model -- fixed population density of susceptible patches

In the two-dimensional model, we have only two state variables: fraction of infected patches $p$ and average population density of infected patches $N$.

The equations are much less complicated than the three-dimensional one, and we get everything analytically in this version. We have now got the disease free equilibrium (DFE), the endemic equilibrium (EE) and the threshold condition for the stability of DFE.

condition for persistence is $L_0>1$ where $L_0$ is the landscape reproduction number: $L_0=(1+B)P_1|_{N=\frac{cS}{1-(1+r)(1-zD)(1-c)}}$

[Detailed results](./fixedns)

## transmission proportional to size of within-patch epidemic

In previous models, burned-out patches are set to be susceptible in the next season and have no chance to spread disease. But it is very likely that a patch with a large epidemic can infect other patches before it burns out.

We can model this by assuming that the transmission for next season is proportional to the size of the epidemic in a patch.

We have $L_0=P_1 + azP_{nf}$, where $a$ is a scaling constant, $P_{nf}=1-\frac{1}{R_0}$ is the “non-fizzle” probability and $z$ is the final size. $P_1$ is evaluated at $N=\frac{cS}{1-(1+r)(1-c)(1-zD)}$

[Details in documents](./within_season_transmission)

$P_1$ is non-monotonous with $R_0$, while $zP_{nf}$ is increasing with $R_0$, so there can be some region where $L_0$ can decrease with $R_0$, and we can find such scenario that decreasing $R_0$ destabilizes DFE and makes disease persist.

$L_0$ is also increasing with $D$, so a lower death probability can also destabilize the DFE.

## keeping track of average population density of all patches

In previous version $N_{DFE}$ is a function of many parameters (c,r,S,R0,D) and moving probability c might be difficult to estimate. And if a patch gets infected and then burns out, its population size should be $(1-zD)N$, but our assumption makes it suddenly go up immediately to S.

To solve these we consider modeling average population density of all patches instead of infected patches only, so moving probability c is no longer necessary, and it’s changing continuously, which makes things easier and more reasonable.

In this version we get $L_0=P_1|_{N=K} + B$ and persistence condition is $L_0>1$

## Stochastic simulation

In previous models we were all dealing with average population size, but this might not be very reasonable and leads to a jump discontinuity in population size when a patch gets infected.

We want to do a simulation to see the more realistic situation.

## something about parameters

The reference in file [parameters.md](./parameters/parameters.md) gives some estimation of bubonic plague R0 (around 1.5). The paper estimated R0 by two methods:

-   using epidemic doubling time and the infectious period
-   using SIR model assuming human-flea-human

They didn't mention rats in their models or data, so it might not be appropriate to say R0 is the same value for rats.

If we assume human cases are caused by spillover from rats and are proportional to rat cases, we can use the same formula $R_0=1+(D/T_d)ln2$ and say that epidemic doubling-time $T_d=3.714$ is the same for humans and rats, but the infectious period $D$ should be the infectious periods of rats instead of humans (for example, 8 days), then we get $R_0=2.5$ (instead of 1.5)

This is close to the value $R_0=2.75$ we get by estimating from simplified ODEs considering spread between rats. I think this might be a more reasonable estimation although there is no reference giving this number...

# Update

Update since 11 Sep 2025:

-   Equations for stochastic simulation

Update since 21 Aug 2025:

-   A version keeping track of average population density of all patches.

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

# Points and assumptions we might need to return to

-   We are assuming there is always only 1 initially infected individual
-   We are dealing with "average" population density in patches. This means that we are neglecting variation in population density within each group.
-   We are assuming that the population density of susceptible patches is constant.
-   We are taking use of burnout probability from a SIR model, while our model of bubonic plague is slightly different (spread by fleas/different birth-and-death process) so some parameters might have a different meaning
-   We are assuming that $R_0$ is independent of host population size $N$ (frequency-dependent transmission). We might consider $R_0$ scaling linearly with $N$ (density-dependent transmission) or something between.
-   We are assuming that a burned-out patch is susceptible in the next season as if it was never infected, and a patch that persists the first wave can still go extinct in the next wave like a newly infected patch.
