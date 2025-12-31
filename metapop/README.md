---
title: "overview of plague metapopulation model"
bibliography: "../virulence.bib"
---

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

We want to do a simulation to see the more realistic situation.[simulation](./stochastic)

## something about parameters

The reference in file [parameters.md](./parameters/parameters.md) gives some estimation of bubonic plague R0 (around 1.5). The paper estimated R0 by two methods:

-   using epidemic doubling time and the infectious period
-   using SIR model assuming human-flea-human

They didn't mention rats in their models or data, so it might not be appropriate to say R0 is the same value for rats.

If we assume human cases are caused by spillover from rats and are proportional to rat cases, we can use the same formula $R_0=1+(D/T_d)ln2$ and say that epidemic doubling-time $T_d=3.714$ is the same for humans and rats, but the infectious period $D$ should be the infectious periods of rats instead of humans (for example, 8 days), then we get $R_0=2.5$ (instead of 1.5)

This is close to the value $R_0=2.75$ we get by estimating from simplified ODEs considering spread between rats. I think this might be a more reasonable estimation although there is no reference giving this number...

## Thoughts about evolutionary models

The existence of a peak in metapopulation occupancy for intermediate values of R0 is exciting, and suggestive that intermediate values of R0 might be an evolutionary optimum, but is far from proving it. (Maximizing (quasi)equilibrium density or metapopulation occupancy is a sufficient condition for evolutionary optimality in **some** simple ecological/epidemiological models, but you definitely can't take this property for granted.)

One paper that briefly discusses eco-evolutionary models: Abrams 2001. “Modelling the Adaptive Dynamics of Traits Involved in Inter- and Intraspecific Interactions: An Assessment of Three Methods.” Ecology Letters 4 (2): 166–75. https://doi.org/10.1046/j.1461-0248.2001.00199.x 

This is particularly tricky in our case because we make a lot of assumptions in collapsing the full system to a discrete-time occupancy model. Some questions:

* in a mixed epidemic (starting from initial conditions $\{S, I_1(0), I_2(0), R\}$), where there are two parasite strains with (slightly) different $R_0$ values, what are (1) final sizes and (2) probabilities of burnout? (fizzle is easier, as these will be independent for each strain)
* how do we model colonization of multiple strains? Since initial conditions will be important, the suggestion is that we pick *Poisson-distributed* numbers of colonists per strain (based on the total population size in infected patches and the colonization-proportion parameter) rather than a hazard-based
* we would also need to pick some rule for starting conditions (number of infected individuals at the start of each season) for patches where infection persists (rather than assuming that infection resets to 1 starting individual every time)

To construct a *pairwise invasibility plot*, we would need to measure the initial growth rate of a small number of strain 2 invading a monoculture of strain 1. More specifically, if the equilibrium of the strain 1 monoculture is $\{S_1^*, I_1^*, 0, R_1^*\}$ [where the zero denotes that strain 2 is absent], then we want to know the short-term growth rate of $I_2$ (where $R_{0,2} = R_{0,1} + \delta$ from a starting condition of $\{S_1^*, I_1^*, I_2^i, R_1^*\}$ where typically we assume $\delta \ll 1$ (small mutations) and $I_2^i \ll 1$ (small invading population). If we are doing a stochastic simulation we also need $\delta$, $I_2^i$ to be large _enough_ so that we can detect a signal and the chance of fizzling isn't too high.

Doing this in a metapopulation context (1) means there are more possible ways to set up initial conditions (do we assume the invaders all start in a single patch, or are multiple patches invaded?) and (2) makes brute-force solutions harder (larger state space).

We also have to make some assumptions about the two-strain infection process. The simplest case is that there is no coinfection and perfect cross-immunity and that the generation times are equal i.e. (in a single patch/well-mixed population, with $\gamma = 1$)

$$
\begin{split}
\dot S & = \sum_j -R_{0,j} S I_j \\
\dot I_j & = R_{0,j} S I_j - I_j \\
\dot R & = \sum I_j
\end{split}
$$

If we can't figure out a good way to get analytical solutions and/or approximations for the final-size and burnout probabilities in a two-strain model, we might be able to use `odin` to quickly spin up a discrete-time stochastic *tau-leaping* model for a single population: see [here](https://mrc-ide.github.io/odin/articles/discrete.html#stochastic-processes), [here](https://mrc-ide.github.io/odin/articles/discrete.html#stochastic-sir-model)

* [example (complicated) multi-strain odin model](https://github.com/abhisheksena/covid_multi_strain/blob/main/inst/odin/covid_multi_strain.R)

```{r}
odin::odin("metapop/odin_twostrain0.R")
```

We might also be able to use that framework for the full stochastic metapop model (discrete-time measured in *epidemic periods* rather than some small Δt): it might well run faster than even a vectorized R-based simulator ...

References (see `../virulence.bib`): @brannstromHitchhikersGuideAdaptive2013a; @diekmannBeginnersGuideAdaptive2002; @abramsModelling2001

### To do (odin)

* two-strain model is more or less working
    * could be extended to arrays to run many replicates simultaneously/efficiently (see https://mrc-ide.github.io/odin/articles/discrete.html)
	* is there a way to do early stopping/stop when extinct?
	* is there a way to store/retrieve the compiled model? (Does it use some kind of caching/make rule anyway?)
* how do we estimate burnout probs from stoch anyway? (How was this done for the burnout paper?) Need vital dynamics ...
* most efficient form of solving the invasion problem would simulate a one-strain stoch model to quasi-equilibrium, then use the final state as the starting condition (+ a few invaders) to look at invasion/growth rate [could also simulate the two-strain model with strain 2 absent, which would be a little bit inefficient]
* shouldn't need a stochastic model for final size(s) of two-strain epidemic, but: how do we get values if there's no closed form solution? A big lookup table? Gaussian process emulator?  (Is this actually solvable directly using some of the same methods as for the regular SIR, i.e. can we solve for the phase-space trajectory without knowing $t$?)
* more generally, how big a lookup table (or multi-dimensional emulator) do we need [what dimensions?] for final size, burnout ... ? Am I overcomplicating this?

* true brute force would run the whole metapop model in (approx) continuous time, i.e. not breaking it up into epidemic episodes ...

# Change log

Update since October 2025:

- Quasi-equilibrium
- Simulations for stochastic model

Update since 11 Sep 2025:

-   Equations for stochastic simulation
-   Hazard model and colonization in infection probability

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
