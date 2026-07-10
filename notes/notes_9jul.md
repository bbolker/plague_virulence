## brain dump on metapopulation dynamics

<!-- LaTeX formatting on GitHub is ugly, need various workarounds that I probably won't bother with. See e.g. https://nschloe.github.io/2022/06/27/math-on-github-follow-up.html . Back-ticks inside $ might help, i.e. $`...`$. For non-whitespace-adjacent stuff, the hacky workaround is $\text{higher-}R_0$ instead of higher-$R_0$ ... -->

Figured out a lot of things (I think) talking with Yuyang today. The take-home messages, I think, are:

* our current approach to simulating/testing evolutionary invasion is probably wrong
* matching our intuition (or my intuition, anyway) for how a low (within-patch) $R_0$ strain can invade a higher-$ R_0$ resident will probably need a scenario where:
   * patch dynamics are asynchronous
   * recurrent exinctions are frequent
   
We may need add seasonality to make this work (we were thinking about seasonality anyway, to increase realism, but it may also be required to get the kind of dynamics we have in mind)

## dynamical scenario

My mental model for how a low-$R_0$ mutant strain can invade a high-$R_0$ wild-type/resident strain in a metapopulation context is that there are enough empty patches (where both strains have gone extinct) that the mutant often finds itself colonizing empty patches alone. Its low $R_0$ allows it to persist (avoid burnout), where a wild-type colonizer would burn out. (If the two strains co-colonize then the mutant will be caught in the wild-type-induced burnout ...)


## what's wrong with our current invasion sims

We simulate a monoculture of the wild-type for a long time (`strain2_delay` disease generations, typically 100: see e.g. `odin/euler_twostrain_singlepatchintro_examples.R`), then we introduce a small number of the mutant [in the most recent runs, 10 individuals in one randomly chosen patch]. Our default initial conditions for the monoculture/transient part of the sim are to start the host population size at the patch carrying capacity and pick a Poisson($\lambda=10$) deviate for each patch (this is an SID model, so $S(0)=K-I(0)$) (see `?plagueMetapop::discrete_run`). In other words, we are simulating an *initial invasion* of the wild-type.

This is a problem because, for large $R_0$/small $K$ parameters, with logistic host dynamics, the entire metapopulation is likely to burn out (intermediate connectivity [$\alpha$] *might* help, via "rescue effect", but low or high $\alpha$ will both be bad). When the mutant arrives many generations later, it will have no competition, so will be able to persist as long as its $R_0$ is large enough not to fizzle and small enough not to burn out. **This doesn't match my concept of the historical scenario**: we don't think that wild-type plague disappeared entirely before the Pla-depleted mutant arrived on the scene ...

Indeed, a high-$R_0$ mutant actually benefits from having a wild-type strain at equilibrium before it invades. The wild-type strain reduces the host population to its equilibrium ($K/R_{0,1}$), which makes the burnout probability negligible. If the mutant $R_0$ is higher than the wild-type, it can then easily invade occupied patches. This is a sort of reverse founder effect (!!)

My first thought was "oh, we should start the wild-type at $I^*$ rather than at very low prevalence". This might not work either; in this case, if the metapopulation persists, then the wild-type will probably achieve close to equilibrium density in every patch through repeated colonization.

I think the only way in which our dynamical scenario above can actually work is in a *fadeout regime*: the resident strain easily persists at the metapopulation scale, but frequently goes locally extinct in patches. We may have to read/think more about the classical *critical community size* [CCS] literature, which addresses exactly this phenomenon -- not the probability of burnout in a single patch on initial introduction, but the probability of 'fadeout' (not necessarily on initial invasion; I don't know the answer offhand, but we could look at historical England and Wales measles data to see the distribution of inter-fadeout durations [the 'seven cities' data on my web site might be sufficient for this, or we could dig up Bryan Grenfell's larger UK data set from wherever it's living on GitHub ... ??]). 

We could do single-strain experiments varying $K$, $R_0$, $\alpha$, and the distribution of initial conditions to see if there is a realistic, not tiny, parameter regime that leads to fadeout-type dynamics. Estimating CCS as a function of these parameters is really hard; I don't know if there are good, fast heuristic or more rigorous approximations in the literature that would let us figure this out without using brute force. For what it's worth I'm pretty sure that all of the existing CCS literature assumes linear demography ($dN/dt = r(1-N/K)$ rather than logistic ($rN(1-N/K)$), which could make a big difference ...) 

A *distribution* of patch sizes might also be important -- I think an implicit assumption in the current CCS literature is that there is some patch in the metapopulation above the CCS, from which faded-out patches can be recolonized. This was certainly the case in the UK in the periods people have considered -- several cities were larger than 250K. I don't know how much anyone has explored pure rescue-effect metapopulations, i.e. metapopulations where all patches are below the CCS and persistence is maintained by colonization between patches with asynchronous dynamics. This is certainly speculated about when people wonder how measles ever established in the first place, but I don't know if there is anything more than speculation.

## recurrent epidemics/fadeouts

As mentioned above, we may need to find a regime with asynchronous fadeouts (i.e. spontaneous local extinctions, not necessarily after a 'virgin soil' epidemic [i.e. one starting from low prevalence, $S(0) \approx K$). The endemic equilibrium prevalence in a population of size $K$, with scaled disease generation time $\epsilon$ (i.e. relative to replenishment rate/host generation time) is $I^* = K \epsilon(1-1/R_0)$; stochastic fluctuations are of order $\sqrt{I^{*}}$ (or, the coefficient of variance is of order $1/\sqrt{I^{*}}$). If the CV is much smaller than 1/2, then we're unlikely ever to get spontaneous stochastic extinction from an endemic equilibrium.

However, two phenomena -- *Bartlett cycles* (interaction between stochastic fluctuations and the damped oscillations of the deterministic model) and *seasonal forcing* can both maintain a much higher degree of variation (and hence higher extinction probability) than the endemic equilibrium would suggest (e.g. see pp. 24 and following of https://bbolker.github.io/math4mb/notes/epi1.pdf ).

We had previously thought that we might need to worry about seasonality and synchrony/asynchrony: @krauerInfluenceTemperatureSeasonality2021 (see `../virulence.bib`) show that European plague epidemics during the second pandemic typically occurred Aug-Nov (and @bacaerModelKermackMcKendrick2012 points out that the 1906 Bombay plague epidemic was probably seasonal). Seasonality is two-edged here; it increases variation (making local extinction more likely) but also increases synchrony (making recolonization/rescue effect less likely) ...  It would be easy enough to add (e.g.) sinusoidal variation in transmission rate (or host demography, or ...), at the cost of one more parameter (amplitude), but we should think before we do it. We have also considered some kind of hybrid model where we simulate until the end of an epidemic "season" and then jump to the beginning of the next -- but we would need to decide what assumptions to make about the dynamics of that interseason jump ...

(Biological note: what if recurrent plague epidemics represent recolonization events from some disease focus [i.e. a persistent reservoir] rather than metapopulation dynamics? Then all of this is out the window ...)

