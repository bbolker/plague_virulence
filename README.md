---
title: "speculative modeling of plague virulence evolution"
bibliography: virulence.bib
---

## overview

This repo is highly organic, reflecting evolution from a collection of ideas about how to *verbally* explain the evolution of attenuated virulence noted by Sidhu *et al* 2025, to a series of models trying to instantiate and explore those verbal models. A lot of this builds on Parsons *et al* 2025 on "burnout" probabilities.

* the documents in the head directory of the repo (`abstract2.qmd`, `main.qmd`, `notes_*`) are discussion of the verbal model; `virulence.bib` is a collection of useful references on evolution of virulence and math models of plague.
* `metapop` represents our attempts to develop analytical (Levins-style) and stochastic-simulation models of the metapopulation dynamics. For the single-species Levins-style model, we considered both 2-dimensional (proportion of patches infected, average patch population size) and 3-D (average pop size for both susceptible and infectious patches) model.

Stuff in `metapop/`, from newest to oldest:

* `odin`: mostly continuous-time simulations with [odin](https://mrc-ide.github.io/odin/) (and also comparable implementations with [macpan2](https://canmod.github.io/macpan2/)
* `plagueMetapop`: mini-R package implementing the functionality used in the `odin` runs
* `parameters`: discussion of estimating/guessing orders of magnitude for parameters
* `poisson`: a single-strain, stochastic model driven by rat colonization limitation (colonization uses Poisson deviates)
* `two_strain`: a two-strain version of the Poisson model (PIP computations)
* `csee_talk`: BMB talk for Canadian Society of Ecology and Evolution, May 2026
* `stochastic`: stochastic model with burnout 
* `two_strain_burnout_probability`: calculations of burnout probability in two-strain (coinfected) patches
* `3d`: Levins model, state variables {average host pop size in S patches; average host pop size in I patches; fraction of patches infected}
`within_season_transmission`: ?
`report_for_summer`: YZ's report
`fixedns`: : 2D model Levins model, state variables: {average host pop size; fraction of patches infected}
`averaging_all_patches`: ?
`outputs`

Parsons, Todd L., Benjamin M. Bolker, Jonathan Dushoff, and David J. D. Earn. 2024. “The Probability of Epidemic Burnout in the Stochastic SIR Model with Vital Dynamics.” Proceedings of the National Academy of Sciences 121 (5): e2313708120. https://doi.org/10.1073/pnas.2313708120.

Sidhu, Ravneet Kaur, Guillem Mas Fiol, Pierre Lê-Bury, et al. 2025. “Attenuation of Virulence in Yersinia Pestis across Three Plague Pandemics.” Science 388 (6750): eadt3880. https://doi.org/10.1126/science.adt3880.


## general thoughts/brain dump

The observed pattern is that an avirulent strain appears to emerge
repeatedly during plague pandemics (reduced copy number of a plasmid
gene coding for plasminogen). Why?

@LensMay1994, @KeelGill2000

From @Lens1988:

> The dramatic declines of the human population in Europe during the great plague epidemics of past centuries were presumably accompanied by comparable declines in the population of susceptible rodents. Not only might these epidemics have been triggered by the appearance of hypervirulent strains of *Y. pestis*, as Rosqvist, Skurnik and Wolf-Watz hypothesize, but the declining populations of susceptible hosts may in turn have favoured less virulent strains.

Lenski and May paper uses logistic growth in host, with possible reproduction by infected hosts at a lower rate than susceptibles. They find the stable equilibrium of the SIR + vital dynamics model. They use a quadratic virulence-transmission tradeoff, and (eventually) find the eco-evolutionary ESS (find the equilibrium density of susceptible hosts for a specified transmission rate $b$, then find the minimum host population size by solving $dH^*/db = 0$ for $b$.

LM94 discuss lots of variations (density-dep vs independent host growth, density vs freq-dependent transmission, effects of migration, immunity, recovery, timescales, etc.).

Points of interest:

* according to HP the PLA gene does *not* decrease case mortality, it just delays it. This departs from the usual "virulence as rate of mortality" model. In principle (???) this should just slow down epidemics, not make them smaller? Also, in the case of bubonic plague, transmission should be strongly tied to host mortality (as fleas don't leave for a new host until their original host dies). (I'm still not sure that I completely understand the effects of PLA copy number on epidemiological parameters ... is transmission really not reduced?)
* How does the previous paragraph interact with differential selection in epidemic and endemic phases? What is the relationship between $r$ and $R_0$ in this case?
* Given basic demographic and epidemiological parameters, what would the expected time scale of a LM94-type scenario? That is, how fast would we expect host (rat) densities to decline from their pre-epidemic level, and how soon would we expect to hit a threshold where an avirulent (or differently virulent) strain would become competitively dominant?

## further thoughts (12 March meeting)

* population structure: what is the minimal model for working out reasonable hypotheses about the effects of population density and structure on invasibility/ESS? Relevant papers are @bootsSmall1999 (pair approximations and stoch-sim: mixture of nearest-neighbour + global dispersal); @claessenEvolution1995 (also NN?); @keelingReinterpreting2000b (patch-level moment approximation results); @maySuperinfection1994a (patch-level competition-colonization tradeoff à la Tilman). Do we need stoch sims to account for local fadeout? (*Not* interested in evolution of resistance, dispersal, local adaptation, ...)
* DJDE's verbal argument:
   * populaton size of rats depleted by epidemics
   * when a given rat dies, fewer hosts
   * inverse density-dependence -- more fleas per host
   * more fleas -> increased prob of infection & mortality
   * less secondary transmission, but tertiary inf is higher
* tangentially: all else equal, a slower generation time is advantageous/selected during the *decreasing* phase of an epidemic. Does this interact with the other mechanisms?
