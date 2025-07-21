---
title: "modeling abstract"
date: today
bibliography: "virulence.bib"
author: David Earn, Ben Bolker, Jonathan Dushoff
format:
  html:
    embed-resources: true
---

$$
\newcommand{\rzero}{{\cal R}_0}
$$

## Introduction

We (DJDE and BMB, with some input from JD) are exploring possible dynamical/evolutionary explanations for the replacement of wild-type by PLA-depleted strains of *Y. pestis* in the late stages ($\approx$ 100 years from onset) of multiple plague pandemics. The main epidemiologically relevant phenotypic differences between wild-type and PLA-depleted strains are:

* case mortality of infected hosts (rats) is lower
* infected hosts that do die take longer to die

Because transmission occurs only (mostly?) after the death of an infected rat, when its fleas leave to find other hosts, this means that PLA-depleted strains will have lower $\rzero$ *and* longer generation time $G$. While there is a huge literature on the evolution of virulence [e.g. @Cress2016] including well-known cases of evolution toward lower virulence (e.g. myxomatosis, syphilis), (nearly?) all theoretical explanations of declining virulence assume some trade-off between transmission rate and infectious period, possibly mediated by ecological conditions [@Aliz+2009a]. That is, the models depend on the lower-virulence strain having *higher* $\rzero$ under some conditions, due to their longer infectious period (part of $G$). In this case, because transmission depends on host mortality, $\rzero$ is (unusually) *independent* of $G$ and strictly proportional to virulence (here equated to case mortality). This rules out most of the standard theory, as evoked by @Lens1988:

> The dramatic declines of the human population in Europe during the great plague epidemics of past centuries were presumably accompanied by comparable declines in the population of susceptible rodents. Not only might these epidemics have been triggered by the appearance of hypervirulent strains of *Y. pestis* [...]  but the declining populations of susceptible hosts may in turn have favoured less virulent strains.

@LensMay1994 published one of the first models of the virulence-transmission trade-off. However, even though plague is the first example they cite in their introduction, they assumed (as is typical of the theoretical literature) that (1) hosts transmit disease throughout their infectious periods, and (2) transmission *rate* is an accelerating function of the *rate* of host mortality (their definition of 'virulence'). Thus their model formulation doesn't match the natural history of plague.

What other mechanisms could explain the emergence of PLA-depleted strains late in successive pandemics?

## "Fast-slow" mechanisms

All other things being equal, pathogen strains with shorter $G$ have a competitive advantage during the growth phase of an epidemic; equally clear from the theory but less widely appreciated[^1], strains with *longer* $G$ will have an advantage during the decline phase. Although the total opportunity for transmission over the course of an epidemic (i.e., the integral of prevalence over time) is larger for the faster strain than the slower strain, the slower strain will dominate the pathogen population later in the epidemic. If new outbreaks are initiated by sampling from late-stage populations, there *might* be a gradual shift toward slower (lower-$G$) strains. If the between-epidemic selection for slower strains is strong enough, it could potentially outweigh a disadvantage in $\rzero$.

This mechanism is fairly easy to model; we need a two-strain SIR model and some decision about how to model the initiation of new epidemics. The questions are:

* what is the biology behind the initiation of new epidemics, and what is the best way to translate that into part of a dynamical model? (We have some ideas.)
* if this mechanism is so prevalent, why don't we see evolution of slow strains all the time? i.e., what might be special about plague relative to all the other pathogens that have outbreak-type dynamics where slower strains don't seem to emerge over time?
* what would account for the long delay (hundreds of times a typical epidemic time scale) in the emergence of the PLA-depleted phenotype?

## Metapopulation/group-selective mechanisms

In spatial predator-prey or host-pathogen systems with patchy populations, phenotypes that underexploit their resource (aka "prudent predators") can evolve if spatial mixing rates are low enough that less prudent (or higher-$\rzero$/more virulent) strains have little opportunity to get into new patches and exploit their unused resources (hosts) [@bootsSmall1999; @claessenEvolution1995; @keelingReinterpreting2000b; @maySuperinfection1994a; @Pels+2002].

Some version of this phenomenon could apply to plague: if rat populations became sparser and patchier over time (due to some combination of disease mortality and environmental change), then the system could switch at some point to an ecological regime where the slower/lower-$\rzero$ strain would dominate.

Prudence-selecting mechanisms may be either *deterministic* (slower strains win by maintaining higher local host population densities, giving them an advantage in colonizing new patches/initiating new epidemics) or *stochastic* (slower strains win due to a lower probability of patch extinction). We are especially interested in the stochastic case because a recent paper of ours [@Pars+2024] derives an accurate approximation for the probability that a pathogen will "burn out" in an epidemic trough, given its $\rzero$ and $G$. In particular, lowering $\rzero$ and increasing $G$ could both increase the probability of persistence. It *might* be (??) relatively easy for us to make a simple patch-level model of extinction and colonization.[^2] (@King+2009 present a related mechanism for the evolution of acuteness in *Bordetella*.)

* the main thing that has held us back from exploring these hypotheses so far is that spatial models are more complicated, harder to program, and computationally slower, than the non-spatial models we need to explore fast-slow mechanisms ... [^3]
* could we take a shortcut by building a patch-level model of strain competition based on the results of @Pars+2024?
* as before: what explains the long time to emergence of PLA- strains? How might we model long-term changes in the rat population[^4]? How would we incorporate such trends in a model in a plausible way?

---

[^1]: this is "obvious"/easy to show, but we don't know of any published paper that states it

[^2]: we just thought of this mechanism today (10 April), so we haven't gone into the details - it might be harder than we think ...

[^3]: @KeelGill2000 is the canonical model of plague metapopulations, but (a) it would be a very large effort to reconstruct/re-program it and (b) there are aspects of it that we find fishy, although TBH I don't remember what they were - their model of transmission *is* reasonably realistic, with an input rate of infectious, searching ("free") fleas proportional to the death rate of infected rats ...

[^4]:  Is there any *data* on rat population/proxies [coalescent models???] @Davis1986 mostly covers historical/archaeological evidence for the *presence* of *Rattus rattus*, not anything (on a quick skim) close to a temporal trend in density ...
