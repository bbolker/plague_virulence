---
title: "speculative modeling of plague virulence evolution"
bibliography: virulence.bib
---

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

