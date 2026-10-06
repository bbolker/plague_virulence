# Do Bramanti et al.'s refs 26–28 support the virulence–transmission trade-off?

*Notes compiled 2026-10-06 (Claude Code session).*

## The question posed

In Bramanti et al., "Assessing the origins of the European Plagues following the
Black Death: A synthesis of genomic, historical, and ecological information", the
authors write:

> One of the optimized survival strategies for an emerging pathogen is to balance
> its virulence to the main host with its transmission strategy. This trade-off
> hypothesis was previously demonstrated for *Y. pestis* (26, 27). This mechanism
> would allow the bacterium to reduce virulence and enhance the time of survival
> of the host and, consequently, of the pathogen (28). After experiencing the
> Black Death and successive waves, the *pla* decay strains might have attempted
> to acquire a fitness advantage, reducing their virulence.

The cited references are:

26. B. J. Hinnebusch, C. O. Jarrett, D. M. Bland, "Fleaing" the plague:
    Adaptations of *Yersinia pestis* to its insect vector that lead to
    transmission. *Annu. Rev. Microbiol.* **71**, 215–232 (2017).
27. Y. Cui et al., Evolutionary selection of biofilm-mediated extended phenotypes
    in *Yersinia pestis* in response to a fluctuating environment. *Nat. Commun.*
    **11**, 281 (2020).
28. W. W. Lathem, P. A. Price, V. L. Miller, W. E. Goldman, A plasminogen-activating
    protease specifically controls the development of primary pneumonic plague.
    *Science* **315**, 509–513 (2007).

**Question:** what information in those papers directly speaks to the trade-off
hypothesis of the evolution of virulence — i.e. that prolonging the infectious
period can increase the reproductive number even if it comes with a reduced
transmission rate?

## Headline conclusion

None of the three papers tests, or even addresses, the trade-off hypothesis in
that sense. None of them measures or models an infectious period, a transmission
rate, or a reproductive number. The claim that "this trade-off hypothesis was
previously demonstrated for *Y. pestis* (26, 27)" is not supported by the cited
work.

### Access status at time of writing

| Ref | Access | Basis for notes below |
|---|---|---|
| 27 Cui et al. 2020 | open ([PMC6962365](https://pmc.ncbi.nlm.nih.gov/articles/PMC6962365/)) | full text read |
| 26 Hinnebusch et al. 2017 | paywalled (Annual Reviews, no PMC deposit) | abstract + related lab literature only |
| 28 Lathem et al. 2007 | paywalled (*Science*) | abstract + secondary coverage only |

Refs 26 and 28 still want a full-text check — in particular whether Hinnebusch
discusses chronic/subclinical infection in resistant rodents, which is the one
place in that review where a duration-of-infectiousness argument could
plausibly live.

## Ref 27 — Cui et al. 2020

The paper does report a trade-off, but it is the wrong trade-off: it is between
**two vector transmission modes**, mediated by the amount of biofilm.

> "The two transmission modes are in a trade-off with each other. Increased
> levels of biofilm formation lead to a better ability of the bacterium to
> maintain itself in the foregut of the flea... However, these increased levels
> of biofilm formation decreases the efficiency of early-phase transmission."

That is blockage-dependent vs. early-phase transmission, both occurring in the
flea. The gene under selection is *rpoZ* (RNA polymerase ω-subunit); variants
show higher in-vitro biofilm production and are enriched during colder, drier
periods in the Guertu focus.

There is **no measurement or discussion of mammalian virulence, host survival
time, attenuation, or R₀** anywhere in the paper. Citing it as a demonstration
of the virulence–transmission trade-off conflates "a trade-off exists among
*Y. pestis* phenotypes" with "the virulence–transmission trade-off hypothesis".

## Ref 26 — Hinnebusch et al. 2017

From the abstract and the surrounding Hinnebusch-lab literature, this review
argues roughly the **opposite** of the use Bramanti et al. make of it. Flea
transmission requires host bacteremia of ~10⁸–10⁹ CFU/ml, essentially a terminal
condition, and transmission probability is *positively* correlated with that
lethal bacteremia level.

The companion experimental paper (Bland et al. 2020, *PLoS Pathog*,
[ppat.1009092](https://journals.plos.org/plospathogens/article?id=10.1371%2Fjournal.ppat.1009092))
makes the point sharply: early-phase transmission more often yields sublethal,
immunizing, **non-productive** infections — hosts that survive do not reach
transmissible bacteremia — whereas blocked-flea transmission yields terminal
disease and onward transmission.

Under that biology, reduced virulence costs transmission rather than buying extra
infectious days. This is a case where virulence and transmissibility are
positively coupled, which is the classic circumstance under which the trade-off
argument *fails*.

## Ref 28 — Lathem et al. 2007

Purely mechanistic and route-specific: Pla is required for *Y. pestis* to
establish fulminant primary pneumonic plague; without Pla expression,
inflammation aborts and lung repair is activated. Pla is less important for
dissemination in pneumonic than in bubonic plague.

No transmission experiment, no co-infection, no host-survival-vs-transmission
comparison. It supports "*pla* loss would impair pneumonic disease" — which, for
Bramanti et al.'s own argument, implies *reduced* onward pneumonic transmission,
not a compensating gain.

## Why the citation chain doesn't do the work

1. **No one has measured a *Y. pestis* transmission–virulence trade-off curve.**
   There is no estimated β(α) relationship, so there is no basis for claiming an
   interior R₀ optimum, let alone that *pla*-decay strains moved toward it.
2. **The mechanism is backwards for flea-borne transmission.** Transmission is
   gated on a bacteremia threshold reached only near death, so β is a steeply
   increasing (near-threshold) function of α. With a hard threshold, lowering
   virulence can drive β toward 0 — reducing virulence cannot increase R₀.
3. **"Prolonging the infectious period" has no demonstrated correlate in these
   papers.** The one duration-like quantity in Cui et al. is persistence of the
   bacterium in the *flea foregut*, not the mammalian infectious period — a
   different compartment entirely.
4. **Refs 26 and 27 concern the enzootic rodent–flea cycle**, whereas the
   Bramanti claim is about human epidemic waves, where (if pneumonic or
   louse-borne transmission matters) the relevant β–α relation is different
   again.

## A defensible reading of the same three papers

*pla* decay is expected to reduce pneumonic transmissibility (28), while
biofilm-regulator variation tunes the balance between two vector transmission
routes under climate fluctuation (27), in a system where host bacteremia high
enough to infect fleas is near-lethal (26). That is a story about
**route-switching and reservoir ecology**, not about attenuation for
host-survival benefit.

## Sources consulted

- Cui et al. 2020, *Nat. Commun.* 11:281 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC6962365/>
- Bland et al. 2020, *PLoS Pathog.* — transmission efficiency and plague progression
  by the two flea mechanisms —
  <https://journals.plos.org/plospathogens/article?id=10.1371%2Fjournal.ppat.1009092>
- Hinnebusch et al. 2017 (record) —
  <https://www.researchgate.net/publication/319613290_Fleaing_the_Plague_Adaptations_of_Yersinia_pestis_to_Its_Insect_Vector_That_Lead_to_Transmission>
- Lathem et al. 2007, summary coverage —
  <https://source.washu.edu/2007/01/disabling-key-protein-may-give-physicians-time-to-treat-pneumonic-plague/>
- Plasminogen activator Pla — <https://en.wikipedia.org/wiki/Plasminogen_activator_Pla>
