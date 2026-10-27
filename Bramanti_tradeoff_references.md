# Do Bramanti et al.'s refs 26–28 support the virulence–transmission trade-off?

*Notes compiled 2026-10-06, from the full texts of all three cited papers.*

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
that sense. None measures a transmission rate or a reproductive number, and none
relates transmission to virulence within a host species. The claim that "this
trade-off hypothesis was previously demonstrated for *Y. pestis* (26, 27)" is not
supported by the cited work.

One clause of the passage does hold up. Ref 28 measures host survival time and
shows that losing Pla prolongs it, so "reduce virulence and enhance the time of
survival of the host" is supported for *pla* in a mouse pneumonic model. What
fails is the step that follows — "and, consequently, of the pathogen" — for
which the same paper supplies evidence to the contrary.

## Ref 26 — Hinnebusch et al. 2017

A review of *Y. pestis* biofilm biology in the flea gut. Mammalian virulence
enters only as a bacteremia threshold that transmission is gated on, and on that
point the review argues roughly the **opposite** of the use Bramanti et al. make
of it.

The decisive passage states that the transmission window is bounded by host
*death*, not extended by host survival:

> "Transmission by this mode is rare unless the infectious blood meal contains
> at least 10⁸ *Y. pestis*/mL, and reported early-phase transmission efficiency
> values are based on blood bacteremia levels of ∼10⁹ *Y. pestis*/mL or higher.
> With the notable exception of mice, terminal bacteremias may not routinely
> reach 10⁹/mL in most mammals, or may occur only shortly before death. **Thus,
> if it occurs at all, there is a very brief interval between the early-phase
> threshold bacteremia level and death.**"

That is a near-threshold β(α) relationship: lowering virulence moves a host
*down* through the threshold and toward β = 0.

Sublethal infections are transmission dead ends:

> "Intermittent challenges from just a few fleas at a time would frequently in
> effect vaccinate animals and remove them from the susceptible population,
> rather than cause the septicemic plague necessary to complete the transmission
> cycle and drive epizootic spread."

The authors conclude that late-stage, blockage-dependent transmission — the mode
depending on terminal bacteremia — "provides the foundation for ecologically
stable plague transmission cycles."

The review contains no discussion of chronic or subclinical infection in
resistant rodents. Its only use of "chronic" refers to chronic proventricular
biofilm infection *in fleas*, and every "persistent" refers to colonization of
the flea gut. There is no mammalian duration-of-infectiousness argument anywhere
in it.

### The one real trade-off here is in the vector, not the host

The paper does contain a duration-versus-intensity trade-off, but it concerns
*flea* survival:

> "Bacot stressed the importance of transmission by partially blocked fleas and
> considered them to be more efficient transmitters than completely blocked
> fleas. They also have a shorter EIP and live longer."

Complete blockage starves the flea, which then "will spend the last few days of
its life trying to obtain a blood meal, making repeated probing attempts, each
one potentially resulting in transmission." Higher biofilm "virulence" toward
the vector thus buys intense but brief infectiousness, while partial blockage
buys a longer infectious period.

This is a genuine virulence–transmission trade-off, structurally the one
Bramanti et al. want — but its host is the *flea*, and the trait under selection
is biofilm production, not *pla*. It cannot be transferred to mammalian virulence
or to *pla* decay without an argument the authors do not make.

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

## Ref 28 — Lathem et al. 2007

This is the only one of the three that measures anything resembling a component
of the trade-off, and it measures the duration side alone.

**Losing Pla prolongs host survival.** Wild-type CO92 killed mice synchronously,
whereas "only 50% of the mice infected with the ∆pla strain developed terminal
plague after 7 days"; the authors state that "the lack of Pla substantially
delayed the time to death." With a tetracycline-inducible construct, mean time
to death was 3.1 days in the *pla*-induced state versus 5.1 days when repressed.
The abstract frames this therapeutically: "inhibition of Pla expression prolonged
the survival of animals with the disease."

So the narrow clause Bramanti et al. attribute to this paper — reduced virulence
lengthening host survival — is supported, for *pla*, in a mouse intranasal model.
That is the mutation at issue, so the citation is apt as far as it goes.

**But the survival is bought by a collapse in pulmonary bacterial load:**

- 100- to 1000-fold fewer bacteria recovered from ∆pla lungs at 24 h; over the
  next two days ∆pla numbers "did not substantially change, whereas wild-type
  bacteria increased by almost 6 logs."
- ∆pla lungs "showed no change in weight, even after 7 days," so the mice died
  of systemic infection rather than pneumonia.
- Inflammation aborted and lung repair was activated (PCNA-positive cells).

For a respiratory pathogen, infectious output is bacterial load in the airway.
The ∆pla host lives longer while being, as far as the lung is concerned, close to
non-infectious. Longer duration times near-zero transmission rate is not a
compensating gain — it is the branch of the trade-off where R₀ falls.

Absent from the paper: any transmission experiment, co-infection, or measurement
of onward infection, and anything about flea-borne transmission, which is the
route relevant to the Black Death reservoir argument.

## Why the citation chain doesn't do the work

1. **No one has measured a *Y. pestis* transmission–virulence trade-off curve.**
   There is no estimated β(α) relationship, so there is no basis for claiming an
   interior R₀ optimum, let alone that *pla*-decay strains moved toward it.
2. **The mechanism is backwards for flea-borne transmission.** Transmission is
   gated on a bacteremia threshold reached only near death, so β is a steeply
   increasing (near-threshold) function of α. With a hard threshold, lowering
   virulence can drive β toward 0 — reducing virulence cannot increase R₀. Ref 26
   states this outright: "there is a very brief interval between the early-phase
   threshold bacteremia level and death."
3. **Longer survival is demonstrated; longer *infectiousness* is not.** Ref 28
   shows prolonged host survival under *pla* loss, but with 100–1000× fewer
   bacteria in the lung — duration up, infectious output down. No paper here
   shows the product increasing. The only other duration-like quantity, in Cui
   et al., is persistence of the bacterium in the *flea foregut*, a different
   compartment entirely.
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

- Hinnebusch, Jarrett & Bland 2017, *Annu. Rev. Microbiol.* 71:215–232 —
  full text, doi:10.1146/annurev-micro-090816-093521
- Cui et al. 2020, *Nat. Commun.* 11:281 — full text,
  doi:10.1038/s41467-019-14099-w
  ([PMC6962365](https://pmc.ncbi.nlm.nih.gov/articles/PMC6962365/))
- Lathem, Price, Miller & Goldman 2007, *Science* 315:509–513 — full text,
  doi:10.1126/science.1137195
- Bland et al. 2020, *PLoS Pathog.* 16:e1009092 — transmission efficiency and
  plague progression by the two flea mechanisms, doi:10.1371/journal.ppat.1009092
- Plasminogen activator Pla — <https://en.wikipedia.org/wiki/Plasminogen_activator_Pla>
