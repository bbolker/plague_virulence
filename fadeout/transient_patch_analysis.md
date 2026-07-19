# Transient Patch Contribution Analysis

## Purpose

The analytical patch-occupancy approximation replaces the metapopulation
infected population by the endemic infected population in persistent patches.
This diagnostic asks how much of the actual infection pressure instead comes
from patches in infection episodes that do not persist.

The analysis directly tests the first mechanism proposed in
[`notes/notes_18jul.md`](../notes/notes_18jul.md): ignoring transiently infected
patches may cause the analytical approximation to underestimate colonization
pressure. The note also proposes a separate mechanism, variation in the number
of infected hosts introduced during recolonization and hence in
`P1(R0,K,I(0))`. That second mechanism is not tested here.

## Episode classification

The analysis uses the same model, 21 one-factor-at-a-time scenarios, seeds,
virgin-soil initial conditions, and simulation horizon as the established
occupancy workflow. No odin equation or simulation rule is changed.

For each patch, a maximal contiguous interval with `I > 0` is one infection
episode. The established criterion uses the fixed forward window `tau = 50`
disease generations (`500` stored steps at `dt = 0.1`). If any time point in an
episode has `I > 0` throughout its complete forward window, the entire episode
is classified as persistent from its first infected observation. All other
infected episodes are classified as transient. Classification is retrospective
and episode-level, not recalculated independently at every time point.

An episode that begins within 50 generations of the simulation boundary and
never has a complete qualifying window is classified as transient. This
right-censoring rule can conservatively overestimate the transient contribution
near the end of the trajectory.

At each time, `transient_fraction_I` is the infected hosts in transient episodes
divided by all infected hosts. Because the odin between-patch force of infection
is linear in `sum(I)`, this is also the instantaneous fraction of mean-field
immigration pressure attributable to transient episodes. It does not count
realized successful colonization events.

## Main results

Transient infection episodes can make a large contribution, especially near
the retrospective persistent-occupancy minimum. Among the 20 scenarios with a
defined fraction at that minimum, 13 have at least 50% of infected hosts in
transient episodes. Across the 21 scenario entries, the median defined fraction
at the minimum is 70.5%. The median transient fraction over the first 500
generations after the minimum is 11.4%, and the median over the complete
trajectory is 3.43%. Thus the contribution is concentrated around and shortly
after burnout, rather than being uniformly large through time.

For `K = 1000`, total infection is zero at the persistent-occupancy minimum, so
the fraction there is mathematically undefined. The minimum figure marks this
case separately rather than treating it as zero.

### Baseline R0 versus R0 = 3

The poor analytical fit at `R0 = 3` is associated with a much larger transient
contribution:

| R0 | At occupancy minimum | Mean, first 500 | Mean, full trajectory | max(transient I) / max(persistent I) |
|---:|---:|---:|---:|---:|
| 2.5 | 70.5% | 12.7% | 3.43% | 2.42 |
| 3.0 | 98.0% | 57.0% | 41.7% | 15.10 |

This strongly supports transient transmission as a plausible contributor to
the underprediction at `R0 = 3`. It is an association from one stochastic
trajectory per scenario, not a causal decomposition of the residual.

### Effect of r

The effect of host growth rate is very strong and systematic. As `r` increases
from `0.05` to `0.4`, the transient fraction at the minimum changes from 100%,
94.5%, 70.5%, 0.5%, to 0%. The full-trajectory mean changes from 100% at
`r = 0.05` to zero at `r = 0.4`. The analytical approximation's sensitivity to
`r` is therefore accompanied by major changes in the transient contribution.

### Effects of K and alpha

`alpha` has a comparatively weak effect at the occupancy minimum: the transient
fraction ranges from 70.5% to 77.0% across the five values. Its first-500 mean
is more variable (4.2% to 12.7%) but remains far less changed than in the `r`
series.

`K` does not have a weak effect. The full-trajectory mean transient fraction
declines from 17.0% at `K = 1000` and 16.2% at `K = 3000` to 0.031% at
`K = 100000`. The at-minimum pattern is non-monotone, in part because the
`K = 1000` fraction is undefined at global extinction, but the post-minimum and
overall contribution decreases strongly at large `K`.

## Interpretation

The results support the first mechanism proposed in `notes/notes_18jul.md` in an
important but qualified way.
At `R0 = 3`, low `r`, and around the post-burnout minimum, transient episodes
account for a substantial share of `sum(I)` and hence of mean-field infection
pressure. An approximation using only `I_star` times persistent patches can
therefore underestimate colonization pressure in exactly these cases.

However, transient infection is not a universal missing correction: it is
negligible in some scenarios, is highly dependent on `R0`, `r`, and `K`, and
this analysis does not show that transient immigration establishes new patches.
In particular, the note describes the analytical mismatch as relatively
insensitive to `K`, whereas the transient contribution changes strongly with
`K`. The transient mechanism therefore does not reproduce the complete pattern
of approximation error by itself; compensating dynamics or the proposed
`I(0)`-dependent establishment probability may also matter. The classification
is retrospective and sensitive to the finite horizon.

## Running and outputs

Run from the repository root:

```bash
Rscript fadeout/run_transient_patch_analysis.R
```

All outputs are isolated under `fadeout/output/transient_patch_analysis/`.
The `data/` directory contains the time series as CSV and RDS plus the
one-row-per-scenario summary. The `figures/` directory contains PDF figures for
the four parameter series, the occupancy-minimum comparison, and the dedicated
`R0 = 2.5` versus `R0 = 3` comparison.

## Validation

For all 420,021 scenario-time rows,
`transient_I + persistent_I == total_I` exactly. Whenever `total_I > 0`,
`transient_fraction_I + persistent_fraction_I == 1` exactly within stored
floating-point precision.
