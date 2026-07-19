# Transient Source-Pressure Audit

## Scope and audit outcome

This audit clarifies the transient-patch diagnostic proposed in
[`notes/notes_18jul.md`](../notes/notes_18jul.md). The estimable quantity is the
contribution of transient infection episodes to pooled metapopulation
colonization pressure. It is not source-target transmission ancestry.

The pre-audit implementation correctly:

1. identified each contiguous `I > 0` interval as an independent infection
   episode, so one patch could have multiple transient and persistent episodes;
2. classified a qualifying persistent episode as persistent from the beginning
   of that episode.

It incorrectly classified a non-qualifying episode still infected at the end of
the simulation as transient. The corrected analysis classifies such episodes as
censored and excludes them from the transient-versus-persistent denominator.
The corrected outputs are stored directly under
`output/transient_patch_analysis/`. Superseded pre-audit outputs were removed
after the audit was completed.

## Corrected episode classification

The persistence window is `tau = 50` disease generations, or 500 stored steps
at `dt = 0.1`.

For each independent infection episode:

- **persistent:** the episode contains at least one time point for which `I > 0`
  throughout the complete forward persistence window; the whole episode is
  persistent from its first infected observation;
- **transient:** the episode returns to `I = 0` before ever satisfying the
  criterion;
- **censored:** the episode is still infected at the simulation boundary and
  has not had enough observed future time to satisfy the criterion.

At every recorded time,

\[
I_{\mathrm{classifiable}}=I_{\mathrm{transient}}+I_{\mathrm{persistent}},
\]

and

\[
I_{\mathrm{total}}=I_{\mathrm{classifiable}}+I_{\mathrm{censored}}.
\]

The primary instantaneous diagnostic is

\[
\text{transient source share}(t)=
\frac{I_{\mathrm{transient}}(t)}
{I_{\mathrm{transient}}(t)+I_{\mathrm{persistent}}(t)}.
\]

Censored infected hosts are excluded from this denominator and are reported
separately. `transient_share_of_total_I` is retained as a secondary quantity.
For clarity, all time-series figures omit the final `tau = 50` time units and
end at time 1950. Every episode visible in that plotting interval has enough
future observation to be classified, so `censored_I = 0` there. The saved data
retain the complete trajectory and the explicit censored component.

## What colonization quantity is estimable?

The odin model calculates

```r
foi[] <- alpha * sum(I[,i]) / n_patch
immig[,] <- rpois(foi[j] * dt)
```

Thus the source pool is aggregated before immigration is drawn. No source patch
identifier is stored for an immigration event. Because the instantaneous
immigration rate is linear in pooled `sum(I)`, `transient_source_share` is the
expected fraction of classifiable instantaneous colonization pressure generated
by transient source episodes.

The current output cannot identify exactly:

- the fraction of new infection episodes caused by transient source patches;
- the fraction of newly persistent episodes caused by transient source patches.

Those are realized ancestry quantities and would require changing the simulator
to record source-target events, followed by new simulations. This audit does not
infer them from pooled immigration counts and does not modify the odin model.

## Integrated source pressure

For each scenario the analysis evaluates trapezoidal time integrals:

\[
C_T=\int \frac{\alpha I_T(t)}{n_{\mathrm{patch}}}\,dt,
\qquad
C_P=\int \frac{\alpha I_P(t)}{n_{\mathrm{patch}}}\,dt,
\]

and

\[
\overline{q}_T=\frac{C_T}{C_T+C_P}.
\]

Within a scenario, `alpha/n_patch` is common to both terms. The pressure-based
share therefore equals the ratio based on integrated transient and persistent
`I`. Across all 21 scenarios, the maximum numerical difference between these
two calculations is `1.11e-16`.

## Results

The scientific pattern from the pre-audit analysis is retained. Among the 20
scenarios with a defined source share at the persistent-occupancy minimum, 13
have a transient source share of at least 50%; the median is 70.5%. The median
mean share over the first 500 generations after that minimum is 11.4%. The
median cumulative transient source-pressure share is 3.42%.

The dedicated `R0` comparison is:

| R0 | Source share at minimum | Mean source share, first 500 | Mean source share, full trajectory | Cumulative source-pressure share |
|---:|---:|---:|---:|---:|
| 2.5 | 70.5% | 12.7% | 3.43% | 3.42% |
| 3.0 | 98.0% | 57.0% | 41.6% | 37.4% |

For `R0 = 3`, cumulative transient and persistent pressures are 9.26 and 15.47,
respectively, in the model's integrated per-patch pressure units. This supports
the hypothesis that ignoring transient sources can substantially underestimate
pooled colonization pressure in the poorly fitted `R0 = 3` scenario.

The `r` effect remains strong: cumulative transient source-pressure shares for
`r = 0.05, 0.1, 0.125, 0.2, 0.4` are approximately 100%, 26.5%, 3.42%, 0.0059%,
and 0%. The cumulative share declines strongly with `K`, from 96.0% at
`K = 1000` to 0.052% at `K = 100000`. Across the `alpha` series it changes more
moderately, from 5.77% at `alpha = 1e-5` to 1.69% at `alpha = 1e-3`.

## Censoring impact

Five of the 21 scenario entries contain censored infected episodes. Censored
infection occurs in 1,916 of 420,021 scenario-time rows; the largest censored
infected population at one time is 14,627. Censoring is concentrated near the
simulation boundary. It does not affect the baseline scenario and changes the
`R0 = 3` full-trajectory mean source share only slightly relative to the
pre-audit calculation, but separating it is necessary for a correct estimand.

## Validation

For every one of 420,021 scenario-time rows:

- `transient_I + persistent_I == classifiable_I` exactly;
- `classifiable_I + censored_I == total_I` exactly;
- when `classifiable_I > 0`, transient and persistent source shares sum to one
  exactly within stored precision.

The 21 reproduced raw trajectories also match the existing occupancy outputs
point by point. No analytical approximation, odin equation, or previous output
was modified.

## Outputs

Run from the repository root with:

```bash
Rscript fadeout/run_transient_patch_analysis.R
```

Corrected files are in `fadeout/output/transient_patch_analysis/`:

- `data/transient_source_pressure_timeseries.csv`
- `data/transient_source_pressure_timeseries.rds`
- `data/transient_source_pressure_summary.csv`
- `figures/transient_source_share_R0_2p5_vs_3.pdf`
- `figures/infected_components_R0_2p5_vs_3.pdf`

Only the representative `R0 = 2.5` and `R0 = 3` comparison is plotted. Both
figures show `t = 0--1950`; the final 50-unit censored window is omitted. The
component figure uses linear-scale stacked areas, with `total_I` as the black
upper boundary.
