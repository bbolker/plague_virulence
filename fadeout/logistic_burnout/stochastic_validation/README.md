# Stochastic validation of logistic burnout

This directory validates the existing first-trough logistic burnout
approximation against complete stochastic single-patch simulations. It is a
validation workflow, not a redesign of the approximation in
`../logistic_burnout_functions.R`.

## Scientific target

Parsons et al. (2024) validate an approximation that separates early stochastic
fizzle, a deterministic major epidemic, and later stochastic burnout near the
post-epidemic boundary layer. The existing logistic approximation in this
repository already validates the boundary-layer birth-death calculation, but
that alone does not test the complete stochastic susceptible-infected process
from the original introduction.

Here the normalized deterministic model is

$$
\frac{dx}{dt}=rx(1-x)-R_0xy,\qquad
\frac{dy}{dt}=(R_0x-1)y,
$$

with endemic state

$$
x^*=\frac{1}{R_0},\qquad
y^*=\frac{r(R_0-1)}{R_0^2}.
$$

The approximation uses the first post-peak downward crossing of

$$
y_{\mathrm{BL}}=y^*
$$

and the existing `logistic_burnout_probability()` function. No boundary
multiplier is introduced.

## Stochastic process

The primary validation simulates a complete non-seasonal, single-patch
logistic S-I process from

$$
S(0)=K-I_0,\qquad I(0)=I_0.
$$

The fixed-step tau-leap diagnostic engine follows the event conventions in
`plagueMetapop/inst/odin/euler_odin_def.R`:

| Event | Implementation |
|---|---|
| Infection | `rbinom(S, 1-exp(-R0*I*dt/K))` |
| Infected removal | `rbinom(I, 1-exp(-dt))` |
| Positive logistic susceptible growth | `rpois(r*S*(1-S/K)*dt)` |
| Negative logistic susceptible growth | `-rbinom(S, -r*S*(1-S/K)*dt/S)` |

The primary validation engine is a hybrid adaptive tau-leap implementation:
large-population phases use adaptive Poisson tau steps, while low-infection or
low-rate phases use exact Gillespie events. This mirrors the reason Parsons et
al. used adaptive tau-leaping for large stochastic simulations while keeping
the first-trough and extinction phases close to the continuous-time Markov
process. The fixed-step engine remains available only as a diagnostic because
it matches the repository's odin/euler implementation. The exact reference
engine is a Gillespie simulation of the corresponding continuous-time
transition rates and is used for representative checks.

## Fizzle criterion

For \(R_0>1\), the early fizzle threshold time is

$$
\tau_\delta
=\frac{1}{R_0-1}
\log
\left[
\frac{(1-\delta)^{-1/I_0}-1/R_0}
{(1-\delta)^{-1/I_0}-1}
\right].
$$

The default is \(\delta=10^{-6}\). A trajectory that reaches \(I=0\) before
\(\tau_\delta\) is classified as `fizzle`.

## Outcome state machine

Each full trajectory is classified into one of five outcomes:

| Outcome | Definition |
|---|---|
| `fizzle` | \(I=0\) before \(\tau_\delta\). |
| `late_extinction_before_boundary` | Escapes fizzle but reaches \(I=0\) before a valid first downward entry into \(I\le I_{\mathrm{BL}}\). |
| `burnout` | Escapes fizzle, rises above \(I_{\mathrm{BL}}\), enters downward into \(I\le I_{\mathrm{BL}}\), then reaches \(I=0\) before upward exit. |
| `persistence` | After first downward boundary entry, leaves the boundary layer upward before extinction. |
| `unresolved` | Reaches the maximum simulation horizon without one of the above outcomes. |

The boundary count uses the same convention as the analytical code:

$$
I_{\mathrm{BL}}=m_{\mathrm{used}}
=\max\{1,\mathrm{round}(Ky^*)\}.
$$

A one-step fluctuation at the boundary is not counted as a later wave unless a
confirmed downward entry has already occurred.

## Analytical and simulation estimands

The analytical code returns:

| Quantity | Meaning |
|---|---|
| \(q_1\) | One-lineage extinction probability in the boundary layer. |
| \(Q_{\mathrm{approx}}=q_1^{m_{\mathrm{used}}}\) | Conditional whole-patch burnout probability at the first trough. |
| \(P_{\mathrm{cond,approx}}=1-Q_{\mathrm{approx}}\) | Conditional persistence through the first trough. |
| \(p_{\mathrm{est,approx}}=1-(1/R_0)^{I_0}\) | Early branching approximation for escaping fizzle. |
| \(P_{\mathrm{uncond,approx}}=p_{\mathrm{est,approx}}P_{\mathrm{cond,approx}}\) | Approximate unconditional persistence from the original introduction. |

The full simulation reports two conditional estimates:

| Quantity | Conditioning |
|---|---|
| `P_cond_sim_established` | Persistence divided by all trajectories that escaped fizzle. |
| `P_cond_sim_boundary` | Persistence divided by trajectories classified as persistence or burnout after boundary entry. |

These differ when `late_extinction_before_boundary` is non-negligible.

## Boundary-start validation

A separate validation starts the stochastic model at

$$
S_0=\mathrm{round}(Kx_{\mathrm{in}}),\qquad I_0=m_{\mathrm{used}},
$$

and treats this as first boundary entry. It estimates the probability of
extinction before upward boundary exit and compares it directly with
\(Q_{\mathrm{approx}}\). This isolates the boundary-layer approximation from
early fizzle and first-wave stochasticity, while still retaining stochastic
susceptible dynamics.

## Parameter grids

The smoke test uses \(R_0=\{1.1,2.5,5\}\), \(r=0.1\),
\(K=\{1000,10000\}\), \(I_0=1\), and 100-200 simulations per cell by default.
It is only a code check.

The main validation fixes \(r=0.1\), uses
\(K=\{1000,3000,10000,30000\}\), and uses a dense grid in \(R_0-1\) from
near 0.02 to 4, including explicit values near 1.05, 1.1, 1.5, 2, 2.5, 3, 4,
and 5.

The secondary sensitivity grid uses
\(r=\{0.05,0.1,0.125,0.2\}\),
\(K=\{3000,10000\}\), and
\(R_0=\{1.05,1.1,1.5,2,2.5,3,4,5\}\).

Simulation counts are adaptive. The default target is a 95% Wilson interval
width of 0.01 for unconditional persistence, subject to minimum and maximum
simulation counts.

## Outputs

| File | Purpose |
|---|---|
| `stochastic_validation_functions.R` | Shared analytical extraction, stochastic engines, outcome classification, confidence intervals, and checks. |
| `run_smoke_test.R` | Short local code-path test. |
| `run_stochastic_validation.R` | Resumable adaptive validation runner. |
| `run_engine_checks.R` | Fixed-step tau-leap, adaptive tau-leap, and Gillespie comparison on representative cells. |
| `plot_stochastic_validation.R` | Generates the focused validation figures from available CSV outputs. |
| `outputs/full_stochastic_validation.csv` | Full simulation estimates and analytical quantities. |
| `outputs/boundary_start_validation.csv` | Boundary-start estimates compared with \(Q_{\mathrm{approx}}\). |
| `outputs/simulation_engine_comparison.csv` | Gillespie, adaptive tau-leap, and fixed-step tau-leap comparison. |
| `outputs/dt_sensitivity.csv` | Tau-leap timestep sensitivity. |
| `outputs/outcome_counts.csv` | Outcome count table. |
| `outputs/validation_status_summary.csv` | Completed-cell and error summary. |
| `outputs/smoke_test_results.csv` | Smoke-test full simulation results. |
| `figures/validation_unconditional_probability_scale_legend.{png,pdf}` | Fig.-4-style unconditional persistence comparison, with explicit legend. |
| `figures/validation_simulator_sensitivity.{png,pdf}` | Tau-leap and Gillespie comparison. |

Additional CSV outputs retain conditional, boundary-start, error, and
outcome-count diagnostics. They are not plotted by default so that the primary
outputs remain easy to inspect.

## Reproducibility

From the repository root:

```bash
Rscript fadeout/logistic_burnout/stochastic_validation/run_smoke_test.R
Rscript fadeout/logistic_burnout/stochastic_validation/run_engine_checks.R
Rscript fadeout/logistic_burnout/stochastic_validation/run_stochastic_validation.R --mode=main
Rscript fadeout/logistic_burnout/stochastic_validation/run_stochastic_validation.R --mode=r_sensitivity --min_sim=2000 --max_sim=10000
Rscript fadeout/logistic_burnout/stochastic_validation/plot_stochastic_validation.R
```

The main validation is resumable: completed cells in
`outputs/full_stochastic_validation.csv` are skipped unless `--overwrite=TRUE`
is supplied.

## Limitations

The stochastic logistic susceptible process is a modeling choice matched to
the repository's odin transition convention. The analytical approximation
uses deterministic \(x_{\mathrm{in}}\) and neglects susceptible stochasticity
inside the boundary layer. Near \(R_0=1\), the distinction among fizzle, a
major epidemic, and burnout is weak. Small \(Ky^*\) is expected to reduce
accuracy. Finite simulations cannot resolve arbitrarily small probabilities.

## Generated results

The current generated results use the main \(r=0.1\) grid with the default
minimum of 5,000 tau-leap realizations per cell, but the run was capped at
5,000 per cell rather than adaptively continuing to 50,000. The settings were
`dt = 0.02`, `delta = 1e-6`, `tmax = 200`, `r = 0.1`, and \(I_0=1\). The run
completed 124 \(R_0,K\) cells and 620,000 full stochastic realizations.

For the 84 cells with defined analytical values, the mean absolute
unconditional persistence error was 0.0446 and the maximum absolute error was
0.1493. The largest discrepancies include \(R_0=2.26426,K=30000\), where the
simulation estimate was 0.325 and the approximation was 0.4743, and
\(R_0\simeq3,K=1000\), where the simulation estimate was about 0.12 while the
approximation was near zero.

The approximation is therefore not uniformly accurate in this preliminary
full-trajectory validation. It is closer through part of the intermediate
\(R_0\) range, especially for larger \(K\), but it tends to miss nonzero
stochastic persistence in some high-\(R_0\) cells and can overestimate
persistence in some larger-\(K\), intermediate-\(R_0\) cells. These are
scientific discrepancies, not changes to the analytical formula.

Late extinction before formal boundary entry reached 0.3454,
showing that the deterministic establishment-to-boundary step can fail in some
parameter regions. The maximum unresolved fraction in full trajectories was
0.052.

The boundary-start run completed all 124 cells, but 40 near-critical cells
were unresolved at `tmax = 200` because no deterministic boundary-start
comparison was identifiable there. Among the 84 cells with a defined
boundary-start comparison, the mean absolute error in \(Q\) was 0.1143 and the
maximum absolute error was 0.5942. This indicates that susceptible
stochasticity and boundary-exit classification can still matter even after
starting at the deterministic boundary state.

The engine check used five representative parameter cells. Tau-leap and
Gillespie agreed well near \(R_0=1.5,K=10000\), with estimates around
0.32-0.35 and overlapping Monte Carlo intervals. Near \(R_0=2.5\), tau-leap
was sensitive to `dt` and lower than Gillespie: for \(K=10000\), tau-leap
estimates were 0.040, 0.061, and 0.103 for `dt = 0.05`, 0.02, and 0.01,
whereas Gillespie estimated 0.255. This should be treated as a serious
timestep/simulator diagnostic before relying on the tau-leap preview in that
region.
