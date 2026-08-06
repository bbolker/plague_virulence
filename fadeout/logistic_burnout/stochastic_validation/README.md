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
| `figures/validation_conditional_probability_round_vs_continuous.{png,pdf}` | Persistence conditional on escaping fizzle (`P_cond_sim_established` vs. `P_cond_approx`), with the rounded-`m` and continuous-`m` analytical curves overlaid. |
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

Two follow-up checks (not in the original committed run) narrow down where
the small-\(K\)/high-\(R_0\) discrepancy comes from:

- **Rounding \(m=Ky^*\) to an integer is not the cause.** `analytical_validation_quantities()`
  now also returns a continuous-\(m\) analytical curve (\(q_1^{m_{\mathrm{raw}}}\),
  matching Eq. 22/26/27 of Parsons et al. 2024, which do not round), stored in
  `Q_approx_continuous`/`P_cond_approx_continuous`/`P_uncond_approx_continuous`
  and plotted in `validation_conditional_probability_round_vs_continuous`. The
  rounded and continuous curves are visually indistinguishable at every \(K\)
  tested; the discrepancy with simulation must come from elsewhere (most
  likely the neglected susceptible-recovery stochasticity noted above).
- **The plotted conditional-persistence curve pins the discrepancy to the
  post-peak, high-\(R_0\) transition specifically** (roughly \(R_0>2\)), not to
  the fizzle-dominated low-\(R_0\) region or the near-1 plateau, both of which
  match simulation closely once fizzle-escape uncertainty (wide CIs from small
  samples) is accounted for.
- Near \(R_0=1\), a nontrivial fraction of full-trajectory simulations end the
  run still classified `unresolved` rather than `persistence`/`fizzle`,
  because near-critical growth from \(I_0=1\) can take
  \(\sim\!\ln(K)/(R_0-1)\) disease generations to reach the boundary-layer
  scale — longer than `tmax=200` for \(R_0\) within roughly \(\ln(K)/t_{\max}\)
  of 1. This inflates the apparent gap between simulated and reference
  unconditional persistence right at the left edge of each panel; it is a
  simulation-horizon artifact, not a tau-leap step-size (`dt`) error. Rerunning
  the affected cells (\(R_0\lesssim1.06\) at \(K=10{,}000\) and \(30{,}000\))
  with \(t_{\max}\) set to roughly \(4\times\ln(K)/(R_0-1)\) reduces the
  maximum unresolved fraction across the whole grid from 0.052 to 0.0016.
- Separately, the deterministic ODE used to locate the epidemic peak can be
  genuinely **overdamped** rather than merely slow when \(R_0\) is close
  enough to 1 (specifically when \(r(R_0-1) < r^2R_0^2/4\), the same
  discriminant as the \(T_{\mathrm{osc}}\) closure elsewhere in this
  repository): \(x(t)\) then approaches \(x_*\) monotonically and never
  crosses it, however long the ODE is integrated (checked out to \(t=25{,}600\)
  for \(R_0=1.02\), \(K\in\{1000,30000\}\)). `logistic_burnout_probability()`'s
  analytical curve is therefore undefined there for a structural reason, not
  merely a short integration horizon — extending `initial_tmax`/`maximum_tmax`
  recovers the analytical curve only for cells above this overdamped threshold
  (confirmed for \(R_0=1.05\), where the peak appears once the horizon reaches
  ~200 generations, since the current retry loop only re-extends the horizon
  when a peak was found but no post-peak crossing was, not when no peak was
  found at all).
- With the unresolved-fraction artifact removed, unconditional persistence in
  the fully resolved near-\(R_0=1\) simulations still comes in systematically
  **below** the early-establishment reference \(p_{\mathrm{est}}=1-(1/R_0)^{I_0}\)
  (e.g. 0.0010 vs. 0.0196 at \(R_0=1.02\), \(K=30{,}000\); 0.0426 vs. 0.0476 at
  \(R_0=1.05\); gap shrinking quickly as \(R_0\) moves away from 1). This is
  consistent with \(p_{\mathrm{est}}\) implicitly assuming a permanently
  supercritical branching environment at the DFE, whereas \(x(t)\) actually
  drifts down toward \(x_*\) (and, in the overdamped regime, never back up),
  so the effective local reproduction number \(R_0x(t)\) relaxes toward
  exactly 1 — a genuinely critical, not supercritical, branching process,
  which is null-recurrent (extinguishes with probability 1 given enough time).
  \(p_{\mathrm{est}}\) is therefore an overestimate of true escape probability
  near \(R_0=1\), not just an approximation with unmodelled noise.
- **The genuinely-\(0\)-looking simulated points at small \(K\) (e.g.
  \(K=1000\), \(R_0\lesssim1.1\)) are real, resolved results, not an
  artifact** — but they reveal that comparing against the *first-trough-only*
  analytical curve (what this module validates) is the wrong comparison at
  small \(K\), not that the simulator is broken. E.g. at \(K=1000\),
  \(R_0=1.1\): \(m_{\mathrm{used}}=8\) is small enough that each boundary-layer
  encounter only has \(\approx35\%\) conditional persistence probability
  (`logistic_multitrough_probabilities()`), and successive encounters recur
  roughly every 75 generations; cumulative persistence through 4 troughs
  drops to 0.016, vs. 0.346 after the first trough alone — a 20-fold gap
  documented in `../validate_multitrough_burnout.R`'s own output, not
  discovered here for the first time. The single-trough
  \(P_{\mathrm{uncond,approx}}=0.0315\) reported by this module's
  `analytical_validation_quantities()` is the *first-trough* quantity; the
  simulation (which runs until true extinction or `tmax`, i.e. implicitly
  through as many troughs as occur) is closer to the multi-trough cumulative
  quantity. This module does not yet compare against the multi-trough
  cumulative persistence — doing so would likely close most of the apparent
  small-\(K\) gap and is a natural extension (`logistic_multitrough_probabilities()`
  already exists in `../logistic_burnout_functions.R`; it would need wiring
  into `analytical_validation_quantities()`/the plotting scripts here).

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
