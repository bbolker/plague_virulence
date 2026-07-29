# Fadeout and patch-occupancy analyses

This directory contains the non-seasonal single-strain patch-occupancy
workflow, its analytical approximations, and older related fadeout scripts.
Occupancy means `I > 0` unless explicitly described as established occupancy.

The separate exploratory logistic boundary-layer burnout approximation is
documented in [`logistic_burnout/README.md`](logistic_burnout/README.md).

The current occupancy analyses vary one of `R0`, `K`, `alpha`, or `r` at a
time, with all other parameters held at their baseline values. Figures are
saved as PDF. Run commands from the repository root.

## Current main workflow

### 1. Stochastic patch occupancy

```bash
Rscript fadeout/run_stochastic_occupancy.R
```

This runs one existing-seed stochastic trajectory for each retained parameter
value, using 200 patches, `dt=0.1`, `t_max=2000`, and virgin-soil initial
conditions `S(0)=K-10`, `I(0)=10` in every patch. Each simulation saves its
full patch-level infected history under `data/full_trajectories/`. Raw
occupancy and the episode/establishment analysis are derived from this common
realization, so the metapopulation model is not run twice. Legacy summary-only
caches are still readable, but running the script normally will regenerate
them once to create the full trajectory. Use `--plot-only` to rebuild figures
from existing summary data without rerunning simulations.

| File | Purpose |
|---|---|
| [`run_stochastic_occupancy.R`](run_stochastic_occupancy.R) | Defines the baseline and the `R0`, `K`, `alpha`, and `r` one-factor-at-a-time scenarios, runs each stochastic realization once, and caches its full infected history. |
| [`occupancy_functions.R`](occupancy_functions.R) | Simulation, raw occupancy, established occupancy, episode classification, summary, and absorbing-display helpers. |
| [`output/stochastic_patch_occupancy/`](output/stochastic_patch_occupancy/) | Full stochastic trajectories, direct raw `I>0` occupancy results, compact CSV files, and four parameter-scan PDFs. |

Main stochastic figures:

- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_R0.pdf`
- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_K.pdf`
- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_alpha.pdf`
- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_r.pdf`

### 2. Original one-state approximation

```bash
Rscript fadeout/compare_one_state_occupancy.R
```

The one-state model tracks persistent occupancy $p(t)$:

$$
\frac{dp}{dt}=\alpha I^*P_1p(1-p).
$$

It assumes that (1) the initial epidemic and first burnout are rapid relative
to recolonization, (2) patches surviving that burnout persist thereafter,
(3) persistent patches rapidly reach the deterministic endemic infected
abundance $I^*$, and (4) empty patches rapidly recover to $S=K$. Thus the
infection pressure on an empty patch is

$$
c(t)=\alpha I^*p(t),
$$

and a new colonization establishes persistence with probability $P_1$. With
$p(0)=P_1$ and $\lambda=\alpha I^*P_1$, the logistic solution is

$$
p(t)=\frac{1}{
1+\dfrac{1-P_1}{P_1}\exp\!\left(-\alpha I^*P_1t\right)
}.
$$

`P1` is estimated as a function of `R0`, `K`, and `r` from the single-patch
demography-grid extinction data. Exact grid values are preferred; an
`R0`-by-`K` GAM within the matching observed `r` slice is used only when
interpolation is required. `I_star` comes from `plagueMetapop::ode_eq()`. The
stochastic comparison target is established occupancy with `tau=50`, while
the common post-burnout time origin is the raw occupancy minimum.

| File | Purpose |
|---|---|
| [`compare_one_state_occupancy.R`](compare_one_state_occupancy.R) | Produces the one-state comparisons and diagnostics without rerunning metapopulation simulations. |
| [`output/one_state_occupancy/`](output/one_state_occupancy/) | Four `tau=50` PDFs and the diagnostic CSV. |

Main one-state figures:

- `output/one_state_occupancy/one_state_established_tau50_compare_R0.pdf`
- `output/one_state_occupancy/one_state_established_tau50_compare_K.pdf`
- `output/one_state_occupancy/one_state_established_tau50_compare_alpha.pdf`
- `output/one_state_occupancy/one_state_established_tau50_compare_r.pdf`

### 3. Current two-state approximation

```bash
Rscript fadeout/compare_two_state_occupancy.R
```

The current model distinguishes persistent occupancy $p(t)$ from transient
occupancy $q(t)$:

$$
B(t)=\alpha\left[I^*p(t)+\bar I_Tq(t)\right]
\left[1-p(t)-q(t)\right],
$$

$$
\frac{dp}{dt}=P_1B(t),
$$

$$
\frac{dq}{dt}=(1-P_1)B(t)-\frac{q(t)}{T_T}.
$$

Persistent source patches contribute $I^*$. Transient source patches
contribute the deterministic invasion-trajectory average

$$
\bar I_T=\frac{1}{T_T}\int_0^{T_T}I(t)\,dt,
$$

where

$$
T_T=cT_{\mathrm{osc}},\qquad c=0.5.
$$

where $T_{\textrm{osc}}$ is the (approximate) oscillation period of the SIR model with vital dynamics.
The provisional oscillation closure uses `mu=r` and `rho=gamma`; this is an approximation
because we're modeling logistic susceptible recruitment rather than
standard constant-turnover SIR demography (i.e., linear demography). 

$p(0)=P_1$; $q(0)$ is the observed
episode-classified transient occupancy at the aligned raw minimum.
The deterministic invasion begins at `S(0)=K-10`, `I(0)=10`. 

| File | Purpose |
|---|---|
| [`two_state_functions.R`](two_state_functions.R) | Computes `T_osc`, `T_T`, `Ibar_T`, and solves the current two-state ODE. |
| [`compare_two_state_occupancy.R`](compare_two_state_occupancy.R) | Compares stochastic and deterministic persistent, transient, and total occupancy for 21 retained scenarios. Its default mode uses the original half-period duration; `--hybrid-duration` uses the separate $I=1$-or-first-trough duration. |
| [`output/two_state_occupancy/`](output/two_state_occupancy/) | Deterministic invasion data, comparison curves, diagnostics, and four PDFs. |
| [`output/two_state_occupancy_hybrid_TT/`](output/two_state_occupancy_hybrid_TT/) | Separate curves, diagnostics, and four PDFs generated with the hybrid transient duration; the original results are not overwritten. |

Main current two-state figures:

- `output/two_state_occupancy/figures/two_state_occupancy_compare_R0.pdf`
- `output/two_state_occupancy/figures/two_state_occupancy_compare_K.pdf`
- `output/two_state_occupancy/figures/two_state_occupancy_compare_alpha.pdf`
- `output/two_state_occupancy/figures/two_state_occupancy_compare_r.pdf`

The figures show only the stochastic trajectory and the current deterministic
approximation. The original one-state curve is retained in the two-state CSV
for quantitative comparison but is not drawn.

### Alternative transient-duration closure

The existing two-state implementation remains unchanged and uses the
provisional closure

$$
T_T=0.5T_{\mathrm{osc}}.
$$

[`compute_transient_outbreak_summary_I1()`](two_state_functions.R) provides a
separate deterministic hybrid alternative. It solves the normalized logistic model

$$
\frac{dS}{dt}=rS\left(1-\frac{S}{K}\right)-R_0\frac{SI}{K},
\qquad
\frac{dI}{dt}=R_0\frac{SI}{K}-I,
$$

from $S(0)=K-I_0$, $I(0)=I_0$. The first peak and following trough are
identified by the downward and upward crossings, respectively, of
$S=K/R_0$. If the first-trough infected count is at most one, $T_T$ is the
post-peak downward crossing of $I=1$; otherwise $T_T$ is the first trough:

$$
T_T=
\begin{cases}
t_{\downarrow,I=1}, & I(t_{\mathrm{trough}})\leq 1,\\
t_{\mathrm{trough}}, & I(t_{\mathrm{trough}})>1.
\end{cases}
$$

All event times are linearly interpolated between solver outputs, and the
integration horizon is extended adaptively when needed. If no post-peak trough
is found, the function returns `T_T=NA` with a clear status.

Both closures use the same infected-load calculation

$$
\bar I_T=\frac{1}{T_T}\int_0^{T_T}I(t)\,dt.
$$

The $I=1$ boundary approximates stochastic fade-out; it is not literal
deterministic extinction. The trough branch prevents a deterministic trajectory
whose first trough remains above one from having an undefined duration.
Susceptible population growth remains logistic, not linear demographic
turnover, and all times are in disease generations.

Run the comparison with:

```bash
Rscript fadeout/compare_transient_duration_methods.R
```

Outputs are saved under `output/transient_duration_methods/`, including the
baseline trajectory and comparison table, three one-parameter-at-a-time PDF
diagnostics, and an endpoint-type summary. Each OFAT PDF shows the old and
hybrid $T_T$, the old and hybrid $\bar I_T$, the first-trough infected count,
and the selected hybrid endpoint.
No existing occupancy script uses this new closure by default.

To run the same two-state occupancy comparison with the hybrid closure without
overwriting the original outputs:

```bash
Rscript fadeout/compare_two_state_occupancy.R --hybrid-duration
```

### Initial-infection sensitivity of P1

```bash
Rscript fadeout/check_P1_initial_I.R
```

[`check_P1_initial_I.R`](check_P1_initial_I.R) checks the baseline
single-patch persistence probability at several Poisson mean initial infected
counts from 1 to 300. It uses the
same logistic Euler model and extinction horizon ($t=200$) as the existing
$P_1$ calibration data. Results and their 95% binomial intervals are cached in
`output/P1_initial_I/P1_initial_I_results.csv`; future runs simulate only
missing initial-count values. The figure is
`output/P1_initial_I/P1_by_initial_I_to300.pdf`.

## Supporting established and episode occupancy

```bash
Rscript fadeout/run_stochastic_episode_occupancy.R
```

[`run_stochastic_episode_occupancy.R`](run_stochastic_episode_occupancy.R)
does not run the odin model. It reads the full trajectories saved by
`run_stochastic_occupancy.R` and derives established occupancy plus persistent,
transient, and censored episode classifications for the 21 retained
scenarios. Patch $i$ is
established at time $t$ when

$$
I_i(s)>0
\quad\text{for every recorded }s\in[t,t+\tau],
\qquad \tau=50.
$$

This is an offline, forward-looking definition: the final 50 time units are
`NA`, and with `dt=0.1` the window contains 500 stored steps. A qualifying
episode is persistent from its beginning; an episode ending before it qualifies
is transient; an unresolved terminal episode is censored.

Let $I_T(t)$ and $I_P(t)$ be infected hosts in transient and persistent
episodes. The primary source-pressure diagnostic is

$$
\textrm{transient\ source\ share}(t)
=\frac{I_T(t)}{I_T(t)+I_P(t)}.
$$

Censored infected hosts are excluded from this denominator. Because odin uses
the pooled immigration rate `alpha*sum(I)/n_patch`, this share is the expected
transient contribution to pooled colonization pressure, not realized
source-target ancestry. The cumulative share is

$$
C_T=\int \frac{\alpha I_T(t)}{n_{\mathrm{patch}}} \, dt,
\qquad
C_P=\int \frac{\alpha I_P(t)}{n_{\mathrm{patch}}} \, dt,
$$

with cumulative transient share

$$
\frac{C_T}{C_T+C_P}.
$$

Outputs are under `output/stochastic_episode_occupancy/`; figures omit the final
censored window and show times 0--1950. At the occupancy minimum, transient
source shares are 70.5% for `R0=2.5` and 98.0% for `R0=3`; cumulative shares are
3.42% and 37.4%, respectively. These single-trajectory results support a
transient-pressure explanation for part of the `R0=3` one-state
underprediction, but do not identify transmission ancestry.

## Extinction time versus oscillation period

```bash
Rscript fadeout/analyze_extinction_Tosc_relationship.R
```

[`analyze_extinction_Tosc_relationship.R`](analyze_extinction_Tosc_relationship.R)
compares the current two-state oscillation-period approximation

$$
T_{\mathrm{osc}}=
\frac{2\pi}{
\sqrt{r\gamma(R_0-1)-r^2R_0^2/4}
}
$$

with `mean_ext_time.I1` from the existing 2600-cell one-patch
$R_0\times K\times r$ stochastic grid. A period is defined only when the
quantity under the square root is positive. The current formula depends on
$R_0$ and $r$, not $K$, so all 13 $K$ values at fixed $R_0,r$ have the same
$T_{\mathrm{osc}}$ and remain as separate vertically aligned points. The first
regression is an unweighted fit of

$$
\log_{10}\!\left(\overline T_{\mathrm{ext}}\right)
=a+b\log_{10}\!\left(T_{\mathrm{osc}}\right)
$$

over all valid grid cells; $K$ is shown by point colour but is not included as
a predictor or averaged out. The regression and figures are restricted to
$0<T_{\mathrm{osc}}<50$ disease generations. Longer periods are not used
because extinction is observed only through `t_max=200`, so their conditional
mean extinction times are especially vulnerable to finite-window selection.
All 2600 grid cells remain in the cleaned output table, with a flag identifying
whether they satisfy this period-window restriction.

`sumfun_discrete()` stores the mean row index of first extinction among runs
that became extinct. Because row 1 is time zero and rows are `dt=0.1` apart,
the analysis uses

$$
\overline T_{\mathrm{ext}}
=\left(\mathtt{mean\_ext\_time.I1}-1\right)dt.
$$

This is the conditional mean extinction time among runs that went extinct
within `t_max=200`, not a replicate-level extinction-time distribution. The
combined data contain only grid-cell summaries. This first analysis does not
separate early fizzles from post-outbreak burnout, which is especially
important near $R_0=1$.

Outputs are under `output/extinction_Tosc_relationship/`: the full 2600-row
cleaned grid, regression and ratio-summary CSV files, and two PDF figures in
the `figures/` subdirectory.

## One-patch extinction: one parameter at a time

```bash
Rscript fadeout/analyze_onepatch_extinction_sensitivity.R
```

[`analyze_onepatch_extinction_sensitivity.R`](analyze_onepatch_extinction_sensitivity.R)
uses the existing one-patch, one-strain logistic continuous-Euler extinction
grid to show how the stochastic summaries change around the baseline
$R_0=2.5$, $K=10^4$, and $r=0.125$. Exactly three sweeps are retained: vary
$R_0$, vary $K$, or vary $r$, while holding the other two parameters fixed.
There is no $\alpha$ sweep because this is a one-patch dataset with
$\alpha=0$.

Each two-panel figure shows the conditional mean extinction time above and the
probability of extinction by `t_max=200` below. The latter indicates where
conditioning and finite-window censoring make the mean harder to interpret.
As above, the time conversion is

$$
\overline T_{\mathrm{ext}}
=\left(\mathtt{mean\_ext\_time.I1}-1\right)dt,
\qquad dt=0.1,
$$

in disease generations. Outputs are under
`output/onepatch_extinction_sensitivity/`: a cleaned full-grid CSV, three sweep CSV
files, and exactly three PDF figures.

## Seasonal recurrent fade-out

The reusable seasonal workflow is documented in
[`seasonal/README.md`](seasonal/README.md). It extracts infection episodes,
classifies post-burnout fade-out using a configurable multiple of the
non-seasonal intrinsic period, calculates annual occupancy diagnostics, and
produces full and conditional Kaplan-Meier episode-survival curves. It supports
explicitly seeded replicates and parameter-grid CSV files; outputs are written
under `output/seasonal_fadeout/<run-id>/`.
