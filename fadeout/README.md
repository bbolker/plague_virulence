# Fadeout and patch-occupancy analyses

This directory contains the non-seasonal single-strain patch-occupancy
workflow, its analytical approximations, and older related fadeout scripts.
Occupancy means `I > 0` unless explicitly described as established occupancy.

The current occupancy analyses fix `r = 0.125` and vary one of `R0`, `K`, or
`alpha` at a time. Figures are saved as PDF. Run commands from the repository
root.

## Current main workflow

### 1. Stochastic patch occupancy

```bash
Rscript fadeout/run_stochastic_occupancy.R
```

This runs one existing-seed stochastic trajectory for each retained parameter
value, using 200 patches, `dt=0.1`, `t_max=2000`, and virgin-soil initial
conditions `S(0)=K-10`, `I(0)=10` in every patch. Use `--plot-only` to rebuild
the saved data and figures from existing trajectories without rerunning the
simulations.

| File | Purpose |
|---|---|
| [`run_stochastic_occupancy.R`](run_stochastic_occupancy.R) | Defines the baseline and the `R0`, `K`, and `alpha` one-factor-at-a-time scenarios. |
| [`occupancy_functions.R`](occupancy_functions.R) | Simulation, raw occupancy, established occupancy, episode classification, summary, and absorbing-display helpers. |
| [`output/stochastic_patch_occupancy/`](output/stochastic_patch_occupancy/) | Direct raw `I>0` stochastic patch-occupancy results, compact CSV files, and three parameter-scan PDFs. |

Main stochastic figures:

- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_R0.pdf`
- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_K.pdf`
- `output/stochastic_patch_occupancy/figures/stochastic_patch_occupancy_compare_alpha.pdf`

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

`P1` is predicted as a function of `R0` and `K` by the existing GAM fitted to
the single-patch extinction data; those data use `r=0.125`. `I_star` comes from
`plagueMetapop::ode_eq()`. The stochastic comparison target is established
occupancy with `tau=50`, while the common post-burnout time origin is the raw
occupancy minimum.

| File | Purpose |
|---|---|
| [`compare_one_state_occupancy.R`](compare_one_state_occupancy.R) | Produces the one-state comparisons and diagnostics without rerunning metapopulation simulations. |
| [`output/one_state_occupancy/`](output/one_state_occupancy/) | Three `tau=50` PDFs and the diagnostic CSV. |

Main one-state figures:

- `output/one_state_occupancy/one_state_established_tau50_compare_R0.pdf`
- `output/one_state_occupancy/one_state_established_tau50_compare_K.pdf`
- `output/one_state_occupancy/one_state_established_tau50_compare_alpha.pdf`

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

The deterministic invasion begins at `S(0)=K-10`, `I(0)=10`. The provisional
oscillation closure uses `mu=r` and `rho=gamma`; this is an approximation
because the repository model has logistic susceptible recruitment rather than
standard constant-turnover SIR demography. $p(0)=P_1$; $q(0)$ is the observed
episode-classified transient occupancy at the aligned raw minimum.

| File | Purpose |
|---|---|
| [`two_state_functions.R`](two_state_functions.R) | Computes `T_osc`, `T_T`, `Ibar_T`, and solves the current two-state ODE. |
| [`compare_two_state_occupancy.R`](compare_two_state_occupancy.R) | Compares stochastic and current deterministic persistent, transient, and total occupancy for 16 retained scenarios. |
| [`output/two_state_occupancy/`](output/two_state_occupancy/) | Deterministic invasion data, comparison curves, diagnostics, and three PDFs. |

Main current two-state figures:

- `output/two_state_occupancy/figures/two_state_occupancy_compare_R0.pdf`
- `output/two_state_occupancy/figures/two_state_occupancy_compare_K.pdf`
- `output/two_state_occupancy/figures/two_state_occupancy_compare_alpha.pdf`

The figures show only the stochastic trajectory and the current deterministic
approximation. The original one-state curve is retained in the two-state CSV
for quantitative comparison but is not drawn.

### Initial-infection sensitivity of P1

```bash
Rscript fadeout/check_P1_initial_I.R
```

[`check_P1_initial_I.R`](check_P1_initial_I.R) checks the baseline
single-patch persistence probability at several Poisson mean initial infected
counts, from 10 to the endemic equilibrium infected count $I^*$. It uses the
same logistic Euler model and extinction horizon ($t=200$) as the existing
$P_1$ calibration data. The single figure is saved directly as
`output/P1_initial_I_effect.pdf`.

## Supporting established and episode occupancy

```bash
Rscript fadeout/run_stochastic_episode_occupancy.R
```

[`run_stochastic_episode_occupancy.R`](run_stochastic_episode_occupancy.R)
runs the 16 retained scenarios and produces established occupancy plus
persistent, transient, and censored episode classifications. Patch $i$ is
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
\operatorname{transient\ source\ share}(t)
=\frac{I_T(t)}{I_T(t)+I_P(t)}.
$$

Censored infected hosts are excluded from this denominator. Because odin uses
the pooled immigration rate `alpha*sum(I)/n_patch`, this share is the expected
transient contribution to pooled colonization pressure, not realized
source-target ancestry. The cumulative share is

$$
C_T=\int \frac{\alpha I_T(t)}{n_{\mathrm{patch}}}\,dt,
\qquad
C_P=\int \frac{\alpha I_P(t)}{n_{\mathrm{patch}}}\,dt,
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
