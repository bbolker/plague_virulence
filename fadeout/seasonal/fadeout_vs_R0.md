# Why does post-establishment fade-out increase with $R_0$?

2026-10-06

## The observation

In the batch runs (`occupancy_factorial_100y_batch.R`), a patch counts
as *established* once its episode has lasted $1.5\,T_0$, where $T_0$ is
the intrinsic period of the damped endemic oscillation. Under constant
transmission, the hazard of fade-out after establishment rises steeply
with $R_0$ at every $K$:

|  R0 | K = 1000 | K = 3000 | K = 5000 | K = 10000 |
|----:|---------:|---------:|---------:|----------:|
| 1.2 |  1.05834 |  0.02829 |  0.00079 |   0.00000 |
| 1.5 |  0.55282 |  0.00334 |  0.00001 |   0.00000 |
| 2.0 |  0.64936 |  0.00634 |  0.00002 |   0.00000 |
| 2.5 |  0.96347 |  0.03310 |  0.00088 |   0.00001 |
| 3.0 |  1.39292 |  0.10871 |  0.00807 |   0.00029 |
| 3.5 |  1.67660 |  0.23269 |  0.03590 |   0.00111 |
| 5.0 |  2.62388 |  0.85070 |  0.33116 |   0.04172 |

Fade-outs per established patch-year (constant transmission, alpha =
1e-4/day, 5 replicates pooled)

## Candidate explanations

1.  **Transient ringing after establishment.** An established patch is
    still recovering from the initial overshoot. At higher $R_0$ the
    overshoot is larger and is damped more slowly, so the 2nd and 3rd
    troughs are deeper. This is the $P_2, P_3, \ldots$ effect: the
    probability of surviving the *second and later* troughs, not just
    the first. It predicts fade-outs concentrated at young ages (a few
    $T_0$).
2.  **Stochastic amplification near equilibrium.** Long after the
    transient has gone, demographic noise keeps exciting the weakly
    damped oscillatory mode (quasi-cycles; Alonso, McKane & Pascual
    2007). Weaker damping means larger sustained fluctuations and
    occasional troughs deep enough to lose the infection. This predicts
    a roughly constant hazard at all ages.
3.  **Smaller endemic prevalence.** A lower $I^*$ at higher $R_0$ would
    raise the hazard. This is ruled out:
    $I^*/K = (r/\gamma)(1 - 1/R_0)/R_0$ is nearly flat across the range
    of $R_0$ (table below).

A minor additional factor is that $T_0$ is shorter at high $R_0$, so
there are more troughs per year. That accounts for at most a factor of
about 1.5, not orders of magnitude.

## Linear analysis

Linearizing the logistic SI model about its endemic equilibrium gives a
damped oscillation with damping rate $r/(2R_0)$ and natural frequency
$\omega_0^2 = \gamma r (R_0 - 1)/R_0 - r^2/(4R_0^2)$. Damping *per
cycle* weakens as $R_0$ increases, because the decay rate falls while
the frequency rises:

| R0 | damping rate (/day) | omega0 (/day) | T0 (days) | amplitude retained per cycle | I\*/K |
|---:|---:|---:|---:|---:|---:|
| 1.2 | 0.0104 | 0.0269 | 233 | 0.09 | 0.0174 |
| 1.5 | 0.0083 | 0.0400 | 157 | 0.27 | 0.0278 |
| 2.0 | 0.0063 | 0.0496 | 127 | 0.45 | 0.0312 |
| 2.5 | 0.0050 | 0.0545 | 115 | 0.56 | 0.0300 |
| 3.0 | 0.0042 | 0.0576 | 109 | 0.63 | 0.0278 |
| 3.5 | 0.0036 | 0.0597 | 105 | 0.69 | 0.0255 |
| 5.0 | 0.0025 | 0.0632 | 99 | 0.78 | 0.0200 |

Both explanation 1 and explanation 2 follow from this weakening damping.
A slowly damped transient rings for longer, and the variance of
noise-driven quasi-cycles scales inversely with the damping rate.

## Young versus older fade-outs

To separate the two mechanisms, split the time after establishment into
a young window ($1.5$–$3\,T_0$, around the second trough) and everything
later. Then compute the hazard in each window (constant transmission;
$K$ = 5000 and 10000; all $\alpha$ pooled). Time at risk in the young
window comes from the episode-duration histograms ($0.1\,T_0$ bins,
using bin midpoints).

| R0 | hazard, 1.5-3 T0 (/yr) | hazard, \> 3 T0 (/yr) | ratio | fraction of fade-outs \> 3 T0 | total fade-outs |
|---:|---:|---:|---:|---:|---:|
| 1.2 | 0.00029 | 0.000364 | 1 | 0.99 | 503 |
| 1.5 | 0.00011 | 0.000010 | 11 | 0.93 | 15 |
| 2.0 | 0.00205 | 0.000020 | 105 | 0.64 | 42 |
| 2.5 | 0.02936 | 0.000394 | 75 | 0.73 | 740 |
| 3.0 | 0.11774 | 0.003154 | 37 | 0.81 | 4775 |
| 3.5 | 0.27207 | 0.010928 | 25 | 0.80 | 13337 |
| 5.0 | 0.64425 | 0.062551 | 10 | 0.74 | 36010 |

- Per unit time, the risk in the young window is 20–80 times the later
  risk, so transient ringing (explanation 1) is real.
- The later hazard also rises by more than two orders of magnitude
  between $R_0 = 2$ and $3.5$. Patches spend much more time at older
  ages, so most fade-outs happen there. Stochastic amplification
  (explanation 2) therefore contributes most of the total.
- The $R_0 = 1.5$ row rests on very few fade-outs and is mostly noise.

## Deterministic-trough calculation of $P_2$, $P_3$

The `burnout` package (Earn & Bolker) computes $\mathcal{P}_m$
(`Pm_prob()`) only for the SIR model with births and deaths. Our
within-patch model is an SI model with disease-induced death and
logistic host growth, so we use the same method as implemented for that
model in `fadeout/logistic_burnout/logistic_burnout_functions.R`
(`logistic_multitrough_probabilities()`). The method has three steps:

1.  Follow the deterministic trajectory from a single infective in a
    patch at $S = K$.
2.  At the entry to each trough’s boundary layer, compute the Kendall
    extinction probability $q$ of one lineage while susceptibles recover
    logistically.
3.  Raise $q$ to the number of infectives at entry.

$Q_j$ is the probability of burnout in trough $j$, given survival
through trough $j-1$. The table also shows when each trough is entered,
in units of $T_0$ since introduction.

To compare with the simulations, take colonization episodes that are
still infected at $1.5\,T_0$ (constant transmission, all $\alpha$
pooled). The comparable simulated quantity is the proportion of these
episodes that go extinct before $3\,T_0$, a window that contains the
second trough.

| R0 | K | trough 2 entry (T0) | trough 3 entry (T0) | Q1 | Q2 | Q3 | simulated P(extinct 1.5-3 T0) | episodes at risk |
|---:|---:|---:|---:|:---|:---|:---|---:|---:|
| 2.0 | 3000 | 1.698 | 2.745 | 0.18 | 0.00036 | 4e-06 | 0.028 | 5553 |
| 2.5 | 3000 | 1.792 | 2.900 | 0.76 | 0.023 | 0.00017 | 0.088 | 15237 |
| 3.0 | 3000 | 1.946 | 3.141 | 0.97 | 0.25 | 0.0072 | 0.141 | 16844 |
| 3.5 | 3000 | 2.126 | 3.429 | 1 | 0.66 | 0.085 | 0.196 | 13940 |
| 2.0 | 5000 | 1.718 | 2.765 | 0.059 | 2e-06 | 1.1e-09 | 0.005 | 575 |
| 2.5 | 5000 | 1.807 | 2.915 | 0.63 | 0.0018 | 5.1e-07 | 0.037 | 2954 |
| 3.0 | 5000 | 1.957 | 3.153 | 0.95 | 0.098 | 0.00026 | 0.088 | 7090 |
| 3.5 | 5000 | 2.136 | 3.438 | 1 | 0.5 | 0.017 | 0.146 | 13465 |
| 2.0 | 10000 | 1.746 | 2.792 | 0.0034 | 3.8e-12 | 1.2e-18 | 0.000 | 120 |
| 2.5 | 10000 | 1.827 | 2.935 | 0.4 | 3.2e-06 | 2.6e-13 | 0.007 | 1504 |
| 3.0 | 10000 | 1.973 | 3.169 | 0.9 | 0.0097 | 6.8e-08 | 0.036 | 3124 |
| 3.5 | 10000 | 2.149 | 3.451 | 0.99 | 0.25 | 0.00029 | 0.084 | 4501 |

- The deterministic-trough theory reproduces the qualitative pattern:
  $Q_2$ rises by many orders of magnitude with $R_0$. At $K = 5000$, for
  example, it goes from about $10^{-6}$ at $R_0 = 2$ to 0.5 at
  $R_0 = 3.5$. Trough 2 starts at about $1.7$–$2.2\,T_0$, just after the
  establishment threshold, so explanation 1 explains the steep young-age
  hazard.
- Quantitatively the fit is poor in both directions:
  - At low $R_0$ the theory is far too small. It has no noise, so it
    misses fade-outs from explanation 2 and from lineages lingering
    after the first trough.
  - At high $R_0$ it is too large, probably because colonizations in the
    metapopulation often reach patches whose hosts have not yet
    recovered to $K$ after an earlier epidemic. The overshoot is then
    smaller than from $S = K$.
- The “Q” values describe a single deterministic pass and ignore the
  subsequent stationary hazard entirely. They therefore cannot account
  for the 70–80% of fade-outs that happen after $3\,T_0$.

## Summary

The fade-out hazard after establishment rises with $R_0$ because damping
per cycle of the within-patch oscillation weakens as $R_0$ rises: the
amplitude retained per cycle goes from 0.27 at $R_0 = 1.5$ to 0.69 at
$R_0 = 3.5$. Weaker damping does two things:

- the post-invasion transient keeps ringing, giving deep second and
  third troughs (high $Q_2$, $Q_3$);
- demographic noise sustains larger quasi-cycles around equilibrium.

The first effect dominates the hazard per unit time early on. The second
dominates the total number of fade-outs.
