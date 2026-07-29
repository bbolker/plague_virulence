# Logistic epidemic burnout

This directory is a self-contained exploratory extension of the epidemic
burnout approximation of [Parsons et al. (2024)](https://doi.org/10.1073/pnas.2313708120)
and the associated [`davidearn/burnout`](https://github.com/davidearn/burnout)
code. It applies the deterministic-to-stochastic decomposition to the
single-patch model used in this repository, in which susceptible hosts recover
logistically rather than through linear vital dynamics.

The primary quantity is a **conditional post-epidemic burnout probability**.
It is not the unconditional probability that one introduced infection
eventually becomes extinct.

## 1. Relationship to the Parsons et al. method

The retained structure is

$$
\text{deterministic major epidemic}
\longrightarrow
\text{entry into a low-infection boundary layer}
\longrightarrow
\text{time-dependent birth-death process}.
$$

Outside the boundary layer, a large epidemic is approximated by a deterministic
trajectory. Inside the boundary layer, infected-host counts are small enough
that demographic stochasticity matters, so the infection process is replaced
by a time-inhomogeneous linear birth-death process.

Parsons et al. use an SIR model with linear host vital dynamics. Their
susceptible trajectory inside the boundary layer consequently has a linear
recovery form, and their paper develops further asymptotic approximations.
This module retains the boundary-layer and Kendall birth-death argument but
changes the susceptible environment to the explicit logistic recovery
appropriate to the present repository. The numerical deterministic entry,
logistic recovery formula, and numerical outer integral are therefore a new
adaptation; code from `davidearn/burnout` is not copied.

The paper distinguishes:

1. **fizzle**: early stochastic extinction following introduction;
2. a deterministic **major epidemic** conditional on escaping fizzle;
3. **burnout** during the first post-epidemic low-infection phase;
4. longer persistence near the endemic state.

This module estimates item 3. Starting the deterministic trajectory at
$I(0)=1$ represents the major-outbreak path from that initial state. It cannot
represent the stochastic chance that this infection disappears before the
major epidemic. For one initially infected host and $x(0)\simeq1$, the separate
early-extinction branching approximation is $1/R_0$; it is not combined with
the burnout probability here.

## 2. Normalized logistic model

In host counts,

$$
\frac{dS}{dt}
=rS\left(1-\frac{S}{K}\right)-R_0\frac{SI}{K},
\qquad
\frac{dI}{dt}
=R_0\frac{SI}{K}-I.
$$

Time is measured in disease generations, so the infectious removal rate is
one. Define

$$
x=\frac{S}{K},\qquad y=\frac{I}{K}.
$$

Dividing the count equations by $K$ gives

$$
\frac{dx}{dt}=rx(1-x)-R_0xy,
\qquad
\frac{dy}{dt}=(R_0x-1)y.
$$

At an endemic equilibrium with $y^*>0$, the second equation gives

$$
x^*=\frac{1}{R_0}.
$$

Substitution into the first gives

$$
y^*=\frac{r(R_0-1)}{R_0^2}.
$$

Therefore,

$$
S^*=\frac{K}{R_0},
\qquad
I^*=Ky^*=\frac{rK(R_0-1)}{R_0^2}.
$$

## 3. Boundary-layer definition

The boundary scale is the logistic endemic infected density,

$$
y_{\mathrm{BL}}=y^*
=\frac{r(R_0-1)}{R_0^2}.
$$

This follows Parsons et al.'s choice to define the boundary by the endemic
infection scale. In the logistic model it remains an exploratory analogue,
not a rigorous matched-asymptotic derivation equivalent to theirs.

The full deterministic system is solved from

$$
S(0)=K-I_0,\qquad I(0)=I_0.
$$

Because

$$
\frac{dI}{dt}=(R_0x-1)I,
$$

the first downward crossing of $x=1/R_0$ locates the first epidemic peak.
After that peak, the first downward crossing of
$y=y_{\mathrm{BL}}$ defines

$$
t_{\mathrm{in}},\qquad
x_{\mathrm{in}}=x(t_{\mathrm{in}}),\qquad
y_{\mathrm{in}}=y(t_{\mathrm{in}}).
$$

Both crossings are linearly interpolated between ODE output times. A valid
entry requires

$$
t_{\mathrm{in}}>t_{\mathrm{peak}},\qquad
\left.\frac{dy}{dt}\right|_{t_{\mathrm{in}}}<0,\qquad
y_{\mathrm{in}}\simeq y_{\mathrm{BL}}.
$$

The solver horizon begins at 50 disease generations and doubles adaptively to
a maximum of 1600. A missing peak or boundary crossing returns an explicit
failure status; the first trough is not substituted.

## 4. Logistic susceptible recovery inside the boundary layer

Reset boundary-layer entry to local time zero. When infection density is low,
neglect its feedback on susceptible dynamics:

$$
R_0xy\simeq0,
\qquad
\frac{dx}{dt}=rx(1-x),
\qquad
x(0)=x_{\mathrm{in}}.
$$

Separating variables,

$$
\frac{dx}{x(1-x)}=r\,dt.
$$

Since

$$
\frac{1}{x(1-x)}=\frac{1}{x}+\frac{1}{1-x},
$$

integration gives

$$
\log\left(\frac{x}{1-x}\right)=rt+C.
$$

Using $x(0)=x_{\mathrm{in}}$,

$$
x(t)=
\frac{x_{\mathrm{in}}e^{rt}}
{1-x_{\mathrm{in}}+x_{\mathrm{in}}e^{rt}}
=
\frac{1}
{1+\left(\frac{1-x_{\mathrm{in}}}{x_{\mathrm{in}}}\right)e^{-rt}}.
$$

Define

$$
A=\frac{1-x_{\mathrm{in}}}{x_{\mathrm{in}}},
$$

so that $x(t)=[1+Ae^{-rt}]^{-1}$. Since $r>0$, $x(t)\to1$.

## 5. Birth and death rates for an infected lineage

The total infection-event rate in the count model is

$$
R_0\frac{SI}{K}=R_0xI.
$$

Conditional on the deterministic environment $x(t)$, each infected individual
therefore creates new infections at per-capita rate

$$
\lambda(t)=R_0x(t).
$$

Here “birth” means creation of a new infection, not demographic birth of a
host. The normalized per-infected removal rate is

$$
\delta(t)=1.
$$

The boundary-layer approximation is thus

$$
I\longrightarrow I+1
\quad\text{at rate }R_0x(t)I,
$$

$$
I\longrightarrow I-1
\quad\text{at rate }I.
$$

Its mean satisfies

$$
\frac{d\,\mathbb E[I]}{dt}
=[R_0x(t)-1]\mathbb E[I],
$$

which agrees with the low-density deterministic infected equation.

## 6. Single-lineage extinction probability

For a time-inhomogeneous linear birth-death process beginning with one infected
individual at boundary-layer time zero, Kendall's representation gives

$$
J=
\int_0^\infty
\delta(t)
\exp\left[
-\int_0^t\{\lambda(s)-\delta(s)\}\,ds
\right]dt.
$$

With $\delta(t)=1$ and $\lambda(t)=R_0x(t)$,

$$
J=
\int_0^\infty
\exp\left[
-\int_0^t\{R_0x(s)-1\}\,ds
\right]dt,
$$

and the eventual extinction probability of one lineage is

$$
q_1=\frac{J}{1+J}.
$$

For a constant supercritical rate $\lambda>\delta$, the integral is
$J=\delta/(\lambda-\delta)$ and hence $q_1=\delta/\lambda$, providing a
consistency check.

## 7. Analytical inner integral

Using $x(t)=[1+Ae^{-rt}]^{-1}$,

$$
\int_0^t x(s)\,ds
=\frac{1}{r}
\log\left(\frac{e^{rt}+A}{1+A}\right).
$$

Therefore,

$$
H(t)
=\int_0^t[R_0x(s)-1]\,ds
=\frac{R_0}{r}
\log\left(\frac{e^{rt}+A}{1+A}\right)-t,
$$

and

$$
J=\int_0^\infty e^{-H(t)}\,dt
=
\int_0^\infty
\exp\left[
t-\frac{R_0}{r}
\log\left(\frac{e^{rt}+A}{1+A}\right)
\right]dt.
$$

The inner integral is analytical; the outer integral is numerical. An
incomplete-beta transformation is deliberately not used.

## 8. Numerical integration

Direct evaluation of $e^{rt}$ can overflow. The implementation instead
computes

$$
\log(e^{rt}+A)
$$

with a log-space addition,

$$
\operatorname{logadd}(a,b)
=m+\log[e^{a-m}+e^{b-m}],
\qquad
m=\max(a,b),
$$

where $a=rt$ and $b=\log A$. The cases $A=0$ ($x_{\mathrm{in}}=1$),
$x_{\mathrm{in}}$ close to zero or one, small $r$, and $R_0$ close to one
remain on this stable scale.

`stats::integrate()` evaluates $J$ with an infinite upper limit. The function
returns its absolute error, subdivisions, message, any recorded warnings, and
a convergence flag. The stable transformation

$$
q_1=\operatorname{logit}^{-1}(\log J)
$$

avoids explicitly forming $J/(1+J)$ when $J$ is large. Failures remain missing;
they are never replaced by zero.

## 9. Multiple lineages at boundary entry

The deterministic boundary count is

$$
m_{\mathrm{raw}}=Ky_{\mathrm{BL}}.
$$

If the lineages are treated as independent, an integer count $m$ gives

$$
P_{\mathrm{post\ burnout}}=q_1^m,
\qquad
\log P_{\mathrm{post\ burnout}}=m\log q_1.
$$

The default is

$$
m=\max[1,\operatorname{round}(m_{\mathrm{raw}})].
$$

`floor`, `ceiling`, and a continuous diagnostic
$m=m_{\mathrm{raw}}$ are also implemented. The continuous option computes
$m_{\mathrm{raw}}\log q_1$ but does not claim that a fractional number of
independent lineages is literal. When $m_{\mathrm{raw}}<1$, integer methods use
one lineage and return a status that records this adjustment.

## 10. Role of the initial infected host

The main surface uses $I_0=1$ in the full deterministic trajectory. This is the
infection count at epidemic introduction. It is not the number of lineages at
post-epidemic boundary entry, which is approximately

$$
Ky_{\mathrm{BL}}.
$$

The deterministic $I_0=1$ trajectory is conditioned conceptually on following
the major-outbreak path. Early stochastic extinction is absent from the ODE and
from the plotted probability.

## 11. Implementation and validation

| File | Purpose |
|---|---|
| `logistic_burnout_functions.R` | Deterministic trajectory, interpolated boundary entry, logistic recovery, stable $H(t)$, Kendall integral, burnout wrapper, and exact thinning simulator. |
| `validate_logistic_burnout.R` | Formula checks, integration convergence checks, and direct nonhomogeneous branching-process validation. |
| `plot_burnout_surface.R` | Generates the $R_0$-$r$ burnout surface and boundary-entry diagnostic. |
| `outputs/validation_results.csv` | Detailed results of all validation checks. |
| `outputs/logistic_burnout_R0_r_grid.csv` | Complete 41 by 41 parameter grid, including failures and integration diagnostics. |
| `outputs/logistic_burnout_status_summary.csv` | Counts of successful and failed grid cells. |
| `figures/validation_q1.png` | Analytical $q_1$ versus thinning simulation. |
| `figures/logistic_burnout_R0_r_heatmap.{png,pdf}` | Main post-epidemic burnout surface. |
| `figures/logistic_boundary_xin_heatmap.{png,pdf}` | Diagnostic surface for $x_{\mathrm{in}}$. |

The validation suite compares the explicit logistic recovery with numerical
ODE solutions, compares $H(t)$ with direct quadrature, tightens the improper
integral tolerances, and uses exact thinning for the nonhomogeneous branching
process. Thinning proposes events from the bound
$(R_0+1)I$, then accepts events according to the current $x(t)$; it does not use
a homogeneous Gillespie approximation. Surviving lineages are followed until
$I=300$ or 500 disease generations. This validates the boundary-layer
birth-death calculation, not the full stochastic logistic S-I approximation.

## 12. Current exploratory surface

The current grid fixes

$$
K=10{,}000,\qquad I_0=1,
$$

and uses 41 evenly spaced values in each range

$$
1.05\le R_0\le5,\qquad 0.01\le r\le0.5.
$$

$K=10{,}000$ is the common baseline in the existing `fadeout` comparisons.
The plotted fill is exactly

$$
P_{\mathrm{post\ burnout}}
=q_1^{\max[1,\operatorname{round}(Ky_{\mathrm{BL}})]}.
$$

It excludes early extinction. Grey cells indicate no valid deterministic
first peak or numerical failure. In the current grid, 1613 of 1681 cells
succeed. The 68 failures are `no_epidemic_peak`: all 41 cells at $R_0=1.05$
and the 27 cells with $R_0=1.14875$ and $r\ge0.1815$. These trajectories do
not exhibit the required first overshoot and therefore do not define the
requested post-epidemic downward boundary entry.

The surface should be treated as a documented exploratory logistic extension.
Agreement with the boundary-layer branching simulation does not establish that
the boundary scale $y_{\mathrm{BL}}=y^*$ is asymptotically correct or that the
approximation agrees with full stochastic logistic S-I simulations.

## Reproducibility

Required R packages are `deSolve` and `ggplot2`.

From the repository root:

```bash
Rscript fadeout/logistic_burnout/validate_logistic_burnout.R
Rscript fadeout/logistic_burnout/plot_burnout_surface.R
```

To load only the functions:

```r
source("fadeout/logistic_burnout/logistic_burnout_functions.R")
```

On the development machine, validation takes approximately 20 seconds and the
41 by 41 surface approximately 20 seconds. A typical laptop should finish each
script within a few minutes. Generated files are confined to `outputs/` and
`figures/` in this directory.
