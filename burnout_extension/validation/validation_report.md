# Numerical validation of the general-recruitment burnout asymptotics

## Summary

The refined first-order approximation is accurate through most of the fixed-
`R0`, small-`rho` regime, and the second-order action-resummed approximation is
usually substantially better. On the predefined main validity subset
`1.2 <= R0 <= 6`, `rho <= 0.03` (341 valid-entry points), median relative
`x_in` error falls from 0.0282% at first order to 0.00129% at second order;
95th-percentile error falls from 0.449% to 0.0942%. Median absolute action error
falls from 0.00821 to 0.000430. Second order improves both `x_in` and action at
97.7% of these points.

The results support exact Kendall as the default semi-analytical probability
step. Laplace is excellent for large action but not uniformly accurate at small
action. `NO_ENTRY` is common at larger `rho` and must be classified rather than
assigned a fabricated entry value.

## Methods and reproducibility

The implementation follows the formulas in the unified theory note directly.
`h_theta`, `F`, the final-size root, definite action differences, the regularized
`C_theta` integral, inverse outer coefficients, `D_theta`, both master equations,
the deterministic ODE, exact Kendall integral, and Laplace approximation are in
the R sources documented in [README.md](README.md).

The canonical ODE starts at `y0=1e-10` on the leading unstable direction
`1-x = R0*y/(R0-1+rho)`. Events are detected in order: downward `x=x*`, downward
`y=y_BL`, and upward `x=x*`. A trough preceding the boundary-layer crossing is
reported as `NO_ENTRY`. ODE tolerances are `rtol=2e-10`, with componentwise
absolute tolerances `2e-12, 2e-14`. Selected trajectories were repeated at
`y0=1e-8,1e-10,1e-12`, with a finer output mesh, and with the finite-K-like
initialization `(1-1/K,1/K)`; raw results are in `data/ode_sensitivity.csv`.

An initial pilot grid was used during implementation and then discarded after
the numerical issues it exposed were corrected. The retained dense grid contains 720 points:
five theta values, 16 R0 values from 1.05 through 10, nine rho values from
0.0025 through 0.1, and K=10000. A 108-point K-sensitivity grid uses K=1000,
10000,100000. Retained outputs are `data/dense_results.csv` and
`data/K_sensitivity.csv`.

## Internal consistency checks

Finite-difference tests of the action primitive pass for theta 0, 0.5, and 1.
The action obtained from the second master equation agrees with the action after
root solving to the test tolerance. At theta=0.5, explicit regularized-integral
values of C are 0.07968025, 0.1406151, and 0.1657055 for R0=2,3,5. The independent
inverse-outer values are 0.07968028, 0.1406152, and 0.165706, respectively.
These reproduce the supplied regression checks without hardcoding.

## Stability of D_theta

`D_theta` is estimated by extrapolating the regularized inverse coefficient over
several windows and two remainder bases. The full retained table is
`data/Dtheta_stability.csv`. Extremely small-y windows are demonstrably worse
in double precision because X2 is obtained after severe cancellation. The
production window is geometrically spaced from 2e-3 to 3e-5 and automatically
shrinks below the zero-order peak near R0=1.

For theta=0.5, production estimates are -0.0478, -0.07923, and -0.02567 at
R0=2,3,5, close to the supplied -0.04866, -0.07962, and -0.02574 checks. The
remaining R0=2 window sensitivity is retained as numerical uncertainty rather
than hidden. Near threshold, D is less reliable because the available inverse-
outer window collapses.

The installed R environment lacks an arbitrary-precision package (`Rmpfr`), so
the requested double-versus-high-precision comparison could not be completed in
pure R. Window, basis, point-count, and tolerance sensitivity were completed;
the missing arbitrary-precision cross-check remains the main numerical-quality
limitation of this validation.

## Deterministic matching and action

Of 720 dense-grid points, 513 have valid entries and 204 are `NO_ENTRY`.
Initially, 35 near-threshold points exceeded the integration window; after an
adaptive-horizon rerun, all 35 were physically classified as `NO_ENTRY`. In the
main validity subset,
absolute relative x errors have (median, 90th, 95th, maximum):

- first order: (0.000282, 0.00217, 0.00449, 0.0799);
- second order: (0.0000129, 0.000371, 0.000942, 0.0138).

Absolute action errors have:

- first order: (0.00821, 0.0319, 0.0441, 0.0899);
- second order: (0.000430, 0.00377, 0.00695, 0.0195).

The maxima are concentrated at validity boundaries. The maps in `figures/`
show near-threshold deterioration, high-R0 stress, and the no-entry boundary.

The current error maps use a dedicated 1875-point grid over five rho values,
five K values, five theta values, and 15 values uniform in `log(R0-1)`. Each
multi-page PDF contains one page per fixed rho, with absolute and relative error
shown in separate panels and with independent colour scales. Files cover first-
and second-order `x_in`, first- and second-order action, first- and second-order exact-Kendall B,
and Laplace-only B. The separate `fig08_P1_stochastic_vs_analytic.pdf` shows the five Gillespie estimates and their
confidence intervals against analytical predictions.

Grey cells are primarily physical `NO_ENTRY` classifications, not missing
computations: 1141 of 1875 points are `NO_ENTRY`, 727 are valid, and only seven
have an ODE entry but no physical second-order master-equation root. Because
`y_BL` is proportional to `K^(-1/2)`, increasing K lowers the matching line and
can increase the no-entry region. Light grey cells in every retained map show
this classification directly, so a separate classification figure is omitted.

## Burnout exponent and Laplace decomposition

Using exact Kendall, the median relative B error in the main validity subset is
0.819% at first order and 0.0432% at second order; the 95th percentiles are 4.54%
and 0.821%. Thus most probability error removed by second order is traceable to
deterministic matching/action error.

With the ODE entry held exact, the Laplace-only median relative B error is 0.108%
in the same subset, but its 95th percentile is 7.17% and maximum 29.1%. Among
81 valid points with Lambda<2, the median Laplace-only error is 9.25%. This
confirms that a bad fully asymptotic result at small action must not be blamed on
`x_in`; exact Kendall should be retained there.

## Internal truncation diagnostic

Spearman correlation between `|Lambda2-Lambda1|` and the actual second-order
action error is 0.572 over valid dense-grid points. The order gap is therefore a
useful qualitative warning but not a universal error estimator. It successfully
identifies many deteriorating regions, but the scatter is too wide to justify a
single hard threshold.

## Stochastic check

The exact count CTMC uses infection, removal, and recruitment rates specified in
the task. Simulations condition on the first downward threshold crossing. A
sqrt(K)-count hysteresis is required before arming the upward recrossing event;
without it, a one-step discrete bounce at an integer threshold is falsely counted
as persistence. Each selected point uses 600 established waves and Wilson 95%
intervals. Results are in `data/stochastic_results.csv`. At
`(R0,rho,theta,K)=(3,0.03,0.5,10000)`, simulation gives 1 persistence in 600
established waves (`P1_hat=0.00167`, interval 0.000294--0.00938), versus
0.000611 from second-order + exact Kendall. At `(5,0.1,0.5,10000)`, simulation
gives 0.9967 (0.9879--0.9991), versus 0.9813. Three low-probability checks saw
zero persistence; their 95% upper bound is 0.00636, so they are compatible with
the tiny analytical predictions but do not resolve them quantitatively.

These finite-K results should be interpreted as a check on the complete
deterministic-entry branching approximation, not as another test of the ODE
matching formula. Discrepancies incorporate random entry location and the
discrete event definition. They are explicitly kept separate from the matching
and Laplace error decompositions above.

## Anomalies and iterations

Three implementation issues were found and corrected during actual runs:

1. Uniform x-grid interpolation destroyed the inverse coefficient near x_f;
   coefficient ODEs are now evaluated exactly at requested inverse roots.
2. Overly small y windows made D appear divergent through floating-point
   cancellation; multiple windows exposed this and motivated the production
   window and stability flagging.
3. Near R0=1 the default y window exceeded the zero-order epidemic peak, and at
   R0=10 the original final-size bracket excluded x_f. Both domains are now
   selected adaptively.
4. The initial discrete stochastic recrossing logic double-counted a downward
   landing on an integer threshold; side-aware crossings plus hysteresis fix it.

No theory formula was altered to improve agreement.

## Scientific recommendation

For fixed R0 away from 1 and small rho, use second-order `x_in` plus the exact
one-dimensional Kendall integral. Refined first order is already strong when a
faster estimate is needed. Use Laplace only after checking that action is large;
the data do not support a universal cutoff, although Lambda below about 2 is a
clear warning region here. Always classify `NO_ENTRY`. Treat large rho,
near-threshold R0, and the high-R0 distinguished regime as outside or near the
edge of the present expansion, and propagate D-window instability as a warning.
