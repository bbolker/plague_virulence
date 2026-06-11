## To do

* update euler examples; thin properly, summarize before storing
* discuss trough-run output
    * interpolate between leaky-bucket and logistic equation?
    * `r*S*(1-S/K)` vs `r*(K-S)` 
    * (**not** theta-logistic)
    * non-dimensionally, this is `x*(1-x)` vs `(1-x)`; dimensionally, `r*K*(S/K)*(1-S/K)` vs `r*K*(1-S/K)`. So `r*K*(S/K)^theta*(1-S/K)` recovers leaky-bucket as `theta → 0` and logistic as `theta → 1` (need to be a little careful implementing this with the special case of `S=0, theta = 0` (`0^0 = 1` in IEEE 754, which is actually what we want in this case ...)
* seasonality!  (see @krauerInfluenceTemperatureSeasonality2021)
* keep working on single-patch-introduction PIP runs
* direct solution of deterministic metapop model??? does it collapse?

* implement non-extinction variant, compare with extinction-allowed (strength of burnout mechanism)
* implement power-law R0(N) ?
* think about an invasibility expression: if we know patch occupancy and S/I distributions of resident, can we calculate invasion, or do extinction dynamics etc. make expectations of instantaneous growth rate irrelevant?
* sensitivity to $r$?
* test ODE with euler integrator (pass through to deSolve)? allow switch to disable hazard correction in stoch models?

### code

* refactor two-strain delay to be part of discrete_run (chunked run), rather than part of model definition (helps with trough calculations as well)
* refactor stop_condition to allow calculations based on more than last row (current state) alone: parameters (e.g. -> equilibrium values), derivative, complete trajectory? (for trough runs, want to know if we're beyond `t_Seq`)
* lightweight S3 code for odin runs (print, plot methods ...); store parameter values as attributes. Make sure this propagates e.g. for multi-sample runs to summary?
* prototype and benchmark hybrid deterministic (above threshold)/stochastic (below threshold) model
* create dedicated single-strain model (doesn't help that much)
* instrument runs to record r(t), foi(t), etc. (to understand invasion criteria etc.)
* makestuff/shellpipes?
   * sort out working-directory stuff (package??)
   * interference between shellpipes/makestuff and command-line args?
* profiling macpan2 runner
* lay out/test mechanisms: dens-indep R0, fixed N, extinction-independent ... ?
* better stochastic (+/-) dynamics (cleaner, avoid clamping in macpan2)

Note effects of [1] generation interval/infectious period distribution; [2] demography (SIR w/ constant (leaky-bucket) vital dynamics vs SIRS vs SID w/ logistic growth) [3] ??

  Chose r=0.125 on the basis of intrinsic rat pop growth rates of 3-5/year; generation time of plague 10-20 days; `3/365*15 ~ 0.125`.

  Linear growth rate would be equal to death rate, in disease generation times (== epsilon); death rate is about 1 (lifespan == 1 year), so approx 1/3 of r == 0.04 (maybe less, 15 is on the high end)

## old

* tweak poisson/sim.R; allow sim prefix
* allow burnout in poisson sim?
* check Levine and Earn for two-strain final info
* run odin two-strain for burnout calcs? (with vital dynamics, with some rubric for finding first trough)
* compare with YZ computations?

* Go back and compare analytical results ...
* miscellaneous thoughts
    * flea travel/colonization decoupled from rat movement?
    * estimate coupling (somehow) from plague movement? (Fisher equation etc)
    * how strongly can we assume that transmission is decoupled from infectious period/rat survival?
    * look over Keeling and Gilligan
    * are we off-base in not considering the sylvatic component? https://journals-asm-org.libaccess.lib.mcmaster.ca/cms/10.1128/aem.01658-25/asset/48d7afed-77ed-4bea-9292-91ec846f371d/assets/images/large/aem.01658-25.f001.jpg
