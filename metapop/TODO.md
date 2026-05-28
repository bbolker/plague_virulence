## To do

* trouble-shoot PIP runs
   * retrieve two-strain mini example.
   * strain 2 is always going extinct in the mini-example. Is that expected?

* switch from job array to META runs? (ugh, nasty)
   * switch back to job arrays?
* odin ODE version; quantify trough depth vs R0 for leaky-bucket vs logistic demography
* compare deterministic models with leaky-bucket vs logistic demography (trough depth vs R0)
* re-run one-strain/one-patch/extinct with more simulations
* implement non-extinction variant, compare with extinction-allowed (strength of burnout mechanism)
* implement power-law R0(N) ?

* create dedicated single-strain model (doens't help that much)

Note effects of [1] generation interval/infectious period distribution; [2] demography (SIR w/ constant (leaky-bucket) vital dynamics vs SIRS vs SID w/ logistic growth) [3] ??

  Chose r=0.125 on the basis of intrinsic rat pop growth rates of 3-5/year; generation time of plague 10-20 days; `3/365*15 ~ 0.125`.

* instrument runs to record r(t), foi(t), etc. (to understand invasion criteria etc.)
* store parameter values as attributes
* makestuff/shellpipes?
   * sort out working-directory stuff (package??)
   * interference between shellpipes/makestuff and command-line args?
* summary info:
    * n occupied, quasi-eq, etc.
    * YZ PIP invasion characteristics
* profiling macpan2 runner
* lay out/test mechanisms: dens-indep R0, fixed N, extinction-independent ... ?
* invasibility equations
* better stochastic (+/-) dynamics (cleaner, avoid clamping in macpan2)

## old

* update CSEE slides?
* tweak poisson/sim.R; allow sim prefix
* allow burnout in poisson sim?
* re-run PIP sims, with ODE final sizes - on SHARCNET?
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
