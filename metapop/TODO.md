## To do

* switch from job array to META runs?
* retrieve batch2 onestrain results, plot
* run two-strain mini example
* set up plotting (à la YZ) for the two-strain PIP plots
* compare deterministic models with leaky-bucket vs logistic demography (trough depth vs R0)
* re-run one-strain/one-patch/extinct with more simulations
  
* create a model that has sensible/expected burnout behaviour for one patch, one strain (note effects of [1] generation interval/infectious period distribution; [2] demography (SIR w/ constant (leaky-bucket) vital dynamics vs SIRS vs SID w/ logistic growth) [3] ??)

  Chose r=0.125 on the basis of intrinsic rat pop growth rates of 3-5/year; generation time of plague 10-20 days; `3/365*15 ~ 0.125`. However, this leads to near-certain burnout for R0 = 4 (with dt = 0.2 disease generations) at K=1e6. 0.125 is *way* above the range shown in Parsons et al (which goes only up to 0.02, reasonable for human diseases).  What do the burnout contours look like for epsilon in the range of 0.02 to 0.2 ... ?

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
* implement density-dependent R0?
* re-run baseline poisson sim example?
* re-run PIP sims, with ODE final sizes - on SHARCNET?
* check Levine and Earn for two-strain final info
* run odin two-strain for burnout calcs? (with vital dynamics, with some rubric for finding first trough)
* compare with YZ computations?

* investigate details of runs in 'blackout' range of new middle PIP
* Go back and compare analytical results ...
* miscellaneous thoughts
    * flea travel/colonization decoupled from rat movement?
    * estimate coupling (somehow) from plague movement? (Fisher equation etc)
    * how strongly can we assume that transmission is decoupled from infectious period/rat survival?
    * look over Keeling and Gilligan
    * are we off-base in not considering the sylvatic component? https://journals-asm-org.libaccess.lib.mcmaster.ca/cms/10.1128/aem.01658-25/asset/48d7afed-77ed-4bea-9292-91ec846f371d/assets/images/large/aem.01658-25.f001.jpg
