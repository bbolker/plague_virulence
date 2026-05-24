## To do

* makestuff/shellpipes?
* summary info:
    * n occupied, quasi-eq, etc.
    * YZ PIP invasion characteristics
* single-strain runs	
* run PIPs
* 1 30x30 PIP array (200 sims per parameter set) should take about 10 hours in serial -- 
* job array/META runs on Nibi?
* profiling macpan2 runner
* set up tinytest machinery for basic tests
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
