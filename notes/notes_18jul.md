## notes on analytic approximation

* get $P_1(R0, K)$ [persistence through first trough] from GAM fit (tensor product smooth) to grid of one-patch, one-strain stochastic sim extinction probabilities
   * could in principle get it by repeating the Parsons et al perturbation analysis for logistic demography?
* assume $P_i \approx 1$ for $i > 1$ (i.e., permanent persistence conditional on non-fizzle/non-burnout in the first trough)
* assume that patches that persist rapidly reach deterministic $I^*$ (joint demog/epi equilibrium)
* assume that extinct patches rapidly reach $S=K$ (disease-free eq.)
* calculate colonization probability $c(t)$ from patch occupancy, $I^*$ , $\alpha$
* calculate *successful* colonization prob as $c(t) \cdot P_1(R0,K)$
* approximate burnout as instantaneous?

YZ's results from `fadeout/run_occupancy_exploration.R`. Baseline parameters

```{r}
R0 = 2.5, K = 1e4, r = 0.125, alpha = 1e-4, gamma=1, n_patch = 200, I_outbreak = 10
```

(plus `dt = 0.1 t_max = 2000`)

results from single-parameter (one-at-a-time) sensitivity analysis

* bad as soon as $R_0$ goes from the baseline 2.5 value up to 3 (approximation strongly underpredicts)
* approximation is pretty good for varying $\alpha$ (starts to diverge for $\alpha=1e-5$)
* terrible as $r$ deviates from the baseline value in either direction (variation in the right direction, but much more sensitive)
* insensitive to $K$

Two possible mechanisms:

* not accounting for contribution of transiently occupied patches to colonization 
    * extend approximation based on expected time until first trough/burnout, which should give the expected number of transiently occupied patches at any instant? What should we assume about $I$ in transiently infected patches? is $I^*$ a sufficiently good first guess?
    * use results from sims to compute proportion of infections coming from transiently vs permanently infected pops (i.e. $\sum I_j$ for each)?
* not accounting for variation in initial conditions (i.e. $I_0$ is currently assumed to equal 10 for all invasions)
    * extend fitted $P_1(R_0, K)$ to $P_1(R_0, K, I(0))$ ?
    * force simulation to use $I(0)=10$ always
