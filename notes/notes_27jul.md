## notes from JD/BB/YZ meeting

* two-state (persistently infected, transiently infected patches, uninfected patches) approximations are working reasonably well, across a one-parameter-at-a-time sweep around the baseline parameters (i.e. varying $R_0$, $K$, $r$, $\alpha$ one at a time)
* the estimates of average infectious pressure and lifetime of transiently infected patches are a bit *ad hoc*; scaling lifetime with $T_{\textrm{osc}}$ and computing average (cross-patch) force of infection from an ODE simulation
* from existing one-strain-one-patch simulations, make plots exploring the actual relationship between (log) time-until-burnout and (log) $T_{\textrm{osc}}$; is it approximately constant? What's the range of ratios? (Comparing on a log scale is a good idea since it will automatically give ratios and may give hints if there is a power-law scaling.) (Draw plots with facets/colours/point shapes etc. reflecting different values of parameters)
* use final size instead of $\int I \, dt$ to estimate average
* maybe use David Earn's new SIR approximations to compute an approximate burnout time? (i.e. $t$ such that $I(t) = 1$ after the first peak? (I have a draft ms.; I will ask if we can share it with YZ)
* go back to the Parsons et al. burnout paper and see what we would have to do to re-do the derivation for logistic dynamics (can a LLM do this ???)
* start writing an outline of what we've discovered/would like to say so far? (BMB)
* consider when/whether the transient phase is so long that we might get a plausible case for persistence even without adding other dynamical details (such as seasonality)?
* continue working on the seasonal version of the model
