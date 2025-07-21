## notes on metapop model etc.

* Hendrik et al. think the generation interval might *not* be longer for PLA-depleted strains? (that would make the "fast-slow" hypotheses irrelevant ...)
* what do we need to do to set up a metapopulation model with burnout?
    * burnout equation depends on (N, R0, epsilon)
	* my first thought was to make N exogenous, make per-patch colonization rates identical. Then strain with lower extinction rate will win. But this doesn't quite work, doesn't give any colonization advantage for higher R0?
	* also doesn't take account of population feedbacks. Allow N (rat population) to be dynamic, (logistic growth - final size of epidemic)?
	* in discrete time, not clear what a time step should be?

## further (April 16)

* suppose a patch-occupancy model with coinfection
* within-patch exclusion on the basis of $\rzero$
* patch extinction based on burnout probabilities
* patch recovery based on logistic rat pop growth?


