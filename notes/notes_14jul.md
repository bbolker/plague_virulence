## meeting notes

* YZ, DJDE, JD, BMB

Becoming clear that the basic dynamics of patch occupancy in a single-strain metapopulation, *if* we start all patches in an "invasion" scenario [i.e. $I(0) = I_0, S(0) = K-I_0, I_0 \ll K$], is (1) an initial burnout of most patches followed by (2) a gradual increase in the patch occupancy (because patches that survive burnout will persist more or less indefinitely). In addition to this overall trend, there is a background of stochastic colonization and burnout.

We should be able to characterize a lot of this; for example, data for the persistence probability (and other summary statistics) of a single patch under logistic demography, with continuous dynamics, as a function of $R_0$ and $K$. The data are in:

```
odin/sharcnet/outputs/euler_onepatch_onestrain_extinct_logistic_continuous.rds
```

And this draws the picture:

```bash
Rscript odin/euler_onepatch_onestrain_extinct_plot.R --combo logistic_continuous 
```

This is not a simple relationship, but we use a lookup table or emulate it with a spline/Gaussian process/etc. fit. (Might be simpler to fit if we look only at $R_0 > 2$, which is the non-fizzle regime ...)

This should give us $P_1(R_0, K)$. Assuming that patches that persist rapidly get close to their deterministic equilibrium, we should be able to compute colonization probability etc. etc.. and get a pretty good characterization of the slope of the patch-occupancy curve.

Adding seasonality will make this a whole different ball game.
