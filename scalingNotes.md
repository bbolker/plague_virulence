The form we should use is

$$S' = \rho(K-S)\left(\frac{S}{K}\right)^\theta.$$

Thus, when $\theta=0$,

$$S' = \rho(K-S),$$

and when $\theta=1$,

$$S' = \rho S(K-S)/K,$$

so $\rho$ matches both $\varepsilon$ and $r$

There is an infelicity here, which I think we should ignore for now.

In the linear case $\theta=0$, we have explicit deaths, whereas for $\theta=1$, we explicitly have _no_ deaths. In between, we're kind of phenomenological (it is not clear how best to parse the S' given above into births and deaths). This seems fine for our purposes but is worth noting.

For the record, there is a pretty nice form with explicit accounting of mortality, but it requires separate parameters for $r$ and $\mu$. The $\mu$ is kind of a nuisance parameter -- in particular, it's the one that disappears when $\theta=0$
$$S' = r(K-(r-\mu)S/r)\left(\frac{S}{K}\right)^\theta - \mu x.$$

Still slightly weird (our model is the limit $\mu\to 0$, which is cool but the classic linear case has the interpretations that births go to zero at equilibrium (rather than that deaths match births).
