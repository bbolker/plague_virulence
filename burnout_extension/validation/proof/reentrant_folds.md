# Three certified local folds at one fixed theta,K

Fix theta=11/20 and K=exp(15/4)>1. Retain the exact amplitude integral C_theta
and the exact action integral A in the unconditional BI expression; no
large-R or small-B replacement is used in these certificates.

Let lambda=1/rho, t=log(1-1/R), q=R(R-1),
B=K C sqrt[(R-1)R^(-theta)/(2 pi lambda)] exp(-lambda A),
s=U+lambda V, H(B)=(exp B-1)/B, D(B)=B H'(B).
The R derivative has the sign of E=H(B)-s, so a fold is the system
F(R,lambda)=(E,E_t)=0. On the whole certification box,

E_t=-sD-s_t,
E_lambda=-(A+1/(2lambda))D-V,
E_tt=-s_t D+s^2 B D'-s_tt,
E_tlambda=-V D+s(A+1/(2lambda)) B D'-V1.

The Jacobian in (R,lambda) coordinates is
[[E_t/q,E_lambda],[E_tt/q,E_tlambda]].

The script `reentrant_folds.py` encloses the amplitude and all its needed
derivatives using validated Arb contour quadrature and parameter Taylor
bounds. The action and its derivatives are also enclosed. With a fixed
nonsingular approximate inverse C, it evaluates the Krawczyk operator

c-C F(c)+(I-C J(X))(X-c).

For each of the three disjoint boxes X, the operator is strictly inside X.
The weighted row-sum bounds of I-CJ(X) are below 0.001, so the map is also a
strict contraction on X. Thus each box contains exactly one zero of F.
The certificates are stored in `reentrant_fold_certificates.json`.

| R near | certified rho enclosure | sign of -E_lambda E_tt | local change as rho increases |
|---|---|---|---|
| 6.791443 | 0.3111753521 +/- 1.88e-11 | positive, 35.60 +/- .00392 | pair annihilation |
| 7.806470 | 0.33337104603 +/- 9.19e-12 | negative, -37.05 +/- .00713 | pair creation |
| 9.555305 | 0.33716252287 +/- 4.33e-12 | positive, 186.0 +/- .0521 | pair annihilation |

The orientation follows because E_rho=-lambda^2 E_lambda and, at a fold,
E_RR=E_tt/q^2. A positive product E_rho E_RR means two nearby stationary
points exist for rho below the fold and disappear above it. A negative
product gives the reverse. Both factors are nonzero in every box, so these
are ordinary local folds, not degenerate higher-order singularities.

This proves three distinct local bifurcations for the same theta,K. It does
not prove that these are the only folds or establish a global number of R
extrema. A fine numerical scan suggests total counts 3,1,3,1 across the four
rho intervals, but those total counts are not certified by the local boxes.
In particular, no global upper bound has been inferred from these examples.

The statement concerns the BI expression itself. The rho values are not in
the small-recruitment limit, so these local bifurcations are not asserted to
be features of the full stochastic epidemic model.
