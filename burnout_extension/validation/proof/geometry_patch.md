# Continuous-parameter Taylor certificates for the full geometry

These routines retain the exact amplitude and action integrals. They are not
large-R replacements of the model. A patch covers a real rectangle in
S=-log(x_f) and theta, with a separately validated holomorphic polydisk bound.

## Analytic enclosure

The amplitude contour identity is differentiated in the S coordinate. A
finite Arb enclosure on a complex parameter square encloses the corresponding
polydisk. The integrand denominators and logarithm branches are checked by the
validated quadrature callback. The fixed s contour contains 0 and S for every
parameter in the polydisk. It has no other enclosed F zeros at the real center;
none can enter while its boundary remains free of zeros. This supplies the
holomorphic continuation used by Cauchy's coefficient estimate.

The second contour kernel evaluation uses log(s)-log(s-S) on the horizontal
and right edges, and log(-s)-log(S-s) on the left edge. The arguments remain
in common half-planes, so these differences agree with log[s/(s-S)] without
crossing its branch. The contour height is pi. In the far-right part of the
contour, factoring F=-u[1-(f+s/R)/u] avoids interval cancellation.

Complex exponentials are bounded componentwise using monotonic real
exponential bounds, endpoint/critical-point trigonometric bounds, and vertex
bounds for real interval multiplication. This is still outward-rounded Arb
arithmetic; it reduces excess width that otherwise obscures nonzero
integrand denominators. The original and rewritten contour evaluations have
been checked to overlap through fourth derivatives at S=1,3,9.

`complex_geometry.py` bounds the analytic functions U,V and their needed t
derivatives on the polydisk. Its action integral uses the exact identity

A_R=exp(-delta S)-S int_{log(R)/S}^1
      exp[-(delta+1)S y]/[1-exp(-S y)] dy.

The logarithmic x coordinate avoids an exponentially small moving endpoint.
The callback verifies that its denominator stays away from zero.

## Taylor coefficients and remainders

`bivariate_jet.py` implements formal bivariate Taylor arithmetic with Arb
coefficients. The amplitude's central mixed derivatives come from validated
contour quadrature. Central theta derivatives of A_R are direct integrals;
S derivatives follow from the exact formula for A_RR. Full geometry functions
are composed before interval evaluation, retaining the algebraic
cancellations between derivatives.

If |Z|<=M_Z on polydisk radii r_S,r_theta, and the target rectangle satisfies
kappa=max(h_S/r_S,h_theta/r_theta)<1, every coefficient obeys
|z_ij|<=M_Z/(r_S^i r_theta^j). A total-degree N Taylor polynomial therefore
has remainder at most

M_Z (N+2) kappa^(N+1)
 /[1-kappa (N+3)/(N+2)],

provided the denominator is positive. This is an analytic Cauchy remainder,
not an error estimated by comparing evaluations at points.

## Certified patch

The degree-14 certificate covers
S in [8.95,9.05], theta in [.545,.555],
using Cauchy radii (.5,.1), with the full exact theory. It proves

V1<0, qa=V1 V3-V2^2<0,
P=-Wt V1+2W V2>0,
W=U1 V2-U2 V1, Wt=U1 V3-U3 V1.

For every lambda=1/rho with s_t<0, put w=-U1/V1 and z=lambda-w>0. Then

Q=s_t s_ttt-s_tt^2
 =-(U2+wV2)^2 +(P/V1)z +qa z^2<0.

Thus log(-s_t) is strictly concave wherever s_t<0 throughout this entire
rectangle, simultaneously for all rho>0. The same patch has U1>0 and Wt>0,
so any W zero in it is an upward crossing, the orientation needed for the
connectedness proof.

This is one continuous patch, not a cover of the entire parameter domain.
The remaining global proof requires an atlas/analytic argument covering the
other regions and the noncompact parameter limits.
