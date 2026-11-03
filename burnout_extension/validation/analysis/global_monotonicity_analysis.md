# Analytical sign criterion for the unconditional BI probability

Working analysis, not manuscript text. R denotes R0; rho, theta, and K are fixed when differentiating in R. This classifies the BI formula, not the stochastic epidemic outside the approximation's validity.

Write b(R)=B_BI(R)/K, L=b'/b, and H(B)=(exp(B)-1)/B. Direct differentiation gives

P'(R) = exp(-B) B / R^2 * [R(R-1)L + H(B)].

All factors outside brackets are positive for R>1 and K>0. H is strictly increasing on (0,infinity), with H(0+)=1 and H(infinity)=infinity. Indeed H'(B)=((B-1)exp(B)+1)/B^2>0, since the numerator is zero at zero and its derivative is B exp(B)>0.

Let s(R)=-R(R-1)L(R).

* If s<=1, P'>0 for every K>0.
* If s>1, there is exactly one positive B_c solving H(B_c)=s. Define K_c(R;rho,theta)=B_c/b(R). Then P'<0 iff K<K_c, P'=0 iff K=K_c, and P'>0 iff K>K_c.

The positive solution can also be written using the lower real Lambert W branch:
B_c=-W_{-1}(-exp(-1/s)/s)-1/s.
The zero solution of exp(B)=1+sB is extraneous to H(B)=s when s>1.

This is a global, exact sign characterization of the BI expression. It does NOT by itself bound the number of intersections of K_c(R) with a fixed K. A complete root-count theorem needs bounds on the connected components of {s>1} and the turning points of K_c; it cannot be obtained from a finite scan.

At a smooth boundary,
(d/dR) log K_c = s' / (B_c H'(B_c)) - L,
s'=-(2R-1)L-R(R-1)L'.
A tangency K=K_c and K_c'=0 is the condition for a pair of simple stationary points to merge. This gives a mathematical interpretation of P'=P''=0.

Near R=1, the established amplitude asymptotic C_theta ~ 2 exp(-2)(R-1)^2 and action O(R-1) give P_uncond ~ 2 exp(-2) K sqrt(rho/(2pi)) (R-1)^(7/2) at fixed rho>0. With the corresponding differentiated expansion, P'>0 sufficiently near one. This may be an extremely narrow, very low-probability region for small rho, and is not a statement of uniform stochastic accuracy.

At infinity, the BI asymptotics imply P'<0 eventually: for theta<1, b~c R^(-(1+theta)/2); for theta=1, log b has derivative -1/rho+O(1/R). Thus every fixed finite-K BI curve rises near one and falls eventually. Generic simple-zero sign sequences start with + and end with -. This requires an odd number of sign-changing roots, but does not exclude five or more without additional analysis.

Approximation diagnostics used in the preview: flag if any of rho*(R-1)*R^(-theta), rho*R^(1-theta)/(R-1), rho*L_BL^2, or 1/(K*y_star) is >=1. These are warning diagnostics, not certified accuracy boundaries. The first two test necessary distances relative to the Laplace width; the actual finite-entry lower-end condition can be stricter. A flag-free point is not automatically accurate.

Validation: 142 finite-difference comparisons of the original probability formula, evaluated with stable quadrature, agree with the threshold sign prediction. An additional audit found that the legacy C_explicit endpoint polynomial is unreliable at extremely small x_f (examples R=50 and 300); those failures are preserved in analytic_monotonicity_legacy_evaluation_checks.csv. This does not affect the displayed validation curves at R<=6, or the heatmaps at R<=8; the large-R analysis uses stable transformed quadrature. No manuscript text has been modified during this analytical investigation.


## Further expansion: exact structure of the R-derivative sign as rho varies

At fixed R,theta,K let q=R(R-1), A=DeltaA_f,
c=(log C_theta)' + 1/[2(R-1)] - theta/(2R),
B=K C_theta sqrt(rho*(R-1)*R^(-theta)/(2pi)) exp(-A/rho).
The sign of the R derivative is the sign of E(rho)=c-A'/rho+H(B)/q.
This is a parameter-slice analysis of the R derivative, not a different research question about dP/drho.

Because B_rho=B*(1/(2rho)+A/rho^2),
rho^2 E_rho=A' + (A+rho/2) D(B)/q,
where D(B)=B H'(B)=((B-1)exp(B)+1)/B.
D>0 and D'>0 for B>0: the numerator of D' is ((B^2-B+1)exp(B)-1), whose derivative is B(B+1)exp(B)>0. Hence rho^2 E_rho is strictly increasing in rho.

If A'>0, E is strictly increasing from minus infinity to plus infinity. There is exactly one positive recruitment threshold rho_c, with P_R<0 below it and P_R>0 above it.
If A'<0, E decreases to one unique minimum then increases, and E tends to plus infinity at both rho endpoints. Thus there are zero or two positive rho roots (one double root at tangency), and P_R<0 on at most one interval in rho.
If A'=0, E is strictly increasing; its lower endpoint is c+1/q, which determines whether there is a threshold.

This gives an exact bound on the boundary multiplicity in the rho direction for every fixed R,theta,K. It is not yet a bound on the number of R roots at fixed rho,theta,K. The remaining question is the shape of the threshold curves as R varies.

For theta=0 the action and derivative simplify to elementary expressions:
A=1-R*x_f+(R-1)log[(1-1/R)/(1-x_f)],
A'=1/R+log[(1-1/R)/(1-x_f)].
For theta=1:
A'=1+log[(1-1/R)/(1-x_f)] > 1-log 2 >0.
To prove the bound, put f=x_f, R=-log(f)/(1-f); the inequality -log(f)>2(1-f)/(1+f) for 0<f<1 implies R(1+f)>2, hence (1-1/R)/(1-f)>1/2. Thus the theta=1 endpoint has exactly one rho threshold for the R derivative at every R>1 and K>0.

Numerical checks of the rho-derivative identity against direct rho differences of E have maximum scaled error 7.8e-7 on the checked cases. Action derivative zeros found numerically are recorded separately; their uniqueness for all theta<1 has not been proved here.


## Exact amplitude differentiation and independent simplification

With f=x_f, xstar=1/R, a=(1-f)f^(theta-1), F(u)=1-u+log(u)/R,
q(u)=(xstar-u)h_theta(u)/(u^2 F(u))-a/(u-f), J=int_f^1 q(u)du,

f'=f(f-1)/(1-Rf), a'=a f'[(theta-1)/f-1/(1-f)],
q(f)=a[(theta-2)/f-1/(1-f)-1/(xstar-f)+1/(2f(1-Rf))],
q_R=h_theta(u)/(R^2u^2)[-1/F+(xstar-u)log(u)/F^2]-a'/(u-f)-a f'/(u-f)^2,
J'=-q(f)f'+int_f^1 q_R du,
(log C)'=-f'/(1-f)+(-R^(-2)-f')/(xstar-f)-f'/f+J'/a-Ja'/a^2.
All derivatives on the right have explicit formulas; the cancellations must be preserved when evaluating q_R.

An independent integration-by-parts identity removes the lower-end subtraction:
w(u)=(1-u)u^(theta-1), a=w(f),
log C=log(1-f)+1/a * int_f^1 [-w'(u) log(F(u)/(u-f))+(w(u)-a)/(u-f)]du,
where w'(u)=u^(theta-2)[theta-1-theta u]. Its lower endpoint is removable, and its upper logarithmic singularity is integrable. This is an exact identity, not a new approximation.

## Small-rho local peak and why finite plots can appear to start decreasing

Put e=R-1. At fixed theta, C~2 exp(-2)e^2, A~(1-log 2)e. For rho->0, e=O(rho), fixed K and theta, B is small and
P_uncond ~ [2 exp(-2)/sqrt(2pi)] K sqrt(rho) e^(7/2) exp[-(1-log 2)e/rho].
Consequently e_peak~[7/(2(1-log 2))]rho=11.406...rho,
and P_peak~16.34192275 K rho^4 (coefficient evaluated separately).
This is a formal BI endpoint/small-rho result. It explains an invisible near-threshold local peak, not a claim of uniform stochastic accuracy. Fixed finite K is essential for this height expansion.

## Broader numerical exploration

The broader scan covers rho=1e-6 to 1e3 (125 log-spaced values), K=1 to 1e15 (11 log-spaced values), and 26 theta values including .99, .995, .999, .9999, and 1. In 35,750 parameter combinations only +,- and +,-,+,- sign sequences were observed. This is evidence, not an exhaustive continuous-parameter root-count theorem.
Near R=1 a local asymptotic expression is used below R-1=2e-4 to avoid coalescing-root loss of accuracy. At large R the exp[-(1-theta)R] contribution is retained through the nonuniform theta->1 crossover. It is dropped only once (1-theta)R is large; logarithmic R coordinates avoid overflow at enormous remote roots.
The uniform large-R evaluation agrees with direct quadrature on R=50,100,300,600 and theta=0,.5,.975,.99,.999,1: maximum log-amplitude error 5.1e-13, action absolute error 1.2e-13, and scaled amplitude derivative error 3.4e-9.

The missing analytical step remains an upper bound on the number of R roots or, equivalently, the number of turns/components of the critical boundary curves. None of the calculations here proves that five or more sign-changing roots are impossible at every continuous parameter value.


## Additional proof attempt: theta=0 action has exactly one maximum

For theta=0, A''=1/[R^2(R-1)]-f/(1-Rf). Its sign is opposite to D-1 where D=Rf(R^2-R+1). Put z=-log f>0, R=z/(1-exp(-z)). Multiplying D-1 by (exp(z)-1)^3 gives
N(z)=-exp(3z)+(z^3-z^2+z+3)exp(2z)+(z^2-2z-3)exp(z)+z+1.
The coefficient of z^n/n! is zero for n=0,...,4; for n>=2 it is
c_n=-3^n+2^(n-3)(n^3-5n^2+8n+24)+n^2-3n-3.
It is positive for n=5,...,12 and negative at n=13.
For n>=13,
c_(n+1)-3c_n=2^(n-3)(-n^3+11n^2-22n-16)-2n^2+8n+4<0,
so all later coefficients are negative. Thus the positive and negative Taylor coefficients have one sign change. Dividing N(z) by z^12 makes every nonconstant term strictly decreasing: positive coefficients have powers <=-1 and negative coefficients have powers >=1. Hence N(z)/z^12 decreases strictly from positive infinity to negative infinity, and N has exactly one positive zero.
Therefore A0'' is negative then positive. Since A0'(1+)=1-log 2>0 and A0'(infinity)=0 from below, A0' has exactly one zero; the action has exactly one global maximum, at the numerical root R~2.572629848.

This proves the action's shape, not the shape of P: the amplitude and establishment factors must still be included.

## Theta=1 amplitude and a failed sufficient-condition route

Integration by parts gives the exact identity
log C1 = [1/(1-f)] int_f^1 log F(u;R) du.
A tempting sufficient condition for a unique probability maximum is pointwise log-concavity in R of F(f+(1-f)t;R), followed by integration. That pointwise claim is false very close to t=0. For z=10 and t=1e-12, high-precision evaluation gives a positive second R-derivative of log F (about 0.00922585). Thus this particular route cannot be used as a proof. Numerical log-concavity of the integrated log C1 remains a conjectural intermediate statement, not an established result.

The endpoint simplification, critical-K criterion, rho-slice classification, and theta=0 action proof are established algebraic results. No root-count proof for the full BI probability over all rho,theta,K has been obtained. Do not promote the two observed curve types to an exhaustive theorem.


## Collected fully expanded first derivative of the amplitude

Combining the endpoint terms and differentiating the normalized integrand removes J and J' from the final formula. Let d=1-Rf, a=(1-f)f^(theta-1), v=f(f-1)/d, and g=(1-theta+theta*f)/d. Then

(log C)'=-(1-theta)(1-f)/d -1/(R*d)+(1-f)/(2*d^2)+int_f^1 V(u)du,

V(u)=(1-u)u^theta/[a*u^2*F(u)] * [(1/R-u)*(log(u)/(R^2*F(u))-g)-1/R^2] -v/(u-f)^2.

The combined integrand has a removable lower singularity; its separated terms must not be integrated independently. There are no uncomputed derivatives in this formula. Numerical quadrature of the combined expression agrees with automatic derivatives at 25 combinations (R=1.2,2,3,5,10; theta=0,.25,.5,.75,1), maximum absolute difference 3.45e-7, dominated by endpoint extrapolation at R=1.2.

The remaining mathematical problem is to bound the zeros of the resulting integral expression after including action and establishment. This collected derivative is not, by itself, a root-count theorem. Exploratory critical-K geometry checks found no curve with more than one peak or one trough, but this is not a proof on the continuous parameter domain.


## A sufficient theorem reducing the full root-count problem

Change independent variable to t=log(1-1/R), so d/dt=q d/dR, q=R(R-1). Set s(t)=-d(log b)/dt. Stationarity is H(B)=s. Only s>1 can contain stationary points. The critical curve is Kc=H^{-1}(s)/b.
Let D(B)=B H'(B)>0 and define F(s)=s D(H^{-1}(s)). Then
(d/dt) log Kc = [F(s)+s_t]/D.
Where s_t>=0 this is positive. On an interval with s_t<0, put h(t)=-s_t and Y(s)=h(t(s)). If log h is concave in t, Y is concave in s, since Y''(s)=(log h)''/h<=0.

F is strictly convex. Indeed with H_j=d^jH/dB^j, mu=H_2/H_1,
F_ss=2+2B mu+(H/H_1)mu+(HB/H_1)mu',
where mu'=(H_3 H_1-H_2^2)/H_1^2>=0 by Cauchy-Schwarz under the positive measure x exp(Bx)dx on [0,1]. All remaining terms are positive.
Thus Y-F is strictly concave, and its positive set is an interval. Kc can therefore turn down then up at most once on a single descending interval of s. If {s>1,s_t<0} has at most one relevant interval, Kc has at most one peak and one trough across its domain components, and a fixed K has at most three intersections. Endpoint derivative signs then restrict the generic probability shapes to one maximum or maximum-minimum-maximum.

The unproved model-specific hypotheses are now:
(1) {s>1,s_t<0} has at most one interval;
(2) log(-s_t) is concave on that interval.
These are sufficient, not known necessary conditions. They must be proved from C_theta and A, not assumed from numerical plots.

Because s=U+lambda V, lambda=1/rho, U=-q*c, V=q*A', the numerator of the log-concavity condition is a quadratic in lambda:
Q(lambda)=(U1+lambda V1)(U3+lambda V3)-(U2+lambda V2)^2,
where subscripts denote t derivatives. The domain s>1,s_t<0 is a lambda interval described by two linear inequalities. Thus for each R,theta, verifying the condition for every rho reduces to maximizing one explicit quadratic on that interval.
Fourth-order Taylor automatic differentiation was implemented to evaluate these quantities without finite-difference derivatives. At 592 sampled R,theta pairs with a nonempty lambda interval, the quadratic maximum is negative. This check covers all rho at each of those sampled R,theta pairs, but is NOT an interval proof between them.

Several stronger easy hypotheses fail: s itself need not be concave on its descending part; -V1 need not be log-concave everywhere; the ratio -U1/V1 need not be convex; its Wronskian derivative need not be positive everywhere. The corresponding exploratory data preserve counterexamples so those assumptions are not reused.

A potentially weaker route for interval connectedness is: on V1<0, show U1>0 and w=-U1/V1 has exactly one minimum. Then s_t<0 is a sublevel set of w and is connected. The sampled w curves have one minimum; proving this remains open. The Wronskian W=U1 V2-U2 V1 changes sign at its stationary points; W_t=U1 V3-U3 V1 is not positive everywhere, so monotonicity of W cannot be used without further restriction.

## Validated local curvature certificates (2026-11-02)

The new `validation/proof` routines use Arb real/complex balls, rather than
floating-point error estimates. They remain local certificates and do not cover
the entire continuous parameter domain.

Let S=-log(f), u=f exp(s). An exact amplitude representation is

log C = log(1-f)+log(1/R-f)-log(f) + int_0^S g(s) ds,

g(s)=(1/R-u)(1-u)exp((theta-1)s)/[(1-f)(f-u+s/R)]
      -1/[1-exp(-s)].

The singularities of the combined function g at 0 and S are removable. There
are no other zeros of F=f-f exp(s)+s/R in |Im s|<2 pi. To see this, a nonreal
zero s=x+iy with 0<y<pi would imply simultaneously
x> -log(Rf)>1-Rf (imaginary part) and
x=y cot(y)-Rf<1-Rf (real part). For pi<=y<2 pi the imaginary-part equation
has the wrong sign. Conjugation excludes the lower half-strip. Strict
concavity of F on the real axis leaves only 0 and S as its real zeros.
Consequently, a rectangle enclosing [0,S] with imaginary heights +/-4 gives

int_0^S g(s)ds = 1/(2 pi i) contour_int g(z) log[z/(z-S)] dz.

This identity permits validated differentiation without numerically subtracting
endpoint poles on the real integration path. Taylor jets implicitly solve for
f, and each differentiated contour integral is enclosed by Arb quadrature.
The implementation computes R derivatives or derivatives in
 t=log(1-1/R). Mixed theta derivatives follow by multiplying the first term of
g by s^m. A parameter Taylor enclosure uses center values/gradients and
full-rectangle second-derivative bounds; this reduces interval dependency.

For the all-rho curvature test, the lambda quadratic is centered inside its
feasible interval before maximization. Also, if s_ttt>0 throughout the interval,
then s_t*s_ttt-s_tt^2<0 follows immediately because feasibility requires s_t<0.

`local_certificate_checks.py` saves reproducible certificates, including:
- R in [1.9999,2.0001], theta in [-.0001,.0001]: Q upper < -715.
- R in [3.9999,4.0001], theta in [.4999,.5001]: Q upper < -18090.
- R in [8.9999,9.0001], theta in [.5499,.5501]: s_ttt lower >15752.
Each applies to every feasible lambda=1/rho, not just selected rho values.
The theta=0 box may be intersected with the model domain theta>=0.
Other tested boxes remain unresolved because interval enclosures are too wide.
Neither these local certificates nor the previous broad scans establish the
required global curvature or interval-connectedness hypotheses.

Arb implementation detail: raising a real ball containing zero to a general
power can return indeterminate even for exponent 2. The new code uses explicit
multiplication for such squares; this enlarges bounds safely and avoids NaNs.

## Exact fold-surface reduction: at most two fold points at fixed R,theta

This is an additional analytic result, independent of interval quadrature.
Retain t=log(1-1/R), lambda=1/rho, b=B/K, s=U+lambda V,
and H(B)=(exp(B)-1)/B. Put D(B)=B H'(B) and
F(s)=s D(H^{-1}(s)), s>1. The earlier argument proves F''(s)>0.
A degenerate stationary point of the probability (a fold of the stationary-point
surface) satisfies both

H(B)=s,  s_t=-s D(B)=-F(s).

Indeed, the exact derivative sign factor is H(B)-s; its t derivative at a zero
is -s D(B)-s_t, since B_t=-s B. Thus these two equations are necessary and
sufficient for P_R=P_RR=0 at an interior point.

If V != 0, eliminate lambda=(s-U)/V to obtain the single equation

G(s)=V F(s)+V1 s+V U1-U V1=0,

on s>1 with (s-U)/V>0. That admissible domain is an interval. Because
G''(s)=V F''(s), G is strictly convex for V>0 and strictly concave for V<0.
It has at most two roots. Each gives exactly one positive rho=V/(s-U),
B=H^{-1}(s), and K=B/b(R,rho,theta). The model restriction K>=1 is then
applied. Therefore, at each fixed R,theta with V!=0, there are at most two
fold points in the positive (rho,K) plane. This is a proved bound on the
parameter boundary, not yet a bound on the number of R extrema for fixed
rho,K.

If V=0, stationarity requires s=U>1 and B=H^{-1}(U). If V1!=0, the fold has
lambda=(-F(U)-U1)/V1 and hence at most one admissible point. If also V1=0,
there is no solution unless U1=-F(U); that exceptional case must be ruled out
or classified separately before claiming a universal two-sheet boundary.

For third-order degeneracy, put E=H(B)-s and differentiate the exact expression
E_t=-s D(B)-s_t. At a fold,
E_tt=s D(B)^2+s^2 B D'(B)-s_tt=F'(s)F(s)-s_tt.
Consequently the additional exact equation is
s_tt=F'(s) F(s)>0.
This locates possible cusps on the recovery side of the descending segment of
s. Proving their total count still requires global information about U,V.

## Uniform near-threshold theorem, not a numerical scan

`validation/proof/threshold_exclusion.md` proves, and
`threshold_exclusion.py` reproduces with outward-rounded arithmetic, that
for every 0<S=-log(f)<=0.01, all theta in [0,1], and every rho>0,

s>1 implies s_t>2.72859.

Thus {s>1,s_t<0} is empty on this whole region, independently of K. This
removes the coalescing-final-size-root endpoint from the two remaining
model-specific sufficient hypotheses. The proof uses an exact regularized
amplitude identity and a Cauchy remainder, not samples or fitted coefficients.

The key exact rescaling is
F(f exp(Sv))=(1-v)exp(-S)+v-exp[-S(1-v)]
           =S^2 v(1-v) E(S,1-v),
where E is an entire divided difference of the exponential. The regular
amplitude L=log(C/S^2) is analytic at S=0, uniformly in theta. An explicit
complex-disk bound |L|<120 on |S|<=1/2 supplies rigorous bounds for every
Taylor remainder. The first terms, with delta=1-theta, are

L=-2-log 2 -S/2+[-delta^2/12+delta/36+5/216]S^2
  +[delta^3/36-delta^2/36]S^3+O(S^4).

The coefficients through order 6 are exact rational polynomials in delta,
generated by `threshold_amplitude_series.py` and stored in
`validation/data/threshold_amplitude_series.json`. Independent validated
contour values agree with the expected order-7 discrepancy. The remainder
certificate itself comes from Cauchy's estimate, not that comparison.

In particular, U=-5/2-(1/6+delta/4)S+O(S^2), so U_t is negative immediately
above threshold. Therefore any global argument that assumes U_t>0 everywhere
is invalid. The stronger fact needed here is that s_t>0 whenever s>1; the
large positive lambda V_t term enforces it uniformly near S=0.

## Additional exact theta=1 amplitude monotonicity

For theta=1, the previously proved identity is
log C1=(1/(1-f)) int_f^1 log F(u) du. Put L=1-exp(-S) and u=1-L z,
0<z<1. Define g_z(S)=-log(1-z+z exp(-S)). It satisfies g_z(0)=0 and
is strictly concave in S because g_z'' is minus a positive Bernoulli variance.
Hence g_z(S)/S decreases strictly. At the same quantile z,

F(1-L z)=L [z-g_z(S)/S].

Both factors are strictly positive and strictly increasing in S. Averaging
the logarithm over z in (0,1) proves C1 strictly increases with S and R.
Therefore U=-1/2-q(log C1)_R<-1/2 for theta=1. This is an exact monotonicity
result for the amplitude, not yet a proof that its full persistence curve has
only one maximum. A sufficient further inequality would be (log C1)_RR<=0;
that inequality is numerically supported but not proved here.

A new exploratory check finds (log C)_RR+2(log C)_R/R<0 and the corresponding
log-amplitude-b curvature in x=1/R negative at 7,380 sampled R,theta pairs.
This inverse-R concavity remains unproved; no argument relies on it yet.
A stronger attempted action-integral log-convexity fails near threshold for
theta close to 1 and must not be assumed globally.

The uniform threshold theorem also gives an exact local root-count result:
for 1<R<=0.01/(1-exp(-0.01)) (~1.00500833332), every stationary point of the
unconditional BI probability is a strict nondegenerate maximum, and there is
at most one. At a stationary point E=H(B)-s=0, E_t=-sD(B)-s_t<0, so every
zero is a downward crossing. Two such zeros cannot occur without an intervening
upward crossing. This holds for every theta,rho,K in the objective, including
arbitrarily small rho and arbitrarily large fixed finite K.

A further useful exact form for the curvature check centers lambda at
w=-U1/V1 when V1<0. With z=lambda-w>0,

Q=-(U2+w V2)^2 +[-Wt+2W V2/V1]z +(V1 V3-V2^2)z^2,
W=U1 V2-U2 V1, Wt=U1 V3-U3 V1.

This keeps the constant term nonpositive exactly, even in parameter boxes
crossing W=0. If both the linear coefficient and the quadratic coefficient
are nonpositive (one strictly), Q<0 follows for z>0 without loss from
interval cancellation. The certificate code now uses this sufficient test.
Log-convexity of w itself is false in some sampled regions and is not assumed.

## Three certified local folds at fixed theta,K

A new inverse parameterization of the fold equation exposed a narrow region
near theta=.55 that was not resolved by the earlier logarithmic rho/K scans.
For theta=11/20 and K=exp(15/4), a validated two-variable Krawczyk calculation
proves three distinct local folds, using the full C and A integrals:

- R~6.791443, rho in .3111753521 +/-1.88e-11: pair annihilation;
- R~7.806470, rho in .33337104603 +/-9.19e-12: pair creation;
- R~9.555305, rho in .33716252287 +/-4.33e-12: pair annihilation.

Every weighted contraction bound is below .001. The orientation products
-E_lambda E_tt are respectively strictly positive, negative, positive.
See `proof/reentrant_folds.md`, `proof/reentrant_folds.py`, and the generated
`data/reentrant_fold_certificates.json`. This is a rigorous local result for
fixed parameters in the requested domain, not an upper-bound theorem.
A dense numerical check at theta=.55,K=exp(3.75) finds counts 3,1,3,1 in the
four intervening rho ranges. The local certificates do not certify those
total counts, because other roots/folds outside the boxes remain possible.
Any proposed parameter classification by just one critical rho is therefore
inadequate to describe the three local bifurcation boundaries that actually
exist. The full classification must accommodate these additional branches.
These rho values are not small; the result is about the formal BI function.

The diagnostic scan of 15 theta values from .53 to .60, and 1601 rho values
from .24 to .40, using 9000 R points per theta, found no more than two turns
of Kc(R) in the tested window. This does not exclude narrow structures between
points and is not being used as a proof of the global three-root bound.

## A numerical blind spot in fixed-R fold scans

For small B, the exact fold equation tends to
M(R,theta)=V U1+(1-U)V1=0, with lambda=(1-U)/V.
At fixed theta, a regular root of M generates a small-B fold branch. Its
K coordinate is K=B/b, which can be finite even when B is exponentially tiny.
For example, the numerical theta=.99 limiting branch has log b about -104,
so K near 1 requires B around exp(-104). Such a fold lies exponentially
close in R to the B=0 branch and cannot be captured on an ordinary R mesh.
Thus the minimum logK values in `fold_surface_scan.csv` are NOT minima of
the full fold curve. They are only minima among its sampled fixed-R points.
The inverse B parameterization, or a validated expansion around B=0, is needed
for those branches. This discovery changes the search strategy; the earlier
finite scans cannot be used to rule out them.

## Exact polynomial tail certificate for the reduced model

For a=delta R>=4, delta=1-theta, consider the leading tail geometry
U=pR, p=1-delta/2,
V=R^2 exp(-a)-R^(1-delta)/(1+delta), d/dt~R^2 d/dR.
Let z=(1+delta)R^(1+delta)exp(-a), and scale lambda so that
s/(pR)=1+eta(z-1). Then
s_t/(pR^2)=1+eta[z(2-a)-(1-delta)],
s_tt/(pR^3)=2+eta[z(a^2-6a+6)-(1-delta)(2-delta)],
s_ttt/(pR^4)=6+eta[z(-a^3+12a^2-36a+24)
                         -(1-delta)(2-delta)(3-delta)].
The negative-slope and positive-s constraints can be imposed exactly by
x=a-2, eta=(1+h)/(xz+1-delta),
z=(delta+h+w)/(x+1+h), with h,w>0.
After a=4+X and delta=alpha/(1+alpha), the numerator of the curvature
Q=s_t s_ttt-s_tt^2, multiplied by the positive denominator-clearing factor,
has 263 strictly negative nonzero coefficients in X,alpha,h,w. The w^2
coefficient is -4. The delta=1 endpoint has a separate negative-coefficient
certificate. `proof/reduced_tail_certificate.py` performs the exact rational
symbolic verification. Hence Q<0 in the reduced model throughout this whole
semi-infinite tail region, not just at points.

This is not yet a certificate for the full theory. The differences between
the exact amplitude/action derivatives and the reduced tail expressions must
be bounded in a way that preserves these signs. The global objective still
requires those bounds and the exact compact-region connectedness/curvature
certificates.

One useful additional tail identity retains the entire action power series:
V_tail=R(R-1)exp(-delta R)-R^(1-delta)/(1+delta)
       +sum_{n>=1} R^(1-delta-n)/[(delta+n)(delta+n+1)].
For R>1 and delta R>=4 its second R derivative is strictly positive.
The exponential term has second derivative
exp(-a)[a^2-4a+2+delta(2-a)] >=exp(-a)(a-1)(a-4)>=0;
the n=0 negative power term is convex for 0<=delta<=1, and every higher
positive term is convex since its exponent is <=0. The series may be
differentiated locally uniformly for R>1. For the pure logarithmic tail
amplitude U=p(R-1)-1/2, W=U1 V2-U2 V1=p q^3 V_RR>0, so the slope-threshold
ratio increases in that tail model. The exact exponentially small amplitude
correction remains to be incorporated before using this for the full W.

## Full-geometry Cauchy/Taylor patch (not a reduced tail model)

The new bivariate Taylor machinery composes the full geometry before taking
interval ranges, preserving cancellations that made earlier derivative-wise
bounds excessively wide. `proof/geometry_patch.py` gives a degree-14
certificate for S in [8.95,9.05], theta in [.545,.555], proving V1<0, qa<0,
and P=-Wt V1+2W V2>0. Consequently the negative-slope log-curvature inequality
holds on this entire rectangle for every rho>0. It also proves U1>0,Wt>0,
so every W zero there has the required upward orientation. The polydisk radii
are .5 and .1 and the remainder is bounded by Cauchy's estimate. See
`proof/geometry_patch.md` and `data/geometry_patch_S9_thetap55_degree14.json`.
The polynomial center values also overlap an independent scalar conversion.
This is genuine continuous coverage but is not yet a full-domain cover.

During the audit the threshold proof's written transformed action kernel
was corrected: its numerator is exp[-(delta+1)S z], not exp(-delta S z).
The associated conservative bound for |d_S A_R| changes from 2.015 to 2.522,
which remains below 3. The certified |A_RR|<6, the V_t/V lower bound, and the
threshold root-count result are unchanged and were rerun successfully.
`threshold_initial_signs.py` further proves uniformly that U_t/S<-0.1133
and V_t>0 throughout 0<S<=.01, all theta in [0,1].

The interval topology has a further important detail: near the theta=0
V1 zero, U1 is still negative. Thus relevant s_t<0 can occur even on V1>0.
The proposed shortcut that excludes all V1>=0 points is FALSE and must not
be used. The initial quadrant U1<0,V1>0 needs its own ratio analysis.
Exploratory checks in that quadrant find W zeros with Wt<0, and M zeros
with Mt<0, where M=V U1+(1-U)V1 and Mt=V U2+(1-U)V2. These orientations
would provide a more precise connectedness proof if established globally.
The sampled orientations alone are not a proof.

## Additional endpoint and proof-method audit

The threshold Cauchy/series argument now also gives the uniform signs
U1/S<-0.1133 and V1>0 for 0<S<=.01, all theta. It uses 20 theta intervals
only to sharpen the enclosure; the underlying bounds and coefficients are
analytic and valid continuously on every interval.

A proposed stronger initial-ratio orientation fails: there are W zeros with
Wt>0 in U1<0,V1>0 arbitrarily near threshold, despite the downward W zeros
observed farther out. Exact rational-polynomial series for F=W/S^3 through
degree 11, with the analytic complex bound |F|<3e8 on |S|<=.25, prove that
all its zeros on 0<S<=.01 have F_S>0, hence Wt>0. Existence follows from
opposite certified F signs at S=.005,delta=.80 and .82. See
`proof/threshold_W_bound.md`, `threshold_W_orientation.py`, and the generated
certificate. These zeros are outside the relevant set by threshold exclusion;
a global orientation hypothesis must therefore be conditioned on relevance.

One more precise candidate for the initial quadrant is to prove
V2<0 and Mt=V U2+(1-U)V2<0 whenever U1<0,V1>0,M<=0,
M=V U1+(1-U)V1. It implies W>0 there, because M<=0 gives
-U1/V1 >= (1-U)/V and V>0 in the initial increasing-V interval. At M=0,
W=-V1 Mt/V, so it also forces every M zero to be downward. Unlike an
unrestricted Wt<0 assertion, this is consistent with the rigorous threshold
counterexample. It is not yet established globally.

`proof/conditional_shape_classification.md` records the complete abstract
classification under the two remaining global connectedness/curvature
hypotheses: one or three ordinary extrema, a maximum plus a stationary
inflection on a fold boundary, or a quartic maximum at a cusp. The parameter
partition is K_L<K<K_H, with K_L=0 for a two-component s>1 domain. This
conditional theorem is not being promoted to a full-theory theorem.

## An integral-free forcing for the action-slope topology

Let g=V_R, q=R(R-1), qp=2R-1, f=x_f and delta=1-theta. Since
V=q A_R, one can eliminate the sole unevaluated integral A_R from

T=g'-2g/qp=2(3q+1) A_RR/qp+q A_RRR.

The exact symbolic reduction is
T=R^(-delta-1) P/[(R-1)qp]-f^delta Q/[qp(1-Rf)^3],
where P=2(1-delta)R^2+(3delta-2)R+1-delta is strictly positive for R>1,
delta in [0,1]. The explicit polynomial Q is stored by
`proof/action_slope_symbolic.py` in `data/action_slope_forcing_polynomials.json`.
Its signs and zero count have not been proved. Where Q>0, T<0 iff
q (Rf)^delta Q/[P(1-Rf)^3]>1; where Q<=0, T>0. This may permit a
logarithmic-derivative/positive-coefficient argument for the unique V1 zero,
without differentiating an integral numerically. It remains a research route.
