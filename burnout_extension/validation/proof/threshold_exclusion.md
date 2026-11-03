# Uniform exclusion of a descending critical segment near R=1

This proves a specific piece of the global problem. It does not prove that
there are at most three probability extrema on the entire positive R axis.

Write S=-log f, R=S/(1-exp(-S)), delta=1-theta, q=R(R-1),
t=log(1-1/R), and b=B/K. Define U,V by

s=-d(log b)/dt=U+lambda V, lambda=1/rho,
U=-q d_R log[C sqrt((R-1)R^(-theta))], V=q A_R.

**Proposition.** For every 0<S<=0.01, every theta in [0,1] and rho>0,
s>1 implies s_t>2.72. This bound is independent of K.
In particular, the set {s>1,s_t<0} is empty in this region. The critical
K curve is strictly increasing on every part of this region where it exists.

**Root-count consequence.** For every fixed finite rho,K>0 and theta in [0,1],
there is at most one stationary point of the unconditional BI probability on
this near-threshold interval. If it exists, it is a strict nondegenerate
maximum. Indeed the derivative sign factor is E=H(B)-s. At a stationary point,
s=H(B)>1 and E_t=-s D(B)-s_t<0. All its zeros must therefore be downward
crossings, so two distinct zeros are impossible. This also excludes degenerate
stationary points in the interval. The endpoint R is
0.01/(1-exp(-0.01)), approximately 1.00500833332.

## 1. A regular amplitude at S=0

Put phi(z)=(1-exp(-z))/z, extended by phi(0)=1, and w=1-v. Set

E(S,w)=sum_{k>=0} (-S)^k [sum_{j=0}^k w^j]/(k+2)!,
Psi(S,w)=(phi(S)-exp(-Sw))/S,

with their entire extensions at S=0. The Jensen gap of the exponential is

w exp(-S)+(1-w)-exp(-Sw)=S^2 w(1-w) E(S,w).

The exact amplitude identity in the main analysis, after s=Sv, gives

L(S,delta)=log(C/S^2)
 =log(phi(S))+log(E(S,1))+S+int_0^1 G(S,v,delta)dv,

G=(Psi(S,1-v) phi(S(1-v)) exp(-delta Sv)
      /[phi(S) E(S,1-v)] -1/phi(Sv))/v.

The numerator vanishes at v=0 identically; the quotient is removable.
At S=0 the integrand is -2 and L(0,delta)=-2-log(2).

## 2. Uniform analytic bound and a rigorous Taylor remainder

For |S|<=1/2, real delta in [0,1], and v,w in [0,1], the power series give

|phi(S)-1|<1/3,  |E(S,w)-1/2|<1/4,
|phi(Sw)|<4/3, |Psi(S,w)|<2,
|Psi_v|<5/3, |(phi(Sw))_v|<5/6, |E_v|<5/12.

For the E bound, sum the absolute coefficients at w=1:
sum_{k>=0} (k+1)(1/2)^k/(k+2)!=4-2 exp(1/2)<3/4.
All inequalities follow, for example, from 13/8<exp(1/2)<5/3.
The denominators are therefore nonzero and their logarithms have analytic
branches agreeing with their positive values at S=0.

Let A be the first term in the numerator defining G. Product differentiation
and the preceding bounds give

|A_v| < 10(20/9+5/3+4/3+40/9)=290/3,
|(1/phi(Sv))_v| <15/8.

Since the numerator vanishes at v=0, G=int_0^1 N_v(tv)dt, and hence
|G|<290/3+15/8. Also |log phi|<1, |log E(S,1)|<2.
It follows that |L(S,delta)|<120 throughout this complex S disk, uniformly
in delta. Consequently every Taylor coefficient satisfies |c_n(delta)|<=120*2^n.

The script `threshold_amplitude_series.py` computes c_0,...,c_6 exactly as
rational polynomials in delta. For j=1,2, the derivative remainder at
0<=S<=epsilon=0.01 is bounded by

120*2^7*(7!/(7-j)!)*epsilon^(7-j)
  /[1-2*epsilon*8/(8-j)].

This is a geometric majorant for the derivative tail. In particular, this is
an analytic remainder bound, not an estimate from comparison with samples.

## 3. Bounds for U and U_t

Let chi(S)=(1-phi(S))/S=int_0^1(1-z)exp(-Sz)dz. The regular part of log b is

g=L+(1/2)log chi-(delta/2)log phi,
log b=(5/2)log S+g+constant,
d/dt=S M d/dS, M=-chi/phi'.

Therefore
U=-M(5/2+S g'),
U_t=-S M[M'(5/2+S g')+M(g'+S g'')].

The integral formulas for phi,chi and their derivatives give endpoint bounds
uniformly on [0,epsilon]; e.g. exp(-epsilon)<=phi<=1,
-1/2<=phi'<=-exp(-epsilon)/2, and
exp(-epsilon)/2<=chi<=1/2. Combining them with the Taylor bounds above in
`threshold_exclusion.py` yields the outward-rounded enclosures

U in [-2.532,-2.468],  |U_t|<0.00479.

These somewhat loose enclosures suffice for the proposition.

## 4. Bounds for V_t/V

Let a=log(R)/S, extended by a(0)=1/2, and
k(S,z)=exp(-(delta+1)Sz)/phi(Sz).
An exact transformed action derivative is

A_R=exp(-delta S)-int_a^1 k(S,z)/z dz.

Because (log R)''=-Var_S(Z) for the exponentially tilted uniform variable
Z in [0,1], -1/4<=(log R)''<=0. Thus

1/2-S/8<=a<=1/2,  -1/8<=a'<=0.

It follows that
A_R>=exp(-epsilon)-exp(epsilon)log[1/(1/2-epsilon/8)]>0.2874.
Differentiating the transformed integral gives

|d_S A_R| <= 1+exp(epsilon)/(8 a_min)
              +(1-a_min)exp(epsilon)[2+exp(epsilon)/2] <3.

Moreover R_S>=1/2, which is equivalent to sinh S>=S. Hence |A_RR|<6.
Since V=q A_R and d/dt=q d/dR,

V_t/V=(2R-1)+q A_RR/A_R
 >=1-6 exp(epsilon)[exp(epsilon)-1]/(A_R lower bound)>0.78808.

## 5. The exclusion

Here V>0 and V_t>0. If s=U+lambda V>1, then
lambda>(1-U)/V. Therefore

s_t=U_t+lambda V_t
 >U_t+(1-U)V_t/V>2.72859.

All arithmetic bounds are reproduced with Arb balls by
`threshold_exclusion.py`, whose output is stored in
`validation/data/threshold_exclusion_certificate.json`.
