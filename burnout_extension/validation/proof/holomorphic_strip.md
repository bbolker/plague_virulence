# A common holomorphic strip for the exact geometry

Write S=-log f, delta=1-theta, R=S/(1-exp(-S)), with the removable value
R(0)=1. The regularized log-amplitude L=log(C/S^2), the action derivative A_R,
and U,V and all of their t derivatives extend holomorphically to

                 |Im S|<pi,   theta in C.

The statement is joint holomorphy in S and theta. It does not assert that B,
which contains S^(5/2) near zero, is holomorphic at S=0. The geometry removes
that known fractional-power singularity.

## Positive Laplace transforms

Let phi(S)=int_0^1 exp(-Sz)dz, chi(S)=int_0^1(1-z)exp(-Sz)dz, and
-phi'(S)=int_0^1 z exp(-Sz)dz. The divided difference in the rescaled
amplitude has the exact representation

E(S,w)=int_{a,b>=0,a+b<=1} exp[-S(a+bw)] da db, 0<=w<=1.

Expanding the exponential and integrating a^p b^q gives the previously used
series sum (-S)^n(1+w+...+w^n)/(n+2)!, so the representation is exact.
All four are positive-measure Laplace transforms supported on [0,1].
For S=sigma+i tau, rotation by exp(i tau/2) gives

Re[exp(i tau/2) F(S)]
 >=cos(|tau|/2) int exp(-sigma z) dmu(z)>0, |tau|<pi.

Thus none of the transforms vanishes in the strip. For tau>0 their arguments
lie strictly between -tau and 0, with the analogous statement for tau<0.
Their principal logarithms are therefore holomorphic there. In particular,
R=1/phi is holomorphic and R_S=-phi'/phi^2 never vanishes.

## Exact amplitude

The regular identity is
L=log phi(S)+log E(S,1)+S+int_0^1 N(S,v,delta)/v dv,

N=Psi(S,1-v)phi(S(1-v))exp(-delta S v)/[phi(S)E(S,1-v)]
       -1/phi(Sv),
Psi=(phi(S)-exp[-S(1-v)])/S.

Psi is entire, and N(S,0,delta)=0 identically. The quotient by v is removable.
On every compact subset of the strip and every compact theta set, all
denominators have uniform positive lower bounds, including uniformly in v.
Integration over v therefore preserves holomorphy. This proves the amplitude
statement throughout the strip, not only on a small disk around S=0.

## The normalized logarithmic endpoint

Set mu(S)=1/S-1/(exp S-1), with mu(0)=1/2. Then
mu(-S)=1-mu(S). On Im S=+/-pi,

Re mu(x+i pi)=x/(x^2+pi^2)+1/(exp x+1),

which lies strictly between zero and one. For x>=0 both terms are positive
and their sum is at most 1/(2pi)+1/2<1; negative x follows by symmetry.
On vertical lines x=+M for sufficiently large M, the bounds

M/(M^2+pi^2)-1/(exp M-1)>0,
1/M+1/(exp M-1)<1

hold uniformly in Im S in [-pi,pi]. The negative vertical boundary follows
by symmetry. Applying the harmonic maximum/minimum principle on arbitrarily
large such rectangles proves

                        0<Re mu(S)<1.

Since log R=-log phi and (log R)'=mu,

a(S)=log R/S=int_0^1 mu(tS)dt

has 0<Re a<1 and is holomorphic, with a(0)=1/2. Thus log a is holomorphic.
Also Im log R lies between 0 and Im S, by the Laplace-transform argument
bound. The straight segment from log R to S stays inside |Im z|<pi.

## Exact action derivative

Define the entire-at-zero function

g_delta(z)=exp(-delta z)/(exp z-1)-1/z,
             g_delta(0)=-delta-1/2.

It is holomorphic on |Im z|<2pi and entire in delta. Logarithmic change of
variable in the exact action integral gives

A_R=exp(-delta S)+log a(S)
       -int_{log R(S)}^S g_delta(z) dz.

The preceding endpoint and path properties prove its holomorphic extension
on the whole strip. At zero this gives A_R(0)=1-log 2. The action A itself
extends by integrating A_R R_S from zero; the strip is simply connected.

## Root-count geometry

Since R-1=S chi/phi,
m=dS/dt=-S chi/phi' is holomorphic. If

g=L+(1/2)log chi-(delta/2)log phi,

then log b0=(5/2)log S+g, up to a constant, and

U=(chi/phi')[5/2+S g_S],   V=R(R-1) A_R.

Both are holomorphic at zero and throughout the strip; so are all t
derivatives, because d/dt=m d/dS. This supplies a common analytic domain for
certified polynomial approximation of every function in the global proof.
It is not, by itself, a sign or zero-count theorem for those functions.
