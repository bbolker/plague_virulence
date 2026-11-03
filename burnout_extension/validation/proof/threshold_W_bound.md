# Threshold Wronskian series and its analytic remainder

This supplements the threshold exclusion theorem. It is useful for checking
which global ratio-shape hypotheses are false; these W zeros occur outside
the relevant negative-slope set.

Write delta=1-theta, phi=(1-exp(-S))/S and
chi=(1-phi)/S=int_0^1(1-z)exp(-Sz)dz. Let
m=dS/dt=-S chi/phi'. The preceding threshold proof gives |L|<120 on
|S|<=1/2 for L=log(C/S^2), uniformly in real delta in [0,1].
On that disk, the exponential series imply

|phi|>=2/3, |phi'|>=1/4, |chi|<=5/8, |chi|>=3/8,
|phi'|<=3/4, |chi'|<=5/18, |m|<=5/4.

For example, |phi'+1/2| is bounded by the positive-coefficient sum
4-2 exp(1/2)-1/2<1/4. Bounds for chi follow from its integral or series.

The normalized lower endpoint a=log R/S satisfies |a-1/2|<.03. Indeed
phi=exp(-S/2)psi, where psi is the moment generating function of a uniform
variable on [-1/2,1/2]. Therefore
|psi-1|<=|S|^2 exp(1/4)/24<|S|^2/18,
and |log psi|/|S|<4|S|/71<=2/71<.03.
The exact transformed action derivative is

A_R=exp(-delta S)-int_a^1
     exp[-(delta+1)S z]/[z phi(Sz)] dz.

The straight path from a to 1 stays in Re z>=.47, |z|<=1, with length<=.53.
The integrand modulus is at most 3/(.47*(2/3)), giving |A_R|<10. Since
|q|=|R(R-1)|<= (3/2)*(1/2)*(5/8)/(2/3)<1, we have |V=q A_R|<10.
These functions are holomorphic at S=0 and on the disk.

On the intermediate disk |S|<=3/8, Cauchy's estimate gives |L'|<=960.
For g=L+(1/2)log chi-(delta/2)log phi, the preceding bounds yield |g'|<962.
Thus U=-(m/S)(5/2+Sg') has |U|<1000 there.

On |S|<=1/4, Cauchy bounds on the surrounding disks give

|U'|<=8000, |U''|<=128000,
|V'|<=40, |V''|<=320, |m|<=5/8, |m'|<=5.

Consequently |U1|<=5000, |U2|<=75000, |V1|<=25, |V2|<=250, where
subscripts denote t derivatives. Hence

|W|=|U1 V2-U2 V1|<=3125000.

The first three Taylor coefficients of W vanish identically, as verified by
exact rational-polynomial arithmetic in `threshold_geometry_series.py`.
Therefore F=W/S^3 is holomorphic at zero. On |S|=1/4, |F|<=200000000;
the maximum modulus principle gives the same bound inside. The implementation
uses the looser bound 300000000 for its Cauchy remainder.

The exact coefficients of F through degree 11 are rational polynomials in
delta and ell=log 2. At 0<S<=.01 the value remainder is at most 5.243e-9
and its S derivative remainder at most 6.314e-6. Interval subdivision of delta
proves that every zero of F in this region has F_S>0. Because
W_t=m S^3 F_S at a W zero, these are upward W crossings.

There is such a zero: at S=.005 the certified F values at delta=.80 and .82
have opposite signs. Here U1<0 and V1>0 by the initial-sign certificate.
Thus the tempting universal assertion that all W zeros in this initial
quadrant have Wt<0 is false. The threshold exclusion theorem already proves
that the whole region has no s>1,s_t<0; the extra ratio turn does not produce
an additional probability extremum. Future global inequalities must therefore
be restricted to the relevant parameter set rather than imposed on all
initial-quadrant points.
