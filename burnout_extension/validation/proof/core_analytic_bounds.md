# Explicit majorants for the compact-core interpolation

The calculation in `core_analytic_bounds.py` uses the holomorphic-strip
theorem, not point sampling. The S interval [.01,4] has a Bernstein ellipse
of parameter 2.5 with extrema
Re S in [-.88775,4.89775], |Im S|<=2.09475.
It lies at distance greater than .5 from the boundary of the rectangle
Re S in [-1.5,5.5], |Im S|<=2.8. Theta uses the Bernstein ellipse parameter
10 for [0,1], contained in Re theta in [-2.025,3.025], |Im theta|<=2.475.

Let sigma_-,sigma_+,tau be the outer S bounds and c=cos(tau/2)>0. Positive
Laplace transforms give explicit uniform lower bounds

|phi|>=c phi(sigma_+), |chi|>=c chi(sigma_+),
|-phi'|>=c[-phi'(sigma_+)], |E(S,w)|>=c[-phi'(sigma_+)].

For the last inequality E(sigma_+,w) decreases with w, so its minimum is
E(sigma_+,1)=-phi'(sigma_+). Upper bounds are the corresponding positive
Laplace integrals at sigma_-. The exact elementary expressions are used.

If B bounds |S| and e bounds |exp(-delta S v)|, product differentiation of
N in the amplitude identity gives

|N_v| <= e/(phi_low E_low) [Psi_v_up phi_up
   +Psi_up phi_v_up+|delta| B Psi_up phi_up
   +Psi_up phi_up E_v_up/E_low] +B phi_prime_up/phi_low^2.

Here Psi_up<=exp(-sigma_-)/2, Psi_v_up<=exp(-sigma_-),
phi_v_up<=B phi_prime_up, E_v_up<=B exp(-sigma_-)/6.
Because N(0)=0, this bounds N(v)/v. Adding the explicit logarithm and S
bounds yields a majorant for L and then for g.

For the action, the harmonic argument on the rectangle |Re S|<=6,
|Im S|<=pi gives Re a>=1/(exp 6+1). We also use |a|<10: if |S|>=1 it
follows from the bound on log phi; if |S|<=1, write
phi=exp(-S/2)psi and use |psi-1|<=|S|^2 exp(1/2)/24.
These bound log a. For H(z)=exp[-(delta+1)z]/phi(z),

|g_delta(z)|=|(H(z)-1)/z|
 <=sup |H'| <= exp_bound[(|delta+1|/phi_low)
                            +phi_prime_up/phi_low^2].

The segment from log R to S stays in the outer rectangle; the script checks
its real bounds explicitly. This bounds A_R and V. The formula for U bounds
it on the rectangle inset by .25. Cauchy circles of radius .25 there, and
radius .5 for m,V on the original rectangle, bound all needed derivatives.
All constants and the strict containment inequalities are checked with Arb.

## Tensor Chebyshev error

If |Z|<=M on the product of Bernstein ellipses with parameters rS,rT, its
bivariate Chebyshev coefficients obey |a_ij|<=4 M rS^(-i)rT^(-j).
At Lobatto nodes, a Chebyshev polynomial aliases to another one of norm at
most one. Summing the omitted/aliased coefficients therefore gives the
following uniform interpolation error for degrees nS,nT:

8 M [ rT/(rT-1) * rS^(-nS)/(rS-1)
        +rS/(rS-1) * rT^(-nT)/(rT-1) ].

This bound is independent of the sampled values. Sampling and interpolation
must still be performed with validated arithmetic, and the resulting
polynomial inequalities must still be proved. The majorant file alone is
not a compact-core sign certificate.
