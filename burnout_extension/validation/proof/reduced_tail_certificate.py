"""Exact polynomial sign certificate for a reduced tail model.
The full amplitude's remainder is NOT yet bounded by this certificate.
"""
import sympy as s,json
from pathlib import Path
x,d,z,h,w=s.symbols('x delta z h w');alpha,X=s.symbols('alpha X')
k=1-d
v1=-x*z-k
v2=z*(x*x-2*x-2)-k*(k+1)
v3=z*(-x**3+6*x*x-8)-k*(k+1)*(k+2)
eta=(1+h)/(-v1)
Q=(1+eta*v1)*(6+eta*v3)-(2+eta*v2)**2
# a=x+2. The new h>0 is minus the normalized s_t, while w>0
# enforces positivity of the normalized s itself.
Z=(d+h+w)/(x+1+h)
num,den=s.fraction(s.factor(Q.subs(z,Z)))
assert s.factor(den)==(-d*h-d+h*x+h+w*x+x+1)**2
transformed=s.cancel(num.subs(d,alpha/(1+alpha))*(1+alpha)**4)
poly=s.Poly(transformed.subs(x,X+2).expand(),X,alpha,h,w)
assert all(c<0 for c in poly.coeffs())
# A strictly negative w^2 term remains for all finite alpha>=0.
assert poly.coeff_monomial(w**2)<0
endpoint=s.Poly(num.subs({d:1,x:X+2}).expand(),X,h,w)
assert all(c<0 for c in endpoint.coeffs())
result={'model':'reduced large-R U=pR, V=R^2 exp(-delta R)-R^(1-delta)/(1+delta)',
        'conditions':'a=delta R>=4, 0<=delta<=1, normalized s>0 and s_t<0',
        'transformation':'a=4+X; delta=alpha/(1+alpha); z=(delta+h+w)/(a-1+h)',
        'polynomial_terms':len(poly.coeffs()),'positive_coefficients':0,
        'strict_negative_w_squared_coefficient':str(poly.coeff_monomial(w**2)),
        'delta_one_terms':len(endpoint.coeffs()),
        'status':'proved for reduced model only; full-theory tail remainder remains'}
print(json.dumps(result,indent=2))
path=Path(__file__).resolve().parents[1]/'data'/'reduced_tail_polynomial_certificate.json';path.write_text(json.dumps(result,indent=2))
