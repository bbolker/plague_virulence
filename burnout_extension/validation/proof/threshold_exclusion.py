"""Rigorous near-threshold exclusion of {s>1,s_t<0}, uniform in theta,rho,K.
Uses the analytic Cauchy bound proved in threshold_exclusion.md.
"""
import json
from pathlib import Path
from math import factorial
import sympy as sp
from flint import arb,ctx
ctx.dps=50
N=6;M=arb(120);r=arb('.5');endpoint=arb('.01')
S=arb('.005 +/- .005');delta=arb('.5 +/- .5')
path=Path(__file__).resolve().parents[1]/'data'/'threshold_amplitude_series.json'
expressions=json.loads(path.read_text())['coefficients']
d=sp.Symbol('delta')
def power(x,n):
    v=arb(1)
    for _ in range(n):v*=x
    return v
def polynomial(text):
    p=sp.Poly(sp.sympify(text),d);out=arb(0)
    for (n,),coef in p.terms():
        rational=arb(int(sp.numer(coef)))/arb(int(sp.denom(coef)))
        out+=rational*power(delta,n)
    return out
coeff=[None]+[polynomial(x) for x in expressions[1:N+1]]
def remainder_bound(j):
    n=N+1
    first=M*power(1/r,n)*factorial(n)/factorial(n-j)*power(endpoint,n-j)
    ratio=endpoint/r*arb(n+1)/arb(n+1-j)
    return (first/(1-ratio)).upper()
def regular_derivative(j):
    val=arb(0)
    for n in range(max(1,j),N+1):
        val+=coeff[n]*arb(factorial(n))/arb(factorial(n-j))*power(S,n-j)
    return val+arb(0,remainder_bound(j))
cp,cpp=regular_derivative(1),regular_derivative(2)
e=(-endpoint).exp()
phi=e.union(arb(1));phip=arb('-.5').union(-e/2);phipp=(e/3).union(arb(1)/3)
chi=(e/2).union(arb('.5'));chip=arb(-1)/6;chip=chip.union(-e/6)
chipp=(e/12).union(arb(1)/12)
gp=cp+chip/(2*chi)-delta*phip/(2*phi)
gpp=cpp+(chipp/chi-chip*chip/(chi*chi))/2-delta*(phipp/phi-phip*phip/(phi*phi))/2
m=-chi/phip
mp=-chip/phip+chi*phipp/(phip*phip)
U=-m*(arb('2.5')+S*gp)
U1=-S*m*(mp*(arb('2.5')+S*gp)+m*(gp+S*gpp))
amin=arb('.5')-endpoint/8
A1min=e-endpoint.exp()*(1/amin).log()
A1Sbound=1+endpoint.exp()/(8*amin)+(1-amin)*endpoint.exp()*(2+endpoint.exp()/2)
assert A1Sbound<3
qmax=endpoint.exp()*(endpoint.exp()-1)
# |A1_S|<3 and R_S>=1/2 give |A_RR|<6.
V1ratio=1-6*qmax/A1min
lower=U1.lower()+(1-U.upper())*V1ratio.lower()
result={'S_upper':str(endpoint),'theta_domain':'[0,1]', 'rho_domain':'(0,infinity)',
        'regular_C_prime':str(cp),'regular_C_second':str(cpp),
        'U':str(U),'U_t':str(U1),'A_R_lower':str(A1min),'A_RS_abs_upper':str(A1Sbound),
        'V_t_over_V_lower':str(V1ratio),'s_t_lower_when_s_gt_1':str(lower),
        'status':'proved' if A1min>0 and V1ratio>0 and lower>0 else 'unresolved'}
print(json.dumps(result,indent=2))
assert result['status']=='proved'
(path.parent/'threshold_exclusion_certificate.json').write_text(json.dumps(result,indent=2))
