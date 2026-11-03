"""Uniform initial signs U_t<0,V_t>0, supplementing threshold exclusion."""
import json
from math import factorial
from pathlib import Path
import sympy as sp
from flint import arb
import threshold_exclusion as t

lower=None;records=[]
for i in range(20):
    delta=(arb(i)/20).union(arb(i+1)/20)
    coeff=[None]
    for expression in t.expressions[1:7]:
        poly=sp.Poly(sp.sympify(expression),t.d);v=arb(0)
        for (n,),c in poly.terms():v+=arb(int(sp.numer(c)))/arb(int(sp.denom(c)))*t.power(delta,n)
        coeff.append(v)
    derivatives=[]
    for j in [1,2]:
        value=sum(coeff[n]*arb(factorial(n))/arb(factorial(n-j))*t.power(t.S,n-j) for n in range(j,7))
        derivatives.append(value+arb(0,t.remainder_bound(j)))
    cp,cpp=derivatives
    gp=cp+t.chip/(2*t.chi)-delta*t.phip/(2*t.phi)
    gpp=cpp+(t.chipp/t.chi-t.chip*t.chip/(t.chi*t.chi))/2-delta*(t.phipp/t.phi-t.phip*t.phip/(t.phi*t.phi))/2
    bracket=t.mp*(arb('2.5')+t.S*gp)+t.m*(gp+t.S*gpp)
    lower=bracket.lower() if lower is None else lower.min(bracket.lower())
    records.append({'delta':str(delta),'bracket_lower':str(bracket.lower())})
assert lower>0 and t.m>0 and t.A1min>0 and t.V1ratio>0
result={'status':'proved','S_domain':'(0,.01]','theta_domain':'[0,1]',
        'U_t_sign':'negative','V_t_sign':'positive',
        'positive_factor_in_minus_U_t_over_S':str((t.m.lower()*lower).lower()),
        'theta_subintervals':records}
print(json.dumps(result,indent=2))
path=Path(__file__).resolve().parents[1]/'data'/'threshold_initial_sign_certificate.json';path.write_text(json.dumps(result,indent=2))
