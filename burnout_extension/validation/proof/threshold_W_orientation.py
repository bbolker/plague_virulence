"""Uniform W-zero orientation near threshold, with an analytic Cauchy tail."""
from flint import arb,ctx
from math import factorial
from pathlib import Path
import json,sympy as sp
ctx.dps=60
root=Path(__file__).resolve().parents[1]/'data'
series=json.loads((root/'threshold_W_series.json').read_text())
d,ell=sp.symbols('delta ell');polys=[sp.Poly(sp.sympify(x),d,ell) for x in series['coefficients']]
N=series['order'];eps=arb('.01');r=arb('.25');MW=arb(300000000)
S=arb(0).union(eps);L=arb(2).log()
def power(v,n):
    out=arb(1)
    for _ in range(n):out*=v
    return out
def coefficient(p,delta):
    out=arb(0)
    for (i,j),v in p.terms():out+=arb(int(sp.numer(v)))/arb(int(sp.denom(v)))*power(delta,i)*power(L,j)
    return out
def rem(j):
    n=N+1;first=MW*power(1/r,n)*arb(factorial(n))/arb(factorial(n-j))*power(eps,n-j)
    ratio=eps/r*arb(n+1)/arb(n+1-j)
    return (first/(1-ratio)).upper()
uncertain=[];excluded=0;queue=[(arb(i)/100,arb(i+1)/100) for i in range(100)]
def evaluate(delta,Sbox,j):
    c=[coefficient(p,delta) for p in polys];value=arb(0)
    if j==0:
        for n in range(N,-1,-1):value=value*Sbox+c[n]
    else:
        for n in range(N,0,-1):value=value*Sbox+n*c[n]
    return value+arb(0,rem(j))
while queue:
    left,right=queue.pop();delta=left.union(right)
    val=evaluate(delta,S,0);deriv=evaluate(delta,S,1)
    if not val.contains(0):excluded+=1;continue
    if deriv>0:
        uncertain.append({'delta_lower':str(left),'delta_upper':str(right),'F_S_lower':str(deriv.lower()),'F_S_upper':str(deriv.upper())})
        continue
    assert right-left>arb('1e-6'),(delta,val,deriv)
    mid=(left+right)/2;queue.extend([(left,mid),(mid,right)])
left_value=evaluate(arb('.80'),arb('.005'),0)
right_value=evaluate(arb('.82'),arb('.005'),0)
assert left_value<0 and right_value>0
result={'status':'proved','S_domain':'(0,.01]','theta_domain':'[0,1]',
        'statement':'W=0 implies W_t>0', 'normalized_function':'F=W/S^3',
        'complex_radius':str(r),'complex_bound_MW':str(MW),'degree':N,
        'value_remainder':str(rem(0)),'derivative_remainder':str(rem(1)),
        'excluded_delta_boxes':excluded,'possible_zero_boxes':uncertain,
        'existence_at_S_005':{'delta_080_value':str(left_value),'delta_082_value':str(right_value)},
        'relevance':'These W zeros are outside s>1,s_t<0 by the threshold exclusion theorem'}
print(json.dumps(result,indent=2))
(root/'threshold_W_orientation_certificate.json').write_text(json.dumps(result,indent=2))
