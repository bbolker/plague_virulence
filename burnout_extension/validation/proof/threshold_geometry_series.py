"""Exact threshold geometry series using polynomial coefficient arithmetic.
No fitted coefficients or floating-point algebra are used.
"""
import sympy as sp,json,time
from pathlib import Path
from math import factorial
D=14
d,ell=sp.symbols('delta ell');P=lambda x:sp.Poly(x,d,ell,domain=sp.QQ)
z=P(0);one=P(1)
def const(v):return [P(v)]+[z]*D
def add(a,b):return [u+v for u,v in zip(a,b)]
def sub(a,b):return [u-v for u,v in zip(a,b)]
def scale(a,v):return [x*v for x in a]
def mul(a,b):return [sum((a[j]*b[n-j] for j in range(n+1)),z) for n in range(D+1)]
def inv(a):
    base=a[0].as_expr();assert not base.free_symbols
    out=[P(1/base)]
    for n in range(1,D+1):out.append(-sum((a[j]*out[n-j] for j in range(1,n+1)),z)/base)
    return out
def div(a,b):return mul(a,inv(b))
def deriv(a):return [(n+1)*a[n+1] for n in range(D)]+[z]
def log(a):
    v=mul(deriv(a),inv(a));base=a[0].as_expr()
    assert base in [1,sp.Rational(1,2)]
    return [z if base==1 else P(-ell)]+[v[n-1]/n for n in range(1,D+1)]
def shift(a,n):return [z]*n+a[:D+1-n]
def unshift(a,n):return a[n:]+[z]*n

t=time.time();root=Path(__file__).resolve().parents[1]/'data'
raw=json.loads((root/'threshold_amplitude_series_order14.json').read_text())['coefficients']
C=[P(sp.sympify(x).xreplace({sp.log(2):ell})) for x in raw]
phi=[P(sp.Rational((-1)**n,factorial(n+1))) for n in range(D+1)]
R=inv(phi);Rp=deriv(R);q=mul(R,sub(R,const(1)));m=div(q,Rp)
chi=unshift(sub(const(1),phi),1)
g=add(add(C,scale(log(chi),sp.Rational(1,2))),scale(log(phi),P(-d/2)))
U=sub(scale(unshift(m,1),-sp.Rational(5,2)),mul(m,deriv(g)))
a=unshift(log(R),1)
k=mul([P((-d-1)**n/factorial(n)) for n in range(D+1)],R)
A1=add([P((-d)**n/factorial(n)) for n in range(D+1)],log(a));ap=const(1)
for n in range(1,D+1):
    ap=mul(ap,a)
    A1=sub(A1,shift(scale(sub(const(1),ap),k[n]/n),n))
V=mul(q,A1);Dt=lambda f:mul(m,deriv(f))
U1=Dt(U);U2=Dt(U1);V1=Dt(V);V2=Dt(V1)
W=sub(mul(U1,V2),mul(U2,V1));assert all(x.is_zero for x in W[:3])
# W coefficient at degree D uses U,V only through degree D-1, all of which
# are known despite truncating the auxiliary shifted series at degree D.
coeff=[str(x.as_expr()) for x in W[3:]]
result={'series':'W/S^3','coordinate':'S=-log(x_f)','delta':'1-theta','ell':'log(2)',
        'order':len(coeff)-1,'coefficients':coeff,
        'U_linear':str(U[1].as_expr()),'V_linear':str(V[1].as_expr())}
(root/'threshold_W_series.json').write_text(json.dumps(result,indent=2))
print('W/S^3 first two coefficients',coeff[:2]);print('order',len(coeff)-1,'seconds',time.time()-t)
