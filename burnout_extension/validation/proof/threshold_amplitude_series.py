"""Exact formal series at S=-log(x_f)=0; coefficients are rational polynomials.
The constant of log(C/S^2) is -2-log(2). No numerical fit is used.
"""
import sys,json
from pathlib import Path
import sympy as sp
from math import factorial
N=int(sys.argv[1]) if len(sys.argv)>1 else 6
w,d=sp.symbols('w delta'); P=lambda x:sp.Poly(x,w,d,domain=sp.QQ).as_expr()
def mul(a,b):
    return [P(sum(a[j]*b[k-j] for j in range(k+1))) for k in range(N+1)]
def inv(a):
    out=[1/a[0]]
    for k in range(1,N+1):out.append(P(-sum(a[j]*out[k-j] for j in range(1,k+1))/a[0]))
    return out
def logcoef(a):
    b=inv(a)
    return [sp.log(a[0])]+[P(sum(j*a[j]*b[k-j] for j in range(1,k+1))/k) for k in range(1,N+1)]
phi=[sp.Rational((-1)**k,factorial(k+1)) for k in range(N+1)]
phiw=[phi[k]*w**k for k in range(N+1)]
phiv=[phi[k]*(1-w)**k for k in range(N+1)]
E=[sp.Rational((-1)**k,factorial(k+2))*sum(w**j for j in range(k+1)) for k in range(N+1)]
Psi=[(-1)**(k+1)*(sp.Rational(1,factorial(k+2))-w**(k+1)/factorial(k+1)) for k in range(N+1)]
expd=[(-d*(1-w))**k/factorial(k) for k in range(N+1)]
a=mul(mul(mul(mul(Psi,phiw),expd),inv(phi)),inv(E)); b=inv(phiv)
ints=[]
for k in range(N+1):
    q,r=sp.div(P(a[k]-b[k]),1-w,w)
    assert sp.expand(r)==0
    integral=P(sp.integrate(q,(w,0,1)))
    ints.append(integral)
    print('integral',k,integral,flush=True)
logs1=logcoef(phi);logs2=logcoef([x.subs(w,1) if hasattr(x,'subs') else x for x in E])
coeff=[sp.expand(logs1[k]+logs2[k]+ints[k]+(1 if k==1 else 0)) for k in range(N+1)]
result={'coordinate':'S=-log(x_f)','delta':'1-theta','series':'log(C/S^2)','order':N,'coefficients':[str(x) for x in coeff]}
print(json.dumps(result,indent=2))
filename='threshold_amplitude_series.json' if N==6 else f'threshold_amplitude_series_order{N}.json'
path=Path(__file__).resolve().parents[1]/'data'/filename;path.write_text(json.dumps(result,indent=2))
