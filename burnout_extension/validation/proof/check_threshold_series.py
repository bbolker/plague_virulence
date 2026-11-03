"""Compare exact Taylor coefficients to independent validated contour values.
This checks implementation agreement; it is not a Taylor remainder certificate.
"""
import json
from pathlib import Path
from flint import arb
from interval_amplitude import amplitude_contour
path=Path(__file__).resolve().parents[1]/'data'/'threshold_amplitude_series.json'
series=json.loads(path.read_text())['coefficients']
for ss in ['.1','.2']:
    for tt in ['0','.5','1']:
        S,T=arb(ss),arb(tt)
        # Preserve exact rational coefficients through the Arb conversion.
        import sympy as sp
        d=sp.Symbol('delta');coef=[]
        for expression in series:
            e=sp.sympify(expression).subs(d,sp.Rational(1)-sp.Rational(tt))
            if e.has(sp.log):v=-2-arb(2).log()
            else:v=arb(int(sp.numer(e)))/arb(int(sp.denom(e)))
            coef.append(v)
        approx=2*S.log()+sum(coef[n]*S**n for n in range(len(coef)))
        exact,_=amplitude_contour(S,T,order=0,tol='1e-14',coordinate='S')
        print(ss,tt,'residual / S^7',(exact[0].real-approx)/S**7,flush=True)
