"""Expand the integral-free forcing for zeros of V_R; no sign theorem yet."""
import sympy as s,json
from pathlib import Path
R,f,d=s.symbols('R f delta',positive=True)
q=R*(R-1);qp=2*R-1;v=f*(f-1)/(1-R*f)
D=lambda z:s.diff(z,R)+v*s.diff(z,f)
b0=(d+(1-d)*f)/(1-R*f)
coef=2*(3*q+1)/qp
Q=s.factor((coef*b0+q*(d*v/f*b0+D(b0)))*qp*(1-R*f)**3)
P=2*(1-d)*R**2+(3*d-2)*R+1-d
print('P=',P)
print('Q=',s.expand(Q))
result={'forcing':'T=V_RR-2 V_R/(2R-1)',
        'identity':'T=R^(-delta-1)*P/[(R-1)(2R-1)]-f^delta*Q/[(2R-1)(1-R f)^3]',
        'P':str(P),'Q':str(s.expand(Q)),'status':'exact algebra only; global signs unproved'}
path=Path(__file__).resolve().parents[1]/'data'/'action_slope_forcing_polynomials.json';path.write_text(json.dumps(result,indent=2))
