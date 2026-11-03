"""Validated point checks of the algebra used to eliminate rho from folds."""
from flint import arb
from interval_amplitude import amplitude_contour
from curvature_certificate import coefficients

for r,theta in [('2','0'),('4','.5'),('9','.55')]:
    R,T=arb(r),arb(theta)
    C,_=amplitude_contour(R,T,order=4,tol='1e-13',coordinate='t')
    U,V,U1,V1,*_=coefficients(R,T,[x.real for x in C],coordinate='t')
    for b in ['.2','1','3']:
        B=arb(b);H=(B.exp()-1)/B;D=((B-1)*B.exp()+1)/B
        lam=(H-U)/V
        direct=H*D+U1+lam*V1
        eliminated=V*H*D+V1*H+V*U1-U*V1
        residual=eliminated-V*direct
        assert residual.contains(0)
    print(r,theta,'fold elimination identities enclosed',flush=True)
