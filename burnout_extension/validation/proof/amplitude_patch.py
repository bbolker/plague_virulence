"""Bivariate Taylor patches with a validated complex-polydisk Cauchy remainder.
Input coordinate is S=-log(x_f). A patch is only valid if its complex contour
bound is finite and the real target box is strictly within the Cauchy radii.
"""
from math import factorial,comb
from flint import arb,acb,arb_series,ctx
from interval_amplitude import amplitude_contour
from curvature_certificate import coefficients


def pow_int(x,n):
    out=arb(1)
    for _ in range(n):out*=x
    return out

class AmplitudePatch:
    def __init__(self,S,theta,rS='.25',rtheta='.1',degree=12,tol='1e-17'):
        self.S=arb(S);self.theta=arb(theta);self.rS=arb(rS);self.rt=arb(rtheta)
        self.N=degree;self.calls=0
        sb=acb(self.S+arb(0,self.rS),arb(0,self.rS))
        tb=acb(self.theta+arb(0,self.rt),arb(0,self.rt))
        bound,calls=amplitude_contour(sb,tb,order=0,tol='1e-10',coordinate='S',use_rectangular_exp=True,contour_height='pi')
        self.calls+=calls
        if not bound[0].is_finite():raise ValueError('complex polydisk not certified')
        self.M=abs(bound[0]).upper()
        self.derivatives=[]
        for m in range(degree+1):
            values,calls=amplitude_contour(self.S,self.theta,order=4+degree-m,tol=tol,coordinate='S',theta_order=m)
            self.calls+=calls
            if any(not v.is_finite() for v in values):raise ValueError('center jet not certified')
            self.derivatives.append([v.real for v in values])

    def enclose(self,S,theta):
        S,theta=arb(S),arb(theta);ds=S-self.S;dt=theta-self.theta
        kappa=(abs(ds).upper()/self.rS.lower()).max(abs(dt).upper()/self.rt.lower())
        if not kappa<1:raise ValueError('box outside polydisk')
        out=[];errors=[]
        for j in range(5):
            value=arb(0)
            # Horner evaluation retains dependence on ds and dt within each
            # polynomial, while all numerical derivative errors stay enclosed.
            for m in range(self.N,-1,-1):
                part=arb(0)
                for i in range(self.N-m,-1,-1):
                    c=self.derivatives[m][j+i]/(factorial(i)*factorial(m))
                    part=part*ds+c
                value=value*dt+part
            ratio=kappa*(self.N+j+3)/(self.N+2)
            if not ratio<1:raise ValueError('Cauchy tail majorant needs subdivision')
            error=self.M*factorial(j)/pow_int(self.rS,j)*comb(self.N+j+2,j+1)*pow_int(kappa,self.N+1)/(1-ratio)
            error=error.upper();errors.append(error)
            out.append(value+arb(0,error))
        return out,errors


def geometry_from_S(S,theta,C):
    """Convert enclosed S derivatives to U,V and t derivatives, exactly."""
    S,theta=arb(S),arb(theta);ctx.cap=5
    SJ=arb_series([S,1]);RR=SJ/(1-(-SJ).exp());R=RR[0]
    Rp=arb_series([k*RR[k] for k in range(1,5)])
    m=RR*(RR-1)/Rp
    G=arb_series([C[n]/factorial(n) for n in range(5)])+(RR-1).log()/2-theta*RR.log()/2
    U=[]
    for _ in range(4):
        derivative=arb_series([k*G[k] for k in range(1,5)])
        G=m*derivative;U.append(-G[0])
    _,V,_,V1,_,V2,_,V3=coefficients(R,theta,[arb(0)]*5)
    return R,(U[0],V,U[1],V1,U[2],V2,U[3],V3)

if __name__=='__main__':
    import time
    ctx.dps=70
    start=time.time();p=AmplitudePatch('9','.55')
    C,err=p.enclose('9 +/- .01','.55 +/- .005')
    R,z=geometry_from_S(arb('9 +/- .01'),arb('.55 +/- .005'),C)
    U,V,U1,V1,U2,V2,U3,V3=z
    qa=V1*V3-V2*V2;W=U1*V2-U2*V1;Wt=U1*V3-U3*V1
    beta=-Wt+2*W*V2/V1
    print('M',p.M,'calls',p.calls,'seconds',time.time()-start)
    print('R',R)
    print('C derivatives',[str(v) for v in C])
    print('remainder bounds',[str(v) for v in err])
    print('qa',qa,'beta',beta)
