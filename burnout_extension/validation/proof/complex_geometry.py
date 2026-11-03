"""Complex-polydisk bounds for analytic geometry functions (no division by V1)."""
from math import factorial
from flint import acb,arb,acb_series,ctx
from interval_amplitude import amplitude_contour,rectangular_exp as exp


def geometry_values(S,theta,tol='1e-10'):
    ctx.cap=5
    values,calls=amplitude_contour(S,theta,order=4,tol=tol,coordinate='S',use_rectangular_exp=True,contour_height='pi')
    if any(not v.is_finite() for v in values):raise ValueError('nonfinite amplitude on polydisk')
    SJ=acb_series([S,1]);R=SJ/(1-exp(-SJ));f=exp(-SJ);de=1-theta
    Rp=acb_series([n*R[n] for n in range(1,5)]);q=R*(R-1);m=q/Rp
    # Logarithmic x coordinate avoids an exponentially small moving endpoint.
    def fun(y,analytic):
        return S*exp(-(de+1)*S*y)/(1-exp(-S*y))
    lower=R[0].log()/S
    A1=exp(-de*S)-acb.integral(fun,lower,acb(1),abs_tol=arb(tol),rel_tol=arb(tol),eval_limit=50000)
    A2=exp((theta-2)*R.log())/(R-1)-exp(-de*SJ)*(de+theta*f)/(1-R*f)
    derivative=Rp*A2
    A1series=acb_series([A1]+[derivative[i-1]/i for i in range(1,5)])
    def Dt(G):return m*acb_series([n*G[n] for n in range(1,5)])
    G=acb_series([values[n]/factorial(n) for n in range(5)])+(R-1).log()/2-theta*R.log()/2
    U=[]
    for _ in range(4):G=Dt(G);U.append(-G[0])
    G=q*A1series;V=[G[0]]
    for _ in range(3):G=Dt(G);V.append(G[0])
    u,u1,u2,u3=U;v,v1,v2,v3=V
    W=u1*v2-u2*v1;Wt=u1*v3-u3*v1
    vals={'U':u,'V':v,'U1':u1,'V1':v1,'U2':u2,'V2':v2,'U3':u3,'V3':v3,
          'W':W,'Wt':Wt,'qa':v1*v3-v2*v2,'P':-Wt*v1+2*W*v2,
          'M':v*u1+(1-u)*v1,'Mt':v*u2+(1-u)*v2}
    if any(not x.is_finite() for x in vals.values()):raise ValueError('nonfinite geometry')
    return vals,calls

def geometry_bounds(S,theta):
    values,calls=geometry_values(S,theta)
    return {name:abs(value).upper() for name,value in values.items()},calls

if __name__=='__main__':
    import time
    ctx.dps=60
    S=acb(arb('9 +/- .25'),arb('0 +/- .25'));T=acb(arb('.55 +/- .1'),arb('0 +/- .1'))
    start=time.time();bounds,calls=geometry_bounds(S,T)
    print({k:str(v) for k,v in bounds.items()},'calls',calls,'seconds',time.time()-start)
