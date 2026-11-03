"""Uniform full-amplitude tail expansion and Cauchy derivative majorants.
The sector inequalities are proved in tail_amplitude_bounds.md. Bounds are
not fitted to the comparison points below.
"""
from math import factorial
from flint import arb,arb_series,ctx
from interval_amplitude import amplitude_contour
import json,time
from pathlib import Path

def coefficient(n,theta):
    if n==0:return 1/theta
    def moment(p,q):
        coeff=[arb(1)]
        for j in range(1,q):
            new=[arb(0)]*(len(coeff)+1)
            for k,v in enumerate(coeff):new[k]+=v*(j-theta);new[k+1]+=v
            coeff=new
        return arb(factorial(p))/arb(factorial(q-1))*sum(v*arb(p+1-k).zeta(theta) for k,v in enumerate(coeff))
    return moment(n,n)-sum(moment(n-1,l+1) for l in range(n-1))

def constants(N):
    assert N>=2
    kappa=arb(2).sqrt();c=1-arb(-5).exp()
    def T(n):return arb(2)**(1-n)*sum(arb(factorial(n))/arb(factorial(n-j))*(arb(2)/5)**j for j in range(n+1))
    A=kappa*(arb('2.5')+arb('6.25')*N)*arb('.5').exp()*factorial(N+1)*arb(2)**(N+2)
    B=kappa*(sum(T(n)/c**n for n in range(N+1))+sum(arb(n-1)/10*T(n-1)/c**(n-1) for n in range(2,N+1)))
    # Checks for the elementary uniform denominator estimates at Re S>=10.
    f=arb(-10).exp();ru=kappa*10/(1-f)*arb(-5).exp();zz=2*kappa/(1-f)*arb(-5).exp()
    assert ru<arb('.1') and zz<arb('.02') and arb('1.5')*f<arb('.0001')
    Aprime=(arb('1.12')/(1-f)/(1-zz)+(1+ru)*(1+arb(-5).exp())/(2*(1-f)*(1-zz)**2))
    assert Aprime<3
    return A.upper(),B.upper()

def derivative_bounds(S,theta,N=10):
    assert S>=14 and theta>=arb('.5') and theta<=1
    A,B=constants(N);delta=1-theta;lo=3*S/4;radius=S/4
    Rmax=(5*S/4)/(1-(-lo).exp())
    M=A*(-delta*lo).exp()/lo**(N+1)+(9*Rmax+B)*(-(1+delta)*lo/2).exp()+(5*S+40)*(-lo).exp()
    return [(arb(factorial(j))*M/radius**j).upper() for j in range(5)]

def reference(S,theta,N=10):
    ctx.cap=5;SJ=arb_series([S,1]);delta=1-theta
    K0=arb(0) if delta==0 else -arb.const_euler()-(delta*S).log()+(-delta*S).ei()
    deriv=((-delta*SJ).exp()-1)/SJ
    K=arb_series([K0]+[deriv[j-1]/j for j in range(1,5)])
    inv=1/SJ;power=arb_series([1]);H=arb_series([0])
    for n in range(N+1):
        if n:power*=inv
        H+=coefficient(n,theta)*power
    C=K-(-delta*SJ).exp()*H
    return [C[j]*factorial(j) for j in range(5)]

if __name__=='__main__':
    ctx.dps=100;N=10;out=[]
    for ss,tt in [('50','.5'),('100','.75'),('200','.99'),('200','1'),('1000','1')]:
        S,T=arb(ss),arb(tt);start=time.time()
        exact,_=amplitude_contour(S,T,order=4,tol='1e-45',coordinate='S')
        ref=reference(S,T,N);bounds=derivative_bounds(S,T,N)
        residuals=[(exact[j].real-ref[j]) for j in range(5)]
        valid=[abs(residuals[j]).upper()<bounds[j].lower() for j in range(5)]
        result={'S':ss,'theta':tt,'enclosed_by_bound':valid,'residuals':[str(v) for v in residuals],'bounds':[str(v) for v in bounds],'seconds':time.time()-start}
        out.append(result);print(json.dumps(result),flush=True)
    assert all(all(z['enclosed_by_bound']) for z in out)
    path=Path(__file__).resolve().parents[1]/'data'/'tail_amplitude_bound_checks.json';path.write_text(json.dumps({'order':N,'A_N':str(constants(N)[0]),'B_N':str(constants(N)[1]),'checks':out},indent=2))
