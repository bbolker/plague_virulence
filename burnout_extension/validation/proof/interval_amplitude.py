"""Validated contour quadrature for the amplitude and its R derivatives.
Uses Arb balls and Taylor arithmetic. Not yet a domain-wide certificate.
Run with PYTHONPATH=/Users/yuyang/.cache/burnout-extension/python-proof.
"""
from flint import arb, acb, arb_series, acb_series, ctx
from math import factorial
import math
ctx.dps=40
ctx.cap=5

def tight_real_product(a,b):
    values=[x*y for x in [a.lower(),a.upper()] for y in [b.lower(),b.upper()]]
    out=values[0]
    for value in values[1:]:out=out.union(value)
    return out

def tight_real_trig(x,cosine=False):
    lo,hi=x.lower(),x.upper();pi=arb.pi()
    if not hi-lo<2*pi:return arb(-1).union(arb(1))
    shift=arb(0) if cosine else pi/2
    # Floating approximations only seed the integer search. The while tests
    # establish, with Arb bounds, that both sentinels lie outside the interval.
    start=math.floor(float((lo-shift)/pi))-2
    end=math.ceil(float((hi-shift)/pi))+2
    while not (shift+start*pi).upper()<lo:start-=1
    while not (shift+end*pi).lower()>hi:end+=1
    fn=lambda z:z.cos() if cosine else z.sin()
    out=fn(lo).union(fn(hi))
    for k in range(start,end+1):
        critical=shift+k*pi
        if not (critical.upper()<lo or critical.lower()>hi):out=out.union(arb(1 if k%2==0 else -1))
    return out

def rectangular_exp(x):
    """Componentwise complex-rectangle enclosure, including Taylor jets."""
    if isinstance(x,acb_series):
        coeff=[rectangular_exp(x[0])]
        for n in range(1,ctx.cap):
            coeff.append(sum(k*x[k]*coeff[n-k] for k in range(1,n+1))/n)
        return acb_series(coeff)
    if not isinstance(x,acb):return x.exp()
    mag=x.real.lower().exp().union(x.real.upper().exp())
    return acb(tight_real_product(mag,tight_real_trig(x.imag,True)),tight_real_product(mag,tight_real_trig(x.imag)))

def root_log_series(R, order=4):
    ctx.cap=order+1
    def fp(r): return -(-r*(-r).exp()).lambertw()/r
    Rball=R[0]
    f=fp(Rball.lower()).union(fp(Rball.upper()))
    z0=f.log(); coeffs=[z0]; slope=1/Rball-f
    for n in range(1,order+1):
        z=arb_series(coeffs); F=1-z.exp()+z/R
        coeffs.append(-F[n]/slope)
    return arb_series(coeffs)

def root_log_jet(Rball, order=4):
    return root_log_series(arb_series([Rball,1]),order)

def amplitude_contour(Rball, theta, order=4, tol='1e-14', eval_limit=50000, coordinate='R', theta_order=0, use_rectangular_exp=False, contour_height=None):
    ctx.cap=order+1
    exp=rectangular_exp if use_rectangular_exp else lambda x:x.exp()
    if coordinate=='S':
        SJ=acb_series([Rball,1]) if isinstance(Rball,acb) else arb_series([Rball,1])
        RR=SJ/(1-exp(-SJ))
        z=acb_series(-SJ)
    elif coordinate=='t':
        J=arb_series([0,1]); RR=Rball/(1-(Rball-1)*(J.exp()-1))
    elif coordinate=='R':RR=arb_series([Rball,1])
    else:raise ValueError(coordinate)
    R=acb_series(RR)
    if coordinate!='S':z=acb_series(root_log_series(RR,order))
    f=exp(z); S=-z
    iR=(1-f)/S if coordinate=='S' else 1/R
    pref=(1-f).log()+(iR-f).log()-z
    Smid=float(S[0].real.mid()); Srad=float(S[0].real.rad())
    margin=arb(4); height=arb(4) if contour_height is None else arb.pi() if contour_height=='pi' else arb(contour_height)
    left=-margin; right=arb(str(Smid+Srad+4))
    pts=[acb(left,-height),acb(right,-height),acb(right,height),acb(left,height)]
    integrals=[];calls=[0]
    for n in range(order+1):
        def fun(s, analytic):
            calls[0]+=1
            try:
                if use_rectangular_exp and s.real.lower()>S[0].real.upper()+1:
                    v=exp(-z-s)
                    first=-exp(theta*s+z)*(v*iR-1)*(v-1)/((1-f)*(1-(f+s*iR)*v))
                else:
                    u=exp(z+s); F=f-u+s*iR
                    first=(iR-u)*(1-u)*exp((theta-1)*s)/((1-f)*F)
                second=1/(1-exp(-s)) if use_rectangular_exp else 1/(-(-s).expm1())
                g=first-second if theta_order==0 else first*s**theta_order
                if use_rectangular_exp:
                    a=-s if s.real.upper()<0 else s
                    b=S-s if s.real.upper()<0 else s-S
                    la=a.log(analytic=analytic);lb=b[0].log(analytic=analytic)
                    if not la.is_finite() or not lb.is_finite():return acb('nan')
                    kernel=la-b.log()
                else:
                    rat=acb_series(s)/(s-S)
                    check=rat[0].log(analytic=analytic)
                    if not check.is_finite(): return acb('nan')
                    kernel=rat.log()
                return (g*kernel)[n]
            except (ValueError, ZeroDivisionError):
                return acb("nan")
        val=acb(0)
        for i in range(4):
            val+=acb.integral(fun,pts[i],pts[(i+1)%4],abs_tol=arb(tol),rel_tol=arb(tol),eval_limit=eval_limit,depth_limit=30)
        integrals.append(val/(2*acb.pi()*acb(0,1)))
    coeffs=[(pref[n] if theta_order==0 else 0)+integrals[n] for n in range(order+1)]
    return [c*factorial(n) for n,c in enumerate(coeffs)],calls[0]

def amplitude_taylor_box(Rball,theta,order=4,tol='1e-13'):
    """Enclose t derivatives on a parameter rectangle using Taylor's theorem.

    The central value and gradient are computed at the midpoint. All second
    parameter derivatives are enclosed on the full rectangle. This preserves
    correlations lost by direct interval evaluation of the fourth derivative.
    """
    rc,tc=Rball.mid(),theta.mid(); dr,dt=Rball-rc,theta-tc
    centers=[]; boxes=[];calls=0
    for m in range(3):
        values,n=amplitude_contour(rc,tc,order=order+2-m,tol=tol,coordinate='t',theta_order=m)
        centers.append([v.real for v in values]);calls+=n
        values,n=amplitude_contour(Rball,theta,order=order+2-m,tol=tol,coordinate='t',theta_order=m)
        boxes.append([v.real for v in values]);calls+=n
    q=Rball*(Rball-1);p=2*Rball-1;qc=rc*(rc-1)
    out=[]
    for n in range(order+1):
        RR=(boxes[0][n+2]-p*boxes[0][n+1])/q**2
        RT=boxes[1][n+1]/q;TT=boxes[2][n]
        value=centers[0][n]+dr*centers[0][n+1]/qc+dt*centers[1][n]
        value+=dr*dr*RR/2+dr*dt*RT+dt*dt*TT/2
        out.append(acb(value))
    return out,calls

if __name__=='__main__':
    import sys,time
    r=arb(sys.argv[1] if len(sys.argv)>1 else '2')
    th=arb(sys.argv[2] if len(sys.argv)>2 else '.5')
    t=time.time();v,calls=amplitude_contour(r,th,order=4,tol='1e-10')
    for n,x in enumerate(v):print(n,x,flush=True)
    print('calls',calls,'seconds',time.time()-t,flush=True)
