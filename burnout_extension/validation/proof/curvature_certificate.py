"""Local, rigorous all-rho curvature certificate; no global coverage claim."""
from flint import arb,acb,arb_series,ctx
from interval_amplitude import amplitude_contour,root_log_jet,amplitude_taylor_box


def coefficients(R,theta,C,coordinate='R'):
    ctx.cap=5
    Z=root_log_jet(R,4);RR=arb_series([R,1]);f=Z.exp();de=1-theta
    def integrand(x,analytic):
        logx=x.log(analytic=analytic)
        if not logx.is_finite():return acb('nan')
        return (de*logx).exp()/(1-x)
    integ=acb.integral(integrand,acb(f[0]),acb(1/R),abs_tol=arb('1e-14'),rel_tol=arb('1e-12'),eval_limit=50000)
    A1=(de*Z[0]).exp()-integ.real
    A2=RR.pow(theta-2)/(RR-1)-(-de*RR).exp()*(de+theta*f)/(1-RR*f) if False else (RR.log()*(theta-2)).exp()/(RR-1)-(de*Z).exp()*(de+theta*f)/(1-RR*f)
    A2,A3,A4=A2[0],A2[1],2*A2[2]
    q=R*(R-1);qp=2*R-1
    c0=C[1]+1/(2*(R-1))-theta/(2*R)
    c1=C[2]-1/(2*(R-1)**2)+theta/(2*R**2)
    c2=C[3]+1/(R-1)**3-theta/R**3
    c3=C[4]-3/(R-1)**4+3*theta/R**4
    U=-q*c0; Up=-qp*c0-q*c1; Upp=-2*c0-2*qp*c1-q*c2; Uppp=-6*c1-3*qp*c2-q*c3
    V=q*A1;Vp=qp*A1+q*A2;Vpp=2*A1+2*qp*A2+q*A3;Vppp=6*A2+3*qp*A3+q*A4
    U1=q*Up;V1=q*Vp;U2=q*(qp*Up+q*Upp);V2=q*(qp*Vp+q*Vpp)
    U3=q*(qp**2+2*q)*Up+3*q*q*qp*Upp+q**3*Uppp
    V3=q*(qp**2+2*q)*Vp+3*q*q*qp*Vpp+q**3*Vppp
    if coordinate=='t':
        T=arb_series([R,q,q*qp/2,q*(qp**2+2*q)/6,q*qp*(qp**2+8*q)/24])
        elem=(T-1).log()/2-theta*T.log()/2
        U,U1,U2,U3=[-C[n]-elem[n]*__import__('math').factorial(n) for n in range(1,5)]
    return U,V,U1,V1,U2,V2,U3,V3


def certify(R,theta,tol='1e-10',coordinate='t',taylor=False):
    C,calls=amplitude_taylor_box(R,theta,order=4,tol=tol) if taylor else amplitude_contour(R,theta,order=4,tol=tol,coordinate=coordinate)
    if any(not v.is_finite() for v in C):return {'status':'unresolved','reason':'integral enclosure nonfinite'}
    C=[v.real for v in C]
    U,V,U1,V1,U2,V2,U3,V3=coefficients(R,theta,C,coordinate=coordinate)
    lo=arb(0);hi=None
    # Necessary feasibility inequalities using coefficient upper bounds; these
    # deliberately overapproximate the union over all parameter values in box.
    for a,b in [(U-1,V),(-U1,-V1)]:
        au=a.upper();bu=b.upper()
        if bu>0:
            if au<0:lo=lo.max((-au/bu).lower())
        elif bu<0:
            if au<=0:return {'status':'empty','calls':calls}
            bound=(-au/bu).upper();hi=bound if hi is None else hi.min(bound)
        elif bu==0:
            if au<=0:return {'status':'empty','calls':calls}
        else:return {'status':'unresolved','reason':'constraint sign'}
    if hi is not None and hi<=lo:return {'status':'empty','calls':calls}
    if V1<0:
        # Center at the parameter-dependent zero of s_t, not at a fixed
        # numerical lambda. This preserves the exact nonpositive constant
        # -(U2+w*V2)^2, including rectangles crossing W=0.
        w=-U1/V1
        qa0=V1*V3-V2*V2
        beta0=-(U1*V3-U3*V1)+2*(U1*V2-U2*V1)*V2/V1
        if qa0<0 and beta0<=0:
            return {'status':'proved','reason':'negative centered coefficients',
                    'quadratic_upper':str(qa0.upper()),'linear_upper':str(beta0.upper()),'calls':calls}
    # On the feasible set S1=s_t<0. If s_ttt>0, the curvature
    # numerator S1*S3-S2^2 is strictly negative without any cancellation.
    v3lo=V3.lower();u3lo=U3.lower()
    if v3lo>=0:
        third_lower=(u3lo+lo*v3lo).lower()
    elif hi is not None:
        third_lower=(u3lo+hi*v3lo).lower()
    else:third_lower=None
    if third_lower is not None and third_lower>0:
        return {'status':'proved','reason':'positive third t derivative','third_lower':str(third_lower),'lambda_lo':str(lo),'lambda_hi':str(hi),'calls':calls}
    center=(lo+hi)/2 if hi is not None else lo
    S1,S2,S3=U1+center*V1,U2+center*V2,U3+center*V3
    qa=V1*V3-V2*V2;qb=S1*V3+V1*S3-2*S2*V2;qc=S1*S3-S2*S2
    a,c=qa.upper(),qc.upper()
    upper=None
    intervals=[(lo-center,arb(0),qb.lower()),(arb(0),None if hi is None else hi-center,qb.upper())]
    for lower,higher,b in intervals:
        candidates=[lower]
        if higher is not None:candidates.append(higher)
        elif not a<0:return {'status':'unresolved','reason':'unbounded quadratic'}
        if a<0:
            vertex=-b/(2*a)
            if vertex.upper()>=lower and (higher is None or vertex.lower()<=higher):candidates.append(vertex)
        for x in candidates:
            val=(a*x*x+b*x+c).upper()
            upper=val if upper is None else upper.max(val)
    return {'status':'proved' if upper<0 else 'unresolved','Q_upper':str(upper),'lambda_lo':str(lo),'lambda_hi':str(hi),'calls':calls}

if __name__=='__main__':
    import sys,time
    R=arb(sys.argv[1] if len(sys.argv)>1 else '2 +/- 0.000001')
    theta=arb(sys.argv[2] if len(sys.argv)>2 else '.5 +/- .000001')
    t=time.time();print(certify(R,theta));print('seconds',time.time()-t)
