"""Explicit Bernstein-ellipse majorants from the holomorphic-strip theorem.
No sampled function values enter these bounds.
"""
from flint import arb,ctx
from pathlib import Path
import json
ctx.dps=60

def maximum(values):
    v=values[0].upper()
    for x in values[1:]:v=v.max(x.upper())
    return v

def majorants(theta_rho='10'):
    # The S Bernstein ellipse for [.01,4] with rho=2.5 lies at distance >.5
    # from the boundary of this outer rectangle. A .25 intermediate inset
    # is used for bounds on U and its derivatives.
    sl,su,tau=arb('-1.5'),arb('5.5'),arb('2.8')
    B=(su*su+tau*tau).sqrt();cos=(tau/2).cos();assert cos>0
    assert tau<arb.pi()
    half=(arb(4)-arb('.01'))/2;center=(arb(4)+arb('.01'))/2;rs=arb('2.5')
    sa=half*(rs+1/rs)/2;sb=half*(rs-1/rs)/2
    assert center-sa-sl>arb('.5') and su-center-sa>arb('.5') and tau-sb>arb('.5')
    def phi(x):return -(-x).expm1()/x
    def chi(x):return (x-1+(-x).exp())/(x*x)
    def j1(x):return (1-(1+x)*(-x).exp())/(x*x)
    pl=cos*phi(su);cl=cos*chi(su);jl=cos*j1(su);El=jl
    pu=phi(sl);cu=chi(sl);ju=j1(sl)
    a=(arb(theta_rho)+1/arb(theta_rho))/4
    b=(arb(theta_rho)-1/arb(theta_rho))/4
    dl,du=arb('.5')-a,arb('.5')+a
    dab=(maximum([abs(dl),abs(du)])**2+b*b).sqrt()
    def expmax(left,right):return (maximum([-x*y for x in [left,right] for y in [sl,su]])+b*tau).exp()
    e=expmax(dl,du);ep=expmax(dl+1,du+1)
    em=(-sl).exp();psi=em/2;psiv=em;phiv=B*ju;Ev=B*em/6
    Av=e/(pl*El)*(psiv*pu+psi*phiv+dab*B*psi*pu+psi*pu*Ev/El)
    MG=Av+B*ju/(pl*pl)
    logbound=lambda low,up:maximum([abs(low.log()),abs(up.log())])+tau
    lp=logbound(pl,pu);lc=logbound(cl,cu);le=logbound(El,ju)
    ML=lp+le+B+MG
    Mg=ML+lc/2+dab*lp/2
    M0=cu/jl;Mm=B*M0
    MU=M0*(arb('2.5')+B*Mg/arb('.25'))
    # Re a(S)>1/(exp(6)+1) follows by the harmonic maximum principle.
    ca=1/(arb(6).exp()+1)
    assert 6/(36+arb.pi()*arb.pi())-1/(arb(6).exp()-1)>ca
    assert arb('.5')+1/(2*arb.pi())<1-ca
    # |a|<=|log phi| for |S|>=1; for |S|<=1 use the centered-uniform
    # moment generating function phi=exp(-S/2) psi.
    assert lp<10
    small=(arb('.5').exp()/24)
    assert arb('.5')+small/(1-small)<1
    loga=maximum([abs(ca.log()),arb(10).log()])+arb.pi()/2
    dplus=(maximum([abs(dl+1),abs(du+1)])**2+b*b).sqrt()
    hprime=ep*(dplus/pl+ju/(pl*pl))
    MA=e+loga+(B+lp)*hprime
    # The integration path from log R to S stays in the outer rectangle.
    assert -pu.log()>sl and -pl.log()<su
    MV=(1/pl)*(1/pl+1)*MA
    mp=Mm/arb('.5');mpp=2*Mm/(arb('.5')**2)
    def derivatives(bound,radius):
        v1=bound/radius;v2=2*bound/(radius*radius);v3=6*bound/(radius*radius*radius)
        return [bound,Mm*v1,Mm*(mp*v1+Mm*v2),Mm*((mp*mp+Mm*mpp)*v1+3*Mm*mp*v2+Mm*Mm*v3)]
    U,U1,U2,U3=derivatives(MU,arb('.25'))
    V,V1,V2,V3=derivatives(MV,arb('.5'))
    W=U1*V2+U2*V1;Wt=U1*V3+U3*V1
    bounds={'U':U,'V':V,'U1':U1,'V1':V1,'U2':U2,'V2':V2,'U3':U3,'V3':V3,
            'W':W,'Wt':Wt,'qa':V1*V3+V2*V2,'P':Wt*V1+2*W*V2,
            'M':V*U1+(1+U)*V1,'Mt':V*U2+(1+U)*V2}
    return {name:value.upper() for name,value in bounds.items()},{'L':ML.upper(),'g':Mg.upper(),'A_R':MA.upper(),'m':Mm.upper(),'phi_lower':pl.lower(),'E_lower':El.lower()}

def interpolation_error(M,ns,nt,rs=arb('2.5'),rt=arb(10)):
    return (8*M*(rt/(rt-1)*rs**(-ns)/(rs-1)+rs/(rs-1)*rt**(-nt)/(rt-1))).upper()

if __name__=='__main__':
    b,a=majorants();out={'status':'analytic bounds; not yet an interpolation certificate',
        'S_interval':['.01','4'],'theta_interval':['0','1'],'S_ellipse_rho':'2.5','theta_ellipse_rho':'10',
        'majorants':{k:str(v) for k,v in b.items()},'auxiliary':{k:str(v) for k,v in a.items()},
        'errors_at_128_64':{k:str(interpolation_error(v,128,64)) for k,v in b.items()}}
    print(json.dumps(out,indent=2))
    p=Path(__file__).resolve().parents[1]/'data'/'core_analytic_majorants.json';p.write_text(json.dumps(out,indent=2))
