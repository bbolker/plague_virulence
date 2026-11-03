"""Taylor patch for the complete analytic geometry, preserving cancellation."""
from flint import arb,acb,ctx
from math import factorial,comb
from amplitude_patch import AmplitudePatch,pow_int,geometry_from_S
from bivariate_jet import BiJet
from complex_geometry import geometry_bounds
import time,json
from pathlib import Path


def build_jets(p,degree=12):
    BiJet.order=degree+4;BiJet.theta_order=degree
    S=BiJet({(0,0):p.S,(1,0):arb(1)})
    T=BiJet({(0,0):p.theta,(0,1):arb(1)})
    R=S/(1-(-S).exp());f=(-S).exp();de=1-T;Rp=R.Ds();q=R*(R-1);m=q/Rp
    C=BiJet({(i,j):v/(factorial(i)*factorial(j)) for j,a in enumerate(p.derivatives) for i,v in enumerate(a) if i+j<=degree+4})
    A2=((T-2)*R.log()).exp()/(R-1)-(-de*S).exp()*(de+T*f)/(1-R*f)
    Ad=Rp*A2
    S0=p.S;T0=p.theta;f0=(-S0).exp();R0=S0/(1-f0);delta=1-T0
    seeds={}
    for j in range(degree+1):
        def integrand(x,analytic):
            lx=x.log(analytic=analytic)
            if not lx.is_finite():return acb('nan')
            return (delta*lx).exp()*(-lx)**j/(1-x)
        integral=acb.integral(integrand,acb(f0),acb(1/R0),abs_tol=arb('1e-17'),rel_tol=arb('1e-15'),eval_limit=50000).real
        seeds[(0,j)]=(pow_int(S0,j)*(-delta*S0).exp()-integral)/factorial(j)
    seeds.update({(i+1,j):v/(i+1) for (i,j),v in Ad.c.items() if i+1+j<=degree+4})
    A1=BiJet(seeds)
    def Dt(G):return m*G.Ds()
    G=C+(R-1).log()/2-T*R.log()/2;us=[]
    for _ in range(4):G=Dt(G);us.append(-G)
    G=q*A1;vs=[G]
    for _ in range(3):G=Dt(G);vs.append(G)
    U,U1,U2,U3=us;V,V1,V2,V3=vs
    W=U1*V2-U2*V1;Wt=U1*V3-U3*V1
    return {'U':U,'V':V,'U1':U1,'V1':V1,'U2':U2,'V2':V2,
            'W':W,'Wt':Wt,'qa':V1*V3-V2*V2,'P':-Wt*V1+2*W*V2,
            'M':V*U1+(1-U)*V1}


def certify_patch(S='9',theta='.55',halfS='.05',halfTheta='.005',degree=14,rS='.5',rtheta='.1'):
    start=time.time();p=AmplitudePatch(S,theta,degree=degree,rS=rS,rtheta=rtheta)
    sb=acb(p.S+arb(0,p.rS),arb(0,p.rS));tb=acb(p.theta+arb(0,p.rt),arb(0,p.rt))
    bounds,calls=geometry_bounds(sb,tb)
    jets=build_jets(p,degree=degree)
    # Independent scalar conversion checks the formal geometry at the center.
    pointC=[p.derivatives[0][j] for j in range(5)]
    _,zz=geometry_from_S(p.S,p.theta,pointC)
    u,v,u1,v1,u2,v2,u3,v3=zz
    ww=u1*v2-u2*v1;wwt=u1*v3-u3*v1
    scalar={'U':u,'V':v,'U1':u1,'V1':v1,'U2':u2,'V2':v2,'W':ww,'Wt':wwt,
            'qa':v1*v3-v2*v2,'P':-wwt*v1+2*ww*v2,'M':v*u1+(1-u)*v1}
    assert all(scalar[name].overlaps(jets[name].at(0,0)) for name in jets)
    ds=arb(0,arb(halfS));dt=arb(0,arb(halfTheta))
    k=(abs(ds).upper()/p.rS.lower()).max(abs(dt).upper()/p.rt.lower())
    ratio=k*(degree+3)/(degree+2)
    assert ratio<1
    factor=(degree+2)*pow_int(k,degree+1)/(1-ratio)
    enclosures={}
    for name,jet in jets.items():
        error=(bounds[name]*factor).upper()
        enclosures[name]=jet.evaluate(ds,dt,degree)+arb(0,error)
    v=enclosures
    curvature=v['V1']<0 and v['qa']<0 and v['P']>0
    result={'S_center':S,'theta_center':theta,'S_halfwidth':halfS,'theta_halfwidth':halfTheta,
            'degree':degree,'Cauchy_S_radius':rS,'Cauchy_theta_radius':rtheta,'center_check':'passed',
            'curvature_status':'proved' if curvature else 'unresolved',
            'W_zero_orientation':'upward' if v['U1']>0 and v['V1']<0 and v['Wt']>0 else 'unresolved',
            'enclosures':{name:str(val) for name,val in v.items()},
            'complex_bounds':{name:str(val) for name,val in bounds.items()},
            'seconds':time.time()-start,'amplitude_calls':p.calls+calls}
    label=f'S{S}_theta{theta}_degree{degree}'.replace('.','p')
    path=Path(__file__).resolve().parents[1]/'data'/f'geometry_patch_{label}.json'
    path.write_text(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    ctx.dps=70
    print(json.dumps(certify_patch(),indent=2),flush=True)
