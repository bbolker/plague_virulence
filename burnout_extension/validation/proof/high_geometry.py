"""Fourth t derivatives for local degeneracy analysis (full amplitude)."""
from flint import arb,acb,acb_series,ctx
from math import factorial
from interval_amplitude import amplitude_contour,rectangular_exp as exp

def geometry4(S,theta,tol='1e-22'):
    ctx.cap=6
    values,calls=amplitude_contour(S,theta,order=5,tol=tol,coordinate='S',use_rectangular_exp=True,contour_height='pi')
    SJ=acb_series([S,1]);R=SJ/(1-exp(-SJ));f=exp(-SJ);de=1-theta
    Rp=acb_series([n*R[n] for n in range(1,6)]);q=R*(R-1);m=q/Rp
    def fun(y,analytic):return S*exp(-(de+1)*S*y)/(1-exp(-S*y))
    lower=R[0].log()/S
    A1=exp(-de*S)-acb.integral(fun,lower,acb(1),abs_tol=arb(tol),rel_tol=arb(tol),eval_limit=50000)
    A2=exp((theta-2)*R.log())/(R-1)-exp(-de*SJ)*(de+theta*f)/(1-R*f)
    derivative=Rp*A2
    A1series=acb_series([A1]+[derivative[i-1]/i for i in range(1,5)])
    def Dt(G):return m*acb_series([n*G[n] for n in range(1,6)])
    G=acb_series([values[n]/factorial(n) for n in range(6)])+(R-1).log()/2-theta*R.log()/2
    U=[]
    for _ in range(5):G=Dt(G);U.append(-G[0])
    G=q*A1series;V=[G[0]]
    for _ in range(4):G=Dt(G);V.append(G[0])
    out={**{f'U{i}':U[i] for i in range(5)},**{f'V{i}':V[i] for i in range(5)}}
    out['Wtt']=U[2]*V[3]+U[1]*V[4]-U[4]*V[1]-U[3]*V[2]
    return out

if __name__=='__main__':
    import json,time
    from pathlib import Path
    ctx.dps=70;start=time.time()
    S=arb('1.41765549195593')+arb('0 +/- .00001')
    theta=arb('.139353662293934')+arb('0 +/- .00001')
    out=geometry4(acb(S),acb(theta))
    assert all(z.imag.contains(0) and z.is_finite() for z in out.values())
    result={'S_box':str(S),'theta_box':str(theta),'values':{k:str(v.real) for k,v in out.items()},'seconds':time.time()-start}
    result['signs']={k:('positive' if out[k].real>0 else 'negative' if out[k].real<0 else 'unresolved') for k in ['U2','U3','V2','Wtt']}
    print(json.dumps(result,indent=2))
    (Path(__file__).resolve().parents[1]/'data'/'joint_slope_local_convexity.json').write_text(json.dumps(result,indent=2))
