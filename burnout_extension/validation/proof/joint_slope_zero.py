"""Krawczyk isolation of the common U_t=V_t=0 degeneracy.
Independent of the running compact-core node cache and its source hash.
"""
from flint import arb,acb,acb_series,ctx
from math import factorial
from pathlib import Path
import csv,json,time
from interval_amplitude import amplitude_contour,rectangular_exp as exp
from complex_geometry import geometry_values

def evaluate(S,theta,tol='1e-20'):
    ctx.cap=5
    z,calls=geometry_values(acb(S),acb(theta),tol=tol)
    Ct,n=amplitude_contour(S,theta,order=2,theta_order=1,tol=tol,
                          coordinate='S',use_rectangular_exp=True,contour_height='pi')
    SJ=acb_series([S,1]);R=SJ/(1-exp(-SJ));f=exp(-SJ);d=1-theta
    Rp=acb_series([n*R[n] for n in range(1,5)]);q=R*(R-1);m=q/Rp
    def Dt(g):return m*acb_series([n*g[n] for n in range(1,5)])
    g=acb_series([Ct[n]/factorial(n) for n in range(3)])-R.log()/2
    U1theta=-Dt(Dt(g))[0]
    lower=R[0].log()/S
    def kernel(y,analytic):return S*S*y*exp(-(d+1)*S*y)/(1-exp(-S*y))
    A1theta=S*exp(-d*S)-acb.integral(kernel,lower,acb(1),abs_tol=arb(tol),rel_tol=arb(tol),eval_limit=50000)
    r=R[0];ff=f[0];qq=q[0]
    A2theta=exp((theta-2)*r.log())*r.log()/(r-1)-exp(-d*S)*(S*(d+theta*ff)+ff-1)/(1-r*ff)
    V1theta=qq*((2*r-1)*A1theta+qq*A2theta)
    vals=[z['U1'],z['V1']]
    J=[[z['U2']/m[0],U1theta],[z['V2']/m[0],V1theta]]
    assert all(a.is_finite() and a.imag.contains(0) for a in vals+[b for row in J for b in row])
    return [a.real for a in vals],[[a.real for a in row] for row in J],{k:str(z[k].real) for k in ['U','V','U2','V2','U3','V3']}

def certify(center,radius='1e-7'):
    rad=[arb(radius)]*2;box=[c+arb(0,r) for c,r in zip(center,rad)]
    F,J0,_=evaluate(*center);_,J,details=evaluate(*box)
    M0=[[a.mid() for a in row] for row in J0]
    det=M0[0][0]*M0[1][1]-M0[0][1]*M0[1][0]
    inv=[[M0[1][1]/det,-M0[0][1]/det],[-M0[1][0]/det,M0[0][0]/det]]
    inv=[[a.mid() for a in row] for row in inv]
    M=[[arb(i==j)-sum(inv[i][k]*J[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    off=[-sum(inv[i][k]*F[k] for k in range(2))+sum(M[i][j]*arb(0,rad[j]) for j in range(2)) for i in range(2)]
    inclusion=[(abs(off[i]).upper()/rad[i].lower()).upper() for i in range(2)]
    contraction=[sum(abs(M[i][j]).upper()*rad[j]/rad[i] for j in range(2)).upper() for i in range(2)]
    success=all(a<1 for a in inclusion+contraction)
    return {'status':'proved local unique common zero' if success else 'unresolved','S_box':str(box[0]),'theta_box':str(box[1]),
            'S_enclosure':str(center[0]+off[0]),'theta_enclosure':str(center[1]+off[1]),
            'inclusion_ratios':[str(a) for a in inclusion],'contraction_bounds':[str(a) for a in contraction],
            'geometry_on_box':details,'scope':'Local existence and uniqueness only; not exclusion of other common zeros.'}

if __name__=='__main__':
    ctx.dps=65;data=Path(__file__).resolve().parents[1]/'data'
    row=next(csv.DictReader((data/'joint_slope_zero_center.csv').open()))
    start=time.time();out=certify([arb(row['S']),arb(row['theta'])]);out['seconds']=time.time()-start
    print(json.dumps(out,indent=2),flush=True)
    (data/'joint_slope_zero_certificate.json').write_text(json.dumps(out,indent=2))
