"""Certify three distinct fold points at theta=.55, K=exp(15/4).
These are local bifurcation certificates, not a global R root-count theorem.
"""
from flint import arb,acb,ctx
from interval_amplitude import amplitude_contour,amplitude_taylor_box,root_log_jet
from curvature_certificate import coefficients
from pathlib import Path
import csv,json,time
ctx.dps=50

def action_value(R,theta):
    z=root_log_jet(R,0);f=z[0].exp()
    def integrand(x,analytic):
        lx=x.log(analytic=analytic)
        if not lx.is_finite():return acb('nan')
        return (1-R*x)/(1-x)*(-theta*lx).exp()
    return acb.integral(integrand,acb(f),acb(1/R),abs_tol=arb('1e-17'),rel_tol=arb('1e-15'),eval_limit=50000).real

def evaluate(R,la,taylor=False):
    theta=arb('.55');K=(arb(15)/4).exp()
    C,_=amplitude_taylor_box(R,theta,order=4,tol='1e-14') if taylor else amplitude_contour(R,theta,order=4,tol='1e-14',coordinate='t')
    C=[x.real for x in C]
    U,V,U1,V1,U2,V2,U3,V3=coefficients(R,theta,C,coordinate='t')
    A=action_value(R,theta)
    logB=(arb(15)/4)+C[0]+(-la.log()+(R-1).log()-theta*R.log()-(2*arb.pi()).log())/2-A*la
    B=logB.exp();H=B.expm1()/B;D=B.exp()-H;Dp=B.exp()-D/B
    s=U+la*V;st=U1+la*V1;stt=U2+la*V2
    E=H-s;Et=-s*D-st;El=-(A+1/(2*la))*D-V
    Ett=-st*D+s*s*B*Dp-stt
    Etl=-V*D+s*(A+1/(2*la))*B*Dp-V1
    q=R*(R-1)
    J=[[Et/q,El],[Ett/q,Etl]]
    return [E,Et],J,{'Ett':Ett,'Elambda':El,'orientation':-El*Ett,'B':B}


def certify(row,radius='1e-6'):
    center=[arb(row['R0']),arb(row['lambda'])];rad=[arb(radius),arb(radius)]
    box=[c+arb(0,r) for c,r in zip(center,rad)]
    F0,J0,_=evaluate(*center)
    _,J,details=evaluate(*box,taylor=True)
    mid=[[v.mid() for v in row] for row in J0]
    det=mid[0][0]*mid[1][1]-mid[0][1]*mid[1][0]
    assert not det.contains(0)
    inverse=[[mid[1][1]/det,-mid[0][1]/det],[-mid[1][0]/det,mid[0][0]/det]]
    inverse=[[v.mid() for v in row] for row in inverse]
    M=[[arb(int(i==j))-sum(inverse[i][k]*J[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    offset=[-sum(inverse[i][k]*F0[k] for k in range(2))+sum(M[i][j]*arb(0,rad[j]) for j in range(2)) for i in range(2)]
    ratios=[(abs(v).upper()/r.lower()).upper() for v,r in zip(offset,rad)]
    contractions=[sum(abs(M[i][j]).upper()*rad[j]/rad[i] for j in range(2)).upper() for i in range(2)]
    success=all(v<1 for v in ratios+contractions)
    orientation=details['orientation']
    result={'status':'proved' if success and not orientation.contains(0) else 'unresolved',
            'theta':'.55','K':'exp(15/4)','R_box':str(box[0]),'lambda_box':str(box[1]),
            'R_enclosure':str(center[0]+offset[0]),
            'rho_enclosure':str(1/(center[1]+offset[1])),
            'Krawczyk_inclusion_ratios':[str(v) for v in ratios],
            'contraction_bounds':[str(v) for v in contractions],
            'orientation_product':str(orientation),
            'local_transition_with_increasing_rho':'annihilation' if orientation>0 else 'creation' if orientation<0 else 'unresolved',
            'B_enclosure':str(details['B'])}
    return result

if __name__=='__main__':
    data=Path(__file__).resolve().parents[1]/'data'
    rows=list(csv.DictReader((data/'reentrant_fold_centers.csv').open()))
    out=[]
    for row in rows:
        start=time.time();result=certify(row);result['seconds']=time.time()-start
        print(json.dumps(result),flush=True);out.append(result)
    (data/'reentrant_fold_certificates.json').write_text(json.dumps(out,indent=2))
