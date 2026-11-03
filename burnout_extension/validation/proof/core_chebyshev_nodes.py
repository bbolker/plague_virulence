"""Validated Lobatto node generation with resumable provenance-aware caching.
This produces input enclosures only. It does not certify polynomial signs.
"""
import argparse,hashlib,json,sqlite3,time,math
from fractions import Fraction
from pathlib import Path
from flint import arb,acb,ctx
from complex_geometry import geometry_values
parser=argparse.ArgumentParser();parser.add_argument('--max-new',type=int,default=0);parser.add_argument('--ns',type=int,default=128);parser.add_argument('--nt',type=int,default=64)
args=parser.parse_args();ctx.dps=80;tol='1e-24'
base=Path(__file__).resolve().parent;data=base.parent/'data'
hashobj=hashlib.sha256()
for name in ['complex_geometry.py','interval_amplitude.py']:
    hashobj.update((base/name).read_bytes())
hashobj.update(b'S=[.01,4];theta=[0,1];dps=80;tol=1e-24;v1')
version=hashobj.hexdigest();db=sqlite3.connect(data/'core_validated_nodes.sqlite')
db.execute('PRAGMA journal_mode=WAL')
db.execute('CREATE TABLE IF NOT EXISTS nodes(version TEXT,nodekey TEXT,payload TEXT,seconds REAL,PRIMARY KEY(version,nodekey))')
db.execute('CREATE TABLE IF NOT EXISTS versions(version TEXT PRIMARY KEY,metadata TEXT)')
meta={'version':version,'S_interval':['.01','4'],'theta_interval':['0','1'],'dps':80,'tol':tol,'source_files':['complex_geometry.py','interval_amplitude.py']}
db.execute('INSERT OR IGNORE INTO versions VALUES (?,?)',(version,json.dumps(meta)));db.commit()
ns,nt=args.ns,args.nt
points=[(i,j) for i in range(ns+1) for j in range(nt+1)]
def level(i,n):return Fraction(i,n).denominator.bit_length()-1
points.sort(key=lambda ij:(max(level(ij[0],ns),level(ij[1],nt)),level(ij[0],ns)+level(ij[1],nt),ij))
start=time.time();new=0;done=0;cached=0
progress=data/'core_node_progress.json'
status={**meta,'ns':ns,'nt':nt,'done':0,'total':len(points),'new_this_run':0,'cached':0,
        'elapsed_this_run':0,'remaining_nodes':len(points),'state':'running'}
for i,j in points:
    key=f'{Fraction(i,ns)}|{Fraction(j,nt)}'
    if db.execute('SELECT 1 FROM nodes WHERE version=? AND nodekey=?',(version,key)).fetchone():
        done+=1;cached+=1;continue
    before=time.time()
    S=(arb(4)+arb('.01'))/2+(arb(4)-arb('.01'))/2*(arb.pi()*i/ns).cos()
    theta=arb('.5')+arb('.5')*(arb.pi()*j/nt).cos()
    vals,calls=geometry_values(acb(S),acb(theta),tol=tol)
    assert all(v.is_finite() and v.imag.contains(0) for v in vals.values())
    payload={'S':str(S),'theta':str(theta),'i_ratio':str(Fraction(i,ns)),'j_ratio':str(Fraction(j,nt)),
             'values':{name:str(v.real) for name,v in vals.items()},'calls':calls}
    # The text format is an enclosure; check that serialization preserves it.
    assert all(arb(payload['values'][name]).contains(value.real) for name,value in vals.items())
    elapsed=time.time()-before
    db.execute('INSERT INTO nodes VALUES (?,?,?,?)',(version,key,json.dumps(payload),elapsed));db.commit()
    new+=1;done+=1
    status={**meta,'ns':ns,'nt':nt,'done':done,'total':len(points),'new_this_run':new,'cached':cached,
            'elapsed_this_run':time.time()-start,'last_node':key,'last_node_seconds':elapsed,
            'remaining_nodes':len(points)-done,'state':'running'}
    progress.write_text(json.dumps(status,indent=2))
    if new%25==0 or new<=3:print(json.dumps({k:status[k] for k in ['done','total','new_this_run','elapsed_this_run','last_node_seconds']}),flush=True)
    if args.max_new and new>=args.max_new:break
else:
    status.update(state='completed',done=done,cached=cached,remaining_nodes=0,elapsed_this_run=time.time()-start)
    progress.write_text(json.dumps(status,indent=2))
    print('all requested nodes completed',flush=True)
if args.max_new and new>=args.max_new:
    status['state']='bounded_run_completed';progress.write_text(json.dumps(status,indent=2))
    print('bounded run completed',new,'nodes; seconds',time.time()-start,flush=True)
db.close()
