"""Reproduce local rectangle certificates, not an exhaustive domain cover."""
import json,time
from pathlib import Path
from flint import arb
from curvature_certificate import certify
cases=[('2 +/- .0001','.0 +/- .0001'),('2 +/- .001','.25 +/- .001'),('4 +/- .0001','.5 +/- .0001'),('9 +/- .0001','.55 +/- .0001'),('10 +/- .001','.8 +/- .001'),('30 +/- .001','.9 +/- .001')]
results=[]
for R,theta in cases:
    start=time.time()
    result=certify(arb(R),arb(theta),tol='1e-12',taylor=True)
    result.update(R=R,theta=theta,seconds=time.time()-start)
    results.append(result)
    print(json.dumps(result),flush=True)
path=Path(__file__).resolve().parents[1]/'data'/'local_curvature_certificates.json'
path.write_text(json.dumps(results,indent=2))
