"""Finalize corrected figures only after full scan and independent checks finish."""
import csv,os,subprocess,sys,time
from pathlib import Path

pid=int(sys.argv[1]);data=Path('validation/data')
while True:
    ready=False
    try:
        rows=list(csv.DictReader((data/'fig1_high_R_guarded_results.csv').open()))
        checks=list(csv.DictReader((data/'tau_bias_audit_20261104/cross_check_results.csv').open()))
        ready=(len(rows)==180 and all(int(r['attempts'])==10000 and
               int(r['unresolved'])==0 for r in rows) and len(checks)==20 and
               all(int(r['unresolved'])==0 for r in checks))
    except (OSError,ValueError,KeyError):
        pass
    if ready:
        subprocess.run(['/usr/local/bin/Rscript','validation/finalize_high_R_figures.R','--require-guarded'],check=True)
        print('Corrected validation figures updated; raw previous scan retained.',flush=True)
        break
    try:
        os.kill(pid,0)
    except ProcessLookupError:
        raise SystemExit('Corrected scan stopped before all completion checks passed; inspect logs. Figures not replaced.')
    time.sleep(10)
