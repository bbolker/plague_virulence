"""Wait for the explicitly launched scan, then rebuild both figures.
Does not launch or restart simulations and does not edit manuscript prose.
"""
import csv,os,subprocess,sys,time
from pathlib import Path

pid=int(sys.argv[1])
results=Path('validation/data/fig1_high_R_results.csv')
while True:
    ready=False
    try:
        rows=list(csv.DictReader(results.open()))
        ready=(len(rows)==180 and all(int(r['attempts'])==10000 and
               int(r['unresolved'])==0 for r in rows))
    except (OSError,ValueError,KeyError):
        pass
    if ready:
        subprocess.run(['/usr/local/bin/Rscript','validation/finalize_high_R_figures.R'],check=True)
        print('Automatic figure update completed.',flush=True)
        break
    try:
        os.kill(pid,0)
    except ProcessLookupError:
        raise SystemExit('Scan process ended before all results completed; figures left unchanged. Inspect fig1_high_R_scan.log.')
    time.sleep(10)
