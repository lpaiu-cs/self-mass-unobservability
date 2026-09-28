from pathlib import Path
import time,resource,json
import numpy as np
import solve_native_incident_reciprocal as solve
start=time.monotonic();cpu=time.process_time();error=None
solve.base.drive.native.deadline(30)
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
out=solve.OUT
try:
    solve.write(out/'material-range-plan.json',dict(classification='Counterexample candidate',claim='Identify why the saved new material early-state Richardson comparison fails before production. Reuse that state; vary arithmetic probes and locate dominant cells. No physical amplitude, clock or tolerance change.',seconds=30,total_cap=5600))
    solve.initialize(1);m=solve.Material(128,64);d=np.load(solve.MATERIAL/'pilot-64.npz');z=d['delta_scaled'];t=float(d['time']);rates={};eps={}
    for p in [8.,16.,32.,64.,128.]:
        m.min_probe=np.inf;rates[p]=m.rhs(t,z,p)[0];eps[p]=float(m.min_probe)
    rows=[]
    for p in [16.,32.,64.,128.]:
        diff=rates[p]-rates[p/2];norm=np.maximum(np.sum(abs(rates[p]),axis=1),1.)
        ids=np.argsort(abs(diff[1]))[-8:][::-1]
        rows.append(dict(probe=p,relative=(np.sum(abs(diff),axis=1)/norm).astype(float).tolist(),largest_momentum_cells=ids.tolist(),momentum_difference=np.asarray(diff[1,ids],float).tolist(),momentum_rate=np.asarray(rates[p][1,ids],float).tolist()))
    solve.write(out/'material-range-result.json',dict(classification='Counterexample candidate',precise=m.precise,rows=rows,eps=eps,physical_branch=m.physical_branch_ratio))
    np.savez_compressed(out/'material-range-rates.npz',**{str(k):v for k,v in rates.items()})
    print(json.dumps(rows),flush=True)
except Exception as exc:error=repr(exc);raise
finally:
    solve.write(out/'material-range-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=solve.sha(__file__)))

