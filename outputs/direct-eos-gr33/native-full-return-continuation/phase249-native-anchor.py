"""Check the actual late high anchor before extending its coupled return."""
from pathlib import Path
import hashlib,json,os,resource,sys,time
import numpy as np
sys.path.insert(0,'verification')
import read_full_captured_history as complete
base=complete.prior;OUT=Path('native-full-return249-work/anchor');OUT.mkdir(parents=True)
resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3))
base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(600)
start=time.monotonic()
for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/part).mkdir(parents=True)
previous=Path('native-complete-return236-work')
for p in list((previous/'sweep-0').rglob('*.npz'))+[previous/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
    dst=OUT/p.relative_to(previous);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
complete.bind(base.endpoint.initialize,OUT=OUT)();m=base.base.run.owner.Model(64)
path=complete.saved(64);p=np.load(path);rows=[]
for k in [222,237]:
    t=p['joint_stage_times'][k];g=base.recovery.prior.restored_gas(m,p['joint_stage_conserved_scaled'][k])
    high=m.native(t,g);saved=p['joint_native_rates_scaled'][k]/m.units
    relative=(np.sum(abs(high-saved)*m.units,axis=0)/np.maximum(np.sum(abs(saved)*m.units,axis=0),base.LD('1e-290'))).astype(float).tolist()
    rows.append(dict(stage_index=k,time=float(t),native_relative=relative))
r=dict(classification='Counterexample candidate',passed=max(v for row in rows for v in row['native_relative'])<1e-12,
    rows=rows,seconds=time.monotonic()-start,peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    scope='Original returned-stage high native anchor on the first/last added coarse interval stages; not actual return acceptance.',
    bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),path]},final_charge_conclusion='unadjudicated')
(OUT/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));assert r['passed'],r
