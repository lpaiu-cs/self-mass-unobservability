"""Apply the equivalent, pressure-regular canonical form without grid expansion."""
from pathlib import Path
import json
import signal
import time
import resource
import numpy as np
import def_gr_energy_modes as modes
import def_gr_canonical_regular as regular

OUT=modes.fem.OUT/'regular';assert not OUT.exists();OUT.mkdir()
paths=[Path(__file__),Path(regular.__file__),Path(modes.__file__),Path(modes.fem.__file__),modes.OUT/'result.json']
modes.write(OUT/'plan.json',dict(classification='Counterexample candidate',
    claim='Remove the free-surface pressure cancellation by a proved canonical momentum gauge, then apply the identical input to actual GR evolution and repeat the frozen acceptance gates.',
    method='Same native nodes, heat input,128/256/512 rational spaces and horizon. New canonical gauge has no reciprocal-pressure coefficients; no grid, basis, time or path-count expansion.',
    decision='If four projection histories pass relative2% and order1.5, run coefficient, exterior3R and every-other-node contrasts; no automatic expansion.',
    budget=dict(maximum_models=4,hard_seconds=120,CPU_threads=1,memory_GB=3,new_EOS_calls=0),
    forecast='Prior four equivalent-size models took22.60s; allow30-60s including the new symbolic cancellation, hard120s.',
    bindings={str(p):modes.digest(p) for p in paths}))
signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
_,symbolic=regular.symbolic();modes.write(OUT/'symbolic.json',symbolic)
modes.fem.canonical=regular;modes.OUT=OUT
model=modes.fem.Model(modes.evolution.BANK/'fine-bank.npz');P=modes.Projection(model)
_,Q,orth=P.basis(512);K,skew=modes.energy_matrix(model,Q)
cases={str(n):modes.series(model,K[:n,:n],Q[:,:n],'fine-'+str(n)) for n in [128,256,512]};comparisons={}
for field in modes.evolution.FIELDS:
    a,b,c=[np.array([v[field] for v in cases[str(n)]['history']]) for n in [128,256,512]]
    norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
    comparisons[field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
passed_projection=all(v['last']<.02 and v['order']>1.5 for v in comparisons.values())
if passed_projection:
    for label,bank,outer,coarse in [('coefficient','coarse-bank.npz',2,False),('outer','fine-bank.npz',3,False),('spatial','fine-bank.npz',2,True)]:
        other=modes.fem.Model(modes.evolution.BANK/bank,outer,coarse);proj=modes.Projection(other);_,pq,_=proj.basis(512);kr,_=modes.energy_matrix(other,pq)
        cases[label]=modes.series(other,kr,pq,label+'-512')
        for field in modes.evolution.FIELDS:
            c=np.array([v[field] for v in cases['512']['history']]);d=np.array([v[field] for v in cases[label]['history']])
            comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
passed=passed_projection and all(v.get('coefficient',1)<.02 and v.get('outer',1)<.002 and v.get('spatial',1)<.02 for v in comparisons.values())
result=dict(classification='Counterexample candidate',passed=passed,projection_passed=passed_projection,
    comparisons=comparisons,paths=list(cases),seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
    endpoint=cases['512']['history'][-1],orthogonality=orth,stiffness_skew=skew,linear_residual=P.error,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
modes.write(OUT/'result.json',result);signal.alarm(0);print('REGULAR',json.dumps(result),flush=True)
