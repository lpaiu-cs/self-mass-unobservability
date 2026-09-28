"""Same spatial candidates, with square-root energy arithmetic and a fixed budget."""
from pathlib import Path
import json
import signal
import resource
import time
import numpy as np
import def_gr_spatial_repair as task
import def_gr_energy_factor as stable

PARENT=task.OUT;OUT=PARENT/'factor';assert not OUT.exists();OUT.mkdir();task.OUT=OUT
task.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='95d99987',
    claim='Test whether square-root energy arithmetic repairs the unstable projected response, then finish the same fixed degree1/2/4 spatial comparison.',
    unchanged='Same polynomial cells, composite quadrature, heat input,128/256/512 spaces, readouts and acceptance gates. No modal clipping, degree or grid expansion.',
    formula='K_r+shift*M_r=X.T X, evaluated by QR followed by SVD instead of forming the normal matrix. Shift is derived from local potential positivity before response evaluation.',
    decision='First run degree4 with75s cap. Only if all four projection comparisons pass relative2% and order1.5, continue degrees2,1. Only if degree2/4 spatial differences all decrease and are below2%, execute the three original coefficient/outer/quadrature contrasts. Total cap300s; stop on failure without expansion.',
    forecast='Previous degree4 case48.65s included a cold symbolic setup. QR instead of normal-matrix formation is unmeasured; first-case bound75s. Reuse saved raw case as the arithmetic contrast. Entire remaining workflow capped300s.',
    budget=dict(first_case_seconds=75,total_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
    bindings={str(p):task.digest(p) for p in [Path(__file__),Path(stable.__file__),Path(task.__file__),PARENT/'plan.json',PARENT/'measured-budget.json']}))
task.write(OUT/'control.json',stable.control())
signal.alarm(75);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic();cases={};propagation={};checks=[]
prop_pass=True
for degree in [4,2,1]:
    began=time.monotonic();model=task.Model(degree,task.BANK/'fine-bank.npz');P=task.Projection(model)
    setup=time.monotonic()-began;began_basis=time.monotonic();_,Q,orth=P.basis(512);basis_seconds=time.monotonic()-began_basis
    began_factor=time.monotonic();R,shift,meta=stable.factor(model,Q);factor_seconds=time.monotonic()-began_factor
    rows={n:stable.series(model,R[:n,:n],shift,Q[:,:n],f'p{degree}-{n}') for n in [128,256,512]}
    cases[str(degree)]=rows[512];propagation[str(degree)]={}
    for field in task.FIELDS:
        a,b,c=[np.array([h[field] for h in rows[n]['history']]) for n in [128,256,512]]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        propagation[str(degree)][field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
    checks.append(dict(degree=degree,seconds=time.monotonic()-began,setup_seconds=setup,basis_seconds=basis_seconds,factor_seconds=factor_seconds,
        orthogonality=orth,linear_residual=P.error,**meta));del model,P,Q,R
    prop_pass=all(v['last']<.02 and v['order']>1.5 for p in propagation.values() for v in p.values())
    if not prop_pass:break
    signal.alarm(max(1,int(300-(time.monotonic()-start))))
comparisons={};spatial_pass=False
if len(cases)==3:
    for field in task.FIELDS:
        a,b,c=[np.array([h[field] for h in cases[str(n)]['history']]) for n in [1,2,4]]
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
        comparisons[field]=dict(previous=float(d1),last=float(d2),decreased=bool(d2<d1))
    spatial_pass=all(v['last']<.02 and v['decreased'] for v in comparisons.values())
if prop_pass and spatial_pass:
    for label,bank,outer,nquad in [('coefficient','coarse-bank.npz',2,6),('outer','fine-bank.npz',3,6),('quadrature','fine-bank.npz',2,8)]:
        model=task.Model(4,task.BANK/bank,outer,nquad);P=task.Projection(model);_,Q,orth=P.basis(512);R,shift,meta=stable.factor(model,Q)
        cases[label]=stable.series(model,R,shift,Q,label+'-p4-512')
        for field in task.FIELDS:
            c=np.array([h[field] for h in cases['4']['history']]);d=np.array([h[field] for h in cases[label]['history']])
            comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
        checks.append(dict(label=label,orthogonality=orth,linear_residual=P.error,**meta));del model,P,Q,R
passed=prop_pass and spatial_pass and all(v.get('coefficient',1)<.02 and v.get('outer',1)<.002 and v.get('quadrature',1)<.002 for v in comparisons.values())
result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,passed=passed,propagation_passed=prop_pass,
    spatial_passed=spatial_pass,propagation=propagation,comparisons=comparisons,checks=checks,paths=list(cases),
    seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
task.write(OUT/'result.json',result);signal.alarm(0);print('FACTOR RESULT',json.dumps(result),flush=True)
