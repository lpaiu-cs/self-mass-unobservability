"""Read-only intervention: does the native branch history cause the defect?"""
from pathlib import Path
import inspect,json,resource,time
import numpy as np
import resolve_returned_native_time as work
prior=work.prior;OUT=work.OUT;LD=work.LD
read,write,sha=work.read,work.write,work.sha
start=time.monotonic();cap=40
assert not (OUT/'branch-receipt.json').exists()
assert read(OUT/'diagnose-receipt.json')['seconds']+cap<75
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));prior.joint.previous.original.inf.incident.native.deadline(cap)
files=[Path(__file__),Path(work.__file__),Path(prior.__file__),OUT/'diagnosis.json',OUT/'decomposition.npz']
write(OUT/'branch-plan.json',dict(classification='Conjectural',
    claim='Separate changing native directional branches from ordinary smooth time quadrature by freezing the real endpoint branch tape in a read-only RHS intervention.',
    decision='If this removes the baryon quadrature defect, resolve actual branch events before any new returned evolution. A frozen branch is NOT an admitted physical equation.',
    cap_seconds=cap,total_diagnostic_cap_seconds=75,new_physical_steps=0,
    stop='No new evolution, clock, source or criterion change. Stop on owner/decomposition/budget mismatch.',
    bindings={str(p):sha(p) for p in files}))
error=None
try:
    prior.OUT=OUT;prior.initialize();m=prior.Model(128)
    coarse,fine=[dict(np.load(work.OLD/f'sweep-1/photons/return-{n}.npz')) for n in [64,128]]
    t=float(fine['joint_stage_times'][-1]);ref=work.gas(m,fine['joint_stage_conserved_scaled'][-1])
    m.native(t,ref);fn,reset=m.selected[t]
    reference=inspect.getclosurevars(fn.__globals__['select']).nonlocals
    frozen=[v.copy() for v in reference['decisions']];events=reference['events']
    integrals=[];traces=[];rates=[]
    for path in [coarse,fine]:
        rows=[]
        for j,t in enumerate(path['joint_stage_times']):
            t=float(t);g=work.dense(m,coarse,t);reset(False)
            rows.append(m.native(t,g,tangent=fn)*m.units)
            if path is fine:
                m.native(t,g);actual=m.selected[t][0]
                tape=inspect.getclosurevars(actual.__globals__['select']).nonlocals
                changed=[]
                for k,(a,b) in enumerate(zip(tape['decisions'],frozen)):
                    different=a!=b
                    if np.any(different):changed.append(dict(event=k,kind=events[k][0],shape=list(a.shape),
                        count=int(np.count_nonzero(different)),coordinates=np.argwhere(different).tolist()))
                traces.append(dict(time=t,changed=changed))
        integrals.append(np.sum(path['joint_stage_weights'][:,None,None].astype(LD)*np.array(rows),axis=0,dtype=LD))
        rates.append(np.array(rows))
    before=dict(np.load(OUT/'decomposition.npz'));defect=integrals[0][:,2]-integrals[1][:,2]
    norm=np.sum(abs(before['baryon_difference']),dtype=LD)
    result=dict(classification='Counterexample candidate',passed=True,
        frozen_branch_quadrature_over_original_difference=float(np.sum(abs(defect))/norm),
        original_quadrature_over_original_difference=float(np.sum(abs(before['quadrature_and_base_history']))/norm),
        branch_events_per_stage=[sum(v['count'] for v in r['changed']) for r in traces],
        scope='Read-only counterfactual source evaluation; no frozen-branch evolution is accepted.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'branch-result.json',result);write(OUT/'branch-traces.json',dict(classification='Counterexample candidate',traces=traces))
    np.savez_compressed(OUT/'frozen-branch-rates.npz',coarse=rates[0],fine=rates[1],baryon_defect=defect)
    print(json.dumps(result))
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'branch-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__)))
