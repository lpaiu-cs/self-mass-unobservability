"""Apply the accepted dense same-history GR input to actual coupled stages.

Counterexample candidate: one high/low feedback iterate on the existing T/8
prefix. No separately evolved charge addition or full nonlinear closure claim.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import read_complete_radau_history as bridge
import evolve_same_solution_gr_return as original

OUT=Path('native-dense-return225-work');INPUT=Path('native-complete-radau224-check-work')
geometry=original.prior;owner=original.owner;joint=original.joint
LD,read,write,sha=original.LD,original.read,original.write,original.sha
CAPS=dict(prepare=180,metric=900,coarse=1800,fine=3600,audit=180)
saved=lambda n:INPUT/'input'/f'material-{n}.npz'


def prepare():
    assert read(INPUT/'result.json')['GR_return_admitted']
    assert read(INPUT/'sources.json')['representation_controls_passed'] and read(INPUT/'driver-polynomial-audit.json')['passed']
    assert not OUT.exists();OUT.mkdir();files=[]
    for folder in ['sweep-0','sweep-1/photons','sweep-1/material','metric','gr']:(OUT/folder).mkdir(parents=True)
    for src in list((INPUT/'sweep-0').rglob('*.npz'))+list((INPUT/'gr').glob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(INPUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    files += [saved(n) for n in [64,128]]+[INPUT/n for n in ['result.json','sources.json','fields-receipt.json','source-receipt.json','driver-polynomial-audit.json','tested-consumer.py']]
    files += [original.OUT/'result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the accepted224dense GR field and its own actual angular packets to the same coupled matter/photon equations as a high/low feedback iterate over the original acceptedT/8prefix.',
        decision='Only original metric, whole stage, physical moment, native anchor, constitutive, balance, branch and2percent paired-time gates admit this iterate. Read the original failedT/64subinterval separately; a longer global norm never erases188failure.',
        model='Reuse188compensated same-stage equations and actual high/low branch choices. Restore exact conserved gas coordinates. High component is the actual saved joint solution; low component is its same-equation increment, never an independent zero-state solution added to a charge. Reject nonzero high-branch crossings. No converged self-GR or fully nonlinear EOS claim.',
        scope='Same531cells,64/128clocks, original front split rule,15/29actualT/8steps. Fixed fine224returned metric drives both clocks. Actual scalar exterior completion, continuum errors and full period remain open.',
        gates=dict(metric_time=.02,metric_quadrature=.002,ray=1e-10,packet=1e-12,stage=1e-12,physical_stage=1e-13,native_anchor=1e-12,constitutive=.002,balance=1e-8,time=.02,max_Newton=3,max_branch_passes=8),
        budgets=CAPS,CPU_threads=1,virtual_GiB=8,
        forecast='188two/fouractualreturnsteps57/97seconds.15/29steps roughly6..10/12..20minutes if marginal cost holds; later branch costs unmeasured. Allow30/60minutes plus15minutes metric. Preserve accepted checkpoints and rejected-stage proposals; no automatic new grid/period or relaxed gate.',
        stop='Any original gate, nonzero branch crossing or cap.225owned paths run sequentially and stop at failure.218,223,224controllers and their source/plan stay unchanged.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(original.OUT/'symbolic.json'))


def initialize_metric():
    geometry.prior.saved=saved
    FunctionType(bridge.endpoint.initialize.__code__,dict(bridge.endpoint.initialize.__globals__,OUT=OUT))()


def metric():
    initialize_metric();m=geometry.Lapse()
    scope=dict(geometry.metric.Lapse.run.__globals__,OUT=OUT/'metric',
        prior=SimpleNamespace(OUT=OUT/'gr',new=SimpleNamespace(OUT=OUT/'gr'),centers=geometry.constraints.centers))
    fn=FunctionType(geometry.metric.Lapse.run.__code__,scope);rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'));coarse=dict(np.load(OUT/'metric/metric-64-g8.npz'));low=dict(np.load(OUT/'metric/metric-128-g4.npz'))
    aligned=bridge.endpoint.aligned
    controls={name:{k:aligned(other,fine,k) for k in geometry.KEYS if k!='delta_lambda_rate'} for name,other in [('time',coarse),('quadrature',low)]}
    for name,other in [('time',coarse),('quadrature',low)]:controls[name]['delta_lambda_rate_integral']=aligned(other,fine,'delta_lambda')
    result=dict(classification='Counterexample candidate',passed=max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002 and max(r['ray_invariant'] for r in rows)<1e-10,
        controls=controls,rows=rows,same_energy_and_actual_packets=True,interval_derivative_pointwise_convergence_unadjudicated=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'metric-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def initialize():
    geometry.OUT=OUT;geometry.prior.saved=saved
    s=inspect.getsource(original.initialize)
    old="return k,np.column_stack([(q[2]-self.kappa*q[0])/self.eu,q[3]/self.nu,q[0]/self.bu,q[1]/self.su])"
    assert s.count(old)==1;s=s.replace(old,'return k,restored_gas(self,q)')
    ns=dict(original.initialize.__globals__,OUT=OUT,restored_gas=bridge.recovery.prior.restored_gas)
    exec(compile(s,__file__,'exec'),ns);ns['initialize']();(OUT/'expanded-increment-initialize.py').write_text(s)
    return ns['Model']


def evolve(n):
    assert read(OUT/'metric-result.json')['passed'];Model=initialize();m=Model(n)
    row=m.run(n,f'return-{n}',n//8)
    p=OUT/f'sweep-1/photons/return-{n}.npz';z=np.load(p);anchor=np.load(saved(n))
    assert np.array_equal(z['actual_step_edges'],anchor['actual_step_edges']),('Original physical clocks changed',n)
    assert np.array_equal(z['joint_stage_times'],anchor['joint_stage_times'])
    audit,_,_,_=geometry.prior.run.verify(p,p)
    row.update(audit=audit,anchor_checks=m.anchor_checks,branch_checks=m.branch_checks,
        maximum_true_stage=max(v[-1]['relative'] for v in m.newton_iterations),maximum_true_physical_stage=max(max(v[-1]['moments']) for v in m.newton_iterations),
        same_saved_stage_equation=True,actual_return_time_evolved=True,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'run-{n}.json',row);assert row['passed']
    write(OUT/f'checks-{n}.json',dict(classification='Counterexample candidate',newton=m.newton_iterations,stages=m.stage_log))


def audit():
    # Original comparison is retained, including its independent stage ledgers.
    FunctionType(original.audit.__code__,dict(original.audit.__globals__,OUT=OUT))()
    r=read(OUT/'result.json');first=[]
    for n,count in [(64,2),(128,4)]:
        z=np.load(OUT/f'sweep-1/photons/return-{n}.npz');first.append(z['joint_stage_conserved_scaled'][2*count-1])
    errors=(np.sum(abs(first[0]-first[1]),axis=-1)/np.maximum(np.sum(abs(first[1]),axis=-1),LD('1e-290'))).astype(float).tolist()
    r.update(original188failure_preserved=True,original_T64_endpoint_conserved_relative=errors,
        limited_to_one_feedback_iterate=True,source_and_metric='same accepted224T/8trajectory')
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
