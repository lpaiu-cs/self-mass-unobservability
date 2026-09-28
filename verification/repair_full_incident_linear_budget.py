"""Counterexample candidate: bounded refinement of the same final interval.

Only the linear iterative-refinement budget changes, from four to seven solves.
The physical RHS, stage equations, Newton limit, time grid and gates are frozen.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,resource,sys,time
import numpy as np
import finish_full_incident_horizon as capture

OUT=capture.OUT;prior=capture.prior;owner=prior.owner;joint=owner.joint
read,write,sha=capture.read,capture.write,capture.sha
CAPS=dict(admit=15,coarse=450,fine=700,audit=25)


def admit():
    assert not (OUT/'repair-plan.json').exists()
    failure=read(OUT/'linear-failure.json');trace=read(OUT/'linear-calls.json')['calls'][failure['call_begin']:failure['call_end']]
    assert len(trace)==4 and trace[-1]['true_relative']<1e-3
    assert read(OUT/'probe-receipt.json')['error']==failure['error']
    files=[Path(__file__),Path(capture.__file__),OUT/'linear-failure.json',OUT/'linear-calls.json',OUT/'probe-receipt.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'repair-plan.json',dict(classification='Conjectural',checkpoint='7fa047343',
        failure='The original failed linear system was reproduced exactly. Four float64 Krylov solves with long-double residual refinement stop at7.22270554e-8 vector and1.40074188e-4 physical error. The last correction reduces its own true residual by0.000394328; it has not reached a precision plateau.',
        repair='Allow at most seven total inner solves instead of four, keeping each original GMRES restart20/maxiter5, extended residual evaluation, rtol, physical scales, three Newton proposals and all physical gates. No preconditioner, method, source, grid, amplitude or horizon change.',
        prediction='If the last measured correction factor persists, three further corrections reduce the worst physical residual to about8.6e-15. This is a hypothesis and must be checked against the actual equations.',
        decision='Finish coarse only from accepted60/64 to64, then fine only from120/128 to128 if coarse passes. Require actual unchanged-prefix, stage, balance, packet and original2percent full-history comparisons before accepting the full horizon.',
        reuse='All accepted15/16history is restored exactly; no accepted prefix replay. Original185and reproduction failures remain failed. Only last8/16existing substeps are new.',
        gates=dict(linear_vector=1e-14,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02,max_Newton=3,max_inner_solves=7),
        forecast=dict(coarse_seconds=[200,450],fine_seconds=[300,700],
            basis='Original terminal attempt91.63s; reproduction92.88s. Each final correction10.5s, so three additional corrections cost about31.5s at the observed system. Previous interval8/16substeps cost136/259s before constructor overhead. Subsequent correction rates are unmeasured assumptions.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Any original physical/time/prefix gate, seven-solve limit or450/700s cap stops. No further solver, grid, horizon, tolerance or iteration ladder; no automatic failed-path restart.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def initialize():
    prior.prior.OUT=OUT;prior.prior.initialize()
    source=inspect.getsource(joint.solve)
    for old,new in [('range(4)','range(7)'),('assert k<3,','assert k<6,')]:
        assert source.count(old)==1;source=source.replace(old,new)
    ns=dict(joint.solve.__globals__);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-seven-solve.py').write_text(source)
    stage=owner.Model.run.__globals__['stages'];stage=FunctionType(stage.__code__,dict(stage.__globals__,solve=ns['solve']))
    run=owner.Model.run;owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=stage),argdefs=run.__defaults__)


def evolve(n):
    if n==128:assert read(OUT/'repair-64.json')['passed']
    initialize();m=owner.Model(n);start=time.monotonic()
    row=m.run(n,f'complete-{n}',n,restart=f'interval-15-{n}')
    file=OUT/f'sweep-1/photons/complete-{n}.npz';prefix=OUT/f'sweep-1/photons/interval-15-{n}.npz'
    audit,_,_,_=prior.verify(file,prefix)
    row.update(audit=audit,continuation_seconds=time.monotonic()-start,accepted15of16preserved=True,
        maximum_new_inner_solves=1+max(v['corrections'] for v in joint.radau.prior.owner.reuse.LINEAR),
        self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'repair-{n}.json',row);assert row['passed']


def audit():
    rows=[read(OUT/f'repair-{n}.json') for n in [64,128]];arrays=[];clocks=[]
    for n in [64,128]:
        report,ph,gas,clock=prior.verify(OUT/f'sweep-1/photons/complete-{n}.npz',OUT/f'sweep-1/photons/interval-15-{n}.npz')
        arrays.append((ph,gas));clocks.append(clock)
    assert np.array_equal(*clocks)
    relative=joint.previous.run.c.relative;errors=relative(arrays[0][0],arrays[1][0])+relative(arrays[0][1],arrays[1][1])
    for k in range(1,16):assert read(capture.OLD/f'comparison-{k:02d}.json')['passed']
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02,rows=rows,time_relative=errors,
        full_horizon_seconds=float(clocks[0][-1]),full_horizon_photon_material_completed=True,
        original_failure_preserved=True,accepted15of16preserved=True,physical_equation_and_gates_unchanged=True,
        maximum_inner_solves_changed_from=4,maximum_inner_solves_changed_to=7,
        self_GR_return_closed=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    assert result['passed'],errors


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'repair-{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='admit':
            for p,h in read(OUT/'repair-plan.json')['bindings'].items():assert sha(p)==h,p
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
