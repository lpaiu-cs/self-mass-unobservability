"""Finish the preregistered near-radial direct ray after the phase-263 one-hour stop.

Phase 263 completed and saved packet 0 (near-grazing) before its registered
3600 s budget stopped packet 31 (near-radial). Its bindings are unchanged, so
packet 0 is reused as saved. Only packet 31 is computed, with the same code,
ODE tolerances, cohort and 0.2 percent gate. Two rays remain a representative
trajectory control, not a continuum certificate.
"""
from pathlib import Path
import json,resource,time
import numpy as np
import propagate_primitive_characteristics as current
b=current.b;p=current.p;root=current.OUT;old=root/'direct-reference';out=root/'direct-reference-264'
read,write,sha=current.read,current.write,current.sha
BUDGET=21600
start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(BUDGET)
try:
    assert not out.exists();out.mkdir()
    registered=read(old/'plan.json');stopped=read(old/'receipt.json');saved=read(old/'progress.json')
    assert 'TimeoutError' in stopped['error'] and not (old/'result.json').exists()
    for q,h in registered['bindings'].items():assert sha(q)==h,q
    assert registered['packets']==[0,31] and registered['emission_cell']==0
    assert [r['packet_index'] for r in saved['rows']]==[0] and read(old/'packet-0.json')==saved['rows'][0]
    write(out/'plan.json',dict(classification='Conjectural',packets=registered['packets'],emission_cell=registered['emission_cell'],
        claim=registered['claim'],gate=registered['gate'],rtol=registered['rtol'],atol=registered['atol'],maximum_seconds=BUDGET,
        scope=registered['scope'],
        continuation='Phase263 stopped at its registered 3600 s budget inside packet 31 after packet 0 had completed and been saved. Reuse that identical packet 0 under the unchanged registered bindings and compute only packet 31.',
        budget_reason='Measured packet 0: 1820.997 s. Packet 31 ran about 1771 s without completing, so its total cost is unmeasured. Allow 6 hours for this ray alone; if it does not finish, stop and re-plan instead of extending.',
        original_receipt=stopped,reused={str(q):sha(q) for q in [old/'packet-0.npz',old/'packet-0.json']},
        bindings={str(q):sha(q) for q in [Path(__file__),Path(current.__file__),Path(b.__file__),Path(p.__file__),
            root/'plan.json',root/'check.json',root/'pilot-0.npz',old/'plan.json',old/'receipt.json',old/'progress.json']}))
    p.OUT=out; p.Metric=b.Metric; b.Moments=current.base.previous.Moments
    p.initialize();m=p.Photons(8,8);original=m.cohorts;index=31
    def cohorts(now,order,cells=None):
        values=original(now,order,[0]);return tuple(v[index:index+1] for v in values)
    m.cohorts=cohorts
    # Observation only: the same rhs values reach the same solver.
    solve=p.solve_ivp;live=dict(rhs_evaluations=0,s=0.,written=time.monotonic())
    def observed(fun,span,y0,**options):
        def rhs(s,y):
            live['rhs_evaluations']+=1;live['s']=max(live['s'],float(s));now=time.monotonic()
            if now-live['written']>60:
                live['written']=now
                write(out/'live-progress.json',dict(packet_index=index,s_maximum_evaluated=live['s'],
                    rhs_evaluations=live['rhs_evaluations'],seconds=now-start))
            return fun(s,y)
        return solve(rhs,span,y0,**options)
    p.solve_ivp=observed
    try:z,row=m.propagate(m.d.T,8,[0])
    finally:p.solve_ivp=solve
    row['packet_index']=index;row['work_identity_independent']=True
    np.savez_compressed(out/f'packet-{index}.npz',**z);write(out/f'packet-{index}.json',row)
    rows=[saved['rows'][0],row];write(out/'progress.json',dict(rows=rows,reused_packet_indices=[0]))
    fine=np.load(root/'pilot-0.npz');errors={}
    for packet,folder in [(0,old),(31,out)]:
        direct=np.load(folder/f'packet-{packet}.npz')
        for key in ['emission_t','owner','background_packet_energy_erg','radius_cm','direction']:
            assert np.array_equal(direct[key],fine[key][packet:packet+1]),(packet,key)
        for key in ['delta_radius_cm','delta_direction','delta_log_H','delta_arrival_seconds','integrated_log_H_work']:
            a=direct[key];v=fine[key][packet:packet+1]
            errors[f'{packet}:{key}']=float(np.max(abs(a-v))/max(np.max(abs(a)),1e-290))
    passed=max(errors.values())<registered['gate']
    write(out/'result.json',dict(classification='Counterexample candidate',passed=passed,relative=errors,rows=rows,
        reused_packet_indices=[0],original_energy_integral_checked=True,physical_final_charge_solved=False));assert passed,errors
except BaseException as exc:error=repr(exc);raise
finally:
    if out.exists():write(out/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
