"""Counterexample candidate: transfer the FAILED source discrepancy to GR.

No admission of206or new physical evolution. Isolate source-state and source
interpolation effects using its two saved histories and the original operator.
"""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import return_resolved_joint_history as prior

OUT=Path('native-early-gr207-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,measure=600)


def prepare():
    assert not OUT.exists();OUT.mkdir();assert not read(OLD/'original-source-norm.json')['passed']
    for name in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','gr']:(OUT/name).mkdir(parents=True,exist_ok=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/name for name in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for p in files:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    for n in [64,128]:
        p=OLD/'gr'/f'source-{n}.npz';os.link(p,OUT/'gr'/p.name);files.append(p)
    files += [OLD/'original-source-norm.json',OLD/'sources.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Determine how the measured early same-solution source error transfers through the original compact GR operator, separating stored source-state differences from interpolation between the same existing edges.',
        decision='If temporal representation dominates, repair that actual input before any coupled return. If state differences remain dominant, the high-component evolution needs repair first.206source failure is never reclassified by this diagnostic.',
        method='Read206coarse and fine sources. Make a reduced fine source at the original coarse edges, preserving its own source and cumulative ports. Apply the identical retarded GR operator to all three. At common times, coarse-minus-fine equals coarse-minus-reduced plus reduced-minus-fine. No new fluid trajectory or charge superposition.',
        controls='Original4/8quadrature and independent direct GR control; exact decomposition at common times. The2percent comparison is diagnostic and cannot override the failed source gate.',
        scope='Original188T/64only. Compact scalar response is not infinity-normalized final charge, self-GR, uniform time error or EOS certification.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,forecast='206source already saved.186three characteristic field reads took6.22s. Four existing-grid GR reads plus setup expected10..60s; allow10minutes. No new physical solve.',
        stop='Input binding, quadrature, independent direct integral or decomposition failure; no new grid/period or altered206verdict.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',prior.base.gr.check())


def measure():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    a=dict(np.load(OUT/'gr/source-64.npz'));b=dict(np.load(OUT/'gr/source-128.npz'))
    ids=np.array([np.argmin(abs(b['t']-t)) for t in a['t']]);assert np.max(abs(b['t'][ids]-a['t']))<1e-18
    reduced={k:(v[ids] if v.ndim and len(v)==len(b['t']) else v.copy()) for k,v in b.items()}
    # 0 labels a stored-data projection, never an extra physical time clock.
    np.savez_compressed(OUT/'gr/source-0.npz',**reduced)
    m=prior.base.gr.Response();fn=FunctionType(m.run.__func__.__code__,dict(m.run.__func__.__globals__,OUT=OUT/'gr'))
    rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(0,8),(128,4)]]
    fine=dict(np.load(OUT/'gr/fields-128-g8.npz'));coarse=dict(np.load(OUT/'gr/fields-64-g8.npz'));projected=dict(np.load(OUT/'gr/fields-0-g8.npz'));low=dict(np.load(OUT/'gr/fields-128-g4.npz'))
    controls=dict(quadrature=prior.aligned(low,fine,'U'));direct,_=prior.base.charge.independent.direct(m,b,8)
    controls['independent_GR']=abs(direct-rows[0]['endpoint_direct'])/max(abs(direct),1e-290)
    result={}
    for key in ['U','U_t','U_x']:
        f=fine[key][ids];d=coarse[key]-f;state=coarse[key]-projected[key];representation=projected[key]-f
        scale=max(float(np.max(abs(f))),1e-290);difference=max(float(np.max(abs(d))),1e-290)
        result[key]=dict(total_relative=float(np.max(abs(d))/scale),state_relative=float(np.max(abs(state))/scale),representation_relative=float(np.max(abs(representation))/scale),decomposition_relative=float(np.max(abs(d-state-representation))/difference))
    endpoints=[r['endpoint_compact_with_metric'] for r in rows[:3]]
    row=dict(classification='Counterexample candidate',diagnostic_controls_passed=controls['quadrature']<.002 and controls['independent_GR']<1e-9 and max(v['decomposition_relative'] for v in result.values())<1e-12,
        controls=controls,common_time_decomposition=result,endpoint_compact_values=endpoints,rows=rows,
        source_admission_passed=False,GR_return_evolution_admitted=False,new_physical_steps=0,physical_horizon_seconds=float(b['t'][-1]),final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',row);print(json.dumps(row),flush=True);assert row['diagnostic_controls_passed'],row


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
