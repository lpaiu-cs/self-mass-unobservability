"""Join actual fine photon captures and read the same complete material source.

Counterexample candidate. This consumes the existing two physical paths only.
It does not integrate a new path or claim a closed GR return or final charge.
"""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import read_complete_radau_history as prior
import complete_coarse_photon_history as coarse

OUT=Path('native-full-captured244-work')
FULL=Path('native-common-arithmetic239-work')
COMPLETE=Path('native-preimage-photon243-work')
INPUT=OUT
read,write,sha,LD,AMP=prior.read,prior.write,prior.sha,prior.LD,prior.AMP
CAPS=dict(check=600,prepare=600,assemble=600,endpoints=2700,source=5400)
saved=lambda n:FULL/f'sweep-1/photons/complete-{n}.npz'


def bind(fn,**values):
    return FunctionType(fn.__code__,dict(fn.__globals__,**values),argdefs=fn.__defaults__)


def joined(z,prefix,captures):
    d={k:np.array(v,copy=True) for k,v in prefix.items()};start=len(d['times'])
    assert start+len(captures)==len(z['joint_stage_times'])
    for key,other in [('times','joint_stage_times'),('weights','joint_stage_weights')]:
        assert np.array_equal(d[key],z[other][:start]),key
    for i,c in enumerate(captures,start):
        assert c['time']==z['joint_stage_times'][i] and c['weight']==z['joint_stage_weights'][i]
        assert np.array_equal(c['collision_rates'],z['joint_collision_rates_scaled'][i])
    for key in ['photon_moments','collision_rates','radial_ports']:
        d[key]=np.concatenate([d[key],np.array([c[key] for c in captures])])
    d['angular']=np.concatenate([d['angular'],z['accepted_angular_luminosity'][start:]])
    d['times']=z['joint_stage_times'];d['weights']=z['joint_stage_weights']
    d['endpoint_occupation']=z['photon_history_scaled_occupation'][-1]
    for key in ['times','weights','photon_moments','collision_rates','radial_ports','angular']:
        assert np.array_equal(d[key][:start],prefix[key]),key
    port=np.sum(d['weights'][:,None,None]*d['radial_ports'],axis=0,dtype=LD)
    radial=float(np.max(abs(port-z['radial_ports'][-1])/np.maximum(abs(z['radial_ports'][-1]),LD('1e-290'))))
    angular=float(np.sum(abs(d['angular']-z['accepted_angular_luminosity']))/max(np.sum(abs(z['accepted_angular_luminosity'])),LD('1e-290')))
    q=z['conserved_material_history'][-1]/AMP;k=z['energy_offset_reference'][-1]/AMP
    actual=np.column_stack([q[2]-k,q[3],q[0],q[1]])+z['material_floor_discard_history_scaled'][-1]
    target=np.sum(d['weights'][:,None,None]*(z['joint_native_rates_scaled']+d['collision_rates']),axis=0,dtype=LD)
    ledger=(np.sum(abs(actual-target),axis=0)/np.maximum(np.sum(abs(actual)+abs(target),axis=0),LD('1e-290'))).astype(float).tolist()
    row=dict(classification='Counterexample candidate',passed=radial<1e-12 and angular<1e-12 and max(ledger)<1e-8,
        reused_prefix_steps=start//2,actual_captured_tail_steps=len(captures)//2,
        radial_port_relative=radial,angular_relative=angular,material_ledger=ledger,
        new_photon_solves=0,new_material_steps=0,final_charge_conclusion='unadjudicated')
    assert row['passed'],row
    return d,row


def check():
    assert read(COMPLETE/'result.json')['passed']
    z=dict(np.load(coarse.saved(64)));prefix=dict(np.load(COMPLETE/'conditional-recovered-64.npz'))
    captures=[dict(np.load(folder/f'captured-64-{i:03d}.npz'))
              for folder,ids in [(coarse.CAPTURE,[234,235]),(coarse.FULL,[236,237])] for i in ids]
    actual,row=joined(z,prefix,captures);expected=dict(np.load(COMPLETE/'recovered-64.npz'))
    assert actual.keys()==expected.keys()
    for key in expected:assert np.array_equal(actual[key],expected[key]),key
    # A displaced capture must fail before the source consumes it.
    bad=[dict(c) for c in captures];bad[0]['time']=np.nextafter(bad[0]['time'],np.inf)
    try:joined(z,prefix,bad)
    except AssertionError:rejected=True
    else:raise AssertionError('Wrong stage time was accepted')
    write(OUT/'regression.json',dict(classification='Counterexample candidate',passed=True,
        every_completed_coarse_array_exact=True,wrong_capture_time_rejected=rejected,coarse_controls=row,
        existing_consumer=prior.check()))


def prepare():
    assert read(OUT/'regression.json')['passed']
    assert read(FULL/'controller-status.json')['state']=='completed'
    result=read(FULL/'result.json');assert result['passed']
    assert result['same_corrected_native_arithmetic'] and result['common_saved_prefix_material_ledger_revalidated']
    assert read(FULL/'audit-receipt.json')['error'] is None
    assert read(FULL/'audit-receipt.json')['source_sha256']==sha(Path('verification/complete_common_arithmetic_fine.py'))
    assert read(COMPLETE/'controller-status.json')['state']=='completed'
    assert not (OUT/'plan.json').exists();files=[]
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr']:(OUT/part).mkdir(parents=True,exist_ok=True)
    for src in list((FULL/'sweep-0').rglob('*.npz'))+[FULL/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(FULL);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for n,count in [(64,119),(128,231)]:
        z=np.load(saved(n));assert len(z['actual_step_edges'])==count+1
        assert read(FULL/f'path-{n}.json')['passed'];files += [saved(n),FULL/f'path-{n}.json']
    z=np.load(saved(128));original=np.load(coarse.archive.saved(128))
    for key in ['joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled']:
        assert np.array_equal(z[key][:430],original[key]),key
    assert abs(float(z['actual_step_edges'][-1])-float(np.load(saved(64))['actual_step_edges'][-1]))<1e-18
    files += [FULL/f'captured-128-{i:03d}.npz' for i in range(430,462)]
    files += [FULL/n for n in ['result.json','controller-status.json','fine-receipt.json','audit-receipt.json','prefix-result.json']]
    files += [COMPLETE/f'recovered-64.{ext}' for ext in ['npz','json']]
    files += [COMPLETE/'result.json',COMPLETE/'assemble-receipt.json',coarse.recovery.OUT/'recovered-128.npz',coarse.recovery.OUT/'recovered-128.json',coarse.archive.saved(128)]
    files += [prior.prior.OUT/n for n in ['expanded-source.py','expanded-fields.py']]
    files += [OUT/'regression.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='After actual239paired-time admission, supply its complete119/231same-solution material/photon/port histories to the existing driver-aware GR source reader.',
        reuse='Use243complete coarse photons and223fine215prefix plus239actual final16step captures. No conditional photon solve, new material step, replayed prefix, finer grid or new parameter path.',
        decision='Original source-time2percent, pressure0.2percent, mapping/dense/polynomial1e-12 and ledger1e-8 gates decide full-period source admission. Field and actual GR return are separate required steps; this source is the primary trajectory, not a substitute for236returned-solution charge.',
        gates=read(prior.OUT/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,CPU_affinity=2,virtual_GiB=16,
        forecast='224measured endpoints185.92s andsource570.36s for326steps.350steps is7.4percent larger; roughly4..8and11..20minutes if per-step costs hold, later branches unmeasured. Allow45and90minutes,16GiB, without relaxing numerical gates.',
        stop='Dependency failure, source SHA change, prefix/capture mismatch or any original representation gate; retain a failed source-time comparison without asserting physical acceptance. No automatic GR-field dispatch.',
        original_return_horizon_seconds=float(z['actual_step_edges'][-1]),full_declared_period=True,
        self_GR_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    write(OUT/'symbolic.json',prior.prior.symbolic())


def assemble():
    for ext in ['npz','json']:os.link(COMPLETE/f'recovered-64.{ext}',OUT/f'recovered-64.{ext}')
    z=dict(np.load(saved(128)));prefix=dict(np.load(coarse.recovery.OUT/'recovered-128.npz'))
    captures=[dict(np.load(FULL/f'captured-128-{i:03d}.npz')) for i in range(430,462)]
    d,row=joined(z,prefix,captures)
    row.update(clock=128,recovered_steps=231,horizon_seconds=float(z['actual_step_edges'][-1]),paired_time_admitted=True)
    np.savez_compressed(OUT/'recovered-128.npz',**d);write(OUT/'recovered-128.json',row)
    write(OUT/'assembly-result.json',dict(classification='Counterexample candidate',passed=True,
        complete_same_solution_photon_histories=True,rows=[read(OUT/'recovered-64.json'),row],
        physical_period_seconds=row['horizon_seconds'],new_material_steps=0,new_photon_solves=0,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))


recovered=bind(prior.recovered,INPUT=OUT)
endpoints=bind(prior.endpoints,OUT=OUT,saved=saved,recovered=recovered)


def source():
    bind(prior.source,OUT=OUT,INPUT=OUT,saved=saved)()
    r=read(OUT/'sources.json')
    r.update(full_declared_period=True,missing_final_interval=False,new_material_steps=0,
        actual_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;OUT.mkdir(exist_ok=True)
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{2});resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    prior.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action not in ['check','prepare']:
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
