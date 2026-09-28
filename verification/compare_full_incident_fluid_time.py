"""Counterexample candidate: original64/128 comparison for the full-input solve.

Reuse181coarse bytes; evolve only128 over the same T/32. No extra resolution.
"""
from pathlib import Path
import json,os,resource,sys,time
import numpy as np
import couple_full_incident_fluid as prior

OUT=Path('native-full-time182-work');BEFORE=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha;LD=prior.LD
CAPS=dict(prepare=15,run=250,audit=15)
CAPS['retry']=245


def prepare():
    assert not OUT.exists();OUT.mkdir();reuse={}
    coarse=read(BEFORE/'pilot-result.json');assert coarse['passed'] and read(BEFORE/'saved-audit.json')['passed']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    files=list((BEFORE/'sweep-0').rglob('*.npz'))
    files.extend(BEFORE/name for name in ['normalization.json','photon-conservation-plan.json','check-result.json',
        'sweep-1/photons/pilot-64.npz','sweep-1/photons/pilot-64.json'])
    for p in files:
        dst=OUT/p.relative_to(BEFORE);os.link(p,dst);reuse[str(dst.relative_to(OUT))]=dict(path=str(p),sha256=sha(p))
    write(OUT/'reuse.json',reuse)
    modules=[Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None)
        and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    files+=modules+[Path(__file__),BEFORE/'pilot-receipt.json',BEFORE/'saved-audit.json',prior.drive.FIELDS/'born-g8.npz']
    forecast=2*coarse['stepping_seconds']+coarse['operator_point_seconds']+read(BEFORE/'pilot-receipt.json')['seconds']-coarse['seconds']+20
    assert forecast<CAPS['run']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='9c72e5190',
        previous_turn='Progress: actual full incident geometry was coupled to photon/B/S/Etilde/H stages and the original64T/32path passed true equations and same-solution ledgers. Its own GR and final charge remain open.',
        claim='Determine whether all original six photon and four material channels remain within2percent between64 and128 over the SAME full-input interval before continuing toward GR.',
        decision='A failed time channel prevents full-horizon admission. A pass warrants a separately bounded same-solution extension/GR return, not a final charge conclusion.',
        reuse='Hard-link accepted64state, actual stages and ledgers byte-for-byte. Reuse EOS/background/incident data. Only the missing original128four-macro T/32path is new. No accepted64replay, new clock, source, splitting rule, EOS or nonlinear background.',
        equations='Unchanged181full incident source and geometric native operator, direct physical forcing, four simultaneous material unknowns, max3Newton, same original physical and constitutive gates. Previous gross-recovery derivative failure remains unresolved and is not relabeled as passed.',
        forecast=dict(seconds=forecast,assumed_range_seconds=[170,250],hard_cap_seconds=250,
            basis='Twice measured coarse stepping97.678s plus measured point/setup overhead and20s reserve. Fine-step Krylov/branch changes remain an assumption.'),
        gates=dict(time=.02,stage=1e-12,physical_stage=1e-13,conservation=1e-8,constitutive=.002,port=1e-12,max_Newton_solves=3),
        budget=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Gate, iteration or cap stops this path. Preserve failure. No automatic finer clock, split, waveform, method, full horizon or GR dispatch.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    # Preserve the exact symbolic identities used by this unchanged equation.
    for name in ['symbolic.json','geometric-face-symbolic.json']:
        p=BEFORE/name;assert read(p)['passed'];os.link(p,OUT/name)


def compare():
    paths=[OUT/f'sweep-1/photons/pilot-{n}.npz' for n in [64,128]]
    rows=[];photons=[];material=[]
    for path in paths:
        rows.append(read(path.with_suffix('.json')));assert rows[-1]['passed']
        with np.load(path) as p:
            photons.append(p['moments'][:,[0,1,2,3,5,6]]);material.append(p['conserved_material_history'])
            if len(rows)==1:times=p['t'].copy()
            else:assert np.array_equal(times,p['t'])
    ph=prior.joint.previous.run.c.relative(*photons);mh=prior.joint.previous.run.c.relative(*material)
    result=dict(classification='Counterexample candidate',passed=max(ph+mh)<.02,photon_time=ph,material_time=mh,rows=rows,
        same_horizon_seconds=float(times[-1]),coarse_bytes_reused=True,full_incident_input_applied=True,
        self_GR_return_closed=False,full_horizon_completed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'time-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def run():
    prior.OUT=OUT;prior.initialize()
    pilot,source=prior.clone(prior.pilot,[('initialize();m=Model(64)','m=Model(128)'),
        ("m.run(64,'pilot-64',2)","m.run(128,'pilot-128',4)"),
        ("paths(1)[0]/'pilot-64.npz'","paths(1)[0]/'pilot-128.npz'")])
    (OUT/'expanded-fine-dispatch.py').write_text(source);pilot();compare()


def audit():
    rows=[]
    for rel,v in read(OUT/'reuse.json').items():assert sha(OUT/rel)==sha(v['path'])==v['sha256'],rel
    for n in [64,128]:
        file=OUT/f'sweep-1/photons/pilot-{n}.npz';p=dict(np.load(file));w=p['joint_stage_weights'].astype(LD)
        assert np.array_equal(w,p['accepted_angular_quadrature_weights'])
        assert np.array_equal(p['joint_stage_times'],p['accepted_angular_times'])
        q=p['conserved_material_history'][-1]/prior.AMP
        final=np.column_stack([q[2]-p['energy_offset_reference'][-1]/prior.AMP,q[3],q[0],q[1]])
        expected=np.sum(w[:,None,None]*(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled']),axis=0,dtype=LD)
        actual=final+p['material_floor_discard_scaled']
        err=np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))
        port=prior.joint.previous.run.packets(file)[2];assert max(err)<1e-8 and port<1e-12
        rows.append(dict(clock=n,local_material_balance=err.astype(float).tolist(),angular_port=port,actual_stages=len(w)))
    write(OUT/'audit-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        final_charge_conclusion='unadjudicated',time_verdict=read(OUT/'time-result.json')['passed']))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                target=OUT/'initial-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(target)==h,p
            if action=='retry':
                failed=read(OUT/'run-receipt.json');assert failed['seconds']+CAPS[action]<250
                write(OUT/'dispatch-repair.json',dict(classification='Counterexample candidate',failure=failed,
                    repair='Initialize the original owner before cloning its pilot globals so the newly constructed Model is bound. Remove the duplicate initialize call. No Model instance or physical step existed in the failed dispatch.',
                    cap=CAPS[action],source_sha256=sha(__file__)))
            elif action=='audit':assert sha(__file__)==read(OUT/'dispatch-repair.json')['source_sha256']
        if action=='retry':run()
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
