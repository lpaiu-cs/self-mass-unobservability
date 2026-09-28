"""Restore archived conserved coordinates exactly, then resume photon recovery.

Counterexample candidate: preserve213's native bit-identity and original gates.
No rate fitting: only the stored coordinate map determines the inverse.
"""
from pathlib import Path
from types import FunctionType
import gc,inspect,json,os,resource,sys,time
import numpy as np
import recover_remaining_joint_photons as prior

OUT=Path('native-coordinate-photon214-work');OLD=prior.OUT
base=prior.base;LD=prior.LD;read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,check=600,check_interval=1200,coarse=7200,fine=10800,audit=300)


def restored_gas(m,q):
    g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
    assert np.all(np.isfinite(q)) and np.all(m.units>0)
    # Restore the three independent products first; E then uses exact stored B.
    # Eight neighboring floats is a bounded inversion, not a changed tolerance.
    for component,coordinate in [(2,0),(3,1),(1,3),(0,2)]:
        for attempt in range(9):
            actual=m.conserved(g)[coordinate];mask=actual!=q[coordinate]
            if not np.any(mask):break
            assert attempt<8,('Coordinate inverse has no nearby exact preimage',component,np.flatnonzero(mask).tolist())
            target=np.where(actual[mask]<q[coordinate,mask],LD('inf'),LD('-inf'))
            g[mask,component]=np.nextafter(g[mask,component],target)
    assert np.array_equal(m.conserved(g),q)
    return g


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[];reused={}
    assert read(OLD/'controller-status.json')['state']=='failed'
    probe=read(OLD/'inverse-probe.json');assert probe['preimages'][0]['native_max_difference']=='0.0'
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64','clock-128']:(OUT/folder).mkdir(parents=True)
    for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);files.append(p);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['controller-status.json','fine-receipt.json','native-replay.json','inverse-probe.json','rejected-original-128.npz','expanded-recovery-64.py','expanded-recovery-128.py','restart-check.json','plan.json']]
    for n in [64,128]:files += [OLD/f'accepted-{n}.npz',prior.saved(n)]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='af6bc01aa',
        claim='Continue the SAME213photon history after exactly restoring its archived conserved coordinates. Keep every original213gate, including native bit identity.',
        cause='213fine step9 passed the original full equation5.365e-17 but native S differed by1024 (4ULPs). The stored conserved E round-trip at cell15 differed by0.5. The nearest lower normalized E exactly restores every conserved value and every archived native value. Native(t,g) factors through conserved(g); no rate or photon moment is used to select the inverse.',
        method='Invert B,S,H products and then Eref addition with at most8nearest representable moves per coordinate. Require exact conserved equality. Check ALL658stored stages and exact native identity before dispatch. Resume213accepted6/9steps bit-identically; reuse the failed fine photon pair as an equation-tested proposal. Apply the same inverse to stage gas and the previous initial gas. No native-fluid reintegration.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads_per_path=1,virtual_GiB_per_path=6,maximum_parallel_paths=2,
        forecast='Reuse all accepted6/9steps and the retained fine pair; remaining105/206conditional steps. Earlier20..45s/step suggests35..79min coarse and69..155min fine. Later speed unmeasured. Allow2/3hours as already authorized.',
        stop='Any exact inverse/native identity, unchanged equation/endpoint/ledger/port,12linear corrections or wall cap fails. Preserve213failure and all accepted214checkpoints. Do not expand grid, period, physical paths or retry without a diagnosed correction.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))


def seed(m,z,n):
    p=dict(np.load(OLD/f'accepted-{n}.npz'));begin=int(p['step']);logs=json.loads(str(p['logs']))
    assert begin==len(logs) and begin in [6,9]
    for key in ['moments','collisions','ports','packets']:assert len(p[key])==2*begin
    return begin,p['x'].copy(),logs,list(p['moments']),list(p['collisions']),list(p['packets']),list(p['ports']),None


def interval_starts(z):
    return np.array([int(np.argmin(abs(z['actual_step_edges']-t))) for t in z['t'][:-1]])


def interval_model(z,n,step):
    start=max(v for v in interval_starts(z) if v<=step)
    m=base.base.owner.Model(n)
    # The actual producer recreates its model at each canonical continuation.
    # Replay only prior native calls inside that interval, never an evolution.
    for i in range(2*start,2*step):
        g=restored_gas(m,z['joint_stage_conserved_scaled'][i])
        assert np.array_equal(m.native(z['joint_stage_times'][i],g)*m.units,z['joint_native_rates_scaled'][i]),('Warm interval native',n,i)
    return m


def check():
    import sympy as sp
    E,B,k,u,v=sp.symbols('E B k u v',nonzero=True)
    assert sp.simplify(((E-k*B)/u)*u+k*((B/v)*v)-E)==0
    FunctionType(base.base.prior.initialize.__code__,dict(base.base.prior.initialize.__globals__,OUT=OUT))()
    rows=[]
    for n in [64,128]:
        z=dict(np.load(prior.saved(n)));m=None;changes=0;starts=interval_starts(z);models=0
        for i,(t,q) in enumerate(zip(z['joint_stage_times'],z['joint_stage_conserved_scaled'])):
            if m is None or i in 2*starts:
                if m is not None:del m;gc.collect()
                m=interval_model(z,n,i//2);models+=1
                write(OUT/'check-progress.json',dict(clock=n,stage=i,total=len(z['joint_stage_times']),model_instances=models))
            naive=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
            g=restored_gas(m,q);changes+=np.count_nonzero(naive!=g)
            native=m.native(t,g)*m.units;old=z['joint_native_rates_scaled'][i]
            assert np.array_equal(native,old),('Original native identity',n,i,str(np.max(abs(native-old))))
        begin,x,logs,moments,collisions,packets,ports,_=seed(m,z,n)
        checkpoint=FunctionType(prior.checkpoint.__code__,dict(prior.checkpoint.__globals__,OUT=OUT))
        checkpoint(n,x=x,step=begin,t=z['actual_step_edges'][begin],moments=moments,collisions=collisions,packets=packets,ports=ports,logs=np.array(json.dumps(logs)))
        original=dict(np.load(OLD/f'accepted-{n}.npz'));copied=dict(np.load(OUT/f'accepted-{n}.npz'))
        for key in original:assert np.array_equal(original[key],copied[key]),key
        rows.append(dict(clock=n,stages=len(z['joint_stage_times']),model_instances=models,changed_normalized_coordinates=int(changes),every_conserved_value_exact=True,every_native_value_exact=True,reused_steps=begin,checkpoint_value_identity=True));del m,z;gc.collect()
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Exact-real inverse of B normalization and Eref mapping; finite-precision inverse and original native rates are independently checked at all stored stages.'))
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,new_physical_steps=0,new_photon_solves=0))


def run(n):
    check_source=inspect.getsource(prior.prior.original_equation)
    old="initial=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])"
    assert check_source.count(old)==1;check_source=check_source.replace(old,'initial=restored_gas(m,q)')
    ns=dict(prior.prior.original_equation.__globals__,OUT=OUT/f'clock-{n}',restored_gas=restored_gas);exec(compile(check_source,__file__,'exec'),ns);check_original=ns['original_equation']
    def original_equation(m,z,step,t,h,x,gas,photons,times,cs,ss):
        try:return check_original(m,z,step,t,h,x,gas,photons,times,cs,ss)
        except BaseException:
            np.savez_compressed(OUT/f'rejected-original-{n}.npz',step=step,time=t,step_size=h,x_initial=x,photon_stage_solution=photons,gas=gas,stage_times=times);raise
    source=(OLD/f'expanded-recovery-{n}.py').read_text()
    old="m=base.owner.Model(n);z=dict(np.load(saved(n)))"
    new="z=dict(np.load(saved(n)));m=interval_model(z,n,int(np.load(OLD/f'accepted-{n}.npz')['step']))"
    assert source.count(old)==1;source=source.replace(old,new)
    old='for step in range(begin,stop):'
    new="for step in range(begin,stop):\n        if step!=begin and step in interval_starts(z):\n            units=m.units.copy();scale=m.scale;del m;gc.collect();m=interval_model(z,n,step)\n            assert np.array_equal(m.units,units) and m.scale==scale"
    assert source.count(old)==1;source=source.replace(old,new)
    old='gas.append(np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su]))'
    assert source.count(old)==1;source=source.replace(old,'gas.append(restored_gas(m,q))')
    old='sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]'
    new="sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]\n        if n==128 and step==begin:\n            proposal=dict(np.load(OLD/'rejected-original-128.npz'));assert int(proposal['step'])==step and np.array_equal(proposal['x_initial'],x)\n            sol=proposal['photon_stage_solution'].copy().ravel()"
    assert source.count(old)==1;source=source.replace(old,new)
    checkpoint=FunctionType(prior.checkpoint.__code__,dict(prior.checkpoint.__globals__,OUT=OUT))
    snapshot=FunctionType(prior.snapshot.__code__,dict(prior.snapshot.__globals__,OUT=OUT))
    ns=dict(base.run.__globals__,OUT=OUT,OLD=OLD,saved=prior.saved,seed_folder=prior.seed_folder,seed=seed,checkpoint=checkpoint,snapshot=snapshot,restored_gas=restored_gas,original_equation=original_equation,interval_model=interval_model,interval_starts=interval_starts,gc=gc)
    exec(compile(source,__file__,'exec'),ns);(OUT/f'expanded-recovery-{n}.py').write_text(source);(OUT/f'expanded-original-{n}.py').write_text(check_source);ns['run'](n)


def audit():
    FunctionType(prior.audit.__code__,dict(prior.audit.__globals__,OUT=OUT))()
    result=read(OUT/'result.json');result.update(original213failure_preserved=True,exact_conserved_coordinate_inversion=True,all_native_bit_identity_gates_retained=True);write(OUT/'result.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():
                bound=OUT/'failed-check-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(bound)==h,p
        if action=='check_interval':
            assert read(OUT/'check-receipt.json')['error'].startswith("AssertionError(('Original native identity', 64, 16,")
            write(OUT/'interval-amendment.json',dict(classification='Conjectural',
                cause='A single model across the entire archive disagrees from coarse stage16. A fresh model at that original interval start reproduces its original native rate exactly; clearing only thermodynamic jets does not. The original185producer recreated the model at every canonical continuation. Reproduce that owner lifecycle in both the precheck and actual recovery. This is distinct from the already-repaired coordinate inverse.',
                scope='All652archived stages (222coarse/430fine);658in the initial plan referred to a different later prefix and is corrected here. No extra physical step or photon solve in this check.',
                check_cap_seconds=CAPS[action],forecast='At14..16seconds/model, up to32original interval model constructions cost8..9minutes; allow20minutes. All actual photon recovery caps and every scientific gate are unchanged.',
                frozen_failure_sha256=sha(OUT/'failed-check-producer.py'),source_sha256=sha(__file__),original_gates_unchanged=True))
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];run(64 if action=='coarse' else 128)
        else:globals()['check' if action=='check_interval' else action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
