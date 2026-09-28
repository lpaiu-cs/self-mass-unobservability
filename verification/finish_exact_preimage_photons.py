"""Select an equivalent conserved preimage using the original B equation."""
from pathlib import Path
import gc,json,os,resource,sys,time
import numpy as np
import finish_stored_time_photons as prior

OUT=Path('native-preimage-photon243-work');OLD=prior.OUT;SEED=prior.OLD
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
base,archive,flux=prior.base,prior.archive,prior.flux
seed_path=lambda n:SEED/f'accepted-{n}.npz'
CAPS=dict(prepare=300,check=600,recover=3600,assemble=300)


def prepare():
    failure=read(OLD/'clock-64/original-equation-115.json')
    assert not failure['passed'] and max(failure['native_identity_absolute'])==0
    assert not OUT.exists();OUT.mkdir();files=[]
    for part in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64']:(OUT/part).mkdir(parents=True)
    for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);files += [p,dst]
    assert int(np.load(seed_path(64))['step'])==115
    for p in (SEED/'clock-64').glob('original-equation-*.json'):
        if read(p)['passed']:os.link(p,OUT/'clock-64'/p.name);files.append(p)
    files += [seed_path(64),SEED/'rejected-original-64.npz',SEED/'expanded-recovery-64.py',OLD/'check-receipt.json',
        OLD/'actual-defect.npz',OLD/'actual-defect-location.json',OLD/'actual-time-types.json']
    files += [Path(p) for p in read(OLD/'plan.json')['bindings']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Restore the normalized B coordinate of the SAME archived material states, then finish the missing coarse photons at all original gates.',
        evidence='242exactly reproduces archived native rates but116whole-vector3.79152e-11fails. Its dominant second-stage B defect is0.0312499996at cell261; reconstructed gas B is5.04107e17. Conserved-coordinate inversion need not recover the original normalized floating-point preimage. The post-floor endpoint erases that inactive cell.',
        method='Use independently evaluated original20260-digit B rates and unchanged Radau weights to select at most8adjacent B floats nearer the original stage equation. Accept a coordinate move ONLY if EVERY conserved component remains bit-identical to the archived state. No state/rate/collision fitting or material integration. Exact native identity, unchanged photon RHS and full coupled/physical/endpoint/port/ledger gates remain mandatory. Failure preserves all data.',
        scope='Only original116/117 pairs under the original rounded conserved-input native owner; do not apply this decoder to238exact-conserved-input arithmetic. Reuse115accepted photons and retained116pair.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,CPU_affinity=2,virtual_GiB=8,
        forecast='242saved-pair check31.67s. Decoder adds at most8small coordinate maps and two native B calls; remaining117conditional solve measured order30s. Allow10minutes check and1hour completion.',
        stop='No equivalent preimage, changed conserved value, native identity, whole/physical/photon/endpoint/port/ledger gate or cap. Preserve240..242failures and236/239live jobs.',
        new_material_steps=0,bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(OLD/'symbolic.json'))


interval_model=bind(prior.interval_model,OUT=OUT)
original_equation=bind(prior.original_equation,OUT=OUT)


def restore_pair(m,z,step,gas,t,h,times):
    assert step in [115,116];gas=np.array(gas,copy=True);before=gas.copy()
    initial=prior.prior.prior.restored_gas(m,z['joint_stage_conserved_scaled'][2*step-1]);initial[~m.material.active(t)]=0
    expected=z['joint_stage_conserved_scaled'][2*step:2*step+2]
    for g,q in zip(gas,expected):assert np.array_equal(m.conserved(g),q)
    with flux.precision.mp.workdps(60):
        rates=np.array([flux.precision.native_B(m,now,g,m.precise_tangent) for now,g in zip(times,gas)])
        target=flux.precision.cast(flux.precision.hp(initial[:,2])+flux.precision.hp(h)*(flux.precision.hp(flux.joint.A)@rates))
    counts=np.zeros(gas.shape[:2],int)
    for _ in range(8):
        changes=0
        for j,g in enumerate(gas):
            trial=g.copy();trial[:,2]=np.nextafter(g[:,2],np.where(target[j]>g[:,2],np.longdouble('inf'),np.longdouble('-inf')))
            same=np.all(m.conserved(trial)==expected[j],axis=0)
            better=abs(trial[:,2]-target[j])<abs(g[:,2]-target[j]);mask=same&better
            g[mask,2]=trial[mask,2];counts[j,mask]+=1;changes+=np.count_nonzero(mask)
        if not changes:break
    for g,q in zip(gas,expected):assert np.array_equal(m.conserved(g),q),'Decoder changed an archived physical state'
    record=dict(classification='Counterexample candidate',step=step,all_conserved_components_exact=True,
        changed_coordinates=int(np.count_nonzero(counts)),maximum_adjacent_moves=int(np.max(counts)),
        changed_indices=np.argwhere(counts).tolist())
    write(OUT/f'preimage-{step}.json',record)
    np.savez_compressed(OUT/f'preimage-{step}.npz',before=before,after=gas,target_B=target,conserved=expected)
    return gas


def check():
    bind(base.base.prior.initialize,OUT=OUT)()
    z=np.load(prior.prior.old.saved(64));p=np.load(SEED/'rejected-original-64.npz');m=interval_model(z,64,115)
    t,h=p['time'][()],p['step_size'][()];times=p['stage_times']
    gas=restore_pair(m,z,115,p['gas'],t,h,times)
    row=original_equation(m,z,115,t,h,p['x_initial'],gas,p['photon_stage_solution'],times,[m.local(v) for v in times],[m.source(v) for v in times])
    write(OUT/'saved-pair-check.json',dict(classification='Counterexample candidate',passed=row['passed'],actual_equation=row,new_photon_solves=0,new_material_steps=0))


def recover():
    assert read(OUT/'saved-pair-check.json')['passed'];source=(SEED/'expanded-recovery-64.py').read_text()
    marker='        cs=[m.local(now) for now in times];ss=[m.source(now) for now in times]'
    assert source.count(marker)==1;source=source.replace(marker,'        gas=restore_pair(m,z,step,gas,t,h,times)\n'+marker)
    ns=dict(base.run.__globals__,OUT=OUT,OLD=SEED,ANCHOR117=prior.prior.old.ANCHOR117,saved=prior.prior.old.saved,
        seed_folder=lambda n:prior.prior.old.recovery.OUT,seed=bind(prior.prior.prior.seed,seed_path=seed_path),seed_path=seed_path,RESUME=True,
        checkpoint=bind(archive.checkpoint,OUT=OUT),snapshot=bind(archive.snapshot,OUT=OUT),
        restored_gas=prior.prior.prior.restored_gas,original_equation=original_equation,restore_pair=restore_pair,
        interval_model=interval_model,interval_starts=prior.prior.prior.interval_starts,gc=gc)
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-recovery-64.py').write_text(source);ns['run'](64)


def assemble():
    bind(prior.prior.old.assemble,OUT=OUT)()
    r=read(OUT/'result.json');r.update(original240_242_failures_preserved=True,actual_stored_Radau_times_used=True,
        exact_conserved_preimage_selected=True,reused_accepted_steps=115,reused_failed_photon_proposal=True)
    write(OUT/'result.json',r);write(OUT/'recovered-64.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    os.sched_setaffinity(0,{2});resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    base.joint.previous.original.inf.incident.native.deadline(CAPS[action]);start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
