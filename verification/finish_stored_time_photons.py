"""Use archived Radau evaluation times consistently in the precise B check."""
from pathlib import Path
from types import SimpleNamespace
import gc,inspect,json,os,resource,sys,time
import numpy as np
import resume_versioned_photon_history as prior

OUT=Path('native-stored-time-photon242-work');OLD=prior.OUT
base,archive,flux=prior.base,prior.archive,prior.flux
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
seed_path=lambda n:OLD/f'accepted-{n}.npz'
CAPS=dict(prepare=300,check=600,recover=3600,assemble=300)


def prepare():
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert read(OLD/'recover-receipt.json')['error'].startswith('KeyError(')
    assert not OUT.exists();OUT.mkdir();files=[]
    for part in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64']:
        (OUT/part).mkdir(parents=True)
    for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);files += [p,dst]
    a=np.load(seed_path(64));z=np.load(prior.old.saved(64));p=np.load(OLD/'rejected-original-64.npz')
    assert int(a['step'])==int(p['step'])==115 and np.array_equal(a['x'],p['x_initial'])
    difference=p['stage_times']-(p['time']+flux.joint.C*p['step_size'])
    assert 0<max(abs(difference))<1e-18
    for p in (OLD/'clock-64').glob('original-equation-*.json'):
        if read(p)['passed']:os.link(p,OUT/'clock-64'/p.name);files.append(p)
    files += [seed_path(64),OLD/'rejected-original-64.npz',OLD/'expanded-recovery-64.py',OLD/'recover-receipt.json',OLD/'saved-pair-check.json']
    files += [Path(p) for p in read(OLD/'plan.json')['bindings']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Finish the accepted coarse photon history using its actual stored evaluation times throughout the same precise B RHS and defect.',
        cause='241restored the199primitive and passed114/115with exact native identity.116photon pair was saved, but B defect lookup reconstructed t+c*h and missed the actual recorded native-time key by a floating-point ULP. Photon coefficients and native calls already used the stored times.',
        method='Use those same archived Radau times in both high-precision B RHS evaluations and cached-rate lookup. Keep actual step weights, times, gas, branches and every numerical gate unchanged. Check clock coincidence<1e-18and exact native identity. Reuse115accepted checkpoint and the saved116pair; only117still needs a new conditional solve.',
        stored_minus_reconstructed_stage_times=difference.astype(float).tolist(),
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,CPU_affinity=2,virtual_GiB=8,
        forecast='241check plus two new pairs112.51seconds. Reuse116pair and solve117only; roughly1..3minutes plus initialization,1hour cap. No material replay or additional physical path.',
        stop='Any original exact native,whole-equation,physical,photon,endpoint,port,ledger gate or cap. Preserve240/241failures and original236/239jobs.',
        new_material_steps=0,reused_accepted_photon_steps=115,reused_failed_proposal=True,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(OLD/'symbolic.json'))


interval_model=bind(prior.interval_model,OUT=OUT)


def original_equation(m,z,step,t,h,x,gas,photons,times,cs,ss):
    assert max(abs(times-(t+flux.joint.C*h)))<1e-18
    functions={}
    for name,changes in {
        'precise_rhs':[('for c,g,(J,_) in zip(joint.C,guides,maps)','for now,g,(J,_) in zip(stage_times,guides,maps)'),('m,t+c*h,g','m,now,g')],
        'precise_defect':[('for c,g in zip(joint.C,pairs)','for now,g in zip(stage_times,pairs)'),('m.precise_values[t+c*h]','m.precise_values[now]')]
    }.items():
        fn=getattr(flux,name);s=inspect.getsource(fn)
        for a,b in changes:assert s.count(a)==1,(name,a);s=s.replace(a,b)
        ns=dict(fn.__globals__,stage_times=times);exec(compile(s,__file__,'exec'),ns);functions[name]=ns[name]
    fixed=SimpleNamespace(**dict(vars(flux),**functions))
    return bind(prior.old.original_equation,OUT=OUT,flux=fixed)(m,z,step,t,h,x,gas,photons,times,cs,ss)


def check():
    bind(base.base.prior.initialize,OUT=OUT)()
    z=np.load(prior.old.saved(64));p=np.load(OLD/'rejected-original-64.npz');m=interval_model(z,64,115)
    times=p['stage_times'];gas=p['gas'];photons=p['photon_stage_solution']
    for g,q in zip(gas,z['joint_stage_conserved_scaled'][230:232]):assert np.array_equal(m.conserved(g),q)
    row=original_equation(m,z,115,p['time'][()],p['step_size'][()],p['x_initial'],gas,photons,times,[m.local(v) for v in times],[m.source(v) for v in times])
    write(OUT/'saved-pair-check.json',dict(classification='Counterexample candidate',passed=row['passed'],actual_equation=row,new_photon_solves=0,new_material_steps=0))


def recover():
    assert read(OUT/'saved-pair-check.json')['passed'];source=(OLD/'expanded-recovery-64.py').read_text()
    ns=dict(base.run.__globals__,OUT=OUT,OLD=OLD,ANCHOR117=prior.old.ANCHOR117,saved=prior.old.saved,
        seed_folder=lambda n:prior.old.recovery.OUT,seed=bind(prior.prior.seed,seed_path=seed_path),seed_path=seed_path,RESUME=True,
        checkpoint=bind(archive.checkpoint,OUT=OUT),snapshot=bind(archive.snapshot,OUT=OUT),
        restored_gas=prior.prior.restored_gas,original_equation=original_equation,
        interval_model=interval_model,interval_starts=prior.prior.interval_starts,gc=gc)
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-recovery-64.py').write_text(source);ns['run'](64)


def assemble():
    bind(prior.old.assemble,OUT=OUT)()
    r=read(OUT/'result.json');r.update(original240_and241_failures_preserved=True,original199_primitive_reproduced=True,
        actual_stored_Radau_times_used=True,reused_accepted_steps=115,reused_failed_photon_proposal=True)
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
