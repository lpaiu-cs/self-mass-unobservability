"""Resume240 with the actual199 primitive owner from accepted stage114."""
from pathlib import Path
import gc,json,os,resource,sys,time
import numpy as np
import complete_coarse_photon_history as old
import repair_native_affine_precision as primitive

OUT=Path('native-versioned-photon241-work');OLD=old.OUT
read,write,sha,bind=old.read,old.write,old.sha,old.bind
prior,archive,base,flux=old.prior,old.archive,old.base,old.flux
seed_path=lambda n:OLD/f'accepted-{n}.npz'
CAPS=dict(prepare=300,check=600,recover=3600,assemble=300)


def prepare():
    assert read(OLD/'controller-status.json')['state']=='failed'
    failure=read(OLD/'clock-64/original-equation-113.json')
    assert not failure['passed'] and failure['native_identity_absolute'][0]>0
    assert not OUT.exists();OUT.mkdir();files=[]
    for part in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64']:
        (OUT/part).mkdir(parents=True)
    for p in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);files += [p,dst]
    a=np.load(seed_path(64));z=np.load(old.saved(64));proposal=np.load(OLD/'rejected-original-64.npz')
    assert int(a['step'])==113 and int(proposal['step'])==113
    assert np.array_equal(a['x'],proposal['x_initial']) and a['t']==z['actual_step_edges'][113]
    for p in (OLD/'clock-64').glob('original-equation-*.json'):
        if read(p)['passed']:os.link(p,OUT/'clock-64'/p.name);files.append(p)
    files += [seed_path(64),OLD/'rejected-original-64.npz',OLD/'expanded-recovery-64.py',OLD/'plan.json',OLD/'recover-receipt.json',
        primitive.OUT/'primitive-result.json',Path('native-precise-continuation199-work/plan.json')]
    files += [Path(p) for p in read(OLD/'plan.json')['bindings']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the actual archived primitive-precision transition to the saved114photon pair, then finish the same original coarse history without replaying accepted photon/material steps.',
        cause='240recovered112/113 exactly but114used the pre199primitive.199accepted114only after promoting five primitive outputs and casts to long double.202B precision starts at115and remains required. Reuse the actual199installer; do not alter native archival identity or full-equation gates.',
        method='Reuse240accepted113checkpoint and its failed114photon proposal. A fresh model with the199native primitive must exactly reproduce both archived114native rates and pass the original whole equation before recovery. Re-evaluate that saved proposal against the current photon RHS; compute only115..117missing pairs, then reuse actual118/119captures.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,CPU_affinity=2,virtual_GiB=8,
        forecast='240three conditional pairs and checks95.21seconds including setup. Reused114pair and three remaining photon solves should take2..4minutes plus setup; allow1hour. Check has10minute cap. No native evolution, extra grid,period or path.',
        stop='Any original exact native/whole-equation/physical/photon/endpoint/port/ledger gate or cap; preserve240failure and rejected pair. Original236/239jobs remain untouched.',
        new_material_steps=0,reused_accepted_photon_steps=113,reused_failed_proposal=True,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(OLD/'symbolic.json'))


def interval_model(z,n,step):
    assert n==64 and step>=113
    bind(primitive.install_primitive_precision,OUT=OUT)()
    m=base.base.owner.Model(n)
    with flux.precision.mp.workdps(60):m.precise_tangent=flux.precision.build(flux.joint,flux.owner)
    m.precise_values={}
    return m


original_equation=bind(old.original_equation,OUT=OUT)


def check():
    bind(base.base.prior.initialize,OUT=OUT)()
    z=np.load(old.saved(64));p=np.load(OLD/'rejected-original-64.npz');m=interval_model(z,64,113)
    t,h=p['time'][()],p['step_size'][()];times=p['stage_times'];gas=p['gas'];photons=p['photon_stage_solution']
    for g,q in zip(gas,z['joint_stage_conserved_scaled'][226:228]):assert np.array_equal(m.conserved(g),q)
    row=original_equation(m,z,113,t,h,p['x_initial'],gas,photons,times,[m.local(v) for v in times],[m.source(v) for v in times])
    write(OUT/'saved-pair-check.json',dict(classification='Counterexample candidate',passed=row['passed'],
        original199_primitive_applied=True,actual_equation=row,new_photon_solves=0,new_material_steps=0))


def recover():
    assert read(OUT/'saved-pair-check.json')['passed']
    source=(OLD/'expanded-recovery-64.py').read_text()
    a='sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]'
    b="sol=np.tile(x[None],(2,1,1,1)).ravel();calls=[]\n        if step==begin:\n            proposal=np.load(OLD/'rejected-original-64.npz');assert np.array_equal(x,proposal['x_initial'])\n            sol=proposal['photon_stage_solution'].copy().ravel()"
    assert source.count(a)==1;source=source.replace(a,b)
    ns=dict(base.run.__globals__,OUT=OUT,OLD=OLD,ANCHOR117=old.ANCHOR117,saved=old.saved,seed_folder=lambda n:old.recovery.OUT,
        seed=bind(prior.seed,seed_path=seed_path),seed_path=seed_path,RESUME=True,
        checkpoint=bind(archive.checkpoint,OUT=OUT),snapshot=bind(archive.snapshot,OUT=OUT),
        restored_gas=prior.restored_gas,original_equation=original_equation,
        interval_model=interval_model,interval_starts=prior.interval_starts,gc=gc)
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-recovery-64.py').write_text(source);ns['run'](64)


def assemble():
    bind(old.assemble,OUT=OUT)()
    r=read(OUT/'result.json');r.update(original240_failure_preserved=True,original199_primitive_reproduced=True,
        reused_240_accepted_steps=2,reused_failed_photon_proposal=True)
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
