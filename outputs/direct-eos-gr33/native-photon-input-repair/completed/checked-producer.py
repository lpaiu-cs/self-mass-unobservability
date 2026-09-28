"""Counterexample candidate: reproduce the original checkpoint input, without port fitting."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import resume_coordinate_exact_photons as prior

OUT=Path('native-archival-photon220-work')
INPUT=Path('native-front-continuation184-work/sweep-1/photons/input-128.npz')
SEED=Path('native-equation-recovery208-work')
read,write,sha=prior.read,prior.write,prior.sha
LD,AMP=prior.LD,prior.base.AMP
CAPS=dict(prepare=180,fine=900)


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[]
    assert read(prior.OUT/'restart-check.json')['passed']
    assert read(SEED/'result.json')['original_endpoint_and_ledger_passed']
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-128']:(OUT/folder).mkdir(parents=True)
    for p in list((prior.OUT/'sweep-0').rglob('*.npz'))+[prior.OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        os.link(p,OUT/p.relative_to(prior.OUT));files.append(p)
    files += [INPUT,prior.prior.saved(128),SEED/'recovered-128.npz',SEED/'recovered-128.json',SEED/'result.json',prior.OUT/'restart-check.json',Path('native-port-order219-work/result.json')]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Test the actual failed conditional-recovery interval with the exact photon initial state used by its original producer, before changing solver accuracy or gates.',
        reason='219reproduces the failed outer photon-number port4.74071e-12; original accumulation still fails and all boundary timestamps match. The current recovery carries its own approximate endpoint into the next original model interval, whereas the original producer reads its stored restart_x.',
        method='Reuse208first8fine stage moments/ports/checks, but start the next8conditional photon blocks from184exact restart_x at the same stored checkpoint. Verify all archived material/time prefixes and checkpoint physical state. Do not force any final or stage port. Keep exact conserved-coordinate inversion and original per-interval model lifetime.',
        controls='Report the old approximate versus original initial photon difference and inherited prefix ports separately. Test every original coupled equation/native identity and the unchanged final snapshot/ledger/port gates. A changed starting photon approximation is not native-fluid reintegration.',
        gates=read(prior.prior.OUT/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='219one block including model construction45.95s. Eight conditional blocks at20..45s implyabout3..6minutes; allow15minutes, no physical steps, one existing fine path only.',
        stop='At the original T/16snapshot whether pass or fail, or any original gate/cap failure. Preserve each last pair before snapshot testing and accepted checkpoint after each pass. No automatic longer recovery or tighter solve.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    np.savez_compressed(OUT/'seed-128.npz',step=8)


def seed(m,z,n):
    assert n==128;p=dict(np.load(SEED/'recovered-128.npz'));old=dict(np.load(INPUT));r=read(SEED/'recovered-128.json');begin=r['recovered_steps'];assert begin==8
    assert old['t'][-1]==z['t'][1]==z['actual_step_edges'][begin]
    for key in ['joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled']:
        assert np.array_equal(old[key],z[key][:2*begin]),key
    assert np.array_equal(old['photon_history_scaled_occupation'],z['photon_history_scaled_occupation'][:2])
    x=old['restart_x'].copy();previous=p['endpoint_occupation']/(m.scale*AMP)
    relative=float(np.sum(abs(x-previous))/max(np.sum(abs(x)),LD('1e-290')))
    photon_identity=float(np.sum(abs(x*m.scale*AMP-z['photon_history_scaled_occupation'][1]))/max(np.sum(abs(z['photon_history_scaled_occupation'][1])),LD('1e-290')))
    ports=np.sum(p['weights'][:,None,None]*p['radial_ports'],axis=0,dtype=LD);target=z['radial_ports'][1]
    prefix_error=float(np.max(abs(ports-target)/np.maximum(abs(target),LD('1e-290'))))
    assert photon_identity==0 and prefix_error<1e-12
    write(OUT/'seed-identity.json',dict(classification='Counterexample candidate',exact_original_checkpoint_photon_input=True,
        all_material_and_time_prefixes_exact=True,original_photon_snapshot_relative=photon_identity,old_recovered_initial_relative=relative,
        inherited_prefix_port_relative=prefix_error,reused_steps=begin,new_physical_steps=0))
    return begin,x,r['rows'],list(p['photon_moments']),list(p['collision_rates']),list(p['angular']),list(p['radial_ports']),relative


def fine():
    source=inspect.getsource(prior.run)
    mark="    exec(compile(source,__file__,'exec'),ns);"
    changes=[("target=float(z['t'][-1])","target=float(z['t'][2])"),
        ("expected=z['photon_history_scaled_occupation'][-1]","expected=z['photon_history_scaled_occupation'][2]"),
        ("        row['snapshot']=snapshot(","        np.savez_compressed(OUT/'current-pair.npz',step=step,time=t,step_size=h,x_initial=x,solution=pairs,gas=gas,rhs=rhs,residual=residual)\n        row['snapshot']=snapshot(")]
    inject=''
    for a,b in changes:inject+=f'    assert source.count({a!r})==1\n    source=source.replace({a!r},{b!r})\n'
    assert source.count(mark)==1;source=source.replace(mark,inject+mark)
    ns=dict(prior.run.__globals__,OUT=OUT,RESUME=True,seed_path=lambda n:OUT/'seed-128.npz',seed=seed)
    exec(compile(source,__file__,'exec'),ns);ns['run'](128)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.base.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
