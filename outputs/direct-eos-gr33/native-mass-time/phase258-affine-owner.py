"""Separate native energy forcing on each accepted physical branch."""
from pathlib import Path
import gc,json,os,resource,time
import numpy as np
import finish_returned_material_accuracy as producer
import resolve_returned_native_time as previous

OUT=Path('native-mass-time258-work');OLD=producer.OUT
read,write,sha,bind=producer.read,producer.write,producer.sha,producer.bind
LD=np.longdouble;AMP=producer.base.actual.AMP
assert not (OUT/'canonical-affine-start.json').exists()
assert sha(OUT/'failed-affine-producer.py')==read(OUT/'affine-receipt.json')['source_sha256']
assert 'Saved-stage anchor precision' in read(OUT/'affine-receipt.json')['error']
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2)
producer.joint.previous.original.inf.incident.native.deadline(3600)
start=time.monotonic();error=None
write(OUT/'canonical-affine-start.json',dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),source_sha256=sha(__file__)))
write(OUT/'canonical-affine-plan.json',dict(classification='Conjectural',
    claim='Separate affine forcing from state-dependent native thermal work on each actually accepted branch, without applying a foreign coarse direction to the fine high/low branch selector.',
    reason='Both prior diagnostics retained a model across canonical intervals, contrary to the original producer and the already established214/226cache lifetime. The affine replay reproduced the exact226error3.8747131e-12. Restore a fresh model at each original canonical interval and require all stored rates below the unchanged1e-12gate. The earlier foreign-direction explanation was an unconfirmed hypothesis. Evaluate zero gas on each verified actual branch only; not a replacement physical law.',
    decision='If forcing dominates, trace its actual geometry/energy owner; if state evolution dominates, retain the failure and assess the physical time solver. No diagnostic energy is added to final charge.',
    cap_seconds=3600,virtual_GiB=8,CPU_affinity=2,physical_steps=0,new_clock_paths=0,
    forecast='Measured full native call0.484..0.518seconds.700accepted stages plus branch replays and32original canonical constructors assumed10..25minutes; retain one-hour headroom. No extra physical solve. Stop on any1e-12saved-native mismatch, branch replay mismatch or resource cap.',
    original_mass_time_gate=.02,final_charge_conclusion='unadjudicated',full_goal_complete=False,
    bindings={str(p):sha(p) for p in [Path(__file__),Path(previous.__file__),Path('verification/continue_dense_gr_return.py'),OUT/'affine-receipt.json',OUT/'failed-affine-producer.py',OLD/'result.json']}))
try:
    Model=bind(producer.base.initialize,OUT=OUT)();rows=[];integrals={}
    for n in [64,128]:
        m=Model(n)
        with np.load(OLD/f'sweep-1/photons/return-{n}.npz') as z:
            times=z['joint_stage_times'];weights=z['joint_stage_weights'];states=z['joint_stage_conserved_scaled'];stored=z['joint_native_rates_scaled']
        forcing=[];errors=[];canonical=0
        for i,t in enumerate(times):
            t=float(t);interval=int(np.clip(np.searchsorted(m.t,t-1e-18,side='right')-1,0,15))
            if interval!=canonical:
                assert interval==canonical+1,(n,i,interval,canonical)
                del m;gc.collect();m=Model(n);canonical=interval
            g=previous.gas(m,states[i]);actual=m.native(t,g)*m.units
            den=np.maximum(np.sum(abs(stored[i]),axis=0),LD('1e-290'))
            err=(np.sum(abs(actual-stored[i]),axis=0)/den).astype(float).tolist()
            assert max(err)<1e-12,(n,i,t,err)
            tangent,reset=m.selected[t];reset(False)
            replay=m.native(t,g,tangent=tangent)*m.units
            assert np.array_equal(replay,actual),(n,i,'branch replay')
            reset(False);zero=m.native(t,np.zeros_like(g),tangent=tangent)*m.units
            forcing.append(zero);errors.append(err)
            if i%64==0:
                np.savez_compressed(OUT/f'canonical-affine-prefix-{n}.npz',forcing=forcing,times=times[:i+1],weights=weights[:i+1])
                write(OUT/'canonical-affine-progress.json',dict(clock=n,completed=i+1,total=len(times),canonical=canonical,seconds=time.monotonic()-start))
        forcing=np.array(forcing);state=stored-forcing
        values=[np.sum(weights[:,None,None].astype(LD)*v,axis=0,dtype=LD)*AMP for v in [stored,forcing,state]]
        integrals[n]=values
        np.savez_compressed(OUT/f'canonical-affine-{n}.npz',times=times,weights=weights,forcing=forcing,state_dependent=state,integrated=values)
        rows.append(dict(clock=n,max_saved_native_relative=float(np.max(errors)),thermal_erg=[float(np.sum(v[:,0])) for v in values]))
        del m
    differences=[a-b for a,b in zip(integrals[64],integrals[128])]
    norm=np.sum(abs(differences[0][:,0]));identity=float(np.sum(abs(differences[0][:,0]-differences[1][:,0]-differences[2][:,0]))/norm)
    assert identity<1e-8,identity
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,
        labels=['native','branch_affine_forcing','state_dependent'],thermal_difference_erg=[float(np.sum(v[:,0])) for v in differences],
        thermal_difference_L1=[float(np.sum(abs(v[:,0]))) for v in differences],identity_relative=identity,
        scope='Branchwise exact saved-history split. The affine part can depend on the accepted branch; no global linearity, unique cause or repaired physical solution is claimed.',
        physical_steps=0,new_clock_paths=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'canonical-affine-result.json',result);print(json.dumps(result),flush=True)
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'canonical-affine-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
