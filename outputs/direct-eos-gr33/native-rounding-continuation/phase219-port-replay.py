"""Replay only the unsaved conditional photon block and inspect its original port."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,time
import numpy as np
import resume_coordinate_exact_photons as prior

OUT=Path('native-port-order219-work');OLD=prior.OUT
LD,AMP=prior.LD,prior.base.AMP
read,write,sha=prior.read,prior.write,prior.sha
assert not OUT.exists();OUT.mkdir();start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.base.joint.previous.original.inf.incident.native.deadline(900)
files=[]
for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-128']:(OUT/folder).mkdir(parents=True)
for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
    os.link(src,OUT/src.relative_to(OLD));files.append(src)
files += [OLD/n for n in ['accepted-128.npz','restart-check.json','clock-128/snapshot-02.json','resume-controller-status.json']]
files += [Path(v.__file__) for v in list(__import__('sys').modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
write(OUT/'plan.json',dict(classification='Conjectural',
    claim='Determine whether the actual original radial-port failure comes from a mismatched accumulation or boundary timestamp before continuing same-history recovery.',
    method='Reuse214accepted15fine conditional steps and652native identity controls. Recompute only its unsaved16th photon block; preserve the original equation/native/conditional tests and snapshot failure. Compare original unscaled binary64 Radau accumulation and saved angular timestamps. No fluid steps or fitted port.',
    budgets=dict(seconds=900,CPU_threads=1,virtual_GiB=6),forecast='214conditional stages about20..45s plus14s model; allow15minutes for one block and source-time checks.',
    stop='After this original snapshot whether it passes or fails; no automatic history or grid expansion.',
    bindings={str(p):sha(p) for p in files+[Path(__file__),prior.prior.saved(128)]},final_charge_conclusion='unadjudicated'))


class ProbeDone(Exception):pass


def snapshot(n,m,z,step,photons,collisions,ports,packets):
    original=FunctionType(prior.prior.snapshot.__code__,dict(prior.prior.snapshot.__globals__,OUT=OUT))
    try:row=original(n,m,z,step,photons,collisions,ports,packets)
    except AssertionError:row=read(OUT/f'clock-{n}/snapshot-02.json')
    assert step==15 and len(ports)==32
    values=np.array(ports);weights=z['joint_stage_weights'][:32];target=z['radial_ports'][2]
    def replay(p):
        value=np.zeros((2,2),float)
        for k in range(16):value+=float(4*weights[2*k+1])*(.75*np.asarray(p[2*k]/AMP,float)+.25*np.asarray(p[2*k+1]/AMP,float))
        return value*AMP
    def relative(p):return np.max(abs(p-target)/np.maximum(abs(target),LD('1e-290')))
    old_times=np.array(z['accepted_angular_times'][:32]);new_times=np.array(z['joint_stage_times'][:32],float)
    delta=[];zero=np.zeros_like(photons)
    for a,b in zip(old_times,new_times):delta.append((m.boundary_ports(a,zero)-m.boundary_ports(b,zero))*AMP)
    delta=np.array(delta)
    direct=m.boundary_ports(old_times[-1],photons)-m.boundary_ports(new_times[-1],photons)
    zero_direct=m.boundary_ports(old_times[-1],zero)-m.boundary_ports(new_times[-1],zero)
    result=dict(classification='Counterexample candidate',original_snapshot=row,
        original_order_relative=float(relative(replay(values))),timestamp_change_count=int(np.count_nonzero(old_times!=new_times)),
        maximum_timestamp_change=float(np.max(abs(old_times-new_times))),time_shifted_original_order_relative=float(relative(replay(values+delta))),
        last_photon_zero_time_difference_equal=bool(np.array_equal(direct,zero_direct)),
        original_port=target.astype(float).tolist(),recovered_port=replay(values).astype(float).tolist(),
        new_physical_steps=0,new_conditional_steps=1,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'snapshot-input.npz',ports=values,weights=weights,target=target,old_times=old_times,new_times=new_times,delta=delta,photons=photons,collisions=collisions,packets=packets)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);raise ProbeDone()


try:
    source=inspect.getsource(prior.run)
    mark="snapshot=FunctionType(prior.snapshot.__code__,dict(prior.snapshot.__globals__,OUT=OUT))"
    assert source.count(mark)==1;source=source.replace(mark,'snapshot=probe_snapshot')
    seed_path=lambda n:OLD/f'accepted-{n}.npz'
    seed=FunctionType(prior.seed.__code__,dict(prior.seed.__globals__,seed_path=seed_path))
    ns=dict(prior.run.__globals__,OUT=OUT,RESUME=True,seed_path=seed_path,seed=seed,probe_snapshot=snapshot)
    exec(compile(source,__file__,'exec'),ns);ns['run'](128)
except ProbeDone:pass
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
