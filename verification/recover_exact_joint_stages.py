"""Counterexample candidate: recover with the exact archived Radau step size."""
from pathlib import Path
import inspect,json,os,resource,sys,time
import numpy as np
from types import FunctionType
import recover_joint_photon_stages as prior

OUT=Path('native-exact-stage205-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=120,fine=1800,coarse=1200,audit=60)


def prepare():
    assert 'collision_relative' in read(OLD/'fine-receipt.json')['error']
    assert read(OLD/'recovered-64.json')['passed']
    fn=FunctionType(prior.prepare.__code__,dict(prior.prepare.__globals__,OUT=OUT));fn()
    plan=read(OUT/'plan.json');plan.update(checkpoint='6437fd1c5',
        correction='204recomputed h by subtracting adjacent rounded edges. The original solver used a common h, whose exact bits are retained by four times the final Radau weight. Use those bits and the stored stage times directly. Preserve204fine collision-H failure1.563538e-12. This identified equation mismatch is corrected before drawing a reconstruction verdict; it is not yet proven to explain the whole failure.',
        output='Save reconstructed stage photon moments, both radial energy/number ports, actual angular packets and accepted recovery checkpoints. Preserve a rejected photon proposal and its prior accepted checkpoint on an identity failure. Never accept a moment-fitted state.',
        forecast='204coarse4steps took91.97s including setup; fine failed at6steps after140.23s. Fine8steps expected160..240s and coarse4steps90..150s under the same20/30minute caps. No native-fluid solve or full-period replay.',
        budgets=CAPS,old_failure_preserved=True)
    for n in ['fine-receipt.json','recovered-64.json','progress-128.json']:
        p=OLD/n;plan['bindings'][str(p)]=sha(p)
    write(OUT/'plan.json',plan)


def run(n):
    source=inspect.getsource(prior.run)
    changes=[
        ("if n==128:assert read(OUT/'recovered-64.json')['passed']", "if n==128:assert read(OLD/'recovered-64.json')['passed']"),
        ('timings=[]','timings=[];ports=[]'),
        ('started=time.monotonic();t,h=edges[step],edges[step+1]-edges[step];times=t+C*h',
         "started=time.monotonic();t=edges[step];h=4*z['joint_stage_weights'][2*step+1];times=z['joint_stage_times'][2*step:2*step+2]\n        assert np.array_equal(h*B,z['joint_stage_weights'][2*step:2*step+2])\n        np.savez_compressed(OUT/f'accepted-{n}.npz',x=x,step=step,t=t,moments=moments,collisions=collisions,packets=packets,ports=ports,logs=np.array(json.dumps(logs)))"),
        ('m.boundary_ports(float(times[j]),xx);packets.append(m.angular[-1])',
         'ports.append(m.boundary_ports(float(times[j]),xx)*AMP);packets.append(m.angular[-1])'),
        ('assert np.max(errors)<1e-12 and packet_error<1e-12,row',
         "if np.max(errors)>=1e-12 or packet_error>=1e-12:\n            np.savez_compressed(OUT/f'rejected-{n}.npz',x_initial=x,photon_stage_solution=pairs,gas=gas,time=t,step_size=h,stage_times=times,collision_rates=collisions[-2:],photon_moments=moments[-2:],ports=ports[-2:],angular=packets[-2:])\n            raise AssertionError(row)"),
        ('angular=np.array(packets),endpoint_occupation=actual)',
         'angular=np.array(packets),radial_ports=np.array(ports),endpoint_occupation=actual)')]
    for old,new in changes:assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
    ns=dict(prior.run.__globals__,OUT=OUT,OLD=OLD);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-recovery.py').write_text(source);ns['run'](n)


def audit():
    FunctionType(prior.audit.__code__,dict(prior.audit.__globals__,OUT=OUT))()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    started=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:run(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-started,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
