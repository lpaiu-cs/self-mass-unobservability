"""Reuse the exact-vacuum optical ODE inside actual photon propagation."""
from pathlib import Path
import json, resource, sys, time
import numpy as np
from scipy.integrate import solve_ivp
import propagate_physical_exterior as base

OLD=base.OUT
base.OUT=OUT=Path('native-dynamic-vacuum261-work')
SlowMetric=base.Metric


class Metric(SlowMetric):
    def __init__(self,order):
        super().__init__(order);d=self.d;span=2*base.C*d.T
        def speed(y):return float(d.bg.metric(np.array([(d.r0+span*y)/d.model.m.R]))[3][0])
        radial=solve_ivp(lambda s,y:[speed(y[0])],[0,1],[0.],method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True)
        optical=solve_ivp(lambda s,y:[1/speed(s)],[0,1],[0.],method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True)
        assert radial.success and optical.success
        def inverse(x):
            x=np.asarray(x);assert np.min(x,initial=0)>=-1e-7 and np.max(x,initial=0)<=span*(1+1e-12)
            return d.r0+span*radial.sol(np.clip(x.ravel()/span,0,1))[0].reshape(x.shape)
        def at(r):
            r=np.asarray(r);z=(r-d.r0)/span
            assert np.min(z,initial=0)>-1e-12 and np.max(z,initial=0)<1+1e-12
            return span*optical.sol(np.clip(z.ravel(),0,1))[0].reshape(r.shape)
        self.original_optical=d.optical;d.inverse=inverse;d.optical=at


base.Metric=Metric
Photons=base.Photons;initialize=base.initialize
read,write,sha=base.read,base.write,base.sha


def prepare():
    assert read(OLD/'pipeline-status.json')['state']=='failed'
    assert read(OLD/'cost-stop.json')['accepted_physical_steps']==0
    base.prepare();p=read(OUT/'plan.json')
    p['cost_repair']='Replace repeated full-envelope inverse calls by dr/dx=N*sqrt(b) and dx/dr=1/(N*sqrt(b)) on the identical exterior vacuum. This is a numerical cache of the same optical map, not a new physical grid, trajectory clock, metric or tolerance.'
    p['original_attempt']='Original source, plan, check and cost-stop remain under native-dynamic-photon261-work; intentional pilot cost termination, not a physical acceptance failure.'
    for f in [Path(__file__),OLD/'cost-stop.json',OLD/'plan.json',OLD/'check.json']:
        p['bindings'][str(f)]=sha(f)
    write(OUT/'plan.json',p)


def check():
    base.check();m=Metric(8);d=m.d;span=2*base.C*d.T
    x=np.unique(np.r_[d.ex,m.nodes.ravel(),np.linspace(0,span,129)])
    r=d.inverse(x)
    independent=float(np.max(abs(m.original_optical(r)-x))/span)
    inverse=float(np.max(abs(d.optical(r)-x))/span)
    assert max(independent,inverse)<1e-10,(independent,inverse)
    slow=SlowMetric(8);now=np.linspace(.05*d.D,.9*d.D,32)
    rr=d.inverse(np.linspace(0,.3*base.C*d.D,32))
    t=time.monotonic();old=slow.at(now,rr);old_seconds=time.monotonic()-t
    t=time.monotonic();new=m.at(now,rr);new_seconds=time.monotonic()-t
    errors={k:float(np.max(abs(new[k]-v))/max(np.max(abs(v)),1e-290)) for k,v in old.items()}
    assert max(errors.values())<1e-9,errors
    write(OUT/'optical-check.json',dict(classification='Counterexample candidate',passed=True,
        independent_original_optical_fraction=independent,inverse_fraction=inverse,
        identical_metric_relative=errors,query_packets=len(now),old_seconds=old_seconds,new_seconds=new_seconds,
        actual_wave_values_and_derivatives_preserved_to_tested_tolerance=True,
        scope='Finite numerical same-map comparison and analytic vacuum optical ODE; not a continuum certificate.'))


if __name__=='__main__':
    action=sys.argv[1];start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    if action=='prepare':cap=120
    else:
        p=read(OUT/'plan.json');cap=p['budgets'].get(action,7200)
        for path,h in p['bindings'].items():assert sha(path)==h,path
    base.incident.native.deadline(cap)
    try:
        if action=='prepare':prepare()
        elif action=='check':check()
        elif action in base.SETTINGS:base.production(action)
        else:getattr(base,action)()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
