"""Evaluate isolated small source moments without extended-precision powers.

The moment construction remains longdouble. No high/low states are combined.
"""
from pathlib import Path
import json, resource, sys, time
import numpy as np
import propagate_returned_exterior as base

OLD=base.OUT
base.OUT=base.p.OUT=OUT=Path('native-returned-moments262-work')
SlowMoments=base.Moments


class Moments(SlowMoments):
    def __init__(self,*args):
        super().__init__(*args)
        self.a=np.asarray(self.a,float); self.b=np.asarray(self.b,float)
        self.q=np.asarray(self.q,float)

    def primitive(self,s,bins,k=0):
        v=np.clip(np.asarray(s)/self.T,0,self.edges[-1])
        j=np.clip(np.searchsorted(self.edges,v,side='right')-1,0,len(self.edges)-2)
        lo=self.edges[j]; a=self.a[j,bins]; b=self.b[j,bins]
        return self.q[k,j,bins]+self.T*((a-b*lo)*(v**(k+1)-lo**(k+1))/(k+1)+b*(v**(k+2)-lo**(k+2))/(k+2))


def prepare():
    base.prepare(); plan=base.read(OUT/'plan.json')
    plan['arithmetic']='Preserve longdouble construction of cumulative moments, then evaluate the isolated returned metric in float64. No large background subtraction or high/low state addition. Require same-query whole metric and derivative comparison <1e-9 against the unchanged longdouble producer.'
    plan['bindings'][str(Path(__file__))]=base.sha(__file__)
    base.write(OUT/'plan.json',plan)


def check():
    base.p.initialize(); base.Moments=SlowMoments; old=base.Metric(8)
    base.Moments=Moments; new=base.Metric(8); d=new.d
    t=np.linspace(.01,.99,256)*d.T; r=d.inverse(.2*base.C*t)
    start=time.monotonic(); a=old.at(t,r); slow=time.monotonic()-start
    start=time.monotonic(); b=new.at(t,r); fast=time.monotonic()-start
    errors={k:float(np.max(abs(b[k]-v))/max(np.max(abs(v)),1e-290)) for k,v in a.items()}
    assert max(errors.values())<1e-9,errors
    base.write(OUT/'arithmetic-check.json',dict(classification='Counterexample candidate',passed=True,
        same_queries=256,old_seconds=slow,new_seconds=fast,relative=errors,scope='Finite metric and derivative comparisons; not a universal rounding-error certificate.'))
    base.check()


if __name__=='__main__':
    action=sys.argv[1]; start=time.monotonic(); error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2); base.p.incident.native.deadline(1800 if action=='pilot' else 600)
    try:
        if action!='prepare':
            for f,h in base.read(OUT/'plan.json')['bindings'].items():assert base.sha(f)==h,f
        base.Moments=Moments
        if action in ['prepare','check']:globals()[action]()
        else:getattr(base,action)()
    except BaseException as exc:error=repr(exc); raise
    finally:
        if OUT.exists():base.write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=base.sha(__file__),adapter_sha256=base.sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
