"""Exact-input native memoization; retain the rejected original forecast."""
from functools import lru_cache
from types import FunctionType
from pathlib import Path
import inspect,json,signal,sys,time
import numpy as np
import def_native_conserved_history_charge as run

OUT=run.OUT;write=run.write;sha=run.sha
Original=run.chem.old.Native;original_setup=run.chem.setup


class CachedNative(Original):
    def __init__(self,cap=24000):
        # Four production jobs plus pilot calls remain below80000 in total.
        super().__init__(min(cap,18000))
        self.state=lru_cache(maxsize=16)(self.state)


def setup(native,d,j):
    original_setup(native,d,j)
    native.state.cache_clear()  # Anchor/inventories are implicit state inputs.


def configure():
    run.chem.old.Native=CachedNative;run.chem.setup=setup


def prepare():
    original=json.loads((OUT/'pilot.json').read_text());assert not original['eligible']
    assert not (OUT/'cache-plan.json').exists()
    write(OUT/'cache-plan.json',dict(classification='Counterexample candidate',
        failure='Original conservative pilot forecast1303.61s wall/3430.05 worker-seconds exceeds900/2400. Full production was not dispatched.',
        intervention='Memoize only EXACT native state(x,logT,y) calls within one fixed anchor. Clear the cache whenever chem.setup changes its implicit inventories. Store at most16 returned states per instance. Native calls with changed coordinates are never approximated or merged.',
        decision='Compare cached root against saved uncached roots at the same actual conserved states, keeping original inverse tolerances. Forecast remaining work at twice the slowest measured root in each region. No production unless original900s wall/2400s aggregate worker cap passes.',
        budget=dict(cached_pilot_seconds=int(90-original['seconds']),original_total_pilot_seconds=90,production_wall_seconds=900,
                    production_worker_seconds=2400,readout_seconds=90,native_calls_per_production_job=18000,workers=3),
        limits='Exact-call reuse accelerates the same native inverse; it does not turn sampled EOS states into a uniform enclosure.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),OUT/'plan.json',OUT/'pilot.json',OUT/'pilot-samples.json']}))


def pilot():
    assert not (OUT/'cached-pilot.json').exists();plan=json.loads((OUT/'cache-plan.json').read_text())
    start=time.monotonic();signal.alarm(plan['budget']['cached_pilot_seconds']);configure()
    model,ds,ats=run.load_inputs();old=json.loads((OUT/'pilot-samples.json').read_text());by={(r['kind'],r['it'],r['cell']):r for r in old}
    n=CachedNative(cap=1200);n.prefix=dict(n.prefix);durations={'atmosphere':[],'deep':[]};checks=[];hits=0
    for it in [0,8,16]:
        ids=ats[it]['active']
        for j in [int(ids[0]),int(ids[len(ids)//2]),int(ids[-1])]:
            tick=time.monotonic();r=run.atmosphere.root(n,model.flow,ats[it],j);durations['atmosphere'].append(time.monotonic()-tick)
            previous=by['atmosphere',it,j];err=abs(r['p']/previous['p']-1);assert err<1e-10
            checks.append(dict(kind='atmosphere',it=it,cell=j,uncached_pressure_relative=err))
    hits+=n.state.cache_info().hits
    for it,j in [(0,0),(8,9),(16,18)]:
        setup(n,model.bulk.eos.base.d,j);tick=time.monotonic();r=run.deep.root(n,ds[it],j);durations['deep'].append(time.monotonic()-tick)
        previous=by['deep',it,j];err=abs(r['raw'][1]/previous['raw'][1]-1);assert err<1e-10
        checks.append(dict(kind='deep',it=it,cell=j,uncached_pressure_relative=err));hits+=n.state.cache_info().hits
    prior=json.loads((OUT/'pilot.json').read_text());remaining=prior['remaining'];work={k:2*max(durations[k])*remaining[k] for k in remaining}
    forecast=work['deep']+work['atmosphere']/3+30
    result=dict(prior,eligible=forecast<900 and sum(work.values())<2400,root_seconds=durations,
        forecast_upper_wall_seconds=forecast,forecast_upper_worker_seconds=sum(work.values()),cache_hits=hits,
        equivalence_checks=checks,cached_pilot_native_calls=n.ion.calls,cached_pilot_seconds=time.monotonic()-start,
        cache_source_sha256=sha(__file__))
    write(OUT/'cached-pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def production():
    cached=json.loads((OUT/'cached-pilot.json').read_text());assert cached['eligible'] and cached['cache_source_sha256']==sha(__file__)
    configure()
    # Reuse the original production body with only its explicit forecast input
    # routed to the new accepted pilot. The failed pilot stays unchanged.
    source=inspect.getsource(run.production).replace("OUT/'pilot.json'","OUT/'cached-pilot.json'")
    scope=vars(run);exec(compile(source,__file__,'exec'),scope);scope['production']()


def audit():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.alarm(30)
    model,ds,ats=run.load_inputs();rows=json.loads((OUT/'production-samples.json').read_text());source=dict(np.load(OUT/'source-knots.npz'))
    checks=[];direct=source['trace_erg'].copy()*0;stress=direct.copy();pressure=direct.copy()
    unit=model.flow.eos.rho0*run.C**2*4*np.pi*model.m.RJ**2*model.m.vol
    for r in rows:
        it,j=r['it'],r['cell'];raw=np.array(r['raw'])
        if r['kind']=='deep':
            error=abs(raw[2]/ds[it]['target'][j]-1);assert error<1.01e-12
            dp=(raw[1]-model.bulk.eos.inventory[0,j]*ds[it]['xi'][j]-ds[it]['p'][j])*model.bulk.volume[j]
            pressure[it,j]=dp;direct[it,j]=-3*dp;stress[it,j]=-dp
        else:
            d=ats[it];D,S=d['U'][:2,j];v=r['v'];rho=r['rho'];p=r['p'];root=np.sqrt(1-v*v)
            E=model.flow.eos.cx*D*v*v/(root*(1+root))+(rho*raw[2]/run.C**2+p)/(1-v*v)-p
            error=abs(E/d['tau'][j]-1);assert error<1.01e-12
            momentum=abs(v*(model.flow.eos.cx*D+d['tau'][j]+p)-S)/max(abs(S),1e-300);assert momentum<1.01e-12
            dp=(raw[1]/(model.flow.eos.rho0*run.C**2)-d['p'][j])*unit[j];vv=d['v'][j]*v;k=19+j
            pressure[it,k]=dp;direct[it,k]=(vv-3)*dp;stress[it,k]=(vv-1)*dp
        checks.append(error)
    errors={k:float(np.max(abs(actual-source[k]))/max(np.max(abs(source[k])),1e-300))
        for k,actual in [('pressure_volume_erg',pressure),('trace_erg',direct),('metric_stress_erg',stress)]}
    assert max(errors.values())<1e-12
    initial_norm=float(np.sum(abs(source['trace_erg'][0])));assert initial_norm>0
    result=dict(classification='Counterexample candidate',passed=True,states=len(rows),maximum_rebuilt_energy_relative=max(checks),
        raw_native_source_reconstruction_relative=errors,initial_trace_correction_L1_erg=initial_norm,
        original_initial_offset_retained=True,seconds=time.monotonic()-start,new_native_calls=0,
        continuous_EOS_error_certified=False,final_charge_solved=False)
    write(OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    cap=3*1024**3;run.resource.setrlimit(run.resource.RLIMIT_AS,(cap,cap));signal.signal(signal.SIGALRM,run.flow.old.optical.timeout)
    try:globals()[sys.argv[1]]()
    except Exception as exc:
        write(OUT/f'cached-{sys.argv[1]}-failure.json',dict(error=repr(exc)));raise
