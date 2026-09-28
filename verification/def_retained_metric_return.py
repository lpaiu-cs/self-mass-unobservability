"""Apply the current retained GR fields to photons and free material."""
from pathlib import Path
from types import FunctionType
import inspect,json,resource,signal,sys,time
import numpy as np
import def_retained_metric_fields as fields
import def_native_updated_material_return as matter
import def_native_stage_aligned_material as aligned

OUT=fields.OUT;PHOTON=fields.PHOTON;MATERIAL=fields.MATERIAL;GR=fields.GR
read,write,sha=fields.read,fields.write,fields.sha;AMP=fields.AMP;C=fields.C


def prepare():
    assert read(fields.METRIC/'result.json')['passed']
    assert not (OUT/'response-plan.json').exists()
    files=[Path(__file__),Path(fields.__file__),Path(matter.__file__),Path(aligned.__file__),
           Path(matter.branch.__file__),Path(matter.branch.base.__file__),OUT/'plan.json',fields.METRIC/'result.json']
    files += list((fields.NATIVE/'photons/bank-128').glob('*.npz'))
    for n in [64,128]:files += [fields.METRIC/f'corrected/metric-{n}-g8.npz',fields.NATIVE/f'photons/steps-{n}-reference-128.json',fields.BEFORE/f'material/steps-{n}-reference-128.json']
    write(OUT/'response-plan.json',dict(classification='Counterexample candidate',
        claim='Consume the current full-source metric in the simultaneous SDIRK photon/thermal/H equation and the free conservative material force/work equation.',
        source='Original radial/angular/frequency GR forcing plus geometric collision volume/lapse terms. No duplicate native constitutive defect; the applied metric already includes its source.',
        material='Same branch-aware directional conservative owner, corrected SSP left closing interval, deep radial kinetic stress2T and canonical geometry subtraction. Retain small-state1e-6, donor1percent, directional0.2percent and original support gates.',
        baseline='First geometric variation on the saved corrected-EOS background. Existing native pressure/collision response and Phase151 finite motion remain the accepted baseline for charge addition, but their native Jacobian is not substituted for the existing table derivative.',
        output='Save separate geometric photon and free-material increments, accepted angular packets and spectral-work ledger. Apply that increment once to the Phase151 direct-readout charge.',
        budgets=fields.CAPS,gates=read(OUT/'plan.json')['gates'],
        limitations='One actual metric return, not a proof of convergence or nonlinear Einstein closure. No new native evaluation, bank, clock, background or domain.',
        bindings={str(p):sha(p) for p in files}))


def initialize():
    global Response,Material
    fields.initialize();prior=fields.previous.prior;old=prior.old
    old.prior.metric.OUT=fields.METRIC;matter.branch.base.photons.OUT=PHOTON
    # Keep the existing collision banks on the current physical background.
    assert old.OUT==fields.NATIVE/'photons'
    class ResponseImpl(prior.updated.Response):
        def __init__(self,steps):
            super().__init__(128);self.angular=[];self.angular_times=[]
            self.g=dict(np.load(fields.METRIC/f'corrected/metric-{steps}-g8.npz'))
            assert np.array_equal(self.g['t'],self.t) and np.any(self.g['delta_lambda'])
        def boundary_ports(self,t,x):
            j=np.clip(np.searchsorted(self.t,t,side='left')-1,0,15);f=(t-self.t[j])/(self.t[j+1]-self.t[j])
            I=(1-f)*self.I[j]+f*self.I[j+1]
            speed=((1-f)*self.g['delta_log_speed'][j]+f*self.g['delta_log_speed'][j+1])/AMP
            actual=x*self.scale+I*speed[:,None,None]
            self.angular_times.append(t)
            self.angular.append(AMP*2*np.pi*C*self.area[-1]*(actual[-1,self.mu>0]@(self.num*self.E)))
            return super().boundary_ports(t,x)
        run=FunctionType(prior.Response.run.__code__,dict(prior.Response.run.__globals__,OUT=PHOTON),argdefs=prior.Response.run.__defaults__)
    Response=ResponseImpl
    matter.photon_path=lambda n,r:PHOTON/f'steps-{n}-reference-{r}.npz'
    class MaterialImpl(matter.Material):
        closing_stage=False
        interval_fields=aligned.ns['fields']
        fields=aligned.Material.fields
        end_rhs=aligned.Material.end_rhs
        directional_rhs=aligned.ns['rhs']
        rhs=aligned.Material.rhs
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps)
            self.metric=dict(np.load(fields.METRIC/f'corrected/metric-{steps}-g8.npz'))
            assert np.array_equal(self.metric['t'],self.t) and np.any(self.metric['delta_lambda'])
        run=FunctionType(aligned.Material.run.__code__,dict(aligned.Material.run.__globals__,OUT=MATERIAL),argdefs=aligned.Material.run.__defaults__)
    Material=MaterialImpl


def photon(pilot):
    start=time.monotonic();rows=[]
    if not pilot:
        p=read(PHOTON/'pilot.json');revision=read(OUT/'response-execution-plan.json')
        assert p['passed'] and revision['eligible'] and p['upper_remaining_seconds']<revision['stage_caps']['photon_production']
    for n in [64,128]:
        m=Response(n);label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        rows.append(m.run(n,label,n//16 if pilot else None,restart=None if pilot else f'pilot-{n}'))
    keys=[0,1,2,3,5,6]
    a=np.load(PHOTON/('pilot-64.npz' if pilot else 'steps-64-reference-128.npz'))['moments'][:,keys]
    b=np.load(PHOTON/('pilot-128.npz' if pilot else 'steps-128-reference-128.npz'))['moments'][:,keys]
    error=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-300)).tolist()
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(error)<.02,
        rows=rows,time_comparison=error,seconds=time.monotonic()-start,current_GR_applied=True)
    if pilot:
        forecast=[]
        for r in rows:
            p=read(fields.NATIVE/f"photons/steps-{r['steps']}-reference-128.json")
            step=max(r['stepping_seconds']/r['new_steps'],p['stepping_seconds']/p['new_steps'])
            forecast.append(17*r['operator_point_seconds']/r['operator_points']+(r['steps']-r['completed_steps'])*step+20)
        result.update(upper_remaining_seconds=2*sum(forecast),eligible=result['passed'] and 2*sum(forecast)<fields.CAPS['photon_production'])
    write(PHOTON/('pilot.json' if pilot else 'result.json'),result);print(json.dumps(result),flush=True)
    assert result['eligible' if pilot else 'passed']


def material(pilot):
    assert read(PHOTON/'result.json')['passed'];start=time.monotonic();rows=[]
    if not pilot:assert read(MATERIAL/'pilot.json')['eligible']
    for n in [64,128]:
        mark=time.monotonic();m=Material(128,n);label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        row=m.run(n,label,2 if pilot else None,restart=None if pilot else f'pilot-{n}')
        row.update(worker_wall_seconds=time.monotonic()-mark,physical_branch_ratio=m.physical_branch_ratio,
            maximum_owner_error=max(p['owner_error'] for p in m.cache.values()))
        row['passed']=row['passed'] and row['physical_branch_ratio']<.01 and row['maximum_owner_error']<1e-8
        write(MATERIAL/f'{label}.json',row);rows.append(row);assert row['passed'],row
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start)
    if pilot:
        forecasts=[]
        for r in rows:
            p=read(fields.BEFORE/f"material/steps-{r['steps']}-reference-128.json")
            # Directional Richardson uses two states at each knot; the prior finite owner used one.
            calls=2*p['raw_owner_calls']
            forecasts.append(calls*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+5)
        result.update(upper_remaining_seconds=2*sum(forecasts),eligible=2*sum(forecasts)<fields.CAPS['material_production'])
    write(MATERIAL/('pilot.json' if pilot else 'production.json'),result);print(json.dumps(result),flush=True)
    if pilot:assert result['eligible']


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    else:
        receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        def timeout(*_):raise TimeoutError('Registered response budget')
        revision=read(OUT/'response-execution-plan.json') if (OUT/'response-execution-plan.json').exists() else None
        cap=(revision['stage_caps'] if revision else fields.CAPS)[action]
        signal.signal(signal.SIGALRM,timeout);signal.alarm(cap);start=time.monotonic();cpu=time.process_time();error=None
        try:
            for p,h in read(OUT/'response-plan.json')['bindings'].items():
                target=OUT/'registered-response-producer.py' if revision and Path(p).name==Path(__file__).name else p
                assert sha(target)==h,p
            if revision:
                assert sha(__file__)==revision['source_sha256']
                for p,h in revision['bindings'].items():assert sha(p)==h,p
            initialize();globals()[action.split('_')[0]](action.endswith('pilot'))
        except Exception as exc:error=repr(exc);raise
        finally:
            signal.alarm(0);write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=sha(__file__)))
