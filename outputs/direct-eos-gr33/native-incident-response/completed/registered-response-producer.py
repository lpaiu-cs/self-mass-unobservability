"""Apply the external scalar input to the existing retained coupled owners."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,resource,sys,time
import numpy as np
import def_native_incident_drive as drive
import def_retained_metric_return as old

OUT=drive.OUT;PHOTON=drive.PHOTON;MATERIAL=drive.MATERIAL;GR=drive.GR
read,write,sha=drive.read,drive.write,drive.sha;AMP=old.AMP;C=drive.C;LD=np.longdouble


def initialize():
    global Response,Material
    for owner in [old,old.fields]:
        owner.OUT=OUT;owner.PHOTON=PHOTON;owner.MATERIAL=MATERIAL;owner.GR=GR
    old.fields.METRIC=drive.METRIC
    old.initialize();PhotonBase=old.Response;MaterialBase=old.Material
    class DrivenPhoton(PhotonBase):
        def __init__(self,steps):
            super().__init__(steps);self.driver=drive.Driver(8);self.drive_scale=1.
        def set_stage(self,t):
            self.g={k:v*self.drive_scale for k,v in self.driver.view(t,self.t).items()}
        def source(self,t):
            self.set_stage(t);return super().source(t)
        def local(self,t):
            self.set_stage(t);return super().local(t)
        def boundary_ports(self,t,x):
            self.set_stage(t);return super().boundary_ports(t,x)
        def frequency(self,I,omega):
            a=super().frequency(I,omega);b=super().frequency(I,-omega)
            return (*[(x-y)/2 for x,y in zip(a[:4],b[:4])],max(a[4],b[4]))
    Response=DrivenPhoton
    class DrivenMaterial(MaterialBase):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps);self.driver=drive.Driver(8)
        def fields(self,t):
            k=int(np.argmin(abs(self.t-t)))
            if abs(self.t[k]-t)<=16*np.finfo(float).eps*self.t[-1]:t=self.t[k]
            self.metric=self.driver.view(t,self.t,side='left' if self.closing_stage else 'right')
            return self.interval_fields(t)
        def rhs(self,t,z,probe=1.):return super().rhs(t,z,probe*8)
    Material=DrivenMaterial


def prepare():
    assert read(drive.METRIC/'result.json')['passed']
    assert not (OUT/'response-plan.json').exists()
    files=[Path(__file__),Path(drive.__file__),Path(old.__file__),Path(old.matter.__file__),Path(old.aligned.__file__),
           drive.METRIC/'result.json',OUT/'lapse-repair-plan.json']
    files += [drive.FIELDS/f'born-g{q}.npz' for q in [4,8]]
    write(OUT/'response-plan.json',dict(classification='Counterexample candidate',
        claim='Consume the declared incoming scalar/GR fields at actual photon and material stage times, save accepted angular packets and actual free-material response, and return their charge source.',
        null='Zero drive gives exactly zero compensated source on the same evolving retained background. Positive scaling and sign reversal check the declared linear source; no new nonlinear unforced background run is needed.',
        spectral_change='Odd conservative frequency forcing at the zero baseline frequency drift, explicitly replacing the old one-sided upwind artificial-diffusion direction. Keep number, energy, endpoint and metric-work ledgers. This is not a spectral continuum certificate.',
        material='Use the established8x Richardson arithmetic probe for both clocks; verify its4/8 and8/16 directions on actual saved endpoint states. These probes do not change eta=1e-30 or the physical incoming packet.',
        budgets=drive.CAPS,gates=read(OUT/'plan.json')['gates'],
        scope='Fine retained background and existing native-anchored interpolation/Jacobian. No uniform native derivative certificate, full nonlinear driven star, stationary transfer function or companion-matched claim.',
        bindings={str(p):sha(p) for p in files}))


def source_check():
    initialize();m=Response(128);rows=[];arrays={};start=time.monotonic()
    for t in [m.t[-1]/3,2*m.t[-1]/3]:
        samples={}
        for scale in [0.,1.,-1.,2.]:
            m.drive_scale=scale;s,l,e=m.source(t);c=m.local(t)
            values=[s,l,c['q'],c['qb'],c['qe']]
            samples[scale]=values
            if scale==0:assert all(np.count_nonzero(v)==0 for v in values)
            assert e<1e-12
        errors=[]
        for j in range(5):
            norm=max(np.max(abs(samples[1.][j])),1e-290)
            errors += [float(np.max(abs(samples[1.][j]+samples[-1.][j]))/norm),float(np.max(abs(2*samples[1.][j]-samples[2.][j]))/norm)]
        rows.append(dict(time=t,null_exact=True,odd_and_scaling=max(errors)));assert max(errors)<1e-12
    write(OUT/'source-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        seconds=time.monotonic()-start,scope='The declared compensated finite source is homogeneous in the incident amplitude; not a nonlinear response or native-Jacobian certificate.'))


def photon(pilot):
    assert read(OUT/'source-check.json')['passed'];initialize();start=time.monotonic();rows=[]
    if not pilot:assert read(PHOTON/'execution-plan.json')['eligible']
    for n in [64,128]:
        m=Response(n);label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        rows.append(m.run(n,label,n//16 if pilot else None,restart=None if pilot else f'pilot-{n}'))
    keys=[0,1,2,3,5,6]
    a=np.load(PHOTON/('pilot-64.npz' if pilot else 'steps-64-reference-128.npz'))['moments'][:,keys]
    b=np.load(PHOTON/('pilot-128.npz' if pilot else 'steps-128-reference-128.npz'))['moments'][:,keys]
    errors=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1e-290)).tolist()
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(errors)<.02,
        rows=rows,time_comparison=errors,seconds=time.monotonic()-start,actual_external_GR_input_consumed=True,
        additional_material_motion=False,full_goal_complete=False)
    if pilot:
        forecasts=[]
        for r in rows:
            prev=read(old.fields.NATIVE/f"photons/steps-{r['steps']}-reference-128.json")
            step=max(r['stepping_seconds']/r['new_steps'],prev['stepping_seconds']/prev['new_steps'])
            forecasts.append(17*r['operator_point_seconds']/r['operator_points']+(r['steps']-r['completed_steps'])*step+25)
        result['upper_remaining_seconds']=2*sum(forecasts);result['eligible']=result['passed'] and result['upper_remaining_seconds']<1000
    write(PHOTON/('pilot.json' if pilot else 'result.json'),result);print(json.dumps(result),flush=True)
    if pilot and result['eligible']:
        write(PHOTON/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,cap_seconds=1000,
            forecast=result['upper_remaining_seconds'],source_sha256=sha(__file__),driver_sha256=sha(drive.__file__)))
    assert result['eligible' if pilot else 'passed'],result


def material(pilot):
    assert read(PHOTON/'result.json')['passed'];initialize();start=time.monotonic();rows=[]
    if not pilot:assert read(MATERIAL/'execution-plan.json')['eligible']
    for n in [64,128]:
        mark=time.monotonic();m=Material(128,n);label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        row=m.run(n,label,2 if pilot else None,restart=None if pilot else f'pilot-{n}')
        row.update(worker_wall_seconds=time.monotonic()-mark,physical_branch_ratio=m.physical_branch_ratio,
            maximum_owner_error=max(p['owner_error'] for p in m.cache.values()))
        if not pilot:
            z=np.load(MATERIAL/f'{label}.npz')['delta_scaled'];rates=[m.rhs(m.t[-1],z,v)[0] for v in [.5,1.,2.]]
            norm=np.maximum(np.sum(abs(rates[1]),axis=1),1.)
            row['endpoint_probe_4_8_16']=[(np.sum(abs(x-rates[1]),axis=1)/norm).tolist() for x in [rates[0],rates[2]]]
            row['passed']=row['passed'] and np.max(row['endpoint_probe_4_8_16'])<.002
        row['passed']=row['passed'] and row['physical_branch_ratio']<.01 and row['maximum_owner_error']<1e-8
        write(MATERIAL/f'{label}.json',row);rows.append(row);assert row['passed'],row
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start)
    if pilot:
        forecasts=[]
        for r in rows:
            path=Path('retained-metric-return152-work')/('material-resolved' if r['steps']==128 else 'material')/f"steps-{r['steps']}-reference-128.json"
            prev=read(path);forecasts.append(prev['raw_owner_calls']*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+5)
        result['upper_remaining_seconds']=2*sum(forecasts);result['eligible']=result['upper_remaining_seconds']<500
    write(MATERIAL/('pilot.json' if pilot else 'production.json'),result);print(json.dumps(result),flush=True)
    if pilot:
        if result['eligible']:write(MATERIAL/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,
            cap_seconds=500,upper_remaining_seconds=result['upper_remaining_seconds'],source_sha256=sha(__file__),driver_sha256=sha(drive.__file__)))
        assert result['eligible'],result


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','source_check','photon_pilot','photon_production','material_pilot','material_production']
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));drive.native.deadline(drive.CAPS.get(action,60))
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action.startswith(('photon_','material_')):globals()[action.split('_')[0]](action.endswith('pilot'))
        else:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
        write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
