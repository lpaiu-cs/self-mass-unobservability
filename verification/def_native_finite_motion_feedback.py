"""Return the actual finite native-collision material motion to photons and GR."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,resource,signal,sys,time
import numpy as np
import sympy as sp
import def_native_compensated_finite_return as matter
import def_native_updated_joint_feedback as joint
import verify_native_compensated_finite_return as audit

native=matter.old.old.repaired;prior=native.response;base=joint.base
OUT=matter.OUT/'motion-feedback';RETURN=OUT/'return';GR=RETURN/'gr'
AMP=matter.AMP;C=matter.C;write=matter.write;sha=matter.sha;configure=matter.configure


class Response(joint.Response):
    def __init__(self,reference=128,steps=128):
        assert reference==128;base.mono.Response.__init__(self,reference)
        self.material=matter.Material(reference,steps)
        for key in self.g:
            if key.startswith('delta_'):self.g[key]=np.zeros_like(self.g[key])
        d=dict(np.load(matter.OUT/f'steps-{steps}-reference-{reference}.npz'))
        ids=[int(np.argmin(abs(d['t']-t))) for t in self.t];assert np.max(abs(d['t'][ids]-self.t))<1e-18
        self.motion=d['history_scaled'][ids].copy()
        self.energy_offset=(self.a-self.model.m.a0)*self.model.cx*C**2*self.motion[:,0]
        self.motion[:,2]-=self.energy_offset
        self.mechanical=self.motion[:,[2,3]]-self.material.transfer[:,[2,3]]
        self.xi=np.array([self.material.model.mech.xi@np.r_[0.,-np.cumsum(z[0,:self.nb])] for z in self.motion])
        self.velocity_jet_error=0.;self.mapping_error=0.;self.map_checks=[];self.angular=[];self.angular_times=[]
        self.defects={};self.zero=np.zeros_like(self.I[0])

    defect=prior.Response.defect

    def local(self,t):
        c=super().local(t);k=max(0,min(np.searchsorted(self.t,t,side='right')-1,15));f=(t-self.t[k])/(self.t[k+1]-self.t[k])
        a,b=self.defect(k),self.defect(k+1)
        for key,source in [('q','photon'),('qb','bound'),('qe','escape')]:
            delta=((1-f)*a[source]+f*b[source])/AMP
            c[key]+=delta if key=='qe' else delta/self.scale
        return c

    source=prior.Response.source
    run=FunctionType(joint.Response.run.__code__,dict(joint.Response.run.__globals__,OUT=OUT),argdefs=joint.Response.run.__defaults__)


def prepare():
    assert not OUT.exists();OUT.mkdir();RETURN.mkdir();GR.mkdir()
    assert json.loads((matter.GR/'result.json').read_text())['passed']
    paths=[Path(__file__),Path(matter.__file__),Path(joint.__file__),Path(base.__file__),Path(native.__file__),
           matter.OUT/'production.json',matter.GR/'result.json',prior.OUT/'result.json',native.OUT/'result.json']
    paths+=[matter.OUT/f'steps-{n}-reference-128.npz' for n in [64,128]]
    paths+=list(native.OUT.glob('point-*.npz'))+list((prior.OUT/'bank-128').glob('point-*.npz'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='550de4ab6',
        claim='Apply the actual finite native-collision material baryon/momentum/inventory and noncollisional transport to simultaneous photon/energy/H stages, retain the same native forcing exactly once, then return the new paired transfers to finite matter and actual GR.',
        decision='Measure whether feedback changes or invalidates the positive conditional charge, and how much photon/material mismatch remains. A small change is not a coupled error bound.',
        reuse='Same531 cells,8 angles,152 frequencies,17 saved corrected backgrounds,3.434431ms and64/128 paths on reference128. Match each material path to its photon clock. No EOS roots, background replay, third path, finer mesh or horizon.',
        scope='Actual finite matter supplies motion to the existing linear conservative radiation primitive maps. Native collision forcing stays fixed on the saved background. Zero extra geometric forcing. This is one waveform sweep, not fully nonlinear radiation or certified EOS derivatives.',
        budgets=dict(check_seconds=45,photon_pilot_seconds=90,photon_production_seconds=900,material_pilot_seconds=60,material_production_seconds=780,source_seconds=90,GR_seconds=90,CPU_threads=1,virtual_GiB=3),
        forecast='Measure equal-horizon4/8 photon prefixes and use the previous full photon mean-step cost as a late-cost floor; require2x summed forecast below900s. Material two-step prefixes use completed finite raw-call counts and2x summed forecast below780s. Resume accepted prefixes. Point/host cost changes remain extrapolated.',
        gates=dict(owner=1e-10,velocity=.0001,mapping=1e-10,linear=1e-12,balance=1e-8,time=.02,finite_resolution=.002,pressure_resolution=.002,angular=1e-12,GR_quadrature=.002,GR_independent=1e-9),
        stop='One full sweep only. Stop on original numerical gate or budget. No automatic tolerance relaxation, further waveform iteration, rescaling, refinement or new native support.',
        bindings={str(p):sha(p) for p in paths}))
    (OUT/'expanded-run.py').write_text(joint.runner)
    a,b,c=sp.symbols('native motion unknown');assert sp.diff(a+b+c,a)==1
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='The affine collision sum contains one copy of the fixed native defect; the paired gas owner uses its identical photon/bound/escape moments. No nonlinear or discretization theorem.'))


def check():
    start=time.monotonic();signal.alarm(45);configure();m=Response();rows=[]
    for k in [0,8,16]:
        t=m.t[k];c=m.local(t);motion=base.Response.local(m,t);d=m.defect(k);errors=[]
        for key,src in [('q','photon'),('qb','bound'),('qe','escape')]:
            expected=d[src]/AMP
            if key!='qe':expected=expected/m.scale
            errors.append(float(np.max(abs(c[key]-motion[key]-expected))/max(np.max(abs(expected)),1.)))
        gas=m.gas(c['q'],c['qb'],c['qe'])
        energy=float(abs(np.sum(c['q']*m.Eweight)+np.sum(c['qe'][1])+np.sum(gas[:,0]*m.eu))/max(np.sum(abs(c['q'])*m.Eweight),np.sum(abs(gas[:,0])*m.eu),1.))
        species=float(abs(np.sum(c['qb']*m.Nweight)-np.sum(gas[:,1]*m.nu))/max(np.sum(abs(c['qb'])*m.Nweight),1.))
        assert max(errors+[energy,species])<1e-10
        rows.append(dict(k=k,single_native_source=errors,energy=energy,species=species))
    assert m.velocity_jet_error<1e-4 and m.mapping_error<1e-10
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,velocity=m.velocity_jet_error,mapping=m.mapping_error,seconds=time.monotonic()-start));signal.alarm(0)


def photon_worker(n,label,limit=None,restart=None):
    started=time.monotonic();cpu=time.process_time();m=Response(128,n);r=m.run(n,label,limit,restart)
    r.update(worker_wall_seconds=time.monotonic()-started,worker_CPU_seconds=time.process_time()-cpu,
             peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    write(OUT/f'{label}.json',r);assert r['passed'];return r


def photon_pilot():
    assert json.loads((OUT/'check.json').read_text())['passed'];assert not (OUT/'pilot.json').exists()
    configure();started=time.monotonic();signal.alarm(90)
    rows=[photon_worker(n,f'pilot-{n}',k) for n,k in [(64,4),(128,8)]]
    before=json.loads((prior.OUT/'result.json').read_text())['paths'];estimates=[]
    for r in rows:
        old=next(p for p in before if p['steps']==r['steps'])
        cost=max(r['stepping_seconds']/r['new_steps'],old['stepping_seconds']/old['new_steps'])
        estimates.append(17*r['operator_point_seconds']/r['operator_points']+(r['steps']-r['completed_steps'])*cost+20)
    a=np.load(OUT/'pilot-64.npz')['moments'][-1,[0,1,2,3,5,6]];b=np.load(OUT/'pilot-128.npz')['moments'][-1,[0,1,2,3,5,6]]
    errors=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist()
    result=dict(classification='Counterexample candidate',rows=rows,equal_horizon=errors,forecast_each_seconds=estimates,upper_seconds=2*sum(estimates),
                eligible=max(errors)<.02 and 2*sum(estimates)<900,seconds=time.monotonic()-started)
    write(OUT/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def photon_production():
    assert json.loads((OUT/'pilot.json').read_text())['eligible'];assert not (OUT/'result.json').exists()
    configure();started=time.monotonic();signal.alarm(900)
    rows=[photon_worker(n,f'steps-{n}-reference-128',restart=f'pilot-{n}') for n in [64,128]]
    histories=[];angular=[]
    for n in [64,128]:
        d=dict(np.load(OUT/f'steps-{n}-reference-128.npz'));histories.append(d['moments'][:,[0,1,2,3,5,6]])
        gamma=1-1/np.sqrt(2);h=d['t'][-1]/n
        expected=np.ravel(h*(np.arange(n)[:,None]+np.array([gamma,1.])))
        assert np.max(abs(d['accepted_angular_times']-expected))<1e-18
        flux=(d['accepted_angular_luminosity']@((np.arange(4)*2+1)/32)).reshape(n,2)
        angular.append(float(abs(h*np.sum(flux*[1-gamma,gamma])-d['radial_ports'][-1,1,1])/max(h*np.sum(abs(flux)),1.)))
    errors=joint.stage.compare(*histories)
    result=dict(classification='Counterexample candidate',passed=bool(max(errors)<.02 and max(angular)<1e-12),paths=rows,time_comparison=errors,
                angular_port_relative=angular,seconds=time.monotonic()-started,actual_finite_motion_in_photons=True,full_nonlinear_radiation=False,coupled_fixed_point_verified=False,final_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


def photon_path(n,r):return OUT/f'steps-{n}-reference-{r}.npz'


class Material(matter.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps);p=dict(np.load(photon_path(steps,reference)))
        ids=[int(np.argmin(abs(p['t']-t))) for t in self.t];assert np.max(abs(p['t'][ids]-self.t))<1e-18
        c=p['collision_transfer'][ids]
        self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/AMP
    run=FunctionType(matter.Material.run.__code__,dict(matter.Material.run.__globals__,OUT=RETURN),argdefs=matter.Material.run.__defaults__)


material_worker=FunctionType(matter.worker.__code__,dict(matter.worker.__globals__,OUT=RETURN,Material=Material),argdefs=matter.worker.__defaults__)


def material_pilot():
    assert json.loads((OUT/'result.json').read_text())['passed'];assert not (RETURN/'pilot.json').exists()
    started=time.monotonic();signal.alarm(60);rows=[material_worker(n,f'pilot-{n}',2) for n in [64,128]]
    before=json.loads((matter.OUT/'production.json').read_text())['rows'];estimates=[]
    for r in rows:
        old=next(p for p in before if p['steps']==r['steps'])
        estimates.append(max(old['raw_owner_calls']*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+15,old['worker_wall_seconds']))
    result=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=estimates,upper_seconds=2*sum(estimates),eligible=2*sum(estimates)<780,seconds=time.monotonic()-started)
    write(RETURN/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


material_production=FunctionType(matter.production.__code__,dict(matter.production.__globals__,OUT=RETURN,worker=material_worker))
sources=FunctionType(matter.old.sources.__code__,dict(matter.old.sources.__globals__,OUT=RETURN,GR=GR,Material=Material,
    old=SimpleNamespace(previous=matter.old.old.previous,photon_path=photon_path)))


def charge():
    assert json.loads((RETURN/'sources.json').read_text())['passed'];configure();started=time.monotonic();signal.alarm(90)
    m=matter.old.old.previous.GRResponse()
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));wave=audit.charge.read(m,d,q)
        np.savez_compressed(GR/f'wave-{n}-g{q}.npz',**wave)
    fn=FunctionType(audit.main.__code__,dict(audit.main.__globals__,OUT=RETURN,GR=GR));fn()
    r=json.loads((GR/'result.json').read_text());old=json.loads((matter.GR/'result.json').read_text())
    r.update(material_motion_returned_to_photons=True,new_photon_transfer_returned_to_finite_material=True,
             previous_collision_endpoint=old['endpoint_collision_charge'],
             endpoint_feedback_change=r['endpoint_collision_charge']-old['endpoint_collision_charge'],seconds=time.monotonic()-started)
    write(GR/'result.json',r);print(json.dumps(r),flush=True);signal.alarm(0)


if __name__=='__main__':
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));signal.signal(signal.SIGALRM,native.forcing.history.flow.old.optical.timeout)
    action=sys.argv[1];started=time.monotonic()
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc),seconds=time.monotonic()-started));raise
