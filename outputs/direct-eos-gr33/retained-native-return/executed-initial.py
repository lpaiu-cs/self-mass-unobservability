"""Counterexample candidate: current native collisions AND pressure forces.

One constitutive-defect sweep on the saved retained-material background.
The atmosphere uses actual native face p,u, with the declared table acoustic
closure. This is neither native derivative certification nor a fixed point.
"""
from pathlib import Path
from types import FunctionType, MethodType, SimpleNamespace
import hashlib, inspect, io, json, resource, signal, sys, textwrap, time

ROOT=Path('outputs/direct-eos-gr33'); OUT=ROOT/'retained-native-return'
CURRENT=ROOT/'native-retained-completion'; EV=CURRENT/'evolution'
EOS=ROOT/'retained-eos-source'; FORCE=OUT/'collision'; HYDRO=OUT/'hydro'
PHOTON=OUT/'photons'; MATERIAL=OUT/'material'; GR=OUT/'gr'; METRIC=OUT/'metric'
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def replace(s,a,b):assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)
CAPS=dict(collision=90,hydro_pilot=90,hydro_production=600,bank=90,
          photon_pilot=90,photon_production=900,material_pilot=60,
          material_production=780,readout=180)
NO_ALARM=SimpleNamespace(alarm=lambda n:0,signal=signal.signal,SIGALRM=signal.SIGALRM)


def prepare():
    assert not OUT.exists()
    for p in [OUT,FORCE,HYDRO,PHOTON,MATERIAL,GR,METRIC/'corrected']:p.mkdir(parents=True,exist_ok=True)
    paths=[Path(__file__),EOS/'result.json',EOS/'production-samples.json',EOS/'plan.json',
        EV/'coupled-128.npz',EV/'source-128.npz',EV/'result.json',
        *[Path('verification')/n for n in ['evolve_native_retained_tail.py','def_retained_eos_source.py',
          'def_native_material_response.py','def_native_compensated_finite_return.py',
          'def_native_finite_collision_return.py','def_native_history_collision_forcing.py',
          'repair_native_collision_defect.py','def_native_collision_defect_response.py',
          'verify_native_updated_gr_return.py','def_native_collision_response.py']]]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='8450ad5d5',
        claim='Apply current native collision and conservative pressure-flux differences to actual photon/energy/H evolution and finite free material, then the represented GR charge.',
        decision='Measure whether native pressure forces or returned collision transfer change the surviving charge, and identify an actual failure before any longer or finer calculation.',
        equation='d(deltaQ)/dt=Ftable(Qsaved+deltaQ)-Ftable(Qsaved)+[Fnative(Qsaved)-Ftable(Qsaved)]+paired photon collision transfer.',
        scope='Existing531 cells,17 fine-background knots,64/128 response clocks and3.4344311179287023ms. Zero additional metric forcing for this isolated constitutive correction. Preserve the initial native pressure offset.',
        pressure='Reuse exact conserved cell roots; evaluate native EOS at reconstructed active faces. Reuse the SAME reconstruction, HLL face, donor and gravity operators. The shared face is evaluated once and enters both sides. Native p,u replace table values; the existing table gamma and deep face_K remain an explicitly declared acoustic closure, not certified native derivatives.',
        reuse='No background replay, new EOS bank, finer grid, smaller floor or longer horizon. Collision roots are reused exactly. Resume accepted response prefixes. Save accepted outgoing angular increments for later infinity return.',
        budgets=CAPS,total_max_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,new_native_face_call_cap=16000,
        forecast='Native face pilot at knots0,8,16: twice worst whole-knot cost for14 remaining plus measured setup and20s must fit600s. Photon equal-horizon4/8 prefixes and previous late cost floor must fit900s with2x margin. Finite material two-step prefixes and prior full raw-call counts must fit780s with2x margin.',
        gates=dict(native_pressure=.002,native_energy=.002,face_density=1e-10,conservation=1e-8,shared=1e-12,
            collision_positive=.002,collision_net=.02,number=1e-12,derivative_samples=1e-4,
            flux_resolution=.002,response_time=.02,pressure_resolution=.002,GR_quadrature=.002,independent_GR=1e-9),
        stop='Stop on the existing physical, arithmetic, support, time-comparison or measured-cost gate. Preserve failures. Do not weaken gates, reduce amplitude, add clocks, or replay completed work.',
        limitations='First constitutive-defect sweep on a saved finite history; subsequent reciprocal material-to-photon/GR feedback, uniform EOS derivatives, spatial convergence and final charge remain unproved.',
        bindings={p.as_posix():sha(p) for p in paths},native_bindings={**read(EOS/'plan.json')['inherited_bindings'],**read(EOS/'plan.json')['native_bindings']}))


def initialize():
    global np, current, retained, repaired, collision, forcing, matter, base, old, updated, C, AMP, Material, Response
    import numpy as np
    import evolve_native_retained_tail as current
    import def_retained_eos_source as retained
    import repair_native_collision_defect as repaired
    import def_native_compensated_finite_return as matter
    collision=repaired.response;forcing=repaired.forcing;updated=collision.updated;old=collision.old
    base=matter.old.old.physical.branch.base;C=matter.C;AMP=matter.AMP
    loader=np.load;saver=np.savez_compressed;cached={}
    def buffered(path,*args,**kw):
        if isinstance(path,(str,Path)) and Path(path).suffix=='.npz':
            p=Path(path)
            if p==EV/'coupled-128.npz':
                if str(p) not in cached:cached[str(p)]=p.read_bytes()
                data=cached[str(p)]
            else:data=p.read_bytes()
            return loader(io.BytesIO(data),*args,**kw)
        return loader(path,*args,**kw)
    def atomic(path,*args,**kw):
        if not isinstance(path,(str,Path)):return saver(path,*args,**kw)
        p=Path(path);s=io.BytesIO();saver(s,*args,**kw);tmp=p.with_suffix('.saving');tmp.write_bytes(s.getbuffer());tmp.replace(p)
    np.load=buffered;np.savez_compressed=atomic
    updated.configure();old.flow.OUT=EV;old.flow.Coupled=current.Capture;old.OUT=PHOTON
    original_metric=old.prior.metric.OUT/'corrected/metric-128-g8.npz'
    if not (METRIC/'corrected/metric-128-g8.npz').exists():
        z=dict(np.load(original_metric))
        for key in z:
            if key.startswith('delta_'):z[key]=np.zeros_like(z[key])
        assert np.max(abs(z['t']-np.load(EV/'coupled-128.npz')['snapshot_t']))<1e-18
        np.savez_compressed(METRIC/'corrected/metric-128-g8.npz',**z)
        write(METRIC/'binding.json',dict(classification='Imported from prior work',path=str(original_metric),sha256=sha(original_metric),zero_additional_geometry=True))
    old.prior.metric.OUT=METRIC;base.photons.OUT=PHOTON;forcing.OUT=FORCE
    retained.np=np;retained.run=current;retained.history=forcing.history;retained.load=np.load
    forcing.inputs=inputs
    # Reuse the material operator without loading obsolete photon histories.
    s=textwrap.dedent(inspect.getsource(base.Material.__init__))
    begin=s.index('    x=np.load(photons.OUT/');end=s.index('    self.rest=',begin)
    s=s[:begin]+'    self.transfer=np.zeros((17,4,self.n))\n'+s[end:]
    ns=dict(vars(base));exec(compile(s,__file__,'exec'),ns);base_init=ns['__init__']
    s=textwrap.dedent(inspect.getsource(matter.Material.__init__))
    s=replace(s,'super().__init__(reference,steps);m=self.model',
        'base_init(self,reference);self.steps=steps;self.finite_resolution=0.;self.finite_calls=0;m=self.model')
    ns=dict(vars(matter),base_init=base_init);exec(compile(s,__file__,'exec'),ns);material_init=ns['__init__']
    class MaterialImpl(matter.Material):
        def __init__(self,reference=128,steps=128,driven=True):
            material_init(self,reference,steps)
            self.hydro={}
            if driven:
                p=np.load(PHOTON/f'steps-{steps}-reference-128.npz');ids=[int(np.argmin(abs(p['t']-t))) for t in self.t]
                assert np.max(abs(p['t'][ids]-self.t))<1e-18;c=p['collision_transfer'][ids]
                self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/AMP
        def rhs(self,t,z,probe=1.):
            r,l,dt=super().rhs(t,z,probe);j,f,_,_=self.fields(t)
            for k,w in [(j,1-f),(j+1,f)]:
                if not w:continue
                if k not in self.hydro:self.hydro[k]=dict(np.load(HYDRO/f'point-{k}.npz'))
                h=self.hydro[k];r+=w*h['rate']/AMP;l+=w*h['ledger']/AMP;dt=min(dt,float(h['dt']))
            return r,l,dt
        run=FunctionType(matter.Material.run.__code__,dict(matter.Material.run.__globals__,OUT=MATERIAL),argdefs=matter.Material.run.__defaults__)
    Material=MaterialImpl
    # Original nonautonomous SDIRK stages and actual accepted outgoing ports.
    runner=updated.runner
    runner=replace(runner,"times=z['t'].tolist();", "self.angular_times=list(z['accepted_angular_times']);self.angular=list(z['accepted_angular_luminosity']);times=z['t'].tolist();")
    runner=runner.replace('transfer_history=transfer_history)', 'transfer_history=transfer_history,accepted_angular_times=self.angular_times,accepted_angular_luminosity=self.angular)')
    ns=dict(vars(updated.mono),OUT=PHOTON);exec(compile(runner,__file__,'exec'),ns)
    class ResponseImpl(collision.Response):
        def __init__(self,reference=128):super().__init__(reference);self.angular=[];self.angular_times=[]
        def boundary_ports(self,t,x):
            actual=x[-1]*self.scale
            self.angular_times.append(t);self.angular.append(AMP*2*np.pi*C*self.area[-1]*(actual[self.mu>0]@(self.num*self.E)))
            return super().boundary_ports(t,x)
        run=ns['run']
    Response=ResponseImpl
    (OUT/'expanded-material-init.py').write_text(s);(OUT/'expanded-photon-run.py').write_text(runner)


def inputs():
    m,ds,ats=retained.load_inputs();z=dict(np.load(EV/'coupled-128.npz'))
    rows=read(EOS/'production-samples.json');by={(r['kind'],r['it'],r['cell']):r for r in rows}
    return m,ds,ats,z,by


def collision_sources():
    fn=FunctionType(repaired.repair.__code__,dict(vars(repaired),OUT=FORCE,signal=NO_ALARM));fn()
    points=[dict(np.load(FORCE/f'point-{k}.npz')) for k in range(17)]
    positive=max(max(read(FORCE/f'bound-repaired/point-{k}.json')['positive_relative']) for k in range(17))
    t=np.array([p['t'] for p in points]);a=np.array([p['original_moments'] for p in points]);b=np.array([p['defect_moments'] for p in points])
    ratios={name:(np.trapezoid(np.sum(abs(b[:,:,sl]),axis=2),t,axis=0)/np.maximum(np.trapezoid(np.sum(abs(a[:,:,sl]),axis=2),t,axis=0),1.)).tolist()
        for name,sl in [('deep',slice(0,19)),('atmosphere',slice(19,None)),('all',slice(None))]}
    passed=positive<.002 and max(v for r in ratios.values() for v in r)<.02
    write(FORCE/'audit.json',dict(classification='Counterexample candidate',passed=passed,positive=positive,integrated_net=ratios));assert passed


def hydro_sources(pilot):
    start=time.monotonic();m=Material(driven=False);model=m.model;f=model.flow;by={(r['kind'],r['it'],r['cell']):r for r in read(EOS/'production-samples.json')}
    import verify_native_conserved_history_charge as cached
    high=cached.CachedNative(16000);high.prefix=dict(high.prefix);low=current.support.fast_native(16000);low.prefix=dict(low.prefix)
    rows=[];durations=[];states={};new=[];maxp=0.;maxu=0.;maxrho=0.
    for r in by.values():
        if r['kind']=='atmosphere':states[(float(np.log(r['rho'])),r['lt'],r['y'])]=np.array(r['raw'])
    native_calls0=high.ion.calls+low.ion.calls
    # ponytail: exact-coordinate cache is local to this saved history; no EOS interpolation is inferred from it.
    class NativeEOS:
        def __init__(self,parent):self.parent=parent
        def __getattr__(self,key):return getattr(self.parent,key)
        @property
        def y(self):return self.parent.y
        @y.setter
        def y(self,value):self.parent.y=value
        def evaluate(self,rho,lt):
            nonlocal maxp,maxu,maxrho
            values=np.array(self.parent.evaluate(rho,lt));y=np.broadcast_to(self.y,rho.shape)
            for j in np.flatnonzero(rho>=self.floor):
                key=(float(np.log(rho[j])),float(lt[j]),float(y[j]))
                if key not in states:
                    n=low if rho[j]<self.parent.original.floor else high;s=n.state(*key);states[key]=s['raw'].copy()
                    new.append(dict(coordinates=key,raw=s['raw'].tolist(),owner='low' if n is low else 'high',population_error=float(s['population_error'])))
                    assert s['population_error']<1e-12
                raw=states[key];p=raw[1]/(self.rho0*C*C);u=raw[2]/C**2
                maxp=max(maxp,abs(p/values[0,j]-1));maxu=max(maxu,abs(u-values[1,j])/max(abs(u),1e-300));maxrho=max(maxrho,abs(raw[0]/(rho[j]*self.rho0)-1))
                assert maxp<.002 and maxu<.002 and maxrho<1e-10,(key,maxp,maxu,maxrho)
                values[0,j]=p;values[1,j]=u
            assert high.ion.calls+low.ion.calls-native_calls0<=16000
            return tuple(values)
        def __call__(self,rho,lt):return self.evaluate(rho,lt)[:5]
    indices=[0,8,16] if pilot else [k for k in range(17) if not (HYDRO/f'point-{k}.json').exists()]
    setup=time.monotonic()-start
    if not pilot:
        p=read(HYDRO/'pilot.json');assert p['eligible']
        upper=2*max(p['point_seconds'])*len(indices)+setup+20
        write(HYDRO/'dispatch.json',dict(classification='Counterexample candidate',upper_seconds=upper,remaining=indices,eligible=upper<CAPS['hydro_production']))
        assert upper<CAPS['hydro_production']
    try:
        for k in indices:
            tick=time.monotonic();p=m.point(k);before=m.raw(k,np.zeros_like(p['Q']),np.zeros((5,m.n)),0.)
            theta,eta=before[3]['theta'],before[3]['eta'];prim=before[3]['primitive'].copy();pg=np.array([by['deep',k,j]['pressure'] for j in range(m.nb)])
            nt=theta+np.array([by['deep',k,j]['lt']-by['deep',k,j]['old_lt'] for j in range(m.nb)])
            for j in np.flatnonzero(p['active'][m.nb:]):
                r=by['atmosphere',k,int(j)];prim[:,j]=[r['rho'],r['v'],r['lt'],r['y']]
            oldprimitive=f.primitive;oldeos=f.eos;gas=model.bulk.eos.gas
            def native_gas(t,e):
                v=list(gas(t,e));v[0]=pg.copy();return v
            model.bulk.eos.gas=native_gas;f.eos=NativeEOS(oldeos);f.primitive=lambda U:tuple(prim.copy())
            f.join_state[2]+=nt[-1]-theta[-1]
            try:
                af,ag,dt,_=m.atmosphere(m.d['snapshot_U'][m.ids[k]],m.t[k]);df,dg,ddt=m.deep(m.d['snapshot_u'][m.ids[k]],nt,eta)
                factor=4*np.pi*model.m.RJ**2*f.eos.rho0*C*np.array([1,C*C,C*C,f.eos.nH]);af=af*factor[:,None];ag=ag*m.V[m.nb:]*f.eos.rho0*C*C
                shared=float(abs(df[1,-1]+C*model.base_momentum_flux[-1]-af[1,0])/max(abs(af[1,0]),1.));assert shared<1e-12
                af[1,0]=df[1,-1];ag[0]+=C*model.base_momentum_flux[-1]
                allshared=float(np.max(abs(df[:,-1]-af[:,0])/np.maximum(abs(af[:,0]),1.)));assert allshared<1e-12
                F=np.c_[df[:,:-1],af];G=np.r_[dg,ag]
            finally:f.primitive=oldprimitive;f.eos=oldeos;model.bulk.eos.gas=gas
            dF=F.astype(np.longdouble)-before[0].astype(np.longdouble);dG=G.astype(np.longdouble)-before[1].astype(np.longdouble)
            rate=-np.diff(dF,axis=1);rate[1]+=dG;ledger=dF[:,0]-dF[:,-1];ledger[1]+=dG.sum()
            balance=float(np.max(abs(rate.sum(1)-ledger)/np.maximum(np.sum(abs(rate),axis=1),1.)));assert balance<1e-8
            np.savez_compressed(HYDRO/f'point-{k}.npz',t=m.t[k],flux=dF,gravity=dG,rate=np.asarray(rate,float),ledger=np.asarray(ledger,float),dt=min(dt,ddt))
            duration=time.monotonic()-tick;durations.append(duration)
            row=dict(classification='Counterexample candidate',k=k,seconds=duration,shared=allshared,balance=balance,owner=p['owner_error'],rate_L1=np.asarray(np.sum(abs(rate),axis=1),float).tolist());assert row['owner']<1e-8
            write(HYDRO/f'point-{k}.json',row);rows.append(row)
            write(HYDRO/f'{"pilot" if pilot else "production"}-native-states.json',new)
        if pilot:
            upper=2*max(durations)*14+setup+20;result=dict(eligible=upper<600,upper_seconds=upper,point_seconds=durations)
        else:result=dict(complete=all((HYDRO/f'point-{k}.json').exists() for k in range(17)))
        result.update(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start,setup_seconds=setup,
            new_native_states=len(new),native_calls=high.ion.calls+low.ion.calls-native_calls0,maximum_pressure_relative=maxp,maximum_energy_relative=maxu,maximum_density_relative=maxrho,native_acoustic_derivatives_certified=False)
        write(HYDRO/f'{"pilot" if pilot else "production"}.json',result)
    finally:write(HYDRO/f'{"pilot" if pilot else "production"}-native-states.json',new)


def bank():
    s=inspect.getsource(old.bank).replace('for reference in [128,64]:','for reference in [128]:')
    ns=dict(vars(old),OUT=PHOTON,Response=collision.Coefficients,signal=NO_ALARM);exec(compile(s,__file__,'exec'),ns);ns['bank']()


def photon(pilot):
    fn=collision.pilot if pilot else collision.production
    ns=dict(fn.__globals__,OUT=PHOTON,Response=Response,configure=lambda:None,signal=NO_ALARM)
    ns['forcing_check']=FunctionType(collision.forcing_check.__code__,dict(vars(collision),OUT=PHOTON))
    FunctionType(fn.__code__,ns)()


def material(pilot):
    assert read(HYDRO/'production.json')['complete'] and read(PHOTON/'result.json')['passed']
    worker=FunctionType(matter.old.worker.__code__,dict(matter.old.worker.__globals__,OUT=MATERIAL,Material=Material,configure=lambda:None),argdefs=matter.old.worker.__defaults__)
    fn=matter.old.pilot if pilot else matter.old.production
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=MATERIAL,worker=worker,signal=NO_ALARM))()


def template():
    d=dict(np.load(EV/'source-128.npz'));nt=len(d['t'])
    return {k:(v[::8] if v.ndim and v.shape[0]==nt else v.copy()) for k,v in d.items()}


def readout():
    import verify_native_stage_energy_charge as charge
    import sympy as sp
    owner=matter.old.old.previous.previous.old.source_owner
    s=owner.source.replace('[[64,128],[128,128],[128,64]]','[[64,128],[128,128]]')
    s=s.replace('background=compare(histories[2],histories[1]),','').replace(',stress_background=compare(allstress[2],allstress[1])','')
    s=replace(s,"base=dict(np.load(wave.base.OUT/f'source-{ref}.npz'))",'base=template()')
    ns=dict(matter.old.old.previous.source_scope,OUT=MATERIAL,GR=GR,Material=Material,
        stress=SimpleNamespace(pressure=matter.old.pressure),signal=NO_ALARM,template=template,
        photons=SimpleNamespace(path=lambda n,r:PHOTON/f'steps-{n}-reference-{r}.npz'),
        material_path=lambda n,r:MATERIAL/f'steps-{n}-reference-{r}.npz')
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-source.py').write_text(s);ns['sources']()
    assert read(MATERIAL/'sources.json')['passed']
    model=charge.gr.Response();waves={}
    for n,q in [(128,8),(64,8),(128,4)]:
        d=dict(np.load(GR/f'source-{n}-reference-128.npz'));waves[n,q]=charge.read(model,d,q)
        np.savez_compressed(GR/f'wave-{n}-g{q}.npz',**waves[n,q])
    fine=waves[128,8]['free_scalar'];norm=max(np.max(abs(fine)),1e-300)
    terr=float(np.max(abs(fine-waves[64,8]['free_scalar']))/norm);qerr=float(np.max(abs(fine-waves[128,4]['free_scalar']))/norm)
    d=dict(np.load(GR/'source-128-reference-128.npz'));direct,coordinate=charge.independent.direct(model,d,8)
    ierr=float(abs(direct-waves[128,8]['direct_scalar'][-1])/max(np.max(abs(waves[128,8]['direct_scalar'])),1e-300))
    total=template();p=dict(np.load(EOS/'source.npz'));keys=forcing.history.physical.capture.KEYS+['metric_stress_erg']
    for key in keys:total[key]=total[key]+p[key][::8]+d[key]
    for key in ['inner_cumulative_energy_erg','outer_cumulative_energy_erg']:total[key]=total[key]+d[key]
    applied=charge.read(model,total,8);baseline=template()
    for key in keys:baseline[key]=baseline[key]+p[key][::8]
    before=charge.read(model,baseline,8);linear=float(np.max(abs(applied['free_scalar']-before['free_scalar']-fine))/max(np.max(abs(before['free_scalar'])),1e-300))
    np.savez_compressed(GR/'applied-source.npz',**total);np.savez_compressed(GR/'applied-charge.npz',**applied,previous=before['free_scalar'],correction=fine)
    # A conservative face enters adjacent cells with exactly opposite signs.
    left,shared,right,g0,g1=sp.symbols('left shared right g0 g1')
    assert sp.expand((left-shared+g0)+(shared-right+g1)-(left-right+g0+g1))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Shared-face telescoping only; no nonlinear or uniform EOS error theorem.'))
    passed=terr<.02 and qerr<.002 and ierr<1e-9 and linear<1e-10
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=passed,
        native_return_compact_endpoint=float(fine[-1]),previous_same_cadence_endpoint=float(before['free_scalar'][-1]),
        applied_compact_endpoint=float(applied['free_scalar'][-1]),time_relative=terr,quadrature_relative=qerr,
        independent_direct_relative=ierr,inverse_radius_residual=float(coordinate),linear_application_relative=linear,
        native_collision_applied=True,native_pressure_force_applied=True,free_material_evolved=True,represented_GR_applied=True,
        actual_material_motion_returned_to_photons=False,new_angular_emission_returned_to_infinity=False,
        native_acoustic_derivatives_certified=False,coupled_fixed_point_verified=False,final_charge_solved=False,full_goal_complete=False))
    assert passed


def execute(action):
    plan=read(OUT/'plan.json')
    for p,h in {**plan['bindings'],**plan['native_bindings']}.items():assert sha(p)==h,p
    initialize()
    if action=='collision':collision_sources()
    elif action.startswith('hydro_'):hydro_sources(action.endswith('pilot'))
    elif action=='bank':bank()
    elif action.startswith('photon_'):photon(action.endswith('pilot'))
    elif action.startswith('material_'):material(action.endswith('pilot'))
    elif action=='readout':readout()
    else:raise ValueError(action)


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    else:
        cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));start=time.monotonic()
        def timeout(*_):raise TimeoutError(action+' wall budget')
        signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,CAPS[action])
        try:execute(action)
        except Exception as exc:
            write(OUT/f'{action}-failure.json',dict(classification='Counterexample candidate',error=repr(exc),seconds=time.monotonic()-start));raise
        finally:
            signal.setitimer(signal.ITIMER_REAL,0)
            write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,cap_seconds=CAPS[action],peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
