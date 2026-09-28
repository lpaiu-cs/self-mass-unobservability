"""Apply native low-density material support to the actual coupled owner."""
from pathlib import Path
from types import FunctionType,MethodType
import io,json,sys,time,signal
import numpy as np
import native_tail_supported_temperature as support
import def_native_updated_gr_return as previous
import def_native_cold_coupling as cold

OUT=support.OUT;BANK=support.BANK;EV=support.EV;GR=support.GR
write,read=support.write,support.read
def sha(p):return support.sha(Path(p))
C=previous.C;flow=previous.flow;physical=previous.physical
BaseEOS=cold.prior.optical.ex.EOS;BaseSpectrum=cold.prior.optical.Spectrum


class EOS(BaseEOS):
    def __init__(self,old):
        self.__dict__.update(old.__dict__);self.original=old
        self.low=object.__new__(BaseEOS);fn=BaseEOS.__init__
        FunctionType(fn.__code__,dict(fn.__globals__,OUT=BANK))(self.low)
        self.floor=self.low.floor

    def evaluate(self,rho,lt):
        low=rho<self.original.floor;y=np.broadcast_to(self.y,rho.shape);result=np.zeros((7,len(rho)))
        for mask,model in [(low,self.low),(~low,self.original)]:
            if mask.any():model.y=y[mask];result[:,mask]=model.evaluate(rho[mask],lt[mask])
        return tuple(result)

    def reactions(self,*_):raise AssertionError('Use actual coupled photons')


class Spectrum(BaseSpectrum):
    def __init__(self,old):
        self.__dict__.update(old.__dict__);self.original=old
        self.low=object.__new__(BaseSpectrum);self.low.__dict__.update(old.__dict__)
        c=self.low;c.d=d=dict(np.load(BANK/'spectrum-bank.npz'));c.f=[];c.rev=[]
        for k,y in enumerate(self.ys):
            T=np.exp(d['lt'])[None,:,None];lf=d['log_fraction'][k]
            pref=lf+(self.binding[0]-self.binding)[None,None,:]/(cold.prior.optical.ex.K*T)
            reverse=lf+np.log(y/(1-y))+d['affinity'][k,:,:,None]-self.binding[None,None,:]/(cold.prior.optical.ex.K*T)
            c.f.append([cold.prior.optical.RectBivariateSpline(d['x'],d['lt'],pref[:,:,j]) for j in range(10)])
            c.rev.append([cold.prior.optical.RectBivariateSpline(d['x'],d['lt'],reverse[:,:,j]) for j in range(10)])

    def levels(self,rho,lt,y):
        low=rho<self.rho0*np.exp(self.low.d['x'][-1]);f=np.zeros((len(rho),10));r=np.zeros_like(f)
        for mask,model in [(low,self.low),(~low,self.original)]:
            if mask.any():f[mask],r[mask]=model.levels(rho[mask],lt[mask],y[mask])
        return f,r


def atomic_npz(path,**values):
    # One mounted-filesystem write replaces many ZIP header seeks. A killed
    # write cannot replace the previous complete checkpoint or angular port.
    stream=io.BytesIO();np.savez_compressed(stream,**values);tmp=path.with_suffix('.saving')
    tmp.write_bytes(stream.getbuffer());tmp.replace(path)


def load_npz(path):return np.load(io.BytesIO(path.read_bytes()))


runner=previous.runner
anchor="""        tmp=path.with_suffix('.tmp')
        with tmp.open('wb') as stream:np.savez_compressed(stream,**values)
        tmp.replace(path)"""
# The original save() is nested inside run_capture().
anchor='\n'.join('    '+line for line in anchor.splitlines())
assert runner.count(anchor)==1
runner=runner.replace(anchor,'            atomic_npz(path,**values)')
for key in ['moments','t','ports','discard','escape']:
    runner=runner.replace("list(saved['"+key+"'])","list(saved['"+key+"'][:begin+1])")
runner=runner.replace('z=np.load(checkpoint)','z=load_npz(checkpoint)').replace('saved=np.load(sidecar)','saved=load_npz(sidecar)')
runner=runner.replace("np.savez_compressed(OUT/f'source-{steps}.npz',**source)","atomic_npz(OUT/f'source-{steps}.npz',**source)")
scope=dict(previous.capture_run.__globals__,OUT=EV,atomic_npz=atomic_npz,load_npz=load_npz)
exec(compile(runner,__file__,'exec'),scope);capture_run=scope['run_capture']


class Capture(previous.Capture):
    def __init__(self):
        super().__init__();self.flow.eos=EOS(self.flow.eos);self.spectrum=Spectrum(self.spectrum);self.tail=[]
        fn=self.flow.primitive.__func__
        self.flow.primitive=MethodType(FunctionType(fn.__code__,dict(fn.__globals__,OUT=EV),argdefs=fn.__defaults__,closure=fn.__closure__),self.flow)

    def run(self,steps,count=None,resume=False):
        self.steps=steps;self.port_file=EV/f'accepted-ports-{steps}.npz'
        self.angular=list(load_npz(self.port_file)['angular_luminosity']) if resume else []
        if resume:
            begin=int(load_npz(EV/f'checkpoint-{steps}.npz')['completed'])
            assert len(self.angular)>=begin
            self.angular=self.angular[:begin];self.tail=list(read(EV/f'tail-{steps}.json')['rows'])[:begin+1]
        try:
            result=capture_run(self,steps,count,resume);result['cost']=self.cost.copy()
            write(EV/f'cost-{steps}-{result["completed"]}.json',result)
            if result['complete']:self.export()
            return result
        finally:write(EV/f'tail-{steps}.json',dict(classification='Counterexample candidate',rows=self.tail))

    def moment(self,U,I,xb,u,theta,eta):
        start=time.monotonic();value=super(previous.Capture,self).moment(U,I,xb,u,theta,eta);k=len(self.times)
        self.cost['physics']+=time.monotonic()-start
        if k%(self.steps//16)==0:
            start=time.monotonic()
            atomic_npz(EV/f'state-{self.steps}-{k}.npz',U=U,I=I,bulk_I=xb*self.bulk.scale,u=u,
                theta=theta,eta=eta,h=self.h,Pi=self.Pi,mass=self.mass,j=self.j)
            self.cost['state_write']+=time.monotonic()-start;self.cost['state_count']+=1
        f=self.flow;m=self.m;seed=f.seed.copy();oldy=np.copy(f.eos.y)
        try:
            rho,v,lt,y=f.primitive(U);p=f.eos(rho,lt)[0];mask=(rho>=f.eos.floor)&(rho<f.eos.original.floor)
            self.tail.append(dict(step=len(self.times),cells=int(mask.sum()),
                mass_g=float(np.sum(U[0,mask]*m.vol[mask])*self.gas_scale/C**2),
                Killing_energy_erg=float(np.sum((U[2,mask]+m.a0*f.eos.cx*U[0,mask])*m.vol[mask])*self.gas_scale),
                pressure_volume_erg=float(np.sum(p[mask]*m.vol[mask])*self.gas_scale),
                maximum_speed_c=float(np.max(abs(v[mask]),initial=0)),
                minimum_temperature_K=float(np.min(np.exp(lt[mask]),initial=np.inf)) if mask.any() else None))
        finally:f.seed=seed;f.eos.y=oldy
        write(EV/f'tail-{self.steps}.json',dict(classification='Counterexample candidate',rows=self.tail))
        return value

    def radiate(self,*args):
        start=time.monotonic();result=super(previous.Capture,self).radiate(*args);self.cost['physics']+=time.monotonic()-start
        start=time.monotonic();I=result[1].sum(0)
        luminosity=2*np.pi*C*self.area[-1]*(I[-1,self.mu>0]@(self.number*self.E))
        self.angular.append(luminosity)
        atomic_npz(self.port_file,angular_luminosity=self.angular,h=flow.old.END/self.steps)
        self.cost['port_write']+=time.monotonic()-start;self.cost['port_count']+=1
        return result

    def export(self):
        ks=range(0,self.steps+1,self.steps//16);states=[dict(load_npz(EV/f'state-{self.steps}-{k}.npz')) for k in ks]
        atomic_npz(EV/f'coupled-{self.steps}.npz',snapshot_t=np.linspace(0,flow.old.END,17),
            **{'snapshot_'+key:np.array([r[key] for r in states]) for key in states[0]})
        current=load_npz(EV/f'source-{self.steps}.npz');mu=self.mu[self.mu>0];w=self.w[self.mu>0]
        energy=flow.old.END/self.steps*np.cumsum(np.asarray(self.angular)@(2*w*mu))
        reference=current['outer_cumulative_energy_erg'][1:]
        err=float(np.max(abs(energy-reference))/max(abs(reference[-1]),1.));assert err<1e-12
        d=dict(current);selected=np.arange(0,self.steps+1,self.steps//16)
        for key in d:
            if key=='t' or key in physical.capture.KEYS+['metric_stress_erg','inner_cumulative_energy_erg','outer_cumulative_energy_erg']:d[key]=d[key][selected]
        atomic_npz(GR/f'source-{self.steps}.npz',**d)
        write(EV/f'capture-{self.steps}.json',dict(classification='Counterexample candidate',passed=True,
            canonical_states=17,accepted_angular_port_relative=err,changed_physics=True,identical_old_path_replay=False))


def prepare():
    assert read(BANK/'result.json')['passed'] and not (OUT/'coupling-plan.json').exists()
    paths=[Path(__file__),BANK/'repaired-bank.npz',BANK/'spectrum-bank.npz',OUT/'reallocated-bank-plan.json',
        Path(previous.__file__),Path(physical.__file__),Path(cold.__file__)]
    write(OUT/'coupling-plan.json',dict(classification='Counterexample candidate',
        claim='Apply the retained native tail to the unchanged actual paired photon/material SSP evolution; record nonzero retained pressure/motion, remaining deletion, and complete source histories.',
        clocks=[64,128],horizon_seconds=flow.old.END,controls_seconds=60,controls_native_calls=220,
        pilot_seconds=70,production_seconds=650,GR_seconds=10,
        controls='Withheld native constitutive and integrated spectral values, two sampled native derivative comparisons, old-support bitwise reuse and the seam. No uniform EOS or derivative certificate is inferred.',
        controls_states=[[-23.7,90.,2e-14],[-22.7,310.,1e-10],[-21.,1900.,1e-8],[-19.2,11200.,.0003],[-18.2,18200.,.0007]],
        gates=read(OUT/'plan.json')['gates'],
        stop='Stop on unsupported actual conserved state; retain its U and root failure. No cold clamp or further floor reduction. Preserve old trajectory and all failed precursor probes.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths}))


def controls():
    assert not (OUT/'coupling-controls.json').exists();start=time.monotonic();support.owner.deadline(start,60)
    plan=read(OUT/'coupling-plan.json');n=support.fast_native(220)
    original=cold.ColdEOS();oldspec=cold.ColdSpectrum();e=EOS(original);s=Spectrum(oldspec)
    d=np.load(cold.prior.optical.photons.OUT/'bank-16-8.npz');a=d['face_a'][-1];E=d['Einf']/a;num=d['num']/a**3
    rows=[];failure=None
    try:
        for x,T,y in plan['controls_states']:
            t=np.log(T);z=n.state(x,t,y);e.y=np.array([y]);p,u,g,_,_,cv,_=e.evaluate(np.array([np.exp(x)]),np.array([t]))
            ce=max(abs(p[0]*e.rho0*C**2/z['raw'][1]-1),abs(u[0]*C**2/z['raw'][2]-1))
            ab=np.zeros_like(E);em=np.zeros_like(E)
            for j in range(10):
                sigma=s.cross(E,j+1);ab+=np.exp(np.log(y)+z['log_fraction'][j])*sigma;good=sigma>0
                em[good]+=np.exp(np.log(y)+z['log_fraction'][j]+z['affinity']-E[good]/(cold.prior.optical.ex.K*T))*sigma[good]
            aa,ee=s.coefficients(np.array([e.rho0*np.exp(x)]),np.array([t]),np.array([y]),E[None,:]);errors=[]
            for Trad in [13400,20000,80000]:
                occ=1/np.expm1(E/(cold.prior.optical.ex.K*Trad))
                for true,est,field in [(ab,aa[0],occ),(em,ee[0],occ),(em,ee[0],np.ones_like(E))]:
                    for factor in [num,num*E]:
                        weight=factor*field;errors.append(float(np.sum(abs(est-true)*weight)/max(float(true@weight),1e-280)))
            derivative=None
            if x in [-22.7,-19.2]:
                h=1e-4;ar=(n.state(x+h,t,y)['raw']-n.state(x-h,t,y)['raw'])/(2*h)
                at=(n.state(x,t+h,y)['raw']-n.state(x,t-h,y)['raw'])/(2*h);raw=z['raw']
                gamma=ar[1]/raw[1]+at[1]/raw[1]*(raw[1]/raw[0]-ar[2])/at[2]
                derivative=float(max(abs(g[0]/gamma-1),abs(cv[0]*C**2/at[2]-1)))
            rows.append(dict(x=x,T=T,y=y,constitutive=float(ce),spectrum=max(errors),derivative=derivative))
        rho=np.exp(np.array([-17.,-6.,-.1]));lt=np.log([180.,1900.,13400.]);yy=np.array([1e-10,1e-8,.0003]);original.y=yy;e.y=yy
        assert all(np.array_equal(a,b) for a,b in zip(e.evaluate(rho,lt),original.evaluate(rho,lt)))
        assert all(np.array_equal(a,b) for a,b in zip(s.levels(rho*e.rho0,lt,yy),oldspec.levels(rho*e.rho0,lt,yy)))
        rr=np.exp(np.array([-18.-1e-8,-18.]));seam=[]
        for T in [90.,310.,1900.,18200.]:
            e.y=np.full(2,1e-7);v=e.evaluate(rr,np.full(2,np.log(T)))
            seam.append(float(max(abs(v[0][0]/v[0][1]-1),abs(v[1][0]/v[1][1]-1),abs(v[2][0]/v[2][1]-1),abs(v[5][0]/v[5][1]-1))))
        assert max(max(r['constitutive'],r['spectrum'],r['derivative'] or 0) for r in rows)<.002
        assert max(seam)<.002
        support.owner.deadline(start,60)
    except Exception as exc:failure=repr(exc)
    result=dict(classification='Counterexample candidate',passed=failure is None,failure=failure,checks=rows,
        native_calls=n.ion.calls+n.variant_initial_calls,seconds=time.monotonic()-start,
        old_support_bitwise=failure is None,seam_relative=seam if 'seam' in locals() else None,uniform_derivative_certificate=False)
    write(OUT/'coupling-controls.json',result);signal.setitimer(signal.ITIMER_REAL,0.);print(json.dumps(result),flush=True)
    assert result['passed']


def pilot():
    assert read(OUT/'coupling-controls.json')['passed'] and not (OUT/'evolution-pilot.json').exists()
    start=time.monotonic();support.owner.deadline(start,70);rows=[]
    try:
        for n in [64,128]:
            row=Capture().run(n,4);row['new_steps']=4;rows.append(row);assert row['passed']
        state=max(r['cost']['state_write']/max(r['cost']['state_count'],1) for r in rows);parts=[]
        for r in rows:
            c=r['cost'];remaining=r['steps']-r['completed'];fixed=max(r['seconds']-c['physics']-c['state_write']-c['port_write'],0.)
            parts.append(dict(setup=c['setup'],steps=(c['physics']+c['port_write'])/4*remaining,
                states=state*int(np.ceil(remaining/(r['steps']//16))),checkpoints=fixed*(int(np.ceil(remaining/16))+1)))
        estimate=sum(sum(p.values()) for p in parts)+30;upper=1.5*max(estimate,353.30)+20;support.owner.deadline(start,70)
        aggregate=read(OUT/'aggregate-budget-plan.json')
        probe_spent=sum(read(p)['seconds'] for p in support.FAILURES)+read(OUT/'warm-probe.json')['seconds']
        spent=probe_spent+aggregate['bank_already_spent_seconds']+read(BANK/'result.json')['seconds']+read(OUT/'coupling-controls.json')['seconds']+(time.monotonic()-start)
        cap=min(650,1095-spent-10)
        result=dict(classification='Counterexample candidate',paths=rows,cost_parts=parts,forecast_seconds=estimate,upper_seconds=upper,
            eligible=upper<cap,seconds=time.monotonic()-start,aggregate_spent_seconds=spent,production_cap_seconds=cap)
        write(OUT/'evolution-pilot.json',result);print(json.dumps(result),flush=True)
        if result['eligible']:write(OUT/'execution-plan.json',dict(eligible=True,hard_seconds=cap,source_sha256=sha(__file__),upper_seconds=upper))
    except Exception as exc:write(OUT/'pilot-failure.json',dict(classification='Counterexample candidate',error=repr(exc),paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def production():
    assert not (EV/'result.json').exists();plan=read(OUT/'execution-plan.json');assert plan['eligible'] and sha(__file__)==plan['source_sha256']
    start=time.monotonic();cap=plan['hard_seconds'];support.owner.deadline(start,cap);rows=[]
    try:
        for n in [64,128]:
            row=Capture().run(n,resume=True);rows.append(row);assert row['passed'];support.owner.deadline(start,cap)
        result=dict(classification='Counterexample candidate',passed=True,paths=rows,seconds=time.monotonic()-start,
            actual_native_low_density_material_applied=True,nonlinear_GR_closed=False,final_charge_solved=False)
        write(EV/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:write(EV/'failure.json',dict(classification='Counterexample candidate',error=repr(exc),paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def resume_pilot():
    assert not (OUT/'resumed-pilot-plan.json').exists()
    write(OUT/'resumed-pilot-plan.json',dict(classification='Counterexample candidate',
        prior='The70s externally bounded pilot exited with process code1 after accepted64/4 and an incomplete128 prefix. No final exception receipt was written, so the exact exit mechanism is not established. Charge the full70s; preserve its original producer, complete64 checkpoint and orphaned fine angular port.',
        decision='Continue64 from accepted step4 for one measured step and run only two fine steps. Persist each path cost immediately. Both are prefixes of the original approved64/128 trajectories. No accepted stored step is replayed.',
        hard_seconds=60,budget='Use unspent native-bank/control allocation under the same aggregate1095s cap. Original failed pilot remains failed. Production must fit the newly calculated aggregate remainder.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'executed-evolution-pilot.py',EV/'checkpoint-64.npz',EV/'orphaned-pilot-ports-128.npz']}))
    start=time.monotonic();support.owner.deadline(start,60);rows=[]
    try:
        for n,count,resume in [(64,5,True),(128,2,False)]:
            model=Capture();row=model.run(n,count,resume);row['new_steps']=1 if resume else 2;rows.append(row);assert row['passed']
            del model
            import gc
            gc.collect()
        state=max(r['cost']['state_write']/max(r['cost']['state_count'],1) for r in rows);parts=[]
        for r in rows:
            c=r['cost'];remaining=r['steps']-r['completed'];fixed=max(r['seconds']-c['physics']-c['state_write']-c['port_write'],0.)
            parts.append(dict(setup=c['setup'],steps=(c['physics']+c['port_write'])/r['new_steps']*remaining,
                states=state*int(np.ceil(remaining/(r['steps']//16))),checkpoints=fixed*(int(np.ceil(remaining/16))+1)))
        estimate=sum(sum(p.values()) for p in parts)+30;upper=1.5*max(estimate,353.30)+20
        aggregate=read(OUT/'aggregate-budget-plan.json');probe_spent=sum(read(p)['seconds'] for p in support.FAILURES)+read(OUT/'warm-probe.json')['seconds']
        spent=probe_spent+aggregate['bank_already_spent_seconds']+read(BANK/'result.json')['seconds']+read(OUT/'coupling-controls.json')['seconds']+70+(time.monotonic()-start)
        cap=min(650,1095-spent-10);support.owner.deadline(start,60)
        result=dict(classification='Counterexample candidate',paths=rows,cost_parts=parts,forecast_seconds=estimate,upper_seconds=upper,
            eligible=upper<cap,seconds=time.monotonic()-start,aggregate_spent_seconds=spent,production_cap_seconds=cap,original_pilot_passed=False)
        write(OUT/'resumed-pilot.json',result);print(json.dumps(result),flush=True)
        if result['eligible']:write(OUT/'execution-plan.json',dict(eligible=True,hard_seconds=cap,source_sha256=sha(__file__),upper_seconds=upper))
    except Exception as exc:write(OUT/'resumed-pilot-failure.json',dict(error=repr(exc),paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def production_plan():
    assert not (OUT/'execution-plan.json').exists();p=read(OUT/'resumed-pilot.json');assert not p['eligible']
    paths=[Path(__file__),OUT/'executed-resumed-pilot.py',OUT/'resumed-pilot.json']
    paths += [EV/f'preproduction-{kind}-{n}.npz' for n in [64,128] for kind in ['checkpoint','history','source']]
    write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=p['upper_seconds']<900,
        hard_seconds=900,source_sha256=sha(__file__),forecast_seconds=p['forecast_seconds'],upper_seconds=p['upper_seconds'],
        original650s_forecast_passed=False,original_aggregate_budget_seconds=1095,
        revised_total_scientific_action_cap_seconds=1300,already_spent_scientific_seconds=p['aggregate_spent_seconds'],
        reassessment='The accepted native constitutive controls and real retained material now permit the decisive full-horizon comparison. Measured original-clock prefixes forecast572.97s and the unchanged1.5x margin gives879.45s. Authorize one900s continuation of the original64/128 paths, from saved5/2 prefixes. This explicitly raises the prior allocation; rejected dispatch forecasts remain false. No extra physical path, resolution, support or accuracy relaxation.',
        reason='An early prefix establishes material retention but cannot determine its full causal source/charge effect. Previous failed setup attempts are not counted as full evolution. The cheaper saved tables, sparse native controls, constructor/import cache and successful prefixes are all reused.',
        stop='One measured continuation only. Stop on actual support, conservation, deadline or later acceptance failure. No automatic second900s run or larger numerical path.',
        limitations='Later SSP/native-root costs remain unmeasured; forecast is not a completion guarantee. Full EOS/derivative, nonlinear GR, same-inventory static comparison and observational closure remain open.',
        bindings={(v.relative_to(Path.cwd()) if v.is_absolute() else v).as_posix():sha(v) for v in paths}))


def checkpoint_control():
    import tempfile,zipfile
    with tempfile.TemporaryDirectory() as folder:
        p=Path(folder)/'checkpoint.npz';values=dict(completed=16,U=np.arange(128,dtype=np.longdouble).reshape(4,32))
        atomic_npz(p,**values);before=p.read_bytes();p.with_suffix('.saving').write_bytes(b'interrupted')
        assert p.read_bytes()==before
        with load_npz(p) as z:assert int(z['completed'])==16 and np.array_equal(z['U'],values['U'])
        history=np.arange(33);assert np.array_equal(history[:int(values['completed'])+1],np.arange(17))
        atomic_npz(p,completed=32,U=values['U']+1)
        with zipfile.ZipFile(p) as z:assert z.testzip() is None
        with load_npz(p) as z:assert int(z['completed'])==32 and np.array_equal(z['U'],values['U']+1)
    assert "list(saved['moments'][:begin+1])" in runner and 'atomic_npz(path,**values)' in runner
    result=dict(passed=True,original_checkpoint_survives_interrupted_write=True,complete_replacement=True,
        history_truncated_to_checkpoint=True,scientific_re_evolution=False,source_sha256=sha(__file__))
    write(OUT/'checkpoint-control.json',result);print(json.dumps(result))


if __name__=='__main__':
    import resource
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    globals()[sys.argv[1]]()
