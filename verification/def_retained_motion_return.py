"""Counterexample candidate: return actual native-driven material to photons.

One finite deferred correction. Frozen native constitutive defects stay fixed;
neither a native Jacobian nor a coupled fixed-point bound is inferred.
"""
from pathlib import Path
from types import FunctionType
import inspect, json, os, resource, signal, sys, textwrap, time
import def_retained_native_return as prior

ROOT=prior.ROOT
BEFORE=Path(os.environ.get('RETAINED_PREVIOUS_OUTPUT','retained-native-return150-work'))
OUT=Path(os.environ.get('RETAINED_MOTION_OUTPUT','retained-motion-return151-work'))
SOURCE=OUT/'collision'; PHOTON=OUT/'increment'; TOTAL=OUT/'photons'; MATERIAL=OUT/'material'; GR=OUT/'gr'
write,read,sha=prior.write,prior.read,prior.sha
CAPS=dict(source_pilot=120,source_production=480,photon_pilot=90,photon_production=1200,
          material_pilot=30,material_production=180,readout=180,audit=60)


def prepare():
    assert not OUT.exists()
    for p in [OUT,SOURCE,PHOTON,TOTAL,MATERIAL,GR]:p.mkdir(parents=True)
    paths=[Path(__file__),Path(prior.__file__),Path('verification/repair_native_finite_collision_remainder.py'),
           prior.EV/'coupled-128.npz',prior.EV/'source-128.npz',BEFORE/'infinity/result.json',BEFORE/'result.json']
    paths+=list((BEFORE/'photons/bank-128').glob('*.npz'))
    for folder in ['photons','material']:
        paths += [BEFORE/folder/f'steps-{n}-reference-128.npz' for n in [64,128]]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='93150dc9d',
        claim='Apply the saved actual finite native-driven material motion to photon collisions, return the paired increment to free material, and compare the represented GR/infinity charge.',
        decision='Resolve whether the missing reciprocal material-to-photon step screens the surviving nominal charge; quantify the remaining state mismatch before claiming a fixed point.',
        equation='R=Ctable(Qb+Qm,Ib+x0)-Ctable(Qb,Ib)-L0(x0,g0); delta_xdot=L0_photon(delta_x,delta_g)+R. The new material radiation baseline is gm, not g0. Its known mechanical transport cancels between the equations; do not add it or the fixed native defect twice.',
        scope='Same531 cells,17 saved knots,8 angles,152 frequencies,64/128 matched response clocks,3.4344311179287023ms. Evaluate finite collisions at actual saved material for EACH clock. Old tangent is a correction operator, not a native derivative certificate.',
        reuse='Reuse completed background, native force/collision defects and coefficient banks. Propagate only the new photon correction; recompute only its required free-material return. No microscopic state calls, bank rebuild, extra clocks or longer horizon.',
        method='Extended conserved recovery and spline evaluation including low-density branches, coefficient-first subtraction, exact scattering-number projection at zero energy/momentum. Half amplitude is a diagnostic at0/8/16, never a replacement result.',
        gates=dict(base_owner=1e-9,photon_negativity=1e-12,arithmetic_over_linear=.002,number=1e-10,
                   projection_over_linear=1e-8,linear=1e-12,conservation=1e-8,time=.02,quadrature=.002,independent_GR=1e-9),
        budgets=CAPS,total_max_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,
        forecast='Reuse source pilot points; twice the slowest representative point for the remaining count plus setup must fit480s. Photon production uses twice the greater measured prefix or previous late step cost. Material uses previous actual raw-call count and fresh two-step costs. Stop if any forecast exceeds its registered cap.',
        stop='Preserve and resolve the first physical/arithmetic/support/time/budget failure. No automatic finer clocks, relaxed gates, lower amplitude or larger aggregate budget.',
        limitations='One finite constitutive correction and free-material sweep; native defect derivatives, new-state remainder, metric recurrence, spatial/floor errors and final physical charge remain open.',
        bindings={str(p):sha(p) for p in paths}))
    import sympy as sp
    L0,R,D,H=sp.symbols('L0 R D H'); gm_dot=L0+D+H
    assert sp.expand((L0+R+D+H)-gm_dot-R)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='Subtract the baseline gm equation: fixed native D and known mechanical H cancel, leaving only finite collision residual R. This proves no nonlinear convergence estimate.'))


def initialize():
    global np,LD,AMP,C,State,Response,Material
    import numpy as np
    import repair_native_finite_collision_remainder as extended
    LD=np.longdouble;AMP=1e-26
    prior.OUT=OUT;prior.PHOTON=BEFORE/'photons';prior.MATERIAL=BEFORE/'material'
    prior.HYDRO=BEFORE/'hydro';prior.FORCE=BEFORE/'collision';prior.METRIC=BEFORE/'metric'
    # Reuse the exact immutable cache; prior.initialize verifies its SHA.
    alias=OUT/'immutable-coupled-128.npz'
    if not alias.exists():alias.symlink_to((BEFORE/'immutable-coupled-128.npz').resolve())
    prior.initialize();C=prior.C

    class StateImpl:
        def __init__(self):
            self.m=m=prior.Response();m.material=prior.Material(driven=False);m.model=m.material.model
            self.owners=[];self.splines=[];seen=set();spline_cache={}
            # The retained tail adds low/original owners to the existing cold/warm pair.
            todo=[m.model.flow.eos,m.model.spectrum]
            while todo:
                obj=todo.pop()
                if id(obj) in seen:continue
                seen.add(id(obj))
                todo += [getattr(obj,k) for k in ['low','original','cold','warm'] if hasattr(obj,k)]
                for key in ['f','rev']:
                    if not hasattr(obj,key):continue
                    groups=getattr(obj,key)
                    if not isinstance(groups,list):continue
                    after=[]
                    for group in groups:
                        row=[]
                        for s in group:
                            if id(s) not in spline_cache:
                                spline_cache[id(s)]=extended.Spline(s);self.splines.append(spline_cache[id(s)])
                            row.append(spline_cache[id(s)])
                        after.append(row)
                    self.owners.append((obj,key,groups,after))
                for key,anchor in [('evaluate','np.zeros((7,len(rho)))'),('levels','np.zeros((len(rho),10))')]:
                    if hasattr(obj,key):
                        fn=getattr(obj,key)
                        if anchor in inspect.getsource(fn):
                            self.owners.append((obj,key,fn,extended.promote(fn,[(anchor,anchor[:-1]+',dtype=LD)')])))
            s=m.model.spectrum;cross=s.cross
            self.owners.append((s,'cross',cross,lambda energy,level:cross(np.asarray(energy,float),level).astype(LD)))
            for key,edits in [('coefficients',[('np.zeros((n,self.q,self.nf))','np.zeros((n,self.q,self.nf),dtype=LD)'),('np.zeros(n)','np.zeros(n,dtype=LD)')]),
                              ('scattering_matrix',[('np.zeros((3,n,q,nf))',"np.zeros((3,n,q,nf),dtype=c['beta'].dtype)")])]:
                fn=getattr(m,key);self.owners.append((m,key,fn,extended.promote(fn,edits)))
        precision=extended.State.precision
        coefficients=extended.State.coefficients
        def select(self,n):
            self.n=n;self.folder=SOURCE/str(n);self.folder.mkdir(exist_ok=True)
            self.p=dict(np.load(BEFORE/f'photons/steps-{n}-reference-128.npz'))
            z=np.load(BEFORE/f'material/steps-{n}-reference-128.npz');ids=[np.argmin(abs(z['t']-t)) for t in self.m.t]
            assert np.max(abs(z['t'][ids]-self.m.t))<1e-18
            self.m.motion=z['history_scaled'][ids]
        def point(self,k,half=False):
            m=self.m;p=self.p;mark=time.monotonic();t=m.t[k];assert p['t'][k]==t
            self.precision(False);m.model.flow.seed=np.asarray(m.model.flow.seed,float)
            target=m.point(k);point=m.point;m.point=lambda index:target
            try:c=m.local(t)
            finally:m.point=point
            x=p['photon_history_scaled_occupation'][k]/(m.scale*AMP);g=p['material_history'][k]/AMP
            linear,_,escape,bound=m.collision(c,x,g,False)
            linear=linear*m.scale*AMP;bound=bound*m.scale*AMP;escape=escape*AMP
            # Full actual free-material state, including its noncollisional E/H.
            z=m.motion[k];self.precision(True);a=self.coefficients(k,np.zeros_like(z),0.)
            sa,ea=m.scattering_matrix(a);bank=np.load(BEFORE/f'photons/bank-128/point-{k}.npz')
            owner=max(float(np.max(abs(a[q]-bank[q]))/max(np.max(abs(bank[q])),1e-300)) for q in ['emit','loss','sc','rho','beta','u','p'])
            I=m.I[k].astype(LD);delta=p['photon_history_scaled_occupation'][k].astype(LD);X=I/m.scale;dx=delta/m.scale
            weight=m.Eweight/m.scale;N=m.Nweight/m.scale;rows=[];full=None
            for factor in ([1.,.5] if half else [1.]):
                b=self.coefficients(k,z,factor);sb,eb=m.scattering_matrix(b)
                db=(b['emit']-a['emit'])-(b['loss']-a['loss'])*I-b['loss']*delta*factor
                scatter=((sb-sa)@X.ravel()).reshape(X.shape)*m.scale+(sb@(dx*factor).ravel()).reshape(X.shape)*m.scale
                de=np.einsum('knqf,nqf->kn',eb-ea,X)+np.einsum('knqf,nqf->kn',eb,dx*factor)
                residual=db+scatter-factor*linear;br=db-factor*bound;er=de-factor*escape
                norm=max(np.sum(abs(factor*linear)*weight),1.)
                rounding=16*np.finfo(LD).eps*np.sum((abs(a['emit'])+abs(b['emit'])+(abs(a['loss'])+abs(b['loss']))*abs(I))*weight)
                if not np.any(z) and not np.any(delta):rounding=0.
                original=residual.copy();elo,ehi=LD(m.E[0]),LD(m.E[-1])
                for _ in range(2):
                    number=np.sum((residual-br)*N,axis=(1,2))+er[0]
                    residual[:,0,0]-=number*ehi/(ehi-elo)/N[:,0,0]
                    residual[:,0,-1]+=number*elo/(ehi-elo)/N[:,0,-1]
                projection=float(np.sum(abs(residual-original)*weight)/norm)
                number=np.sum((residual-br)*N,axis=(1,2))+er[0]
                ns=np.maximum(np.sum((abs(residual)+abs(br))*N,axis=(1,2))+abs(er[0]),1.)
                negative=float(np.sum(np.maximum(-(I+factor*delta),0)*weight)/max(np.sum(abs(I)*weight),1.))
                rows.append(dict(factor=factor,residual_over_linear=float(np.sum(abs(residual)*weight)/norm),
                    rounding_over_linear=float(rounding/norm),number_relative=float(np.max(abs(number)/ns)),
                    projection_over_linear=projection,negative_photon_energy_relative=negative))
                if factor==1.:full=dict(t=t,photon=residual,bound=br,escape=er)
            row=dict(classification='Counterexample candidate',steps=self.n,k=k,t=float(t),base_owner_relative=owner,rows=rows,
                seconds=time.monotonic()-mark,passed=bool(owner<1e-9 and all(r['rounding_over_linear']<.002 and r['number_relative']<1e-10 and r['projection_over_linear']<1e-8 and r['negative_photon_energy_relative']<1e-12 for r in rows)))
            np.savez_compressed(self.folder/f'point-{k}.npz',**full);write(self.folder/f'point-{k}.json',row)
            self.precision(False);print(json.dumps(row),flush=True);return row
    State=StateImpl

    class ResponseImpl(prior.Response):
        def __init__(self,reference=128,steps=128):super().__init__(reference);self.steps=steps
        def defect(self,k):
            if k not in self.defects:
                d=np.load(SOURCE/str(self.steps)/f'point-{k}.npz')
                self.defects[k]={key:np.asarray(d[key],float) for key in ['photon','bound','escape']}
                for j in list(self.defects):
                    if j not in [k-1,k,k+1]:del self.defects[j]
            return self.defects[k]
        run=FunctionType(prior.Response.run.__code__,dict(prior.Response.run.__globals__,OUT=PHOTON),argdefs=prior.Response.run.__defaults__)
    Response=ResponseImpl
    class MaterialImpl(prior.Material):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps,driven=False);p=np.load(TOTAL/f'steps-{steps}-reference-128.npz');ids=[np.argmin(abs(p['t']-t)) for t in self.t]
            assert np.max(abs(p['t'][ids]-self.t))<1e-18;c=p['collision_transfer'][ids]
            self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/AMP
        run=FunctionType(prior.Material.run.__code__,dict(prior.Material.run.__globals__,OUT=MATERIAL),argdefs=prior.Material.run.__defaults__)
    Material=MaterialImpl


def source(pilot):
    if not pilot:assert read(SOURCE/'pilot.json')['eligible']
    mark=time.monotonic();s=State();rows=[]
    for n in [64,128]:
        s.select(n)
        for k in ([0,8,16] if pilot else range(17)):
            if (s.folder/f'point-{k}.json').exists():row=read(s.folder/f'point-{k}.json')
            else:row=s.point(k,half=pilot)
            rows.append(row)
            assert row['passed'],row
    if pilot:
        forecast=2*max(r['seconds'] for r in rows)*28+15
        result=dict(rows=rows,upper_remaining_seconds=forecast,eligible=forecast<CAPS['source_production'])
    else:result=dict(rows=rows,passed=True)
    write(SOURCE/('pilot.json' if pilot else 'result.json'),dict(classification='Counterexample candidate',seconds=time.monotonic()-mark,**result))
    if pilot:assert result['eligible'],result['upper_remaining_seconds']


def photon(pilot):
    assert read(SOURCE/'result.json')['passed'];mark=time.monotonic();rows=[]
    if not pilot:assert read(PHOTON/'pilot.json')['eligible']
    for n in [64,128]:
        m=Response(steps=n);label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        rows.append(m.run(n,label,n//16 if pilot else None,restart=None if pilot else f'pilot-{n}'))
    keys=[0,1,2,3,5,6];a=np.load(PHOTON/('pilot-64.npz' if pilot else 'steps-64-reference-128.npz'))['moments'][:,keys]
    b=np.load(PHOTON/('pilot-128.npz' if pilot else 'steps-128-reference-128.npz'))['moments'][:,keys]
    errors=(np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1.)).tolist()
    result=dict(classification='Counterexample candidate',rows=rows,time_comparison=errors,passed=all(r['passed'] for r in rows) and max(errors)<.02,seconds=time.monotonic()-mark)
    if pilot:
        previous=read(BEFORE/'photons/result.json')['paths'];forecast=[]
        for r in rows:
            old=next(v for v in previous if v['steps']==r['steps'])
            step=max(r['stepping_seconds']/r['new_steps'],old['stepping_seconds']/old['new_steps'])
            forecast.append(17*r['operator_point_seconds']/r['operator_points']+(r['steps']-r['completed_steps'])*step+20)
        result.update(upper_remaining_seconds=2*sum(forecast),eligible=result['passed'] and 2*sum(forecast)<CAPS['photon_production'])
    write(PHOTON/('pilot.json' if pilot else 'result.json'),result);print(json.dumps(result),flush=True)
    assert result['eligible' if pilot else 'passed']


def combine():
    assert read(PHOTON/'result.json')['passed'];m=Response();rows=[]
    for n in [64,128]:
        a=dict(np.load(BEFORE/f'photons/steps-{n}-reference-128.npz'));b=dict(np.load(PHOTON/f'steps-{n}-reference-128.npz'))
        assert np.array_equal(a['t'],b['t']) and np.array_equal(a['accepted_angular_times'],b['accepted_angular_times'])
        z=np.load(BEFORE/f'material/steps-{n}-reference-128.npz');ids=[np.argmin(abs(z['t']-t)) for t in a['t']];motion=z['history_scaled'][ids]*AMP
        total=dict(a)
        for key in ['moments','delta_packet_scaled_occupation','ledger','escape','collision_transfer','radial_ports','photon_history_scaled_occupation','accepted_angular_luminosity']:
            total[key]=a[key]+b[key]
        total['material_history']=np.moveaxis(motion[:,[2,3]],1,2)/np.stack([m.eu,m.nu],axis=-1)+b['material_history']
        total['delta_material']=total['material_history'][-1]
        total['moments'][:,[1,2]]=motion[:,[2,3]]+b['moments'][:,[1,2]]
        total['moments'][:,4]=np.sum(abs(total['photon_history_scaled_occupation'])*(m.Eweight/m.scale),axis=(2,3))
        # ponytail: moment6 is diagnostic tangent pressure; GR uses the actual free-material finite pressure below.
        np.savez_compressed(TOTAL/f'steps-{n}-reference-128.npz',**total)
        rows.append(dict(steps=n,actual_material_baseline=True,collision_transfer_separate_from_mechanical_energy=True))
    write(TOTAL/'result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        old_photons_replayed=False,new_state_remainder_certified=False))


def material(pilot):
    assert read(TOTAL/'result.json')['passed'];mark=time.monotonic()
    if not pilot:assert read(MATERIAL/'pilot.json')['eligible']
    fn=prior.matter.old.worker
    worker=FunctionType(fn.__code__,dict(fn.__globals__,OUT=MATERIAL,Material=Material,configure=lambda:None),argdefs=fn.__defaults__)
    rows=[worker(n,f'pilot-{n}' if pilot else f'steps-{n}-reference-128',2 if pilot else None,None if pilot else f'pilot-{n}') for n in [64,128]]
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-mark)
    if pilot:
        previous=[read(BEFORE/f'material/steps-{n}-reference-128.json') for n in [64,128]]
        forecast=[old['raw_owner_calls']*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+5 for r,old in zip(rows,previous)]
        result.update(upper_remaining_seconds=2*sum(forecast),eligible=2*sum(forecast)<CAPS['material_production'])
    write(MATERIAL/('pilot.json' if pilot else 'production.json'),result);print(json.dumps(result),flush=True)
    if pilot:assert result['eligible']


def readout():
    ns=dict(vars(prior),OUT=OUT,PHOTON=TOTAL,MATERIAL=MATERIAL,GR=GR,Material=Material)
    FunctionType(prior.readout.__code__,ns)()
    compact=read(OUT/'result.json');previous=read(BEFORE/'result.json')
    compact.update(actual_material_motion_returned_to_photons=True,
        motion_return_compact_endpoint=compact['native_return_compact_endpoint']-previous['native_return_compact_endpoint'])
    write(OUT/'result.json',compact)
    ns['CAPS']=dict(prior.CAPS,infinity=120)
    FunctionType(prior.infinity.__code__,ns)()
    r=read(OUT/'infinity/result.json');old=read(BEFORE/'infinity/result.json')
    # Subtract the already separated constitutive corrections, not the much
    # larger baseline charge. Both use the same preserved Phase149 baseline.
    fine=np.load(OUT/'infinity/charge-128-a8-r8.npz')['native_return'];before=np.load(BEFORE/'infinity/charge-128-a8-r8.npz')['native_return'];inc=fine-before
    changes={}
    for key,n,a,b in [('time',64,8,8),('angular',128,4,8),('radial',128,8,4)]:
        now=np.load(OUT/f'infinity/charge-{n}-a{a}-r{b}.npz')['native_return'];was=np.load(BEFORE/f'infinity/charge-{n}-a{a}-r{b}.npz')['native_return']
        changes[key]=float(np.max(abs(now-was-inc))/max(np.max(abs(inc)),1e-300))
    r.update(actual_material_motion_returned_to_photons=True,previous_phase_endpoint=old['endpoint_normalized'],
        motion_return_endpoint=float(inc[-1]),motion_return_controls=changes,
        coupled_fixed_point_verified=False,final_charge_solved=False,full_goal_complete=False)
    r['passed']=r['passed'] and changes['time']<.02 and max(changes['angular'],changes['radial'])<.002
    write(OUT/'infinity/result.json',r);assert r['passed'],changes


def execute(action):
    plan=read(OUT/'plan.json')
    revision=read(OUT/'execution-plan.json') if (OUT/'execution-plan.json').exists() else None
    for p,h in plan['bindings'].items():
        target=OUT/'registered-producer.py' if revision and Path(p).name==Path(__file__).name else p
        assert sha(target)==h,p
    if revision:assert sha(__file__)==revision['source_sha256']
    initialize()
    if action.startswith('source_'):source(action.endswith('pilot'))
    elif action.startswith('photon_'):photon(action.endswith('pilot'))
    elif action=='combine':combine()
    elif action.startswith('material_'):material(action.endswith('pilot'))
    elif action=='readout':readout()
    else:raise ValueError(action)


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    else:
        assert not (OUT/f'{action}-receipt.json').exists(),'Preserve completed attempts; do not overwrite receipts'
        spent=sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))
        remaining=sum(CAPS.values())-spent;assert remaining>0
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        def timeout(*_):raise TimeoutError('Registered stage budget')
        signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,min(CAPS.get(action,60),remaining))
        start=time.monotonic();cpu=time.process_time();error=None
        try:execute(action)
        except Exception as exc:error=repr(exc);raise
        finally:
            signal.alarm(0);write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=sha(__file__)))
