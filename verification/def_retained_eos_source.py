"""Re-use exact native roots and apply changed retained-history EOS sources."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
from functools import lru_cache
import hashlib,io,json,signal,sys,time,subprocess

ROOT=Path('outputs/direct-eos-gr33');CURRENT=ROOT/'native-retained-completion';EV=CURRENT/'evolution'
THERMAL=ROOT/'def-native-conservative-rates/thermal-refined';PRIOR=THERMAL/'updated-gr-return/conserved-history-charge'
OUT=ROOT/'retained-eos-source';read=lambda p:json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path('verification/evolve_native_retained_tail.py'),Path('verification/def_native_conserved_history_charge.py'),
        Path('verification/verify_native_conserved_history_charge.py'),Path('verification/def_native_atmosphere_inverse.py'),
        Path('verification/def_native_conservative_eos_readout.py'),PRIOR/'production-samples.json',PRIOR/'production.json',
        PRIOR/'audit.json',EV/'coupled-128.npz',EV/'source-128.npz',EV/'result.json',CURRENT/'gr/result.json',
        CURRENT/'infinity/completed/result.json',CURRENT/'infinity/completed/retained-128-a8-r8.npz',
        THERMAL/'updated-gr-return/background/coupled-128.npz',THERMAL/'bank.npz',
        CURRENT/'gr/wave-128-g8.npz',PRIOR/'plan.json',
        Path('verification/native_tail_supported_temperature.py'),Path('verification/def_native_cold_coupling.py'),
        ROOT/'def-native-cold-coupling/logarithmic-build.json',
        ROOT/'native-retained-tail/supported-temperature/support/repaired-bank.npz']
    inherited=read(PRIOR/'plan.json')['bindings']
    native=read(ROOT/'def-native-cold-coupling/logarithmic-build.json')['bindings']
    for p,h in {**inherited,**native}.items():assert sha(p)==h,p
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='fb4ade4f0',previous_goal_turn='progress',
        claim='Apply native conserved EOS pressure/trace/stress replacement on the actual retained-material history and propagate it into its compact and mass-normalized null-infinity charge.',
        decision='Determine whether the EOS source difference materially changes the surviving scalar residue. A pressure0.2percent or charge2percent failure requires repair rather than physical interpretation.',
        reuse='Reuse old native roots only at exact native conserved coordinates with the same anchor/library, or exactly equal atmospheric U in old support. Recompute changed/new states, including the newly retained tail. Preserve initial pressure offsets. No fluid replay, new grid, clock, horizon or EOS table.',
        scope='All17 saved fine snapshots. Native pressure source replacement holds conserved trajectory, photons and initial GR coefficients fixed; not a uniform EOS/derivative or coupled-gain certificate. Do not inherit old native collisions.',
        budget=dict(pilot_seconds=120,production_seconds=900,worker_seconds=2400,readout_seconds=120,total_seconds=1140,
            workers=3,threads_per_worker=1,memory_GiB_per_process=3,native_calls=80000),
        forecast='Prior3340 native states took348.31s with3 isolated workers and exact-call cache. Require twice measured maximum per-root cost for actual missing high/low/deep groups plus measured setup and60s margin within900s and2400 worker-seconds. No dispatch if ineligible.',
        gates=dict(energy=1e-12,momentum=1e-12,density=1e-10,reseed=1e-10,pressure=.002,charge=.02,cadence=.002,quadrature=.002,linear=1e-10),
        stop='Stop on native/root/support, pressure, forecast or resource failure. Reuse successful new roots; no automatic budget or support expansion.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths},
        inherited_bindings=inherited,native_bindings=native))
    print(json.dumps(dict(prepared=True,new_physical_paths=0,maximum_workers=3,cap_seconds=1140)),flush=True)


def load_inputs():
    model=run.Capture();f=model.flow;e=model.bulk.eos
    z=load(EV/'coupled-128.npz');names=['snapshot_t','snapshot_U','snapshot_Pi','snapshot_h','snapshot_j','snapshot_theta','snapshot_eta','snapshot_u']
    z={key:z[key] for key in names};ds=[];ats=[]
    for it,t in enumerate(z['snapshot_t']):
        history.restore(model,z,it);theta=z['snapshot_theta'][it];eta=z['snapshot_eta'][it];p,u,ut,*_=e.gas(theta,eta)
        ds.append(dict(t=float(t),x=np.log1p(e.density_shift)+np.log1p(e.x),xi=e.xi.copy(),
            lt=np.log(e.base.d['T'])+e.theta0+theta,y=e.base.d['y0']*(1+eta),
            target=z['snapshot_u'][it]+e.inventory[1]*e.xi,p=p.copy(),u=z['snapshot_u'][it],ut=ut.copy()))
        U=z['snapshot_U'][it];rho,v,lt,y=f.primitive(U);p,u,_,_,_,cv,_=f.eos.evaluate(rho,lt)
        active=np.flatnonzero(rho>=f.eos.floor);assert np.array_equal(active,np.flatnonzero(U[0]>=f.eos.floor))
        ats.append(dict(U=U,rho=rho,v=v,lt=lt,y=y,p=p,u=u,cv=cv,
            tau=(U[2]-(model.m.a-model.m.a0)*f.eos.cx*U[0])/model.m.a,active=active))
    return model,ds,ats


def roots(kind,items,label,cap=1200):
    model,ds,ats=INPUTS;mark=time.monotonic();rows=[];durations=[]
    n=run.support.fast_native(cap) if kind=='low' else cached.CachedNative(cap)
    if kind=='low':n.state=lru_cache(maxsize=16)(n.state)
    n.prefix=dict(n.prefix);setup=time.monotonic()-mark;anchor=None;reseed=[]
    try:
        for it,j in items:
            if kind=='deep' and anchor!=j:cached.setup(n,model.bulk.eos.base.d,j);anchor=j
            begin=time.monotonic()
            r=history.deep.root(n,ds[it],j) if kind=='deep' else history.atmosphere.root(n,model.flow,ats[it],j)
            durations.append(time.monotonic()-begin)
            r=history.decorate('deep' if kind=='deep' else 'atmosphere',it,j,r,model,ds,ats)
            rows.append(dict(r,native_owner=kind,reused=False))
            if len(rows)%64==0:write(OUT/f'{label}-partial.json',rows)
        if label.startswith('pilot-'):
            other=history.deep.root(n,ds[it],j,.001) if kind=='deep' else history.atmosphere.root(n,model.flow,ats[it],j,.001)
            error=abs(other['raw'][1]/r['raw'][1]-1);assert error<1e-10
            reseed.append(dict(it=it,cell=j,relative=error,offset_logT=.001,root=other))
        write(OUT/f'{label}-samples.json',rows)
        result=dict(label=label,kind=kind,states=len(rows),setup_seconds=setup,root_seconds=durations,
            native_calls=n.ion.calls+getattr(n,'variant_initial_calls',0),reseed=reseed,seconds=time.monotonic()-mark)
        write(OUT/f'{label}.json',result);return result
    except Exception as exc:
        write(OUT/f'{label}-partial.json',rows);write(OUT/f'{label}-failure.json',dict(error=repr(exc),kind=kind,
            it=it,cell=j,completed=len(rows),native_calls=n.ion.calls+getattr(n,'variant_initial_calls',0),seconds=time.monotonic()-mark));raise


def worker(job):return roots(*job)


def evaluate():
    global np,run,history,cached,load,INPUTS
    start=time.monotonic();plan=read(OUT/'plan.json')
    for group in ['bindings','inherited_bindings','native_bindings']:
        for p,h in plan[group].items():
            assert sha(p)==h,p
            if '/verification/' in p:assert sha(Path('verification')/Path(p).name)==h,p
    import numpy as np
    import evolve_native_retained_tail as run
    import def_native_conserved_history_charge as history
    import verify_native_conserved_history_charge as cached
    import multiprocessing as mp
    import sympy as sp
    load=lambda p:np.load(io.BytesIO(p.read_bytes()))
    run.support.owner.deadline(start,120);INPUTS=load_inputs();model,ds,ats=INPUTS
    assert len(ds)==17 and all(abs(v['t']-ds[-1]['t']*it/16)<1e-17 for it,v in enumerate(ds))
    previous={(r['kind'],r['it'],r['cell']):r for r in read(PRIOR/'production-samples.json')}
    oldU=load(THERMAL/'updated-gr-return/background/coupled-128.npz')['snapshot_U']
    rows=[];missing={'deep':[],'high':[],'low':[]};reused={'deep':0,'atmosphere':0}
    try:
        for it in range(17):
            for j in range(19):
                r=previous['deep',it,j]
                if all(float(ds[it][a][j])==r[b] for a,b in [('x','x'),('y','y'),('target','target_u')]):
                    rows.append(dict(history.decorate('deep',it,j,r,model,ds,ats),reused=True));reused['deep']+=1
                else:missing['deep'].append((it,j))
            for j in map(int,ats[it]['active']):
                low=ats[it]['rho'][j]<model.flow.eos.original.floor;r=previous.get(('atmosphere',it,j))
                if not low and r is not None and np.array_equal(oldU[it,:,j],ats[it]['U'][:,j]):
                    rows.append(dict(history.decorate('atmosphere',it,j,r,model,ds,ats),reused=True));reused['atmosphere']+=1
                else:missing['low' if low else 'high'].append((it,j))
        write(OUT/'reuse.json',dict(classification='Counterexample candidate',exact_reused=reused,
            missing={k:len(v) for k,v in missing.items()},active_atmosphere_per_time=[len(v['active']) for v in ats],
            no_near_match_reuse=True,seconds=time.monotonic()-start))
        write(OUT/'reused-samples.json',rows);print((OUT/'reuse.json').read_text(),flush=True)
        pilots={};pilot_calls=0
        for kind,items in missing.items():
            if not items:continue
            choices=[items[i] for i in sorted(set([0,len(items)//2,len(items)-1]))]
            result=roots(kind,choices,'pilot-'+kind);pilots[kind]=result;pilot_calls+=result['native_calls']
            rows+=read(OUT/f'pilot-{kind}-samples.json');missing[kind]=[v for v in items if v not in choices]
        work={k:2*max(pilots[k]['root_seconds'])*len(v) if v else 0 for k,v in missing.items()}
        setup=max((r['setup_seconds'] for r in pilots.values()),default=0)
        upper=work['deep']+work['low']+work['high']/3+2*setup+60
        eligible=upper<900 and sum(work.values())<2400
        pilot=dict(classification='Counterexample candidate',passed=True,eligible=eligible,regions=pilots,
            exact_reused=reused,remaining={k:len(v) for k,v in missing.items()},upper_wall_seconds=upper,
            upper_worker_seconds=sum(work.values()),seconds=time.monotonic()-start,native_calls=pilot_calls)
        write(OUT/'pilot.json',pilot);run.support.owner.deadline(start,120);print(json.dumps(pilot),flush=True);assert eligible
        mark=time.monotonic();run.support.owner.deadline(mark,900)
        jobs=[]
        if missing['deep']:jobs.append(('deep',sorted(missing['deep'],key=lambda v:(v[1],v[0])),'deep'))
        for group in range(3):
            if missing['high'][group::3]:jobs.append(('high',missing['high'][group::3],f'high-{group}'))
        if missing['low']:jobs.append(('low',missing['low'],'low'))
        # Native Fortran state is isolated per process, never shared by threads.
        results=[]
        if jobs:
            cap=min(18000,(80000-pilot_calls-1)//len(jobs))
            jobs=[(*job,cap) for job in jobs]
            with mp.get_context('fork').Pool(min(3,len(jobs))) as pool:results=pool.map(worker,jobs)
        for r in results:rows+=read(OUT/f"{r['label']}-samples.json")
        calls=pilot_calls+sum(r['native_calls'] for r in results);worker_time=sum(r['seconds'] for r in results)
        assert calls<80000 and worker_time<2400
        keys={(r['kind'],r['it'],r['cell']) for r in rows}
        expected={('deep',it,j) for it in range(17) for j in range(19)}|{('atmosphere',it,int(j)) for it,v in enumerate(ats) for j in v['active']}
        assert keys==expected and len(rows)==len(keys)
        write(OUT/'production-samples.json',rows)
        production=dict(classification='Counterexample candidate',passed=True,states=len(rows),exact_reused=reused,
            native_calls=calls,worker_seconds=worker_time,seconds=time.monotonic()-mark,
            maximum_pressure_relative=max(abs(r['pressure_relative']) for r in rows),
            maximum_energy_relative=max(abs(r.get('energy_relative',r.get('energy_residual'))) for r in rows),
            maximum_momentum_relative=max(r.get('momentum_relative',0.) for r in rows),
            maximum_density_relative=max(r['density_relative'] for r in rows),jobs=results)
        write(OUT/'production.json',production);run.support.owner.deadline(mark,900);print(json.dumps(production),flush=True)
        readout(model,ds,ats,rows)
        run.support.owner.deadline(start,1140)
    except Exception as exc:write(OUT/'failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def readout(model,ds,ats,rows):
    start=time.monotonic();run.support.owner.deadline(start,120);charge=history.charge
    physical=SimpleNamespace(**dict(vars(history.physical),EV=EV))
    source=FunctionType(history.source.__code__,dict(history.source.__globals__,physical=physical))
    m=charge.gr.Response();waves={}
    for label,stride,region,order in [('fine',1,'all',8),('nine',2,'all',8),('g4',1,'all',4),('deep',1,'deep',8)]:
        d,knots=source(model,ds,ats,rows,stride,region);waves[label]=charge.read(m,d,order)
        run.atomic_npz(OUT/f'wave-{label}.npz',**waves[label])
        if label=='fine':run.atomic_npz(OUT/'source.npz',**d);run.atomic_npz(OUT/'source-knots.npz',**knots)
    d=dict(load(OUT/'source.npz'));total=dict(load(EV/'source-128.npz'))
    for key in physical.capture.KEYS+['metric_stress_erg']:total[key]=total[key]+d[key]
    applied=charge.read(m,total,8);before=load(CURRENT/'gr/wave-128-g8.npz');a=waves['fine']['free_scalar'];norm=float(np.max(abs(before['free_scalar'])))
    linear=float(np.max(abs(applied['free_scalar']-before['free_scalar']-a)))/norm
    cadence=float(np.max(abs(a-waves['nine']['free_scalar'])))/norm;quad=float(np.max(abs(a-waves['g4']['free_scalar'])))/norm
    direct,error=charge.independent.direct(m,d,8);agreement=float(abs(direct-waves['fine']['direct_scalar'][-1]))/max(float(np.max(abs(waves['fine']['direct_scalar']))),1e-300)
    materiality=float(np.max(abs(a)))/norm
    # Rebuild conserved native energy and sources from the recorded raw states.
    clone=SimpleNamespace(**dict(vars(history),OUT=OUT,load_inputs=lambda:(model,ds,ats)))
    # The inherited audit shares this action's absolute deadline.
    audit_signal=SimpleNamespace(alarm=lambda _:run.support.owner.deadline(start,120))
    audit=FunctionType(cached.audit.__code__,dict(cached.audit.__globals__,OUT=OUT,run=clone,signal=audit_signal))
    audit();run.support.owner.deadline(start,120)
    import sympy as sp
    S,Q,p,dp=sp.symbols('S Q p dp');v=S/(Q+p);vn=S/(Q+p+dp)
    assert sp.factor(-S*(vn-v)-3*dp-(v*vn-3)*dp)==0
    assert sp.factor(-S*(vn-v)-dp-(v*vn-1)*dp)==0
    infinity=dict(load(CURRENT/'infinity/completed/retained-128-a8-r8.npz'));ids=np.arange(0,129,8)
    correction=a[ids]/(1-infinity['epsilon']);updated=infinity['normalized']+correction
    run.atomic_npz(OUT/'applied-charge.npz',**applied,previous=before['free_scalar'],correction=a)
    run.atomic_npz(OUT/'normalized-charge.npz',t=infinity['t'],previous=infinity['normalized'],correction=correction,normalized=updated)
    passed=linear<1e-10 and cadence<.002 and quad<.002 and agreement<1e-9 and materiality<.02
    result=dict(classification='Counterexample candidate',passed=bool(passed),old_compact_endpoint=float(before['free_scalar'][-1]),
        new_compact_endpoint=float(applied['free_scalar'][-1]),correction_endpoint=float(a[-1]),
        old_normalized_endpoint=float(infinity['normalized'][-1]),new_normalized_endpoint=float(updated[-1]),
        deep_correction_endpoint=float(waves['deep']['free_scalar'][-1]),atmosphere_correction_endpoint=float(a[-1]-waves['deep']['free_scalar'][-1]),
        maximum_correction_relative=materiality,cadence_relative_to_signal=cadence,quadrature_relative_to_signal=quad,
        direct_application_relative=linear,independent_direct_relative=agreement,inverse_radius_residual=float(error),
        original_initial_offset_retained=True,symbolic_conserved_trace_checked=True,seconds=time.monotonic()-start,
        actual_new_native_source_applied=True,photons_and_conserved_trajectory_held_fixed=True,
        native_correction_time_convergence_tested=False,uniform_EOS_derivative_bound=False,full_floor_feedback_enclosed=False,
        coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert passed


def launch():
    assert not (OUT/'receipt.json').exists();start=time.monotonic()
    with (OUT/'external.log').open('w') as stream:
        result=subprocess.run(['timeout','--signal=TERM','--kill-after=5s','1140s',sys.executable,__file__,'evaluate'],stdout=stream,stderr=subprocess.STDOUT)
    write(OUT/'receipt.json',dict(returncode=result.returncode,seconds=time.monotonic()-start,external_cap_seconds=1140))
    print((OUT/'external.log').read_text(),end='');print((OUT/'receipt.json').read_text(),flush=True);sys.exit(result.returncode)


if __name__=='__main__':
    if sys.platform!='win32':
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    globals()[sys.argv[1]]()
