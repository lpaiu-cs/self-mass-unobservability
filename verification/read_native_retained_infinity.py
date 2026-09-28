"""Read the retained-material emission at infinity within the unused budget."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import hashlib,io,json,signal,subprocess,sys,time

ROOT=Path('outputs/direct-eos-gr33');RUN=ROOT/'native-retained-completion';FAILED=RUN/'infinity';OUT=FAILED/'completed'
BEFORE=ROOT/'native-retained-tail/supported-temperature/evolution';EV=RUN/'evolution'
THERMAL=ROOT/'def-native-conservative-rates/thermal-refined';OLDPORT=THERMAL/'updated-gr-return/background'
read=lambda p:json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def prepare():
    assert read(EV/'result.json')['passed'] and read(RUN/'gr/result.json')['passed']
    assert not OUT.exists();OUT.mkdir()
    failed=read(FAILED/'receipt.json');assert failed['returncode']==1
    spent=sum(read(RUN/f'{name}-receipt.json')['seconds'] for name in ['evolve','charge'])+failed['seconds']
    remaining=read(RUN/'plan.json')['budget']['total_scientific_seconds']-spent
    previous=read(ROOT/'native-charge-null-infinity/production.json')['seconds']
    setup=read(BEFORE/'cost-64-64.json')['cost']['setup'];upper=2*(previous+setup)
    cap=min(120.-failed['seconds'],remaining);assert upper<cap,('Unspent aggregate budget gate',upper,cap)
    paths=[Path(__file__),Path('verification/def_native_global_scalar_closure.py'),Path('verification/verify_native_stage_energy_charge.py'),
        RUN/'gr/result.json',ROOT/'def-native-global-scalar-closure/bound.json',ROOT/'def-native-global-scalar-closure/bound-plan.json',
        FAILED/'plan.json',FAILED/'receipt.json',FAILED/'failure.json',FAILED/'executed-source.py']
    for n,folder in [(64,BEFORE),(128,EV)]:
        paths += [folder/f'accepted-ports-{n}.npz',folder/f'source-{n}.npz',folder/f'history-{n}.npz',
            RUN/f'gr/wave-{n}-g8.npz',OLDPORT/f'accepted-ports-{n}.npz',THERMAL/f'evolution/source-{n}.npz',THERMAL/f'evolution/gr/wave-{n}-g8.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',claim='Propagate the actual lower-floor emission to null infinity and pair exterior stress, body energy debit and arrived-photon mass normalization. Compare to the same thermal background without later native corrections.',
        budget_seconds=cap,prior_actions_seconds=spent,aggregate_cap_seconds=1650,
        prior_failure='Exact manifest lookup used the Linux code-cache alias rather than the original input path. No kernel or emission evaluation was reached. Resolve both to the same existing file, retain original SHA comparison, and use only the remainder of the same120s allocation.',
        measured_upper_seconds=upper,allocation='Only unused1650s Phase148 allocation. Prior exact kernel/readout17.334s plus a32.204s setup allowance, doubled. No new physical evolution or automatic extension.',
        settings=[[128,8,8],[128,4,8],[128,8,4],[64,8,8]],matched_reference_settings=[[128,8,8],[64,8,8]],
        gates=dict(port=1e-12,quadrature=.002,time=.02,identity=1e-12),
        bounds='Rebuild original conditional compact-source envelope on this source/history. Retain its larger exterior envelope even though exterior forcing is explicitly computed. Includes frozen potential, specified port and direct omitted positive-matter terms; does not bound later floor feedback, microscopic or nonlinear errors.',
        stop='Stop on failed compact readout, budget or original quadrature/time controls; no new clock, EOS state, kernel order or physical horizon.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths}))
    print(json.dumps(dict(prepared=True,cap_seconds=cap,upper_seconds=upper,aggregate_remaining=remaining)),flush=True)


def compute():
    start=time.monotonic();plan=read(OUT/'plan.json')
    for p,h in plan['bindings'].items():assert sha(Path(p))==h,p
    import numpy as np
    import sympy as sp
    import evolve_native_retained_tail as run
    import def_native_global_scalar_closure as original
    import verify_native_stage_energy_charge as bound_owner
    run.support.owner.deadline(start,plan['budget_seconds']);LD=np.longdouble;C,G=original.C,original.G
    load=lambda p:np.load(io.BytesIO(p.read_bytes()))
    try:
        m=original.Exterior();cache={};rows=[]
        geom=read(original.OUT/'bound-plan.json')['bindings'];bg=original.flow.INPUT/'balanced-20.npz'
        matches=[name for name in geom if Path(name).exists() and Path(name).resolve()==bg.resolve()]
        assert len(matches)==1;assert sha(bg)==geom[matches[0]]
        bound_input=Path(matches[0]).parent
        def inputs(n,baseline):
            key=(n,baseline)
            if key in cache:return cache[key]
            folder=OLDPORT if baseline else (BEFORE if n==64 else EV)
            z=load(folder/f'accepted-ports-{n}.npz');lum=z['angular_luminosity'].astype(LD);h=float(z['h'])
            assert lum.shape==(n,4) and np.min(lum)>=0 and abs(n*h-m.T)<1e-18
            q=np.vstack([np.zeros(4,LD),np.cumsum(LD(h)*lum,axis=0)])
            r=np.vstack([np.zeros(4,LD),np.cumsum(LD(h)*q[:-1]+LD(h)**2*lum/2,axis=0)])
            sf=THERMAL/'evolution' if baseline else folder;source=load(sf/f'source-{n}.npz')
            aw=np.arange(1,8,2,dtype=LD)/32;port=source['outer_cumulative_energy_erg']
            error=float(np.max(abs(q@aw-port))/max(abs(port[-1]),1.));assert error<1e-12
            def primitive(at):
                x=np.clip(np.asarray(at),0,m.T);j=np.clip((x/h).astype(int),0,n-1);s=(x-j*h).astype(LD)[...,None]
                return q[j]+s*lum[j],r[j]+s*q[j]+s*s*lum[j]/2
            cache[key]=(primitive,error);return cache[key]
        def evaluate(n,a,r,baseline=False):
            mark=time.monotonic();primitive,port=inputs(n,baseline);k=m.kernel(a,r);mass=[];stress=[];arrival=[]
            for t in m.t:
                h,hh=primitive(np.maximum(t-k['delay'],0.));ids=np.arange(len(h));bins=k['node_bins']
                mass.append(float(-G/C**3*np.sum(k['weights']*k['mass']*hh[ids,bins],dtype=LD)))
                stress.append(float(-G/(2*C**4)*np.sum(k['weights']*k['stress']*h[ids,bins],dtype=LD)))
                h,_=primitive(np.maximum(t-k['infinity'],0.));arrival.append(float(np.sum(k['mw']*k['mu']*h[np.arange(len(h)),k['bins']],dtype=LD)))
            data=dict(t=m.t,mass_U=np.asarray(mass),stress_U=np.asarray(stress),
                exterior=-(np.asarray(mass)+stress)/m.M,arrived=np.asarray(arrival))
            wave=load((THERMAL/'evolution/gr' if baseline else RUN/'gr')/f'wave-{n}-g8.npz')
            ids=np.array([int(np.argmin(abs(wave['t']-t))) for t in m.t]);assert max(abs(wave['t'][ids]-m.t))<1e-18
            compact=wave['free_scalar'][ids];eps=G/C**4*data['arrived']/m.M;alpha=-m.K/m.M
            assert np.all(eps>=0) and max(eps)<1
            scalar=compact+data['exterior'];data.update(compact=compact,epsilon=eps,scalar=scalar,
                mass_term=alpha*eps,normalized=(scalar+alpha*eps)/(1-eps))
            label=f'{"baseline" if baseline else "retained"}-{n}-a{a}-r{r}'
            run.atomic_npz(OUT/f'{label}.npz',**data)
            rows.append(dict(label=label,port_relative=port,seconds=time.monotonic()-mark,
                endpoint=float(data['normalized'][-1]),delay_inverse_seconds=k['inversion_error']))
            run.support.owner.deadline(start,plan['budget_seconds']);return data
        paths=[evaluate(*v) for v in plan['settings']];refs={n:evaluate(n,a,r,True) for n,a,r in plan['matched_reference_settings']}
        fine=paths[0];controls={}
        for field in ['exterior','arrived','scalar','normalized']:
            norm=max(float(np.max(abs(fine[field]))),1e-300)
            controls[field]={key:float(np.max(abs(d[field]-fine[field])))/norm for key,d in zip(['angular','radial','time'],paths[1:])}
        delta=fine['normalized']-refs[128]['normalized'];delta64=paths[-1]['normalized']-refs[64]['normalized']
        delta_time=float(np.max(abs(delta-delta64)))/max(float(np.max(abs(delta))),1e-300)
        # Same existing file and SHA; only the path spelling differs in cache.
        bound_flow=SimpleNamespace(**dict(vars(bound_owner.flow),INPUT=bound_input))
        enclosure=FunctionType(bound_owner.enclosure.__code__,dict(bound_owner.enclosure.__globals__,OUT=OUT,flow=bound_flow))
        bound=enclosure(m.response,dict(load(EV/'source-128.npz')),load(EV/'history-128.npz'))
        error=float(np.nextafter(bound['total_conditional_normalized_bound']/(1-fine['epsilon'][-1]),np.inf))
        interval=[float(np.nextafter(fine['normalized'][-1]-error,-np.inf)),float(np.nextafter(fine['normalized'][-1]+error,np.inf))]
        A,s,e=sp.symbols('A s e');assert sp.simplify((A+s)/(1-e)-A-(s+A*e)/(1-e))==0
        passed=all(v['angular']<.002 and v['radial']<.002 and v['time']<.02 for v in controls.values())
        run.support.owner.deadline(start,plan['budget_seconds'])
        result=dict(classification='Counterexample candidate',passed=bool(passed),controls=controls,rows=rows,
            endpoint_normalized=float(fine['normalized'][-1]),endpoint_scalar=float(fine['scalar'][-1]),
            endpoint_compact=float(fine['compact'][-1]),endpoint_exterior=float(fine['exterior'][-1]),
            endpoint_mass_term=float(fine['mass_term'][-1]),matched_baseline_endpoint=float(refs[128]['normalized'][-1]),
            retained_material_change=float(delta[-1]),retained_change_time_relative=delta_time,
            retained_change_resolved_on_time_gate=delta_time<.02,conditional_interval=interval,conditional_envelope=bound,
            actual_current_emission_at_null_infinity=True,actual_emission_mass_normalization=True,
            normalization_identity_symbolic=True,seconds=time.monotonic()-start,
            old_native_collision_corrections_inherited=False,uniform_EOS_derivative_bound=False,
            full_floor_feedback_enclosed=False,coupled_fixed_point_verified=False,nonlinear_GR=False,
            final_charge_solved=False,full_goal_complete=False)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:write(OUT/'failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def launch():
    cap=read(OUT/'plan.json')['budget_seconds'];start=time.monotonic();assert not (OUT/'receipt.json').exists()
    with (OUT/'external.log').open('w') as stream:
        result=subprocess.run(['timeout','--signal=TERM','--kill-after=5s',str(cap)+'s',sys.executable,__file__,'compute'],stdout=stream,stderr=subprocess.STDOUT)
    write(OUT/'receipt.json',dict(returncode=result.returncode,seconds=time.monotonic()-start,external_cap_seconds=cap))
    print((OUT/'external.log').read_text(),end='');print((OUT/'receipt.json').read_text(),flush=True);sys.exit(result.returncode)


if __name__=='__main__':
    if sys.platform!='win32':
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    globals()[sys.argv[1]]()
