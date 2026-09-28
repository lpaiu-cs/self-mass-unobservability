"""Fixed-inventory native acoustic derivatives on the retained trajectory.

Counterexample candidate. Sampled native derivatives and their actual flux
application; never a uniform native-EOS derivative certificate.
"""
from pathlib import Path
import inspect, json, os, signal, sys, time, multiprocessing as mp
from types import FunctionType, SimpleNamespace
import numpy as np
import sympy as sp
os.environ['RETAINED_NATIVE_OUTPUT']='retained-native-return150-work'
import def_retained_native_return as prior
import verify_native_conserved_history_charge as cached

OUT=Path('retained-native-acoustic153-work')
HYDRO=OUT/'hydro'
read,write,sha=prior.read,prior.write,prior.sha


def deadline(seconds):
    def stop(*_):raise TimeoutError('Registered native acoustic budget')
    signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,seconds)


def prepare():
    assert not OUT.exists();HYDRO.mkdir(parents=True)
    paths=[Path(__file__),Path(prior.__file__),prior.EOS/'production-samples.json',
           Path('verification/def_native_hydrogen_exchange.py'),
           Path('verification/def_native_material_join.py'),
           Path('verification/def_native_energy_flow.py'),
           prior.HYDRO/'production.json',Path('retained-metric-return152-work/final-result.json')]
    paths=[p for p in paths if p.exists()]
    paths+=sorted(prior.HYDRO.glob('*-native-states.json'))
    paths+=sorted(prior.HYDRO.glob('point-*.npz'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Replace the inherited atmosphere gamma and deep face_K by fixed-inventory native thermodynamic derivatives, apply their conservative shared-face force difference, and measure the actual retained-history charge response.',
        decision='A changed acoustic force that changes charge by2percent requires constitutive repair before interpreting the residue. Passing only derivative samples does not complete this lever: apply the force to the existing response paths.',
        scope='Same conserved background, physical amplitude,531cells,17source knots,64/128 clocks and3.4344311179287023ms. No background replay, new mesh, longer horizon or full nonlinear Einstein claim.',
        derivative='Differentiate native constrained pressure, energy and entropy in log density and log temperature while holding neutral H and all other declared inventories fixed. Never borrow equilibrium derivative columns. Centered h=1e-4 and h/2 with Richardson extrapolation; preserve both estimates.',
        budget=dict(pilot_seconds=100,native_production_seconds=1000,native_worker_seconds=2600,
                    native_calls=350000,force_seconds=120,material_seconds=500,
                    photon_seconds=1200,readout_seconds=180,total_seconds=3100,workers=3,
                    BLAS_threads=1,virtual_GiB_per_process=3),
        gates=dict(derivative=.0002,first_law=.0001,conservation=1e-8,pressure=.002,
                   time=.02,quadrature=.002,independent_GR=1e-9),
        forecast='Representative native h/h2 states first; dispatch only if2x worst measured cost plus setup fits1000s wall and2600 worker-seconds. Response paths require measured prefixes and previous full late-CFL counts within their caps.',
        stop='Preserve failures; no automatic support, grid, clock, horizon, amplitude or gate changes. Reuse completed native centers and accepted prefixes.',
        limits='Finite difference convergence is sampled evidence, not a continuous enclosure. The atmospheric derivative uses current fixed-H chemistry; deep additive advected-inventory terms remain explicit. Current native collision Jacobians and full coupled gain remain separate boundaries.',
        bindings={str(p):sha(p) for p in paths}))
    px,pt,ux,ut,p,rho=sp.symbols('Px PT ux uT P rho',nonzero=True)
    tx=(p/rho-ux)/ut;K=px+pt*tx
    assert sp.simplify(ux+ut*tx-p/rho)==0
    assert sp.simplify((px+pt*tx)/p-K/p)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        premise='Differentiable constrained EOS at fixed specific chemical inventories; reversible adiabatic first law du=P/rho dlnrho.',
        temperature_slope='(P/rho-u_x)/u_T',bulk_modulus='P_x+P_T*(P/rho-u_x)/u_T',
        gamma='bulk_modulus/P',sound_speed_squared_over_c_squared='bulk_modulus/(rho*(cx*c^2+u)+P)'))


def derivative(native,key,raw,physical_pressure=None):
    x,t,y=key;h=1e-4;values=[];begin=time.monotonic();calls=native.ion.calls
    for d in [h,h/2]:
        values.append(np.array([native.state(x+d,t,y)['raw'],native.state(x-d,t,y)['raw'],
                                native.state(x,t+d,y)['raw'],native.state(x,t-d,y)['raw']]))
    a=np.array(values);dr=np.array([(v[0]-v[1])/(2*d) for v,d in zip(a,[h,h/2])])
    dt=np.array([(v[2]-v[3])/(2*d) for v,d in zip(a,[h,h/2])])
    p=float(raw[1] if physical_pressure is None else physical_pressure);r=float(raw[0]);T=np.exp(t)
    gam=(dr[:,1]+dt[:,1]*(p/r-dr[:,2])/dt[:,2])/p
    dx=(4*dr[1]-dr[0])/3;dT=(4*dt[1]-dt[0])/3
    K=dx[1]+dT[1]*(p/r-dx[2])/dT[2]
    errors=[abs(gam[1]/gam[0]-1),abs(dt[1,2]/dt[0,2]-1)]
    # Native first-law check uses native P, before explicit inventory offsets.
    laws=[(T*dT[3]-dT[2])/dT[2],(T*dx[3]-dx[2]+raw[1]/r)/(raw[1]/r)]
    assert dT[2]>0 and K>0 and max(errors)<.0002,(key,errors)
    assert max(abs(np.array(laws)))<.0001,(key,laws)
    return dict(coordinates=list(key),raw=list(map(float,raw)),samples=a.tolist(),
                dx=dx.tolist(),dT=dT.tolist(),gamma=float(K/p),K=float(K),
                physical_pressure=p,comparison=list(map(float,errors)),first_law=list(map(float,laws)),
                calls=native.ion.calls-calls,seconds=time.monotonic()-begin)


def warm_start(native):
    original=native.state
    def state(x,t,y):
        result=original(x,t,y)
        fields=native.ion.fields.copy()
        fields[0]-=np.log((1-y)/y)-np.log((1-native.y0)/native.y0)
        native.prefix=dict(T=np.array([np.exp(t)]),log_density_ratio=np.array([x]),fields=fields[None])
        return result
    # The converged fields are only an initial guess. Full target and native
    # residual checks still run at every coordinate; no nearby reuse of values.
    native.state=state
    return native


def warm_pilot():
    old=read(OUT/'pilot.json');assert old['passed'] and not old['eligible']
    assert not (OUT/'warm-pilot.json').exists();deadline(85);start=time.monotonic();prior.initialize()
    model=prior.current.Capture();data=model.bulk.eos.base.d
    high=cached.CachedNative(12000);high.prefix=dict(high.prefix)
    low=warm_start(prior.current.support.fast_native(12000));low.prefix=dict(low.prefix)
    rows=[];setup_calls=high.ion.calls+low.ion.calls
    write(OUT/'warm-plan.json',dict(classification='Counterexample candidate',
        failure='Original science pilot passed but2x5080worker/1723wall forecast exceeds2600/1000; no production dispatched.',
        repair='Use the immediately preceding converged constrained chemical fields solely as the next native solve initial guess. Preserve exact coordinates, target inventories, eight calls per derivative and every original gate. Compare all12 original derivative probes.',
        remaining_pilot_seconds=85,production_seconds=1000,worker_seconds=2600,new_budget=False,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'registered-pilot-producer.py',OUT/'pilot.json',OUT/'pilot-samples.json']}))
    try:
        for r in read(OUT/'pilot-samples.json'):
            if r['kind']=='deep':
                # Warm prefix belongs to one chemical anchor only.
                high=cached.CachedNative(12000);cached.setup(high,data,r['cell']);setup_calls+=high.ion.calls
                n=warm_start(high)
            else:
                if not rows or rows[-1]['kind']=='deep':
                    high=warm_start(cached.CachedNative(12000));setup_calls+=high.ion.calls
                n=low if np.exp(r['coordinates'][0])<model.flow.eos.original.floor else high
            v=derivative(n,tuple(r['coordinates']),np.array(r['raw']),r['physical_pressure'])
            raw_error=float(np.max(abs(np.array(v['samples'])[:,:,[0,1,2,3]]-np.array(r['samples'])[:,:,[0,1,2,3]])/np.maximum(abs(np.array(r['samples'])[:,:,[0,1,2,3]]),1.)))
            assert raw_error<1e-10 and abs(v['gamma']/r['gamma']-1)<.0002
            rows.append(dict(v,kind=r['kind'],it=r['it'],cell=r['cell'],cold_start_raw_relative=raw_error))
            write(OUT/'warm-pilot-samples.json',rows)
        costs={k:max(r['seconds'] for r in rows if r['kind']==k) for k in old['counts']}
        work=2*sum(costs[k]*old['counts'][k] for k in costs);wall=work/3+30
        result=dict(classification='Counterexample candidate',passed=True,eligible=work<2600 and wall<1000,
                    costs=costs,estimated_worker_seconds=work,upper_wall_seconds=wall,seconds=time.monotonic()-start,
                    calls=sum(r['calls'] for r in rows),setup_calls=setup_calls,
                    cold_start_max_relative=max(r['cold_start_raw_relative'] for r in rows))
        write(OUT/'warm-pilot.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'warm-failure.json',dict(error=repr(exc),completed=len(rows),seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def job(arg):
    kind,label,rows=arg;start=time.monotonic();cpu=time.process_time();deadline(500)
    n=prior.current.support.fast_native(16000) if kind=='low' else cached.CachedNative(16000)
    if kind=='deep':cached.setup(n,prior.current.Capture().bulk.eos.base.d,rows[0]['cell'])
    n=warm_start(n);setup_calls=n.ion.calls;results=[]
    try:
        for r in rows:
            key=tuple(r['coordinates'])
            results.append(derivative(n,key,np.array(r['raw']),r.get('pressure')))
        # Dense native probes are numerical artifacts, stored once without
        # decimal JSON inflation. Metadata stays inspectable and SHA-bound.
        arrays={k:np.array([r[k] for r in results]) for k in results[0]}
        np.savez_compressed(HYDRO/f'{label}.npz',**arrays)
        report=dict(kind=kind,label=label,states=len(results),seconds=time.monotonic()-start,
                    CPU_seconds=time.process_time()-cpu,calls=n.ion.calls,setup_calls=setup_calls,
                    derivative=max(max(r['comparison']) for r in results),
                    first_law=max(max(abs(np.array(r['first_law']))) for r in results))
        write(HYDRO/f'{label}.json',report);return report
    except Exception as exc:
        write(HYDRO/f'{label}-failure.json',dict(error=repr(exc),completed=len(results),calls=n.ion.calls,seconds=time.monotonic()-start));raise


def native():
    assert not (OUT/'native.json').exists();p=read(OUT/'warm-pilot.json');assert p['passed']
    prior.initialize();faces={}
    for path in prior.HYDRO.glob('*-native-states.json'):
        for r in read(path):faces[tuple(r['coordinates'])]=r
    roots=read(prior.EOS/'production-samples.json');jobs=[]
    for kind in ['high','low']:
        rows=sorted([r for r in faces.values() if r['owner']==kind],key=lambda r:r['coordinates'])
        for k in range(0,len(rows),200):jobs.append((kind,f'{kind}-{k//200}',rows[k:k+200]))
    for j in range(19):
        rows=[dict(r,coordinates=[r['x'],r['lt'],r['y']]) for r in roots if r['kind']=='deep' and r['cell']==j]
        jobs.append(('deep',f'deep-{j}',rows))
    forecast=p['upper_wall_seconds'];assert forecast<1500 and p['estimated_worker_seconds']<4500
    write(OUT/'native-execution-plan.json',dict(classification='Counterexample candidate',
        reassessment='The completed warm pilot gives the SAME eight-probe derivatives. Reserve1500s wall/4500worker for the complete native face/deep calculation by moving500s of unstarted response reserve; total wall3100s and scientific gates stay unchanged. After completion, only measured unused time may be returned to response work.',
        acoustic_scope='Native derivatives at all6392 reconstructed face states and323 deep states. Cell-centered atmospheric gamma is used only by the legacy CFL estimate, not the force; supplement that estimate with a conservative cap from the maximum actual native face speed. No extra3191 center probes that cannot change the semidiscrete force.',
        total_wall_seconds=3100,native_wall_seconds=1500,native_worker_seconds=4500,forecast_wall_seconds=forecast,
        native_call_cap=350000,workers=3,jobs=[dict(kind=k,label=l,states=len(r)) for k,l,r in jobs],
        bindings={str(f):sha(f) for f in [Path(__file__),OUT/'registered-warm-producer.py',OUT/'warm-pilot.json']}))
    start=time.monotonic();deadline(1500)
    with mp.get_context('fork').Pool(3) as pool:reports=pool.map(job,jobs)
    seconds=time.monotonic()-start;cpu=sum(r['CPU_seconds'] for r in reports);calls=sum(r['calls'] for r in reports)
    assert seconds<1500 and sum(r['seconds'] for r in reports)<4500 and calls<350000
    result=dict(classification='Counterexample candidate',passed=True,seconds=seconds,CPU_seconds=cpu,calls=calls,
                states=sum(r['states'] for r in reports),maximum_derivative=max(r['derivative'] for r in reports),
                maximum_first_law=max(r['first_law'] for r in reports),rows=reports,uniform_derivative_enclosure=False)
    write(OUT/'native.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    signal.setitimer(signal.ITIMER_REAL,0.)


def force():
    assert read(OUT/'native.json')['passed'] and not (OUT/'force.json').exists()
    start=time.monotonic();deadline(120);prior.initialize();lookup={};deep={}
    roots=read(prior.EOS/'production-samples.json')
    for r in read(OUT/'native.json')['rows']:
        data=np.load(HYDRO/f"{r['label']}.npz")
        if r['kind']=='deep':
            j=int(r['label'].split('-')[1]);bykey={tuple(k):v for k,v in zip(data['coordinates'],data['K'])}
            states=sorted([v for v in roots if v['kind']=='deep' and v['cell']==j],key=lambda v:v['it'])
            deep[j]=dict(K=np.array([bykey[v['x'],v['lt'],v['y']] for v in states]))
        else:
            for key,gamma in zip(data['coordinates'],data['gamma']):lookup[tuple(key)]=float(gamma)
    # Reuse the original native force owner, with its EXACT native p/u states.
    # Only face gamma and deep K change; the old table-centered CFL estimate
    # is supplemented by the largest actual native face characteristic speed.
    code=inspect.getsource(prior.hydro_sources)
    code=prior.replace(code,"    native_calls0=high.ion.calls+low.ion.calls", "    native_calls0=high.ion.calls+low.ion.calls\n    acoustic=dict(face_calls=0,center_calls=0,max_speed=0.,max_gamma_change=0.,max_deep_K_change=0.)\n    baselineK=model.face_K.copy()")
    code=prior.replace(code,'            return tuple(values)',
        "            for j in np.flatnonzero(rho>=self.floor):\n                key=(float(np.log(rho[j])),float(lt[j]),float(y[j]))\n                if len(rho)==m.n-m.nb:\n                    acoustic['center_calls']+=1\n                    continue\n                gamma=lookup[key];acoustic['face_calls']+=1\n                acoustic['max_gamma_change']=max(acoustic['max_gamma_change'],abs(gamma/values[2,j]-1))\n                values[2,j]=gamma\n                cs=np.sqrt(gamma*values[0,j]/(rho[j]*(self.cx+values[1,j])+values[0,j]))\n                acoustic['max_speed']=max(acoustic['max_speed'],float(cs))\n            return tuple(values)")
    code=prior.replace(code,"            oldprimitive=f.primitive;oldeos=f.eos;gas=model.bulk.eos.gas",
        "            oldprimitive=f.primitive;oldeos=f.eos;gas=model.bulk.eos.gas\n            model.face_K=np.interp(model.edge,model.bulk.d['r'],np.array([deep[j]['K'][k] for j in range(m.nb)]))\n            acoustic['max_deep_K_change']=max(acoustic['max_deep_K_change'],float(np.max(abs(model.face_K/baselineK-1))))")
    code=prior.replace(code,'            finally:f.primitive=oldprimitive;f.eos=oldeos;model.bulk.eos.gas=gas',
        '            finally:f.primitive=oldprimitive;f.eos=oldeos;model.bulk.eos.gas=gas;model.face_K=baselineK.copy()')
    code=prior.replace(code,"            dF=F.astype(np.longdouble)-before[0].astype(np.longdouble);dG=G.astype(np.longdouble)-before[1].astype(np.longdouble)",
        "            previous=np.load(OLDHYDRO/f'point-{k}.npz')\n            dF=F.astype(np.longdouble)-before[0].astype(np.longdouble)-previous['flux']\n            dG=G.astype(np.longdouble)-before[1].astype(np.longdouble)-previous['gravity']\n            dt=min(dt,.35*np.min(np.diff(model.m.rf)/(C*model.m.a/model.m.B*(np.max(abs(prim[1]))+acoustic['max_speed']+1e-100))))")
    code=prior.replace(code,"        write(HYDRO/f'{label}.json',result);return result",
        "        result['acoustic']=acoustic\n        write(HYDRO/f'{label}.json',result);return result")
    ns=dict(vars(prior),HYDRO=HYDRO,OLDHYDRO=prior.HYDRO,lookup=lookup,deep=deep,
            FACE_INPUTS=list(prior.HYDRO.glob('*-native-states.json')))
    write(OUT/'force-plan.json',dict(classification='Counterexample candidate',seconds=120,
        derivative_order='Match deep states by exact native coordinates and then canonical time, never production serialization order.',
        owner='Restore inherited face_K after each knot; only the new native force evaluation gets the replacement.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'native.json',OUT/'registered-native-producer.py']}))
    exec(compile(code,__file__,'exec'),ns);(OUT/'expanded-native-force.py').write_text(code)
    result=ns['hydro_sources'](False,list(range(17)),'force')
    assert result['new_native_states']==0 and result['native_calls']==0
    result.update(seconds=time.monotonic()-start,complete=True)
    write(OUT/'force.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    signal.setitimer(signal.ITIMER_REAL,0.)


def pilot():
    assert not (OUT/'pilot.json').exists();deadline(100);begin=time.monotonic();prior.initialize()
    model=prior.current.Capture();data=model.bulk.eos.base.d
    high=cached.CachedNative(12000);high.prefix=dict(high.prefix)
    low=prior.current.support.fast_native(12000);low.prefix=dict(low.prefix)
    rows=read(prior.EOS/'production-samples.json');deep=[r for r in rows if r['kind']=='deep']
    at=[r for r in rows if r['kind']=='atmosphere'];samples=[]
    choices=[next(r for r in deep if r['it']==it and r['cell']==j) for it,j in [(0,0),(8,9),(16,18)]]
    for it in [0,8,16]:
        group=[r for r in at if r['it']==it]
        choices += [group[i] for i in sorted({0,len(group)//2,len(group)-1})]
    try:
        for r in choices:
            if r['kind']=='deep':
                cached.setup(high,data,r['cell']);n=high;key=(r['x'],r['lt'],r['y']);physical=r['pressure']
            else:
                # A changed deep anchor must never leak into atmospheric calls.
                if not samples or samples[-1]['kind']=='deep':
                    high=cached.CachedNative(12000);high.prefix=dict(high.prefix)
                n=low if r['rho']<model.flow.eos.original.floor else high
                key=(float(np.log(r['rho'])),r['lt'],r['y']);physical=None
            value=derivative(n,key,np.array(r['raw']),physical)
            samples.append(dict(value,kind=r['kind'],it=r['it'],cell=r['cell']))
            write(OUT/'pilot-samples.json',samples)
        counts=dict(atmosphere=6392,deep=19*17)
        costs={kind:max(v['seconds'] for v in samples if v['kind']==kind) for kind in counts}
        work=2*sum(costs[k]*counts[k] for k in counts);wall=work/3+30
        result=dict(classification='Counterexample candidate',passed=True,eligible=wall<1000 and work<2600,
                    estimated_worker_seconds=work,upper_wall_seconds=wall,costs=costs,counts=counts,
                    seconds=time.monotonic()-begin,calls=sum(r['calls'] for r in samples),
                    maximum_derivative_comparison=max(max(r['comparison']) for r in samples),
                    maximum_first_law=max(max(abs(np.array(r['first_law']))) for r in samples))
        write(OUT/'pilot.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'pilot-failure.json',dict(error=repr(exc),completed=len(samples),seconds=time.monotonic()-begin));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


if __name__=='__main__':globals()[sys.argv[1]]()
