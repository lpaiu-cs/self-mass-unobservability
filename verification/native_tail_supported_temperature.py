"""Apply lower-density native support without assuming unneeded cold states.

Cold corner failures are retained. Actual conservative inversion must remain
inside the original80K temperature support or stop; no temperature clamping.
"""
from pathlib import Path
from types import FunctionType,MethodType
import json,sys,time,signal
import def_native_retained_tail as owner

OUT=owner.OUT/'supported-temperature';BANK=OUT/'support';EV=OUT/'evolution';GR=OUT/'gr'
read,write,sha=owner.read,owner.write,owner.sha
FAILURES=[owner.OUT/'probe-failure.json',owner.OUT/'density-inverse/probe-failure.json',
    owner.OUT/'population-continuation/probe-failure.json',owner.OUT/'forced-population-continuation/probe-failure.json']


def prepare():
    assert not OUT.exists()
    for p in [OUT,BANK,EV,GR]:p.mkdir(exist_ok=True)
    original=read(owner.OUT/'plan.json');spent=sum(read(p)['seconds'] for p in FAILURES);calls=sum(read(p)['native_calls'] for p in FAILURES)
    paths=[Path(__file__),Path(owner.__file__),owner.OUT/'plan.json',*FAILURES,
        owner.COLD/'repaired-bank.npz',owner.COLD/'spectrum-bank.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='7178d21af',
        rejected='The originally proposed1K support has not passed. Native error104 is the electron-inventory underflow guard; direct/electron-variable and two continuation attempts failed. Those are not accepted support.',
        decision='The actual trajectory has not yet requested1K. Test the same complete lower-density64/128 coupled trajectories using independently checked native support at the original80K lower temperature. If a conservative state requests colder temperature, stop and solve that actual state; never clamp it or count this branch as a full1K table.',
        support=dict(x=[-24.,-22.,-20.,-18.],temperature_K=[80.,20000.],new_states=186,reused_boundary_states=62),
        budgets=dict(warm_probe_seconds=45-spent,warm_probe_native_calls=160-calls,bank_seconds=180,bank_native_calls=2000,
            controls_seconds=60,controls_native_calls=220,pilot_seconds=70,production_seconds=650,GR_seconds=90,CPU_threads=1,virtual_GiB=3),
        gates=original['gates'],
        budget_scope='Remaining original probe allocation, original unstarted table/control/evolution allocations; no larger total allocation. External timeout includes imports and SHA checks. No second resolution, horizon or further floor reduction.',
        stop='Stop at unsupported conserved state, native/control failure, clock gate or budget. Original cold-support failure remains false even if the actual warmer trajectory succeeds.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths}))


def warm_probe():
    assert not (OUT/'warm-probe.json').exists();plan=read(OUT/'plan.json');start=time.monotonic()
    owner.deadline(start,plan['budgets']['warm_probe_seconds'])
    import numpy as np
    import def_native_cold_coupling as cold
    n=cold.logarithmic_native(plan['budgets']['warm_probe_native_calls']);setup=time.monotonic()-start;rows=[];states=[]
    for T,y in [(80.,.001),(20000.,1e-16),(20000.,.001)]:
        begin=time.monotonic();s=n.state(-24.,float(np.log(T)),y);states.append(s)
        assert s['population_error']<1e-12 and np.all(np.isfinite(s['raw'])) and s['raw'][1]>0
        rows.append(dict(T=T,y=y,seconds=time.monotonic()-begin,population_error=s['population_error']))
        owner.deadline(start,plan['budgets']['warm_probe_seconds'])
    np.savez_compressed(OUT/'warm-probe.npz',raw=[s['raw'] for s in states],log_fraction=[s['log_fraction'] for s in states],affinity=[s['affinity'] for s in states],
        x=[s['x'] for s in states],lt=[s['lt'] for s in states],y=[s['y'] for s in states])
    estimate=setup+max(r['seconds'] for r in rows)*183+10;upper=2*estimate+20
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,native_calls=n.ion.calls+n.variant_initial_calls,
        seconds=time.monotonic()-start,bank_forecast_seconds=estimate,bank_upper_seconds=upper,eligible=upper<180,
        colder_support_accepted=False,actual_coupled_retention_applied=False)
    write(OUT/'warm-probe.json',result);signal.setitimer(signal.ITIMER_REAL,0.);print(json.dumps(result),flush=True)


def cached_bank():
    assert not (OUT/'cached-bank-plan.json').exists()
    write(OUT/'cached-bank-plan.json',dict(classification='Counterexample candidate',
        original_forecast_accepted=False,repair='Materialize the repeatedly decompressed native seed NPZ arrays once, preserving their dtype and values. Recheck the same three native states before measuring the remaining identical bank. No EOS, grid, gate or budget changes.',
        budget_seconds=180,budget_scope='Setup, cache controls, native table construction and output all share the original unstarted180s bank allocation.',
        native_call_cap=2000,gate='All three cached outputs agree bitwise with the saved native outputs; twice the new maximum per-state forecast plus20s must fit the remaining bank allocation.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'executed-probe.py',OUT/'warm-probe.json',OUT/'warm-probe.npz']}))
    start=time.monotonic();owner.deadline(start,180)
    import numpy as np
    import def_native_cold_coupling as cold
    n=cold.logarithmic_native(2000);n.prefix=dict(n.prefix)
    z=np.load(OUT/'warm-probe.npz');rows=[]
    for j in range(len(z['x'])):
        begin=time.monotonic();a=n.state(float(z['x'][j]),float(z['lt'][j]),float(z['y'][j]))
        assert all(np.array_equal(a[key],z[key][j]) for key in ['raw','log_fraction','affinity']), 'Cached native seed changed output'
        rows.append(time.monotonic()-begin)
    elapsed=time.monotonic()-start;upper=2*(max(rows)*183+10)+20
    result=dict(classification='Counterexample candidate',bitwise=True,seconds=elapsed,per_state_seconds=rows,
        remaining_upper_seconds=upper,remaining_seconds=180-elapsed,eligible=upper<180-elapsed)
    write(OUT/'cached-bank-pilot.json',result);print(json.dumps(result),flush=True)
    assert result['eligible'],'Cached bank forecast exceeds original remaining allocation'
    bank(n,start)


def reallocated_bank():
    assert not (OUT/'reallocated-bank-plan.json').exists()
    previous=read(OUT/'cached-bank-pilot.json');cap=260-previous['seconds']
    write(OUT/'reallocated-bank-plan.json',dict(classification='Counterexample candidate',
        original_forecasts_accepted=False,decision='Keep both rejected dispatch forecasts. Transfer80s from the unstarted90s GR readout to the bank, preserving the total action allocation1195s, all physical states and all scientific gates. Bank setup already spent39.095s is deducted. This targets the actual coupled material bottleneck; GR readout has only10s left in this phase.',
        bank_total_seconds=260,bank_already_spent_seconds=previous['seconds'],hard_seconds=cap,GR_seconds=10,
        total_original_seconds=1195,total_reallocated_seconds=1195,native_call_cap=1983,
        forecast='Use the already measured cached maximum per-state cost, doubled plus20s and10s output allowance. Add the actual new setup cost before dispatch; stop if it no longer fits.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'executed-cache-pilot.py',OUT/'cached-bank-pilot.json']}))
    start=time.monotonic();owner.deadline(start,cap)
    import def_native_cold_coupling as cold
    n=cold.logarithmic_native(1983);n.prefix=dict(n.prefix)
    setup=time.monotonic()-start;upper=setup+previous['remaining_upper_seconds']
    result=dict(eligible=upper<cap,forecast_upper_seconds=upper,hard_seconds=cap,setup_seconds=setup)
    write(OUT/'reallocated-bank-dispatch.json',result);print(json.dumps(result),flush=True)
    assert result['eligible'],'Reallocated bank forecast'
    bank(n,start,cap,'reallocated-bank-dispatch.json')


def bank(native=None,started=None,cap=180,dispatch=None):
    assert read(OUT/(dispatch or ('cached-bank-pilot.json' if native is not None else 'warm-probe.json')))['eligible'] and not (BANK/'result.json').exists()
    start=time.monotonic() if started is None else started;owner.deadline(start,cap)
    import numpy as np
    import def_native_cold_coupling as cold
    d=dict(np.load(owner.COLD/'repaired-bank.npz'));s=dict(np.load(owner.COLD/'spectrum-bank.npz'))
    xs=np.array([-24.,-22.,-20.,-18.]);lt=d['lt'];ys=d['ys'];raw=np.zeros((2,4,len(lt),21));lf=np.zeros((*raw.shape[:3],10));af=np.zeros(raw.shape[:3]);done=np.zeros(af.shape,bool)
    raw[:,-1]=d['raw'][:,0];lf[:,-1]=s['log_fraction'][:,0];af[:,-1]=s['affinity'][:,0];done[:,-1]=True
    z=np.load(OUT/'warm-probe.npz')
    for j in range(len(z['x'])):
        iy=int(np.argmin(abs(ys-z['y'][j])));it=int(np.argmin(abs(lt-z['lt'][j])))
        assert xs[0]==z['x'][j] and ys[iy]==z['y'][j] and abs(lt[it]-z['lt'][j])<1e-14
        raw[iy,0,it]=z['raw'][j];lf[iy,0,it]=z['log_fraction'][j];af[iy,0,it]=z['affinity'][j];done[iy,0,it]=True
    n=cold.logarithmic_native(2000) if native is None else native;initial_calls=n.ion.calls+n.variant_initial_calls
    try:
        for iy,y in enumerate(ys):
            for ix,x in enumerate(xs[:-1]):
                for it,t in enumerate(lt):
                    if done[iy,ix,it]:continue
                    a=n.state(float(x),float(t),float(y));assert a['population_error']<1e-12
                    raw[iy,ix,it]=a['raw'];lf[iy,ix,it]=a['log_fraction'];af[iy,ix,it]=a['affinity'];done[iy,ix,it]=True
                owner.deadline(start,cap)
                print(json.dumps(dict(bank_y=float(y),bank_x=float(x),done=int(done.sum()),seconds=time.monotonic()-start)),flush=True)
        values={k:d[k] for k in ['y0','rho0','cx','sunit','nH','s0']};values.update(x=xs,lt=lt,ys=ys,raw=raw,rates=np.ones((*raw.shape[:3],3,2)))
        np.savez_compressed(BANK/'repaired-bank.npz',**values)
        np.savez_compressed(BANK/'spectrum-bank.npz',x=xs,lt=lt,ys=ys,log_fraction=lf,affinity=af,binding=s['binding'],rho0=s['rho0'],nH=s['nH'])
        owner.deadline(start,cap)
        result=dict(classification='Counterexample candidate',passed=True,reused_old_boundary_states=62,reused_probe_states=3,new_states=183,
            native_calls=n.ion.calls+n.variant_initial_calls,initial_native_calls=initial_calls,seconds=time.monotonic()-start,
            prescribed_bath_rates_used=False,colder_support_accepted=False)
        write(BANK/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        np.savez_compressed(BANK/'partial.npz',raw=raw,log_fraction=lf,affinity=af,done=done,x=xs,lt=lt,ys=ys)
        write(BANK/'failure.json',dict(error=repr(exc),completed=int(done.sum()),native_calls=n.ion.calls+n.variant_initial_calls,seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def fast_native(cap):
    # Only rho,T,X,cx of Fan are consumed by this native state owner. Reading
    # their saved values avoids constructing unused stellar/opacity models.
    from types import SimpleNamespace
    import numpy as np
    import def_native_cold_coupling as cold
    module=cold.COLD.old.native;old=module.Fan
    env=np.load(module.old.prior.OUT/'final-envelope.npz');core=np.load(module.old.prior.OUT/'final-core.npz')
    d=np.load(owner.COLD/'repaired-bank.npz')
    fan=SimpleNamespace(rho=float(env['rho'][-1]),T=float(env['T'][-1]),X=core['X'][0],cx=float(d['cx']))
    assert fan.rho==float(d['rho0'])
    module.Fan=lambda *args,**kwargs:fan
    try:n=cold.logarithmic_native(cap)
    finally:module.Fan=old
    n.prefix=dict(n.prefix);return n


def reused_constructor_bank():
    assert not (OUT/'constructor-reuse-plan.json').exists()
    spent=read(OUT/'cached-bank-pilot.json')['seconds']+read(OUT/'reallocated-bank-dispatch.json')['setup_seconds']
    cap=260-spent
    write(OUT/'constructor-reuse-plan.json',dict(classification='Counterexample candidate',
        repair='Native Ions constructed Fan twice, which loads a complete unused stellar envelope and opacity model. Read only the exact saved rho,T,X,cx actually consumed by Native.state. Restore the Fan class immediately after construction. Keep the actual native library and all constraint/optical equations.',
        validation='Require bitwise equality at the three saved native corners; remeasure their state cost and keep the2x maximum forecast plus40s output/margin within the same remaining bank allocation.',
        hard_seconds=cap,bank_total_seconds=260,already_spent_seconds=spent,GR_seconds=10,native_call_cap=1981,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'executed-reallocation.py',OUT/'cached-bank-pilot.json',OUT/'reallocated-bank-dispatch.json',OUT/'warm-probe.npz']}))
    start=time.monotonic();owner.deadline(start,cap)
    import numpy as np
    n=fast_native(1981);setup=time.monotonic()-start;z=np.load(OUT/'warm-probe.npz');cost=[]
    for j in range(len(z['x'])):
        begin=time.monotonic();a=n.state(float(z['x'][j]),float(z['lt'][j]),float(z['y'][j]))
        assert all(np.array_equal(a[key],z[key][j]) for key in ['raw','log_fraction','affinity'])
        cost.append(time.monotonic()-begin)
    elapsed=time.monotonic()-start;upper=elapsed+2*max(cost)*183+40
    result=dict(classification='Counterexample candidate',bitwise=True,setup_seconds=setup,seconds=elapsed,per_state_seconds=cost,
        forecast_upper_seconds=upper,hard_seconds=cap,eligible=upper<cap)
    write(OUT/'constructor-reuse-dispatch.json',result);print(json.dumps(result),flush=True)
    assert result['eligible'],'Constructor reuse budget'
    bank(n,start,cap,'constructor-reuse-dispatch.json')


if __name__=='__main__':globals()[sys.argv[1]]()
