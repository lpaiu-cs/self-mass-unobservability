"""Native atmospheric inverse at the actual updated conserved endpoints.

No time history is fabricated from endpoints. Preserve S and Killing energy,
including the pressure-induced change of velocity in the scalar trace.
"""
from pathlib import Path
import json,signal,sys,time
import numpy as np
import sympy as s
import def_native_refined_thermochemistry as run

OUT=run.flow.OUT.parent/'def-native-atmosphere-inverse'
C=run.C;write=run.write;sha=run.sha


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(run.__file__),run.EV/'final-64.npz',run.EV/'final-128.npz',run.GR/'result.json',
        run.flow.INPUT/'balanced-20.npz',run.chem.old.old.cold.CACHE/'stable/gas.so',run.chem.old.LEVELS]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='51f42a5f4',
        claim='At the actual updated initial state and both final states, invert the native atmosphere EOS at fixed conserved D,S,K,H and insert its pressure/velocity into the exact trace and radial metric stress.',
        decision='A0.2percent constitutive failure or2percent atmospheric endpoint source change requires repair of the EOS owner before another coupled path. Small endpoint corrections remain limited evidence, not a final-charge or whole-trajectory certificate.',
        reuse='All128 initial and223+223 final active atmospheric cells saved by Phase130; no fluid replay or enlarged mesh. Native ion/molecule constraints and H definition remain those of the existing atmospheric owner.',
        gates=dict(conservative_energy=1e-12,conservative_momentum=1e-12,density=1e-10,root_logT_shift=.02,
            root_reseed=1e-10,constitutive=.002,endpoint_source=.02),
        budget=dict(pilot_seconds=40,pilot_native_calls=1000,production_seconds=180,production_native_calls=18000,
            CPU_threads=1,memory_GB=3,new_fluid_steps=0),
        forecast='Six actual cells span initial/deep atmosphere, pressure fan and dilute front. Reuse pilot states; require twice slowest measured root time for remaining states plus10s to fit180s. Do not extend the budget or support automatically.',
        stop='Preserve failed native root, constitutive or source gate. Endpoints do not determine the intervening retarded charge; do not interpolate an invented correction history or launch a full replay for diagnostics.',
        bindings={str(p):sha(p) for p in files}))
    S,Q,p,dp,k=s.symbols('S Q p dp k')
    dv=s.factor(S/(Q+p+dp)-S/(Q+p))
    change=s.factor(-S*dv-k*dp)
    assert s.simplify(change-(S*S/((Q+p)*(Q+p+dp))-k)*dp)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        premise='Fixed D,S,tau,H and lapse; Q=cx*D+tau, v=S/(Q+p).',
        trace_change='(v_old*v_new-3)*delta_p',radial_stress_change='(v_old*v_new-1)*delta_p',
        total_energy_change=0,baryon_change=0,scope='Exact conservative perfect-fluid identity; not a time-history or GR feedback bound.'))


def inputs(model):
    f=model.flow;m=model.m;rows={}
    for steps in [0,64,128]:
        U=f.initial.copy() if steps==0 else np.load(run.EV/f'final-{steps}.npz')['U']
        rho,v,lt,y=f.primitive(U);p,u,_,_,_,cv,_=f.eos.evaluate(rho,lt)
        tau=(U[2]-(m.a-m.a0)*f.eos.cx*U[0])/m.a
        rows[steps]=dict(U=U,rho=rho,v=v,lt=lt,y=y,p=p,u=u,cv=cv,tau=tau,active=np.flatnonzero(U[0]>=f.eos.floor))
    return rows


def root(native,f,data,j,offset=0.):
    D,S=float(data['U'][0,j]),float(data['U'][1,j]);tau=float(data['tau'][j]);y=float(data['y'][j])
    lt=float(data['lt'][j]+offset);p=float(data['p'][j]);start=native.ion.calls
    for iteration in range(8):
        for momentum_iteration in range(6):
            v=S/(f.eos.cx*D+tau+p);assert abs(v)<1
            r=np.sqrt(1-v*v);rho=D*r;state=native.state(float(np.log(rho)),lt,y)
            newp=float(state['raw'][1]/(f.eos.rho0*C*C))
            errp=abs(newp/p-1);p=newp
            if errp<1e-12:break
        else:raise AssertionError(('Native momentum inverse',j,errp))
        u=float(state['raw'][2]/C**2);wm=v*v/(r*(1+r))
        energy=f.eos.cx*D*wm+(rho*u+p)/(r*r)-p
        error=(energy-tau)/max(abs(tau),1e-300)
        if abs(error)<1e-12:break
        lt-=(energy-tau)/(rho*float(data['cv'][j])/(r*r))
        assert abs(lt-data['lt'][j])<.02,('Native temperature displacement',j,lt)
    else:raise AssertionError(('Native energy inverse',j,error))
    momentum=abs(v*(f.eos.cx*D+tau+p)-S)/max(abs(S),1e-300)
    density=abs(state['raw'][0]/(rho*f.eos.rho0)-1)
    assert momentum<1e-12 and density<1e-10
    return dict(cell=j,lt=lt,rho=rho,v=v,p=p,u=u,energy_relative=error,momentum_relative=momentum,
        density_relative=density,pressure_relative=p/float(data['p'][j])-1,
        old_pressure=float(data['p'][j]),old_velocity=float(data['v'][j]),old_lt=float(data['lt'][j]),
        y=y,raw=state['raw'].tolist(),fraction=state['fraction'].tolist(),affinity=float(state['affinity']),
        population_error=float(state['population_error']),native_calls=native.ion.calls-start,iterations=iteration+1)


def sample(pilot=False):
    name='pilot' if pilot else 'production';assert not (OUT/f'{name}.json').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    if not pilot:assert json.loads((OUT/'pilot.json').read_text())['eligible']
    signal.signal(signal.SIGALRM,run.flow.old.optical.timeout);signal.alarm(40 if pilot else 180);start=time.monotonic()
    model=run.Coupled();f=model.flow;data=inputs(model);native=run.chem.old.Native(cap=1000 if pilot else 18000)
    assert abs(np.exp(native.lr)/f.eos.rho0-1)<1e-10
    rows=[] if pilot else json.loads((OUT/'pilot-samples.json').read_text());done={(r['steps'],r['cell']) for r in rows}
    choices=[(0,0),(0,127),(128,0),(128,128),(128,222),(64,222)] if pilot else [(n,int(j)) for n,d in data.items() for j in d['active']]
    durations=[];reseed=[]
    try:
        for n,j in choices:
            if (n,j) in done:continue
            tick=time.monotonic();r=root(native,f,data[n],j);durations.append(time.monotonic()-tick)
            r['steps']=n;rows.append(r);done.add((n,j))
            if pilot:
                other=root(native,f,data[n],j,.001);diff=abs(other['p']/r['p']-1);assert diff<1e-10
                reseed.append(dict(steps=n,cell=j,pressure_relative=diff))
        write(OUT/f'{name}-samples.json',rows)
    except Exception as exc:
        write(OUT/f'{name}-partial.json',rows);write(OUT/f'{name}-failure.json',dict(error=repr(exc),states=len(rows),native_calls=native.ion.calls,seconds=time.monotonic()-start));raise
    remaining=sum(len(d['active']) for d in data.values())-len(rows);forecast=2*max(durations)*remaining+10
    by={(r['steps'],r['cell']):r for r in rows};source=[]
    if not pilot:
        vol=4*np.pi*model.m.RJ**2*model.m.vol;unit=f.eos.rho0*C*C*vol
        init=data[0];initial_trace=init['tau']-init['U'][1]*init['v']-3*init['p']
        for n in [0,64,128]:
            d=data[n];active=d['active'];dp=np.zeros(f.n);dv=np.zeros(f.n)
            for j in active:dp[j]=by[n,int(j)]['p']-d['p'][j];dv[j]=by[n,int(j)]['v']-d['v'][j]
            trace=(-d['U'][1]*dv-3*dp)*unit;stress=(-d['U'][1]*dv-dp)*unit
            identity=np.array([(by[n,int(j)]['v']*d['v'][j]-3)*dp[j]*unit[j] for j in active])
            err=float(np.max(abs(trace[active]-identity))/max(np.max(abs(trace)),1e-300));assert err<1e-8
            # Retain the original initial pressure reference. No subtraction
            # of the native initial offset from the corrected source.
            baseline=(d['tau']-d['U'][1]*d['v']-3*d['p']-initial_trace)*unit
            scale=max(np.sum(abs(baseline)),1.)
            source.append(dict(steps=n,trace_correction_L1_erg=float(np.sum(abs(trace))),trace_correction_sum_erg=float(np.sum(trace)),
                old_response_L1_erg=float(np.sum(abs(baseline))),relative_to_response=None if n==0 else float(np.sum(abs(trace))/scale),identity_relative=err))
            np.savez_compressed(OUT/f'source-{n}.npz',U=d['U'],rho=d['rho'],v=d['v'],lt=d['lt'],y=d['y'],old_p=d['p'],delta_p=dp,delta_v=dv,trace_correction_erg=trace,stress_correction_erg=stress,old_trace_response_erg=baseline)
    worst=max(abs(r['pressure_relative']) for r in rows)
    source_ok=all(r['relative_to_response'] is None or r['relative_to_response']<.02 for r in source)
    result=dict(classification='Counterexample candidate',passed=bool(worst<.002 and source_ok),states=len(rows),
        maximum_pressure_relative=worst,maximum_energy_relative=max(abs(r['energy_relative']) for r in rows),
        maximum_momentum_relative=max(r['momentum_relative'] for r in rows),maximum_density_relative=max(r['density_relative'] for r in rows),
        maximum_logT_change=max(abs(r['lt']-r['old_lt']) for r in rows),native_calls=native.ion.calls,seconds=time.monotonic()-start,
        forecast_upper_seconds=forecast,eligible=bool(forecast<180 and worst<.002),reseed=reseed,source=source,
        native_endpoint_inverse_complete=not pilot,full_atmosphere_history_audited=False,uniform_EOS_error_bound=False,
        actual_atmosphere_EOS_re_evolution=False,final_charge_solved=False)
    write(OUT/f'{name}.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1]
    if action=='prepare':prepare()
    elif action=='pilot':sample(True)
    elif action=='production':sample()
    else:raise ValueError(action)
