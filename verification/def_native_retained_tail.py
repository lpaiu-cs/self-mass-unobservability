"""Extend the same native EOS into the formerly projected atmosphere tail.

Counterexample candidate: one declared support change, original clocks and
physical horizon. No cold clamp, ideal-gas replacement or automatic refinement.
"""
from pathlib import Path
import json,sys,time,signal
import def_native_charge_null_infinity as saved

OUT=saved.ROOT/'native-retained-tail';BANK=OUT/'support';EV=OUT/'evolution';GR=OUT/'gr'
COLD=saved.ROOT/'def-native-cold-coupling'
read,write,sha=saved.read,saved.write,saved.sha


def prepare():
    assert not OUT.exists()
    for p in [OUT,BANK,EV,GR]:p.mkdir(exist_ok=True)
    paths=[Path(__file__),Path('verification/def_native_cold_coupling.py'),Path('verification/def_native_hydrogen_exchange.py'),
        Path('verification/def_native_material_join.py'),Path('verification/def_native_refined_thermochemistry.py'),
        Path('verification/def_native_updated_gr_return.py'),COLD/'repaired-bank.npz',COLD/'spectrum-bank.npz',
        saved.THERMAL/'bank.npz',saved.THERMAL/'evolution/result.json',saved.ROOT/'native-discarded-radiation/audit.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='7178d21af',
        previous_goal_turn='Progress: passive radiative floor uncertainty quantified; actual material exchange remained missing and production budget failure was preserved.',
        claim='Retain the formerly deleted atmosphere with the same native EOS, actual pressure/momentum/species and paired photon transfers, then compare its same-solution GR charge with the accepted trajectory.',
        decision='Does applying native low-density support to actual coupled64/128 evolution change the scalar residue or expose a physical support/energy failure? Passive-tail bounds alone do not answer this.',
        support=dict(log_density=[-24.,-22.,-20.,-18.],added_temperature_K=[1.,3.,10.,30.],
            old_minimum_log_density=-18.,old_minimum_temperature_K=80.,neutral_planes=[1e-16,.001]),
        reuse='Reuse62 old boundary EOS/optical states and all old constitutive owners above the old density floor. Preserve original531 cells,8 angles,152 frequencies,3.434431ms and64/128 clocks; resume measured new prefixes, never replay an accepted new prefix.',
        changes='Lower the numerical projection threshold only within the registered native support. This is an explicit physical-input test, not permission to extrapolate. Record retained sub-old-floor material, actual remaining discard, accepted angular ports, all-step GR moments and17 full states in this run.',
        budgets=dict(probe_seconds=45,probe_native_calls=160,bank_seconds=180,bank_native_calls=2000,
            controls_seconds=60,controls_native_calls=220,pilot_seconds=70,production_seconds=650,GR_seconds=90,
            CPU_threads=1,virtual_GiB=3),
        enforcement='Use an external timeout process around scientific actions in addition to an absolute monotonic deadline. SHA checks and WSL startup are reported separately; constructor and output costs count in action timers.',
        forecast='Probe the eight native cold/hot density corners first. Require2x measured native per-state forecast plus20s within180s. For evolution compare measured prefix physics/I-O costs with the prior353s full run and require1.5x the larger estimate plus20s below650s. Late CFL and low-density root costs are unmeasured.',
        gates=dict(native_population=1e-12,relative_H=1e-8,constitutive=.002,spectrum=.002,primitive=2e-14,
            energy=1e-8,port=1e-10,time=.02,quadrature=.002),
        stop='One support change only. Stop on unsupported T/rho/y, failed native or coupled gate, external timeout or forecast. No automatic further support, floor, grid, clock or horizon change; do not weaken criteria or rename partial success as full closure.',
        limitations='A successful lower-floor comparison is not a uniform continuum EOS derivative, complete floor-limit or nonlinear-GR certificate. Static comparison and observational completion requirements remain intact.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths}))


def deadline(start,cap):
    remaining=cap-(time.monotonic()-start)
    if remaining<=0:raise TimeoutError('Phase147 action deadline')
    def timeout(*_):raise TimeoutError('Phase147 action deadline')
    signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,remaining)


def probe():
    assert not (OUT/'probe.json').exists() and not (OUT/'probe-failure.json').exists()
    plan=read(OUT/'plan.json')
    for p,h in plan['bindings'].items():assert sha(Path(p))==h,p
    start=time.monotonic();deadline(start,plan['budgets']['probe_seconds'])
    import resource
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    import numpy as np
    import def_native_cold_coupling as cold
    n=cold.logarithmic_native(plan['budgets']['probe_native_calls']);deadline(start,45)
    initial=time.monotonic()-start;rows=[];states=[]
    try:
        for x,T in [(-24.,1.),(-24.,20000.),(-18.,1.),(-18.,80.)]:
            for y in [1e-16,.001]:
                before=time.monotonic();s=n.state(x,float(np.log(T)),y)
                assert np.all(np.isfinite(s['raw'])) and s['raw'][1]>0 and np.all(np.isfinite(s['log_fraction']))
                assert s['population_error']<1e-12
                states.append(s);rows.append(dict(x=x,T=T,y=y,seconds=time.monotonic()-before,
                    population_error=s['population_error'],pressure=float(s['raw'][1]),energy=float(s['raw'][2])))
                deadline(start,45)
        # Native and optical outputs are reused verbatim by the table action.
        np.savez_compressed(OUT/'probe-states.npz',raw=[s['raw'] for s in states],log_fraction=[s['log_fraction'] for s in states],
            affinity=[s['affinity'] for s in states],x=[s['x'] for s in states],lt=[s['lt'] for s in states],y=[s['y'] for s in states])
        per_state=sum(r['seconds'] for r in rows)/len(rows);remaining=218-6
        forecast=initial+per_state*remaining+10;upper=2*forecast+20;deadline(start,45)
        result=dict(classification='Counterexample candidate',passed=True,rows=rows,native_calls=n.ion.calls+n.variant_initial_calls,
            seconds=time.monotonic()-start,setup_seconds=initial,remaining_new_states=remaining,
            measured_bank_forecast_seconds=forecast,bank_upper_seconds=upper,bank_eligible=upper<180,
            actual_coupled_retention_applied=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/'probe.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        if states:np.savez_compressed(OUT/'partial-probe.npz',raw=[s['raw'] for s in states],x=[s['x'] for s in states],lt=[s['lt'] for s in states],y=[s['y'] for s in states])
        write(OUT/'probe-failure.json',dict(classification='Counterexample candidate',error=repr(exc),rows=rows,
            failed_state=dict(x=x,T=T,y=y),native_calls=n.ion.calls+n.variant_initial_calls,seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


if __name__=='__main__':globals()[sys.argv[1]]()
