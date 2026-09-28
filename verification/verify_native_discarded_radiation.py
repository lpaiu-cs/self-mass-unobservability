"""Independent compact-data audit; retain Phase146's production budget failure."""
from pathlib import Path
from fractions import Fraction as Q
from types import SimpleNamespace
import ast,json,sys,time,signal
import def_native_discarded_radiation as owner

OUT=owner.OUT;read,write,sha=owner.read,owner.write,owner.sha


def prepare():
    assert not (OUT/'audit-plan.json').exists()
    paths=[Path(__file__),Path(owner.__file__),OUT/'executed-producer.py',OUT/'plan.json',OUT/'result.json',
        owner.prior.APPLIED/'return/gr/combined-source.npz',owner.prior.OUT/'production.json',
        owner.prior.ROOT/'def-native-global-scalar-closure/bound.json']
    for n in [64,128]:paths += [OUT/f'inputs-{n}.npz',owner.prior.APPLIED/f'return/steps-{n}-reference-128.npz']
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',budget_seconds=20,
        scope='Compact ledger verification only; no new physical trajectory or photon-history replay. Exact rational dose audit, add signed finite-return discard envelopes, use the conservative global lapse in exterior coefficients, test the monotonic deadline repair.',
        acceptance='Preserve the measured111.666857171s production as a failure of its90s cap. This audit does not retroactively enlarge that cap or turn the phase into an overall pass.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths}))


def run():
    assert not (OUT/'audit.json').exists()
    plan=read(OUT/'audit-plan.json')
    for p,h in plan['bindings'].items():assert sha(Path(p))==h,p
    original=read(OUT/'plan.json');assert sha(OUT/'executed-producer.py')==original['bindings']['verification/def_native_discarded_radiation.py']
    started=time.monotonic()
    def timeout(*_):raise TimeoutError('Phase146 compact audit20s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(plan['budget_seconds'])
    import numpy as np
    import mpmath as mp
    iv=mp.iv;iv.dps=40
    def B(x):
        n,d=float(x).as_integer_ratio();return iv.mpf(n)/iv.mpf(d)
    def I(x):return iv.mpf(x.numerator)/iv.mpf(x.denominator)
    def up(x):return float(np.nextafter(float(x.b),np.inf))
    prod=read(OUT/'result.json');last=read(owner.prior.OUT/'production.json')
    geo=read(owner.prior.ROOT/'def-native-global-scalar-closure/bound.json');bb=geo['frozen_polynomial_bounds']
    with np.load(owner.prior.APPLIED/'return/gr/combined-source.npz') as d:
        M=B(float(d['M_cm']));K=B(abs(float(d['K_cm'])));cx=B(float(d['cx']));r0=B(float(d['edges'][-1]));T=B(float(d['t'][-1]))
    C=B(29979245800.);G=B(6.6743e-8);phi=B(bb['phi']['absolute']);N=B(bb['lapse']['lower']);amin=N*iv.exp(-2*phi**2)
    rmin=B(geo['causal_radius_min_cm']);Phi=B(bb['Phi']['absolute']);eta=B(geo['global_potential_contraction'])
    c0=amin*iv.sqrt(1-2*M/r0);b0=1-2*M/r0;rk=M/(r0*b0)+K*K/(2*c0*c0*r0*r0)
    opacity=B(prod['opacity_bound_cm2_g']);alpha=B(last['alpha0']);eps0=B(last['epsilon']);s=B(abs(last['endpoint_scalar']))
    rows=[]
    for row in prod['paths']:
        n=row['steps'];z=np.load(OUT/f'inputs-{n}.npz');tt=z['original_t'];mass=z['baryon_discard_g'];ucap=z['interval_radiation_cap']
        # Separate exact sum in the stored binary values, not the producer's
        # interval accumulation or its transient Python objects.
        dose=sum((Q(float(mass[j+1]))*Q(float(u))*Q(float(tt[j+1]-tt[j])) for j,u in enumerate(ucap)),Q(0))
        A=2*C*opacity*I(dose)
        rel=abs(up(A)/row['incident_background_energy_bound_erg']-1);assert rel<1e-14
        with np.load(owner.prior.APPLIED/f'return/steps-{n}-reference-128.npz') as d:
            assert np.max(abs(d['t']-tt))<1e-18
            signed=d['discards_scaled'].astype(np.longdouble)*np.longdouble(float(d['amplitude']))
        dm=B(float(np.max(abs(signed[:,0]))))
        df=B(float(np.max(abs(signed[:,2]))))+cx*C*C*dm  # a0<=1, no cancellation assumed
        da=2*C*opacity*dm*T*B(float(max(z['radiation_cap'])))
        F=B(row['discarded_Killing_erg'])+df;A+=da;E=F+A;P=F+2*A
        direct=G/(2*C**3)*T*4*phi/rmin*2*E/amin/M
        stress=G/(2*C**3)*T*Phi*(E+P)/amin/M
        ext=(G/C**4*P*K*iv.pi/(4*c0*c0*r0*iv.sqrt(1-rk))+C*T/2*G/C**4*P*K/(c0*b0*r0*r0))/M
        length=(r0-rmin)/(N*iv.sqrt(1-2*B(bb['mass']['upper'])/rmin))
        constraint=C*T/2*length*B(geo['mass_source_coefficient_absolute_per_cm2'])*G/C**4/amin*(E+P)/M
        fixed=(direct+stress+ext+constraint)/(1-eta);de=G/C**4*P/M
        norm=fixed/(1-eps0-de)+(alpha+s)*de/((1-eps0-de)*(1-eps0))
        rows.append(dict(steps=n,exact_rational_dose_relative=rel,
            native_discard_mass_absolute_g=up(dm),native_discard_Killing_energy_absolute_erg=up(df),
            native_discard_additional_dose_bound_erg=up(da),
            fixed_mass_scalar_bound=up(fixed),normalized_charge_bound=up(norm),
            fixed_scalar_fraction=up(fixed)/float(s.a),normalized_over_scalar=up(norm)/float(s.a),
            passive_normalized_interval=[float(np.nextafter(last['endpoint_normalized']-up(norm),-np.inf)),
                                         float(np.nextafter(last['endpoint_normalized']+up(norm),np.inf))]))
    # Regression: a cancelled/overwritten alarm cannot remove the monotonic
    # budget check, and rearming never resets the original deadline.
    tree=ast.parse(Path(owner.__file__).read_text());run=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
    nodes=[x for x in run.body if isinstance(x,ast.FunctionDef) and x.name in ['timeout','deadline']]
    ticks=[80.];calls=[]
    scope=dict(started=0.,plan={'budget':{'seconds':90}},time=SimpleNamespace(monotonic=lambda:ticks[0]),
        signal=SimpleNamespace(SIGALRM=14,ITIMER_REAL=0,signal=lambda *a:calls.append(('handler',a)),setitimer=lambda *a:calls.append(('timer',a))))
    exec(compile(ast.Module(body=nodes,type_ignores=[]),'<deadline-check>','exec'),scope)
    scope['deadline']();assert calls[-1][1][1]==10.
    calls.clear();ticks[0]=89.;scope['deadline']();assert calls[-1][1][1]==1.
    ticks[0]=91.
    try:scope['deadline']();raise AssertionError('Over-budget execution accepted')
    except TimeoutError:pass
    budget=prod['seconds']<=original['budget']['seconds'];assert not budget
    elapsed=time.monotonic()-started;assert elapsed<plan['budget_seconds']
    result=dict(classification='Counterexample candidate',passed=True,paths=rows,
        numerical_subproblem_passed=all(x['normalized_over_scalar']<.02 for x in rows),production_budget_passed=budget,
        overall_phase_passed=False,deadline_repair_tested=True,production_repeated=False,
        conservative_boundary_lapse_used=True,signed_native_discard_envelope_included=True,
        interpretation='This closes the passive radiative energy-growth subproblem of the declared compact reservoir. It does not verify confinement or the incident envelope in a full live coupled trajectory, or bound pressure-work and changed-photon amplification by retained matter. Preserve the resource failure.',
        full_floor_feedback_enclosed=False,full_source_error_enclosed=False,coupled_fixed_point_verified=False,
        nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,seconds=elapsed)
    write(OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
