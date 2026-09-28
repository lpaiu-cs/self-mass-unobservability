"""Independent primitive, mass-normalization and emission-pairing audit."""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_right,bisect_left
import ast,json,signal,time
import def_native_charge_null_infinity as run


def main():
    out=run.OUT;assert not (out/'audit.json').exists()
    paths=[Path(__file__),Path(run.__file__),out/'plan.json',out/'production.json',out/'charge.npz',
           out/'exterior-128-a8-r8.npz',run.APPLIED/'return/gr/combined-source.npz',
           run.ROOT/'def-native-global-scalar-closure/bound.json',run.ROOT/'def-native-global-scalar-closure/check.json']
    run.write(out/'audit-plan.json',dict(classification='Counterexample candidate',budget_seconds=60,
        claim='Check the emission primitives by direct sums, the normalized charge by high-precision subtraction, and all between-knot body-debit mismatches using exact rational piecewise histories.',
        scope='Exact rational audit of the declared stored step/packet/linear-port representation. It does not certify the continuous physical emission, EOS, discarded matter, coupled closure, or numerical quadrature.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():run.sha(p) for p in paths}))
    start=time.monotonic()
    def timeout(*_):raise TimeoutError('Phase145 independent audit cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    import numpy as np
    import mpmath as mp
    from types import SimpleNamespace
    LD=np.longdouble
    raw=run.read(out/'production.json');assert raw['passed']
    wave=dict(np.load(out/'charge.npz'));source=dict(np.load(run.APPLIED/'return/gr/combined-source.npz'))
    env=dict(vars(run),np=np,LD=LD,SimpleNamespace=SimpleNamespace,emission={},m=SimpleNamespace(T=float(wave['t'][-1])))
    node=next(n for n in ast.walk(ast.parse(Path(run.__file__).read_text())) if isinstance(n,ast.FunctionDef) and n.name=='inputs')
    exec(compile(ast.Module(body=[node],type_ignores=[]),'saved-owner-primitives','exec'),env)
    primitive=[];fine=None
    for n in [64,128]:
        with np.load(run.UPDATED/f'background/accepted-ports-{n}.npz') as z:L=z['angular_luminosity'];h=float(z['h'])
        with np.load(run.APPLIED/f'total-photons/steps-{n}-reference-128.npz') as z:t=z['accepted_angular_times'];delta=z['accepted_angular_luminosity']
        gamma=1-1/np.sqrt(2);w=np.tile([1-gamma,gamma],n).astype(LD)*LD(h);packets=w[:,None]*delta
        query=np.unique(np.r_[np.linspace(0,wave['t'][-1],33),t,np.nextafter(t,-np.inf),np.nextafter(t,np.inf)])
        query=np.clip(query,0,wave['t'][-1]);returned=env['inputs'](n).primitives(query)
        dt=np.maximum(query[:,None].astype(LD)-np.arange(n,dtype=LD)*LD(h),0)
        direct=[np.minimum(dt,LD(h))@L,(dt*dt-np.maximum(dt-LD(h),0)**2)/2@L,
                (query[:,None]>=t)@packets,np.maximum(query[:,None].astype(LD)-t,0)@packets]
        errors=[float(np.max(abs(a-b))/max(np.max(abs(b)),LD('1e-300'))) for a,b in zip(returned,direct)]
        assert max(errors)<1e-12,errors;primitive.append(dict(steps=n,relative=errors,queries=len(query)))
        if n==128:fine=(L,h,t,delta,gamma)
    def fraction(v):return F(*v.as_integer_ratio())
    L,h,t,delta,gamma=fine;hf=F(h);T=F(float(wave['t'][-1]));aw=[F(i,32) for i in [1,3,5,7]]
    lum=[sum((fraction(v)*w for v,w in zip(row,aw)),F(0)) for row in L]
    bp=[i*hf for i in range(len(L)+1)];bq=[F(0)]
    for value in lum:bq.append(bq[-1]+hf*value)
    pt=[F(float(v)) for v in t];pack=[];absolute=F(0)
    weights=[F(float(1-gamma)),F(float(gamma))]
    for i,row in enumerate(delta):
        pk=[hf*weights[i%2]*fraction(v)*w for v,w in zip(row,aw)]
        pack.append(sum(pk,F(0)));absolute+=sum(map(abs,pk),F(0))
    pq=[F(0)]
    for v in pack:pq.append(pq[-1]+v)
    st=[F(float(v)) for v in source['t']];sv=[fraction(v) for v in source['outer_cumulative_energy_erg']]
    knots=sorted(v for v in set(bp+pt+st+[F(0),T]) if 0<=v<=T)
    mismatches=[]
    for q in knots:
        j=max(0,min(bisect_right(bp,q)-1,len(lum)-1));bg=bq[j]+(q-bp[j])*lum[j]
        k=max(0,min(bisect_right(st,q)-1,len(st)-2));expected=sv[k]+(q-st[k])*(sv[k+1]-sv[k])/(st[k+1]-st[k])
        for idx in [bisect_left(pt,q),bisect_right(pt,q)]:mismatches.append(abs(bg+pq[idx]-expected))
    mismatch=max(mismatches);absolute+=bq[-1]
    assert absolute>0 and mismatch>=0
    mp.mp.dps=80;mp.iv.dps=40;I=mp.iv.mpf
    def B(v):
        q=v if isinstance(v,F) else fraction(v)
        return I(q.numerator)/I(q.denominator)
    def up(v):return float(np.nextafter(float(v.b),np.inf))
    prior=run.read(run.ROOT/'def-native-global-scalar-closure/bound.json')
    geom=run.read(run.ROOT/'def-native-global-scalar-closure/check.json')
    # Same finite source coefficients, with a new exact total-variation energy.
    C=B(29979245800.);G=B(6.6743e-8);M=B(float(source['M_cm']));K=abs(B(float(source['K_cm'])))
    r0=B(float(geom['outer_radius_cm']));rmin=I(prior['causal_radius_min_cm']);bounds=prior['frozen_polynomial_bounds']
    Nmin=I(bounds['lapse']['lower']);bmin=1-2*I(bounds['mass']['upper'])/rmin
    b0=1-2*M/r0;c0=Nmin*mp.iv.sqrt(b0);kap=M/(r0*b0)+K*K/(2*c0*c0*r0*r0)
    assert float(kap.b)<1
    energy=G/C**4*B(absolute)
    es=energy*K*mp.iv.pi/(4*c0*c0*r0*mp.iv.sqrt(1-kap));em=C*B(T)/2*energy*K/(c0*b0*r0*r0)
    integral=I(0)
    # The stored coefficient maxima include every cell. Keeping deeper cells
    # is conservative and avoids inferring the original face selection.
    for row,half in zip(wave['new_source_coefficient_maxima'],wave['source_optical_half']):integral+=2*B(float(half))*sum((B(float(v)) for v in row),I(0))
    norm=C*B(T)/2*integral;eta=I(prior['global_potential_contraction'])
    ke=I(prior['mass_source_coefficient_absolute_per_cm2']);length=(r0-rmin)/(Nmin*mp.iv.sqrt(bmin))
    pair=C*B(T)/2*length*ke*G/C**4/Nmin*B(mismatch)
    potential=eta/(1-eta)*(norm+es+em+pair)
    epsilon=B(float(wave['epsilon'][-1]));uncertainty=up((potential+pair)/M/(1-epsilon))
    value=float(wave['normalized'][-1]);interval=[float(np.nextafter(value-uncertainty,-np.inf)),float(np.nextafter(value+uncertainty,np.inf))]
    alpha=mp.mpf(raw['alpha0']);eps=mp.mpf(float(wave['epsilon'][-1]));scalar=mp.mpf(float(wave['scalar'][-1]))
    direct=(alpha+scalar)/(1-eps)-alpha;identity=float(abs(direct-mp.mpf(value))/abs(direct));assert identity<1e-12 and interval[0]>0
    exterior=dict(np.load(out/'exterior-128-a8-r8.npz'))
    assert np.max(abs(exterior['normalized_exterior_parts'].sum(1)-wave['exterior']))<1e-45
    increment_arrival=float(exterior['arrived_parts_erg'][-1,1])
    result=dict(classification='Counterexample candidate',passed=True,primitive_checks=primitive,
        normalization_identity_relative=identity,exact_rational_pairing_knots=len(knots),
        exact_discrete_maximum_port_mismatch_erg=up(B(mismatch)),absolute_emission_energy_upper_erg=up(B(absolute)),
        source_pairing_mass_effect_bound=up(pair/M),all_orders_potential_bound=up(potential/M),
        conditional_stored_source_interval=interval,conditional_bound=uncertainty,
        mass_normalization_fraction=raw['endpoint_mass_normalization_term']/raw['endpoint_normalized'],
        native_increment_arrived_energy_erg=increment_arrival,
        native_increment_mass_normalization=raw['alpha0']*6.6743e-8/29979245800.**4*increment_arrival/float(source['M_cm']),
        source_pairing_scope='Every endpoint and both sides of each signed packet in the exact rational step/packet/linear-port representation; only its mass-constraint consequence is bounded.',
        interval_scope='Same frozen first-variation source and coefficient model, absolute signed emission and all potential repetitions. Not a complete physical error interval. All source cells retained conservatively.',
        missing='EOS/derivative reconstruction, continuous coupling, emission/charge quadrature and time error, missing/discarded matter feedback, nonlinear GR and static comparator/observability.',
        static_mass_redefinition_is_not_dynamic_novelty=True,uniform_EOS_derivative_bound=False,full_source_error_enclosed=False,
        coupled_fixed_point_verified=False,nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False,
        seconds=time.monotonic()-start)
    run.write(out/'audit.json',result);print(json.dumps(result),flush=True);signal.alarm(0)


if __name__=='__main__':
    try:main()
    except Exception as exc:run.write(run.OUT/'audit-failure.json',dict(error=repr(exc)));raise
