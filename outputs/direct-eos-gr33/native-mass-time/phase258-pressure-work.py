"""Test the known metric-pressure work primitive on the same saved clocks.

This locates a possible forcing quadrature defect. It does not patch a charge
or accept an unevolved correction as a physical solution.
"""
from pathlib import Path
import json,resource,time
import numpy as np
import sympy as s
import finish_returned_material_accuracy as producer

OUT=Path('native-mass-time258-work');OLD=producer.OUT
read,write,sha,bind=producer.read,producer.write,producer.sha,producer.bind
LD=np.longdouble
assert read(OUT/'canonical-affine-result.json')['passed']
assert not (OUT/'pressure-work.json').exists()
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2)
producer.joint.previous.original.inf.incident.native.deadline(900)
start=time.monotonic();error=None
write(OUT/'pressure-plan.json',dict(classification='Conjectural',
    claim='Test whether a right-sided derivative at closing Radau endpoints causes the native pressure-work clock difference.',
    method='Reproduce the applied mass-constraint derivative from the SAME original source polynomials, then evaluate their left derivative at knots. All source values and scalar fields remain identical. Compare its pressure work with direct work and integration by parts on the SAME stage/end values; jumps remain separate.',
    decision='Only a dominant, independently checked derivative-quadrature defect warrants a primitive-based repair of actual coupled stages. No diagnostic correction enters final charge.',
    budget_seconds=900,virtual_GiB=8,physical_steps=0,new_clock_paths=0,
    forecast='One original constructor measured about16seconds;17background pressure rows and700small array evaluations assumedunder2minutes,15minutes allowed.',
    final_charge_conclusion='unadjudicated',full_goal_complete=False,
    bindings={str(p):sha(p) for p in [Path(__file__),OUT/'canonical-affine-result.json',OLD/'metric/metric-128-g8.npz']}))
try:
    x=s.symbols('x');p0,p1,q0,q1,q2,q3=s.symbols('p0 p1 q0 q1 q2 q3')
    p=p0+p1*x;q=q0+q1*x+q2*x*x+q3*x**3
    assert s.expand(s.integrate(p*s.diff(q,x),(x,0,1))-(p.subs(x,1)*q.subs(x,1)-p.subs(x,0)*q.subs(x,0)-s.integrate(s.diff(p,x)*q,(x,0,1))))==0
    Model=bind(producer.base.initialize,OUT=OUT)();m=Model(128)
    pressures=np.array([[m.material.point(k)[name].copy() for name in ['Pr','Pg']] for k in range(17)],LD)
    lapse=producer.base.Lapse();controls=[]
    def evaluate(poly,t,side):
        ids=np.clip(np.searchsorted(poly.x,t,side=side)-1,0,len(poly.x)-2)
        dt=(t-poly.x[ids]).reshape((-1,)+(1,)*(poly.c.ndim-2))
        v=np.zeros((len(t),)+poly.c.shape[2:])
        for co in poly.c:v=v*dt+co[ids]
        return v
    for n,q in [(64,8),(128,4),(128,8)]:
        path=OLD/f'metric/metric-{n}-g{q}.npz';z=dict(np.load(path))
        d=dict(np.load(OLD/f'gr/source-{n}.npz'));field=np.load(OLD/f'gr/fields-{n}-g{q}.npz')
        _,poly=producer.base.representation(n)
        derivative=dict(delta_phi=field['U_t']/field['radius_E'],delta_Phi=np.zeros_like(field['U_t']))
        rates=[]
        for side in ['right','left']:
            dot=dict(d);dot.update({k:evaluate(p.derivative(),d['t'],side) for k,p in poly.items()})
            c=producer.base.geometry.constraints.centers(lapse.response,dot,derivative,q)
            rates.append(c['delta_lambda'][:,:-1])
        denom=max(np.max(np.sum(abs(z['actual_delta_lambda_rate']),axis=1)),1e-290)
        err=float(np.max(np.sum(abs(rates[0]-z['actual_delta_lambda_rate']),axis=1))/denom)
        assert err<1e-12,(n,q,'Original right derivative reproduction',err)
        # Keep the original causally spliced history and apply only the exact
        # difference in the source-polynomial derivative convention.
        z['actual_delta_lambda_rate']=z['actual_delta_lambda_rate']+(rates[1]-rates[0])
        np.savez_compressed(OUT/f'metric-left-{n}-g{q}.npz',**z)
        controls.append(dict(clock=n,order=q,right_reproduction=err,
            derivative_change_relative=float(np.max(np.sum(abs(rates[1]-rates[0]),axis=1))/denom)))
    metric=dict(np.load(OLD/'metric/metric-128-g8.npz'))
    left_metric=np.load(OUT/'metric-left-128-g8.npz')
    clock=metric['t'];a=m.material.a.astype(LD)
    def at(times):
        ids=np.array([np.argmin(abs(clock-t)) for t in times])
        assert np.max(abs(clock[ids]-times))<1e-18
        return [metric[k][ids].astype(LD) for k in ['delta_u','delta_lambda','delta_u_t','actual_delta_lambda_rate']]
    def pressure(times):
        j=np.clip(np.searchsorted(m.t,times,side='right')-1,0,15)
        w=((times-m.t[j])/(m.t[j+1]-m.t[j])).astype(LD)
        return (1-w[:,None,None])*pressures[j]+w[:,None,None]*pressures[j+1]
    rows=[];values={}
    for n in [64,128]:
        with np.load(OLD/f'sweep-1/photons/return-{n}.npz') as z:
            times=z['joint_stage_times'];weights=z['joint_stage_weights'];edges=z['actual_step_edges']
        pp=pressure(times);u,lam,ut,lt=at(times)
        direct=-a*((pp[:,0]+2*pp[:,1])*ut+pp[:,0]*lt)
        pe=pressure(edges);ue,le,*_=at(edges);h=np.diff(edges).astype(LD)
        assert np.all(h>0) and np.array_equal(times[1::2],edges[1:])
        b=np.array([LD(3)/4,LD(1)/4]);iv=lambda v:h[:,None]*np.einsum('j,kjn->kn',b,v.reshape(-1,2,m.n))
        P=pe[:,0]+2*pe[:,1];R=pe[:,0]
        primitive=-a*(P[1:]*ue[1:]-P[:-1]*ue[:-1]+R[1:]*le[1:]-R[:-1]*le[:-1]
            -np.diff(P,axis=0)/h[:,None]*iv(u)-np.diff(R,axis=0)/h[:,None]*iv(lam))
        original=iv(direct);delta=primitive-original
        ids=np.array([np.argmin(abs(clock-t)) for t in times])
        left_change=iv(-a*pp[:,0]*(left_metric['actual_delta_lambda_rate'][ids].astype(LD)-lt))
        values[n]=[v.sum(0,dtype=LD) for v in [original,primitive,delta,left_change]]
        rows.append(dict(clock=n,thermal_erg=[float(v.sum()) for v in values[n]],
            unweighted_u_rate_integral_error=float(np.max(abs(np.sum(weights[:,None]*ut,axis=0)-(ue[-1]-ue[0])))/max(np.max(abs(ue[-1]-ue[0])),LD('1e-290')))))
        np.savez_compressed(OUT/f'pressure-work-{n}.npz',edges=edges,direct=original,by_parts=primitive,correction=delta,left_change=left_change)
    diff=[x-y for x,y in zip(values[64],values[128])]
    native=read(OUT/'canonical-affine-result.json')['thermal_difference_erg'][0]
    result=dict(classification='Counterexample candidate',symbolic_identity_passed=True,rows=rows,derivative_controls=controls,
        labels=['direct_pressure_work','by_parts_pressure_work','by_parts_minus_direct','left_minus_right'],
        coarse_minus_fine_erg=[float(v.sum()) for v in diff],native_difference_erg=native,
        estimated_native_difference_after_pressure_primitive=native+float(diff[2].sum()),
        estimated_native_difference_after_left_derivative=native+float(diff[3].sum()),
        endpoint_metric_jumps_separately_unadjudicated=True,physical_steps=0,
        scope='Known explicit pressure term only; retained value quadrature and jumps are not certified. Not an evolved energy correction.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'pressure-work.json',result);print(json.dumps(result),flush=True)
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'pressure-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
