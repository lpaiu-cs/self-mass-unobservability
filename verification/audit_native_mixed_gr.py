"""Independent exact-source check for the analytic reciprocal derivative."""
from pathlib import Path
import json, resource, time
import numpy as np
from numpy.polynomial import legendre as leg
import couple_native_mixed_gr as run


def derivative_check():
    errors=[];C=run.C
    for order in [4,8]:
        m=run.wave.Response.__new__(run.wave.Response);m.order=order;m.t=np.array([0.,.17,.43,.79,1.])/C
        m.xfaces=np.array([-1.,0.,1.]);m.mid=np.array([-.5,.5]);m.half=np.array([.5,.5]);gx,gw=leg.leggauss(order)
        m.x=(m.mid[:,None]+m.half[:,None]*gx).ravel();m.ids=np.repeat(np.arange(2),order)
        m.dx=(m.half[:,None]*np.broadcast_to(gw,(2,order))).ravel();m.tx=np.array([-.75,0.,.3,2.])
        m.inverse=np.linalg.inv(leg.legvander(np.broadcast_to(gx,(2,order)),order-1));m.gx,m.gw=leg.leggauss((order+3)//2)
        m.distance=abs(m.tx[:,None]-m.x[None,:])/C;m.sign=np.sign(m.tx[:,None]-m.x[None,:])
        # Exact S(t,x)=2+x^2+c*t*(1+x) on [-1,1]. Nonzero initial
        # source and both causal fronts distinguish it from a slope-only test.
        source=m.dx*(2+m.x*m.x+(C*m.t)[:,None]*(1+m.x))
        actual=run.second_derivative(m,source);exact=[]
        for t in C*m.t:
            row=[]
            for x in m.tx:
                lo=max(-1,x-t);hi=min(1,x+t)
                integral=(hi+hi*hi/2-lo-lo*lo/2)/2 if hi>lo else 0.
                front=sum((2+y*y)/2 for y in [x-t,x+t] if -1<=y<=1)
                local=2+x*x+t*(1+x) if -1<=x<=1 else 0.
                row.append(integral+front-local)
            exact.append(row)
        errors.append(float(np.max(abs(actual-exact))))
    assert max(errors)<1e-11,errors
    return dict(classification='Proven',passed=True,maximum_absolute_errors=errors,
        scope='Analytic Uxx for the stated polynomial source including nonzero initial fronts. No physical continuum error is inferred.')


def main():
    out=run.OUT;assert not (out/'audit-plan.json').exists()
    files=[Path(__file__),Path(run.__file__),out/'result.json',out/'repair-plan.json',out/'field-128-g8.npz']
    run.write(out/'audit-plan.json',dict(classification='Counterexample candidate',budget_seconds=60,
        claim='Independently verify the analytic second derivative with nonzero initial source, mixed identities, frozen inputs and the actual applied source sum.',
        gates=dict(derivative=1e-11,source_sum=1e-12),
        limits='This audit does not certify physical source interpolation, exterior mixed terms, material return or nonlinear closure.',
        bindings={str(p):run.sha(p) for p in files}))
    check=derivative_check();symbolic=run.symbolic();z=np.load(out/'field-128-g8.npz')
    density=sum(z['density_'+k] for k in ['trace','static_gradient','reciprocal','constraint'])
    relative=float(np.max(abs(z['source']/z['source_dx']-density))/max(np.max(abs(density)),1e-290));assert relative<1e-12
    plan=run.read(out/'plan.json')
    for p,h in plan['bindings'].items():
        assert run.sha(out/'registered-producer.py' if p==str(Path(run.__file__)) else p)==h,p
    result=run.read(out/'result.json');assert result['passed'],result['controls']
    report=dict(classification='Counterexample candidate',passed=True,derivative=check,mixed_symbolic=symbolic,
        applied_source_sum_relative=relative,original_inputs_preserved=True,
        reciprocal_compact_source_applied=True,full_goal_complete=False)
    run.write(out/'audit.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.inc.native.deadline(60)
    started=time.monotonic();cpu=time.process_time();error=None;receipt=run.OUT/'audit-receipt.json';assert not receipt.exists()
    try:main()
    except Exception as exc:error=repr(exc);raise
    finally:run.write(receipt,dict(seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=run.sha(__file__)))
