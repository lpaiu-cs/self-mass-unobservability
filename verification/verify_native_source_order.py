"""Independent order-verdict replay and exact static-condensation control."""
import json
import time
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csc_matrix
import def_native_source_order as task


def condensation_control():
    rng=np.random.default_rng(95)
    model=SimpleNamespace(size=10,n=2,degree=2,grid=np.arange(5)/2,
        cells=np.arange(3),indices=np.arange(10).reshape(5,2))
    model.permutation=np.argsort(np.r_[np.repeat(model.grid,2),[1.,2.]],kind='stable')
    a=rng.normal(size=(12,12));a=(a+a.T)/2
    # The two separate element-interior pairs must not couple directly.
    a[np.ix_([2,3],[6,7])]=0;a[np.ix_([6,7],[2,3])]=0
    a[np.diag_indices(12)]=np.sum(abs(a),axis=1)+1
    p=model.permutation;matrix=a[p,:][:,p];rhs=rng.normal(size=12)
    factor=task.condense_factor(csc_matrix(matrix),model)
    calculated=factor.solve(rhs);reference=np.linalg.solve(matrix,rhs)
    error=float(max(abs(calculated-reference))/max(abs(reference)))
    residual=float(max(abs(matrix@calculated-rhs))/max(abs(rhs)))
    assert error<1e-12 and residual<1e-12
    import sympy as s
    B,C,D,A,x,y,rb,rk=s.symbols('B C D A x y rb rk',nonzero=True)
    recovered=(rb-C*x)/B
    assert s.expand((D*recovered+A*x-rk)-((A-D*C/B)*x-(rk-D*rb/B)))==0
    return dict(classification='Proven',passed=True,
        identity='Eliminate B*y+C*x=rb, solve (A-D*B^-1*C)*x=rk-D*B^-1*rb, and reconstruct y. No basis coefficient is omitted.',
        numerical_classification='Counterexample candidate',dense_control_relative=error,residual=residual)


def main():
    start=time.monotonic();out=task.OUT;result=json.loads((out/'order-result.json').read_text())
    plan=json.loads((out/'order-budget-review.json').read_text())
    for name,h in plan['bindings'].items():assert task.task.old.photons.digest(task.task.old.ROOT/name)==h,name
    a=np.load(out/'consistent-p6-32.npz');b=np.load(out/'consistent-p6-64.npz');c=np.load(out/'consistent-p4-64.npz')
    replay_time=task.compare(a,{f:b[f][::2] for f in ['temperature','velocity','scalar']});space=task.compare(c,b)
    passed=True;worst={}
    for f,g in plan['gates'].items():
        assert abs(replay_time[f]-result['time_relative'][f])<1e-14
        assert abs(space[f]-result['space_relative'][f])<1e-14
        passed &= replay_time[f]<g and space[f]<g
        difference=abs(c[f]-b[f]);j,i=np.unravel_index(np.argmax(difference),difference.shape)
        worst[f]=dict(native_index=int(i),radius_fraction=float(b['radius'][i]),p4=float(c[f][j,i]),p6=float(b[f][j,i]),reference_maximum=float(max(abs(b[f]).ravel())))
    assert bool(passed)==result['passed']
    max_residual=max_heat=0.;memory=0.
    for label in ['consistent-p6-pilot','condensed-p6-pilot','consistent-p6-32','consistent-p6-64']:
        d=np.load(out/(label+'.npz'));report=json.loads((out/(label+'.json')).read_text())
        assert all(np.all(np.isfinite(d[f])) for f in d.files)
        assert all(np.max(abs(d[f][0]))==0 for f in plan['gates'])
        max_residual=max(max_residual,report['max_linear_residual']);max_heat=max(max_heat,report['max_heat_identity'])
        memory=max(memory,report['memory_GB'])
    assert max_residual<1e-9 and max_heat<1e-9
    pilot=task.compare(np.load(out/'condensed-p6-pilot.npz'),np.load(out/'consistent-p6-pilot.npz'))
    assert max(pilot.values())<1e-7
    checks=condensation_control()
    report=dict(classification='Counterexample candidate',artifact_checks_passed=True,
        consistent_mass_order_passed=bool(passed),worst=worst,condensation=checks,
        full_banded_pilot_relative=pilot,maximum_original_equation_residual=max_residual,
        maximum_original_heat_identity=max_heat,memory_GB=memory,
        seconds=time.monotonic()-start,physical_source_profile_certified=False,
        moving_surface_solved=False,final_dynamic_charge_solved=False,full_goal_complete=False)
    task.task.write(out/'order-audit.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':main()
