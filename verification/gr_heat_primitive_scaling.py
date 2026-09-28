"""Separate native primitive-Jacobian checks from MINPACK variable scaling."""
import json,shutil,sys
from types import FunctionType
import numpy as np
import gr_heat_primitive_inverse as original

g=original.g;OUT=g.OUT/'gr-heat-primitive-scaling';diagnostics=[]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.verify()
    result=json.loads((original.OUT/'result.json').read_text())
    assert [(r['cell'],r['case']) for r in result['rows'] if not r['passed']]==[(2972,0),(2972,2)]
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='a2d50b8',
        change='Reuse the exact prior run code object, target generation, EOS/analytic Jacobian and every acceptance criterion. Only the root invocation adds explicit unit variable scaling diag=[1,1]. Independently check the initial analytic Jacobian by symmetric native differences at two fixed resolutions before solving.',
        finite_steps=[[1e-4,1e-6],[5e-5,5e-7]],finite_Jacobian_row_relative_tolerance=1e-5,
        original_failed_cases=[[2972,0],[2972,2]],
        scope='Post-failure solver diagnosis and separately registered scaling candidate. Original failure retained; no physical, global-root or finite-time GR certification.')
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_heat_primitive_scaling.py',original.OUT/'manifest.json',original.OUT/'result.json']})
    save('plan.json',plan);shutil.copy2(original.OUT/'symbolic.json',OUT/'symbolic.json')


def root(fun,x0,jac,method,options):
    plan=json.loads((OUT/'plan.json').read_text());x0=np.asarray(x0,float)
    residual=fun(x0);matrix=jac(x0);differences=[]
    for steps in plan['finite_steps']:
        fd=np.empty((2,2))
        for k,h in enumerate(steps):
            change=np.eye(2)[k]*h;fd[:,k]=(fun(x0+change)-fun(x0-change))/(2*h)
        error=float(np.max(np.sum(abs(fd-matrix),axis=1)/np.maximum(np.sum(abs(matrix),axis=1),1e-100)))
        differences.append(dict(steps=steps,finite_matrix=fd.tolist(),row_relative_error=error,
            passed=error<plan['finite_Jacobian_row_relative_tolerance']))
    solution=original.root(fun,x0,jac=jac,method=method,options=dict(options,diag=[1.,1.]))
    final=fun(solution.x);record=dict(case_index=len(diagnostics),initial_residual=residual.tolist(),
        initial_analytic_Jacobian=matrix.tolist(),initial_Jacobian_condition=float(np.linalg.cond(matrix)),
        finite_differences=differences,unit_scale_solution=solution.x.tolist(),solver_success=bool(solution.success),
        final_residual=final.tolist())
    diagnostics.append(record);save('root-diagnostics.json',dict(classification='Counterexample candidate',rows=diagnostics))
    return solution


def run():
    # A private function globals mapping changes no historical module or
    # running scientific source. The nonlinear equations remain identical.
    function=FunctionType(original.run.__code__,dict(original.run.__globals__,OUT=OUT,save=save,root=root,verify=verify))
    function()


def verify():
    FunctionType(original.verify.__code__,dict(original.verify.__globals__,OUT=OUT))()
    plan=json.loads((OUT/'plan.json').read_text())
    for i in plan['cells']:
        for case in range(len(plan['probes'])):
            a=np.load(OUT/f'cell-{i}-case-{case}.npz');b=np.load(original.OUT/f'cell-{i}-case-{case}.npz')
            assert np.array_equal(a['conserved'],b['conserved']) and np.array_equal(a['true_primitive'],b['true_primitive'])
    print('PASS identical nonlinear primitive targets and separate unit-scaling results',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
