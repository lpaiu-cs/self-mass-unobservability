"""Preserve failed hybrid solves; audit a residual-decreasing coupled Newton path."""
import json,shutil,sys
from types import FunctionType,SimpleNamespace
import numpy as np
import gr_heat_primitive_scaling as scaling

original=scaling.original;g=original.g;OUT=g.OUT/'gr-heat-primitive-newton';traces=[]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();scaling.verify()
    result=json.loads((scaling.OUT/'result.json').read_text());assert not result['all_passed']
    audits=json.loads((scaling.OUT/'root-diagnostics.json').read_text())['rows']
    assert all(d['passed'] for a in audits for d in a['finite_differences'])
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='a2d50b8',max_iterations=30,max_backtracks=24,Newton_scaled_residual=1e-12,
        logT_step_cap=.1,velocity_step_cap=.01,
        change='Keep every original target, nonlinear equation, analytic Jacobian, starting point and acceptance gate. Replace the hybrid iteration by coupled Newton with a recorded strict decrease of maximum scaled residual. A separately stricter 1e-12 solver stopping criterion is used; no physical/root gate is relaxed.',
        scope='Post-failure numerical root repair. Full and rejected Newton trials are recorded; local convergence is not global primitive uniqueness, certified native root error, conserved cell closure or GR evolution.')
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_heat_primitive_newton.py',original.OUT/'manifest.json',scaling.OUT/'manifest.json']})
    save('plan.json',plan);shutil.copy2(original.OUT/'symbolic.json',OUT/'symbolic.json')


def root(fun,x0,jac,method,options):
    plan=json.loads((OUT/'plan.json').read_text());x=np.asarray(x0,float);history=[];success=False;message='iteration limit'
    for iteration in range(plan['max_iterations']):
        f=fun(x);J=jac(x);merit=float(abs(f).max())
        row=dict(iteration=iteration,x=x.tolist(),residual=f.tolist(),merit=merit,
            Jacobian_condition=float(np.linalg.cond(J)),trials=[])
        history.append(row)
        if merit<=plan['Newton_scaled_residual']:success=True;message='scaled residual converged';break
        delta=np.linalg.solve(J,-f)
        delta*=min(1.,plan['logT_step_cap']/max(abs(delta[0]),1e-300),
            plan['velocity_step_cap']/max(abs(delta[1]),1e-300))
        for backtrack in range(plan['max_backtracks']):
            proposed=x+delta*(.5**backtrack)
            if abs(proposed[1])>=.5:
                row['trials'].append(dict(backtrack=backtrack,outside_velocity_bracket=True));continue
            residual=fun(proposed);new_merit=float(abs(residual).max())
            row['trials'].append(dict(backtrack=backtrack,x=proposed.tolist(),residual=residual.tolist(),merit=new_merit))
            if new_merit<merit:
                x=proposed;row['accepted_backtrack']=backtrack;break
        else:message='no residual-decreasing step';break
    traces.append(dict(case_index=len(traces),history=history,success=success,message=message))
    save('Newton-traces.json',dict(classification='Counterexample candidate',rows=traces))
    return SimpleNamespace(x=x,success=success,message=message)


def run():
    FunctionType(original.run.__code__,dict(original.run.__globals__,OUT=OUT,save=save,root=root,verify=verify))()


def verify():
    FunctionType(original.verify.__code__,dict(original.verify.__globals__,OUT=OUT))()
    plan=json.loads((OUT/'plan.json').read_text())
    for i in plan['cells']:
        for case in range(len(plan['probes'])):
            a=np.load(OUT/f'cell-{i}-case-{case}.npz');b=np.load(original.OUT/f'cell-{i}-case-{case}.npz')
            assert np.array_equal(a['conserved'],b['conserved']) and np.array_equal(a['true_primitive'],b['true_primitive'])
    for row in json.loads((OUT/'Newton-traces.json').read_text())['rows']:
        history=row['history']
        assert all(b['merit']<a['merit'] for a,b in zip(history,history[1:]))
    print('PASS identical primitive targets and decreasing Newton residual histories',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
