"""Reuse the measured six-variable Jacobian after the precision gate failed."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_native_whole_star_match as task


def main():
    out=task.OUT;old=json.loads((out/'result.json').read_text())
    assert not old['passed'] and not (out/'corrected-result.json').exists()
    history=json.loads((out/'progress.json').read_text())['history']
    pilot=json.loads((out/'pilot.json').read_text())
    # The original run evaluated exactly the same finite-difference Jacobian
    # twice; reuse the first 12 calls and retain both copies as evidence.
    columns=[]
    for i in range(6):
        a,b=history[1+2*i],history[2+2*i]
        step=a['parameters'][i]-b['parameters'][i]
        columns.append((np.array(a['residual'])-b['residual'])/step)
    jac=np.array(columns).T
    forecast=old['total_compute_seconds']+pilot['objective_seconds']*12+100
    assert forecast<600
    task.write(out/'correction-plan.json',dict(classification='Counterexample candidate',
        failed_result_preserved=True,cause='The root is accurate for the coarse ODE but not invariant under the registered tighter envelope/core evaluation. The dense envelope temperature error drives the residual, so root status is insufficient.',
        method='Reuse the already measured 6x6 Jacobian; at most3 Newton corrections at core sub8/native envelope tolerance2e-9. One final independent sub16/tolerance2e-10 evaluation. No new finite-difference Jacobian, unchanged equations, input inventory and original gates.',
        budget=dict(prior_compute_seconds=old['total_compute_seconds'],forecast_total_seconds=forecast,
            total_cap_seconds=600,hard_remaining_seconds=int(600-old['total_compute_seconds']),newton_corrections_max=3,new_full_jacobians=0,automatic_expansion=False),
        bindings={str(p.relative_to(task.ROOT)):task.envelope.digest(p) for p in [Path(__file__),Path(task.__file__),out/'result.json',out/'progress.json',out/'plan.json']},
        jacobian=jac.tolist()))
    signal.alarm(int(600-old['total_compute_seconds']));start=time.monotonic()
    model=task.Match(sub=8,tol=2e-9);x=np.array(old['parameters']);rows=[]
    for k in range(4):
        err,*_=model.branches(x)
        rows.append(dict(parameters=x.tolist(),residual=err.tolist()))
        task.write(out/'correction-progress.json',dict(classification='Counterexample candidate',history=rows))
        print('CORRECT',k,err.tolist(),flush=True)
        if max(abs(err))<1e-9:break
        assert k<3,'Registered correction cap'
        x-=np.linalg.solve(jac,err)
    matched=model.save(x,'corrected')
    fine=task.Match(sub=16,tol=2e-10);control=fine.save(x,'control')
    residual=np.array(control['residual'])
    correction=-np.linalg.solve(jac,residual)
    accepted=max(abs(residual))<1e-8 and abs(residual[-1])<1e-7
    result=dict(classification='Counterexample candidate',passed=bool(accepted),matched=matched,control=control,
        fixed_parameter_control_residual=residual.tolist(),
        linearized_parameter_correction_estimate=correction.tolist(),
        sensitivity_scope='Inverse cached Jacobian times residual is an estimate, not a rigorous root/charge error enclosure.',
        native_calls=model.env.calls+fine.env.calls,seconds=time.monotonic()-start,
        total_compute_seconds=old['total_compute_seconds']+time.monotonic()-start,
        corrected_mechanical_material_match=bool(accepted),steady_thermal_match=False,
        physical_zero_pressure_surface=False,final_dynamic_charge_solved=False,full_goal_complete=False)
    assert result['total_compute_seconds']<600
    task.write(out/'corrected-result.json',result);signal.alarm(0)
    print('FINAL',json.dumps(result),flush=True)


def final_repair():
    out=task.OUT;old=json.loads((out/'corrected-result.json').read_text())
    assert not old['passed'] and not (out/'final-result.json').exists()
    plan=json.loads((out/'correction-plan.json').read_text());jac=np.array(plan['jacobian'])
    charged=old['total_compute_seconds']+4
    task.write(out/'final-repair-plan.json',dict(classification='Counterexample candidate',
        cause='Measured scalar-coordinate subtraction loss and unscaled absolute ODE tolerances survive relative-tolerance refinement. Preserve both failed results. Repair the increment-to-core coordinate and absolute tolerance floor, without changing equations or gates.',
        method='One cached-Jacobian correction at sub8/tolerance2e-10, save that root and one independent sub16/tolerance2e-11 control. No repeated root search or extra paths after this finite repair.',
        budget=dict(prior_compute_seconds=charged,remaining_seconds=600-charged,
            forecast_additional_seconds_range=[90,165],total_cap_seconds=600,automatic_expansion=False),
        original_production_binding='The original plan binds original-source.py, preserved byte-for-byte; current source has the two diagnosed numerical corrections.',
        bindings={str(p.relative_to(task.ROOT)):task.envelope.digest(p) for p in
            [Path(__file__),Path(task.__file__),out/'original-source.py',out/'corrected-result.json',out/'scalar-coordinate-diagnosis.json']}))
    signal.alarm(int(600-charged));start=time.monotonic()
    model=task.Match(sub=8,tol=2e-10);x=np.array(old['matched']['parameters'])
    err,*_=model.branches(x);x-=np.linalg.solve(jac,err)
    matched=model.save(x,'final')
    control=task.Match(sub=16,tol=2e-11);reference=control.save(x,'final-control')
    errors=[max(abs(np.array(r['residual']))) for r in [matched,reference]]
    accepted=max(errors)<1e-8 and max(abs(r['residual'][-1]) for r in [matched,reference])<1e-7
    correction=-np.linalg.solve(jac,reference['residual'])
    value=dict(classification='Counterexample candidate',passed=bool(accepted),matched=matched,control=reference,
        junction_errors=errors,linearized_parameter_correction_estimate=correction.tolist(),
        native_calls=model.env.calls+control.env.calls,seconds=time.monotonic()-start,
        total_compute_seconds=charged+time.monotonic()-start,
        same_material_mechanical_match=bool(accepted),steady_thermal_match=False,
        physical_zero_pressure_surface=False,final_dynamic_charge_solved=False,full_goal_complete=False)
    task.write(out/'final-result.json',value);signal.alarm(0);print(json.dumps(value),flush=True)


if __name__=='__main__':
    import sys
    final_repair() if len(sys.argv)>1 and sys.argv[1]=='final' else main()
