"""Measure warm kernel costs, then finish only the two registered gas paths."""
from pathlib import Path
import json
import signal
import time
import runpy
import numpy as np
import def_native_energy_flow as task


def main():
    out=task.OUT;assert not (out/'result.json').exists();start=time.monotonic();signal.alarm(125)
    task.write(out/'flow-completion-plan.json',dict(classification='Counterexample candidate',
        decision='Reuse the completed224-cell pilot and native trajectory cache. The original cell-squared extrapolation included one-time EOS construction. Measure actual896/1792-cell RHS costs after adding only a missing-domain fast path, before any further evolution.',
        source_sha256=task.bank.old.photons.digest(Path(task.__file__)),EOS_source_sha256=task.bank.old.photons.digest(Path(task.temperature.__file__)),
        remaining_total_seconds=125,additional_native_calls='Unchanged3000-call runtime cap minus already saved calls',
        prior_benchmark_seconds=36.987579548,repair='Skip fixed-temperature pressure iterations only after relative pressure change is below1e-12; compare recovered states against the frozen three-iteration producer. Save accepted intermediate states for any subsequent stop.',
        benchmark='Three actual RHS calls at initial and interpolated pilot endpoint states on each registered grid. Offset the inverse seed by0.03 in logT so recovery cost is included. No time steps are advanced by this benchmark.',
        forecast='2 RHS per step, final/initial measured CFL with50percent step margin,15percent output overhead, plus0.022s for every remaining allowed native call.',
        stop='Do not start remaining evolution if the measured forecast exceeds remaining time. No new grid, horizon, gate or native-call extension.'))
    pilot=np.load(out/'pilot-224.npz');rows=[];end=.0034344311179287023
    before=runpy.run_path(str(out/'flow-before-pressure-convergence.py'))['Flow'](224)
    after=task.Flow(224);before.seed[:]-=.03;after.seed[:]-=.03
    a=before.primitive(pilot['U']);b=after.primitive(pilot['U']);error=float(np.max(abs(a-b)/np.maximum(abs(a),1e-30)))
    task.write(out/'pressure-convergence-control.json',dict(classification='Counterexample candidate',maximum_relative=error,passed=error<1e-11))
    assert error<1e-11,'Pressure fixed-point optimization changed the recovered state'
    for n in [896,1792]:
        m=task.Flow(n);b=m.base
        rho=np.interp(b.x,pilot['x_cm'],pilot['rho'])/m.eos.rho0
        v=np.interp(b.x,pilot['x_cm'],pilot['velocity_cm_s'])/task.C
        theta=np.interp(b.x,pilot['x_cm'],np.log(pilot['T']))
        final=m.conserved(rho,v,theta,b.a)[0]
        costs=[];steps=[]
        for U,t in [(m.initial,0.),(final,end)]:
            m.primitive(U);seed=m.seed.copy();sample=[]
            for _ in range(3):
                m.seed=seed-.03;clock=time.monotonic();_,_,dt=m.rhs(U,t);sample.append(time.monotonic()-clock)
            costs.append(max(sample));steps.append(int(np.ceil(end/dt*1.5)))
        rows.append(dict(cells=n,maximum_RHS_seconds=max(costs),forecast_steps=max(steps),forecast_seconds=2*max(costs)*max(steps)*1.15))
    native=np.load(out/'runtime-columns.npz');spent=int(native['runtime_native_calls']);native_allowance=(3000-spent)*.022
    elapsed=time.monotonic()-start;forecast=sum(r['forecast_seconds'] for r in rows)+native_allowance
    passed=elapsed+forecast<125
    task.write(out/'measured-flow-budget.json',dict(classification='Counterexample candidate',passed=passed,rows=rows,benchmark_seconds=elapsed,
        runtime_native_calls_already_spent=spent,native_call_time_allowance_seconds=native_allowance,forecast_remaining_seconds=forecast,hard_remaining_seconds=125-elapsed))
    assert passed,'Measured warm forecast still exceeds remaining budget'
    results=[task.Flow(n).run(f'cells-{n}') for n in [896,1792]]
    errors={k:abs(results[0][k]/results[1][k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    row=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in results) and max(errors.values())<.02,refinement=errors,
        seconds=time.monotonic()-start,conservative_energy_evolved=True,energy_defect_added_as_heat=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    task.write(out/'result.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':main()
