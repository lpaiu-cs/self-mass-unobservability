"""One explicitly reassessed completion run of the registered1792-cell path."""
from pathlib import Path
import json
import signal
import time
import def_native_energy_flow as task


def main():
    out=task.OUT;assert not (out/'cells-1792.npz').exists();start=time.monotonic();signal.alarm(90)
    budget=json.loads((out/'final-path-budget.json').read_text());assert not budget['allowed']
    task.write(out/'final-path-reassessment.json',dict(classification='Counterexample candidate',
        reason='The actual896-cell path completed and closed energy. The final-path conservative forecast is83.55seconds, with about5.53seconds measured setup, versus80.15seconds left. A single additional10seconds of production budget finishes the originally registered convergence decision; no new grid or physical tolerance is introduced.',
        original_production_budget_seconds=240,reassessed_total_production_budget_seconds=250,
        final_path_cap_seconds=90,forecast_evolution_and_native_seconds=budget['forecast_last_path_seconds'],measured_previous_setup_seconds=budget['elapsed_seconds']-budget['coarse_actual_seconds'],
        unchanged_native_cap=3000,remaining_native_calls=3000-budget['runtime_native_calls'],
        unchanged=['1792cells','original3.434ms horizon','all conservation/EOS/refinement gates'],
        stop='No further budget increase in this phase. On timeout or a native/domain failure retain accepted progress for a later scientific and resource reassessment.',
        producer_sha256=task.bank.old.photons.digest(Path(task.__file__)),EOS_sha256=task.bank.old.photons.digest(Path(task.temperature.__file__))))
    coarse=json.loads((out/'cells-896.json').read_text());fine=task.Flow(1792).run('cells-1792')
    errors={k:abs(coarse[k]/fine[k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    result=dict(classification='Counterexample candidate',passed=coarse['passed'] and fine['passed'] and max(errors.values())<.02,refinement=errors,
        final_path_including_setup_seconds=time.monotonic()-start,conservative_energy_evolved=True,energy_defect_added_as_heat=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    task.write(out/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
