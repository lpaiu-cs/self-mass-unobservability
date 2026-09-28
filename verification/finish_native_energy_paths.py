"""Use the affordable registered896-cell path to decide the final path budget."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import def_native_energy_flow as task


def main():
    out=task.OUT;assert not (out/'result.json').exists();start=time.monotonic();signal.alarm(100)
    task.write(out/'sequential-path-plan.json',dict(classification='Counterexample candidate',remaining_total_seconds=100,
        reason='The combined conservative forecast exceeds the remaining budget by9.23seconds, but the registered896-cell path itself is affordable and provides an actual energy-evolution speed. Execute that path first, preserve it, then decide the1792-cell path from its measured cost.',
        first_path_forecast_seconds=21.584329309201415,
        last_path_rule='Use the larger of5.2 times the896-cell evolution cost and measured1792 RHS cost times2.5 times the actual896 steps times2 RHS times1.15; add the entire remaining native-call allowance. The2.5 step ratio includes25percent margin above2x spatial scaling and exceeds the prior entropy-evolution1335/567 ratio. The new fine-grid speed is still an assumption.',
        paths=[896,1792],unchanged=['horizon','conservation gates','EOS native-call cap'],stop='Do not start1792 if its forecast exceeds the remaining100second hard clock; persist accepted intermediate states.',
        source_sha256=task.bank.old.photons.digest(Path(task.__file__)),EOS_source_sha256=task.bank.old.photons.digest(Path(task.temperature.__file__))))
    coarse=task.Flow(896).run('cells-896');elapsed=time.monotonic()-start
    measured=json.loads((out/'measured-flow-budget.json').read_text());stage=measured['rows'][1]['maximum_RHS_seconds'];calls=int(np.load(out/'runtime-columns.npz')['runtime_native_calls'])
    forecast=max(5.2*coarse['seconds'],stage*2.5*coarse['steps']*2*1.15)+(3000-calls)*.022
    allowed=elapsed+forecast<100
    task.write(out/'final-path-budget.json',dict(classification='Counterexample candidate',allowed=allowed,coarse_actual_seconds=coarse['seconds'],coarse_actual_steps=coarse['steps'],elapsed_seconds=elapsed,forecast_last_path_seconds=forecast,hard_remaining_seconds=100-elapsed,runtime_native_calls=calls))
    assert allowed,'Final registered path forecast exceeds remaining budget'
    fine=task.Flow(1792).run('cells-1792');errors={k:abs(coarse[k]/fine[k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    result=dict(classification='Counterexample candidate',passed=coarse['passed'] and fine['passed'] and max(errors.values())<.02,refinement=errors,
        seconds=time.monotonic()-start,conservative_energy_evolved=True,energy_defect_added_as_heat=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    task.write(out/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
