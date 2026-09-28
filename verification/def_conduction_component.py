"""Check the added heat response in its own norm, not the larger old source."""
from pathlib import Path
import json
import time
import numpy as np
import def_conductive_cauchy as task

PARENT=task.OUT
OUT=PARENT/'component'


def main():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(task.__file__),PARENT/'plan.json',PARENT/'result.json',PARENT/'heat-only-64.npz',PARENT/'heat-only-64.json']
    task.core.ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(task.h.ROOT).as_posix():task.h.digest(p) for p in paths},
        reason='Combined readout is dominated by the previously solved reactive source. Passing its norm cannot validate the much smaller added conduction response. Isolate the already declared linear source component without increasing grid resolution, horizon, collision samples or relaxing the old criteria.',
        cases=['heat-16','heat-32','coarse-bank-64','outer-64'],
        gates=dict(time_relative=.02,time_order=1.5,coefficient_relative=.02,outer_relative=.002,heat_balance=2e-13),
        budget=dict(steps=176,paths=4,hard_seconds=60,forecast_seconds=52,
                    basis='Measured setup 3.71 s/path and 0.1143 s/step, times 1.5 allowance. Reuse the already saved fine heat-only 64 path.',automatic_expansion=False)))
    task.OUT=OUT;start=time.monotonic();cases={}
    for steps in [16,32]:cases[str(steps)]=task.solve(steps,task.core.OUT/'bank.npz',False,label=f'heat-{steps}')
    cases['64']=json.loads((PARENT/'heat-only-64.json').read_text())
    cases['coarse']=task.solve(64,task.core.OUT/'coarse-bank.npz',False,label='coarse-bank-64')
    cases['outer']=task.solve(64,task.core.OUT/'bank.npz',False,3,label='outer-64')
    comparisons={}
    for field in ['velocity_mass_RMS_m_s','scalar_mass_RMS']:
        series=lambda name:np.array([x[field] for x in cases[name]['history']])
        a,b,c=series('16'),series('32'),series('64');norm=max(abs(c).max(),1e-100)
        d1=np.max(abs(a-b[::2]))/norm;d2=np.max(abs(b-c[::2]))/norm
        comparisons[field]=dict(time_previous=float(d1),time_last=float(d2),order=float(np.log2(d1/d2)),
            coefficients=float(np.max(abs(c-series('coarse')))/norm),outer=float(np.max(abs(c-series('outer')))/norm))
    # Check endpoint superposition separately for each state component, in the
    # smaller heat-only norm. Do not divide by the much larger combined state.
    total=np.load(PARENT/'combined-64.npz');heat=np.load(PARENT/'heat-only-64.npz');old=np.load(task.coupled.OUT/'fine-64.npz')
    errors={k:[float(np.max(abs(total[k][:,j]-old[k][:,j]-heat[k][:,j]))/max(abs(heat[k][:,j]).max(),1e-100)) for j in range(4)] for k in ['response','velocity']}
    passed=all(x['time_last']<.02 and x['order']>1.5 and x['coefficients']<.02 and x['outer']<.002 for x in comparisons.values())
    passed=passed and max(x['heat_telescoping'] for x in cases.values())<2e-13
    result=dict(classification='Counterexample candidate',passed=passed,comparisons=comparisons,
        heat_norm_superposition_components=errors,seconds=time.monotonic()-start,
        full_temperature_feedback=False,core_interface_is_physical_boundary=False,continuum_physical_error_certified=False,full_dynamic_charge_solved=False)
    task.core.ex.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
