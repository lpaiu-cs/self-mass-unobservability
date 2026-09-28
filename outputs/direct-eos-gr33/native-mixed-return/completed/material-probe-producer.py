"""Inspect arithmetic probe dependence on the failed actual material state."""
from pathlib import Path
import json,resource,time
import numpy as np
import return_native_mixed_gr as run

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha


def main():
    files=[Path(__file__),Path(run.__file__),OUT/'material_pilot-1-receipt.json',
        OUT/'sweep-1/material/pilot-64.npz',OUT/'sweep-1/material/pilot-64.json',OUT/'sweep-1/photons/result.json']
    assert not (OUT/'material-probe-plan.json').exists()
    write(OUT/'material-probe-plan.json',dict(classification='Counterexample candidate',
        failure='Actual first material prefix fails the unchanged0.2percent half/nominal directional comparison.',
        question='Does the failure shrink on reducing the artificial finite-difference probe, or on increasing it? Inspect the same saved state; do not change physical amplitude, EOS, grid, clock or donor rule.',
        factors=[1/64,1/32,1/16,1/8,1/4,1/2,1,2],base_arithmetic_probe=512,
        selection='Choose the smallest central factor whose half/nominal/double pairwise L1 directional discrepancies are below0.2percent. This only admits rerunning the failed two-step prefix; all future trajectory checks remain required.',
        budget_seconds=30,total_action_seconds=run.TOTAL,new_evolution_steps=0,
        bindings={str(p):sha(p) for p in files}))
    run.initialize(1);m=run.reuse.coupled.Material(128,64);d=np.load(OUT/'sweep-1/material/pilot-64.npz')
    t=float(d['time']);z=d['delta_scaled'];factors=read(OUT/'material-probe-plan.json')['factors'];values=[];rows=[]
    for factor in factors:
        m.min_probe=np.inf;m.probe_error=0.;mark=time.monotonic();v,_,cfl=m.rhs(t,z,factor);values.append(v)
        rows.append(dict(factor=factor,actual_epsilon=m.min_probe,forward_indicator=m.probe_error,cfl_seconds=cfl,seconds=time.monotonic()-mark))
    controls=[]
    for i in range(1,len(values)-1):
        errors=[(np.sum(abs(values[j]-values[i]),axis=1)/np.maximum(np.sum(abs(values[i]),axis=1),1.)).tolist() for j in [i-1,i+1]]
        controls.append(dict(factor=factors[i],half_nominal_double=errors,maximum=float(np.max(errors))))
    choices=[r for r in controls if r['maximum']<.002]
    np.savez_compressed(OUT/'material-probe-rates.npz',t=t,state=z,factors=factors,rates=values)
    result=dict(classification='Counterexample candidate',passed=bool(choices),rows=rows,controls=controls,
        selected_factor=choices[0]['factor'] if choices else None,
        selected_nominal_arithmetic_probe=512*choices[0]['factor'] if choices else None,
        physical_input_unchanged=True,trajectory_certified=False)
    write(OUT/'material-probe-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.mixed.inc.native.deadline(30)
    start=time.monotonic();cpu=time.process_time();error=None;receipt=OUT/'material-probe-receipt.json';assert not receipt.exists()
    try:
        assert sum(read(p)['seconds'] for p in OUT.rglob('*-receipt.json'))+30<=run.TOTAL
        main()
    except Exception as exc:error=repr(exc);raise
    finally:write(receipt,dict(action='material_probe',seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
