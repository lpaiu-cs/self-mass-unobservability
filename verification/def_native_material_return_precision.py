"""Counterexample candidate: resolve measured deep-force probe roundoff."""
from pathlib import Path
from types import FunctionType
import json,signal,sys,time
import numpy as np
import def_native_feedback_return as base

OUT=base.OUT/'probe-scaled';write=base.write;sha=base.sha
class Material(base.Material):
    def rhs(self,t,z,probe=1.):return super().rhs(t,z,8*probe)
    run=FunctionType(base.old.base.Material.run.__code__,dict(vars(base.old.base),OUT=OUT),argdefs=base.old.base.Material.run.__defaults__)


def run():
    assert not OUT.exists();OUT.mkdir();prior=json.loads((base.OUT/'production.json').read_text());scan=json.loads((base.OUT/'probe.json').read_text())
    assert not prior['passed'] and all(r['passed'] for r in prior['paths'][:2])
    assert max(max(r['relative']) for r in scan['rows'][-2:])<.002
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        failure='Completed reference64 material return fails endpoint momentum directional comparison0.30497percent against unchanged0.2percent. Errors are concentrated in deep cells0-8 and increase as probes shrink:1.34percent at0.125/0.25,0.305percent at0.5/1,0.0341percent at4/8,0.0135percent at8/16.',
        repair='The deep nearly balanced force subtraction reaches its arithmetic floor. Use8x finite probe within the same Richardson directional derivative, retaining actual donor branches and original physical amplitude. Check both4/8 and8/16 on the new evolved state; do not change the0.002 gate.',
        reuse='Retain both accepted reference128 material paths and all accepted photons. Recompute only failed128/reference64 using identical SSP,CFL,geometry and physical fluxes. No new mesh, clocks,EOS banks or horizon.',
        budget=dict(pilot_seconds=25,production_seconds=250,CPU_threads=1,memory_GB=3),
        forecast='Use the measured79.88s failed full path, scaled by new/old prefix time per raw call, plus5s. Require2x forecast below250s.',
        stop='Stop on prefix gate, forecast, completed directional gate or250s. Do not automatically search further probe scales or physical resolutions.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),base.OUT/'production.json',base.OUT/'probe.json',base.OUT/'production-producer.py']}))
    start=time.monotonic();signal.signal(signal.SIGALRM,base.old.base.flow.old.optical.timeout);signal.alarm(25)
    m=Material(64,128);pilot=m.run(128,'pilot-128',2)
    oldpilot=json.loads((base.OUT/'pilot.json').read_text())['rows'][2]
    rate=(pilot['seconds']/pilot['raw_owner_calls'])/(oldpilot['seconds']/oldpilot['raw_owner_calls'])
    forecast=prior['paths'][2]['seconds']*rate+5
    p=dict(classification='Counterexample candidate',row=pilot,forecast_seconds=forecast,upper_seconds=2*forecast,eligible=pilot['passed'] and 2*forecast<250,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',p);print(json.dumps(p),flush=True);signal.alarm(0)
    if not p['eligible']:return
    write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',budget_seconds=250,forecast_seconds=forecast,upper_seconds=2*forecast,
        bindings={str(q):sha(q) for q in [Path(__file__),Path(base.__file__),OUT/'plan.json',OUT/'pilot.json']}))
    start=time.monotonic();signal.alarm(250)
    try:
        m=Material(64,128);row=m.run(128,'steps-128-reference-64',restart='pilot-128');row['physical_branch_ratio']=m.physical_branch_ratio
        d=np.load(OUT/'steps-128-reference-64.npz');z=d['delta_scaled'];t=float(d['time']);a=m.rhs(t,z)[0];b=m.rhs(t,z,2.)[0]
        reverse=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist()
        result=dict(classification='Counterexample candidate',passed=row['passed'] and max(reverse)<.002,path=row,reverse_probe_8_16_relative=reverse,seconds=time.monotonic()-start)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
        if result['passed']:
            audit=dict(prior,passed=True,paths=prior['paths'][:2]+[row],original_material_return_verdict=False,
                precision_repair_seconds=time.monotonic()-start,seconds=prior['seconds']+time.monotonic()-start,
                reverse_probe_8_16_relative=reverse,accepted_paths={f'{s}/{r}':str(base.material_path(s,r)) for s,r in [(64,128),(128,128),(128,64)]})
            write(base.OUT/'material-audit.json',audit)
    finally:signal.alarm(0)

if __name__=='__main__':run()
