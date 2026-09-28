"""Bounded correction solves; preserve the original three-refinement baseline."""
from pathlib import Path
import argparse
import inspect
import json
import signal
import textwrap
import time
import numpy as np
import def_gr_temperature_feedback as prior

OUT=prior.OUT
source=textwrap.dedent(inspect.getsource(prior.Problem.transform_pair))
source=source.replace('    def solve(E):','    solves=0\n    def solve(E):\n        nonlocal solves\n        refinements=3 if solves==0 else 1\n        solves+=1')
source=source.replace('for _ in range(3):','for _ in range(refinements):')
assert source.count('for _ in range(refinements):')==1
namespace=dict(vars(prior));exec(compile(source,'fast-transform-source.py','exec'),namespace)


class Problem(prior.Problem):
    transform_pair=namespace['transform_pair']


def prepare():
    assert not (OUT/'fast-pilot.json').exists();signal.alarm(90);start=time.monotonic()
    (OUT/'fast-transform-source.py').write_text(source)
    prior.write(OUT/'fast-plan.json',dict(classification='Counterexample candidate',
        reason='Original forecast639.05s exceeded the unchanged600s production cap; no full path ran. Preserve baseline three refinements; use one extended residual refinement for each tiny feedback correction, retaining the original1e-9 solve and1e-10 loop residual gates.',
        decision='Run only if the same16 frequencies agree within2e-8 relative in each native correction field and heat increment, and the new40-percent-margin forecast is below600s. Keep total720s including both pilots. No additional optimization or automatic budget increase after failure.',
        bindings={str(p):prior.go.task.digest(p) for p in [Path(__file__),Path(prior.__file__),OUT/'fast-transform-source.py',OUT/'plan.json',OUT/'pilot.json',OUT/'pilot.npz']}))
    prior.patch.install();tick=time.monotonic();p=Problem();setup=time.monotonic()-tick
    saved=np.load(OUT/'pilot.npz');rows=[];tick=time.monotonic()
    results=[p.transform_pair(prior.go.contour(int(k),12)[0]) for k in saved['ids']]
    elapsed=time.monotonic()-tick
    for i,result in enumerate(results):
        q,dq,E,dE=result;ref=saved['correction_q'][i];m=p.model
        errors=[]
        for field in range(2):
            refv=m.nativeV[field]@ref;delta=m.nativeV[field]@(dq-ref)
            for mask in ([np.ones(len(refv),bool),*m.original.masks] if field==0 else [np.ones(len(refv),bool)]):
                errors.append(float(np.linalg.norm(delta[mask])/max(np.linalg.norm(refv[mask]),1e-100)))
        errors.append(float(np.linalg.norm(dE-saved['correction_E'][i])/max(np.linalg.norm(saved['correction_E'][i]),1e-100)))
        rows.append(max(errors));assert rows[-1]<2e-8
        assert np.array_equal(q,saved['base_q'][i]) and np.array_equal(E,saved['base_E'][i])
    old=json.loads((OUT/'pilot.json').read_text());cost=json.loads((prior.patch.PRIOR/'pilot-budget.json').read_text())
    forecast=1.4*(setup+elapsed/16*2048+cost['inversion_forecast_seconds']+30)
    row=dict(old,forecast_seconds=forecast,setup_seconds=setup,transfer16_seconds=elapsed,
        seconds=old['seconds']+time.monotonic()-start,original_pilot_seconds=old['seconds'],
        correction_native_max_relative=max(rows),unchanged_baseline_bitwise=True,
        linear_residual=max(old['linear_residual'],p.error),feedback_residual=max(old['feedback_residual'],p.loop_residual))
    prior.write(OUT/'fast-pilot.json',row);signal.alarm(0);print('FAST FEEDBACK PILOT',json.dumps(row),flush=True)


def run():
    for path,h in json.loads((OUT/'fast-plan.json').read_text())['bindings'].items():assert prior.go.task.digest(Path(path))==h,path
    src=inspect.getsource(prior.run).replace("OUT/'pilot.json'","OUT/'fast-pilot.json'")
    (OUT/'fast-run-source.py').write_text(src)
    ns=dict(vars(prior),Problem=Problem);exec(compile(src,'fast-run-source.py','exec'),ns);ns['run']()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
