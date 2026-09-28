"""Counterexample candidate: stage-time nonautonomous photon/material return.

Preserve the rejected midpoint-frozen waveform sweep. Reuse its operator and
change only when SDIRK evaluates the actual coefficients and paired sources.
"""
from pathlib import Path
import json,signal,sys,time
import numpy as np
import sympy as sp
import def_native_matter_photon_feedback as base

OUT=base.OUT/'stage-time';write=base.write;sha=base.sha
runner=base.runner
runner=base.replace(runner,'t=k*h;c=self.local(t+h/2);inverse=self.inverse(c,gamma*h)',
    't=k*h;c=self.local(t+gamma*h);mechanical=c["mechanical"];inverse=self.inverse(c,gamma*h)')
runner=base.replace(runner,'s2,l2,e2=self.source(t+h);',
    'c=self.local(t+h);c["mechanical"]=mechanical;inverse=self.inverse(c,gamma*h)\n        qg=self.gas(c["q"],c["qb"],c["qe"])+mechanical\n        s2,l2,e2=self.source(t+h);')
# All macro steps end on or inside the prescribed 17-knot intervals. At a
# knot use the derivative from this step's left interval, not the next one.
namespace=dict(vars(base.mono),OUT=OUT);exec(compile(runner,__file__,'exec'),namespace)
class Response(base.Response):
    run=namespace['run']


def prepare():
    assert not OUT.exists();OUT.mkdir();prior=json.loads((base.OUT/'result.json').read_text());assert not prior['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d6ccddff9',
        failure='Full midpoint-frozen return passed invariants but neutral H64/128 difference was6.6447percent, above the unchanged2percent gate.',
        repair='Evaluate the actual nonautonomous collision operator and affine drive at each SDIRK stage time. The old midpoint freeze has an O(h) displaced equilibrium in the stiff limit. Mechanical interval derivatives retain their left value at aligned endpoint knots; paired transfers and ledgers use the same stages.',
        reuse='Same531 cells,8 angles,152 frequencies,17 input knots,64/128 clocks and3.434ms. Reuse accepted background/EOS/material histories. Preserve all failed midpoint paths. No finer clocks or native EOS replay.',
        decision='Accept only if the same three full paths pass2percent time/background and existing conservation/linear gates, then return collision histories to matter and GR.',
        budget=dict(pilot_seconds=60,production_seconds=900,CPU_threads=1,memory_GB=3),
        forecast='Measure equal-horizon4/8-step prefixes and17-point preparation. Retain2x forecast within900s; full failed run cost343.37s provides a late-cost comparison.',
        stop='Stop on failed pilot or path, cap, or unchanged2percent comparison. No automatic finer clock, longer horizon or gate relaxation.',
        limits='A single waveform return on stored GR and prescribed baryon/momentum; not a coupled fixed point or final charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),base.OUT/'result.json',base.OUT/'execution-plan.json']}))
    (OUT/'expanded-run.py').write_text(runner)
    h,L,t,y=sp.symbols('h L t y',positive=True);g=1-1/sp.sqrt(2)
    def step(a,b):
        stage=(y+g*h*L*a)/(1+g*h*L)
        return (y+(1-g)*h*L*(a-stage)+g*h*L*b)/(1+g*h*L)
    midpoint=sp.simplify(sp.limit(step(t+h/2,t+h/2),L,sp.oo))
    actual=sp.simplify(sp.limit(step(t+g*h,t+h),L,sp.oo))
    assert sp.simplify(midpoint-t-h/2)==0 and sp.simplify(actual-t-h)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Scalar stiff control y_prime=L*(t-y), L to infinity. Frozen midpoint SDIRK tends to t+h/2 while actual stage-time SDIRK tends to t+h. This identifies an algorithmic defect, not by itself the entire observed full-system error.',midpoint=str(midpoint),stage_time=str(actual)))


def compare(a,b):
    return (np.max(np.sum(abs(a-b),axis=2),axis=0)/np.maximum(np.max(np.sum(abs(b),axis=2),axis=0),1.)).tolist()


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,base.flow.old.optical.timeout);signal.alarm(60)
    rows=[Response(128).run(64,'pilot-64',4),Response(128).run(128,'pilot-128',8)]
    a=np.load(OUT/'pilot-64.npz')['moments'][-1:];b=np.load(OUT/'pilot-128.npz')['moments'][-1:];dif=compare(a[:,[0,1,2,3,5,6]],b[:,[0,1,2,3,5,6]])
    point=max(r['operator_point_seconds']/r['operator_points'] for r in rows)
    forecast=point*51+rows[0]['stepping_seconds']/4*60+rows[1]['stepping_seconds']/8*248+18
    result=dict(classification='Counterexample candidate',rows=rows,equal_horizon_comparison=dif,forecast_seconds=forecast,upper_seconds=2*forecast,eligible=all(r['passed'] for r in rows) and max(dif)<.02 and 2*forecast<900,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);print(json.dumps(result),flush=True);signal.alarm(0)
    if result['eligible']:
        write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=900,forecast_seconds=forecast,upper_seconds=2*forecast,
            paths=[[64,128,'pilot-64'],[128,128,'pilot-128'],[128,64,None]],
            bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),OUT/'plan.json',OUT/'symbolic.json',OUT/'pilot.json',base.matter.OUT/'source-audit.json']}))


from types import FunctionType
production=FunctionType(base.production.__code__,dict(vars(base),OUT=OUT,Response=Response,__file__=__file__))

if __name__=='__main__':globals()[sys.argv[1]]()
