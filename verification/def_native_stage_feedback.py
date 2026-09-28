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


def path(steps,reference):
    folder=OUT/'linear-refined' if (steps,reference)==(128,64) else OUT
    return folder/f'steps-{steps}-reference-{reference}.npz'


def refine():
    folder=OUT/'linear-refined';assert not folder.exists();folder.mkdir()
    prior=json.loads((OUT/'result.json').read_text());assert not prior['passed'] and all(r['passed'] for r in prior['paths'][:2])
    precise=base.replace(runner,'rtol=1e-10','rtol=1e-12');precise=base.replace(precise,'info==0 and err<1e-10','info==0 and err<1e-12')
    ns=dict(vars(base.mono),OUT=folder);exec(compile(precise,__file__,'exec'),ns)
    class Precise(Response):run=ns['run']
    write(folder/'plan.json',dict(classification='Counterexample candidate',
        failure='Only reference64 stage-time path exceeded species conservation:2.803864e-8 vs1e-8. Global unweighted Krylov residual below1e-10 need not bound the differently weighted H invariant.',
        repair='Require full stage linear residual below1e-12 and GMRES rtol1e-12. Reuse both accepted reference128 paths byte-for-byte. Recompute only the failed original128/reference64 path, with unchanged equations, clocks and physical gates.',
        budgets=dict(pilot_seconds=30,production_seconds=500,CPU_threads=1,memory_GB=3),
        forecast='Failed full path cost130.03s. Measure four tightened stages, scale its observed late stepping cost by new/old pilot stage cost, add all17 operator points and15s. Require2x total under500s; no automatic increase.',
        stop='Stop on failed prefix, forecast cap, completed path gate, or500s. Preserve the failed full path; no finer clock or physical-gate relaxation.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),OUT/'result.json',OUT/'producer.py']}))
    (folder/'expanded-run.py').write_text(precise)
    start=time.monotonic();signal.signal(signal.SIGALRM,base.flow.old.optical.timeout);signal.alarm(30)
    row=Precise(64).run(128,'pilot-128',4)
    previous=json.loads((OUT/'pilot.json').read_text())['rows'][1]
    ratio=(row['stepping_seconds']/4)/(previous['stepping_seconds']/8)
    forecast=prior['paths'][2]['stepping_seconds']*ratio+17*row['operator_point_seconds']/row['operator_points']+15
    pilot=dict(classification='Counterexample candidate',row=row,forecast_seconds=forecast,upper_seconds=2*forecast,eligible=row['passed'] and 2*forecast<500,seconds=time.monotonic()-start)
    write(folder/'pilot.json',pilot);print(json.dumps(pilot),flush=True);signal.alarm(0)
    if not pilot['eligible']:return
    write(folder/'execution-plan.json',dict(classification='Counterexample candidate',budget_seconds=500,forecast_seconds=forecast,upper_seconds=2*forecast,
        bindings={str(p):sha(p) for p in [Path(__file__),folder/'plan.json',folder/'pilot.json']}))
    start=time.monotonic();signal.alarm(500)
    try:
        row=Precise(64).run(128,'steps-128-reference-64',restart='pilot-128')
        write(folder/'result.json',dict(classification='Counterexample candidate',passed=row['passed'],path=row,seconds=time.monotonic()-start))
        if not row['passed']:return
        def history(steps,ref):
            d=np.load(path(steps,ref));t=np.linspace(0,base.flow.old.END,17);ids=[np.argmin(abs(d['t']-s)) for s in t];assert np.max(abs(d['t'][ids]-t))<1e-18
            return d['moments'][ids][:,[0,1,2,3,5,6]]
        fine=history(128,128);comparisons=dict(time=compare(history(64,128),fine),background=compare(history(128,64),fine))
        result=dict(prior,passed=max(max(v) for v in comparisons.values())<.02,comparisons=comparisons,paths=prior['paths'][:2]+[row],seconds=prior['seconds']+time.monotonic()-start,
            original_stage_time_verdict=False,linear_refinement_seconds=time.monotonic()-start,accepted_paths={f'{s}/{r}':str(path(s,r)) for s,r in [(64,128),(128,128),(128,64)]})
        write(OUT/'audit.json',result);print(json.dumps(result),flush=True)
    finally:signal.alarm(0)

if __name__=='__main__':globals()[sys.argv[1]]()
