"""One fixed scale change, retaining contour and coefficient budgets."""
from pathlib import Path
import json
import signal
import resource
import time
import numpy as np
import def_gr_full_weeks as go

ROOT=go.OUT;OUT=ROOT/'beta1024';assert not OUT.exists();OUT.mkdir();go.OUT=OUT;go.BETA=1024
lu=go.splu;go.splu=lambda A:lu(A,permc_spec='NATURAL')
go.write(OUT/'plan.json',dict(classification='Counterexample candidate',
    claim='Resolve the remaining global-velocity tail through a single broader Laguerre frequency scale, without adding contour nodes or coefficients.',
    evidence='The unreduced beta512 trajectory passes scalar,old interface,new boundary and every contour check; global order1.442 is still below1.5. Local failed velocity indicators are about200-870 and the beta512 control leaves omega2000 unresolved.',
    method='One beta1024 candidate. Keep4096 contour nodes,512/1024/2048 coefficient counts,sigma12,same degree4,input,horizon and65 readouts. Extended transform and NATURAL pivoted LU unchanged. Compare with saved beta512 finest histories as an additional independent scale contrast; do not try other beta values after failure.',
    gates=dict(propagation_relative=.02,propagation_order=1.5,contour=.0002,scale_contrast=.002,linear_residual=1e-9,heat_balance=2e-13),
    budget=dict(pilot_seconds=90,first_case_seconds=450,parent_total_seconds=1200,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
    decision='Control on fixed known oscillators and16 representative measured resolvents first. Require forecast<450s. Stop any original propagation/contour gate; no automatic counts or scale expansion.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),ROOT/'p4-result.json',ROOT/'natural-pilot.json',ROOT/'control-extended.json',go.task.BANK/'fine-bank.npz']}))
signal.alarm(90)
L=go.laguerre(12);errors={}
for omega in [3.,50.,300.,600.,800.,1000.,1500.,2000.]:
    half=np.zeros((go.COUNT//2+1,1),dtype=np.clongdouble)
    for k in range(1,len(half)):
        z,f=go.contour(k,12);half[k,0]=f/(z*(z*z+omega*omega))
    half[-1]=half[-1].real;a=go.coefficients(half)[:,0];exact=(1-np.cos(omega*go.TIMES))/omega**2
    errors[str(omega)]=[float(max(abs((L[:,:n]@a[:n])[1:]-exact[1:]))*omega**2) for n in go.DEGREES]
assert max(errors[str(w)][-1] for w in [3.,50.,300.,600.,800.,1000.])<1e-6
go.write(OUT/'control.json',dict(classification='Counterexample candidate',known_oscillator_scaled_errors=errors,
    scope='Fixed known oscillators,including unresolved cases when present. No uniform all-mode claim or replacement of actual four-component gates.'))
start=time.monotonic();p=go.Problem();setup=time.monotonic()-start
ids=np.linspace(1,go.COUNT//2,16,dtype=int);start=time.monotonic();answers=np.array([p.transform(go.contour(k,12)[0]) for k in ids]);elapsed=time.monotonic()-start
np.savez_compressed(OUT/'pilot.npz',ids=ids,answers=answers)
old=json.loads((ROOT/'pilot-budget.json').read_text());forecast=1.4*(setup+elapsed/16*2048+old['inversion_forecast_seconds']+20)
go.write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,resolvent16_seconds=elapsed,
    first_case_forecast_seconds=forecast,linear_residual=p.error,assumption='Same native size and long inversion cost; actual new16 mapped points scaled to2048 with40 percent margin.'))
signal.alarm(0);print('BETA1024 FORECAST',forecast,flush=True)
assert forecast<450
signal.alarm(450);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
result=go.solve(reuse=True);contrasts={}
old=json.loads((ROOT/'p4-2048.json').read_text());new=json.loads((OUT/'p4-2048.json').read_text())
for key in go.task.FIELDS:
    a=np.array([r[key] for r in old['history']]);b=np.array([r[key] for r in new['history']])
    contrasts[key]=float(max(abs(a-b))/max(abs(b).max(),1e-100))
go.write(OUT/'scale-contrast.json',dict(classification='Counterexample candidate',differences=contrasts,scale_contrast_passed=max(contrasts.values())<.002,
    propagation_passed=result['propagation_passed'],spatial_comparison_started=False,original_failure_resolved=False))
signal.alarm(0);print('SCALE CONTRAST',json.dumps(contrasts),flush=True)
