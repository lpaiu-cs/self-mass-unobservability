"""Use the identical Linux code cache within the remaining bank allocation."""
from pathlib import Path
import json,time,sys
sys.path.insert(0,'/home/lpaiu/work/native-retained-tail-runtime/verification')
import native_tail_supported_temperature as s

spent=sum([s.read(s.OUT/'cached-bank-pilot.json')['seconds'],
    s.read(s.OUT/'reallocated-bank-dispatch.json')['setup_seconds'],s.read(s.OUT/'constructor-reuse-dispatch.json')['seconds']])
cap=260-spent
assert not (s.OUT/'runtime-bank-plan.json').exists()
s.write(s.OUT/'runtime-bank-plan.json',dict(classification='Counterexample candidate',hard_seconds=cap,already_spent_seconds=spent,
    repair='Byte-identical Linux import cache; all691 source hashes verified. Native constructor reuse already passed bitwise controls. Remaining bank allocation only; no physical resolution, gate or total science-budget increase.',
    engineering_import_profile_seconds=s.read(s.owner.OUT/'runtime-import.json')['seconds'],
    prior_failed_forecasts_accepted=False,GR_seconds=10,
    bindings={str(p):s.sha(p) for p in [Path(__file__),s.owner.OUT/'runtime-source.json',s.OUT/'constructor-reuse-dispatch.json']}))
start=time.monotonic();s.owner.deadline(start,cap)
import numpy as np
n=s.fast_native(1964);z=np.load(s.OUT/'warm-probe.npz');cost=[]
for j in range(len(z['x'])):
    begin=time.monotonic();a=n.state(float(z['x'][j]),float(z['lt'][j]),float(z['y'][j]))
    assert all(np.array_equal(a[key],z[key][j]) for key in ['raw','log_fraction','affinity'])
    cost.append(time.monotonic()-begin)
elapsed=time.monotonic()-start;upper=elapsed+2*max(cost)*183+40
result=dict(classification='Counterexample candidate',bitwise=True,seconds=elapsed,per_state_seconds=cost,
    upper_seconds=upper,hard_seconds=cap,eligible=upper<cap)
s.write(s.OUT/'runtime-bank-dispatch.json',result);print(json.dumps(result),flush=True)
assert result['eligible'],'Runtime cache bank forecast'
s.bank(n,start,cap,'runtime-bank-dispatch.json')
