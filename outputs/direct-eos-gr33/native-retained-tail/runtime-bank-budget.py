"""Finish the same bank with a bounded advance from unused production time."""
from pathlib import Path
import sys,time
sys.path.insert(0,'/home/lpaiu/work/native-retained-tail-runtime/verification')
import native_tail_supported_temperature as s

assert not (s.OUT/'aggregate-budget-plan.json').exists()
old=s.read(s.OUT/'runtime-bank-plan.json');probe=s.read(s.OUT/'runtime-bank-dispatch.json')
spent=old['already_spent_seconds']+probe['seconds'];cap=340-spent
s.write(s.OUT/'aggregate-budget-plan.json',dict(classification='Counterexample candidate',
    correction='Earlier aggregate prose mistakenly added the original stage caps to1195s; the correct sum is1095s. Use1095s, not the erroneous larger sum.',
    original_total_seconds=1095,maximum_total_seconds=1095,bank_cap_seconds=cap,bank_already_spent_seconds=spent,
    advance_from_unstarted_production_seconds=80,maximum_production_seconds=650,
    decision='Advance up to80s from the unused production allocation to finish the actual native table after the verified Linux cache reduced setup to8.51s. Production may use at most650s and must also fit the aggregate remainder after all actual bank/control/pilot costs and reserved10s GR. This is a budget reallocation, not new resolution or a larger total allocation.',
    forecast='Reuse the measured identical-cache state maximum, doubled plus40s; add the actual setup before dispatch. Stop if this fails. Preserve every earlier rejected forecast.',
    bindings={str(p):s.sha(p) for p in [Path(__file__),s.OUT/'runtime-bank-plan.json',s.OUT/'runtime-bank-dispatch.json',s.owner.OUT/'runtime-source.json']}))
start=time.monotonic();s.owner.deadline(start,cap)
n=s.fast_native(1947);setup=time.monotonic()-start;upper=setup+2*max(probe['per_state_seconds'])*183+40
result=dict(eligible=upper<cap,setup_seconds=setup,upper_seconds=upper,hard_seconds=cap)
s.write(s.OUT/'aggregate-bank-dispatch.json',result);print(result,flush=True)
assert result['eligible'],'Aggregate bank budget'
s.bank(n,start,cap,'aggregate-bank-dispatch.json')
