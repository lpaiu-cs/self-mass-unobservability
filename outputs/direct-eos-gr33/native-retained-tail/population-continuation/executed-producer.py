"""Cold native population continuation; retain equations and all residual gates."""
from pathlib import Path
from types import FunctionType
import sys
import def_native_retained_tail as owner

OUT=owner.OUT/'population-continuation';read,write,sha=owner.read,owner.write,owner.sha


def prepare():
    assert not OUT.exists();OUT.mkdir()
    plan=read(owner.OUT/'plan.json')
    failures=[owner.OUT/'probe-failure.json',owner.OUT/'density-inverse/probe-failure.json']
    for p in failures:
        row=read(p);plan['budgets']['probe_seconds']-=row['seconds'];plan['budgets']['probe_native_calls']-=row['native_calls']
        plan['bindings'][p.as_posix()]=sha(p)
    plan['bindings'][Path(__file__).relative_to(Path.cwd()).as_posix()]=sha(Path(__file__))
    source=Path('/home/lpaiu/work/direct-eos-gr33/molecular-spectral/source/src')
    code=(source/'eos_calc.f90').read_text();codes=(source/'mod_info_data.f90').read_text()
    assert 'info_offset_eos_calc = 100' in codes
    guard=code[code.index('     ne = hne + h_ion + h2plus'):];guard=guard[:guard.index('     if(ifnr03)')]
    assert 'info = info_offset_eos_calc + 4' in guard and 'ne.lt.underflow_limit' in guard
    (OUT/'native-error-104.txt').write_text(guard)
    plan['repair']=dict(classification='Conjectural',
        correction='104 is the native electron-inventory underflow guard (100+4), not density Newton error4. Both direct density and electron-variable initial evaluations failed.',
        method='For a requested cold state, start with the same constrained native state at80K and follow converged affinity seeds downward. Each next temperature is at least0.7 times the previous and its inverse-temperature increment is at most0.1 per kelvin. These are solver seeds, not new physical states or a temperature clamp.',
        invariants='Preserve original population1e-12 and relativeH1e-8 gates, chemistry, partitions, density target and final requested temperature. Cache successful affinity seeds in process memory only. The original prefix and libraries are immutable.',
        stop='No extra probe budget: retain both failures and use only remaining45s/160calls. Stop on another failure, insufficient call/time budget or failed warm controls.',
        source_hashes={str(source/p):sha(source/p) for p in ['eos_calc.f90','mod_info_data.f90']})
    write(OUT/'plan.json',plan)


def install(native):
    import numpy as np
    base=native.state;d={key:native.prefix[key].copy() for key in ['T','log_density_ratio','fields']};native.prefix=d
    rows=[];lowest={}
    def sample(x,T,y):
        state=base(x,float(np.log(T)),y)
        fields=native.ion.fields.copy()
        # The original owner adds this H offset relative to y0. Undo it in
        # the saved seed so the next call applies the requested offset once.
        fields[0]-=np.log((1-y)/y)-np.log((1-native.y0)/native.y0)
        d['T']=np.r_[d['T'],T];d['log_density_ratio']=np.r_[d['log_density_ratio'],x];d['fields']=np.vstack([d['fields'],fields])
        rows.append(dict(x=x,T=T,y=y,population_error=state['population_error'],native_calls=native.ion.calls+native.variant_initial_calls))
        lowest[x,y]=min(lowest.get((x,y),80.),T)
        return state
    def state(x,lt,y):
        target=float(np.exp(lt))
        if target<80:
            key=(x,y)
            if key not in lowest:sample(x,80.,y)
            current=lowest[key]
            while current>target*(1+1e-14):
                next_T=max(target,.7*current,1/(1/current+.1))
                value=sample(x,next_T,y);current=next_T
            if abs(current-target)<1e-13*target:return value
        return sample(x,target,y)
    native.state=state;native.continuation_rows=rows;return native


def run():
    import def_native_cold_coupling as cold
    plan=read(OUT/'plan.json');remaining=plan['budgets']['probe_seconds'];original=cold.logarithmic_native;instances=[]
    def build(cap):
        n=install(original(cap));instances.append(n);return n
    cold.logarithmic_native=build
    def deadline(start,cap):return owner.deadline(start,min(cap,remaining))
    probe=FunctionType(owner.probe.__code__,dict(owner.probe.__globals__,OUT=OUT,deadline=deadline))
    try:probe()
    finally:
        cold.logarithmic_native=original
        write(OUT/'continuation.json',dict(classification='Counterexample candidate',rows=[r for n in instances for r in n.continuation_rows],
            original_failures_preserved=True,equations_and_gates_unchanged=True,physical_tail_evolved=False))


if __name__=='__main__':globals()[sys.argv[1]]()
