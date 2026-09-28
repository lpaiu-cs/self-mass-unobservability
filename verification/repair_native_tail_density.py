"""Recover the identical native density root through its electron variable."""
from pathlib import Path
from types import FunctionType
import json,sys,time
import def_native_retained_tail as owner

OUT=owner.OUT/'density-inverse';read,write,sha=owner.read,owner.write,owner.sha


def prepare():
    assert not OUT.exists();OUT.mkdir()
    original=read(owner.OUT/'plan.json');failure=read(owner.OUT/'probe-failure.json')
    assert failure['failed_state']==dict(x=-24.,T=1.,y=1e-16)
    original['budgets']['probe_seconds']-=failure['seconds']
    original['budgets']['probe_native_calls']-=failure['native_calls']
    original['bindings'].update({Path(__file__).relative_to(Path.cwd()).as_posix():sha(Path(__file__)),
        (owner.OUT/'probe-failure.json').as_posix():sha(owner.OUT/'probe-failure.json')})
    original['repair']=dict(classification='Conjectural',
        failure='Density mode2 returned104 at the first1K state. The existing fallback handles native error4 only.',
        test='Try the same native EOS in electron-variable mode0, using the already existing density Newton update and original2e-12 log-density residual. Require finite native output and all original population/H gates. No physical equation, EOS floor or tolerance change.',
        acceptance='If mode0 also fails or the same density root is not verified, preserve failure and stop; a native error is not permission to replace the EOS.',
        budget='Use only the remaining original45s/160calls probe budget; external timeout35s. Original failure is immutable.')
    write(OUT/'plan.json',original)


def install(native):
    import numpy as np
    ion=native.ion;previous=ion.gas.call;rows=[]
    def call(mode,value,t,eps,out,info):
        previous(mode,value,t,eps,out,info)
        if mode!=2 or info._obj.value!=104:return
        electron=float(eps@ion.Z)
        ne=np.exp(value)*6.02214076e23*electron
        thermal=2*(2*np.pi*9.1093837015e-28*1.380649e-16*np.exp(t)/(6.62607015e-27)**2)**1.5
        fl=np.log(ne/thermal)-(2-np.log(4.))
        for iteration in range(12):
            assert ion.calls<ion.cap,'Native density fallback budget'
            ion.calls+=1;previous(0,float(fl),t,eps,out,info)
            assert info._obj.value==0 and np.all(np.isfinite(out)),('Same native electron root',info._obj.value)
            residual=np.log(out[0])-value
            assert out[7]>0 and np.isfinite(residual),'Same native monotone density equation'
            if abs(residual)<2e-12:break
            fl-=np.clip(residual/out[7],-2.,2.)
        else:raise AssertionError('Same native electron density root did not converge')
        rows.append(dict(log_native_rho=value,logT=t,iterations=iteration+1,residual=float(residual)))
    ion.gas.call=call;native.repaired_density_rows=rows;return native


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
        write(OUT/'density-roots.json',dict(classification='Counterexample candidate',rows=[row for n in instances for row in n.repaired_density_rows],
            equation='Same native density match; electron-variable solution, original2e-12 residual.',
            original_failure_preserved=True,full_EOS_certified=False))


if __name__=='__main__':globals()[sys.argv[1]]()
