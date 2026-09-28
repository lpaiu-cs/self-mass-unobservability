"""Repair the saved cold native state before continuing its coupled trajectory."""
from pathlib import Path
import ctypes
import json
import signal
import sys
import time
import numpy as np
import def_native_two_way_atmosphere as prior

OUT=prior.OUT.parent/'def-native-cold-coupling'
write=prior.write;sha=prior.sha


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1297eba71',
        claim='Repair the exact native EOS state blocking the real two-way photon/moving-atmosphere fine trajectory, then apply the validated input to its remaining original interval.',
        reuse='Original448/64 and448/128 completed paths, fine104 complete snapshot and fine108 history, original918 native states and the saved failed conservative cell. No full-path replay or added spatial/time grid.',
        first_decision='Use exposed native molecular diagnostics and constrained states to determine why finite requested density underflows internally at160K. Preserve failed attempts. No ideal-EOS substitution or temperature clamp.',
        first_budget=dict(native_calls=60,seconds=20,CPU_threads=1,memory_GB=1,fluid_steps=0),
        forecast='Phase112 seven native endpoint calls took1.75s; this bounded single-cell diagnostic is expected1-10s. If native/build/support changes are needed, assess and register those before execution.',
        gates=dict(constitutive=.002,population=1e-12,relative_H=1e-8,primitive=2e-11,energy=1e-8,space_trace=.02,space_mass=.02),
        stop='Do not enlarge the horizon, add full flow paths, relax gates or silently extrapolate an EOS or spectrum. A recovered state must actually enter the original coupled evolution before calling its blocking condition solved.',
        bindings={str(p):sha(p) for p in [Path(__file__),prior.OUT/'root-failure-state.npz',prior.OUT/'cells-896-steps-128.npz',prior.OUT/'result.json',prior.OUT/'spectrum-bank.npz',prior.optical.ex.OUT/'repaired-bank.npz']}))


def diagnose():
    assert not (OUT/'diagnostic.json').exists();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(20);start=time.monotonic()
    native=prior.optical.ex.Native(cap=60);f=prior.Flow(896);m=f.base;U=np.load(prior.OUT/'root-failure-state.npz')['U'][:,416]
    D=U[0];tau=(U[2]-(m.a[416]-m.a0)*f.eos.cx*D)/m.a[416];v=U[1]/(f.eos.cx*D+tau);rho=D*np.sqrt(1-v*v);x=float(np.log(rho));y=float(U[3]/D)
    rows=[];fields=[];lib=native.ion.gas.gas_lib
    def array(name,n):return np.ctypeslib.as_array((ctypes.c_double*n).in_dll(lib,'__mod_nuvar_MOD_'+name)).copy().tolist()
    try:
        for T in [240.,200.,170.,160.]:
            failure=None
            try:state=native.state(x,float(np.log(T)),y)
            except Exception as exc:failure=repr(exc)
            row=dict(T=T,failure=failure,native_calls=native.ion.calls,molecular_log_populations=array('mol_logs',3),molecular_equilibria=array('mol_eq',2),molecular_fields=array('mol_dv',2),fields_H_H2_H2plus=native.ion.fields[[0,316,317]].tolist())
            if failure is None:row.update(pressure=float(state['raw'][1]),energy=float(state['raw'][2]),population_error=state['population_error'])
            rows.append(row);fields.append(native.ion.fields.copy());print(json.dumps(row),flush=True)
    finally:
        np.savez_compressed(OUT/'diagnostic-states.npz',fields=fields,temperatures=[r['T'] for r in rows],x=x,y=y,rho=rho,U=U)
        write(OUT/'diagnostic.json',dict(classification='Counterexample candidate',checks=rows,native_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=sha(__file__)))
        signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
