"""Repair the saved cold native state before continuing its coupled trajectory."""
from pathlib import Path
import ctypes
import json
import signal
import sys
import time
import shutil
import subprocess
from types import FunctionType
import numpy as np
import def_native_two_way_atmosphere as prior

OUT=prior.OUT.parent/'def-native-cold-coupling'
write=prior.write;sha=prior.sha
COLD=prior.optical.ex.old.cold
CACHE=COLD.CACHE.parent/'native-cold-coupling'


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


def diagnostic_build():
    assert not CACHE.exists();CACHE.mkdir()
    write(OUT/'partition-build-plan.json',dict(classification='Counterexample candidate',
        evidence='At200K the molecular equilibrium and internal H2plus field become infinite while its external constrained field remains finite. The failure is upstream of density recovery.',
        decision='Read the exponential and partition operands at their owner; correct only their demonstrated arithmetic cause and test old supported states before using new support.',
        budget=dict(build_seconds=60,diagnostic_native_calls=8,diagnostic_seconds=20),
        forecast='The same excitation module previously compiled in3.59s. Reuse all unchanged native objects and build a separate diagnostic library; no old library is overwritten.',
        source_sha256=sha(__file__)))
    for name in ['mod_excitation.f90','excitation_pi.f90','excitation_sum.f90']:shutil.copyfile(COLD.CACHE/name,CACHE/name)
    p=CACHE/'excitation_pi.f90';src=p.read_text();anchor='                 mu(ion) = exparg*qstar(nmin_s,izqstar)';assert src.count(anchor)==1
    src=src.replace(anchor,"                 if(tl.lt.log(240._fp_kind)) write(*,*) 'COLD_EXP',iz,nmin_s,izqstar,exparg,qstar(nmin_s,izqstar),c2t*bion(ion),plop(ion),qh2plus\n"+anchor);p.write_text(src)
    build_library('diagnostic')


def build_library(label):
    old=json.loads((COLD.OUT/'stable-build.json').read_text())['receipts'];start=time.monotonic();receipts=[]
    library='free_eos_native_cold_coupling_'+label
    compile_cmd=old[0]['command'].copy();compile_cmd=[v.replace(str(COLD.CACHE),str(CACHE)).replace('mod_excitation-stable.o','mod_excitation-'+label+'.o') for v in compile_cmd]
    link_cmd=old[1]['command'].copy();link_cmd=[v.replace(str(COLD.CACHE/'libfree_eos_native_cold_stable.so'),str(CACHE/('lib'+library+'.so'))).replace('libfree_eos_native_cold_stable.so','lib'+library+'.so').replace('mod_excitation-stable.o','mod_excitation-'+label+'.o') for v in link_cmd]
    (CACHE/label).mkdir(exist_ok=True)
    bridge=COLD.old.native.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    bridge_cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(COLD.old.BUILD),str(bridge),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-l'+library,'-o',str(CACHE/label/'gas.so')]
    for i,cmd in enumerate([compile_cmd,link_cmd,bridge_cmd]):
        row=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=max(1,60-(time.monotonic()-start)))
        (OUT/f'{label}-build-{i}.log').write_text(row.stdout+row.stderr);receipts.append(dict(command=cmd,returncode=row.returncode));assert row.returncode==0,row.stderr
    write(OUT/(label+'-build.json'),dict(classification='Counterexample candidate',seconds=time.monotonic()-start,receipts=receipts,
        bindings={str(p):sha(p) for p in [CACHE/'excitation_pi.f90',CACHE/'excitation_sum.f90',CACHE/(label+'/gas.so'),CACHE/('lib'+library+'.so')]}))


def native_variant(label,cap):
    native=prior.optical.ex.Native(cap=cap)
    ion=object.__new__(prior.optical.ex.old.InventoryIons)
    init=COLD.old.Ions.__init__
    FunctionType(init.__code__,dict(init.__globals__,CACHE=CACHE/label))(ion,cap-1)
    # The physical baseline is unchanged; retain the initial constructor call
    # in resource accounting and verify the replacement surface evaluation.
    baseline=ion.snapshot(native.lr,np.log(native.fan.T),np.zeros(318))
    assert np.max(abs(baseline['eos']-native.base['eos'])/np.maximum(abs(native.base['eos']),1.))<1e-10
    native.ion=ion;native.variant_initial_calls=1
    return native


def diagnostic_native():
    signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(20);start=time.monotonic();native=native_variant('diagnostic',8)
    d=np.load(OUT/'diagnostic-states.npz');failure=None
    try:native.state(float(d['x']),np.log(200.),float(d['y']))
    except Exception as exc:failure=repr(exc)
    write(OUT/'partition-diagnostic.json',dict(classification='Counterexample candidate',failure=failure,native_calls=native.ion.calls+native.variant_initial_calls,seconds=time.monotonic()-start));signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
