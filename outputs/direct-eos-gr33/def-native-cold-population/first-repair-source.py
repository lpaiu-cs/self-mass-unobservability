"""Counterexample candidate: repair the exact stopped constrained EOS state."""
from pathlib import Path
from types import FunctionType
import ctypes
import json
import re
import signal
import subprocess
import sys
import time
import numpy as np
import def_native_ion_closure as old

OUT=old.OUT.parent/'def-native-cold-population'
CACHE=old.CACHE.parent/'native-cold-population'
write=old.write;sha=old.sha


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='14e7524ef',
        claim='Repair the precise786K fixed-inventory native-EOS failure, then complete the original49node expansion without changing its free energy, populations, domain or scientific gates.',
        decision='Only a successful same-state recovery with native inventory/thermodynamic checks can replace the blocked EOS input. Retain the original LTE charge as conditional until actual chemistry and fluid are evolved together.',
        reuse='Original libraries, saved34node prefix, exact failed next state, original density nodes and entropy root. No new stellar/fluid integration while its EOS is not evaluable.',
        budget=dict(native_build_seconds=60,diagnostic_calls=2,diagnostic_seconds=20,repair_EOS_evaluations=1000,repair_seconds=90,CPU_threads=1,memory_GB=1),
        forecast='Previous native build3.35s,499 constrained state evaluations15.15s including failed inversions. Frozen-population strong fields may need a different initializer or evaluation route; its cost is unmeasured and bounded by the caps.',
        gates=dict(unbiased_bitwise=True,old_supported_state_relative=1e-10,inventory=1e-10,constrained_population=1e-12,entropy_root=2e-10,first_law=1e-4,coarse_fine_speed=.002),
        stop='No relaxation rate, EOS/temperature extension, hidden bound population change, relaxed gate or long fluid job. Diagnose the exact native equations and preserve each executed version.',
        inputs={str(p):sha(p) for p in [Path(__file__),old.LIB,old.CACHE/'gas.so',old.OUT/'accepted-adiabat-prefix.npz',old.OUT/'adiabat-failure.json',old.OUT/'final-adiabat-source.py']}))
    original=(old.native.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90').read_text()
    bridge=original.replace('call free_eos(0,3,11,-2,kif','call free_eos(2,3,11,-2,kif')
    assert bridge!=original
    (OUT/'diagnostic-bridge.f90').write_text(bridge)
    cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(old.BUILD),str(OUT/'diagnostic-bridge.f90'),
         '-L'+str(old.CACHE),'-Wl,-rpath,'+str(old.CACHE),'-lfree_eos_native_ions','-o',str(CACHE/'gas.so')]
    start=time.monotonic();p=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
    (OUT/'diagnostic-build.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    write(OUT/'diagnostic-build.json',dict(command=cmd,seconds=time.monotonic()-start,sha256=sha(CACHE/'gas.so')))
    d=np.load(old.OUT/'accepted-adiabat-prefix.npz');fields=d['fields'][-1].copy()
    ion=old.Ions(cap=2);lr0=np.log(ion.fan.rho);x=float(d['log_density_ratio'][-1]-.125)
    lt0=np.log(d['T'][-1]);lt=lt0-2*.125/3
    raw=(old.SOURCE/'mod_ionization_data.f90').read_text();body=raw[raw.index('monatomic_ip(nions) ='):];body=body[:body.index(']')]
    body='\n'.join(line.split('!')[0] for line in body.splitlines())
    potentials=np.array([float(v.replace('d','e')) for v in re.findall(r'([-+]?[0-9]+\.[0-9]*(?:[edED][-+]?[0-9]+)?)_fp_kind',body)])
    assert len(potentials)==316
    c2=float(json.loads((old.native.OUT.parent/'def-photon-shared-atomic/catalog.json').read_text())['constants'][2])
    energy=np.concatenate([np.cumsum(potentials[ion.starts[e]:ion.starts[e+1]]) for e in range(24)])*c2
    charges=np.concatenate([np.arange(1,z+1) for z in ion.Z]);target=d['number_fractions'][0]
    active=np.concatenate([target[e,1:z+1]>target[e].sum()*1e-18 for e,z in enumerate(ion.Z)])
    dt=lt-lt0;dr=-.125;dinv=np.exp(-lt)-np.exp(-lt0)
    fields[:316]+=np.where(active,energy*dinv+charges*(dr-1.5*dt),0.)
    diss=float(re.search(r'h2diss = ([0-9.]+)_fp_kind',raw).group(1))*c2
    fields[316]+=-diss*dinv-dr+1.5*dt
    fields[317]+=(diss+potentials[0]*c2-21375.95*c2)*dinv+dr-1.5*dt
    np.savez_compressed(OUT/'failed-state.npz',lrho=lr0+x,logT=lt,fields=fields,X=ion.fan.X,target=target,
        last_fields=d['fields'][-1],last_lrho=lr0+d['log_density_ratio'][-1],last_logT=lt0)
    with (OUT/'diagnostic.log').open('w') as f:
        p=subprocess.run([sys.executable,__file__,'diagnose'],stdout=f,stderr=subprocess.STDOUT,text=True,timeout=20)
    write(OUT/'diagnostic-receipt.json',dict(returncode=p.returncode,log_sha256=sha(OUT/'diagnostic.log')))
    assert p.returncode==0
    print('PREPARED and reproduced exact cold failure',flush=True)


def diagnose(bridge='gas.so',label='diagnostic'):
    module=old.native.split;gas=object.__new__(module.GasEOS);init=module.GasEOS.__init__
    FunctionType(init.__code__,dict(init.__globals__,BRIDGE=CACHE/bridge),closure=init.__closure__)(gas)
    fields=np.ctypeslib.as_array((ctypes.c_double*318).in_dll(gas.gas_lib,'__mod_eos_calc_MOD_population_field'))
    d=np.load(OUT/'failed-state.npz');fields[:]=d['fields']
    try:gas(2,float(d['lrho']),float(d['logT']),d['X'])
    except AssertionError as exc:
        write(OUT/(label+'.json'),dict(classification='Counterexample candidate',reproduced=True,error=repr(exc),T=float(np.exp(d['logT'])),source_sha256=sha(__file__)))
    else:raise AssertionError('Expected exact original density failure did not reproduce')


def jacobian():
    import shutil
    start=time.monotonic();assert not (OUT/'jacobian-build.json').exists()
    owner=(old.SOURCE/'mod_free_eos_detailed.f90').read_text();(CACHE/'mod_free_eos_detailed.f90').write_text(owner)
    original=(old.SOURCE/'free_eos_detailed.f90').read_text()
    anchor='        jacobian_save(1:njacobian,1:njacobian) = jacobian(1:njacobian, 1:njacobian)'
    assert original.count(anchor)==1
    added='''
        if(t.lt.1000._fp_kind.and.ioncount.eq.1) then
           write(stderr,*) 'COLD_MATRIX',njacobian,n_partial_aux
           do ijacobian=1,njacobian
              write(stderr,*) 'ROW',ijacobian,jacobian(ijacobian,1:njacobian),'RHS',rhs1(ijacobian)
           enddo
           do ijacobian=1,n_partial_aux
              iaux=partial_aux(ijacobian)
              write(stderr,*) 'AUX',ijacobian,iaux,aux_old(iaux),aux(iaux)
           enddo
        endif'''
    (CACHE/'free_eos_detailed.f90').write_text(original.replace(anchor,anchor+added))
    for name in ['eos_cold_start.f90','eos_tqft.f90']:shutil.copyfile(old.SOURCE/name,CACHE/name)
    objects=[p for p in (old.BUILD/'CMakeFiles/free_eos.dir').glob('*.f90.o') if p.name not in ['mod_eos_calc.f90.o','mod_free_eos_detailed.f90.o']]
    library=CACHE/'libfree_eos_native_cold_debug.so'
    commands=[['gfortran','-cpp','-DUSINGDLL','-Dfree_eos_EXPORTS','-O3','-DNDEBUG','-fPIC','-I'+str(CACHE),'-I'+str(old.BUILD),'-I'+str(old.SOURCE),'-c','mod_free_eos_detailed.f90','-o','mod_free_eos_detailed.o'],
        ['gfortran','-shared','-Wl,-soname,'+library.name,str(old.CACHE/'mod_eos_calc.o'),'mod_free_eos_detailed.o',*map(str,objects),'-llapack','-lblas','-o',str(library)],
        ['gfortran','-O2','-fPIC','-shared','-I'+str(old.BUILD),str(OUT/'diagnostic-bridge.f90'),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-lfree_eos_native_cold_debug','-o','gas-debug.so']]
    receipts=[]
    for i,cmd in enumerate(commands):
        p=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=30)
        (OUT/f'jacobian-build-{i}.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
        receipts.append(dict(command=cmd,returncode=p.returncode))
    write(OUT/'jacobian-build.json',dict(seconds=time.monotonic()-start,receipts=receipts,source_sha256=sha(CACHE/'free_eos_detailed.f90')))
    with (OUT/'jacobian.log').open('w') as f:
        p=subprocess.run([sys.executable,__file__,'inspect'],stdout=f,stderr=subprocess.STDOUT,text=True,timeout=20)
    assert p.returncode==0
    print('JACOBIAN captured',flush=True)


def inspect():diagnose('gas-debug.so','jacobian')


def repair():
    import difflib
    start=time.monotonic();assert not (OUT/'repair-build.json').exists();(CACHE/'fixed').mkdir()
    source=(old.SOURCE/'excitation_pi.f90').read_text()
    before='mu(ion)*qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)'
    after='mu(ion)*(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))'
    assert source.count(before)==1;fixed=source.replace(before,after)
    line='                       mu_aux(4,ion) = -occ_const*mu(ion)*r_ion3(ion) + '+after
    assert fixed.count(line)==1
    # Operand logging is restricted to the separate verbose diagnostic bridge.
    fixed=fixed.replace(line,"                       if(verbosity.ge.2.and.tl.lt.log(1000._fp_kind)) &\n                            write(stderr,*) 'SCALED_PRODUCT',ion,mu(ion),qstarx(5,nmin_s,izqstar),qstar(nmin_s,izqstar)\n"+line)
    (CACHE/'excitation_pi.f90').write_text(fixed)
    (CACHE/'mod_excitation.f90').write_text((old.SOURCE/'mod_excitation.f90').read_text())
    (OUT/'repair.patch').write_text(''.join(difflib.unified_diff(source.splitlines(True),fixed.splitlines(True),fromfile='native/excitation_pi.f90',tofile='stable/excitation_pi.f90')))
    objects=[p for p in (old.BUILD/'CMakeFiles/free_eos.dir').glob('*.f90.o') if p.name not in ['mod_eos_calc.f90.o','mod_excitation.f90.o']]
    library=CACHE/'libfree_eos_native_cold_fixed.so';bridge=old.native.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    commands=[['gfortran','-cpp','-DUSINGDLL','-Dfree_eos_EXPORTS','-O3','-DNDEBUG','-O3','-fPIC','-I'+str(CACHE),'-I'+str(old.BUILD),'-I'+str(old.SOURCE),'-c','mod_excitation.f90','-o','mod_excitation-fixed.o'],
        ['gfortran','-shared','-Wl,-soname,'+library.name,str(old.CACHE/'mod_eos_calc.o'),'mod_excitation-fixed.o',*map(str,objects),'-llapack','-lblas','-o',str(library)]]
    for path,target in [(bridge,CACHE/'fixed/gas.so'),(OUT/'diagnostic-bridge.f90',CACHE/'gas-fixed-verbose.so')]:
        commands.append(['gfortran','-O2','-fPIC','-shared','-I'+str(old.BUILD),str(path),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-lfree_eos_native_cold_fixed','-o',str(target)])
    receipts=[]
    for i,cmd in enumerate(commands):
        p=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=30)
        (OUT/f'repair-build-{i}.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
        receipts.append(dict(command=cmd,returncode=p.returncode))
    write(OUT/'repair-build.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-start,receipts=receipts,
        source_sha256=sha(CACHE/'excitation_pi.f90'),library_sha256=sha(library),bridge_sha256=sha(CACHE/'fixed/gas.so'),
        intervention='Parenthesize the molecular mu*qstarx/qstar derivative like the already-stable sibling atomic derivatives. Do not alter partitions, affinities, derivatives or solver thresholds.'))
    with (OUT/'repaired-state.log').open('w') as f:
        p=subprocess.run([sys.executable,__file__,'recovered'],stdout=f,stderr=subprocess.STDOUT,text=True,timeout=20)
    assert p.returncode==0,(OUT/'repaired-state.log').read_text()[-2000:]
    print('REPAIRED exact failed state',flush=True)


class Ions(old.Ions):
    __init__=FunctionType(old.Ions.__init__.__code__,dict(old.Ions.__init__.__globals__,CACHE=CACHE/'fixed'))
    save=FunctionType(old.Ions.save.__code__,dict(old.Ions.save.__globals__,OUT=OUT))


def recovered():
    module=old.native.split;gas=object.__new__(module.GasEOS);init=module.GasEOS.__init__
    FunctionType(init.__code__,dict(init.__globals__,BRIDGE=CACHE/'gas-fixed-verbose.so'),closure=init.__closure__)(gas)
    fields=np.ctypeslib.as_array((ctypes.c_double*318).in_dll(gas.gas_lib,'__mod_eos_calc_MOD_population_field'))
    d=np.load(OUT/'failed-state.npz');fields[:]=d['fields']
    a=gas(2,float(d['lrho']),float(d['logT']),d['X']);assert np.all(np.isfinite(a))
    np.savez_compressed(OUT/'repaired-state.npz',eos=a,**dict(d))
    write(OUT/'repaired-state.json',dict(classification='Counterexample candidate',same_failed_input_evaluated=True,
        T=float(np.exp(d['logT'])),rho=float(a[0]),P=float(a[1]),native_info=0,source_sha256=sha(__file__)))


if __name__=='__main__':globals()[sys.argv[1]]()
