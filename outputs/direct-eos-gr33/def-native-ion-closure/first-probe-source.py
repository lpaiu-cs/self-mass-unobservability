"""Counterexample candidate: constrained ions in the actual atmosphere EOS.

Reuse the external-affinity hook, but not Phase70's different atomic catalog.
No biased thermodynamic value is physical until the first-law audit passes.
"""
from pathlib import Path
from types import FunctionType
import ctypes
import difflib
import json
import re
import shutil
import signal
import subprocess
import sys
import time
import numpy as np
import def_native_vacuum_release as native
import eos_species_inventory as inventory

OUT=native.OUT.parent/'def-native-ion-closure'
CACHE=Path('/home/lpaiu/work/direct-eos-gr33/native-ion-closure')
SOURCE=native.split.model.CACHE/'source/src'
BUILD=native.split.model.CACHE/'build/src'
LIB=CACHE/'libfree_eos_native_ions.so'
write=native.write


def sha(p):
    import hashlib
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir()
    paths=[Path(__file__),Path(native.__file__),Path(inventory.__file__),
           native.split.model.LIB,native.split.BRIDGE,native.OUT/'fine.npz',
           native.OUT.parent/'def-native-metric-release/fine-1792.npz']
    # The actual saved flow filename is bound by enumerating its NPZ inputs.
    paths=[p for p in paths if p.exists()]
    paths+=sorted((native.OUT.parent/'def-native-metric-release').glob('*.npz'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='ca1099fc6',
        claim='Connect finite ion populations to the original molecular-spectral atmosphere EOS without substituting the Phase66 finite atomic catalog.',
        decision='Measure the equilibrium versus constrained-ion thermal/pressure response on saved release states. Only a thermodynamically audited constrained EOS may drive a later conservative evolution.',
        model='Conjugate constant dimensionless fields for native ion stages. Original partitions, molecules, Coulomb and photon-free gas EOS are retained. No invented relaxation rate and no opacity transfer from another EOS.',
        budget=dict(build_seconds=60,probe_seconds=90,EOS_calls=240,cpu_threads=1,memory_GB=1,new_fluid_steps=0),
        forecast='Prior related module build was a few seconds; current native calls unmeasured. First probe measures call speed. Stop on call/time/domain/identity failure; do not enlarge the budget or silently adopt biased energy.',
        gates=dict(unbiased_bitwise=True,inventory=1e-10,charge=1e-10,constrained_population=1e-8,first_law=1e-4),
        bindings={str(p):sha(p) for p in paths},
        boundary='Internal levels and molecules remain conditional equilibrium. Finite collision rates, incident spectrum, nonlinear kinetics, full GR and final charge are not certified by this probe.'))
    original=(SOURCE/'eos_calc.f90').read_text();s=original
    owner=(SOURCE/'mod_eos_calc.f90').read_text()
    anchor='  public eos_calc';assert owner.count(anchor)==1
    owner=owner.replace(anchor,anchor+'\n  real(fp_kind), save, public :: population_field(318)=0._fp_kind')
    (CACHE/'mod_eos_calc.f90').write_text(owner)
    anchor='  ! Local variables';assert anchor in s
    s=s.replace(anchor,'  real(fp_kind) :: population_dv(size(dv))\n\n'+anchor,1)
    start=s.index('  nionsp2 = ')
    s=s[:start]+'  population_dv=dv+population_field(:size(dv))\n'+re.sub(r'\bdv\b','population_dv',s[start:])
    (CACHE/'eos_calc.f90').write_text(s)
    shutil.copyfile(SOURCE/'ionize.f90',CACHE/'ionize.f90')
    (OUT/'source.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='native/eos_calc.f90',tofile='constrained/eos_calc.f90')))
    shutil.copyfile(CACHE/'mod_eos_calc.f90',OUT/'mod_eos_calc.f90')
    objects=sorted((BUILD/'CMakeFiles/free_eos.dir').glob('*.f90.o'))
    assert len(objects)>25
    objects=[str(p) for p in objects if p.name!='mod_eos_calc.f90.o']
    bridge=native.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    commands=[['gfortran','-cpp','-DUSINGDLL','-Dfree_eos_EXPORTS','-O3','-DNDEBUG','-O3','-fPIC',
               '-I'+str(CACHE),'-I'+str(BUILD),'-I'+str(SOURCE),'-c','mod_eos_calc.f90','-o','mod_eos_calc.o'],
              ['gfortran','-shared','-Wl,-soname,'+LIB.name,'mod_eos_calc.o',*objects,'-llapack','-lblas','-o',str(LIB)],
              ['gfortran','-O2','-fPIC','-shared','-I'+str(CACHE),'-I'+str(BUILD),str(bridge),
               '-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-lfree_eos_native_ions','-o','gas.so']]
    begin=time.monotonic();receipts=[]
    for i,cmd in enumerate(commands):
        result=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=max(1,60-(time.monotonic()-begin)))
        (OUT/f'build-{i}.log').write_text(result.stdout+result.stderr)
        receipts.append(dict(command=cmd,returncode=result.returncode))
        assert result.returncode==0,result.stderr
    write(OUT/'build.json',dict(seconds=time.monotonic()-begin,receipts=receipts,
        bindings={str(p):sha(p) for p in [LIB,CACHE/'gas.so',CACHE/'eos_calc.f90',*map(Path,objects)]}))
    print('BUILD',time.monotonic()-begin,flush=True)


class Ions:
    def __init__(self,cap=240):
        self.fan=native.Fan(reuse=True)
        module=native.split;gas=object.__new__(module.GasEOS);init=module.GasEOS.__init__
        FunctionType(init.__code__,dict(init.__globals__,BRIDGE=CACHE/'gas.so'),closure=init.__closure__)(gas)
        gas.inventory_lib=gas.gas_lib;call=gas.gas_lib.ionization_inventory
        def capture(mode,value,t,eps,out,info):
            raw=np.full(24,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
        gas.call=capture;self.gas=gas
        self.fields=np.ctypeslib.as_array((ctypes.c_double*318).in_dll(gas.gas_lib,'__mod_eos_calc_MOD_population_field'))
        self.fields[:]=0.;self.calls=0;self.cap=cap
        self.Z=inventory.g.d.CHARGES;self.starts=np.r_[0,np.cumsum(self.Z)]
        assert self.starts[-1]==316
        self.states=[]

    def snapshot(self,lr,lt,fields=None):
        assert self.calls<self.cap,'EOS call budget'
        self.calls+=1
        if fields is not None:self.fields[:]=fields
        a=inventory.InventoryEOS.snapshot(self.gas,float(lr),float(lt),self.fan.X)
        row=inventory.check(a,self.fan.X,float(lr),self.gas)
        assert row['inventory_error']<1e-10 and row['charge_error']<1e-10,row
        assert not row['missing_nonzero_elements'] and np.all(np.isfinite(a['eos']))
        self.states.append(dict(lrho=lr,logT=lt,fields=self.fields.copy(),**a))
        return a

    def constrain(self,lr,lt,target,guess=None):
        if guess is not None:self.fields[:]=guess
        # Chemical fields enforce stage inventories, not physical rate laws.
        # Electron degeneracy/nonideal terms are iterated by the native EOS.
        active=target>target.sum(1)[:,None]*1e-18
        for iteration in range(16):
            a=self.snapshot(lr,lt);current=a['number_fractions']
            err=float(np.max(abs(current-target))/target.sum())
            if err<1e-10:return a,iteration+1,err
            for e,z in enumerate(self.Z):
                if target[e].sum()==0:continue
                anchor=int(np.argmax(target[e]))
                ratios=np.log(np.maximum(target[e,:z+1],1e-290)/np.maximum(current[e,:z+1],1e-290))
                ratios-=ratios[anchor]
                # Fix neutral affinity to zero. All ionic affinities use the
                # same neutral reference, including a tiny neutral inventory.
                ratios-=ratios[0]
                self.fields[self.starts[e]:self.starts[e+1]]+=ratios[1:]
        raise AssertionError(('Constrained ion solve',lr,lt,err))

    def save(self,name):
        np.savez_compressed(OUT/name,**{k:np.array([s[k] for s in self.states]) for k in self.states[0]})


def probe():
    assert not (OUT/'probe.json').exists();begin=time.monotonic();signal.alarm(90)
    ion=Ions();fan=ion.fan;lr=np.log(fan.rho);lt=np.log(fan.T)
    base=ion.snapshot(lr,lt,np.zeros(318));baseline=fan.call(lr,lt)
    same=bool(np.array_equal(base['eos'],baseline));assert same,'Unbiased current native EOS replay'
    target=base['number_fractions'];rows=[]
    saved=dict(np.load(native.OUT/'fine.npz'))
    try:
        for x in [0.,-1.,-2.,-4.]:
            j=int(np.argmin(abs(saved['log_density_ratio']-x)))
            r=lr+float(saved['log_density_ratio'][j]);t=float(np.log(saved['T'][j]))
            eq=ion.snapshot(r,t,np.zeros(318))
            fixed,n,error=ion.constrain(r,t,target,np.zeros(318))
            fields=ion.fields.copy();contrasts=[]
            for dr,dt in [(1e-4,0),(-1e-4,0),(0,1e-4),(0,-1e-4)]:
                a,_,_=ion.constrain(r+dr,t+dt,target,fields)
                contrasts.append(a['eos'])
            ar,br,at,bt=contrasts
            ur=(ar[2]-br[2])/2e-4;ut=(at[2]-bt[2])/2e-4
            sr=(ar[3]-br[3])/2e-4;st=(at[3]-bt[3])/2e-4
            a=fixed['eos'];T=np.exp(t);rho=np.exp(r)
            residual=[(T*st-ut)/max(abs(ut),1.),(T*sr-ur+a[1]/rho)/max(abs(a[1]/rho),1.)]
            rows.append(dict(log_density_ratio=r-lr,T=T,iterations=n,population_error=error,
                equilibrium_H_ion=float(eq['eos'][14]),fixed_H_ion=float(a[14]),
                equilibrium_gamma=float(eq['eos'][4]),raw_fixed_pressure=float(a[1]),
                pressure_ratio=float(a[1]/eq['eos'][1]),raw_fixed_u=float(a[2]),
                fixed_cvT=float(ut),equilibrium_cvT=float(eq['eos'][10]),
                raw_entropy_first_law_residual=list(map(float,residual)),
                molecular_H_fractions=fixed['molecular_H_fractions'].tolist(),max_affinity=float(max(abs(fields)))))
            write(OUT/'probe-progress.json',dict(rows=rows,EOS_calls=ion.calls,seconds=time.monotonic()-begin))
            print('PROBE',rows[-1],flush=True)
        ion.fields[:]=0.;replay=ion.snapshot(lr,lt)
        assert all(np.array_equal(base[k],replay[k]) for k in base)
        result=dict(classification='Counterexample candidate',unbiased_bitwise=same,history_bitwise=True,
            EOS_calls=ion.calls,seconds=time.monotonic()-begin,rows=rows,
            raw_thermodynamic_values_accepted=bool(max(abs(v) for r in rows for v in r['raw_entropy_first_law_residual'])<1e-4),
            full_kinetics=False,new_fluid_steps=0,full_goal_complete=False)
        write(OUT/'probe.json',result);print('PROBE COMPLETE',result['EOS_calls'],result['seconds'],flush=True)
    except Exception as exc:
        write(OUT/'probe-failure.json',dict(error=repr(exc),EOS_calls=ion.calls,seconds=time.monotonic()-begin,rows=rows))
        raise
    finally:ion.save('probe-states.npz')
    signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
