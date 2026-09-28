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
        self.density_fallbacks=[]
        def capture(mode,value,t,eps,out,info):
            raw=np.full(24,np.nan);call(mode,value,t,eps,raw,info)
            if mode==2 and info._obj.value==4:
                # Native density Newton can start from an effectively neutral
                # cold gas despite imposed ion affinities. Solve the identical
                # density equation with its exposed electron variable instead.
                electron=float(eps@inventory.g.d.CHARGES)
                ne=np.exp(value)*6.02214076e23*electron
                thermal=2*(2*np.pi*9.1093837015e-28*1.380649e-16*np.exp(t)/(6.62607015e-27)**2)**1.5
                fl=np.log(ne/thermal)-(2-np.log(4.))
                for iteration in range(12):
                    assert self.calls<self.cap,'EOS density-fallback budget'
                    self.calls+=1
                    call(0,float(fl),t,eps,raw,info)
                    assert info._obj.value==0,('Electron-variable native solve',info._obj.value)
                    residual=np.log(raw[0])-value
                    assert raw[7]>0 and np.isfinite(residual),'Native density monotonicity'
                    if abs(residual)<2e-12:break
                    fl-=np.clip(residual/raw[7],-2.,2.)
                else:raise AssertionError('Electron-variable density root')
                self.density_fallbacks.append(dict(log_native_rho=value,logT=t,iterations=iteration+1,residual=float(residual)))
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

    def constrain(self,lr,lt,target,guess=None,target_molecules=None,tolerance=1e-10):
        if guess is not None:self.fields[:]=guess
        # Chemical fields enforce stage inventories, not physical rate laws.
        # Electron degeneracy/nonideal terms are iterated by the native EOS.
        active=target>target.sum(1)[:,None]*1e-18
        for iteration in range(16):
            a=self.snapshot(lr,lt);current=a['number_fractions']
            err=float(np.max(abs(current-target))/target.sum())
            if err<tolerance:return a,iteration+1,err
            for e,z in enumerate(self.Z):
                if target[e].sum()==0:continue
                anchor=int(np.argmax(target[e]))
                ratios=np.log(np.maximum(target[e,:z+1],1e-290)/np.maximum(current[e,:z+1],1e-290))
                ratios-=ratios[anchor]
                # Fix neutral affinity to zero. All ionic affinities use the
                # same neutral reference, including a tiny neutral inventory.
                ratios-=ratios[0]
                # Do not invert native underflow zeros into finite ions. Only
                # enforce initially populated coordinates; the remaining
                # coordinates equilibrate and their total is checked below.
                take=active[e,1:z+1]
                self.fields[self.starts[e]:self.starts[e+1]]+=np.where(take,ratios[1:],0.)
            if target_molecules is not None:
                now=a['molecular_H_fractions']
                correction=np.log(np.maximum(target_molecules,1e-290)/np.maximum(now,1e-290))
                correction-=2*np.log(target[0,0]/max(current[0,0],1e-290))
                self.fields[316:]+=np.where(now>1e-18,correction,0.)
        raise AssertionError(('Constrained ion solve',lr,lt,err))

    def save(self,name):
        np.savez_compressed(OUT/name,**{k:np.array([s[k] for s in self.states]) for k in self.states[0]})


def probe():
    assert not (OUT/'probe.json').exists();begin=time.monotonic();signal.alarm(90)
    previous=json.loads((OUT/'first-probe-failure.json').read_text()) if (OUT/'first-probe-failure.json').exists() else dict(EOS_calls=0,seconds=0)
    write(OUT/'probe-repair-plan.json',dict(classification='Counterexample candidate',
        failure='Attempting to invert underflow-zero trace ions amplified populations with initial elemental fractions below1e-18. The saved nitrogen trace cycles caused the full-vector convergence failure.',
        repair='Constrain populated coordinates only; retain native conditional equilibrium of smaller coordinates and check the original total population error. No EOS equation or numerical acceptance gate is relaxed.',
        source_sha256=sha(__file__),prior_calls=previous['EOS_calls'],remaining_calls=240-previous['EOS_calls'],
        remaining_seconds=90-previous['seconds'],minimum_initial_stage_fraction=1e-18))
    signal.alarm(max(1,int(90-previous['seconds'])))
    ion=Ions(cap=240-previous['EOS_calls']);fan=ion.fan;lr=np.log(fan.rho);lt=np.log(fan.T)
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


def adiabat():
    from scipy.interpolate import PchipInterpolator
    from scipy.integrate import cumulative_simpson
    assert not (OUT/'adiabat.json').exists()
    assert not (OUT/'adiabat-failure.json').exists(),'Registered cold closure stopped; a new physical/numerical design is required before another attempt'
    failures=[json.loads(p.read_text()) for p in OUT.glob('*-adiabat-failure.json')]
    previous={k:sum(d[k] for d in failures) for k in ['EOS_calls','seconds']}
    begin=time.monotonic();signal.alarm(max(1,int(60-previous['seconds'])))
    write(OUT/'adiabat-plan.json',dict(classification='Counterexample candidate',
        claim='Solve the same-native nonlinear adiabatic expansion at fixed initially populated ion-stage inventories. This is the finite-reaction limiting constitutive branch required by a population-conserving fluid solver.',
        decision='Supply a physically explicit alternative to LTE and quantify its temperature, pressure and characteristic-speed change; do not relabel a fixed-density replay as a new self-consistent GR evolution.',
        state='Original surface density, temperature, isotope mixture, gas entropy, ion and molecular inventories; all initially significant ion stages constrained, H2/H2+ fixed at initial negligible inventories; smaller stages and internal levels conditionally equilibrated.',
        domain='ln(rho/rho0) from0 to-6; native T>=100K. No extrapolated cold EOS or fitted relaxation time.',
        budget=dict(seconds=60,EOS_calls=1200,cpu_threads=1,new_whole_star_steps=0),
        forecast='Current probe60 calls1.43s.49 new density nodes with up to6 entropy iterations and16 population iterations have an explicit1200call cap. Expected5-30s; stop60s.',
        gates=dict(entropy_scaled=1e-8,population=1e-8,positive_compressibility=True,coarse_fine_speed=.002),
        stop='No density/domain extension, gate relaxation or inherited whole-flow/charge acceptance.',
        repair='Preserve entropy-root noise and the cold Newton-electron NaN at density-4.25. Retain the stricter1e-12 inventory solve; reuse34 accepted nodes. Predict native cumulative ion affinities with the exact native ionization energies before each colder call, so the internal electron solve starts near the constrained rather than neutral gas. Original equations and gates remain unchanged; debit both attempts.',
        prior_calls=previous['EOS_calls'],prior_seconds=previous['seconds'],source_sha256=sha(__file__)))
    ion=Ions(cap=1200-previous['EOS_calls']);fan=ion.fan;lr=np.log(fan.rho);lt=np.log(fan.T)
    base=ion.snapshot(lr,lt,np.zeros(318));target=base['number_fractions'];s0=base['eos'][3]
    values=[];affinities=[];counts=[];errors=[]
    density=np.linspace(0,-6,49);lastT=lt;lastx=0.;fields=np.zeros(318)
    solvedT=[];start=0
    if (OUT/'second-adiabat-native-states.npz').exists():
        cache=np.load(OUT/'second-adiabat-native-states.npz');start=failures[-1]['nodes']
        # The second attempt's accepted count is authoritative, independent
        # of directory enumeration order.
        start=json.loads((OUT/'second-adiabat-failure.json').read_text())['nodes']
        for x in density[:start]:
            j=np.flatnonzero(abs(cache['lrho']-(lr+x))<1e-12)[-1]
            row=cache['eos'][j];theta=float(cache['logT'][j]);T=np.exp(theta)
            values.append(row);affinities.append(cache['fields'][j]);counts.append(cache['number_fractions'][j])
            errors.append([(row[3]-s0)*T/(1.5*row[1]/row[0]),float(np.max(abs(cache['number_fractions'][j]-target))/target.sum())]);solvedT.append(T)
        fields=affinities[-1].copy();lastT=np.log(solvedT[-1]);lastx=density[start-1]
    raw=(SOURCE/'mod_ionization_data.f90').read_text()
    body=raw[raw.index('monatomic_ip(nions) ='):];body=body[:body.index(']')]
    body='\n'.join(line.split('!')[0] for line in body.splitlines())
    potentials=np.array([float(v.replace('d','e')) for v in re.findall(r'([-+]?[0-9]+\.[0-9]*(?:[edED][-+]?[0-9]+)?)_fp_kind',body)])
    assert len(potentials)==316,('Native ionization constants',len(potentials))
    c2=float(json.loads((native.OUT.parent/'def-photon-shared-atomic/catalog.json').read_text())['constants'][2])
    energy=np.concatenate([np.cumsum(potentials[ion.starts[e]:ion.starts[e+1]]) for e in range(24)])*c2
    charges=np.concatenate([np.arange(1,z+1) for z in ion.Z])
    active=np.concatenate([target[e,1:z+1]>target[e].sum()*1e-18 for e,z in enumerate(ion.Z)])
    diss=float(re.search(r'h2diss = ([0-9.]+)_fp_kind',raw).group(1))*c2
    molecule_ion=diss+potentials[0]*c2-21375.95*c2
    try:
        for i,x in enumerate(density):
            if i<start:continue
            theta=lastT+2*(x-lastx)/3
            fields[:316]+=np.where(active,energy*(np.exp(-theta)-np.exp(-lastT))+charges*((x-lastx)-1.5*(theta-lastT)),0.)
            fields[316]+=-diss*(np.exp(-theta)-np.exp(-lastT))-(x-lastx)+1.5*(theta-lastT)
            fields[317]+=molecule_ion*(np.exp(-theta)-np.exp(-lastT))+(x-lastx)-1.5*(theta-lastT)
            for iteration in range(8):
                assert np.exp(theta)>=100,'Registered cold domain boundary'
                a,n,error=ion.constrain(lr+x,theta,target,fields,target_molecules=base['molecular_H_fractions'],tolerance=1e-12);fields=ion.fields.copy()
                row=a['eos'];T=np.exp(theta);scale=1.5*row[1]/row[0]
                residual=(row[3]-s0)*T/scale
                if abs(residual)<2e-10:break
                assert abs(residual)<.15,'Entropy root left bounded local step'
                theta-=residual
            else:raise AssertionError('Fixed-ion entropy solve')
            values.append(row);affinities.append(fields.copy());counts.append(a['number_fractions']);errors.append([residual,error])
            solvedT.append(T)
            lastT=theta;lastx=x
            write(OUT/'adiabat-progress.json',dict(nodes=i+1,EOS_calls=ion.calls,seconds=time.monotonic()-begin,T=T))
        rows=np.array(values)
        # Retain the solved temperatures independently of any EOS output index.
        T=np.array(solvedT)
        curves=[]
        for stride in [2,1]:
            d=-density[::stride];a=rows[::stride]
            gamma=-PchipInterpolator(d,np.log(a[:,1])).derivative()(d)
            enthalpy=fan.cx*fan.c**2+a[:,2]+a[:,1]/a[:,0]
            cs=np.sqrt(gamma*a[:,1]/a[:,0]/enthalpy)
            rapidity=cumulative_simpson(cs,x=d,initial=0);v=np.tanh(rapidity);xi=(v-cs)/(1-v*cs)
            assert np.all(gamma>1) and np.all(np.diff(xi)>0)
            curves.append(dict(gamma=gamma,cs=cs,rapidity=rapidity,velocity=v,xi=xi))
        contrast=float(max(abs(curves[0]['rapidity']-curves[1]['rapidity'][::2]))/max(curves[1]['rapidity']))
        eq=dict(np.load(native.OUT/'fine.npz'));comparisons=[]
        for x in [-1.,-2.,-4.,-6.]:
            i=int(np.argmin(abs(density-x)));j=int(np.argmin(abs(eq['log_density_ratio']-x)))
            comparisons.append(dict(log_density_ratio=x,fixed_T=float(T[i]),LTE_T=float(eq['T'][j]),
                pressure_ratio=float(rows[i,1]/eq['raw'][j,1]),fixed_gamma=float(curves[1]['gamma'][i]),
                LTE_gamma=float(eq['raw'][j,4]),fixed_H_ion=float(rows[i,14]),LTE_H_ion=float(eq['raw'][j,14])))
        np.savez_compressed(OUT/'adiabat.npz',log_density_ratio=density,T=T,raw=rows,fields=affinities,
            number_fractions=counts,errors=errors,**curves[1])
        result=dict(classification='Counterexample candidate',passed=bool(contrast<.002 and np.max(abs(errors),axis=0)[0]<1e-8),
            nodes=len(rows),EOS_calls=ion.calls,seconds=time.monotonic()-begin,
            coarse_fine_rapidity_relative=contrast,maximum_entropy_scaled=float(np.max(abs(errors),axis=0)[0]),
            maximum_population_error=float(np.max(abs(errors),axis=0)[1]),comparisons=comparisons,
            physical_frozen_ion_limit_certified=False,full_fluid_evolution=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/'adiabat.json',result);assert result['passed'];print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'adiabat-failure.json',dict(error=repr(exc),EOS_calls=ion.calls,seconds=time.monotonic()-begin,nodes=len(values)))
        raise
    finally:
        ion.save('adiabat-native-states.npz')
        write(OUT/'density-fallbacks.json',dict(rows=ion.density_fallbacks,source_sha256=sha(__file__)))
    signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
