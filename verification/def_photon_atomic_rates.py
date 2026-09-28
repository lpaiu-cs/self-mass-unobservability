"""Actual atomic amplitudes with the common EOS states, no stellar integration.

Counterexample candidate: finite atom data, translated thresholds and a frozen
thermodynamic reverse-rate closure. This is not a microscopic plasma opacity.
"""
from pathlib import Path
from types import FunctionType
import argparse
import ctypes
import difflib
import json
import resource
import shutil
import signal
import subprocess
import time
import numpy as np
import sympy as sp
import def_photon_shared_atomic as a
import def_photon_profile_integral as profile

OUT=a.OUT.parent/'def-photon-atomic-rates'
CACHE=a.CACHE.parent/'photon-atomic-rates'
write=a.p.a.write
digest=a.p.a.digest
# SIGK's unsuffixed F77 literal is rounded to binary32 before promotion.
NATIVE_H=float(np.float32(6.6256e-27))


def continuum_rates(bound_density,parent_density,inverse_factor,sigma,u,T,h,c,k):
    """Pointwise volume coefficients; parent is the next ground-state density.

    The frozen inverse_factor comes from EOS chemical potentials and level
    weights. Recompute it for another bath; it is not a non-LTE kinetic model.
    """
    forward=bound_density*sigma
    reverse=parent_density*inverse_factor*sigma
    stimulated=reverse*np.exp(-u)
    frequency=u*k*T/h
    spontaneous=stimulated*(2*h*frequency**3/c**2)
    return forward,reverse,stimulated,spontaneous,forward-stimulated


def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='69e3dfa9',
        claim='Use parsed atomic oscillator strengths and native photoionization cross sections with the common finite EOS levels; derive reverse-rate factors from native chemical stationarity, not a fitted opacity or population ratio.',
        decision='Accept the reusable pointwise rate provider only if the same EOS state is recovered, chemical ratios agree, amplitudes are finite/nonnegative, subthreshold continuum is zero, and forward/reverse LTE balance passes. Do not run photon or GR evolution on this partial opacity.',
        model='Fixed stored state. Preserve imported excitation differences and translate photoionization energy by the EOS ground threshold shift at fixed outgoing electron energy. Use vacuum photon frequencies. Bound-bound joint survival remains min(w_i,w_j); continuum inverse is a declared thermodynamic closure with frozen native activities, not a microscopic Milne derivation. No guessed dissolved continuum or missing transitions.',
        budget=dict(cpu_threads=1,total_compute_seconds=240,build_seconds=120,build_attempts=2,native_EOS_calls=1,extraction_seconds=30,analysis_seconds=60,native_extractions=2,new_stellar_steps=0,coupled_runs=0,memory_GB=2),
        forecast='Prior EOS 12 calls took 0.410 s; this one call is capped at10 s. Existing native spectra took about6-13 s; parser plus cross-section sampling has not been measured and is capped at30 s per extraction. Two builds total capped120 s. No unmeasured cost is asserted as an ETA.',
        gates=dict(same_state_bitwise=True,chemical_log_ratio=1e-9,kirchhoff_relative=1e-9,threshold_negative=0),
        bindings={str(p):digest(p) for p in [a.OUT/'catalog.json',a.OUT/'eos-state.npz',a.OUT/'optical-levels.npz',a.LIB,a.CACHE/'gas.so',profile.BUILD/'synspec54.f']}))
    for p in profile.BUILD.glob('*.FOR'):shutil.copyfile(p,CACHE/p.name)
    patch_source()
    work=CACHE/'run';work.mkdir()
    previous=profile.CACHE/'balanced-fine'
    for name in ['fort.5','fort.55','tas']:shutil.copyfile(previous/name,work/name)
    for name in ['data','atom']:(work/name).symlink_to((previous/name).resolve(),target_is_directory=True)


def patch_source():
    source=(profile.BUILD/'synspec54.f').read_text();patched=source
    anchor='      CALL START\n'
    assert patched.count(anchor)==1
    patched=patched.replace(anchor,anchor+'''C     Standalone atomic provider: no atmosphere or spectral evolution.
      CALL ATOMIC_REQUESTS
      STOP
''',1)
    start=patched.index('      SUBROUTINE RDATA(');end=patched.index('      SUBROUTINE NSTPAR(',start)
    block=patched[start:end]
    anchor='''      READ(IUNIT,*,END=30,ERR=30) II,JJ,MODE,IFANCY,ICOLIS,
     *                              IFRQ0,IFRQ1,OSC,CPARAM'''
    assert block.count(anchor)==1
    block=block.replace(anchor,anchor+'''
C     Export even unsupported records to preserve the coverage boundary.
      WRITE(92,'(7i8,1p2e26.17)') ION,II,JJ,
     *   NFIRST(ION)+II-1,NFIRST(ION)+JJ-1,MODE,IFANCY,OSC,CPARAM''')
    anchor='      IBF(II)=IFANCY'
    assert block.count(anchor)==1
    block=block.replace(anchor,'''      WRITE(98,'(6i8,1pe26.17)') ION,II,JJ,
     *   NNEXT(ION),MODE,IFANCY,OSC
'''+anchor)
    patched=patched[:start]+block+patched[end:]
    # SIGK was designed for wavelength spectra and shifts edges to air. This
    # provider uses physical vacuum photon energy throughout its interface.
    start=patched.index('      FUNCTION SIGK(');end=patched.index('      END',start)
    block=patched[start:end];anchor='      IF(WL0.GT.vaclim) THEN'
    assert block.count(anchor)==1
    block=block.replace(anchor,'      IF(.FALSE.) THEN')
    patched=patched[:start]+block+patched[end:]
    patched+='''
      SUBROUTINE ATOMIC_REQUESTS
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      COMMON/TOPCS/CTOP(MFIT,MCROSS),XTOP(MFIT,MCROSS)
      INTEGER NREQ,IREQ,ILEV
      DO ILEV=1,NLEVEL
         WRITE(93,'(3i8,1p2e26.17)') ILEV,IBF(ILEV),
     *      NQUANT(ILEV),ENION(ILEV),G(ILEV)
      END DO
      READ(90,*) NREQ
      DO IREQ=1,NREQ
         READ(90,*) ILEV,FR
         IF(ILEV.LT.1.OR.ILEV.GT.NLEVEL.OR.FR.LE.0)
     *      STOP 'INVALID ATOMIC REQUEST'
         SIG=SIGK(FR,ILEV,0)
         WRITE(91,'(2i8,1p2e26.17)') IREQ,ILEV,FR,SIG
      END DO
      END
'''
    (CACHE/'synspec54.f').write_text(patched)
    (OUT/'source.patch').write_text(''.join(difflib.unified_diff(source.splitlines(True),patched.splitlines(True),fromfile='phase64/synspec54.f',tofile='phase67/synspec54.f')))


def build():
    assert not (OUT/'build.json').exists()
    cmd=['gfortran','-std=legacy','-fno-automatic','-mcmodel=medium','-g','-fcheck=all','-fbacktrace','synspec54.f','-o','atomic-rates']
    begin=time.monotonic();p=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=90)
    (OUT/'build.log').write_text(p.stdout+p.stderr)
    write(OUT/'build.json',dict(command=cmd,returncode=p.returncode,seconds=time.monotonic()-begin,
        source_sha256=digest(CACHE/'synspec54.f'),binary_sha256=digest(CACHE/'atomic-rates') if p.returncode==0 else None))
    assert p.returncode==0,p.stderr
    print('BUILT',time.monotonic()-begin,flush=True)


def station():
    assert not (OUT/'station.npz').exists();signal.alarm(10);begin=time.monotonic()
    prior=a.p.a.previous;gas_module=prior.old.previous.matter.old
    gas=object.__new__(gas_module.GasEOS);init=gas_module.GasEOS.__init__
    FunctionType(init.__code__,dict(init.__globals__,BRIDGE=a.CACHE/'gas.so'),closure=init.__closure__)(gas)
    gas.inventory_lib=gas.gas_lib;native=gas.gas_lib.ionization_inventory
    def call(mode,value,t,eps,out,info):
        raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info)
        out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
    gas.call=call
    d,_=prior.base.inputs();r=float(d['lnd'][0]);X=d['X'][0]
    T=float(np.load(a.p.OUT/'eos-state.npz')['T'][0])
    snap=prior.inventory_reader.InventoryEOS.snapshot(gas,r,float(np.log(T)),X)
    saved=np.load(a.OUT/'eos-state.npz')
    assert all(np.array_equal(saved[k][0],snap[k]) for k in snap)
    values={}
    for key,n in [('dv',316),('zero',24),('ce',316),('binding',316),('plop',318),('tc2',1),('flags',3),('active',24)]:
        dtype=ctypes.c_int if key in ['flags','active'] else ctypes.c_double
        values[key]=np.ctypeslib.as_array((dtype*n).in_dll(gas.gas_lib,'__mod_nuvar_MOD_station_'+key)).copy()
    values['T']=np.array(T);values['rho']=np.array(np.exp(r))
    np.savez_compressed(OUT/'station.npz',**values)
    write(OUT/'station.json',dict(classification='Counterexample candidate',same_state_bitwise=True,native_EOS_calls=1,seconds=time.monotonic()-begin))
    signal.alarm(0);print('STATION',time.monotonic()-begin,flush=True)


def extract():
    assert not (OUT/'cross-sections.npz').exists()
    catalog=json.loads((a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants']
    T=float(np.load(OUT/'station.npz')['T']);queries=[];metadata=[]
    # Threshold-centered probes plus a fixed excess-energy grid. These are
    # pointwise operator checks, never a resonance integration error estimate.
    excess=np.r_[-1e-6,0.,1e-6,np.geomspace(1e-5,60.,257)]
    for stage,row in enumerate(catalog['rows']):
        for local,index in enumerate(np.flatnonzero(row['included'])):
            new=row['binding_cm_inverse'][index]*erg
            old=new-row['anchor_shift_erg']
            for x in excess:
                energy=new+x*k*T;native_energy=old+x*k*T
                assert energy>0 and native_energy>0
                queries.append([row['global_levels'][index],native_energy/NATIVE_H])
                metadata.append([stage,local,energy/(k*T),x])
    work=CACHE/'run'
    with (work/'fort.90').open('w') as f:
        f.write(str(len(queries))+'\n');np.savetxt(f,queries,fmt=['%d','%.17e'])
    begin=time.monotonic()
    with (work/'fort.5').open() as stdin:
        p=subprocess.run([str(CACHE/'atomic-rates')],cwd=work,stdin=stdin,capture_output=True,text=True,timeout=30)
    (OUT/'native.log').write_text(p.stdout+p.stderr)
    write(OUT/'extraction.json',dict(returncode=p.returncode,seconds=time.monotonic()-begin,requests=len(queries),maxrss_kB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss))
    assert p.returncode==0,p.stderr
    raw=np.loadtxt(work/'fort.91');assert raw.shape==(len(queries),4)
    assert np.array_equal(raw[:,0],np.arange(1,len(queries)+1))
    assert np.array_equal(raw[:,1],np.array(queries)[:,0])
    assert np.allclose(raw[:,2],np.array(queries)[:,1],rtol=1e-15,atol=0)
    assert np.isfinite(raw).all() and (raw[:,3]>=0).all()
    np.savez_compressed(OUT/'cross-sections.npz',query=np.array(queries),metadata=np.array(metadata),sigma=raw[:,3])
    for name in ['fort.92','fort.93','fort.98']:shutil.copyfile(work/name,OUT/name)
    print('EXTRACTED',len(queries),'in',time.monotonic()-begin,flush=True)


def analyze():
    assert not (OUT/'result.json').exists();signal.alarm(60);begin=time.monotonic()
    catalog=json.loads((a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants'];h=erg/c
    st=np.load(OUT/'station.npz');T=float(st['T']);rho=float(st['rho'])
    saved=np.load(a.OUT/'eos-state.npz');optical=np.load(a.OUT/'optical-levels.npz')
    inv=saved['number_fractions'][0];nav=6.02214076e23
    # Compile-time parameter (not a dynamic-library symbol); bind its source.
    constants=a.p.SOURCE/'mod_free_eos_constants.f90'
    assert 'avogadro = 6.02214076e23_fp_kind' in constants.read_text()
    conversion=rho*nav
    coeff=[];ratio_errors=[];level_map={};rows=[];coefficient_rows=[];parent_fractions=[]
    for stage,row in enumerate(catalog['rows']):
        element=row['element_index']-1;q=row['charge'];i=row['ion_index']-1
        lower=-st['zero'][element] if q==0 else st['ce'][i-1]+st['dv'][i-1]-st['binding'][i-1]*st['tc2'][0]
        upper=st['ce'][i]+st['dv'][i]-st['binding'][i]*st['tc2'][0]
        ln_ratio=upper-lower
        ratio_errors.append(float(abs(np.log(inv[element,q+1]/inv[element,q])-ln_ratio)))
        # This uses native free-energy stationarity coefficients, independently
        # of the saved ion population ratio. Electron activity is in dv.
        parent_match=np.flatnonzero(np.all(saved['ids'][0]==[element+1,q+1],axis=1))
        assert len(parent_match)<=1
        parent_L=0. if not len(parent_match) else saved['value'][0,parent_match[0],3]/saved['scale'][0,parent_match[0]]
        parent_fraction=np.exp(-parent_L);parent_fractions.append(parent_fraction)
        # Native photoionization channels end in the next ground state.
        # Make its LTE fraction explicit instead of treating the whole stage
        # as a ground state. The factor cancels only under this declared LTE.
        K=np.exp(-np.longdouble(ln_ratio))*optical[f'fraction_{stage}']/parent_fraction
        coeff.append(K)
        for local,index in enumerate(np.flatnonzero(row['included'])):
            level=row['global_levels'][index];level_map[level]=(stage,local)
            coefficient_rows.append([level,stage,local,float(K[local]),row['binding_cm_inverse'][index]*erg,row['anchor_shift_erg'],float(parent_fraction)])
        rows.append(dict(Z=row['Z'],charge=q,log_parent_over_stage=float(ln_ratio),parent_ground_fraction=float(parent_fraction),chemical_log_ratio_error=ratio_errors[-1]))
    raw=np.loadtxt(OUT/'fort.92');lines=[];excluded=[];bb_errors=[]
    echarge=1.602176634e-19*.1*c
    me=ryd*(h/echarge**2)*(h/echarge)**2*c/(2*np.pi**2)
    sigma_integral=np.pi*echarge**2/(me*c)
    for index,line in enumerate(raw):
        lo,hi,mode=int(line[3]),int(line[4]),int(line[5]);f=float(line[7])
        reason=None
        if abs(mode)!=1:reason='disabled_or_unsupported_mode'
        elif lo not in level_map or hi not in level_map:reason='outside_common_retained_levels'
        elif f<=0:reason='nonpositive_oscillator_strength'
        if reason:excluded.append(dict(row=index,reason=reason));continue
        stage,l=level_map[lo];stage2,u=level_map[hi];assert stage==stage2
        weight=optical[f'weight_{stage}'];lw=optical[f'logw_{stage}'];pop=optical[f'populations_{stage}']*conversion
        energy=(optical[f'excitation_cm_inverse_{stage}'][u]-optical[f'excitation_cm_inverse_{stage}'][l])*erg
        if energy<=0:excluded.append(dict(row=index,reason='nonpositive_photon_energy'));continue
        nu=energy/h;ratio=np.exp(lw[u]-lw[l]);up=min(np.longdouble(1),ratio);down=min(np.longdouble(1),1/ratio)
        amplitude=sigma_integral*f
        abs0=pop[l]*amplitude*up
        stim=pop[u]*amplitude*weight[l]/weight[u]*down
        spontaneous=stim*(2*h*nu**3/c**2)
        planck=2*h*nu**3/c**2/np.expm1(energy/(k*T))
        net=abs0-stim
        err=abs(spontaneous-net*planck)/max(abs(spontaneous),abs(net*planck),1e-300)
        bb_errors.append(float(err));assert net>0
        lines.append([stage,l,u,float(nu),f,float(abs0),float(stim),float(spontaneous),float(net),float(up),float(down)])
    cross=np.load(OUT/'cross-sections.npz');meta=cross['metadata'];sigma=cross['sigma'];forward=np.zeros(len(sigma));reverse=np.zeros(len(sigma))
    continuum=np.loadtxt(OUT/'fort.98');channels={}
    for record in continuum:
        level=int(record[1])
        if level not in level_map:continue
        assert level not in channels,'Multiple photoionization channels need an explicit sum'
        assert record[2]==record[3],'An excited parent channel needs its own threshold'
        channels[level]=record
    assert set(channels)==set(level_map),'Never invent a hydrogenic channel for an unlisted level'
    for stage,row in enumerate(catalog['rows']):
        mask=meta[:,0]==stage;local=meta[mask,1].astype(int)
        pop=optical[f'populations_{stage}']*conversion
        parent=inv[row['element_index']-1,row['charge']+1]*conversion*parent_fractions[stage]
        forward[mask],reverse[mask],_,_,_=continuum_rates(pop[local],parent,coeff[stage][local],sigma[mask],meta[mask,2],T,h,c,k)
    u=meta[:,2];stim=reverse*np.exp(-u);net=forward-stim
    nu=u*k*T/h;spontaneous=stim*(2*h*nu**3/c**2)
    planck=2*h*nu**3/c**2/np.expm1(u)
    positive=forward>0
    bf_error=float(np.max(abs(spontaneous[positive]-net[positive]*planck[positive])/np.maximum(abs(spontaneous[positive]),abs(net[positive]*planck[positive]))))
    assert np.all(sigma[meta[:,3]<0]==0) and np.all(net>=0)
    assert all(np.isfinite(x).all() for x in [forward,reverse,stim,net,spontaneous])
    np.savez_compressed(OUT/'rates.npz',lines=np.array(lines),continuum_coefficients=np.array(coefficient_rows),bf_forward=forward,bf_reverse_prefactor=reverse,bf_stimulated=stim,bf_spontaneous=spontaneous,bf_net=net)
    write(OUT/'coverage.json',dict(classification='Counterexample candidate',atomic_stages=rows,parsed_lines=len(raw),retained_lines=len(lines),excluded=excluded,
        threshold_zero_samples=int(np.sum(meta[:,3]<0)),continuum_positive_samples=int(positive.sum()),sampled_levels=len(level_map),
        unresolved=['Missing H/He II provider amplitudes in this metal/He I extraction','Missing ion stages and source model transitions','Dissolved bound states and oscillator-strength redistribution','Line broadening, finite-profile redistribution and resonant continuum integration','Microscopic nonideal recombination activity correction','Opacity derivatives and off-equilibrium population kinetics']))
    passed=bool(max(ratio_errors)<1e-9 and max(bb_errors)<1e-9 and bf_error<1e-9)
    result=dict(classification='Counterexample candidate',passed=passed,decision='COMMON_FINITE_ATOMIC_POINTWISE_RATES_IMPLEMENTED_PHYSICAL_COMPLETENESS_OPEN',
        same_state_bitwise=True,native_EOS_calls=1,atomic_stages=len(rows),common_levels=len(level_map),bound_bound_lines=len(lines),bound_free_samples=len(sigma),
        chemical_log_ratio_max=max(ratio_errors),bound_bound_balance_relative=max(bb_errors),bound_free_balance_relative=bf_error,
        line_net_strength_sum_cm_inverse_Hz=float(np.array(lines)[:,8].sum()),seconds=time.monotonic()-begin,avogadro=nav,
        new_stellar_steps=0,coupled_runs=0,physical_opacity_certified=False,full_dynamic_charge_solved=False,
        bindings={str(p):digest(p) for p in [a.OUT/'eos-state.npz',a.OUT/'optical-levels.npz',a.OUT/'catalog.json',OUT/'station.npz',OUT/'cross-sections.npz',OUT/'fort.92',OUT/'fort.93',OUT/'fort.98',CACHE/'atomic-rates',a.LIB,a.CACHE/'gas.so',constants]})
    write(OUT/'result.json',result)
    n,K,s,z,b=sp.symbols('n K s z b',positive=True)
    A=n*K*s;R=n*K*s*z;J=R*b;B=b*z/(1-z)
    assert sp.simplify(J-(A-R)*B)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        identity='If the native thermodynamic coefficient K independently implies n_bound=n_parent_ground K, absorption A=n_bound sigma, stimulated R=n_parent_ground K sigma exp(-u), and spontaneous J=R 2h nu^3/c^2 obey J=(A-R) Bnu. Actual sigma is retained in both directions.',
        limitation='Conditional fixed-state equilibrium closure. It does not derive microscopic nonideal recombination or off-equilibrium activities.'))
    signal.alarm(0);print(json.dumps(result,indent=2),flush=True);assert passed


def verify():
    result=json.loads((OUT/'result.json').read_text());assert result['passed']
    for record in [json.loads((OUT/'plan.json').read_text()),result]:
        for name,sha in record['bindings'].items():assert digest(Path(name))==sha,name
    for name,sha in json.loads((OUT/'provider-inputs.json').read_text())['sha256'].items():assert digest(Path(name))==sha,name
    rates=np.load(OUT/'rates.npz');cross=np.load(OUT/'cross-sections.npz')
    assert np.all(cross['sigma']>=0) and np.all(cross['sigma'][cross['metadata'][:,3]<0]==0)
    assert np.all(rates['bf_net']>=0) and np.all(rates['lines'][:,8]>0)
    catalog=json.loads((a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants'];h=erg/c
    st=np.load(OUT/'station.npz');T=float(st['T'])
    native=np.loadtxt(OUT/'fort.93');levels={int(x[0]):x for x in native}
    for row in catalog['rows']:
        for index in np.flatnonzero(row['included']):
            data=levels[row['global_levels'][index]]
            assert data[4]==row['weight'][index]
            assert np.isclose(data[3],row['binding_cm_inverse'][index]*erg-row['anchor_shift_erg'],rtol=1e-12,atol=0)
    # Positive control: the implementation must retain independent absorption
    # and emission, rather than forcing j=chi B after changing populations.
    f,r,s,j,net=continuum_rates(2.,3.,2/3,1e-18,1.,T,h,c,k)
    f2,r2,s2,j2,net2=continuum_rates(2.2,3.,2/3,1e-18,1.,T,h,c,k)
    assert j==j2 and net2>net and f2>f
    assert abs(net*j2/(net2*j)-1)>.05
    assert not result['physical_opacity_certified'] and result['new_stellar_steps']==0
    if (OUT/'manifest.json').exists():
        for name,sha in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert digest(Path(name))==sha,name
    print('PASS atomic mapping, positive rates, threshold zeros and independent-emission positive control; physical completeness remains open',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','build','station','extract','analyze','verify'])
    globals()[p.parse_args().action]()
