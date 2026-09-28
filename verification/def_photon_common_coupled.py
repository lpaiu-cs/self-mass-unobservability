"""Apply the finite common atomic input to the existing coupled time equation.

Counterexample candidate. Only listed atomic channels are present. This is an
actual local photon/material/spatial evolution, not a full stellar GR solve.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import inspect
import json
import resource
import shutil
import signal
import time
import numpy as np
import def_photon_mhd_spectrum as a
import def_photon_atomic_rates as rates
import def_photon_profile_integral as profile
import def_photon_actual_opacity_validate as v

OUT=a.OUT/'coupled';RATES=OUT/'rates';write=a.write;digest=a.digest
candidate=SimpleNamespace(**dict(vars(a.old),OUT=a.OUT,CACHE=a.CACHE,LIB=a.LIB))


def prepare():
    assert not OUT.exists();OUT.mkdir();RATES.mkdir()
    paths=[Path(__file__),Path(a.__file__),Path(rates.__file__),Path(profile.__file__),
        Path(v.__file__),Path(v.c.__file__),Path(v.c.m.__file__),
        Path(v.c.m.previous.__file__),Path(v.old.__file__),
        a.OUT/'eos-state.npz',a.OUT/'optical-levels.npz',a.OUT/'result.json',
        a.OUT/'catalog.json',a.LIB,a.CACHE/'gas.so',rates.OUT/'cross-sections.npz',
        rates.OUT/'fort.92',rates.OUT/'fort.93',rates.OUT/'fort.98',v.old.OUT/'bank.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the newly implemented common finite EOS/rates to the existing photon/material/spatial time equation. Test the original numerical gates without changing the initial perturbation, duration, cells or time paths.',
        decision='Report separately actual coupled execution, numerical acceptance, same-state closure for retained channels, and whether the ORIGINAL complete-opacity failure is resolved. Missing H/He II channels prohibit the last claim even if this finite model passes.',
        model='Replace legacy absorption completely with the listed 566 levels and 3674 lines. No inconsistent old H/He background is carried over. Maxwell/RPA inventory, electron density and material heat capacity use the new EOS. Thermal Doppler plus natural widths from retained Einstein coefficients only; pressure widths and missing transitions are not certified. Frequency-dependent thermal reverse profile implements the declared fixed-bath LTE closure.',
        controls='Initial material perturbation1 K, zero photons, duration H/c, kH1, existing24134 cells and8 angular modes; fine16/32/64, half source knots64, Gauss4 profiles64 and zero-atomic-input64. No grid or duration expansion.',
        budget=dict(additional_native_EOS_calls=1,native_extractions=0,input_seconds=180,
            coupled_seconds=480,coupled_paths=6,cpu_threads=1,memory_GB=4,new_stellar_steps=0),
        budget_reason='The initial12 EOS attempts are exhausted and preserved. One additional identical-state call captures chemical coefficients and mass normalization; it adds no state or derivative sweep. Prior10 calls took0.365s. Reuse all cross sections. Prior coupled runs took about156s; new RPA table/opacity cost is unmeasured, forecast180-400s, hard480s.',
        gates=dict(time_order=1.8,time_error=.001,source_grid_plus_two_tails=.001,
            line_quadrature=.0001,energy_balance=1e-9,solver_residual=1e-11,
            entropy_growth=1e-10,energy_number_null=1e-10,chemical_ratio=1e-9),
        bindings={str(p):digest(p) for p in paths}))
    for name in ['cross-sections.npz','fort.92','fort.93','fort.98']:
        shutil.copyfile(rates.OUT/name,RATES/name)
    print('PREPARED actual coupled experiment',flush=True)


def inputs():
    assert not (OUT/'bank.npz').exists();signal.alarm(180);start=time.monotonic()
    for p,sha in json.loads((OUT/'plan.json').read_text())['bindings'].items():
        assert digest(Path(p))==sha,p
    # Preserve Phase67 source and its raw verdict; add the native mass scale
    # where normalized EOS fractions become volume densities.
    src=inspect.getsource(rates.station)
    anchor="    values['T']=np.array(T);values['rho']=np.array(np.exp(r))"
    assert src.count(anchor)==1
    src=src.replace(anchor,anchor+"\n    values['mass_scale']=np.array(float(((X/prior.inventory_reader.g.c.A)@gas.mapping)@gas.weights))\n    values['weights']=gas.weights.copy()")
    ns=dict(vars(rates),OUT=RATES,a=candidate)
    exec(compile(src,str(OUT/'station-source.py'),'exec'),ns)
    if not (RATES/'station.npz').exists():
        (OUT/'station-source.py').write_text(src);ns['station']()
    st=np.load(RATES/'station.npz');cx=float(st['mass_scale'])
    src=inspect.getsource(rates.analyze)
    assert src.count('conversion=rho*nav')==1
    src=src.replace('conversion=rho*nav',"conversion=rho*nav*float(st['mass_scale'])")
    (OUT/'rates-source.py').write_text(src)
    exec(compile(src,str(OUT/'rates-source.py'),'exec'),ns)
    if not (RATES/'result.json').exists():ns['analyze']()
    else:
        for p,sha in json.loads((RATES/'result.json').read_text())['bindings'].items():assert digest(Path(p))==sha,p
    signal.alarm(max(1,int(180-(time.monotonic()-start))))
    s=np.load(a.OUT/'eos-state.npz');b=dict(np.load(v.old.OUT/'bank.npz'))
    T=float(st['T']);rho=float(st['rho']);eos=s['eos'][0];base=v.c.m.base
    reference,_=v.c.m.base.inputs()
    assert T==float(b['T']) and rho==float(np.exp(reference['lnd'][0]))
    native_density_change=rho*cx/float(b['rho'])-1
    b['rho']=rho*cx;b['material_rho']=rho
    b['Cm']=rho*eos[10]/T;oldne=float(b['ne']);ne=float(eos[13])*base.Avogadro*1e6
    for key in ['rate_e','rate_C']:b[key]=b[key]*ne/oldne
    b['ne']=ne
    ni=s['number_fractions'][0]*rho*cx*base.Avogadro*1e6
    strength=(ni*np.arange(29)**2).sum(axis=1)/ne
    mass=st['weights']*v.c.m.atomic_mass;active=strength>0
    strength=strength[active];mass=mass[active]
    # The native molecular fraction counts H nuclei, two per H2+ ion.
    epsH=s['number_fractions'][0,0].sum()/(1-s['molecular_H_fractions'][0].sum())
    molecular=rho*cx*base.Avogadro*1e6*epsH*s['molecular_H_fractions'][0,1]/2
    if molecular>0:strength=np.r_[strength,molecular/ne];mass=np.r_[mass,2*st['weights'][0]*v.c.m.atomic_mass]
    charge_error=abs(((ni*np.arange(29)).sum()+molecular)/ne-1)
    assert charge_error<1e-10
    ld=np.sqrt(base.k*T/(4*np.pi*base.E2*ne));sigma=np.sqrt(base.k*T/base.m_e)
    np.savez_compressed(OUT/'inventory.npz',ion_strength=strength,ion_mass=mass,ne=ne,
        plasma_u=base.hbar*sigma/ld/(base.k*T),wave_lambda=base.k*T/(base.hbar*base.c)*ld)
    catalog=json.loads((a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants'];h=erg/c
    data=np.load(RATES/'rates.npz');lines=data['lines'];opt=np.load(a.OUT/'optical-levels.npz')
    echarge=1.602176634e-19*.1*c;me=ryd*(h/echarge**2)*(h/echarge)**2*c/(2*np.pi**2)
    lifetimes={}
    for row in lines:
        stage,lo,hi=map(int,row[:3]);g=opt[f'weight_{stage}']
        spontaneous=8*np.pi**2*echarge**2*row[3]**2/(me*c**3)*float(g[lo]/g[hi])*row[4]*row[10]
        lifetimes[stage,hi]=lifetimes.get((stage,hi),0.)+spontaneous
    profiles=[]
    for index,row in enumerate(lines):
        stage,lo,hi=map(int,row[:3]);mass_g=st['weights'][catalog['rows'][stage]['element_index']-1]/base.Avogadro
        width=row[3]*np.sqrt(2*k*T/mass_g)/c
        damping=(lifetimes.get((stage,lo),0.)+lifetimes[stage,hi])/(4*np.pi*width)
        profiles.append([index,row[3],row[4],0,row[5]/(np.sqrt(np.pi)*width),1/width,damping,1])
    profiles=np.array(profiles);np.save(OUT/'profiles.npy',profiles)
    # Reuse the same quadrature, with the common EOS constants (u=h nu/kT).
    src=inspect.getsource(profile.integrate)
    src=src.replace('stimulated_scale=4.79928144e-11/(T*scale)','stimulated_scale=1.')
    (OUT/'profile-source.py').write_text(src);pn=dict(vars(profile));exec(compile(src,str(OUT/'profile-source.py'),'exec'),pn)
    q4,i4=pn['integrate'](b,profiles,4);q8,i8=pn['integrate'](b,profiles,8)
    cross=np.load(RATES/'cross-sections.npz');meta=cross['metadata'];n=len(q8)
    edges=b['edges_u'][:n+1];g,w=np.polynomial.legendre.leggauss(8)
    capacity0=15*float(b['arad'])*T**3/np.pi**4
    def continuum(stride):
        total=np.zeros(n)
        for level in np.unique(cross['query'][:,0]):
            ids=np.flatnonzero(cross['query'][:,0]==level)
            # Below-edge probes are excluded, never interpolated across an edge.
            ids=ids[meta[ids,3]>=0];ids=ids[np.unique(np.r_[np.arange(0,len(ids),stride),len(ids)-1])]
            u=meta[ids,2];amp=data['bf_forward'][ids]
            if u[0]>=edges[-1]:continue
            upper=min(edges[-1],u[-1]);mesh=np.unique(np.r_[u[(u>=edges[0])&(u<=upper)],edges[(edges>=u[0])&(edges<=upper)],upper])
            mid=(mesh[:-1]+mesh[1:])/2;half=np.diff(mesh)/2;uq=mid[:,None]+half[:,None]*g
            density=capacity0*uq**4*np.exp(-uq)/(-np.expm1(-uq))**2
            # A and reverse activities have already been independently checked.
            absorb=np.interp(uq,u,amp)*(-np.expm1(-uq))
            amount=(density*absorb*half[:,None]*w).sum(axis=1)
            bins=np.searchsorted(edges,mid,side='right')-1
            np.add.at(total,bins,amount)
        return total
    bf=continuum(1);bc=continuum(2)
    for name,amount in [('fine',bf+q8),('coarse',bc+q8),('quadrature4',bf+q4)]:
        b['rate_a_'+name]=np.zeros_like(b['rate_a']);b['rate_a_'+name][:n]=c*amount/b['Ci'][:n]
    b['rate_a']=b['rate_a_fine'];np.savez_compressed(OUT/'bank.npz',**b)
    write(OUT/'input.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
        mass_scale=cx,native_density_change_from_legacy_request=native_density_change,
        charge_relative=charge_error,Cm=float(b['Cm']),ne_m3=ne,
        line_integration=[i4,i8],retained_channels_common_EOS=True,legacy_background_retained=False,
        physical_opacity_certified=False,original_complete_input_failure_resolved=False,
        limitations=json.loads((RATES/'coverage.json').read_text())['unresolved'],
        bindings={str(p):digest(p) for p in [OUT/'bank.npz',OUT/'inventory.npz',RATES/'rates.npz']}))
    signal.alarm(0);print('INPUT READY',time.monotonic()-start,flush=True)


def run():
    assert not (OUT/'response.json').exists();signal.alarm(480);start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    info=json.loads((OUT/'input.json').read_text())
    for p,sha in {**json.loads((OUT/'plan.json').read_text())['bindings'],**info['bindings']}.items():assert digest(Path(p))==sha,p
    # Rebind only the data directory; reuse the unchanged scattering and SDIRK solver.
    m=SimpleNamespace(**dict(vars(v.c.m),OUT=OUT))
    spec=FunctionType(v.c.spectrum.__code__,dict(vars(v.c),OUT=OUT,m=m))
    coef=FunctionType(v.c.raw_coefficients.__code__,dict(v.c.raw_coefficients.__globals__,OUT=OUT,spectrum=spec),argdefs=v.c.raw_coefficients.__defaults__)
    class Operator(v.c.Operator):
        __init__=FunctionType(v.c.Operator.__init__.__code__,dict(v.c.Operator.__init__.__globals__,coefficients=coef),argdefs=v.c.Operator.__init__.__defaults__,closure=v.c.Operator.__init__.__closure__)
    b=dict(np.load(OUT/'bank.npz'));b.update(velocity_order=256,q_points=257)
    op=Operator(b,8,1,24);duration=float(b['H']/b['c']);Cm=float(b['Cm']);rows=[]
    def evolve(name,steps):
        result=op.evolve(duration,steps)
        row=dict(name=name,algebra=v.algebra(op,b),**v.diagnostics(op,result));rows.append(row)
        write(OUT/(name+'.json'),row)
        np.savez_compressed(OUT/(name+'.npz'),T=result['T'],E=result['E'],u=op.u,Ci=op.Ci,Cm=Cm,duration=duration)
        print('COUPLED',name,'seconds',time.monotonic()-start,'balance',row['balance'],flush=True)
        return result
    runs=[evolve('fine-'+str(n),n) for n in [16,32,64]];fine=runs[-1]
    errors=[v.old.compare(x,y,op) for x,y in zip(runs[:-1],runs[1:])];order=float(np.log2(errors[0]/errors[1]))
    differences={}
    for name in ['coarse','quadrature4','zero-atomic']:
        b['rate_a']=b['rate_a_'+name] if name!='zero-atomic' else np.zeros_like(b['rate_a'])
        v.set_absorption(op,b['rate_a']);other=evolve(name+'-64',64)
        differences[name]=v.old.compare(fine,other,op)
    tail=info['line_integration'][-1]['omitted_wing_response_bound']
    passed=bool(order>=1.8 and errors[-1]<.001 and differences['coarse']+2*tail<.001
        and differences['quadrature4']<.0001 and all(r['balance']<1e-9 and r['energy_residual']<1e-9
        and r['solver_residual']<1e-11 and r['entropy_growth']<1e-10
        and max(r['algebra']['energy_number_null_relative'])<1e-10 for r in rows))
    write(OUT/'response.json',dict(classification='Counterexample candidate',actual_coupled_evolution=True,
        numerical_response_passed=passed,time_order=order,time_differences_initial=errors,
        differences_initial=differences,source_grid_plus_two_tails=differences['coarse']+2*tail,
        temperature_perturbation_K=float((fine['T']/np.sqrt(Cm)).real),duration_seconds=duration,
        rows=rows,coefficients=op.info,seconds=time.monotonic()-start,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        same_EOS_retained_channels=True,physical_opacity_certified=False,
        original_complete_input_failure_resolved=False,full_GR_photon_feedback_evolved=False,
        full_dynamic_charge_solved=False,
        decision='The new finite common input was actually evolved. This establishes or rejects its local numerical behavior; absent H/HeII and other channels prevent declaring the original complete-input physical failure fixed. No stellar GR trajectory is relabeled.'))
    signal.alarm(0);print('COUPLED FINISHED',passed,order,differences,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','inputs','run']);globals()[p.parse_args().action]()
