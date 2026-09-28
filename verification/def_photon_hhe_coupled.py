"""Put same-EOS H/HeII bound and charged-ion free-free channels into evolution.

Counterexample candidate: native finite hydrogenic/thermal continuum model.
Pressure profiles, dissolved strength and missing atomic channels remain open.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import ctypes
import difflib
import inspect
import json
import shutil
import signal
import subprocess
import time
import numpy as np
import def_photon_common_coupled as old

OUT=old.a.OUT.parent/'def-photon-hhe-coupled'
CACHE=old.a.CACHE.parent/'photon-hhe-coupled'
LIB=CACHE/'libfree_eos_hhe_metadata.so'
write=old.write;digest=old.digest


def prepare():
    assert not (OUT/'station.npz').exists() and not (OUT/'build.json').exists()
    OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
    if (OUT/'plan.json').exists():(OUT/'plan.json').rename(OUT/'plan-first-setup.json')
    paths=[Path(__file__),Path(old.__file__),old.a.LIB,old.a.CACHE/'gas.so',
        old.a.OUT/'eos-state.npz',old.a.OUT/'catalog.json',old.OUT/'bank.npz',
        old.OUT/'inventory.npz',old.OUT/'spectrum-256-257.npz',old.OUT/'input.json',
        old.RATES/'station.npz',old.rates.CACHE/'synspec54.f',old.rates.CACHE/'PARAMS.FOR']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='13e7c84a',
        claim='Add the missing H/HeII amplitudes and all charged-ion free-free absorption from the same saved EOS to the actual coupled solver. Recover the skipped native H ground metadata without changing EOS equations.',
        decision='Accept this channel connection only after bitwise same-EOS replay, independent native partition/chemical checks, actual same-interval coupled evolution and unchanged numerical gates. Never call remaining missing atomic/pressure/dissolved channels or full GR complete.',
        model='Native n<=10 MHD hydrogenic EOS levels; native STARK0 oscillator strengths, GAUNT hydrogenic photoionization and GFREE thermal charged-ion absorption. Translate each bound-free threshold at fixed electron energy. No merged-level opacity or modified free-free tail outside the EOS. Thermal Doppler plus retained natural widths, same declared LTE reciprocal closure as Phase68.',
        budget=dict(build_attempts=2,build_seconds=120,new_EOS_calls=1,input_seconds=120,
            coupled_seconds=480,coupled_paths=6,cpu_threads=1,memory_GB=4,new_stellar_steps=0),
        forecast='Prior two native builds totaled8.4s; one same-state EOS call0.110s; input26.84s and six coupled paths196.20s. New channel quadrature is unmeasured and capped120s. Reuse the exact inventory/RPA table and existing metal rates. Coupled forecast150-350s, hard480s; no automatic density/grid/duration expansion.',
        setup_revision='Initial text anchor matched more than the H branch and preparation stopped. An erroneous dependent build then reached absent atomic files. No EOS calls or evolution occurred. Restrict the patch to the H molecular branch and permit one completed build under the same120s aggregate cap; preserve the setup failure.',
        gates=dict(same_EOS_bitwise=True,partition_relative=1e-12,chemical_log_ratio=1e-9,
            balance_relative=1e-9,time_order=1.8,time_error=.001,source_plus_two_tails=.001,
            line_quadrature=.0001,energy_balance=1e-9,solver=1e-11),
        bindings={str(p):digest(p) for p in paths}))
    for p in old.a.CACHE.glob('*.f90'):shutil.copyfile(p,CACHE/p.name)
    for p in old.a.CACHE.glob('*.mod'):shutil.copyfile(p,CACHE/p.name)
    p=CACHE/'excitation_sum.f90';before=p.read_text()
    start=before.index('     ! now do neutral monatomic H')
    point=before.index('        exparg = exp(exparg)\n',start)+len('        exparg = exp(exparg)\n')
    assert '        ion = 1\n' in before[start:point]
    after=before[:point]+'''        shared_ground_logw(ion)=plop(ion)
        if(ifmhd_logical) shared_ground_logw(ion)=plop(ion)+ln_ground_occ
'''+before[point:]
    p.write_text(after)
    (OUT/'metadata.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='phase68/excitation_sum.f90',tofile='phase69/excitation_sum.f90')))
    # Compile the existing three atomic functions unchanged, with a thin array
    # bridge. No START/atmosphere/population override or fitted normalization.
    source=(old.rates.CACHE/'synspec54.f').read_text();parts=[]
    for start,end in [('      FUNCTION GAUNT(', '      FUNCTION GNTK('),
        ('      FUNCTION GFREE(', '      FUNCTION SFFHMI_old('),
        ('      SUBROUTINE STARK0(', '      FUNCTION STARKA(')]:
        i=source.index(start);parts.append(source[i:source.index(end,i)])
    (CACHE/'atomic.f').write_text(''.join(parts));shutil.copyfile(old.rates.CACHE/'PARAMS.FOR',CACHE/'PARAMS.FOR')
    (CACHE/'bridge.f90').write_text('''subroutine hydrogen_cross(n,z,level,t,fr,bf,ff) bind(C)
 use iso_c_binding
 implicit none
 integer(c_int),value::n,z,level
 real(c_double),value::t
 real(c_double),intent(in)::fr(n)
 real(c_double),intent(out)::bf(n),ff(n)
 real(c_double),external::gaunt,gfree
 integer::i
 do i=1,n
  bf(i)=2.815d29*real(z,c_double)**4/real(level,c_double)**5/fr(i)**3*gaunt(level,fr(i)/real(z*z,c_double))
  ff(i)=gfree(t,fr(i)/real(z*z,c_double))
 enddo
end subroutine
subroutine hydrogen_line(i,j,z,f) bind(C)
 use iso_c_binding
 implicit none
 integer(c_int),value::i,j,z
 real(c_double),intent(out)::f
 real(c_double)::x,w,f1,f0
 call stark0(i,j,z,x,w,f1,f0)
 f=f1+f0
end subroutine
''')
    for name in ['atomic.f','bridge.f90','PARAMS.FOR']:shutil.copyfile(CACHE/name,OUT/name)
    print('PREPARED H/He actual coupling',flush=True)


def build():
    assert not (OUT/'build.json').exists();begin=time.monotonic();receipts=[]
    def call(cmd,name):
        p=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=60)
        (OUT/(name+'.log')).write_text(p.stdout+p.stderr);receipts.append(dict(command=cmd,returncode=p.returncode))
        assert p.returncode==0,p.stderr
    call(['gfortran','-cpp','-O2','-fPIC','-fcheck=all','-ffree-line-length-none','-I'+str(CACHE),'-I'+str(old.a.old.PARENT),'-I'+str(old.a.old.p.SOURCE),'-c','mod_excitation.f90','-o','mod_excitation.o'],'excitation')
    objects=json.loads((old.a.OUT/'build.json').read_text())['reused_objects']
    call(['gfortran','-shared','-Wl,-soname,'+LIB.name,'mod_excitation.o',str(old.a.CACHE/'mod_free_eos_detailed.o'),str(old.a.CACHE/'mod_free_eos.o'),*objects,'-llapack','-lblas','-o',str(LIB)],'link')
    source=old.a.old.p.a.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    call(['gfortran','-O2','-fPIC','-shared','-I'+str(CACHE),'-I'+str(old.a.old.PARENT),str(source),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-lfree_eos_hhe_metadata','-o','gas.so'],'gas')
    call(['gfortran','-std=legacy','-O2','-fPIC','-shared','-fcheck=all','-ffree-line-length-none','atomic.f','bridge.f90','-o','atomic.so'],'atomic')
    write(OUT/'build.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-begin,
        receipts=receipts,runtime={str(p):digest(p) for p in [LIB,CACHE/'gas.so',CACHE/'atomic.so']},
        source={str(p):digest(p) for p in CACHE.glob('*.f*')}))
    assert time.monotonic()-begin<120;print('BUILT',time.monotonic()-begin,flush=True)


def inputs():
    assert not (OUT/'bank.npz').exists();signal.alarm(120);begin=time.monotonic()
    for p,sha in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert digest(Path(p))==sha,p
    # Exactly one same-state EOS call; expose H metadata on its molecular path.
    candidate=SimpleNamespace(**dict(vars(old.candidate),CACHE=CACHE,LIB=LIB))
    src=inspect.getsource(old.rates.station)
    anchor="    values['T']=np.array(T);values['rho']=np.array(np.exp(r))"
    src=src.replace(anchor,anchor+"\n    values['ground_logw']=np.ctypeslib.as_array((ctypes.c_double*318).in_dll(gas.gas_lib,'__mod_excitation_MOD_shared_ground_logw')).copy()")
    ns=dict(vars(old.rates),OUT=OUT,a=candidate)
    (OUT/'station-source.py').write_text(src);exec(compile(src,str(OUT/'station-source.py'),'exec'),ns);ns['station']()
    signal.alarm(max(1,int(120-(time.monotonic()-begin))))
    saved=np.load(old.a.OUT/'eos-state.npz');st=np.load(OUT/'station.npz');prior=np.load(old.RATES/'station.npz')
    assert all(np.array_equal(st[k],prior[k]) for k in st.files if k!='ground_logw')
    T=float(st['T']);rho=float(st['rho']);cx=float(prior['mass_scale']);NA=6.02214076e23
    catalog=json.loads((old.a.OUT/'catalog.json').read_text());erg,ryd,c2,c,k=catalog['constants'];h=erg/c
    _,_,electron_mass,h_mass=np.loadtxt(old.a.old.p.OUT/'levels-repaired/constants.txt')
    lib=ctypes.CDLL(str(LIB));array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
    fn=lib.hydrogen_terms;fn.argtypes=[ctypes.c_int,ctypes.c_double,array,array];fn.restype=None
    atom=ctypes.CDLL(str(CACHE/'atomic.so'));cross=atom.hydrogen_cross
    cross.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,array,array,array];cross.restype=None
    line=atom.hydrogen_line;line.argtypes=[ctypes.c_int]*3+[ctypes.POINTER(ctypes.c_double)];line.restype=None
    def atomic(z,n,frequency):
        frequency=np.ascontiguousarray(frequency.ravel());bf=np.zeros_like(frequency);ff=np.zeros_like(frequency)
        cross(len(frequency),z,n,T,frequency,bf,ff);assert np.all(bf>=0) and np.all(ff>0)
        return bf,ff
    levels=[];line_rows=[];checks=[]
    echarge=1.602176634e-19*.1*c;me=ryd*(h/echarge**2)*(h/echarge)**2*c/(2*np.pi**2)
    integral=np.pi*echarge**2/(me*c)
    native=np.loadtxt(old.rates.OUT/'fort.93')
    for z,element,charge,ion,global_ground in [(1,0,0,0,1),(2,1,1,2,35)]:
        out=np.zeros(30);fn(z,float(np.log(T)),np.ascontiguousarray(saved['x'][0]),out)
        tail=out.reshape(3,10)[0];terms=tail-np.r_[tail[1:],0.]
        bind=ryd*z*z/(1+electron_mass/(h_mass*(1 if z==1 else 4)))/np.arange(1,11)**2
        bind[0]=st['binding'][ion]-(st['binding'][ion-1] if charge else 0)
        g=2*np.arange(1,11,dtype=float)**2;logw=np.log(terms/g)-c2*bind/T
        logw[0]=st['ground_logw'][ion]
        ratios=terms[1:]*np.exp(-c2*bind[0]/T-logw[0])/g[0]
        fraction=np.r_[1.,ratios];fraction/=fraction.sum()
        j=np.flatnonzero(np.all(saved['ids'][0]==[element+1,charge],axis=1));assert len(j)==1
        L=float(saved['value'][0,j[0],3]/saved['scale'][0,j[0]])
        partition_error=abs(np.log1p(ratios.sum())/L-1);assert partition_error<1e-12
        inv=saved['number_fractions'][0];pop=inv[element,charge]*rho*cx*NA*fraction
        lower=-st['zero'][element] if charge==0 else st['ce'][ion-1]+st['dv'][ion-1]-st['binding'][ion-1]*st['tc2'][0]
        upper=st['ce'][ion]+st['dv'][ion]-st['binding'][ion]*st['tc2'][0]
        K=np.exp(lower-upper)*fraction;parent=inv[element,charge+1]*rho*cx*NA
        chemical_error=float(np.max(abs(parent*K/pop-1)));assert chemical_error<1e-9
        oldground=native[global_ground-1,3];rows=[]
        for lo in range(10):
            for hi in range(lo+1,10):
                f=ctypes.c_double();line(lo+1,hi+1,z,ctypes.byref(f))
                nu=(bind[lo]-bind[hi])*erg/h;ratio=np.exp(logw[hi]-logw[lo]);up=min(1.,ratio);down=min(1.,1/ratio)
                A=pop[lo]*integral*f.value*up;R=pop[hi]*integral*f.value*g[lo]/g[hi]*down
                err=abs(R*np.expm1(h*nu/(k*T))/(A-R)-1);assert err<1e-9
                rows.append([z,lo,hi,nu,f.value,A,R,up,down]);line_rows.append(rows[-1])
        levels.append(dict(z=z,element=element,bind=bind,logw=logw,pop=pop,K=K,parent=parent,oldground=oldground,weight=g,rows=rows))
        checks.append(dict(Z=z,partition_relative=float(partition_error),chemical_relative=chemical_error,
            native_ground_logw=float(logw[0]),previous_ground_logw=float(saved['ground_logw'][0,ion]),
            stage_population_cm3=float(pop.sum()),positive_lines=len(rows)))
    # Thermal/natural line profiles; population and energy enter explicitly.
    profiles=[]
    for level in levels:
        decay=np.zeros(10);g=level['weight']
        for z,lo,hi,nu,f,A,R,up,down in level['rows']:
            decay[hi]+=8*np.pi**2*echarge**2*nu**2/(me*c**3)*g[lo]/g[hi]*f*down
        for z,lo,hi,nu,f,A,R,up,down in level['rows']:
            mass=prior['weights'][level['element']]/NA;width=nu*np.sqrt(2*k*T/mass)/c
            profiles.append([len(profiles),nu,f,0,A/(np.sqrt(np.pi)*width),1/width,(decay[lo]+decay[hi])/(4*np.pi*width),1])
    b=dict(np.load(old.OUT/'bank.npz'));n=int((b['u']<=60).sum());edges=b['edges_u'][:n+1]
    src=(old.OUT/'profile-source.py').read_text();pn=dict(vars(old.profile));exec(compile(src,'phase68/profile-source.py','exec'),pn)
    l4,li4=pn['integrate'](b,np.array(profiles),4);l8,li8=pn['integrate'](b,np.array(profiles),8)
    cap=15*float(b['arad'])*T**3/np.pi**4;ne=float(b['ne'])/1e6
    ni=saved['number_fractions'][0]*rho*cx*NA;charge_density=ni.sum(axis=0)
    epsH=ni[0].sum()/(1-saved['molecular_H_fractions'][0].sum())
    charge_density[1]+=epsH*saved['molecular_H_fractions'][0,1]/2
    assert abs(charge_density@np.arange(29)/ne-1)<1e-10
    def continuum(order):
        gx,gw=np.polynomial.legendre.leggauss(order);bf_total=np.zeros(n);ff_total=np.zeros(n)
        for level in levels:
            for local in range(10):
                threshold=level['bind'][local]*erg/(k*T)
                mesh=np.unique(np.r_[edges[edges>threshold],max(edges[0],threshold)])
                if len(mesh)<2:continue
                mid=(mesh[:-1]+mesh[1:])/2;half=np.diff(mesh)/2;u=mid[:,None]+half[:,None]*gx
                native_energy=level['oldground']/(local+1)**2+(u-threshold)*k*T
                sigma,_=atomic(level['z'],local+1,native_energy/old.rates.NATIVE_H);sigma=sigma.reshape(u.shape)
                A,R0,R,J,net=old.rates.continuum_rates(level['pop'][local],level['parent'],level['K'][local],sigma,u,T,h,c,k)
                B=2*h*(u*k*T/h)**3/c**2/np.expm1(u)
                assert np.max(abs(J/(net*B)-1))<1e-9
                density=cap*u**4*np.exp(-u)/(-np.expm1(-u))**2
                amount=(density*net*half[:,None]*gw).sum(1)
                ids=np.searchsorted(edges,mid,side='right')-1;np.add.at(bf_total,ids,amount)
        mid=(edges[:-1]+edges[1:])/2;half=np.diff(edges)/2;u=mid[:,None]+half[:,None]*gx;nu=u*k*T/h
        amp=np.zeros_like(u)
        for z in np.flatnonzero(charge_density[1:]>0)+1:
            _,gaunt=atomic(int(z),1,nu)
            amp+=float(np.float32(3.694e8))/np.sqrt(T)*ne*charge_density[z]*z*z/nu**3*gaunt.reshape(u.shape)
        net=amp*(-np.expm1(-u));density=cap*u**4*np.exp(-u)/(-np.expm1(-u))**2
        ff_total=(density*net*half[:,None]*gw).sum(1)
        return bf_total,ff_total
    bf4,ff4=continuum(4);bf8,ff8=continuum(8)
    for name,added in [('fine',l8+bf8+ff8),('coarse',l8+bf4+ff4),('quadrature4',l4+bf8+ff8)]:
        b['rate_a_'+name][:n]+=c*added/b['Ci'][:n]
    b['rate_a']=b['rate_a_fine'];np.savez_compressed(OUT/'bank.npz',**b)
    np.savez_compressed(OUT/'channels.npz',lines=np.array(line_rows),profiles=np.array(profiles),bf4=bf4,bf8=bf8,ff4=ff4,ff8=ff8,line4=l4,line8=l8,
        populations=np.array([x['pop'] for x in levels]),binding_cm_inverse=np.array([x['bind'] for x in levels]),logw=np.array([x['logw'] for x in levels]))
    for name in ['inventory.npz','spectrum-256-257.npz','spectrum-256-257.json']:shutil.copyfile(old.OUT/name,OUT/name)
    old_tail=json.loads((old.OUT/'input.json').read_text())['line_integration'][-1]['omitted_wing_response_bound']
    li8['omitted_wing_response_bound']+=old_tail
    write(OUT/'input.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-begin,
        same_EOS_bitwise=True,new_EOS_calls=1,hydrogen_helium_checks=checks,line_integration=[li4,li8],
        lines_added=90,bound_free_levels_added=20,charged_ion_free_free=True,
        physical_opacity_certified=False,original_complete_input_failure_resolved=False,
        remaining='Pressure widths, dissolved/omitted oscillator strength, missing metal transitions, neutral/molecular absorption and bound-electron redistribution; atmosphere, moving GR material and charge/observations.',
        bindings={str(p):digest(p) for p in [OUT/'bank.npz',OUT/'inventory.npz',OUT/'spectrum-256-257.npz',OUT/'channels.npz']}))
    signal.alarm(0);print('HHE INPUT',time.monotonic()-begin,checks,flush=True)


def run():
    # The original numerical gates and all six paths remain exactly the same.
    fn=FunctionType(old.run.__code__,dict(vars(old),OUT=OUT),argdefs=old.run.__defaults__)
    fn()
    result=json.loads((OUT/'response.json').read_text())
    result.update(hydrogen_helium_coupled=True,charged_ion_free_free_coupled=True,
        decision='Same-EOS H/HeII bound and charged-ion free-free inputs were actually evolved. Their previous absence is fixed. Other missing channels, pressure/dissolved-state physics and full GR feedback remain; the original complete-input failure is still open.')
    write(OUT/'response.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','build','inputs','run']);globals()[p.parse_args().action]()
