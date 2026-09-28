"""Common finite atomic states in chemical equilibrium and excitation free energy.

Counterexample candidate: the inherited PL/MHD model is explicit, not certified
plasma physics. Frozen gas libraries and stellar histories are never modified.
"""
from pathlib import Path
import argparse
import ctypes
import difflib
import json
import re
import shutil
import subprocess
import time
import numpy as np
import def_photon_eos_populations as p

OUT=p.OUT.parent/'def-photon-shared-atomic'
CACHE=Path('/home/lpaiu/work/direct-eos-gr33/photon-shared-atomic')
PARENT=p.SOURCE.parent.parent/'build/src'
NAME='free_eos_shared_atomic'
LIB=CACHE/('lib'+NAME+'.so')


def command(args, name, timeout=90):
    start=time.monotonic()
    done=subprocess.run(args,cwd=CACHE,capture_output=True,text=True,timeout=timeout)
    (OUT/(name+'.log')).write_text(done.stdout+done.stderr)
    assert done.returncode==0,done.stderr
    return dict(command=args,seconds=time.monotonic()-start)


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    p.a.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='8309f9b0b',
        claim='One finite atomic level catalog enters both chemical potentials and excitation free energy, including analytic derivatives, before the optical populations are exported.',
        decision='Accept the same-state EOS/optical provider only if inventory, total thermodynamic derivatives, free-energy identities and detailed balance pass. No stellar or photon time integration is licensed by this local check.',
        budget=dict(cpu_threads=1,build_seconds=120,evaluation_seconds=60,native_EOS_calls=12,new_stellar_steps=0,coupled_runs=0),
        forecast='Previous four same-state calls took 0.209 s. New source build and EOS root speed are unmeasured; hard caps are 120 s and 60 s, no automatic extension.',
        gates=dict(inventory_relative=1e-10,charge_relative=1e-10,derivative_relative=1e-5,free_energy_relative=1e-5,level_sum_relative=1e-12,kirchhoff_relative=1e-10,history_bitwise=True),
        model='Retain native H/He II and molecular states. Replace He I and add the frozen explicit metal catalogs; retain EOS ground energies/weights, preserving imported level excitation energies. Retain native n_eff cutoff 10 and PL/MHD weights, with the native rounded-n_eff radius approximation, no added Rydberg tail. Unknown stages and dissolved continua remain unresolved.',
        original_runtime=json.loads((p.a.OUT.parent/'gr-radiation-eos-split/build.json').read_text())['runtime']))
    constants='''program export_constants
  use mod_free_eos_constants, only: ergspercmm1,rydberg,c2,clight,boltzmann
  use mod_statistical_weight_data, only: iqneutral,iqion
  use mod_ionization_data, only: bi
  implicit none
  print *,ergspercmm1,rydberg,c2,clight,boltzmann
  print *,iqneutral
  print *,iqion
  print *,bi
end program
'''
    (CACHE/'constants.f90').write_text(constants);(OUT/'constants.f90').write_text(constants)
    rec=command(['gfortran','-I'+str(PARENT),'constants.f90','-L'+str(PARENT),'-Wl,-rpath,'+str(PARENT),'-lfree_eos_direct24_molecular_spectral','-o','constants'],'constants-build')
    raw=subprocess.check_output([str(CACHE/'constants')],text=True);(OUT/'constants.txt').write_text(raw)
    c,ground,parent,bi=[np.fromstring(line,sep=' ') for line in raw.splitlines()]
    assert len(c)==5 and len(ground)==24 and len(parent)==316 and len(bi)==318
    import direct_ion_eos as d
    explicit=p.a.OUT/'profile-integral/balanced-fine'
    models=[]
    for line in (explicit/'fort.5').read_text().splitlines()[104:]:
        if "'" not in line:continue
        nums=line.split("'",1)[0].split();z,q,n=map(int,nums[:3])
        if z==0:break
        models.append((z,q,n))
    levels=np.array([[float(x) for x in line.split()[1:]] for line in (explicit/'fort.28').read_text().splitlines() if line.startswith('LEVEL')])
    assert sum(n for z,q,n in models)==len(levels)
    rows=[];offset=0
    for z,q,n in models:
        block=levels[offset:offset+n];offset+=n
        if (z,q) in [(1,0),(2,1)] or n==1:continue
        element=int(np.flatnonzero(d.CHARGES==z)[0]);ion=int(d.CHARGES[:element].sum())+q+1
        g0=ground[element] if q==0 else parent[ion-2]
        assert block[0,4]==g0,(z,q,block[0,4],g0)
        excitation=(block[0,5]-block[:,5])/c[0]
        binding=bi[ion-1]-excitation
        assert np.all(excitation>=0) and np.all(binding>0),(z,q)
        neff=np.sqrt(c[1]*(q+1)**2/binding)
        included=np.floor(neff+0.5)<=10
        assert included[0]
        rows.append(dict(Z=z,charge=q,element_index=element+1,ion_index=ion,
            ground_weight=float(g0),parent_weight=float(parent[ion-1]),ground_cm_inverse=float(bi[ion-1]),
            imported_ground_erg=float(block[0,5]),anchor_shift_erg=float(bi[ion-1]*c[0]-block[0,5]),
            global_levels=block[:,0].astype(int).tolist(),weight=block[:,4].tolist(),
            excitation_cm_inverse=excitation.tolist(),binding_cm_inverse=binding.tolist(),
            neff=neff.tolist(),included=included.tolist()))
    assert len(rows)==16,len(rows)
    p.a.write(OUT/'catalog.json',dict(classification='Counterexample candidate',constants=c.tolist(),rows=rows,
        energy_rule='Imported excitation energies in erg are retained; the stage threshold is anchored to the original EOS ionization energy. This explicit shift is not an atomic-data calibration.',
        radius_rule='qmhd_calc ignores its ell array and uses rounded neff and ellapprox=n-1. Reuse exactly; do not infer LS term L is a one-electron ell.',
        bindings={str(f):p.a.digest(f) for f in [explicit/'fort.5',explicit/'fort.28',p.a.OUT/'inputs.tar.gz',OUT/'constants.txt']}))
    p.a.write(OUT/'prepare.json',dict(constants_build=rec,stages=len(rows),excited_levels=sum(sum(r['included'][1:]) for r in rows)))
    print('CATALOG',[(r['Z'],r['charge'],len(r['weight']),sum(r['included'])) for r in rows],flush=True)


def build():
    assert not (OUT/'build.json').exists()
    catalog=json.loads((OUT/'catalog.json').read_text());rows=catalog['rows'];begin=time.monotonic()
    starts=np.zeros(318,dtype=int);ends=starts.copy();weights=[];neff=[]
    for r in rows:
        i=r['ion_index']-1;starts[i]=len(weights)+1
        for j in range(1,len(r['weight'])):
            if r['included'][j]:weights.append(r['weight'][j]);neff.append(r['neff'][j])
        ends[i]=len(weights)
    def arr(name,values,integer=False):
        typ='integer' if integer else 'real(fp_kind)'
        vals=[str(int(x)) if integer else f'{x:.17e}_fp_kind' for x in values]
        return f'  {typ}, parameter :: {name}({len(vals)}) = [ &\n'+', &\n'.join('    '+x for x in vals)+' ]\n'
    declaration=arr('shared_start',starts,True)+arr('shared_end',ends,True)+arr('shared_weight',weights)+arr('shared_neff',neff)
    helper='''logical function shared_has(ion)
  integer,intent(in) :: ion
  shared_has=shared_start(ion)>0
end function

subroutine shared_atomic(ion,tl,x,q,qt,qx,qtt,qtx,qxx)
  use mod_ionization_data, only: nion
  use mod_statistical_weight_data, only: iqion
  integer,intent(in) :: ion
  real(fp_kind),intent(in) :: tl,x(5)
  real(fp_kind),intent(out) :: q,qt,qx(5),qtt,qtx(5),qxx(5,5)
  integer :: lo,hi
  lo=shared_start(ion);hi=shared_end(ion)
  if(lo<1.or.hi<lo) error stop 'shared atomic catalog missing'
  qx=0; qtx=0; qxx=0
  ! Finite imported levels only. Native radius approximation ignores ell.
  call qmhd_calc(1,1._fp_kind,.true.,.false.,1._fp_kind,11,10,tl,nion(ion),1._fp_kind, &
       shared_weight(lo:hi),shared_neff(lo:hi),0._fp_kind*shared_neff(lo:hi),x,q,qt,qx,qtt,qtx,qxx)
  q=q/iqion(ion);qt=qt/iqion(ion);qx=qx/iqion(ion)
  qtt=qtt/iqion(ion);qtx=qtx/iqion(ion);qxx=qxx/iqion(ion)
end subroutine

subroutine shared_terms(ion,tl,x,values) bind(C)
  use iso_c_binding
  use mod_ionization_data, only: nion
  integer(c_int),value :: ion
  real(c_double),value :: tl
  real(c_double),intent(in) :: x(5)
  real(c_double),intent(out) :: values(38,700)
  real(fp_kind) :: q,qt,qtt,qx(5),qtx(5),qxx(5,5)
  integer :: j,k
  values=0
  if(.not.shared_has(ion)) error stop 'shared optical catalog missing'
  do j=shared_start(ion),shared_end(ion)
     qx=0;qtx=0;qxx=0
     call qmhd_calc(1,1._fp_kind,.true.,.false.,1._fp_kind,11,10,tl,nion(ion),1._fp_kind, &
          shared_weight(j:j),shared_neff(j:j),0._fp_kind*shared_neff(j:j),x,q,qt,qx,qtt,qtx,qxx)
     k=j-shared_start(ion)+1
     values(:,k)=[q,qt,qtt,qx,qtx,reshape(qxx,[25])]
  end do
end subroutine
'''
    (OUT/'shared-data.f90').write_text(declaration);(OUT/'shared-functions.f90').write_text(helper)
    module=(p.SOURCE/'mod_excitation.f90').read_text().replace('  private','  private\n'+declaration).replace('        excitation_sum','        excitation_sum, shared_terms').replace('contains','contains\n'+helper)
    (CACHE/'mod_excitation.f90').write_text(module);(OUT/'mod_excitation.f90').write_text(module)
    diffs=[]
    files=re.findall(r'#include "([^"]+)"',module)
    init='''
           if(shared_has(ion)) then
              if(.not.ifmhd_logical.or..not.ifpl_logical.or.ifpi_fit.ne.1.or.nmax.ne.10) &
                   error stop 'shared atomic requires frozen PL/MHD options'
              call shared_atomic(ion,tl,x,qa,qat,qax,qatt,qatx,qaxx)
           else
              qa=qstar(nmin_s,iz);qat=qstart(nmin_s,iz);qatt=qstart2(nmin_s,iz)
              qax=qstarx(:,nmin_s,iz);qatx=qstartx(:,nmin_s,iz);qaxx=qstarx2(:,:,nmin_s,iz)
           endif
'''
    for name in files:
        before=(p.SOURCE/name).read_text();after=before
        if name in ['excitation_sum.f90','excitation_pi.f90']:
            after=after.replace('  ! Local variables','  real(fp_kind) :: qa,qat,qatt,qax(5),qatx(5),qaxx(5,5)\n\n  ! Local variables',1)
            start=after.index('  ! This do loop is over atomic species only')
            end=after.index('     ! do all the molecular stuff below',start)
            section=after[start:end]
            section=section.replace('ifhe1_special.and.ion.eq.2','ifhe1_special.and.ion.eq.2.and..not.shared_has(ion)')
            for old,new in [('qstar(nmin_s,iz)','qa'),('qstart(nmin_s,iz)','qat'),('qstart2(nmin_s,iz)','qatt')]:section=section.replace(old,new)
            for old,new in [('qstarx','qax'),('qstartx','qatx'),('qstarx2','qaxx')]:
                section=re.sub(old+r'\(([^()]*)?,nmin_s,iz\)',lambda m:new+'('+m[1]+')',section)
            if name=='excitation_sum.f90':
                section=section.replace('(ielement.le.2.or.mod(ifexcited,10).eq.3)) then','(ielement.le.2.or.mod(ifexcited,10).eq.3.or.any(shared_start>0))) then',1)
                section=section.replace('         do ion = ion_start, ion_end(ielement)','         do ion = ion_start, ion_end(ielement)\n            if(ielement.gt.2.and.mod(ifexcited,10).ne.3.and..not.shared_has(ion)) cycle',1)
            else:
                section=section.replace('if(ielement.le.2.or.mod(ifexcited,10).eq.3) then','if(ielement.le.2.or.mod(ifexcited,10).eq.3.or.shared_has(ion)) then',1)
            anchor=next(line for line in section.splitlines() if "invalid nmin or nmin_max'" in line)
            section=section.replace(anchor,anchor+'\n'+init,1)
            after=after[:start]+section+after[end:]
            assert 'call shared_atomic' in after and 'real(fp_kind) :: qa' in after
        (CACHE/name).write_text(after);(OUT/name).write_text(after)
        if before!=after:diffs.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original/'+name,tofile='candidate/'+name))
    (OUT/'source.patch').write_text(''.join(diffs))
    compile_info=command(['gfortran','-cpp','-O2','-fPIC','-fcheck=all','-ffree-line-length-none','-I'+str(PARENT),'-c','mod_excitation.f90','-o','mod_excitation.o'],'module-build')
    objects=sorted(PARENT.glob('CMakeFiles/*free_eos.dir/*.o'));assert len(objects)>=25
    objects=[q for q in objects if q.name!='mod_excitation.f90.o']
    link_info=command(['gfortran','-shared','-Wl,-soname,'+LIB.name,'mod_excitation.o',*map(str,objects),'-llapack','-lblas','-o',str(LIB)],'library-link')
    bridge=p.a.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    bridge_info=command(['gfortran','-O2','-fPIC','-shared','-I'+str(PARENT),str(bridge),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-l'+NAME,'-o','gas.so'],'bridge-build')
    assert time.monotonic()-begin<120
    p.a.write(OUT/'build.json',dict(classification='Counterexample candidate',compile=compile_info,link=link_info,bridge=bridge_info,
        seconds=time.monotonic()-begin,runtime={str(f):p.a.digest(f) for f in [LIB,CACHE/'gas.so']},
        reused_objects={str(f):p.a.digest(f) for f in objects},
        source_bindings={str(p.SOURCE/n):p.a.digest(p.SOURCE/n) for n in files+['mod_excitation.f90']},
        candidate_source_bindings={str(OUT/n):p.a.digest(OUT/n) for n in files+['mod_excitation.f90']}))
    print('BUILT',time.monotonic()-begin,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','build'])
    globals()[parser.parse_args().action]()
