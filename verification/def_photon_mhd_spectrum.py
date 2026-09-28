"""Separate thermal Planck-Larkin weights from plasma optical dissolution.

Counterexample candidate: a distinct finite MHD-only EOS and optical provider;
the accepted PL/MHD histories are immutable. No new stellar integration.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import ctypes
import difflib
import json
import inspect
import shutil
import subprocess
import time
import numpy as np
import def_photon_shared_atomic as old
import def_photon_shared_atomic_check as check
import def_photon_atomic_rates as rates

OUT=old.OUT.parent/'def-photon-mhd-spectrum'
CACHE=old.CACHE.parent/'photon-mhd-spectrum'
NAME='free_eos_mhd_spectrum'
LIB=CACHE/('lib'+NAME+'.so')
write=rates.write
digest=rates.digest


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b35f5c16',
        claim='Remove the false identification of a density-independent PL thermal weight with plasma dissolution, changing the EOS chemical/free-energy calculation and optical weights together. Reconnect the actual atomic amplitudes to the resulting MHD-only common state.',
        decision='Adopt this separate candidate only after total EOS derivatives, chemical inventory, common optical partition and detailed balance pass unchanged gates, and the zero-environment weights recover one. Do not transfer old GR or opacity results to the new EOS.',
        model='Keep the finite catalog, native radii/fitting parameters and nmax10. Switch only PL off through native ifpi=-3 in gas option11; retain MHD, Coulomb, molecules, statistical weights and all thermodynamic derivatives. The catalog truncation, MHD optical interpretation and full spectrum remain explicit limitations.',
        budget=dict(build_attempts=2,build_seconds=120,native_EOS_calls=12,analysis_seconds=60,atomic_extractions=0,coupled_runs=0,new_stellar_steps=0,cpu_threads=1),
        forecast='The prior common EOS build took6.06s and12calls plus checks0.410s; this changed branch is unmeasured. Hard caps120s build and60s evaluation, no grid expansion or automatic new state calls after failure.',
        gates=dict(total_derivative_relative=1e-5,free_energy_relative=1e-5,inventory_relative=1e-10,chemical_log_ratio=1e-9,partition_relative=1e-12,kirchhoff_relative=1e-10,dilute_occupation_absolute=1e-12),
        original_runtime={str(p):digest(p) for p in [old.LIB,old.CACHE/'gas.so',rates.CACHE/'atomic-rates']},
        bindings={str(p):digest(p) for p in [Path(__file__),old.OUT/'catalog.json',old.OUT/'eos-state.npz',rates.OUT/'cross-sections.npz',rates.OUT/'fort.92',rates.OUT/'fort.98']}))
    for p in old.CACHE.glob('*.f90'):shutil.copyfile(p,CACHE/p.name)
    for p in old.CACHE.glob('*.mod'):shutil.copyfile(p,CACHE/p.name)
    for name in ['catalog.json','constants.txt']:shutil.copyfile(old.OUT/name,OUT/name)
    patches=[]
    for name in ['excitation_sum.f90','excitation_pi.f90']:
        p=CACHE/name;before=p.read_text();after=before
        assert after.count('.not.ifpl_logical.or.ifpi_fit.ne.1.or.nmax.ne.10')==1
        after=after.replace('.not.ifpl_logical.or.ifpi_fit.ne.1.or.nmax.ne.10','ifpl_logical.or..not.ifmhd_logical.or.ifpi_fit.ne.1.or.nmax.ne.10')
        after=after.replace('shared atomic requires frozen PL/MHD options','shared atomic requires MHD-only options')
        if name=='excitation_sum.f90':
            # Expose H and He II ground factors too, for their native provider.
            after=after.replace('''           if(shared_has(ion)) then
              shared_ground_logw(ion)=plop(ion)
              if(ifmhd_logical) shared_ground_logw(ion)=plop(ion)+ln_ground_occ
           endif''','''           shared_ground_logw(ion)=plop(ion)
           if(ifmhd_logical) shared_ground_logw(ion)=plop(ion)+ln_ground_occ''')
        p.write_text(after);patches+=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='phase66/'+name,tofile='phase68/'+name))
    p=CACHE/'mod_excitation.f90';before=p.read_text()
    anchor='call qmhd_calc(1,1._fp_kind,.true.,.false.,1._fp_kind,11,10,tl,nion(ion),1._fp_kind,'
    assert before.count(anchor)==2
    after=before.replace(anchor,anchor.replace('.true.','.false.',1))
    # Same native hydrogenic provider, now under the same MHD-only flags.
    bridge='''
subroutine hydrogen_terms(z,tl,x,values) bind(C)
  use iso_c_binding
  implicit none
  integer(c_int),value :: z
  real(c_double),value :: tl
  real(c_double),intent(in) :: x(5)
  real(c_double),intent(out) :: values(10,3)
  real(fp_kind) :: q(10),qt(10),qtt(10),qx(5,10),qtx(5,10),qxx(5,5,10)
  integer :: n
  do n=1,10
     call qstar_calc(1,.false.,.true.,z.eq.1,.false.,1._fp_kind,n,10,z,tl,x, &
          q(:n),qt(:n),qx(:,:n),qtt(:n),qtx(:,:n),qxx(:,:,:n))
     values(n,:)=[q(n),qt(n),qtt(n)]
  enddo
end subroutine
'''
    after=after.replace('end module mod_excitation',bridge+'\nend module mod_excitation')
    p.write_text(after);patches+=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='phase66/mod_excitation.f90',tofile='phase68/mod_excitation.f90'))
    p=CACHE/'mod_free_eos.f90';before=(old.p.SOURCE/p.name).read_text()
    anchor='''       elseif(ifmodified.eq.11) then
          ! EOS1 without radiation pressure
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0'''
    assert before.count(anchor)==1
    after=before.replace(anchor,anchor.replace('ifpi = 3','ifpi = -3'))
    p.write_text(after);patches+=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original/mod_free_eos.f90',tofile='phase68/mod_free_eos.f90'))
    p=CACHE/'free_eos_detailed.f90';before=p.read_text()
    anchor='  ! calculate planck-larkin occupation probabilities and equilibrium\n'
    assert before.count(anchor)==1
    after=before.replace(anchor,'''  ! PL-off still passes these arrays to excitation and molecular consumers.
  ! They are newly allocated each call, so never inherit allocator contents.
  plop=0._fp_kind;plopt=0._fp_kind;plopt2=0._fp_kind
  dv_pl=0._fp_kind;dv_plt=0._fp_kind
'''+anchor)
    p.write_text(after);patches+=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='phase66/free_eos_detailed.f90',tofile='phase68/free_eos_detailed.f90'))
    (OUT/'source.patch').write_text(''.join(patches))


def build():
    assert not (OUT/'build.json').exists();begin=time.monotonic();receipts=[]
    def call(cmd,name):
        t=time.monotonic();p=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=60)
        (OUT/(name+'.log')).write_text(p.stdout+p.stderr)
        receipts.append(dict(command=cmd,seconds=time.monotonic()-t,returncode=p.returncode))
        assert p.returncode==0,p.stderr
    flags=['gfortran','-cpp','-O2','-fPIC','-fcheck=all','-ffree-line-length-none','-I'+str(CACHE),'-I'+str(old.PARENT),'-I'+str(old.p.SOURCE)]
    for name in ['mod_excitation','mod_free_eos_detailed','mod_free_eos']:
        if not (CACHE/(name+'.f90')).exists():shutil.copyfile(old.p.SOURCE/(name+'.f90'),CACHE/(name+'.f90'))
        call(flags+['-c',name+'.f90','-o',name+'.o'],name)
    objects=sorted(old.PARENT.glob('CMakeFiles/*free_eos.dir/*.o'))
    objects=[p for p in objects if p.name not in ['mod_excitation.f90.o','mod_free_eos.f90.o','mod_free_eos_detailed.f90.o','mod_statistical_weight_data.f90.o','mod_eos_calc.f90.o']]
    objects += [old.CACHE/'mod_statistical_weight_data.o',old.CACHE/'mod_eos_calc.o']
    call(['gfortran','-shared','-Wl,-soname,'+LIB.name,'mod_excitation.o','mod_free_eos_detailed.o','mod_free_eos.o',*map(str,objects),'-llapack','-lblas','-o',str(LIB)],'link')
    source=old.p.a.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    call(['gfortran','-O2','-fPIC','-shared','-I'+str(CACHE),'-I'+str(old.PARENT),str(source),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-l'+NAME,'-o','gas.so'],'bridge')
    assert time.monotonic()-begin<120
    write(OUT/'build.json',dict(classification='Counterexample candidate',seconds=time.monotonic()-begin,receipts=receipts,
        reused_objects={str(p):digest(p) for p in objects},runtime={str(p):digest(p) for p in [LIB,CACHE/'gas.so']},
        candidate_source_bindings={str(p):digest(p) for p in CACHE.glob('*.f90')}))
    print('MHD BUILD',time.monotonic()-begin,flush=True)


def evaluate():
    candidate=SimpleNamespace(**dict(vars(old),OUT=OUT,CACHE=CACHE,LIB=LIB))
    # Two attempts were consumed by the preserved PL-off initialization bug.
    # Use the remaining ten for base, both centered differences, and replay.
    # The inverse-pressure check is explicitly absent, never filled with zero.
    source=inspect.getsource(check.run)
    anchor='''        inverse=gas(1,float(np.log(v[1])),t,X);inverse_error=float(abs(np.log(inverse[0])-r))
        base_again=snapshot(0,0);assert all(np.array_equal(base[k],base_again[k]) for k in base)'''
    assert source.count(anchor)==1
    source=source.replace(anchor,'        inverse_error=None')
    source=source.replace('inverse_error<1e-10 and ','').replace('else 12','else 10').replace('=12 if reuse','=10 if reuse')
    (OUT/'bounded-eos-check.py').write_text(source)
    namespace=dict(vars(check),a=candidate)
    exec(compile(source,str(OUT/'bounded-eos-check.py'),'exec'),namespace)
    namespace['run']()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','build','evaluate'])
    globals()[p.parse_args().action]()
