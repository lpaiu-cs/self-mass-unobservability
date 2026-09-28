"""Reuse the gas EOS's own PL/MHD hydrogenic level sums, without fitting."""
from pathlib import Path
import ctypes
import difflib
import json
import shutil
import subprocess
import time
import numpy as np
import def_photon_eos_populations as p

OUT=p.OUT/'levels-repaired';CACHE=Path('/home/lpaiu/work/direct-eos-gr33/photon-eos-levels-repaired')


def run():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    source=p.SOURCE;parent=source.parent.parent/'build/src'
    files=['qstar_calc.f90','qryd_calc.f90','qryd_approx.f90','plsum.f90','plsum_approx.f90']
    for name in files:
        shutil.copyfile(source/name,OUT/name);shutil.copyfile(source/name,CACHE/name)
    module="""module mod_photon_eos_levels
  use mod_free_eos_types, only: fp_kind
  implicit none
contains
#include "qstar_calc.f90"
#include "qryd_calc.f90"
#include "qryd_approx.f90"
#include "plsum.f90"
#include "plsum_approx.f90"
end module

subroutine photon_levels(z,tl,x,values) bind(C)
  use iso_c_binding
  use mod_photon_eos_levels, only: qstar_calc
  implicit none
  integer(c_int),value :: z
  real(c_double),value :: tl
  real(c_double),intent(in) :: x(5)
  real(c_double),intent(out) :: values(10,3)
  real(c_double) :: q(10),qt(10),qtt(10),qx(5,10),qtx(5,10),qxx(5,5,10)
  integer :: n
  ! The native vector assignment is not a recursive cumulative sum.
  ! Its last slot is the exact tail: request that slot for each n.
  do n=1,10
     call qstar_calc(1,.true.,.true.,z.eq.1,.false.,1.d0,n,10,z,tl,x, &
          q(:n),qt(:n),qx(:,:n),qtt(:n),qtx(:,:n),qxx(:,:,:n))
     values(n,:)=[q(n),qt(n),qtt(n)]
  end do
end subroutine
"""
    (OUT/'bridge.F90').write_text(module);(CACHE/'bridge.F90').write_text(module)
    p.a.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Expose the SAME finite PL/MHD level sum used by the gas EOS for H and He II. Read the exact final cumulative slot for each n; native slice assignment is nonrecursive. Preserve the first failed extraction. Optical survival interpretation remains separate.',
        budget=dict(build_seconds=90,evaluation_seconds=30,native_EOS_calls=0,new_stellar_steps=0,coupled_runs=0),
        gates=dict(native_cumulative_relative=1e-12,positive_levels=True,derivative_relative=1e-5),
        bindings={str(q):p.a.digest(q) for q in [Path(__file__),p.OUT/'eos-state.npz',OUT/'bridge.F90']+[OUT/n for n in files]}))
    cmd=['gfortran','-cpp','-O2','-fPIC','-shared','-fcheck=all','-I'+str(parent),'bridge.F90',
        '-L'+str(parent),'-Wl,-rpath,'+str(parent),'-lfree_eos_direct24_molecular_spectral','-o','levels.so']
    begin=time.monotonic();done=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=90)
    (OUT/'build.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr
    p.a.write(OUT/'build.json',dict(command=cmd,seconds=time.monotonic()-begin,binary_sha256=p.a.digest(CACHE/'levels.so')))
    lib=ctypes.CDLL(str(CACHE/'levels.so'));fn=lib.photon_levels
    array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
    fn.argtypes=[ctypes.c_int,ctypes.c_double,array,array];fn.restype=None
    data=np.load(p.OUT/'eos-state.npz');T=float(data['T'][0]);x=np.ascontiguousarray(data['x'][0])
    rows=[];level_values=[]
    for z in [1,2]:
        def evaluate(lnT):
            out=np.zeros(30);fn(z,lnT,x,out);return out.reshape(3,10)
        base=evaluate(np.log(T));h=2e-4;minus=evaluate(np.log(T)-h);plus=evaluate(np.log(T)+h)
        cumulative=base[0];weights=cumulative-np.r_[cumulative[1:],0.]
        assert np.all(weights>0)
        # Only n>=2 native qstar slots are source-defined at this EOS state.
        native=data['qstar'][0,z-1,1:3]
        native_error=float(np.max(abs(cumulative[1:3]/native-1)))
        d1=(plus[0]-minus[0])/(2*h);d2=(plus[0]-2*base[0]+minus[0])/h**2
        e1=float(np.max(abs(d1/base[1]-1)));e2=float(np.max(abs(d2/base[2]-1)))
        assert native_error<1e-12 and max(e1,e2)<1e-5
        rows.append(dict(Z=z,native_cumulative_relative=native_error,first_derivative_relative=e1,second_derivative_relative=e2))
        level_values.append(base)
    np.savez_compressed(OUT/'level-sums.npz',T=T,x=x,cumulative=np.array(level_values))
    p.a.write(OUT/'result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        physical_optical_occupations_certified=False,
        boundary='This exports exact source-defined thermodynamic weights and partial temperature derivatives at frozen environmental moments. It does not replace total EOS derivatives or license unmodified Einstein absorption/emission with these weights.'))
    print('LEVEL PROVIDER',rows)


if __name__=='__main__':run()
