module mod_photon_eos_levels
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
