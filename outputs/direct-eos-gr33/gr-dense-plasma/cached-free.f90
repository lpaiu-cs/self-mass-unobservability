subroutine cached_free(res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,cpe
  use mod_master_coulomb_data, only: fcoulomb
  use mod_master_exchange_data, only: iforder,fexprime,n_e,t,p_e,pstarprime,psiprime,psi
  implicit none
  real(c_double),intent(out) :: res(5)
  real(c_double) :: fex2
  fex2=boltzmann*t*(psiprime-psi)*n_e-(cpe*pstarprime(1)-p_e)
  res=[fcoulomb,fexprime+fex2,n_e,t,real(iforder,c_double)]
end subroutine cached_free
