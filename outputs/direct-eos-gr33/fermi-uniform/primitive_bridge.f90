subroutine primitive(nu,eta,beta,answer) bind(C)
  use iso_c_binding
  use mod_fermi_dirac, only: fermi_dirac_direct
  implicit none
  real(c_double), value :: nu,eta,beta
  real(c_double), intent(out) :: answer(9)
  call fermi_dirac_direct(nu,eta,beta,answer)
end subroutine primitive
