! Request29: expose actual FreeEOS auxiliaries without rebuilding the EOS.
! Imported from prior work: Request21 bridge and unmodified FreeEOS 3.0.0.
subroutine ionization_inventory(kif, match_value, tl, eps, res, info) bind(C)
  use iso_c_binding
  use mod_free_eos, only: free_eos
  use mod_free_eos_constants, only: c2, cr
  use mod_ionization_data, only: h2diss
  implicit none
  integer(c_int), value :: kif
  real(c_double), value :: match_value, tl
  real(c_double), intent(in) :: eps(24)
  real(c_double), intent(out) :: res(24)
  integer(c_int), intent(out) :: info
  integer :: iterations
  real(c_double) :: fl,t,rho,rl,p,pl,cf,cp,s,sf,st,grada,rtp,qe,qv,rmue, &
       fh2,fhe2,fhe3,xmu1,xmu3,eta,gamma1,gamma2,gamma3,h2rat,h2plusrat, &
       lambda,gamma_e,sound2,pressure(3),density(3),energy(3),entropy(3)
  call free_eos(0,3,11,-2,kif,eps,match_value,tl,fl,t,rho,rl,p,pl, &
       cf,cp,s,sf,st,grada,rtp,qe,qv,rmue,fh2,fhe2,fhe3,xmu1,xmu3,eta, &
       gamma1,gamma2,gamma3,h2rat,h2plusrat,lambda,gamma_e,sound2, &
       iterations,info,pressure=pressure,density=density,energy=energy,entropy=entropy)
  if (info /= 0) return
  res(:12) = [rho,p,qe-0.5d0*c2*cr*h2diss*eps(1),s,gamma1, &
       pressure(2),pressure(3),density(2),density(3),energy(2),energy(3), &
       0.5d0*c2*cr*h2diss*eps(1)]
  res(13:22) = [eta,rmue,fh2,fhe2,fhe3,xmu1,xmu3,lambda,gamma_e,sound2]
  res(23:24) = [h2rat,h2plusrat]
end subroutine ionization_inventory

subroutine photon_constants(values) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: prad_const,clight
  implicit none
  real(c_double),intent(out) :: values(2)
  values=[3._c_double*prad_const,clight]
end subroutine photon_constants
