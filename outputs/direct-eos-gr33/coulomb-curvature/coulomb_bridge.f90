! Request29: expose actual FreeEOS auxiliaries without rebuilding the EOS.
! Imported from prior work: Request21 bridge and unmodified FreeEOS 3.0.0.
subroutine coulomb_inventory(kif, match_value, tl, eps, res, info) bind(C)
  use iso_c_binding
  use mod_free_eos, only: free_eos
  use mod_free_eos_constants, only: c2, cr
  use mod_ionization_data, only: h2diss
  implicit none
  integer(c_int), value :: kif
  real(c_double), value :: match_value, tl
  real(c_double), intent(in) :: eps(24)
  real(c_double), intent(out) :: res(25)
  integer(c_int), intent(out) :: info
  integer :: iterations
  real(c_double) :: fl,t,rho,rl,p,pl,cf,cp,s,sf,st,grada,rtp,qe,qv,rmue, &
       fh2,fhe2,fhe3,xmu1,xmu3,eta,gamma1,gamma2,gamma3,h2rat,h2plusrat, &
       lambda,gamma_e,sound2,pressure(3),density(3),energy(3),entropy(3)
  call free_eos(0,3,1,-2,kif,eps,match_value,tl,fl,t,rho,rl,p,pl, &
       cf,cp,s,sf,st,grada,rtp,qe,qv,rmue,fh2,fhe2,fhe3,xmu1,xmu3,eta, &
       gamma1,gamma2,gamma3,h2rat,h2plusrat,lambda,gamma_e,sound2, &
       iterations,info,pressure=pressure,density=density,energy=energy,entropy=entropy)
  if (info /= 0) return
  res(:12) = [rho,p,qe-0.5d0*c2*cr*h2diss*eps(1),s,gamma1, &
       pressure(2),pressure(3),density(2),density(3),energy(2),energy(3), &
       0.5d0*c2*cr*h2diss*eps(1)]
  res(13:22) = [eta,rmue,fh2,fhe2,fhe3,xmu1,xmu3,lambda,gamma_e,sound2]
  res(23:24) = [h2rat,h2plusrat]
  res(25) = fl
end subroutine coulomb_inventory

subroutine native_constants(res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: avogadro,boltzmann,ct,c_e,cpe
  real(c_double), intent(out) :: res(5)
  res=[avogadro,boltzmann,ct,c_e,cpe]
end subroutine native_constants

subroutine coulomb_probe(fl,tl,sum0,sum2,res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,ct,c_e,cpe
  use mod_fermi_dirac, only: fermi_dirac
  use mod_coulomb, only: master_coulomb
  use mod_master_coulomb_data, only: fcoulomb
  implicit none
  real(c_double), value :: fl,tl,sum0,sum2
  real(c_double), intent(out) :: res(19)
  real(c_double) :: r(9),p(9),ss(3),uu(3),f,t,ne,pe,kt,lam,gamma_e,den, &
       dve,dvef,dvet,dve0,dve2,dv0,dv0f,dv0t,dv2,dv2f,dv2t,dv00,dv02,dv22
  call fermi_dirac(0,fl,tl+log(ct),r,p,ss,uu,21)
  f=exp(fl);t=exp(tl);ne=c_e*r(1);pe=cpe*p(1);kt=boltzmann*t
  call master_coulomb(3,r,f,sum0,0.d0,0.d0,0.d0,sum2,0.d0,0.d0,0.d0, &
       ne,t,pe,p,lam,gamma_e,5,0,0,dve,dvef,dvet,dve0,dve2, &
       dv0,dv0f,dv0t,dv2,dv2f,dv2t,dv00,dv02,dv22)
  den=ne*r(2)
  res=[fcoulomb,-kt*dve,-kt*dv0,-kt*dv2,-kt*dvef/den, &
       -kt*dve0,-kt*dve2,-kt*dv00,-kt*dv02,-kt*dv22, &
       ne,pe,r(2),sqrt(1.d0+f),lam,kt,-kt*dv0f/den,-kt*dv2f/den,r(2)/sqrt(1.d0+f)]
end subroutine coulomb_probe
