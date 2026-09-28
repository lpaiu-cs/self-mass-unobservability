! Extended arithmetic candidate. No physical EOS or native error certificate.
subroutine direct_ion_eos(kif, vh, vl, th, tlow, cx, eps, hi, lo, info) bind(C)
  use iso_c_binding
  use mod_free_eos_types, only: fp_kind
  use mod_free_eos, only: free_eos
  use mod_free_eos_constants, only: c2, cr
  use mod_ionization_data, only: h2diss
  implicit none
  integer(c_int), value :: kif
  real(c_double), value :: vh, vl, th, tlow, cx
  real(c_double), intent(in) :: eps(24)
  real(c_double), intent(out) :: hi(22),lo(22)
  integer(c_int), intent(out) :: info
  integer :: iterations
  real(fp_kind) :: match_value,tl,epsq(24),res(22),cxq
  real(fp_kind) :: fl,t,rho,rl,p,pl,cf,cp,s,sf,st,grada,rtp,qe,qv,rmue, &
       fh2,fhe2,fhe3,xmu1,xmu3,eta,gamma1,gamma2,gamma3,h2rat,h2plusrat, &
       lambda,gamma_e,sound2,pressure(3),density(3),energy(3),entropy(3)
  match_value=real(vh,fp_kind)+real(vl,fp_kind)
  tl=real(th,fp_kind)+real(tlow,fp_kind)
  epsq=real(eps,fp_kind);cxq=real(cx,fp_kind)
  ! The disabled diffraction output is padding, discarded by the Python caller.
  gamma_e=0._fp_kind
  call free_eos(0,3,1,-2,kif,epsq,match_value,tl,fl,t,rho,rl,p,pl, &
       cf,cp,s,sf,st,grada,rtp,qe,qv,rmue,fh2,fhe2,fhe3,xmu1,xmu3,eta, &
       gamma1,gamma2,gamma3,h2rat,h2plusrat,lambda,gamma_e,sound2, &
       iterations,info,pressure=pressure,density=density,energy=energy,entropy=entropy)
  if (info /= 0) return
  res(:12) = [rho,p,qe-0.5_fp_kind*c2*cr*h2diss*epsq(1),s,gamma1, &
       pressure(2),pressure(3),density(2),density(3),energy(2),energy(3), &
       0.5_fp_kind*c2*cr*h2diss*epsq(1)]
  res(13:) = [eta,rmue,fh2,fhe2,fhe3,xmu1,xmu3,lambda,gamma_e,sound2]
  res(1)=res(1)/cxq;res(3:4)=res(3:4)*cxq;res(10:11)=res(10:11)*cxq
  hi=real(res,c_double);lo=real(res-real(hi,fp_kind),c_double)
end subroutine direct_ion_eos

subroutine precision_info(values) bind(C)
  use iso_c_binding
  use mod_free_eos_types, only: fp_kind,lapack_fp_kind
  implicit none
  integer(c_int),intent(out) :: values(5)
  values=[precision(1._fp_kind),digits(1._fp_kind),range(1._fp_kind), &
       storage_size(1._fp_kind),digits(1._lapack_fp_kind)]
end subroutine precision_info
