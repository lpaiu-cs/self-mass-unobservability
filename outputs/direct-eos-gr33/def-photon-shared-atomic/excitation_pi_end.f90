!*******************************************************************************
!       Copyright (C) 1996-2022 Alan W. Irwin
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!*******************************************************************************

!> This excitation_pi_end subroutine calculates the change to the
!> pressure, entropy per unit mass, and energy per unit mass due to
!> the excited state pressure ionization component of the free energy
!> as well as partial derivatives of those first two quantities wrt fl
!> and tl.
!>
!> \param[in] t PARAMETERS NEED DOCUMENTATION
!>
subroutine excitation_pi_end(t, rho, rf, rt,&
     pexcited, pexcitedf, pexcitedt,&
     sexcited, sexcitedf, sexcitedt, uexcited)

  use mod_free_eos_constants, only: cr
  ! These variables all calculated by excitation_sum.
  use mod_excitation_block, only:&
       psum, psumf, psumt,&
       ssum, ssumf, ssumt, usum,&
       ifnr03, excitation_sum_called

  ! Arguments
  real(fp_kind), intent(in) :: t, rho, rf, rt
  real(fp_kind), intent(out) ::&
       pexcited, pexcitedf, pexcitedt,&
       sexcited, sexcitedf, sexcitedt, uexcited

  ! No internal variables

  ! Sanity check
  if(.not.excitation_sum_called) error stop 'excitation_pi_end: excitation_sum must be called first'

  ! free energy  = -kTV sum n(species) delta ln Z
  ! psum is kept in n/(rho*avogardro) form.
  pexcited = cr*t*psum*rho
  ! ssum is kept in n/(rho*avogardro) form above so sexcited
  ! is per unit mass.
  sexcited = cr*ssum
  ! E = -T^2 partial (F/T) wrt T -->
  ! energy/unit mass = R T sum (n/(rho*avogadro))
  ! partial delta ln Z wrt ln t
  uexcited = cr*t*usum
  if(ifnr03) then
     pexcitedf = cr*t*psumf*rho + pexcited*rf
     pexcitedt = cr*t*psumt*rho + pexcited*(rt+1._fp_kind)
     sexcitedf = cr*ssumf
     sexcitedt = cr*ssumt
  endif
end subroutine excitation_pi_end

!> This excitation_pi_pressure_free subroutine calculates the change to the
!> pressure and free energy per unit mass due to
!> the excited state pressure ionization component of the free energy
!> as well as partial derivatives of those quantities
!> wrt fl, tl, and dv.
!>
!> \param[in] t PARAMETERS NEED DOCUMENTATION
!>
subroutine excitation_pi_pressure_free(&
       t, rho, rf, rt, r_dv,&
       pexcited, pexcitedf, pexcitedt, pexcited_dv,&
       free_excited, free_excitedf, free_excited_dv)

  use mod_free_eos_constants, only: cr
  ! These variables all calculated by excitation_sum.
  use mod_excitation_block, only:&
       psum, psumf, psumt, psum_dv,&
       free_sum, free_sumf, free_sum_dv,&
       ifnr03, ifnr13, excitation_sum_called

  ! Arguments
  real(fp_kind), intent(in) :: t, rho, rf, rt, r_dv(:)
  real(fp_kind), intent(out) ::&
       pexcited,  pexcitedf, pexcitedt, pexcited_dv(:),&
       free_excited,  free_excitedf, free_excited_dv(:)

  ! Internal variable
  integer max_index

  max_index = size(pexcited_dv)

  ! Sanity checks
  if(.not.excitation_sum_called) error stop 'excitation_pi_pressure_free: excitation_sum must be called first'
  if(&
       max_index.gt.size(psum_dv).or.&
       max_index.gt.size(free_sum_dv).or.&
       max_index.ne.size(r_dv).or.&
       max_index.ne.size(free_excited_dv))&
       error stop 'excitation_pi_pressure_free: inconsistent sizes for pexcited_dv, psum_dv, r_dv, free_excited_dv,or free_sum_dv'

  pexcited = cr*t*psum*rho
  ! free_sum kept in n/(rho*avogardro) form so free_excited
  ! is per unit mass.
  free_excited = -cr*t*free_sum
  if(ifnr03) then
     pexcitedf = cr*t*psumf*rho + pexcited*rf
     pexcitedt = cr*t*psumt*rho + pexcited*(rt+1._fp_kind)
     free_excitedf = -cr*t*free_sumf
  endif
  if(ifnr13) then
     pexcited_dv = cr*t*rho*(psum_dv(1:max_index) + psum*r_dv)
     free_excited_dv = -cr*t*free_sum_dv(1:max_index)
  endif
end subroutine excitation_pi_pressure_free
