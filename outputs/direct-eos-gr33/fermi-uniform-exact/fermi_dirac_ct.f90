!*******************************************************************************
!    Copyright (C) 1996-2022 Alan W. Irwin
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

! calculate nderiv derivative of fermi-dirac integral of order 3/2.
! derived from cody and thacher 1967 math. comp. 21, 30.
! n.b. see errata for table iic 1967 math. comp. 21, 525!!!!!!!!!
! according to c-t paper, worst precision of these approximations is one
! part in 10^(8.58).  Errors in derivatives will be worse (probably
! something like an order of magnitude per order of derivative).
! This routine also tested against
! eggleton, faulkner, and flannery, generalized
! fermi-dirac integrals and also tables in back of cox and guili vol. 2.
! if(ifsimple == 1) use low degeneracy approximation.

!> This fermi_dirac_ct function calculates the nderiv derivative of
!> the fermi-dirac integral of order 3/2 using the Cody-Thacher
!> approximation (1967 math. comp. 21, 30 plus errata in
!> math. comp. 21, 525).
!>
!> \param[in] eta PARAMETERS NEED DOCUMENTATION
!>
function fermi_dirac_ct(eta, nderiv, ifsimple)
  use mod_free_eos_constants, only: pi
  real(fp_kind) fermi_dirac_ct
  real(fp_kind), intent(in) :: eta
  integer, intent(in) :: nderiv, ifsimple
  real(fp_kind), save :: last_eta = -huge(1._fp_kind), cached(0:5) = 0._fp_kind
  logical, save :: have_zero = .false., have_middle = .false., have_five = .false.
  real(fp_kind) answer(7)
  if(nderiv.lt.0.or.nderiv.gt.5) error stop 'direct CT replacement: bad derivative order'
  if(ifsimple.eq.1) then
     fermi_dirac_ct = 0.75_fp_kind*sqrt(pi)*exp(eta)
     return
  endif
  ! The existing EOS library is serial/stateful; this cache has the same scope.
  if(eta.ne.last_eta) then
     last_eta = eta
     have_zero = .false.
     have_middle = .false.
     have_five = .false.
  endif
  if(nderiv.eq.0.and..not.have_zero) then
     call fermi_dirac_direct(1.5_fp_kind,eta,0._fp_kind,answer(:1))
     cached(0) = answer(1)
     have_zero = .true.
  elseif(1.le.nderiv.and.nderiv.le.4.and..not.have_middle) then
     call fermi_dirac_direct(0.5_fp_kind,eta,0._fp_kind,answer)
     cached(1:4) = 1.5_fp_kind*[answer(1),answer(2),answer(4),answer(7)]
     have_middle = .true.
  elseif(nderiv.eq.5.and..not.have_five) then
     call fermi_dirac_direct(-0.5_fp_kind,eta,0._fp_kind,answer)
     cached(5) = 0.75_fp_kind*answer(7)
     have_five = .true.
  endif
  fermi_dirac_ct = cached(nderiv)
end function fermi_dirac_ct

! For eta greater than 4 calculate 2.5 F_3/2 - eta*F'3/2 with analytical
! cancellation of differences.

!> This fermi_dirac_ct_diff function calculates the nderiv derivative
!> of 2.5 F_3/2 - eta*F'3/2 with analytical cancellation of
!> differences where F_3/2 is the fermi-dirac integral of order 3/2
!> that is calculated using the Cody-Thacher approximation (1967
!> math. comp. 21, 30 plus errata in math. comp. 21, 525).
!>
!> \param[in] eta PARAMETERS NEED DOCUMENTATION
!>
function fermi_dirac_ct_diff(eta, nderiv)

  use mod_poly_sum, only: poly_sum
  use mod_fermi_dirac_ct_data, only: maxderiv, p, q
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Function and arguments
  real(fp_kind) fermi_dirac_ct_diff
  real(fp_kind), intent(in) :: eta
  integer, intent(in) :: nderiv

  ! Internal variables
  integer itab, ideriv

  real(fp_kind) y, exp4, const4

  real(fp_kind), allocatable ::&
       top(:),&
       bot(:),&
       ratio1(:),&
       ratio2(:),&
       ratio3(:),&
       ratio4(:),&
       ratio5(:)

  ! Sanity checks
  if(eta.le.4._fp_kind) error stop 'fermi_dirac_ct_diff: bad eta value'
  if(nderiv.lt.0.or.nderiv.gt.maxderiv-1) error stop 'illegal nderiv argument to fermi_dirac_ct_diff'

  allocate(&
       top(0:maxderiv),&
       bot(0:maxderiv),&
       ratio1(0:maxderiv),&
       ratio2(0:maxderiv),&
       ratio3(0:maxderiv),&
       ratio4(0:maxderiv),&
       ratio5(0:maxderiv))

  if(if_taint_allocated_real) then
     call taint_allocated_real(top)
     call taint_allocated_real(bot)
     call taint_allocated_real(ratio1)
     call taint_allocated_real(ratio2)
     call taint_allocated_real(ratio3)
     call taint_allocated_real(ratio4)
     call taint_allocated_real(ratio5)
  endif

  ! Use cody and thacher table ivc
  itab = 3
  y = 1._fp_kind/(eta*eta)
  top(0) = poly_sum(y,p(0:4,0,itab))
  bot(0) = poly_sum(y,q(0:4,0,itab))
  ratio2(0) = deriv_ratio(top(:0),bot(:0))
  do ideriv = 0,nderiv
     if(ideriv.eq.0) then
        ratio1(ideriv) = y
        exp4 = -0.25_fp_kind
        const4 = 2._fp_kind
     elseif(ideriv.eq.1) then
        ratio1(ideriv) = 1._fp_kind
        const4 = const4*exp4
        exp4 = exp4 - 1._fp_kind
     else
        ratio1(ideriv) = 0._fp_kind
        const4 = const4*exp4
        exp4 = exp4 - 1._fp_kind
     endif
     ratio4(ideriv) = const4*y**exp4
     top(ideriv+1) = poly_sum(y,p(0:3-ideriv,ideriv+1,itab))
     bot(ideriv+1) = poly_sum(y,q(0:3-ideriv,ideriv+1,itab))
     ratio2(ideriv+1) = deriv_ratio(top(:ideriv+1),bot(:ideriv+1))
     ratio3(ideriv) = ratio2(ideriv) + deriv_product(ratio1(:ideriv),ratio2(1:ideriv+1))
     ratio5(ideriv) = deriv_product(ratio4(:ideriv),ratio3(:ideriv))
  enddo
  ! ratio1 is y and its x derivatives.
  ratio1(0) = y
  do ideriv = 1,nderiv
     ratio1(ideriv) = -real((ideriv+1),fp_kind)*ratio1(ideriv-1)/eta
  enddo
  ! transform from y derivatives to x derivatives.
  ! fermi_dirac_ct_diff is *always* defined below (on assumption that 0 <= nderiv <= maxderiv-1 = 4 which
  ! is confirmed by the sanity check above).
  ! Unfortunately, gfortran is unable to figure out this
  ! initialization logic so it emits a spurious -Wmaybe-uninitialized
  ! warning when -Wall is an option.
  ! Suppress that warning with the following line of code.
  fermi_dirac_ct_diff = 0._fp_kind
  if(nderiv.eq.0) then
     fermi_dirac_ct_diff = ratio5(0)
  elseif(nderiv.eq.1) then
     fermi_dirac_ct_diff = ratio5(1)*ratio1(1)
  elseif(nderiv.eq.2) then
     fermi_dirac_ct_diff =&
          ratio5(2)*ratio1(1)*ratio1(1) + ratio5(1)*ratio1(2)
  elseif(nderiv.eq.3) then
     fermi_dirac_ct_diff = ratio5(3)*ratio1(1)*ratio1(1)*ratio1(1) +&
          3._fp_kind*ratio5(2)*ratio1(1)*ratio1(2) + ratio5(1)*ratio1(3)
  elseif(nderiv.eq.4) then
     fermi_dirac_ct_diff =&
          ratio5(4)*ratio1(1)*ratio1(1)*ratio1(1)*ratio1(1) +&
          6._fp_kind*ratio5(3)*ratio1(1)*ratio1(1)*ratio1(2) +&
          4._fp_kind*ratio5(2)*ratio1(1)*ratio1(3) +&
          3._fp_kind*ratio5(2)*ratio1(2)*ratio1(2) + ratio5(1)*ratio1(4)
  endif
end function fermi_dirac_ct_diff
