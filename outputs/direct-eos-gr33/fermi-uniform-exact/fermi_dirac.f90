!*******************************************************************************
!       Copyright (C) 1989-2022 Alan W. Irwin
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

! subroutine to calculate *scaled* n_e, p_e, s_e, and u_e and the
! derivative of the *ln* of these quantities wrt fl, and tcl.  value
! and the fl and tcl derivatives are stored
! in the arrays rhostar(9), pstar(9), sstar(3), and ustar(3).
! The six extra derivatives of ln rhostar and pstar are in
! the order flfl, fltcl, tcltcl, flflfl, flfltcl, fltcltcl.
! fl is ln f and tcl is ln tc, where tc = kt/mc^2 equiv beta. (see Eggleton,
! Faulkner and Flannery paper and vdb notes).
! scaling:
! n_e (number density of electrons) = 8 pi/compton^3 rhostar = c_e*rhostar
!   ==> rhostar = sqrt(2) beta^3/2 [F(1/2 + beta F(3/2)] (C.G. 24.98)
! p_e (electron pressure) = 8 pi/compton^3 m c^2 pstar = cpe*pstar
!   ==> pstar = (2/3)sqrt(2) beta^5/2 [F(3/2 + (1/2)beta F(5/2)] (C.G. 24.99)
! s_e (entropy per unit volume of the free electrons) CG. 24.76b
!   = (u_e + p_e)/T - psi n_e k
!   = n_e k sstar
! note u_e (see later) is the internal energy per unit volume of the electrons.
! internally we calculate:
! qstar, where the convention for qstar in this programme is
! (1+f)/g times the qstar defined in the paper.
! Therefore, the following transformations apply using this convention:
! sstar = qstar/(rhostar*sqrt(1+f)) + 2 sqrt(1+f) - psi
!   = qstar/(rhostar*sqrt(1+f)) + ln((1+sqrt(1+f))^2/f)
! finally we calculate
! ustar = sstar + psi - pstar sqrt(1+f)/(g rhostar)
! this is related to u_e via
! u_e = n_e k T  ustar = n_e m c^2 beta ustar
! ==> ustar = [F(3/2 + beta F(5/2)]/[F(1/2 + beta F(3/2)] (C.G. 24.100)
!   = s_e T + psi n_e k T - pstar mc^2 n_e/rhostar
!   = s_e T + psi n_e k T - p_e (checks with previous).
! n.b. the ustar in this routine is different from the ustar in
! the vdb notes.
! more relation to CG:
! sstar = (u_e + p_e)/n_e k T - psi
!   = ustar + pstar/(rhostar beta) - psi
!   = ustar + (2/3) [F(3/2 + (1/2)beta F(5/2)]/[F(1/2 + beta F(3/2)] - psi
!   = [(5/3)F(3/2 + (4/3)beta F(5/2)]/[F(1/2 + beta F(3/2)] - psi
! morder_in = 3, 5, or 8  (or 13, 15, or 18) means 3rd, 5th, or 8th order
! fermi_dirac integral approximation following EFF fit (or modified
! version of EFF fit which reduces to Cody-Thacher approximation for
! low relativistic correction).
! morder_in = -3, -5, or -8 (or -13, -15, or -18) means use above
! approximations in non-relativistic limit.
! morder_in = -1 or 1 uses Cody-Thacher approximation directly.
! morder_in = 21 calculate Fermi-Dirac integrals with slow, but precise
! (~1.d-9 relative errors) numerical integration.
! morder_in = -21 is same as 21 in non-relativistic limit.
! morder_in = 23 means original 3rd order eff result
! morder_in = -23 is same as 23 in non-relativistic limit.

!> This fermi_dirac subroutine calculates scaled n_e, p_e, s_e, and
!> u_e (the number density, partial pressure, partial entropy and
!> partial energy of free electrons) as well as several different orders of
!> the (mixed) partial derivatives of the *ln* of these quantities wrt
!> fl, and tcl.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine fermi_dirac(verbosity, fl_in, tcl, rhostar, pstar, sstar, ustar, morder_in)

  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real
  use mod_flow_data, only: ln_overflow_limit, ln_underflow_limit

  ! Arguments
  integer, intent(in) :: verbosity, morder_in
  real(fp_kind), intent(in) :: fl_in, tcl
  real(fp_kind), intent(out) :: rhostar(:), pstar(:), sstar(:), ustar(:)

  ! Internal variables
  integer, parameter :: mderiv = 2
  integer, parameter :: mordermax = 8
  integer, parameter :: maxfd_direct = 9

  integer morder_old
  ! must be invalid value!
  data morder_old/0/

  integer mstar, morder, norder
  integer i, ideriv, nderiv, ioffset, moffset

  real(fp_kind) rhostarcon, pstarcon
  real(fp_kind) fl, f, wf, tc, g, vf, vg, uf, duf, duf2, ug, dug, dug2,&
       fdf, psi, qconst, pconst, sumadd, sumaddconst,&
       power, powerufm1, poweruf,&
       argm1, logarg, dlogarg

  ! n.b ccoeff in order of rhocoeff, pcoeff, qcoeff
  ! n.b. dimensions on ccoeff large enough for morder = 8 call to fermi_dirac_coeff.
  ! ccoeff not allocatable because it must be saved.
  real(fp_kind) ccoeff(3*(mordermax+1)*(mordermax+1))

  real(fp_kind), allocatable ::&
       acoeff(:, :, :),&
       fd(:),&
       fd_diff(:),&
       fd_direct(:,:),&
       sum(:)

  save morder_old, ccoeff

  mstar = size(sstar)
  ! Sanity check
  if(mstar+6.ne.size(rhostar).or.mstar+6.ne.size(pstar).or.mstar.ne.size(ustar))&
       error stop 'fermi_dirac: inconsistent sizes of rhostar, pstar, sstar, or ustar'

  ! As some protection against an fl iteration that has gone really
  ! bad, pin the fl value used in the calculations below to within a
  ! very wide dynamic range which nevertheless has sensible margins
  ! that make it unlikely the results below will generate any serious
  ! floating-point exceptions.
  ! Maintenance, 2021.  Make these limits all the same in
  ! master_exchange, exchange_gcpf, and fermi_dirac.

  fl = max(min(fl_in, ln_overflow_limit), ln_underflow_limit)

  f = exp(fl)
  ! d psi d fl
  wf = sqrt(1._fp_kind+f)
  if(f.gt.1._fp_kind) then
     psi = 2._fp_kind*wf + log((wf-1._fp_kind)/(wf+1._fp_kind))
  else
     psi = fl + 2._fp_kind*(wf - log(wf+1._fp_kind))
  endif
  ! argument of log = (1+wf)/(1-wf) = (1+wf)^2/f.
  ! subtract 1 from that argument in a way that does not generate
  ! large significance loss.
  argm1 = 2._fp_kind*(1._fp_kind+wf)/f
  if(argm1.gt.1.e-3_fp_kind) then
     logarg = log(1._fp_kind+argm1)
  else
     ! alternating series so relative error is less than first
     ! missing term which is argm1^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
     logarg = argm1*&
          (1._fp_kind      - argm1*&
          (1._fp_kind/2._fp_kind - argm1*&
          (1._fp_kind/3._fp_kind - argm1*&
          (1._fp_kind/4._fp_kind - argm1*&
          (1._fp_kind/5._fp_kind)))))
  endif
  ! Calculate (straightforward) derivative of log ((1+wf)/(1-wf)) rather
  ! than using (more complex) argm1 form.
  ! no large significance loss at any f.
  dlogarg = -1._fp_kind/wf
  ! a.k.a. beta = kt/mc^2
  tc = exp(tcl)
  morder = abs(morder_in)
  if(morder.eq.1) then
     ! Start of logic block taken independent of all else to use
     ! Cody-Thacher approximations to Fermi-Dirac integral in *non*
     ! relativistic limit.
     allocate(fd(0:4))

     if(if_taint_allocated_real) then
        call taint_allocated_real(fd)
     endif

     do ideriv = 0, 4
        ! F(3/2) and 4 derivatives
        fd(ideriv) = fermi_dirac_ct(psi,ideriv,0)
     enddo
     ! N.B. F 1/2 = (2/3) d F 3/2 d psi so that is where the factor
     ! of (2/3) comes from.
     rhostarcon = sqrt(2._fp_kind) * exp(1.5_fp_kind*tcl) * (2._fp_kind/3._fp_kind)
     rhostar(1) = rhostarcon*fd(1)
     ! ln derivatives in order function, f, t, ff, ft, tt, fff, fft, ftt
     rhostar(2) = fd(2)/fd(1)*wf
     rhostar(3) = 1.5_fp_kind
     rhostar(4) = (-fd(2)*fd(2)*(1._fp_kind+f)/fd(1) +&
          fd(3)*(1._fp_kind+f) + 0.5_fp_kind*fd(2)*f/wf)/fd(1)
     rhostar(5) = 0._fp_kind
     rhostar(6) = 0._fp_kind
     rhostar(7) = -rhostar(4)*(fd(2)/fd(1))*wf + (&
          -fd(2)*fd(2)*f/fd(1) -&
          2._fp_kind*wf*fd(2)*fd(3)*(1._fp_kind+f)/fd(1) +&
          wf*fd(2)*(fd(2)/fd(1))*(fd(2)/fd(1))*(1._fp_kind+f) +&
          fd(4)*(1._fp_kind+f)*wf + 1.5_fp_kind*fd(3)*f +&
          0.5_fp_kind*fd(2)*(f/wf)*(1._fp_kind - 0.5_fp_kind*f/(1._fp_kind+f)))/fd(1)
     rhostar(8) = 0._fp_kind
     rhostar(9) = 0._fp_kind
     pstarcon = sqrt(2._fp_kind) * exp(2.5_fp_kind*tcl) * (2._fp_kind/3._fp_kind)
     pstar(1) = pstarcon*fd(0)
     pstar(2) = fd(1)/fd(0)*wf
     pstar(3) = 2.5_fp_kind
     pstar(4) = (-fd(1)*fd(1)*(1._fp_kind+f)/fd(0) +&
          fd(2)*(1._fp_kind+f) + 0.5_fp_kind*fd(1)*f/wf)/fd(0)
     pstar(5) = 0._fp_kind
     pstar(6) = 0._fp_kind
     pstar(7) = -pstar(4)*(fd(1)/fd(0))*wf + (&
          -fd(1)*fd(1)*f/fd(0) -&
          2._fp_kind*wf*fd(1)*fd(2)*(1._fp_kind+f)/fd(0) +&
          wf*fd(1)*(fd(1)/fd(0))*(fd(1)/fd(0))*(1._fp_kind+f) +&
          fd(3)*(1._fp_kind+f)*wf + 1.5_fp_kind*fd(2)*f +&
          0.5_fp_kind*fd(1)*(f/wf)*(1._fp_kind - 0.5_fp_kind*f/(1._fp_kind+f)))/fd(0)
     pstar(8) = 0._fp_kind
     pstar(9) = 0._fp_kind
     ustar(1) = 1.5_fp_kind*fd(0)/fd(1)
     ustar(2) = wf*(fd(1)/fd(0) - fd(2)/fd(1))
     ustar(3) = 0._fp_kind
     ! N.B. must be same large psi limit as in fermi_dirac_ct_diff.
     if(psi.le.4._fp_kind) then
        sstar(1) = 2.5_fp_kind*fd(0)/fd(1) - psi
        sstar(2) = wf*(1.5_fp_kind - 2.5_fp_kind*(fd(0)/fd(1))*(fd(2)/fd(1)))/&
             sstar(1)
     else
        ! calculate 2.5 F_3/2 - psi F'_3/2 and derivative wrt psi
        ! without incurring large signficance loss.
        allocate(fd_diff(0:1))

        if(if_taint_allocated_real) then
           call taint_allocated_real(fd_diff)
        endif

        fd_diff(0) = fermi_dirac_ct_diff(psi, 0)
        fd_diff(1) = fermi_dirac_ct_diff(psi, 1)
        sstar(1) = fd_diff(0)/fd(1)
        sstar(2) = wf*(fd(1)*fd_diff(1) - fd(2)*fd_diff(0))/(fd(1)*fd(1)*sstar(1))
     endif
     sstar(3) = 0._fp_kind
     ! End of logic block taken independent of all else to use
     ! Cody-Thacher approximations to Fermi-Dirac integral in *non*
     ! relativistic limit.
     return
  elseif(abs(morder_in).eq.21) then
     ! Start of logic block for direct, near-exact, but slow
     ! calculation of Fermi-Dirac integrals using numerical
     ! integration.
     if(morder_in.eq.-21) then
        allocate(fd_direct(7,2))

        if(if_taint_allocated_real) then
           call taint_allocated_real(fd_direct)
        endif

        call fermi_dirac_direct(0.5_fp_kind, psi, 0._fp_kind, fd_direct(:,1))
        call fermi_dirac_direct(1.5_fp_kind, psi, 0._fp_kind, fd_direct(:,2))
        ! derivatives in order function, f, t, ff, ft, tt, fff, fft, ftt
        rhostar(1) = fd_direct(1,1)
        rhostar(2) = fd_direct(2,1)
        rhostar(4) = fd_direct(4,1)
        rhostar(7) = fd_direct(7,1)
        ! Zero beta derivatives.
        rhostar(3) = 0._fp_kind
        rhostar(5) = 0._fp_kind
        rhostar(6) = 0._fp_kind
        rhostar(8) = 0._fp_kind
        rhostar(9) = 0._fp_kind
        pstar(1) = fd_direct(1,2)
        pstar(2) = fd_direct(2,2)
        pstar(4) = fd_direct(4,2)
        pstar(7) = fd_direct(7,2)
        ! Zero beta derivatives.
        pstar(3) = 0._fp_kind
        pstar(5) = 0._fp_kind
        pstar(6) = 0._fp_kind
        pstar(8) = 0._fp_kind
        pstar(9) = 0._fp_kind
        ustar(1) = fd_direct(1,2)
        ustar(2) = fd_direct(2,2)
        ! Zero beta derivatives.
        ustar(3) = 0._fp_kind
     else !if(morder_in.eq.-21) then
        allocate(fd_direct(maxfd_direct,3))

        if(if_taint_allocated_real) then
           call taint_allocated_real(fd_direct)
        endif

        call fermi_dirac_direct(0.5_fp_kind, psi, tc, fd_direct(:,1))
        call fermi_dirac_direct(1.5_fp_kind, psi, tc, fd_direct(:,2))
        call fermi_dirac_direct(2.5_fp_kind, psi, tc, fd_direct(:,3))

        ! derivatives in order function, f, t, ff, ft, tt, fff, fft, ftt
        rhostar(1:maxfd_direct) = fd_direct(1:maxfd_direct,1) + tc*fd_direct(1:maxfd_direct,2)
        pstar(1:maxfd_direct) = fd_direct(1:maxfd_direct,2) + 0.5_fp_kind*tc*fd_direct(1:maxfd_direct,3)
        ustar(1:mstar) = fd_direct(1:mstar,2) + tc*fd_direct(1:mstar,3)

        ! Transform to total beta=tc derivative of above rhostar, pstar.
        rhostar(3) = rhostar(3) + fd_direct(1,2)
        rhostar(5) = rhostar(5) + fd_direct(2,2)
        rhostar(6) = rhostar(6) + 2._fp_kind*fd_direct(3,2)
        rhostar(8) = rhostar(8) + fd_direct(4,2)
        rhostar(9) = rhostar(9) + 2._fp_kind*fd_direct(5,2)
        pstar(3) = pstar(3) + 0.5_fp_kind*fd_direct(1,3)
        pstar(5) = pstar(5) + 0.5_fp_kind*fd_direct(2,3)
        pstar(6) = pstar(6) + fd_direct(3,3)
        pstar(8) = pstar(8) + 0.5_fp_kind*fd_direct(4,3)
        pstar(9) = pstar(9) + fd_direct(5,3)
        ustar(3) = ustar(3) + fd_direct(1,3)
        ! Transform to tcl = ln beta derivatives from tc=beta derivatives
        rhostar(9) = rhostar(9)*tc*tc + rhostar(5)*tc
        rhostar(8) = rhostar(8)*tc
        rhostar(6) = rhostar(6)*tc*tc + rhostar(3)*tc
        rhostar(5) = rhostar(5)*tc
        rhostar(3) = rhostar(3)*tc
        pstar(9) = pstar(9)*tc*tc + pstar(5)*tc
        pstar(8) = pstar(8)*tc
        pstar(6) = pstar(6)*tc*tc + pstar(3)*tc
        pstar(5) = pstar(5)*tc
        pstar(3) = pstar(3)*tc
        ustar(3) = ustar(3)*tc
     endif !f(morder_in.eq.-21) then
     ! Transform to fl derivatives from psi derivatives.  wf = dpsi/dfl
     rhostar(9) = rhostar(9)*wf
     rhostar(8) = rhostar(8)*wf*wf + rhostar(5)*(0.5_fp_kind*f/wf)
     rhostar(7) = rhostar(7)*wf*wf*wf + rhostar(4)*(1.5_fp_kind*f) +&
          rhostar(2)*0.5_fp_kind*f/wf*(1._fp_kind - 0.5_fp_kind*(f/wf)/wf)
     rhostar(5) = rhostar(5)*wf
     rhostar(4) = rhostar(4)*wf*wf + rhostar(2)*(0.5_fp_kind*f/wf)
     rhostar(2) = rhostar(2)*wf
     pstar(9) = pstar(9)*wf
     pstar(8) = pstar(8)*wf*wf + pstar(5)*(0.5_fp_kind*f/wf)
     pstar(7) = pstar(7)*wf*wf*wf + pstar(4)*(1.5_fp_kind*f) + pstar(2)*0.5_fp_kind*f/wf*(1._fp_kind - 0.5_fp_kind*(f/wf)/wf)
     pstar(5) = pstar(5)*wf
     pstar(4) = pstar(4)*wf*wf + pstar(2)*(0.5_fp_kind*f/wf)
     pstar(2) = pstar(2)*wf
     ustar(2) = ustar(2)*wf
     ! Transform to derivatives of ln of function, but do not transform
     ! function itself (which is the convention used in this routine).
     rhostar(9) = rhostar(9)/rhostar(1) -&
          (rhostar(6)/rhostar(1))*(rhostar(2)/rhostar(1)) -&
          2._fp_kind*(rhostar(5)/rhostar(1))*(rhostar(3)/rhostar(1)) +&
          2._fp_kind*(rhostar(3)/rhostar(1))*(rhostar(3)/rhostar(1))*&
          (rhostar(2)/rhostar(1))
     rhostar(8) = rhostar(8)/rhostar(1) -&
          (rhostar(4)/rhostar(1))*(rhostar(3)/rhostar(1)) -&
          2._fp_kind*(rhostar(5)/rhostar(1))*(rhostar(2)/rhostar(1)) +&
          2._fp_kind*(rhostar(3)/rhostar(1))*(rhostar(2)/rhostar(1))*&
          (rhostar(2)/rhostar(1))
     rhostar(7) = rhostar(7)/rhostar(1) -&
          3._fp_kind*(rhostar(4)/rhostar(1))*(rhostar(2)/rhostar(1)) +&
          2._fp_kind*(rhostar(2)/rhostar(1))*(rhostar(2)/rhostar(1))*&
          (rhostar(2)/rhostar(1))
     rhostar(6) = rhostar(6)/rhostar(1) -&
          (rhostar(3)/rhostar(1))*(rhostar(3)/rhostar(1))
     rhostar(5) = rhostar(5)/rhostar(1) -&
          (rhostar(2)/rhostar(1))*(rhostar(3)/rhostar(1))
     rhostar(4) = rhostar(4)/rhostar(1) -&
          (rhostar(2)/rhostar(1))*(rhostar(2)/rhostar(1))
     rhostar(3) = rhostar(3)/rhostar(1)
     rhostar(2) = rhostar(2)/rhostar(1)
     pstar(9) = pstar(9)/pstar(1) -&
          (pstar(6)/pstar(1))*(pstar(2)/pstar(1)) -&
          2._fp_kind*(pstar(5)/pstar(1))*(pstar(3)/pstar(1)) +&
          2._fp_kind*(pstar(3)/pstar(1))*(pstar(3)/pstar(1))*&
          (pstar(2)/pstar(1))
     pstar(8) = pstar(8)/pstar(1) -&
          (pstar(4)/pstar(1))*(pstar(3)/pstar(1)) -&
          2._fp_kind*(pstar(5)/pstar(1))*(pstar(2)/pstar(1)) +&
          2._fp_kind*(pstar(3)/pstar(1))*(pstar(2)/pstar(1))*&
          (pstar(2)/pstar(1))
     pstar(7) = pstar(7)/pstar(1) -&
          3._fp_kind*(pstar(4)/pstar(1))*(pstar(2)/pstar(1)) +&
          2._fp_kind*(pstar(2)/pstar(1))*(pstar(2)/pstar(1))*&
          (pstar(2)/pstar(1))
     pstar(6) = pstar(6)/pstar(1) -&
          (pstar(3)/pstar(1))*(pstar(3)/pstar(1))
     pstar(5) = pstar(5)/pstar(1) -&
          (pstar(2)/pstar(1))*(pstar(3)/pstar(1))
     pstar(4) = pstar(4)/pstar(1) -&
          (pstar(2)/pstar(1))*(pstar(2)/pstar(1))
     pstar(3) = pstar(3)/pstar(1)
     pstar(2) = pstar(2)/pstar(1)
     ustar(3) = ustar(3)/ustar(1)
     ustar(2) = ustar(2)/ustar(1)
     ! N.B. at this stage we have the following:
     ! rhostar = F 1/2 + beta F 3/2
     ! pstar = F 3/2 + 0.5 beta F 5/2
     ! ustar = F 3/2 + beta F 5/2
     !   plus partial derivatives of the ln of these quantities.
     !
     ! Transform to final forms:
     ! ustar to ratio of present ustar and rhostar.
     ustar(1) = ustar(1)/rhostar(1)
     ustar(2) = ustar(2) - rhostar(2)
     ustar(3) = ustar(3) - rhostar(3)
     ! Multiply by factor to convert to rhostar (or add 1.5 tcl)
     ! to ln rhostar for partial-derivative purposes).
     ! N.B. no factor of 2/3 in this version since we are dealing with
     ! FD integrals rather than their derivatives as above for the
     ! non-releativistic Cody-Thacher case.
     rhostarcon = sqrt(2._fp_kind) * exp(1.5_fp_kind*tcl)
     rhostar(1) = rhostarcon*rhostar(1)
     rhostar(3) = rhostar(3) + 1.5_fp_kind
     ! Multiply by factor to convert to pstar (or add 2.5 tcl)
     ! to ln pstar for partial-derivative purposes).
     pstarcon = sqrt(2._fp_kind) * exp(2.5_fp_kind*tcl) * (2._fp_kind/3._fp_kind)
     pstar(1) = pstarcon*pstar(1)
     pstar(3) = pstar(3) + 2.5_fp_kind
     ! finally we calculate
     ! sstar = ustar - psi + pstar sqrt(1+f)/(g rhostar)
     ! significance loss for large degeneracy to be fixed later
     ! if ever. (only loses 3 figures at psi = 100)
     sstar(1) = ustar(1) - psi + pstar(1)/(tc*rhostar(1))
     sstar(2) = (ustar(1)*ustar(2) - wf + (pstar(1)/(tc*rhostar(1)))*(pstar(2)-rhostar(2)))/sstar(1)
     sstar(3) = (ustar(1)*ustar(3) + (pstar(1)/(tc*rhostar(1)))*(pstar(3)-rhostar(3)-1._fp_kind))/sstar(1)
     ! End of logic block for direct, near-exact, but slow
     ! calculation of Fermi-Dirac integrals using numerical
     ! integration.
     return
  endif
  ! At this stage have dealt with mord_in = -1, 1, -21, or 21 result.
  ! Now deal with remaining possibilities (which include special case of morder_in = -23 or 23).
  ! N.B. morder = abs(morder_in) (see above).
  if(mstar.ne.mderiv+1) error stop 'invalid mstar input to fermi_dirac'
  if(abs(morder_in).gt.20) then
     morder = morder - 20
  elseif(abs(morder_in).gt.10) then
     morder = morder - 10
  endif

  if(abs(morder_in).ne.morder_old) then
     morder_old = abs(morder_in)
     allocate(acoeff(0:morder, 0:morder, 3))

     if(if_taint_allocated_real) then
        call taint_allocated_real(acoeff)
     endif

     ! fill coefficients with thermodynamically consistent set.
     if(abs(morder_in).gt.20) then
        call fermi_dirac_original_coeff(verbosity, acoeff)
     elseif(abs(morder_in).gt.10) then
        call fermi_dirac_minusnr_coeff(acoeff)
     else
        call fermi_dirac_coeff(acoeff)
     endif
     ccoeff(1:3*(morder+1)*(morder+1)) = reshape(acoeff, shape(ccoeff(1:3*(morder+1)*(morder+1))))
  endif
  if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then
     ! Collect zero-temperature Cody-Thacher results that were subtracted
     ! in the fit used to generate the fermi_dirac_minusnr_coeff results.
     allocate(fd(0:4))

     if(if_taint_allocated_real) then
        call taint_allocated_real(fd)
     endif

     do ideriv = 0, 4
        ! F(3/2) and 4 derivatives
        fd(ideriv) = fermi_dirac_ct(psi,ideriv,0)
     enddo
  endif
  ! f = exp(fl)
  ! wf = sqrt(1.d0+f)
  vf = 1._fp_kind/(1._fp_kind+f)
  ! partial of ln (1+f) wrt ln f
  uf = f*vf
  ! partial uf wrt ln f.  N.B. 1-uf = 1/(1+f)
  duf = uf*vf
  duf2 = uf*(1._fp_kind-f)*vf*vf
  if(morder_in.gt.0) then
     g = tc*wf
     vg = 1._fp_kind/(1._fp_kind+g)
     fdf = g*g+g
     fdf = uf*fdf*sqrt(fdf)*(vf*vg)**morder
     ! partial of ln (1+g) wrt ln g
     ug = g*vg
     ! partial ug wrt ln g
     dug = ug*vg
     dug2 = ug*(1._fp_kind-g)*vg*vg
  else
     ! negative morder_in flags NR limit.
     g = 0._fp_kind
     vg = 1._fp_kind
     ! d ln(1+g)/dln g in limit where 1+g replaced by 1.
     ug = 0._fp_kind
     fdf = uf*(tc*wf)*sqrt(tc*wf)*vf**morder
     dug = 0._fp_kind
     dug2 = 0._fp_kind
  endif
  if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then
     ! When non-relativistic limit is treated separately, lowest order
     ! (in g) coefficients are zero.  Ignore these coefficients in
     ! sum and (later) multiply series results by an additional factor
     ! of g.
     norder = morder - 1
     moffset = 1 + (morder+1)
  else
     norder = morder
     moffset = 1
  endif
  do i = 1,3
     ! i == 1 corresponds to calculating rhostar
     ! i == 2 corresponds to calculating pstar
     ! i == 3 corresponds to calculating qstar which is stored in sstar for now
     if(i.le.2) then
        nderiv = mderiv + 7
     else
        nderiv = mderiv
     endif
     allocate(sum(nderiv+1))

     if(if_taint_allocated_real) then
        call taint_allocated_real(sum)
     endif

     ! access rhostar, pstar, or qstar coefficients.
     ioffset = (i-1)*(morder+1)*(morder+1) + moffset
     ! Note for 10< abs(morder_in)<= 20, the zero-order g coefficients are
     ! zero so that norder and moffset are manipulated so that
     ! effsum_calc just calculates polynomial multiplier for g.
     call effsum_calc(f, g, reshape(ccoeff(ioffset:ioffset+(morder+1)*(norder+1)-1),[morder+1,norder+1]), sum)

     if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then
        ! Factor to convert from fd to non-relativistic pstar and rhostar.
        pstarcon = sqrt(2._fp_kind)*(2._fp_kind/3._fp_kind)

        ! Multiply sum by g and adjust sum derivatives accordingly since
        ! for this case effsum_calc (with the appropriate norder and
        ! moffset) just calculates polynomial multiplier of g.
        if(i.le.2) then
           ! ggg
           sum(10) = g*(sum(1) + 3._fp_kind*sum(3) + 3._fp_kind*sum(6) + sum(10))
           ! fgg
           sum(9) = g*(sum(2) + 2._fp_kind*sum(5) + sum(9))
           ! ffg
           sum(8) = g*(sum(4) + sum(8))
           ! fff
           sum(7) = g*sum(7)
           ! gg
           sum(6) = g*(sum(1) + 2._fp_kind*sum(3) + sum(6))
           ! fg
           sum(5) = g*(sum(2) + sum(5))
           ! ff
           sum(4) = g*sum(4)
        endif
        ! g
        sum(3) = g*(sum(1) + sum(3))
        ! f
        sum(2) = g*sum(2)
        sum(1) = g*sum(1)
        ! Add in zero-temperature Cody-Thacher results that were subtracted
        ! in the fit used to generate the fermi_dirac_minusnr_coeff results.
        ! N.B. results will be multiplied afterward by fdf factor where
        ! fdf = f/(1+f) (g*(1+g))^3/2 ((1+f)*(1+g))^(-morder)
        if(i.eq.1) then
           ! Start of logic block for i.eq.1, i.e., rhostar calculation
           power = real((morder),fp_kind)+0.25_fp_kind
           powerufm1 = power*uf - 1._fp_kind
           sumadd = pstarcon*(1._fp_kind+g)*fd(1)*((1._fp_kind+f)**power)/f
           sum(1) = sum(1) + sumadd
           sum(2) = sum(2) + sumadd*(fd(2)*wf/fd(1) + powerufm1)
           sum(3) = sum(3) + ug*sumadd
           sum(4) = sum(4) + sumadd*(&
                (fd(3)*wf + 0.5_fp_kind*fd(2)*uf)*wf/fd(1) + power*duf +&
                powerufm1*(2._fp_kind*fd(2)*wf/fd(1) + powerufm1))
           sum(5) = sum(5) + ug*sumadd*(fd(2)*wf/fd(1) + powerufm1)
           sum(6) = sum(6) + ug*sumadd
           ! sum(7) = sum(7) + sumadd*(&
           !   ((fd(4)*wf + fd(3)*uf)*wf + 0.5d0*fd(2)*duf)*wf/fd(1) + &
           !   (fd(3)*wf + 0.5d0*fd(2)*uf)*wf*&
           !   (0.5d0*uf - fd(2)*wf/fd(1))/fd(1) + power*duf2 +&
           !   power*duf*(2.d0*fd(2)*wf/fd(1) + 2.d0*powerufm1) +&
           !   powerufm1*2.d0*(&
           !   fd(3)*wf + 0.5d0*fd(2)*uf -&
           !   fd(2)*fd(2)*wf/fd(1))*wf/fd(1) + (&
           !   (fd(3)*wf + 0.5d0*fd(2)*uf)*wf/fd(1) + power*duf +&
           !   powerufm1*(2.d0*fd(2)*wf/fd(1) + powerufm1))*(&
           !   fd(2)*wf/fd(1) + powerufm1))
           ! commented out expression above is straight derivative which has
           ! been tested to be correct.  Expression below has terms
           ! consolidated from above.  This is also tested to be correct,
           ! but it a lot less understandable than expression above.
           sum(7) = sum(7) + sumadd*(&
                (&
                (fd(4)*wf + 1.5_fp_kind*fd(3)*uf)*wf + 0.5_fp_kind*fd(2)*duf&
                )*wf/fd(1) +&
                (0.5_fp_kind*fd(2)*uf)*wf*(0.5_fp_kind*uf)/fd(1) +&
                power*duf2 +&
                power*duf*(3._fp_kind*fd(2)*wf/fd(1)) +&
                powerufm1*(&
                3._fp_kind*(fd(3)*wf + 0.5_fp_kind*fd(2)*uf)*wf/fd(1) +&
                3._fp_kind*power*duf + powerufm1*(&
                3._fp_kind*fd(2)*wf/fd(1) + powerufm1))&
                )
           sum(8) = sum(8) + ug*sumadd*(&
                (fd(3)*wf + 0.5_fp_kind*fd(2)*uf)*wf/fd(1) +&
                power*duf +&
                powerufm1*(2._fp_kind*fd(2)*wf/fd(1) +&
                powerufm1))
           sum(9) = sum(9) + ug*sumadd*(fd(2)*wf/fd(1) + powerufm1)
           sum(10) = sum(10) + ug*sumadd
           ! add in second term
           power = real((morder),fp_kind)-1.25_fp_kind
           poweruf = power*uf
           sumadd = 0.5_fp_kind*(2.5_fp_kind - real((morder),fp_kind))*pstarcon*&
                g*fd(0)*((1._fp_kind+f)**power)
           sum(1) = sum(1) + sumadd
           sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + poweruf)
           sum(3) = sum(3) + sumadd
           sum(4) = sum(4) + sumadd*(&
                (fd(2)*wf + 0.5_fp_kind*fd(1)*uf)*wf/fd(0) + power*duf +&
                poweruf*(2._fp_kind*fd(1)*wf/fd(0) + poweruf))
           sum(5) = sum(5) + sumadd*(fd(1)*wf/fd(0) + poweruf)
           sum(6) = sum(6) + sumadd
           ! sum(7) = sum(7) + sumadd*(&
           !   ((fd(3)*wf + fd(2)*uf)*wf + 0.5d0*fd(1)*duf)*wf/fd(0) +&
           !   (fd(2)*wf + 0.5d0*fd(1)*uf)*wf*&
           !   (0.5d0*uf - fd(1)*wf/fd(0))/fd(0) + power*duf2 +&
           !   power*duf*(2.d0*fd(1)*wf/fd(0) + 2.d0*poweruf) +&
           !   poweruf*2.d0*(&
           !   fd(2)*wf + 0.5d0*fd(1)*uf -&
           !   fd(1)*fd(1)*wf/fd(0))*wf/fd(0) + (&
           !   (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +&
           !   poweruf*(2.d0*fd(1)*wf/fd(0) + poweruf))*(&
           !   fd(1)*wf/fd(0) + poweruf))
           ! commented out expression above is straight derivative which has
           ! been tested to be correct.  Expression below has terms
           ! consolidated from above.  This is also tested to be correct,
           ! but it a lot less understandable than expression above.
           sum(7) = sum(7) + sumadd*(&
                (&
                (fd(3)*wf + 1.5_fp_kind*fd(2)*uf)*wf + 0.5_fp_kind*fd(1)*duf&
                )*wf/fd(0) +&
                (0.5_fp_kind*fd(1)*uf)*wf*(0.5_fp_kind*uf)/fd(0) +&
                power*duf2 +&
                power*duf*(3._fp_kind*fd(1)*wf/fd(0)) +&
                poweruf*(&
                3._fp_kind*(fd(2)*wf + 0.5_fp_kind*fd(1)*uf)*wf/fd(0) +&
                3._fp_kind*power*duf + poweruf*(&
                3._fp_kind*fd(1)*wf/fd(0) + poweruf))&
                )
           sum(8) = sum(8) + sumadd*(&
                (fd(2)*wf + 0.5_fp_kind*fd(1)*uf)*wf/fd(0) + power*duf +&
                poweruf*(2._fp_kind*fd(1)*wf/fd(0) + poweruf))
           sum(9) = sum(9) + sumadd*(fd(1)*wf/fd(0) + poweruf)
           sum(10) = sum(10) + sumadd
           ! End of logic block for i.eq.1
        elseif(i.eq.2) then
           ! Start of logic block for i.eq.2, i.e., pstar calculation
           power = real((morder),fp_kind)-0.25_fp_kind
           powerufm1 = power*uf - 1._fp_kind
           sumadd = pstarcon*(1._fp_kind+g)*fd(0)*((1._fp_kind+f)**power)/f
           sum(1) = sum(1) + sumadd
           sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + powerufm1)
           sum(3) = sum(3) + ug*sumadd
           sum(4) = sum(4) + sumadd*(&
                (fd(2)*wf + 0.5_fp_kind*fd(1)*uf)*wf/fd(0) + power*duf +&
                powerufm1*(2._fp_kind*fd(1)*wf/fd(0) + powerufm1))
           sum(5) = sum(5) + ug*sumadd*(fd(1)*wf/fd(0) + powerufm1)
           sum(6) = sum(6) + ug*sumadd
           ! sum(7) = sum(7) + sumadd*(&
           !   ((fd(3)*wf + fd(2)*uf)*wf + 0.5d0*fd(1)*duf)*wf/fd(0) +&
           !   (fd(2)*wf + 0.5d0*fd(1)*uf)*wf*&
           !   (0.5d0*uf - fd(1)*wf/fd(0))/fd(0) + power*duf2 +&
           !   power*duf*(2.d0*fd(1)*wf/fd(0) + 2.d0*powerufm1) +&
           !   powerufm1*2.d0*(&
           !   fd(2)*wf + 0.5d0*fd(1)*uf -&
           !   fd(1)*fd(1)*wf/fd(0))*wf/fd(0) + (&
           !   (fd(2)*wf + 0.5d0*fd(1)*uf)*wf/fd(0) + power*duf +&
           !   powerufm1*(2.d0*fd(1)*wf/fd(0) + powerufm1))*(&
           !   fd(1)*wf/fd(0) + powerufm1))
           ! commented out expression above is straight derivative which has
           ! been tested to be correct.  Expression below has terms
           ! consolidated from above.  This is also tested to be correct,
           ! but it a lot less understandable than expression above.
           sum(7) = sum(7) + sumadd*(&
                (&
                (fd(3)*wf + 1.5_fp_kind*fd(2)*uf)*wf + 0.5_fp_kind*fd(1)*duf&
                )*wf/fd(0) +&
                (0.5_fp_kind*fd(1)*uf)*wf*(0.5_fp_kind*uf)/fd(0) +&
                power*duf2 +&
                power*duf*(3._fp_kind*fd(1)*wf/fd(0)) +&
                powerufm1*(&
                3._fp_kind*(fd(2)*wf + 0.5_fp_kind*fd(1)*uf)*wf/fd(0) +&
                3._fp_kind*power*duf + powerufm1*(&
                3._fp_kind*fd(1)*wf/fd(0) + powerufm1))&
                )
           sum(8) = sum(8) + ug*sumadd*(&
                (fd(2)*wf + 0.5_fp_kind*fd(1)*uf)*wf/fd(0) +&
                power*duf +&
                powerufm1*(2._fp_kind*fd(1)*wf/fd(0) +&
                powerufm1))
           sum(9) = sum(9) + ug*sumadd*(fd(1)*wf/fd(0) + powerufm1)
           sum(10) = sum(10) + ug*sumadd
           ! End of logic block for i.eq.2
        elseif(i.eq.3) then
           ! Start of logic block for i.eq.3, i.e., qstar calculation (temporarily stored in ustar)
           ! N.B. must be same large psi limit as in fermi_dirac_ct_diff.
           if(psi.le.4._fp_kind) then
              ! add in first term.
              power = real((morder),fp_kind)+0.75_fp_kind
              powerufm1 = power*uf - 1._fp_kind
              sumadd = 2.5_fp_kind*pstarcon*(1._fp_kind+g)*fd(0)*((1._fp_kind+f)**power)/f
              sum(1) = sum(1) + sumadd
              sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + powerufm1)
              sum(3) = sum(3) + ug*sumadd
              ! add in second term.
              power = real((morder),fp_kind)+1.25_fp_kind
              powerufm1 = power*uf - 1._fp_kind
              sumadd = -2._fp_kind*pstarcon*(1._fp_kind+g)*fd(1)*((1._fp_kind+f)**power)/f
              sum(1) = sum(1) + sumadd
              sum(2) = sum(2) + sumadd*(fd(2)*wf/fd(1) + powerufm1)
              sum(3) = sum(3) + ug*sumadd
           else
              ! add in combination of first and second terms with analytical
              ! subtraction to reduce significance loss.
              ! calculate 2.5 F_3/2 - psi F'_3/2 and derivative wrt psi
              ! without incurring large significance loss.
              allocate(fd_diff(0:1))

              if(if_taint_allocated_real) then
                 call taint_allocated_real(fd_diff)
              endif

              fd_diff(0) = fermi_dirac_ct_diff(psi, 0)
              fd_diff(1) = fermi_dirac_ct_diff(psi, 1)
              ! psi = 2 wf + log((wf-1.d0)/(wf+1.d0)) = 2 wf - logarg
              ! Therefore multiplier is
              ! 2.5 F_3/2 - 2 wf F'_3/2 = fd_diff(0) + fd(0)*(psi - 2wf)
              ! = fd_diff(0) + fd(1)*(psi - 2 wf)
              ! = fd_diff(0) - fd(1)*logarg
              ! N.B. fd_diff(0) goes as pi^2/2 psi^{1/2} for large psi while
              ! fd(1)*logarg goes as psi^{3/2} 2/psi = 2 psi^{1/2}.
              ! Therefore, no large significance loss from subtraction of
              ! the two terms.
              power = real((morder),fp_kind)+0.75_fp_kind
              powerufm1 = power*uf - 1._fp_kind
              sumaddconst = pstarcon*(1._fp_kind+g)*((1._fp_kind+f)**power)/f
              sumadd = sumaddconst*(fd_diff(0) - fd(1)*logarg)
              sum(1) = sum(1) + sumadd
              sum(2) = sum(2) + sumadd*powerufm1 + sumaddconst*(wf*(fd_diff(1) - fd(2)*logarg) - fd(1)*dlogarg)
              sum(3) = sum(3) + ug*sumadd
           endif
           ! add in third term.
           power = real((morder),fp_kind)-0.25_fp_kind
           powerufm1 = power*uf - 1._fp_kind
           sumadd = (2.5_fp_kind - real((morder),fp_kind))*pstarcon*g*fd(0)*((1._fp_kind+f)**power)/f
           sum(1) = sum(1) + sumadd
           sum(2) = sum(2) + sumadd*(fd(1)*wf/fd(0) + powerufm1)
           sum(3) = sum(3) + sumadd
           ! End of logic block for i.eq.3
        endif
     endif !if((.not.abs(morder_in).gt.20).and.abs(morder_in).gt.10) then

     ! Transform to derivatives of ln of function, but do not transform
     ! function itself (which is the convention used in this routine).
     if(i.le.2) then
        sum(10) = sum(10)/sum(1) -&
             3._fp_kind*(sum(6)/sum(1))*(sum(3)/sum(1)) +&
             2._fp_kind*(sum(3)/sum(1))*(sum(3)/sum(1))*&
             (sum(3)/sum(1))
        sum(9) = sum(9)/sum(1) -&
             (sum(6)/sum(1))*(sum(2)/sum(1)) -&
             2._fp_kind*(sum(5)/sum(1))*(sum(3)/sum(1)) +&
             2._fp_kind*(sum(3)/sum(1))*(sum(3)/sum(1))*&
             (sum(2)/sum(1))
        sum(8) = sum(8)/sum(1) -&
             (sum(4)/sum(1))*(sum(3)/sum(1)) -&
             2._fp_kind*(sum(5)/sum(1))*(sum(2)/sum(1)) +&
             2._fp_kind*(sum(3)/sum(1))*(sum(2)/sum(1))*&
             (sum(2)/sum(1))
        sum(7) = sum(7)/sum(1) -&
             3._fp_kind*(sum(4)/sum(1))*(sum(2)/sum(1)) +&
             2._fp_kind*(sum(2)/sum(1))*(sum(2)/sum(1))*&
             (sum(2)/sum(1))
        sum(6) = sum(6)/sum(1) -&
             (sum(3)/sum(1))*(sum(3)/sum(1))
        sum(5) = sum(5)/sum(1) -&
             (sum(2)/sum(1))*(sum(3)/sum(1))
        sum(4) = sum(4)/sum(1) -&
             (sum(2)/sum(1))*(sum(2)/sum(1))
     endif
     sum(3) = sum(3)/sum(1)
     sum(2) = sum(2)/sum(1)
     ! Multiply sum by fdf where
     ! fdf = f/(1+f) (g*(1+g))^3/2 ((1+f)*(1+g))^(-morder)
     if(i.le.2) then
        sum(10) = sum(10) + (1.5_fp_kind-real((morder),fp_kind))*dug2
        sum(7) = sum(7) - real((morder+1),fp_kind)*duf2
        sum(6) = sum(6) + (1.5_fp_kind-real((morder),fp_kind))*dug
        sum(4) = sum(4) - real((morder+1),fp_kind)*duf
     endif
     sum(3) = sum(3) + 1.5_fp_kind + (1.5_fp_kind-real((morder),fp_kind))*ug
     sum(2) = sum(2) + 1._fp_kind - real((morder+1),fp_kind)*uf
     sum(1) = fdf*sum(1)
     ! Transform to tcl derivative from ln g
     ! derivative recalling that dlng/dlntc = 1, and dlng/dlnf = 0.5*uf.
     if(i.le.2) then
        sum(4) = sum(4) +&
             0.5_fp_kind*(duf*sum(3) + uf*sum(5)) +&
             0.5_fp_kind*uf*(sum(5) + 0.5_fp_kind*uf*sum(6))
        sum(7) = sum(7) +&
             0.5_fp_kind*(duf2*sum(3) + 2._fp_kind*duf*sum(5) + uf*sum(8)) +&
             0.5_fp_kind*(duf*sum(5) + uf*sum(8) +&
             0.5_fp_kind*uf*(2._fp_kind*duf*sum(6) + uf*sum(9))) +&
             0.5_fp_kind*uf*(sum(8) +&
             0.5_fp_kind*(duf*sum(6) + uf*sum(9)) +&
             0.5_fp_kind*uf*(sum(9) + 0.5_fp_kind*uf*sum(10)))
        ! do sum(5) later than sum(7) so that rhs of sum(7) is undisturbed.
        sum(5) = sum(5) + 0.5_fp_kind*uf*sum(6)
        sum(8) = sum(8) +&
             0.5_fp_kind*(duf*sum(6) + uf*sum(9)) +&
             0.5_fp_kind*uf*(sum(9) + 0.5_fp_kind*uf*sum(10))
        sum(9) = sum(9) + 0.5_fp_kind*uf*sum(10)
     endif
     sum(2) = sum(2) + 0.5_fp_kind*uf*sum(3)
     if(i.eq.1) then
        rhostar(1:maxfd_direct) = sum(1:maxfd_direct)
     elseif(i.eq.2) then
        pstar(1:maxfd_direct) = sum(1:maxfd_direct)
     elseif(i.eq.3) then
        sstar(1:mstar) = sum(1:mstar)
     endif
     ! Just in case the compiler is too crude to implicitly deallocate sum when it falls out
     ! of scope for each iteration of this do loop.
     deallocate(sum)
  enddo !do i = 1,3

  ! psi = 2 wf + log((wf-1.d0)/(wf+1.d0)) = 2 wf - logarg
  ! transform to new functions of fl, tcl.
  ! N.B. up to now sstar has stored qstar and its derivatives, transform
  ! to derivative of ln se wrt ln f and ln tc.
  ! N.B. qstar in this programme is defined with different convention
  ! than paper, see initial commentary.
  qconst = sstar(1)/(rhostar(1)*wf)
  ! sstar(1) = qconst + 2 wf - psi
  ! = qconst + logarg
  sstar(1) = qconst + logarg
  sstar(2) = (qconst*(sstar(2)-rhostar(2)-0.5_fp_kind*uf)+dlogarg)/sstar(1)
  sstar(3) = (qconst*(sstar(3)-rhostar(3)))/sstar(1)
  pconst = pstar(1)*wf/rhostar(1)

  ustar(1) = sstar(1)+psi-pconst
  ustar(2) = (sstar(1)*sstar(2)+wf-pconst*(pstar(2)+0.5_fp_kind*uf-rhostar(2)))/ustar(1)
  ustar(3) = (sstar(1)*sstar(3)-pconst*(pstar(3)-rhostar(3)))/ustar(1)

  ! transform pstar = pstar*g
  pstar(1) = (tc*wf)*pstar(1)
  pstar(2) = pstar(2) + 0.5_fp_kind*uf
  pstar(3) = pstar(3) + 1._fp_kind
  pstar(4) = pstar(4) + 0.5_fp_kind*duf
  pstar(7) = pstar(7) + 0.5_fp_kind*duf2
end subroutine fermi_dirac
