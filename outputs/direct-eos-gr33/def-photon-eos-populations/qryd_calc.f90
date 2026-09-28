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

! calculate Rydberg partition function sums from the range
! of minimum principal quantum numbers (given by nmin to nmin_max)
! to nmax of the sum 2 n^2 exp(alpha/n^2) times
! planck_larkin occupation probability (ifpl true) times
! mhd occupation probability (ifmhd true).
! The planck-larkin occupation probability is given by (1 - exp(-a)*(1+a)).
! The *ln* of the mhd occupation probability is given by
!   -[b(1) + b(2)*n^2*(1+g(n)) + b(3)*n^4*(1+g(n))^2 +
!   b(4)*n^6*(1+g(n))^3 + b(5)*n^(15/2)*(1+h(n))].
!   g(n) = 1/(2n).
!   h(n) = (16/(3*n*K_n))^(3/2) - 1
!   = (16/(3*n))^(3/2) - 1 [for n <= 3]
!   = ((n+1)/n)^3*((n^2 + n + 1/2)/(n^2 + 7/6*n))^(3/2) - 1 [for n >=3]

! input quantities:
! ifpl (logical) controls whether to use planck-larkin occupation probability.
! ifmhd (logical) controls whether to use mhd occupation probability.
! ifneutral (logical) controls when b(1) through b(4) are employed
!   in the occupation probability calculation.
! eps_factor is a factor used to help terminate the principal
!   quantum number sum. eps_factor = exp(-c2*min(chi)/t), where
!   min(chi) = the minimum ionization potential for all species
!   with the same a value (i.e. charge).  This definition means
!   that at the limit summand contributes approximately eps (see
!   parameter below) to ln(1+q_excited/q_ground).  With exponential
!   cutoff (MHD), this insures maximum cutoff errors are of order eps.
!   if not MHD (i.e., only Planck-Larkin occupation probability) the
!   remainder of the sum is n*summand.  With nmax of order 300000
!   this means total relative error is of order 1.d-4 which is still
!   fine.
! nmin to nmin_max is the range of minimum principal quantum numbers required.
! nmax is the maximum principal quantum number.  In Planck-Larkin case
!    obtain relative errors of 1.d-5 for nmin of 3 if nmax = 300000.
! a = c2*Z^2*R/T.
! b(nb), where nb = 5 determines the mhd occupation probability.

! output quantities:
! qryd(nmin_max), qryda(nmin_max), qrydb(nb, nmin_max),
! qryda2(nmin_max), qrydab(nb, nmin_max), qrydb2(nb, nb, nmin_max) are
!   the resulting sum plus partial derivatives wrt a and b.
! nmax_reached is the returned maximum principal quantum number
!   that is used in the sum taking account of the convergence criteria.

!> This qryd_calc subroutine calculates (using explicit summation
!> without an approximation) occupation-weighted partition function
!> sums and their first and second (mixed) partial derivatives wrt to
!> b (a convenient transformation of the subset of auxiliary variables
!> that are relevant) and a (a convenient transformation of t) for
!> Rydberg excited states.
!>
!> \param[in] ifpl PARAMETERS NEED DOCUMENTATION
!>
subroutine qryd_calc(ifpl, ifmhd, ifneutral, eps_factor, nmin, nmax, a, b,&
     qryd, qryda, qrydb, qryda2, qrydab, qrydb2, nmax_reached)

  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  logical, intent(in) :: ifpl, ifmhd, ifneutral
  integer, intent(in) :: nmin, nmax
  integer, intent(out) :: nmax_reached
  real(fp_kind), intent(in) :: eps_factor, a, b(:)
  real(fp_kind), intent(out) :: qryd(:), qryda(:), qrydb(:, :), qryda2(:), qrydab(:, :), qrydb2(:, :, :)

  ! Internal variables

  real(fp_kind), parameter :: eps=1.e-10_fp_kind

  integer nmin_max, nb
  integer n
  real(fp_kind) rn2, arg, summand, summanda, summanda2,&
       kn, hprime, gprime, occupation, lnoccupation

  real(fp_kind), allocatable :: dlnoccupation(:), dsummand(:), dsummanda(:), dsummand2(:,:)


  nb = size(b)
  nmin_max = size(qryd)

  ! Sanity checking.
  if(&
       nb.ne.size(qrydb,1).or.&
       nb.ne.size(qrydab,1).or.&
       nb.ne.size(qrydb2,1).or.&
       nb.ne.size(qrydb2,2))&
       error stop 'qryd_calc: inconsistent nb dimensions for b, qrydb, qrydab, or qrydb2'
  if(&
       nmin_max.ne.size(qryda).or.&
       nmin_max.ne.size(qrydb,2).or.&
       nmin_max.ne.size(qryda2).or.&
       nmin_max.ne.size(qrydab,2).or.&
       nmin_max.ne.size(qrydb2,3))&
       error stop 'qryd_calc: inconsistent nmin_max dimensions for qryd, qryda, qrydb, qryda2, qrydab, or qrydb2'
  if(nmin.lt.1.or.nmin.gt.nmin_max) error stop 'qryd_calc: bad nmin'

  allocate(dlnoccupation(nb), dsummand(nb), dsummanda(nb), dsummand2(nb,nb))

  if(if_taint_allocated_real) then
     call taint_allocated_real(dlnoccupation)
     call taint_allocated_real(dsummand)
     call taint_allocated_real(dsummanda)
     call taint_allocated_real(dsummand2)
  endif

  nmax_reached = 0
  ! define to allow first time through while loop
  lnoccupation = 1._fp_kind
  summand = 1._fp_kind
  ! n.b. note that returns zero for all values if nmax < nmin, which
  ! is proper thing to do.
  n = nmin
  do while(n.le.nmax.and.&
       (n.eq.nmin.or.(summand.gt.eps/eps_factor.or.&
       (ifmhd.and.lnoccupation.ge.-1._fp_kind))))
     ! avoid integer overflow by calculating this value using real(n,fp_kind)
     rn2 = real(n,fp_kind)*real(n,fp_kind)
     arg = a/rn2
     if(ifpl) then
        call plsummand_normalized(rn2, arg, summand, summanda, summanda2)
     else
        summand = 2._fp_kind*rn2*exp(arg)
        summanda = 2._fp_kind*exp(arg)
        summanda2 = 2._fp_kind*exp(arg)/rn2
     endif
     if(ifmhd) then
        ! quantum correction K_n see Hummer and Mihalas eq. 4.24
        if(n.le.3) then
           kn = 1._fp_kind
        else
           kn =&
                (16._fp_kind*rn2*(real((n),fp_kind) + 7._fp_kind/6._fp_kind))/&
                (real((3*(n+1)),fp_kind)*real((n+1),fp_kind)*(rn2 + real((n),fp_kind) + 0.5_fp_kind))
        endif
        hprime = (16._fp_kind/(3._fp_kind*real((n),fp_kind)*kn))*sqrt(16._fp_kind/(3._fp_kind*real((n),fp_kind)*kn))
        ! mhd occupation probability for neutral-ion and ion-ion interactions
        dlnoccupation(5) = -(real((n),fp_kind))**7*sqrt(real((n),fp_kind))*hprime
        lnoccupation = b(5)*dlnoccupation(5)
        if(ifneutral) then
           ! mhd occupation probability for neutral-neutral interactions
           ! l = l_max = n-1 --> 1/2 * (3n^2 - l(l+1)) = n^2*(1+1/2n)
           gprime = (1._fp_kind + 0.5_fp_kind/real((n),fp_kind))*rn2
           ! l = l_min = 0 --> 1/2 * (3n^2 - l(l+1)) = n^2*(3/2) (test case)
           ! gprime = (1.5d0)*rn2
           lnoccupation = lnoccupation - (b(1) + gprime*(b(2) + gprime*(b(3) + gprime*b(4))))
           dlnoccupation(1) = -1._fp_kind
           dlnoccupation(2) = dlnoccupation(1)*gprime
           dlnoccupation(3) = dlnoccupation(2)*gprime
           dlnoccupation(4) = dlnoccupation(3)*gprime
        endif
        occupation = exp(lnoccupation)
        summand = summand*occupation
        summanda = summanda*occupation
        summanda2 = summanda2*occupation
        dsummand(5) = summand*dlnoccupation(5)
        dsummanda(5) = summanda*dlnoccupation(5)
        dsummand2(5,5) = summand*dlnoccupation(5)*dlnoccupation(5)
        if(ifneutral) then
           dsummand(1) = summand*dlnoccupation(1)
           dsummand(2) = summand*dlnoccupation(2)
           dsummand(3) = summand*dlnoccupation(3)
           dsummand(4) = summand*dlnoccupation(4)
           dsummanda(1) = summanda*dlnoccupation(1)
           dsummanda(2) = summanda*dlnoccupation(2)
           dsummanda(3) = summanda*dlnoccupation(3)
           dsummanda(4) = summanda*dlnoccupation(4)
           dsummand2(1,1) = summand*dlnoccupation(1)*dlnoccupation(1)
           dsummand2(2,1) = summand*dlnoccupation(2)*dlnoccupation(1)
           dsummand2(3,1) = summand*dlnoccupation(3)*dlnoccupation(1)
           dsummand2(4,1) = summand*dlnoccupation(4)*dlnoccupation(1)
           dsummand2(5,1) = summand*dlnoccupation(5)*dlnoccupation(1)
           dsummand2(2,2) = summand*dlnoccupation(2)*dlnoccupation(2)
           dsummand2(3,2) = summand*dlnoccupation(3)*dlnoccupation(2)
           dsummand2(4,2) = summand*dlnoccupation(4)*dlnoccupation(2)
           dsummand2(5,2) = summand*dlnoccupation(5)*dlnoccupation(2)
           dsummand2(3,3) = summand*dlnoccupation(3)*dlnoccupation(3)
           dsummand2(4,3) = summand*dlnoccupation(4)*dlnoccupation(3)
           dsummand2(5,3) = summand*dlnoccupation(5)*dlnoccupation(3)
           dsummand2(4,4) = summand*dlnoccupation(4)*dlnoccupation(4)
           dsummand2(5,4) = summand*dlnoccupation(5)*dlnoccupation(4)
        endif
        if(n.le.nmin_max) then
           qrydb(5,n) = dsummand(5)
           qrydab(5,n) = dsummanda(5)
           qrydb2(5,5,n) = dsummand2(5,5)
           if(ifneutral) then
              qrydb(1,n) = dsummand(1)
              qrydb(2,n) = dsummand(2)
              qrydb(3,n) = dsummand(3)
              qrydb(4,n) = dsummand(4)
              qrydab(1,n) = dsummanda(1)
              qrydab(2,n) = dsummanda(2)
              qrydab(3,n) = dsummanda(3)
              qrydab(4,n) = dsummanda(4)
              qrydb2(1,1,n) = dsummand2(1,1)
              qrydb2(2,1,n) = dsummand2(2,1)
              qrydb2(3,1,n) = dsummand2(3,1)
              qrydb2(4,1,n) = dsummand2(4,1)
              qrydb2(5,1,n) = dsummand2(5,1)
              qrydb2(2,2,n) = dsummand2(2,2)
              qrydb2(3,2,n) = dsummand2(3,2)
              qrydb2(4,2,n) = dsummand2(4,2)
              qrydb2(5,2,n) = dsummand2(5,2)
              qrydb2(3,3,n) = dsummand2(3,3)
              qrydb2(4,3,n) = dsummand2(4,3)
              qrydb2(5,3,n) = dsummand2(5,3)
              qrydb2(4,4,n) = dsummand2(4,4)
              qrydb2(5,4,n) = dsummand2(5,4)
           endif
        else
           qrydb(5,nmin_max) = qrydb(5,nmin_max) +&
                dsummand(5)
           qrydab(5,nmin_max) = qrydab(5,nmin_max) +&
                dsummanda(5)
           qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +&
                dsummand2(5,5)
           if(ifneutral) then
              qrydb(1,nmin_max) = qrydb(1,nmin_max) +&
                   dsummand(1)
              qrydb(2,nmin_max) = qrydb(2,nmin_max) +&
                   dsummand(2)
              qrydb(3,nmin_max) = qrydb(3,nmin_max) +&
                   dsummand(3)
              qrydb(4,nmin_max) = qrydb(4,nmin_max) +&
                   dsummand(4)
              qrydab(1,nmin_max) = qrydab(1,nmin_max) +&
                   dsummanda(1)
              qrydab(2,nmin_max) = qrydab(2,nmin_max) +&
                   dsummanda(2)
              qrydab(3,nmin_max) = qrydab(3,nmin_max) +&
                   dsummanda(3)
              qrydab(4,nmin_max) = qrydab(4,nmin_max) +&
                   dsummanda(4)
              qrydb2(1,1,nmin_max) = qrydb2(1,1,nmin_max) +&
                   dsummand2(1,1)
              qrydb2(2,1,nmin_max) = qrydb2(2,1,nmin_max) +&
                   dsummand2(2,1)
              qrydb2(3,1,nmin_max) = qrydb2(3,1,nmin_max) +&
                   dsummand2(3,1)
              qrydb2(4,1,nmin_max) = qrydb2(4,1,nmin_max) +&
                   dsummand2(4,1)
              qrydb2(5,1,nmin_max) = qrydb2(5,1,nmin_max) +&
                   dsummand2(5,1)
              qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +&
                   dsummand2(2,2)
              qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +&
                   dsummand2(3,2)
              qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +&
                   dsummand2(4,2)
              qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +&
                   dsummand2(5,2)
              qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +&
                   dsummand2(3,3)
              qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +&
                   dsummand2(4,3)
              qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +&
                   dsummand2(5,3)
              qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +&
                   dsummand2(4,4)
              qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +&
                   dsummand2(5,4)
           endif
        endif
     endif
     if(n.le.nmin_max) then
        qryd(n) = summand
        qryda(n) = summanda
        qryda2(n) = summanda2
     else
        qryd(nmin_max) = qryd(nmin_max) + summand
        qryda(nmin_max) = qryda(nmin_max) + summanda
        qryda2(nmin_max) = qryda2(nmin_max) + summanda2
     endif
     n = n + 1
  enddo
  nmax_reached = max(nmax_reached, n-1)

  ! zero remainder of qryd if necessary.
  qryd(n:nmin_max) = 0._fp_kind
  qryda(n:nmin_max) = 0._fp_kind
  qryda2(n:nmin_max) = 0._fp_kind
  if(ifmhd) then
     qrydb(5,n:nmin_max) = 0._fp_kind
     qrydab(5,n:nmin_max) = 0._fp_kind
     qrydb2(5,5,n:nmin_max) = 0._fp_kind
     if(ifneutral) then
        qrydb(1,n:nmin_max) = 0._fp_kind
        qrydb(2,n:nmin_max) = 0._fp_kind
        qrydb(3,n:nmin_max) = 0._fp_kind
        qrydb(4,n:nmin_max) = 0._fp_kind
        qrydab(1,n:nmin_max) = 0._fp_kind
        qrydab(2,n:nmin_max) = 0._fp_kind
        qrydab(3,n:nmin_max) = 0._fp_kind
        qrydab(4,n:nmin_max) = 0._fp_kind
        qrydb2(1,1,n:nmin_max) = 0._fp_kind
        qrydb2(2,1,n:nmin_max) = 0._fp_kind
        qrydb2(3,1,n:nmin_max) = 0._fp_kind
        qrydb2(4,1,n:nmin_max) = 0._fp_kind
        qrydb2(5,1,n:nmin_max) = 0._fp_kind
        qrydb2(2,2,n:nmin_max) = 0._fp_kind
        qrydb2(3,2,n:nmin_max) = 0._fp_kind
        qrydb2(4,2,n:nmin_max) = 0._fp_kind
        qrydb2(5,2,n:nmin_max) = 0._fp_kind
        qrydb2(3,3,n:nmin_max) = 0._fp_kind
        qrydb2(4,3,n:nmin_max) = 0._fp_kind
        qrydb2(5,3,n:nmin_max) = 0._fp_kind
        qrydb2(4,4,n:nmin_max) = 0._fp_kind
        qrydb2(5,4,n:nmin_max) = 0._fp_kind
     endif
  endif

  ! convert from summand to sum for lower indices
  qryd(nmin_max-1:nmin:-1) = qryd(nmin_max-1:nmin:-1) + qryd(nmin_max:nmin+1:-1)
  qryda(nmin_max-1:nmin:-1) = qryda(nmin_max-1:nmin:-1) + qryda(nmin_max:nmin+1:-1)
  qryda2(nmin_max-1:nmin:-1) = qryda2(nmin_max-1:nmin:-1) + qryda2(nmin_max:nmin+1:-1)
  if(ifmhd) then
     qrydb(5,nmin_max-1:nmin:-1) = qrydb(5,nmin_max-1:nmin:-1) + qrydb(5,nmin_max:nmin+1:-1)
     qrydab(5,nmin_max-1:nmin:-1) = qrydab(5,nmin_max-1:nmin:-1) + qrydab(5,nmin_max:nmin+1:-1)
     qrydb2(5,5,nmin_max-1:nmin:-1) = qrydb2(5,5,nmin_max-1:nmin:-1) + qrydb2(5,5,nmin_max:nmin+1:-1)
     if(ifneutral) then
        qrydb(1,nmin_max-1:nmin:-1) = qrydb(1,nmin_max-1:nmin:-1) + qrydb(1,nmin_max:nmin+1:-1)
        qrydb(2,nmin_max-1:nmin:-1) = qrydb(2,nmin_max-1:nmin:-1) + qrydb(2,nmin_max:nmin+1:-1)
        qrydb(3,nmin_max-1:nmin:-1) = qrydb(3,nmin_max-1:nmin:-1) + qrydb(3,nmin_max:nmin+1:-1)
        qrydb(4,nmin_max-1:nmin:-1) = qrydb(4,nmin_max-1:nmin:-1) + qrydb(4,nmin_max:nmin+1:-1)
        qrydab(1,nmin_max-1:nmin:-1) = qrydab(1,nmin_max-1:nmin:-1) + qrydab(1,nmin_max:nmin+1:-1)
        qrydab(2,nmin_max-1:nmin:-1) = qrydab(2,nmin_max-1:nmin:-1) + qrydab(2,nmin_max:nmin+1:-1)
        qrydab(3,nmin_max-1:nmin:-1) = qrydab(3,nmin_max-1:nmin:-1) + qrydab(3,nmin_max:nmin+1:-1)
        qrydab(4,nmin_max-1:nmin:-1) = qrydab(4,nmin_max-1:nmin:-1) + qrydab(4,nmin_max:nmin+1:-1)
        qrydb2(1,1,nmin_max-1:nmin:-1) = qrydb2(1,1,nmin_max-1:nmin:-1) + qrydb2(1,1,nmin_max:nmin+1:-1)
        qrydb2(2,1,nmin_max-1:nmin:-1) = qrydb2(2,1,nmin_max-1:nmin:-1) + qrydb2(2,1,nmin_max:nmin+1:-1)
        qrydb2(3,1,nmin_max-1:nmin:-1) = qrydb2(3,1,nmin_max-1:nmin:-1) + qrydb2(3,1,nmin_max:nmin+1:-1)
        qrydb2(4,1,nmin_max-1:nmin:-1) = qrydb2(4,1,nmin_max-1:nmin:-1) + qrydb2(4,1,nmin_max:nmin+1:-1)
        qrydb2(5,1,nmin_max-1:nmin:-1) = qrydb2(5,1,nmin_max-1:nmin:-1) + qrydb2(5,1,nmin_max:nmin+1:-1)
        qrydb2(2,2,nmin_max-1:nmin:-1) = qrydb2(2,2,nmin_max-1:nmin:-1) + qrydb2(2,2,nmin_max:nmin+1:-1)
        qrydb2(3,2,nmin_max-1:nmin:-1) = qrydb2(3,2,nmin_max-1:nmin:-1) + qrydb2(3,2,nmin_max:nmin+1:-1)
        qrydb2(4,2,nmin_max-1:nmin:-1) = qrydb2(4,2,nmin_max-1:nmin:-1) + qrydb2(4,2,nmin_max:nmin+1:-1)
        qrydb2(5,2,nmin_max-1:nmin:-1) = qrydb2(5,2,nmin_max-1:nmin:-1) + qrydb2(5,2,nmin_max:nmin+1:-1)
        qrydb2(3,3,nmin_max-1:nmin:-1) = qrydb2(3,3,nmin_max-1:nmin:-1) + qrydb2(3,3,nmin_max:nmin+1:-1)
        qrydb2(4,3,nmin_max-1:nmin:-1) = qrydb2(4,3,nmin_max-1:nmin:-1) + qrydb2(4,3,nmin_max:nmin+1:-1)
        qrydb2(5,3,nmin_max-1:nmin:-1) = qrydb2(5,3,nmin_max-1:nmin:-1) + qrydb2(5,3,nmin_max:nmin+1:-1)
        qrydb2(4,4,nmin_max-1:nmin:-1) = qrydb2(4,4,nmin_max-1:nmin:-1) + qrydb2(4,4,nmin_max:nmin+1:-1)
        qrydb2(5,4,nmin_max-1:nmin:-1) = qrydb2(5,4,nmin_max-1:nmin:-1) + qrydb2(5,4,nmin_max:nmin+1:-1)
     endif
  endif
end subroutine qryd_calc
