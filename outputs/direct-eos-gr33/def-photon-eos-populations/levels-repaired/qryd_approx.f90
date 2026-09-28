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
!
!*******************************************************************************

! use explicit summations and calls to approximation routines  to:
! calculate Rydberg partition function sums from the range
! of minimum principal quantum numbers (given by nmin to nmin_max)
! to nmax of 2 n^2 exp(a/n^2) times
! planck_larkin occupation probability (ifpl true) times
! mhd occupation probability (ifmhd true).
! The planck-larkin occupation probability is given by (1 - exp(-a)*(1+a)).
! The *ln* of the mhd occupation probability is given by
!   -[b(1) + b(2)*n^2*(1+g(n)) + b(3)*n^4*(1+g(n))^2 +
!   b(4)*n^6*(1+g(n))^3 + b(5)*n^(15/2)*(1+h(n))].
!   g(n) = 1/(2n).
!   h(n) = (16/(3*n*K_n))^(3/2) - 1
!        = (16/(3*n))^(3/2) - 1 [for n <= 3]
!        = ((n+1)/n)^3*((n^2 + n + 1/2)/(n^2 + 7/6*n))^(3/2) - 1 [for n >=3]

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
! nmin to nmin_max is the range of required minimum principal quantum numbers.
! nmax is the maximum principal quantum number.  In Planck-Larkin case
!   obtain relative errors of 1.d-5 for nmin of 3 if nmax = 300000.
! a = c2*Z^2*R/T.
! b(nb), (nb = 5) determines the mhd occupation probability.

! output quantities:
! qryd(nmin_max), qryda(nmin_max), qrydb(nb, nmin_max),
! qryda2(nmin_max), qrydab(nb, nmin_max),
! qrydab2(nb,nb,nmin_max) (lower triangle in nb) are
!   the resulting sum plus partial derivatives wrt a and b.
! nmax_reached is the returned maximum principal quantum number
!   that is used in the sum taking account of the convergence criteria.

!> This qryd_approx subroutine calculates (using an approximation)
!> occupation-weighted partition function sums and their first and
!> second (mixed) partial derivatives wrt to b (a convenient
!> transformation of the subset of auxiliary variables that are
!> relevant) and a (a convenient transformation of t) for Rydberg
!> excited states.
!>
!> \param[in] ifpl PARAMETERS NEED DOCUMENTATION
!>
subroutine qryd_approx(ifpl, ifmhd, ifneutral, eps_factor,&
     nmin, nmax, a, b,&
     qryd, qryda, qrydb, qryda2, qrydab, qrydb2, nmax_reached)

  use mod_free_eos_constants, only: pi
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! arguments
  logical, intent(in) :: ifpl, ifmhd, ifneutral
  integer, intent(in) :: nmin, nmax
  integer, intent(out) :: nmax_reached
  real(fp_kind), intent(in) :: eps_factor, a, b(:)
  real(fp_kind), intent(out) :: qryd(:), qryda(:), qrydb(:,:), qryda2(:), qrydab(:,:), qrydb2(:,:,:)

  ! internal variables:

  real(fp_kind), parameter :: eps=1.e-10_fp_kind
  ! this is a standard value for most of the approximations
  ! which is used to limit their range of applicability
  real(fp_kind), parameter ::exp_max = 1._fp_kind
  ! limits on continuous interpolation between exact sum
  ! and approximation
  real(fp_kind), parameter :: lim_sum = -1.e-3_fp_kind
  real(fp_kind), parameter :: lim_approx = -1.e-4_fp_kind

  integer nmin_max, nb
  integer n, nmaxa

  real(fp_kind) rn2, arg, summand, summanda, summanda2,&
       kn, hprime, gprime, occupation,&
       lnoccupation, lnoccupationa
  real(fp_kind) expmb1
  real(fp_kind) sarg, carg, dsarg, d2sarg,&
       wapprox, wsum

  real(fp_kind), allocatable ::&
       dlnoccupation(:),&
       dsummand(:),&
       dsummanda(:),&
       dsummand2(:,:),&
       lqryd(:),&
       lqryda(:),&
       lqrydb(:),&
       lqryda2(:),&
       lqrydab(:),&
       lqrydb2(:,:),&
       dlnoccupationa(:),&
       dwapprox(:),&
       dwapprox2(:,:),&
       dwsum(:),&
       dwsum2(:,:)

  nb = size(b)
  nmin_max = size(qryd)

  if(&
       nb.ne.size(qrydb,1).or.&
       nb.ne.size(qrydab,1).or.&
       nb.ne.size(qrydb2,1).or.&
       nb.ne.size(qrydb2,2))&
       error stop 'qryd_approx: inconsistent nb dimensions for b, qrydb, qrydab, or qrydb2'
  if(&
       nmin_max.ne.size(qryda).or.&
       nmin_max.ne.size(qrydb,2).or.&
       nmin_max.ne.size(qryda2).or.&
       nmin_max.ne.size(qrydab,2).or.&
       nmin_max.ne.size(qrydb2,3))&
       error stop 'qryd_approx: inconsistent nmin_max dimensions for qryd, qryda, qrydb, qryda2, qrydab, or qrydb2'
  if(nmin.lt.1.or.nmin.gt.nmin_max) error stop 'qryd_approx: bad nmin'

  allocate(&
       dlnoccupation(nb),&
       dsummand(nb),&
       dsummanda(nb),&
       dsummand2(nb,nb),&
       lqryd(1),&
       lqryda(1),&
       lqrydb(nb),&
       lqryda2(1),&
       lqrydab(nb),&
       lqrydb2(nb,nb),&
       dlnoccupationa(nb),&
       dwapprox(nb),&
       dwapprox2(nb,nb),&
       dwsum(nb),&
       dwsum2(nb,nb))

  if(if_taint_allocated_real) then
       call taint_allocated_real(dlnoccupation)
       call taint_allocated_real(dsummand)
       call taint_allocated_real(dsummanda)
       call taint_allocated_real(dsummand2)
       call taint_allocated_real(lqryd)
       call taint_allocated_real(lqryda)
       call taint_allocated_real(lqrydb)
       call taint_allocated_real(lqryda2)
       call taint_allocated_real(lqrydab)
       call taint_allocated_real(lqrydb2)
       call taint_allocated_real(dlnoccupationa)
       call taint_allocated_real(dwapprox)
       call taint_allocated_real(dwapprox2)
       call taint_allocated_real(dwsum)
       call taint_allocated_real(dwsum2)
  endif

  nmax_reached = 0
  if(.not.ifmhd) then
     if(.not.ifpl) error stop 'qryd_approx: bad values of ifmhd or ifpl'
     call plsum(nmin, nmin_max, 0, a, qryd, qryda, qryda2)
     if(nmax.lt.nmin_max) then
        error stop 'qryd_approx: this case not programmed'
     elseif(nmax.le.300000) then
        call plsum(nmax+1, nmax+1, nmax, a, lqryd, lqryda, lqryda2)
        qryd(nmin:nmin_max) = qryd(nmin:nmin_max) - lqryd(1)
        qryda(nmin:nmin_max) = qryda(nmin:nmin_max) - lqryda(1)
        qryda2(nmin:nmin_max) = qryda2(nmin:nmin_max) - lqryda2(1)
     endif
  else
     error stop 'qryd_approx: excitation approximation is invalid'
     ! n.b. the code below works in many cases, but there seems to
     ! be bad significance loss (or else partial derivative troubles)
     ! plaguing it and it is slow.  Thus, should be using truncation
     ! approximation instead (see paper).
     ! n.b. to reduce the size of the resulting executable, I have
     ! commented out calls to approximation routines and removed them
     ! from this directory.
     ! adopt minimum limit since approximations
     ! only work for nmin_max.ge.3 because of change in definition of
     ! kn at n = 3.
     if(nmin_max.lt.3)&
          error stop 'qryd_approx: nmin_max must be 3 or greater for ifmhd'
     ! n must be >= to this number in order to use approximations
     ! max(nmin_max,... assures at least one loop for negligible a.
     nmaxa = max(nmin_max,int(sqrt(a/exp_max))+1)
     n = nmin
     ! define these quantities to clear out any undefined garbage
     summand = 0._fp_kind
     lnoccupation = 0._fp_kind
     lnoccupationa = 0._fp_kind
     ! n.b. nmax and summand limits used below only for *very* low
     ! temperatures where nmaxa can get large.  In particular nmax
     ! is usually a no-op for the MHD case where you are using
     ! approximations.  One could re-programme to make nmax be
     ! effective for part of approximation where not in cutoff region,
     ! but I have judged this is not worth it since ordinarily ignore
     ! all approximations for MHD case.
     do while(.not.(n.gt.nmaxa.and.lnoccupation.gt.lim_approx).and.&
          n.le.nmax.and.&
          (n.eq.nmin.or.summand.gt.eps/eps_factor))
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
        ! quantum correction K_n see Hummer and Mihalas eq. 4.24
        if(n.le.3) then
           kn = 1._fp_kind
        else
           kn =&
                (16._fp_kind*rn2*(real((n),fp_kind) + 7._fp_kind/6._fp_kind))/&
                (real((3*(n+1)),fp_kind)*real((n+1),fp_kind)*(rn2 + real((n),fp_kind) + 0.5_fp_kind))
        endif
        hprime = (16._fp_kind/(3._fp_kind*real((n),fp_kind)*kn))*sqrt(16._fp_kind/(3._fp_kind*real((n),fp_kind)*kn))
        ! mhd ln occupation probability for neutral-ion and
        ! ion-ion interactions
        dlnoccupation(5) = -(real((n),fp_kind))**7*sqrt(real((n),fp_kind))*hprime
        lnoccupation = b(5)*dlnoccupation(5)
        if(ifneutral) then
           ! mhd ln occupation probability for neutral-neutral interactions
           gprime = (1._fp_kind + 0.5_fp_kind/real((n),fp_kind))*rn2
           ! n.b. b(1) applied later.
           lnoccupation = lnoccupation -&
                (gprime*(b(2) + gprime*(b(3) + gprime*b(4))))
           dlnoccupation(2) = -gprime
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
           dsummand(2) = summand*dlnoccupation(2)
           dsummand(3) = summand*dlnoccupation(3)
           dsummand(4) = summand*dlnoccupation(4)
           dsummanda(2) = summanda*dlnoccupation(2)
           dsummanda(3) = summanda*dlnoccupation(3)
           dsummanda(4) = summanda*dlnoccupation(4)
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
        if(n.eq.nmaxa) then
           lnoccupationa = lnoccupation
           dlnoccupationa(5) = dlnoccupation(5)
           if(ifneutral) then
              dlnoccupationa(2) = dlnoccupation(2)
              dlnoccupationa(3) = dlnoccupation(3)
              dlnoccupationa(4) = dlnoccupation(4)
           endif
        endif
        if(n.ge.nmaxa) then
           ! n.b. by definition of nmaxa, n >= nmin_max.
           ! n.n.b. qryd contains sum from nmin_max to nmaxa-1 of summand,
           ! thus, must zero it in this  special case.
           if(nmaxa.eq.nmin_max) then
              qryd(nmin_max) = 0._fp_kind
              qryda(nmin_max) = 0._fp_kind
              qryda2(nmin_max) = 0._fp_kind
              qrydb(5,nmin_max) = 0._fp_kind
              qrydab(5,nmin_max) = 0._fp_kind
              qrydb2(5,5,nmin_max) = 0._fp_kind
              if(ifneutral) then
                 qrydb(2,nmin_max) = 0._fp_kind
                 qrydb(3,nmin_max) = 0._fp_kind
                 qrydb(4,nmin_max) = 0._fp_kind
                 qrydab(2,nmin_max) = 0._fp_kind
                 qrydab(3,nmin_max) = 0._fp_kind
                 qrydab(4,nmin_max) = 0._fp_kind
                 qrydb2(2,2,nmin_max) = 0._fp_kind
                 qrydb2(3,2,nmin_max) = 0._fp_kind
                 qrydb2(4,2,nmin_max) = 0._fp_kind
                 qrydb2(5,2,nmin_max) = 0._fp_kind
                 qrydb2(3,3,nmin_max) = 0._fp_kind
                 qrydb2(4,3,nmin_max) = 0._fp_kind
                 qrydb2(5,3,nmin_max) = 0._fp_kind
                 qrydb2(4,4,nmin_max) = 0._fp_kind
                 qrydb2(5,4,nmin_max) = 0._fp_kind
              endif
           endif
           if(n.eq.nmaxa) then
              lqryd(1) = summand
              lqryda(1) = summanda
              lqryda2(1) = summanda2
              lqrydb(5) = dsummand(5)
              lqrydab(5) = dsummanda(5)
              lqrydb2(5,5) = dsummand2(5,5)
              if(ifneutral) then
                 lqrydb(2) = dsummand(2)
                 lqrydb(3) = dsummand(3)
                 lqrydb(4) = dsummand(4)
                 lqrydab(2) = dsummanda(2)
                 lqrydab(3) = dsummanda(3)
                 lqrydab(4) = dsummanda(4)
                 lqrydb2(2,2) = dsummand2(2,2)
                 lqrydb2(3,2) = dsummand2(3,2)
                 lqrydb2(4,2) = dsummand2(4,2)
                 lqrydb2(5,2) = dsummand2(5,2)
                 lqrydb2(3,3) = dsummand2(3,3)
                 lqrydb2(4,3) = dsummand2(4,3)
                 lqrydb2(5,3) = dsummand2(5,3)
                 lqrydb2(4,4) = dsummand2(4,4)
                 lqrydb2(5,4) = dsummand2(5,4)
              endif
           else
              lqryd(1) = lqryd(1) +&
                   summand
              lqryda(1) = lqryda(1) +&
                   summanda
              lqryda2(1) = lqryda2(1) +&
                   summanda2
              lqrydb(5) = lqrydb(5) +&
                   dsummand(5)
              lqrydab(5) = lqrydab(5) +&
                   dsummanda(5)
              lqrydb2(5,5) = lqrydb2(5,5) +&
                   dsummand2(5,5)
              if(ifneutral) then
                 lqrydb(2) = lqrydb(2) +&
                      dsummand(2)
                 lqrydb(3) = lqrydb(3) +&
                      dsummand(3)
                 lqrydb(4) = lqrydb(4) +&
                      dsummand(4)
                 lqrydab(2) = lqrydab(2) +&
                      dsummanda(2)
                 lqrydab(3) = lqrydab(3) +&
                      dsummanda(3)
                 lqrydab(4) = lqrydab(4) +&
                      dsummanda(4)
                 lqrydb2(2,2) = lqrydb2(2,2) +&
                      dsummand2(2,2)
                 lqrydb2(3,2) = lqrydb2(3,2) +&
                      dsummand2(3,2)
                 lqrydb2(4,2) = lqrydb2(4,2) +&
                      dsummand2(4,2)
                 lqrydb2(5,2) = lqrydb2(5,2) +&
                      dsummand2(5,2)
                 lqrydb2(3,3) = lqrydb2(3,3) +&
                      dsummand2(3,3)
                 lqrydb2(4,3) = lqrydb2(4,3) +&
                      dsummand2(4,3)
                 lqrydb2(5,3) = lqrydb2(5,3) +&
                      dsummand2(5,3)
                 lqrydb2(4,4) = lqrydb2(4,4) +&
                      dsummand2(4,4)
                 lqrydb2(5,4) = lqrydb2(5,4) +&
                      dsummand2(5,4)
              endif
           endif
        else
           if(n.le.nmin_max) then
              qryd(n) = summand
              qryda(n) = summanda
              qryda2(n) = summanda2
              qrydb(5,n) = dsummand(5)
              qrydab(5,n) = dsummanda(5)
              qrydb2(5,5,n) = dsummand2(5,5)
              if(ifneutral) then
                 qrydb(2,n) = dsummand(2)
                 qrydb(3,n) = dsummand(3)
                 qrydb(4,n) = dsummand(4)
                 qrydab(2,n) = dsummanda(2)
                 qrydab(3,n) = dsummanda(3)
                 qrydab(4,n) = dsummanda(4)
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
              qryd(nmin_max) = qryd(nmin_max) +&
                   summand
              qryda(nmin_max) = qryda(nmin_max) +&
                   summanda
              qryda2(nmin_max) = qryda2(nmin_max) +&
                   summanda2
              qrydb(5,nmin_max) = qrydb(5,nmin_max) +&
                   dsummand(5)
              qrydab(5,nmin_max) = qrydab(5,nmin_max) +&
                   dsummanda(5)
              qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +&
                   dsummand2(5,5)
              if(ifneutral) then
                 qrydb(2,nmin_max) = qrydb(2,nmin_max) +&
                      dsummand(2)
                 qrydb(3,nmin_max) = qrydb(3,nmin_max) +&
                      dsummand(3)
                 qrydb(4,nmin_max) = qrydb(4,nmin_max) +&
                      dsummand(4)
                 qrydab(2,nmin_max) = qrydab(2,nmin_max) +&
                      dsummanda(2)
                 qrydab(3,nmin_max) = qrydab(3,nmin_max) +&
                      dsummanda(3)
                 qrydab(4,nmin_max) = qrydab(4,nmin_max) +&
                      dsummanda(4)
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
        n = n + 1
     enddo
     nmax_reached = max(nmax_reached, n-1)
     if(n-1.lt.nmaxa) then
        ! in this case, sum terminated by nmax or summand conditions.
        ! only remaining thing to do is zero remainder of qryd if necessary.
        qryd(n:nmin_max) = 0._fp_kind
        qryda(n:nmin_max) = 0._fp_kind
        qryda2(n:nmin_max) = 0._fp_kind
        qrydb(5,n:nmin_max) = 0._fp_kind
        qrydab(5,n:nmin_max) = 0._fp_kind
        qrydb2(5,5,n:nmin_max) = 0._fp_kind
        if(ifneutral) then
           qrydb(2:4,n:nmin_max) = 0._fp_kind
           qrydab(2:4,n:nmin_max) = 0._fp_kind
           qrydb2(2:5,2,n:nmin_max) = 0._fp_kind
           qrydb2(3:5,3,n:nmin_max) = 0._fp_kind
           qrydb2(4:5,4,n:nmin_max) = 0._fp_kind
        endif
        ! n.b. in all other cases n-1 reached nmaxa >= nmin_max, so
        ! (1) lnoccupationa defined
        ! (2) zero (when nmaxa = nmin_max) or
        !   sum from nmin_max to nmaxa-1 of summand stored in qryd(nmin_max)
        ! (3) sum from nmaxa to n-1 of summand stored in lqryd(1)
     elseif(lnoccupationa.le.lim_sum) then
        ! in this case sum terminated by nmax or summand conditions and
        ! to finish simply add lqryd(1) to qrd
        qryd(nmin_max) = qryd(nmin_max) + lqryd(1)
        qryda(nmin_max) = qryda(nmin_max) + lqryda(1)
        qryda2(nmin_max) = qryda2(nmin_max) + lqryda2(1)
        qrydb(5,nmin_max) = qrydb(5,nmin_max) + lqrydb(5)
        qrydab(5,nmin_max) = qrydab(5,nmin_max) + lqrydab(5)
        qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) + lqrydb2(5,5)
        if(ifneutral) then
           qrydb(2,nmin_max) = qrydb(2,nmin_max) + lqrydb(2)
           qrydb(3,nmin_max) = qrydb(3,nmin_max) + lqrydb(3)
           qrydb(4,nmin_max) = qrydb(4,nmin_max) + lqrydb(4)
           qrydab(2,nmin_max) = qrydab(2,nmin_max) + lqrydab(2)
           qrydab(3,nmin_max) = qrydab(3,nmin_max) + lqrydab(3)
           qrydab(4,nmin_max) = qrydab(4,nmin_max) + lqrydab(4)
           qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +&
                lqrydb2(2,2)
           qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +&
                lqrydb2(3,2)
           qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +&
                lqrydb2(4,2)
           qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +&
                lqrydb2(5,2)
           qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +&
                lqrydb2(3,3)
           qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +&
                lqrydb2(4,3)
           qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +&
                lqrydb2(5,3)
           qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +&
                lqrydb2(4,4)
           qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +&
                lqrydb2(5,4)
        endif
     elseif(lnoccupationa.le.lim_approx) then
        ! lim_sum < lnoccupationa <= lim_approx
        ! interpolate between
        ! complete approximation at lnoccupationa = lim_approx and
        ! complete summation at lnoccupationa = lim_sum using interpolation
        ! calculate interpolating function with zero first derivatives
        ! at each limit.
        sarg = 0.5_fp_kind*sin(pi*0.5_fp_kind*(&
             (lnoccupationa-lim_sum) -&
             (lim_approx-lnoccupationa)&
             )/(lim_approx-lim_sum))
        carg = 0.5_fp_kind*cos(pi*0.5_fp_kind*(&
             (lnoccupationa-lim_sum) -&
             (lim_approx-lnoccupationa)&
             )/(lim_approx-lim_sum))
        dsarg = carg*pi/(lim_approx-lim_sum)
        d2sarg = -sarg*pi*pi/&
             ((lim_approx-lim_sum)*(lim_approx-lim_sum))
        wapprox = 0.5_fp_kind + sarg
        wsum  = 0.5_fp_kind - sarg
        dwapprox(5) = dsarg*dlnoccupationa(5)
        dwsum(5) = -dwapprox(5)
        dwapprox2(5,5) = d2sarg*dlnoccupationa(5)*dlnoccupationa(5)
        dwsum2(5,5) = -dwapprox2(5,5)
        qryd(nmin_max) = qryd(nmin_max) + lqryd(1)*wsum
        qryda(nmin_max) = qryda(nmin_max) + lqryda(1)*wsum
        qryda2(nmin_max) = qryda2(nmin_max) + lqryda2(1)*wsum
        qrydb(5,nmin_max) = qrydb(5,nmin_max) +&
             lqrydb(5)*wsum + lqryd(1)*dwsum(5)
        qrydab(5,nmin_max) = qrydab(5,nmin_max) +&
             lqrydab(5)*wsum + lqryda(1)*dwsum(5)
        qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +&
             lqrydb2(5,5)*wsum +&
             lqrydb(5)*dwsum(5) + lqrydb(5)*dwsum(5) +&
             lqryd(1)*dwsum2(5,5)
        if(ifneutral) then
           dwapprox(2) = dsarg*dlnoccupationa(2)
           dwapprox(3) = dsarg*dlnoccupationa(3)
           dwapprox(4) = dsarg*dlnoccupationa(4)
           dwsum(2) = -dwapprox(2)
           dwsum(3) = -dwapprox(3)
           dwsum(4) = -dwapprox(4)
           dwapprox2(2,2) = d2sarg*dlnoccupationa(2)*dlnoccupationa(2)
           dwapprox2(3,2) = d2sarg*dlnoccupationa(3)*dlnoccupationa(2)
           dwapprox2(4,2) = d2sarg*dlnoccupationa(4)*dlnoccupationa(2)
           dwapprox2(5,2) = d2sarg*dlnoccupationa(5)*dlnoccupationa(2)
           dwapprox2(3,3) = d2sarg*dlnoccupationa(3)*dlnoccupationa(3)
           dwapprox2(4,3) = d2sarg*dlnoccupationa(4)*dlnoccupationa(3)
           dwapprox2(5,3) = d2sarg*dlnoccupationa(5)*dlnoccupationa(3)
           dwapprox2(4,4) = d2sarg*dlnoccupationa(4)*dlnoccupationa(4)
           dwapprox2(5,4) = d2sarg*dlnoccupationa(5)*dlnoccupationa(4)
           dwsum2(2,2) = -dwapprox2(2,2)
           dwsum2(3,2) = -dwapprox2(3,2)
           dwsum2(4,2) = -dwapprox2(4,2)
           dwsum2(5,2) = -dwapprox2(5,2)
           dwsum2(3,3) = -dwapprox2(3,3)
           dwsum2(4,3) = -dwapprox2(4,3)
           dwsum2(5,3) = -dwapprox2(5,3)
           dwsum2(4,4) = -dwapprox2(4,4)
           dwsum2(5,4) = -dwapprox2(5,4)
           qrydb(2,nmin_max) = qrydb(2,nmin_max) +&
                lqrydb(2)*wsum + lqryd(1)*dwsum(2)
           qrydb(3,nmin_max) = qrydb(3,nmin_max) +&
                lqrydb(3)*wsum + lqryd(1)*dwsum(3)
           qrydb(4,nmin_max) = qrydb(4,nmin_max) +&
                lqrydb(4)*wsum + lqryd(1)*dwsum(4)
           qrydab(2,nmin_max) = qrydab(2,nmin_max) +&
                lqrydab(2)*wsum + lqryda(1)*dwsum(2)
           qrydab(3,nmin_max) = qrydab(3,nmin_max) +&
                lqrydab(3)*wsum + lqryda(1)*dwsum(3)
           qrydab(4,nmin_max) = qrydab(4,nmin_max) +&
                lqrydab(4)*wsum + lqryda(1)*dwsum(4)
           qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +&
                lqrydb2(2,2)*wsum +&
                lqrydb(2)*dwsum(2) + lqrydb(2)*dwsum(2) +&
                lqryd(1)*dwsum2(2,2)
           qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +&
                lqrydb2(3,2)*wsum +&
                lqrydb(3)*dwsum(2) + lqrydb(2)*dwsum(3) +&
                lqryd(1)*dwsum2(3,2)
           qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +&
                lqrydb2(4,2)*wsum +&
                lqrydb(4)*dwsum(2) + lqrydb(2)*dwsum(4) +&
                lqryd(1)*dwsum2(4,2)
           qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +&
                lqrydb2(5,2)*wsum +&
                lqrydb(5)*dwsum(2) + lqrydb(2)*dwsum(5) +&
                lqryd(1)*dwsum2(5,2)
           qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +&
                lqrydb2(3,3)*wsum +&
                lqrydb(3)*dwsum(3) + lqrydb(3)*dwsum(3) +&
                lqryd(1)*dwsum2(3,3)
           qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +&
                lqrydb2(4,3)*wsum +&
                lqrydb(4)*dwsum(3) + lqrydb(3)*dwsum(4) +&
                lqryd(1)*dwsum2(4,3)
           qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +&
                lqrydb2(5,3)*wsum +&
                lqrydb(5)*dwsum(3) + lqrydb(3)*dwsum(5) +&
                lqryd(1)*dwsum2(5,3)
           qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +&
                lqrydb2(4,4)*wsum +&
                lqrydb(4)*dwsum(4) + lqrydb(4)*dwsum(4) +&
                lqryd(1)*dwsum2(4,4)
           qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +&
                lqrydb2(5,4)*wsum +&
                lqrydb(5)*dwsum(4) + lqrydb(4)*dwsum(5) +&
                lqryd(1)*dwsum2(5,4)
        endif
        ! lnoccupationa > lim_sum so may calculate approximation
        ! lnoccupationa evaluated for n = nmaxa which satisfies lnoccupation
        ! criterion for approximation (note that both hprime and gprime
        ! are always greater than unity.)
        ! also note that n-1 >= to nmaxa which satisfies "a" criterion
        ! call pi_totalsum_approx(ifpl, ifneutral, nmaxa,&
        !   a, b, nb, lqryd(1), lqryda(1), lqrydb,&
        !   lqryda2(1), lqrydab, lqrydb2)
        qryd(nmin_max) = qryd(nmin_max) + lqryd(1)*wapprox
        qryda(nmin_max) = qryda(nmin_max) + lqryda(1)*wapprox
        qryda2(nmin_max) = qryda2(nmin_max) + lqryda2(1)*wapprox
        qrydb(5,nmin_max) = qrydb(5,nmin_max) +&
             lqrydb(5)*wapprox + lqryd(1)*dwapprox(5)
        qrydab(5,nmin_max) = qrydab(5,nmin_max) +&
             lqrydab(5)*wapprox + lqryda(1)*dwapprox(5)
        qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) +&
             lqrydb2(5,5)*wapprox +&
             lqrydb(5)*dwapprox(5) + lqrydb(5)*dwapprox(5) +&
             lqryd(1)*dwapprox2(5,5)
        if(ifneutral) then
           qrydb(2,nmin_max) = qrydb(2,nmin_max) +&
                lqrydb(2)*wapprox + lqryd(1)*dwapprox(2)
           qrydb(3,nmin_max) = qrydb(3,nmin_max) +&
                lqrydb(3)*wapprox + lqryd(1)*dwapprox(3)
           qrydb(4,nmin_max) = qrydb(4,nmin_max) +&
                lqrydb(4)*wapprox + lqryd(1)*dwapprox(4)
           qrydab(2,nmin_max) = qrydab(2,nmin_max) +&
                lqrydab(2)*wapprox + lqryda(1)*dwapprox(2)
           qrydab(3,nmin_max) = qrydab(3,nmin_max) +&
                lqrydab(3)*wapprox + lqryda(1)*dwapprox(3)
           qrydab(4,nmin_max) = qrydab(4,nmin_max) +&
                lqrydab(4)*wapprox + lqryda(1)*dwapprox(4)
           qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +&
                lqrydb2(2,2)*wapprox +&
                lqrydb(2)*dwapprox(2) + lqrydb(2)*dwapprox(2) +&
                lqryd(1)*dwapprox2(2,2)
           qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +&
                lqrydb2(3,2)*wapprox +&
                lqrydb(3)*dwapprox(2) + lqrydb(2)*dwapprox(3) +&
                lqryd(1)*dwapprox2(3,2)
           qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +&
                lqrydb2(4,2)*wapprox +&
                lqrydb(4)*dwapprox(2) + lqrydb(2)*dwapprox(4) +&
                lqryd(1)*dwapprox2(4,2)
           qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +&
                lqrydb2(5,2)*wapprox +&
                lqrydb(5)*dwapprox(2) + lqrydb(2)*dwapprox(5) +&
                lqryd(1)*dwapprox2(5,2)
           qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +&
                lqrydb2(3,3)*wapprox +&
                lqrydb(3)*dwapprox(3) + lqrydb(3)*dwapprox(3) +&
                lqryd(1)*dwapprox2(3,3)
           qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +&
                lqrydb2(4,3)*wapprox +&
                lqrydb(4)*dwapprox(3) + lqrydb(3)*dwapprox(4) +&
                lqryd(1)*dwapprox2(4,3)
           qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +&
                lqrydb2(5,3)*wapprox +&
                lqrydb(5)*dwapprox(3) + lqrydb(3)*dwapprox(5) +&
                lqryd(1)*dwapprox2(5,3)
           qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +&
                lqrydb2(4,4)*wapprox +&
                lqrydb(4)*dwapprox(4) + lqrydb(4)*dwapprox(4) +&
                lqryd(1)*dwapprox2(4,4)
           qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +&
                lqrydb2(5,4)*wapprox +&
                lqrydb(5)*dwapprox(4) + lqrydb(4)*dwapprox(5) +&
                lqryd(1)*dwapprox2(5,4)
        endif
     else
        ! lnoccupationa > lim_approx so may calculate approximation
        ! lnoccupationa evaluated for n = nmaxa which satisfies lnoccupation
        ! criterion for approximation (note that both hprime and gprime
        ! are always greater than unity.)
        ! also note that n-1 >= to nmaxa which satisfies"a" criterion
        ! call pi_totalsum_approx(ifpl, ifneutral, nmaxa,&
        !   a, b, nb, lqryd(1), lqryda(1), lqrydb,&
        !   lqryda2(1), lqrydab, lqrydb2)
        !  use pure approximation
        qryd(nmin_max) = qryd(nmin_max) + lqryd(1)
        qryda(nmin_max) = qryda(nmin_max) + lqryda(1)
        qryda2(nmin_max) = qryda2(nmin_max) + lqryda2(1)
        qrydb(5,nmin_max) = qrydb(5,nmin_max) + lqrydb(5)
        qrydab(5,nmin_max) = qrydab(5,nmin_max) + lqrydab(5)
        qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max) + lqrydb2(5,5)
        if(ifneutral) then
           qrydb(2,nmin_max) = qrydb(2,nmin_max) + lqrydb(2)
           qrydb(3,nmin_max) = qrydb(3,nmin_max) + lqrydb(3)
           qrydb(4,nmin_max) = qrydb(4,nmin_max) + lqrydb(4)
           qrydab(2,nmin_max) = qrydab(2,nmin_max) + lqrydab(2)
           qrydab(3,nmin_max) = qrydab(3,nmin_max) + lqrydab(3)
           qrydab(4,nmin_max) = qrydab(4,nmin_max) + lqrydab(4)
           qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max) +&
                lqrydb2(2,2)
           qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max) +&
                lqrydb2(3,2)
           qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max) +&
                lqrydb2(4,2)
           qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max) +&
                lqrydb2(5,2)
           qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max) +&
                lqrydb2(3,3)
           qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max) +&
                lqrydb2(4,3)
           qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max) +&
                lqrydb2(5,3)
           qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max) +&
                lqrydb2(4,4)
           qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max) +&
                lqrydb2(5,4)
        endif
     endif
     ! calculate constant occupation probability factor and apply it.
     if(ifneutral) then
        expmb1 = exp(-b(1))
        dlnoccupation(1) = -1._fp_kind
     else
        expmb1 = 1._fp_kind
     endif
     qryd(nmin_max) = qryd(nmin_max)*expmb1
     qryda(nmin_max) = qryda(nmin_max)*expmb1
     qryda2(nmin_max) = qryda2(nmin_max)*expmb1
     qrydb(5,nmin_max) = qrydb(5,nmin_max)*expmb1
     qrydab(5,nmin_max) = qrydab(5,nmin_max)*expmb1
     qrydb2(5,5,nmin_max) = qrydb2(5,5,nmin_max)*expmb1
     if(ifneutral) then
        qrydb(1,nmin_max) = qryd(nmin_max)*dlnoccupation(1)
        qrydb(2,nmin_max) = qrydb(2,nmin_max)*expmb1
        qrydb(3,nmin_max) = qrydb(3,nmin_max)*expmb1
        qrydb(4,nmin_max) = qrydb(4,nmin_max)*expmb1
        qrydab(1,nmin_max) = qryda(nmin_max)*dlnoccupation(1)
        qrydab(2,nmin_max) = qrydab(2,nmin_max)*expmb1
        qrydab(3,nmin_max) = qrydab(3,nmin_max)*expmb1
        qrydab(4,nmin_max) = qrydab(4,nmin_max)*expmb1
        qrydb2(1,1,nmin_max) = qrydb(1,nmin_max)*dlnoccupation(1)
        qrydb2(2,1,nmin_max) = qrydb(2,nmin_max)*dlnoccupation(1)
        qrydb2(3,1,nmin_max) = qrydb(3,nmin_max)*dlnoccupation(1)
        qrydb2(4,1,nmin_max) = qrydb(4,nmin_max)*dlnoccupation(1)
        qrydb2(5,1,nmin_max) = qrydb(5,nmin_max)*dlnoccupation(1)
        qrydb2(2,2,nmin_max) = qrydb2(2,2,nmin_max)*expmb1
        qrydb2(3,2,nmin_max) = qrydb2(3,2,nmin_max)*expmb1
        qrydb2(4,2,nmin_max) = qrydb2(4,2,nmin_max)*expmb1
        qrydb2(5,2,nmin_max) = qrydb2(5,2,nmin_max)*expmb1
        qrydb2(3,3,nmin_max) = qrydb2(3,3,nmin_max)*expmb1
        qrydb2(4,3,nmin_max) = qrydb2(4,3,nmin_max)*expmb1
        qrydb2(5,3,nmin_max) = qrydb2(5,3,nmin_max)*expmb1
        qrydb2(4,4,nmin_max) = qrydb2(4,4,nmin_max)*expmb1
        qrydb2(5,4,nmin_max) = qrydb2(5,4,nmin_max)*expmb1
     endif
     ! convert from summand to sum for lower indices
     qryd(nmin_max-1:nmin:-1) = qryd(nmin_max-1:nmin:-1)*expmb1
     qryda(nmin_max-1:nmin:-1) = qryda(nmin_max-1:nmin:-1)*expmb1
     qryda2(nmin_max-1:nmin:-1) = qryda2(nmin_max-1:nmin:-1)*expmb1
     qrydb(5,nmin_max-1:nmin:-1) = qrydb(5,nmin_max-1:nmin:-1)*expmb1
     qrydab(5,nmin_max-1:nmin:-1) = qrydab(5,nmin_max-1:nmin:-1)*expmb1
     qrydb2(5,5,nmin_max-1:nmin:-1) = qrydb2(5,5,nmin_max-1:nmin:-1)*expmb1
     if(ifneutral) then
        qrydb(1,nmin_max-1:nmin:-1) = qryd(nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydb(2,nmin_max-1:nmin:-1) = qrydb(2,nmin_max-1:nmin:-1)*expmb1
        qrydb(3,nmin_max-1:nmin:-1) = qrydb(3,nmin_max-1:nmin:-1)*expmb1
        qrydb(4,nmin_max-1:nmin:-1) = qrydb(4,nmin_max-1:nmin:-1)*expmb1
        qrydab(1,nmin_max-1:nmin:-1) = qrydab(1,nmin_max:nmin+1:-1) + qryda(nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydab(2,nmin_max-1:nmin:-1) = qrydab(2,nmin_max:nmin+1:-1) + qrydab(2,nmin_max-1:nmin:-1)*expmb1
        qrydab(3,nmin_max-1:nmin:-1) = qrydab(3,nmin_max:nmin+1:-1) + qrydab(3,nmin_max-1:nmin:-1)*expmb1
        qrydab(4,nmin_max-1:nmin:-1) = qrydab(4,nmin_max:nmin+1:-1) + qrydab(4,nmin_max-1:nmin:-1)*expmb1
        qrydb2(1,1,nmin_max-1:nmin:-1) = qrydb2(1,1,nmin_max:nmin+1:-1) + qrydb(1,nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydb2(2,1,nmin_max-1:nmin:-1) = qrydb2(2,1,nmin_max:nmin+1:-1) + qrydb(2,nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydb2(3,1,nmin_max-1:nmin:-1) = qrydb2(3,1,nmin_max:nmin+1:-1) + qrydb(3,nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydb2(4,1,nmin_max-1:nmin:-1) = qrydb2(4,1,nmin_max:nmin+1:-1) + qrydb(4,nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydb2(5,1,nmin_max-1:nmin:-1) = qrydb2(5,1,nmin_max:nmin+1:-1) + qrydb(5,nmin_max-1:nmin:-1)*dlnoccupation(1)
        qrydb(1,nmin_max-1:nmin:-1) = qrydb(1,nmin_max:nmin+1:-1) + qrydb(1,nmin_max-1:nmin:-1)
        qrydb(2,nmin_max-1:nmin:-1) = qrydb(2,nmin_max:nmin+1:-1) + qrydb(2,nmin_max-1:nmin:-1)
        qrydb(3,nmin_max-1:nmin:-1) = qrydb(3,nmin_max:nmin+1:-1) + qrydb(3,nmin_max-1:nmin:-1)
        qrydb(4,nmin_max-1:nmin:-1) = qrydb(4,nmin_max:nmin+1:-1) + qrydb(4,nmin_max-1:nmin:-1)
        qrydb2(2,2,nmin_max-1:nmin:-1) = qrydb2(2,2,nmin_max:nmin+1:-1) + qrydb2(2,2,nmin_max-1:nmin:-1)*expmb1
        qrydb2(3,2,nmin_max-1:nmin:-1) = qrydb2(3,2,nmin_max:nmin+1:-1) + qrydb2(3,2,nmin_max-1:nmin:-1)*expmb1
        qrydb2(4,2,nmin_max-1:nmin:-1) = qrydb2(4,2,nmin_max:nmin+1:-1) + qrydb2(4,2,nmin_max-1:nmin:-1)*expmb1
        qrydb2(5,2,nmin_max-1:nmin:-1) = qrydb2(5,2,nmin_max:nmin+1:-1) + qrydb2(5,2,nmin_max-1:nmin:-1)*expmb1
        qrydb2(3,3,nmin_max-1:nmin:-1) = qrydb2(3,3,nmin_max:nmin+1:-1) + qrydb2(3,3,nmin_max-1:nmin:-1)*expmb1
        qrydb2(4,3,nmin_max-1:nmin:-1) = qrydb2(4,3,nmin_max:nmin+1:-1) + qrydb2(4,3,nmin_max-1:nmin:-1)*expmb1
        qrydb2(5,3,nmin_max-1:nmin:-1) = qrydb2(5,3,nmin_max:nmin+1:-1) + qrydb2(5,3,nmin_max-1:nmin:-1)*expmb1
        qrydb2(4,4,nmin_max-1:nmin:-1) = qrydb2(4,4,nmin_max:nmin+1:-1) + qrydb2(4,4,nmin_max-1:nmin:-1)*expmb1
        qrydb2(5,4,nmin_max-1:nmin:-1) = qrydb2(5,4,nmin_max:nmin+1:-1) + qrydb2(5,4,nmin_max-1:nmin:-1)*expmb1
     endif
     qryd(nmin_max-1:nmin:-1) = qryd(nmin_max:nmin+1:-1) + qryd(nmin_max-1:nmin:-1)
     qryda(nmin_max-1:nmin:-1) = qryda(nmin_max:nmin+1:-1) + qryda(nmin_max-1:nmin:-1)
     qryda2(nmin_max-1:nmin:-1) = qryda2(nmin_max:nmin+1:-1) + qryda2(nmin_max-1:nmin:-1)
     qrydb(5,nmin_max-1:nmin:-1) = qrydb(5,nmin_max:nmin+1:-1) + qrydb(5,nmin_max-1:nmin:-1)
     qrydab(5,nmin_max-1:nmin:-1) = qrydab(5,nmin_max:nmin+1:-1) + qrydab(5,nmin_max-1:nmin:-1)
     qrydb2(5,5,nmin_max-1:nmin:-1) = qrydb2(5,5,nmin_max:nmin+1:-1) + qrydb2(5,5,nmin_max-1:nmin:-1)
  endif
end subroutine qryd_approx
