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

! calculate planck-larkin sum where
! sum(nmin) = sum from nmin to infinity (nmin = nlow,...,nhigh) of
! 2 n^2 (exp(a/n^2) -(1+a/n^2))
! dsum is the derivative of sum wrt to a
! dsum2 is the second derivative of sum wrt to a
! use direct summation until can use plsum_approx
! if nlow too large for plsum_approx, use scaling rule

!> This plsum subroutine calculates (using a combination of explicit
!> summation and an approximation) the Planck-Larkin sums and the
!> first and second derivatives of those quantities wrt to the
!> independent variable a (a convenient transformation of t).  These
!> sums are defined as the sum from n = nmin (where nmin takes on the
!> range of values from nlow to nigh) to n = infinity of 2 n^2
!> (exp(a/n^2) -(1+a/n^2)).
!>
!> \param[in] nlow PARAMETERS NEED DOCUMENTATION
!>
subroutine plsum(nlow, nhigh, ioffset, a, sum, dsum, dsum2)

  ! arguments
  integer, intent(in) :: nlow, nhigh, ioffset
  real(fp_kind), intent(in) :: a
  real(fp_kind), intent(out) :: sum(nhigh-ioffset), dsum(nhigh-ioffset), dsum2(nhigh-ioffset)

  ! internal variables

  ! at this value or greater use scaling law
  ! have proved through tests that the scaling error goes like
  ! nlow_scale^3 and the relative error at this value is a maximum
  ! of 2.d-5
  integer, parameter :: nlow_scale=31
  ! this is a standard value for most of the approximations
  ! which is used to limit their range of applicability
  real(fp_kind), parameter :: exp_max = 1._fp_kind

  integer n, nstore, nstorep1, nstart

  real(fp_kind) summand, summanda, summanda2
  real(fp_kind) arg, rn2, savearg, scale

  ! Sanity checks
  if(nlow-ioffset.lt.1) error stop 'plsum: bad ioffset value'
  if(a.le.0._fp_kind) error stop 'plsum: bad a value'
  if(nlow.le.0.or.nlow.gt.nhigh) error stop 'plsum: bad nlow or nhigh'

  ! full sum results via scaling or by approximation.
  ! minimum principal quantum number that can be done this way.
  nstart = max(2,nlow,int(sqrt(a/exp_max))+1)
  n = nstart
  do while(n.eq.nstart.or.n.le.nhigh)
     ! first part of while assures sum calculated for n = nstart,
     ! and next statement assures this result stored in nhigh, if
     ! this do loop only done once because nstart > nhigh.
     nstore = min(n,nhigh)-ioffset
     if(n.ge.nlow_scale) then
        ! scaling rule.
        scale = real((n),fp_kind)/real((nlow_scale-1),fp_kind)
        call plsum_approx(nlow_scale-1, a/scale/scale,&
             sum(nstore), dsum(nstore), dsum2(nstore))
        dsum(nstore) = dsum(nstore)/(scale*scale)
        dsum2(nstore) = dsum2(nstore)/(scale*scale*scale*scale)
        ! transform results to integral from nlow_scale -1 of
        ! f(a/scale/scale,t) according to trapezoidal rule,
        ! transform to integral from n of f(a,t) by multiplying by scale^3,
        ! and finally transform to sum by trapezoidal rule
        ! n.b. argprime = a/scale/scale/(n_low_scale-1)^2 = arg
        arg = a/(real((n),fp_kind)*real((n),fp_kind))
        call plsummand(arg, summand, summanda, summanda2)
        savearg = real((n),fp_kind)*real((n),fp_kind)
        ! sum(nstore) = (sum(nstore) - savearg/scale/scale*(exparg-(1.d0+arg)))*&
        !   scale*scale*scale + savearg*(exparg-(1.d0+arg))
        ! dsum(nstore) = (dsum(nstore) - 1.d0/scale/scale*(exparg-1.d0))*&
        !   scale*scale*scale + (exparg-1.d0)
        ! dsum2(nstore) = (dsum2(nstore) - 1.d0/scale/scale/savearg*exparg)*&
        !   scale*scale*scale + exparg
        sum(nstore) = sum(nstore)*scale*scale*scale +&
             savearg*summand*(1._fp_kind-scale)
        dsum(nstore) = dsum(nstore)*scale*scale*scale +&
             1._fp_kind*summanda*(1._fp_kind-scale)
        dsum2(nstore) = dsum2(nstore)*scale*scale*scale +&
             1._fp_kind/savearg*summanda2*(1._fp_kind-scale)
     else
        call plsum_approx(n, a, sum(nstore), dsum(nstore), dsum2(nstore))
     endif
     n = n + 1
  enddo
  ! at this point have calculated full sum results for indexes from
  ! nstart to nhigh (n.b. or just nstart and stored in nhigh index
  ! if nhigh smaller than nstart).
  ! do explicit summation for remaining minimum principal
  ! quantum numbers
  do n = nstart-1,nlow,-1
     rn2 = real((n),fp_kind)*real((n),fp_kind)
     arg = a/rn2
     call plsummand_normalized(rn2, arg, summand, summanda, summanda2)
     nstore = min(nhigh,n)-ioffset
     nstorep1 = min(nhigh,n+1)-ioffset
     sum(nstore) = sum(nstorep1) + summand
     dsum(nstore) = dsum(nstorep1) + summanda
     dsum2(nstore) = dsum2(nstorep1) + summanda2
  enddo
end subroutine plsum

! Calculate Planck-Larkin summand factor
! exp(arg) - (1. + arg)
! and its derivatives with respect to arg.

!> This plsummand subroutine calculates the Planck-Larkin summand
!> factor, exp(arg) - (1. + arg) and its first two derivatives wrt
!> arg.
!>
!> \param[in] arg PARAMETERS NEED DOCUMENTATION
!>

subroutine plsummand(arg, summand, summanda, summanda2)

  ! arguments
  real(fp_kind), intent(in) :: arg
  real(fp_kind), intent(out) :: summand, summanda, summanda2

  ! internal variables
  real(fp_kind) exparg

  if(arg.gt.0.01_fp_kind) then
     ! lose a maximum of 4 significant digits
     exparg = exp(arg)
     summand = exparg - (1._fp_kind + arg)
     summanda = exparg - 1._fp_kind
     summanda2 = exparg
  else
     summand = arg*arg/2._fp_kind*(&
          1._fp_kind + arg/3._fp_kind*(&
          1._fp_kind + arg/4._fp_kind*(&
          1._fp_kind + arg/5._fp_kind*(&
          1._fp_kind + arg/6._fp_kind*(&
          1._fp_kind + arg/7._fp_kind*(&
          1._fp_kind + arg/8._fp_kind*(&
          1._fp_kind + arg/9._fp_kind*(&
          1._fp_kind + arg/10._fp_kind*(&
          1._fp_kind + arg/11._fp_kind)))))))))
     summanda = arg*(&
          1._fp_kind + arg/2._fp_kind*(&
          1._fp_kind + arg/3._fp_kind*(&
          1._fp_kind + arg/4._fp_kind*(&
          1._fp_kind + arg/5._fp_kind*(&
          1._fp_kind + arg/6._fp_kind*(&
          1._fp_kind + arg/7._fp_kind*(&
          1._fp_kind + arg/8._fp_kind*(&
          1._fp_kind + arg/9._fp_kind*(&
          1._fp_kind + arg/10._fp_kind)))))))))
     summanda2 = (&
          1._fp_kind + arg*(&
          1._fp_kind + arg/2._fp_kind*(&
          1._fp_kind + arg/3._fp_kind*(&
          1._fp_kind + arg/4._fp_kind*(&
          1._fp_kind + arg/5._fp_kind*(&
          1._fp_kind + arg/6._fp_kind*(&
          1._fp_kind + arg/7._fp_kind*(&
          1._fp_kind + arg/8._fp_kind*(&
          1._fp_kind + arg/9._fp_kind)))))))))
  endif
end subroutine plsummand

!> This plsummand_normalized subroutine calculates the normalized Planck-Larkin summand
!> factor, 2.*rn2*(exp(arg) - (1. + arg)) and its first two derivatives wrt
!> a where arg = a/rn2.
!>
!> \param[in] arg PARAMETERS NEED DOCUMENTATION
!>

subroutine plsummand_normalized(rn2, arg, summand, summanda, summanda2)

  ! arguments
  real(fp_kind), intent(in) :: rn2, arg
  real(fp_kind), intent(out) :: summand, summanda, summanda2

  ! internal variables
  real(fp_kind) exparg

  if(arg.gt.0.01_fp_kind) then
     ! lose a maximum of 4 significant digits
     exparg = exp(arg)
     summand = 2._fp_kind*rn2*(exparg - (1._fp_kind + arg))
     summanda = 2._fp_kind*(exparg - 1._fp_kind)
     summanda2 = 2._fp_kind*exparg/rn2
  else
     summand = rn2*arg*arg*(&
          1._fp_kind + arg/3._fp_kind*(&
          1._fp_kind + arg/4._fp_kind*(&
          1._fp_kind + arg/5._fp_kind*(&
          1._fp_kind + arg/6._fp_kind*(&
          1._fp_kind + arg/7._fp_kind*(&
          1._fp_kind + arg/8._fp_kind*(&
          1._fp_kind + arg/9._fp_kind*(&
          1._fp_kind + arg/10._fp_kind*(&
          1._fp_kind + arg/11._fp_kind)))))))))
     summanda = 2._fp_kind*arg*(&
          1._fp_kind + arg/2._fp_kind*(&
          1._fp_kind + arg/3._fp_kind*(&
          1._fp_kind + arg/4._fp_kind*(&
          1._fp_kind + arg/5._fp_kind*(&
          1._fp_kind + arg/6._fp_kind*(&
          1._fp_kind + arg/7._fp_kind*(&
          1._fp_kind + arg/8._fp_kind*(&
          1._fp_kind + arg/9._fp_kind*(&
          1._fp_kind + arg/10._fp_kind)))))))))
     summanda2 = 2._fp_kind*(&
          1._fp_kind + arg*(&
          1._fp_kind + arg/2._fp_kind*(&
          1._fp_kind + arg/3._fp_kind*(&
          1._fp_kind + arg/4._fp_kind*(&
          1._fp_kind + arg/5._fp_kind*(&
          1._fp_kind + arg/6._fp_kind*(&
          1._fp_kind + arg/7._fp_kind*(&
          1._fp_kind + arg/8._fp_kind*(&
          1._fp_kind + arg/9._fp_kind)))))))))/rn2
  endif
end subroutine plsummand_normalized
