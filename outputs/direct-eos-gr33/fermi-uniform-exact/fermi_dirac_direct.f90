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

! Calculate Fermi-Dirac integral (C.G. eq. 24,97) (answer(1)) and
! its partial derivative with respect to e=eta (answer(2)),
! b=beta (answer(3)), ee (answer(4)), eb (answer(5)), bb (answer(6)),
! eee (answer(7)), eeb (answer(8)), and ebb (answer(9)).
! N.B.
! only answer(1) through answer(min(1, max(9,ifderivative)))
! is returned.

!> This fermi_dirac_direct subroutine calculates the Fermi-Dirac
!> integral defined by C.G. eq. 24.97 and its partial derivatives
!> using numerical quadrature.
!>
!> \param[in] akindex PARAMETERS NEED DOCUMENTATION
!>
subroutine fermi_dirac_direct(akindex, eta, beta, answer)

  !  Follow
  !  <https://www.fortran90.org/src/best-practices.html#ii-general-structure>
  !  for implementing the required callback functions and their data and
  !  context arguments and implementing the (integration) routine that calls that callback function.
  use mod_fermi_dirac_integrator, only: fermi_dirac_context, fermi_dirac_integrator

  ! Arguments
  real(fp_kind), intent(in) :: akindex, eta, beta
  real(fp_kind), intent(out) :: answer(:)

  ! Internal variables

  ! Requested accuracy of integration.
  real(fp_kind), parameter :: fderr = 1.e-12_fp_kind

  integer nanswer

  real(fp_kind) amp, ampmid
  real(fp_kind) anslow, ansmid
  type(fermi_dirac_context) :: fddata

  nanswer = size(answer)

  ! Sanity check
  if(nanswer.lt.1.or.nanswer.gt.9) error stop 'fermi_dirac_direct: invalid size of answer'

  ! Prepare context data to be passed (via fermi_dirac_integrator) to callback routines.
  fddata%akindex = akindex
  fddata%eta = eta
  fddata%beta = beta

  ! range of integration.
  ! 1/exp(40) ~ 4.d-18
  amp= 40._fp_kind+abs(eta)
  if(eta.le.0._fp_kind.or.amp-80._fp_kind.le.0._fp_kind) then
     ! divide integral into two ranges:
     ! first range where exp factor is variable,
     ! but not overwhelmingly so, and upper range where exponential
     ! cutoff factor is most active.
     if(eta.le.0._fp_kind) then
        ! Most of the integral determined below this value.
        ampmid = 1._fp_kind
     else
        ! N.B. eta > 0 ==> amp = eta + 40 <= 80 ==> eta <= 40 ==> amp - 80 < eta
        ampmid = 1._fp_kind + eta
     endif
     call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,fermi_fun,fddata,fderr)
     call fermi_dirac_integrator(answer(1),ampmid,amp,fermi_fun,fddata,fderr)
     answer(1) = ansmid + answer(1)
     if(nanswer.ge.2) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_de,fddata,fderr)
        call fermi_dirac_integrator(answer(2),ampmid,amp,dfermi_fun_de,fddata,fderr)
        answer(2) = ansmid + answer(2)
     endif
     if(nanswer.ge.3) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_db,fddata,fderr)
        call fermi_dirac_integrator(answer(3),ampmid,amp,dfermi_fun_db,fddata,fderr)
        answer(3) = ansmid + answer(3)
     endif
     if(nanswer.ge.4) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_dee,fddata,fderr)
        call fermi_dirac_integrator(answer(4),ampmid,amp,dfermi_fun_dee,fddata,fderr)
        answer(4) = ansmid + answer(4)
     endif
     if(nanswer.ge.5) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_deb,fddata,fderr)
        call fermi_dirac_integrator(answer(5),ampmid,amp,dfermi_fun_deb,fddata,fderr)
        answer(5) = ansmid + answer(5)
     endif
     if(nanswer.ge.6) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_dbb,fddata,fderr)
        call fermi_dirac_integrator(answer(6),ampmid,amp,dfermi_fun_dbb,fddata,fderr)
        answer(6) = ansmid + answer(6)
     endif
     if(nanswer.ge.7) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_deee,fddata,fderr)
        call fermi_dirac_integrator(answer(7),ampmid,amp,dfermi_fun_deee,fddata,fderr)
        answer(7) = ansmid + answer(7)
     endif
     if(nanswer.ge.8) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_deeb,fddata,fderr)
        call fermi_dirac_integrator(answer(8),ampmid,amp,dfermi_fun_deeb,fddata,fderr)
        answer(8) = ansmid + answer(8)
     endif
     if(nanswer.ge.9) then
        call fermi_dirac_integrator(ansmid,0._fp_kind,ampmid,dfermi_fun_debb,fddata,fderr)
        call fermi_dirac_integrator(answer(9),ampmid,amp,dfermi_fun_debb,fddata,fderr)
        answer(9) = ansmid + answer(9)
     endif
  else
     ! divide integral into three ranges: lower range where cutoff factor
     ! 1/(exp(x-eta)+1) is unity, middle range where exp factor is variable,
     ! but not overwhelmingly so,and upper range where exponential
     ! cutoff factor is most active.
     ! N.B. eta > 0 ==> amp = eta + 40 > 80 ==> eta > 40 ==> amp - 80 < eta
     call fermi_dirac_integrator(anslow,0._fp_kind,amp-80._fp_kind,fermi_fun0,fddata,fderr)
     call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,fermi_fun,fddata,fderr)
     call fermi_dirac_integrator(answer(1),eta,amp,fermi_fun,fddata,fderr)
     answer(1) = anslow + ansmid + answer(1)
     if(nanswer.ge.2) then
        ! note: lower range integral is zero so just have to calculate
        ! middle and upper range.
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_de,fddata,fderr)
        call fermi_dirac_integrator(answer(2),eta,amp,dfermi_fun_de,fddata,fderr)
        answer(2) = ansmid + answer(2)
     endif
     if(nanswer.ge.3) then
        call fermi_dirac_integrator(anslow,0._fp_kind,amp-80._fp_kind,dfermi_fun0_db,fddata,fderr)
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_db,fddata,fderr)
        call fermi_dirac_integrator(answer(3),eta,amp,dfermi_fun_db,fddata,fderr)
        answer(3) = anslow + ansmid + answer(3)
     endif
     if(nanswer.ge.4) then
        ! note: lower range integral is zero so just have to calculate
        ! middle and upper range.
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_dee,fddata,fderr)
        call fermi_dirac_integrator(answer(4),eta,amp,dfermi_fun_dee,fddata,fderr)
        answer(4) = ansmid + answer(4)
     endif
     if(nanswer.ge.5) then
        ! note: lower range integral is zero so just have to calculate
        ! middle and upper range.
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_deb,fddata,fderr)
        call fermi_dirac_integrator(answer(5),eta,amp,dfermi_fun_deb,fddata,fderr)
        answer(5) = ansmid + answer(5)
     endif
     if(nanswer.ge.6) then
        call fermi_dirac_integrator(anslow,0._fp_kind,amp-80._fp_kind,dfermi_fun0_dbb,fddata,fderr)
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_dbb,fddata,fderr)
        call fermi_dirac_integrator(answer(6),eta,amp,dfermi_fun_dbb,fddata,fderr)
        answer(6) = anslow + ansmid + answer(6)
     endif
     if(nanswer.ge.7) then
        ! note: lower range integral is zero so just have to calculate
        ! middle and upper range.
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_deee,fddata,fderr)
        call fermi_dirac_integrator(answer(7),eta,amp,dfermi_fun_deee,fddata,fderr)
        answer(7) = ansmid + answer(7)
     endif
     if(nanswer.ge.8) then
        ! note: lower range integral is zero so just have to calculate
        ! middle and upper range.
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_deeb,fddata,fderr)
        call fermi_dirac_integrator(answer(8),eta,amp,dfermi_fun_deeb,fddata,fderr)
        answer(8) = ansmid + answer(8)
     endif
     if(nanswer.ge.9) then
        ! note: lower range integral is zero so just have to calculate
        ! middle and upper range.
        call fermi_dirac_integrator(ansmid,amp-80._fp_kind,eta,dfermi_fun_debb,fddata,fderr)
        call fermi_dirac_integrator(answer(9),eta,amp,dfermi_fun_debb,fddata,fderr)
        answer(9) = ansmid + answer(9)
     endif
  endif
end subroutine fermi_dirac_direct

!> This fermi_fun function calculates the integrand of the Fermi-Dirac
!> integral defined by C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function fermi_fun(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) fermi_fun
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  fermi_fun=x**akindex
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  a2=exp(x-eta) + 1._fp_kind
  fermi_fun=fermi_fun*a1/a2
end function fermi_fun

!> This dfermi_fun_de function calculates the partial derivative wrt
!> eta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_de(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_de
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2, a3, a4
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_de = x**akindex
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  ! 1.d300 is roughly exp(690).  So guard against under and overflow.
  if(x-eta.le.-690._fp_kind) then
     ! underflow. a2 ~ 1, a4 gets very small
     dfermi_fun_de = 0._fp_kind
  elseif(x-eta.le.690._fp_kind) then
     a3 = exp(x-eta)
     a2 = a3 + 1._fp_kind
     a4 = a3/a2
     dfermi_fun_de = dfermi_fun_de*a1*a4/a2
  else
     ! overflow.  a4 ~ 1, a2 gets very large
     dfermi_fun_de = 0._fp_kind
  endif
end function dfermi_fun_de

!> This dfermi_fun_db function calculates the partial derivative wrt
!> beta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_db(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Functions and arguments
  real(fp_kind) dfermi_fun_db
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_db = 0.25_fp_kind*x**(akindex+1._fp_kind)
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  a2=exp(x-eta) + 1._fp_kind
  dfermi_fun_db = dfermi_fun_db/(a1*a2)
end function dfermi_fun_db

!> This dfermi_fun_dee function calculates the second partial derivative wrt
!> eta, eta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_dee(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_dee
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2, a3, a4
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_dee = x**akindex
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  ! 1.d300 is roughly exp(690).  So guard against under and overflow.
  if(x-eta.le.-690._fp_kind) then
     ! underflow. a2 ~ 1, a4 gets very small
     dfermi_fun_dee = 0._fp_kind
  elseif(x-eta.le.690._fp_kind) then
     a3 = exp(x-eta)
     a2 = a3 + 1._fp_kind
     a4 = a3/a2
     dfermi_fun_dee = dfermi_fun_dee*a1*&
          (2._fp_kind*a4-1._fp_kind)*a4/a2
  else
     ! overflow.  a4 ~ 1, a2 gets very large
     dfermi_fun_dee = 0._fp_kind
  endif
end function dfermi_fun_dee

!> This dfermi_fun_deb function calculates the second partial derivative wrt
!> eta, beta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_deb(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_deb
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2, a3, a4
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_deb = 0.25_fp_kind*x**(akindex+1._fp_kind)
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  !      1.d300 is roughly exp(690).  So guard against under and overflow.
  if(x-eta.le.-690._fp_kind) then
     ! underflow. a2 ~ 1, a4 gets very small
     dfermi_fun_deb = 0._fp_kind
  elseif(x-eta.le.690._fp_kind) then
     a3 = exp(x-eta)
     a2 = a3 + 1._fp_kind
     a4 = a3/a2
     dfermi_fun_deb = dfermi_fun_deb*a4/(a1*a2)
  else
     ! overflow.  a4 ~ 1, a2 gets very large
     dfermi_fun_deb = 0._fp_kind
  endif
end function dfermi_fun_deb

!> This dfermi_fun_dbb function calculates the second partial derivative wrt
!> beta, beta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_dbb(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_dbb
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_dbb = -0.0625_fp_kind*x**(akindex+2._fp_kind)
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)*(1._fp_kind+0.5_fp_kind*beta*x)
  a2=exp(x-eta) + 1._fp_kind
  dfermi_fun_dbb = dfermi_fun_dbb/(a1*a2)
end function dfermi_fun_dbb

!> This dfermi_fun_deee function calculates the third partial derivative wrt
!> eta, eta, eta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_deee(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_deee
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2, a3, a4
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_deee = x**akindex
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  ! 1.d300 is roughly exp(690).  So guard against under and overflow.
  if(x-eta.le.-690._fp_kind) then
     ! underflow. a2 ~ 1, a4 gets very small
     dfermi_fun_deee = 0._fp_kind
  elseif(x-eta.le.690._fp_kind) then
     a3 = exp(x-eta)
     a2 = a3 + 1._fp_kind
     a4 = a3/a2
     ! 1 - a4 = (a2 - a3)/a2 = 1/a2
     dfermi_fun_deee = dfermi_fun_deee*a1*&
          (1._fp_kind - 6._fp_kind*a4/a2)*a4/a2
  else
     ! overflow.  a4 ~ 1, a2 gets very large
     dfermi_fun_deee = 0._fp_kind
  endif
end function dfermi_fun_deee

!> This dfermi_fun_deeb function calculates the third partial derivative wrt
!> eta, eta, beta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_deeb(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_deeb
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2, a3, a4
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_deeb = 0.25_fp_kind*x**(akindex+1._fp_kind)
  a1=sqrt(1._fp_kind+0.5_fp_kind*beta*x)
  ! 1.d300 is roughly exp(690).  So guard against under and overflow.
  if(x-eta.le.-690._fp_kind) then
     ! underflow. a2 ~ 1, a4 gets very small
     dfermi_fun_deeb = 0._fp_kind
  elseif(x-eta.le.690._fp_kind) then
     a3 = exp(x-eta)
     a2 = a3 + 1._fp_kind
     a4 = a3/a2
     dfermi_fun_deeb = dfermi_fun_deeb*&
          (2._fp_kind*a4-1._fp_kind)*a4/(a1*a2)
  else
     ! overflow.  a4 ~ 1, a2 gets very large
     dfermi_fun_deeb = 0._fp_kind
  endif
end function dfermi_fun_deeb

!> This dfermi_fun_debb function calculates the third partial derivative wrt
!> eta, beta, beta of the integrand of the Fermi-Dirac integral defined by
!> C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun_debb(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun_debb
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) a1, a2, a3, a4
  real(fp_kind) akindex, eta, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  eta = fddata%eta
  beta = fddata%beta

  dfermi_fun_debb = -0.0625_fp_kind*x**(akindex+2._fp_kind)
  a1 = sqrt(1._fp_kind+0.5_fp_kind*beta*x)*(1._fp_kind+0.5_fp_kind*beta*x)
  a3 = exp(x-eta)
  a2 = a3  + 1._fp_kind
  a4 = a3/a2
  dfermi_fun_debb = dfermi_fun_debb*a4/(a1*a2)
end function dfermi_fun_debb

!> This fermi_fun0 function calculates (in the x-eta >> 1 limit where
!> denominator = 1) the integrand of the Fermi-Dirac integral defined
!> by C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function fermi_fun0(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) fermi_fun0
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) akindex, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  beta = fddata%beta

  fermi_fun0=(x**akindex)*sqrt(1._fp_kind+0.5_fp_kind*beta*x)
end function fermi_fun0

!> This dfermi_fun0_db function calculates (in the x-eta >> 1 limit
!> where denominator = 1) the partial derivative wrt beta of the
!> integrand of the Fermi-Dirac integral defined by C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun0_db(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun0_db
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) akindex, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  beta = fddata%beta

  dfermi_fun0_db = (0.25_fp_kind*x**(akindex+1._fp_kind))/&
       sqrt(1._fp_kind+0.5_fp_kind*beta*x)
end function dfermi_fun0_db

!> This dfermi_fun0_dbb function calculates (in the x-eta >> 1 limit
!> where denominator = 1) the second partial derivative wrt beta, beta of the
!> integrand of the Fermi-Dirac integral defined by C.G. eq. 24.97.
!>
!> \param[in] x PARAMETERS NEED DOCUMENTATION
!>
function dfermi_fun0_dbb(x, fddata)
  use mod_fermi_dirac_integrator, only: fermi_dirac_context

  ! Function and arguments
  real(fp_kind) dfermi_fun0_dbb
  real(fp_kind), intent(in) :: x
  type(fermi_dirac_context), intent(inout) :: fddata

  ! Internal variables
  real(fp_kind) akindex, beta

  ! Prepare local context data that has been passed (via
  ! fermi_dirac_integrator and fddata) to this callback routine.
  akindex = fddata%akindex
  beta = fddata%beta

  dfermi_fun0_dbb = (-0.0625_fp_kind*x**(akindex+2._fp_kind))/&
       (sqrt(1._fp_kind+0.5_fp_kind*beta*x)*(1._fp_kind+0.5_fp_kind*beta*x))
end function dfermi_fun0_dbb
