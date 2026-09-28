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

!> This module provides variables that are used to help communicate
!> between the subset of the mod_pi module procedures that implement
!> the MDH pressure ionization scheme.
module mod_mdh_pi_data
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public mdh_pi_called, quad, pi_trace_state, pi_trace_extra
  real(fp_kind), save :: pi_trace_state(12)=0._fp_kind, pi_trace_extra(9,3)=0._fp_kind

  logical mdh_pi_called
  data mdh_pi_called /.false./

  real(fp_kind) quad

  save mdh_pi_called, quad
end module mod_mdh_pi_data

!> This module provides variables that are used to help communicate
!> between the subset of the mod_pi module procedures that implement
!> the PTEH pressure ionization scheme.
module mod_pteh_pi_data
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public pteh_pi_called, hion,&
       c1pi, c2pi, c3pi, c4pi, c5pi, c6pi,&
       onethird, onethirdm1, onethirdm2,&
       da, dax, day, daf, dat, daxf, daxt, dayf, dayt, damday, dafmdayf, datmdayt

  ! This is the exact value used by PTEH for H2 ionization in K^-1
  ! this is a fairly inaccurate value, but we use it for
  ! our own PTEH pressure ionization implementation to help mimic their
  ! implementation as closely as possible.
  real(fp_kind), parameter :: hion=11605.0_fp_kind*13.6_fp_kind

  real(fp_kind), parameter :: c5pi = 1.07654_fp_kind
  real(fp_kind), parameter :: c6pi = 0.61315_fp_kind

  real(fp_kind), parameter :: onethird = 1._fp_kind/3._fp_kind
  real(fp_kind), parameter :: onethirdm1 = onethird-1._fp_kind
  real(fp_kind), parameter :: onethirdm2 = onethird-2._fp_kind

  logical pteh_pi_called
  data pteh_pi_called /.false./

  real(fp_kind) c1pi, c2pi, c3pi, c4pi

  real(fp_kind) da, dax, day, daf, dat, daxf, daxt, dayf, dayt, damday, dafmdayf, datmdayt

  save pteh_pi_called,&
       c1pi, c2pi, c3pi, c4pi,&
       da, dax, day, daf, dat, daxf, daxt, dayf, dayt, damday, dafmdayf, datmdayt
end module mod_pteh_pi_data

!> This module provides variables that are used to help communicate
!> between the subset of the mod_pi module procedures that implement
!> the FJS pressure ionization scheme.
module mod_fjs_pi_data
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public&
       ! Parameters used by subroutine fjs_pi and other fjs-related routines
       maxfjs_aux,&

       ! Parameters used just by subroutine fjs_pi
       api1_orig, api2_orig, api3_orig, api4_orig,&
       chipi1_orig, chipi2_orig, chipi3_orig, chipi4_orig,&
       api1_mod, api2_mod, api3_mod, api4_mod,&
       chipi1_mod, chipi2_mod, chipi3_mod, chipi4_mod,&
       fudge_ln,&

       ! Variables calculated by subroutine fjs_pi and saved by that routine (for
       ! the case when ifmodified changes) and potentially for other fjs-related
       ! routines.
       ifmodified_old,&
       api1, api2, api3, api4,&
       chipi1, chipi2, chipi3, chipi4,&
       omht,&
       omhe1t,&
       omhe2t,&
       omet,&
       chpi1,&
       chpi2,&
       chpi3,&
       chpi4,&

       ! Variables calculated by subroutine fjs_pi and saved by that routine (strictly
       ! for use by other fjs-related routines).
       fjs_pi_called,&
       omh,&
       omhe1,&
       omhe2,&
       ome

  ! aux(1:maxfjs_aux) are the aux variables that are relevant to fjs
  ! pressure ionization.  N.B. h2plus = aux(5) may affect the
  ! if_mc.eq.1 approximations for the Coulomb sums that are often used
  ! in combination with fjs pressure ionization, but that auxiliary
  ! variable is not counted here since it does not affect fjs pressure
  ! ionization.
  integer, parameter :: maxfjs_aux = 4

  real(fp_kind), parameter::&
       ! original fjs set.
       api1_orig = 1.2e-8_fp_kind,&
       api2_orig = 0.45e-8_fp_kind,&
       api3_orig = 0.3e-8_fp_kind,&
       api4_orig = 0.57e-8_fp_kind,&
       chipi1_orig = 0.5_fp_kind,&
       chipi2_orig = 36.0_fp_kind,&
       chipi3_orig = 4.0_fp_kind,&
       chipi4_orig = 8.0_fp_kind,&
       ! working set = api*_org + fudge that gives reasonable fit to EOS1 for EOS2 case.
       api1_mod = 1.2e-8_fp_kind,&
       api2_mod = 0.45e-8_fp_kind,&
       api3_mod = 0.3e-8_fp_kind,&
       api4_mod = 0.57e-8_fp_kind,&
       ! Fixme: I originally got peculiar fit result that the best
       ! fits of EOS2 to EOS1 were for negative values of these parameters (which is
       ! contrary to the idea these have something to do with
       ! ionization potentials).  So I derived the best rough fit I
       ! could with these fixed at zero.  However, may want to try
       ! something non-zero for the next fit in conjunction with
       ! different fudge factors or give up all together with regard
       ! to physical meaning and simply adopt negative values of these
       ! that give the best fit.
       chipi1_mod = 0._fp_kind,&
       chipi2_mod = 0._fp_kind,&
       chipi3_mod = 0._fp_kind,&
       chipi4_mod = 0._fp_kind,&
       fudge_ln(4) = [&
       ! n.b. (3) must be sufficiently larger than (2) otherwise, He+ is
       ! never pressure ionized to He++.
       !old  data fudge_ln /1.d0, 4.5d0, 5.5d0, -15.0d0/
       ! small adjustment made to fit pure He opal tables.  Didn't degrade
       ! rest of table fit or extension fit.
       ! data fudge_ln /1.d0, 5.0d0, 5.5d0, -5.0d0/
       ! now adjust again against best fit to EOS1 rather than opal especially for
       ! zero chipi*_mod case
       1._fp_kind, 2.5_fp_kind, 3.0_fp_kind, -5.0_fp_kind&
       ]

  ! Saved variables calculated by fjs_pi that are used by that routine
  ! and the other fjs-related routines.

  logical fjs_pi_called
  data fjs_pi_called /.false./

  integer ifmodified_old
  ! Initialize to invalid value to start.
  data ifmodified_old/-1000/

  real(fp_kind)&
       api1, api2, api3, api4,&
       chipi1, chipi2, chipi3, chipi4,&
       omht,&
       omhe1t,&
       omhe2t,&
       omet,&
       chpi1,&
       chpi2,&
       chpi3,&
       chpi4,&
       omh,&
       omhe1,&
       omhe2,&
       ome

  save&
       fjs_pi_called,&
       ifmodified_old,&
       api1, api2, api3, api4,&
       chipi1, chipi2, chipi3, chipi4,&
       omht,&
       omhe1t,&
       omhe2t,&
       omet,&
       chpi1,&
       chpi2,&
       chpi3,&
       chpi4,&
       omh,&
       omhe1,&
       omhe2,&
       ome

end module mod_fjs_pi_data

!> This module provides public module procedures that implement the
!> chosen non-ideal pressure-ionization component of the free-energy
!> and all relevant derivatives of same.  One of the MDH, PTEH, or FJS
!> approximations for pressure ionization are chosen for this
!> component depending on option suite.
module mod_pi
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public&
       mdh_pi, mdh_pi_pressure_free, mdh_pi_end,&
       fjs_pi, fjs_pi_free, fjs_pi_end,&
       pteh_pi, pteh_pi_pressure, pteh_pi_free, pteh_pi_end
contains
#include "fjs_pi.f90"
#include "mdh_pi.f90"
#include "pteh_pi.f90"
end module mod_pi
