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

!> This module provides the public molecular_hydrogen module procedure that calculates
!> ideal internal partition functions of H_2 and H_2+ as a function of temperature
!> according to a number of different approximations for these partition
!> functions that are available to the user.
!>
module mod_molecular_hydrogen
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public molecular_hydrogen
contains

  ! calculate partition function eos quantities for molecular hydogen and
  ! its positive ion.
  ! ifh2 = 0  no h2
  ! ifh2 = 1  vdb h2
  ! ifh2 = 2  st h2
  ! ifh2 = 3  irwin h2 (recommended)
  ! ifh2 = 4  pteh h2
  ! ifh2plus = 0  no h2plus
  ! ifh2plus = 1  st h2plus
  ! ifh2plus = 2  irwin h2plus (recommended)
  ! qh2 is the ln partition function
  ! qh2t is dlnq/dlnt
  ! qh2tt is d2lnq/(dlnt)^2
  ! qh2plus is the ln partition function
  ! qh2plust is dlnq/dlnt
  ! qh2plustt is d2lnq/(dlnt)^2

  !> This molecular_hydrogen subroutine calculates ideal internal
  !> partition functions of H_2 and H_2+ and their first and second
  !> derivatives wrt temperatureas a function of temperature according
  !> to a number of different approximations for these partition
  !> functions that have been implemented.
  !>
  !> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
  !>
  subroutine molecular_hydrogen(verbosity, ifh2, ifh2plus, tl,&
       qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt)
    use mod_free_eos_constants, only: boltzmann, ergsperev, ln10
    use mod_poly_sum, only: poly_sum
    use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

    ! Arguments
    integer, intent(in) :: verbosity, ifh2, ifh2plus
    real(fp_kind), intent(in) :: tl
    real(fp_kind), intent(out) :: qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt

    ! Internal variables
    ! partition function coefficients
    real(fp_kind) qh2vdb(2)
    data qh2vdb/4.874e-3_fp_kind, 1.371_fp_kind/
    ! parameter(nqh2i87 = 8)
    ! dimension qh2i87(nqh2i87), dqh2i87(nqh2i87-1), d2qh2i87(nqh2i87-2)
    ! data qh2i87/
    ! 11.69179D+00, -1.7227D+00, 7.98033D-01, -1.57089D-01, -5.35313D-01,
    ! 21.75818D+00, -2.63895D+00, 1.35708D+00/
    ! parameter(nqh2plusi89 = 10)
    ! dimension qh2plusi89(nqh2plusi89), dqh2plusi89(nqh2plusi89-1),&
    !   d2qh2plusi89(nqh2plusi89-2)
    ! data qh2plusi89/&
    !   2.564408220D+00, -2.15158152D+00, 4.046884D-01, 1.882055D+00,&
    !   -2.439623D+00, -3.47104D+00, 9.45846D+00, 6.971D-02,&
    !   -1.34950D+01, 8.59387D+00/
    integer nqh2i90
    parameter(nqh2i90 = 9)
    real(fp_kind) qh2i90(nqh2i90), dqh2i90(nqh2i90-1), d2qh2i90(nqh2i90-2)
    data qh2i90/&
         1.6918292822_fp_kind,-1.72246845_fp_kind,7.9631758e-01_fp_kind,&
         -1.706903e-01_fp_kind,-4.574240e-01_fp_kind,1.825633_fp_kind,&
         -3.52468_fp_kind,2.91102_fp_kind,-8.40929e-01_fp_kind/
    integer nqh2plusi90
    parameter(nqh2plusi90 = 11)
    real(fp_kind) qh2plusi90(nqh2plusi90), dqh2plusi90(nqh2plusi90-1),&
         d2qh2plusi90(nqh2plusi90-2)
    data qh2plusi90/&
         2.5644276281_fp_kind,-2.152026894_fp_kind,4.0264504e-01_fp_kind,&
         1.9297089_fp_kind,-2.432017_fp_kind,-4.784831_fp_kind,&
         1.208765e+01_fp_kind,8.31227_fp_kind,-4.94136e+01_fp_kind,5.46883e+01_fp_kind,-2.03073e+01_fp_kind/
    integer nqh2i90_hi
    parameter(nqh2i90_hi = 9)
    real(fp_kind) qh2i90_hi(nqh2i90_hi), dqh2i90_hi(nqh2i90_hi-1),&
         d2qh2i90_hi(nqh2i90_hi-2)
    data qh2i90_hi/&
         1.6936496106_fp_kind,-1.785296204_fp_kind,-3.3412934e-01_fp_kind,&
         -7.18046192_fp_kind,-2.13716663e+01_fp_kind,-2.73213417e+01_fp_kind,&
         -1.86763490e+01_fp_kind,-6.74105728_fp_kind,-1.01516624_fp_kind/
    integer nqh2plusi90_hi
    parameter(nqh2plusi90_hi = 8)
    real(fp_kind) qh2plusi90_hi(nqh2plusi90_hi),&
         dqh2plusi90_hi(nqh2plusi90_hi-1),&
         d2qh2plusi90_hi(nqh2plusi90_hi-2)
    data qh2plusi90_hi/&
         2.56934142_fp_kind,-2.0517792_fp_kind,1.239529_fp_kind,&
         5.653001_fp_kind,6.798597_fp_kind,4.355492_fp_kind,&
         1.508621_fp_kind,2.226748e-01_fp_kind/
    integer nqh2st
    parameter(nqh2st = 4)
    real(fp_kind) qh2st(nqh2st), dqh2st(nqh2st-1), d2qh2st(nqh2st-2)
    data qh2st/&
         1.6498_fp_kind, -1.6265_fp_kind, 0.7472_fp_kind, -0.2751_fp_kind/
    integer nqh2plusst
    parameter(nqh2plusst = 5)
    real(fp_kind) qh2plusst(nqh2plusst), dqh2plusst(nqh2plusst-1),&
         d2qh2plusst(nqh2plusst-2)
    data qh2plusst/&
         2.5410_fp_kind, -2.4336_fp_kind, 1.4979_fp_kind, 0.0192_fp_kind, -0.7483_fp_kind/
    integer nqh2pteh
    parameter (nqh2pteh = 5)
    real(fp_kind) qh2pteh(nqh2pteh)
    ! from PTEH paper with DH2 adjusted from 4.48 to 4.477 ev (Pols, private
    ! communication, 1996).
    data qh2pteh /&
         6608.8_fp_kind, 0.448_fp_kind, 0.1562_fp_kind, 0.0851_fp_kind, 4.477_fp_kind/
    integer iffirst, ifprint
    ! flag for first entry into routine
    data iffirst/1/
    ! flag for one printout for h2 or h2+ extrapolation beyond one dex
    ! above range.
    data ifprint/1/
    integer iorder, ifextrapolate, ifextrapolate_dex
    real(fp_kind) t, dtl, ti, thetal, eh2,&
         argzeta, zeta, dzeta, dzeta2, tmin, tmax
    ! Drop tmin slightly more than delta ln T = 0.1 below actual
    ! limit to allow test routines to evaluate derivatives numerically
    ! with fairly coarse step sizes right at T = 1.d3 without the
    ! extrapolated warning message that is otherwise emitted below.
    parameter (tmin=0.90e3_fp_kind, tmax=1.000000000000001e5_fp_kind)

    ! Most/all Fortran compilers specify the save attribute for all variables
    ! intialized by data statements.  But just in case...
    save qh2vdb, qh2i90, qh2plusi90, qh2i90_hi, qh2plusi90_hi, qh2st,  qh2plusst, qh2pteh, iffirst, ifprint
    ! Additional variables that are initialized below for iffirst.eq.1
    save dqh2i90, d2qh2i90, dqh2plusi90, d2qh2plusi90,  dqh2i90_hi, d2qh2i90_hi, dqh2plusi90_hi,&
         d2qh2plusi90_hi,  dqh2st, d2qh2st, dqh2plusst, d2qh2plusst

    if(iffirst.eq.1) then
       iffirst = 0
       do iorder = 1,nqh2i90-1
          dqh2i90(iorder) = real((iorder),fp_kind)*qh2i90(iorder+1)
       enddo
       do iorder = 1,nqh2i90-2
          d2qh2i90(iorder) = real((iorder),fp_kind)*dqh2i90(iorder+1)
       enddo
       do iorder = 1,nqh2plusi90-1
          dqh2plusi90(iorder) = real((iorder),fp_kind)*qh2plusi90(iorder+1)
       enddo
       do iorder = 1,nqh2plusi90-2
          d2qh2plusi90(iorder) = real((iorder),fp_kind)*dqh2plusi90(iorder+1)
       enddo
       do iorder = 1,nqh2i90_hi-1
          dqh2i90_hi(iorder) = real((iorder),fp_kind)*qh2i90_hi(iorder+1)
       enddo
       do iorder = 1,nqh2i90_hi-2
          d2qh2i90_hi(iorder) = real((iorder),fp_kind)*dqh2i90_hi(iorder+1)
       enddo
       do iorder = 1,nqh2plusi90_hi-1
          dqh2plusi90_hi(iorder) = real((iorder),fp_kind)*qh2plusi90_hi(iorder+1)
       enddo
       do iorder = 1,nqh2plusi90_hi-2
          d2qh2plusi90_hi(iorder) =&
               real((iorder),fp_kind)*dqh2plusi90_hi(iorder+1)
       enddo
       do iorder = 1,nqh2st-1
          dqh2st(iorder) = real((iorder),fp_kind)*qh2st(iorder+1)
       enddo
       do iorder = 1,nqh2st-2
          d2qh2st(iorder) = real((iorder),fp_kind)*dqh2st(iorder+1)
       enddo
       do iorder = 1,nqh2plusst-1
          dqh2plusst(iorder) = real((iorder),fp_kind)*qh2plusst(iorder+1)
       enddo
       do iorder = 1,nqh2plusst-2
          d2qh2plusst(iorder) = real((iorder),fp_kind)*dqh2plusst(iorder+1)
       enddo
       qh2pteh(1) = log(qh2pteh(1))
       do iorder = 2, nqh2pteh
          qh2pteh(iorder) = (ergsperev/boltzmann)*qh2pteh(iorder)
       enddo
    endif

    ! Sanity checks
    if(ifh2.lt.0.or.ifh2.gt.4) error stop 'molecular hydrogen: ifh2 must be in range from 0 to 4'
    if(.not.(ifh2plus.eq.ifh2-1.or.ifh2plus.eq.0.or.(ifh2plus.eq.2.and.ifh2.eq.4)))&
         error stop 'molecular_hydrogen: bad ifh2plus for input ifh2'

    ! n.b. this is a local t and should not affect anything outside this
    ! routine.
    t = exp(tl)

    ! warning: the ifh2 = 1 or 2 options are simply compatibility modes
    ! with old programmes.  Some effort has been made to keep these
    ! options, current, but they are no longer tested, and they might
    ! not work.

    ! The code inside the following nested if statements clearly
    ! initializes both ifextrapolate and ifextrapolate_dex for *all*
    ! cases, but gfortran apparently does not understand that logic so
    ! it emits spurious [-Wmaybe-uninitialized] ifextrapolate and
    ! ifextrapolate_dex warnings, and the following redundant
    ! initializations of ifextrapolate and ifextrapolate_dex are
    ! required to suppress those warnings.
    ifextrapolate = 0
    ifextrapolate_dex = 0

    if(ifh2.eq.1) then
       if(t.lt.tmin.or.t.gt.8.e3_fp_kind) then
          ifextrapolate = 1
          t=min(8.e3_fp_kind,max(tmin,t))
       else
          ifextrapolate = 0
       endif
       ifextrapolate_dex = ifextrapolate
    elseif(2.le.ifh2.and.ifh2.le.4) then
       if(ifh2.gt.2) then
          ! H2 and H2+ are negligible at 1.d6, but second-order
          ! extrapolation is ill-behaved beyond this temperature
          ! (underflows) so warn.
          ! n.b. get this test done *before* t is modified below.
          if(t.lt.tmin.or.t.gt.10._fp_kind*tmax) then
             ifextrapolate_dex = 1
          else
             ifextrapolate_dex = 0
          endif
       endif
       if(t.lt.tmin.or.t.gt.tmax) then
          ifextrapolate = 1
          t=min(tmax,max(tmin,t))
       else
          ifextrapolate = 0
       endif
       if(ifh2.eq.2) then
          ! this mode has tremendous errors at T = 1.d5 = tmax
          ! (original polynomial fit only to 9000 K) so warn.
          ! n.b. this test must be made *after* ifextrapolate is
          ! calculated.
          ifextrapolate_dex = ifextrapolate
       endif
    endif
    thetal = log10(5040._fp_kind/t)
    ti=11605.5_fp_kind/t
    if(ifh2.le.0.and.ifh2plus.gt.0) error stop 'molecular_hydrogen: invalid ifh2, ifh2plus switches'
    if(ifh2.eq.1) then
       ! vdb h2 partition function.
       eh2 = qh2vdb(2)/ti
       qh2 = log(qh2vdb(1)*t) + eh2
       qh2t = 1._fp_kind+eh2
       qh2tt = eh2
    elseif(ifh2.eq.2) then
       ! st h2 partition function.
       qh2 = poly_sum(thetal,qh2st)*ln10
       qh2t = -poly_sum(thetal,dqh2st)
       qh2tt = poly_sum(thetal,d2qh2st)/ln10
    elseif(ifh2.eq.3) then
       ! irwin (1990) h2 partition function (published in 1996)
       if(t.lt.9.e3_fp_kind) then
          qh2 = poly_sum(thetal,qh2i90)*ln10
          qh2t = -poly_sum(thetal,dqh2i90)
          qh2tt = poly_sum(thetal,d2qh2i90)/ln10
       else
          qh2 = poly_sum(thetal,qh2i90_hi)*ln10
          qh2t = -poly_sum(thetal,dqh2i90_hi)
          qh2tt = poly_sum(thetal,d2qh2i90_hi)/ln10
       endif
    elseif(ifh2.eq.4) then
       argzeta = qh2pteh(5)/t
       zeta = 1._fp_kind - exp(-argzeta)*(1._fp_kind + argzeta)
       qh2 = qh2pteh(1) + (qh2pteh(2) - (qh2pteh(3)*qh2pteh(3) -&
            qh2pteh(4)*qh2pteh(4)*qh2pteh(4)/t)/t)/t -&
            2.5_fp_kind*log(qh2pteh(5)/t) + log(zeta)
       dzeta = -exp(-argzeta)*argzeta*argzeta
       qh2t = -(qh2pteh(2) - (2._fp_kind*qh2pteh(3)*qh2pteh(3) -&
            3._fp_kind*qh2pteh(4)*qh2pteh(4)*qh2pteh(4)/t)/t)/t + 2.5_fp_kind +&
            dzeta/zeta
       dzeta2 = -exp(-argzeta)*argzeta*argzeta*(argzeta-2._fp_kind)
       qh2tt = (qh2pteh(2) - (4._fp_kind*qh2pteh(3)*qh2pteh(3) -&
            9._fp_kind*qh2pteh(4)*qh2pteh(4)*qh2pteh(4)/t)/t)/t +&
            (zeta*dzeta2 - dzeta*dzeta)/(zeta*zeta)
    endif
    if(ifh2.gt.0.and.ifextrapolate.eq.1) then
       if(ifh2.eq.3.or.ifh2.eq.4) then
          ! Taylor series approach for ln q
          ! n.b. t has been adjusted to maximum or minimum
          dtl = tl - log(t)
          qh2 = qh2 + dtl*qh2t + 0.5_fp_kind*dtl*dtl*qh2tt
          qh2t = qh2t + dtl*qh2tt
       else
          ! don't bother with old crummy partition functions
          qh2t = 0._fp_kind
          qh2tt = 0._fp_kind
       endif
    endif
    if(ifh2.gt.0.and.ifextrapolate_dex.eq.1) then
       if(ifprint.eq.1) then
          ifprint = 0
          if(verbosity.ge.2) then
             if(tl - log(t).gt.0._fp_kind) then
                write(stderr,'(a)')&
                     'WARNING: h2 p.f.extrapolated 1 dex above valid T range'
                if(ifh2plus.gt.0) write(stderr,'(a)')&
                     'WARNING: h2+ p.f. extrapolated 1 dex above valid T '//&
                     'range'
             else
                write(stderr,'(a)')&
                     'WARNING: h2 p.f.extrapolated below valid T range'
                if(ifh2plus.gt.0) write(stderr,'(a)')&
                     'WARNING: h2+ p.f. extrapolated below valid T range'
             endif
          endif
       endif
    endif
    if(ifh2plus.eq.1) then
       ! st h2plus partition function.
       qh2plus = poly_sum(thetal,qh2plusst)*ln10
       qh2plust = -poly_sum(thetal,dqh2plusst)
       qh2plustt = poly_sum(thetal,d2qh2plusst)/ln10
    elseif(ifh2plus.eq.2) then
       ! irwin (1990) h2plus partition function. (published in 1996)
       if(t.lt.9.e3_fp_kind) then
          qh2plus = poly_sum(thetal,qh2plusi90)*ln10
          qh2plust = -poly_sum(thetal,dqh2plusi90)
          qh2plustt = poly_sum(thetal,d2qh2plusi90)/ln10
       else
          qh2plus = poly_sum(thetal,qh2plusi90_hi)*ln10
          qh2plust = -poly_sum(thetal,dqh2plusi90_hi)
          qh2plustt = poly_sum(thetal,d2qh2plusi90_hi)/ln10
       endif
    endif
    if(ifh2plus.gt.0.and.ifextrapolate.eq.1) then
       if(ifh2plus.eq.2) then
          ! Taylor series approach for ln q
          ! n.b. t has been adjusted to maximum or minimum
          dtl = tl - log(t)
          qh2plus = qh2plus + dtl*qh2plust +&
               0.5_fp_kind*dtl*dtl*qh2plustt
          qh2plust = qh2plust + dtl*qh2plustt
       else
          ! don't bother with old crummy partition functions
          qh2plust = 0._fp_kind
          qh2plustt = 0._fp_kind
       endif
    endif
  end subroutine molecular_hydrogen
end module mod_molecular_hydrogen
