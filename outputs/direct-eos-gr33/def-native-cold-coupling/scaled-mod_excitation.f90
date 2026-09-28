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
!> This module provides module procedures to help determine the
!> non-ideal excitation component of the free-energy and all relevant
!> derivatives of same.
module mod_excitation
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public&
       excitation_pi, excitation_pi_pressure_free, excitation_pi_end,&
       excitation_sum
  real(fp_kind), save :: qstar_logscale(10,29)=0._fp_kind, qryd_logscale=0._fp_kind
contains
#include "excitation_pi.f90"
#include "excitation_pi_end.f90"
#include "excitation_sum.f90"
#include "qmhd_calc.f90"
#include "qstar_calc.f90"
#include "qryd_approx.f90"
#include "qryd_calc.f90"
#include "plsum.f90"
#include "plsum_approx.f90"

  function partition_product(logfactor, partition, logscale) result(value)
    real(fp_kind), intent(in) :: logfactor, partition, logscale
    real(fp_kind) :: value
    if(partition.le.0._fp_kind) then
       value=0._fp_kind
    elseif(logscale.eq.0._fp_kind.and.abs(logfactor).lt.700._fp_kind) then
       value=exp(logfactor)*partition
    else
       value=exp(logfactor+logscale+log(partition))
    endif
  end function partition_product

  subroutine scaled_summand(rn2,arg,lnocc,value,first,second)
    real(fp_kind), intent(in) :: rn2,arg,lnocc
    real(fp_kind), intent(out) :: value,first,second
    real(fp_kind) :: factor,em
    if(arg.gt.40._fp_kind) then
       factor=exp(arg+lnocc-qryd_logscale);em=exp(-arg)
       value=2._fp_kind*rn2*factor*(1._fp_kind-(1._fp_kind+arg)*em)
       first=2._fp_kind*factor*(1._fp_kind-em)
       second=2._fp_kind*factor/rn2
    else
       call plsummand_normalized(rn2,arg,value,first,second)
       factor=exp(lnocc-qryd_logscale)
       value=value*factor;first=first*factor;second=second*factor
    endif
  end subroutine scaled_summand
end module mod_excitation
