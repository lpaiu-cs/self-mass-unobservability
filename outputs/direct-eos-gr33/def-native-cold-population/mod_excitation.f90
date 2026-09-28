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
end module mod_excitation
