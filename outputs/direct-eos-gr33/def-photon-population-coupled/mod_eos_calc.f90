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
!> This module provides the public module procedure eos_calc.  That
!> procedure calculates ionization and dissociation fractions and a
!> new set of auxiliary variables as a function of non-ideal
!> equilibrium constants that are calculated from the old set of
!> auxiliary variables.
module mod_eos_calc
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public eos_calc
  real(fp_kind), save, public :: population_field(318)=0._fp_kind
contains
#include "eos_calc.f90"
#include "ionize.f90"
end module mod_eos_calc
