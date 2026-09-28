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

!> This module provides parameters describing constant ideal internal
!> partition functions (i.e., average statistical weights) of the
!> lower (i.e., non-hydrogenic) energy levels of monatomic neutral and
!> ion species.  These data are corrected elsewhere for hydrogenic
!> levels and non-ideal effects.
!>
!> FIXME(2021) The constant ideal partition functions given
!> here for the lower energy levels should be replaced by the correct
!> temperature-dependent ideal partition functions of those levels.

module mod_statistical_weight_data
  implicit none
  private
  public nelements_stat, iqneutral, nions_stat, iqion

  !> Number of elements
  integer, parameter :: nelements_stat = 24
  !> Constant ideal partition functions of lower electronic species for neutral species in element order.
  integer, parameter :: iqneutral(nelements_stat) = [&
       2, &
       1, &
       9, &
       4, &
       9, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       1, &
       21, &
       7, &
       6, &
       25, &
       21, &
       2, &
       1, &
       6, &
       6 ]

  !> Number of ionized species for the 20 elements
  integer, parameter :: nions_stat = 316
  !> Constant ideal partition functions of lower electronic species
  !> for ionized species in ion species order.
  integer, parameter :: iqion(nions_stat) = [&
       1, &
       2, &
       1, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       9, &
       2, &
       1, &
       2, &
       1, &
       2, &
       1, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       28, &
       21, &
       10, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       6, &
       25, &
       28, &
       21, &
       10, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       7, &
       6, &
       25, &
       28, &
       21, &
       10, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       30, &
       25, &
       6, &
       25, &
       28, &
       21, &
       10, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       10, &
       15, &
       30, &
       25, &
       6, &
       25, &
       28, &
       21, &
       10, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       6, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1, &
       1, &
       2, &
       1, &
       2, &
       1, &
       2, &
       1, &
       1, &
       2, &
       1, &
       2, &
       1, &
       9, &
       4, &
       9, &
       6, &
       1, &
       2, &
       1, &
       2, &
       1 ]
end module mod_statistical_weight_data
