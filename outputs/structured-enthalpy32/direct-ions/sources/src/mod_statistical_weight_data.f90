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
       ! H
       2,&
       ! He
       1,&
       ! C
       9,&
       ! N
       4,&
       ! O
       9,&
       ! Ne
       1,&
       ! Na
       2,&
       ! Mg
       1,&
       ! Al
       6,&
       ! Si
       9,&
       ! P
       4,&
       ! S
       9,&
       ! Cl
       6,&
       ! A
       1,&
       ! Ca
       1,&
       ! Ti
       21,&
       ! Cr
       7,&
       ! Mn
       6,&
       ! Fe
       25,&
       ! Ni
       21,&
       2,1,6,6]

  !> Number of ionized species for the 20 elements
  integer, parameter :: nions_stat = 316
  !> Constant ideal partition functions of lower electronic species
  !> for ionized species in ion species order.
  integer, parameter :: iqion(nions_stat) = [&
       ! H
       1,&
       ! He
       2,1,&
       ! C
       6,1,2,1,2,1,&
       ! N
       9,6,1,2,1,2,1,&
       ! O
       4,9,6,1,2,1,2,1,&
       ! Ne
       6,9,4,9,6,1,2,1,2,1,&
       ! Na
       1,6,9,4,9,6,1,2,1,2,1,&
       ! Mg
       2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Al
       1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Si
       6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! P
       9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! S
       4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Cl
       9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! A
       6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Ca
       2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Ti
       28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Cr
       6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Mn
       7,6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Fe
       30,25,6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Ni
       10,15,30,25,6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,&
       ! Li
       1,2,1,&
       ! Be
       2,1,2,1,&
       ! B
       1,2,1,2,1,&
       ! F
       9,4,9,6,1,2,1,2,1]
end module mod_statistical_weight_data
