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

!> This module provides number fraction data for all atomic and
!> molecular species used in the free-energy model for a particular
!> FreeEOS calculation.  These arrays are used to help communicate
!> these essential data between mod_eos_calc and mod_excitation module
!> routines and are also available for external use.

module mod_nuvar
  use mod_free_eos_types, only: fp_kind
  implicit none
  public audit_logweight, audit_zero, audit_reused, audit_ifnr, audit_calls
  real(fp_kind), save :: audit_logweight(29,24)=0._fp_kind
  integer, save :: audit_zero(29,24)=-1, audit_reused=0, audit_ifnr=-1, audit_calls=0
  private
  public nuvar, nuvarf, nuvart, nuvar_dv, nuvar_index_element, nuvar_atomic_number, nuvar_nelements, maxionstage_nuvar

  !     mod_nuvar is a module containing saved nuvar = n/(rho*navogadro) and related
  !       module variable information.
  !       nuvar is a function of dv, tl, but dv is a function of fl, and tl.
  !       for ifnr = 0:
  !       nuvarf is the fl derivative of nuvar(fl, tl)
  !       nuvart is the tl derivative of nuvar(fl, tl).
  !       for ifnr = 3:
  !       nuvarf is the fl derivative of nuvar(dv, tl) = 0
  !       nuvart is the tl derivative of nuvar(dv, tl).
  !       for ifnr = 1 or ifnr = 3:
  !       nuvar_dv is the dv derivative of nuvar(dv, tl).

  integer maxionstage_nuvar, nelements_nuvar
  !       must increase this if add element to mix with atomic number greater
  !       than 28
  parameter(maxionstage_nuvar = 28)
  !       must be the same as nelements in awieos_detailed.f
  parameter(nelements_nuvar = 24)
  ! nuvar_nelements is identical to npartial_elements (assigned
  !  to that value in ionize.)  Like npartial_elements, it is
  !  the number of different elements saved with the compact
  !  elements index.

  ! nuvar_index_element is identical to partial_elements (assigned
  !  to the same values as that array in ionize.)  Like that array
  !  it transforms from the compact element index to the element
  !  index.

  ! nuvar_atomic_number is the atomic number (or the number of neutral
  !  and ionized species up to next to bare ion).  It is assigned in
  !  ionize to be the same values as iatomic_number although it
  !  uses the compact element index rather than the iatomic_number
  !  ordinary element index.

  ! normally just use npartial_elements, partial_elements, and
  !  iatomic_number, but nuvar_nelements, nuvar_index_element, and
  !  nuvar_atomic_number assigned in ionize and saved here for
  !  the case of external programmes which will not have
  !  npartial_elements, partial_elements, and
  !  iatomic_number available.

  ! N.B. nelements_nuvar index is compact indexed to relevant elements.

  integer nuvar_index_element(nelements_nuvar),&
       nuvar_atomic_number(nelements_nuvar), nuvar_nelements
  ! N.B. first index refers to neutral, first ion, etc., up to next to
  !   bare ion (usually).  However, nuvar is also used in
  !   eos_free_calc where the bare ion is required (but no
  !   derivatives) so leave space for the bare ion value just in the
  !   nuvar case.
  real(fp_kind) nuvar(maxionstage_nuvar+1,nelements_nuvar),&
       nuvarf(maxionstage_nuvar,nelements_nuvar),&
       nuvart(maxionstage_nuvar,nelements_nuvar),&
       nuvar_dv(maxionstage_nuvar,maxionstage_nuvar,nelements_nuvar)

  save nuvar, nuvarf, nuvart, nuvar_dv, nuvar_index_element, nuvar_atomic_number, nuvar_nelements
end module mod_nuvar
