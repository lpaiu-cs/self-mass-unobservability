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
!
!*******************************************************************************

! calculate effective radii to be used in an MDH-style pressure ionization
! ifpi_fit = 2, use best fit to Saumon table
! ifpi_fit = 1, use best fit to opal table + extensions
! ifpi_fit = 0, use best fit to original MDH table.
! bi(nions+2) are the ionization potentials in cm^-1 (last two H2 and H2+).
! nion(nions+2) are the ionic charge of the parent ion (last two H2 and H2+).
! r_neutral(nelements+2) are the radii (cm) of the ground state of
!   the "neutral" species.  (These are actually species with
!   non-zero hard-sphere radii which happen to all be neutral in the
!   MHD model except for H2+.)  Aside from a proportionality fitting
!   factor these are the MDH values, i.e., aside from some
!   adopted values for H2, H2+, H, and He, the remaining values are
!   calculated from r = pi_fit*3 n^2 a0/2 Z, (hydrogenic values), where
!   n is the principal quantum number calculated from chi = Z^2 R/n^2.
!   The values are in order of the elements encountered for
!   bi, with the H2 and H2+ ground state radii appended as the
!   nelements+1 and nelements+2 values.
! r_ion3(nions+2) are the *cube of the*
!   effective radius (cm) of the ion-ion products.  aside from
!   a proportionality fitting factor this is the MDH value for the
!   ground state assuming K_n = unity.  The values are ordered
!   the same as bi, but the reference is to the lower ionization
!   state, with the bare nucleii skipped, e.g, h, he, he+, c-c+++++, etc.
!   the nions+1 and nions+2 values are for H2 and H2+

!> This effective_radius subroutine calculates effective radii
!> characterizing the MDH form of pressure ionization.
!>
!> \param[in] ifpi_fit PARAMETERS NEED DOCUMENTATION
!>
subroutine effective_radius(ifpi_fit, bi, nion, r_ion3, r_neutral)

  use mod_free_eos_constants, only: bohr, echarge, ergspercmm1, rydberg

  integer, intent(in) :: ifpi_fit, nion(:)
  real(fp_kind), intent(in) :: bi(:)
  real(fp_kind), intent(out) :: r_ion3(:), r_neutral(:)

  ! Local variables
  real(fp_kind), parameter :: rconst = 16._fp_kind**(1._fp_kind/3._fp_kind)*echarge*echarge/ergspercmm1

  integer ion, ielement, nionsp2, nions, nelements

  real(fp_kind) principal2

  nionsp2 = size(nion)
  nions = nionsp2 - 2
  nelements = size(r_neutral) - 2

  ! sanity checks
  if(nelements+2.ne.nelements_pi_fitp2) error stop 'effective_radius: bad size of r_neutral'
  if(nionsp2.ne.nions_pi_fitp2) error stop 'effective_radius: bad size of nion'
  if(nionsp2.ne.size(bi).or.nionsp2.ne.size(r_ion3)) error stop 'effective_radius: inconsistent nions sizes'

  ! calculate hydrogenic radii for all neutrals
  ielement = 0
  do ion = 1, nions
     if(nion(ion).eq.1) then
        ! only if first ion
        ielement = ielement + 1
        if(ielement.gt.nelements) error stop 'effective_radius: inconsistent input'
        principal2 = rydberg/bi(ion)
        r_neutral(ielement) = principal2*bohr*(1._fp_kind + 0.5_fp_kind/sqrt(principal2))
     endif
  enddo

  ! One more sanity check.
  if(ielement.ne.nelements) error stop 'effective_radius: inconsistent input'

  ! adopt MHD values (multiplied by fitting factor later) for
  ! H, He, H2, and H2+ radii
  ! r_neutral(1) = bohr (see Table 3 of MHD II paper for hydrogen
  ! species radii and text of same paper Section b) other elements
  ! for adopted helium ground-state radius).
  r_neutral(1) = 0.529e-8_fp_kind
  r_neutral(2) = 0.5e-8_fp_kind
  r_neutral(nelements+1) = 1.45e-8_fp_kind
  r_neutral(nelements+2) = 1.56e-8_fp_kind
  if(ifpi_fit.eq.0) then
     r_neutral = pi_fit_neutral_original*r_neutral
  elseif(ifpi_fit.eq.1) then
     r_neutral = pi_fit_neutral*r_neutral
  elseif(ifpi_fit.eq.2) then
     r_neutral = pi_fit_neutral_saumon*r_neutral
  endif

  ! calculate ion perturbation radii for all monatomic species
  ! (excluding bare nucleii) and for H2 and H2+.  Assume K_n = 1
  if(ifpi_fit.eq.0) then
     r_ion3 = (pi_fit_ion_original*sqrt(real(nion,fp_kind))*rconst/bi)**3
  elseif(ifpi_fit.eq.1) then
     r_ion3 = (pi_fit_ion*sqrt(real(nion,fp_kind))*rconst/bi)**3
  elseif(ifpi_fit.eq.2) then
     r_ion3 = (pi_fit_ion_saumon*sqrt(real(nion,fp_kind))*rconst/bi)**3
  else
     error stop 'effective_radius: bad ifpi_fit value input'
  endif

end subroutine effective_radius
