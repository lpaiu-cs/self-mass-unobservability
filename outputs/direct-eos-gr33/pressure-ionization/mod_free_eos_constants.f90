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
!> This module provides physical (and mathematical) constants required by FreeEOS.
!>
!> FIXME(2021).  These physical constants and certain parts of the
!> FreeEOS code that depend on them should be adjusted to be
!> consistent with the 2019 redefinition of the mole following the
!> discussion in <https://en.wikipedia.org/wiki/Molar_mass_constant>.
!> This redefinition should make negligible (O(10{-10})) relative
!> changes in FreeEOS results, but neverthless it would be good to get
!> this right.
module mod_free_eos_constants
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public &
       cr, compton, avogadro, pi, ln10, c_e, cd, clight,&
       electron_mass, h_mass, ct, cpe, c2, alpha20, alpha2, logalpha,&
       boltzmann, echarge, ergsperev, stefan_boltzmann, prad_const,&
       planck, rydberg, bohr, ergspercmm1, proton_mass
  ! math constant: 3.14159265358979323846264338327950288419... from
  ! <https://oeis.org/A000796/constant>, but 20 figures should be more than enough
  ! since even 80-bit (Intel) floating point has only ~18 decimal digits of precision.
  real(fp_kind), parameter :: pi = 3.1415926535897932385_fp_kind
  real(fp_kind), parameter :: ln10 = log(10._fp_kind)
  ! primary physical constants from 2018 CODATA
  ! <https://physics.nist.gov/cuu/pdf/wall_2018.pdf> which have been converted
  ! to CGS electrostatic units since FreeEOS has always used those units.

  ! exact values (so these should never change from now on).
  real(fp_kind), parameter :: clight = 2.99792458e10_fp_kind    !c (exact)
  real(fp_kind), parameter :: planck = 6.62607015e-27_fp_kind   !h (exact
  real(fp_kind), parameter :: avogadro = 6.02214076e23_fp_kind  !Avogadro number N_A (exact)
  real(fp_kind), parameter :: boltzmann = 1.380649e-16_fp_kind  !k (exact)
  ! real(fp_kind), parameter :: echarge = 4.8032068d-10    !Fritz's value in esu's
  ! 2018 CODATA (exact_ value transformed to electrostatic units.
  real(fp_kind), parameter :: echarge = 1.602176634e-19_fp_kind*1.e-1_fp_kind*clight ! e (exact)

  ! MAINTENANCE 2020-03.
  ! inexact (i.e., measured) values so these should be the only two values
  ! to change from now on as CODATA values evolve.
  real(fp_kind), parameter :: rydberg = 109737.31568160_fp_kind  !R(infinity)
  real(fp_kind), parameter :: proton_mass = 1.007276466621_fp_kind ! units are AMU

  ! secondary constants (all calculated from primary constants, although
  ! I have often specified the rounded (and therefore possibly
  ! thermodynamically inconsistent) CODATA value as a commented out
  ! value.
  ! real(fp_kind), parameter :: cr = 8.314462618d7   !R (gas constant) = k*N_A
  real(fp_kind), parameter :: cr = boltzmann*avogadro   !R (gas constant) = k*N_A
  ! a Joule is a Coulomb-volt, thus an electron volt is the number of
  ! Coulombs in a fundamental charge Joules, or 10^7 times that
  ! quantity ergs.
  ! real(fp_kind), parameter :: ergsperev = 1.60217653d-19*1.d7
  real(fp_kind), parameter :: ergsperev = (echarge/(1.e-1_fp_kind*clight))*1.e7_fp_kind
  ! 2002 CODATA electron mass in amu.
  ! real(fp_kind), parameter :: electron_mass = 5.4857990945d-04
  ! real(fp_kind), parameter :: electron_mass = 9.1093826d-28*avogadro
  real(fp_kind), parameter :: electron_mass = avogadro*rydberg*&
       (planck/(echarge*echarge))*&
       (planck/echarge)*(planck/echarge)*clight/(2._fp_kind*pi*pi)
  ! real(fp_kind), parameter :: h_mass = 1.007825035d0    !amu Wapstra and Audi, 1985
  real(fp_kind), parameter :: h_mass = proton_mass + electron_mass
  ! Compton wavelength of the electron.
  ! real(fp_kind), parameter :: compton = 2.426310238d-10  !2002 CODATA.
  real(fp_kind), parameter :: compton = planck*avogadro/(electron_mass*clight)
  real(fp_kind), parameter :: c_e = 8._fp_kind*pi/(compton*compton*compton)  ! ==> n_e = c_e*rhostar
  real(fp_kind), parameter :: cd = c_e/avogadro  !8 pi/(avogadro*compton^3)
  real(fp_kind), parameter :: ct = cr/(electron_mass*clight*clight)  !k/(m c^2)
  real(fp_kind), parameter :: cpe = cr*cd/ct  !8 pi m c^2/compton^3 ==> p_e = cpe*pstar
  ! second radiation constant
  ! real(fp_kind), parameter :: c2 = 1.4387752d0  !2002 CODATA.
  real(fp_kind), parameter :: c2 = planck*clight/boltzmann
  real(fp_kind), parameter :: alpha20 = c2*c2*cr/(2._fp_kind*pi*clight)
  ! alpha2 = (2piHk/h^2)^-3/H^2.
  real(fp_kind), parameter :: alpha2 = alpha20*alpha20*alpha20*&
       (avogadro/(clight*clight))*(avogadro/clight)
  real(fp_kind), parameter :: logalpha = 0.5_fp_kind*log(alpha2)
  real(fp_kind), parameter :: stefan_boltzmann = 2._fp_kind*pi*pi*pi*pi*pi/&
       (15._fp_kind*c2*c2*c2)*clight*boltzmann
  real(fp_kind), parameter :: prad_const = 4._fp_kind*stefan_boltzmann/(3._fp_kind*clight)
  ! real(fp_kind), parameter :: bohr = 0.5291772108d-8  !Bohr radius 2002 CODATA
  real(fp_kind), parameter :: bohr = 0.5_fp_kind*echarge*echarge/(rydberg*clight*planck)
  ! number of ergs per cm^-1 = planck*clight = c2*boltzmann
  real(fp_kind), parameter :: ergspercmm1 = c2*boltzmann
end module mod_free_eos_constants
