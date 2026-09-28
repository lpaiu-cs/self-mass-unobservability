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
!> This module provides parameters and variables that are used to help communicate between
!> the module procedures defined by the mod_excitation module.
module mod_excitation_block
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public c2t, x, x_old_excitation,&
       qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,&
       qstar, qstart, qstarx, qstart2, qstartx, qstarx2,&
       qmhd_he1, qmhd_he1t, qmhd_he1x,&
       qmhd_he1t2, qmhd_he1tx, qmhd_he1x2,&
       psum, psumf, psumt, psum_dv,&
       ssum, ssumf, ssumt, usum,&
       free_sum, free_sumf, free_sum_dv,&
       tl_old_excitation, ifpi_fit_old_excitation,&
       ifh2_old_excitation, ifh2plus_old_excitation,&
       ifpl_logical_old_excitation,&
       ifmhd_logical_old_excitation,&
       ifapprox_old_excitation, ifdiff_x_excitation,&
       ifhe1_special, nx, max_nmin_max, nions_excitationp2, max_izhi,&
       ifnr03, ifnr13, excitation_sum_called

  logical, parameter :: ifhe1_special = .true.

  ! number of auxiliary variables that qstar depends on
  integer, parameter :: nx = 5
  ! maximum value of nmin_max
  integer, parameter :: max_nmin_max = 10
  ! maximum value of izhi corresponds to nickel (28) plus room for one
  ! extra corresponding to special H2+ calculation.
  integer, parameter :: max_izhi = 29
  ! must be same as nionsp2
  integer, parameter:: nions_excitationp2 = 318

  real(fp_kind) c2t, x(nx), x_old_excitation(nx),&
       qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,&
       qstar(max_nmin_max, max_izhi),&
       qstart(max_nmin_max, max_izhi),&
       qstarx(nx, max_nmin_max, max_izhi),&
       qstart2(max_nmin_max, max_izhi),&
       qstartx(nx, max_nmin_max, max_izhi),&
       qstarx2(nx, nx, max_nmin_max, max_izhi),&
       qmhd_he1, qmhd_he1t, qmhd_he1x(nx),&
       qmhd_he1t2, qmhd_he1tx(nx), qmhd_he1x2(nx,nx),&
       psum, psumf, psumt, psum_dv(nions_excitationp2),&
       ssum, ssumf, ssumt, usum,&
       free_sum, free_sumf, free_sum_dv(nions_excitationp2)
  ! Variables for keeping track of conditions for h2/h2+ partition
  ! function calculations and qstar_calc calculations for
  ! excitation_pi and excitation_sum.
  real(fp_kind) tl_old_excitation
  integer ifpi_fit_old_excitation,&
       ifh2_old_excitation, ifh2plus_old_excitation
  logical ifpl_logical_old_excitation,&
       ifmhd_logical_old_excitation,&
       ifapprox_old_excitation, ifdiff_x_excitation,&
       ifnr03, ifnr13, excitation_sum_called
  data excitation_sum_called/.false./

  save c2t, x, x_old_excitation,&
       qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,&
       qstar, qstart, qstarx, qstart2, qstartx, qstarx2,&
       qmhd_he1, qmhd_he1t, qmhd_he1x,&
       qmhd_he1t2, qmhd_he1tx, qmhd_he1x2,&
       psum, psumf, psumt, psum_dv,&
       ssum, ssumf, ssumt, usum,&
       free_sum, free_sumf, free_sum_dv,&
       tl_old_excitation, ifpi_fit_old_excitation,&
       ifh2_old_excitation, ifh2plus_old_excitation,&
       ifpl_logical_old_excitation,&
       ifmhd_logical_old_excitation,&
       ifapprox_old_excitation, ifdiff_x_excitation,&
       ifnr03, ifnr13, excitation_sum_called
end module mod_excitation_block
