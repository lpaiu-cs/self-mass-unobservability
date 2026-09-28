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

!> This module provides isotopic mass data from AME2020 in convenient form.
module mod_isotopic_mass_data
  use mod_free_eos_types, only: fp_kind

  implicit none
  private
  public nelements_iso, isotopic_mass

  integer, parameter :: nelements_iso = 24

  ! These AME2020 isotopic mass data *for every known element* are derived
  ! as follows:

  ! This raw isotopic data was converted to equivalent parameter
  ! statements using the Fortran 2008 utility, generate_mod_isotopic_mass_data.  That utility
  ! is built by building the target of the same name, and that utility is run from the build tree
  ! using
  !
  ! # From the top directory in the build tree:
  ! cd test_FreeEOS/utils/
  ! sed -e 's?[*] ?0.?g' -e 's?#?.?g' <../../../free_eos.git/src/mass_1.mas20.txt |./generate_mod_isotopic_mass_data
  !
  ! where the two sed scripts convert the raw AME2020 data to a form that
  ! can be read by that utility, and the entire command generates the two
  ! files isotopic_mass_2d.data and isotopic_mass.data.  Those two files
  ! have been appropriately inserted into src/mod_isotopic_mass_data.f90
  ! to provide parameter arrays containing all AME2020 isotopic mass data
  ! for all elements, and selected AME2020 isotopic mass data for just the
  ! most abundant isotope of the 20 most abundant elements.  Note that for
  ! now, FreeEOS only uses the latter data, but that should change in the
  ! near future.
  ! START OF AME2020 DATA (inserted from the file isotopic_mass_2d.data documented above)
  integer, parameter :: min_nmz =  -8
  integer, parameter :: max_nmz =  64
  integer, parameter :: min_z =   0
  integer, parameter :: max_z = 118
  ! Z =   0
  real(fp_kind), parameter :: isotopic_mass_000(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   1.008664915900_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   1
  real(fp_kind), parameter :: isotopic_mass_001(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   1.007825031898_fp_kind,&
         2.014101777844_fp_kind,   3.016049281320_fp_kind,   4.026431867000_fp_kind,   5.035311492000_fp_kind,&
         6.044955437000_fp_kind,   7.052749000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   2
  real(fp_kind), parameter :: isotopic_mass_002(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   3.016029321970_fp_kind,&
         4.002603254130_fp_kind,   5.012057224000_fp_kind,   6.018885889000_fp_kind,   7.027990652000_fp_kind,&
         8.033934388000_fp_kind,   9.043946414000_fp_kind,  10.052815306000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   3
  real(fp_kind), parameter :: isotopic_mass_003(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   3.030775000000_fp_kind,   4.027185561000_fp_kind,   5.012537800000_fp_kind,&
         6.015122887420_fp_kind,   7.016003434260_fp_kind,   8.022486244000_fp_kind,   9.026790191000_fp_kind,&
        10.035483453000_fp_kind,  11.043723581000_fp_kind,  12.052613942000_fp_kind,  13.061171503000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   4
  real(fp_kind), parameter :: isotopic_mass_004(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   5.039870000000_fp_kind,   6.019726409000_fp_kind,   7.016928714000_fp_kind,&
         8.005305102000_fp_kind,   9.012183062000_fp_kind,  10.013534692000_fp_kind,  11.021661080000_fp_kind,&
        12.026922082000_fp_kind,  13.036134506000_fp_kind,  14.042892920000_fp_kind,  15.053490215000_fp_kind,&
        16.061672036000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   5
  real(fp_kind), parameter :: isotopic_mass_005(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         6.050800000000_fp_kind,   7.029712000000_fp_kind,   8.024607315000_fp_kind,   9.013329645000_fp_kind,&
        10.012936862000_fp_kind,  11.009305166000_fp_kind,  12.014352638000_fp_kind,  13.017779981000_fp_kind,&
        14.025404010000_fp_kind,  15.031087023000_fp_kind,  16.039841045000_fp_kind,  17.046931399000_fp_kind,&
        18.055601683000_fp_kind,  19.064166000000_fp_kind,  20.074505644000_fp_kind,  21.084147485000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   6
  real(fp_kind), parameter :: isotopic_mass_006(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         8.037643039000_fp_kind,   9.031037202000_fp_kind,  10.016853217000_fp_kind,  11.011432597000_fp_kind,&
        12.000000000000_fp_kind,  13.003354835340_fp_kind,  14.003241988620_fp_kind,  15.010599256000_fp_kind,&
        16.014701255000_fp_kind,  17.022578650000_fp_kind,  18.026751930000_fp_kind,  19.034797594000_fp_kind,&
        20.040261732000_fp_kind,  21.049000000000_fp_kind,  22.057553990000_fp_kind,  23.068890000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   7
  real(fp_kind), parameter :: isotopic_mass_007(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
        10.041653540000_fp_kind,  11.026157593000_fp_kind,  12.018613180000_fp_kind,  13.005738609000_fp_kind,&
        14.003074004250_fp_kind,  15.000108898270_fp_kind,  16.006101925000_fp_kind,  17.008448876000_fp_kind,&
        18.014077563000_fp_kind,  19.017022389000_fp_kind,  20.023367295000_fp_kind,  21.027087573000_fp_kind,&
        22.034100918000_fp_kind,  23.039421000000_fp_kind,  24.050390000000_fp_kind,  25.060100000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   8
  real(fp_kind), parameter :: isotopic_mass_008(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  11.051249828000_fp_kind,&
        12.034367726000_fp_kind,  13.024815435000_fp_kind,  14.008596706000_fp_kind,  15.003065636000_fp_kind,&
        15.994914619260_fp_kind,  16.999131755950_fp_kind,  17.999159612140_fp_kind,  19.003577969000_fp_kind,&
        20.004075357000_fp_kind,  21.008654948000_fp_kind,  22.009965744000_fp_kind,  23.015696686000_fp_kind,&
        24.019861000000_fp_kind,  25.029338919000_fp_kind,  26.037210155000_fp_kind,  27.047955000000_fp_kind,&
        28.055910000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =   9
  real(fp_kind), parameter :: isotopic_mass_009(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  13.045121000000_fp_kind,&
        14.034315196000_fp_kind,  15.017785139000_fp_kind,  16.011460278000_fp_kind,  17.002095237000_fp_kind,&
        18.000937324000_fp_kind,  18.998403162070_fp_kind,  19.999981252000_fp_kind,  20.999948893000_fp_kind,&
        22.002998812000_fp_kind,  23.003526875000_fp_kind,  24.008099370000_fp_kind,  25.012167727000_fp_kind,&
        26.020048065000_fp_kind,  27.026981897000_fp_kind,  28.035860448000_fp_kind,  29.043103000000_fp_kind,&
        30.052561000000_fp_kind,  31.061023000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  10
  real(fp_kind), parameter :: isotopic_mass_010(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  15.043172977000_fp_kind,&
        16.025750860000_fp_kind,  17.017713962000_fp_kind,  18.005708696000_fp_kind,  19.001880906000_fp_kind,&
        19.992440175250_fp_kind,  20.993846685000_fp_kind,  21.991385113000_fp_kind,  22.994466905000_fp_kind,&
        23.993610649000_fp_kind,  24.997814797000_fp_kind,  26.000516496000_fp_kind,  27.007569462000_fp_kind,&
        28.012130767000_fp_kind,  29.019753000000_fp_kind,  30.024992235000_fp_kind,  31.033474816000_fp_kind,&
        32.039720000000_fp_kind,  33.049523000000_fp_kind,  34.056728000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  11
  real(fp_kind), parameter :: isotopic_mass_011(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  17.037273000000_fp_kind,&
        18.026879388000_fp_kind,  19.013880264000_fp_kind,  20.007354301000_fp_kind,  20.997654459000_fp_kind,&
        21.994437547000_fp_kind,  22.989769281950_fp_kind,  23.990963012000_fp_kind,  24.989953974000_fp_kind,&
        25.992634649000_fp_kind,  26.994076408000_fp_kind,  27.998939000000_fp_kind,  29.002877091000_fp_kind,&
        30.009097931000_fp_kind,  31.013146654000_fp_kind,  32.020011024000_fp_kind,  33.025529000000_fp_kind,&
        34.034010000000_fp_kind,  35.040614000000_fp_kind,  36.049279000000_fp_kind,  37.057042000000_fp_kind,&
        38.066458000000_fp_kind,  39.075123000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  12
  real(fp_kind), parameter :: isotopic_mass_012(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  19.034179920000_fp_kind,&
        20.018763075000_fp_kind,  21.011705764000_fp_kind,  21.999570597000_fp_kind,  22.994123768000_fp_kind,&
        23.985041689000_fp_kind,  24.985836966000_fp_kind,  25.982592972000_fp_kind,  26.984340647000_fp_kind,&
        27.983875426000_fp_kind,  28.988607163000_fp_kind,  29.990465454000_fp_kind,  30.996648232000_fp_kind,&
        31.999110138000_fp_kind,  33.005327862000_fp_kind,  34.008935455000_fp_kind,  35.016790000000_fp_kind,&
        36.021879000000_fp_kind,  37.030286265000_fp_kind,  38.036580000000_fp_kind,  39.045921000000_fp_kind,&
        40.053194000000_fp_kind,  41.062373000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  13
  real(fp_kind), parameter :: isotopic_mass_013(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  21.029082000000_fp_kind,&
        22.019540000000_fp_kind,  23.007244351000_fp_kind,  23.999947598000_fp_kind,  24.990428308000_fp_kind,&
        25.986891876000_fp_kind,  26.981538408000_fp_kind,  27.981910009000_fp_kind,  28.980453164000_fp_kind,&
        29.982969171000_fp_kind,  30.983949754000_fp_kind,  31.988084338000_fp_kind,  32.990877685000_fp_kind,&
        33.996781924000_fp_kind,  34.999759816000_fp_kind,  36.006388000000_fp_kind,  37.010531000000_fp_kind,&
        38.017681000000_fp_kind,  39.023070000000_fp_kind,  40.030940000000_fp_kind,  41.037134000000_fp_kind,&
        42.045078000000_fp_kind,  43.051820000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  14
  real(fp_kind), parameter :: isotopic_mass_014(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  22.036114000000_fp_kind,  23.025711000000_fp_kind,&
        24.011535430000_fp_kind,  25.004108798000_fp_kind,  25.992333818000_fp_kind,  26.986704687000_fp_kind,&
        27.976926534420_fp_kind,  28.976494664340_fp_kind,  29.973770137000_fp_kind,  30.975363196000_fp_kind,&
        31.974151538000_fp_kind,  32.977976964000_fp_kind,  33.978538045000_fp_kind,  34.984550111000_fp_kind,&
        35.986649271000_fp_kind,  36.992945191000_fp_kind,  37.995523000000_fp_kind,  39.002491000000_fp_kind,&
        40.006083641000_fp_kind,  41.014171000000_fp_kind,  42.018078000000_fp_kind,  43.026119000000_fp_kind,&
        44.031466000000_fp_kind,  45.039818000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  15
  real(fp_kind), parameter :: isotopic_mass_015(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  24.036522000000_fp_kind,  25.021675000000_fp_kind,&
        26.011780000000_fp_kind,  26.999292499000_fp_kind,  27.992326460000_fp_kind,  28.981800368000_fp_kind,&
        29.978313490000_fp_kind,  30.973761997680_fp_kind,  31.973907643000_fp_kind,  32.971725692000_fp_kind,&
        33.973645886000_fp_kind,  34.973314045000_fp_kind,  35.978259610000_fp_kind,  36.979606942000_fp_kind,&
        37.984303105000_fp_kind,  38.986285865000_fp_kind,  39.991262221000_fp_kind,  40.994654000000_fp_kind,&
        42.001172140000_fp_kind,  43.005411000000_fp_kind,  44.011927000000_fp_kind,  45.017134000000_fp_kind,&
        46.024520000000_fp_kind,  47.030929000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  16
  real(fp_kind), parameter :: isotopic_mass_016(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  26.029716000000_fp_kind,  27.018777000000_fp_kind,&
        28.004372762000_fp_kind,  28.996678000000_fp_kind,  29.984906770000_fp_kind,  30.979557002000_fp_kind,&
        31.972071173540_fp_kind,  32.971458908620_fp_kind,  33.967867011000_fp_kind,  34.969032321000_fp_kind,&
        35.967080692000_fp_kind,  36.971125500000_fp_kind,  37.971163300000_fp_kind,  38.975133850000_fp_kind,&
        39.975482561000_fp_kind,  40.979593451000_fp_kind,  41.981065100000_fp_kind,  42.986907635000_fp_kind,&
        43.990118846000_fp_kind,  44.996414000000_fp_kind,  46.000687000000_fp_kind,  47.007730000000_fp_kind,&
        48.013301000000_fp_kind,  49.021891000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  17
  real(fp_kind), parameter :: isotopic_mass_017(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  28.030349000000_fp_kind,  29.015053000000_fp_kind,&
        30.005018333000_fp_kind,  30.992448097000_fp_kind,  31.985684605000_fp_kind,  32.977451988000_fp_kind,&
        33.973762490000_fp_kind,  34.968852694000_fp_kind,  35.968306822000_fp_kind,  36.965902573000_fp_kind,&
        37.968010408000_fp_kind,  38.968008151000_fp_kind,  39.970415466000_fp_kind,  40.970684525000_fp_kind,&
        41.973342000000_fp_kind,  42.974063700000_fp_kind,  43.978014918000_fp_kind,  44.980394353000_fp_kind,&
        45.985254926000_fp_kind,  46.989715000000_fp_kind,  47.995405000000_fp_kind,  49.000794000000_fp_kind,&
        50.008266000000_fp_kind,  51.015341000000_fp_kind,  52.024004000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  18
  real(fp_kind), parameter :: isotopic_mass_018(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  29.040761000000_fp_kind,  30.023694000000_fp_kind,  31.012158000000_fp_kind,&
        31.997637824000_fp_kind,  32.989925545000_fp_kind,  33.980270092000_fp_kind,  34.975257719000_fp_kind,&
        35.967545106000_fp_kind,  36.966776301000_fp_kind,  37.962732102000_fp_kind,  38.964313037000_fp_kind,&
        39.962383122040_fp_kind,  40.964500570000_fp_kind,  41.963045737000_fp_kind,  42.965636056000_fp_kind,&
        43.964923814000_fp_kind,  44.968039731000_fp_kind,  45.968039244000_fp_kind,  46.972767112000_fp_kind,&
        47.976001000000_fp_kind,  48.981685000000_fp_kind,  49.985797000000_fp_kind,  50.993033000000_fp_kind,&
        51.998519000000_fp_kind,  53.007290000000_fp_kind,  54.013484000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  19
  real(fp_kind), parameter :: isotopic_mass_019(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  31.036780000000_fp_kind,  32.023607000000_fp_kind,  33.008095000000_fp_kind,&
        33.998690000000_fp_kind,  34.988005406000_fp_kind,  35.981301887000_fp_kind,  36.973375890000_fp_kind,&
        37.969081114000_fp_kind,  38.963706484820_fp_kind,  39.963998165000_fp_kind,  40.961825256110_fp_kind,&
        41.962402305000_fp_kind,  42.960734701000_fp_kind,  43.961586984000_fp_kind,  44.960691491000_fp_kind,&
        45.961981584000_fp_kind,  46.961661612000_fp_kind,  47.965341184000_fp_kind,  48.968210753000_fp_kind,&
        49.972380015000_fp_kind,  50.975828664000_fp_kind,  51.981602000000_fp_kind,  52.986800000000_fp_kind,&
        53.994471000000_fp_kind,  55.000505000000_fp_kind,  56.008567000000_fp_kind,  57.015169000000_fp_kind,&
        58.023543000000_fp_kind,  59.030864000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  20
  real(fp_kind), parameter :: isotopic_mass_020(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  33.033312000000_fp_kind,  34.015985000000_fp_kind,  35.005572000000_fp_kind,&
        35.993074388000_fp_kind,  36.985897849000_fp_kind,  37.976319223000_fp_kind,  38.970710811000_fp_kind,&
        39.962590850000_fp_kind,  40.962277905000_fp_kind,  41.958617780000_fp_kind,  42.958766381000_fp_kind,&
        43.955481489000_fp_kind,  44.956186270000_fp_kind,  45.953687726000_fp_kind,  46.954541134000_fp_kind,&
        47.952522654000_fp_kind,  48.955662625000_fp_kind,  49.957499215000_fp_kind,  50.960995663000_fp_kind,&
        51.963213646000_fp_kind,  52.968451000000_fp_kind,  53.972989000000_fp_kind,  54.979978000000_fp_kind,&
        55.985496000000_fp_kind,  56.992958000000_fp_kind,  57.998357000000_fp_kind,  59.006237000000_fp_kind,&
        60.011809000000_fp_kind,  61.020408000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  21
  real(fp_kind), parameter :: isotopic_mass_021(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  35.029093000000_fp_kind,  36.017338000000_fp_kind,  37.004058000000_fp_kind,&
        37.995438000000_fp_kind,  38.984784953000_fp_kind,  39.977967275000_fp_kind,  40.969251163000_fp_kind,&
        41.965516686000_fp_kind,  42.961150425000_fp_kind,  43.959402818000_fp_kind,  44.955907051000_fp_kind,&
        45.955167034000_fp_kind,  46.952402444000_fp_kind,  47.952222903000_fp_kind,  48.950013159000_fp_kind,&
        49.952187437000_fp_kind,  50.953568838000_fp_kind,  51.956496170000_fp_kind,  52.958379173000_fp_kind,&
        53.963029359000_fp_kind,  54.966889637000_fp_kind,  55.972607611000_fp_kind,  56.977048000000_fp_kind,&
        57.983382000000_fp_kind,  58.988374000000_fp_kind,  59.995115000000_fp_kind,  61.000537000000_fp_kind,&
        62.007848000000_fp_kind,  63.014031000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  22
  real(fp_kind), parameter :: isotopic_mass_022(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  37.027021000000_fp_kind,  38.012206000000_fp_kind,  39.002684000000_fp_kind,&
        39.990345146000_fp_kind,  40.983148000000_fp_kind,  41.973049369000_fp_kind,  42.968528420000_fp_kind,&
        43.959689936000_fp_kind,  44.958120758000_fp_kind,  45.952626356000_fp_kind,  46.951757491000_fp_kind,&
        47.947940677000_fp_kind,  48.947864391000_fp_kind,  49.944785622000_fp_kind,  50.946609468000_fp_kind,&
        51.946883509000_fp_kind,  52.949670714000_fp_kind,  53.950892000000_fp_kind,  54.955091000000_fp_kind,&
        55.957677675000_fp_kind,  56.963068098000_fp_kind,  57.966808519000_fp_kind,  58.972217000000_fp_kind,&
        59.976275000000_fp_kind,  60.982426000000_fp_kind,  61.986903000000_fp_kind,  62.993709000000_fp_kind,&
        63.998411000000_fp_kind,  65.005593000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  23
  real(fp_kind), parameter :: isotopic_mass_023(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  39.024230000000_fp_kind,  40.013387000000_fp_kind,  41.000333000000_fp_kind,&
        41.991820000000_fp_kind,  42.980766000000_fp_kind,  43.974440977000_fp_kind,  44.965768498000_fp_kind,&
        45.960197389000_fp_kind,  46.954903558000_fp_kind,  47.952250900000_fp_kind,  48.948510509000_fp_kind,&
        49.947156681000_fp_kind,  50.943957664000_fp_kind,  51.944773636000_fp_kind,  52.944334940000_fp_kind,&
        53.946432009000_fp_kind,  54.947262000000_fp_kind,  55.950420082000_fp_kind,  56.952297000000_fp_kind,&
        57.956595985000_fp_kind,  58.959623343000_fp_kind,  59.964479215000_fp_kind,  60.967603529000_fp_kind,&
        61.972932556000_fp_kind,  62.976661000000_fp_kind,  63.982480000000_fp_kind,  64.986999000000_fp_kind,&
        65.993237000000_fp_kind,  66.998128000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  24
  real(fp_kind), parameter :: isotopic_mass_024(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  41.021911000000_fp_kind,  42.007579000000_fp_kind,  42.997885000000_fp_kind,&
        43.985591000000_fp_kind,  44.979050000000_fp_kind,  45.968360969000_fp_kind,  46.962894995000_fp_kind,&
        47.954029431000_fp_kind,  48.951333720000_fp_kind,  49.946042209000_fp_kind,  50.944765388000_fp_kind,&
        51.940504714000_fp_kind,  52.940646304000_fp_kind,  53.938877359000_fp_kind,  54.940836637000_fp_kind,&
        55.940648977000_fp_kind,  56.943612112000_fp_kind,  57.944184501000_fp_kind,  58.948345426000_fp_kind,&
        59.949641656000_fp_kind,  60.954378130000_fp_kind,  61.956142920000_fp_kind,  62.961161000000_fp_kind,&
        63.963886000000_fp_kind,  64.969608000000_fp_kind,  65.973011000000_fp_kind,  66.979313000000_fp_kind,&
        67.983156000000_fp_kind,  68.989662000000_fp_kind,  69.993945000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  25
  real(fp_kind), parameter :: isotopic_mass_025(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  43.018647000000_fp_kind,  44.008009000000_fp_kind,  44.994654000000_fp_kind,&
        45.986669000000_fp_kind,  46.975774000000_fp_kind,  47.968548760000_fp_kind,  48.959613350000_fp_kind,&
        49.954238157000_fp_kind,  50.948208770000_fp_kind,  51.945559090000_fp_kind,  52.941287497000_fp_kind,&
        53.940355772000_fp_kind,  54.938043040000_fp_kind,  55.938902816000_fp_kind,  56.938285944000_fp_kind,&
        57.940066643000_fp_kind,  58.940391111000_fp_kind,  59.943136574000_fp_kind,  60.944452541000_fp_kind,&
        61.947907384000_fp_kind,  62.949664672000_fp_kind,  63.953849369000_fp_kind,  64.956019749000_fp_kind,&
        65.960546833000_fp_kind,  66.963950000000_fp_kind,  67.968953000000_fp_kind,  68.972775000000_fp_kind,&
        69.978046000000_fp_kind,  70.982158000000_fp_kind,  71.988009000000_fp_kind,  72.992807000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  26
  real(fp_kind), parameter :: isotopic_mass_026(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  45.015467000000_fp_kind,  46.001299000000_fp_kind,  46.992346000000_fp_kind,&
        47.980667000000_fp_kind,  48.973429000000_fp_kind,  49.962988000000_fp_kind,  50.956855137000_fp_kind,&
        51.948113364000_fp_kind,  52.945305629000_fp_kind,  53.939608189000_fp_kind,  54.938291158000_fp_kind,&
        55.934935537000_fp_kind,  56.935391950000_fp_kind,  57.933273575000_fp_kind,  58.934873492000_fp_kind,&
        59.934070249000_fp_kind,  60.936746241000_fp_kind,  61.936791809000_fp_kind,  62.940272698000_fp_kind,&
        63.940987761000_fp_kind,  64.945015323000_fp_kind,  65.946249958000_fp_kind,  66.950930000000_fp_kind,&
        67.952875000000_fp_kind,  68.957918000000_fp_kind,  69.960397000000_fp_kind,  70.965722000000_fp_kind,&
        71.968599000000_fp_kind,  72.974246000000_fp_kind,  73.977821000000_fp_kind,  74.984219000000_fp_kind,&
        75.988631000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  27
  real(fp_kind), parameter :: isotopic_mass_027(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,  47.011401000000_fp_kind,  48.001857000000_fp_kind,  48.989501000000_fp_kind,&
        49.981117000000_fp_kind,  50.970647000000_fp_kind,  51.963130224000_fp_kind,  52.954203278000_fp_kind,&
        53.948459075000_fp_kind,  54.941996416000_fp_kind,  55.939838032000_fp_kind,  56.936289819000_fp_kind,&
        57.935751292000_fp_kind,  58.933193524000_fp_kind,  59.933815536000_fp_kind,  60.932476031000_fp_kind,&
        61.934058198000_fp_kind,  62.933599630000_fp_kind,  63.935810176000_fp_kind,  64.936462071000_fp_kind,&
        65.939442943000_fp_kind,  66.940609625000_fp_kind,  67.944559401000_fp_kind,  68.945909000000_fp_kind,&
        69.950053400000_fp_kind,  70.952366923000_fp_kind,  71.956736000000_fp_kind,  72.959238000000_fp_kind,&
        73.963993000000_fp_kind,  74.967192000000_fp_kind,  75.972453000000_fp_kind,  76.976479000000_fp_kind,&
        77.983553000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  28
  real(fp_kind), parameter :: isotopic_mass_028(min_nmz:max_nmz) = [&
        48.019515000000_fp_kind,  49.009157000000_fp_kind,  49.996286000000_fp_kind,  50.987493000000_fp_kind,&
        51.975781000000_fp_kind,  52.968190000000_fp_kind,  53.957833000000_fp_kind,  54.951329846000_fp_kind,&
        55.942127761000_fp_kind,  56.939791394000_fp_kind,  57.935341650000_fp_kind,  58.934345442000_fp_kind,&
        59.930785129000_fp_kind,  60.931054819000_fp_kind,  61.928344753000_fp_kind,  62.929669021000_fp_kind,&
        63.927966228000_fp_kind,  64.930084585000_fp_kind,  65.929139333000_fp_kind,  66.931569413000_fp_kind,&
        67.931868787000_fp_kind,  68.935610267000_fp_kind,  69.936431300000_fp_kind,  70.940518962000_fp_kind,&
        71.941785924000_fp_kind,  72.946206681000_fp_kind,  73.947718000000_fp_kind,  74.952506000000_fp_kind,&
        75.954707000000_fp_kind,  76.959903000000_fp_kind,  77.962555000000_fp_kind,  78.969769000000_fp_kind,&
        79.975051000000_fp_kind,  80.982727000000_fp_kind,  81.988492000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  29
  real(fp_kind), parameter :: isotopic_mass_029(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  51.997982000000_fp_kind,  52.985894000000_fp_kind,&
        53.977198000000_fp_kind,  54.966038000000_fp_kind,  55.958529278000_fp_kind,  56.949211686000_fp_kind,&
        57.944532283000_fp_kind,  58.939496713000_fp_kind,  59.937363787000_fp_kind,  60.933457375000_fp_kind,&
        61.932594803000_fp_kind,  62.929597119000_fp_kind,  63.929764001000_fp_kind,  64.927789476000_fp_kind,&
        65.928868804000_fp_kind,  66.927729490000_fp_kind,  67.929610887000_fp_kind,  68.929429267000_fp_kind,&
        69.932392078000_fp_kind,  70.932676831000_fp_kind,  71.935820306000_fp_kind,  72.936674376000_fp_kind,&
        73.939874860000_fp_kind,  74.941523817000_fp_kind,  75.945268974000_fp_kind,  76.947543599000_fp_kind,&
        77.951916524000_fp_kind,  78.954473100000_fp_kind,  79.960623000000_fp_kind,  80.965743000000_fp_kind,&
        81.972378000000_fp_kind,  82.978110000000_fp_kind,  83.985271000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  30
  real(fp_kind), parameter :: isotopic_mass_030(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  53.993879000000_fp_kind,  54.984681000000_fp_kind,&
        55.972743000000_fp_kind,  56.965056000000_fp_kind,  57.954590296000_fp_kind,  58.949311886000_fp_kind,&
        59.941841317000_fp_kind,  60.939506964000_fp_kind,  61.934333359000_fp_kind,  62.933211140000_fp_kind,&
        63.929141776000_fp_kind,  64.929240534000_fp_kind,  65.926033639000_fp_kind,  66.927127422000_fp_kind,&
        67.924844232000_fp_kind,  68.926550360000_fp_kind,  69.925319175000_fp_kind,  70.927719578000_fp_kind,&
        71.926842806000_fp_kind,  72.929582580000_fp_kind,  73.929407260000_fp_kind,  74.932840244000_fp_kind,&
        75.933114956000_fp_kind,  76.936887197000_fp_kind,  77.938289204000_fp_kind,  78.942638067000_fp_kind,&
        79.944552929000_fp_kind,  80.950402617000_fp_kind,  81.954574097000_fp_kind,  82.961041000000_fp_kind,&
        83.965829000000_fp_kind,  84.973054000000_fp_kind,  85.978463000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  31
  real(fp_kind), parameter :: isotopic_mass_031(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  55.995878000000_fp_kind,  56.983457000000_fp_kind,&
        57.974729000000_fp_kind,  58.963757000000_fp_kind,  59.957498000000_fp_kind,  60.949398861000_fp_kind,&
        61.944189639000_fp_kind,  62.939294194000_fp_kind,  63.936840366000_fp_kind,  64.932734424000_fp_kind,&
        65.931589766000_fp_kind,  66.928202276000_fp_kind,  67.927980161000_fp_kind,  68.925573528000_fp_kind,&
        69.926021914000_fp_kind,  70.924702554000_fp_kind,  71.926367452000_fp_kind,  72.925174680000_fp_kind,&
        73.926945725000_fp_kind,  74.926504484000_fp_kind,  75.928827624000_fp_kind,  76.929154299000_fp_kind,&
        77.931610854000_fp_kind,  78.932851582000_fp_kind,  79.936420773000_fp_kind,  80.938133841000_fp_kind,&
        81.943176531000_fp_kind,  82.947120300000_fp_kind,  83.952663000000_fp_kind,  84.957333000000_fp_kind,&
        85.963757000000_fp_kind,  86.969007000000_fp_kind,  87.975963000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  32
  real(fp_kind), parameter :: isotopic_mass_032(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  57.991863000000_fp_kind,  58.982426000000_fp_kind,&
        59.970445000000_fp_kind,  60.963725000000_fp_kind,  61.954761000000_fp_kind,  62.949628000000_fp_kind,&
        63.941689912000_fp_kind,  64.939368136000_fp_kind,  65.933862124000_fp_kind,  66.932716999000_fp_kind,&
        67.928095305000_fp_kind,  68.927964467000_fp_kind,  69.924248542000_fp_kind,  70.924952120000_fp_kind,&
        71.922075824000_fp_kind,  72.923458954000_fp_kind,  73.921177760000_fp_kind,  74.922858370000_fp_kind,&
        75.921402725000_fp_kind,  76.923549843000_fp_kind,  77.922852911000_fp_kind,  78.925359506000_fp_kind,&
        79.925350773000_fp_kind,  80.928832941000_fp_kind,  81.929774031000_fp_kind,  82.934539100000_fp_kind,&
        83.937575090000_fp_kind,  84.942969658000_fp_kind,  85.946967000000_fp_kind,  86.953204000000_fp_kind,&
        87.957574000000_fp_kind,  88.964530000000_fp_kind,  89.969436000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  33
  real(fp_kind), parameter :: isotopic_mass_033(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  59.993945000000_fp_kind,  60.981535000000_fp_kind,&
        61.973784000000_fp_kind,  62.964036000000_fp_kind,  63.957560000000_fp_kind,  64.949611000000_fp_kind,&
        65.944148778000_fp_kind,  66.939251110000_fp_kind,  67.936774127000_fp_kind,  68.932246289000_fp_kind,&
        69.930934642000_fp_kind,  70.927113594000_fp_kind,  71.926752291000_fp_kind,  72.923829086000_fp_kind,&
        73.923928596000_fp_kind,  74.921594562000_fp_kind,  75.922392011000_fp_kind,  76.920647555000_fp_kind,&
        77.921827771000_fp_kind,  78.920948419000_fp_kind,  79.922474440000_fp_kind,  80.922132288000_fp_kind,&
        81.924738731000_fp_kind,  82.925206900000_fp_kind,  83.929303290000_fp_kind,  84.932163658000_fp_kind,&
        85.936701532000_fp_kind,  86.940291716000_fp_kind,  87.945840000000_fp_kind,  88.950048000000_fp_kind,&
        89.955995000000_fp_kind,  90.960816000000_fp_kind,  91.967386000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  34
  real(fp_kind), parameter :: isotopic_mass_034(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  62.981911000000_fp_kind,&
        63.971165000000_fp_kind,  64.964552000000_fp_kind,  65.955276000000_fp_kind,  66.949994000000_fp_kind,&
        67.941825236000_fp_kind,  68.939414845000_fp_kind,  69.933515521000_fp_kind,  70.932209431000_fp_kind,&
        71.927140506000_fp_kind,  72.926754881000_fp_kind,  73.922475933000_fp_kind,  74.922522870000_fp_kind,&
        75.919213702000_fp_kind,  76.919914150000_fp_kind,  77.917309244000_fp_kind,  78.918499252000_fp_kind,&
        79.916521761000_fp_kind,  80.917993019000_fp_kind,  81.916699531000_fp_kind,  82.919118604000_fp_kind,&
        83.918466761000_fp_kind,  84.922260758000_fp_kind,  85.924311732000_fp_kind,  86.928688616000_fp_kind,&
        87.931417490000_fp_kind,  88.936669058000_fp_kind,  89.940096000000_fp_kind,  90.945700000000_fp_kind,&
        91.949840000000_fp_kind,  92.956135000000_fp_kind,  93.960490000000_fp_kind,  94.967300000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  35
  real(fp_kind), parameter :: isotopic_mass_035(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  64.982297000000_fp_kind,&
        65.974697000000_fp_kind,  66.965078000000_fp_kind,  67.958356000000_fp_kind,  68.950338410000_fp_kind,&
        69.944792321000_fp_kind,  70.939342153000_fp_kind,  71.936594606000_fp_kind,  72.931673441000_fp_kind,&
        73.929910279000_fp_kind,  74.925810566000_fp_kind,  75.924541574000_fp_kind,  76.921379193000_fp_kind,&
        77.921145858000_fp_kind,  78.918337574000_fp_kind,  79.918529784000_fp_kind,  80.916288197000_fp_kind,&
        81.916801752000_fp_kind,  82.915175285000_fp_kind,  83.916496417000_fp_kind,  84.915645758000_fp_kind,&
        85.918805432000_fp_kind,  86.920674016000_fp_kind,  87.924083290000_fp_kind,  88.926704558000_fp_kind,&
        89.931292848000_fp_kind,  90.934398617000_fp_kind,  91.939631595000_fp_kind,  92.943220000000_fp_kind,&
        93.948846000000_fp_kind,  94.952925000000_fp_kind,  95.958980000000_fp_kind,  96.963499000000_fp_kind,&
        97.969887000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  36
  real(fp_kind), parameter :: isotopic_mass_036(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  66.983305000000_fp_kind,&
        67.972489000000_fp_kind,  68.965496000000_fp_kind,  69.955877000000_fp_kind,  70.950265695000_fp_kind,&
        71.942092406000_fp_kind,  72.939289193000_fp_kind,  73.933084016000_fp_kind,  74.930945744000_fp_kind,&
        75.925910743000_fp_kind,  76.924669999000_fp_kind,  77.920366341000_fp_kind,  78.920082919000_fp_kind,&
        79.916377940000_fp_kind,  80.916589703000_fp_kind,  81.913481153680_fp_kind,  82.914126516000_fp_kind,&
        83.911497727080_fp_kind,  84.912527260000_fp_kind,  85.910610624680_fp_kind,  86.913354759000_fp_kind,&
        87.914447879000_fp_kind,  88.917835449000_fp_kind,  89.919527929000_fp_kind,  90.923806309000_fp_kind,&
        91.926173092000_fp_kind,  92.931147172000_fp_kind,  93.934140452000_fp_kind,  94.939710922000_fp_kind,&
        95.943014473000_fp_kind,  96.949088782000_fp_kind,  97.952635000000_fp_kind,  98.958776000000_fp_kind,&
        99.962995000000_fp_kind, 100.969318000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  37
  real(fp_kind), parameter :: isotopic_mass_037(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  70.965335000000_fp_kind,  71.958851000000_fp_kind,  72.950604506000_fp_kind,&
        73.944265867000_fp_kind,  74.938573200000_fp_kind,  75.935073031000_fp_kind,  76.930401599000_fp_kind,&
        77.928141866000_fp_kind,  78.923990095000_fp_kind,  79.922516442000_fp_kind,  80.918993900000_fp_kind,&
        81.918209023000_fp_kind,  82.915114181000_fp_kind,  83.914375223000_fp_kind,  84.911789736040_fp_kind,&
        85.911167443000_fp_kind,  86.909180529000_fp_kind,  87.911315590000_fp_kind,  88.912278136000_fp_kind,&
        89.914797557000_fp_kind,  90.916537261000_fp_kind,  91.919728477000_fp_kind,  92.922039334000_fp_kind,&
        93.926394819000_fp_kind,  94.929263849000_fp_kind,  95.934133398000_fp_kind,  96.937177117000_fp_kind,&
        97.941632317000_fp_kind,  98.945119190000_fp_kind,  99.950331532000_fp_kind, 100.954302000000_fp_kind,&
       101.960008000000_fp_kind, 102.964401000000_fp_kind, 103.970531000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  38
  real(fp_kind), parameter :: isotopic_mass_038(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  72.965700000000_fp_kind,  73.956170000000_fp_kind,  74.949952767000_fp_kind,&
        75.941762760000_fp_kind,  76.937945454000_fp_kind,  77.932179979000_fp_kind,  78.929704692000_fp_kind,&
        79.924517538000_fp_kind,  80.923211393000_fp_kind,  81.918399845000_fp_kind,  82.917554372000_fp_kind,&
        83.913419118000_fp_kind,  84.912932041000_fp_kind,  85.909260724730_fp_kind,  86.908877494540_fp_kind,&
        87.905612253000_fp_kind,  88.907450808000_fp_kind,  89.907727870000_fp_kind,  90.910195942000_fp_kind,&
        91.911038222000_fp_kind,  92.914024314000_fp_kind,  93.915355641000_fp_kind,  94.919358282000_fp_kind,&
        95.921719045000_fp_kind,  96.926375621000_fp_kind,  97.928692636000_fp_kind,  98.932883604000_fp_kind,&
        99.935783270000_fp_kind, 100.940606264000_fp_kind, 101.944004679000_fp_kind, 102.949243000000_fp_kind,&
       103.953022000000_fp_kind, 104.959001000000_fp_kind, 105.963177000000_fp_kind, 106.969672000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  39
  real(fp_kind), parameter :: isotopic_mass_039(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  74.965840000000_fp_kind,  75.958937000000_fp_kind,  76.950146000000_fp_kind,&
        77.943990000000_fp_kind,  78.937946000000_fp_kind,  79.934354750000_fp_kind,  80.929454283000_fp_kind,&
        81.926930189000_fp_kind,  82.922484026000_fp_kind,  83.920671060000_fp_kind,  84.916433039000_fp_kind,&
        85.914886095000_fp_kind,  86.910876100000_fp_kind,  87.909501274000_fp_kind,  88.905838156000_fp_kind,&
        89.907141749000_fp_kind,  90.907298048000_fp_kind,  91.908945752000_fp_kind,  92.909578434000_fp_kind,&
        93.911592062000_fp_kind,  94.912819697000_fp_kind,  95.915909305000_fp_kind,  96.918286702000_fp_kind,&
        97.922394841000_fp_kind,  98.924160839000_fp_kind,  99.927727678000_fp_kind, 100.930160817000_fp_kind,&
       101.934328471000_fp_kind, 102.937243796000_fp_kind, 103.941943000000_fp_kind, 104.945711000000_fp_kind,&
       105.950842000000_fp_kind, 106.954943000000_fp_kind, 107.960515000000_fp_kind, 108.965131000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  40
  real(fp_kind), parameter :: isotopic_mass_040(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  76.966076000000_fp_kind,  77.956146000000_fp_kind,  78.949790000000_fp_kind,&
        79.941213000000_fp_kind,  80.938245000000_fp_kind,  81.931707497000_fp_kind,  82.929240926000_fp_kind,&
        83.923325663000_fp_kind,  84.921443199000_fp_kind,  85.916296814000_fp_kind,  86.914817338000_fp_kind,&
        87.910220715000_fp_kind,  88.908879751000_fp_kind,  89.904698755000_fp_kind,  90.905640205000_fp_kind,&
        91.905035336000_fp_kind,  92.906470661000_fp_kind,  93.906312523000_fp_kind,  94.908040276000_fp_kind,&
        95.908277615000_fp_kind,  96.910963802000_fp_kind,  97.912740448000_fp_kind,  98.916675081000_fp_kind,&
        99.918010499000_fp_kind, 100.921458454000_fp_kind, 101.923154181000_fp_kind, 102.927204054000_fp_kind,&
       103.929449193000_fp_kind, 104.934021832000_fp_kind, 105.936930000000_fp_kind, 106.942007000000_fp_kind,&
       107.945303000000_fp_kind, 108.950907000000_fp_kind, 109.954675000000_fp_kind, 110.960837000000_fp_kind,&
       111.965196000000_fp_kind, 112.971723000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  41
  real(fp_kind), parameter :: isotopic_mass_041(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  78.966022000000_fp_kind,  79.958754000000_fp_kind,  80.950230000000_fp_kind,&
        81.944380000000_fp_kind,  82.938150000000_fp_kind,  83.934305711000_fp_kind,  84.928845836000_fp_kind,&
        85.925781536000_fp_kind,  86.920692473000_fp_kind,  87.918226476000_fp_kind,  88.913444696000_fp_kind,&
        89.911259201000_fp_kind,  90.906990256000_fp_kind,  91.907188580000_fp_kind,  92.906373170000_fp_kind,&
        93.907279001000_fp_kind,  94.906831110000_fp_kind,  95.908101586000_fp_kind,  96.908101622000_fp_kind,&
        97.910332645000_fp_kind,  98.911609377000_fp_kind,  99.914340578000_fp_kind, 100.915306508000_fp_kind,&
       101.918090447000_fp_kind, 102.919453416000_fp_kind, 103.922907728000_fp_kind, 104.924942577000_fp_kind,&
       105.928928505000_fp_kind, 106.931589685000_fp_kind, 107.936075604000_fp_kind, 108.939141000000_fp_kind,&
       109.943843000000_fp_kind, 110.947439000000_fp_kind, 111.952689000000_fp_kind, 112.956833000000_fp_kind,&
       113.962469000000_fp_kind, 114.966849000000_fp_kind, 115.972914000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  42
  real(fp_kind), parameter :: isotopic_mass_042(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  80.966226000000_fp_kind,  81.956661000000_fp_kind,  82.950252000000_fp_kind,&
        83.941846000000_fp_kind,  84.938260736000_fp_kind,  85.931174092000_fp_kind,  86.928196198000_fp_kind,&
        87.921967779000_fp_kind,  88.919468149000_fp_kind,  89.913931270000_fp_kind,  90.911745190000_fp_kind,&
        91.906807153000_fp_kind,  92.906808772000_fp_kind,  93.905083586000_fp_kind,  94.905837436000_fp_kind,&
        95.904674770000_fp_kind,  96.906016903000_fp_kind,  97.905403609000_fp_kind,  98.907707299000_fp_kind,&
        99.907467982000_fp_kind, 100.910337648000_fp_kind, 101.910293725000_fp_kind, 102.913091954000_fp_kind,&
       103.913747443000_fp_kind, 104.916981989000_fp_kind, 105.918273231000_fp_kind, 106.922119770000_fp_kind,&
       107.924047508000_fp_kind, 108.928438318000_fp_kind, 109.930717956000_fp_kind, 110.935651966000_fp_kind,&
       111.938293000000_fp_kind, 112.943478000000_fp_kind, 113.946666000000_fp_kind, 114.952174000000_fp_kind,&
       115.955759000000_fp_kind, 116.961686000000_fp_kind, 117.965249000000_fp_kind, 118.971465000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  43
  real(fp_kind), parameter :: isotopic_mass_043(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  82.966377000000_fp_kind,  83.959527000000_fp_kind,  84.950778000000_fp_kind,&
        85.944637000000_fp_kind,  86.938067185000_fp_kind,  87.933794211000_fp_kind,  88.927648649000_fp_kind,&
        89.924073919000_fp_kind,  90.918424972000_fp_kind,  91.915269777000_fp_kind,  92.910245147000_fp_kind,&
        93.909652319000_fp_kind,  94.907652281000_fp_kind,  95.907866675000_fp_kind,  96.906360720000_fp_kind,&
        97.907211206000_fp_kind,  98.906249681000_fp_kind,  99.907652715000_fp_kind, 100.907305271000_fp_kind,&
       101.909207239000_fp_kind, 102.909173960000_fp_kind, 103.911433718000_fp_kind, 104.911662024000_fp_kind,&
       105.914356674000_fp_kind, 106.915458437000_fp_kind, 107.918493493000_fp_kind, 108.920254107000_fp_kind,&
       109.923741263000_fp_kind, 110.925898966000_fp_kind, 111.929941658000_fp_kind, 112.932569032000_fp_kind,&
       113.937090000000_fp_kind, 114.940100000000_fp_kind, 115.945020000000_fp_kind, 116.948320000000_fp_kind,&
       117.953526000000_fp_kind, 118.956876000000_fp_kind, 119.962426000000_fp_kind, 120.966140000000_fp_kind,&
       121.971760000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  44
  real(fp_kind), parameter :: isotopic_mass_044(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,  84.967117000000_fp_kind,  85.957305000000_fp_kind,  86.950907000000_fp_kind,&
        87.941664000000_fp_kind,  88.937337849000_fp_kind,  89.930344378000_fp_kind,  90.926741530000_fp_kind,&
        91.920234373000_fp_kind,  92.917104442000_fp_kind,  93.911342860000_fp_kind,  94.910404415000_fp_kind,&
        95.907588910000_fp_kind,  96.907545776000_fp_kind,  97.905286709000_fp_kind,  98.905930284000_fp_kind,&
        99.904210460000_fp_kind, 100.905573086000_fp_kind, 101.904340312000_fp_kind, 102.906314846000_fp_kind,&
       103.905425312000_fp_kind, 104.907745478000_fp_kind, 105.907328181000_fp_kind, 106.909969837000_fp_kind,&
       107.910185793000_fp_kind, 108.913323707000_fp_kind, 109.914038501000_fp_kind, 110.917567566000_fp_kind,&
       111.918806922000_fp_kind, 112.922846729000_fp_kind, 113.924614430000_fp_kind, 114.929033049000_fp_kind,&
       115.931219191000_fp_kind, 116.936135000000_fp_kind, 117.938808000000_fp_kind, 118.944090000000_fp_kind,&
       119.946623000000_fp_kind, 120.952098000000_fp_kind, 121.955147000000_fp_kind, 122.960762000000_fp_kind,&
       123.963940000000_fp_kind, 124.969544000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  45
  real(fp_kind), parameter :: isotopic_mass_045(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  87.960429000000_fp_kind,  88.950992000000_fp_kind,&
        89.944569000000_fp_kind,  90.937123000000_fp_kind,  91.932367692000_fp_kind,  92.925912778000_fp_kind,&
        93.921730450000_fp_kind,  94.915897893000_fp_kind,  95.914451705000_fp_kind,  96.911327872000_fp_kind,&
        97.910707734000_fp_kind,  98.908121241000_fp_kind,  99.908114147000_fp_kind, 100.906158903000_fp_kind,&
       101.906834282000_fp_kind, 102.905494081000_fp_kind, 103.906645309000_fp_kind, 104.905687787000_fp_kind,&
       105.907285879000_fp_kind, 106.906747975000_fp_kind, 107.908715304000_fp_kind, 108.908749555000_fp_kind,&
       109.911079745000_fp_kind, 110.911643164000_fp_kind, 111.914405199000_fp_kind, 112.915440212000_fp_kind,&
       113.918721680000_fp_kind, 114.920311649000_fp_kind, 115.924062060000_fp_kind, 116.926036291000_fp_kind,&
       117.930341116000_fp_kind, 118.932556951000_fp_kind, 119.937069000000_fp_kind, 120.939613000000_fp_kind,&
       121.944305000000_fp_kind, 122.947192000000_fp_kind, 123.952002000000_fp_kind, 124.955094000000_fp_kind,&
       125.960064000000_fp_kind, 126.963789000000_fp_kind, 127.970649000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  46
  real(fp_kind), parameter :: isotopic_mass_046(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  89.957370000000_fp_kind,  90.950435000000_fp_kind,&
        91.941192225000_fp_kind,  92.936680426000_fp_kind,  93.929036286000_fp_kind,  94.924888506000_fp_kind,&
        95.918213739000_fp_kind,  96.916471985000_fp_kind,  97.912698335000_fp_kind,  98.911773073000_fp_kind,&
        99.908520438000_fp_kind, 100.908284824000_fp_kind, 101.905632292000_fp_kind, 102.906111074000_fp_kind,&
       103.904030393000_fp_kind, 104.905079479000_fp_kind, 105.903480287000_fp_kind, 106.905128058000_fp_kind,&
       107.903891806000_fp_kind, 108.905950576000_fp_kind, 109.905172878000_fp_kind, 110.907690358000_fp_kind,&
       111.907330557000_fp_kind, 112.910261912000_fp_kind, 113.910369430000_fp_kind, 114.913659333000_fp_kind,&
       115.914297872000_fp_kind, 116.917955584000_fp_kind, 117.919067273000_fp_kind, 118.923341138000_fp_kind,&
       119.924551745000_fp_kind, 120.928950342000_fp_kind, 121.930631693000_fp_kind, 122.935126000000_fp_kind,&
       123.937305000000_fp_kind, 124.942072000000_fp_kind, 125.944401000000_fp_kind, 126.949307000000_fp_kind,&
       127.952345000000_fp_kind, 128.959334000000_fp_kind, 129.964863000000_fp_kind, 130.972367000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  47
  real(fp_kind), parameter :: isotopic_mass_047(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  91.959710000000_fp_kind,  92.950188000000_fp_kind,&
        93.943744000000_fp_kind,  94.935688000000_fp_kind,  95.930743903000_fp_kind,  96.923881400000_fp_kind,&
        97.921559970000_fp_kind,  98.917645766000_fp_kind,  99.916115443000_fp_kind, 100.912683951000_fp_kind,&
       101.911704538000_fp_kind, 102.908960558000_fp_kind, 103.908623715000_fp_kind, 104.906525604000_fp_kind,&
       105.906663499000_fp_kind, 106.905091509000_fp_kind, 107.905950245000_fp_kind, 108.904755778000_fp_kind,&
       109.906110724000_fp_kind, 110.905296827000_fp_kind, 111.907048548000_fp_kind, 112.906572865000_fp_kind,&
       113.908823029000_fp_kind, 114.908767445000_fp_kind, 115.911386809000_fp_kind, 116.911774086000_fp_kind,&
       117.914595484000_fp_kind, 118.915570309000_fp_kind, 119.918784765000_fp_kind, 120.920125279000_fp_kind,&
       121.923664446000_fp_kind, 122.925315060000_fp_kind, 123.928899227000_fp_kind, 124.930735000000_fp_kind,&
       125.934814000000_fp_kind, 126.937037000000_fp_kind, 127.941266000000_fp_kind, 128.944315000000_fp_kind,&
       129.950727000000_fp_kind, 130.956253000000_fp_kind, 131.963070000000_fp_kind, 132.968781000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  48
  real(fp_kind), parameter :: isotopic_mass_048(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  93.956586000000_fp_kind,  94.949483000000_fp_kind,&
        95.940341000000_fp_kind,  96.934799343000_fp_kind,  97.927389315000_fp_kind,  98.924925845000_fp_kind,&
        99.920348829000_fp_kind, 100.918586209000_fp_kind, 101.914481797000_fp_kind, 102.913416922000_fp_kind,&
       103.909856228000_fp_kind, 104.909463893000_fp_kind, 105.906459791000_fp_kind, 106.906612049000_fp_kind,&
       107.904183588000_fp_kind, 108.904986697000_fp_kind, 109.903007470000_fp_kind, 110.904183776000_fp_kind,&
       111.902763896000_fp_kind, 112.904408105000_fp_kind, 113.903364998000_fp_kind, 114.905437426000_fp_kind,&
       115.904763230000_fp_kind, 116.907226039000_fp_kind, 117.906921956000_fp_kind, 118.909847052000_fp_kind,&
       119.909868065000_fp_kind, 120.912963660000_fp_kind, 121.913459050000_fp_kind, 122.916892460000_fp_kind,&
       123.917659772000_fp_kind, 124.921257590000_fp_kind, 125.922430290000_fp_kind, 126.926203291000_fp_kind,&
       127.927816778000_fp_kind, 128.932235597000_fp_kind, 129.934387563000_fp_kind, 130.940727740000_fp_kind,&
       131.945823136000_fp_kind, 132.952614000000_fp_kind, 133.957638000000_fp_kind, 134.964766000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  49
  real(fp_kind), parameter :: isotopic_mass_049(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,  95.959109000000_fp_kind,  96.949125000000_fp_kind,&
        97.942129000000_fp_kind,  98.934110000000_fp_kind,  99.931101929000_fp_kind, 100.926414025000_fp_kind,&
       101.924105911000_fp_kind, 102.919878830000_fp_kind, 103.918214538000_fp_kind, 104.914502322000_fp_kind,&
       105.913463596000_fp_kind, 106.910287497000_fp_kind, 107.909693654000_fp_kind, 108.907149679000_fp_kind,&
       109.907170674000_fp_kind, 110.905107236000_fp_kind, 111.905538718000_fp_kind, 112.904060451000_fp_kind,&
       113.904916405000_fp_kind, 114.903878772000_fp_kind, 115.905259992000_fp_kind, 116.904515729000_fp_kind,&
       117.906356705000_fp_kind, 118.905851622000_fp_kind, 119.907967489000_fp_kind, 120.907852778000_fp_kind,&
       121.910282458000_fp_kind, 122.910435252000_fp_kind, 123.913184873000_fp_kind, 124.913673841000_fp_kind,&
       125.916468202000_fp_kind, 126.917466040000_fp_kind, 127.920353637000_fp_kind, 128.921808534000_fp_kind,&
       129.924952257000_fp_kind, 130.926972839000_fp_kind, 131.932998444000_fp_kind, 132.938067000000_fp_kind,&
       133.944208000000_fp_kind, 134.949425000000_fp_kind, 135.956017000000_fp_kind, 136.961535000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  50
  real(fp_kind), parameter :: isotopic_mass_050(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,  98.948495000000_fp_kind,&
        99.938648944000_fp_kind, 100.935259252000_fp_kind, 101.930289525000_fp_kind, 102.927973000000_fp_kind,&
       103.923105195000_fp_kind, 104.921268421000_fp_kind, 105.916957394000_fp_kind, 106.915713649000_fp_kind,&
       107.911894290000_fp_kind, 108.911292857000_fp_kind, 109.907844835000_fp_kind, 110.907741143000_fp_kind,&
       111.904824894000_fp_kind, 112.905175857000_fp_kind, 113.902780130000_fp_kind, 114.903344695000_fp_kind,&
       115.901742825000_fp_kind, 116.902954036000_fp_kind, 117.901606630000_fp_kind, 118.903311266000_fp_kind,&
       119.902202557000_fp_kind, 120.904243488000_fp_kind, 121.903445494000_fp_kind, 122.905727065000_fp_kind,&
       123.905279619000_fp_kind, 124.907789370000_fp_kind, 125.907658958000_fp_kind, 126.910391726000_fp_kind,&
       127.910507828000_fp_kind, 128.913482440000_fp_kind, 129.913974531000_fp_kind, 130.917053067000_fp_kind,&
       131.917823898000_fp_kind, 132.923913753000_fp_kind, 133.928680430000_fp_kind, 134.934908603000_fp_kind,&
       135.939699000000_fp_kind, 136.946162000000_fp_kind, 137.951143000000_fp_kind, 138.957799000000_fp_kind,&
       139.962973000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  51
  real(fp_kind), parameter :: isotopic_mass_051(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       101.945142000000_fp_kind, 102.939162000000_fp_kind, 103.936344000000_fp_kind, 104.931276547000_fp_kind,&
       105.928637979000_fp_kind, 106.924150621000_fp_kind, 107.922226731000_fp_kind, 108.918141203000_fp_kind,&
       109.916854283000_fp_kind, 110.913218187000_fp_kind, 111.912399903000_fp_kind, 112.909374664000_fp_kind,&
       113.909289155000_fp_kind, 114.906598000000_fp_kind, 115.906792732000_fp_kind, 116.904841519000_fp_kind,&
       117.905532194000_fp_kind, 118.903944062000_fp_kind, 119.905080308000_fp_kind, 120.903811353000_fp_kind,&
       121.905169335000_fp_kind, 122.904215292000_fp_kind, 123.905937065000_fp_kind, 124.905254264000_fp_kind,&
       125.907253158000_fp_kind, 126.906925557000_fp_kind, 127.909146121000_fp_kind, 128.909146623000_fp_kind,&
       129.911662686000_fp_kind, 130.911989339000_fp_kind, 131.914508013000_fp_kind, 132.915272128000_fp_kind,&
       133.920537334000_fp_kind, 134.925184354000_fp_kind, 135.930749009000_fp_kind, 136.935522519000_fp_kind,&
       137.941331000000_fp_kind, 138.946269000000_fp_kind, 139.952345000000_fp_kind, 140.957552000000_fp_kind,&
       141.963918000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  52
  real(fp_kind), parameter :: isotopic_mass_052(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       103.946723408000_fp_kind, 104.943304516000_fp_kind, 105.937498521000_fp_kind, 106.934882000000_fp_kind,&
       107.929380469000_fp_kind, 108.927304532000_fp_kind, 109.922458102000_fp_kind, 110.921000587000_fp_kind,&
       111.916727848000_fp_kind, 112.915891000000_fp_kind, 113.912087820000_fp_kind, 114.911902000000_fp_kind,&
       115.908465558000_fp_kind, 116.908646227000_fp_kind, 117.905860104000_fp_kind, 118.906405699000_fp_kind,&
       119.904065779000_fp_kind, 120.904945065000_fp_kind, 121.903044708000_fp_kind, 122.904271022000_fp_kind,&
       123.902818341000_fp_kind, 124.904431178000_fp_kind, 125.903312144000_fp_kind, 126.905226993000_fp_kind,&
       127.904461237000_fp_kind, 128.906596419000_fp_kind, 129.906222745000_fp_kind, 130.908522210000_fp_kind,&
       131.908546713000_fp_kind, 132.910963330000_fp_kind, 133.911396376000_fp_kind, 134.916554715000_fp_kind,&
       135.920101180000_fp_kind, 136.925599354000_fp_kind, 137.929472452000_fp_kind, 138.935367191000_fp_kind,&
       139.939487057000_fp_kind, 140.945604000000_fp_kind, 141.950027000000_fp_kind, 142.956489000000_fp_kind,&
       143.961116000000_fp_kind, 144.967783000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  53
  real(fp_kind), parameter :: isotopic_mass_053(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       105.953516000000_fp_kind, 106.946935000000_fp_kind, 107.943348000000_fp_kind, 108.938086022000_fp_kind,&
       109.935085102000_fp_kind, 110.930269236000_fp_kind, 111.928004548000_fp_kind, 112.923650062000_fp_kind,&
       113.922018900000_fp_kind, 114.918048000000_fp_kind, 115.916885513000_fp_kind, 116.913645649000_fp_kind,&
       117.913074000000_fp_kind, 118.910060910000_fp_kind, 119.910093729000_fp_kind, 120.907411492000_fp_kind,&
       121.907590094000_fp_kind, 122.905589753000_fp_kind, 123.906210297000_fp_kind, 124.904630610000_fp_kind,&
       125.905624205000_fp_kind, 126.904472592000_fp_kind, 127.905809355000_fp_kind, 128.904983643000_fp_kind,&
       129.906670168000_fp_kind, 130.906126375000_fp_kind, 131.907993511000_fp_kind, 132.907828400000_fp_kind,&
       133.909775660000_fp_kind, 134.910059355000_fp_kind, 135.914604693000_fp_kind, 136.918028178000_fp_kind,&
       137.922726392000_fp_kind, 138.926493400000_fp_kind, 139.931715914000_fp_kind, 140.935666081000_fp_kind,&
       141.941166595000_fp_kind, 142.945475000000_fp_kind, 143.951336000000_fp_kind, 144.955845000000_fp_kind,&
       145.961846000000_fp_kind, 146.966505000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  54
  real(fp_kind), parameter :: isotopic_mass_054(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       107.954232285000_fp_kind, 108.950434955000_fp_kind, 109.944258759000_fp_kind, 110.941470000000_fp_kind,&
       111.935559068000_fp_kind, 112.933221663000_fp_kind, 113.927980329000_fp_kind, 114.926293943000_fp_kind,&
       115.921580955000_fp_kind, 116.920358758000_fp_kind, 117.916178678000_fp_kind, 118.915410641000_fp_kind,&
       119.911784267000_fp_kind, 120.911453012000_fp_kind, 121.908367655000_fp_kind, 122.908482235000_fp_kind,&
       123.905885174000_fp_kind, 124.906387640000_fp_kind, 125.904297422000_fp_kind, 126.905183636000_fp_kind,&
       127.903530753410_fp_kind, 128.904780857420_fp_kind, 129.903509346000_fp_kind, 130.905084128080_fp_kind,&
       131.904155083460_fp_kind, 132.905910748000_fp_kind, 133.905393030000_fp_kind, 134.907231441000_fp_kind,&
       135.907214474000_fp_kind, 136.911557771000_fp_kind, 137.914146268000_fp_kind, 138.918792200000_fp_kind,&
       139.921645814000_fp_kind, 140.926787181000_fp_kind, 141.929973095000_fp_kind, 142.935369550000_fp_kind,&
       143.938945076000_fp_kind, 144.944719631000_fp_kind, 145.948518245000_fp_kind, 146.954482000000_fp_kind,&
       147.958508000000_fp_kind, 148.964573000000_fp_kind, 149.968878000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  55
  real(fp_kind), parameter :: isotopic_mass_055(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 110.953945000000_fp_kind, 111.950172000000_fp_kind, 112.944428484000_fp_kind,&
       113.941292244000_fp_kind, 114.935910000000_fp_kind, 115.933395000000_fp_kind, 116.928616723000_fp_kind,&
       117.926559517000_fp_kind, 118.922377327000_fp_kind, 119.920677277000_fp_kind, 120.917227235000_fp_kind,&
       121.916108144000_fp_kind, 122.912996060000_fp_kind, 123.912247366000_fp_kind, 124.909725953000_fp_kind,&
       125.909445821000_fp_kind, 126.907417527000_fp_kind, 127.907748452000_fp_kind, 128.906065910000_fp_kind,&
       129.906709281000_fp_kind, 130.905468457000_fp_kind, 131.906437740000_fp_kind, 132.905451958000_fp_kind,&
       133.906718501000_fp_kind, 134.905976907000_fp_kind, 135.907311431000_fp_kind, 136.907089296000_fp_kind,&
       137.911017119000_fp_kind, 138.913363822000_fp_kind, 139.917283707000_fp_kind, 140.920045279000_fp_kind,&
       141.924299514000_fp_kind, 142.927347346000_fp_kind, 143.932075402000_fp_kind, 144.935528927000_fp_kind,&
       145.940621867000_fp_kind, 146.944261512000_fp_kind, 147.949639026000_fp_kind, 148.953516000000_fp_kind,&
       149.959023000000_fp_kind, 150.963199000000_fp_kind, 151.968728000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  56
  real(fp_kind), parameter :: isotopic_mass_056(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 112.957370000000_fp_kind, 113.950718489000_fp_kind, 114.947482000000_fp_kind,&
       115.941621000000_fp_kind, 116.938316403000_fp_kind, 117.933226000000_fp_kind, 118.930659683000_fp_kind,&
       119.926044997000_fp_kind, 120.924052286000_fp_kind, 121.919904000000_fp_kind, 122.918781060000_fp_kind,&
       123.915093627000_fp_kind, 124.914471840000_fp_kind, 125.911250202000_fp_kind, 126.911091272000_fp_kind,&
       127.908352446000_fp_kind, 128.908683409000_fp_kind, 129.906326002000_fp_kind, 130.906946315000_fp_kind,&
       131.905061231000_fp_kind, 132.906007443000_fp_kind, 133.904508249000_fp_kind, 134.905688447000_fp_kind,&
       135.904575800000_fp_kind, 136.905827207000_fp_kind, 137.905247059000_fp_kind, 138.908841164000_fp_kind,&
       139.910608231000_fp_kind, 140.914403653000_fp_kind, 141.916432904000_fp_kind, 142.920625149000_fp_kind,&
       143.922954821000_fp_kind, 144.927518400000_fp_kind, 145.930363200000_fp_kind, 146.935303900000_fp_kind,&
       147.938223000000_fp_kind, 148.943284000000_fp_kind, 149.946441100000_fp_kind, 150.951755000000_fp_kind,&
       151.955330000000_fp_kind, 152.960848000000_fp_kind, 153.964659000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  57
  real(fp_kind), parameter :: isotopic_mass_057(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 115.957005000000_fp_kind, 116.950326000000_fp_kind,&
       117.946731000000_fp_kind, 118.940934000000_fp_kind, 119.938196000000_fp_kind, 120.933236000000_fp_kind,&
       121.930710000000_fp_kind, 122.926300000000_fp_kind, 123.924574275000_fp_kind, 124.920815931000_fp_kind,&
       125.919512667000_fp_kind, 126.916375083000_fp_kind, 127.915592123000_fp_kind, 128.912695592000_fp_kind,&
       129.912369413000_fp_kind, 130.910070000000_fp_kind, 131.910119047000_fp_kind, 132.908218000000_fp_kind,&
       133.908514011000_fp_kind, 134.906984427000_fp_kind, 135.907634962000_fp_kind, 136.906450438000_fp_kind,&
       137.907124041000_fp_kind, 138.906362927000_fp_kind, 139.909487285000_fp_kind, 140.910971155000_fp_kind,&
       141.914090760000_fp_kind, 142.916079482000_fp_kind, 143.919645589000_fp_kind, 144.921808065000_fp_kind,&
       145.925688017000_fp_kind, 146.928417800000_fp_kind, 147.932679400000_fp_kind, 148.935351259000_fp_kind,&
       149.939547500000_fp_kind, 150.942769000000_fp_kind, 151.947085000000_fp_kind, 152.950553000000_fp_kind,&
       153.955416000000_fp_kind, 154.959280000000_fp_kind, 155.964519000000_fp_kind, 156.968792000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  58
  real(fp_kind), parameter :: isotopic_mass_058(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 118.952957000000_fp_kind,&
       119.946613000000_fp_kind, 120.943435000000_fp_kind, 121.937870000000_fp_kind, 122.935280000000_fp_kind,&
       123.930310000000_fp_kind, 124.928440000000_fp_kind, 125.923971000000_fp_kind, 126.922727000000_fp_kind,&
       127.918911000000_fp_kind, 128.918102000000_fp_kind, 129.914736000000_fp_kind, 130.914429465000_fp_kind,&
       131.911466226000_fp_kind, 132.911520402000_fp_kind, 133.908928142000_fp_kind, 134.909160662000_fp_kind,&
       135.907129256000_fp_kind, 136.907762416000_fp_kind, 137.905994180000_fp_kind, 138.906647029000_fp_kind,&
       139.905448433000_fp_kind, 140.908285991000_fp_kind, 141.909250208000_fp_kind, 142.912391953000_fp_kind,&
       143.913652763000_fp_kind, 144.917265113000_fp_kind, 145.918812294000_fp_kind, 146.922689900000_fp_kind,&
       147.924424186000_fp_kind, 148.928426900000_fp_kind, 149.930384032000_fp_kind, 150.934272200000_fp_kind,&
       151.936682000000_fp_kind, 152.941052000000_fp_kind, 153.943940000000_fp_kind, 154.948706000000_fp_kind,&
       155.951884000000_fp_kind, 156.957133000000_fp_kind, 157.960773000000_fp_kind, 158.966355000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  59
  real(fp_kind), parameter :: isotopic_mass_059(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 120.955393000000_fp_kind,&
       121.951927000000_fp_kind, 122.946076000000_fp_kind, 123.942940000000_fp_kind, 124.937659000000_fp_kind,&
       125.935240000000_fp_kind, 126.930710000000_fp_kind, 127.928791000000_fp_kind, 128.925095000000_fp_kind,&
       129.923590000000_fp_kind, 130.920234960000_fp_kind, 131.919240000000_fp_kind, 132.916330558000_fp_kind,&
       133.915696729000_fp_kind, 134.913111772000_fp_kind, 135.912677470000_fp_kind, 136.910679183000_fp_kind,&
       137.910757495000_fp_kind, 138.908932700000_fp_kind, 139.909085600000_fp_kind, 140.907659604000_fp_kind,&
       141.910051640000_fp_kind, 142.910822624000_fp_kind, 143.913310682000_fp_kind, 144.914517987000_fp_kind,&
       145.917687630000_fp_kind, 146.919007438000_fp_kind, 147.922129992000_fp_kind, 148.923736100000_fp_kind,&
       149.926676391000_fp_kind, 150.928309066000_fp_kind, 151.931552900000_fp_kind, 152.933903511000_fp_kind,&
       153.937885165000_fp_kind, 154.940509193000_fp_kind, 155.944766900000_fp_kind, 156.948003100000_fp_kind,&
       157.952603000000_fp_kind, 158.956232000000_fp_kind, 159.961138000000_fp_kind, 160.965121000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  60
  real(fp_kind), parameter :: isotopic_mass_060(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       123.951873000000_fp_kind, 124.948395000000_fp_kind, 125.942694000000_fp_kind, 126.939978000000_fp_kind,&
       127.935018000000_fp_kind, 128.933038000000_fp_kind, 129.928506000000_fp_kind, 130.927248020000_fp_kind,&
       131.923321237000_fp_kind, 132.922348000000_fp_kind, 133.918790207000_fp_kind, 134.918181318000_fp_kind,&
       135.914976061000_fp_kind, 136.914563099000_fp_kind, 137.911950938000_fp_kind, 138.911951208000_fp_kind,&
       139.909546130000_fp_kind, 140.909616690000_fp_kind, 141.907728824000_fp_kind, 142.909819815000_fp_kind,&
       143.910092798000_fp_kind, 144.912579151000_fp_kind, 145.913122459000_fp_kind, 146.916105969000_fp_kind,&
       147.916899027000_fp_kind, 148.920154583000_fp_kind, 149.920901322000_fp_kind, 150.923839363000_fp_kind,&
       151.924691242000_fp_kind, 152.927717868000_fp_kind, 153.929597404000_fp_kind, 154.933135598000_fp_kind,&
       155.935370358000_fp_kind, 156.939351074000_fp_kind, 157.942205620000_fp_kind, 158.946619085000_fp_kind,&
       159.949839172000_fp_kind, 160.954664000000_fp_kind, 161.958121000000_fp_kind, 162.963414000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  61
  real(fp_kind), parameter :: isotopic_mass_061(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       125.957327000000_fp_kind, 126.951358000000_fp_kind, 127.948234000000_fp_kind, 128.942909000000_fp_kind,&
       129.940451000000_fp_kind, 130.935834000000_fp_kind, 131.933840000000_fp_kind, 132.929782000000_fp_kind,&
       133.928326000000_fp_kind, 134.924785000000_fp_kind, 135.923595949000_fp_kind, 136.920479519000_fp_kind,&
       137.919576119000_fp_kind, 138.916799228000_fp_kind, 139.916035918000_fp_kind, 140.913555081000_fp_kind,&
       141.912890982000_fp_kind, 142.910938068000_fp_kind, 143.912596208000_fp_kind, 144.912755748000_fp_kind,&
       145.914702240000_fp_kind, 146.915144944000_fp_kind, 147.917481091000_fp_kind, 148.918341507000_fp_kind,&
       149.920990014000_fp_kind, 150.921216613000_fp_kind, 151.923505185000_fp_kind, 152.924156252000_fp_kind,&
       153.926712791000_fp_kind, 154.928136951000_fp_kind, 155.931114059000_fp_kind, 156.933121298000_fp_kind,&
       157.936546948000_fp_kind, 158.939286409000_fp_kind, 159.943215272000_fp_kind, 160.946229837000_fp_kind,&
       161.950574000000_fp_kind, 162.953881000000_fp_kind, 163.958819000000_fp_kind, 164.962780000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  62
  real(fp_kind), parameter :: isotopic_mass_062(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       127.957971000000_fp_kind, 128.954557000000_fp_kind, 129.948792000000_fp_kind, 130.946022000000_fp_kind,&
       131.940805000000_fp_kind, 132.938560000000_fp_kind, 133.934110000000_fp_kind, 134.932520000000_fp_kind,&
       135.928275553000_fp_kind, 136.927007959000_fp_kind, 137.923243988000_fp_kind, 138.922296631000_fp_kind,&
       139.918994714000_fp_kind, 140.918481545000_fp_kind, 141.915209415000_fp_kind, 142.914634848000_fp_kind,&
       143.912006285000_fp_kind, 144.913417157000_fp_kind, 145.913046835000_fp_kind, 146.914904401000_fp_kind,&
       147.914829233000_fp_kind, 148.917191211000_fp_kind, 149.917281993000_fp_kind, 150.919938859000_fp_kind,&
       151.919738646000_fp_kind, 152.922103576000_fp_kind, 153.922215756000_fp_kind, 154.924646645000_fp_kind,&
       155.925538191000_fp_kind, 156.928418598000_fp_kind, 157.929949262000_fp_kind, 158.933217130000_fp_kind,&
       159.935337032000_fp_kind, 160.939160062000_fp_kind, 161.941621687000_fp_kind, 162.945679085000_fp_kind,&
       163.948550061000_fp_kind, 164.953290000000_fp_kind, 165.956575000000_fp_kind, 166.962072000000_fp_kind,&
       167.966033000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  63
  real(fp_kind), parameter :: isotopic_mass_063(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       129.964022000000_fp_kind, 130.957634000000_fp_kind, 131.954696000000_fp_kind, 132.949290000000_fp_kind,&
       133.946537000000_fp_kind, 134.941870000000_fp_kind, 135.939620000000_fp_kind, 136.935430719000_fp_kind,&
       137.933709000000_fp_kind, 138.929792307000_fp_kind, 139.928087633000_fp_kind, 140.924931734000_fp_kind,&
       141.923446719000_fp_kind, 142.920298678000_fp_kind, 143.918819481000_fp_kind, 144.916272659000_fp_kind,&
       145.917210852000_fp_kind, 146.916752440000_fp_kind, 147.918091288000_fp_kind, 148.917936875000_fp_kind,&
       149.919707092000_fp_kind, 150.919856606000_fp_kind, 151.921750980000_fp_kind, 152.921236789000_fp_kind,&
       153.922985699000_fp_kind, 154.922899847000_fp_kind, 155.924762976000_fp_kind, 156.925432556000_fp_kind,&
       157.927782192000_fp_kind, 158.929099512000_fp_kind, 159.931836982000_fp_kind, 160.933663991000_fp_kind,&
       161.936958329000_fp_kind, 162.939265510000_fp_kind, 163.942852943000_fp_kind, 164.945540070000_fp_kind,&
       165.949813000000_fp_kind, 166.953011000000_fp_kind, 167.957863000000_fp_kind, 168.961717000000_fp_kind,&
       169.966870000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  64
  real(fp_kind), parameter :: isotopic_mass_064(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 132.961288000000_fp_kind, 133.955416000000_fp_kind, 134.952496000000_fp_kind,&
       135.947300000000_fp_kind, 136.945020000000_fp_kind, 137.940247000000_fp_kind, 138.938130000000_fp_kind,&
       139.933674000000_fp_kind, 140.932126000000_fp_kind, 141.928116000000_fp_kind, 142.926750678000_fp_kind,&
       143.922963000000_fp_kind, 144.921710051000_fp_kind, 145.918318513000_fp_kind, 146.919101014000_fp_kind,&
       147.918121414000_fp_kind, 148.919347666000_fp_kind, 149.918663949000_fp_kind, 150.920354922000_fp_kind,&
       151.919798414000_fp_kind, 152.921756945000_fp_kind, 153.920872974000_fp_kind, 154.922629356000_fp_kind,&
       155.922130120000_fp_kind, 156.923967424000_fp_kind, 157.924111200000_fp_kind, 158.926395822000_fp_kind,&
       159.927061202000_fp_kind, 160.929676267000_fp_kind, 161.930991812000_fp_kind, 162.934096640000_fp_kind,&
       163.935916193000_fp_kind, 164.939317080000_fp_kind, 165.941630413000_fp_kind, 166.945490012000_fp_kind,&
       167.948309000000_fp_kind, 168.952882000000_fp_kind, 169.956146000000_fp_kind, 170.961127000000_fp_kind,&
       171.964605000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  65
  real(fp_kind), parameter :: isotopic_mass_065(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 134.964516000000_fp_kind, 135.961460000000_fp_kind, 136.956020000000_fp_kind,&
       137.953193000000_fp_kind, 138.948330000000_fp_kind, 139.945805048000_fp_kind, 140.941448000000_fp_kind,&
       141.939280858000_fp_kind, 142.935137332000_fp_kind, 143.933045000000_fp_kind, 144.928717001000_fp_kind,&
       145.927252739000_fp_kind, 146.924054620000_fp_kind, 147.924275476000_fp_kind, 148.923253792000_fp_kind,&
       149.923664799000_fp_kind, 150.923108970000_fp_kind, 151.924081855000_fp_kind, 152.923441694000_fp_kind,&
       153.924683681000_fp_kind, 154.923509511000_fp_kind, 155.924754209000_fp_kind, 156.924031888000_fp_kind,&
       157.925419942000_fp_kind, 158.925353707000_fp_kind, 159.927174553000_fp_kind, 160.927576806000_fp_kind,&
       161.929275400000_fp_kind, 162.930653609000_fp_kind, 163.933327561000_fp_kind, 164.934955198000_fp_kind,&
       165.937939727000_fp_kind, 166.940007046000_fp_kind, 167.943337074000_fp_kind, 168.945807000000_fp_kind,&
       169.949855000000_fp_kind, 170.953011000000_fp_kind, 171.957391000000_fp_kind, 172.960805000000_fp_kind,&
       173.965679000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  66
  real(fp_kind), parameter :: isotopic_mass_066(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 137.962500000000_fp_kind, 138.959527000000_fp_kind,&
       139.954020000000_fp_kind, 140.951280000000_fp_kind, 141.946194000000_fp_kind, 142.943994332000_fp_kind,&
       143.939269512000_fp_kind, 144.937473992000_fp_kind, 145.932844526000_fp_kind, 146.931082712000_fp_kind,&
       147.927149944000_fp_kind, 148.927327516000_fp_kind, 149.925593068000_fp_kind, 150.926191279000_fp_kind,&
       151.924725274000_fp_kind, 152.925771729000_fp_kind, 153.924428920000_fp_kind, 154.925758049000_fp_kind,&
       155.924283593000_fp_kind, 156.925469555000_fp_kind, 157.924414817000_fp_kind, 158.925745938000_fp_kind,&
       159.925203578000_fp_kind, 160.926939425000_fp_kind, 161.926804507000_fp_kind, 162.928737221000_fp_kind,&
       163.929180819000_fp_kind, 164.931709402000_fp_kind, 165.932812810000_fp_kind, 166.935682415000_fp_kind,&
       167.937134977000_fp_kind, 168.940315231000_fp_kind, 169.942340000000_fp_kind, 170.946312000000_fp_kind,&
       171.948728000000_fp_kind, 172.953043000000_fp_kind, 173.955845000000_fp_kind, 174.960569000000_fp_kind,&
       175.963918000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  67
  real(fp_kind), parameter :: isotopic_mass_067(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 139.968526000000_fp_kind, 140.963108000000_fp_kind,&
       141.960010000000_fp_kind, 142.954860000000_fp_kind, 143.952109712000_fp_kind, 144.947267392000_fp_kind,&
       145.944993503000_fp_kind, 146.940142293000_fp_kind, 147.937743925000_fp_kind, 148.933820457000_fp_kind,&
       149.933498353000_fp_kind, 150.931698176000_fp_kind, 151.931717618000_fp_kind, 152.930206671000_fp_kind,&
       153.930606776000_fp_kind, 154.929103363000_fp_kind, 155.929641634000_fp_kind, 156.928251974000_fp_kind,&
       157.928944910000_fp_kind, 158.927718683000_fp_kind, 159.928735538000_fp_kind, 160.927861815000_fp_kind,&
       161.929102543000_fp_kind, 162.928740260000_fp_kind, 163.930240548000_fp_kind, 164.930329116000_fp_kind,&
       165.932291209000_fp_kind, 166.933140254000_fp_kind, 167.935523766000_fp_kind, 168.936879890000_fp_kind,&
       169.939626548000_fp_kind, 170.941472713000_fp_kind, 171.944730000000_fp_kind, 172.947020000000_fp_kind,&
       173.950757000000_fp_kind, 174.953516000000_fp_kind, 175.957713000000_fp_kind, 176.961052000000_fp_kind,&
       177.965507000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  68
  real(fp_kind), parameter :: isotopic_mass_068(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 141.970016000000_fp_kind, 142.966548000000_fp_kind,&
       143.960700000000_fp_kind, 144.957874000000_fp_kind, 145.952418357000_fp_kind, 146.949964456000_fp_kind,&
       147.944735026000_fp_kind, 148.942306000000_fp_kind, 149.937915524000_fp_kind, 150.937448567000_fp_kind,&
       151.935050347000_fp_kind, 152.935086350000_fp_kind, 153.932790799000_fp_kind, 154.933215710000_fp_kind,&
       155.931065926000_fp_kind, 156.931922652000_fp_kind, 157.929893474000_fp_kind, 158.930690790000_fp_kind,&
       159.929077193000_fp_kind, 160.930003530000_fp_kind, 161.928787299000_fp_kind, 162.930039908000_fp_kind,&
       163.929207739000_fp_kind, 164.930733482000_fp_kind, 165.930301067000_fp_kind, 166.932056192000_fp_kind,&
       167.932378282000_fp_kind, 168.934598444000_fp_kind, 169.935471933000_fp_kind, 170.938037372000_fp_kind,&
       171.939363461000_fp_kind, 172.942400000000_fp_kind, 173.944230000000_fp_kind, 174.947770000000_fp_kind,&
       175.949940000000_fp_kind, 176.953990000000_fp_kind, 177.956779000000_fp_kind, 178.961267000000_fp_kind,&
       179.964380000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  69
  real(fp_kind), parameter :: isotopic_mass_069(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 143.976211000000_fp_kind, 144.970389000000_fp_kind,&
       145.966661000000_fp_kind, 146.961379887000_fp_kind, 147.958384026000_fp_kind, 148.952828000000_fp_kind,&
       149.950090000000_fp_kind, 150.945494433000_fp_kind, 151.944476000000_fp_kind, 152.942058023000_fp_kind,&
       153.941570062000_fp_kind, 154.939209576000_fp_kind, 155.938985746000_fp_kind, 156.936973000000_fp_kind,&
       157.936979525000_fp_kind, 158.934975000000_fp_kind, 159.935264177000_fp_kind, 160.933549000000_fp_kind,&
       161.934001211000_fp_kind, 162.932658282000_fp_kind, 163.933538019000_fp_kind, 164.932441843000_fp_kind,&
       165.933562136000_fp_kind, 166.932857206000_fp_kind, 167.934178457000_fp_kind, 168.934218956000_fp_kind,&
       169.935807093000_fp_kind, 170.936435162000_fp_kind, 171.938406959000_fp_kind, 172.939606630000_fp_kind,&
       173.942174061000_fp_kind, 174.943842310000_fp_kind, 175.946997707000_fp_kind, 176.948932000000_fp_kind,&
       177.952506000000_fp_kind, 178.955018000000_fp_kind, 179.959023000000_fp_kind, 180.961954000000_fp_kind,&
       181.966194000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  70
  real(fp_kind), parameter :: isotopic_mass_070(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       147.967547000000_fp_kind, 148.964219000000_fp_kind, 149.958314000000_fp_kind, 150.955402453000_fp_kind,&
       151.950326699000_fp_kind, 152.949372000000_fp_kind, 153.946395696000_fp_kind, 154.945783216000_fp_kind,&
       155.942817096000_fp_kind, 156.942651368000_fp_kind, 157.939871202000_fp_kind, 158.940060257000_fp_kind,&
       159.937559210000_fp_kind, 160.937912384000_fp_kind, 161.935779342000_fp_kind, 162.936345406000_fp_kind,&
       163.934500743000_fp_kind, 164.935270241000_fp_kind, 165.933876439000_fp_kind, 166.934954069000_fp_kind,&
       167.933891297000_fp_kind, 168.935184208000_fp_kind, 169.934767242000_fp_kind, 170.936331515000_fp_kind,&
       171.936386654000_fp_kind, 172.938216211000_fp_kind, 173.938867545000_fp_kind, 174.941281907000_fp_kind,&
       175.942574706000_fp_kind, 176.945263846000_fp_kind, 177.946669400000_fp_kind, 178.949930000000_fp_kind,&
       179.951991000000_fp_kind, 180.955890000000_fp_kind, 181.958239000000_fp_kind, 182.962426000000_fp_kind,&
       183.965002000000_fp_kind, 184.969425000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  71
  real(fp_kind), parameter :: isotopic_mass_071(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       149.973407000000_fp_kind, 150.967471000000_fp_kind, 151.964120000000_fp_kind, 152.958802248000_fp_kind,&
       153.957416000000_fp_kind, 154.954326005000_fp_kind, 155.953086606000_fp_kind, 156.950144807000_fp_kind,&
       157.949315620000_fp_kind, 158.946635615000_fp_kind, 159.946033000000_fp_kind, 160.943572000000_fp_kind,&
       161.943282776000_fp_kind, 162.941179000000_fp_kind, 163.941339000000_fp_kind, 164.939406758000_fp_kind,&
       165.939859000000_fp_kind, 166.938243000000_fp_kind, 167.938729798000_fp_kind, 168.937645845000_fp_kind,&
       169.938479230000_fp_kind, 170.937918591000_fp_kind, 171.939091320000_fp_kind, 172.938935722000_fp_kind,&
       173.940342840000_fp_kind, 174.940777211000_fp_kind, 175.942691711000_fp_kind, 176.943763570000_fp_kind,&
       177.945960065000_fp_kind, 178.947332985000_fp_kind, 179.949890744000_fp_kind, 180.951908000000_fp_kind,&
       181.955158000000_fp_kind, 182.957363000000_fp_kind, 183.961030000000_fp_kind, 184.963542000000_fp_kind,&
       185.967450000000_fp_kind, 186.970188000000_fp_kind, 187.974428000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  72
  real(fp_kind), parameter :: isotopic_mass_072(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 152.970692000000_fp_kind, 153.964863000000_fp_kind, 154.963167000000_fp_kind,&
       155.959399083000_fp_kind, 156.958288000000_fp_kind, 157.954801217000_fp_kind, 158.953995837000_fp_kind,&
       159.950682728000_fp_kind, 160.950277927000_fp_kind, 161.947215526000_fp_kind, 162.947107211000_fp_kind,&
       163.944370709000_fp_kind, 164.944567000000_fp_kind, 165.942180000000_fp_kind, 166.942600000000_fp_kind,&
       167.940568000000_fp_kind, 168.941259000000_fp_kind, 169.939609000000_fp_kind, 170.940492000000_fp_kind,&
       171.939449716000_fp_kind, 172.940513000000_fp_kind, 173.940048377000_fp_kind, 174.941511424000_fp_kind,&
       175.941409797000_fp_kind, 176.943230187000_fp_kind, 177.943708322000_fp_kind, 178.945825705000_fp_kind,&
       179.946559537000_fp_kind, 180.949110834000_fp_kind, 181.950563684000_fp_kind, 182.953533203000_fp_kind,&
       183.955448507000_fp_kind, 184.958862000000_fp_kind, 185.960897000000_fp_kind, 186.964573000000_fp_kind,&
       187.966903000000_fp_kind, 188.970853000000_fp_kind, 189.973376000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  73
  real(fp_kind), parameter :: isotopic_mass_073(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 154.974248000000_fp_kind, 155.972087000000_fp_kind, 156.968227445000_fp_kind,&
       157.966593000000_fp_kind, 158.963028046000_fp_kind, 159.961541678000_fp_kind, 160.958369489000_fp_kind,&
       161.957292907000_fp_kind, 162.954337194000_fp_kind, 163.953534000000_fp_kind, 164.950780287000_fp_kind,&
       165.950512000000_fp_kind, 166.948093000000_fp_kind, 167.948047000000_fp_kind, 168.946011000000_fp_kind,&
       169.946175000000_fp_kind, 170.944476000000_fp_kind, 171.944895000000_fp_kind, 172.943750000000_fp_kind,&
       173.944454000000_fp_kind, 174.943737000000_fp_kind, 175.944857000000_fp_kind, 176.944481940000_fp_kind,&
       177.945680000000_fp_kind, 178.945939050000_fp_kind, 179.947467589000_fp_kind, 180.947998528000_fp_kind,&
       181.950154612000_fp_kind, 182.951375380000_fp_kind, 183.954009958000_fp_kind, 184.955561317000_fp_kind,&
       185.958553036000_fp_kind, 186.960391000000_fp_kind, 187.963596000000_fp_kind, 188.965690000000_fp_kind,&
       189.969168000000_fp_kind, 190.971530000000_fp_kind, 191.975201000000_fp_kind, 192.977660000000_fp_kind,&
       193.981610000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  74
  real(fp_kind), parameter :: isotopic_mass_074(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 156.978862000000_fp_kind, 157.974565000000_fp_kind, 158.972696000000_fp_kind,&
       159.968513946000_fp_kind, 160.967249000000_fp_kind, 161.963500341000_fp_kind, 162.962524251000_fp_kind,&
       163.958952445000_fp_kind, 164.958280663000_fp_kind, 165.955031952000_fp_kind, 166.954811080000_fp_kind,&
       167.951805459000_fp_kind, 168.951778689000_fp_kind, 169.949231235000_fp_kind, 170.949451000000_fp_kind,&
       171.947292000000_fp_kind, 172.947689000000_fp_kind, 173.946079000000_fp_kind, 174.946717000000_fp_kind,&
       175.945634000000_fp_kind, 176.946643000000_fp_kind, 177.945885791000_fp_kind, 178.947079378000_fp_kind,&
       179.946713304000_fp_kind, 180.948218733000_fp_kind, 181.948205636000_fp_kind, 182.950224416000_fp_kind,&
       183.950933180000_fp_kind, 184.953421206000_fp_kind, 185.954365140000_fp_kind, 186.957161249000_fp_kind,&
       187.958488325000_fp_kind, 188.961557000000_fp_kind, 189.963103542000_fp_kind, 190.966531000000_fp_kind,&
       191.968202000000_fp_kind, 192.971884000000_fp_kind, 193.973795000000_fp_kind, 194.977735000000_fp_kind,&
       195.979882000000_fp_kind, 196.984036000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  75
  real(fp_kind), parameter :: isotopic_mass_075(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 158.984106000000_fp_kind, 159.981880000000_fp_kind, 160.977624313000_fp_kind,&
       161.975896000000_fp_kind, 162.972085434000_fp_kind, 163.970507122000_fp_kind, 164.967085831000_fp_kind,&
       165.965821216000_fp_kind, 166.962604000000_fp_kind, 167.961572607000_fp_kind, 168.958765979000_fp_kind,&
       169.958234844000_fp_kind, 170.955716000000_fp_kind, 171.955376165000_fp_kind, 172.953243000000_fp_kind,&
       173.953115000000_fp_kind, 174.951381000000_fp_kind, 175.951623000000_fp_kind, 176.950328000000_fp_kind,&
       177.950989000000_fp_kind, 178.949989686000_fp_kind, 179.950791568000_fp_kind, 180.950061507000_fp_kind,&
       181.951211560000_fp_kind, 182.950821306000_fp_kind, 183.952528073000_fp_kind, 184.952958320000_fp_kind,&
       185.954989172000_fp_kind, 186.955752217000_fp_kind, 187.958113658000_fp_kind, 188.959227764000_fp_kind,&
       189.961800064000_fp_kind, 190.963123322000_fp_kind, 191.966088000000_fp_kind, 192.967545000000_fp_kind,&
       193.970735000000_fp_kind, 194.972560000000_fp_kind, 195.975996000000_fp_kind, 196.978153000000_fp_kind,&
       197.981760000000_fp_kind, 198.984187000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  76
  real(fp_kind), parameter :: isotopic_mass_076(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 160.989054000000_fp_kind, 161.984434000000_fp_kind, 162.982462000000_fp_kind,&
       163.978073158000_fp_kind, 164.976654000000_fp_kind, 165.972698135000_fp_kind, 166.971552304000_fp_kind,&
       167.967799050000_fp_kind, 168.967017521000_fp_kind, 169.963579273000_fp_kind, 170.963180402000_fp_kind,&
       171.960017309000_fp_kind, 172.959808387000_fp_kind, 173.957063192000_fp_kind, 174.956945126000_fp_kind,&
       175.954770315000_fp_kind, 176.954957902000_fp_kind, 177.953253334000_fp_kind, 178.953815985000_fp_kind,&
       179.952381665000_fp_kind, 180.953247188000_fp_kind, 181.952110154000_fp_kind, 182.953125028000_fp_kind,&
       183.952492919000_fp_kind, 184.954045969000_fp_kind, 185.953837569000_fp_kind, 186.955749569000_fp_kind,&
       187.955837292000_fp_kind, 188.958145949000_fp_kind, 189.958445442000_fp_kind, 190.960928105000_fp_kind,&
       191.961478765000_fp_kind, 192.964149637000_fp_kind, 193.965179407000_fp_kind, 194.968318000000_fp_kind,&
       195.969643261000_fp_kind, 196.973076000000_fp_kind, 197.974664000000_fp_kind, 198.978239000000_fp_kind,&
       199.980086000000_fp_kind, 200.984069000000_fp_kind, 201.986548000000_fp_kind, 202.992195000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  77
  real(fp_kind), parameter :: isotopic_mass_077(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 162.994299000000_fp_kind, 163.991966000000_fp_kind, 164.987552000000_fp_kind,&
       165.985716000000_fp_kind, 166.981671973000_fp_kind, 167.979960978000_fp_kind, 168.976281743000_fp_kind,&
       169.975113000000_fp_kind, 170.971645520000_fp_kind, 171.970607035000_fp_kind, 172.967505477000_fp_kind,&
       173.966949939000_fp_kind, 174.964149519000_fp_kind, 175.963626261000_fp_kind, 176.961301500000_fp_kind,&
       177.961079395000_fp_kind, 178.959117594000_fp_kind, 179.959229446000_fp_kind, 180.957634691000_fp_kind,&
       181.958076296000_fp_kind, 182.956841231000_fp_kind, 183.957476000000_fp_kind, 184.956698000000_fp_kind,&
       185.957946754000_fp_kind, 186.957542000000_fp_kind, 187.958834999000_fp_kind, 188.958722602000_fp_kind,&
       189.960543374000_fp_kind, 190.960591455000_fp_kind, 191.962602414000_fp_kind, 192.962923753000_fp_kind,&
       193.965075703000_fp_kind, 194.965976898000_fp_kind, 195.968399669000_fp_kind, 196.969657217000_fp_kind,&
       197.972399000000_fp_kind, 198.973807097000_fp_kind, 199.976844000000_fp_kind, 200.978701000000_fp_kind,&
       201.982136000000_fp_kind, 202.984573000000_fp_kind, 203.989726000000_fp_kind, 204.993988000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  78
  real(fp_kind), parameter :: isotopic_mass_078(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 164.999658000000_fp_kind, 165.994866000000_fp_kind, 166.992750000000_fp_kind,&
       167.988180196000_fp_kind, 168.986619000000_fp_kind, 169.982502087000_fp_kind, 170.981248868000_fp_kind,&
       171.977341059000_fp_kind, 172.976449922000_fp_kind, 173.972820431000_fp_kind, 174.972400593000_fp_kind,&
       175.968938162000_fp_kind, 176.968469541000_fp_kind, 177.965649288000_fp_kind, 178.965358742000_fp_kind,&
       179.963038010000_fp_kind, 180.963089946000_fp_kind, 181.961171605000_fp_kind, 182.961595895000_fp_kind,&
       183.959921929000_fp_kind, 184.960613659000_fp_kind, 185.959350845000_fp_kind, 186.960616646000_fp_kind,&
       187.959397521000_fp_kind, 188.960848485000_fp_kind, 189.959949823000_fp_kind, 190.961676261000_fp_kind,&
       191.961042667000_fp_kind, 192.962984546000_fp_kind, 193.962683498000_fp_kind, 194.964794325000_fp_kind,&
       195.964954648000_fp_kind, 196.967343030000_fp_kind, 197.967896718000_fp_kind, 198.970597022000_fp_kind,&
       199.971444609000_fp_kind, 200.974513305000_fp_kind, 201.975639000000_fp_kind, 202.979055000000_fp_kind,&
       203.981084000000_fp_kind, 204.986237000000_fp_kind, 205.990080000000_fp_kind, 206.995556000000_fp_kind,&
       207.999463000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  79
  real(fp_kind), parameter :: isotopic_mass_079(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 168.002716000000_fp_kind, 168.998080000000_fp_kind,&
       169.996024000000_fp_kind, 170.991881533000_fp_kind, 171.989996704000_fp_kind, 172.986224263000_fp_kind,&
       173.984908000000_fp_kind, 174.981316375000_fp_kind, 175.980116925000_fp_kind, 176.976869701000_fp_kind,&
       177.976056714000_fp_kind, 178.973173666000_fp_kind, 179.972489738000_fp_kind, 180.970079102000_fp_kind,&
       181.969614433000_fp_kind, 182.967588106000_fp_kind, 183.967451523000_fp_kind, 184.965798871000_fp_kind,&
       185.965952703000_fp_kind, 186.964542147000_fp_kind, 187.965247966000_fp_kind, 188.963948286000_fp_kind,&
       189.964751746000_fp_kind, 190.963716452000_fp_kind, 191.964817615000_fp_kind, 192.964138442000_fp_kind,&
       193.965419051000_fp_kind, 194.965037823000_fp_kind, 195.966571213000_fp_kind, 196.966570103000_fp_kind,&
       197.968243714000_fp_kind, 198.968766573000_fp_kind, 199.970756558000_fp_kind, 200.971657678000_fp_kind,&
       201.973856000000_fp_kind, 202.975154492000_fp_kind, 203.978110000000_fp_kind, 204.980064000000_fp_kind,&
       205.984766000000_fp_kind, 206.988577000000_fp_kind, 207.993655000000_fp_kind, 208.997606000000_fp_kind,&
       210.002877000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  80
  real(fp_kind), parameter :: isotopic_mass_080(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 170.005814000000_fp_kind, 171.003585000000_fp_kind,&
       171.998860581000_fp_kind, 172.997143000000_fp_kind, 173.992870575000_fp_kind, 174.991444451000_fp_kind,&
       175.987348670000_fp_kind, 176.986284590000_fp_kind, 177.982484756000_fp_kind, 178.981821759000_fp_kind,&
       179.978260180000_fp_kind, 180.977819368000_fp_kind, 181.974689173000_fp_kind, 182.974444652000_fp_kind,&
       183.971717709000_fp_kind, 184.971890696000_fp_kind, 185.969362061000_fp_kind, 186.969813540000_fp_kind,&
       187.967580738000_fp_kind, 188.968194776000_fp_kind, 189.966322250000_fp_kind, 190.967158301000_fp_kind,&
       191.965634263000_fp_kind, 192.966653395000_fp_kind, 193.965449108000_fp_kind, 194.966705809000_fp_kind,&
       195.965833445000_fp_kind, 196.967213715000_fp_kind, 197.966769177000_fp_kind, 198.968280994000_fp_kind,&
       199.968326941000_fp_kind, 200.970303054000_fp_kind, 201.970643604000_fp_kind, 202.972872396000_fp_kind,&
       203.973494037000_fp_kind, 204.976073151000_fp_kind, 205.977513837000_fp_kind, 206.982300000000_fp_kind,&
       207.985759000000_fp_kind, 208.990757000000_fp_kind, 209.994310000000_fp_kind, 210.999581000000_fp_kind,&
       212.003242000000_fp_kind, 213.008803000000_fp_kind, 214.012636000000_fp_kind, 215.018368000000_fp_kind,&
       216.022459000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  81
  real(fp_kind), parameter :: isotopic_mass_081(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 176.000627731000_fp_kind, 176.996414252000_fp_kind,&
       177.995047000000_fp_kind, 178.991122185000_fp_kind, 179.989918950000_fp_kind, 180.986259978000_fp_kind,&
       181.985692649000_fp_kind, 182.982192843000_fp_kind, 183.981874973000_fp_kind, 184.978789189000_fp_kind,&
       185.978654787000_fp_kind, 186.975904740000_fp_kind, 187.976020886000_fp_kind, 188.973573525000_fp_kind,&
       189.973841771000_fp_kind, 190.971784093000_fp_kind, 191.972225000000_fp_kind, 192.970501994000_fp_kind,&
       193.971081408000_fp_kind, 194.969774052000_fp_kind, 195.970481189000_fp_kind, 196.969560492000_fp_kind,&
       197.970446669000_fp_kind, 198.969877000000_fp_kind, 199.970963608000_fp_kind, 200.970820235000_fp_kind,&
       201.972108874000_fp_kind, 202.972344098000_fp_kind, 203.973863420000_fp_kind, 204.974427318000_fp_kind,&
       205.976110108000_fp_kind, 206.977418605000_fp_kind, 207.982018006000_fp_kind, 208.985351713000_fp_kind,&
       209.990072942000_fp_kind, 210.993475000000_fp_kind, 211.998335000000_fp_kind, 213.001915000000_fp_kind,&
       214.006940000000_fp_kind, 215.010768000000_fp_kind, 216.015964000000_fp_kind, 217.020032000000_fp_kind,&
       218.025454000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  82
  real(fp_kind), parameter :: isotopic_mass_082(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 178.003836171000_fp_kind, 179.002202492000_fp_kind,&
       179.997916177000_fp_kind, 180.996660600000_fp_kind, 181.992673537000_fp_kind, 182.991862527000_fp_kind,&
       183.988135634000_fp_kind, 184.987610000000_fp_kind, 185.984239409000_fp_kind, 186.983910842000_fp_kind,&
       187.980879079000_fp_kind, 188.980843658000_fp_kind, 189.978081872000_fp_kind, 190.978216455000_fp_kind,&
       191.975789598000_fp_kind, 192.976135914000_fp_kind, 193.974011788000_fp_kind, 194.974516167000_fp_kind,&
       195.972787552000_fp_kind, 196.973434737000_fp_kind, 197.972015450000_fp_kind, 198.972912620000_fp_kind,&
       199.971818546000_fp_kind, 200.972870431000_fp_kind, 201.972151613000_fp_kind, 202.973390617000_fp_kind,&
       203.973043506000_fp_kind, 204.974481682000_fp_kind, 205.974465210000_fp_kind, 206.975896821000_fp_kind,&
       207.976652005000_fp_kind, 208.981089978000_fp_kind, 209.984188381000_fp_kind, 210.988735288000_fp_kind,&
       211.991895891000_fp_kind, 212.996560796000_fp_kind, 213.999803521000_fp_kind, 215.004661591000_fp_kind,&
       216.008062000000_fp_kind, 217.013162000000_fp_kind, 218.016779000000_fp_kind, 219.022136000000_fp_kind,&
       220.025905000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  83
  real(fp_kind), parameter :: isotopic_mass_083(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 184.001347000000_fp_kind, 184.997600000000_fp_kind,&
       185.996623169000_fp_kind, 186.993147272000_fp_kind, 187.992276064000_fp_kind, 188.989195139000_fp_kind,&
       189.988624828000_fp_kind, 190.985786972000_fp_kind, 191.985470077000_fp_kind, 192.982947220000_fp_kind,&
       193.982798581000_fp_kind, 194.980648759000_fp_kind, 195.980666509000_fp_kind, 196.978864927000_fp_kind,&
       197.979201316000_fp_kind, 198.977672841000_fp_kind, 199.978131290000_fp_kind, 200.976995017000_fp_kind,&
       201.977723042000_fp_kind, 202.976892077000_fp_kind, 203.977835687000_fp_kind, 204.977385182000_fp_kind,&
       205.978498843000_fp_kind, 206.978470551000_fp_kind, 207.979742060000_fp_kind, 208.980398599000_fp_kind,&
       209.984120237000_fp_kind, 210.987268715000_fp_kind, 211.991285030000_fp_kind, 212.994383570000_fp_kind,&
       213.998710909000_fp_kind, 215.001749095000_fp_kind, 216.006305985000_fp_kind, 217.009372000000_fp_kind,&
       218.014188000000_fp_kind, 219.017520000000_fp_kind, 220.022501000000_fp_kind, 221.025980000000_fp_kind,&
       222.031079000000_fp_kind, 223.034611000000_fp_kind, 224.039796000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  84
  real(fp_kind), parameter :: isotopic_mass_084(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 186.004403174000_fp_kind, 187.003031482000_fp_kind,&
       187.999415586000_fp_kind, 188.998473425000_fp_kind, 189.995101731000_fp_kind, 190.994558494000_fp_kind,&
       191.991340274000_fp_kind, 192.991062421000_fp_kind, 193.988186058000_fp_kind, 194.988065781000_fp_kind,&
       195.985540722000_fp_kind, 196.985621939000_fp_kind, 197.983388753000_fp_kind, 198.983640445000_fp_kind,&
       199.981812355000_fp_kind, 200.982263799000_fp_kind, 201.980738934000_fp_kind, 202.981416072000_fp_kind,&
       203.980310078000_fp_kind, 204.981190006000_fp_kind, 205.980473662000_fp_kind, 206.981593334000_fp_kind,&
       207.981246035000_fp_kind, 208.982430361000_fp_kind, 209.982873686000_fp_kind, 210.986653171000_fp_kind,&
       211.988867982000_fp_kind, 212.992857154000_fp_kind, 213.995201287000_fp_kind, 214.999418385000_fp_kind,&
       216.001913416000_fp_kind, 217.006316145000_fp_kind, 218.008971234000_fp_kind, 219.013614000000_fp_kind,&
       220.016386000000_fp_kind, 221.021228000000_fp_kind, 222.024140000000_fp_kind, 223.029070000000_fp_kind,&
       224.032110000000_fp_kind, 225.037123000000_fp_kind, 226.040310000000_fp_kind, 227.045390000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  85
  real(fp_kind), parameter :: isotopic_mass_085(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 191.004148081000_fp_kind, 192.003140912000_fp_kind, 192.999927725000_fp_kind,&
       193.999230816000_fp_kind, 194.996274480000_fp_kind, 195.995799034000_fp_kind, 196.993177353000_fp_kind,&
       197.992797864000_fp_kind, 198.990527715000_fp_kind, 199.990351099000_fp_kind, 200.988417058000_fp_kind,&
       201.988625686000_fp_kind, 202.986942904000_fp_kind, 203.987251393000_fp_kind, 204.986060546000_fp_kind,&
       205.986645768000_fp_kind, 206.985799715000_fp_kind, 207.986613011000_fp_kind, 208.986168701000_fp_kind,&
       209.987147423000_fp_kind, 210.987496226000_fp_kind, 211.990737301000_fp_kind, 212.992936593000_fp_kind,&
       213.996372331000_fp_kind, 214.998651002000_fp_kind, 216.002422643000_fp_kind, 217.004717794000_fp_kind,&
       218.008695941000_fp_kind, 219.011160587000_fp_kind, 220.015433000000_fp_kind, 221.018017000000_fp_kind,&
       222.022494000000_fp_kind, 223.025151000000_fp_kind, 224.029749000000_fp_kind, 225.032528000000_fp_kind,&
       226.037209000000_fp_kind, 227.040183000000_fp_kind, 228.044960000000_fp_kind, 229.048191000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  86
  real(fp_kind), parameter :: isotopic_mass_086(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 193.009707973000_fp_kind, 194.006145636000_fp_kind, 195.005421703000_fp_kind,&
       196.002120431000_fp_kind, 197.001621446000_fp_kind, 197.998679197000_fp_kind, 198.998325436000_fp_kind,&
       199.995705335000_fp_kind, 200.995590511000_fp_kind, 201.993263982000_fp_kind, 202.993361155000_fp_kind,&
       203.991443729000_fp_kind, 204.991723228000_fp_kind, 205.990195409000_fp_kind, 206.990730224000_fp_kind,&
       207.989634513000_fp_kind, 208.990401389000_fp_kind, 209.989688862000_fp_kind, 210.990600767000_fp_kind,&
       211.990703946000_fp_kind, 212.993885147000_fp_kind, 213.995362650000_fp_kind, 214.998745037000_fp_kind,&
       216.000271942000_fp_kind, 217.003927632000_fp_kind, 218.005601123000_fp_kind, 219.009478683000_fp_kind,&
       220.011392443000_fp_kind, 221.015535637000_fp_kind, 222.017576017000_fp_kind, 223.021889283000_fp_kind,&
       224.024095803000_fp_kind, 225.028485572000_fp_kind, 226.030861380000_fp_kind, 227.035304393000_fp_kind,&
       228.037835415000_fp_kind, 229.042257272000_fp_kind, 230.045271000000_fp_kind, 231.049973000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  87
  real(fp_kind), parameter :: isotopic_mass_087(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 197.011008086000_fp_kind,&
       198.010282081000_fp_kind, 199.007269384000_fp_kind, 200.006584666000_fp_kind, 201.003852491000_fp_kind,&
       202.003329637000_fp_kind, 203.000940867000_fp_kind, 204.000651972000_fp_kind, 204.998593854000_fp_kind,&
       205.998661441000_fp_kind, 206.996941450000_fp_kind, 207.997139082000_fp_kind, 208.995939701000_fp_kind,&
       209.996410596000_fp_kind, 210.995555189000_fp_kind, 211.996225420000_fp_kind, 212.996184410000_fp_kind,&
       213.998971193000_fp_kind, 215.000341534000_fp_kind, 216.003189523000_fp_kind, 217.004631980000_fp_kind,&
       218.007578620000_fp_kind, 219.009250664000_fp_kind, 220.012326789000_fp_kind, 221.014253714000_fp_kind,&
       222.017582615000_fp_kind, 223.019734241000_fp_kind, 224.023348096000_fp_kind, 225.025572466000_fp_kind,&
       226.029544512000_fp_kind, 227.031865413000_fp_kind, 228.035839433000_fp_kind, 229.038291443000_fp_kind,&
       230.042390787000_fp_kind, 231.045175353000_fp_kind, 232.049461219000_fp_kind, 233.052517833000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  88
  real(fp_kind), parameter :: isotopic_mass_088(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 201.012814699000_fp_kind, 202.009742305000_fp_kind, 203.009233907000_fp_kind,&
       204.006506855000_fp_kind, 205.006230692000_fp_kind, 206.003827842000_fp_kind, 207.003772420000_fp_kind,&
       208.001855012000_fp_kind, 209.001994902000_fp_kind, 210.000475406000_fp_kind, 211.000893049000_fp_kind,&
       211.999786619000_fp_kind, 213.000370971000_fp_kind, 214.000099560000_fp_kind, 215.002718208000_fp_kind,&
       216.003533534000_fp_kind, 217.006322676000_fp_kind, 218.007134297000_fp_kind, 219.010084715000_fp_kind,&
       220.011027542000_fp_kind, 221.013917293000_fp_kind, 222.015373371000_fp_kind, 223.018500648000_fp_kind,&
       224.020210361000_fp_kind, 225.023610502000_fp_kind, 226.025408186000_fp_kind, 227.029176205000_fp_kind,&
       228.031068574000_fp_kind, 229.034956703000_fp_kind, 230.037054776000_fp_kind, 231.041027085000_fp_kind,&
       232.043475267000_fp_kind, 233.047594570000_fp_kind, 234.050382100000_fp_kind, 235.054890000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  89
  real(fp_kind), parameter :: isotopic_mass_089(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 205.015144152000_fp_kind,&
       206.014476477000_fp_kind, 207.011965967000_fp_kind, 208.011552251000_fp_kind, 209.009495375000_fp_kind,&
       210.009408625000_fp_kind, 211.007668846000_fp_kind, 212.007836442000_fp_kind, 213.006592665000_fp_kind,&
       214.006906400000_fp_kind, 215.006474061000_fp_kind, 216.008749101000_fp_kind, 217.009342325000_fp_kind,&
       218.011648860000_fp_kind, 219.012420425000_fp_kind, 220.014754527000_fp_kind, 221.015599721000_fp_kind,&
       222.017844232000_fp_kind, 223.019135982000_fp_kind, 224.021722249000_fp_kind, 225.023228601000_fp_kind,&
       226.026096999000_fp_kind, 227.027750594000_fp_kind, 228.031019685000_fp_kind, 229.032947000000_fp_kind,&
       230.036327000000_fp_kind, 231.038393000000_fp_kind, 232.042034000000_fp_kind, 233.044346000000_fp_kind,&
       234.048139000000_fp_kind, 235.050840000000_fp_kind, 236.054988000000_fp_kind, 237.057993000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  90
  real(fp_kind), parameter :: isotopic_mass_090(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       208.017915348000_fp_kind, 209.017601000000_fp_kind, 210.015093515000_fp_kind, 211.014896923000_fp_kind,&
       212.013001570000_fp_kind, 213.013011470000_fp_kind, 214.011481480000_fp_kind, 215.011724640000_fp_kind,&
       216.011055933000_fp_kind, 217.013103443000_fp_kind, 218.013276248000_fp_kind, 219.015526432000_fp_kind,&
       220.015769866000_fp_kind, 221.018185757000_fp_kind, 222.018468220000_fp_kind, 223.020811083000_fp_kind,&
       224.021466137000_fp_kind, 225.023950975000_fp_kind, 226.024903699000_fp_kind, 227.027702546000_fp_kind,&
       228.028739741000_fp_kind, 229.031761357000_fp_kind, 230.033132267000_fp_kind, 231.036302764000_fp_kind,&
       232.038053606000_fp_kind, 233.041580126000_fp_kind, 234.043599801000_fp_kind, 235.047255000000_fp_kind,&
       236.049657000000_fp_kind, 237.053629000000_fp_kind, 238.056388000000_fp_kind, 239.060655000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  91
  real(fp_kind), parameter :: isotopic_mass_091(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 211.023674036000_fp_kind, 212.023184819000_fp_kind, 213.021099644000_fp_kind,&
       214.020891055000_fp_kind, 215.019113955000_fp_kind, 216.019134633000_fp_kind, 217.018309024000_fp_kind,&
       218.020021133000_fp_kind, 219.019949909000_fp_kind, 220.021769753000_fp_kind, 221.021873393000_fp_kind,&
       222.023687064000_fp_kind, 223.023980414000_fp_kind, 224.025617286000_fp_kind, 225.026147927000_fp_kind,&
       226.027948217000_fp_kind, 227.028803586000_fp_kind, 228.031050758000_fp_kind, 229.032095585000_fp_kind,&
       230.034539717000_fp_kind, 231.035882500000_fp_kind, 232.038590205000_fp_kind, 233.040246535000_fp_kind,&
       234.043305555000_fp_kind, 235.045399000000_fp_kind, 236.048668000000_fp_kind, 237.051023000000_fp_kind,&
       238.054637000000_fp_kind, 239.057260000000_fp_kind, 240.061203000000_fp_kind, 241.064134000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  92
  real(fp_kind), parameter :: isotopic_mass_092(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 215.026719774000_fp_kind,&
       216.024762829000_fp_kind, 217.024660000000_fp_kind, 218.023504877000_fp_kind, 219.025009233000_fp_kind,&
       220.024706000000_fp_kind, 221.026323297000_fp_kind, 222.026057957000_fp_kind, 223.027960754000_fp_kind,&
       224.027635913000_fp_kind, 225.029385050000_fp_kind, 226.029338669000_fp_kind, 227.031181124000_fp_kind,&
       228.031368959000_fp_kind, 229.033505976000_fp_kind, 230.033940114000_fp_kind, 231.036292180000_fp_kind,&
       232.037154765000_fp_kind, 233.039634294000_fp_kind, 234.040950296000_fp_kind, 235.043928117000_fp_kind,&
       236.045566130000_fp_kind, 237.048728309000_fp_kind, 238.050786936000_fp_kind, 239.054291989000_fp_kind,&
       240.056592411000_fp_kind, 241.060330000000_fp_kind, 242.062931000000_fp_kind, 243.067075000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  93
  real(fp_kind), parameter :: isotopic_mass_093(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 219.031601865000_fp_kind, 220.032716280000_fp_kind, 221.032110000000_fp_kind,&
       222.033574706000_fp_kind, 223.032913340000_fp_kind, 224.034388030000_fp_kind, 225.033943422000_fp_kind,&
       226.035230364000_fp_kind, 227.034975012000_fp_kind, 228.036313000000_fp_kind, 229.036287269000_fp_kind,&
       230.037828060000_fp_kind, 231.038243598000_fp_kind, 232.040107000000_fp_kind, 233.040739421000_fp_kind,&
       234.042893245000_fp_kind, 235.044061518000_fp_kind, 236.046568296000_fp_kind, 237.048171640000_fp_kind,&
       238.050944603000_fp_kind, 239.052937538000_fp_kind, 240.056163778000_fp_kind, 241.058309671000_fp_kind,&
       242.061639548000_fp_kind, 243.064204000000_fp_kind, 244.067891000000_fp_kind, 245.070693000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  94
  real(fp_kind), parameter :: isotopic_mass_094(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 221.038572000000_fp_kind, 222.037638000000_fp_kind, 223.038777000000_fp_kind,&
       224.037875000000_fp_kind, 225.038970000000_fp_kind, 226.038250000000_fp_kind, 227.039474000000_fp_kind,&
       228.038763325000_fp_kind, 229.040145099000_fp_kind, 230.039648313000_fp_kind, 231.041125946000_fp_kind,&
       232.041182133000_fp_kind, 233.042997411000_fp_kind, 234.043317489000_fp_kind, 235.045284609000_fp_kind,&
       236.046056661000_fp_kind, 237.048407888000_fp_kind, 238.049558175000_fp_kind, 239.052161596000_fp_kind,&
       240.053811740000_fp_kind, 241.056849651000_fp_kind, 242.058740979000_fp_kind, 243.062002068000_fp_kind,&
       244.064204401000_fp_kind, 245.067824554000_fp_kind, 246.070204172000_fp_kind, 247.074300000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  95
  real(fp_kind), parameter :: isotopic_mass_095(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 223.045840000000_fp_kind, 224.046442000000_fp_kind, 225.045508000000_fp_kind,&
       226.046130000000_fp_kind, 227.045282000000_fp_kind, 228.046001000000_fp_kind, 229.045282534000_fp_kind,&
       230.046025000000_fp_kind, 231.045529000000_fp_kind, 232.046613000000_fp_kind, 233.046468000000_fp_kind,&
       234.047731000000_fp_kind, 235.047906478000_fp_kind, 236.049427000000_fp_kind, 237.049995000000_fp_kind,&
       238.051982531000_fp_kind, 239.053022729000_fp_kind, 240.055298374000_fp_kind, 241.056827343000_fp_kind,&
       242.059547358000_fp_kind, 243.061379889000_fp_kind, 244.064282892000_fp_kind, 245.066452827000_fp_kind,&
       246.069774000000_fp_kind, 247.072092000000_fp_kind, 248.075752000000_fp_kind, 249.078480000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  96
  real(fp_kind), parameter :: isotopic_mass_096(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 231.050746000000_fp_kind,&
       232.049740000000_fp_kind, 233.050771485000_fp_kind, 234.050158568000_fp_kind, 235.051545000000_fp_kind,&
       236.051372112000_fp_kind, 237.052868988000_fp_kind, 238.053081606000_fp_kind, 239.054908519000_fp_kind,&
       240.055528233000_fp_kind, 241.057651218000_fp_kind, 242.058834187000_fp_kind, 243.061387329000_fp_kind,&
       244.062750622000_fp_kind, 245.065491047000_fp_kind, 246.067222016000_fp_kind, 247.070352678000_fp_kind,&
       248.072349086000_fp_kind, 249.075953992000_fp_kind, 250.078357541000_fp_kind, 251.082284988000_fp_kind,&
       252.084870000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  97
  real(fp_kind), parameter :: isotopic_mass_097(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 233.056652000000_fp_kind,&
       234.057322000000_fp_kind, 235.056651000000_fp_kind, 236.057479000000_fp_kind, 237.057123000000_fp_kind,&
       238.058204000000_fp_kind, 239.058239000000_fp_kind, 240.059758000000_fp_kind, 241.060098000000_fp_kind,&
       242.061999000000_fp_kind, 243.063005905000_fp_kind, 244.065178969000_fp_kind, 245.066359814000_fp_kind,&
       246.068671300000_fp_kind, 247.070305889000_fp_kind, 248.073141689000_fp_kind, 249.074983118000_fp_kind,&
       250.078317195000_fp_kind, 251.080760555000_fp_kind, 252.084310000000_fp_kind, 253.086880000000_fp_kind,&
       254.090600000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  98
  real(fp_kind), parameter :: isotopic_mass_098(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 237.062199272000_fp_kind, 238.061490000000_fp_kind, 239.062482000000_fp_kind,&
       240.062253447000_fp_kind, 241.063690000000_fp_kind, 242.063754544000_fp_kind, 243.065475000000_fp_kind,&
       244.065999447000_fp_kind, 245.068046755000_fp_kind, 246.068803685000_fp_kind, 247.070971348000_fp_kind,&
       248.072182905000_fp_kind, 249.074850428000_fp_kind, 250.076404494000_fp_kind, 251.079587171000_fp_kind,&
       252.081626507000_fp_kind, 253.085133723000_fp_kind, 254.087323575000_fp_kind, 255.091046000000_fp_kind,&
       256.093442000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z =  99
  real(fp_kind), parameter :: isotopic_mass_099(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 239.068310000000_fp_kind, 240.068949000000_fp_kind, 241.068592000000_fp_kind,&
       242.069567000000_fp_kind, 243.069508000000_fp_kind, 244.070881000000_fp_kind, 245.071192000000_fp_kind,&
       246.072806474000_fp_kind, 247.073621929000_fp_kind, 248.075469000000_fp_kind, 249.076409000000_fp_kind,&
       250.078611000000_fp_kind, 251.079991431000_fp_kind, 252.082979173000_fp_kind, 253.084821241000_fp_kind,&
       254.088024337000_fp_kind, 255.090273504000_fp_kind, 256.093597000000_fp_kind, 257.095979000000_fp_kind,&
       258.099520000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 100
  real(fp_kind), parameter :: isotopic_mass_100(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 241.074311000000_fp_kind, 242.073430000000_fp_kind, 243.074414000000_fp_kind,&
       244.074036000000_fp_kind, 245.075354000000_fp_kind, 246.075353334000_fp_kind, 247.076944000000_fp_kind,&
       248.077185451000_fp_kind, 249.078926042000_fp_kind, 250.079519765000_fp_kind, 251.081545130000_fp_kind,&
       252.082466019000_fp_kind, 253.085180945000_fp_kind, 254.086852424000_fp_kind, 255.089963495000_fp_kind,&
       256.091771699000_fp_kind, 257.095105419000_fp_kind, 258.097077000000_fp_kind, 259.100596000000_fp_kind,&
       260.102809000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 101
  real(fp_kind), parameter :: isotopic_mass_101(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 244.081157000000_fp_kind, 245.080864000000_fp_kind,&
       246.081713000000_fp_kind, 247.081520000000_fp_kind, 248.082607000000_fp_kind, 249.082857155000_fp_kind,&
       250.084164934000_fp_kind, 251.084774287000_fp_kind, 252.086385000000_fp_kind, 253.087143000000_fp_kind,&
       254.089590000000_fp_kind, 255.091081702000_fp_kind, 256.093888000000_fp_kind, 257.095537343000_fp_kind,&
       258.098433634000_fp_kind, 259.100445000000_fp_kind, 260.103650000000_fp_kind, 261.105828000000_fp_kind,&
       262.109144000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 102
  real(fp_kind), parameter :: isotopic_mass_102(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       248.086623000000_fp_kind, 249.087802000000_fp_kind, 250.087565000000_fp_kind, 251.088942000000_fp_kind,&
       252.088966070000_fp_kind, 253.090562780000_fp_kind, 254.090954211000_fp_kind, 255.093196439000_fp_kind,&
       256.094281912000_fp_kind, 257.096884203000_fp_kind, 258.098205000000_fp_kind, 259.100998364000_fp_kind,&
       260.102641000000_fp_kind, 261.105696000000_fp_kind, 262.107463000000_fp_kind, 263.110714000000_fp_kind,&
       264.112734000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 103
  real(fp_kind), parameter :: isotopic_mass_103(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 251.094289000000_fp_kind, 252.095048000000_fp_kind, 253.095033850000_fp_kind,&
       254.096238813000_fp_kind, 255.096562399000_fp_kind, 256.098494024000_fp_kind, 257.099480000000_fp_kind,&
       258.101753000000_fp_kind, 259.102900000000_fp_kind, 260.105504000000_fp_kind, 261.106879000000_fp_kind,&
       262.109615000000_fp_kind, 263.111293000000_fp_kind, 264.114198000000_fp_kind, 265.116193000000_fp_kind,&
       266.119874000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 104
  real(fp_kind), parameter :: isotopic_mass_104(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 253.100528000000_fp_kind, 254.100055000000_fp_kind, 255.101267000000_fp_kind,&
       256.101151464000_fp_kind, 257.102916796000_fp_kind, 258.103429895000_fp_kind, 259.105601000000_fp_kind,&
       260.106440000000_fp_kind, 261.108769591000_fp_kind, 262.109923000000_fp_kind, 263.112461000000_fp_kind,&
       264.113876000000_fp_kind, 265.116683000000_fp_kind, 266.118236000000_fp_kind, 267.121787000000_fp_kind,&
       268.123968000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 105
  real(fp_kind), parameter :: isotopic_mass_105(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 255.106919000000_fp_kind, 256.107674000000_fp_kind, 257.107520042000_fp_kind,&
       258.108972995000_fp_kind, 259.109491859000_fp_kind, 260.111297000000_fp_kind, 261.111979000000_fp_kind,&
       262.114067000000_fp_kind, 263.114987000000_fp_kind, 264.117297000000_fp_kind, 265.118500000000_fp_kind,&
       266.121032000000_fp_kind, 267.122399000000_fp_kind, 268.125669000000_fp_kind, 269.127911000000_fp_kind,&
       270.131399000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 106
  real(fp_kind), parameter :: isotopic_mass_106(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 258.113040000000_fp_kind, 259.114353000000_fp_kind,&
       260.114383435000_fp_kind, 261.115948135000_fp_kind, 262.116338978000_fp_kind, 263.118299000000_fp_kind,&
       264.118930000000_fp_kind, 265.121089000000_fp_kind, 266.121973000000_fp_kind, 267.124323000000_fp_kind,&
       268.125389000000_fp_kind, 269.128495000000_fp_kind, 270.130362000000_fp_kind, 271.133782000000_fp_kind,&
       272.135825000000_fp_kind, 273.139475000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 107
  real(fp_kind), parameter :: isotopic_mass_107(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 260.121443000000_fp_kind, 261.121395733000_fp_kind,&
       262.122654688000_fp_kind, 263.122916000000_fp_kind, 264.124486000000_fp_kind, 265.124955000000_fp_kind,&
       266.126790000000_fp_kind, 267.127499000000_fp_kind, 268.129584000000_fp_kind, 269.130411000000_fp_kind,&
       270.133366000000_fp_kind, 271.135115000000_fp_kind, 272.138259000000_fp_kind, 273.140294000000_fp_kind,&
       274.143599000000_fp_kind, 275.145766000000_fp_kind, 276.149169000000_fp_kind, 277.151477000000_fp_kind,&
       278.154988000000_fp_kind]
  ! Z = 108
  real(fp_kind), parameter :: isotopic_mass_108(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 263.128479000000_fp_kind,&
       264.128356330000_fp_kind, 265.129791744000_fp_kind, 266.130048783000_fp_kind, 267.131678000000_fp_kind,&
       268.132011000000_fp_kind, 269.133649000000_fp_kind, 270.134313000000_fp_kind, 271.137082000000_fp_kind,&
       272.138492000000_fp_kind, 273.141458000000_fp_kind, 274.143217000000_fp_kind, 275.146530000000_fp_kind,&
       276.148348000000_fp_kind, 277.151772000000_fp_kind, 278.153753000000_fp_kind, 279.157274000000_fp_kind,&
       280.159335000000_fp_kind]
  ! Z = 109
  real(fp_kind), parameter :: isotopic_mass_109(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 265.135937000000_fp_kind,&
       266.137062253000_fp_kind, 267.137189000000_fp_kind, 268.138649000000_fp_kind, 269.138809000000_fp_kind,&
       270.140322000000_fp_kind, 271.140741000000_fp_kind, 272.143298000000_fp_kind, 273.144695000000_fp_kind,&
       274.147343000000_fp_kind, 275.148972000000_fp_kind, 276.151705000000_fp_kind, 277.153525000000_fp_kind,&
       278.156487000000_fp_kind, 279.158439000000_fp_kind, 280.161579000000_fp_kind, 281.163608000000_fp_kind,&
       282.166888000000_fp_kind]
  ! Z = 110
  real(fp_kind), parameter :: isotopic_mass_110(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind, 267.143726000000_fp_kind,&
       268.143477000000_fp_kind, 269.144750965000_fp_kind, 270.144586620000_fp_kind, 271.145951000000_fp_kind,&
       272.146091000000_fp_kind, 273.148455000000_fp_kind, 274.149434000000_fp_kind, 275.152085000000_fp_kind,&
       276.153022000000_fp_kind, 277.155763000000_fp_kind, 278.157007000000_fp_kind, 279.159984000000_fp_kind,&
       280.161375000000_fp_kind, 281.164545000000_fp_kind, 282.166174000000_fp_kind, 283.169437000000_fp_kind,&
       284.171187000000_fp_kind]
  ! Z = 111
  real(fp_kind), parameter :: isotopic_mass_111(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind, 272.153273000000_fp_kind, 273.153393000000_fp_kind,&
       274.155247000000_fp_kind, 275.156088000000_fp_kind, 276.158226000000_fp_kind, 277.159322000000_fp_kind,&
       278.161590000000_fp_kind, 279.162880000000_fp_kind, 280.165204000000_fp_kind, 281.166757000000_fp_kind,&
       282.169343000000_fp_kind, 283.171101000000_fp_kind, 284.173882000000_fp_kind, 285.175771000000_fp_kind,&
       286.178756000000_fp_kind]
  ! Z = 112
  real(fp_kind), parameter :: isotopic_mass_112(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       276.161418000000_fp_kind, 277.163535000000_fp_kind, 278.164083000000_fp_kind, 279.166422000000_fp_kind,&
       280.167102000000_fp_kind, 281.169563000000_fp_kind, 282.170507000000_fp_kind, 283.173202000000_fp_kind,&
       284.174360000000_fp_kind, 285.177227000000_fp_kind, 286.178691000000_fp_kind, 287.181826000000_fp_kind,&
       288.183501000000_fp_kind]
  ! Z = 113
  real(fp_kind), parameter :: isotopic_mass_113(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       278.170725000000_fp_kind, 279.171187000000_fp_kind, 280.173098000000_fp_kind, 281.173710000000_fp_kind,&
       282.175770000000_fp_kind, 283.176666000000_fp_kind, 284.178843000000_fp_kind, 285.180106000000_fp_kind,&
       286.182456000000_fp_kind, 287.184064000000_fp_kind, 288.186764000000_fp_kind, 289.188461000000_fp_kind,&
       290.191429000000_fp_kind]
  ! Z = 114
  real(fp_kind), parameter :: isotopic_mass_114(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
       284.181192000000_fp_kind, 285.183503000000_fp_kind, 286.184226000000_fp_kind, 287.186720000000_fp_kind,&
       288.187781000000_fp_kind, 289.190517000000_fp_kind, 290.191875000000_fp_kind, 291.194848000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 115
  real(fp_kind), parameter :: isotopic_mass_115(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 287.190820000000_fp_kind, 288.192879000000_fp_kind, 289.193971000000_fp_kind,&
       290.196235000000_fp_kind, 291.197725000000_fp_kind, 292.200323000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 116
  real(fp_kind), parameter :: isotopic_mass_116(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 289.198023000000_fp_kind, 290.198635000000_fp_kind, 291.201014000000_fp_kind,&
       292.201969000000_fp_kind, 293.204583000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 117
  real(fp_kind), parameter :: isotopic_mass_117(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 291.205748000000_fp_kind, 292.207861000000_fp_kind, 293.208727000000_fp_kind,&
       294.210840000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! Z = 118
  real(fp_kind), parameter :: isotopic_mass_118(min_nmz:max_nmz) = [&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind, 293.213423000000_fp_kind, 294.213979000000_fp_kind, 295.216178000000_fp_kind,&
         0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,   0.000000000000_fp_kind,&
         0.000000000000_fp_kind]
  ! END OF AME2020 DATA (inserted from the file isotopic_mass_2d.data documented above)
  ! Form final result from AME2020 data for each Z value.
  real(fp_kind), parameter :: isotopic_mass_2d(min_nmz:max_nmz, min_z:max_z) = reshape([&
       isotopic_mass_000, isotopic_mass_001, isotopic_mass_002, isotopic_mass_003, isotopic_mass_004,&
       isotopic_mass_005, isotopic_mass_006, isotopic_mass_007, isotopic_mass_008, isotopic_mass_009,&
       isotopic_mass_010, isotopic_mass_011, isotopic_mass_012, isotopic_mass_013, isotopic_mass_014,&
       isotopic_mass_015, isotopic_mass_016, isotopic_mass_017, isotopic_mass_018, isotopic_mass_019,&
       isotopic_mass_020, isotopic_mass_021, isotopic_mass_022, isotopic_mass_023, isotopic_mass_024,&
       isotopic_mass_025, isotopic_mass_026, isotopic_mass_027, isotopic_mass_028, isotopic_mass_029,&
       isotopic_mass_030, isotopic_mass_031, isotopic_mass_032, isotopic_mass_033, isotopic_mass_034,&
       isotopic_mass_035, isotopic_mass_036, isotopic_mass_037, isotopic_mass_038, isotopic_mass_039,&
       isotopic_mass_040, isotopic_mass_041, isotopic_mass_042, isotopic_mass_043, isotopic_mass_044,&
       isotopic_mass_045, isotopic_mass_046, isotopic_mass_047, isotopic_mass_048, isotopic_mass_049,&
       isotopic_mass_050, isotopic_mass_051, isotopic_mass_052, isotopic_mass_053, isotopic_mass_054,&
       isotopic_mass_055, isotopic_mass_056, isotopic_mass_057, isotopic_mass_058, isotopic_mass_059,&
       isotopic_mass_060, isotopic_mass_061, isotopic_mass_062, isotopic_mass_063, isotopic_mass_064,&
       isotopic_mass_065, isotopic_mass_066, isotopic_mass_067, isotopic_mass_068, isotopic_mass_069,&
       isotopic_mass_070, isotopic_mass_071, isotopic_mass_072, isotopic_mass_073, isotopic_mass_074,&
       isotopic_mass_075, isotopic_mass_076, isotopic_mass_077, isotopic_mass_078, isotopic_mass_079,&
       isotopic_mass_080, isotopic_mass_081, isotopic_mass_082, isotopic_mass_083, isotopic_mass_084,&
       isotopic_mass_085, isotopic_mass_086, isotopic_mass_087, isotopic_mass_088, isotopic_mass_089,&
       isotopic_mass_090, isotopic_mass_091, isotopic_mass_092, isotopic_mass_093, isotopic_mass_094,&
       isotopic_mass_095, isotopic_mass_096, isotopic_mass_097, isotopic_mass_098, isotopic_mass_099,&
       isotopic_mass_100, isotopic_mass_101, isotopic_mass_102, isotopic_mass_103, isotopic_mass_104,&
       isotopic_mass_105, isotopic_mass_106, isotopic_mass_107, isotopic_mass_108, isotopic_mass_109,&
       isotopic_mass_110, isotopic_mass_111, isotopic_mass_112, isotopic_mass_113, isotopic_mass_114,&
       isotopic_mass_115, isotopic_mass_116, isotopic_mass_117, isotopic_mass_118], shape(isotopic_mass_2d))

  ! The most precise way to treat isotopes in EOS calculations is to
  ! have a separate abundance constraint equations for each isotope,
  ! but to very good approximation the energy levels that determine
  ! the partition functions (and therefore isotopic equilibrium
  ! constants) are independent of isotopic species for a given
  ! element.

  ! I adopt that approximation for FreeEOS.  In which case, all the
  ! isotopic abundance constraint equations can be analytically added
  ! (taking proper account of different isotopic nuclear-spin
  ! statistics for homonuclear molecules) to form an exact summed
  ! abundance constraint equation where the various mean equilbrium
  ! constants used are identical (as a result of the stated
  ! approximation) to the equivalent equilibrium constant for any
  ! single choice of isotope per element.  And to improve that
  ! approximation (which would be perfect for isotopic mixes that only
  ! have one isotope per element in them), I choose the most abundant
  ! isotope for each element to calculate the relevant mean
  ! equilibrium constants.  And that calculation requires the masses
  ! of the most abundant isotopes for each element which are provided
  ! by this data module.

  ! N.B. These most abundant isotopic weight data *should not be
  ! confused* with the mean atomic weights required for abundance and
  ! density calculations.  The element order is
  ! H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni. We use the
  ! atomic weight scale where (neutral monatomic) C(12) has a weight
  ! of 12.00000000....  All weights are for the neutral monatomic
  ! species of the most abundant isotope.
  ! START OF AME2020 DATA (inserted from the file isotopic_mass.data documented above)
  ! Start of element section.
  real(fp_kind), parameter :: isotopic_mass(nelements_iso) = [&
       ! Z =  1, A =  1, N = A - Z =  0, N - Z = -1
       isotopic_mass_2d(-1, 1),&
       ! Z =  2, A =  4, N = A - Z =  2, N - Z =  0
       isotopic_mass_2d( 0, 2),&
       ! Z =  6, A = 12, N = A - Z =  6, N - Z =  0
       isotopic_mass_2d( 0, 6),&
       ! Z =  7, A = 14, N = A - Z =  7, N - Z =  0
       isotopic_mass_2d( 0, 7),&
       ! Z =  8, A = 16, N = A - Z =  8, N - Z =  0
       isotopic_mass_2d( 0, 8),&
       ! Z = 10, A = 20, N = A - Z = 10, N - Z =  0
       isotopic_mass_2d( 0,10),&
       ! Z = 11, A = 23, N = A - Z = 12, N - Z =  1
       isotopic_mass_2d( 1,11),&
       ! Z = 12, A = 24, N = A - Z = 12, N - Z =  0
       isotopic_mass_2d( 0,12),&
       ! Z = 13, A = 27, N = A - Z = 14, N - Z =  1
       isotopic_mass_2d( 1,13),&
       ! Z = 14, A = 28, N = A - Z = 14, N - Z =  0
       isotopic_mass_2d( 0,14),&
       ! Z = 15, A = 31, N = A - Z = 16, N - Z =  1
       isotopic_mass_2d( 1,15),&
       ! Z = 16, A = 32, N = A - Z = 16, N - Z =  0
       isotopic_mass_2d( 0,16),&
       ! Z = 17, A = 35, N = A - Z = 18, N - Z =  1
       isotopic_mass_2d( 1,17),&
       ! Z = 18, A = 40, N = A - Z = 22, N - Z =  4
       isotopic_mass_2d( 4,18),&
       ! Z = 20, A = 40, N = A - Z = 20, N - Z =  0
       isotopic_mass_2d( 0,20),&
       ! Z = 22, A = 48, N = A - Z = 26, N - Z =  4
       isotopic_mass_2d( 4,22),&
       ! Z = 24, A = 52, N = A - Z = 28, N - Z =  4
       isotopic_mass_2d( 4,24),&
       ! Z = 25, A = 55, N = A - Z = 30, N - Z =  5
       isotopic_mass_2d( 5,25),&
       ! Z = 26, A = 56, N = A - Z = 30, N - Z =  4
       isotopic_mass_2d( 4,26),&
       ! Z = 28, A = 58, N = A - Z = 30, N - Z =  2
       isotopic_mass_2d( 2,28),&
       isotopic_mass_2d(1,3),&
       isotopic_mass_2d(1,4),&
       isotopic_mass_2d(1,5),&
       isotopic_mass_2d(1,9)]
  ! END OF AME2020 DATA (inserted from the file isotopic_mass.data documented above)

end module mod_isotopic_mass_data
