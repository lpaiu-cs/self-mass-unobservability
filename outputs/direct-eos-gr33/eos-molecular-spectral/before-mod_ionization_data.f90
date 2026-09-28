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
!> This module provides physical data describing the ion charge and
!> ionization potential of the 295 neutral and ionized (other than
!> bare ion) species of the chosen 20 elements and H_2 and H_2+.  In
!> addition this module provides the H_2 dissociation energy.
module mod_ionization_data
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public nions, nion, bi, h2diss

  integer, parameter :: nions = 316
  ! charge of parent ion in ion order.
  integer, parameter ::  nion(nions+2) = [&
       ! H
       1,&
       ! He
       1,2,&
       ! C
       1,2,3,4,5,6,&
       ! N
       1,2,3,4,5,6,7,&
       ! O
       1,2,3,4,5,6,7,8,&
       ! Ne
       1,2,3,4,5,6,7,8,9,10,&
       ! Na
       1,2,3,4,5,6,7,8,9,10,11,&
       ! Mg
       1,2,3,4,5,6,7,8,9,10,11,12,&
       ! Al
       1,2,3,4,5,6,7,8,9,10,11,12,13,&
       ! Si
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,&
       ! P
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,&
       ! S
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,&
       ! Cl
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,&
       ! A
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,&
       ! Ca
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
       ! Ti
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
       21,22,&
       ! Cr
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
       21,22,23,24,&
       ! Mn
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
       21,22,23,24,25,&
       ! Fe
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
       21,22,23,24,25,26,&
       ! Ni
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
       21,22,23,24,25,26,27,28,&
       ! Li
       1,2,3,&
       ! Be
       1,2,3,4,&
       ! B
       1,2,3,4,5,&
       ! F
       1,2,3,4,5,6,7,8,9,&
       ! H2
       1,&
       ! H2+
       2]

  ! The energy of all relevant species in units of cm**{-1} (the unit
  ! used for spectroscopic energy [inverse wavelength] measures, and thus
  ! the unit of energy least susceptible to systematic errors due to
  ! changes in adopted physical constants.)
  ! To convert this energy unit to ergs multiply by hc.
  ! Thus, E (in ergs)/kt = c2 E (in cm^{-1})/T where c2 = hc/k = 1.4....
  ! Note the energies are organized in the following order: energy of
  ! first ion relative to neutral, energy of second ion relative to first
  ! ion, energy of third ion relative to second, ...., for each element.
  ! These relative energies are what are measured spectroscopically, but
  ! for EOS purposes we need the energy relative to the neutral monatomic
  ! state (or ultimately the lowest ro-vibrational state of H_2 in the
  ! case of the hydrogen element) and these conversions are done later.

  ! On 2020-03-11 these data were cut and pasted from the NIST form at
  ! <https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html> using all
  ! defaults other than energy in cm^-1, and using the search terms:
  ! H-Ni (i.e., all elements from atomic number 1 through 28 which
  ! includes data from Li, Be, B, F, K, Sc, V, Co that will be filtered
  ! out later).  From data included with these results, these data
  ! should be referenced as follows:
  ! Kramida, A., Ralchenko, Yu.,
  ! Reader, J., and NIST ASD Team (2019). NIST Atomic Spectra Database
  ! (ver. 5.7.1), [Online]. Available: https://physics.nist.gov/asd
  ! [2020, March 11]. National Institute of Standards and Technology,
  ! Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F

  ! Further processing of those cut and pasted data was done as follows:
  ! sed -f free_eos.git/utils/NIST_ionization.sed <20100311_NIST_ionization |\
  ! python3 free_eos.git/utils/NIST_ionization.py >| /tmp/test_irwin
  ! where the sed script parsed the data for each neutral and ion into
  ! lines containing <iatomic>, <ion>, <ip>, where iatomic is the atomic number,
  ! ion is the ion number (starting at zero for neutrals), and ip is the
  ! ionization potential in cm^-1, and the Python script reformatted those data
  ! into the Fortran form used here.
  real(fp_kind), parameter :: monatomic_ip(nions) = [&
       ! H
       109678.77174307_fp_kind,&
       ! He
       198310.66637_fp_kind, 438908.8785_fp_kind,&
       ! missing Li
       ! missing Be
       ! missing B
       ! C
       90820.348_fp_kind, 196663.40_fp_kind, 386241.0_fp_kind, 520175.3_fp_kind, 3162423.3_fp_kind,&
       3952061.67_fp_kind,&
       ! N
       117225.7_fp_kind, 238750.2_fp_kind, 382672._fp_kind, 624866._fp_kind, 789537.2_fp_kind,&
       4452723.3_fp_kind, 5380089.8_fp_kind,&
       ! O
       109837.02_fp_kind, 283270.9_fp_kind, 443085.0_fp_kind, 624382.0_fp_kind, 918657._fp_kind,&
       1114004._fp_kind, 5963073.0_fp_kind, 7028394.7_fp_kind,&
       ! missing F
       ! Ne
       173929.75_fp_kind, 330388.6_fp_kind, 511543.5_fp_kind, 783890._fp_kind, 1018250._fp_kind,&
       1273820._fp_kind, 1671750._fp_kind, 1928447._fp_kind, 9644840.7_fp_kind, 10986877.2_fp_kind,&
       ! Na
       41449.451_fp_kind, 381390.2_fp_kind, 577654._fp_kind, 797970._fp_kind, 1116300._fp_kind,&
       1389110._fp_kind, 1681700._fp_kind, 2130850._fp_kind, 2418500._fp_kind, 11817106.7_fp_kind,&
       13297680._fp_kind,&
       ! Mg
       61671.05_fp_kind, 121267.64_fp_kind, 646402._fp_kind, 881285._fp_kind, 1139940._fp_kind,&
       1506300._fp_kind, 1814870._fp_kind, 2144820._fp_kind, 2645380._fp_kind, 2964000._fp_kind,&
       14209914.7_fp_kind, 15829950._fp_kind,&
       ! Al
       48278.480_fp_kind, 151862.5_fp_kind, 229445.71_fp_kind, 967804._fp_kind, 1240684._fp_kind,&
       1536400._fp_kind, 1949900._fp_kind, 2295760._fp_kind, 2663340._fp_kind, 3215300._fp_kind,&
       3565010._fp_kind, 16824539.3_fp_kind, 18584143._fp_kind,&
       ! Si
       65747.76_fp_kind, 131838.14_fp_kind, 270139.3_fp_kind, 364093.1_fp_kind, 1345070._fp_kind,&
       1655690._fp_kind, 1988730._fp_kind, 2448640._fp_kind, 2833300._fp_kind, 3237350._fp_kind,&
       3841400._fp_kind, 4221630._fp_kind, 19661038.9_fp_kind, 21560631._fp_kind,&
       ! P
       84580.83_fp_kind, 159451.7_fp_kind, 243600.7_fp_kind, 414922.8_fp_kind, 524462.9_fp_kind,&
       1777891._fp_kind, 2125800._fp_kind, 2497100._fp_kind, 3002900._fp_kind, 3423000._fp_kind,&
       3866950._fp_kind, 4521700._fp_kind, 4934020._fp_kind, 22719901.6_fp_kind, 24759942._fp_kind,&
       ! S
       83559.1_fp_kind, 188232.7_fp_kind, 281130._fp_kind, 380870._fp_kind, 585514.1_fp_kind,&
       710194.7_fp_kind, 2266050._fp_kind, 2651900._fp_kind, 3063600._fp_kind, 3611300._fp_kind,&
       4069500._fp_kind, 4552250._fp_kind, 5258400._fp_kind, 5702290._fp_kind, 26001545.1_fp_kind,&
       28182526._fp_kind,&
       ! Cl
       104591.01_fp_kind, 192070.0_fp_kind, 320970._fp_kind, 429430._fp_kind, 545840._fp_kind,&
       781900._fp_kind, 921096._fp_kind, 2809280._fp_kind, 3233080._fp_kind, 3683400._fp_kind,&
       4274500._fp_kind, 4771400._fp_kind, 5293400._fp_kind, 6051000._fp_kind, 6526620._fp_kind,&
       29506532.5_fp_kind, 31828983._fp_kind,&
       ! Ar
       127109.842_fp_kind, 222848.3_fp_kind, 328550._fp_kind, 480560._fp_kind, 603660._fp_kind,&
       736300._fp_kind, 1003450._fp_kind, 1157056._fp_kind, 3408480._fp_kind, 3869500._fp_kind,&
       4358900._fp_kind, 4992200._fp_kind, 5528700._fp_kind, 6090500._fp_kind, 6899800._fp_kind,&
       7407190._fp_kind, 33235410._fp_kind, 35699895._fp_kind,&
       ! missing K
       ! Ca
       49305.9240_fp_kind, 95751.87_fp_kind, 410642.3_fp_kind, 542595._fp_kind, 680230._fp_kind,&
       877400._fp_kind, 1026000._fp_kind, 1187600._fp_kind, 1520640._fp_kind, 1704049._fp_kind,&
       4771550._fp_kind, 5309100._fp_kind, 5876800._fp_kind, 6591300._fp_kind, 7210400._fp_kind,&
       7853400._fp_kind, 8766000._fp_kind, 9337690._fp_kind, 41367028._fp_kind, 44117409._fp_kind,&
       ! missing Sc
       ! Ti
       55072.5_fp_kind, 109494._fp_kind, 221735.6_fp_kind, 348973.3_fp_kind, 800900._fp_kind,&
       964100._fp_kind, 1134660._fp_kind, 1374900._fp_kind, 1549000._fp_kind, 1741500._fp_kind,&
       2137900._fp_kind, 2351109._fp_kind, 6353000._fp_kind, 6968700._fp_kind, 7618000._fp_kind,&
       8408200._fp_kind, 9116000._fp_kind, 9842400._fp_kind, 10858900._fp_kind, 11495470._fp_kind,&
       50401766._fp_kind, 53440740._fp_kind,&
       ! missing V
       ! Cr
       54575.6_fp_kind, 132971.02_fp_kind, 249700._fp_kind, 396500._fp_kind, 560200._fp_kind,&
       731020._fp_kind, 1292830._fp_kind, 1490230._fp_kind, 1690070._fp_kind, 1972200._fp_kind,&
       2184000._fp_kind, 2393400._fp_kind, 2860500._fp_kind, 3098480._fp_kind, 8159300._fp_kind,&
       8849700._fp_kind, 9582000._fp_kind, 10443000._fp_kind, 11247400._fp_kind, 12059000._fp_kind,&
       13180000._fp_kind, 13882280._fp_kind, 60345293._fp_kind, 63675850._fp_kind,&
       ! Mn
       59959.560_fp_kind, 126145.0_fp_kind, 271550._fp_kind, 413000._fp_kind, 584000._fp_kind,&
       771100._fp_kind, 961440._fp_kind, 1576600._fp_kind, 1789630._fp_kind, 2005400._fp_kind,&
       2307900._fp_kind, 2536000._fp_kind, 2771000._fp_kind, 3250000._fp_kind, 3509900._fp_kind,&
       9144100._fp_kind, 9872800._fp_kind, 10649000._fp_kind, 11541000._fp_kind, 12398500._fp_kind,&
       13253400._fp_kind, 14426800._fp_kind, 15162200._fp_kind, 65659877._fp_kind, 69137430._fp_kind,&
       ! Fe
       63737.704_fp_kind, 130655.4_fp_kind, 247220._fp_kind, 442900._fp_kind, 604900._fp_kind,&
       798370._fp_kind, 1008000._fp_kind, 1218380._fp_kind, 1884000._fp_kind, 2114000._fp_kind,&
       2346000._fp_kind, 2668200._fp_kind, 2912000._fp_kind, 3163000._fp_kind, 3679500._fp_kind,&
       3946570._fp_kind, 10184000._fp_kind, 10951200._fp_kind, 11772000._fp_kind, 12708000._fp_kind,&
       13606800._fp_kind, 14505300._fp_kind, 15731400._fp_kind, 16500160._fp_kind, 71204137._fp_kind,&
       74829550._fp_kind,&
       ! missing Co
       ! Ni
       61619.77_fp_kind, 146541.56_fp_kind, 283800._fp_kind, 443000._fp_kind, 613500._fp_kind,&
       871000._fp_kind, 1065000._fp_kind, 1307000._fp_kind, 1558000._fp_kind, 1812000._fp_kind,&
       2577200._fp_kind, 2836060._fp_kind, 3101600._fp_kind, 3462700._fp_kind, 3732400._fp_kind,&
       3995400._fp_kind, 4606000._fp_kind, 4895950._fp_kind, 12429400._fp_kind, 13274000._fp_kind,&
       14183000._fp_kind, 15166000._fp_kind, 16196300._fp_kind, 17183300._fp_kind, 18515300._fp_kind,&
       19351330._fp_kind, 82985464._fp_kind, 86909350._fp_kind,&
       ! Li
       4.34871141979026288e+04_fp_kind,&
       6.10078525778856245e+05_fp_kind,&
       9.87661013882954605e+05_fp_kind,&
       ! Be
       7.51926383991815528e+04_fp_kind,&
       1.46882830474657094e+05_fp_kind,&
       1.24125660321880155e+06_fp_kind,&
       1.75601880998812616e+06_fp_kind,&
       ! B
       6.69280368374585669e+04_fp_kind,&
       2.02887386601550068e+05_fp_kind,&
       3.05930840214578668e+05_fp_kind,&
       2.09199545004716655e+06_fp_kind,&
       2.74410793310331134e+06_fp_kind,&
       ! F
       1.40524520222526597e+05_fp_kind,&
       2.82058604579691193e+05_fp_kind,&
       5.05773967912415625e+05_fp_kind,&
       7.03113792738417513e+05_fp_kind,&
       9.21480329298210097e+05_fp_kind,&
       1.26760596903544711e+06_fp_kind,&
       1.49363227201710106e+06_fp_kind,&
       7.69370865041271877e+06_fp_kind,&
       8.89724290788824484e+06_fp_kind]

  ! H2 I.P. from Huber and Herzberg, footnote b
  real(fp_kind), parameter :: h2_ip = 124417.2_fp_kind
  ! H2 dissociation energy cm^{-1} taken from Huber and Herzberg
  ! footnote a.
  real(fp_kind), parameter :: h2diss = 36118.3_fp_kind

  ! Calculate the ionization potential of H2+ that is
  ! self-consistent with the H2 dissociation energy, h2diss.

  ! Note that H2 dissociates to H + H and ionizes to H2+ which implies
  ! (if
  ! D(H2) is the dissociation energy of H2,
  ! IP(H2) is the ionization energy of H2,
  ! D(H2+) is the dissociation energy of H2+,
  ! IP(H) is the ionization energy of monatomic H,
  ! E(H2) = ground-state energy of H2
  ! E(H + H) energy of dissociated H2
  ! E(H2+) = ground-state energy of H2+
  ! E(H + H+) energy of dissociated H2+
  ! )
  ! that
  ! 1.  D(H2) = E(H + H) - E(H2),
  ! 2.  IP(H2) = E(H2+) - E(H2),
  ! 3.  D(H2+) = E(H + H+) - E(H2+), and
  ! 4.  IP(H) = E(H + H+) - E(H + H).
  ! Equations 4 and 1 imply
  ! 5. E(H + H+) = IP(H) + D(H2) + E(H2),
  ! Equations 2 and 3 imply
  ! 6. D(H2+) = E(H + H+) - (IP(H2) + E(H2)),
  ! and the combination of Equations 5 and 6 imply
  ! 7. D(H2+) =  D(H2) + IP(H) - IP(H2)
  ! Furthermore, note that H2+ dissociates to H + H+ which implies
  ! 8. IP(H2+) = D(H2+) + IP(H)
  ! which in combination with equation 7 implies
  ! 9. IP(H2+) = D(H2) + 2*IP(H) - IP(H2)
  real(fp_kind), parameter :: h2plus_ip = h2diss + 2._fp_kind*monatomic_ip(1) - h2_ip

  integer i
  real(fp_kind), parameter :: bi(nions+2) = [(monatomic_ip(i), i=1,nions), h2_ip, h2plus_ip]

  ! N.B. The experimental change below was tried long ago when bi was
  ! initialized with a data statement, but if you wanted to try it
  ! again (or some variation of this idea) it is easy to do the same
  ! thing for the present implementation where bi is specified with
  ! the parameter attribute.

  ! Experimental change to bi to get energy scale onto the opal system, but
  ! this did not quite work (perhaps because c2 inconsistent with what Forrest
  ! assumed?) so comment out.
  ! Forrest Rogers believed he used 1/k = 11605.4 (where k expressed
  ! in ev/K) and he also used 1 Rydberg = 13.6058ev
  ! put Helium on Rydberg scale
  !use mod_free_eos_constants, only: c2
  !bi(3) = bi(3)/bi(1)
  !bi(2) = bi(2)/bi(1)
  ! 1 Rydberg in cm^-1 under opal energy system
  !bi(1) = 11605.4._fp_kind*13.6058._fp_kind/c2
  !bi(2) = bi(2)*bi(1)
  !bi(3) = bi(3)*bi(1)
end module mod_ionization_data
