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

! calculate approximation to planck-larkin sum where
! sum = approximation to sum from nlow to infinity of
! 2 n^2 (exp(a/n^2) -(1+a/n^2))
! dsum is the derivative of sum wrt to a
! dsum2 is the second derivative of sum wrt to a

!> This plsum_approx subroutine calculates an approximation to the
!> Planck-Larkin sum and the first and second derivative of that
!> quantity wrt to the independent variable a (a convenient
!> transformation of t).  This sum is defined as the sum from n = nlow
!> to n = infinity of 2 n^2 (exp(a/n^2) -(1+a/n^2)).
!>
!> \param[in] nlow PARAMETERS NEED DOCUMENTATION
!>
subroutine plsum_approx(nlow, a, sum, dsum, dsum2)

  ! arguments
  integer, intent(in) :: nlow
  real(fp_kind), intent(in) :: a
  real(fp_kind), intent(out) :: sum, dsum, dsum2

  ! internal variables

  integer, parameter :: nblock = 29
  integer, parameter :: nc = 174

  ! this is a standard value for most of the approximations
  ! which is used to limit their range of applicability
  real(fp_kind),parameter :: exp_max = 1._fp_kind

  integer, parameter :: istart(nblock) = [&
       1,&
       7,&
       13,&
       19,&
       25,&
       31,&
       37,&
       43,&
       49,&
       55,&
       61,&
       67,&
       73,&
       79,&
       85,&
       91,&
       97,&
       103,&
       109,&
       115,&
       121,&
       127,&
       133,&
       139,&
       145,&
       151,&
       157,&
       163,&
       169]
  integer, parameter :: istop(nblock) = [&
       6,&
       12,&
       18,&
       24,&
       30,&
       36,&
       42,&
       48,&
       54,&
       60,&
       66,&
       72,&
       78,&
       84,&
       90,&
       96,&
       102,&
       108,&
       114,&
       120,&
       126,&
       132,&
       138,&
       144,&
       150,&
       156,&
       162,&
       168,&
       174]

  real(fp_kind), parameter :: coeff(nc) = [&
       0.64493379440_fp_kind,&
       0.027442692300_fp_kind,&
       0.0014429412090_fp_kind,&
       0.000069399073200_fp_kind,&
       0.0000023183767810_fp_kind,&
       0.00000016231026170_fp_kind,&
       0.39493394420_fp_kind,&
       0.0066080674370_fp_kind,&
       0.00014296613790_fp_kind,&
       0.0000029087927880_fp_kind,&
       0.000000042230264560_fp_kind,&
       0.0000000012771898410_fp_kind,&
       0.28382288500_fp_kind,&
       0.0024926229090_fp_kind,&
       0.000028822472370_fp_kind,&
       0.00000031735823540_fp_kind,&
       0.0000000025425126030_fp_kind,&
       4.1909882160e-11_fp_kind,&
       0.22132290890_fp_kind,&
       0.0011904792530_fp_kind,&
       0.0000085047618610_fp_kind,&
       0.000000058209565780_fp_kind,&
       0.00000000029397810190_fp_kind,&
       3.0126577460e-12_fp_kind,&
       0.18132292190_fp_kind,&
       0.00065712383250_fp_kind,&
       0.0000031780560560_fp_kind,&
       0.000000014773016170_fp_kind,&
       5.1199129740e-11_fp_kind,&
       3.5536844680e-13_fp_kind,&
       0.15354515200_fp_kind,&
       0.00039991247700_fp_kind,&
       0.0000013940112180_fp_kind,&
       0.0000000046789925720_fp_kind,&
       1.1801540130e-11_fp_kind,&
       5.8924487680e-14_fp_kind,&
       0.13313699390_fp_kind,&
       0.00026107658150_fp_kind,&
       0.00000068646906340_fp_kind,&
       0.0000000017399165930_fp_kind,&
       3.3342434700e-12_fp_kind,&
       1.2521406030e-14_fp_kind,&
       0.11751199740_fp_kind,&
       0.00017969370620_fp_kind,&
       0.00000036891145780_fp_kind,&
       0.00000000073055796580_fp_kind,&
       1.0991856350e-12_fp_kind,&
       3.2128765430e-15_fp_kind,&
       0.10516632100_fp_kind,&
       0.00012888686620_fp_kind,&
       0.00000021226333330_fp_kind,&
       0.00000000033733901520_fp_kind,&
       4.0894862060e-13_fp_kind,&
       9.5590406600e-16_fp_kind,&
       0.095166323060_fp_kind,&
       0.000095552545750_fp_kind,&
       0.00000012901127490_fp_kind,&
       1.6813362860e-10_fp_kind,&
       1.6769544050e-13_fp_kind,&
       3.2040707030e-16_fp_kind,&
       0.086901861790_fp_kind,&
       0.000072784782960_fp_kind,&
       0.000000082016323360_fp_kind,&
       8.9222806310e-11_fp_kind,&
       7.4488930220e-14_fp_kind,&
       1.1845535660e-16_fp_kind,&
       0.079957418560_fp_kind,&
       0.000056709238600_fp_kind,&
       0.000000054133913460_fp_kind,&
       4.9893950530e-11_fp_kind,&
       3.5374324490e-14_fp_kind,&
       4.7534713320e-17_fp_kind,&
       0.074040259790_fp_kind,&
       0.000045038001250_fp_kind,&
       0.000000036884805220_fp_kind,&
       2.9168052560e-11_fp_kind,&
       1.7779035170e-14_fp_kind,&
       2.0449386150e-17_fp_kind,&
       0.068938219800_fp_kind,&
       0.000036360836240_fp_kind,&
       0.000000025827066500_fp_kind,&
       1.7714305280e-11_fp_kind,&
       9.3816209760e-15_fp_kind,&
       9.3391113500e-18_fp_kind,&
       0.064493776040_fp_kind,&
       0.000029776307570_fp_kind,&
       0.000000018517463670_fp_kind,&
       1.1120075040e-11_fp_kind,&
       5.1642625690e-15_fp_kind,&
       4.4922784960e-18_fp_kind,&
       0.060587526620_fp_kind,&
       0.000024689918680_fp_kind,&
       0.000000013554656220_fp_kind,&
       7.1857933550e-12_fp_kind,&
       2.9500334770e-15_fp_kind,&
       2.2613901900e-18_fp_kind,&
       0.057127319510_fp_kind,&
       0.000020698809590_fp_kind,&
       0.000000010105132760_fp_kind,&
       4.7638023940e-12_fp_kind,&
       1.7412365120e-15_fp_kind,&
       1.1850307570e-18_fp_kind,&
       0.054040900200_fp_kind,&
       0.000017523404920_fp_kind,&
       0.0000000076570711610_fp_kind,&
       3.2308497310e-12_fp_kind,&
       1.0581184140e-15_fp_kind,&
       6.4358508690e-19_fp_kind,&
       0.051270817480_fp_kind,&
       0.000014965558120_fp_kind,&
       0.0000000058872066350_fp_kind,&
       2.2362913510e-12_fp_kind,&
       6.5998646280e-16_fp_kind,&
       3.6089599170e-19_fp_kind,&
       0.048770817810_fp_kind,&
       0.000012882176620_fp_kind,&
       0.0000000045861820630_fp_kind,&
       1.5765494200e-12_fp_kind,&
       4.2143940870e-16_fp_kind,&
       2.0829390140e-19_fp_kind,&
       0.046503244410_fp_kind,&
       0.000011168174040_fp_kind,&
       0.0000000036153297030_fp_kind,&
       1.1300522580e-12_fp_kind,&
       2.7489649090e-16_fp_kind,&
       1.2339678830e-19_fp_kind,&
       0.044437128980_fp_kind,&
       0.0000097451972920_fp_kind,&
       0.0000000028809239860_fp_kind,&
       8.2233337340e-13_fp_kind,&
       1.8281106310e-16_fp_kind,&
       7.4858107840e-20_fp_kind,&
       0.042546770050_fp_kind,&
       0.0000085540181700_fp_kind,&
       0.0000000023184427720_fp_kind,&
       6.0671858680e-13_fp_kind,&
       1.2373952120e-16_fp_kind,&
       4.6407612160e-20_fp_kind,&
       0.040810659150_fp_kind,&
       0.0000075493021010_fp_kind,&
       0.0000000018827200040_fp_kind,&
       4.5333472780e-13_fp_kind,&
       8.5123757470e-17_fp_kind,&
       2.9347665050e-20_fp_kind,&
       0.039210659350_fp_kind,&
       0.0000066959501480_fp_kind,&
       0.0000000015416526840_fp_kind,&
       3.4269265320e-13_fp_kind,&
       5.9438683070e-17_fp_kind,&
       1.8901821500e-20_fp_kind,&
       0.037731369590_fp_kind,&
       0.0000059665014670_fp_kind,&
       0.0000000012721008190_fp_kind,&
       2.6185260980e-13_fp_kind,&
       4.2079208910e-17_fp_kind,&
       1.2381417260e-20_fp_kind,&
       0.036359627640_fp_kind,&
       0.0000053392625460_fp_kind,&
       0.0000000010571683890_fp_kind,&
       2.0208332240e-13_fp_kind,&
       3.0172045180e-17_fp_kind,&
       8.2381910580e-21_fp_kind,&
       0.035084117580_fp_kind,&
       0.0000047969422520_fp_kind,&
       0.00000000088437055360_fp_kind,&
       1.5740462100e-13_fp_kind,&
       2.1892057470e-17_fp_kind,&
       5.5616564710e-21_fp_kind,&
       0.033895057080_fp_kind,&
       0.0000043256439000_fp_kind,&
       0.00000000074437948750_fp_kind,&
       1.2366354190e-13_fp_kind,&
       1.6060516410e-17_fp_kind,&
       3.8058628800e-21_fp_kind]

  integer i

  ! Sanity checks
  if(nlow.lt.2.or.nlow.gt.nblock+1) error stop 'plsum_approx: bad nlow value'
  ! limit of approximation
  if(a.gt.exp_max*real((nlow),fp_kind)*real((nlow),fp_kind)) error stop 'plsum_approx: invalid a value'

  ! derivatives evaluated using NR p. 137 synthetic division trick.
  sum = coeff(istop(nlow-1))
  dsum = 0._fp_kind
  dsum2 = 0._fp_kind

  do i = istop(nlow-1)-1, istart(nlow-1), -1
     dsum2 = dsum2*a + dsum
     dsum = dsum*a + sum
     sum = sum*a + coeff(i)
  enddo

  ! multiply power series by a^2
  dsum2 = dsum2*a + dsum
  dsum = dsum*a + sum
  sum = sum*a
  dsum2 = dsum2*a + dsum
  dsum = dsum*a + sum
  sum = sum*a
  ! multiply nth derivative by n factorial
  dsum2 = 2._fp_kind*dsum2
end subroutine plsum_approx
