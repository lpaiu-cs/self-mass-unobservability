! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!

! ***********************************************************************

      module net_derivs_support
      use net_def, only: Net_General_Info, Net_Info
      use const_def
      use crlibm_lib
      use chem_def
      use rates_def

      implicit none

#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif
      
      real(dp), parameter :: r_min = 1d-99

      logical, parameter :: checking_deriv_flags = .false.

      
      logical, parameter :: show_rate = .false.
      logical, parameter :: show_jac = .false.
      logical, parameter :: show_neuQs = .false.
      
      real(dp), parameter :: show_dydt_y = 7.0097206738032283D-02
      logical, parameter :: show_dydt = .false.
      
      logical, parameter :: show_d_dydt_dRho = .false.
      logical, parameter :: show_d_dydt_dT = .false.
      
      logical, parameter :: show_eps_nuc = .false.
      logical, parameter :: show_d_eps_nuc_dy = .false.
      
      logical, parameter :: checkQs = .false.
      
      real(dp), parameter :: checkQ_frac = 1d-4

      logical, parameter :: call_reaclib = .false.
      !logical, parameter :: call_reaclib = .true.
      
      
      contains
      
      
      real(dp) function isoB(ci)
         integer, intent(in) :: ci
         real(dp), parameter :: & ! new names for values defined in chem_def
               deltap     = del_Mp,  &
               deltan     = del_Mn
         isoB = chem_isos% binding_energy(ci) -  &
                chem_isos% Z(ci)*deltap - chem_isos% N(ci)*deltan
      end function isoB


      subroutine do_two_two( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, r_in, c3, i3, c4, i4,  &
            idr1, dr1, idr2, dr2, deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 -> c3 * i3 + c4 * i4
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4 ! net isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         
         real(dp) :: rvs(num_rvs), d, d1, d2, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, cid3, cid4, reaction_id
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0
         
         g => n% g
         chem_id => g% chem_id
         
         d = n% rate_screened(i)
         d1  = dr1 * d
         d2  = dr2 * d
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)

         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_lhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)             
         lhs = lhs + c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c2*rvs(i_rate), -c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_rhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c3*rvs(i_rate), c3, rvs(i_rate), r

         call check(i4, 'i4')
         cid4 = chem_id(i4)
         call do_rhs_iso(n, dydt, i, c4, i4, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         rhs = rhs + c4*(chem_isos% Z(cid4) + chem_isos% N(cid4))
         if (cid4 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c4*rvs(i_rate), c4, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = n% reaction_Qs(reaction_id) - n% reaction_neuQs(reaction_id)
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + n% reaction_neuQs(reaction_id)*rvs(i_rate)
         if (n% reaction_neuQs(reaction_id)*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu',  &
                  n% reaction_neuQs(reaction_id)*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), n% reaction_Qs(reaction_id), &
                  n% reaction_neuQs(reaction_id)
     
         n% d_eps_nuc_dy(idr2) = n% d_eps_nuc_dy(idr2) + d2*Q
         if (chem_id(idr2) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d2*Q, n% d_eps_nuc_dy(idr2), n% reaction_Qs(reaction_id), &
                  n% reaction_neuQs(reaction_id)
                  
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c3*isoB(cid3) + c4*isoB(cid4) - c1*isoB(cid1) - c2*isoB(cid2)
            if (abs(n% reaction_Qs(reaction_id) - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_two_two checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  n% reaction_Qs(reaction_id), checkQ
               !stop
            end if
         end if
         
         if (.false.) then
            write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               n% reaction_Qs(reaction_id), n% reaction_neuQs(reaction_id), rvs(1)*Q, r, rvs(1)
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_two_two: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_two_two


      subroutine do_two_three( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, r_in, c3, i3, c4, i4, c5, i5,  &
            idr1, dr1, idr2, dr2, deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 -> c3 * i3 + c4 * i4 + c5 * i5
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4, i5 ! net isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4, c5 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         
         real(dp) :: rvs(num_rvs), d, d1, d2, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, cid3, cid4, cid5, reaction_id
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))

         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0
         
         g => n % g
         chem_id => g% chem_id
         
         d = n% rate_screened(i)
         d1  = dr1 * d
         d2  = dr2 * d
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)

         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_lhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)             
         lhs = lhs + c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c2*rvs(i_rate), -c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_rhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c3*rvs(i_rate), c3, rvs(i_rate), r

         call check(i4, 'i4')
         cid4 = chem_id(i4)
         call do_rhs_iso(n, dydt, i, c4, i4, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         rhs = rhs + c4*(chem_isos% Z(cid4) + chem_isos% N(cid4))
         if (cid4 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c4*rvs(i_rate), c4, rvs(i_rate), r

         call check(i5, 'i5')
         cid5 = chem_id(i5)
         call do_rhs_iso(n, dydt, i, c5, i5, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         rhs = rhs + c5*(chem_isos% Z(cid5) + chem_isos% N(cid5))
         if (cid5 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c5*rvs(i_rate), c5, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = n% reaction_Qs(reaction_id) - n% reaction_neuQs(reaction_id)
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + n% reaction_neuQs(reaction_id)*rvs(i_rate)
         if (n% reaction_neuQs(reaction_id)*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu',  &
                  n% reaction_neuQs(reaction_id)*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), n% reaction_Qs(reaction_id), &
                  n% reaction_neuQs(reaction_id)
     
         n% d_eps_nuc_dy(idr2) = n% d_eps_nuc_dy(idr2) + d2*Q
         if (chem_id(idr2) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d2*Q, n% d_eps_nuc_dy(idr2), n% reaction_Qs(reaction_id),&
                  n% reaction_neuQs(reaction_id)
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c3*isoB(cid3) + c4*isoB(cid4) + &
                  c5*isoB(cid5) - c1*isoB(cid1) - c2*isoB(cid2)
            if (abs(n% reaction_Qs(reaction_id) - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_two_three checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  n% reaction_Qs(reaction_id), checkQ
               !stop
            end if
         end if
         
         if (.false.) then
            write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               n% reaction_Qs(reaction_id), n% reaction_neuQs(reaction_id), rvs(1)*Q, r, rvs(1)
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_two_three: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_two_three


      subroutine do_three_one( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, c3, i3, r_in, c4, i4,  &
            idr1, dr1, idr2, dr2, deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 + c3 * i3 -> c4 * i4
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         integer :: j
         j = n% g% reaction_id(i)
         call do_three_one_neu( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, c3, i3, r_in, c4, i4, idr1, dr1, idr2, dr2,   &
            n% reaction_Qs(j), n% reaction_neuQs(j),  &
            deriv_flgs, &
            symbolic, just_dydt)
      end subroutine do_three_one


      subroutine do_three_one_neu( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, c3, i3, r_in, c4, i4,  &
            idr1, dr1, idr2, dr2, reaction_Q, Qneu,  &
            deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 + c3 * i3 -> c4 * i4
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         real(dp), intent(in) :: reaction_Q, Qneu
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt

         real(dp) :: rvs(num_rvs), d, d1, d2, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, cid3, cid4, reaction_id
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0

         g => n % g
         chem_id => g% chem_id
         
         d = n% rate_screened(i)
         d1  = dr1 * d
         d2  = dr2 * d
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)
         
         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_lhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)             
         lhs = lhs + c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c2*rvs(i_rate), -c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_lhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)       
         lhs = lhs + c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c3*rvs(i_rate), -c3, rvs(i_rate), r

         call check(i4, 'i4')
         cid4 = chem_id(i4)
         call do_rhs_iso(n, dydt, i, c4, i4, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         rhs = c4*(chem_isos% Z(cid4) + chem_isos% N(cid4))
         if (cid4 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c4*rvs(i_rate), c4, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         
         Q = reaction_Q - Qneu
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + Qneu*rvs(i_rate)
         if (n% reaction_neuQs(reaction_id)*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu', Qneu*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), reaction_Q, Qneu

         n% d_eps_nuc_dy(idr2) = n% d_eps_nuc_dy(idr2) + d2*Q
         if (chem_id(idr2) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d2*Q, n% d_eps_nuc_dy(idr2), reaction_Q, Qneu
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c4*isoB(cid4) - c1*isoB(cid1) - c2*isoB(cid2) - c3*isoB(cid3)
            if (abs(n% reaction_Qs(reaction_id) - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_three_one checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  n% reaction_Qs(reaction_id), checkQ
               !stop
            end if
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_three_one: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_three_one_neu


      subroutine do_three_two( &
               n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, c3, i3, r_in, c4, i4, c5, i5, &
               idr1, dr1, idr2, dr2, idr3, dr3, deriv_flgs, &
               symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 + c3 * i3 -> c4 * i4 + c5 * i5
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4, i5 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4, c5 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         integer, intent(in) :: idr3 ! isotope number for dr3
         real(dp), intent(in) :: dr3 ! coefficient for Jacobian entries d_dydt_dy(idr3)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt

         real(dp) :: rvs(num_rvs), d, d1, d2, d3, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, cid3, cid4, cid5, reaction_id
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0

         g => n % g
         chem_id => g% chem_id
         
         d = n% rate_screened(i)
         d1  = dr1 * d
         d2  = dr2 * d
         d3  = dr3 * d
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)
         
         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, idr2, d2, idr3, d3, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_lhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, idr2, d2, idr3, d3, &
               symbolic, just_dydt)             
         lhs = lhs + c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c2*rvs(i_rate), -c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_lhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, idr2, d2, idr3, d3, &
               symbolic, just_dydt)       
         lhs = lhs + c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c3*rvs(i_rate), -c3, rvs(i_rate), r


         call check(i4, 'i4')
         cid4 = chem_id(i4)
         call do_rhs_iso(n, dydt, i, c4, i4, rvs, idr1, d1, idr2, d2, idr3, d3, &
               symbolic, just_dydt)
         rhs = c4*(chem_isos% Z(cid4) + chem_isos% N(cid4))
         if (cid4 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c4*rvs(i_rate), c4, rvs(i_rate), r

         call check(i5, 'i5')
         cid5 = chem_id(i5)
         call do_rhs_iso(n, dydt, i, c5, i5, rvs, idr1, d1, idr2, d2, idr3, d3, &
               symbolic, just_dydt)
         rhs = rhs + c5*(chem_isos% Z(cid5) + chem_isos% N(cid5))
         if (cid5 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c5*rvs(i_rate), c5, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = n% reaction_Qs(reaction_id) - n% reaction_neuQs(reaction_id)
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + n% reaction_neuQs(reaction_id)*rvs(i_rate)
         if (n% reaction_neuQs(reaction_id)*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu',  &
                           n% reaction_neuQs(reaction_id)*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), n% reaction_Qs(reaction_id), &
                  n% reaction_neuQs(reaction_id)

         n% d_eps_nuc_dy(idr2) = n% d_eps_nuc_dy(idr2) + d2*Q
         if (chem_id(idr2) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d2*Q, n% d_eps_nuc_dy(idr2), n% reaction_Qs(reaction_id), &
                  n% reaction_neuQs(reaction_id)
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = &
                  c4*isoB(cid4) + c5*isoB(cid5) - &
                  c1*isoB(cid1) - c2*isoB(cid2) - c3*isoB(cid3)
            if (abs(n% reaction_Qs(reaction_id) - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_three_two checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  n% reaction_Qs(reaction_id), checkQ
               !stop
            end if
         end if
         
         if (.false.) then
            write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               n% reaction_Qs(reaction_id), n% reaction_neuQs(reaction_id), rvs(1)*Q, r, rvs(1)
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_three_two: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_three_two


      subroutine do_two_one( &
               n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, r_in, c3, i3, idr1, dr1, idr2, dr2,  &
               deriv_flgs, &
               symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 -> c3 * i3
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         integer :: j
         j = n% g% reaction_id(i)
         call do_two_one_neu( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, r_in, c3, i3, idr1, dr1, idr2, dr2,  &
            n% reaction_Qs(j), n% reaction_neuQs(j),  &
            deriv_flgs, &
            symbolic, just_dydt)
      end subroutine do_two_one
      

      subroutine do_two_one_neu( &
            n, dydt, eps_nuc_MeV, i, c1, i1, c2, i2, r_in, c3, i3, idr1, dr1, idr2, dr2,  &
            reaction_Q, Qneu, deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il + c2 * i2 -> c3 * i3
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! coefficient for Jacobian entries d_dydt_dy(idr1)
         integer, intent(in) :: idr2 ! isotope number for dr2
         real(dp), intent(in) :: dr2 ! coefficient for Jacobian entries d_dydt_dy(idr2)
         real(dp), intent(in) :: reaction_Q, Qneu
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         
         real(dp) :: rvs(num_rvs), d, d1, d2, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: ir, j, cid1, icat, cid2, cid3, reaction_id
         
         include 'formats.dek'
         
         r = r_in
         ir = n% g% reaction_id(i)
         icat = reaction_categories(ir)
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0

         g => n % g
         chem_id => g% chem_id
         
         d = n% rate_screened(i)
         d1  = dr1 * d
         d2  = dr2 * d
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)
         
         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_lhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)             
         lhs = lhs + c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c2*rvs(i_rate), -c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_rhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, idr2, d2, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c3*rvs(i_rate), c3, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = reaction_Q - Qneu
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + Qneu*rvs(i_rate)

         if (Qneu*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu', Qneu*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), reaction_Q, Qneu

         n% d_eps_nuc_dy(idr2) = n% d_eps_nuc_dy(idr2) + d2*Q
         if (chem_id(idr2) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d2*Q, n% d_eps_nuc_dy(idr2), reaction_Q, Qneu
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c3*isoB(cid3) - c1*isoB(cid1) - c2*isoB(cid2)
            if (abs(n% reaction_Qs(reaction_id) - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_two_one checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  n% reaction_Qs(reaction_id), checkQ
               !stop
            end if
         end if
         
         
         if (.false. .and. ir == irn14ag_lite) then
            write(*,1) trim(reaction_Name(ir)) // ' epsnuc',  &
               Q*rvs(i_rate), reaction_Q*rvs(i_rate), Qneu, reaction_Q
         end if
         
         if (.false. .and. reaction_name(ir) == 'r_c12_ag_o16' .and. &
               n% y(i1) == 7.6682722655358818D-02) then ! .and. &
!$omp critical               

               write(*,*)
               write(*,2) trim(reaction_name(ir))
               !write(*,2) trim(reaction_name(r_ir))
               write(*,*)
               write(*,1) 'c1', c1
               write(*,1) 'c2', c2
               write(*,1) 'c3', c3
               write(*,*)
               write(*,1) 'd1', d1
               write(*,1) 'd2', d2
               write(*,*)
               write(*,1) 'c1*d1', c1*d1
               write(*,1) 'c1*d2', c1*d2
               write(*,*)
               write(*,1) 'c2*d1', c2*d1
               write(*,1) 'c2*d2', c2*d2
               write(*,*)
               write(*,1) 'c3*d1', c3*d1
               write(*,1) 'c3*d2', c3*d2
               write(*,*)
               write(*,1) 'y1', n% y(i1)
               write(*,1) 'y2', n% y(i2)
               write(*,1) 'y3', n% y(i3)
               write(*,*)
               write(*,1) 'n% rate_screened(i)', n% rate_screened(i)
               write(*,*)


               write(*,*)
               stop 'do_two_one_neu'
!$omp end critical

         end if
         
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_two_one: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check

      end subroutine do_two_one_neu


      subroutine do_one_three( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, c3, i3, c4, i4, idr1, dr1,  &
               deriv_flgs, &
               symbolic, just_dydt)
         ! reaction of form c1 * il -> c2 * i2 + c3 * i3 + c4 * i4
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! rate info for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! for Jacobian entries d_dydt_dy(idr1)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         integer :: j
         j = n% g% reaction_id(i)
         call do_one_three_neu( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, c3, i3, c4, i4, idr1, dr1,  &
               n% reaction_Qs(j), n% reaction_neuQs(j),  &
               deriv_flgs, &
               symbolic, just_dydt)
      end subroutine do_one_three


      subroutine do_one_three_neu( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, c3, i3, c4, i4, idr1, dr1,  &
               reaction_Q, Qneu, deriv_flgs, &
               symbolic, just_dydt)
         ! reaction of form c1 * il -> c2 * i2 + c3 * i3 + c4 * i4
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3, i4 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3, c4 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! rate info for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! for Jacobian entries d_dydt_dy(idr1)
         real(dp), intent(in) :: reaction_Q, Qneu
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         
         real(dp) :: rvs(num_rvs), d1, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, cid3, cid4, reaction_id
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0

         g => n % g
         chem_id => g% chem_id

         d1  = dr1 * n% rate_screened(i)
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)
         
         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_rhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)
         rhs = c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c2*rvs(i_rate), c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_rhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = rhs + c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c3*rvs(i_rate), c3, rvs(i_rate), r

         call check(i4, 'i4')
         cid4 = chem_id(i4)
         call do_rhs_iso(n, dydt, i, c4, i4, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = rhs + c4*(chem_isos% Z(cid4) + chem_isos% N(cid4))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c4*rvs(i_rate), c4, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = reaction_Q - Qneu
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + Qneu*rvs(i_rate)
         if (Qneu*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu', Qneu*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), reaction_Q, Qneu
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c2*isoB(cid2) + c3*isoB(cid3) + c4*isoB(cid4) - c1*isoB(cid1)
            if (abs(reaction_Q - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) ' checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  reaction_Q, checkQ
               !stop
            end if
         end if
         
         if (.false.) then
            write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               reaction_Q, Qneu, rvs(1)*Q, r, rvs(1)
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_one_three_neu


      subroutine do_one_two( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, c3, i3, idr1, dr1,  &
               deriv_flgs, &
               symbolic, just_dydt)
         ! reaction of form c1 * il -> c2 * i2 + c3 * i3
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! rate info for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! for Jacobian entries d_dydt_dy(idr1)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         integer :: j
         j = n% g% reaction_id(i)
         call do_one_two_neu( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, c3, i3, idr1, dr1,  &
               n% reaction_Qs(j), n% reaction_neuQs(j),  &
               deriv_flgs, &
               symbolic, just_dydt)
      end subroutine do_one_two
      

      subroutine do_one_two_neu( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, c3, i3, idr1, dr1,  &
               reaction_Q, Qneu, deriv_flgs, &
               symbolic, just_dydt)
         ! reaction of form c1 * il -> c2 * i2 + c3 * i3
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2, i3 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2, c3 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! rate info for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! for Jacobian entries d_dydt_dy(idr1)
         real(dp), intent(in) :: reaction_Q, Qneu
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         
         real(dp) :: rvs(num_rvs), d1, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, cid3, reaction_id, ir
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0

         g => n % g
         chem_id => g% chem_id

         d1  = dr1 * n% rate_screened(i)
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)
         
         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_rhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)
         rhs = c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c2*rvs(i_rate), c2, rvs(i_rate), r

         call check(i3, 'i3')
         cid3 = chem_id(i3)
         call do_rhs_iso(n, dydt, i, c3, i3, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = rhs + c3*(chem_isos% Z(cid3) + chem_isos% N(cid3))
         if (cid3 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c3*rvs(i_rate), c3, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = reaction_Q - Qneu
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + Qneu*rvs(i_rate)
         if (Qneu*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu', Qneu*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), reaction_Q, Qneu
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c2*isoB(cid2) + c3*isoB(cid3) - c1*isoB(cid1)
            if (abs(reaction_Q - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_one_two checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                  n% reaction_Qs(reaction_id), checkQ
               !stop
            end if
         end if
         
         ir = n% g% reaction_id(i)
         if (.false. .and. reaction_name(ir) == 'r_o16_ga_c12' .and. &
               n% y(i1) == 6.8716855628023080D-04) then ! .and. &
               !n% y(i1) > 1d-20 .and. &
               !n% y(i2) > 1d-20 .and. &
               !n% y(i3) > 1d-20) then
               
               write(*,*)
               write(*,2) trim(reaction_name(ir))
               write(*,*)
               write(*,1) 'c1', c1
               write(*,1) 'c2', c2
               write(*,1) 'c3', c3
               write(*,*)
               write(*,1) 'd1', d1
               write(*,*)
               write(*,1) 'y1',  n% y(i1)
               write(*,1) 'y2',  n% y(i2)
               write(*,1) 'y3',  n% y(i3)
               write(*,*)
               write(*,1) 'c1*d1', c1*d1
               write(*,*)
               write(*,1) 'c2*d1', c2*d1
               write(*,*)
               write(*,1) 'c3*d1', c3*d1
               write(*,*)

               stop 'do_one_two_neu'

         end if
         
                  
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_one_two: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_one_two_neu


      subroutine do_one_one( &
            n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, idr1, dr1,  &
            deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il -> c2 * i2
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! rate info for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! for Jacobian entries d_dydt_dy(idr1)
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         integer :: j
         j = n% g% reaction_id(i)
         call do_one_one_neu( &
               n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, idr1, dr1,  &
               n% reaction_Qs(j), n% reaction_neuQs(j),  &
               deriv_flgs, &
               symbolic, just_dydt)
      end subroutine do_one_one


      subroutine do_one_one_neu( &
            n, dydt, eps_nuc_MeV, i, c1, i1, r_in, c2, i2, idr1, dr1,  &
            reaction_Q, Qneu, deriv_flgs, &
            symbolic, just_dydt)
         ! reaction of form c1 * il -> c2 * i2
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
         integer, intent(in) :: i ! the reaction number
         integer, intent(in) :: i1, i2 ! isotope numbers for the reaction
         real(dp), intent(in) :: c1, c2 ! isotope coefficients in reaction equation
         real(dp), intent(in) :: r_in ! rate info for the reaction
         integer, intent(in) :: idr1 ! isotope number for dr1
         real(dp), intent(in) :: dr1 ! for Jacobian entries d_dydt_dy(idr1)
         real(dp), intent(in) :: reaction_Q, Qneu
         logical, pointer :: deriv_flgs(:)
         logical, intent(in) :: symbolic, just_dydt
         
         real(dp) :: rvs(num_rvs), d1, lhs, rhs, r, Q, checkQ
         type (Net_General_Info), pointer  :: g
         integer, pointer :: chem_id(:)
         integer :: j, cid1, icat, cid2, reaction_id
         logical :: dbg
         
         include 'formats.dek'
         
         r = r_in
         icat = reaction_categories(n% g% reaction_id(i))
      
         if (r < r_min .or. n% rate_screened(i) < r_min) r = 0
         
         g => n % g
         chem_id => g% chem_id
         
      	dbg = .false. ! (reaction_Name(g% reaction_id(i)) == 'r_h1_h1_wk_h2')

         d1  = dr1 * n% rate_screened(i)
         rvs(i_rate) = r * n% rate_screened(i)
         rvs(i_rate_dT) = r * n% rate_screened_dT(i)
         rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)
         
         call check(i1, 'i1')
         cid1 = chem_id(i1)
         call do_lhs_iso(n, dydt, i, c1, i1, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)
         lhs = c1*(chem_isos% Z(cid1) + chem_isos% N(cid1))
         if (cid1 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               -c1*rvs(i_rate), -c1, rvs(i_rate), r

         call check(i2, 'i2')
         cid2 = chem_id(i2)
         call do_rhs_iso(n, dydt, i, c2, i2, rvs, idr1, d1, 0, 0d0, 0, 0d0, &
               symbolic, just_dydt)       
         rhs = c2*(chem_isos% Z(cid2) + chem_isos% N(cid2))
         if (cid2 == ihe4 .and. show_rate) write(*,1) trim(reaction_Name(g% reaction_id(i))),  &
               c2*rvs(i_rate), c2, rvs(i_rate), r
         
         reaction_id = g% reaction_id(i)
         Q = reaction_Q - Qneu
         do j=1,num_rvs
            eps_nuc_MeV(j) = eps_nuc_MeV(j) + rvs(j)*Q
         end do
         if (show_eps_nuc .and. abs(rvs(1)*Q) > 1d2) &
            write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q
         n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*Q
         n% eps_neu_total = n% eps_neu_total + Qneu*rvs(i_rate)
         if (Qneu*rvs(i_rate) /= 0 .and. show_neuQs .and. &
                  abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu', Qneu*rvs(:)
         
         n% d_eps_nuc_dy(idr1) = n% d_eps_nuc_dy(idr1) + d1*Q
         if (chem_id(idr1) == ini56 .and. show_d_eps_nuc_dy)  &
               write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                  d1*Q, n% d_eps_nuc_dy(idr1), reaction_Q, Qneu
         
         call check_balance(n, i, lhs, rhs)
         if (checking_deriv_flags) deriv_flgs(i) = .true.
         
         if (checkQs) then
            checkQ = c2*isoB(cid2) - c1*isoB(cid1)
            if (abs(n% reaction_Qs(reaction_id) - checkQ) > checkQ_frac*abs(checkQ)) then
               write(*,1) 'do_one_one_neu checkQ ' // trim(reaction_Name(g% reaction_id(i))), &
                  n% reaction_Qs(reaction_id), checkQ, c1, isoB(cid1), c2, isoB(cid2)
               !stop
            end if
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: check
#endif         
         subroutine check(ii, str)
            integer, intent(in) :: ii
            character (len=*), intent(in) :: str
            if (ii <= 0) then
               write(*,*)  &
                  'do_one_one_neu: bad iso num for ' // trim(str) // &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
               stop 1
            end if
         end subroutine check
      
      end subroutine do_one_one_neu
      
         
      subroutine do_lhs_iso( &
            n, dydt, i, c, i1, rvs, i2, dr2, i3, dr3, i4, dr4, symbolic, just_dydt)
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         integer, intent(in) :: i, i1, i2, i3, i4
         real(dp), intent(in) :: c, rvs(:), dr2, dr3, dr4
         logical, intent(in) :: symbolic, just_dydt
         
         ! i1, i2, i3, and 14 are isotope numbers
         ! dr2 = dr/dy(i2); dr3 = dr/dy(i3); dr4 = dr/dy(i4)
         ! -c * r   = dydt(i1)
         ! -c * dr2 = d_dydt(i1)_dy(i2)
         ! -c * dr3 = d_dydt(i1)_dy(i3)
         ! -c * dr4 = d_dydt(i1)_dy(i4)
         
         integer :: j
         integer, pointer :: chem_id(:)
         chem_id => n% g% chem_id
         
         include 'formats.dek'
         
         if (symbolic) then
            n% d_dydt_dy(i1, i2) = 1
            if (i3 <= 0) return
            n% d_dydt_dy(i1, i3) = 1
            if (i4 <= 0) return
            n% d_dydt_dy(i1, i4) = 1
            return
         end if

         ! update the dydt terms for i1
         do j=1,num_rvs
            dydt(j,i1) = dydt(j,i1) - c * rvs(j)
         end do
         if (chem_id(i1) == io16 .and. show_dydt .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) 'lhs ' // trim(reaction_Name(n% g% reaction_id(i))), &
                  -c * rvs(i_rate), dydt(i_rate,i1), &
                  n% rate_screened(i), &
                  n% rate_raw(i), &
                  n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dRho .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dRho ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  -c * rvs(i_rate_dRho), n% rate_screened(i), &
                  n% rate_screened_dRho(i), n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dT .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dT ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  -c * rvs(i_rate_dT), n% rate_screened(i), &
                  n% rate_screened_dT(i), n% y(i1)
         
         if (just_dydt) return
         
         ! update the Jacobian for d_dydt(i1)_dy(i2)
         n% d_dydt_dy(i1, i2) = n% d_dydt_dy(i1, i2)  - c * dr2
         
         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr2 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy l2 dr2 '  &
                  // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i2))),  &
                  - c * dr2, n% d_dydt_dy(i1, i2)
            
         if (i3 <= 0) return
            
         ! update the Jacobian for d_dydt(i1)_dy(i3)
         n% d_dydt_dy(i1, i3) = n% d_dydt_dy(i1, i3)  - c * dr3
         
         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr3 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy l3 dr3 '  &
           // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i3))),  &
           - c * dr3, n% d_dydt_dy(i1, i3)
            
         if (i4 <= 0) return
            
         !write(*,4) trim(reaction_Name(n% g% reaction_id(i))) // ' do_lhs_iso', i, i1, i4
         ! update the Jacobian for d_dydt(i1)_dy(i4)
         n% d_dydt_dy(i1, i4) = n% d_dydt_dy(i1, i4)  - c * dr4
         
         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr4 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy l4 dr4 '  &
           // trim(chem_isos% name(chem_id(i1))) // ' ' // &
                  trim(chem_isos% name(chem_id(i4))),  &
                  - c * dr4, n% d_dydt_dy(i1, i4)
         
      end subroutine do_lhs_iso
            

      subroutine do_rhs_iso( &
            n, dydt, i, c, i1, rvs, i2, dr2, i3, dr3, i4, dr4, symbolic, just_dydt)
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         integer, intent(in) :: i, i1, i2, i3, i4
         real(dp), intent(in) :: c, rvs(:), dr2, dr3, dr4
         logical, intent(in) :: symbolic, just_dydt
         
         ! i1, i2, i3, and 14 are isotope numbers
         ! dr2 = dr/dy(i2); dr3 = dr/dy(i3); dr4 = dr/dy(i4)
         ! c * r   = dydt(i1)
         ! c * dr2 = d_dydt(i1)_dy(i2)
         ! c * dr3 = d_dydt(i1)_dy(i3)
         ! c * dr4 = d_dydt(i1)_dy(i4)
         
         integer :: j
         integer, pointer :: chem_id(:)
         chem_id => n% g% chem_id
         
         include 'formats.dek'
         
         if (symbolic) then
            n% d_dydt_dy(i1, i2) = 1
            if (i3 <= 0) return
            n% d_dydt_dy(i1, i3) = 1
            if (i4 <= 0) return
            n% d_dydt_dy(i1, i4) = 1
            return
         end if

         ! update the dydt terms for i1
         do j=1,num_rvs
            dydt(j,i1) = dydt(j,i1) + c * rvs(j)
         end do
         if (chem_id(i1) == io16 .and. show_dydt .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) 'rhs ' // trim(reaction_Name(n% g% reaction_id(i))), &
               c * rvs(i_rate), dydt(i_rate,i1), &
               n% rate_screened(i), n% rate_raw(i), n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dRho .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dRho ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  c * rvs(i_rate_dRho), n% rate_screened(i), &
                  n% rate_screened_dRho(i), n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dT .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dT ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  c * rvs(i_rate_dT), n% rate_screened(i), &
                  n% rate_screened_dT(i), n% y(i1)
         
         if (just_dydt) return

         ! update the Jacobian for d_dydt(i1)_dy(i2)
         n% d_dydt_dy(i1, i2) = n% d_dydt_dy(i1, i2)  + c * dr2
         
         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr2 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy r2 dr2 ' &
                  // trim(chem_isos% name(chem_id(i1))) // ' ' // &
                  trim(chem_isos% name(chem_id(i2))),  &
                  c * dr2, n% d_dydt_dy(i1, i2)
            
         if (i3 <= 0) return
            
         ! update the Jacobian for d_dydt(i1)_dy(i3)
         n% d_dydt_dy(i1, i3) = n% d_dydt_dy(i1, i3)  + c * dr3
         
         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr3 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy r3 dr3 ' &
            // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i3))),  &
           c * dr3, n% d_dydt_dy(i1, i3)
            
         if (i4 <= 0) return
            
         ! update the Jacobian for d_dydt(i1)_dy(i4)
         n% d_dydt_dy(i1, i4) = n% d_dydt_dy(i1, i4)  + c * dr4
         
         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr4 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy r4 dr4 ' &
            // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i4))),  &
           c * dr4, n% d_dydt_dy(i1, i4)

      end subroutine do_rhs_iso
      
         
      subroutine check_balance(n, i, lhs, rhs) ! check conservation of nucleons
         type (Net_Info), pointer :: n
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         integer, intent(in) :: i
         real(dp), intent(in) :: lhs, rhs
         if (lhs == rhs) return
         if (abs(lhs-rhs) < 1d-6) return
         write(*,'(2a)') 'non-conservation of nucleons in reaction ',  &
               reaction_Name(n% g% reaction_id(i))
         write(*,*) 'lhs aion sum', lhs
         write(*,*) 'rhs aion sum', rhs
         stop
      end subroutine check_balance


#ifdef offload
      !dir$ end options
#endif

      end module net_derivs_support



