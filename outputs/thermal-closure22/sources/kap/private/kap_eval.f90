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

      module kap_eval
      use utils_lib,only:is_bad_num
      use kap_def
      use const_def, only: dp
      use crlibm_lib
      
      implicit none
      
#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif            
      
            
      contains
      
      
      real(dp) function get_logT_Compton_blend_hi()
         get_logT_Compton_blend_hi = kap_z_tables(1)% x_tables(1)% logT_max - 0.01d0
      end function get_logT_Compton_blend_hi


      subroutine Get_kap_Results( &
            rq, zbar, X, dXC, dXO, Zbase, Z, frac_Type1, &
            Rho_in, logRho_in, T_in, logT_in, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use const_def
         use kap_eval_co
         use kap_eval_fixed
         use utils_lib, only: is_bad_num, is_bad_real
         
         type (Kap_General_Info), pointer :: rq
         real, intent(in) :: X, dXC, dXO, Zbase, Z, zbar
         real(dp), intent(in) :: frac_Type1
         real, intent(in) :: Rho_in, logRho_in, T_in, logT_in
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         real :: logR, logT, T, logRho, Rho
         real(dp) :: logT_Compton_blend_hi, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            beta, kap_beta, dlnkap_beta_dlnRho, dlnkap_beta_dlnT, &
            alfa, kap_alfa, dlnkap_alfa_dlnRho, dlnkap_alfa_dlnT
         
         logical, parameter :: dbg = .false.

         include 'formats.dek'
         
         logRho = logRho_in; Rho = Rho_in
         logT = logT_in; T = T_in
         logR = logRho - 3*logT + 18
         if (dbg) write(*,1) 'logRho', logRho
         if (dbg) write(*,1) 'logT', logT
         if (dbg) write(*,1) 'logR', logR

            
         if (is_bad_real(dXC) .or. is_bad_real(dXO)) then
            write(*,1) 'Get_kap_Results dXC', dXC
            write(*,1) 'Get_kap_Results dXO', dXO
            write(*,1) 'X', X
            write(*,1) 'Z', Z
            write(*,1) 'Zbase', Zbase
            write(*,1) 'zbar', zbar
            write(*,1) 'logRho', logRho
            write(*,1) 'logT', logT
            write(*,1) 'logR', logR
            write(*,*)
            ierr = -1
            return
            stop 'Get_kap_Results'
         end if

         beta = frac_Type1
         
         if (logT <= kap_type2_logT_lower_bdy) beta = 1 ! all lowT means all type1
         
         alfa = 1 - beta
         if (dbg) write(*,1) 'alfa', alfa
         if (dbg) write(*,1) 'beta', beta
         
         logT_Compton_blend_hi = get_logT_Compton_blend_hi()
         if (logT >= logT_Compton_blend_hi) then ! just use compton

            kap_rad = 0
            dlnkap_rad_dlnRho = 0
            dlnkap_rad_dlnT = 0

         else

            if (beta > 0) then ! get value from fixed metal tables
               call Get_kap_fixed_metal_Results( &
                  rq, zbar, Z, X, Rho, logRho, T, logT, &
                  kap_beta, dlnkap_beta_dlnRho, dlnkap_beta_dlnT, ierr)
               if (dbg) then
                  write(*,1) 'Get_kap_fixed_metal_Results'
                  write(*,1) 'zbar', zbar
                  write(*,1) 'Z', Z
                  write(*,1) 'X', X
                  write(*,1) 'Rho', Rho
                  write(*,1) 'logRho', logRho
                  write(*,1) 'T', T
                  write(*,1) 'logT', logT
                  write(*,1) 'kap_beta', kap_beta
                  write(*,*)
               end if
            else
               kap_beta = 0; dlnkap_beta_dlnRho = 0; dlnkap_beta_dlnT = 0
            end if
            
            if (alfa > 0) then ! get value from C/O enhanced tables
               call Get_kap_CO_Results( &
                  rq, zbar, Zbase, X, dXC, dXO, Rho, logRho, T, logT, &
                  kap_alfa, dlnkap_alfa_dlnRho, dlnkap_alfa_dlnT, ierr)
            else
               kap_alfa = 0; dlnkap_alfa_dlnRho = 0; dlnkap_alfa_dlnT = 0
            end if
         
            kap_rad = alfa*kap_alfa + beta*kap_beta
            if (dbg) then
               write(*,1) 'alfa', alfa
               write(*,1) 'kap_alfa (Type2)', kap_alfa
               write(*,1) 'beta', beta
               write(*,1) 'kap_beta (Type1)', kap_beta
               write(*,1) 'kap_rad', kap_rad
            end if
            if (kap_rad < 1e-30) then
               kap_rad = 1e-30
               dlnkap_rad_dlnRho = 0
               dlnkap_rad_dlnT = 0
            else
               dlnkap_rad_dlnRho = alfa*dlnkap_alfa_dlnRho + beta*dlnkap_beta_dlnRho
               dlnkap_rad_dlnT = alfa*dlnkap_alfa_dlnT + beta*dlnkap_beta_dlnT
            end if
         
         end if
         
         call combine_rad_with_compton_and_conduction( &
            rq, Rho, logRho, T, logT, zbar, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine Get_kap_Results
      
      
      subroutine combine_rad_with_compton_and_conduction( &
            rq, Rho, logRho, T, logT, zbar, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         use condint, only: do_electron_conduction
         type (Kap_General_Info), pointer :: rq
         real, intent(in) :: Rho, logRho, T, logT, zbar
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(inout) :: kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
      
         real(dp) :: &
            dkap_dlnRho, dkap_dlnT, dkap_rad_dlnRho, dkap_rad_dlnT, &
            kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, dkap_ec_dlnRho, dkap_ec_dlnT, &
            kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT, kap_highT, &
            dkap_compton_dlnRho, dkap_compton_dlnT, fac
         
         real(dp) :: logT_Compton_blend_lo, logT_Compton_blend_hi
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         logT_Compton_blend_hi = get_logT_Compton_blend_hi()
         logT_Compton_blend_lo = logT_Compton_blend_hi - 0.50d0
         if (logT_Compton_blend_hi < 7.5d0) then
            write(*,1) 'Get_kap_Results: logT_Compton_blend_hi', logT_Compton_blend_hi
            write(*,1) 'bogus??? expect 7.99 for OP and 8.69 for OPAL'
            stop 1
         end if
         
         if (logT > logT_Compton_blend_lo) then ! combine kap_rad with rad_compton
            
            call Compton_Opacity( &
               dble(Rho), dble(T), lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT, ierr)
            if (ierr /= 0) return
            if (dbg) write(*,1) 'kap_compton', kap_compton
            if (logT >= logT_Compton_blend_hi) then
               kap_rad = kap_compton
               dlnkap_rad_dlnRho = dlnkap_compton_dlnRho
               dlnkap_rad_dlnT = dlnkap_compton_dlnT
            else
               fac = (logT - logT_Compton_blend_lo) / (logT_Compton_blend_hi - logT_Compton_blend_lo)
               !fac = 0.5 * (1.0 - cos(fac * pi))
               kap_highT = fac*kap_compton + (1-fac)*kap_rad
               dkap_compton_dlnRho = kap_compton*dlnkap_compton_dlnRho
               dkap_rad_dlnRho = kap_rad*dlnkap_rad_dlnRho
               dlnkap_rad_dlnRho = (fac*dkap_compton_dlnRho + (1-fac)*dkap_rad_dlnRho)/kap_highT
               dkap_compton_dlnT = kap_compton*dlnkap_compton_dlnT
               dkap_rad_dlnT = kap_rad*dlnkap_rad_dlnT
               dlnkap_rad_dlnT = (fac*dkap_compton_dlnT + (1-fac)*dkap_rad_dlnT)/kap_highT
               kap_rad = kap_highT
            end if
            
         end if
         
         if (.not. rq% include_electron_conduction) then
            kap = kap_rad
            dlnkap_dlnRho = dlnkap_rad_dlnRho
            dlnkap_dlnT = dlnkap_rad_dlnT
            return
         end if
         
         call do_electron_conduction( &
            zbar, logRho, logT, &
            kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
         if (ierr /= 0) return
         
         if (is_bad_num(kap_ec)) then
            write(*,*) 'kap_ec', kap_ec
            stop 'Get_kap_Results'
         end if
         
         kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)
         if (dbg) write(*,1) 'kap_rad', kap_rad
         if (dbg) write(*,1) 'kap_ec', kap_ec
         if (dbg) write(*,1) 'kap', kap
         if (dbg) write(*,1) 'log10_cr(kap)', log10_cr(kap)
         
         if (is_bad_num(kap)) then
            ierr = -1; return
            write(*,1) 'kap', kap
            stop 'Get_kap_Results'
         end if
         
         dkap_rad_dlnRho = kap_rad*dlnkap_rad_dlnRho
         dkap_ec_dlnRho = kap_ec*dlnkap_ec_dlnRho
         dkap_dlnRho = (kap_rad*kap_rad*dkap_ec_dlnRho + kap_ec*kap_ec*dkap_rad_dlnRho)/&
            ((kap_rad+kap_ec)*(kap_rad+kap_ec))
         dlnkap_dlnRho = dkap_dlnRho/kap

         if (is_bad_num(dlnkap_dlnRho)) then
            ierr = -1; return
            write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,1) 'kap', kap
            write(*,1) 'dkap_dlnRho', dkap_dlnRho
            write(*,1) 'dkap_ec_dlnRho', dkap_ec_dlnRho
            write(*,1) 'dkap_rad_dlnRho', dkap_rad_dlnRho
            write(*,1) 'kap_rad', kap_rad
            write(*,1) 'kap_ec', kap_ec
            stop 'Get_kap_Results'
         end if
         
         dkap_rad_dlnT = kap_rad*dlnkap_rad_dlnT
         dkap_ec_dlnT = kap_ec*dlnkap_ec_dlnT
         dkap_dlnT = (kap_rad*kap_rad*dkap_ec_dlnT + kap_ec*kap_ec*dkap_rad_dlnT)/&
            ((kap_rad+kap_ec)*(kap_rad+kap_ec))
         dlnkap_dlnT = dkap_dlnT/kap
         
         if (is_bad_num(dlnkap_dlnT)) then
            ierr = -1; return
            write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
            stop 'Get_kap_Results'
         end if
      
      
      end subroutine combine_rad_with_compton_and_conduction
      

      subroutine Compton_Opacity( &
            Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use eos_lib
         use eos_def
         use const_def

         ! same as Thompson at low T and low RHO.
         ! high T and high degeneracy both lengthen mean free path and thus decrease opacity
         ! formula for approximating this effect from Buchler & Yueh, 1976, Apj, 210:440-446.

         ! this approximation breaks down for high degeneracy (eta > 4), but
         ! by then conduction should be dominant anyway, so things should be okay.

         real(dp), intent(in) :: Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         
         real(dp) :: T, rho, edensity, efermi, eta, psi, psi2, Gbar_inverse, kT_mc2, free_e, tmp, &
            d_free_e_dlnRho, d_free_e_dlnT, d_edensity_dlnRho, d_edensity_dlnT, &
            edensity_one_third, d_efermi_dlnRho, d_efermi_dlnT, d_kT_mc2_dlnT, &
            d_eta_dlnRho, d_eta_dlnT, &
            d_psi_dlnRho, d_psi_dlnT, d_psi2_dlnRho, d_psi2_dlnT, &
            d_Gbar_inverse_dlnRho, d_Gbar_inverse_dlnT
         real(dp), parameter :: sigma_Thompson = 6.6524616d-25
         real(dp), parameter :: rho_limit = 5d7
         real(dp), parameter :: T_limit = 2d10
         real(dp), parameter :: pow_3o8pi = 0.242430689510993d0 ! pow_cr(3/(8*pi),2d0/3d0)
         
         logical :: hit_rho_limit, hit_T_limit
         
         include 'formats.dek'
         
         ierr = 0
         free_e = exp_cr(lnfree_e)
         d_free_e_dlnRho = free_e*d_lnfree_e_dlnRho
         d_free_e_dlnT = free_e*d_lnfree_e_dlnT
         
         T = T_in
         hit_T_limit = (T > T_limit)
         if (hit_T_limit) T = T_limit
         rho = rho_in
         hit_rho_limit = (rho > rho_limit)
         if (hit_rho_limit) rho = rho_limit

         edensity = free_e*rho/amu
         if (edensity <= 0) then
            ierr = -1; return
         end if
         d_edensity_dlnRho = d_free_e_dlnRho*rho/amu + edensity
         d_edensity_dlnT = d_free_e_dlnT*rho/amu
         
         tmp = planck_h*planck_h/(2*me)*pow_3o8pi
         edensity_one_third = pow_cr(edensity,1d0/3d0)
         efermi = tmp*edensity_one_third*edensity_one_third
         d_efermi_dlnRho = (2d0/3d0)*tmp*d_edensity_dlnRho/edensity_one_third
         d_efermi_dlnT = (2d0/3d0)*tmp*d_edensity_dlnT/edensity_one_third
         
         kT_mc2 = T * (kerg / (me * clight*clight))
         d_kT_mc2_dlnT = kT_mc2

         eta = (efermi - me*clight*clight)/(kerg*T)
         d_eta_dlnRho = d_efermi_dlnRho/(kerg*T)
         d_eta_dlnT = d_efermi_dlnT/(kerg*T) - eta

         if (eta > 4) then
            eta = 4
            d_eta_dlnRho = 0
            d_eta_dlnT = 0
         else if (eta < -10) then  
            eta = -10
            d_eta_dlnRho = 0
            d_eta_dlnT = 0
         end if

         psi = exp_cr(eta*(0.8168 - eta*0.05522)) ! coefficients from B&Y
         d_psi_dlnRho = psi*(0.8168 - 2*0.05522*eta)*d_eta_dlnRho
         d_psi_dlnT = psi*(0.8168 - 2*0.05522*eta)*d_eta_dlnT
         
         psi2 = psi*psi
         d_psi2_dlnRho = 2*psi*d_psi_dlnRho
         d_psi2_dlnT = 2*psi*d_psi_dlnT

         ! formula for Gbar_inverse from B&Y (incorrectly given as Gbar in their paper!)
         Gbar_inverse = 1.129+0.2965*psi-0.005594*psi2 &
            +(11.47+0.3570*psi+0.1078*psi2)*kT_mc2 &
            +(0.1678*psi-3.249-0.04706*psi2)*kT_mc2*kT_mc2
         d_Gbar_inverse_dlnRho = &
            1.129+0.2965*d_psi_dlnRho-0.005594*d_psi2_dlnRho &
            +(11.47+0.3570*d_psi_dlnRho+0.1078*d_psi2_dlnRho)*kT_mc2 &
            +(0.1678*d_psi_dlnRho-3.249-0.04706*d_psi2_dlnRho)*kT_mc2*kT_mc2
         d_Gbar_inverse_dlnT = &
            1.129+0.2965*d_psi_dlnT-0.005594*d_psi2_dlnT &
            +(11.47+0.3570*d_psi_dlnT+0.1078*d_psi2_dlnT)*kT_mc2 &
            +(0.1678*d_psi_dlnT-3.249-0.04706*d_psi2_dlnT)*kT_mc2*kT_mc2 &
            +(11.47+0.3570*psi+0.1078*psi2)*d_kT_mc2_dlnT &
            +(0.1678*psi-3.249-0.04706*psi2)*2*kT_mc2*d_kT_mc2_dlnT
         
         kap = sigma_Thompson*edensity/(rho*Gbar_inverse)
         if (kap <= 0.0) then
            write(*,*) 'lnfree_e', lnfree_e
            write(*,*) 'Rho', Rho
            write(*,*) 'T', T
            write(*,*) 'edensity', edensity
            write(*,*) 'efermi', efermi
            write(*,*) 'sigma_Thompson', sigma_Thompson
            write(*,*) 'Gbar_inverse', Gbar_inverse
            write(*,*) 'eta',eta
            write(*,*) 'psi', psi
            write(*,*) 'bad compton kap', kap
            stop 1
         end if
         
         dlnkap_dlnRho = kap* &
            (d_edensity_dlnRho/edensity - d_Gbar_inverse_dlnRho/Gbar_inverse - 1)
         dlnkap_dlnT = kap* &
            (d_edensity_dlnT/edensity - d_Gbar_inverse_dlnT/Gbar_inverse)
            
         if (hit_rho_limit) dlnkap_dlnRho = 0
         if (hit_T_limit) dlnkap_dlnT = 0

      end subroutine Compton_Opacity

#ifdef offload
      !dir$ end options
#endif

      end module kap_eval
      