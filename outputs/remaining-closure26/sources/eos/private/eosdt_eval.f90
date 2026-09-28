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
!
! ***********************************************************************

      module eosDT_eval
      use eos_def
      use const_def, only: avo, crad, ln10, arg_not_provided, mp, kerg, dp
      use utils_lib, only: is_bad_num, is_bad_real
      use crlibm_lib

      implicit none

#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif      
			
         
      integer, parameter :: i_doing_Rho = 1
      integer, parameter :: i_which_other = 2
      integer, parameter :: i_handle = 3
      integer, parameter :: i_count = 4
      integer, parameter :: i_species = 5
      
      integer, parameter :: eos_lipar = 5

      integer, parameter :: r_other_value = 1
      integer, parameter :: r_Z = 2
      integer, parameter :: r_X = 3
      integer, parameter :: r_abar = 4
      integer, parameter :: r_zbar = 5
      integer, parameter :: r_rho = 6
      integer, parameter :: r_T = 7
      integer, parameter :: r_the_other_log = 8

      integer, parameter :: eos_lrpar = 8
      
      

      logical, parameter :: dbg = .false.
      
      
      
      contains


      subroutine Get_eosDT_Results(rq,  &
               Z_in, X_in, abar, zbar,   &
               species, chem_id, net_iso, xa, &
               arho, alogrho, atemp, alogtemp, &
               basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         use const_def
         use eosDT_load_tables, only: Load_eosDT_Table
         
         ! INPUT
         
         type (EoS_General_Info), pointer :: rq

         real(dp), intent(in) :: Z_in ! the desired Z
         real(dp), intent(in) :: X_in ! the desired X
            
         real(dp), intent(in) :: abar 
            ! mean atomic number (nucleons per nucleus; grams per mole)
         real(dp), intent(in) :: zbar 
            ! mean charge per nucleus
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: arho, alogrho ! the density
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            
         real(dp), intent(in) :: atemp, alogtemp ! the temperature
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided

         logical, intent(in) :: basic_flag ! if true, only return the basic results, skip the partials.

         
         ! OUTPUT    
              
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 
         integer, intent(out) :: ierr
         
         real(dp), dimension(num_eos_basic_results) ::  &
               res1, d_dlnRho_c_T1, d_dlnT_c_Rho1, d_dabar_c_TRho1, d_dzbar_c_TRho1
         real(dp), dimension(num_eos_basic_results) ::  &
               res2, d_dlnRho_c_T2, d_dlnT_c_Rho2
         real(dp) :: X, Z, Rho, logRho, T, logT
         integer :: iregion
         character (len=256) :: message
         
         integer, parameter :: pure_helm = 1
         integer, parameter :: pure_opal_scvh = 2
         integer, parameter :: blend_in_x = 3
         integer, parameter :: blend_in_y = 4
         integer, parameter :: blend_corner_out = 5
         integer, parameter :: blend_corner_in = 6
         
         real(dp), parameter :: logRho1_max = 3.71d0

         real(dp), parameter :: logQ1 = 5.69d0 ! SCVH full off for logQ > this
         real(dp), parameter :: logQ2 = 5.19d0 ! must have logQ2 < logQ1

         real(dp), parameter :: logRho5 = -9.0d0
         real(dp), parameter :: logRho6 = -9.99d0
         real(dp), parameter :: logRho7 = -13.5d0
         real(dp), parameter :: logRho8 = -13.9d0
         
         real(dp), parameter :: logT5   = 3.6d0
         real(dp), parameter :: logT6   = 3.5d0

         real(dp) :: logRho1
         real(dp) :: logT3
         real(dp) :: logT4
         real(dp) :: logRho2
         
         
         real(dp), parameter :: logT_HELM   = 3d0  ! smallest logT for HELM
         real(dp), parameter :: T_HELM   = 1d3  ! smallest T for HELM
         real(dp), parameter :: logQmin = -9.8d0
         real(dp), parameter :: logQmax = logQ1
         
         real(dp) :: alfa, beta, c_dx, c_dy, logRho_lo, logRho_hi,  &
            logT1, logT2, logT7, logT8, logRho3, logRho4
         real(dp) :: logQ, A, B, sin_pi_alfa, dA_dlnT, dA_dlnRho, dB_dlnT, dB_dlnRho
         real(dp) :: d_dx_dlogT, d_dx_dlogRho, d_dy_dlogT, d_dy_dlogRho
         real(dp) :: &
            lnPgas, Pgas, dlnPgas_dlnRho, dlnPgas_dlnT, dPgas_dlnRho, dPgas_dlnT,  &
            Prad, dPrad_dlnT, P, dP_dlnRho, dP_dlnT, dlnP_dlnRho, dlnP_dlnT,  &
            lnE, energy, dlnE_dlnRho, dlnE_dlnT,  &
            lnS, entropy, dlnS_dlnRho, dlnS_dlnT,  &
            d_alfa_dlogT, d_alfa_dlogRho, d_beta_dlogT, d_beta_dlogRho
         real(dp), parameter :: dZ_transition = 0.01d0 ! must be > 0
         real(dp), parameter :: logT_margin = 0.1d0
         real(dp), parameter :: tiny = 1d-20
         
         logical :: debug
         
         include 'formats.dek'         
         
         ierr = 0
         debug = dbg
         
         if (.not. eosDT_is_initialized) then
            call Load_eosDT_Table(ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'Load_eosDT_Table ierr', ierr
               return
            end if
         end if
         
         if (is_bad_num(X_in) .or. is_bad_num(Z_in)) then
            ierr = -1
            return
         end if
         
         X = X_in; Z = Z_in
         if (X < tiny) X = 0
         if (Z < tiny) Z = 0
         
         !..get temp and rho args
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (alogtemp == arg_not_provided) logT = log10_cr(T)
         if (atemp == arg_not_provided) T = exp10_cr(logT)
         
         if (T <= 0) then
            ierr = -1
            return
         end if
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -3; return
         end if
         if (alogrho == arg_not_provided) logRho = log10_cr(Rho)
         if (arho == arg_not_provided) Rho = exp10_cr(logRho)
         
         if (Rho <= 0) then
            ierr = -1
            return
         end if
         
         if (is_bad_num(Rho) .or. is_bad_num(T)) then
            ierr = -1
            return
         end if
         
         logRho1 = rq% logRho1_OPAL_SCVH_limit
         if (logRho1 > logRho1_max) then
            write(*,*) 'sorry: value for logRho1_OPAL_SCVH_limit is too large.  max allowed is ', &
               logRho1_max
            ierr = -1
            return
         end if
         logRho2 = rq% logRho2_OPAL_SCVH_limit
         
         logT1 = rq% logT_all_HELM
         
         if (logT >= logT1) then ! just use HELM
         
            alfa = 1
            beta = 0
         
         else
         
            logT2 = rq% logT_all_OPAL
            logT3 = (logRho1 - logQ1 + 12)/2
            logT4 = (logRho2 - logQ2 + 12)/2
            logT7 = rq% logT_low_all_SCVH
            logT8 = rq% logT_low_all_HELM

            logRho3 = logQ1 + 2*logT7 - 12
            logRho4 = logQ2 + 2*logT7 - 12
         
            ! check validity of Rho's and T's for region boundaries
            if (logRho1 <= logRho2 .or. logRho2 <= logRho3 .or. &
                logRho3 <= logRho4 .or. logRho4 <= logRho5 .or. &
                logRho5 <= logRho6 .or. logRho6 <= logRho7 .or. &
                logRho7 <= logRho8 .or. &
                logT1 <= logT2 .or. logT2 <= logT3 .or. &
                logT3 <= logT4 .or. logT4 <= logT5 .or. &
                logT5 <= logT6 .or. logT6 <= logT7 .or. &
                logT7 <= logT8) then
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,'(a)') 'must have strictly decreasing values for eos logT + logRho region boundaries'
               write(*,*)
               write(*,1) 'logRho1', logRho1
               write(*,1) 'logRho2', logRho2
               write(*,1) 'logRho3', logRho3
               write(*,1) 'logRho4', logRho4
               write(*,1) 'logRho5', logRho5
               write(*,1) 'logRho6', logRho6
               write(*,1) 'logRho7', logRho7
               write(*,1) 'logRho8', logRho8
               write(*,*)
               write(*,1) 'logT1', logT1
               write(*,1) 'logT2', logT2
               write(*,1) 'logT3', logT3
               write(*,1) 'logT4', logT4
               write(*,1) 'logT5', logT5
               write(*,1) 'logT6', logT6
               write(*,1) 'logT7', logT7
               write(*,1) 'logT8', logT8
               write(*,*)
               write(*,1) 'logQ1', logQ1
               write(*,1) 'logQ2', logQ2
               write(*,*)
               write(*,*)
               if (logT3 <= logT4) then
                  write(*,'(a)') 'must have logRho1 > logRho2 + logQ1 - logQ2'
                  write(*,1) 'must have logRho1 > ', logRho2 + logQ1 - logQ2
                  write(*,1) 'logRho1', logRho1
                  write(*,'(a)') 'logRho1_OPAL_SCVH_limit sets logRho1'
                  write(*,'(a)') 'logRho2_OPAL_SCVH_limit sets logRho2'
                  write(*,1) 'max allowed logRho1 is', logRho1_max
               end if
               write(*,*)
               ierr = -1
               return
            end if
         
            call determine_region
         
            call set_alfa_beta
            if (ierr /= 0) return
            
         end if
         
         if (alfa /= 0d0 .and. logT < logT_HELM .and. logQ >= logQ2) then ! clip T
            if (dbg) then
               write(*, *) 'logT < logT_HELM', logT, logT_HELM
            end if
            logT = logT_HELM
            T = T_HELM
         end if
         
         if (beta /= 0D0) then
            ! before calling OPAL_SCVH, check that Z isn't too large. 
            if (Z >= rq% Z_all_HELM) then
               beta = 0d0
            else if (Z > eos_Zs(num_eosDT_Zs)) then
               beta = beta*(1-(Z-eos_Zs(num_eosDT_Zs))/dZ_transition)
            end if
            alfa = 1d0 - beta
         end if
         
         if (beta > 0.9999) then
            beta = 1; alfa = 0
         else if (alfa > 0.9999) then
            alfa = 1; beta = 0
         end if
               
         if (beta /= 0d0) then
            call Get_OPAL_SCVH_Results( &
                     rq, Z, X, abar, zbar, Rho, logRho, T, logT,  &
                     basic_flag, res2, d_dlnRho_c_T2, d_dlnT_c_Rho2, ierr)
            if (ierr /= 0) then               
               ierr = 0
               beta = 0
               alfa = 1
            end if
         end if
      
         if (alfa /= 0D0) then
            call Get_HELM_PC_Results( &
                     rq, Z, X, abar, zbar,  &
                     species, chem_id, net_iso, xa, &
                     Rho, logRho, T, logT, basic_flag,  &
                     res1, d_dlnRho_c_T1, d_dlnT_c_Rho1,  &
                     d_dabar_c_TRho1, d_dzbar_c_TRho1, debug, ierr)
            if (ierr /= 0) then
               return
            end if
         end if
               
         if (beta == 0d0) then
            res = res1
            if (.not. basic_flag) then
               d_dlnRho_c_T = d_dlnRho_c_T1
               d_dlnT_c_Rho = d_dlnT_c_Rho1
               d_dabar_c_TRho = d_dabar_c_TRho1
               d_dzbar_c_TRho = d_dzbar_c_TRho1
            end if
            return
         end if
         
         if (alfa == 0d0) then
            res = res2
            if (.not. basic_flag) then
               d_dlnRho_c_T = d_dlnRho_c_T2
               d_dlnT_c_Rho = d_dlnT_c_Rho2
               d_dabar_c_TRho = 0d0
               d_dzbar_c_TRho = 0d0
            end if
            return
         end if
         
         A = alfa ! (1 - cospi_cr(alfa))/2 ! smooth the transitions
         B = 1 - A
         
         res = A*res1 + B*res2
         
         if (basic_flag) return
         
         ! just use approximate partials for these
         d_dlnRho_c_T = A*d_dlnRho_c_T1 + B*d_dlnRho_c_T2
         d_dlnT_c_Rho = A*d_dlnT_c_Rho1 + B*d_dlnT_c_Rho2
         d_dabar_c_TRho = A*d_dabar_c_TRho1
         d_dzbar_c_TRho = A*d_dzbar_c_TRho1

         
         contains
         
         
#ifdef offload
         !dir$ attributes offload : mic :: determine_region
#endif      
         subroutine determine_region
            logical, parameter :: dbg = .false.
            logQ = logRho - 2*logT + 12
            d_dx_dlogRho=0; d_dx_dlogT=0
            d_dy_dlogT=0; d_dy_dlogRho=0
            if (logT >= logT1 .or. logT <= logT8 .or. logRho >= logRho1) then
               if (dbg .and. logT >= logT1) write(*,*) 'logT >= logT1', logT >= logT1
               if (dbg .and. logT <= logT8) write(*,*) 'logT <= logT8', logT <= logT8
               if (dbg .and. logRho >= logRho1) write(*,*) 'logRho >= logRho1', logRho >= logRho1
               iregion = pure_helm
            else if (logQ < logQmin) then          
               if (dbg) write(*,*) 'logQ < logQmin', logQ < logQmin
               iregion = pure_helm
            else if (logQ > logQmax) then          
               if (dbg) write(*,*) 'logQ > logQmax', logQ > logQmax
               iregion = pure_helm
            else if (logT >= logT2) then         
               if (dbg) write(*,*) 'logT >= logT2', logT >= logT2
               d_dy_dlogT = 1/(logT1 - logT2)
               c_dy = (logT - logT2)*d_dy_dlogT
               if (logRho > logRho2) then
                  d_dx_dlogRho = 1/(logRho1 - logRho2)
                  c_dx = (logRho - logRho2)*d_dx_dlogRho 
                  iregion = blend_corner_out
               else
                  iregion = blend_in_y
               end if        
                
            else if (logT >= logT3) then  ! NOTE: this assumes logT3 > logT4
               if (dbg) write(*,*) 'logT >= logT3', logT >= logT3
               if (logRho > logRho2) then
                  d_dx_dlogRho = 1/(logRho1 - logRho2)
                  c_dx = (logRho - logRho2)*d_dx_dlogRho
                  iregion = blend_in_x
               else
                  iregion = pure_opal_scvh
               end if    
               
               
            else if (logT >= logT4) then         
               if (dbg) write(*,*) 'logT >= logT4', logT >= logT4
               logRho_hi = logQ1 + 2*logT - 12
               if (logRho >= logRho_hi) then
                  iregion = pure_helm
               else if (logRho > logRho2) then
                  d_dx_dlogRho = 1/(logRho_hi - logRho2)
                  c_dx = (logRho - logRho2)*d_dx_dlogRho
                  iregion = blend_in_x
               else ! logRho <= logRho2
                  iregion = pure_opal_scvh
               end if                   
            else if (logRho > logRho4) then         
               if (dbg) write(*,*) 'logRho > logRho4', logRho > logRho4
               if (logT > logT7) then
                  A = ((logQ1+2*logT4-12) - logRho3)/(logT4-logT7)
                  logRho_hi = logRho3 + (logT-logT7)*A
                  B = (logRho2-logRho4)/(logT4-logT7)
                  logRho_lo = logRho4 + (logT-logT7)*B
                  if (logRho >= logRho_hi) then
                     iregion = pure_helm
                  else if (logRho >= logRho_lo) then
                     c_dx = (logRho - logRho_lo)/(logRho_hi - logRho_lo)
                     d_dx_dlogRho = 1/(logRho3 - logRho4 + (A - B)*(logT - logT7))
                     d_dx_dlogT = (B*(logRho - logRho3) - &
                        A*(logRho - logRho4))*d_dx_dlogRho*d_dx_dlogRho
                     iregion = blend_in_x               
                  else
                     iregion = pure_opal_scvh               
                  end if            
               else ! logT is > logT8            
               if (dbg) write(*,*) 'logT > logT8', logT > logT8
                  if (logRho > logRho3) then
                     iregion = pure_helm
                  else ! logRho is > logRho4
                     d_dx_dlogRho = 1/(logRho3 - logRho4)
                     c_dx = (logRho - logRho4)*d_dx_dlogRho
                     d_dy_dlogT = 1/(logT8 - logT7)
                     c_dy = (logT - logT7)*d_dy_dlogT
                     iregion = blend_corner_out
                  end if            
               end if         
            else ! logRho <= logRho4         
               if (dbg) write(*,*) 'logRho <= logRho4', logRho <= logRho4
               if (logT >= logT5) then
                  if (logRho > logRho7) then
                     iregion = pure_opal_scvh
                  else if (logRho > logRho8) then
                     d_dx_dlogRho = 1/(logRho8 - logRho7)
                     c_dx = (logRho - logRho7)*d_dx_dlogRho
                     iregion = blend_in_x
                  else
                     iregion = pure_helm
                  end if
               else if (logT >= logT6) then            
                  if (logRho >= logRho5) then
                     iregion = pure_opal_scvh
                  else 
                     if (logRho > logRho6) then
                        d_dx_dlogRho = 1/(logRho5 - logRho6)
                        c_dx = (logRho - logRho6)*d_dx_dlogRho
                        d_dy_dlogT = 1/(logT5 - logT6)
                        c_dy = (logT - logT6)*d_dy_dlogT
                        iregion = blend_corner_in
                     else if (logRho > logRho7) then
                        d_dy_dlogT = 1/(logT6 - logT5)
                        c_dy = (logT - logT5)*d_dy_dlogT
                        iregion = blend_in_y
                     else if (logRho > logRho8) then
                        d_dx_dlogRho = 1/(logRho8 - logRho7)
                        c_dx = (logRho - logRho7)*d_dx_dlogRho
                        d_dy_dlogT = 1/(logT6 - logT5)
                        c_dy = (logT - logT5)*d_dy_dlogT
                        iregion = blend_corner_out
                     else
                        iregion = pure_helm
                     end if               
                  end if               
               else if (logT >= logT7) then               
                  if (logRho >= logRho5) then
                     iregion = pure_opal_scvh
                  else if (logRho > logRho6) then
                     d_dx_dlogRho = 1/(logRho6 - logRho5)
                     c_dx = (logRho - logRho5)*d_dx_dlogRho
                     iregion = blend_in_x
                  else ! logRho <= logRho6
                     iregion = pure_helm
                  end if            
               else ! logT > logT8            
                  if (logRho >= logRho5) then
                     d_dy_dlogT = 1/(logT8 - logT7)
                     c_dy = (logT - logT7)*d_dy_dlogT
                     iregion = blend_in_y
                  else if (logRho >= logRho6) then
                     d_dx_dlogRho = 1/(logRho6 - logRho5)
                     c_dx = (logRho - logRho5)*d_dx_dlogRho
                     d_dy_dlogT = 1/(logT8 - logT7)
                     c_dy = (logT - logT7)*d_dy_dlogT
                     iregion = blend_corner_out
                  else ! logRho < logRho6
                     iregion = pure_helm
                  end if            
               end if         
            end if
            
         end subroutine determine_region


#ifdef offload
         !dir$ attributes offload : mic :: set_alfa_beta
#endif      
         subroutine set_alfa_beta
            logical, parameter :: dbg = .false.
            include 'formats.dek'
            
            d_alfa_dlogT=0; d_alfa_dlogRho=0; d_beta_dlogT=0; d_beta_dlogRho=0
            if (iregion == pure_helm) then
               if (dbg) write(*,*) 'iregion == pure_helm'
               alfa = 1
               beta = 0
            else if (iregion == pure_opal_scvh) then
               if (dbg) write(*,*) 'iregion == pure_opal_scvh'
               alfa = 0
               beta = 1
            else if (iregion == blend_in_y) then
               if (dbg) write(*,*) 'iregion == blend_in_y'
               alfa = c_dy
               d_alfa_dlogT = d_dy_dlogT
               d_alfa_dlogRho = d_dy_dlogRho
               beta = 1 - alfa
               d_beta_dlogT = -d_alfa_dlogT
               d_beta_dlogRho = -d_alfa_dlogRho
            else if (iregion == blend_in_x) then
               if (dbg) write(*,*) 'iregion == blend_in_x'
               alfa = c_dx
               d_alfa_dlogT = d_dx_dlogT
               d_alfa_dlogRho = d_dx_dlogRho
               beta = 1 - alfa
               d_beta_dlogT = -d_alfa_dlogT
               d_beta_dlogRho = -d_alfa_dlogRho
            else if (iregion == blend_corner_out) then
               if (dbg) write(*,*) 'iregion == blend_corner_out'
               alfa = sqrt(c_dx*c_dx + c_dy*c_dy)
               if (alfa >= 1d0) then
                  alfa = 1
                  d_alfa_dlogT = 0
                  d_alfa_dlogRho = 0
               else if (alfa < 1d-10) then
                  alfa = 0
                  d_alfa_dlogT = 0
                  d_alfa_dlogRho = 0
               else
                  d_alfa_dlogT = (c_dx*d_dx_dlogT + c_dy*d_dy_dlogT)/alfa
                  d_alfa_dlogRho = (c_dx*d_dx_dlogRho + c_dy*d_dy_dlogRho)/alfa
               end if
               beta = 1 - alfa
               d_beta_dlogT = -d_alfa_dlogT
               d_beta_dlogRho = -d_alfa_dlogRho
            else if (iregion == blend_corner_in) then
               if (dbg) write(*,*) 'iregion == blend_corner_in'
               beta = sqrt(c_dx*c_dx + c_dy*c_dy)
               if (beta >= 1d0) then
                  beta = 1
                  d_beta_dlogT = 0
                  d_beta_dlogRho = 0
               else if (beta < 1d-10) then
                  beta = 0
                  d_beta_dlogT = 0
                  d_beta_dlogRho = 0
               else
                  d_beta_dlogT = (c_dx*d_dx_dlogT + c_dy*d_dy_dlogT)/beta
                  d_beta_dlogRho = (c_dx*d_dx_dlogRho + c_dy*d_dy_dlogRho)/beta
               end if
               alfa = 1 - beta
               d_alfa_dlogT = -d_beta_dlogT
               d_alfa_dlogRho = -d_beta_dlogRho
            else
               ierr = -1
               return
            end if
         end subroutine set_alfa_beta


      end subroutine Get_eosDT_Results
      

      subroutine Get_HELM_Results( &
               X, abar, zbar, arho, alogrho, atemp, alogtemp, &
               include_radiation, always_skip_elec_pos, res, ierr)
         use const_def
         use helm

         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: X, abar, zbar
         real(dp), intent(in) :: arho, alogrho
         real(dp), intent(in) :: atemp, alogtemp 
         logical, intent(in) :: include_radiation, always_skip_elec_pos
         real(dp), intent(out) :: res(:) ! (num_helm_results)
         integer, intent(out) :: ierr ! 0 means AOK.

         real(dp) :: Rho, logRho, T, logT, dse, dpe, dsp
         
         logical, parameter :: clip_to_table_boundaries = .true.
         
         include 'formats.dek'
         
         ierr = 0

         !..get temp and rho args
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (atemp == arg_not_provided) T = exp10_cr(logT)
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -3; return
         end if
         if (arho == arg_not_provided) Rho = exp10_cr(logRho)
         
         call helmeos2(T, logT, Rho, logRho, X, abar, zbar, res,  &
                  clip_to_table_boundaries, include_radiation, &
                  always_skip_elec_pos, ierr)
         res(h_valid) = 1

         
      end subroutine Get_HELM_Results

      
      subroutine Get_OPAL_SCVH_Results( &
               rq, Z, X, abar, zbar, Rho, logRho, T, logT, basic_flag,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         use chem_def
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar, Rho, logRho, T, logT
         logical, intent(in) :: basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         integer, intent(out) :: ierr

         real(dp), dimension(num_eos_basic_results, num_eosDT_Zs) ::  &
            res_zx, d_dlnRho_c_T_zx, d_dlnT_c_Rho_zx, d_res_dX
         real(dp), dimension(num_eos_basic_results) :: d_dX, d_dZ
         real(dp) :: d_res_dZ(num_eos_basic_results), dZ, denom, c(3), &
            dP_dZ, dlnP_dZ
         real(dp), parameter :: tiny = 1d-6

         character (len=256) :: message
         integer :: iz, j, ci
         logical, parameter :: OPAL_SCVH_dbg = dbg
         
         include 'formats.dek'

         ierr = 0
         
         if (num_eosDT_Zs < 3) then
            write(*, *) 'error: Get_OPAL_SCVH_Results assumes num_eosDT_Zs >= 3'
            stop 1
         end if
         
         if (eos_Zs(1) /= 0) then
            write(*, *) 'error: Get_OPAL_SCVH_Results assumes eos_Zs(1) == 0'
            stop 1
         end if
         
         if (abs(eos_Zs(1) - 2*eos_Zs(2) + eos_Zs(3)) > tiny) then
            write(*, *) 'error: Get_OPAL_SCVH_Results assumes equal spaced Zs(1:3)'
            stop 1
         end if
         
         if (Z < 1d-20) then
            !write(*,1) 'Get_OPAL_SCVH_for_X Z', Z
            call Get_OPAL_SCVH_for_X( &
                  rq, 1, X, abar, zbar, Rho, logRho, T, logT,  &
                  basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
                  ierr)
            return
         end if
         
         if (Z > eos_Zs(3) .and. num_eosDT_Zs == 5) then
         
            if (Z <= eos_Zs(4)) then
               call do_interp2(3,4,ierr)
               if (ierr /= 0) return
            else if (Z <= eos_Zs(5) - tiny) then
               call do_interp2(4,5,ierr)
               if (ierr /= 0) return
            else
               call Get_OPAL_SCVH_for_X( &
                  rq, 5, X, abar, zbar, Rho, logRho, T, logT,  &
                  basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
                  ierr)
               return
            end if
            
         else
         
            do iz = 1, 3
               if (.false.) write(*,2) 'Get_OPAL_SCVH_for_X iz', iz, eos_Zs(iz)
               call Get_OPAL_SCVH_for_X( &
                  rq, iz, X, abar, zbar, Rho, logRho, T, logT,  &
                  basic_flag, res_zx(:, iz), d_dlnRho_c_T_zx(:, iz), d_dlnT_c_Rho_zx(:, iz),  &
                  ierr)
               if (ierr /= 0) return
            end do

            dZ = eos_Zs(2) - eos_Zs(1)
            denom = 2*dZ*dZ
            c(1) = (2*dZ*dZ - 3*dZ*Z + Z*Z)/denom
            c(2) = 2*(2*dZ-Z)*Z/denom
            c(3) = Z*(Z-dZ)/denom         

         end if
         
         res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)
         if (basic_flag) return
         
         d_dlnRho_c_T(:) =  &
           c(1)*d_dlnRho_c_T_zx(:, 1) +  &
           c(2)*d_dlnRho_c_T_zx(:, 2) +  &
           c(3)*d_dlnRho_c_T_zx(:, 3)
     
         d_dlnT_c_Rho(:) = &
           c(1)*d_dlnT_c_Rho_zx(:, 1) +  &
           c(2)*d_dlnT_c_Rho_zx(:, 2) +  &
           c(3)*d_dlnT_c_Rho_zx(:, 3)
     
     
         contains
         
         
#ifdef offload
         !dir$ attributes offload : mic :: do_interp2
#endif      
         subroutine do_interp2(iz1, iz2, ierr)
            use const_def, only: pi
            integer, intent(in) :: iz1, iz2
            integer, intent(out) :: ierr
            real*8 :: Z1, Z2, alfa, beta
            ierr = 0
            Z1 = eos_Zs(iz1)
            Z2 = eos_Zs(iz2)
            alfa = (Z - Z1) / (Z2 - Z1)
            !alfa = (1 - cospi_cr(alfa))/2 ! smooth the transitions
            beta = 1 - alfa
            call Get_OPAL_SCVH_for_X( &
                  rq, iz1, X, abar, zbar, Rho, logRho, T, logT,  &
                  basic_flag, res_zx(:, 1), d_dlnRho_c_T_zx(:, 1), d_dlnT_c_Rho_zx(:, 1),  &
                  ierr)
            if (ierr /= 0) return
            call Get_OPAL_SCVH_for_X( &
                  rq, iz2, X, abar, zbar, Rho, logRho, T, logT,  &
                  basic_flag, res_zx(:, 2), d_dlnRho_c_T_zx(:, 2), d_dlnT_c_Rho_zx(:, 2),  &
                  ierr)
            if (ierr /= 0) return
            c(1) = beta
            c(2) = alfa
            c(3) = 0
            res_zx(:,3) = 0
            d_dlnRho_c_T_zx(:,3) = 0
            d_dlnT_c_Rho_zx(:,3) = 0
         end subroutine do_interp2
     
      end subroutine Get_OPAL_SCVH_Results

      
      subroutine Get_OPAL_SCVH_for_X( &
               rq, iz, X, abar, zbar, Rho, logRho, T, logT,  &
               basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: iz ! the index in eos_Zs
         real(dp), intent(in) :: X, abar, zbar, Rho, logRho, T, logT
         logical, intent(in) :: basic_flag
         real(dp), dimension(:, :), pointer :: d_dx ! (num_eos_basic_results, species)
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         integer, intent(out) :: ierr

         real(dp), dimension(num_eos_basic_results, 4) ::  &
               res_zx, d_dlnRho_c_T_zx, d_dlnT_c_Rho_zx
         real(dp) :: dX, c(4), denom, delX, coef, dmu_dX, dlnmu_dX, dlnP_dX
         character (len=256) :: message
         integer :: ix, ix_lo, ix_hi, j, num_Xs
         real(dp), parameter :: tiny = 1d-6
         logical, parameter :: dbg_for_X = dbg ! .or. .true.
         
         include 'formats.dek'
         
         ierr = 0
         
         if (num_eos_Xs /= 6) then
            write(*, *) 'error: Get_OPAL_SCVH_for_X assumes num_eos_Xs == 6'
            stop 1
         end if
         
         if (eos_Xs(1) /= 0) then
            write(*, *) 'error: Get_OPAL_SCVH_for_X assumes eos_Xs(1) == 0'
            stop 1
         end if
         
         num_Xs = num_eos_Xs_for_Z(iz)
         
         if (X < tiny .or. num_Xs == 1) then
            !write(*,1) 'call Get_eosDT_XTable_Results X', X
            call Get_eosDT_XTable_Results(rq, 1, iz, Rho, logRho, T, logT,  &
               basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
            return
         end if
         
         dX = eos_Xs(2)-eos_Xs(1)
         
         do ix = 3, num_Xs
            if (abs(dX - (eos_Xs(ix) - eos_Xs(ix-1))) > tiny) then
               write(*, *) 'error: Get_OPAL_SCVH_for_X assumes equal spaced Xs'
               stop 1
            end if
         end do
         
         ix_hi = -1
         if (X <= eos_Xs(2)) then
            ix_lo = 1; ix_hi = 3
         else if (X >= eos_Xs(num_Xs-1)) then
            ix_lo = num_Xs-2; ix_hi = num_Xs
         else
            do ix = 3, num_Xs-1
               if (X <= eos_Xs(ix)) then
                  ix_lo = ix-2; ix_hi = ix+1; exit
               end if
            end do
         end if
         
         if (ix_hi < 0) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
            write(*, *) 'error: Get_OPAL_SCVH_for_X logic bug'
            stop 1
         end if
         
         if (dbg_for_X) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
         end if
         
         do ix=ix_lo, ix_hi
            j = ix-ix_lo+1
            
            if (.false.) write(*,2) 'Get_eosDT_XTable_Results X', ix, eos_Xs(ix)
            
            call Get_eosDT_XTable_Results(rq, ix, iz, Rho, logRho, T, logT,  &
               basic_flag, res_zx(:, j), d_dlnRho_c_T_zx(:, j), d_dlnT_c_Rho_zx(:, j), ierr)
            if (ierr /= 0) return
         end do
         
         delX = X - eos_Xs(ix_lo)
         if (ix_hi-ix_lo==2) then
         
            denom = 2*dX*dX
            c(1) = (2*dX*dX - 3*dX*delX + delX*delX)/denom
            c(2) = 2*(2*dX-delX)*delX/denom
            c(3) = delX*(delX-dX)/denom
            res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)

            if (dbg_for_X) then
               write(*, 1) 'res(i_grad_ad)', res(i_grad_ad)
               write(*, 1) 'res_zx(i_grad_ad, 1)', res_zx(i_grad_ad, 1)
               write(*, 1) 'res_zx(i_grad_ad, 2)', res_zx(i_grad_ad, 2)
               write(*, 1) 'res_zx(i_grad_ad, 3)', res_zx(i_grad_ad, 3)
               write(*, 1) 'c(1)*res_zx(i_grad_ad, 1)', c(1)*res_zx(i_grad_ad, 1)
               write(*, 1) 'c(2)*res_zx(i_grad_ad, 2)', c(2)*res_zx(i_grad_ad, 2)
               write(*, 1) 'c(3)*res_zx(i_grad_ad, 3)', c(3)*res_zx(i_grad_ad, 3)
            end if

            if (basic_flag) return
            
            d_dlnRho_c_T(:) =  &
               c(1)*d_dlnRho_c_T_zx(:,1) +  &
               c(2)*d_dlnRho_c_T_zx(:,2) +  &
               c(3)*d_dlnRho_c_T_zx(:,3)
            d_dlnT_c_Rho(:) = &
               c(1)*d_dlnT_c_Rho_zx(:,1) +  &
               c(2)*d_dlnT_c_Rho_zx(:,2) +  &
               c(3)*d_dlnT_c_Rho_zx(:,3)
     
            c(1) = (2*delX-3*dX)/denom
            c(2) = 4*(dX-delX)/denom
            c(3) = (2*delX-dX)/denom
            
            dlnP_dX = c(1)*res_zx(i_lnPgas, 1)  &
                  + c(2)*res_zx(i_lnPgas, 2) + c(3)*res_zx(i_lnPgas, 3)
            dmu_dX = c(1)*res_zx(i_mu, 1)  &
                  + c(2)*res_zx(i_mu, 2) + c(3)*res_zx(i_mu, 3)
            
         else
         
            coef = (X-eos_Xs(ix_lo+1))/dX 
            ! coef = fractional location of X between 2nd and 3rd X's for fit.
            ! coef is weight for the quadratic based on points 2, 3, 4 of fit.
            ! (1-coef) is weight for quadratic based on points 1, 2, 3 of fit.
            if (coef < 0 .or. coef > 1) then
               ierr = -1
               return
            end if
            c(1) = -coef*(coef-1)*(coef-1)/2
            c(2) = (2 - coef*coef*(5 - 3*coef))/2
            c(3) = coef*(1 + coef*(4 - 3*coef))/2
            c(4) = coef*coef*(coef-1)/2
            res(:) = c(1)*res_zx(:, 1) + &
                        (c(2)*res_zx(:, 2) + &
                           (c(3)*res_zx(:, 3) + &
                              c(4)*res_zx(:, 4)))
     
            if (basic_flag) return
            
            d_dlnRho_c_T(:) =  &
               c(1)*d_dlnRho_c_T_zx(:, 1) +  &
                  (c(2)*d_dlnRho_c_T_zx(:, 2) +  &
                     (c(3)*d_dlnRho_c_T_zx(:, 3) +  &
                           c(4)*d_dlnRho_c_T_zx(:, 4)))
            d_dlnT_c_Rho(:) = &
               c(1)*d_dlnT_c_Rho_zx(:, 1) +  &
                  (c(2)*d_dlnT_c_Rho_zx(:, 2) +  &
                     (c(3)*d_dlnT_c_Rho_zx(:, 3) +  &
                           c(4)*d_dlnT_c_Rho_zx(:, 4)))
     
            denom = 2*dX
            c(1) = (-1 + coef*(3 - 2*coef))/denom
            c(2) = coef*(-7 + 6*coef)/denom
            c(3) = (1 + coef*(5 - 6*coef))/denom
            c(4) = coef*(-1 + 2*coef)/denom
            
            dlnP_dX = &
               c(1)*res_zx(i_lnPgas, 1) + &
                  (c(2)*res_zx(i_lnPgas, 2) + &
                     (c(3)*res_zx(i_lnPgas, 3) + &
                        c(4)*res_zx(i_lnPgas, 4)))
            dmu_dX = &
               c(1)*res_zx(i_mu, 1) + &
                  (c(2)*res_zx(i_mu, 2) + &
                     (c(3)*res_zx(i_mu, 3) + &
                        c(4)*res_zx(i_mu, 4)))
     
         end if
         
         
      end subroutine Get_OPAL_SCVH_for_X
         
      
      subroutine Get_eosDT_XTable_Results(rq, ix, iz, Rho, logRho_in, T, logT_in,  &
               basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: ix, iz
         real(dp), intent(in) :: Rho, logRho_in, T, logT_in
         logical, intent(in) :: basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
            ! only does partials wrt lnRho and lnT
         integer, intent(out) :: ierr
         
         integer, parameter :: nv = num_eos_vals
         real :: fval(nv), df_dx(nv), df_dy(nv), df_dlogRho_c_T(nv), df_dlogT_c_Rho(nv)
         real :: logQ0, logQ1, logT0, logT1
         integer :: iQ, jtemp
         type (EosDT_XZ_Info), pointer :: ep
         
         real :: logRho, logT, logQ
         
         include 'formats.dek'

         logRho = logRho_in
         logT = logT_in
         logQ = logRho - 2*logT + 12

         ierr = 0
         ep => eosDT_XZ_data(ix,iz)     
         
         call Locate_logQ(rq, ep, logQ, iQ, logQ0, logQ1, ierr)
         if (ierr /= 0) then
            write(*,1) 'failed in Locate_logQ', logQ
            return
         end if
      
         call Locate_logT(rq, ep, logT, jtemp, logT0, logT1, ierr)
         if (ierr /= 0) then
            write(*,1) 'failed in Locate_logT', logT
            return
         end if
         
         call Do_EoS_Interpolations( &
                  nv, ep% num_logQs, ep% logQs, ep% num_logTs, ep% logTs,  &
                  ep% tbl1, iQ, jtemp, logQ0, logQ, logQ1, logT0, logT, logT1,  &
                  fval, df_dx, df_dy, ierr)
         if (ierr /= 0) then
            write(*,1) 'failed in Do_EoS_Interpolations'
            return
         end if
     
         res(i_lnPgas) = fval(eosDT_ilnPgas)
         res(i_lnE) = fval(eosDT_ilnE)
         res(i_lnS) = fval(eos_ilnS)
         
         res(i_grad_ad) = fval(eos_igrad_ad)
         res(i_chiRho) = fval(eos_ichiRho)
         res(i_chiT) = fval(eos_ichiT)
         
         res(i_Cp) = fval(eos_iCp)
         res(i_Cv) = fval(eos_iCv)
         
         res(i_dE_dRho) = fval(eos_idE_dRho)
         res(i_dS_dT) = fval(eos_idS_dT)
         res(i_dS_dRho) = fval(eos_idS_dRho)
         
         res(i_mu) = fval(eos_imu)
         res(i_lnfree_e) = fval(eos_ilnfree_e)
         res(i_gamma1) = fval(eos_igamma1)
         res(i_gamma3) = fval(eos_igamma3)
         res(i_eta) = fval(eos_ieta)
         
         if (dbg) write(*, 4) 'res(i_eta)', ix, iQ, jtemp, res(i_eta),  &
               logRho_in, logT_in, eos_Xs(ix), eos_Zs(iz), res(i_gamma1)
         
         if (basic_flag) return
         
         ! convert df_dx and df_dy to df_dlogRho_c_T and df_dlogT_c_Rho
         
         ! df_dx is df_dlogQ at const T
         ! df_dy is df_dlogT_c_Rho at const Q
         ! logQ = logRho - 2*logT + 12
         
         ! f = f(logQ(logRho,logT),logT)
         ! df/dlogRho|T = df/dlogQ|T * dlogQ/dlogRho|T = df_dx
         ! df/dlogT|Rho = df/dlogT|Q + df/dlogQ|T * dlogQ/dlogT|Rho = df_dy - 2*df_dx
         
         df_dlogRho_c_T = df_dx
         df_dlogT_c_Rho = df_dy - 2*df_dx

         d_dlnRho_c_T(i_lnPgas) = df_dlogRho_c_T(eosDT_ilnPgas)/ln10
         d_dlnRho_c_T(i_lnE) = df_dlogRho_c_T(eosDT_ilnE)/ln10
         d_dlnRho_c_T(i_lnS) = df_dlogRho_c_T(eos_ilnS)/ln10
         d_dlnRho_c_T(i_grad_ad) = df_dlogRho_c_T(eos_igrad_ad)/ln10
         d_dlnRho_c_T(i_chiRho) = df_dlogRho_c_T(eos_ichiRho)/ln10
         d_dlnRho_c_T(i_chiT) = df_dlogRho_c_T(eos_ichiT)/ln10
         
         d_dlnRho_c_T(i_Cp) = df_dlogRho_c_T(eos_iCp)/ln10
         d_dlnRho_c_T(i_Cv) = df_dlogRho_c_T(eos_iCv)/ln10
         d_dlnRho_c_T(i_dE_dRho) = df_dlogRho_c_T(eos_idE_dRho)/ln10
         d_dlnRho_c_T(i_dS_dT) = df_dlogRho_c_T(eos_idS_dT)/ln10
         d_dlnRho_c_T(i_dS_dRho) = df_dlogRho_c_T(eos_idS_dRho)/ln10
         d_dlnRho_c_T(i_mu) = df_dlogRho_c_T(eos_imu)/ln10
         d_dlnRho_c_T(i_lnfree_e) = df_dlogRho_c_T(eos_ilnfree_e)/ln10
         d_dlnRho_c_T(i_gamma1) = df_dlogRho_c_T(eos_igamma1)/ln10
         d_dlnRho_c_T(i_gamma3) = df_dlogRho_c_T(eos_igamma3)/ln10
         d_dlnRho_c_T(i_eta) = df_dlogRho_c_T(eos_ieta)/ln10
      
         d_dlnT_c_Rho(i_lnPgas) = df_dlogT_c_Rho(eosDT_ilnPgas)/ln10
         d_dlnT_c_Rho(i_lnE) = df_dlogT_c_Rho(eosDT_ilnE)/ln10
         d_dlnT_c_Rho(i_lnS) = df_dlogT_c_Rho(eos_ilnS)/ln10
         d_dlnT_c_Rho(i_grad_ad) = df_dlogT_c_Rho(eos_igrad_ad)/ln10
         d_dlnT_c_Rho(i_chiRho) = df_dlogT_c_Rho(eos_ichiRho)/ln10
         d_dlnT_c_Rho(i_chiT) = df_dlogT_c_Rho(eos_ichiT)/ln10
         d_dlnT_c_Rho(i_Cp) = df_dlogT_c_Rho(eos_iCp)/ln10
         d_dlnT_c_Rho(i_Cv) = df_dlogT_c_Rho(eos_iCv)/ln10
         d_dlnT_c_Rho(i_dE_dRho) = df_dlogT_c_Rho(eos_idE_dRho)/ln10
         d_dlnT_c_Rho(i_dS_dT) = df_dlogT_c_Rho(eos_idS_dT)/ln10
         d_dlnT_c_Rho(i_dS_dRho) = df_dlogT_c_Rho(eos_idS_dRho)/ln10
         d_dlnT_c_Rho(i_mu) = df_dlogT_c_Rho(eos_imu)/ln10
         d_dlnT_c_Rho(i_lnfree_e) = df_dlogT_c_Rho(eos_ilnfree_e)/ln10
         d_dlnT_c_Rho(i_gamma1) = df_dlogT_c_Rho(eos_igamma1)/ln10
         d_dlnT_c_Rho(i_gamma3) = df_dlogT_c_Rho(eos_igamma3)/ln10
         d_dlnT_c_Rho(i_eta) = df_dlogT_c_Rho(eos_ieta)/ln10

      end subroutine Get_eosDT_XTable_Results
      
      
      subroutine Locate_logQ(rq, ep, logQ, iQ, logQ0, logQ1, ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosDT_XZ_Info), pointer :: ep
         real, intent(inout) :: logQ
         integer, intent(out) :: iQ
         real, intent(out) :: logQ0, logQ1
         integer, intent(out) :: ierr
      
         ierr = 0
         iQ = int((logQ - ep% logQ_min) / ep% del_logQ + 1d-4) + 1
         
         if (iQ < 1 .or. iQ >= ep% num_logQs) then
            
            if (iQ < 1) then
               iQ = 1
               logQ0 = ep% logQ_min
               logQ1 = logQ0 + ep% del_logQ
               logQ = logQ0
            else
               iQ = ep% num_logQs-1
               logQ0 = ep% logQ_min + (iQ-1) * ep% del_logQ
               logQ1 = logQ0 + ep% del_logQ
               logQ = logQ1
            end if
            
         else
         
            logQ0 = ep% logQ_min + (iQ-1) * ep% del_logQ
            logQ1 = logQ0 + ep% del_logQ

         end if

      end subroutine Locate_logQ
      
      
      subroutine Locate_logT(rq, ep, logT, iT, logT0, logT1, ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosDT_XZ_Info), pointer :: ep
         real, intent(inout) :: logT
         integer, intent(out) :: iT
         real, intent(out) :: logT0, logT1
         integer, intent(out) :: ierr
      
         ierr = 0
         iT = int((logT - ep% logT_min) / ep% del_logT + 1d-4) + 1
         
         if (iT < 1 .or. iT >= ep% num_logTs) then
            
            if (iT < 1) then
               iT = 1
               logT0 = ep% logT_min
               logT1 = logT0 + ep% del_logT
               logT = logT0
            else
               iT = ep% num_logTs-1
               logT0 = ep% logT_min + (iT-1) * ep% del_logT
               logT1 = logT0 + ep% del_logT
               logT = logT1
            end if
            
         else
         
            logT0 = ep% logT_min + (iT-1) * ep% del_logT
            logT1 = logT0 + ep% del_logT

         end if

      end subroutine Locate_logT
      
      
      subroutine Do_EoS_Interpolations( &
                   nv, nx, x, ny, y, fin1, i, j,  &
                   x0, xget, x1, y0, yget, y1,  &
                   fval, df_dx, df_dy, ierr)
         use interp_2d_lib_sg, only: interp_evbipm_sg
         integer, intent(in) :: nv, nx, ny
         real, intent(in) :: x(:) ! (nx)
         real, intent(in) :: y(:) ! (ny)
         real, intent(in), pointer :: fin1(:) ! =(4,nv,nx,ny)
         integer, intent(in) :: i, j           ! target cell in f
         real, intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
         real, intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
         real, intent(out) :: fval(:), df_dx(:), df_dy(:) ! (nv)
         integer, intent(out) :: ierr
   
         real, parameter :: sixth = 1.0/6.0
         real, parameter :: z36th = 1.0/36.0
         real :: xp, xpi, xp2, xpi2, ax, axbar, bx, bxbar, cx, cxi, hx2, cxd, cxdi, hx, hxi
         real :: yp, ypi, yp2, ypi2, ay, aybar, by, bybar, cy, cyi, hy2, cyd, cydi, hy, hyi
         real :: sixth_hx2, sixth_hy2, z36th_hx2_hy2
         real :: sixth_hx, sixth_hxi_hy2, z36th_hx_hy2
         real :: sixth_hx2_hyi, sixth_hy, z36th_hx2_hy
         integer :: k, ip1, jp1
         real, pointer :: fin(:,:,:,:)
         
         include 'formats.dek'
         
         ierr = 0
         
         fin(1:4,1:nv,1:nx,1:ny) => fin1(1:4*nv*nx*ny)
         
         hx=x1-x0
         hxi=1.0/hx
         hx2=hx*hx
   
         xp=(xget-x0)*hxi
         xpi=1.0-xp
         xp2=xp*xp
         xpi2=xpi*xpi

         ax=xp2*(3.0-2.0*xp)
         axbar=1.0-ax
         
         bx=-xp2*xpi
         bxbar=xpi2*xp
   
         cx=xp*(xp2-1.0)
         cxi=xpi*(xpi2-1.0)
         cxd=3.0*xp2-1.0
         cxdi=-3.0*xpi2+1.0
   
         hy=y1-y0
         hyi=1.0/hy
         hy2=hy*hy
   
         yp=(yget-y0)*hyi
         ypi=1.0-yp
         yp2=yp*yp
         ypi2=ypi*ypi

         ay=yp2*(3.0-2.0*yp)
         aybar=1.0-ay
         
         by=-yp2*ypi
         bybar=ypi2*yp
   
         cy=yp*(yp2-1.0)
         cyi=ypi*(ypi2-1.0)
         cyd=3.0*yp2-1.0
         cydi=-3.0*ypi2+1.0
                  
         sixth_hx2 = sixth*hx2
         sixth_hy2 = sixth*hy2
         z36th_hx2_hy2 = z36th*hx2*hy2
         
         sixth_hx = sixth*hx
         sixth_hxi_hy2 = sixth*hxi*hy2
         z36th_hx_hy2 = z36th*hx*hy2
         
         sixth_hx2_hyi = sixth*hx2*hyi
         sixth_hy = sixth*hy
         z36th_hx2_hy = z36th*hx2*hy
         
         ip1 = i+1
         jp1 = j+1
         
         do k=1,nv
            ! bicubic spline interpolation
            fval(k) = &
                  xpi*( &
                     ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                     +xp*(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1)) &
                  +sixth_hx2*( &
                     cxi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                     cx*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
                  +sixth_hy2*( &
                     xpi*(cyi*fin(3,k,i,j) +cy*fin(3,k,i,jp1))+ &
                     xp*(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
                  +z36th_hx2_hy2*( &
                     cxi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                     cx*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

            ! derivatives of bicubic splines
            df_dx(k) = &
                  hxi*( &
                     -(ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                     +(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1))) &
                  +sixth_hx*( &
                     cxdi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                     cxd*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
                  +sixth_hxi_hy2*( &
                     -(cyi*fin(3,k,i,j)  +cy*fin(3,k,i,jp1)) &
                     +(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
                  +z36th_hx_hy2*( &
                     cxdi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                     cxd*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

            df_dy(k) = &
                  hyi*( &
                     xpi*(-fin(1,k,i,j) +fin(1,k,i,jp1))+ &
                     xp*(-fin(1,k,ip1,j)+fin(1,k,ip1,jp1))) &
                  +sixth_hx2_hyi*( &
                     cxi*(-fin(2,k,i,j) +fin(2,k,i,jp1))+ &
                     cx*(-fin(2,k,ip1,j)+fin(2,k,ip1,jp1))) &
                  +sixth_hy*( &
                     xpi*(cydi*fin(3,k,i,j) +cyd*fin(3,k,i,jp1))+ &
                     xp*(cydi*fin(3,k,ip1,j)+cyd*fin(3,k,ip1,jp1))) &
                  +z36th_hx2_hy*( &
                     cxi*(cydi*fin(4,k,i,j) +cyd*fin(4,k,i,jp1))+ &
                     cx*(cydi*fin(4,k,ip1,j)+cyd*fin(4,k,ip1,jp1)))
         end do
         
      end subroutine Do_EoS_Interpolations


      subroutine Get_HELM_PC_Results( &
            rq, Z, X, abar, zbar,  &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, basic_flag,  &
            res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
            d_dabar_c_TRho, d_dzbar_c_TRho, debug, ierr)
         use helm
         use const_def
         
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(inout) :: Rho, logRho, T, logT
         logical, intent(in) :: basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)
         logical, intent(in) :: debug
         integer, intent(out) :: ierr
         
         !real(dp), parameter :: logG0 = 6.67173d0 ! for pure o16
         real(dp) :: helm_res(num_helm_results)
         real(dp), dimension(num_eos_basic_results) ::  &
            res1, d_dlnRho_c_T1, d_dlnT_c_Rho1, d_dabar_c_TRho1, d_dzbar_c_TRho1, &
            res2, d_dlnRho_c_T2, d_dlnT_c_Rho2
         real(dp) :: logG0, logG, logG_lo, logG_hi, alfa, beta, sin_pi_alfa, &
            A, B, dA_dlnT, dA_dlnRho, dB_dlnT, dB_dlnRho, dlogG_dlogT, dlogG_dlogRho, &
            d_alfa_dlogT, d_alfa_dlogRho, d_beta_dlogT, d_beta_dlogRho,  &
            logT_lo, logT_hi, logRho_lo, logRho_hi
     
         logical :: dbg
         
         include 'formats.dek'
         
         dbg = debug
         
         
         !gamma = (zbar*qe)**2*((4.0d0/3.0d0)*pi*avo*Rho/abar)**one_third/(kerg*T)
         !logG = logG0 + logRho/3 - logT
         ! where G0 = (zbar*qe)**2*((4.0d0/3.0d0)*pi*avo/abar)**one_third/kerg
         logG0 = log10_cr( &
            (zbar*qe)*(zbar*qe)*pow_cr((4.0d0/3.0d0)*pi*avo/abar,one_third)/kerg)
         
         !logG0 = 6.67173d0 ! for pure o16
         
         logG = logG0 + logRho/3 - logT
         logG_lo = rq% log_Gamma_all_HELM ! HELM for logG <= this
         logG_hi = rq% log_Gamma_all_PC ! PC for logG >= this
         logRho_lo = rq% logRho2_PC_limit ! don't use PC for logRho < this
         logRho_hi = rq% logRho1_PC_limit ! okay for pure PC for logRho > this
         logT_lo = logG0 + logRho_lo/3 - logG_lo
         logT_hi = logG0 + logRho_hi/3 - logG_hi
         
         if (dbg) then
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,1) 'logG', logG
            write(*,1) 'logG_lo', logG_lo
            write(*,1) 'logG_hi', logG_hi
            write(*,1) 'logRho_lo', logRho_lo
            write(*,1) 'logRho_hi', logRho_hi
            write(*,1) 'logT_lo', logT_lo
            write(*,1) 'logT_hi', logT_hi
            write(*,*)
            !stop
         end if
         
         d_alfa_dlogT=0
         d_alfa_dlogRho=0
         if (Z < rq% PC_min_Z .or. logRho <= logRho_lo .or. logG <= logG_lo) then
            alfa = 1 ! pure HELM
         else if (logRho >= logRho_hi .and. logG >= logG_hi) then
            alfa = 0 ! pure PC
         else if (logT >= logT_hi) then ! blend in logG
            alfa = (logG_hi - logG)/(logG_hi - logG_lo)
            dlogG_dlogT = -1
            dlogG_dlogRho = 1d0/3d0
            d_alfa_dlogT = -dlogG_dlogT/(logG_hi - logG_lo)
            d_alfa_dlogRho = -dlogG_dlogRho/(logG_hi - logG_lo)
         else ! blend in logRho
            if (logT >= logT_lo) logRho_lo = (logT_lo + logG_lo - logG0)*3
            alfa = (logRho_hi - logRho)/(logRho_hi - logRho_lo)
            d_alfa_dlogRho = -1/(logRho_hi - logRho_lo)
         end if

         beta = 1 - alfa
         d_beta_dlogT = -d_alfa_dlogT
         d_beta_dlogRho = -d_alfa_dlogRho
         
         if (alfa /= 0) then
            call Get_HELMEOS_Results( &
               rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
               rq% include_radiation, rq% always_skip_elec_pos, basic_flag,  &
               res1, d_dlnRho_c_T1, d_dlnT_c_Rho1, &
               d_dabar_c_TRho1, d_dzbar_c_TRho1, helm_res, ierr)
            if (ierr /= 0) return
         end if

         if (beta /= 0) then            
            call Get_PC_Results( &
               rq, Z, X, abar, zbar,  &
               species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, basic_flag,  &
               res2, d_dlnRho_c_T2, d_dlnT_c_Rho2, ierr)
            if (ierr /= 0) then
               return
               write(*,*) 'ierr from Get_PC_Results in Get_HELM_PC_Results'
               write(*,1) 'logG_lo', logG_lo
               write(*,1) 'logG', logG
               write(*,1) 'logG_hi', logG_hi
               write(*,1) 'beta', beta
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               !write(*,1) '', 
               stop 'debug eos'
            end if
         end if
      
         if (beta == 0d0) then
            res = res1
            if (.not. basic_flag) then
               d_dlnRho_c_T = d_dlnRho_c_T1
               d_dlnT_c_Rho = d_dlnT_c_Rho1
               d_dabar_c_TRho = d_dabar_c_TRho1
               d_dzbar_c_TRho = d_dzbar_c_TRho1
            end if
            return
         end if
         
         if (alfa == 0d0) then
            res = res2
            if (.not. basic_flag) then
               d_dlnRho_c_T = d_dlnRho_c_T2
               d_dlnT_c_Rho = d_dlnT_c_Rho2
               d_dabar_c_TRho = 0d0
               d_dzbar_c_TRho = 0d0
            end if
            return
         end if
         
         A = alfa ! (1 - cospi_cr(alfa))/2 ! smooth the transitions
         B = 1 - A
                  
         res = A*res1 + B*res2
         
         !write(*,1) 'Get_HELM_PC_Results grad_ad', res(i_grad_ad), res1(i_grad_ad), res2(i_grad_ad), A, B
         
         if (basic_flag) return
         
         ! just use approximate partials for these
         d_dlnRho_c_T = A*d_dlnRho_c_T1 + B*d_dlnRho_c_T2
         d_dlnT_c_Rho = A*d_dlnT_c_Rho1 + B*d_dlnT_c_Rho2
         d_dabar_c_TRho = A*d_dabar_c_TRho1
         d_dzbar_c_TRho = A*d_dzbar_c_TRho1
         
      end subroutine Get_HELM_PC_Results


      subroutine Get_PC_Results( &
            rq, Z, X, abar, zbar,  &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, basic_flag,  &
            res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         use pc_eos
         use chem_def, only: chem_isos
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: Rho, logRho, T, logT
         logical, intent(in) :: basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         integer, intent(out) :: ierr
         
         real(dp) :: Ni, N, lnT, lnRho, dlnT, dlnRho, Rho2, T3
         integer :: i, j
         real(dp), parameter :: lnT_frac = 1d-7, lnRho_frac = 1d-7
         real(dp), dimension(species) :: AY, AZion, ACMI
         real(dp), dimension(num_eos_basic_results) :: &
            res2, d_dlnRho_c_T2, d_dlnRho_c_Rho2, &
            res3, d_dlnRho_c_T3, d_dlnRho_c_Rho3
         
         include 'formats.dek'
         
         ierr = 0

         Ni = rho/(abar*mp)
         N = Ni/rho
         AZion(1:species) = chem_isos% Z(chem_id(1:species))
         ACMI(1:species) = chem_isos% W(chem_id(1:species)) ! this really is atomic weight.
         do j=1,species
            if (xa(j) < rq% mass_fraction_limit_for_PC) then
               AY(j) = 0
            else
               AY(j) = xa(j)/ACMI(j)
            end if
         end do
         
         lnT = logT*ln10
         lnRho = logRho*ln10
         
         dlnT = lnT*lnT_frac
         dlnRho = lnRho*lnRho_frac
         
         call do1(.false.,Rho,T,res,ierr)
         if (ierr /= 0) return
         
         !write(*,1) 'Get_PC_Results grad_ad', res(i_grad_ad), log10_cr(Rho), log10_cr(T)
         
         Rho2 = exp_cr(lnRho+dlnRho)
         call do1(.false.,Rho2,T,res2,ierr)
         if (ierr /= 0) return
         
         T3 = exp_cr(lnT+dlnT)
         call do1(.false.,Rho,T3,res3,ierr)
         if (ierr /= 0) return
         
         d_dlnRho_c_T(:) = (res2(:)-res(:))/dlnRho
         d_dlnT_c_Rho(:) = (res3(:)-res(:))/dlnT

         contains
         
         
#ifdef offload
         !dir$ attributes offload : mic :: do1
#endif      
         subroutine do1(show,RHO,T,res,ierr)
            logical, intent(in) :: show
            real(dp), intent(in) :: RHO,T
            real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
            integer, intent(out) :: ierr

            integer :: j, LIQSOL
            real(dp) :: AY_priv(species), dsp, dse, dpe, &
               DENS,Zmean,CMImean,Z2mean,GAMI,CHI,TPT,GAMImean,TEMP, &
               PnkT,UNkT,SNk,CVNkt,CV,CHIR,CHIT,Tnk,P,Cp,gamma1,gamma3,grad_ad, &
               Prad,Pgas,PRADnkT,dE_dRho,dS_dT,dS_dRho,mu,lnfree_e
            real(dp), parameter :: UN_T6=0.3157746d0
            
            include 'formats.dek'
            
            ierr = 0
            
            AY_priv(1:species) = AY(1:species)
         
            TEMP=T*1d-6/UN_T6 ! T [au]

            if (show) then
               write(*,1) 'RHO', RHO
               write(*,1) 'TEMP', TEMP
               write(*,1) 'AZion(1:species)', AZion(1:species)
               write(*,1) 'ACMI(1:species)', ACMI(1:species)
               write(*,1) 'AY(1:species)', AY_priv(1:species)
               write(*,1) 'xa(1:species)', xa(1:species)
            end if

            call MELANGE9( &
               species,AY_priv,AZion,ACMI,RHO,TEMP,PRADnkT, &
               DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL, &
               PnkT,UNkT,SNk,CVNkt,CHIR,CHIT,ierr)
            
            if (ierr /= 0) then
               return
               write(*,1) 'RHO', RHO
               write(*,1) 'T', T
               write(*,1) 'logRho', log10_cr(RHO)
               write(*,1) 'logT', log10_cr(T)
               write(*,*) 'ierr from MELANGE9'
               stop 'debug eos'
            end if
            
            if (show) then
               write(*,1) 'PRADnkT', PRADnkT
               write(*,1) 'DENS', DENS
               write(*,1) 'Zmean', Zmean
               write(*,1) 'CMImean', CMImean
               write(*,1) 'Z2mean', Z2mean
               write(*,1) 'GAMI', GAMImean
               write(*,1) 'CHI', CHI
               write(*,1) 'TPT', TPT
               write(*,1) 'PnkT', PnkT
               write(*,1) 'UNkT', UNkT
               write(*,1) 'SNk', SNk
               write(*,1) 'CV', CVNkt
               write(*,1) 'CHIR', CHIR
               write(*,1) 'CHIT', CHIT
               write(*,*)
            end if
            
            Tnk=8.31447d7/CMImean*RHO*T ! n_i kT [erg/cc]
            Pgas = PnkT*Tnk
            if (rq% include_radiation) then
               ! results from MELANGE9 do not include radiation.  add it now.
               CHIT=CHIT*(PnkT/(PnkT+PRADnkT))  + 4.*PRADnkT/(PnkT+PRADnkT)
               PnkT=PnkT+PRADnkT
               UNkT=UNkT+3.*PRADnkT
               SNk=SNk+4.*PRADnkT
               CVNkt=CVNkt+12.*PRADnkT
            end if
            P = PnkT*Tnk
            
            Ni = rho/(abar*mp)
            N = Ni/rho
            CV = CVNkt*N*kerg
            gamma3 = 1d0 + P/rho * chit/(T*CV)
            gamma1 = chit*(gamma3-1d0) + chir
            grad_ad = (gamma3-1d0)/gamma1
            Cp = CV * gamma1/chir
            dE_dRho = (1d0-chiT)*P/(rho*rho)
            dS_dT = CV/T
            dS_dRho = -P*chiT/(rho*rho * T)
            mu = abar / (1 + zbar)
            lnfree_e = log_cr(zbar/abar) ! for complete ionization
            
            res(i_lnPgas) = log_cr(Pgas)
            res(i_lnE) = safe_log_cr(UNkt*N*kerg*T)
            res(i_lnS) = safe_log_cr(SNk*N*kerg)

            res(i_grad_ad) = grad_ad
            res(i_chiRho) = CHIR
            res(i_chiT) = CHIT
            res(i_Cp) = Cp
            res(i_Cv) = CV
            res(i_dE_dRho) = dE_dRho
            res(i_dS_dT) = dS_dT
            res(i_dS_dRho) = dS_dRho
            res(i_mu) = mu
            res(i_lnfree_e) = lnfree_e
            res(i_gamma1) = gamma1
            res(i_gamma3) = gamma3
            res(i_eta) = CHI
            
            if (.false.) then ! check consistency
               dse = T*res(i_dS_dT)/res(i_Cv) - 1d0
               dpe = rho*rho*res(i_dE_dRho)/P + res(i_chiT) - 1d0
               dsp = -rho*rho*res(i_dS_dRho)/(res(i_chiT)*P/T) - 1d0
               if (abs(dse) > 1d-6 .or. abs(dpe) > 1d-6 .or. abs(dsp) > 1d-6) then
                  write(*,1) 'dse', dse
                  write(*,1) 'dpe', dpe
                  write(*,1) 'dsp', dsp
                  stop 'debug eos'
               end if
            end if
            
            
            if (UNkt < 0 .or. is_bad_num(res(i_lnE))) then
               if (.false. .and. RHO < 5d2) then
!$omp critical (PC_eos_problem)
                  write(*,*) 'PC eos returned negative energy -- perhaps means logRho limit too small'
                  if (.false.) then
                     write(*,1) 'res(i_lnE)', res(i_lnE)
                     write(*,1) 'T', T
                     write(*,1) 'RHO', RHO
                     write(*,1) 'logT', log10_cr(T)
                     write(*,1) 'logRho', log10_cr(RHO)
                     write(*,1) 'N', N
                     write(*,1) 'PRADnkT', PRADnkT
                     write(*,1) 'DENS', DENS
                     write(*,1) 'Zmean', Zmean
                     write(*,1) 'CMImean', CMImean
                     write(*,1) 'Z2mean', Z2mean
                     write(*,1) 'GAMI', GAMImean
                     write(*,1) 'CHI', CHI
                     write(*,1) 'TPT', TPT
                     write(*,1) 'PnkT', PnkT
                     write(*,1) 'UNkT', UNkT
                     write(*,1) 'SNk', SNk
                     write(*,1) 'CV', CVNkt
                     write(*,1) 'CHIR', CHIR
                     write(*,1) 'CHIT', CHIT
                     write(*,*)
                  end if
!$omp end critical (PC_eos_problem)
               end if
               ierr = -1
               return
               stop 'Get_PC_Results' 
            end if

         end subroutine do1
     
      end subroutine Get_PC_Results

         
      subroutine get_eosDT_HELMEOS_Results( &
               rq, Z_in, X_in, abar, zbar, &
               arho, alogrho, atemp, alogtemp, &
               include_radiation, always_skip_elec_pos, basic_flag,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, helm_res, ierr)
         use helm
         use const_def
         use eosDT_load_tables, only: Load_eosDT_Table
         
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z_in, X_in, abar, zbar
         real(dp), intent(in) :: arho, alogrho, atemp, alogtemp
         logical, intent(in) :: include_radiation, always_skip_elec_pos, basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: helm_res(:) ! (num_helm_results)
         integer, intent(out) :: ierr
         
         real(dp) :: Z, X, Rho, logRho, T, logT
         real(dp), parameter :: tiny = 1d-20
         
         !logical, parameter :: clip_to_table_boundaries = .false.
         logical, parameter :: clip_to_table_boundaries = .true.

         include 'formats'
         
         if (.not. eosDT_is_initialized) then
            call Load_eosDT_Table(ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'Load_eosDT_Table ierr', ierr
               return
            end if
         end if
         
         if (is_bad_num(X_in) .or. is_bad_num(Z_in)) then
            ierr = -1
            return
         end if
         
         X = X_in; Z = Z_in
         if (X < tiny) X = 0
         if (Z < tiny) Z = 0
         
         !..get temp and rho args
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (alogtemp == arg_not_provided) logT = log10_cr(T)
         if (atemp == arg_not_provided) T = exp10_cr(logT)
         
         if (T <= 0) then
            ierr = -1
            return
         end if
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -3; return
         end if
         if (alogrho == arg_not_provided) logRho = log10_cr(Rho)
         if (arho == arg_not_provided) Rho = exp10_cr(logRho)
         
         if (Rho <= 0) then
            ierr = -1
            return
         end if
         
         if (is_bad_num(Rho) .or. is_bad_num(T)) then
            ierr = -1
            return
         end if
         
         !write(*,1) 'Rho', Rho
         !write(*,1) 'T', T
         !stop
         
         call helmeos2( &
               T, logT, Rho, logRho, X, abar, zbar, helm_res,  &
               clip_to_table_boundaries, include_radiation, &
               always_skip_elec_pos, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in helmeos2'
            return
         end if
         
         call do_convert_helm_results( &
               helm_res, Z, X, abar, zbar, Rho, T, basic_flag,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in do_convert_helm_results'
            return
         end if

!         call Get_HELMEOS_Results( &
!               rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
!               include_radiation, always_skip_elec_pos, basic_flag,  &
!               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
!               d_dabar_c_TRho, d_dzbar_c_TRho, helm_res, ierr)
               
      end subroutine get_eosDT_HELMEOS_Results
         
         
      subroutine Get_HELMEOS_Results( &
               rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
               include_radiation, always_skip_elec_pos, basic_flag,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, helm_res, ierr)
         use helm
         use const_def
         
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         real(dp), intent(inout) :: Rho, logRho, T, logT
         logical, intent(in) :: include_radiation, always_skip_elec_pos, basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: helm_res(:) ! (num_helm_results)
         integer, intent(out) :: ierr
         
         logical, parameter :: clip_to_table_boundaries = .true.
         
         include 'formats.dek'

         ierr = 0
         
         call helmeos2( &
               T, logT, Rho, logRho, X, abar, zbar, helm_res,  &
               clip_to_table_boundaries, include_radiation, &
               always_skip_elec_pos, ierr)
         if (ierr /= 0) then
            return
         end if
         
         call do_convert_helm_results( &
               helm_res, Z, X, abar, zbar, Rho, T, basic_flag,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in do_convert_helm_results'
            return
         end if
                  
      end subroutine Get_HELMEOS_Results


      subroutine do_convert_helm_results( &
               helm_res, Z, X, abar, zbar, Rho, T, basic_flag,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         use helm
         use chem_def
         real(dp), intent(in) :: helm_res(:) ! (num_helm_results)
         real(dp), intent(in) :: Z, X, abar, zbar, Rho, T
         logical, intent(in) :: basic_flag
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)
         integer, intent(out) :: ierr

         real(dp) :: mu, P, Pgas, free_e, dse, dpe, dsp
         integer :: j, ci
         
         include 'formats.dek'
         
         ierr = 0

         ! store the basic results first

         P = helm_res(h_ptot)
         
         Pgas = helm_res(h_pgas)
         res(i_lnPgas) = log_cr(Pgas)
         res(i_lnE) = log_cr(helm_res(h_etot))
         res(i_lnS) = log_cr(helm_res(h_stot))

         res(i_grad_ad) = helm_res(h_nabad)
         res(i_chiRho) = helm_res(h_chid)
         res(i_chiT) = helm_res(h_chit)
         res(i_Cp) = helm_res(h_cp)
         res(i_Cv) = helm_res(h_cv)
         res(i_dE_dRho) = helm_res(h_ded)
         res(i_dS_dT) = helm_res(h_dst)
         res(i_dS_dRho) = helm_res(h_dsd)
         mu = abar / (1 + zbar)
         res(i_mu) = mu
         free_e = max(1d-99, helm_res(h_xne)) / (avo * Rho) ! assuming complete ionization
         res(i_lnfree_e) = log_cr(free_e)
         res(i_gamma1) = helm_res(h_gam1)
         res(i_gamma3) = helm_res(h_gam3)
         res(i_eta) = helm_res(h_etaele)
         
         if (basic_flag) return
         
         d_dlnRho_c_T(i_lnPgas) = helm_res(h_dpgasd) * Rho / Pgas
         d_dlnRho_c_T(i_lnE) = helm_res(h_ded) * Rho / helm_res(h_etot)
         d_dlnRho_c_T(i_lnS) = helm_res(h_dsd) * Rho / helm_res(h_stot)
         d_dlnRho_c_T(i_grad_ad) = helm_res(h_dnabdd) * Rho
         d_dlnRho_c_T(i_chiRho) = helm_res(h_dchiddd) * Rho
         d_dlnRho_c_T(i_chiT) = helm_res(h_dchitdd) * Rho
         d_dlnRho_c_T(i_Cp) = helm_res(h_dcpdd) * Rho
         d_dlnRho_c_T(i_Cv) = helm_res(h_dcvdd) * Rho
         d_dlnRho_c_T(i_dE_dRho) = helm_res(h_dedd) * Rho
         d_dlnRho_c_T(i_dS_dT) = helm_res(h_dsdt) * Rho
         d_dlnRho_c_T(i_dS_dRho) = helm_res(h_dsdd) * Rho
         d_dlnRho_c_T(i_mu) = 0
         d_dlnRho_c_T(i_lnfree_e) = (helm_res(h_dxned) / avo - res(i_lnfree_e)) / free_e
         d_dlnRho_c_T(i_gamma1) = helm_res(h_dgam1dd) * Rho
         d_dlnRho_c_T(i_gamma3) = helm_res(h_dgam3dd) * Rho
         d_dlnRho_c_T(i_eta) = helm_res(h_detad) * Rho
           
         d_dlnT_c_Rho(i_lnPgas) = helm_res(h_dpgast) * T / Pgas
         d_dlnT_c_Rho(i_lnE) = helm_res(h_det) * T / helm_res(h_etot)
         d_dlnT_c_Rho(i_lnS) = helm_res(h_dst) * T / helm_res(h_stot)
         d_dlnT_c_Rho(i_grad_ad) = helm_res(h_dnabdt) * T
         d_dlnT_c_Rho(i_chiRho) = helm_res(h_dchiddt) * T
         d_dlnT_c_Rho(i_chiT) = helm_res(h_dchitdt) * T
         d_dlnT_c_Rho(i_Cp) = helm_res(h_dcpdt) * T
         d_dlnT_c_Rho(i_Cv) = helm_res(h_dcvdt) * T
         d_dlnT_c_Rho(i_dE_dRho) = helm_res(h_dedt) * T
         d_dlnT_c_Rho(i_dS_dT) = helm_res(h_dstt) * T
         d_dlnT_c_Rho(i_dS_dRho) = helm_res(h_dsdt) * T
         d_dlnT_c_Rho(i_mu) = 0
         d_dlnT_c_Rho(i_lnfree_e) = (helm_res(h_dxnet) * T / (avo * Rho)) / free_e
         d_dlnT_c_Rho(i_gamma1) = helm_res(h_dgam1dt) * T
         d_dlnT_c_Rho(i_gamma3) = helm_res(h_dgam3dt) * T
         d_dlnT_c_Rho(i_eta) = helm_res(h_detat) * T
           
         d_dabar_c_TRho(i_lnPgas) = helm_res(h_dpgasa) / Pgas
         d_dabar_c_TRho(i_lnE) = helm_res(h_dea) / helm_res(h_etot)
         d_dabar_c_TRho(i_lnS) = helm_res(h_dsa) / helm_res(h_stot)
         d_dabar_c_TRho(i_grad_ad) = helm_res(h_dnabda)
         d_dabar_c_TRho(i_chiRho) = helm_res(h_dchidda)
         d_dabar_c_TRho(i_chiT) = helm_res(h_dchitda)
         d_dabar_c_TRho(i_Cp) = helm_res(h_dcpda)
         d_dabar_c_TRho(i_Cv) = helm_res(h_dcvda)
         d_dabar_c_TRho(i_dE_dRho) = helm_res(h_deda)
         d_dabar_c_TRho(i_dS_dT) = helm_res(h_dsta)
         d_dabar_c_TRho(i_dS_dRho) = helm_res(h_dsda)
         d_dabar_c_TRho(i_mu) = 0
         d_dabar_c_TRho(i_lnfree_e) = (helm_res(h_dxnea) / (avo * Rho)) / free_e
         d_dabar_c_TRho(i_gamma1) = helm_res(h_dgam1da)
         d_dabar_c_TRho(i_gamma3) = helm_res(h_dgam3da)
         d_dabar_c_TRho(i_eta) = helm_res(h_detaa)
           
         d_dzbar_c_TRho(i_lnPgas) = helm_res(h_dpgasz) / Pgas
         d_dzbar_c_TRho(i_lnE) = helm_res(h_dez) / helm_res(h_etot)
         d_dzbar_c_TRho(i_lnS) = helm_res(h_dsz) / helm_res(h_stot)
         d_dzbar_c_TRho(i_grad_ad) = helm_res(h_dnabdz)
         d_dzbar_c_TRho(i_chiRho) = helm_res(h_dchiddz)
         d_dzbar_c_TRho(i_chiT) = helm_res(h_dchitdz)
         d_dzbar_c_TRho(i_Cp) = helm_res(h_dcpdz)
         d_dzbar_c_TRho(i_Cv) = helm_res(h_dcvdz)
         d_dzbar_c_TRho(i_dE_dRho) = helm_res(h_dedz)
         d_dzbar_c_TRho(i_dS_dT) = helm_res(h_dstz)
         d_dzbar_c_TRho(i_dS_dRho) = helm_res(h_dsdz)
         d_dzbar_c_TRho(i_mu) = 0
         d_dzbar_c_TRho(i_lnfree_e) = (helm_res(h_dxnez) / (avo * Rho)) / free_e
         d_dzbar_c_TRho(i_gamma1) = helm_res(h_dgam1dz)
         d_dzbar_c_TRho(i_gamma3) = helm_res(h_dgam3dz)
         d_dzbar_c_TRho(i_eta) = helm_res(h_detaz)

      end subroutine do_convert_helm_results
      

      subroutine get_T( &
               handle, Z, X, abar, zbar,  &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess,  &
               logT_bnd1, logT_bnd2,  other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar
            ! mean atomic number (nucleons per nucleus; grams per mole)
         real(dp), intent(in) :: zbar ! mean charge per nucleus
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: logRho ! log10 of density
         integer, intent(in) :: which_other ! from eos_def.  e.g., i_P for pressure
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logT_tol
         integer, intent(in) :: max_iter ! max number of iterations        

         real(dp), intent(in) :: logT_guess
         real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            ! set to arg_not_provided if do not know bounds
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in c_def)
         
         real(dp), intent(out) :: logT_result
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)
         
         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         logical, parameter :: doing_Rho = .false.
         
         call do_safe_get_Rho_T( &
               handle, Z, X, abar, zbar,  &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, doing_Rho,  &
               logT_guess, logT_result, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_tol, other_tol, max_iter, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
      
      end subroutine get_T
      

      subroutine get_Rho( &
               handle, Z, X, abar, zbar,  &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess,  &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
     
         use const_def
         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar
            ! mean atomic number (nucleons per nucleus; grams per mole)
         real(dp), intent(in) :: zbar ! mean charge per nucleus
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: logT ! log10 of temperature

         integer, intent(in) :: which_other ! from eos_def.
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logRho_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logRho_guess
         real(dp), intent(in) :: logRho_bnd1, logRho_bnd2 ! bounds for logrho
            ! set to arg_not_provided if do not know bounds
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in c_def)
            
         real(dp), intent(out) :: logRho_result
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         logical, parameter :: doing_Rho = .true.
         real(dp) :: Prad
         
         call do_safe_get_Rho_T( &
               handle, Z, X, abar, zbar,  &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, doing_Rho,  &
               logRho_guess, logRho_result, logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_tol, other_tol, max_iter, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)

      end subroutine get_Rho
      
      
      subroutine do_safe_get_Rho_T( &
               handle, Z, XH1, abar, zbar,  &
               species, chem_id, net_iso, xa, &
               the_other_log, which_other, other_value, doing_Rho,  &
               initial_guess, x, xbnd1, xbnd2, other_at_bnd1, other_at_bnd2,  &
               xacc, yacc, ntry,  &
               res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
         use const_def
         use chem_def, only: num_chem_isos
         use utils_lib, only: is_bad_num
         use num_lib, only: safe_root_without_brackets
         use eosDT_load_tables, only: Load_eosDT_Table
         integer, intent(in) :: handle
         real(dp), intent(in) :: Z, XH1, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         integer, intent(in) :: which_other
         real(dp), intent(in) :: other_value
         logical, intent(in) :: doing_Rho
         real(dp), intent(in) :: initial_guess ! for x
         real(dp), intent(out) :: x ! if doing_Rho, then logRho, else logT
         real(dp), intent(in) :: the_other_log
         real(dp), intent(in) :: xbnd1, xbnd2, other_at_bnd1, other_at_bnd2
         real(dp), intent(in) :: xacc, yacc ! tolerances
         integer, intent(in) :: ntry ! max number of iterations        
         real(dp), intent(out) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dabar_c_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(out) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results)
         integer, intent(out) :: eos_calls, ierr
         
         integer :: i, j
         integer, parameter :: lrpar = 0, lipar = 0, newt_imax = 6
         real(dp), parameter :: dx = 0.1d0
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)
         real(dp) :: the_other_val, logRho, logT, rho, T
         type (EoS_General_Info), pointer :: rq
         logical, parameter :: basic_flag = .false.

         include 'formats'
         
         ierr = 0

         if (.not. eosDT_is_initialized) then
            call Load_eosDT_Table(ierr)
            if (ierr /= 0) return
         end if
			
         call get_eos_ptr(handle, rq, ierr)
         if (ierr /= 0) then
            write(*, *) 'get_eos_ptr returned ierr', ierr
            return
         end if
         
         eos_calls = 0
			the_other_val = exp10_cr(the_other_log)
         nullify(ipar, rpar)

         x = safe_root_without_brackets( &
            f, initial_guess, dx, min(ntry,newt_imax), ntry, xacc, yacc, &
            lrpar, rpar, lipar, ipar, ierr)
         
         
         contains
         
         
         real(dp) function f(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
            ! returns with ierr = 0 if was able to evaluate f and df/dx at x
            ! if df/dx not available, it is okay to set it to 0
            use const_def, only: dp
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(in) :: x
            real(dp), intent(out) :: dfdx
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr
            
            include 'formats'
            ierr = 0
            eos_calls = eos_calls + 1
            f = 0; dfdx = 0
			
   			if (doing_Rho) then
   			   logRho = x
   			   rho = exp10_cr(logRho)
   			   logT = the_other_log
   			   T = the_other_val
   			else
   			   logT = x
   			   T = exp10_cr(logT)
   			   logRho = the_other_log
   			   rho = the_other_val
   			end if
   			
            
            call Get_eosDT_Results(rq, Z, XH1, abar, zbar,  &
                  species, chem_id, net_iso, xa, &
                  rho, logRho, T, logT, &
                  basic_flag, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
                  d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
            if (ierr /= 0) then
               if (.false.) then
                  write(*,*) 'Get_eosDT_Results returned ierr', ierr
                  write(*,1) 'Z', Z
                  write(*,1) 'XH1', XH1
                  write(*,1) 'abar', abar
                  write(*,1) 'zbar', zbar
                  write(*,1) 'rho', rho
                  write(*,1) 'logRho', logRho
                  write(*,1) 'T', T
                  write(*,1) 'logT', logT
                  write(*,*)
                  stop
               end if
               return
            end if
         
   			f = res(which_other) - other_value
			
   			if (doing_Rho) then
   			   dfdx = d_dlnRho_c_T(which_other)*ln10
   			else
   			   dfdx = d_dlnT_c_Rho(which_other)*ln10
   			end if
			   !write(*,2) 'logT logRho val/ln10 target/ln10 f dfdx', &
			   !   eos_calls, logT, logRho, res(which_other)/ln10, other_value/ln10, f, dfdx
            
         end function f
         
      end subroutine do_safe_get_Rho_T

      
      real(dp) function eval_theta_e(eta, d_theta_e_deta)
         use interp_1d_lib, only : interp_values_and_slopes
         real(dp), intent(in) :: eta
         real(dp), intent(out) :: d_theta_e_deta
         integer, parameter :: nv = 1
         integer :: ierr
         real(dp) :: x(nv), vals(nv), slopes(nv)
         
         x(1) = eta
         call interp_values_and_slopes( &
                  x_theta_e, theta_e_nx, f_theta_e1, nv, x, vals, slopes, ierr)
         eval_theta_e = vals(1)
         d_theta_e_deta = slopes(1)
         
      end function eval_theta_e 

#ifdef offload
      !dir$ end options
#endif

      end module eosDT_eval
      
