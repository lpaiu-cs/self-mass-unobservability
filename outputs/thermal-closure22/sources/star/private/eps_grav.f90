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

      module eps_grav
      
      use star_private_def
      use const_def
      use chem_def, only: chem_isos
      use utils_lib, only:is_bad_num
      
      implicit none


      contains
      
      
      subroutine prepare_for_eval_eps_grav_and_partials(s, ierr)
         use alloc
         use interp_1d_def
         use interp_1d_lib
         
         type (star_info), pointer :: s 
         integer, intent(out) :: ierr
         
         real(dp) :: max_gam
         integer :: nwork, prev_nz, sz, nz, k
         real(dp), pointer, dimension(:) :: xm, prev_xq, xq, work, f1
         real(dp), pointer :: f(:,:)
         logical :: do_mu, do_lnS
         
         include 'formats'
         ierr = 0
         
         if (s% E_flag) then
            write(*,*) 'eps_grav not allowed with E_flag'
            stop 1
         end if
         
         !write(*,2) 's% k_below_Eulerian_eps_grav', s% k_below_Eulerian_eps_grav
         !write(*,2) 's% k_Lagrangian_eps_grav', s% k_Lagrangian_eps_grav
         !write(*,2) 's% nz', s% nz
         
         do_mu = (s% include_dmu_dt_in_eps_grav .and. .not. s% have_prev_mu)
         
         do_lnS = .false.
         if (.not. s% have_prev_lnS) then
            max_gam = maxval(s% gam(1:s% nz))
            if (max_gam > s% Gamma_lnS_eps_grav_full_off .or. &
                s% use_lnS_for_eps_grav) do_lnS = .true.
         end if
         
         if ((.not. do_mu) .and. (.not. do_lnS)) return
         
         nwork = pm_work_size
         prev_nz = s% prev_mesh_nz
         nz = s% nz
         sz = nwork*prev_nz
         
         nullify(xm, prev_xq, xq, work, f1)
            
         call non_crit_get_work_array( &
            s, xm, nz, 0, 'prepare_for_eval_eps_grav_and_partials', ierr)
         if (failed('xm')) return
         
         call non_crit_get_work_array( &
            s, work, sz, 0, 'prepare_for_eval_eps_grav_and_partials', ierr)
         if (failed('work')) return
         
         call non_crit_get_work_array( &
            s, f1, 4*prev_nz, 0, 'prepare_for_eval_eps_grav_and_partials', ierr)
         if (failed('f1')) return
         
         f(1:4,1:prev_nz) => f1(1:4*prev_nz)
         
         xm(1) = s% xmstar_old - s% xmstar
         do k=2,nz
            xm(k) = xm(k-1) + s% dm(k-1)
         end do
         
         if (do_mu) then   

            do k=1,prev_nz
               f(1,k) = s% prev_mesh_mu(k)
            end do
                
            call non_crit_do1_alloc_if_necessary( &
               s, s% prev_mu, nz, 'prepare_for_eval_eps_grav_and_partials', ierr)
            if (failed('prev_mu')) return
            
            call interp_pm(s% prev_mesh_xm, prev_nz, f1, &
               nwork, work, 'prepare_for_eval_eps_grav_and_partials', ierr)
            if (failed('interp_pm mu')) return    
                   
            call interp_values( &
               s% prev_mesh_xm, prev_nz, f1, nz, xm, s% prev_mu, ierr)
            if (failed('interp_values mu')) return
            
            s% have_prev_mu = .true. 
                    
         end if
         
         
         if (do_lnS) then
         
            do k=1,prev_nz
               f(1,k) = s% prev_mesh_lnS(k)
            end do
            
            if (s% k_below_Eulerian_eps_grav < s% nz + 1) then ! some cells use Lagrangian
               ! do prev_lnS
               !write(*,*) 'set prev_lnS'
            
               call non_crit_do1_alloc_if_necessary( &
                  s, s% prev_lnS, nz, 'prepare_for_eval_eps_grav_and_partials', ierr)
               if (failed('prev_lnS')) return
            
               if (s% k_const_mass > 1) then            
                  call interp_pm(s% prev_mesh_xm, prev_nz, f1, &
                     nwork, work, 'prepare_for_eval_eps_grav_and_partials', ierr)
                  if (failed('interp_pm lnS')) return                       
                  call interp_values( &
                     s% prev_mesh_xm, prev_nz, f1, s% k_const_mass, xm, s% prev_lnS, ierr)
                  if (failed('interp_values lnS')) return              
               end if
               
               do k=s% k_const_mass,nz
                  s% prev_lnS(k) = s% lnS_pre(k)
               end do
               
            end if
            
            if (s% k_Lagrangian_eps_grav > 1) then ! some cells use Eulerian
               ! do prev_lnS_const_q
               !write(*,*) 'prev_lnS_const_q'
            
               call non_crit_get_work_array( &
                  s, xq, nz, 0, 'prepare_for_eval_eps_grav_and_partials', ierr)
               if (failed('xq')) return
            
               call non_crit_get_work_array( &
                  s, prev_xq, prev_nz, 0, 'prepare_for_eval_eps_grav_and_partials', ierr)
               if (failed('prev_xq')) return
            
               do k=1,prev_nz
                  prev_xq(k) = s% prev_mesh_xm(k)/s% xmstar_old
               end do
               xq(1) = 0
               do k=2,nz
                  xq(k) = xq(k-1) + s% dq(k-1)
               end do
            
               call non_crit_do1_alloc_if_necessary( &
                  s, s% prev_lnS_const_q, nz, 'prepare_for_eval_eps_grav_and_partials', ierr)
               if (failed('prev_lnS_const_q')) return
               
               call interp_pm(prev_xq, prev_nz, f1, &
                  nwork, work, 'prepare_for_eval_eps_grav_and_partials', ierr)
               if (failed('interp_pm lnS_const_q')) return 
                         
               call interp_values( &
                  prev_xq, prev_nz, f1, nz, xq, s% prev_lnS_const_q, ierr)
               if (failed('interp_values lnS')) return
            
            end if
            
            s% have_prev_lnS = .true.      
               
         end if
         
         call dealloc
         
         
         contains
         
         
         subroutine dealloc
            if (associated(xm)) call non_crit_return_work_array( &
               s, xm, 'prepare_for_eval_eps_grav_and_partials')
            if (associated(work)) call non_crit_return_work_array( &
               s, work, 'prepare_for_eval_eps_grav_and_partials')
            if (associated(f1)) call non_crit_return_work_array( &
               s, f1, 'prepare_for_eval_eps_grav_and_partials')
            if (associated(xq)) call non_crit_return_work_array( &
               s, xq, 'prepare_for_eval_eps_grav_and_partials')
            if (associated(prev_xq)) call non_crit_return_work_array( &
               s, prev_xq, 'prepare_for_eval_eps_grav_and_partials')
         end subroutine dealloc
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            failed = .true.
            call dealloc
            if (s% report_ierr) &
               write(*,'(a)') 'prepare_for_eval_eps_grav_and_partials: ' // trim(str)
         end function failed

         
      end subroutine prepare_for_eval_eps_grav_and_partials
      
      
      subroutine eval_eps_grav_and_partials(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: f
         include 'formats'
         ierr = 0
         
         if (s% dt <= 0) then
            call zero_eps_grav_and_partials(s, k)
            return
         end if

         call eval1_eps_grav_and_partials(s, k, ierr)
         if (ierr /= 0) return
         
         if (s% use_other_eps_grav) then
            ! note: call this after 1st doing the standard calculation
            call s% other_eps_grav(s% id, k, s% dt, ierr)
            if (ierr /= 0) return
         end if
         
         f = s% eps_grav_factor
         if (abs(f - 1d0) < 1d-12) return
         
         s% eps_grav(k) = f*s% eps_grav(k)
         s% d_eps_grav_dlndm1(k) = f*s% d_eps_grav_dlndm1(k)
         s% d_eps_grav_dlnd00(k) = f*s% d_eps_grav_dlnd00(k)
         s% d_eps_grav_dlndp1(k) = f*s% d_eps_grav_dlndp1(k)
         s% d_eps_grav_dlnTm1(k) = f*s% d_eps_grav_dlnTm1(k)
         s% d_eps_grav_dlnT00(k) = f*s% d_eps_grav_dlnT00(k)
         s% d_eps_grav_dlnTp1(k) = f*s% d_eps_grav_dlnTp1(k)
         s% d_eps_grav_dlnR00(k) = f*s% d_eps_grav_dlnR00(k)
         s% d_eps_grav_dL00(k) = f*s% d_eps_grav_dL00(k)
         s% d_eps_grav_dlnPgas00_const_T(k) = f*s% d_eps_grav_dlnPgas00_const_T(k)
         s% d_eps_grav_dlnPgasm1_const_T(k) = f*s% d_eps_grav_dlnPgasm1_const_T(k)
         s% d_eps_grav_dlnPgasp1_const_T(k) = f*s% d_eps_grav_dlnPgasp1_const_T(k)
         s% d_eps_grav_dlnTm1_const_Pgas(k) = f*s% d_eps_grav_dlnTm1_const_Pgas(k)
         s% d_eps_grav_dlnT00_const_Pgas(k) = f*s% d_eps_grav_dlnT00_const_Pgas(k)
         s% d_eps_grav_dlnTp1_const_Pgas(k) = f*s% d_eps_grav_dlnTp1_const_Pgas(k)
         s% d_eps_grav_dlnRp1(k) = f*s% d_eps_grav_dlnRp1(k)
         s% d_eps_grav_dv00(k) = f*s% d_eps_grav_dv00(k)
         s% d_eps_grav_dvp1(k) = f*s% d_eps_grav_dvp1(k)
         
      end subroutine eval_eps_grav_and_partials
      
   
      subroutine eval1_eps_grav_and_partials(s, k, ierr)
         use eos_def
         use eos_lib
         use mlt_def
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         include 'formats'         
         ierr = 0         
            
         if (k < s% k_below_just_added .and. &
               s% zero_eps_grav_in_just_added_material) then
            call zero_eps_grav_and_partials(s, k)
            return
         end if
                  
         if (s% use_lnS_for_eps_grav) then         
            call do_eps_grav_with_lnS(s, k, ierr)            
         else if (s% use_dEdRho_form_for_eps_grav) then         
            call do_dEdRho_form_for_eps_grav(s, k, ierr)            
         else if (s% use_PdVdt_form_for_eps_grav) then         
            call do_PdVdt_form_for_eps_grav(s, k, ierr)            
         else if (s% use_dlnd_dt_form_for_eps_grav) then         
            call do_dlnd_dt_form_for_eps_grav(s, k, ierr)            
         else ! check to see if need to use lnS form         
            if (s% gam_start(k) >= s% Gamma_lnS_eps_grav_full_on) then
               call do_eps_grav_with_lnS(s, k, ierr)
            else if (s% gam_start(k) > s% Gamma_lnS_eps_grav_full_off) then
               call blend_with_lnS_form(s, k, ierr)
            else if (s% lnPgas_flag) then
               call do_eps_grav_with_lnPgas(s, k, ierr)
            else
               call do_eps_grav_with_lnd(s, k, ierr)
            end if            
         end if
         
         if (.false. .and. k == s% nz) &
               write(*,3) 'dlnd/dt, dlnT/dt, eps_grav', &
                  k, s% model_number, s% dlnd_dt(k), s% dlnT_dt(k), s% eps_grav(k)

         if (ierr /= 0 .or. is_bad_num(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) &
                  'failed in eval_eps_grav_and_partials', k, s% eps_grav(k)
            end if
            return
         end if  

      end subroutine eval1_eps_grav_and_partials
      
      
      real(dp) function get_Eulerian_fraction_for_eps_grav(s, k, ierr) result(alfa)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: q_Eulerian, q_Lagrangian
         include 'formats'
         ! alfa is fraction of Eulerian form in result
         
         ierr = 0

         if (k < s% k_below_Eulerian_eps_grav .or. &
               s% k_below_Eulerian_eps_grav > s% nz) then
            alfa = 1d0 ! pure Eulerian
         else if (k >= s% k_Lagrangian_eps_grav) then
            alfa = 0d0 ! pure Lagrangian
         else
            if (s% k_below_Eulerian_eps_grav == 1) then
               q_Eulerian = 1d0
            else
               q_Eulerian = s% q(s% k_below_Eulerian_eps_grav-1)
            end if
            if (s% k_Lagrangian_eps_grav > s% nz) then
               q_Lagrangian = 0d0
            else
               q_Lagrangian = s% q(s% k_Lagrangian_eps_grav)
            end if
            alfa = max(0d0, min(1d0, (s% q(k) - q_Lagrangian)/(q_Eulerian - q_Lagrangian)))
            if (is_bad_num(alfa)) then
               write(*,2) 's% k_below_Eulerian_eps_grav', s% k_below_Eulerian_eps_grav, q_Eulerian
               write(*,2) 's% q(k)', k, s% q(k)
               write(*,2) 's% k_Lagrangian_eps_grav', s% k_Lagrangian_eps_grav, q_Lagrangian
               write(*,2) 'alfa', k, alfa 
               write(*,*) 'failed in get_Eulerian_fraction_for_eps_grav'
               stop 1
            end if   
         end if
         
         if (k == s% trace_k) then
            write(*,5) 'Eulerian_fraction_for_eps_grav', &
               k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, alfa
         end if
                        
      end function get_Eulerian_fraction_for_eps_grav


      subroutine do_dlnd_dt_form_for_eps_grav(s, k, ierr)
         use eos_def, only: i_lnE
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, dedt, P_div_rho_dlnddt, P, rho
         
         include 'formats'
         
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         if (alfa /= 0d0) then
            stop 'no support for Eulerian form of do_dlnd_dt_form_for_eps_grav'
         end if
                  
         dedt = (s% energy(k) - s% energy_start(k))*s% dVARDOT_dVAR
         
         rho = 0.5*(s% rho(k) + s% rho_start(k)) ! time centering to match PdvAdm
         P = s% P(k) ! no time centering for P
         P_div_rho_dlnddt = P/rho*s% dlnd_dt(k)
         
         ! eps_grav := -dedt + P/rho*dlnd_dt, with time centering for rho
         s% eps_grav(k) = -dedt + P_div_rho_dlnddt
         
         s% d_eps_grav_dlnd00(k) = &
            -s% energy(k)*s% d_eos_dlnd(i_lnE,k)*s% dVARDOT_dVAR &
            + s% chiRho(k)*P_div_rho_dlnddt &
            - 0.5d0*s% rho(k)*P/(rho*rho)*s% dlnd_dt(k) &
            + P/rho*s% dVARDOT_dVAR
            
         s% d_eps_grav_dlnT00(k) = &
            - s% energy(k)*s% d_eos_dlnT(i_lnE,k)*s% dVARDOT_dVAR &
            + s% chiT(k)*P_div_rho_dlnddt

         call include_dmu_dt_in_eps_grav(s, k)
         
      end subroutine do_dlnd_dt_form_for_eps_grav


      subroutine do_dEdRho_form_for_eps_grav(s, k, ierr)
         use eos_def, only: i_Cv, i_dE_dRho
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, eps_grav1, d_dlnd1, d_dlnT1, &
            rho, cvT, rho_dE_dRho, P, P_div_rho, dlnT_dt, dlnd_dt, &
            d_cvT_dlnd, d_cvT_dlnT, d_rho_dE_dRho_dlnd, d_rho_dE_dRho_dlnT, &
            d_P_div_rho_dlnd, d_P_div_rho_dlnT
         
         include 'formats'
         
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         if (alfa /= 0d0) then
            stop 'no support for Eulerian form of do_dEdRho_form_for_eps_grav'
         end if
         
         call do_dlnd_dt_form_for_eps_grav(s, k, ierr)
         eps_grav1 = s% eps_grav(k)
         d_dlnd1 = s% d_eps_grav_dlnd00(k)
         d_dlnT1 = s% d_eps_grav_dlnT00(k)
         
         rho = 0.5d0*(s% rho(k) + s% rho_start(k))
         cvT = 0.5d0*(s% cv(k)*s% T(k) + s% cv_start(k)*s% T_start(k))
         rho_dE_dRho = 0.5d0*(s% rho(k)*s% dE_dRho(k) + &
                              s% rho_start(k)*s% dE_dRho_start(k))
         P = s% P(k) ! no time centering for P in momentum equation or here
         P_div_rho = P/rho
         
         dlnT_dt = s% dlnT_dt(k)
         dlnd_dt = s% dlnd_dt(k)
         
         ! eps_grav := -(cv*T*dlnT_dt + (rho*dE_dRho - P/rho)*dlnd_dt)
         s% eps_grav = -(cvT*dlnT_dt + (rho_dE_dRho - P_div_rho)*dlnd_dt)
         
         ! eps_grav = -cvT*dlnT_dt - rho_dE_dRho*dlnd_dt + P_div_rho*dlnd_dt
         
         d_cvT_dlnd = 0.5d0*s% d_eos_dlnd(i_Cv,k)*s% T(k)
         d_cvT_dlnT = 0.5d0*(s% d_eos_dlnT(i_Cv,k) + s% cv(k))*s% T(k)
         
         d_rho_dE_dRho_dlnd = &
            0.5d0*s% rho(k)*s% dE_dRho(k)*dlnd_dt + &
            0.5d0*s% rho(k)*s% d_eos_dlnd(i_dE_dRho,k)*dlnd_dt + &
            0.5d0*s% rho(k)*s% dE_dRho(k)*s% dVARDOT_dVAR
         d_rho_dE_dRho_dlnT = &
            0.5d0*s% rho(k)*s% d_eos_dlnT(i_dE_dRho,k)*dlnd_dt
            
         d_P_div_rho_dlnd = &
            s% chiRho(k)*P_div_rho - P_div_rho*0.5d0*s% rho(k)/rho
         d_P_div_rho_dlnT = s% chiT(k)*P_div_rho
         
         s% d_eps_grav_dlnd00(k) = &     
            - d_cvT_dlnd*dlnT_dt &  
            - d_rho_dE_dRho_dlnd*dlnd_dt &
            + d_P_div_rho_dlnd*dlnd_dt
            
         s% d_eps_grav_dlnT00(k) = &
            - d_cvT_dlnT*dlnT_dt &  
            - d_rho_dE_dRho_dlnT*dlnd_dt &
            + d_P_div_rho_dlnT*dlnd_dt &
            - cvT*s% dVARDOT_dVAR
            
         !if (abs(s% eps_grav(k) - eps_grav1) > 1d-5*abs(eps_grav1)) then
         if (k == 167 .and. abs(s% eps_grav(k)) > 0) then
!$omp critical
            write(*,2) 's% eps_grav(k) - eps_grav1', k, s% eps_grav(k) - eps_grav1
            write(*,2) 's% eps_grav(k)', k, s% eps_grav(k)
            write(*,2) 'eps_grav1', k, eps_grav1
            write(*,*)
            stop 'do_dEdRho_form_for_eps_grav'
!$omp end critical
         end if
            
         call include_dmu_dt_in_eps_grav(s, k)
         
      end subroutine do_dEdRho_form_for_eps_grav


      subroutine do_PdVdt_form_for_eps_grav(s, k, ierr)
         use eos_def, only: i_lnE
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, dedt, P_dVdt
         
         include 'formats'
         
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         if (alfa /= 0d0) then
            stop 'no support for Eulerian form of do_PdVdt_form_for_eps_grav'
         end if
         
         ! eps_grav := -dedt + P*d(1/rho)/dt
         
         dedt = (s% energy(k) - s% energy_start(k))*s% dVARDOT_dVAR
         P_dVdt = s% P(k)*(1/s% rho(k) - 1/s% rho_start(k))*s% dVARDOT_dVAR
         
         s% eps_grav(k) = -(dedt + P_dVdt)

         s% d_eps_grav_dlnd00(k) = &
            -s% energy(k)*s% d_eos_dlnd(i_lnE,k)*s% dVARDOT_dVAR &
            + s% chiRho(k)*P_dVdt &
            - s% P(k)/s% rho(k)*s% dVARDOT_dVAR
            
         s% d_eps_grav_dlnT00(k) = &
            - s% energy(k)*s% d_eos_dlnT(i_lnE,k)*s% dVARDOT_dVAR &
            + s% chiT(k)*P_dVdt

         call include_dmu_dt_in_eps_grav(s, k)
         
      end subroutine do_PdVdt_form_for_eps_grav


      subroutine do_eps_grav_with_lnd(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa
         include 'formats'
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         call combine_two_eps_gravs( &
            s, k, alfa, 1d0 - alfa, &
            do_eps_grav_with_lnd_Eulerian, &
            do_eps_grav_with_lnd_Lagrangian, ierr)
         if (ierr /= 0) return
         call include_dmu_dt_in_eps_grav(s, k)
      end subroutine do_eps_grav_with_lnd
      

      subroutine do_eps_grav_with_lnd_Lagrangian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call do_lnd_eps_grav(s, k, s% dlnd_dt(k), s% dlnT_dt(k), ierr)
      end subroutine do_eps_grav_with_lnd_Lagrangian


      subroutine do_eps_grav_with_lnd_Eulerian(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         include 'formats'
         call combine_two_eps_gravs( &
            s, k, 1d0, 1d0, do_spatial_term, do_dt_at_const_q_with_lnd, ierr)
      end subroutine do_eps_grav_with_lnd_Eulerian


      subroutine do_dt_at_const_q_with_lnd(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call do_lnd_eps_grav( &
            s, k, s% dlnd_dt_const_q(k), s% dlnT_dt_const_q(k), ierr)
         s% eps_grav_dt_term_const_q(k) = s% eps_grav(k)
      end subroutine do_dt_at_const_q_with_lnd
      
      
      ! this uses the given args to calculate -T*ds/dt
      subroutine do_lnd_eps_grav(s, k, dlnd_dt, dlnT_dt, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(in) :: dlnd_dt, dlnT_dt
         integer, intent(out) :: ierr

         real(dp) :: dT1_dlnTdot, a1, da1_dlnd, da1_dlnT, &
            T1, dT1_dlnd, dT1_dlnT00, dT1_dlnd_dt, dT1_d_dlnTdt, &
            a2, da2_dlnd, da2_dlnT, &
            T2, dT2_dlnT, dT2_dlnd00, &
            T3, dT3_dlnd, dT3_dlnT

         include 'formats'
         ierr = 0
         
         call zero_eps_grav_and_partials(s, k)
      
         !s% eps_grav(k) = -s% T(k)*s% cp(k)* &
         !      ((1-s% grada(k)*s% chiT(k))*dlnT_dt &
         !        - s% grada(k)*s% chiRho(k)*dlnd_dt)

         a1 = 1 - s% grada(k)*s% chiT(k)
         da1_dlnd = -(s% d_eos_dlnd(i_grad_ad,k)*s% chiT(k) + s% grada(k)*s% d_eos_dlnd(i_chiT,k))
         da1_dlnT = -(s% d_eos_dlnT(i_grad_ad,k)*s% chiT(k) + s% grada(k)*s% d_eos_dlnT(i_chiT,k))

         T1 = dlnT_dt*a1
         dT1_dlnd = dlnT_dt*da1_dlnd

         dT1_d_dlnTdt = a1
         dT1_dlnT00 = s% dVARDOT_dVAR*a1 + dlnT_dt*da1_dlnT
         
         a2 = s% grada(k)*s% chiRho(k)
         da2_dlnd = s% d_eos_dlnd(i_grad_ad,k)*s% chiRho(k) + s% grada(k)*s% d_eos_dlnd(i_chiRho,k)
         da2_dlnT = s% d_eos_dlnT(i_grad_ad,k)*s% chiRho(k) + s% grada(k)*s% d_eos_dlnT(i_chiRho,k)

         T2 = dlnd_dt*a2
         dT2_dlnT = dlnd_dt*da2_dlnT
         
         dT2_dlnd00 = s% dVARDOT_dVAR*a2 + dlnd_dt*da2_dlnd

         T3 = -s% T(k)*s% cp(k)
         dT3_dlnd = -s% T(k)*s% d_eos_dlnd(i_Cp,k)
         dT3_dlnT = -s% T(k)*(s% cp(k) + s% d_eos_dlnT(i_Cp,k))

         ! eps_grav = T3*(T1-T2)
         s% eps_grav(k) = T3*(T1-T2)
      
         if (is_bad_num(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'do_lnd_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
            return
         end if
      
         s% d_eps_grav_dlndm1(k) = 0
         s% d_eps_grav_dlndp1(k) = 0
         s% d_eps_grav_dlnd00(k) = (T3*(dT1_dlnd - dT2_dlnd00) + dT3_dlnd*(T1-T2))
         
         s% d_eps_grav_dlnTm1(k) = 0
         s% d_eps_grav_dlnTp1(k) = 0
         s% d_eps_grav_dlnT00(k) = (T3*(dT1_dlnT00 - dT2_dlnT) + dT3_dlnT*(T1-T2))
         
         if (k == s% trace_k) then
            write(*,5) 'do_lnd_eps_grav', &
               k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,2) 's% rho(k)', k, s% rho(k)
            write(*,2) 'dlnd_dt', k, dlnd_dt
            write(*,2) 'dlnT_dt', k, dlnT_dt
            write(*,2) 's% cp(k)', k, s% cp(k)
            write(*,2) 's% grada(k)', k, s% grada(k)
            write(*,2) 's% chiT(k)', k, s% chiT(k)
            write(*,2) 'a1', k, a1
            write(*,2) 'da1_dlnd', k, da1_dlnd
            write(*,2) 'da1_dlnT', k, da1_dlnT
            write(*,2) 's% chiRho(k)', k, s% chiRho(k)
            write(*,2) 'a2', k, a2
            write(*,2) 'da2_dlnd', k, da2_dlnd
            write(*,2) 'da2_dlnT', k, da2_dlnT
            write(*,2) 'T1', k, T1
            write(*,2) 'T2', k, T2
            write(*,2) 'T3', k, T3
            write(*,2) 'dT1_dlnd', k, dT1_dlnd
            write(*,2) 'dT2_dlnd00', k, dT2_dlnd00
            write(*,2) 'dT3_dlnd', k, dT3_dlnd
            write(*,2) 'dT1_dlnT00', k, dT1_dlnT00
            write(*,2) 'dT2_dlnT', k, dT2_dlnT
            write(*,2) 'dT3_dlnT', k, dT3_dlnT
         end if
         
      end subroutine do_lnd_eps_grav


      subroutine do_eps_grav_with_lnPgas(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         call combine_two_eps_gravs( &
            s, k, alfa, 1d0 - alfa, &
            do_eps_grav_with_lnPgas_Eulerian, &
            do_eps_grav_with_lnPgas_Lagrangian, ierr)
         if (ierr /= 0) return
         call include_dmu_dt_in_eps_grav(s, k)
      end subroutine do_eps_grav_with_lnPgas


      subroutine do_eps_grav_with_lnPgas_Lagrangian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call do_lnPgas_eps_grav(s, k, s% dlnPgas_dt(k), s% dlnT_dt(k), ierr)
      end subroutine do_eps_grav_with_lnPgas_Lagrangian


      subroutine do_eps_grav_with_lnPgas_Eulerian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call combine_two_eps_gravs( &
            s, k, 1d0, 1d0, do_spatial_term, do_dt_at_const_q_with_lnPgas, ierr)
      end subroutine do_eps_grav_with_lnPgas_Eulerian


      subroutine do_dt_at_const_q_with_lnPgas(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call do_lnPgas_eps_grav( &
            s, k, s% dlnPgas_dt_const_q(k), s% dlnT_dt_const_q(k), ierr)
         s% eps_grav_dt_term_const_q(k) = s% eps_grav(k)
      end subroutine do_dt_at_const_q_with_lnPgas
      

      ! this uses the given args to calculate -T*ds/dt
      ! for Lagrangian eps_grav, call it with constant mass time derivatives.
      ! for Eulerian non-spatial_term term, call it with constant q time derivatives.
      subroutine do_lnPgas_eps_grav(s, k, dlnPgas_dt, dlnT_dt, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(in) :: dlnPgas_dt, dlnT_dt
         integer, intent(out) :: ierr
      
         !s% eps_grav(k) = -s% T(k)*s% cp(k)*(dlnT_dt - s% grada(k)*dlnP_dt)
         ! dlnP_dt = (dPgas_dt + dPrad_dt)/P
         !         = (Pgas*dlnPgas_dt + 4*Prad*dlnT_dt)/P
         !s% eps_grav(k) = -s% T(k)*s% cp(k)* &
         !     ((1 - s% grada(k)*4*Prad/P)*dlnT_dt - s% grada(k)*Pgas/P*dlnPgas_dt)
         
         real(dp) :: &
            P2, a1, da1_dlnd, da1_dlnT, &
            T1, dT1_dlnd, dT1_dlnT00, dT1_dlnd_dt, dT1_d_dlnTdt, &
            a2, da2_dlnd, da2_dlnT, &
            T2, dT2_dlnT, dT2_dlnd00, &
            T3, dT3_dlnd, dT3_dlnT, &
            dlnd_dlnPgas00, dgrada_dlnPgas, &
            dT1_dlnPgas, dT2_dlnPgas00, dT3_dlnPgas, &
            d_eps_grav_dlnPgas00, &
            dP_dlnPgas, dPinv_dlnPgas, dPrad_dlnT, &
            da1_dlnPgas, da2_dlnPgas, dP_dlnT_const_Pgas, dPinv_dlnT_const_Pgas, &
            da1_dlnPgas_const_T, da1_dlnT_const_Pgas, da2_dlnPgas_const_T, &
            da2_dlnT_const_Pgas, dcp_dlnPgas_const_T, dcp_dlnT_const_Pgas, &
            dgrada_dlnPgas_const_T, dgrada_dlnT_const_Pgas, dPinv_dlnPgas_const_T, &
            dPgas_dlnPgas_const_T, dPrad_dlnT_const_Pgas, dP_dlnPgas_const_T, &
            dT1_dlnPgas_const_T, dT1_dlnT00_const_Pgas, &
            dT2_dlnPgas00_const_T, &
            dT2_dlnT_const_Pgas, &
            dT3_dlnT_const_Pgas, dT3_dlnPgas_const_T

         include 'formats'
         ierr = 0
         call zero_eps_grav_and_partials(s, k)
         
         dgrada_dlnPgas_const_T = s% d_eos_dlnd(i_grad_ad,k)*s% dlnRho_dlnPgas_const_T(k)
         dgrada_dlnT_const_Pgas = &
            s% d_eos_dlnT(i_grad_ad,k) + s% d_eos_dlnd(i_grad_ad,k)*s% dlnRho_dlnT_const_Pgas(k)
         
         dgrada_dlnPgas_const_T = s% d_eos_dlnd(i_grad_ad,k)*s% dlnRho_dlnPgas_const_T(k)
         dgrada_dlnT_const_Pgas = &
            s% d_eos_dlnT(i_grad_ad,k) + s% d_eos_dlnd(i_grad_ad,k)*s% dlnRho_dlnT_const_Pgas(k)

         dcp_dlnPgas_const_T = s% d_eos_dlnd(i_cp,k)*s% dlnRho_dlnPgas_const_T(k)
         dcp_dlnT_const_Pgas = &
            s% d_eos_dlnT(i_cp,k) + s% d_eos_dlnd(i_cp,k)*s% dlnRho_dlnT_const_Pgas(k)
         
         dPgas_dlnPgas_const_T = s% Pgas(k)
         dP_dlnPgas_const_T = dPgas_dlnPgas_const_T
         P2 = s% P(k)*s% P(k)
         dPinv_dlnPgas_const_T = -dP_dlnPgas_const_T/P2
         
         dPrad_dlnT_const_Pgas = 4*s% Prad(k)
         dP_dlnT_const_Pgas = dPrad_dlnT_const_Pgas
         dPinv_dlnT_const_Pgas = -dP_dlnT_const_Pgas/P2
         
         a1 = 1 - s% grada(k)*4*s% Prad(k)/s% P(k)
         da1_dlnPgas_const_T = &
            - dgrada_dlnPgas_const_T*4*s% Prad(k)/s% P(k) &
            - s% grada(k)*4*s% Prad(k)*dPinv_dlnPgas_const_T
         da1_dlnT_const_Pgas = &
            - dgrada_dlnT_const_Pgas*4*s% Prad(k)/s% P(k) &
            - s% grada(k)*4*dPrad_dlnT_const_Pgas/s% P(k) &
            - s% grada(k)*4*s% Prad(k)*dPinv_dlnT_const_Pgas
         
         T1 = dlnT_dt*a1
         dT1_dlnPgas_const_T = dlnT_dt*da1_dlnPgas_const_T
         dT1_dlnT00_const_Pgas = s% dVARDOT_dVAR*a1 + dlnT_dt*da1_dlnT_const_Pgas
         
         a2 = s% grada(k)*s% Pgas(k)/s% P(k)
         da2_dlnPgas_const_T = &
            dgrada_dlnPgas_const_T*s% Pgas(k)/s% P(k) + &
            s% grada(k)*dPgas_dlnPgas_const_T/s% P(k) + &
            s% grada(k)*s% Pgas(k)*dPinv_dlnPgas_const_T
         da2_dlnT_const_Pgas = &
            dgrada_dlnT_const_Pgas*s% Pgas(k)/s% P(k) + &
            s% grada(k)*s% Pgas(k)*dPinv_dlnT_const_Pgas
         
         T2 = dlnPgas_dt*a2
         dT2_dlnT_const_Pgas = dlnPgas_dt*da2_dlnT_const_Pgas
         dT2_dlnPgas00_const_T = s% dVARDOT_dVAR*a2 + dlnPgas_dt*da2_dlnPgas_const_T
         
         T3 = -s% T(k)*s% cp(k)
         dT3_dlnPgas_const_T = -s% T(k)*dcp_dlnPgas_const_T
         dT3_dlnT_const_Pgas =  -(s% T(k)*s% cp(k) + s% T(k)*dcp_dlnT_const_Pgas)

         s% eps_grav(k) = T3*(T1-T2)
      
         if (is_bad_num(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'do_lnPgas_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
            return
         end if

         s% d_eps_grav_dlnPgasm1_const_T(k) = 0
         s% d_eps_grav_dlnPgasp1_const_T(k) = 0
         s% d_eps_grav_dlnPgas00_const_T(k) = &
            (T3*(dT1_dlnPgas_const_T - dT2_dlnPgas00_const_T) + dT3_dlnPgas_const_T*(T1-T2))
         
         s% d_eps_grav_dlnTm1_const_Pgas(k) = 0
         s% d_eps_grav_dlnTp1_const_Pgas(k) = 0
         s% d_eps_grav_dlnT00_const_Pgas(k) = &
            (T3*(dT1_dlnT00_const_Pgas - dT2_dlnT_const_Pgas) + dT3_dlnT_const_Pgas*(T1-T2))
         
         if (k == s% trace_k) then
            write(*,5) 'do_lnPgas_eps_grav', &
               k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,2) 's% grada(k)', k, s% grada(k)
            write(*,2) 's% Prad(k)', k, s% Prad(k)
            write(*,2) 's% P(k)', k, s% P(k)
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,2) 's% cp(k)', k, s% cp(k)
            write(*,2) 'dlnT_dt', k, dlnT_dt
            write(*,2) 'dlnPgas_dt', k, dlnPgas_dt
         end if

      end subroutine do_lnPgas_eps_grav


      subroutine do_eps_grav_with_lnS(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa
         include 'formats'
         alfa = get_Eulerian_fraction_for_eps_grav(s, k, ierr)
         if (ierr /= 0) return
         call combine_two_eps_gravs( &
            s, k, alfa, 1d0 - alfa, &
            do_eps_grav_with_lnS_Eulerian, &
            do_eps_grav_with_lnS_Lagrangian, ierr)
         if (ierr /= 0) return
         call include_dmu_dt_in_eps_grav(s, k)
      end subroutine do_eps_grav_with_lnS


      subroutine do_eps_grav_with_lnS_Lagrangian(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad, i_chiRho, i_chiT
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         include 'formats'
         if (.not. associated(s% prev_lnS)) &
            stop 'prev_lnS missing do_eps_grav_with_lnS_Lagrangian'
         call do_lnS_eps_grav(s, k, s% prev_lnS(k), ierr)
      end subroutine do_eps_grav_with_lnS_Lagrangian


      subroutine do_eps_grav_with_lnS_Eulerian(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         call combine_two_eps_gravs( &
            s, k, 1d0, 1d0, do_spatial_term, do_dt_at_const_q_with_lnS, ierr)
      end subroutine do_eps_grav_with_lnS_Eulerian


      subroutine do_dt_at_const_q_with_lnS(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         if (.not. associated(s% prev_lnS_const_q)) &
            stop 'prev_lnS_const_q missing do_eps_grav_with_lnS_Lagrangian'
         call do_lnS_eps_grav(s, k, s% prev_lnS_const_q(k), ierr)
         s% eps_grav_dt_term_const_q(k) = s% eps_grav(k)
      end subroutine do_dt_at_const_q_with_lnS


      subroutine do_lnS_eps_grav(s, k, prev_lnS, ierr)
         use eos_def, only: i_lnS
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(in) :: prev_lnS
         integer, intent(out) :: ierr
         
         real(dp) :: entropy, S_prev, T, dS_dlnT, dS_dlnd
         
         include 'formats'
         ierr = 0
         call zero_eps_grav_and_partials(s, k)

         entropy = exp_cr(s% lnS(k))
         T = s% T(k)
         S_prev = exp_cr(prev_lnS)
                     
         s% eps_grav(k) = -T*(entropy - S_prev)*s% dVARDOT_dVAR
      
         if (is_bad_num(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'do_lnS_eps_grav -- bad value for eps_grav', k, s% eps_grav(k)
            return
         end if
         
         dS_dlnT = entropy*s% d_eos_dlnT(i_lnS,k)
         dS_dlnd = entropy*s% d_eos_dlnd(i_lnS,k)
         s% d_eps_grav_dlnT00(k) = -T*dS_dlnT*s% dVARDOT_dVAR + s% eps_grav(k)
         s% d_eps_grav_dlnd00(k) = -T*dS_dlnd*s% dVARDOT_dVAR

         if (s% lnPgas_flag) then
            s% d_eps_grav_dlnPgas00_const_T(k) = &
               s% d_eps_grav_dlnd00(k)*s% dlnRho_dlnPgas_const_T(k)
            s% d_eps_grav_dlnT00_const_Pgas(k) = &
               s% d_eps_grav_dlnT00(k) + &
               s% d_eps_grav_dlnd00(k)*s% dlnRho_dlnT_const_Pgas(k)
         end if
         
         if (k == s% trace_k) then
            write(*,5) 'do_lnS_eps_grav', &
               k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,2) 'entropy', k, entropy
            write(*,2) 'T', k, T
            write(*,2) 'S_prev', k, S_prev
         end if
            
      end subroutine do_lnS_eps_grav
      
      
      subroutine do_spatial_term(s, k, ierr)
         use eos_def, only: i_Cp, i_grad_ad
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
      
         ! Townsley & Bildsten, The Astrophysical Journal, 600:390–403, 2004 January 1
         ! ds/dm = (ds/dP)*(dP/dm)
         ! dP/dm = -Gm/(4 pi r^4) 
         ! ds/dP = (gradT - grada)*Cp/P
         ! gradT = (dlnT/dlnP)_actual
         ! grada = (dlnT/dlnP)_adiabatic
         ! in mesa/star m = M_center + q*xmstar,
         ! so (dm/dt)_q = q*Mdot since M_center is constant during timesteps
         ! and d(xmstar)/dt = Mdot.
         ! Putting it all together,
         ! we get eps_grav_h = (-G m q Mdot Cp T)*(gradT - grada)/(4 pi r^4 P)
         
         real(dp) :: &
            P2, rmid3, rmid4, inv_rmid4, inv_4pi_r4, P, inv_4pi_r4_P, &
            gradT_mid, dg, mmid, qmid, G_m_q_mdot, CpT, &
            d_inv_rmid4_dlnR00, d_inv_rmid4_dlnRp1, &
            d_inv_4pi_r4_P_dlnR00, d_inv_4pi_r4_P_dlnRp1, &
            d_dg_dlnR00, d_dg_dlnRp1, d_dg_dL00, d_dg_dLp1, &
            dgradT_mid_dlndm1, dgradT_mid_dlnd00, dgradT_mid_dlndp1, &
            d_dg_dlndm1, d_dg_dlnd00, d_dg_dlndp1, &
            dgradT_mid_dlnTm1, dgradT_mid_dlnT00, dgradT_mid_dlnTp1, &
            d_dg_dlnTm1, d_dg_dlnT00, d_dg_dlnTp1, d_CpT_dlnd, d_CpT_dlnT, &
            d_inv_4pi_r4_P_dlnd, d_inv_4pi_r4_P_dlnT, &
            dP_dlnPgas_const_T, dP_dlnT_const_Pgas, &
            dPinv_dlnPgas_const_T, dPinv_dlnT_const_Pgas, &
            d_inv_4pi_r4_P_dlnPgas_const_T, d_inv_4pi_r4_P_dlnT_const_Pgas, &
            d_CpT_dlnPgas_const_T, d_CpT_dlnT_const_Pgas, &
            d_dg_dlnPgas00_const_T, d_dg_dlnT00_const_Pgas, &
            d_dg_dlnPgasp1_const_T, d_dg_dlnTp1_const_Pgas, &
            d_dg_dlnPgasm1_const_T, d_dg_dlnTm1_const_Pgas

         include 'formats'
         ierr = 0
         
         call zero_eps_grav_and_partials(s, k)
         s% eps_grav_dm_term_const_q(k) = 0
         s% eps_grav_dt_term_const_q(k) = 0
         if (k == s% nz) return
         
         rmid3 = s% rmid(k)*s% rmid(k)*s% rmid(k)
         rmid4 = rmid3*s% rmid(k)
         inv_rmid4 = 1d0/rmid4
         inv_4pi_r4 = inv_rmid4/(4*pi)
         P = s% P(k)
         inv_4pi_r4_P = inv_4pi_r4/P
         gradT_mid = 0.5d0*(s% gradT(k) + s% gradT(k+1))
         dg = s% grada(k) - gradT_mid
         mmid = 0.5d0*(s% m(k) + s% m(k+1))
         qmid = 0.5d0*(s% q(k) + s% q(k+1))
         G_m_q_mdot = s% cgrav(k)*mmid*qmid*s% mstar_dot
         CpT = s% Cp(k)*s% T(k)
         
         s% eps_grav(k) = G_m_q_mdot*CpT*dg*inv_4pi_r4_P
         s% eps_grav_dm_term_const_q(k) = s% eps_grav(k)
         
         if (is_bad_num(s% eps_grav(k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'do_spatial_term -- bad value for eps_grav', k, s% eps_grav(k)
            return
         end if
         
         d_inv_rmid4_dlnR00 = -2*inv_rmid4*s% r(k)*s% r(k)*s% r(k)/rmid3
         d_inv_rmid4_dlnRp1 = -2*inv_rmid4*s% r(k+1)*s% r(k+1)*s% r(k+1)/rmid3
         d_inv_4pi_r4_P_dlnR00 = d_inv_rmid4_dlnR00/(4*pi*P)
         d_inv_4pi_r4_P_dlnRp1 = d_inv_rmid4_dlnRp1/(4*pi*P)
         
         d_dg_dlnR00 = -0.5d0*s% d_gradT_dlnR(k)
         d_dg_dlnRp1 = -0.5d0*s% d_gradT_dlnR(k+1)

         s% d_eps_grav_dlnR00(k) = G_m_q_mdot*CpT* &
            (d_dg_dlnR00*inv_4pi_r4_P + dg*d_inv_4pi_r4_P_dlnR00)
         s% d_eps_grav_dlnRp1(k) = G_m_q_mdot*CpT* &
            (d_dg_dlnRp1*inv_4pi_r4_P + dg*d_inv_4pi_r4_P_dlnRp1)

         d_dg_dL00 = -0.5d0*s% d_gradT_dL(k)
         d_dg_dLp1 = -0.5d0*s% d_gradT_dL(k+1)
         
         s% d_eps_grav_dL00(k) = d_dg_dL00*G_m_q_mdot*CpT*inv_4pi_r4_P
         s% d_eps_grav_dLp1(k) = d_dg_dLp1*G_m_q_mdot*CpT*inv_4pi_r4_P
         
         dgradT_mid_dlndm1 = 0.5d0*s% d_gradT_dlndm1(k)
         dgradT_mid_dlnd00 = 0.5d0*(s% d_gradT_dlnd00(k) + s% d_gradT_dlndm1(k+1))
         dgradT_mid_dlndp1 = 0.5d0*s% d_gradT_dlnd00(k+1)
      
         d_dg_dlndm1 = -dgradT_mid_dlndm1
         d_dg_dlnd00 = s% d_eos_dlnd(i_grad_ad,k) - dgradT_mid_dlnd00
         d_dg_dlndp1 = -dgradT_mid_dlndp1
         
         dgradT_mid_dlnTm1 = 0.5d0*s% d_gradT_dlnTm1(k)
         dgradT_mid_dlnT00 = 0.5d0*(s% d_gradT_dlnT00(k) + s% d_gradT_dlnTm1(k+1))
         dgradT_mid_dlnTp1 = 0.5d0*s% d_gradT_dlnT00(k+1)

         d_dg_dlnTm1 = -dgradT_mid_dlnTm1
         d_dg_dlnT00 = s% d_eos_dlnT(i_grad_ad,k) - dgradT_mid_dlnT00
         d_dg_dlnTp1 = -dgradT_mid_dlnTp1

         d_CpT_dlnd = s% d_eos_dlnd(i_Cp,k)*s% T(k)
         d_CpT_dlnT = s% d_eos_dlnT(i_Cp,k)*s% T(k) + CpT 

         if (.not. s% lnPgas_flag) then

            d_inv_4pi_r4_P_dlnd = -s% chiRho(k)*inv_4pi_r4_P
            d_inv_4pi_r4_P_dlnT = -s% chiT(k)*inv_4pi_r4_P

            s% d_eps_grav_dlndm1(k) = d_dg_dlndm1*G_m_q_mdot*CpT*inv_4pi_r4_P
            s% d_eps_grav_dlndp1(k) = d_dg_dlndp1*G_m_q_mdot*CpT*inv_4pi_r4_P
            s% d_eps_grav_dlnd00(k) = G_m_q_mdot*( &
               d_dg_dlnd00*CpT*inv_4pi_r4_P + &
               dg*d_CpT_dlnd*inv_4pi_r4_P + &
               dg*CpT*d_inv_4pi_r4_P_dlnd)

            s% d_eps_grav_dlnTm1(k) = d_dg_dlnTm1*G_m_q_mdot*CpT*inv_4pi_r4_P
            s% d_eps_grav_dlnTp1(k) = d_dg_dlnTp1*G_m_q_mdot*CpT*inv_4pi_r4_P
            s% d_eps_grav_dlnT00(k) = G_m_q_mdot*( &
               d_dg_dlnT00*CpT*inv_4pi_r4_P + &
               dg*d_CpT_dlnT*inv_4pi_r4_P + &
               dg*CpT*d_inv_4pi_r4_P_dlnT)         
         
         else

            dP_dlnPgas_const_T = s% Pgas(k)
            dP_dlnT_const_Pgas = 4*s% Prad(k)
            
            P2 = s% P(k)*s% P(k)
            dPinv_dlnPgas_const_T = -dP_dlnPgas_const_T/P2
            dPinv_dlnT_const_Pgas = -dP_dlnT_const_Pgas/P2
            
            d_inv_4pi_r4_P_dlnPgas_const_T = inv_4pi_r4*dPinv_dlnPgas_const_T
            d_inv_4pi_r4_P_dlnT_const_Pgas = inv_4pi_r4*dPinv_dlnT_const_Pgas
            
            d_CpT_dlnPgas_const_T = d_CpT_dlnd*s% dlnRho_dlnPgas_const_T(k)
            d_CpT_dlnT_const_Pgas = &
               d_CpT_dlnT + d_CpT_dlnd*s% dlnRho_dlnT_const_Pgas(k)
            
            d_dg_dlnPgas00_const_T = d_dg_dlnd00*s% dlnRho_dlnPgas_const_T(k)
            d_dg_dlnT00_const_Pgas = &
               d_dg_dlnT00 + d_dg_dlnd00*s% dlnRho_dlnT_const_Pgas(k)
            
            d_dg_dlnPgasp1_const_T = d_dg_dlndp1*s% dlnRho_dlnPgas_const_T(k+1)
            d_dg_dlnTp1_const_Pgas = &
               d_dg_dlnTp1 + d_dg_dlndp1*s% dlnRho_dlnT_const_Pgas(k+1)
            
            if (k == 1) then
               d_dg_dlnPgasm1_const_T = 0
               d_dg_dlnTm1_const_Pgas = 0
            else
               d_dg_dlnPgasm1_const_T = &
                  d_dg_dlndm1*s% dlnRho_dlnPgas_const_T(k-1)
               d_dg_dlnTm1_const_Pgas = &
                  d_dg_dlnTm1 + d_dg_dlndm1*s% dlnRho_dlnT_const_Pgas(k-1)
            end if
            
            s% d_eps_grav_dlnPgasm1_const_T(k) = &
               d_dg_dlnPgasm1_const_T*G_m_q_mdot*CpT*inv_4pi_r4_P
            s% d_eps_grav_dlnPgas00_const_T(k) = G_m_q_mdot*( &
               d_dg_dlnPgas00_const_T*CpT*inv_4pi_r4_P + &
               dg*d_CpT_dlnPgas_const_T*inv_4pi_r4_P + &
               dg*CpT*d_inv_4pi_r4_P_dlnPgas_const_T)
            s% d_eps_grav_dlnPgasp1_const_T(k) = &
               d_dg_dlnPgasp1_const_T*G_m_q_mdot*CpT*inv_4pi_r4_P
            
            s% d_eps_grav_dlnTm1_const_Pgas(k) = &
               d_dg_dlnTm1_const_Pgas*G_m_q_mdot*CpT*inv_4pi_r4_P
            s% d_eps_grav_dlnT00_const_Pgas(k) = G_m_q_mdot*( &
               d_dg_dlnT00_const_Pgas*CpT*inv_4pi_r4_P + &
               dg*d_CpT_dlnT_const_Pgas*inv_4pi_r4_P + &
               dg*CpT*d_inv_4pi_r4_P_dlnT_const_Pgas)
            s% d_eps_grav_dlnTp1_const_Pgas(k) = &
               d_dg_dlnTp1_const_Pgas*G_m_q_mdot*CpT*inv_4pi_r4_P
         
         end if
         
         if (k == s% trace_k) then
            write(*,5) 'do_spatial_term eps_grav', &
               k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,*) 'gradT_mid', k, gradT_mid
            write(*,*) 's% gradT(k)', k, s% gradT(k)
            write(*,*) 's% gradT(k+1)', k+1, s% gradT(k+1)
            write(*,*) 'dg', k, dg
            write(*,*) 'G_m_q_mdot', k, G_m_q_mdot
            write(*,*) 'CpT', k, CpT
            write(*,*) 'qmid', k, qmid
            write(*,*) 'mmid', k, mmid
            write(*,*) 's% grada(k)', k, s% grada(k)
            write(*,*) 'inv_4pi_r4_P', k, inv_4pi_r4_P
            write(*,*) 'P', k, P
            write(*,*) 'rmid3', k, rmid3
            write(*,*) 'rmid4', k, rmid4
            write(*,*) 'inv_rmid4', k, inv_rmid4
            write(*,*) 'inv_4pi_r4', k, inv_4pi_r4
         end if
         
      end subroutine do_spatial_term


      subroutine blend_with_lnS_form(s, k, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: alfa, Gamma     
         include 'formats'
         ierr = 0
         Gamma = s% gam_start(k)
         alfa = (Gamma - s% Gamma_lnS_eps_grav_full_off) / &
            (s% Gamma_lnS_eps_grav_full_on - s% Gamma_lnS_eps_grav_full_off)
         ! alfa is fraction of lnS form in result
         if (alfa >= 1 .or. alfa <= 0) then
            ierr = -1
            if (s% report_ierr) &
               write(*,*) 'eval_eps_grav_and_partials -- error in blend_with_lnS_form'
            return
         end if         
         if (s% lnPgas_flag) then
            call combine_two_eps_gravs( &
               s, k, alfa, 1d0 - alfa, do_eps_grav_with_lnS, do_eps_grav_with_lnPgas, ierr)
         else
            call combine_two_eps_gravs( &
               s, k, alfa, 1d0 - alfa, do_eps_grav_with_lnS, do_eps_grav_with_lnd, ierr)
         end if      
      end subroutine blend_with_lnS_form


      recursive subroutine combine_two_eps_gravs( &
            s, k, alfa, beta, eps_grav_proc1, eps_grav_proc2, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(in) :: alfa, beta
         interface
            subroutine eps_grav_proc1(s, k, ierr)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s 
               integer, intent(in) :: k
                     integer, intent(out) :: ierr
            end subroutine eps_grav_proc1
            subroutine eps_grav_proc2(s, k, ierr)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s 
               integer, intent(in) :: k
                     integer, intent(out) :: ierr
            end subroutine eps_grav_proc2
         end interface
         integer, intent(out) :: ierr
         
         real(dp) :: &
            eps_grav, d_eps_grav_dlndm1, d_eps_grav_dlnd00, d_eps_grav_dlndp1, &
            d_eps_grav_dlnTm1, d_eps_grav_dlnT00, d_eps_grav_dlnTp1, &
            d_eps_grav_dlnR00, d_eps_grav_dlnRp1, d_eps_grav_dL00, d_eps_grav_dLp1, &
            d_eps_grav_dlnPgas00_const_T, &
            d_eps_grav_dlnPgasm1_const_T, d_eps_grav_dlnPgasp1_const_T, &
            d_eps_grav_dlnTm1_const_Pgas, d_eps_grav_dlnT00_const_Pgas, &
            d_eps_grav_dlnTp1_const_Pgas, d_eps_grav_dv00, d_eps_grav_dvp1

         include 'formats'
         ierr = 0

         ! alfa is multiplier of result from calling eps_grav_proc1
         ! beta is multiplier of result from calling eps_grav_proc2
         ! i.e., eps_grav = alfa*eps_grav1 + beta*eps_grav2

         if (alfa > 1d0 .or. alfa < 0d0 .or. beta > 1d0 .or. beta < 0d0) then
            if (s% report_ierr) &
               write(*,2) 'combine_two_eps_gravs: alfa beta', k, alfa, beta
            ierr = -1
            return
         end if
         
         ! result is alfa*eps_grav_proc1 + beta*eps_grav_proc2
         
         if (alfa > 0d0) then
            call eps_grav_proc1(s, k, ierr)
            if (ierr /= 0) return
            if (beta == 0d0) return
            ! save results
            eps_grav = s% eps_grav(k)
            d_eps_grav_dlndm1 = s% d_eps_grav_dlndm1(k)
            d_eps_grav_dlnd00 = s% d_eps_grav_dlnd00(k)
            d_eps_grav_dlndp1 = s% d_eps_grav_dlndp1(k)
            d_eps_grav_dlnTm1 = s% d_eps_grav_dlnTm1(k)
            d_eps_grav_dlnT00 = s% d_eps_grav_dlnT00(k)
            d_eps_grav_dlnTp1 = s% d_eps_grav_dlnTp1(k)
            d_eps_grav_dlnR00 = s% d_eps_grav_dlnR00(k)
            d_eps_grav_dlnRp1 = s% d_eps_grav_dlnRp1(k)
            d_eps_grav_dL00 = s% d_eps_grav_dL00(k)
            d_eps_grav_dLp1 = s% d_eps_grav_dLp1(k)
            d_eps_grav_dlnPgas00_const_T = s% d_eps_grav_dlnPgas00_const_T(k)
            d_eps_grav_dlnPgasm1_const_T = s% d_eps_grav_dlnPgasm1_const_T(k)
            d_eps_grav_dlnPgasp1_const_T = s% d_eps_grav_dlnPgasp1_const_T(k)
            d_eps_grav_dlnTm1_const_Pgas = s% d_eps_grav_dlnTm1_const_Pgas(k)
            d_eps_grav_dlnT00_const_Pgas = s% d_eps_grav_dlnT00_const_Pgas(k)
            d_eps_grav_dlnTp1_const_Pgas = s% d_eps_grav_dlnTp1_const_Pgas(k)
            d_eps_grav_dv00 = s% d_eps_grav_dv00(k)
            d_eps_grav_dvp1 = s% d_eps_grav_dvp1(k)
         else ! not needed, but to keep the compiler happy we set these to 0
            eps_grav = 0
            d_eps_grav_dlndm1 = 0
            d_eps_grav_dlnd00 = 0
            d_eps_grav_dlndp1 = 0
            d_eps_grav_dlnTm1 = 0
            d_eps_grav_dlnT00 = 0
            d_eps_grav_dlnTp1 = 0
            d_eps_grav_dlnR00 = 0
            d_eps_grav_dlnRp1 = 0
            d_eps_grav_dL00 = 0
            d_eps_grav_dLp1 = 0
            d_eps_grav_dlnPgas00_const_T = 0
            d_eps_grav_dlnPgasm1_const_T = 0
            d_eps_grav_dlnPgasp1_const_T = 0
            d_eps_grav_dlnTm1_const_Pgas = 0
            d_eps_grav_dlnT00_const_Pgas = 0
            d_eps_grav_dlnTp1_const_Pgas = 0
            d_eps_grav_dv00 = 0
            d_eps_grav_dvp1 = 0
         end if
         
         call eps_grav_proc2(s, k, ierr)
         if (ierr /= 0) return
         if (alfa == 0d0) return
         
         ! combine results
         s% eps_grav(k) = alfa*eps_grav + beta*s% eps_grav(k)
         
         s% d_eps_grav_dlndm1(k) = alfa*d_eps_grav_dlndm1 + beta*s% d_eps_grav_dlndm1(k)
         s% d_eps_grav_dlnd00(k) = alfa*d_eps_grav_dlnd00 + beta*s% d_eps_grav_dlnd00(k)
         s% d_eps_grav_dlndp1(k) = alfa*d_eps_grav_dlndp1 + beta*s% d_eps_grav_dlndp1(k)
         
         s% d_eps_grav_dlnTm1(k) = alfa*d_eps_grav_dlnTm1 + beta*s% d_eps_grav_dlnTm1(k)
         s% d_eps_grav_dlnT00(k) = alfa*d_eps_grav_dlnT00 + beta*s% d_eps_grav_dlnT00(k)
         s% d_eps_grav_dlnTp1(k) = alfa*d_eps_grav_dlnTp1 + beta*s% d_eps_grav_dlnTp1(k)
         
         s% d_eps_grav_dlnPgas00_const_T(k) = &
            alfa*d_eps_grav_dlnPgas00_const_T + beta*s% d_eps_grav_dlnPgas00_const_T(k)
         s% d_eps_grav_dlnPgasm1_const_T(k) = &
            alfa*d_eps_grav_dlnPgasm1_const_T + beta*s% d_eps_grav_dlnPgasm1_const_T(k)
         s% d_eps_grav_dlnPgasp1_const_T(k) = &
            alfa*d_eps_grav_dlnPgasp1_const_T + beta*s% d_eps_grav_dlnPgasp1_const_T(k)
            
         s% d_eps_grav_dlnTm1_const_Pgas(k) = &
            alfa*d_eps_grav_dlnTm1_const_Pgas + beta*s% d_eps_grav_dlnTm1_const_Pgas(k)
         s% d_eps_grav_dlnT00_const_Pgas(k) = &
            alfa*d_eps_grav_dlnT00_const_Pgas + beta*s% d_eps_grav_dlnT00_const_Pgas(k)
         s% d_eps_grav_dlnTp1_const_Pgas(k) = &
            alfa*d_eps_grav_dlnTp1_const_Pgas + beta*s% d_eps_grav_dlnTp1_const_Pgas(k)
         
         s% d_eps_grav_dlnR00(k) = alfa*d_eps_grav_dlnR00 + beta*s% d_eps_grav_dlnR00(k)
         s% d_eps_grav_dlnRp1(k) = alfa*d_eps_grav_dlnRp1 + beta*s% d_eps_grav_dlnRp1(k)
         
         s% d_eps_grav_dL00(k) = alfa*d_eps_grav_dL00 + beta*s% d_eps_grav_dL00(k)
         s% d_eps_grav_dLp1(k) = alfa*d_eps_grav_dLp1 + beta*s% d_eps_grav_dLp1(k)
            
         s% d_eps_grav_dv00(k) = alfa*d_eps_grav_dv00 + beta*s% d_eps_grav_dv00(k)
         s% d_eps_grav_dvp1(k) = alfa*d_eps_grav_dvp1 + beta*s% d_eps_grav_dvp1(k)
      
      end subroutine combine_two_eps_gravs
      
      
      subroutine include_dmu_dt_in_eps_grav(s, k)
         ! Lagrangian only since want to avoid doing a numerical difference for dmu/dm
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp) :: mu, mu_prev, dmu_dt, dE_dmu, d_dEdmu_dlnT
         include 'formats'
         if (.not. s% include_dmu_dt_in_eps_grav) return
         if (k < s% k_below_just_added) return
         mu = s% mu(k)
         mu_prev = s% prev_mu(k) ! interpolated to same mass coordinate
         dmu_dt = (mu - mu_prev)*s% dVARDOT_dVAR
         dE_dmu = -1.5d0*cgas*s% T(k)/(mu*mu)
         d_dEdmu_dlnT = dE_dmu
         s% eps_grav_composition_term(k) = -dE_dmu*dmu_dt
         s% eps_grav(k) = s% eps_grav(k) + s% eps_grav_composition_term(k)
         s% d_eps_grav_dlnT00(k) = s% d_eps_grav_dlnT00(k) - d_dEdmu_dlnT*dmu_dt
         
         if (k == s% trace_k) then
            write(*,5) 'include_dmu_dt_in_eps_grav', &
               k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, s% eps_grav(k)
            write(*,2) 'mu', k, mu
            write(*,2) 'mu_prev', k, mu_prev
            write(*,2) 'dmu_dt', k, dmu_dt
            write(*,2) 'dE_dmu', k, dE_dmu
         end if

      end subroutine include_dmu_dt_in_eps_grav
      
      
      subroutine zero_eps_grav_and_partials(s, k)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         s% eps_grav(k) = 0
         s% d_eps_grav_dlndm1(k) = 0
         s% d_eps_grav_dlnd00(k) = 0
         s% d_eps_grav_dlndp1(k) = 0
         s% d_eps_grav_dlnTm1(k) = 0
         s% d_eps_grav_dlnT00(k) = 0
         s% d_eps_grav_dlnTp1(k) = 0
         s% d_eps_grav_dlnPgasm1_const_T(k) = 0
         s% d_eps_grav_dlnPgas00_const_T(k) = 0
         s% d_eps_grav_dlnPgasp1_const_T(k) = 0
         s% d_eps_grav_dlnTm1_const_Pgas(k) = 0
         s% d_eps_grav_dlnT00_const_Pgas(k) = 0
         s% d_eps_grav_dlnTp1_const_Pgas(k) = 0
         s% d_eps_grav_dlnR00(k) = 0
         s% d_eps_grav_dlnRp1(k) = 0
         s% d_eps_grav_dL00(k) = 0
         s% d_eps_grav_dLp1(k) = 0
         s% d_eps_grav_dv00(k) = 0
         s% d_eps_grav_dvp1(k) = 0
      end subroutine zero_eps_grav_and_partials

         
      end module eps_grav

