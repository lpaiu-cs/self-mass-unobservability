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

      module star_private_def
      
      use star_def
      use crlibm_lib

      implicit none


      ! phases of evolution for logs
      integer, parameter :: phase_starting = 0
      integer, parameter :: phase_early_main_seq = 1
      integer, parameter :: phase_mid_main_seq = 2
      integer, parameter :: phase_wait_for_he = 3
      integer, parameter :: phase_he_ignition_over = 4
      integer, parameter :: phase_he_igniting = 5
      integer, parameter :: phase_helium_burning = 6
      integer, parameter :: phase_carbon_burning = 7

         
      real(dp), parameter :: del_cntr_rho = 1d0
      real(dp), parameter :: min_cntr_rho = 3d0
      real(dp), parameter :: no_he_ignition_limit = 0.75d0
      real(dp), parameter :: no_cntr_T_drops_limit = 6.5d0
      
      real(dp), parameter :: center_h_gone = 1d-3
      real(dp), parameter :: center_h_going = 1d0/3d0
      real(dp), parameter :: center_he_going = 5d-2


         

      contains


      subroutine star_private_def_init
         use num_def
         integer :: i      
         logical :: okay
         
         include 'formats'
         
         okay = .true.
         
         termination_code_str(1:num_termination_codes) = ''
         
         termination_code_str(t_max_age) = 'max_age'
         termination_code_str(t_max_omega_div_omega_crit) = 'max_omega_div_omega_crit'
         termination_code_str(t_peak_burn_vconv_div_cs_limit) = 'peak_burn_vconv_div_cs_limit'
         termination_code_str(t_max_model_number) = 'max_model_number'
         termination_code_str(t_eta_center_limit) = 'eta_center_limit'
         termination_code_str(t_log_center_temp_limit) = 'log_center_temp_limit'
         termination_code_str(t_log_center_temp_lower_limit) = 'log_center_temp_lower_limit'
         termination_code_str(t_center_entropy_limit) = 'center_entropy_limit'
         termination_code_str(t_center_entropy_lower_limit) = 'center_entropy_lower_limit'
         termination_code_str(t_max_entropy_limit) = 'max_entropy_limit'
         termination_code_str(t_max_entropy_lower_limit) = 'max_entropy_lower_limit'
         termination_code_str(t_log_center_density_limit) = 'log_center_density_limit'
         termination_code_str(t_log_center_density_lower_limit) = 'log_center_density_lower_limit'
         termination_code_str(t_gamma_center_limit) = 'gamma_center_limit'
         termination_code_str(t_log_max_temp_upper_limit) = 'log_max_temp_upper_limit'
         termination_code_str(t_log_max_temp_lower_limit) = 'log_max_temp_lower_limit'
         termination_code_str(t_HB_limit) = 'HB_limit'
         termination_code_str(t_stop_at_TP) = 'stop_at_TP'
         termination_code_str(t_star_mass_min_limit) = 'star_mass_min_limit'
         termination_code_str(t_star_mass_max_limit) = 'star_mass_max_limit'
         termination_code_str(t_star_H_mass_max_limit) = 'star_H_mass_max_limit'
         termination_code_str(t_star_H_mass_min_limit) = 'star_H_mass_min_limit'
         termination_code_str(t_star_He_mass_max_limit) = 'star_He_mass_max_limit'
         termination_code_str(t_star_He_mass_min_limit) = 'star_He_mass_min_limit'
         termination_code_str(t_star_C_mass_max_limit) = 'star_C_mass_max_limit'
         termination_code_str(t_star_C_mass_min_limit) = 'star_C_mass_min_limit'
         termination_code_str(t_xmstar_min_limit) = 'xmstar_min_limit'
         termination_code_str(t_xmstar_max_limit) = 'xmstar_max_limit'
         termination_code_str(t_envelope_mass_limit) = 'envelope_mass_limit'
         termination_code_str(t_envelope_fraction_left_limit) = 'envelope_fraction_left_limit'
         
         termination_code_str(t_he_core_mass_limit) = 'he_core_mass_limit'
         termination_code_str(t_c_core_mass_limit) = 'c_core_mass_limit'
         termination_code_str(t_o_core_mass_limit) = 'o_core_mass_limit'
         termination_code_str(t_si_core_mass_limit) = 'si_core_mass_limit'
         termination_code_str(t_fe_core_mass_limit) = 'fe_core_mass_limit'
         termination_code_str(t_neutron_rich_core_mass_limit) = 'neutron_rich_core_mass_limit'
         
         termination_code_str(t_he_layer_mass_lower_limit) = 'he_layer_mass_lower_limit'
         termination_code_str(t_abs_diff_lg_LH_lg_Ls_limit) = 'abs_diff_lg_LH_lg_Ls_limit'
         termination_code_str(t_Teff_lower_limit) = 'Teff_lower_limit'
         termination_code_str(t_Teff_upper_limit) = 'Teff_upper_limit'
         termination_code_str(t_delta_nu_lower_limit) = 'delta_nu_lower_limit'
         termination_code_str(t_delta_nu_upper_limit) = 'delta_nu_upper_limit'
         termination_code_str(t_delta_Pg_lower_limit) = 'delta_Pg_lower_limit'
         termination_code_str(t_delta_Pg_upper_limit) = 'delta_Pg_upper_limit'
         termination_code_str(t_photosphere_r_lower_limit) = 'photosphere_r_lower_limit'
         termination_code_str(t_photosphere_r_upper_limit) = 'photosphere_r_upper_limit'
         termination_code_str(t_log_Teff_lower_limit) = 'log_Teff_lower_limit'
         termination_code_str(t_log_Teff_upper_limit) = 'log_Teff_upper_limit'
         termination_code_str(t_log_Tsurf_lower_limit) = 'log_Tsurf_lower_limit'
         termination_code_str(t_log_Tsurf_upper_limit) = 'log_Tsurf_upper_limit'
         termination_code_str(t_log_Psurf_lower_limit) = 'log_Psurf_lower_limit'
         termination_code_str(t_log_Psurf_upper_limit) = 'log_Psurf_upper_limit'
         termination_code_str(t_log_Dsurf_lower_limit) = 'log_Dsurf_lower_limit'
         termination_code_str(t_log_Dsurf_upper_limit) = 'log_Dsurf_upper_limit'
         termination_code_str(t_log_L_lower_limit) = 'log_L_lower_limit'
         termination_code_str(t_log_L_upper_limit) = 'log_L_upper_limit'
         termination_code_str(t_log_g_lower_limit) = 'log_g_lower_limit'
         termination_code_str(t_log_g_upper_limit) = 'log_g_upper_limit'
         termination_code_str(t_power_nuc_burn_upper_limit) = 'power_nuc_burn_upper_limit'
         termination_code_str(t_power_h_burn_upper_limit) = 'power_h_burn_upper_limit'
         termination_code_str(t_power_he_burn_upper_limit) = 'power_he_burn_upper_limit'
         termination_code_str(t_power_c_burn_upper_limit) = 'power_c_burn_upper_limit'
         termination_code_str(t_power_nuc_burn_lower_limit) = 'power_nuc_burn_lower_limit'
         termination_code_str(t_power_h_burn_lower_limit) = 'power_h_burn_lower_limit'
         termination_code_str(t_power_he_burn_lower_limit) = 'power_he_burn_lower_limit'
         termination_code_str(t_power_c_burn_lower_limit) = 'power_c_burn_lower_limit'
         termination_code_str(t_phase_of_evolution_stop) = 'phase_of_evolution_stop'
         termination_code_str(t_center_Ye_lower_limit) = 'center_Ye_lower_limit'
         termination_code_str(t_fe_core_infall_limit) = 'fe_core_infall_limit'
         termination_code_str(t_non_fe_core_infall_limit) = 'non_fe_core_infall_limit'
         termination_code_str(t_v_div_csound_surf_limit) = 'v_div_csound_surf_limit'
         termination_code_str(t_Pgas_div_P_limit) = 'Pgas_div_P_limit'
         termination_code_str(t_Lnuc_div_L_lower_limit) = 'Lnuc_div_L_lower_limit'
         termination_code_str(t_Lnuc_div_L_upper_limit) = 'Lnuc_div_L_upper_limit'
         termination_code_str(t_v_surf_div_v_kh_lower_limit) = 'v_surf_div_v_kh_lower_limit'
         termination_code_str(t_v_surf_div_v_kh_upper_limit) = 'v_surf_div_v_kh_upper_limit'
         termination_code_str(t_v_surf_div_v_esc_limit) = 'v_surf_div_v_esc_limit'
         termination_code_str(t_Lnuc_div_L_zams_limit) = 'Lnuc_div_L_zams_limit'
         termination_code_str(t_logQ_limit) = 'logQ_limit'
         termination_code_str(t_xa_central_lower_limit) = 'xa_central_lower_limit'
         termination_code_str(t_xa_central_upper_limit) = 'xa_central_upper_limit'
         termination_code_str(t_xa_surface_lower_limit) = 'xa_surface_lower_limit'
         termination_code_str(t_xa_surface_upper_limit) = 'xa_surface_upper_limit'
         termination_code_str(t_xa_average_lower_limit) = 'xa_average_lower_limit'
         termination_code_str(t_xa_average_upper_limit) = 'xa_average_upper_limit'
         termination_code_str(t_surface_accel_div_grav_limit) = 'surface_accel_div_grav_limit'
         termination_code_str(t_stop_because_he_ignited) = 'stop_because_he_ignited'
         termination_code_str(t_adjust_mesh_failed) = 'adjust_mesh_failed'
         termination_code_str(t_dt_is_zero) = 'dt_is_zero'
         termination_code_str(t_min_timestep_limit) = 'min_timestep_limit'
         termination_code_str(t_failed_prepare_for_new_try) = 'failed_prepare_for_new_try'
         termination_code_str(t_negative_total_angular_momentum) = 'negative_total_angular_momentum'
         termination_code_str(t_max_number_retries) = 'max_number_retries'
         termination_code_str(t_max_number_backups) = 'max_number_backups'
         termination_code_str(t_max_backups_in_a_row) = 'max_backups_in_a_row'
         termination_code_str(t_solve_burn) = 'solve_burn'
         termination_code_str(t_solve_hydro) = 'solve_hydro'
         termination_code_str(t_solve_mix) = 'solve_mix'
         termination_code_str(t_solve_omega_mix) = 'solve_omega_mix'
         termination_code_str(t_timestep_controller) = 'timestep_controller'
         termination_code_str(t_relax_finished_okay) = 'relax_finished_okay'
         termination_code_str(t_done_with_piston) = 'done with piston'
         termination_code_str(t_done_with_center_flash) = 'done with center flash'

         termination_code_str(t_extras_check_model) = 'extras_check_model'
         termination_code_str(t_extras_finish_step) = 'extras_finish_step'

         termination_code_str(t_xtra1) = 'customize by setting termination_code_str(t_xtra1)'
         termination_code_str(t_xtra2) = 'customize by setting termination_code_str(t_xtra2)'
         termination_code_str(t_xtra3) = 'customize by setting termination_code_str(t_xtra3)'
         termination_code_str(t_xtra4) = 'customize by setting termination_code_str(t_xtra4)'
         termination_code_str(t_xtra5) = 'customize by setting termination_code_str(t_xtra5)'
         termination_code_str(t_xtra6) = 'customize by setting termination_code_str(t_xtra6)'
         termination_code_str(t_xtra7) = 'customize by setting termination_code_str(t_xtra7)'
         termination_code_str(t_xtra8) = 'customize by setting termination_code_str(t_xtra8)'
         termination_code_str(t_xtra9) = 'customize by setting termination_code_str(t_xtra9)'
         
         do i=1,num_termination_codes
            if (len_trim(termination_code_str(i)) == 0) then
               if (i > 1) then
                  write(*,2) 'missing termination_code_str following ' // &
                     trim(termination_code_str(i-1)), i
               else
                  write(*,2) 'missing termination_code_str 1'
               end if
               okay = .false.
            end if
         end do

         dt_why_str(1:numTlim) = ''
         
         dt_why_str(Tlim_struc) = 'varcontrol'
         dt_why_str(Tlim_max_timestep_factor) = 'max increase'
         dt_why_str(Tlim_min_timestep_factor) = 'max decrease'
         dt_why_str(Tlim_num_newton_iterations) = 'newton iters'
         dt_why_str(Tlim_num_rotation_steps) = 'rotation steps'
         dt_why_str(Tlim_num_diff_solver_steps) = 'diff steps'
         dt_why_str(Tlim_num_diff_solver_iters) = 'diff iters'
         dt_why_str(Tlim_burn_max_num_substeps) = 'burn steps'
         dt_why_str(Tlim_burn_max_num_iters) = 'burn iters'
         dt_why_str(Tlim_max_mix_fixup) = 'mix fixup'
         dt_why_str(Tlim_dH) = 'dH'
         dt_why_str(Tlim_dHe) = 'dHe'
         dt_why_str(Tlim_dHe3) = 'dHe3'
         dt_why_str(Tlim_dX) = 'dX'
         dt_why_str(Tlim_dH_div_H) = 'dH/H'
         dt_why_str(Tlim_dHe_div_He) = 'dHe/He'
         dt_why_str(Tlim_dHe3_div_He3) = 'dHe3/He3'
         dt_why_str(Tlim_dX_div_X) = 'dX/X'
         dt_why_str(Tlim_dL_div_L) = 'dL/L'
         dt_why_str(Tlim_dlgP) = 'lgP'
         dt_why_str(Tlim_dlgRho) = 'lgRho'
         dt_why_str(Tlim_dlgT) = 'lgT'
         dt_why_str(Tlim_dlgE) = 'lgE'
         dt_why_str(Tlim_dlgR) = 'lgR'
         dt_why_str(Tlim_d_deltaR_shrink) = 'd_delR_shrink'
         dt_why_str(Tlim_d_deltaR_grow) = 'd_delR_grow'
         dt_why_str(Tlim_dlgL_nuc_cat) = 'Lnuc_cat'
         dt_why_str(Tlim_dlgL_H) = 'Lnuc_H'
         dt_why_str(Tlim_dlgL_He) = 'Lnuc_He'
         dt_why_str(Tlim_dlgL_z) = 'Lnuc_z'
         dt_why_str(Tlim_dlgL_photo) = 'Lnuc_photo'
         dt_why_str(Tlim_dlgL_nuc) = 'Lnuc'
         dt_why_str(Tlim_dlgTeff) = 'lgTeff'
         dt_why_str(Tlim_dlgRho_cntr) = 'lgRho_cntr'
         dt_why_str(Tlim_dlgT_cntr) = 'lgT_cntr'
         dt_why_str(Tlim_dlgT_max) = 'lgT_max'
         dt_why_str(Tlim_dlgRho_max) = 'lgRho_max'
         dt_why_str(Tlim_dlog_eps_nuc_cntr) = 'eps_nuc_cntr'
         dt_why_str(Tlim_dlog_eps_nuc) = 'log_eps_nuc'
         dt_why_str(Tlim_lg_XH_cntr) = 'lg_XH_cntr'
         dt_why_str(Tlim_dmstar) = 'delta_mstar'
         dt_why_str(Tlim_del_mdot) = 'delta mdot'
         dt_why_str(Tlim_deltaJ) = 'delta total J'
         dt_why_str(Tlim_dt_div_dt_cell_collapse) = 'dt_collapse'
         dt_why_str(Tlim_dt_div_dt_Courant) = 'dt_Courant'
         dt_why_str(Tlim_dt_div_dt_thermal) = 'dt_thermal'
         dt_why_str(Tlim_dt_div_dt_dynamic) = 'dt_dynamic'
         dt_why_str(Tlim_dt_div_dt_acoustic) = 'dt_acoustic'
         dt_why_str(Tlim_dt_div_dt_mass_loss) = 'dt_mass_loss'
         dt_why_str(Tlim_lgL) = 'lgL'
         dt_why_str(Tlim_max_timestep) = 'max_dt'
         dt_why_str(Tlim_timestep_hold) = 'hold'
         dt_why_str(Tlim_dX_nuc_drop) = 'dX_nuc_drop'
         dt_why_str(Tlim_lg_XHe_cntr) = 'lg_XHe_cntr'
         dt_why_str(Tlim_lg_XC_cntr) = 'lg_XC_cntr'
         dt_why_str(Tlim_lg_XNe_cntr) = 'lg_XNe_cntr'
         dt_why_str(Tlim_lg_XO_cntr) = 'lg_XO_cntr'
         dt_why_str(Tlim_lg_XSi_cntr) = 'lg_XSi_cntr'
         dt_why_str(Tlim_CpT_absMdot_div_L) = 'CpT_absMdot_div_L'
         dt_why_str(Tlim_neg_X) = 'neg_mass_frac'
         dt_why_str(Tlim_bad_Xsum) = 'bad_X_sum'
         dt_why_str(Tlim_delta_HR) = 'delta_HR'
         dt_why_str(Tlim_delta_Ye) = 'delta Ye'
         dt_why_str(Tlim_delta_Ye_highT) = 'highT del Ye'
         dt_why_str(Tlim_error_rate_energy_conservation) = 'error rate'
         dt_why_str(Tlim_avg_E_residual) = 'avg E resid'
         dt_why_str(Tlim_max_E_residual) = 'max E resid'
         dt_why_str(Tlim_error_in_energy_conservation) = 'error enrg'
         dt_why_str(Tlim_truncation_ratio) = 'trnc ratio'
         dt_why_str(Tlim_v_div_v_crit) = 'v/v_crit'
         dt_why_str(Tlim_binary) = 'binary'
         dt_why_str(Tlim_error_other) = 'other'
         
         do i=1,numTlim
            if (len_trim(dt_why_str(i)) == 0) then
               if (i > 1) then
                  write(*,2) 'missing dt_why_str following ' // trim(dt_why_str(i-1)), i
               else
                  write(*,2) 'missing dt_why_str 1'
               end if
               okay = .false.
            end if
         end do
         
         if (.not. okay) stop 'star_private_def_init'
                  
      end subroutine star_private_def_init

      
      integer function alloc_star(ierr)
         use rates_def, only: rates_NACRE_if_available
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s
         
         ierr = 0
         alloc_star = -1
!$omp critical (star_handle)
         if (.not. have_initialized_star_handles) then
            do i = 1, max_star_handles
               star_handles(i)% id = i
               star_handles(i)% in_use = .false.
            end do
            have_initialized_star_handles = .true.
         end if
         do i = 1, max_star_handles
            if (.not. star_handles(i)% in_use) then
               star_handles(i)% in_use = .true.
               alloc_star = i
               exit
            end if
         end do
!$omp end critical (star_handle)
         if (alloc_star == -1) then
            ierr = -1
            return
         end if
         if (star_handles(alloc_star)% id /= alloc_star) then
            ierr = -1
            return
         end if
         s => star_handles(alloc_star)
         
      end function alloc_star
      
      
      subroutine free_star(s)
         type (star_info), pointer :: s
         star_handles(s% id)% in_use = .false.
      end subroutine free_star
      
      
      subroutine stardata_init( &
            my_mesa_dir, chem_isotopes_filename, &
            kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
            kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, kappa_type2_logT_lower_bdy, &
            eos_file_prefix, eosDT_Z1_suffix, eosPT_Z1_suffix, &
            net_reaction_filename, jina_reaclib_filename, &
            use_special_weak_rates, special_weak_states_file, special_weak_transitions_file, &
            reaclib_min_T9_in, &
            rate_tables_dir, rates_cache_suffix, &
            ionization_file_prefix, ionization_Z1_suffix, &
            eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir, &
            ionization_cache_dir, kap_cache_dir, rates_cache_dir, &
            ierr)
         use colors_lib
         use eos_lib, only: eos_init, set_eos_Z1_suffix
         use kap_lib, only: kap_init
         use rates_lib, only: rates_init
         use rates_def, only: reaclib_min_T9
!         use nse_lib, only: nse_init
         use net_lib, only: net_init
         use ionization_lib, only: ionization_init
         use atm_lib
         use chem_lib
         use const_lib
         use const_def, only: mesa_data_dir
         use utils_lib
         use star_history_def, only: history_column_names_init
         use star_profile_def, only: profile_column_names_init
         character (len=*), intent(in) :: &
            my_mesa_dir, chem_isotopes_filename, net_reaction_filename, &
            jina_reaclib_filename, rate_tables_dir, &
            special_weak_states_file, special_weak_transitions_file, &
            kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
            eosDT_Z1_suffix, eosPT_Z1_suffix, &
            eos_file_prefix, rates_cache_suffix, &
            ionization_file_prefix, ionization_Z1_suffix, &
            eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir, &
            ionization_cache_dir, kap_cache_dir, rates_cache_dir
         real(dp), intent(in) :: kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
            kappa_type2_logT_lower_bdy, reaclib_min_T9_in
         logical, intent(in) :: use_special_weak_rates
         integer, intent(out) :: ierr

         logical, parameter :: use_cache = .true.
         character (len=strlen) :: fname
         integer :: iounit
         
         logical, parameter :: dbg = .false.
         
         include 'formats'

         ierr = 0
         
         rate_tables_dir_for_star = rate_tables_dir
         rates_cache_suffix_for_star = rates_cache_suffix
         
         if (dbg) write(*,*) 'call const_init'
         call const_init(my_mesa_dir,ierr)
      	if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call crlibm_init'
         call crlibm_init
           
         call star_private_def_init
         call result_reason_init

         call history_column_names_init(ierr)
      	if (ierr /= 0) return

         call profile_column_names_init(ierr)
      	if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call chem_init'
         call chem_init(chem_isotopes_filename, ierr)
      	if (ierr /= 0) return

         if (dbg) write(*,*) 'call colors_init'
         call colors_init(ierr)
      	if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call eos_init'
         call eos_init( &
            eos_file_prefix, &
            eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir, &
            use_cache, ierr)
      	if (ierr /= 0) return
         
         call set_eos_Z1_suffix(eosDT_Z1_suffix, eosPT_Z1_suffix) 
         
         if (dbg) write(*,*) 'call kap_init'
         call kap_init( &
            kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
            kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
            kappa_type2_logT_lower_bdy, use_cache, kap_cache_dir, &
            ierr)
      	if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call rates_init'
         call rates_init( &
            net_reaction_filename, jina_reaclib_filename, &
            use_special_weak_rates, &
            special_weak_states_file, &
            special_weak_transitions_file, &
            rates_cache_dir, ierr)
         
      	if (ierr /= 0) return
      	
      	if (reaclib_min_T9_in > 0 .and. reaclib_min_T9_in /= reaclib_min_T9) then
      	   reaclib_min_T9 = reaclib_min_T9_in
      	   write(*,*)
      	   write(*,*)
      	   write(*,*)
      	   write(*,*)
      	   write(*,1) 'change reaclib_min_T9', reaclib_min_T9
      	   write(*,1) 'must clear data/rates_data/cache of old reaclib rates'
      	   write(*,*)
      	   write(*,*)
      	   write(*,*)
      	   write(*,*)
      	   write(*,*)
      	end if
                  
         if (dbg) write(*,*) 'call net_init'
         call net_init(ierr)
      	if (ierr /= 0) return
                  
!         if (dbg) write(*,*) 'call nse_init'
!         call nse_init(ierr)   
!      	if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call atm_init'
         call atm_init(use_cache, ierr)
      	if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call ionization_init'
         call ionization_init( &
            ionization_file_prefix, ionization_Z1_suffix, &
            ionization_cache_dir, use_cache, ierr)
      	if (ierr /= 0) return
         
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         
         version_number = -1
         fname = trim(mesa_data_dir) // '/version_number'
         open(iounit, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open ' // trim(fname)
            call free_iounit(iounit)
            return
         end if                  
         read(iounit, *, iostat=ierr) version_number
         if (ierr /= 0) then
            close(iounit)
            call free_iounit(iounit)
            return
         end if
         write(*,2) 'version_number', version_number
         
      end subroutine stardata_init


      end module star_private_def

