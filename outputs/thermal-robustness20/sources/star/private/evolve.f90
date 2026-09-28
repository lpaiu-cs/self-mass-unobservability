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

      module evolve

      use star_private_def
      use const_def
      use star_utils

      implicit none

      contains
      
      integer function do_evolve_step_part1(id, first_try, just_did_backup)
         use winds, only: set_mdot
         use alloc, only: check_sizes
         use do_one_utils, only: write_terminal_header
         logical, intent(in) :: first_try, just_did_backup
         integer, intent(in) :: id

         type (star_info), pointer :: s            
         integer :: ierr, j, k, time0, clock_rate
         logical :: trace
         
         logical, parameter :: dbg = .false.

         include 'formats'

         do_evolve_step_part1 = terminate
         
         ierr = 0         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return         

         s% termination_code = 0
         s% retry_message = ''

         call system_clock(s% system_clock_at_start_of_step, clock_rate)
         
         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do j=1,s% species
               write(*,4) 'evolve starting xa(j)', &
                  s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
            end do
            if (s% rotation_flag) then
               do k=1,s% nz
                  write(*,3) 'evolve starting j_rot(k)', &
                     s% model_number, k, s% j_rot(k)
               end do
            end if
         end if

         trace = s% trace_evolve
         s% need_to_reset_cumulative_energy_info = .false.
         s% just_did_backup = just_did_backup
         s% doing_newton_iterations = .false.
         s% num_rotation_solver_steps = 0
         s% have_mixing_info = .false.
         dt_why_str(Tlim_max_timestep_factor) = 'max increase'   
               
         if (dbg) then
            call check_sizes(s, ierr)
            if (ierr /= 0) then
               write(*,*) 'do_evolve_step_part1: check_sizes returned ierr', ierr
               return
            end if
         end if

         if (s% doing_first_model_of_run) then
            if (s% do_history_file) then
               if (first_try) then
                  call write_terminal_header(s)
               else
                  write(*,1) '1st model retry log10(dt/yr)', log10_cr(s% dt/secyer)
               end if
            end if
            call system_clock(time0,clock_rate)
            s% starting_system_clock_time = time0
            s% system_clock_rate = clock_rate
            s% initial_timestep = s% dt_next
            s% last_backup = -s% backup_hold
            s% initial_L_center = s% L_center
            s% initial_R_center = s% R_center
            s% initial_v_center = s% v_center
            s% doing_center_flash = (s% center_flash_total_time > 0d0)
            if (s% doing_center_flash .and. s% v_center /= 0) then
               s% v_center = 0
               write(*,*) 'automatically set s% v_center = 0 at start of center flash'
            end if
            if (s% use_piston .and. s% time == 0 .and. s% piston_period <= 0) then   
               !s% piston_inward_time = (s% piston_Rmin - s% R_center)/s% v_center
               s% time = -s% piston_inward_time
               s% star_age = s% time/secyer
               s% profile_age = s% star_age
               s% post_he_age = s% star_age
               !write(*,1) 'set starting time to piston_inward_time', s% time
               if (s% piston_inward_time > 0) then
                  s% piston_vfinal_inward = &
                     2*(s% piston_Rmin - s% R_center)/s% piston_inward_time - s% v_center
                  if (.false.) then
                     write(*,1) 's% piston_Rmin', s% piston_Rmin
                     write(*,1) 's% R_center', s% R_center
                     write(*,1) 's% piston_Rmin - s% R_center', s% piston_Rmin - s% R_center
                     write(*,1) 'piston_vfinal_inward', s% piston_vfinal_inward
                     write(*,1) 'initial_v_center', s% initial_v_center
                     !stop
                  end if
               end if
               if (s% piston_Rmax > 0 .and. s% piston_Rmin > 0) then
                  s% piston_alpha = -s% piston_v0*s% piston_v0 / &
                     (standard_cgrav*s% M_center* &
                        (1/s% piston_Rmax - 1/s% piston_Rmin))
                  !write(*,2) 'piston_alpha', s% model_number, s% piston_alpha
                  !stop
               end if
               !stop
            end if
            if (s% min_ebdf_order > 0) then
               s% ebdf_order = s% min_ebdf_order
               s% startup_increment_ebdf_order = .true.
            end if
            s% timestep_hold = -111
            if (first_try) s% model_number_old = s% model_number
         end if

         call debug('before prepare_for_new_step', s)
         
         if (first_try) then
            if (trace) write(*,'(/,a,i8)') 'call prepare_for_new_step', s% model_number
            do_evolve_step_part1 = prepare_for_new_step(s)
            if (do_evolve_step_part1 /= keep_going) return
         end if

         call debug('before prepare_for_new_try', s)
         if (trace) write(*,'(/,a,i8)') 'call prepare_for_new_try', s% model_number
         do_evolve_step_part1 = prepare_for_new_try(s)
         if (do_evolve_step_part1 /= keep_going) return

         call debug('before set_mdot', s)
         if (trace) write(*,'(/,a,i8)') 'call set_mdot', s% model_number

         ! set mdot for the step
         call set_mdot(s, s% L_phot*Lsun, s% mstar, s% Teff, ierr)
         if (ierr /= 0) then
            do_evolve_step_part1 = retry
            s% result_reason = nonzero_ierr
            if (s% report_ierr) write(*, *) 'do_evolve_step_part1: set_mdot ierr'
            return
         end if
         
         if (s% use_other_adjust_mdot) then
            call s% other_adjust_mdot(s% id, ierr)
            if (ierr /= 0) then
               do_evolve_step_part1 = retry
               s% result_reason = nonzero_ierr
               if (s% report_ierr) write(*, *) 'do_evolve_step_part1: other_adjust_mdot ierr'
               return
            end if
         end if

      end function do_evolve_step_part1


      integer function do_evolve_step_part2(id, first_try, just_did_backup)
         use num_def
         use report, only: do_report
         use adjust_mass, only: do_adjust_mass         
         use utils_lib, only: is_bad_num, has_bad_num
         use element_diffusion, only: do_element_diffusion
         use evolve_support, only: set_current_to_old
         use solve_hydro, only: set_L_burn_by_category
         use struct_burn_mix, only: do_struct_burn_mix
         use hydro_vars, only: set_vars, set_final_vars
         use hydro_mtx, only: dump_struct
         use star_utils, only: start_time, update_time, &
            get_total_energy_integral
         use solve_omega_mix, only: do_solve_omega_mix
         use mix_info, only: set_mixing_info
         use hydro_rotation, only: set_rotation_info
         use profile
         
         logical, intent(in) :: first_try, just_did_backup
         integer, intent(in) :: id

         type (star_info), pointer :: s            
         integer :: ierr, time0, clock_rate, &
            j, k, j_cnt, mdot_redo_cnt, max_mdot_redo_cnt, cnt, max_cnt, nz
         logical :: okay, trace, skip_global_corr_coeff_limit, &
            have_too_large_wind_mdot, have_too_small_wind_mdot, &
            ignored_first_step, was_in_implicit_wind_limit
         real(dp) :: J_tot1, J_tot2, rel_error, piston_displacement, piston_t_end, &
            w_div_w_crit, w_div_w_crit_prev, mstar_dot, mstar_dot_prev, abs_mstar_delta, &
            explicit_mdot, max_wind_mdot, wind_mdot, r_phot, kh_timescale, dmskhf, dmsfac, &
            too_large_wind_mdot, too_small_wind_mdot, boost, mstar_dot_nxt, total, &
            surf_w_div_w_crit_limit, dt, r_nz, v_nz, time, max_dt, dv, total_energy, &
            new_R_center, time_new, amplitude, flash_max
         
         logical, parameter :: dbg = .false.
         
         include 'formats'

         ierr = 0         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return         

         time0 = s% starting_system_clock_time
         clock_rate = s% system_clock_rate
         trace = s% trace_evolve
         ignored_first_step = .false.
         
         mstar_dot = 0
         w_div_w_crit = -1
         surf_w_div_w_crit_limit = s% surf_w_div_w_crit_limit
         mdot_redo_cnt = 0
         max_mdot_redo_cnt = s% max_mdot_redo_cnt
         
         max_wind_mdot = 10*Msun/secyer
         have_too_large_wind_mdot = .false.
         have_too_small_wind_mdot = .false.
         too_large_wind_mdot = 0
         too_small_wind_mdot = 0

         explicit_mdot = s% mstar_dot
         
         was_in_implicit_wind_limit = s% was_in_implicit_wind_limit
         if (was_in_implicit_wind_limit .and. &
             s% generations >= 2 .and. &
               abs(s% mstar_dot_old) > 0 .and. &
               abs((s% mstar_dot-s% mstar_dot_old)/s% mstar_dot_old)+1 > &
                  s% mdot_revise_factor) then
             write(*,*) "Skipping first step in implicit mdot"
             s% mstar_dot = s% mstar_dot_old
             mdot_redo_cnt = 1
             ignored_first_step = .true.
         end if
         
         if (s% dt <= 0) then
            do_evolve_step_part2 = terminate
            s% termination_code = t_dt_is_zero
            s% result_reason = dt_is_zero
            return
         end if
      
         abs_mstar_delta = 0 
         
         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do j=1,s% species
               write(*,4) 'evolve after remesh xa(j)', &
                  s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
            end do
         end if
         
         if (s% done_with_piston .and. s% v_center /= 0d0) s% v_center = 0d0
         
      implicit_mdot_loop: do

            time_new = s% time_old + s% dt ! haven't updated s% time yet
            
            if (s% doing_center_flash) then
               if (time_new < s% center_flash_total_time) then ! boost L_center
                  flash_max = s% center_flash_total_ergs / &
                     (s% center_flash_total_time - 0.5d0*( &
                        s% center_flash_ramp_up_duration + &
                        s% center_flash_ramp_down_duration))
                  if (time_new < s% center_flash_ramp_up_duration) then
                     s% L_center = s% initial_L_center + &
                        flash_max*(time_new/s% center_flash_ramp_up_duration) 
                     write(*,2) 'ramp up: flash ergs', s% model_number, &
                        flash_max*0.5d0*time_new*time_new/s% center_flash_ramp_up_duration                        
                  else if (time_new > &
                        s% center_flash_total_time - s% center_flash_ramp_down_duration) then
                     s% L_center = s% initial_L_center + &
                        flash_max*(s% center_flash_total_time - time_new)/ &
                           s% center_flash_ramp_down_duration
                     write(*,2) 'ramp down: flash ergs', s% model_number, &
                        s% center_flash_total_ergs - &
                        flash_max* &
                           0.5d0*(s% center_flash_total_time - time_new)* &
                              (s% center_flash_total_time - time_new)/ &
                                 s% center_flash_ramp_down_duration
                  else ! during the flat center period
                     s% L_center = s% initial_L_center + flash_max
                     write(*,2) 'plateau: flash ergs', s% model_number, &
                        flash_max*(time_new - 0.5d0*s% center_flash_ramp_up_duration)
                  end if
               else
                  s% L_center = s% initial_L_center
                  s% done_with_center_flash = .true.
               end if
            end if
            
            if (.not. s% done_with_piston) then
            
               if (s% use_piston) then
                  if (s% piston_period > 0) then ! periodic piston
                     piston_t_end = s% piston_period*s% periodic_piston_number_of_cycles
                     amplitude = max(0d0,min(1d0,(piston_t_end - time_new)/piston_t_end))* &
                        s% periodic_piston_max_displacement
                     if (s% dt > 0) then
                        if (time_new < piston_t_end) then
                           piston_displacement = &
                              amplitude*sin_cr(2*pi*time_new/s% piston_period)
                        else
                           piston_displacement = &
                              amplitude*sin_cr(2*pi*s% periodic_piston_number_of_cycles)
                        end if
                        new_R_center = s% initial_R_center + piston_displacement
                        if (new_R_center < 0d0) then
                           write(*,2) 'new_R_center', s% model_number, new_R_center
                           do_evolve_step_part2 = terminate
                           return
                        endif
                        s% v_center = (new_R_center - s% R_center_old)/s% dt
                        if (s% v_center < 0 .and. s% stop_when_piston_v_goes_negative) then
                           do_evolve_step_part2 = terminate
                           s% termination_code = t_done_with_piston
                           return
                        end if
                     end if
                  else if (s% v_center < 0 .and. s% R_center > s% piston_Rmin) then
                     ! accelerate until reach piston_vfinal_inward
                     s% v_center = max(s% piston_vfinal_inward, s% initial_v_center + &
                        ((s% piston_inward_time + s% time)/s% piston_inward_time)* &                        
                           (s% piston_vfinal_inward - s% initial_v_center))
                  else ! moving outward
                     if (s% R_center <= s% piston_Rmin) then
                        s% v_center = s% piston_v0
                     else ! piston slows because of gravitation
                        s% v_center = sqrt(max(0d0, s% piston_v0*s% piston_v0 + &
                           s% piston_alpha*standard_cgrav*s% M_center* &
                              (1/s% R_center - 1/s% piston_Rmin)))
                     end if
                  end if
               end if
               
               if (s% v_center /= 0d0) then ! adjust R_center
            
                  r_nz = exp_cr(s% xh(s% i_lnR,s% nz))
                  if (s% v_flag) then
                     v_nz = s% xh(s% i_v,s% nz)
                  else
                     v_nz = 0
                  end if
                  if (v_nz > 0d0) v_nz = 0.99d0*v_nz ! safety factor
                  dv = s% v_center - v_nz
                  
                  s% R_center = s% R_center_old + s% dt*s% v_center
                  
                  if (s% v_center < 0 .and. s% piston_period <= 0 .and. &
                        s% R_center < s% piston_Rmin) then ! fix it
                     s% R_center = s% piston_Rmin
                     s% v_center = (s% R_center - s% R_center_old)/s% dt
                     write(*,2) 'fix v_center to hit piston_Rmin exactly', &
                        s% model_number, s% v_center, s% R_center
                  end if
                  if (s% R_center < 0) then
                     write(*,2) 's% R_center', s% model_number, s% R_center
                     do_evolve_step_part2 = retry
                     return
                  end if
                                
               end if
            
            end if

            s% time = s% time_old + s% dt

            call debug('before do_adjust_mass', s)

            if (trace) write(*,'(/,a)') 'call do_adjust_mass'         
            if (s% doing_timing) call start_time(s, time0, total)
            call do_adjust_mass(s, s% species, ierr)
            if (s% doing_timing) call update_time(s, time0, total, s% time_adjust_mass)         
            if (ierr /= 0) then
               do_evolve_step_part2 = retry
               s% result_reason = adjust_mass_failed
               if (s% report_ierr) write(*, *) 'do_evolve_step_part2: do_adjust_mass ierr'
               return
            end if

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do j=1,s% species
                  write(*,4) 'evolve after do_adjust_mass xa(j)', &
                     s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
               end do
            end if

            dt = s% dt
            
            call debug('before do_set_vars', s)
            
            if (.not. s% do_element_diffusion) then

               call do_set_vars(0, dt, ierr)
               if (ierr /= 0) return
               do k=1,s% nz ! save in case needed by eps_grav
                  s% lnS_pre(k) = s% lnS(k)
               end do    

            else
            
               call do_set_vars(1, dt, ierr)
               if (ierr /= 0) return
               do k=1,s% nz ! save in case needed by eps_grav
                  s% lnS_pre(k) = s% lnS(k)
               end do    
               
               call debug('before do_element_diffusion', s)
               if (trace) write(*,'(/,a)') 'call do_element_diffusion'
               if (s% doing_timing) call start_time(s, time0, total)
               okay = do_element_diffusion(s, s% dt)
               if (s% doing_timing) &
                  call update_time(s, time0, total, s% time_element_diffusion)         
               if (.not. okay) then
                  if (s% report_ierr) then
                     write(*, *) 'element diffusion failed: retry', s% model_number
                  end if
                  do_evolve_step_part2 = retry
                  s% result_reason = diffusion_failed
                  return
               end if
                  
               call do_set_vars(2, dt, ierr)
               if (ierr /= 0) return
               
            end if
            
            ! get starting energy info for later conservation tests
            nz = s% nz
                        
            call eval_total_energy_integrals(s, &
               s% total_internal_energy_start, &
               s% total_gravitational_energy_start, &
               s% total_linear_kinetic_energy_start, &
               s% total_rotational_kinetic_energy_start, &
               s% total_energy_start)          
            
            if (.not. s% have_initial_energy_integrals) then
               s% total_internal_energy_initial = &
                  s% total_internal_energy_start
               s% total_gravitational_energy_initial = &
                  s% total_gravitational_energy_start
               s% total_linear_kinetic_energy_initial = &
                  s% total_linear_kinetic_energy_start
               s% total_rotational_kinetic_energy_initial = &
                  s% total_rotational_kinetic_energy_start
               s% total_energy_initial = s% total_energy_start
               s% have_initial_energy_integrals = .true.
            end if     
               
            if (s% need_to_reset_cumulative_energy_info) then
               s% need_to_reset_cumulative_energy_info = .false.
               s% cumulative_visc_heat_added = 0
               s% cumulative_eps_grav = 0
               s% cumulative_acoustic_L = 0
               s% cumulative_acoustic_L_center = 0
               s% cumulative_extra_heating = 0
               s% cumulative_irradiation_heating = 0
               s% cumulative_nuclear_heating = 0
               s% cumulative_non_nuc_neu_cooling = 0
               s% cumulative_sources_and_sinks = 0
               s% cumulative_energy_error = 0
               s% cumulative_visc_heat_added_old = 0
               s% cumulative_eps_grav_old = 0
               s% cumulative_acoustic_L_old = 0
               s% cumulative_acoustic_L_center_old = 0
               s% cumulative_extra_heating_old = 0
               s% cumulative_irradiation_heating_old = 0
               s% cumulative_nuclear_heating_old = 0
               s% cumulative_non_nuc_neu_cooling_old = 0
               s% cumulative_sources_and_sinks_old = 0
               s% cumulative_energy_error_old = 0
               s% cumulative_visc_heat_added_older = 0
               s% cumulative_eps_grav_older = 0
               s% cumulative_acoustic_L_older = 0
               s% cumulative_acoustic_L_center_older = 0
               s% cumulative_extra_heating_older = 0
               s% cumulative_irradiation_heating_older = 0
               s% cumulative_nuclear_heating_older = 0
               s% cumulative_non_nuc_neu_cooling_older = 0
               s% cumulative_sources_and_sinks_older = 0
               s% cumulative_energy_error_older = 0
            end if

            s% need_to_adjust_J_lost = .true.
            if (s% rotation_flag .and. s% premix_omega) then
               do_evolve_step_part2 = do_solve_omega_mix(s, 0.5d0*dt)
               if (do_evolve_step_part2 /= keep_going) return
               call set_rotation_info(s, ierr)
               if (ierr /= 0) return
               call set_mixing_info(s, ierr)
               if (ierr /= 0) return            
            end if
            
            call save_start_values(s, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'save_start_values failed'
               return
            end if
            
            call set_Eulerian_Lagrangian_for_eps_grav

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do j=1,s% species
                  write(*,4) 'evolve before do_struct_burn_mix xa(j)', &
                     s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
               end do
            end if

            call debug('before do_struct_burn_mix', s)
            if (trace) write(*,'(/,a,i8)') 'call do_struct_burn_mix', s% model_number            
            skip_global_corr_coeff_limit = (first_try .or. &
                (just_did_backup .and. s% number_of_backups_in_a_row == 1) .or. &
                (s% model_number_for_last_retry /= s% model_number &
                .and. .not. just_did_backup)) ! last alternative is for redo's

            if (s% doing_timing) call start_time(s, time0, total)
            
            do_evolve_step_part2 = do_struct_burn_mix( &
               s, skip_global_corr_coeff_limit, dt)
            if (s% doing_timing) &
               call update_time(s, time0, total, s% time_struct_burn_mix)  
                   
            if (do_evolve_step_part2 /= keep_going) return
            call debug('after do_struct_burn_mix', s)
        
            if (.not. s% rotation_flag) exit implicit_mdot_loop
            if (s% mstar_dot == 0) exit implicit_mdot_loop
            if (max_mdot_redo_cnt <= 0) exit implicit_mdot_loop
                        
            mstar_dot_prev = mstar_dot
            mstar_dot = s% mstar_dot
            wind_mdot = -s% mstar_dot 
               ! for a wind, wind_mdot > 0, but we don't require this in the following.
               ! i.e., the following does not assume that mstar_dot < 0
            
            if (mdot_redo_cnt == 1 .or. ignored_first_step) then 
               ! this is the 1st correction to mdot
               r_phot = sqrt(s% L(1)/(pi*crad*clight*pow4(s% Teff)))
               kh_timescale = eval_kh_timescale(s% cgrav(1), s% mstar, r_phot, s% L(1))
               dmskhf = s% rotational_mdot_kh_fac
               dmsfac = s% rotational_mdot_boost_fac
               max_wind_mdot = dmskhf*s% mstar/kh_timescale
               if (wind_mdot > 0) max_wind_mdot = min(max_wind_mdot, wind_mdot*dmsfac)
            end if
            
            w_div_w_crit_prev = w_div_w_crit
            ! check the new w_div_w_crit to make sure not too large
            call set_surf_avg_rotation_info(s)
            w_div_w_crit = s% w_div_w_crit_avg_surf
            
            !write(*,2) 'log wind_mdot', mdot_redo_cnt, log10(abs(wind_mdot)/(Msun/secyer))
            
            if (wind_mdot >= max_wind_mdot) then
               if (mdot_redo_cnt == 0) then
                  write(*,*) 'cannot fix omega >= omega_crit -- mass loss already at max'
               else
                  write(*,2) 'retry: at max wind mass loss', s% model_number, &
                     log10_cr(max_wind_mdot/(Msun/secyer))
                  do_evolve_step_part2 = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               write(*,*)
               if (w_div_w_crit > surf_w_div_w_crit_limit) then
                  write(*,1) 'retry: w_div_w_crit > surf_w_div_w_crit_limit', &
                     w_div_w_crit, surf_w_div_w_crit_limit
                  do_evolve_step_part2 = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               exit implicit_mdot_loop
            end if
               
            ! NOTE: we assume that if surface omega/omega_crit (w_div_w_crit) is too large,
            ! then mass loss needs to be made larger to fix the problem.
            ! if that assumption is wrong,
            ! i.e. if bigger mass loss makes w_div_w_crit worse,
            ! then in an unstable situation and will remove mass until regain stability.
            
            if (w_div_w_crit <= surf_w_div_w_crit_limit &
                  .and. mdot_redo_cnt == 0) then
               s% was_in_implicit_wind_limit = .false.
               exit implicit_mdot_loop 
            end if
               ! normal case; no problem; no redo required.
            
            if (w_div_w_crit <= surf_w_div_w_crit_limit &
                  .and. s% mstar_dot == explicit_mdot) exit implicit_mdot_loop 
               ! implicit scheme reached the limit setted by the explicit_mdot;
               ! no problem; no redo required.

            s% was_in_implicit_wind_limit = .true.
               
            if (s% dt/secyer < s% min_years_dt_for_redo_mdot) then
               if (.true.) write(*,1) &
                  'dt too small for fix to fix w > w_crit; min_years_dt_for_redo_mdot', &
                  s% dt/secyer, s% min_years_dt_for_redo_mdot
               exit implicit_mdot_loop
            end if
            
            ! if get here, need to revise mdot to fix w_div_w_crit
            
            mdot_redo_cnt = mdot_redo_cnt + 1
            
            if (mdot_redo_cnt == 1) then ! this is the 1st correction to mdot

               call set_current_to_old(s)
               do_evolve_step_part2 = prepare_for_new_try(s)
               if (do_evolve_step_part2 /= keep_going) return
               
               have_too_small_wind_mdot = .true.
               too_small_wind_mdot = wind_mdot
               if (s% mstar_dot < 0) then
                  s% mstar_dot = mstar_dot*s% mdot_revise_factor
               else
                  s% mstar_dot = mstar_dot/s% mdot_revise_factor
               end if
               
               if (-s% mstar_dot > max_wind_mdot) s% mstar_dot = -max_wind_mdot
                  
               write(*,3) 'w > w_crit: revise mdot and redo', &
                  s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))

               !abs_mstar_delta = max(abs(s% mstar_dot), 1d-6*Msun/secyer)
               abs_mstar_delta = abs(s% mstar_dot)
               
               cycle implicit_mdot_loop
               
            else if (mdot_redo_cnt == 2 .and. ignored_first_step) then
               abs_mstar_delta = abs(s% mstar_dot_old)
            end if
            
            ! have already done at least one correction -- check if okay now
            if (w_div_w_crit <= surf_w_div_w_crit_limit .and. &
                  have_too_small_wind_mdot .and. &
                  abs((wind_mdot-too_small_wind_mdot)/wind_mdot) < &
                     s% surf_w_div_w_crit_tol) then
               write(*,3) 'OKAY', s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))
               write(*,*)
               exit implicit_mdot_loop ! in bounds so accept it
            end if

            if (mdot_redo_cnt >= max_mdot_redo_cnt) then
               if (max_mdot_redo_cnt > 0) then
                  write(*,3) 'failed to fix w > w_crit: too many tries', &
                     s% model_number, mdot_redo_cnt, w_div_w_crit, &
                     log10_cr(abs(s% mstar_dot)/(Msun/secyer))
                  do_evolve_step_part2 = retry
                  s% result_reason = nonzero_ierr
                  return
               end if
               exit implicit_mdot_loop
            end if
            
            if (w_div_w_crit > surf_w_div_w_crit_limit &
                  .and. w_div_w_crit_prev >= surf_w_div_w_crit_limit &
                  .and. -mstar_dot >= max_wind_mdot) then
               write(*,3) 'failed to fix w > w_crit', &
                  s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))
               write(*,*)
               do_evolve_step_part2 = retry
               s% result_reason = nonzero_ierr
               return
            end if
            
            if (w_div_w_crit >= surf_w_div_w_crit_limit) then ! wind too small
               !write(*,*) "entering too small wind mdot"
               if (.not. have_too_small_wind_mdot) then
                  !write(*,*) "setting too small wind mdot"
                  too_small_wind_mdot = wind_mdot
                  have_too_small_wind_mdot = .true.
               else if (wind_mdot > too_small_wind_mdot) then
                  !write(*,*) "changing too small wind mdot"
                  too_small_wind_mdot = wind_mdot
               end if
            else if (w_div_w_crit < surf_w_div_w_crit_limit) then ! wind too large
               !write(*,*) "entering too large wind mdot"
               if (.not. have_too_large_wind_mdot) then
                  !write(*,*) "setting too large wind mdot"
                  too_large_wind_mdot = wind_mdot
                  have_too_large_wind_mdot = .true.
               else if (wind_mdot < too_large_wind_mdot) then
                  !write(*,*) "changing too large wind mdot"
                  too_large_wind_mdot = wind_mdot
               end if
            end if

            call set_current_to_old(s)
            do_evolve_step_part2 = prepare_for_new_try(s)
            if (do_evolve_step_part2 /= keep_going) return
            
            if (have_too_large_wind_mdot .and. have_too_small_wind_mdot) then
               if (abs((too_large_wind_mdot-too_small_wind_mdot)/too_large_wind_mdot) &
                   < s% surf_w_div_w_crit_tol) then
                  write(*,*) "too_large_wind_mdot good enough, using it"
                  s% mstar_dot = -too_large_wind_mdot
               else
                  ! have bracketing mdots; bisect for next one.
                  s% mstar_dot = -0.5d0*(too_large_wind_mdot + too_small_wind_mdot)
                  write(*,3) 'fix w > w_crit: bisect mdots and redo', &
                     s% model_number, mdot_redo_cnt, w_div_w_crit, &
                     log10_cr(abs(s% mstar_dot)/(Msun/secyer)), &
                     log10_cr(abs(too_large_wind_mdot)/(Msun/secyer)), &
                     log10_cr(abs(too_small_wind_mdot)/(Msun/secyer))
               end if
                  
            else ! still have wind too small so boost it again
               if (have_too_small_wind_mdot) then
                  if (mod(mdot_redo_cnt,2) == 1) then
                     boost = s% implicit_mdot_boost
                     ! increase mass loss
                     mstar_dot_nxt = mstar_dot - boost*abs_mstar_delta
                  else
                     if (mstar_dot < 0) then ! increase mass loss
                        mstar_dot_nxt = mstar_dot*s% mdot_revise_factor
                     else ! decrease mass gain
                        mstar_dot_nxt = mstar_dot/s% mdot_revise_factor
                     end if
                  end if
               else
                  if (mod(mdot_redo_cnt,2) == 1) then
                     boost = s% implicit_mdot_boost
                     ! decrease mass loss
                     mstar_dot_nxt = mstar_dot + boost*abs_mstar_delta
                  else
                     if (mstar_dot < 0) then ! decrease mass loss
                        mstar_dot_nxt = mstar_dot/s% mdot_revise_factor
                     else ! increase mass gain
                        mstar_dot_nxt = mstar_dot*s% mdot_revise_factor
                     end if
                  end if
               end if
               if (mstar_dot_prev /= explicit_mdot) &
                  mstar_dot_nxt = min(mstar_dot_nxt, explicit_mdot)
               if (mstar_dot_nxt == explicit_mdot) &
                  write(*,*) "implicit mdot: reached explicit_mdot"
               s% mstar_dot = mstar_dot_nxt
               if (-s% mstar_dot > max_wind_mdot) s% mstar_dot = -max_wind_mdot
               !abs_mstar_delta = max(abs_mstar_delta, abs(s% mstar_dot))
               write(*,3) 'fix w > w_crit: change mdot and redo', &
                  s% model_number, mdot_redo_cnt, w_div_w_crit, &
                  log10_cr(abs(s% mstar_dot)/(Msun/secyer))
            end if

         end do implicit_mdot_loop

         call debug('before set_final_vars', s)
         if (trace) write(*,'(/,a)') 'call set_final_vars'
         call set_final_vars(s, s% dt, ierr) ! will use these for next step
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'do_evolve_step_part2: set_final_vars ierr'
            do_evolve_step_part2 = retry
            s% result_reason = nonzero_ierr
            return
         end if
         
         if (s% max_timestep_hi_T_limit > 0 .and. &
               s% max_years_for_timestep /= s% hi_T_max_years_for_timestep) then
            if (maxval(s% T(1:s% nz)) >= s% max_timestep_hi_T_limit) then
               write(*,1) 'switch to high T max timesteps'
               s% max_years_for_timestep = s% hi_T_max_years_for_timestep
               s% max_timestep = secyer*s% max_years_for_timestep
            end if
         end if
         
         call debug('before do_report', s)

         if (trace) write(*,'(/,a)') 'call do_report'
         call system_clock(time0,clock_rate)
         s% current_system_clock_time = time0
         s% total_elapsed_time = &
            dble(time0 - s% starting_system_clock_time)/dble(clock_rate)
         call do_report(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'do_evolve_step_part2: do_report ierr'
            do_evolve_step_part2 = retry
            s% result_reason = nonzero_ierr
            return
         end if
         
         call set_L_burn_by_category(s) ! final values for use in selecting timestep
         
         if (.not. okay_energy_conservation()) return
         
         !if (s% model_number == 11) stop
         
         s% total_angular_momentum = total_angular_momentum(s)
         !write(*,2) 'total_angular_momentum after do_evolve_step_part2', &
         !   s% model_number, s% total_angular_momentum
         
         if (trace) write(*,'(/,a)') 'done do_evolve_step_part2'
         
         call debug('done do_evolve_step_part2', s)
         
         !call show_debug
         
       
         contains
         
         
         subroutine show_debug
            integer :: k
            real(dp) :: alfa, beta, gamma1, Cv, chiRho, chiT, Cp, grada, &
               Pgas, Prad, P, opacity
            include 'formats'
            k = 1205
            
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1 - alfa

            gamma1 = alfa*s% gamma1(k) + beta*s% gamma1(k-1)
            Cv = alfa*s% Cv(k) + beta*s% Cv(k-1)
            chiRho = alfa*s% chiRho(k) + beta*s% chiRho(k-1)
            chiT = alfa*s% chiT(k) + beta*s% chiT(k-1)
            Cp = alfa*s% Cp(k) + beta*s% Cp(k-1)
            grada = alfa*s% grada(k) + beta*s% grada(k-1)
            Pgas = alfa*s% Pgas(k) + beta*s% Pgas(k-1)
            Prad = alfa*s% Prad(k) + beta*s% Prad(k-1)
            P = alfa*s% P(k) + beta*s% P(k-1)
            opacity = alfa*s% opacity(k) + beta*s% opacity(k-1)
            
            write(*,2) 'at end of step', s% model_number
            write(*,2) 'gamma1', k, gamma1
            write(*,2) 'Cv', k, Cv
            write(*,2) 'chiRho', k, chiRho
            write(*,2) 'chiT', k, chiT
            write(*,2) 'Cp', k, Cp
            write(*,2) 'grada', k, grada
            write(*,2) 'Pgas', k, Pgas
            write(*,2) 'Prad', k, Prad
            write(*,2) 'P', k, P
            write(*,2) 'opacity', k, opacity
            write(*,2) 'L', k, s% L(k)
            write(*,2) 'gradr', k, s% gradr(k)
            write(*,2) 'gradr/grada', k, s% gradr(k)/grada
            write(*,3) 'mixing_type', k, s% mixing_type(k)
            write(*,*)
         
         end subroutine show_debug
         
         
         logical function okay_energy_conservation()
            use star_utils, only: set_acoustic_L
            integer :: nz
            include 'formats'
            
            okay_energy_conservation = .false.
         
            nz = s% nz
         
            call eval_total_energy_integrals(s, &
               s% total_internal_energy, &
               s% total_gravitational_energy, &
               s% total_linear_kinetic_energy, &
               s% total_rotational_kinetic_energy, &
               s% total_energy)

            s% total_eps_grav = dt*s% eps_grav_dot_product_dm_average
               
            if (s% v_flag .and. s% use_artificial_viscosity) then
               s% total_visc_heat_added = &
                  s% eps_visc_factor*dt*s% eps_visc_dot_product_dm_average
            else
               s% total_visc_heat_added = 0
            end if

            s% cumulative_eps_grav = &
               s% cumulative_eps_grav_old + s% total_eps_grav

            s% cumulative_visc_heat_added = &
               s% cumulative_visc_heat_added_old + s% total_visc_heat_added
            
            call set_acoustic_L(s)

            s% cumulative_acoustic_L = &
               s% cumulative_acoustic_L_old + dt*s% acoustic_L_average
            s% cumulative_acoustic_L_center = &
               s% cumulative_acoustic_L_center_old + dt*s% acoustic_L_center
         
            if (s% rotation_flag .and. &
                  (s% use_other_torque .or. s% use_other_torque_implicit .or. &
                     associated(s% binary_other_torque))) then
               ! keep track of rotational kinetic energy
            end if
               
            s% total_nuclear_heating = dt*s% eps_nuc_dot_product_dm_average
            s% cumulative_nuclear_heating = &
               s% cumulative_nuclear_heating_old + s% total_nuclear_heating

            s% total_non_nuc_neu_cooling = dt*0.5d0*( &
               s% non_nuc_neu_dot_product_dm_average + &
               sum(s% non_nuc_neu_start(1:nz)*s% dm(1:nz)))
            s% cumulative_non_nuc_neu_cooling = &
               s% cumulative_non_nuc_neu_cooling_old + s% total_non_nuc_neu_cooling

            s% total_irradiation_heating = dt*s% irradiation_heat_dot_product_dm_average
            s% cumulative_irradiation_heating = &
               s% cumulative_irradiation_heating_old + s% total_irradiation_heating

            s% total_extra_heating = dt*s% extra_heat_dot_product_dm_average
            s% cumulative_extra_heating = &
               s% cumulative_extra_heating_old + s% total_extra_heating

            s% cumulative_L_center = &
               s% cumulative_L_center_old + dt*s% L_center
            
            s% total_energy_sources_and_sinks = &
                 s% total_nuclear_heating &
               - s% total_non_nuc_neu_cooling &
               + s% total_irradiation_heating &
               + s% total_extra_heating &
               + dt*(s% acoustic_L_center - s% acoustic_L_average) &
               + dt*s% L_center

            if (s% L_flag) then
               s% cumulative_L_surf = &
                  s% cumulative_L_surf_old + dt*s% surface_L_average
               s% total_energy_sources_and_sinks = &
                  s% total_energy_sources_and_sinks - dt*s% surface_L_average
            end if

            s% cumulative_sources_and_sinks = &
               s% cumulative_sources_and_sinks_old + s% total_energy_sources_and_sinks
         
            s% error_in_energy_conservation = &
               ! positive error means ended up with too much energy in the model.
               (s% total_energy - s% total_energy_start) - s% total_energy_sources_and_sinks
               
            s% cumulative_energy_error = s% cumulative_energy_error_old + &
               s% error_in_energy_conservation
               
               
            
            
            if (.false.) then
            
            

            write(*,1) 'set total_internal_energy', s% total_internal_energy
            write(*,1) 'set total_gravitational_energy', s% total_gravitational_energy
            write(*,1) 'set total_linear_kinetic_energy', s% total_linear_kinetic_energy
            write(*,1) 'set total_rotational_kinetic_energy', s% total_rotational_kinetic_energy
            write(*,1) 'set total_energy', s% total_energy
               
            write(*,1) 'step delta total_internal_energy', s% total_internal_energy - &
               s% total_internal_energy_start
            write(*,1) 'step delta total_gravitational_energy', s% total_gravitational_energy - &
               s% total_gravitational_energy_start
            write(*,1) 'step delta total_linear_kinetic_energy', s% total_linear_kinetic_energy - &
               s% total_linear_kinetic_energy_start
            write(*,1) 'step delta total_rotational_kinetic_energy', s% total_rotational_kinetic_energy - &
               s% total_rotational_kinetic_energy_start
            write(*,1) 'step delta total_energy', s% total_energy - &
               s% total_energy_start
            
            write(*,2) 'sum dEdt_expected', &
               s% model_number, &
               sum(s% d_IE_dt_expected(1:nz) + &
                     s% d_KE_dt_expected(1:nz) + &
                     s% d_PE_dt_expected(1:nz))
            write(*,2) 'sum dEdt_actual', &
               s% model_number, &
               sum(s% d_IE_dt_actual(1:nz) + &
                     s% d_KE_dt_actual(1:nz) + &
                     s% d_PE_dt_actual(1:nz))
            write(*,2) 'sum sources_and_sinks/dt', &
               s% model_number, &
               s% total_energy_sources_and_sinks/dt
            write(*,*)
            write(*,2) 'dt*sum dEdt_actual', &
               s% model_number, &
               dt*sum(s% d_IE_dt_actual(1:nz) + &
                     s% d_KE_dt_actual(1:nz) + &
                     s% d_PE_dt_actual(1:nz))
            write(*,2) 'total_energy_sources_and_sinks', &
               s% model_number, &
               s% total_energy_sources_and_sinks
            write(*,2) 'dt*sum dEdt_actual - total_energy_sources_and_sinks', &
               s% model_number, &
               dt*sum(s% d_IE_dt_actual(1:nz) + &
                     s% d_KE_dt_actual(1:nz) + &
                     s% d_PE_dt_actual(1:nz)) - s% total_energy_sources_and_sinks
            write(*,*)
            write(*,2) 'total_internal_energy_start', &
               s% model_number, s% total_internal_energy_start
            write(*,2) 'total_internal_energy', &
               s% model_number, s% total_internal_energy
            write(*,*)
            write(*,2) 'total_internal_energy_change', &
               s% model_number, s% total_internal_energy - s% total_internal_energy_start
            write(*,*)
            
            write(*,2) 's% total_nuclear_heating', &
               s% model_number, s% total_nuclear_heating
            write(*,2) 's% total_non_nuc_neu_cooling', &
               s% model_number, s% total_non_nuc_neu_cooling
            write(*,2) 's% total_irradiation_heating', &
               s% model_number, s% total_irradiation_heating
            write(*,2) 's% total_extra_heating', &
               s% model_number, s% total_extra_heating
            write(*,2) 's% acoustic_L_center', &
               s% model_number, s% extra_heat_dot_product_dm_average
            write(*,2) 's% extra_heat_dot_product_dm_average', &
               s% model_number, s% acoustic_L_center
            write(*,2) 's% acoustic_L_average', &
               s% model_number, s% acoustic_L_average
            write(*,2) 'dt', &
               s% model_number, dt
            write(*,*)
            
            write(*,2) 's% total_visc_heat_added', &
               s% model_number, s% total_visc_heat_added
            write(*,2) 'total_kinetic_energy_change', &
               s% model_number, s% total_linear_kinetic_energy - s% total_linear_kinetic_energy_start
            write(*,2) 'total_internal_energy_change', &
               s% model_number, s% total_internal_energy - s% total_internal_energy_start
            write(*,2) 'total_gravitational_energy_change', &
               s% model_number, s% total_gravitational_energy - s% total_gravitational_energy_start
            write(*,2) 's% total_energy - s% total_energy_start', &
               s% model_number, s% total_energy - s% total_energy_start
            write(*,2) 's% total_energy_sources_and_sinks', &
               s% model_number, s% total_energy_sources_and_sinks
            write(*,2) 's% error_in_energy_conservation', &
               s% model_number, s% error_in_energy_conservation
            write(*,*)

            if (abs(s% error_in_energy_conservation/s% total_energy) > 1d0) then
               write(*,2) 'error in energy', s% model_number, &
                  s% error_in_energy_conservation/s% total_energy, &
                  s% error_in_energy_conservation, s% total_energy
               stop 'evolve'
            end if
            
            stop 'evolve'

            end if
         
            okay_energy_conservation = .true.
         
         end function okay_energy_conservation

         
         subroutine do_set_vars(which_case, dt, ierr)
            use hydro_vars, only: set_vars, set_vars_before_diffusion, &
               reset_vars_after_diffusion
            integer, intent(in) :: which_case
            real(dp), intent(in) :: dt
            integer, intent(out) :: ierr
            integer :: nz, k
            include 'formats'
            nz = s% nz
            do k=1,nz
               s% T_start(k) = -1d99
               s% r_start(k) = -1d99
               s% v_start(k) = -1d99
               s% csound_start(k) = -1d99
            end do
            ierr = 0
            if (trace) write(*,'(/,a)') 'call set_vars'
         
            select case(which_case)
            case (0)
               call set_vars(s, dt, ierr)
            case (1)
               call set_vars_before_diffusion(s, dt, ierr)
            case (2)
               call reset_vars_after_diffusion(s, dt, ierr)
            end select
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*, *) 'do_evolve_step_part2: set_vars ierr: retry', s% model_number
               end if
               do_evolve_step_part2 = retry
               s% result_reason = nonzero_ierr
               return
            end if
         end subroutine do_set_vars
         
         
         subroutine set_Eulerian_Lagrangian_for_eps_grav
            use star_utils, only: set_k_CpTMdot_lt_L
            real(dp) :: dxm_CpTMdot_lt_L, dxm_kA, dxm_target
            integer :: k, kA, kB, nz, transition
            
            include 'formats'
            
            if (s% mstar_dot == 0d0) then
               ! all Lagrangian
               s% k_below_Eulerian_eps_grav = 1
               s% k_Lagrangian_eps_grav = 1
            else if (s% mstar_dot > 0) then
               ! lock Ds/dt form to match type of grid at each point
               ! limits are handled in grid selection
               s% k_below_Eulerian_eps_grav = s% k_below_const_q ! pure "Eulerian" for k < this
               s% k_Lagrangian_eps_grav = s% k_const_mass ! pure Lagrangian for k >= this
            else !  ( s% mstar_dot < 0 )
               ! might want to match this to grid in future
               nz = s% nz
               kA = 1
               dxm_kA = 0

               if( s% min_dxm_Eulerian_div_dxm_removed > 0) then
                  dxm_target = (s% xmstar_old - s% xmstar)*s% min_dxm_Eulerian_div_dxm_removed
                  do k = kA, nz
                     if (dxm_kA >= dxm_target) exit
                     kA = k
                     dxm_kA = dxm_kA + s% dm(k)
                  end do
               end if
            
               ! ensure some consistency
               if (kA == 1) then ! pure Lagrangian
                  kB = 1
               else if (kA >= nz) then ! pure Eulerian
                  kA = nz+1
                  kB = kA
               else ! transition zone between
                  transition = max(0, s% min_cells_for_Eulerian_to_Lagrangian_transition)
                  kB = kA
                  kA = kB-transition
                  if (kA < 1) then
                     kA = 1
                     kB = kA+transition
                  end if
               end if

               s% k_below_Eulerian_eps_grav = kA ! pure Eulerian for k < this
               s% k_Lagrangian_eps_grav = kB ! pure Lagrangian for k >= this
            end if
             
            !  the CpTMdot control has been superceded and should be removed
            !if (s% min_dxm_Eulerian_div_dxm_CpTMdot_lt_L > 0) then
            !   call set_k_CpTMdot_lt_L(s)
            !   dxm_CpTMdot_lt_L = sum(s% dm(1:s% k_CpTMdot_lt_L))
            !   dxm_target = s% min_dxm_Eulerian_div_dxm_CpTMdot_lt_L*dxm_CpTMdot_lt_L
            !   do k = kA, nz
            !      if (dxm_kA >= dxm_target) exit
            !      kA = k
            !      dxm_kA = dxm_kA + s% dm(k)
            !   end do
            !end if
           
     
         end subroutine set_Eulerian_Lagrangian_for_eps_grav  


      end function do_evolve_step_part2
         
      subroutine debug(str, s) 
         use chem_def
         character (len=*), intent(in) :: str
         type(star_info), pointer :: s
         integer :: k, j
         include 'formats'
         
         !write(*,2) trim(str) // ' s% xh(s% i_lnd,1)', s% model_number, s% xh(s% i_lnd,1)

         !write(*,2) 's% xtra_coef_czb_full_on', s% model_number, s% xtra_coef_czb_full_on
         !write(*,2) 's% xtra_coef_czb_full_off', s% model_number, s% xtra_coef_czb_full_off
         return
         
         if (.not. s% rotation_flag) return
         k = 1
         write(*,3) trim(str) // ' s% omega(k)', k, s% model_number, s% omega(k)
         return
         j = 2
         !do j=1,1 !s% species
            if (.true. .or. s% xa(j,k) > 1d-9) &
               write(*,1) trim(str) // ' xin(net_iso(i' // &
                  trim(chem_isos% name(s% chem_id(j))) // '))= ', &
                  s% xa(j,k), s% abar(k)
         !end do
      end subroutine debug
                     

      subroutine save_start_values(s, ierr)
         use solve_hydro, only: set_L_burn_by_category
         use chem_def, only: num_categories
         use star_utils, only: get_Lrad
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, j
         real(dp), pointer :: p1(:), p2(:,:)
         integer, pointer :: ip(:)
         real(dp) :: L_rad
         include 'formats'    
         ierr = 0         
         call set_L_burn_by_category(s)
         do k=1,s% nz
            !s% T_start(k) set elsewhere
            !s% r_start(k) set elsewhere
            !s% v_start(k) set elsewhere
            !s% csound_start(k) set elsewhere
            s% lnd_start(k) = s% lnd(k)
            s% lnP_start(k) = s% lnP(k)
            s% P_start(k) = s% P(k)
            s% lnPgas_start(k) = s% lnPgas(k)
            s% lnT_start(k) = s% lnT(k)
            s% lnE_start(k) = s% lnE(k)
            s% energy_start(k) = s% energy(k)
            s% lnR_start(k) = s% lnR(k)
            s% v_start(k) = s% v(k)
            s% L_start(k) = s% L(k)
            s% omega_start(k) = s% omega(k)
            s% ye_start(k) = s% ye(k)
            s% Z_start(k) = min(1d0, max(0d0, 1d0 - (s% X(k) + s% Y(k))))
            s% i_rot_start(k) = s% i_rot(k)
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% non_nuc_neu_start(k) = s% non_nuc_neu(k)
            s% mass_correction_start(k) = s% mass_correction(k)
            s% P_div_rho_start(k) = s% P(k)/s% rho(k)
            do j=1,s% species
               s% dxdt_nuc_start(j,k) = s% dxdt_nuc(j,k)
            end do
            do j=1,num_categories
               s% luminosity_by_category_start(j,k) = &
                  s% luminosity_by_category(j,k)
            end do
            s% scale_height_start(k) = s% scale_height(k)
            s% gradL_start(k) = s% gradL(k)
            s% grada_start(k) = s% grada(k)
            s% gradr_start(k) = s% gradr(k)
            s% grada_at_face_start(k) = s% grada_at_face(k)
            s% chiT_start(k) = s% chiT(k)
            s% chiRho_start(k) = s% chiRho(k)
            s% cp_start(k) = s% cp(k)
            s% cv_start(k) = s% cv(k)
            s% dE_dRho_start(k) = s% dE_dRho_start(k)
            s% gam_start(k) = s% gam(k)
            s% entropy_start(k) = s% entropy(k)
            s% rho_start(k) = s% rho(k)
            s% eta_start(k) = s% eta(k)
            s% abar_start(k) = s% abar(k)
            s% zbar_start(k) = s% zbar(k)
            s% mu_start(k) = s% mu(k)
            s% eps_nuc_start(k) = s% eps_nuc(k)
            s% opacity_start(k) = s% opacity(k)            
            s% mlt_mixing_length_start(k) = s% mlt_mixing_length(k)
            s% mlt_mixing_type_start(k) = s% mlt_mixing_type(k)
            s% mlt_D_start(k) = s% mlt_D(k)
            s% mlt_vc_start(k) = s% mlt_vc(k)
            s% mlt_Gamma_start(k) = s% mlt_Gamma(k)
            s% mlt_cdc_start(k) = s% mlt_cdc(k)            
            s% burn_num_iters(k) = 0
            s% burn_num_substeps(k) = 0
            if (s% qmin_freeze_non_radiative_luminosity < 1d0) then
               L_rad = get_Lrad(s,k)
               s% L_non_rad_start(k) = s% L(k) - L_rad
            end if
         end do
         s% have_start_values = .true.         
      end subroutine save_start_values

      
      integer function prepare_for_new_step(s)
         use utils_lib, only: has_bad_num, is_bad_num
         use evolve_support, only: new_generation
         use chem_def

         type (star_info), pointer :: s
         
         integer :: ierr, nz, k
         real(dp) :: total_energy, piston_done_time
         logical :: trace
         
         include 'formats'

         ierr = 0
         trace = s% trace_evolve
         nz = s% nz

         prepare_for_new_step = keep_going

         if (s% dt_next <= 0) then
            prepare_for_new_step = terminate
            if ((s% time >= s% max_age*secyer .and. s% max_age > 0) .or. &
                (s% time >= s% max_age_in_seconds .and. s% max_age_in_seconds > 0)) then
               s% result_reason = result_reason_normal
               s% termination_code = t_max_age
            else
               s% result_reason = dt_is_zero
               s% termination_code = t_dt_is_zero
            end if
            return
         end if

         if (s% dt_next < s% min_timestep_limit) then
            write(*, *) 's% dt_next', s% dt_next
            write(*, *) 's% min_timestep_limit', s% min_timestep_limit
            prepare_for_new_step = terminate
            s% termination_code = t_min_timestep_limit
            s% result_reason = timestep_limits
            return
         end if
         
         if (s% stop_when_done_with_center_flash .and. &
               s% done_with_center_flash) then
            prepare_for_new_step = terminate
            s% termination_code = t_done_with_center_flash
            return
         end if

         if (s% doing_center_flash .and. &
               s% time >= s% center_flash_total_time) then
            s% L_center = s% initial_L_center
            s% doing_center_flash = .false.
            write(*,*) 'done with center flash'
         end if
         
         if (s% use_piston .and. .not. s% done_with_piston) then
            if (s% piston_period > 0) then ! periodic piston
               piston_done_time = s% piston_period* &
                  (s% periodic_piston_number_of_cycles + s% periodic_piston_done_delay)
               if (s% time >= piston_done_time) then
                  s% v_center = 0d0        
                  s% done_with_piston = .true.
                  write(*,*)
                  write(*,2) 'turn off piston', s% model_number
                  write(*,*)
               end if
            else if (s% R_center >= s% piston_Rmax) then
               s% v_center = 0d0        
               s% done_with_piston = .true.          
               write(*,*)
               write(*,2) 'turn off piston', s% model_number, s% R_center, &
                  s% piston_Rmax
               write(*,*)
            end if
            if (s% done_with_piston .and. s% stop_when_done_with_piston) then
               prepare_for_new_step = terminate
               s% termination_code = t_done_with_piston
               return
            end if
            if (s% done_with_piston .and. &
                  s% reset_total_energy_initial_when_done_with_piston) then
               s% need_to_reset_cumulative_energy_info = .true.
            end if
         end if
            
         if (trace) write(*,*) 'call set_start_of_step_info'
         call set_start_of_step_info(s, ierr)
         if (failed('set_start_of_step_info ierr')) return

         if (.not. s% doing_first_model_of_run) then
            if (trace) write(*,*) 'call do_mesh'
            prepare_for_new_step = do_mesh(s)
            if (prepare_for_new_step /= keep_going) then
               write(*,*) 'failed in do_mesh'
               prepare_for_new_step = terminate
               return
            end if
            if (s% ebdf_order == s% max_ebdf_order) then
               s% startup_increment_ebdf_order = .false.
            else if (s% ebdf_order < s% max_ebdf_order .and. &
                  (s% startup_increment_ebdf_order .or. &
                   s% model_number >= s% ebdf_hold)) then
               s% ebdf_order = s% ebdf_order + 1
               write(*,3) 'increase edbf order', s% ebdf_order, s% model_number
               s% ebdf_hold = s% model_number + &
                  max(1, s% steps_before_try_higher_ebdf_order)
            end if
         end if

         if (trace) write(*,*) 'call new_generation'
         call new_generation(s, ierr)         
         if (failed('new_generation ierr')) return

         s% dt = s% dt_next 
         s% dt_start = s% dt         
         s% retry_cnt = 0
         s% redo_cnt = 0
         s% generations = min(max_generations, s% generations+1)
         
         s% need_to_save_profiles_now = .false.
         s% need_to_update_history_now = s% doing_first_model_of_run
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            failed = .true.
            prepare_for_new_step = terminate
            if (s% report_ierr) write(*, *) 'prepare_for_new_step: ' // trim(str)
            s% result_reason = nonzero_ierr
         end function failed


      end function prepare_for_new_step


      integer function do_mesh(s)
         use adjust_mesh, only: remesh
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s      
         logical, parameter :: okay_to_merge = .true.
         integer :: time0, clock_rate
         real(dp) :: total
         include 'formats'
         if (s% doing_timing) call start_time(s, time0, total)
         do_mesh = remesh(s, okay_to_merge)         
         if (s% doing_timing) call update_time(s, time0, total, s% time_remesh)         
         if (do_mesh /= keep_going) then
            s% result_reason = adjust_mesh_failed
            if (s% report_ierr) write(*, *) 'do_mesh: remesh failed'
            return
         end if
         s% have_start_values = .false.
      end function do_mesh
      

      integer function prepare_for_new_try(s) 
         ! return keep_going, terminate, retry, or backup
         ! if don't return keep_going, then set result_reason to say why.
         use utils_lib, only: is_bad_num, has_bad_num
         use net_lib, only: clean_up_fractions
         use net, only: get_screening_mode
         use star_utils, only: save_for_d_dt

         type (star_info), pointer :: s
         
         integer :: ierr, i, j, k, nz, nvar, nvar_hydro
         real(dp), parameter :: max_sum_abs = 10d0, xsum_tol = 1d-2
         real(dp) :: r00, r003, rp13, rm13, r_in, r_out, screening 
         logical :: okay

         include 'formats'
         
         ierr = 0
         nvar = s% nvar
         nvar_hydro = s% nvar_hydro
         nz = s% nz

         s% result_reason = result_reason_normal
         s% model_number = s% model_number_old + 1
         s% termination_code = 0
         s% newton_iter = 0
         s% newton_adjust_iter = 0
         
         screening = get_screening_mode(s,ierr)
         if (ierr /= 0) then
            write(*,*) 'bad value for screening_mode ' // trim(s% screening_mode)
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         
         ! check dimensions
         if (size(s% xh_old,dim=1) /= nvar_hydro .or. size(s% xh_old,dim=2) < nz) then
            write(*,*) 'bad dimensions for xh_old', size(s% xh_old,dim=1), nvar_hydro, &
               size(s% xh_old,dim=2), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         if (size(s% xa_old,dim=1) /= s% species .or. size(s% xa_old,dim=2) < nz) then
            write(*,*) 'bad dimensions for xa_old', size(s% xa_old,dim=1), s% species, &
               size(s% xa_old,dim=2), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         if (size(s% q_old,dim=1) < nz) then
            write(*,*) 'bad dimensions for q_old', size(s% q_old,dim=1), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         if (size(s% dq_old,dim=1) < nz) then
            write(*,*) 'bad dimensions for dq_old', size(s% dq_old,dim=1), nz
            prepare_for_new_try = terminate
            s% termination_code = t_failed_prepare_for_new_try
            return
         end if
         
         do k = 1, nz
            do j=1,nvar_hydro
               s% xh(j,k) = s% xh_old(j,k) ! start from copy of old structure
            end do
            do j=1,s% species
               s% xa(j,k) = s% xa_old(j,k) ! start from copy of old composition
            end do
            s% q(k) = s% q_old(k) ! start with same q's
            s% dq(k) = s% dq_old(k) ! start with same dq's
         end do

         call set_m_and_dm(s)
         call set_dm_bar(nz, s% dm, s% dm_bar)
         
         if (s% rotation_flag) then
            okay = .true.
            do k=1,nz
               s% omega(k) = s% omega_old(k)
               if (is_bad_num(s% omega(k)) .or. &
                     s% omega(k) > 1d50 .or. &
                     (s% omega_old(nz) /= 0 .and. s% omega(k) < 1d-50)) then
                  write(*,2) 's% omega(k)', k, s% omega(k)
                  okay = .false.
               end if
            end do
            if (.not. okay) then
               write(*,2) 'model_number', s% model_number
               stop 'prepare_for_new_try'
            end if
            call use_xh_to_update_i_rot_and_j_rot(s)
            s% total_angular_momentum = total_angular_momentum(s)
            !write(*,2) 'total_angular_momentum after use_xh_to_update_i_rot_and_j_rot', &
            !   s% model_number, s% total_angular_momentum
            if (s% total_angular_momentum < 0) then
               write(*,*) 'ERROR: doing rotation, but total_angular_momentum < 0'
               prepare_for_new_try = terminate
               s% termination_code = t_negative_total_angular_momentum
               return
               stop 'prepare_for_new_try'
            end if
         end if
         
         if (s% just_did_backup) then
            !write(*,*) 'just_did_backup call save_for_d_dt'
            call save_for_d_dt(s)
            !write(*,*) 'just_did_backup call set_start_of_step_info'
            call set_start_of_step_info(s, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) &
                  write(*, *) 'prepare_for_new_try: set_start_of_step_info ierr'
               s% result_reason = nonzero_ierr
               prepare_for_new_try = retry; return
            end if
         end if
         
         prepare_for_new_try = keep_going
         
         
         !write(*,2) 'done prepare_for_new_try omega(1)', s% model_number, s% omega(1)
         
         
      end function prepare_for_new_try
      
      
      integer function pick_next_timestep(id)
         ! determine what we want for the next timestep
         ! if don't return keep_going, then set result_reason to say why.
         use timestep, only: timestep_controller
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: i, j, n
         real(dp) :: max_timestep, remaining_years, min_max, prev_max_years
         include 'formats'
         
         pick_next_timestep = terminate         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% trace_evolve) write(*,'(/,a)') 'pick_next_timestep'

         if (s% max_years_for_timestep > 0) then
            max_timestep = secyer*s% max_years_for_timestep
         else
            max_timestep = s% max_timestep
         end if
         pick_next_timestep = timestep_controller(s, max_timestep)
         if (pick_next_timestep /= keep_going) then
            if (s% trace_evolve) &
               write(*,*) 'pick_next_timestep: timestep_controller /= keep_going'
            return
         end if
         s% dt_next_unclipped = s% dt_next 
               ! write out the unclipped timestep in saved models
         if (s% time < 0 .and. s% time + s% dt_next > 0) then
            s% dt_next = -s% time            
         else if ((s% time + s% dt_next) > s% max_age*secyer .and. s% max_age > 0) then
            s% dt_next = max(0d0, s% max_age*secyer - s% time)            
         else if ((s% time + s% dt_next) > s% max_age_in_seconds &
                  .and. s% max_age_in_seconds > 0) then
            s% dt_next = max(0d0, s% max_age_in_seconds - s% time)            
         else if (s% num_adjusted_dt_steps_before_max_age > 0 .and. &
                  s% max_years_for_timestep > 0) then
            if (s% max_age > 0) then
               remaining_years = s% max_age - s% star_age
            else if (s% max_age_in_seconds > 0) then
               remaining_years = (s% max_age_in_seconds - s% time)/secyer
            else
               remaining_years = 1d99
            end if
            if (s% using_revised_max_yr_dt) &
               s% max_years_for_timestep = s% revised_max_yr_dt
            n = floor(remaining_years/s% max_years_for_timestep + 1d-6)
            j = s% num_adjusted_dt_steps_before_max_age
            if (remaining_years <= s% max_years_for_timestep) then
               s% max_years_for_timestep = remaining_years
               s% using_revised_max_yr_dt = .true.
               s% revised_max_yr_dt = s% max_years_for_timestep
               s% dt_next = s% max_years_for_timestep*secyer
               write(*,3) 'remaining steps and years until max age', &
                  s% model_number, 1, remaining_years
            else if (n <= j) then
               prev_max_years = s% max_years_for_timestep
               i = floor(remaining_years/s% dt_years_for_steps_before_max_age + 1d-6)
               if ((i+1d-9)*s% dt_years_for_steps_before_max_age < remaining_years) then
                  s% max_years_for_timestep = remaining_years/(i+1)
               else
                  s% max_years_for_timestep = remaining_years/i
               end if
               min_max = prev_max_years*s% reduction_factor_for_max_timestep
               if (s% max_years_for_timestep < min_max) &
                  s% max_years_for_timestep = min_max
               if (.not. s% using_revised_max_yr_dt) then
                  s% using_revised_max_yr_dt = .true.
                  write(*,2) 'begin reducing max timestep prior to max age', &
                     s% model_number, remaining_years
               else if (s% revised_max_yr_dt > s% max_years_for_timestep) then
                  write(*,2) 'reducing max timestep prior to max age', &
                     s% model_number, remaining_years
               else if (s% max_years_for_timestep <= &
                     s% dt_years_for_steps_before_max_age) then
                  i = floor(remaining_years/s% max_years_for_timestep + 1d-6)
                  write(*,3) 'remaining steps and years until max age', &
                     s% model_number, i, remaining_years
               else 
                  write(*,2) 'remaining_years until max age', &
                     s% model_number, remaining_years
               end if
               s% revised_max_yr_dt = s% max_years_for_timestep
               if (s% dt_next/secyer > s% max_years_for_timestep) &
                  s% dt_next = s% max_years_for_timestep*secyer
            end if
            
         end if
         
      end function pick_next_timestep
      
      
      integer function prepare_to_redo(id)
         use evolve_support, only: set_current_to_old
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr
         include 'formats'        
         ierr = 0 
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            prepare_to_redo = terminate
            return
         end if         
                  
         s% redo_cnt = s% redo_cnt + 1
         if (s% redo_limit > 0 .and. s% redo_cnt > s% redo_limit) then
            s% dt_start = sqrt(s% dt*s% dt_start)
            prepare_to_redo = backup
            write(*,*) 'have reached redo limit so now backup'
            !stop
            return
         end if
         prepare_to_redo = keep_going         
         if (s% trace_evolve) write(*,'(/,a)') 'prepare_to_redo'        
         call set_current_to_old(s)             
      end function prepare_to_redo
      
      
      integer function prepare_to_retry(id)
         use evolve_support, only: set_current_to_old
         integer, intent(in) :: id
         
         real(dp) :: retry_factor
         type (star_info), pointer :: s
         integer :: ierr, k
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            prepare_to_retry = terminate
            return
         end if
         
         if (s% trace_evolve) write(*,'(/,a)') 'prepare_to_retry'
                  
         s% retry_cnt = s% retry_cnt + 1
         if (s% retry_limit > 0 .and. s% retry_cnt > s% retry_limit) then
            s% dt_start = sqrt(s% dt*s% dt_start)
            prepare_to_retry = backup
            write(*,*) 'have reached retry limit so now backup'
            !stop
            return
         end if

         prepare_to_retry = keep_going
         
         retry_factor = s% timestep_factor_for_retries
         s% dt = s% dt*retry_factor
         if (len_trim(s% retry_message) > 0) write(*,'(a)') trim(s% retry_message)
         write(*,'(a50,2i6,2f16.6)') 'retry cnt, step, log10(dt/yr), retry_factor', &
            s% retry_cnt, s% model_number, &
            log10_cr(s% dt*retry_factor/secyer), retry_factor
         if (s% dt <= max(s% min_timestep_limit,0d0)) then
            write(*,1) 'dt', s% dt
            write(*,1) 'min_timestep_limit', s% min_timestep_limit
            call report_convergence_problems(s, 'dt < min_timestep_limit')
            prepare_to_retry = terminate
            s% termination_code = t_min_timestep_limit
            s% result_reason = timestep_limits
            return
         end if

         if (s% max_years_for_timestep > 0) &
            s% max_timestep = secyer*s% max_years_for_timestep
         if (s% max_timestep > 0) s% dt = min(s% dt, s% max_timestep)
         
         call set_current_to_old(s)
         s% num_retries = s% num_retries+1
         !write(*,2) 'prepare_to_retry s% num_retries', s% num_retries
         !write(*,2) 's% max_number_retries', s% max_number_retries
         if (s% num_retries > s% max_number_retries .and. s% max_number_retries >= 0) then
            write(*,2) 'num_retries', s% num_retries
            write(*,2) 'max_number_retries', s% max_number_retries
            call report_convergence_problems(s, '-- too many retries')
            s% termination_code = t_max_number_retries
            prepare_to_retry = terminate; return
         end if

         s% model_number_for_last_retry = s% model_number
         if (s% why_Tlim == Tlim_neg_X) &
            s% timestep_hold = s% model_number + s% neg_mass_fraction_hold         
            
      end function prepare_to_retry
      
      
      subroutine report_convergence_problems(s,str)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: str
         write(*,*)
         write(*,*) 'stopping because of convergence problems ' // trim(str)
         write(*,*)
      end subroutine report_convergence_problems
      
      
      integer function do1_backup(id)
         ! returns either keep_going or terminate
         ! at end of this, must have same vars set as at end of prepare_for_new_step.
         
         use evolve_support, only: restore_older, set_current_to_old
         use alloc, only: free_star_info_arrays, allocate_star_info_arrays, &
            set_var_info, set_chem_names
         
         integer, intent(in) :: id
         
         integer :: ierr
         type (star_info), pointer :: s
         
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            s% result_reason = nonzero_ierr
            do1_backup = terminate
            return
         end if
         
         if (s% trace_evolve) write(*,'(/,a)') 'do1_backup'
         
         do1_backup = keep_going
         
         s% dt = s% dt_start * &
            (s% timestep_factor_for_backups**(1+2*s% number_of_backups_in_a_row))         
         if (s% doing_first_model_of_run) s% dt = 0.01d0*s% dt

         if (s% dt < s% min_timestep_limit) then
            write(*,1) 'dt', s% dt
            write(*,1) 'min_timestep_limit', s% min_timestep_limit
            call report_convergence_problems(s, 'dt < min_timestep_limit')
            do1_backup = terminate
            s% termination_code = t_min_timestep_limit
            s% result_reason = timestep_limits
            return
         end if
         
         if (s% generations > 2) then        
            !write(*,2) 'do1_backup omega_older(1)', s% model_number, s% omega_older(1)
            !write(*,2) 'do1_backup omega_old(1)', s% model_number, s% omega_old(1)
            call restore_older(s) ! set old = older
            s% generations = s% generations-1        
            call set_var_info(s, ierr)
            if (ierr /= 0) then
               write(*, *) 'do1_backup: set_var_info ierr', ierr
               s% result_reason = nonzero_ierr
               do1_backup = terminate
               return
            end if            
            call set_chem_names(s)            
            call set_current_to_old(s)
            call free_star_info_arrays(s)
            call allocate_star_info_arrays(s, ierr)
            if (ierr /= 0) then
               write(*, *) 'do1_backup: allocate_star_info_arrays ierr'
               s% result_reason = nonzero_ierr
               do1_backup = terminate
               return
            end if            
         else
            call set_current_to_old(s)         
         end if
      
         s% dt_limit_ratio_old = 0 ! don't use predictive timestep control after backup
         s% dt_start = s% dt                   
         s% last_backup = s% model_number               
         s% num_backups = s% num_backups + 1
         if (s% num_backups > s% max_number_backups .and. s% max_number_backups >= 0) then
            write(*,2) 'num_backups', s% num_backups
            write(*,2) 'max_number_backups', s% max_number_backups
            call report_convergence_problems(s, 'num_backups > max_number_backups')
            s% termination_code = t_max_number_backups
            s% result_reason = nonzero_ierr
            do1_backup = terminate
            return
         end if

         s% number_of_backups_in_a_row = s% number_of_backups_in_a_row + 1
         if (s% number_of_backups_in_a_row > s% max_backups_in_a_row &
               .and. s% max_backups_in_a_row > 0) then
            write(*,2) 'number_of_backups_in_a_row', s% number_of_backups_in_a_row
            write(*,2) 'max_backups_in_a_row', s% max_backups_in_a_row
            call report_convergence_problems(s, 'too many backups in a row')
            s% termination_code = t_max_backups_in_a_row
            s% result_reason = nonzero_ierr
            do1_backup = terminate
            return
         end if
         
         if (s% why_Tlim == Tlim_neg_X) then
            s% timestep_hold = s% model_number + &
               max(s% backup_hold, s% neg_mass_fraction_hold)
         else         
            s% timestep_hold = s% model_number + s% backup_hold 
         end if       
         s% retry_cnt = 0
         s% redo_cnt = 0
         s% model_number_for_last_retry = s% model_number
         s% have_start_values = .false.

         if (s% report_ierr) write(*, *) 'backup model_number', &
            s% model_number, s% num_backups
         
      end function do1_backup
      
      
      subroutine set_start_of_step_info(s, ierr)
         use report, only: do_report
         use hydro_vars, only: set_vars
         use mlt_info, only: set_gradT_excess_alpha
         use star_utils, only: dt_Courant
         use alloc, only: non_crit_get_work_array, non_crit_return_work_array

         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         logical :: trace
         integer :: nz, k
         
         include 'formats'
         
         ierr = 0
         trace = s% trace_evolve
         nz = s% nz
         
         do k=1,nz
            s% T_start(k) = -1d99
            s% r_start(k) = -1d99
            s% v_start(k) = -1d99
         end do
                  
         if (trace) write(*,*) 'call set_vars'
         call set_vars(s, s% dt, ierr) ! this does set_mixing_info too
         if (failed('set_vars')) return
         
         if (trace) write(*,*) 'call do_report'
         call do_report(s, ierr) ! set values in case used during step
         if (failed('do_report ierr')) return
         
         ! save a few things from start of step that will need later
         s% dt_Courant_start = dt_Courant(s,k) 
         s% prev_Lmax = maxval(abs(s% L(1:nz)))
         if (s% rotation_flag) then
            s% surf_r_equatorial = s% r_equatorial(1)
         else
            s% surf_r_equatorial = s% r(1)
         end if
         s% starting_T_center = s% T(nz)
         s% surf_opacity = s% opacity(1)
         s% surf_csound = s% csound(1)
         s% surf_rho = s% rho(1)
         s% prev_Ledd = eval_Ledd(s)
         
         if (s% generations == 1) then
            s% surf_accel_grav_ratio = 0
         else
            s% surf_accel_grav_ratio = &
               (s% v_surf - s% v_surf_old)/(s% dt*s% grav(1))
         end if
         
         if (trace) write(*,*) 'call set_gradT_excess_alpha'
         call set_gradT_excess_alpha(s, ierr)
         if (failed('set_gradT_excess_alpha ierr')) return
         
         if (trace) write(*,*) 'call save_prev_mesh_info'
         call save_prev_mesh_info(ierr)
         if (failed('save_prev_mesh_info ierr')) return
         
         
         
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            failed = .true.
            if (s% report_ierr) write(*, *) 'set_start_of_step_info: ' // trim(str)
            s% result_reason = nonzero_ierr
         end function failed
         
         
         subroutine save_prev_mesh_info(ierr)
            integer, intent(out) :: ierr
            real(dp) :: xm
            integer :: k           
            ierr = 0         
            call do1_alloc(s% prev_mesh_xm, nz, ierr)
            if (failed('non_crit_get_work_array')) return         
            call do1_alloc(s% prev_mesh_lnS, nz, ierr)
            if (failed('non_crit_get_work_array')) return         
            call do1_alloc(s% prev_mesh_mu, nz, ierr)
            if (failed('non_crit_get_work_array')) return         
            xm = 0
            do k=1,nz
               s% prev_mesh_xm(k) = xm
               xm = xm + s% dm(k)
               s% prev_mesh_lnS(k) = s% lnS(k)
               s% prev_mesh_mu(k) = s% mu(k)
            end do
            s% prev_mesh_nz = nz
            s% have_prev_lnS = .false.
            s% have_prev_mu = .false.         
         end subroutine save_prev_mesh_info
         
         
         subroutine do1_alloc(p, sz, ierr)
            use alloc, only: non_crit_do1_alloc_if_necessary
            real(dp), pointer :: p(:)
            integer, intent(in) :: sz
            integer, intent(out) :: ierr
            call non_crit_do1_alloc_if_necessary( &
               s, p, sz, 'prepare_for_new_step', ierr)
         end subroutine do1_alloc

      
      end subroutine set_start_of_step_info


      integer function finish_step( &
            id, id_extra, do_photo, &
            how_many_extra_profile_columns, data_for_extra_profile_columns, &
            how_many_extra_history_columns, data_for_extra_history_columns, ierr)
         ! returns keep_going or terminate
         ! if don't return keep_going, then set result_reason to say why.
         use evolve_support, only: output
         use profile, only: do_save_profiles
         use history, only: write_history_info
         use utils_lib, only: free_iounit, number_iounits_allocated
         use alloc, only: size_work_arrays

         integer, intent(in) :: id, id_extra
         logical, intent(in) :: do_photo ! if true, then save "photo" for restart
         interface
            include 'extra_profile_cols.inc'
            include 'extra_history_cols.inc'
         end interface
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer, parameter :: nvals = 1, n_ivals = 0
         integer :: j, k, nz, &
            current_num_iounits_in_use, prev_num_iounits_in_use
         integer :: ivals(n_ivals)
         real(dp) :: vals(nvals), total_energy
         logical :: trace
         
         include 'formats'
         
         finish_step = terminate
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         nz = s% nz
         trace = s% trace_evolve
         
         if (trace) write(*,*) 'start finish_step'     
         
         s% h1_czb_mass_prev = s% h1_czb_mass
         prev_num_iounits_in_use = number_iounits_allocated()
         
         finish_step = keep_going
         s% result_reason = result_reason_normal
         
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         
         if (s% need_to_save_profiles_now .and. s% write_profiles_flag) then
            if (trace) write(*,*) 'call do_save_profiles'
            call do_save_profiles( &
               s, id_extra, how_many_extra_profile_columns, &
               data_for_extra_profile_columns, ierr)
            s% need_to_save_profiles_now = .false.
         end if
         
         call check(1)
         
         if (s% need_to_update_history_now .and. s% do_history_file) then
            if (trace) write(*,*) 'call write_history_info'
            
            call write_history_info( &
               s, id_extra, how_many_extra_history_columns, &
               data_for_extra_history_columns, ierr)
            if (ierr /= 0) then
               finish_step = terminate
               if (s% report_ierr) write(*, *) 'finish_step: write_history_info ierr', ierr
               s% result_reason = nonzero_ierr
               return
            end if
            s% need_to_update_history_now = .false.
         end if
         
         call check(2)
         
         if (do_photo .or. &
               (s% photostep > 0 .and. mod(s% model_number, s% photostep) == 0)) then
               
            if (trace) write(*,*) 'call output'
            call output(id, ierr)

            if (ierr /= 0) then
               finish_step = terminate
               if (s% report_ierr) write(*, *) 'finish_step: output ierr', ierr
               s% result_reason = nonzero_ierr
               return
            end if
         
            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do j=1,s% species
                  write(*,4) 'finish_step after save photo xa(j)', &
                     s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
               end do
            end if

            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do k=1,s% nz
                  write(*,3) 'lnr', s% model_number, k, s% xh(s% i_lnR, k)
                  if (s% rotation_flag) then
                     write(*,3) 'i_rot', s% model_number, k, s% i_rot(k)
                     write(*,3) 'j_rot', s% model_number, k, s% j_rot(k)
                     write(*,3) 'omega', s% model_number, k, s% omega(k)
                  end if
               end do
            end if
            
            
         end if
         
         call check(3)

         s% screening_mode_value = -1 ! force a new lookup for next step         
         s% doing_first_model_of_run = .false.
         s% number_of_backups_in_a_row = 0    
         
         if (trace) write(*,*) 'done finish_step'     
         
                  
         contains
         
         
         subroutine check(i)
            integer, intent(in) :: i
            include 'formats'
            !return
            
            current_num_iounits_in_use = number_iounits_allocated()
            if (current_num_iounits_in_use > 3 .and. &
                  current_num_iounits_in_use > prev_num_iounits_in_use) then
               write(*,2) 's% model_number', s% model_number
               write(*,2) 'prev_num_iounits_in_use', prev_num_iounits_in_use
               write(*,2) 'current_num_iounits_in_use', current_num_iounits_in_use
               write(*,2) 'i', i
               stop 'finish_step' 
            end if
            prev_num_iounits_in_use = current_num_iounits_in_use
         end subroutine check
         
         
      end function finish_step
      
      
      subroutine set_age(id, age, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: age
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% time = age*secyer
         s% star_age = age
         s% profile_age = age
         s% post_he_age = age
      end subroutine set_age



      end module evolve


