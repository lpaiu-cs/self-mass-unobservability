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
 
      module run_star_support

      use star_lib
      use star_def
      use chem_def
      use chem_lib
      use const_def
      use crlibm_lib
      use eos_lib
      use net_def
      use net_lib
      use atm_lib
      use rates_lib, only: set_which_rate_1212
     

      implicit none
      
       
      integer :: id_from_read_star_job = 0
      
      contains 

      subroutine run1_star( &
            do_alloc_star, do_free_star, okay_to_restart, &
            id, restart, &
            extras_controls, &
            extras_startup, &
            extras_check_model, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns, &
            extras_finish_step, &
            extras_after_evolve, &
            ierr, &
            inlist_fname_arg)
         
         logical, intent(in) :: do_alloc_star, do_free_star, okay_to_restart
         integer, intent(inout) :: id ! input if not do_alloc_star
         logical, intent(inout) :: restart ! input if not do_alloc_star
         character (len=*) :: inlist_fname_arg
         integer, intent(out) :: ierr
         optional inlist_fname_arg
         
         interface

            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
     
            integer function extras_startup(id, restart, ierr)
               integer, intent(in) :: id
               logical, intent(in) :: restart
               integer, intent(out) :: ierr
            end function extras_startup
      
            integer function extras_check_model(id, id_extra)
               integer, intent(in) :: id, id_extra
            end function extras_check_model
      
            integer function extras_finish_step(id, id_extra)
               integer, intent(in) :: id, id_extra
               integer :: ierr
            end function extras_finish_step     
      
            subroutine extras_after_evolve(id, id_extra, ierr)
               integer, intent(in) :: id, id_extra
               integer, intent(out) :: ierr
            end subroutine extras_after_evolve

            include 'extra_history_cols.inc'

            include 'extra_profile_cols.inc'

         end interface
         
         integer :: id_extra, result, model_number, source, nsteps, j, ci, nz
         logical :: continue_evolve_loop, first_try, just_did_backup
         type (star_info), pointer :: s
         character (len=strlen) :: inlist_fname
            
         logical, parameter :: dbg = .false.
         
         1 format(a35, 99(1pe26.16))
         2 format(a35, i7, 1pe26.16)
         3 format(a15, 2x, f15.6)
         4 format(a15, 2x, e15.6)

         11 format(a35, f20.10)

         ierr = 0

         call GET_ENVIRONMENT_VARIABLE( &
            "MESA_INLIST", inlist_fname, status=ierr, trim_name=.true.)
         if (ierr /= 0 .or. len_trim(inlist_fname) == 0) then
            if (present(inlist_fname_arg)) then
               inlist_fname = inlist_fname_arg
            else
               inlist_fname = 'inlist'
            end if
         ierr=0
         end if
         
         ! star is initiated here
         call before_evolve_loop(do_alloc_star, okay_to_restart, restart, &
              extras_startup, null_binary_controls, extras_controls, &
              how_many_extra_profile_columns, data_for_extra_profile_columns, &
              id_from_read_star_job, inlist_fname, "restart_photo", &
              dbg, 0, id, id_extra, ierr)
         if (failed('before_evolve_loop',ierr)) return

         call star_ptr(id, s, ierr)
         if (failed('star_ptr',ierr)) return

         continue_evolve_loop = .true.
         s% doing_timing = .false.
         s% job% check_before_step_timing = 0
         s% job% check_step_loop_timing = 0
         s% job% check_after_step_timing = 0
         s% job% time0_initial = 0
         
         if (dbg) write(*,*) 'start evolve_loop'
         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop

            call before_step_loop(s, ierr)
            if (failed('before_step_loop',ierr)) return

            first_try = .true.
            just_did_backup = .false.
         
            step_loop: do ! may need to repeat this loop
            
               if (stop_is_requested(s)) then
                  continue_evolve_loop = .false.
                  result = terminate
                  exit
               end if
            
               result = star_evolve_step(id, first_try, just_did_backup)
               if (result == keep_going) result = star_check_model(id)
               if (result == keep_going) result = extras_check_model(id, id_extra)
               if (result == keep_going) result = star_pick_next_timestep(id)            
               if (result == keep_going) exit step_loop
               
               model_number = get_model_number(id, ierr)
               if (failed('get_model_number',ierr)) return
                              
               if (result == retry .and. s% job% report_retries) then
                  write(*,'(i6,3x,a,/)') model_number, &
                     'retry reason ' // trim(result_reason_str(s% result_reason))
               else if (result == backup .and. s% job% report_backups) then
                  write(*,'(i6,3x,a,/)') model_number, &
                     'backup reason ' // trim(result_reason_str(s% result_reason))
               end if
               
               if (result == redo) then
                  result = star_prepare_to_redo(id)
               end if
               if (result == retry) then
                  result = star_prepare_to_retry(id)
               end if
               if (result == backup) then
                  result = star_do1_backup(id)
                  just_did_backup = .true.
               else
                  just_did_backup = .false.
               end if
               if (result == terminate) then
                  continue_evolve_loop = .false.
                  exit step_loop
               end if
               first_try = .false.
               
            end do step_loop
            
            ! once we get here, the only options are keep_going or terminate.
            ! redo, retry, or backup must be done inside the step_loop
            
            call after_step_loop(s, extras_finish_step, inlist_fname, &
                how_many_extra_history_columns, data_for_extra_history_columns, &
                how_many_extra_profile_columns, data_for_extra_profile_columns, &
                id_extra, dbg, result, ierr)
            if (failed('after_step_loop',ierr)) return
               
            if (result /= keep_going) then
               if (result /= terminate) then
                  write(*,2) 'ERROR in result value in run_star_extras: model', &
                     s% model_number
                  write(*,2) 'result', result
                  exit evolve_loop
               end if
               if (s% result_reason == result_reason_normal) then

                  call terminate_normal_evolve_loop(s, &
                      how_many_extra_history_columns, data_for_extra_history_columns, &
                      how_many_extra_profile_columns, data_for_extra_profile_columns, &
                      id_extra, dbg, result, ierr)
                  if (failed('terminate_normal_evolve_loop',ierr)) return

               end if
               exit evolve_loop
            end if
            
            call do_saves( &
               id, id_extra, s, &
               how_many_extra_history_columns, &
               data_for_extra_history_columns, &
               how_many_extra_profile_columns, &
               data_for_extra_profile_columns)

            if (s% doing_timing) then
               call system_clock(s% job% time1_extra,s% job% clock_rate)
               s% job% after_step_timing = s% job% after_step_timing + &
                  dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
               s% job% check_time_end = eval_total_times(s% id, ierr)
               s% job% check_after_step_timing = s% job% check_after_step_timing + &
                  (s% job% check_time_end - s% job% check_time_start)
            end if
            
         end do evolve_loop

         call after_evolve_loop(s, extras_after_evolve, &
             how_many_extra_history_columns, data_for_extra_history_columns, &
             how_many_extra_profile_columns, data_for_extra_profile_columns, &
             id_extra, do_free_star, ierr)
         if (failed('after_evolve_loop',ierr)) return

      end subroutine run1_star     

      subroutine null_binary_controls(id, binary_id, ierr)
         integer, intent(in) :: id, binary_id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_binary_controls

      ! Binary requires to set some controls in here, which is why
      ! binary_controls and binary_id are arguments. These do nothing
      ! for the case of single star evolution.
      subroutine before_evolve_loop(do_alloc_star, okay_to_restart, restart, &
              extras_startup, binary_controls, extras_controls, &
              how_many_extra_profile_columns, data_for_extra_profile_columns, &
              id_from_read_star_job, inlist_fname, restart_filename, &
              dbg, binary_id, id, id_extra, ierr)
         use se_support, only: se_startup
         logical, intent(in) :: do_alloc_star, okay_to_restart
         logical :: restart
         interface
            integer function extras_startup(id, restart, ierr)
               integer, intent(in) :: id
               logical, intent(in) :: restart
               integer, intent(out) :: ierr
            end function extras_startup
            subroutine binary_controls(id, binary_id, ierr)
               integer, intent(in) :: id, binary_id
               integer, intent(out) :: ierr
            end subroutine binary_controls     
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
            include 'extra_profile_cols.inc'
         end interface
         integer :: id_from_read_star_job
         character (len=*) :: inlist_fname, restart_filename
         logical, intent(in) :: dbg
         integer, intent(in) :: binary_id
         integer, intent(out) :: id, id_extra, ierr

         type (star_info), pointer :: s         

         if (do_alloc_star) then           
            if (id_from_read_star_job /= 0) then 
               ! already allocated by read_star_job
               id = id_from_read_star_job
               id_from_read_star_job = 0
            else
               id = alloc_star(ierr)
               if (failed('alloc_star',ierr)) return
            end if         
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
         else
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return
            call init_starting_star_data(s, ierr)
            if (failed('init_starting_star_data',ierr)) return
         end if
         
         if (dbg) write(*,*) 'call starlib_init'
         call starlib_init(s, ierr) ! okay to do extra calls on this
         if (failed('star_init',ierr)) return
         
         if (dbg) write(*,*) 'call star_set_kap_and_eos_handles'
         call star_set_kap_and_eos_handles(id, ierr)
         if (failed('set_star_kap_and_eos_handles',ierr)) return
         
         if (dbg) write(*,*) 'call star_setup'
         call star_setup(id, inlist_fname, ierr)
         if (failed('star_setup',ierr)) return
         
         if (len_trim(s% op_mono_data_path) == 0) &
            call get_environment_variable( &
               "MESA_OP_MONO_DATA_PATH", s% op_mono_data_path)
         
         if (len_trim(s% op_mono_data_cache_filename) == 0) &
            call get_environment_variable( &
               "MESA_OP_MONO_DATA_CACHE_FILENAME", s% op_mono_data_cache_filename)         
         if (restart_filename /= "restart_photo") &
            restart_filename  = trim(s% photo_directory) // '/' // trim(restart_filename)

         if (okay_to_restart) then
            restart = doing_a_restart(restart_filename)
         else
            restart = .false.
         end if
         
         if (s% job% show_log_description_at_start .and. .not. restart) then
            write(*,*)
            call show_log_description(id, ierr)
            if (failed('show_log_description',ierr)) return
         end if

         call extras_controls(id, ierr)
         if (ierr /= 0) return

         call binary_controls(id, binary_id, ierr)
         if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call do_star_job_controls_before'
         call do_star_job_controls_before(id, s, restart, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call do_load1_star'
         call do_load1_star(id, s, restart, restart_filename, ierr)
         if (failed('do_load1_star',ierr)) return
         
         if (dbg) write(*,*) 'call do_star_job_controls_after'
         call do_star_job_controls_after(id, s, restart, ierr)
         if (failed('do_star_job_controls_after',ierr)) return

         write(*,*)
         write(*,*)
         
         if (.not. restart) then
            call before_evolve(id, ierr)
            if (failed('before_evolve',ierr)) return
            call start_new_run_for_pgstar(s, ierr)
            if (failed('start_new_run_for_pgstar',ierr)) return
         else
            call show_terminal_header(id, ierr)
            if (failed('show_terminal_header',ierr)) return
            call restart_run_for_pgstar(s, ierr)
            if (failed('restart_run_for_pgstar',ierr)) return
         end if
         
         id_extra = extras_startup(id, restart, ierr)
         if (failed('extras_startup',ierr)) return
         
         call se_startup(s, id, restart, s% job% use_se_output, ierr)
         if (failed('se_startup',ierr)) return

         if (s% job% profile_starting_model) then
            write(*, '(a, i12)') 'save profile for model number ', s% model_number
            !call star_set_vars(id, 1, s% nz, ierr)
            if (failed('star_set_vars',ierr)) return
            call save_profile(id, id_extra, &
               how_many_extra_profile_columns, data_for_extra_profile_columns, &
               3, ierr)
            if (failed('save_profile',ierr)) return
         end if

      end subroutine before_evolve_loop

      subroutine before_step_loop(s, ierr)
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr

         integer :: id, k, model_number
         real(dp) :: gamma1_integral, integral_norm
         
         1 format(a35, 99(1pe26.16))
         2 format(a35, i7, 1pe26.16)
         3 format(a15, 2x, f15.6)
         4 format(a15, 2x, e15.6)

         11 format(a35, f20.10)

         s% result_reason = result_reason_normal
         id = s% id
         
         if (s% job% first_model_for_timing >= 0 .and. &
               s% model_number >= s% job% first_model_for_timing .and. &
               .not. s% doing_timing) then
            s% doing_timing = .true.
            write(*,*) 'start timing'
            write(*,*)
            call system_clock(s% job% time0, s% job% clock_rate)
            s% job% time0_initial = s% job% time0
            s% job% step_loop_timing = 0
            s% job% after_step_timing = 0
            s% job% before_step_timing = 0
         end if
         
         if (s% doing_timing) then
            call system_clock(s% job% time0_extra,s% job% clock_rate)
            s% job% check_time_start = eval_total_times(s% id, ierr)
         end if
         
         if (s% use_other_adjust_net) then
            call s% other_adjust_net(s% id, ierr)
            if (failed('other_adjust_net',ierr)) return
         end if
         
         if (s% job% enable_adaptive_network) then
            call star_adjust_net(s% id, &
               s% job% min_x_for_keep, &
               s% job% min_x_for_n, &
               s% job% min_x_for_add, &
               s% job% max_Z_for_add, &
               s% job% max_N_for_add, &
               s% job% max_A_for_add, &
               ierr)
            if (failed('star_adjust_net',ierr)) return
         end if
         
         if (s% job% auto_extend_net) then
            call extend_net(s, ierr)
            if (failed('extend_net',ierr)) return
         end if
         
         if (s% center_ye <= s% job% center_ye_limit_for_v_flag &
               .and. .not. s% v_flag) then
            write(*,1) 'have reached center ye limit', &
               s% center_ye, s% job% center_ye_limit_for_v_flag
            write(*,1) 'set v_flag true'
            call star_set_v_flag(id, .true., ierr)
            if (failed('star_set_v_flag',ierr)) return
            if (ierr /= 0) return
         end if
         
         if (s% log_max_temperature > 9d0 .and. .not. s% v_flag) then 
            ! thanks go to Roni Waldman for this
            gamma1_integral = 0
            integral_norm = 0
            do k=1,s% nz
               integral_norm = integral_norm + s% P(k)*s% dm(k)/s% rho(k)
               gamma1_integral = gamma1_integral + &
                  (s% gamma1(k)-4.d0/3.d0)*s% P(k)*s% dm(k)/s% rho(k)
            end do
            gamma1_integral = gamma1_integral/max(1d-99,integral_norm)
            if (gamma1_integral <= s% job% gamma1_integral_for_v_flag) then
               write(*,1) 'have reached gamma1 integral limit', gamma1_integral
               write(*,1) 'set v_flag true'
               call star_set_v_flag(id, .true., ierr)
               if (failed('star_set_v_flag',ierr)) return
               if (ierr /= 0) return
            end if
         end if
         
         if (s% job% report_mass_not_fe56) call do_report_mass_not_fe56(s)
         if (s% job% report_cell_for_xm > 0) call do_report_cell_for_xm(s)
         
         model_number = get_model_number(id, ierr)
         if (failed('get_model_number',ierr)) return
         
         if (s% doing_timing) then
         
            call system_clock(s% job% time1_extra, s% job% clock_rate)
            s% job% before_step_timing = &
               s% job% before_step_timing + &
                  dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
            
            s% job% check_time_end = eval_total_times(s% id, ierr)
            s% job% check_before_step_timing = &
               s% job% check_before_step_timing + &
                  (s% job% check_time_end - s% job% check_time_start)

            s% job% time0_extra = s% job% time1_extra
            s% job% check_time_start = s% job% check_time_end

         end if

      end subroutine before_step_loop


      subroutine after_step_loop(s, extras_finish_step, inlist_fname, &
             how_many_extra_history_columns, data_for_extra_history_columns, &
             how_many_extra_profile_columns, data_for_extra_profile_columns, &
             id_extra, dbg, result, ierr)
         use se_support, only: se_finish_step
         type (star_info), pointer :: s         
         character (len=*) :: inlist_fname
         integer, intent(in) :: id_extra
         logical, intent(in) :: dbg
         integer, intent(out) :: result, ierr
         interface
            integer function extras_finish_step(id, id_extra)
               integer, intent(in) :: id, id_extra
               integer :: ierr
            end function extras_finish_step     
            include 'extra_history_cols.inc'
            include 'extra_profile_cols.inc'
         end interface

         integer :: id

         id = s% id

         if (s% doing_timing) then            
            call system_clock(s% job% time1_extra,s% job% clock_rate)
            s% job% step_loop_timing = s% job% step_loop_timing + &
               dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate               
            s% job% check_time_end = eval_total_times(s% id, ierr)
            s% job% check_step_loop_timing = s% job% check_step_loop_timing + &
                (s% job% check_time_end - s% job% check_time_start)
            s% job% time0_extra = s% job% time1_extra
            s% job% check_time_start = s% job% check_time_end               
         end if
                     
         if (result == keep_going) then              
            call adjust_tau_factor(s)
            if (s% L_nuc_burn_total/s% L_phot >= s% Lnuc_div_L_zams_limit &
                  .and. .not. s% rotation_flag) then  
               call do_rotation_near_zams(s,ierr)
               if (ierr /= 0) return
            end if           
            if (s% rotation_flag) then      
               call do_rotation(s,ierr)
               if (ierr /= 0) return
            end if 
            ! if you have data that needs to be saved and restored for restarts, 
            ! save it in s% extra_iwork and s% extra_work
            ! before calling star_finish_step
            if (s% job% pgstar_flag) call read_pgstar_inlist(s, inlist_fname, ierr) 
               ! do this before call extras_finish_step
            if (failed('read_pgstar_controls',ierr)) return
            result = extras_finish_step(id, id_extra)               
         end if
                     
         if (result == keep_going) then
            if (dbg) write(*,*) 'call se_finish_step'
            result = se_finish_step(s, id, s% job% use_se_output, &
               how_many_extra_history_columns, data_for_extra_history_columns, &
               how_many_extra_profile_columns, data_for_extra_profile_columns)
         end if
                     
         if (result == keep_going) then
            if (dbg) write(*,*) 'call star_finish_step'
            result = star_finish_step(id, id_extra, .false., &
               how_many_extra_profile_columns, data_for_extra_profile_columns, &
               how_many_extra_history_columns, data_for_extra_history_columns, ierr)
            if (failed('star_finish_step',ierr)) return
         end if
                     
         if (result == keep_going .and. s% job% pgstar_flag) then
            if (dbg) write(*,*) 'call update_pgstar_plots'
            call update_pgstar_plots( &
               s, .false., id_extra, &
               how_many_extra_history_columns, &
               data_for_extra_history_columns, &
               how_many_extra_profile_columns, &
               data_for_extra_profile_columns, &
               ierr)
            if (failed('update_pgstar_plots',ierr)) return
         end if

      end subroutine after_step_loop

      subroutine terminate_normal_evolve_loop(s, &
             how_many_extra_history_columns, data_for_extra_history_columns, &
             how_many_extra_profile_columns, data_for_extra_profile_columns, &
             id_extra, dbg, result, ierr)
         use se_support, only: se_finish_step
         type (star_info), pointer :: s         
         integer, intent(in) :: id_extra
         logical, intent(in) :: dbg
         integer, intent(out) :: result, ierr
         interface
            include 'extra_history_cols.inc'
            include 'extra_profile_cols.inc'
         end interface

         integer :: id

         id = s% id

         if (dbg) write(*,*) 'call star_pick_next_timestep'
         result = star_pick_next_timestep(id) ! for saved model if any  
         if (dbg) write(*,*) 'call se_finish_step'
         result = se_finish_step(s, id, s% job% use_se_output, &
            how_many_extra_history_columns, data_for_extra_history_columns, &
            how_many_extra_profile_columns, data_for_extra_profile_columns)
         if (dbg) write(*,*) 'call save_profile'
         call save_profile(id, id_extra, &
            how_many_extra_profile_columns, data_for_extra_profile_columns, &
            3, ierr)
         s% need_to_save_profiles_now = .false.
         s% need_to_update_history_now = .true.
         if (dbg) write(*,*) 'call star_finish_step'
         result = star_finish_step( &
            id, id_extra, s% job% save_photo_when_terminate, &
            how_many_extra_profile_columns, data_for_extra_profile_columns, &
            how_many_extra_history_columns, data_for_extra_history_columns, ierr)
         if (failed('star_finish_step',ierr)) return
         if (s% job% save_model_when_terminate) &
            s% job% save_model_number = s% model_number                  
         if (s% job% save_pulsation_info_when_terminate) &
            s% job% save_pulsation_info_for_model_number = s% model_number
         if (s% job% write_profile_when_terminate .and. &
               len_trim(s% job% filename_for_profile_when_terminate) > 0) then
            call star_write_profile_info( &
               id, s% job% filename_for_profile_when_terminate, id_extra, &
               how_many_extra_profile_columns, data_for_extra_profile_columns, ierr)
            if (failed('star_write_profile_info',ierr)) return
         end if
         call do_saves( &
            id, id_extra, s, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns)

      end subroutine terminate_normal_evolve_loop

      subroutine after_evolve_loop(s, extras_after_evolve, &
             how_many_extra_history_columns, data_for_extra_history_columns, &
             how_many_extra_profile_columns, data_for_extra_profile_columns, &
             id_extra, do_free_star, ierr)
         use se_support, only: se_after_evolve
         type (star_info), pointer :: s         
         integer, intent(in) :: id_extra
         logical, intent(in) :: do_free_star
         integer, intent(out) :: ierr
         interface
            subroutine extras_after_evolve(id, id_extra, ierr)
               integer, intent(in) :: id, id_extra
               integer, intent(out) :: ierr
            end subroutine extras_after_evolve
            include 'extra_history_cols.inc'
            include 'extra_profile_cols.inc'
         end interface

         integer :: id
         id = s% id

         if (s% doing_timing) then
            call system_clock(s% job% time1,s% job% clock_rate)
            s% job% elapsed_time = &
                dble(s% job% time1 - s% job% time0_initial) / s% job% clock_rate
            call show_times(id,s)
         end if
         
         if (s% result_reason /= result_reason_normal) then
            write(*, *) 
            write(*, '(a)') 'terminated evolution: ' // &
               trim(result_reason_str(s% result_reason))
            write(*, *)
         end if
         
         if (s% termination_code > 0 .and. s% termination_code <= num_termination_codes) then
            write(*, '(a)') 'termination code: ' // &
               trim(termination_code_str(s% termination_code))
         end if
         
         if (s% job% pause_before_terminate) then
            write(*,'(a)') 'pause_before_terminate: hit RETURN to continue'
            read(*,*)
         end if
         
         if (s% result_reason == result_reason_normal) then
         
            if (s% job% pgstar_flag) &
               call update_pgstar_plots( &
                  s, s% job% save_pgstar_files_when_terminate, id_extra, &
                  how_many_extra_history_columns, &
                  data_for_extra_history_columns, &
                  how_many_extra_profile_columns, &
                  data_for_extra_profile_columns, &
                  ierr)
            if (failed('update_pgstar_plots',ierr)) return

            call extras_after_evolve(id, id_extra, ierr)
            if (failed('after_evolve_extras',ierr)) return

            call se_after_evolve(s, id, ierr)
            if (failed('se_after_evolve',ierr)) return

            call star_after_evolve(id, ierr)
            if (failed('star_after_evolve',ierr)) return

            call write_terminal_summary(id, ierr)
            if (failed('write_terminal_summary',ierr)) return
         
         end if
         
         if (do_free_star) then
            call free_star(id, ierr)
            if (failed('free_star',ierr)) return
         end if

      end subroutine after_evolve_loop
         
      
      subroutine adjust_tau_factor(s)
         type (star_info), pointer :: s         
         include 'formats'
         if (s% job% set_tau_factor_after_core_He_burn > 0 .and. &
               abs(s% tau_factor - s% job% set_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% set_to_this_tau_factor)) then
            if (check_for_after_He_burn(s, s% job% set_tau_factor_after_core_He_burn)) then
               s% tau_factor = s% job% set_to_this_tau_factor
               write(*,1) 'set_tau_factor_after_core_He_burn', s% tau_factor
            end if
         end if
   
         if (s% job% set_tau_factor_after_core_C_burn > 0 .and. &
               abs(s% tau_factor - s% job% set_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% set_to_this_tau_factor)) then
            if (check_for_after_C_burn(s, s% job% set_tau_factor_after_core_C_burn)) then
               s% tau_factor = s% job% set_to_this_tau_factor
               write(*,1) 'set_tau_factor_after_core_C_burn', s% tau_factor
            end if
         end if
   
         if (s% job% relax_tau_factor_after_core_He_burn > 0 .and. &
               abs(s% tau_factor - s% job% relax_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% relax_to_this_tau_factor)) then
            if (check_for_after_He_burn(s, s% job% relax_tau_factor_after_core_He_burn)) &
               call relax_tau_factor(s)
         end if
   
         if (s% job% relax_tau_factor_after_core_C_burn > 0 .and. &
               abs(s% tau_factor - s% job% relax_to_this_tau_factor) > &
                  1d-6*max(s% tau_factor, s% job% relax_to_this_tau_factor)) then
            write(*,*) 'call check_for_after_C_burn'
            if (check_for_after_C_burn(s, s% job% relax_tau_factor_after_core_C_burn)) &
               call relax_tau_factor(s)
         end if
      end subroutine adjust_tau_factor

            
      subroutine do_rotation(s,ierr)
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         
         if (s% model_number <= s% job% set_surf_rotation_v_step_limit) then
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/(s% photosphere_r*Rsun)
            write(*,2) 'surface_rotation_v', s% model_number, s% job% new_surface_rotation_v
            write(*,2) 'omega', s% model_number, s% job% new_omega
            call star_set_uniform_omega(s% id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         
         else if (s% model_number <= s% job% set_omega_step_limit) then
            write(*,2) 'omega', s% model_number, s% job% new_omega
            if (failed('star_surface_omega_crit',ierr)) return
            call star_set_uniform_omega(s% id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         
         else if (s% model_number <= s% job% set_omega_div_omega_crit_step_limit) then
            s% job% new_omega = &
               s% job% new_omega_div_omega_crit*star_surface_omega_crit(s% id, ierr)
            write(*,2) 'omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
            write(*,2) 'omega', s% model_number, s% job% new_omega
            if (failed('star_surface_omega_crit',ierr)) return
            call star_set_uniform_omega(s% id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if    
      end subroutine do_rotation  
                     
      
      subroutine do_rotation_near_zams(s,ierr)
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
                      
         if (s% job% set_near_zams_surface_rotation_v_steps > 0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            s% job% set_surf_rotation_v_step_limit = &
               s% model_number + s% job% set_near_zams_surface_rotation_v_steps - 1
            write(*,2) 'near zams: set_surf_rotation_v_step_limit', &
               s% job% set_surf_rotation_v_step_limit

         else if (s% job% set_near_zams_omega_steps > 0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            s% job% set_omega_step_limit = &
               s% model_number + s% job% set_near_zams_omega_steps - 1
            write(*,2) 'near zams: set_omega_step_limit', s% job% set_omega_step_limit

         else if (s% job% set_near_zams_omega_div_omega_crit_steps > 0) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            s% job% set_omega_div_omega_crit_step_limit = &
               s% model_number + s% job% set_near_zams_omega_div_omega_crit_steps - 1
            write(*,2) 'near zams: set_omega_div_omega_crit_step_limit', &
               s% job% set_omega_div_omega_crit_step_limit

         else if (s% job% near_zams_relax_omega) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            write(*,2) 'new_omega', s% model_number, s% job% new_omega
            call star_relax_uniform_omega( &
               s% id, 0, s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return

         else if (s% job% near_zams_relax_omega_div_omega_crit) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
            call star_relax_uniform_omega( &
               s% id, 1, s% job% new_omega_div_omega_crit, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return

         else if (s% job% near_zams_relax_initial_surface_rotation_v) then
            s% job% new_rotation_flag = .true.
            call star_set_rotation_flag(s% id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, s% job% new_surface_rotation_v
            call star_relax_uniform_omega( &
               s% id, 2, s% job% new_surface_rotation_v, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return

         end if  
      end subroutine do_rotation_near_zams   

         
      subroutine relax_tau_factor(s)
         type (star_info), pointer :: s         
         real(dp) :: next
         include 'formats.inc'
         write(*,*) 'relax_to_this_tau_factor < s% tau_factor', &
            s% job% relax_to_this_tau_factor < s% tau_factor
         write(*,1) 'relax_to_this_tau_factor', s% job% relax_to_this_tau_factor
         write(*,1) 's% tau_factor', s% tau_factor
         if (s% job% relax_to_this_tau_factor < s% tau_factor) then
            next = 10**(safe_log10_cr(s% tau_factor) - s% job% dlogtau_factor)
            if (next < s% job% relax_to_this_tau_factor) &
               next = s% job% relax_to_this_tau_factor
         else
            next = 10**(safe_log10_cr(s% tau_factor) + s% job% dlogtau_factor)
            if (next > s% job% relax_to_this_tau_factor) &
               next = s% job% relax_to_this_tau_factor
         end if
         s% tau_factor = next
         write(*,1) 'relax_tau_factor', next, s% job% relax_to_this_tau_factor
      end subroutine relax_tau_factor
      
      
      logical function stop_is_requested(s)
         use utils_lib
         type (star_info), pointer :: s         
         integer :: iounit, ierr
         stop_is_requested = .false.
         if (len_trim(s% job% stop_if_this_file_exists) == 0) return
         ierr = 0
         iounit = alloc_iounit(ierr); if (ierr /= 0) return
         open(unit=iounit, file=trim(s% job% stop_if_this_file_exists), &
            status='old', action='read', iostat=ierr)
         call free_iounit(iounit)
         if (ierr /= 0) return
         close(iounit)
         write(*,*) 'stopping because found file ' // &
            trim(s% job% stop_if_this_file_exists)
         stop_is_requested = .true.
      end function stop_is_requested
      
      
      logical function failed(str,ierr)
         character (len=*), intent(in) :: str
         integer, intent(in) :: ierr
         failed = (ierr /= 0)
         if (failed) write(*, *) trim(str) // ' ierr', ierr
      end function failed
      
      
      subroutine show_times(id, s)
         use utils_lib, only: utils_OMP_GET_MAX_THREADS
         use num_lib, only: qsort
         
         integer, intent(in) :: id
         type (star_info), pointer :: s

         integer, parameter :: max_num_items = 50
         character(len=60) :: item_names(max_num_items)
         real(dp) :: item_values(max_num_items)
         integer, target :: index_arry(max_num_items) 
         integer, pointer :: index(:) 
         integer :: item_order(max_num_items)
         integer :: ierr, omp_num_threads, item_num, num_items, i, j
         real(dp) :: total, misc, tmp
         include 'formats.inc'
         ierr = 0
         omp_num_threads = utils_OMP_GET_MAX_THREADS()
         s% time_total = s% job% check_before_step_timing + &
             s% job% check_step_loop_timing + s% job% check_after_step_timing
         
         write(*,*)
         write(*,'(a50,i8)') 'nz', s% nz
         write(*,'(a50,i8)') 'nvar', s% nvar
         write(*,'(a50,i8)') trim(s% net_name) // ' species', s% species
         write(*,'(a50,i8)') 'total_num_newton_iterations', &
            s% total_num_newton_iterations
         write(*,*)
         write(*,'(a50,i8)') 'threads', omp_num_threads
         total = 0
         item_num = 0
         call save1('remesh', s% time_remesh, total)
         call save1('adjust_mass', s% time_adjust_mass, total)
         call save1('element_diffusion', s% time_element_diffusion, total)
         call save1('burn', s% time_solve_burn, total)
         call save1('mix', s% time_solve_mix, total)
         call save1('semi_explicit_hydro', s% time_semi_explicit_hydro, total)
         call save1('solve', s% time_struct_burn_mix, total)
         call save1('matrix', s% time_newton_matrix, total)
         call save1('omega_mix', s% time_solve_omega_mix, total)
         call save1('eos', s% time_eos, total)
         call save1('neu_and_kap', s% time_neu_kap, total)
         call save1('net', s% time_nonburn_net, total)
         call save1('mlt', s% time_mlt, total)
         call save1('hydro_vars', s% time_set_hydro_vars, total)
         call save1('evolve_step', s% time_evolve_step, total)
         call save1('run1_star', s% job% elapsed_time - total, total)
         call save1('total', total, tmp)
         
         num_items = item_num
         index(1:num_items) => index_arry(1:num_items)
         call qsort(index, num_items, item_values)
         
         write(*,*)
         write(*,*)
         do i=1,num_items
            j = index(num_items+1-i)
            if (item_values(j) == 0d0) cycle
            write(*,'(a50,2f9.3)') trim(item_names(j)), &
               item_values(j), item_values(j)/total
            if (j == num_items) write(*,*)
         end do
      
         
         if (s% job% step_loop_timing/s% job% elapsed_time < 0.9) then
            write(*,*)
            write(*,*)
            write(*,1) 'before_step', s% job% before_step_timing/s% job% elapsed_time
            write(*,1) 'step_loop', s% job% step_loop_timing/s% job% elapsed_time
            write(*,1) 'after_step', s% job% after_step_timing/s% job% elapsed_time
            write(*,*)
         end if
         write(*,*)
         write(*,*)
         
         
         contains
         
         
         subroutine save1(name, value, total)
            character (len=*), intent(in) :: name
            real(dp), intent(in) :: value
            real(dp), intent(inout) :: total
            item_num = item_num + 1
            item_names(item_num) = name
            item_values(item_num) = value
            total = total + value
         end subroutine save1
         

      end subroutine show_times

      
      
      subroutine do_saves( &
            id, id_extra, s, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns)
         integer, intent(in) :: id, id_extra
         type (star_info), pointer :: s
         interface
            include 'extra_profile_cols.inc'
            include 'extra_history_cols.inc'
         end interface

         integer :: ierr
         ierr = 0
      
         if (s% model_number == s% job% save_model_number) then
            call star_write_model(id, s% job% save_model_filename, ierr)
            if (failed('star_write_model',ierr)) return
            write(*, *) 'saved to ' // trim(s% job% save_model_filename)
         end if
         
         if (s% model_number == s% job% save_pulsation_info_for_model_number) then
            call star_write_pulsation_info(id, &
               s% add_center_point_to_pulse_info, &
               s% keep_surface_point_for_pulse_info, &
               s% add_atmosphere_to_pulse_info, &
               s% pulse_info_format, &
               s% job% save_pulsation_info_filename, ierr)
            if (failed('star_write_pulsation_info',ierr)) return
            write(*, *) 'pulsation data saved to ' // &
               trim(s% job% save_pulsation_info_filename)
         end if
         
         if (s% model_number == s% job% save_short_format_for_model_number) then
            call star_write_short_format(id, s% job% save_short_format_filename, ierr)
            if (failed('star_write_short_format',ierr)) return
            write(*, *) 'short_format saved to ' // &
               trim(s% job% save_short_format_filename)
         end if
         
         if (s% model_number == s% job% profile_model_number) then
            write(*, '(a, i7)') 'save profile for model number', s% model_number
            call save_profile(id, id_extra, &
               how_many_extra_profile_columns, data_for_extra_profile_columns, &
               3, ierr)
            if (failed('save_profile',ierr)) return
         end if
         
         if (s% job% internals_num >= 0) then
            write(*, '(a, i7)') 'write internals for model number', s% model_number
            call std_write_internals(id, s% job% internals_num)
            stop 'finished std_write_internals'
         end if
         
      end subroutine do_saves

                  
      subroutine write_colors_info(id, s, ierr)
         use colors_lib
         use colors_def
         use chem_def, only: zsol
         use utils_lib, only: alloc_iounit, free_iounit
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: io, i, j
         character (len=strlen) :: fname
			real(dp)  :: log_Teff ! log10 of surface temp
			real(dp)  :: log_L ! log10 of luminosity in solar units
			real(dp)  :: mass ! mass in solar units
			real(dp)  :: Fe_H ! [Fe/H]
			! output
			real(dp) :: results(n_colors)
			real(dp) :: log_g
         
         ierr = 0
         
         io = alloc_iounit(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in alloc_iounit'
            return
         end if
         
         fname = 'colors.log'
         !if (s% doing_first_model_of_run) then
         if (.false.) then
            open(unit=io, file=trim(fname), action='write', status='replace', iostat=ierr)
            ! write column numbers
            j = 1
            write(io,fmt='(i10)',advance='no') j
            j = j+1
            do i=1,4+n_colors
               write(io,fmt='(i25)',advance='no') j
               j = j+1
            end do
            write(io,fmt='(i25)') j
            ! write column labels
            write(io,fmt='(a10)',advance='no') 'model'
            write(io,fmt='(a25)',advance='no') 'log_Teff'
            write(io,fmt='(a25)',advance='no') 'log_L'
            write(io,fmt='(a25)',advance='no') 'mass'
            write(io,fmt='(a25)',advance='no') 'Fe_H'
            do i=1,n_colors
               write(io,fmt='(a25)',advance='no') trim(colors_name(i))
            end do
            write(io,fmt='(a25)') 'log_g'
         else
            open(unit=io, file=trim(fname), action='write', position='append', iostat=ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed to open colors.log'
            call free_iounit(io)
            return
         end if
         
         log_Teff = log10(s% Teff)
         log_L = s% log_surface_luminosity
         mass = s% star_mass
         Fe_H = safe_log10_cr(get_current_z_at_point(id, 1, ierr) / zsol)
         if (ierr /= 0) then
            write(*,*) 'failed in get_current_z_at_point'
            call cleanup
            return
         end if
         
         call colors_get(log_Teff, log_L, mass, Fe_H, results, log_g, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in colors_get'
            call cleanup
            return
         end if
         
         1 format(1x,f24.12)
         write(io,fmt='(i10)',advance='no') s% model_number
         write(io,fmt=1,advance='no') log_Teff
         write(io,fmt=1,advance='no') log_L
         write(io,fmt=1,advance='no') mass
         write(io,fmt=1,advance='no') Fe_H
         do i=1,n_colors
            write(io,fmt=1,advance='no') results(i)
         end do
         write(io,1) log_g
         
         call cleanup
         
         contains
         
         subroutine cleanup
            close(io)
            call free_iounit(io)
         end subroutine cleanup
      
      end subroutine write_colors_info
      
      
      subroutine read_masses(filename, masses, nmasses, ierr)
         character (len=*), intent(in) :: filename
         real(dp), pointer, intent(out) :: masses(:)
         integer, intent(out) :: nmasses, ierr
         call read_items(filename, masses, nmasses, 'masses', ierr)
      end subroutine read_masses
      
      
      subroutine read_items(filename, items, nitems, name, ierr)
         use utils_lib
         use utils_def
         character (len=*), intent(in) :: filename, name
         real(dp), pointer, intent(out) :: items(:)
         integer, intent(out) :: nitems, ierr
         
         integer :: iounit, n, i, t, capacity
         character (len=strlen) :: buffer, string
         
         nitems = 0
         if (.not. associated(items)) then
            capacity = 10
            allocate(items(capacity))
         else
            capacity = size(items,dim=1)
         end if
         
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return

         open(unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            call free_iounit(iounit)
            write(*,*) 'failed to open file ' // trim(filename)
            return
         end if
         
         n = 0
         i = 0
         
         do
            t = token(iounit, n, i, buffer, string)
            select case(t)
               case(name_token)
                  if (string == name) then
                     call do_read_items(ierr)
                     if (ierr /= 0) then
                        call free_iounit(iounit)
                        return
                     end if
                     exit ! for now, nothing else to be read
                  end if
                  call error; return
               case(eof_token)
                  exit
               case default
                  call error; return
            end select
            
         end do
         
         close(iounit)
         call free_iounit(iounit)
         
         contains
         
         
         subroutine error
            ierr = -1
            write(*,*) 'error in reading file' // trim(filename)
            close(iounit)
            call free_iounit(iounit)
         end subroutine error
         
         
         subroutine do_read_items(ierr)
            integer, intent(out) :: ierr
            real(dp) :: mass
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
         mass_loop: do
               t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  call error; return
               end if
               read(string,fmt=*,iostat=ierr) mass
               if (ierr /= 0) then
                  call error; return
               end if
               nitems = nitems+1
               if (nitems > capacity) then
                  capacity = capacity + 10
                  call realloc_double(items,capacity,ierr)
                  if (ierr /= 0) then
                     call error; return
                  end if
               end if
               items(nitems) = mass
               t = token(iounit, n, i, buffer, string)
               if (t == right_paren_token) exit mass_loop
               if (t /= comma_token) then
                  call error; return
               end if
            end do mass_loop
         end subroutine do_read_items
         
      
      end subroutine read_items
      
      
      subroutine do_report_mass_not_fe56(s)
         use const_def
         type (star_info), pointer :: s
         integer :: k, fe56
         real(dp) :: sumdq
         include 'formats.inc'
         fe56 = s% net_iso(ife56)
         if (fe56 == 0) return
         sumdq = 0
         do k = 1, s% nz
            sumdq = sumdq + s% dq(k)*(1-s% xa(fe56,k))
         end do
         write(*,1) 'R', s% r(1)
         write(*,1) 'g', s% cgrav(1)*s% mstar/(s% r(1)**2)
         write(*,1) 'mass non fe56', s% xmstar*sumdq, sumdq
         write(*,1) 'M_center (Msun)', s% M_center/Msun
         write(*,1) 'xmstar (g)', s% xmstar
         do k=1,s% nz
            if (fe56 == maxloc(s% xa(:,k),dim=1)) then
               write(*,2) 'mass exterior to fe56 (g)', k, (1d0 - s% q(k))*s% xmstar
               write(*,2) 'mass coord top of fe56 (g)', k, s% q(k)*s% xmstar
               return
            end if
         end do
      end subroutine do_report_mass_not_fe56
      
      
      subroutine do_report_cell_for_xm(s)
         use const_def
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: sumdq, dq
         include 'formats.inc'
         dq = s% job% report_cell_for_xm/s% xmstar
         if (dq > 1) then
            write(*,2) 'report_cell_for_xm > xmstar', s% nz
            return
         end if
         sumdq = 0
         do k = 1, s% nz
            sumdq = sumdq + s% dq(k)
            if (sumdq >= dq) then
               write(*,*)
               write(*,2) 'total mass in cells from 1 to k', k, sumdq*s% xmstar
               write(*,2) 'logT(k)', k, s% lnT(k)/ln10
               write(*,2) 'logRho(k)', k, s% lnd(k)/ln10
               write(*,2) 'entropy(k)', k, exp(s% lnS(k))*amu/kerg
               write(*,2) 'xmstar*q(k)', k, s% xmstar*s% q(k)
               write(*,2) 'q(k)', k, s% q(k)
               write(*,*)
               return
            end if
         end do
         write(*,2) 'total mass in cells from 1 to nz', s% nz, s% xmstar
      end subroutine do_report_cell_for_xm
      
      
      subroutine set_which_rates(id, ierr)
         use rates_def
         use rates_lib
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: which_rate
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% job% set_rates_preference) then
            write(*,*) 'change rates preference to', s% job% new_rates_preference
            s% which_rates(:) = s% job% new_rates_preference
         else
            s% which_rates(:) = rates_NACRE_if_available
         end if
         
         if (len_trim(s% job% set_rate_c12ag) > 0) then
            if (s% job% set_rate_c12ag == 'NACRE') then
               which_rate = use_rate_c12ag_NACRE
            else if (s% job% set_rate_c12ag == 'jina reaclib') then
               which_rate = use_rate_c12ag_JR
            else if (s% job% set_rate_c12ag == 'Kunz') then
               which_rate = use_rate_c12ag_Kunz
            else if (s% job% set_rate_c12ag == 'CF88') then
               which_rate = use_rate_c12ag_CF88
            else if (s% job% set_rate_c12ag == 'Buchmann') then
               write(*,*) 'Buchmann rate for c12ag is not in the current jina reaclib'
               write(*,*) 'to use it, switch to the old jina file '
               write(*,*) 'and use set_rate_c12ag == "jina reaclib"'
               write(*,*) '.'
               ierr = -1
               return
            else
               write(*,*) 'invalid string for set_rate_c12ag ' // trim(s% job% set_rate_c12ag)
               write(*,*) 'options are NACRE, jina reaclib, Kunz, CF88'
               ierr = -1
               return
            end if
            call set_which_rate_c12ag(s% which_rates, which_rate)
         end if
         
         if (len_trim(s% job% set_rate_n14pg) > 0) then
            if (s% job% set_rate_n14pg == 'NACRE') then
               which_rate = use_rate_n14pg_NACRE
            else if (s% job% set_rate_n14pg == 'jina reaclib') then
               which_rate = use_rate_n14pg_JR
            else if (s% job% set_rate_n14pg == 'CF88') then
               which_rate = use_rate_n14pg_CF88
            else if (s% job% set_rate_n14pg == 'Imbriani') then
               write(*,*) 'Imbriani rate for n14pg is not in the current jina reaclib'
               write(*,*) 'to use it, switch to the old jina file '
               write(*,*) 'and use set_rate_n14pg == "jina reaclib"'
               write(*,*) '.'
               ierr = -1
               return
            else
               write(*,*) 'invalid string for set_rate_n14pg ' // trim(s% job% set_rate_n14pg)
               write(*,*) 'options are NACRE, jina reaclib, CF88'
               ierr = -1
               return
            end if
            call set_which_rate_n14pg(s% which_rates, which_rate)
         end if
         
         if (len_trim(s% job% set_rate_3a) > 0) then
            if (s% job% set_rate_3a == 'NACRE') then
               which_rate = use_rate_3a_NACRE
            else if (s% job% set_rate_3a == 'jina reaclib') then
               which_rate = use_rate_3a_JR
            else if (s% job% set_rate_3a == 'CF88') then
               which_rate = use_rate_3a_CF88
            else if (s% job% set_rate_3a == 'FL87') then
               which_rate = use_rate_3a_FL87
            else
               write(*,*) 'invalid string for set_rate_3a ' // trim(s% job% set_rate_3a)
               write(*,*) 'options are NACRE, jina reaclib, CF88, FL87'
               ierr = -1
               return
            end if
            call set_which_rate_3a(s% which_rates, which_rate)
         end if
         
         if (len_trim(s% job% set_rate_1212) > 0) then
            if (s% job% set_rate_1212 == 'CF88_basic_1212') then
               which_rate = use_rate_1212_CF88_basic
            else if (s% job% set_rate_1212 == 'CF88_multi_1212') then
               which_rate = use_rate_1212_CF88_multi
            else if (s% job% set_rate_1212 == 'G05') then
               which_rate = use_rate_1212_G05
            else
               write(*,*) 'invalid string for set_rate_1212 ' // trim(s% job% set_rate_1212)
               ierr = -1
               return
            end if
            call set_which_rate_1212(s% which_rates, which_rate)
         end if

      end subroutine set_which_rates
      
      
      subroutine set_rate_factors(id, ierr)
         use net_lib, only: get_net_reaction_table_ptr
         use rates_lib, only: rates_reaction_id
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, i, ir
         integer, pointer :: net_reaction_ptr(:) 
         
         include 'formats.inc'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% rate_factors(:) = 1
         if (s% job% num_special_rate_factors <= 0) return
         
         call get_net_reaction_table_ptr(s% net_handle, net_reaction_ptr, ierr)
         if (ierr /= 0) return
         
         do i=1,s% job% num_special_rate_factors
            if (len_trim(s% job% reaction_for_special_factor(i)) == 0) cycle
            ir = rates_reaction_id(s% job% reaction_for_special_factor(i))
            j = 0
            if (ir > 0) j = net_reaction_ptr(ir)
            if (j <= 0) cycle
            s% rate_factors(j) = s% job% special_rate_factor(i)
            write(*,2) 'set special rate factor for ' // &
                  trim(s% job% reaction_for_special_factor(i)), &
                  j, s% job% special_rate_factor(i)
         end do
         
      end subroutine set_rate_factors

      

      subroutine do_star_job_controls_before(id, s, restart, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr

         integer :: which_atm

         include 'formats.inc'
      
         ierr = 0

         s% set_which_rates => set_which_rates ! will be called after net is defined
         s% set_rate_factors => set_rate_factors ! will be called after net is defined
         
         if (s% job% set_tau_factor .or. &
               (s% job% set_initial_tau_factor .and. .not. restart)) then
            write(*,1) 'set_tau_factor', s% job% set_to_this_tau_factor
            s% tau_factor = s% job% set_to_this_tau_factor
         end if
         
         if (s% job% set_other_HELM_flags) then
            write(*,*) 'set_other_HELM_flags'
            write(*,*) 'HELM_include_radiation', s% job% HELM_include_radiation
            write(*,*) 'HELM_always_skip_elec_pos', s% job% HELM_always_skip_elec_pos
            call eos_set_HELM_flags(s% eos_handle, &
               s% job% HELM_include_radiation, s% job% HELM_always_skip_elec_pos, ierr)
            if (failed('eos_set_HELM_flags',ierr)) return
         end if

         if (s% job% set_Z_all_HELM) then
            write(*,*) 'set_Z_all_HELM'
            write(*,1) 'Z_all_HELM', s% job% Z_all_HELM
            call eos_set_Z_all_HELM(s% eos_handle, s% job% Z_all_HELM, ierr)
            if (failed('eos_set_Z_all_HELM',ierr)) return
         end if

         if (s% job% set_logRho_OPAL_SCVH_limits) then
            write(*,*) 'set_logRho_OPAL_SCVH_limits'
            write(*,1) 'logRho1_OPAL_SCVH_limit', s% job% logRho1_OPAL_SCVH_limit
            write(*,1) 'logRho2_OPAL_SCVH_limit', s% job% logRho2_OPAL_SCVH_limit
            call eos_set_logRhos_OPAL_SCVH( s% eos_handle, &
               s% job% logRho1_OPAL_SCVH_limit, s% job% logRho2_OPAL_SCVH_limit, ierr)
            if (failed('eos_set_logRhos_OPAL_SCVH',ierr)) return
         end if

         if (s% job% set_HELM_SCVH_lgTs) then
            write(*,*) 'set_HELM_SCVH_lgTs'
            write(*,1) 'logT_low_all_HELM', s% job% logT_low_all_HELM
            write(*,1) 'logT_low_all_SCVH', s% job% logT_low_all_SCVH
            if (s% job% logT_low_all_HELM < 2.1d0) then
               write(*,*) 'for current eos tables, min allowed logT_low_all_HELM is 2.1'
               ierr = -1
               return
            end if
            if (s% job% logT_low_all_HELM > s% job% logT_low_all_SCVH) then
               write(*,*) 'logT_low_all_HELM must be <= logT_low_all_SCVH'
               ierr = -1
               return
            end if
            call eos_set_HELM_SCVH_lgTs( &
               s% eos_handle, s% job% logT_low_all_HELM, s% job% logT_low_all_SCVH, ierr)
            if (failed('eos_set_HELM_SCVH_lgTs',ierr)) return
         end if

         if (s% job% set_HELM_OPAL_lgTs) then
            write(*,*) 'set_HELM_OPAL_lgTs'
            write(*,1) 'logT_all_HELM', s% job% logT_all_HELM
            write(*,1) 'logT_all_OPAL', s% job% logT_all_OPAL
            if (s% job% logT_all_HELM > 7.7d0) then
               write(*,*) 'for current eos tables, max allowed logT_all_HELM is 7.7'
               ierr = -1
               return
            end if
            if (s% job% logT_all_HELM < s% job% logT_all_OPAL) then
               write(*,*) 'logT_all_HELM must be >= logT_all_OPAL'
               ierr = -1
               return
            end if
            call eos_set_HELM_OPAL_lgTs( &
               s% eos_handle, s% job% logT_all_HELM, s% job% logT_all_OPAL, ierr)
            if (failed('eos_set_HELM_OPAL_lgTs',ierr)) return
         end if

         if (s% job% set_eos_PC_parameters) then
            write(*,*) 'set_eos_PC_parameters'
            write(*,1) 'mass_fraction_limit_for_PC', s% job% mass_fraction_limit_for_PC
            write(*,1) 'logRho1_PC_limit', s% job% logRho1_PC_limit
            write(*,1) 'logRho2_PC_limit', s% job% logRho2_PC_limit
            write(*,1) 'log_Gamma_all_HELM', s% job% log_Gamma_all_HELM
            write(*,1) 'log_Gamma_all_PC', s% job% log_Gamma_all_PC
            write(*,1) 'PC_min_Z', s% job% PC_min_Z
            call eos_set_PC_parameters(s% eos_handle, &
               s% job% mass_fraction_limit_for_PC, s% job% logRho1_PC_limit, &
               s% job% logRho2_PC_limit, &
               s% job% log_Gamma_all_HELM, s% job% log_Gamma_all_PC, s% job% PC_min_Z, ierr)
            if (failed('star_set_eos_PC_params',ierr)) return
         end if
         
         which_atm = atm_option(s% which_atm_option, ierr)
         if (failed('atm_option',ierr)) return
         s% tau_base = atm_tau_base(which_atm, ierr)
         if (failed('atm_tau_base',ierr)) return
         

      end subroutine do_star_job_controls_before

      
      
      ! in a perfect world, we'd pass s as an arg to this routine.
      ! but for backward compatibility for a large number of users
      ! we do it this strange way instead.
      subroutine do_read_star_job(filename, ierr)
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr  

         integer :: id
         type (star_info), pointer :: s
         
         include 'formats'
         ierr = 0   

         if (id_from_read_star_job /= 0) then
            write(*,2) 'id_from_read_star_job', id_from_read_star_job
            ierr = -1
            return
         end if
         
         id = alloc_star(ierr)
         if (ierr /= 0) then
            write(*,*) 'do_read_star_job failed in alloc_star'
            return
         end if
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_read_star_job failed in star_ptr'
            return
         end if
         
         call read_star_job(s, filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr from read_star_job ' // trim(filename)
            return
         end if
         
         id_from_read_star_job = id

         if (s% job% save_star_job_namelist) then
            call write_star_job(s, s% job% star_job_namelist_name, ierr)
            if (ierr /= 0) then
               write(*,*) 'ierr from write_star_job ' // &
                  trim(s% job% star_job_namelist_name)
               return
            end if
         end if

      end subroutine do_read_star_job
      
      
      subroutine do_load1_star(id, s, restart, restart_filename, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         character (len=*), intent(in) :: restart_filename
         integer, intent(out) :: ierr
      
         if (restart) then
            call star_load_restart_photo(id, restart_filename, ierr)
            if (failed('star_load_restart_photo',ierr)) return
         else if (s% job% load_saved_model) then
            if (s% job% create_pre_main_sequence_model) then
               write(*,*) 'you have both load_saved_model and ' // &
                  'create_pre_main_sequence_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            if (s% job% create_initial_model) then
               write(*,*) 'you have both load_saved_model and create_initial_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            write(*,'(a)') 'load saved model ' // trim(s% job% saved_model_name)
            write(*,*)
            call star_read_model(id, s% job% saved_model_name, ierr)
            if (failed('star_read_model',ierr)) return
         else if (s% job% create_pre_main_sequence_model) then
            if (.not. restart) write(*, *) 'create pre-main-sequence model'
            if (s% job% create_initial_model) then
               write(*,*) 'you have both create_pre_main_sequence_model ' // &
                  'and create_initial_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            call star_create_pre_ms_model( &
               id, s% job% pre_ms_T_c, s% job% pre_ms_guess_rho_c, &
               s% job% pre_ms_d_log10_P, s% job% pre_ms_logT_surf_limit, &
               s% job% pre_ms_logP_surf_limit, s% job% initial_zfracs, &
               s% job% dump_missing_metals_into_heaviest, &
               (s% job% change_net .or. (s% job% change_initial_net .and. .not. restart)), &
               s% job% new_net_name, s% job% pre_ms_relax_num_steps, ierr)
            if (failed('star_create_pre_ms_model',ierr)) return
         else if (s% job% create_initial_model) then
            if (.not. restart) write(*, *) 'create initial model'
            if (s% job% create_pre_main_sequence_model) then
               write(*,*) 'you have both create_initial_model and ' // &
                  'create_pre_main_sequence_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            call star_create_initial_model(id, &
               s% job% radius_in_cm_for_create_initial_model, &
               s% job% mass_in_gm_for_create_initial_model, &
               s% job% center_logP_1st_try_for_create_initial_model, &
               s% job% entropy_1st_try_for_create_initial_model, &
               s% job% max_tries_for_create_initial_model, &
               s% job% abs_e01_tolerance_for_create_initial_model, &
               s% job% abs_e02_tolerance_for_create_initial_model, &
               s% job% initial_zfracs, &
               s% job% dump_missing_metals_into_heaviest, &
               (s% job% change_net .or. (s% job% change_initial_net .and. .not. restart)), &
               s% job% new_net_name, s% job% initial_model_relax_num_steps, &
               s% job% initial_model_eps, &
               ierr)
            if (failed('star_create_initial_model',ierr)) return
         else
            call star_load_zams(id, ierr)
            if (failed('star_load_zams',ierr)) return
         end if

      end subroutine do_load1_star
      

      subroutine extend_net(s, ierr)
         use net_def
         use chem_def
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp), parameter :: tiny = 1d-10, small = 1d-2
         
         real(dp) :: cntr_h, cntr_he
         
         include 'formats.inc'
         
         ierr = 0
         
         if (s% net_name == s% job% adv_net) return

         if (s% net_name == s% job% co_net) then
            if (s% log_max_temperature > 8.8d0 .or. s% log_center_density > 9d0) then
               call change_net(s% job% adv_net)
               if (len_trim(s% job% profile_columns_file) > 0) &
                  write(*,*) 'read ' // trim(s% job% profile_columns_file)
               call star_set_profile_columns( &
                  s% id, s% job% profile_columns_file, .true., ierr)
            end if
            return
         end if
         
         if (s% net_name == s% job% h_he_net) then
            cntr_h = current_abundance_at_point(s% id, ih1, s% nz, ierr)
            if (ierr /= 0) return
            if (cntr_h > tiny) return
            cntr_he = current_abundance_at_point(s% id, ihe4, s% nz, ierr)
            if (ierr /= 0) return
            if (cntr_he > small) return
            if (s% log_max_temperature > 8.3d0 .or. s% log_center_density > 8.5d0) then
               call change_net(s% job% co_net)
               if (len_trim(s% job% profile_columns_file) > 0) &
                  write(*,*) 'read ' // trim(s% job% profile_columns_file)
               call star_set_profile_columns( &
                  s% id, s% job% profile_columns_file, .true., ierr)
            end if
         end if

         
         contains
                  
         
         subroutine change_net(net_name)
            use const_def
            character (len=*), intent(in) :: net_name
            integer :: j
            
            include 'formats'
            
            call star_change_to_new_net( &
               s% id, s% job% adjust_abundances_for_new_isos, net_name, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_change_to_new_net ' // trim(net_name)
               stop 'change_net'
               return
            end if
            
            if (net_name /= s% net_name) then
               write(*,*) '   new net_name ', trim(net_name)
               write(*,*) 'old s% net_name ', trim(s% net_name)
               write(*,*) 'failed to change'
               stop 'change_net'
            end if

            write(*,'(a)') ' new net = ' // trim(s% net_name)
            !do j=1,s% species
            !   write(*,fmt='(a,x)',advance='no') trim(chem_isos% name(s% chem_id(j)))
            !end do
            !write(*,*)
            s% dt_next = s% dt_next/5
            write(*,1) 'reduce timestep', log10(s% dt_next/secyer)
            write(*,*)
         end subroutine change_net
         
         
      end subroutine extend_net         


      subroutine before_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine before_evolve         

       

      subroutine do_star_job_controls_after(id, s, restart, ierr)
         use const_def
         use omp_lib
         use rates_def
         use rates_lib

         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         
         real(dp) :: log_m, log_lifetime, max_dt, minq, maxq
         integer :: i, j, k, nzlo, nzhi, chem_id, chem_id1, chem_id2

         include 'formats.inc'
         
         if (s% job% set_initial_age .and. .not. restart) then
            write(*,1) 'set_initial_age', s% job% initial_age ! in years
            call star_set_age(id, s% job% initial_age, ierr)
            if (failed('star_set_age',ierr)) return
         end if

         if (s% job% set_initial_dt .and. .not. restart) then
            write(*,1) 'set_initial_dt', s% job% years_for_initial_dt
            s% dt_next = s% job% years_for_initial_dt*secyer
         end if

         if (s% job% limit_initial_dt .and. .not. restart) then
            write(*,1) 'limit_initial_dt', s% job% years_for_initial_dt
            s% dt_next = min(s% dt_next, s% job% years_for_initial_dt*secyer)
         end if

         if (s% job% set_initial_model_number .and. .not. restart) then
            write(*,2) 'set_initial_model_number', s% job% initial_model_number
            s% model_number = s% job% initial_model_number
         end if

         if (s% job% steps_to_take_before_terminate > 0) then
            s% max_model_number = s% model_number + s% job% steps_to_take_before_terminate
            write(*,2) 'steps_to_take_before_terminate', &
               s% job% steps_to_take_before_terminate
            write(*,2) 'max_model_number', s% max_model_number
         end if

         if (s% job% change_net .or. (s% job% change_initial_net .and. .not. restart)) then         
            write(*,*) 'change to ' // trim(s% job% new_net_name)
            call star_change_to_new_net( &
               id, s% job% adjust_abundances_for_new_isos, s% job% new_net_name, ierr)
            if (failed('star_change_to_new_net',ierr)) return
            write(*,*) 'number of species =', s% species
         end if
         
         if (abs(s% job% T9_weaklib_full_off - T9_weaklib_full_off) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_off', s% job% T9_weaklib_full_off
            T9_weaklib_full_off = s% job% T9_weaklib_full_off
         end if
         
         if (abs(s% job% T9_weaklib_full_on - T9_weaklib_full_on) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_on', s% job% T9_weaklib_full_on
            T9_weaklib_full_on = s% job% T9_weaklib_full_on
         end if
         
         if (s% job% weaklib_blend_hi_Z /= weaklib_blend_hi_Z) then
            write(*,1) 'set weaklib_blend_hi_Z', s% job% weaklib_blend_hi_Z
            weaklib_blend_hi_Z = s% job% weaklib_blend_hi_Z
         end if
         
         if (abs(s% job% T9_weaklib_full_off_hi_Z - T9_weaklib_full_off_hi_Z) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_off_hi_Z', s% job% T9_weaklib_full_off_hi_Z
            T9_weaklib_full_off_hi_Z = s% job% T9_weaklib_full_off_hi_Z
         end if
         
         if (abs(s% job% T9_weaklib_full_on_hi_Z - T9_weaklib_full_on_hi_Z) > 1d-6) then
            write(*,1) 'set T9_weaklib_full_on_hi_Z', s% job% T9_weaklib_full_on_hi_Z
            T9_weaklib_full_on_hi_Z = s% job% T9_weaklib_full_on_hi_Z
         end if

         ! set up coulomb corrections for the special weak rates
         which_mui_coulomb = get_mui_value(s% job% ion_coulomb_corrections)
         which_vs_coulomb = get_vs_value(s% job% electron_coulomb_corrections)

         if (s% job% change_lnPgas_flag .or. &
               (s% job% change_initial_lnPgas_flag .and. .not. restart)) then
            write(*,*) 'new_lnPgas_flag =', s% job% new_lnPgas_flag
            call star_set_lnPgas_flag(id, s% job% new_lnPgas_flag, ierr)
            if (failed('star_set_lnPgas_flag',ierr)) return
         end if

         if (s% job% change_v_flag .or. &
               (s% job% change_initial_v_flag .and. .not. restart)) then
            write(*,*) 'new_v_flag =', s% job% new_v_flag
            call star_set_v_flag(id, s% job% new_v_flag, ierr)
            if (failed('star_set_v_flag',ierr)) return
         end if

         if (s% job% change_L_flag .or. &
               (s% job% change_initial_L_flag .and. .not. restart)) then
            write(*,*) 'new_L_flag =', s% job% new_L_flag
            call star_set_L_flag(id, s% job% new_L_flag, ierr)
            if (failed('star_set_L_flag',ierr)) return
         end if

         if (s% job% change_E_flag .or. &
               (s% job% change_initial_E_flag .and. .not. restart)) then
            write(*,*) 'new_E_flag =', s% job% new_E_flag
            call star_set_E_flag(id, s% job% new_E_flag, ierr)
            if (failed('star_set_E_flag',ierr)) return
         end if

         if (s% job% change_rotation_flag .or. &
               (s% job% change_initial_rotation_flag .and. .not. restart)) then
            write(*,*) 'new_rotation_flag =', s% job% new_rotation_flag
            call star_set_rotation_flag(id, s% job% new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_omega) then
            write(*,1) 'new_omega =', s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_initial_omega .and. .not. restart) then
            write(*,1) 'new_omega =', s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_surface_rotation_v) then
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/s% r(1)
            write(*,1) 'new_surface_rotation_v =', &
               s% job% new_surface_rotation_v, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. &
             s% job% set_initial_surface_rotation_v .and. .not. restart) then
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/s% r(1)
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, s% job% new_surface_rotation_v, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% set_omega_div_omega_crit) then
            s% job% new_omega = &
               s% job% new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. &
             s% job% set_initial_omega_div_omega_crit .and. .not. restart) then
            s% job% new_omega = &
               s% job% new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit, s% job% new_omega
            call star_set_uniform_omega(id, s% job% new_omega, ierr)
            if (failed('star_set_uniform_omega',ierr)) return
         end if
         
         if (s% job% set_to_xa_for_accretion .or. &
               (s% job% set_initial_to_xa_for_accretion .and. .not. restart)) then
            write(*,*) 'set_to_xa_for_accretion'
            call change_to_xa_for_accretion(id, s% job% set_nzlo, s% job% set_nzhi, ierr)
            if (failed('set_to_xa_for_accretion',ierr)) return
         end if
         
         if (s% job% first_model_for_timing > 0) &
            write(*,2) 'first_model_for_timing', s% job% first_model_for_timing
         
         if (s% job% set_uniform_initial_composition .and. .not. restart) then
            write(*,*)
            write(*,1) 'set_uniform_initial_composition'
            write(*,1) 'initial_h1', s% job% initial_h1
            write(*,1) 'initial_h2', s% job% initial_h2
            write(*,1) 'initial_he3', s% job% initial_he3
            write(*,1) 'initial_he4', s% job% initial_he4
            select case(s% job% initial_zfracs)
               case (AG89_zfracs)
                  write(*,1) 'metals AG89'
               case (GN93_zfracs)
                  write(*,1) 'metals GN93'
               case (GS98_zfracs)
                  write(*,1) 'metals GS98'
               case (L03_zfracs)
                  write(*,1) 'metals L03'
               case (AGS05_zfracs)
                  write(*,1) 'metals AGS05'
               case (AGSS09_zfracs)
                  write(*,1) 'metals AGSS09'
               case (L09_zfracs)
                  write(*,1) 'metals L09'
               case (A09_Prz_zfracs)
                  write(*,1) 'metals A09_Prz'
               case default
                  write(*,2) 'unknown value for initial_zfracs', s% job% initial_zfracs
            end select
            call star_set_standard_composition( &
               id, s% job% initial_h1, s% job% initial_h2, &
               s% job% initial_he3, s% job% initial_he4, s% job% initial_zfracs, &
               s% job% dump_missing_metals_into_heaviest, ierr)
            if (failed('set_uniform_initial_composition',ierr)) return
         end if
         
         if (s% job% relax_initial_composition .and. .not. restart) then
            call do_relax_initial_composition(ierr)
            if (failed('do_relax_initial_composition',ierr)) return
         end if
         
         if (s% job% relax_initial_to_xaccrete .and. .not. restart) then
            call star_relax_to_xaccrete(id, s% job% num_steps_to_relax_composition, ierr)
            if (failed('star_relax_to_xaccrete',ierr)) return
         end if

         if (s% job% set_uniform_xa_from_file) then
            call star_uniform_xa_from_file(id, s% job% file_for_uniform_xa, ierr)
            if (failed('star_uniform_xa_from_file',ierr)) return
         end if

         if (s% job% mix_initial_envelope_down_to_T > 0d0 .and. .not. restart) then
            call uniform_mix_envelope_down_to_T(id, s% job% mix_initial_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T',ierr)) return
         end if

         if (s% job% mix_envelope_down_to_T > 0d0) then
            call uniform_mix_envelope_down_to_T(id, s% job% mix_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T',ierr)) return
         end if

         if (s% job% mix_initial_envelope_down_to_T > 0d0) then
            call uniform_mix_envelope_down_to_T(id, s% job% mix_initial_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T',ierr)) return
         end if

         if (s% job% set_uniform_initial_xa_from_file .and. .not. restart) then
            call star_uniform_xa_from_file(id, s% job% file_for_uniform_xa, ierr)
            if (failed('star_uniform_xa_from_file',ierr)) return
         end if
         
         ! do change Z before change Y since changing Z can change Y
         if (s% job% change_Z) then
            write(*,1) 'new_Z', s% job% new_Z
            call star_set_z(id, s% job% new_Z, ierr)
            if (failed('star_set_z',ierr)) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z',ierr)) return
         end if

         if (s% job% change_initial_Z .and. .not. restart) then
            write(*,1) 'new_Z', s% job% new_Z
            call star_set_z(id, s% job% new_Z, ierr)
            if (failed('star_set_z',ierr)) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z',ierr)) return
         end if

         if (s% job% change_Y) then
            write(*,1) 'new_Y', s% job% new_Y
            call star_set_y(id, s% job% new_Y, ierr)
            if (failed('change_Y',ierr)) return
            write(*, 1) 'new y', get_current_abundance(id, ihe4, ierr)
            if (failed('get_current_abundance',ierr)) return
         end if

         if (s% job% change_initial_Y .and. .not. restart) then
            write(*,1) 'new_Y', s% job% new_Y
            call star_set_y(id, s% job% new_Y, ierr)
            if (failed('change_initial_Y',ierr)) return
            write(*, 1) 'new y', get_current_abundance(id, ihe4, ierr)
            if (failed('get_current_abundance',ierr)) return
         end if

         if (s% job% set_abundance .or. &
               (s% job% set_initial_abundance .and. .not. restart)) then
            nzlo = s% job% set_abundance_nzlo
            nzhi = s% job% set_abundance_nzhi
            if (nzhi <= 0) nzhi = s% nz
            if (nzlo <= 0) nzlo = 1
            write(*, *) 'set_abundance of ', &
               trim(s% job% chem_name), s% job% new_frac, nzlo, nzhi 
            chem_id = get_nuclide_index(s% job% chem_name)
            if (chem_id <= 0) then
               write(*,*) 'failed to find ' // trim(s% job% chem_name)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            call set_abundance_in_section(id, chem_id, s% job% new_frac, nzlo, nzhi, ierr)
            if (failed('set_abundance_in_section',ierr)) return
         end if
         
         if (s% job% replace_element .or. &
               (s% job% replace_initial_element .and. .not. restart)) then
            write(*, *) 'replace_element ', &
               trim(s% job% chem_name1), ' by ', trim(s% job% chem_name2)
            chem_id1 = get_nuclide_index(s% job% chem_name1)
            chem_id2 = get_nuclide_index(s% job% chem_name2)
            if (chem_id1 <= 0) then
               write(*,*) 'failed to find ' // trim(s% job% chem_name1)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            if (chem_id2 <= 0) then
               write(*,*) 'failed to find ' // trim(s% job% chem_name2)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            nzhi = s% job% replace_element_nzlo
            nzlo = s% job% replace_element_nzhi
            if (nzhi <= 0) nzhi = s% nz
            if (nzlo <= 0) nzlo = 1
            write(*, *) 'in section', nzlo, nzhi
            call replace_element_in_section( &
               id, chem_id1, chem_id2, nzlo, nzhi, ierr)
            if (failed('replace_element_in_section',ierr)) return
         end if

         if (s% job% set_irradiation .or. &
               (s% job% set_initial_irradiation .and. .not. restart)) then
            write(*,2) 'set_irradiation'
            s% irradiation_flux = s% job% set_to_this_irrad_flux
            s% column_depth_for_irradiation = s% job% irrad_col_depth
         end if
         
         if (s% job% do_special_test) then
            write(*, *) 'do_special_test'
            call star_special_test(id, ierr)
            if (failed('star_special_test',ierr)) return
         end if
         
         if (s% job% remove_center_by_temperature > 0) then
            write(*, 1) 'remove_center_by_temperature', s% job% remove_center_by_temperature
            call star_remove_center_by_temperature( &
               id, s% job% remove_center_by_temperature, ierr)
            if (failed('star_remove_center_by_temperature',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_temperature > 0 .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_temperature', &
               s% job% remove_initial_center_by_temperature
            call star_remove_center_by_temperature( &
               id, s% job% remove_initial_center_by_temperature, ierr)
            if (failed('star_remove_center_by_temperature',ierr)) return
         end if
         
         if (s% job% remove_center_by_radius_cm > s% R_center .and. &
               s% job% remove_center_by_radius_cm < s% r(1)) then
            write(*, 1) 'remove_center_by_radius_cm', &
               s% job% remove_center_by_radius_cm
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_center_by_radius_cm, ierr)
            if (failed('star_remove_center_by_radius_cm',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_radius_cm > s% R_center .and. &
               s% job% remove_initial_center_by_radius_cm < s% r(1) .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_radius_cm', &
               s% job% remove_initial_center_by_radius_cm
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_initial_center_by_radius_cm, ierr)
            if (failed('star_remove_center_by_radius_cm',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_mass_fraction_q > 0d0 .and. &
               s% job% remove_initial_center_by_mass_fraction_q < 1d0 &
                  .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_mass_fraction_q', &
               s% job% remove_initial_center_by_mass_fraction_q
            call star_remove_center_by_mass_fraction_q( &
               id, s% job% remove_initial_center_by_mass_fraction_q, ierr)
            if (failed('star_remove_initial_center_by_mass_fraction_q',ierr)) return
         end if
         
         if (s% job% remove_center_by_mass_fraction_q > 0d0 .and. &
               s% job% remove_center_by_mass_fraction_q < 1d0) then
            write(*, 1) 'remove_center_by_mass_fraction_q', &
               s% job% remove_center_by_mass_fraction_q
            call star_remove_center_by_mass_fraction_q( &
               id, s% job% remove_center_by_mass_fraction_q, ierr)
            if (failed('star_remove_center_by_mass_fraction_q',ierr)) return
         end if
         
         if (s% job% remove_center_by_mass_gm > s% M_center .and. &
               s% job% remove_center_by_mass_gm < s% m(1)) then
            write(*, 1) 'remove_center_by_mass_gm', &
               s% job% remove_center_by_mass_gm
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_center_by_mass_gm, ierr)
            if (failed('star_remove_center_by_mass_gm',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_mass_gm > s% M_center .and. &
               s% job% remove_initial_center_by_mass_gm < s% m(1) .and. .not. restart) then
            write(*, 1) 'remove_initial_center_by_mass_gm', &
               s% job% remove_initial_center_by_mass_gm
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_initial_center_by_mass_gm, ierr)
            if (failed('star_remove_center_by_mass_gm',ierr)) return
         end if
         
         if (s% job% remove_center_by_radius_Rsun > s% R_center/Rsun .and. &
               s% job% remove_center_by_radius_Rsun < s% r(1)/Rsun) then
            write(*, 1) 'remove_center_by_radius_Rsun', &
               s% job% remove_center_by_radius_Rsun
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_center_by_radius_Rsun*Rsun, ierr)
            if (failed('star_remove_center_by_radius_Rsun',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_radius_Rsun > s% R_center/Rsun .and. &
               s% job% remove_initial_center_by_radius_Rsun < s% r(1)/Rsun .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_center_by_radius_Rsun', &
               s% job% remove_initial_center_by_radius_Rsun
            call star_remove_center_by_radius_cm( &
               id, s% job% remove_initial_center_by_radius_Rsun*Rsun, ierr)
            if (failed('star_remove_center_by_radius_Rsun',ierr)) return
         end if
         
         if (s% job% remove_center_by_mass_Msun > s% M_center/Msun .and. &
               s% job% remove_center_by_mass_Msun < s% m(1)/Msun) then
            write(*, 1) 'remove_center_by_mass_Msun', &
               s% job% remove_center_by_mass_Msun
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_center_by_mass_Msun*Msun, ierr)
            if (failed('star_remove_center_by_mass_Msun',ierr)) return
         end if
         
         if (s% job% remove_initial_center_by_mass_Msun > s% M_center/Msun .and. &
               s% job% remove_initial_center_by_mass_Msun < s% m(1)/Msun .and. &
               .not. restart) then
            write(*, 1) 'remove_initial_center_by_mass_Msun', &
               s% job% remove_initial_center_by_mass_Msun
            call star_remove_center_by_mass_gm( &
               id, s% job% remove_initial_center_by_mass_Msun*Msun, ierr)
            if (failed('star_remove_center_by_mass_Msun',ierr)) return
         end if

         if (s% job% remove_initial_center_at_cell_k > 0 .and. .not. restart .and. &
               s% job% remove_initial_center_at_cell_k <= s% nz) then
            write(*, 2) 'remove_initial_center_at_cell_k', s% job% remove_initial_center_at_cell_k
            call star_remove_center_at_cell_k( &
               id, s% job% remove_initial_center_at_cell_k, ierr)
            if (failed('star_remove_center_at_cell_k',ierr)) return
         end if

         if (s% job% remove_center_at_cell_k > 0 .and. &
               s% job% remove_center_at_cell_k <= s% nz) then
            write(*, 2) 'remove_center_at_cell_k', s% job% remove_center_at_cell_k
            call star_remove_center_at_cell_k(id, s% job% remove_center_at_cell_k, ierr)
            if (failed('star_remove_center_at_cell_k',ierr)) return
         end if
         
         if (s% job% set_v_center .or. &
               (s% job% set_initial_v_center .and. .not. restart)) then
            write(*, 1) 'set_v_center', s% job% new_v_center
            s% v_center = s% job% new_v_center
         end if
         
         if (s% job% set_L_center .or. &
               (s% job% set_initial_L_center .and. .not. restart)) then
            write(*, 1) 'set_L_center', s% job% new_L_center
            s% L_center = s% job% new_L_center
         end if

         ! do "set" before "relax"
         
         ! must do relax Z before relax Y since relax Z can change Y
         ! (Warrick Ball pointed out this requirement)
         if (s% job% relax_initial_Z .and. .not. restart) then
            minq = 0
            maxq = 1
            write(*,1) 'relax_initial_Z', s% job% new_Z
            call star_relax_Z(id, s% job% new_Z, s% relax_dlnZ, minq, maxq, ierr)
            if (failed('star_relax_Z',ierr)) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z',ierr)) return
         end if

         if (s% job% relax_Z) then
            minq = 0
            maxq = 1
            write(*,1) 'relax_Z', s% job% new_Z
            call star_relax_Z(id, s% job% new_Z, s% relax_dlnZ, minq, maxq, ierr)
            if (failed('star_relax_Z',ierr)) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z',ierr)) return
         end if

         if (s% job% relax_initial_Y .and. .not. restart) then
            write(*,1) 'relax_initial_Y', s% job% new_Y
            call star_relax_Y(id, s% job% new_Y, s% relax_dY, ierr)
            if (failed('star_relax_Y',ierr)) return
            write(*, 1) 'new y', get_current_y(id, ierr)
            if (failed('get_current_y',ierr)) return
         end if

         if (s% job% relax_Y) then
            write(*,1) 'relax_Y', s% job% new_Y
            call star_relax_Y(id, s% job% new_Y, s% relax_dY, ierr)
            if (failed('star_relax_Y',ierr)) return
            write(*, 1) 'new y', get_current_y(id, ierr)
            if (failed('get_current_y',ierr)) return
         end if

         if (s% job% relax_mass) then
            write(*, 1) 'relax_mass', s% job% new_mass
            call star_relax_mass(id, s% job% new_mass, s% job% lg_max_abs_mdot, ierr)
            if (failed('star_relax_mass',ierr)) return
         end if

         if (s% job% relax_dxdt_nuc_factor .or. &
               (s% job% relax_initial_dxdt_nuc_factor .and. .not. restart)) then
            write(*, 1) 'relax_dxdt_nuc_factor', s% job% new_dxdt_nuc_factor
            call star_relax_dxdt_nuc_factor( &
               id, s% job% new_dxdt_nuc_factor, s% job% dxdt_nuc_factor_multiplier, ierr)
            if (failed('star_relax_dxdt_nuc_factor',ierr)) return
         end if

         if (s% job% relax_eps_nuc_factor .or. &
               (s% job% relax_initial_eps_nuc_factor .and. .not. restart)) then
            write(*, 1) 'relax_eps_nuc_factor', s% job% new_eps_nuc_factor
            call star_relax_eps_nuc_factor( &
               id, s% job% new_eps_nuc_factor, s% job% eps_nuc_factor_multiplier, ierr)
            if (failed('star_relax_eps_nuc_factor',ierr)) return
         end if

         if (s% job% relax_opacity_max .or. &
               (s% job% relax_initial_opacity_max .and. .not. restart)) then
            write(*, 1) 'relax_opacity_max', s% job% new_opacity_max
            call star_relax_opacity_max( &
               id, s% job% new_opacity_max, s% job% opacity_max_multiplier, ierr)
            if (failed('star_relax_opacity_max',ierr)) return
         end if

         if (s% job% relax_max_surf_dq .or. &
               (s% job% relax_initial_max_surf_dq .and. .not. restart)) then
            write(*, 1) 'relax_max_surf_dq', s% job% new_max_surf_dq
            call star_relax_max_surf_dq( &
               id, s% job% new_max_surf_dq, s% job% max_surf_dq_multiplier, ierr)
            if (failed('star_relax_max_surf_dq',ierr)) return
         end if

         if (s% job% relax_initial_mass .and. .not. restart) then
            write(*, 1) 'relax_initial_mass to new_mass =', s% job% new_mass
            call star_relax_mass(id, s% job% new_mass, s% job% lg_max_abs_mdot, ierr)
            if (failed('relax_initial_mass',ierr)) return
         end if

         if (s% job% relax_mass_scale .or. &
               (s% job% relax_initial_mass_scale .and. .not. restart)) then
            write(*, 1) 'relax_mass_scale', s% job% new_mass
            call star_relax_mass_scale( &
               id, s% job% new_mass, s% job% dlgm_per_step, &
               s% job% change_mass_years_for_dt, ierr)
            if (failed('star_relax_mass_scale',ierr)) return
         end if

         if (s% job% relax_core .or. &
               (s% job% relax_initial_core .and. .not. restart)) then
            write(*, 1) 'relax_core', s% job% new_core_mass
            call star_relax_core( &
               id, s% job% new_core_mass, s% job% dlg_core_mass_per_step, &
               s% job% relax_core_years_for_dt, &
               s% job% core_avg_rho, s% job% core_avg_eps, ierr)
            if (failed('star_relax_core',ierr)) return
         end if
         
         if (s% job% relax_M_center .or. &
               (s% job% relax_initial_M_center .and. .not. restart)) then
            write(*, 1) 'relax_M_center', s% job% new_mass
            call star_relax_M_center( &
               id, s% job% new_mass, s% job% dlgm_per_step, s% job% relax_M_center_dt, ierr)
            if (failed('star_relax_M_center',ierr)) return
         end if
         
         if (s% job% relax_R_center .or. &
               (s% job% relax_initial_R_center .and. .not. restart)) then
            write(*, 1) 'relax_R_center', s% job% new_R_center
            call star_relax_R_center( &
               id, s% job% new_R_center, s% job% dlgR_per_step, s% job% relax_R_center_dt, ierr)
            if (failed('star_relax_R_center',ierr)) return
         end if
         
         if (s% job% relax_v_center .or. &
               (s% job% relax_initial_v_center .and. .not. restart)) then
            write(*, 1) 'relax_v_center', s% job% new_v_center
            call star_relax_v_center( &
               id, s% job% new_v_center, s% job% dv_per_step, s% job% relax_v_center_dt, ierr)
            if (failed('star_relax_v_center',ierr)) return
         end if
         
         if (s% job% relax_L_center .or. &
               (s% job% relax_initial_L_center .and. .not. restart)) then
            write(*, 1) 'relax_L_center', s% job% new_L_center
            call star_relax_L_center( &
               id, s% job% new_L_center, s% job% dlgL_per_step, s% job% relax_L_center_dt, ierr)
            if (failed('star_relax_L_center',ierr)) return
         end if
         
         if (s% job% relax_tau_factor .or. &
               (s% job% relax_initial_tau_factor .and. .not. restart)) then
            write(*,1) 'relax_tau_factor', s% job% relax_to_this_tau_factor
            call star_relax_tau_factor( &
               id, s% job% relax_to_this_tau_factor, s% job% dlogtau_factor, ierr)
            if (failed('star_relax_tau_factor',ierr)) return
         end if

         if (s% job% relax_irradiation .or. &
               (s% job% relax_initial_irradiation .and. .not. restart)) then
            write(*,2) 'relax_irradiation -- min steps', s% job% relax_irradiation_min_steps
            write(*,1) 'relax_irradiation -- max yrs dt', s% job% relax_irradiation_max_yrs_dt
            call star_relax_irradiation(id, &
               s% job% relax_irradiation_min_steps, &
               s% job% relax_to_this_irrad_flux, s% job% irrad_col_depth, &
               s% job% relax_irradiation_max_yrs_dt, ierr)
            if (failed('star_relax_irradiation',ierr)) return
         end if

         if (s% job% relax_mass_change .or. &
               (s% job% relax_initial_mass_change .and. .not. restart)) then
            write(*,2) 'relax_mass_change -- min steps', &
               s% job% relax_mass_change_min_steps
            write(*,1) 'relax_mass_change -- max yrs dt', &
               s% job% relax_mass_change_max_yrs_dt
            write(*,1) 'relax_mass_change -- initial_mass_change', &
               s% job% relax_mass_change_init_mdot
            write(*,1) 'relax_mass_change -- final_mass_change', &
               s% job% relax_mass_change_final_mdot
            call star_relax_mass_change(id, &
               s% job% relax_mass_change_min_steps, &
               s% job% relax_mass_change_init_mdot, &
               s% job% relax_mass_change_final_mdot, &
               s% job% relax_mass_change_max_yrs_dt, ierr)
            if (failed('star_relax_mass_change',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% relax_omega) then
            write(*,1) 'new_omega =', s% job% new_omega
            call star_relax_uniform_omega( &
               id, 0, s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
         end if

         if (s% rotation_flag .and. s% job% relax_initial_omega .and. .not. restart) then
            call star_relax_uniform_omega( &
               id, 0, s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,1) 'new_omega =', s% job% new_omega
         end if

         if (s% rotation_flag .and. s% job% relax_surface_rotation_v) then
            call star_relax_uniform_omega( &
               id, 1, s% job% new_surface_rotation_v, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            s% job% new_omega = s% job% new_surface_rotation_v*1d5/s% r(1)
            write(*,1) 'new_surface_rotation_v =', &
               s% job% new_surface_rotation_v, s% job% new_omega
         end if

         if (s% rotation_flag .and. &
               s% job% relax_initial_surface_rotation_v .and. .not. restart) then
            write(*,1) 'new_omega', s% job% new_omega
            write(*,*) 'call star_relax_uniform_omega'
            call star_relax_uniform_omega( &
               id, 1, s% job% new_surface_rotation_v, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, s% job% new_surface_rotation_v
         end if

         if (s% rotation_flag .and. s% job% relax_omega_div_omega_crit) then
            if (failed('star_surface_omega_crit',ierr)) return
            call star_relax_uniform_omega( &
               id, 2, s% job% new_omega_div_omega_crit, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
         end if

         if (s% rotation_flag .and. &
               s% job% relax_initial_omega_div_omega_crit .and. .not. restart) then
            if (failed('star_surface_omega_crit',ierr)) return
            call star_relax_uniform_omega( &
               id, 2, s% job% new_omega_div_omega_crit, &
               s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, s% job% new_omega_div_omega_crit
         end if

         if (s% job% set_max_dt_to_frac_lifetime) then
            log_m = log10(s% star_mass) ! in Msun units
            log_lifetime = 9.921 - 3.6648*log_m + 1.9697*log_m**2 - 0.9369*log_m**3
            ! Iben & Laughlin (1989) as quoted in H&K (eqn 2.3)
            max_dt = s% job% max_frac_of_lifetime_per_step*secyer*10**(log_lifetime)
            if (max_dt < s% max_timestep) then
               s% max_timestep = max_dt
               write(*, *) 'set_max_dt_to_frac_lifetime: lg(maxdt/secyer)', &
                  log10(s% max_timestep/secyer)
            end if
         end if
         
         if (len_trim(s% job% history_columns_file) > 0) &
            write(*,*) 'read ' // trim(s% job% history_columns_file)
         call star_set_history_columns(id, s% job% history_columns_file, .true., ierr)
         if (failed('star_set_history_columns',ierr)) return
         
         if (len_trim(s% job% profile_columns_file) > 0) &
            write(*,*) 'read ' // trim(s% job% profile_columns_file)
         call star_set_profile_columns(id, s% job% profile_columns_file, .true., ierr)
         if (failed('star_set_profile_columns',ierr)) return
         
         ! print out info about selected non-standard parameter settings
         
         write(*,*) 'net name ' // trim(s% net_name)
         
         if (len_trim(s% extra_terminal_output_file) > 0) &
            write(*,*) 'extra_terminal_output_file: ' // trim(s% extra_terminal_output_file)
         
         if (s% do_element_diffusion) &
            write(*,*) 'do_element_diffusion', s% do_element_diffusion
         
         if (s% lnPgas_flag) &
            write(*,*) 'lnPgas_flag =', s% lnPgas_flag
         
         if (s% E_flag) &
            write(*,*) 'E_flag =', s% E_flag
         
         if (.not. s% L_flag) &
            write(*,*) 'L_flag =', s% L_flag
         
         if (s% v_flag) &
            write(*,*) 'v_flag =', s% v_flag
         
         if (s% rotation_flag) &
            write(*,*) 'rotation_flag =', s% rotation_flag
         
         if (s% mix_factor /= 1d0) &
            write(*,1) 'mix_factor =', s% mix_factor

         if (s% hydro_numerical_jacobian) &
            write(*,*) 'hydro_numerical_jacobian', s% hydro_numerical_jacobian
            
         if (abs(s% tau_base - 2d0/3d0) > 1d-4) &
            write(*,1) 'tau_base', s% tau_base
            
         if (abs(s% tau_factor - 1) > 1d-4) &
            write(*,1) 'tau_factor', s% tau_factor
            
         if (s% eps_grav_factor /= 1) &
            write(*,1) 'eps_grav_factor', s% eps_grav_factor
            
         if (s% dxdt_nuc_factor /= 1) &
            write(*,1) 'dxdt_nuc_factor', s% dxdt_nuc_factor
            
         if (s% which_atm_option /= 'simple_photosphere') &
            write(*,1) 'which_atm_option: ' // trim(s% which_atm_option)
           
         if (s% M_center /= 0) then
            write(*,1) 'xmstar/mstar', s% xmstar/s% mstar
            write(*,1) 'xmstar (g)', s% xmstar
            write(*,1) 'M_center (g)', s% M_center
            write(*,1) 'xmstar/Msun', s% xmstar/Msun
            write(*,1) 'M_center/Msun', s% M_center/Msun
         end if
            
         if (s% R_center /= 0) then
            write(*,1) 'R_center (cm)', s% R_center
            write(*,1) 'R_center/Rsun', s% R_center/Rsun
            write(*,1) 'core density', &
               s% M_center/(4*pi/3*s% R_center*s% R_center*s% R_center)
         end if
            
         if (s% v_center /= 0) &
            write(*,1) 'v_center (cm/s)', s% v_center
            
         if (s% L_center /= 0) &
            write(*,1) 'L_center/Lsun', s% L_center/Lsun
                     
         if (s% opacity_max > 0) &
            write(*,1) 'opacity_max', s% opacity_max
            
         if (s% split_mixing_choice /= 0) &
            write(*,2) 'split_mixing_choice', s% split_mixing_choice
         
         if (s% job% show_net_reactions_info) then
            write(*,'(a)') ' net reactions '
            call show_net_reactions_and_info(s% net_handle, 6, ierr)
            if (failed('show_net_reactions_and_info',ierr)) return
         end if
         
         if (s% job% list_net_reactions) then
            write(*,'(a)') ' net reactions '
            call show_net_reactions(s% net_handle, 6, ierr)
            if (failed('show_net_reactions',ierr)) return
         end if
         
         if (s% hydro_decsol_switch < 10000) then
            if (s% species > s% hydro_decsol_switch) then
               write(*,3) 'use large_mtx_decsol ' // trim(s% large_mtx_decsol), &
                  s% species, s% hydro_decsol_switch
            else
               write(*,3) 'use small_mtx_decsol ' // trim(s% small_mtx_decsol), &
                  s% species, s% hydro_decsol_switch
            end if
         end if
         
         if (s% job% show_net_species_info) then
            write(*,'(a)') ' species'
            do j=1,s% species
               write(*,'(i6,3x,a)') j, chem_isos% name(s% chem_id(j))
            end do
            write(*,*)
         end if
         
         if (s% job% show_eqns_and_vars_names) then
            do i=1,s% nvar
               write(*,*) i, s% nameofvar(i), s% nameofequ(i)
            end do
            write(*,*)
         end if         
         
         write(*,'(a)') ' kappa_file_prefix ' // trim(s% job% kappa_file_prefix)
         write(*,'(a)') ' kappa_lowT_prefix ' // trim(s% job% kappa_lowT_prefix)
         
         write(*,'(a)') '   eos_file_prefix ' // trim(s% job% eos_file_prefix)

         write(*,2) 'OMP_NUM_THREADS', omp_get_max_threads()

         
         contains


         subroutine do_relax_initial_composition(ierr)
            use utils_lib
            integer, intent(out) :: ierr
            real(dp), pointer :: xq(:), xa(:,:)
            integer :: num_pts, num_species, i, iounit
            include 'formats.inc'
            
            write(*,*)
            write(*,1) 'relax_initial_composition'

            iounit = alloc_iounit(ierr)
            if (ierr /= 0) return
            open(unit=iounit, file=trim(s% job% relax_composition_filename), &
                  status='old', action='read', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'open failed', ierr, iounit
               write(*, '(a)') 'failed to open ' // trim(s% job% relax_composition_filename)
               call free_iounit(iounit)
               return
            end if
            read(iounit, *, iostat=ierr) num_pts, num_species
            if (ierr /= 0) then
               close(iounit)
               call free_iounit(iounit)
               write(*, '(a)') 'failed while trying to read 1st line of ' // &
                  trim(s% job% relax_composition_filename)
               return
            end if
            allocate(xq(num_pts), xa(num_species,num_pts))
            do i = 1, num_pts
               read(iounit,*,iostat=ierr) xq(i), xa(1:num_species,i)
               if (ierr /= 0) then
                  close(iounit)
                  call free_iounit(iounit)
                  write(*, '(a)') &
                     'failed while trying to read ' // trim(s% job% relax_composition_filename)
                  write(*,*) 'line', i+1
                  write(*,*) 'perhaps wrong info in 1st line?'
                  write(*,*) '1st line must have num_pts and num_species in that order'
                  deallocate(xq,xa)
                  return
               end if
            end do
            close(iounit)
            call free_iounit(iounit)
            call star_relax_composition( &
               id, s% job% num_steps_to_relax_composition, num_pts, num_species, xa, xq, ierr)
            deallocate(xq,xa)
         end subroutine do_relax_initial_composition
         

      end subroutine do_star_job_controls_after
      

      end module run_star_support
      
      
      
      
