! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      
      implicit none
      
      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 1 ! minutes

      integer, parameter :: restart_info_alloc = 1
      integer, parameter :: restart_info_get = 2
      integer, parameter :: restart_info_put = 3
      
      
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine extras_controls
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         integer :: restart_time, prev_time_used
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         if (.not. restart) then
            call system_clock(time0,clock_rate)
            call alloc_restart_info(s)
         else
            call unpack_restart_info(s)
            call system_clock(restart_time,clock_rate)
            prev_time_used = time1 - time0
            time1 = restart_time
            time0 = time1 - prev_time_used
         end if
         s% xtra1 = log10(abs(b% mtransfer_rate/(Msun/secyer)))  
         write(*,*) "central hydrogen ", s% center_h1 
         if (s% center_h1 >0.25) then 
            s% lxtra2 =.false.
         else 
           s% lxtra2 =.true. 
         end if 
         if (s% center_h1 > 0.0001) then 
            s% lxtra1 = .false.
         else 
           s% lxtra1 = .true.
        end if 
         
        if (safe_log10_cr(s% L(1)/Lsun) < -2.5) then
            s% lxtra3 = .true.
        else
           s% lxtra3 = .false.
       
        end if
           
           if (s% Teff <10000) then
          s% which_atm_option = 'WD_tau_25_tables'
      end if
 
           extras_startup = keep_going
      end function extras_startup
      

      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         ierr = 0
         if (n /= 1) then
            stop 'bad n for data_for_extra_history_columns'
         end if
         dt = dble(time1 - time0) / clock_rate / 60
         names(1) = 'runtime_minutes'
         vals(1) = dt
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call system_clock(time1,clock_rate)
         if (log10(abs(b% mtransfer_rate/(Msun/secyer))) >-15.0 .and.  &
            log10(abs(b% mtransfer_rate_older/(Msun/secyer))) <-15.0) then  

         write(*,*) "begining of mass transfer"
         s% need_to_save_profiles_now = .true.
         s% save_profiles_model_priority = 6   
         end if

       !  write (*,*) "mass tarnsfer rate ",  log10(abs(b% mtransfer_rate)/(Msun/secyer)),&
       !  log10(abs(b% mtransfer_rate_older)/(Msun/secyer))
         if (log10(abs(b% mtransfer_rate)/(Msun/secyer)) <-12.0 .and. log10(abs(b% mtransfer_rate_older)/(Msun/secyer)) >-12.0) then
            write(*,*) "end of mass transfer"
            s% need_to_save_profiles_now = .true.
            s% lxtra4 =.true.
         !   s% xtra1=s %xtra1+1 
            s% save_profiles_model_priority = 7   
         end if          
           
         if ((s% center_h1 < 0.0001) .and. (s% lxtra1.eqv..false.)) then 
            write(*,*) "end of main sequence"
            s% need_to_save_profiles_now = .true.
            s% save_profiles_model_priority = 5   
             s% lxtra1 =.true.
          end if  
        
         if ((s% center_h1 < 0.25) .and.(s% lxtra2.eqv..false.))  then 
            write(*,*) "0.25 hydrogen  main sequence"
            s% need_to_save_profiles_now = .true.
            s% save_profiles_model_priority = 4   
             s% lxtra2 =.true.
         end if  
         
         if ((safe_log10_cr(s% L(1)/Lsun) < -2.5) .and.(s% lxtra3.eqv..false.))  then
            write(*,*) "cooling track : log L =2.5"
            s% need_to_save_profiles_now = .true.
            s% save_profiles_model_priority = 8
             s% lxtra3 =.true.
         end if
         

         if (safe_log10_cr(s% power_h_burn) > 3) then

             s% tol_correction_norm = 3d-4
             s% tol_max_correction = 3d-2
             write(*,*) "relax tolerance"
             else
                s% tol_correction_norm = 3d-5
                s% tol_max_correction = 3d-3
                write(*,*) "normal tolerance"
             end if
         call store_restart_info(s)
       if (s% model_number > 1)     then 
          s% do_element_diffusion = .true.
          write (*,*) "element diffusion on "
      end if 
        if (s% Teff <10000) then 
          s% which_atm_option = 'WD_tau_25_tables'
      end if    
 end function extras_finish_step
     
 
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         dt = dble(time1 - time0) / clock_rate / 60
         if (dt > 10*expected_runtime) then
            write(*,'(/,a30,2f18.6,a,/)') '>>>>>>> EXCESSIVE runtime', &
               dt, expected_runtime, '   <<<<<<<<<  ERROR'
         else
            write(*,'(/,a50,2f18.6,99i10/)') 'runtime, retries, backups, steps', &
               dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
         end if
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring data so can do restarts

      
      subroutine alloc_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_alloc)
      end subroutine alloc_restart_info
      
      
      subroutine unpack_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_get)
      end subroutine unpack_restart_info
      
      
      subroutine store_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_put)
      end subroutine store_restart_info
      
      
      subroutine move_restart_info(s,op)
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg 
         call move_int(time0)
         call move_int(time1)
         
         num_ints = i
         
         i = 0
         ! call move_dbl 
         
         num_dbls = i
         
         if (op /= restart_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (restart_info_get)
               dbl = s% extra_work(i)
            case (restart_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            include 'formats'
            i = i+1
            select case (op)
            case (restart_info_get)
               !write(*,3) 'restore int', i, s% extra_iwork(i)
               int = s% extra_iwork(i)
            case (restart_info_put)
               !write(*,3) 'save int', i, int
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (restart_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (restart_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_restart_info
      
      


      end module run_star_extras
      
