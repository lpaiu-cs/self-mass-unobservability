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

      module read_model

      use star_private_def
      use const_def
      use utils_lib, only: is_bad_num

      implicit none
      
      integer, parameter :: bit_for_zams_file = 0
      integer, parameter :: bit_for_lnPgas = 1 ! includes lnPgas variables in place of lnd
      integer, parameter :: bit_for_2models = 2 ! file contains two consecutive models
      integer, parameter :: bit_for_velocity = 3 ! includes velocity variables
      integer, parameter :: bit_for_rotation = 4 ! includes rotation variables
      integer, parameter :: bit_for_conv_vel = 5 ! includes convection velocities
      integer, parameter :: bit_for_E = 6 ! includes internal energy variables in place of lnT
      integer, parameter :: bit_for_no_L = 7 ! no luminosity variables

      integer, parameter :: mesa_zams_file_type = 2**bit_for_zams_file
      
      character (len=100000) :: buf


      contains
      
      
      subroutine finish_load_model(s, ierr)
         use hydro_vars, only: set_vars
         use star_utils, only: use_xh_to_update_i_rot_and_j_rot, &
            set_m_and_dm, set_dm_bar, total_angular_momentum, save_for_d_dt
         use report, only: do_report
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, j

         include 'formats'
         
         s% brunt_B(1:s% nz) = 0 ! temporary brunt_B for set_vars
         
         if (s% rotation_flag) then
            call set_m_and_dm(s)
            call set_dm_bar(s% nz, s% dm, s% dm_bar)
            call use_xh_to_update_i_rot_and_j_rot(s)
            s% total_angular_momentum = total_angular_momentum(s)
            !write(*,1) 'finish_load_model log J tot', log10(s% total_angular_momentum)


            if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
               do k=1,s% nz
                  write(*,3) 'lnr', s% model_number, k, s% xh(s% i_lnR, k)
                  write(*,3) 'i_rot', s% model_number, k, s% i_rot(k)
                  write(*,3) 'j_rot', s% model_number, k, s% j_rot(k)
                  write(*,3) 'omega', s% model_number, k, s% omega(k)
               end do
            end if


         end if
         
         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do j=1,s% species
               write(*,4) 'finish_load_model before set_vars xa(j)', &
                  s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
            end do
         end if
         
         do k=1,s% nz
            s% T_start(k) = -1d99
            s% r_start(k) = -1d99
            s% v_start(k) = -1d99
            s% csound_start(k) = -1d99
         end do
         
         call save_for_d_dt(s)
         
         call set_vars(s, s% dt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*, *) 'load_model: set_vars ierr'
            end if
            write(*,*) 'finish_load_model: failed in set_vars'
            return
         end if      
         
         if (s% trace_k > 0 .and. s% trace_k <= s% nz) then
            do j=1,s% species
               write(*,4) 'finish_load_model after set_vars xa(j)', &
                  s% model_number, s% trace_k, j, s% xa(j,s% trace_k)
            end do
         end if

         ! actual brunt_B is set as part of do_report
         call do_report(s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'load_model: do_report ierr'
            write(*,*) 'finish_load_model: failed in do_report'
            return
         end if

      end subroutine finish_load_model
      
      
      subroutine set_zero_age_params(s)
         type (star_info), pointer :: s
         s% model_number = 0
         s% time = 0
         s% dt = 0
         s% dt_old = 0
         s% total_num_newton_iterations = 0
         s% num_retries = 0
         s% num_backups = 0
         s% mesh_call_number = 0
         s% hydro_call_number = 0
         s% diffusion_call_number = 0
         s% model_number_for_last_retry = 0
         s% Teff = -1 ! need to calculate it
      end subroutine set_zero_age_params


      subroutine read1_model( &
            s, species, nvar_hydro, nz, iounit, &
            xh, xa, q, dq, omega, lnT, conv_vel, perm, ierr)
         use utils_lib, only: is_bad_num
         use star_utils, only: set_qs
         use crlibm_lib, only: str_to_vector
         use chem_def   
         type (star_info), pointer :: s
         integer, intent(in) :: species, nvar_hydro, nz, iounit, perm(:)
         real(dp), dimension(:,:), intent(out) :: xh, xa
         real(dp), dimension(:), intent(out) :: q, dq, omega, lnT, conv_vel
         integer, intent(out) :: ierr

         integer :: i, j, n, i_lnd, i_lnT, i_lnR, i_lum, &
            i_lnPgas, i_E, i_v, ii
         real(dp), target :: vec_ary(species + nvar_hydro + 10)
         real(dp), pointer :: vec(:)
         integer :: nvec
            
         include 'formats'
         
         ierr = 0
         vec => vec_ary
      
         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_lnPgas = s% i_lnPgas
         i_E = s% i_E
         i_v = s% i_v
         n = species + nvar_hydro + 1 ! + 1 is for dq
         if (i_lnPgas /= 0) n = n+1 ! read both lnPgas and lnd 
         if (s% rotation_flag) n = n+1 ! read omega 
         if (s% have_previous_conv_vel) n = n+1 ! read conv_vel 

!$omp critical (read1_model_loop)
! make this a critical section to so don't have to dynamically allocate buf 
         do i = 1, nz
            read(iounit,'(a)',iostat=ierr) buf
            if (ierr == 0) call str_to_vector(buf, vec, nvec, ierr)
            if (ierr /= 0) exit 
            j = int(vec(1))           
            if (j /= i) then
               ierr = -1
               write(*, *) 'error in reading model data   j /= i'
               write(*, *) 'species', species
               write(*, *) 'j', j
               write(*, *) 'i', i
               exit
            end if
            j = 1
            j=j+1; if (i_lnd /= 0) xh(i_lnd,i) = vec(j)
            j=j+1
            if (i_lnT /= 0) then
               xh(i_lnT,i) = vec(j)
            else
               lnT(i) = vec(j)
            end if
            j=j+1; xh(i_lnR,i) = vec(j)
            if (i_lum /= 0) then
               j=j+1; xh(i_lum,i) = vec(j)
            end if
            j=j+1; dq(i) = vec(j)
            if (i_lnPgas /= 0) then
               j=j+1; xh(i_lnPgas,i) = vec(j)
            end if
            if (i_E /= 0) then
               j=j+1; xh(i_E,i) = vec(j)
            end if
            if (i_v /= 0) then
               j=j+1; xh(i_v,i) = vec(j)
            end if
            if (s% rotation_flag) then
               j=j+1; omega(i) = vec(j)
            end if
            if (s% have_previous_conv_vel) then
               j=j+1; conv_vel(i) = vec(j)
            end if
            if (j+species > nvec) then
               ierr = -1
               write(*, *) 'error in reading model data  j+species > nvec'
               write(*, *) 'j+species', j+species
               write(*, *) 'nvec', nvec
               write(*, *) 'j', j
               write(*, *) 'species', species
               exit
            end if
            do ii=1,species
               xa(perm(ii),i) = vec(j+ii)
            end do
         end do
!$omp end critical (read1_model_loop)

         if (ierr /= 0) then
            !if (s% report_ierr) &
               write(*,*) 'read1_model_loop failed'
            return
         end if
         
         if (i_lum /= 0) then
            s% prev_Lmax = maxval(abs(xh(i_lum,1:nz)))
         else
            s% prev_Lmax = 0
         end if
         
         call set_qs(nz, q, dq, ierr)
         if (ierr /= 0) then
            !if (s% report_ierr) &
               write(*,*) 'set_qs failed in read1_model sum(dq)', sum(dq(1:nz))
            return
         end if
         
      end subroutine read1_model

      
      subroutine do_read_saved_model(s, filename, ierr)
         use utils_lib
         use utils_def
         use chem_def
         use net, only: set_net
         use alloc, only: set_var_info, &
            free_star_info_arrays, allocate_star_info_arrays, set_chem_names
         use star_utils, only: yrs_for_init_timestep, set_phase_of_evolution
         type (star_info), pointer :: s
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr   
         
         integer :: iounit, n, i, t, file_type, year_month_day_when_created, nz, species, nvar
         logical :: do_read_prev
         real(dp) :: initial_mass, initial_z, initial_y, tau_factor, mixing_length_alpha
         character (len=strlen) :: buffer, string, message
         character (len=net_name_len) :: net_name
      	character(len=iso_name_length), pointer :: names(:) ! (species)
      	integer, pointer :: perm(:) ! (species)
      	
      	include 'formats'
         
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         open(unit=iounit, file=trim(filename), status='old', action='read', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'open failed', ierr, iounit
            write(*, '(a)') 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if
         
         ! use token to get file_type so can have comments at start of file
         n = 0
         i = 0
         t = token(iounit, n, i, buffer, string)
         if (t == eof_token) then
            write(*, '(a)') 'failed to find file type at start of ' // trim(filename)
            call free_iounit(iounit)
            return
         end if
         if (t /= name_token) then
            write(*, '(a)') 'failed to find file type at start of ' // trim(filename)
            call free_iounit(iounit)
            return
         end if
         read(string,fmt=*,iostat=ierr) file_type
         if (ierr /= 0) then
            write(*, '(a)') 'failed to find file type at start of ' // trim(filename)
            call free_iounit(iounit)
            return
         end if 
         
         read(iounit, *, iostat=ierr) ! skip the blank line after the file type
         if (ierr /= 0) then
            call free_iounit(iounit)
            return
         end if
         
         s% model_number = 0
         s% star_age = 0
         s% xmstar = -1
         
         tau_factor = s% tau_factor
         mixing_length_alpha = s% mixing_length_alpha
         s% free_fall_mdot_Bondi = 0
         s% free_fall_entropy = 0

         call read_properties(iounit, &
            net_name, species, nz, year_month_day_when_created, &
            initial_mass, initial_z, initial_y, mixing_length_alpha, &
            s% model_number, s% star_age, &
            tau_factor, s% xmstar, s% R_center, s% L_center, s% v_center, &
            s% free_fall_mdot_Bondi, s% free_fall_entropy, ierr)
      
         if (ierr /= 0 .or. initial_mass < 0 .or. nz < 0 &
               .or. initial_z < 0 .or. species < 0) then
            ierr = -1
            write(*, *) 'do_read_model: missing required properties'
            write(*,*) 'initial_mass', initial_mass
            write(*,*) 'initial_z', initial_z
            write(*,*) 'nz', nz
            write(*,*) 'species', species
            call free_iounit(iounit)
            return
         end if
         
         s% time = s% star_age*secyer
         if (s% free_fall_mdot_Bondi > 0) then
            s% using_free_fall_surface_PT = .true.
            write(*,*) 'using_free_fall_surface_PT'
         end if
         
         if (abs(tau_factor - s% tau_factor) > tau_factor*1d-9) then
            write(*,*)
            write(*,1) 'WARNING: model saved with tau_factor =', tau_factor
            write(*,1) 'but current setting for tau_factor =', s% tau_factor
            write(*,*)
         end if
         
         if (abs(mixing_length_alpha - s% mixing_length_alpha) > mixing_length_alpha*1d-9) then
            write(*,*)
            write(*,1) 'WARNING: model saved with mixing_length_alpha =', mixing_length_alpha
            write(*,1) 'but current setting for mixing_length_alpha =', s% mixing_length_alpha
            write(*,*)
         end if
         
         s% net_name = trim(net_name)
         s% species = species
         s% lnPgas_flag = BTEST(file_type, bit_for_lnPgas)
         s% E_flag = BTEST(file_type, bit_for_E)
         s% L_flag = .not. BTEST(file_type, bit_for_no_L)
         s% v_flag = BTEST(file_type, bit_for_velocity)
         s% rotation_flag = BTEST(file_type, bit_for_rotation)
         s% have_previous_conv_vel = BTEST(file_type, bit_for_conv_vel)
         s% initial_z = initial_z

         s% mstar = initial_mass*Msun
         if (s% xmstar < 0) then
            s% M_center = 0
            s% xmstar = s% mstar
         else
            s% M_center = s% mstar - s% xmstar
         end if

         call set_net(s, s% net_name, ierr)
         if (ierr /= 0) then
            write(*,*) &
               'do_read_saved_model failed in set_net for net_name = ' // trim(s% net_name)
            call free_iounit(iounit)
            return
         end if

         call set_var_info(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_read_saved_model failed in set_var_info'
            call free_iounit(iounit)
            return
         end if

         ! fixup chem names now that have nvar_hydro
         call set_chem_names(s)

         s% nz = nz
         call free_star_info_arrays(s)
         call allocate_star_info_arrays(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_read_saved_model failed in allocate_star_info_arrays'
            call free_iounit(iounit)
            return
         end if

         allocate(names(species), perm(species))
         call get_chem_col_names(s, iounit, species, names, perm, ierr)
         if (ierr /= 0) then
            deallocate(names, perm)
            write(*,*) 'do_read_saved_model failed in get_chem_col_names'
            call free_iounit(iounit)
            return
         end if
         
         nvar = s% nvar
         call read1_model( &
               s, s% species, s% nvar_hydro, nz, iounit, &
               s% xh, s% xa, s% q, s% dq, s% omega, s% lnT, &
               s% prev_conv_vel_from_file, perm, ierr)
         deallocate(names, perm)
         if (ierr /= 0) then
            write(*,*) 'do_read_saved_model failed in read1_model'
            call free_iounit(iounit)
            return
         end if
         
         s% use_previous_conv_vel_from_file = s% have_previous_conv_vel

         do_read_prev = BTEST(file_type, bit_for_2models)
         if (ierr == 0) then
            if (do_read_prev) then
               call read_prev
            else
               s% generations = 1
            end if
         end if
                  
         close(iounit)
         call free_iounit(iounit)   

         !call dump_4_debug(s); stop 'dump_4_debug'
         
         
         contains
         
         
         subroutine read_prev
            integer :: k
            
            do k=1, 3
               read(iounit, *, iostat=ierr) 
               if (ierr /= 0) return
            end do
            call read_prev_properties
            if (ierr /= 0) return
            
            ! we do read_prev_properties to set initial timestep,
            ! but we don't use the previous model
            ! because we need to have other info about that isn't saved
            ! such as conv_vel and mixing_type
            
            s% generations = 1
                  
         end subroutine read_prev      
         
         
         subroutine read_prev_properties
            character (len=132) :: line
            real(dp) :: tmp, skip_val
            integer :: i
            include 'formats'
            ierr = 0
            s% dt = -1
            s% mstar_old = -1
            s% dt_next = -1
            s% nz_old = -1

            do ! until reach a blank line
               read(iounit, fmt='(a)', iostat=ierr) line
               if (ierr /= 0) return
               
               if (len_trim(line) == 0) exit ! blank line
               
               if (match_keyword('previous n_shells', line, tmp)) then
                  s% nz_old = int(tmp)
                  cycle
               end if
               
               if (match_keyword('timestep (seconds)', line, s% dt)) then
                  cycle
               end if
               
               if (match_keyword('previous mass (grams)', line, s% mstar_old)) then
                  cycle
               end if
               
               if (match_keyword('dt_next (seconds)', line, s% dt_next)) then
                  cycle
               end if

               if (match_keyword('year_month_day_when_created', line, skip_val)) cycle

            end do
            if (s% dt < 0) then
               ierr = -1
               write(*, *) 'missing dt for previous model'
            end if
            if (s% mstar_old < 0) then
               ierr = -1
               write(*, *) 'missing mstar_old for previous model'
            end if
            if (s% dt_next < 0) then
               ierr = -1
               write(*, *) 'missing dt_next for previous model'
            end if
            
         end subroutine read_prev_properties


      end subroutine do_read_saved_model
      
      
      subroutine do_read_saved_model_number(fname, model_number, ierr)
         character (len=*), intent(in) :: fname
         integer, intent(inout) :: model_number
         integer, intent(out) :: ierr
         character (len=strlen) :: net_name
         integer :: species, n_shells, &
            year_month_day_when_created
         real(dp) :: m_div_msun, initial_z, &
            mixing_length_alpha, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy
         call do_read_saved_model_properties(fname, &
            net_name, species, n_shells, year_month_day_when_created, &
            m_div_msun, initial_z, mixing_length_alpha, &
            model_number, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy, &
            ierr)
      end subroutine do_read_saved_model_number
      
      
      subroutine do_read_saved_model_age(fname, star_age, ierr)
         character (len=*), intent(in) :: fname
         real(dp), intent(inout) :: star_age
         integer, intent(out) :: ierr
         character (len=strlen) :: net_name
         integer :: species, n_shells, model_number, &
            year_month_day_when_created
         real(dp) :: m_div_msun, initial_z, &
            mixing_length_alpha, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy
         call do_read_saved_model_properties(fname, &
            net_name, species, n_shells, year_month_day_when_created, &
            m_div_msun, initial_z, mixing_length_alpha, &
            model_number, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy, &
            ierr)
      end subroutine do_read_saved_model_age
      
      
      subroutine do_read_saved_model_properties(fname, &
            net_name, species, n_shells, year_month_day_when_created, &
            m_div_msun, initial_z, mixing_length_alpha, &
            model_number, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy, &
            ierr)
         use utils_lib
         character (len=*), intent(in) :: fname
         character (len=*), intent(inout) :: net_name
         integer, intent(inout) :: species, n_shells, &
            year_month_day_when_created, model_number
         real(dp), intent(inout) :: m_div_msun, initial_z, &
            mixing_length_alpha, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy
         integer, intent(out) :: ierr
         integer :: iounit
         real(dp) :: initial_y
         ierr = 0
         iounit = alloc_iounit(ierr); if (ierr /= 0) return
         open(iounit, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open ' // trim(fname)
            call free_iounit(iounit)
            return
         end if                  
         read(iounit, *, iostat=ierr) 
         if (ierr /= 0) then
            close(iounit)
            call free_iounit(iounit)
            return
         end if
         read(iounit, *, iostat=ierr) 
         if (ierr /= 0) then
            close(iounit)
            call free_iounit(iounit)
            return
         end if
         call read_properties(iounit, &
            net_name, species, n_shells, year_month_day_when_created, &
            m_div_msun, initial_z, initial_y, mixing_length_alpha, &
            model_number, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy, &
            ierr)
         close(iounit)
         call free_iounit(iounit)
      end subroutine do_read_saved_model_properties
      
      
      subroutine read_properties(iounit, &
            net_name, species, n_shells, year_month_day_when_created, &
            m_div_msun, initial_z, initial_y, mixing_length_alpha, &
            model_number, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy, &
            ierr)
         integer, intent(in) :: iounit
         character (len=*), intent(inout) :: net_name
         integer, intent(inout) :: species, n_shells, &
            year_month_day_when_created, model_number
         real(dp), intent(inout) :: m_div_msun, initial_z, initial_y, &
            mixing_length_alpha, star_age, &
            tau_factor, xmstar, R_center, L_center, v_center, &
            free_fall_mdot_Bondi, free_fall_entropy
         integer, intent(out) :: ierr
         character (len=132) :: line
         real(dp) :: tmp
         ierr = 0
         do ! until reach a blank line
            read(iounit, fmt='(a)', iostat=ierr) line
            if (ierr /= 0) return
            if (len_trim(line) == 0) return ! blank line
            if (match_keyword_for_string('net_name', line, net_name)) then; cycle; end if
            if (match_keyword('species', line, tmp)) then; species = int(tmp); cycle; end if
            if (match_keyword('n_shells', line, tmp)) then; n_shells = int(tmp); cycle; end if
            if (match_keyword('model_number', line, tmp)) then; model_number = int(tmp); cycle; end if
            if (match_keyword('M/Msun', line, m_div_msun)) cycle
            if (match_keyword('star_age', line, star_age)) cycle
            if (match_keyword('initial_z', line, initial_z)) cycle
            if (match_keyword('initial_y', line, initial_y)) cycle
            if (match_keyword('mixing_length_alpha', line, mixing_length_alpha)) cycle
            if (match_keyword('tau_factor', line, tau_factor)) cycle
            if (match_keyword('xmstar', line, xmstar)) cycle
            if (match_keyword('R_center', line, R_center)) cycle
            if (match_keyword('L_center', line, L_center)) cycle
            if (match_keyword('v_center', line, v_center)) cycle
            if (match_keyword('year_month_day_when_created', line, tmp)) then
               year_month_day_when_created = int(tmp); cycle; end if
            if (match_keyword('tau_factor', line, tmp)) cycle
            if (match_keyword('tau_photosphere', line, tmp)) cycle
            if (match_keyword('free_fall_mdot_Bondi', line, free_fall_mdot_Bondi)) cycle
            if (match_keyword('free_fall_entropy', line, free_fall_entropy)) cycle
         end do
      end subroutine read_properties
         
         
      logical function match_keyword(key, txt, value)
         ! returns true if leading non-blank part of txt is same as key.
         ! i.e., skips leading blanks in txt before testing equality.
         character (len=*), intent(in) :: key, txt
         real(dp), intent(inout) :: value
         integer :: i, j, k, ierr
         i = len(key)
         k = len(txt)
         j = 1
         do while (j <= k .and. txt(j:j) == ' ')
            j = j+1
         end do
         match_keyword = (txt(j:j+i-1) == key)
         ierr = 0
         if (match_keyword) then
            read(txt(j+i:k), fmt=*, iostat=ierr) value
            if (ierr /= 0) match_keyword = .false.
         end if
      end function match_keyword
         
         
      logical function match_keyword_for_string(key, txt, value)
         ! returns true if leading non-blank part of txt is same as key.
         ! i.e., skips leading blanks in txt before testing equality.
         character (len=*), intent(in) :: key, txt
         character (len=*), intent(inout) :: value
         integer :: i, j, k, str_len
         logical, parameter :: dbg = .false.
         i = len(key)
         k = len(txt)
         j = 1
         do while (j <= k .and. txt(j:j) == ' ')
            j = j+1
         end do
         match_keyword_for_string = (txt(j:j+i-1) == key)
         if (.not. match_keyword_for_string) return
         if (dbg) then
            write(*,*) 'matching ' // trim(key)
            write(*,*) 'txt ' // trim(txt)
         end if
         j = j+i
         do while (j <= k .and. txt(j:j) == ' ')
            j = j+1
         end do
         if (j > k) then
            match_keyword_for_string = .false.
            if (dbg) write(*,*) 'j > k'
            return
         end if
         if (txt(j:j) /= '''') then
            match_keyword_for_string = .false.
            if (dbg) write(*,*) 'no leading quote'
            return
         end if
         j = j+1
         i = 1
         str_len = len(value)
         do while (j <= k .and. txt(j:j) /= '''')
            value(i:i) = txt(j:j)
            i = i+1
            j = j+1
         end do
         do while (i <= str_len)
            value(i:i) = ' '
            i = i+1
         end do
         if (dbg) write(*,*) 'value <' // trim(value) // ">"
      end function match_keyword_for_string
      
      
      
      subroutine get_chem_col_names(s, iounit, species, names, perm, ierr)
         use chem_def, only: iso_name_length
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         integer, intent(in) :: iounit, species
      	character(len=iso_name_length), intent(out) :: names(species)
      	integer, intent(out) :: perm(species)
         integer, intent(out) :: ierr
         
         character (len=50000) :: buffer
         character (len=20) :: string
         integer :: n, i, j1, j2, str_len, l, indx, j, num_found
         
         ierr = 0
         read(iounit,fmt='(a)',iostat=ierr) buffer
         if (ierr /= 0) return
         
         n = len_trim(buffer)         
         i = 0
         num_found = 0
       token_loop: do ! have non-empty buffer
            i = i+1
            if (i > n) then
               write(*,*) 'get_chem_col_names: failed to find all of the names'
               ierr = -1
               return
            end if
            if (buffer(i:i) == char(9)) cycle token_loop ! skip tabs
            select case(buffer(i:i))
               case (' ')
                  cycle token_loop ! skip spaces
               case default
                  j1 = i; j2 = i
                  name_loop: do
                     if (i+1 > n) exit
                     if (buffer(i+1:i+1) == ' ') exit
                     if (buffer(i+1:i+1) == '(') exit
                     if (buffer(i+1:i+1) == ')') exit
                     if (buffer(i+1:i+1) == ',') exit
                     i = i+1
                     j2 = i
                  end do name_loop
                  str_len = len(string)
                  l = j2-j1+1
                  if (l > str_len) then
                     l = str_len
                     j2 = l+j1-1
                  end if
                  string(1:l) = buffer(j1:j2)
                  do j = l+1, str_len
                     string(j:j) = ' '
                  end do

                  indx = chem_get_iso_id(string)
                  
                  if (indx > 0) then
                     num_found = num_found+1
                     names(num_found) = trim(string)
                     perm(num_found) = s% net_iso(indx)
                     !write(*,*) trim(string), num_found, perm(num_found)
                     if (num_found == species) return
                  end if

            end select
         end do token_loop
      
      end subroutine get_chem_col_names
      
      
      subroutine dump_4_debug(s)
         type (star_info), pointer :: s
         integer :: k, j
         include 'formats'
         write(*,*) 's% E_flag', s% E_flag
         write(*,*) 's% lnPgas_flag', s% lnPgas_flag
         write(*,*) 's% v_flag ', s% v_flag 
         write(*,*) 's% rotation_flag', s% rotation_flag
         write(*,*) 's% have_previous_conv_vel', s% have_previous_conv_vel
         write(*,2) 's% species', s% species
         write(*,2) 's% nz', s% nz
         write(*,1) 's% star_mass', s% star_mass
         write(*,1) 's% star_age', s% star_age
         write(*,1) 's% initial_mass', s% initial_mass
         write(*,1) 's% initial_z', s% initial_z
         write(*,1) 's% xmstar', s% xmstar
         write(*,1) 's% R_center', s% R_center
         write(*,1) 's% L_center', s% L_center
         write(*,1) 's% v_center', s% v_center
         write(*,1) 's% tau_factor', s% tau_factor
         write(*,*) 's% net_name ', trim(s% net_name)
         write(*,2) 's% nz_old', s% nz_old
         write(*,1) 's% mstar_old', s% mstar_old
         write(*,1) 's% dt', s% dt
         write(*,1) 's% dt_next', s% dt_next
         write(*,1) 's% prev_Lmax', s% prev_Lmax
         do k=1,s% nz
            if (s% E_flag) then
               write(*,2) 's% xh(i_E,k)', k, s% xh(s% i_E,k)
            else
               write(*,2) 's% xh(i_lnT,k)', k, s% xh(s% i_lnT,k)
            end if
            write(*,2) 's% xh(i_lnR,k)', k, s% xh(s% i_lnR,k)
            if (s% L_flag) then
               write(*,2) 's% xh(i_lum,k)', k, s% xh(s% i_lum,k)
            end if
            write(*,2) 's% dq(k)', k, s% dq(k)
            if (s% lnPgas_flag) then
               write(*,2) 's% xh(i_lnPgas,k)', k, s% xh(s% i_lnPgas,k)
            else
               write(*,2) 's% xh(i_lnd,k)', k, s% xh(s% i_lnd,k)
            end if
            if (s% v_flag) write(*,2) 's% xh(i_v,k)', k, s% xh(s% i_v,k)
            if (s% rotation_flag) write(*,2) 's% omega(k)', k, s% omega(k)
            if (s% have_previous_conv_vel) &
               write(*,2) 's% conv_vel(k)', k, s% conv_vel(k)
            do j=1,s% species
               write(*,3) 's% xa(j,k)', j, k, s% xa(j,k)
            end do
         end do
      end subroutine dump_4_debug
      

      end module read_model
