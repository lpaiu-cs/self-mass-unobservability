! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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

      module write_model

      use star_private_def
      use const_def
      use read_model

      implicit none
      

      contains
      
      
     
      subroutine do_write_model(id, filename, ierr)
         use utils_lib  
         use chem_def
         integer, intent(in) :: id
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr 
             
         integer :: iounit, i, k, nvar_hydro, nz, species, file_type
         integer, pointer :: chem_id(:)
         type (star_info), pointer :: s
         logical :: lnPgas_flag, E_flag, v_flag, L_flag, &
            prev_flag, rotation_flag, write_conv_vel
         integer :: time_vals(8)
         
         1 format(a32, 2x, 1pd26.16)
         11 format(a32, 2x, 1pd26.16, 2x, a, 2x, 99(1pd26.16))
         2 format(a32, 2x, i30)
         3 format(a32, 3x, a8)
         4 format(a32, 3x, a)

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         chem_id => s% chem_id
         nvar_hydro = s% nvar_hydro
         nz = s% nz
         lnPgas_flag = s% lnPgas_flag
         E_flag = s% E_flag
         v_flag = s% v_flag
         L_flag = s% L_flag
         rotation_flag = s% rotation_flag
         write_conv_vel = L_flag .and. s% have_previous_conv_vel .and. &
            s% min_T_for_acceleration_limited_conv_velocity < 1d12
         species = s% species
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         open(unit=iounit, file=trim(filename), action='write', status='replace')
         write(iounit,'(a)') '! note: initial lines of file can contain comments'
         write(iounit,'(a)') '!'
         prev_flag = (s% nz_old == s% nz .and. s% generations > 1)
         file_type = 0
         if (lnPgas_flag) file_type = file_type + 2**bit_for_lnPgas
         if (E_flag) file_type = file_type + 2**bit_for_E
         if (prev_flag) file_type = file_type + 2**bit_for_2models
         if (v_flag) file_type = file_type + 2**bit_for_velocity
         if (rotation_flag) file_type = file_type + 2**bit_for_rotation
         if (write_conv_vel) file_type = file_type + 2**bit_for_conv_vel
         if (.not. L_flag) file_type = file_type + 2**bit_for_no_L
         write(iounit, '(i14)', advance='no') file_type
         write(iounit,'(a)',advance='no') ' -- model for mesa/star'
         if (BTEST(file_type, bit_for_lnPgas)) &
            write(iounit,'(a)',advance='no') ', with lnPgas variables'
         if (BTEST(file_type, bit_for_E)) &
            write(iounit,'(a)',advance='no') ', with internal energy variables'
         if (BTEST(file_type, bit_for_velocity)) &
            write(iounit,'(a)',advance='no') ', with velocity variables (v)'
         if (BTEST(file_type, bit_for_rotation)) &
            write(iounit,'(a)',advance='no') ', with angular velocity variables (omega)'
         if (BTEST(file_type, bit_for_conv_vel)) &
            write(iounit,'(a)',advance='no') ', with convection velocities (conv_vel)'
         if (BTEST(file_type, bit_for_no_L)) then
            write(iounit,'(a)',advance='no') ', without luminosity variables'
            write(iounit,'(a)',advance='no') &
               '. cgs units. lnd=ln(density), lnT=ln(temperature), lnR=ln(radius)'
         else
            write(iounit,'(a)',advance='no') &
               '. cgs units. lnd=ln(density), lnT=ln(temperature), lnR=ln(radius), L=luminosity'
         end if
         write(iounit,'(a)',advance='no') &
            ', dq=fraction of total mstar in cell; remaining cols are mass fractions.'
         write(iounit,*)
         ! write property list
         write(iounit, *) ! blank line before start of property list
         write(iounit, 1) 'M/Msun', s% star_mass
         write(iounit, 2) 'model_number', s% model_number
         write(iounit, 1) 'star_age', s% star_age
         write(iounit, 1) 'initial_z', s% initial_z
         write(iounit, 2) 'n_shells', nz
         write(iounit, 4) 'net_name', "'" // trim(s% net_name) // "'"
         write(iounit, 2) 'species', species
         if (s% M_center /= 0) then
            write(iounit, 11) 'xmstar', s% xmstar, &
               '! mass above core (g).  mass of core: grams, Msun:', &
               s% M_center, s% M_center/Msun
         end if
         if (s% R_center /= 0) then
            write(iounit, 11) 'R_center', s% R_center, &
               '! radius of core (cm).  R/Rsun, avg core density (g/cm^3):', &
                  s% R_center/Rsun, s% M_center/(4*pi/3*pow3(s% R_center))
         end if
         if (s% v_center /= 0) then
            write(iounit, 11) 'v_center', s% v_center, &
               '! velocity of outer edge of core (cm/s)'
         end if
         if (s% L_center /= 0) then
            write(iounit, 11) 'L_center', s% L_center, &
               '! luminosity of core (erg/s). L/Lsun, avg core eps (erg/g/s):', &
                  s% L_center/Lsun, s% L_center/max(1d0,s% M_center)
         end if
         if (s% tau_factor /= 1) then
            write(iounit, 1) 'tau_factor', s% tau_factor
         end if
         if (s% using_free_fall_surface_PT) then
            write(iounit, 1) 'free_fall_mdot_Bondi', s% free_fall_mdot_Bondi
            write(iounit, 1) 'free_fall_entropy', s% free_fall_entropy
         end if
         write(iounit, *) ! blank line for end of property list

         call header
         do k=1, nz
            if (ierr /= 0) exit
            write(iounit, fmt='(i5, 1x)', advance='no') k
            call write1(s% lnd(k),ierr); if (ierr /= 0) exit
            call write1(s% lnT(k),ierr); if (ierr /= 0) exit
            call write1(s% lnR(k),ierr); if (ierr /= 0) exit
            if (L_flag) then
               call write1(s% L(k),ierr); if (ierr /= 0) exit
            end if
            call write1(s% dq(k),ierr); if (ierr /= 0) exit
            if (lnPgas_flag) then
               call write1(s% lnPgas(k),ierr); if (ierr /= 0) exit
            end if
            if (E_flag) then
               call write1(s% energy(k),ierr); if (ierr /= 0) exit
            end if
            if (v_flag) then
               call write1(s% v(k),ierr); if (ierr /= 0) exit
            end if
            if (rotation_flag) then
               call write1(s% omega(k),ierr); if (ierr /= 0) exit
            end if
            if (write_conv_vel) then
               call write1(s% conv_vel(k),ierr); if (ierr /= 0) exit
            end if
            do i=1, species
               call write1(s% xa(i,k),ierr); if (ierr /= 0) exit
            end do
            write(iounit, *)
         end do

         if (prev_flag) then         
            write(iounit, *)
            write(iounit, '(a)') '        previous model'
            ! write prev properties
            write(iounit, *)
            write(iounit, 2) 'previous n_shells', s% nz_old
            write(iounit, 1, advance='no') 'previous mass (grams)'
            call write1_eol(s% mstar_old,ierr)
            write(iounit, 1, advance='no') 'timestep (seconds)'
            call write1_eol(s% dt,ierr)
            write(iounit, 1, advance='no') 'dt_next (seconds)'
            call write1_eol(s% dt_next_unclipped,ierr)
            write(iounit, *)
         end if
         close(iounit)
         call free_iounit(iounit)   
         
         contains
         
         subroutine write1_eol(val,ierr)
            real(dp), intent(in) :: val
            integer, intent(out) :: ierr
            call write1(val,ierr)
            write(iounit,*)
         end subroutine write1_eol
         
         subroutine write1(val,ierr)
            use crlibm_lib, only: double_to_str_1pd26pt16
            real(dp), intent(in) :: val
            integer, intent(out) :: ierr
            integer, parameter :: str_len = 26
            character (len=str_len) :: string
            ierr = 0
            call double_to_str_1pd26pt16(val,string)
            write(iounit, fmt='(a26,1x)', advance='no') string 
         end subroutine write1
         
         subroutine header
            write(iounit, fmt='(10x, a9, 1x, 99(a26, 1x))', advance='no') &
                  'lnd', 'lnT', 'lnR'
            if (L_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'L'
            write(iounit, fmt='(a26, 1x)', advance='no') 'dq'
            if (lnPgas_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'lnPgas'
            if (E_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'energy'
            if (v_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'v'
            if (rotation_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'omega'
            if (write_conv_vel) &
               write(iounit, fmt='(a26, 1x)', advance='no') 'conv_vel'
            do i=1, species
               write(iounit, fmt='(a26, 1x)', advance='no') chem_isos% name(chem_id(i))
            end do
            write(iounit, *)
         end subroutine header
         
      end subroutine do_write_model


      end module write_model
