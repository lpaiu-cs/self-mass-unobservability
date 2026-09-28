! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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


! when lnPgas_flag is true, convert partials for f(lnd,lnT) by
! f = f(lnPgas(lnd,lnT),lnT)
! df_dlnPgas_const_T = df_dlnd_const_T*dlnRho_dlnPgas_const_T
! df_dlnT_const_Pgas = df_dlnT_const_Rho + df_dlnd_const_T*sdlnRho_dlnT_const_Pgas

! when E_flag is true, convert partials for f(lnd,lnT) by
! f = f(lnd,lnT(lnd,lnE))
! df_dlnd_const_E = df_dlnd_const_T + df_dlnT_const_Rho*dlnT_dlnd_const_E
! df_dlnE_const_Rho = df_dlnT_const_Rho*dlnT_dlnE_const_Rho
! dlnT_dlnd_const_E = -Rho*dE_dRho/(Cv*T)
! dlnT_dlnE_const_Rho = energy/(Cv*T)


      module hydro_eqns

      use star_private_def
      use const_def
      use utils_lib, only: is_bad_num
      use star_utils, only: em1, e00, ep1
      use hydro_vars, only: eval_dlnPdm_qhse, do_dlnP_dm
      use crlibm_lib
      

      implicit none

      real(dp), parameter :: one = 1, zero = 0
      
      logical, parameter :: dbg = .false.   
      
      integer, parameter :: dbg_cell = -1
      

      contains


      subroutine eval_equ(s, nvar, skip_partials, &
            convert_ODE_to_DAE_form, xscale, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nvar
         logical, intent(in) :: skip_partials, convert_ODE_to_DAE_form
         real(dp), pointer, dimension(:,:) :: xscale ! (nvar, nz)
         integer, intent(out) :: ierr
         call eval_equ_for_solver(s, nvar, 1, s% nz, skip_partials, &
            convert_ODE_to_DAE_form, xscale, ierr)
      end subroutine eval_equ
            
            
      integer function equ_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         type (star_info), pointer :: s
         integer :: ierr
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         equ_extra_profile_columns = s% nvar_hydro
      end function equ_extra_profile_columns
      
      
      subroutine equ_data_for_extra_profile_columns( &
            id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: maxlen_profile_column_name, star_info
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: i, k
         type (star_info), pointer :: s
         real(dp), dimension(:, :), pointer :: equ
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% nvar_hydro /= n) then
            write(*,3) 'nvar_hydro /= n', s% nvar_hydro, n
            stop 'equ_data_for_extra_profile_columns'
         end if
         equ(1:n,1:nz) => s% equ1(1:n*nz)
         do i=1,n
            do k=1,nz
               vals(k,i) = equ(i,k)
            end do
            names(i) = s% nameofequ(i)
         end do
      end subroutine equ_data_for_extra_profile_columns
      

      subroutine eval_equ_for_solver( &
            s, nvar, nzlo, nzhi, skip_partials, &
            convert_ODE_to_DAE_form, xscale, ierr)
         use chem_def
         use mesh_functions
         use hydro_chem_eqns, only: do_chem_eqns, do1_chem_eqns
         use eps_grav, only: prepare_for_eval_eps_grav_and_partials
         use profile, only: do_save_profiles
         use star_utils, only: &
            no_extra_profile_columns, no_data_for_extra_profile_columns

         type (star_info), pointer :: s         
         integer, intent(in) :: nvar, nzlo, nzhi
         logical, intent(in) :: skip_partials, convert_ODE_to_DAE_form
         real(dp), pointer, dimension(:,:) :: xscale ! (nvar, nz)
         integer, intent(out) :: ierr

         integer :: &
            i_dv_dt, i_dlnT_dm, i_dlnd_dt, i_dE_dt, i_dlnR_dt, equchem1, &
            i, k, j, nvar_hydro, nvar_chem, nz, op_err, matrix_type
         integer :: &
            i_lnd, i_lnPgas, i_lnR, i_lnT, i_lum, i_v, i_E, &
            i_chem1, i_xh1, i_xhe4, kmax_equ(nvar), species
         real(dp) :: max_equ(nvar)
         real(dp) :: dVARdot_dVAR, L_phot_old
         real(dp), dimension(:), pointer :: &
            L, lnR, lnP, lnT, energy
         real(dp), dimension(:, :), pointer :: equ
         logical :: v_flag, L_flag, dump_for_debug, do_visc, do_chem, do_mix
         
         include 'formats'
         
         ierr = 0
         
         if (s% do_struct_hydro .and. s% use_dvdt_form_of_momentum_eqn) then         
            if (s% lnPgas_flag) &
               stop 'must not call do1_dvdt_eqn when using lnPgas vars'
            if (.not. s% v_flag) &
               stop 'must not call do1_dvdt_eqn unless using velocity vars'
            if (s% conv_dP_term_factor > 0) & 
               stop 'conv dP not supported in do1_dvdt'
            if (s% use_mass_corrections) &
               stop 'use_mass_corrections dP not supported in do1_dvdt'
            if (s% use_gr_factors) &
               stop 'use_gr_factors dP not supported in do1_dvdt'
            if (s% rotation_flag .and. s% use_dP_dm_rotation_correction) &
               stop 'fp_rot corrections not supported in do1_dvdt'
            if (s% tau(2) < 2d0/3d0) then
               stop 'tau < 2/3 corrections not supported in do1_dvdt'
            end if
         end if

         dump_for_debug = .false. !(s% hydro_call_number == 3)
         !dump_for_debug = .true.         
         
         do_visc = (s% v_flag .and. s% use_artificial_viscosity)
         do_mix = s% do_mix
         do_chem = (do_mix .or. s% do_burn)
         
         call unpack
         
         if (.not. (s% do_struct_hydro .or. s% do_struct_thermo)) then 
         
            if (.not. do_chem) then
               write(*,*) 'bug: do_chem, do_struct_hydro, and do_struct_thermo'
               stop 1
            end if
            ! hold structure constant while solve burn and/or mix
            do j=1,nvar_hydro
               call null_eqn(j)
            end do
            if (ierr == 0) &
               call do_chem_eqns( &
                  s, xscale, nvar, equchem1, nvar_chem, &
                  skip_partials, convert_ODE_to_DAE_form, equ, ierr)
               
         else ! solving structure equations
            
            if (s% use_other_energy_implicit) then
               call s% other_energy_implicit(s% id, ierr)
               if (ierr /= 0) then
                  if (s% report_ierr) &
                     write(*,*) 'ierr from other_energy_implicit'
                  return
               end if
            end if
            
            if (.not. s% use_dedt_form_of_energy_eqn) then
               call prepare_for_eval_eps_grav_and_partials(s,ierr)
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,*) &
                     'ierr from prepare_for_eval_eps_grav_and_partials'
                  return
               end if
            end if

!$OMP PARALLEL DO PRIVATE(op_err,k)
            do k = nzlo, nzhi     
               if (.not. skip_partials .and. &
                     matrix_type == block_tridiag_dble_matrix_type) then
                  s% dblk(:,:,k) = 0
                  s% ublk(:,:,k) = 0
                  s% lblk(:,:,k) = 0
               end if
               op_err = 0
               if (s% do_struct_hydro) then
                  call do1_dlnd_dt_eqn( &
                     s, k, xscale, equ, &
                     skip_partials, convert_ODE_to_DAE_form, &
                     nvar, op_err)
                  if (op_err /= 0) ierr = op_err            
                  if (k > 1) then
                     if (s% use_dvdt_form_of_momentum_eqn) then
                        call do1_dvdt_eqn( &
                           s, k, -1d0, xscale, equ, do_visc, &
                           skip_partials, convert_ODE_to_DAE_form, &
                           nvar, op_err)
                     else
                        call do1_dPdm_eqn( &
                           s, k, -1d0, xscale, equ, skip_partials, &
                           do_visc, nvar, op_err)
                     end if
                  end if
                  if (op_err /= 0) ierr = op_err 
                  if (i_dlnR_dt /= 0) &
                     call do1_dlnR_dt_eqn( &
                        s, k, xscale, equ, &
                        skip_partials, convert_ODE_to_DAE_form, &
                        nvar, op_err)
                  if (op_err /= 0) ierr = op_err 
               end if
               if (s% do_struct_thermo) then
                  if (s% use_dedt_form_of_energy_eqn) then
                     call do1_dedt_eqn( &
                        s, k, xscale, equ, &
                        skip_partials, convert_ODE_to_DAE_form, &
                        do_chem, do_visc, nvar, op_err)
                  else
                     call do1_dLdm_eqn( &
                        s, k, xscale, equ, skip_partials, do_chem, do_visc, &
                        nvar, op_err)
                  end if
                  if (op_err /= 0) ierr = op_err
                  if (L_flag .and. k > 1) then
                     call do1_dlnT_dm_eqn( &
                        s, k, xscale, equ, skip_partials, nvar, op_err)
                     if (op_err /= 0) ierr = op_err
                  end if
               end if                                
               if (do_chem) then
                  call do1_chem_eqns( &
                     s, xscale, k, nvar, equchem1, nvar_chem, &
                     skip_partials, convert_ODE_to_DAE_form, equ, op_err)
                  if (op_err /= 0) ierr = op_err
               end if
            end do      
!$OMP END PARALLEL DO

            if (ierr == 0 .and. nzlo == 1 .and. &
                  (s% do_struct_hydro .or. s% do_struct_thermo)) then
               if (dbg) write(*,*) 'call PT_eqns_surf'
               call PT_eqns_surf( &
                  s, xscale, equ, skip_partials, do_visc, &
                  convert_ODE_to_DAE_form, nvar, ierr)
               if (dbg) write(*,*) 'done PT_eqns_surf'
            end if

            if (.not. s% do_struct_hydro) then
               call dummy_eqn(i_dlnd_dt,i_lnR,nzlo,nzhi) 
               if (s% lnPgas_flag) then
                  call dummy_eqn(i_dv_dt,i_lnPgas,max(2,nzlo),nzhi)
               else
                  call dummy_eqn(i_dv_dt,i_lnd,max(2,nzlo),nzhi)
               end if
               if (i_dlnR_dt /= 0) call dummy_eqn(i_dlnR_dt,i_v,nzlo,nzhi)
            end if
            
            if (L_flag .and. .not. s% do_struct_thermo) then
               call dummy_eqn(i_dE_dt,i_lum,nzlo,nzhi)
               call dummy_eqn(i_dlnT_dm,i_lnT,max(2,nzlo),nzhi)
            end if
            
         end if

         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'ierr after eval equations'
            return
         end if

         if (s% hydro_check_everything .and. &
               s% hydro_call_number == s% hydro_dump_call_number) then
            call check_everything
         end if
         
         if (dump_for_debug .and. .not. skip_partials) call dump_equ_and_partials 
         
         
         return
         
         if (skip_partials) return
         
         
         !do i=1,nvar
         !   do k=1,nz
         !      s% profile_extra(k,i) = s% equ(i,k)
         !   end do
         !end do
         s% profile_extra(1,1:10) = 0
         call do_save_profiles( &
            s, 0, equ_extra_profile_columns, &
            equ_data_for_extra_profile_columns, ierr)
            
         stop 'equations'
         
         
         
         
         contains
         
         
         subroutine unpack
         
            include 'formats'
            
            nz = s% nz
            species = s% species
            nvar_hydro = s% nvar_hydro
            matrix_type = s% hydro_matrix_type
            
            lnT => s% lnT
            lnR => s% lnR
            L => s% L
            lnP => s% lnP
            energy => s% energy
            equ(1:nvar,1:nz) => s% equ1(1:nvar*nz)
            dVARdot_dVAR = s% dVARdot_dVAR
            
            i_dv_dt = s% i_dv_dt
            i_dlnT_dm = s% i_dlnT_dm
            i_dlnd_dt = s% i_dlnd_dt
            i_dE_dt = s% i_dE_dt
            i_dlnR_dt = s% i_dlnR_dt
            
            equchem1 = s% equchem1
            
            i_lnd = s% i_lnd
            i_lnPgas = s% i_lnPgas
            i_lnT = s% i_lnT
            i_lnR = s% i_lnR
            i_lum = s% i_lum
            i_v = s% i_v
            i_E = s% i_E

            i_chem1 = s% i_chem1
            i_xh1 = i_chem1-1+s% net_iso(ih1)
            i_xhe4 = i_chem1-1+s% net_iso(ihe4)
            
            L_phot_old = s% L_phot_old
            v_flag = s% v_flag
            L_flag = s% L_flag
            
            nvar_chem = s% nvar_chem
            
         end subroutine unpack
         
         
         subroutine null_eqn(j)
            integer, intent(in) :: j
            integer :: k
            do k=nzlo,nzhi
               equ(j,k) = 0 ! s% xs(j,k) - s% xs_pre_pass(j,k)
               if (.not. skip_partials) call e00(s, xscale,j,j,k,nvar,one)
            end do
         end subroutine null_eqn
         
         
         subroutine dummy_eqn(j,i,nzlo,nzhi)
            integer, intent(in) :: j, i, nzlo, nzhi
            integer :: k
            do k=nzlo,nzhi
               equ(j,k) = 0
               if (.not. skip_partials) call e00(s,xscale,j,i,k,nvar,one)
            end do
         end subroutine dummy_eqn


         subroutine check(j_equ, k)
            integer, intent(in) :: j_equ, k
            integer :: j
            logical, parameter :: dump_it = .false.
            if (abs(equ(j_equ, k)) > max_equ(j_equ)) then
               max_equ(j_equ) = abs(equ(j_equ, k))
               kmax_equ(j_equ) = k
            end if
            
            8 format(a30, 2i6, 1pd26.16, 3x, a)
            9 format(a30, 3i6, 1pd26.16, 3x, a)
                        
            if (dump_it .or. is_bad_num(equ(j_equ, k))) then
               write(*,8) 'equ(j_equ, k)', j_equ, k, &
                  equ(j_equ, k), trim(s% nameofequ(j_equ))
               write(*,*) 'skip_partials', skip_partials
            end if
         end subroutine check
         
         
         subroutine dump_i_dE_dt_info
            integer :: k, j, k0, k1
            include 'formats'
            write(*,*) 'dump_i_dE_dt_info'
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 's% energy(k)', k, s% energy(k)
                  write(*,2) 's% P(k)', k, s% P(k)
                  write(*,2) 's% rho(k)', k, s% rho(k)
                  write(*,2) 's% T(k)', k, s% T(k)
                  write(*,2) 's% Cv(k)', k, s% Cv(k)
                  write(*,2) 's% dE_dRho(k)', k, s% dE_dRho(k)
                  write(*,2) 's% chiT(k)', k, s% chiT(k)
                  write(*,2) 's% chiRho(k)', k, s% chiRho(k)
                  write(*,2) 's% d_eps_grav_dlnT00(k)', k, s% d_eps_grav_dlnT00(k)
                  write(*,2) 's% d_eps_grav_dlnd00(k)', k, s% d_eps_grav_dlnd00(k)
                  write(*,2) 's% eps_grav(k)', k, s% eps_grav(k)
               end do
               write(*,*)
            end do
            stop 'dump_i_dE_dt_info'      
         end subroutine dump_i_dE_dt_info


         subroutine dump_equ
            integer :: k, j, k0, k1
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 'equ ' // trim(s% nameofequ(j)), k, equ(j, k)
               end do
               write(*,*)
            end do
         end subroutine dump_equ


         subroutine dump_equ_and_partials
            integer :: k, j, i
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,2) 'equ ' // trim(s% nameofequ(j)), k, equ(j, k)
                  do i=1,nvar
                     if (s% dblk(i,j,k) /= 0d0) &
                        write(*,2) 'e00 ' // trim(s% nameofequ(j)) // ' wrt ' // &
                           trim(s% nameofvar(i)), k, s% dblk(i,j,k)
                     if (k < s% nz .and. s% ublk(i,j,k) /= 0d0) &
                        write(*,2) 'ep1 ' // trim(s% nameofequ(j)) // ' wrt ' // &
                           trim(s% nameofvar(i)), k, s% ublk(i,j,k)
                     if (k > 1 .and. s% lblk(i,j,k) /= 0d0) &
                        write(*,2) 'em1 ' // trim(s% nameofequ(j)) // ' wrt ' // &
                           trim(s% nameofvar(i)), k, s% lblk(i,j,k)
                  end do
               end do
               write(*,*)
            end do
            stop 'dump_equ'      
         end subroutine dump_equ_and_partials


         subroutine dump_some_equ
            integer :: k, j, i
            include 'formats'
            do k=536, 538
               do j=1,nvar
                  write(*,2) 'equ ' // trim(s% nameofequ(j)), k, equ(j, k)
               end do
               write(*,*)
            end do
            !stop 'dump_some_equ'      
         end subroutine dump_some_equ


         subroutine check_everything
            use star_utils, only: std_write_internals_to_file
            integer :: k, j, k0, k1
            include 'formats'
            write(*,*) 'check everything', nvar, nz
            max_equ(:) = 0; kmax_equ(:) = 0
            do k=nzlo,nzhi
               do j=1, nvar
                  call check(j, k)
               end do
            end do
            write(*,*)
            do j=1,nvar
               write(*,5) trim(s% nameofequ(j)) // ' kmax var iter model', &
                  kmax_equ(j), j, &
                  s% newton_iter, s% model_number, max_equ(j)
            end do
            write(*,*)
            ! set k0 and k1 to locations where want more info
            k0 = 1; k1 = 0
            if (k0 <= k1) then
               do k=k0,k1
                  do j=1,nvar
                     write(*,4) 'equ ' // trim(s% nameofequ(j)) // ' k iter model', &
                        k, s% newton_iter, s% model_number, equ(j, k)
                  end do
                  write(*,*)
               end do
               write(*,*)
               write(*,*)
               do k=k0,k1
                  write(*,4) '1-q(k)', k, s% newton_iter, s% model_number, 1-s% q(k)
               end do
               write(*,*)
               write(*,*)
               do k=k0,k1
                  write(*,4) 'lnd(k)/ln10', &
                     k, s% newton_iter, s% model_number, s% lnd(k)/ln10
               end do
               write(*,*)
               write(*,*)
               do k=k0,k1
                  write(*,4) 'lnT(k)/ln10', &
                     k, s% newton_iter, s% model_number, s% lnT(k)/ln10
               end do
               write(*,*)
               write(*,*)
               do k=k0,k1
                  do j=1,nvar_chem
                     write(*,4) 'xa_old(j,k) ' // &
                        trim(chem_isos% name(s% chem_id(j))), &
                        j, s% xa_old(j,k)
                  end do
               end do
               write(*,*)
               write(*,*)
               do k=k0,k1
                  do j=1,nvar_chem
                     write(*,3) 'xa(j,k) ' // trim(chem_isos% name(s% chem_id(j))), &
                        j, s% xa(j,k)
                  end do
               end do
            end if
            write(*,*) 'nz', nz
            write(*,*)
            write(*,*) 'finished check everything'
            !stop 'check_everything'
            write(*,*)
            write(*,*)
            write(*,*)
            call std_write_internals_to_file(s% id, s% newton_iter)
            
            if (skip_partials) return
            
         end subroutine check_everything
         

      end subroutine eval_equ_for_solver
      

      subroutine do1_dvdt_eqn( &
            s, k, P_surf, xscale, equ, do_visc, &
            skip_partials, convert_ODE_to_DAE_form, &
            nvar, ierr)
         use hydro_vars, only: do1_dvdt
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: &
            do_visc, skip_partials, convert_ODE_to_DAE_form
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
      
         integer :: j, i_dv_dt, i_v
         logical :: dbg
         real(dp) :: &
            dvdt_expected, scal, dPdm_scal, energy, &
            d_dvdt_dlnd00, d_dvdt_dlndm1, &
            d_dvdt_dlnT00, d_dvdt_dlnTm1, &
            d_dvdt_dlnR00, d_dvdt_dlnRm1, d_dvdt_dlnRp1, &
            d_dvdt_dv00, d_dvdt_dvm1, d_dvdt_dvp1, &
            dlnT_dlnd_const_E, dlnT_dlnE_const_Rho, &
            d_dlnd00_const_E, d_dlnE00_const_Rho, &
            d_dlndm1_const_E, d_dlnEm1_const_Rho

         
         real(dp), parameter :: &   ! DEBUGGING
            lnR00 =  2.6083030127653043D+01, &
              L00 =  3.6037336523777561D+39, &
              v00 = -2.2031511576631115D+00, &
            lnE00 =  3.5914870255149779D+01, &
            lnd00 = -7.8472080075830164D-01, &
            lnRm1 =  2.6083346490237602D+01, &
            lndm1 = -7.8680089580703194D-01, &
            lnEm1 =  3.5913655274076980D+01, &
              Lm1 =  3.6120639103861149D+39, &
              vm1 = -2.2152957245440734D+00, &
            lnRp1 =  2.6082318695086247D+01, &
            lndp1 = -7.8185523563292814D-01, &
            lnEp1 =  3.5916559212315519D+01, &
              Lp1 =  3.5849522633647453D+39, &
              vp1 = -2.1747690066868754D+00, &
              equ_test = -3.6696185099721572D-06


         include 'formats'

         ierr = 0
         i_dv_dt = s% i_dv_dt
         i_v = s% i_v
   
         dbg = .false.
      
         call do1_dvdt( &
            s, k, P_surf, do_visc, skip_partials, dvdt_expected, &
            d_dvdt_dlnd00, d_dvdt_dlndm1, &
            d_dvdt_dlnT00, d_dvdt_dlnTm1, &
            d_dvdt_dlnR00, d_dvdt_dlnRm1, d_dvdt_dlnRp1, &
            d_dvdt_dv00, d_dvdt_dvm1, d_dvdt_dvp1, &
            ierr)
         if (ierr /= 0) return            
         
         s% d_KE_dt_actual(k) = s% dm_bar(k)*s% v(k)*s% dv_dt(k)
         s% d_KE_dt_expected(k) = s% dm_bar(k)*s% v(k)*dvdt_expected
         
         if (convert_ODE_to_DAE_form) then      
            ! make residual units be relative difference in energy
            if (k == 1) then
               energy = s% dm_bar(k)*s% energy_start(k)
            else if (k < s% nz) then
               energy = &
                  0.5d0*(s% dm(k-1)*s% energy_start(k-1) + s% dm(k)*s% energy_start(k))
            else
               energy = &
                  0.5d0*s% dm(k-1)*s% energy_start(k-1) + s% dm(k)*s% energy_start(k)
            end if
            scal = s% dm_bar(k)*max(abs(s% v_start(k)),1d0)/(energy*s% dVARdot_dVAR)
            equ(i_dv_dt, k) = (dvdt_expected - s% dv_dt(k))*scal
         else
            scal = 1d0
            equ(i_dv_dt, k) = dvdt_expected
         end if
         s% v_residual(k) = equ(i_dv_dt, k)
         
         if (s% using_ode_form) s% ode(i_dv_dt,k) = dvdt_expected

         if (skip_partials) return

         call e00(s, xscale, i_dv_dt, i_v, k, nvar, d_dvdt_dv00*scal)
         if (convert_ODE_to_DAE_form) &
            call e00(s, xscale, i_dv_dt, i_v, k, nvar, -s% dVARdot_dVAR*scal)

         if (k > 1) &
            call em1(s, xscale, i_dv_dt, i_v, k, nvar, d_dvdt_dvm1*scal)
         if (k < s% nz) &
            call ep1(s, xscale, i_dv_dt, i_v, k, nvar, d_dvdt_dvp1*scal)
      
         call e00(s, xscale, i_dv_dt, s% i_lnR, k, nvar, d_dvdt_dlnR00*scal)
         if (k > 1) &
            call em1(s, xscale, i_dv_dt, s% i_lnR, k, nvar, d_dvdt_dlnRm1*scal)
         if (k < s% nz) &
            call ep1(s, xscale, i_dv_dt, s% i_lnR, k, nvar, d_dvdt_dlnRp1*scal)

         if (s% E_flag) then
            
            dlnT_dlnd_const_E = s% dlnT_dlnd_c_E(k)
            dlnT_dlnE_const_Rho = s% dlnT_dlnE_c_Rho(k)
         
            d_dlnd00_const_E = d_dvdt_dlnd00 + d_dvdt_dlnT00*dlnT_dlnd_const_E
            d_dlnE00_const_Rho = d_dvdt_dlnT00*dlnT_dlnE_const_Rho
            
            if (k > 1) then
               dlnT_dlnd_const_E = s% dlnT_dlnd_c_E(k-1)
               dlnT_dlnE_const_Rho = s% dlnT_dlnE_c_Rho(k-1)
            else
               dlnT_dlnd_const_E = 0
               dlnT_dlnE_const_Rho = 0
            end if
         
            d_dlndm1_const_E = d_dvdt_dlndm1 + d_dvdt_dlnTm1*dlnT_dlnd_const_E
            d_dlnEm1_const_Rho = d_dvdt_dlnTm1*dlnT_dlnE_const_Rho
         
            call e00(s, xscale, i_dv_dt, s% i_lnd, k, nvar, d_dlnd00_const_E*scal)
            if (k > 1) &
               call em1(s, xscale, i_dv_dt, s% i_lnd, k, nvar, d_dlndm1_const_E*scal)
         
            if (s% do_struct_thermo) then
               call e00(s, xscale, i_dv_dt, s% i_E, k, nvar, &
                  d_dlnE00_const_Rho*scal/s% energy(k))
               if (k > 1) &
                  call em1(s, xscale, i_dv_dt, s% i_E, k, nvar, &
                     d_dlnEm1_const_Rho*scal/s% energy(k-1))
            end if


            
            if (k == -1) then
               write(*,2) 'dequ_dlndm1', k, d_dlndm1_const_E*scal
               write(*,2) 'dequ_dlnd00', k, d_dlnd00_const_E*scal
               write(*,2) 'dequ_dEm1', k, d_dlnEm1_const_Rho*scal/s% energy(k)
               write(*,2) 'dequ_dE00', k, d_dlnE00_const_Rho*scal/s% energy(k-1)
               write(*,2) 'dequ_dlnR00', k, d_dvdt_dlnR00*scal
               write(*,2) 'dequ_dv00', k, d_dvdt_dv00*scal-s% dVARdot_dVAR*scal
               write(*,2) 'dequ_dvm1', k, d_dvdt_dvm1*scal
               write(*,2) 'dequ_dvp1', k, d_dvdt_dvp1*scal
               write(*,2) 'dequ_dlnRm1', k, d_dvdt_dlnRm1*scal
               write(*,2) 'dequ_dlnRp1', k, d_dvdt_dlnRp1*scal
            end if



         
         else
         
            call e00(s, xscale, i_dv_dt, s% i_lnd, k, nvar, d_dvdt_dlnd00*scal)
            if (k > 1) &
               call em1(s, xscale, i_dv_dt, s% i_lnd, k, nvar, d_dvdt_dlndm1*scal)
         
            if (s% do_struct_thermo) then
               call e00(s, xscale, i_dv_dt, s% i_lnT, k, nvar, d_dvdt_dlnT00*scal)
               if (k > 1) &
                  call em1(s, xscale, i_dv_dt, s% i_lnT, k, nvar, d_dvdt_dlnTm1*scal)
            end if
         
         end if
      
      end subroutine do1_dvdt_eqn
      
            
      subroutine do1_dlnd_dt_eqn( &
            s, k, xscale, equ, &
            skip_partials, convert_ODE_to_DAE_form, &
            nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: &
            skip_partials, convert_ODE_to_DAE_form
         integer, intent(out) :: ierr

         real(dp) :: rp13, dm, rho, dlnd_dt, dr3, dequ_ddr3, ddr3_dlnd, &
            ddr3_dlnPgas_const_T, ddr3_dlnT_const_Pgas, res, &
            rho_div_rho_start, rho_dvA_dm, R200, R2p1, v00, vp1
         integer :: nz, i_dlnd_dt, i_lnd
            
         include 'formats'

         ierr = 0
         nz = s% nz
         i_dlnd_dt = s% i_dlnd_dt
         i_lnd = s% i_lnd
         
         dm = s% dm(k)
         rho = s% rho(k)
         
         R200 = s% R2(k)
         v00 = s% vc(k)
         if (k < nz) then
            R2p1 = s% R2(k+1)
            vp1 = s% vc(k+1)
         else
            R2p1 = s% R_center*s% R_center
            vp1 = s% v_center
         end if
         dlnd_dt = -4*pi*rho*(R200*v00 - R2p1*vp1)/dm
         s% lnd_residual(k) = (dlnd_dt - s% dlnd_dt(k))/s% dVARdot_dVAR
            
         if (s% use_dlnddt_form_of_density_eqn) then
         
            stop 'use_dlnddt_form_of_density_eqn is bad for energy conservation residuals'
         
            if (convert_ODE_to_DAE_form) then               
               equ(i_dlnd_dt,k) = (dlnd_dt - s% dlnd_dt(k))
            else
               equ(i_dlnd_dt,k) = dlnd_dt
            end if
         
            if (s% using_ode_form) s% ode(i_dlnd_dt,k) = dlnd_dt
            
            if (skip_partials) return
            
            call e00(s, xscale, i_dlnd_dt, i_lnd, k, nvar, &
               -4*pi*rho*(R200*v00 - R2p1*vp1)/dm)
            if (convert_ODE_to_DAE_form) then
               call e00(s, xscale, i_dlnd_dt, i_lnd, k, nvar, &
                  -s% dVARdot_dVAR)
            end if

            call e00(s, xscale, i_dlnd_dt, s% i_v, k, nvar, &
               -4*pi*rho*R200*s% d_vc_dv/dm)
            call e00(s, xscale, i_dlnd_dt, s% i_lnR, k, nvar, &
               -4*pi*rho*s% d_R2_dlnR(k)*v00/dm)
            
            if (k < nz) then
               call ep1(s, xscale, i_dlnd_dt, s% i_v, k, nvar, &
                  4*pi*rho*R2p1*s% d_vc_dv/dm)
               call ep1(s, xscale, i_dlnd_dt, s% i_lnR, k, nvar, &
                  4*pi*rho*s% d_R2_dlnR(k+1)*vp1/dm)
            end if
            
            return
            
         end if
         
         dr3 = (dm/rho)/(pi4/3)
      
         if (k < nz) then
            rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
         else
            rp13 = s% R_center*s% R_center*s% R_center
         end if
         ! dm = (pi4/3)*(r(k)**3 - rp13)*rho
         ! r(k)**3 = rp13 + (dm/rho)/(pi4/3) = rp13 + dr3      
         
         res = log_cr(rp13 + dr3)
         
         equ(i_dlnd_dt, k) = s% lnR(k) - res/3d0
         s% lnd_residual(k) = equ(i_dlnd_dt, k)

         if (k == s% trace_k) then
            write(*,5) 'i_dlnd_dt', k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, equ(i_dlnd_dt, k)
         end if

         if (s% using_ode_form) s% ode(i_dlnd_dt,k) = 0
      
         if (skip_partials) return

         call e00(s, xscale, i_dlnd_dt, s% i_lnR, k, nvar, one)
         if (k < nz) call ep1( &
            s, xscale, i_dlnd_dt, s% i_lnR, k, nvar, -rp13/(rp13 + dr3))

         dequ_ddr3 = -one_third/(rp13 + dr3)
         ddr3_dlnd = -dr3
         
         if (s% lnPgas_flag) then
            if (s% i_lnT == 0) then
               write(*,*) 'when not using lnT, then must use_dlnddt_form_of_density_eqn'
               stop 1
            end if
            ddr3_dlnPgas_const_T = ddr3_dlnd*s% dlnRho_dlnPgas_const_T(k)
            call e00(s, xscale, i_dlnd_dt, s% i_lnPgas, &
               k, nvar, dequ_ddr3*ddr3_dlnPgas_const_T)
            ddr3_dlnT_const_Pgas = ddr3_dlnd*s% dlnRho_dlnT_const_Pgas(k)
            call e00(s, xscale, i_dlnd_dt, s% i_lnT, &
               k, nvar, dequ_ddr3*ddr3_dlnT_const_Pgas)
         else
            call e00(s, xscale, i_dlnd_dt, s% i_lnd, k, nvar, dequ_ddr3*ddr3_dlnd)
         end if
         
      end subroutine do1_dlnd_dt_eqn
         
         
      subroutine do1_dlnR_dt_eqn( &
            s, k, xscale, equ, &
            skip_partials, convert_ODE_to_DAE_form, &
            nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: &
            skip_partials, convert_ODE_to_DAE_form
         integer, intent(out) :: ierr
         
         real(dp) :: cs, r, vc, scal, PE_scal
         integer :: nz, i_dlnR_dt, i_v, i_lnR
         logical :: force_zero_v
         
         include 'formats'
         ierr = 0
         nz = s% nz
         i_dlnR_dt = s% i_dlnR_dt
         i_v = s% i_v
         i_lnR = s% i_lnR
         r = s% r(k)
         vc = s% vc(k)
         
         PE_scal = s% cgrav(k)*s% m(k)*s% dm_bar(k)/r ! convert dlnR/dt to PE
         s% d_PE_dt_expected(k) = (vc/r)*PE_scal
         s% d_PE_dt_actual(k) = s% dlnR_dt(k)*PE_scal
         
         force_zero_v = (s% q(k) > s% velocity_q_upper_bound)
         if (s% i_lnT /= 0 .and. .not. force_zero_v) &
            force_zero_v = &
               (s% xh_old(s% i_lnT,k)/ln10 < s% velocity_logT_lower_bound .and. &
                  s% dt < secyer*s% max_dt_yrs_for_velocity_logT_lower_bound)
         if (force_zero_v) then
            cs = s% csound_start(k)
            equ(i_dlnR_dt, k) = s% v(k)/cs ! this makes v(k) => 0
            if (skip_partials) return
            call e00(s, xscale, i_dlnR_dt, i_v, k, nvar, one/cs)     
            return
         end if
         
         if (convert_ODE_to_DAE_form) then
            scal = max(1d-10,s% dt)
            equ(i_dlnR_dt, k) = (vc/r - s% dlnR_dt(k))*scal
         else
            scal = 1d0
            equ(i_dlnR_dt, k) = vc/r
         end if
         
         s% lnR_residual(k) = equ(i_dlnR_dt, k)
     
         if (skip_partials) return

         call e00(s, xscale, i_dlnR_dt, i_lnR, k, nvar, -vc/r*scal)
         if (convert_ODE_to_DAE_form) &
            call e00(s, xscale, i_dlnR_dt, i_lnR, k, nvar, -s% dVARdot_dVAR*scal)
      
         call e00(s, xscale, i_dlnR_dt, i_v, k, nvar, s% d_vc_dv/r*scal)
            
      end subroutine do1_dlnR_dt_eqn
      
      
      subroutine do1_dedt_eqn( &
            s, k, xscale, equ, &
            skip_partials, convert_ODE_to_DAE_form, do_chem, do_visc, &
            nvar, ierr)
         use hydro_vars, only: do1_dedt
         use eos_def, only: i_lnE
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: &
            skip_partials, convert_ODE_to_DAE_form, do_chem, do_visc
         integer, intent(out) :: ierr
         
         real(dp) :: &
            dedt, energy, scal, dedt_expected, &
            d_dedt_dlnd, d_dedt_dlnT, &
            d_dedt_dlnR00, d_dedt_dlnRp1, &
            d_dedt_dv00, d_dedt_dvp1, &
            d_dedt_dL00, d_dedt_dLp1, &
            d_dedt_dlnd_c_E, d_dedt_dE_c_Rho
         integer :: j, i_dE_dt
         
         include 'formats'
         
         ierr = 0
         i_dE_dt = s% i_dE_dt
         
         !write(*,2) 'call do1_dedt', k

         call do1_dedt( &
            s, k, do_visc, skip_partials, dedt_expected, &
            d_dedt_dlnd, d_dedt_dlnT, &
            d_dedt_dlnR00, d_dedt_dlnRp1, &
            d_dedt_dv00, d_dedt_dvp1, &
            d_dedt_dL00, d_dedt_dLp1, &
            d_dedt_dlnd_c_E, d_dedt_dE_c_Rho, &
            ierr)
         if (ierr /= 0) return
         
         energy = s% energy(k)
         
         if (s% E_flag) then
            dedt = s% de_dt(k)
         else
            dedt = (s% energy(k) - s% energy_start(k))/s% dt  
            if (s% ebdf_order > 1) then
!$omp critical
               write(*,*) 'cannot use use_dedt_form_of_energy_eqn with higher order'
               write(*,*) 'unless switch to using lnE instead of lnT as basic variable.'
               stop 'do1_dlnEdt_eqn'
!$omp end critical
            end if             
         end if

         s% d_IE_dt_actual(k) = s% dm(k)*dedt
         s% d_IE_dt_expected(k) = s% dm(k)*dedt_expected
         
         if (convert_ODE_to_DAE_form .or. .not. s% E_flag) then
            scal = 1d0/(s% energy_start(k)*s% dVARdot_dVAR) ! ~ dt/energy
            equ(i_dE_dt, k) = (dedt_expected - dedt)*scal
         else
            scal = 1d0
            equ(i_dE_dt, k) = dedt_expected
         end if
         s% E_residual(k) = equ(i_dE_dt, k)

         if (s% using_ode_form) s% ode(i_dE_dt,k) = dedt_expected
         
         if (is_bad_num(equ(i_dE_dt, k))) then
!$omp critical
            write(*,2) 'equ(i_dE_dt, k)', k, equ(i_dE_dt, k)
            write(*,2) 'dedt_expected', k, dedt_expected
            write(*,2) 'dedt', k, dedt
            write(*,2) 'energy', k, energy
            stop 'do1_dedt_eqn'
!$omp end critical
         end if
  
         if (skip_partials) return
         
         if (do_chem .and. & 
               (s% dxdt_nuc_factor > 0d0 .or. s% mix_factor > 0d0)) then
            do j=1,s% nvar_chem
               call e00(s, xscale, i_dE_dt, s% i_chem1+j-1, k, nvar, &
                        s% d_epsnuc_dx(j,k)*scal)
            end do
         end if
         
         if (s% E_flag) then         
            call e00(s, xscale, i_dE_dt, s% i_lnd, k, nvar, d_dedt_dlnd_c_E*scal)              
            call e00(s, xscale, i_dE_dt, s% i_E, k, nvar, d_dedt_dE_c_Rho*scal)
         else         
            call e00(s, xscale, i_dE_dt, s% i_lnd, k, nvar, d_dedt_dlnd*scal)
            call e00(s, xscale, i_dE_dt, s% i_lnT, k, nvar, d_dedt_dlnT*scal)            
         end if
            
         if (convert_ODE_to_DAE_form) then
            if (s% E_flag) then
               call e00(s, xscale, i_dE_dt, s% i_E, k, nvar, -s% dVARdot_dVAR*scal)
            else
               call e00(s, xscale, i_dE_dt, s% i_lnd, k, nvar, &
                  -s% d_eos_dlnd(i_lnE,k)/energy*s% dVARdot_dVAR*scal)
               call e00(s, xscale, i_dE_dt, s% i_lnT, k, nvar, &
                  -s% d_eos_dlnT(i_lnE,k)/energy*s% dVARdot_dVAR*scal)
            end if
         end if
            
         call e00(s, xscale, i_dE_dt, s% i_lnR, k, nvar, d_dedt_dlnR00*scal)
         if (s% i_lum /= 0) &
            call e00(s, xscale, i_dE_dt, s% i_lum, k, nvar, d_dedt_dL00*scal)
         if (s% i_v /= 0) &
            call e00(s, xscale, i_dE_dt, s% i_v, k, nvar, d_dedt_dv00*scal)
         
         if (k < s% nz) then
            call ep1(s, xscale, i_dE_dt, s% i_lnR, k, nvar, d_dedt_dlnRp1*scal)
            if (s% i_lum /= 0) &
               call ep1(s, xscale, i_dE_dt, s% i_lum, k, nvar, d_dedt_dLp1*scal)
            if (s% i_v /= 0) &
               call ep1(s, xscale, i_dE_dt, s% i_v, k, nvar, d_dedt_dvp1*scal)
         end if

      end subroutine do1_dedt_eqn
      
      
      ! just relate L_rad to T gradient.
      ! d_P_rad/dm = -<opacity_at_face>*L_rad/(clight*area^2) -- see, e.g., K&W (5.12)
      ! P_rad = (1/3)*crad*T^4
      ! d_P_rad/dm = (crad/3)*(T(k-1)^4 - T(k)^4)/dm_bar
      ! L_rad = L - L_non_rad, L_non_rad = L_start - L_rad_start
      ! L_rad_start = (-d_P_rad/dm_bar*clight*area^2/<opacity_at_face>)_start
      subroutine do1_alt_dlnT_dm_eqn( &
            s, k, xscale, equ, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, beta, opacity_face, area, area2, T4_m1, T4_00, &
            L, P_rad_start, P_rad_m1, P_rad_00, scale, dm_bar, L_rad, L_non_rad, &
            d_P_rad_expected, d_P_rad_actual, r, dr_dL00, &
            d_kap_dlnT00, d_expected_dlnT00, d_actual_dlnT00, dr_dlnT00, &
            d_kap_dlnTm1, d_expected_dlnTm1, d_actual_dlnTm1, dr_dlnTm1, &
            d_kap_dlnd00, d_expected_dlnd00, dr_dlnd00, &
            d_kap_dlndm1, d_expected_dlndm1, dr_dlndm1, &
            dr_dlnT00_const_Pgas, dr_dlnTm1_const_Pgas, dr_dlnR00, &
            dr_dlnPgas00_const_T, dr_dlnPgasm1_const_T, d_Lrad_dL00, &
            d_Lrad_dlnT00, d_Lrad_dlnTm1, d_Lrad_dlnd00, d_Lrad_dlndm1, &
            d_Lrad_dlnR, gradr2
         integer :: i_dlnT_dm
         
         include 'formats'
         ierr = 0
         i_dlnT_dm = s% i_dlnT_dm

         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         beta = 1d0 - alfa
         !P_rad_start = (crad/3)*s% T_start(k)*s% T_start(k)*s% T_start(k)*s% T_start(k)
         !scale = P_rad_start ! divide by this to scale the equation
         scale = s% energy_start(k)*s% rho_start(k)
         dm_bar = s% dm_bar(k)
         L = s% L(k)

         if (s% q(k) > s% qmin_freeze_non_radiative_luminosity) then
            L_non_rad = s% L_non_rad_start(k)
            L_rad = L - L_non_rad
         else if (s% q(k) <= s% qmax_zero_non_radiative_luminosity) then
            L_non_rad = 0d0
            L_rad = L
         else if (s% use_dPrad_dm_form_of_T_gradient_eqn) then
            if (s% gradT(k) < s% gradr(k)) then
               L_rad = L*s% gradT(k)/s% gradr(k) ! C&G 14.109
            else
               L_rad = L
            end if
            L_non_rad = L - L_rad
         else
            ierr = -1
            return
         end if
         
         ! calculate expected d_P_rad from current L_rad
         opacity_face = alfa*s% opacity(k) + beta*s% opacity(k-1)         
         area = s% area(k)
         area2 = area*area
         d_P_rad_expected = -dm_bar*opacity_face*L_rad/(clight*area2)
         
         ! calculate actual d_P_rad in current model
         T4_m1 = s% T(k-1)*s% T(k-1)*s% T(k-1)*s% T(k-1)
         T4_00 = s% T(k)*s% T(k)*s% T(k)*s% T(k)
         P_rad_m1 = (crad/3)*T4_m1
         P_rad_00 = (crad/3)*T4_00
         d_P_rad_actual = P_rad_m1 - P_rad_00
         
         r = (d_P_rad_expected - d_P_rad_actual)/scale
         equ(i_dlnT_dm, k) = r
         if (k == -1) then
            write(*,2) 'lnE/ln10', k, s% lnE(k)/ln10
            write(*,2) 'lnd/ln10', k, s% lnd(k)/ln10
            write(*,2) 'lnT/ln10', k, s% lnT(k)/ln10
            write(*,2) 'opacity_face', k, opacity_face
            write(*,2) 'L_rad', k, L_rad
            write(*,2) 'area', k, area
            write(*,2) 'd_P_rad_expected', k, d_P_rad_expected
            write(*,2) 'd_P_rad_actual', k, d_P_rad_actual
            write(*,2) 'scale', k, scale
            write(*,2) 'r', k, r
            write(*,*) 's% gradT(k) < s% gradr(k)', s% gradT(k) < s% gradr(k)
            write(*,*)
         end if

         s% dlnTdm_residual(k) = equ(i_dlnT_dm,k)
         if (s% using_ode_form) s% ode(i_dlnT_dm,k) = 0
      
         if (skip_partials) return
         
         if (s% use_dPrad_dm_form_of_T_gradient_eqn .and. &
               s% gradT(k) < s% gradr(k)) then ! L_rad = L*s% gradT(k)/s% gradr(k)
            gradr2 = s% gradr(k)*s% gradr(k)
            d_Lrad_dL00 = s% gradT(k)/s% gradr(k) + L*( &
               s% d_gradT_dL(k)/s% gradr(k) &
               - s% gradT(k)*s% d_gradr_dL(k)/gradr2)
            d_Lrad_dlnR = L*( &
               s% d_gradT_dlnR(k)/s% gradr(k) &
               - s% gradT(k)*s% d_gradr_dlnR(k)/gradr2)
            d_Lrad_dlnT00 = L*( &
               s% d_gradT_dlnT00(k)/s% gradr(k) &
               - s% gradT(k)*s% d_gradr_dlnT00(k)/gradr2)
            d_Lrad_dlnTm1 = L*( &
               s% d_gradT_dlnTm1(k)/s% gradr(k) &
               - s% gradT(k)*s% d_gradr_dlnTm1(k)/gradr2)
            d_Lrad_dlnd00 = L*( &
               s% d_gradT_dlnd00(k)/s% gradr(k) &
               - s% gradT(k)*s% d_gradr_dlnd00(k)/gradr2)
            d_Lrad_dlndm1 = L*( &
               s% d_gradT_dlndm1(k)/s% gradr(k) &
               - s% gradT(k)*s% d_gradr_dlndm1(k)/gradr2)               
         else ! Lrad = L
            d_Lrad_dL00 = 1d0
            d_Lrad_dlnR = 0d0
            d_Lrad_dlnT00 = 0d0
            d_Lrad_dlnTm1 = 0d0
            d_Lrad_dlnd00 = 0d0
            d_Lrad_dlndm1 = 0d0
         end if
         
         dr_dL00 = -dm_bar*opacity_face*d_Lrad_dL00/(clight*area2)/scale
         call e00(s, xscale, i_dlnT_dm, s% i_lum, k, nvar, dr_dL00) ! d/d_L00
         
         d_kap_dlnT00 = alfa*s% d_opacity_dlnT(k)            
         d_expected_dlnT00 = -dm_bar/(clight*area2)*( &
            d_kap_dlnT00*L_rad + opacity_face*d_Lrad_dlnT00)
         d_actual_dlnT00 = -4*P_rad_00
         
         dr_dlnT00 = (d_expected_dlnT00 - d_actual_dlnT00)/scale
         
         if (k == -1) then
            write(*,2) 'dr_dlnT00', k, dr_dlnT00
            write(*,2) 'd_expected_dlnT00', k, d_expected_dlnT00
            write(*,2) 'd_actual_dlnT00', k, d_actual_dlnT00
            write(*,2) 'scale', k, scale
            write(*,2) 'P_rad_00', k, P_rad_00
            write(*,2) 'd_kap_dlnT00', k, d_kap_dlnT00
            write(*,2) 'd_Lrad_dlnT00', k, d_Lrad_dlnT00
            write(*,2) 'opacity_face', k, opacity_face
            write(*,2) 'dm_bar', k, dm_bar
            write(*,2) 's% T_start(k)', k, s% T_start(k)
            write(*,*)
         end if
         
         d_kap_dlnTm1 = beta*s% d_opacity_dlnT(k-1)
         d_expected_dlnTm1 = -dm_bar/(clight*area2)*( &
            d_kap_dlnTm1*L_rad + opacity_face*d_Lrad_dlnTm1)
         d_actual_dlnTm1 = 4*P_rad_m1
         
         dr_dlnTm1 = (d_expected_dlnTm1 - d_actual_dlnTm1)/scale
         
         if (s% do_struct_hydro) then ! partials wrt lnd 
                   
            d_kap_dlnd00 = alfa*s% d_opacity_dlnd(k)
            d_expected_dlnd00 = -dm_bar/(clight*area2)*( &
               d_kap_dlnd00*L_rad + opacity_face*d_Lrad_dlnd00)
            dr_dlnd00 = d_expected_dlnd00/scale   
                        
            d_kap_dlndm1 = beta*s% d_opacity_dlnd(k-1)
            d_expected_dlndm1 = -dm_bar/(clight*area2)*( &
               d_kap_dlndm1*L_rad + opacity_face*d_Lrad_dlndm1)
            dr_dlndm1 = d_expected_dlndm1/scale    
                       
         else
            dr_dlnd00 = 0
            dr_dlndm1 = 0
         end if
                     
         if (s% lnPgas_flag) then ! d/dlnT00_const_Pgas, d/dlnTm1_const_Pgas
            dr_dlnT00_const_Pgas = &
               dr_dlnT00 + dr_dlnd00*s% dlnRho_dlnT_const_Pgas(k)
            call e00(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnT00_const_Pgas)
            dr_dlnTm1_const_Pgas = &
               dr_dlnTm1 + dr_dlndm1*s% dlnRho_dlnT_const_Pgas(k-1)
            call em1(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnTm1_const_Pgas)
         else if (s% E_flag) then ! d_dlnE_const_Rho
            call e00(s, xscale, i_dlnT_dm, s% i_E, &
               k, nvar, dr_dlnT00*s% dlnT_dlnE_c_Rho(k)/s% energy(k))
            call em1(s, xscale, i_dlnT_dm, s% i_E, &
               k, nvar, dr_dlnTm1*s% dlnT_dlnE_c_Rho(k-1)/s% energy(k-1))
            if (k == -1) then
               write(*,2) 'd_dlnE_const_Rho', k, dr_dlnT00*s% dlnT_dlnE_c_Rho(k)
               write(*,2) 'dr_dlnT00', k, dr_dlnT00
               write(*,2) 's% dlnT_dlnE_c_Rho(k)', k, s% dlnT_dlnE_c_Rho(k)
            end if
         else ! d/dlnT00, d/dlnTm1
            call e00(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnT00)
            call em1(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnTm1)
         end if

         if (s% do_struct_hydro) then ! partials wrt lnR and (lnd or lnPgas)
            dr_dlnR00 = &
               -4*d_P_rad_expected/scale &
               -dm_bar*opacity_face*d_Lrad_dlnR/(clight*area2)/scale
            call e00(s, xscale, i_dlnT_dm, s% i_lnR, k, nvar, dr_dlnR00) ! d/dlnR
            if (s% lnPgas_flag) then ! d/d_lnPgas00, d/d_lnPgasm1
               dr_dlnPgas00_const_T = dr_dlnd00*s% dlnRho_dlnPgas_const_T(k)
               call e00(s, xscale, i_dlnT_dm, s% i_lnPgas, &
                  k, nvar, dr_dlnPgas00_const_T)
               dr_dlnPgasm1_const_T = dr_dlndm1*s% dlnRho_dlnPgas_const_T(k-1)
               call em1(s, xscale, i_dlnT_dm, s% i_lnPgas, &
                  k, nvar, dr_dlnPgasm1_const_T)
            else if (s% E_flag) then ! d_lnd_const_E
               call e00(s, xscale, i_dlnT_dm, s% i_lnd, k, nvar, &
                  dr_dlnd00 + dr_dlnT00*s% dlnT_dlnd_c_E(k))
               call em1(s, xscale, i_dlnT_dm, s% i_lnd, k, nvar, &
                  dr_dlndm1 + dr_dlnTm1*s% dlnT_dlnd_c_E(k-1))
            else ! d/d_lnd00, d/d_lndm1
               call e00(s, xscale, i_dlnT_dm, s% i_lnd, k, nvar, dr_dlnd00)
               call em1(s, xscale, i_dlnT_dm, s% i_lnd, k, nvar, dr_dlndm1)
            end if
         end if
         
      end subroutine do1_alt_dlnT_dm_eqn


      subroutine do1_dlnT_dm_eqn( &
            s, k, xscale, equ, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
                     
         real(dp) :: alfa, beta, r, dr_dL00, dr_dlnR00, dr_dlnd00, dr_dlnT00, &
            dr_dlndm1, dr_dlnTm1, dr_lnq00, dr_lndq00, dr_lndqm1, m, dlnPdm, &
            d_dlnPdm_dlnR, d_dlnPdm_dL, d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            dP_dlnPgas00_const_T, dP_dlnPgasm1_const_T, &
            dP_dlnT00_const_Pgas, dP_dlnTm1_const_Pgas, &
            grad_star, d_grad_star_dL, d_grad_star_dlnR, &       
            d_grad_star_dlnd00, d_grad_star_dlndm1, &    
            d_grad_star_dlnT00, d_grad_star_dlnTm1, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            d_dlnTdm_dLum, d_dlnTdm_dlnR, d_dlnTdm_dlnd00, &
            d_dlnTdm_dlnT00, d_dlnTdm_dlndm1, d_dlnTdm_dlnTm1, &
            delm, d_delm_dqm1, d_delm_dq00, T00, Tm1, dT, Tpoint, &
            lnTdiff, d_lnTdiff_dlnT00, d_lnTdiff_dlnTm1, d_dlnTdm_dlnq, &
            d_delm_dlndq00, d_delm_dlndqm1, diff, diff_old, &
            d_grad_star_dlnT00_const_Pgas, d_grad_star_dlnTm1_const_Pgas, &
            d_grad_star_dlnPgas00_const_T, d_grad_star_dlnPgasm1_const_T, &
            d_dlnTdm_dlnPgas00_const_T, d_dlnTdm_dlnPgasm1_const_T, &
            d_dlnTdm_dlnT00_const_Pgas, d_dlnTdm_dlnTm1_const_Pgas, &
            dr_dlnPgas00_const_T, dr_dlnPgasm1_const_T, &
            dr_dlnT00_const_Pgas, dr_dlnTm1_const_Pgas, &
            other, other_equ, FL, other_FL
         integer :: i_dlnT_dm
         
         include 'formats'
         ierr = 0
         i_dlnT_dm = s% i_dlnT_dm
         
         if (.not. s% L_flag) then
            write(*,*) 'error to call T gradient eqn when not using L vars'
            stop 1
         end if

         if (s% q(k) > s% qmin_freeze_non_radiative_luminosity .or. &
             s% q(k) <= s% qmax_zero_non_radiative_luminosity .or. &
             s% use_dPrad_dm_form_of_T_gradient_eqn) then
            call do1_alt_dlnT_dm_eqn(s, k, xscale, equ, skip_partials, nvar, ierr)
            return
         end if
         
         if (s% E_flag) then
            write(*,*) 'must use_dPrad_dm_form_of_T_gradient_eqn when have E_flag'
            stop 1
         end if
         
         ! dT/dm = dP/dm * T/P * grad_T, grad_T = dlnT/dlnP from MLT.
         ! but use hydrostatic value for dP/dm in this.
         ! this is because of limitations of MLT for calculating grad_T.
         ! (MLT assumes hydrostatic equilibrium)
         ! see comment in K&W chpt 9.1.

         call eval_dlnPdm_qhse(s, k, m, &
            dlnPdm, d_dlnPdm_dlnR, d_dlnPdm_dL, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            ierr)
         if (ierr /= 0) return
      
         call eval_grad_star_info( &
            s, k, grad_star, d_grad_star_dL, d_grad_star_dlnR, &
            d_grad_star_dlnd00, d_grad_star_dlndm1, &
            d_grad_star_dlnT00, d_grad_star_dlnTm1, &
            d_grad_star_dlnT00_const_Pgas, d_grad_star_dlnTm1_const_Pgas, &
            d_grad_star_dlnPgas00_const_T, d_grad_star_dlnPgasm1_const_T, &
            ierr)
         if (ierr /= 0) return
                  
         s% dlnT_dm_expected(k) = dlnPdm*grad_star 
      
         delm = (s% dm(k) + s% dm(k-1))/2
         Tm1 = s% T(k-1)
         alfa = s% dm(k-1)/(s% dm(k-1) + s% dm(k))
         beta = 1 - alfa

         T00 = s% T(k)
         dT = Tm1 - T00
         Tpoint = alfa*T00 + beta*Tm1
         lnTdiff = dT/Tpoint ! use this in place of lnT(k-1)-lnT(k)
      
         r = delm*s% dlnT_dm_expected(k) - lnTdiff
         equ(i_dlnT_dm, k) = r

         if (k == s% trace_k) then
            write(*,5) 'i_dlnT_dm', k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, equ(i_dlnT_dm, k)
         end if
         
         if (is_bad_num(equ(i_dlnT_dm, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_dlnT_dm, k)', k, equ(i_dlnT_dm, k)
            return
            write(*,2) 'equ(i_dlnT_dm, k)', k, equ(i_dlnT_dm, k)
            write(*,2) 'lnTdiff', k, lnTdiff
            write(*,2) 'delm', k, delm
            write(*,2) 'dlnT_dm_expected(k)', k, s% dlnT_dm_expected(k)
            write(*,2) 'dlnPdm', k, dlnPdm
            write(*,2) 'grad_star', k, grad_star
            stop 'i_dlnT_dm'
         end if
         
         other = 8.3276707878306060D+40
         other_equ = 1.4606961473706548D-12
         other_FL = 9.4222987519533120D+01
         
         s% dlnTdm_residual(k) = equ(i_dlnT_dm,k)
         if (s% using_ode_form) s% ode(i_dlnT_dm,k) = 0
      
         if (skip_partials) return
      
         d_lnTdiff_dlnTm1 = T00*Tm1/(Tpoint*Tpoint)
         d_lnTdiff_dlnT00 = -d_lnTdiff_dlnTm1

         d_dlnTdm_dLum = dlnPdm*d_grad_star_dL + d_dlnPdm_dL*grad_star
         dr_dL00 = delm*d_dlnTdm_dLum

         d_dlnTdm_dlnR = dlnPdm*d_grad_star_dlnR + d_dlnPdm_dlnR*grad_star
         dr_dlnR00 = delm*d_dlnTdm_dlnR
         
         if (s% i_lum /= 0) &
            call e00(s, xscale, i_dlnT_dm, s% i_lum, k, nvar, dr_dL00) ! d/d_L00
         
         if (s% lnPgas_flag) then ! d/dlnT00_const_Pgas, d/dlnTm1_const_Pgas
            d_dlnTdm_dlnT00_const_Pgas = &
               dlnPdm*d_grad_star_dlnT00_const_Pgas + &
                  d_dlnPdm_dlnT00_const_Pgas*grad_star
            dr_dlnT00_const_Pgas = delm*d_dlnTdm_dlnT00_const_Pgas - d_lnTdiff_dlnT00
            call e00(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnT00_const_Pgas)
            d_dlnTdm_dlnTm1_const_Pgas =  &
               dlnPdm*d_grad_star_dlnTm1_const_Pgas + &
                  d_dlnPdm_dlnTm1_const_Pgas*grad_star
            dr_dlnTm1_const_Pgas = delm*d_dlnTdm_dlnTm1_const_Pgas - d_lnTdiff_dlnTm1
            call em1(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnTm1_const_Pgas)
         else ! d/dlnT00, d/dlnTm1
            d_dlnTdm_dlnT00 = dlnPdm*d_grad_star_dlnT00 + d_dlnPdm_dlnT00*grad_star
            dr_dlnT00 = delm*d_dlnTdm_dlnT00 - d_lnTdiff_dlnT00
            call e00(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnT00)
            d_dlnTdm_dlnTm1 = dlnPdm*d_grad_star_dlnTm1 + d_dlnPdm_dlnTm1*grad_star
            dr_dlnTm1 = delm*d_dlnTdm_dlnTm1 - d_lnTdiff_dlnTm1
            call em1(s, xscale, i_dlnT_dm, s% i_lnT, k, nvar, dr_dlnTm1)
         end if

         if (s% do_struct_hydro) then
            call e00(s, xscale, i_dlnT_dm, s% i_lnR, k, nvar, dr_dlnR00) ! d/dlnR
            if (s% lnPgas_flag) then ! d/d_lnPgas00, d/d_lnPgasm1
               d_dlnTdm_dlnPgas00_const_T = &
                  dlnPdm*d_grad_star_dlnPgas00_const_T + &
                     d_dlnPdm_dlnPgas00_const_T*grad_star
               dr_dlnPgas00_const_T = delm*d_dlnTdm_dlnPgas00_const_T
               call e00(s, xscale, i_dlnT_dm, s% i_lnPgas, k, nvar, dr_dlnPgas00_const_T)
               d_dlnTdm_dlnPgasm1_const_T = &
                  dlnPdm*d_grad_star_dlnPgasm1_const_T + &
                     d_dlnPdm_dlnPgasm1_const_T*grad_star
               dr_dlnPgasm1_const_T = delm*d_dlnTdm_dlnPgasm1_const_T
               call em1(s, xscale, i_dlnT_dm, s% i_lnPgas, k, nvar, dr_dlnPgasm1_const_T)
            else ! d/d_lnd00, d/d_lndm1
               d_dlnTdm_dlnd00 = dlnPdm*d_grad_star_dlnd00 + d_dlnPdm_dlnd00*grad_star
               dr_dlnd00 = delm*d_dlnTdm_dlnd00
               call e00(s, xscale, i_dlnT_dm, s% i_lnd, k, nvar, dr_dlnd00)
               d_dlnTdm_dlndm1 = dlnPdm*d_grad_star_dlndm1 + d_dlnPdm_dlndm1*grad_star
               dr_dlndm1 = delm*d_dlnTdm_dlndm1
               call em1(s, xscale, i_dlnT_dm, s% i_lnd, k, nvar, dr_dlndm1)
            end if
         end if

      end subroutine do1_dlnT_dm_eqn
         
         
      subroutine do1_dLdm_eqn( & ! energy conservation
            s, k, xscale, equ, skip_partials, do_chem, do_visc, nvar, ierr)
         use rates_def
         use eps_grav, only: eval_eps_grav_and_partials
         use chem_def, only: iprot
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, do_chem, do_visc
         integer, intent(out) :: ierr
         
         real(dp) :: L00, Lp1, Lx, dA_dlnR, diff, diff_old, &
            dm_max, dm, L_scale, dLdm, partial, &
            d_dLdM_dlndm1, d_dLdM_dlnTm1, &
            d_dLdM_dlnd00, d_dLdM_dlnT00, &
            d_dLdM_dlndp1, d_dLdM_dlnTp1, &
            d_dLdM_dlnR00, d_dLdM_dlnRp1, &
            d_dLdM_dL00, d_dLdM_dLp1, &
            d_faclnd_dlnd, d_faclnT_dlnT, &
            eps_burn, d_eps_burn_dlnd, d_eps_burn_dlnT, &
            d_epsnuc_dlnd, d_epsnuc_dlnT, &
            non_nuc_neu, d_nonnucneu_dlnd, d_nonnucneu_dlnT, faclnd, faclnT, &
            Rho_new, Rho_start, T_new, T_start, dln_dlnd, dln_dlnT, &
            dL_expected, dL_actual, eps_nuc, Flx, dLx, d_dLx_epsnuc, &
            d_del_dlnR00, d_del_dlnRp1, d_del_dvel00, d_del_dvelp1, &
            d_dLAV_dlnd, d_dLAV_dlnR00, d_dLAV_dlnRp1, &
            d_dLAV_dvel00, d_dLAV_dvelp1, dLAV, del, qmid, eps_extra, &
            d_eps_burn_dlnPgas_const_T, d_eps_burn_dlnT_const_Pgas, &
            d_dLdM_dlnPgasm1_const_T, d_dLdM_dlnPgas00_const_T, &
            d_dLdM_dlnPgasp1_const_T, &
            d_dLdM_dlnTm1_const_Pgas, d_dLdM_dlnT00_const_Pgas, &
            d_dLdM_dlnTp1_const_Pgas, &
            d_dLAV_dlnPgas_const_T, &
            d_dLdM_dv00, d_dLdM_dvp1, &
            dLx_dLp1, dLx_dlnT00_const_Pgas, dLx_dlnTm1_const_Pgas, &
            dLx_dlnTp1_const_Pgas, &
            dLx_dlnT00, dLx_dlnTm1, dLx_dlnTp1, dLx_dlnd00, dLx_dx, &
            dLx_dlnPgas00_const_T, dLx_dlnR00, dLx_dv00, dLx_dL00, &
            dLx_dlnRp1, dLx_dvp1, dLx_dlnPgasp1_const_T, &
            dLx_dlndp1, dLx_dlnPgasm1_const_T, dLx_dlndm1, &
            other_equ, other_f, visc_factor
            
         integer :: j, i_dE_dt, i_lum
         logical :: dbg
            
         include 'formats'
         
         i_dE_dt = s% i_dE_dt
         i_lum = s% i_lum
         
         ierr = 0
         call eval_eps_grav_and_partials(s, k, ierr) ! get eps_grav info
         if (ierr /= 0) return

         dm = s% dm(k)
         
         eps_nuc = s% eps_nuc(k)
         d_epsnuc_dlnd = s% d_epsnuc_dlnd(k)
         d_epsnuc_dlnT = s% d_epsnuc_dlnT(k)

         non_nuc_neu = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))
         d_nonnucneu_dlnd = 0.5d0*s% d_nonnucneu_dlnd(k)
         d_nonnucneu_dlnT = 0.5d0*s% d_nonnucneu_dlnT(k)
         
         eps_burn = eps_nuc - non_nuc_neu + s% extra_heat(k) + s% irradiation_heat(k)
         dLdm = s% eps_grav(k) + eps_burn
         
         visc_factor = s% eps_visc_factor
         if (do_visc) dLdm = dLdm + s% eps_visc(k)*visc_factor

         s% dL_dm_expected(k) = dLdm
         dL_expected = dm*dLdm
         
         if (i_lum /= 0) then
            L00 = s% L(k)
            if (k < s% nz) then
               Lp1 = s% L(k+1)
            else
               Lp1 = s% L_center
            end if
         else
            L00 = 0
            Lp1 = 0
         end if
         dL_actual = L00 - Lp1
         
         dLx = dLdm*dm ! dL expected
         Lx = Lp1 + dLx ! Lx = L00 expected = Lp1 + dLdm*dm
         
         if (s% L_flag) then
            L_scale = max(1d-3*Lsun, abs(s% xh_pre(i_lum,k)))
         else
            L_scale = Lsun
         end if
         !L_scale = s% energy_start(k)*dm*s% dVARdot_dVAR
            ! this form of L_scale looks great, but it sometimes causes trouble.
            ! e.g., test cases irradiated_planet and 20M_si_burn_qp_bcyclic fail.
         
         ! energy conservation (luminosity equation)
         equ(i_dE_dt, k) = (Lx - L00)/L_scale
         
         dbg = .false.
         !dbg = (k == 1920) !.and. (.not. skip_partials)
   
         if (is_bad_num(equ(i_dE_dt, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_dE_dt, k)', k, equ(i_dE_dt, k)
            return
         end if
         
         if (k >= s% min_k_for_max_dE_residual) then
            s% E_residual(k) = &
               abs(s% dt/s% energy(k)* &
                  (s% dL_dm_expected(k) - dL_actual/s% dm(k)))
         else
            s% E_residual(k) = 0
         end if

         if (s% using_ode_form) s% ode(i_dE_dt,k) = 0

         if (skip_partials) return
         
         
         if (.false. .and. k == s% nz) then
         
            !write(*,2) 's% dL_dm_expected(k)', k, s% dL_dm_expected(k)
            !write(*,2) 'dL_dm_actual', k, dL_actual/s% dm(k)
            write(*,2) 'dL_actual', k, dL_actual
            write(*,2) 'dL expected', k, s% dL_dm_expected(k)*s% dm(k)
            write(*,2) 's% eps_grav(k)*s% dm(k)', k, s% eps_grav(k) *s% dm(k)
            write(*,2) 's% eps_nuc(k)*s% dm(k)', k, s% eps_nuc(k)*s% dm(k)
            write(*,2) 'non_nuc_neu*s% dm(k)', k, non_nuc_neu*s% dm(k)
            !write(*,2) 's% extra_heat(k)', k, s% extra_heat(k)
            !write(*,2) 's% irradiation_heat(k)', k, s% irradiation_heat(k)
            !if (do_visc) write(*,2) 's% eps_visc(k)', k, s% eps_visc(k)
            write(*,*)
         
         end if
         
         
         
         
         
         
      
         d_dLdM_dlndm1 = 0
         d_dLdM_dlnd00 = 0
         d_dLdM_dlndp1 = 0
         d_dLdM_dlnTm1 = 0
         d_dLdM_dlnT00 = 0
         d_dLdM_dlnTp1 = 0
         
         d_eps_burn_dlnd = d_epsnuc_dlnd - d_nonnucneu_dlnd + s% d_extra_heat_dlnd(k)
         d_eps_burn_dlnT = d_epsnuc_dlnT - d_nonnucneu_dlnT + s% d_extra_heat_dlnT(k)
   
         d_dLdM_dlnR00 = s% d_eps_grav_dlnR00(k)
         d_dLdM_dlnRp1 = s% d_eps_grav_dlnRp1(k)
         
         d_dLdM_dL00 = s% d_eps_grav_dL00(k)
         d_dLdM_dLp1 = s% d_eps_grav_dLp1(k)

         d_dLdM_dv00 = s% d_eps_grav_dv00(k)
         d_dLdM_dvp1 = s% d_eps_grav_dvp1(k)
         
         if (s% lnPgas_flag) then
            
            d_eps_burn_dlnPgas_const_T = d_eps_burn_dlnd*s% dlnRho_dlnPgas_const_T(k)
            d_eps_burn_dlnT_const_Pgas = &
               d_eps_burn_dlnT + d_eps_burn_dlnd*s% dlnRho_dlnT_const_Pgas(k)
   
            d_dLdM_dlnPgasm1_const_T = s% d_eps_grav_dlnPgasm1_const_T(k)
            d_dLdM_dlnPgas00_const_T = &
               s% d_eps_grav_dlnPgas00_const_T(k) + d_eps_burn_dlnPgas_const_T
            d_dLdM_dlnPgasp1_const_T = s% d_eps_grav_dlnPgasp1_const_T(k)
         
            d_dLdM_dlnTm1_const_Pgas = s% d_eps_grav_dlnTm1_const_Pgas(k)
            d_dLdM_dlnT00_const_Pgas = &
               s% d_eps_grav_dlnT00_const_Pgas(k) + d_eps_burn_dlnT_const_Pgas
            d_dLdM_dlnTp1_const_Pgas = s% d_eps_grav_dlnTp1_const_Pgas(k)

         else
   
            d_dLdM_dlndm1 = s% d_eps_grav_dlndm1(k)
            d_dLdM_dlnd00 = s% d_eps_grav_dlnd00(k) + d_eps_burn_dlnd
            d_dLdM_dlndp1 = s% d_eps_grav_dlndp1(k)
         
            d_dLdM_dlnTm1 = s% d_eps_grav_dlnTm1(k)
            d_dLdM_dlnT00 = s% d_eps_grav_dlnT00(k) + d_eps_burn_dlnT
            d_dLdM_dlnTp1 = s% d_eps_grav_dlnTp1(k)
            
            d_dLdM_dlnPgasm1_const_T = 0
            d_dLdM_dlnPgas00_const_T = 0
            d_dLdM_dlnPgasp1_const_T = 0
         
            d_dLdM_dlnTm1_const_Pgas = 0
            d_dLdM_dlnT00_const_Pgas = 0
            d_dLdM_dlnTp1_const_Pgas = 0

         end if

         dLx_dlnT00 = d_dLdm_dlnT00*dm
         dLx_dlnTm1 = d_dLdM_dlnTm1*dm
         dLx_dlnTp1 = d_dLdM_dlnTp1*dm
         dLx_dlnd00 = d_dLdm_dlnd00*dm
         dLx_dlnR00 = d_dLdM_dlnR00*dm
         dLx_dv00 = d_dLdM_dv00*dm
         dLx_dL00 = d_dLdM_dL00*dm
         dLx_dLp1 = d_dLdM_dLp1*dm
         dLx_dlnRp1 = d_dLdM_dlnRp1*dm
         dLx_dvp1 = d_dLdM_dvp1*dm
         dLx_dlndp1 = d_dLdM_dlndp1*dm
         dLx_dlndm1 = d_dLdM_dlndm1*dm
         
         if (s% lnPgas_flag) then
            dLx_dlnT00_const_Pgas = d_dLdm_dlnT00_const_Pgas*dm
            dLx_dlnTm1_const_Pgas = d_dLdM_dlnTm1_const_Pgas*dm
            dLx_dlnTp1_const_Pgas = d_dLdM_dlnTp1_const_Pgas*dm
            dLx_dlnPgasp1_const_T = d_dLdM_dlnPgasp1_const_T*dm
            dLx_dlnPgas00_const_T = d_dLdm_dlnPgas00_const_T*dm
            dLx_dlnPgasm1_const_T = d_dLdM_dlnPgasm1_const_T*dm
         else
            dLx_dlnT00_const_Pgas = 0
            dLx_dlnTm1_const_Pgas = 0
            dLx_dlnTp1_const_Pgas = 0
            dLx_dlnPgasp1_const_T = 0
            dLx_dlnPgas00_const_T = 0
            dLx_dlnPgasm1_const_T = 0
         end if
         
         if (i_lum /= 0) then
            call e00(s, xscale, i_dE_dt, i_lum, k, nvar, dLx_dL00/L_scale - 1/L_scale)
            if (k < s% nz) call ep1( &
               s, xscale, i_dE_dt, i_lum, k, nvar, dLx_dLp1/L_scale + 1/L_scale)
         end if
         
         if (s% lnPgas_flag) then
            call e00(s, xscale, i_dE_dt, s% i_lnT, &
               k, nvar, dLx_dlnT00_const_Pgas/L_scale)
            if (k > 1) then
               call em1(s, xscale, i_dE_dt, s% i_lnT, &
               k, nvar, dLx_dlnTm1_const_Pgas/L_scale)
            end if
            if (k < s% nz) then
               call ep1(s, xscale, i_dE_dt, s% i_lnT, &
               k, nvar, dLx_dlnTp1_const_Pgas/L_scale)
            end if
         else
            call e00(s, xscale, i_dE_dt, s% i_lnT, k, nvar, dLx_dlnT00/L_scale)
            if (k > 1) call em1(s, xscale, i_dE_dt, s% i_lnT, &
               k, nvar, dLx_dlnTm1/L_scale)
            if (k < s% nz) call ep1(s, xscale, i_dE_dt, s% i_lnT, &
               k, nvar, dLx_dlnTp1/L_scale)
         end if

         if (s% do_struct_hydro) then
            if (s% lnPgas_flag) then
               call e00(&
                  s, xscale, i_dE_dt, s% i_lnPgas, &
               k, nvar, dLx_dlnPgas00_const_T/L_scale)
            else
               call e00(&
                  s, xscale, i_dE_dt, s% i_lnd, k, nvar, dLx_dlnd00/L_scale)
            end if
            
            call e00(s, xscale, i_dE_dt, s% i_lnR, k, nvar, dLx_dlnR00/L_scale)
            
            if (s% i_v /= 0) call e00(s, xscale, i_dE_dt, s% i_v, &
               k, nvar, dLx_dv00/L_scale)
            if (k < s% nz) then
               call ep1(s, xscale, i_dE_dt, s% i_lnR, k, nvar, dLx_dlnRp1/L_scale)
               if (s% i_v /= 0) &
                  call ep1(s, xscale, i_dE_dt, s% i_v, k, nvar, dLx_dvp1/L_scale)
               if (s% lnPgas_flag) then
                  call ep1(&
                     s, xscale, i_dE_dt, s% i_lnPgas, &
                     k, nvar, dLx_dlnPgasp1_const_T/L_scale)
               else
                  call ep1(&
                     s, xscale, i_dE_dt, s% i_lnd, &
                     k, nvar, dLx_dlndp1/L_scale)
               end if
            end if
            if (k > 1) then
               if (s% lnPgas_flag) then
                  call em1(&
                     s, xscale, i_dE_dt, s% i_lnPgas, &
                     k, nvar, dLx_dlnPgasm1_const_T/L_scale)
               else
                  call em1(s, xscale, i_dE_dt, s% i_lnd, &
                     k, nvar, dLx_dlndm1/L_scale)
               end if
            end if
         end if
         
         if (do_chem .and. & 
               (s% dxdt_nuc_factor > 0d0 .or. s% mix_factor > 0d0)) then
            do j=1,s% nvar_chem
               dLx_dx = s% d_epsnuc_dx(j,k)*dm
               call e00(s, xscale, i_dE_dt, s% i_chem1+j-1, k, nvar, dLx_dx/L_scale)
            end do
         end if

         if (.not. do_visc) return
         
         ! include partials of eps_visc
      
         ! dLAV = eps_visc*dm*visc_factor            
         d_dLAV_dlnd = s% d_eps_visc_dlnd(k)*dm*visc_factor
         d_dLAV_dlnR00 = s% d_eps_visc_dlnR00(k)*dm*visc_factor
         d_dLAV_dvel00 = s% d_eps_visc_dvel00(k)*dm*visc_factor
         d_dLAV_dlnRp1 = s% d_eps_visc_dlnRp1(k)*dm*visc_factor
         d_dLAV_dvelp1 = s% d_eps_visc_dvelp1(k)*dm*visc_factor

         if (s% do_struct_hydro) then
         
            if (s% lnPgas_flag) then   
               d_dLAV_dlnPgas_const_T = d_dLAV_dlnd*s% dlnRho_dlnPgas_const_T(k)
               call e00(s, xscale, i_dE_dt, s% i_lnPgas, k, nvar, &
                        d_dLAV_dlnPgas_const_T/L_scale)
            else
               call e00(s, xscale, i_dE_dt, s% i_lnd, k, nvar, d_dLAV_dlnd/L_scale)
            end if
            
            call e00(s, xscale, i_dE_dt, s% i_lnR, k, nvar, d_dLAV_dlnR00/L_scale)

            call e00(s, xscale, i_dE_dt, s% i_v, k, nvar, d_dLAV_dvel00/L_scale)
            if (k < s% nz) then
               call ep1(s, xscale, i_dE_dt, s% i_lnR, k, nvar, d_dLAV_dlnRp1/L_scale)
               call ep1(s, xscale, i_dE_dt, s% i_v, k, nvar, d_dLAV_dvelp1/L_scale)
            end if
            
         end if

      end subroutine do1_dLdm_eqn
      
            
      subroutine do1_dPdm_eqn( &
            s, k, P_surf, xscale, equ, skip_partials, do_visc, nvar, ierr)
         use chem_def, only: chem_isos
         use eos_def, only: i_lnPgas

         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, do_visc
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
      
         real(dp) :: &
            lnP_surf, delm, d_delm_dqm1, d_delm_dq00, P00, Pm1, &
            dPr, Ppoint, Ppoint_inv, Ppoint_inv2, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            d_Ppoint_inv_dlndm1, d_Ppoint_inv_dlnd00, &
            d_Ppoint_inv_dlnTm1, d_Ppoint_inv_dlnT00, &
            d_Ppoint_inv_dlnPgas00_const_T, d_Ppoint_inv_dlnPgasm1_const_T, &
            d_Ppoint_inv_dlnT00_const_Pgas, d_Ppoint_inv_dlnTm1_const_Pgas, &
            lnPdiff, d_lnPdiff_dlnP00, d_lnPdiff_dlnPm1, var_test, equ_test, &
            d_dlnPdm_dlnRp1, d_dlnPdm_dlnR00, d_dlnPdm_dlnRm1, &
            d_dlnPdm_dvelp1, d_dlnPdm_dvel00, d_dlnPdm_dvelm1, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlndm1, &
            d_dlnPdm_dlnT00, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, & 
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, & 
            d_delm_dlndq00, d_delm_dlndqm1, diff, diff_old, &
            dA_dlnR, delQ, dlnPQ, Ainv, &
            d_delQ_dlnd00, d_delQ_dlnT00, d_delQ_dlnR00, d_delQ_dvel00, &
            d_delQ_dlndm1, d_delQ_dlnTm1, d_delQ_dlnRm1, d_delQ_dvelm1, &
            d_delQ_dlnRp1, d_delQ_dvelp1, d_Ainv_dlnR00, &
            d_delQ_dlnPgas00_const_T, d_delQ_dlnPgasm1_const_T, &
            d_delQ_dlnT00_const_Pgas, d_delQ_dlnTm1_const_Pgas, &
            d_dlnPQ_dlnd00, d_dlnPQ_dlnT00, d_dlnPQ_dlnR00, d_dlnPQ_dvel00, &
            d_dlnPQ_dlndm1, d_dlnPQ_dlnTm1, d_dlnPQ_dlnRm1, d_dlnPQ_dvelm1, &
            d_dlnPQ_dlnRp1, d_dlnPQ_dvelp1, &
            r003, rp13, T1, lnT1, lnP1, &
            dlnT1_dL, dlnT1_dlnR, &
            dlnP1_dL, dlnP1_dlnR, &
            dlnP1_dlnm, dlnP1_dlnkap, &
            dlnkap_dlnd, dlnkap_dlnT, &
            dlnP1_dlnd, dlnP1_dlnT, &
            dlnP1_dlnR00, dlnP1_dlnRp1, d_dlnPdm_dL00, &
            dlnP00_dlnPgas_const_T, dlnP00_dlnT_const_Pgas, &
            dlnPm1_dlnPgas_const_T, dlnPm1_dlnT_const_Pgas, &
            d_lnPdiff_dlnTm1_const_Pgas, d_lnPdiff_dlnT00_const_Pgas, &
            d_lnPdiff_dlnPpgas00_const_T, d_lnPdiff_dlnPgasm1_const_T, &
            d_dlnPQ_dlnPgas00_const_T, d_dlnPQ_dlnPgasm1_const_T, &
            d_dlnPQ_dlnT00_const_Pgas, d_dlnPQ_dlnTm1_const_Pgas
         integer :: j, i_dv_dt
         logical :: dbg 

         include 'formats'
         
         if (s% E_flag) then
            write(*,*) 'must not call pressure eqn when E_flag is true'
            stop 1
         end if
         
         ierr = 0
         i_dv_dt = s% i_dv_dt

         dbg = .false.
         
         if (k > 1) then
            delm = (s% dm(k) + s% dm(k-1))/2
            Pm1 = s% P(k-1)
         else
            delm = s% dm(k)/2
            Pm1 = P_surf
         end if
      
         call do_dlnP_dm(s, k, do_visc, delm, &
            d_dlnPdm_dlnRp1, d_dlnPdm_dlnR00, d_dlnPdm_dlnRm1, &
            d_dlnPdm_dvelp1, d_dlnPdm_dvel00, d_dlnPdm_dvelm1, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlndm1, &
            d_dlnPdm_dlnT00, d_dlnPdm_dlnTm1, &         
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, & 
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, & 
            d_dlnPdm_dL00, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            ierr)
         if (ierr /= 0) return            
      
         P00 = s% P(k)
         dPr = Pm1 - P00
         Ppoint_inv = 1/Ppoint
   
         ! basic eqn is dP = -G m / (4 pi r^4)
         ! divide by <P> to make it unitless
         ! simple average is adequate for <P> since is only for normalizing the equation.
         ! however, be careful to use same <P> for both sides of equation..... 
      
         lnPdiff = dPr/Ppoint ! use this in place of lnP(k-1)-lnP(k)
      
         equ(i_dv_dt, k) = delm*s% dlnP_dm_expected(k) - lnPdiff
         s% v_residual(k) = equ(i_dv_dt, k)

         if (k == s% trace_k) then
            write(*,5) 'i_dv_dt', k, s% newton_iter, s% newton_adjust_iter, &
               s% model_number, equ(i_dv_dt, k), s% dlnP_dm_expected(k)
         end if
                  
         if (is_bad_num(equ(i_dv_dt, k))) then
            ierr = -1
            if (s% report_ierr) &
               write(*,2) 'P_eqn: is_bad_num(equ(i_dv_dt, k))', k, equ(i_dv_dt,k)
            !stop
            return
         end if

         if (s% using_ode_form) s% ode(i_dv_dt,k) = 0
         
         if (skip_partials) return

         d_lnPdiff_dlnPm1 = P00*Pm1/(Ppoint*Ppoint)
         d_lnPdiff_dlnP00 = -d_lnPdiff_dlnPm1
      
         call e00(s, xscale, i_dv_dt, s% i_lnR, k, nvar, delm*d_dlnPdm_dlnR00)
         if (s% i_lum /= 0) &
            call e00(s, xscale, i_dv_dt, s% i_lum, k, nvar, delm*d_dlnPdm_dL00)
         
         if (s% lnPgas_flag) then
            dlnP00_dlnPgas_const_T = s% Pgas(k)/s% P(k)
            dlnP00_dlnT_const_Pgas = 4*s% Prad(k)/s% P(k)
            d_lnPdiff_dlnPpgas00_const_T = d_lnPdiff_dlnP00*dlnP00_dlnPgas_const_T
            d_lnPdiff_dlnT00_const_Pgas = d_lnPdiff_dlnP00*dlnP00_dlnT_const_Pgas
            call e00(s, xscale, i_dv_dt, s% i_lnPgas, k, nvar, &
                  delm*d_dlnPdm_dlnPgas00_const_T - d_lnPdiff_dlnPpgas00_const_T)
            if (s% do_struct_thermo) &
               call e00(s, xscale, i_dv_dt, s% i_lnT, k, nvar, &
                     delm*d_dlnPdm_dlnT00_const_Pgas - d_lnPdiff_dlnT00_const_Pgas)
         else
            call e00(s, xscale, i_dv_dt, s% i_lnd, k, nvar, &
                  delm*d_dlnPdm_dlnd00 - d_lnPdiff_dlnP00*s% chiRho(k))
            if (s% do_struct_thermo) &
               call e00(s, xscale, i_dv_dt, s% i_lnT, k, nvar, &
                     delm*d_dlnPdm_dlnT00 - d_lnPdiff_dlnP00*s% chiT(k))
         end if
         
         if (k < s% nz) then
            call ep1(s, xscale, i_dv_dt, s% i_lnR, k, nvar, delm*d_dlnPdm_dlnRp1)
         end if
         
         if (k > 1) then
            call em1(s, xscale, i_dv_dt, s% i_lnR, k, nvar, delm*d_dlnPdm_dlnRm1)
            if (s% lnPgas_flag) then
               dlnPm1_dlnPgas_const_T = s% Pgas(k-1)/s% P(k-1)
               dlnPm1_dlnT_const_Pgas = 4*s% Prad(k-1)/s% P(k-1)
               d_lnPdiff_dlnPgasm1_const_T = d_lnPdiff_dlnPm1*dlnPm1_dlnPgas_const_T
               d_lnPdiff_dlnTm1_const_Pgas = d_lnPdiff_dlnPm1*dlnPm1_dlnT_const_Pgas
               call em1(s, xscale, i_dv_dt, s% i_lnPgas, k, nvar, &
                  delm*d_dlnPdm_dlnPgasm1_const_T - d_lnPdiff_dlnPgasm1_const_T)
               if (s% do_struct_thermo) &
                  call em1(s, xscale, i_dv_dt, s% i_lnT, k, nvar, &
                     delm*d_dlnPdm_dlnTm1_const_Pgas - d_lnPdiff_dlnTm1_const_Pgas)
            else
               call em1(s, xscale, i_dv_dt, s% i_lnd, k, nvar, &
                  delm*d_dlnPdm_dlndm1 - d_lnPdiff_dlnPm1*s% chiRho(k-1))
               if (s% do_struct_thermo) &
                  call em1(s, xscale, i_dv_dt, s% i_lnT, k, nvar, &
                     delm*d_dlnPdm_dlnTm1 - d_lnPdiff_dlnPm1*s% chiT(k-1))
            end if
         end if

         if (s% v_flag) then
            call e00(s, xscale, i_dv_dt, s% i_v, k, nvar, delm*d_dlnPdm_dvel00)
            if (k < s% nz) &
               call ep1(s, xscale, i_dv_dt, s% i_v, k, nvar, delm*d_dlnPdm_dvelp1)
            if (k > 1) &
               call em1(s, xscale, i_dv_dt, s% i_v, k, nvar, delm*d_dlnPdm_dvelm1)
         end if
         
         if (.false.) then
            write(*,2) 'do1_dPdm_eqn equ(i_dv_dt, k)', k, equ(i_dv_dt, k)
            write(*,2) 'do1_dPdm_eqn vel00', k, delm*d_dlnPdm_dvel00
            write(*,2) 'do1_dPdm_eqn velp1', k, delm*d_dlnPdm_dvelp1
            write(*,2) 'do1_dPdm_eqn velm1', k, delm*d_dlnPdm_dvelm1
            write(*,2) 'do1_dPdm_eqn lnR00', k, delm*d_dlnPdm_dlnR00
            write(*,2) 'do1_dPdm_eqn lnRp1', k, delm*d_dlnPdm_dlnRp1
            write(*,2) 'do1_dPdm_eqn lnRm1', k, delm*d_dlnPdm_dlnRm1
            write(*,2) 'do1_dPdm_eqn lnd00', k, &
               delm*d_dlnPdm_dlnd00 - d_lnPdiff_dlnP00*s% chiRho(k)
            write(*,2) 'do1_dPdm_eqn lndm1', k, &
               delm*d_dlnPdm_dlndm1 - d_lnPdiff_dlnPm1*s% chiRho(k-1)
            write(*,2) 'do1_dPdm_eqn lnT00', k, &
               delm*d_dlnPdm_dlnT00 - d_lnPdiff_dlnP00*s% chiT(k)
            write(*,2) 'do1_dPdm_eqn lnTm1', k, &
               delm*d_dlnPdm_dlnTm1 - d_lnPdiff_dlnPm1*s% chiT(k-1)
            write(*,2) 'do1_dPdm_eqn L00', k, delm*d_dlnPdm_dL00
         end if

      
      end subroutine do1_dPdm_eqn


      subroutine PT_eqns_surf( &
            s, xscale, equ, skip_partials, do_visc, &
            convert_ODE_to_DAE_form, nvar, ierr)
         use hydro_vars, only: set_Teff_info_for_eqns
         use chem_def
         use atm_def
         use atm_lib, only: atm_option
         
         type (star_info), pointer :: s         
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, do_visc, convert_ODE_to_DAE_form
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         
         integer :: i_lnd, i_lnPgas, i_lnT, i_E, i_lnR, i_lum, i_v, i_dv_dt, i_dlnT_dm
         real(dp) :: r, L, Teff, &
            lnT_surf, dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
            lnP_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap
         real(dp) :: &
            dlnT_bc_dlnd, dlnT_bc_dlnT, dlnT_bc_dlnR, &
            dlnT_bc_dL, dlnP_bc_dlnd, dlnP_bc_dlnT, dlnP_bc_dlnR, dlnP_bc_dL, &
            dlnkap_dlnd, dlnkap_dlnT, dPinv_dlnd, dPinv_dlnT, dP0, dT0, &
            P_surf, T_surf, dlnP_bc_dlnPsurf, &
            dlnT_bc_dlnTsurf, P_bc, T_bc, lnT_bc, lnP_bc, &
            dP0_dlnR, dT0_dlnR, dT0_dlnT, dT0_dlnd, dT0_dL, dlnP_bc_dP0, dlnT_bc_dT0, &
            dlnP_bc_dlnPgas_const_T, dlnP_bc_dlnT_const_Pgas, dlnP_dlnPgas_const_T, &
            dlnP_dlnT_const_Pgas, dlnT_bc_dlnPgas_const_T, dlnT_bc_dlnT_const_Pgas, &
            dlnT_bc_dE_const_Rho, dlnT_dE_const_Rho, dlnP_dE_c_Rho, &
            dlnP_bc_dE_c_Rho, dlnT_bc_dlnd_c_E, dlnP_bc_dlnd_c_E
         
         include 'formats'
         
         ierr = 0
         s% dlnP_dm_expected(1) = s% dlnP_dm_expected(2) ! not used except for output
         s% dlnT_dm_expected(1) = s% dlnT_dm_expected(2) ! not used except for output

         i_dv_dt = s% i_dv_dt
         i_dlnT_dm = s% i_dlnT_dm

         i_lnPgas = s% i_lnPgas
         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_E = s% i_E
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_v = s% i_v
         
         if (s% use_hydrodynamic_surface_BCs) then
            call set_hydrodynamic_surface_BCs( &
               s, xscale, equ, skip_partials, convert_ODE_to_DAE_form, nvar, ierr)
            return
         end if
         
         if (i_dlnT_dm /= 0 .and. i_lum == 0) then
            write(*,*) 'error in PT_eqns_surf: i_dlnT_dm /= 0 .and. i_lum == 0'
            write(*,*) 'set use_ODE_var_eqn_pairing = .true.'
         end if

         if (s% lnPgas_flag .and. s% E_flag) then
            write(*,*) 'sorry: for surf BCs, cannot have both E_flag and lnPgas_flag'
            ierr = -1
            return
         end if
         
         call set_Teff_info_for_eqns(s, skip_partials, r, L, Teff, &
            lnT_surf, dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
            lnP_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap, &
            ierr)
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*) 'P_eqn_surf: ierr from set_Teff_info_for_eqns'
            end if
            return
         end if
         
         ! P_surf and T_surf are at outer boundary of cell 1
         P_surf = exp_cr(lnP_surf)
         T_surf = exp_cr(lnT_surf)
         s% P_surf = P_surf
         s% T_surf = T_surf
         
         if (s% use_atm_PT_at_center_of_surface_cell .or. &
             s% using_free_fall_surface_PT) then
            dP0 = 0
            dT0 = 0
         else ! offset P and T from outer edge of cell 1 to center of cell 1
            dP0 = s% cgrav(1)*s% mstar*s% dm(1)/(8*pi*r*r*r*r)
            dT0 = dP0*s% gradT(1)*s% T(1)/s% P(1)
         end if
         
         P_bc = P_surf + dP0
         T_bc = T_surf + dT0

         lnP_bc = log_cr(P_bc)
         lnT_bc = log_cr(T_bc)
         
         if (s% do_struct_hydro) then
            if (.not. s% use_momentum_eqn_outer_BC) then
               equ(i_dv_dt, 1) = lnP_bc - s% lnP(1)
            else if (s% use_dvdt_form_of_momentum_eqn) then
               call do1_dvdt_eqn( &
                  s, 1, P_surf, xscale, equ, do_visc, &
                  skip_partials, convert_ODE_to_DAE_form, &
                  nvar, ierr)
            else
               call do1_dPdm_eqn( &
                  s, 1, P_surf, xscale, equ, do_visc, &
                  skip_partials, nvar, ierr)
            end if
         else
            equ(i_dv_dt,1) = 0
         end if
         
         s% v_residual(1) = equ(i_dv_dt,1)
         if (s% using_ode_form) s% ode(i_dv_dt,1) = 0
         
         if (i_dlnT_dm /= 0) then
            if (s% do_struct_thermo) then
               equ(i_dlnT_dm, 1) = lnT_bc - s% lnT(1)        
            else
               equ(i_dlnT_dm,1) = 0
            end if
            if (s% using_ode_form) s% ode(i_dlnT_dm,1) = 0
         end if

         if (skip_partials) return

         !dT0 = dP0*s% gradT(1)*s% T(1)/s% P(1)
         
         dP0_dlnR = -4*dP0         
         dT0_dlnR = -4*dT0 + dP0*s% d_gradT_dlnR(1)*s% T(1)/s% P(1)
         
         dPinv_dlnT = -s% chiT(1)/s% P(1)
         dT0_dlnT = &
            dT0 + &
            dP0*s% d_gradT_dlnT00(1)*s% T(1)/s% P(1) + &
            dP0*s% gradT(1)*s% T(1)*dPinv_dlnT
         dPinv_dlnd = -s% chiRho(1)/s% P(1)
         dT0_dlnd = &
            dP0*s% d_gradT_dlnd00(1)*s% T(1)/s% P(1) + &
            dP0*s% gradT(1)*s% T(1)*dPinv_dlnd
            
         dT0_dL = dP0*s% d_gradT_dL(1)*s% T(1)/s% P(1)
         
         dlnP_bc_dP0 = 1/P_bc
         dlnT_bc_dT0 = 1/T_bc
         
         dlnP_bc_dlnPsurf = P_surf/P_bc
         dlnT_bc_dlnTsurf = T_surf/T_bc

         dlnkap_dlnd = s% d_opacity_dlnd(1)/s% opacity(1)
         dlnkap_dlnT = s% d_opacity_dlnT(1)/s% opacity(1)

         dlnP_bc_dlnd = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnd
         dlnP_bc_dlnT = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnT

         dlnT_bc_dlnT = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnT &
               + dlnT_bc_dT0*dT0_dlnT
         dlnT_bc_dlnd = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnd &
               + dlnT_bc_dT0*dT0_dlnd
         
         if (i_dlnT_dm /= 0) then
         
            if (s% do_struct_thermo) then ! temperature eqn
               ! equ(i_dlnT_dm, 1) = lnT_bc - s% lnT(1)        
         
               if (s% lnPgas_flag) then
                  dlnT_bc_dlnT_const_Pgas = &
                     dlnT_bc_dlnT + dlnT_bc_dlnd*s% dlnRho_dlnT_const_Pgas(1)
                  call e00(s, xscale, i_dlnT_dm, i_lnT, &
                     1, nvar, dlnT_bc_dlnT_const_Pgas - 1)
               else if (s% E_flag) then
                  dlnT_dE_const_Rho = s% dlnT_dlnE_c_Rho(1)/s% energy(1)
                  dlnT_bc_dE_const_Rho = dlnP_bc_dlnT*dlnT_dE_const_Rho
                  call e00(s, xscale, i_dlnT_dm, i_E, &
                     1, nvar, dlnT_bc_dE_const_Rho - dlnT_dE_const_Rho)
               else
                  call e00(s, xscale, i_dlnT_dm, i_lnT, 1, nvar, dlnT_bc_dlnT - 1)
               end if
               
               dlnT_bc_dL = dlnT_bc_dlnTsurf*dlnTsurf_dL + dlnT_bc_dT0*dT0_dL
               call e00(s, xscale, i_dlnT_dm, i_lum, 1, nvar, dlnT_bc_dL)
            
               if (s% do_struct_hydro) then
                  ! partial of temperature eqn wrt (lnR or v) and (lnd or lnPgas)
            
                  dlnT_bc_dlnR = dlnT_bc_dlnTsurf*dlnTsurf_dlnR + dlnT_bc_dT0*dT0_dlnR
                  call e00(s, xscale, i_dlnT_dm, i_lnR, 1, nvar, dlnT_bc_dlnR)

                  if (s% lnPgas_flag) then
                     dlnT_bc_dlnPgas_const_T = dlnT_bc_dlnd*s% dlnRho_dlnPgas_const_T(1)
                     call e00(s, xscale, i_dlnT_dm, i_lnPgas, &
                        1, nvar, dlnT_bc_dlnPgas_const_T)
                  else if (s% E_flag) then
                     dlnT_bc_dlnd_c_E = dlnP_bc_dlnd + dlnP_bc_dlnT*s% dlnT_dlnd_c_E(1)
                     call e00(s, xscale, i_dlnT_dm, i_lnd, 1, nvar, &
                        dlnT_bc_dlnd_c_E - s% dlnT_dlnd_c_E(1))
                  else
                     call e00(s, xscale, i_dlnT_dm, i_lnd, 1, nvar, dlnT_bc_dlnd)
                  end if
               
               end if
            
            else ! dummy eqn
         
               call e00(s,xscale,i_dlnT_dm,i_lnT,1,nvar,one)
            
            end if
            
         end if

         if (s% do_struct_hydro .and. .not. s% use_momentum_eqn_outer_BC) then ! pressure eqn
            ! equ(i_dv_dt, 1) = lnP_bc - s% lnP(1)
            
            if (s% lnPgas_flag) then
               dlnP_bc_dlnPgas_const_T = dlnP_bc_dlnd*s% dlnRho_dlnPgas_const_T(1)
               dlnP_dlnPgas_const_T = s% Pgas(1)/s% P(1)
               call e00(s, xscale, i_dv_dt, i_lnPgas, 1, nvar, &
                  dlnP_bc_dlnPgas_const_T - dlnP_dlnPgas_const_T)
            else if (s% E_flag) then
               dlnP_bc_dlnd_c_E = dlnP_bc_dlnd + dlnP_bc_dlnT*s% dlnT_dlnd_c_E(1)
               call e00(s, xscale, i_dv_dt, i_lnd, &
                  1, nvar, dlnP_bc_dlnd_c_E - s% dlnP_dlnd_c_E(1))
            else
               call e00(s, xscale, i_dv_dt, i_lnd, 1, nvar, dlnP_bc_dlnd - s% chiRho(1))
            end if
            
            dlnP_bc_dlnR = dlnP_bc_dlnPsurf*dlnPsurf_dlnR + dlnP_bc_dP0*dP0_dlnR
            
            call e00(s, xscale, i_dv_dt, i_lnR, 1, nvar, dlnP_bc_dlnR)
            
            if (s% do_struct_thermo) then 
               
               if (s% lnPgas_flag) then
                  dlnP_bc_dlnT_const_Pgas = &
                     dlnP_bc_dlnT + dlnP_bc_dlnd*s% dlnRho_dlnT_const_Pgas(1)
                  dlnP_dlnT_const_Pgas = 4*s% Prad(1)/s% P(1)
                  call e00(s, xscale, i_dv_dt, i_lnT, 1, &
                     nvar, dlnP_bc_dlnT_const_Pgas - dlnP_dlnT_const_Pgas)
               else if (s% E_flag) then
                  dlnP_bc_dE_c_Rho = dlnP_bc_dlnT*s% dlnT_dlnE_c_Rho(1)/s% energy(1)
                  dlnP_dE_c_Rho = s% dlnP_dlnE_c_Rho(1)/s% energy(1)
                  call e00(s, xscale, i_dv_dt, i_E, 1, &
                     nvar, dlnP_bc_dE_c_Rho - dlnP_dE_c_Rho)
               else
                  call e00(s, xscale, i_dv_dt, i_lnT, 1, &
                     nvar, dlnP_bc_dlnT - s% chiT(1))
               end if
               
               if (i_lum /= 0) then
                  dlnP_bc_dL = dlnP_bc_dlnPsurf*dlnPsurf_dL
                  call e00(s, xscale, i_dv_dt, i_lum, 1, nvar, dlnP_bc_dL)
               end if
               
            end if
            
         else ! dummy eqn
            
            if (s% lnPgas_flag) then
               call e00(s,xscale,i_dv_dt,i_lnPgas,1,nvar,one)
            else
               call e00(s,xscale,i_dv_dt,i_lnd,1,nvar,one)
            end if
            
         end if

      end subroutine PT_eqns_surf


      subroutine eval_grad_star_info( &
            s, k, grad_star, d_grad_star_dL, d_grad_star_dlnR, &
            d_grad_star_dlnd00, d_grad_star_dlndm1, &
            d_grad_star_dlnT00, d_grad_star_dlnTm1, &
            d_grad_star_dlnT00_const_Pgas, d_grad_star_dlnTm1_const_Pgas, &
            d_grad_star_dlnPgas00_const_T, d_grad_star_dlnPgasm1_const_T, &
            ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k        
         real(dp), intent(out) :: grad_star, d_grad_star_dL, d_grad_star_dlnR, &
            d_grad_star_dlnd00, d_grad_star_dlndm1, &
               d_grad_star_dlnT00, d_grad_star_dlnTm1, &
            d_grad_star_dlnT00_const_Pgas, d_grad_star_dlnTm1_const_Pgas, &
            d_grad_star_dlnPgas00_const_T, d_grad_star_dlnPgasm1_const_T
         integer, intent(out) :: ierr         
         include 'formats'
              
         ierr = 0      
         
         grad_star = s% gradT(k)  
          
         d_grad_star_dL = s% d_gradT_dL(k)   
         d_grad_star_dlnR = s% d_gradT_dlnR(k)

         d_grad_star_dlnd00 = s% d_gradT_dlnd00(k)
         d_grad_star_dlndm1 = s% d_gradT_dlndm1(k)   
         d_grad_star_dlnT00 = s% d_gradT_dlnT00(k)
         d_grad_star_dlnTm1 = s% d_gradT_dlnTm1(k)         
         
         if (s% lnPgas_flag) then
            d_grad_star_dlnPgas00_const_T = &
               d_grad_star_dlnd00*s% dlnRho_dlnPgas_const_T(k)
            d_grad_star_dlnT00_const_Pgas = &
               d_grad_star_dlnT00 + d_grad_star_dlnd00*s% dlnRho_dlnT_const_Pgas(k)
            if (k > 1) then
               d_grad_star_dlnPgasm1_const_T = &
                  d_grad_star_dlndm1*s% dlnRho_dlnPgas_const_T(k-1)
               d_grad_star_dlnTm1_const_Pgas = &
                  d_grad_star_dlnTm1 + d_grad_star_dlndm1*s% dlnRho_dlnT_const_Pgas(k-1)
            else
               d_grad_star_dlnPgasm1_const_T = 0
               d_grad_star_dlnTm1_const_Pgas = 0
            end if
         end if

      end subroutine eval_grad_star_info


      subroutine set_hydrodynamic_surface_BCs( &
            s, xscale, equ, skip_partials, convert_ODE_to_DAE_form, nvar, ierr)
         use star_utils, only: get_r_phot
         type (star_info), pointer :: s         
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, convert_ODE_to_DAE_form
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         
         real(dp) :: L_surf, r_phot, Teff, T1, dT1_dL, dT1_dlnR, &
            Lmid, rmid, dT1_dlnR00, dT1_dlnRp1, dT1_dL00, dT1_dLp1, Tscale
         integer :: i_lnd, i_lnT, i_E, i_lnR, i_lum, i_dv_dt, i_dlnT_dm
         
         include 'formats'

         ierr = 0

         i_dv_dt = s% i_dv_dt
         i_dlnT_dm = s% i_dlnT_dm

         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_E = s% i_E
         i_lnR = s% i_lnR
         i_lum = s% i_lum

         if (s% lnPgas_flag) then
            write(*,*) 'cannot use hydrodynamic_surface_BCs with lnPgas_flag'
            ierr = -1
            return
         end if
         
         if (i_lum /= 0) then
         
            L_surf = s% L(1)
            if (L_surf <= 0) then
               ierr = -1
               if (s% report_ierr) &
                  write(*,1) 'L_surf <= 0 for hydrodynamic_surface_BCs', L_surf
               return
            end if
         
            r_phot = get_r_phot(s)
            Teff = pow_cr(L_surf/(4*pi*r_phot*r_phot*boltz_sigma), 0.25d0)
            s% Teff = Teff
            s% L_phot = L_surf/Lsun
            s% photosphere_L = s% L_phot
            s% photosphere_r = r_phot/Rsun
            s% P_surf = s% P(1)
            s% T_surf = s% T(1)
            
         else
         
            s% L_phot = 0
            s% photosphere_L = 0
            s% photosphere_r = s% r(1)
            s% P_surf = s% P(1)
            s% T_surf = s% T(1)
            s% Teff = s% T_surf
            
         end if
         
         if (s% do_struct_hydro) then
            ! gradient of compression vanishes fixes density for cell 1
               ! d_dt(1/rho(1)) = d_dt(1/rho(2))
               ! Grott, Chernigovski, Glatzel, 2005.
            equ(i_dv_dt, 1) = s% rho_start(1)*( &
               (1/s% rho(1) - 1/s% rho_start(1)) - (1/s% rho(2) - 1/s% rho_start(2)))
            if (is_bad_num(equ(i_dv_dt, 1))) then
               write(*,1) 'equ(i_dv_dt, 1)', equ(i_dv_dt, 1)
               stop 'hydrodynamic_surface_BCs'
            end if
            if (s% using_ode_form) s% ode(i_dv_dt,1) = 0
            s% d_KE_dt_actual(1) = 0 ! clear these for plotting
            s% d_KE_dt_expected(1) = 0
            s% v_residual(1) = 0
         end if
         
         if (i_dlnT_dm /= 0) then
            
            if (.not. s% L_flag) &
               stop 'set_hydrodynamic_surface_BCs: must have L_flag when i_dlnT_dm /= 0'
               
            rmid = s% rmid(1)
            Lmid = (s% L(1) + s% L(2))/2
            if (Lmid <= 0) then
               if (s% report_ierr) &
                  write(*,2) 'Lmid <= 0 for hydrodynamic_surface_BCs', s% model_number
               ierr = -1
               return
            end if
            
            T1 = pow_cr(Lmid/(4*pi*rmid*rmid*boltz_sigma), 0.25d0)
            dT1_dL = T1/(4*Lmid)
            dT1_dL00 = dT1_dL/2
            dT1_dLp1 = dT1_dL/2
            dT1_dlnR = -T1/2
            dT1_dlnR00 = dT1_dlnR*s% drmid_dlnR00(1)/rmid
            dT1_dlnRp1 = dT1_dlnR*s% drmid_dlnRp1(1)/rmid
         
            if (s% do_struct_thermo) then         
               Tscale = 1d6*s% T_start(1) ! 1d6 to reduce equ size compared to k > 1 cases
               equ(i_dlnT_dm, 1) = (T1 - s% T(1))/Tscale 
               if (is_bad_num(equ(i_dlnT_dm, 1))) then
                  write(*,1) 'equ(i_dlnT_dm, 1)', equ(i_dlnT_dm, 1)
                  write(*,1) 'L_surf', L_surf
                  write(*,1) 's% r(1)', s% r(1)
                  write(*,1) 'T1', T1
                  write(*,1) 'Lmid', Lmid
                  write(*,1) 'rmid', rmid
                  write(*,1) 's% L(1)', s% L(1)
                  write(*,1) 's% L(2)', s% L(2)
                  write(*,1) 's% r(1)', s% r(1)
                  write(*,1) 's% r(2)', s% r(2)
                  write(*,1) 's% T(1)', s% T(1)
                  write(*,1) 's% T_start(1)', s% T_start(1)
                  write(*,1) 'Tscale', Tscale
                  stop 'hydrodynamic_surface_BCs'
               end if
            else
               equ(i_dlnT_dm,1) = 0
            end if
            s% dlnTdm_residual(1) = equ(i_dlnT_dm,1)
            if (s% using_ode_form) s% ode(i_dlnT_dm,1) = 0
            
         else
         
            dT1_dL = 0
            dT1_dL00 = 0
            dT1_dLp1 = 0
            dT1_dlnR = 0
            dT1_dlnR00 = 0
            dT1_dlnRp1 = 0
            
         end if

         if (skip_partials) return
         
         if (s% do_struct_hydro) then
            call e00(s, xscale, i_dv_dt, i_lnd, 1, nvar, -s% rho_start(1)/s% rho(1))
            call ep1(s, xscale, i_dv_dt, i_lnd, 1, nvar,  s% rho_start(1)/s% rho(2))
         end if

         if (i_dlnT_dm /= 0) then
            if (s% do_struct_thermo) then ! uses T(1), L(1), L(2), r(1), r(2)            
               if (s% E_flag) then ! d(-s% T(1)/s% T_start(1))/dE              
                  call e00(s, xscale, i_dlnT_dm, i_lnd, 1, nvar, &
                     -s% dlnT_dlnd_c_E(1)/Tscale)
                  call e00(s, xscale, i_dlnT_dm, i_E, 1, nvar, &
                     -s% dlnT_dlnE_c_Rho(1)/s% energy(1)/Tscale)                 
               else
                  call e00(s, xscale, i_dlnT_dm, i_lnT, 1, nvar, -s% T(1)/Tscale)
               end if               
               if (i_lum /= 0) then
                  call e00(s, xscale, i_dlnT_dm, i_lum, 1, nvar, dT1_dL00/Tscale)
                  call ep1(s, xscale, i_dlnT_dm, i_lum, 1, nvar, dT1_dLp1/Tscale)
               end if
               
               if (s% do_struct_hydro) then
                  call e00(s, xscale, i_dlnT_dm, i_lnR, 1, nvar, dT1_dlnR00/Tscale)
                  call ep1(s, xscale, i_dlnT_dm, i_lnR, 1, nvar, dT1_dlnRp1/Tscale)
               end if    
            else       
               call e00(s, xscale, i_dlnT_dm, i_lnT, 1, nvar, one)            
            end if
         end if

      end subroutine set_hydrodynamic_surface_BCs


      end module hydro_eqns

