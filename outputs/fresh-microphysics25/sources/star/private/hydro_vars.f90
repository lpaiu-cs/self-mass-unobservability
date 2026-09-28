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

      module hydro_vars
      
      use star_private_def
      use const_def
      use chem_def, only: chem_isos
      use utils_lib, only:is_bad_num
      
      implicit none

      logical, parameter :: dbg = .false.      

      contains
      
      
      subroutine set_vars(s, dt, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         logical, parameter :: &
            skip_time_derivatives = .false., &
            skip_m_grav_and_grav = .false., &
            skip_net = .false., &
            skip_neu = .false., &
            skip_kap = .false., &
            skip_grads = .false., &
            skip_rotation = .false., &
            skip_brunt = .false., &
            skip_mixing_info = .false.
         if (dbg) write(*, *) 'set_vars'
         call set_some_vars( &
            s, skip_time_derivatives, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, dt, ierr)
      end subroutine set_vars
      
      
      subroutine set_vars_before_diffusion(s, dt, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         call set_vars(s, dt, ierr)
      end subroutine set_vars_before_diffusion
      
      
      subroutine reset_vars_after_diffusion(s, dt, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         logical, parameter :: &
            skip_time_derivatives = .true., &
            skip_m_grav_and_grav = .true., &
            skip_net = .false., &
            skip_neu = .false., &
            skip_kap = .false., &
            skip_grads = .false., &
            skip_rotation = .false., &
            skip_brunt = .false., &
            skip_mixing_info = .true.
         if (dbg) write(*, *) 'reset_vars_after_diffusion'
         call set_some_vars( &
            s, skip_time_derivatives, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, dt, ierr)
      end subroutine reset_vars_after_diffusion
      
      
      subroutine set_final_vars(s, dt, ierr)
         use alloc, only: non_crit_get_work_array, non_crit_return_work_array
         use rates_def, only: num_rvs
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         
         logical, parameter :: &
            skip_time_derivatives = .true., &
            skip_m_grav_and_grav = .false., &
            skip_net = .false., &
            skip_neu = .false., &
            skip_kap = .false., &
            skip_grads = .false., &
            skip_rotation = .true., &
            skip_brunt = .false., &
            skip_mixing_info = .true.
         integer :: nz, ierr1, k  
         real(dp), pointer, dimension(:) :: &
            d_mx, cv, nu_st, d_st, d_dsi, d_sh, d_ssi, d_es, d_gsf
         
         ierr = 0
         nz = s% nz
         
         nullify(d_mx, cv, nu_st, d_st, d_dsi, d_sh, d_ssi, d_es, d_gsf)
                  
         ! save and restore mixing coeffs needed for time smoothing
         if (s% rotation_flag) then
            call get_cpy(s% nu_ST, nu_st, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_ST, d_st, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_DSI, d_dsi, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_SH, d_sh, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_SSI, d_ssi, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_ES, d_es, ierr1)
            if (ierr1 /= 0) ierr = -1
            call get_cpy(s% D_GSF, d_gsf, ierr1)
            if (ierr1 /= 0) ierr = -1
         end if

         call get_cpy(s% D_mix, d_mx, ierr1)
         if (ierr1 /= 0) ierr = -1
         
         if (s% min_T_for_acceleration_limited_conv_velocity < 1d12) then
            call get_cpy(s% conv_vel, cv, ierr1)
            if (ierr1 /= 0) ierr = -1
         end if
         
         if (ierr == 0) then
            if (dbg) write(*, *) 'set_vars'
            
            ! save values used during the step before they are replaced
            ! by values at end of step
            
            s% n_conv_regions_older = s% n_conv_regions_old
            s% n_conv_regions_old = s% n_conv_regions
         
            s% cz_bot_mass_older(:) = s% cz_bot_mass_old(:)
            s% cz_bot_mass_old(:) = s% cz_bot_mass(:)
         
            s% cz_top_mass_older(:) = s% cz_top_mass_old(:)
            s% cz_top_mass_old(:) = s% cz_top_mass(:)

            call set_some_vars( &
               s, skip_time_derivatives, skip_m_grav_and_grav, &
               skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, skip_brunt, &
               skip_mixing_info, dt, ierr)
            
         end if
         
         if (ierr == 0) then
         
            if (s% rotation_flag) then
               call restore(s% nu_ST, nu_st)
               call restore(s% D_ST, d_st)
               call restore(s% D_DSI, d_dsi)
               call restore(s% D_SH, d_sh)
               call restore(s% D_SSI, d_ssi)
               call restore(s% D_ES, d_es)
               call restore(s% D_GSF, d_gsf)         
               s% have_previous_rotation_info = .true.
            end if
         
            if (s% min_T_for_acceleration_limited_conv_velocity < 1d12) then
               call restore(s% conv_vel, cv)
               s% have_previous_conv_vel = .true.
               s% use_previous_conv_vel_from_file = .false.
            end if
         
            call restore(s% D_mix, d_mx)
            s% have_previous_D_mix = .true.
            
         end if
            
            
         contains

         
         subroutine get_cpy(src,cpy,ierr)
            real(dp), pointer, dimension(:) :: src, cpy
            integer, intent(out) :: ierr
            integer :: k
            ierr = 0
            call non_crit_get_work_array( &
               s, cpy, nz, nz_alloc_extra, 'set_final_vars', ierr)
            if (ierr /= 0) return
            do k=1,nz
               cpy(k) = src(k)
            end do
         end subroutine get_cpy
         
         
         subroutine restore(src,cpy)
            real(dp), pointer, dimension(:) :: src, cpy
            integer :: k
            if (.not. associated(cpy)) return
            do k=1,nz
               src(k) = cpy(k)
            end do
            call non_crit_return_work_array(s, cpy, 'set_final_vars')
         end subroutine restore
         
         
      end subroutine set_final_vars

      
      subroutine set_some_vars( &
            s, skip_time_derivatives, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_mixing_info, dt, ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         logical, intent(in) :: &
            skip_time_derivatives, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, &
            skip_rotation, skip_brunt, skip_mixing_info
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         integer :: time0, clock_rate
         real(dp) :: total
         logical, parameter :: skip_other_cgrav = .false.
         logical, parameter :: skip_basic_vars = .false.
         logical, parameter :: skip_micro_vars = .false.
         
         include 'formats'

         call update_vars(s, &
            skip_time_derivatives, skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, &
            skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, dt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr .or. dbg) &
               write(*, *) 'set_some_vars: update_vars returned ierr', ierr
            return
         end if

      end subroutine set_some_vars
      
      
      subroutine update_vars(s, &
            skip_time_derivatives, skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, &
            skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, dt, ierr)
         use utils_lib, only: has_bad_num, is_bad_num
         use mlt_def, only: no_mixing
         use star_utils, only: &
            eval_irradiation_heat, set_dm_bar, set_m_and_dm
         use mix_info, only: set_mixing_info
         type (star_info), pointer :: s 
         logical, intent(in) :: &
            skip_time_derivatives, skip_basic_vars, skip_micro_vars, &
            skip_m_grav_and_grav, skip_net, skip_neu, skip_kap, skip_grads, &
            skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         
         integer :: i_lnd, i_lnPgas, i_lnT, i_E, i_lnR, i_lum, i_v, &
            j, k, species, nvar_chem, nzlo, nz, k_below_just_added
         real(dp) :: dt_inv, del_t, v, start_time, end_time, &
            q00, qp1, qmin, qmax, qtop, qbot, extra
         logical, parameter :: skip_mlt = .false.
         
         include 'formats'
         
         ierr = 0
         nz = s% nz
         nzlo = 1
         k_below_just_added = s% k_below_just_added
         species = s% species
         nvar_chem = s% nvar_chem
         i_lnd = s% i_lnd
         i_lnPgas = s% i_lnPgas
         i_lnT = s% i_lnT
         i_E = s% i_E
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_v = s% i_v    
         
         do j=1,s% nvar_hydro
            if (j == i_lnd) then
               do k=nzlo,nz
                  s% lnd(k) = s% xh(i_lnd,k)
               end do
            else if (j == i_lnPgas) then
               do k=nzlo,nz
                  s% lnPgas(k) = s% xh(i_lnPgas,k)
               end do
            else if (j == i_lnT) then
               do k=nzlo,nz
                  s% lnT(k) = s% xh(i_lnT,k)
               end do
            else if (j == i_E) then
               do k=nzlo,nz
                  s% energy(k) = s% xh(i_E,k)
               end do
            else if (j == i_lnR) then
               do k=nzlo,nz
                  s% lnR(k) = s% xh(i_lnR,k)
               end do
            else if (j == i_lum) then
               do k=nzlo,nz
                  s% L(k) = s% xh(i_lum, k)
               end do
            else if (j == i_v) then
               do k=nzlo,nz
                  s% v(k) = s% xh(i_v,k)
               end do
            end if
         end do
         
         if (s% trace_k > 0 .and. s% trace_k <= nz) then
            k = s% trace_k
            if (i_lnd /= 0) &
               write(*,5) 'update_vars: lnd', k, s% newton_iter, s% newton_adjust_iter, &
                        s% model_number, s% lnd(k)
            if (i_lnPgas /= 0) &
               write(*,5) 'update_vars: lnPgas', k, s% newton_iter, s% newton_adjust_iter, &
                        s% model_number, s% lnPgas(k)
            if (i_lnT /= 0) &
               write(*,5) 'update_vars: lnT', k, s% newton_iter, s% newton_adjust_iter, &
                           s% model_number, s% lnT(k)
            if (i_E /= 0) &
               write(*,5) 'update_vars: E', k, s% newton_iter, s% newton_adjust_iter, &
                           s% model_number, s% energy(k)
            if (s% L_flag) &
               write(*,5) 'update_vars: L', k, s% newton_iter, s% newton_adjust_iter, &
                        s% model_number, s% L(k)
            write(*,5) 'update_vars: lnR', k, s% newton_iter, s% newton_adjust_iter, &
                     s% model_number, s% lnR(k)
            if (i_v /= 0) &
               write(*,5) 'update_vars: v', k, s% newton_iter, s% newton_adjust_iter, &
                        s% model_number, s% v(k)
         end if
                  
         if (i_v == 0) then
            do k=nzlo,nz
               s% v(k) = 0d0
            end do
         end if
         
         call set_m_and_dm(s)
         call set_dm_bar(s% nz, s% dm, s% dm_bar)

         if (.not. skip_mixing_info) then
            s% mixing_type(nzlo:nz) = no_mixing
            s% adjust_mlt_gradT_fraction(nzlo:nz) = -1
         end if
         
         if (.not. skip_time_derivatives) then
         
            ! time derivatives at constant q
            if (s% generations < 2 .or. dt <= 0 .or. s% nz /= s% nz_old) then
         
               if (i_lnd /= 0) s% dlnd_dt_const_q(nzlo:nz) = 0
               if (i_lnPgas /= 0) s% dlnPgas_dt_const_q(nzlo:nz) = 0
               if (i_lnT /= 0) s% dlnT_dt_const_q(nzlo:nz) = 0
               if (i_E /= 0) s% dE_dt_const_q(nzlo:nz) = 0
            
            else
         
               dt_inv = 1/dt
               if (i_lnd /= 0) then
                  do k=1,nz 
                     s% dlnd_dt_const_q(k) = &
                        (s% xh(i_lnd,k) - s% lnd_for_d_dt_const_q(k))*dt_inv
                     if (is_bad_num(s% dlnd_dt_const_q(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnd_dt_const_q', k, &
                              s% dlnd_dt_const_q(k), &
                              s% lnd_for_d_dt_const_q(k), s% xh(i_lnd,k)
                        return
                        stop 'update_vars'
                     end if
                  end do
               end if
               
               if (i_lnPgas /= 0) then
                  do k=1,nz 
                     s% dlnPgas_dt_const_q(k) = &
                        (s% xh(i_lnPgas,k) - s% lnPgas_for_d_dt_const_q(k))*dt_inv
                     if (is_bad_num(s% dlnPgas_dt_const_q(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnPgas_dt_const_q', k, &
                              s% dlnPgas_dt_const_q(k), &
                              s% lnPgas_for_d_dt_const_q(k), s% xh(i_lnPgas,k)
                        return
                     end if
                  end do
               end if
               
               if (i_lnT /= 0) then
                  do k=1,nz 
                     s% dlnT_dt_const_q(k) = &
                        (s% xh(i_lnT,k) - s% lnT_for_d_dt_const_q(k))*dt_inv
                     if (is_bad_num(s% dlnT_dt_const_q(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnT_dt_const_q', k, &
                              s% dlnT_dt_const_q(k), &
                              s% lnT_for_d_dt_const_q(k), s% xh(i_lnT,k)
                        return
                     end if
                  end do 
               end if              
               
               if (i_E /= 0) then
                  do k=1,nz 
                     s% dE_dt_const_q(k) = &
                        (s% xh(i_E,k) - s% E_for_d_dt_const_q(k))*dt_inv
                     if (is_bad_num(s% dE_dt_const_q(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dE_dt_const_q', k, &
                              s% dE_dt_const_q(k), &
                              s% E_for_d_dt_const_q(k), s% xh(i_E,k)
                        return
                     end if
                  end do 
               end if              
               
            end if
         
            ! time derivatives at constant mass
            if (s% generations < 2 .or. dt <= 0 .or. s% nz /= s% nz_old) then
         
               if (i_lnd /= 0) s% dlnd_dt(nzlo:nz) = 0
               if (i_lnPgas /= 0) s% dlnPgas_dt(nzlo:nz) = 0
               if (i_lnT /= 0) s% dlnT_dt(nzlo:nz) = 0
               if (i_E /= 0) s% dE_dt(nzlo:nz) = 0
               s% dlnR_dt(nzlo:nz) = 0
               if (s% v_flag) s% dv_dt(nzlo:nz) = 0
               s% dVARdot_dVAR = 0
            
            else
         
               dt_inv = 1/dt
               s% dVARdot_dVAR = dt_inv
            
               do k=1,k_below_just_added-1         
                  del_t = s% del_t_for_just_added(k)
                  if (k == 1 .or. del_t <= 0d0) then
                     s% dlnR_dt(k) = 0
                     if (i_v /= 0) s% dv_dt(k) = 0       
                     if (i_lnd /= 0) s% dlnd_dt(k) = 0        
                     if (i_lnPgas /= 0) s% dlnPgas_dt(k) = 0        
                     if (i_lnT /= 0) s% dlnT_dt(k) = 0  
                     if (i_E /= 0) s% dE_dt(k) = 0  
                  else
                     s% dlnR_dt(k) = (s% lnR(k) - s% lnR(1))/del_t
                     if (i_v /= 0) s% dv_dt(k) = (s% v(k) - s% v(1))/del_t            
                     if (i_lnd /= 0) s% dlnd_dt(k) = (s% lnd(k) - s% lnd(1))/del_t            
                     if (i_lnPgas /= 0) s% dlnPgas_dt(k) = &
                        (s% lnPgas(k) - s% lnPgas(1))/del_t            
                     if (i_lnT /= 0) s% dlnT_dt(k) = (s% lnT(k) - s% lnT(1))/del_t            
                     if (i_E /= 0) s% dE_dt(k) = (s% energy(k) - s% energy(1))/del_t            
                  end if           
               end do
               
               if (i_lnd /= 0) then
                  do k=k_below_just_added,nz 
                     s% dlnd_dt(k) = &
                        (s% xh(i_lnd,k) - s% lnd_for_d_dt_const_m(k))*dt_inv
                     if (is_bad_num(s% dlnd_dt(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnd_dt', k, &
                              s% dlnd_dt(k), s% lnd_for_d_dt_const_m(k), s% xh(i_lnd,k)
                        return
                        stop 'update_vars'
                     end if
                  end do
               end if
               
               if (i_lnPgas /= 0) then
                  do k=k_below_just_added,nz 
                     s% dlnPgas_dt(k) = &
                        (s% xh(i_lnPgas,k) - s% lnPgas_for_d_dt_const_m(k))*dt_inv
                     if (is_bad_num(s% dlnPgas_dt(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnPgas_dt', k, &
                              s% dlnPgas_dt(k), &
                              s% lnPgas_for_d_dt_const_m(k), s% xh(i_lnPgas,k)
                     end if
                  end do
               end if
               
               if (i_lnT /= 0) then
                  do k=k_below_just_added,nz 
                     s% dlnT_dt(k) = (s% xh(i_lnT,k) - s% lnT_for_d_dt_const_m(k))*dt_inv
                     if (is_bad_num(s% dlnT_dt(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dlnT_dt', k, s% dlnT_dt(k)
                     end if
                  end do
               end if
               
               if (i_E /= 0) then
                  do k=k_below_just_added,nz 
                     s% dE_dt(k) = (s% xh(i_E,k) - s% E_for_d_dt_const_m(k))*dt_inv
                     if (is_bad_num(s% dE_dt(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dE_dt', k, s% dE_dt(k)
                     end if
                  end do
               end if
               
               if (s% v_flag) then
                  do k=k_below_just_added,nz 
                     s% dv_dt(k) = (s% xh(i_v,k) - s% v_for_d_dt_const_m(k))*dt_inv
                     if (is_bad_num(s% dv_dt(k))) then
                        ierr = -1
                        if (s% report_ierr) &
                           write(*,2) 'update_vars: bad dv_dt', k, s% dv_dt(k)
                     end if
                  end do
               end if
               
               do k=k_below_just_added,nz
                  s% dlnR_dt(k) = (s% xh(i_lnR,k) - s% lnR_for_d_dt_const_m(k))*dt_inv
                  if (is_bad_num(s% dlnR_dt(k))) then
                     ierr = -1
                     if (s% report_ierr) &
                        write(*,2) 'update_vars: bad dlnR_dt', k, s% dlnR_dt(k)
                  end if
               end do
            
            end if
            
         end if

         call set_hydro_vars( &
            s, nzlo, nz, skip_basic_vars, &
            skip_micro_vars, skip_m_grav_and_grav, skip_net, skip_neu, &
            skip_kap, skip_grads, skip_rotation, skip_brunt, skip_other_cgrav, &
            skip_mixing_info, skip_mlt, ierr)
         if (ierr /= 0) return         

         if (s% Teff <= 0) then
            call set_Teff(s, ierr) 
            if (ierr /= 0) return
         end if

         if (s% use_other_momentum) then
            call s% other_momentum(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*, *) 'set_vars: other_momentum returned ierr', ierr
               return
            end if
         end if
         
         if (.not. s% use_other_energy_implicit) then
            s% extra_heat(nzlo:nz) = s% extra_power_source
            s% d_extra_heat_dlnd(nzlo:nz) = 0
            s% d_extra_heat_dlnT(nzlo:nz) = 0         
            if (s% use_other_energy) then
               call s% other_energy(s% id, ierr)
               if (ierr /= 0) then
                  if (s% report_ierr .or. dbg) &
                     write(*, *) 'set_vars: other_energy returned ierr', ierr
                  return
               end if
            else 
            
               if (s% inject_uniform_extra_heat /= 0) then 
                  qp1 = 0d0
                  qmin = s% min_q_for_uniform_extra_heat
                  qmax = s% max_q_for_uniform_extra_heat
                  extra = s% inject_uniform_extra_heat
                  do k=nz,nzlo,-1
                     q00 = s% q(k)
                     if (qp1 >= qmin .and. q00 <= qmax) then ! all inside of region
                        s% extra_heat(k) = s% extra_heat(k) + extra
                     else
                        qtop = min(q00, qmax)
                        qbot = max(qp1, qmin)
                        if (qtop > qbot) then ! overlaps region
                           s% extra_heat(k) = s% extra_heat(k) + extra*(qtop - qbot)/s% dq(k)
                        end if
                     end if
                     qp1 = q00
                  end do
               end if
               
               start_time = s% start_time_for_inject_extra_ergs_sec
               if (s% duration_for_inject_extra_ergs_sec > 0) then
                  end_time = start_time + s% duration_for_inject_extra_ergs_sec
               else
                  end_time = 1d99
               end if
               if (s% inject_extra_ergs_sec /= 0 .and. dt > 0 .and. &
                     s% total_mass_for_inject_extra_ergs_sec > 0 .and. &
                     s% time >= start_time .and. &
                     s% time_old < end_time .and. &
                     s% total_energy_start - s% total_energy_initial < &
                        s% inject_until_reach_delta_total_energy) then 
                  qp1 = 0d0
                  qmin = max(0d0, Msun*s% base_of_inject_extra_ergs_sec - s% M_center)/s% xmstar
                  qmax = min(1d0, qmin + Msun*s% total_mass_for_inject_extra_ergs_sec/s% xmstar)
                  extra = s% inject_extra_ergs_sec/(s% xmstar*(qmax - qmin))
                  if (s% time > end_time .or. s% time_old < start_time) &
                     extra = extra*(min(s% time, end_time) - max(s% time_old, start_time))/dt
                  do k=nz,nzlo,-1
                     q00 = s% q(k)
                     if (qp1 >= qmin .and. q00 <= qmax) then ! all inside of region
                        s% extra_heat(k) = s% extra_heat(k) + extra
                     else
                        qtop = min(q00, qmax)
                        qbot = max(qp1, qmin)
                        if (qtop > qbot) then ! overlaps region
                           s% extra_heat(k) = s% extra_heat(k) + extra*(qtop - qbot)/s% dq(k)
                        end if
                     end if
                     qp1 = q00
                  end do
               end if
               
            end if
            
         end if
         
         if (s% irradiation_flux /= 0) then
            do k=nzlo,nz
               s% irradiation_heat(k) = eval_irradiation_heat(s,k)
            end do
         else
            s% irradiation_heat(nzlo:nz) = 0
         end if

      end subroutine update_vars

         
      subroutine set_adjust_mlt_gradT_fraction(s,ierr)
         use eos_def
         type (star_info), pointer :: s 
         integer, intent(out) :: ierr   
         integer :: k     
         include 'formats'
         ierr = 0
         do k=2,s% nz
            if (abs(s% grada_at_face(k) - s% gradr(k)) < 1d-50) then
               s% adjust_mlt_gradT_fraction(k) = 0d0
            else
               s% adjust_mlt_gradT_fraction(k) = &
                  (s% gradT(k) - s% gradr(k))/ &
                     (s% grada_at_face(k) - s% gradr(k))
            end if
         end do
         s% adjust_mlt_gradT_fraction(1) = 0d0 ! not used
      end subroutine set_adjust_mlt_gradT_fraction
      
      
      subroutine set_Teff(s, ierr)
         use atm_lib, only: atm_option
         use star_utils, only: get_r_phot
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: r_phot, L_surf
         logical, parameter :: skip_partials = .true.
         real(dp) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         ierr = 0         
         if (s% use_hydrodynamic_surface_BCs) then
            r_phot = get_r_phot(s)
            L_surf = s% L(1)
            Teff = pow_cr(L_surf/(4*pi*r_phot*r_phot*boltz_sigma), 0.25d0)
            s% Teff = Teff
            s% L_phot = L_surf/Lsun
            s% photosphere_L = s% L_phot
            s% photosphere_r = r_phot/Rsun
            return
         end if
         call set_Teff_info_for_eqns(s, skip_partials, r_phot, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
      end subroutine set_Teff
      
      
      subroutine set_Teff_info_for_eqns(s, skip_partials, r_phot, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
         use atm_lib, only: atm_option
         use star_utils, only: get_r_phot
         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: r_phot, L_surf, Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

         integer :: which_atm_option, off_table_option
         real(dp) :: Tm1, T00, T4_m1, T4_00, P_rad_m1, P_rad_00, &
            alfa, beta, A, opacity_face
         
         include 'formats'
         
         r_phot = get_r_phot(s)
         if (s% L_flag) then
            L_surf = s% L(1)
         else ! assume L(1) = L(2), and L(2) = L_rad(2)
            ! get L_rad(2) from temperature gradient and opacity
            ! L_rad = -d_P_rad/dm_bar*clight*area^2/<opacity_at_face>
            ! d_P_rad = P_rad(1) - P_rad(2)
            ! use values at start of step to avoid having partials
            Tm1 = s% T_start(1)
            T00 = s% T_start(2)
            T4_m1 = Tm1*Tm1*Tm1*Tm1
            T4_00 = T00*T00*T00*T00
            P_rad_m1 = (crad/3)*T4_m1
            P_rad_00 = (crad/3)*T4_00
            alfa = s% dq(1)/(s% dq(1) + s% dq(2))
            beta = 1d0 - alfa
            opacity_face = alfa*s% opacity_start(2) + beta*s% opacity_start(1)               
            A = 4*pi*s% r_start(2)*s% r_start(2)
            L_surf = -(P_rad_m1 - P_rad_00)/s% dm_bar(2)*clight*A*A/opacity_face
         end if
         
         if (L_surf <= 0) then
            if (s% report_ierr) &
               write(*,2) 'L_surf <= 0', s% model_number, L_surf
            ierr = -1
            return
         end if
         
         which_atm_option = atm_option(s% which_atm_option, ierr)
         if (ierr /= 0) then
            write(*,*) 'unknown value for which_atm_option ' // trim(s% which_atm_option)
            return
         end if
      
         off_table_option = atm_option(s% which_atm_off_table_option, ierr)  
         if (ierr /= 0) then
            write(*,*) 'unknown value for which_atm_off_table_option ' // &
               trim(s% which_atm_off_table_option)
            return
         end if
         
         if (s% use_other_surface_PT) then
            call s% other_surface_PT( &
               s% id, which_atm_option, off_table_option, r_phot, L_surf, skip_partials, &
               Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         else if (s% using_free_fall_surface_PT) then
            call get_free_fall_surf_PT( &
               s, which_atm_option, off_table_option, r_phot, L_surf, skip_partials, &
               Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         else
            call get_surf_PT( &
               s, which_atm_option, off_table_option, r_phot, L_surf, skip_partials, &
               Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
         end if
         
         if (ierr /= 0) then
            return
         end if

         s% Teff = Teff
         s% L_phot = L_surf/Lsun
         s% photosphere_L = s% L_phot
         s% photosphere_r = r_phot/Rsun
         
      end subroutine set_Teff_info_for_eqns
      
      
      subroutine set_hydro_vars( &
            s, nzlo, nzhi, skip_basic_vars, skip_micro_vars, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_grads, skip_rotation, &
            skip_brunt, skip_other_cgrav, skip_mixing_info, skip_mlt, ierr)
         use micro, only: set_micro_vars
         use mlt_info, only: set_mlt_vars, set_grads
         use star_utils, only: set_k_CpTMdot_lt_L, start_time, update_time, &
            set_m_grav_and_grav, set_scale_height, get_tau
         use hydro_rotation, only: set_rotation_info
         use brunt, only: do_brunt_N2
         use mix_info, only: set_mixing_info

         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: &
            skip_basic_vars, skip_micro_vars, skip_m_grav_and_grav, &
            skip_net, skip_neu, skip_kap, skip_brunt, &
            skip_grads, skip_rotation, skip_other_cgrav, &
            skip_mixing_info, skip_mlt
         integer, intent(out) :: ierr
         
         integer :: nz, num_nse, k, time0
         logical :: do_Qvisc
         logical, parameter :: skip_eos = .false.
         real(dp) :: total
         
         include 'formats'

         if (dbg) write(*, *) 'set_hydro_vars', nzlo, nzhi
         if (s% doing_timing) call start_time(s, time0, total)

         ierr = 0
         nz = s% nz
         
         if (.not. skip_basic_vars) then
            if (dbg) write(*, *) 'call set_basic_vars'
            call set_basic_vars( &
               s, nzlo, nzhi, skip_other_cgrav, skip_rotation, ierr)
            if (failed('set_basic_vars')) return
         end if

         if (.not. skip_micro_vars) then
            if (dbg) write(*, *) 'call set_micro_vars'
            call set_micro_vars( &
               s, nzlo, nzhi, skip_eos, skip_net, skip_neu, skip_kap, ierr)
            if (failed('set_micro_vars')) return
         end if
         
         call get_tau(s, s% tau)
         
         if (.not. skip_m_grav_and_grav) then 
            ! don't change m_grav or grav during newton iteratons
            if (dbg) write(*, *) 'call set_m_grav_and_grav'
            call set_m_grav_and_grav(s)
         end if
         
         if (dbg) write(*, *) 'call set_scale_height'
         call set_scale_height(s)
         
         if (dbg) write(*, *) 'call set_k_CpTMdot_lt_L'
         call set_k_CpTMdot_lt_L(s)

         if (s% rotation_flag .and. .not. skip_rotation) then
            if (dbg) write(*, *) 'call set_rotation_info'
            call set_rotation_info(s, ierr)
            if (failed('set_rotation_info')) return
         end if
         
         do_Qvisc = (s% v_flag .and. s% use_artificial_viscosity .and. &
            .not. (s% hold_artificial_viscosity_const_during_step &
                        .and. skip_mixing_info))
         if (do_Qvisc) then
            if (dbg) write(*, *) 'call set_Qvisc'
            call set_Qvisc(s, nzlo, nzhi, &
               s% hold_artificial_viscosity_const_during_step, ierr)
            if (failed('set_Qvisc')) return
         end if
         
         if (.not. skip_brunt) then
            if (dbg) write(*,*) 'call do_brunt_N2'
            call do_brunt_N2(s, nzlo, nzhi, ierr)
            if (failed('do_brunt_N2')) return
         end if

         if (.not. skip_grads) then            
            if (dbg) write(*, *) 'call set_grads'
            call set_grads(s, ierr)
            if (failed('set_grads')) return
         end if
         
         if (s% L_flag) then ! only do mixing when have luminosity variables
         
            if (.not. skip_mixing_info) then
               if (dbg) write(*, *) 'call other_adjust_mlt_gradT_fraction'
               call s% other_adjust_mlt_gradT_fraction(s% id,ierr)
               if (failed('other_adjust_mlt_gradT_fraction')) return
               if (dbg) write(*,*) 'call other_after_set_mixing_info'
               call s% other_after_set_mixing_info(s% id, ierr)
               if (failed('other_after_set_mixing_info')) return
            end if
         
            if (.not. skip_mlt) then 
               if (dbg) write(*, *) 'call set_mlt_vars'
               call set_mlt_vars(s, nzlo, nzhi, ierr)
               if (failed('set_mlt_vars')) return
            end if
         
            if (.not. skip_mixing_info) then         
               if (skip_net) then
                  write(*,*) 'ERROR: calling set_mixing_info without doing net 1st'
                  stop 1
               end if           
               if (dbg) write(*, *) 'call set_mixing_info'
               call set_mixing_info(s, ierr)
               if (ierr /= 0) return
               if (s% op_split_gradT_fraction_grada) then
                  call set_adjust_mlt_gradT_fraction(s, ierr)
                  if (ierr /= 0) return
               end if           
            end if
         
         end if
         
         if (s% doing_timing) &
            call update_time(s, time0, total, s% time_set_hydro_vars)         

         
         contains
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            if (ierr == 0) then
               failed = .false.
               return
            end if
            if (s% report_ierr .or. dbg) &
               write(*,*) 'set_hydro_vars failed in call to ' // trim(str)
            failed = .true.
         end function failed

      end subroutine set_hydro_vars
      
      
      subroutine check_rs(s, ierr)
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         integer :: k
         logical :: okay
         include 'formats'
         ierr = 0
         okay = .true.
         do k=2, s% nz
            if (s% r(k) > s% r(k-1)) then
               if (s% report_ierr) then
                  write(*,2) 's% r(k) > s% r(k-1)', k, &
                     s% r(k)/Rsun, s% r(k-1)/Rsun, s% r(k)/Rsun-s% r(k-1)/Rsun
               end if
               okay = .false.
            end if
         end do
         if (okay) return
         ierr = -1
         if (s% report_ierr) write(*,*)
      end subroutine check_rs
      
      
      subroutine set_basic_vars( &
            s, nzlo, nzhi, skip_other_cgrav, skip_rotation, ierr)
         use star_utils, only: set_rv_info
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: skip_other_cgrav, skip_rotation
         integer, intent(out) :: ierr
         integer :: j, k, species
         real(dp) :: twoGmrc2, r2, alfa, beta, sum_xa, v
         
         include 'formats'
         
         if (dbg) write(*,4) 'enter set_basic_vars: nzlo, nzhi, nz', nzlo, nzhi, s% nz
         ierr = 0
         species = s% species
         s% L_phot = s% L(1)/Lsun
         if (.not. s% use_energy_conservation_form) then
            s% d_vc_dv = 1d0
         else
            s% d_vc_dv = 0.5d0
         end if
         
!$OMP PARALLEL DO PRIVATE(j,k,twoGmrc2,sum_xa,v)
         do k=nzlo, nzhi
            if (s% lnPgas_flag) then
               s% Pgas(k) = exp_cr(s% lnPgas(k))
            else
               s% rho(k) = exp_cr(s% lnd(k))
            end if
            if (s% E_flag) then
               s% lnE(k) = log_cr(s% energy(k))
            else
               s% T(k) = exp_cr(s% lnT(k))
               if (s% T_start(k) < 0) s% T_start(k) = s% T(k)
            end if
            if (s% v_start(k) < -1d90) s% v_start(k) = s% v(k) ! v set above
            s% r(k) = exp_cr(s% lnR(k))
            if (s% r_start(k) < 0) s% r_start(k) = s% r(k)
            call set_rv_info(s,k)
            do j=1,species 
               s% xa(j,k) = max(0d0, min(1d0, s% xa(j,k)))
            end do
            sum_xa = sum(s% xa(1:species,k))
            if (abs(sum_xa - 1d0) > 1d-12) then
               do j=1,species
                  s% xa(j,k) = s% xa(j,k)/sum_xa
               end do
            end if
         end do    
!$OMP END PARALLEL DO

         call set_rmid_and_Amid(s, nzlo, nzhi, ierr)
         
         if (.not. skip_other_cgrav) then
            if (s% use_other_cgrav) then
               call s% other_cgrav(s% id, ierr)
               if (ierr /= 0) then
                  if (s% report_ierr .or. dbg) &
                     write(*, *) 'other_cgrav returned ierr', ierr
                  return
               end if
            else
               s% cgrav(nzlo:nzhi) = standard_cgrav
            end if
         end if
         
      end subroutine set_basic_vars

      
      subroutine set_rmid_and_Amid(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         
         integer :: k, nz         
         include 'formats'         
         ierr = 0
         nz = s% nz
!$OMP PARALLEL DO PRIVATE(k)
         do k=nzlo, nzhi
            call do1(k)
         end do
!$OMP END PARALLEL DO
         
         contains
         
         subroutine do1(k)
            integer, intent(in) :: k
            real(dp) :: r003, rp13, rmid, rmid2
            r003 = s% r(k)*s% r(k)*s% r(k)
            if (k < nz) then
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
            else
               rp13 = s% R_center*s% R_center*s% R_center
            end if
            rmid = pow_cr(0.5d0*(r003 + rp13),1d0/3d0)
            s% rmid(k) = rmid
            rmid2 = rmid*rmid
            s% drmid_dlnR00(k) = 0.5d0*r003/rmid2
            s% drmid_dlnRp1(k) = 0.5d0*rp13/rmid2
            s% drmid2_dlnR00(k) = r003/rmid
            s% drmid2_dlnRp1(k) = rp13/rmid
            s% Amid(k) = 4*pi*rmid2
         end subroutine do1         

      end subroutine set_rmid_and_Amid
      
      
      subroutine set_Qvisc(s, nzlo, nzhi, skip_partials, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         
         integer :: k, nz, op_err
         
         include 'formats'

         ierr = 0   
         nz = s% nz
!$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=nzlo, nzhi
            call do1_Qvisc(s,k,skip_partials,nz,op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO

         
      end subroutine set_Qvisc

         
      subroutine zero_Qvisc(s,k)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         s% eta_visc(k) = 0
         s% d_eta_visc_dlnd(k) = 0
         s% d_eta_visc_dlnR00(k) = 0
         s% d_eta_visc_dvel00(k) = 0
         s% d_eta_visc_dlnRp1(k) = 0
         s% d_eta_visc_dvelp1(k) = 0
         s% Qvisc(k) = 0
         s% dQvisc_dlnd(k) = 0
         s% dQvisc_dlnR00(k) = 0
         s% dQvisc_dvel00(k) = 0
         s% dQvisc_dlnRp1(k) = 0
         s% dQvisc_dvelp1(k) = 0
         s% eps_visc(k) = 0
         s% d_eps_visc_dlnd(k) = 0
         s% d_eps_visc_dlnR00(k) = 0
         s% d_eps_visc_dvel00(k) = 0
         s% d_eps_visc_dlnRp1(k) = 0
         s% d_eps_visc_dvelp1(k) = 0
      end subroutine zero_Qvisc


      
      ! see Weaver, Zimmerman, and Woosley, ApJ, 225: 1021-1029, 1978.

      subroutine do1_Qvisc(s,k,skip_partials,nz,ierr) ! set Qvisc and eps_visc
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nz
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         
         real(dp) :: A00, Ap1, v00, vp1, r00, rp1, &
            l1, l22, tmp, d_vc_dv, cs_start, dAv, del, d_A00_dlnR00, &
            d_Ap1_dlnRp1, d_dAv_dlnR00, d_dAv_dlnRp1, &
            d_dAv_dvel00, d_dAv_dvelp1, dm, dm3, rho, Amid, &
            rho2, rmid, rmid2, rmid3, rmid4, rmid6, &
            d_rmid_dlnR00, d_rmid_dlnRp1, inv_rmid3, &
            d_rmid2_dlnR00, d_rmid2_dlnRp1, &
            d_rmid3_dlnR00, d_rmid3_dlnRp1, &
            d_rmid4_dlnR00, d_rmid4_dlnRp1, &
            d_inv_rmid3_dlnR00, d_inv_rmid3_dlnRp1, &
            d_Amid_dlnR00, d_Amid_dlnRp1, eta_visc, &
            d_drdivr_dlnR00, d_drdivr_dlnRp1, &
            dl1_dlnR00, dl1_dlnRp1, &
            dl22_dlnR00, dl22_dlnRp1, &
            Qvisc, &
            d_del_dlnRp1, d_del_dvel00, d_del_dvelp1, &
            d_d_v_div_r_dm_dlnR00, d_d_v_div_r_dm_dlnRp1, &
            d_d_v_div_r_dm_dvel00, d_d_v_div_r_dm_dvelp1, &
            d_d_v_div_r_dr_dlnd, d_d_v_div_r_dr_dlnR00, &
            d_d_v_div_r_dr_dlnRp1, d_d_v_div_r_dr_dvel00, Hq, fac, &
            d_d_v_div_r_dr_dvelp1, d_v_div_r_dm, d_v_div_r_dr, &
            dist_to_inner, dist_to_outer, v_div_cs_00, v_div_cs_m1, &
            v_div_cs_max, v_div_cs_min, v_div_cs_p1
         integer :: kk, k_inner, k_outer, k_lo, k_hi
         
         include 'formats'
         
         ierr = 0
         if (k == 1 .or. &
               s% lnd_start(k)/ln10 <= &
               s% art_visc_full_off_logRho_le_this) then
            call zero_Qvisc(s,k)
            return
         end if
         
         dm = s% dm(k)
         dm3 = dm*dm*dm
         rho = s% rho(k)
         rho2 = rho*rho
         rmid = s% rmid(k)
         cs_start = s% csound_start(k)
         
         r00 = s% r(k)
         if (k == s% nz) then
            rp1 = s% R_center
            if (rp1 <= 0) then
               call zero_Qvisc(s,k)
               return
            end if
            vp1 = s% v_center
         else
            rp1 = s% r(k+1)
            vp1 = s% vc(k+1)
         end if
         v00 = s% vc(k)

         A00 = 4*pi*r00*r00
         Ap1 = 4*pi*rp1*rp1
         
         dAv = Ap1*vp1 - A00*v00
         if (dAv < 0d0) dAv = 0d0
         
         rmid2 = rmid*rmid
         rmid3 = rmid2*rmid
         rmid4 = rmid3*rmid
         rmid6 = rmid4*rmid2
         
         inv_rmid3 = 1d0/rmid3
         
         Amid = 4*pi*rmid2
         
         l1 = s% shock_spread_linear
         l22 = s% shock_spread_quadratic*s% shock_spread_quadratic
         
         if (.false. .and. s% ebdf_stage > 1) then
            eta_visc = s% eta_visc(k)
         else
            eta_visc = 0.75d0*(l1*rmid*rho*cs_start + l22*rmid2*rho2*dAv/dm) ! [g cm^-1 s^-1]
            s% eta_visc(k) = eta_visc
         end if

         if (.false. .and. s% ebdf_stage > 1) then
            d_v_div_r_dm = s% d_v_div_r_dm(k)
            d_v_div_r_dr = s% d_v_div_r_dr(k)
         else
            d_v_div_r_dm = (v00/r00 - vp1/rp1)/dm ! [g^-1 s^-1]            
            d_v_div_r_dr = 4*pi*rmid2*rho*d_v_div_r_dm ! [cm^-1 s^-1]
            s% d_v_div_r_dm(k) = d_v_div_r_dm
            s% d_v_div_r_dr(k) = d_v_div_r_dr
         end if

         fac = 4d0/3d0
         dist_to_inner = 1d99
         dist_to_outer = 1d99
         k_inner = 0
         k_outer = 0
         if (s% pre_shock_viscosity_decay_factor > 0) then      
            ! search inward for next Mach 1 location
            v_div_cs_00 = s% v_start(k)/s% csound_start(k)
            do kk = k+1,nz
               v_div_cs_m1 = v_div_cs_00
               v_div_cs_00 = s% v_start(kk)/s% csound_start(kk)
               v_div_cs_max = max(v_div_cs_00, v_div_cs_m1)
               v_div_cs_min = min(v_div_cs_00, v_div_cs_m1)
               if ((v_div_cs_max >= 1d0 .and. v_div_cs_min < 1d0) .or. &
                   (v_div_cs_min <= -1d0 .and. v_div_cs_max > -1d0)) then
                  dist_to_inner = s% r_start(k) - s% r_start(kk)
                  k_inner = kk
                  exit
               end if
            end do     
         end if
         if (s% post_shock_viscosity_decay_factor > 0) then      
            ! search outward for next Mach 1 location
            v_div_cs_00 = s% v(k)/s% csound(k)
            do kk = k-1,1,-1
               v_div_cs_p1 = v_div_cs_00
               v_div_cs_00 = s% v(kk)/s% csound(kk)
               v_div_cs_max = max(v_div_cs_00, v_div_cs_p1)
               v_div_cs_min = min(v_div_cs_00, v_div_cs_p1)
               if ((v_div_cs_max >= 1d0 .and. v_div_cs_min < 1d0) .or. &
                   (v_div_cs_min <= -1d0 .and. v_div_cs_max > -1d0)) then
                  dist_to_outer = s% r_start(kk) - s% r_start(k)
                  k_outer = kk
                  exit
               end if
            end do
         end if
         if (k_inner > 0 .and. dist_to_inner < dist_to_outer) then ! pre shock
            Hq = s% r_start(k)*s% pre_shock_viscosity_decay_factor
            fac = fac*exp_cr(-dist_to_inner/Hq)
         else if (k_outer > 0) then ! post shock
            Hq = s% r_start(k)*s% post_shock_viscosity_decay_factor
            fac = fac*exp_cr(-(dist_to_outer/Hq)*(dist_to_outer/Hq))
         end if
         
         Qvisc = fac*eta_visc*rmid4*d_v_div_r_dr ! [g cm^2 s^-2] = [ergs]            
         Qvisc = Qvisc + max(0d0, s% artificial_viscosity_Q_shift)
         if (Qvisc >= 0) then
            call zero_Qvisc(s,k)
            return
         end if
         if (s% lnd_start(k)/ln10 < s% art_visc_full_on_logRho_ge_this) &
            Qvisc = Qvisc * &
               (s% lnd_start(k)/ln10 - s% art_visc_full_off_logRho_le_this) / &
               (s% art_visc_full_on_logRho_ge_this - &
                  s% art_visc_full_off_logRho_le_this)            
         
         s% Qvisc(k) = Qvisc
         s% eps_visc(k) = 4*pi*Qvisc*d_v_div_r_dm ! [ergs g^-1 s^-1]
         if (k == -1) then
            write(*,2) 'Qvisc', k, Qvisc
            write(*,2) 'rho', k, rho
            write(*,2) 'fac', k, fac
            write(*,2) 'eta_visc', k, eta_visc
            write(*,2) 'rmid4', k, rmid4
            write(*,2) 'd_v_div_r_dr', k, d_v_div_r_dr
            write(*,*)
         end if

         if (is_bad_num(s% eps_visc(k))) then
            ierr = -1
            return
!$OMP critical
            write(*,2) 's% eps_visc(k)', k, s% eps_visc(k)
            write(*,2) 'Qvisc', k, Qvisc
            write(*,2) 'd_v_div_r_dm', k, d_v_div_r_dm
            write(*,2) 'fac', k, fac
            write(*,2) 'eta_visc', k, eta_visc
            write(*,2) 'rmid4', k, rmid4
            write(*,2) 'd_v_div_r_dr', k, d_v_div_r_dr
            write(*,2) 'cs_start', k, cs_start
            write(*,2) 'dAv', k, dAv
            write(*,2) 'l22', k, l22
            write(*,2) 'l1', k, l1
            write(*,2) 'rmid', k, rmid
            write(*,2) 'rho2', k, rho2
            write(*,2) 'dm', k, dm
            stop 'Qvisc'
!$OMP end critical
         end if
      
         s% d_eta_visc_dlnd(k) = 0
         s% d_eta_visc_dlnR00(k) = 0
         s% d_eta_visc_dlnRp1(k) = 0
         s% d_eta_visc_dvel00(k) = 0
         s% d_eta_visc_dvelp1(k) = 0

         s% dQvisc_dlnd(k) = 0
         s% dQvisc_dlnR00(k) = 0
         s% dQvisc_dlnRp1(k) = 0
         s% dQvisc_dvel00(k) = 0
         s% dQvisc_dvelp1(k) = 0
      
         s% d_eps_visc_dlnd(k) = 0
         s% d_eps_visc_dlnR00(k) = 0
         s% d_eps_visc_dlnRp1(k) = 0
         s% d_eps_visc_dvel00(k) = 0
         s% d_eps_visc_dvelp1(k) = 0
         
         if (skip_partials) return
         
         d_rmid_dlnR00 = s% drmid_dlnR00(k)
         d_rmid_dlnRp1 = s% drmid_dlnRp1(k)
         d_vc_dv = s% d_vc_dv

         d_A00_dlnR00 = 2*A00
         d_Ap1_dlnRp1 = 2*Ap1
         
         if (dAv <= 0d0) then
            d_dAv_dlnR00 = 0
            d_dAv_dlnRp1 = 0
            d_dAv_dvel00 = 0
            d_dAv_dvelp1 = 0
         else
            d_dAv_dlnR00 = -d_A00_dlnR00*v00
            d_dAv_dlnRp1 = d_Ap1_dlnRp1*vp1
            d_dAv_dvel00 = -d_vc_dv*A00
            d_dAv_dvelp1 = d_vc_dv*Ap1
         end if
                     
         d_rmid2_dlnR00 = 2*rmid*d_rmid_dlnR00
         d_rmid2_dlnRp1 = 2*rmid*d_rmid_dlnRp1
         d_rmid3_dlnR00 = 3*rmid2*d_rmid_dlnR00
         d_rmid3_dlnRp1 = 3*rmid2*d_rmid_dlnRp1
         d_rmid4_dlnR00 = 4*rmid3*d_rmid_dlnR00
         d_rmid4_dlnRp1 = 4*rmid3*d_rmid_dlnRp1
         
         d_inv_rmid3_dlnR00 = -3*d_rmid_dlnR00/rmid4
         d_inv_rmid3_dlnRp1 = -3*d_rmid_dlnRp1/rmid4
         
         d_Amid_dlnR00 = 4*pi*d_rmid2_dlnR00
         d_Amid_dlnRp1 = 4*pi*d_rmid2_dlnRp1

         if (s% ebdf_stage > 1 .and. .false.) then
            s% d_eta_visc_dlnd(k) = 0
            s% d_eta_visc_dvel00(k) = 0
            s% d_eta_visc_dvelp1(k) = 0
            s% d_eta_visc_dlnR00(k) = 0
            s% d_eta_visc_dlnRp1(k) = 0
         else
            s% d_eta_visc_dlnd(k) = 0.75d0*(l1*rmid*rho*cs_start + 2*l22*rmid2*rho2*dAv/dm)
            s% d_eta_visc_dvel00(k) = 0.75d0*l22*rmid2*rho2*d_dAv_dvel00/dm
            s% d_eta_visc_dvelp1(k) = 0.75d0*l22*rmid2*rho2*d_dAv_dvelp1/dm
            s% d_eta_visc_dlnR00(k) = 0.75d0*( &
               l1*d_rmid_dlnR00*rho*cs_start + &
               l22*rho2/dm*(d_rmid2_dlnR00*dAv + rmid2*d_dAv_dlnR00))
            s% d_eta_visc_dlnRp1(k) = 0.75d0*( &
               l1*d_rmid_dlnRp1*rho*cs_start + &
               l22*rho2/dm*(d_rmid2_dlnRp1*dAv + rmid2*d_dAv_dlnRp1))
         end if
         
         if (.false. .and. s% ebdf_stage > 1) then
            d_d_v_div_r_dm_dlnR00 = 0
            d_d_v_div_r_dm_dlnRp1 = 0
            d_d_v_div_r_dm_dvel00 = 0
            d_d_v_div_r_dm_dvelp1 = 0
            d_d_v_div_r_dr_dlnd = 0
            d_d_v_div_r_dr_dlnR00 = 0
            d_d_v_div_r_dr_dlnRp1 = 0
            d_d_v_div_r_dr_dvel00 = 0
            d_d_v_div_r_dr_dvelp1 = 0
         else
         
            d_d_v_div_r_dm_dlnR00 = -v00/r00/dm
            d_d_v_div_r_dm_dlnRp1 = vp1/rp1/dm
            d_d_v_div_r_dm_dvel00 = d_vc_dv/r00/dm
            d_d_v_div_r_dm_dvelp1 = -d_vc_dv/rp1/dm
            
            d_d_v_div_r_dr_dlnd = d_v_div_r_dr
            d_d_v_div_r_dr_dlnR00 = 4*pi*rho*( &
               d_rmid2_dlnR00*d_v_div_r_dm + rmid2*d_d_v_div_r_dm_dlnR00)
            d_d_v_div_r_dr_dlnRp1 = 4*pi*rho*( &
               d_rmid2_dlnRp1*d_v_div_r_dm + rmid2*d_d_v_div_r_dm_dlnRp1) 
            d_d_v_div_r_dr_dvel00 = 4*pi*rmid2*rho*d_d_v_div_r_dm_dvel00
            d_d_v_div_r_dr_dvelp1 = 4*pi*rmid2*rho*d_d_v_div_r_dm_dvelp1
            
         end if
         
         if (.false. .and. s% ebdf_stage > 1) then
         
            s% dQvisc_dlnd(k) = 0
            s% dQvisc_dlnR00(k) = 0
            s% dQvisc_dlnRp1(k) = 0
            s% dQvisc_dvel00(k) = 0
            s% dQvisc_dvelp1(k) = 0
         
            s% d_eps_visc_dlnd(k) = 0
            s% d_eps_visc_dlnR00(k) = 0
            s% d_eps_visc_dlnRp1(k) = 0
            s% d_eps_visc_dvel00(k) = 0
            s% d_eps_visc_dvelp1(k) = 0
         else
         
            s% dQvisc_dlnd(k) = fac*rmid4*( &
               s% d_eta_visc_dlnd(k)*d_v_div_r_dr + &
               eta_visc*d_d_v_div_r_dr_dlnd)
            
            if (k == -1) then
               write(*,2) 's% dQvisc_dlnd(k)', k, s% dQvisc_dlnd(k)
               write(*,2) 'eta_visc', k, eta_visc
               write(*,2) 's% d_eta_visc_dlnd(k)', k, s% d_eta_visc_dlnd(k)
               write(*,2) 'd_v_div_r_dr', k, d_v_div_r_dr
               write(*,2) 'd_d_v_div_r_dr_dlnd', k, d_d_v_div_r_dr_dlnd
               write(*,*)
            end if

            s% dQvisc_dlnR00(k) = fac*( &
               s% d_eta_visc_dlnR00(k)*rmid4*d_v_div_r_dr + &
               eta_visc*d_rmid4_dlnR00*d_v_div_r_dr + &
               eta_visc*rmid4*d_d_v_div_r_dr_dlnR00)
            s% dQvisc_dlnRp1(k) = fac*( &
               s% d_eta_visc_dlnRp1(k)*rmid4*d_v_div_r_dr + &
               eta_visc*d_rmid4_dlnRp1*d_v_div_r_dr + &
               eta_visc*rmid4*d_d_v_div_r_dr_dlnRp1) 
            s% dQvisc_dvel00(k) = fac*rmid4*( &
               s% d_eta_visc_dvel00(k)*d_v_div_r_dr + eta_visc*d_d_v_div_r_dr_dvel00) 
            s% dQvisc_dvelp1(k) = fac*rmid4*( &
               s% d_eta_visc_dvelp1(k)*d_v_div_r_dr + eta_visc*d_d_v_div_r_dr_dvelp1)  
         
            s% d_eps_visc_dlnd(k) = 4*pi*s% dQvisc_dlnd(k)*d_v_div_r_dm
            s% d_eps_visc_dlnR00(k) = 4*pi*( &
               s% dQvisc_dlnR00(k)*d_v_div_r_dm + Qvisc*d_d_v_div_r_dm_dlnR00)
            s% d_eps_visc_dlnRp1(k) = 4*pi*( &
               s% dQvisc_dlnRp1(k)*d_v_div_r_dm + Qvisc*d_d_v_div_r_dm_dlnRp1) 
            s% d_eps_visc_dvel00(k) = 4*pi*( &
               s% dQvisc_dvel00(k)*d_v_div_r_dm + Qvisc*d_d_v_div_r_dm_dvel00) 
            s% d_eps_visc_dvelp1(k) = 4*pi*( &
               s% dQvisc_dvelp1(k)*d_v_div_r_dm + Qvisc*d_d_v_div_r_dm_dvelp1) 
               
         end if
         

      end subroutine do1_Qvisc
      
      
      real(dp) function get_Teff(s, ierr) result(Teff)
         use atm_lib, only: atm_option
         use star_utils, only: get_r_phot
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: which_atm_option, off_table_option
         real(dp) :: r_phot, L_surf
         logical :: skip_partials
         real(dp) :: &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
            
            
         include 'formats'
            
         ierr = 0
         
         r_phot = get_r_phot(s)
         L_surf = s% L(1)
         skip_partials = .true.
         
         which_atm_option = atm_option(s% which_atm_option, ierr)
         if (ierr /= 0) then
            write(*,*) 'unknown value for which_atm_option ' // trim(s% which_atm_option)
            return
         end if
      
         off_table_option = atm_option(s% which_atm_off_table_option, ierr)  
         if (ierr /= 0) then
            write(*,*) 'unknown value for which_atm_off_table_option ' // &
               trim(s% which_atm_off_table_option)
            return
         end if
         
         call get_surf_PT( &
            s, which_atm_option, off_table_option, r_phot, L_surf, skip_partials, &
            Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
         if (ierr /= 0) then
            return
         end if

      end function get_Teff

      
      subroutine get_surf_PT( &
            s, which_atm_option, off_table_option, r_phot, L_surf, skip_partials, &
            Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
         use atm_lib
         use atm_def
         use chem_def
         use create_atm, only: get_Paczynski_atm_surf_PT
         use opacities, only: fraction_of_op_mono
         use eos_lib, only: Radiation_Pressure

         type (star_info), pointer :: s
         integer, intent(in) :: which_atm_option, off_table_option
         real(dp), intent(in) :: r_phot, L_surf
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr
         
         real(dp) :: opacity, X, Y, Z, M, tau_surface, tau_center_cell1, &
            kap_surf, err, Teff4, T_surf4, T_surf, P_surf_atm, P_surf
         integer :: iters, num_atm_structure_points
         real(dp), pointer :: atm_structure_results(:,:)
         logical :: save_atm_structure_info
         
         include 'formats'
         
         if (s% L(1) <= 0) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) 'get_surf_PT L <= 0', ierr, s% L(1)
               if (s% surface_extra_Pgas > 0) &
                  write(*,1) 'surface_extra_Pgas', s% surface_extra_Pgas
            end if
            return
         end if
      
         opacity = s% opacity(1)
         X = s% X(1); Y = s% Y(1); Z = 1 - (X+Y)
         tau_surface = s% tau_factor*s% tau_base ! tau at outer edge of cell 1
         M = s% m_grav(1)
         
         !write(*,1) 'get_surf_PT tau_surface', tau_surface

         !write(*,*) 's% dtau1_start > tau_surface', s% dtau1_start > tau_surface
         if (s% use_other_atm) then
         
            call s% other_atm( &
               s% id, M, r_phot, L_surf, X, Z, opacity, Teff, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               which_atm_option, s% atm_switch_to_grey_as_backup, ierr)
               
         else if (which_atm_option == atm_grey_and_kap) then
         
            if (fraction_of_op_mono(s,1) > 0d0) then
               call atm_get_grey_and_op_mono_kap( &
                  s% Pextra_factor, tau_surface, opacity, &
                  s% cgrav(1), M, r_phot, L_surf, X, Z, s% abar(1), s% zbar(1), & 
                  s% species, s% chem_id, s% net_iso, s% xa(:,1), &
                  s% atm_grey_and_kap_max_tries, s% atm_grey_and_kap_atol, &
                  s% atm_grey_and_kap_rtol, &
                  s% eos_handle, s% kap_handle, &
                  s% use_op_mono_alt_get_kap, s% op_mono_min_X_to_include, &
                  s% op_mono_data_path, s% op_mono_data_cache_filename, &
                  lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, & 
                  lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, & 
                  kap_surf, Teff, iters, err, ierr)
            else
               call atm_get_grey_and_kap( &
                  s% Pextra_factor, tau_surface, opacity, &
                  s% cgrav(1), M, r_phot, L_surf, X, Z, s% abar(1), s% zbar(1), & 
                  s% species, s% chem_id, s% net_iso, s% xa(:,1), &
                  s% atm_grey_and_kap_max_tries, s% atm_grey_and_kap_atol, &
                  s% atm_grey_and_kap_rtol, &
                  s% eos_handle, s% kap_handle, &
                  lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, & 
                  lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, & 
                  kap_surf, Teff, iters, err, ierr)
            end if
            
            dlnT_dlnkap = 0
            dlnP_dlnkap = 0
            if (s% trace_atm_grey_and_kap .and. .not. skip_partials) then
               write(*,2) 'grey&kap: itrs, err, kap_s, kap1', &
                  iters, err, kap_surf, opacity
            end if    
            
         else if (which_atm_option == atm_Paczynski_grey) then
         
            call get_Paczynski_atm_surf_PT( &
               s, r_phot, L_surf, skip_partials, Teff, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
            if (s% trace_atm_Paczynski_grey .and. .not. skip_partials) then
               write(*,*) 
               write(*,'(99a20)') &
                  'Paczynski_grey', 'Teff', 'lnT', 'dlnT_dL', &
                  'dlnT_dlnR', 'lnP', 'dlnP_dL', 'dlnP_dlnR'
               write(*,'(i20,99e20.10)') &
                  s% model_number, Teff, lnT_surf, dlnT_dL, dlnT_dlnR, &
                  lnP_surf, dlnP_dL, dlnP_dlnR
            end if 
               
         else if (which_atm_option == atm_grey_irradiated) then
         
            if (s% atm_grey_irradiated_simple_kap_th .and. &
                  s% atm_grey_irr_kap_v_div_kap_th > 0) &
               opacity = s% atm_grey_irradiated_kap_v/s% atm_grey_irr_kap_v_div_kap_th
            call atm_grey_irradiated_get( &
               s% atm_grey_irradiated_T_eq, s% atm_grey_irradiated_kap_v, &
               opacity, s% atm_grey_irradiated_simple_kap_th, &
               s% atm_grey_irr_kap_v_div_kap_th, s% atm_grey_irradiated_P_surf, &
               s% cgrav(1), M, r_phot, L_surf, X, Z, s% abar(1), s% zbar(1), & 
               s% species, s% chem_id, s% net_iso, s% xa(:,1), &
               s% atm_grey_irradiated_max_tries, s% atm_grey_irradiated_atol, &
               s% atm_grey_irradiated_rtol, &
               s% eos_handle, s% kap_handle, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, & 
               kap_surf, tau_surface, Teff, iters, err, ierr)
            s% tau_factor = tau_surface/s% tau_base
            lnP_surf = log_cr(s% atm_grey_irradiated_P_surf)
            dlnP_dL = 0
            dlnP_dlnR = 0
            dlnP_dlnM = 0
            dlnT_dlnkap = 0
            dlnP_dlnkap = 0
            if (s% trace_atm_grey_irradiated .and. .not. skip_partials) then
               write(*,2) &
                  'grey_irradiated: itrs, err, kap_s, kap1, tau_surface, Teff, T', &
                  iters, err, kap_surf, opacity, tau_surface, Teff, exp_cr(lnT_surf)
            end if    
            
         else if (atm_by_integration(which_atm_option)) then
         
            save_atm_structure_info = .false.
            nullify(atm_structure_results)
            call atm_get_int_T_tau( &
               s% atm_int_errtol, s% cgrav(1), M, r_phot, L_surf, &
               X, Z, s% abar(1), s% zbar(1), & 
               s% species, s% chem_id, s% net_iso, s% xa(:,1), &
               which_atm_option, s% eos_handle, s% kap_handle, save_atm_structure_info, &
               tau_surface, skip_partials, Teff, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               num_atm_structure_points, atm_structure_results, ierr)
               
         else if (atm_by_table(which_atm_option)) then
         
            if (s% tau_factor /= 1) then
               write(*,*) 'cannot simultaneously have atm by table and tau_factor /= 1'
               ierr = -1
               return
            end if
            s% tau_base = atm_tau_base(which_atm_option, ierr)
            if (ierr == 0) call atm_get_table( &
               which_atm_option, off_table_option, &
               s% cgrav(1), M, r_phot, L_surf, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, & 
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, & 
               s% Pextra_factor, opacity, X, Z, s% abar(1), s% zbar(1), & 
               s% species, s% chem_id, s% net_iso, s% xa(:,1), &
               s% atm_grey_and_kap_max_tries, &
               s% atm_grey_and_kap_atol, s% atm_grey_and_kap_rtol, &
               s% eos_handle, s% kap_handle, Teff, ierr)
            dlnT_dlnkap = 0
            dlnP_dlnkap = 0
               
         else if (which_atm_option == atm_fixed_Teff) then ! use atm_fixed_Teff
            ! set Tsurf from Eddington T-tau relation
            !     for current surface tau and Teff = `atm_fixed_Teff`.
            ! set Psurf = Radiation_Pressure(Tsurf)

            Teff = s% atm_fixed_Teff
            Teff4 = Teff*Teff*Teff*Teff
            T_surf = pow_cr(3d0/4d0*Teff4*(tau_surface + 2d0/3d0), 0.25d0)
            lnT_surf = log_cr(T_surf)
            lnP_surf = Radiation_Pressure(T_surf)
            
            dlnT_dL = 0; dlnT_dlnR = 0; dlnT_dlnM = 0; dlnT_dlnkap = 0
            dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
               
         else if (which_atm_option == atm_fixed_Tsurf) then ! use atm_fixed_Tsurf
            ! set Teff from Eddington T-tau relation for given Tsurf and tau=2/3
            ! set Psurf = Radiation_Pressure(Tsurf)
            
            T_surf = s% atm_fixed_Tsurf
            lnT_surf = log_cr(T_surf)
            T_surf4 = T_surf*T_surf*T_surf*T_surf
            Teff = pow_cr(4d0/3d0*T_surf4/(tau_surface + 2d0/3d0), 0.25d0)
            lnP_surf = Radiation_Pressure(T_surf)
            
            dlnT_dL = 0; dlnT_dlnR = 0; dlnT_dlnM = 0; dlnT_dlnkap = 0
            dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
               
         else if (which_atm_option == atm_fixed_Psurf) then ! use atm_fixed_Psurf
            ! set Tsurf from L and R using L = 4*pi*R^2*boltz_sigma*T^4.
            ! set Teff using Eddington T-tau relation for tau=2/3 and T=Tsurf.
            
            lnP_surf = safe_log_cr(s% atm_fixed_Psurf)
            T_surf4 = L_surf/(4*pi*r_phot*r_phot*boltz_sigma)
            T_surf = pow_cr(T_surf4, 0.25d0)
            lnT_surf = log_cr(T_surf)
            dlnT_dlnR = -0.5d0
            dlnT_dL = 1d0/(4d0*L_surf)
            Teff = pow_cr(4d0/3d0*T_surf4/(tau_surface + 2d0/3d0), 0.25d0)
            
            dlnT_dlnM = 0; dlnT_dlnkap = 0
            dlnP_dL = 0; dlnP_dlnR = 0; dlnP_dlnM = 0; dlnP_dlnkap = 0
            
         else
         
            call atm_get_grey( &
               s% Pextra_factor, tau_surface, s% cgrav(1), &
               M, r_phot, L_surf, opacity, Teff, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
               ierr)
               
         end if
         
         if (s% surface_extra_Pgas /= 0d0) then
            P_surf_atm = exp_cr(lnP_surf)
            P_surf = P_surf_atm + s% surface_extra_Pgas
            if (P_surf < 1d-50) then
               lnP_surf = -50*ln10
               dlnP_dL = 0
               dlnP_dlnR = 0
               dlnP_dlnM = 0
               dlnP_dlnkap = 0
            else
               lnP_surf = log_cr(P_surf)
               dlnP_dL = dlnP_dL*P_surf_atm/P_surf
               dlnP_dlnR = dlnP_dlnR*P_surf_atm/P_surf
               dlnP_dlnM = dlnP_dlnM*P_surf_atm/P_surf
               dlnP_dlnkap = dlnP_dlnkap*P_surf_atm/P_surf
            end if
         end if
         
         !if (s% model_number == 2) then
         !   call show
         !   stop 'get PT surf'
         !end if
         
         if (is_bad_num(lnT_surf) .or. is_bad_num(lnP_surf)) then
            call show
            ierr = -1
            return
            !stop 'get_surf_PT'
         end if

         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*, *) 'atm_get returned ierr', ierr
               write(*, *)
            end if
            return
         end if
         
         !write(*,1) 'get_surf_PT log(Teff) R L', &
         !   log10_cr(Teff), r_phot/Rsun, L_surf/Lsun
         
         contains
         
         subroutine show
            include 'formats'
            write(*,1) 'M =', M
            write(*,1) 'R =', r_phot
            write(*,1) 'L =', L_surf
            write(*,1) 'X =', X
            write(*,1) 'Z =', Z
            write(*,1) 'opacity =', opacity
            write(*,2) 'which_atm_option =', which_atm_option
            write(*,*) trim(s% which_atm_option)
            write(*,*) 's% atm_switch_to_grey_as_backup', s% atm_switch_to_grey_as_backup
            write(*,*)
            write(*,1) 's% tau_factor', s% tau_factor
            write(*,*)
            write(*,1) 'Teff', Teff
            write(*,1) 'logT_surf', lnT_surf/ln10
            write(*,1) 'logP_surf', lnP_surf/ln10
            write(*,*)
            write(*,1) 'lnT_surf', lnT_surf
            write(*,1) 'lnP_surf', lnP_surf
            write(*,*)
            write(*,1) 'dlnT_dL', dlnT_dL
            write(*,1) 'dlnT_dlnR', dlnT_dlnR
            write(*,1) 'dlnT_dlnM', dlnT_dlnM
            write(*,1) 'dlnT_dlnkap', dlnT_dlnkap
            write(*,*)
            write(*,1) 'dlnP_dL', dlnP_dL
            write(*,1) 'dlnP_dlnR', dlnP_dlnR
            write(*,1) 'dlnP_dlnM', dlnP_dlnM
            write(*,1) 'dlnP_dlnkap', dlnP_dlnkap
            write(*,*)
         end subroutine show
         
      end subroutine get_surf_PT
            

      real(dp) function eval_cell_collapse_timescale(s) result(dt)
         type (star_info), pointer :: s
         
         integer :: nz, k
         real(dp) :: dr, dv, cell_dt
         include 'formats'
         
         dt = 1d99
         if (.not. s% v_flag) return
         
         nz = s% nz
         do k=1,nz-1
            dv = s% v(k+1) - s% v(k)
            if (dv <= 0) cycle
            dr = s% r(k) - s% r(k+1)
            cell_dt = dr/dv
            if (cell_dt < dt) dt = cell_dt
         end do
         
      end function eval_cell_collapse_timescale
      

      real(dp) function eval_chem_timescale(s) result(dt)
         type (star_info), pointer :: s

         real(dp), pointer, dimension(:) :: sig
         integer :: nz, k, j, species, min_k, min_j
         real(dp) :: dt0, dm, sig00, sigp1, dxdt_nuc, dxdt_chem, &
            x00, xm1, xp1, dx00, dxp1, flux00, fluxp1, dxdt_mix
         
         include 'formats'
         
         dt = 1d99
         min_k = 0
         min_j = 0
         nz = s% nz
         species = s% species
         sig => s% sig
         
         do k = 1, s% nz
            dm = s% dm(k)
            sig00 = sig(k)   
            if (k < nz) then
               sigp1 = sig(k+1)
            else
               sigp1 = 0
            end if
            do j=1,species
               if (.not. s% do_burn) then
                  dxdt_nuc = 0
               else
                  dxdt_nuc = s% dxdt_nuc(j,k)
               end if           
               x00 = s% xa(j,k)
               if (x00 < 1d-10) cycle
               if (.not. s% do_mix) then
                  dxdt_mix = 0
               else
                  if (k > 1) then
                     xm1 = s% xa(j,k-1)
                     dx00 = xm1 - x00
                     flux00 = -sig00*dx00
                  else
                     flux00 = 0
                  end if         
                  if (k < nz) then
                     xp1 = s% xa(j,k+1)
                     dxp1 = x00 - xp1
                     fluxp1 = -sigp1*dxp1
                  else
                     fluxp1 = 0
                  end if
                  dxdt_mix = (fluxp1 - flux00)/dm
               end if
               dxdt_chem = dxdt_nuc + dxdt_mix
               if (dxdt_chem >= 0) cycle   
               dt0 = max(1d-20, x00)/max(1d-50, abs(dxdt_chem))
               if (dt0 < dt) then
                  dt = dt0
                  min_k = k
                  min_j = j
               end if
            end do
         end do
         
         return
         
         j = min_j
         k = min_k
         write(*,3) 'chem timescale', j, k, dt, s% xa(j,k), s% dxdt_nuc(j,k) 
         
      end function eval_chem_timescale
         
         
      subroutine get_gr_gravity_factor(s, k, gr_factor, d_gr_factor_dlnR00)
         ! note: this uses gravitational mass, m_grav
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(out) :: gr_factor, d_gr_factor_dlnR00            
         real(dp) :: twoGmrc2, invfac, dtwoGmrc2_dlnR, d_invfac_dlnR
         twoGmrc2 = 2*s% cgrav(k)*s% m_grav(k)/(s% r(k)*clight*clight)
         invfac = sqrt(1d0 - twoGmrc2)
         gr_factor = 1d0/invfac         
         dtwoGmrc2_dlnR = -twoGmrc2
         d_invfac_dlnR = -dtwoGmrc2_dlnR/(2*invfac)
         d_gr_factor_dlnR00 = -d_invfac_dlnR/(invfac*invfac)            
      end subroutine get_gr_gravity_factor


      subroutine do_dlnP_dm(s, k, do_visc, delm, &
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
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: do_visc
         real(dp), intent(in) :: delm
         real(dp), intent(out) :: &
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
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas
         integer, intent(out) :: ierr
         
         real(dp) :: &
            m, dv_dt, dlnP_dm, fac, vfac, &
            r, R2, rPterm, accel, d_rPterm_dlnR, dlnP_dm_accel_term, d_accel_dlnR, &
            d_rPterm_dlnd00, d_rPterm_dlnT00, d_rPterm_dlndm1, d_rPterm_dlnTm1, &
            d_rPterm_dlnPgas00_const_T, d_rPterm_dlnT00_const_Pgas, &
            d_rPterm_dlnPgasm1_const_T, d_rPterm_dlnTm1_const_Pgas, &
            d_conv_dP_term_dlnPgas00_const_T, d_conv_dP_term_dlnPgasm1_const_T, &
            d_conv_dP_term_dlnT00_const_Pgas, d_conv_dP_term_dlnTm1_const_Pgas, &
            gr_factor, d_gr_factor_dlnR00, &
            d_accel_dlnR00, d_accel_dlnRp1, d_accel_dlnRm1, &
            dfac_dlnd00, dfac_dlnT00, dfac_dlndm1, dfac_dlnTm1, dfac_dlnR, dfac_dL, &
            dfac_dlnPgas00_const_T, dfac_dlnPgasm1_const_T, &
            dfac_dlnT00_const_Pgas, dfac_dlnTm1_const_Pgas, &
            d_dlnPdm_visc_dlnd00, d_dlnPdm_visc_dlndm1, &
            d_dlnPdm_visc_dlnTm1, d_dlnPdm_visc_dlnT00, &
            d_dlnPdm_visc_dlnRm1, d_dlnPdm_visc_dlnR00, &
            d_dlnPdm_visc_dlnRp1, d_dlnPdm_visc_dvelm1, &
            d_dlnPdm_visc_dvel00, d_dlnPdm_visc_dvelp1, &
            d_dlnPdm_visc_dlnPgas00_const_T, d_dlnPdm_visc_dlnPgasm1_const_T, &
            d_dlnPdm_visc_dlnT00_const_Pgas, d_dlnPdm_visc_dlnTm1_const_Pgas, &
            d_dQ_dlnd00, d_dQ_dlndm1, d_dQ_dlnR00, d_dQ_dlnRm1, d_dQ_dlnRp1, &
            d_dQ_dvel00, d_dQ_dvelm1, d_dQ_dvelp1, visc, d_visc_dlnd00, &
            d_visc_dlndm1, d_visc_dlnR00, d_visc_dlnRm1, &
            d_visc_dlnRp1, d_visc_dvel00, d_visc_dvelm1, &
            d_visc_dvelp1, dlnP_dm_visc, dm_bar, dQ, dvdt_visc, &
            inv_Ppoint, d_inv_Ppoint_dlnd00, d_inv_Ppoint_dlnT00, &
            d_inv_Ppoint_dlndm1, d_inv_Ppoint_dlnTm1, dvfac_dlnd00, &
            dvfac_dlnT00, dvfac_dlndm1, dvfac_dlnTm1, dvfac_dlnR00
                        
         logical :: local_v_flag, dbg

         include 'formats'
         
         ierr = 0
         
         dbg = .false.

         s% dlnP_dm_expected(k) = 0
         dlnP_dm = 0
         d_dlnPdm_dlnRp1 = 0
         d_dlnPdm_dlnR00 = 0
         d_dlnPdm_dlnRm1 = 0
         d_dlnPdm_dvelp1 = 0
         d_dlnPdm_dvel00 = 0
         d_dlnPdm_dvelm1 = 0
         d_dlnPdm_dlnd00 = 0
         d_dlnPdm_dlndm1 = 0
         d_dlnPdm_dlnT00 = 0
         d_dlnPdm_dlnTm1 = 0
         d_dlnPdm_dlnPgas00_const_T = 0
         d_dlnPdm_dlnT00_const_Pgas = 0
         d_dlnPdm_dlnPgasm1_const_T = 0
         d_dlnPdm_dlnTm1_const_Pgas = 0
         d_dlnPdm_dL00 = 0

         call eval_dlnPdm_qhse(s, k, m, &
            dlnP_dm, &
            d_dlnPdm_dlnR00, d_dlnPdm_dL00, &
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
         s% dlnP_dm_expected(k) = dlnP_dm

         d_dlnPdm_dlnRm1 = 0

         if (is_bad_num(s% dlnP_dm_expected(k))) then
            write(*,*) 'k', k
            write(*,1) 'dlnP_dm_expected(k)', s% dlnP_dm_expected(k)
            write(*,*)
            ierr = -1
            return
         end if
         
         r = s% r(k)
         
         if (s% conv_dP_term(k) /= 0 .and.  s% conv_dP_term_factor > 0) then 
            ! include effect of convective turbulence on pressure gradient
            
            fac = 1d0/(1d0 + s% conv_dP_term_factor*s% conv_dP_term(k))

            s% dlnP_dm_expected(k) = s% dlnP_dm_expected(k)*fac
            
            dfac_dlnd00 = -fac*fac*s% conv_dP_term_factor*s% d_conv_dP_term_dlnd00(k)
            dfac_dlnT00 = -fac*fac*s% conv_dP_term_factor*s% d_conv_dP_term_dlnT00(k)
            dfac_dlndm1 = -fac*fac*s% conv_dP_term_factor*s% d_conv_dP_term_dlndm1(k)
            dfac_dlnTm1 = -fac*fac*s% conv_dP_term_factor*s% d_conv_dP_term_dlnTm1(k)
            dfac_dlnR = -fac*fac*s% conv_dP_term_factor*s% d_conv_dP_term_dlnR(k)
            dfac_dL = -fac*fac*s% conv_dP_term_factor*s% d_conv_dP_term_dL(k)
                                
            if (s% lnPgas_flag) then
               d_conv_dP_term_dlnPgas00_const_T = &
                  dfac_dlnd00*s% dlnRho_dlnPgas_const_T(k)
               d_conv_dP_term_dlnT00_const_Pgas = &
                  dfac_dlnT00 + dfac_dlnd00*s% dlnRho_dlnT_const_Pgas(k)
               
               if (k > 1) then
                  d_conv_dP_term_dlnPgasm1_const_T = &
                     dfac_dlndm1*s% dlnRho_dlnPgas_const_T(k-1)
                  d_conv_dP_term_dlnTm1_const_Pgas = &
                     dfac_dlnTm1 + dfac_dlndm1*s% dlnRho_dlnT_const_Pgas(k-1)
               else
                  d_conv_dP_term_dlnPgasm1_const_T = 0
                  d_conv_dP_term_dlnTm1_const_Pgas = 0
               end if
                  
               dfac_dlnPgas00_const_T = &
                  -fac*fac*s% conv_dP_term_factor*d_conv_dP_term_dlnPgas00_const_T
               dfac_dlnPgasm1_const_T = &
                  -fac*fac*s% conv_dP_term_factor*d_conv_dP_term_dlnPgasm1_const_T
               dfac_dlnT00_const_Pgas = &
                  -fac*fac*s% conv_dP_term_factor*d_conv_dP_term_dlnT00_const_Pgas
               dfac_dlnTm1_const_Pgas = &
                  -fac*fac*s% conv_dP_term_factor*d_conv_dP_term_dlnTm1_const_Pgas
                  
               d_dlnPdm_dlnPgas00_const_T = d_dlnPdm_dlnPgas00_const_T*fac + &
                  s% dlnP_dm_expected(k)*dfac_dlnPgas00_const_T
               d_dlnPdm_dlnT00_const_Pgas = d_dlnPdm_dlnT00_const_Pgas*fac + &
                  s% dlnP_dm_expected(k)*dfac_dlnT00_const_Pgas
               d_dlnPdm_dlnPgasm1_const_T = d_dlnPdm_dlnPgasm1_const_T*fac + &
                  s% dlnP_dm_expected(k)*dfac_dlnPgasm1_const_T
               d_dlnPdm_dlnTm1_const_Pgas = d_dlnPdm_dlnTm1_const_Pgas*fac + &
                  s% dlnP_dm_expected(k)*dfac_dlnTm1_const_Pgas
            else
               d_dlnPdm_dlnd00 = d_dlnPdm_dlnd00*fac + s% dlnP_dm_expected(k)*dfac_dlnd00
               d_dlnPdm_dlnT00 = d_dlnPdm_dlnT00*fac + s% dlnP_dm_expected(k)*dfac_dlnT00
               d_dlnPdm_dlndm1 = d_dlnPdm_dlndm1*fac + s% dlnP_dm_expected(k)*dfac_dlndm1
               d_dlnPdm_dlnTm1 = d_dlnPdm_dlnTm1*fac + s% dlnP_dm_expected(k)*dfac_dlnTm1
            end if
            
            d_dlnPdm_dlnR00 = d_dlnPdm_dlnR00*fac + s% dlnP_dm_expected(k)*dfac_dlnR
            d_dlnPdm_dL00 = d_dlnPdm_dL00*fac + s% dlnP_dm_expected(k)*dfac_dL

         end if
         
         if (s% use_gr_factors) then ! GR gravity factor = 1/sqrt(1-2Gm/(rc^2))
            
            call get_gr_gravity_factor(s, k, gr_factor, d_gr_factor_dlnR00)
            
            if (dbg) then
               write(*,2) 'std s% dlnP_dm_expected(k)', k, s% dlnP_dm_expected(k)
               write(*,2) 'gr_factor', k, gr_factor
               write(*,2) 'new s% dlnP_dm_expected(k)', k, s% dlnP_dm_expected(k)*gr_factor
               !write(*,2) '', k, 
               !write(*,2) '', k, 
               !write(*,2) '', k, 
            end if

            s% dlnP_dm_expected(k) = s% dlnP_dm_expected(k)*gr_factor                        
            d_dlnPdm_dlnR00 = &
               d_dlnPdm_dlnR00*gr_factor + s% dlnP_dm_expected(k)*d_gr_factor_dlnR00
            
         end if
         
         if (s% v_flag) then ! include velocity term
            
            if (s% i_lnT == 0) then
               local_v_flag = .true.
            else
               local_v_flag = &
                  (s% xh_old(s% i_lnT,k)/ln10 >= s% velocity_logT_lower_bound)
            end if
            
            fac = s% accel_factor
            if (s% use_mass_corrections) fac = fac*s% mass_correction_start(k)
            
            if (local_v_flag) then
               dv_dt = s% dv_dt(k) ! = (v(k) - vstart(k))/dt
            else ! assume vstart(k) = 0 and
               ! constant acceleration dv_dt so vfinal(k) = dv_dt*dt
               ! v(k) = dr/dt = average velocity =
               !      (vstart + vfinal)/2 = dv_dt*dt/2 when vstart = 0
               ! so (1/2)*dv_dt*dt = v(k) 
               dv_dt = 2d0*s% v(k)*s% dVARDOT_dVAR
            end if
            accel = dv_dt
            
            R2 = s% R2(k) ! for time weighting
            rPterm = 1/(4*pi*R2*Ppoint)          
            d_rPterm_dlnR = -rPterm*s% d_R2_dlnR(k)/R2

            dlnP_dm_accel_term = -fac*accel*rPterm
            
            s% dlnP_dm_expected(k) = s% dlnP_dm_expected(k) + dlnP_dm_accel_term
            
            if (s% lnPgas_flag) then
               d_rPterm_dlnPgas00_const_T = -rPterm*dPpoint_dlnPgas00_const_T/Ppoint
               d_rPterm_dlnT00_const_Pgas = -rPterm*dPpoint_dlnT00_const_Pgas/Ppoint
               d_rPterm_dlnPgasm1_const_T = -rPterm*dPpoint_dlnPgasm1_const_T/Ppoint
               d_rPterm_dlnTm1_const_Pgas = -rPterm*dPpoint_dlnTm1_const_Pgas/Ppoint
               d_rPterm_dlnd00 = 0
               d_rPterm_dlnT00 = 0
               d_rPterm_dlndm1 = 0
               d_rPterm_dlnTm1 = 0
            else
               d_rPterm_dlnd00 = -rPterm*dPpoint_dlnd00/Ppoint
               d_rPterm_dlnT00 = -rPterm*dPpoint_dlnT00/Ppoint
               d_rPterm_dlndm1 = -rPterm*dPpoint_dlndm1/Ppoint
               d_rPterm_dlnTm1 = -rPterm*dPpoint_dlnTm1/Ppoint
               d_rPterm_dlnPgas00_const_T = 0
               d_rPterm_dlnT00_const_Pgas = 0
               d_rPterm_dlnPgasm1_const_T = 0
               d_rPterm_dlnTm1_const_Pgas = 0
            end if

            if (s% lnPgas_flag) then
               d_dlnPdm_dlnPgas00_const_T = &
                  d_dlnPdm_dlnPgas00_const_T - fac*accel*d_rPterm_dlnPgas00_const_T
               d_dlnPdm_dlnT00_const_Pgas = &
                  d_dlnPdm_dlnT00_const_Pgas - fac*accel*d_rPterm_dlnT00_const_Pgas
               d_dlnPdm_dlnPgasm1_const_T = &
                  d_dlnPdm_dlnPgasm1_const_T - fac*accel*d_rPterm_dlnPgasm1_const_T
               d_dlnPdm_dlnTm1_const_Pgas = &
                  d_dlnPdm_dlnTm1_const_Pgas - fac*accel*d_rPterm_dlnTm1_const_Pgas
            else
               d_dlnPdm_dlnd00 = d_dlnPdm_dlnd00 - fac*accel*d_rPterm_dlnd00
               d_dlnPdm_dlnT00 = d_dlnPdm_dlnT00 - fac*accel*d_rPterm_dlnT00
               d_dlnPdm_dlndm1 = d_dlnPdm_dlndm1 - fac*accel*d_rPterm_dlndm1
               d_dlnPdm_dlnTm1 = d_dlnPdm_dlnTm1 - fac*accel*d_rPterm_dlnTm1
            end if

            d_dlnPdm_dlnR00 = d_dlnPdm_dlnR00 - fac*accel*d_rPterm_dlnR
            
            if (local_v_flag) then
               d_dlnPdm_dvel00 = -fac*s% dVARDOT_dVAR*rPterm
            else
               d_dlnPdm_dvel00 = -fac*rPterm*2d0*s% dVARDOT_dVAR
            end if

            if (do_visc) then ! include artificial viscosity
            
               if (k > 1) then
                  dQ = s% Qvisc(k-1) - s% Qvisc(k)
                  d_dQ_dlndm1 = s% dQvisc_dlnd(k-1)
                  d_dQ_dlnRm1 = s% dQvisc_dlnR00(k-1)
                  d_dQ_dlnR00 = s% dQvisc_dlnRp1(k-1) - s% dQvisc_dlnR00(k)
                  d_dQ_dvelm1 = s% dQvisc_dvel00(k-1)
                  d_dQ_dvel00 = s% dQvisc_dvelp1(k-1) - s% dQvisc_dvel00(k)
               else
                  dQ = -s% Qvisc(k)
                  d_dQ_dlndm1 = 0
                  d_dQ_dlnRm1 = 0
                  d_dQ_dlnR00 = -s% dQvisc_dlnR00(k)
                  d_dQ_dvelm1 = 0
                  d_dQ_dvel00 = -s% dQvisc_dvel00(k)
               end if
               
               d_dQ_dlnd00 = -s% dQvisc_dlnd(k)
               d_dQ_dlnRp1 = -s% dQvisc_dlnRp1(k)
               d_dQ_dvelp1 = -s% dQvisc_dvelp1(k)
                              
               dm_bar = s% dm_bar(k)
               s% dvdt_visc(k) = (4*pi/r)*(dQ/dm_bar)
               
               if (.false.) then
                  write(*,2) 's% dvdt_visc(k)', k, s% dvdt_visc(k)
                  write(*,2) 'dQ', k, dQ
                  write(*,2) 'dm_bar', k, dm_bar
                  write(*,2) 'Ppoint', k, Ppoint                  
                  write(*,*)
               end if
               
               fac = s% dvdt_visc_factor
               if (s% use_mass_corrections) fac = fac*s% mass_correction_start(k)
               
               inv_Ppoint = 1/Ppoint
               d_inv_Ppoint_dlnd00 = -dPpoint_dlnd00/(Ppoint*Ppoint)
               d_inv_Ppoint_dlnT00 = -dPpoint_dlnT00/(Ppoint*Ppoint)
               d_inv_Ppoint_dlndm1 = -dPpoint_dlndm1/(Ppoint*Ppoint)
               d_inv_Ppoint_dlnTm1 = -dPpoint_dlnTm1/(Ppoint*Ppoint)

               vfac = fac*inv_Ppoint/(r*R2*dm_bar)
               dvfac_dlnd00 = fac*d_inv_Ppoint_dlnd00/(r*R2*dm_bar)
               dvfac_dlnT00 = fac*d_inv_Ppoint_dlnT00/(r*R2*dm_bar)
               dvfac_dlndm1 = fac*d_inv_Ppoint_dlndm1/(r*R2*dm_bar)
               dvfac_dlnTm1 = fac*d_inv_Ppoint_dlnTm1/(r*R2*dm_bar)
               dvfac_dlnR00 = -vfac*(1d0 + s% d_R2_dlnR(k)/R2)

               dlnP_dm_visc = vfac*dQ
               s% dlnP_dm_visc(k) = dlnP_dm_visc

               d_dlnPdm_visc_dlndm1 = dvfac_dlndm1*dQ + vfac*d_dQ_dlndm1
               d_dlnPdm_visc_dlnd00 = dvfac_dlnd00*dQ + vfac*d_dQ_dlnd00
               
               d_dlnPdm_visc_dlnTm1 = dvfac_dlnTm1*dQ
               d_dlnPdm_visc_dlnT00 = dvfac_dlnT00*dQ
               
               d_dlnPdm_visc_dlnRm1 = vfac*d_dQ_dlnRm1
               d_dlnPdm_visc_dlnR00 = vfac*d_dQ_dlnR00 + dvfac_dlnR00*dQ
               d_dlnPdm_visc_dlnRp1 = vfac*d_dQ_dlnRp1
               
               d_dlnPdm_visc_dvelm1 = vfac*d_dQ_dvelm1
               d_dlnPdm_visc_dvel00 = vfac*d_dQ_dvel00
               d_dlnPdm_visc_dvelp1 = vfac*d_dQ_dvelp1
                             
               s% dlnP_dm_expected(k) = s% dlnP_dm_expected(k) + dlnP_dm_visc

               d_dlnPdm_dlnR00 = d_dlnPdm_dlnR00 + d_dlnPdm_visc_dlnR00
               d_dlnPdm_dlnRp1 = d_dlnPdm_dlnRp1 + d_dlnPdm_visc_dlnRp1
               d_dlnPdm_dlnRm1 = d_dlnPdm_dlnRm1 + d_dlnPdm_visc_dlnRm1

               d_dlnPdm_dvel00 = d_dlnPdm_dvel00 + d_dlnPdm_visc_dvel00
               d_dlnPdm_dvelp1 = d_dlnPdm_dvelp1 + d_dlnPdm_visc_dvelp1
               d_dlnPdm_dvelm1 = d_dlnPdm_dvelm1 + d_dlnPdm_visc_dvelm1
               
               if (s% lnPgas_flag) then    
               
                  d_dlnPdm_visc_dlnPgas00_const_T = &
                     d_dlnPdm_visc_dlnd00*s% dlnRho_dlnPgas_const_T(k)                 
               
                  d_dlnPdm_visc_dlnT00_const_Pgas = &
                     d_dlnPdm_visc_dlnT00 + &
                     d_dlnPdm_visc_dlnd00*s% dlnRho_dlnT_const_Pgas(k)
                  
                  if (k > 1) then
                     d_dlnPdm_visc_dlnPgasm1_const_T = &
                        d_dlnPdm_visc_dlndm1*s% dlnRho_dlnPgas_const_T(k-1) 
                     d_dlnPdm_visc_dlnTm1_const_Pgas = &
                        d_dlnPdm_visc_dlnTm1 + &
                        d_dlnPdm_visc_dlndm1*s% dlnRho_dlnT_const_Pgas(k-1)
                  else
                     d_dlnPdm_visc_dlnPgasm1_const_T = 0
                     d_dlnPdm_visc_dlnTm1_const_Pgas = 0
                  end if
                         
                  d_dlnPdm_dlnPgas00_const_T = d_dlnPdm_dlnPgas00_const_T + &     
                     d_dlnPdm_visc_dlnPgas00_const_T
                     
                  d_dlnPdm_dlnPgasm1_const_T = d_dlnPdm_dlnPgasm1_const_T + &   
                     d_dlnPdm_visc_dlnPgasm1_const_T
                     
                  d_dlnPdm_dlnT00_const_Pgas = d_dlnPdm_dlnT00_const_Pgas + &
                     d_dlnPdm_visc_dlnT00_const_Pgas
                     
                  d_dlnPdm_dlnTm1_const_Pgas = d_dlnPdm_dlnTm1_const_Pgas + &
                     d_dlnPdm_visc_dlnTm1_const_Pgas
                    
               else         
                    
                  d_dlnPdm_dlnd00 = d_dlnPdm_dlnd00 + d_dlnPdm_visc_dlnd00
                  d_dlnPdm_dlndm1 = d_dlnPdm_dlndm1 + d_dlnPdm_visc_dlndm1          
                  d_dlnPdm_dlnT00 = d_dlnPdm_dlnT00 + d_dlnPdm_visc_dlnT00
                  d_dlnPdm_dlnTm1 = d_dlnPdm_dlnTm1 + d_dlnPdm_visc_dlnTm1        
                    
               end if
                        
            end if
                        
         end if            
               
         if (.false.) then
            write(*,2) 's% dlnP_dm_expected(k)', k, s% dlnP_dm_expected(k)
            write(*,*)
         end if

      end subroutine do_dlnP_dm
      
      
      subroutine eval_dlnPdm_qhse(s, k, &
            m, dlnP_dm_qhse, d_dlnPdm_dlnR, d_dlnPdm_dL, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(out) :: &
            m, dlnP_dm_qhse, d_dlnPdm_dlnR, d_dlnPdm_dL, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, r4, d_r4_dlnR, dlnq_dq, P00, Pm1, Ppoint_start, &
            dlnP_dm_qhse0, fs, dfs_dlnT, dfs_dlnR, dfs_dlnd, dfs_dL, &
            dfs_dlnPgas00_const_T, dfs_dlnT00_const_Pgas, correction_factor, &
            rtheta, d_rtheta_dlnR
         logical :: lnPgas_flag

         include 'formats'
         
         ierr = 0
         lnPgas_flag = s% lnPgas_flag
         
         ! basic eqn is dP/dm = -G m / (4 pi r^4)
         ! divide by <P> to make it unitless
         ! simple average is adequate for <P> since is only for normalizing the equation.
         ! however, be careful to use same <P> for both sides of equation..... 
         
         ! for rotation, multiply by factor fp
         
         ! for tau < 2/3, multiply by Paczynski factor for dilution of radiation
            ! B. Paczynski, 1969, Acta Astr., vol. 19.  eqn 13
                  
         ! dlnP_dm_qhse = -G m / (4 pi r^4 <P>)
         
         if (.not. s% use_energy_conservation_form) then ! for proper time weighting
            rtheta = s% r(k)
            d_rtheta_dlnR = s% r(k)
         else
            rtheta = s% r_start(k)
            d_rtheta_dlnR = 0
         end if
         r4 = s% R2(k)*s% r(k)*rtheta
         d_r4_dlnR = (s% d_R2_dlnR(k) + s% R2(k))*s% r(k)*rtheta + &
            s% R2(k)*s% r(k)*d_rtheta_dlnR
         
         P00 = s% P(k)
         if (k == 1) then
            alfa = 1
            Pm1 = 0d0
            Ppoint = alfa*P00
            dPpoint_dlndm1 = 0
            dPpoint_dlnTm1 = 0
            dPpoint_dlnPgasm1_const_T = 0
            dPpoint_dlnTm1_const_Pgas = 0
         else
            Pm1 = s% P(k-1)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*P00 + (1-alfa)*Pm1
            if (lnPgas_flag) then
               dPpoint_dlnPgasm1_const_T = (1-alfa)*s% Pgas(k-1)
               dPpoint_dlnTm1_const_Pgas = (1-alfa)*4*s% Prad(k-1)
            else
               dPpoint_dlndm1 = (1-alfa)*s% P(k-1)*s% chiRho(k-1)
               dPpoint_dlnTm1 = (1-alfa)*s% P(k-1)*s% chiT(k-1)
            end if
         end if
         
         if (lnPgas_flag) then
            dPpoint_dlnPgas00_const_T = alfa*s% Pgas(k)
            dPpoint_dlnT00_const_Pgas = alfa*4*s% Prad(k)
         else
            dPpoint_dlnd00 = alfa*s% P(k)*s% chiRho(k)
            dPpoint_dlnT00 = alfa*s% P(k)*s% chiT(k)
         end if
         

         if (s% use_mass_corrections) then ! skip time weighting of the corrections
            correction_factor = &
               s% mass_correction_start(k) + &
               s% P_div_rho_start(k)/(clight*clight)
               
            if (correction_factor > 1.05 .or. correction_factor < 0.95) then
               write(*,2) 'bad correction factor?', k, correction_factor
               write(*,2) 's% mass_correction_start(k)', k, s% mass_correction_start(k)
               write(*,2) 's% P_div_rho_start(k)/clight**2', &
                  k, s% P_div_rho_start(k)/(clight*clight)
               stop 1
            end if
            
            m = s% m_grav(k)
         else
            correction_factor = 1
            m = s% m(k)
         end if
         
         dlnP_dm_qhse = -s% cgrav(k)*m*correction_factor/(pi4*r4*Ppoint)
         
         if (s% rotation_flag .and. s% use_dP_dm_rotation_correction) then
            dlnP_dm_qhse = dlnP_dm_qhse*s% fp_rot(k)
         end if
            
         if (s% use_other_momentum) then
            if (k == 1) then
               Ppoint_start = s% P_start(1)
            else
               Ppoint_start = alfa*s% P_start(k) + (1-alfa)*s% P_start(k-1)
            end if
            dlnP_dm_qhse = dlnP_dm_qhse + s% extra_dPdm(k)/Ppoint_start
         end if
            
         dlnP_dm_qhse0 = dlnP_dm_qhse
         
         if (s% tau(k) >= 2d0/3d0) then
            fs = 0
            dfs_dlnT = 0
            dfs_dlnR = 0
            dfs_dlnd = 0
            dfs_dL = 0
         else ! (tau < 2d0/3d0) then ! B. Paczynski, 1969, Acta Astr., vol. 19, 1.
            fs = (1d0 - 1.5d0*s% tau(k))*(2*crad*pow3(s% T(k))*sqrt(s% r(k)))/ &
                  (3d0*s% cgrav(k)*m*s% rho(k))* &
                  pow_cr(s% L(k)/(8d0*pi*boltz_sigma),0.25d0)  ! eqn 15
            dfs_dlnT = 3d0*fs
            dfs_dlnR = 0.5d0*fs
            dfs_dlnd = -fs
            dfs_dL = 0.25d0*fs/s% L(k)
            dlnP_dm_qhse = dlnP_dm_qhse0*(1 + fs) ! eqn 13
         end if

         d_dlnPdm_dlnR = -(d_r4_dlnR/r4)*dlnP_dm_qhse + dlnP_dm_qhse0*dfs_dlnR
         dlnq_dq = 1/s% q(k)
         
         if (lnPgas_flag) then      
            dfs_dlnPgas00_const_T = dfs_dlnd*s% dlnRho_dlnPgas_const_T(k)
            dfs_dlnT00_const_Pgas = dfs_dlnT + dfs_dlnd*s% dlnRho_dlnT_const_Pgas(k)
            d_dlnPdm_dlnPgas00_const_T = &
               -(dPpoint_dlnPgas00_const_T/Ppoint)*dlnP_dm_qhse + &
               dlnP_dm_qhse0*dfs_dlnPgas00_const_T
            d_dlnPdm_dlnT00_const_Pgas = &
               -(dPpoint_dlnT00_const_Pgas/Ppoint)*dlnP_dm_qhse + &
               dlnP_dm_qhse0*dfs_dlnT00_const_Pgas
            d_dlnPdm_dlnPgasm1_const_T = -(dPpoint_dlnPgasm1_const_T/Ppoint)*dlnP_dm_qhse
            d_dlnPdm_dlnTm1_const_Pgas = -(dPpoint_dlnTm1_const_Pgas/Ppoint)*dlnP_dm_qhse
         else
            d_dlnPdm_dlnd00 = &
               -(dPpoint_dlnd00/Ppoint)*dlnP_dm_qhse + dlnP_dm_qhse0*dfs_dlnd
            d_dlnPdm_dlnT00 = &
               -(dPpoint_dlnT00/Ppoint)*dlnP_dm_qhse + dlnP_dm_qhse0*dfs_dlnT
            d_dlnPdm_dlndm1 = -(dPpoint_dlndm1/Ppoint)*dlnP_dm_qhse
            d_dlnPdm_dlnTm1 = -(dPpoint_dlnTm1/Ppoint)*dlnP_dm_qhse
         end if
         
         d_dlnPdm_dL = dlnP_dm_qhse0*dfs_dL

         if (is_bad_num(dlnP_dm_qhse)) then
            ierr = -1
            if (s% report_ierr) then
!$OMP critical
               write(*,*) 'eval_dlnPdm_qhse: is_bad_num(dlnP_dm_qhse)'
               write(*,2) 'dlnP_dm_qhse', k, dlnP_dm_qhse
               write(*,2) 's% tau(k)', k, s% tau(k)
               write(*,2) 's% fp_rot(k)', k, s% fp_rot(k)
               write(*,2) 'r4', k, r4
               write(*,2) 'Ppoint', k, Ppoint
               write(*,2) 'correction_factor', k, correction_factor
               write(*,2) 'm', k, m
               write(*,2) 's% cgrav(k)', k, s% cgrav(k)
               stop
!$OMP end critical
            end if
            return
         end if

      end subroutine eval_dlnPdm_qhse

      
      subroutine get_free_fall_surf_PT(s, &
            which_atm_option, off_table_option, r_phot, L_surf, skip_partials, &
            Teff, lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
         use crlibm_lib
         type (star_info), pointer :: s
         integer, intent(in) :: which_atm_option, off_table_option
         real(dp), intent(in) :: r_phot, L_surf
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr
         
         real(dp) :: mdot_Bondi, infall_entropy, infall_lnS, mdot, &
            G, M, z, lnR, dlnR, P, lnP_alt, P_alt, lnT_alt, T_alt
         include 'formats'

         ierr = 0
         Teff = 0
         lnT_surf = 0
         dlnT_dL = 0
         dlnT_dlnR = 0
         dlnT_dlnM = 0
         dlnT_dlnkap = 0
         lnP_surf = 0
         dlnP_dL = 0
         dlnP_dlnR = 0
         dlnP_dlnM = 0
         dlnP_dlnkap = 0
         
         mdot_Bondi = s% free_fall_mdot_Bondi
         !write(*,1) 'mdot_Bondi', mdot_Bondi
         infall_entropy = s% free_fall_entropy
         !write(*,1) 'infall_entropy', infall_entropy
         infall_lnS = log_cr(infall_entropy*avo*kerg)
         mdot = mdot_Bondi*(Msun/secyer)
         G = s% cgrav(1)
         M = s% mstar
         z = 1d0 - max(0d0,min(1d0,s% X(1)+s% Y(1)))
                  
         ! use rmid(1) instead of r_phot
         lnR = log_cr(s% rmid(1))
         call eval_PT(s% rmid(1), lnP_surf, P, lnT_surf, Teff, ierr)
         if (ierr /= 0) then
            write(*,1) 'infall_surface_PT failed in eval_PT', r_phot
            stop 'infall_surface_PT'
            return
         end if
         
         dlnR = 1d-7*lnR
         call eval_PT( &
            exp_cr(lnR+dlnR), lnP_alt, P_alt, lnT_alt, T_alt, ierr)
         if (ierr /= 0) then
            write(*,1) 'infall_surface_PT failed in eval_PT', exp_cr(lnR+dlnR)
            stop 'infall_surface_PT'
            return
         end if
         
         dlnT_dlnR = (lnT_alt - lnT_surf)/dlnR
         dlnP_dlnR = (lnP_alt - lnP_surf)/dlnR
         
         
         contains
         
         subroutine eval_PT(r, lnP, P, lnT, T, ierr)
            use eos_def, only: &
               i_lnS, i_lnPgas, num_eos_basic_results
            use eos_lib, only: Radiation_Pressure, eosDT_get_T
            use star_utils, only: get_XYZ
            real(dp), intent(in) :: r
            real(dp), intent(out) :: lnP, P, lnT, T
            integer, intent(out) :: ierr
            
            real(dp) :: logT_guess,logT_result, &
               log10Rho,dlnRho_dlnPgas_const_T,dlnRho_dlnT_const_Pgas, &
               v_ff, rho_ff, X, Y, Z, Prad, Pgas, logRho
            real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_const_TRho, d_dzbar_const_TRho
            integer, parameter :: max_iter = 100
            integer :: eos_calls
            real(dp), parameter :: logT_tol = 1.d-12, other_tol = 1.d-12, &
               logT_bnd1= arg_not_provided, logT_bnd2= arg_not_provided, &
               other_at_bnd1= arg_not_provided, other_at_bnd2= arg_not_provided
               
            include 'formats'
               
            ierr = 0
            lnP = 0
            lnT = 0
            T = 0
            
            ! use r to get rho_ff at surface
            v_ff = sqrt(2*G*M/r)
            !write(*,1) 'v_ff', v_ff
            !write(*,1) 'r/v_ff', r/v_ff
            rho_ff = mdot/(4*pi*r*r*v_ff)
            
            ! use rho_ff and infall entropy to get T and P at surface
            logRho = log10_cr(rho_ff)
            logT_guess = s% lnT(1)/ln10
            call get_XYZ(s, s% xa(:,1), X, Y, Z)
            call eosDT_get_T( &
               s% eos_handle, Z, X, s% abar(1), s% zbar(1), &
               s% species, s% chem_id, s% net_iso, s% xa(:,1), &
               logRho, i_lnS, infall_lnS, &
               logT_tol, other_tol, max_iter, logT_guess,  &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho,  &
               d_dabar_const_TRho, d_dzbar_const_TRho, eos_calls, ierr)
            if (ierr /= 0) then
               write(*,1) 'r', r
               write(*,1) 'mdot/(Msun/secyer)', mdot/(Msun/secyer)
               write(*,1) 'eosDT_get_T ierr'
               write(*,1) 'v_ff', v_ff
               write(*,1) 'rho_ff', rho_ff
               write(*,1) 'logRho', logRho
               write(*,1) 'logT_guess', logT_guess
               write(*,1) 'logS', infall_lnS/ln10
               write(*,1) 'z', z
               write(*,1) 'X', s% X(1)
               write(*,1) 'abar', s% abar(1)
               write(*,1) 'zbar', s% zbar(1)
               write(*,*)
               return
            end if
            
            lnT = ln10*logT_result
            Pgas = exp_cr(res(i_lnPgas))
            T = exp_cr(lnT)
            Prad = Radiation_Pressure(T)
            P = Prad + Pgas
            lnP = log_cr(P)

         end subroutine eval_PT
         
      end subroutine get_free_fall_surf_PT


      subroutine do1_dedt( &
            s, k, do_visc, skip_partials, dedt_expected, &
            d_dedt_dlnd, d_dedt_dlnT, &
            d_dedt_dlnR00, d_dedt_dlnRp1, &
            d_dedt_dv00, d_dedt_dvp1, &
            d_dedt_dL00, d_dedt_dLp1, &
            d_dedt_dlnd_c_E, d_dedt_dE_c_Rho, &
            ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: do_visc, skip_partials
         real(dp), intent(out) :: dedt_expected, &
            d_dedt_dlnd, d_dedt_dlnT, &
            d_dedt_dlnR00, d_dedt_dlnRp1, &
            d_dedt_dv00, d_dedt_dvp1, &
            d_dedt_dL00, d_dedt_dLp1, &
            d_dedt_dlnd_c_E, d_dedt_dE_c_Rho
         integer, intent(out) :: ierr
         
         real(dp) :: &
            dm, energy, eps_nuc, non_nuc_neu, theta, &
            Lp1, dR2v, P, PdVdt, rho, L00, visc_factor, &
            d_PdvAdm_dlnR00, d_PdvAdm_dv00, &
            d_PdvAdm_dlnRp1, d_PdvAdm_dvp1, d_PdVdt_dlnP
         integer :: nz
            
         ! de_dt_expected = 
         !   - (L(k) - L(k+1))/dm(k)
         !   - P*dVdt
         !   + eps_nuc
         !   - non_nuc_neutrinos
         !   + extra_heat
         !   + irradiation_heat
         !   + eps_visc

         
         ! dVdt = dvAdm = 4*pi*(R2(k)*vc(k) - R2(k+1)*vc(k+1))/dm(k)
         ! or
         ! dVdt = d(1/rho)/dt = (1/rho - 1/rho_start)/dt
         ! or
         ! dVdt = d(1/rho)/dt = -1/rho*dlnd_dt
         
         include 'formats'

         ierr = 0
         nz = s% nz
            
         energy = s% energy(k)
         eps_nuc = s% eps_nuc(k)
         non_nuc_neu = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))            
         s% eps_heat(k) = eps_nuc - non_nuc_neu + &
               s% extra_heat(k) + s% irradiation_heat(k)

         dm = s% dm(k)
         rho = s% rho(k)
         theta = s% theta_P
         P = theta*s% P(k) + (1d0-theta)*s% P_start(k)
         if (s% L_flag) then
            L00 = s% L(k)
            if (k < nz) then
               Lp1 = s% L(k+1)
            else
               Lp1 = s% L_center
            end if
         else
            L00 = 0
            if (k == nz) then
               Lp1 = s% L_center ! allow injection of energy from center
            else
               Lp1 = 0
            end if
         end if
         
         if (k < nz) then
            dR2v = s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1)
         else
            dR2v = s% R2(k)*s% vc(k) - s% r_center*s% r_center*s% v_center
         end if
         PdVdt = 4*pi*P*dR2v/dm

         dedt_expected = -(L00 - Lp1)/dm - PdVdt + s% eps_heat(k)

         if (do_visc) then
            visc_factor = s% eps_visc_factor
            dedt_expected = dedt_expected + s% eps_visc(k)*visc_factor
         else
            visc_factor = 0
         end if
         
         d_dedt_dlnd = 0
         d_dedt_dlnT = 0
         d_dedt_dlnR00 = 0
         d_dedt_dlnRp1 = 0
         d_dedt_dv00 = 0
         d_dedt_dvp1 = 0
         d_dedt_dL00 = 0
         d_dedt_dLp1 = 0
         d_dedt_dlnd_c_E = 0
         d_dedt_dE_c_Rho = 0
         
         if (is_bad_num(dedt_expected)) then
!$OMP critical
            write(*,2) 'dedt_expected', k, dedt_expected
            write(*,2) 's% eps_visc(k)', k, s% eps_visc(k)
            write(*,2) 'visc_factor', k, visc_factor
            write(*,2) 'eps_heat', k, s% eps_heat(k)
            write(*,2) 'PdVdt', k, PdVdt
            write(*,2) 'L00 - Lp1', k, L00 - Lp1
            stop 'do1_dedt'
!$OMP end critical
         end if
         
         if (skip_partials) return
         
         !dedt_expected = -(L00 - Lp1)/dm - PdVdt + eps_heat
         
         if (s% L_flag) then
            d_dedt_dL00 = -1/dm
            if (k < nz) d_dedt_dLp1 = 1/dm
         else
            d_dedt_dL00 = 0
            d_dedt_dLp1 = 0
         end if

         d_dedt_dlnT = &
            s% d_epsnuc_dlnT(k) - 0.5d0*s% d_nonnucneu_dlnT(k) + s% d_extra_heat_dlnT(k)
         d_dedt_dlnd = &
            s% d_epsnuc_dlnd(k) - 0.5d0*s% d_nonnucneu_dlnd(k) + s% d_extra_heat_dlnd(k)
            
         if (do_visc) &
            d_dedt_dlnd = d_dedt_dlnd + s% d_eps_visc_dlnd(k)*visc_factor
            
         !PdvAdm = 4*pi*P*(s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1))/dm
         d_PdvAdm_dlnR00 = 4*pi*P*s% d_R2_dlnR(k)*s% vc(k)/dm
         d_PdvAdm_dv00 = 4*pi*P*s% R2(k)*s% d_vc_dv/dm
         d_dedt_dlnR00 = -d_PdvAdm_dlnR00
         d_dedt_dv00 = -d_PdvAdm_dv00
         if (k < nz) then     
            d_PdvAdm_dlnRp1 = -4*pi*P*s% d_R2_dlnR(k+1)*s% vc(k+1)/dm
            d_PdvAdm_dvp1 = -4*pi*P*s% R2(k+1)*s% d_vc_dv/dm
            d_dedt_dlnRp1 = -d_PdvAdm_dlnRp1
            d_dedt_dvp1 = -d_PdvAdm_dvp1
         end if    
         
         if (k == -672) then
            write(*,*)
            write(*,2) 's% d_R2_dlnR(k)', k, s% d_R2_dlnR(k)
            write(*,2) 'd_PdvAdm_dlnR00', k, d_PdvAdm_dlnR00
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            !write(*,2) 'delta PdVdt', k, PdVdt - 4.7477338060905886D+03
         end if
         
         if (do_visc) then
            d_dedt_dlnR00 = d_dedt_dlnR00 + s% d_eps_visc_dlnR00(k)*visc_factor
            d_dedt_dv00 = d_dedt_dv00 + s% d_eps_visc_dvel00(k)*visc_factor
            if (k < nz) then
               d_dedt_dlnRp1 = d_dedt_dlnRp1 + s% d_eps_visc_dlnRp1(k)*visc_factor
               d_dedt_dvp1 = d_dedt_dvp1 + s% d_eps_visc_dvelp1(k)*visc_factor
            end if
         end if
                  
         d_PdVdt_dlnP = 4*pi*theta*s% P(k)*dR2v/dm
         
         if (s% E_flag) then
            ! 1st convert the previously set values
            d_dedt_dlnd_c_E = d_dedt_dlnd + d_dedt_dlnT*s% dlnT_dlnd_c_E(k)
            d_dedt_dE_c_Rho = d_dedt_dlnT*s% dlnT_dlnE_c_Rho(k)/energy
            ! then add the PdVdt terms
            d_dedt_dlnd_c_E = d_dedt_dlnd_c_E - d_PdVdt_dlnP*s% dlnP_dlnd_c_E(k)
            d_dedt_dE_c_Rho = d_dedt_dE_c_Rho - d_PdVdt_dlnP*s% dlnP_dlnE_c_Rho(k)/energy
         end if
         
         d_dedt_dlnT = d_dedt_dlnT - d_PdVdt_dlnP*s% chiT(k)
         d_dedt_dlnd = d_dedt_dlnd - d_PdVdt_dlnP*s% chiRho(k)
         
         if (k == -345) then
            write(*,2) 'dedt_expected', k, dedt_expected
            write(*,2) 'd_dedt_dlnd', k, d_dedt_dlnd
            write(*,2) 'd_dedt_dlnT', k, d_dedt_dlnT
            write(*,2) 'd_dedt_dlnR00', k, d_dedt_dlnR00
            write(*,2) 'd_dedt_dlnRp1', k, d_dedt_dlnRp1
            write(*,2) 'd_dedt_dv00', k, d_dedt_dv00
            write(*,2) 'd_dedt_dvp1', k, d_dedt_dvp1
            write(*,2) 'd_dedt_dL00', k, d_dedt_dL00
            write(*,2) 'd_dedt_dLp1', k, d_dedt_dLp1
            write(*,2) 'd_dedt_dlnd_c_E', k, d_dedt_dlnd_c_E
            write(*,2) 'd_dedt_dE_c_Rho', k, d_dedt_dE_c_Rho
            write(*,2) 'P', k, P
            write(*,2) 'rho', k, rho
         end if

      end subroutine do1_dedt


      subroutine do1_dvdt( &
            s, k, P_surf, do_visc, skip_partials, dvdt_expected, &
            d_dvdt_dlnd00, d_dvdt_dlndm1, &
            d_dvdt_dlnT00, d_dvdt_dlnTm1, &
            d_dvdt_dlnR00, d_dvdt_dlnRm1, d_dvdt_dlnRp1, &
            d_dvdt_dv00, d_dvdt_dvm1, d_dvdt_dvp1, &
            ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         logical, intent(in) :: do_visc, skip_partials
         real(dp), intent(out) :: dvdt_expected, &
            d_dvdt_dlnd00, d_dvdt_dlndm1, &
            d_dvdt_dlnT00, d_dvdt_dlnTm1, &
            d_dvdt_dlnR00, d_dvdt_dlnRm1, d_dvdt_dlnRp1, &
            d_dvdt_dv00, d_dvdt_dvm1, d_dvdt_dvp1
         integer, intent(out) :: ierr

         ! dv_dt_expected(k) = 
         !   - cgrav(k)*m(k)/(r(k)*rtheta(k))
         !   - (4*pi*R2(k))*(P(k-1) - P(k))/dm_bar(k)
         !   + (4*pi/r(k))*(Q(k-1) - Q(k))/dm_bar(k)
         
         ! someday add various correction factors to this.
         
         real(dp) :: &
            theta, rtheta, d_rtheta_dlnR, G, Pm1, P00, dvdt_grav, &
            r, R2, m, delm, Qm1, Q00, dQ, pi4_div_r_dm, dvdt_AdPdm, &
            d_dQ_dlndm1, d_dQ_dlnd00, &
            d_dQ_dlnRm1, d_dQ_dlnR00, d_dQ_dlnRp1, &
            d_dQ_dvelm1, d_dQ_dvel00, d_dQ_dvelp1
         
         include 'formats'
         
         ierr = 0
         if (.not. s% use_energy_conservation_form) then
            rtheta = s% r(k)
            d_rtheta_dlnR = s% r(k)
         else
            rtheta = s% r_start(k)
            d_rtheta_dlnR = 0
         end if
         
         G = s% cgrav(k)
         
         theta = s% theta_P
         if (k == 1) then
            Pm1 = P_surf
         else
            Pm1 = theta*s% P(k-1) + (1d0-theta)*s% P_start(k-1)
         end if
         P00 = theta*s% P(k) + (1d0-theta)*s% P_start(k)
         r = s% r(k)
         R2 = s% R2(k)
         m = s% m(k)
         !delm = 0.5d0*(s% dm(k-1) + s% dm(k)) ! difference between cell centers
         !   ! note delm is dm_bar for all k except nz
         delm = s% dm_bar(k)
         
         dvdt_grav = -G*m/(r*rtheta)
         dvdt_AdPdm = -4*pi*R2*(Pm1 - P00)/delm
         if (s% use_other_momentum) &
            dvdt_AdPdm = dvdt_AdPdm - 4*pi*R2*s% extra_dPdm(k)

         dvdt_expected = dvdt_grav + dvdt_AdPdm
            
         pi4_div_r_dm = 4*pi/(r*delm)
         if (do_visc) then
            if (k == 1) then
               Qm1 = 0
            else
               Qm1 = s% Qvisc(k-1)
            endif
            Q00 = s% Qvisc(k)
            dQ = Qm1 - Q00
            s% dvdt_visc(k) = pi4_div_r_dm*dQ
            dvdt_expected = dvdt_expected + s% dvdt_visc(k)
         else
            s% dvdt_visc(k) = 0
            dQ = 0
         end if
         
         if (skip_partials) then
            s% d_dvdt_visc_dlndm1(k) = 0
            s% d_dvdt_visc_dlnRm1(k) = 0
            s% d_dvdt_visc_dvm1(k) = 0
            s% d_dvdt_visc_dlnd00(k) = 0
            s% d_dvdt_visc_dlnR00(k) = 0
            s% d_dvdt_visc_dv00(k) = 0
            s% d_dvdt_visc_dlnRp1(k) = 0
            s% d_dvdt_visc_dvp1(k) = 0
            d_dvdt_dlnd00 = 0
            d_dvdt_dlndm1 = 0
            d_dvdt_dlnT00 = 0
            d_dvdt_dlnTm1 = 0
            d_dvdt_dlnR00 = 0
            d_dvdt_dlnRm1 = 0
            d_dvdt_dlnRp1 = 0
            d_dvdt_dv00 = 0
            d_dvdt_dvm1 = 0
            d_dvdt_dvp1 = 0
            return
         end if
         
         ! dvdt_AdPdm = -4*pi*R2*(Pm1 - P00)/delm
         d_dvdt_dlnd00 = 4*pi*R2*theta*s% P(k)*s% chiRho(k)/delm
         d_dvdt_dlnT00 = 4*pi*R2*theta*s% P(k)*s% chiT(k)/delm
         
         if (k > 1) then
            d_dvdt_dlndm1 = -4*pi*R2*theta*s% P(k-1)*s% chiRho(k-1)/delm
            d_dvdt_dlnTm1 = -4*pi*R2*theta*s% P(k-1)*s% chiT(k-1)/delm
         end if
         
         d_dvdt_dlnR00 = -4*pi*s% d_R2_dlnR(k)*(Pm1 - P00)/delm + &
            G*m*(1d0 + d_rtheta_dlnR/rtheta)/(r*rtheta)

         if (do_visc) then  
                      
            d_dQ_dlnd00 = -s% dQvisc_dlnd(k)           

            if (k > 1) then
               d_dQ_dlndm1 = s% dQvisc_dlnd(k-1)
               d_dQ_dlnRm1 = s% dQvisc_dlnR00(k-1)
               d_dQ_dvelm1 = s% dQvisc_dvel00(k-1)
               d_dQ_dlnR00 = s% dQvisc_dlnRp1(k-1) - s% dQvisc_dlnR00(k)
               d_dQ_dvel00 = s% dQvisc_dvelp1(k-1) - s% dQvisc_dvel00(k)
            else
               d_dQ_dlndm1 = 0
               d_dQ_dlnRm1 = 0
               d_dQ_dvelm1 = 0
               d_dQ_dlnR00 = -s% dQvisc_dlnR00(k)
               d_dQ_dvel00 = -s% dQvisc_dvel00(k)
            end if
            
            d_dQ_dlnRp1 = -s% dQvisc_dlnRp1(k)           
            d_dQ_dvelp1 = -s% dQvisc_dvelp1(k)
            
            s% d_dvdt_visc_dlndm1(k) = pi4_div_r_dm*d_dQ_dlndm1
            s% d_dvdt_visc_dlnRm1(k) = pi4_div_r_dm*d_dQ_dlnRm1
            s% d_dvdt_visc_dvm1(k) = pi4_div_r_dm*d_dQ_dvelm1
            s% d_dvdt_visc_dlnd00(k) = pi4_div_r_dm*d_dQ_dlnd00
            s% d_dvdt_visc_dlnR00(k) = pi4_div_r_dm*(d_dQ_dlnR00 - dQ)
            s% d_dvdt_visc_dv00(k) = pi4_div_r_dm*d_dQ_dvel00
            s% d_dvdt_visc_dlnRp1(k) = pi4_div_r_dm*d_dQ_dlnRp1
            s% d_dvdt_visc_dvp1(k) = pi4_div_r_dm*d_dQ_dvelp1
            
            d_dvdt_dlndm1 = d_dvdt_dlndm1 + s% d_dvdt_visc_dlndm1(k)
            d_dvdt_dlnRm1 = s% d_dvdt_visc_dlnRm1(k)
            d_dvdt_dvm1 = s% d_dvdt_visc_dvm1(k)

            d_dvdt_dlnd00 = d_dvdt_dlnd00 + s% d_dvdt_visc_dlnd00(k)
            d_dvdt_dlnR00 = d_dvdt_dlnR00 + s% d_dvdt_visc_dlnR00(k)
            d_dvdt_dv00 = s% d_dvdt_visc_dv00(k)

            d_dvdt_dlnRp1 = s% d_dvdt_visc_dlnRp1(k)
            d_dvdt_dvp1 = s% d_dvdt_visc_dvp1(k)
            
         else
            s% d_dvdt_visc_dlndm1(k) = 0
            s% d_dvdt_visc_dlnRm1(k) = 0
            s% d_dvdt_visc_dvm1(k) = 0
            s% d_dvdt_visc_dlnd00(k) = 0
            s% d_dvdt_visc_dlnR00(k) = 0
            s% d_dvdt_visc_dv00(k) = 0
            s% d_dvdt_visc_dlnRp1(k) = 0
            s% d_dvdt_visc_dvp1(k) = 0
            d_dvdt_dlnRm1 = 0
            d_dvdt_dlnRp1 = 0
            d_dvdt_dv00 = 0
            d_dvdt_dvm1 = 0
            d_dvdt_dvp1 = 0
         end if

         
      end subroutine do1_dvdt

         
      end module hydro_vars

