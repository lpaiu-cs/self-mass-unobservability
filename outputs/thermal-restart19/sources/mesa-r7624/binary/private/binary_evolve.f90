! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton and Pablo Marchant
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


      module binary_evolve

      use const_def
      use crlibm_lib
      use star_lib
      use star_def
      use crlibm_lib
      use binary_def

      implicit none

      contains

      subroutine binarydata_init(b)
         use utils_lib, only: is_bad_num
         type (binary_info), pointer :: b
         logical :: evolve_both_stars
         logical :: trace_binary_rlo
         integer :: finish_step_result
         integer :: i
         include 'formats.inc'

         b% doing_first_model_of_run = .true.

         b% s_donor => b% s1
         b% s_accretor => b% s2
         b% d_i = 1
         b% a_i = 2

         b% max_timestep = 1d99
         b% change_factor = b% initial_change_factor

         initial_mass(1) = b% s1% mstar / Msun
         if (b% evolve_both_stars) then
            initial_mass(2) = b% s2% mstar / Msun
         else
            initial_mass(2) = b% m2
         end if

         b% m(1) = initial_mass(1)*Msun
         b% m(2) = initial_mass(2)*Msun
         b% r(1) = Rsun*b% s1% photosphere_r
         if (b% evolve_both_stars) then
            b% r(2) = Rsun*b% s2% photosphere_r
         else
            b% r(2) = 0
         end if
         
         b% eccentricity = b% initial_eccentricity

         if (b% initial_period_in_days <= 0) then ! calculate from initial_separation_in_Rsuns
            b% separation = b% initial_separation_in_Rsuns*Rsun
            b% period = &
               (2*pi)*sqrt(b% separation*b% separation*b% separation/&
                     (standard_cgrav*(b% m(1)+b% m(2))))
            b% initial_period_in_days = b% period / (24d0*60d0*60d0)
         else
            b% period = b% initial_period_in_days*(24d0*60d0*60d0)
            b% separation = &
               pow_cr((b% s1% cgrav(1)*(b% m(1)+b% m(2)))*(b% period/(2*pi))**2,1d0/3d0)
         end if
         b% angular_momentum_j = b% m(1) * b% m(2) * sqrt( b% s1% cgrav(1) *&
            b% separation * (1 - b% eccentricity**2) / (b% m(1) + b% m(2)) )

         b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1) * (1 - b% eccentricity) ) / &
             b% rl(1) / (1 - b% eccentricity) ! gap < 0 means out of contact 
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2) * (1 - b% eccentricity) ) / &
             b% rl(2) / (1 - b% eccentricity) ! gap < 0 means out of contact

         if (is_bad_num(b% rl_relative_gap(1))) stop 'binarydata_init'
         if (is_bad_num(b% rl_relative_gap(2))) stop 'binarydata_init'
         ! these will be adjusted properly by check_radiative_core
         b% have_radiative_core(1) = .true.
         b% have_radiative_core(2) = .true.
         
         ! Set all parameters nessessary for integration over the binary orbit
         ! 1) true anomaly = polar angle from periastron 0 -> 2pi
         do i = 1,b% anomaly_steps 
            b% theta_co(i) = (i-1) * (2 * pi) / b% anomaly_steps
         end do
         
         ! 2) time between periastron and polar angle theta 0 -> 1 (fraction of the
         !    orbital period)
         do i = 1,b% anomaly_steps ! time between periastron and polar angle theta
            b% time_co(i) = ( 2 * atan( sqrt( (1-b% eccentricity)/(1 + b% eccentricity) ) * &
                            tan(b% theta_co(i)/2d0) ) - b% eccentricity * &
                            sqrt(1 - b% eccentricity**2) * sin(b% theta_co(i)) / &
                            (1 + b% eccentricity * cos(b% theta_co(i)) ) ) /2. /pi
            if (i > b% anomaly_steps/2+1) then
               b% time_co(i) = b% time_co(i) + b% time_co(b% anomaly_steps/2+1) * 2
            end if
         end do
         
         write(*,*)
         write(*,1) 'm2', b% m2
         write(*,1) 'm1', b% m1
         write(*,1) 'initial_period_in_days', b% initial_period_in_days
         write(*,1) 'initial_separation_in_Rsun', b% separation/Rsun
         write(*,1) 'jdot_multiplier', b% jdot_multiplier
         write(*,1) 'fr', b% fr
         write(*,*)

         b% have_radiative_core = .false.
         just_evolve = .false.
         b% s1% mesh_delta_coeff_pre_ms = 1
         min_binary_period = b% period
         b% min_binary_separation = b% separation

         b% num_tries = 0

         finish_step_result = binary_finish_step(b)
          
      end subroutine

      subroutine set_donor_star(b)
         type (binary_info), pointer :: b
         logical :: switch_donor
         real(dp) :: mdot_hi_temp
         include 'formats.inc'

         switch_donor = .false.

         if (b% keep_donor_fixed .and. b% mdot_scheme /= "contact") return

         if (b% mdot_scheme == "roche_lobe" .and. &
            abs(b% mtransfer_rate/(Msun/secyer)) < b% mdot_limit_donor_switch .and. &
            b% rl_relative_gap_old(b% a_i) > b% rl_relative_gap_old(b% d_i)) then
            switch_donor = .true.
         else if (b% mtransfer_rate > 0d0) then
            switch_donor = .true.
            b% mtransfer_rate = - b% mtransfer_rate
            mdot_hi_temp = b% mdot_hi
            b% mdot_hi = - b% mdot_lo
            b% mdot_lo = - mdot_hi_temp
            if (.not. b% have_mdot_lo) then
               b% have_mdot_hi = .false.  
            end if
            b% have_mdot_lo = .true.
            b% fixed_delta_mdot = b% fixed_delta_mdot / 2.0
         else if (b% mdot_scheme == "contact" .and. &
            b% rl_relative_gap_old(b% a_i) > b% rl_relative_gap_old(b% d_i) .and. &
            b% rl_relative_gap_old(b% a_i) < - b% implicit_scheme_tolerance .and. &
            b% rl_relative_gap_old(b% d_i) < - b% implicit_scheme_tolerance .and. &
            abs(b% mtransfer_rate/(Msun/secyer)) < b% mdot_limit_donor_switch) then
            switch_donor = .true.
         end if

         if (switch_donor) then
            write(*,*) "switching donor"
            if (b% d_i == 2) then
               b% d_i = 1
               b% a_i = 2
               b% s_donor => b% s1
               b% s_accretor => b% s2
            else
               b% d_i = 2
               b% a_i = 1
               b% s_donor => b% s2
               b% s_accretor => b% s1
            end if
         end if
      end subroutine

      subroutine binary_evolve_step(b)
         use utils_lib, only: is_bad_num
         use binary_jdot, only: get_jdot
         use binary_edot, only: get_edot
         type(binary_info), pointer :: b
         integer :: i
         
         include 'formats.inc'

         b% m(1) = b% s1% mstar
         if (b% evolve_both_stars) then
            b% m(2) = b% s2% mstar
         else
            b% m(2) = b% m(2) - (b% xfer_fraction*b% mtransfer_rate + &
               b% mdot_wind_transfer(b% d_i)) * b% s1% dt
         end if
         
         b% r(1) = b% s1% photosphere_r*Rsun ! radius at photosphere in cm
         if (b% evolve_both_stars) b% r(2) = b% s2% photosphere_r*Rsun ! radius at photosphere in cm

         ! solve the winds in the system for jdot calculation,
         ! these don't include mass lost due to mass_transfer_efficiency < 1.0
         b% mdot_system_wind(b% d_i) = b% s_donor% mstar_dot - b% mtransfer_rate + &
            b% mdot_wind_transfer(b% a_i)
         if (b% evolve_both_stars) then
            b% mdot_system_wind(b% a_i) = b% s_accretor% mstar_dot + &
                b% mtransfer_rate * b% xfer_fraction + b% mdot_wind_transfer(b% d_i)
         else
            b% mdot_system_wind(b% a_i) = 0.0d0
         end if

         ! get jdot and update orbital J
         b% jdot = get_jdot(b, b% mtransfer_rate, b% xfer_fraction)
         b% angular_momentum_j = b% angular_momentum_j + b% jdot*b% s1% time_step*secyer

         if (b% angular_momentum_j <= 0) then
            stop 'bad angular_momentum_j'
         end if
         
         ! update the eccentricity (ignore in first step)
         if (.not. b% doing_first_model_of_run) then
            b% eccentricity = b% eccentricity + get_edot(b) *b% s1% time_step*secyer
            if (b% eccentricity < b% min_eccentricity) b% eccentricity = b% min_eccentricity
            if (b% eccentricity > b% max_eccentricity) b% eccentricity = b% max_eccentricity
         end if
         
         !use new eccentricity to calculate new time coordinate
         do i = 1,b% anomaly_steps ! time between periastron and polar angle theta
            b% time_co(i) = ( 2 * atan( sqrt( (1-b% eccentricity)/(1 + b% eccentricity) ) * &
                            tan(b% theta_co(i)/2d0) ) - b% eccentricity * &
                            sqrt(1 - b% eccentricity**2) * sin(b% theta_co(i)) / &
                            (1 + b% eccentricity * cos(b% theta_co(i)) ) ) /2. /pi
            if (i > b% anomaly_steps/2+1) then
               b% time_co(i) = b% time_co(i) + b% time_co(b% anomaly_steps/2+1) * 2
            end if
         end do
         
         ! use the new j to calculate new separation
         b% separation = ((b% angular_momentum_j/(b% m(1)*b% m(2)))**2) *&
             (b% m(1)+b% m(2)) / b% s1% cgrav(1) * 1 / (1 - b% eccentricity**2)
         if (b% separation < b% min_binary_separation) &
            b% min_binary_separation = b% separation
         
         b% period = 2*pi*sqrt(b% separation*b% separation*b% separation/&
               (b% s1% cgrav(1)*(b% m(1)+b% m(2)))) 
         if (b% period < min_binary_period) min_binary_period = b% period
         
         ! use the new separation to calculate the new roche lobe radius
         
         b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1) * (1 - b% eccentricity) ) / &
             b% rl(1) / (1 - b% eccentricity) ! gap < 0 means out of contact 
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2) * (1 - b% eccentricity) ) / &
             b% rl(2) / (1 - b% eccentricity) ! gap < 0 means out of contact

         if (is_bad_num(b% rl_relative_gap(1)) .or. is_bad_num(b% rl_relative_gap(2))) then
            stop 'error solving rl_rel_gap'
         end if

      end subroutine

      integer function binary_check_model(b)
         use binary_mdot, only: rlo_mdot, check_implicit_rlo
         use binary_irradiation
         type (binary_info), pointer :: b

         integer :: i, j, ierr, id
         logical :: implicit_rlo
         real(dp) :: new_mdot, q


         include 'formats.inc'

         binary_check_model = retry
         ierr = 0
         
         implicit_rlo = (b% max_tries_to_achieve > 0 .and. b% implicit_scheme_tolerance > 0d0)
         
         binary_check_model = keep_going
                  
         if (.not. just_evolve) then
            if (implicit_rlo) then ! check agreement between new r and new rl
               b% s_donor% min_abs_mdot_for_change_limits = 1d99
               binary_check_model = check_implicit_rlo(b, new_mdot)
               if (binary_check_model == keep_going) then
                  b% num_tries = 0
               end if
            else
               if (.not. b% use_other_rlo_mdot) then
                  call rlo_mdot(b% binary_id, new_mdot, ierr) ! grams per second
                  if (ierr /= 0) then
                     write(*,*) 'failed in rlo_mdot'
                     binary_check_model = retry
                     return
                  end if
               else
                  call b% other_rlo_mdot(b% binary_id, new_mdot, ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed in other rlo_mdot'
                     binary_check_model = retry
                     return
                  end if
               end if
               if (new_mdot > 0) then
                  new_mdot = 0.0d0
                  write(*,*) "WARNING: explicit computation of mass transfer results in accreting donor"
                  write(*,*) "Not transfering mass"
               end if
               ! smooth out the changes in mdot
               new_mdot = b% cur_mdot_frac*b% mtransfer_rate + (1-b% cur_mdot_frac)*new_mdot
               if (-new_mdot/(Msun/secyer) > b% max_explicit_abs_mdot) new_mdot = -b% max_explicit_abs_mdot*Msun/secyer 
            end if
            b% mtransfer_rate = new_mdot
            call adjust_irradiation(b, b% mtransfer_rate, b% xfer_fraction)
         end if

         if (b% evolve_both_stars .and. b% rl_relative_gap(b% a_i) >= 0.0d0) then
            if (b% rl_relative_gap(b% a_i) >= b% accretor_overflow_terminate) then
               binary_check_model = terminate
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = &
                   "Terminate because accretor (r-rl)/rl > accretor_overflow_terminate"
            end if
            write(*,'(a)') 'accretor photosphere has reached its Roche Lobe'
            write(*,2) 'accretor (r-rl)/rl', b% s_accretor% model_number, b% rl_relative_gap(b% a_i)
         end if
         if (b% doing_first_model_of_run .and. b% terminate_if_initial_overflow) then
            if (b% rl_relative_gap(1) >= 0.0d0 .or. (b% evolve_both_stars .and. b% rl_relative_gap(2) >= 0.0d0)) then
               binary_check_model = terminate
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = &
                   "Terminate because of overflowing initial model"
            end if
         end if
         if (b% evolve_both_stars .and. b% mdot_scheme == "contact") then
            if (b% m(1) > b% m(2)) then
               q = b% m(2) / b% m(1)
               id = 2
            else
               q = b% m(1) / b% m(2)
               id = 1
            end if
            if (b% rl_relative_gap(id) > 0.29858997d0*atan_cr(1.83530121d0*pow_cr(q,0.39661426d0))) then
               binary_check_model = terminate
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = &
                   "Terminate because of L2 overflow"
            end if
         end if

      end function binary_check_model

      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow_cr(m1/m2,one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
      end function eval_rlobe

      integer function binary_finish_step(b)
         type (binary_info), pointer :: b
         real(dp) :: spin_period

         if (b% do_tidal_sync) then
            if (mod(b% s1% model_number, b% s1% terminal_interval) == 0) then
               spin_period = 2*pi/b% s1% omega_avg_surf
               write(*,*) 'star_1_spin_period/binary_period', &
                  spin_period/b% period, spin_period/(60*60*24), b% period/(60*60*24)
               if (b% evolve_both_stars) then
                  spin_period = 2*pi/b% s2% omega_avg_surf
                  write(*,*) 'star_2_spin_period/binary_period', &
                     spin_period/b% period, spin_period/(60*60*24), b% period/(60*60*24)
               end if
            end if
         end if

         binary_finish_step = keep_going
         ! update change factor in case mtransfer_rate has changed
         if(b% mtransfer_rate_old /= b% mtransfer_rate .and. &
             b% mtransfer_rate /= 0 .and. b% mtransfer_rate_old/=0.0) then
            if(b% mtransfer_rate < b% mtransfer_rate_old) then
               b% change_factor = b% change_factor*(1.0-b% implicit_lambda) + b% implicit_lambda* &
                  (1+b% change_factor_fraction*(b% mtransfer_rate/b% mtransfer_rate_old-1))
            else
               b% change_factor = b% change_factor*(1.0-b% implicit_lambda) + b% implicit_lambda* &
                   (1+b% change_factor_fraction*(b% mtransfer_rate_old/b% mtransfer_rate-1))
            end if
            if(b% change_factor > b% max_change_factor) b% change_factor = b% max_change_factor
            if(b% change_factor < b% min_change_factor) b% change_factor = b% min_change_factor
         end if

         ! store all variables into "old" and "older"
         b% mtransfer_rate_older = b% mtransfer_rate_old
         b% angular_momentum_j_older = b% angular_momentum_j_old
         b% separation_older = b% separation_old
         b% eccentricity_older = b% eccentricity_old
         b% dt_older = b% dt_old
         b% env_older = b% env_old
         b% xfer_fraction_older = b% xfer_fraction_old
         b% sum_div_qloc_older(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc_older(2) = b% sum_div_qloc_old(2)
         b% period_older = b% period_old
         b% rl_relative_gap_older(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap_older(2) = b% rl_relative_gap_old(2)
         b% r_older(1) = b% r_old(1)
         b% r_older(2) = b% r_old(2)
         b% rl_older(1) = b% rl_old(1)
         b% rl_older(2) = b% rl_old(2)
         b% m_older(1) = b% m_old(1)
         b% m_older(2) = b% m_old(2)
         b% have_radiative_core_older = b% have_radiative_core_old
         b% max_timestep_older = b% max_timestep_old
         b% change_factor_older = b% change_factor_old

         b% d_i_older = b% d_i_old
         b% a_i_older = b% a_i_old

         b% dt_why_reason_older = b% dt_why_reason_old

         b% mtransfer_rate_old = b% mtransfer_rate
         b% angular_momentum_j_old = b% angular_momentum_j
         b% separation_old = b% separation
         b% eccentricity_old = b% eccentricity
         b% dt_old = b% dt
         b% env_old = b% env
         b% xfer_fraction_old = b% xfer_fraction
         b% sum_div_qloc_old(1) = b% sum_div_qloc(1)
         b% sum_div_qloc_old(2) = b% sum_div_qloc(2)
         b% period_old = b% period
         b% rl_relative_gap_old(1) = b% rl_relative_gap(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap(2)
         b% r_old(1) = b% r(1)
         b% r_old(2) = b% r(2)
         b% rl_old(1) = b% rl(1)
         b% rl_old(2) = b% rl(2)
         b% m_old(1) = b% m(1)
         b% m_old(2) = b% m(2)
         b% have_radiative_core_old = b% have_radiative_core
         b% max_timestep_old = b% max_timestep
         b% change_factor_old = b% change_factor

         b% d_i_old = b% d_i
         b% a_i_old = b% a_i

         b% dt_why_reason_old = b% dt_why_reason

      end function binary_finish_step

      integer function binary_prepare_to_redo(b)
         type (binary_info), pointer :: b

         binary_prepare_to_redo = redo
         call binary_set_current_to_old(b)

      end function binary_prepare_to_redo

      integer function binary_prepare_to_retry(b)
         type (binary_info), pointer :: b

         b% num_tries = 0
         ! this call takes care of restoring variables
         call binary_set_current_to_old(b)
         binary_prepare_to_retry = retry

      end function binary_prepare_to_retry

      integer function binary_do1_backup(b)
         type (binary_info), pointer :: b

         binary_do1_backup = backup
         if (b% s1% generations > 2) then
            call binary_restore_older(b)
         end if
         call binary_set_current_to_old(b)

         b% num_tries = 0

      end function binary_do1_backup

      subroutine binary_set_current_to_old(b)
         type (binary_info), pointer :: b
         ! restore variables
         ! do not restore mtransfer_rate during implicit rlo
         if (b% num_tries == 0) b% mtransfer_rate = b% mtransfer_rate_old
         b% angular_momentum_j = b% angular_momentum_j_old
         b% separation = b% separation_old
         b% eccentricity = b% eccentricity_old
         b% dt = b% dt_old
         b% env = b% env_old
         b% xfer_fraction = b% xfer_fraction_old
         b% sum_div_qloc(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc(2) = b% sum_div_qloc_old(2)
         b% period = b% period_old
         b% rl_relative_gap(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap(2) = b% rl_relative_gap_old(2)
         b% r(1) = b% r_old(1)
         b% r(2) = b% r_old(2)
         b% rl(1) = b% rl_old(1)
         b% rl(2) = b% rl_old(2)
         b% m(1) = b% m_old(1)
         b% m(2) = b% m_old(2)
         b% have_radiative_core = b% have_radiative_core_old
         b% max_timestep = b% max_timestep_old
         b% change_factor = b% change_factor_old

         ! do not restore donor and accretor ids during implicit mdot
         if (b% num_tries == 0) b% d_i = b% d_i_old
         if (b% num_tries == 0) b% a_i = b% a_i_old

         b% dt_why_reason = b% dt_why_reason_old
      end subroutine binary_set_current_to_old

      subroutine binary_restore_older(b)
         type (binary_info), pointer :: b

         b% mtransfer_rate_old = b% mtransfer_rate_older
         b% angular_momentum_j_old = b% angular_momentum_j_older
         b% separation_old = b% separation_older
         b% eccentricity_old = b% eccentricity_older
         b% dt_old = b% dt_older
         b% env_old = b% env_older
         b% xfer_fraction_old = b% xfer_fraction_older
         b% sum_div_qloc_old(1) = b% sum_div_qloc_older(1)
         b% sum_div_qloc_old(2) = b% sum_div_qloc_older(2)
         b% period_old = b% period_older
         b% rl_relative_gap_old(1) = b% rl_relative_gap_older(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap_older(2)
         b% r_old(1) = b% r_older(1)
         b% r_old(2) = b% r_older(2)
         b% rl_old(1) = b% rl_older(1)
         b% rl_old(2) = b% rl_older(2)
         b% m_old(1) = b% m_older(1)
         b% m_old(2) = b% m_older(2)
         b% have_radiative_core_old = b% have_radiative_core_older
         b% max_timestep_old = b% max_timestep_older
         b% change_factor_old = b% change_factor_older

         b% d_i_old = b% d_i_older
         b% a_i_old = b% a_i_older

         b% dt_why_reason_old = b% dt_why_reason_older

      end subroutine binary_restore_older

      end module binary_evolve
