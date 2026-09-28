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

      module photo_out

      use star_private_def

      implicit none     
      
      contains


      subroutine output_star_photo(s,iounit)
         use rates_def,only:num_rvs
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         
         integer :: ierr, part_number, nz, nz_old, nz_older, len_history_col_spec
         integer :: k

         include 'formats'

         ierr = 0
         nz = s% nz
         nz_old = s% nz_old
         nz_older = s% nz_older
         part_number = 0 ! part_numbers are just a consistency check on the data file

         write(iounit) star_def_version
         
         call write_part_number(iounit)
         write(iounit) &
            s% generations, s% total_num_newton_iterations, &
            s% nz, s% nz_old, s% nz_older, &
            s% nvar_hydro, s% nvar_chem, s% nvar, &
            s% lnPgas_flag, s% v_flag, s% L_flag, s% E_flag, s% rotation_flag, &
            s% prev_Lmax, s% species, s% num_reactions, &
            s% model_number, s% model_number_old, s% model_number_older, &
            s% mstar, s% mstar_old, s% mstar_older, &
            s% xmstar, s% xmstar_old, s% xmstar_older, &
            s% M_center, s% M_center_old, s% M_center_older, &
            s% v_center, s% v_center_old, s% v_center_older, &
            s% R_center, s% R_center_old, s% R_center_older, &
            s% L_center, s% L_center_old, s% L_center_older, &
            s% total_radiation, s% total_radiation_old, s% total_radiation_older, &
            s% piston_vfinal_inward, s% piston_vfinal_inward_old, s% piston_vfinal_inward_older, &
            s% piston_alpha, s% piston_alpha_old, s% piston_alpha_older, &
            s% time, s% time_old, s% time_older, &
            s% ebdf_order, s% ebdf_order_old, s% ebdf_order_older, &
            s% ebdf_hold, s% ebdf_hold, s% ebdf_hold, &
            s% total_angular_momentum, s% total_angular_momentum_old, &
            s% total_angular_momentum_older, s% prev_create_atm_R0_div_R, &
            s% dt, s% dt_old, &
            s% have_previous_rotation_info, s% have_previous_conv_vel, &
            s% have_previous_D_mix, &
            s% initial_mass, s% initial_z, s% initial_y, &
            s% was_in_implicit_wind_limit, &
            s% using_revised_net_name, &
            s% revised_net_name, &
            s% revised_net_name_old, &
            s% revised_net_name_older, &
            s% using_revised_max_yr_dt,  &
            s% revised_max_yr_dt, &
            s% revised_max_yr_dt_old, &
            s% revised_max_yr_dt_older, &
            s% astero_using_revised_max_yr_dt, &
            s% astero_revised_max_yr_dt, &
            s% astero_revised_max_yr_dt_old, &
            s% astero_revised_max_yr_dt_older, &
            s% cumulative_extra_heating, &
            s% cumulative_extra_heating_old, &
            s% cumulative_extra_heating_older, &            
            s% cumulative_irradiation_heating, &
            s% cumulative_irradiation_heating_old, &
            s% cumulative_irradiation_heating_older, &            
            s% cumulative_nuclear_heating, &
            s% cumulative_nuclear_heating_old, &
            s% cumulative_nuclear_heating_older, &            
            s% cumulative_L_surf, &
            s% cumulative_L_surf_old, &
            s% cumulative_L_surf_older, &            
            s% cumulative_L_center, &
            s% cumulative_L_center_old, &
            s% cumulative_L_center_older, &            
            s% cumulative_non_nuc_neu_cooling, &
            s% cumulative_non_nuc_neu_cooling_old, &
            s% cumulative_non_nuc_neu_cooling_older, &            
            s% cumulative_sources_and_sinks, &
            s% cumulative_sources_and_sinks_old, &
            s% cumulative_sources_and_sinks_older, &
            s% cumulative_visc_heat_added, &
            s% cumulative_visc_heat_added_old, &
            s% cumulative_visc_heat_added_older, &
            s% cumulative_eps_grav, &
            s% cumulative_eps_grav_old, &
            s% cumulative_eps_grav_older, &
            s% cumulative_energy_error, &
            s% cumulative_energy_error_old, &
            s% cumulative_energy_error_older, &
            s% cumulative_acoustic_L, &
            s% cumulative_acoustic_L_old, &
            s% cumulative_acoustic_L_older, &
            s% cumulative_acoustic_L_center, &
            s% cumulative_acoustic_L_center_old, &
            s% cumulative_acoustic_L_center_older, &
            s% have_initial_energy_integrals, &
            s% total_internal_energy_initial, &
            s% total_gravitational_energy_initial, &
            s% total_linear_kinetic_energy_initial, &
            s% total_rotational_kinetic_energy_initial, &
            s% total_energy_initial, &
            s% free_fall_mdot_Bondi, s% free_fall_entropy, &
            s% using_free_fall_surface_PT
 
         write(iounit) s% net_name
         
         call write_part_number(iounit)
         write(iounit) &
            s% dq(1:nz), s% q(1:nz), s% xa(:,1:nz), s% xh(:,1:nz), s% omega(1:nz), &
            s% dlnd_dt(1:nz), s% dlnPgas_dt(1:nz), s% dlnT_dt(1:nz), s% dE_dt(1:nz), &
            s% eps_grav(1:nz), s% conv_vel(1:nz), s% lnT(1:nz)
         
         call write_part_number(iounit)
         if (s% generations > 1) &
            write(iounit) &
               s% lnT_old(1:nz_old), &
               s% conv_vel_old(1:nz_old), &
               s% nu_ST_old(1:nz_old), &
               s% D_ST_old(1:nz_old), &
               s% D_DSI_old(1:nz_old), &
               s% D_SH_old(1:nz_old), &
               s% D_SSI_old(1:nz_old), &
               s% D_ES_old(1:nz_old), &
               s% D_GSF_old(1:nz_old), &
               s% D_mix_old(1:nz_old), &
               s% omega_old(1:nz_old), &
               s% dq_old(1:nz_old), s% q_old(1:nz_old), &
               s% xa_old(:,1:nz_old), s% xh_old(:,1:nz_old)
            
         call write_part_number(iounit)
         if (s% generations > 2) &
            write(iounit) &
               s% lnT_older(1:nz_older), &
               s% conv_vel_older(1:nz_older), &
               s% nu_ST_older(1:nz_older), &
               s% D_ST_older(1:nz_older), &
               s% D_DSI_older(1:nz_older), &
               s% D_SH_older(1:nz_older), &
               s% D_SSI_older(1:nz_older), &
               s% D_ES_older(1:nz_older), &
               s% D_GSF_older(1:nz_older), &
               s% D_mix_older(1:nz_older), &
               s% omega_older(1:nz_older), &
               s% dq_older(1:nz_older), s% q_older(1:nz_older), &
               s% xa_older(:,1:nz_older), s% xh_older(:,1:nz_older)

         call write_part_number(iounit)
         write(iounit) &
            s% mstar_dot, s% mstar_dot_old, s% mstar_dot_older, &
            s% v_surf, s% v_surf_old, s% v_surf_older, &
            s% L_nuc_burn_total, s% L_nuc_burn_total_old, s% L_nuc_burn_total_older, &
            s% L_by_category, s% L_by_category_old, s% L_by_category_older, &
            s% gradT_excess_alpha, s% gradT_excess_alpha_old, s% gradT_excess_alpha_older, &
            s% dt_limit_ratio, s% dt_limit_ratio_old, s% dt_limit_ratio_older, &
            s% L_phot, s% L_phot_old, s% L_phot_older, s% T_surf, s% P_surf, &
            s% h1_czb_mass, s% h1_czb_mass_old, s% h1_czb_mass_older, s% h1_czb_mass_prev, &
            s% he_core_mass, s% he_core_mass_old, s% he_core_mass_older, &
            s% c_core_mass, s% c_core_mass_old, s% c_core_mass_older, &
            s% tau_base, s% Teff, s% Teff_old, s% Teff_older, &
            s% center_eps_nuc, s% center_eps_nuc_old, s% center_eps_nuc_older, &
            s% Lrad_div_Ledd_avg_surf, s% Lrad_div_Ledd_avg_surf_old, s% Lrad_div_Ledd_avg_surf_older, &
            s% w_div_w_crit_avg_surf, s% w_div_w_crit_avg_surf_old, s% w_div_w_crit_avg_surf_older, &
            s% have_done_TP, s% TP_state, s% TP_state_old, s% TP_state_older, &
            s% TP_count, s% TP_count_old, s% TP_count_older, &
            s% TP_M_H_on, s% TP_M_H_on_old, s% TP_M_H_on_older, &
            s% TP_M_H_min, s% TP_M_H_min_old, s% TP_M_H_min_older, &
            s% n_conv_regions, s% n_conv_regions_old, s% n_conv_regions_older, &
            s% cz_bot_mass(:), s% cz_bot_mass_old(:), s% cz_bot_mass_older(:), &
            s% cz_top_mass(:), s% cz_top_mass_old(:), s% cz_top_mass_older(:)

         write(iounit) &
            s% dt_next, s% i_lnd, s% i_lnT, s% i_E, s% i_lnR, &
            s% i_lum, s% i_lnPgas, s% i_v, &
            s% i_dv_dt, s% i_dlnT_dm, s% i_dlnd_dt, s% i_dE_dt, s% i_dlnR_dt
         
         write(iounit) &
            s% model_profile_filename, s% model_controls_filename, s% model_data_filename, &
            s% most_recent_profile_filename, s% most_recent_controls_filename, &
            s% most_recent_model_data_filename
         
         call write_part_number(iounit)
         write(iounit) &
            s% helium_ignition, s% carbon_ignition, &
            s% recent_log_header, s% phase_of_evolution, &
            s% prev_Tcntr1, s% prev_age1, s% prev_Tcntr2, s% prev_age2, s% prev_Tsurf, &
            s% prv_log_luminosity, s% prv_log_surface_temp, s% prv_log_center_temp, s% prv_log_center_density, &
            s% profile_age, s% post_he_age, s% prev_luminosity, s% ignition_center_xhe, s% he_luminosity_limit, &
            s% prev_cntr_rho, s% next_cntr_rho

         call write_part_number(iounit)
         write(iounit) &
            s% hydro_matrix_type, s% done_with_piston, s% done_with_center_flash, &
            s% num_newton_iterations, s% num_retries, s% num_backups, &
            s% number_of_backups_in_a_row, s% last_backup, s% initial_L_center, s% doing_center_flash, &
            s% initial_R_center, s% initial_v_center, &
            s% timestep_hold, s% model_number_for_last_retry, &
            s% mesh_call_number, s% hydro_call_number, s% diffusion_call_number, &
            s% initial_timestep, s% why_Tlim, s% result_reason, &
            s% need_to_update_history_now, s% need_to_save_profiles_now, s% save_profiles_model_priority, &
            s% doing_flash_wind, s% doing_rlo_wind, s% doing_nova_wind, s% most_recent_photo_name

         call write_part_number(iounit)
         write(iounit) s% len_extra_iwork, s% len_extra_work
         if (s% len_extra_iwork > 0) then
            write(iounit) s% extra_iwork(1:s% len_extra_iwork)
            write(iounit) s% extra_iwork_old(1:s% len_extra_iwork)
            write(iounit) s% extra_iwork_older(1:s% len_extra_iwork)
         end if
         if (s% len_extra_work > 0) then
            write(iounit) s% extra_work(1:s% len_extra_work)
            write(iounit) s% extra_work_old(1:s% len_extra_work)
            write(iounit) s% extra_work_older(1:s% len_extra_work)
         end if

         write(iounit) &
            s% ixtra1, s% ixtra2, s% ixtra3, s% ixtra4, s% ixtra5, &
            s% ixtra6, s% ixtra7, s% ixtra8, s% ixtra9, s% ixtra10, &
            s% ixtra11, s% ixtra12, s% ixtra13, s% ixtra14, s% ixtra15, &
            s% ixtra16, s% ixtra17, s% ixtra18, s% ixtra19, s% ixtra20, &
            s% ixtra21, s% ixtra22, s% ixtra23, s% ixtra24, s% ixtra25, &
            s% ixtra26, s% ixtra27, s% ixtra28, s% ixtra29, s% ixtra30
         write(iounit) &
            s% xtra1, s% xtra2, s% xtra3, s% xtra4, s% xtra5, &
            s% xtra6, s% xtra7, s% xtra8, s% xtra9, s% xtra10, &
            s% xtra11, s% xtra12, s% xtra13, s% xtra14, s% xtra15, &
            s% xtra16, s% xtra17, s% xtra18, s% xtra19, s% xtra20, &
            s% xtra21, s% xtra22, s% xtra23, s% xtra24, s% xtra25, &
            s% xtra26, s% xtra27, s% xtra28, s% xtra29, s% xtra30
         write(iounit) &
            s% lxtra1, s% lxtra2, s% lxtra3, s% lxtra4, s% lxtra5, &
            s% lxtra6, s% lxtra7, s% lxtra8, s% lxtra9, s% lxtra10, &
            s% lxtra11, s% lxtra12, s% lxtra13, s% lxtra14, s% lxtra15, &
            s% lxtra16, s% lxtra17, s% lxtra18, s% lxtra19, s% lxtra20, &
            s% lxtra21, s% lxtra22, s% lxtra23, s% lxtra24, s% lxtra25, &
            s% lxtra26, s% lxtra27, s% lxtra28, s% lxtra29, s% lxtra30

         write(iounit) &
            s% ixtra1_old, s% ixtra2_old, s% ixtra3_old, s% ixtra4_old, s% ixtra5_old, &
            s% ixtra6_old, s% ixtra7_old, s% ixtra8_old, s% ixtra9_old, s% ixtra10_old, &
            s% ixtra11_old, s% ixtra12_old, s% ixtra13_old, s% ixtra14_old, s% ixtra15_old, &
            s% ixtra16_old, s% ixtra17_old, s% ixtra18_old, s% ixtra19_old, s% ixtra20_old, &
            s% ixtra21_old, s% ixtra22_old, s% ixtra23_old, s% ixtra24_old, s% ixtra25_old, &
            s% ixtra26_old, s% ixtra27_old, s% ixtra28_old, s% ixtra29_old, s% ixtra30_old
         write(iounit) &
            s% xtra1_old, s% xtra2_old, s% xtra3_old, s% xtra4_old, s% xtra5_old, &
            s% xtra6_old, s% xtra7_old, s% xtra8_old, s% xtra9_old, s% xtra10_old, &
            s% xtra11_old, s% xtra12_old, s% xtra13_old, s% xtra14_old, s% xtra15_old, &
            s% xtra16_old, s% xtra17_old, s% xtra18_old, s% xtra19_old, s% xtra20_old, &
            s% xtra21_old, s% xtra22_old, s% xtra23_old, s% xtra24_old, s% xtra25_old, &
            s% xtra26_old, s% xtra27_old, s% xtra28_old, s% xtra29_old, s% xtra30_old
         write(iounit) &
            s% lxtra1_old, s% lxtra2_old, s% lxtra3_old, s% lxtra4_old, s% lxtra5_old, &
            s% lxtra6_old, s% lxtra7_old, s% lxtra8_old, s% lxtra9_old, s% lxtra10_old, &
            s% lxtra11_old, s% lxtra12_old, s% lxtra13_old, s% lxtra14_old, s% lxtra15_old, &
            s% lxtra16_old, s% lxtra17_old, s% lxtra18_old, s% lxtra19_old, s% lxtra20_old, &
            s% lxtra21_old, s% lxtra22_old, s% lxtra23_old, s% lxtra24_old, s% lxtra25_old, &
            s% lxtra26_old, s% lxtra27_old, s% lxtra28_old, s% lxtra29_old, s% lxtra30_old
         
         write(iounit) &
            s% ixtra1_older, s% ixtra2_older, s% ixtra3_older, s% ixtra4_older, s% ixtra5_older, &
            s% ixtra6_older, s% ixtra7_older, s% ixtra8_older, s% ixtra9_older, s% ixtra10_older, &
            s% ixtra11_older, s% ixtra12_older, s% ixtra13_older, s% ixtra14_older, s% ixtra15_older, &
            s% ixtra16_older, s% ixtra17_older, s% ixtra18_older, s% ixtra19_older, s% ixtra20_older, &
            s% ixtra21_older, s% ixtra22_older, s% ixtra23_older, s% ixtra24_older, s% ixtra25_older, &
            s% ixtra26_older, s% ixtra27_older, s% ixtra28_older, s% ixtra29_older, s% ixtra30_older
         write(iounit) &
            s% xtra1_older, s% xtra2_older, s% xtra3_older, s% xtra4_older, s% xtra5_older, &
            s% xtra6_older, s% xtra7_older, s% xtra8_older, s% xtra9_older, s% xtra10_older, &
            s% xtra11_older, s% xtra12_older, s% xtra13_older, s% xtra14_older, s% xtra15_older, &
            s% xtra16_older, s% xtra17_older, s% xtra18_older, s% xtra19_older, s% xtra20_older, &
            s% xtra21_older, s% xtra22_older, s% xtra23_older, s% xtra24_older, s% xtra25_older, &
            s% xtra26_older, s% xtra27_older, s% xtra28_older, s% xtra29_older, s% xtra30_older
         write(iounit) &
            s% lxtra1_older, s% lxtra2_older, s% lxtra3_older, s% lxtra4_older, s% lxtra5_older, &
            s% lxtra6_older, s% lxtra7_older, s% lxtra8_older, s% lxtra9_older, s% lxtra10_older, &
            s% lxtra11_older, s% lxtra12_older, s% lxtra13_older, s% lxtra14_older, s% lxtra15_older, &
            s% lxtra16_older, s% lxtra17_older, s% lxtra18_older, s% lxtra19_older, s% lxtra20_older, &
            s% lxtra21_older, s% lxtra22_older, s% lxtra23_older, s% lxtra24_older, s% lxtra25_older, &
            s% lxtra26_older, s% lxtra27_older, s% lxtra28_older, s% lxtra29_older, s% lxtra30_older
                  
         if (associated(s% history_column_spec)) then
            len_history_col_spec = size(s% history_column_spec)
            write(iounit) len_history_col_spec
            write(iounit) s% history_column_spec(1:len_history_col_spec)
         else
            write(iounit) 0 ! len_log_col_spec
         end if
         
         write(iounit)  &
            s% number_of_history_columns, s% model_number_of_history_values, &
            s% need_to_set_history_names_etc
         if (s% number_of_history_columns > 0) then
            write(iounit) s% history_value_is_integer(1:s% number_of_history_columns)
            do k=1,s% number_of_history_columns
               write(iounit) s% history_names(k)
            end do
         end if

         call write_part_number(iounit)

         call s% other_photo_write(s% id, iounit)
         
         call write_part_number(iounit)

         
         contains
         
         subroutine write_part_number(iounit)
            integer, intent(in) :: iounit
            part_number = part_number + 1
            write(iounit) part_number
         end subroutine write_part_number
         
      
      end subroutine output_star_photo


      end module photo_out
