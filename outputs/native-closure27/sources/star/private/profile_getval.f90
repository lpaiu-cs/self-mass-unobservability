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

      module profile_getval
      
      use star_private_def
      use star_profile_def
      use const_def
      use utils_lib, only: is_bad_num
      use star_utils, only: interp_val_to_pt

      implicit none
      
      
      integer, parameter :: idel = 10000
      
      integer, parameter :: add_abundances = idel
      integer, parameter :: add_log_abundances = add_abundances + 1
      integer, parameter :: category_offset = add_log_abundances + 1
      integer, parameter :: abundance_offset = category_offset + idel
      integer, parameter :: log_abundance_offset = abundance_offset + idel
      integer, parameter :: xadot_offset = log_abundance_offset + idel
      integer, parameter :: xaprev_offset = xadot_offset + idel
      integer, parameter :: ionization_offset = xaprev_offset + idel
      integer, parameter :: typical_charge_offset = ionization_offset + idel
      integer, parameter :: edv_offset = typical_charge_offset + idel
      integer, parameter :: v_rad_offset = edv_offset + idel
      integer, parameter :: log_g_rad_offset = v_rad_offset + idel
      integer, parameter :: log_concentration_offset = log_g_rad_offset + idel
      integer, parameter :: diffusion_dX_offset = log_concentration_offset + idel
      integer, parameter :: diffusion_D_offset = diffusion_dX_offset + idel
      integer, parameter :: extra_offset = diffusion_D_offset + idel
      
      integer, parameter :: max_profile_offset = extra_offset + idel
      
      
      contains
      
      
      subroutine getval_for_profile(s, c, k, val, int_flag, int_val)
         use chem_def
         use rates_def
         use mlt_def
         use nse_def, only: min_T_for_NSE
         use ionization_def
         use ionization_lib, only: eval_typical_charge
         use star_utils, only: get_log_concentration
         use star_utils, only: interp_val_to_pt, omega_crit, &
            get_Ledd, get_Lrad, get_phi_Joss, get_Lrad_div_Ledd, &
            get_total_energy_integral, eval_cell_collapse_time
         type (star_info), pointer :: s
         integer, intent(in) :: c, k
         real(dp), intent(out) :: val
         integer, intent(out) :: int_val
         logical, intent(inout) :: int_flag 
         
         real(dp) :: cno, z, x, frac, eps, eps_alt, L_rad, L_edd, P_face, rho_face, dr, &
            r00, rp1, v00, vp1, A00, Ap1, Amid, r00_start, rp1_start, dr3, dr3_start
         integer :: j, nz, ionization_k, klo, khi, i, ii, kk
         real(dp) :: ionization_res(num_ion_vals)
         real(dp) :: f, lgT, full_on, full_off, prefact         
         include 'formats'

         if (s% rotation_flag) then
            full_on = s% D_mix_rotation_max_logT_full_on
            full_off = s% D_mix_rotation_min_logT_full_off
            lgT = s% lnT(k)/ln10
            if (lgT <= full_on) then
               f = 1d0
            else if (lgT >= full_off) then
               f = 0d0
            else                   ! lgT > full_on and < full_off
               f = (lgT - full_on) / (full_off - full_on)
            end if
            prefact = f*s% am_D_mix_factor
         else
            prefact = 1d0
         end if

         val = 0; int_val = 0; int_flag = .false.
         nz = s% nz
         ionization_k = 0

         int_flag = .false.

         if (c > extra_offset) then
            i = c - extra_offset
            val = s% profile_extra(k,i)
         else if (c > diffusion_D_offset) then
            i = c - diffusion_D_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% diffusion_D_self(ii,k)
         else if (c > diffusion_dX_offset) then
            i = c - diffusion_dX_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% diffusion_dX(ii,k)
         else if (c > log_concentration_offset) then
            i = c - log_concentration_offset
            ii = s% net_iso(i)
            if (ii > 0) val = get_log_concentration(s,ii,k)
         else if (c > log_g_rad_offset) then
            i = c - log_g_rad_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = safe_log10_cr(s% g_rad(ii,k))
         else if (c > v_rad_offset) then
            i = c - v_rad_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% v_rad(ii,k)
         else if (c > edv_offset) then
            i = c - edv_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% edv(ii,k)
         else if (c > typical_charge_offset) then
            i = c - typical_charge_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% typical_charge(ii,k)
         else if (c > ionization_offset) then
            i = c - ionization_offset
            ii = s% net_iso(i)
            val = eval_typical_charge( &
               i, s% abar(k), exp_cr(s% lnfree_e(k))*s% abar(k), &
               s% T(k), s% lnT(k)/ln10, s% rho(k), s% lnd(k)/ln10)
         else if (c > xaprev_offset) then
            i = c - xaprev_offset
            ii = s% net_iso(i)
            if (ii > 0) val = s% xa_pre(ii,k)
         else if (c > xadot_offset) then
            i = c - xadot_offset
            ii = s% net_iso(i)
            if (ii > 0) val = s% xa(ii,k) - s% xa_pre(ii,k)
         else if (c > log_abundance_offset) then
            i = c - log_abundance_offset
            ii = s% net_iso(i)
            if (ii > 0) then
               val = safe_log10_cr(s% xa(ii,k))
            else
               val = -99d0
            end if
         else if (c > abundance_offset) then
            i = c - abundance_offset
            ii = s% net_iso(i)
            if (ii > 0) val = s% xa(ii,k)
         else if (c > category_offset) then
            i = c - category_offset
            val = s% eps_nuc_categories(i,k)
         else
   
         select case(c)
            case (p_zone)
               val = dble(k)
               int_val = k
               int_flag = .true.
            case (p_conv_L_div_L)
               val = s% L_conv_div_L(k)
            case (p_log_conv_L_div_L)
               val = safe_log10_cr(s% L_conv_div_L(k))
            case (p_luminosity)
               val = s% L(k)/Lsun
            case (p_lum_rad)
               L_rad = get_Lrad(s,k)
               val = L_rad/Lsun
            case (p_lum_conv)
               L_rad = get_Lrad(s,k)
               val = (s% L(k) - L_rad)/Lsun
               
            case (p_lum_rad_div_L_Edd_sub_fourPrad_div_PchiT)
               val = get_Lrad_div_Ledd(s,k) - 4*s% Prad(k)/(s% P(k)*s% chiT(k))
               
            case (p_lum_rad_div_L_Edd)
               val = get_Lrad_div_Ledd(s,k)
               
            case (p_lum_rad_div_L_Edd_div_phi_Joss)
               val = get_Lrad_div_Ledd(s,k)/get_phi_Joss(s,k)
               
            case (p_lum_conv_div_lum_Edd)
               L_rad = get_Lrad(s,k)
               L_edd = get_Ledd(s,k)
               val = (s% L(k) - L_rad)/L_edd
               
            case (p_lum_conv_div_lum_rad)
               L_rad = get_Lrad(s,k)
               val = (s% L(k) - L_rad)/L_rad

            case (p_lum_rad_div_L)
               L_rad = get_Lrad(s,k)
               val = L_rad/max(1d0,s% L(k))
            case (p_lum_conv_div_L)
               L_rad = get_Lrad(s,k)
               val = (s% L(k) - L_rad)/max(1d0,s% L(k))

            case (p_log_Lrad)
               L_rad = get_Lrad(s,k)
               val = safe_log10_cr(L_rad/Lsun)
            case (p_log_Lconv)
               L_rad = get_Lrad(s,k)
               val = safe_log10_cr((s% L(k) - L_rad)/Lsun)
               
            case (p_log_Lrad_div_L)
               L_rad = get_Lrad(s,k)
               val = safe_log10_cr(L_rad/s% L(k))
            case (p_log_Lrad_div_Ledd)
               val = safe_log10_cr(get_Lrad_div_Ledd(s,k))
               
            case (p_log_g)
               val = safe_log10_cr(s% grav(k))
            case (p_grav)
               val = s% grav(k)
            case (p_g_div_r)
               val = s% grav(k)/s% r(k)
            case (p_r_div_g)
               val = s% r(k)/s% grav(k)
            case (p_signed_log_eps_grav)
               val = s% eps_grav(k)
               val = sign(1d0,val)*log10_cr(max(1d0,abs(val)))
            case (p_net_nuclear_energy)
               val = s% eps_nuc(k) - s% eps_nuc_neu_total(k) - s% non_nuc_neu(k)
               val = sign(1d0,val)*log10_cr(max(1d0,abs(val)))
            case (p_eps_nuc_plus_nuc_neu)
               val = s% eps_nuc(k) + s% eps_nuc_neu_total(k)
            case (p_eps_nuc_minus_non_nuc_neu)
               val = s% eps_nuc(k) - s% non_nuc_neu(k)
            case (p_log_abs_eps_nuc_div_non_nuc_neu)
               val = log10_cr(max(1d0, abs(s% eps_nuc(k))) / max(1d-20,s% non_nuc_neu(k)))
            case (p_net_energy)
               val = s% eps_nuc(k) - s% non_nuc_neu(k) + s% eps_grav(k)
               val = sign(1d0,val)*log10_cr(max(1d0,abs(val)))
            case (p_signed_log_power)
               val = s% L(k)
               val = sign(1d0,val)*log10_cr(max(1d0,abs(val)))
            case (p_logL)
               val = log10_cr(max(1d-12,s% L(k)/Lsun))
            case (p_log_Ledd)
               val = log10_cr(get_Ledd(s,k)/Lsun)
            case (p_log_L_div_Ledd)
               val = log10_cr(max(1d-12,s% L(k)/get_Ledd(s,k)))
            case (p_velocity)
               val = s% v(k)
            case (p_vel_km_per_s)
               val = s% v(k)*1d-5
            case (p_v_div_r)
               val = s% v(k)/s% r(k)
            case (p_radius)
               val = s% r(k)/Rsun
            case (p_radius_cm)
               val = s% r(k)
            case (p_rmid)
               val = s% rmid(k)/Rsun
            case (p_logR_cm)
               val = safe_log10_cr(s% r(k))
            case (p_logR)
               val = safe_log10_cr(s% r(k)/Rsun)
               
            case (p_grav_gr_factor)
               val = 1d0/sqrt(1d0 - &
                        2*s% cgrav(k)*s% m_grav(k)/(s% r(k)*clight*clight))
            case (p_log_grav_gr_factor)
               val = safe_log10_cr(1d0/sqrt(1d0 - &
                        2*s% cgrav(k)*s% m_grav(k)/(s% r(k)*clight*clight)))

            case (p_avg_charge_H)
               val = get_ion_info(ion_iZ_H,k)
            case (p_avg_charge_He)
               val = get_ion_info(ion_iZ_He,k)
            case (p_avg_charge_C)
               val = get_ion_info(ion_iZ_C,k)
            case (p_avg_charge_N)
               val = get_ion_info(ion_iZ_N,k)
            case (p_avg_charge_O)
               val = get_ion_info(ion_iZ_O,k)
            case (p_avg_charge_Ne)
               val = get_ion_info(ion_iZ_Ne,k)
            case (p_avg_charge_Mg)
               val = get_ion_info(ion_iZ_Mg,k)
            case (p_avg_charge_Si)
               val = get_ion_info(ion_iZ_Si,k)
            case (p_avg_charge_Fe)
               val = get_ion_info(ion_iZ_Fe,k)
 
            case (p_neutral_fraction_H)
               val = get_ion_info(ion_ifneut_H,k)
            case (p_neutral_fraction_He)
               val = get_ion_info(ion_ifneut_He,k)
            case (p_neutral_fraction_C)
               val = get_ion_info(ion_ifneut_C,k)
            case (p_neutral_fraction_N)
               val = get_ion_info(ion_ifneut_N,k)
            case (p_neutral_fraction_O)
               val = get_ion_info(ion_ifneut_O,k)
            case (p_neutral_fraction_Ne)
               val = get_ion_info(ion_ifneut_Ne,k)
            case (p_neutral_fraction_Mg)
               val = get_ion_info(ion_ifneut_Mg,k)
            case (p_neutral_fraction_Si)
               val = get_ion_info(ion_ifneut_Si,k)
            case (p_neutral_fraction_Fe)
               val = get_ion_info(ion_ifneut_Fe,k)
 
            case (p_q)
               val = s% q(k)
            case (p_log_q)
               val = safe_log10_cr(s% q(k))
            case (p_dq)
               val = s% dq(k)
            case (p_logtau)
               val = safe_log10_cr(s% tau(k))
            case (p_mass)
               val = s% m(k)/Msun
            case (p_mass_grams)
               val = s% m(k)
            case (p_mmid)
               val = (s% M_center + s% xmstar*(s% q(k) - s% dq(k)/2))/Msun               

            case (p_dm)
               val = s% dm(k)
            case (p_dm_bar)
               val = s% dm_bar(k)

            case (p_m_div_r)
               val = s% m(k)/s% r(k)
            case (p_dmbar_m_div_r)
               val = s% dm_bar(k)*s% m(k)/s% r(k)
            case (p_log_dmbar_m_div_r)
               val = log10_cr(s% dm_bar(k)*s% m(k)/s% r(k))
               
            case (p_m_grav)
               val = s% m_grav(k)/Msun
            case (p_m_grav_div_m_baryonic)
               val = s% m_grav(k)/s% m(k)
            case (p_mass_correction_factor)
               val = s% mass_correction(k)                              
               
            case (p_x)
               val = s% X(k)
            case (p_log_x)
               val = safe_log10_cr(s% X(k))
            case (p_y)
               val = s% Y(k)
            case (p_log_y)
               val = safe_log10_cr(s% Y(k))
            case (p_z)
               val = 1d0 - (s% X(k) + s% Y(k))
            case (p_log_z)
               val = safe_log10_cr(1d0 - (s% X(k) + s% Y(k)))
            case (p_xm)
               val = s% xmstar*sum(s% dq(1:k-1))/Msun
            case (p_logxm)
               val = safe_log10_cr(s% xmstar*sum(s% dq(1:k-1))/Msun)
            case (p_logxq)
               val = safe_log10_cr(sum(s% dq(1:k-1)))
            case (p_logdq)
               val = safe_log10_cr(s% dq(k))
            case (p_log_column_depth)
               val = safe_log10_cr(s% xmstar*sum(s% dq(1:k-1))/(4*pi*s% r(k)*s% r(k)))
            case (p_log_radial_depth)
               val = safe_log10_cr(s% r(1) - s% r(k))
               
            case (p_rmid_div_R)
               val = s% rmid(k)/s% r(1)
            case (p_r_div_R)
               val = s% r(k)/s% r(1)
            case (p_log_dr)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10_cr(val)
            case (p_dr_div_rmid)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = val/s% rmid(k)
            case (p_log_dr_div_rmid)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10_cr(val/s% rmid(k))
            case (p_log_acoustic_radius)
               val = safe_log10_cr(sum(s% dr_div_csound(k:nz)))
            case (p_acoustic_radius)
               val = sum(s% dr_div_csound(k:nz))
            case (p_log_acoustic_depth)
               if (k > 1) &
                  val = sum(s% dr_div_csound(1:k-1))
               val = safe_log10_cr(val)
            case (p_acoustic_depth)
               if (k > 1) &
                  val = sum(s% dr_div_csound(1:k-1))
            case (p_acoustic_r_div_R_phot)
               val = sum(s% dr_div_csound(k:nz))/s% photosphere_acoustic_r

            case (p_lnR_residual)
               val = s% lnR_residual(k)
            case (p_lnd_residual)
               val = s% lnd_residual(k)
            case (p_dlnTdm_residual)
               val = s% dlnTdm_residual(k)
            case (p_E_residual)
               val = s% E_residual(k)
            case (p_v_residual)
               val = s% v_residual(k)

            case (p_sum_lnR_residual)
               val = sum(s% lnR_residual(k:s% nz))
            case (p_sum_lnd_residual)
               val = sum(s% lnd_residual(k:s% nz))
            case (p_sum_E_residual)
               val = sum(s% E_residual(k:s% nz))
            case (p_sum_v_residual)
               val = sum(s% v_residual(k:s% nz))
            case (p_sum_structure_residuals)
               val = sum(s% lnR_residual(k:s% nz) + &
                         s% lnd_residual(k:s% nz) + &
                         s% E_residual(k:s% nz) + &
                         s% v_residual(k:s% nz))

            case (p_log_lnR_residual)
               val = safe_log10_cr(abs(s% lnR_residual(k)))
            case (p_log_lnd_residual)
               val = safe_log10_cr(abs(s% lnd_residual(k)))
            case (p_log_dlnTdm_residual)
               val = safe_log10_cr(abs(s% dlnTdm_residual(k)))
            case (p_log_E_residual)
               val = safe_log10_cr(abs(s% E_residual(k)))
            case (p_log_v_residual)
               val = safe_log10_cr(abs(s% v_residual(k)))

            case (p_d_IE_dt_residual)
               val = s% d_IE_dt_actual(k) - s% d_IE_dt_expected(k)
            case (p_d_KE_dt_residual)
               val = s% d_KE_dt_actual(k) - s% d_KE_dt_expected(k)
            case (p_d_PE_dt_residual)
               val = s% d_PE_dt_actual(k) - s% d_PE_dt_expected(k)
            case (p_sum_from_cntr_dEdt_residuals)
               val = sum( &
                     s% d_IE_dt_actual(k:nz) - s% d_IE_dt_expected(k:nz) + &
                     s% d_KE_dt_actual(k:nz) - s% d_KE_dt_expected(k:nz) + &
                     s% d_PE_dt_actual(k:nz) - s% d_PE_dt_expected(k:nz))
               
            case (p_log_d_IE_dt_residual)
               val = s% d_IE_dt_actual(k) - s% d_IE_dt_expected(k)
               val = safe_log10_cr(abs(val))
            case (p_log_d_KE_dt_residual)
               val = s% d_KE_dt_actual(k) - s% d_KE_dt_expected(k)
               val = safe_log10_cr(abs(val))
            case (p_log_d_PE_dt_residual)
               val = s% d_PE_dt_actual(k) - s% d_PE_dt_expected(k)
               val = safe_log10_cr(abs(val))
            case (p_log_sum_from_cntr_dEdt_residuals)
               val = sum( &
                     s% d_IE_dt_actual(k:nz) - s% d_IE_dt_expected(k:nz) + &
                     s% d_KE_dt_actual(k:nz) - s% d_KE_dt_expected(k:nz) + &
                     s% d_PE_dt_actual(k:nz) - s% d_PE_dt_expected(k:nz))
               val = safe_log10_cr(abs(val))

            case (p_dEdt_expected)
               val = &
                  s% d_IE_dt_expected(k) + &
                  s% d_KE_dt_expected(k) + &
                  s% d_PE_dt_expected(k)

            case (p_check_dEdt_expected)
               val = 0
               ! this assumes v_center = 0 and eps_heat = 0 and probably other stuff too
               if (k > 1) &
                  val = val + &
                     4*pi*s% v(k)*(s% R2(k)*s% P(k-1) - s% Qvisc(k-1)/s% r(k))
               if (k < nz) &
                  val = val + &
                     4*pi*s% v(k+1)*(s% Qvisc(k)/s% r(k+1) - s% P(k)*s% R2(k+1))
               if (k == nz) &
                  val = val - s% L_center
               if (.false. .and. k == 1) then
                  write(*,1) 'val(1)', val
                  write(*,1) '-4*pi*s% v(k+1)*s% P(k)*s% R2(k+1)', &
                     -4*pi*s% v(k+1)*s% P(k)*s% R2(k+1)
                  write(*,1) '4*pi*s% v(k+1)*s% Qvisc(k)/s% r(k+1)', &
                     4*pi*s% v(k+1)*s% Qvisc(k)/s% r(k+1)
                  write(*,1) 's% d_KE_dt_expected(1)', s% d_KE_dt_expected(1)
                  write(*,1) 's% d_IE_dt_expected(1)', s% d_IE_dt_expected(1)
                  write(*,1) 's% d_PE_dt_expected(1)', s% d_PE_dt_expected(1)
                  write(*,1) 'sum', &
                     s% d_IE_dt_expected(k) + &
                     s% d_KE_dt_expected(k) + &
                     s% d_PE_dt_expected(k)
               end if
               val = val + &
                  s% d_IE_dt_expected(k) + &
                  s% d_KE_dt_expected(k) + &
                  s% d_PE_dt_expected(k)

            case (p_sum_check_dEdt_expected)
               val = 0
               ! this assumes v_center = 0 and eps_heat = 0 and probably other stuff too
               do kk=k,nz
                  if (kk > 1) &
                     val = val + &
                        4*pi*s% v(kk)*(s% R2(kk)*s% P(kk-1) - s% Qvisc(kk-1)/s% r(kk))
                  if (kk < nz) &
                     val = val + &
                        4*pi*s% v(kk+1)*(s% Qvisc(kk)/s% r(kk+1) - s% P(kk)*s% R2(kk+1))
                  if (kk == nz) &
                     val = val - s% L_center
                  val = val + &
                     s% d_IE_dt_expected(kk) + &
                     s% d_KE_dt_expected(kk) + &
                     s% d_PE_dt_expected(kk)
               end do

            case (p_dEdt_actual)
               val = &
                  s% d_IE_dt_actual(k) + &
                  s% d_KE_dt_actual(k) + &
                  s% d_PE_dt_actual(k)

            case (p_check_dEdt_actual)
               val = 0
               ! this assumes v_center = 0 and eps_heat = 0 and probably other stuff too
               if (k > 1) &
                  val = val + &
                     4*pi*s% v(k)*(s% R2(k)*s% P(k-1) - s% Qvisc(k-1)/s% r(k))
               if (k < nz) &
                  val = val + &
                     4*pi*s% v(k+1)*(s% Qvisc(k)/s% r(k+1) - s% P(k)*s% R2(k+1))
               if (k == nz) &
                  val = val - s% L_center
               val = val + &
                  s% d_IE_dt_actual(k) + &
                  s% d_KE_dt_actual(k) + &
                  s% d_PE_dt_actual(k)

            case (p_sum_check_dEdt_actual)
               val = 0
               ! this assumes v_center = 0 and eps_heat = 0 and probably other stuff too
               do kk=k,nz
                  if (kk > 1) &
                     val = val + &
                        4*pi*s% v(kk)*(s% R2(kk)*s% P(kk-1) - s% Qvisc(kk-1)/s% r(kk))
                  if (kk < nz) &
                     val = val + &
                        4*pi*s% v(kk+1)*(s% Qvisc(kk)/s% r(kk+1) - s% P(kk)*s% R2(kk+1))
                  if (kk == nz) &
                     val = val - s% L_center
                  val = val + &
                     s% d_IE_dt_actual(kk) + &
                     s% d_KE_dt_actual(kk) + &
                     s% d_PE_dt_actual(kk)
               end do
               
            case (p_dr_div_cs)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = val/s% csound(k)
            case (p_log_dr_div_cs)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10_cr(val/s% csound(k))
            case (p_dr_div_cs_yr)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = val/s% csound(k)/secyer
            case (p_log_dr_div_cs_yr)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10_cr(val/s% csound(k)/secyer)

            case (p_pgas_div_ptotal)
               val = s% Pgas(k)/s% P(k)

            case (p_pgas_div_p)
               val = s% Pgas(k)/s% P(k)

            case (p_cell_collapse_time)
               if (s% v_flag) then
                  if (k == s% nz) then
                     rp1 = s% R_center
                     vp1 = s% v_center
                  else
                     rp1 = s% r(k+1)
                     vp1 = s% v(k+1)
                  end if
                  r00 = s% r(k)
                  v00 = s% v(k)
                  if (vp1 > v00) val = (r00 - rp1)/(vp1 - v00)
               end if
               
            case (p_log_cell_collapse_time)
               if (s% v_flag) then
                  if (k == s% nz) then
                     rp1 = s% R_center
                     vp1 = s% v_center
                  else
                     rp1 = s% r(k+1)
                     vp1 = s% v(k+1)
                  end if
                  r00 = s% r(k)
                  v00 = s% v(k)
                  if (vp1 > v00) val = (r00 - rp1)/(vp1 - v00)
               end if
               val = safe_log10_cr(val)
               
            case (p_dlnT_dlnE_c_Rho)
               val = s% dlnT_dlnE_c_Rho(k)
            case (p_dlnPgas_dlnE_c_Rho)
               val = s% dlnPgas_dlnE_c_Rho(k)
            case (p_dlnP_dlnE_c_Rho)
               val = s% dlnP_dlnE_c_Rho(k)
            case (p_dlnT_dlnd_c_E)
               val = s% dlnT_dlnd_c_E(k)
            case (p_dlnPgas_dlnd_c_E)
               val = s% dlnPgas_dlnd_c_E(k)
            case (p_dlnP_dlnd_c_E)
               val = s% dlnP_dlnd_c_E(k)
               
            case (p_pturb_div_pgas_plus_prad)
               val = s% gamma1(k)/3*pow2(s% conv_vel(k)/s% csound(k))

            case (p_fourPrad_div_PchiT)
               val = 4*s% Prad(k)/(s% P(k)*s% chiT(k))

            case (p_dq_ratio)
               if (k == 1 .or. k == s% nz) then
                  val = 1
               else
                  val = s% dq(k-1)/s% dq(k)
               end if
            case (p_tau)
               val = s% tau(k)
            case (p_kap_frac_Type2)
               val = s% kap_frac_Type2(k)
            case (p_log_opacity)
               val = safe_log10_cr(s% opacity(k))
            case (p_energy)
               val = s% energy(k)
            case (p_logM)
               val = safe_log10_cr(s% m(k)/Msun)
            case (p_temperature)
               val = s% T(k)
            case (p_logT)
               val = s% lnT(k)/ln10
            case (p_density)
               val = s% rho(k)
            case (p_rho)
               val = s% rho(k)
            case (p_logRho)
               val = s% lnd(k)/ln10
            case (p_pgas)
               val = s% Pgas(k)
            case (p_logPgas)
               val = s% lnPgas(k)/ln10
            case (p_prad)
               val = s% Prad(k)
            case (p_pressure)
               val = s% P(k)
            case (p_logP)
               val = s% lnP(k)/ln10
            case (p_logE)
               val = s% lnE(k)/ln10
            case (p_grada)
               val = s% grada(k)
            case (p_dE_dRho)
               val = s% dE_dRho(k)
            case (p_cv)
               val = s% cv(k)
            case (p_cp)
               val = s% cp(k)

            case (p_thermal_time_to_surface)
               if (s% L(1) > 0) &
                  val = sum(s% dm(1:k)*s% cp(1:k)*s% T(1:k))/s% L(1)
            case (p_log_thermal_time_to_surface)
               if (s% L(1) > 0) then
                  val = sum(s% dm(1:k)*s% cp(1:k)*s% T(1:k))/s% L(1)
                  val = safe_log10_cr(val)
               end if
               
            case (p_log_CpT)
               val = safe_log10_cr(s% cp(k)*s% T(k))
            case (p_log_CpT_absMdot_div_L)
               val = safe_log10_cr(s% cp(k)*s% T(k)*abs(s% mstar_dot)/max(1d-99,s% L(k)))
            case (p_logS)
               val = s% lnS(k)/ln10
            case (p_logS_per_baryon)
               val = s% lnS(k)/ln10 + log10_cr(amu)
            case (p_gamma1)
               val = s% gamma1(k)
            case (p_gamma3)
               val = s% gamma3(k)
            case (p_eta)
               val = s% eta(k)
            case (p_theta_e)
               val = s% theta_e(k)
            case (p_gam)
               val = s% gam(k)
            case (p_mu)
               val = s% mu(k)
            case (p_v_div_cs)
               val = s% v_div_csound(k)
            case (p_v_div_csound)
               val = s% v_div_csound(k)
            case (p_csound)
               val = s% csound(k)
            case (p_csound_at_face)
               val = s% csound_at_face(k)
            case (p_scale_height)
               val = s% scale_height(k)/Rsun
            case (p_scale_height_start)
               val = s% scale_height_start(k)/Rsun
            case (p_entropy)
               val = s% entropy(k)
            case (p_free_e)
               val = exp_cr(s% lnfree_e(k))
            case (p_logfree_e)
               val = s% lnfree_e(k)/ln10
            case (p_chiRho)
               val = s% chiRho(k)
            case (p_chiT)
               val = s% chiT(k)

            case (p_dlnRho_dlnT_const_Pgas)
               if (s% lnPgas_flag) then
                  val = s% dlnRho_dlnT_const_Pgas(k)
               else
                  val = 0
               end if
            case (p_dlnRho_dlnPgas_const_T)
               if (s% lnPgas_flag) then
                  val = s% dlnRho_dlnPgas_const_T(k)
               else
                  val = 0
               end if

            case (p_x_mass_fraction_H)
               val = s% X(k)
            case (p_y_mass_fraction_He)
               val = s% Y(k)
            case (p_z_mass_fraction_metals)
               val = 1d0 - (s% X(k) + s% Y(k))

            case (p_abar)
               val = s% abar(k)
            case (p_zbar)
               val = s% zbar(k)
            case (p_z2bar)
               val = s% z2bar(k)
            case (p_ye)
               val = s% ye(k)
            case (p_opacity)
               val = s% opacity(k)
            case (p_dkap_dlnrho_at_face)
               val = interp_val_to_pt(s% d_opacity_dlnd,k,nz,s% dq,'p_dkap_dlnrho_at_face')
            case (p_dkap_dlnt_at_face)
               val = interp_val_to_pt(s% d_opacity_dlnT,k,nz,s% dq,'p_dkap_dlnt_at_face')
            case (p_eps_nuc)
               val = s% eps_nuc(k)
            case (p_log_abs_eps_nuc)
               val = safe_log10_cr(abs(s% eps_nuc(k)))
            case (p_d_epsnuc_dlnd)
               val = s% d_epsnuc_dlnd(k)
            case (p_d_lnepsnuc_dlnd)
               val = s% d_epsnuc_dlnd(k)/max(1d0,abs(s% eps_nuc(k)))
            case (p_d_epsnuc_dlnT)
               val = s% d_epsnuc_dlnT(k)
            case (p_d_lnepsnuc_dlnT)
               val = s% d_epsnuc_dlnT(k)/max(1d0,abs(s% eps_nuc(k)))

            case (p_deps_dlnd_at_face)
               val = interp_val_to_pt(s% d_epsnuc_dlnd,k,nz,s% dq,'p_deps_dlnd_at_face')
            case (p_deps_dlnT_at_face)
               val = interp_val_to_pt(s% d_epsnuc_dlnT,k,nz,s% dq,'p_deps_dlnT_at_face')
            case (p_eps_nuc_neu_total)
               val = s% eps_nuc_neu_total(k)
            case (p_non_nuc_neu)
               val = s% non_nuc_neu(k)
            case (p_nonnucneu_plas)
               val = s% nonnucneu_plas(k)
            case (p_nonnucneu_brem)
               val = s% nonnucneu_brem(k)
            case (p_nonnucneu_phot)
               val = s% nonnucneu_phot(k)
            case (p_nonnucneu_pair)
               val = s% nonnucneu_pair(k)
            case (p_nonnucneu_reco)
               val = s% nonnucneu_reco(k)
            case (p_d_nonnucneu_dlnd)
               val = s% d_nonnucneu_dlnd(k)
            case (p_d_nonnucneu_dlnT)
               val = s% d_nonnucneu_dlnT(k)
               
            case (p_log_irradiation_heat)
               val = safe_log10_cr(s% irradiation_heat(k))
            case (p_cgrav_factor)
               val = s% cgrav(k)/standard_cgrav
               
            case (p_extra_jdot)
               val = s% extra_jdot(k)
            case (p_extra_omegadot)
               val = s% extra_omegadot(k)
            case (p_extra_heat)
               val = s% extra_heat(k)
            case (p_extra_dPdm)
               val = s% extra_dPdm(k)
            case (p_extra_L)
               val = dot_product(s% dm(k:s% nz),s% extra_heat(k:s% nz))/Lsun
            case (p_log_extra_L)
               val = safe_log10_cr( &
                  dot_product(s% dm(k:s% nz),s% extra_heat(k:s% nz))/Lsun)

            case (p_log_abs_eps_grav_dm_div_L)
               val = safe_log10_cr( &
                  abs(s% eps_grav(k))*s% dm(k)/max(1d0,abs(s% L(k))))

            case (p_eps_grav_composition_term)
               if (s% include_dmu_dt_in_eps_grav) &
                  val = s% eps_grav_composition_term(k)
               
            case (p_area)
               val = s% area(k)
            case (p_area_times_v)
               if (s% v_flag) val = s% area(k)*s% v(k)
            case (p_div_v)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     Ap1 = 4*pi*s% R_center*s% R_center
                  else
                     vp1 = s% v(k+1)
                     Ap1 = s% area(k+1)
                  end if
                  val = (s% area(k)*s% v(k) - Ap1*vp1)*s% rho(k)/s% dm(k)
               end if
                                             
            case (p_d_v_div_r_dm)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     rp1 = s% R_center
                  else
                     vp1 = s% v(k+1)
                     rp1 = s% r(k+1)
                  end if
                  v00 = s% v(k)
                  r00 = s% r(k)
                  if (rp1 > 0) then
                     val = (v00/r00 - vp1/rp1)/s% dm(k)
                  end if
               end if
                                             
            case (p_d_v_div_r_dr)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     rp1 = s% R_center
                  else
                     vp1 = s% v(k+1)
                     rp1 = s% r(k+1)
                  end if
                  v00 = s% v(k)
                  r00 = s% r(k)
                  if (rp1 > 0) then
                     val = 4*pi*s% rmid(k)*s% rmid(k)*s% rho(k)* &
                           (v00/r00 - vp1/rp1)/s% dm(k)
                  end if
               end if
               
            case (p_dm_div_4pi_rmid2_rho_dr)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     rp1 = s% R_center
                  else
                     vp1 = s% v(k+1)
                     rp1 = s% r(k+1)
                  end if
                  v00 = s% v(k)
                  r00 = s% r(k)
                  val = s% dm(k)/(4*pi*s% rmid(k)*s% rmid(k)*s% rho(k)*(r00 - rp1))
               end if
            
            case (p_rmid3_rho)
               val = s% rmid(k)*s% rmid(k)*s% rmid(k)*s% rho(k)
            case (p_log_rmid3_rho)
               val = safe_log10_cr(s% rmid(k)*s% rmid(k)*s% rmid(k)*s% rho(k))
            case (p_eta_visc)
               if (s% v_flag .and. s% use_artificial_viscosity) val = s% eta_visc(k)
            case (p_dP_dm_visc)
               if (s% v_flag .and. s% use_artificial_viscosity) val = s% dP_dm_visc(k)
            case (p_eps_visc)
               if (s% v_flag .and. s% use_artificial_viscosity) val = s% eps_visc(k)
            case (p_qvisc)
               if (s% v_flag .and. s% use_artificial_viscosity) val = s% Qvisc(k)

            case (p_dvdt_visc)
               if (s% v_flag .and. s% use_artificial_viscosity) val = s% dvdt_visc(k)
            case (p_dvdt_grav)
               val = -s% cgrav(k)*s% m(k)/(s% r(k)*s% r(k))
            case (p_dvdt_dPdm)
               if (k > 1) val = -4*pi*s% r(k)*s% r(k)*(s% P(k-1) - s% P(k))/s% dm_bar(k)
               
               
            case (p_dm_eps_grav)
               val = s% eps_grav(k)*s% dm(k)
            case (p_eps_grav)
               val = s% eps_grav(k)
            case (p_eps_grav_dm_term_const_q)
               val = s% eps_grav_dm_term_const_q(k)
            case (p_eps_grav_dt_term_const_q)
               val = s% eps_grav_dt_term_const_q(k)
               
            case (p_env_eps_grav)
               val = -s% gradT_sub_grada(k)*s% grav(k)*s% mstar_dot*s% Cp(k)*s% T(k) / &
                        (4*pi*s% r(k)*s% r(k)*s% P(k))
                        
            case (p_mlt_mixing_type)
               int_val = s% mlt_mixing_type(k)
               val = dble(int_val)
               int_flag = .true.
            case (p_mlt_mixing_length)
               val = s% mlt_mixing_length(k)
            case (p_mlt_Gamma)
               val = s% mlt_Gamma(k)
            case (p_mlt_Zeta)
               if (abs(s% gradr(k) - s% grada_at_face(k)) > 1d-20) &
                  val = (s% gradr(k) - s% gradT(k))/(s% gradr(k) - s% grada_at_face(k))
            case (p_adjust_mlt_gradT_fraction)
               val = s% adjust_mlt_gradT_fraction(k)

            case (p_burn_num_iters)
               if (s% split_mixing_choice /= 0) then
                  int_val = s% burn_num_iters(k)
               else
                  int_val = 0
               end if
               val = dble(int_val)
               int_flag = .true.
            case (p_burn_num_substeps)
               if (s% split_mixing_choice /= 0) then
                  int_val = s% burn_num_substeps(k)
               else
                  int_val = 0
               end if
               val = dble(int_val)
               int_flag = .true.

            case (p_num_steps)
               val = dble(s% num_steps(k))
               int_val = s% num_steps(k)
               int_flag = .true.
            case (p_mtx_solve)
               if (s% split_mixing_choice /= 0) then
                  val = dble(s% mtx_solve(k))
                  int_val = s% mtx_solve(k)
               else
                  int_val = 0
               end if
               int_flag = .true.
            case (p_mtx_factor)
               val = dble(s% mtx_factor(k))
               int_val = s% mtx_factor(k)
               int_flag = .true.
            case (p_avg_order)
               val = s% avg_order(k)
               
            case (p_grad_density)
               val = s% grad_density(k)
            case (p_grad_temperature)
               val = s% grad_temperature(k)
               
            case (p_gradL_sub_gradr)
               val = s% gradL(k) - s% gradr(k)
            case (p_grada_sub_gradr)
               val = s% grada_at_face(k) - s% gradr(k)
               
            case (p_gradL)
               val = s% gradL(k)
            case (p_sch_stable)
               if (s% grada(k) > s% gradr(k)) val = 1
            case (p_ledoux_stable)
               if (s% gradL(k) > s% gradr(k)) val = 1

            case (p_eps_nuc_start)
               val = s% eps_nuc_start(k)

            case (p_dominant_isoA_for_thermohaline)
               int_val = chem_isos% Z_plus_N(s% dominant_iso_for_thermohaline(k))
               int_flag = .true.
            case (p_dominant_isoZ_for_thermohaline)
               int_val = chem_isos% Z(s% dominant_iso_for_thermohaline(k))
               int_flag = .true.
            case (p_gradL_composition_term)
               val = s% gradL_composition_term(k)
               
            case (p_log_D_conv)
               if (s% mixing_type(k) == convective_mixing) then
                  val = safe_log10_cr(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_soft)
               if (s% mixing_type(k) == softened_convective_mixing) then
                  val = safe_log10_cr(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_semi)
               if (s% mixing_type(k) == semiconvective_mixing) then
                  val = safe_log10_cr(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_ovr)
               if (s% mixing_type(k) == overshoot_mixing) then
                  val = safe_log10_cr(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_anon)
               if (s% mixing_type(k) == anonymous_mixing) then
                  val = safe_log10_cr(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_thrm)
               if (s% mixing_type(k) == thermohaline_mixing) then
                  val = safe_log10_cr(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if

            case (p_log_D_minimum)
               if (s% mixing_type(k) == minimum_mixing) then
                  val = safe_log10_cr(s% D_mix(k))
               else
                  val = -99
               end if

            case (p_log_D_mix_non_rotation)
               val = safe_log10_cr(s% D_mix_non_rotation(k))
            case (p_log_D_mix)
               val = safe_log10_cr(s% D_mix(k))
            case (p_sig_raw_div_sig)
               if (s% sig(k) < 1d0) then
                  val = 1d0
               else
                  val = s% sig_raw(k)/s% sig(k)
               end if
            case (p_log_sig_mix)
               val = safe_log10_cr(s% sig(k))

            case (p_d_ln_sig_dlnd00)
               if (s% sig(k) > 0) val = s% d_sig_dlnd00(k)/s% sig(k)
            case (p_d_ln_sig_dlnT00)
               if (s% sig(k) > 0) val = s% d_sig_dlnT00(k)/s% sig(k)
            case (p_d_ln_sig_dlndm1)
               if (s% sig(k) > 0) val = s% d_sig_dlndm1(k)/s% sig(k)
            case (p_d_ln_sig_dlnTm1)
               if (s% sig(k) > 0) val = s% d_sig_dlnTm1(k)/s% sig(k)
            case (p_d_ln_sig_dlnR)
               if (s% sig(k) > 0) val = s% d_sig_dlnR(k)/s% sig(k)
            case (p_d_ln_sig_dL)
               if (s% sig(k) > 0) val = s% d_sig_dL(k)/s% sig(k)

            case (p_d_gradT_dlnd00)
               val = s% d_gradT_dlnd00(k)
            case (p_d_gradT_dlnT00)
               val = s% d_gradT_dlnT00(k)
            case (p_d_gradT_dlndm1)
               val = s% d_gradT_dlndm1(k)
            case (p_d_gradT_dlnTm1)
               val = s% d_gradT_dlnTm1(k)
            case (p_d_gradT_dlnR)
               val = s% d_gradT_dlnR(k)
            case (p_d_gradT_dL)
               val = s% d_gradT_dL(k)

            case (p_d_ln_D_mix_dlnd00)
               if (s% D_mix(k) > 0) val = s% d_D_mix_dlnd00(k)/s% D_mix(k)
            case (p_d_ln_D_mix_dlnT00)
               if (s% D_mix(k) > 0) val = s% d_D_mix_dlnT00(k)/s% D_mix(k)
            case (p_d_ln_D_mix_dlndm1)
               if (s% D_mix(k) > 0) val = s% d_D_mix_dlndm1(k)/s% D_mix(k)
            case (p_d_ln_D_mix_dlnTm1)
               if (s% D_mix(k) > 0) val = s% d_D_mix_dlnTm1(k)/s% D_mix(k)
            case (p_d_ln_D_mix_dlnR)
               if (s% D_mix(k) > 0) val = s% d_D_mix_dlnR(k)/s% D_mix(k)
            case (p_d_ln_D_mix_dL)
               if (s% D_mix(k) > 0) val = s% d_D_mix_dL(k)/s% D_mix(k)

            case (p_d_ln_mlt_D_dlnd00)
               if (s% mlt_D(k) > 0) val = s% d_mlt_D_dlnd00(k)/s% mlt_D(k)
            case (p_d_ln_mlt_D_dlnT00)
               if (s% mlt_D(k) > 0) val = s% d_mlt_D_dlnT00(k)/s% mlt_D(k)
            case (p_d_ln_mlt_D_dlndm1)
               if (s% mlt_D(k) > 0) val = s% d_mlt_D_dlndm1(k)/s% mlt_D(k)
            case (p_d_ln_mlt_D_dlnTm1)
               if (s% mlt_D(k) > 0) val = s% d_mlt_D_dlnTm1(k)/s% mlt_D(k)
            case (p_d_ln_mlt_D_dlnR)
               if (s% mlt_D(k) > 0) val = s% d_mlt_D_dlnR(k)/s% mlt_D(k)
            case (p_d_ln_mlt_D_dL)
               if (s% mlt_D(k) > 0) val = s% d_mlt_D_dL(k)/s% mlt_D(k)

            case (p_rp_div_GMrho_face)
               P_face = interp_val_to_pt(s% P,k,s% nz,s% dq,'p_rp_div_GMrho_face')
               rho_face = interp_val_to_pt(s% rho,k,s% nz,s% dq,'p_rp_div_GMrho_face')
               val = s% r(1)*P_face/(s% cgrav(k)*s% m(1)*rho_face)
               
            case (p_log_conv_vel_old)
               val = safe_log10_cr(s% conv_vel_old(k))
            case (p_conv_vel)
               val = s% conv_vel(k)
            case (p_log_conv_vel)
               val = safe_log10_cr(s% conv_vel(k))
            case (p_conv_vel_div_L_vel)
               val = s% conv_vel(k)/max(1d0,get_L_vel(k))
            case (p_conv_vel_div_csound)
               val = s% conv_vel(k)/s% csound(k)
            case (p_log_tau_conv_yrs)
               if (s% conv_vel(k) > 1d-99) then
                  val = safe_log10_cr(s% mlt_mixing_length(k)/(4*s% conv_vel(k)*secyer))
               else
                  val = -99
               end if
            case (p_newly_nonconvective)
               if (s% newly_nonconvective(k)) then
                  int_val = 1; val = 1
               else
                  int_val = 0; val = 0
               end if
               int_flag = .true.
            case (p_mixing_type)
               val = dble(s% mixing_type(k))
               int_val = s% mixing_type(k)
               int_flag = .true.
            case (p_conv_mixing_type) ! OBSOLETE
               val = dble(s% mixing_type(k))
               int_val = s% mixing_type(k)
               int_flag = .true.
            case (p_conv_dP_term)
               val = s% conv_dP_term_factor*s% conv_dP_term(k)
            case (p_log_mlt_D_mix)
               val = safe_log10_cr(s% mlt_D(k))
            case (p_dt_div_conv_tau)
               val = min(1d0,s% dt*s% vc(k)/(s% scale_height(k)/2))
            case (p_pressure_scale_height)
               val = s% P(k)/(s% rho(k)*s% grav(k))/Rsun
            case (p_actual_gradT)
               val = s% actual_gradT(k)
            case (p_gradT_sub_actual_gradT)
               val = s% gradT(k) - s% actual_gradT(k)
            case (p_gradT)
               val = s% gradT(k)
            case (p_gradr)
               val = s% gradr(k)
            case (p_grada_sub_actual_gradT)
               val = s% grada_at_face(k) - s% actual_gradT(k)
            case (p_grada_sub_gradT)
               val = s% grada_at_face(k) - s% gradT(k)

            case (p_log_dr_div_scale_height)
               if (k == 1) then
                  val = 2*(s% r(1) - s% rmid(1))
               else
                  val = s% rmid(k-1) - s% rmid(k)
               end if
               val = safe_log10_cr(val/(s% P(k)/(s% rho(k)*s% grav(k))))
               
            case (p_omega)
               val = if_rot(s% omega,k)
            case (p_omega_start)
               val = if_rot(s% omega_start,k)

            case (p_log_omega)
               val = safe_log10_cr(if_rot(s% omega,k))
            case (p_log_j_rot)
               val = safe_log10_cr(if_rot(s% j_rot,k))
            case (p_log_J_inside)
               if (s% rotation_flag) then
                  val = safe_log10_cr(dot_product(s% j_rot(k:s% nz), s% dm(k:s% nz)))
               else
                  val = -99
               end if
            case (p_log_J_div_M53)
               if (s% rotation_flag) then
                  val = safe_log10_cr( &
                     dot_product(s% j_rot(k:s% nz), &
                                 s% dm(k:s% nz))*1d-50/pow_cr(s% m(k)/Msun,5d0/3d0))
               else
                  val = -99
               end if

            case (p_shear)
               if (s% rotation_flag .and. s% omega(k) > 0) val = s% omega_shear(k)
            case (p_log_abs_shear)
               if (s% rotation_flag .and. s% omega(k) > 0) then
                  val = safe_log10_cr(s% omega_shear(k))
               else
                  val = -99
               end if
            case (p_i_rot)
               val = if_rot(s% i_rot,k)
            case (p_j_rot)
               val = if_rot(s% j_rot,k)
            case (p_v_rot)
               val = if_rot(s% omega,k)*s% r(k)*1d-5 ! km/sec
            case (p_fp_rot)
               val = if_rot(s% fp_rot,k)
            case (p_ft_rot)
               val = if_rot(s% ft_rot,k)
            case (p_fp_rot_minus_1)
               val = if_rot(s% fp_rot,k) - 1d0
            case (p_ft_rot_minus_1)
               val = if_rot(s% ft_rot,k) - 1d0
            case (p_ft_rot_div_fp_rot)
               if (s% rotation_flag) val = s% ft_rot(k)/s% fp_rot(k)
               
            case (p_log_am_nu_non_rot)
               val = safe_log10_cr(if_rot(s% am_nu_non_rot,k))
            case (p_log_am_nu_rot)
               val = safe_log10_cr(if_rot(s% am_nu_rot,k))
            case (p_log_am_nu)
               val = safe_log10_cr(if_rot(s% am_nu_rot,k) + if_rot(s% am_nu_non_rot,k))

            case (p_r_polar)
               val = if_rot(s% r_polar,k)/Rsun
            case (p_log_r_polar)
               val = safe_log10_cr(if_rot(s% r_polar,k)/Rsun)
            case (p_r_equatorial)
               val = if_rot(s% r_equatorial,k)/Rsun
            case (p_log_r_equatorial)
               val = safe_log10_cr(if_rot(s% r_equatorial,k)/Rsun)
            case (p_r_e_div_r_p)
               if (s% rotation_flag) val = s% r_equatorial(k)/s% r_polar(k)
            case (p_omega_crit)
               if (s% rotation_flag) val = omega_crit(s,k)
            case (p_omega_div_omega_crit)
               if (s% rotation_flag) then
                  val = omega_crit(s,k)
                  if (val < 1d-50) then
                     val = 0
                  else
                     val = s% omega(k)/val
                  end if
               end if
            
            case (p_e_field)
               if (s% do_element_diffusion) val = s% E_field(k)
            case (p_log_e_field)
               if (s% do_element_diffusion) val = safe_log10_cr(s% E_field(k))
            
            case (p_g_field_element_diffusion)
               if (s% do_element_diffusion) val = s% g_field_element_diffusion(k)
            case (p_log_g_field_element_diffusion)
               if (s% do_element_diffusion) &
                  val = safe_log10_cr(s% g_field_element_diffusion(k))
            
            case (p_e_div_g_element_diffusion)
               if (s% do_element_diffusion) val = s% E_field(k)/s% g_field_element_diffusion(k)
            case (p_log_e_div_g_element_diffusion)
               if (s% do_element_diffusion) &
                  val = safe_log10_cr(s% E_field(k)/s% g_field_element_diffusion(k))

            case (p_richardson_number)
               val = if_rot(s% richardson_number,k)
            case (p_am_domega_dlnR)
               val = if_rot(s% domega_dlnR,k)
               
            case (p_am_log_sig) ! == am_log_sig_omega
               val = safe_log10_cr(if_rot(s% am_sig_omega,k))
            case (p_am_log_sig_omega)
               val = safe_log10_cr(if_rot(s% am_sig_omega,k))
            case (p_am_log_sig_j)
               val = safe_log10_cr(if_rot(s% am_sig_j,k))

            case (p_am_log_nu_omega)
               val = safe_log10_cr(if_rot(s% am_nu_omega,k))
            case (p_am_log_nu_j)
               val = safe_log10_cr(if_rot(s% am_nu_j,k))

            case (p_am_log_nu_rot)
               val = safe_log10_cr(if_rot(s% am_nu_rot,k))
            case (p_am_log_nu_non_rot)
               val = safe_log10_cr(if_rot(s% am_nu_non_rot,k))
                  
            case (p_am_log_D_visc)
               val = safe_log10_cr(prefact*s% D_visc_factor*if_rot(s% D_visc,k))
            case (p_am_log_D_DSI)
               val = safe_log10_cr(prefact*s% D_DSI_factor*if_rot(s% D_DSI,k))
            case (p_am_log_D_SH)
               val = safe_log10_cr(prefact*s% D_SH_factor*if_rot(s% D_SH,k))
            case (p_am_log_D_SSI)
               val = safe_log10_cr(prefact*s% D_SSI_factor*if_rot(s% D_SSI,k))
            case (p_am_log_D_ES)
               val = safe_log10_cr(prefact*s% D_ES_factor*if_rot(s% D_ES,k))
            case (p_am_log_D_GSF)
               val = safe_log10_cr(prefact*s% D_GSF_factor*if_rot(s% D_GSF,k))
            case (p_am_log_D_ST)
               val = safe_log10_cr(prefact*s% D_ST_factor*if_rot(s% D_ST,k))
            case (p_am_log_nu_ST)
               val = safe_log10_cr(prefact*s% am_nu_ST_factor*if_rot(s% nu_ST,k))

            case (p_dynamo_log_B_r)
               val = safe_log10_cr(if_rot(s% dynamo_B_r,k))
            case (p_dynamo_log_B_phi)
               val = safe_log10_cr(if_rot(s% dynamo_B_phi,k))

            case (p_gradr_start)
               val = s% gradr_start(k)
            case (p_grada_at_face_start)
               val = s% grada_at_face_start(k)
            case (p_gradL_start)
               val = s% gradL_start(k)

            case (p_gradr_sub_gradL_start)
               val = s% gradr_start(k) - s% gradL_start(k)
               
            case (p_grada_at_face)
               val = s% grada_at_face(k)
            case (p_gradr_div_grada)
               val = s% gradr(k)/s% grada_at_face(k)
            case (p_gradr_sub_grada)
               val = s% gradr(k) - s% grada_at_face(k)
            case (p_gradT_sub_a)
               val = s% gradT(k) - s% grada_at_face(k)
            case (p_gradT_sub_grada)
               val = s% gradT(k) - s% grada_at_face(k)
            case (p_gradT_div_grada)
               val = s% gradT(k) / s% grada_at_face(k)
            case (p_gradr_sub_gradT)
               val = s% gradr(k) - s% gradT(k)
            case (p_gradT_sub_gradr)
               val = s% gradT(k) - s% gradr(k)

            case (p_gradT_div_gradr)
               if (abs(s% gradr(k)) < 1d-99) then
                  val = 1d0
               else
                  val = s% gradT(k) / s% gradr(k)
               end if
            case (p_log_gradT_div_gradr)
               if (abs(s% gradr(k)) < 1d-99) then
                  val = 0d0
               else
                  val = safe_log10_cr(s% gradT(k) / s% gradr(k))
               end if
               
            case (p_log_mlt_Gamma)
               val = safe_log10_cr(s% mlt_Gamma(k))

            case (p_mlt_mixing_length_start)
               val = s% mlt_mixing_length_start(k)
            case (p_mlt_mixing_type_start)
               int_val = s% mlt_mixing_type_start(k)
               val = dble(int_val)
               int_flag = .true.
            case (p_mlt_D_start)
               val = s% mlt_D_start(k)
            case (p_mlt_vc_start)
               val = s% mlt_vc_start(k)
            case (p_mlt_Gamma_start)
               val = s% mlt_Gamma_start(k)
            case (p_mlt_cdc_start)
               val = s% mlt_cdc_start(k)
          
            case (p_delta_r)
               val = s% r(k) - s% r_start(k)
            case (p_delta_L)
               val = s% L(k) - s% L_start(k)
            case (p_delta_cell_vol)
               if (k == s% nz) then
                  rp1 = s% R_center
                  rp1_start = s% R_center_old
               else
                  rp1 = s% r(k+1)
                  rp1_start = s% r_start(k+1)
               end if
               r00 = s% r(k)
               r00_start = s% r_start(k)
               dr3 = r00*r00*r00 - rp1*rp1*rp1
               dr3_start = r00_start*r00_start*r00_start - rp1_start*rp1_start*rp1_start
               val = 4d0/3d0*pi*(dr3 - dr3_start)
            case (p_delta_entropy)
               val = s% entropy(k) - s% entropy_start(k)
            case (p_delta_T)
               val = s% T(k) - s% T_start(k)
            case (p_delta_rho)
               val = s% rho(k) - exp_cr(s% lnd_start(k))
            case (p_delta_eps_nuc)
               val = s% eps_nuc(k) - s% eps_nuc_start(k)
            case (p_delta_mu)
               val = s% mu(k) - s% mu_start(k)

            case (p_dlnR_dk)
               if (k >= 1) val = s% lnR(k) - s% lnR(k+1)

            case (p_super_ad)
               val = max(0d0, s% gradT(k) - s% grada_at_face(k))
               
            case (p_dlnP_dm_visc)
               val = s% dlnP_dm_visc(k)
            case (p_accel_div_grav)
               if (s% v_flag) val = s% dv_dt(k)/s% grav(k)
                              
            case (p_log_maxabs_dt_d_dxdt_dt)
               val = maxval(abs(s% d_dxdt_dx(1:s% species,1:s% species,k)))
               val = log10_cr(s% dt*val)
                              
            case (p_dlnPgas_dt_const_q)
               if (s% lnPgas_flag) val = s% dlnPgas_dt_const_q(k)
            case (p_dlnd_dt_const_q)
               if (.not. s% lnPgas_flag) val = s% dlnd_dt_const_q(k)
            case (p_dlnT_dt_const_q)
               val = s% dlnT_dt_const_q(k)
                              
            case (p_dlnPgas)
               if (s% lnPgas_flag) val = s% dt*s% dlnPgas_dt(k)
            case (p_dlnd)
               if (.not. s% lnPgas_flag) val = s% dt*s% dlnd_dt(k)
            case (p_dlnT)
               val = s% dt*s% dlnT_dt(k)
            case (p_dlnR)
               val = s% dt*s% dlnR_dt(k)
                              
            case (p_dlnPgas_dt)
               if (s% lnPgas_flag) val = s% dlnPgas_dt(k)
            case (p_dlnd_dt)
               if (.not. s% lnPgas_flag) val = s% dlnd_dt(k)
            case (p_dlnT_dt)
               val = s% dlnT_dt(k)
            case (p_dlnR_dt)
               val = s% dlnR_dt(k)
            case (p_dr_dt)
               val = s% dlnR_dt(k)*s% r(k)
            case (p_dv_dt)
               if (s% v_flag) val = s% dv_dt(k)

            case (p_del_entropy)
               val = s% entropy(k) - s% entropy_start(k)
            case (p_ds_from_eps_grav)
               val = -s% dt*s% eps_grav(k)/s% T(k)/(avo*kerg)
            
            
            case (p_energy_check)
               if (k == 1 .or. k == s% nz) then
                  val = 0d0
               else
                  val = &
                     -0.5d0*s% dm_bar(k)*(s% v(k)*s% v(k) - s% v_start(k)*s% v_start(k)) &
                     -s% cgrav(k)*s% m(k)*s% dm_bar(k)*(1d0/s% r_start(k) - 1d0/s% r(k)) &
                     -4*pi*s% R2(k)*s% vc(k)*(s% P(k-1) - s% P(k))*s% dt &
                     +4*pi*(s% Qvisc(k-1) - s% Qvisc(k))*s% vc(k)/s% r(k)*s% dt &
                     -(s% energy(k) - s% energy_start(k))*s% dm(k) &
                     -s% dt*(s% L(k) - s% L(k+1)) &
                     -s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1)) &
                     +s% dt*4*pi*s% Qvisc(k)*(s% vc(k)/s% r(k) - s% vc(k+1)/s% r(k+1)) &
                     +s% dt*s% eps_nuc(k)*s% dm(k) &
                     -s% dt*0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))*s% dm(k)
               end if
            
            
            case (p_vcheck)
               if (k == 1) then
                  val = 0d0
               else
                  val = &
                     -0.5d0*s% dm_bar(k)*(s% v(k)*s% v(k) - s% v_start(k)*s% v_start(k)) &
                     -s% cgrav(k)*s% m(k)*s% dm_bar(k)*(1d0/s% r_start(k) - 1d0/s% r(k)) &
                     -4*pi*s% R2(k)*s% vc(k)*(s% P(k-1) - s% P(k))*s% dt &
                     +4*pi*(s% Qvisc(k-1) - s% Qvisc(k))*s% vc(k)/s% r(k)*s% dt     
               end if               
            
            case (p_echeck)
               if (k == s% nz) then
                  val = 0d0
               else
                  val = &
                     -(s% energy(k) - s% energy_start(k))*s% dm(k) &
                     -s% dt*(s% L(k) - s% L(k+1)) &
                     -s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1)) &
                     +s% dt*4*pi*s% Qvisc(k)*(s% vc(k)/s% r(k) - s% vc(k+1)/s% r(k+1)) &
                     +s% dt*s% eps_nuc(k)*s% dm(k) &
                     -s% dt*0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))*s% dm(k)
               end if
            
            
            case (p_dm_eps_grav_plus_eps_visc)
               val = (s% eps_grav(k)+s% eps_visc(k))*s% dm(k)
            case(p_minus_dm_eps_visc)
               val = -s% eps_visc(k)*s% dm(k)
            case(p_dm_eps_visc)
               val = s% eps_visc(k)*s% dm(k)
               
            case(p_dt_dm_eps_grav)
               val = s% eps_grav(k)*s% dm(k)*s% dt
            
            case(p_minus_dm_de_plus_dtPdvA)
               if (k < s% nz .and. s% v_flag) val = &
                  - s% dm(k)*(s% energy(k) - s% energy_start(k)) &
                  - s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1))
                  
            
            case (p_sum_echeck_from_surf)
               kk = min(k,s% nz-1)
               val = sum( &
                  -(s% energy(1:kk) - s% energy_start(1:kk))*s% dm(1:kk) &
                  -s% dt*(s% L(1:kk) - s% L(2:kk+1)) &
                  -s% dt*4*pi*s% P(1:kk)*(s% R2(1:kk)*s% vc(1:kk) - s% R2(2:kk+1)*s% vc(2:kk+1)) &
                  +s% dt*4*pi*s% Qvisc(1:kk)*(s% vc(1:kk)/s% r(1:kk) - s% vc(2:kk+1)/s% r(2:kk+1)) &
                  +s% dt*s% eps_nuc(1:kk)*s% dm(1:kk) &
                  -s% dt*0.5d0*(s% non_nuc_neu_start(1:kk) + s% non_nuc_neu(1:kk))*s% dm(1:kk))
               if (k == s% nz .and. s% r_center > 0) val = val + (&
                  -(s% energy(k) - s% energy_start(k))*s% dm(k) &
                  -s% dt*(s% L(k) - s% L_center) &
                  -s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% r_center*s% r_center*s% v_center) &
                  +s% dt*4*pi*s% Qvisc(k)*(s% vc(k)/s% r(k) - s% v_center/s% r_center) &
                  +s% dt*s% eps_nuc(k)*s% dm(k) &
                  -s% dt*0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))*s% dm(k))
                  
                  
            case(p_sum_check_eps_grav_from_surf)
               kk = min(k,s% nz-1)
               val = sum( &
                  -s% eps_grav(1:kk)*s% dm(1:kk)*s% dt &
                  - s% dm(1:kk)*(s% energy(1:kk) - s% energy_start(1:kk)) &
                  - s% dt*4*pi*s% P(1:kk)*(s% R2(1:kk)*s% vc(1:kk) - s% R2(2:kk+1)*s% vc(2:kk+1)))
               if (k == s% nz .and. s% r_center > 0) val = val + (&
                  -s% eps_grav(k)*s% dm(k)*s% dt &
                  - s% dm(k)*(s% energy(k) - s% energy_start(k)) &
                  - s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% r_center*s% r_center*s% v_center))
                  
            case(p_check_eps_grav)
               if (k < s% nz) val = -s% eps_grav(k)*s% dm(k)*s% dt &
                  - s% dm(k)*(s% energy(k) - s% energy_start(k)) &
                  - s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1))
            
            case(p_dm_de)
               val = s% dm(k)*(s% energy(k) - s% energy_start(k))
            case(p_dt_dL)
               if (k < s% nz) val = s% dt*(s% L(k) - s% L(k+1))
            case(p_dt_PdvA)
               if (k < s% nz) &
                  val = s% dt*4*pi*s% P(k)*(s% R2(k)*s% vc(k) - s% R2(k+1)*s% vc(k+1))

            case (p_lnd_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_lnd)
            case (p_lnT_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_lnT)
            case (p_e_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_E)
            case (p_lnR_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_lnR)
            case (p_lnPgas_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_lnPgas)
            case (p_v_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_v)
            case (p_lum_truncation_ratio)
               val = eval_truncation_ratio_term(s,k,s% i_lum)
            
            case (p_signed_dlnd)
               if (.not. s% lnPgas_flag) then
                  val = 1d6*s% dlnd_dt(k)*s% dt
                  val = sign(1d0,val)*log10_cr(max(1d0,abs(val)))
               end if
            case (p_signed_dlnT)
               val = 1d6*s% dlnT_dt(k)*s% dt
               val = sign(1d0,val)*log10_cr(max(1d0,abs(val)))
            case (p_cno_div_z)
               cno = s% xa(s% net_iso(ic12),k) + &
                     s% xa(s% net_iso(in14),k) + s% xa(s% net_iso(io16),k)
               z = 1 - (s% xa(s% net_iso(ih1),k) + s% xa(s% net_iso(ihe4),k))
               if (z > 1d-50) then
                  val = cno/z
               else
                  val = 0
               end if
            case (p_dE)
               val = s% energy(k) - s% energy_start(k)
            case (p_dr)
               if (k < s% nz) then
                  val = s% r(k) - s% r(k+1)
               else
                  val = s% r(k) - s% R_center
               end if
            case (p_dv)
               if (.not. s% v_flag) then
                  val = 0
               else if (k < s% nz) then
                  val = s% v(k+1) - s% v(k)
               else
                  val = -s% v(k)
               end if
            case (p_dt_dv_div_dr)
               if (.not. s% v_flag) then
                  val = 0
               else if (k < s% nz) then
                  val = s% dt*(s% v(k+1) - s% v(k))/(s% r(k) - s% r(k+1))
               else
                  val = -s% dt*s% v(k)/s% r(k)
               end if

            case (p_dlnH1_dlnP)
               val = get_dlogX_dlogP(ih1, k)
            case (p_dlnHe3_dlnP)
               val = get_dlogX_dlogP(ihe3, k)
            case (p_dlnHe4_dlnP)
               val = get_dlogX_dlogP(ihe4, k)
            case (p_dlnC12_dlnP)
               val = get_dlogX_dlogP(ic12, k)
            case (p_dlnC13_dlnP)
               val = get_dlogX_dlogP(ic13, k)
            case (p_dlnN14_dlnP)
               val = get_dlogX_dlogP(in14, k)
            case (p_dlnO16_dlnP)
               val = get_dlogX_dlogP(io16, k)
            case (p_dlnNe20_dlnP)
               val = get_dlogX_dlogP(ine20, k)
            case (p_dlnMg24_dlnP)
               val = get_dlogX_dlogP(img24, k)
            case (p_dlnSi28_dlnP)
               val = get_dlogX_dlogP(isi28, k)

            case (p_dlog_pp_dlogP)
               val = get_dlog_eps_dlogP(ipp, k)
            case (p_dlog_cno_dlogP)
               val = get_dlog_eps_dlogP(icno, k)
            case (p_dlog_3alf_dlogP)
               val = get_dlog_eps_dlogP(i3alf, k)

            case (p_dlog_burn_c_dlogP)
               val = get_dlog_eps_dlogP(i_burn_c, k)
            case (p_dlog_burn_n_dlogP)
               val = get_dlog_eps_dlogP(i_burn_n, k)
            case (p_dlog_burn_o_dlogP)
               val = get_dlog_eps_dlogP(i_burn_o, k)

            case (p_dlog_burn_ne_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ne, k)
            case (p_dlog_burn_na_dlogP)
               val = get_dlog_eps_dlogP(i_burn_na, k)
            case (p_dlog_burn_mg_dlogP)
               val = get_dlog_eps_dlogP(i_burn_mg, k)

            case (p_dlog_cc_dlogP)
               val = get_dlog_eps_dlogP(icc, k)
            case (p_dlog_co_dlogP)
               val = get_dlog_eps_dlogP(ico, k)
            case (p_dlog_oo_dlogP)
               val = get_dlog_eps_dlogP(ioo, k)

            case (p_dlog_burn_si_dlogP)
               val = get_dlog_eps_dlogP(i_burn_si, k)
            case (p_dlog_burn_s_dlogP)
               val = get_dlog_eps_dlogP(i_burn_s, k)
            case (p_dlog_burn_ar_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ar, k)
            case (p_dlog_burn_ca_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ca, k)
            case (p_dlog_burn_ti_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ti, k)
            case (p_dlog_burn_cr_dlogP)
               val = get_dlog_eps_dlogP(i_burn_cr, k)
            case (p_dlog_burn_fe_dlogP)
               val = get_dlog_eps_dlogP(i_burn_fe, k)
            case (p_dlog_pnhe4_dlogP)
               val = get_dlog_eps_dlogP(ipnhe4, k)
            case (p_dlog_photo_dlogP)
               val = get_dlog_eps_dlogP(iphoto, k)
            case (p_dlog_other_dlogP)
               val = get_dlog_eps_dlogP(iother, k)
               
            case (p_total_energy) ! specific total energy at k
               if (k == 1) then
                  val = s% energy(k)
               else if (k == s% nz) then
                  val = (s% dm(k)*s% energy(k) + &
                              0.5d0*s% dm(k-1)*s% energy(k-1))/ &
                        (s% dm(k) + 0.5d0*s% dm(k-1))
               else
                  val = (s% dm(k)*s% energy(k) + &
                              s% dm(k-1)*s% energy(k-1))/ &
                        (s% dm(k) + s% dm(k-1))
               end if
               val = val - s% cgrav(k)*s% m_grav(k)/s% r(k) + &
                           0.5d0*s% v(k)*s% v(k)
            case (p_total_energy_integral)
               val = get_total_energy_integral(s,k)                     
                     
            case (p_dlnX_dr)
               klo = max(1,k-1)
               khi = min(nz,k+1)
               val = log_cr(max(1d-99,max(1d-99,s% X(klo))/max(1d-99,s% X(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
            case (p_dlnY_dr)
               klo = max(1,k-1)
               khi = min(nz,k+1)
               val = log_cr(max(1d-99,max(1d-99,s% Y(klo))/max(1d-99,s% Y(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
            case (p_dlnRho_dr)
               klo = max(1,k-1)
               khi = min(nz,k+1)
               val = (s% lnd(klo) - s% lnd(khi))/(s% rmid(klo) - s% rmid(khi))

            case (p_brunt_g_r_div_cs2)
               val = s% grav(k)*s% r(k)/pow2(max(1d-20,s% csound_at_face(k)))
               
            case (p_brunt_B)
               val = s% brunt_B(k)
            case (p_brunt_nonB)
               if (s% calculate_Brunt_N2) val = -s% gradT_sub_grada(k)
            case (p_log_brunt_B)
               val = log_cr(max(1d-99,s% brunt_B(k)))
            case (p_log_brunt_nonB)
               if (s% calculate_Brunt_N2) val = log_cr(max(1d-99,-s% gradT_sub_grada(k)))
               
            case (p_brunt_N2)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)
            case (p_brunt_N2_composition_term)
               if (s% calculate_Brunt_N2) val = s% brunt_N2_composition_term(k)
            case (p_brunt_N2_structure_term)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k) - s% brunt_N2_composition_term(k)
            case (p_log_brunt_N2_composition_term)
               if (s% calculate_Brunt_N2) val = &
                  safe_log10_cr(s% brunt_N2_composition_term(k))
            case (p_log_brunt_N2_structure_term)
               if (s% calculate_Brunt_N2) val = &
                  safe_log10_cr(s% brunt_N2(k) - s% brunt_N2_composition_term(k))
               
            case (p_brunt_A)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)*s% r(k)/s% grav(k)
            case (p_brunt_A_div_x2)
               x = s% r(k)/s% r(1)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)*s% r(k)/s% grav(k)/x/x
            case (p_log_brunt_N2_dimensionless)
               if (s% calculate_Brunt_N2) val = &
                  safe_log10_cr(s% brunt_N2(k)/(3*s% cgrav(1)*s% m_grav(1)/pow3(s% r(1))))
            case (p_brunt_N2_dimensionless)
               if (s% calculate_Brunt_N2) val = &
                  s% brunt_N2(k)/(3*s% cgrav(1)*s% m_grav(1)/pow3(s% r(1)))
            case (p_brunt_N_dimensionless)
               if (s% calculate_Brunt_N2) val = &
                  sqrt(max(0d0,s% brunt_N2(k))/(3*s% cgrav(1)*s% m_grav(1)/pow3(s% r(1))))
            case (p_brunt_N)
               if (s% calculate_Brunt_N2) val = sqrt(max(0d0,s% brunt_N2(k)))
            case (p_brunt_frequency) ! cycles per day
               if (s% calculate_Brunt_N2) val = &
                  (24d0*60d0*60d0/(2*pi))*sqrt(max(0d0,s% brunt_N2(k)))
            case (p_log_brunt_N)
               if (s% calculate_Brunt_N2) val = safe_log10_cr(sqrt(max(0d0,s% brunt_N2(k))))
            case (p_log_brunt_N2)
               if (s% calculate_Brunt_N2) val = safe_log10_cr(s% brunt_N2(k))

            case (p_brunt_nu) ! micro Hz
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)
               val = (1d6/(2*pi))*sqrt(max(0d0,val))
            case (p_log_brunt_nu) ! micro Hz
               if (s% calculate_Brunt_N2) &
                  val = safe_log10_cr((1d6/(2*pi))*sqrt(max(0d0,s% brunt_N2(k))))

            case (p_lamb_S)
               val = sqrt(2d0)*s% csound_at_face(k)/s% r(k) ! for l=1
            case (p_lamb_S2)
               val = 2d0*pow2(s% csound_at_face(k)/s% r(k)) ! for l=1               
                              
            case (p_lamb_Sl1)
               val = (1d6/(2*pi))*sqrt(2d0)*s% csound_at_face(k)/s% r(k) ! microHz
            case (p_lamb_Sl2)
               val = (1d6/(2*pi))*sqrt(6d0)*s% csound_at_face(k)/s% r(k) ! microHz
            case (p_lamb_Sl3)
               val = (1d6/(2*pi))*sqrt(12d0)*s% csound_at_face(k)/s% r(k) ! microHz
               
            case (p_log_lamb_Sl1)
               val = safe_log10_cr((1d6/(2*pi))*sqrt(2d0)*s% csound_at_face(k)/s% r(k)) ! microHz
            case (p_log_lamb_Sl2)
               val = safe_log10_cr((1d6/(2*pi))*sqrt(6d0)*s% csound_at_face(k)/s% r(k)) ! microHz
            case (p_log_lamb_Sl3)
               val = safe_log10_cr((1d6/(2*pi))*sqrt(12d0)*s% csound_at_face(k)/s% r(k)) ! microHz
            case (p_log_lamb_Sl10)
               val = safe_log10_cr((1d6/(2*pi))*sqrt(110d0)*s% csound_at_face(k)/s% r(k)) ! microHz
               
            case (p_brunt_N_div_r_integral)
               if (s% calculate_Brunt_N2) val = get_brunt_N_div_r_integral(k)
            case (p_sign_brunt_N2)
               if (s% calculate_Brunt_N2) val = sign(1d0,s% brunt_N2(k))
               
            case (p_k_r_integral)
               if (s% calculate_Brunt_N2) val = get_k_r_integral(k,1,1d0)
               
            case (p_brunt_N2_sub_omega2)
               if (s% calculate_Brunt_N2) then
                  if (.false.) write(*,2) 'profile N2 - omega2', k, &
                     s% brunt_N2(k) - pow2(2*pi*s% nu_max/1d6), &
                     s% brunt_N2(k), pow2(2*pi*s% nu_max/1d6)
                  val = s% brunt_N2(k) - pow2(2*pi*s% nu_max/1d6)
                  if (val > 0d0) then
                     val = 1
                  else
                     val = 0
                  end if
               end if
            case (p_sl2_sub_omega2)
               if (s% calculate_Brunt_N2) then
                  val = 2*pow2(s% csound_at_face(k)/s% r(k)) - pow2(2*pi*s% nu_max/1d6)
                  if (val >= 0d0) then
                     val = 1
                  else
                     val = 0
                  end if
               end if
               
            case (p_cs_at_cell_bdy)
               val = s% csound_at_face(k)
            case (p_log_mdot_cs) ! log10_cr(4 Pi r^2 csound rho / (Msun/year))
               val = safe_log10_cr(4*pi*s% r(k)*s% r(k)*s% csound(k)*s% rho(k)/(Msun/secyer))
            case (p_log_mdot_v) ! log10_cr(4 Pi r^2 v rho / (Msun/year))
               val = safe_log10_cr(4*pi*s% r(k)*s% r(k)*s% v(k)*s% rho(k)/(Msun/secyer))
            case (p_log_L_div_CpTMdot)
               if (s% star_mdot == 0) then
                  val = 0
               else
                  val = safe_log10_cr(s% L(k)/(s% cp(k)*s% T(k)*abs(s% star_mdot)*(Msun/secyer)))
               end if
            case (p_logQ)
               val = s% lnd(k)/ln10 - 2*s% lnT(k)/ln10 + 12
               
            case default
               write(*,*) 'FATAL ERROR in profile columns specifications', c, k
               write(*,*) 'between ' // trim(profile_column_name(c-1)) // ' and ' // &
                  trim(profile_column_name(c+1))
               val = 0
               stop 'profile_getval'
               
         end select
         
         end if
         
         
         contains
         
         
         real(dp) function get_L_vel(k) result(v) ! velocity if L carried by convection
            integer, intent(in) :: k
            real(dp) :: rho_face
            integer :: j
            if (k == 1) then
               j = 2
            else
               j = k
            end if
            rho_face = interp_val_to_pt(s% rho,j,nz,s% dq,'profile get_L_vel')
            v = pow_cr(max(1d0,s% L(k))/(4*pi*s% r(k)*s% r(k)*rho_face),1d0/3d0)
         end function get_L_vel
         

         real(dp) function get_k_r_integral(k_in, el, nu_factor)
            integer, intent(in) :: k_in
            integer, intent(in) :: el
            real(dp), intent(in) :: nu_factor
            real(dp) :: integral, integral_for_k, &
               cs2, r2, n2, sl2, omega2, L2, kr2, dr
            integer :: k, k1, k_inner, k_outer         
            include 'formats'

            if (k_in == 1) then
               get_k_r_integral = 1
               return
            end if
            
            get_k_r_integral = 0
            L2 = el*(el+1)
            omega2 = pow2(1d-6*2*pi*s% nu_max*nu_factor)
            
            ! k_inner and k_outer are bounds of evanescent region
            
            ! k_outer is outermost k where Sl2 <= omega2 at k-1 and Sl2 > omega2 at k
            ! 1st find outermost where Sl2 <= omega2
            k1 = 0
            do k = 2, s% nz
               r2 = s% r(k)*s% r(k)
               cs2 = s% csound_at_face(k)*s% csound_at_face(k)
               sl2 = L2*cs2/r2
               if (sl2 <= omega2) then
                  k1 = k; exit
               end if
            end do
            if (k1 == 0) return
            ! then find next k where Sl2 >= omega2
            k_outer = 0
            do k = k1+1, s% nz
               r2 = s% r(k)*s% r(k)
               cs2 = s% csound_at_face(k)*s% csound_at_face(k)
               sl2 = L2*cs2/r2
               if (sl2 > omega2) then
                  k_outer = k; exit
               end if
            end do
            if (k_outer == 0) return
            if (k_in <= k_outer) then
               get_k_r_integral = 1
               return
            end if
            
            ! k_inner is next k where N2 >= omega2 at k+1 and N2 < omega2 at k
            k_inner = 0
            do k = k_outer+1, s% nz
               if (s% brunt_N2(k) >= omega2) then
                  k_inner= k; exit
               end if
            end do
            if (k_inner == 0) return
            if (k_in > k_inner) then
               get_k_r_integral = 1
               return
            end if

            integral = 0; integral_for_k = 0
            get_k_r_integral = 0
            do k = k_inner, k_outer, -1
               r2 = s% r(k)*s% r(k)
               cs2 = s% csound_at_face(k)*s% csound_at_face(k)
               n2 = s% brunt_N2(k)
               sl2 = L2*cs2/r2
               kr2 = (1 - n2/omega2)*(1 - Sl2/omega2)/cs2
               dr = s% rmid(k-1) - s% rmid(k)
               if (kr2 < 0 .and. omega2 < Sl2) integral = integral + sqrt(-kr2)*dr
               if (k == k_in) integral_for_k = integral
            end do
            if (integral < 1d-99) return
            get_k_r_integral = integral_for_k/integral
            
            if (is_bad_num(get_k_r_integral)) then
               write(*,2) 'get_k_r_integral', k_in, integral_for_k, integral
               stop 'get_k_r_integral'
            end if
            
         end function get_k_r_integral
         

         real(dp) function get_brunt_N_div_r_integral(k_in)
            integer, intent(in) :: k_in
            real(dp) :: integral, integral_for_k, dr
            integer :: k
            integral = 0
            integral_for_k = 0
            get_brunt_N_div_r_integral = 1
            if (k_in == 1) return
            get_brunt_N_div_r_integral = 0
            do k = s% nz, 2, -1
               dr = s% rmid(k-1) - s% rmid(k)
               if (s% brunt_N2(k) > 0) &
                  integral = integral + sqrt(s% brunt_N2(k))*dr/s% r(k)
               if (k == k_in) integral_for_k = integral
            end do
            if (integral < 1d-99) return
            get_brunt_N_div_r_integral = integral_for_k/integral
         end function get_brunt_N_div_r_integral

         
         real(dp) function get_dlogX_dlogP(j, k)
            integer, intent(in) :: j, k
            integer :: ii, i
            real(dp) :: val, x00, xm1, dlogP, dlogX
            include 'formats'
            get_dlogx_dlogp = 0
            if (k > 1) then
               ii = k
            else
               ii = 2
            end if
            i = s% net_iso(j)
            if (i == 0) return
            x00 = s% xa(i,ii)
            xm1 = s% xa(i,ii-1)
            if (x00 < 1d-20 .or. xm1 < 1d-20) return
            dlogP = (s% lnP(ii) - s% lnP(ii-1))/ln10
            if (dlogP <= 0d0) return
            dlogX = log10_cr(x00/xm1)
            get_dlogX_dlogP = dlogX/dlogP
         end function get_dlogX_dlogP

         
         real(dp) function get_dlog_eps_dlogP(cat, k)
            integer, intent(in) :: cat, k
            integer :: ii
            real(dp) :: val, eps, epsm1, dlogP, dlog_eps
            get_dlog_eps_dlogP = 0
            if (k > 1) then
               ii = k
            else
               ii = 2
            end if
            eps = s% eps_nuc_categories(cat,ii)
            epsm1 = s% eps_nuc_categories(cat,ii-1)
            if (eps < 1d-3 .or. epsm1 < 1d-3) return
            dlogP = (s% lnP(ii) - s% lnP(ii-1))/ln10
            if (dlogP <= 0d0) return
            dlog_eps = log10_cr(eps/epsm1)
            get_dlog_eps_dlogP = dlog_eps/dlogP
         end function get_dlog_eps_dlogP

         
         real(dp) function pt(v,k)
            integer, intent(in) :: k
            real(dp), pointer :: v(:)
            if (k == 1) then
               pt = v(k)
            else
               pt = (v(k)*s% dq(k-1) + v(k-1)*s% dq(k))/(s% dq(k-1) + s% dq(k))
            endif
         end function pt

         
         real(dp) function if_rot(v,k)
            integer, intent(in) :: k
            real(dp), pointer :: v(:)
            if (s% rotation_flag) then
               if_rot = v(k)
            else
               if_rot = 0
            endif
         end function if_rot
               
         real(dp) function get_ion_info(id,k)
            use ionization_lib, only: eval_ionization
            integer, intent(in) :: id, k
            integer :: ierr
            if (ionization_k /= k) then
               ierr = 0
               call eval_ionization( &
                  1d0 - (s% X(k) + s% Y(k)), s% X(k), s% Rho(k), s% lnd(k)/ln10, &
                  s% T(k), s% lnT(k)/ln10, ionization_res, ierr)
               if (ierr /= 0) ionization_res = 0
               ionization_k = k
            end if
            get_ion_info = ionization_res(id)
         end function get_ion_info
         
      end subroutine getval_for_profile

      
      real(dp) function eval_truncation_ratio_term(s,k,i) result(err)
         type (star_info), pointer :: s
         integer, intent(in) :: k, i
         if (i == 0 .or. .not. associated(s% xh_compare)) then
            err = 0
         else
            err = abs(s% xh(i,k) - s% xh_compare(i,k)) / &
               (s% truncation_ratio_xh_atol + &
                s% truncation_ratio_xh_rtol* &
                  max(abs(s% xh(i,k)), abs(s% xh_compare(i,k))))
         end if
      end function eval_truncation_ratio_term

      end module profile_getval
      
