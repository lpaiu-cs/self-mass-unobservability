! ***********************************************************************
!
!  Copyright (C) 2010  Bill Paxton
!
!  MESA is free software; you can use it and/or modify
!  it under the combined terms and restrictions of the MESA MANIFESTO
!  and the GNU General Library Public License as published
!  by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  You should have received a copy of the MESA MANIFESTO along with
!  this software; if not, it is available at the mesa website:
!  http://mesa.sourceforge.net/
!
!  MESA is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!  See the GNU Library General Public License for more details.
!
!  You should have received a copy of the GNU Library General Public License
!  along with this software; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module rates_def_mic
      use const_def, only: dp
      use utils_def, only: integer_dict
      use chem_def, only: iso_name_length
      
      
      
      ! routines that run on the mic should use rates_def_mic rather than rates_def.
      
      ! everything defined in rates_def_mic must be present on the MIC as well as the host.
      
      implicit none

#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif


         ! reactions information
         
         integer, parameter :: maxlen_reaction_Name = 32
         character (len=maxlen_reaction_Name), pointer :: reaction_Name(:) ! (rates_reaction_id_max)
            
         integer, parameter :: maxlen_reaction_Info = 72
         character (len=maxlen_reaction_Info), pointer :: reaction_Info(:) ! (rates_reaction_id_max)
         
         real(dp), pointer :: std_reaction_Qs(:) ! (rates_reaction_id_max) 
            ! set at initialization; read-only afterwards.
            ! avg energy including neutrinos
         
         real(dp), pointer :: std_reaction_neuQs(:) ! (rates_reaction_id_max) 
            ! set at initialization; read-only afterwards.
            ! avg neutrino loss
         
         real(dp), pointer :: weak_lowT_rate(:) ! (rates_reaction_id_max) 
            ! these are from reaclib or weak_info.list
            ! set at initialization; read-only afterwards.
            
         integer, pointer :: reaction_screening_info(:,:) !(3,rates_reaction_id_max)
            ! reaction_screen_info(1:2,i) = [chem_id1, chem_id2] for screening.  0's if no screening.
            
         integer, pointer :: weak_reaction_info(:,:) ! (2,rates_reaction_id_max)
            ! weak_reaction_info(1:2,i) = [chem_id_in, chem_id_out].  0's if not a weak reaction.

         integer, pointer :: reaction_ye_rho_exponents(:,:) ! (2,rates_reaction_id_max)
            ! multiply T dependent rate by Ye^a(i) * Rho^b(i)
            ! reaction_ye_rho_coeffs(1,i) is a(i)
            ! reaction_ye_rho_coeffs(2,i) is b(i)
            ! (0,0) for photodisintegrations and decays
            ! (0,1) for standard 2 body reactions
            ! (0,2) for 3 body reactions such as triple alpha
            ! (1,1) for 2 body electron captures
            ! (1,2) for 3 body electron captures (e.g., pep)      
      
         integer, parameter :: max_num_reaction_inputs = 3
         integer, pointer :: reaction_inputs(:,:) ! (2*max_num_reaction_inputs,rates_reaction_id_max)
            ! up to max_num_reaction_inputs pairs of coefficients and chem id's, terminated by 0's.
            ! e.g.,  o16(p,g)f17 would be (/ 1, io16, 1, ih1, 0 /)
            ! triple alpha would be (/ 3, ihe4, 0 /)
            ! he3(he4, g)be7(e-,nu)li7(p,a)he4 would be (/ 1, ihe3, 1, ihe4, i, ih1, 0 /)

         integer, parameter :: max_num_reaction_outputs = 4
         integer, pointer :: reaction_outputs(:,:) ! (2*max_num_reaction_outputs,rates_reaction_id_max)
            ! up to max_num_reaction_outputs pairs of coefficients and chem id's, terminated by 0's.
            ! e.g.,  o16(p,g)f17 would be (/ 1, if17, 0 /)
            ! c12(a, p)n15 would be (/ 1, in15, 1, ih1, 0 /)
         
      ! weak_info_list
      
         integer :: num_weak_info_list_reactions
         real(dp), pointer :: weak_info_list_halflife(:), weak_info_list_Qneu(:)
         type (integer_dict), pointer :: weak_info_list_dict

      ! weaklib
      
         integer, parameter :: weak_num_T9 = 12, weak_num_lYeRho = 11
         real :: weak_reaction_T9s(weak_num_T9) = &
            (/ 0.01, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0 /)
         real :: weak_reaction_lYeRhos(weak_num_lYeRho) = &
            (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 /)
            ! NOTE: code assumes these are evenly spaced.

         real(dp) :: &
            T9_weaklib_full_off = 0.01d0, & ! use pure reaclib for T <= this
            T9_weaklib_full_on = 0.02d0 ! use pure weaklib for T >= this
            ! blend for intermediate temperatures

         ! for high Z elements, switch to reaclib at temp where no longer fully ionized
         ! as rough approximation for this, we switch at Fe to higher values of T9
         integer :: weaklib_blend_hi_Z = 26
         ! if input element has Z >= weaklib_blend_hi_Z, then use the following T9 limits
         real(dp) :: &
            T9_weaklib_full_off_hi_Z = 0.063d0, & ! use pure reaclib for T <= this
            T9_weaklib_full_on_hi_Z = 0.073d0 ! use pure weaklib for T >= this

         integer, parameter :: i_ldecay = 1, i_lcapture = 2, i_lneutrino = 3
         integer :: num_weak_reactions
         real, pointer :: weak_reactions_data1(:) 
         real, pointer :: weak_reactions_data(:,:,:,:,:) 
            ! (4, weak_num_T9, weak_num_lYeRho, 3, num_weak_reactions)

         integer, pointer, dimension(:) :: & ! (num_weak_reactions)
            weak_lhs_nuclide_id, weak_rhs_nuclide_id, weak_reaclib_id
         character(len=iso_name_length), dimension(:), pointer :: &
            weak_lhs_nuclide_name, weak_rhs_nuclide_name ! (num_weak_reactions)
         type (integer_dict), pointer :: weak_reactions_dict

         logical :: weak_bicubic = .false.  
            ! true means do bicubic splines for interpolation
            ! false means just do bilinear
            ! bilinear is safe; bicubic can overshoot near jumps
      
         
      ! ecapture
      
      logical :: do_ecapture = .false.

      type (integer_dict), pointer :: ecapture_states_number_dict
      type (integer_dict), pointer :: ecapture_states_offset_dict

      type (integer_dict), pointer :: ecapture_transitions_number_dict
      type (integer_dict), pointer :: ecapture_transitions_offset_dict

      integer, parameter :: max_num_ecapture_nuclei = 10
      integer, parameter :: max_num_ecapture_states = 10

      integer, parameter :: max_num_ecapture_reactions = 10
      integer, parameter :: max_num_ecapture_transitions = 25

      integer, parameter :: num_transitions_data = 2
      integer, parameter :: i_Si = 1, i_Sf = 2
      integer :: num_ecapture_transitions
      integer, pointer :: ecapture_transitions_data(:,:)
      real(dp), pointer :: ecapture_logft_data(:)

      integer, parameter :: num_states_data = 2
      integer, parameter :: i_E = 1, i_J = 2
      integer :: num_ecapture_states
      real(dp), pointer :: ecapture_states_data(:,:)

      integer :: num_ecapture_nuclei, num_ecapture_reactions
      integer, pointer, dimension(:) :: ecapture_nuclide_id, &
             ecapture_lhs_nuclide_id, ecapture_rhs_nuclide_id ! (num_ecapture_reactions)
      character(len=iso_name_length), dimension(:), pointer :: ecapture_nuclide_name, &
             ecapture_lhs_nuclide_name, ecapture_rhs_nuclide_name ! (num_ecapture_reactions)
      type (integer_dict), pointer :: ecapture_reactions_dict
      

      integer, pointer :: reaction_categories(:) ! (rates_reaction_id_max) set by net using reactions.list info
      

      integer, pointer,dimension(:) :: &
         reaction_is_reverse, reaction_reaclib_lo, reaction_reaclib_hi, reverse_reaction_id
         ! caches for get_reaclib_rate_and_dlnT (in raw_rates.f)
         ! for all of these, 0 means "cache entry not yet set -- don't have the information"



         
      ! for tabular evaluation of the raw reaction rates
      real(dp) :: rattab_thi != 10.301029995664d0 ! log10(highest temp = 2e10)
      real(dp) :: rattab_tlo != 5.30102999566398d0 ! log10(lowest temp = 2e5)
      real(dp) :: rattab_temp_hi != 10**rattab_thi
      real(dp) :: rattab_temp_lo != 10**rattab_tlo
      
      integer :: rattab_points_per_decade = 2000
      integer :: nrattab ! number of reaction rate table temperatures
         ! nrattab = <points per decade>*(rattab_thi - rattab_tlo) + 1
               
      real(dp) :: rattab_tstp != (rattab_thi-rattab_tlo)/(nrattab-1)! step size
               
      ! reactions for hardwired nets and reactions with multiple choices for rates
         integer, parameter :: ir1212 = 1
         integer, parameter :: ir1216 = ir1212+1
         integer, parameter :: ir1216_to_mg24 = ir1216+1
         integer, parameter :: ir1216_to_si28 = ir1216_to_mg24+1
         integer, parameter :: ir1616 = ir1216_to_si28+1
         integer, parameter :: ir1616a = ir1616+1
         integer, parameter :: ir1616g = ir1616a+1
         integer, parameter :: ir1616p_aux = ir1616g+1
         integer, parameter :: ir1616ppa = ir1616p_aux+1
         integer, parameter :: ir1616ppg = ir1616ppa+1
         integer, parameter :: ir_he3_ag_be7 = ir1616ppg+1
         integer, parameter :: ir34_pp2 = ir_he3_ag_be7+1
         integer, parameter :: ir34_pp3 = ir34_pp2+1
         integer, parameter :: ir_al27_pa_mg24 = ir34_pp3+1
         integer, parameter :: ir_ar36_ag_ca40 = ir_al27_pa_mg24+1
         integer, parameter :: ir_ar36_ga_s32 = ir_ar36_ag_ca40+1
         integer, parameter :: ir_b8_gp_be7 = ir_ar36_ga_s32+1
         integer, parameter :: ir_be7_pg_b8 = ir_b8_gp_be7+1
         integer, parameter :: ir_c12_ag_o16 = ir_be7_pg_b8+1
         integer, parameter :: ir_c12_ap_n15 = ir_c12_ag_o16+1
         integer, parameter :: ir_c12_pg_n13 = ir_c12_ap_n15+1
         integer, parameter :: ir_c12_to_he4_he4_he4 = ir_c12_pg_n13+1
         integer, parameter :: ir_c13_an_o16 = ir_c12_to_he4_he4_he4+1
         integer, parameter :: ir_c13_pg_n14 = ir_c13_an_o16+1
         integer, parameter :: ir_ca40_ag_ti44 = ir_c13_pg_n14+1
         integer, parameter :: ir_ca40_ga_ar36 = ir_ca40_ag_ti44+1
         integer, parameter :: ir_cr48_ag_fe52 = ir_ca40_ga_ar36+1
         integer, parameter :: ir_cr48_ga_ti44 = ir_cr48_ag_fe52+1
         integer, parameter :: ir_f17_ap_ne20 = ir_cr48_ga_ti44+1
         integer, parameter :: ir_f17_gp_o16 = ir_f17_ap_ne20+1
         integer, parameter :: ir_f17_pa_o14 = ir_f17_gp_o16+1
         integer, parameter :: ir_f18_gp_o17 = ir_f17_pa_o14+1
         integer, parameter :: ir_f18_pa_o15 = ir_f18_gp_o17+1
         integer, parameter :: ir_f19_ap_ne22 = ir_f18_pa_o15+1
         integer, parameter :: ir_f19_gp_o18 = ir_f19_ap_ne22+1
         integer, parameter :: ir_f19_pa_o16 = ir_f19_gp_o18+1
         integer, parameter :: ir_f19_pg_ne20 = ir_f19_pa_o16+1
         integer, parameter :: ir_fe52_ag_ni56 = ir_f19_pg_ne20+1
         integer, parameter :: ir_fe52_ga_cr48 = ir_fe52_ag_ni56+1
         integer, parameter :: ir_h2_be7_to_h1_he4_he4 = ir_fe52_ga_cr48+1
         integer, parameter :: ir_h2_h2_to_he4 = ir_h2_be7_to_h1_he4_he4+1
         integer, parameter :: ir_h2_he3_to_h1_he4 = ir_h2_h2_to_he4+1

         integer, parameter :: ir_he3_be7_to_h1_h1_he4_he4 = ir_h2_he3_to_h1_he4+1
         integer, parameter :: ir_he3_he3_to_h1_h1_he4 = ir_he3_be7_to_h1_h1_he4_he4+1
         integer, parameter :: ir_he4_he4_he4_to_c12 = ir_he3_he3_to_h1_h1_he4+1
         integer, parameter :: ir_li7_pa_he4 = ir_he4_he4_he4_to_c12+1
         integer, parameter :: ir_mg24_ag_si28 = ir_li7_pa_he4+1
         integer, parameter :: ir_mg24_ap_al27 = ir_mg24_ag_si28+1
         integer, parameter :: ir_mg24_ga_ne20 = ir_mg24_ap_al27+1
         integer, parameter :: ir_n13_ap_o16 = ir_mg24_ga_ne20+1
         integer, parameter :: ir_n13_gp_c12 = ir_n13_ap_o16+1
         integer, parameter :: ir_n13_pg_o14 = ir_n13_gp_c12+1
         integer, parameter :: ir_n14_ag_f18 = ir_n13_pg_o14+1
         integer, parameter :: ir_n14_ap_o17 = ir_n14_ag_f18+1
         integer, parameter :: ir_n14_gp_c13 = ir_n14_ap_o17+1
         integer, parameter :: ir_n14_pg_o15 = ir_n14_gp_c13+1
         integer, parameter :: ir_n15_ag_f19 = ir_n14_pg_o15+1
         integer, parameter :: ir_n15_ap_o18 = ir_n15_ag_f19+1
         integer, parameter :: ir_n15_pa_c12 = ir_n15_ap_o18+1
         integer, parameter :: ir_n15_pg_o16 = ir_n15_pa_c12+1
         integer, parameter :: ir_na23_pa_ne20 = ir_n15_pg_o16+1
         integer, parameter :: ir_ne18_gp_f17 = ir_na23_pa_ne20+1
         integer, parameter :: ir_ne19_ga_o15 = ir_ne18_gp_f17+1
         integer, parameter :: ir_ne19_gp_f18 = ir_ne19_ga_o15+1
         integer, parameter :: ir_ne20_ag_mg24 = ir_ne19_gp_f18+1
         integer, parameter :: ir_ne20_ap_na23 = ir_ne20_ag_mg24+1
         integer, parameter :: ir_ne20_ga_o16 = ir_ne20_ap_na23+1
         integer, parameter :: ir_ne20_gp_f19 = ir_ne20_ga_o16+1
         integer, parameter :: ir_ne22_ag_mg26 = ir_ne20_gp_f19+1
         integer, parameter :: ir_ne22_pg_na23 = ir_ne22_ag_mg26+1
         integer, parameter :: ir_ni56_ga_fe52 = ir_ne22_pg_na23+1
         integer, parameter :: ir_o14_ag_ne18 = ir_ni56_ga_fe52+1
         integer, parameter :: ir_o14_ap_f17 = ir_o14_ag_ne18+1
         integer, parameter :: ir_o14_gp_n13 = ir_o14_ap_f17+1
         integer, parameter :: ir_o15_ag_ne19 = ir_o14_gp_n13+1
         integer, parameter :: ir_o15_ap_f18 = ir_o15_ag_ne19+1
         integer, parameter :: ir_o15_gp_n14 = ir_o15_ap_f18+1
         integer, parameter :: ir_o16_ag_ne20 = ir_o15_gp_n14+1
         integer, parameter :: ir_o16_ap_f19 = ir_o16_ag_ne20+1
         integer, parameter :: ir_o16_ga_c12 = ir_o16_ap_f19+1
         integer, parameter :: ir_o16_gp_n15 = ir_o16_ga_c12+1
         integer, parameter :: ir_o16_pg_f17 = ir_o16_gp_n15+1
         integer, parameter :: ir_o17_pa_n14 = ir_o16_pg_f17+1
         integer, parameter :: ir_o17_pg_f18 = ir_o17_pa_n14+1
         integer, parameter :: ir_o18_ag_ne22 = ir_o17_pg_f18+1
         integer, parameter :: ir_o18_pa_n15 = ir_o18_ag_ne22+1
         integer, parameter :: ir_o18_pg_f19 = ir_o18_pa_n15+1
         integer, parameter :: ir_s32_ag_ar36 = ir_o18_pg_f19+1
         integer, parameter :: ir_s32_ga_si28 = ir_s32_ag_ar36+1
         integer, parameter :: ir_si28_ag_s32 = ir_s32_ga_si28+1
         integer, parameter :: ir_si28_ga_mg24 = ir_si28_ag_s32+1
         integer, parameter :: ir_ti44_ag_cr48 = ir_si28_ga_mg24+1
         integer, parameter :: ir_ti44_ga_ca40 = ir_ti44_ag_cr48+1
         integer, parameter :: iral27pa_aux = ir_ti44_ga_ca40+1
         integer, parameter :: iral27pg_aux = iral27pa_aux+1
         integer, parameter :: irar36ap_aux = iral27pg_aux+1
         integer, parameter :: irar36ap_to_ca40 = irar36ap_aux+1
         integer, parameter :: irar36gp_aux = irar36ap_to_ca40+1
         integer, parameter :: irar36gp_to_s32 = irar36gp_aux+1
         integer, parameter :: irbe7ec_li7_aux = irar36gp_to_s32+1
         integer, parameter :: irbe7pg_b8_aux = irbe7ec_li7_aux+1
         integer, parameter :: irc12_to_c13 = irbe7pg_b8_aux+1
         integer, parameter :: irc12_to_n14 = irc12_to_c13+1
         integer, parameter :: irc12ap_aux = irc12_to_n14+1
         integer, parameter :: irc12ap_to_o16 = irc12ap_aux+1
         integer, parameter :: irca40ap_aux = irc12ap_to_o16+1
         integer, parameter :: irca40ap_to_ti44 = irca40ap_aux+1
         integer, parameter :: irca40gp_aux = irca40ap_to_ti44+1
         integer, parameter :: irca40gp_to_ar36 = irca40gp_aux+1
         integer, parameter :: ircl35pa_aux = irca40gp_to_ar36+1
         integer, parameter :: ircl35pg_aux = ircl35pa_aux+1
         integer, parameter :: irco55gprot_aux = ircl35pg_aux+1
         integer, parameter :: irco55pg_aux = irco55gprot_aux+1
         integer, parameter :: irco55protg_aux = irco55pg_aux+1
         integer, parameter :: ircr48ap_aux = irco55protg_aux+1
         integer, parameter :: ircr48ap_to_fe52 = ircr48ap_aux+1
         integer, parameter :: ircr48gp_aux = ircr48ap_to_fe52+1
         integer, parameter :: ircr48gp_to_ti44 = ircr48gp_aux+1
         integer, parameter :: irf19pg_aux = ircr48gp_to_ti44+1
         integer, parameter :: irfe52ap_aux = irf19pg_aux+1
         integer, parameter :: irfe52ap_to_ni56 = irfe52ap_aux+1
         integer, parameter :: irfe52aprot_aux = irfe52ap_to_ni56+1
         integer, parameter :: irfe52aprot_to_fe54 = irfe52aprot_aux+1
         integer, parameter :: irfe52aprot_to_ni56 = irfe52aprot_to_fe54+1
         integer, parameter :: irfe52gp_aux = irfe52aprot_to_ni56+1
         integer, parameter :: irfe52gp_to_cr48 = irfe52gp_aux+1
         integer, parameter :: irfe52neut_to_fe54 = irfe52gp_to_cr48+1
         integer, parameter :: irfe52ng_aux = irfe52neut_to_fe54+1
         integer, parameter :: irfe53gn_aux = irfe52ng_aux+1
         integer, parameter :: irfe53ng_aux = irfe53gn_aux+1
         integer, parameter :: irfe54a_to_ni56 = irfe53ng_aux+1
         integer, parameter :: irfe54an_aux = irfe54a_to_ni56+1
         integer, parameter :: irfe54an_to_ni56 = irfe54an_aux+1
         integer, parameter :: irfe54aprot_to_fe56 = irfe54an_to_ni56+1
         integer, parameter :: irfe54g_to_fe52 = irfe54aprot_to_fe56+1
         integer, parameter :: irfe54ng_aux = irfe54g_to_fe52+1
         integer, parameter :: irfe54ng_to_fe56 = irfe54ng_aux+1
         integer, parameter :: irfe54prot_to_fe52 = irfe54ng_to_fe56+1
         integer, parameter :: irfe54prot_to_ni56 = irfe54prot_to_fe52+1
         integer, parameter :: irfe54protg_aux = irfe54prot_to_ni56+1
         integer, parameter :: irfe55gn_aux = irfe54protg_aux+1
         integer, parameter :: irfe55ng_aux = irfe55gn_aux+1
         
         integer, parameter :: irfe56ec_fake_to_mn56 = irfe55ng_aux+1
         integer, parameter :: irfe56ec_fake_to_mn57 = irfe56ec_fake_to_mn56+1
         integer, parameter :: irfe56ec_fake_to_cr56 = irfe56ec_fake_to_mn57+1
         integer, parameter :: irfe56ec_fake_to_cr57 = irfe56ec_fake_to_cr56+1
         integer, parameter :: irfe56ec_fake_to_cr58 = irfe56ec_fake_to_cr57+1
         integer, parameter :: irfe56ec_fake_to_cr59 = irfe56ec_fake_to_cr58+1
         integer, parameter :: irfe56ec_fake_to_cr60 = irfe56ec_fake_to_cr59+1
         integer, parameter :: irfe56ec_fake_to_cr61 = irfe56ec_fake_to_cr60+1
         integer, parameter :: irfe56ec_fake_to_cr62 = irfe56ec_fake_to_cr61+1
         integer, parameter :: irfe56ec_fake_to_cr63 = irfe56ec_fake_to_cr62+1
         integer, parameter :: irfe56ec_fake_to_cr64 = irfe56ec_fake_to_cr63+1
         integer, parameter :: irfe56ec_fake_to_cr65 = irfe56ec_fake_to_cr64+1
         integer, parameter :: irfe56ec_fake_to_cr66 = irfe56ec_fake_to_cr65+1
         
         integer, parameter :: irfe56ee_to_ni56 = irfe56ec_fake_to_cr66+1
         integer, parameter :: irfe56gn_aux = irfe56ee_to_ni56+1
         integer, parameter :: irfe56gn_to_fe54 = irfe56gn_aux+1
         integer, parameter :: irfe56prot_to_fe54 = irfe56gn_to_fe54+1
         integer, parameter :: irh2_protg_aux = irfe56prot_to_fe54+1
         integer, parameter :: irh2g_neut_aux = irh2_protg_aux+1
         integer, parameter :: irhe3_neutg_aux = irh2g_neut_aux+1
         integer, parameter :: irhe3gprot_aux = irhe3_neutg_aux+1
         integer, parameter :: irhe4_breakup = irhe3gprot_aux+1
         integer, parameter :: irhe4_rebuild = irhe4_breakup+1
         integer, parameter :: irhe4g_neut_aux = irhe4_rebuild+1
         integer, parameter :: irk39pa_aux = irhe4g_neut_aux+1
         integer, parameter :: irk39pg_aux = irk39pa_aux+1
         integer, parameter :: irmg24ap_aux = irk39pg_aux+1
         integer, parameter :: irmg24ap_to_si28 = irmg24ap_aux+1
         integer, parameter :: irmg24gp_aux = irmg24ap_to_si28+1
         integer, parameter :: irmg24gp_to_ne20 = irmg24gp_aux+1
         integer, parameter :: irmn51pg_aux = irmg24gp_to_ne20+1
         integer, parameter :: irn14_to_c12 = irmn51pg_aux+1
         integer, parameter :: irn14_to_n15 = irn14_to_c12+1
         integer, parameter :: irn14_to_o16 = irn14_to_n15+1
         integer, parameter :: irn14ag_lite = irn14_to_o16+1
         integer, parameter :: irn14gc12 = irn14ag_lite+1
         integer, parameter :: irn14pg_aux = irn14gc12+1
         integer, parameter :: irn15pa_aux = irn14pg_aux+1
         integer, parameter :: irn15pg_aux = irn15pa_aux+1
         integer, parameter :: irna23pa_aux = irn15pg_aux+1
         integer, parameter :: irna23pg_aux = irna23pa_aux+1
         integer, parameter :: irne18ag_to_mg24 = irna23pg_aux+1
         integer, parameter :: irne18ap_to_mg22 = irne18ag_to_mg24+1
         integer, parameter :: irne18ap_to_mg24 = irne18ap_to_mg22+1
         integer, parameter :: irne19pg_to_mg22 = irne18ap_to_mg24+1
         integer, parameter :: irne19pg_to_mg24 = irne19pg_to_mg22+1
         integer, parameter :: irne20ap_aux = irne19pg_to_mg24+1
         integer, parameter :: irne20ap_to_mg24 = irne20ap_aux+1
         integer, parameter :: irne20gp_aux = irne20ap_to_mg24+1
         integer, parameter :: irne20gp_to_o16 = irne20gp_aux+1
         integer, parameter :: irne20pg_to_mg22 = irne20gp_to_o16+1
         integer, parameter :: irne20pg_to_mg24 = irne20pg_to_mg22+1
         integer, parameter :: irneut_to_prot = irne20pg_to_mg24+1
         integer, parameter :: irni56ec_to_fe54 = irneut_to_prot+1
         integer, parameter :: irni56ec_to_fe56 = irni56ec_to_fe54+1
         integer, parameter :: irni56ec_to_co56 = irni56ec_to_fe56+1
         integer, parameter :: irco56ec_to_fe56 = irni56ec_to_co56+1
         integer, parameter :: irni56gp_aux = irco56ec_to_fe56+1
         integer, parameter :: irni56gp_to_fe52 = irni56gp_aux+1
         integer, parameter :: irni56gprot_aux = irni56gp_to_fe52+1
         integer, parameter :: irni56gprot_to_fe52 = irni56gprot_aux+1
         integer, parameter :: irni56gprot_to_fe54 = irni56gprot_to_fe52+1
         integer, parameter :: irni56ng_to_fe54 = irni56gprot_to_fe54+1
         integer, parameter :: irni57na_aux = irni56ng_to_fe54+1
         integer, parameter :: iro16_to_n14 = irni57na_aux+1
         integer, parameter :: iro16_to_o17 = iro16_to_n14+1
         integer, parameter :: iro16ap_aux = iro16_to_o17+1
         integer, parameter :: iro16ap_to_ne20 = iro16ap_aux+1
         integer, parameter :: iro16gp_aux = iro16ap_to_ne20+1
         integer, parameter :: iro16gp_to_c12 = iro16gp_aux+1
         integer, parameter :: iro17_to_o18 = iro16gp_to_c12+1
         integer, parameter :: irp31pa_aux = iro17_to_o18+1
         integer, parameter :: irp31pg_aux = irp31pa_aux+1
         integer, parameter :: irpep_to_he3 = irp31pg_aux+1
         integer, parameter :: irpp_to_he3 = irpep_to_he3+1
         integer, parameter :: irprot_neutg_aux = irpp_to_he3+1
         integer, parameter :: irprot_to_neut = irprot_neutg_aux+1
         integer, parameter :: irs32ap_aux = irprot_to_neut+1
         integer, parameter :: irs32ap_to_ar36 = irs32ap_aux+1
         integer, parameter :: irs32gp_aux = irs32ap_to_ar36+1
         integer, parameter :: irs32gp_to_si28 = irs32gp_aux+1
         integer, parameter :: irsc43pa_aux = irs32gp_to_si28+1
         integer, parameter :: irsc43pg_aux = irsc43pa_aux+1
         integer, parameter :: irsi28ap_aux = irsc43pg_aux+1
         integer, parameter :: irsi28ap_to_s32 = irsi28ap_aux+1
         integer, parameter :: irsi28gp_aux = irsi28ap_to_s32+1
         integer, parameter :: irsi28gp_to_mg24 = irsi28gp_aux+1
         integer, parameter :: irti44ap_aux = irsi28gp_to_mg24+1
         integer, parameter :: irti44ap_to_cr48 = irti44ap_aux+1
         integer, parameter :: irti44gp_aux = irti44ap_to_cr48+1
         integer, parameter :: irti44gp_to_ca40 = irti44gp_aux+1
         integer, parameter :: irv47pa_aux = irti44gp_to_ca40+1
         integer, parameter :: irv47pg_aux = irv47pa_aux+1
         integer, parameter :: ir_h1_h1_wk_h2 = irv47pg_aux+1
         integer, parameter :: ir_h1_h1_ec_h2 = ir_h1_h1_wk_h2+1
         integer, parameter :: irn14ag_to_ne22 = ir_h1_h1_ec_h2+1
         integer, parameter :: irf19pa_aux = irn14ag_to_ne22+1
         integer, parameter :: ir_b8_wk_he4_he4 = irf19pa_aux+1
         integer, parameter :: irmn51pa_aux = ir_b8_wk_he4_he4+1
         integer, parameter :: irfe54gn_aux = irmn51pa_aux+1
         integer, parameter :: irco55pa_aux = irfe54gn_aux+1
         integer, parameter :: irco55prota_aux = irco55pa_aux+1
         integer, parameter :: irn14ag_to_o18 = irco55prota_aux+1
         integer, parameter :: ir_h1_he3_wk_he4 = irn14ag_to_o18+1
         integer, parameter :: ir_be7_wk_li7 = ir_h1_he3_wk_he4+1
         ! rates added for approx21
         integer, parameter :: ir_al27_pg_si28 = ir_be7_wk_li7+1
         integer, parameter :: ir_si28_gp_al27 = ir_al27_pg_si28+1
         integer, parameter :: ir_si28_ap_p31 = ir_si28_gp_al27+1
         integer, parameter :: ir_p31_pa_si28 = ir_si28_ap_p31+1
         integer, parameter :: ir_p31_pg_s32 = ir_p31_pa_si28+1
         integer, parameter :: ir_s32_gp_p31 = ir_p31_pg_s32+1
         integer, parameter :: ir_s32_ap_cl35 = ir_s32_gp_p31+1
         integer, parameter :: ir_cl35_pa_s32 = ir_s32_ap_cl35+1
         integer, parameter :: ir_cl35_pg_ar36 = ir_cl35_pa_s32+1
         integer, parameter :: ir_ar36_gp_cl35 = ir_cl35_pg_ar36+1
         integer, parameter :: ir_ar36_ap_k39 = ir_ar36_gp_cl35+1
         integer, parameter :: ir_k39_pa_ar36 = ir_ar36_ap_k39+1
         integer, parameter :: ir_k39_pg_ca40 = ir_k39_pa_ar36+1
         integer, parameter :: ir_ca40_gp_k39 = ir_k39_pg_ca40+1
         integer, parameter :: ir_ca40_ap_sc43 = ir_ca40_gp_k39+1
         integer, parameter :: ir_sc43_pa_ca40 = ir_ca40_ap_sc43+1
         integer, parameter :: ir_sc43_pg_ti44 = ir_sc43_pa_ca40+1
         integer, parameter :: ir_ti44_gp_sc43 = ir_sc43_pg_ti44+1
         integer, parameter :: ir_ti44_ap_v47 = ir_ti44_gp_sc43+1
         integer, parameter :: ir_v47_pa_ti44 = ir_ti44_ap_v47+1
         integer, parameter :: ir_v47_pg_cr48 = ir_v47_pa_ti44+1
         integer, parameter :: ir_cr48_gp_v47 = ir_v47_pg_cr48+1
         integer, parameter :: ir_cr48_ap_mn51 = ir_cr48_gp_v47+1
         integer, parameter :: ir_mn51_pa_cr48 = ir_cr48_ap_mn51+1
         integer, parameter :: ir_mn51_pg_fe52 = ir_mn51_pa_cr48+1
         integer, parameter :: ir_fe52_gp_mn51 = ir_mn51_pg_fe52+1
         integer, parameter :: ir_fe52_ap_co55 = ir_fe52_gp_mn51+1
         integer, parameter :: ir_co55_pa_fe52 = ir_fe52_ap_co55+1
         integer, parameter :: ir_co55_pg_ni56 = ir_co55_pa_fe52+1
         integer, parameter :: ir_ni56_gp_co55 = ir_co55_pg_ni56+1
         integer, parameter :: ir_fe52_ng_fe53 = ir_ni56_gp_co55+1
         integer, parameter :: ir_fe53_gn_fe52 = ir_fe52_ng_fe53+1
         integer, parameter :: ir_fe53_ng_fe54 = ir_fe53_gn_fe52+1
         integer, parameter :: ir_fe54_gn_fe53 = ir_fe53_ng_fe54+1
         integer, parameter :: ir_fe54_pg_co55 = ir_fe54_gn_fe53+1
         integer, parameter :: ir_co55_gp_fe54 = ir_fe54_pg_co55+1
         integer, parameter :: ir_he3_ng_he4 = ir_co55_gp_fe54+1
         integer, parameter :: ir_he4_gn_he3 = ir_he3_ng_he4+1
         integer, parameter :: ir_h1_ng_h2 = ir_he4_gn_he3+1
         integer, parameter :: ir_h2_gn_h1 = ir_h1_ng_h2+1
         integer, parameter :: ir_h2_pg_he3 = ir_h2_gn_h1+1
         integer, parameter :: ir_he3_gp_h2 = ir_h2_pg_he3+1
         integer, parameter :: ir_fe54_ng_fe55 = ir_he3_gp_h2+1
         integer, parameter :: ir_fe55_gn_fe54 = ir_fe54_ng_fe55+1
         integer, parameter :: ir_fe55_ng_fe56 = ir_fe55_gn_fe54+1
         integer, parameter :: ir_fe56_gn_fe55 = ir_fe55_ng_fe56+1
         integer, parameter :: ir_fe54_ap_co57 = ir_fe56_gn_fe55+1
         integer, parameter :: ir_co57_pa_fe54 = ir_fe54_ap_co57+1
         integer, parameter :: ir_fe56_pg_co57 = ir_co57_pa_fe54+1
         integer, parameter :: ir_co57_gp_fe56 = ir_fe56_pg_co57+1
         
         integer, parameter :: ir_c12_c12_to_h1_na23 = ir_co57_gp_fe56+1
         integer, parameter :: ir_he4_ne20_to_c12_c12 = ir_c12_c12_to_h1_na23+1
         integer, parameter :: ir_c12_c12_to_he4_ne20 = ir_he4_ne20_to_c12_c12+1
         integer, parameter :: ir_he4_mg24_to_c12_o16 = ir_c12_c12_to_he4_ne20+1

         integer, parameter :: num_predefined_reactions = ir_he4_mg24_to_c12_o16
         integer :: rates_reaction_id_max


      
      ! for mazurek's ni56 electron capture rate interpolation
         real(dp) :: tv(7),rv(6),rfdm(4),rfd0(4),rfd1(4),rfd2(4),tfdm(5),tfd0(5),tfd1(5),tfd2(5)


      type T_Factors
         real(dp) :: lnT9
         real(dp) :: T9
         real(dp) :: T92
         real(dp) :: T93
         real(dp) :: T94
         real(dp) :: T95
         real(dp) :: T96
         real(dp) :: T912
         real(dp) :: T932
         real(dp) :: T952
         real(dp) :: T972
         real(dp) :: T913
         real(dp) :: T923
         real(dp) :: T943
         real(dp) :: T953
         real(dp) :: T973
         real(dp) :: T9113
         real(dp) :: T914
         real(dp) :: T934
         real(dp) :: T954
         real(dp) :: T974
         real(dp) :: T915
         real(dp) :: T935
         real(dp) :: T945
         real(dp) :: T965
         real(dp) :: T917
         real(dp) :: T927
         real(dp) :: T947
         real(dp) :: T918
         real(dp) :: T938
         real(dp) :: T958
         real(dp) :: T9i
         real(dp) :: T9i2
         real(dp) :: T9i3
         real(dp) :: T9i12
         real(dp) :: T9i32
         real(dp) :: T9i52
         real(dp) :: T9i72
         real(dp) :: T9i13
         real(dp) :: T9i23
         real(dp) :: T9i43
         real(dp) :: T9i53
         real(dp) :: T9i14
         real(dp) :: T9i34
         real(dp) :: T9i54
         real(dp) :: T9i15
         real(dp) :: T9i35
         real(dp) :: T9i45
         real(dp) :: T9i65
         real(dp) :: T9i17
         real(dp) :: T9i27
         real(dp) :: T9i47
         real(dp) :: T9i18
         real(dp) :: T9i38
         real(dp) :: T9i58
         real(dp) :: T916
         real(dp) :: T976
         real(dp) :: T9i76
      end type T_Factors
      
      
      
      ! rate results components
      
      integer, parameter :: i_rate = 1        
      integer, parameter :: i_rate_dT = 2  
      integer, parameter :: i_rate_dRho = 3 
      integer, parameter :: num_rvs = 3
      
      
      
      ! screening
      
      integer, parameter :: no_screening = 0
      integer, parameter :: classic_screening = 1
         ! DeWitt, Graboske, Cooper, "Screening Factors for Nuclear Reactions. 
         !    I. General Theory", ApJ, 181:439-456, 1973.
         ! Graboske, DeWitt, Grossman, Cooper, "Screening Factors for Nuclear Reactions. 
         !    II. Intermediate Screening and Astrophysical Applications", ApJ, 181:457-474, 1973.
      integer, parameter :: extended_screening = 2
         ! based on code from Frank Timmes
         ! extends the Graboske method using results from Alastuey and Jancovici (1978),
         ! along with plasma parameters from Itoh et al (1979) for strong screening.
      integer, parameter :: salpeter_screening = 3
         ! weak screening only.  following Salpeter (1954),
         ! with equations (4-215) and (4-221) of Clayton (1968).

      type Screen_Info
         real(dp) :: temp
         real(dp) :: den
         real(dp) :: logT
         real(dp) :: logRho
         real(dp) :: theta_e
         real(dp) :: zbar
         real(dp) :: abar
         real(dp) :: z2bar
         real(dp) :: zbar13
         real(dp) :: zbar0pt28
         real(dp) :: z1pt58bar
         real(dp) :: ztilda ! sqrt(z2bar + zbar*theta_e)  ! (Dewitt eqn 4)
         real(dp) :: ztilda0pt58
         real(dp) :: Lambda0 ! = 1.88d8*sqrt(rho/(abar*T**3)) ! (Graboske eqn 19; mu_I = abar)
         real(dp) :: Lambda0b ! Lambda0**0.86
         real(dp) :: Lambda0_23 ! Lambda0**(2d0/3d0)
         real(dp) :: ytot
         real(dp) :: rr
         real(dp) :: tempi
         real(dp) :: dtempi
         real(dp) :: deni
         real(dp) :: pp
         real(dp) :: dppdt
         real(dp) :: dppdd
         real(dp) :: qlam0z
         real(dp) :: qlam0zdt
         real(dp) :: qlam0zdd
         real(dp) :: taufac
         real(dp) :: taufacdt
         real(dp) :: xni
         real(dp) :: dxnidd
         real(dp) :: aa
         real(dp) :: daadt
         real(dp) :: daadd
         integer :: num_calls, num_cache_hits
      end type Screen_Info
      
      
      real(dp) :: reaclib_min_T9 ! for T9 < this, return 0 for reaclib strong rates
      

      ! integers to 1/3 (for graboske screening)
      integer, parameter :: num_one_thirds = 60
      real(dp) :: one_third_power(num_one_thirds)
      ! integers to 1.86 (for graboske screening)
      integer, parameter :: num_pow_186 = 60
      real(dp) :: pow_186(num_pow_186)


      type (integer_dict), pointer :: reaction_names_dict
      
      logical :: have_finished_initialization = .false.

      
      contains


      subroutine do_start_rates_def_init(ierr)
         use crlibm_lib
         integer, intent(out) :: ierr
         integer :: i         
         ierr = 0
         do i=1,num_one_thirds
            one_third_power(i) = pow_cr(dble(i),1d0/3d0)
         end do
         do i=1,num_pow_186
            pow_186(i) = pow_cr(dble(i),1.86d0)
         end do
         call set_rattab_range(5.30102999566398d0, 10.301029995664d0)
         
         reaclib_min_T9 = 1d-2 
            ! need <= 2d-3 for pre-ms li7 burning
            ! pre-ms deuterium burning needs much lower (4d-4)
            ! but that seems to cause problems during advanced burning.
                        
      end subroutine do_start_rates_def_init
      
      
      subroutine set_rattab_range(tlo, thi)
         use crlibm_lib
         real(dp), intent(in) :: tlo, thi         
         if (abs(thi - tlo) < 1d-6) then
            rattab_tlo = tlo
            rattab_temp_lo = exp10_cr(rattab_tlo)
            rattab_thi = rattab_tlo
            rattab_temp_hi = rattab_temp_lo
            nrattab = 1
            rattab_tstp = 0
            return
         end if
         rattab_thi = thi
         rattab_tlo = tlo
         rattab_temp_hi = exp10_cr(rattab_thi)
         rattab_temp_lo = exp10_cr(rattab_tlo)
         nrattab = rattab_points_per_decade*(rattab_thi - rattab_tlo) + 1
         if (nrattab <= 1) then
            rattab_thi = rattab_tlo
            nrattab = 1
            rattab_tstp = 0
         else
            rattab_tstp = (rattab_thi-rattab_tlo)/(nrattab-1)
         end if
      end subroutine set_rattab_range
      
      
		integer function get_rates_reaction_id(reaction_name) result(value)
		   use utils_lib, only: integer_dict_lookup
			character (len=*), intent(in)  :: reaction_name 
			integer :: ierr
			integer :: indx
			ierr = 0
			call integer_dict_lookup(reaction_names_dict, reaction_name, value, ierr)
			if (ierr /= 0) value = 0
		end function get_rates_reaction_id


      subroutine create_ecapture_dict_key(ecapture_lhs, ecapture_rhs, key)
         character(len=iso_name_length), intent(in) :: ecapture_lhs, ecapture_rhs
         character(len=2*iso_name_length+1), intent(out) :: key
         key = trim(ecapture_lhs) // ' ' // trim(ecapture_rhs)
      end subroutine create_ecapture_dict_key
      
      
      subroutine create_weak_dict_key(weak_lhs, weak_rhs, key)
	      character(len=iso_name_length), intent(in) :: weak_lhs, weak_rhs
	      character(len=2*iso_name_length+1), intent(out) :: key
         key = trim(weak_lhs) // ' ' // trim(weak_rhs)
      end subroutine create_weak_dict_key
      

		integer function do_get_weak_rate_id(lhs, rhs) ! returns 0 if reaction not found
		   use utils_lib
			character (len=*), intent(in)  :: lhs, rhs 
			integer :: ierr, i
	      character(len=2*iso_name_length+1) :: key
			character (len=iso_name_length) :: lhs_name, rhs_name
			ierr = 0
			do_get_weak_rate_id = 0
			lhs_name = adjustl(lhs)
			rhs_name = adjustl(rhs)
			call create_weak_dict_key(lhs_name, rhs_name, key)
			call integer_dict_lookup(weak_reactions_dict, key, i, ierr)
		   if (ierr /= 0) then
		      !write(*,*) 'failed in integer_dict_lookup for key ' // trim(key)
		      return
		   end if
		   do_get_weak_rate_id = i
		end function do_get_weak_rate_id
      

		integer function do_get_weak_info_list_id(lhs, rhs) ! returns 0 if reaction not found
		   ! value can be used to index weak_info_list_halflife and weak_info_list_Qneu
		   use utils_lib
			character (len=*), intent(in)  :: lhs, rhs ! names as in weak_info.list file
			integer :: ierr, i
	      character(len=2*iso_name_length+1) :: key
			character (len=iso_name_length) :: lhs_name, rhs_name
			ierr = 0
			do_get_weak_info_list_id = 0
			lhs_name = adjustl(lhs)
			rhs_name = adjustl(rhs)
			call create_weak_dict_key(lhs_name, rhs_name, key)
			call integer_dict_lookup(weak_info_list_dict, key, i, ierr)
		   if (ierr /= 0) then
		      !write(*,'(a)') 'get_weak_info_list_id failed for ' // trim(key)
		      return
		   end if
		   do_get_weak_info_list_id = i
		end function do_get_weak_info_list_id
		
		
		integer function get_num_reaction_inputs(ir)
		   integer, intent(in) :: ir
		   integer :: j
		   include 'formats.dek'
		   if (max_num_reaction_inputs == 3) then
            if (reaction_inputs(5,ir) /= 0) then
               get_num_reaction_inputs = 3
            else if (reaction_inputs(3,ir) /= 0) then
               get_num_reaction_inputs = 2
            else if (reaction_inputs(1,ir) /= 0) then
               get_num_reaction_inputs = 1
            else 
               get_num_reaction_inputs = 0
            end if
            return
		   end if
         get_num_reaction_inputs = max_num_reaction_inputs
         do j = 1, 2*max_num_reaction_inputs-1, 2
            if (reaction_inputs(j,ir) == 0) then
               get_num_reaction_inputs = (j-1)/2
               exit
            end if
         end do
		end function get_num_reaction_inputs
		
		
		integer function get_num_reaction_outputs(ir)
		   integer, intent(in) :: ir
		   integer :: j
		   if (max_num_reaction_outputs == 3) then
            if (reaction_outputs(5,ir) /= 0) then
               get_num_reaction_outputs = 3
            else if (reaction_outputs(3,ir) /= 0) then
               get_num_reaction_outputs = 2
            else if (reaction_outputs(1,ir) /= 0) then
               get_num_reaction_outputs = 1
            else 
               get_num_reaction_outputs = 0
            end if
            return
		   end if
         get_num_reaction_outputs = max_num_reaction_outputs
         do j = 1, 2*max_num_reaction_outputs-1, 2
            if (reaction_outputs(j,ir) == 0) then
               get_num_reaction_outputs = (j-1)/2
               exit
            end if
         end do
		end function get_num_reaction_outputs
      
      
#ifdef offload
      !dir$ attributes offload: mic :: do_copy_weak_info_to_coprocessor
#endif
      subroutine do_copy_weak_info_to_coprocessor( & ! on mic
            num_info_list_reactions, info_list_halflife, info_list_Qneu, &
            info_list_dict_keys, info_list_dict_values, &
            full_off, full_on, blend_hi_Z, full_off_hi_Z, full_on_hi_Z, &
            num_reactions, reactions_data1, &
            lhs_nuclide_id, rhs_nuclide_id, reaclib_id, &
            lhs_nuclide_name, rhs_nuclide_name, &
            reactions_dict_keys, reactions_dict_values, bicubic, &
            ierr)
         use utils_lib, only: integer_dict_define, integer_dict_create_hash
         use utils_def, only: maxlen_key_string
         integer, intent(in) :: num_info_list_reactions
         real(dp), pointer, intent(in) :: &
            info_list_halflife(:), info_list_Qneu(:)
         character (len=maxlen_key_string), pointer, intent(in) :: &
            info_list_dict_keys(:)
         integer, pointer, intent(in) :: info_list_dict_values(:)
         real(dp), intent(in) :: full_off, full_on
         integer, intent(in) :: blend_hi_Z
         real(dp), intent(in) :: full_off_hi_Z, full_on_hi_Z
         integer, intent(in) :: num_reactions
         real, pointer, intent(in) :: reactions_data1(:)
         integer, pointer, intent(in) :: &
            lhs_nuclide_id(:), rhs_nuclide_id(:), reaclib_id(:)
         character(len=iso_name_length), pointer, intent(in) :: &
            lhs_nuclide_name(:), rhs_nuclide_name(:)
         character (len=maxlen_key_string), pointer, intent(in) :: &
            reactions_dict_keys(:)
         integer, pointer, intent(in) :: reactions_dict_values(:)
         logical, intent(in) :: bicubic
         integer, intent(out) :: ierr
      
         integer :: sz, i
         
         ierr = 0

         num_weak_info_list_reactions = num_info_list_reactions
         sz = size(weak_info_list_halflife,dim=1)
         allocate(weak_info_list_halflife(sz), weak_info_list_Qneu(sz))
         do i=1,sz
            weak_info_list_halflife(i) = info_list_halflife(i)
            weak_info_list_Qneu(i) = info_list_Qneu(i)
         end do
         
         ! copy weak_info_list_dict
         nullify(weak_info_list_dict)
         sz = size(info_list_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(weak_info_list_dict, &
               info_list_dict_keys(i), info_list_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(weak_info_list_dict, ierr)
         if (ierr /= 0) return
         
         T9_weaklib_full_off = full_off
         T9_weaklib_full_on = full_on
         weaklib_blend_hi_Z = blend_hi_Z
         T9_weaklib_full_off_hi_Z = full_off_hi_Z
         T9_weaklib_full_on_hi_Z = full_on_hi_Z
         num_weak_reactions = num_reactions
         sz = size(reactions_data1,dim=1)
         allocate(weak_reactions_data1(sz))
         do i=1,sz
            weak_reactions_data1(i) = reactions_data1(i)
         end do
         weak_reactions_data( &
            1:4, 1:weak_num_T9, 1:weak_num_lYeRho, 1:3, 1:num_weak_reactions) => &
               weak_reactions_data1( &
                  1:4*weak_num_T9*weak_num_lYeRho*3*num_weak_reactions)
         sz = num_weak_reactions
         allocate( &
            weak_lhs_nuclide_id(sz), &
            weak_rhs_nuclide_id(sz), &
            weak_reaclib_id(sz), &
            weak_lhs_nuclide_name(sz), &
            weak_rhs_nuclide_name(sz))
         do i=1,sz
            weak_lhs_nuclide_id(i) = lhs_nuclide_id(i)
            weak_rhs_nuclide_id(i) = rhs_nuclide_id(i)
            weak_reaclib_id(i) = reaclib_id(i)
            weak_lhs_nuclide_name(i) = lhs_nuclide_name(i)
            weak_rhs_nuclide_name(i) = rhs_nuclide_name(i)
         end do
         
         ! copy weak_reactions_dict
         nullify(weak_reactions_dict)
         sz = size(reactions_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(weak_reactions_dict, &
               reactions_dict_keys(i), reactions_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(weak_reactions_dict, ierr)
         if (ierr /= 0) return
         
         weak_bicubic = bicubic
         
         !assume not necessary to copy weak_reaction_T9s or weak_reaction_lYeRhos
      
      end subroutine do_copy_weak_info_to_coprocessor
      
#ifdef offload
      !dir$ attributes offload: mic :: do_copy_ecapture_info_to_coprocessor
#endif
      subroutine do_copy_ecapture_info_to_coprocessor( & ! on mic
            do_ecapture_flag, &
            num_nuclei, num_reactions, &
            nuclide_id, lhs_nuclide_id, rhs_nuclide_id, &
            nuclide_name, lhs_nuclide_name, rhs_nuclide_name, &
            reactions_dict_keys, reactions_dict_values, &
            ecapture_states, states_data, &
            states_number_dict_keys, states_offset_dict_keys, &
            states_number_dict_values, states_offset_dict_values, &
            ecapture_transitions, transitions_data, logft_data, &
            transitions_number_dict_keys, transitions_offset_dict_keys, &
            transitions_number_dict_values, transitions_offset_dict_values, &
            ierr)
         use utils_lib, only: integer_dict_define, integer_dict_create_hash
         use utils_def, only: maxlen_key_string
         logical, intent(in) :: do_ecapture_flag
         integer, intent(in) :: num_nuclei, num_reactions
         integer, pointer, dimension(:), intent(in) :: &
            nuclide_id, lhs_nuclide_id, rhs_nuclide_id
         character(len=iso_name_length), dimension(:), pointer, intent(in) :: &
            nuclide_name, lhs_nuclide_name, rhs_nuclide_name
         character (len=maxlen_key_string), pointer, intent(in) :: &
            reactions_dict_keys(:)
         integer, pointer, intent(in) :: reactions_dict_values(:)
         integer, intent(in) :: ecapture_states
         real(dp), pointer, intent(in) :: states_data(:,:)
         character (len=maxlen_key_string), pointer, intent(in) :: &
            states_number_dict_keys(:), states_offset_dict_keys(:)
         integer, pointer, intent(in) :: &
            states_number_dict_values(:), states_offset_dict_values(:)
         integer, intent(in) :: ecapture_transitions
         integer, pointer, intent(in) :: transitions_data(:,:)
         real(dp), pointer, intent(in) :: logft_data(:)
         character (len=maxlen_key_string), pointer, intent(in) :: &
            transitions_number_dict_keys(:), transitions_offset_dict_keys(:)
         integer, pointer, intent(in) :: &
            transitions_number_dict_values(:), transitions_offset_dict_values(:)
         integer, intent(out) :: ierr

         integer :: sz, i, j, sz1, sz2
         
         ierr = 0

         do_ecapture = do_ecapture_flag
         num_ecapture_nuclei = num_nuclei
         num_ecapture_reactions = num_reactions
         sz = num_ecapture_reactions
         allocate( &
            ecapture_nuclide_id(sz), &
            ecapture_lhs_nuclide_id(sz), &
            ecapture_rhs_nuclide_id(sz), &
            ecapture_nuclide_name(sz), &
            ecapture_lhs_nuclide_name(sz), &
            ecapture_rhs_nuclide_name(sz))
         do i=1,sz
            ecapture_nuclide_id(i) = nuclide_id(i)
            ecapture_lhs_nuclide_id(i) = lhs_nuclide_id(i)
            ecapture_rhs_nuclide_id(i) = rhs_nuclide_id(i)
            ecapture_nuclide_name(i) = nuclide_name(i)
            ecapture_lhs_nuclide_name(i) = lhs_nuclide_name(i)
            ecapture_rhs_nuclide_name(i) = rhs_nuclide_name(i)
         end do
         
         ! ecapture_reactions_dict
         nullify(ecapture_reactions_dict)
         sz = size(reactions_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(ecapture_reactions_dict, &
               reactions_dict_keys(i), reactions_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(ecapture_reactions_dict, ierr)
         if (ierr /= 0) return

         num_ecapture_states = ecapture_states
         sz1 = size(states_data,dim=1)
         sz2 = size(states_data,dim=2)
         allocate(ecapture_states_data(sz1,sz2))
         do j=1,sz2
            do i=1,sz1
               ecapture_states_data(i,j) = states_data(i,j)
            end do
         end do

         ! ecapture_states_number_dict
         nullify(ecapture_states_number_dict)
         sz = size(states_number_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(ecapture_states_number_dict, &
               states_number_dict_keys(i), states_number_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(ecapture_states_number_dict, ierr)
         if (ierr /= 0) return

         ! ecapture_states_offset_dict
         nullify(ecapture_states_offset_dict)
         sz = size(states_offset_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(ecapture_states_offset_dict, &
               states_offset_dict_keys(i), states_offset_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(ecapture_states_offset_dict, ierr)
         if (ierr /= 0) return

         num_ecapture_transitions = ecapture_transitions
         sz1 = size(transitions_data,dim=1)
         sz2 = size(transitions_data,dim=2)
         allocate(ecapture_transitions_data(sz1,sz2))
         do j=1,sz2
            do i=1,sz1
               ecapture_transitions_data(i,j) = transitions_data(i,j)
            end do
         end do
         
         sz = size(logft_data,dim=1)
         allocate(ecapture_logft_data(sz))
         do i=1,sz
            ecapture_logft_data(i) = logft_data(i)
         end do

         ! ecapture_transitions_number_dict
         nullify(ecapture_transitions_number_dict)
         sz = size(transitions_number_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(ecapture_transitions_number_dict, &
               transitions_number_dict_keys(i), transitions_number_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(ecapture_transitions_number_dict, ierr)
         if (ierr /= 0) return

         ! ecapture_transitions_offset_dict
         nullify(ecapture_transitions_offset_dict)
         sz = size(transitions_offset_dict_values,dim=1)
         do i=1,sz
            call integer_dict_define(ecapture_transitions_offset_dict, &
               transitions_offset_dict_keys(i), transitions_offset_dict_values(i), ierr)
            if (ierr /= 0) return
         end do
         call integer_dict_create_hash(ecapture_transitions_offset_dict, ierr)
         if (ierr /= 0) return
         
      end subroutine do_copy_ecapture_info_to_coprocessor     

         
#ifdef offload
      !dir$ attributes offload: mic :: do_copy1_reaction_info_to_coprocessor
#endif         
      subroutine do_copy1_reaction_info_to_coprocessor( & ! runs on mic
            ir, rates_max, &
            reaction_Name_in, &
            reaction_Info_in, &
            reaction_category_in, &
            reaction_is_reverse_in, &
            reaction_reaclib_lo_in, &
            reaction_reaclib_hi_in, &
            reverse_reaction_id_in, &
            std_reaction_Q_in, &
            std_reaction_neuQ_in, &
            weak_lowT_rate_in, &
            reaction_screening_info_in, &
            weak_reaction_info_in, &
            reaction_ye_rho_exponents_in, &
            reaction_inputs_in, &
            reaction_outputs_in, &
            ierr)
         use utils_lib, only: integer_dict_define

         integer, intent(in) :: ir, rates_max
         character (len=maxlen_reaction_Name), intent(in) :: reaction_Name_in
         character (len=maxlen_reaction_Info), intent(in) :: reaction_Info_in
         integer, intent(in) :: reaction_category_in
         integer, intent(in) :: reaction_is_reverse_in
         integer, intent(in) :: reaction_reaclib_lo_in
         integer, intent(in) :: reaction_reaclib_hi_in
         integer, intent(in) :: reverse_reaction_id_in
         integer, intent(in) :: reaction_screening_info_in(3) 
         integer, intent(in) :: weak_reaction_info_in(2) 
         integer, intent(in) :: reaction_ye_rho_exponents_in(2) 
         integer, intent(in) :: reaction_inputs_in(2*max_num_reaction_inputs) 
         integer, intent(in) :: reaction_outputs_in(2*max_num_reaction_outputs) 
         real(dp), intent(in) :: std_reaction_Q_in
         real(dp), intent(in) :: std_reaction_neuQ_in
         real(dp), intent(in) :: weak_lowT_rate_in
         integer, intent(out) :: ierr
         
         integer :: i
         
         include 'formats'
         
         ierr = 0

         
         if (ir > rates_max) then
            write(*,*) 'do_copy1_reaction_info_to_coprocessor failed', ir, rates_max
            stop 1
         end if

         rates_reaction_id_max = rates_max
         reaction_Name(ir) = reaction_Name_in
         reaction_Info(ir) = reaction_Info_in
         reaction_categories(ir) = reaction_category_in
         reaction_is_reverse(ir) = reaction_is_reverse_in
         reaction_reaclib_lo(ir) = reaction_reaclib_lo_in
         reaction_reaclib_hi(ir) = reaction_reaclib_hi_in         
         reverse_reaction_id(ir) = reverse_reaction_id_in         
         std_reaction_Qs(ir) = std_reaction_Q_in
         std_reaction_neuQs(ir) = std_reaction_neuQ_in         
         weak_lowT_rate(ir) = weak_lowT_rate_in
         reaction_screening_info(:,ir) = reaction_screening_info_in(:)
         weak_reaction_info(:,ir) = weak_reaction_info_in(:)
         reaction_ye_rho_exponents(:,ir) = reaction_ye_rho_exponents_in(:)
         reaction_inputs(:,ir) = reaction_inputs_in(:)
         reaction_outputs(:,ir) = reaction_outputs_in(:)

         call integer_dict_define( &
            reaction_names_dict, reaction_Name_in, ir, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_copy1_reaction_info_to_coprocessor failed'
            return
         end if
            
         if (reaction_screening_info(1,2) < 0) stop 1
         
      end subroutine do_copy1_reaction_info_to_coprocessor

#ifdef offload
      !dir$ end options
#endif

      end module rates_def_mic

