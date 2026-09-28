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

      module rates_def
      use utils_def, only: integer_dict
      use const_def, only: dp
      use chem_def, only: iso_name_length, nuclide_data, npart
      
      use rates_def_mic
      
      
      implicit none

      ! weaklib
      
      character (len=256) :: weak_data_dir         

      ! ecapture

      character (len=1000) :: ecapture_states_file
      character (len=1000) :: ecapture_transitions_file
      
      ! reaclib


      integer, parameter :: max_nreaclib=100000
      integer, parameter :: max_species_per_reaction=6
      integer, parameter :: ncoefficients=7
      integer, parameter :: nchapters=11
      integer, parameter :: ninverse_coeff = 2
      integer, parameter :: max_terms_per_rate = 20
      integer, parameter :: max_id_length = 36


      ! i/o parameters
      integer, parameter :: internal_format = 0
      integer, parameter :: pretty_print_format = 1
      integer, parameter :: short_format = 2


      ! flags for reaction types -- values are chapter numbers
      integer, parameter :: &
         r_one_one = 1, &
         r_one_two = 2, &
         r_one_three = 3, &
         r_two_one = 4, &
         r_two_two = 5, &
         r_two_three = 6, &
         r_two_four = 7, &
         r_three_one = 8, &
         r_three_two = 9, &
         r_four_two = 10, &
         r_one_four = 11


      integer, dimension(nchapters) :: Nin = (/1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 1/)
      integer, dimension(nchapters) :: Nout = (/1, 2, 3, 1, 2, 3, 4, 1, 2, 2, 4/)


      type reaclib_data
      	integer, dimension(:), pointer :: chapter
      	character(len=iso_name_length), dimension(:,:), pointer :: species
      	character(len=iso_name_length), dimension(:), pointer :: label
      	character, dimension(:), pointer :: reaction_flag
      	character, dimension(:), pointer :: reverse_flag
      	real(dp), dimension(:), pointer :: Qvalue
      	real(dp), dimension(:,:), pointer :: coefficients
         ! following are for 1D allocation of 2D arrays
      	character(len=iso_name_length), dimension(:), pointer :: species1
      	real(dp), dimension(:), pointer :: coefficients1
      end type reaclib_data


      type reaction_data
      
      	integer :: nreactions
      	integer :: nchapters_included
      	integer, dimension(nchapters) :: chapters_present
      	type(nuclide_data), pointer :: nuclides
      	
      	integer :: num_from_weaklib
      	integer, dimension(:), pointer :: weaklib_ids ! (num_from_weaklib)
      	logical, dimension(:), pointer :: also_in_reaclib ! (num_from_weaklib)
      	
      	integer, dimension(2, nchapters) :: bookmarks
      	character(len=max_id_length), dimension(:), pointer :: reaction_handle, reverse_handle
      	integer, dimension(:), pointer :: category
      	integer, dimension(:), pointer :: chapter
      	integer, dimension(:, :), pointer :: pspecies
      	character, dimension(:), pointer :: reaction_flag
      	real(dp), dimension(:), pointer :: weight
      	real(dp), dimension(:), pointer :: weight_reverse
      	real(dp), dimension(:, :), pointer :: coefficients
      	real(dp), dimension(:), pointer :: weak_mask
      	real(dp), dimension(:, :), pointer :: inverse_coefficients
      	integer, dimension(:), pointer :: inverse_exp
      	real(dp), dimension(:, :), pointer :: inverse_part
         real(dp), dimension(:), pointer :: Q
         real(dp), dimension(:), pointer :: Qneu
         type (integer_dict), pointer :: reaction_dict, reverse_dict
         
      end type reaction_data
      
      
      character (len=256) :: reaclib_dir, reaclib_filename




      ! choices for various rates
         ! NOTE: if change these, must edit raw_rates to match.
         
         ! NACRE = Angulo et al. 1999 Nucl. Phys. A, 656, 3
         ! JR = jina reaclib -- (Sakharuk et al. 2006)
         ! CF88 = Frank Timmes' version of 
            ! Caughlin, G. R. & Fowler, W. A. 1988, Atom. Data and Nuc. Data Tables, 40, 283
         
         
         ! when possible, NACRE is 1, jina reaclib is 2
         integer, parameter :: rates_NACRE_if_available = 1
         integer, parameter :: rates_JR_if_available = 2
         
         ! the jina reaclib rates are not valid below T = 10^7.
         ! if we have the option, we automatically blend over to nacre for low temperatures.
         ! the following values determine the blend region.
         real(dp) :: JR_T_full_off = 1.0d7 ! don't use JR below this
         real(dp) :: JR_T_full_on = 1.1d7 ! don't need to blend above this
         
                  
         ! triple alpha
         integer, parameter :: use_rate_3a_NACRE = 1
         integer, parameter :: use_rate_3a_JR = 2 
         integer, parameter :: use_rate_3a_CF88 = 3
         integer, parameter :: use_rate_3a_FL87 = 4 ! Fushiki and Lamb, Apj, 317, 368-388, 1987
            ! note: use_rate_3a_FL87 is a special case. see eval_FL_epsnuc_3alf in rate_lib
         
         ! c12(a,g)o16
         integer, parameter :: use_rate_c12ag_NACRE = 1
         integer, parameter :: use_rate_c12ag_JR = 2 
         integer, parameter :: use_rate_c12ag_Kunz = 3 ! Kunz et al. (2002)
         integer, parameter :: use_rate_c12ag_CF88 = 4
         
         ! c12 + c12
         integer, parameter :: use_rate_1212_CF88_multi = 1
            ! combines the rates for the n, p, and a channels.
            ! using neutron branching from dayras switkowski and woosley 1976
            ! and an estimate of proton branching.
         integer, parameter :: use_rate_1212_CF88_basic = 2
            ! the single rate approximation from CF88
         integer, parameter :: use_rate_1212_G05 = 3 ! Gasques, et al. Phys Review C, 72, 025806 (2005)
            ! NOTE: Gasques option for c12+c12 is implemented in net, not in rates.
      
         ! n14(p,g)o15
         integer, parameter :: use_rate_n14pg_NACRE = 1
         integer, parameter :: use_rate_n14pg_JR = 2
         integer, parameter :: use_rate_n14pg_CF88 = 3
         
         ! o14(a,p)f17
         integer, parameter :: use_rate_o14ap_CF88 = 1 ! (NOT IN NACRE)
         integer, parameter :: use_rate_o14ap_JR = 2
         
         ! o15(a,g)ne19
         integer, parameter :: use_rate_o15ag_CF88 = 1 ! (NOT IN NACRE)
         integer, parameter :: use_rate_o15ag_JR = 2
               
         ! f17(p,g)ne18
         integer, parameter :: use_rate_f17pg_wk82 = 1 ! (NOT IN NACRE)
            ! wiescher and kettner, ap. j., 263, 891 (1982)
         integer, parameter :: use_rate_f17pg_JR = 2
      
         ! f18(p,g)ne19
         integer, parameter :: use_rate_f18pg_wk82 = 1 ! (NOT IN NACRE)
            ! wiescher and kettner, ap. j., 263, 891 (1982)
         integer, parameter :: use_rate_f18pg_JR = 2
      
         ! f18(p,a)o15
         integer, parameter :: use_rate_f18pa_wk82 = 1 ! (NOT IN NACRE) 
            ! wiescher and kettner, ap. j., 263, 891 (1982)
         integer, parameter :: use_rate_f18pa_JR = 2
         
         ! o16 + o16
         integer, parameter :: use_rate_1616_CF88 = 1
         integer, parameter :: use_rate_1616_reaclib = 2
         
         
         

      ! info for rates being evaluated using tables (rate_list.txt)
      type rate_table_info
         logical :: use_rate_table
         logical :: need_to_read
         character (len=132) :: rate_fname
         integer :: nT8s
         real(dp), pointer :: T8s(:) ! (nT8s)
         real(dp), pointer :: f1(:) ! =(4,nT8s)
      end type rate_table_info
      
      type (rate_table_info), pointer :: raw_rates_records(:)
      character (len=164) :: rates_table_dir

      type (integer_dict), pointer :: skip_warnings_dict
   	
      type (reaction_data), target :: reaclib_rates
         
      character (len=1000) :: rates_dir, rates_cache_dir
      


      ! coulomb corrections for weak reactions
      integer :: which_vs_coulomb = 0
      integer :: which_mui_coulomb = 0

      type Coulomb_Info
         real(dp) :: temp
         real(dp) :: den
         real(dp) :: logT
         real(dp) :: logRho
         real(dp) :: theta_e
         real(dp) :: zbar
         real(dp) :: abar
         real(dp) :: z2bar
         real(dp) :: ye
         real(dp) :: z52bar
         real(dp) :: zbar13
         real(dp) :: abari
         real(dp) :: rr
         real(dp) :: tempi
         real(dp) :: dtempi
         real(dp) :: deni
         real(dp) :: pp
         real(dp) :: rs
         real(dp) :: gamma_e
      end type Coulomb_Info

      
      contains
      
      
      subroutine set_rates_cache_dir(rates_cache_dir_in,ierr)
         use const_def, only: mesa_data_dir, mesa_caches_dir
         use utils_lib, only : mkdir
         character (len=*), intent(in) :: rates_cache_dir_in
         integer, intent(out) :: ierr
         ierr = 0    
         rates_dir = trim(mesa_data_dir) // '/rates_data'
         if (len_trim(rates_cache_dir_in) > 0) then
            rates_cache_dir = rates_cache_dir_in
         else if (len_trim(mesa_caches_dir) > 0) then
            rates_cache_dir = trim(mesa_caches_dir) // '/rates_cache'
         else
            rates_cache_dir = trim(rates_dir) // '/cache'
         end if
         call mkdir(rates_cache_dir)
      end subroutine set_rates_cache_dir


      subroutine start_rates_def_init(ierr)
         use utils_lib, only: integer_dict_define
         use crlibm_lib
         integer, intent(out) :: ierr
         
         integer :: i
         ierr = 0    
         call create_skip_warnings_dict(ierr)  
         if (ierr /= 0) return
         nullify(reaction_names_dict)
         do i=1,rates_reaction_id_max
            call integer_dict_define(reaction_names_dict, reaction_Name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: rates_def_init failed in integer_dict_define'
               return
            end if
         end do
         call do_start_rates_def_init(ierr)
         
      end subroutine start_rates_def_init
      
      
      subroutine create_skip_warnings_dict(ierr)
         use utils_lib
         use utils_def
         integer, intent(out) :: ierr
         
         integer :: iounit, n, i, t, id, read_int
         character (len=256) :: buffer, string, filename, list_filename
         
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         
         nullify(skip_warnings_dict)

         list_filename = 'skip_warnings.list'
         ! first try the local directory
         filename = trim(list_filename)
         open(unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in rates_data
            filename = trim(rates_dir) // '/' // trim(list_filename)
            ierr = 0
            open(unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               call free_iounit(iounit)
               write(*,*) 'failed to open file ' // trim(list_filename)
               return
            end if
         end if
         
         n = 0
         i = 0
         
      reaction_loop: do
            t = token(iounit, n, i, buffer, string)
            if (t == eof_token) exit
            if (t /= name_token) then
               call error; return
            end if
            call integer_dict_define(skip_warnings_dict, string, 1, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: create_skip_warnings_dict failed in integer_dict_define'
               return
            end if
         end do reaction_loop
         
         close(iounit)
         call free_iounit(iounit)

         call integer_dict_create_hash(skip_warnings_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: create_skip_warnings_dict failed'
            return
         end if

         contains
         
#ifdef offload
      !dir$ attributes offload: mic :: error
#endif         
         subroutine error
            ierr = -1
            close(iounit)
            call free_iounit(iounit)
         end subroutine error

      end subroutine create_skip_warnings_dict
      
      
      integer function reaclib_index(handle) result(indx)
      	use utils_lib, only: integer_dict_lookup
      	character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         integer :: ierr
         ierr = 0
      	call integer_dict_lookup(reaclib_rates% reaction_dict, handle, indx, ierr)
      	if (ierr /= 0) indx = 0
      end function reaclib_index
      
      
      integer function reaclib_reverse(handle) result(indx)
      	use utils_lib, only: integer_dict_lookup
      	character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         integer :: ierr
         ierr = 0
      	call integer_dict_lookup(reaclib_rates% reverse_dict, handle, indx, ierr)
      	if (ierr /= 0) indx = 0
      end function reaclib_reverse
   
   
      subroutine weaklib_init(ierr)
         use const_def, only: mesa_data_dir
         integer, intent(out) :: ierr 
         integer :: i
         ierr = 0
         weak_data_dir = trim(mesa_data_dir) // '/rates_data'
         nullify(weak_reactions_dict)
         nullify(weak_info_list_dict)
      end subroutine weaklib_init
      
      
      subroutine free_weak_info
         use utils_lib, only: integer_dict_free
         deallocate( &
            weak_reactions_data, weak_lhs_nuclide_id, weak_rhs_nuclide_id, &
            weak_lhs_nuclide_name, weak_rhs_nuclide_name, weak_reaclib_id)
         call integer_dict_free(weak_reactions_dict)
      end subroutine free_weak_info


      subroutine ecapture_init(ierr)
         use const_def, only: mesa_data_dir
         integer, intent(out) :: ierr 
         ierr = 0
         nullify(ecapture_reactions_dict)
         nullify(ecapture_transitions_number_dict)
         nullify(ecapture_transitions_offset_dict)
         nullify(ecapture_states_number_dict)
         nullify(ecapture_states_offset_dict)
      end subroutine ecapture_init


      subroutine free_ecapture_info
         use utils_lib, only: integer_dict_free
         deallocate(ecapture_transitions_data, ecapture_states_data, ecapture_logft_data, &
                ecapture_nuclide_id, ecapture_lhs_nuclide_id, ecapture_rhs_nuclide_id, &
                ecapture_nuclide_name, ecapture_lhs_nuclide_name, ecapture_rhs_nuclide_name)
         call integer_dict_free(ecapture_reactions_dict)
         call integer_dict_free(ecapture_transitions_number_dict)
         call integer_dict_free(ecapture_transitions_offset_dict)
         call integer_dict_free(ecapture_states_number_dict)
         call integer_dict_free(ecapture_states_offset_dict)
      end subroutine free_ecapture_info
      
      
      subroutine reaclib_init(jina_reaclib_filename)
         use const_def, only: mesa_data_dir
         character (len=*), intent(in) :: jina_reaclib_filename
         reaclib_dir = trim(mesa_data_dir) // '/rates_data'
         !reaclib_filename = 'jina_reaclib_results05301331'
         reaclib_filename = jina_reaclib_filename
         if (len_trim(reaclib_filename) == 0) &
            reaclib_filename = 'jina_reaclib_results_20130213default2'
      end subroutine reaclib_init
      
      
      subroutine allocate_reaclib_data(r, n, ierr)
      	type(reaclib_data), intent(inout) :: r
      	integer, intent(in) :: n
      	integer, intent(out) :: ierr
      	ierr = 0
      	allocate( &
      	   r% chapter(n), r% species1(max_species_per_reaction*n), &
      	   r% label(n), r% reaction_flag(n), r% reverse_flag(n), &
      		r% Qvalue(n), r% coefficients1(ncoefficients*n), stat=ierr)
      	r% species(1:max_species_per_reaction,1:n) => r% species1(1:max_species_per_reaction*n)
      	r% coefficients(1:ncoefficients,1:n) => r% coefficients1(1:ncoefficients*n)
      end subroutine allocate_reaclib_data
      

      subroutine free_reaclib_data(reaclib, ierr)
      	type(reaclib_data), intent(inout) :: reaclib
      	integer, intent(out) :: ierr
      	ierr = 0
      	if (associated(reaclib% chapter)) & 
      		deallocate( &
      		   reaclib% chapter, reaclib% species1, reaclib% label, reaclib% reaction_flag, &
      		   reaclib% reverse_flag, reaclib% Qvalue, reaclib% coefficients1, stat=ierr)
      end subroutine free_reaclib_data
      

      subroutine allocate_reaction_data(r, n, nweak, ierr)
      	type(reaction_data), intent(out) :: r
      	integer, intent(in) :: n ! number of rates
      	integer, intent(in) :: nweak ! number of weaklib rates
      	integer, intent(out) :: ierr
      	allocate( &
      	   r% reaction_handle(n), r% reverse_handle(n), r% category(n), r% chapter(n), &
      	   r% weaklib_ids(nweak), r% also_in_reaclib(nweak), &
      	   r% pspecies(max_species_per_reaction, n), r% reaction_flag(n), &
      		r% weight(n), r% weight_reverse(n), r% coefficients(ncoefficients, n), &
      		r% weak_mask(n), r% inverse_coefficients(ninverse_coeff, n), &
      		r% inverse_exp(n), r% inverse_part(npart, n), r% Q(n), r% Qneu(n), stat=ierr)
      	nullify(r% reaction_dict)
      	nullify(r% reverse_dict)
      end subroutine allocate_reaction_data
      

      subroutine free_reaction_data(r)
         use utils_lib, only: integer_dict_free
      	type(reaction_data), intent(inout) :: r
      	if (associated(r% chapter)) then
      		deallocate( &
      		   r% reaction_handle, r% reverse_handle, r% category, r% chapter, &
      		   r% weaklib_ids, r% also_in_reaclib, r% pspecies, &
      		   r% reaction_flag, r% weight, r% weight_reverse, r% coefficients, r% weak_mask, &
      		   r% inverse_coefficients, r% inverse_exp, r% inverse_part, r% Q, r% Qneu)
      		if (associated(r% reaction_dict)) call integer_dict_free(r% reaction_dict)
      		if (associated(r% reverse_dict)) call integer_dict_free(r% reverse_dict)
      	end if
      end subroutine free_reaction_data
      

      end module rates_def

