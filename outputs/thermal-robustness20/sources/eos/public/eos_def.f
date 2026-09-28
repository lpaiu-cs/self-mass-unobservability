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
!
! ***********************************************************************

      module eos_def
      
      use const_def, only: dp
      
      implicit none

#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif
       
      ! cgs units
      
      ! the basic eos results
      
      integer, parameter :: i_lnPgas = 1
            ! gas pressure (total pressure minus radiation pressure)
      integer, parameter :: i_lnE = 2 
            ! internal energy per gram
      integer, parameter :: i_lnS = 3 
            ! entropy per gram
      integer, parameter :: i_grad_ad = 4 
            ! dlnT_dlnP at constant S
      integer, parameter :: i_chiRho = 5 
            ! dlnP_dlnRho at constant T
      integer, parameter :: i_chiT = 6 
            ! dlnP_dlnT at constant Rho
      integer, parameter :: i_Cp = 7 
            ! dh_dT at constant P, specific heat at constant total pressure
            ! where h is enthalpy, h = E + P/Rho
      integer, parameter :: i_Cv = 8 
            ! dE_dT at constant Rho, specific heat at constant volume
      integer, parameter :: i_dE_dRho = 9 
            ! at constant T
      integer, parameter :: i_dS_dT = 10 
            ! at constant Rho
      integer, parameter :: i_dS_dRho = 11 
            ! at constant T
      integer, parameter :: i_mu = 12 
            ! mean molecular weight per gas particle (ions + free electrons)
      integer, parameter :: i_lnfree_e = 13
            ! free_e := total combined number per nucleon of free electrons and positrons
      integer, parameter :: i_gamma1 = 14 
            ! dlnP_dlnRho at constant S
      integer, parameter :: i_gamma3 = 15 
            ! gamma3 - 1 = dlnT_dlnRho at constant S
      integer, parameter :: i_eta = 16 
            ! electron degeneracy parameter (eta > 1 for significant degeneracy)
            ! eta = ratio of electron chemical potential to kT
      
      integer, parameter :: num_eos_basic_results = 16
      
      
      character (len=20) :: eosDT_result_names(num_eos_basic_results)

      
      ! NOTE: the calculation of eta is based on the following equation for ne, 
      ! the mean number of free electrons per cm^3,
      ! assuming non-relativistic electrons (okay for T < 10^9 at least)
      !
      !  ne = 4 Pi / h^3 (2 me k T)^1.5 F[1/2,eta]   -- see, for example, Clayton, eqn 2-57
      !        where F is the fermi-dirac integral: 
      !        F[beta,eta] := Integrate[(u^beta)/(1+E^(u-eta)),{u,0,Infinity}]
      ! 
      ! CAVEAT: when free_e, the mean number of free electrons per nucleon gets really small, 
      ! eta isn't very interesting because there aren't a lot of free electrons to be degenerate!
      ! our calculation of eta gets flaky at this point as well.
      ! we sweep this problem under the rug by making eta tend to a fairly large negative value 
      ! when free_e < 0.01 or so. this roughly corresponds to T < 10^4 or less.

      
      integer, parameter :: max_num_eos_Zs = 5
      real(dp), parameter :: eos_Zs(max_num_eos_Zs) = (/ 0.00d0, 0.02d0, 0.04d0, 0.20d0, 1.00d0 /)
      integer, parameter :: num_eos_Xs_for_Z(max_num_eos_Zs) = (/ 6, 5, 5, 5, 1 /)

      integer, parameter :: num_eosDT_Zs = 5
      integer, parameter :: num_eosPT_Zs = 5
      integer, parameter :: num_eosDE_Zs = 3
      
      integer, parameter :: num_eos_Xs = 6
      real(dp), parameter :: eos_Xs(num_eos_Xs) =
     >      (/ 0.0d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0, 1.0d0 /)

      integer, parameter :: sz_per_eos_point = 4 ! for bicubic spline interpolation

      type EoS_General_Info
         logical :: include_radiation, always_skip_elec_pos
         real(dp) :: mass_fraction_limit_for_PC ! skip any species with abundance < this
         real(dp) :: logRho1_OPAL_SCVH_limit ! don't use OPAL_SCVH for logRho > this
         real(dp) :: logRho2_OPAL_SCVH_limit ! full OPAL_SCVH okay for logRho < this
         real(dp) :: logRho1_PC_limit ! okay for pure PC for logRho > this
         real(dp) :: logRho2_PC_limit ! don't use PC for logRho < this (>= 2.8)
         ! transition log_Gamma for PC to HELM
         real(dp) :: log_Gamma_all_HELM ! HELM for log_Gamma <= this
         real(dp) :: Gamma_all_HELM ! 10**log_Gamma_all_HELM
         real(dp) :: log_Gamma_all_PC ! PC for log_Gamma >= this
         real(dp) :: PC_min_Z ! don't use PC for Z < this
         ! transition Z for OPAL to HELM
         real(dp) :: Z_all_HELM ! HELM for Z >= this
         ! transition temperature zone for OPAL to HELM at high T
         real(dp) :: logT_all_HELM ! HELM for lgT >= this
         real(dp) :: logT_all_OPAL ! OPAL for lgT <= this
         ! transition temperature zone for SCVH to HELM at very low T
         real(dp) :: logT_low_all_HELM ! HELM for lgT <= this
         real(dp) :: logT_low_all_SCVH ! SCVH for lgT >= this
         ! transition energy zone for OPAL to HELM (for eosDE)
         real(dp) :: logE_all_HELM ! HELM for lgE >= this
         real(dp) :: logE_all_OPAL ! OPAL for lgE <= this
         ! bookkeeping
         integer :: handle
         logical :: in_use
      end type EoS_General_Info

      
      include 'helm_def.dek'


      ! THE FOLLOWING ARE PRIVATE DEFS -- NOT FOR USE BY CLIENTS
      
      
      ! data tables -- read only after they are initialized
      
         type (Helm_Table), pointer :: eos_ht

         type EosDT_XZ_Info
            real :: logQ_min ! logQ = logRho - 2*logT + 12
            real :: logQ_max
            real :: del_logQ ! spacing for the logQs
            integer :: num_logQs
            real :: logT_min
            real :: logT_max
            real :: del_logT ! spacing for the logTs
            integer :: num_logTs
            real, pointer :: logQs(:), logTs(:)
            real, pointer :: tbl1(:) ! =(sz_per_eos_point, num_eos_vals, num_logQs, num_logTs)
            integer :: version
         end type EosDT_XZ_Info

         type EosPT_XZ_Info
            real :: logW_min ! logW = logPgas - 4*logT
            real :: logW_max
            real :: del_logW ! spacing for the logWs
            integer :: num_logWs
            real :: logT_min
            real :: logT_max
            real :: del_logT ! spacing for the logTs
            integer :: num_logTs
            real, pointer :: logWs(:), logTs(:)
            real, pointer :: tbl1(:) ! =(sz_per_eos_point, num_eos_vals, num_logWs, num_logTs)
            integer :: version
         end type EosPT_XZ_Info

         type EosDE_XZ_Info
            real :: logV_min ! logV = logRho - 0.7*logE + 20
            real :: logV_max
            real :: del_logV ! spacing for the logVs
            integer :: num_logVs
            real :: logE_min
            real :: logE_max
            real :: del_logE ! spacing for the logEs
            integer :: num_logEs
            real, pointer :: logVs(:), logEs(:)
            real, pointer :: tbl1(:) ! =(sz_per_eos_point, num_eos_vals, num_logVs, num_logEs)
            integer :: version
         end type EosDE_XZ_Info
      
         type (EosDT_XZ_Info), target :: eosDT_XZ_data(num_eos_Xs, num_eosDT_Zs)
         type (EosPT_XZ_Info), target :: eosPT_XZ_data(num_eos_Xs, num_eosPT_Zs)
         type (EosDE_XZ_Info), target :: eosDE_XZ_data(num_eos_Xs, num_eosDE_Zs)
      

      ! parameters for eosDT

      ! internal storage of table info
      integer, parameter :: eosDT_ilnE = 1
      integer, parameter :: eosDT_ilnPgas = 2
      ! the rest are the same for all of the cases
      integer, parameter :: eos_ilnS = 3
      integer, parameter :: eos_igrad_ad = 4
      integer, parameter :: eos_ichiRho = 5
      integer, parameter :: eos_ichiT = 6
      integer, parameter :: eos_iCp = 7
      integer, parameter :: eos_iCv = 8
      integer, parameter :: eos_idE_dRho = 9
      integer, parameter :: eos_idS_dT = 10
      integer, parameter :: eos_idS_dRho = 11
      integer, parameter :: eos_imu = 12
      integer, parameter :: eos_ilnfree_e = 13
      integer, parameter :: eos_igamma1 = 14
      integer, parameter :: eos_igamma3 = 15
      integer, parameter :: eos_ieta = 16

      integer, parameter :: num_eos_vals = 16

      ! we are storing more than we'd need to store if we had unlimited numerical precision.
      ! for example, Cv, the S derivatives, and the gammas can all be derived from the other stuff.
      ! better to be safe, we'll take the storage hit and do them independently.


      ! parameters for eosPT
      integer, parameter :: eosPT_ilnE = 1
      integer, parameter :: eosPT_ilnRho = 2

      ! parameters for eosDE
      integer, parameter :: eosDE_ilnT = 1
      integer, parameter :: eosDE_ilnPgas = 2


      ! interpolation info for theta_e
      integer :: theta_e_nx
      real(dp), pointer :: f_theta_e1(:), f_theta_e(:,:)
      real(dp), pointer :: x_theta_e(:)
      
      
      integer, parameter :: max_eos_handles = 10
      type (EoS_General_Info), target :: eos_handles(max_eos_handles)
      
      
      logical :: use_cache_for_eos = .true.
      logical :: eos_root_is_initialized = .false.
      logical :: eosDT_is_initialized = .false.
      logical :: eosPT_is_initialized = .false.
      logical :: eosDE_is_initialized = .false.
      
      character(len=128) :: eosDT_file_prefix, eosDT_Z1_suffix
      character(len=128) :: eosPT_file_prefix, eosPT_Z1_suffix
      character(len=128) :: eosDE_file_prefix
      character(len=1000) :: eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir

      
      contains
      
      
      subroutine eos_def_init
         integer :: i
         do i=1,max_eos_handles
            eos_handles(i)% handle = i
            eos_handles(i)% in_use = .false.
         end do
         eos_root_is_initialized = .false.
         eosDT_is_initialized = .false.
         eosPT_is_initialized = .false.
         eosDE_is_initialized = .false.
         use_cache_for_eos = .true.
         eosDT_file_prefix = 'mesa'
         eosPT_file_prefix = 'mesa'
         eosDE_file_prefix = 'mesa'
         eosDT_Z1_suffix = '_CO_1'
         eosPT_Z1_suffix = '_CO_1'
         eosDT_result_names(i_lnPgas) = 'lnPgas'
         eosDT_result_names(i_lnE) = 'lnE'
         eosDT_result_names(i_lnS) = 'lnS'
         eosDT_result_names(i_grad_ad) = 'grad_ad'
         eosDT_result_names(i_chiRho) = 'chiRho'
         eosDT_result_names(i_chiT) = 'chiT'
         eosDT_result_names(i_Cp) = 'Cp'
         eosDT_result_names(i_Cv) = 'Cv'
         eosDT_result_names(i_dE_dRho) = 'dE_dRho'
         eosDT_result_names(i_dS_dT) = 'dS_dT'
         eosDT_result_names(i_dS_dRho) = 'dS_dRho'
         eosDT_result_names(i_mu) = 'mu'
         eosDT_result_names(i_lnfree_e) = 'lnfree_e'
         eosDT_result_names(i_gamma1) = 'gamma1'
         eosDT_result_names(i_gamma3) = 'gamma3'
         eosDT_result_names(i_eta) = 'eta'
      end subroutine eos_def_init

      
      integer function do_alloc_eos(ierr)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do_alloc_eos = -1
!$omp critical (eos_handle)
         do i = 1, max_eos_handles
            if (.not. eos_handles(i)% in_use) then
               eos_handles(i)% in_use = .true.
               do_alloc_eos = i
               exit
            end if
         end do
!$omp end critical (eos_handle)
         if (do_alloc_eos == -1) then
            ierr = -1
            return
         end if
         if (eos_handles(do_alloc_eos)% handle /= do_alloc_eos) then
            ierr = -1
            return
         end if
         call init_eos_handle_data(do_alloc_eos)
#ifdef offload
         !dir$ offload target(mic) in(do_alloc_eos)
         call init_eos_handle_data(do_alloc_eos)
#endif
      end function do_alloc_eos


      subroutine get_result_names(names)
         character (len=8) :: names(num_eos_basic_results)
         names(i_lnPgas) = 'lnPgas'  
         names(i_lnE) = 'lnE' ! internal energy per gram
         names(i_lnS) = 'lnS' ! entropy per gram
         names(i_grad_ad) = 'grad_ad' ! dlnT_dlnP at constant S
         names(i_chiRho) = 'chiRho' ! dlnP_dlnRho at constant T      
         names(i_chiT) = 'chiT' ! dlnP_dlnT at constant Rho      
         names(i_Cp) = 'Cp' ! dE_dT at constant P, specific heat at constant pressure     
         names(i_Cv) = 'Cv' ! dE_dT at constant Rho, specific heat at constant volume
         names(i_dE_dRho) = 'dE_dRho' ! at constant T      
         names(i_dS_dT) = 'dS_dT' ! at constant Rho      
         names(i_dS_dRho) = 'dS_dRho'  ! at constant T      
         names(i_mu) = 'mu'  
         names(i_lnfree_e) = 'lnfree_e'       
         names(i_gamma1) = 'gamma1'  ! dlnP_dlnRho at constant S      
         names(i_gamma3) = 'gamma3'  ! gamma3 - 1) = '' dlnT_dlnRho at constant S            
         names(i_eta) = 'eta'  
      end subroutine get_result_names
      
      
      subroutine init_eos_handle_data(handle)
         use crlibm_lib
         integer, intent(in) :: handle
         type (EoS_General_Info), pointer :: rq
         rq => eos_handles(handle)
         rq% in_use = .true.
         rq% handle = handle
         
         ! set defaults
         rq% include_radiation = .true.
         rq% always_skip_elec_pos = .false.
         
         rq% logT_all_HELM = 7.7d0 ! HELM for lgT >= this
         rq% logT_all_OPAL = 7.6d0 ! OPAL for lgT <= this
         
         rq% logT_low_all_HELM = 2.2d0 ! HELM for lgT <= this
         rq% logT_low_all_SCVH = 2.3d0 ! SCVH for lgT >= this

         rq% mass_fraction_limit_for_PC = 1d-2
         
         rq% logRho1_OPAL_SCVH_limit = 3.6d0  ! must be <= 3.7
         rq% logRho2_OPAL_SCVH_limit = 3.0d0
         rq% logRho1_PC_limit = 2.999d0 
            ! keep < logRho2_OPAL_SCVH_limit so can go from table to PC without HELM
         rq% logRho2_PC_limit = 2.8d0 ! must be > 2.8 or so to avoid NaN's from PC

         !rq% Z_all_HELM = 0.040000001d0 ! old
         rq% Z_all_HELM = 1.000000001d0 ! new
         
         rq% log_Gamma_all_HELM = log10_cr(40d0)
         rq% Gamma_all_HELM = exp10_cr(rq% log_Gamma_all_HELM)
         rq% log_Gamma_all_PC = log10_cr(80d0)
         rq% PC_min_Z = 0.999d0
         
         rq% logE_all_HELM = 16.0d0 ! eosDE, HELM for lgE >= this
         rq% logE_all_OPAL = 15.9d0 ! eosDE, OPAL for lgE <= this
         
      end subroutine init_eos_handle_data
            
      
      subroutine do_free_eos_handle(handle)
         integer, intent(in) :: handle
         type (EoS_General_Info), pointer :: rq
         if (handle >= 1 .and. handle <= max_eos_handles) then
            rq => eos_handles(handle)
            eos_handles(handle)% in_use = .false.
         end if
      end subroutine do_free_eos_handle
      

      subroutine get_eos_ptr(handle,rq,ierr)
         integer, intent(in) :: handle
         type (EoS_General_Info), pointer :: rq
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_eos_handles) then
            ierr = -1
            return
         end if
         rq => eos_handles(handle)
         ierr = 0
      end subroutine get_eos_ptr
      
      
      subroutine do_set_HELM_flags( 
     >      handle, include_radiation, always_skip_elec_pos, ierr)
         integer, intent(in) :: handle
         logical, intent(in) :: include_radiation, always_skip_elec_pos
         integer, intent(out) :: ierr ! 0 means AOK.
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% include_radiation = include_radiation
         rq% always_skip_elec_pos = always_skip_elec_pos
      end subroutine do_set_HELM_flags

      
      subroutine do_set_Z_all_HELM(handle, Z_all_HELM, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: Z_all_HELM ! HELM for Z >= this
         integer, intent(out) :: ierr ! 0 means AOK.
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% Z_all_HELM = Z_all_HELM
      end subroutine do_set_Z_all_HELM

      
      subroutine do_set_logRhos_OPAL_SCVH( 
     >      handle, logRho1_OPAL_SCVH_limit, logRho2_OPAL_SCVH_limit, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: logRho1_OPAL_SCVH_limit 
            ! don't use OPAL_SCVH for logRho > this
         real(dp), intent(in) :: logRho2_OPAL_SCVH_limit 
            ! full OPAL_SCVH okay for logRho < this
         integer, intent(out) :: ierr ! 0 means AOK.
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% logRho1_OPAL_SCVH_limit = logRho1_OPAL_SCVH_limit
         rq% logRho2_OPAL_SCVH_limit = logRho2_OPAL_SCVH_limit
      end subroutine do_set_logRhos_OPAL_SCVH

      
      ! transition temperature zone for OPAL to HELM at high T
      subroutine do_set_HELM_OPAL_lgTs( 
     >      handle, logT_all_HELM, logT_all_OPAL, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: logT_all_HELM ! HELM for lgT >= this
         real(dp), intent(in) :: logT_all_OPAL ! OPAL/SCVH for lgT <= this
         integer, intent(out) :: ierr ! 0 means AOK.
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% logT_all_HELM = logT_all_HELM
         rq% logT_all_OPAL = logT_all_OPAL
      end subroutine do_set_HELM_OPAL_lgTs

      
      ! transition temperature zone for SCVH to HELM at low T
      subroutine do_set_HELM_SCVH_lgTs( 
     >      handle, logT_low_all_HELM, logT_low_all_SCVH, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: logT_low_all_HELM ! HELM for lgT <= this
         real(dp), intent(in) :: logT_low_all_SCVH ! SCVH for lgT >= this
         integer, intent(out) :: ierr ! 0 means AOK.
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% logT_low_all_HELM = logT_low_all_HELM
         rq% logT_low_all_SCVH = logT_low_all_SCVH
      end subroutine do_set_HELM_SCVH_lgTs


      subroutine do_set_PC_parameters(
     >      handle, mass_fraction_limit_for_PC,
     >      logRho1_PC_limit, logRho2_PC_limit,
     >      log_Gamma_all_HELM, log_Gamma_all_PC, PC_min_Z, ierr)
         use crlibm_lib
         integer, intent(in) :: handle
         real(dp), intent(in) :: mass_fraction_limit_for_PC ! skip species if abundance < this
         real(dp), intent(in) :: logRho1_PC_limit ! okay for pure PC for logRho > this
         real(dp), intent(in) :: logRho2_PC_limit ! don't use PC for logRho < this (>= 2.8)
         real(dp), intent(in) :: log_Gamma_all_HELM ! HELM for log_Gamma <= this
         real(dp), intent(in) :: log_Gamma_all_PC ! PC for log_Gamma >= this
         real(dp), intent(in) :: PC_min_Z ! don't use PC for Z < this
         integer, intent(out) :: ierr ! 0 means AOK.
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% mass_fraction_limit_for_PC = mass_fraction_limit_for_PC
         rq% logRho1_PC_limit = logRho1_PC_limit
         rq% logRho2_PC_limit = logRho2_PC_limit
         rq% log_Gamma_all_HELM = log_Gamma_all_HELM
         rq% Gamma_all_HELM = exp10_cr(log_Gamma_all_HELM)
         rq% log_Gamma_all_PC = log_Gamma_all_PC
         rq% PC_min_Z = PC_min_Z
      end subroutine do_set_PC_parameters
    
#ifdef offload
      
      subroutine init_target_eos(ierr) ! this is executed on the host
         integer, intent(out) :: ierr
         
         real(dp), pointer :: f_theta(:), x_theta(:)
         integer :: nx_theta
         
         ! rest is for Helm_Table
         logical :: with_coulomb_corrections
      
         integer :: imax
         integer :: jmax

         real(dp) :: logtlo ! log10 temp
         real(dp) :: logthi ! log10 temp
         real(dp) :: templo ! 10**logtlo
         real(dp) :: temphi ! 10**logthi
         real(dp) :: logtstp
         real(dp) :: logtstpi
         real(dp) :: logdlo ! log10 rho
         real(dp) :: logdhi ! log10 rho
         real(dp) :: denlo ! 10**logdlo
         real(dp) :: denhi ! 10**logdhi
         real(dp) :: logdstp
         real(dp) :: logdstpi
         real(dp), dimension(:), pointer :: d ! (imax) 
         real(dp), dimension(:), pointer :: t ! (jmax) 

         real(dp), dimension(:,:), pointer :: f ! (imax,jmax) 
         real(dp), dimension(:,:), pointer :: fd ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: ft ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: fdd ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: ftt ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: fdt ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: fddt ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: fdtt ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: fddtt ! (imax,jmax)

         real(dp), dimension(:,:), pointer :: dpdf ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: dpdfd ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: dpdft ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: dpdfdt ! (imax,jmax)

         real(dp), dimension(:,:), pointer :: ef ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: efd ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: eft ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: efdt ! (imax,jmax)

         real(dp), dimension(:,:), pointer :: xf ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: xfd ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: xft ! (imax,jmax)
         real(dp), dimension(:,:), pointer :: xfdt ! (imax,jmax)

         real(dp), dimension(:), pointer :: dt_sav ! (jmax)
         real(dp), dimension(:), pointer :: dt2_sav ! (jmax)
         real(dp), dimension(:), pointer :: dti_sav ! (jmax)
         real(dp), dimension(:), pointer :: dt2i_sav ! (jmax)
         real(dp), dimension(:), pointer :: dt3i_sav ! (jmax)
         real(dp), dimension(:), pointer :: dd_sav ! (imax)
         real(dp), dimension(:), pointer :: dd2_sav ! (imax)
         real(dp), dimension(:), pointer :: ddi_sav ! (imax)
         real(dp), dimension(:), pointer :: dd2i_sav ! (imax)
         real(dp), dimension(:), pointer :: dd3i_sav ! (jmax)

         ierr = 0
         
         !dir$ offload target(mic)
         call eos_def_init
         
         nx_theta = theta_e_nx
         if (.not. associated(f_theta_e1)) then
            write(*,*) 'init_target_eos: (.not. associated(f_theta_e1))'
            stop 1
         end if
         if (.not. associated(x_theta_e)) then
            write(*,*) 'init_target_eos: (.not. associated(x_theta_e))'
            stop 1
         end if
         f_theta => f_theta_e1
         x_theta => x_theta_e
         !dir$ offload target(mic)
         call copy_theta_e_info(nx_theta, f_theta, x_theta)
         
         ! now copy the Helm_Table
         with_coulomb_corrections = eos_ht% with_coulomb_corrections
         imax = eos_ht% imax
         jmax = eos_ht% jmax
         logtlo = eos_ht% logtlo ! log10 temp
         logthi = eos_ht% logthi ! log10 temp
         templo = eos_ht% templo ! 10**logtlo
         temphi = eos_ht% temphi ! 10**logthi
         logtstp = eos_ht% logtstp
         logtstpi = eos_ht% logtstpi
         logdlo = eos_ht% logdlo ! log10 rho
         logdhi = eos_ht% logdhi ! log10 rho
         denlo = eos_ht% denlo ! 10**logdlo
         denhi = eos_ht% denhi ! 10**logdhi
         logdstp = eos_ht% logdstp
         logdstpi = eos_ht% logdstpi
         d => eos_ht% d ! (imax) 
         t => eos_ht% t ! (jmax) 
         f => eos_ht% f ! (imax,jmax) 
         fd => eos_ht% fd ! (imax,jmax)
         ft => eos_ht% ft ! (imax,jmax)
         fdd => eos_ht% fdd ! (imax,jmax)
         ftt => eos_ht% ftt ! (imax,jmax)
         fdt => eos_ht% fdt ! (imax,jmax)
         fddt => eos_ht% fddt ! (imax,jmax)
         fdtt => eos_ht% fdtt ! (imax,jmax)
         fddtt => eos_ht% fddtt ! (imax,jmax)
         dpdf => eos_ht% dpdf ! (imax,jmax)
         dpdfd => eos_ht% dpdfd ! (imax,jmax)
         dpdft => eos_ht% dpdft ! (imax,jmax)
         dpdfdt => eos_ht% dpdfdt ! (imax,jmax)
         ef => eos_ht% ef ! (imax,jmax)
         efd => eos_ht% efd ! (imax,jmax)
         eft => eos_ht% eft ! (imax,jmax)
         efdt => eos_ht% efdt ! (imax,jmax)
         xf => eos_ht% xf ! (imax,jmax)
         xfd => eos_ht% xfd ! (imax,jmax)
         xft => eos_ht% xft ! (imax,jmax)
         xfdt => eos_ht% xfdt ! (imax,jmax)
         dt_sav => eos_ht% dt_sav ! (jmax)
         dt2_sav => eos_ht% dt2_sav ! (jmax)
         dti_sav => eos_ht% dti_sav ! (jmax)
         dt2i_sav => eos_ht% dt2i_sav ! (jmax)
         dt3i_sav => eos_ht% dt3i_sav ! (jmax)
         dd_sav => eos_ht% dd_sav ! (imax)
         dd2_sav => eos_ht% dd2_sav ! (imax)
         ddi_sav => eos_ht% ddi_sav ! (imax)
         dd2i_sav => eos_ht% dd2i_sav ! (imax)
         dd3i_sav => eos_ht% dd3i_sav ! (jmax)
 
         !dir$ offload target(mic) in(
     >      with_coulomb_corrections, ! 
     >      imax, ! 
     >      jmax, ! 
     >      logtlo, ! log10 temp
     >      logthi, ! log10 temp
     >      templo, ! 10**logtlo
     >      temphi, ! 10**logthi
     >      logtstp, ! 
     >      logtstpi, ! 
     >      logdlo, ! log10 rho
     >      logdhi, ! log10 rho
     >      denlo, ! 10**logdlo
     >      denhi, ! 10**logdhi
     >      logdstp, ! 
     >      logdstpi, ! 
     >      d, ! (imax) 
     >      t, ! (jmax) 
     >      f, ! (imax,jmax) 
     >      fd, ! (imax,jmax)
     >      ft, ! (imax,jmax)
     >      fdd, ! (imax,jmax)
     >      ftt, ! (imax,jmax)
     >      fdt, ! (imax,jmax)
     >      fddt, ! (imax,jmax)
     >      fdtt, ! (imax,jmax)
     >      fddtt, ! (imax,jmax)
     >      dpdf, ! (imax,jmax)
     >      dpdfd, ! (imax,jmax)
     >      dpdft, ! (imax,jmax)
     >      dpdfdt, ! (imax,jmax)
     >      ef, ! (imax,jmax)
     >      efd, ! (imax,jmax)
     >      eft, ! (imax,jmax)
     >      efdt, ! (imax,jmax)
     >      xf, ! (imax,jmax)
     >      xfd, ! (imax,jmax)
     >      xft, ! (imax,jmax)
     >      xfdt, ! (imax,jmax)
     >      dt_sav, ! (jmax)
     >      dt2_sav, ! (jmax)
     >      dti_sav, ! (jmax)
     >      dt2i_sav, ! (jmax)
     >      dt3i_sav, ! (jmax)
     >      dd_sav, ! (imax)
     >      dd2_sav, ! (imax)
     >      ddi_sav, ! (imax)
     >      dd2i_sav, ! (imax)
     >      dd3i_sav) ! (jmax)
         call do_copy_Helm_Table(
     >      with_coulomb_corrections, ! 
     >      imax, ! 
     >      jmax, ! 
     >      logtlo, ! log10 temp
     >      logthi, ! log10 temp
     >      templo, ! 10**logtlo
     >      temphi, ! 10**logthi
     >      logtstp, ! 
     >      logtstpi, ! 
     >      logdlo, ! log10 rho
     >      logdhi, ! log10 rho
     >      denlo, ! 10**logdlo
     >      denhi, ! 10**logdhi
     >      logdstp, ! 
     >      logdstpi, ! 
     >      d, ! (imax) 
     >      t, ! (jmax) 
     >      f, ! (imax,jmax) 
     >      fd, ! (imax,jmax)
     >      ft, ! (imax,jmax)
     >      fdd, ! (imax,jmax)
     >      ftt, ! (imax,jmax)
     >      fdt, ! (imax,jmax)
     >      fddt, ! (imax,jmax)
     >      fdtt, ! (imax,jmax)
     >      fddtt, ! (imax,jmax)
     >      dpdf, ! (imax,jmax)
     >      dpdfd, ! (imax,jmax)
     >      dpdft, ! (imax,jmax)
     >      dpdfdt, ! (imax,jmax)
     >      ef, ! (imax,jmax)
     >      efd, ! (imax,jmax)
     >      eft, ! (imax,jmax)
     >      efdt, ! (imax,jmax)
     >      xf, ! (imax,jmax)
     >      xfd, ! (imax,jmax)
     >      xft, ! (imax,jmax)
     >      xfdt, ! (imax,jmax)
     >      dt_sav, ! (jmax)
     >      dt2_sav, ! (jmax)
     >      dti_sav, ! (jmax)
     >      dt2i_sav, ! (jmax)
     >      dt3i_sav, ! (jmax)
     >      dd_sav, ! (imax)
     >      dd2_sav, ! (imax)
     >      ddi_sav, ! (imax)
     >      dd2i_sav, ! (imax)
     >      dd3i_sav) ! (jmax)
            
         !dir$ offload target(mic)
         call done_init_target_eos
         
      end subroutine init_target_eos
      
      
      subroutine done_init_target_eos ! this is executed on the target
         eosDT_is_initialized = .true.
         eosPT_is_initialized = .true.
         eosDE_is_initialized = .true.
      end subroutine done_init_target_eos


      subroutine copy_theta_e_info(nx, f, x) ! this is executed on the target
         integer, intent(in) :: nx
         real(dp), pointer, intent(in) :: f(:), x(:)
         integer :: i
         theta_e_nx = nx
         allocate(f_theta_e1(4*nx), x_theta_e(nx))
         f_theta_e(1:4,1:nx) => f_theta_e1(1:4*nx)
         do i=1,nx
            x_theta_e(i) = x(i)
         end do
         do i=1,4*nx
            f_theta_e1(i) = f(i)
         end do
      end subroutine copy_theta_e_info
      
      
      subroutine do_copy_Helm_Table( ! this is executed on the target
     >      with_coulomb_corrections, ! 
     >      imax, ! 
     >      jmax, ! 
     >      logtlo, ! log10 temp
     >      logthi, ! log10 temp
     >      templo, ! 10**logtlo
     >      temphi, ! 10**logthi
     >      logtstp, ! 
     >      logtstpi, ! 
     >      logdlo, ! log10 rho
     >      logdhi, ! log10 rho
     >      denlo, ! 10**logdlo
     >      denhi, ! 10**logdhi
     >      logdstp, ! 
     >      logdstpi, ! 
     >      d, ! (imax) 
     >      t, ! (jmax) 
     >      f, ! (imax,jmax) 
     >      fd, ! (imax,jmax)
     >      ft, ! (imax,jmax)
     >      fdd, ! (imax,jmax)
     >      ftt, ! (imax,jmax)
     >      fdt, ! (imax,jmax)
     >      fddt, ! (imax,jmax)
     >      fdtt, ! (imax,jmax)
     >      fddtt, ! (imax,jmax)
     >      dpdf, ! (imax,jmax)
     >      dpdfd, ! (imax,jmax)
     >      dpdft, ! (imax,jmax)
     >      dpdfdt, ! (imax,jmax)
     >      ef, ! (imax,jmax)
     >      efd, ! (imax,jmax)
     >      eft, ! (imax,jmax)
     >      efdt, ! (imax,jmax)
     >      xf, ! (imax,jmax)
     >      xfd, ! (imax,jmax)
     >      xft, ! (imax,jmax)
     >      xfdt, ! (imax,jmax)
     >      dt_sav, ! (jmax)
     >      dt2_sav, ! (jmax)
     >      dti_sav, ! (jmax)
     >      dt2i_sav, ! (jmax)
     >      dt3i_sav, ! (jmax)
     >      dd_sav, ! (imax)
     >      dd2_sav, ! (imax)
     >      ddi_sav, ! (imax)
     >      dd2i_sav, ! (imax)
     >      dd3i_sav) ! (jmax)
            
         ! controls
         logical :: with_coulomb_corrections
      
         ! sizes of the arrays
         integer :: imax
         integer :: jmax

         !..density and temperature
         real(dp), intent(in) :: logtlo ! log10 temp
         real(dp), intent(in) :: logthi ! log10 temp
         real(dp), intent(in) :: templo ! 10**logtlo
         real(dp), intent(in) :: temphi ! 10**logthi
         real(dp), intent(in) :: logtstp
         real(dp), intent(in) :: logtstpi
         real(dp), intent(in) :: logdlo ! log10 rho
         real(dp), intent(in) :: logdhi ! log10 rho
         real(dp), intent(in) :: denlo ! 10**logdlo
         real(dp), intent(in) :: denhi ! 10**logdhi
         real(dp), intent(in) :: logdstp
         real(dp), intent(in) :: logdstpi
         real(dp), intent(in), dimension(:), pointer :: d ! (imax) 
         real(dp), intent(in), dimension(:), pointer :: t ! (jmax) 

         !..for the helmholtz free energy tables
         real(dp), intent(in), dimension(:,:), pointer :: f ! (imax,jmax) 
         real(dp), intent(in), dimension(:,:), pointer :: fd ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: ft ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: fdd ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: ftt ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: fdt ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: fddt ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: fdtt ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: fddtt ! (imax,jmax)

         !..for the pressure derivative with density tables
         real(dp), intent(in), dimension(:,:), pointer :: dpdf ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: dpdfd ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: dpdft ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: dpdfdt ! (imax,jmax)

         !..for chemical potential tables
         real(dp), intent(in), dimension(:,:), pointer :: ef ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: efd ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: eft ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: efdt ! (imax,jmax)

         !..for the number density tables
         real(dp), intent(in), dimension(:,:), pointer :: xf ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: xfd ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: xft ! (imax,jmax)
         real(dp), intent(in), dimension(:,:), pointer :: xfdt ! (imax,jmax)

         !..for storing the differences
         real(dp), intent(in), dimension(:), pointer :: dt_sav ! (jmax)
         real(dp), intent(in), dimension(:), pointer :: dt2_sav ! (jmax)
         real(dp), intent(in), dimension(:), pointer :: dti_sav ! (jmax)
         real(dp), intent(in), dimension(:), pointer :: dt2i_sav ! (jmax)
         real(dp), intent(in), dimension(:), pointer :: dt3i_sav ! (jmax)
         real(dp), intent(in), dimension(:), pointer :: dd_sav ! (imax)
         real(dp), intent(in), dimension(:), pointer :: dd2_sav ! (imax)
         real(dp), intent(in), dimension(:), pointer :: ddi_sav ! (imax)
         real(dp), intent(in), dimension(:), pointer :: dd2i_sav ! (imax)
         real(dp), intent(in), dimension(:), pointer :: dd3i_sav ! (jmax)
         
         integer :: i, j
         
         allocate(eos_ht)
         
         allocate(
     >      eos_ht% d(imax),
     >      eos_ht% t(jmax),
     >      eos_ht% f(imax,jmax),
     >      eos_ht% fd(imax,jmax),
     >      eos_ht% ft(imax,jmax),
     >      eos_ht% fdd(imax,jmax),
     >      eos_ht% ftt(imax,jmax),
     >      eos_ht% fdt(imax,jmax),
     >      eos_ht% fddt(imax,jmax),
     >      eos_ht% fdtt(imax,jmax),
     >      eos_ht% fddtt(imax,jmax),
     >      eos_ht% dpdf(imax,jmax),
     >      eos_ht% dpdfd(imax,jmax),
     >      eos_ht% dpdft(imax,jmax),
     >      eos_ht% dpdfdt(imax,jmax),
     >      eos_ht% ef(imax,jmax),
     >      eos_ht% efd(imax,jmax),
     >      eos_ht% eft(imax,jmax),
     >      eos_ht% efdt(imax,jmax),
     >      eos_ht% xf(imax,jmax),
     >      eos_ht% xfd(imax,jmax),
     >      eos_ht% xft(imax,jmax),
     >      eos_ht% xfdt(imax,jmax),
     >      eos_ht% dt_sav(jmax),
     >      eos_ht% dt2_sav(jmax),
     >      eos_ht% dti_sav(jmax),
     >      eos_ht% dt2i_sav(jmax),
     >      eos_ht% dt3i_sav(jmax),
     >      eos_ht% dd_sav(imax),
     >      eos_ht% dd2_sav(imax),
     >      eos_ht% ddi_sav(imax),
     >      eos_ht% dd2i_sav(imax),
     >      eos_ht% dd3i_sav(jmax))

         eos_ht% with_coulomb_corrections = with_coulomb_corrections
         eos_ht% imax = imax
         eos_ht% jmax = jmax
         eos_ht% logtlo = logtlo ! log10 temp
         eos_ht% logthi = logthi ! log10 temp
         eos_ht% templo = templo ! 10**logtlo
         eos_ht% temphi = temphi ! 10**logthi
         eos_ht% logtstp = logtstp
         eos_ht% logtstpi = logtstpi
         eos_ht% logdlo = logdlo ! log10 rho
         eos_ht% logdhi = logdhi ! log10 rho
         eos_ht% denlo = denlo ! 10**logdlo
         eos_ht% denhi = denhi ! 10**logdhi
         eos_ht% logdstp = logdstp
         eos_ht% logdstpi = logdstpi
         
         do i=1,imax
            eos_ht% d(i) = d(i)
            eos_ht% dd_sav(i) = dd_sav(i)
            eos_ht% dd2_sav(i) = dd2_sav(i)
            eos_ht% ddi_sav(i) = ddi_sav(i)
            eos_ht% dd2i_sav(i) = dd2i_sav(i)
         end do
         do j=1,jmax
            eos_ht% t(j) = t(j)
            eos_ht% dt_sav(j) = dt_sav(j)
            eos_ht% dt2_sav(j) = dt2_sav(j)
            eos_ht% dti_sav(j) = dti_sav(j)
            eos_ht% dt2i_sav(j) = dt2i_sav(j)
            eos_ht% dt3i_sav(j) = dt3i_sav(j)
            eos_ht% dd3i_sav(j) = dd3i_sav(j)
         end do
         
         do j=1,jmax
            do i=1,imax
               eos_ht% f(i,j) = f(i,j) 
               eos_ht% fd(i,j) = fd(i,j)
               eos_ht% ft(i,j) = ft(i,j)
               eos_ht% fdd(i,j) = fdd(i,j)
               eos_ht% ftt(i,j) = ftt(i,j)
               eos_ht% fdt(i,j) = fdt(i,j)
               eos_ht% fddt(i,j) = fddt(i,j)
               eos_ht% fdtt(i,j) = fdtt(i,j)
               eos_ht% fddtt(i,j) = fddtt(i,j)
               eos_ht% dpdf(i,j) = dpdf(i,j)
               eos_ht% dpdfd(i,j) = dpdfd(i,j)
               eos_ht% dpdft(i,j) = dpdft(i,j)
               eos_ht% dpdfdt(i,j) = dpdfdt(i,j)
               eos_ht% ef(i,j) = ef(i,j)
               eos_ht% efd(i,j) = efd(i,j)
               eos_ht% eft(i,j) = eft(i,j)
               eos_ht% efdt(i,j) = efdt(i,j)
               eos_ht% xf(i,j) = xf(i,j)
               eos_ht% xfd(i,j) = xfd(i,j)
               eos_ht% xft(i,j) = xft(i,j)
               eos_ht% xfdt(i,j) = xfdt(i,j)
            end do
         end do
         
      end subroutine do_copy_Helm_Table
      
      
      subroutine copy_eosDT_table(ix, iz, ierr)
         integer, intent(in) :: ix, iz
         integer, intent(out) :: ierr

         real :: logQ_min ! logQ = logRho - 2*logT + 12
         real :: logQ_max
         real :: del_logQ ! spacing for the logQs
         integer :: num_logQs
         real :: logT_min
         real :: logT_max
         real :: del_logT ! spacing for the logTs
         integer :: num_logTs
         real, pointer :: logQs(:), logTs(:)
         real, pointer :: tbl1(:)
         integer :: version
         type (EosDT_XZ_Info), pointer :: ep
         
         ep => eosDT_XZ_data(ix,iz)     
         logQ_min = ep% logQ_min
         logQ_max = ep% logQ_max
         del_logQ = ep% del_logQ
         num_logQs = ep% num_logQs
         logT_min = ep% logT_min
         logT_max = ep% logT_max
         del_logT = ep% del_logT
         num_logTs = ep% num_logTs
         logQs => ep% logQs
         logTs => ep% logTs
         tbl1 => ep% tbl1
         version = ep% version

         !dir$ offload target(mic) out(ierr) in( 
     >      ix, iz, 
     >      logQ_min, logQ_max, del_logQ, num_logQs, 
     >      logT_min, logT_max, del_logT, num_logTs, 
     >      logQs, logTs, tbl1, version)
         call do_copy_eosDT_table( 
     >      ix, iz, 
     >      logQ_min, logQ_max, del_logQ, num_logQs, 
     >      logT_min, logT_max, del_logT, num_logTs, 
     >      logQs, logTs, tbl1, version, ierr)
         
      end subroutine copy_eosDT_table
      
      
      subroutine do_copy_eosDT_table( 
     >      ix, iz, 
     >      logQ_min, logQ_max, del_logQ, num_logQs, 
     >      logT_min, logT_max, del_logT, num_logTs, 
     >      logQs, logTs, tbl1, version, ierr)
         integer, intent(in) :: ix, iz
         real, intent(in) :: logQ_min ! logQ = logRho - 2*logT + 12
         real, intent(in) :: logQ_max
         real, intent(in) :: del_logQ ! spacing for the logQs
         integer, intent(in) :: num_logQs
         real, intent(in) :: logT_min
         real, intent(in) :: logT_max
         real, intent(in) :: del_logT ! spacing for the logTs
         integer, intent(in) :: num_logTs
         real, pointer, intent(in) :: logQs(:), logTs(:)
         real, pointer, intent(in) :: tbl1(:) ! =(sz_per_eos_point, num_eos_vals, num_logQs, num_logTs)
         integer, intent(in) :: version
         integer, intent(out) :: ierr
         
         type (EosDT_XZ_Info), pointer :: ep
         integer :: i
         
         ierr = 0
         ep => eosDT_XZ_data(ix,iz) 
         allocate( 
     >      ep% tbl1(sz_per_eos_point*num_eos_vals*num_logQs*num_logTs), 
     >      ep% logQs(num_logQs), ep% logTs(num_logTs), stat=ierr)
         if (ierr /= 0) return
             
         ep% logQ_min = logQ_min
         ep% logQ_max = logQ_max
         ep% del_logQ = del_logQ
         ep% num_logQs = num_logQs
         ep% logT_min = logT_min
         ep% logT_max = logT_max
         ep% del_logT = del_logT
         ep% num_logTs = num_logTs
         ep% version = version
         do i=1,num_logQs
            ep% logQs(i) = logQs(i)
         end do
         do i=1,num_logTs
            ep% logTs(i) = logTs(i)
         end do
         do i=1,size(tbl1,dim=1)
            ep% tbl1(i) = tbl1(i)
         end do
         
      end subroutine do_copy_eosDT_table
      
      
      subroutine copy_eosPT_table(ix, iz, ierr)
         integer, intent(in) :: ix, iz
         integer, intent(out) :: ierr

         real :: logW_min ! logW = logRho - 2*logT + 12
         real :: logW_max
         real :: del_logW ! spacing for the logWs
         integer :: num_logWs
         real :: logT_min
         real :: logT_max
         real :: del_logT ! spacing for the logTs
         integer :: num_logTs
         real, pointer :: logWs(:), logTs(:)
         real, pointer :: tbl1(:) ! =(sz_per_eos_point, num_eos_vals, num_logWs, num_logTs)
         integer :: version
         type (EosPT_XZ_Info), pointer :: ep
         
         ierr = 0
         
         ep => eosPT_XZ_data(ix,iz)     
         logW_min = ep% logW_min
         logW_max = ep% logW_max
         del_logW = ep% del_logW
         num_logWs = ep% num_logWs
         logT_min = ep% logT_min
         logT_max = ep% logT_max
         del_logT = ep% del_logT
         num_logTs = ep% num_logTs
         logWs => ep% logWs
         logTs => ep% logTs
         tbl1 => ep% tbl1
         version = ep% version

         !dir$ offload target(mic) out(ierr) in(
     >      ix, iz, 
     >      logW_min, logW_max, del_logW, num_logWs, 
     >      logT_min, logT_max, del_logT, num_logTs, 
     >      logWs, logTs, tbl1, version)
         call do_copy_eosPT_table( 
     >      ix, iz, 
     >      logW_min, logW_max, del_logW, num_logWs, 
     >      logT_min, logT_max, del_logT, num_logTs, 
     >      logWs, logTs, tbl1, version, ierr)
         
      end subroutine copy_eosPT_table
      
      
      subroutine do_copy_eosPT_table( 
     >      ix, iz, 
     >      logW_min, logW_max, del_logW, num_logWs, 
     >      logT_min, logT_max, del_logT, num_logTs, 
     >      logWs, logTs, tbl1, version, ierr)
         integer, intent(in) :: ix, iz
         real, intent(in) :: logW_min
         real, intent(in) :: logW_max
         real, intent(in) :: del_logW ! spacing for the logWs
         integer, intent(in) :: num_logWs
         real, intent(in) :: logT_min
         real, intent(in) :: logT_max
         real, intent(in) :: del_logT ! spacing for the logTs
         integer, intent(in) :: num_logTs
         real, pointer, intent(in) :: logWs(:), logTs(:)
         real, pointer, intent(in) :: tbl1(:)
         integer, intent(in) :: version
         integer, intent(out) :: ierr
         
         type (EosPT_XZ_Info), pointer :: ep
         integer :: i
         
         ierr = 0
         ep => eosPT_XZ_data(ix,iz) 
         allocate( 
     >      ep% tbl1(sz_per_eos_point*num_eos_vals*num_logWs*num_logTs), 
     >      ep% logWs(num_logWs), ep% logTs(num_logTs), stat=ierr)
         if (ierr /= 0) return
             
         ep% logW_min = logW_min
         ep% logW_max = logW_max
         ep% del_logW = del_logW
         ep% num_logWs = num_logWs
         ep% logT_min = logT_min
         ep% logT_max = logT_max
         ep% del_logT = del_logT
         ep% num_logTs = num_logTs
         ep% version = version
         do i=1,num_logWs
            ep% logWs(i) = logWs(i)
         end do
         do i=1,num_logTs
            ep% logTs(i) = logTs(i)
         end do
         do i=1,size(tbl1,dim=1)
            ep% tbl1(i) = tbl1(i)
         end do
         
      end subroutine do_copy_eosPT_table
      
      
#endif
      
   
#ifdef offload
      !dir$ end options
#endif

      end module eos_def
      
