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


      module kap_lib
      ! library for calculating opacities

      ! the data interface for the library is defined in kap_def
      
      use const_def, only: dp
      
      implicit none


      contains ! the procedure interface for the library
      ! client programs should only call these routines.
            
      
      ! call this routine to initialize the kap module. 
      ! only needs to be done once at start of run.
      ! Reads data from the 'kap' directory in the data_dir.
      ! If use_cache is true and there is a 'kap/cache' directory, it will try that first.
      ! If it doesn't find what it needs in the cache, 
      ! it reads the data and writes the cache for next time.    
      subroutine kap_init( &
            kap_file_prefix, CO_prefix, lowT_prefix, &
            blend_logT_upper_bdy, blend_logT_lower_bdy, &
            type2_logT_lower_bdy, use_cache, kap_cache_dir, info)      
         use kap_def, only : kap_def_init, kap_is_initialized, kap_init_on_coprocessor
         use load_kap, only : Setup_Kap_Tables
         character(*), intent(in) :: &
            kap_file_prefix, CO_prefix, lowT_prefix
         real(dp), intent(in) :: blend_logT_upper_bdy, blend_logT_lower_bdy
         real(dp), intent(in) :: type2_logT_lower_bdy
         logical, intent(in) :: use_cache
         character (len=*), intent(in) :: kap_cache_dir ! '' means use default
         integer, intent(out) :: info ! 0 means AOK.         
#ifdef offload
         logical, parameter :: load_on_demand = .false.
#else
         logical, parameter :: load_on_demand = .true.
#endif
         info = 0    
         if (kap_is_initialized) return
         call kap_def_init(lowT_prefix, kap_cache_dir)
         call Setup_Kap_Tables( &
            kap_file_prefix, CO_prefix, lowT_prefix, &
            blend_logT_upper_bdy, blend_logT_lower_bdy, &
            type2_logT_lower_bdy, use_cache, load_on_demand, info)
#ifdef offload
         if (info /= 0) return
         call kap_init_on_coprocessor(info)
#endif
         if (info /= 0) return
         kap_is_initialized = .true.
      end subroutine kap_init

      
      subroutine kap_shutdown
         use kap_def
         use load_kap, only : Free_Kap_Tables
         call Free_Kap_Tables  
         kap_is_initialized = .false.
      end subroutine kap_shutdown

      
      ! after kap_init has finished, you can allocate a "handle".
      
      integer function alloc_kap_handle(ierr)
         use kap_def,only: do_alloc_kap,kap_is_initialized
         integer, intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         alloc_kap_handle = do_alloc_kap(ierr)
      end function alloc_kap_handle      
      
      subroutine free_kap_handle(handle)
         ! frees the handle and all associated data
         use kap_def,only: Kap_General_Info,do_free_kap,get_kap_ptr
         integer, intent(in) :: handle
         call do_free_kap(handle)
      end subroutine free_kap_handle

      
      
      
      
      subroutine kap_get_op_mono( &
            handle, zbar, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            ! args for op_mono
            use_op_mono_alt_get_kap, &
            nel, izzp, fap, screening, umesh, ff, rs, ss, &
            ! output
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use kap_eval, only: combine_rad_with_compton_and_conduction
         use crlibm_lib
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: log10_rho ! the density
         real(dp), intent(in) :: log10_T ! the temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         ! args for op_mono_get_kap
         logical, intent(in) :: use_op_mono_alt_get_kap
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(:) ! (nel)
         real(dp), intent(in) :: fap(:) ! (nel)
         logical, intent(in) :: screening
         ! work arrays
         real, pointer :: umesh(:), ff(:,:,:,:), rs(:,:,:), ss(:,:,:,:)
            ! umesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! rs(nptot, 4, 4)
            ! ss(nptot, nrad, 4, 4)
         ! output
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         real(dp) :: g1, gp1(nel), kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         real :: logRho, logT, Rho, T
         type (Kap_General_Info), pointer :: rq

         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return

         logRho = real(log10_rho); Rho = real(exp10_cr(log10_rho))
         logT = real(log10_T); T = real(exp10_cr(log10_T))
         
         if (use_op_mono_alt_get_kap) then
            call op_mono_alt_get_kap( &
               nel, izzp, fap, log10_T, log10_rho, screening, &
               g1, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, &
               umesh, ff, rs, ierr)
         else
            call op_mono_get_kap( &
               nel, izzp, fap, log10_T, log10_rho, screening, &
               g1, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, gp1, &
               umesh, ff, rs, ss, ierr)
         end if
         if (ierr /= 0) return
         
         kap_rad = exp10_cr(g1)
         
         call combine_rad_with_compton_and_conduction( &
            rq, rho, logRho, T, logT, real(zbar), &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine kap_get_op_mono
      
      
      
      ! interface to OP routines as modified by Haili Hu for radiative levitation in diffusion
      
      ! ref: Hu et al MNRAS 418 (2011)
      
      subroutine load_op_mono_data(op_mono_data_path, op_mono_data_cache_filename, ierr)
         use kap_def
         use op_load, only: op_dload
         character (len=*), intent(in) :: op_mono_data_path, op_mono_data_cache_filename
         integer, intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call op_dload(op_mono_data_path, op_mono_data_cache_filename, ierr)
      end subroutine load_op_mono_data
      
      
      ! sizes for work arrays
      subroutine get_op_mono_params(op_nptot, op_ipe, op_nrad)
         use kap_def
         use op_def, only: nptot, ipe, nrad
         integer, intent(out) :: op_nptot, op_ipe, op_nrad
         op_nptot = nptot
         op_ipe = ipe
         op_nrad = nrad
      end subroutine get_op_mono_params
      
      
! HH: Based on "op_ax.f"
! Input:   kk = number of elements to calculate g_rad for   
!          iz1(kk) = charge of element to calculate g_rad for   
!          nel = number of elements in mixture   
!          izzp(nel) = charge of elements   
!          fap(nel) = number fractions of elements   
!          flux = local radiative flux (Lrad/4*pi*r^2)   
!          fltp = log10 T   
!          flrhop = log10 rho   
!          screening   if true, use screening corrections   
! Output: g1 = log10 kappa   
!         gx1 = d(log kappa)/d(log T)   
!         gy1 = d(log kappa)/d(log rho)   
!         gp1(kk) = d(log kappa)/d(log xi)    
!         grl1(kk) = log10 grad   
!         fx1(kk) = d(log grad)/d(log T)    
!         fy1(kk) = d(log grad)/d(log rho)   
!         grlp1(kk) = d(log grad)/d(log chi),
!              chi is the fraction with which the number fraction is varied, i.e.:
!                 chi = nf_new/nf_previous
!                 where nf is the number fraction
!         meanZ(nel) = average ionic charge of elements   
!         zetx1(nel) = d(meanZ)/d(log T)    
!         zety1(nel) = d(meanZ)/d(log rho)   
!         ierr = 0 for correct use
!         ierr = 101 for rho out of range for this T 
!         ierr = 102 for T out of range
      subroutine op_mono_get_radacc( &
            kk, izk, nel, izzp, fap, flux, fltp, flrhop, screening, &
            g1, gx1, gy1, gp1, grl1, fx1, fy1, grlp1, meanZ, zetx1, zety1, &
            umesh, ff, ta, rs, s, ierr)
         use kap_def
         use op_eval, only: eval_op_radacc
         integer, intent(in) :: kk, nel
         integer, intent(in) :: izk(kk), izzp(nel)
         real(dp), intent(in) :: fap(nel)
         real(dp), intent(in) :: flux, fltp
         real(dp), intent(in) :: flrhop
         logical, intent(in) :: screening
         real(dp), intent(out) :: g1, gx1, gy1
         real(dp), intent(out) :: &
            grl1(kk), meanZ(nel), grlp1(kk), gp1(kk), &
            fx1(kk), fy1(kk), zetx1(nel), zety1(nel)
         ! work arrays 
         real, pointer :: umesh(:), ff(:,:,:,:), ta(:,:,:,:), rs(:,:,:), s(:,:,:,:)
            ! umesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! ta(nptot, nrad, 4, 4), 
            ! rs(nptot, 4, 4)
            ! s(nptot, nrad, 4, 4)
         integer,intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call eval_op_radacc( &
            kk, izk, nel, izzp, fap, flux, fltp, flrhop, screening, &
            g1, gx1, gy1, gp1, grl1, fx1, fy1, grlp1, meanZ, zetx1, zety1, &
            umesh, ff, ta, rs, s, ierr)
      end subroutine op_mono_get_radacc
      

! note: for op mono, elements must come from the set given in op_mono_element_Z in kap_def.
      
! HH: Based on "op_mx.f", opacity calculations to be used for stellar evolution calculations 
! Input:   nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fltp = log10 (temperature)
!          flrhop = log10 (mass density) 
!          screening   if true, use screening corrections
! Output: g1 = log10 kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         ierr = 0 for correct use
!         ierr = 101 for rho out of range for this T 
!         ierr = 102 for T out of range
      subroutine op_mono_get_kap( &
            nel, izzp, fap, fltp, flrhop, screening, &
            g1, gx1, gy1, gp1, &
            umesh, ff, rs, s, ierr)
         use kap_def
         use op_eval, only: eval_op_ev
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(nel)
         real(dp), intent(in) :: fap(nel)
         real(dp), intent(in) :: fltp, flrhop
         logical, intent(in) :: screening
         real(dp), intent(out) :: g1, gx1, gy1, gp1(nel)
         ! work arrays
         real, pointer :: umesh(:), ff(:,:,:,:), rs(:,:,:), s(:,:,:,:)
            ! umesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! rs(nptot, 4, 4)
            ! s(nptot, nrad, 4, 4)
         integer,intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call eval_op_ev( &
            nel, izzp, fap, fltp, flrhop, screening, &
            g1, gx1, gy1, gp1, &
            umesh, ff, rs, s, ierr)
      end subroutine op_mono_get_kap
      
      
! HH: Based on "op_mx.f", opacity calculations to be used for non-adiabatic pulsation calculations
! Special care is taken to ensure smoothness of opacity derivatives
! Input:   nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fltp = log10 (temperature)
!          flrhop = log10 (mass density) 
!          screening   if true, use screening corrections
! Output: g1 = log10 kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         ierr = 0 for correct use
!         ierr = 101 for rho out of range for this T 
!         ierr = 102 for T out of range
      subroutine op_mono_alt_get_kap( &
            nel, izzp, fap, fltp, flrhop, screening, &
            g1, gx1, gy1, &
            umesh, ff, rs, ierr)
         use kap_def
         use op_eval, only: eval_alt_op
         implicit none
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(nel)
         real(dp), intent(in) :: fap(nel)
         real(dp), intent(in) :: fltp, flrhop
         logical, intent(in) :: screening
         real(dp), intent(out) :: g1, gx1, gy1
         ! work arrays
         real, pointer :: umesh(:), ff(:,:,:,:), rs(:,:,:)
            ! umesh(nptot)
            ! ff(nptot, ipe, 0:5, 0:5)
            ! rs(nptot, 0:5, 0:5)
         integer,intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call eval_alt_op( &
            nel, izzp, fap, fltp, flrhop, screening, &
            g1, gx1, gy1, &
            umesh, ff, rs, ierr)
      end subroutine op_mono_alt_get_kap
      
      
		subroutine get_op_mono_args( &
		      species, X, min_X_to_include, chem_id, &
		      nel, izzp, fap, ierr)
		   use chem_def, only: chem_isos
         use kap_def
         integer, intent(in) :: species, chem_id(:)
		   real(dp), intent(in) :: X(:) ! mass fractions (assumed baryonic)
		   real(dp), intent(in) :: min_X_to_include ! skip iso if X < this
         integer, intent(out) :: nel
         integer, intent(out) :: izzp(:)
         real(dp), intent(out) :: fap(:)
         integer,intent(out) :: ierr
         
         integer :: i, cid, j, Z, iel
         real(dp) :: tot
         
         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         
         nel = 0
         izzp(:) = 0
         fap(:) = 0d0
         
         do i=1,species
            if (X(i) < min_X_to_include) cycle
            cid = chem_id(i)
            Z = chem_isos% Z(cid)
            if (Z == 0) cycle
            ! change Z if necessary so that in op set
            do j=num_op_mono_elements,1,-1
               if (Z >= op_mono_element_Z(j)) then
                  Z = op_mono_element_Z(j)
                  exit
               end if
            end do
            iel = 0
            do j=1,nel
               if (izzp(j) == Z) then
                  iel = j
                  exit
               end if
            end do
            if (iel == 0) then
               nel = nel+1
               iel = nel
               izzp(nel) = Z
            end if
            fap(iel) = fap(iel) + X(i)/dble(chem_isos% Z_plus_N(cid))
         end do
         
         tot = sum(fap(1:nel))
         if (tot <= 0d0) then
            ierr = -1
            return
         end if
         
         do j=1,nel
            fap(j) = fap(j)/tot ! number fractions
         end do
		   
		end subroutine get_op_mono_args
		
      
      
      
      
#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif            
      
      
      subroutine kap_ptr(handle,rq,ierr)
         use kap_def,only:Kap_General_Info,get_kap_ptr,kap_is_initialized
         integer, intent(in) :: handle ! from alloc_kap_handle
         type (Kap_General_Info), pointer :: rq
         integer, intent(out):: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call get_kap_ptr(handle,rq,ierr)
      end subroutine kap_ptr
      

      ! kap evaluation
      ! you can call these routines after you've setup the tables for the handle.
      ! NOTE: the structures referenced via the handle are read-only
      ! for the evaulation routines, so you can do multiple evaluations in parallel
      ! using the same handle. 
      
      


      subroutine kap_set_choices( &
            handle, cubic_interpolation_in_X, cubic_interpolation_in_Z, &
            include_electron_conduction, &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
            ierr)
         use kap_def, only: do_kap_set_choices,kap_is_initialized
         integer, intent(in) :: handle ! from alloc_kap_handle
         logical, intent(in) :: &
            cubic_interpolation_in_X, cubic_interpolation_in_Z, include_electron_conduction
         real(dp), intent(in) :: &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ
         integer, intent(out) :: ierr ! 0 means AOK.
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call do_kap_set_choices( &
            handle, cubic_interpolation_in_X, cubic_interpolation_in_Z, &
            include_electron_conduction, &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
            ierr)
#ifdef offload
         !dir$ offload target(mic) out(ierr) &
            in(handle, cubic_interpolation_in_X, cubic_interpolation_in_Z, &
            include_electron_conduction, &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ)
         call do_kap_set_choices( &
            handle, cubic_interpolation_in_X, cubic_interpolation_in_Z, &
            include_electron_conduction, &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
            ierr)
#endif
      end subroutine kap_set_choices
           

      
      subroutine kap_get_Type1( &
            handle, zbar, X, Z, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def, only : kap_is_initialized
         ! INPUT
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: X ! the hydrogen mass fraction
         real(dp), intent(in) :: Z ! the metallicity
         real(dp), intent(in) :: log10_rho ! the density
         real(dp), intent(in) :: log10_T ! the temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         
         ! OUTPUT
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         real(dp), parameter :: max_frac_Type2 = 0d0, &
            XC = 0d0, XN = 0d0, XO = 0d0, XNe = 0d0
         real(dp) :: frac_Type2

         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call kap_get_blend_1_2( &
            max_frac_Type2, handle, zbar, X, Z, Z, XC, XN, XO, XNe, &
            log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine kap_get_Type1
           

      ! Type2 opacities use a base metallicity and enhancements to C and O
      ! from the 3 definitions
         ! Z = Zbase + dXC + dXO ! total metals
         ! XC = dXC + base_fC*Zbase ! total mass fraction of carbon
         ! XO = dXO + base_fO*Zbase ! total mass fraction of oxygen
      ! we get expressions for the 3 parameters, Zbase, dXC, and dXO
      ! using the base values for fractions of C and O and
      ! the base_fC and base_fO provided in the call to kap_set_choices.
      subroutine kap_get_Type2( &
            handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def, only : kap_is_initialized
         ! INPUT
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: X, Z, Zbase, XC, XN, XO, XNe ! composition    
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         
         ! OUTPUT
         real(dp), intent(out) :: frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         real(dp), parameter :: max_frac_Type2 = 1d0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call kap_get_blend_1_2( &
            max_frac_Type2, handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
            log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
                     
      end subroutine kap_get_Type2


      subroutine kap_get_blend_1_2( &
            max_frac_Type2, handle, zbar, &
            X, Z, Zbase_in, XC, XN, XO, XNe, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use kap_eval, only: Get_kap_Results
         use crlibm_lib
         use utils_lib, only: is_bad_num

         ! INPUT
         real(dp), intent(in) :: max_frac_Type2
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: X, Z, Zbase_in, XC, XN, XO, XNe ! composition    
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         
         ! OUTPUT
         real(dp), intent(out) :: frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         type (Kap_General_Info), pointer :: rq
         real :: logRho, logT, Rho, T
         real(dp) :: frac_Type1, Zbase, dZ, frac_X, frac_dZ, &
            dXC, dXO, fC, fN, fO, fNe, fHeavy, ZHeavy, dXsum
         logical, parameter :: dbg = .false.
         include 'formats'
         
         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         
         logRho = real(log10_rho); Rho = real(exp10_cr(log10_rho))
         logT = real(log10_T); T = real(exp10_cr(log10_T))
         
         Zbase = real(Zbase_in)
         dXC = 0
         dxO = 0
         frac_Type2 = max_frac_Type2
         
         if (frac_Type2 > 0d0) then
            
            fC = kap_co_z_tables(1)% Zfrac_C
            fN = kap_co_z_tables(1)% Zfrac_N
            fO = kap_co_z_tables(1)% Zfrac_O
            fNe = kap_co_z_tables(1)% Zfrac_Ne
         
            if (Zbase < 0d0) then
               write(*,1) 'must supply Zbase for Type2 opacities', Zbase_in
               write(*,1) 'frac_Type2', frac_Type2
               ierr = -1
               return
               fHeavy = 1d0 - fC - fN - fO - fNe
               ZHeavy = Z - XC - XN - XO - XNe
               Zbase = ZHeavy/fHeavy
            end if
            
            dXC = max(0d0, xC - fC*Zbase)
            dXO = max(0d0, xO - fO*Zbase)
            
            if (X >= rq% kap_Type2_full_off_X) then
               frac_X = 0
            else if (X <= rq% kap_Type2_full_on_X) then
               frac_X = 1
            else ! blend
               frac_X = (rq% kap_Type2_full_off_X - X) / &
                  (rq% kap_Type2_full_off_X - rq% kap_Type2_full_on_X)
            end if
            
            dZ = Z - Zbase
            if (dZ <= rq% kap_Type2_full_off_dZ) then
               frac_dZ = 0
            else if (dZ >= rq% kap_Type2_full_on_dZ) then
               frac_dZ = 1
            else ! blend
               frac_dZ = (rq% kap_Type2_full_off_dZ - dZ) / &
                  (rq% kap_Type2_full_off_dZ - rq% kap_Type2_full_on_dZ)
            end if
            
            frac_Type2 = frac_Type2*frac_X*frac_dZ
            if (frac_Type2 == 0d0) then
               dXC = 0
               dXO = 0
            end if
            
            if (is_bad_num(dXC) .or. is_bad_num(dXO)) then
               write(*,1) 'X', X
               write(*,1) 'Z', Z
               write(*,1) 'Zbase_in', Zbase_in
               write(*,1) 'XC', XC
               write(*,1) 'XN', XN
               write(*,1) 'XO', XO
               write(*,1) 'XNe', XNe
               write(*,1) 'log10_rho', log10_rho
               write(*,1) 'log10_T', log10_T
               write(*,1) 'zbar', zbar
               write(*,1) 'lnfree_e', lnfree_e
               write(*,1) 'd_lnfree_e_dlnRho', d_lnfree_e_dlnRho
               write(*,1) 'd_lnfree_e_dlnT', d_lnfree_e_dlnT
               write(*,1) 'dXC', dXC
               write(*,1) 'dXO', dXO
               write(*,1) 'frac_X', frac_X
               write(*,1) 'frac_dZ', frac_dZ
               write(*,1) 'frac_Type2', frac_Type2
               ierr = -1
               return
               stop 'kap_get_blend_1_2'
            end if
         
         end if

         frac_Type1 = 1d0 - frac_Type2
         if (dbg) then
            write(*,1) 'max_frac_Type2', max_frac_Type2
            write(*,1) 'frac_Type2', frac_Type2
            write(*,1) 'frac_Type1', frac_Type1
            write(*,*)
            write(*,1) 'X', X
            write(*,1) 'dXC', dXC
            write(*,1) 'dXO', dXO
            write(*,1) 'Zbase', Zbase
         end if
         
         call Get_kap_Results( &
            rq, real(zbar), real(X), real(dXC), real(dXO), real(Zbase), real(Z), &
            frac_Type1, Rho, logRho, T, logT,  &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
                     
      end subroutine kap_get_blend_1_2
      
      
      subroutine kap_get_elect_cond_opacity( &
            zbar, log10_rho, log10_T,  &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use condint, only: do_electron_conduction
         use kap_def, only : kap_is_initialized
         ! INPUT
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: log10_rho ! the density
         real(dp), intent(in) :: log10_T ! the temperature
         
         ! OUTPUT
         real(dp), intent(out) :: kap ! electron conduction opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         real(dp) :: kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT

         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call do_electron_conduction( &
            real(zbar), real(log10_rho), real(log10_T), &
            kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
         kap = dble(kap_ec)
         dlnkap_dlnRho = dble(dlnkap_ec_dlnRho)
         dlnkap_dlnT = dble(dlnkap_ec_dlnT)
         
      end subroutine kap_get_elect_cond_opacity
      
      
      subroutine kap_get_compton_opacity( &
            Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_eval, only: Compton_Opacity
         use kap_def, only : kap_is_initialized
         ! INPUT
         real(dp), intent(in) :: Rho, T
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         
         ! OUTPUT
         real(dp), intent(out) :: kap ! electron conduction opacity
         real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.

         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call Compton_Opacity(Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine kap_get_compton_opacity
      
#ifdef offload
      !dir$ end options
#endif
		
      end module kap_lib

