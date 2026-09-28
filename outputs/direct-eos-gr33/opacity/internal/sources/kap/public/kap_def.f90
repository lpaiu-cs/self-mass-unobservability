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

      module kap_def
      use const_def, only: dp
      
      implicit none
      
#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif
      
      ! info about op_mono elements
      integer, parameter :: num_op_mono_elements = 17
      integer :: op_mono_element_Z(num_op_mono_elements) 
      character(len=2) :: op_mono_element_name(num_op_mono_elements) 
      real(dp) :: op_mono_element_mass(num_op_mono_elements)
      
      
      integer, parameter :: kap_table_fixed_metal_form = 1
      integer, parameter :: kap_table_co_enhanced_form = 2
      
      
      ! for fixed metal tables (no enhancements)
      type Kap_Z_Table ! holds pointers to all the X tables for a particular Z
         logical :: lowT_flag
         real :: Z ! the Z for this table
         integer :: num_Xs ! number of X's for this Z
         type (Kap_X_Table), dimension(:), pointer :: x_tables ! stored in order of increasing X
      end type Kap_Z_Table

      type (Kap_Z_Table), dimension(:), pointer :: kap_z_tables ! stored in order of increasing Z
      type (Kap_Z_Table), dimension(:), pointer :: kap_lowT_z_tables ! stored in order of increasing Z
      
      
      ! note:  logR = logRho - 3*logT + 18
      
      type Kap_X_Table
         logical :: not_loaded_yet
         real :: X
         real :: Z
         real :: logR_min
         real :: logR_max
         integer :: num_logRs
         integer :: ili_logRs ! =1 if logRs are evenly spaced
         real, pointer :: logRs(:) ! indexed from 1 to num_logRs
         real :: logT_min
         real :: logT_max
         integer :: num_logTs
         integer :: ili_logTs ! =1 if logTs are evenly spaced
         real, pointer :: logTs(:) ! indexed from 1 to num_logTs
         real, pointer :: kap1(:)
      end type Kap_X_Table
      
      integer, parameter :: num_kap_Xs = 10
      ! 0,  .1, .2, .35, .5, .7, .8, .9, .95, 1-Z
      real(dp), parameter :: kap_Xs(1:num_kap_Xs) =  &
            (/ 0.00d0, 0.10d0, 0.20d0, 0.35d0, 0.50d0, 0.70d0, 0.80d0, 0.90d0, 0.95d0, -1d0 /)
         ! -1 means X = 1-Z for whatever Z
      
      integer, parameter :: num_kap_Zs = 13
      real(dp), parameter :: kap_Zs(1:num_kap_Zs) =  &
            (/ 0.000d0, 0.0001d0, 0.0003d0, 0.001d0, 0.002d0, 0.004d0, &
            0.01d0, 0.02d0, 0.03d0, 0.04d0, 0.06d0, 0.08d0, 0.100d0 /)
         
      integer, parameter :: num_kap_Xs_for_this_Z(1:num_kap_Zs) =  &
            (/ num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, &
            num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, &
            num_kap_Xs-2, num_kap_Xs-2, num_kap_Xs-2 /)      
      
      integer, parameter :: num_kap_Freedman_Zs = 7
      real(dp), parameter :: kap_Freedman_Zs(1:num_kap_Freedman_Zs) =  &
            (/ 0.01d0, 0.02d0, 0.04d0, 0.100d0, 0.200d0, 0.63d0, 1.00d0 /)
      
      ! for C/O enhanced tables
      type Kap_CO_Z_Table
         real :: Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne
         real :: log10_Zbase ! log10(Zbase)
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables ! stored in order of increasing X
            ! the X tables need not be equally spaced
      end type Kap_CO_Z_Table

      type (Kap_CO_Z_Table), dimension(:), pointer :: kap_co_z_tables ! stored in order of increasing Z

      integer, parameter :: num_kap_CO_Zs = 8
      real(dp), parameter :: kap_CO_Zs(1:num_kap_CO_Zs) =  &
         (/ 0.000, 0.001, 0.004, 0.010, 0.020, 0.030, 0.050, 0.100 /)
      
      integer, parameter :: num_kap_CO_Xs = 5
      real(dp), parameter :: kap_CO_Xs(1:num_kap_CO_Xs) =  &
         (/ 0.00d0, 0.03d0, 0.10d0, 0.35d0, 0.70d0 /)
      
      integer, parameter :: num_kap_CO_dXs = 8
      real(dp), parameter :: kap_CO_dXs(num_kap_CO_dXs) =  &
         (/ 0.00d0, 0.01d0, 0.03d0, 0.10d0, 0.20d0, 0.40d0, 0.60d0, 1.0d0 /)
      
      type Kap_CO_Table
         integer :: table_num ! the table number from the data file
         real :: X
         real :: Z
         real :: dXC
         real :: dXO
         real :: dXC_lookup
         real :: dXO_lookup
         real, dimension(:), pointer :: kap1
      end type Kap_CO_Table
      
      
      integer, parameter :: max_num_CO_tables = 70
      
      ! standard number of CO tables for each X+Z combo
      !           X  
      ! Z         0.0   0.1   0.35  0.7
      ! 0.0       58    58    51    32
      ! 0.001     58    58    51    30
      ! 0.004     58    58    51    30
      ! 0.010     58    58    51    30
      ! 0.020     58    58    51    30
      ! 0.030     58    58    49    30
      ! 0.100     58    53    43    26
      
      ! 1362 tables total
      
      
      
      type Kap_CO_X_Table
      
         logical :: not_loaded_yet
         real :: X
         real :: Z
         real :: logR_min
         real :: logR_max
         integer :: num_logRs
         integer :: ili_logRs ! =1 if logRs are evenly spaced
         real, dimension(:), pointer :: logRs ! indexed from 1 to num_logRs
         real :: logT_min
         real :: logT_max
         integer :: num_logTs
         integer :: ili_logTs ! =1 if logTs are evenly spaced
         real, dimension(:), pointer :: logTs ! indexed from 1 to num_logTs
         
         integer :: num_CO_tables
         ! the tables are in 3 groups
         ! 1) tables with dXC > dXO, ordered by increasing dXO, and by increasing dXC within same dXO.
         ! 2) tables with dXC = dXO, ordered by increasing value.
         ! 3) tables with dXC < dXO, ordered by increasing dXC, and by increasing dXO within same dXC.
         ! the spacing of dXC's is the same as dXO's, so there are as many tables in 3) as in 1).
         integer :: num_dXC_gt_dXO ! the number of tables with dXC > dXO
         integer :: CO_table_numbers(num_kap_CO_dXs,num_kap_CO_dXs) 
            ! entry (i,j) is the co_index for table with dXC=Xs(i) and dXO=Xs(j), or -1 if no such table.
         integer :: next_dXO_table(max_num_CO_tables) 
            ! entry (i) is the co_index for the table with same dXC and next larger dXO, or -1 if none such.
         integer :: next_dXC_table(max_num_CO_tables) 
            ! entry (i) is the co_index for the table with same dXO and next larger dXC, or -1 if none such.
         type (Kap_CO_Table), dimension(:), pointer :: co_tables
          
      end type Kap_CO_X_Table


      type (Kap_Z_Table), pointer :: z_tables_copy(:)
      type (Kap_X_Table), pointer :: x_tables_copy(:)      
      type (Kap_CO_X_Table), pointer :: co_x_tables_copy(:)
      type (Kap_CO_Table), pointer :: co_tables_copy(:)
      

      type Kap_General_Info
      
         logical :: cubic_interpolation_in_X
         logical :: cubic_interpolation_in_Z
         logical :: include_electron_conduction
         
         ! switch to Type1 if X too large
         real(dp) :: kap_Type2_full_off_X ! Type2 full off for X >= this
         real(dp) :: kap_Type2_full_on_X ! Type2 full on for X <= this
   
         ! switch to Type1 if dZ too small (dZ = Z - Zbase)
         real(dp) :: kap_Type2_full_off_dZ ! Type2 is full off for dZ <= this
         real(dp) :: kap_Type2_full_on_dZ ! Type2 can be full on for dZ >= this

         ! for logR > 1, we extrapolate the radiative opacities.
         ! for low T, this can run into problems, so we need to clip logT when logR > 1.
         real(dp) :: min_logT_for_logR_gt_1 ! 3.5 is the default for this
         ! this clipping doesn't apply to special low T tables
         ! i.e. skip it if min logT in table is < 2.5
         
         ! bookkeeping
         integer :: handle
         logical :: in_use
         
      end type Kap_General_Info



      ! NOTE: in the following, "log" means base 10, "ln" means natural log, and units are cgs.

      integer, parameter :: sz_per_kap_point = 4
      !
      ! function f(x,y) with samples f(i,j) has bicubic spline fit s(x,y).
      ! compact representation of spline fit uses 4 entries as follows:
      !
      ! d(1,i,j) = s(i,j)
      ! d(2,i,j) = d2s_dx2(i,j)
      ! d(3,i,j) = d2s_dy2(i,j)
      ! d(4,i,j) = d4s_dx2_dy2(i,j)
      ! 
      ! given f(i,j), the spline fitting code can compute the other entries
      !
      ! given d(1:4,i,j), spline interpolation code can compute s(x,y)
      ! and also the partials ds_dx(x,y) and ds_dy(x,y)
      !

      
      logical :: kap_is_initialized = .false.
      
      character (len=1000) :: kap_dir, kap_cache_dir
      character (len=256) :: kap_prefix = 'gn93'
      character (len=256) :: kap_CO_prefix = 'gn93_co'
      character (len=256) :: kap_lowT_prefix = 'lowT_fa05_gn93'
      logical :: kap_use_cache = .true.

      ! blending in T is done between the following limits
      real(dp) :: kap_blend_logT_upper_bdy = 3.88 ! old value was 4.1d0
      real(dp) :: kap_blend_logT_lower_bdy = 3.80 ! old value was 4.0d0
      ! last time I looked, the table bottom for the higher T tables was logT = 3.75
      ! while max logT for the lower T Freeman tables was 4.5
      ! so for those, you need to keep kap_blend_logT_upper_bdy < 4.5
      ! and kap_blend_logT_lower_bdy > 3.75 
      ! it is also probably a good idea to keep the blend away from H ionization
      ! logT upper of about 3.9 or a bit less will do that.

      real(dp) :: kap_type2_logT_lower_bdy = 3.8d0
      ! below this, turn off type2

      logical :: kap_using_lowT_Freedman = .false.
         ! this is set true if kap_lowT_prefix starts with 'lowT_Freedman'
      
      logical :: clip_to_kap_table_boundaries ! typically, this should be set true.
         ! if this is set true, then temperature and density args are
         ! clipped to the boundaries of the table.
      
      integer, parameter :: max_kap_handles = 10
      type (Kap_General_Info), target :: kap_handles(max_kap_handles)
      

      contains
      
      
      subroutine kap_def_init(lowT_prefix, kap_cache_dir_in)
         use chem_def
         use utils_lib, only : mkdir
         use const_def, only: mesa_data_dir, mesa_caches_dir
         
         character (*), intent(in) :: lowT_prefix, kap_cache_dir_in
         integer :: i
         logical :: lowT_Freedman
         
         if (len_trim(kap_cache_dir_in) > 0) then
            kap_cache_dir = kap_cache_dir_in
         else if (len_trim(mesa_caches_dir) > 0) then
            kap_cache_dir = trim(mesa_caches_dir) // '/kap_cache'
         else
            kap_cache_dir = trim(mesa_data_dir) // '/kap_data/cache'
         end if
         call mkdir(kap_cache_dir)
         
         if (len_trim(lowT_prefix) < 13) then
            lowT_Freedman = .false.
         else
            lowT_Freedman = lowT_prefix(1:13) == 'lowT_Freedman'
         end if
         
         call do_kap_def_init(lowT_Freedman)

      end subroutine kap_def_init
      
      
      subroutine do_kap_def_init(lowT_Freedman)
         logical, intent(in) :: lowT_Freedman
         integer :: i
         
         kap_using_lowT_Freedman = lowT_Freedman
         clip_to_kap_table_boundaries = .true.
         
         do i=1,max_kap_handles
            kap_handles(i)% handle = i
            kap_handles(i)% in_use = .false.
         end do
         
         op_mono_element_Z(1:num_op_mono_elements) = (/ &
            1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 24, 25, 26, 28 /)
         op_mono_element_name(1:num_op_mono_elements) = (/ &
            'h ', 'he', 'c ', 'n ', 'o ', 'ne', 'na', 'mg', 'al', &
            'si', 's ', 'ar', 'ca', 'cr', 'mn', 'fe', 'ni' /)
         op_mono_element_mass(1:num_op_mono_elements) = (/ &
            1.0080d0, 4.0026d0, 12.0111d0, 14.0067d0, 15.9994d0, 20.179d0, &
            22.9898d0, 24.305d0, 26.9815d0, 28.086d0, 32.06d0, 39.948d0, &
            40.08d0, 51.996d0, 54.9380d0, 55.847d0, 58.71d0 /)
      
      end subroutine do_kap_def_init


      
      integer function do_alloc_kap(ierr)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do_alloc_kap = -1
!$omp critical (kap_handle)
         do i = 1, max_kap_handles
            if (.not. kap_handles(i)% in_use) then
               kap_handles(i)% in_use = .true.
               do_alloc_kap = i
               exit
            end if
         end do
!$omp end critical (kap_handle)
         if (do_alloc_kap == -1) then
            ierr = -1
            return
         end if
         if (kap_handles(do_alloc_kap)% handle /= do_alloc_kap) then
            ierr = -1
            return
         end if
         call init_kap_handle_data(do_alloc_kap)
#ifdef offload
         !dir$ offload target(mic) in(do_alloc_kap)
         call init_kap_handle_data(do_alloc_kap)
#endif
      end function do_alloc_kap
      
      
      subroutine init_kap_handle_data(i) ! set defaults
         integer, intent(in) :: i
         kap_handles(i)% cubic_interpolation_in_X = .true.
         kap_handles(i)% cubic_interpolation_in_Z = .false. 
         kap_handles(i)% include_electron_conduction = .true.
         kap_handles(i)% min_logT_for_logR_gt_1 = 3.5
         kap_handles(i)% kap_Type2_full_off_X = 0.71d0 ! Type2 full off for X >= this
         kap_handles(i)% kap_Type2_full_on_X = 0.70d0 ! Type2 full on for X <= this
         kap_handles(i)% kap_Type2_full_off_dZ = 0.001d0 ! Type2 is full off for dZ <= this
         kap_handles(i)% kap_Type2_full_on_dZ = 0.01d0 ! Type2 can be full on for dZ >= this
      end subroutine init_kap_handle_data
      
      
      subroutine do_free_kap(handle)
         integer, intent(in) :: handle
         if (handle >= 1 .and. handle <= max_kap_handles) then
            kap_handles(handle)% in_use = .false.
         end if
      end subroutine do_free_kap
      

      subroutine get_kap_ptr(handle,rq,ierr)
         integer, intent(in) :: handle
         type (Kap_General_Info), pointer :: rq
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_kap_handles) then
            ierr = -1
            return
         end if
         rq => kap_handles(handle)
         ierr = 0
      end subroutine get_kap_ptr
      
      
      subroutine do_kap_set_choices( &
            handle, cubic_interpolation_in_X, cubic_interpolation_in_Z, &
            include_electron_conduction, &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
            ierr)
         integer, intent(in) :: handle ! from alloc_kap_handle
         logical, intent(in) :: &
            cubic_interpolation_in_X, cubic_interpolation_in_Z, include_electron_conduction
         real(dp), intent(in) :: &
            kap_Type2_full_off_X, kap_Type2_full_on_X, &
            kap_Type2_full_off_dZ, kap_Type2_full_on_dZ
         integer, intent(out) :: ierr ! 0 means AOK.
         type (Kap_General_Info), pointer :: rq
         ierr = 0
         call get_kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return
         rq% cubic_interpolation_in_X = cubic_interpolation_in_X
         rq% cubic_interpolation_in_Z = cubic_interpolation_in_Z
         rq% include_electron_conduction = include_electron_conduction
         if (kap_Type2_full_off_X > 0) &
            rq% kap_Type2_full_off_X = kap_Type2_full_off_X
         if (kap_Type2_full_on_X > 0) &
            rq% kap_Type2_full_on_X = kap_Type2_full_on_X
         if (kap_Type2_full_off_dZ > 0) &
            rq% kap_Type2_full_off_dZ = kap_Type2_full_off_dZ
         if (kap_Type2_full_on_dZ > 0) &
            rq% kap_Type2_full_on_dZ = kap_Type2_full_on_dZ
      end subroutine do_kap_set_choices
      
      
      subroutine kap_init_on_coprocessor(ierr) ! this routine runs on host
         integer, intent(out):: ierr         
         logical :: lowT_Freedman         
         ierr = 0         
#ifdef offload
         lowT_Freedman = kap_using_lowT_Freedman
         !dir$ offload target(mic) in(lowT_Freedman)
         call do_kap_def_init(lowT_Freedman)
         call copy_z_tables(.false.,ierr)
         if (ierr /= 0) return
         call copy_z_tables(.true.,ierr)
         if (ierr /= 0) return
         call copy_co_z_tables(ierr)
         if (ierr /= 0) return
         !dir$ offload target(mic)
         call done_kap_init_on_coprocessor
#endif
      end subroutine kap_init_on_coprocessor

#ifdef offload
      
      subroutine done_kap_init_on_coprocessor ! on mic
         kap_is_initialized = .true.
      end subroutine done_kap_init_on_coprocessor
      
      subroutine copy_z_tables(lowT_flag,ierr) ! on host
         logical, intent(in) :: lowT_flag
         integer, intent(out) :: ierr         
         integer :: num, i
         type (Kap_Z_Table), pointer :: zt(:)   
         ierr = 0
         if (lowT_flag) then
            zt => kap_lowT_z_tables
         else
            zt => kap_z_tables
         end if
         if (.not. associated(zt)) return
         num = size(zt, dim=1)
         !dir$ offload target(mic)
         call start_copy_z_tables(lowT_flag, num, ierr)
         if (ierr /= 0) return
         do i=1,num
            call copy_z_table(zt, i, ierr)
            if (ierr /= 0) return
         end do         
      end subroutine copy_z_tables
      
      subroutine start_copy_z_tables(lowT_flag,num,ierr) ! on mic
         logical, intent(in) :: lowT_flag
         integer, intent(in) :: num
         integer, intent(out) :: ierr
         ierr = 0
         allocate(z_tables_copy(num),stat=ierr)
         if (lowT_flag) then
            kap_lowT_z_tables => z_tables_copy
         else
            kap_z_tables => z_tables_copy
         end if
      end subroutine start_copy_z_tables
      
      subroutine copy_z_table(zt,i,ierr) ! on host
         type (Kap_Z_Table), intent(in), pointer :: zt(:)
         integer, intent(in) :: i
         integer, intent(out) :: ierr         
         integer :: j, num, num_Xs
         real :: Z
         logical :: lowT_flag
         ierr = 0
         num = zt(i)% num_Xs
         !dir$ offload target(mic) out(ierr) in(num)
         call start_copy_x_tables(num,ierr)
         if (ierr /= 0) return         
         do j=1,num
            call copy_x_table(zt(i)% x_tables, j, ierr)
            if (ierr /= 0) return
         end do         
         lowT_flag = zt(i)% lowT_flag
         Z = zt(i)% Z
         num_Xs = zt(i)% num_Xs
         !dir$ offload target(mic) out(ierr) in(i, lowT_flag, Z, num_Xs)
         call do_copy_z_table(i, lowT_flag, Z, num_Xs, ierr)            
      end subroutine copy_z_table
      
      subroutine do_copy_z_table( & ! on mic
            i, lowT_flag, Z, num_Xs, ierr)
         integer, intent(in) :: i, num_Xs
         logical, intent(in) :: lowT_flag
         real, intent(in) :: Z
         integer, intent(out) :: ierr
         ierr = 0
         z_tables_copy(i)% lowT_flag = lowT_flag
         z_tables_copy(i)% Z = Z
         z_tables_copy(i)% num_Xs = num_Xs
         z_tables_copy(i)% x_tables => x_tables_copy
         nullify(x_tables_copy)
      end subroutine do_copy_z_table     
      
      subroutine start_copy_x_tables(num,ierr) ! on mic
         integer, intent(in) :: num
         integer, intent(out) :: ierr
         allocate(x_tables_copy(num),stat=ierr)
      end subroutine start_copy_x_tables
      
      subroutine copy_x_table(xt,i,ierr) ! on host
         type (Kap_X_Table), intent(in), pointer :: xt(:)
         integer, intent(in) :: i
         integer, intent(out) :: ierr
         real :: X, Z
         real :: logR_min, logR_max
         integer :: num_logRs, ili_logRs
         real, pointer :: logRs(:) ! indexed from 1 to num_logRs
         real :: logT_min, logT_max
         integer :: num_logTs, ili_logTs
         real, pointer :: logTs(:) ! indexed from 1 to num_logTs
         real, pointer :: kap1(:)
         if (xt(i)% not_loaded_yet) then
            ierr = -1
            return
         end if
         ierr = 0
         X = xt(i)% X
         Z = xt(i)% Z
         logR_min = xt(i)% logR_min
         logR_max = xt(i)% logR_max
         num_logRs = xt(i)% num_logRs
         ili_logRs = xt(i)% ili_logRs
         logT_min = xt(i)% logT_min
         logT_max = xt(i)% logT_max
         num_logTs = xt(i)% num_logTs
         ili_logTs = xt(i)% ili_logTs
         logRs => xt(i)% logRs
         logTs => xt(i)% logTs
         kap1 => xt(i)% kap1
         !dir$ offload target(mic) out(ierr) in( &
            i, X, Z, &
            logR_min, logR_max, num_logRs, ili_logRs, logRs, &
            logT_min, logT_max, num_logTs, ili_logTs, logTs, &
            kap1)
         call do_copy_x_table( &
            i, X, Z, &
            logR_min, logR_max, num_logRs, ili_logRs, logRs, &
            logT_min, logT_max, num_logTs, ili_logTs, logTs, &
            kap1, ierr)
      end subroutine copy_x_table

      subroutine do_copy_x_table( & ! on mic
            i, X, Z, &
            logR_min, logR_max, num_logRs, ili_logRs, logRs, &
            logT_min, logT_max, num_logTs, ili_logTs, logTs, &
            kap1, ierr)
         integer, intent(in) :: i
         real, intent(in) :: X, Z
         real, intent(in) :: logR_min, logR_max
         integer, intent(in) :: num_logRs, ili_logRs
         real, pointer, intent(in) :: logRs(:) ! indexed from 1 to num_logRs
         real, intent(in) :: logT_min, logT_max
         integer, intent(in) :: num_logTs, ili_logTs
         real, pointer, intent(in) :: logTs(:) ! indexed from 1 to num_logTs
         real, pointer, intent(in) :: kap1(:)
         integer, intent(out) :: ierr         
         integer :: j, sz
         ierr = 0
         x_tables_copy(i)% not_loaded_yet = .false.
         x_tables_copy(i)% X = X
         x_tables_copy(i)% Z = Z
         x_tables_copy(i)% logR_min = logR_min
         x_tables_copy(i)% logR_max = logR_max
         x_tables_copy(i)% num_logRs = num_logRs
         x_tables_copy(i)% ili_logRs = ili_logRs
         allocate(x_tables_copy(i)% logRs(num_logRs), stat=ierr)
         if (ierr /= 0) return
         do j = 1, num_logRs
            x_tables_copy(i)% logRs(j) = logRs(j)
         end do
         x_tables_copy(i)% logT_min = logT_min
         x_tables_copy(i)% logT_max = logT_max
         x_tables_copy(i)% num_logTs = num_logTs
         x_tables_copy(i)% ili_logTs = ili_logTs
         allocate(x_tables_copy(i)% logTs(num_logTs), stat=ierr)
         if (ierr /= 0) return
         do j = 1, num_logTs
            x_tables_copy(i)% logTs(j) = logTs(j)
         end do
         sz = size(kap1,dim=1)
         allocate(x_tables_copy(i)% kap1(sz), stat=ierr)
         if (ierr /= 0) return
         do j = 1, sz
            x_tables_copy(i)% kap1(j) = kap1(j)
         end do
      end subroutine do_copy_x_table

      subroutine copy_co_z_tables(ierr) ! on host
         integer, intent(out) :: ierr         
         integer :: num, i
         ierr = 0
         if (.not. associated(kap_co_z_tables)) return
         num = size(kap_co_z_tables, dim=1)
         !dir$ offload target(mic)
         call start_copy_co_z_tables(num,ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in start_copy_co_z_tables'
            stop 1
            return
         end if
         do i=1,num
            call copy_co_z_table(i,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in copy_co_z_table', i
               stop 1
               return
            end if
         end do         
      end subroutine copy_co_z_tables
      
      subroutine start_copy_co_z_tables(num,ierr) ! on mic
         integer, intent(in) :: num
         integer, intent(out) :: ierr
         ierr = 0
         allocate(kap_co_z_tables(num),stat=ierr)
      end subroutine start_copy_co_z_tables
      
      subroutine copy_co_z_table(i,ierr) ! on host
         integer, intent(in) :: i
         integer, intent(out) :: ierr         
         integer :: j, num
         real :: Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, log10_Zbase
         ierr = 0
         num = size(kap_co_z_tables(i)% x_tables,dim=1)
         !dir$ offload target(mic) out(ierr) in(num)
         call start_copy_co_x_tables(num, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in start_copy_co_x_tables', i
            stop 1
            return
         end if
         do j=1,num
            call copy_co_x_table(kap_co_z_tables(i)% x_tables, j, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in copy_co_x_table', i, j
               stop 1
               return
            end if
         end do      
         Zbase = kap_co_z_tables(i)% Zbase
         Zfrac_C = kap_co_z_tables(i)% Zfrac_C
         Zfrac_N = kap_co_z_tables(i)% Zfrac_N
         Zfrac_O = kap_co_z_tables(i)% Zfrac_O
         Zfrac_Ne = kap_co_z_tables(i)% Zfrac_Ne
         log10_Zbase = kap_co_z_tables(i)% log10_Zbase
         !dir$ offload target(mic) out(ierr) in( &
            i, Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, log10_Zbase)
         call do_copy_co_z_table( &
            i, Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, log10_Zbase, ierr)            
      end subroutine copy_co_z_table
      
      subroutine do_copy_co_z_table( & ! on mic
            i, Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, log10_Zbase, ierr)
         integer, intent(in) :: i
         real, intent(in) :: Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne, log10_Zbase
         integer, intent(out) :: ierr
         ierr = 0
         kap_co_z_tables(i)% Zbase = Zbase
         kap_co_z_tables(i)% Zfrac_C = Zfrac_C
         kap_co_z_tables(i)% Zfrac_N = Zfrac_N
         kap_co_z_tables(i)% Zfrac_O = Zfrac_O
         kap_co_z_tables(i)% Zfrac_Ne = Zfrac_Ne
         kap_co_z_tables(i)% log10_Zbase = log10_Zbase
         kap_co_z_tables(i)% x_tables => co_x_tables_copy
         nullify(co_x_tables_copy)
      end subroutine do_copy_co_z_table
      
      subroutine start_copy_co_x_tables(num,ierr) ! on mic
         integer, intent(in) :: num
         integer, intent(out) :: ierr
         ierr = 0
         allocate(co_x_tables_copy(num),stat=ierr)
      end subroutine start_copy_co_x_tables
      
      subroutine copy_co_x_table(xt,i,ierr) ! on host
         type (Kap_CO_X_Table), pointer :: xt(:)
         integer, intent(in) :: i
         integer, intent(out) :: ierr
         integer :: num, j
         real :: X, Z
         real :: logR_min, logR_max
         integer :: num_logRs, ili_logRs
         real, pointer :: logRs(:) ! indexed from 1 to num_logRs
         real :: logT_min, logT_max
         integer :: num_logTs, ili_logTs
         real, pointer :: logTs(:) ! indexed from 1 to num_logTs
         integer :: num_CO_tables, num_dXC_gt_dXO
         integer :: CO_table_numbers(num_kap_CO_dXs,num_kap_CO_dXs) 
         integer :: next_dXO_table(max_num_CO_tables) 
         integer :: next_dXC_table(max_num_CO_tables) 
         ierr = 0
         num = size(xt(i)% co_tables, dim=1)
         !dir$ offload target(mic) out(ierr) in(num)
         call start_copy_co_tables(num, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in start_copy_co_tables'
            stop 1
            return
         end if
         do j=1,num
            call copy_co_table(xt(i)% co_tables, j, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in copy_co_table', i, j
               stop 1
               return
            end if
         end do         
         X = xt(i)% X
         Z = xt(i)% Z
         logR_min = xt(i)% logR_min
         logR_max = xt(i)% logR_max
         num_logRs = xt(i)% num_logRs
         ili_logRs = xt(i)% ili_logRs
         logRs => xt(i)% logRs
         logT_min = xt(i)% logT_min
         logT_max = xt(i)% logT_max
         num_logTs = xt(i)% num_logTs
         ili_logTs = xt(i)% ili_logTs
         logTs => xt(i)% logTs
         num_CO_tables = xt(i)% num_CO_tables
         num_dXC_gt_dXO = xt(i)% num_dXC_gt_dXO
         CO_table_numbers = xt(i)% CO_table_numbers
         next_dXO_table = xt(i)% next_dXO_table
         next_dXC_table = xt(i)% next_dXC_table
         !dir$ offload target(mic) out(ierr) in( &
            i, X, Z, &
            logR_min, logR_max, num_logRs, ili_logRs, logRs, &
            logT_min, logT_max, num_logTs, ili_logTs, logTs, &
            num_CO_tables, num_dXC_gt_dXO, CO_table_numbers, &
            next_dXO_table, next_dXC_table)
         call do_copy_co_x_table( &
            i, X, Z, &
            logR_min, logR_max, num_logRs, ili_logRs, logRs, &
            logT_min, logT_max, num_logTs, ili_logTs, logTs, &
            num_CO_tables, num_dXC_gt_dXO, CO_table_numbers, &
            next_dXO_table, next_dXC_table, ierr)            
      end subroutine copy_co_x_table

      subroutine do_copy_co_x_table( & ! on mic
            i, X, Z, &
            logR_min, logR_max, num_logRs, ili_logRs, logRs, &
            logT_min, logT_max, num_logTs, ili_logTs, logTs, &
            num_CO_tables, num_dXC_gt_dXO, CO_table_numbers, &
            next_dXO_table, next_dXC_table, &
            ierr)
         integer, intent(in) :: i
         real, intent(in) :: X, Z
         real, intent(in) :: logR_min, logR_max
         integer, intent(in) :: num_logRs, ili_logRs
         real, pointer, intent(in) :: logRs(:) ! indexed from 1 to num_logRs
         real, intent(in) :: logT_min, logT_max
         integer, intent(in) :: num_logTs, ili_logTs
         real, pointer, intent(in) :: logTs(:) ! indexed from 1 to num_logTs
         integer, intent(in) :: num_CO_tables, num_dXC_gt_dXO
         integer, intent(in) :: CO_table_numbers(num_kap_CO_dXs,num_kap_CO_dXs) 
         integer, intent(in) :: next_dXO_table(max_num_CO_tables) 
         integer, intent(in) :: next_dXC_table(max_num_CO_tables) 
         integer, intent(out) :: ierr         
         integer :: j
         ierr = 0
         co_x_tables_copy(i)% not_loaded_yet = .false.
         co_x_tables_copy(i)% X = X
         co_x_tables_copy(i)% Z = Z
         co_x_tables_copy(i)% logR_min = logR_min
         co_x_tables_copy(i)% logR_max = logR_max
         co_x_tables_copy(i)% num_logRs = num_logRs
         co_x_tables_copy(i)% ili_logRs = ili_logRs
         allocate(co_x_tables_copy(i)% logRs(num_logRs))
         do j = 1, num_logRs
            co_x_tables_copy(i)% logRs(j) = logRs(j)
         end do
         co_x_tables_copy(i)% logT_min = logT_min
         co_x_tables_copy(i)% logT_max = logT_max
         co_x_tables_copy(i)% num_logTs = num_logTs
         co_x_tables_copy(i)% ili_logTs = ili_logTs
         allocate(co_x_tables_copy(i)% logTs(num_logTs))
         do j = 1, num_logTs
            co_x_tables_copy(i)% logTs(j) = logTs(j)
         end do         
         co_x_tables_copy(i)% num_CO_tables = num_CO_tables
         co_x_tables_copy(i)% num_dXC_gt_dXO = num_dXC_gt_dXO
         co_x_tables_copy(i)% CO_table_numbers = CO_table_numbers
         co_x_tables_copy(i)% next_dXO_table = next_dXO_table
         co_x_tables_copy(i)% next_dXC_table = next_dXC_table
         co_x_tables_copy(i)% co_tables => co_tables_copy
         nullify(co_tables_copy)
      end subroutine do_copy_co_x_table
      
      subroutine start_copy_co_tables(num,ierr) ! on mic
         integer, intent(in) :: num
         integer, intent(out) :: ierr
         allocate(co_tables_copy(num),stat=ierr)
      end subroutine start_copy_co_tables

      subroutine copy_co_table(co, i, ierr) ! on host
         type (Kap_CO_Table), pointer :: co(:)
         integer, intent(in) :: i
         integer, intent(out) :: ierr
         integer :: table_num
         real :: X, Z, dXC, dXO, dXC_lookup, dXO_lookup
         real, pointer :: kap1(:)
         ierr = 0
         table_num = co(i)% table_num
         X = co(i)% X
         Z = co(i)% Z
         dXC = co(i)% dXC
         dXO = co(i)% dXO
         dXC_lookup = co(i)% dXC_lookup
         dXO_lookup = co(i)% dXO_lookup
         kap1 => co(i)% kap1
         !dir$ offload target(mic) out(ierr) in( &
            i, table_num, X, Z, dXC, dXO, dXC_lookup, dXO_lookup, kap1)
         call do_copy_co_table( &
            i, table_num, X, Z, dXC, dXO, dXC_lookup, dXO_lookup, kap1, ierr)
      end subroutine copy_co_table

      subroutine do_copy_co_table( & ! on mic
            i, table_num, X, Z, dXC, dXO, dXC_lookup, dXO_lookup, kap1, ierr)
         integer, intent(in) :: i, table_num
         real, intent(in) :: X, Z, dXC, dXO, dXC_lookup, dXO_lookup
         real, pointer :: kap1(:)
         integer, intent(out) :: ierr
         integer :: sz, j
         ierr = 0
         co_tables_copy(i)% table_num = table_num
         co_tables_copy(i)% X = X
         co_tables_copy(i)% Z = Z
         co_tables_copy(i)% dXC = dXC
         co_tables_copy(i)% dXO = dXO
         co_tables_copy(i)% dXC_lookup = dXC_lookup
         co_tables_copy(i)% dXO_lookup = dXO_lookup
         sz = size(kap1,dim=1)
         allocate(co_tables_copy(i)% kap1(sz),stat=ierr)
         if (ierr /= 0) return
         do j=1,sz
            co_tables_copy(i)% kap1(j) = kap1(j)
         end do
      end subroutine do_copy_co_table
#endif

#ifdef offload
      !dir$ end options
#endif

      end module kap_def
      
