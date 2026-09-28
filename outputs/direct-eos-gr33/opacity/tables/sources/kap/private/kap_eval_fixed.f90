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

      module kap_eval_fixed
      use utils_lib,only:is_bad_num
      use kap_eval_support
      use const_def, only: dp
      use crlibm_lib
      
      implicit none
      
      contains
      
#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif            
      
      subroutine Get_kap_fixed_metal_Results( &
               rq, zbar, Z, X, Rho, logRho, T, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use const_def, only: pi
         
         ! INPUT
         type (Kap_General_Info), pointer :: rq
         real, intent(in) :: zbar, Z, X
         real, intent(inout) :: Rho, logRho ! can be modified to clip to table boundaries
         real, intent(inout) :: T, logT

         ! OUTPUT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
         
         integer :: lowT_num_Zs
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         real(dp) :: alfa, beta, &
            kap1, dlnkap1_dlnRho, dlnkap1_dlnT, &
            kap2, dlnkap2_dlnRho, dlnkap2_dlnT, &
            lower_bdy, upper_bdy
         
         if (kap_using_lowT_Freedman) then
            lowT_num_Zs = num_kap_Freedman_Zs
         else
            lowT_num_Zs = num_kap_Zs
         end if

         lower_bdy = max(kap_blend_logT_lower_bdy, kap_z_tables(1)% x_tables(1)% logT_min)
         if (logT <= lower_bdy) then ! all lowT
            if (dbg) write(*,*) 'all lowT call Get1_kap_fixed_metal_Results'
            call Get1_kap_fixed_metal_Results( &
               kap_lowT_z_tables, lowT_num_Zs, rq, zbar, Z, X, Rho, logRho, T, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            return
         end if
         
         upper_bdy = min(kap_blend_logT_upper_bdy, kap_lowT_z_tables(1)% x_tables(1)% logT_max)
         if (logT >= upper_bdy) then ! no lowT
            if (dbg) write(*,*) 'no lowT call Get1_kap_fixed_metal_Results'
            call Get1_kap_fixed_metal_Results( &
               kap_z_tables, num_kap_Zs, rq, zbar, Z, X, Rho, logRho, T, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (dbg) write(*,*) 'no lowT kap', kap, log10_cr(kap)
            return
         end if
         
         ! in blend region

         call Get1_kap_fixed_metal_Results( &
            kap_lowT_z_tables, lowT_num_Zs, rq, zbar, Z, X, Rho, logRho, T, logT, &
            kap2, dlnkap2_dlnRho, dlnkap2_dlnT, ierr)
         if (ierr /= 0) return
         
         call Get1_kap_fixed_metal_Results( &
            kap_z_tables, num_kap_Zs, rq, zbar, Z, X, Rho, logRho, T, logT, &
            kap1, dlnkap1_dlnRho, dlnkap1_dlnT, ierr)
         if (ierr /= 0) return

         alfa = (logT - lower_bdy) / (upper_bdy - lower_bdy)
         !alfa = 0.5 * (1.0 - cos(alfa * pi))
         beta = 1 - alfa
         
         if (dbg) then
            write(*,1) 'alfa', alfa
            write(*,1) 'beta', beta
            write(*,1) 'kap1', kap1
            write(*,1) 'kap2', kap2
            write(*,1) 'kap', kap
         end if
         
         kap = alfa*kap1 + beta*kap2
         dlnkap_dlnRho = (alfa*kap1*dlnkap1_dlnRho + beta*kap2*dlnkap2_dlnRho)/kap
         dlnkap_dlnT = (alfa*kap1*dlnkap1_dlnT + beta*kap2*dlnkap2_dlnT)/kap
         
      end subroutine Get_kap_fixed_metal_Results
      
      
      
      subroutine Get1_kap_fixed_metal_Results( &
               z_tables, num_Zs, rq, zbar, Z, X, Rho, logRho, T, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use const_def
         
         ! INPUT
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         integer, intent(in) :: num_Zs
         type (Kap_General_Info), pointer :: rq
         real, intent(in) :: zbar, Z, X
         real, intent(inout) :: Rho, logRho ! can be modified to clip to table boundaries
         real, intent(inout) :: T, logT

         ! OUTPUT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
      
         integer :: iz, i
         real :: Z0, Z1, alfa, beta, lnZ, lnZ0, lnZ1
         real :: K0, logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real(dp) :: res
         character (len=256) :: message
      
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'

         ierr = 0
         kap = 0
         dlnkap_dlnRho = 0
         dlnkap_dlnT = 0
         
         if (num_Zs > 1) then
            if (z_tables(1)% Z >= z_tables(2)% Z) then
               ierr = -3
               return
            end if
         end if         
         
         if (num_Zs == 1 .or. Z >= z_tables(num_Zs)% Z) then ! use the largest Z
            call Get_Kap_for_X( &
               z_tables, rq, num_Zs, X, logRho, logT, &
               logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
            if (ierr /= 0) return
            kap = exp10_cr(dble(logK0))
            dlnkap_dlnRho = dble(dlogK0_dlogRho)
            dlnkap_dlnT = dble(dlogK0_dlogT)
            return
         end if
         
         if (Z <= z_tables(1)% Z) then ! use the smallest Z
            if (dbg) then
               write(*,*) 'use the smallest Z'
               write(*,*) 'Z', Z
               write(*,*) 'z_tables(1)% Z', z_tables(1)% Z
            end if
            call Get_Kap_for_X( &
               z_tables, rq, 1, X, logRho, logT, &
               logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
            if (ierr /= 0) return
            kap = exp10_cr(dble(logK0))
            dlnkap_dlnRho = DBLE(dlogK0_dlogRho)
            dlnkap_dlnT = DBLE(dlogK0_dlogT)
            return
         end if
         
         if (dbg) then
            do iz = 1, num_Zs
               write(*,*) 'Z(iz)', iz, z_tables(iz)% Z
            end do
         end if

         do iz = 1, num_Zs-1
            if (Z < z_tables(iz+1)% Z) exit
         end do
      
         Z0 = z_tables(iz)% Z
         Z1 = z_tables(iz+1)% Z
         
         if (dbg) then
            write(*,*) 'Z0', Z0
            write(*,*) 'Z ', Z
            write(*,*) 'Z1', Z1
         end if
         
         if (Z1 <= Z0) then
            ierr = 1
            return
         end if
         
         if (Z <= Z0) then ! use the Z0 table
            if (dbg) then
               write(*,*) 'use the Z0 table', iz, Z0, Z, Z-Z0
               do i = 1, z_tables(iz)% num_Xs
                  write(*,*) 'X', i, z_tables(iz)% x_tables(i)% X
               end do
            end if
            call Get_Kap_for_X( &
               z_tables, rq, iz, X, logRho, logT, &
               logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
            if (ierr /= 0) return
            kap = exp10_cr(dble(logK0))
            dlnkap_dlnRho = DBLE(dlogK0_dlogRho)
            dlnkap_dlnT = DBLE(dlogK0_dlogT)
            return
         end if
         
         if (Z >= Z1) then ! use the Z1 table
            if (dbg) write(*,*) 'use the Z1 table', Z1
            call Get_Kap_for_X( &
               z_tables, rq, iz+1, X, logRho, logT, &
               logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
            if (ierr /= 0) return
            kap = exp10_cr(dble(logK0))
            dlnkap_dlnRho = DBLE(dlogK0_dlogRho)
            dlnkap_dlnT = DBLE(dlogK0_dlogT)
            return
         end if
         
         if (num_Zs >= 4 .and. rq% cubic_interpolation_in_Z) then
            if (dbg) write(*,*) 'call Get_Kap_for_Z_cubic'
            call Get_Kap_for_Z_cubic(z_tables, num_Zs, rq, iz, Z, X, logRho, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in Get_Kap_for_Z_cubic'
               return
            end if
         else ! linear
            if (dbg) write(*,*) 'call Get_Kap_for_Z_linear'
            call Get_Kap_for_Z_linear(z_tables, rq, iz, Z, Z0, Z1, X, logRho, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in Get_Kap_for_Z_linear'
               return
            end if
         end if
         
         if (dbg) then
            write(*,1) 'final logK at X and Z', log10_cr(kap), logT, logRho, X, Z
            write(*,*)
         end if

      end subroutine Get1_kap_fixed_metal_Results
      
      
      ! use tables iz-1 to iz+2 to do piecewise monotonic cubic interpolation in Z
      subroutine Get_Kap_for_Z_cubic( &
            z_tables, num_Zs, rq, iz, Z, X, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector, interp_pm

         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: num_Zs, iz
         real, intent(in) :: Z, X
         real, intent(inout) :: logRho, logT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         
         integer, parameter :: n_old = 4, n_new = 1
         real, dimension(n_old) :: logKs, dlogKs_dlogRho, dlogKs_dlogT
         real(dp) :: logK, z_old(n_old), z_new(n_new)
         integer :: i, i1, izz
         real(dp), target :: work_ary(n_old*pm_work_size)
         real(dp), pointer :: work(:)
         
         logical, parameter :: dbg = .false.
         
         11 format(a40,e20.10)
         
         ierr = 0
         work => work_ary
         
         if (iz+2 > num_Zs) then
            i1 = num_Zs-2
         else if (iz == 1) then
            i1 = 2
         else
            i1 = iz
         end if
         
         if (dbg) write(*,*) 'iz', iz
         
         do i=1,n_old
            izz = i1-2+i
            if (dbg) write(*,*) 'izz', izz
            z_old(i) = z_tables(izz)% Z
            call Get_Kap_for_X( &
               z_tables, rq, izz, X, logRho, logT, &
               logKs(i), dlogKs_dlogRho(i), dlogKs_dlogT(i), ierr)
            if (ierr /= 0) then
               if (dbg) write(*,11) 'logRho', logRho
               if (dbg) write(*,11) 'logT', logT
               return
            end if
         end do
         z_new(1) = dble(Z)
         
         call interp1(logKs, logK, ierr)
         if (ierr /= 0) then
            stop 'failed in interp1 for logK'
            return
         end if
         kap = exp10_cr(logK)
         
         call interp1(dlogKs_dlogRho, dlnkap_dlnRho, ierr)
         if (ierr /= 0) then
            stop 'failed in interp1 for dlogK_dlogRho'
            return
         end if
                  
         call interp1(dlogKs_dlogT, dlnkap_dlnT, ierr)
         if (ierr /= 0) then
            stop 'failed in interp1 for dlogK_dlogT'
            return
         end if
         
         if (dbg) then
         
            do i=1,n_old
               write(*,*) 'z_old(i)', z_old(i)
            end do
            write(*,*) 'z_new(1)', z_new(1)
            write(*,*)
         
            do i=1,n_old
               write(*,*) 'logKs(i)', logKs(i)
            end do
            write(*,*) 'logK', logK
            write(*,*)

            do i=1,n_old
               write(*,*) 'dlogKs_dlogRho(i)', dlogKs_dlogRho(i)
            end do
            write(*,*) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,*)
         
            do i=1,n_old
               write(*,*) 'dlogKs_dlogT(i)', dlogKs_dlogT(i)
            end do
            write(*,*) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,*)
            
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: interp1
#endif      
         subroutine interp1(old, new, ierr)
            real, intent(in) :: old(n_old)
            real(dp), intent(out) :: new
            integer, intent(out) :: ierr
            real(dp) :: v_old(n_old), v_new(n_new)
            v_old(:) = dble(old(:))
            call interpolate_vector( &
                  n_old, z_old, n_new, z_new, v_old, v_new, interp_pm, pm_work_size, work, &
                  'Get_Kap_for_Z_cubic', ierr)
            new = v_new(1)
         end subroutine interp1
      
      end subroutine Get_Kap_for_Z_cubic
      
      
      subroutine Get_Kap_for_Z_linear(z_tables, rq, iz, Z, Z0, Z1, X, logRho, logT, &
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real, intent(in) :: Z, Z0, Z1, X
         real, intent(inout) :: logRho, logT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr

         real :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real :: alfa, beta
      
         logical, parameter :: dbg = .false.

         ierr = 0
         call Get_Kap_for_X( &
            z_tables, rq, iz, X, logRho, logT, &
            logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'logK0', logK0
      
         call Get_Kap_for_X( &
            z_tables, rq, iz+1, X, logRho, logT, &
            logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'logK1', logK1

         ! Z0 result in logK0, Z1 result in logK1
         beta = (Z - Z1) / (Z0 - Z1) ! beta -> 1 as Z -> Z0
         alfa = 1 - beta
         
         kap = exp10_cr(dble(beta * logK0 + alfa * logK1))
         dlnkap_dlnRho = DBLE(beta * dlogK0_dlogRho + alfa * dlogK1_dlogRho)
         dlnkap_dlnT = DBLE(beta * dlogK0_dlogT + alfa * dlogK1_dlogT)

      end subroutine Get_Kap_for_Z_linear
      
      
      subroutine Get_Kap_for_X(z_tables, rq, iz, X, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz
         real, intent(in) :: X
         real, intent(inout) :: logRho, logT
         real, intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         real :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real :: X0, X1, alfa, beta
         integer :: ix, i, num_Xs
      
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         ierr = 0
         x_tables => z_tables(iz)% x_tables
         num_Xs = z_tables(iz)% num_Xs
                  
         if (X < 0 .or. X > 1) then
            ierr = -3
            return
         end if
         
         if (num_Xs > 1) then
            if (x_tables(1)% X >= x_tables(2)% X) then
               ierr = -3
               write(*,*) 'x_tables must have increasing X values for Get_Kap_for_X'
               write(*,*) 'lowT', z_tables(iz)% lowT_flag
               write(*,2) 'Z', iz, z_tables(iz)% Z
               write(*,2) 'num_Xs', num_Xs
               do i=1,num_Xs
                  write(*,2) 'X', i, x_tables(i)% X
               end do
               stop
               
               return
            end if
         end if
         
         if (num_Xs == 1 .or. X <= x_tables(1)% X) then ! use the first X
            !write(*,*) 'use the first X'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, 1, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (X >= x_tables(num_Xs)% X) then ! use the last X
            if (dbg) write(*,*) 'use the last X: call Get_Kap_for_logRho_logT'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, num_Xs, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return    
         end if

         ! search for the X
         !write(*,*) 'search for the X'
         ix = num_Xs
         do i = 1, num_Xs-1
            if (X < x_tables(i+1)% X) then
               ix = i; exit
            end if
         end do
         
         if (ix == num_Xs) then
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, num_Xs, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         X0 = x_tables(ix)% X
         X1 = x_tables(ix+1)% X
         
         if (X1 <= X0) then
            ierr = 1
            return
         end if
         
         if (X0 >= X) then ! use the X0 table
            !write(*,*) 'use the X0 table'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (X1 <= X) then ! use the X1 table
            !write(*,*) 'use the X1 table'
            call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix+1, &
                  logRho, logT, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if

         if (num_Xs >= 4 .and. rq% cubic_interpolation_in_X) then
            !write(*,*) 'call Get_Kap_for_X_cubic'
            call Get_Kap_for_X_cubic( &
               z_tables, rq, iz, ix, num_Xs, x_tables, X, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in Get_Kap_for_X_cubic'
               return
            end if
         else ! linear
            !write(*,*) 'call Get_Kap_for_X_linear'
            call Get_Kap_for_X_linear( &
               z_tables, rq, iz, ix, x_tables, X, X0, X1, logRho, logT, &
               logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in Get_Kap_for_X_linear'
               return
            end if
         end if
         
         if (.false.) then
            write(*,1) 'logK at X for Z', logK, logT, logRho, X, z_tables(iz)% Z
            write(*,*)
         end if
         
      end subroutine Get_Kap_for_X
      
      
      ! use tables ix and ix+1 to do linear interpolation in X
      subroutine Get_Kap_for_X_linear( &
            z_tables, rq, iz, ix, x_tables, X, X0, X1, logRho, logT, &
            logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz, ix
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         real, intent(in) :: X, X0, X1
         real, intent(inout) :: logRho, logT
         real, intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         real :: logK0, dlogK0_dlogRho, dlogK0_dlogT, logK1, dlogK1_dlogRho, dlogK1_dlogT
         real :: alfa, beta
         
         ierr = 0
         
         call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix, &
                  logRho, logT, logK0, dlogK0_dlogRho, dlogK0_dlogT, ierr)
         if (ierr /= 0) return
      
         call Get_Kap_for_logRho_logT( &
                  z_tables, rq, iz, x_tables, ix+1, &
                  logRho, logT, logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) return
         
         ! X0 result in logK0, X1 result in logK1
         beta = (X - X1) / (X0 - X1) ! beta -> 1 as X -> X0
         alfa = 1 - beta
         
         logK = beta * logK0 + alfa * logK1
         dlogK_dlogRho = beta * dlogK0_dlogRho + alfa * dlogK1_dlogRho
         dlogK_dlogT = beta * dlogK0_dlogT + alfa * dlogK1_dlogT
      
      end subroutine Get_Kap_for_X_linear
      
      
      ! use tables ix-1 to ix+2 to do piecewise monotonic cubic interpolation in X
      subroutine Get_Kap_for_X_cubic( &
            z_tables, rq, iz, ix, num_Xs, x_tables, X, logRho, logT, &
            logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_def
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interpolate_vector, interp_pm

         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         integer, intent(in) :: iz, ix, num_Xs
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         real, intent(in) :: X
         real, intent(inout) :: logRho, logT
         real, intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         integer, parameter :: n_old = 4, n_new = 1
         real, dimension(n_old) :: logKs, dlogKs_dlogRho, dlogKs_dlogT
         real(dp) :: x_old(n_old), x_new(n_new)
         real(dp), target :: work_ary(n_old*pm_work_size)
         real(dp), pointer :: work(:)
         integer :: i, i1, ixx
         
         logical, parameter :: dbg = .false.
         
         11 format(a40,e20.10)
         
         ierr = 0
         work => work_ary
         
         if (ix+2 > num_Xs) then
            i1 = num_Xs-2
         else if (ix == 1) then
            i1 = 2
         else
            i1 = ix
         end if
         
         if (dbg) write(*,*) 'ix', ix
         
         do i=1,n_old
            ixx = i1-2+i
            if (dbg) write(*,*) 'ixx', ixx
            x_old(i) = x_tables(ixx)% X
            call Get_Kap_for_logRho_logT( &
                     z_tables, rq, iz, x_tables, ixx, &
                     logRho, logT, logKs(i), dlogKs_dlogRho(i), dlogKs_dlogT(i), ierr)
            if (ierr /= 0) then
               if (dbg) write(*,11) 'logRho', logRho
               if (dbg) write(*,11) 'logT', logT
               return
            end if
         end do
         x_new(1) = dble(X)
         
         call interp1(logKs, logK, ierr)
         if (ierr /= 0) then
            stop 'failed in interp1 for logK'
            return
         end if
         
         call interp1(dlogKs_dlogRho, dlogK_dlogRho, ierr)
         if (ierr /= 0) then
            stop 'failed in interp1 for dlogK_dlogRho'
            return
         end if
                  
         call interp1(dlogKs_dlogT, dlogK_dlogT, ierr)
         if (ierr /= 0) then
            stop 'failed in interp1 for dlogK_dlogT'
            return
         end if
         
         if (dbg) then
         
            do i=1,n_old
               write(*,*) 'x_old(i)', x_old(i)
            end do
            write(*,*) 'x_new(1)', x_new(1)
            write(*,*)
         
            do i=1,n_old
               write(*,*) 'logKs(i)', logKs(i)
            end do
            write(*,*) 'logK', logK
            write(*,*)

            do i=1,n_old
               write(*,*) 'dlogKs_dlogRho(i)', dlogKs_dlogRho(i)
            end do
            write(*,*) 'dlogK_dlogRho', dlogK_dlogRho
            write(*,*)
         
            do i=1,n_old
               write(*,*) 'dlogKs_dlogT(i)', dlogKs_dlogT(i)
            end do
            write(*,*) 'dlogK_dlogT', dlogK_dlogT
            write(*,*)
            
         end if
         
         contains
         
#ifdef offload
         !dir$ attributes offload: mic :: interp1
#endif      
         subroutine interp1(old, new, ierr)
            real, intent(in) :: old(n_old)
            real, intent(out) :: new
            integer, intent(out) :: ierr
            real(dp) :: v_old(n_old), v_new(n_new)
            v_old(:) = dble(old(:))
            call interpolate_vector( &
                  n_old, x_old, n_new, x_new, v_old, v_new, interp_pm, pm_work_size, work, &
                  'Get_Kap_for_X_cubic', ierr)
            new = real(v_new(1))
         end subroutine interp1
      
      end subroutine Get_Kap_for_X_cubic

      
      subroutine Get_Kap_for_logRho_logT( &
               z_tables, rq, iz, x_tables, ix, &
               logRho, logT_in, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use load_kap, only: load_one
         use kap_def
         type (Kap_Z_Table), dimension(:), pointer :: z_tables
         type (Kap_General_Info), pointer :: rq
         type (Kap_X_Table), dimension(:), pointer :: x_tables
         integer :: iz, ix
         real, intent(inout) :: logRho, logT_in
         real, intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr

         real :: logR0, logR1, logT, logT0, logT1, logR, df_dx, df_dy
         real :: logR_blend1, logR_blend2, fac, logKap_es, logR_in
         real, parameter :: logR_blend_width = 0.5
         real :: Zbase, X, dXC, dXO, logKcond, dlogT, alfa, beta
         integer :: iR, jtemp, i, num_logRs, num_logTs
         logical :: clipped_logT
         logical, parameter :: read_later = .false.
      
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         ierr = 0
         if (x_tables(ix)% not_loaded_yet) then ! avoid doing critical section if possible
!$omp critical (load_kap_x_table)
            if (x_tables(ix)% not_loaded_yet) then
               call load_one( &
                  z_tables, iz, ix, dble(x_tables(ix)% X), dble(x_tables(ix)% Z), read_later, ierr)
            end if
!$omp end critical (load_kap_x_table)
         end if
         if (ierr /= 0) return
         
         logT = logT_in
         clipped_logT = .false.
         
         logR = logRho -3*logT + 18
         logR_in = logR
         
         num_logRs = x_tables(ix)% num_logRs
         num_logTs = x_tables(ix)% num_logTs
         
         if (logR > 1 .and. logT < rq% min_logT_for_logR_gt_1 .and. &
               x_tables(ix)% logT_min >= 2.5) then
            logT = rq% min_logT_for_logR_gt_1
            clipped_logT = .true.
         end if
         
         logR_blend1 = -20 !x_tables(ix)% logR_min
         logR_blend2 = logR_blend1 + logR_blend_width
         
         if (logR_in >= logR_blend1) then ! evaluate using table
         
            if (num_logRs <= 0) then
               write(*,*) 'num_logRs', num_logRs
               write(*,*) 'ix', ix
               write(*,*) 'x_tables(ix)% not_loaded_yet', x_tables(ix)% not_loaded_yet
               stop 'Get_Kap_for_logRho_logT'
            end if
          
            call Locate_logR( &
               rq, num_logRs, x_tables(ix)% logR_min, x_tables(ix)% logR_max, &
               x_tables(ix)% ili_logRs, x_tables(ix)% logRs, logR, iR, logR0, logR1, ierr)
            if (ierr /= 0) return
        
            if (num_logTs <= 0) then
               write(*,*) 'num_logTs', num_logRs
               stop 'Get_Kap_for_logRho_logT'
            end if
         
            call Locate_logT( &
               rq, num_logTs, x_tables(ix)% logT_min, x_tables(ix)% logT_max, &
               x_tables(ix)% ili_logTs, x_tables(ix)% logTs, logT, jtemp, logT0, logT1, ierr)
            if (ierr /= 0) return
            
            call Do_Kap_Interpolations( &
                  x_tables(ix)% kap1, num_logRs, num_logTs, &
                  iR, jtemp, logR0, logR, logR1, logT0, logT, logT1, &
                  logK, df_dx, df_dy)
            
            if (dbg) write(*,1) 'Do_Kap_Interpolations: logK', logK
            ! convert df_dx and df_dy to dlogK_dlogRho, dlogK_dlogT
            dlogK_dlogRho = df_dx
            if (clipped_logT) then
               dlogK_dlogT = 0   
            else
               dlogK_dlogT = df_dy - 3*df_dx
            end if      
            
         end if
         
         if (logR_in <= logR_blend2) then ! evaluate electron scattering
         
            logKap_es = log10_cr(dble(0.2*(1 + x_tables(ix)% X)))
            if (dbg) write(*,1) 'logKap_es', logKap_es

            !write(*,*) '    logR_in', logR_in
            !write(*,*) 'logR_blend1', logR_blend1
            !write(*,*) 'logR_blend2', logR_blend2
            !write(*,*)
            
            if (logR_in <= logR_blend1) then
               
               logK = logKap_es
               dlogK_dlogRho = 0
               dlogK_dlogT = 0
               
            else
               
               fac = (logR_in - logR_blend1) / (logR_blend2 - logR_blend1)
               if (dbg) then
                  write(*,1) 'logR_in', logR_in
                  write(*,1) 'logR_blend1', logR_blend1
                  write(*,1) 'logR_blend2', logR_blend2
                  write(*,1) 'fac', fac
               end if
               fac = 0.5*(1 - real(cospi_cr(dble(fac))))
               
               !write(*,*) '      fac', fac
               !write(*,*) 'logKap_es', logKap_es
               !write(*,*) '     logK', logK
               !write(*,*) '      new', fac*logK + (1-fac)*logKap_es
               !write(*,*)

               logK = fac*logK + (1-fac)*logKap_es
               dlogK_dlogRho = fac*dlogK_dlogRho
               dlogK_dlogT = fac*dlogK_dlogRho
               if (dbg) then
                  write(*,1) 'fac', fac
                  write(*,1) 'fac*logK', fac*logK
                  write(*,1) '(1-fac)', (1-fac)
                  write(*,1) '(1-fac)*logKap_es', (1-fac)*logKap_es
               end if
               
            end if
         
         end if
         
         if (dbg) then
            write(*,1) 'logK', logK, logT, logRho, x_tables(ix)% X, z_tables(iz)% Z
         end if
         
      end subroutine Get_Kap_for_logRho_logT
      
#ifdef offload
      !dir$ end options
#endif

      end module kap_eval_fixed
