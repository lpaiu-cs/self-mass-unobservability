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

      module opacities
      
      use star_private_def
      use utils_lib,only:is_bad_num
      use const_def
      
      implicit none
      
      logical, parameter :: dbg = .false.
      

      contains
      
      
      subroutine set_kap_params(s,ierr)
         use kap_lib, only: kap_set_choices
         type (star_info), pointer :: s         
         integer, intent(out) :: ierr
         call kap_set_choices( &
            s% kap_handle, s% cubic_interpolation_in_X, s% cubic_interpolation_in_Z, &
            s% include_electron_conduction, &
            s% kap_Type2_full_off_X, s% kap_Type2_full_on_X, &
            s% kap_Type2_full_off_dZ, s% kap_Type2_full_on_dZ, &
            ierr)
      end subroutine set_kap_params
      
      
      subroutine do_opacities(s,nzlo,nzhi,ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         call set_kap_params(s,ierr)
         if (ierr /= 0) return
         if (dbg) write(*,*) 'call do_kap'
         call do_kap(s,nzlo,nzhi,ierr)
      end subroutine do_opacities
      
      
      subroutine do_kap(s,nzlo,nzhi,ierr)
         use star_utils, only: foreach_cell
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         logical, parameter :: use_omp = .false.
         if (s% use_other_opacity_factor) then
            call s% other_opacity_factor(s% id, ierr)
            if (ierr /= 0) return
         else
            s% extra_opacity_factor(1:s% nz) = s% opacity_factor
         end if
         call foreach_cell(s, nzlo, nzhi, use_omp, do_kap_for_cell, ierr)
      end subroutine do_kap
      
      
      subroutine do_kap_for_cell(s,k,ierr)
         use const_def,only:ln10
         use utils_lib,only:is_bad_num
         use net_def,only:net_general_info
         use rates_def, only:i_rate
         use chem_def
         use kap_def
         use kap_lib
         use eos_def, only : i_lnfree_e
         use eos_lib, only: Radiation_Pressure
         use star_utils, only: get_XYZ
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer, intent(out) :: ierr

         integer, pointer :: net_iso(:)
         integer :: i, iz, kh
         real(dp) :: &
               log10_rho, log10_T, dlnkap_dlnd, dlnkap_dlnT, &
               opacity_max, opacity_max0, opacity_max1, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               log_r, log_r_in, log_r_out, log_r_frac, frac, min_cno_for_kap_limit, &
               P, Prad, Pgas, Pgas_div_P, Ledd_factor, Ledd_kap, Ledd_log, &
               a, b, da_dlnd, da_dlnT, db_dlnd, db_dlnT, opacity_factor
         character (len=100) :: message
         real(dp), pointer :: xa(:)
         
         include 'formats'
                  
         ierr = 0
         
         log10_rho = s% lnd(k)/ln10
         log10_T = s% lnT(k)/ln10
         
         lnfree_e = s% lnfree_e(k)
         d_lnfree_e_dlnRho = s% d_eos_dlnd(i_lnfree_e,k)
         d_lnfree_e_dlnT = s% d_eos_dlnT(i_lnfree_e,k)
         
         Pgas_div_P = 1
         xa(1:s% species) => s% xa(1:s% species,k)
         s% opacity(k) = get1_kap( &
            s, k, s% zbar(k), xa, s% q(k), Pgas_div_P, &
            log10_rho, log10_T, fraction_of_op_mono(s,k), &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            s% kap_frac_Type2(k), dlnkap_dlnd, dlnkap_dlnT, ierr)
            
         if (.false.) then
!$omp critical (star_kap_get)
            call show_stuff
            stop 'debug: do_kap_for_cell'
!$omp end critical (star_kap_get)
         end if
         
         if (ierr /= 0 .or. is_bad_num(s% opacity(k))) then
            if (s% report_ierr) then
               return
!$omp critical (star_kap_get)
               write(*,*) 'do_kap_for_cell: kap_get failure for cell ', k
               call show_stuff
               stop 'debug: do_kap_for_cell'
!$omp end critical (star_kap_get)
            end if
            ierr = -1
            return
         end if
         
         if (s% opacity(k) < 1d-99) then
            s% opacity(k) = 1d-99
            dlnkap_dlnd = 0
            dlnkap_dlnT = 0
         end if
         
         opacity_factor = s% extra_opacity_factor(k)
         if (s% min_logT_for_opacity_factor_off > 0) then
            if (log10_T >= s% max_logT_for_opacity_factor_off .or. &
                log10_T <= s% min_logT_for_opacity_factor_off) then
               opacity_factor = 1
            else if (log10_T > s% max_logT_for_opacity_factor_on) then
               opacity_factor = 1 + (opacity_factor-1)* &
                  (log10_T - s% max_logT_for_opacity_factor_off)/ &
                  (s% max_logT_for_opacity_factor_on - s% max_logT_for_opacity_factor_off)
            else if (log10_T < s% min_logT_for_opacity_factor_on) then
               opacity_factor = 1 + (opacity_factor-1)* &
                  (log10_T - s% min_logT_for_opacity_factor_off)/ &
                  (s% min_logT_for_opacity_factor_on - s% min_logT_for_opacity_factor_off)
            end if
         end if
         
         s% opacity(k) = s% opacity(k)*opacity_factor
         s% d_opacity_dlnd(k) = s% opacity(k)*dlnkap_dlnd
         s% d_opacity_dlnT(k) = s% opacity(k)*dlnkap_dlnT

         if (s% opacity(k) > s% opacity_max .and. s% opacity_max > 0) then
            s% opacity(k) = s% opacity_max
            s% d_opacity_dlnd(k) = 0
            s% d_opacity_dlnT(k) = 0
         end if         
         
         if (.false. .and. k > 365 .and. k < 374 .and. s% model_number == 761) then
            write(*,2) 'logKap = ', k, log10_cr(s% opacity(k))
            !call show_stuff
            if (k == 373) write(*,*)
            !if (k > 374) stop 'debug: do_kap_for_cell'
         end if
         
         if (k == -773) then
            call show_stuff
            !stop 'kap'
         end if

         
         contains
         
         subroutine show_stuff
            include 'formats'
            real(dp) :: xc, xo, xh, xhe, Z
            integer :: i, iz
            call get_XYZ(s, s% xa(:,k), xh, xhe, Z)
            xc = 0; xo = 0
            do i=1, s% species
               iz = floor(chem_isos% Z(s% chem_id(i)) + 1d-10)
               select case(iz)
                  case (6)
                     xc = xc + s% xa(i,k)
                  case (8)
                     xo = xo + s% xa(i,k)
               end select
            end do
            write(*,2) 'show opacity info'
            write(*,2) 'logT = ', k, log10_T
            write(*,2) 'logRho = ', k, log10_rho
            write(*,2) 'z = ', k, z
            write(*,2) 'xh = ', k, xh
            write(*,2) 'xc = ', k, xc
            write(*,2) 'xo = ', k, xo
            write(*,2) 'lnfree_e = ', k, lnfree_e
            write(*,2) 'd_lnfree_e_dlnRho = ', k, d_lnfree_e_dlnRho
            write(*,2) 'd_lnfree_e_dlnT = ', k, d_lnfree_e_dlnT
            write(*,2) 'abar = ', k, s% abar(k)
            write(*,2) 'zbar = ', k, s% zbar(k)
            write(*,*)
            write(*,2) 'rho = ', k, s% rho(k)
            write(*,2) 'lnrho = ', k, s% lnd(k)
            write(*,2) 'T = ', k, s% T(k)
            write(*,2) 'lnT = ', k, s% lnT(k)
            write(*,2) 'logKap = ', k, log10_cr(s% opacity(k))
            write(*,2) 'opacity = ', k, s% opacity(k)
            write(*,2) 'dlnkap_dlnd = ', k, dlnkap_dlnd
            write(*,2) 'dlnkap_dlnT = ', k, dlnkap_dlnT
            write(*,2) 'd_opacity_dlnd = ', k, s% d_opacity_dlnd(k)
            write(*,2) 'd_opacity_dlnT = ', k, s% d_opacity_dlnT(k)
            write(*,*)
         end subroutine show_stuff
         
      end subroutine do_kap_for_cell
      
      
      subroutine eval_kap_Type2( &
            s, k, zbar, X, Z, Zbase, XC, XN, XO, XNe, log10_rho, log10_T, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         use kap_lib
         type (star_info), pointer :: s   
         integer, intent(in) :: k      
         real(dp), intent(in) :: zbar, X, Z, Zbase, XC, XN, XO, XNe, &
            log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(in), pointer :: xa(:)
         real(dp), intent(out) :: frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
         if (s% use_other_kap) then
            call s% other_kap_get_Type2( &
               s% id, k, s% kap_handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
               log10_rho, log10_T, s% species, s% chem_id, s% net_iso, xa, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         else
            call kap_get_Type2( &
               s% kap_handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, & 
               log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         end if
      end subroutine eval_kap_Type2
      
      
      subroutine eval_kap_Type1( &
            s, k, zbar, X_in, Z_in, q, Pgas_div_P, log10_rho, log10_T, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         use kap_lib
         type (star_info), pointer :: s         
         integer, intent(in) :: k      
         real(dp), intent(in) :: zbar, X_in, Z_in, log10_rho, log10_T, xa(:), &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, q, Pgas_div_P
         real(dp), intent(out) :: kap, dln_kap_dlnRho, dln_kap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
         real(dp) :: X, Z
         include 'formats'
         X = X_in
         Z = Z_in
         if (s% use_other_kap) then
            call s% other_kap_get_Type1( &
               s% id, k, s% kap_handle, zbar, X, Z, log10_rho, log10_T, & 
               s% species, s% chem_id, s% net_iso, xa, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         else
            call kap_get_Type1( &
               s% kap_handle, zbar, X, Z, log10_rho, log10_T, & 
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
            if (k == -1) then
               write(*,2) 'zbar', k, zbar
               write(*,2) 'X', k, X
               write(*,2) 'Z', k, Z
               write(*,2) 'log10_rho', k, log10_rho
               write(*,2) 'log10_T', k, log10_T
               write(*,2) 'lnfree_e', k, lnfree_e
               write(*,2) 'd_lnfree_e_dlnRho', k, d_lnfree_e_dlnRho
               write(*,2) 'd_lnfree_e_dlnT', k, d_lnfree_e_dlnT
               write(*,2) 'kap', k, kap
               write(*,2) 'dln_kap_dlnRho', k, dln_kap_dlnRho
               write(*,2) 'dln_kap_dlnT', k, dln_kap_dlnT
            end if
         end if
      end subroutine eval_kap_Type1
      
      
      real(dp) function fraction_of_op_mono(s, k) result(beta)
         type (star_info), pointer :: s
         integer, intent(in) :: k         
         real(dp) :: log10_T, alfa
         if (k <= 0 .or. k > s% nz) then
            beta = 0d0
            return
         end if
         log10_T = s% lnT(k)/ln10         
         ! alfa is fraction standard opacity
         if (log10_T >= s% high_logT_op_mono_full_off .or. &
             log10_T <= s% low_logT_op_mono_full_off) then
            alfa = 1d0
         else if (log10_T <= s% high_logT_op_mono_full_on .and. &
                  log10_T >= s% low_logT_op_mono_full_on) then
            alfa = 0d0
         else if (log10_T > s% high_logT_op_mono_full_on) then
            alfa = (log10_T - s% high_logT_op_mono_full_on) / &
               (s% high_logT_op_mono_full_off - s% high_logT_op_mono_full_on)
         else
            alfa = (log10_T - s% low_logT_op_mono_full_off) / &
               (s% low_logT_op_mono_full_on - s% low_logT_op_mono_full_off)
         end if
         beta = 1d0 - alfa ! beta is fraction of op mono      
      end function fraction_of_op_mono
      
      
      real(dp) function get1_kap( &
            s, k, zbar, xa, q, Pgas_div_P, log10_rho, log10_T, frac_op_mono, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, dlnkap_dlnd, dlnkap_dlnT, ierr)
         use utils_lib
         use num_lib
         use kap_lib, only: &
            load_op_mono_data, get_op_mono_params, &
            get_op_mono_args, kap_get_op_mono
         use chem_def, only: ih1, ihe3, ihe4, chem_isos
         use star_utils, only: get_XYZ
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: &
            zbar, q, Pgas_div_P, log10_rho, log10_T, frac_op_mono, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(in), pointer :: xa(:)
         real(dp), intent(out) :: frac_Type2, dlnkap_dlnd, dlnkap_dlnT
         integer, intent(out) :: ierr

         integer :: i, iz, nptot, ipe, nrad, thread_num, sz, offset
         real(dp) :: opacity, alfa, beta, lnkap, &
            kap_op, dlnkap_op_dlnRho, dlnkap_op_dlnT, &
            Z, xh, xhe, xc, xn, xo, xne, xheavy
         integer, pointer :: net_iso(:)
         real, pointer :: &
            umesh(:), ff(:,:,:,:), rs(:,:,:), ss(:,:,:,:)
         integer :: nel, izzp(s% species)
         real(dp) :: fap(s% species), gp1(s% species)
         logical :: screening
         real(dp), parameter :: &
            eps = 1d-6, minus_eps = -eps, one_plus_eps = 1d0 + eps
         
         include 'formats'
         
         get1_kap = -1d99
         frac_Type2 = 0
         alfa = 1d0
         beta = 0d0
         
         net_iso => s% net_iso
         xc = 0; xn = 0; xo = 0; xne = 0
         do i=1, s% species
            iz = floor(chem_isos% Z(s% chem_id(i)) + 1d-10)
            select case(iz)
               case (6)
                  xc = xc + xa(i)
               case (7)
                  xn = xn + xa(i)
               case (8)
                  xo = xo + xa(i)
               case (10)
                  xne = xne + xa(i)
            end select
         end do
         
         if (xc < minus_eps .or. xn < minus_eps .or. &
             xo < minus_eps .or. xne < minus_eps .or. &
             xc > one_plus_eps .or. xn > one_plus_eps .or. &
             xo > one_plus_eps .or. xne > one_plus_eps) then
            ierr = -1
            if (s% report_ierr) then
               if (xc < 0d0 .or. xc > 1d0) write(*,2) 'get1_kap: xc', k, xc
               if (xn < 0d0 .or. xn > 1d0) write(*,2) 'get1_kap: xn', k, xn
               if (xo < 0d0 .or. xo > 1d0) write(*,2) 'get1_kap: xo', k, xo
               if (xne < 0d0 .or. xne > 1d0) write(*,2) 'get1_kap: xne', k, xne
            end if
            return
            stop 'bad mass fraction: get1_kap'
         end if
         
         call get_XYZ(s, xa, xh, xhe, Z)
         
         if (xh < 0d0 .or. xhe < 0d0 .or. Z < 0d0) then
            ierr = -1
            if (s% report_ierr) then
               if (xh < 0d0) write(*,2) 'xh', k, xh
               if (xhe < 0d0) write(*,2) 'xhe', k, xhe
               if (Z < 0d0) write(*,2) 'Z', k, Z
            end if
            return
            stop 'negative mass fraction: get1_kap'
         end if

         if (s% use_simple_es_for_kap) then
            get1_kap = 0.2d0*(1 + xh)
            dlnkap_dlnd = 0
            dlnkap_dlnT = 0
            return
         end if
         
         beta = frac_op_mono
         alfa = 1d0 - beta
                  
         if (beta > 0d0) then
            
            call get_op_mono_args( &
   		      s% species, xa, s% op_mono_min_X_to_include, s% chem_id, &
   		      nel, izzp, fap, ierr)
            if (ierr /= 0) then
               write(*,*) 'error in get_op_mono_args, ierr = ',ierr
               return
            end if
               
            if (associated(s% op_mono_umesh1)) then
            
               thread_num = utils_OMP_GET_THREAD_NUM() ! in range 0 to op_mono_n-1
               if (thread_num < 0) then
                  write(*,3) 'thread_num < 0', thread_num, s% op_mono_n
                  ierr = -1
                  return
               end if
               if (thread_num >= s% op_mono_n) then
                  write(*,3) 'thread_num >= s% op_mono_n', thread_num, s% op_mono_n
                  ierr = -1
                  return
               end if
               nptot = s% op_mono_nptot
               ipe = s% op_mono_ipe
               nrad = s% op_mono_nrad
               
               sz = nptot; offset = thread_num*sz
               umesh(1:nptot) => s% op_mono_umesh1(offset+1:offset+sz)
               sz = nptot*ipe*4*4; offset = thread_num*sz
               ff(1:nptot,1:ipe,1:4,1:4) => s% op_mono_ff1(offset+1:offset+sz)
               sz = nptot*4*4; offset = thread_num*sz
               rs(1:nptot,1:4,1:4) => s% op_mono_rs1(offset+1:offset+sz)
               sz = nptot*nrad*4*4; offset = thread_num*sz
               ss(1:nptot,1:nrad,1:4,1:4) => s% op_mono_s1(offset+1:offset+sz)
               
            else

               call load_op_mono_data( &
                  s% op_mono_data_path, s% op_mono_data_cache_filename, ierr)
               if (ierr /= 0) then
                  write(*,*) 'error while loading OP data, ierr = ',ierr
                  return
               end if

               call get_op_mono_params(nptot, ipe, nrad)
               allocate( &
                  umesh(nptot), ff(nptot,ipe,4,4), &
                  rs(nptot,4,4), ss(nptot,nrad,4,4), stat=ierr)
               if (ierr /= 0) return
               
            end if

            screening = .true.
            if (s% use_other_kap) then
               call s% other_kap_get_op_mono( &
                  s% kap_handle, zbar, log10_rho, log10_T, &
                  lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  s% use_op_mono_alt_get_kap, &
                  nel, izzp, fap, screening, umesh, ff, rs, ss, &
                  kap_op, dlnkap_op_dlnRho, dlnkap_op_dlnT, ierr)
            else
               call kap_get_op_mono( &
                  s% kap_handle, zbar, log10_rho, log10_T, &
                  lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  s% use_op_mono_alt_get_kap, &
                  nel, izzp, fap, screening, umesh, ff, rs, ss, &
                  kap_op, dlnkap_op_dlnRho, dlnkap_op_dlnT, ierr)
            end if
         
            if (.not. associated(s% op_mono_umesh1)) deallocate(umesh, ff, rs, ss)
         
            if (ierr /= 0) then             
               if (s% report_ierr) write(*,*) 'error in op_mono kap, ierr = ', ierr
               alfa = 1d0
               beta = 0d0
               ierr = 0               
            else if (beta == 1d0) then
               get1_kap = kap_op
               dlnkap_dlnT = dlnkap_op_dlnT
               dlnkap_dlnd = dlnkap_op_dlnRho
               return
            end if

         end if

         if (s% use_Type2_opacities) then
            call eval_kap_Type2( &
               s, k, zbar, xh, Z, s% Zbase, XC, XN, XO, XNe, log10_rho, log10_T, xa, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, opacity, dlnkap_dlnd, dlnkap_dlnT, ierr)
         
            if (ierr /= 0) then
               return
               write(*,*)
               write(*,*)
               write(*,2) 'ierr from eval_kap_Type2', ierr
               write(*,2) 'zbar', k, zbar
               write(*,2) 'xh', k, xh
               write(*,2) 'Z', k, Z
               write(*,2) 's% Zbase', k, s% Zbase
               write(*,2) 'XC', k, XC
               write(*,2) 'XN', k, XN
               write(*,2) 'XO', k, XO
               write(*,2) 'XNe', k, XNe
               write(*,2) 'log10_rho', k, log10_rho
               write(*,2) 'log10_T', k, log10_T
               write(*,*)
               stop 'get1_kap'
            end if

         
            if (.false.) then
               write(*,2) 'opacity', k, opacity
               write(*,2) 'zbar', k, zbar
               write(*,2) 'xh', k, xh
               write(*,2) 'Z', k, Z
               write(*,2) 's% Zbase', k, s% Zbase
               write(*,2) 'XC', k, XC
               write(*,2) 'XN', k, XN
               write(*,2) 'XO', k, XO
               write(*,2) 'XNe', k, XNe
               write(*,2) 'lnfree_e', k, lnfree_e
               write(*,2) 'log10_rho', k, log10_rho
               write(*,2) 'log10_T', k, log10_T
               write(*,2) 'beta', k, beta
               write(*,2) 'frac_Type2', k, frac_Type2
               write(*,*)
               stop 'get1_kap'
            end if


         else
            Z = min( max(1 - (xh+xhe), 0d0), 0.1d0 )
            call eval_kap_Type1( &
               s, k, zbar, xh, Z, q, Pgas_div_P, log10_rho, log10_T, xa, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               opacity, dlnkap_dlnd, dlnkap_dlnT, ierr)

         
            if (.false.) then
               write(*,2) 'opacity', k, opacity
               write(*,2) 'zbar', k, zbar
               write(*,2) 'xh', k, xh
               write(*,2) 'Z', k, Z
               write(*,2) 's% Zbase', k, s% Zbase
               write(*,2) 'log10_rho', k, log10_rho
               write(*,2) 'log10_T', k, log10_T
               write(*,2) 'beta', k, beta
               write(*,*)
               stop 'get1_kap'
            end if


         end if
         
         if (beta == 0d0) then
            get1_kap = opacity
            return
         end if
         
         lnkap = alfa*log_cr(opacity) + beta*log_cr(kap_op)
         opacity = exp_cr(lnkap)
         dlnkap_dlnT = alfa*dlnkap_dlnT + beta*dlnkap_op_dlnT
         dlnkap_dlnd = alfa*dlnkap_dlnd + beta*dlnkap_op_dlnRho         
         
         get1_kap = opacity

         if (k == s% trace_k) then
            write(*,5) 'opacity', &
               k, s% newton_iter, s% model_number, s% newton_adjust_iter, &
                        opacity
         end if
      
      end function get1_kap
      
      

      end module opacities

