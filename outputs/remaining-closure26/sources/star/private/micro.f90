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

      module micro
      
      use star_private_def
      use utils_lib,only:is_bad_num
      use const_def
      use star_utils, only: foreach_cell
      
      implicit none

      logical, parameter :: dbg = .false.
      


      contains
      
      
      subroutine set_micro_vars( &
            s, nzlo, nzhi, skip_eos, skip_net, skip_neu, skip_kap, ierr)
         use opacities, only: set_kap_params, do_kap_for_cell
         use star_utils, only: start_time, update_time
         use net_lib, only: net_work_size
         use net, only: do_net
         !use solve_burn, only: do_one_zone_burns
         use chem_def, only: icno, ipp
         use rates_def, only: i_rate
         
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: skip_eos, skip_net, skip_neu, skip_kap
         integer, intent(out) :: ierr
         
         integer :: j, k, op_err, k_bad, res
         integer :: time0, clock_rate
         real(dp) :: total, alfa, beta
         
         include 'formats'

         ierr = 0
         if (dbg) then
            write(*,*) 
            write(*,*) 'set_micro_vars'
            write(*,*) 'nzlo', nzlo
            write(*,*) 'nzhi', nzhi
            write(*,*) 'skip_net', skip_net
            write(*,*) 'skip_kap', skip_kap
         end if
         
         if (.not. skip_eos) then
            call set_eos(ierr)
            if (ierr /= 0) return
         end if
         
         do k=nzlo, nzhi
            if (k == 1) then
               s% rho_face(k) = s% rho(k)
               s% csound_at_face(1) = s% csound(1)
            else
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
               beta = 1 - alfa
               s% rho_face(k) = alfa*s% rho(k) + beta*s% rho(k-1)
               s% csound_at_face(k) = alfa*s% csound(k) + beta*s% csound(k-1)
            end if
         end do

         if (.not. (skip_kap .and. skip_neu)) then
         
            if (s% doing_timing) call start_time(s, time0, total)
            
            if (.not. skip_kap) then
               call set_kap_params(s,ierr)
               if (ierr /= 0) return
               if (s% use_other_opacity_factor) then
                  call s% other_opacity_factor(s% id, ierr)
                  if (ierr /= 0) return
               else
                  s% extra_opacity_factor(1:s% nz) = s% opacity_factor
               end if
               call setup_for_op_mono(s, ierr)
               if (ierr /= 0) return
            end if
               
!$OMP PARALLEL DO PRIVATE(k,op_err)
            do k = nzlo, nzhi
               op_err = 0
               call do1_neu_kap(s,k,op_err)
               if (op_err /= 0) ierr = op_err
            end do
!$OMP END PARALLEL DO
            if (s% doing_timing) call update_time(s, time0, total, s% time_neu_kap)         
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'do1_neu_kap returned ierr', ierr
               return
            end if
            
         end if
         
         if (ierr /= 0) return
         
         if (.not. skip_net) then
         
            
            if (dbg) then
               write(*,2) 's% split_mixing_choice', s% split_mixing_choice
               write(*,*) 's% doing_newton_iterations', s% doing_newton_iterations
            end if
            
            if (s% split_mixing_choice < 0 .and. &
                  (.not. s% split_mix_do_burn) .and. &
                  s% dt > 0d0 .and. s% doing_newton_iterations) then
                  
                  stop 'split burn no longer supported'
!               if (dbg) write(*,*) 'micro: call do_one_zone_burns'
!               res = do_one_zone_burns(s, dt, 1, 1, k_bad)
!               if (res /= keep_going) then
!                  ierr = -1
!                  if (s% report_ierr) &
!                     write(*,3) 'do_one_zone_burns failed: model, k_bad, T', &
!                        s% model_number, k_bad, s% T(k_bad)
!                  return
!               end if
            
            else

               if (dbg) write(*,*) 'micro: call do_net'
               if (s% doing_timing) call start_time(s, time0, total)
               
               call do_net(s, nzlo, nzhi, .false., ierr)
               
               if (s% doing_timing) call update_time(s, time0, total, s% time_nonburn_net)         
               if (dbg) write(*,*) 'micro: done do_net'
               if (ierr /= 0) then
                  if (s% report_ierr) write(*,*) 'do_net returned ierr', ierr
                  return
               end if
            
            end if
         
            
         end if 

         
         contains
         
         
         subroutine set_eos(ierr)
            integer, intent(out) :: ierr
            integer :: k
            real(dp) :: alfa
            include 'formats'
            ierr = 0
            if (dbg) write(*,*) 'call do_eos'
            if (s% doing_timing) call start_time(s, time0, total)
            call do_eos(s,nzlo,nzhi,ierr)
            if (s% doing_timing) call update_time(s, time0, total, s% time_eos)         
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'do_eos returned ierr', ierr
               return
            end if
         end subroutine set_eos
         
         
         subroutine debug(str) 
            use chem_def
            character (len=*), intent(in) :: str
            integer :: k, j
            include 'formats'
            k = 1469
            do j=1,1 !s% species
               if (.true. .or. s% xa(j,k) > 1d-9) &
                  write(*,1) trim(str) // ' xin(net_iso(i' // &
                     trim(chem_isos% name(s% chem_id(j))) // '))= ', s% xa(j,k)
            end do
         end subroutine debug
         
         
         subroutine do1_neu_kap(s,k,ierr)
            use neu_def,only:Tmin_neu
            use neu, only: do_neu_for_cell, do_clear_neu_for_cell
            type (star_info), pointer :: s         
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            if (s% T(k) >= Tmin_neu .and. s% non_nuc_neu_factor > 0d0) then
               call do_neu_for_cell(s,k,ierr)
            else
               call do_clear_neu_for_cell(s,k,ierr)
            end if
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'do_neu returned ierr', ierr
               return
            end if      
            if (skip_kap) return
            call do_kap_for_cell(s,k,ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'do_kap returned ierr', ierr
               return
            end if
         end subroutine do1_neu_kap


      end subroutine set_micro_vars
      
      
      subroutine setup_for_op_mono(s, ierr)
         use opacities, only: fraction_of_op_mono
         use kap_lib, only: load_op_mono_data, get_op_mono_params
         use utils_lib, only: utils_OMP_GET_MAX_THREADS

         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: k, nptot, ipe, nrad, n
         logical :: need_op
         real(dp) :: beta
         
         include 'formats'
         
         ierr = 0
         
         if (s% high_logT_op_mono_full_off < 0d0 .or. &
             s% low_logT_op_mono_full_off > 99d0) return

         if (s% op_mono_n > 0) return ! already setup
         
         need_op = .false.
         do k=1,s% nz
            beta = fraction_of_op_mono(s,k)
            if (beta > 0d0) then
               need_op = .true.
               exit
            end if
         end do
         if (.not. need_op) return

         call load_op_mono_data( &
            s% op_mono_data_path, s% op_mono_data_cache_filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'error while loading OP data, ierr = ',ierr
            return
         end if

         call get_op_mono_params(nptot, ipe, nrad)

         n = utils_OMP_GET_MAX_THREADS()
         
         if (n /= s% op_mono_n .or. &
             nptot /= s% op_mono_nptot .or. &
             ipe /= s% op_mono_ipe .or. &
             nrad /= s% op_mono_nrad) then
            if (associated(s% op_mono_umesh1)) deallocate(s% op_mono_umesh1)
            if (associated(s% op_mono_ff1)) deallocate(s% op_mono_ff1)
            if (associated(s% op_mono_rs1)) deallocate(s% op_mono_rs1)
            if (associated(s% op_mono_s1)) deallocate(s% op_mono_s1)
            allocate( &
               s% op_mono_umesh1(nptot*n), s% op_mono_ff1(nptot*ipe*4*4*n), &
               s% op_mono_rs1(nptot*4*4*n), s% op_mono_s1(nptot*nrad*4*4*n), &
               stat=ierr)
            if (ierr /= 0) return
            s% op_mono_n = n
            s% op_mono_nptot = nptot
            s% op_mono_ipe = ipe
            s% op_mono_nrad = nrad
         end if
      
      end subroutine setup_for_op_mono
      
      
      subroutine do_eos(s,nzlo,nzhi,ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         logical, parameter :: use_omp = .true.
         include 'formats'
         integer :: k
         ierr = 0
         if (dbg) write(*,*) 'do_eos call foreach_cell', nzlo, nzhi
         
         !write(*,*) 'omp off for eos'
         !call foreach_cell(s,nzlo,nzhi,.false.,do_eos_for_cell,ierr)
         !return
         
         call foreach_cell(s,nzlo,nzhi,use_omp,do_eos_for_cell,ierr)

      end subroutine do_eos
      
      
      subroutine do_eos_for_cell(s,k,ierr)
         use const_def
         use chem_def
         use chem_lib
         use eos_def
         use eos_lib
         use net_def,only:net_general_info
         use utils_lib,only:is_bad_num

         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         real(dp), dimension(num_eos_basic_results) :: &
            res, res_a, res_b, d_dlnd, d_dlnT
         real(dp) :: &
            z, xh, xhe3, xhe4, xhe, sumx, dx, dxh_a, dxh_b, &
            Rho, log10Rho, lnd, lnE, log10T, T, energy
         integer, pointer :: net_iso(:)
         integer :: j, h1, h2, he3, he4, species
         real(dp), parameter :: epsder = 1d-4, Z_limit = 0.5d0
         
         logical, parameter :: testing = .false.
            
         include 'formats'
         
         !write(*,3) 'do_eos_for_cell', k, s% nz
         
         ierr = 0
         
         net_iso => s% net_iso
         species = s% species
         call basic_composition_info( &
            species, s% chem_id, s% xa(1:species,k), s% X(k), s% Y(k), &
            s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), s% mass_correction(k), sumx)
         
         if (.false. .and. k == 1469) then
            j = 2
            write(*,1) 'do_eos_for_cell xin(net_iso(i' // &
               trim(chem_isos% name(s% chem_id(j))) // '))= ', &
                     s% xa(j,k), s% abar(k)
         end if
         
         h1 = net_iso(ih1)
         h2 = net_iso(ih2)
         he3 = net_iso(ihe3)
         he4 = net_iso(ihe4)
         
         if (h1 /= 0) then
            xh = s% xa(h1,k)
         else
            xh = 0
         end if
         if (h2 /= 0) xh = xh + s% xa(h2,k)
         
         if (he3 /= 0) then
            xhe3 = s% xa(he3,k)
         else
            xhe3 = 0
         end if
         
         if (he4 /= 0) then
            xhe4 = s% xa(he4,k)
         else
            xhe4 = 0
         end if
         
         xhe = xhe3 + xhe4
         z = max(0d0,1d0-(xh+xhe))
         
         if (.false. .and. z < s% initial_z - 1d-6) then
            write(*,2) 'z dropped', k, z
            stop 'do_eos_for_cell'
         end if
         
         if (s% lnPgas_flag) then
            if (s% Pgas(k) == 0) then
               ierr = -1
               return               
               write(*,2) 's% lnPgas(k)/ln10', k, s% lnPgas(k)/ln10
               stop 'do_eos_for_cell s% Pgas(k) == 0'
            end if        
            call eval_eosPT( &
               s, k, z, xh, s% abar(k), s% zbar(k), s% xa(:,k), &
               s% Pgas(k), s% lnPgas(k)/ln10, s% T(k), s% lnT(k)/ln10, & 
               Rho, log10Rho, s% dlnRho_dlnPgas_const_T(k), s% dlnRho_dlnT_const_Pgas(k), &
               res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
               s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), ierr)
            s% lnd(k) = log10Rho*ln10
            s% rho(k) = Rho
            call store_stuff(ierr)
            return            
         end if
         
         if (s% E_flag) then            
            lnd = s% lnd(k)
            rho = exp_cr(lnd)
            lnE = s% lnE(k)
            energy = s% energy(k)            
            call eval_eosDE( &
               s, k, z, xh, s% abar(k), s% zbar(k), s% xa(:,k), &
               rho, lnd/ln10, lnE/ln10, & 
               T, log10T, res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
               s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), &
               s% dlnT_dlnE_c_Rho(k), s% dlnT_dlnd_c_E(k), &
               s% dlnPgas_dlnE_c_Rho(k), s% dlnPgas_dlnd_c_E(k), &
               ierr)
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,*) 'do_eos_for_cell: eval_eosDE ierr', ierr
                  !stop 'do_eos_for_cell'
               end if
               return
            end if
            s% T(k) = T
            s% lnT(k) = log10T*ln10
            s% lnPgas(k) = res(i_lnPgas)
            s% Pgas(k) = exp_cr(s% lnPgas(k))
            call store_stuff(ierr)
            s% lnE(k) = lnE
            s% energy(k) = energy
            call set_dlnP_dlnd_c_E(s,k)
            call set_dlnP_dlnE_c_Rho(s,k)
            return            
         end if
         
         call eos_get( &
            s, k, z, xh, s% abar(k), s% zbar(k), s% xa(:,k), &
            s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
            res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
            s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), ierr)
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,*) 'do_eos_for_cell: eos_get ierr', ierr
               stop 'do_eos_for_cell'
            end if
            return
         end if
         s% lnPgas(k) = res(i_lnPgas)
         s% Pgas(k) = exp_cr(s% lnPgas(k))
         call store_stuff(ierr)
         
         
         contains
         
         subroutine store_stuff(ierr)  
            integer, intent(out) :: ierr
            include 'formats'
            call store_eos_for_cell(s, k, res, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,2) 'store_eos_for_cell failed', k
   !               return
                  call write_info
                  write(*,*) 's% E_flag', s% E_flag
                  stop 'debug: eos'
               end if
               return
            end if
            if (k == s% trace_k) then
               write(*,5) 'grada', k, s% newton_iter, s% newton_adjust_iter, &
                           s% model_number, s% grada(k)
            end if
            if (s% model_number == -1) then
               write(*,4) 'grada', k, s% newton_iter, s% model_number, s% grada(k)
            end if
         end subroutine store_stuff
         
         subroutine write_info
            integer :: j
            include 'formats'
            write(*,*) 'do_eos_for_cell k', k
            write(*,*) 's% lnPgas_flag', s% lnPgas_flag
            write(*,1) 'z = ', z
            write(*,1) 'x = ', xh
            write(*,1) 'Xhe3 =', s% xa(he3,k)
            write(*,1) 'Xhe4 =', s% xa(he4,k)
            write(*,1) 'sum(xa) =', sum(s% xa(:,k))
            write(*,1) 'abar = ', s% abar(k)
            write(*,1) 'zbar = ', s% zbar(k)
            write(*,1) 'Pgas = ', s% Pgas(k)
            write(*,1) 'logPgas = ', s% lnPgas(k)/ln10
            write(*,1) 'rho = ', s% rho(k)
            write(*,1) 'logRho = ', s% lnd(k)/ln10
            write(*,1) 'T = ', s% T(k)
            write(*,1) 'logT = ', s% lnT(k)/ln10
            write(*,1) 'logQ = ', s% lnd(k)/ln10 - 2*s% lnT(k)/ln10 + 12
            if (s% E_flag) write(*,1) 'logV = ', s% lnd(k)/ln10 - 0.7d0*lnE/ln10 + 20
            if (s% E_flag) write(*,1) 'logE = ', lnE/ln10
            if (s% E_flag) write(*,1) 'lgT_old = ', s% lnT_old(k)/ln10
            if (s% lnPgas_flag) write(*,1) 'logW = ', s% lnPgas(k)/ln10 - 4*s% lnT(k)/ln10
            write(*,*)
            if (ierr /= 0) return
            write(*,1) 'P = ', s% P(k)
            write(*,1) 'lnP = ', s% lnP(k)
            write(*,1) 'Prad = ', s% Prad(k)
            write(*,1) 'logS = ', s% lnS(k)/ln10
            write(*,1) 'logE = ', s% lnE(k)/ln10
            write(*,1) 'grada = ', s% grada(k)
            write(*,1) 'cv = ', s% cv(k)
            write(*,1) 'cp = ', s% cp(k)
            write(*,1) 'gamma1 = ', s% gamma1(k)
            write(*,1) 'gamma3 = ', s% gamma3(k)
            write(*,1) 'eta = ', s% eta(k)
            write(*,1) 'gam = ', s% gam(k)
            write(*,1) 'mu = ', s% mu(k)
            write(*,1) 'log_free_e = ', s% lnfree_e(k)/ln10
            write(*,1) 'chiRho = ', s% chiRho(k)
            write(*,1) 'chiT = ', s% chiT(k)
            write(*,*)
            do j=1,s% species
               write(*,3) 'xa(j,k)', j, k, s% xa(j,k)
            end do            
         end subroutine write_info
         
      end subroutine do_eos_for_cell
   
      
      subroutine eos_get( &
            s, k, z, xh, abar, zbar, xa, &
            Rho, log10Rho, T, log10T, & 
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         use eos_lib, only: eosDT_get, eosDT_ideal_gas_get, eosDT_HELMEOS_get
         use eos_def, only: num_eos_basic_results, num_helm_results, i_Cp, i_Cv
         type (star_info), pointer :: s         
         integer, intent(in) :: k ! 0 means not being called for a particular cell
         real(dp), intent(in) :: &
            z, xh, abar, zbar, xa(:), Rho, log10Rho, T, log10T
         real(dp), dimension(num_eos_basic_results), intent(out) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer, intent(out) :: ierr
         
         real(dp) :: eos_z, eos_x
         real(dp) :: helm_res(num_helm_results)

         include 'formats'
         if (s% use_fixed_XZ_for_eos) then
            eos_x = s% fixed_X_for_eos
            eos_z = s% fixed_Z_for_eos
         else
            eos_x = xh
            eos_z = z
         end if
         if (s% use_other_eos) then
            call s% other_eosDT_get( &
               s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               Rho, log10Rho, T, log10T, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         else if (s% use_eosDT_ideal_gas) then
            call eosDT_ideal_gas_get( &
               s% eos_handle, eos_z, eos_x, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               Rho, log10Rho, T, log10T, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         else if (s% use_eosDT_HELMEOS) then
            call eosDT_HELMEOS_get( &
               s% eos_handle, eos_z, eos_x, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               Rho, log10Rho, T, log10T, &
               s% eosDT_HELMEOS_include_radiation, &
               s% eosDT_HELMEOS_always_skip_elec_pos, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, helm_res, ierr)
         else
            call eosDT_get( &
               s% eos_handle, eos_z, eos_x, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               Rho, log10Rho, T, log10T, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         end if
         if (ierr /= 0) then
            if (s% report_ierr) then
!$OMP critical (star_eos_get)
               write(*,*) 'eos_get ierr', ierr
               write(*,2) 'k', k
               write(*,1) 'z', z
               write(*,1) 'xh', xh
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'log10Rho', log10Rho
               write(*,1) 'log10T', log10T
!$OMP end critical (star_eos_get)

               stop 'eos_get'
               
            end if
            return
         end if
         if (.false. .and. (res(i_Cp) < 0 .or. res(i_Cv) < 0)) then
            write(*,*) 'eos_get ierr', ierr
            write(*,2) 'k', k
            write(*,1) 'res(i_Cp)', res(i_Cp)
            write(*,1) 'res(i_Cv)', res(i_Cv)
            write(*,1) 'z', z
            write(*,1) 'xh', xh
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'log10Rho', log10Rho
            write(*,1) 'log10T', log10T
            stop 'debug: eos_get'
         end if
         
      end subroutine eos_get
      
      
      subroutine check_eos_res(s, k, res, z, xh, ierr)
         use eos_def
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp) :: res(num_eos_basic_results)
         real(dp) :: z, xh
         integer, intent(out) :: ierr
         1 format(a40,1pd26.16)
         ierr = 0
         if (abs(s% lnPgas(k) - res(i_lnPgas)) > 1d-5*res(i_lnPgas)) then
            write(*,*) 'failure for check_eos_res for lnPgas k', k
            write(*,1) 'res(i_lnPgas)', res(i_lnPgas)
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% lnS(k) - res(i_lnS)) > 1d-5*res(i_lnS)) then
            write(*,*) 'failure for check_eos_res for entropy k', k
            write(*,1) 'res(i_lnS)', res(i_lnS)
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% lnE(k) - res(i_lnE)) > 1d-5*res(i_lnE)) then
            write(*,*) 'failure for check_eos_res for energy k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% grada(k) - res(i_grad_ad)) > 1d-5*res(i_grad_ad)) then
            write(*,*) 'failure for check_eos_res for grada k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% dE_drho(k) - res(i_dE_drho)) > 1d-5*abs(res(i_dE_drho))) then
            write(*,*) 'failure for check_eos_res for dE_dRho k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% cv(k) - res(i_Cv)) > 1d-5*res(i_Cv)) then
            write(*,*) 'failure for check_eos_res for cv k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% cp(k) - res(i_cp)) > 1d-5*res(i_cp)) then
            write(*,*) 'failure for check_eos_res for cp k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% gamma1(k) - res(i_gamma1)) > 1d-5*res(i_gamma1)) then
            write(*,*) 'failure for check_eos_res for gamma1 k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% gamma3(k) - res(i_gamma3)) > 1d-5*res(i_gamma3)) then
            write(*,*) 'failure for check_eos_res for gamma3 k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% chiRho(k) - res(i_chiRho)) > 1d-5*abs(res(i_chiRho))) then
            write(*,*) 'failure for check_eos_res for chiRho k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
         if (abs(s% chiT(k) - res(i_chiT)) > 1d-5*abs(res(i_chiT))) then
            write(*,*) 'failure for check_eos_res for chiT k', k
            call show_res_prev_diffs(s, k, res, z, xh)
            ierr = -1; return
         end if
      end subroutine check_eos_res
   
      
      subroutine eval_eosPT( &
            s, k, z, xh, abar, zbar, xa, &
            Pgas, log10Pgas, T, log10T, & 
            Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         use eos_lib, only: eosPT_get
         use eos_def, only: num_eos_basic_results, i_Cp, i_Cv
         type (star_info), pointer :: s         
         integer, intent(in) :: k ! 0 means not being called for a particular cell
         real(dp), intent(in) :: &
            z, xh, abar, zbar, xa(:), Pgas, log10Pgas, T, log10T
            
         real(dp), intent(out) :: &
            Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         real(dp), dimension(num_eos_basic_results), intent(out) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer, intent(out) :: ierr
         
         real(dp) :: eos_z, eos_x
         include 'formats'
         if (s% use_fixed_XZ_for_eos) then
            eos_x = s% fixed_X_for_eos
            eos_z = s% fixed_Z_for_eos
         else
            eos_x = xh
            eos_z = z
         end if
         ierr = 0
         if (s% use_other_eos) then
            call s% other_eosPT_get( &
               s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               Pgas, log10Pgas, T, log10T, &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         else
            call eosPT_get( &
               s% eos_handle, eos_z, eos_x, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               Pgas, log10Pgas, T, log10T, &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         end if
         if (ierr /= 0) then
            if (s% report_ierr) then
!$OMP critical (star_eosPT_get)
               write(*,*) 'eval_eosPT ierr', ierr
               write(*,2) 'k', k
               write(*,1) 'z', z
               write(*,1) 'xh', xh
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'log10Pgas', log10Pgas
               write(*,1) 'log10T', log10T
!$OMP end critical (star_eosPT_get)

               !stop 'eval_eosPT'
               
            end if
            return
         end if
         
      end subroutine eval_eosPT
      
      
      subroutine eval_eosDE( &
            s, k, z, xh, abar, zbar, xa, &
            rho, log10Rho, log10E, & 
            T, log10T, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
            ierr)
         use eos_def, only: num_eos_basic_results
         type (star_info), pointer :: s         
         integer, intent(in) :: k ! 0 means not being called for a particular cell
         real(dp), intent(in) :: &
            z, xh, abar, zbar, xa(:), rho, log10Rho, log10E
         real(dp), intent(out) :: T, log10T, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         real(dp), dimension(num_eos_basic_results), intent(out) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer, intent(out) :: ierr
         
         integer :: k_off
         real(dp) :: logT_guess
         real(dp), parameter :: logT_tol = 1d-9, logE_tol = 1d-9
         
         include 'formats'
         
         ierr = 0
         
         if (k /= 0) then
            if (s% generations > 1) then
               logT_guess = s% lnT_old(k)/ln10
            else
               logT_guess = s% lnT(k)/ln10
            end if
            if (is_bad_num(logT_guess)) then
               do k_off = 1, 10
                  if (k-k_off >= 1) then
                     if (s% generations > 1) then
                        logT_guess = s% lnT_old(k-k_off)/ln10
                     else
                        logT_guess = s% lnT(k-k_off)/ln10
                     end if
                     if (.not. is_bad_num(logT_guess)) exit
                  end if
                  if (k+k_off <= s% nz) then
                     if (s% generations > 1) then
                        logT_guess = s% lnT_old(k+k_off)/ln10
                     else
                        logT_guess = s% lnT(k+k_off)/ln10
                     end if
                     if (.not. is_bad_num(logT_guess)) exit
                  end if
               end do
               if (is_bad_num(logT_guess)) logT_guess = 7
            end if
         else
            logT_guess = 7
         end if

         call eval_eosDE_with_logT_guess( &
            s, k, z, xh, abar, zbar, xa, &
            rho, log10Rho, log10E, logT_guess, logT_tol, logE_tol, & 
            T, log10T, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
            ierr)

      end subroutine eval_eosDE     
      
      
      subroutine eval_eosDE_with_logT_guess( &
            s, k, z, xh, abar, zbar, xa, &
            rho, log10Rho, log10E, logT_guess, logT_tol, logE_tol, & 
            T, log10T, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
            ierr)
         ! note: dlnT_dlnE_const_Rho = 1/(dE_dT_const_rho*T/E) = E/(Cv*T)
         ! dlnT_dlnRho_const_E = -dlnE_dlnRho_const_T/dlnE_dnT_const_Rho
         !     = -(dE_dRho_const_T*Rho/E)/(dE_dT_const_Rho*T/E) = -Rho*dE_dRho/(Cv*T)
         use eos_def
         use eos_lib, only: eosDT_get_T, eosDE_get
         type (star_info), pointer :: s        
         integer, intent(in) :: k ! 0 indicates not for a particular cell. 
         real(dp), intent(in) :: &
            z, xh, abar, zbar, xa(:), rho, log10Rho, log10E, &
            logT_guess, logT_tol, logE_tol
         real(dp), intent(out) :: T, log10T, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         real(dp), dimension(num_eos_basic_results), intent(out) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer, intent(out) :: ierr

         integer :: which_other, max_iter, eos_calls, species
         real(dp) :: other_value, other_tol, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
            logT_result, sumx, eos_x, eos_z, energy, Cv
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         species = s% species
         
         if (s% use_fixed_XZ_for_eos) then
            eos_x = s% fixed_X_for_eos
            eos_z = s% fixed_Z_for_eos
         else
            eos_x = xh
            eos_z = z
         end if

         if (s% use_eosDE_get .and. s% job% eos_file_prefix == 'mesa') then
            
            energy = exp10_cr(log10E)
            call eosDE_get( &
               s% eos_handle, eos_z, eos_x, abar, zbar, &
               species, s% chem_id, s% net_iso, xa, &
               energy, log10E, rho, log10Rho, logT_guess, &
               T, log10T, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
               dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
               ierr)
            if (ierr /= 0 .and. s% report_ierr) then
!$OMP critical (star_eosDE_get)
               write(*,3) 'eosDE_get ierr', k, ierr
               write(*,1) 'eos_z', eos_z
               write(*,1) 'eos_x', eos_x
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'energy', energy
               write(*,1) 'log10E', log10E
               write(*,1) 'rho', rho
               write(*,1) 'log10Rho', log10Rho
               write(*,1) 'logT_guess', logT_guess
!$OMP end critical (star_eosDE_get)
               stop 'eval_eosDE_with_logT_guess'               
            end if
               
         else
         
            max_iter = 100
            which_other = i_lnE
            other_value = log10E*ln10
            other_tol = logE_tol*ln10
            logT_bnd1 = arg_not_provided
            other_at_bnd1 = arg_not_provided
            logT_bnd2 = arg_not_provided
            other_at_bnd2 = arg_not_provided
      
            if (s% use_other_eos) then
               call s% other_eosDT_get_T( &
                  s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
                  species, s% chem_id, s% net_iso, xa, &
                  log10Rho, which_other, other_value, &
                  logT_tol, other_tol, max_iter, logT_guess, & 
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  logT_result, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
                  eos_calls, ierr)               
            else
               call eosDT_get_T( &
                  s% eos_handle, eos_z, eos_x, abar, zbar, &
                  species, s% chem_id, s% net_iso, xa, &
                  log10Rho, which_other, other_value, &
                  logT_tol, other_tol, max_iter, logT_guess,  &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  logT_result, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
                  eos_calls, ierr)               
            end if

            log10T = logT_result
            T = exp10_cr(log10T)
            energy = exp10_cr(log10E)
            Cv = res(i_Cv)

            dlnT_dlnd_c_E = -rho*res(i_dE_dRho)/(T*Cv)
            dlnPgas_dlnd_c_E = &
               d_dlnd(i_lnPgas) + d_dlnT(i_lnPgas)*dlnT_dlnd_c_E

            dlnT_dlnE_c_Rho = energy/(T*Cv)
            dlnPgas_dlnE_c_Rho = d_dlnT(i_lnPgas)*dlnT_dlnE_c_Rho
            
         end if
         
         if (ierr /= 0) then
            if (s% report_ierr) then
!$OMP critical (star_eosDE_get)
               write(*,3) 'eval_eosDE ierr', k, ierr
               write(*,1) 'z', z
               write(*,1) 'xh', xh
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'log10Rho', log10Rho
               write(*,1) 'log10E', log10E
               write(*,1) 'logT_guess', logT_guess
!$OMP end critical (star_eosDE_get)
               stop 'eval_eosDE_with_logT_guess'               
            end if
         end if
         
         if (k == s% trace_k) then
            write(*,1) 'z', z
            write(*,1) 'xh', xh
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'rho', rho
            write(*,2) 's% lnd(k)', k, s% lnd(k)
            write(*,2) 's% xh(s% i_lnd,k)', k, s% xh(s% i_lnd,k)
            write(*,1) 'log10Rho', log10Rho
            write(*,1) 'log10E', log10E
            write(*,1) 'logT_guess', logT_guess
            write(*,1) 'log10T', log10T
         end if
         
      end subroutine eval_eosDE_with_logT_guess      
      
      
      real(dp) function eval_csound_and_rho(s, k, rho, ierr) result(cs)
         ! doesn't use any vars other than xh and xa
         use chem_lib
         use chem_def
         use eos_lib
         use eos_def
         
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(out) :: rho
         integer, intent(out) :: ierr
         
         integer :: species
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp) :: xsum, mass_correction, abar, zbar, z2bar, &
            ye, z, y, xh, lnd, T, lnT, log10T, &
            gamma1, log10Pgas, Pgas, Prad, P, lnE, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
            log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            dlnT_dlnE_const_Rho, dlnT_dlnRho_const_E

         ierr = 0
         cs = -1
         species = s% species
         call basic_composition_info( &
            species, s% chem_id, s% xa(1:species,k), &
            xh, y, abar, zbar, z2bar, ye, mass_correction, xsum)  
         z = 1 - (xh + y)
         
         if (s% lnPgas_flag) then
            log10Pgas = s% xh(s% i_lnPgas,k)/ln10
            Pgas = exp10_cr(log10Pgas)
            lnT = s% xh(s% i_lnT, k)
            T = exp_cr(lnT)
            call eval_eosPT( &
               s, k, z, xh, abar, zbar, s% xa(:,k), &
               Pgas, log10Pgas, T, lnT/ln10, & 
               rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,*) 'eval_csound: eval_eosPT ierr', ierr
               end if
               return
            end if
         else if (s% E_flag) then
            lnd = s% xh(s% i_lnd, k)
            rho = exp_cr(lnd)
            lnE = log_cr(s% xh(s% i_E, k))
            call eval_eosDE( &
               s, k, z, xh, abar, zbar, s% xa(:,k), &
               rho, lnd/ln10, lnE/ln10, & 
               T, log10T, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
               dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
               ierr)
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,*) 'eval_csound: eval_eosDE ierr', ierr
               end if
               return
            end if
            lnT = log10T*ln10
         else
            lnd = s% xh(s% i_lnd, k)
            rho = exp_cr(lnd)
            lnT = s% xh(s% i_lnT, k)
            T = exp_cr(lnT)
            call eos_get( &
               s, k, z, xh, abar, zbar, s% xa(:,k), &
               rho, lnd/ln10, T, lnT/ln10, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)         
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,*) 'eval_csound: eos_get ierr', ierr
               end if
               return
            end if
            Pgas = exp_cr(res(i_lnPgas))
         end if
         
         gamma1 = res(i_gamma1)
         Prad = crad * T*T*T*T / 3
         P = Pgas + Prad
         
         if (s% use_sr_sound_speed) then
            cs = sqrt(gamma1*P/(P + exp_cr(res(i_lnE)) + clight*clight*rho))
         else
            cs = sqrt(gamma1*P/rho)
         end if
         if (is_bad_num(cs)) then
            ierr = -1
            if (s% report_ierr) then
               write(*,*) 'eval_csound: bad num', cs
            end if
         end if

      end function eval_csound_and_rho
      
      
      subroutine eval_rho( &
            s, z, x, xa, abar, zbar, T, lnT, Pgas, &
            rho, res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar, ierr)
         use eos_lib, only: eosPT_get, eosDT_get_Rho, Radiation_Pressure
         use eos_def
         
         type (star_info), pointer :: s
         real(dp), intent(in) :: &
            z, x, xa(:), abar, zbar, T, lnT, Pgas
         real(dp), intent(out) :: rho
         real(dp), intent(out) :: &
            res(num_eos_basic_results), &
            d_eos_dlnd(num_eos_basic_results), &
            d_eos_dlnT(num_eos_basic_results), &
            d_eos_dabar(num_eos_basic_results), &
            d_eos_dzbar(num_eos_basic_results)
         integer, intent(out) :: ierr
         
         real(dp) :: logPgas, logT, &
            logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         integer :: which_other, max_iter, eos_calls
         real(dp) :: other, other_tol, logRho_tol, logRho_guess, &
            logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2
            
         include 'formats'
         
         ierr = 0

         logPgas = log10_cr(Pgas)
         logT = lnT/ln10
         
         if (s% use_other_eos) then
            call s% other_eosPT_get( &
               s% id, 0, s% eos_handle, z, x, abar, zbar,  &
               s% species, s% chem_id, s% net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar, &
               ierr)
         else
            call eosPT_get( &
               s% eos_handle, z, x, abar, zbar,  &
               s% species, s% chem_id, s% net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar, &
               ierr)
         end if
         which_other = i_lnPgas
         max_iter = 100
         other = logPgas*ln10
         logRho_tol = 1d-6
         other_tol = 1d-6
         logRho_guess = logRho ! use the eosPT result as guess
         logRho_bnd1 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         logRho_bnd2 = arg_not_provided
         other_at_bnd2 = arg_not_provided
         if (s% use_other_eos) then
            call s% other_eosDT_get_Rho( &
               s% id, 0, s% eos_handle, Z, X, abar, zbar, &
               s% species, s% chem_id, s% net_iso, xa, &
               logT, which_other, other, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho, res, d_eos_dlnd, d_eos_dlnT, &
               d_eos_dabar, d_eos_dzbar, eos_calls, ierr)
         else
            call eosDT_get_Rho( &
               s% eos_handle, Z, X, abar, zbar,  &
               s% species, s% chem_id, s% net_iso, xa, &
               logT, which_other, other, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho, res, d_eos_dlnd, d_eos_dlnT, &
               d_eos_dabar, d_eos_dzbar, eos_calls, ierr)
         end if
         Rho = exp10_cr(logRho)
         
      end subroutine eval_rho
      
      
      subroutine eval_lnPgas( &
            s, z, x, xa, abar, zbar, T, lnT, Rho, lnd, &
            dbg, lnPgas, ierr)
         use eos_lib, only: eosPT_get_Pgas_for_Rho
         use eos_def
         
         type (star_info), pointer :: s
         real(dp), intent(in) :: &
            z, x, xa(:), abar, zbar, T, lnT, Rho, lnd
         logical, intent(in) :: dbg
         real(dp), intent(out) :: lnPgas
         integer, intent(out) :: ierr
         
         real(dp) :: logPgas_guess, logPgas_result, logPgas_tol, logRho_tol, &
            logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2, rho_result, logRho_result, &
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_eos_dlnd, d_eos_dlnT, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_eos_dabar, d_eos_dzbar
            
         integer :: max_iter, eos_calls
            
         include 'formats'
         
         ierr = 0
         
         if (is_bad_num(z)) then
            ierr = -1
            return
            write(*,1) 'Z', Z
            write(*,1) 'X', X
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'lnT/ln10', lnT/ln10
            write(*,1) 'lnd/ln10', lnd/ln10
            stop 'eval_lnPgas'
         end if

         call eos_get( &
            s, 0, z, x, abar, zbar, xa, &
            Rho, lnd/ln10, T, lnT/ln10, &
            res, d_eos_dlnd, d_eos_dlnT, &
            d_eos_dabar, d_eos_dzbar, &
            ierr)
         if (ierr /= 0) return
         
         logPgas_guess = res(i_lnPgas)/ln10
         max_iter = 100
         logPgas_tol = 1d-6
         logRho_tol = 1d-6
         logPgas_bnd1 = arg_not_provided
         logRho_at_bnd1 = arg_not_provided
         logPgas_bnd2 = arg_not_provided
         logRho_at_bnd2 = arg_not_provided
         
         if (dbg) then
            write(*,*) 'call eosPT_get_Pgas_for_Rho'
            write(*,1) 'lnT/ln10', lnT/ln10
            write(*,1) 'lnd/ln10', lnd/ln10
            write(*,1) 'logPgas_guess', logPgas_guess
            write(*,1) 'Z', Z
            write(*,1) 'X', X
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'guess eta', res(i_eta)
         end if
         
         if (s% use_other_eos) then
            call s% other_eosPT_get_Pgas_for_Rho( &
               s% id, 0, s% eos_handle, Z, X, abar, zbar,  &
               s% species, s% chem_id, s% net_iso, xa, &
               lnT/ln10, lnd/ln10, &
               logPgas_tol, logRho_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2, &
               logPgas_result, rho_result, logRho_result, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, & 
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_eos_dabar, d_eos_dzbar, eos_calls, ierr)
         else
            call eosPT_get_Pgas_for_Rho( &
               s% eos_handle, Z, X, abar, zbar,  &
               s% species, s% chem_id, s% net_iso, xa, &
               lnT/ln10, lnd/ln10, &
               logPgas_tol, logRho_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2, &
               logPgas_result, rho_result, logRho_result, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, & 
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_eos_dabar, d_eos_dzbar, eos_calls, ierr)
         end if
         lnPgas = logPgas_result*ln10
         
         if (dbg) then
            write(*,1) 'done eosPT_get_Pgas_for_Rho', lnPgas/ln10
            write(*,2) 'ierr', ierr
            stop 'eval_lnPgas'
         end if
         
      end subroutine eval_lnPgas
      
      
      subroutine eval_eos_rho_e(s,k,ierr)
         ! use lnd(k) and lnE(k) to set lnT(k)
         use chem_lib, only: basic_composition_info
         use eos_def
         use eos_lib, only: eosDT_get_T
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         real(dp) :: sumx, z, T, log10T, &
            res(num_eos_basic_results), &
            d_dlnd(num_eos_basic_results), &
            d__dlnT(num_eos_basic_results), &
            d_dabar(num_eos_basic_results), &
            d_dzbar(num_eos_basic_results)
         integer :: species
         
         include 'formats'
         
         if (s% use_eosDT_ideal_gas) then
            call eval_eosDT_rho_e(s,k,ierr)
            return
         end if
         
         ierr = 0
         species = s% species
         call basic_composition_info( &
            species, s% chem_id, s% xa(1:species,k), s% X(k), s% Y(k), s% abar(k), &
            s% zbar(k), s% z2bar(k), s% ye(k), s% mass_correction(k), sumx)
         
         z = max(0d0,1d0-(s% X(k)+s% Y(k)))

         call eval_eosDE( &
               s, k, z, s% X(k), s% abar(k), s% zbar(k), s% xa(:,k), &
               s% rho(k), s% lnd(k)/ln10, s% lnE(k)/ln10, & 
               T, log10T, res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
               s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), &
               s% dlnT_dlnE_c_Rho(k), s% dlnT_dlnd_c_E(k), &
               s% dlnPgas_dlnE_c_Rho(k), s% dlnPgas_dlnd_c_E(k), &
               ierr)
         if (ierr /= 0) return

         s% lnT(k) = log10T*ln10
         s% xh(s% i_lnT,k) = s% lnT(k)
         s% T(k) = T
         s% lnPgas(k) = res(i_lnPgas)
         s% Pgas(k) = exp_cr(s% lnPgas(k))
         if (s% E_flag) then
            call set_dlnP_dlnd_c_E(s,k)
            call set_dlnP_dlnE_c_Rho(s,k)
         end if

      end subroutine eval_eos_rho_e

      
      subroutine set_dlnP_dlnd_c_E(s,k)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp) :: Prad, dPrad_dlnd_c_E, Pgas, dPgas_dlnd_c_E, dP_dlnd_c_E         
         include 'formats'         
         if (.not. s% E_flag) then
            write(*,2) 'eval_dlnP_dlnd_c_E requires E_flag', k
            stop 'eval_dlnP_dlnd_c_E'
         end if        
         Prad = s% Prad(k) ! = crad * T*T*T*T / 3
         dPrad_dlnd_c_E = 4*Prad*s% dlnT_dlnd_c_E(k)        
         Pgas = s% Pgas(k)
         dPgas_dlnd_c_E = Pgas*s% dlnPgas_dlnd_c_E(k)       
         dP_dlnd_c_E = dPgas_dlnd_c_E + dPrad_dlnd_c_E
         s% dlnP_dlnd_c_E(k) = dP_dlnd_c_E/s% P(k)      
      end subroutine set_dlnP_dlnd_c_E
      
      
      subroutine set_dlnP_dlnE_c_Rho(s,k)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp) :: Prad, dPrad_dlnE_c_Rho, Pgas, dPgas_dlnE_c_Rho, dP_dlnE_c_Rho         
         include 'formats'         
         if (.not. s% E_flag) then
            write(*,2) 'eval_dlnP_dlnE_c_Rho requires E_flag', k
            stop 'eval_dlnP_dlnE_c_Rho'
         end if        
         Prad = s% Prad(k) ! = crad * T*T*T*T / 3
         dPrad_dlnE_c_Rho = 4*Prad*s% dlnT_dlnE_c_Rho(k)        
         Pgas = s% Pgas(k)
         dPgas_dlnE_c_Rho = Pgas*s% dlnPgas_dlnE_c_Rho(k)       
         dP_dlnE_c_Rho = dPgas_dlnE_c_Rho + dPrad_dlnE_c_Rho
         s% dlnP_dlnE_c_Rho(k) = dP_dlnE_c_Rho/s% P(k)  
         if (.false.) then ! k == 878 .or. k == 879 .or. k == 880) then
!$OMP critical
            write(*,2) 'T', k, s% T(k)
            write(*,2) 'Prad', k, Prad
            write(*,2) 'dPrad_dlnE_c_Rho', k, dPrad_dlnE_c_Rho
            write(*,2) 'Pgas', k, Pgas
            write(*,2) 'dPgas_dlnE_c_Rho', k, dPgas_dlnE_c_Rho
            write(*,2) 'dP_dlnE_c_Rho', k, dP_dlnE_c_Rho
            write(*,2) 'P', k, s% P(k)
            write(*,2) 's% dlnP_dlnE_c_Rho(k)', k, s% dlnP_dlnE_c_Rho(k)
            write(*,2) 'long form', k, &
                s% chiT(k)*(s% energy(k)/(s% Cv(k)*s% T(k)))
            write(*,2) 's% chiT(k)', k, s% chiT(k)
            write(*,2) 's% energy(k)', k, s% energy(k)
            write(*,2) 's% Cv(k)', k, s% Cv(k)
            write(*,2) 's% T(k)', k, s% T(k)
            write(*,2) 'logP', k, s% lnP(k)/ln10
            write(*,2) 'logT', k, s% lnT(k)/ln10
            write(*,2) 'logE', k, s% lnE(k)/ln10
            write(*,2) 'logRho', k, s% lnd(k)/ln10
!$OMP end critical
         end if   
      end subroutine set_dlnP_dlnE_c_Rho
      
      
      subroutine eval_eosDT_rho_e(s,k,ierr)
         ! use lnd(k) and lnE(k) to set lnT(k)
         use chem_lib, only: basic_composition_info
         use eos_def
         use eos_lib, only: eosDT_get_T
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         real(dp) :: &
            res(num_eos_basic_results), &
            d_eos_dlnd(num_eos_basic_results), &
            d_eos_dlnT(num_eos_basic_results), &
            d_eos_dabar(num_eos_basic_results), &
            d_eos_dzbar(num_eos_basic_results)
         integer :: which_other, max_iter, eos_calls, species
         real(dp) :: other_value, other_tol, logT_tol, logT_guess, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
            z, logRho, logT_result, sumx
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         species = s% species
         call basic_composition_info( &
            species, s% chem_id, s% xa(1:species,k), s% X(k), s% Y(k), &
            s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), s% mass_correction(k), sumx)
         
         z = max(0d0,1d0-(s% X(k)+s% Y(k)))
         
         max_iter = 100
         logRho = s% lnd(k)/ln10
         which_other = i_lnE
         other_value = s% lnE(k)
         if (s% generations > 1) then
            logT_guess = s% lnT_old(k)/ln10
         else
            logT_guess = s% lnT(k)/ln10
         end if
         logT_tol = 1d-8
         other_tol = 1d-8
         logT_bnd1 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         logT_bnd2 = arg_not_provided
         other_at_bnd2 = arg_not_provided
      
         if (s% use_other_eos) then
            call s% other_eosDT_get_T( &
               s% id, k, s% eos_handle, z, s% X(k), s% abar(k), s% zbar(k), &
               species, s% chem_id, s% net_iso, s% xa(:,k), &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, & 
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
               s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), eos_calls, ierr)               
         else
            call eosDT_get_T( &
               s% eos_handle, z, s% X(k), s% abar(k), s% zbar(k), &
               species, s% chem_id, s% net_iso, s% xa(:,k), &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess,  &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k),  &
               s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), eos_calls, ierr)               
         end if
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,2) 'eval_eos_rho_e: eosDT_get_T failed', k
            if (dbg) then
               write(*,2) 'z', k, z
               write(*,2) 's% X(k)', k, s% X(k)
               write(*,2) 's% abar(k)', k, s% abar(k)
               write(*,2) 's% zbar(k)', k, s% zbar(k)
               write(*,2) 'logRho', k, logRho
               write(*,2) 'logT_guess', k, logT_guess
               write(*,2) 'other_value', k, other_value
               stop 'eval_eos_rho_e'
            end if
            return
         end if

         s% lnT(k) = logT_result*ln10
         s% xh(s% i_lnT,k) = s% lnT(k)
         s% T(k) = exp_cr(s% lnT(k))
         s% lnPgas(k) = res(i_lnPgas)
         s% Pgas(k) = exp_cr(s% lnPgas(k))

         call store_eos_for_cell(s, k, res, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,2) 'eval_eos_rho_e: store_eos_for_cell failed', k
            return
         end if

      end subroutine eval_eosDT_rho_e


      subroutine show_res_prev_diffs(s, k, res, z, xh)
         use eos_def
         use chem_def
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp) :: res(num_eos_basic_results)
         real(dp) :: z, xh
         integer :: he3, he4
         he3 = s% net_iso(ihe3)
         he4 = s% net_iso(ihe4)
         
         1 format(a40,1pd26.16)
         write(*,*)
         write(*,*) 'show_res_prev_diffs: k', k
         write(*,*)
         write(*,1) 'z = ', z
         write(*,1) 'x = ', xh
         write(*,1) 'Xhe3 =', s% xa(he3,k)
         write(*,1) 'Xhe4 =', s% xa(he4,k)
         write(*,1) 'sum(xa) =', sum(s% xa(:,k))
         write(*,1) 'abar = ', s% abar(k)
         write(*,1) 'zbar = ', s% zbar(k)
         write(*,1) 'rho = ', s% rho(k)
         write(*,1) 'log10_rho = ', s% lnd(k)/ln10
         write(*,1) 'T = ', s% T(k)
         write(*,1) 'log10_T = ', s% lnT(k)/ln10
         write(*,*)
         write(*,1) 's% lnPgas(k)', s% lnPgas(k)
         write(*,1) 'res(i_lnPgas)',res(i_lnPgas)
         write(*,1) 's% lnPgas(k)-res(i_lnPgas)',s% lnPgas(k)-res(i_lnPgas)
         write(*,*) 
         write(*,1) 's% lnS(k)', s% lnS(k)
         write(*,1) 'res(i_lnS)', res(i_lnS)
         write(*,1) 's% lnS(k)-res(i_lnS)', s% lnS(k)-res(i_lnS)
         write(*,*)
         write(*,1) 's% lnE(k)', s% lnE(k)
         write(*,1) 'res(i_lnE)', res(i_lnE)
         write(*,1) 's% lnE(k)-res(i_lnE)', s% lnE(k)-res(i_lnE)
         write(*,*)
         write(*,1) 's% grada(k)', s% grada(k)
         write(*,1) 'res(i_grad_ad)' , res(i_grad_ad)
         write(*,1) 's% grada(k)-res(i_grad_ad)' , s% grada(k)-res(i_grad_ad)
         write(*,*)

         write(*,1) 's% dE_drho(k)', s% dE_drho(k)
         write(*,1) 'res(i_dE_drho)' , res(i_dE_drho)
         write(*,1) 's% dE_drho(k)-res(i_dE_drho)' , s% dE_drho(k)-res(i_dE_drho)
         write(*,*)

         write(*,1) 's% cv(k)', s% cv(k)
         write(*,1) 'res(i_Cv)' , res(i_Cv)
         write(*,1) 's% cv(k)-res(i_Cv)' , s% cv(k)-res(i_Cv)
         write(*,*)

         write(*,1) 's% cp(k)', s% cp(k)
         write(*,1) 'res(i_cp)' , res(i_cp)
         write(*,1) 's% cp(k)-res(i_cp)' , s% cp(k)-res(i_cp)
         write(*,*)

         write(*,1) 's% gamma1(k)', s% gamma1(k)
         write(*,1) 'res(i_gamma1)' , res(i_gamma1)
         write(*,1) 's% gamma1(k)-res(i_gamma1)' , s% gamma1(k)-res(i_gamma1)
         write(*,*)

         write(*,1) 's% gamma3(k)', s% gamma3(k)
         write(*,1) 'res(i_gamma3)' , res(i_gamma3)
         write(*,1) 's% gamma3(k)-res(i_gamma3)' , s% gamma3(k)-res(i_gamma3)
         write(*,*)

         write(*,1) 's% eta(k)', s% eta(k)
         write(*,1) 'res(i_eta)' , res(i_eta)
         write(*,1) 's% eta(k)-res(i_eta)' , s% eta(k)-res(i_eta)
         write(*,*)

         write(*,1) 's% mu(k)', s% mu(k)
         write(*,1) 'res(i_mu)' , res(i_mu)
         write(*,1) 's% mu(k)-res(i_mu)' , s% mu(k)-res(i_mu)
         write(*,*)

         write(*,1) 's% lnfree_e(k)', s% lnfree_e(k)
         write(*,1) 'res(i_lnfree_e)' , res(i_lnfree_e)
         write(*,1) 's% lnfree_e(k)-res(i_lnfree_e)' , s% lnfree_e(k)-res(i_lnfree_e)
         write(*,*)

         write(*,1) 's% chiRho(k)', s% chiRho(k)
         write(*,1) 'res(i_chiRho)' , res(i_chiRho)
         write(*,1) 's% chiRho(k)-res(i_chiRho)' , s% chiRho(k)-res(i_chiRho)
         write(*,*)

         write(*,1) 's% chiT(k)', s% chiT(k)
         write(*,1) 'res(i_chiT)' , res(i_chiT)
         write(*,1) 's% chiT(k)-res(i_chiT)' , s% chiT(k)-res(i_chiT)
         write(*,*)

      end subroutine show_res_prev_diffs
      
      
      subroutine store_eos_for_cell(s, k, res, ierr)
         use eos_def
         use eos_lib
         use star_utils, only: eval_csound
         type (star_info), pointer :: s
         integer, intent(in) :: k      
         real(dp), intent(in) ::res(num_eos_basic_results)
         integer, intent(out) :: ierr
         
         real(dp) :: d_theta_e_deta, T
         
         include 'formats'
         ierr = 0
         T = s% T(k)
         s% Prad(k) = crad * T*T*T*T / 3
         s% P(k) = s% Prad(k) + s% Pgas(k)
         s% lnP(k) = log_cr(s% P(k))
         s% lnS(k) = res(i_lnS)
         s% lnE(k) = res(i_lnE)
         s% energy(k) = exp_cr(s% lnE(k))
         s% grada(k) = res(i_grad_ad)
         s% dE_drho(k) = res(i_dE_drho)
         s% cv(k) = res(i_Cv)
         s% cp(k) = res(i_cp)
         s% gamma1(k) = res(i_gamma1)
         s% gamma3(k) = res(i_gamma3)
         s% eta(k) = res(i_eta)
         s% gam(k) = pow2(s% zbar(k)*qe) * &
            pow_cr((4.0d0/3.0d0)*pi*avo*s% rho(k)/s% abar(k),one_third) / (kerg*T)
         s% mu(k) = res(i_mu)
         s% lnfree_e(k) = res(i_lnfree_e)
         s% chiRho(k) = res(i_chiRho)
         s% chiT(k) = res(i_chiT)
         if (s% screening_mode == 'classic') then
            s% theta_e(k) = eos_theta_e(res(i_eta), d_theta_e_deta)
         else
            s% theta_e(k) = 0d0
         end if
         if (s% T_start(k) < 0) s% T_start(k) = s% T(k)
         s% csound(k) = eval_csound(s,k,ierr)
         if (s% csound_start(k) < 0) s% csound_start(k) = s% csound(k)
         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,2) 'csound', k, s% csound(k)
               write(*,2) 's% gamma1(k)', k, s% gamma1(k)
               write(*,2) 's% P(k)', k, s% P(k)
               write(*,2) 's% rho(k)', k, s% rho(k)
            end if
            ierr = -1
         end if
         
      end subroutine store_eos_for_cell
      
      
      subroutine shutdown_microphys
         use eos_lib
         use kap_lib
         use net_lib
         call eos_shutdown
         call kap_shutdown
         call net_shutdown
      end subroutine shutdown_microphys
      
      
      subroutine adjust1_for_new_xa(s,k,ierr)
         ! adjust T for change in composition at constant P and Rho
         ! i.e., restore P(k) = P_prev
         ! keeping lnd(k) fixed and changing lnT(k) to compensate
         ! for changes in xa(:,k)
         use eos_def
         use chem_lib, only: basic_composition_info

         type (star_info), pointer :: s         
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         
         integer :: species, j
         real(dp) :: xsum, P_prev, fact
         real(dp), dimension(num_eos_basic_results) :: res
         
         include 'formats'
         
         ierr = 0
         species = s% species

         call basic_composition_info( &
            species, s% chem_id, s% xa(1:species,k), s% X(k), s% Y(k), &
            s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), &
            s% mass_correction(k), xsum)
      
         P_prev = s% P(k)
         if (k == -1) write(*,3) 'old T P', 0, k, s% T(k), P_prev
         do j = 1, 20 ! max tries for fixing T
         
            call eos_get( &
               s, k, 1d0-max(0d0,min(1d0,s% X(k) + s% Y(k))), &
               s% X(k), s% abar(k), s% zbar(k), s% xa(:,k), &
               s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
               res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
               s% d_eos_dabar(:,k), s% d_eos_dzbar(:,k), ierr)
            if (ierr /= 0) then
               write(*,3) 'check mix P failed in eos_get', j, k, s% T(k)
               exit
            end if
            s% lnPgas(k) = res(i_lnPgas)
            s% Pgas(k) = exp_cr(s% lnPgas(k))
            call store_eos_for_cell(s, k, res, ierr)
            if (ierr /= 0) then
               write(*,3) 'adjust1_for_new_xa failed store_eos_for_cell', &
                  j, k, s% T(k)
               exit
            end if
            fact = s% P(k) - P_prev
            if (j > 10) fact = fact/dble(j*j) ! undercorrect and relax tolerance
            s% T(k) = s% T(k) - fact/(s% chiT(k)*s% P(k)/s% T(k))
            s% lnT(k) = log_cr(s% T(k))
            s% xh(s% i_lnT,k) = s% lnT(k)

            if (abs(fact) < 1d-7*P_prev) then
               if (k == -1) write(*,3) 'OKAY: fix P after mix', j, k
               exit
            end if
            if (k == -1) write(*,3) 'new T P fact', &
               j, k, s% T(k), s% P(k), fact/P_prev
            
         end do


      end subroutine adjust1_for_new_xa
      


      end module micro

