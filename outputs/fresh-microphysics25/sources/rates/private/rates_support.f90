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

      module rates_support
      use const_def, only: dp
      use crlibm_lib
      use rates_def
      
      implicit none

      integer, parameter :: cache_version = 4
      
      
      contains
      
#ifdef offload
      !dir$ attributes offload: mic :: do_get_raw_rates
#endif         
      subroutine do_get_raw_rates( &
            num_reactions, reaction_id, rattab, rattab_f1, nT8s, &
            ye, logtemp, btemp, bden, raw_rate_factor, logttab, &
            rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
         use const_def, only : missing_value
         integer, intent(in) :: num_reactions, reaction_id(:), nT8s
         real(dp), intent(in) ::  &
            ye, logtemp, btemp, bden, raw_rate_factor(:),  &
            rattab(:,:), logttab(:)
         real(dp), pointer, intent(in) :: rattab_f1(:)
         real(dp), intent(out), dimension(:) :: rate_raw, rate_raw_dT, rate_raw_dRho
         integer, intent(out) :: ierr
         
         integer :: imax, iat0, iat, ir, i, j, irho
         integer, parameter :: mp = 4
         real(dp) :: dtab(num_reactions), ddtab(num_reactions), fac
         real(dp), pointer :: rattab_f(:,:,:)
         
         include 'formats.dek'
         
         ierr = 0
         
         rattab_f(1:4,1:nT8s,1:num_reactions) => rattab_f1(1:4*nT8s*num_reactions)

         do i = 1,num_reactions
         
            ir = reaction_id(i)
            !dtab(i) = ye**reaction_ye_rho_exponents(1,ir)
            select case(reaction_ye_rho_exponents(1,ir))
            case (0)
               dtab(i) = 1d0
            case (1)
               dtab(i) = ye
            case (2)
               dtab(i) = ye*ye
            end select

            !dtab(i) = dtab(i)*bden**reaction_ye_rho_exponents(2,ir)
            irho = reaction_ye_rho_exponents(2,ir)
            select case(irho)
            case (1)
               dtab(i) = dtab(i)*bden
            case (2)
               dtab(i) = dtab(i)*bden*bden
            case (3)
               dtab(i) = dtab(i)*bden*bden*bden
            case (4)
               dtab(i) = dtab(i)*bden*bden*bden*bden
            end select

            ddtab(i) = irho*dtab(i)/bden
            
         end do

         if (nrattab > 1) then
            imax = nrattab
            if (logtemp > rattab_thi) then
               ierr = -1
               return
            end if
            iat0 = int((logtemp - rattab_tlo)/rattab_tstp) + 1
            iat = max(1, min(iat0 - mp/2 + 1, imax - mp + 1))
            call get_rates_from_table(1, num_reactions)
         else ! table only has a single temperature
            do i=1,num_reactions
               rate_raw(i) = rattab(i,1)*dtab(i)
               rate_raw_dT(i) = 0
               rate_raw_dRho(i) = rate_raw(i)*ddtab(i)/dtab(i)
            end do
         end if
         
         do i=1,num_reactions
            fac = raw_rate_factor(i)
            rate_raw(i) = rate_raw(i)*fac
            rate_raw_dT(i) = rate_raw_dT(i)*fac
            rate_raw_dRho(i) = rate_raw_dRho(i)*fac
         end do

         
         contains
         
#ifdef offload
      !dir$ attributes offload: mic :: get_rates_from_table
#endif         
         subroutine get_rates_from_table(r1, r2)
            use const_def, only: ln10
            integer, intent(in) :: r1, r2
            
            integer :: i, k, cnt
            real(dp) :: denom, am1, a00, ap1, ap2, cm1, c00, cp1, cp2,  &
                  rate, dr_dT, dx, dt, old_rate, old_dr_dT
            
            include 'formats.dek'
               
            k = iat+1 ! starting guess for search
            do while (logtemp < logttab(k) .and. k > 1)
               k = k-1
            end do
            do while (logtemp > logttab(k+1) .and. k+1 < nrattab)
               k = k+1
            end do
            dt = logtemp - logttab(k)
            
            do i = r1,r2
            
               rate_raw(i) =  &
                     (rattab_f(1,k,i) + dt*(rattab_f(2,k,i) +   &
                           dt*(rattab_f(3,k,i) + dt*rattab_f(4,k,i))) &
                              ) * dtab(i)

               rate_raw_dRho(i) = rate_raw(i) * ddtab(i) / dtab(i)

               rate_raw_dT(i) =  &
                     (rattab_f(2,k,i) + 2*dt*(rattab_f(3,k,i) +   &
                           1.5d0*dt*rattab_f(4,k,i)) &
                              ) * dtab(i) / (btemp * ln10)

            end do
            
         end subroutine get_rates_from_table
         
      
      end subroutine do_get_raw_rates

      
      subroutine do_make_rate_tables( &
           num_reactions, cache_suffix, net_reaction_id, which_rates,  &
           rattab, rattab_f1, nT8s, ttab, logttab, ierr)  
         use const_def
         use interp_1d_lib, only: interp_pm, interp_m3q
         use interp_1d_def, only: pm_work_size, mp_work_size
         use utils_lib
         integer, intent(in) :: nT8s, num_reactions, net_reaction_id(:), which_rates(:)
         character (len=*), intent(in) :: cache_suffix
         real(dp) :: rattab(:,:), ttab(:), logttab(:)
         real(dp), pointer :: rattab_f1(:)
         integer, intent(out) :: ierr
         
         integer :: i, j, operr, ir, io_unit, num_to_add_to_cache
         real(dp) :: logT, btemp
         real(dp), pointer :: work(:), work1(:), f1(:), rattab_f(:,:,:)
         integer, pointer :: reaction_id(:)

         logical :: all_okay, a_okay, all_in_cache
         
         include 'formats.dek'
         
         ierr = 0
         io_unit = alloc_iounit(ierr)
         if (ierr /= 0) return
         
         rattab_f(1:4,1:nrattab,1:num_reactions) =>  &
               rattab_f1(1:4*nrattab*num_reactions)
         
         allocate(reaction_id(num_reactions))
         reaction_id(:) = net_reaction_id(:)
         
         num_to_add_to_cache = 0
         if (nrattab == 1) then
            all_in_cache = .false.
         else
            all_in_cache = .true.
            do i=1, num_reactions
               if (read_reaction_from_cache( &
                        net_reaction_id, cache_suffix, which_rates, i, io_unit, rattab)) then
                  reaction_id(i) = 0
                  cycle
               end if
               all_in_cache = .false.
               num_to_add_to_cache = num_to_add_to_cache + 1
               write(*,'(a)') 'create rate data for ' //  &
                     trim(reaction_Name(net_reaction_id(i)))
               !stop
            end do
         end if
         
         if (all_in_cache) then
         
!$OMP PARALLEL DO PRIVATE(i, logT, btemp)
            do i=1, nrattab
               logT = rattab_tlo + float(i-1)*rattab_tstp
               btemp = exp10_cr(logT)
               ttab(i) = btemp
               logttab(i) = logT
            end do
!$OMP END PARALLEL DO

         else 
            
            if (num_to_add_to_cache > 20) then
               write(*,2) 'number being added to net cache:', num_to_add_to_cache
               if (num_to_add_to_cache > 100) write(*,*) 'this will take some time .....'
            end if
            all_okay = .true.
!$OMP PARALLEL DO PRIVATE(i, operr, logT, btemp, a_okay, j)
            do i=1, nrattab
               logT = rattab_tlo + float(i-1)*rattab_tstp
               btemp = exp10_cr(logT)
               ttab(i) = btemp
               logttab(i) = logT
               do j=1, num_reactions
                  if (reaction_id(j) <= 0) cycle
                  rattab(j, i) = missing_value ! so can check
               end do
               operr = 0
               !write(*,2) 'logT', i, logT
               call get_net_rates_for_tables(which_rates, reaction_id, logT, btemp, num_reactions,  &
                     rattab(1:num_reactions, i), operr)
               if (operr /= 0) then
                  ierr = -1
                  cycle
               end if
               a_okay = .true.
               do j = 1, num_reactions
                  if (reaction_id(j) <= 0) cycle
                  if (rattab(j, i) == missing_value) then
                     write(*, '(a,i4,2x,a)') 'missing raw rate for ',  &
                        j, trim(reaction_Name(reaction_id(j)))
                     a_okay = .false.
                  end if
               end do
               if (.not. a_okay) all_okay = .false.
            end do
!$OMP END PARALLEL DO
            if (.not. all_okay) stop 'make_rate_tables'
            if (ierr /= 0) then
               write(*,*) 'make_rate_tables failed'
               return
            end if
         end if

         if (nrattab > 1) then ! create interpolants
            allocate(work(nrattab*mp_work_size*num_reactions), stat=ierr)
            if (ierr /= 0) return
!$OMP PARALLEL DO PRIVATE(i,operr,work1,f1)
            do i=1,num_reactions
               work1(1:nrattab*mp_work_size) =>  &
                     work(nrattab*mp_work_size*(i-1)+1:nrattab*mp_work_size*i)
               rattab_f(1,1:nrattab,i) = rattab(i,1:nrattab)
               f1(1:4*nT8s) => rattab_f1(1+4*nT8s*(i-1):4*nT8s*i)
               call interp_m3q(logttab, nrattab, f1, mp_work_size, work1,  &
                        'rates do_make_rate_tables', operr)
               if (operr /= 0) ierr = -1
            end do
!$OMP END PARALLEL DO
            deallocate(work)
         end if
         
         if (ierr == 0 .and. nrattab > 1 .and. .not. all_in_cache) then
            do i=1, num_reactions
               if (reaction_id(i) <= 0) cycle
               call write_reaction_to_cache(reaction_id, cache_suffix, which_rates, i, io_unit, rattab) 
            end do
         end if
         
         call free_iounit(io_unit)
         deallocate(reaction_id)

      end subroutine do_make_rate_tables
      
      
      subroutine reaction_filename(ir, cache_suffix, which, cache_filename, ierr)
         integer, intent(in) :: ir, which
         character (len=*), intent(in) :: cache_suffix
         character (len=*), intent(out) :: cache_filename
         integer, intent(out) :: ierr
         character (len=64) :: suffix
         ierr = 0
         if (which == 0 .and. len_trim(cache_suffix) > 0) then
            suffix = cache_suffix
         else
            if (which < 0) then
               ierr = -1
               suffix = '?'
            else if (which >= 100) then
               write(suffix,'(i3)') which
            else if (which >= 10) then
               write(suffix,'(i2)') which
            else
               write(suffix,'(i1)') which
            end if
         end if
         write(cache_filename,'(a)')  &
               trim(rates_cache_dir) // '/' //  &
               trim(reaction_Name(ir)) // '_' // trim(suffix) // '.bin'
      end subroutine reaction_filename
      

      
      logical function read_reaction_from_cache(reaction_id, cache_suffix, which_rates, i, io_unit, rattab) 
         integer, intent(in) :: i, io_unit, reaction_id(:), which_rates(:)
         character (len=*), intent(in) :: cache_suffix
         real(dp),intent(out) :: rattab(:,:)
         
         integer :: file_version, file_nrattab, file_which
         real(dp) :: file_rattab_thi, file_rattab_tlo, file_rattab_tstp
         character (len=256) :: cache_filename
         integer :: ios, ir, which, j, ierr
         real(dp), parameter :: tiny = 1e-6
         character (len=maxlen_reaction_Name) :: name
         real(dp), pointer :: rates(:)
         
         logical, parameter :: show_read_cache = .false.

         read_reaction_from_cache = .false.
         ierr = 0
         
         ir = reaction_id(i)
         which = which_rates(ir)
         if (raw_rates_records(ir)% use_rate_table) which = 0
         
         call reaction_filename(reaction_id(i), cache_suffix, which, cache_filename, ierr)
         if (ierr /= 0) then
            if (show_read_cache) write(*,*) 'read cache -- bad reaction_filename ' // trim(cache_filename)
            return
         end if   
         
         ios = 0
         open(unit=io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
         if (ios /= 0) then
            if (show_read_cache) write(*,*) 'read cache failed for open ' // trim(cache_filename)
            return
         end if       
         
         read(io_unit, iostat=ios)  &
               name, file_which, file_version, file_nrattab,  &
               file_rattab_thi, file_rattab_tlo, file_rattab_tstp
         if (ios /= 0) then
            if (show_read_cache) write(*,*) 'read cache failed for read header ' // trim(cache_filename)
            close(io_unit)
            return
         end if       
         
         if (name /= reaction_Name(ir)) then
            if (show_read_cache) write(*,*) 'read cache failed for name'
            close(io_unit)
            return
         end if
         
         if (which /= file_which) then
            if (show_read_cache) write(*,*) 'read cache failed for which reaction'
            close(io_unit)
            return
         end if
         
         if (cache_version /= file_version) then
            if (show_read_cache) write(*,*) 'read cache failed for version'
            close(io_unit)
            return
         end if

         if (abs(rattab_thi-file_rattab_thi) > tiny) then
            if (show_read_cache) write(*,*) 'read cache failed for rattab_thi'
            close(io_unit)
            return
         end if

         if (abs(rattab_tlo-file_rattab_tlo) > tiny) then
            if (show_read_cache) write(*,*) 'read cache failed for rattab_tlo'
            close(io_unit)
            return
         end if

         if (abs(rattab_tstp-file_rattab_tstp) > tiny) then
            if (show_read_cache) write(*,*) 'read cache failed for rattab_tstp'
            close(io_unit)
            return
         end if
         
         do j = 1, nrattab
            read(io_unit, iostat=ios) rattab(i,j)
            if (ios /= 0) then
               if (show_read_cache) write(*,*) 'read cache failed for reaction'
               close(io_unit)
               return
            end if
         end do         
         
         close(io_unit)
         
         read_reaction_from_cache = .true.
      
      end function read_reaction_from_cache
      
      
      
      subroutine write_reaction_to_cache(reaction_id, cache_suffix, which_rates, i, io_unit, rattab) 
         integer, intent(in) :: i, io_unit
         character (len=*), intent(in) :: cache_suffix
         integer, intent(in) :: reaction_id(:), which_rates(:)
         real(dp), intent(in) :: rattab(:,:)

         character (len=256) :: cache_filename
         integer :: ios, ir, which, ierr, j
         
         logical, parameter :: show_write_cache = .true.
         
         ierr = 0
         ir = reaction_id(i)
         which = which_rates(ir)
         if (raw_rates_records(ir)% use_rate_table) which = 0
         
         call reaction_filename(reaction_id(i), cache_suffix, which, cache_filename, ierr)
         if (ierr /= 0) return
         
         ios = 0
         open(unit=io_unit, file=trim(cache_filename), iostat=ios,  &
               action='write', form='unformatted')
         if (ios /= 0) then
            if (show_write_cache) write(*,*) 'write_cache failed to open ', trim(cache_filename)
            return
         end if

         if (show_write_cache) write(*,'(a)') 'write ' // trim(cache_filename)
         
         write(io_unit)  &
            reaction_Name(ir), which, cache_version, nrattab,  &
            rattab_thi, rattab_tlo, rattab_tstp

         do j = 1, nrattab
            write(io_unit) rattab(i,j)
         end do         

         close(io_unit)
      
      end subroutine write_reaction_to_cache

      
      subroutine do_show_reaction_from_cache(cache_filename, ierr) 
         character (len=*) :: cache_filename
         integer, intent(out) :: ierr
         
         integer :: version, nrattab, which
         real(dp) :: rattab_thi, rattab_tlo, rattab_tstp, rate, T8, logT
         integer :: ios, ir, i, j, io_unit
         real(dp), parameter :: tiny = 1e-6
         character (len=maxlen_reaction_Name) :: name

         ierr = 0
         io_unit = 37
         ios = 0
         open(unit=io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
         if (ios /= 0) then
            write(*,*) 'read cache failed for open ' // trim(cache_filename)
            return
         end if       
         
         read(io_unit, iostat=ios)  &
               name, which, version, nrattab,  &
               rattab_thi, rattab_tlo, rattab_tstp
         if (ios /= 0) then
            write(*,*) 'read cache failed for read header ' // trim(cache_filename)
            close(io_unit)
            return
         end if       
         
         write(*,'(a)') '#    T8     rate'
         write(*,*)
         write(*,*) nrattab
         
         do j = 1, nrattab
            read(io_unit, iostat=ios) rate
            if (ios /= 0) then
               write(*,*) 'read cache failed for reaction data', j
               close(io_unit)
               return
            end if
            logT = rattab_tlo + dble(j-1)*rattab_tstp
            T8 = exp10_cr(logT - 8d0)
            write(*,'(1pd26.16,3x,1pd26.16)') T8, rate
         end do         
         write(*,*)
         
         close(io_unit)
      
      end subroutine do_show_reaction_from_cache


      subroutine get_net_rates_for_tables( &
               which_rates, reaction_id, logT, btemp, num_reactions, rates, ierr)
         use ratelib, only: tfactors
         use raw_rates, only: set_raw_rates
         use utils_lib, only: is_bad_num

         real(dp), intent(in) :: logT, btemp
         integer, intent(in) :: num_reactions, reaction_id(:), which_rates(:)
         real(dp), intent(out) :: rates(:)
         integer, intent(out) :: ierr
      
         integer :: i, ir
         type (T_Factors), target ::  tfs
         type (T_Factors), pointer :: tf

         include 'formats.dek'
         
         ierr = 0

         tf => tfs
         call tfactors(tf, logT, btemp)
         call set_raw_rates( &
               num_reactions, reaction_id, which_rates, btemp, tf, rates, ierr)
         if (ierr /= 0) return

         do i = 1, num_reactions
            ir = reaction_id(i)
            if (ir <= 0) cycle
            if (is_bad_num(rates(i))) then
               write(*,2) trim(reaction_Name(ir)) // ' rates', i, rates(i)
               stop 'get_net_rates_for_tables'
            end if
         end do
               
      end subroutine get_net_rates_for_tables
      
      
      subroutine do_eval_reaclib_21( &
            ir, temp, den, rate_raw, reverse_rate_raw, ierr)
         use raw_rates, only: get_reaclib_rate_and_dlnT
         integer, intent(in) :: ir ! reaction_id
         real(dp), intent(in) :: temp, den
         real(dp), intent(out) :: rate_raw(:), reverse_rate_raw(:)
         integer, intent(out) :: ierr
         
         real(dp) :: lambda, dlambda_dlnT, rlambda, drlambda_dlnT
         
         include 'formats.dek'
         
         ierr = 0
         call get_reaclib_rate_and_dlnT( &
            ir, temp, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
         if (ierr /= 0) return
                  
         if (reaction_ye_rho_exponents(2,ir) /= 1) then
            ierr = -1
            return
         end if
            
         rate_raw(i_rate) = lambda*den
         rate_raw(i_rate_dT) = dlambda_dlnT*den/temp
         rate_raw(i_rate_dRho) = lambda
            
         reverse_rate_raw(i_rate) = rlambda
         reverse_rate_raw(i_rate_dT) = drlambda_dlnT/temp
         reverse_rate_raw(i_rate_dRho) = 0d0
         
         return

!$omp critical
         write(*,1) 'do_eval_reaclib_21 ' // trim(reaction_Name(ir))
         write(*,*)
         write(*,1) 'den', den
         write(*,1) 'temp', temp
         write(*,*)
         write(*,1) 'lambda', lambda
         write(*,1) 'dlambda_dlnT', dlambda_dlnT
         write(*,1) 'rate_raw', rate_raw(1:num_rvs)
         write(*,*)
         write(*,1) 'rlambda', rlambda
         write(*,1) 'drlambda_dlnT', drlambda_dlnT
         write(*,1) 'reverse_rate_raw', reverse_rate_raw(1:num_rvs)
         write(*,*)
!$omp end critical
         
      end subroutine do_eval_reaclib_21

      
      subroutine do_eval_reaclib_22( &
            ir, temp, den, rate_raw, reverse_rate_raw, ierr)
         use raw_rates, only: get_reaclib_rate_and_dlnT
         integer, intent(in) :: ir ! reaction_id
         real(dp), intent(in) :: temp, den
         real(dp), intent(out) :: rate_raw(:), reverse_rate_raw(:)
         integer, intent(out) :: ierr
         
         real(dp) :: lambda, dlambda_dlnT, rlambda, drlambda_dlnT
         
         include 'formats.dek'
         
         ierr = 0
         call get_reaclib_rate_and_dlnT( &
            ir, temp, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
         if (ierr /= 0) return
                  
         if (reaction_ye_rho_exponents(2,ir) /= 1) then
            ierr = -1
            return
         end if
            
         rate_raw(i_rate) = lambda*den
         rate_raw(i_rate_dT) = dlambda_dlnT*den/temp
         rate_raw(i_rate_dRho) = lambda
            
         reverse_rate_raw(i_rate) = rlambda*den
         reverse_rate_raw(i_rate_dT) = drlambda_dlnT*den/temp
         reverse_rate_raw(i_rate_dRho) = rlambda
         
         return

!$omp critical
         write(*,1) 'do_eval_reaclib_22 ' // trim(reaction_Name(ir))
         write(*,*)
         write(*,1) 'den', den
         write(*,1) 'temp', temp
         write(*,*)
         write(*,1) 'lambda', lambda
         write(*,1) 'dlambda_dlnT', dlambda_dlnT
         write(*,1) 'rate_raw', rate_raw(1:num_rvs)
         write(*,*)
         write(*,1) 'rlambda', rlambda
         write(*,1) 'drlambda_dlnT', drlambda_dlnT
         write(*,1) 'reverse_rate_raw', reverse_rate_raw(1:num_rvs)
         write(*,*)
         !stop 'do_eval_reaclib_22'
!$omp end critical
         
      end subroutine do_eval_reaclib_22
      

      end module rates_support

