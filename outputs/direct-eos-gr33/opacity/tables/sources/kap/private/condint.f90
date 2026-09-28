! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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


      module condint
      
      use const_def, only: dp
      use crlibm_lib
      
      implicit none

#ifdef offload
      !dir$ options /offload_attribute_target=mic
#endif            

      integer, parameter :: num_logTs=29, num_logRhos=71, num_logzs=15
      !!! NB: These parameters must be consistent with the table "condtabl.d"!
      logical :: initialized = .false.
      
      real :: logTs(num_logTs), logRhos(num_logRhos), logzs(num_logzs)
      real, target :: f_ary(4*num_logRhos*num_logTs*num_logzs) ! for bicubic splines
      real, pointer :: f(:,:,:,:)
      integer :: ilinx(num_logzs), iliny(num_logzs)
      
      
      contains
      
		
		subroutine init_potekhin(ierr)
   	   use kap_def, only: kap_dir
		   use interp_2d_lib_sg, only: interp_mkbicub_sg
		   integer, intent(out) :: ierr
		   
         character (len=256) :: filename
         integer :: read_err, iz, it, ir, shift
         integer :: ibcxmin                   ! bc flag for x=xmin
         real :: bcxmin(num_logTs)               ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real :: bcxmax(num_logTs)               ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real :: bcymin(num_logRhos)               ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real :: bcymax(num_logRhos)               ! bc data vs. x at y=ymax
         real :: Z
         real, pointer :: f1(:)
         
         ierr = 0
         if (initialized) return
         
         shift = 4*num_logRhos*num_logTs
         f(1:4,1:num_logRhos,1:num_logTs,1:num_logzs) => f_ary(1:shift*num_logzs)
                  
         filename = trim(kap_dir) // '/condtabl.data'
         open(1,file=trim(filename),status='OLD',iostat=ierr)
         if (ierr /= 0) then
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*) 'NOTICE: missing ' // trim(filename)
            write(*,*) 'Please remove the directory mesa/data/kap_data,'
            write(*,*) 'and rerun the mesa ./install script.'
            write(*,*)
            stop 1
         end if
         !print*,'Reading thermal conductivity data...'
         read_err = 0
         read(1,'(A)',iostat=read_err) ! skip the first line
         if (read_err /= 0) ierr = read_err
         do iz = 1, num_logzs
            read (1,*,iostat=read_err) z, (logTs(it),it=1,num_logTs)
            if (read_err /= 0) ierr = read_err
            logzs(iz) = log10_cr(dble(z))
            do ir = 1, num_logRhos
               read(1,*,iostat=read_err) logRhos(ir), (f(1,ir,it,iz),it=1,num_logTs)
               if (read_err /= 0) ierr = read_err
            end do
         end do
         close(1)
         if (ierr /= 0) then
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*) 'NOTICE: error trying to read ' // trim(filename)
            write(*,*) 'Please remove the directory mesa/data/kap_data,'
            write(*,*) 'and rerun the mesa ./install script.'
            write(*,*)
            stop 1
         end if
         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(1:num_logTs) = 0
         ibcxmax = 0; bcxmax(1:num_logTs) = 0
         ibcymin = 0; bcymin(1:num_logRhos) = 0
         ibcymax = 0; bcymax(1:num_logRhos) = 0
         do iz = 1, num_logzs
            f1(1:shift) => f_ary(1+(iz-1)*shift:iz*shift) 
            call interp_mkbicub_sg( &
               logRhos, num_logRhos, logTs, num_logTs, f1, num_logRhos, &
               ibcxmin, bcxmin, ibcxmax, bcxmax, &
               ibcymin, bcymin, ibcymax, bcymax, &
               ilinx(iz), iliny(iz), ierr)
            if (ierr /= 0) then
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*) 'NOTICE: error in ' // trim(filename)
               write(*,*) 'Please report the problem.'
               write(*,*)
               stop 1
            end if
         end do
         initialized = .true.
#ifdef offload
         f1 => f_ary
         call init_potekhin_on_coprocessor( &
		      logTs, logRhos, logzs, f1, ilinx, iliny)
#endif
		end subroutine init_potekhin
		
#ifdef offload
		subroutine init_potekhin_on_coprocessor( & ! on host
		      logTs, logRhos, logzs, f1, ilinx, iliny)
         real, intent(in) :: logTs(num_logTs), logRhos(num_logRhos), logzs(num_logzs)
         real, pointer, intent(in) :: f1(:)
         integer, intent(in) :: ilinx(num_logzs), iliny(num_logzs)
         !dir$ offload target(mic) in(logTs, logRhos, logzs, f1, ilinx, iliny)
         call do_init_potekhin_on_coprocessor( &
		      logTs, logRhos, logzs, f1, ilinx, iliny)
		end subroutine init_potekhin_on_coprocessor

		subroutine do_init_potekhin_on_coprocessor( & ! on mic
		      logTs_in, logRhos_in, logzs_in, f1, ilinx_in, iliny_in)
         real, intent(in) :: &
            logTs_in(num_logTs), logRhos_in(num_logRhos), logzs_in(num_logzs)
         real, pointer, intent(in) :: f1(:)
         integer, intent(in) :: ilinx_in(num_logzs), iliny_in(num_logzs)
         integer :: sz, i
         logTs = logTs_in
         logRhos = logRhos_in
         logzs = logzs_in
         ilinx = ilinx_in
         iliny = iliny_in
         sz = size(f1,dim=1)
         do i=1,sz
            f_ary(i) = f1(i)
         end do
         f(1:4,1:num_logRhos,1:num_logTs,1:num_logzs) => &
            f_ary(1:num_logRhos*num_logTs*num_logzs)
         initialized = .true.
      end subroutine do_init_potekhin_on_coprocessor
#endif
      
      subroutine do_electron_conduction( &
            zbar, logRho_in, logT_in, kap, dlogK_dlogRho, dlogK_dlogT, ierr)
         use utils_lib, only: is_bad_num
         real, intent(in) :: zbar, logRho_in, logT_in
         real(dp), intent(out) :: kap, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         integer :: iz, iz1, iz2, shift
         real :: zlog, logRho, logT
         real(dp) :: alfa, beta, logK, &
            logK1, kap1, dlogK1_dlogRho, dlogK1_dlogT, &
            logK2, kap2, dlogK2_dlogRho, dlogK2_dlogT
         real, pointer :: f1(:)
            
         include 'formats.dek'
         
         ierr = 0
         shift = 4*num_logRhos*num_logTs

         logRho = max(logRhos(1),min(logRhos(num_logRhos),logRho_in))
         logT = max(logTs(1),min(logTs(num_logTs),logT_in))
         zlog = max(logzs(1),min(logzs(num_logzs),log10_cr(max(1d-30,zbar))))
         
         if (zlog <= logzs(1)) then ! use 1st
            call get1(1, logK, kap, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         if (zlog >= logzs(num_logzs)) then ! use last
            call get1(num_logzs, logK, kap, dlogK_dlogRho, dlogK_dlogT, ierr)
            return
         end if
         
         iz1 = -1
         do iz = 2, num_logzs
            if (zlog >= logzs(iz-1) .and. zlog <= logzs(iz)) then
               iz1 = iz-1; iz2 = iz; exit
            end if
         end do
         if (iz1 < 0) then
            write(*,2) 'num_logzs', num_logzs
            do iz = 1, num_logzs
               write(*,2) 'logzs(iz)', iz, logzs(iz)
            end do
            write(*,1) 'zlog', zlog
            write(*,*) 'confusion in do_electron_conduction'
            stop 1
         end if
         
         call get1(iz1, logK1, kap1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) then
            write(*,*) 'interp failed for iz1 in do_electron_conduction', iz1, logRho, logT
            stop 1
         end if
         
         call get1(iz2, logK2, kap2, dlogK2_dlogRho, dlogK2_dlogT, ierr)
         if (ierr /= 0) then
            write(*,*) 'interp failed for iz2 in do_electron_conduction', iz2, logRho, logT
            stop 1
         end if
         
         ! linear interpolation in zlog
         alfa = (zlog - logzs(iz1)) / (logzs(iz2) - logzs(iz1))
         beta = 1-alfa
         logK = alfa*logK2 + beta*logK1
         kap = exp10_cr(logK)
         dlogK_dlogRho = alfa*dlogK2_dlogRho + beta*dlogK1_dlogRho
         dlogK_dlogT = alfa*dlogK2_dlogT + beta*dlogK1_dlogT
         
         
         contains
         
         
#ifdef offload
         !dir$ attributes offload: mic :: get1
#endif      
         subroutine get1(iz, logK, kap, dlogK_dlogRho, dlogK_dlogT, ierr)
            use interp_2d_lib_sg, only: interp_evbicub_sg
            integer, intent(in) :: iz
            real(dp), intent(out) :: logK, kap, dlogK_dlogRho, dlogK_dlogT
            integer, intent(out) :: ierr
            integer :: ict(6)                    ! code specifying output desired
            real :: fval(6)                      ! output data
            logical, parameter :: dbg = .false.
            include 'formats'
            ierr = 0
            ict(1:3) = 1
            ict(4:6) = 0
            f1(1:shift) => f_ary(1+(iz-1)*shift:iz*shift) 
            call interp_evbicub_sg( &
               logRho, logT, logRhos, num_logRhos, logTs, num_logTs, &
               ilinx(iz), iliny(iz), f1, num_logRhos, ict, fval, ierr)
            if (ierr /= 0) return
            ! fval(1) = CK; fval(2) = DRK = dCK/dlogRho; fval(3) = DTK = dCK/dlogT
            ! chi = thermal conductivity, = 10**CK (cgs units)
   		   ! conduction opacity kappa = 16*boltz_sigma*T^3 / (3*rho*chi)
   		   ! logK = 3*logT - logRho - CK + log10(16*boltz_sigma/3)
   		   ! log10(16*boltz_sigma/3) = -3.51937938116756d0
            logK = 3*dble(logT) - dble(logRho) - fval(1) - 3.51937938116756d0
            kap = exp10_cr(logK)
            if (dbg) write(*,2) 'do_electron_conduction', iz, dble(logK)
            dlogK_dlogRho = - 1 - dble(fval(2))
            dlogK_dlogT = 3 - dble(fval(3))         
         end subroutine get1


      end subroutine do_electron_conduction

#ifdef offload
      !dir$ end options
#endif

      end module condint
