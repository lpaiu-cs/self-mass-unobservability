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

      module eval_weak
      
      use rates_def_mic
      use const_def, only: dp
      use crlibm_lib
      
      implicit none


      contains
		
#ifdef offload
      !dir$ attributes offload: mic :: do_eval_weak_reaction_info
#endif         
		subroutine do_eval_weak_reaction_info( &
		      n, ids, T9, YeRho, &
		      eta, d_eta_dlnT, d_eta_dlnRho, &
		      ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
		      lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
		      lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
		      lambda, dlambda_dlnT, dlambda_dlnRho, &
		      Q, dQ_dlnT, dQ_dlnRho, &
		      Qneu, dQneu_dlnT, dQneu_dlnRho, &
		      ierr)
         use const_def, only: ln10, kerg, mev_to_ergs
         use chem_def
         use interp_1d_def
         use utils_lib, only: is_bad_real, is_bad_num
		   integer, intent(in) :: n, ids(:)
		   real(dp), intent(in) :: T9, YeRho, eta, d_eta_dlnT, d_eta_dlnRho
		   real(dp), dimension(:), intent(out), pointer :: &
		      ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
		      lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
		      lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
		      lambda, dlambda_dlnT, dlambda_dlnRho, &
		      Q, dQ_dlnT, dQ_dlnRho, &
		      Qneu, dQneu_dlnT, dQneu_dlnRho
		   integer, intent(out) :: ierr
		   
         real(dp) :: alfa, beta, alfa_hi_Z, beta_hi_Z
         integer :: i, ir, cid

         include 'formats.dek'
		   
	      call do_eval_weaklib_reaction_info( &
		      n, ids, T9, YeRho, &
		      eta, d_eta_dlnT, d_eta_dlnRho, &
		      ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
		      lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
		      lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
		      lambda, dlambda_dlnT, dlambda_dlnRho, &
		      Q, dQ_dlnT, dQ_dlnRho, &
		      Qneu, dQneu_dlnT, dQneu_dlnRho, &
		      ierr)
	      if (ierr /= 0) then
	         !write(*,*) 'failed in do_eval_weaklib_reaction_info'
	         return
	      end if

	      if (T9 >= max(T9_weaklib_full_on, T9_weaklib_full_on_hi_Z)) return
	      
	      ! revise lambda using reaclib rate for low T		   
         ! alfa is fraction from weaklib
	      if (T9 >= T9_weaklib_full_on) then
	         alfa = 1d0
		   else if (T9 > T9_weaklib_full_off) then
		      alfa = (T9 - T9_weaklib_full_off)/ &
		             (T9_weaklib_full_on - T9_weaklib_full_off)
		   else
		      alfa = 0d0
		   end if
		   beta = 1d0 - alfa ! beta is fraction from reaclib

	      if (T9 >= T9_weaklib_full_on_hi_Z) then
	         alfa_hi_Z = 1d0
		   else if (T9 > T9_weaklib_full_off_hi_Z) then
		      alfa_hi_Z = (T9 - T9_weaklib_full_off_hi_Z)/ &
		             (T9_weaklib_full_on_hi_Z - T9_weaklib_full_off_hi_Z)
		   else
		      alfa_hi_Z = 0d0
		   end if
		   beta_hi_Z = 1d0 - alfa_hi_Z ! beta is fraction from reaclib

			do i = 1, n
			   if (ids(i) <= 0) cycle
			   ir = weak_reaclib_id(ids(i))
			   if (ir <= 0) cycle
	         if (weak_lowT_rate(ir) <= 0d0) cycle
	         cid = weak_reaction_info(1,ir)
	         if (cid <= 0) cycle
	         if (chem_isos% Z(cid) >= weaklib_blend_hi_Z) then
      	      lambda(i) = alfa_hi_Z*lambda(i) + beta_hi_Z*weak_lowT_rate(ir)
         	   dlambda_dlnT(i) = alfa_hi_Z*dlambda_dlnT(i)
         	   dlambda_dlnRho(i) = alfa_hi_Z*dlambda_dlnRho(i)
	         else
      	      lambda(i) = alfa*lambda(i) + beta*weak_lowT_rate(ir)
         	   dlambda_dlnT(i) = alfa*dlambda_dlnT(i)
         	   dlambda_dlnRho(i) = alfa*dlambda_dlnRho(i)
	         end if
      	end do
		   
		end subroutine do_eval_weak_reaction_info

		
#ifdef offload
      !dir$ attributes offload: mic :: do_eval_weaklib_reaction_info
#endif         
		subroutine do_eval_weaklib_reaction_info( &
		      n, ids, T9_in, YeRho_in, &
		      eta, d_eta_dlnT, d_eta_dlnRho, &
		      ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
		      lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
		      lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
		      lambda, dlambda_dlnT, dlambda_dlnRho, &
		      Q, dQ_dlnT, dQ_dlnRho, &
		      Qneu, dQneu_dlnT, dQneu_dlnRho, &
		      ierr)
         use const_def, only: ln10, kerg, mev_to_ergs
         use chem_def
         use interp_1d_def
         use utils_lib, only: is_bad_real, is_bad_num
		   integer, intent(in) :: n, ids(:)
		   real(dp), intent(in) :: T9_in, YeRho_in, eta, d_eta_dlnT, d_eta_dlnRho
		   real(dp), dimension(:), intent(out), pointer :: &
		      ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
		      lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
		      lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, &
		      lambda, dlambda_dlnT, dlambda_dlnRho, &
		      Q, dQ_dlnT, dQ_dlnRho, &
		      Qneu, dQneu_dlnT, dQneu_dlnRho
		   integer, intent(out) :: ierr
		   
         logical, parameter :: dbg = .false.
		   
		   real :: T, T9, YeRho, lYeRho, f, dfdx, dfdy, conv
		   real :: delta_T9, dT9, dlYeRho, delta_lYeRho, y_alfa, y_beta, x_alfa, x_beta
		   integer :: iT9, ilYeRho, i, ir, in, out

         real, parameter :: sixth = 1e0 / 6e0
         real, parameter :: z36th = 1e0 / 36e0

         integer :: ix, jy          ! target cell in the spline data
         real :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(ix), x1 = xs(ix+1)
         real :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(jy), y1 = ys(jy+1)
         real :: xp, xpi, xp2, xpi2, cx, cxi, hx2, cxd, cxdi, hx, hxi
         real :: yp, ypi, yp2, ypi2, cy, cyi, hy2, cyd, cydi, hy, hyi
         
         logical :: neg
         real(dp) :: decay, capture, Qx, Qn, mue, d_mue_dlnRho, d_mue_dlnT
	      character(len=iso_name_length) :: weak_lhs, weak_rhs
         integer, parameter :: nwork = pm_work_size

         include 'formats.dek'

         ierr = 0
			
			T9 = real(T9_in)
			YeRho = real(YeRho_in)
			lYeRho = real(log10_cr(YeRho_in))
			if (is_bad_real(lYeRho)) then
			   ierr = -1
			   return
			   
		      write(*,1) 'lYeRho', lYeRho
		      write(*,1) 'YeRho_in', YeRho_in
		      write(*,1) 'log10_cr(YeRho_in)', log10_cr(YeRho_in)
			   !stop 'weak lYeRho'
			end if
			
			if (n == 0) then
		      ldecay=-100; d_ldecay_dT9=0; d_ldecay_dlYeRho=0
		      lcapture=-100; d_lcapture_dT9=0; d_lcapture_dlYeRho=0
		      lneutrino=-100; d_lneutrino_dT9=0; d_lneutrino_dlYeRho=0
		      write(*,*) 'problem in eval_weak_reaction_info: n == 0'
		      write(*,2) 'n', n
		      write(*,*) 'T9 < weak_reaction_T9s(1)', T9 < weak_reaction_T9s(1)
		      write(*,*) 'lYeRho < weak_reaction_lYeRhos(1)', lYeRho < weak_reaction_lYeRhos(1)
		      write(*,1) 'weak_reaction_T9s(1)', weak_reaction_T9s(1)
		      write(*,1) 'T9', T9
		      write(*,1) 'lYeRho', lYeRho
		      write(*,1) 'weak_reaction_lYeRhos(1)', weak_reaction_lYeRhos(1)
		      
			   return
			end if
			
			! clip small values to edge of table
			if (T9 < weak_reaction_T9s(1)) &
			   T9 = weak_reaction_T9s(1)
			if (lYeRho < weak_reaction_lYeRhos(1)) &
			   lYeRho = weak_reaction_lYeRhos(1)
			
			! clip large values to edge of table
			if (T9 > weak_reaction_T9s(weak_num_T9)) &
			   T9 = weak_reaction_T9s(weak_num_T9)
			if (lYeRho > weak_reaction_lYeRhos(weak_num_lYeRho)) &
			   lYeRho = weak_reaction_lYeRhos(weak_num_lYeRho)
			
			xget = T9
			yget = lYeRho
			if (weak_bicubic) then
			   call setup_for_bicubic_interpolations
			else
			   call setup_for_linear_interp
			endif
			
			! convert to MeV
			conv = real(kerg/mev_to_ergs)
			T = T9*1e9
      	mue = eta*conv*T
      	d_mue_dlnRho = d_eta_dlnRho*conv*T
      	if (d_eta_dlnT == 0) then
      	   d_mue_dlnT = 0
      	else
      	   d_mue_dlnT = d_eta_dlnT*conv*T + mue
      	end if

			do i = 1, n
			
			   ir = ids(i)
			   if (ir <= 0) cycle
			   if (ir > size(weak_reactions_data,dim=5)) then
			      call show_stuff
			      write(*,*) 'bad ir', ir, i
			      stop 1
			      ierr = -1
			      return
			   end if

			   neg = (mod(ir,2)==1) ! neg is true for electron capture and positron emission
			   
			   if (weak_bicubic) then
			   
   			   call do_bicubic_interpolations( &
   			      weak_reactions_data(1:4,1:weak_num_T9,1:weak_num_lYeRho,i_ldecay,ir), f, dfdx, dfdy)
   			   ldecay(i)=dble(f); d_ldecay_dT9(i)=dble(dfdx); d_ldecay_dlYeRho(i)=dble(dfdy)
			   
   			   call do_bicubic_interpolations( &
   			      weak_reactions_data(1:4,1:weak_num_T9,1:weak_num_lYeRho,i_lcapture,ir), f, dfdx, dfdy)
   			   lcapture(i)=dble(f); d_lcapture_dT9(i)=dble(dfdx); d_lcapture_dlYeRho(i)=dble(dfdy)
			   
   			   call do_bicubic_interpolations( &
   			      weak_reactions_data(1:4,1:weak_num_T9,1:weak_num_lYeRho,i_lneutrino,ir), f, dfdx, dfdy)
   			   lneutrino(i)=dble(f); d_lneutrino_dT9(i)=dble(dfdx); d_lneutrino_dlYeRho(i)=dble(dfdy)
   			
   			else
   			
   			   call do_linear_interp( &
   			      weak_reactions_data(1:4,1:weak_num_T9,1:weak_num_lYeRho,i_ldecay,ir), &
   			      ldecay(i), d_ldecay_dT9(i), d_ldecay_dlYeRho(i), ierr)
			   
   			   call do_linear_interp( &
   			      weak_reactions_data(1:4,1:weak_num_T9,1:weak_num_lYeRho,i_lcapture,ir), &
   			      lcapture(i), d_lcapture_dT9(i), d_lcapture_dlYeRho(i), ierr)
			   
   			   call do_linear_interp( &
   			      weak_reactions_data(1:4,1:weak_num_T9,1:weak_num_lYeRho,i_lneutrino,ir), &
   			      lneutrino(i), d_lneutrino_dT9(i), d_lneutrino_dlYeRho(i), ierr)
   			
   			end if
			   
	         decay = exp10_cr(ldecay(i))
	         capture = exp10_cr(lcapture(i))
      	   lambda(i) = decay + capture
      	   ! d_dlnT = T9*d_dT9
      	   dlambda_dlnT(i) = ln10*T9*(decay*d_ldecay_dT9(i) + capture*d_lcapture_dT9(i))
      	   ! d_dlnRho = rho*d_drho = ye*rho*d_d(ye*rho) = d_dlnYeRho = ln10*d_dlYeRho
      	   dlambda_dlnRho(i) = ln10*(decay*d_ldecay_dlYeRho(i) + capture*d_lcapture_dlYeRho(i))
      	   
      	   in = weak_lhs_nuclide_id(ir)
      	   out = weak_rhs_nuclide_id(ir)
      	   Qx = chem_isos% mass_excess(in) - chem_isos% mass_excess(out)
      	   
      	   if (neg) then ! electron capture and positron emission
      	      Q(i) = Qx + mue 
      	      dQ_dlnT(i) = d_mue_dlnT
      	      dQ_dlnRho(i) = d_mue_dlnRho
      	   else ! positron capture and electron emission
      	      Q(i) = Qx - mue 
      	      dQ_dlnT(i) = -d_mue_dlnT
      	      dQ_dlnRho(i) = -d_mue_dlnRho
      	   end if
      	   
      	   if (lambda(i) < 1e-30) then
         	   Qneu(i) = 0
         	   dQneu_dlnT(i) = 0
         	   dQneu_dlnRho(i) = 0
      	   else
         	   Qn = exp10_cr(lneutrino(i))/lambda(i)
         	   Qneu(i) = Qn
         	   dQneu_dlnT(i) = ln10*T9*Qn*d_lneutrino_dT9(i) - dlambda_dlnT(i)*Qn/lambda(i)
         	   dQneu_dlnRho(i) = ln10*Qn*d_lneutrino_dlYeRho(i) - dlambda_dlnRho(i)*Qn/lambda(i)
      	   end if

            if (is_bad_num(Qneu(i))) then
               ierr = -1
               return

               write(*,2) 'lambda', i, lambda(i)
               write(*,2) 'lneutrino', i, lneutrino(i)
               write(*,2) 'Qneu', i, Qneu(i)
               call show_stuff
            end if
            
			end do
						
         if (is_bad_real(lYeRho)) then
            ierr = -1
            return
            call show_stuff
         end if
      
      
         contains
         
#ifdef offload
      !dir$ attributes offload: mic :: show_stuff
#endif         
         subroutine show_stuff
            include 'formats.dek'
   	      write(*,1) 'T9', T9
   	      write(*,1) 'lYeRho', lYeRho
   	      write(*,1) 'eta', eta
   	      write(*,1) 'mue', mue
   	      write(*,*)
   			do i = 1, n
			      ir = ids(i)
         	   in = weak_lhs_nuclide_id(i)
         	   out = weak_rhs_nuclide_id(i)
         	   if (.true.) then
         	      write(*,4) 'ir, i, n', ir, i, n
         	      if (ir <= 0 .or. ir > ubound(weak_lhs_nuclide_id,dim=1)) cycle
      		      weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
      		      weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
      		      write(*,'(a30,3i5)') weak_lhs // weak_rhs, ir, i, n
         	      !write(*,1) 'chem_isos% mass_excess(in)', chem_isos% mass_excess(in)
         	      !write(*,1) 'chem_isos% mass_excess(out)', chem_isos% mass_excess(out)
         	      write(*,2) 'Qx', i, chem_isos% mass_excess(in) - chem_isos% mass_excess(out)
         	      write(*,2) 'Q', i, Q(i)
         	      write(*,2) 'Qneu', i, Qneu(i)
         	      write(*,*)
         	   end if
   			end do
         	stop 1
         end subroutine show_stuff         
         
#ifdef offload
      !dir$ attributes offload: mic :: find_location
#endif         
         subroutine find_location ! set ix, jy; x is T9; y is lYeRho
            integer i, j
            real :: del
            include 'formats.dek'                       
            ! x0 <= T9 <= x1
            ix = weak_num_T9-1 ! since weak_num_T9 is small, just do a linear search
            do i = 2, weak_num_T9-1
               if (T9 > weak_reaction_T9s(i)) cycle
               ix = i-1
               exit
            end do           
            ! assume lYeRho values are evenly spaced
            ! y0 <= lYeRho <= y1
            del = (weak_reaction_lYeRhos(weak_num_lYeRho) - weak_reaction_lYeRhos(1))/(weak_num_lYeRho-1)
            jy = int((lYeRho - weak_reaction_lYeRhos(1)) / del) + 1
            if (jy <= 0) jy = 1
            if (jy == weak_num_lYeRho) jy = weak_num_lYeRho-1            
            x0 = weak_reaction_T9s(ix); x1 = weak_reaction_T9s(ix+1)
            y0 = weak_reaction_lYeRhos(jy); y1 = weak_reaction_lYeRhos(jy+1)
         end subroutine find_location
                  
#ifdef offload
      !dir$ attributes offload: mic :: setup_for_bicubic_interpolations
#endif         
         subroutine setup_for_bicubic_interpolations
            integer i, j
            real :: del
            
            include 'formats.dek'
            
            call find_location
            
            ! set factors for interpolation
            
            hx=x1-x0
            hxi=1.0/hx
            hx2=hx*hx
   
            xp=(xget-x0)*hxi
            xpi=1.0-xp
            xp2=xp*xp
            xpi2=xpi*xpi

            cx=xp*(xp2-1.0)
            cxi=xpi*(xpi2-1.0)
            cxd=3.0*xp2-1.0
            cxdi=-3.0*xpi2+1.0
   
            hy=y1-y0
            hyi=1.0/hy
            hy2=hy*hy
   
            yp=(yget-y0)*hyi
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
            
            if (dbg) then
               write(*,2) 'T9', ix, x0, T9, x1
               write(*,2) 'lYeRho', jy, y0, lYeRho, y1
               write(*,1) 'xpi', xpi
               write(*,1) 'ypi', ypi
               write(*,*)
            end if

         end subroutine setup_for_bicubic_interpolations         
      
#ifdef offload
      !dir$ attributes offload: mic :: do_bicubic_interpolations
#endif         
         subroutine do_bicubic_interpolations(fin, fval, df_dx, df_dy)
            ! derived from routines in the PSPLINE package written by Doug McCune 
            real, dimension(:,:,:) :: fin ! the spline data array, dimensions (4, nx, ny)
            real, intent(out) :: fval, df_dx, df_dy
         
            ! bicubic spline interpolation
            fval = &
               xpi*( &
                  ypi*fin(1,ix,jy)  +yp*fin(1,ix,jy+1)) &
                  +xp*(ypi*fin(1,ix+1,jy)+yp*fin(1,ix+1,jy+1)) &
               +sixth*hx2*( &
                  cxi*(ypi*fin(2,ix,jy) +yp*fin(2,ix,jy+1))+ &
                  cx*(ypi*fin(2,ix+1,jy)+yp*fin(2,ix+1,jy+1))) &
               +sixth*hy2*( &
                  xpi*(cyi*fin(3,ix,jy) +cy*fin(3,ix,jy+1))+ &
                  xp*(cyi*fin(3,ix+1,jy)+cy*fin(3,ix+1,jy+1))) &
               +z36th*hx2*hy2*( &
                  cxi*(cyi*fin(4,ix,jy) +cy*fin(4,ix,jy+1))+ &
                  cx*(cyi*fin(4,ix+1,jy)+cy*fin(4,ix+1,jy+1)))

            include 'formats.dek'
            if (dbg) then
               write(*,1) 'fin(1,ix,jy)', fin(1,ix,jy)
               write(*,1) 'fin(1,ix,jy+1)', fin(1,ix,jy+1)
               write(*,1) 'fin(1,ix+1,jy)', fin(1,ix+1,jy)
               write(*,1) 'fin(1,ix+1,jy+1)', fin(1,ix+1,jy+1)
               write(*,1) 'fval', fval
               
               write(*,*)
               stop 'debug: do_bicubic_interpolations'
            end if

            ! derivatives of bicubic splines
            df_dx = &
               hxi*( &
                  -(ypi*fin(1,ix,jy)  +yp*fin(1,ix,jy+1)) &
                  +(ypi*fin(1,ix+1,jy)+yp*fin(1,ix+1,jy+1))) &
               +sixth*hx*( &
                  cxdi*(ypi*fin(2,ix,jy) +yp*fin(2,ix,jy+1))+ &
                  cxd*(ypi*fin(2,ix+1,jy)+yp*fin(2,ix+1,jy+1))) &
               +sixth*hxi*hy2*( &
                  -(cyi*fin(3,ix,jy)  +cy*fin(3,ix,jy+1)) &
                  +(cyi*fin(3,ix+1,jy)+cy*fin(3,ix+1,jy+1))) &
               +z36th*hx*hy2*( &
                  cxdi*(cyi*fin(4,ix,jy) +cy*fin(4,ix,jy+1))+ &
                  cxd*(cyi*fin(4,ix+1,jy)+cy*fin(4,ix+1,jy+1)))

            df_dy = &
               hyi*( &
                  xpi*(-fin(1,ix,jy) +fin(1,ix,jy+1))+ &
                  xp*(-fin(1,ix+1,jy)+fin(1,ix+1,jy+1))) &
               +sixth*hx2*hyi*( &
                  cxi*(-fin(2,ix,jy) +fin(2,ix,jy+1))+ &
                  cx*(-fin(2,ix+1,jy)+fin(2,ix+1,jy+1))) &
               +sixth*hy*( &
                  xpi*(cydi*fin(3,ix,jy) +cyd*fin(3,ix,jy+1))+ &
                  xp*(cydi*fin(3,ix+1,jy)+cyd*fin(3,ix+1,jy+1))) &
               +z36th*hx2*hy*( &
                  cxi*(cydi*fin(4,ix,jy) +cyd*fin(4,ix,jy+1))+ &
                  cx*(cydi*fin(4,ix+1,jy)+cyd*fin(4,ix+1,jy+1)))

         end subroutine do_bicubic_interpolations

#ifdef offload
      !dir$ attributes offload: mic :: setup_for_linear_interp
#endif         
         subroutine setup_for_linear_interp
            include 'formats.dek'
         
            call find_location
            
            dT9 = T9 - x0
            delta_T9 = x1 - x0
            x_beta = dT9 / delta_T9 ! fraction of x1 result
            x_alfa = 1 - x_beta ! fraction of x0 result
            if (x_alfa < 0 .or. x_alfa > 1) then
               write(*,1) 'weaklib: x_alfa', x_alfa
               write(*,1) 'T9', T9
               write(*,1) 'x0', x0
               write(*,1) 'x1', x1
               stop 1
            end if          
            
            dlYeRho = lYeRho - y0
            delta_lYeRho = y1 - y0
            y_beta = dlYeRho / delta_lYeRho ! fraction of y1 result
            y_alfa = 1 - y_beta ! fraction of y0 result     
            if (is_bad_real(y_alfa) .or. y_alfa < 0 .or. y_alfa > 1) then
               write(*,1) 'weaklib: y_alfa', y_alfa
               write(*,1) 'T9', T9
               write(*,1) 'x0', x0
               write(*,1) 'dT9', dT9
               write(*,1) 'delta_T9', delta_T9
               write(*,1) 'lYeRho', lYeRho
               write(*,1) 'y0', y0
               write(*,1) 'dlYeRho', dlYeRho
               write(*,1) 'y1', y1
               write(*,1) 'delta_lYeRho', delta_lYeRho
               write(*,1) 'y_beta', y_beta
               !stop 'weak setup_for_linear_interp'
            end if          
                   
         end subroutine setup_for_linear_interp         
      
#ifdef offload
      !dir$ attributes offload: mic :: do_linear_interp
#endif         
         subroutine do_linear_interp(f, fval, df_dx, df_dy, ierr)
            use interp_1d_lib
            use utils_lib, only: is_bad_real, is_bad_num
            real, dimension(:,:,:) :: f ! (4, nx, ny)         
            real(dp), intent(out) :: fval, df_dx, df_dy
            integer, intent(out) :: ierr
            
            real(dp) :: fx0, fx1, fy0, fy1
            
            include 'formats'
            
            ierr = 0
            
            fx0 = y_alfa*f(1,ix,jy) + y_beta*f(1,ix,jy+1)
            fx1 = y_alfa*f(1,ix+1,jy) + y_beta*f(1,ix+1,jy+1)
            
            fy0 = x_alfa*f(1,ix,jy) + x_beta*f(1,ix+1,jy)
            fy1 = x_alfa*f(1,ix,jy+1) + x_beta*f(1,ix+1,jy+1)
            
            fval = x_alfa*fx0 + x_beta*fx1
            df_dx = (fx1 - fx0)/(x1 - x0)
            df_dy = (fy1 - fy0)/(y1 - y0)    
            
            if (is_bad_num(fval)) then
               ierr = -1
               return

               write(*,1) 'x_alfa', x_alfa
               write(*,1) 'x_beta', x_beta
               write(*,1) 'fx0', fx0
               write(*,1) 'fx1', fx1
               write(*,1) 'y_alfa', y_alfa
               write(*,1) 'y_beta', y_beta
               write(*,1) 'f(1,ix,jy)', f(1,ix,jy)
               write(*,1) 'f(1,ix,jy+1)', f(1,ix,jy+1)
               !stop 'weak do_linear_interp'
            end if        
            
         end subroutine do_linear_interp

         
		end subroutine do_eval_weaklib_reaction_info
      



      end module eval_weak

