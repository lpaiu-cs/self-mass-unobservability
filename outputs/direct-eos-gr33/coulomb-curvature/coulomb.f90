!*******************************************************************************
!       Copyright (C) 1996-2022 Alan W. Irwin
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!*******************************************************************************

! routine to calculate EOS quantitities relevant to the
! Coulomb interaction using the PTEH approximation (ifcoulomb = 3, 4
! pteh theta calculated outside, otherwise fermi-dirac theta_e calculated
! outside) or DH smoothly joined by cubic to modified OCP
! (preferred, ifcoulomb = 5)
! or the Debye-Huckel limit (ifcoulomb = 1) or that limit corrected
! by tau(x) (ifcoulomb = 2) for the Coulomb
! free energy calculated for the partial ionization case.
! input variables:
! ifcoulomb controls whether Debye-Huckel (1), Debye-Huckel with tau(x)
!   correction (2), PTEH approximation (3 or 4 or 9),
!   DH smoothly joined using cubic with modified OCP (preferred) (5)
!   or alternative version of cubic with different end points (6),
!   DH for Gamma <1 and OCP for Gamma >= 1 (for illustrative purposes)
!   (ifcoulomb = 7 or 8, with different theta_e values calculated by the calling routine)
!   ifcoulomb = 3 means theta_e by pteh approximation outside
!   ifcoulomb != 3 means theta_e by usual fermi-dirac formula.
!   ifcouloumb <= 4 means sum0 calculated outside with
!   PTEH definition otherwise calculated outside with DeWitt definition.
!   N.B. this distinction is meaningless for ifcoulomb = 0, 1, or 2.
! if_dc = 0 means apply no diffraction correction to Lambda
! if_dc = 1 means apply diffraction correction to Lambda
!   n.b. tau(x) correction and diffraction correction are not
!   simultaneously allowed.
! sum0 = sum_ions nion, where nion is the number per unit volume,
!   and the sum taken over positive ions when PTEH definition is
!   used outside.  Otherwise n_e thetae has been added outside
!   (DeWitt definition).
! sum2 = n_e thetae + sum_ions Z^2 nion,
!   n_e is the electron number density, thetae is the degeneracy
!   correction (see PTEH or notes), Z is the charge on the ions.
! t = temperature in K
! if if_dc = 0, then
!   x = thetax*n_e = (n_e k T/P_e)*n_e, where n_e is the
!   number density of free electrons and P_e is associated pressure.
!     n.b. thetax*n_e is only relevant when the tau(x) correction is
!     being applied to the Debye-Huckle limit, i.e., ifcoulomb = 2
! if if_dc = 1, then
!   x = n_e thetae or possibly n_e
! output variables:
! fcoulomb = free energy per unit volume as a function of
!   sum0, sum2, t, and x.
! fcoulomb0 = partial fcoulomb(sum0, sum2, t, x) wrt sum0
! fcoulomb2 = partial fcoulomb(sum0, sum2, t, x) wrt sum2
! fcoulombt = partial fcoulomb(sum0, sum2, t, x) wrt t
! fcoulombx = partial fcoulomb(sum0, sum2, t, x) wrt x
! fcoulomb00, etc., higher order derivatives
! lambda = plasma interaction parameter (defined by PTEH and consistent
! with the simpler CG 15.60 criterion for single component plasmas.)
! gamma_e ~ gamma_e_i of DeWitt, paper II, with a correction
! that is only good (currently) for 1 >> gamma_e > Lambda.

!> This coulomb subroutine calculates EOS quantitities relevant to the
!> Coulomb interaction for a large variety of different free-energy models
!> for the Coulomb effect.
!>
!> \param[in] ifcoulomb PARAMETERS NEED DOCUMENTATION
!>
subroutine coulomb(ifcoulomb, if_dc, sum0, sum2, t, x,&
     fcoulomb, fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,&
     fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x,&
     fcoulomb22, fcoulomb2t, fcoulomb2x,&
     fcoulombtt, fcoulombtx, fcoulombxx, lambda, gamma_e)
  use mod_pi_fit, only: xdh10, xmocp10
  use mod_free_eos_constants, only: alpha2, avogadro, boltzmann, echarge, electron_mass, pi, ln10
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  integer, intent(in) :: ifcoulomb, if_dc
  real(fp_kind), intent(in) :: sum0, sum2, t, x
  real(fp_kind), intent(out) :: fcoulomb,&
       fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,&
       fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x,&
       fcoulomb22, fcoulomb2t, fcoulomb2x,&
       fcoulombtt, fcoulombtx, fcoulombxx,&
       lambda, gamma_e

  ! Internal variables
  real(fp_kind) gamma, zeta2,&
       dc_lambda,&
       dc_lambda2, dc_lambdat, dc_lambdax,&
       dc_lambda22, dc_lambda2t, dc_lambda2x,&
       dc_lambdatt, dc_lambdatx, dc_lambdaxx,&
       lambda2, lambdat, lambdax,&
       lambda22, lambda2t, lambda2x,&
       lambdatt, lambdatx, lambdaxx,&
       dh, ddh, ddh2, ocp, docp, docp2,&
       gc, dgc, dgc2, tau,&
       xvalue, xdh, xmocp
  real(fp_kind) dscmccon, dscmdcon,&
       dx, acap, bcap, ccap, dcap, ylo, yhi, dylo2, dyhi, dyhi2
  real(fp_kind), allocatable :: dtau(:), d2tau(:,:)

  logical ifdh_limit
  integer ifcoulomb_old
  ! Must be different from 5 or 6.
  data ifcoulomb_old/-1/

  ! size of dtau and d2tau dimensions
  integer, parameter :: ntau = 4

  ! PTEH values
  real(fp_kind), parameter :: acon = 0.89752_fp_kind
  ! n.b. changing 0.208 to 0.5 gives good LMS fit, but bad solar
  ! fit (at least by solar standards)
  ! pteh values give even worse solar fit so went to dh for
  ! solar conditions (log10(gamma) < -0.4) smoothly joined by cubic
  ! to modified ocp result.
  real(fp_kind), parameter :: dcon = 0.208_fp_kind*acon
  real(fp_kind), parameter :: gcon = -1._fp_kind/0.768_fp_kind
  ! fitting coefficients from 1995, DeWitt, Slattery, and Chabrier
  real(fp_kind), parameter :: a = -0.899126_fp_kind
  real(fp_kind), parameter :: b = 0.60712_fp_kind
  real(fp_kind), parameter :: c = -0.27998_fp_kind
  real(fp_kind), parameter :: s = 0.321308_fp_kind
  real(fp_kind), parameter :: f1 = -0.436484_fp_kind
  ! n.b. f1 is from HNC calculation
  real(fp_kind), parameter :: d = f1 -(a +b/s)
  real(fp_kind), parameter :: dscacon = -a
  real(fp_kind), parameter :: dscbcon = -b/s
  real(fp_kind), parameter :: dscscon = s
  real(fp_kind), parameter :: dscccon = -c
  real(fp_kind), parameter :: dscdcon = -d

  real(fp_kind), parameter :: lambda_const = 2._fp_kind*echarge*echarge*echarge*sqrt(pi/boltzmann)/boltzmann
  ! use eq. (36) of DeWitt, 1966 (J. Math. Phys. 7, 616), but
  ! use different definition of gamma_e, assume m_e/m_i << 1
  ! so that gamma_e_e = sqrt(2) gamma_e, gamma_e_i ~ gamma_e.
  ! (2 pi k/avogadro/h^2)^(-0.5)
  ! = (alpha2/avogadro^5)**(1/6)*avogdro**(3/6)
  ! = alpha2**(1/6)/avogadro**(1/3)
  real(fp_kind), parameter :: gamma_e_const = -3._fp_kind*sqrt(pi)/8._fp_kind
  ! include lambda as part of definition (see below)
  real(fp_kind), parameter :: dc_const = gamma_e_const*alpha2**(1._fp_kind/6._fp_kind)/&
       avogadro**(1._fp_kind/3._fp_kind)/sqrt(4._fp_kind*pi*electron_mass)*&
       (4._fp_kind*pi*lambda_const)**(1._fp_kind/3._fp_kind)*lambda_const
  real(fp_kind), parameter :: dc_const1 = 1._fp_kind/sqrt(2._fp_kind) - 1._fp_kind

  ! These variables obviously must be saved according to the ifcoulomb_old logic below
  ! when noting that coulomb_adjust outputs dscmccon and dscmdcon
  save ifcoulomb_old, dscmccon, dscmdcon, dx, ylo, dylo2, yhi, dyhi, dyhi2

  if(ifcoulomb.eq.2.and.if_dc.eq.1) error stop 'coulomb: ifcoulomb = 2 and if_dc = 1 not allowed'

  ! PTEH, eq. 19
  zeta2 = sum2/sum0
  lambda = lambda_const*zeta2*sqrt(sum2/(t*t*t))
  ! lambda0 = -lambda/sum0
  lambda2 = 1.5_fp_kind*lambda/sum2
  lambdat = -1.5_fp_kind*lambda/t
  ! lambda00 = -2.d0*lambda0/sum0
  ! lambda02 = -lambda2/sum0
  ! lambda0t = -lambdat/sum0
  lambda22 = 0.5_fp_kind*lambda2/sum2
  lambda2t = -1.5_fp_kind*lambda2/t
  lambdatt = -2.5_fp_kind*lambdat/t
  if(if_dc.eq.1) then
     ! I have long since disabled this option because of its known
     ! invalidity for cooler temperatures and higher densities right
     ! where the Coulomb interaction has its strongest effect on EOS
     ! results.  Note, the implementation of this option was a quick
     ! hack before I decided to drop it.  Thus, if you ever want to try
     ! this option again, the implementation should be carefully checked
     ! against the original literature on the diffraction correction.
     ! Also, the implementation should be checked for correct
     ! partial derivatives and thermodynamic consistency.
     error stop 'coulomb: the diffraction correction is disabled.'
     ! delta lambda due to diffraction correction
     !   = dc_const*lambda*(sum2 + dc_const1*x)*x/(T * sum2^1.5d0)
     !   = dc_const*lambda_const*(sum2 + dc_const1*x)/(t^2.5)*(x/sum0)
     dc_lambda = dc_const*(sum2 + dc_const1*x)/(t*t*sqrt(t))*(x/sum0)
     ! by definition of gamma_e
     gamma_e = dc_lambda/lambda/gamma_e_const
     ! dc_lambda0 = -dc_lambda/sum0
     dc_lambda2 = dc_const/(t*t*sqrt(t))*(x/sum0)
     dc_lambdat = -2.5_fp_kind*dc_lambda/t
     dc_lambdax = dc_lambda/x + dc_const*dc_const1/(t*t*sqrt(t))*(x/sum0)
     ! dc_lambda00 = -2.d0*dc_lambda0/sum0
     ! dc_lambda02 = -dc_lambda2/sum0
     ! dc_lambda0t = -2.5d0*dc_lambda0/t
     ! dc_lambda0x = -dc_lambdax/sum0
     dc_lambda22 = 0._fp_kind
     dc_lambda2t = -2.5_fp_kind*dc_lambda2/t
     dc_lambda2x = dc_lambda2/x
     dc_lambdatt = -3.5_fp_kind*dc_lambdat/t
     dc_lambdatx = -2.5_fp_kind*dc_lambdax/t
     dc_lambdaxx = 2._fp_kind*dc_const*dc_const1/(t*t*sqrt(t))/sum0
     lambda = lambda + dc_lambda
     ! lambda0 = lambda0 + dc_lambda0
     lambda2 = lambda2 + dc_lambda2
     lambdat = lambdat + dc_lambdat
     lambdax = dc_lambdax
     ! lambda00 = lambda00 + dc_lambda00
     ! lambda02 = lambda02 + dc_lambda02
     ! lambda0t = lambda0t + dc_lambda0t
     ! lambda0x = dc_lambda0x
     lambda22 = lambda22 + dc_lambda22
     lambda2t = lambda2t + dc_lambda2t
     lambda2x = dc_lambda2x
     lambdatt = lambdatt + dc_lambdatt
     lambdatx = dc_lambdatx
     lambdaxx = dc_lambdaxx
  else
     lambdax = 0._fp_kind
     ! lambda0x = 0.d0
     lambda2x = 0._fp_kind
     lambdatx = 0._fp_kind
     lambdaxx = 0._fp_kind
  endif
  ! just before PTEH eq. 23
  gamma = (lambda*lambda/3._fp_kind)**(1._fp_kind/3._fp_kind)
  xvalue = log(gamma)
  ! For this case must define this to an arbitrary value to avoid
  ! Boolean .and. innocuous but uninitialized logic possibility below
  ! that causes a nagfor warning.
  if(.not.(ifcoulomb.eq.5.or.ifcoulomb.eq.6)) ifdh_limit = .false.
  if(ifcoulomb.ge.3.and.ifcoulomb.le.9) then
     ! reduces to PTEH eq. 25 for original parameters
     ! dh = debye-huckel limit
     dh = gamma*sqrt(gamma)/sqrt(3._fp_kind)
     ddh = 1.5_fp_kind*dh/gamma
     ddh2 = 0.5_fp_kind*ddh/gamma
     if(ifcoulomb.eq.5.or.ifcoulomb.eq.6) then
        if(ifcoulomb.eq.5) then
           ! (For solar conditions, Gamma < 0.367 or log10(gamma) < -0.435
           ! for proper Lambda and Gamma normalization.)
           ! Final attempt at correction to DH so that Coulomb treatment is
           ! essentially DH for solar conditions, and for larger Gammas there
           ! is a rapid, but smooth transition so this results corrects DH to
           ! the OCP result.  In this latter case our mean
           ! Gamma is defined in such a way that we obtain a good
           ! approximation to the multi-component plasma results (see paper).
           xdh = xdh10*ln10
           xmocp = xmocp10*ln10
        elseif(ifcoulomb.eq.6) then
           ! alternative variation to show how sensitive results are
           ! to approximation for intermediate Gamma values.
           xdh = -1._fp_kind*ln10
           xmocp = 0._fp_kind
        endif
        ifdh_limit = xvalue.le.xdh
        if(ifcoulomb.ne.ifcoulomb_old) then
           ifcoulomb_old = ifcoulomb
           ! find modified coefficients dscmccon, dscmdcon so that modified
           ! OCP result is second-order continuous (via a cubic polynomial)
           ! with the DH result.
           call coulomb_adjust(xdh, xmocp,&
                dscacon, dscbcon, dscscon, dscccon, dscdcon,&
                dscmccon, dscmdcon)
           dx = xmocp - xdh
           ylo = 1.5_fp_kind*xdh - 0.5_fp_kind*log(3._fp_kind)
           ! dylo = 1.5d0
           dylo2 = 0._fp_kind
           yhi = dscacon*exp(xmocp) +&
                (dscbcon*exp(dscscon*xmocp) + dscmccon*xmocp +&
                dscmdcon)
           dyhi = dscacon*exp(xmocp) +&
                (dscscon*dscbcon*exp(dscscon*xmocp) +&
                dscmccon)
           dyhi2 = dscacon*exp(xmocp) +&
                dscscon*dscscon*dscbcon*exp(dscscon*xmocp)
           ! transform to ln yhi.
           ! First do commented form using untransformed variables on RHS
           ! yhi = log(yhi)
           ! dyhi = dyhi/yhi
           ! dyhi2 = dyhi2/yhi - dyhi*dyhi/(yhi*yhi)
           dyhi = dyhi/yhi
           dyhi2 = dyhi2/yhi - dyhi*dyhi
           yhi = log(yhi)
        endif
        if(ifdh_limit) then
           ! These old expressions work out to the D-H limit once
           ! the derivatives transformed to Lambda derivatives, but
           ! that transformation has significance loss for dgc2.
           ! gc = dh
           ! dgc = ddh
           ! dgc2 = ddh2
           ! Debye-Huckel limit with derivatives already transformed
           ! to Lambda derivatives to eliminate significance loss.
           gc = lambda/3._fp_kind
           dgc = 1._fp_kind/3._fp_kind
           dgc2 = 0._fp_kind
        elseif(xvalue.lt.xmocp) then
           ! cubic which by design is second-order continuous with
           ! ln(DH) at xdh, and second-order continuous with ln(modified
           ! ocp) at xmocp, where coulomb_adjust above makes sure
           ! of that continuity by adjusting the coefficients that
           ! define the modiefied ocp result.
           ! acap, bcap, ccap, and dcap from numerical recipes 3.3.2, 3.3.4
           acap = (xmocp-xvalue)/dx
           bcap = (xvalue-xdh)/dx
           ! ccap = acap*(acap*acap-1.d0)*dx*dx/6.d0
           ccap = (xmocp-xvalue)*(acap*acap-1._fp_kind)*dx/6._fp_kind
           ! dcap = bcap*(bcap*bcap-1.d0)*dx*dx/6.d0
           dcap = (xvalue-xdh)*(bcap*bcap-1._fp_kind)*dx/6._fp_kind
           ! interpolated cubic, derivative, and second derivative
           ! from numerical recipes 3.3.3, 3.3.5, 3.3.6
           gc = acap*ylo + bcap*yhi + ccap*dylo2 + dcap*dyhi2
           dgc = (yhi-ylo)/dx -&
                (3._fp_kind*acap*acap-1._fp_kind)*dx*dylo2/6._fp_kind +&
                (3._fp_kind*bcap*bcap-1._fp_kind)*dx*dyhi2/6._fp_kind
           dgc2 = acap*dylo2 + bcap*dyhi2
           ! *****convert from ln gc to gc.
           gc = exp(gc)
           ! commented out versions have RHS expressed in untransformed
           ! variables.
           ! dgc2 = exp(gc)*(dgc*dgc + dgc2)
           ! the dgc on the RHS below is untransformed so this must be done
           ! before dgc is transformed below.  The gc on the RHS below
           ! is transformed so must be done after the gc transform above.
           dgc2 = gc*(dgc*dgc + dgc2)
           ! dgc = exp(gc)*dgc
           ! must be done after gc transform above.
           dgc = gc*dgc
           ! *****convert from X = ln gamma independent variable to gamma.
           ! second derivative done fir so can use untransformed
           ! quantities on RHS.
           dgc2 = (dgc2-dgc)/(gamma*gamma)
           dgc = dgc/gamma
        else
           ! MODIFIED one-component plasma (OCP) approximation taken from
           ! DeWitt et al and modified to be second-order continuous
           ! (via a cubic polynomial) with DH result.
           gc = dscacon*gamma +&
                (dscbcon*gamma**dscscon + dscmccon*xvalue +&
                dscmdcon)
           dgc = dscacon +&
                (dscscon*dscbcon*gamma**(-1._fp_kind+dscscon) +&
                dscmccon/gamma)
           dgc2 =&
                ((-1._fp_kind+dscscon)*dscscon*dscbcon*&
                gamma**(-2._fp_kind+dscscon) -&
                dscmccon/(gamma*gamma))
        endif
        ! end of "if(ifcoulomb.eq.5.or.ifcoulomb.eq.6) then" block
     elseif(ifcoulomb.eq.7.or.ifcoulomb.eq.8) then
        ! Use same limit as indicated by extraordinarily good DH solar fit.
        if(log10(gamma).lt.-0.4_fp_kind) then
           dgc2 = ddh2
           dgc = ddh
           gc = dh
        else
           ! one-component plasma (OCP) approximation taken from
           ! DeWitt et al.
           gc = dscacon*gamma +&
                (dscbcon*gamma**dscscon + dscccon*xvalue +&
                dscdcon)
           dgc = dscacon +&
                (dscscon*dscbcon*gamma**(-1._fp_kind+dscscon) +&
                dscccon/gamma)
           dgc2 =&
                ((-1._fp_kind+dscscon)*dscscon*dscbcon*&
                gamma**(-2._fp_kind+dscscon) -&
                dscccon/(gamma*gamma))
        endif
        ! end of "elseif(ifcoulomb.eq.7.or.ifcoulomb.eq.8) then" block
     elseif(ifcoulomb.eq.3.or.ifcoulomb.eq.4.or.ifcoulomb.eq.9)&
          then
        ! ocp = simplified one-component plasma approximation taken from
        ! pteh (works in conjunction with smooth transition from dh to
        ! provide what is actually quite good approximation for ocp
        ! above gamma = 1.
        ocp = acon*gamma + dcon
        docp = acon
        docp2 = 0._fp_kind
        ! smooth transition between dh and ocp
        gc = (dh**gcon + ocp**gcon)**(1._fp_kind/gcon)
        dgc = (ddh*dh**(-1._fp_kind+gcon) + docp*ocp**(-1._fp_kind+gcon))*&
             (dh**gcon + ocp**gcon)**(-1._fp_kind + 1._fp_kind/gcon)
        dgc2 = (ddh2*dh**(-1._fp_kind+gcon) + docp2*ocp**(-1._fp_kind+gcon) +&
             ddh*ddh*(-1._fp_kind+gcon)*dh**(-2._fp_kind+gcon) +&
             docp*docp*(-1._fp_kind+gcon)*ocp**(-2._fp_kind+gcon))*&
             (dh**gcon + ocp**gcon)**(-1._fp_kind + 1._fp_kind/gcon) +&
             (1._fp_kind - gcon)*&
             (ddh*dh**(-1._fp_kind+gcon) + docp*ocp**(-1._fp_kind+gcon))*&
             (ddh*dh**(-1._fp_kind+gcon) + docp*ocp**(-1._fp_kind+gcon))*&
             (dh**gcon + ocp**gcon)**(-2._fp_kind + 1._fp_kind/gcon)
        ! end of "elseif(ifcoulomb.eq.3.or.ifcoulomb.eq.4.or.ifcoulomb.eq.9)&" block
     else
        error stop 'coulomb: logic error'
     endif
     ! convert derivatives to lambda except for case where that has
     ! been done already.

     if(.not.((ifcoulomb.eq.5.or.ifcoulomb.eq.6).and.ifdh_limit)) then
        dgc = dgc*(2._fp_kind/3._fp_kind)*gamma/lambda
        dgc2 = dgc2*(4._fp_kind/9._fp_kind)*(gamma/lambda)*(gamma/lambda) -&
             dgc/(3._fp_kind*lambda)
     endif
     ! end of "if(ifcoulomb.ge.3.and.ifcoulomb.le.9) then" block
  elseif(ifcoulomb.eq.1.or.ifcoulomb.eq.2) then
     ! Debye-Huckel limit
     gc = lambda/3._fp_kind
     dgc = 1._fp_kind/3._fp_kind
     dgc2 = 0._fp_kind
  else
     error stop 'coulomb: ifcoulomb must be in range from 1-9'
  endif
  fcoulomb = -sum0*boltzmann*t*gc
  if(ifcoulomb.eq.1.or.ifcoulomb.eq.2) then
     ! in this case gc proportional lambda proportional to sum0^-1 and
     ! fcoulomb is completely independent of sum0.
     fcoulomb0 = 0._fp_kind
  else
     fcoulomb0 = -boltzmann*t*(gc - dgc*lambda)
  endif
  fcoulomb2 = fcoulomb*dgc/gc*lambda2
  fcoulombt = fcoulomb*(1._fp_kind/t + dgc/gc*lambdat)
  fcoulombx = fcoulomb*dgc/gc*lambdax
  ! use the fact that lambda proportional to sum0^-1
  ! (diffraction corrected or not)
  fcoulomb00 = -boltzmann*t*dgc2*lambda*lambda/sum0
  fcoulomb02 = -boltzmann*t*(-dgc2*lambda2*lambda)
  fcoulomb0t = fcoulomb0/t - boltzmann*t*(-dgc2*lambdat*lambda)
  fcoulomb0x = -boltzmann*t*(-dgc2*lambdax*lambda)
  fcoulomb22 = fcoulomb2*(dgc2*lambda2/dgc + lambda22/lambda2)
  fcoulomb2t = fcoulomb2*&
       (1._fp_kind/t + dgc2*lambdat/dgc + lambda2t/lambda2)
  fcoulomb2x = fcoulomb2*(dgc2*lambdax/dgc + lambda2x/lambda2)
  fcoulombtt = -sum0*boltzmann*(2._fp_kind*dgc*lambdat +&
       t*(dgc2*lambdat*lambdat + dgc*lambdatt))
  fcoulombtx = -sum0*boltzmann*(dgc*lambdax +&
       t*(dgc2*lambdat*lambdax + dgc*lambdatx))
  fcoulombxx = -sum0*boltzmann*t*&
       (dgc2*lambdax*lambdax + dgc*lambdaxx)
  if(ifcoulomb.eq.2) then
     ! prior to this, ifcoulomb = 1 and 2 logic is identical.
     ! with this flag on correct DH results by tau(x) correction
     ! also note that prior x derivatives are zero because if_dc
     ! must be 0 if ifcoulomb = 2.

     allocate(dtau(ntau), d2tau(ntau,ntau))

     if(if_taint_allocated_real) then
        call taint_allocated_real(dtau)
        call taint_allocated_real(d2tau)
     endif

     call tau_calc(sum0, sum2, t, x, tau, dtau, d2tau)
     ! second partials first to use untransformed first partials
     ! and fcoulomb
     fcoulomb00 = fcoulomb*d2tau(1,1)
     fcoulomb02 = fcoulomb2*dtau(1) + fcoulomb*d2tau(1,2)
     fcoulomb0t = fcoulombt*dtau(1) + fcoulomb*d2tau(1,3)
     fcoulomb0x = fcoulomb*d2tau(1,4)
     fcoulomb22 = fcoulomb22*tau + 2._fp_kind*fcoulomb2*dtau(2) +&
          fcoulomb*d2tau(2,2)
     fcoulomb2t = fcoulomb2t*tau + fcoulomb2*dtau(3) +&
          fcoulombt*dtau(2) + fcoulomb*d2tau(2,3)
     fcoulomb2x = fcoulomb2*dtau(4) + fcoulomb*d2tau(2,4)
     fcoulombtt = fcoulombtt*tau + 2._fp_kind*fcoulombt*dtau(3) +&
          fcoulomb*d2tau(3,3)
     fcoulombtx = fcoulombt*dtau(4) + fcoulomb*d2tau(3,4)
     fcoulombxx = fcoulomb*d2tau(4,4)
     ! first partials next to use untransformed fcoulomb
     ! original fcoulomb0 and fcoulombx are zero in DH limit.
     fcoulomb0 = fcoulomb*dtau(1)
     fcoulomb2 = fcoulomb2*tau + fcoulomb*dtau(2)
     fcoulombt = fcoulombt*tau + fcoulomb*dtau(3)
     fcoulombx = fcoulomb*dtau(4)
     fcoulomb = fcoulomb*tau
  endif
end subroutine coulomb
