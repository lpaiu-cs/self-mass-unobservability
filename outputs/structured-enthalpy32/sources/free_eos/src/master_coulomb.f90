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

! input quantities:
! ifnr_in = 0, simple iteration, calculate dve, dv0, dv2, and all f and
!   t derivatives with no other variable fixed.
! ifnr_in = 1, NR iteration, calculate dve, dv0, dv2, and all sum0 and
!   sum2 derivatives.
! ifnr_in = 2, calculate dve, dv0, dv2, and all f and
!   t derivatives *keeping* sum0 and sum2 fixed.
! ifnr_in = 3, is combination of ifnr_in = 1 and 2.
! rhostar(9), 1=re proportional to n_e,
! the rest are partials of ln n_e wrt fl and tl,
! 2 = f, 3 = t, 4 = ff, 5 = ft, 6 = tt, 7 = fff, 8 = fft 9 = ftt
! sum0 = sum_ions nion, where nion is the number per unit volume,
!   and the sum taken over positive ions, but excluding n_e.
! sum2 = sum_ions Z^2 nion,
!   Z is the charge on the ions.
!   outside this routine, sum0 and sum2 are calculated as a function
!     of n_e *or* f and T, but for programming convenience sum0
!     and sum2 are viewed as functions of n_e, f, and T.
!     The derivatives calculated outside are:
!       sum0ne = partial sum0(n_e, f, T) wrt n_e
!       sum0f = partial sum0(n_e, f, T) wrt fl
!       sum0t = partial sum0(n_e, f, T) wrt tl
!       sum2ne = partial sum2(n_e, f, T) wrt n_e
!       sum2f = partial sum2(n_e, f, T) wrt fl
!       sum2t = partial sum2(n_e, f, T) wrt tl
!       Since the actual functional dependence is n_e *or* f, T, it is the
!         calling routine's responsibility to calculate one set of
!         derivatives and set the the other set to zero as appropriate.
! n_e = electron number density per unit volume
! t = temperature (Kelvin)
! p_e = electron pressure (cgs)
! pstar(9), 1=proportional to p_e,
! the rest are partials of ln p_e wrt fl and tl,
! 2 = f, 3 = t, 41 = ff, 5 = ft, 6 = tt, 7 = fff, 8 = fft 9 = ftt
! ifcoulomb = 0 ignore the effects of the Coulomb interaction.
! ifcoulomb = 1 treat Coulomb interaction in the Debye-Huckel
!   approximation.
! ifcoulomb = 2 treat Coulomb interaction in the Debye-Huckel
!   approximation corrected by tau(x).
! ifcoulomb = 3 treat Coulomb interaction in the PTEH approximation
!   with pteh theta_e.
! ifcoulomb = 4 treat Coulomb interaction in the PTEH approximation
!   with fermi-dirac theta_e.
! ifcoulomb = 9 treat Coulomb interaction in the PTEH approximation
!   with fermi-dirac theta_e and DeWitt definition of lambda
!   (using sum0a = sum0 + ne*theta_e).
! ifcoulomb = 5 DH smoothly connected to modified OCP and DeWitt
!   definition of lambda.
! ifcouloumb = 6 same as 5 with alternative smooth connection
! ifcoulomb = 7, DH (Gamma < 1) or OCP using new DeWitt lambda.
! ifcoulomb = 8, same as 7 with theta_e = 0.
! if_dc = 0, ignore effects of the diffraction correction
! if_dc = 1, apply diffraction correction to Lambda
! if_pteh = 0, sum0 and sum2 are functions of n_i only
! if_pteh = 1, sum0 and sum2 are functions of n_e only (PTEH
!   approximation)
! output quantities:
! lambda = plasma interaction parameter (defined by PTEH *or* defined
!   by DeWitt if ifcoulomb > 4 and consistent
!   with the simpler CG 15.60 criterion for single component plasmas.)
!   n.b., if_dc = 1 means this is the diffraction-corrected quantity.
! gamma_e =  DeWitt paper II gamma_e = diffraction correction parameter.
! dve, dvef, dvet, dve0, dve2:
! dv0, dv0f, dv0t, dv00, dv02:
! dv2, dv2f, dv2t, dv22:
!   N.B. the f, t, 0, and 2 suffixes for dve, dv0, and dv2 type of
!   output variables refer to fl, tl, sum0, and sum2 partial derivatives.
!   dve is the change in the electron component of the equilibrium constant.
!   dv0 is the change in the sum0 component of the equilibrium constant.
!   dv2 is the change in the sum2 component of the equilibrium constant.

!> This master_coulomb subroutine calculates quantities required to
!> help determine the non-ideal Coulomb component of equilibrium
!> constants that are iteratively used to help converge the EOS
!> solution.  These calculated results include intermediate quantities
!> (saved in the mod_master_coulomb_data module) that are also needed
!> for other purposes (e.g., calculation of the non-ideal component of
!> the pressure due to the Coulomb effect).
!>
!> \param[in] ifnr_in PARAMETERS NEED DOCUMENTATION
!>
subroutine master_coulomb(ifnr_in, rhostar, f,&
     sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,&
     n_e, t, p_e, pstar, lambda, gamma_e,&
     ifcoulomb, if_dc, if_pteh,&
     dve, dvef, dvet, dve0, dve2,&
     dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
  use mod_free_eos_constants, only: boltzmann, cd
  use mod_master_coulomb_data, only:&
       dmucoulomb0, dmucoulomb2, dmucoulombf, dmucoulombt, dmucoulomb,&
       fcft0f, fcft0t, fcft2f, fcft2t, fcftxf, fcftxt,&
       fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0,&
       fcoulomb22, fcoulomb2t, fcoulomb2,&
       fcoulombtt, fcoulombtx, fcoulombt, fcoulombx, fcoulomb,&
       ifstart,&
       sum0af, sum0ane, sum0atl,&
       sum2af, sum2ane, sum2atl,&
       theta_eff, theta_eft, theta_ef, thetatf, thetatt,&
       xf, xt, x
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) :: ifnr_in, ifcoulomb, if_dc, if_pteh
  real(fp_kind), intent(in) :: rhostar(:), f,&
       sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,&
       n_e, t, p_e, pstar(:)
  real(fp_kind), intent(out) :: lambda, gamma_e,&
       dve, dvef, dvet, dve0, dve2,&
       dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22

  ! Local variables
  integer, parameter :: nstar = 9
  integer ifnr
  real(fp_kind) dpsidf, dpsidff, dpsidfff
  real(fp_kind) theta_e, theta_et, theta_ett,&
       thetane, thetanef, thetanet, thetat, theta_pow,&
       sum0a, sum0anef, sum0anet,&
       sum2a, sum2anef, sum2anet,&
       dvcoulomb, dvcoulombf, dvcoulombt,&
       fcoulomb0x, fcoulomb2x, fcoulombxx,&
       thetax, xe, xef, xet
  logical ifnr13, ifnr23

  ! Sanity check
  if(nstar.ne.size(rhostar).or.nstar.ne.size(pstar))&
       error stop 'master_coulomb: bad sizes for rhostar or pstar'

  ! Mark this routine as called.
  if(ifstart.eq.1) ifstart = 0

  if(ifcoulomb.lt.0.or.ifcoulomb.gt.9) then
     write(stderr,*) 'master_coulomb: ifcoulomb = ', ifcoulomb
     error stop 'master_coulomb: invalid ifcoulomb'
  elseif(ifcoulomb.eq.0) then
     ! For ifcount.eq.0, zero all intent(out) results.
     lambda = 0._fp_kind
     gamma_e = 0._fp_kind
     dve = 0._fp_kind
     dvef = 0._fp_kind
     dvet = 0._fp_kind
     dve0 = 0._fp_kind
     dve2 = 0._fp_kind
     dv0 = 0._fp_kind
     dv0f = 0._fp_kind
     dv0t = 0._fp_kind
     dv2 = 0._fp_kind
     dv2f = 0._fp_kind
     dv2t = 0._fp_kind
     dv00 = 0._fp_kind
     dv02 = 0._fp_kind
     dv22 = 0._fp_kind
  else
     if(ifnr_in.gt.3) error stop 'master_coulomb: bad ifnr_in'

     ! must zero NR derivative calculations when if_pteh = 1
     if(if_pteh.eq.1) then
        ifnr = min(ifnr_in,0)
        dve0 = 0._fp_kind
        dve2 = 0._fp_kind
        dv00 = 0._fp_kind
        dv02 = 0._fp_kind
        dv22 = 0._fp_kind
     else
        ifnr = ifnr_in
     endif
     ifnr13 = ifnr.eq.1.or.ifnr.eq.3
     ifnr23 = ifnr.eq.2.or.ifnr.eq.3
     ! don't want tau(x) correction and diffraction correction
     ! simultaneously
     if(ifcoulomb.eq.2.and.if_dc.eq.1)&
          error stop 'master_coulomb: bad combination of ifcoulomb and if_dc'
     ! set up coulomb interaction treatment
     ! partial psi/d ln f and higher derivatives
     dpsidf = sqrt(1._fp_kind+f)
     dpsidff = 0.5_fp_kind*f*dpsidf/(1._fp_kind+f)
     dpsidfff = dpsidff*(1._fp_kind - 0.5_fp_kind*f/(1._fp_kind+f))
     ! Decide power of original theta_e
     if(ifcoulomb.gt.5) then
        ! Increase this to get more rapid quenching of electronic
        ! component as the degeneracy is increased.
        ! Currently experimental, but at some point this will also
        ! be used for standard Coulomb treament (ifcoulomb = 5).
        ! theta_pow = 1.3d0
        theta_pow = 1._fp_kind
     else
        ! Traditional value for MDH, etc.
        theta_pow = 1._fp_kind
     endif
     ! partial ln n_e(T, psi)/partial psi
     if(ifcoulomb.ne.3.and.ifcoulomb.ne.8) then
        theta_e = rhostar(2)/dpsidf
        theta_ef = (rhostar(4)-rhostar(2)*dpsidff/dpsidf)/dpsidf
        theta_et = rhostar(5)/dpsidf
        theta_eff = (rhostar(7)-(2._fp_kind*rhostar(4)*dpsidff +&
             rhostar(2)*(dpsidfff-2._fp_kind*dpsidff/dpsidf*dpsidff))/&
             dpsidf)/dpsidf
        theta_eft = (rhostar(8)-rhostar(5)*dpsidff/dpsidf)/dpsidf
        theta_ett = rhostar(9)/dpsidf
        if(theta_pow.ne.1._fp_kind) then
           ! transform to power of original theta_e
           theta_eff =&
                theta_eff*theta_pow*theta_e**(theta_pow-1._fp_kind) +&
                theta_ef*theta_ef*theta_pow*(theta_pow-1._fp_kind)*&
                theta_e**(theta_pow-2._fp_kind)
           theta_eft =&
                theta_eft*theta_pow*theta_e**(theta_pow-1._fp_kind) +&
                theta_ef*theta_et*theta_pow*(theta_pow-1._fp_kind)*&
                theta_e**(theta_pow-2._fp_kind)
           theta_ett =&
                theta_ett*theta_pow*theta_e**(theta_pow-1._fp_kind) +&
                theta_et*theta_et*theta_pow*(theta_pow-1._fp_kind)*&
                theta_e**(theta_pow-2._fp_kind)
           theta_ef = theta_ef*theta_pow*theta_e**(theta_pow-1._fp_kind)
           theta_et = theta_et*theta_pow*theta_e**(theta_pow-1._fp_kind)
           theta_e = theta_e**theta_pow
        endif
        ! n.b. partial theta_e(f(ne,T), T)/partial ln ne =
        ! theta_ef/rhostar(2)
        thetane = theta_ef/rhostar(2)
        thetanef = (theta_eff - theta_ef*rhostar(4)/&
             rhostar(2))/rhostar(2)
        thetanet = (theta_eft - theta_ef*rhostar(5)/&
             rhostar(2))/rhostar(2)
        ! n.b. partial theta_e(ne, T)/partial ln T =
        ! theta_et - theta_ef*rhostar(3)/rhostar(2)
        thetat = theta_et - theta_ef*rhostar(3)/rhostar(2)
        thetatf = theta_eft - theta_eff*rhostar(3)/rhostar(2) -&
             theta_ef*(rhostar(5) - rhostar(3)*rhostar(4)/rhostar(2))/&
             rhostar(2)
        thetatt = theta_ett - theta_eft*rhostar(3)/rhostar(2) -&
             theta_ef*(rhostar(6) - rhostar(3)*rhostar(5)/rhostar(2))/&
             rhostar(2)
     elseif(ifcoulomb.eq.3) then
        ! use simplified pteh approximation for theta_e and
        ! all subsequent derivatives based on their eq. 29.
        call pteh_theta(cd*rhostar(1), t, rhostar,&
             theta_e, theta_ef, theta_et,&
             theta_eff, theta_eft, theta_ett,&
             thetane, thetanef, thetanet,&
             thetat, thetatf, thetatt)
     elseif(ifcoulomb.eq.8) then
        theta_e = 0._fp_kind
        theta_ef = 0._fp_kind
        theta_et = 0._fp_kind
        theta_eff = 0._fp_kind
        theta_eft = 0._fp_kind
        theta_ett = 0._fp_kind
        thetane = 0._fp_kind
        thetanef = 0._fp_kind
        thetanet = 0._fp_kind
        thetat = 0._fp_kind
        thetatf = 0._fp_kind
        thetatf = 0._fp_kind
     endif
     if(ifcoulomb.gt.4) then
        ! DeWitt definition
        sum0a = sum0 + theta_e*n_e
        sum0af = theta_ef*n_e + theta_e*n_e*rhostar(2)
     else
        ! PTEH definition
        sum0a = sum0
        sum0af = 0._fp_kind
     endif
     sum2a = sum2 + theta_e*n_e
     sum2af = theta_ef*n_e + theta_e*n_e*rhostar(2)
     ! sum0af and sum2af are the partials of sum?a wrt fl for fixed
     ! auxiliary variables.  Ordinarily sum0 and sum2 *are* auxiliary
     ! variables (or are transformed later to a different set of
     ! auxiliary variables) and we are done.  However, for special case
     ! of if_pteh = 1, sum0 and sum2 are not auxiliary variables, and
     ! we must do more.
     if(if_pteh.eq.1) then
        sum0af = sum0af + sum0ne*rhostar(2)
        sum2af = sum2af + sum2ne*rhostar(2)
     endif
     ! n.b. x and its derivatives and the partials of fcoulomb
     ! wrt x returned by coulomb below *only* relevant when
     ! if_dc.eq.1 or ifcoulomb.eq.2.
     ! n.b.
     ! xf = partial x(f,T)/partial fl
     ! xt = partial x(f,T)/partial ft
     ! xe = partial x(n_e,T)/partial n_e
     ! xef = partial xe(f,T)/partial fl
     ! xet = partial xe(f,T)/partial tl
     if(if_dc.eq.1) then
        x = theta_e*n_e
        xf = theta_ef*n_e + x*rhostar(2)
        xt = theta_et*n_e + x*rhostar(3)
        xe = theta_e + thetane
        xef = theta_ef + thetanef
        xet = theta_et + thetanet
     elseif(ifcoulomb.eq.2) then
        ! thetax = (n_e k T/p_e) ~ 1
        thetax = (n_e*boltzmann*t/p_e)
        x = thetax*n_e
        xf = x*(2._fp_kind*rhostar(2) - pstar(2))
        xt = x*(2._fp_kind*rhostar(3) - pstar(3) + 1._fp_kind)
        xe = thetax*(2._fp_kind - pstar(2)/rhostar(2))
        xef = thetax*((2._fp_kind - pstar(2)/rhostar(2))*&
             (rhostar(2)-pstar(2)) +&
             (-pstar(4) + pstar(2)*rhostar(4)/rhostar(2))/rhostar(2))
        xet = thetax*((2._fp_kind - pstar(2)/rhostar(2))*&
             (rhostar(3) - pstar(3) + 1._fp_kind) +&
             (-pstar(5) + pstar(2)*rhostar(5)/rhostar(2))/rhostar(2))
     else
        ! should define since used below (although x derivatives
        ! set to zero by coulomb).
        x = 0._fp_kind
        xf = 0._fp_kind
        xt = 0._fp_kind
        xe = 0._fp_kind
        xef = 0._fp_kind
        xet = 0._fp_kind
     endif
     call coulomb(ifcoulomb, if_dc, sum0a, sum2a, t, x,&
          fcoulomb, fcoulomb0, fcoulomb2, fcoulombt, fcoulombx,&
          fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0x,&
          fcoulomb22, fcoulomb2t, fcoulomb2x,&
          fcoulombtt, fcoulombtx, fcoulombxx, lambda, gamma_e)
     ! fcoulomb (free energy per unit volume) is a function of all
     ! n_i (where i ranges over positive ions), n_e, T.
     ! Calculate additions to chemical potential of electrons,
     ! partial (V*fcoulomb(n_i, n_e, T)) wrt n_e * partial n_e wrt (V n_e)
     !   = partial fcoulomb(sum0a(n_i,n_e), sum2a(n_i,n_e), n_e, T) wrt n_e
     ! n.b. sum0 and sum2 are functions of n_e alone or are
     !   independent of n_e.  (depending on ifcoulomb for
     !   sum0) sum0a and sum2a have the extra factor
     !   theta_e(f(n_e,T),T)*n_e
     if(ifcoulomb.gt.4) then
        ! DeWitt definition
        ! sum0a = sum0 + theta_e*n_e
        ! calculate partial sum0a(n_e, n_i, T) wrt n_e
        ! n.b., when n_i is fixed, sum0 is fixed except in pteh
        ! approximation when it varies with n_e.
        sum0ane = sum0ne + theta_e + thetane
        ! calculate partial sum0a(n_e, n_i, T) wrt ln T
        ! (ignore sum0 dependence on t, since fixed n_i means sum0 is fixed
        ! or depends on fixed n_e in initial approximation)
        sum0atl = n_e*thetat
     else
        ! PTEH definition
        ! sum0a = sum0
        sum0ane = sum0ne
        sum0atl = 0._fp_kind
     endif
     ! calculate partial sum2a(n_e, n_i, T) wrt n_e
     ! n.b., when n_i is fixed, sum2 is fixed except in pteh
     ! approximation when it varies with n_e.
     sum2ane = sum2ne + theta_e + thetane
     ! calculate partial sum2a(n_e, n_i, T) wrt ln T
     ! (ignore sum2 dependence on t, since fixed n_i means sum2 is fixed
     ! or depends on fixed n_e in initial approximation)
     sum2atl = n_e*thetat
     ! partial fcoulomb(sum0, sum2, t, ne) wrt ne
     dmucoulomb = fcoulomb0*sum0ane + fcoulomb2*sum2ane +&
          fcoulombx*xe
     ! partials of dmucoulomb wrt sum0 and sum2
     ! note that sum[02]ane and xe are completely independent of sum[02].
     dmucoulomb0 = fcoulomb00*sum0ane + fcoulomb02*sum2ane +&
          fcoulomb0x*xe
     dmucoulomb2 = fcoulomb02*sum0ane + fcoulomb22*sum2ane +&
          fcoulomb2x*xe
     ! change in equilibrium constant related to negative chemical
     ! potential/kT.
     dvcoulomb = -dmucoulomb/(boltzmann*t)
     if(ifnr.eq.0.or.ifnr23) then
        ! following partial derivatives of fcoulomb0, fcoulomb2, and
        ! fcoulombx are wrt ln f, ln t holding nothing else fixed
        ! if(ifnr.eq.0) or holding input sum0 and sum2 fixed (ifnr23)
        ! f derivatives first....
        fcft0f =&
             n_e*rhostar(2)*(fcoulomb00*sum0ane +&
             fcoulomb02*sum2ane) + fcoulomb0x*xe*n_e*rhostar(2)
        fcft2f =&
             n_e*rhostar(2)*(fcoulomb02*sum0ane +&
             fcoulomb22*sum2ane) + fcoulomb2x*xe*n_e*rhostar(2)
        fcftxf =&
             n_e*rhostar(2)*(fcoulomb0x*sum0ane +&
             fcoulomb2x*sum2ane) + fcoulombxx*xe*n_e*rhostar(2)
        ! we recognize that mixed second partial derivatives of (input)
        ! sum0 and sum2 wrt ne, f = 0.
        if(ifcoulomb.gt.4) then
           ! DeWitt definition
           ! sum0a = sum0 + theta_e*n_e
           sum0anef = theta_ef + thetanef
        else
           ! PTEH definition
           ! sum0a = sum0
           sum0anef = 0._fp_kind
        endif
        sum2anef = theta_ef + thetanef
        ! now t derivatives....
        fcft0t =&
             t*fcoulomb0t +&
             fcoulomb00*sum0atl + fcoulomb02*sum2atl +&
             n_e*rhostar(3)*(fcoulomb00*sum0ane + fcoulomb02*sum2ane) +&
             fcoulomb0x*xt
        fcft2t =&
             t*fcoulomb2t +&
             fcoulomb02*sum0atl + fcoulomb22*sum2atl +&
             n_e*rhostar(3)*(fcoulomb02*sum0ane + fcoulomb22*sum2ane) +&
             fcoulomb2x*xt
        fcftxt =&
             t*fcoulombtx +&
             fcoulomb0x*sum0atl + fcoulomb2x*sum2atl +&
             n_e*rhostar(3)*(fcoulomb0x*sum0ane + fcoulomb2x*sum2ane) +&
             fcoulombxx*xt
        ! we recognize that the second mixed partial derivatives of
        ! (input) sum0 and sum2 wrt ne, t = 0.
        if(ifcoulomb.gt.4) then
           ! DeWitt definition
           ! sum0a = sum0 + theta_e*n_e
           sum0anet = theta_et + thetanet
        else
           ! PTEH definition
           ! sum0a = sum0
           sum0anet = 0._fp_kind
        endif
        sum2anet = theta_et + thetanet
        if(ifnr.eq.0) then
           ! account for sum0, sum2 dependence on f, t
           fcft0f = fcft0f + fcoulomb00*sum0f + fcoulomb02*sum2f
           fcft2f = fcft2f + fcoulomb02*sum0f + fcoulomb22*sum2f
           fcftxf = fcftxf + fcoulomb0x*sum0f + fcoulomb2x*sum2f
           fcft0t = fcft0t + fcoulomb00*sum0t + fcoulomb02*sum2t
           fcft2t = fcft2t + fcoulomb02*sum0t + fcoulomb22*sum2t
           fcftxt = fcftxt + fcoulomb0x*sum0t + fcoulomb2x*sum2t
        endif
        ! dmucoulomb = fcoulomb0*sum0ane + fcoulomb2*sum2ane + fcoulombx*xe
        dmucoulombf = fcft0f*sum0ane + fcft2f*sum2ane +&
             fcoulomb0*sum0anef + fcoulomb2*sum2anef +&
             fcftxf*xe + fcoulombx*xef
        dvcoulombf = -dmucoulombf/(boltzmann*t)
        dmucoulombt = fcft0t*sum0ane + fcft2t*sum2ane +&
             fcoulomb0*sum0anet + fcoulomb2*sum2anet +&
             fcftxt*xe + fcoulombx*xet
        dvcoulombt = (dmucoulomb - dmucoulombt)/(boltzmann*t)
     endif
     ! calculate negative of electron component (dvcoulomb) appropriately....
     dve = dvcoulomb
     if(ifnr.eq.0) then
        ! derivative wrt f, t holding no other variable fixed.
        dvef = dvcoulombf
        dvet = dvcoulombt
     elseif(ifnr23) then
        ! derivative wrt f, t holding sum0 and sum2 fixed.
        dvef = dvcoulombf
        dvet = dvcoulombt
     endif
     if(ifnr13) then
        ! derivative wrt sum0, sum2, holding f, t fixed.
        dve0 = -dmucoulomb0/(boltzmann*t)
        dve2 = -dmucoulomb2/(boltzmann*t)
     endif
     ! remaining equilibrium constant changes are absolute,
     ! *not* relative to next lower ion.
     if(if_pteh.eq.1) then
        ! in this approximation, sum0 and sum2 are independent
        ! of n_i.
        dv0 = 0._fp_kind
        dv0f = 0._fp_kind
        dv0t = 0._fp_kind
        dv2 = 0._fp_kind
        dv2f = 0._fp_kind
        dv2t = 0._fp_kind
     else
        ! calculate dv0 = negative of chemical potential/kT due to
        ! sum0 term = -partial f(sum0, sum2, T, x(f,T))/partial sum0 * partial sum0/partial n_i /kT
        dv0 = -fcoulomb0/(boltzmann*t)
        ! calculate dv2 = negative of chemical potential/kT/Z^2 due to
        ! sum2 term = -partial f(sum0, sum2, T, x(f,T))/partial sum2 * partial sum2/partial n_i /kT
        dv2 = -fcoulomb2/(boltzmann*t)
        if(ifnr.eq.0) then
           ! derivative wrt f, t holding no other variable fixed.
           dv0f = -fcft0f/(boltzmann*t)
           dv0t = -fcft0t/(boltzmann*t) -dv0
           dv2f = -fcft2f/(boltzmann*t)
           dv2t = -fcft2t/(boltzmann*t) -dv2
        elseif(ifnr23) then
           ! derivative wrt f, t holding sum0 and sum2 fixed.
           dv0f = -fcft0f/(boltzmann*t)
           dv0t = -fcft0t/(boltzmann*t) -dv0
           dv2f = -fcft2f/(boltzmann*t)
           dv2t = -fcft2t/(boltzmann*t) -dv2
        endif
        if(ifnr13) then
           ! derivative wrt sum0, sum2, holding f, t fixed.
           dv00 = -fcoulomb00/(boltzmann*t)
           dv02 = -fcoulomb02/(boltzmann*t)
           dv22 = -fcoulomb22/(boltzmann*t)
        endif
     endif
  endif
end subroutine master_coulomb

!> This master_coulomb_pressure subroutine calculates the change to
!> the pressure due to the Coulomb effect and the partial derivatives
!> of that change wrt ln f, ln T, sum0 and sum2.
!>
!> \param[in] ifcoulomb PARAMETERS NEED DOCUMENTATION
!>
subroutine master_coulomb_pressure(&
       ifcoulomb, if_dc, if_pteh, sum0, sum2,&
       t, n_e, rhostar, pstar,&
       dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2)
  use mod_master_coulomb_data, only:&
       dmucoulomb0, dmucoulomb2, dmucoulombf, dmucoulombt, dmucoulomb,&
       fcft0f, fcft0t, fcft2f, fcft2t,&
       fcoulomb00, fcoulomb02, fcoulomb0,&
       fcoulomb22, fcoulomb2,&
       fcoulombt, fcoulombx, fcoulomb,&
       ifstart,&
       sum0atl,&
       sum2atl,&
       theta_ef,&
       x
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) :: ifcoulomb, if_dc, if_pteh
  real(fp_kind), intent(in) :: sum0, sum2, t, n_e, rhostar(:), pstar(:)
  real(fp_kind), intent(out) :: dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2

  ! Internal variables
  real(fp_kind) dscoulomb_local, xfixedt

  ! Sanity check
  if(ifstart.eq.1) error stop 'master_coulomb_pressure: master_coulomb must be called first'

  if(ifcoulomb.lt.0.or.ifcoulomb.gt.9) then
     write(stderr,*) 'master_coulomb_pressure: ifcoulomb = ', ifcoulomb
     error stop 'master_coulomb_pressure: invalid ifcoulomb'
  elseif(ifcoulomb.eq.0) then
     ! For ifcount.eq.0, zero all intent(out) results.
     dpcoulomb = 0._fp_kind
     dpcoulombf = 0._fp_kind
     dpcoulombt = 0._fp_kind
     dpcoulomb0 = 0._fp_kind
     dpcoulomb2 = 0._fp_kind
  else
     ! don't want tau(x) correction and diffraction correction
     ! simultaneously
     if(ifcoulomb.eq.2.and.if_dc.eq.1)&
          error stop 'master_coulomb_pressure: bad ifcoulomb or if_dc'
     ! calculate change to thermodynamic quantities
     ! n.b. fcoulomb is per unit volume so is dscoulomb, ducoulomb
     ! n.b. fcoulomb is considered to be a function of n_e, n_i, T
     ! because sum0a(n_e, n_i), sum2a(n_e, n_i, T), and x(n_e, T),
     ! where x = thetax*n_e (ifcoulomb=2) or theta_e*n_e (if_dc = 1)
     ! partial -fcoulomb(n_i, n_e, t) wrt t
     ! n.b.
     ! xfixedt = partial x(ne,T)/partial T
     ! also note that
     ! partial fl(ne, T)/partial ln T = -rhostar(3)/rhostar(2)
     if(if_dc.eq.1) then
        xfixedt = (-theta_ef*rhostar(3)/rhostar(2))*n_e/t
     elseif(ifcoulomb.eq.2) then
        xfixedt = x*(1._fp_kind + pstar(2)*rhostar(3)/rhostar(2) -&
             pstar(3))/t
     else
        xfixedt = 0._fp_kind
     endif
     ! non-returned version needed for this subroutine
     dscoulomb_local = -fcoulombt -&
          fcoulomb0*sum0atl/t - fcoulomb2*sum2atl/t -&
          fcoulombx*xfixedt
     ! dpcoulomb = - partial(V*fcoulomb(n_e, n_i, T)) wrt V
     !   = -fcoulomb + n_e*partial fcoulomb wrt n_e + sum_i n_i* partial fcoulomb wrt n_i
     !   = -fcoulomb + n_e*dmucoulomb + sum0*partial fcoulomb
     !     wrt sum0a + sum2*partial fcoulomb wrt sum2a
     !   n.b. the latter two terms are dropped if sum0 and sum2
     !   depend only on n_e.  This occurs if
     !   the PTEH approximation being used for sum0, sum2
     if(if_pteh.eq.1) then
        dpcoulomb = -fcoulomb + n_e*dmucoulomb
        ! dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
        ! fcoulomb (for the fully ionized case) is considered
        ! to be a function of (sum0(n_e), sum2(n_e)+
        ! theta_e(n_e,t)*n_e, t, x(n_e,T)) = function(n_e, T),
        ! completely independent of n_i.
        ! also recall dmucoulomb = partial fcoulomb(n_e, T) wrt n_e
        ! in the fully ionized case.
        ! There is a useful cancellation of terms when using
        ! this viewpoint....
        dpcoulombf = n_e*dmucoulombf
        ! dscoulomb = partial -fcoulomb(n_e, t) wrt t
        dpcoulombt = t*dscoulomb_local + n_e*dmucoulombt
        dpcoulomb0 = 0._fp_kind
        dpcoulomb2 = 0._fp_kind
     else
        dpcoulomb = -fcoulomb + n_e*dmucoulomb +&
             sum0*fcoulomb0 + sum2*fcoulomb2
        ! dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
        ! fcoulomb is considered a function of (sum0(f,t), sum2(f,t)+
        ! theta_e(n_e,t)*n_e, t, x(n_e,T)).  In other words, the previous
        ! n_i dependence is replaced with sum0 and sum2 dependence.
        ! There is a useful cancellation of terms when using
        ! this viewpoint....
        dpcoulombf = n_e*dmucoulombf +&
             sum0*fcft0f + sum2*fcft2f
        ! dscoulomb = partial -fcoulomb(sum0, sum2, n_e, t) wrt t
        dpcoulombt = t*dscoulomb_local + n_e*dmucoulombt +&
             sum0*fcft0t + sum2*fcft2t
        dpcoulomb0 = n_e*dmucoulomb0 +&
             sum0*fcoulomb00 + sum2*fcoulomb02
        dpcoulomb2 = n_e*dmucoulomb2 +&
             sum0*fcoulomb02 + sum2*fcoulomb22
     endif
  endif
end subroutine master_coulomb_pressure

!> This master_coulomb_free subroutine calculates the change to the
!> free energy per unit volume due to the Coulomb effect and the
!> partial derivatives of that change wrt fl, sum0, and sum2.
!>
!> \param[in] ifcoulomb PARAMETERS NEED DOCUMENTATION
!>
subroutine master_coulomb_free(ifcoulomb, if_pteh, free_coulomb, free_coulombf, free_coulomb0, free_coulomb2)
  use mod_master_coulomb_data, only:&
       fcoulomb0,&
       fcoulomb2,&
       fcoulombx, fcoulomb,&
       ifstart,&
       sum0af,&
       sum2af,&
       xf
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) :: ifcoulomb, if_pteh
  real(fp_kind), intent(out) :: free_coulomb, free_coulombf, free_coulomb0, free_coulomb2

  ! Sanity check
  if(ifstart.eq.1) error stop 'master_coulomb_free: master_coulomb must be called first'

  if(ifcoulomb.lt.0.or.ifcoulomb.gt.9) then
     write(stderr,*) 'master_coulomb_free: ifcoulomb = ', ifcoulomb
     error stop 'master_coulomb_free: invalid ifcoulomb'
  elseif(ifcoulomb.eq.0) then
     ! For ifcount.eq.0, zero all intent(out) results.
     free_coulomb = 0._fp_kind
     free_coulombf = 0._fp_kind
     free_coulomb0 = 0._fp_kind
     free_coulomb2 = 0._fp_kind
  else
     free_coulomb = fcoulomb
     free_coulombf = fcoulomb0*sum0af + fcoulomb2*sum2af +&
          fcoulombx*xf
     ! n.b. The fcoulomb0 and fcoulomb2 derivatives are actually with
     ! respect to sum0a and sum2a with x = thetax*n_e fixed (if there is a
     ! tau correction), but with f, and t fixed, thetax*n_e is fixed
     ! and also the partials of sum0a with respect to sum0 and sum2a wrt
     ! sum2 are both unity.
     if(if_pteh.eq.1) then
        ! if the pteh Coulomb approximation is used, then sum0 and
        ! sum2 are just functions of f and t.
        free_coulomb0 = 0._fp_kind
        free_coulomb2 = 0._fp_kind
     else
        free_coulomb0 = fcoulomb0
        free_coulomb2 = fcoulomb2
     endif
  endif
end subroutine master_coulomb_free

!> This master_coulomb_end subroutine calculates the change to the
!> pressure, entropy per unit volume, and energy per unit volume due
!> to the Coulomb effect as well as partial derivatives of those
!> those first two quantities wrt fl and tl.
!>
!> \param[in] rhostar PARAMETERS NEED DOCUMENTATION
!>
subroutine master_coulomb_end(rhostar,&
       sum0, sum0ne, sum0f, sum0t,&
       sum2, sum2ne, sum2f, sum2t,&
       n_e, t, pstar,&
       ifcoulomb, if_dc, if_pteh,&
       dpcoulomb, dpcoulombf, dpcoulombt,&
       dscoulomb, dscoulombf, dscoulombt, ducoulomb)
  use mod_master_coulomb_data, only:&
       dmucoulombf, dmucoulombt, dmucoulomb,&
       fcft0f, fcft0t, fcft2f, fcft2t, fcftxf, fcftxt,&
       fcoulomb0t, fcoulomb0,&
       fcoulomb2t, fcoulomb2,&
       fcoulombtt, fcoulombtx, fcoulombt, fcoulombx, fcoulomb,&
       ifstart,&
       sum0ane, sum0atl,&
       sum2ane, sum2atl,&
       theta_eff, theta_eft, theta_ef, thetatf, thetatt,&
       xf, xt, x
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) :: ifcoulomb, if_dc, if_pteh
  real(fp_kind), intent(in) :: rhostar(:),&
       sum0, sum0ne, sum0f, sum0t,&
       sum2, sum2ne, sum2f, sum2t,&
       n_e, t, pstar(:)
  real(fp_kind), intent(out) ::&
       dpcoulomb, dpcoulombf, dpcoulombt,&
       dscoulomb, dscoulombf, dscoulombt, ducoulomb

  ! Local variables
  real(fp_kind) fcfttf, fcfttt, sum0atlf, sum0atlt, sum2atlf, sum2atlt,&
       xfixedt, xfixedtf, xfixedtt

  ! Sanity check
  if(ifstart.eq.1) error stop 'master_coulomb_end: master_coulomb must be called first'

  if(ifcoulomb.lt.0.or.ifcoulomb.gt.9) then
     write(stderr,*) 'master_coulomb_end: ifcoulomb = ', ifcoulomb
     error stop 'master_coulomb_end: invalid ifcoulomb'
  elseif(ifcoulomb.eq.0) then
     ! For ifcount.eq.0, zero all intent(out) results.
     dpcoulomb = 0._fp_kind
     dpcoulombf = 0._fp_kind
     dpcoulombt = 0._fp_kind
     dscoulomb = 0._fp_kind
     dscoulombf = 0._fp_kind
     dscoulombt = 0._fp_kind
     ducoulomb = 0._fp_kind
  else
     ! don't want tau(x) correction and diffraction correction
     ! simultaneously
     if(ifcoulomb.eq.2.and.if_dc.eq.1)&
          error stop 'master_coulomb_end: bad combination of ifcoulomb and if_dc'
     ! calculate change to thermodynamic quantities
     ! n.b. fcoulomb is per unit volume so is dscoulomb, ducoulomb
     ! n.b. fcoulomb is considered to be a function of n_e, n_i, T
     ! because sum0a(n_e, n_i), sum2a(n_e, n_i, T), and x(n_e, T),
     ! where x = thetax*n_e (ifcoulomb=2) or theta_e*n_e (if_dc = 1)
     ! partial -fcoulomb(n_i, n_e, t) wrt t
     ! n.b.
     ! xfixedt = partial x(ne,T)/partial T
     ! also note that
     ! partial fl(ne, T)/partial ln T = -rhostar(3)/rhostar(2)
     if(if_dc.eq.1) then
        xfixedt = (-theta_ef*rhostar(3)/rhostar(2))*n_e/t
        xfixedtf = xfixedt*(theta_eff/theta_ef +&
             rhostar(5)/rhostar(3) - rhostar(4)/rhostar(2) +&
             rhostar(2))
        xfixedtt = xfixedt*(theta_eft/theta_ef +&
             rhostar(6)/rhostar(3) - rhostar(5)/rhostar(2) +&
             rhostar(3) - 1._fp_kind)
     elseif(ifcoulomb.eq.2) then
        xfixedt = x*(1._fp_kind + pstar(2)*rhostar(3)/rhostar(2) -&
             pstar(3))/t
        xfixedtf = x*((2._fp_kind*rhostar(2) - pstar(2))*&
             (1._fp_kind + pstar(2)*rhostar(3)/rhostar(2) - pstar(3)) +&
             (pstar(4)*rhostar(3) + pstar(2)*&
             (rhostar(5)-rhostar(3)*rhostar(4)/rhostar(2)))/&
             rhostar(2) - pstar(5))/t
        xfixedtt = x*((2._fp_kind*rhostar(3) - pstar(3))*&
             (1._fp_kind + pstar(2)*rhostar(3)/rhostar(2) - pstar(3)) +&
             (pstar(5)*rhostar(3) + pstar(2)*&
             (rhostar(6)-rhostar(3)*rhostar(5)/rhostar(2)))/&
             rhostar(2) - pstar(6))/t
     else
        xfixedt = 0._fp_kind
        xfixedtf = 0._fp_kind
        xfixedtt = 0._fp_kind
     endif
     dscoulomb = -fcoulombt -&
          fcoulomb0*sum0atl/t - fcoulomb2*sum2atl/t -&
          fcoulombx*xfixedt
     ! partial fcoulombt(f,t) wrt ln f
     fcfttf =&
          fcoulomb0t*sum0f + fcoulomb2t*sum2f +&
          n_e*rhostar(2)*(fcoulomb0t*sum0ane +&
          fcoulomb2t*sum2ane) + fcoulombtx*xf
     ! partial fcoulombt(f,t) wrt ln t
     fcfttt =&
          fcoulomb0t*sum0t + fcoulomb2t*sum2t +&
          n_e*rhostar(3)*(fcoulomb0t*sum0ane +&
          fcoulomb2t*sum2ane) + fcoulomb0t*sum0atl +&
          fcoulomb2t*sum2atl + fcoulombtt*t + fcoulombtx*xt
     if(ifcoulomb.gt.4) then
        ! DeWitt definition
        ! sum0a = sum0 + theta_e*n_e
        ! partial sum0atl(f,t) wrt ln f
        sum0atlf = sum0atl*rhostar(2) + n_e*thetatf
        ! partial sum0atl(f,t) wrt ln t
        sum0atlt = sum0atl*rhostar(3) + n_e*thetatt
     else
        ! PTEH definition
        ! sum0a = sum0
        sum0atlf = 0._fp_kind
        sum0atlt = 0._fp_kind
     endif
     ! partial sum2atl(f,t) wrt ln f
     sum2atlf = sum2atl*rhostar(2) + n_e*thetatf
     ! partial sum2atl(f,t) wrt ln t
     sum2atlt = sum2atl*rhostar(3) + n_e*thetatt
     dscoulombf = -fcfttf -&
          fcft0f*sum0atl/t - fcoulomb0*sum0atlf/t -&
          fcft2f*sum2atl/t - fcoulomb2*sum2atlf/t -&
          fcftxf*xfixedt - fcoulombx*xfixedtf
     dscoulombt = -fcfttt -&
          fcft0t*sum0atl/t - fcft2t*sum2atl/t -&
          fcoulomb0*(sum0atlt/t - sum0atl/t) -&
          fcoulomb2*(sum2atlt/t - sum2atl/t) -&
          fcftxt*xfixedt - fcoulombx*xfixedtt
     ! from u = f + Ts
     ducoulomb = fcoulomb + t*dscoulomb
     ! dpcoulomb = - partial(V*fcoulomb(n_e, n_i, T)) wrt V
     !   = -fcoulomb + n_e*partial fcoulomb wrt n_e + sum_i n_i*
     !   partial fcoulomb wrt n_i
     !   = -fcoulomb + n_e*dmucoulomb + sum0*partial fcoulomb
     !   wrt sum0a + sum2*partial fcoulomb wrt sum2a
     ! n.b. the latter two terms are dropped if sum0 and sum2
     ! depend only on n_e.  This occurs if
     ! the PTEH approximation being used for sum0, sum2
     if(if_pteh.eq.1) then
        dpcoulomb = -fcoulomb + n_e*dmucoulomb
        ! dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
        ! fcoulomb (for the fully ionized case) is considered
        ! to be a function of (sum0(n_e), sum2(n_e)+
        ! theta_e(n_e,t)*n_e, t, x(n_e,T)) = function(n_e, T),
        ! completely independent of n_i.
        ! also recall dmucoulomb = partial fcoulomb(n_e, T) wrt n_e
        ! in the fully ionized case.
        ! There is a useful cancellation of terms when using
        ! this viewpoint....
        dpcoulombf = n_e*dmucoulombf
        ! dscoulomb = partial -fcoulomb(n_e, t) wrt t
        dpcoulombt = t*dscoulomb + n_e*dmucoulombt
        ! temporary check confirming definition of dmucoulomb
        ! dpcoulomb = fcoulomb
        ! dpcoulombf = dmucoulomb*n_e*rhostar(2)
        ! dpcoulombt = dmucoulomb*n_e*rhostar(3) - t*dscoulomb
     else
        dpcoulomb = -fcoulomb + n_e*dmucoulomb +&
             sum0*fcoulomb0 + sum2*fcoulomb2
        ! dpcoulomb(f,t)/ ln f and ln t:  to do these derivatives
        ! fcoulomb is considered a function of (sum0(f,t), sum2(f,t)+
        ! theta_e(n_e,t)*n_e, t, x(n_e,T)).  In other words, the previous
        ! n_i dependence is replaced with sum0 and sum2 dependence.
        ! There is a useful cancellation of terms when using
        ! this viewpoint....
        dpcoulombf = n_e*dmucoulombf +&
             sum0*fcft0f + sum2*fcft2f
        ! dscoulomb = partial -fcoulomb(sum0, sum2, n_e, t) wrt t
        dpcoulombt = t*dscoulomb + n_e*dmucoulombt +&
             sum0*fcft0t + sum2*fcft2t
     endif
  endif
end subroutine master_coulomb_end
