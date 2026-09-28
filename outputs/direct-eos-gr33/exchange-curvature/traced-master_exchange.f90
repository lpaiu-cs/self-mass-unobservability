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

! calculate transformed fermi_dirac integrals
! rhostar, pstar, sstar, ustar as a function of fl and tl.
! also calculate exchange if ifexchange_in != 0.
! input data:
! if abs(ifexchange_in) > 100 then use linear transform approximation
!   for exchange treatement.  Otherwise, use numerical transform as
!   described in research note.
! ifexchange = mod(ifexchange_in,100):
! One other wrinkle on deciding ifexchange is that large degeneracy
! approximations (mod(ifexchange,10) = 2) only allowed above psi_lim.
! ifexchange:
! To understand ifexchange, must summarize free-energy model of fex used
! in research note.  In general from kapusta relation,
! fex is proportional to I - J + 2^1.5 pi^2/3 beta^2 K, where
! K and J and known integrals which are functions of psi and beta,
! and I = K^2. (see Kovetz et al 1972, ApJ 174, 109, hereafter KLVH)
! 0  < ifexchange < 10 --> KLVH treatment (i.e., drop K term).
! 10 <= ifexchange < 20 --> Kapusta treatment (i.e., retain K term).
! 20 <= ifexchange < 30 --> special test case, use I term alone.
! 30 <= ifexchange < 40 --> special test case, use J term alone.
! 40 <= ifexchange < 50 --> special test case, use K term alone.
! N.B. if ifexchange is negative, then use non-relativistic limit
!   of corresponding positive ifexchange option.
! ifexchange details:
! ifexchange = 1 --> G(psi) + weak relativistic correction from KLVH
! ifexchange = 2 --> I-J degenerate expression from KLVH (with corrected
!   sign error on a2 and high numerical precision a2 and a3, see
!   paper IV)
! ifexchange = 11 or 12 K term added (both series from CG).
! ifexchange = 21 or 22 I term alone (both series from KLVH)
! ifexchange = 31 or 32 -J term alone (both series from KLVH).
! ifexchange = 41 or 42 K term alone.
! mod(ifexchange,10) = 4 is lowest order fit of J, K
! mod(ifexchange,10) = 5 is next higher order fit of J, K
! mod(ifexchange,10) = 6 is highest order fit of J, K

!> This master_exchange subroutine calculates transformed fermi_dirac
!> integrals (rhostar, pstar, sstar, ustar) as a function of fl
!> (determined from fl_in and iterative solution for an fl that is
!> consistent with the grand canonicial partition function version of
!> exchange) and tl.  In addition this subroutine calculates
!> quantities required to help determine the non-ideal exchange
!> component of equilibrium constants that are iteratively used to
!> help converge the EOS solution.  These calculated results include
!> intermediate quantities (saved in the mod_master_exchange_data
!> module) that are also needed for other purposes (e.g., calculation
!> of the non-ideal component of the pressure due to the exchange
!> effect) once the EOS is converged.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine master_exchange(verbosity, fl_in, tl, rhostar, pstar, sstar, ustar, morder,&
     ifexchange_in, dve, dvef, dvet)
  use mod_free_eos_constants, only: boltzmann, c_e, cpe, ct
  use mod_fermi_dirac, only: fermi_dirac
  use mod_master_exchange_data, only: last_fl_input,last_tl_input
  use mod_master_exchange_data, only:&
       iforder, ifstart,&
       dpsidf, dpsidf2,&
       dpsiprimedf, dpsiprimedf2,&
       fex, fexf, fext, fexf2, fexft, fext2, n_e, t,&
       fexprime, fexprimef, fexprimet, fexprimef2, fexprimeft, fexprimet2,&
       flprimef, flprimet, flprimef2, flprimeft, flprimet2,&
       muex, muexf, muext,&
       muex2, muex2f, muex2t,&
       p_e, pstarprime, psiprime, psi
  use mod_flow_data, only: ln_overflow_limit, ln_underflow_limit
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) :: verbosity, ifexchange_in, morder
  real(fp_kind), intent(in) :: fl_in, tl
  real(fp_kind), intent(out) :: rhostar(:), pstar(:), sstar(:), ustar(:)
  real(fp_kind), intent(out) :: dve, dvef, dvet

  ! Local variables
  real(fp_kind) tcl, fl
  ! variables for full-order approach
  real(fp_kind) psi_lim,&
       dflprime, flprime, n_eprime,&
       fzero, fzerox, fzerof, fzerot, fzerox2, fzeroxf, fzeroxt,&
       fzerof2, fzeroft, fzerot2,&
       fex_psi, fex_psif, fex_psit,&
       fex_psif2, fex_psift, fex_psit2

  real(fp_kind) rhostarprime(9), sstarprime(3), ustarprime(3)

  integer icount, ifexchange, nstar
  logical ifcontinue

  ! Mark this routine as called.
  if(ifstart.eq.1) ifstart = 0

  ! As some protection against an fl iteration that has gone really
  ! bad, pin the fl value used in the calculations below to within a
  ! very wide dynamic range which nevertheless has sensible margins
  ! that make it unlikely the results below will generate any serious
  ! floating-point exceptions.
  ! Maintenance, 2021.  Make these limits all the same in
  ! master_exchange, exchange_gcpf, and fermi_dirac.

  last_fl_input=fl_in;last_tl_input=tl
  fl = max(min(fl_in, ln_overflow_limit), ln_underflow_limit)

  nstar = size(sstar)

  ! Sanity check
  if(nstar.ne.3) error stop 'master_exchange: incorrect size of sstar'
  if(nstar+6.ne.size(rhostar).or.nstar+6.ne.size(pstar).or.nstar.ne.size(ustar))&
       error stop 'master_exchange: incorrect sizes for rhostar, pstar, or ustar'

  ! Silence spurious [-Wmaybe-uninitialized] gfortran warning for
  ! psi_lim. The reason why the warning is spurious is psi_lim is
  ! always initialized below if ifexchange_in is non-zero while all
  ! use cases for psi_lim occur if ifexchange.ne.0.  If ifexchange_in
  ! is non-zero, then ifexchange=mod(ifexchange_in,100) must also be
  ! non-zero so psi_lim is never used uninitialized.  However, gfortran
  ! cannot figure that out and therefore generates the spurious warning.
  psi_lim = 0._fp_kind

  if(ifexchange_in.ne.0) then
     ! iforder = 1 means exchange using first-order transformation from gcpf
     ! iforder = 2 means exchange using numerical transformation from gcpf
     if(abs(ifexchange_in).ge.100) then
        iforder = 1
     else
        iforder = 2
     endif
     ifexchange = mod(ifexchange_in,100)
     ! validity limit for Kovetz et al strong degeneracy series.
     psi_lim = max(3._fp_kind,0.5_fp_kind/(ct*exp(tl)))
  else
     ! signal that exchange not wanted.
     iforder = 0
  endif
  ! always calculate fermi_dirac (and return values) for fl.
  ! later may require primed values for non-linear transform
  ! approach (iforder = 2) which correspond to flprime.
  tcl = tl + log(ct)
  call fermi_dirac(verbosity, fl, tcl, rhostar, pstar, sstar, ustar, morder)

  if(iforder.eq.1) then
     ! first-order transformation from gcpf form
     ! large degeneracy approximation.
     if(mod(ifexchange,10).eq.2) then
        dpsidf = sqrt(1._fp_kind+exp(fl))
        psi = fl + 2._fp_kind*(dpsidf - log(1._fp_kind + dpsidf))
        if(psi.lt.psi_lim) then
           ifexchange = ifexchange -1
        endif
     endif
     call exchange_gcpf(&
          ifexchange, fl, tl,&
          psi, dpsidf, dpsidf2,&
          fex, fexf, fext, fexf2, fexft, fext2,&
          fex_psi, fex_psif, fex_psit,&
          fex_psif2, fex_psift, fex_psit2)
     ! number density of free electrons
     n_e = c_e*rhostar(1)
     t = exp(tl)
     ! electron exchange chemical potential (divided by kt)
     ! = partial fex(T, n_e)/partial n_e/k T
     muex = fexf/(n_e*rhostar(2)*boltzmann*t)
     muexf = muex*(fexf2/fexf - rhostar(2) - rhostar(4)/rhostar(2))
     muext = muex*(fexft/fexf - rhostar(3) - rhostar(5)/rhostar(2)&
          - 1._fp_kind)
     ! dve is *negative* chemical potential of electron/kT
     dve = -muex
     dvef = -muexf
     dvet = -muext
  elseif(iforder.eq.2) then
     ! full order approach (ignoring differential Coulomb effects)
     ! start with non-relativistic degeneracy treatment if considering
     ! option of strong degeneracy treatment
     if(mod(ifexchange,10).eq.2) then
        ifexchange = ifexchange -1
     endif
     call exchange_gcpf(&
          ifexchange, fl, tl,&
          psi, dpsidf, dpsidf2,&
          fex, fexf, fext, fexf2, fexft, fext2,&
          fex_psi, fex_psif, fex_psit,&
          fex_psif2, fex_psift, fex_psit2)
     ! number density of free electrons
     n_e = c_e*rhostar(1)
     t = exp(tl)
     ! needed for exchange_end and exchange_free entries later.
     p_e = cpe*pstar(1)
     ! solve for fl' such that
     ! n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
     ! where fex is the first-order free energy/volume = -kT ln Z_x
     ! first approximation for delta fl = fl' - fl (see below)
     ! function (see below) evaluated at fl' = fl
     fzero = - fex_psi/(boltzmann*t)
     ! first derivative wrt x = flprime (see below) evaluated at fl' = fl
     fzerox = n_e*rhostar(2)&
          - fex_psif/(boltzmann*t)
     ! second derivative wrt x = flprime (see below) evaluated at fl' = fl
     fzerox2 =&
          n_e*(rhostar(2)*rhostar(2) + rhostar(4))&
          - fex_psif2/(boltzmann*t)
     ! first linear estimate
     dflprime = -fzero/fzerox
     ! could solve quadratic, but use instead linear estimate of dflprime
     ! to estimate quadratic correction and store temporarily in fzerox2
     ! variable.
     fzerox2 = -0.5_fp_kind*dflprime*(dflprime*fzerox2/fzerox)
     ! make second-orer correction if significantly less than first-order
     ! estimate.
     if(abs(fzerox2).lt.0.5_fp_kind*abs(dflprime))&
          dflprime = dflprime + fzerox2
     flprime = fl
     ifcontinue = .true.
     icount = 0
     do while(ifcontinue)
        ! go one iteration beyond 1.d-7 which should give ~1.d-14 accuracy
        ! for quadratic convergence
        icount = icount + 1
        if(.false..and.verbosity.ge.4) write(stderr,*) 'icount, flprime, dflprime = ',&
             icount, flprime, dflprime
        ifcontinue = abs(dflprime).gt.1.e-7_fp_kind.and.icount.le.20
        ! fzero is monotonic in flprime so assured convergence so long
        ! as don't change flprime too wildly.  So far experience is
        ! that essentially no limits on change are fine so long as
        ! don't cross discontinuity at fl_limit
        flprime = flprime + max(-1.e1_fp_kind,min(1.e1_fp_kind,dflprime))
        call fermi_dirac(verbosity, flprime, tcl, rhostarprime, pstarprime, sstarprime, ustarprime, morder)
        call exchange_gcpf(&
             ifexchange, flprime, tl,&
             psiprime, dpsiprimedf, dpsiprimedf2,&
             fexprime, fexprimef, fexprimet,&
             fexprimef2, fexprimeft, fexprimet2,&
             fex_psi, fex_psif, fex_psit,&
             fex_psif2, fex_psift, fex_psit2)
        n_eprime = c_e*rhostarprime(1)
        ! n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
        ! where fex is the first-order free energy/volume = -kT ln Z_x
        fzero = n_eprime - n_e - fex_psi/(boltzmann*t)
        !           derivative wrt x = flprime
        fzerox = n_eprime*rhostarprime(2)&
             - fex_psif/(boltzmann*t)
        dflprime = -fzero/fzerox
     enddo
     if(icount.gt.20) then
        write(stderr,*) 'fl, tl = ', fl, tl
        write(stderr,*) 'flprime, dflprime = ', flprime, dflprime
        error stop 'master_exchange: first iteration failed to converge'
     endif
     if(mod(mod(ifexchange_in,100),10).eq.2 .and.psiprime.ge.psi_lim) then
        ifexchange = mod(ifexchange_in,100)
        ! if considering the possibility of a strong degeneracy treatment
        ! and if result for non-relativistic degeneracy is in appropriate
        ! regime with psiprime.ge.psi_lim
        ! then use strong degeneracy treatment
        call exchange_gcpf(&
             ifexchange, fl, tl,&
             psi, dpsidf, dpsidf2,&
             fex, fexf, fext, fexf2, fexft, fext2,&
             fex_psi, fex_psif, fex_psit,&
             fex_psif2, fex_psift, fex_psit2)
        ! number density of free electrons
        n_e = c_e*rhostar(1)
        t = exp(tl)
        ! solve for fl' such that
        ! n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
        ! where fex is the first-order free energy/volume = -kT ln Z_x
        ! first approximation for delta fl = fl' - fl (see below)
        ! function (see below) evaluated at fl' = fl
        fzero = - fex_psi/(boltzmann*t)
        ! first derivative wrt x = flprime (see below) evaluated at fl' = fl
        fzerox = n_e*rhostar(2)&
             - fex_psif/(boltzmann*t)
        ! second derivative wrt x = flprime (see below) evaluated at fl' = fl
        fzerox2 =&
             n_e*(rhostar(2)*rhostar(2) + rhostar(4))&
             - fex_psif2/(boltzmann*t)
        ! first linear estimate
        dflprime = -fzero/fzerox
        ! could solve quadratic, but use instead linear estimate of dflprime
        ! to estimate quadratic correction and store temporarily in fzerox2
        ! variable.
        fzerox2 = -0.5_fp_kind*dflprime*(dflprime*fzerox2/fzerox)
        ! make second-orer correction if significantly less than first-order
        ! estimate.
        if(abs(fzerox2).lt.0.5_fp_kind*abs(dflprime))&
             dflprime = dflprime + fzerox2
        flprime = fl
        ifcontinue = .true.
        icount = 0
        do while(ifcontinue)
           ! go one iteration beyond 1.d-7 which should give ~1.d-14 accuracy
           ! for quadratic convergence
           icount = icount + 1
           if(.false..and.verbosity.ge.4) write(stderr,*) 'icount, flprime, dflprime = ',&
                icount, flprime, dflprime
           ifcontinue = abs(dflprime).gt.1.e-7_fp_kind.and.icount.le.20
           ! fzero is monotonic in flprime so assured convergence so long
           ! as don't change flprime too wildly.  So far experience is
           ! that essentially no limits on change are fine so long as
           ! don't cross discontinuity at fl_limit
           flprime = flprime + max(-1.e1_fp_kind,min(1.e1_fp_kind,dflprime))
           call fermi_dirac(verbosity, flprime, tcl, rhostarprime, pstarprime, sstarprime, ustarprime, morder)
           call exchange_gcpf(&
                ifexchange, flprime, tl,&
                psiprime, dpsiprimedf, dpsiprimedf2,&
                fexprime, fexprimef, fexprimet,&
                fexprimef2, fexprimeft, fexprimet2,&
                fex_psi, fex_psif, fex_psit,&
                fex_psif2, fex_psift, fex_psit2)
           n_eprime = c_e*rhostarprime(1)
           ! n_e(fl') - n_e(fl) - (1/kT) partial fex(psi(fl'),T)/partial psi = 0
           ! where fex is the first-order free energy/volume = -kT ln Z_x
           fzero = n_eprime - n_e - fex_psi/(boltzmann*t)
           ! derivative wrt x = flprime
           fzerox = n_eprime*rhostarprime(2)&
                - fex_psif/(boltzmann*t)
           dflprime = -fzero/fzerox
        enddo
        if(icount.gt.20) then
           write(stderr,*) 'fl, tl = ', fl, tl
           write(stderr,*) 'flprime, dflprime = ', flprime, dflprime
           error stop 'master_exchange: second iteration failed to converge'
        endif
     endif
     ! flprime (and therefore all primed quantities) is an implicit
     ! function of fl and tl.  Find the derivatives from the chain
     ! rule for implicit functions.
     fzerof = -n_e*rhostar(2)
     fzerot = (n_eprime*rhostarprime(3) - n_e*rhostar(3))&
          - (fex_psit - fex_psi)/(boltzmann*t)
     flprimef = -fzerof/fzerox
     flprimet = -fzerot/fzerox
     !         second partial derivatives:
     fzerox2 =&
          n_eprime*(rhostarprime(2)*rhostarprime(2) + rhostarprime(4))&
          - fex_psif2/(boltzmann*t)
     fzeroxf = 0._fp_kind
     fzeroxt =&
          n_eprime*(rhostarprime(2)*rhostarprime(3) + rhostarprime(5))&
          - (fex_psift - fex_psif)/(boltzmann*t)
     fzerof2 = -n_e*(rhostar(2)*rhostar(2) + rhostar(4))
     fzeroft = -n_e*(rhostar(2)*rhostar(3) + rhostar(5))
     fzerot2 =&
          n_eprime*(rhostarprime(3)*rhostarprime(3) + rhostarprime(6))&
          - n_e*(rhostar(3)*rhostar(3) + rhostar(6))&
          - (fex_psit2 - 2._fp_kind*fex_psit + fex_psi)&
          /(boltzmann*t)
     flprimef2 = -(&
          fzerox2*flprimef*flprimef +&
          fzeroxf*flprimef + fzeroxf*flprimef +&
          fzerof2&
          )/fzerox
     flprimeft = -(&
          fzerox2*flprimef*flprimet +&
          fzeroxf*flprimet + fzeroxt*flprimef +&
          fzeroft&
          )/fzerox
     flprimet2 = -(&
          fzerox2*flprimet*flprimet +&
          fzeroxt*flprimet + fzeroxt*flprimet +&
          fzerot2&
          )/fzerox
     ! free energy/volume =
     ! = fex(psi') + kT(psi'-psi) n_e -kT (ln Z_e(psi')-ln Z_e(psi))/V
     ! = fex(psi') + kT(psi'-psi) n_e +
     !   (f_e(psi')-kt ne' psi') - (f_e(psi)-kt ne psi)
     ! = fex(psi') + kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
     ! = fex(psi') + kT(n_e-n_e') psi' + (f_e(psi') - f_e(psi))
     ! where n_e = n_e(psi) and n_e' = n_e(psi')
     ! chemical potential/kt from fex(psi')
     muex = fexprimef*flprimef/(n_e*rhostar(2)*boltzmann*t)
     muexf = muex*(&
          fexprimef2*flprimef/fexprimef&
          + flprimef2/flprimef - rhostar(2) - rhostar(4)/rhostar(2))
     muext = muex*(&
          (fexprimef2*flprimet + fexprimeft)/fexprimef&
          - 1._fp_kind +&
          flprimeft/flprimef - rhostar(3) - rhostar(5)/rhostar(2))
     ! chemical potential/kt from rest of free energy:
     ! n.b. by definition we have
     ! partial f_e(psi')/n_e' = kT psi'
     ! partial f_e(psi)/n_e = kT psi
     ! so there is much cancellation
     ! muex2 = psiprime - psi + (ne-ne') partial psi'/partial ne
     ! first all but psiprime-psi
     muex2 = -dpsiprimedf*flprimef/(n_e*rhostar(2))
     muex2f = muex2*(&
          (dpsiprimedf2/dpsiprimedf)*flprimef +&
          flprimef2/flprimef - rhostar(2) - rhostar(4)/rhostar(2))
     muex2t = muex2*(&
          (dpsiprimedf2/dpsiprimedf)*flprimet +&
          flprimeft/flprimef - rhostar(3) - rhostar(5)/rhostar(2))
     muex2f = (n_eprime*rhostarprime(2)*flprimef - n_e*rhostar(2))*&
          muex2 + (n_eprime - n_e)*muex2f
     muex2t = (n_eprime*(rhostarprime(2)*flprimet + rhostarprime(3))&
          - n_e*rhostar(3))*muex2 + (n_eprime - n_e)*muex2t
     muex2 = (n_eprime - n_e)*muex2
     ! now do psiprime - psi part
     muex2 = muex2 + psiprime - psi
     muex2f = muex2f + dpsiprimedf*flprimef - dpsidf
     muex2t = muex2t + dpsiprimedf*flprimet
     ! dve is *negative* chemical potential of electron/kT
     dve = - muex - muex2
     dvef = - muexf - muex2f
     dvet = - muext - muex2t
  elseif(iforder.eq.0) then
     dve = 0._fp_kind
     dvef = 0._fp_kind
     dvet = 0._fp_kind
  else
     error stop 'master_exchange: should not happen'
  endif
end subroutine master_exchange

!> This exchange_pressure subroutine calculates the change to the
!> pressure due to the exchange effect and the partial derivatives of
!> that change wrt ln f and ln T.
!>
!> \param[in] rhostar PARAMETERS NEED DOCUMENTATION
!>
subroutine exchange_pressure(rhostar, pstar, pex, pexf, pext)
  use mod_free_eos_constants, only: boltzmann, cpe
  use mod_master_exchange_data, only:&
       iforder, ifstart,&
       dpsidf,&
       dpsiprimedf,&
       fex, fexf, fext, n_e, t,&
       fexprime, fexprimef, fexprimet,&
       flprimef, flprimet,&
       muex, muexf, muext,&
       muex2, muex2f, muex2t,&
       p_e, pstarprime, psiprime, psi

  ! Arguments
  real(fp_kind), intent(in) :: rhostar(:), pstar(:)
  real(fp_kind), intent(out) :: pex, pexf, pext

  ! Local variables
  real(fp_kind) p_eprime, fex2, fex2f, fex2t

  ! Sanity checks
  if(ifstart.eq.1) error stop 'exchange_pressure: master_exchange must be called first'
  if(size(rhostar).ne.9.or.size(pstar).ne.9)&
       error stop 'exchange_pressure: incorrect sizes for rhostar or pstar'

  if(iforder.eq.1) then
     ! electron exchange pressure
     ! = - fex + ne (partial fex(T,ne)/partial ne)
     pex = -fex + n_e*boltzmann*t*muex
     pexf = -fexf + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
     pext = -fext + n_e*boltzmann*t*((rhostar(3)+1._fp_kind)*muex + muext)
  elseif(iforder.eq.2) then
     ! full order approach (ignoring differential Coulomb effects)
     ! free energy/volume =
     ! = fex(psi') + kT(psi'-psi) n_e -kT (ln Z_e(psi')-ln Z_e(psi))/V
     ! = fex(psi') + kT(psi'-psi) n_e +
     !   (f_e(psi')-kt ne' psi') - (f_e(psi)-kt ne psi)
     ! = fex(psi') + kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
     ! = fex(psi') + kT(n_e-n_e') psi' + (f_e(psi') - f_e(psi))
     ! where n_e = n_e(psi) and n_e' = n_e(psi')
     ! first do just electron exchange pressure from fex(psi')
     ! pex  = - fex' + ne (partial fex'(T,ne)/partial ne)
     pex = - fexprime + n_e*boltzmann*t*muex
     pexf = -fexprimef*flprimef&
          + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
     pext = -(fexprimef*flprimet + fexprimet)&
          + n_e*boltzmann*t*((rhostar(3)+1._fp_kind)*muex + muext)
     ! remainder of exchange free energy in full-order approximation
     ! fex2  = kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
     p_eprime = cpe*pstarprime(1)
     fex2 = boltzmann*t*(psiprime - psi)*n_e - (p_eprime - p_e)
     fex2f = boltzmann*t*((dpsiprimedf*flprimef - dpsidf)*n_e +&
          (psiprime - psi)*n_e*rhostar(2)) -&
          (p_eprime*pstarprime(2)*flprimef - p_e*pstar(2))
     fex2t = boltzmann*t*(&
          (psiprime - psi + dpsiprimedf*flprimet)*n_e +&
          (psiprime - psi)*n_e*rhostar(3)) -&
          (p_eprime*(pstarprime(2)*flprimet + pstarprime(3)) -&
          p_e*pstar(3))
     ! electron exchange pressure
     ! = - fex2 + ne (partial fex2(T,ne)/partial ne)
     pex = pex - fex2 + n_e*boltzmann*t*muex2
     pexf = pexf - fex2f + n_e*boltzmann*t*&
          (rhostar(2)*muex2 + muex2f)
     pext = pext - fex2t + n_e*boltzmann*t*&
          ((rhostar(3) + 1._fp_kind)*muex2 + muex2t)
  elseif(iforder.eq.0) then
     pex = 0._fp_kind
     pexf = 0._fp_kind
     pext = 0._fp_kind
  else
     error stop 'exchange_pressure: internal error that should not happen'
  endif
end subroutine exchange_pressure

!> This exchange_free subroutine calculates the change to the free
!> energy per unit volume due to the exchange effect and the partial
!> derivative of that change wrt fl.
!>
!> \param[in] rhostar PARAMETERS NEED DOCUMENTATION
!>
subroutine exchange_free(rhostar, pstar, free_ex, free_exf)
  use mod_free_eos_constants, only: boltzmann, cpe
  use mod_master_exchange_data, only:&
       iforder, ifstart,&
       dpsidf,&
       dpsiprimedf,&
       fex, fexf, n_e, t,&
       fexprime, fexprimef,&
       flprimef,&
       p_e, pstarprime, psiprime, psi

  ! Arguments
  real(fp_kind), intent(in) :: rhostar(:), pstar(:)
  real(fp_kind), intent(out) :: free_ex, free_exf

  ! Local variables
  real(fp_kind) p_eprime, fex2, fex2f

  ! Sanity checks
  if(ifstart.eq.1) error stop 'exchange_free: master_exchange must be called first'
  if(size(rhostar).ne.9.or.size(pstar).ne.9)&
       error stop 'exchange_free: incorrect sizes for rhostar or pstar'

  if(iforder.eq.1) then
     free_ex = fex
     free_exf = fexf
  elseif(iforder.eq.2) then
     !         remainder of exchange free energy in full-order approximation
     !         fex2  = kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
     p_eprime = cpe*pstarprime(1)
     fex2 = boltzmann*t*(psiprime - psi)*n_e - (p_eprime - p_e)
     fex2f = boltzmann*t*((dpsiprimedf*flprimef - dpsidf)*n_e +&
          (psiprime - psi)*n_e*rhostar(2)) -&
          (p_eprime*pstarprime(2)*flprimef - p_e*pstar(2))
     free_ex = fexprime + fex2
     free_exf = fexprimef*flprimef + fex2f
  elseif(iforder.eq.0) then
     free_ex = 0._fp_kind
     free_exf = 0._fp_kind
  else
     error stop 'exchange_free: should not happen'
  endif
end subroutine exchange_free

!> This exchange_end subroutine calculates the change to the pressure,
!> entropy per unit volume, and internal energy per unit volume due to
!> the exchange effect as well as the partial derivatives of those
!> first two quantities wrt fl and ln T.
!>
!> \param[in] rhostar PARAMETERS NEED DOCUMENTATION
!>
subroutine exchange_end(rhostar, pstar, pex, pext, pexf, sex, sexf, sext, uex)
  use mod_free_eos_constants, only: boltzmann, cpe
  use mod_master_exchange_data, only:&
       iforder, ifstart,&
       dpsidf, dpsidf2,&
       dpsiprimedf, dpsiprimedf2,&
       fex, fexf, fext, fexf2, fexft, fext2, n_e, t,&
       fexprime, fexprimef, fexprimet, fexprimef2, fexprimeft, fexprimet2,&
       flprimef, flprimet, flprimef2, flprimeft, flprimet2,&
       muex, muexf, muext,&
       muex2, muex2f, muex2t,&
       p_e, pstarprime, psiprime, psi

  ! Arguments
  real(fp_kind), intent(in) :: rhostar(:), pstar(:)
  real(fp_kind), intent(out) :: pex, pext, pexf, sex, sexf, sext, uex

  ! Local variables
  real(fp_kind) p_eprime, fex2, fex2f, fex2t, fex2f2, fex2ft, fex2t2, sex2

  ! Sanity checks
  if(ifstart.eq.1) error stop 'exchange_end: master_exchange must be called first'
  if(size(rhostar).ne.9.or.size(pstar).ne.9)&
       error stop 'exchange_end: incorrect sizes for rhostar or pstar'

  if(iforder.eq.1) then
     ! electron exchange pressure
     ! = - fex + ne (partial fex(T,ne)/partial ne)
     pex = -fex + n_e*boltzmann*t*muex
     pexf = -fexf + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
     pext = -fext + n_e*boltzmann*t*((rhostar(3)+1._fp_kind)*muex + muext)
     ! electron exchange entropy per unit volume
     ! = - partial fex(T,n_e)/partial T
     ! = - (partial fex(tl,fl)/partial tl +
     !   partial fex(tl,fl)/partial fl *
     !   partial fl(tl,n_e)/partial tl)/T, where
     !   partial fl(tl,n_e)/partial tl = -rhostar(3)/rhostar(2)
     sex = -(fext - fexf*rhostar(3)/rhostar(2))/t
     sexf = -(fexft&
          - (fexf2/fexf + rhostar(5)/rhostar(3)&
          - rhostar(4)/rhostar(2))&
          *fexf*rhostar(3)/rhostar(2))/t
     sext = - sex - (fext2&
          - (fexft/fexf + rhostar(6)/rhostar(3)&
          - rhostar(5)/rhostar(2))&
          *fexf*rhostar(3)/rhostar(2))/t
     uex = fex + t*sex
  elseif(iforder.eq.2) then
     ! full order approach (ignoring differential Coulomb effects)
     ! free energy/volume =
     ! = fex(psi') + kT(psi'-psi) n_e -kT (ln Z_e(psi')-ln Z_e(psi))/V
     ! = fex(psi') + kT(psi'-psi) n_e +
     !   (f_e(psi')-kt ne' psi') - (f_e(psi)-kt ne psi)
     ! = fex(psi') + kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
     ! = fex(psi') + kT(n_e-n_e') psi' + (f_e(psi') - f_e(psi))
     ! where n_e = n_e(psi) and n_e' = n_e(psi')
     ! first do just electron exchange pressure from fex(psi')
     ! pex  = - fex' + ne (partial fex'(T,ne)/partial ne)
     pex = - fexprime + n_e*boltzmann*t*muex
     pexf = -fexprimef*flprimef&
          + n_e*boltzmann*t*(rhostar(2)*muex + muexf)
     pext = -(fexprimef*flprimet + fexprimet)&
          + n_e*boltzmann*t*((rhostar(3)+1._fp_kind)*muex + muext)
     ! electron exchange entropy per unit volume
     !   = - partial fex'(T,n_e)/partial T
     !   = - (partial fex'(tl,fl)/partial tl +
     !     partial fex'(tl,fl)/partial fl *
     !     partial fl(tl,n_e)/partial tl)/T, where
     !     partial fl(tl,n_e)/partial tl = -rhostar(3)/rhostar(2)
     sex = -(fexprimef*flprimet + fexprimet&
          - fexprimef*flprimef*rhostar(3)/rhostar(2))/t
     sexf = -(fexprimef2*flprimef*flprimet + fexprimef*flprimeft&
          + fexprimeft*flprimef&
          - (fexprimef2*flprimef/fexprimef + flprimef2/flprimef&
          + rhostar(5)/rhostar(3)&
          - rhostar(4)/rhostar(2))&
          *fexprimef*flprimef*rhostar(3)/rhostar(2))/t
     sext = - sex - (&
          (fexprimef2*flprimet + fexprimeft)*flprimet&
          + fexprimef*flprimet2&
          + fexprimeft*flprimet + fexprimet2&
          - ((fexprimef2*flprimet + fexprimeft)/fexprimef&
          + flprimeft/flprimef&
          + rhostar(6)/rhostar(3)&
          - rhostar(5)/rhostar(2))&
          *fexprimef*flprimef*rhostar(3)/rhostar(2))/t
     uex = fexprime + t*sex
     ! remainder of exchange free energy in full-order approximation
     ! fex2  = kT(psi'-psi) n_e - (P_e(psi') - P_e(psi))
     p_eprime = cpe*pstarprime(1)
     fex2 = boltzmann*t*(psiprime - psi)*n_e - (p_eprime - p_e)
     fex2f = boltzmann*t*((dpsiprimedf*flprimef - dpsidf)*n_e +&
          (psiprime - psi)*n_e*rhostar(2)) -&
          (p_eprime*pstarprime(2)*flprimef - p_e*pstar(2))
     fex2t = boltzmann*t*(&
          (psiprime - psi + dpsiprimedf*flprimet)*n_e +&
          (psiprime - psi)*n_e*rhostar(3)) -&
          (p_eprime*(pstarprime(2)*flprimet + pstarprime(3)) -&
          p_e*pstar(3))
     fex2f2 = boltzmann*t*((dpsiprimedf2*flprimef*flprimef +&
          dpsiprimedf*flprimef2 - dpsidf2)*n_e +&
          2._fp_kind*(dpsiprimedf*flprimef - dpsidf)*n_e*rhostar(2) +&
          (psiprime - psi)*n_e*(rhostar(2)*rhostar(2) + rhostar(4))) -&
          (p_eprime*(pstarprime(2)*pstarprime(2) + pstarprime(4))*&
          flprimef*flprimef +&
          p_eprime*pstarprime(2)*flprimef2 -&
          p_e*(pstar(2)*pstar(2) + pstar(4)))
     fex2ft = boltzmann*t*(&
          (dpsiprimedf*flprimef - dpsidf +&
          dpsiprimedf2*flprimef*flprimet + dpsiprimedf*flprimeft)*n_e +&
          (psiprime - psi + dpsiprimedf*flprimet)*n_e*rhostar(2) +&
          (dpsiprimedf*flprimef - dpsidf)*n_e*rhostar(3) +&
          (psiprime - psi)*n_e*(rhostar(2)*rhostar(3) + rhostar(5))) -&
          (p_eprime*pstarprime(2)*flprimef*&
          (pstarprime(2)*flprimet + pstarprime(3)) +&
          p_eprime*(pstarprime(4)*flprimef*flprimet +&
          pstarprime(2)*flprimeft + pstarprime(5)*flprimef) -&
          p_e*(pstar(2)*pstar(3) + pstar(5)))
     fex2t2 = boltzmann*t*(&
          (psiprime - psi + dpsiprimedf*flprimet)*n_e +&
          (psiprime - psi)*n_e*rhostar(3) +&
          (dpsiprimedf*flprimet +&
          dpsiprimedf2*flprimet*flprimet + dpsiprimedf*flprimet2)*n_e +&
          (psiprime - psi + dpsiprimedf*flprimet)*n_e*rhostar(3) +&
          (dpsiprimedf*flprimet)*n_e*rhostar(3) +&
          (psiprime - psi)*n_e*(rhostar(3)*rhostar(3) + rhostar(6))) -&
          (p_eprime*(pstarprime(2)*flprimet + pstarprime(3))*&
          (pstarprime(2)*flprimet + pstarprime(3)) +&
          p_eprime*(pstarprime(4)*flprimet*flprimet +&
          pstarprime(5)*flprimet +&
          pstarprime(2)*flprimet2 +&
          pstarprime(5)*flprimet + pstarprime(6)) -&
          p_e*(pstar(3)*pstar(3) + pstar(6)))
     ! electron exchange pressure
     ! = - fex2 + ne (partial fex2(T,ne)/partial ne)
     pex = pex - fex2 + n_e*boltzmann*t*muex2
     pexf = pexf - fex2f + n_e*boltzmann*t*&
          (rhostar(2)*muex2 + muex2f)
     pext = pext - fex2t + n_e*boltzmann*t*&
          ((rhostar(3) + 1._fp_kind)*muex2 + muex2t)
     ! electron exchange entropy per unit volume
     ! = - partial fex2(T,n_e)/partial T
     ! = - (partial fex2(tl,fl)/partial tl +
     !   partial fex2(tl,fl)/partial fl *
     !   partial fl(tl,n_e)/partial tl)/T, where
     !     partial fl(tl,n_e)/partial tl = -rhostar(3)/rhostar(2)
     sex2 = -(fex2t - fex2f*rhostar(3)/rhostar(2))/t
     sex = sex + sex2
     sexf = sexf - (fex2ft -(fex2f2*rhostar(3) +&
          fex2f*(rhostar(5) - rhostar(3)*rhostar(4)/rhostar(2)))/&
          rhostar(2))/t
     sext = sext - sex2 -&
          (fex2t2 -(fex2ft*rhostar(3) +&
          fex2f*(rhostar(6) - rhostar(3)*rhostar(5)/rhostar(2)))/&
          rhostar(2))/t
     ! electron exchange internal energy per unit volume
     ! from u = f + T s
     uex = uex + fex2 + t*sex2
  elseif(iforder.eq.0) then
     pex = 0._fp_kind
     pexf = 0._fp_kind
     pext = 0._fp_kind
     sex = 0._fp_kind
     sexf = 0._fp_kind
     sext = 0._fp_kind
     uex = 0._fp_kind
  else
     error stop 'exchange_end: should not happen'
  endif
end subroutine exchange_end
