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
!
!*******************************************************************************

! Also, calculate
! free energy/V in linear inversion approximation. this free
! energy is given by fex = -kT ln Z_x/V, where ln Z_x is the grand
! canonical partition function (gcpf) exchange term described in the
! research note.  fex is used either to get the free energy/V in the
! linear inversion approximation or else is used by the calling routine
! to solve the full inversion problem numerically.

! Also calculate higher order partial derivatives of fex wrt
! fl, tl, fl2, fltl, tl2.

! Also calculate fex_psi = partial fex(t,psi) wrt psi and higher order
! partial derivatives of fex_psi wrt to fl, tl, fl2, fltl, tl2.
!
! input data:
! ifexchange:
!   To understand ifexchange, must summarize free-energy model of fex used
!   in research note.  In general from kapusta relation,
!   fex is proportional to I - J + 2^1.5 pi^2/3 beta^2 K, where
!   K and J and known integrals which are functions of psi and beta,
!   and I = K^2. (see Kovetz et al 1972, ApJ 174, 109, hereafter KLVH)
!    0  < ifexchange < 10 --> KLVH treatment (i.e., drop K term).
!   10 <= ifexchange < 20 --> Kapusta treatment (i.e., retain K term).
!   20 <= ifexchange < 30 --> special test case, use I term alone.
!   30 <= ifexchange < 40 --> special test case, use J term alone.
!   40 <= ifexchange < 50 --> special test case, use K term alone.
!   N.B. if ifexchange is negative, then use non-relativistic limit
!     of corresponding positive ifexchange option.
!   ifexchange details:
!   ifexchange = 1 --> G(psi) + weak relativistic correction from KLVH
!   ifexchange = 2 --> I-J degenerate expression from KLVH (with corrected
!     sign error on a2 and high numerical precision a2 and a3, see
!     paper IV)
!   ifexchange = 11 or 12 K term added (both series from CG).
!   ifexchange = 21 or 22 I term alone (both series from KLVH)
!   ifexchange = 31 or 32 -J term alone (both series from KLVH).
!   ifexchange = 41 or 42 K term alone.
!   mod(ifexchange,10) = 4 is lowest order fit of J, K
!   mod(ifexchange,10) = 5 is next higher order fit of J, K
!   mod(ifexchange,10) = 6 is highest order fit of J, K
!   fl = ln f where f is EFF degeneracy parameter
!   tl = ln T.
!   output data:
!   psi = CG eta degeneracy parameter and derivatives wrt fl.
!   fex (explained above) and various partial derivatives.

!> This exchange_gcpf subroutine calculates the grand canonical partition
!> function form of the component of the free energy per unit volume due to the exchange effect
!> as well as first and second (mixed) partial derivatives of that quantity wrt fl and ln T.
!>
!> \param[in] ifexchange PARAMETERS NEED DOCUMENTATION
!>
subroutine exchange_gcpf(&
     ifexchange, fl_in, tl,&
     psi, dpsidf, dpsidf2,&
     fex, fexf, fext,&
     fexf2, fexft, fext2,&
     fex_psi, fex_psif, fex_psit,&
     fex_psif2, fex_psift, fex_psit2)

  use mod_free_eos_constants, only: pi, c2, clight, cr, ct, echarge, electron_mass
  use mod_fermi_dirac, only: fermi_dirac_ct, effsum_calc
  use mod_flow_data, only: ln_overflow_limit, ln_underflow_limit
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) :: ifexchange
  real(fp_kind), intent(in) :: fl_in, tl
  real(fp_kind), intent(out) ::&
     psi, dpsidf, dpsidf2,&
     fex, fexf, fext,&
     fexf2, fexft, fext2,&
     fex_psi, fex_psif, fex_psit,&
     fex_psif2, fex_psift, fex_psit2

  ! Local variables
  real(fp_kind)&
       con_exchange, fex0,&
       fex0f, fex0g, fex0ff, fex0fg, fex0gg,&
       fex0fff, fex0ffg, fex0fgg, fex0ggg,&
       fex0log, fex0log2, d1fex0log2, d2fex0log2, d3fex0log2,&
       fl, f, t,&
       fpsi, dfpsi, d2fpsi, d3fpsi,&
       dfpsidf, dfpsidf2, ddfpsidf, ddfpsidf2,&
       dpsidf3
  ! variables needed for relativistic corrections
  real(fp_kind) con_nr1_ratio, con_d1_ratio, beta,&
       fd32(6), fd12(5),&
       fex1, fex1f, fex1t,&
       fex1f2, fex1ft, fex1t2,&
       fex1_psi, fex1_psif, fex1_psit,&
       fex1_psif2, fex1_psift, fex1_psit2,&
       fex2, fex2f, fex2t,&
       fex2f2, fex2ft, fex2t2,&
       fex2_psi, fex2_psif, fex2_psit,&
       fex2_psif2, fex2_psift, fex2_psit2
  integer ideriv, maxderiv
  parameter (maxderiv=4)
  real(fp_kind) fd(maxderiv), dfd(maxderiv-1)
  integer iffirst
  data iffirst/1/
  ! large degeneracy case variables with parameters taken from KLVH et
  ! al. paper.
  real(fp_kind) a1, a2, con2, con4,&
       z, x, xz, xz2, xz3,&
       phi, phiz,&
       lnalpha, lnalphaz, lnalphaz2, lnalphaz3
  ! parameter (a1 = 0.449d0)  !eq. 21
  !   this coefficient is negative of that in paper, but gives vastly
  !   improved agreement with G(psi) asympototic series
  ! parameter (a2 = 0.504d0) !eq. 22 (n.b. negated!)
  ! parameter (a1 = 0.44904d0) !Lamb and Graziani (1997, private communication)
  ! parameter (a2 = 0.50475d0) !Lamb and Graziani (1997, private communication)
  ! These two values the result of running
  ! ./exchange_a1_direct
  ! A1 =         0.4490427560
  ! ./exchange_a2_direct
  ! A2 =         0.5047526561
  ! These results are derived using similar numerical integration used
  ! for other exchange integral testing which should be good to
  ! nominal precision of 10^{-9}.  Actually, various tests with
  ! convergence criteria indicates above results are good to 1 digit
  ! in last place.
  parameter (a1 = 0.4490427560_fp_kind)
  parameter (a2 = 0.5047526561_fp_kind)
  parameter (con2 = pi*pi/12._fp_kind)  !eq. 8
  parameter (con4 = 7._fp_kind*pi*pi*pi*pi/720._fp_kind)  !eq.8
  ! real(fp_kind) con6
  ! parameter (con6 = 31.d0*pi*pi*pi*pi*pi*pi/30240.d0)  !CG 24.182
  ! variables needed for EFF-style summations and derivatives.
  integer mordermax, maxc, maxkind,&
       ikind, nderiv
  parameter (mordermax = 9)
  parameter (maxc = (mordermax+1)*(mordermax+1))
  parameter (maxkind = 6)
  parameter (nderiv = 9)
  real(fp_kind) ccoeff(maxc, maxkind)
  integer mforder(maxkind), mgorder(maxkind), n_max
  real(fp_kind) sum(1+nderiv),  rmforder(maxkind), rmgorder(maxkind), rn_max
  real(fp_kind) fexi, fexif, fexit,&
       fexif2, fexift, fexit2,&
       fexi_psi, fexi_psif, fexi_psit,&
       fexi_psif2, fexi_psift, fexi_psit2,&
       fexj, fexjf, fexjt,&
       fexjf2, fexjft, fexjt2,&
       fexj_psi, fexj_psif, fexj_psit,&
       fexj_psif2, fexj_psift, fexj_psit2,&
       fexk, fexkf, fexkt,&
       fexkf2, fexkft, fexkt2,&
       fexk_psi, fexk_psif, fexk_psit,&
       fexk_psif2, fexk_psift, fexk_psit2,&
       fexl, fexlf, fexlt,&
       fexlf2, fexlft, fexlt2,&
       fexl_psi, fexl_psif, fexl_psit,&
       fexl_psif2, fexl_psift, fexl_psit2
  real(fp_kind) wf, vf, uf, duf, duf2, r_uf, r_duf, r_duf2, g, vg, ug, fdf, dug, dug2,&
       power
  real(fp_kind) acon, bcon, ccon, dcon, econ
  integer index_exchange
  logical iffexi, iffexj, iffexk
  data index_exchange /0/
  save iffirst, con_exchange, con_nr1_ratio, con_d1_ratio,&
       index_exchange, ccoeff, mgorder, mforder, rmgorder, rmforder

  if(iffirst.eq.1) then
     iffirst = 0
     ! m_e * e * k/h^2
     ! = (m_e*avogadro) * e * k^2/(h^2 * (k*avogadro))
     con_exchange = electron_mass*echarge*clight*clight/&
          (c2*c2*cr)
     ! in non-relativistic degenerate limit,
     ! free energy/volume = con_exchange*t*t*G(psi)
     ! thus, coefficient of G(psi) is con_exchange*t*t
     con_exchange = -8._fp_kind*pi*con_exchange*con_exchange
     ! coefficient of F_{1/2} is con_nr1_ratio*beta^1.5 times coefficient of
     ! G(psi).  Also, negative of ratio between I-J and K integral.
     con_nr1_ratio = -2._fp_kind*pi*pi*sqrt(2._fp_kind)/3._fp_kind
     ! coefficient of I-J is con_d1_ratio/t^2 times coefficient of G(psi)
     ! where con_d1_ratio = -0.5*(m_e*c^2/k)^2.  This also works out from
     ! eq. 6a of KLVH.
     con_d1_ratio = -0.5_fp_kind*(electron_mass*clight*clight/cr)**2
  endif
  t = exp(tl)
  beta = ct*t

  ! As some protection against an fl iteration that has gone really
  ! bad, pin the fl value used in the calculations below to within a
  ! very wide dynamic range which nevertheless has sensible margins
  ! that make it unlikely the results below will generate any serious
  ! floating-point exceptions.
  ! Maintenance, 2021.  Make these limits all the same in
  ! master_exchange, exchange_gcpf, and fermi_dirac.

  fl = max(min(fl_in, ln_overflow_limit), ln_underflow_limit)

  f = exp(fl)
  ! partial psi/d ln f and higher derivatives
  dpsidf = sqrt(1._fp_kind+f)
  ! dpsidf2 = 0.5d0*f*dpsidf/(1.d0+f)
  dpsidf2 = 0.5_fp_kind*f/dpsidf
  dpsidf3 = dpsidf2*(1._fp_kind - dpsidf2/dpsidf)
  psi = fl + 2._fp_kind*(dpsidf - log(1._fp_kind + dpsidf))
  if(mod(abs(ifexchange),10).eq.1) then
     ! required for most cases so don't bother with logic to exclude
     fd32(1) = fermi_dirac_ct(psi, 0, 0)
     do ideriv = 1,4
        fd32(ideriv+1) = fermi_dirac_ct(psi, ideriv, 0)
        fd12(ideriv) = (2._fp_kind/3._fp_kind)*fd32(ideriv+1)
     enddo
     if(abs(ifexchange).eq.1.or.abs(ifexchange).eq.11.or.&
          abs(ifexchange).eq.31) then
        ! These cases require non-relativistic limit of J
        call f_psi(psi, fpsi, dfpsi, d2fpsi, d3fpsi)
        dfpsidf = dfpsi*dpsidf
        dfpsidf2 = d2fpsi*dpsidf*dpsidf + dfpsi*dpsidf2
        ddfpsidf = d2fpsi*dpsidf
        ddfpsidf2 = d3fpsi*dpsidf*dpsidf + d2fpsi*dpsidf2
        ! definition of fps derived from definitions of fp12 and exct (see their
        ! commentary):
        ! fpsi(psi) =
        ! integral from -inf to psi of square of derivative of fp12(psi') d psi'
        ! fp12(psi) is C&G F(1/2, psi)/Gamma(3/2), DeWitt's, script
        ! capital I(1/2) and Rogers and DeWitt's script capital F(1/2).
        ! Now, Gamma(3/2)^2 = pi/4, and d F(1/2, psi)/d psi = 1/2 F(-1/2, psi).
        ! Thus, fpsi(psi) = G(psi)/pi where G(psi) is defined in paper.
        !
        ! electron exchange free energy per unit volume.
        ! see DeWitt 1969 in "Low Luminosity Stars" ed. S. Kumar, eq. 17
        ! for the canonical ensemble equation we are programming.
        ! see also Rogers and DeWitt, 1973, Phys. Rev. A8, 1061 eq. 45
        ! for the equivalent equation for the corresponding *grand*
        ! canonical ensemble equation.
        ! extra factor of pi needed because fpsi = G(psi)/pi
        fex0 = con_exchange*pi*t*t
        fex = fex0*fpsi
        fexf = fex0*dfpsidf
        fext = fex*2._fp_kind
        fexf2 = fex0*dfpsidf2
        fexft = fexf*2._fp_kind
        fext2 = fext*2._fp_kind
        fex_psi = fex0*dfpsi
        fex_psif = fex0*ddfpsidf
        fex_psit = fex_psi*2._fp_kind
        fex_psif2 = fex0*ddfpsidf2
        fex_psift = fex_psif*2._fp_kind
        fex_psit2 = fex_psit*2._fp_kind
     elseif(ifexchange.eq.-21) then
        ! Term from Non-relativistic I = K^2  integral
        ! extra factor is -4 beta^3/(2 beta^2)
        fex0 = con_exchange*t*t*(-2._fp_kind*beta)
        fex = fex0*fd12(1)*fd12(1)
        fexf = 2._fp_kind*fex0*fd12(1)*fd12(2)*dpsidf
        fext = 3._fp_kind*fex
        fexf2 = 2._fp_kind*fex0*(&
             (fd12(2)*fd12(2) + fd12(1)*fd12(3))*dpsidf*dpsidf +&
             fd12(1)*fd12(2)*dpsidf2)
        fexft = 3._fp_kind*fexf
        fext2 = 3._fp_kind*fext
        fex_psi = 2._fp_kind*fex0*fd12(1)*fd12(2)
        fex_psif = 2._fp_kind*fex0*(&
             fd12(2)*fd12(2) + fd12(1)*fd12(3))*dpsidf
        fex_psit = 3._fp_kind*fex_psi
        fex_psif2 = 2._fp_kind*fex0*(&
             (3._fp_kind*fd12(2)*fd12(3) + fd12(1)*fd12(4))*dpsidf*dpsidf +&
             (fd12(2)*fd12(2) + fd12(1)*fd12(3))*dpsidf2)
        fex_psift = 3._fp_kind*fex_psif
        fex_psit2 = 3._fp_kind*fex_psit
     elseif(ifexchange.eq.-41) then
        ! Term from non-relativistic K integral
        ! extra factor is con_nr1_ratio*beta^2*(2 beta^1.5/(2 beta^2)
        fex0 = con_exchange*t*t*(con_nr1_ratio*beta*sqrt(beta))
        fex = fex0*fd12(1)
        fexf = fex0*fd12(2)*dpsidf
        fext = 3.5_fp_kind*fex
        fexf2 = fex0*(fd12(3)*dpsidf*dpsidf + fd12(2)*dpsidf2)
        fexft = 3.5_fp_kind*fexf
        fext2 = 3.5_fp_kind*fext
        fex_psi = fex0*fd12(2)
        fex_psif = fex0*fd12(3)*dpsidf
        fex_psit = 3.5_fp_kind*fex_psi
        fex_psif2 = fex0*(fd12(4)*dpsidf*dpsidf + fd12(3)*dpsidf2)
        fex_psift = 3.5_fp_kind*fex_psif
        fex_psit2 = 3.5_fp_kind*fex_psit
     endif
     ! N.B. nothing calculated in above blocks
     ! for ifexchange = 21 or ifexchange = 41 cases.
     if(0.lt.ifexchange.and.ifexchange.lt.40) then
        ! ifexchange = 1, 11, 21, 31
        ! second Kapusta term is equivalent to KLVH so use their
        ! weakly relativistic series.
        ! free energy per unit volume in non-relativisitic case is
        ! con_exchange*t*t*G(psi), and here we add on F_{1/2}F_{1/2} and
        ! F_{1/2}F_{3/2} terms from KLVH eqs. 38, 36, and 36.  acon is
        ! the overall multiplier compared to NR case, while bcon is
        ! the ratio of the second to first term for relativistic part of
        ! J-I, -I, and J.
        if(ifexchange.eq.1.or.ifexchange.eq.11) then
           ! relativistic part of J-I
           acon = -1.5_fp_kind
           bcon = 0.75_fp_kind
        elseif(ifexchange.eq.21) then
           ! -I
           acon = -2._fp_kind
           bcon = 0.5_fp_kind
        elseif(ifexchange.eq.31) then
           ! J
           acon = 0.5_fp_kind
           bcon = -0.25_fp_kind
        else
           error stop 'exchange_gcpf: logic error 1'
        endif
        fex0 = acon*con_exchange*t*t*beta
        ! n.b. multiply by fd12(1) later
        fex1 = fex0*(fd12(1) + bcon*beta*fd32(1))
        fex1f = fex0*(fd12(2) + bcon*beta*fd32(2))*dpsidf
        fex1t = 3._fp_kind*fex1 + fex0*bcon*beta*fd32(1)
        fex1f2 = fex0&
             *((fd12(3) + bcon*beta*fd32(3))*dpsidf*dpsidf&
             + (fd12(2) + bcon*beta*fd32(2))*dpsidf2)
        fex1ft = 3._fp_kind*fex1f + fex0*bcon*beta*fd32(2)*dpsidf
        fex1t2 = 3._fp_kind*fex1t + 4._fp_kind*fex0*bcon*beta*fd32(1)
        fex1_psi = fex0*(fd12(2) + bcon*beta*fd32(2))
        fex1_psif = fex0*(fd12(3) + bcon*beta*fd32(3))*dpsidf
        fex1_psit = 3._fp_kind*fex1_psi + fex0*bcon*beta*fd32(2)
        fex1_psif2 = fex0&
             *((fd12(4) + bcon*beta*fd32(4))*dpsidf*dpsidf&
             + (fd12(3) + bcon*beta*fd32(3))*dpsidf2)
        fex1_psift = 3._fp_kind*fex1_psif&
             + fex0*bcon*beta*fd32(3)*dpsidf
        fex1_psit2 = 3._fp_kind*fex1_psit&
             + 4._fp_kind*fex0*bcon*beta*fd32(2)
        ! to get final result must multiply by fd12(1)
        ! n.b. the order of assignment statements is important here
        ! definition of fex_psi is partial of fex1 wrt psi
        ! thus from multiplication transformation we have
        ! fex1_psi = fex1_psi*fd12(1) + fex1*fd12(2)
        ! **with the rhs evaluated before any multiplication transformation
        ! is made**.
        ! apply multiplication transformation to fex1_psi
        fex1_psif2 = fex1_psif2*fd12(1)&
             + 2._fp_kind*fex1_psif*fd12(2)*dpsidf&
             + fex1_psi*(fd12(3)*dpsidf*dpsidf + fd12(2)*dpsidf2)&
             + fex1f2*fd12(2) + 2._fp_kind*fex1f*fd12(3)*dpsidf&
             + fex1*(fd12(4)*dpsidf*dpsidf + fd12(3)*dpsidf2)
        fex1_psift = fex1_psift*fd12(1) + fex1_psit*fd12(2)*dpsidf&
             + fex1ft*fd12(2) + fex1t*fd12(3)*dpsidf
        fex1_psit2 = fex1_psit2*fd12(1)&
             + fex1t2*fd12(2)
        fex1_psif = fex1_psif*fd12(1) + fex1_psi*fd12(2)*dpsidf&
             + fex1f*fd12(2) + fex1*fd12(3)*dpsidf
        fex1_psit = fex1_psit*fd12(1)&
             + fex1t*fd12(2)
        fex1_psi = fex1_psi*fd12(1)&
             + fex1*fd12(2)
        ! apply multiplication transformation to fex1
        fex1f2 = fex1f2*fd12(1) + 2._fp_kind*fex1f*fd12(2)*dpsidf&
             + fex1*(fd12(3)*dpsidf*dpsidf + fd12(2)*dpsidf2)
        fex1ft = fex1ft*fd12(1) + fex1t*fd12(2)*dpsidf
        fex1t2 = fex1t2*fd12(1)
        fex1f = fex1f*fd12(1) + fex1*fd12(2)*dpsidf
        fex1t = fex1t*fd12(1)
        fex1 = fex1*fd12(1)
        if(ifexchange.ne.21) then
           ! Add in previously calculated NR J term for J-I, and J cases.
           fex = fex + fex1
           fexf = fexf + fex1f
           fext = fext + fex1t
           fexf2 = fexf2 + fex1f2
           fexft = fexft + fex1ft
           fext2 = fext2 + fex1t2
           fex_psi = fex_psi + fex1_psi
           fex_psif = fex_psif + fex1_psif
           fex_psit = fex_psit + fex1_psit
           fex_psif2 = fex_psif2 + fex1_psif2
           fex_psift = fex_psift + fex1_psift
           fex_psit2 = fex_psit2 + fex1_psit2
        else
           ! For -I case, no NR term previously calculated.
           fex = fex1
           fexf = fex1f
           fext = fex1t
           fexf2 = fex1f2
           fexft = fex1ft
           fext2 = fex1t2
           fex_psi = fex1_psi
           fex_psif = fex1_psif
           fex_psit = fex1_psit
           fex_psif2 = fex1_psif2
           fex_psift = fex1_psift
           fex_psit2 = fex1_psit2
        endif
     endif
     if(ifexchange.eq.11.or.ifexchange.eq.41) then
        ! - K term.  CG, 24.289
        fex0 = con_exchange*con_nr1_ratio*t*t*beta*sqrt(beta)
        fex2 = fex0*(fd12(1) + 0.25_fp_kind*beta*fd32(1))
        fex2f = fex0*(fd12(2) + 0.25_fp_kind*beta*fd32(2))*dpsidf
        fex2t = 3.5_fp_kind*fex2 + fex0*0.25_fp_kind*beta*fd32(1)
        fex2f2 = fex0&
             *((fd12(3) + 0.25_fp_kind*beta*fd32(3))*dpsidf*dpsidf&
             + (fd12(2) + 0.25_fp_kind*beta*fd32(2))*dpsidf2)
        fex2ft = 3.5_fp_kind*fex2f + fex0*0.25_fp_kind*beta*fd32(2)*dpsidf
        fex2t2 = 3.5_fp_kind*fex2t + 4.5_fp_kind*fex0*0.25_fp_kind*beta*fd32(1)
        fex2_psi = fex0*(fd12(2) + 0.25_fp_kind*beta*fd32(2))
        fex2_psif = fex0*(fd12(3) + 0.25_fp_kind*beta*fd32(3))*dpsidf
        fex2_psit = 3.5_fp_kind*fex2_psi + fex0*0.25_fp_kind*beta*fd32(2)
        fex2_psif2 = fex0&
             *((fd12(4) + 0.25_fp_kind*beta*fd32(4))*dpsidf*dpsidf&
             + (fd12(3) + 0.25_fp_kind*beta*fd32(3))*dpsidf2)
        fex2_psift = 3.5_fp_kind*fex2_psif&
             + fex0*0.25_fp_kind*beta*fd32(3)*dpsidf
        fex2_psit2 = 3.5_fp_kind*fex2_psit&
             + 4.5_fp_kind*fex0*0.25_fp_kind*beta*fd32(2)
        if(ifexchange.eq.11) then
           ! Add in - K term to previously calculated J-I term.
           fex = fex + fex2
           fexf = fexf + fex2f
           fext = fext + fex2t
           fexf2 = fexf2 + fex2f2
           fexft = fexft + fex2ft
           fext2 = fext2 + fex2t2
           fex_psi = fex_psi + fex2_psi
           fex_psif = fex_psif + fex2_psif
           fex_psit = fex_psit + fex2_psit
           fex_psif2 = fex_psif2 + fex2_psif2
           fex_psift = fex_psift + fex2_psift
           fex_psit2 = fex_psit2 + fex2_psit2
        elseif(ifexchange.eq.41) then
           ! K term alone with no reference to previously
           ! calculated NR component.
           fex = fex2
           fexf = fex2f
           fext = fex2t
           fexf2 = fex2f2
           fexft = fex2ft
           fext2 = fex2t2
           fex_psi = fex2_psi
           fex_psif = fex2_psif
           fex_psit = fex2_psit
           fex_psif2 = fex2_psif2
           fex_psift = fex2_psift
           fex_psit2 = fex2_psit2
        else
           error stop 'exchange_gcpf: logic error 2'
        endif
     endif
  elseif(0.lt.ifexchange.and.mod(ifexchange,10).eq.2) then
     ! ifexchange = 2, 12, 22, 32, 42
     ! strong degeneracy
     if(psi.lt.3._fp_kind)&
          error stop 'exchange_gcpf: psi too low for degenerate series'
     z = psi*beta
     ! beta ~ T/5.930d9 so z < 0.01d0 corresponds to
     ! beta = z/psi < 0.01/psi or T < 5.930d7/psi ~ 2.d7 for psi = 3.
     if(z.lt.0.01_fp_kind) then
        write(stderr,*) 'exchange_gcpf: psi beta so low that '//&
             'significance loss occurs for degenerate series.'
        write(stderr,*) 'exchange_gcpf: use weakly relativistic '//&
             'series instead.'
        error stop 'exchange_gcpf: psi beta too low'
     endif
     phi = 1._fp_kind + z
     phiz = 1._fp_kind
     x = sqrt(2._fp_kind*z + z*z)
     xz = phi/x
     ! xz2 = 1.d0/x - phi*xz/x^2 = (x^2 - phi^2)/x^3 = -1.d0/x^3
     xz2 = -1._fp_kind/(x*x*x)
     xz3 = (-3._fp_kind*xz2/x)*xz
     lnalpha = x + phi
     lnalphaz = 1._fp_kind + xz
     lnalphaz2 = xz2
     lnalphaz3 = xz3
     ! order is important here.
     lnalphaz3 = (lnalphaz3 - 3._fp_kind*lnalphaz*lnalphaz2/lnalpha&
          + 2._fp_kind*lnalphaz*lnalphaz*lnalphaz/&
          (lnalpha*lnalpha))/lnalpha
     lnalphaz2 = (lnalphaz2 - lnalphaz*lnalphaz/lnalpha)/lnalpha
     lnalphaz = lnalphaz/lnalpha
     lnalpha = log(lnalpha)
     if(ifexchange.lt.40) then
        ! coefficient of I - J is
        ! con_d1_ratio/t^2 * coefficient of G(psi)
        fex0 = con_d1_ratio*con_exchange
        if(ifexchange.eq.2.or.ifexchange.eq.12) then
           ! I_1 - J_1 term from KLVH eq. 24.
           acon = 1.5_fp_kind
           bcon = -3.0_fp_kind
           ccon = 1.5_fp_kind
           dcon = 0.5_fp_kind
        elseif(ifexchange.eq.22) then
           ! I_1 term from KLVH eq. 10.
           acon = 0.5_fp_kind
           bcon = -1.0_fp_kind
           ccon = 0.5_fp_kind
           dcon = 0.5_fp_kind
        elseif(ifexchange.eq.32) then
           ! -J_1 term from KLVH eq. 17.
           acon = 1.0_fp_kind
           bcon = -2.0_fp_kind
           ccon = 1.0_fp_kind
           dcon = 0.0_fp_kind
        else
           error stop 'exchange_gcpf: logic error 3'
        endif
        fex = fex0*(acon*lnalpha*lnalpha&
             + x*(bcon*phi*lnalpha&
             + x*(ccon + dcon*x*x)))
        ! first store z derivatives in fex_psi, fex_psif, etc.
        fex_psi = fex0*(acon*2._fp_kind*lnalpha*lnalphaz&
             + x*bcon*(phiz*lnalpha + phi*lnalphaz)&
             + xz*(bcon*phi*lnalpha&
             + x*(2._fp_kind*ccon + 4._fp_kind*dcon*x*x)))
        fex_psif = fex0*(&
             acon*2._fp_kind*(lnalphaz*lnalphaz + lnalpha*lnalphaz2)&
             + x*bcon*(2._fp_kind*phiz*lnalphaz + phi*lnalphaz2)&
             + 2._fp_kind*xz*bcon*(phiz*lnalpha + phi*lnalphaz)&
             + xz2*(bcon*phi*lnalpha&
             + x*(2._fp_kind*ccon + 4._fp_kind*dcon*x*x))&
             + xz*xz*(2._fp_kind*ccon + 12._fp_kind*dcon*x*x))
        fex_psif2 = fex0*(&
             acon*2._fp_kind*(3._fp_kind*lnalphaz*lnalphaz2 + lnalpha*lnalphaz3)&
             + x*bcon*(3._fp_kind*phiz*lnalphaz2 + phi*lnalphaz3)&
             + 3._fp_kind*xz*bcon*(2._fp_kind*phiz*lnalphaz + phi*lnalphaz2)&
             + 3._fp_kind*xz2*bcon*(phiz*lnalpha + phi*lnalphaz)&
             + xz3*(bcon*phi*lnalpha&
             + x*(2._fp_kind*ccon + 4._fp_kind*dcon*x*x))&
             + 3._fp_kind*xz*xz2*(2._fp_kind*ccon + 12._fp_kind*dcon*x*x)&
             + xz*xz*xz*24._fp_kind*dcon*x)
        ! convert z derivatives to fl and tl derivatives
        fexf = fex_psi*beta*dpsidf
        fext = fex_psi*z
        fexf2 = fex_psif*beta*beta*dpsidf*dpsidf +&
             fex_psi*beta*dpsidf2
        fexft = fex_psif*z*beta*dpsidf + fex_psi*beta*dpsidf
        fext2 = fex_psif*z*z + fex_psi*z
        ! order of assignment statements is important
        fex_psit2 = fex_psif2*z*z*beta + 3._fp_kind*fex_psif*z*beta&
             + fex_psi*beta
        fex_psift = fex_psif2*z*beta*beta*dpsidf&
             + 2._fp_kind*fex_psif*beta*beta*dpsidf
        fex_psif2 = fex_psif2*beta*beta*beta*dpsidf*dpsidf&
             + fex_psif*beta*beta*dpsidf2
        fex_psit = fex_psif*z*beta + fex_psi*beta
        fex_psif = fex_psif*beta*beta*dpsidf
        fex_psi = fex_psi*beta
        ! note the first term above
        ! is the lowest order term in psi which for small x reduces to the
        ! equivalent of the first term in the asymptotic expression for
        ! G(psi) = 2 psi^2 + ...
        ! now do beta^2 term (which brings in psi^0 terms)
        fex0 = fex0*4._fp_kind*beta*beta
        if(ifexchange.eq.2.or.ifexchange.eq.12) then
           ! I_2 - J_2 term from KLVH eq. 24.
           acon = -1.0_fp_kind
           bcon = 4.0_fp_kind
           ccon = 1.0_fp_kind
           dcon = -3.0_fp_kind
        elseif(ifexchange.eq.22) then
           ! I_2 term from KLVH eq. 10.
           acon = 0.0_fp_kind
           bcon = 0.0_fp_kind
           ccon = 1.0_fp_kind
           dcon = -1.0_fp_kind
        elseif(ifexchange.eq.32) then
           ! -J_2 term from KLVH eq. 20 and 22.
           acon = -1.0_fp_kind
           bcon = 4.0_fp_kind
           ccon = 0.0_fp_kind
           dcon = -2.0_fp_kind
        else
           error stop 'exchange_gcpf: logic error 4'
        endif
        ! I_2-J_2, KLVH eq. 24. excluding ln beta/2) term.
        fex1 = fex0*(&
             + acon*(2._fp_kind*a1 + a2) + con2*(&
             + bcon*log(x) + ccon*(1._fp_kind + x*x) + dcon*phi*lnalpha/x))
        ! first store z derivatives in fex1_psi, fex1_psif, etc.
        fex1_psi = fex0*(&
             con2*(xz*(bcon/x + 2._fp_kind*ccon*x - dcon*phi*lnalpha/(x*x))&
             + dcon*(phiz*lnalpha + phi*lnalphaz)/x&
             ))
        fex1_psif = fex0*(&
             con2*(xz2*(bcon/x + 2._fp_kind*ccon*x - dcon*phi*lnalpha/(x*x))&
             + xz*xz*(-bcon/(x*x) + 2._fp_kind*ccon&
             + 2._fp_kind*dcon*phi*lnalpha/(x*x*x))&
             - 2._fp_kind*dcon*xz*(phiz*lnalpha + phi*lnalphaz)/(x*x)&
             + dcon*(2._fp_kind*phiz*lnalphaz + phi*lnalphaz2)/x&
             ))
        fex1_psif2 = fex0*(&
             con2*(xz3*(bcon/x + 2._fp_kind*ccon*x - dcon*phi*lnalpha/(x*x))&
             + 3._fp_kind*xz*xz2*(-bcon/(x*x) + 2._fp_kind*ccon&
             + 2._fp_kind*dcon*phi*lnalpha/(x*x*x))&
             - 3._fp_kind*dcon*xz2*(phiz*lnalpha + phi*lnalphaz)/(x*x)&
             + xz*xz*xz*(+2._fp_kind*bcon/(x*x*x)&
             - 6._fp_kind*dcon*phi*lnalpha/(x*x*x*x))&
             + 6._fp_kind*dcon*xz*xz*(phiz*lnalpha + phi*lnalphaz)/(x*x*x)&
             - 3._fp_kind*dcon*xz*(2._fp_kind*phiz*lnalphaz + phi*lnalphaz2)/(x*x)&
             + dcon*(3._fp_kind*phiz*lnalphaz2 + phi*lnalphaz3)/x&
             ))
        ! convert z derivatives to fl and tl derivatives
        fex1f = fex1_psi*beta*dpsidf
        if(ifexchange.eq.2.or.ifexchange.eq.12.or.&
             ifexchange.eq.32) then
           ! log(beta/2) term from -J_2, KLVH eq. 20.
           fex1 = fex1 + fex0*con2*(-2._fp_kind*log(beta/2._fp_kind))
           fex1t = fex1_psi*z + 2._fp_kind*fex1 - fex0*2._fp_kind*con2
        elseif(ifexchange.eq.22) then
           fex1t = fex1_psi*z + 2._fp_kind*fex1
        else
           error stop 'exchange_gcpf: logic error 5'
        endif
        fex1f2 = fex1_psif*beta*beta*dpsidf*dpsidf&
             + fex1_psi*beta*dpsidf2
        fex1ft = fex1_psif*z*beta*dpsidf&
             + 3._fp_kind*fex1_psi*beta*dpsidf
        if(ifexchange.eq.2.or.ifexchange.eq.12.or.&
             ifexchange.eq.32) then
           fex1t2 = fex1_psif*z*z + 3._fp_kind*fex1_psi*z&
                + 2._fp_kind*fex1t - 4._fp_kind*fex0*con2
        else
           fex1t2 = fex1_psif*z*z + 3._fp_kind*fex1_psi*z&
                + 2._fp_kind*fex1t
        endif
        ! order of assignment statements is important
        fex1_psit2 = fex1_psif2*z*z*beta + 7._fp_kind*fex1_psif*z*beta&
             + 9._fp_kind*fex1_psi*beta
        fex1_psift = fex1_psif2*z*beta*beta*dpsidf&
             + 4._fp_kind*fex1_psif*beta*beta*dpsidf
        fex1_psif2 = fex1_psif2*beta*beta*beta*dpsidf*dpsidf&
             + fex1_psif*beta*beta*dpsidf2
        fex1_psit = fex1_psif*z*beta + 3._fp_kind*fex1_psi*beta
        fex1_psif = fex1_psif*beta*beta*dpsidf
        fex1_psi = fex1_psi*beta
        fex = fex + fex1
        fexf = fexf + fex1f
        fext = fext + fex1t
        fexf2 = fexf2 + fex1f2
        fexft = fexft + fex1ft
        fext2 = fext2 + fex1t2
        fex_psi = fex_psi + fex1_psi
        fex_psif = fex_psif + fex1_psif
        fex_psit = fex_psit + fex1_psit
        fex_psif2 = fex_psif2 + fex1_psif2
        fex_psift = fex_psift + fex1_psift
        fex_psit2 = fex_psit2 + fex1_psit2
        ! now do beta^4 term (which brings in psi^{-2} terms)
        fex0 = fex0*beta*beta
        if(ifexchange.eq.2.or.ifexchange.eq.12) then
           ! I_3 - J_3 term from KLVH eq. 24.
           acon = 2.0_fp_kind
           bcon = 1.0_fp_kind
           ccon = 2.0_fp_kind
           dcon = 1.0_fp_kind
           econ = 3.0_fp_kind
        elseif(ifexchange.eq.22) then
           ! I_3 term from KLVH eq. 10.
           acon = 2.0_fp_kind
           bcon = 0.0_fp_kind
           ccon = -1.0_fp_kind
           dcon = -1.0_fp_kind
           econ = 1.0_fp_kind
        elseif(ifexchange.eq.32) then
           ! -J_3 term from KLVH eq. 20 and 22.
           acon = 0.0_fp_kind
           bcon = 1.0_fp_kind
           ccon = 3.0_fp_kind
           dcon = 2.0_fp_kind
           econ = 2.0_fp_kind
        else
           error stop 'exchange_gcpf: logic error 6'
        endif
        fex2 = fex0*(&
             con2*con2*(acon + (acon + bcon/(x*x))/(x*x))&
             - 3._fp_kind*con4/(x*x*x*x)*&
             (ccon + dcon*x*x + econ*phi*lnalpha/x))
        ! first store z derivatives in fex2_psi, fex2_psif, etc.
        fex2_psi = fex0*(&
             - 2._fp_kind*con2*con2*xz*(acon + 2._fp_kind*bcon/(x*x))/(x*x*x)&
             + 3._fp_kind*con4*xz/(x*x*x*x*x)&
             *(4._fp_kind*ccon + 2._fp_kind*dcon*x*x + 5._fp_kind*econ*phi*lnalpha/x)&
             - 3._fp_kind*con4/(x*x*x*x)&
             *(econ*(phiz*lnalpha + phi*lnalphaz)/x)&
             )
        fex2_psif = fex0*(&
             - 2._fp_kind*con2*con2*(xz2*(acon + 2._fp_kind*bcon/(x*x))/(x*x*x)&
             + xz*xz*(-3._fp_kind*acon - 10._fp_kind*bcon/(x*x))/(x*x*x*x))&
             + 3._fp_kind*con4*xz2/(x*x*x*x*x)&
             *(4._fp_kind*ccon + 2._fp_kind*dcon*x*x + 5._fp_kind*econ*phi*lnalpha/x)&
             + 3._fp_kind*con4*xz*xz/(x*x*x*x*x*x)&
             *(-20._fp_kind*ccon - 6._fp_kind*dcon*x*x - 30._fp_kind*econ*phi*lnalpha/x)&
             + 3._fp_kind*con4*xz/(x*x*x*x*x)&
             *(10._fp_kind*econ*(phiz*lnalpha + phi*lnalphaz)/x)&
             + 3._fp_kind*con4/(x*x*x*x)&
             *(-econ*(2._fp_kind*phiz*lnalphaz + phi*lnalphaz2)/x)&
             )
        fex2_psif2 = fex0*(&
             - 2._fp_kind*con2*con2*(xz3*(acon + 2._fp_kind*bcon/(x*x))/(x*x*x)&
             + 3._fp_kind*xz*xz2*(-3._fp_kind*acon - 10._fp_kind*bcon/(x*x))/(x*x*x*x)&
             + xz*xz*xz*(12._fp_kind*acon + 60._fp_kind*bcon/(x*x))/(x*x*x*x*x))&
             + 3._fp_kind*con4*xz3/(x*x*x*x*x)&
             *(4._fp_kind*ccon + 2._fp_kind*dcon*x*x + 5._fp_kind*econ*phi*lnalpha/x)&
             + 9._fp_kind*con4*xz*xz2/(x*x*x*x*x*x)&
             *(-20._fp_kind*ccon - 6._fp_kind*dcon*x*x - 30._fp_kind*econ*phi*lnalpha/x)&
             + 3._fp_kind*con4*xz2/(x*x*x*x*x)&
             *(15._fp_kind*econ*(phiz*lnalpha + phi*lnalphaz)/x)&
             + 3._fp_kind*con4*xz*xz/(x*x*x*x*x*x)&
             *(-90._fp_kind*econ*(phiz*lnalpha + phi*lnalphaz)/x)&
             + 3._fp_kind*con4*xz*xz*xz/(x*x*x*x*x*x*x)&
             *(+120._fp_kind*ccon +24._fp_kind*dcon*x*x + 210._fp_kind*econ*phi*lnalpha/x)&
             + 3._fp_kind*con4*xz/(x*x*x*x*x)&
             *(15._fp_kind*econ*(2._fp_kind*phiz*lnalphaz + phi*lnalphaz2)/x)&
             + 3._fp_kind*con4/(x*x*x*x)&
             *(-econ*(3._fp_kind*phiz*lnalphaz2 + phi*lnalphaz3)/x)&
             )
        ! convert z derivatives to fl and tl derivatives
        fex2f = fex2_psi*beta*dpsidf
        fex2t = fex2_psi*z + 4._fp_kind*fex2
        fex2f2 = fex2_psif*beta*beta*dpsidf*dpsidf&
             + fex2_psi*beta*dpsidf2
        fex2ft = fex2_psif*z*beta*dpsidf&
             + 5._fp_kind*fex2_psi*beta*dpsidf
        fex2t2 = fex2_psif*z*z + 5._fp_kind*fex2_psi*z&
             + 4._fp_kind*fex2t
        ! order of assignment statements is important
        fex2_psit2 = fex2_psif2*z*z*beta + 11._fp_kind*fex2_psif*z*beta&
             + 25._fp_kind*fex2_psi*beta
        fex2_psift = fex2_psif2*z*beta*beta*dpsidf&
             + 6._fp_kind*fex2_psif*beta*beta*dpsidf
        fex2_psif2 = fex2_psif2*beta*beta*beta*dpsidf*dpsidf&
             + fex2_psif*beta*beta*dpsidf2
        fex2_psit = fex2_psif*z*beta + 5._fp_kind*fex2_psi*beta
        fex2_psif = fex2_psif*beta*beta*dpsidf
        fex2_psi = fex2_psi*beta
        fex = fex + fex2
        fexf = fexf + fex2f
        fext = fext + fex2t
        fexf2 = fexf2 + fex2f2
        fexft = fexft + fex2ft
        fext2 = fext2 + fex2t2
        fex_psi = fex_psi + fex2_psi
        fex_psif = fex_psif + fex2_psif
        fex_psit = fex_psit + fex2_psit
        fex_psif2 = fex_psif2 + fex2_psif2
        fex_psift = fex_psift + fex2_psift
        fex_psit2 = fex_psit2 + fex2_psit2
     endif
     if(ifexchange.eq.12.or.ifexchange.eq.42) then
        ! add in first Kapusta term (with no KLVH equivalent)
        ! expression for strong degeneracy
        ! coefficient of F_{1/2} is
        ! fex0 = con_exchange*con_nr1_ratio*t*t*beta*sqrt(beta)
        ! but F_{1/2} = (2 beta)^{-3/2}(x*phi + ...), hence coefficient
        ! of x*phi is ...
        fex0 = con_exchange*con_nr1_ratio*t*t/(2._fp_kind*sqrt(2._fp_kind))
        fex1 = fex0*(x*phi - lnalpha)
        ! temporary to show how to deal with severe significance loss
        ! in the degenerate series for small x.  But comment out because
        ! I don't want to deal with this consistently for all series
        ! and all partial derivatives, and
        ! the weakly relativistic series can be used instead for
        ! the large degeneracy case.  In fact, because of the significance
        ! loss problem, I stop this code above now if x < 1.d-1 and an
        ! attempt is made to use the degenerate series.
        ! if(x.lt.1.d-1) then
        !   CG. eq. 24.191.
        !   fex1 = fex0*2.d0*x**3*(&
        !   1.d0/3.d0 - x*x*(&
        !   1.d0/10.d0 - x*x*(&
        !   3.d0/56.d0 - x*x*(&
        !   5.d0/144.d0 - x*x*(&
        !   35.d0/1408.d0)))))
        ! endif
        ! first store z derivatives in fex1_psi, fex1_psif, etc.
        fex1_psi = fex0*(xz*phi + x*phiz - lnalphaz)
        fex1_psif = fex0*(xz2*phi + 2._fp_kind*xz*phiz - lnalphaz2)
        fex1_psif2 = fex0*(xz3*phi + 3._fp_kind*xz2*phiz - lnalphaz3)
        ! convert z derivatives to fl and tl derivatives
        ! copied all conversions from above because that fex0 term was also
        ! just proportional to t**2 (excluding con2 term).
        fex1f = fex1_psi*beta*dpsidf
        fex1t = fex1_psi*z + 2._fp_kind*fex1
        fex1f2 = fex1_psif*beta*beta*dpsidf*dpsidf&
             + fex1_psi*beta*dpsidf2
        fex1ft = fex1_psif*z*beta*dpsidf&
             + 3._fp_kind*fex1_psi*beta*dpsidf
        fex1t2 = fex1_psif*z*z + 3._fp_kind*fex1_psi*z&
             + 2._fp_kind*fex1t
        ! order of assignment statements is important
        fex1_psit2 = fex1_psif2*z*z*beta + 7._fp_kind*fex1_psif*z*beta&
             + 9._fp_kind*fex1_psi*beta
        fex1_psift = fex1_psif2*z*beta*beta*dpsidf&
             + 4._fp_kind*fex1_psif*beta*beta*dpsidf
        fex1_psif2 = fex1_psif2*beta*beta*beta*dpsidf*dpsidf&
             + fex1_psif*beta*beta*dpsidf2
        fex1_psit = fex1_psif*z*beta + 3._fp_kind*fex1_psi*beta
        fex1_psif = fex1_psif*beta*beta*dpsidf
        fex1_psi = fex1_psi*beta
        ! next higher term in beta
        fex0 = fex0*4._fp_kind*con2*beta*beta
        fex2 = fex0*phi/x
        ! first store z derivatives in fex1_psi, fex1_psif, etc.
        fex2_psi = fex0*(phiz - phi*xz/x)/x
        fex2_psif = fex0*(-2._fp_kind*phiz*xz - phi*xz2&
             + 2._fp_kind*phi*xz*xz/x)/(x*x)
        fex2_psif2 = fex0*(3._fp_kind*phiz*(-xz2 + 2._fp_kind*xz*xz/x)&
             + phi*(-xz3 + 6._fp_kind*xz*(xz2 - xz*xz/x)/x))/(x*x)
        ! convert z derivatives to fl and tl derivatives
        fex2f = fex2_psi*beta*dpsidf
        fex2t = fex2_psi*z + 4._fp_kind*fex2
        fex2f2 = fex2_psif*beta*beta*dpsidf*dpsidf&
             + fex2_psi*beta*dpsidf2
        fex2ft = fex2_psif*z*beta*dpsidf&
             + 5._fp_kind*fex2_psi*beta*dpsidf
        fex2t2 = fex2_psif*z*z + 5._fp_kind*fex2_psi*z&
             + 4._fp_kind*fex2t
        ! order of assignment statements is important
        fex2_psit2 = fex2_psif2*z*z*beta + 11._fp_kind*fex2_psif*z*beta&
             + 25._fp_kind*fex2_psi*beta
        fex2_psift = fex2_psif2*z*beta*beta*dpsidf&
             + 6._fp_kind*fex2_psif*beta*beta*dpsidf
        fex2_psif2 = fex2_psif2*beta*beta*beta*dpsidf*dpsidf&
             + fex2_psif*beta*beta*dpsidf2
        fex2_psit = fex2_psif*z*beta + 5._fp_kind*fex2_psi*beta
        fex2_psif = fex2_psif*beta*beta*dpsidf
        fex2_psi = fex2_psi*beta
        ! add in lower order term
        fex2 = fex2 + fex1
        fex2f = fex2f + fex1f
        fex2t = fex2t + fex1t
        fex2f2 = fex2f2 + fex1f2
        fex2ft = fex2ft + fex1ft
        fex2t2 = fex2t2 + fex1t2
        fex2_psi = fex2_psi + fex1_psi
        fex2_psif = fex2_psif + fex1_psif
        fex2_psit = fex2_psit + fex1_psit
        fex2_psif2 = fex2_psif2 + fex1_psif2
        fex2_psift = fex2_psift + fex1_psift
        fex2_psit2 = fex2_psit2 + fex1_psit2
        ! temporary.  Add more K terms to show how to increase precision
        ! at large degeneracy at expense of precision for low degeneracy
        ! for low psi beta.  However, this low psi beta region is handled
        ! by weakly relativistic series in any case because of
        ! significance loss in some of the expressions so leave all this
        ! commented out now that I have proved the point with some plots.
        ! fex0 = fex0/con2*beta*beta
        ! fex2 = fex2 + fex0*(3.d0*con4*phi/x**5 +&
        !   15.d0*con6*beta*beta*phi*(7.d0 + 4.d0*x*x)/x**9)
        if(ifexchange.eq.12) then
           fex = fex + fex2
           fexf = fexf + fex2f
           fext = fext + fex2t
           fexf2 = fexf2 + fex2f2
           fexft = fexft + fex2ft
           fext2 = fext2 + fex2t2
           fex_psi = fex_psi + fex2_psi
           fex_psif = fex_psif + fex2_psif
           fex_psit = fex_psit + fex2_psit
           fex_psif2 = fex_psif2 + fex2_psif2
           fex_psift = fex_psift + fex2_psift
           fex_psit2 = fex_psit2 + fex2_psit2
        elseif(ifexchange.eq.42) then
           ! test mode for K term alone.
           fex = fex2
           fexf = fex2f
           fext = fex2t
           fexf2 = fex2f2
           fexft = fex2ft
           fext2 = fex2t2
           fex_psi = fex2_psi
           fex_psif = fex2_psif
           fex_psit = fex2_psit
           fex_psif2 = fex2_psif2
           fex_psift = fex2_psift
           fex_psit2 = fex2_psit2
        else
           error stop 'exchange_gcpf: logic error 7'
        endif
     endif
  elseif(4.le.abs(mod(ifexchange,10)).and.&
       abs(mod(ifexchange,10)).le.6) then
     ! Must provide fexi in these cases.
     iffexi = (0.le.abs(ifexchange).and.abs(ifexchange).lt.30)
     ! Must provide fexj in all cases except fexi alone or fexk alone.
     iffexj = .not.(&
          (20.le.abs(ifexchange).and.abs(ifexchange).lt.30).or.&
          (40.le.abs(ifexchange).and.abs(ifexchange).lt.50))
     ! Must provide fexk in these cases.
     iffexk = (10.le.abs(ifexchange).and.abs(ifexchange).lt.20).or.&
          (40.le.abs(ifexchange).and.abs(ifexchange).lt.50)

     if(index_exchange.ne.mod(abs(ifexchange),10)-3) then
        index_exchange = mod(abs(ifexchange),10)-3
        ! fill ccoeff, mforder, and mgorder if that hasn't been done before.
        call exchange_coeff(index_exchange, ccoeff, mforder, mgorder)
        rmforder = real(mforder, fp_kind)
        rmgorder = real(mgorder, fp_kind)
     endif
     fexi = 0._fp_kind
     fexif = 0._fp_kind
     fexit = 0._fp_kind
     fexif2 = 0._fp_kind
     fexift = 0._fp_kind
     fexit2 = 0._fp_kind
     fexi_psi = 0._fp_kind
     fexi_psif = 0._fp_kind
     fexi_psit = 0._fp_kind
     fexi_psif2 = 0._fp_kind
     fexi_psift = 0._fp_kind
     fexi_psit2 = 0._fp_kind
     fexj = 0._fp_kind
     fexjf = 0._fp_kind
     fexjt = 0._fp_kind
     fexjf2 = 0._fp_kind
     fexjft = 0._fp_kind
     fexjt2 = 0._fp_kind
     fexj_psi = 0._fp_kind
     fexj_psif = 0._fp_kind
     fexj_psit = 0._fp_kind
     fexj_psif2 = 0._fp_kind
     fexj_psift = 0._fp_kind
     fexj_psit2 = 0._fp_kind
     fexk = 0._fp_kind
     fexkf = 0._fp_kind
     fexkt = 0._fp_kind
     fexkf2 = 0._fp_kind
     fexkft = 0._fp_kind
     fexkt2 = 0._fp_kind
     fexk_psi = 0._fp_kind
     fexk_psif = 0._fp_kind
     fexk_psit = 0._fp_kind
     fexk_psif2 = 0._fp_kind
     fexk_psift = 0._fp_kind
     fexk_psit2 = 0._fp_kind
     wf = sqrt(1._fp_kind+f)
     vf = 1._fp_kind/(1._fp_kind+f)
     uf = f*vf
     ! first and second derivatives of uf wrt ln f.  N.B. 1-uf = (1+f-f)/(1+f) = 1/(1+f).
     ! Note the limits of uf, duf, and duf2 as f approaches 0 are all f.
     ! And the limits of these quantities as f approaches infinity are respectively 1, 1/f, and -1/f
     duf = uf/(1._fp_kind+f)
     duf2 = duf*(1._fp_kind-f)/(1._fp_kind+f)
     if(mod(ifexchange,10).gt.0) then
        g = beta*wf
     else
        g = 0._fp_kind
        if(iffexj) then
           ! if fexj is provided it contains the NR component.
           iffexi = .false.
           iffexk = .false.
        endif
     endif
     vg = 1._fp_kind/(1._fp_kind+g)
     ug = g*vg
     ! first and second derivatives of ug wrt ln g.  Since ug is the same function of
     ! g that uf is of f, follow exactly what was implemented above for duf and duf2
     dug = ug/(1._fp_kind+g)
     dug2 = dug*(1._fp_kind-g)/(1._fp_kind+g)

     if(iffexj) then
        ! Calculate approximation to J integral.
        ! From Paper IV, eqs. 13 and 32
        ! NR component = 2 beta^2 G(eta)/(1+g)^{N_max} =
        !   2 pi beta^2 fpsi/(1+g)^{N_max} =
        !   2 pi g^2/(1+f) fpsi/(1+g)^{N_max},
        ! where N_max is the maximum of N+1, N'+1, N''+3, N'''+1.
        call f_psi(psi, fpsi, dfpsi, d2fpsi, d3fpsi)
        dfpsidf = dfpsi*dpsidf
        dfpsidf2 = d2fpsi*dpsidf*dpsidf + dfpsi*dpsidf2
        ddfpsidf = d2fpsi*dpsidf
        ddfpsidf2 = d3fpsi*dpsidf*dpsidf + d2fpsi*dpsidf2
        n_max = max(mgorder(1)+1, mgorder(2)+1,mgorder(3)+3, mgorder(4)+1)
        rn_max = real(n_max,fp_kind)
        ! From Paper IV, eqs. 13 and 32 and G(eta) = pi*fpsi
        fex0 = 2._fp_kind*pi*beta*beta*vg**n_max
        fexj = fex0*fpsi
        fexjf = fex0*dfpsidf + fexj*(-0.5_fp_kind*rn_max)*ug*uf
        fexjt = fexj*(2._fp_kind - rn_max*ug)
        fexjf2 = fex0*dfpsidf2 +&
             fex0*dfpsidf*(-0.5_fp_kind*rn_max)*ug*uf +&
             fexjf*(-0.5_fp_kind*rn_max)*ug*uf +&
             fexj*(-0.5_fp_kind*rn_max)*(dug*uf*0.5_fp_kind*uf + ug*duf)
        fexjft = fexjf*(2._fp_kind - rn_max*ug) +&
             fexj*(-0.5_fp_kind*rn_max)*(dug*uf)
        fexjt2 = fexjt*(2._fp_kind - rn_max*ug) - fexj*rn_max*dug
        fexj_psi = fex0*dfpsi
        fexj_psif = fexj_psi*(-0.5_fp_kind*rn_max)*ug*uf +&
             fex0*ddfpsidf
        fexj_psit = fexj_psi*(2._fp_kind - rn_max*ug)
        fexj_psif2 = fexj_psi*(-0.5_fp_kind*rn_max)*&
             (dug*uf*0.5_fp_kind*uf + ug*duf) +&
             fexj_psif*(-0.5_fp_kind*rn_max)*ug*uf +&
             fex0*ddfpsidf*(-0.5_fp_kind*rn_max)*ug*uf +&
             fex0*ddfpsidf2
        fexj_psift = fexj_psif*(2._fp_kind - rn_max*ug) +&
             fexj_psi*(-rn_max*dug)*0.5_fp_kind*uf
        fexj_psit2 = fexj_psit*(2._fp_kind - rn_max*ug) +&
             fexj_psi*(-rn_max*dug)
        fexj_psi = fexj_psi +&
             fexj*(-0.5_fp_kind*rn_max)*ug*uf/dpsidf
        fexj_psif = fexj_psif +&
             fexjf*(-0.5_fp_kind*rn_max)*ug*uf/dpsidf +&
             fexj*(-0.5_fp_kind*rn_max)*&
             (dug*uf*0.5_fp_kind*uf + ug*duf - ug*uf*dpsidf2/dpsidf)/dpsidf
        fexj_psit = fexj_psit +&
             fexjt*(-0.5_fp_kind*rn_max)*ug*uf/dpsidf +&
             fexj*(-0.5_fp_kind*rn_max)*dug*uf/dpsidf
        fexj_psif2 = fexj_psif2 +&
             fexjf2*(-0.5_fp_kind*rn_max)*ug*uf/dpsidf +&
             2._fp_kind*fexjf*(-0.5_fp_kind*rn_max)*&
             (dug*uf*0.5_fp_kind*uf + ug*duf - ug*uf*dpsidf2/dpsidf)/dpsidf +&
             fexj*(-0.5_fp_kind*rn_max)*&
             (dug2*uf*0.5_fp_kind*uf*0.5_fp_kind*uf + 1.5_fp_kind*dug*uf*duf + ug*duf2 -&
             2._fp_kind*(dug*uf*0.5_fp_kind*uf + ug*duf)*dpsidf2/dpsidf -&
             ug*uf*(dpsidf3 - 2._fp_kind*dpsidf2*dpsidf2/dpsidf)/dpsidf&
             )/dpsidf
        fexj_psift = fexj_psift +&
             fexjft*(-0.5_fp_kind*rn_max)*ug*uf/dpsidf +&
             fexjt*(-0.5_fp_kind*rn_max)*&
             (dug*uf*0.5_fp_kind*uf + ug*duf - ug*uf*dpsidf2/dpsidf)/dpsidf +&
             fexjf*(-0.5_fp_kind*rn_max)*dug*uf/dpsidf +&
             fexj*(-0.5_fp_kind*rn_max)*&
             (dug2*uf*0.5_fp_kind*uf + dug*duf - dug*uf*dpsidf2/dpsidf&
             )/dpsidf
        fexj_psit2 = fexj_psit2 +&
             fexjt2*(-0.5_fp_kind*rn_max)*ug*uf/dpsidf +&
             2._fp_kind*fexjt*(-0.5_fp_kind*rn_max)*dug*uf/dpsidf +&
             fexj*(-0.5_fp_kind*rn_max)*dug2*uf/dpsidf
        ! fdf = (f/(1+f))^2 g^2 (g/(1+g))
        fdf = uf*uf*g*g*ug
        do ikind = 1,4
           if(mforder(ikind).ge.0.and.mgorder(ikind).ge.0) then
              ! Note, for negative mod(ifexchange,10), g is zero.
              ! Nevertheless, for programming
              ! simplicity effsum_calc grinds through entire sum.
              call effsum_calc(f, g,&
                   reshape(ccoeff(1:(mforder(ikind)+1)*(mgorder(ikind)+1),ikind), [mforder(ikind)+1, mgorder(ikind)+1]), sum)
              ! from Paper IV, equation 32
              fex0 = fdf*vf**mforder(ikind)*vg**mgorder(ikind)
              ! ln fex0 = 2 ln f - 2 ln 1+f + 3 ln g - ln 1+g
              ! -mforder ln 1+f -mgorder ln 1+g
              fex0f = fex0*(2._fp_kind - (2._fp_kind + rmforder(ikind))*uf)
              fex0g = fex0*(3._fp_kind - (1._fp_kind + rmgorder(ikind))*ug)
              fex0ff = fex0f*(2._fp_kind - (2._fp_kind + rmforder(ikind))*uf) +&
                   fex0*(-(2._fp_kind + rmforder(ikind)))*duf
              fex0fg = fex0g*(2._fp_kind - (2._fp_kind + rmforder(ikind))*uf)
              fex0gg = fex0g*(3._fp_kind - (1._fp_kind + rmgorder(ikind))*ug) +&
                   fex0*(-(1._fp_kind + rmgorder(ikind)))*dug
              fex0fff = fex0ff*(2._fp_kind - (2._fp_kind + rmforder(ikind))*uf) +&
                   2._fp_kind*(fex0f*(-(2._fp_kind + rmforder(ikind)))*duf) +&
                   fex0*(-(2._fp_kind + rmforder(ikind)))*duf2
              fex0ffg = fex0fg*(2._fp_kind - (2._fp_kind + rmforder(ikind))*uf) +&
                   fex0g*(-(2._fp_kind + rmforder(ikind))*duf)
              fex0fgg = fex0gg*(2._fp_kind - (2._fp_kind + rmforder(ikind))*uf)
              fex0ggg = fex0gg*(3._fp_kind - (1._fp_kind + rmgorder(ikind))*ug) +&
                   2._fp_kind*(fex0g*(-(1._fp_kind + rmgorder(ikind)))*dug) +&
                   fex0*(-(1._fp_kind + rmgorder(ikind)))*dug2
              if(ikind.eq.2) then
                 if(g.gt.1.e-2_fp_kind) then
                    fex0log = log(1._fp_kind+g)
                 else
                    fex0log = g*(1._fp_kind - g*(1._fp_kind/2._fp_kind - g*(1._fp_kind/3._fp_kind -&
                         g*(1._fp_kind/4._fp_kind - g*(1._fp_kind/5._fp_kind - g*(1._fp_kind/6._fp_kind -&
                         g*(1._fp_kind/7._fp_kind - g*(1._fp_kind/8._fp_kind))))))))
                 endif
                 fex0ggg = fex0log*fex0ggg + 3._fp_kind*ug*fex0gg +&
                      3._fp_kind*dug*fex0g + dug2*fex0
                 fex0fgg = fex0log*fex0fgg + 2._fp_kind*ug*fex0fg +&
                      dug*fex0f
                 fex0ffg = fex0log*fex0ffg + ug*fex0ff
                 fex0fff = fex0log*fex0fff
                 fex0gg = fex0log*fex0gg + 2._fp_kind*ug*fex0g + dug*fex0
                 fex0fg = fex0log*fex0fg + ug*fex0f
                 fex0ff = fex0log*fex0ff
                 fex0g = fex0log*fex0g + ug*fex0
                 fex0f = fex0log*fex0f
                 fex0 = fex0log*fex0
              elseif(ikind.eq.3) then
                 ! N.B. this branch normally not used, but I have tested
                 ! it anyhow with exchange_test and suitable swapping
                 ! of ikind=2 and 3.  There was substantial significance
                 ! loss, but I think all derivatives were okay, but should
                 ! retest if ever seriously use this branch.
                 if(g.gt.1.e-2_fp_kind) then
                    fex0log = log(1._fp_kind+g)
                 else
                    fex0log = g*(1._fp_kind - g*(1._fp_kind/2._fp_kind - g*(1._fp_kind/3._fp_kind -&
                         g*(1._fp_kind/4._fp_kind - g*(1._fp_kind/5._fp_kind - g*(1._fp_kind/6._fp_kind -&
                         g*(1._fp_kind/7._fp_kind - g*(1._fp_kind/8._fp_kind))))))))
                 endif
                 fex0log2 = -fex0log*fex0log*vg*vg
                 d1fex0log2 = -2._fp_kind*(&
                      ug*(fex0log*vg*vg + fex0log2)&
                      )
                 d2fex0log2 = -2._fp_kind*(&
                      dug*(fex0log*vg*vg + fex0log2) +&
                      ug*(ug*vg*vg - 2._fp_kind*fex0log*vg*vg*ug + d1fex0log2)&
                      )
                 d3fex0log2 = -2._fp_kind*(&
                      dug2*(fex0log*vg*vg + fex0log2) +&
                      2._fp_kind*dug*(ug*vg*vg - 2._fp_kind*fex0log*vg*vg*ug +&
                      d1fex0log2) +&
                      ug*(dug*vg*vg - 4._fp_kind*ug*vg*vg*ug -&
                      2._fp_kind*fex0log*(-2._fp_kind*vg*vg*ug*ug + vg*vg*dug) +&
                      d2fex0log2)&
                      )

                 fex0ggg = d3fex0log2*fex0 + 3._fp_kind*d2fex0log2*fex0g +&
                      3._fp_kind*d1fex0log2*fex0gg + fex0log2*fex0ggg
                 fex0fgg = d2fex0log2*fex0f + 2._fp_kind*d1fex0log2*fex0fg +&
                      fex0log2*fex0fgg
                 fex0ffg = d1fex0log2*fex0ff + fex0log2*fex0ffg
                 fex0fff = fex0log2*fex0fff
                 fex0gg = d2fex0log2*fex0 + 2._fp_kind*d1fex0log2*fex0g +&
                      fex0log2*fex0gg
                 fex0fg = d1fex0log2*fex0f + fex0log2*fex0fg
                 fex0ff = fex0log2*fex0ff
                 fex0g = d1fex0log2*fex0 + fex0log2*fex0g
                 fex0f = fex0log2*fex0f
                 fex0 = fex0log2*fex0
              elseif(ikind.eq.4) then
                 if(f.gt.1.e-2_fp_kind) then
                    fex0log = log(1._fp_kind+f)
                 else
                    fex0log = f*(1._fp_kind - f*(1._fp_kind/2._fp_kind - f*(1._fp_kind/3._fp_kind -&
                         f*(1._fp_kind/4._fp_kind - f*(1._fp_kind/5._fp_kind - f*(1._fp_kind/6._fp_kind -&
                         f*(1._fp_kind/7._fp_kind - f*(1._fp_kind/8._fp_kind))))))))
                 endif
                 ! The limits of fex0log as f approaches zero is f and
                 ! the same is true for uf, duf, and duf2, see above.
                 ! Therefore, near but still above the underflow limit
                 ! for f, there is no chance of gitting a divide by
                 ! zero with the following implementation which is
                 ! a great improvement on the previous implementation which
                 ! used factors of (1/fex0log) cubed!
                 r_uf = uf/fex0log
                 r_duf = duf/fex0log
                 r_duf2 = duf2/fex0log
                 fex0log2 = -fex0log/(1._fp_kind+f)
                 d1fex0log2 = fex0log2*(r_uf-uf)
                 d2fex0log2 = d1fex0log2*(r_uf-uf) + fex0log2*((r_duf-duf) - r_uf*r_uf)
                 d3fex0log2 = d2fex0log2*(r_uf-uf) + 2._fp_kind*d1fex0log2*((r_duf-duf) - r_uf*r_uf) +&
                      fex0log2*((r_duf2-duf2) - 3._fp_kind*r_duf*r_uf + 2._fp_kind*r_uf*r_uf*r_uf)

                 fex0ggg = fex0log2*fex0ggg
                 fex0fgg = d1fex0log2*fex0gg + fex0log2*fex0fgg
                 fex0ffg = d2fex0log2*fex0g + 2._fp_kind*d1fex0log2*fex0fg + fex0log2*fex0ffg
                 fex0fff = d3fex0log2*fex0 + 3._fp_kind*d2fex0log2*fex0f + 3._fp_kind*d1fex0log2*fex0ff + fex0log2*fex0fff
                 fex0gg = fex0log2*fex0gg
                 fex0fg = d1fex0log2*fex0g + fex0log2*fex0fg
                 fex0ff = d2fex0log2*fex0 + 2._fp_kind*d1fex0log2*fex0f + fex0log2*fex0ff
                 fex0g = fex0log2*fex0g
                 fex0f = d1fex0log2*fex0 + fex0log2*fex0f
                 fex0 = fex0log2*fex0
              endif
              ! reverse order to preserve values on RHS.
              ! ggg
              sum(10) = fex0ggg*sum(1) + 3._fp_kind*fex0gg*sum(3) +&
                   3._fp_kind*fex0g*sum(6) + fex0*sum(10)
              ! fgg
              sum(9) = fex0fgg*sum(1) + fex0gg*sum(2) +&
                   2._fp_kind*(fex0fg*sum(3) + fex0g*sum(5)) +&
                   fex0f*sum(6) + fex0*sum(9)
              ! ffg
              sum(8) = fex0ffg*sum(1) + fex0ff*sum(3) +&
                   2._fp_kind*(fex0fg*sum(2) + fex0f*sum(5)) +&
                   fex0g*sum(4) + fex0*sum(8)
              ! fff
              sum(7) = fex0fff*sum(1) + 3._fp_kind*fex0ff*sum(2) +&
                   3._fp_kind*fex0f*sum(4) + fex0*sum(7)
              ! gg
              sum(6) = fex0gg*sum(1) + 2._fp_kind*fex0g*sum(3) + fex0*sum(6)
              ! fg
              sum(5) = fex0fg*sum(1) + fex0g*sum(2) + fex0f*sum(3) +&
                   fex0*sum(5)
              ! ff
              sum(4) = fex0ff*sum(1) + 2._fp_kind*fex0f*sum(2) + fex0*sum(4)
              ! g
              sum(3) = fex0g*sum(1) + fex0*sum(3)
              ! f
              sum(2) = fex0f*sum(1) + fex0*sum(2)
              sum(1) = fex0*sum(1)
              ! Transform to tcl derivative from ln g
              ! derivative recalling that dlng/dlntc = 1, and dlng/dlnf = 0.5*uf.
              sum(4) = sum(4) +&
                   0.5_fp_kind*(duf*sum(3) + uf*sum(5)) +&
                   0.5_fp_kind*uf*(sum(5) + 0.5_fp_kind*uf*sum(6))
              sum(7) = sum(7) +&
                   0.5_fp_kind*(duf2*sum(3) + 2._fp_kind*duf*sum(5) + uf*sum(8)) +&
                   0.5_fp_kind*(duf*sum(5) + uf*sum(8) +&
                   0.5_fp_kind*uf*(2._fp_kind*duf*sum(6) + uf*sum(9))) +&
                   0.5_fp_kind*uf*(sum(8) +&
                   0.5_fp_kind*(duf*sum(6) + uf*sum(9)) +&
                   0.5_fp_kind*uf*(sum(9) + 0.5_fp_kind*uf*sum(10)))
              ! do sum(5) later than sum(7) so that rhs of sum(7) is undisturbed.
              sum(5) = sum(5) + 0.5_fp_kind*uf*sum(6)
              sum(8) = sum(8) +&
                   0.5_fp_kind*(duf*sum(6) + uf*sum(9)) +&
                   0.5_fp_kind*uf*(sum(9) + 0.5_fp_kind*uf*sum(10))
              sum(9) = sum(9) + 0.5_fp_kind*uf*sum(10)
              sum(2) = sum(2) + 0.5_fp_kind*uf*sum(3)
              fexj = fexj + sum(1)
              fexjf = fexjf + sum(2)
              fexjt = fexjt + sum(3)
              fexjf2 = fexjf2 + sum(4)
              fexjft = fexjft + sum(5)
              fexjt2 = fexjt2 + sum(6)
              fexj_psi = fexj_psi + sum(2)/dpsidf
              fexj_psif = fexj_psif +&
                   (sum(4) - sum(2)*dpsidf2/dpsidf)/dpsidf
              fexj_psit = fexj_psit + sum(5)/dpsidf
              fexj_psif2 = fexj_psif2 +&
                   (sum(7) - (2._fp_kind*sum(4)*dpsidf2 + sum(2)*dpsidf3 -&
                   2._fp_kind*sum(2)*dpsidf2*dpsidf2/dpsidf)/dpsidf)/dpsidf
              fexj_psift = fexj_psift +&
                   (sum(8) - sum(5)*dpsidf2/dpsidf)/dpsidf
              fexj_psit2 = fexj_psit2 + sum(9)/dpsidf
           endif
        enddo
     endif
     if(iffexi.or.iffexk) then
        ! Must calculate I=K^2 integral and/or K integral.
        ! Calculate approximation to K integral.
        ! From Paper IV, eqs. 17 and 37
        ! NR component = 2 beta^{3/2} fd(1)/(1+g)^{N + 0.5}
        do ideriv = 1, maxderiv
           ! F(1/2) and 3 derivatives
           fd(ideriv) = (2._fp_kind/3._fp_kind)*fermi_dirac_ct(psi,ideriv,0)
        enddo
        ! Order is important.
        dfd(1) = fd(2)
        dfd(2) = fd(3)*dpsidf
        dfd(3) = fd(4)*dpsidf*dpsidf + fd(3)*dpsidf2
        fd(3) = fd(3)*dpsidf*dpsidf + fd(2)*dpsidf2
        fd(2) = fd(2)*dpsidf
        power = real(max(mgorder(5), mgorder(6)),fp_kind) + 0.5_fp_kind
        fex0 = 2._fp_kind*beta*sqrt(beta)*vg**power

        fexl = fex0*fd(1)
        fexlf = fex0*fd(2) + fexl*(-0.5_fp_kind*power)*ug*uf
        fexlt = fexl*(1.5_fp_kind - power*ug)
        fexlf2 = fex0*fd(3) +&
             fex0*fd(2)*(-0.5_fp_kind*power)*ug*uf +&
             fexlf*(-0.5_fp_kind*power)*ug*uf +&
             fexl*(-0.5_fp_kind*power)*(dug*uf*0.5_fp_kind*uf + ug*duf)
        fexlft = fexlf*(1.5_fp_kind - power*ug) +&
             fexl*(-0.5_fp_kind*power)*(dug*uf)
        fexlt2 = fexlt*(1.5_fp_kind - power*ug) - fexl*power*dug

        fexl_psi = fex0*dfd(1)
        fexl_psif = fex0*dfd(2) + fexl_psi*(-0.5_fp_kind*power)*ug*uf
        fexl_psit = fexl_psi*(1.5_fp_kind - power*ug)
        fexl_psif2 = fex0*dfd(3) +&
             fex0*dfd(2)*(-0.5_fp_kind*power)*ug*uf +&
             fexl_psif*(-0.5_fp_kind*power)*ug*uf +&
             fexl_psi*(-0.5_fp_kind*power)*(dug*uf*0.5_fp_kind*uf + ug*duf)
        fexl_psift = fexl_psif*(1.5_fp_kind - power*ug) +&
             fexl_psi*(-0.5_fp_kind*power)*(dug*uf)
        fexl_psit2 = fexl_psit*(1.5_fp_kind - power*ug) -&
             fexl_psi*power*dug

        fexl_psi = fexl_psi +&
             fexl*(-0.5_fp_kind*power)*ug*uf/dpsidf
        fexl_psif = fexl_psif + (&
             fexlf*(-0.5_fp_kind*power)*ug*uf +&
             fexl*(-0.5_fp_kind*power)*(&
             dug*uf*0.5_fp_kind*uf + ug*duf - ug*uf*dpsidf2/dpsidf)&
             )/dpsidf
        fexl_psit = fexl_psit +&
             fexlt*(-0.5_fp_kind*power)*ug*uf/dpsidf +&
             fexl*(-0.5_fp_kind*power)*dug*uf/dpsidf
        fexl_psif2 = fexl_psif2 + (&
             fexlf2*(-0.5_fp_kind*power)*ug*uf +&
             fexlf*(-0.5_fp_kind*power)*(&
             dug*uf*uf + 2._fp_kind*ug*duf - ug*uf*dpsidf2/dpsidf) +&
             fexl*(-0.5_fp_kind*power)*(&
             dug2*uf*0.5_fp_kind*uf*0.5_fp_kind*uf + 1.5_fp_kind*dug*uf*duf + ug*duf2 - (&
             dug*uf*dpsidf2*0.5_fp_kind*uf + ug*duf*dpsidf2 + ug*uf*dpsidf3 -&
             ug*uf*dpsidf2*dpsidf2/dpsidf&
             )/dpsidf&
             ) - dpsidf2*(&
             fexlf*(-0.5_fp_kind*power)*ug*uf +&
             fexl*(-0.5_fp_kind*power)*(&
             dug*uf*0.5_fp_kind*uf + ug*duf - ug*uf*dpsidf2/dpsidf)&
             )/dpsidf&
             )/dpsidf
        fexl_psift = fexl_psift + (&
             fexlft*(-0.5_fp_kind*power)*ug*uf +&
             fexlf*(-0.5_fp_kind*power)*dug*uf +&
             fexlt*(-0.5_fp_kind*power)*(&
             dug*uf*0.5_fp_kind*uf + ug*duf - ug*uf*dpsidf2/dpsidf) +&
             fexl*(-0.5_fp_kind*power)*(&
             dug2*uf*0.5_fp_kind*uf + dug*duf - dug*uf*dpsidf2/dpsidf)&
             )/dpsidf
        fexl_psit2 = fexl_psit2 +&
             fexlt2*(-0.5_fp_kind*power)*ug*uf/dpsidf +&
             2._fp_kind*fexlt*(-0.5_fp_kind*power)*dug*uf/dpsidf +&
             fexl*(-0.5_fp_kind*power)*dug2*uf/dpsidf

        ! fdf = f/(1+f) g^2 sqrt(g/(1+g))
        fdf = uf*g*g*sqrt(ug)
        do ikind = 5, 6
           if(mforder(ikind).ge.0.and.mgorder(ikind).ge.0) then
              ! Note, for negative mod(ifexchange,10), g is zero.
              ! Nevertheless, for programming
              ! simplicity effsum_calc grinds through entire sum.
              call effsum_calc(f, g,&
                   reshape(ccoeff(1:(mforder(ikind)+1)*(mgorder(ikind)+1),ikind), [mforder(ikind)+1, mgorder(ikind)+1]), sum)
              ! from Paper IV, equation 37
              fex0 = fdf*vf**mforder(ikind)*vg**mgorder(ikind)
              ! ln fex0 = ln f + 2.5 ln g
              ! -(mforder+1) ln 1+f -(mgorder+0.5) ln 1+g
              fex0f = fex0*(1._fp_kind - (1._fp_kind + rmforder(ikind))*uf)
              fex0g = fex0*(2.5_fp_kind - (0.5_fp_kind + rmgorder(ikind))*ug)
              fex0ff = fex0f*(1._fp_kind - (1._fp_kind + rmforder(ikind))*uf) +&
                   fex0*(-(1._fp_kind + rmforder(ikind)))*duf
              fex0fg = fex0g*(1._fp_kind - (1._fp_kind + rmforder(ikind))*uf)
              fex0gg = fex0g*(2.5_fp_kind - (0.5_fp_kind + rmgorder(ikind))*ug) +&
                   fex0*(-(0.5_fp_kind + rmgorder(ikind)))*dug
              fex0fff = fex0ff*(1._fp_kind - (1._fp_kind + rmforder(ikind))*uf) +&
                   2._fp_kind*(fex0f*(-(1._fp_kind + rmforder(ikind)))*duf) +&
                   fex0*(-(1._fp_kind + rmforder(ikind)))*duf2
              fex0ffg = fex0fg*(1._fp_kind - (1._fp_kind + rmforder(ikind))*uf) +&
                   fex0g*(-(1._fp_kind + rmforder(ikind))*duf)
              fex0fgg = fex0gg*(1._fp_kind - (1._fp_kind + rmforder(ikind))*uf)
              fex0ggg = fex0gg*(2.5_fp_kind - (0.5_fp_kind + rmgorder(ikind))*ug) +&
                   2._fp_kind*(fex0g*(-(0.5_fp_kind + rmgorder(ikind)))*dug) +&
                   fex0*(-(0.5_fp_kind + rmgorder(ikind)))*dug2

              if(ikind.eq.6) then
                 ! N.B. there seems to be some significance loss here
                 ! that I cannot track down when looked at in isolation.
                 ! But it doesn't matter when combined in the do ikind loop.
                 if(g.gt.1.e-2_fp_kind) then
                    fex0log = log(1._fp_kind+g)
                 else
                    fex0log = g*(1._fp_kind - g*(1._fp_kind/2._fp_kind - g*(1._fp_kind/3._fp_kind -&
                         g*(1._fp_kind/4._fp_kind - g*(1._fp_kind/5._fp_kind - g*(1._fp_kind/6._fp_kind -&
                         g*(1._fp_kind/7._fp_kind - g*(1._fp_kind/8._fp_kind))))))))
                 endif
                 fex0log2 = -fex0log*vg*vg
                 d1fex0log2 =&
                      -ug*(vg*vg + 2._fp_kind*fex0log2)
                 d2fex0log2 =&
                      -dug*(vg*vg + 2._fp_kind*fex0log2) -&
                      2._fp_kind*ug*(-vg*vg*ug + d1fex0log2)
                 d3fex0log2 =&
                      -dug2*(vg*vg + 2._fp_kind*fex0log2) -&
                      4._fp_kind*dug*(-vg*vg*ug + d1fex0log2) -&
                      2._fp_kind*ug*(2._fp_kind*vg*vg*ug*ug - vg*vg*dug + d2fex0log2)

                 fex0ggg = d3fex0log2*fex0 + 3._fp_kind*d2fex0log2*fex0g +&
                      3._fp_kind*d1fex0log2*fex0gg + fex0log2*fex0ggg
                 fex0fgg = d2fex0log2*fex0f + 2._fp_kind*d1fex0log2*fex0fg +&
                      fex0log2*fex0fgg
                 fex0ffg = d1fex0log2*fex0ff + fex0log2*fex0ffg
                 fex0fff = fex0log2*fex0fff
                 fex0gg = d2fex0log2*fex0 + 2._fp_kind*d1fex0log2*fex0g +&
                      fex0log2*fex0gg
                 fex0fg = d1fex0log2*fex0f + fex0log2*fex0fg
                 fex0ff = fex0log2*fex0ff
                 fex0g = d1fex0log2*fex0 + fex0log2*fex0g
                 fex0f = fex0log2*fex0f
                 fex0 = fex0log2*fex0
              endif
              ! reverse order to preserve values on RHS.
              ! ggg
              sum(10) = fex0ggg*sum(1) + 3._fp_kind*fex0gg*sum(3) +&
                   3._fp_kind*fex0g*sum(6) + fex0*sum(10)
              ! fgg
              sum(9) = fex0fgg*sum(1) + fex0gg*sum(2) +&
                   2._fp_kind*(fex0fg*sum(3) + fex0g*sum(5)) +&
                   fex0f*sum(6) + fex0*sum(9)
              ! ffg
              sum(8) = fex0ffg*sum(1) + fex0ff*sum(3) +&
                   2._fp_kind*(fex0fg*sum(2) + fex0f*sum(5)) +&
                   fex0g*sum(4) + fex0*sum(8)
              ! fff
              sum(7) = fex0fff*sum(1) + 3._fp_kind*fex0ff*sum(2) +&
                   3._fp_kind*fex0f*sum(4) + fex0*sum(7)
              ! gg
              sum(6) = fex0gg*sum(1) + 2._fp_kind*fex0g*sum(3) + fex0*sum(6)
              ! fg
              sum(5) = fex0fg*sum(1) + fex0g*sum(2) + fex0f*sum(3) +&
                   fex0*sum(5)
              ! ff
              sum(4) = fex0ff*sum(1) + 2._fp_kind*fex0f*sum(2) + fex0*sum(4)
              ! g
              sum(3) = fex0g*sum(1) + fex0*sum(3)
              ! f
              sum(2) = fex0f*sum(1) + fex0*sum(2)
              sum(1) = fex0*sum(1)
              ! Transform to tcl derivative from ln g
              ! derivative recalling that dlng/dlntc = 1, and dlng/dlnf = 0.5*uf.
              sum(4) = sum(4) +&
                   0.5_fp_kind*(duf*sum(3) + uf*sum(5)) +&
                   0.5_fp_kind*uf*(sum(5) + 0.5_fp_kind*uf*sum(6))
              sum(7) = sum(7) +&
                   0.5_fp_kind*(duf2*sum(3) + 2._fp_kind*duf*sum(5) + uf*sum(8)) +&
                   0.5_fp_kind*(duf*sum(5) + uf*sum(8) +&
                   0.5_fp_kind*uf*(2._fp_kind*duf*sum(6) + uf*sum(9))) +&
                   0.5_fp_kind*uf*(sum(8) +&
                   0.5_fp_kind*(duf*sum(6) + uf*sum(9)) +&
                   0.5_fp_kind*uf*(sum(9) + 0.5_fp_kind*uf*sum(10)))
              ! do sum(5) later than sum(7) so that rhs of sum(7) is undisturbed.
              sum(5) = sum(5) + 0.5_fp_kind*uf*sum(6)
              sum(8) = sum(8) +&
                   0.5_fp_kind*(duf*sum(6) + uf*sum(9)) +&
                   0.5_fp_kind*uf*(sum(9) + 0.5_fp_kind*uf*sum(10))
              sum(9) = sum(9) + 0.5_fp_kind*uf*sum(10)
              sum(2) = sum(2) + 0.5_fp_kind*uf*sum(3)
              fexl = fexl + sum(1)
              fexlf = fexlf + sum(2)
              fexlt = fexlt + sum(3)
              fexlf2 = fexlf2 + sum(4)
              fexlft = fexlft + sum(5)
              fexlt2 = fexlt2 + sum(6)
              fexl_psi = fexl_psi + sum(2)/dpsidf
              fexl_psif = fexl_psif +&
                   (sum(4) - sum(2)*dpsidf2/dpsidf)/dpsidf
              fexl_psit = fexl_psit + sum(5)/dpsidf
              fexl_psif2 = fexl_psif2 +&
                   (sum(7) - (2._fp_kind*sum(4)*dpsidf2 + sum(2)*dpsidf3 -&
                   2._fp_kind*sum(2)*dpsidf2*dpsidf2/dpsidf)/dpsidf)/dpsidf
              fexl_psift = fexl_psift +&
                   (sum(8) - sum(5)*dpsidf2/dpsidf)/dpsidf
              fexl_psit2 = fexl_psit2 + sum(9)/dpsidf
           endif
        enddo
        ! at this point K integral is stored in fexl.
        ! Now,transform fexl depending on iffexi and iffexk values.
        if(iffexi) then
           fexi = fexl*fexl
           fexif = 2._fp_kind*fexl*fexlf
           fexit = 2._fp_kind*fexl*fexlt
           fexif2 = 2._fp_kind*(fexlf*fexlf + fexl*fexlf2)
           fexift = 2._fp_kind*(fexlt*fexlf + fexl*fexlft)
           fexit2 = 2._fp_kind*(fexlt*fexlt + fexl*fexlt2)
           fexi_psi = 2._fp_kind*fexl*fexl_psi
           fexi_psif = 2._fp_kind*(fexlf*fexl_psi + fexl*fexl_psif)
           fexi_psit = 2._fp_kind*(fexlt*fexl_psi + fexl*fexl_psit)
           fexi_psif2 = 2._fp_kind*(fexlf2*fexl_psi + 2._fp_kind*fexlf*fexl_psif +&
                fexl*fexl_psif2)
           fexi_psift = 2._fp_kind*(fexlft*fexl_psi + fexlf*fexl_psit +&
                fexlt*fexl_psif + fexl*fexl_psift)
           fexi_psit2 = 2._fp_kind*(fexlt2*fexl_psi + 2._fp_kind*fexlt*fexl_psit +&
                fexl*fexl_psit2)
        endif
        if(iffexk) then
           fex0 = -con_nr1_ratio*beta*beta
           fexk = fex0*fexl
           fexkf = fex0*fexlf
           fexkt = fex0*(2._fp_kind*fexl + fexlt)
           fexkf2 = fex0*fexlf2
           fexkft = fex0*(2._fp_kind*fexlf + fexlft)
           fexkt2 = fex0*(4._fp_kind*fexl + 4._fp_kind*fexlt + fexlt2)
           fexk_psi = fex0*fexl_psi
           fexk_psif = fex0*fexl_psif
           fexk_psit = fex0*(2._fp_kind*fexl_psi + fexl_psit)
           fexk_psif2 = fex0*fexl_psif2
           fexk_psift = fex0*(2._fp_kind*fexl_psif + fexl_psift)
           fexk_psit2 =&
                fex0*(4._fp_kind*fexl_psi + 4._fp_kind*fexl_psit + fexl_psit2)
        endif
     endif
     ! from Paper IV equation 21:
     ! fex = fex0*(K^2 - J - con_nr1_ratio*beta^2 K)
     ! and fex0 = 4 pi(e m^2 c^2/h^2)^2 when converted to cgs and
     ! multiplied by -kT/V in accordance with fex definition.
     fex0 = con_d1_ratio*con_exchange
     fex = fex0*(fexi - fexj + fexk)
     fexf = fex0*(fexif - fexjf + fexkf)
     fext = fex0*(fexit - fexjt + fexkt)
     fexf2 = fex0*(fexif2 - fexjf2 + fexkf2)
     fexft = fex0*(fexift - fexjft + fexkft)
     fext2 = fex0*(fexit2 - fexjt2 + fexkt2)
     fex_psi = fex0*(fexi_psi - fexj_psi + fexk_psi)
     fex_psif = fex0*(fexi_psif - fexj_psif + fexk_psif)
     fex_psit = fex0*(fexi_psit - fexj_psit + fexk_psit)
     fex_psif2 = fex0*(fexi_psif2 - fexj_psif2 + fexk_psif2)
     fex_psift = fex0*(fexi_psift - fexj_psift + fexk_psift)
     fex_psit2 = fex0*(fexi_psit2 - fexj_psit2 + fexk_psit2)
  else
     error stop 'exchange_gcpf: invalid ifexchange'
  endif
end subroutine exchange_gcpf
