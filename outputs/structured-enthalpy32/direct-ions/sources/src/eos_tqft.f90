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

! eos_tqft should only be called after all auxiliary variable
! derivatives have been converged with NR procedure.

! This routine calls whole eos procedure again to determine fl and tl
! derivatives of all thermodynamic quantities using ifnr = 0 and the
! converged auxiliary variables and their derivatives as input.

! The above documentation of eos_tqft needs expansion concerning the
! arguments of eos_tqft, but as an interim measure use the eos_calc
! argument notes instead since many of the eos_tqft arguments are the
! same as the equivalently named eos_calc arguments.

! eos_calc calls ionize then fixes up the result consistent with
! formation of hydrogen molecules.
! the purpose of the two routines is to calculate ionization fractions
! and weighted sums over those quantities.  the weighted sums are
! ultimately used in an iterative way to calculate chemical potentials
! according to some free energy model.  The appropriate difference in
! chemical potentials are then combined to form the dv quantities,
! the change in the equilibrium constant of an ion relative to the
! un-ionized reference state and the free electron
! dv = (partial F/partial nref - partial F/partial nion - ion*partial F/partial ne)/kT.
! These dv values (held in an array) are required input to eos_calc.
!
! The free energy model:
! it is the responsibility of the calling programme to
! fill in the dv quantities for the various ions in a manner
! consistent with a particular free energy model.  It should be noted that
! at the startup of the iteration, the free-energy terms corresponding
! to pressure ionization and the Coulomb effect are approximated as
! functions of only ne, so the dv array elements are quite similar
! to each other.  (The electron exchange term is always just a function
! of ne.) Later as more complicated models are used for the Coulomb
! effect and pressure ionization, the dv array elements can become
! more varied.
!
! input quantitites:
! ifexcited > 0 means use excited states (must have Planck-Larkin or ifpi = 3 or 4).
!    0 < ifexcited < 10 means use approximation to explicit summation
!   10 < ifexcited < 20 means use explicit summation
!   mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules) and helium.
!   mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
!   mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
! ifnuform: usually set to false so that h2, h2plus, h_ion, he_ion, he_ion2, and
!   extrasum returned in standard number density form.
!   however, on last call to eos_calc for a given fl, tl, this
!   quantity should be set to true so that these quantities
!   returned in nu = n/(rho*avogadro) form.  This reduces
!   the significance loss in the entropy due to pressure ionization.
! ifsame_under = .true. means use same underflow zeroing for each
!   ion as in previous call.  This option is useful for removing
!   small discontinuities caused by variations in the underflow
!   zeroing which sometimes foil the last stages of convergence.
! ifsame_under = .false. means calculate underflow zeroing.
! ifnr = 0
!   calculate derivatives of nuvar (delivered through common
!   block), hne, sum0, sum2, s, extrasum, h_ion, sh, he_ion, he_ion2, h_neutral,
!   hequil, and sumpl1 wrt f, t and sumpl0 wrt f with no other
!   variables fixed, i.e., use input dvf and dvt and chain rule.
!   N.B. chain of calling routines depend on the assertion that
!   ifnr = 0 means no *_dv variables should be read or written.
! ifnr = 1
!   calculate derivatives of nuvar (delivered through common
!   block), hne, sum0, sum2, extrasum, h_ion, he_ion, and he_ion2 wrt dv with
!   f, t fixed.  n.b. list of variables is subset of those listed
!   for ifnr = 0 because we only need dv derivatives for variables
!   used to calculate (output) auxiliary variables.
! ifnr = 2 (should not occur)
! ifnr = 3 same as combination of ifnr = 0 and ifnr = 1.
!   Note this is quite different in detail from ifnr = 3
!   interpretation for many other routines, but the general
!   motivation is the same for all ifnr = 3 results; calculate
!   both f and t derivatives and other derivatives.
! inv_ion(nions+2) maps ion index to contiguous ion index used in NR iteration dot products.
! max_index is the maximum contiguous index used in NR iteration
! partial_elements(n_partial_elements+2) index of elements treated as
!   partially ionized consistent with ifelement.
! ion_end(nelements) keeps track of largest ion index for each element.
! ifionized = 0 means all elements treated as partially ionized
! ifionized = 1 means trace metals fully ionized.
! ifionized = 2 means all elements treated as fully ionized.
! if_pteh controls whether pteh approximation used for Coulomb sums (1)
!   or whether use detailed Coulomb sum (0).
! if_mc controls whether special approximation used for metal part
!   of Coulomb sum for the partially ionized metals.
!   if_mc = 1, use approximation
!   if_mc = 0, don't use approximation
! ifreducedmass = 1 (use reduced mass in equilibrium constant)
! ifreducedmass = 0 (use electron mass in equilibrium constant)
! ifsame_abundances = 1, if call made with same abundances as previous
!   it is the responsibility of the calling programme to set this
!   flag where appropriate.  ifsame_abundances = 0 means slower
!   execution of ionize, but is completely safe against abundance
!   changes.
! ifmtrace just to keep track of changes in ifelement.
!   n.b. if ifionized is different from previous call or
!   ifsame_abundances = 0, or ifmtrace different from previous
!   call, then local variables ionized_hne, ionized_sum0,
!   ionized_sum2, and ionized_extrasum(2) are recalculated.
! iatomic_number(nelements), the atomic number of each element.
!   n.b. the code logic now demands there are iatomic_number
!   ionized states + 1 neutral state for each element.  Thus, the
!   ionization potentials must include all ions, not just first
!   ions (as for trace metals in old code).
! ifpi = 0, 1, 2, 3, 4 if no, pteh, geff, mhd, or Saumon-style pressure ionization
! ifpl = 0, 1 if there is no or if there is planck-larkin
! ifmodified > 0 or not affects ifpi_fit inside excitation_sum.
! ifh2 non-zero implies h2 is included.
! ifh2plus non-zero implies h2plus is included.
! izlo is lowest core charge (1 for H, He, 2 for He+, etc.) for all species.
! izhi is highest core charge for all species.
! bmin(izhi) is the minimum bion for all species with a given core charge.
! nmin(izhi) is the minimum excited principal quantum number for
!   all species with a given core charge.
! nmin_max(izhi) is the largest minimum excited principal
!   quantum number for all species with a given core charge.
! nmin_species(nions+2) is the minimum excited principal quantum number
!   organized by species.
! nmax is maximum principal quantum number included in sum
!   (to be compatible with opal which used nmax = 4 rather than infinity).
!   if(nmax > 300000) then treated as infinity in qryd_approx.
!   otherwise nmax is meant to be used with qryd_calc only (i.e.,
!   case for mhd approximations not programmed).
! eps(nelements) is an array of neps = nelements = 24 values
!   of relative abundance by weight divided by the appropriate
!   atomic weight. We use the atomic weight scale where
!   (un-ionized) C(12) has a weight of 12.00000000....  All
!   weights are for the un-ionized element. The eps value for an
!   element should be the sum of the individual isotopic eps
!   values for that element.  The eps array refers to the elements
!   in the following order:
!     H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
! tl = log(t)
! tc2 = c2/t
! bi(nions+2) ionization potentials (cm**-1) in ion order:
!   h_ion, he_ion2, he_ion2+, c, c+, etc.  H2 and H2+ are last two.
!   n.b. the code logic now demands there are iatomic_number ionized
!     states + 1 neutral state for each element.  Thus, the
!     ionization potentials must include all ions, not just
!     first ions (as for trace metals in old code).
! h2diss = h2 dissociation energy
! plop(nions+2), plopt, plopt2, planck-larkin occupation probability
!   and derivatives in ion order (starting with neutral in each
!   case and avoiding bare nuclei as in bi.
! r_ion3(nions+2), *the cube of the* effective radii (in ion
!   order but going from neutral to next to bare ion for each
!   species) of MDH interaction between all but bare nucleii
!   species and ionized species.  last 2 are H2 and H2+.
! r_neutral(nelements+2) effective radii for MDH neutral-neutral
!    interactions.  Last two (used externally) are H2 and H2+ (the
!    only ionic species in the MHD model with a non-zero
!    hard-sphere radius).
! ifelement(nelements) = 1 if element treated as partially
!   ionized (depending on abundance, and trace metal treatment), 0
!   if element treated as fully ionized or has no abundance.
! dvzero(nelements) is a zero point shift (zero for hydrogen)
!   that must be added to dv in all uses below.  Sometimes (high
!   density, high ionization) this quantity can be larger than
!   1.d9 so if we avoid using it in differences of dv quanties we
!   can gain many significant digits.
! dv(nions) change in equilibrium constant (explained above)
! dvf(nions) = d dv/d ln f
! dvt(nions) = d dv/d ln t
! nion(nions), charge on ion in ion order (must be same order
!   as bi) e.g., for H+, He+, He++, etc.
! re, ref, ret: fermi-dirac integral and derivatives
! output quantities:
! ne, nef, net is nu(e) = n(positive ions)/(Navogadro*rho) and its
!   derivatives calculated from all elements.
! sion, sionf, siont is the ideal entropy/R per unit mass and its
!   derivatives with respect to lnf and ln t.
!   n.b. if ionized, then the contribution from all elements is
!   calculated in ionize.f.  If partially ionized, then ionize separates
!   the hydrogen from the rest of the components.  If molecular, the
!   hydrogen component from ionize is completely ignored, and
!   recalculated in this routine.
!   n.b.  sion returns a component of entropy/R per unit mass,
!   sion = the sum over all non-electron species of
!     -nu_i*[-5/2 - 3/2 ln T + ln(alpha) - ln Na
!     + ln (n_i/[A_i^{3/2} Q_i]) - dln Q_i/d ln T],
!     where nu_i = n_i/(Na rho), Na is the Avogadro number,
!     A_i is the atomic weight,
!     alpha =  (2 pi k/[Na h^2])^(-3/2) Na
!     (see documentation of alpha^2 in the constants module), and
!     Q_i is the internal ideal partition function
!     of the non-Rydberg states.  (Currently, we calculate this
!     partition function by the statistical weights of the ground states
!     of helium and the combined lower states (roughly
!     approximated) of each of the metals.
!     helium is subsequently corrected for detailed excitation
!     of the non-Rydberg states, and this crude approximation for the metals
!     is not currently corrected.  Thus, in all *current* monatomic
!     cases ln Q_i is a constant, and dln Q_i/d ln T is zero, but this
!     will change for the metals eventually.
!     currently, hydrogen is treated exactly for all cases (full ionization,
!     partial ionization, partial molecular formation).
!     From the equilibrium constant approach and the monatomic species
!     treated in this subroutine (molecules treated outside) we have
!     -nu_i ln (n_i/[A_i^{3/2} Q_i]) =
!     -nu_i * [ln (n_neutral/[A_neutral^{3/2} Q_neutral) +
!     (- chi_i/kT + dv_i)]
!     if we sum this term over all species
!     of an element without molecules we obtain
!     s_element = - eps * ln (eps*n_neutral/sum(n))
!       - eps ln (alpha/(A_neutral^{3/2} Q_neutral)) -
!       - sum over all species of the element of nu_i*(-chi/kT + dv(i))
!     where we have ignored the term
!     eps * (-5/2 - 3/2 ln T + ln rho)
!     (taking into account the first 4 terms above).
!     If hydrogen molecules are included,
!     the result is the same except for the addition of the
!     d ln Q_i/d ln T term (which will also appear for the metals eventually)
!     and the 5/2 factor is multiplied by sum over all species which is
!     corrected to eps by subtracting 5/2 (nu(H2) + nu(H2+)) from sion below.
!     note a final correction of the s zero point occurs
!     in free_eos_detailed.f which puts back the ignored terms for both
!     the case of molecules and no molecules.
! uion is the ideal internal energy (cm^-1 per unit mass divided by avogadro).
! h_ion, h_ionf, h_iont, h_ion_dv is n(H+) and its derivatives with respect to ln f, ln t, and dv.
! he_ion, he_ionf, he_iont, he_ion_dv is n(He+) and its derivatives with respect to ln f, ln t, and dv.
! he_ion2, he_ion2f, he_ion2t, he_ion2_dv is n(He++) and its derivatives with respect to ln f, ln t, and dv.
! The following are returned only if .not.( ifpteh.eq.1.or.if_mc.eq.1)
! sum0 and f,t,dv derivatives, sum over positive charge number densities with uniform weights.
! sum2 and f,t,dv derivatives, sum over positive charge number densities weighted by charge^2.
! The following are returned only if ifpi = 3 or 4.
! extrasum(nextrasum) weighted sums over n(i).
!   for iextrasum = 1,nextrasum-2, sum is only over
!   neutral species + H2+ (the only species with
!   non-zero radii according to the MHD model) and
!   weight is r_neutral^{iextrasum-1}.
!   for iextrasum = nextrasum-1, sum is over all ionized species including
!   bare nucleii, but excluding free electrons, weight is Z^1.5.
!   for iextrasum = nextrasum, sum is over all species excluding bare nucleii
!   and free electrons, the weight is rion^3.
! extrasumf(nextrasum) = partial of extrasum/partial ln f
! extrasumt(nextrasum) = partial of extrasum/partial ln t
! The following are returned only if ifpl = 1
! sumpl0, sumpl1 and sumpl2 = weighted sums over non-H nu(i) = n(i)/(rho/H).
!   for sumpl0 sum is over Planck-Larkin occupation probabilities.
!   for sumpl1 sum is over Planck-Larkin occupation probabilities + d ln w/d ln T
!   for sumpl2, sum is over Planck-Larkin d ln w/d ln T
! sumpl0f, sumpl0_dv = derivative of sumpl0 wrt ln f and dv.
! sumpl1f and sumpl1t = derivatives of sumpl1 wrt lnf and lnt.
! rl, rf, rt are the ln mass density and derivatives
! h2, h2f, h2t, h2_dv are n(H2) and derivatives.
! h2plus, h2plusf, h2plust are n(H2+) and derivatives.
! h2plus_dv is the h2plus derivatives wrt dv.  n.b. this vector only
!   returned if molecular hydrogen is calculated.  The calling
!   routine uses it only if if_mc.eq.1.and.ifh2plus.gt.0.
! The following are returned only if ifexcited.gt.0 and ifpi.eq.3.or4.
! xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
!   partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
! xextrasumf(4), xextrasumt(4), xextrasum_dv(nions+2,4) are the
!  fl and tl derivatives (ifnr.eq.0) or dv derivatives (ifnr.eq.1)
!  required so that can store nuh2 and nuh2plus for eos_free_calc.

!> This eos_tqft subroutine calculates for a converged EOS solution
!> the ifnr = 0 partial derivatives of auxiliary variables wrt fl, tl,
!> and dv.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine eos_tqft(&
     verbosity, dv_aux,&
     lambda, gamma_e,&
     nux, nuy, nuz,&
     n_partial_ions, n_partial_aux,&
     sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
     partial_ions, f, eta, wf, t, n_e, pstar,&
     dve_exchange, dve_exchangef, dve_exchanget,&
     full_sum0, full_sum1, full_sum2, charge, charge2,&
     dv_pl, dv_plt,&
     ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
     ifsame_zero_abundances, ifexcited,&
     inv_aux,&
     inv_ion, max_index,&
     partial_elements, ion_end, mion_end,&
     ifionized, if_pteh, if_mc, ifreducedmass,&
     ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,&
     ifpl, ifmodified, ifh2, ifh2plus,&
     izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
     eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
     r_ion3, r_neutral,&
     ifelement, dvzero, dv, dvf, dvt, nion,&
     rhostar,&
     ne, nef, net, sion, sionf, siont, uion,&
     h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t,&
     h_ion_dv, he_ion_dv, he_ion2_dv,&
     sum0, sum0f, sum0t, sum0_dv,&
     sum2, sum2f, sum2t, sum2_dv,&
     extrasum, extrasumf, extrasumt, extrasum_dv,&
     sumpl1, sumpl1f, sumpl1t, sumpl2,&
     rl, rf, rt, r_dv,&
     h2, h2f, h2t, h2plus, h2plusf, h2plust, h2plus_dv,&
     xextrasum, xextrasumf, xextrasumt, xextrasum_dv, info)

  use mod_free_eos_constants, only: avogadro, cpe
  use mod_coulomb, only: master_coulomb
  use mod_eos_calc, only: eos_calc
  use mod_pi, only: mdh_pi, fjs_pi, pteh_pi
  use mod_excitation, only: excitation_pi
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  logical, intent(in) :: ifdvzero

  integer, intent(in) :: nion(:)
  integer, intent(in) :: ion_end(:)
  integer, intent(in) :: verbosity, ifreducedmass, ifmodified, ifexcited, nmax
  integer, intent(in) :: nmin(:), nmin_max(:), izlo, izhi
  integer, intent(in) :: nmin_species(:)
  integer, intent(in) :: iatomic_number(:)
  integer, intent(in) ::&
       ifh2, ifh2plus, ifmtrace,&
       if_pteh, ifcoulomb_mod, if_mc, if_dc,&
       ifionized, max_index,&
       ifpi_local, ifpl,&
       ifsame_abundances,&
       ifsame_zero_abundances
  integer, intent(in) :: ifdv(:), ifelement(:),&
       partial_ions(:), inv_ion(:),&
       partial_elements(:)
  integer, intent(in) :: inv_aux(:)
  integer, intent(in) :: n_partial_aux
  integer, intent(in) :: mion_end, n_partial_ions

  integer, intent(out) :: info

  real(fp_kind), intent(in) :: bi(:), h2diss, t, tl, eps(:), eta
  real(fp_kind), intent(in) :: rhostar(:), pstar(:)
  real(fp_kind), intent(in) ::&
       dve_exchange, dve_exchangef, dve_exchanget,&
       plop(:), plopt(:), plopt2(:),&
       dv_pl(:), dv_plt(:)
  real(fp_kind), intent(in) ::&
       full_sum0, full_sum1, full_sum2,&
       r_ion3(:), r_neutral(:),&
       charge(:), charge2(:),&
       tc2
  real(fp_kind), intent(in) :: f, wf, n_e
  real(fp_kind), intent(in) :: bmin(:)
  real(fp_kind), intent(in) ::&
       sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
       nux, nuy, nuz

  real(fp_kind), intent(out) ::&
       dvzero(:),&
       dv(:), dvf(:), dvt(:),&
       lambda, gamma_e,&
       ne, nef, net,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       sion, sionf, siont, uion,&
       h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t,&
       h2, h2f, h2t,&
       h_ion_dv(:), he_ion_dv(:), he_ion2_dv(:),&
       r_dv(:), sum0_dv(:), sum2_dv(:),&
       h2plus, h2plusf, h2plust, h2plus_dv(:),&
       extrasum_dv(:,:),&
       xextrasum_dv(:,:)

  real(fp_kind), intent(inout) ::&
       sum0, sum0f, sum0t,&
       sum2, sum2f, sum2t,&
       extrasum(:), extrasumf(:), extrasumt(:),&
       rl, rf, rt,&
       xextrasum(:), xextrasumf(:), xextrasumt(:),&
       dv_aux(:,:)

  ! Internal variables

  logical ifsame_under, ifnuform

  integer nionsp2, naux, nelements, n_partial_elements, nextrasum
  integer&
       ion0,&
       index,&
       ion, index_ion,&
       ifnr,&
       ielement

  real(fp_kind)&
       dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,&
       sum0ne, sum2ne,&
       sumion0, sumion2,&
       dve_coulomb, dve_coulombf, dve_coulombt,&
       dve, dvef, dvet, dve0, dve2,&
       dve_pi, dve_pif, dve_pit,&
       sumpl0, sumpl0f,&
       fion, fionf,&
       rho,&
       dvh, dvhf, dvht,&
       dvhe1, dvhe1f, dvhe1t,&
       dvhe2, dvhe2f, dvhe2t,&
       nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt

  real(fp_kind), allocatable ::&
       dve_pi_aux(:),&
       sumpl0_dv(:),&
       fion_dv(:),&
       h2_dv(:),&
       dvh_aux(:),&
       dvhe1_aux(:),&
       dvhe2_aux(:)

  integer, parameter :: maxfjs_aux = 4
  integer, parameter :: nstar = 9
  ! maximum core charge for non-bare ion
  integer, parameter :: maxcorecharge = 28
  integer, parameter :: nxextrasum = 4
  integer, parameter :: nions_localp2 = 318
  ! iextraoff is the number of non-extrasum auxiliary variables
  integer, parameter :: iextraoff = 7

  nionsp2 = size(nion)
  naux = size(inv_aux)
  nelements = size(ion_end) - 2
  n_partial_elements = size(partial_elements) - 2
  nextrasum = size(extrasum)

  ! Sanity checks
  if(nionsp2.ne.nions_localp2) error stop 'eos_tqft: bad nionsp2 value'
  if(&
       nionsp2.ne.size(nmin_species).or.&
       nionsp2.ne.size(ifdv).or.&
       nionsp2.ne.size(partial_ions).or.&
       nionsp2.ne.size(inv_ion).or.&
       nionsp2.ne.size(bi).or.&
       nionsp2.ne.size(plop).or.&
       nionsp2.ne.size(plopt).or.&
       nionsp2.ne.size(plopt2).or.&
       nionsp2.ne.size(dv_pl).or.&
       nionsp2.ne.size(dv_plt).or.&
       nionsp2.ne.size(r_ion3).or.&
       nionsp2.ne.size(charge).or.&
       nionsp2.ne.size(charge2).or.&
       nionsp2.ne.size(dv).or.&
       nionsp2.ne.size(dvf).or.&
       nionsp2.ne.size(dvt).or.&
       nionsp2.ne.size(h_ion_dv).or.&
       nionsp2.ne.size(he_ion_dv).or.&
       nionsp2.ne.size(he_ion2_dv).or.&
       nionsp2.ne.size(r_dv).or.&
       nionsp2.ne.size(sum0_dv).or.&
       nionsp2.ne.size(sum2_dv).or.&
       nionsp2.ne.size(h2plus_dv).or.&
       nionsp2.ne.size(xextrasum_dv,1).or.&
       nionsp2.ne.size(dv_aux,1))&
       error stop 'eos_tqft: inconsistent nionsp2 array sizes'
  if(&
       naux.ne.size(dv_aux,2))&
       error stop 'eos_tqft: inconsistent naux array sizes'
  if(&
       nelements.ne.size(iatomic_number).or.&
       nelements.ne.size(ifelement).or.&
       nelements.ne.size(eps).or.&
       nelements+2.ne.size(r_neutral).or.&
       nelements.ne.size(dvzero))&
       error stop 'eos_tqft: inconsistent nelements array sizes'
  if(&
       nextrasum.ne.size(extrasumf).or.&
       nextrasum.ne.size(extrasumt).or.&
       nextrasum.ne.size(extrasum_dv,2))&
       error stop 'eos_tqft: inconsistent nextrasum array sizes'
  if(&
       nstar.ne.size(rhostar).or.&
       nstar.ne.size(pstar))&
       error stop 'eos_tqft: inconsistent nstar array sizes'

  if(&
       maxcorecharge.ne.size(nmin).or.&
       maxcorecharge.ne.size(nmin_max).or.&
       maxcorecharge.ne.size(bmin))&
       error stop 'eos_tqft: inconsistent maxcorecharge array sizes'
  if(&
       nxextrasum.ne.size(xextrasum_dv,2).or.&
       nxextrasum.ne.size(xextrasum).or.&
       nxextrasum.ne.size(xextrasumf).or.&
       nxextrasum.ne.size(xextrasumt))&
       error stop 'eos_tqft: inconsistent nxextrasum array sizes'

  allocate(&
       dve_pi_aux(maxfjs_aux),&
       sumpl0_dv(nionsp2),&
       fion_dv(nionsp2),&
       h2_dv(nionsp2),&
       dvh_aux(maxfjs_aux),&
       dvhe1_aux(maxfjs_aux),&
       dvhe2_aux(maxfjs_aux)&
       )

  if(if_taint_allocated_real) then
     call taint_allocated_real(dve_pi_aux)
     call taint_allocated_real(sumpl0_dv)
     call taint_allocated_real(fion_dv)
     call taint_allocated_real(h2_dv)
     call taint_allocated_real(dvh_aux)
     call taint_allocated_real(dvhe1_aux)
     call taint_allocated_real(dvhe2_aux)
    endif

  ifnr = 0
  rho = exp(rl)
  do index = 1, n_partial_elements
     ielement = partial_elements(index)
     dvzero(ielement) = 0._fp_kind
  enddo
  do index_ion = 1,max_index
     ion = partial_ions(index_ion)
     dv(ion) = 0._fp_kind
     dvf(ion) = 0._fp_kind
     dvt(ion) = 0._fp_kind
  enddo
  if(ifpi_local.eq.0) then
     dve_pi = 0._fp_kind
     dve_pif = 0._fp_kind
     dve_pit = 0._fp_kind
  elseif(ifpi_local.eq.1) then
     call pteh_pi(ifmodified,&
          rhostar(1), rhostar(2), rhostar(3), t,&
          dve_pi, dve_pif, dve_pit)
  elseif(ifpi_local.eq.2) then
     ! use FJS pressure ionization
     nx = nux*rho*avogadro
     ny = nuy*rho*avogadro
     nz = nuz*rho*avogadro
     nxf = nx*rf
     nxt = nx*rt
     nyf = ny*rf
     nyt = ny*rt
     nzf = nz*rf
     nzt = nz*rt
     call fjs_pi(ifnr, ifmodified,&
          t, nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,&
          h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t,&
          n_e, n_e*rhostar(2), n_e*rhostar(3),&
          dvh, dvhf, dvht, dvh_aux,&
          dvhe1, dvhe1f, dvhe1t, dvhe1_aux,&
          dvhe2, dvhe2f, dvhe2t, dvhe2_aux,&
          dve_pi, dve_pif, dve_pit, dve_pi_aux)
     ! update dv keeping in mind that calculated dvh etc. have
     ! dve_pi included.  the dve_pi quantity
     ! must be subtracted so that further free_eos_detailed logic
     ! which adds dve_pi works properly.
     if(eps(1).gt.0._fp_kind) then
        dv(1) = dv(1) + dvh - dve_pi
        dvf(1) = dvf(1) + dvhf - dve_pif
        dvt(1) = dvt(1) + dvht - dve_pit
     endif
     if(eps(2).gt.0._fp_kind) then
        if(ifdvzero) then
           dv(2) = dv(2) + dvhe1 - dve_pi
           dv(3) = dv(3) + dvhe1 + dvhe2 - 2._fp_kind*dve_pi
        else
           dvzero(2) = dvzero(2) + dvhe1 + dvhe2 - 2._fp_kind*dve_pi
           dv(2) = dv(2) - dvhe2 + dve_pi
        endif
        dvf(2) = dvf(2) + dvhe1f - dve_pif
        dvt(2) = dvt(2) + dvhe1t - dve_pit
        dvf(3) = dvf(3) + dvhe1f + dvhe2f - 2._fp_kind*dve_pif
        dvt(3) = dvt(3) + dvhe1t + dvhe2t - 2._fp_kind*dve_pit
     endif
     ! no change in H2 equilibrium constant because fjs_pi is independent
     ! n(H) or n(H2).
  elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
     ! use mdh-like pressure ionization and dissociation
     call mdh_pi(&
          ifdvzero, ifnr, inv_aux, iextraoff, inv_ion,&
          partial_elements,&
          ion_end(:mion_end),&
          ifmodified,&
          r_ion3, nion, r_neutral,&
          extrasum, extrasumf, extrasumt,&
          ifdv, dvzero, dv, dvf, dvt, dv_aux)
     ! MDH pressure ionization has no N_e dependence
     dve_pi = 0._fp_kind
     dve_pif = 0._fp_kind
     dve_pit = 0._fp_kind
  endif
  ! calculate change in equilibrium constant due to electron degeneracy
  dve = dve_pi - eta
  dvef = dve_pif - wf
  dvet = dve_pit
  ! n.b. for master_coulomb to work properly, sum0 and sum2
  ! are considered to be functions of ne, f, t.
  ! for pteh sum approximation, the sums depend on ne, and
  ! there is no explicit f, t dependence.
  if(if_pteh.eq.1) then
     ! PTEH (full ionization) approximation to sum0, sum2
     sum0ne = full_sum0/full_sum1
     sum0 = sum0ne*n_e
     sum2ne = full_sum2/full_sum1
     sum2 = sum2ne*n_e
     sum0f = 0._fp_kind
     sum0t = 0._fp_kind
     sum2f = 0._fp_kind
     sum2t = 0._fp_kind
  else
     ! derivative wrt ne is zero (this must be
     ! asserted otherwise master_coulomb doesn't work
     ! correctly)
     sum0ne = 0._fp_kind
     sum2ne = 0._fp_kind
  endif
  if(if_mc.eq.1) then
     ! in just this case
     ! ionize and eos_calc have only calculated fully
     ! ionized part/(rho*avogadro).
     sumion0 = sum0_mc*rho*avogadro
     sumion2 = sum2_mc*rho*avogadro
     sum0 = sumion0 + h2plus + h_ion*(1._fp_kind+hcon_mc) +&
          he_ion +  he_ion2
     sum2 = sumion2 + h2plus + h_ion*(1._fp_kind+hcon_mc) +&
          he_ion + he_ion2*(4._fp_kind+hecon_mc)
     sum0f = sumion0*rf + h2plusf + h_ionf*(1._fp_kind+hcon_mc) +&
          he_ionf + he_ion2f
     sum0t = sumion0*rt + h2plust + h_iont*(1._fp_kind+hcon_mc) +&
          he_iont + he_ion2t
     sum2f = sumion2*rf + h2plusf + h_ionf*(1._fp_kind+hcon_mc) +&
          he_ionf + he_ion2f*(4._fp_kind+hecon_mc)
     sum2t = sumion2*rt + h2plust + h_iont*(1._fp_kind+hcon_mc) +&
          he_iont + he_ion2t*(4._fp_kind+hecon_mc)
  endif
  call master_coulomb(ifnr,&
       rhostar, f,&
       sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,&
       n_e, t, cpe*pstar(1), pstar, lambda, gamma_e,&
       ifcoulomb_mod, if_dc, if_pteh,&
       dve_coulomb, dve_coulombf, dve_coulombt, dve0, dve2,&
       dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
  dve = dve + dve_coulomb
  dvef = dvef + dve_coulombf
  dvet = dvet + dve_coulombt
  ! add in pre-calculated exchange effects
  dve = dve + dve_exchange
  dvef = dvef + dve_exchangef
  dvet = dvet + dve_exchanget
  if(if_mc.ne.1) then
     do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if(ifdvzero.or.ielement.eq.1) then
        else
           dvzero(ielement) = dvzero(ielement) +&
                dv0 + charge(ion_end(ielement))*&
                (dve + charge(ion_end(ielement))*dv2)
        endif
        if(ielement.gt.1) then
           ion0 = ion_end(ielement-1)
        else
           ion0 = 0
        endif
        if(ifdvzero.or.ielement.eq.1) then
           dv(ion0+1:ion_end(ielement)) = dv(ion0+1:ion_end(ielement)) +&
                dv0 + charge(ion0+1:ion_end(ielement))*(dve + charge(ion0+1:ion_end(ielement))*dv2)
        else
           dv(ion0+1:ion_end(ielement)) = dv(ion0+1:ion_end(ielement)) +&
                real((nion(ion0+1:ion_end(ielement))-nion(ion_end(ielement))),fp_kind)*&
                (dve + real((nion(ion0+1:ion_end(ielement))+nion(ion_end(ielement))),fp_kind)*dv2)
        endif
        dvf(ion0+1:ion_end(ielement)) = dvf(ion0+1:ion_end(ielement)) +&
             dv0f + charge(ion0+1:ion_end(ielement))*(dvef + charge(ion0+1:ion_end(ielement))*dv2f)
        dvt(ion0+1:ion_end(ielement)) = dvt(ion0+1:ion_end(ielement)) +&
             dv0t + charge(ion0+1:ion_end(ielement))*(dvet + charge(ion0+1:ion_end(ielement))*dv2t)
     enddo
  else
     ! if_mc *is* 1
     do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if(ifdvzero.or.ielement.eq.1) then
        elseif(ielement.eq.2) then
           dvzero(ielement) = dvzero(ielement) +&
                dv0 + charge(ion_end(ielement))*&
                (dve + charge(ion_end(ielement))*dv2)
        else
           dvzero(ielement) = dvzero(ielement) +&
                charge(ion_end(ielement))*dve
        endif
        if(ielement.gt.1) then
           ion0 = ion_end(ielement-1)
        else
           ion0 = 0
        endif
        do ion = ion0+1,ion_end(ielement)
           if(ion.ge.4) then
              if(ifdvzero) then
                 dv(ion) = dv(ion) + charge(ion)*dve
              else
                 dv(ion) = dv(ion) +&
                      real((nion(ion)-nion(ion_end(ielement))),fp_kind)*dve
              endif
              dvf(ion) = dvf(ion) + charge(ion)*dvef
              dvt(ion) = dvt(ion) + charge(ion)*dvet
           elseif(ion.eq.1) then
              dv(ion) = dv(ion) +&
                   charge(ion)*dve + (1._fp_kind + hcon_mc)*&
                   (dv0 + charge2(ion)*dv2)
              dvf(ion) = dvf(ion) +&
                   charge(ion)*dvef + (1._fp_kind + hcon_mc)*&
                   (dv0f + charge2(ion)*dv2f)
              dvt(ion) = dvt(ion) +&
                   charge(ion)*dvet + (1._fp_kind + hcon_mc)*&
                   (dv0t + charge2(ion)*dv2t)
           elseif(ion.eq.2) then
              if(ifdvzero) then
                 dv(ion) = dv(ion) +&
                      dv0 + charge(ion)*(dve + charge(ion)*dv2)
              else
                 dv(ion) = dv(ion) +&
                      real((nion(ion)-nion(ion_end(ielement))),fp_kind)*&
                      (dve + real((nion(ion)+nion(ion_end(ielement))),fp_kind)*dv2)
              endif
              dvf(ion) = dvf(ion) +&
                   dv0f + charge(ion)*(dvef + charge(ion)*dv2f)
              dvt(ion) = dvt(ion) +&
                   dv0t + charge(ion)*(dvet + charge(ion)*dv2t)
           elseif(ion.eq.3) then
              if(ifdvzero) then
                 dv(ion) = dv(ion) + charge(ion)*dve + dv0 +&
                      (hecon_mc + charge2(ion))*dv2
              else
                 dv(ion) = dv(ion) + hecon_mc*dv2
              endif
              dvf(ion) = dvf(ion) + charge(ion)*dvef + dv0f +&
                   (hecon_mc + charge2(ion))*dv2f
              dvt(ion) = dvt(ion) + charge(ion)*dvet + dv0t +&
                   (hecon_mc + charge2(ion))*dv2t
           endif
        enddo
     enddo
  endif
  ! H2+
  if(ifdv(nionsp2).eq.1) then
     dv(nionsp2) = dv(nionsp2) + dv0 + dve + dv2
     dvf(nionsp2) = dvf(nionsp2) + dv0f + dvef + dv2f
     dvt(nionsp2) = dvt(nionsp2) + dv0t + dvet + dv2t
  endif
  ! add in Planck-Larkin occupation probability effect
  if(ifpl.eq.1) then
     do index_ion = 1,max_index
        ion = partial_ions(index_ion)
        dv(ion) = dv(ion) + dv_pl(ion)
        dvt(ion) = dvt(ion) + dv_plt(ion)
     enddo
  endif
  if(ifexcited.gt.0)&
       call excitation_pi(verbosity, ifexcited, ifsame_zero_abundances,&
          ifpl, ifpi_local, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,&
          partial_elements, ion_end,&
          tl, izlo, bmin(:izhi), nmin(:izhi), nmin_max(:izhi), nmin_species, nmax,&
          bi, plop, plopt, plopt2,&
          r_ion3, nion, r_neutral,&
          inv_aux, iextraoff,&
          extrasum, extrasumf, extrasumt,&
          xextrasum, xextrasumf, xextrasumt,&
          dv, dvf, dvt, dv_aux)
  ! last call to eos_calc should calculate nu = n/(rho*avogadro)
  ! form of h_ion, he_ion, he_ion2, and extrasum.
  ifnuform = .true.
  ! use same underflow limits as converged solution
  ifsame_under = .true.
  call eos_calc(verbosity, ifexcited, ifsame_zero_abundances,&
       ifnuform, ifsame_under, ifnr, inv_ion, max_index,&
       partial_elements, ion_end,&
       ifionized, if_pteh, if_mc, ifreducedmass,&
       ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,&
       ifpl, ifmodified, ifh2, ifh2plus,&
       izlo, bmin(:izhi), nmin(:izhi), nmin_max(:izhi), nmin_species, nmax,&
       eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
       r_ion3, r_neutral,&
       ifelement, dvzero, dv, dvf, dvt, nion,&
       rhostar(1), rhostar(2), rhostar(3),&
       ne, nef, net,&
       sion, sionf, siont,&
       uion,&
       h_ion, h_ionf, h_iont, h_ion_dv,&
       he_ion, he_ionf, he_iont, he_ion_dv,&
       he_ion2, he_ion2f, he_ion2t, he_ion2_dv,&
       fion, fionf, fion_dv,&
       sum0, sum0f, sum0t, sum0_dv,&
       sum2, sum2f, sum2t, sum2_dv,&
       extrasum, extrasumf, extrasumt, extrasum_dv,&
       sumpl0, sumpl0f, sumpl0_dv,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       rl, rf, rt, r_dv,&
       h2, h2f, h2t, h2_dv, h2plus, h2plusf, h2plust, h2plus_dv,&
       xextrasum, xextrasumf, xextrasumt, xextrasum_dv, info)
  if(info.ne.0) return
end subroutine eos_tqft
