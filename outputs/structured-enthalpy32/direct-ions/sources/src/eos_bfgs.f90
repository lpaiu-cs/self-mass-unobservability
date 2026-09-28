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

! eos_bfgs is obviously different from eos_calc, but here are the
! eos_calc argument notes as an interim documentation measure that
! assumes that many eos_bfgs arguments are the same as the
! equivalently named eos_calc arguments.

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
! ifnuform: usually set to false so that h2, h2plus, h, hd, he, and
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
!   block), hne, sum0, sum2, s, extrasum, h, sh, hd, he, g,
!   hequil, and sumpl1 wrt f, t and sumpl0 wrt f with no other
!   variables fixed, i.e., use input dvf and dvt and chain rule.
!   N.B. chain of calling routines depend on the assertion that
!   ifnr = 0 means no *_dv variables should be read or written.
! ifnr = 1
!   calculate derivatives of nuvar (delivered through common
!   block), hne, sum0, sum2, extrasum, h, hd, and he wrt dv with
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
!   h, he, he+, c, c+, etc.  H2 and H2+ are last two.
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
! h, hf, ht, h_dv is n(H+) and its derivatives with respect to ln f, ln t, and dv.
! hd, hdf, hdt, hd_dv is n(He+) and its derivatives with respect to ln f, ln t, and dv.
! he, hef, het, he_dv is n(He++) and its derivatives with respect to ln f, ln t, and dv.
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

!> This eos_bfgs subroutine solves the EOS by minimizing the Helmholtz
!> free energy per unit mass as a function of auxiliary variables
!> using the BFGS technique.  N.B. for a given set of auxiliary
!> variables, fl must be adjusted so that rho is matched so this
!> minimization using auxiliary variable adjustment effectively occurs
!> for fixed rho.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine eos_bfgs(&
     verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
     ifrad, match_variable, kif, fl,&
     aux_old, aux, auxf, auxt, aux_dv,&
     njacobian,&
     faux, jacobian, p, pr,&
     dv_aux,&
     lambda, gamma_e,&
     ifnr, nux, nuy, nuz,&
     n_partial_ions,&
     sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
     partial_ions, t, morder, ifexchange_in,&
     full_sum0, full_sum1, full_sum2, charge, charge2,&
     dv_pl, dv_plt,&
     ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
     ifsame_zero_abundances,&
     ifexcited, ifsame_under,&
     inv_aux,&
     inv_ion, max_index,&
     partial_elements, n_partial_elements, ion_end, mion_end,&
     ifionized, if_pteh, if_mc, ifreducedmass,&
     ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,&
     ifpl, ifmodified, ifh2, ifh2plus,&
     izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
     eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
     r_ion3, r_neutral,&
     ifelement, dvzero, dv, dvf, dvt, nion,&
     ne, nef, net, sion, sionf, siont, uion,&
     pnorad, pnoradf, pnoradt, pnorad_aux,&
     free, free_aux,&
     nextrasum,&
     sumpl1, sumpl1f, sumpl1t, sumpl2,&
     h2, h2f, h2t, h2_dv, info)

  use mod_free_eos_constants, only: cr, ln10
  use mod_bfgs, only: bfgs_iterate, string_equal
  use mod_info_data, only: info_offset_eos_bfgs
  use mod_flow_data, only: ln_overflow_limit
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real, if_taint_allocated_integer, taint_allocated_integer
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit


  ! Arguments
  integer, intent(in) :: nion(:)
  integer, intent(in) :: ion_end(:)
  integer, intent(in) :: verbosity, ifreducedmass, ifmodified, ifexcited, nmax
  integer, intent(in) :: nmin(:), nmin_max(:), izlo, izhi
  integer, intent(in) :: nmin_species(:)
  integer, intent(in) :: iatomic_number(:)
  integer, intent(in) ::&
       ifnr,&
       ifh2, ifh2plus, ifmtrace,&
       if_pteh, ifcoulomb_mod, if_mc, if_dc,&
       ifionized, max_index,&
       nextrasum, ifpi_local, ifpl,&
       ifsame_abundances,&
       ifsame_zero_abundances
  integer, intent(in) :: ifdv(:), ifelement(:),&
       partial_ions(:), inv_ion(:),&
       n_partial_elements, partial_elements(:)
  integer, intent(in) :: inv_aux(:)
  integer, intent(in) :: partial_aux(:)
  integer, intent(in) :: mion_end, n_partial_ions
  integer, intent(in) :: kif, njacobian
  integer, intent(in) :: ifrad
  integer, intent(in) :: morder, ifexchange_in

  integer, intent(out) :: info

  logical, intent(in) :: ifsame_under, ifdvzero

  logical, intent(in) :: old_allow_log(:), old_allow_log_neg(:)
  logical, intent(out) :: old_new_allow_log(:), old_new_allow_log_neg(:)

  real(fp_kind), intent(in) :: bi(:), h2diss, t, eps(:)
  real(fp_kind), intent(in) ::&
       plop(:), plopt(:), plopt2(:),&
       dv_pl(:), dv_plt(:)
  real(fp_kind), intent(in) ::&
       full_sum0, full_sum1, full_sum2,&
       r_ion3(:), r_neutral(:),&
       charge(:), charge2(:),&
       tc2
  real(fp_kind), intent(in):: bmin(:)
  real(fp_kind), intent(in):: pr
  real(fp_kind), intent(in)::&
       sum0_mc, sum2_mc, hcon_mc, hecon_mc, nux, nuy, nuz

  real(fp_kind), intent(in):: match_variable, tl
  real(fp_kind), intent(out) :: lambda, gamma_e,&
       dvzero(:),&
       dv(:), dvf(:), dvt(:),&
       ne, nef, net,&
       sion, sionf, siont, uion,&
       h2, h2f, h2t, h2_dv(:),&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       aux(:), auxf(:), auxt(:), aux_dv(:,:),&
       faux(:), jacobian(:,:), p,&
       dv_aux(:,:),&
       pnorad, pnoradf, pnoradt, pnorad_aux(:),&
       free, free_aux(:)

  real(fp_kind), intent(inout) ::&
       aux_old(:),&
       fl

  ! Internal variables
  integer nions, nionsp2, naux, nelements, n_partial_aux
  integer linesearch_count, function_count, gradient_count
  integer index_aux, iaux
  integer nbfgs
  integer index_bfgs

  real(fp_kind) fl_acceptable
  real(fp_kind) fl_restore
  real(fp_kind) func
  real(fp_kind) fold, step_ratio, max_gradient_component

  character(:), allocatable :: bfgs_status
  character(:), allocatable :: bfgs_name

  logical, allocatable ::&
       allow_log_bfgs(:),&
       allow_log_neg_bfgs(:)

  integer, allocatable ::&
       transform_bfgs(:)

  real(fp_kind), allocatable ::&
       aux_acceptable(:),&
       aux_restore(:),&
       x(:),&
       xi(:),&
       dx(:),&
       gradient(:),&
       x_original(:)

  ! maximum core charge for non-bare ion
  integer, parameter :: maxcorecharge = 28
  integer, parameter  :: nxextrasum = 4
  integer, parameter :: nions_local = 316
  integer, parameter :: maxnextrasum_local = 9
  integer, parameter :: nxextrasum_local = 4
  ! each line search represents something like ~3 function and gradient
  ! calculations and therefore represents something like ~10 calls to
  ! eos_calc (due to the embedded fl iteration).  Thus, be careful that
  ! you do not specify too large an integer here since the number of line
  ! searches may often be the best termination criterion for the bfgs
  ! technique.
  integer, parameter :: max_linesearch_count = 20

  ! Setting fbar to a large negative number essentially disables the
  ! Fletcher bracket-limiting test which should be fine for well-behaved
  ! functions.  This parameter is in absolute units of func, but
  ! func = free/(full_sum0*RT) is expected to stay within a few orders
  ! of magnitude of unity.
  real(fp_kind), parameter :: fbar = -1.e200_fp_kind
  ! epsilon is a delta f expected from roundoff error to add
  ! some robustness in the presence of round-off error.  df_criterion
  ! is the function value convergence parameter for the bfgs iteration.
  ! gradient_criterion is the maximum absolute value of all gradient
  ! components before can exit from bfgs iteration.
  ! these three parameters are in absolute units of func, but
  ! func = free/(full_sum0*RT) is expected to stay within a few orders
  ! of magnitude of unity.
  ! These values may be fairly critical so some experimentation may be
  ! required for a large range of physical conditions to find the best
  ! values.
  real(fp_kind), parameter :: epsilon = 1.e-14_fp_kind
  real(fp_kind), parameter :: df_criterion = 1.e-14_fp_kind
  real(fp_kind), parameter :: gradient_criterion = 1.e-12_fp_kind

  if(verbosity.ge.3) then
     write(stderr,'(a,/,1p2e25.15e4)') 'entering eos_bfgs for '//&
          'log rho, log T = ',&
          match_variable/ln10, tl/ln10
  endif

  ! Assume success unless otherwise indicated by error condition below.
  info = 0

  nionsp2 = size(nion)
  nions = nionsp2 - 2
  naux = size(inv_aux)
  nelements = size(ion_end) - 2
  n_partial_aux = size(pnorad_aux) - 2

  ! sanity checks.
  if(nextrasum.gt.maxnextrasum_local) error stop 'eos_bfgs: nextrasum too large'
  if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3)) error stop 'eos_bfgs: ifnr must be 0, 1, or 3'
  if(ifnr.eq.0) error stop 'eos_bfgs: ifnr = 0 is disabled'
  if(nions.ne.nions_local) error stop 'eos_bfgs: nions must be equal to nions_local'
  if(size(partial_elements).gt.nelements+2) error stop 'eos_bfgs: partial_elements too large'
  if(n_partial_aux.gt.naux) error stop 'eos_bfgs: n_partial_aux too large'
  if(nxextrasum.gt.nxextrasum_local) error stop 'eos_bfgs: nxextrasum too large'
  if(ifrad.eq.2.and.kif.ne.1) error stop 'eos_bfgs: for ifrad == 2, kif values other than 1 are disabled.'
  if(kif.ne.2) error stop 'eos_bfgs: kif must be 2'

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
       nionsp2.ne.size(h2_dv).or.&
       nionsp2.ne.size(aux_dv,1).or.&
       nionsp2.ne.size(dv_aux,1))&
       error stop 'eos_bfgs: inconsistent nionsp2 array sizes'
  if(&
       naux.ne.size(partial_aux).or.&
       naux.ne.size(old_allow_log).or.&
       naux.ne.size(old_allow_log_neg).or.&
       naux.ne.size(old_new_allow_log).or.&
       naux.ne.size(old_new_allow_log_neg).or.&
       naux.ne.size(aux).or.&
       naux.ne.size(auxf).or.&
       naux.ne.size(auxt).or.&
       naux.ne.size(aux_dv,2).or.&
       naux.ne.size(faux).or.&
       naux.ne.size(jacobian,1).or.&
       naux.ne.size(jacobian,2).or.&
       naux.ne.size(dv_aux,2).or.&
       naux.ne.size(aux_old))&
       error stop 'eos_bfgs: inconsistent naux array sizes'
  if(&
       nelements.ne.size(iatomic_number).or.&
       nelements.ne.size(ifelement).or.&
       nelements.ne.size(eps).or.&
       nelements+2.ne.size(r_neutral).or.&
       nelements.ne.size(dvzero))&
       error stop 'eos_bfgs: inconsistent nelements array sizes'
  if(&
       maxcorecharge.ne.size(nmin).or.&
       maxcorecharge.ne.size(nmin_max).or.&
       maxcorecharge.ne.size(bmin))&
       error stop 'eos_bfgs: inconsistent maxcorecharge array sizes'
  if(&
       n_partial_aux+2.ne.size(free_aux))&
       error stop 'eos_bfgs: inconsistent n_partial_aux array sizes'

  allocate(&
       allow_log_bfgs(naux),&
       allow_log_neg_bfgs(naux),&
       transform_bfgs(naux),&
       aux_acceptable(naux),&
       aux_restore(naux),&
       x(naux),&
       xi(naux),&
       dx(naux),&
       gradient(naux),&
       x_original(naux))

  if(if_taint_allocated_real) then
     call taint_allocated_real(aux_acceptable)
     call taint_allocated_real(aux_restore)
     call taint_allocated_real(x)
     call taint_allocated_real(xi)
     call taint_allocated_real(dx)
     call taint_allocated_real(gradient)
     call taint_allocated_real(x_original)
  endif

  if(if_taint_allocated_integer) then
     call taint_allocated_integer(transform_bfgs)
  endif

  ! save values for later restoration.
  fl_restore = fl
  do index_aux = 1,n_partial_aux
     iaux = partial_aux(index_aux)
     aux_restore(iaux) = aux_old(iaux)
  enddo

  ! starting solution:
  index_bfgs = 0
  do index_aux = 1,n_partial_aux
     iaux = partial_aux(index_aux)
     ! ignore certain auxiliary variables for BFGS minimization.
     if(old_allow_log(index_aux).or.old_allow_log_neg(index_aux).or.iaux.eq.4) then
        index_bfgs = index_bfgs + 1
        transform_bfgs(index_bfgs) = iaux
        ! Save allow_log_bfgs and allow_log_neg_bfgs here which
        ! refer strictly to aux_old and which are indexed differently
        ! so that the relationship between xi and aux_old remains consistent.
        allow_log_bfgs(index_bfgs) = old_allow_log(index_aux)
        allow_log_neg_bfgs(index_bfgs) = old_allow_log_neg(index_aux)
        if(allow_log_bfgs(index_bfgs)) then
           x(index_bfgs) = log(aux_old(iaux))
        elseif(allow_log_neg_bfgs(index_bfgs)) then
           x(index_bfgs) = log(-aux_old(iaux))
        elseif(iaux.eq.4) then
           x(index_bfgs) = aux_old(iaux)
        else
           error stop 'eos_bfgs: logic error 1'
        endif
        x_original(index_bfgs) = x(index_bfgs)
     endif
  enddo
  nbfgs = index_bfgs
  bfgs_status = 'start'
  linesearch_count = 0
  function_count = 0
  gradient_count = 0
  ! Mark uninitialized fl_acceptable in case no bfgs step completes.
  fl_acceptable = 1.e300_fp_kind
  max_gradient_component = 2._fp_kind*gradient_criterion
  ! Force at least two iterations.
  func = 1.e300_fp_kind
  fold = 2._fp_kind*func
  do while(linesearch_count.le.max_linesearch_count.and.&
       fold.gt.func+df_criterion.and.&
       max_gradient_component.gt.gradient_criterion.and.&
       .not.string_equal(bfgs_status, 'error'))
     fold = func
     ! the combination of this bfgs_iterate and the completion of
     ! the following do while loop corresponds to the complete
     ! initialization of the bfgs technique (for linesearch_count.eq.0)
     ! or one line search followed by a bfgs correction in direction
     ! (linesearch_count.gt.0)
     call bfgs_iterate(verbosity, bfgs_status, bfgs_name, fbar, epsilon,&
          nbfgs, x, xi, func, gradient, step_ratio, dx)
     do while(.not.(string_equal(bfgs_status, 'error').or.string_equal(bfgs_status, 'complete')))
        ! A gradient calculation is always asked for right after a function
        ! calculation done at the same xi.
        if(string_equal(bfgs_status, 'function').or.string_equal(bfgs_status, 'both')) then
           do index_bfgs = 1,nbfgs
              iaux = transform_bfgs(index_bfgs)
              if(allow_log_bfgs(index_bfgs)) then
                 if(xi(index_bfgs).gt.ln_overflow_limit) then
                    if(verbosity.ge.1) write(stderr,*)&
                         "eos_bfgs: positive aux_old(iaux) for iaux = ", iaux, " is too close to overflowing"
                    info = info_offset_eos_bfgs + 1
                    return
                 endif
                 aux_old(iaux) = exp(xi(index_bfgs))
              elseif(allow_log_neg_bfgs(index_bfgs)) then
                 if(xi(index_bfgs).gt.ln_overflow_limit) then
                    if(verbosity.ge.1) write(stderr,*)&
                         "eos_bfgs: negative aux_old(iaux) for iaux = ", iaux, " is too close to overflowing"
                    info = info_offset_eos_bfgs + 2
                    return
                 endif
                 aux_old(iaux) = -exp(xi(index_bfgs))
              else
                 aux_old(iaux) = xi(index_bfgs)
              endif
           enddo
           if(.false..and.verbosity.ge.4) then
              write(stderr,'(a,/,(1p5e25.15e4))')&
                   'eos_bfgs: xi= ',&
                   (xi(index_bfgs), index_bfgs = 1,nbfgs)
              write(stderr,'(a,/,(1p5e25.15e4))')&
                   'eos_bfgs: aux_old= ',&
                   (aux_old(partial_aux(index_aux)),&
                   index_aux = 1,n_partial_aux)
           endif
           call eos_free_calc(&
                verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
                ifrad, match_variable, kif, fl,&
                aux_old, aux, auxf, auxt, aux_dv,&
                njacobian,&
                faux, jacobian, p, pr,&
                dv_aux,&
                lambda, gamma_e,&
                ifnr, nux, nuy, nuz,&
                n_partial_ions,&
                sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
                partial_ions, t, morder, ifexchange_in,&
                full_sum0, full_sum1, full_sum2, charge, charge2,&
                dv_pl, dv_plt,&
                ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
                ifsame_zero_abundances,&
                ifexcited, ifsame_under,&
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
                ne, nef, net, sion, sionf, siont, uion,&
                pnorad, pnoradf, pnoradt, pnorad_aux,&
                free, free_aux,&
                nextrasum,&
                sumpl1, sumpl1f, sumpl1t, sumpl2,&
                h2, h2f, h2t, h2_dv, info)
           if(info.ne.0) return
           ! Calculate func and gradient for the requested xi.
           func = free/(full_sum0*cr*t)
           do index_bfgs = 1,nbfgs
              iaux = transform_bfgs(index_bfgs)
              index_aux = inv_aux(iaux)
              if(index_aux.le.0) error stop 'eos_bfgs: logic error 2'
              gradient(index_bfgs) =&
                   free_aux(index_aux)/(full_sum0*cr*t)
              ! correct xi gradient for differences between
              ! allow_log_bfgs variables used to compute relationship
              ! between xi and aux_old and old_new_allow_log variables
              ! computed inside eos_jacobian to calculate free_aux.
              if(allow_log_bfgs(index_bfgs)) then
                 if(.not.(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)))&
                      gradient(index_bfgs) = gradient(index_bfgs)*aux_old(iaux)
              elseif(allow_log_neg_bfgs(index_bfgs)) then
                 if(.not.(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)))&
                      gradient(index_bfgs) = gradient(index_bfgs)*aux_old(iaux)
              elseif(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)) then
                 ! n.b. for this case aux_old must be non-zero (see old_new_allow_log
                 ! and old_new_allow_log_neg logic in eos_jacobian).
                 gradient(index_bfgs) = gradient(index_bfgs)/aux_old(iaux)
              endif
           enddo
        endif
        if(string_equal(bfgs_status, 'both')) then
           function_count = function_count + 1
           gradient_count = gradient_count + 1
           if(.false..and.verbosity.ge.4) then
              write(stderr,'(a,/,(1p5e25.15e4))') 'xi = ',&
                   (xi(index_aux), index_aux = 1,nbfgs)
              write(stderr,'(a,/,(1p5e25.15e4))') 'func = ', func
              write(stderr,'(a,/,(1p5e25.15e4))') 'gradient = ',&
                   (gradient(index_bfgs), index_bfgs = 1,nbfgs)
           endif
        elseif(string_equal(bfgs_status, 'function')) then
           function_count = function_count + 1
           if(.false..and.verbosity.ge.4) then
              write(stderr,'(a,/,(1p5e25.15e4))') 'xi = ',&
                   (xi(index_bfgs), index_bfgs = 1,nbfgs)
              write(stderr,'(a,/,(1p5e25.15e4))') 'func = ', func
           endif
        elseif(string_equal(bfgs_status, 'gradient')) then
           gradient_count = gradient_count + 1
           if(.false..and.verbosity.ge.4) then
              write(stderr,'(a,/,(1p5e25.15e4))') 'xi = ',&
                   (xi(index_bfgs), index_bfgs = 1,nbfgs)
              write(stderr,'(a,/,(1p5e25.15e4))') 'gradient = ',&
                   (gradient(index_bfgs), index_bfgs = 1,nbfgs)
           endif
        else
           error stop 'eos_bfgs: logic error 3'
        endif
        call bfgs_iterate(verbosity, bfgs_status, bfgs_name, fbar, epsilon,&
             nbfgs, x, xi, func, gradient, step_ratio, dx)
        ! end of do while(.not.(string_equal(bfgs_status, 'error').or.string_equal(bfgs_status, 'complete')))
     enddo
     if(string_equal(bfgs_status, 'complete')) then
        ! save latest acceptable fl and aux_old (signalled by 'complete'
        ! bfgs_status code) found by BFGS technique.
        ! n.b. this logic only works because the last evaluation done
        ! by the bfgs line search is always an acceptable point.
        fl_acceptable = fl
        do index_aux = 1,n_partial_aux
           iaux = partial_aux(index_aux)
           aux_acceptable(iaux) = aux_old(iaux)
        enddo
        if(verbosity.ge.3) then
           write(stderr,'(a,/,3i5)')&
                'linesearch_count, function_count, gradient_count =',&
                linesearch_count, function_count, gradient_count
           write(stderr,'(a,/,(1p5e25.15e4))') 'x = ',&
                (x(index_bfgs), index_bfgs = 1,nbfgs)
           write(stderr,'(a,/,(1p5e25.15e4))') 'x - x_original = ',&
                (x(index_bfgs)-x_original(index_bfgs),&
                index_bfgs = 1,nbfgs)
           write(stderr,'(a,/,(1p5e25.15e4))') 'fl = ', fl
           write(stderr,'(a,/,(1p5e25.15e4))') 'fl - fl_original = ',&
                fl-fl_restore
           if(linesearch_count.gt.0) then
              write(stderr,'(a,/,(1p5e25.15e4))') 'dx = ',&
                   (dx(index_bfgs), index_bfgs = 1,nbfgs)
              write(stderr,'(a,/,1pe25.15e4)') 'step_ratio = ', step_ratio
           else
              write(stderr,'(a,/,(1p5e25.15e4))') 'dx = ',&
                   (0._fp_kind, index_bfgs = 1,nbfgs)
              write(stderr,'(a,/,1pe25.15e4)') 'step_ratio = ', 0._fp_kind
           endif
           write(stderr,'(a,/,1pe25.15e4)') 'func = ', func
        endif
        max_gradient_component = 0._fp_kind
        do index_bfgs = 1, nbfgs
           max_gradient_component = max(max_gradient_component, abs(gradient(index_bfgs)))
        enddo
        if(verbosity.ge.3) then
           write(stderr,'(a,/,1pe25.15e4)') 'func - fold = ', func-fold
           write(stderr,'(a,/,1pe25.15e4)')&
                'maximum absolute value of gradient component = ',&
                max_gradient_component
           write(stderr,'(a,/,(1p5e25.15e4))') 'gradient = ',&
                (gradient(index_bfgs), index_bfgs = 1,nbfgs)
        endif
        ! increment at end of loop since linesearch = 0 corresponds to
        ! bfgs_status = 'start' which does initialization calculations at
        ! initial x.
        linesearch_count = linesearch_count + 1
     else
        if(verbosity.ge.1)&
             write(stderr,*) 'eos_bfgs: unexpected bfgs_status = ', bfgs_status, ', but continue anyway'
     endif
     ! End of
       !do while(linesearch_count.le.max_linesearch_count.and.&
       !     fold.gt.func+df_criterion.and.&
       !     max_gradient_component.gt.gradient_criterion.and.&
       !     .not.string_equal(bfgs_status, 'error'))
  enddo

  ! return last acceptable point to calling routine.
  if(fl_acceptable.ne.1.e300_fp_kind) then
     ! fl_acceptable has been initialized with at least one complete bfgs iteration
     fl = fl_acceptable
     do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
        aux_old(iaux) = aux_acceptable(iaux)
     enddo
  else
     if(verbosity.ge.1) write(stderr,*) 'eos_bfgs: no bfgs iteration completed'
     info = info_offset_eos_bfgs + 3
     return
  endif
  !temporary
  ! !restore input
  ! fl = fl_restore
  ! do index_aux = 1,n_partial_aux
  !   iaux = partial_aux(index_aux)
  !   aux_old(iaux) = aux_restore(iaux)
  ! enddo
  if(verbosity.ge.3) then
     write(stderr,*) 'Successful return from eos_bfgs'
  endif

end subroutine eos_bfgs
