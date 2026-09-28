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

! The purpose of this subroutine is to calculate auxiliary variables
! and their derivatives that are associated with the free energy of
! excited states.
! The excitation free energy (per unit volume) model is
! f = -kT sum n(species) delta ln Z, where
! delta ln Z = ln(1 + gion*exp(-c2 eion/T)*qstar/(g*plop*piop)),
! where:
! gion is the internal partition function of the ion of the species
! g is the internal partition function of the species
! qstar is the excited state partition function (relative to the
!   energy of the ion) in the hydrogenic approximation corrected for
!   the Planck-Larkin and MHD occupation probabilities
!   = sum (from n=nmin to nmax) 2 n^2 plop(n)*piop(n)* exp[arg(n)];
! arg(n) = c2 R iz^2/(n^2 T);
! R is Rydbergs constant
! plop = is the ground state Planck-Larkin occupation probability
! plop(n) = 1 - exp(-arg)*(1 + arg) (see notes for qstar_calc);
! piop is the ground state MHD occupation probability;
! piop(n) is the same thing for the excited states (see notes for qstar_calc).

! input quantities:
! ifexcited > 0 means use excited states (must have Planck-Larkin or ifpi = 3 or 4).
!    0 < ifexcited < 10 means use approximation to explicit summation
!   10 < ifexcited < 20 means use explicit summation
!   mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules) and helium.
!   mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
!   mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
! ifpl = 1, calculate Planck-Larkin occupation probability, otherwise set to unity;
! ifpi = 3 or 4, calculate MHD occupation probability, otherwise set to unity;
! ifmodified > 0 (or not) affects ifpi_fit
! ifnr = 0
!   calculate derivatives of psum and usum (delivered through common
!   block) and xextrasum wrt f, t with no other variables fixed,
!   i.e., use input f and t derivatives and chain rule.
!   N.B. chain of calling routines depend on the assertion that ifnr = 0
!   means no *_dv variables should be read or written.
! ifnr = 1
!   calculate derivatives of xextrasum wrt dv with f, t fixed.
!   n.b. this variable is subset of the list for ifnr = 0
!   because we only need dv derivatives for variables used to calculate
!   (output) auxiliary variables.
! ifnr = 2 (should not occur)
! ifnr = 3 same as combination of ifnr = 0 and ifnr = 1.
!   Note this is quite different in detail from ifnr = 3 interpretation
!   for many other routines, but the general motivation is the same
!   for all ifnr = 3 results; calculate both f and t derivatives and
!   other derivatives.
! inv_ion(nions+2) maps ion index to contiguous ion index used in NR iteration dot products.
! ifh2 = 0  no h2
! ifh2 = 1  vdb h2
! ifh2 = 2  st h2
! ifh2 = 3  irwin h2 (recommended)
! ifh2 = 4  pteh h2
! ifh2plus = 0  no h2plus
! ifh2plus = 1  st h2plus
! ifh2plus = 2  irwin h2plus (recommended)
! partial_elements(n_partial_elements+2) index of elements treated as
!   partially ionized consistent with ifelement.
! ion_end(nelements) keeps track of largest ion index for each element
! tl = log(t).
! izlo is lowest core charge (1 for H, He, 2 for He+, etc.) for all species.
! izhi is highest core charge for all species.
! bmin(izhi) is the minimum bion for all species with a given core charge.
! nmin(izhi) is the minimum excited principal quantum number for
!   all species with a given core charge.
! nmin_max(izhi) is the largest minimum excited principal
!   quantum number for all species with a given core charge.
! nmin_species(nions+2) is the minimum excited principal quantum number organized by species.
! nmax is maximum principal quantum number included in sum
!   (to be compatible with opal which used nmax = 4 rather than infinity).
!   if(nmax > 300000) then treated as infinity in qryd_approx.
!   otherwise nmax is meant to be used with qryd_calc only (i.e.,
!   case for mhd approximations not programmed).
! bion(nions+2) = ionization energy of next higher ion relative to species
!   that is being calculated.  nions+1 refers to H2, nions+2 refers to H2+.
! plop(nions+2), plopt(nions+2), plopt2(nions+2) = *ln* of ground
!   state Planck-Larkin occupation probabilty and first and second tl derivatives
! r_ion3(nions+2) is *the cube of the*
!   effective radii (in ion order but going from neutral
!   to next to bare ion for each species) of MDH interaction between
!   all but bare nucleii species and ionized species.  last 2 are
!   H2 and H2+ (note, only need r_ion3 of particular species)
! nion(nions+2), charge on ion in ion order (must be same order as bi)
!   e.g., for H+, He+, He++, etc.
! r_neutral(nelements+2) effective radii for MDH neutral-neutral
!   interactions.   Last two are H2 and H2+ (the only ionic
!   species in the MHD model with a non-zero hard-sphere radius).
! extrasum(nextrasum = 9) weighted sums over n(i)
!   for iextrasum = 1,nextrasum-2,
!   sum is only over neutral species (and H2+) and
!   weight is r_neutral^{iextrasum-1}
!   for iextrasum = nextrasum-1, sum is over all ionized species including
!   bare nucleii, but excluding free electrons, weight is Z^1.5.
!   for iextrasum = nextrasum, sum is over all species excluding bare nucleii
!   and free electrons, the weight is rion^3.
! extrasumf(nextrasum) = partial of extrasum/partial ln f
! extrasumt(nextrasum) = partial of extrasum/partial ln t
! extrasum_dv(nions+2,nextrasum) = partial derivatives wrt dv
! max_index
! nug
! nuh2
! nuh2plus
! xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
!   partial qratio/partial extrasum(k), k = 1, 2, 3, and nextrasum-1.

!> This excitation_sum subroutine calculates auxiliary variables that
!> are associated with the excited state pressure ionization component
!> of the free energy as well as partial derivatives of those
!> quantities wrt fl, tl, and dv.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine excitation_sum(verbosity, ifexcited, ifsame_zero_abundances,&
     ifpl, ifpi, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,&
     partial_elements, ion_end,&
     tl, izlo, bmin, nmin, nmin_max, nmin_species, nmax,&
     bion, plop, plopt, plopt2,&
     r_ion3, nion, r_neutral,&
     extrasum, extrasumf, extrasumt, extrasum_dv,&
     max_index, nug, nugf, nugt, nug_dv,&
     nuh2, nuh2f, nuh2t, nuh2_dv,&
     nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv,&
     xextrasum, xextrasumf, xextrasumt, xextrasum_dv)

  use mod_excitation_block, only: extrace_count, extrace_tag, extrace_ids, &
       extrace_value, extrace_grad, extrace_hess, extrace_hf, extrace_ht, extrace_scale, &
       extrace_mode, extrace_nr, extrace_moments
  use mod_free_eos_constants, only: pi, c2, electron_mass, h_mass
  use mod_molecular_hydrogen, only: molecular_hydrogen
  use mod_nuvar, only: nuvar, nuvarf, nuvart, nuvar_dv
  use mod_aux_scale, only: xextrasum_scale
  use mod_excitation_block, only: nx, c2t, free_sum, free_sumf, free_sum_dv, ifapprox_old_excitation,&
       ifdiff_x_excitation, ifh2_old_excitation, ifh2plus_old_excitation, ifhe1_special,&
       psum, psumf, psumt, psum_dv, qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,&
       qmhd_he1, qmhd_he1t, qmhd_he1x, qmhd_he1t2, qmhd_he1tx, qmhd_he1x2,&
       ssum, ssumf, ssumt, tl_old_excitation, usum,&
       x, x_old_excitation, qstar, qstart, qstarx, qstart2, qstart2, qstart2, qstartx, qstarx2,&
       ifmhd_logical_old_excitation, ifpi_fit_old_excitation, ifpl_logical_old_excitation,&
       max_nmin_max, nions_excitationp2,&
       ifnr03, ifnr13, excitation_sum_called
  use mod_helium1_data, only: ell_helium1, neff_helium1, weight_helium1
  use mod_statistical_weight_data, only: iqion, iqneutral
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  integer, intent(in) :: verbosity, ifexcited, ifsame_zero_abundances,&
       ifpl, ifpi, ifmodified,&
       ifnr, ifh2, ifh2plus,&
       partial_elements(:),&
       ion_end(:),&
       izlo, nmin(:), nmin_max(:), nmax,&
       nion(:), nmin_species(:),&
       inv_ion(:), max_index
  real(fp_kind), intent(in) :: tl, bmin(:), bion(:),&
       plop(:), plopt(:), plopt2(:),&
       r_ion3(:), r_neutral(:),&
       extrasum(:), extrasumf(:), extrasumt(:), extrasum_dv(:,:),&
       nug, nugf, nugt, nug_dv(:),&
       nuh2, nuh2f, nuh2t, nuh2_dv(:),&
       nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv(:)
  real(fp_kind), intent(out) :: xextrasum(:), xextrasumf(:),xextrasumt(:), xextrasum_dv(:,:)

  ! Parameters
  ! number of xextrasum variables
  integer, parameter :: nxextrasum = 4
  ! minimum value of principal quantum number where approximations are
  ! used
  integer, parameter :: min_nmin_max = 3
  ! minimum rydberg level for neutral helium (lower levels treated with
  ! exact neff, weight, and ell).
  ! 4 is indistinguishable from 10 for solar comparisons.
  integer, parameter :: nmin_rydberg_he1 = 4

  ! For what it is worth, this value corresponds roughly to ~1.d-217
  real(fp_kind), parameter :: exparg_lim = -500._fp_kind
  real(fp_kind), parameter :: exp_lim = exp(exparg_lim)
  ! corresponds to ~1.d-300
  real(fp_kind), parameter :: eps_factor_lim = -690._fp_kind
  ! This seems like a reasonable limit on the ln (ground state occupation
  ! probability) for when the ground state is worth correcting for
  ! excitation.
  real(fp_kind), parameter :: ground_state_lim = -50._fp_kind
  real(fp_kind), parameter :: occ_const = -4._fp_kind*pi/3._fp_kind
  real(fp_kind), parameter :: exparg_shift = 600._fp_kind
  ! multiply qratio and its derivatives by
  ! qratio_scale = exp(exparg_shift)
  ! in order to get more dynamic range.
  real(fp_kind), parameter :: qratio_scale = exp(exparg_shift)

  ! Local variables
  logical ifpl_logical, ifmhd_logical, ifneutral, ifapprox

  integer n_partial_elements
  integer izhi, nions, nextrasum
  integer iz, izqstar, index, ielement, ion_start, ion, nmin_s, ifpi_fit
  integer nelements, nionsp2, ix, nmin_max_max

  integer index_x(nx)
  ! Use data statement rather than parameter because last element needs to be set below at run time.
  data index_x /4, 3, 2, 1, 0/
  integer iffirst
  data iffirst/1/

  real(fp_kind) exparg,&
       eps_factor, rfactor,&
       qratio, qratiot, qratiot2,&
       qratiov, qratiovt,&
       dqratiof, dqratiotf,&
       dqratiot, dqratiott,&
       dqratiovf, dqratiovt
  real(fp_kind) qratio_us, onepq, onepqs
  real(fp_kind) ln_ground_occ

  real(fp_kind), allocatable ::&
       qratio_dx(:),&
       qratiot_dx(:),&
       qratio_dx2(:,:),&
       qratiov_dx(:),&
       dqratio_dxf(:),&
       dqratio_dxt(:),&
       dqratio_dv(:),&
       dqratio_dxdv(:,:),&
       dqratiov_dv(:)

  ! Most/all Fortran compilers specify the save attribute for all variables
  ! intialized by data statements.  But just in case...
  save iffirst, index_x

  excitation_sum_called = .true.

  if(iffirst.eq.1) then
     ! go through this just in case this routine is called without
     ! excitation_pi called first.  In normal case this means we
     ! have an extra evaluation of qh2, qh2plus, and qstar on first
     ! call, but from then on, no extra calls at all.
     iffirst = 0
     ! Something non-astrophysical
     tl_old_excitation = -1.e30_fp_kind
     ! This required so that on first call
     ! x_old_excitation doesn't need to be initialized,
     ! ifdiff_x_excitation = .true., and
     ! if(ifdiff_x_excitation.or.&... test is always
     ! true regardless of whether the *_old_excitation
     ! variables below are initialized or not.
     ifmhd_logical_old_excitation = .false.
     ! Initialization of these "old_excitation" variables
     ! shouldn't matter (see above argument), but on
     ! the principle that the if(ifdiff_x_excitation.or.&... test
     ! below should not depend on uninitialized variables set
     ! these values to anything.
     ifpi_fit_old_excitation = 0
     ifh2_old_excitation = 0
     ifh2plus_old_excitation = 0
     ifpl_logical_old_excitation = .false.
     ifapprox_old_excitation = .false.
  endif

  n_partial_elements = size(partial_elements) - 2
  izhi = size(bmin)
  nelements = size(r_neutral) - 2
  nionsp2 = size(nion)
  nions = nionsp2 - 2
  nextrasum = size(extrasum)
  ifpl_logical = ifpl.eq.1
  ifmhd_logical = ifpi.eq.3.or.ifpi.eq.4

  ! Sanity checks
  if(nionsp2.ne.nions_excitationp2) error stop 'excitation_sum: inconsistent nionsp2 values'
  if(ifnr.lt.0.or.ifnr.eq.2.or.ifnr.gt.3) error stop 'excitation_sum: ifnr must be 0, 1, or 3'
  if(ifexcited.le.0) error stop 'excitation_sum: invalid ifexcited'
  if(.not.(ifpl_logical.or.ifmhd_logical)) error stop 'excitation_sum: invalid combination of ifpl and ifpi'
  if(.not.((ifmhd_logical.and.nextrasum.eq.9).or.(.not.ifmhd_logical.and.nextrasum.eq.0)))&
       error stop 'excitation_sum: inconsistent ifpi and nextrasum'

  if(&
       nx.ne.size(x).or.&
       nx.ne.size(x_old_excitation)&
       )&
       error stop 'excitation_sum: inconsistent sizes for x and x_old_excitation'
  if(nelements+2.ne.size(ion_end)) error stop 'excitation_sum: inconsistent sizes for ion_end and r_neutral'
  if(izhi.ne.size(nmin).or.izhi.ne.size(nmin_max))&
       error stop 'excitation_sum: inconsistent sizes for bmin, nmin, or nmin_max'
  if(&
       nionsp2.ne.size(nmin_species).or.&
       nionsp2.ne.size(inv_ion).or.&
       nionsp2.ne.size(bion).or.&
       nionsp2.ne.size(plop).or.&
       nionsp2.ne.size(plopt).or.&
       nionsp2.ne.size(plopt2).or.&
       nionsp2.ne.size(r_ion3).or.&
       nionsp2.ne.size(extrasum_dv,1).or.&
       nionsp2.ne.size(nug_dv).or.&
       nionsp2.ne.size(nuh2_dv).or.&
       nionsp2.ne.size(nuh2plus_dv).or.&
       nionsp2.ne.size(xextrasum_dv,1))&
       error stop 'excitation_sum: inconsistent nionsp2 sizes for 13 variables'
  if(&
       nextrasum.ne.size(extrasumf).or.&
       nextrasum.ne.size(extrasumt).or.&
       nextrasum.ne.size(extrasum_dv,2))&
       error stop 'excitation_sum: inconsistent sizes for extrasum, extrasumt, extrasumf, or extrasum_dv'
  if(&
       nxextrasum.ne.size(xextrasum).or.&
       nxextrasum.ne.size(xextrasumf).or.&
       nxextrasum.ne.size(xextrasumt).or.&
       nxextrasum.ne.size(xextrasum_dv,2))&
       error stop 'excitation_sum: inconsistent sizes for xextrasum, xextrasumt, xextrasumf, or xextrasum_dv'

  allocate(&
       qratio_dx(nxextrasum),&
       qratiot_dx(nxextrasum),&
       qratio_dx2(nxextrasum,nxextrasum),&
       qratiov_dx(nxextrasum),&
       dqratio_dxf(nxextrasum),&
       dqratio_dxt(nxextrasum),&
       dqratio_dv(nionsp2),&
       dqratio_dxdv(nionsp2,nxextrasum),&
       dqratiov_dv(nionsp2)&
       )

  if(if_taint_allocated_real) then
     call taint_allocated_real(qratio_dx)
     call taint_allocated_real(qratiot_dx)
     call taint_allocated_real(qratio_dx2)
     call taint_allocated_real(qratiov_dx)
     call taint_allocated_real(dqratio_dxf)
     call taint_allocated_real(dqratio_dxt)
     call taint_allocated_real(dqratio_dv)
     call taint_allocated_real(dqratio_dxdv)
     call taint_allocated_real(dqratiov_dv)
  endif

  ifnr03 = ifnr.eq.0.or.ifnr.eq.3
  ifnr13 = ifnr.eq.1.or.ifnr.eq.3
  ifapprox = ifexcited.lt.10
  c2t = c2*exp(-tl)
  ! sort out what fitting factors will be applied to pressure ionization.
  if(ifpi.eq.4) then
     ! factors to fit Saumon results.
     ifpi_fit = 2
  elseif(ifmodified.gt.0) then
     ! factors to fit opal results.
     ifpi_fit = 1
  else
     ! unity factors to mimic mdh results as closely as possible.
     ifpi_fit = 0
  endif

  ! Start of blocks of code that should be identical (aside from error
  ! stop identifications) for excitation_pi and excitation_sum.  The
  ! idea is all the qh2, qh2plus, and qstar values will be taken from
  ! previous calculations (either excitation_pi or excitation_sum) if
  ! nothing has been changed.

  ! ifdiff_x_excitation is true if ifmhd_logical is true and
  ! the old version not or if both true and the old version has
  ! different x values.
  ! also initialize x and x_old_excitation if needed.
  if(ifmhd_logical) then
     ! index_x(nx) is a variable value so cannot be set in the data statement above.
     index_x(nx) = nextrasum-1
     do ix = 1,nx
        x(ix) = extrasum(index_x(ix))
     enddo
     if(ifmhd_logical_old_excitation) then
        ix = 1
        do while(ix.lt.nx.and.x(ix).eq.x_old_excitation(ix))
           ix = ix + 1
        enddo
        ifdiff_x_excitation = x(ix).ne.x_old_excitation(ix)
     else
        ifdiff_x_excitation = .true.
     endif
     x_old_excitation = x
  else
     ifdiff_x_excitation = .false.
  endif

  ! On initial call to this routine, tl_old_excitation is initialized to
  ! something non-astrophysical (see iffirst code above)
  ! so this entire Boolean logic expression should be true
  ! regardless of the values of the other old_excitation variables.
  if(&
       tl.ne.tl_old_excitation.or.&
       ifdiff_x_excitation.or.&
       (ifmhd_logical.neqv.ifmhd_logical_old_excitation).or.&
       ifpi_fit.ne.ifpi_fit_old_excitation.or.&
       (ifpl_logical.neqv.ifpl_logical_old_excitation).or.&
       (ifapprox.neqv.ifapprox_old_excitation).or.&
       ifsame_zero_abundances.ne.1) then
     ifpi_fit_old_excitation = ifpi_fit
     ifpl_logical_old_excitation = ifpl_logical
     ifmhd_logical_old_excitation = ifmhd_logical
     ifapprox_old_excitation = ifapprox

     ! On initial call to this routine, tl_old_excitation is initialized to
     ! something non-astrophysical (see iffirst code above)
     ! so for this case the first Boolean logic block below should be
     ! true regardless of how the other old_excitation variables have been initialized
     ! (on the principle that in Fortran all variables used in if statements should
     ! be initialized).
     if(&
          (&
          tl.ne.tl_old_excitation.or.&
          ifh2.ne.ifh2_old_excitation.or.&
          ifh2plus.ne.ifh2plus_old_excitation&
          ).and.&
          (ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1))&
          call molecular_hydrogen(verbosity, ifh2, ifh2plus, tl,&
          qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2)
     ifh2_old_excitation = ifh2
     ifh2plus_old_excitation = ifh2plus
     tl_old_excitation = tl
     do iz = izlo, izhi
        if(nmin_max(iz).gt.max_nmin_max) error stop 'excitation_sum: nmin_max too large'
        ! eps_factor determines the limit on the excited-state
        ! partition function sum.  (Include summands beyond nmin where
        ! summand > eps/eps_factor where the error due to the
        ! cutoff is roughly eps = 1.d-10.)
        ! N.B. the limit on ln(eps_factor) avoids underflows
        ! and divide by zeros in the summand test.  Also, for such near-zero
        ! eps_factor values, the partition function sum never passes the
        ! summand test so is just given by the minimum nmin value.
        ! N.B. for ground-state occupation probabilities near unity,
        ! eps_factor is always larger than epsarg below so one might
        ! be tempted to avoid the call to qstar_calc altogether if
        ! eps_factor less than exp_factor_lim.  However, that logic does
        ! not work for small ground-state occupation probabilities so
        ! use simple floor logic on eps_factor and always call
        ! qstar_calc.
        eps_factor = max(eps_factor_lim,-c2t*bmin(iz))
        eps_factor = exp(eps_factor)
        ifneutral = iz.eq.1
        nmin_max_max = max(min_nmin_max,nmin_max(iz))
        call qstar_calc(ifpi_fit,&
             ifpl_logical, ifmhd_logical, ifneutral, ifapprox,&
             eps_factor, nmin(iz), nmax, iz, tl, x,&
             qstar(:nmin_max_max,iz), qstart(:nmin_max_max,iz),&
             qstarx(:,:nmin_max_max,iz),&
             qstart2(:nmin_max_max,iz), qstartx(:,:nmin_max_max,iz),&
             qstarx2(:,:,:nmin_max_max,iz),qstar_logscale(:nmin_max_max,iz))
        ! ! zero results (including derivatives below) to provide
        ! ! consistent small value zeroing rather than inconsistent
        ! ! underflow zeroing.
        ! do n = nmin(iz), nmin_max(iz)
        !   if(qstar(n,iz).lt.exp_lim) qstar(n,iz) = 0.d0
        ! enddo
        if(iz.eq.2.and.ifh2plus.gt.0) then
           ifneutral = .true.
           ! special values of qstar and friends calculated including
           ! "neutral" occupation probability for H2+.
           nmin_max_max = max(min_nmin_max,nmin_max(iz))
           call qstar_calc(ifpi_fit,&
                ifpl_logical, ifmhd_logical, ifneutral, ifapprox,&
                eps_factor, nmin(iz), nmax, iz, tl, x,&
                qstar(:nmin_max_max,izhi+1), qstart(:nmin_max_max,izhi+1),&
                qstarx(:,:nmin_max_max,izhi+1),&
                qstart2(:nmin_max_max,izhi+1), qstartx(:,:nmin_max_max,izhi+1),&
                qstarx2(:,:,:nmin_max_max,izhi+1),qstar_logscale(:nmin_max_max,izhi+1))
           ! ! zero results (including derivatives below) to provide
           ! ! consistent small value zeroing rather than inconsistent
           ! ! underflow zeroing.
           ! do n = nmin(iz), nmin_max(iz)
           !   if(qstar(n,izhi+1).lt.exp_lim) qstar(n,izhi+1) = 0.d0
           ! enddo
        ! else ! if(iz.eq.2.and.ifh2plus.gt.0) then
        !   do n = nmin(iz), nmin_max(iz)
        !     qstar(n,iz) = 0.d0
        !   enddo
        !   if(iz.eq.2.and.ifh2plus.gt.0) then
        !     do n = nmin(iz), nmin_max(iz)
        !       qstar(n,izhi+1) = 0.d0
        !     enddo
        !   endif
        endif ! if(iz.eq.2.and.ifh2plus.gt.0) then
     enddo ! do iz = izlo, izhi
  endif !test for previous calculation of qh2, qh2plus, and qstar
  ! End of blocks of code that should be identical (aside from error
  ! stop identifications) for excitation_pi and excitation_sum.

  psum = 0._fp_kind
  ssum = 0._fp_kind
  usum = 0._fp_kind
  ! free_sum = ssum - usum proportional to negative of free energy
  ! associated with excitation.
  free_sum = 0._fp_kind
  extrace_count=0;extrace_ids=0;extrace_value=0._fp_kind
  extrace_grad=0._fp_kind;extrace_hess=0._fp_kind
  extrace_hf=0._fp_kind;extrace_ht=0._fp_kind
  extrace_mode=ifexcited;extrace_nr=ifnr
  extrace_moments(:,1)=extrasum([1,2,3,nextrasum-1])
  extrace_moments(:,2)=extrasumf([1,2,3,nextrasum-1])
  extrace_moments(:,3)=extrasumt([1,2,3,nextrasum-1])
  if(ifmhd_logical) then
     xextrasum(1:nxextrasum) = 0._fp_kind
  endif
  if(ifnr03) then
     psumf = 0._fp_kind
     psumt = 0._fp_kind
     ssumf = 0._fp_kind
     ssumt = 0._fp_kind
     free_sumf = 0._fp_kind
     if(ifmhd_logical) then
        xextrasumf(1:nxextrasum) = 0._fp_kind
        xextrasumt(1:nxextrasum) = 0._fp_kind
     endif
  endif
  if(ifnr13) then
     psum_dv(1:max_index) = 0._fp_kind
     free_sum_dv(1:max_index) = 0._fp_kind
  endif
  if(ifmhd_logical.and.ifnr13) then
     xextrasum_dv(1:max_index,1:nxextrasum) = 0._fp_kind
  endif
  ! This do loop is over atomic species only
  do index = 1, n_partial_elements
     ielement = partial_elements(index)
     if((ielement.gt.1.or.ifh2.eq.0).and.&
          (ielement.le.2.or.mod(ifexcited,10).eq.3)) then
        if(ielement.gt.1) then
           ion_start = ion_end(ielement-1) + 1
        else
           ion_start = 1
        endif
        ! calculate delta ln Z = ln(1 + qratio), where
        ! qratio is the ratio of excited to ground
        ! electronic state partition function
        ! n.b. index goes from neutral to next to bare ion
        do ion = ion_start, ion_end(ielement)
           iz = nion(ion)
           nmin_s = nmin_species(ion)
           if(nmin_s.lt.nmin(iz).or.nmin_s.gt.nmin_max(iz)) error stop 'excitation_sum: invalid nmin or nmin_max'
           exparg = exparg_shift - c2t*bion(ion) - plop(ion)
           ! correct for MHD ground state occupation probability.
           if(ifmhd_logical) then
              if(iz.eq.1) then
                 ln_ground_occ = occ_const*(&
                      x(1) + r_neutral(ielement)*(&
                      x(2)*3._fp_kind + r_neutral(ielement)*(&
                      x(3)*3._fp_kind + r_neutral(ielement)*(&
                      x(4)))))
              else
                 ln_ground_occ = 0._fp_kind
              endif
              ln_ground_occ = ln_ground_occ + occ_const*r_ion3(ion)*x(5)
              ! ignore excitation correction if ground state wiped out by
              ! pressure ionization in any case.
              if(ln_ground_occ.gt.ground_state_lim) then
                 exparg = exparg - ln_ground_occ
              else
                 ! signal to ignore this species
                 exparg = -1000._fp_kind!exparg_lim - 1._fp_kind
              endif
           endif
           ! if(exparg.gt.exparg_lim) then
           ! Keep the Boltzmann factor logarithmic until the partition product.
           ! else
           !   exparg = 0.d0
           ! endif
           ! qratio is first ratio of excited to ground partition function,
           ! but then is transformed to ln(1 + qratio).  Similarly, there
           ! is a subsequent transformation of all derivatives.
           ! N.B. important convention on partial derivative variable names:
           ! names starting with "qratio" are partial derivatives assuming
           ! that qratio is a function of tl and x.
           ! names starting with "dqratio" are the *change* to the partial
           ! derivative caused by x being a function of fl, tl, dv.
           if(iz.eq.1) then
              ! deal with neutral monatomic species of element
              ! ratio of excited to ground state partition functions
              if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
                 ! special for neutral helium
                 ! eps_factor determines the limit on the excited-state
                 ! partition function sum.  (Include summands beyond nmin where
                 ! summand > eps/eps_factor where the error due to the
                 ! cutoff is roughly eps = 1.d-10.)
                 ! N.B. the limit on ln(eps_factor) avoids underflows
                 ! and divide by zeros in the summand test.  Also, for
                 ! such near-zero eps_factor values, the partition
                 ! function sum never passes the
                 ! summand test so is just given by the minimum nmin value.
                 ! N.B. for ground-state occupation probabilities near unity,
                 ! eps_factor is the same as epsarg above so one might
                 ! be tempted to avoid the call to qstar_calc altogether if
                 ! eps_factor less than exp_factor_lim.  However, that logic does
                 ! not work for small ground-state occupation probabilities so
                 ! use simple floor logic on eps_factor and always call
                 ! qstar_calc.
                 eps_factor = max(eps_factor_lim,-c2t*bion(ion))
                 eps_factor = exp(eps_factor)
                 rfactor = 1._fp_kind/(1._fp_kind + electron_mass/(4._fp_kind*h_mass))
                 call qmhd_calc(ifpi_fit, real((iqion(ion)),fp_kind),&
                      ifpl_logical, ifapprox, eps_factor,&
                      nmin_rydberg_he1, nmax,&
                      tl, iz, rfactor, weight_helium1, neff_helium1, ell_helium1, x,&
                      qmhd_he1, qmhd_he1t, qmhd_he1x, qmhd_he1t2, qmhd_he1tx, qmhd_he1x2)
                 ! ! zero results (including derivatives below) to provide
                 ! ! consistent small value zeroing rather than inconsistent
                 ! ! underflow zeroing.
                 ! if(qmhd_he1.lt.exp_lim) qmhd_he1 = 0.d0
                 qratio = partition_product(exparg,qmhd_he1,0._fp_kind)/real((iqneutral(ielement)),fp_kind)
                 ! else
                 !   qratio = 0.d0
                 !  endif
              else
                 ! Not special case for helium
                 qratio = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))*real((iqion(ion)),fp_kind)/real((iqneutral(ielement)),fp_kind)
              endif
           else
              ! ionized atomic case
              qratio = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))*real((iqion(ion)),fp_kind)/real((iqion(ion-1)),fp_kind)
           endif
           if(ion.eq.1.and.qstar(nmin_s,iz).gt.0._fp_kind) native_h_log_ratio= &
                exparg+qstar_logscale(nmin_s,iz)+log(qstar(nmin_s,iz))+ &
                log(real(iqion(ion),fp_kind)/real(iqneutral(ielement),fp_kind))-exparg_shift
           if(qratio.gt.exp_lim) then
              qratio_us = qratio/qratio_scale
              onepq = 1._fp_kind + qratio_us
              onepqs = onepq*qratio_scale
              ! qratiot and qratiot2 are the first and second partial wrt tl
              if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
                 ! special for neutral helium
                 qratiot = qratio*(c2t*bion(ion) - plopt(ion) +&
                      qmhd_he1t/qmhd_he1)
                 qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +&
                      qmhd_he1t/qmhd_he1) +&
                      qratio*(-c2t*bion(ion) - plopt2(ion) +&
                      qmhd_he1t2/qmhd_he1 -&
                      (qmhd_he1t/qmhd_he1)*&
                      (qmhd_he1t/qmhd_he1))
              else
                 qratiot = qratio*(c2t*bion(ion) - plopt(ion) +&
                      qstart(nmin_s,iz)/qstar(nmin_s,iz))
                 qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +&
                      qstart(nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(-c2t*bion(ion) - plopt2(ion) +&
                      qstart2(nmin_s,iz)/qstar(nmin_s,iz) -&
                      (qstart(nmin_s,iz)/qstar(nmin_s,iz))*&
                      (qstart(nmin_s,iz)/qstar(nmin_s,iz)))
              endif
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiot2 = (qratiot2 - qratiot*&
                   (qratiot/onepqs))/onepq
              if(ifmhd_logical) then
                 ! qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are
                 ! partials wrt tl, extrasum(k) and extrasum(l) following
                 ! special convention for k or l = 4.
                 if(ifhe1_special.and.ion.eq.2) then
                    ! special for neutral helium
                    qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +&
                         qmhd_he1x(5)/qmhd_he1)
                    qratiot_dx(4) = qratiot*&
                         (-occ_const*r_ion3(ion) +&
                         qmhd_he1x(5)/qmhd_he1) +&
                         qratio*(qmhd_he1tx(5)/qmhd_he1-(qmhd_he1x(5)/qmhd_he1)*(qmhd_he1t/qmhd_he1))
                    qratio_dx2(4,4) = qratio_dx(4)*&
                         (-occ_const*r_ion3(ion) +&
                         qmhd_he1x(5)/qmhd_he1) +&
                         qratio*(qmhd_he1x2(5,5)/qmhd_he1-(qmhd_he1x(5)/qmhd_he1)*(qmhd_he1x(5)/qmhd_he1))
                 else
                    qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +&
                         qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
                    qratiot_dx(4) = qratiot*&
                         (-occ_const*r_ion3(ion) +&
                         qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +&
                         qratio*(qstartx(5,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                    qratio_dx2(4,4) = qratio_dx(4)*&
                         (-occ_const*r_ion3(ion) +&
                         qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +&
                         qratio*(qstarx2(5,5,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)))
                 endif
                 ! partial wrt ln V, tl, and extrasum(k) following special
                 ! k convention.
                 qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
                 qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
                 qratiov_dx(4) = -qratio_dx2(4,4)*&
                      extrasum(nextrasum-1) -&
                      qratio_dx(4)
                 ! transform second order quantities to derivative of
                 ! ln(1 + qratio)
                 qratiot_dx(4) = (qratiot_dx(4) - qratiot*&
                      (qratio_dx(4)/onepqs))/onepq
                 qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*&
                      (qratio_dx(4)/onepqs))/onepq
                 if(iz.eq.1) then
                    ! N.B. x(1) (or extrasum(4)) dependence divides out of qratio
                    ! n.b. indices are reordered here so that
                    ! qratio_dx(k) refers to derivative wrt extrasum(k) for
                    ! k = 1, 2, 3, and qratio_dx(4) refers to the derivative
                    ! wrt extrasum(nextrasum-1).
                    if(ifhe1_special.and.ion.eq.2) then
                       ! special for neutral helium
                       qratio_dx(1) = qratio*(&
                            -occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qmhd_he1x(4)/qmhd_he1)
                       qratio_dx(2) = qratio*(&
                            -occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qmhd_he1x(3)/qmhd_he1)
                       qratio_dx(3) = qratio*(&
                            -occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qmhd_he1x(2)/qmhd_he1)
                       qratiot_dx(1) = qratiot*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qmhd_he1x(4)/qmhd_he1) +&
                            qratio*(qmhd_he1tx(4)/qmhd_he1-(qmhd_he1x(4)/qmhd_he1)*(qmhd_he1t/qmhd_he1))
                       qratiot_dx(2) = qratiot*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qmhd_he1x(3)/qmhd_he1) +&
                            qratio*(qmhd_he1tx(3)/qmhd_he1-(qmhd_he1x(3)/qmhd_he1)*(qmhd_he1t/qmhd_he1))
                       qratiot_dx(3) = qratiot*&
                            (-occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qmhd_he1x(2)/qmhd_he1) +&
                            qratio*(qmhd_he1tx(2)/qmhd_he1-(qmhd_he1x(2)/qmhd_he1)*(qmhd_he1t/qmhd_he1))
                       qratio_dx2(1,1) = qratio_dx(1)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qmhd_he1x(4)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(4,4)/qmhd_he1-(qmhd_he1x(4)/qmhd_he1)*(qmhd_he1x(4)/qmhd_he1))
                       qratio_dx2(2,1) = qratio_dx(2)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qmhd_he1x(4)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(4,3)/qmhd_he1-(qmhd_he1x(4)/qmhd_he1)*(qmhd_he1x(3)/qmhd_he1))
                       qratio_dx2(3,1) = qratio_dx(3)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qmhd_he1x(4)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(4,2)/qmhd_he1-(qmhd_he1x(4)/qmhd_he1)*(qmhd_he1x(2)/qmhd_he1))
                       qratio_dx2(4,1) = qratio_dx(4)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qmhd_he1x(4)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(5,4)/qmhd_he1-(qmhd_he1x(5)/qmhd_he1)*(qmhd_he1x(4)/qmhd_he1))
                       qratio_dx2(2,2) = qratio_dx(2)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qmhd_he1x(3)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(3,3)/qmhd_he1-(qmhd_he1x(3)/qmhd_he1)*(qmhd_he1x(3)/qmhd_he1))
                       qratio_dx2(3,2) = qratio_dx(3)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qmhd_he1x(3)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(3,2)/qmhd_he1-(qmhd_he1x(3)/qmhd_he1)*(qmhd_he1x(2)/qmhd_he1))
                       qratio_dx2(4,2) = qratio_dx(4)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qmhd_he1x(3)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(5,3)/qmhd_he1-(qmhd_he1x(5)/qmhd_he1)*(qmhd_he1x(3)/qmhd_he1))
                       qratio_dx2(3,3) = qratio_dx(3)*&
                            (-occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qmhd_he1x(2)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(2,2)/qmhd_he1-(qmhd_he1x(2)/qmhd_he1)*(qmhd_he1x(2)/qmhd_he1))
                       qratio_dx2(4,3) = qratio_dx(4)*&
                            (-occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qmhd_he1x(2)/qmhd_he1) +&
                            qratio*(qmhd_he1x2(5,2)/qmhd_he1-(qmhd_he1x(5)/qmhd_he1)*(qmhd_he1x(2)/qmhd_he1))
                    else
                       qratio_dx(1) = qratio*(&
                            -occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                       qratio_dx(2) = qratio*(&
                            -occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                       qratio_dx(3) = qratio*(&
                            -occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                       qratiot_dx(1) = qratiot*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstartx(4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                       qratiot_dx(2) = qratiot*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstartx(3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                       qratiot_dx(3) = qratiot*&
                            (-occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstartx(2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(1,1) = qratio_dx(1)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(4,4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(2,1) = qratio_dx(2)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(4,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(3,1) = qratio_dx(3)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(4,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(4,1) = qratio_dx(4)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*r_neutral(ielement) +&
                            qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(5,4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(2,2) = qratio_dx(2)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(3,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(3,2) = qratio_dx(3)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(3,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(4,2) = qratio_dx(4)*&
                            (-occ_const*r_neutral(ielement)*&
                            r_neutral(ielement)*3._fp_kind +&
                            qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(5,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(3,3) = qratio_dx(3)*&
                            (-occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(2,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                       qratio_dx2(4,3) = qratio_dx(4)*&
                            (-occ_const*r_neutral(ielement)*&
                            3._fp_kind +&
                            qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                            qratio*(qstarx2(5,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                    endif !if(ifhe1_special.and.ion.eq.2) ...
                    ! derivative of qratio wrt ln V
                    qratiov = qratiov -&
                         qratio_dx(1)*extrasum(1) -&
                         qratio_dx(2)*extrasum(2) -&
                         qratio_dx(3)*extrasum(3)
                    qratiovt = qratiovt -&
                         qratiot_dx(1)*extrasum(1) -&
                         qratiot_dx(2)*extrasum(2) -&
                         qratiot_dx(3)*extrasum(3)
                    qratiov_dx(1) = -&
                         qratio_dx2(1,1)*extrasum(1) -&
                         qratio_dx2(2,1)*extrasum(2) -&
                         qratio_dx2(3,1)*extrasum(3) -&
                         qratio_dx2(4,1)*extrasum(nextrasum-1) -&
                         qratio_dx(1)
                    qratiov_dx(2) = -&
                         qratio_dx2(2,1)*extrasum(1) -&
                         qratio_dx2(2,2)*extrasum(2) -&
                         qratio_dx2(3,2)*extrasum(3) -&
                         qratio_dx2(4,2)*extrasum(nextrasum-1) -&
                         qratio_dx(2)
                    qratiov_dx(3) = -&
                         qratio_dx2(3,1)*extrasum(1) -&
                         qratio_dx2(3,2)*extrasum(2) -&
                         qratio_dx2(3,3)*extrasum(3) -&
                         qratio_dx2(4,3)*extrasum(nextrasum-1) -&
                         qratio_dx(3)
                    qratiov_dx(4) = qratiov_dx(4) -&
                         qratio_dx2(4,1)*extrasum(1) -&
                         qratio_dx2(4,2)*extrasum(2) -&
                         qratio_dx2(4,3)*extrasum(3)
                    ! transform second order quantities to derivative of
                    ! ln(1 + qratio)
                    qratiot_dx(1) = (qratiot_dx(1) - qratiot*&
                         (qratio_dx(1)/onepqs))/onepq
                    qratiot_dx(2) = (qratiot_dx(2) - qratiot*&
                         (qratio_dx(2)/onepqs))/onepq
                    qratiot_dx(3) = (qratiot_dx(3) - qratiot*&
                         (qratio_dx(3)/onepqs))/onepq
                    qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*&
                         (qratio_dx(1)/onepqs))/onepq
                    qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*&
                         (qratio_dx(1)/onepqs))/onepq
                    qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*&
                         (qratio_dx(1)/onepqs))/onepq
                    qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*&
                         (qratio_dx(1)/onepqs))/onepq
                    qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*&
                         (qratio_dx(2)/onepqs))/onepq
                    qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*&
                         (qratio_dx(2)/onepqs))/onepq
                    qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*&
                         (qratio_dx(2)/onepqs))/onepq
                    qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*&
                         (qratio_dx(3)/onepqs))/onepq
                    qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*&
                         (qratio_dx(3)/onepqs))/onepq
                    qratiov_dx(1) = (qratiov_dx(1) - qratiov*&
                         (qratio_dx(1)/onepqs))/onepq
                    qratiov_dx(2) = (qratiov_dx(2) - qratiov*&
                         (qratio_dx(2)/onepqs))/onepq
                    qratiov_dx(3) = (qratiov_dx(3) - qratiov*&
                         (qratio_dx(3)/onepqs))/onepq
                    ! transform first order quantities to derivative of
                    ! ln(1 + qratio)
                    qratio_dx(1) = qratio_dx(1)/onepq
                    qratio_dx(2) = qratio_dx(2)/onepq
                    qratio_dx(3) = qratio_dx(3)/onepq
                 endif
                 ! transform second order quantities to derivative of
                 ! ln(1 + qratio)
                 qratiovt = (qratiovt - qratiov*&
                      (qratiot/onepqs))/onepq
                 qratiov_dx(4) = (qratiov_dx(4) - qratiov*&
                      (qratio_dx(4)/onepqs))/onepq
                 ! transform first order quantities to derivative of
                 ! ln(1 + qratio)
                 qratio_dx(4) = qratio_dx(4)/onepq
                 qratiov = qratiov/onepq
              endif
              ! transform first order quantities to derivative of
              ! ln(1 + qratio)
              qratiot = qratiot/onepq
              if(qratio_us.gt.1.e-3_fp_kind) then
                 qratio = qratio_scale*log(onepq)
              else
                 ! alternating series so relative error is less than first
                 ! missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
                 qratio = qratio*&
                      (1._fp_kind      - qratio_us*&
                      (1._fp_kind/2._fp_kind - qratio_us*&
                      (1._fp_kind/3._fp_kind - qratio_us*&
                      (1._fp_kind/4._fp_kind - qratio_us*&
                      (1._fp_kind/5._fp_kind)))))
              endif
              if(ifmhd_logical) then
                 ! calculate derivative *correction* due to x dependence on fl
                 ! and tl.
                 if(ifnr03) then
                    dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                    dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                    dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                    dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                    dqratio_dxf(4) = qratio_dx2(4,4)*extrasumf(nextrasum-1)
                    dqratio_dxt(4) = qratio_dx2(4,4)*extrasumt(nextrasum-1)
                    dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                    dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
                 endif
                 if(ifnr13) then
                    dqratio_dv(1:max_index) = qratio_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                    dqratiov_dv(1:max_index) = qratiov_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,4) = qratio_dx2(4,4)*extrasum_dv(1:max_index,nextrasum-1)
                 endif
                 if(iz.eq.1) then
                    if(ifnr03) then
                       dqratiof = dqratiof + qratio_dx(1)*extrasumf(1) + qratio_dx(2)*extrasumf(2) + qratio_dx(3)*extrasumf(3)
                       dqratiot = dqratiot + qratio_dx(1)*extrasumt(1) + qratio_dx(2)*extrasumt(2) + qratio_dx(3)*extrasumt(3)
                       dqratiotf = dqratiotf + qratiot_dx(1)*extrasumf(1) + qratiot_dx(2)*extrasumf(2) + qratiot_dx(3)*extrasumf(3)
                       dqratiott = dqratiott + qratiot_dx(1)*extrasumt(1) + qratiot_dx(2)*extrasumt(2) + qratiot_dx(3)*extrasumt(3)
                       dqratio_dxf(1) =&
                            qratio_dx2(1,1)*extrasumf(1) +&
                            qratio_dx2(2,1)*extrasumf(2) +&
                            qratio_dx2(3,1)*extrasumf(3) +&
                            qratio_dx2(4,1)*extrasumf(nextrasum-1)
                       dqratio_dxt(1) =&
                            qratio_dx2(1,1)*extrasumt(1) +&
                            qratio_dx2(2,1)*extrasumt(2) +&
                            qratio_dx2(3,1)*extrasumt(3) +&
                            qratio_dx2(4,1)*extrasumt(nextrasum-1)
                       dqratio_dxf(2) =&
                            qratio_dx2(2,1)*extrasumf(1) +&
                            qratio_dx2(2,2)*extrasumf(2) +&
                            qratio_dx2(3,2)*extrasumf(3) +&
                            qratio_dx2(4,2)*extrasumf(nextrasum-1)
                       dqratio_dxt(2) =&
                            qratio_dx2(2,1)*extrasumt(1) +&
                            qratio_dx2(2,2)*extrasumt(2) +&
                            qratio_dx2(3,2)*extrasumt(3) +&
                            qratio_dx2(4,2)*extrasumt(nextrasum-1)
                       dqratio_dxf(3) =&
                            qratio_dx2(3,1)*extrasumf(1) +&
                            qratio_dx2(3,2)*extrasumf(2) +&
                            qratio_dx2(3,3)*extrasumf(3) +&
                            qratio_dx2(4,3)*extrasumf(nextrasum-1)
                       dqratio_dxt(3) =&
                            qratio_dx2(3,1)*extrasumt(1) +&
                            qratio_dx2(3,2)*extrasumt(2) +&
                            qratio_dx2(3,3)*extrasumt(3) +&
                            qratio_dx2(4,3)*extrasumt(nextrasum-1)
                       dqratio_dxf(4) = dqratio_dxf(4) +&
                            qratio_dx2(4,1)*extrasumf(1) +&
                            qratio_dx2(4,2)*extrasumf(2) +&
                            qratio_dx2(4,3)*extrasumf(3)
                       dqratio_dxt(4) = dqratio_dxt(4) +&
                            qratio_dx2(4,1)*extrasumt(1) +&
                            qratio_dx2(4,2)*extrasumt(2) +&
                            qratio_dx2(4,3)*extrasumt(3)
                       dqratiovf = dqratiovf +&
                            qratiov_dx(1)*extrasumf(1) +&
                            qratiov_dx(2)*extrasumf(2) +&
                            qratiov_dx(3)*extrasumf(3)
                       dqratiovt = dqratiovt +&
                            qratiov_dx(1)*extrasumt(1) +&
                            qratiov_dx(2)*extrasumt(2) +&
                            qratiov_dx(3)*extrasumt(3)
                    endif
                    if(ifnr13) then
                       dqratio_dv(1:max_index) = dqratio_dv(1:max_index) +&
                            qratio_dx(1)*extrasum_dv(1:max_index,1) +&
                            qratio_dx(2)*extrasum_dv(1:max_index,2) +&
                            qratio_dx(3)*extrasum_dv(1:max_index,3)
                       dqratiov_dv(1:max_index) = dqratiov_dv(1:max_index) +&
                            qratiov_dx(1)*extrasum_dv(1:max_index,1) +&
                            qratiov_dx(2)*extrasum_dv(1:max_index,2) +&
                            qratiov_dx(3)*extrasum_dv(1:max_index,3)
                       dqratio_dxdv(1:max_index,1) =&
                            qratio_dx2(1,1)*extrasum_dv(1:max_index,1) +&
                            qratio_dx2(2,1)*extrasum_dv(1:max_index,2) +&
                            qratio_dx2(3,1)*extrasum_dv(1:max_index,3) +&
                            qratio_dx2(4,1)*extrasum_dv(1:max_index,nextrasum-1)
                       dqratio_dxdv(1:max_index,2) =&
                            qratio_dx2(2,1)*extrasum_dv(1:max_index,1) +&
                            qratio_dx2(2,2)*extrasum_dv(1:max_index,2) +&
                            qratio_dx2(3,2)*extrasum_dv(1:max_index,3) +&
                            qratio_dx2(4,2)*extrasum_dv(1:max_index,nextrasum-1)
                       dqratio_dxdv(1:max_index,3) =&
                            qratio_dx2(3,1)*extrasum_dv(1:max_index,1) +&
                            qratio_dx2(3,2)*extrasum_dv(1:max_index,2) +&
                            qratio_dx2(3,3)*extrasum_dv(1:max_index,3) +&
                            qratio_dx2(4,3)*extrasum_dv(1:max_index,nextrasum-1)
                       dqratio_dxdv(1:max_index,4) = dqratio_dxdv(1:max_index,4) +&
                            qratio_dx2(4,1)*extrasum_dv(1:max_index,1) +&
                            qratio_dx2(4,2)*extrasum_dv(1:max_index,2) +&
                            qratio_dx2(4,3)*extrasum_dv(1:max_index,3)
                    endif ! if(ifnr03) then.
                 endif ! if(iz.eq.1) then ...
              endif ! if(ifmhd_logical) then.....
              extrace_tag=[ielement,iz-1]
              call exsum_component_add(&
                   ifmhd_logical, iz.eq.1, ifnr03, ifnr13,.false.,&
                   ion_start, ion_end(ielement),&
                   1, max_index,&
                   inv_ion,&
                   nuvar(iz,index), nuvarf(iz,index), nuvart(iz,index), nuvar_dv(:, iz,index),&
                   xextrasum_scale, qratio_scale,&
                   qratio, qratiot, qratio_dx,&
                   qratiot2, qratiot_dx,&
                   qratio_dx2,&
                   qratiov, qratiovt, qratiov_dx,&
                   dqratiof, dqratiotf, dqratio_dxf,&
                   dqratiot, dqratiott, dqratio_dxt,&
                   dqratio_dv,&
                   dqratio_dxdv,&
                   dqratiovf, dqratiovt, dqratiov_dv,&
                   psum, psumf, psumt, psum_dv,&
                   ssum, ssumf, ssumt,&
                   usum,&
                   free_sum, free_sumf, free_sum_dv,&
                   xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
           endif ! if(qratio.gt.exp_lim) then
        enddo ! do ion = ion_start, i... neutral to next to bare.
     endif ! if((ielement.gt
  enddo ! do index = 1,... for atoms only
  if(partial_elements(1).eq.1) then
     ! do all the molecular stuff below only if hydrogen has a non-zero
     ! abundance
     ! now do H2 equilibrium constant change due to neutral monatomic H.
     ! now do neutral monatomic H
     ! (which skipped previously if molecular formation).
     if(ifh2.gt.0) then
        ielement = 1
        iz = 1
        ion = 1
        nmin_s = nmin_species(ion)
        if(nmin_s.lt.nmin(iz).or.&
             nmin_s.gt.nmin_max(iz))&
             error stop 'excitation_sum: invalid nmin or nmin_max'
        exparg = exparg_shift -&
             c2t*bion(ion) - plop(ion)
        ! correct for MHD ground state occupation probability.
        if(ifmhd_logical) then
           if(iz.eq.1) then
              ln_ground_occ = occ_const*(&
                   x(1) + r_neutral(ielement)*(&
                   x(2)*3._fp_kind + r_neutral(ielement)*(&
                   x(3)*3._fp_kind + r_neutral(ielement)*(&
                   x(4)))))
           else
              ln_ground_occ = 0._fp_kind
           endif
           ln_ground_occ = ln_ground_occ +&
                occ_const*r_ion3(ion)*x(5)
           !            ignore excitation correction if ground state wiped out by
           !            pressure ionization in any case.
           if(ln_ground_occ.gt.ground_state_lim) then
              exparg = exparg - ln_ground_occ
           else
              !              signal to ignore this species
              exparg = -1000._fp_kind!exparg_lim - 1._fp_kind
           endif
        endif
        ! if(exparg.gt.exparg_lim) then
        ! Keep the Boltzmann factor logarithmic until the partition product.
        ! else
        !   exparg = 0.d0
        ! endif
        ! qratio is first ratio of excited to ground partition function,
        ! but then is transformed to ln(1 + qratio).  Similarly, there
        ! is a subsequent transformation of all derivatives.
        ! N.B. important convention on partial derivative variable names:
        ! names starting with "qratio" are partial derivatives assuming
        ! that qratio is a function of tl and x.
        ! names starting with "dqratio" are the *change* to the partial
        ! derivative caused by x being a function of fl, tl, dv.
        if(iz.eq.1) then
           ! neutral monatomic H as part of H2 calculation
           qratio = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))*&
                real((iqion(ion)),fp_kind)/real((iqneutral(ielement)),fp_kind)
        else
           ! unused since iz set to 1 above.
           qratio = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))*&
                real((iqion(ion)),fp_kind)/real((iqion(ion-1)),fp_kind)
        endif
        if(qratio.gt.exp_lim) then
           qratio_us = qratio/qratio_scale
           onepq = 1._fp_kind + qratio_us
           onepqs = onepq*qratio_scale
           ! qratiot and qratiot2 are the first and second partial wrt tl
           qratiot = qratio*(c2t*bion(ion) - plopt(ion) +&
                qstart(nmin_s,iz)/qstar(nmin_s,iz))
           qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +&
                qstart(nmin_s,iz)/qstar(nmin_s,iz)) +&
                qratio*(-c2t*bion(ion) - plopt2(ion) +&
                qstart2(nmin_s,iz)/qstar(nmin_s,iz) -&
                (qstart(nmin_s,iz)/qstar(nmin_s,iz))*&
                (qstart(nmin_s,iz)/qstar(nmin_s,iz)))
           ! transform second order quantities to derivative of
           ! ln(1 + qratio)
           qratiot2 = (qratiot2 - qratiot*&
                (qratiot/onepqs))/onepq
           if(ifmhd_logical) then
              ! qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are
              ! partials wrt tl, extrasum(k) and extrasum(l) following
              ! special convention for k or l = 4.
              qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
              qratiot_dx(4) = qratiot*&
                   (-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +&
                   qratio*(qstartx(5,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
              qratio_dx2(4,4) = qratio_dx(4)*&
                   (-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +&
                   qratio*(qstarx2(5,5,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)))
              ! derivative of qratio wrt ln V
              qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
              qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
              qratiov_dx(4) =&
                   -qratio_dx2(4,4)*extrasum(nextrasum-1) -&
                   qratio_dx(4)
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiot_dx(4) = (qratiot_dx(4) - qratiot*&
                   (qratio_dx(4)/onepqs))/onepq
              qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*&
                   (qratio_dx(4)/onepqs))/onepq
              if(iz.eq.1) then
                 ! N.B. x(1) (or extrasum(4)) dependence divides out of qratio
                 ! n.b. indices are reordered here so that
                 ! qratio_dx(k) refers to derivative wrt extrasum(k) for
                 ! k = 1, 2, 3, and qratio_dx(4) refers to the derivative
                 ! wrt extrasum(nextrasum-1).
                 qratio_dx(1) = qratio*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                 qratio_dx(2) = qratio*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                 qratio_dx(3) = qratio*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                 qratiot_dx(1) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstartx(4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                 qratiot_dx(2) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstartx(3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                 qratiot_dx(3) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstartx(2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(1,1) = qratio_dx(1)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(4,4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(2,1) = qratio_dx(2)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(4,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(3,1) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(4,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(4,1) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(5,4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(2,2) = qratio_dx(2)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(3,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(3,2) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(3,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(4,2) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(5,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(3,3) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(2,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(4,3) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(5,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 ! derivative of qratio wrt ln V
                 qratiov = qratiov -&
                      qratio_dx(1)*extrasum(1) -&
                      qratio_dx(2)*extrasum(2) -&
                      qratio_dx(3)*extrasum(3)
                 qratiovt = qratiovt -&
                      qratiot_dx(1)*extrasum(1) -&
                      qratiot_dx(2)*extrasum(2) -&
                      qratiot_dx(3)*extrasum(3)
                 qratiov_dx(1) = -&
                      qratio_dx2(1,1)*extrasum(1) -&
                      qratio_dx2(2,1)*extrasum(2) -&
                      qratio_dx2(3,1)*extrasum(3) -&
                      qratio_dx2(4,1)*extrasum(nextrasum-1) -&
                      qratio_dx(1)
                 qratiov_dx(2) = -&
                      qratio_dx2(2,1)*extrasum(1) -&
                      qratio_dx2(2,2)*extrasum(2) -&
                      qratio_dx2(3,2)*extrasum(3) -&
                      qratio_dx2(4,2)*extrasum(nextrasum-1) -&
                      qratio_dx(2)
                 qratiov_dx(3) = -&
                      qratio_dx2(3,1)*extrasum(1) -&
                      qratio_dx2(3,2)*extrasum(2) -&
                      qratio_dx2(3,3)*extrasum(3) -&
                      qratio_dx2(4,3)*extrasum(nextrasum-1) -&
                      qratio_dx(3)
                 qratiov_dx(4) = qratiov_dx(4) -&
                      qratio_dx2(4,1)*extrasum(1) -&
                      qratio_dx2(4,2)*extrasum(2) -&
                      qratio_dx2(4,3)*extrasum(3)
                 ! transform second order quantities to derivative of
                 ! ln(1 + qratio)
                 qratiot_dx(1) = (qratiot_dx(1) - qratiot*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratiot_dx(2) = (qratiot_dx(2) - qratiot*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratiot_dx(3) = (qratiot_dx(3) - qratiot*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratiov_dx(1) = (qratiov_dx(1) - qratiov*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratiov_dx(2) = (qratiov_dx(2) - qratiov*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratiov_dx(3) = (qratiov_dx(3) - qratiov*&
                      (qratio_dx(3)/onepqs))/onepq
                 ! transform first order quantities to derivative of
                 ! ln(1 + qratio)
                 qratio_dx(1) = qratio_dx(1)/onepq
                 qratio_dx(2) = qratio_dx(2)/onepq
                 qratio_dx(3) = qratio_dx(3)/onepq
              endif ! if(iz.eq.1) then....
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiovt = (qratiovt - qratiov*&
                   (qratiot/onepqs))/onepq
              qratiov_dx(4) = (qratiov_dx(4) - qratiov*&
                   (qratio_dx(4)/onepqs))/onepq
              ! transform first order quantities to derivative of
              ! ln(1 + qratio)
              qratio_dx(4) = qratio_dx(4)/onepq
              qratiov = qratiov/onepq
           endif ! if(ifmhd_logical) then....
           ! transform first order quantities to derivative of
           ! ln(1 + qratio)
           qratiot = qratiot/onepq
           if(qratio_us.gt.1.e-3_fp_kind) then
              qratio = qratio_scale*log(onepq)
           else
              ! alternating series so relative error is less than first
              ! missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              qratio = qratio*&
                   (1._fp_kind      - qratio_us*&
                   (1._fp_kind/2._fp_kind - qratio_us*&
                   (1._fp_kind/3._fp_kind - qratio_us*&
                   (1._fp_kind/4._fp_kind - qratio_us*&
                   (1._fp_kind/5._fp_kind)))))
           endif
           if(ifmhd_logical) then
              ! calculate derivative *correction* due to x dependence on fl
              ! and tl.
              if(ifnr03) then
                 dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                 dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                 dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                 dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                 dqratio_dxf(4) = qratio_dx2(4,4)*extrasumf(nextrasum-1)
                 dqratio_dxt(4) = qratio_dx2(4,4)*extrasumt(nextrasum-1)
                 dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                 dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
              endif
              if(ifnr13) then
                 dqratio_dv(1:max_index) = qratio_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                 dqratiov_dv(1:max_index) = qratiov_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                 dqratio_dxdv(1:max_index,4) = qratio_dx2(4,4)*extrasum_dv(1:max_index,nextrasum-1)
              endif
              if(iz.eq.1) then
                 if(ifnr03) then
                    dqratiof = dqratiof + qratio_dx(1)*extrasumf(1) + qratio_dx(2)*extrasumf(2) + qratio_dx(3)*extrasumf(3)
                    dqratiot = dqratiot + qratio_dx(1)*extrasumt(1) + qratio_dx(2)*extrasumt(2) + qratio_dx(3)*extrasumt(3)
                    dqratiotf = dqratiotf + qratiot_dx(1)*extrasumf(1) + qratiot_dx(2)*extrasumf(2) + qratiot_dx(3)*extrasumf(3)
                    dqratiott = dqratiott + qratiot_dx(1)*extrasumt(1) + qratiot_dx(2)*extrasumt(2) + qratiot_dx(3)*extrasumt(3)
                    dqratio_dxf(1) =&
                         qratio_dx2(1,1)*extrasumf(1) +&
                         qratio_dx2(2,1)*extrasumf(2) +&
                         qratio_dx2(3,1)*extrasumf(3) +&
                         qratio_dx2(4,1)*extrasumf(nextrasum-1)
                    dqratio_dxt(1) =&
                         qratio_dx2(1,1)*extrasumt(1) +&
                         qratio_dx2(2,1)*extrasumt(2) +&
                         qratio_dx2(3,1)*extrasumt(3) +&
                         qratio_dx2(4,1)*extrasumt(nextrasum-1)
                    dqratio_dxf(2) =&
                         qratio_dx2(2,1)*extrasumf(1) +&
                         qratio_dx2(2,2)*extrasumf(2) +&
                         qratio_dx2(3,2)*extrasumf(3) +&
                         qratio_dx2(4,2)*extrasumf(nextrasum-1)
                    dqratio_dxt(2) =&
                         qratio_dx2(2,1)*extrasumt(1) +&
                         qratio_dx2(2,2)*extrasumt(2) +&
                         qratio_dx2(3,2)*extrasumt(3) +&
                         qratio_dx2(4,2)*extrasumt(nextrasum-1)
                    dqratio_dxf(3) =&
                         qratio_dx2(3,1)*extrasumf(1) +&
                         qratio_dx2(3,2)*extrasumf(2) +&
                         qratio_dx2(3,3)*extrasumf(3) +&
                         qratio_dx2(4,3)*extrasumf(nextrasum-1)
                    dqratio_dxt(3) =&
                         qratio_dx2(3,1)*extrasumt(1) +&
                         qratio_dx2(3,2)*extrasumt(2) +&
                         qratio_dx2(3,3)*extrasumt(3) +&
                         qratio_dx2(4,3)*extrasumt(nextrasum-1)
                    dqratio_dxf(4) = dqratio_dxf(4) +&
                         qratio_dx2(4,1)*extrasumf(1) +&
                         qratio_dx2(4,2)*extrasumf(2) +&
                         qratio_dx2(4,3)*extrasumf(3)
                    dqratio_dxt(4) = dqratio_dxt(4) +&
                         qratio_dx2(4,1)*extrasumt(1) +&
                         qratio_dx2(4,2)*extrasumt(2) +&
                         qratio_dx2(4,3)*extrasumt(3)
                    dqratiovf = dqratiovf +&
                         qratiov_dx(1)*extrasumf(1) +&
                         qratiov_dx(2)*extrasumf(2) +&
                         qratiov_dx(3)*extrasumf(3)
                    dqratiovt = dqratiovt +&
                         qratiov_dx(1)*extrasumt(1) +&
                         qratiov_dx(2)*extrasumt(2) +&
                         qratiov_dx(3)*extrasumt(3)
                 endif
                 if(ifnr13) then
                    dqratio_dv(1:max_index) =&
                         dqratio_dv(1:max_index) +&
                         qratio_dx(1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx(2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx(3)*extrasum_dv(1:max_index,3)
                    dqratiov_dv(1:max_index) =&
                         dqratiov_dv(1:max_index) +&
                         qratiov_dx(1)*extrasum_dv(1:max_index,1) +&
                         qratiov_dx(2)*extrasum_dv(1:max_index,2) +&
                         qratiov_dx(3)*extrasum_dv(1:max_index,3)
                    dqratio_dxdv(1:max_index,1) =&
                         qratio_dx2(1,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(2,1)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,1)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,1)*&
                         extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,2) =&
                         qratio_dx2(2,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(2,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,2)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,2)*&
                         extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,3) =&
                         qratio_dx2(3,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(3,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,3)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,3)*&
                         extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,4) =&
                         dqratio_dxdv(1:max_index,4) +&
                         qratio_dx2(4,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(4,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(4,3)*extrasum_dv(1:max_index,3)
                 endif ! if(ifnr03) then...
              endif ! if(iz.eq.1) then...
           endif ! if(ifmhd_logical)...
           ! monatomic H component of excitation sums.
           extrace_tag=[1,0]
              call exsum_component_add(&
                ifmhd_logical, iz.eq.1, ifnr03, ifnr13,.true.,&
                1, max_index,&
                1, max_index,&
                inv_ion,&
                nug, nugf, nugt, nug_dv,&
                xextrasum_scale, qratio_scale,&
                qratio, qratiot, qratio_dx,&
                qratiot2, qratiot_dx,&
                qratio_dx2,&
                qratiov, qratiovt, qratiov_dx,&
                dqratiof, dqratiotf, dqratio_dxf,&
                dqratiot, dqratiott, dqratio_dxt,&
                dqratio_dv,&
                dqratio_dxdv,&
                dqratiovf, dqratiovt, dqratiov_dv,&
                psum, psumf, psumt, psum_dv,&
                ssum, ssumf, ssumt,&
                usum,&
                free_sum, free_sumf, free_sum_dv,&
                xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
        endif ! if(qratio.gt.exp_lim) then
     endif ! if(ifh2.gt.0) then... neutral monatomic H change to H2
     ! now do H2 change to H2 equilibrium constant
     if(ifh2.gt.0.and.ifh2plus.gt.0.and.&
          mod(ifexcited,10).gt.1) then
        ielement = nelements + 1
        iz = 1
        ion = nions + 1
        nmin_s = nmin_species(ion)
        if(nmin_s.lt.nmin(iz).or.&
             nmin_s.gt.nmin_max(iz))&
             error stop 'excitation_sum: invalid nmin or nmin_max'
        exparg = exparg_shift -&
             c2t*bion(ion) - plop(ion) + qh2plus - qh2
        !          correct for MHD ground state occupation probability.
        if(ifmhd_logical) then
           if(iz.eq.1) then
              ln_ground_occ = occ_const*(&
                   x(1) + r_neutral(ielement)*(&
                   x(2)*3._fp_kind + r_neutral(ielement)*(&
                   x(3)*3._fp_kind + r_neutral(ielement)*(&
                   x(4)))))
           else
              ln_ground_occ = 0._fp_kind
           endif
           ln_ground_occ = ln_ground_occ +&
                occ_const*r_ion3(ion)*x(5)
           ! ignore excitation correction if ground state wiped out by
           ! pressure ionization in any case.
           if(ln_ground_occ.gt.ground_state_lim) then
              exparg = exparg - ln_ground_occ
           else
              ! signal to ignore this species
              exparg = -1000._fp_kind!exparg_lim - 1._fp_kind
           endif
        endif
        ! if(exparg.gt.exparg_lim) then
        ! Keep the Boltzmann factor logarithmic until the partition product.
        ! else
        !   exparg = 0.d0
        ! endif
        ! qratio is first ratio of excited to ground partition function,
        ! but then is transformed to ln(1 + qratio).  Similarly, there
        ! is a subsequent transformation of all derivatives.
        ! N.B. important convention on partial derivative variable names:
        ! names starting with "qratio" are partial derivatives assuming
        ! that qratio is a function of tl and x.
        ! names starting with "dqratio" are the *change* to the partial
        ! derivative caused by x being a function of fl, tl, dv.
        ! ratio of excited to ground state partition functions
        ! h2
        qratio = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))
        if(qratio.gt.exp_lim) then
           qratio_us = qratio/qratio_scale
           onepq = 1._fp_kind + qratio_us
           onepqs = onepq*qratio_scale
           ! qratiot and qratiot2 are the first and second partial wrt tl
           qratiot = qratio*(c2t*bion(ion) - plopt(ion) +&
                qh2plust - qh2t +&
                qstart(nmin_s,iz)/qstar(nmin_s,iz))
           qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) +&
                qh2plust - qh2t +&
                qstart(nmin_s,iz)/qstar(nmin_s,iz)) +&
                qratio*(-c2t*bion(ion) - plopt2(ion) +&
                qh2plust2 - qh2t2 +&
                qstart2(nmin_s,iz)/qstar(nmin_s,iz) -&
                (qstart(nmin_s,iz)/qstar(nmin_s,iz))*&
                (qstart(nmin_s,iz)/qstar(nmin_s,iz)))
           ! transform second order quantities to derivative of
           ! ln(1 + qratio)
           qratiot2 = (qratiot2 - qratiot*&
                (qratiot/onepqs))/onepq
           if(ifmhd_logical) then
              ! qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are
              ! partials wrt tl, extrasum(k) and extrasum(l) following
              ! special convention for k or l = 4.
              qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
              qratiot_dx(4) = qratiot*&
                   (-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +&
                   qratio*(qstartx(5,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
              qratio_dx2(4,4) = qratio_dx(4)*&
                   (-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)) +&
                   qratio*(qstarx2(5,5,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz)))
              ! derivative of qratio wrt ln V
              qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
              qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
              qratiov_dx(4) =&
                   -qratio_dx2(4,4)*extrasum(nextrasum-1) -&
                   qratio_dx(4)
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiot_dx(4) = (qratiot_dx(4) - qratiot*&
                   (qratio_dx(4)/onepqs))/onepq
              qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*&
                   (qratio_dx(4)/onepqs))/onepq
              if(iz.eq.1) then
                 ! N.B. x(1) (or extrasum(4)) dependence divides out of qratio
                 ! n.b. indices are reordered here so that
                 ! qratio_dx(k) refers to derivative wrt extrasum(k) for
                 ! k = 1, 2, 3, and qratio_dx(4) refers to the derivative
                 ! wrt extrasum(nextrasum-1).
                 qratio_dx(1) = qratio*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                 qratio_dx(2) = qratio*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                 qratio_dx(3) = qratio*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                 qratiot_dx(1) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstartx(4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                 qratiot_dx(2) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstartx(3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                 qratiot_dx(3) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstartx(2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))*(qstart(nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(1,1) = qratio_dx(1)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(4,4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(2,1) = qratio_dx(2)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(4,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(3,1) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(4,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(4,1) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(5,4,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(2,2) = qratio_dx(2)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(3,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(3,2) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(3,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(4,2) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(5,3,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(3,3) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(2,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 qratio_dx2(4,3) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)) +&
                      qratio*(qstarx2(5,2,nmin_s,iz)/qstar(nmin_s,iz)-(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz)))
                 ! derivative of qratio wrt ln V
                 qratiov = qratiov -&
                      qratio_dx(1)*extrasum(1) -&
                      qratio_dx(2)*extrasum(2) -&
                      qratio_dx(3)*extrasum(3)
                 qratiovt = qratiovt -&
                      qratiot_dx(1)*extrasum(1) -&
                      qratiot_dx(2)*extrasum(2) -&
                      qratiot_dx(3)*extrasum(3)
                 qratiov_dx(1) = -&
                      qratio_dx2(1,1)*extrasum(1) -&
                      qratio_dx2(2,1)*extrasum(2) -&
                      qratio_dx2(3,1)*extrasum(3) -&
                      qratio_dx2(4,1)*extrasum(nextrasum-1) -&
                      qratio_dx(1)
                 qratiov_dx(2) = -&
                      qratio_dx2(2,1)*extrasum(1) -&
                      qratio_dx2(2,2)*extrasum(2) -&
                      qratio_dx2(3,2)*extrasum(3) -&
                      qratio_dx2(4,2)*extrasum(nextrasum-1) -&
                      qratio_dx(2)
                 qratiov_dx(3) = -&
                      qratio_dx2(3,1)*extrasum(1) -&
                      qratio_dx2(3,2)*extrasum(2) -&
                      qratio_dx2(3,3)*extrasum(3) -&
                      qratio_dx2(4,3)*extrasum(nextrasum-1) -&
                      qratio_dx(3)
                 qratiov_dx(4) = qratiov_dx(4) -&
                      qratio_dx2(4,1)*extrasum(1) -&
                      qratio_dx2(4,2)*extrasum(2) -&
                      qratio_dx2(4,3)*extrasum(3)
                 ! transform second order quantities to derivative of
                 ! ln(1 + qratio)
                 qratiot_dx(1) = (qratiot_dx(1) - qratiot*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratiot_dx(2) = (qratiot_dx(2) - qratiot*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratiot_dx(3) = (qratiot_dx(3) - qratiot*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratiov_dx(1) = (qratiov_dx(1) - qratiov*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratiov_dx(2) = (qratiov_dx(2) - qratiov*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratiov_dx(3) = (qratiov_dx(3) - qratiov*&
                      (qratio_dx(3)/onepqs))/onepq
                 ! transform first order quantities to derivative of
                 ! ln(1 + qratio)
                 qratio_dx(1) = qratio_dx(1)/onepq
                 qratio_dx(2) = qratio_dx(2)/onepq
                 qratio_dx(3) = qratio_dx(3)/onepq
              endif
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiovt = (qratiovt - qratiov*&
                   (qratiot/onepqs))/onepq
              qratiov_dx(4) = (qratiov_dx(4) - qratiov*&
                   (qratio_dx(4)/onepqs))/onepq
              ! transform first order quantities to derivative of
              ! ln(1 + qratio)
              qratio_dx(4) = qratio_dx(4)/onepq
              qratiov = qratiov/onepq
           endif
           ! transform first order quantities to derivative of
           ! ln(1 + qratio)
           qratiot = qratiot/onepq
           if(qratio_us.gt.1.e-3_fp_kind) then
              qratio = qratio_scale*log(onepq)
           else
              ! alternating series so relative error is less than first
              ! missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              qratio = qratio*&
                   (1._fp_kind      - qratio_us*&
                   (1._fp_kind/2._fp_kind - qratio_us*&
                   (1._fp_kind/3._fp_kind - qratio_us*&
                   (1._fp_kind/4._fp_kind - qratio_us*&
                   (1._fp_kind/5._fp_kind)))))
           endif
           if(ifmhd_logical) then
              ! calculate derivative *correction* due to x dependence on fl
              ! and tl.
              if(ifnr03) then
                 dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                 dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                 dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                 dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                 dqratio_dxf(4) = qratio_dx2(4,4)*&
                      extrasumf(nextrasum-1)
                 dqratio_dxt(4) = qratio_dx2(4,4)*&
                      extrasumt(nextrasum-1)
                 dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                 dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
              endif
              if(ifnr13) then
                 dqratio_dv(1:max_index) = qratio_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                 dqratiov_dv(1:max_index) = qratiov_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                 dqratio_dxdv(1:max_index,4) = qratio_dx2(4,4)*extrasum_dv(1:max_index,nextrasum-1)
              endif
              if(iz.eq.1) then
                 if(ifnr03) then
                    dqratiof = dqratiof + qratio_dx(1)*extrasumf(1) + qratio_dx(2)*extrasumf(2) + qratio_dx(3)*extrasumf(3)
                    dqratiot = dqratiot + qratio_dx(1)*extrasumt(1) + qratio_dx(2)*extrasumt(2) + qratio_dx(3)*extrasumt(3)
                    dqratiotf = dqratiotf + qratiot_dx(1)*extrasumf(1) + qratiot_dx(2)*extrasumf(2) + qratiot_dx(3)*extrasumf(3)
                    dqratiott = dqratiott + qratiot_dx(1)*extrasumt(1) + qratiot_dx(2)*extrasumt(2) + qratiot_dx(3)*extrasumt(3)
                    dqratio_dxf(1) =&
                         qratio_dx2(1,1)*extrasumf(1) +&
                         qratio_dx2(2,1)*extrasumf(2) +&
                         qratio_dx2(3,1)*extrasumf(3) +&
                         qratio_dx2(4,1)*extrasumf(nextrasum-1)
                    dqratio_dxt(1) =&
                         qratio_dx2(1,1)*extrasumt(1) +&
                         qratio_dx2(2,1)*extrasumt(2) +&
                         qratio_dx2(3,1)*extrasumt(3) +&
                         qratio_dx2(4,1)*extrasumt(nextrasum-1)
                    dqratio_dxf(2) =&
                         qratio_dx2(2,1)*extrasumf(1) +&
                         qratio_dx2(2,2)*extrasumf(2) +&
                         qratio_dx2(3,2)*extrasumf(3) +&
                         qratio_dx2(4,2)*extrasumf(nextrasum-1)
                    dqratio_dxt(2) =&
                         qratio_dx2(2,1)*extrasumt(1) +&
                         qratio_dx2(2,2)*extrasumt(2) +&
                         qratio_dx2(3,2)*extrasumt(3) +&
                         qratio_dx2(4,2)*extrasumt(nextrasum-1)
                    dqratio_dxf(3) =&
                         qratio_dx2(3,1)*extrasumf(1) +&
                         qratio_dx2(3,2)*extrasumf(2) +&
                         qratio_dx2(3,3)*extrasumf(3) +&
                         qratio_dx2(4,3)*extrasumf(nextrasum-1)
                    dqratio_dxt(3) =&
                         qratio_dx2(3,1)*extrasumt(1) +&
                         qratio_dx2(3,2)*extrasumt(2) +&
                         qratio_dx2(3,3)*extrasumt(3) +&
                         qratio_dx2(4,3)*extrasumt(nextrasum-1)
                    dqratio_dxf(4) = dqratio_dxf(4) +&
                         qratio_dx2(4,1)*extrasumf(1) +&
                         qratio_dx2(4,2)*extrasumf(2) +&
                         qratio_dx2(4,3)*extrasumf(3)
                    dqratio_dxt(4) = dqratio_dxt(4) +&
                         qratio_dx2(4,1)*extrasumt(1) +&
                         qratio_dx2(4,2)*extrasumt(2) +&
                         qratio_dx2(4,3)*extrasumt(3)
                    dqratiovf = dqratiovf +&
                         qratiov_dx(1)*extrasumf(1) +&
                         qratiov_dx(2)*extrasumf(2) +&
                         qratiov_dx(3)*extrasumf(3)
                    dqratiovt = dqratiovt +&
                         qratiov_dx(1)*extrasumt(1) +&
                         qratiov_dx(2)*extrasumt(2) +&
                         qratiov_dx(3)*extrasumt(3)
                 endif
                 if(ifnr13) then
                    dqratio_dv(1:max_index) =&
                         dqratio_dv(1:max_index) +&
                         qratio_dx(1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx(2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx(3)*extrasum_dv(1:max_index,3)
                    dqratiov_dv(1:max_index) =&
                         dqratiov_dv(1:max_index) +&
                         qratiov_dx(1)*extrasum_dv(1:max_index,1) +&
                         qratiov_dx(2)*extrasum_dv(1:max_index,2) +&
                         qratiov_dx(3)*extrasum_dv(1:max_index,3)
                    dqratio_dxdv(1:max_index,1) =&
                         qratio_dx2(1,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(2,1)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,1)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,1)*&
                         extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,2) =&
                         qratio_dx2(2,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(2,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,2)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,2)*&
                         extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,3) =&
                         qratio_dx2(3,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(3,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,3)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,3)*&
                         extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,4) =&
                         dqratio_dxdv(1:max_index,4) +&
                         qratio_dx2(4,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(4,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(4,3)*extrasum_dv(1:max_index,3)
                 endif
              endif
           endif
           ! H2 component of excitation sums.
           extrace_tag=[25,0]
              call exsum_component_add(&
                ifmhd_logical, iz.eq.1, ifnr03, ifnr13,.true.,&
                1, max_index,&
                1, max_index,&
                inv_ion,&
                nuh2, nuh2f, nuh2t, nuh2_dv,&
                xextrasum_scale, qratio_scale,&
                qratio, qratiot, qratio_dx,&
                qratiot2, qratiot_dx,&
                qratio_dx2,&
                qratiov, qratiovt, qratiov_dx,&
                dqratiof, dqratiotf, dqratio_dxf,&
                dqratiot, dqratiott, dqratio_dxt,&
                dqratio_dv,&
                dqratio_dxdv,&
                dqratiovf, dqratiovt, dqratiov_dv,&
                psum, psumf, psumt, psum_dv,&
                ssum, ssumf, ssumt,&
                usum,&
                free_sum, free_sumf, free_sum_dv,&
                xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
        endif !if(qratio.gt.exp_lim) then
        ! finished H2 change to H2 equilibrium constant.
     endif
     ! now do H2+
     if(ifh2.gt.0.and.ifh2plus.gt.0.and.&
          mod(ifexcited,10).gt.1) then
        ielement = nelements + 2
        iz = 2
        izqstar = izhi + 1
        ion = nions + 2
        nmin_s = nmin_species(ion)
        if(nmin_s.lt.nmin(iz).or.&
             nmin_s.gt.nmin_max(iz))&
             error stop 'excitation_sum: invalid nmin or nmin_max'
        ! n.b. excited state core is two protons with unity statistical weight
        ! (following usual convention that nuclear spin statistical weights
        ! are divided out).
        exparg = exparg_shift -&
             c2t*bion(ion) - plop(ion) - qh2plus
        ! correct for MHD ground state occupation probability.
        if(ifmhd_logical) then
           ! H2+ part of "neutral" list in this case.
           if(iz.eq.2) then
              ln_ground_occ = occ_const*(&
                   x(1) + r_neutral(ielement)*(&
                   x(2)*3._fp_kind + r_neutral(ielement)*(&
                   x(3)*3._fp_kind + r_neutral(ielement)*(&
                   x(4)))))
           else
              ln_ground_occ = 0._fp_kind
           endif
           ln_ground_occ = ln_ground_occ +&
                occ_const*r_ion3(ion)*x(5)
           ! ignore excitation correction if ground state wiped out by
           ! pressure ionization in any case.
           if(ln_ground_occ.gt.ground_state_lim) then
              exparg = exparg - ln_ground_occ
           else
              ! signal to ignore this species
              exparg = -1000._fp_kind!exparg_lim - 1._fp_kind
           endif
        endif
        ! if(exparg.gt.exparg_lim) then
        ! Keep the Boltzmann factor logarithmic until the partition product.
        ! else
        !   exparg = 0.d0
        ! endif
        ! qratio is first ratio of excited to ground partition function,
        ! but then is transformed to ln(1 + qratio).  Similarly, there
        ! is a subsequent transformation of all derivatives.
        ! N.B. important convention on partial derivative variable names:
        ! names starting with "qratio" are partial derivatives assuming
        ! that qratio is a function of tl and x.
        ! names starting with "dqratio" are the *change* to the partial
        ! derivative caused by x being a function of fl, tl, dv.
        ! ratio of excited to ground state partition functions
        ! h2+
        qratio = partition_product(exparg,qstar(nmin_s,izqstar),qstar_logscale(nmin_s,izqstar))
        if(qratio.gt.exp_lim) then
           qratio_us = qratio/qratio_scale
           onepq = 1._fp_kind + qratio_us
           onepqs = onepq*qratio_scale
           ! qratiot and qratiot2 are the first and second partial wrt tl
           qratiot = qratio*(c2t*bion(ion) - plopt(ion) -&
                qh2plust +&
                qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar))
           qratiot2 = qratiot*(c2t*bion(ion) - plopt(ion) -&
                qh2plust +&
                qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                qratio*(-c2t*bion(ion) - plopt2(ion) -&
                qh2plust2 +&
                qstart2(nmin_s,izqstar)/qstar(nmin_s,izqstar) -&
                (qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar))*&
                (qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)))
           ! transform second order quantities to derivative of
           ! ln(1 + qratio)
           qratiot2 = (qratiot2 - qratiot*&
                (qratiot/onepqs))/onepq
           if(ifmhd_logical) then
              ! qratio_dx(k), qratiot_dx(k), qratio_dx2(k,l) are
              ! partials wrt tl, extrasum(k) and extrasum(l) following
              ! special convention for k or l = 4.
              qratio_dx(4) = qratio*(-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))
              qratiot_dx(4) = qratiot*&
                   (-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                   qratio*(qstartx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)))
              qratio_dx2(4,4) = qratio_dx(4)*&
                   (-occ_const*r_ion3(ion) +&
                   qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                   qratio*(qstarx2(5,5,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
              ! derivative of qratio wrt ln V
              qratiov = -qratio_dx(4)*extrasum(nextrasum-1)
              qratiovt = -qratiot_dx(4)*extrasum(nextrasum-1)
              qratiov_dx(4) =&
                   -qratio_dx2(4,4)*extrasum(nextrasum-1) -&
                   qratio_dx(4)
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiot_dx(4) = (qratiot_dx(4) - qratiot*&
                   (qratio_dx(4)/onepqs))/onepq
              qratio_dx2(4,4) = (qratio_dx2(4,4) - qratio_dx(4)*&
                   (qratio_dx(4)/onepqs))/onepq
              if(iz.eq.2) then
                 ! N.B. x(1) (or extrasum(4)) dependence divides out of qratio
                 ! n.b. indices are reordered here so that
                 ! qratio_dx(k) refers to derivative wrt extrasum(k) for
                 ! k = 1, 2, 3, and qratio_dx(4) refers to the derivative
                 ! wrt extrasum(nextrasum-1).
                 qratio_dx(1) =&
                      qratio*(-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                 qratio_dx(2) =&
                      qratio*(-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                 qratio_dx(3) =&
                      qratio*(-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                 qratiot_dx(1) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstartx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratiot_dx(2) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstartx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratiot_dx(3) = qratiot*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstartx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(1,1) = qratio_dx(1)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(4,4,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(2,1) = qratio_dx(2)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(4,3,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(3,1) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(4,2,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(4,1) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*r_neutral(ielement) +&
                      qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(5,4,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(2,2) = qratio_dx(2)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(3,3,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(3,2) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(3,2,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(4,2) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      r_neutral(ielement)*3._fp_kind +&
                      qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(5,3,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(3,3) = qratio_dx(3)*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(2,2,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 qratio_dx2(4,3) = qratio_dx(4)*&
                      (-occ_const*r_neutral(ielement)*&
                      3._fp_kind +&
                      qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)) +&
                      qratio*(qstarx2(5,2,nmin_s,izqstar)/qstar(nmin_s,izqstar)-(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))*(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar)))
                 !                derivative of qratio wrt ln V
                 qratiov = qratiov -&
                      qratio_dx(1)*extrasum(1) -&
                      qratio_dx(2)*extrasum(2) -&
                      qratio_dx(3)*extrasum(3)
                 qratiovt = qratiovt -&
                      qratiot_dx(1)*extrasum(1) -&
                      qratiot_dx(2)*extrasum(2) -&
                      qratiot_dx(3)*extrasum(3)
                 qratiov_dx(1) = -&
                      qratio_dx2(1,1)*extrasum(1) -&
                      qratio_dx2(2,1)*extrasum(2) -&
                      qratio_dx2(3,1)*extrasum(3) -&
                      qratio_dx2(4,1)*extrasum(nextrasum-1) -&
                      qratio_dx(1)
                 qratiov_dx(2) = -&
                      qratio_dx2(2,1)*extrasum(1) -&
                      qratio_dx2(2,2)*extrasum(2) -&
                      qratio_dx2(3,2)*extrasum(3) -&
                      qratio_dx2(4,2)*extrasum(nextrasum-1) -&
                      qratio_dx(2)
                 qratiov_dx(3) = -&
                      qratio_dx2(3,1)*extrasum(1) -&
                      qratio_dx2(3,2)*extrasum(2) -&
                      qratio_dx2(3,3)*extrasum(3) -&
                      qratio_dx2(4,3)*extrasum(nextrasum-1) -&
                      qratio_dx(3)
                 qratiov_dx(4) = qratiov_dx(4) -&
                      qratio_dx2(4,1)*extrasum(1) -&
                      qratio_dx2(4,2)*extrasum(2) -&
                      qratio_dx2(4,3)*extrasum(3)
                 ! transform second order quantities to derivative of
                 ! ln(1 + qratio)
                 qratiot_dx(1) = (qratiot_dx(1) - qratiot*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratiot_dx(2) = (qratiot_dx(2) - qratiot*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratiot_dx(3) = (qratiot_dx(3) - qratiot*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratio_dx2(1,1) = (qratio_dx2(1,1) - qratio_dx(1)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(2,1) = (qratio_dx2(2,1) - qratio_dx(2)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(3,1) = (qratio_dx2(3,1) - qratio_dx(3)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(4,1) = (qratio_dx2(4,1) - qratio_dx(4)*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratio_dx2(2,2) = (qratio_dx2(2,2) - qratio_dx(2)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(3,2) = (qratio_dx2(3,2) - qratio_dx(3)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(4,2) = (qratio_dx2(4,2) - qratio_dx(4)*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratio_dx2(3,3) = (qratio_dx2(3,3) - qratio_dx(3)*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratio_dx2(4,3) = (qratio_dx2(4,3) - qratio_dx(4)*&
                      (qratio_dx(3)/onepqs))/onepq
                 qratiov_dx(1) = (qratiov_dx(1) - qratiov*&
                      (qratio_dx(1)/onepqs))/onepq
                 qratiov_dx(2) = (qratiov_dx(2) - qratiov*&
                      (qratio_dx(2)/onepqs))/onepq
                 qratiov_dx(3) = (qratiov_dx(3) - qratiov*&
                      (qratio_dx(3)/onepqs))/onepq
                 ! transform first order quantities to derivative of
                 ! ln(1 + qratio)
                 qratio_dx(1) = qratio_dx(1)/onepq
                 qratio_dx(2) = qratio_dx(2)/onepq
                 qratio_dx(3) = qratio_dx(3)/onepq
              endif ! if(iz.eq.2) then...
              ! transform second order quantities to derivative of
              ! ln(1 + qratio)
              qratiovt = (qratiovt - qratiov*&
                   (qratiot/onepqs))/onepq
              qratiov_dx(4) = (qratiov_dx(4) - qratiov*&
                   (qratio_dx(4)/onepqs))/onepq
              ! transform first order quantities to derivative of
              ! ln(1 + qratio)
              qratio_dx(4) = qratio_dx(4)/onepq
              qratiov = qratiov/onepq
           endif ! if(ifmhd_logical) then
           ! transform first order quantities to derivative of
           ! ln(1 + qratio)
           qratiot = qratiot/onepq
           if(qratio_us.gt.1.e-3_fp_kind) then
              qratio = qratio_scale*log(onepq)
           else
              ! alternating series so relative error is less than first
              ! missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
              qratio = qratio*&
                   (1._fp_kind      - qratio_us*&
                   (1._fp_kind/2._fp_kind - qratio_us*&
                   (1._fp_kind/3._fp_kind - qratio_us*&
                   (1._fp_kind/4._fp_kind - qratio_us*&
                   (1._fp_kind/5._fp_kind)))))
           endif
           if(ifmhd_logical) then
              ! calculate derivative *correction* due to x dependence on fl
              ! and tl.
              if(ifnr03) then
                 dqratiof = qratio_dx(4)*extrasumf(nextrasum-1)
                 dqratiot = qratio_dx(4)*extrasumt(nextrasum-1)
                 dqratiotf = qratiot_dx(4)*extrasumf(nextrasum-1)
                 dqratiott = qratiot_dx(4)*extrasumt(nextrasum-1)
                 dqratio_dxf(4) = qratio_dx2(4,4)*&
                      extrasumf(nextrasum-1)
                 dqratio_dxt(4) = qratio_dx2(4,4)*&
                      extrasumt(nextrasum-1)
                 dqratiovf = qratiov_dx(4)*extrasumf(nextrasum-1)
                 dqratiovt = qratiov_dx(4)*extrasumt(nextrasum-1)
              endif
              if(ifnr13) then
                 dqratio_dv(1:max_index) = qratio_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                 dqratiov_dv(1:max_index) = qratiov_dx(4)*extrasum_dv(1:max_index,nextrasum-1)
                 dqratio_dxdv(1:max_index,4) = qratio_dx2(4,4)*extrasum_dv(1:max_index,nextrasum-1)
              endif
              if(iz.eq.2) then
                 if(ifnr03) then
                    dqratiof = dqratiof + qratio_dx(1)*extrasumf(1) + qratio_dx(2)*extrasumf(2) + qratio_dx(3)*extrasumf(3)
                    dqratiot = dqratiot + qratio_dx(1)*extrasumt(1) + qratio_dx(2)*extrasumt(2) + qratio_dx(3)*extrasumt(3)
                    dqratiotf = dqratiotf + qratiot_dx(1)*extrasumf(1) + qratiot_dx(2)*extrasumf(2) + qratiot_dx(3)*extrasumf(3)
                    dqratiott = dqratiott + qratiot_dx(1)*extrasumt(1) + qratiot_dx(2)*extrasumt(2) + qratiot_dx(3)*extrasumt(3)
                    dqratio_dxf(1) =&
                         qratio_dx2(1,1)*extrasumf(1) +&
                         qratio_dx2(2,1)*extrasumf(2) +&
                         qratio_dx2(3,1)*extrasumf(3) +&
                         qratio_dx2(4,1)*extrasumf(nextrasum-1)
                    dqratio_dxt(1) =&
                         qratio_dx2(1,1)*extrasumt(1) +&
                         qratio_dx2(2,1)*extrasumt(2) +&
                         qratio_dx2(3,1)*extrasumt(3) +&
                         qratio_dx2(4,1)*extrasumt(nextrasum-1)
                    dqratio_dxf(2) =&
                         qratio_dx2(2,1)*extrasumf(1) +&
                         qratio_dx2(2,2)*extrasumf(2) +&
                         qratio_dx2(3,2)*extrasumf(3) +&
                         qratio_dx2(4,2)*extrasumf(nextrasum-1)
                    dqratio_dxt(2) =&
                         qratio_dx2(2,1)*extrasumt(1) +&
                         qratio_dx2(2,2)*extrasumt(2) +&
                         qratio_dx2(3,2)*extrasumt(3) +&
                         qratio_dx2(4,2)*extrasumt(nextrasum-1)
                    dqratio_dxf(3) =&
                         qratio_dx2(3,1)*extrasumf(1) +&
                         qratio_dx2(3,2)*extrasumf(2) +&
                         qratio_dx2(3,3)*extrasumf(3) +&
                         qratio_dx2(4,3)*extrasumf(nextrasum-1)
                    dqratio_dxt(3) =&
                         qratio_dx2(3,1)*extrasumt(1) +&
                         qratio_dx2(3,2)*extrasumt(2) +&
                         qratio_dx2(3,3)*extrasumt(3) +&
                         qratio_dx2(4,3)*extrasumt(nextrasum-1)
                    dqratio_dxf(4) = dqratio_dxf(4) +&
                         qratio_dx2(4,1)*extrasumf(1) +&
                         qratio_dx2(4,2)*extrasumf(2) +&
                         qratio_dx2(4,3)*extrasumf(3)
                    dqratio_dxt(4) = dqratio_dxt(4) +&
                         qratio_dx2(4,1)*extrasumt(1) +&
                         qratio_dx2(4,2)*extrasumt(2) +&
                         qratio_dx2(4,3)*extrasumt(3)
                    dqratiovf = dqratiovf +&
                         qratiov_dx(1)*extrasumf(1) +&
                         qratiov_dx(2)*extrasumf(2) +&
                         qratiov_dx(3)*extrasumf(3)
                    dqratiovt = dqratiovt +&
                         qratiov_dx(1)*extrasumt(1) +&
                         qratiov_dx(2)*extrasumt(2) +&
                         qratiov_dx(3)*extrasumt(3)
                 endif
                 if(ifnr13) then
                    dqratio_dv(1:max_index) =&
                         dqratio_dv(1:max_index) +&
                         qratio_dx(1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx(2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx(3)*extrasum_dv(1:max_index,3)
                    dqratiov_dv(1:max_index) =&
                         dqratiov_dv(1:max_index) +&
                         qratiov_dx(1)*extrasum_dv(1:max_index,1) +&
                         qratiov_dx(2)*extrasum_dv(1:max_index,2) +&
                         qratiov_dx(3)*extrasum_dv(1:max_index,3)
                    dqratio_dxdv(1:max_index,1) =&
                         qratio_dx2(1,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(2,1)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,1)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,1)*extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,2) =&
                         qratio_dx2(2,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(2,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,2)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,2)*extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,3) =&
                         qratio_dx2(3,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(3,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(3,3)*extrasum_dv(1:max_index,3) +&
                         qratio_dx2(4,3)*extrasum_dv(1:max_index,nextrasum-1)
                    dqratio_dxdv(1:max_index,4) = dqratio_dxdv(1:max_index,4) +&
                         qratio_dx2(4,1)*extrasum_dv(1:max_index,1) +&
                         qratio_dx2(4,2)*extrasum_dv(1:max_index,2) +&
                         qratio_dx2(4,3)*extrasum_dv(1:max_index,3)
                 endif ! if(ifnr03) then
              endif ! if(iz.eq.2) then
           endif ! if(ifmhd_logical) then
           ! H2plus component of excitation sums.
           extrace_tag=[26,1]
              call exsum_component_add(&
                ifmhd_logical, iz.eq.2, ifnr03, ifnr13,.true.,&
                1, max_index,&
                1, max_index,&
                inv_ion,&
                nuh2plus, nuh2plusf, nuh2plust, nuh2plus_dv,&
                xextrasum_scale, qratio_scale,&
                qratio, qratiot, qratio_dx,&
                qratiot2, qratiot_dx,&
                qratio_dx2,&
                qratiov, qratiovt, qratiov_dx,&
                dqratiof, dqratiotf, dqratio_dxf,&
                dqratiot, dqratiott, dqratio_dxt,&
                dqratio_dv,&
                dqratio_dxdv,&
                dqratiovf, dqratiovt, dqratiov_dv,&
                psum, psumf, psumt, psum_dv,&
                ssum, ssumf, ssumt,&
                usum,&
                free_sum, free_sumf, free_sum_dv,&
                xextrasum, xextrasumf, xextrasumt, xextrasum_dv)
        endif ! if(qratio.gt.exp_lim) then
     endif ! if(ifh2.gt.0.and.ifh2plus.gt.0.and.... now do H2+
  endif ! if(partial_element.... non-zero H abundance
end subroutine excitation_sum

! add components to each excitation sum and its derivatives.

!> This exsum_component_add subroutine increments various
!> number-density sums and their partial derivatives that are
!> associated with the excited state pressure ionization component of
!> the free energy.  This subroutine is a private module procedure
!> called in several places by excitation_sum to reduce code
!> duplication in that routine.
!>
!> \param[in] ifmhd_logical PARAMETERS NEED DOCUMENTATION
!>
subroutine exsum_component_add(&
     ifmhd_logical, iz_logical, ifnr03, ifnr13, ifsimple_index,&
     ion_start1, ion_end1,&
     ion_start2, ion_end2,&
     inv_ion,&
     nuvar, nuvarf, nuvart, nuvar_dv,&
     xextrasum_scale, qratio_scale,&
     qratio, qratiot, qratio_dx,&
     qratiot2, qratiot_dx,&
     qratio_dx2,&
     qratiov, qratiovt, qratiov_dx,&
     dqratiof, dqratiotf, dqratio_dxf,&
     dqratiot, dqratiott, dqratio_dxt,&
     dqratio_dv,&
     dqratio_dxdv,&
     dqratiovf, dqratiovt, dqratiov_dv,&
     psum, psumf, psumt, psum_dv,&
     ssum, ssumf, ssumt,&
     usum,&
     free_sum, free_sumf, free_sum_dv,&
     xextrasum, xextrasumf, xextrasumt, xextrasum_dv)

  use mod_excitation_block, only: extrace_count, extrace_tag, extrace_ids, &
       extrace_value, extrace_grad, extrace_hess, extrace_hf, extrace_ht, extrace_scale
  ! Arguments
  integer, intent(in) :: ion_start1, ion_end1, ion_start2, ion_end2, inv_ion(:)
  logical, intent(in) :: ifmhd_logical, iz_logical, ifnr03, ifnr13, ifsimple_index
  real(fp_kind), intent(in) :: nuvar, nuvarf, nuvart, nuvar_dv(:),&
       xextrasum_scale(:), qratio_scale,&
       qratio, qratiot, qratio_dx(:),&
       qratiot2, qratiot_dx(:),&
       qratio_dx2(:,:),&
       qratiov, qratiovt, qratiov_dx(:),&
       dqratiof, dqratiotf, dqratio_dxf(:),&
       dqratiot, dqratiott, dqratio_dxt(:),&
       dqratio_dv(:),&
       dqratio_dxdv(:,:),&
       dqratiovf, dqratiovt, dqratiov_dv(:)

  real(fp_kind), intent(inout) :: psum, psumf, psumt, psum_dv(:),&
       ssum, ssumf, ssumt,&
       usum,&
       free_sum, free_sumf, free_sum_dv(:),&
       xextrasum(:), xextrasumf(:), xextrasumt(:), xextrasum_dv(:,:)

  integer :: trace_i,trace_j,trace_slot
  ! Local variables:
  integer nxextrasum, ion_index, ion_index1, inv_ion_index, nionsp2
  real(fp_kind) nuvarsum_scale, nuvarfsum_scale,&
       nuvartsum_scale, nuvar_dvsum_scale,&
       nuvar_scale(4), nuvarf_scale(4),&
       nuvart_scale(4), nuvar_dv_scale(4)

  nionsp2 = size(inv_ion)
  nxextrasum = size(xextrasum_scale)

  ! Sanity checks
  if(&
       nionsp2.ne.size(dqratio_dv).or.&
       nionsp2.ne.size(dqratio_dxdv,1).or.&
       nionsp2.ne.size(dqratiov_dv).or.&
       nionsp2.ne.size(psum_dv).or.&
       nionsp2.ne.size(free_sum_dv).or.&
       nionsp2.ne.size(xextrasum_dv,1))&
       error stop 'exsum_component_add: inconsistent nionsp2 sizes for 7 arrays'
  if(&
       nxextrasum.ne.size(qratio_dx).or.&
       nxextrasum.ne.size(qratiot_dx).or.&
       nxextrasum.ne.size(qratio_dx2,1).or.&
       nxextrasum.ne.size(qratio_dx2,2).or.&
       nxextrasum.ne.size(qratiov_dx).or.&
       nxextrasum.ne.size(dqratio_dxf).or.&
       nxextrasum.ne.size(dqratio_dxt).or.&
       nxextrasum.ne.size(dqratio_dxdv,2).or.&
       nxextrasum.ne.size(xextrasum).or.&
       nxextrasum.ne.size(xextrasumf).or.&
       nxextrasum.ne.size(xextrasumt).or.&
       nxextrasum.ne.size(xextrasum_dv,2))&
       error stop 'exsum_component_add: inconsistent nxextrasum sizes for 13 arrays'

  extrace_count=extrace_count+1;trace_slot=extrace_count
  if(trace_slot.gt.318) error stop 'excitation trace capacity exceeded'
  extrace_ids(:,trace_slot)=extrace_tag
  extrace_scale(trace_slot)=qratio_scale
  extrace_value(:,trace_slot)=[nuvar,nuvarf,nuvart,qratio, &
       qratiot,qratiov]
  if(ifmhd_logical) then
    do trace_i=1,4
      if(iz_logical.or.trace_i.eq.4) then
        extrace_grad(trace_i,trace_slot)=qratio_dx(trace_i)
        if(ifnr03) then
          extrace_hf(trace_i,trace_slot)=dqratio_dxf(trace_i)
          extrace_ht(trace_i,trace_slot)=dqratio_dxt(trace_i)
        endif
        do trace_j=1,trace_i
          if(iz_logical.or.trace_j.eq.4) then
            extrace_hess(trace_i,trace_j,trace_slot)=qratio_dx2(trace_i,trace_j)
            extrace_hess(trace_j,trace_i,trace_slot)=extrace_hess(trace_i,trace_j,trace_slot)
          endif
        enddo
      endif
    enddo
  endif
  nuvarsum_scale = nuvar/qratio_scale
  nuvar_scale(1) = nuvar*(xextrasum_scale(1)/qratio_scale)
  nuvar_scale(2) = nuvar*(xextrasum_scale(2)/qratio_scale)
  nuvar_scale(3) = nuvar*(xextrasum_scale(3)/qratio_scale)
  nuvar_scale(4) = nuvar*(xextrasum_scale(4)/qratio_scale)
  if(ifnr03) then
     nuvarfsum_scale = nuvarf/qratio_scale
     nuvarf_scale(1) = nuvarf*(xextrasum_scale(1)/qratio_scale)
     nuvarf_scale(2) = nuvarf*(xextrasum_scale(2)/qratio_scale)
     nuvarf_scale(3) = nuvarf*(xextrasum_scale(3)/qratio_scale)
     nuvarf_scale(4) = nuvarf*(xextrasum_scale(4)/qratio_scale)
     nuvartsum_scale = nuvart/qratio_scale
     nuvart_scale(1) = nuvart*(xextrasum_scale(1)/qratio_scale)
     nuvart_scale(2) = nuvart*(xextrasum_scale(2)/qratio_scale)
     nuvart_scale(3) = nuvart*(xextrasum_scale(3)/qratio_scale)
     nuvart_scale(4) = nuvart*(xextrasum_scale(4)/qratio_scale)
  endif
  ! psum = sum n(species)/(rho*avogadro)
  !   partial delta ln Z wrt ln V
  ! psumf and psumt are derivatives (including dependence of n/rho
  !   on fl and tl) wrt fl and tl of psum.
  ! ssum = sum n(species)/(rho*avogadro) (delta ln Z +
  !   partial delta ln Z wrt ln t)
  ! ssumf and ssumt are derivatives (including dependence of n/rho
  !   on fl and tl) wrt fl and tl of ssum.
  ! usum = sum n(species)/(rho*avogadro)
  !   partial delta ln Z wrt ln t
  ssum = ssum + nuvarsum_scale*(qratio + qratiot)
  usum = usum + nuvarsum_scale*qratiot
  free_sum = free_sum + nuvarsum_scale*qratio
  if(ifmhd_logical) then
     psum = psum + nuvarsum_scale*qratiov
     if(iz_logical) then
        xextrasum(1) = xextrasum(1) - nuvar_scale(1)*qratio_dx(1)
        xextrasum(2) = xextrasum(2) - nuvar_scale(2)*qratio_dx(2)
        xextrasum(3) = xextrasum(3) - nuvar_scale(3)*qratio_dx(3)
     endif
     xextrasum(4) = xextrasum(4) - nuvar_scale(4)*qratio_dx(4)
  endif
  if(ifnr03) then
     ssumf = ssumf + nuvarfsum_scale*(qratio + qratiot)
     ssumt = ssumt + nuvartsum_scale*(qratio + qratiot) + nuvarsum_scale*(qratiot + qratiot2)
     free_sumf = free_sumf + nuvarfsum_scale*qratio
     if(ifmhd_logical) then
        psumf = psumf + nuvarfsum_scale*qratiov
        psumt = psumt + nuvartsum_scale*qratiov + nuvarsum_scale*qratiovt
        if(iz_logical) then
           xextrasumf(1) = xextrasumf(1) - nuvarf_scale(1)*qratio_dx(1) - nuvar_scale(1)*dqratio_dxf(1)
           xextrasumf(2) = xextrasumf(2) - nuvarf_scale(2)*qratio_dx(2) - nuvar_scale(2)*dqratio_dxf(2)
           xextrasumf(3) = xextrasumf(3) - nuvarf_scale(3)*qratio_dx(3) - nuvar_scale(3)*dqratio_dxf(3)
           xextrasumt(1) = xextrasumt(1) - nuvart_scale(1)*qratio_dx(1) -&
                nuvar_scale(1)*(qratiot_dx(1) + dqratio_dxt(1))
           xextrasumt(2) = xextrasumt(2) - nuvart_scale(2)*qratio_dx(2) -&
                nuvar_scale(2)*(qratiot_dx(2) + dqratio_dxt(2))
           xextrasumt(3) = xextrasumt(3) - nuvart_scale(3)*qratio_dx(3) -&
                nuvar_scale(3)*(qratiot_dx(3) + dqratio_dxt(3))
        endif
        xextrasumf(4) = xextrasumf(4) - nuvarf_scale(4)*qratio_dx(4) -&
             nuvar_scale(4)*dqratio_dxf(4)
        xextrasumt(4) = xextrasumt(4) - nuvart_scale(4)*qratio_dx(4) -&
             nuvar_scale(4)*(qratiot_dx(4) + dqratio_dxt(4))
        psumf = psumf + nuvarsum_scale*dqratiovf
        psumt = psumt + nuvarsum_scale*dqratiovt
        ssumf = ssumf + nuvarsum_scale*(dqratiof + dqratiotf)
        ssumt = ssumt + nuvarsum_scale*(dqratiot + dqratiott)
        free_sumf = free_sumf + nuvarsum_scale*dqratiof
     endif
  endif
  if(ifnr13) then
     do ion_index = ion_start1, ion_end1
        if(ifsimple_index) then
           inv_ion_index = ion_index
           ion_index1 = inv_ion_index
        else
           inv_ion_index = inv_ion(ion_index)
           ion_index1 = ion_index - ion_start1 + 1
        endif
        nuvar_dvsum_scale = nuvar_dv(ion_index1)/qratio_scale
        free_sum_dv(inv_ion_index) = free_sum_dv(inv_ion_index) + nuvar_dvsum_scale*qratio
     enddo
  endif
  if(ifmhd_logical.and.ifnr13) then
     do ion_index = ion_start1, ion_end1
        if(ifsimple_index) then
           inv_ion_index = ion_index
           ion_index1 = inv_ion_index
        else
           inv_ion_index = inv_ion(ion_index)
           ion_index1 = ion_index - ion_start1 + 1
        endif
        nuvar_dvsum_scale = nuvar_dv(ion_index1)/qratio_scale
        nuvar_dv_scale(1) = nuvar_dv(ion_index1)*(xextrasum_scale(1)/qratio_scale)
        nuvar_dv_scale(2) = nuvar_dv(ion_index1)*(xextrasum_scale(2)/qratio_scale)
        nuvar_dv_scale(3) = nuvar_dv(ion_index1)*(xextrasum_scale(3)/qratio_scale)
        nuvar_dv_scale(4) = nuvar_dv(ion_index1)*(xextrasum_scale(4)/qratio_scale)
        psum_dv(inv_ion_index) = psum_dv(inv_ion_index) + nuvar_dvsum_scale*qratiov
        if(iz_logical) then
           xextrasum_dv(inv_ion_index,1) = xextrasum_dv(inv_ion_index,1) - nuvar_dv_scale(1)*qratio_dx(1)
           xextrasum_dv(inv_ion_index,2) = xextrasum_dv(inv_ion_index,2) - nuvar_dv_scale(2)*qratio_dx(2)
           xextrasum_dv(inv_ion_index,3) = xextrasum_dv(inv_ion_index,3) - nuvar_dv_scale(3)*qratio_dx(3)
        endif
        xextrasum_dv(inv_ion_index,4) = xextrasum_dv(inv_ion_index,4) - nuvar_dv_scale(4)*qratio_dx(4)
     enddo
     do inv_ion_index = ion_start2, ion_end2
        free_sum_dv(inv_ion_index) = free_sum_dv(inv_ion_index) +&
             nuvarsum_scale*dqratio_dv(inv_ion_index)
        psum_dv(inv_ion_index) = psum_dv(inv_ion_index) +&
             nuvarsum_scale*dqratiov_dv(inv_ion_index)
        if(iz_logical) then
           xextrasum_dv(inv_ion_index,1) = xextrasum_dv(inv_ion_index,1) -&
                nuvar_scale(1)*dqratio_dxdv(inv_ion_index,1)
           xextrasum_dv(inv_ion_index,2) = xextrasum_dv(inv_ion_index,2) -&
                nuvar_scale(2)*dqratio_dxdv(inv_ion_index,2)
           xextrasum_dv(inv_ion_index,3) = xextrasum_dv(inv_ion_index,3) -&
                nuvar_scale(3)*dqratio_dxdv(inv_ion_index,3)
        endif
        xextrasum_dv(inv_ion_index,4) = xextrasum_dv(inv_ion_index,4) -&
             nuvar_scale(4)*dqratio_dxdv(inv_ion_index,4)
     enddo
  endif
end subroutine exsum_component_add
