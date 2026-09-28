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

! The purpose of this subroutine is to calculate the chemical potential/kT
! for each species (neutral and non-bare ions) for the free
! energy associated with excited states.  Afterward, calculate the
! corresponding change in the equilibrium constants dv(nions+2).
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
! R is the Rydberg constant
! plop = is the ground state Planck-Larkin occupation probability
! plop(n) = 1 - exp(-arg)*(1 + arg) (see notes for qstar_calc);
! piop is the ground state MHD occupation probability;
! piop(n) is the same thing for the excited states (see notes for qstar_calc).
! input quantities:
! ifexcited > 0 means use excited states (must have Planck-Larkin or ifpi = 3 or 4).
! 0 < ifexcited < 10 means use approximation to explicit summation
! 10 < ifexcited < 20 means use explicit summation
! mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules) and helium.
! mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
! mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
! ifpl = 1, use Planck-Larkin occupation probability, otherwise set to unity;
! ifpi = 3 or 4, use MHD occupation probability, otherwise set to unity;
! ifmodified > 0 (or not) affects ifpi_fit
! ifnr = 0, calculate fl and tl derivatives of output using chain rule and auxf and auxt
! ifnr = 1, calculate aux derivatives of output with fl, tl fixed
! ifnr = 2, calculate fl and tl derivatives of output with aux fixed.
! ifnr = 3 is combination of ifnr = 1 and ifnr = 2.
! ifsame_zero_abundances = 1 means izlo, izhi, bmin, nmin, nmin_max have not changed since last call.
! inv_aux(naux) is pointer to actual auxiliary variable index
! iextraoff (= 7) is the index offset between extrasum (*not xextrasum*) and aux.
! naux is number of active auxiliary variables
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
! nmax is the maximum principal quantum number included in sum
!   (to be compatible with opal which used nmax = 4 rather than infinity).
! if(nmax > 300000) then treated as infinity in qryd_approx.
!   otherwise nmax is meant to be used with qryd_calc only (i.e.,
!   case for mhd approximations not programmed).
! bion(nions+2) = ionization energy of next higher ion relative to
!   species that is being calculated.  nions+1 refers to H2,
!   nions+2 refers to H2+.
! plop(nions+2), plopt(nions+2), plopt2(nions+2) = *ln* of ground
!   state Planck-Larkin occupation probabilty and first and
!   second tl derivatives
! r_ion3(nions+2) is *the cube of the*
!   effective radii (in ion order but going from neutral
!   to next to bare ion for each species) of MDH interaction between
!   all but bare nucleii species and ionized species.  last 2 are
!   H2 and H2+ (note, only need r_ion3 of particular species)
! nion(nions+2), charge on ion in ion order (must be same order as bi)
!   e.g., for H+, He+, He++, etc.
! r_neutral(nelements+2) effective radii for MDH neutral-neutral
!   interactions.  Last two are for for H2 and H2+ (the only ionic
!   species in the MHD model with a non-zero hard-sphere radius).
! extrasum(nextrasum = 9) weighted sums over n(i)
!   for iextrasum = 1,nextrasum-2,
!   sum is only over neutral species (and H2+) and
!   weight is r_neutral^{iextrasum-1}
!   for iextrasum = nextrasum-1, sum is over all ionized species
!   including bare nucleii, but excluding free electrons,
!   weight is Z^1.5.
!   for iextrasum = nextrasum, sum is over all species
!   excluding bare nucleii and free electrons, the weight is rion^3.
! extrasumf(nextrasum) = partial of extrasum/partial ln f
! extrasumt(nextrasum) = partial of extrasum/partial ln t
! xextrasum(nxextrasum=4) is the *negative* sum nuvar/(1 + qratio)*
!   partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
! xextrasumf(nxextrasum) = partial of xextrasum/partial ln f
! xextrasumt(nxextrasum) = partial of xextrasum/partial ln t
! modified arrays:
! dv(nions+2), dvf(nions+2), dvt(nions+2), dv_aux(nions+2, naux)
!   equilibrium constants and derivatives (depending on ifnr).

!> This excitation_pi subroutine calculates the chemical potential/kT
!> and the partial derivaties of that quantity wrt fl, ln T, and the
!> excitation subset of auxiliary variables for each species (neutral
!> and non-bare ions) of each element for the Helmholz free-energy
!> component associated with excitation.  The corresponding changes
!> are added to the equilibrium constant array, dv, and its partial
!> derivatives, dvf, dvt, and dv_aux.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine excitation_pi(verbosity, ifexcited, ifsame_zero_abundances,&
     ifpl, ifpi, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,&
     partial_elements, ion_end,&
     tl, izlo, bmin, nmin, nmin_max, nmin_species, nmax,&
     bion, plop, plopt, plopt2,&
     r_ion3, nion, r_neutral,&
     inv_aux, iextraoff,&
     extrasum, extrasumf, extrasumt,&
     xextrasum, xextrasumf, xextrasumt,&
     dv, dvf, dvt, dv_aux)

  use mod_free_eos_constants, only: pi, c2, electron_mass, h_mass
  use mod_molecular_hydrogen, only: molecular_hydrogen
  use mod_excitation_block, only: nx, c2t, ifapprox_old_excitation,&
       ifdiff_x_excitation, ifh2_old_excitation, ifh2plus_old_excitation, ifhe1_special,&
       ifmhd_logical_old_excitation, ifpi_fit_old_excitation, ifpl_logical_old_excitation,&
       max_izhi, max_nmin_max,&
       qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,&
       qmhd_he1, qmhd_he1t, qmhd_he1x, qmhd_he1t2, qmhd_he1tx, qmhd_he1x2,&
       tl_old_excitation, x, x_old_excitation,&
       qstar, qstart, qstarx, qstart2, qstartx, qstarx2
  use mod_helium1_data, only: ell_helium1, neff_helium1, weight_helium1
  use mod_statistical_weight_data, only: nions_stat, nelements_stat, iqion, iqneutral
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  integer, intent(in) :: verbosity, ifexcited, ifpl, ifpi, ifmodified,&
       ifnr, ifsame_zero_abundances,&
       inv_aux(:), iextraoff, ifh2, ifh2plus,&
       partial_elements(:),&
       ion_end(:),&
       izlo, nmin(:), nmin_max(:), nmax,&
       nion(:), nmin_species(:),&
       inv_ion(:)
  real(fp_kind), intent(in) ::&
       tl, bmin(:), bion(:),&
       plop(:), plopt(:), plopt2(:),&
       r_ion3(:), r_neutral(:),&
       extrasum(:), extrasumf(:), extrasumt(:),&
       xextrasum(:), xextrasumf(:), xextrasumt(:)
  real(fp_kind), intent(inout) ::&
       dv(:), dvf(:), dvt(:), dv_aux(:,:)

  ! Parameters
  ! number of xextrasum variables
  integer, parameter :: nxextrasum = 4
  integer, parameter :: maxnextrasum = 9
  ! minimum value of principal quantum number where approximations are
  ! used
  integer, parameter :: min_nmin_max = 3
  ! minimum rydberg level for neutral helium (lower levels treated with
  ! exact neff, weight, and ell.
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
  logical ifmhd_logical, ifpl_logical, ifapprox, ifneutral
  logical ifnr13, ifnr0, ifnr023

  integer&
       n_partial_elements, izhi, nelements, nionsp2, nions, nextrasum, naux,&
       ! next line for temporary test
       ! ion_index, ion_index1,&
       inv_ion_index
  integer iz, izqstar, index, ielement, ion_start, ion, nmin_s, ifpi_fit
  integer&
       !n,&
       ix, nmin_max_max
  integer index_aux

  integer index_x(nx)
  ! Use data statement rather than parameter because last element needs to be set below at run time.
  data index_x /4, 3, 2, 1, 0/
  integer iffirst
  data iffirst/1/

  real(fp_kind) exparg, eps_factor
  real(fp_kind) qratio_us, onepq
  real(fp_kind) rfactor, ln_ground_occ

  real(fp_kind), allocatable ::&
       mu(:),&
       muf(:),&
       mut(:),&
       mu_aux(:, :)

  ! Most/all Fortran compilers specify the save attribute for all variables
  ! intialized by data statements.  But just in case...
  save iffirst, index_x

  if(iffirst.eq.1) then
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
  naux = size(inv_aux)
  ifpl_logical = ifpl.eq.1
  ifmhd_logical = ifpi.eq.3.or.ifpi.eq.4

  ! Sanity checks:
  if(nions.ne.nions_stat) error stop 'excitation_pi: inconsistent nions values'
  if(ifnr.lt.0.or.ifnr.gt.3) error stop 'excitation_pi: ifnr must be 0, 1, 2, or 3'
  if(ifexcited.le.0) error stop 'excitation_pi: invalid ifexcited'
  if(.not.(ifpl_logical.or.ifmhd_logical)) error stop 'excitation_pi: invalid combination of ifpl and ifpi'
  if(nelements.ne.nelements_stat) error stop 'excitation_pi: inconsistent nelements values'
  if(izhi+1.gt.max_izhi) error stop 'excitation_pi: izhi too large'

  if(&
       nx.ne.size(x).or.&
       nx.ne.size(x_old_excitation)&
       )&
       error stop 'excitation_pi: inconsistent sizes for x and x_old_excitation'
  if(nelements+2.ne.size(ion_end)) error stop 'excitation_pi: inconsistent sizes for ion_end and r_neutral'
  if(izhi.ne.size(nmin).or.izhi.ne.size(nmin_max))&
       error stop 'excitation_pi: inconsistent sizes for bmin, nmin, or nmin_max'
  if(&
       nionsp2.ne.size(nmin_species).or.&
       nionsp2.ne.size(inv_ion).or.&
       nionsp2.ne.size(bion).or.&
       nionsp2.ne.size(plop).or.&
       nionsp2.ne.size(plopt).or.&
       nionsp2.ne.size(plopt2).or.&
       nionsp2.ne.size(r_ion3).or.&
       nionsp2.ne.size(dv_aux,1))&
       error stop 'excitation_pi: inconsistent nionsp2 sizes for 8 variables'
  if(&
       nextrasum.ne.size(extrasumf).or.&
       nextrasum.ne.size(extrasumt))&
       error stop 'excitation_pi: inconsistent sizes for extrasum, extrasumt, or, extrasumf'
  if(&
       nxextrasum.ne.size(xextrasum).or.&
       nxextrasum.ne.size(xextrasumf).or.&
       nxextrasum.ne.size(xextrasumt))&
       error stop 'excitation_pi: inconsistent sizes for xextrasum, xextrasumt, or xextrasumf'

  if(naux.ne.size(dv_aux,2)) error stop 'excitation_pi: inconsistent naux sizes for inv_aux and dv_aux'

  ifnr13 = ifnr.eq.1.or.ifnr.eq.3
  ifnr0 = ifnr.eq.0
  ifnr023 = ifnr.eq.0.or.ifnr.eq.2.or.ifnr.eq.3
  ifapprox = ifexcited.lt.10

  allocate(&
       mu(nionsp2),&
       muf(nionsp2),&
       mut(nionsp2),&
       mu_aux(2*nxextrasum, nionsp2))

  if(if_taint_allocated_real) then
     call taint_allocated_real(mu)
     call taint_allocated_real(muf)
     call taint_allocated_real(mut)
     call taint_allocated_real(mu_aux)
  endif

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
        if(nmin_max(iz).gt.max_nmin_max) error stop 'excitation_pi: nmin_max too large'
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
        !   qstar(nmin(iz):nmin_max(iz),iz) = 0.d0
        !   if(iz.eq.2.and.ifh2plus.gt.0) then
        !     qstar(nmin(iz):nmin_max(iz),izhi+1) = 0.d0
        !   endif
        endif ! if(iz.eq.2.and.ifh2plus.gt.0) then
     enddo ! do iz = izlo, izhi
  endif !test for previous calculation of qh2, qh2plus, and qstar
  ! End of blocks of code that should be identical (aside from error
  ! stop identifications) for excitation_pi and excitation_sum.

  ! This do loop is over atomic species only
  do index = 1, n_partial_elements
     ielement = partial_elements(index)
     if(ielement.gt.1) then
        ion_start = ion_end(ielement-1) + 1
     else
        ion_start = 1
     endif
     ! calculate chemical potential/kt of neutral through next to bare ion
     do ion = ion_start, ion_end(ielement)
        iz = nion(ion)
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
        endif ! if(ifmhd_logical) then
        ! if(exparg.gt.exparg_lim) then
        ! Keep the Boltzmann factor logarithmic until the partition product.
        ! else
        !   exparg = 0.d0
        ! endif
        if(ielement.le.2.or.mod(ifexcited,10).eq.3) then
           nmin_s = nmin_species(ion)
           if(nmin_s.lt.nmin(iz).or.nmin_s.gt.nmin_max(iz)) error stop 'excitation_pi: (1) invalid nmin or nmin_max'
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
                 mu(ion) = partition_product(exparg,qmhd_he1,0._fp_kind)/real((iqneutral(ielement)),fp_kind)
                 ! else
                 !   mu(ion) = 0.d0
                 ! endif
              else ! if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
                 ! Not special case for helium
                 mu(ion) = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))*real((iqion(ion)),fp_kind)/real((iqneutral(ielement)),fp_kind)
              endif ! if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
           else ! if(iz.eq.1) then
              ! ionized atomic case
              mu(ion) = partition_product(exparg,qstar(nmin_s,iz),qstar_logscale(nmin_s,iz))*real((iqion(ion)),fp_kind)/real((iqion(ion-1)),fp_kind)
           endif ! if(iz.eq.1) then

           if(mu(ion).gt.exp_lim) then
              qratio_us = mu(ion)/qratio_scale
              onepq = 1._fp_kind + qratio_us
              if(ifnr023) then
                 ! calculate fl, tl derivatives with extrasum fixed
                 muf(ion) = 0._fp_kind
                 if(ifhe1_special.and.ion.eq.2.and.ifmhd_logical) then
                    ! special for neutral helium
                    mut(ion) = mu(ion)*(c2t*bion(ion) - plopt(ion) + qmhd_he1t/qmhd_he1)
                 else
                    mut(ion) = mu(ion)*(c2t*bion(ion) - plopt(ion) + qstart(nmin_s,iz)/qstar(nmin_s,iz))
                 endif
              endif
              if(ifmhd_logical) then
                 if(iz.eq.1) then
                    if(ifnr0) then
                       muf(ion) = muf(ion) -&
                            occ_const*mu(ion)*r_neutral(ielement)*(&
                            3._fp_kind*extrasumf(3) + r_neutral(ielement)*(&
                            3._fp_kind*extrasumf(2) + r_neutral(ielement)*(&
                            extrasumf(1))))
                       mut(ion) = mut(ion) -&
                            occ_const*mu(ion)*r_neutral(ielement)*(&
                            3._fp_kind*extrasumt(3) + r_neutral(ielement)*(&
                            3._fp_kind*extrasumt(2) + r_neutral(ielement)*(&
                            extrasumt(1))))
                       if(ifhe1_special.and.ion.eq.2) then
                          ! special for neutral helium
                          muf(ion) = muf(ion) +&
                               mu(ion)/qmhd_he1*(&
                               qmhd_he1x(4)*extrasumf(1) +&
                               qmhd_he1x(3)*extrasumf(2) +&
                               qmhd_he1x(2)*extrasumf(3))
                          mut(ion) = mut(ion) +&
                               mu(ion)/qmhd_he1*(&
                               qmhd_he1x(4)*extrasumt(1) +&
                               qmhd_he1x(3)*extrasumt(2) +&
                               qmhd_he1x(2)*extrasumt(3))
                       else
                          muf(ion) = muf(ion) +&
                               mu(ion)/qstar(nmin_s,iz)*(&
                               qstarx(4,nmin_s,iz)*extrasumf(1) +&
                               qstarx(3,nmin_s,iz)*extrasumf(2) +&
                               qstarx(2,nmin_s,iz)*extrasumf(3))
                          mut(ion) = mut(ion) +&
                               mu(ion)/qstar(nmin_s,iz)*(&
                               qstarx(4,nmin_s,iz)*extrasumt(1) +&
                               qstarx(3,nmin_s,iz)*extrasumt(2) +&
                               qstarx(2,nmin_s,iz)*extrasumt(3))
                       endif ! if(ifhe1_special.and.ion.eq.2) then
                    elseif(ifnr13) then ! if(ifnr0) then
                       mu_aux(3,ion) = -occ_const*mu(ion)*r_neutral(ielement)*3._fp_kind
                       mu_aux(2,ion) = mu_aux(3,ion)*r_neutral(ielement)
                       mu_aux(1,ion) = mu_aux(2,ion)*r_neutral(ielement)/3._fp_kind
                       if(ifhe1_special.and.ion.eq.2) then
                          ! special for neutral helium
                          mu_aux(3,ion) = mu_aux(3,ion) + mu(ion)*(qmhd_he1x(2)/qmhd_he1)
                          mu_aux(2,ion) = mu_aux(2,ion) + mu(ion)*(qmhd_he1x(3)/qmhd_he1)
                          mu_aux(1,ion) = mu_aux(1,ion) + mu(ion)*(qmhd_he1x(4)/qmhd_he1)
                       else
                          mu_aux(3,ion) = mu_aux(3,ion) + mu(ion)*(qstarx(2,nmin_s,iz)/qstar(nmin_s,iz))
                          mu_aux(2,ion) = mu_aux(2,ion) + mu(ion)*(qstarx(3,nmin_s,iz)/qstar(nmin_s,iz))
                          mu_aux(1,ion) = mu_aux(1,ion) + mu(ion)*(qstarx(4,nmin_s,iz)/qstar(nmin_s,iz))
                       endif
                    endif ! if(ifnr0) then
                 endif ! if(iz.eq.1) then
                 if(ifnr0) then
                    muf(ion) = muf(ion) -occ_const*mu(ion)*r_ion3(ion)*extrasumf(nextrasum-1)
                    mut(ion) = mut(ion) - occ_const*mu(ion)*r_ion3(ion)*extrasumt(nextrasum-1)
                    if(ifhe1_special.and.ion.eq.2) then
                       ! special for neutral helium
                       muf(ion) = muf(ion) + mu(ion)/qmhd_he1*(qmhd_he1x(5)*extrasumf(nextrasum-1))
                       mut(ion) = mut(ion) + mu(ion)/qmhd_he1*(qmhd_he1x(5)*extrasumt(nextrasum-1))
                    else
                       muf(ion) = muf(ion) + mu(ion)/qstar(nmin_s,iz)*(qstarx(5,nmin_s,iz)*extrasumf(nextrasum-1))
                       mut(ion) = mut(ion) + mu(ion)/qstar(nmin_s,iz)*(qstarx(5,nmin_s,iz)*extrasumt(nextrasum-1))
                    endif
                 elseif(ifnr13) then ! if(ifnr0) then
                    mu_aux(4,ion) = -occ_const*mu(ion)*r_ion3(ion)
                    if(ifhe1_special.and.ion.eq.2) then
                       ! special for neutral helium
                       mu_aux(4,ion) = mu_aux(4,ion) + mu(ion)*(qmhd_he1x(5)/qmhd_he1)
                    else
                       mu_aux(4,ion) = mu_aux(4,ion) + mu(ion)*(qstarx(5,nmin_s,iz)/qstar(nmin_s,iz))
                    endif
                 endif ! if(ifnr0) then
              endif ! if(ifmhd_logical) then
              if(ifnr023) then
                 ! transform to chemical potential/kt
                 muf(ion) = -muf(ion)/onepq
                 mut(ion) = -mut(ion)/onepq
              endif
              if(ifnr13.and.ifmhd_logical) then
                 if(iz.eq.1) then
                    mu_aux(1,ion) = -mu_aux(1,ion)/onepq
                    mu_aux(2,ion) = -mu_aux(2,ion)/onepq
                    mu_aux(3,ion) = -mu_aux(3,ion)/onepq
                 endif
                 mu_aux(4,ion) = -mu_aux(4,ion)/onepq
              endif
              if(qratio_us.gt.1.e-3_fp_kind) then
                 mu(ion) = -qratio_scale*log(onepq)
              else
                 ! alternating series so relative error is less than first
                 ! missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
                 mu(ion) = -mu(ion)*&
                      (1._fp_kind      - qratio_us*&
                      (1._fp_kind/2._fp_kind - qratio_us*&
                      (1._fp_kind/3._fp_kind - qratio_us*&
                      (1._fp_kind/4._fp_kind - qratio_us*&
                      (1._fp_kind/5._fp_kind)))))
              endif
           else ! if(mu(ion).gt.exp_lim) then
              mu(ion) = 0._fp_kind
              if(ifnr023) then
                 muf(ion) = 0._fp_kind
                 mut(ion) = 0._fp_kind
              endif
              if(ifnr13.and.ifmhd_logical) then
                 if(iz.eq.1) then
                    mu_aux(1,ion) = 0._fp_kind
                    mu_aux(2,ion) = 0._fp_kind
                    mu_aux(3,ion) = 0._fp_kind
                 endif
                 mu_aux(4,ion) = 0._fp_kind
              endif
           endif !if(mu(ion).gt.exp_lim) then
        else ! if(ielement.le.2.or.mod(ifexcited,10).eq.3) then
           mu(ion) = 0._fp_kind
           if(ifnr023) then
              muf(ion) = 0._fp_kind
              mut(ion) = 0._fp_kind
           endif
           if(ifnr13.and.ifmhd_logical) then
              if(iz.eq.1) then
                 mu_aux(1,ion) = 0._fp_kind
                 mu_aux(2,ion) = 0._fp_kind
                 mu_aux(3,ion) = 0._fp_kind
              endif
              mu_aux(4,ion) = 0._fp_kind
           endif
        endif ! if(ielement.le.2.or.mod(ifexcited,10).eq.3) then

        if(ifmhd_logical) then
           if(iz.eq.1) then
              mu(ion) = mu(ion) + qratio_scale*(&
                   xextrasum(1) +&
                   xextrasum(2)*r_neutral(ielement) +&
                   xextrasum(3)*r_neutral(ielement)*r_neutral(ielement))
              if(ifnr0) then
                 muf(ion) = muf(ion) + qratio_scale*(&
                      xextrasumf(1) +&
                      xextrasumf(2)*r_neutral(ielement) +&
                      xextrasumf(3)*r_neutral(ielement)*&
                      r_neutral(ielement))
                 mut(ion) = mut(ion) + qratio_scale*(&
                      xextrasumt(1) +&
                      xextrasumt(2)*r_neutral(ielement) +&
                      xextrasumt(3)*r_neutral(ielement)*&
                      r_neutral(ielement))
              elseif(ifnr13) then
                 mu_aux(5,ion) = qratio_scale
                 mu_aux(6,ion) = qratio_scale*r_neutral(ielement)
                 mu_aux(7,ion) = qratio_scale*r_neutral(ielement)*&
                      r_neutral(ielement)
              endif
           else ! if(iz.eq.1) then
              mu(ion) = mu(ion) + qratio_scale*&
                   xextrasum(4)*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
              if(ifnr0) then
                 muf(ion) = muf(ion) + qratio_scale*&
                      xextrasumf(4)*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
                 mut(ion) = mut(ion) + qratio_scale*&
                      xextrasumt(4)*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
              elseif(ifnr13) then
                 mu_aux(8,ion) = qratio_scale*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
              endif
           endif ! if(iz.eq.1) then
        endif ! if(ifmhd_logical) then
        ! ! temporary test of derivatives
        ! dv(ion) = mu(ion)
        ! inv_ion_index = inv_ion(ion)
        ! if(iz.eq.1) then
        !   index_aux = inv_aux(iextraoff+1)
        !   dv_aux(inv_ion_index,index_aux) = mu_aux(1,ion)
        !   index_aux = inv_aux(iextraoff+2)
        !   dv_aux(inv_ion_index,index_aux) = mu_aux(2,ion)
        !   index_aux = inv_aux(iextraoff+3)
        !   dv_aux(inv_ion_index,index_aux) = mu_aux(3,ion)
        ! else
        !   index_aux = inv_aux(iextraoff+1)
        !   dv_aux(inv_ion_index,index_aux) = 0.d0
        !   index_aux = inv_aux(iextraoff+2)
        !   dv_aux(inv_ion_index,index_aux) = 0.d0
        !   index_aux = inv_aux(iextraoff+3)
        !   dv_aux(inv_ion_index,index_aux) = 0.d0
        ! endif
        ! index_aux = inv_aux(iextraoff+nextrasum-1)
        ! dv_aux(inv_ion_index,index_aux) = mu_aux(4,ion)
     enddo ! do ion = ion_start, ion_end(ielement)

     ! calculate change in equilibrium constant of monatomic ions
     ! relative to neutral monatomic.
     ! n.b. ion now refers to equilibrium constant of h+, he+, he++, etc.,
     ! i.e. offset by 1 from previous ion meaning.
     do ion = ion_start, ion_end(ielement)
        if(ion.lt.ion_end(ielement)) then
           dv(ion) = dv(ion) + (mu(ion_start) - mu(ion+1))/qratio_scale
           if(ifnr023) then
              dvf(ion) = dvf(ion) + (muf(ion_start) - muf(ion+1))/qratio_scale
              dvt(ion) = dvt(ion) + (mut(ion_start) - mut(ion+1))/qratio_scale
           endif
           if(ifnr13.and.ifmhd_logical) then
              inv_ion_index = inv_ion(ion)
              index_aux = inv_aux(iextraoff+1)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(1,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+2)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(2,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+3)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(3,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+nextrasum-1)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   (mu_aux(4,ion_start) - mu_aux(4,ion+1))/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+1)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(5,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+2)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(6,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+3)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(7,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+4)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                   mu_aux(8,ion+1)/qratio_scale
           endif ! if(ifnr13.and.ifmhd_logical) then
        else ! if(ion.lt.ion_end(ielement)) then
           dv(ion) = dv(ion) + mu(ion_start)/qratio_scale
           if(ifmhd_logical) then
              ! n.b. bare ion has a chemical potential of xextrasum(4)*iz^{3/2}
              iz = nion(ion)
              dv(ion) = dv(ion) - xextrasum(4)*real((iz),fp_kind)*sqrt(real((iz),fp_kind))
              if(ifnr0) then
                 dvf(ion) = dvf(ion) - xextrasumf(4)*real((iz),fp_kind)*sqrt(real((iz),fp_kind))
                 dvt(ion) = dvt(ion) - xextrasumt(4)*real((iz),fp_kind)*sqrt(real((iz),fp_kind))
              endif
           endif
           if(ifnr023) then
              dvf(ion) = dvf(ion) + muf(ion_start)/qratio_scale
              dvt(ion) = dvt(ion) + mut(ion_start)/qratio_scale
           endif
           if(ifnr13.and.ifmhd_logical) then
              inv_ion_index = inv_ion(ion)
              index_aux = inv_aux(iextraoff+1)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(1,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+2)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(2,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+3)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(3,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+nextrasum-1)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(4,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+1)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(5,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+2)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(6,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+3)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                   mu_aux(7,ion_start)/qratio_scale
              index_aux = inv_aux(iextraoff+maxnextrasum+4)
              dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                   real((iz),fp_kind)*sqrt(real((iz),fp_kind))
           endif ! if(ifnr13.and.ifmhd_logical) then
        endif ! if(ion.lt.ion_end(ielement)) then
     enddo ! do ion = ion_start, ion_end(ielement)
  enddo ! do index = 1, n_partial_elements
  if(partial_elements(1).eq.1) then
     ! do all the molecular stuff below only if hydrogen has a non-zero
     ! abundance
     ! now do H2 equilibrium constant change due to neutral monatomic H.
     if(ifh2.gt.0) then
        ! calculate change in H2 equilibrium constant
        ! relative to neutral monatomic H
        ! ion refers to H2, ion_start refers to monatomic H.
        ion = nions+1
        ion_start = 1
        dv(ion) = dv(ion) + 2._fp_kind*mu(ion_start)/qratio_scale
        if(ifnr023) then
           dvf(ion) = dvf(ion) + 2._fp_kind*muf(ion_start)/qratio_scale
           dvt(ion) = dvt(ion) + 2._fp_kind*mut(ion_start)/qratio_scale
        endif
        if(ifnr13.and.ifmhd_logical) then
           inv_ion_index = inv_ion(ion)
           index_aux = inv_aux(iextraoff+1)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(1,ion_start)/qratio_scale
           index_aux = inv_aux(iextraoff+2)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(2,ion_start)/qratio_scale
           index_aux = inv_aux(iextraoff+3)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(3,ion_start)/qratio_scale
           index_aux = inv_aux(iextraoff+nextrasum-1)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(4,ion_start)/qratio_scale
           index_aux = inv_aux(iextraoff+maxnextrasum+1)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(5,ion_start)/qratio_scale
           index_aux = inv_aux(iextraoff+maxnextrasum+2)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(6,ion_start)/qratio_scale
           index_aux = inv_aux(iextraoff+maxnextrasum+3)
           dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                2._fp_kind*mu_aux(7,ion_start)/qratio_scale
        endif ! if(ifnr13.and.ifmhd_logical) then
     endif ! if(ifh2.gt.0) then
     ! now do H2 and H2plus chemical potentials/kT
     ! n.b. ifh2, ifh2plus, and ifexcited control whether *xextrasum* has
     ! a molecular component.  However, *only* ifh2 and ifh2plus control
     ! whether *extrasum* has a molecular component.
     if((ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1).or.ifmhd_logical) then
        do iz = 1,2
           if((iz.eq.1.and.ifh2.gt.0).or.(iz.eq.2.and.ifh2plus.gt.0)) then
              if(iz.eq.1) then
                 ielement = nelements + 1
                 ion = nions+1
                 izqstar = iz
              else
                 ielement = nelements + 2
                 ion = nions+2
                 ! special values of qstar and friends calculated including
                 ! "neutral" occupation probability.
                 izqstar = izhi+1
              endif
              nmin_s = nmin_species(ion)
              if(nmin_s.lt.nmin(iz).or.nmin_s.gt.nmin_max(iz)) error stop 'excitation_pi: (2) invalid nmin or nmin_max'
              if(ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1) then
                 if(iz.eq.1) then
                    exparg = exparg_shift - c2t*bion(ion) - plop(ion) + qh2plus - qh2
                 else
                    ! n.b. excited state core is two protons with unity
                    ! statistical weight (following usual convention that
                    ! nuclear spin statistical weights are divided out).
                    exparg = exparg_shift - c2t*bion(ion) - plop(ion) - qh2plus
                 endif
                 ! correct for MHD ground state occupation probability.
                 if(ifmhd_logical) then
                    ln_ground_occ = occ_const*(&
                         x(1) + r_neutral(ielement)*(&
                         x(2)*3._fp_kind + r_neutral(ielement)*(&
                         x(3)*3._fp_kind + r_neutral(ielement)*(&
                         x(4)))))
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
                 mu(ion) = partition_product(exparg,qstar(nmin_s,izqstar),qstar_logscale(nmin_s,izqstar))
              else ! if(ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1) then
                 mu(ion) = 0._fp_kind
              endif ! if(ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1) then
              if(mu(ion).gt.exp_lim) then
                 qratio_us = mu(ion)/qratio_scale
                 onepq = 1._fp_kind + qratio_us
                 if(ifnr023) then
                    ! calculate fl, tl derivatives with extrasum fixed
                    muf(ion) = 0._fp_kind
                    mut(ion) = mu(ion)*(c2t*bion(ion) - plopt(ion) +&
                         qstart(nmin_s,izqstar)/qstar(nmin_s,izqstar))
                    if(iz.eq.1) then
                       mut(ion) = mut(ion) + mu(ion)*(qh2plust - qh2t)
                    else
                       mut(ion) = mut(ion) + mu(ion)*(-qh2plust)
                    endif
                 endif
                 if(ifmhd_logical) then
                    if(ifnr0) then
                       muf(ion) = muf(ion) -&
                            occ_const*mu(ion)*r_neutral(ielement)*(&
                            3._fp_kind*extrasumf(3) + r_neutral(ielement)*(&
                            3._fp_kind*extrasumf(2) + r_neutral(ielement)*(&
                            extrasumf(1)))) +&
                            mu(ion)/qstar(nmin_s,izqstar)*(&
                            qstarx(4,nmin_s,izqstar)*extrasumf(1) +&
                            qstarx(3,nmin_s,izqstar)*extrasumf(2) +&
                            qstarx(2,nmin_s,izqstar)*extrasumf(3))
                       mut(ion) = mut(ion) -&
                            occ_const*mu(ion)*r_neutral(ielement)*(&
                            3._fp_kind*extrasumt(3) + r_neutral(ielement)*(&
                            3._fp_kind*extrasumt(2) + r_neutral(ielement)*(&
                            extrasumt(1)))) +&
                            mu(ion)/qstar(nmin_s,izqstar)*(&
                            qstarx(4,nmin_s,izqstar)*extrasumt(1) +&
                            qstarx(3,nmin_s,izqstar)*extrasumt(2) +&
                            qstarx(2,nmin_s,izqstar)*extrasumt(3))
                       muf(ion) = muf(ion) -&
                            occ_const*mu(ion)*r_ion3(ion)*extrasumf(nextrasum-1) +&
                            mu(ion)/qstar(nmin_s,izqstar)*(qstarx(5,nmin_s,izqstar)*extrasumf(nextrasum-1))
                       mut(ion) = mut(ion) -&
                            occ_const*mu(ion)*r_ion3(ion)*extrasumt(nextrasum-1) +&
                            mu(ion)/qstar(nmin_s,izqstar)*(qstarx(5,nmin_s,izqstar)*extrasumt(nextrasum-1))
                    elseif(ifnr13) then ! if(ifnr0) then
                       mu_aux(3,ion) = -occ_const*mu(ion)*r_neutral(ielement)*3._fp_kind
                       mu_aux(2,ion) = mu_aux(3,ion)*r_neutral(ielement)
                       mu_aux(1,ion) = mu_aux(2,ion)*r_neutral(ielement)/3._fp_kind
                       mu_aux(3,ion) = mu_aux(3,ion) + mu(ion)*(qstarx(2,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                       mu_aux(2,ion) = mu_aux(2,ion) + mu(ion)*(qstarx(3,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                       mu_aux(1,ion) = mu_aux(1,ion) + mu(ion)*(qstarx(4,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                       if(verbosity.ge.2.and.tl.lt.log(1000._fp_kind)) &
                            write(*,*) 'SCALED_PRODUCT',ion,mu(ion),qstarx(5,nmin_s,izqstar),qstar(nmin_s,izqstar)
                       mu_aux(4,ion) = -occ_const*mu(ion)*r_ion3(ion) + mu(ion)*(qstarx(5,nmin_s,izqstar)/qstar(nmin_s,izqstar))
                    endif ! if(ifnr0) then
                 endif ! if(ifmhd_logical) then
                 if(ifnr023) then
                    ! transform to chemical potential/kt
                    muf(ion) = -muf(ion)/onepq
                    mut(ion) = -mut(ion)/onepq
                 endif
                 if(ifnr13.and.ifmhd_logical) then
                    mu_aux(1,ion) = -mu_aux(1,ion)/onepq
                    mu_aux(2,ion) = -mu_aux(2,ion)/onepq
                    mu_aux(3,ion) = -mu_aux(3,ion)/onepq
                    mu_aux(4,ion) = -mu_aux(4,ion)/onepq
                 endif
                 if(qratio_us.gt.1.e-3_fp_kind) then
                    mu(ion) = -qratio_scale*log(onepq)
                 else
                    ! alternating series so relative error is less than first
                    ! missing term which is qratio_us^5/6 < (1.d-3)^-5/6 ~ 2.d-16.
                    mu(ion) = -mu(ion)*&
                         (1._fp_kind      - qratio_us*&
                         (1._fp_kind/2._fp_kind - qratio_us*&
                         (1._fp_kind/3._fp_kind - qratio_us*&
                         (1._fp_kind/4._fp_kind - qratio_us*&
                         (1._fp_kind/5._fp_kind)))))
                 endif
              else ! if(mu(ion).gt.exp_lim) then
                 mu(ion) = 0._fp_kind
                 if(ifnr023) then
                    muf(ion) = 0._fp_kind
                    mut(ion) = 0._fp_kind
                 endif
                 if(ifnr13.and.ifmhd_logical) then
                    mu_aux(1,ion) = 0._fp_kind
                    mu_aux(2,ion) = 0._fp_kind
                    mu_aux(3,ion) = 0._fp_kind
                    mu_aux(4,ion) = 0._fp_kind
                 endif
              endif ! if(mu(ion).gt.exp_lim) then
              if(ifmhd_logical) then
                 mu(ion) = mu(ion) + qratio_scale*(&
                      xextrasum(1) +&
                      xextrasum(2)*r_neutral(ielement) +&
                      xextrasum(3)*r_neutral(ielement)*r_neutral(ielement))
                 if(ifnr0) then
                    muf(ion) = muf(ion) + qratio_scale*(&
                         xextrasumf(1) +&
                         xextrasumf(2)*r_neutral(ielement) +&
                         xextrasumf(3)*r_neutral(ielement)*r_neutral(ielement))
                    mut(ion) = mut(ion) + qratio_scale*(&
                         xextrasumt(1) +&
                         xextrasumt(2)*r_neutral(ielement) +&
                         xextrasumt(3)*r_neutral(ielement)*r_neutral(ielement))
                 elseif(ifnr13) then
                    mu_aux(5,ion) = qratio_scale
                    mu_aux(6,ion) = qratio_scale*r_neutral(ielement)
                    mu_aux(7,ion) = qratio_scale*r_neutral(ielement)*r_neutral(ielement)
                 endif
                 if(iz.eq.2) then
                    mu(ion) = mu(ion) + qratio_scale*xextrasum(4)*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
                    if(ifnr0) then
                       muf(ion) = muf(ion) + qratio_scale*xextrasumf(4)*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
                       mut(ion) = mut(ion) + qratio_scale*xextrasumt(4)*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
                    elseif(ifnr13) then
                       mu_aux(8,ion) = qratio_scale*real((iz-1),fp_kind)*sqrt(real((iz-1),fp_kind))
                    endif
                 endif
              endif ! if(ifmhd_logical) then
              if(iz.eq.1) then
                 ! calculate change in H2 equilibrium constant
                 ! due to chemical potential of H2 (effect of chemical potential
                 ! of neutral H already applied).
                 dv(ion) = dv(ion) - mu(ion)/qratio_scale
                 if(ifnr023) then
                    dvf(ion) = dvf(ion) - muf(ion)/qratio_scale
                    dvt(ion) = dvt(ion) - mut(ion)/qratio_scale
                 endif
                 if(ifnr13.and.ifmhd_logical) then
                    inv_ion_index = inv_ion(ion)
                    index_aux = inv_aux(iextraoff+1)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(1,ion)/qratio_scale
                    index_aux = inv_aux(iextraoff+2)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(2,ion)/qratio_scale
                    index_aux = inv_aux(iextraoff+3)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(3,ion)/qratio_scale
                    index_aux = inv_aux(iextraoff+nextrasum-1)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(4,ion)/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+1)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(5,ion)/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+2)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(6,ion)/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+3)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(7,ion)/qratio_scale
                 endif ! if(ifnr13.and.ifmhd_logical) then
              else ! if(iz.eq.1) then
                 ! calculate change in H2+ equilibrium constant
                 ! relative to H2 and e-.
                 ion_start = nions + 1
                 ! ion_start refers to H2, ion refers to H2+.
                 dv(ion) = dv(ion) + (mu(ion_start) - mu(ion))/qratio_scale
                 if(ifnr023) then
                    dvf(ion) = dvf(ion) + (muf(ion_start) - muf(ion))/qratio_scale
                    dvt(ion) = dvt(ion) + (mut(ion_start) - mut(ion))/qratio_scale
                 endif
                 if(ifnr13.and.ifmhd_logical) then
                    inv_ion_index = inv_ion(ion)
                    index_aux = inv_aux(iextraoff+1)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(1,ion_start) - mu_aux(1,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+2)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(2,ion_start) - mu_aux(2,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+3)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(3,ion_start) - mu_aux(3,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+nextrasum-1)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(4,ion_start) - mu_aux(4,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+1)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(5,ion_start) - mu_aux(5,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+2)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(6,ion_start) - mu_aux(6,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+3)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) +&
                         (mu_aux(7,ion_start) - mu_aux(7,ion))/qratio_scale
                    index_aux = inv_aux(iextraoff+maxnextrasum+4)
                    dv_aux(inv_ion_index,index_aux) = dv_aux(inv_ion_index,index_aux) -&
                         mu_aux(8,ion)/qratio_scale
                 endif ! if(ifnr13.and.ifmhd_logical) then
              endif ! if(iz.eq.1) then
           endif ! if((iz.eq.1.and.ifh2.gt.0).or.(iz.eq.2.and.ifh2plus.gt.0)) then
        enddo ! do iz = 1,2
     endif ! if((ifh2.gt.0.and.ifh2plus.gt.0.and.mod(ifexcited,10).gt.1).or.ifmhd_logical) then
  endif ! if(partial_elements(1).eq.1) then
end subroutine excitation_pi
