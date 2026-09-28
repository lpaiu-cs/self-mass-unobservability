!> This module provides the public effective radius module procedure
!> for calculating effective pressure-ionization interaction radii and
!> also provides parameters that optionally (depending on option
!> suite) adjust the semi-empirical non-ideal Coulomb and pressure
!> ionization components of the free energy

module mod_pi_fit
  use mod_free_eos_types, only: fp_kind
  implicit none
  real(fp_kind), save, public :: pi_trace_neutral(26)=0._fp_kind, pi_trace_ion3(318)=0._fp_kind
  private
  public&
       xdh10, xmocp10,&
       pi_fitx_neutral_original, pi_fitx_neutral, pi_fitx_neutral_saumon,&
       pi_fitx_ion_original,  pi_fitx_ion, pi_fitx_ion_saumon,&
       quad_original, quad_modified,&
       effective_radius
  ! Coulomb fitting factors used to fit various detailed EOS results.
  ! These are the log10(gamma) limits of the DH region and the modified
  ! OCP region.
  real(fp_kind), parameter ::  xdh10 = -0.4_fp_kind, xmocp10 = 0._fp_kind

  ! Pressure-ionization fitting factors used to fit various detailed
  ! EOS results.
  integer, parameter :: nelements_pi_fitp2 = 26

  ! adjust these values for neutral-neutral pressure ionization

  ! best fit to MDH results.
  real(fp_kind), parameter :: pi_fit_neutral_ln_original(nelements_pi_fitp2) =&
       [0._fp_kind, 0._fp_kind, spread(0._fp_kind,1,22), spread(0._fp_kind,1,2)]

  real(fp_kind), parameter :: pi_fit_neutral_ln(nelements_pi_fitp2) =&
  ! best fit to opal results
  ! [0.5d0, 0.5d0, spread(0.5d0,1,18), -50.d0]
  ! attempt to fit Saumon H results
  ! [-0.050d0, 0.5d0, spread(0.5d0,1,18), -2.d0]
  ! keep neutral metals same as original mdh for simplicity
  ! [-0.050d0, 0.5d0, spread(0.d0,1,18), -2.d0]
  ! try neutrals identical to Saumon fit for opal fit, metals hydrogenic
  ! like excited-state pi_fit factors
  ! Force best fit, quad=0 results for Saumon
  ! [0.8d0, 0.d0, spread(0.d0,1,18), -2.7d0]
  ! Same for quad=10 Saumon fit.
  ! [-0.4d0, -0.2d0, spread(0.d0,1,18), 2*-4.0d0]
  ! Same for quad=20 Saumon fit.
  ! best pre eos2001 fit.
  ! [0.d0, -0.2d0, spread(0.d0,1,18), -3.6d0]
  ! eos2001 fit
  ! [-3.d0, -4.d0, spread(-5.d0,1,18), -5.d0, -5.d0]
  ! eos2001 fit (refined solar, Helium fit to start)
  ! [-0.4d0, -0.2d0, spread(0.d0,1,18), 2*-4.0d0]
  ! [-0.4d0, -4.d0, spread(0.d0,1,18), 2*-4.0d0]
  ! Go back to overall eos2001 fit since it gives reasonable fit
  ! for solar case.
  ! [-3.d0, -4.d0, spread(-5.d0,1,18), -5.d0, -5.d0]
  ! Revert back to 1.1.0
       [-0.4_fp_kind, -0.2_fp_kind, spread(0._fp_kind,1,22), -4.0_fp_kind, -40._fp_kind]
  real(fp_kind), parameter :: pi_fit_neutral_ln_saumon(nelements_pi_fitp2) =&
  ! This next version was a successful historical ad hoc experiment to
  ! reduce the influence of the molecule H2+ which currently (for the
  ! above pi_fit_neutral_ln values) causes convergence and
  ! discontinuity issues for the PI region as can be seen in the
  ! current convergence paper figures.  However, hold off on this
  ! change for now for the reasons mentioned in the commit message.
  ![-0.4d0, -0.2d0, spread(0.d0,1,18), -3.5d0, 0.d0]
  ! fit to Saumon H, and He tables, metals don't matter
  ! Best fit, quad=0 results.
  ! [0.8d0, 0.d0, spread(0.d0,1,18), -2.7d0]
  ! Best fit, quad = 10 results.  (Note effective neutral radii
  ! have to be reduced from their quad=0 results.)
  ! [-0.4d0, -0.2d0, spread(0.d0,1,18), -4.0d0]
  ! Best fit, quad = 20 results.  This comparison with Saumon went
  ! to higher density where we had to change coefficients such
  ! that H+ reduced and H2 reduced.
  ! [0.d0, -0.2d0, spread(0.d0,1,18), 2*-3.6d0]
  ! Revert back to 1.1.0.
       [-0.4_fp_kind, -0.2_fp_kind, spread(0._fp_kind,1,22),-4.0_fp_kind, -40._fp_kind]
  real(fp_kind), parameter :: pi_fit_neutral_original(nelements_pi_fitp2) = exp(pi_fit_neutral_ln_original/3._fp_kind)
  real(fp_kind), parameter :: pi_fit_neutral(nelements_pi_fitp2) = exp(pi_fit_neutral_ln/3._fp_kind)
  real(fp_kind), parameter :: pi_fit_neutral_saumon(nelements_pi_fitp2) = exp(pi_fit_neutral_ln_saumon/3._fp_kind)

  integer, parameter :: nions_pi_fitp2 = 318

  real(fp_kind), parameter :: pi_fit_ion_ln_original(nions_pi_fitp2) =&
  ! best fit to MDH results
       [0._fp_kind, 0._fp_kind, 0._fp_kind, spread(0._fp_kind,1,nions_pi_fitp2-5), 0._fp_kind, 0._fp_kind]
  ! best fit to saumon H and He table has negligible ion interactions
  real(fp_kind), parameter :: pi_fit_ion_ln_saumon(nions_pi_fitp2) =&
  [-40._fp_kind, -40._fp_kind, -40._fp_kind, spread(-40._fp_kind,1,nions_pi_fitp2-3)]

  ! best fit to opal results (H2 and H2+ constants not relevant)
  ! data pi_fit_ion_ln/-1.7d0, 0.d0, 0.2d0, 292*-1.d0, -50.d0, -50d0/
  ! attempt to get some rough agreement with other pressure ionization
  ! modes and also Saumon high density H results.
  ! data pi_fit_ion_ln/-2.13d0, 0.d0, 0.2d0, 292*-1.d0, -0.7d0, -1.d0/
  ! try opal fit as close to Saumon fit as possible
  ! (i.e. reduced H, He, low ions of metals, H2, H2+  ion interactions).
  ! further modification to high metal ions to force full metal ionization
  ! at high rho, T
  ! further modification to get monotonic metals for LMS model, change
  ! -1.8 to zero for first two ions
  real(fp_kind), parameter :: pi_fit_ion_ln(nions_pi_fitp2) = [&
  !Attempt to get smoothest LMS fit when unconstrained by eos2001
  ! &  -2.5d0, !H
  !OLD Attempt to solar fit eos2001
  ! &  -40.d0, !H
  !20041211 Attempt to solar fit eos2001 (start with quad=10 result)
  ! &  -1.8d0, !H
  ! back to overall eos2001 fit/
  ! &  -40.d0, !H
  ! &  -2.5d0, -2.5d0, !He
  !Attempt to get smoothest LMS fit when unconstrained by eos2001
  ! &  -1.5d0, -2.5d0, !He
  !Old Attempt to solar fit eos2001
  ! &  -40.d0, -5.d0, !He
  !20041211 Attempt to solar fit eos2001 (start with quad=10 result)
  ! &  -1.8d0, -1.8d0, !He
  ! &  -40.d0, -5.d0, !He
  ! &  0.d0, 0.d0,  4*6.d0, !C
  ! &  -1.5d0, -2.5d0,  4*3.d0, !C
  ! &  -1.5d0, -2.5d0,  5*3.d0, !N
  ! &  -1.5d0, -2.5d0,  6*3.d0, !O
  ! &  -1.5d0, -2.5d0,  8*3.d0, !Ne
  ! &  -1.5d0, -2.5d0,  9*3.d0, !Na
  ! &  -1.5d0, -2.5d0, 10*3.d0, !Mg
  ! &  -1.5d0, -2.5d0, 11*3.d0, !Al
  ! &  -1.5d0, -2.5d0, 12*3.d0, !Si
  ! &  -1.5d0, -2.5d0, 13*3.d0, !P
  ! &  -1.5d0, -2.5d0, 14*3.d0, !S
  ! &  -1.5d0, -2.5d0, 15*3.d0, !Cl
  ! &  -1.5d0, -2.5d0, 16*3.d0, !A
  ! &  -1.5d0, -2.5d0, 18*3.d0, !Ca
  ! &  -1.5d0, -2.5d0, 20*3.d0, !Ti
  ! &  -1.5d0, -2.5d0, 22*3.d0, !Cr
  ! &  -1.5d0, -2.5d0, 23*3.d0, !Mn
  ! &  -1.5d0, -2.5d0, 24*3.d0, !Fe
  ! &  -1.5d0, -2.5d0, 26*3.d0, !Ni
  ! &  2*-1.8d0/ !H2, H2+
  ! quad=20 result used by Cassisi student for thesis work
  ! &  -4.5d0, -4.5d0/ !H2, H2+
  ! quad=20 result adjusted for monotonic H charge fraction (H2+ + H+),
  ! suppression of H2+ above 10^5 K, and avoidance of grada sign change
  ! all for 0.075 model locus.
  ! &  -2.0d0, -2.5d0/ !H2, H2+
  ! and suppress H2 and H2+ some more now I have discontinuity index to worry
  ! about
  ! Best fit that minimizes discontinuity for pre eos2001
  ! &  -4.0d0, -4.5d0/ !H2, H2+
  ! eos_2001 fit.  Start with Saumon
  ! &  -40.d0, -40.d0/ !H2, H2+
  ! Revert back to 1.1.0
  !H
       -1.8_fp_kind, &
  !He
       -1.8_fp_kind, -1.8_fp_kind, &
  !C
       0._fp_kind, 0._fp_kind,  spread(6._fp_kind,1,4), &
  !N
       0._fp_kind, 0._fp_kind,  spread(6._fp_kind,1,5), &
  !O
       0._fp_kind, 0._fp_kind,  spread(6._fp_kind,1,6), &
  !Ne
       0._fp_kind, 0._fp_kind,  spread(6._fp_kind,1,8), &
  !Na
       0._fp_kind, 0._fp_kind,  spread(6._fp_kind,1,9), &
  !Mg
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,10), &
  !Al
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,11), &
  !Si
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,12), &
  !P
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,13), &
  !S
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,14), &
  !Cl
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,15), &
  !A
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,16), &
  !Ca
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,18), &
  !Ti
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,20), &
  !Cr
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,22), &
  !Mn
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,23), &
  !Fe
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,24), &
  !Ni
       0._fp_kind, 0._fp_kind, spread(6._fp_kind,1,26), &
       ! Li
       0._fp_kind,0._fp_kind,spread(6._fp_kind,1,1),&
       ! Be
       0._fp_kind,0._fp_kind,spread(6._fp_kind,1,2),&
       ! B
       0._fp_kind,0._fp_kind,spread(6._fp_kind,1,3),&
       ! F
       0._fp_kind,0._fp_kind,spread(6._fp_kind,1,7),&
  !H2, H2+
       -1.8_fp_kind, -1.8_fp_kind]

  real(fp_kind), parameter :: pi_fit_ion_original(nions_pi_fitp2) = exp(pi_fit_ion_ln_original/3._fp_kind)
  real(fp_kind), parameter :: pi_fit_ion(nions_pi_fitp2) = exp(pi_fit_ion_ln/3._fp_kind)
  real(fp_kind), parameter :: pi_fit_ion_saumon(nions_pi_fitp2) = exp(pi_fit_ion_ln_saumon/3._fp_kind)

  ! pressure-ionization data for excited states.
  real(fp_kind), parameter ::&
       pi_fitx_neutral_original = exp(0._fp_kind/3._fp_kind),&
       pi_fitx_ion_original = exp(0._fp_kind/3._fp_kind)

  ! pre eos2001 fit to original opal.
  ! data pi_fitx_neutral_ln, pi_fitx_ion_ln/0.d0,-1.8d0/
  ! old fit to eos2001, start with Saumon
  ! data pi_fitx_neutral_ln, pi_fitx_ion_ln/-4.d0,-40.d0/
  ! new fit to eos2001, start with solar considerations.  It turns out the
  ! second fudge factor is important (changes of 0.003 in ln P) for
  ! solar conditions near log T = 5.
  ! data pi_fitx_neutral_ln, pi_fitx_ion_ln/0.d0,-1.8d0/
  ! data pi_fitx_neutral_ln, pi_fitx_ion_ln/-4.d0,-40.d0/
  ! Revert back to 1.1.0.
  real(fp_kind), parameter ::&
       pi_fitx_neutral = exp(0._fp_kind/3._fp_kind),&
       pi_fitx_ion = exp(-1.8_fp_kind/3._fp_kind)
  ! values which fit Saumon H table
  ! N.B. Fit of Saumon He table not relevant because it was calculated with no excited states.
  real(fp_kind), parameter ::&
       pi_fitx_neutral_saumon = exp(0._fp_kind/3._fp_kind),&
       pi_fitx_ion_saumon = exp(-40._fp_kind/3._fp_kind)
  ! quad pressure-ionization data.
  real(fp_kind), parameter :: quad_original = 10._fp_kind
  ! all modified forms now use this term (adjusted for best
  ! possible fit and convergence I could find for highest density Saumon H table.)
  ! parameter(quad_modified = 20.d0)
  ! eos2001 fit (value doesn't matter for relatively low density fit done at
  ! this time [rho_lim = -0.5], but revisit later when want more stability
  ! at higher density.)
  ! parameter(quad_modified = 20.d0)
  ! Revert back to 1.1.0 values.
   real(fp_kind), parameter :: quad_modified = 10._fp_kind

contains
#include "effective_radius.f90"
end module mod_pi_fit
