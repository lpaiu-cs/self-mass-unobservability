!*******************************************************************************
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
! ifh2 = 0  no h2
! ifh2 = 1  vdb h2
! ifh2 = 2  st h2
! ifh2 = 3  irwin h2
! ifh2plus = 0  no h2plus
! ifh2plus = 1  st h2plus
! ifh2plus = 2  irwin h2plus
! morder_in = 3, 5, or 8  (or 13, 15, or 18) means 3rd, 5th, or 8th order
! fermi_dirac integral approximation following EFF fit (or modified
! version of EFF fit which reduces to Cody-Thacher approximation for
! low relativistic correction).
! morder_in = -3, -5, or -8 (or -13, -15, or -18) means use above
! approximations in non-relativistic limit.
! morder_in = 1 uses Cody-Thacher approximation directly.
! morder_in = 21 calculate Fermi-Dirac integrals with slow, but precise
! (~1.d-9 relative errors) numerical integration.
! morder_in = -21 is same as 21 in non-relativistic limit.
! morder_in = 23 means original 3rd order eff result
! morder_in = -23 is same as 23 in non-relativistic limit.
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
!  0  < ifexchange < 10 --> KLVH treatment (i.e., drop K term).
! 10 <= ifexchange < 20 --> Kapusta treatment (i.e., retain K term).
! 20 <= ifexchange < 30 --> special test case, use I term alone.
! 30 <= ifexchange < 40 --> special test case, use J term alone.
! 40 <= ifexchange < 50 --> special test case, use K term alone.
! N.B. if ifexchange is negative, then use non-relativistic limit
! of corresponding positive ifexchange option.
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

! n.b. important metals controlled by array iftracemetal below.
!   currently list includes C, N, O, Ne, Mg, Si, S, Fe
! ifmtrace = 0 treat important metals like H and He.
! ifmtrace = 1 treat everything but H, He like trace metals; e.g.,
!   partially ionized (lw.gt.3) or fully ionized (lw.le.3).
! ifcoulomb > 9, do diffraction correction
! ifcoulomb > 0 means do metal Coulomb contribution to sum0 and sum2 exactly
! -10 < ifcoulomb < 0 means do metal Coulomb contribution to sum0 and
!   sum2 using the "metal Coulomb" approximation.
!   n.b. this approximation is internally replaced by PTEH
!   approximation (next line) if either eps(1) or eps(2) are zero.
! ifcoulomb < -9 means use PTEH approximation for sum0 and sum2.
!   n.b. if this approximation is combined with the PTEH
!   approximation for pressure ionization, then the
!   routine is much faster because no ionization fraction
!   iterations are required.
! mod(|ifcoulomb|,10) = 0 ignore Coulomb interaction.
! mod(|ifcoulomb|,10) = 1 use Debye-Huckel Coulomb approximation.
! mod(|ifcoulomb|,10) = 2 use Debye-Huckel Coulomb approximation with tau(x) correction
! = 3 use PTEH Coulomb approximation with their theta_e
! mod(|ifcoulomb|,10) = 4 use PTEH Coulomb approximation with fermi-dirac theta_e
! mod(|ifcoulomb|,10) = 9 same as 4 with DeWitt definition of lambda
!   (using sum0a = sum0 + ne*theta_e).
! mod(|ifcoulomb|,10) = 5 use DH smoothly connected to modified OCP
!   with DeWitt definition of lambda.
! mod(|ifcoulomb|,10) = 6 same as 5 with alternative smooth connection
! mod(|ifcoulomb|,10) = 7 DH (Gamma < 1) or OCP using new DeWitt lambda.
! mod(|ifcoulomb|,10) = 8 same as 7 with theta_e = 0.
!
! ifpi contains meaning of two flags:
! ifpi > 0 means use Planck-Larkin occupation probability, otherwise not.
! remaining meaning in absolute value of ifpi
! |ifpi| = 0, use no pressure ionization
! |ifpi| = 1, use pteh pressure ionization
! |ifpi| = 2, use fjs pressure ionization
! |ifpi| = 4, use Saumon-like variation of MDH pressure ionization
! |ifpi| > 4 same as zero, i.e., use no pressure ionization.
! ifrad = 0, no radiation pressure,
!   input match_variable is consistent (ln P excluding radiation
!   pressure) for kif = 1.
! ifrad = 1, radiation pressure included,
!   input match_variable is consistent (ln P including radiation
!   pressure) for kif = 1.
! ifrad = 2, radiation pressure included,
!   input match_variable is ln(ptotal-prad) for kif = 1, but
!   all output quantities are calculated with radiation included.
!   this feature is used to reduce significance loss in regions which
!   are dominated by radiation pressure.
!   for kif = 2, the input match_variable is ln rho as per normal, but
!   the output pressure(1) is ln(ptotal-prad).
!   n.b. this latter case is only used for some tables which have
!   rho and T as the independent variable and all quantities including
!   pressure derivatives calculated with radiation pressure *except for*
!   the pressure itself.  Also note for this latter case
!   (ifrad = 2, kif =2) that pressure(1) = ln (ptotal-prad) will be
!   inconsistent with pressure(2) and pressure(3) which will be
!   partials of ln ptotal wrt ln rho and ln t.
! lw < 3 every element treated as fully ionized.
! lw = 3 trace metals treated as fully ionized.
! lw > 3 very slow option with all elements treated as partially ionized
!   using all stages of ionization.
! ifexcited > 0 means use excited states (must have Planck-Larkin or |ifpi| = 3 or 4).
!  0 < ifexcited < 10 means use approximation to explicit summation
! 10 < ifexcited < 20 means use explicit summation
! mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules)
!   and helium.
! mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
! mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
! nmax is maximum principal quantum number included in excited state
!   sum for explicit summation.  Also, same interpretation when approximate
!   summation used for pure Planck-Larkin case with special signal that
!   nmax > 300000 means use approximation for infinite Planck-Larkin sum
!   rather than approximation for finite Planck-Larkin sum up to nmax.
! ifreducedmass = 1 (use reduced mass in equilibrium constant)
! ifreducedmass = 0 (use electron mass in equilibrium constant)
! iftc is only currently meaningful with kif = 1 or 2.
! iftc = 1 (use thermodynamic consistency arguments to obtain appropriate
!   derivatives of entropy).
! iftc = 0 (use direct analytical derivatives of entropy.))

! ifmodified <= 0 (original form of pressure ionization for |ifpi| = 1,2,3,4
! ifmodified > 0 (modified form of pressure ionization)
! kif = 0, ln f and ln t are independent variables.
! kif = 1, ln p = match_variable and ln t are independent variables.
! kif = 2, ln rho = match_variable and ln t are independent variables.
! n.b. kif > 0 are much slower options because ln f must be iterated to
!   match match_variable
! eps(nelements_in) is an array of neps = nelements = 20 values of relative
!   abundance by weight divided by the appropriate atomic weight.
!   We use the atomic weight scale where (un-ionized) C(12)
!   has a weight of 12.00000000....  All weights are for the
!   un-ionized element. The eps value for an element
!   should be the sum of the individual isotopic eps values for
!   that element.  The eps array refers to the
!   elements in the following order:
!   H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
! match_variable is matched by iterative adjustment of fl when kif > 0
! fl = ln of EFF degeneracy parameter related to eta by equation below.
!   fl is starting value for iterative adjustment upon input
!   and iteratively adjusted for output when kif > 0
! tl = ln t
! output quantities:
! fm = partial fl(tl, match_variable)/ partial match_variable (kif > 0).
! ft = partial fl(tl, match_variable)/ partial tl  (kif > 0).
! n.b. when ifrad = 2, kif = 1, these derivatives are with respect to
! ln P_gas.  All other output derivatives for this combination of flags
! are taken w.r.t. ln pressure.
! t = temperature
! rho = density
! rlout = ln rho
! p = pressure (total if ifrad = 1 or 2, excluding radiation if ifrad = 0)
! pl = ln p
! cf = cp (- partial ln rho(T,P)/partial T)^{1/2}
! cp = specific heat at constant pressure
! sf = partial entropy(fl, tl)/partial fl
! st = partial entropy(fl, tl)/partial tl
! grada = the adiabatic temperature gradient
! rtp = - partial ln rho(T,P) wrt ln T.  n.b. note negative sign so
!   quantity should be positive for normal EOS.
! rmue = rho/mu_e, where mu_e is the mean molecular
!   weight per electron.
! fh2 = n(H+)/n(all H)
! fhe2 = n(He+)/n(all He)
! fhe3 = n(He++)/n(all He)
! xmu1 = 1/mu, where mu is the mean molecular
!   weight per particle.
! xmu3 = 1/mu_e, where mu is the mean molecular
!   weight per electron.
! eta = degeneracy parameter (Cox and Guili eta), related
! to EFF ln f = fl by
! wf = sqrt(1.d0 + f) = d eta/d fl
! eta = fl+2.d0*(wf-log(1.d0+wf))
! gamma1-gamma3 are the gammas as defined in Cox and Guili
! h2rat = 2*n(H2)/n(all H)
! h2plusrat = 2*n(H2+)/n(all H)
! lambda = Coulomb interaction parameter (two definitions depending on ifcoulomb).
! gamma_e = Coulomb diffraction parameter.
! iteration_count = the total number of ionization fraction loops completed
!   to attempt to determine the solution
! info = status code.  IMPORTANT: anything other than zero means
!   there was an abnormal ending to the FreeEOS calculation and no
!   returned quantities should be considered to be reliable.
! degeneracy, pressure, density, energy, enthalpy, and entropy
!   are all 3-vectors, with the second component being the derivative
!   of the first component wrt match_variable (except for the case where
!   ifrad = 2), and the third component being the derivative of the first
!   component wrt tl.
! definitions:
!   degeneracy(1) is EFF degeneracy parameter ln f defined above.
!   pressure(1) = ln pressure (except for kif=2, ifrad=2).
!   density(1) = ln density.
!   energy(1) = internal energy per unit mass.
!   enthalpy(1) = enthalpy per unit mass = energy(1) + p/rho.
!   entropy(1) = entropy per unit mass.

!> This free_eos_detailed subroutine calculates the EOS for the
!> variety of different free-energy models specified by the large
!> number of input option parameters.  Because of the combinatorial
!> complexity of those options, this subroutine should not be called
!> by external users since many of those combinations are not valid or
!> are untested.  Instead, the generic free_eos wrappers for this
!> subroutine should be called since those wrappers specify the
!> different free-energy models that are tested with FreeEOS with just
!> four parameters, ifoption, ifmodified, ifion, and kif.  See the
!> internals of free_eos_modern (the particular wrapper which directly
!> calls free_eos_detailed) for how those 4 parameters are translated
!> into the large number of free_eos_detailed options.
!>
!> \param[in] verbosity PARAMETERS NEED DOCUMENTATION
!>
subroutine free_eos_detailed(&
     verbosity, ifh2, ifh2plus, morder, ifexchange_in,&
     ifmtrace, ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,&
     ifreducedmass, iftc, ifmodified, kif,&
     eps, match_variable_in, fl, tl_in, fm, ft,&
     t, rho, rlout, p, pl, cf, cp, sf, st, grada, rtp,&
     rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
     gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2,&
     degeneracy, pressure, density, energy, enthalpy, entropy,&
     iteration_count, info)

  use mod_pi_fit, only: effective_radius
  use mod_ionization_data, only: nions, nion, bi, h2diss
  use mod_free_eos_constants, only: avogadro, c2, c_e, cd, cpe, cr, prad_const, clight, ln10
  use mod_aux_scale, only: aux_underflow, ln_aux_underflow
  use mod_diagnostics, only: pr_ratio, ppie_ratio, pc_ratio, pex_ratio
  use mod_coulomb, only: master_coulomb_end
  use mod_lapack, only: solve_linear_svd, dgesvx
  use mod_exchange, only: master_exchange, exchange_end
  use mod_pi, only: mdh_pi_end, fjs_pi_end, pteh_pi_end
  use mod_pl, only: pl_prepare
  use mod_eos_jacobian, only: eos_jacobian
  use mod_eos_bfgs, only: eos_bfgs
  use mod_info_data, only: info_offset_free_eos_detailed, info_offset_lapack
  use mod_excitation, only: excitation_pi_end
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real, if_taint_allocated_integer, taint_allocated_integer
  use, intrinsic :: ieee_arithmetic, only: ieee_support_nan, ieee_is_nan
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  integer, intent(in) ::&
       verbosity,&
       ifh2, ifh2plus, morder, ifexchange_in,&
       ifmtrace, ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,&
       ifreducedmass, iftc, ifmodified, kif

  integer, intent(out) :: iteration_count, info

  real(fp_kind), intent(in) :: eps(:), match_variable_in, tl_in

  ! Initial value of fl provided by the calling routine.  That value is updated
  ! within this routine if kif /= 0.
  real(fp_kind), intent(inout) :: fl

  real(fp_kind), intent(out) ::&
       fm, ft,&
       t, rho, rlout, p, pl, cf, cp, sf, st, grada, rtp,&
       rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
       gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2,&
       degeneracy(:), pressure(:), density(:), energy(:),&
       enthalpy(:), entropy(:)

  ! Internal variables

  ! Parameters

  ! To give some background, I implemented solve_linear_svd
  ! because of the enthusiasm of Numerical Recipes for the SVD
  ! technique, and that routine appears to work reasonably well
  ! (e.g., it appears to be well debugged) from the perspective of
  ! both free_eos convergence and also ../utils/solve_linear_test
  ! results.
  !
  ! Nevertheless, LINPACK developers (see that users' guide) are
  ! not nearly as enthusiastic about SVD as Numerical Recipes and
  ! that leadership is followed by lapack which, for example, does
  ! provide a SVD routine but does not provide practical methods
  ! of using those results similar to solve_linear_svd.
  !
  ! Furthermore, my experience is that solve_linear_svd typically
  ! does not converge quite as well as dgesvd for extreme density
  ! conditions (log rho (cgs) > 2) regardless of tolerance value.
  ! And I got similar results for ../utils/solve_linear_test.
  ! Therefore, my decision at least for now is to set if_svd to
  ! false, i.e., always use dgesvx rather than solve_linear_svd.
  logical, parameter :: if_svd = .false.
  ! All direct calls to eos_jacobian skip the call to free_non_ideal_calc in that routine.
  logical, parameter :: if_free_non_ideal_calc = .false.

  ! Three parameters required for hacked partial derivative tests.
  ! Normally, all of these are .false., but locally set one (and only one) of
  ! these to be .true. to do some hacked tests of the dv_aux, aux_dv, or
  ! jacobian partial derivatives internally calculated by
  ! free_eos_detailed or the routines it calls.
  logical, parameter :: debug_dv_aux = .false.
  logical, parameter :: debug_aux_dv = .false.
  logical, parameter :: debug_jacobian = .false.
  ! Useful Boolean combination of the above three.
  logical, parameter :: debug_any = debug_dv_aux.or.debug_aux_dv.or.debug_jacobian

  integer, parameter :: rl_aux_index = 4
  integer, parameter :: neps_local = 20
  ! maximum number of allowed ln f iterations
  integer, parameter :: ltermax = 200
  ! these two parameters must agree exactly with corresponding parameters
  ! in the ionize subroutine.
  integer, parameter :: nelements = 20
  !integer, parameter :: maxfjs_aux = 4
  ! maximum core charge for non-bare ion
  integer, parameter :: maxcorecharge = 28
  integer, parameter :: maxioncount=2000
  ! mdh-like treatment has a maximum of 9 extra parameters
  integer, parameter :: maxnextrasum = 9
  integer, parameter :: nxextrasum = 4
  ! iextraoff is the number of non-extrasum auxiliary variables
  integer, parameter :: iextraoff = 7
  ! naux is the total number of auxiliary variables.
  ! one extra in case degeneracy parameter is included as one
  ! of the auxiliary variables.
  integer, parameter :: naux = iextraoff + maxnextrasum + nxextrasum + 1
  integer, parameter :: nderivp1 = 3
  integer, parameter :: iatomic_number(nelements) = [&
       1,2,6,7,8,10,11,12,13,14,15,16,17,18,20,22,24,25,26,28]
  integer, parameter :: iftracemetal(nelements) = [&
       ! iftracemetal controls what elements are treated as fully ionized
       ! when lw = 3.
       ! H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni
       ! MDH dela residual ridge < 0.02
       ! 10,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
       ! MDH dela residual ridge < 0.01 if Fe not trace
       ! 10,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1/
       ! MDH dela residual ridge < 0.001 if Mg, Si, and S not trace
       0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1]

  ! This value is only used to decide when to output rcond as a warning
  ! that there are ill-conditioning problems in determining the NR solution
  ! of the EOS.  3 figures of accuracy is probably sufficient to
  ! obtain a good NR solution.  So don't remark on rcond unless you
  ! might have fewer than 3 figures of accuracy in the NR iteration.
  real(lapack_fp_kind), parameter :: rcond_lapack_min = 1.e-13_lapack_fp_kind
  ! This tolerance discussion is only relevant for the if_svd = .true. case.
  ! Throw out parts of SVD solution that don't satisfy this
  ! criterion.  See Press et al. Numerical Recipes Section 2.9
  ! discussion for why it is a good idea to throw out
  ! low-significance parts of the solutions of linear sets of
  ! equations.
  ! N.B. the resulting SVD condition number will be larger than
  ! this tolerance.
  ! Notes:
  ! tolerance = 1.d-7 gave somewhat worse results for the
  ! difficult convergence issues encountered by the plots for the
  ! convergence paper than the dgesvx results.  What were severe
  ! discontinuities before turned into divergences for a
  ! relatively limited set of cases that are far beyond stellar
  ! densities except possibly for the pure H case which diverged
  ! for log rho (cgs) ~ 2 - 2.5 at log T = 6 rather than lower log
  ! T values.
  ! Essentially, there are the same issues (including the log T = 6
  ! divergence cases) for tolerance = 1.d-12 and 1.d-14.
  real(lapack_fp_kind), parameter :: tolerance_lapack = 1.e-12_lapack_fp_kind
  ! Set up fast attack (immediate raise of simple lambda to
  ! simple_lambda_max if simple criteria are satisfied) and slow
  ! decay (multiply simple_lambda by simple_lambda_ratio whenever the
  ! criteria are not met).
  real(fp_kind), parameter :: simple_lambda_max = 1.e3_fp_kind
  ! relatively slow decay (6 iterations to zero simple_lambda)
  real(fp_kind), parameter :: simple_lambda_ratio = 1.e-1_fp_kind
  ! floor value below which simple_lambda is set to zero.
  real(fp_kind), parameter :: simple_lambda_min = 1.1e-2_fp_kind
  ! lapack variables
  real(fp_kind), parameter :: faux_limit_small_start = 2._fp_kind
  ! used 50 with dgeco, and 50 seems to be okay with lapack svd.
  ! 14 of these steps would be close to the over or under flow limit.
  real(fp_kind), parameter :: faux_limit_large = 50._fp_kind
  ! new control logic:
  ! limit used to help control when diag solution used.  That diag solution
  ! only used if the (log) change in auxiliary variable is greater than
  ! this value and diag solution smaller than NR solution.
  real(fp_kind), parameter :: faux_limit_diag = 0.5_fp_kind
  ! limit on smallest value of faux_limit_small (which is ordinarily
  ! decreased when maxfaux changes sign)
  ! minimum faux_limit
  real(fp_kind), parameter :: faux_limit_min = 1.e-10_fp_kind
  real(fp_kind), parameter :: major_crit = 1.e-3_fp_kind

  ! Allocatable arrays
  logical, allocatable ::&
       old_allow_log(:), old_allow_log_neg(:), old_new_allow_log(:), old_new_allow_log_neg(:), ifmajor(:)

  integer, allocatable ::&
       ifaux(:),&
       partial_aux(:),&
       inv_aux(:),&
       ipiv_lapack(:),&
       iwork_lapack(:)

  real(fp_kind), allocatable ::&
       dvzero(:),&
       dv(:),&
       dvf(:),&
       dvt(:),&
       plop(:),&
       plopt(:),&
       plopt2(:),&
       dv_pl(:),&
       dv_plt(:),&
       h2_dv(:),&
       aux_old(:),&
       aux_restore(:),&
       faux(:),&
       faux_nr(:),&
       zerolim(:),&
       dv_aux(:,:),&
       rhs1(:),&
       rhs1_save(:),&
       rhs2(:,:),&
       jacobian(:,:),&
       jacobian_save(:,:),&
       pnorad_aux(:),&
       simple_lambda(:),&
       free_aux(:)

  real(lapack_fp_kind), allocatable ::&
       lu_lapack(:,:),&
       row_lapack(:),&
       col_lapack(:),&
       sol_lapack(:,:),&
       ferr_lapack(:),&
       berr_lapack(:),&
       work_lapack(:),&
       rhs1_lapack(:),&
       rhs2_lapack(:,:),&
       jacobian_lapack(:,:)
  character equed_lapack

  logical&
       ifnoaux_iteration, ifwarm, any_ifsimple,&
       ifsame_under, ifdvzero, ifnuform
  logical ifnr03

  integer ifpi_fit, ifpi_fit_old
  data ifpi_fit_old/-1000/
  integer ifprintrtp, ifprintcp, ifprintgrada, ifprintsound2
  data ifprintrtp/1/  !one warning on negative rtp
  data ifprintcp/1/  !one warning on negative cp
  data ifprintgrada/1/  !one warning on negative grada
  data ifprintsound2/1/  !one warning on negative sound2
  integer mion_end, number_electrons
  integer iz, izlo, izhi, ion_start
  integer ifnr,&
       if_pteh, ifcoulomb_mod, if_mc, if_dc, i, lter,&
       ifionized, ioncount, ion, max_index,&
       index,&
       nextrasum, ifpi_local, ifpl,&
       ifsame_abundances, ifnear_abundances,&
       ifsame_zero_abundances,&
       iflast, lw_old, ifexcited_old, nmax_old,&
       ifh2_old, ifh2plus_old, morder_old, ifmtrace_old, ifcoulomb_old,&
       ifmtrace_olda, lw_olda, if_mc_olda, ifh2_olda, ifh2plus_olda
  data lw_old, ifexcited_old, nmax_old,&
       ifh2_old, ifh2plus_old, morder_old, ifmtrace_old, ifcoulomb_old,&
       ifmtrace_olda, lw_olda, if_mc_olda, ifh2_olda, ifh2plus_olda&
       /13*-1000/
  integer ion_stop, n_partial_ions, n_partial_elements
  integer n_partial_aux, njacobian,&
       index_aux, iaux, jndex_aux, jaux,&
       ielement, index_max,&
       ioncountzero, ifzerocount
  integer info_lapack
  integer iffirst
  data iffirst/1/
  integer ifmodified_old
  ! silly value
  data ifmodified_old/-1000/
  integer itemp1, itemp2, iaux1, iaux2
  integer jtemp1, jtemp2, jtemp3, jtemp4, jtemp5, jtemp6
  integer debug_count
  data debug_count/0/
  integer ijacobian, jjacobian
  integer isimple_bfgs
  integer neps

  ! Non-allocatable arrays due to data statements).
  real(fp_kind) eps_old(neps_local)
  ! n.b. MUST be zero to start to force ifsame_abundances and
  ! ifsame_zero_abundances to be 0 on first entry into routine.
  data eps_old/neps_local*0._fp_kind/
  real(fp_kind) dni(2), dne(2)
  data dni/0.75_fp_kind, 0.166667_fp_kind/
  data dne/0.5_fp_kind, 0.0_fp_kind/

  ! Non-allocatable arrays due to equivalence statements.
  real(fp_kind) rhostar(9), pstar(9), sstar(3), ustar(3)
  real(fp_kind) aux(naux), auxf(naux), auxt(naux)
  real(fp_kind)&
       extrasum(maxnextrasum), extrasumf(maxnextrasum), extrasumt(maxnextrasum),&
       xextrasum(nxextrasum), xextrasumf(nxextrasum), xextrasumt(nxextrasum)
  real(fp_kind) aux_dv(nions+2,naux)
  real(fp_kind)&
       h2plus_dv(nions+2),&
       h_ion_dv(nions+2), he_ion_dv(nions+2), he_ion2_dv(nions+2),&
       r_dv(nions+2), sum0_dv(nions+2), sum2_dv(nions+2),&
       extrasum_dv(nions+2, maxnextrasum),&
       xextrasum_dv(nions+2,nxextrasum)

  ! Non-allocatable arrays due to save statements.
  integer ion_end(nelements+2),&
       nmin(maxcorecharge), nmin_max(maxcorecharge),&
       nmin_species(nions+2),&
       ifdv(nions+2),&
       ifelement(nelements),&
       partial_ions(nions+2), partial_elements(nelements+2),&
       inv_ion(nions+2)

  real(fp_kind) charge(nions+2), charge2(nions+2),&
       r_ion3(nions+2), r_neutral(nelements+2),&
       bmin(maxcorecharge)

  real(fp_kind) match_variablef, match_variablet

  ! required for all debugging tests.
  ! dv_aux refers to the partial of dv wrt old auxiliary variables
  ! and fl and tl.
  ! aux_dv refers to the partial of the new auxiliary variables wrt dv
  ! and fl and tl.
  ! jacobian refers to the partial of the zeroed functions (either
  ! logarithmic form of aux(old) - aux(new) or calculated match_variable
  ! - input match_variable) wrt old auxiliary variables and fl (and ft).
  real(fp_kind) flminus, flplus, paaminus, paaplus
  real(fp_kind) ne, ni, nef, net
  real(fp_kind) re, ref, ret, pe, pef, pet, se, sef, set, ue, uef, uet, pe_cgs
  equivalence&
       (rhostar(1),re), (rhostar(2),ref), (rhostar(3),ret),&
       (pstar(1),pe), (pstar(2),pef), (pstar(3),pet),&
       (sstar(1),se), (sstar(2),sef), (sstar(3),set),&
       (ustar(1),ue), (ustar(2),uef), (ustar(3),uet)
  real(fp_kind) dve_exchange, dve_exchangef, dve_exchanget
  real(fp_kind)&
       full_sum0, full_sum1, full_sum2,&
       sum0, sum0_mc, sum0ne, sum0f, sum0t,&
       sum2, sum2_mc, sum2ne, sum2f, sum2t,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       pr, tc2
  real(fp_kind) paa, pab, pac_nr, pac, lnrho_5, f, wf, n_e,&
       hcon_mc, hecon_mc,&
       dpcoulomb, dpcoulombt, dpcoulombf,&
       dscoulomb, dscoulombf, dscoulombt, ducoulomb,&
       ne_old, lambda_old, s,&
       sion, sionf, siont, u, uion,&
       h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t, rt, rf,&
       h2, h2f, h2t,&
       h2plus, h2plusf, h2plust,&
       p0, pion, pionf, piont, pf, sr, pt, rpt,&
       vx, vz,&
       pex, pexf, pext, sex, sexf, sext, uex,&
       ppi, ppif, ppit, spi, spif, spit, upi,&
       pexcited, pexcitedf, pexcitedt,&
       sexcited, sexcitedf, sexcitedt, uexcited,&
       nux, nuy, nuz,&
       pnorad, pnoradf, pnoradt,&
       match_variableold, tlold, flold
  ! must be wild values.
  data match_variableold, tlold, flold/3*1.e30_fp_kind/
  ! set up equivalenced auxiliary variables that must be iterated to
  ! consistency.
  real(fp_kind) rl, rl_old
  real(fp_kind) fl_restore,&
       eps_aux, faux_limit,&
       faux_limit_small, faux_scale,&
       maxfaux, maxfaux_diag, maxfaux_diag_old

  ! this equivalence organization is *important* for later logic on
  ! deciding ifaux.  The order is h_ion, he_ion, he_ion2, rl, h2plus, sum0, sum2,
  ! extrasum(maxnextrasum), xextrasum(nxextrasum)
  equivalence&
       (aux(1),h_ion),(aux(2),he_ion),(aux(3),he_ion2),&
       (auxf(1),h_ionf),(auxf(2),he_ionf),(auxf(3),he_ion2f),&
       (auxt(1),h_iont),(auxt(2),he_iont),(auxt(3),he_ion2t),&
       (aux(4),rl),(aux(5),h2plus),&
       (auxf(4),rf),(auxf(5),h2plusf),&
       (auxt(4),rt),(auxt(5),h2plust),&
       (aux(6),sum0),(aux(7),sum2),&
       (auxf(6),sum0f),(auxf(7),sum2f),&
       (auxt(6),sum0t),(auxt(7),sum2t),&
       (aux(iextraoff+1),extrasum(1)),&
       (auxf(iextraoff+1),extrasumf(1)),&
       (auxt(iextraoff+1),extrasumt(1)),&
       (aux(iextraoff+maxnextrasum+1),xextrasum(1)),&
       (auxf(iextraoff+maxnextrasum+1),xextrasumf(1)),&
       (auxt(iextraoff+maxnextrasum+1),xextrasumt(1))
  ! equivalence aux_dv and h_ion_dv, he_ion_dv, he_ion2_dv, r_dv, h2plus, extrasum_dv
  ! note aux index of aux_dv cannot be indexed because these arrays
  ! used separately inside eos_calc so that aux position must be fixed.
  ! Note, however, that in all _dv arrays, the dv index *is* indexed.
  equivalence&
       (aux_dv(1,1), h_ion_dv),&
       (aux_dv(1,2), he_ion_dv),&
       (aux_dv(1,3), he_ion2_dv),&
       (aux_dv(1,4), r_dv),&
       (aux_dv(1,5), h2plus_dv),&
       (aux_dv(1,6), sum0_dv),&
       (aux_dv(1,7), sum2_dv),&
       (aux_dv(1,iextraoff+1), extrasum_dv),&
       (aux_dv(1,iextraoff+maxnextrasum+1), xextrasum_dv)
  !temporary
  ! variables associated with free energy calculation:
  real(fp_kind) free_rad, free_e, free_ion, free_pi, free_excited,&
       free_coulomb,&
       free_ex, free_pl
  real(fp_kind) free, free_fp
  real(fp_kind) match_variable, tl
  real(fp_kind) row_norm
  ! end of variables required for determining ifmajor

  ! These variables used to keep match_variable and tl changes
  ! required for hacked partial derivative tests independent of
  ! match_variable_in and tl_in.
  real(fp_kind) match_variable_save, tl_save, fl_save, match_variableold_save, flold_save, tlold_save
  ! These required to store change in independent variables of whatever
  ! hacked partial derivative tests might be run.
  real(fp_kind) delta1, delta2
  real(fp_kind), allocatable :: debug_results(:,:)

  real(lapack_fp_kind) rcond_lapack

  ! Most (if not all) compilers use save attribute for variables
  ! initialized with data statements, but just in case....
  save ifpi_fit_old, ifprintrtp, ifprintcp, ifprintgrada, ifprintsound2,&
       lw_old, ifexcited_old, nmax_old,&
       ifh2_old, ifh2plus_old, morder_old, ifmtrace_old, ifcoulomb_old,&
       ifmtrace_olda, lw_olda, if_mc_olda, ifh2_olda, ifh2plus_olda,&
       iffirst, ifmodified_old, debug_count, eps_old, dni, dne,&
       match_variableold, tlold, flold

  ! In debug mode these variables are saved every 33 calls.
  save match_variable_save, tl_save, fl_save, match_variableold_save, flold_save, tlold_save

  ! ifpi_fit_old-related values.
  save r_ion3, r_neutral

  ! iffirst-related values (except for bi which is saved in
  ! mod_ionization_data).
  save ion_end, charge, charge2, ion_start, nmin_species

  ! These values determined as the result of the following if statement
  ! if(ifsame_zero_abundances.ne.1.or.ifmtrace.ne.ifmtrace_old.or.&
  !     lw.ne.lw_old.or.ifexcited.ne.ifexcited_old.or.&
  !     ifh2.ne.ifh2_olda.or.ifh2plus.ne.ifh2plus_olda) then
  ! N.B. ignore ion_end which is saved above.
  save ion_stop, n_partial_ions, n_partial_elements, izlo, izhi, ifelement,&
       partial_elements, partial_ions, inv_ion, bmin, nmin, nmin_max,&
       max_index, mion_end, ifdv

  ! These values determined as the result of the following if statement
  !  if(ifsame_abundances.ne.1.or.ifmtrace.ne.ifmtrace_olda.or.&
  !     lw.ne.lw_olda.or.if_mc.ne.if_mc_olda) then
  save full_sum0, full_sum1, full_sum2, nux, nuy, nuz, sum0_mc, sum2_mc, hcon_mc, hecon_mc

  ! Logic below depends on saving prior results from aux, auxf, and auxt
  ! and the variables equivalenced to those arrays.
  save aux, auxf, auxt

  allocate(&
       old_allow_log(naux),&
       old_allow_log_neg(naux),&
       old_new_allow_log(naux),&
       old_new_allow_log_neg(naux),&
       ifmajor(naux),&
       ifaux(naux),&
       partial_aux(naux),&
       inv_aux(naux),&
       ipiv_lapack(naux),&
       iwork_lapack(naux),&
       dvzero(nelements),&
       dv(nions+2),&
       dvf(nions+2),&
       dvt(nions+2),&
       plop(nions+2),&
       plopt(nions+2),&
       plopt2(nions+2),&
       dv_pl(nions+2),&
       dv_plt(nions+2),&
       h2_dv(nions+2),&
       aux_old(naux),&
       aux_restore(naux),&
       faux(naux),&
       faux_nr(naux),&
       zerolim(naux),&
       dv_aux(nions+2,naux),&
       rhs1(naux),&
       rhs1_save(naux),&
       rhs1_lapack(naux),&
       rhs2(naux,2),&
       rhs2_lapack(naux,2),&
       jacobian(naux,naux),&
       jacobian_save(naux, naux),&
       jacobian_lapack(naux, naux),&
       pnorad_aux(naux+1),&
       simple_lambda(naux),&
       lu_lapack(naux,naux),&
       row_lapack(naux),&
       col_lapack(naux),&
       sol_lapack(naux,2),&
       ferr_lapack(2),&
       berr_lapack(2),&
       work_lapack(4*naux),&
       free_aux(naux+1))

  if(if_taint_allocated_real) then
     call taint_allocated_real(dvzero)
     call taint_allocated_real(dv)
     call taint_allocated_real(dvf)
     call taint_allocated_real(dvt)
     call taint_allocated_real(plop)
     call taint_allocated_real(plopt)
     call taint_allocated_real(plopt2)
     call taint_allocated_real(dv_pl)
     call taint_allocated_real(dv_plt)
     call taint_allocated_real(h2_dv)
     call taint_allocated_real(aux_old)
     call taint_allocated_real(aux_restore)
     call taint_allocated_real(faux)
     call taint_allocated_real(faux_nr)
     call taint_allocated_real(zerolim)
     call taint_allocated_real(dv_aux)
     call taint_allocated_real(rhs1)
     call taint_allocated_real(rhs1_save)
     call taint_allocated_real(rhs1_lapack)
     call taint_allocated_real(rhs2)
     call taint_allocated_real(rhs2_lapack)
     call taint_allocated_real(jacobian)
     call taint_allocated_real(jacobian_save)
     call taint_allocated_real(jacobian_lapack)
     call taint_allocated_real(pnorad_aux)
     call taint_allocated_real(simple_lambda)
     call taint_allocated_real(lu_lapack)
     call taint_allocated_real(row_lapack)
     call taint_allocated_real(col_lapack)
     call taint_allocated_real(sol_lapack)
     call taint_allocated_real(ferr_lapack)
     call taint_allocated_real(berr_lapack)
     call taint_allocated_real(work_lapack)
     call taint_allocated_real(free_aux)
  endif

  ! These are equivalanced static real arrays
  ! Fixme (2021).  This code has to be changed when the equivalence statements are
  ! eventually dropped in this routine in conformance with best practices.

  if(if_taint_allocated_real) then
     call taint_allocated_real(rhostar)
     call taint_allocated_real(pstar)
     call taint_allocated_real(sstar)
     call taint_allocated_real(ustar)
     ! aux, auxf, and auxt are saved so can only be tainted on first call.
     if(iffirst.eq.1) then
        call taint_allocated_real(aux)
        call taint_allocated_real(auxf)
        call taint_allocated_real(auxt)
     endif
     call taint_allocated_real(aux_dv)
  endif

  if(if_taint_allocated_integer) then
     call taint_allocated_integer(ifaux)
     call taint_allocated_integer(partial_aux)
     call taint_allocated_integer(inv_aux)
     call taint_allocated_integer(ipiv_lapack)
     call taint_allocated_integer(iwork_lapack)
  endif

  ! default good status
  info = 0

  ! Sanity checks
  if(kif.lt.0.or.kif.gt.2) error stop 'free_eos_detailed: kif outside valid range from 0 to 2'
  if(ifrad.eq.2.and.kif.ne.1) error stop 'free_eos_detailed: for ifrad == 2, kif values other than 1 are disabled.'
  ! approximation for excited states has derivative error and/or
  ! significance loss for combination of mhd and pl occupation
  ! probability.  Direct summation works quite well with nmax ~ 100
  ! when mhd occupation probabilities are active so drop excited
  ! approximation in all cases but pl-excitation (where it works
  ! very well and saves enormous amounts of time).
  if(ifexcited.gt.0.and.ifexcited.lt.10.and.(&
       abs(ifpi).eq.3.or.abs(ifpi).eq.4))&
       error stop 'free_eos_detailed: ifexcited approximation not debugged'
  if(ifexcited.gt.0.and..not.&
       (abs(ifpi).eq.3.or.abs(ifpi).eq.4.or.ifpi.gt.0))&
       error stop 'free_eos_detailed: bad ifexcited and ifpi combination'
  if(&
       nderivp1.ne.size(degeneracy).or.&
       nderivp1.ne.size(pressure).or.&
       nderivp1.ne.size(density).or.&
       nderivp1.ne.size(energy).or.&
       nderivp1.ne.size(enthalpy).or.&
       nderivp1.ne.size(entropy))&
       error stop 'free_eos_detailed: inconsistent nderivp1 array sizes'

  ! NR auxiliary variable iteration.  Go one more step than
  ! this criterion so should guarantee near machine precision
  ! for the expected (and often tested) quadratic convergence.
  eps_aux = 1.e-7_fp_kind

  neps = size(eps)
  if(neps.ne.nelements.or.neps.ne.neps_local) error stop 'free_eos_detailed: bad size of eps'

  ! Insulate local debug-only changes of match_variable and tl (if those occur) from input values.

  match_variable = match_variable_in
  tl = tl_in

  if(debug_any) then
     ! Save the original independent variable values used by
     ! free_eos_test logic.  The 33 is a hack (!) that accounts for
     ! free_eos_test making 1 call to free_eos using the original
     ! independent variables plus an *assumed* (check
     ! free_eos_test.stdin is consistent with this) 8 stanzas of +/-
     ! delta1 and +/- delta2 calls to free_eos for a total of 33 calls
     ! to free_eos (and free_eos_detailed) per grid point.
     if(mod(debug_count,33).eq.0) then
        match_variable_save = match_variable
        tl_save = tl
        ! This is the initial (quite wild) value of fl supplied by the calling routine.
        ! But save it so that fl is determined to be the same for each debug call of
        ! free_eos_detailed.
        fl_save = fl
        ! Save these so iterative criteria are the same (in case that makes any difference).
        match_variableold_save = match_variableold
        flold_save = flold
        tlold_save = tlold
     endif
     ! Calculate changes in independent variables used for hacked partial derivative
     ! tests.  Note these changes are logarithmic.
     delta1 = match_variable - match_variable_save
     delta2 = tl - tl_save
     ! for hacked partial derivative tests, match_variable and tl must be the saved
     ! values to insure numerical differences only occur for internal independent variables
     ! for fixed match_variable and tl.
     match_variable = match_variable_save
     tl = tl_save
     ! In addition, fl must be the saved value so that the iterative
     ! adjustment of fl that occurs via the call to eos_cold_start
     ! that is mandated by any debug case always generates the
     ! identical value of fl.
     fl = fl_save
     ! These must be saved values (in case that change in the iterative critera makes any difference).
     match_variableold = match_variableold_save
     flold = flold_save
     tlold = tlold_save
     debug_count = debug_count + 1
  endif
  if(iffirst.eq.1) then
     iffirst = 0
     if(verbosity.ge.2) write(stderr,'(a, i5,a,i5)')&
          "configured fp_kind range, precision = ", range(0._fp_kind), ",", precision(0._fp_kind)
     ion_end(1) = iatomic_number(1)
     do ielement = 2, nelements
        ion_end(ielement) = ion_end(ielement-1) + iatomic_number(ielement)
     enddo
     charge(1:nions+2) = real(nion(1:nions+2),fp_kind)
     charge2(1:nions+2) = charge(1:nions+2)*charge(1:nions+2)
     ion_start = 1
     do ielement = 1, nelements
        ! for each element, and hydrogen molecules calculate nmin
        ! nmin is 2 for 1 - 2-electron systems,
        ! nmin is 3 for 3 - 10-electron systems,
        ! nmin is 4 for 11 - 28-electron systems
        do ion = ion_start, ion_end(ielement)
           number_electrons = ion_end(ielement) - ion + 1
           if(number_electrons.le.2) then
              nmin_species(ion) = 2
           elseif(number_electrons.le.10) then
              nmin_species(ion) = 3
           else
              nmin_species(ion) = 4
           endif
        enddo
        ion_start = ion_end(ielement) + 1
     enddo
     ! H2 is a 2-electron system
     nmin_species(nions+1) = 2
     ! H2+ is a 1-electron system
     nmin_species(nions+2) = 2
  endif
  ! sort out what fitting factors will be applied to pressure ionization.
  if(abs(ifpi).eq.4) then
     ! factors to fit Saumon results.
     ifpi_fit = 2
  elseif(ifmodified.gt.0) then
     ! factors to fit opal results.
     ifpi_fit = 1
  else
     ! unity factors to mimic mdh results as closely as possible.
     ifpi_fit = 0
  endif
  if(ifpi_fit.ne.ifpi_fit_old.and.(abs(ifpi).eq.3.or.abs(ifpi).eq.4)) then
     ifpi_fit_old = ifpi_fit
     ! calculate effective MDH radii for anything-ion interactions
     ! (r_ion) and neutral-neutral (r_neutral) interactions.
     ! n.b. "neutral" here refers to neutral species and H2+, the
     ! only species considered to have non-zero radii in MHD model.
     call effective_radius(ifpi_fit, bi, nion, r_ion3, r_neutral)
  endif
  ! above 10^5 K, hydrogen is ionized and helium is usually at least
  ! first ionized.  In these circumstances it is best to use the
  ! bare ion as the dv pressure ionization zero rather than the neutral
  ! to avoid significance loss.
  if(tl.le.log(1.e5_fp_kind)) then
     ifdvzero = .true.
  else
     ifdvzero = .false.
  endif
  if(ifcoulomb.gt.9) then
     if_dc = 1
     error stop 'free_eos_detailed: diffraction correction disabled'
  else
     if_dc = 0
  endif
  if(ifpi.gt.0) then
     ifpl = 1
  else
     ifpl = 0
  endif
  ifpi_local = abs(ifpi)
  if(ifpi_local.gt.4) ifpi_local = 0
  if(lw.lt.3) then
     ! full ionization approximation...
     ! use pteh (e.g., full ionization) approximation for Coulomb sums.
     if_pteh = 1
     if_mc = 0
     ! use no pressure ionization, since all forms of pressure
     ! ionization have zero effect on pressure and entropy for full
     ! ionization
     ifpi_local = 0
     ! Planck-Larkin reduces to zero for full ionization.
     ifpl = 0
  elseif(ifcoulomb.ge.0) then
     if_mc = 0
     if_pteh = 0
  elseif(ifcoulomb.gt.-10.and.((eps(1).gt.0._fp_kind.and.eps(2).gt.0._fp_kind).or.(lw.eq.3.and.ifmtrace.eq.1))) then
     ! non-zero H, and He or metals fully ionized
     if_mc = 1
     if_pteh = 0
  else
     if_mc = 0
     if_pteh = 1
  endif
  ifnoaux_iteration = (ifcoulomb.eq.0.or.if_pteh.eq.1).and.ifpi_local.le.1
  ifcoulomb_mod = mod(abs(ifcoulomb),10)
  if(ifpi_local.eq.3.or.ifpi_local.eq.4) then
     nextrasum = maxnextrasum
  else
     nextrasum = 0
  endif
  ! each entry to free_eos_detailed could potentially have a different
  ! set of abundances or a change from zero to non-zero abundance.
  ! check this.
  ifsame_abundances = 1
  ifnear_abundances = 1
  ifsame_zero_abundances = 1
  do i = 1,neps
     if(eps(i).ne.eps_old(i)) then
        if(eps(i).eq.0._fp_kind.or.eps_old(i).eq.0._fp_kind) then
           ifsame_zero_abundances = 0
           ifnear_abundances = 0
        elseif(abs(eps(i)-eps_old(i)).gt.0.05_fp_kind*abs(eps(i))) then
           ifnear_abundances = 0
        endif
        ifsame_abundances = 0
        eps_old(i) = eps(i)
     endif
  enddo
  if(ifsame_zero_abundances.ne.1.or.ifmtrace.ne.ifmtrace_old.or.&
       lw.ne.lw_old.or.ifexcited.ne.ifexcited_old.or.&
       ifh2.ne.ifh2_olda.or.ifh2plus.ne.ifh2plus_olda) then
     ifh2_olda = ifh2
     ifh2plus_olda = ifh2plus
     ion_stop = 0
     n_partial_ions = 0
     n_partial_elements = 0
     ! requires ridiculous values.
     izlo = 1000
     izhi = 0
     do i = 1,nelements
        if(eps(i).gt.0._fp_kind.and.((iftracemetal(i).ne.1.and.ifmtrace.ne.1).or.i.le.2.or.lw.ge.4)) then
           ! accumulate data for non-zero eps if any of following conditions are true:
           ! 1) element not treated as trace metal.
           ! 2) element is H or He.
           ! 3) all elements (including trace metals) treated as partially ionized
           ifelement(i) = 1
           n_partial_elements = n_partial_elements + 1
           partial_elements(n_partial_elements) = i
           do ion = ion_stop+1, ion_stop+iatomic_number(i)
              ! ifdv(ion) = 1
              n_partial_ions = n_partial_ions + 1
              partial_ions(n_partial_ions) = ion
              inv_ion(ion) = n_partial_ions
              if(i.le.2.or.mod(ifexcited,10).eq.3) then
                 iz = ion - ion_stop
                 if(iz.lt.izlo.or.iz.gt.izhi) then
                    ! new iz value
                    izlo = min(izlo, iz)
                    izhi = max(izhi, iz)
                    if(iz.gt.maxcorecharge) error stop 'free_eos_detailed: internal logic failure 1'
                    bmin(iz) = bi(ion)
                    nmin(iz) = nmin_species(ion)
                    nmin_max(iz) = nmin_species(ion)
                 else
                    bmin(iz) = min(bmin(iz), bi(ion))
                    nmin(iz) = min(nmin(iz), nmin_species(ion))
                    nmin_max(iz) = max(nmin_max(iz), nmin_species(ion))
                 endif
              endif
           enddo
        else
           ifelement(i) = 0
           do ion = ion_stop+1, ion_stop+iatomic_number(i)
              ! ifdv(ion) = 0
              inv_ion(ion) = 0
           enddo
        endif
        ion_stop = ion_stop + iatomic_number(i)
     enddo
     if(n_partial_elements.eq.0)&
          error stop 'free_eos_detailed: ionized pure metals not implemented'
     max_index = n_partial_ions
     mion_end = partial_elements(n_partial_elements)
     ! set ifdv for H2 and H2+
     if(eps(1).gt.0._fp_kind) then
        if(ifh2.gt.0) then
           mion_end = nelements + 1
           ion_end(mion_end) = ion_end(mion_end-1) + 1
           if(ion_end(mion_end).ne.nions+1)&
                error stop 'free_eos_detailed: internal logic failure 2'
           ifdv(nions+1) = 1
           partial_elements(n_partial_elements+1) = nelements + 1
           max_index = max_index + 1
           partial_ions(n_partial_ions+1) = nions+1
           inv_ion(nions+1) = n_partial_ions+1
           ! H2+ required for excited H2.
           if(ifh2plus.gt.0.and.mod(ifexcited,10).gt.1) then
              iz = 1
              ion = nions + 1
              if(iz.lt.izlo.or.iz.gt.izhi) then
                 ! new iz value
                 izlo = min(izlo, iz)
                 izhi = max(izhi, iz)
                 if(iz.gt.maxcorecharge)&
                      error stop 'free_eos_detailed: internal logic failure 3'
                 bmin(iz) = bi(ion)
                 nmin(iz) = nmin_species(ion)
                 nmin_max(iz) = nmin_species(ion)
              else
                 bmin(iz) = min(bmin(iz), bi(ion))
                 nmin(iz) = min(nmin(iz), nmin_species(ion))
                 nmin_max(iz) = max(nmin_max(iz), nmin_species(ion))
              endif
           endif
        else
           ifdv(nions+1) = 0
           inv_ion(nions+1) = 0
        endif
        if(ifh2plus.gt.0) then
           mion_end = nelements + 2
           ion_end(mion_end) = ion_end(mion_end-1) + 1
           if(ion_end(mion_end).ne.nions+2) error stop 'free_eos_detailed: internal logic failure 4'
           ifdv(nions+2) = 1
           partial_elements(n_partial_elements+2) = nelements + 2
           max_index = max_index + 1
           partial_ions(n_partial_ions+2) = nions+2
           inv_ion(nions+2) = n_partial_ions+2
           if(mod(ifexcited,10).gt.1) then
              iz = 2
              ion = nions + 2
              if(iz.lt.izlo.or.iz.gt.izhi) then
                 ! new iz value
                 izlo = min(izlo, iz)
                 izhi = max(izhi, iz)
                 if(iz.gt.maxcorecharge) error stop 'free_eos_detailed: internal logic failure 5'
                 bmin(iz) = bi(ion)
                 nmin(iz) = nmin_species(ion)
                 nmin_max(iz) = nmin_species(ion)
              else
                 bmin(iz) = min(bmin(iz), bi(ion))
                 nmin(iz) = min(nmin(iz), nmin_species(ion))
                 nmin_max(iz) = max(nmin_max(iz), nmin_species(ion))
              endif
           endif
        else
           ifdv(nions+2) = 0
           inv_ion(nions+2) = 0
        endif
     else
        ifdv(nions+1) = 0
        ifdv(nions+2) = 0
        inv_ion(nions+1) = 0
        inv_ion(nions+2) = 0
     endif
  endif
  ! n.b. the conditions on ifmtrace and lw arise because partial_elements
  ! may change due to these flags.
  if(ifsame_abundances.ne.1.or.ifmtrace.ne.ifmtrace_olda.or.lw.ne.lw_olda.or.if_mc.ne.if_mc_olda) then
     ! these are distingushed from other "old" flags with "a" suffix
     ifmtrace_olda = ifmtrace
     lw_olda = lw
     if_mc_olda = if_mc
     ! calculate total number of electrons available to be freed.
     ! n.b. done once each call in case of abundance changes.
     ! full_sum0 is full ionization approximation to sum0/rho*NA,
     ! where sum0 is the sum over positive ion number densities,
     ! rho is the density, and NA is Avogadro's number.
     full_sum0 = 0._fp_kind
     ! full_sum2 is full ionization approximation to (sum2-ne*thetae)/rho*NA
     ! where sum2-ne*thetae = sum over positive ion number densities times
     ! the charge squared on those ions.
     full_sum2 = 0._fp_kind
     ! nux, nuy, nuz are number/volume of maximum possible ionization
     ! electrons (divided by rho*NA) for hydrogen, helium, and metals.
     nux = eps(1)
     nuy = 2._fp_kind*eps(2)
     nuz = 0._fp_kind
     do i = 1,nelements
        if(eps(i).gt.0._fp_kind) then
           full_sum0 = full_sum0 + eps(i)
           full_sum2 = full_sum2 + eps(i)*&
                real((iatomic_number(i)*iatomic_number(i)),fp_kind)
           if(i.gt.2) nuz = nuz + eps(i)*real((iatomic_number(i)),fp_kind)
        endif
     enddo
     ! define to ridiculous values as kludgey check that these values
     ! never used unless if_mc is 1 (in which case, good values defined).
     sum0_mc = huge(1._fp_kind)
     sum2_mc = huge(1._fp_kind)
     hcon_mc = huge(1._fp_kind)
     hecon_mc = huge(1._fp_kind)
     if(if_mc.eq.1) then
        sum0_mc = 0._fp_kind
        sum2_mc = 0._fp_kind
        ! Account for trace metals that are assumed fully ionized and which
        ! are ignored in the hcon_mc and hecon_mc calculations.
        do ielement = 3, nelements
           if(eps(ielement).gt.0._fp_kind.and.(iftracemetal(ielement).eq.1.or.ifmtrace.eq.1)) then
              sum0_mc = sum0_mc + eps(ielement)
              sum2_mc = sum2_mc + eps(ielement)*real((iatomic_number(ielement)*iatomic_number(ielement)),fp_kind)
           endif
        enddo
        ! n.b. other flags can supersede the hcon_mc and hecon_mc calculation.  Only accumulate those
        ! variables for metals that are allowed to be
        ! partially ionized.  Coulomb effect of fully ionized trace metals
        ! accounted for with sum0_mc and sum2_mc (see above) and Coulomb effect
        ! of completely ionized approximation calculated with separate code
        ! (in which case hcon_mc and hecon_mc will be zero).
        hcon_mc = 0._fp_kind
        hecon_mc = 0._fp_kind
        do index = 1, n_partial_elements
           ielement = partial_elements(index)
           if(ielement.gt.2) then
              hcon_mc = hcon_mc + eps(ielement)
              hecon_mc = hecon_mc + eps(ielement)*real((iatomic_number(ielement)*iatomic_number(ielement)),fp_kind)
           endif
        enddo
        hecon_mc = hecon_mc - hcon_mc
        if(eps(2).gt.0._fp_kind) then
           hecon_mc = hecon_mc/eps(2)
        elseif(hecon_mc.ne.0._fp_kind) then
           error stop 'free_eos_detailed: hecon_mc logic screwup'
        endif
        if(eps(1).gt.0._fp_kind) then
           hcon_mc = hcon_mc/eps(1)
        elseif(hcon_mc.ne.0._fp_kind) then
           error stop 'free_eos_detailed: hcon_mc logic screwup'
        endif
     endif
     ! full_sum1 is number/volume of maximum possible ionization electrons
     ! divided by rho*NA.
     full_sum1 = nux + nuy + nuz
  endif
  ! set up auxiliary variable logic
  ! h_ion, he_ion, he_ion2
  if(if_mc.eq.1.or.ifpi_local.eq.2) then
     if(eps(1).gt.0._fp_kind) then
        ifaux(1) = 1
     else
        ifaux(1) = 0
     endif
     if(eps(2).gt.0._fp_kind) then
        ifaux(2) = 1
        ifaux(3) = 1
     else
        ifaux(2) = 0
        ifaux(3) = 0
     endif
     ! this somewhat redundant when kif=2 since an auxiliary variable is
     ! equal to the match_variable, but it works without any special
     ! programming so I won't try to reduce number of auxiliary variables
     ! in this special case.
     ifaux(4) = 1
  else
     ifaux(1) = 0
     ifaux(2) = 0
     ifaux(3) = 0
     ifaux(4) = 0
  endif
  if(if_mc.eq.1.and.ifh2plus.gt.0) then
     ! need h2plus for if_mc approximation for Coulomb sums.
     ifaux(5) = 1
  else
     ifaux(5) = 0
  endif
  ! sum0, sum2
  if(if_pteh.eq.1.or.if_mc.eq.1) then
     ifaux(6) = 0
     ifaux(7) = 0
  else
     ifaux(6) = 1
     ifaux(7) = 1
  endif
  if(ifpi_local.eq.3.or.ifpi_local.eq.4) then
     ! Note that nextrasum = maxnextrasum for ifpi_local.eq.3.or.ifpi_local.eq.4.
     ifaux(iextraoff+1:iextraoff+maxnextrasum) = 1
  else
     ifaux(iextraoff+1:iextraoff+maxnextrasum) = 0
  endif
  ! xextrasum
  if((ifpi_local.eq.3.or.ifpi_local.eq.4).and.ifexcited.gt.0) then
     ifaux(iextraoff+maxnextrasum+1:iextraoff+maxnextrasum+4) = 1
  else
     ifaux(iextraoff+maxnextrasum+1:iextraoff+maxnextrasum+4) = 0
  endif
  ! fl is nominally the nauxth auxiliary variable, but it is treated
  ! in a different manner than the other auxiliary variables.
  ! N.B. because ifaux(naux).eq.0, then n_partial_aux.le.naux-1 and
  ! n_partial_aux+2.le.naux+1, see pnorad_aux and free_aux
  ! dimensions above.
  ifaux(naux) = 0
  n_partial_aux = 0
  do iaux = 1,naux
     if(ifaux(iaux).eq.1) then
        n_partial_aux = n_partial_aux + 1
        partial_aux(n_partial_aux) = iaux
        inv_aux(iaux) = n_partial_aux
     else
        inv_aux(iaux) = 0
     endif
  enddo
  if(kif.eq.0) then
     njacobian = n_partial_aux
  else
     njacobian = n_partial_aux + 1
  endif
  t = exp(tl)
  tc2 = c2/t
  ! calculate planck-larkin occupation probabilities and equilibrium
  ! constant changes.
  if(ifpl.eq.1)&
       call pl_prepare(partial_ions(:n_partial_ions+2), max_index,&
       nion, ifdv,&
       bi, tc2, plop, plopt, plopt2, dv_pl, dv_plt)
  if(ifrad.ge.1) then
     pr = prad_const*t**4
  else
     pr = 0._fp_kind
  endif
  ! criterion for doing warm start with previous auxiliary
  ! variables.
  ! At this point flold is either a deliberately wild initial
  ! value or else the fl value determined by the last call to
  ! free_eos_detailed.
  ! Cannot get into much trouble if both fl and flold < -10 so
  ! relax delta fl criterion in that case.

  if(ifnear_abundances.eq.1.and.abs(tlold-tl).le.0.05001_fp_kind.and.&
       (abs(flold-fl).le.0.50001_fp_kind.or.(max(flold,fl).lt.-10._fp_kind.and.abs(flold-fl).le.20._fp_kind)).and.&
       lw.eq.lw_old.and.&
       ifexcited.eq.ifexcited_old.and.&
       nmax.eq.nmax_old.and.&
       ifh2.eq.ifh2_old.and.&
       ifh2plus.eq.ifh2plus_old.and.&
       morder.eq.morder_old.and.&
       ifmtrace.eq.ifmtrace_old.and.&
       ifcoulomb.eq.ifcoulomb_old.and.&
       ifmodified.eq.ifmodified_old) then
     ifwarm = .true.
  else
     ifwarm = .false.
  endif
  if(debug_any) then
     ! always use cold start when testing derivatives
     ifwarm = .false.
  endif
  if(verbosity.ge.4.and..not.ifwarm) then
     write(stderr,*) 'ifwarm is .false.'
     write(stderr,*) 'tlold, tl, abs(tlold-tl) = ', tlold, tl, abs(tlold-tl)
     write(stderr,*) 'flold, fl, abs(flold-fl) = ', flold, fl, abs(flold-fl)
     write(stderr,*) 'lw_old, lw = ', lw_old, lw
     write(stderr,*) 'ifexcited_old, ifexcited = ', ifexcited_old, ifexcited
     write(stderr,*) 'nmax_old, nmax = ', nmax_old, nmax
     write(stderr,*) 'ifh2_old, ifh2 = ', ifh2_old, ifh2
     write(stderr,*) 'ifh2plus_old, ifh2plus = ', ifh2plus_old, ifh2plus
     write(stderr,*) 'morder_old, morder = ', morder_old, morder
     write(stderr,*) 'ifmtrace_old, ifmtrace = ', ifmtrace_old, ifmtrace
     write(stderr,*) 'ifcoulomb_old, ifcoulomb = ', ifcoulomb_old, ifcoulomb
     write(stderr,*) 'ifmodified_old, ifmodified = ', ifmodified_old, ifmodified
  endif
  ! Save flags for the next time the above test is done.
  ! flold and tlold saved later.
  lw_old = lw
  ifexcited_old = ifexcited
  nmax_old = nmax
  ifh2_old = ifh2
  ifh2plus_old = ifh2plus
  morder_old = morder
  ifmtrace_old = ifmtrace
  ifcoulomb_old = ifcoulomb
  ifmodified_old = ifmodified
  match_variableold = match_variable
  iteration_count = 0
  if(ifnoaux_iteration.or..not.ifwarm) then
     ! if no auxiliary variable iteration required or cold start
     ! then iterate on fl to make it consistent with match_variable.
     ! Initialization for loop convergence criteria
     paaplus = 1.e30_fp_kind
     paaminus = -1.e30_fp_kind
     ! mark undefined by ridiculous values
     flplus = -1.e30_fp_kind
     flminus = 1.e30_fp_kind
     lter = 0
     iflast = 0
     paa = 1._fp_kind  !assure at least twice (once if kif=0) through loop.
     pab = 1._fp_kind
     do while(iflast.ne.1)
        lter = lter + 1
        ! Newton-Raphson iteration is quadratic, so obtain
        ! machine precision (within significance loss noise of say
        ! 1.d-14) if do one more iteration after 10^-7 convergence.
        if(kif.eq.0.or.lter.ge.ltermax.or.abs(fl-flold).le.1.e-7_fp_kind) iflast = 1
        ! find themodynamically consistent set of fermi-dirac integrals and
        ! put the values and their derivatives into rhostar, pstar,
        ! sstar, and ustar which are equivalenced to re, pe, se, and ue
        ! and their f and t derivatives.  Also calculate exchange (when
        ! ifexchange_in > 0) and its effects on dv.
        call master_exchange(verbosity, fl, tl,&
             rhostar, pstar, sstar, ustar, morder, ifexchange_in,&
             dve_exchange, dve_exchangef, dve_exchanget)
        f = exp(fl)
        wf = sqrt(1._fp_kind + f)
        eta = fl+2._fp_kind*(wf-log(1._fp_kind+wf))
        ! number density of free electrons and electron pressure
        n_e = c_e*re
        ! rho/mu_e = n_e H = cd*re
        rmue = cd*re
        ! set ifionized.
        if(lw.gt.3) then
           ! partial ionization with all stages of ionization
           ! of all elements.
           ifionized = 0
        elseif(lw.eq.3) then
           ! full ionization of all trace metals.
           ifionized = 1
        else
           ! full ionization of all elements.
           ifionized = 2
        endif
        ! All calls of eos_cold_start should calculate nu = n/(rho*avogadro)
        ! form of h2, h2plus, h_ion, he_ion, he_ion2, and extrasum for the
        ! ifnoaux_iteration case or kif=0 case.  This allows the pressure
        ! and its derivatives to be calculated properly during the
        ! eos_cold_start iteration on fl for kif=1, affects nothing
        ! else that is used during the iterations, and gives the right
        ! result for entropy, etc., after the iterations are completed for
        ! the ifnoaux_iteration case.  (For the
        ! .not.ifnoaux_iteration.and.kif.eq.1 case, a final call to
        ! eos_cold_start is required with ifnuform = .false., see below).
        if(ifnoaux_iteration.or.kif.eq.1) then
           ifnuform = .true.
        else
           ifnuform = .false.
        endif
        iteration_count = iteration_count + 1
        ! N.B. The actual eos_cold_start EOS calculation is done adopting
        ! no or PTEH pressure ionization and adopting the PTEH approximation
        ! for the Coulomb sums.  However, the sums resulting from that
        ! EOS calculation are calculated depending on the values of
        ! if_pteh, if_mc, and ifpi_local which depend on the overall
        ! free-energy model, not the special model used for the eos_cold_start
        ! calculation.
        ! if near fl convergence, make sure to fix underflow criterion
        ! in ionize called by eos_calc.
        ! note paa = 1. (and pab = 1.) on first time through so always false
        ! on first time
        if(kif.gt.0.and.max(abs(paa),abs(pab)).le.1.e-2_fp_kind) then
           ifsame_under = .true.
        else
           ifsame_under = .false.
        endif
        call eos_cold_start(&
             verbosity, ifsame_under, lambda, gamma_e,&
             partial_ions, f, eta, wf, t, n_e, pstar,&
             dve_exchange, dve_exchangef, dve_exchanget,&
             full_sum0, full_sum1, full_sum2, charge, dv_pl, dv_plt,&
             ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
             ifsame_zero_abundances, ifexcited, ifnuform,&
             inv_aux,&
             inv_ion, max_index,&
             partial_elements(1:n_partial_elements+2), ion_end,&
             ifionized, if_pteh, if_mc, ifreducedmass,&
             ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,&
             ifpl, ifmodified, ifh2, ifh2plus,&
             izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
             eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
             r_ion3, r_neutral,&
             ifelement, dvzero, dv, dvf, dvt, nion,&
             rhostar,&
             nextrasum,&
             ne, nef, net, sion, sionf, siont, uion,&
             sumpl1, sumpl1f, sumpl1t, sumpl2,&
             h2, h2f, h2t,&
             aux, auxf, auxt, info)
        if(info.ne.0) return

        ! abundances don't change within auxiliary *AND*
        ! match_variable loops.
        ifsame_abundances = 1
        ifnear_abundances = 1
        ! update rho to be consistent with calculated rl
        rho = exp(rl)
        tlold = tl
        flold = fl
        if(kif.eq.0) then
           ! Don't need to iterate fl for this case.
        else
           ! There is an if statement below that uses the construct (kif.eq.1.and.pf.lt.0.d0)
           ! where pf is always defined if kif.eq.1, but Fortran makes no guarantees about
           ! the order of evaluation for such constructs so uninitialized pf could
           ! get compared to 0.d0 before the kif test is done.  To avoid uninitialized
           ! use as a matter of principle and to avoid tainting false alarms and
           ! static analysis false alarms, initialize pf here.
           pf = 0._fp_kind
           ! naively there appears to be a similar issue with kif.eq.2
           ! and rf, but in that case rf is defined above by the call
           ! to eos_cold_start which determines rf in eos_calc for
           ! ifnr = 0, stores that result in auxf and therefore rf
           ! here contains that result because of the auxf, rf
           ! equivalance in this routine.  So it is currently an
           ! effective but really ugly (!) communication of rf that
           ! solves this issue, and (fixme) this ugly communication
           ! method should be smoothed out in future.
           if(kif.eq.1) then
              ! there is the possibility of pressure ionization and
              ! Planck-Larkin terms
              if(ifpi_local.eq.0) then
                 ! there are no p, s, u terms from pressure ionization.
                 ppi = 0._fp_kind
                 ppif = 0._fp_kind
              else
                 ! Use ifnr = 0 to be consistent with what eos_cold_start uses internally.
                 ifnr = 0
                 ifnr03 = .true.
                 call pteh_pi_end(&
                      ifnr, full_sum1, rho, rf, rt, t, ne, nef, net,&
                      ppi, ppif, ppit, spi, spif, spit, upi)
              endif
              if(ifexcited.gt.0.and.ifpl.eq.1.and.(ifpi_local.eq.3.or.ifpi_local.eq.4)) then
                 ! For the above combination of flags, the excitation_pi
                 ! code evaluates non-zero values of psum and its
                 ! derivatives so that excitation_pi_end
                 ! returns non-zero values of pexcited and its derivatives.
                 ! n.b. for .not.ifwarm, the simplified free-energy model
                 ! used to calculate the equilibrium constants, dv, is different
                 ! from the full free-energy model used to calculate psum,
                 ! etc., in excitation_sum.
                 call excitation_pi_end(t, rho, rf, rt,&
                      pexcited, pexcitedf, pexcitedt,&
                      sexcited, sexcitedf, sexcitedt, uexcited)
              else
                 ! otherwise, excitation_pi_end would return zero for
                 ! pexcited and its derivatives so save some time by not
                 ! calling it.
                 pexcited = 0._fp_kind
                 pexcitedf = 0._fp_kind
              endif
              ! All the Coulomb stuff is calculated using if_pteh = 1 for this
              ! case.
              ! PTEH (full ionization) approximation to sum0, sum2
              ! reassert this approximation because eos_calc
              ! calls ionize which messes a bit with sum0 and sum2
              sum0ne = full_sum0/full_sum1
              sum0 = sum0ne*n_e
              sum0f = 0._fp_kind
              sum0t = 0._fp_kind
              sum2ne = full_sum2/full_sum1
              sum2 = sum2ne*n_e
              sum2f = 0._fp_kind
              sum2t = 0._fp_kind
              call master_coulomb_end(rhostar,&
                   sum0, sum0ne, sum0f, sum0t,&
                   sum2, sum2ne, sum2f, sum2t,&
                   n_e, t, pstar,&
                   ifcoulomb_mod, if_dc, 1,&
                   dpcoulomb, dpcoulombf, dpcoulombt,&
                   dscoulomb, dscoulombf, dscoulombt, ducoulomb)
              call exchange_end(rhostar, pstar,&
                   pex, pext, pexf,&
                   sex, sexf, sext, uex)
              p0 = cr*rho*t
              ni = full_sum0 - (h2+h2plus)
              pion = ni*p0
              pionf = ni*p0*rf - (h2f+h2plusf)*p0
              pe_cgs = cpe*pe
              pnorad = pe_cgs + pion + ppi + pexcited +&
                   dpcoulomb + pex
              p = pnorad + pr
              if(p.le.0._fp_kind) then
                 if(verbosity.ge.1) then
                    write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
                    write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
                    write(stderr,*) 'eps ='
                    write(stderr,'(1p5e25.15e4)') eps
                    write(stderr,*) 'free_eos_detailed ERROR: (1) negative p calculated for the eos_cold_start case'
                 endif
                 info = info_offset_free_eos_detailed + 1
                 return
              endif
              ! d ln p(f,t)/d ln f
              pf = (pe_cgs*pef + pionf + ppif + pexcitedf + dpcoulombf + pexf)/p
              if(ifrad.gt.1) then
                 paa = log(pnorad) - match_variable
                 ! the p/pnorad factor converts ln ptotal derivative
                 ! to ln pnorad derivative.
                 pnoradf = pf*p/pnorad
                 pac = -paa/(pnoradf)
              else
                 paa = log(p) - match_variable
                 pac = -paa/pf
              endif
           elseif(kif.eq.2) then
              paa = log(rho) - match_variable
              pac = -paa/rf
           else
              error stop 'free_eos_detailed: kif must be 0, 1, or 2'
           endif
           ! Apply limits to potential fl change.
           if(fl + pac .lt. -10._fp_kind) then
              ! Cannot get into much trouble for derived  fl < -10.
              pab = min(20._fp_kind,max(-20._fp_kind,pac))
           else
              lnrho_5 = log(rho) - 1.5_fp_kind*(tl - log(1.e5_fp_kind))
              if(tl.gt.log(1.e6_fp_kind).or.lnrho_5 + rf*min(10._fp_kind,pac).lt.log(1.e-3_fp_kind)) then
                 pab = min(10._fp_kind,max(-10._fp_kind,pac))
              elseif(ifnoaux_iteration.or.lnrho_5 + rf*min(3._fp_kind,pac).lt.log(1.e-2_fp_kind)) then
                 pab = min(3._fp_kind,max(-3._fp_kind,pac))
              elseif(lnrho_5 + rf*min(1._fp_kind,pac).lt.log(1.e-1_fp_kind)) then
                 pab = min(1._fp_kind,max(-1._fp_kind,pac))
              else
                 pab = min(0.5_fp_kind,max(-0.5_fp_kind,pac))
              endif
           endif
           if(paa.ge.0._fp_kind) then
              paaplus = paa
              flplus = fl
           else
              paaminus = paa
              flminus = fl
           endif
           if(paaminus.gt.-1.e30_fp_kind.and.paaplus.lt.1.e30_fp_kind) then
              ! have bracket already!
              ! make sure step would keep within bracket
              if((flminus.le.fl+pab.and.fl+pab.le.flplus).or.&
                   (flminus.ge.fl+pab.and.fl+pab.ge.flplus)) then
                 fl = fl + pab
              else
                 fl = 0.5_fp_kind*(flminus+flplus)
              endif
           else
              ! if no bracket, yet, then normal step is only allowed for
              ! when pf (or rf) is positive.  that is local paa
              ! is following global trend that paa
              ! and pl (or rl) generally increases with fl.
              ! if pf (or rf) negative with no bracket,
              ! then move out of region in the direction
              ! which should produce (global) bracket.
              if((kif.eq.1.and.pf.lt.0._fp_kind).or.(kif.eq.2.and.rf.lt.0._fp_kind)) then
                 if(fl.eq.flminus) then
                    pab = 0.5_fp_kind
                 else
                    pab = -0.5_fp_kind
                 endif
              endif
              fl = fl + pab
           endif
        endif
        if(ieee_support_nan(1._fp_kind)) then
           if(ieee_is_nan(fl)) then
              if(verbosity.ge.1) then
                 write(stderr,*) 'match_variable/ln(10), tl/ln(10) ='
                 write(stderr,'(1p5e25.15e4)') match_variable/ln10, tl/ln10
                 write(stderr,*) 'paa, pac, pab, fl = ', paa, pac, pab, fl
                 write(stderr,*) 'free_eos_detailed ERROR: (2) NaN value of fl detected during cold start iteration'
              endif
              info = info_offset_free_eos_detailed + 2
              return
           endif
        endif
        if(.false..and.verbosity.ge.4) then
           write(stderr,'(a,/,i5,1p5e25.15e4)') 'lter, flold, paa, pac, pab, fl =', lter, flold, paa, pac, pab, fl
        endif
        ! End of iteration on fl to match match_variable with
        ! eos_cold_start
     enddo
     if(kif.gt.0.and.lter.ge.ltermax) then
        if(verbosity.ge.1) then
           write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
           write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
           write(stderr,*) 'eps ='
           write(stderr,'(1p5e25.15e4)') eps
           write(stderr,*) 'flplus, paa(flplus), flminus, paa(flminus) ='
           write(stderr,'(1p5e25.15e4)') flplus, paaplus, flminus, paaminus
           write(stderr,'(a,/,i5,1p5e25.15e4)') 'lter, flold, paa, pac, pab, fl =', lter, flold, paa, pac, pab, fl
           write(stderr,*) 'free_eos_detailed ERROR: (3) fl iteration did not converge for the eos_cold_start case'
        endif
        info = info_offset_free_eos_detailed + 3
        return
     endif
     ! end of if(ifnoaux_iteration.or..not.ifwarm) then branch
  endif
  if(.not.ifnoaux_iteration) then
     lw_old = lw
     ifexcited_old = ifexcited
     nmax_old = nmax
     ifh2_old = ifh2
     ifh2plus_old = ifh2plus
     morder_old = morder
     ifmtrace_old = ifmtrace
     ifcoulomb_old = ifcoulomb
     ifmodified_old = ifmodified
     match_variableold = match_variable
     ! set ifionized.
     if(lw.gt.3) then
        ! partial ionization with all stages of ionization
        ! of all elements.
        ifionized = 0
     elseif(lw.eq.3) then
        ! full ionization of all trace metals.
        ifionized = 1
     else
        ! full ionization of all elements.
        ifionized = 2
     endif
     if(ifwarm) then
        ! warm start:
        ! N.B. This branch of the code depends on many variables/arrays saved
        ! from a previous call to free_eos_detailed.

        ! Fix ancient bug originally found by J. Christensen-Dalsgaard.
        ! Make rho consistent with saved local variable rl (equivalenced to aux(4)).
        ! This is necessary because rho is not a saved local variable (i.e., it
        ! is one of the arguments), and the calling routine likely does not save rho
        ! (which exposed the bug in the jcd case).
        ! rf and rt are okay because they are already saved local variables (equivalenced
        ! to auxf(4) and auxt(4)).
        rho = exp(rl)

        ! To shut up valgrind, fix inconsequential uninitialized
        ! variable issues for lambda and ne.  This fix implies
        ! lambda_old and ne_old are initialized to rediculous values
        ! rather than "random" unitialized values below for the ifwarm
        ! case so that the initial convergence iteration check of
        ! lambda-lambda_old and ne-ne_old done below consistently
        ! fails (as it does in any case) due to these rediculous (but
        ! properly initialized to shut up valgrind!) values.  N.B. use
        ! the "moderate" ridiculous value of 1.d30 rather than 1.d300
        ! to avoid overflows below when taking the ratio of
        ! lambda_old/lambda or ne_old/ne.
        lambda = 1.e30_fp_kind
        ne = 1.e30_fp_kind

        ! For warm start, previous calculation of eos was done with
        ! eos_tqft which produces nu form of many auxiliary variables.
        do index_aux = 1,n_partial_aux
           ! N.B. iaux points only to auxiliary variables that are relevant
           ! for the current set of flags that describe the free-energy model.
           iaux = partial_aux(index_aux)
           if(iaux.ne.4.and.iaux.ne.6.and.iaux.ne.7) then
              ! for all but rl, sum0, and sum2 ... convert
              ! to n form.
              aux(iaux) = aux(iaux)*rho*avogadro
              auxf(iaux) = auxf(iaux)*rho*avogadro + aux(iaux)*rf
              auxt(iaux) = auxt(iaux)*rho*avogadro + aux(iaux)*rt
           endif
           ! Logarithmic Taylor series unless auxiliary variable is zero.
           ! At this point flold is the fl value determined by the last call to
           ! free_eos_detailed.
           if(iaux.eq.rl_aux_index) then
              ! already in log form.
              faux(index_aux) =&
                   (fl-flold)*auxf(iaux) + (tl-tlold)*auxt(iaux)
              if(abs(faux(index_aux)).le.faux_limit_small_start)&
                   aux(iaux) = aux(iaux) + faux(index_aux)
           elseif(aux(iaux).gt.0._fp_kind) then
              faux(index_aux) =&
                   ((fl-flold)*auxf(iaux) + (tl-tlold)*auxt(iaux))/&
                   aux(iaux)
              if(abs(faux(index_aux)).le.faux_limit_small_start)&
                   aux(iaux) = exp(log(aux(iaux)) + faux(index_aux))
           elseif(iaux.gt.maxnextrasum.and.aux(iaux).lt.0._fp_kind) then
              faux(index_aux) =&
                   ((fl-flold)*auxf(iaux) + (tl-tlold)*auxt(iaux))/&
                   aux(iaux)
              if(abs(faux(index_aux)).le.faux_limit_small_start)&
                   aux(iaux) = -exp(log(-aux(iaux)) + faux(index_aux))
           else
              aux(iaux) = 0._fp_kind
           endif
        enddo
     elseif(kif.eq.1) then
        ! for cold start must recalculate eos_cold_start with correct
        ! ifnuform = .false. for this special case.  See ifnuform
        ! shenanigans for kif=1 above.
        ! N.B. here ifwarm is false so ifsame_under will be identical to
        ! the value used for the above call to eos_cold_start.
        ifnuform = .false.
        iteration_count = iteration_count + 1
        call eos_cold_start(&
             verbosity, ifsame_under, lambda, gamma_e,&
             partial_ions, f, eta, wf, t, n_e, pstar,&
             dve_exchange, dve_exchangef, dve_exchanget,&
             full_sum0, full_sum1, full_sum2, charge, dv_pl, dv_plt,&
             ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
             ifsame_zero_abundances, ifexcited, ifnuform,&
             inv_aux,&
             inv_ion, max_index,&
             partial_elements(1:n_partial_elements+2), ion_end,&
             ifionized, if_pteh, if_mc, ifreducedmass,&
             ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,&
             ifpl, ifmodified, ifh2, ifh2plus,&
             izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
             eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
             r_ion3, r_neutral,&
             ifelement, dvzero, dv, dvf, dvt, nion,&
             rhostar,&
             nextrasum,&
             ne, nef, net, sion, sionf, siont, uion,&
             sumpl1, sumpl1f, sumpl1t, sumpl2,&
             h2, h2f, h2t,&
             aux, auxf, auxt, info)
        if(info.ne.0) return
     endif
     tlold = tl
     flold = fl
     do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
        ! initially only allow small magnitude above aux_underflow
        ! for zero to non-zero auxiliary variable change.
        zerolim(iaux) = 1.e10_fp_kind*aux_underflow
     enddo
     ! save values for later use to start bfgs iteration.
     fl_restore = fl
     do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
        aux_restore(iaux) = aux(iaux)
     enddo
     ! this factor must be unity in order to judge
     ! the initial relative change in ne.
     faux_scale = 1._fp_kind
     ! initialize logic for counting number of zero trys.
     ioncountzero = -1
     ifzerocount = 0
     faux_limit_small = faux_limit_small_start
     maxfaux = 1._fp_kind    !assure at least two iterations
     maxfaux_diag = 1._fp_kind !needed to initialize maxfaux_diag_old
     ioncount = 0
     iflast = 0
     ! if not straight through (iflast.ne.1) then either warm
     ! start (with all auxiliary variables defined from  previous call)
     ! or else cold start (with all auxiliary variables defined by
     ! above call to eos_cold_start using pteh_pi and pteh sum approximation).

     ! For first iteration assume NR solution unless simple iteration
     ! solution criteria are fulfilled.
     any_ifsimple = .false.
     simple_lambda(1:njacobian) = 0._fp_kind
     do while(iflast.ne.1.and.ioncount.lt.maxioncount)
        ! find themodynamically consistent set of fermi-dirac integrals and
        ! put the values and their derivatives into rhostar, pstar,
        ! sstar, and ustar which are equivalenced to re, pe, se, and ue
        ! and their f and t derivatives.  Also calculate exchange (when
        ! ifexchange_in > 0) and its effects on dv.
        call master_exchange(verbosity, fl, tl,&
             rhostar, pstar, sstar, ustar, morder,&
             ifexchange_in, dve_exchange, dve_exchangef, dve_exchanget)
        f = exp(fl)
        wf = sqrt(1._fp_kind + f)
        eta = fl+2._fp_kind*(wf-log(1._fp_kind+wf))
        ! number density of free electrons and electron pressure
        n_e = c_e*re
        ! rho/mu_e = n_e H = cd*re
        rmue = cd*re
        ioncount = ioncount + 1
        maxfaux_diag_old = maxfaux_diag
        ! either cold or warm start may require huge ln changes in
        ! auxiliary variables.  However, these huge ln changes
        ! are often only required of auxiliary variables that
        ! are approaching zero and which have little practical
        ! effect on the equilibrium constants.  If we use a
        ! over-cautious approach when limiting the maximum
        ! ln change, then a large number of auxiliary variable
        ! iterations will be required.  However, with our present
        ! approach of allowing large ln changes, the maximum
        ! number of iterations is usually less than 10.
        ! n.b. use one extra iteration after satisfy ending criteria.
        if(.not.any_ifsimple.and.&
             ((ifwarm.and.ioncount.ge.maxioncount).or.&
             ioncount.ge.maxioncount.or.&
             abs(maxfaux).le.eps_aux))&
             iflast = 1
        ! if near final stages of auxiliary variable convergence,
        ! make sure to fix underflow criterion in ionize called by eos_calc.
        ! if(ioncount.ge.2.and.&
        !   max(abs(rl-rl_old),&
        !   abs(ne-ne_old)/ne,&
        !   abs(lambda-lambda_old)/max(1.d-15,lambda))/&
        !   faux_scale.le.1.d-2) then
        ! N.B. maxfaux = 1. on first time through this loop so ifsame_under
        !   always false on first time.
        if(ioncount.ge.2.and.abs(maxfaux).le.1.e-3_fp_kind) then
           ifsame_under = .true.
        else
           ifsame_under = .false.
        endif
        if(.false..and.verbosity.ge.4) then
           write(stderr,*) 'aux'
           write(stderr,*) (aux(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
        endif
        if(debug_dv_aux.or.debug_jacobian) then
           ! N.B. as the result of the prior debug logic, for every call of free_eos_detailed
           ! that reaches here, match_variable and tl are the same, delta1 and delta2 take
           ! on the expected pattern for the 33 calls (both zero for the first call, and only
           ! one non-zero for the remaining 32 calls), and fl and the set of
           ! auxiliary variables is the same (via the forced debug call of eos_cold_start for every
           ! call of free_eos_detailed).
           ! specify change (if any) in that standard set of auxiliary variables (or fl or tl)
           ! and leave dv zero point as zero.
           ifdvzero = .true.
           ! Choose itemp1 and itemp2 to correspond to auxiliary variables
           ! actually used for particular free-energy model being tested.
           ! For example, EOS1 has 15 auxiliary variables in order
           ! sum0, sum2, 7 neutral auxiliary variables, 2 ion auxiliary
           ! variables, and 4 xextrasum auxiliary variables.

           ! FIXME (2022).  The code below works perfectly for the
           ! debug_dv_aux case for all values of itemp1 and itemp2 and
           ! also works perfectly for the debug_jacobian case for all
           ! values of itemp1 *except for* where fl and tl differences
           ! are used to check the partials of faux (the RHS vector)
           ! wrt fl and tl.  Note, the ifnr = 0 eos_tqft case does
           ! calculate the partial derivatives of faux wrt to fl and
           ! tl correctly (otherwise the normal non-debug
           ! free_eos_test results would show derivative errors) so my
           ! best guess is the source of this debug_jacobian issue is
           ! my attempt to mimic those eos_tqft results with ifnr=3
           ! eos_jacobian calculations has not yet been successful.
           ! Of course, this is a bug in debug code that is normally
           ! not used so it is a low priority to fix this bug.
           if(.false.) then
              itemp1 = 1
              itemp2 = 2
           elseif(.true.) then
              itemp1 = 3
              itemp2 = 4
           elseif(.false.) then
              itemp1 = 5
              itemp2 = 6
           elseif(.false.) then
              itemp1 = 7
              itemp2 = 8
           elseif(.false.) then
              ! repeat 8 for alignment of like results.  8 and 9 are
              ! last two derivatives wrt neutral sums.
              itemp1 = 8
              itemp2 = 9
           elseif(.false.) then
              itemp1 = 10
              itemp2 = 11
           elseif(.false.) then
              itemp1 = 12
              itemp2 = 13
           elseif(.false.) then
              itemp1 = 14
              itemp2 = 15
           elseif(.true.) then
              itemp1 = n_partial_aux + 1
              itemp2 = n_partial_aux + 2
           endif
           if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
              iaux1 = partial_aux(itemp1)
              if(iaux1.eq.4) then
                 aux(iaux1) = aux(iaux1) + delta1 !density
                 rho = exp(rl)
              else
                 aux(iaux1) = aux(iaux1)*(1._fp_kind + delta1)
              endif
           elseif(itemp1.eq.n_partial_aux + 1) then
              fl = fl + delta1
           elseif(itemp1.eq.n_partial_aux + 2) then
              tl = tl + delta1
           else
              error stop 'free_eos_detailed: bad itemp1'
           endif
           if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
              iaux2 = partial_aux(itemp2)
              if(iaux2.eq.4) then
                 aux(iaux2) = aux(iaux2) + delta2 !density
                 rho = exp(rl)
              else
                 aux(iaux2) = aux(iaux2)*(1._fp_kind + delta2)
              endif
           elseif(itemp2.eq.n_partial_aux + 1) then
              fl = fl + delta2
           elseif(itemp2.eq.n_partial_aux + 2) then
              tl = tl + delta2
           else
              error stop 'free_eos_detailed: bad itemp2'
           endif
           ! just in case there was a change in tl
           t = exp(tl)
           tc2 = c2/t
           ! calculate planck-larkin occupation probabilities and equilibrium
           ! constant changes.
           if(ifpl.eq.1)&
                call pl_prepare(partial_ions(:n_partial_ions+2), max_index,&
                nion, ifdv,&
                bi, tc2, plop, plopt, plopt2, dv_pl, dv_plt)
           if(ifrad.ge.1) then
              pr = prad_const*t**4
           else
              pr = 0._fp_kind
           endif
           ! just in case there was a change in fl or tl.
           call master_exchange(verbosity, fl, tl,&
                rhostar, pstar, sstar, ustar, morder,&
                ifexchange_in, dve_exchange, dve_exchangef,&
                dve_exchanget)
           f = exp(fl)
           wf = sqrt(1._fp_kind + f)
           eta = fl+2._fp_kind*(wf-log(1._fp_kind+wf))
           ! number density of free electrons and electron pressure
           n_e = c_e*re
           ! rho/mu_e = n_e H = cd*re
           rmue = cd*re
        elseif(debug_aux_dv) then
           ! choose first two dv or x values to vary
           ! partial derivative of auxiliary variables wrt dv
           ! or partial derivative of y wrt x.
           ! 1 ==> 4 indices correspond to H+, He+, He++, C+, while
           ! max_index-1 ==> max_index correspond to H2 and H2+.
           if(.false.) then
              ! H+, He+
              itemp1 = 1
              itemp2 = 2
           elseif(.false.) then
              ! He++, C+
              itemp1 = 3
              itemp2 = 4
           elseif(.true.) then
              ! H2, H2+
              itemp1 = max_index-1
              itemp2 = max_index
           endif
           ! the next part of this initial step for this debugging mode is
           ! done inside eos_jacobian after consistent dv values are calculated as
           ! a function of aux_old.
        endif
        ! Allocate debug_results with correct size (3 by defined compact size of function vector
        ! whose partial derivative is being tested) if any of the debug_* parameters is set to .true.
        if(debug_dv_aux) then
           ! Function vector is dv.
           allocate(debug_results(3, max_index))
        elseif(debug_aux_dv) then
           ! Function vector is aux.
           allocate(debug_results(3, n_partial_aux))
        elseif(debug_jacobian) then
           ! Function vector is faux.
           allocate(debug_results(3, njacobian))
        endif

        ! rl is a locally saved variable, but ne and lambda are not (they are arguments)
        ! so set ne_old and lambda_old to either the eos_cold_start results or
        ! inconsequential ridiculous values for the first iteration for the ifwarm case.
        rl_old = rl
        ne_old = ne
        lambda_old = lambda
        do index_aux = 1,n_partial_aux
           iaux = partial_aux(index_aux)
           aux_old(iaux) = aux(iaux)
           ! Adjust old_allow_log variables for changed aux_old.

           ! n.b. auxiliary variable corresponding to jaux = rl_aux_index is already in
           ! log form (i.e., rl = log(rho)).
           ! If element of old_allow_log is .true. can take log of
           ! corresponding element of aux_old.
           old_allow_log(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).gt.0._fp_kind
           ! If element of old_allow_log_neg is .true. can take log of
           ! corresponding negative element of aux_old.
           old_allow_log_neg(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).lt.0._fp_kind
        enddo
        flold = fl
        if(kif.eq.0.and.iflast.ne.1) then
           ! No need to calculate fl and tl derivatives until last
           ! iteration for kif = 0
           ifnr = 1
        else
           ! Calculate combined auxiliary variable and fl (and ft)
           ! partial derivatives.
           ifnr = 3
        endif
        ifnr03 = ifnr.eq.0.or.ifnr.eq.3
        if(.true..and.kif.eq.2.and.ioncount.eq.100) then
           if(verbosity.ge.3) write(stderr,*) 'start one-time emergency eos_bfgs per NR cycle'
           if(.true.) then
              fl = fl_restore
              do index_aux = 1,n_partial_aux
                 iaux = partial_aux(index_aux)
                 aux_old(iaux) = aux_restore(iaux)
                 ! Adjust old_allow_log variables for changed aux_old.

                 ! n.b. auxiliary variable corresponding to jaux = rl_aux_index is already in
                 ! log form (i.e., rl = log(rho)).
                 ! If element of old_allow_log is .true. can take log of
                 ! corresponding element of aux_old.
                 old_allow_log(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).gt.0._fp_kind
                 ! If element of old_allow_log_neg is .true. can take log of
                 ! corresponding negative element of aux_old.
                 old_allow_log_neg(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).lt.0._fp_kind
              enddo
              isimple_bfgs = 1
              ! ordinarily do not do any preliminary simple iterations at fixed
              ! fl prior to using the bfgs technique.
              do while(isimple_bfgs.le.0)
                 !temporary.
                 !do while(isimple_bfgs.le.5)
                 ! No debug mode ever reaches here so there are no debug arguments!
                 call eos_jacobian(&
                      verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
                      ifrad, match_variable, kif, fl,&
                      aux_old, aux, auxf, auxt, aux_dv,&
                      njacobian,&
                      rhs1, jacobian, p, pr,&
                      dv_aux,&
                      lambda, gamma_e,&
                      ifnr, nux, nuy, nuz,&
                      n_partial_ions,&
                      sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
                      partial_ions, f, eta, wf, t, n_e, pstar,&
                      dve_exchange, dve_exchangef, dve_exchanget,&
                      full_sum0, full_sum1, full_sum2, charge, charge2,&
                      dv_pl, dv_plt,&
                      ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
                      ifsame_zero_abundances,&
                      ifexcited, ifsame_under, if_free_non_ideal_calc,&
                      inv_aux,&
                      inv_ion, max_index,&
                      partial_elements(1:n_partial_elements+2), ion_end, mion_end,&
                      ifionized, if_pteh, if_mc, ifreducedmass,&
                      ifsame_abundances, ifmtrace,&
                      iatomic_number, ifpi_local,&
                      ifpl, ifmodified, ifh2, ifh2plus,&
                      izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
                      eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
                      r_ion3, r_neutral,&
                      ifelement, dvzero, dv, dvf, dvt, nion,&
                      rhostar,&
                      ne, nef, net, sion, sionf, siont, uion,&
                      pnorad, pnoradf, pnoradt, pnorad_aux(1:n_partial_aux+2),&
                      free_fp, free_aux(1:n_partial_aux+2),&
                      nextrasum,&
                      sumpl1, sumpl1f, sumpl1t, sumpl2,&
                      h2, h2f, h2t, h2_dv, info)
                 if(info.ne.0) return

                 if(.false..and.verbosity.ge.4) then
                    write(stderr,*) 'pre-bfgs isimple_bfgs = ', isimple_bfgs
                    write(stderr,*) 'old pre-bfgs value of aux for frozen fl'
                    write(stderr,'(5(1pe15.5e4,2x))') (aux_old(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
                    write(stderr,*) 'new pre-bfgs value of aux for frozen fl'
                    write(stderr,'(5(1pe15.5e4,2x))') (aux(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
                 endif
                 do index_aux = 1,n_partial_aux
                    iaux = partial_aux(index_aux)
                    aux_old(iaux) = aux(iaux)
                    ! Adjust old_allow_log variables for changed aux_old.

                    ! n.b. auxiliary variable corresponding to jaux = rl_aux_index is already in
                    ! log form (i.e., rl = log(rho)).
                    ! If element of old_allow_log is .true. can take log of
                    ! corresponding element of aux_old.
                    old_allow_log(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).gt.0._fp_kind
                    ! If element of old_allow_log_neg is .true. can take log of
                    ! corresponding negative element of aux_old.
                    old_allow_log_neg(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).lt.0._fp_kind
                 enddo
                 isimple_bfgs = isimple_bfgs + 1
              enddo
           endif
           call eos_bfgs(&
                verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
                ifrad, match_variable, kif, fl,&
                aux_old, aux, auxf, auxt, aux_dv,&
                njacobian,&
                rhs1, jacobian, p, pr,&
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
                partial_elements(1:n_partial_elements+2), n_partial_elements, ion_end, mion_end,&
                ifionized, if_pteh, if_mc, ifreducedmass,&
                ifsame_abundances, ifmtrace,&
                iatomic_number, ifpi_local,&
                ifpl, ifmodified, ifh2, ifh2plus,&
                izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
                eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
                r_ion3, r_neutral,&
                ifelement, dvzero, dv, dvf, dvt, nion,&
                ne, nef, net, sion, sionf, siont, uion,&
                pnorad, pnoradf, pnoradt, pnorad_aux(1:n_partial_aux+2),&
                free_fp, free_aux(1:n_partial_aux+2),&
                nextrasum,&
                sumpl1, sumpl1f, sumpl1t, sumpl2,&
                h2, h2f, h2t, h2_dv, info)
           if(info.ne.0) return
           call master_exchange(verbosity, fl, tl,&
                rhostar, pstar, sstar, ustar, morder,&
                ifexchange_in, dve_exchange, dve_exchangef,&
                dve_exchanget)
           flold = fl
           f = exp(fl)
           wf = sqrt(1._fp_kind + f)
           eta = fl+2._fp_kind*(wf-log(1._fp_kind+wf))
           ! number density of free electrons and electron pressure
           n_e = c_e*re
           ! rho/mu_e = n_e H = cd*re
           rmue = cd*re
           isimple_bfgs = 1
           ! ordinarily do not do any post-BFGS simple iterations at fixed
           ! fl prior to using the NR technique.
           do while(isimple_bfgs.le.0)
              !temporary.
              !do while(isimple_bfgs.le.5)
              ! No debug mode ever reaches here so there are no debug arguments!
              call eos_jacobian(&
                   verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
                   ifrad, match_variable, kif, fl,&
                   aux_old, aux, auxf, auxt, aux_dv,&
                   njacobian,&
                   rhs1, jacobian, p, pr,&
                   dv_aux,&
                   lambda, gamma_e,&
                   ifnr, nux, nuy, nuz,&
                   n_partial_ions,&
                   sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
                   partial_ions, f, eta, wf, t, n_e, pstar,&
                   dve_exchange, dve_exchangef, dve_exchanget,&
                   full_sum0, full_sum1, full_sum2, charge, charge2,&
                   dv_pl, dv_plt,&
                   ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
                   ifsame_zero_abundances,&
                   ifexcited, ifsame_under, if_free_non_ideal_calc,&
                   inv_aux,&
                   inv_ion, max_index,&
                   partial_elements(1:n_partial_elements+2), ion_end, mion_end,&
                   ifionized, if_pteh, if_mc, ifreducedmass,&
                   ifsame_abundances, ifmtrace,&
                   iatomic_number, ifpi_local,&
                   ifpl, ifmodified, ifh2, ifh2plus,&
                   izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
                   eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
                   r_ion3, r_neutral,&
                   ifelement, dvzero, dv, dvf, dvt, nion,&
                   rhostar,&
                   ne, nef, net, sion, sionf, siont, uion,&
                   pnorad, pnoradf, pnoradt, pnorad_aux(1:n_partial_aux+2),&
                   free_fp, free_aux(1:n_partial_aux+2),&
                   nextrasum,&
                   sumpl1, sumpl1f, sumpl1t, sumpl2,&
                   h2, h2f, h2t, h2_dv, info)
              if(info.ne.0) return

              if(.false..and.verbosity.ge.4) then
                 write(stderr,*) 'post-bfgs isimple_bfgs = ', isimple_bfgs
                 write(stderr,*) 'old post-bfgs value of aux for frozen fl'
                 write(stderr,'(5(1pe15.5e4,2x))') (aux_old(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
                 write(stderr,*) 'new post-bfgs value of aux for frozen fl'
                 write(stderr,'(5(1pe15.5e4,2x))') (aux(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
              endif
              do index_aux = 1,n_partial_aux
                 iaux = partial_aux(index_aux)
                 aux_old(iaux) = aux(iaux)
                 ! Adjust old_allow_log variables for changed aux_old.

                 ! n.b. auxiliary variable corresponding to jaux = rl_aux_index is already in
                 ! log form (i.e., rl = log(rho)).
                 ! If element of old_allow_log is .true. can take log of
                 ! corresponding element of aux_old.
                 old_allow_log(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).gt.0._fp_kind
                 ! If element of old_allow_log_neg is .true. can take log of
                 ! corresponding negative element of aux_old.
                 old_allow_log_neg(index_aux) = iaux.ne.rl_aux_index.and.aux_old(iaux).lt.0._fp_kind
              enddo
              isimple_bfgs = isimple_bfgs + 1
           enddo !do while(isimple_bfgs.le.0)
           ! restore initial size of faux_limit_small.
           faux_limit_small = faux_limit_small_start
           ! For first post-bfgs iteration assume NR solution unless
           ! simple iteration solution criteria are fulfilled.
           any_ifsimple = .false.
           simple_lambda(1:njacobian) = 0._fp_kind
           ! force at least two NR iterations after bfgs minimization.
           iflast = 0
           if(verbosity.ge.3) write(stderr,*) 'end one-time emergency eos_bfgs per NR cycle'
        endif !if(.true..and.kif.eq.2.and.ioncount.eq.100) then
        iteration_count = iteration_count + 1
        ! This call has debug arguments just in case debug_any is .true.
        call eos_jacobian(&
             verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
             ifrad, match_variable, kif, fl,&
             aux_old, aux, auxf, auxt, aux_dv,&
             njacobian,&
             rhs1, jacobian, p, pr,&
             dv_aux,&
             lambda, gamma_e,&
             ifnr, nux, nuy, nuz,&
             n_partial_ions,&
             sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
             partial_ions, f, eta, wf, t, n_e, pstar,&
             dve_exchange, dve_exchangef, dve_exchanget,&
             full_sum0, full_sum1, full_sum2, charge, charge2,&
             dv_pl, dv_plt,&
             ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
             ifsame_zero_abundances,&
             ifexcited, ifsame_under, if_free_non_ideal_calc,&
             inv_aux,&
             inv_ion, max_index,&
             partial_elements(1:n_partial_elements+2), ion_end, mion_end,&
             ifionized, if_pteh, if_mc, ifreducedmass,&
             ifsame_abundances, ifmtrace,&
             iatomic_number, ifpi_local,&
             ifpl, ifmodified, ifh2, ifh2plus,&
             izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
             eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
             r_ion3, r_neutral,&
             ifelement, dvzero, dv, dvf, dvt, nion,&
             rhostar,&
             ne, nef, net, sion, sionf, siont, uion,&
             pnorad, pnoradf, pnoradt, pnorad_aux(1:n_partial_aux+2),&
             free_fp, free_aux(1:n_partial_aux+2),&
             nextrasum,&
             sumpl1, sumpl1f, sumpl1t, sumpl2,&
             h2, h2f, h2t, h2_dv, info,&
             debug_aux_dv, debug_dv_aux, debug_jacobian, itemp1, itemp2, delta1, delta2, debug_results)
        if(info.ne.0) return

        if(debug_any) then
           ! Fill in uninitialized intent(out) variables so that they
           ! don't obfuscate detection of any remaining uninitialized
           ! variables for the debug_any case.
           fm = 0._fp_kind
           ft = 0._fp_kind
           rlout = rl
           p = 0._fp_kind
           pl = 0._fp_kind
           cf = 0._fp_kind
           cp = 0._fp_kind
           sf = 0._fp_kind
           st = 0._fp_kind
           grada = 0._fp_kind
           rtp = 0._fp_kind
           rmue = 0._fp_kind
           fh2 = 0._fp_kind
           fhe2 = 0._fp_kind
           fhe3 = 0._fp_kind
           xmu1 = 0._fp_kind
           xmu3 = 0._fp_kind
           gamma1 = 0._fp_kind
           gamma2 = 0._fp_kind
           gamma3 = 0._fp_kind
           h2rat = 0._fp_kind
           h2plusrat = 0._fp_kind
           lambda = 0._fp_kind
           gamma_e = 0._fp_kind
           sound2 = 0._fp_kind

           if(debug_dv_aux) then
              ! partial derivative of dv wrt auxiliary variables, fl, or tl.
              ! N.B. dv, dvf, and dvt have non-compact indices, and the compact indices below
              ! correspond to H+, He+, He++, C+, and H2 and H2+ when H, He, and C have non-zero
              ! abundances.
              jtemp1 = 1
              jtemp2 = 2
              jtemp3 = 3
              jtemp4 = 4
              jtemp5 = max_index-1
              jtemp6 = max_index
           elseif(debug_aux_dv) then
              ! test of aux_dv derivatives
              ! Choose jtemp1 ==> jtemp6 to correspond to auxiliary variables
              ! actually used for particular free-energy model being tested.
              ! reminder, there are only 15 of these variables for EOS1 in
              ! order sum0, sum2, 7 neutral sums, 2 ionized sums, and 4 xextra
              ! sums.
              if(.false.) then
                 jtemp1 = 1
                 jtemp2 = 2
                 jtemp3 = 3
                 jtemp4 = 4
                 jtemp5 = 5
                 jtemp6 = 6
              elseif(.false.) then
                 jtemp1 = 7
                 jtemp2 = 8
                 jtemp3 = 9
                 jtemp4 = 10
                 jtemp5 = 11
                 jtemp6 = 12
              elseif(.true.) then
                 jtemp1 = 10
                 jtemp2 = 11
                 jtemp3 = 12
                 ! prior to this are repeats.
                 jtemp4 = 13
                 jtemp5 = 14
                 jtemp6 = 15
              endif
           elseif(debug_jacobian) then
              ! Choose jtemp[1-6] to correspond to auxiliary variables actually
              ! used for the particular free-energy model being tested.
              ! Reminders: There are only n_partial_aux = 15 of these variables
              ! for EOS1 in order sum0, sum2, 7 neutral sums, 2 ionized sums,
              ! and 4 xextra sums, and njacobian = n_partial_aux + 1 = 16 to
              ! account for extra faux equation to be zeroed assuming kif > 0.
              ! N.B. The first ".true." logical block below sets jtemp[1-6].
              ! Distinct values must be assigned to those variables within that
              ! block (and therefore all blocks that are ever going to become
              ! the ".true." one) because otherwise the auxf and auxt
              ! transformations below will (incorrectly) be applied more than
              ! once for that block.
              if(.true.) then
                 jtemp1 = 1
                 jtemp2 = 2
                 jtemp3 = 3
                 jtemp4 = 4
                 jtemp5 = 5
                 jtemp6 = 6
              elseif(.false.) then
                 jtemp1 = 7
                 jtemp2 = 8
                 jtemp3 = 9
                 jtemp4 = 10
                 jtemp5 = 11
                 jtemp6 = 12
              elseif(.false.) then
                 ! Assume n_partial_aux + 1 = 16
                 ! 11 and 12 repeat to get good alignment

                 jtemp1 = 11
                 jtemp2 = 12
                 jtemp3 = 13
                 jtemp4 = njacobian-2
                 jtemp5 = njacobian-1
                 jtemp6 = njacobian
                 ! FIXME(2022).  faux(njacobian) has different form
                 ! then the rest.  The jacobian matrix is calculated in
                 ! a way that is consistent with that so debugging for
                 ! itemp1 or itemp2 corresponding to auxold should be fine,
                 ! but the special case of partial faux(njacobian) wrt fl or tl
                 ! has not yet been implemented in eos_jacobian.
                 if(itemp1.gt.n_partial_aux.or.itemp2.gt.n_partial_aux)&
                      error stop 'eos_jacobian, debug_jacobian: this combination not yet implemented'
              endif
           endif
           degeneracy(1:3) = debug_results(1:3, jtemp1)
           pressure(1:3) = debug_results(1:3, jtemp2)
           density(1:3) = debug_results(1:3, jtemp3)
           energy(1:3) = debug_results(1:3, jtemp4)
           enthalpy(1:3) = debug_results(1:3, jtemp5)
           entropy(1:3) = debug_results(1:3, jtemp6)
           return
        endif

        ! if aux or aux_old is zero (aside from the log density), this
        ! is the result of underflow.  this problem should
        ! be isolated from the rest of the solution as much as possible
        ! by zeroing all rows and columns of the Jacobian that are
        ! concerned with the auxiliary variable that is zero while
        ! the diagonal element is set to unity.
        ! (the assumption here is that one must be careful with auxiliary
        ! variables that are the result of an error.  We do the most
        ! conservative thing which is a simple iteration.  This should leave
        ! the remainder of the solution largely unaffected since *usually*
        ! a transition to or from the underflow condition means the associated
        ! non-zero variable is small.)
        !   n.b. this transformation has already been effectively done
        !   so comment out.
        ! do jndex_aux = 1, n_partial_aux
        !   jaux = partial_aux(jndex_aux)
        !   do index_aux = 1, n_partial_aux
        !     iaux = partial_aux(index_aux)
        !!    zero rows *or* columns....
        !     if(&
        !       ((jaux.ne.rl_aux_index.and.&
        !       (aux(jaux).eq.0.d0.or.aux_old(jaux).eq.0.d0)).or.&
        !       (iaux.ne.rl_aux_index.and.&
        !       (aux(iaux).eq.0.d0.or.aux_old(iaux).eq.0.d0)))&
        !       ) then
        !         if(iaux.eq.jaux) then
        !           jacobian(jndex_aux,index_aux) = 1.d0
        !       else
        !         jacobian(jndex_aux,index_aux) = 0.d0
        !       endif
        !     endif
        !   enddo
        !   if(njacobian.gt.n_partial_aux) then
        !     if(jaux.ne.rl_aux_index.and.(aux(jaux).eq.0.d0.or.aux_old(jaux).eq.0.d0)) then
        !       jacobian(jndex_aux,njacobian) = 0.d0
        !       jacobian(njacobian,jndex_aux) = 0.d0
        !     endif
        !   endif
        ! enddo
        ! Jacobian(i,j) is negative partial ith equation wrt jth auxiliary
        ! variable in (usually) log form.  Decide on whether an auxiliary
        ! variable is of major importance based on whether any off-diagonal
        ! element of its Jacobian column is greater than major_crit * row
        ! norm.
        ifmajor(1:njacobian) = .false.
        rhs1_save(1:njacobian) = rhs1(1:njacobian)
        jacobian_save(1:njacobian,1:njacobian) = jacobian(1:njacobian, 1:njacobian)
        do ijacobian = 1,njacobian
           row_norm = 0._fp_kind
           do jjacobian = 1,njacobian
              row_norm = max(row_norm, abs(jacobian(ijacobian,jjacobian)))
           enddo
           do jjacobian = 1,njacobian
              if(jjacobian.ne.ijacobian)&
                   ifmajor(jjacobian) = ifmajor(jjacobian).or.&
                   (abs(jacobian(ijacobian,jjacobian)).gt.major_crit*row_norm)
           enddo
        enddo
        if(.false..and.verbosity.ge.4) then
           write(stderr,*) 'njacobian, naux = ', njacobian, naux
           do jndex_aux = 1,njacobian
              if(jndex_aux.le.n_partial_aux) then
                 jaux = partial_aux(jndex_aux)
                 write(stderr,*) 'jndex_aux, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg = '
                 write(stderr,'(i5,4l5)')&
                      jndex_aux, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg
                 write(stderr,*) 'jndex_aux, aux_old, aux, rhs1 ='
                 write(stderr,'(i5,1p3e20.10e4)')&
                      jndex_aux, aux_old(jaux), aux(jaux), rhs1(jndex_aux)
              else
                 write(stderr,*) 'jndex_aux, rhs1 ='
                 write(stderr,'(i5,30x,1pe20.10e4)') jndex_aux, rhs1(jndex_aux)
              endif
              write(stderr,*) 'jacobian(jndex_aux, index_aux) ='
              write(stderr,'(5(1pe15.5e4,2x))') (jacobian(jndex_aux, index_aux), index_aux = 1, njacobian)
           enddo
        endif
        rhs1_lapack(1:njacobian) = real(rhs1(1:njacobian),lapack_fp_kind)
        jacobian_lapack(1:njacobian,1:njacobian) = real(jacobian(1:njacobian,1:njacobian), lapack_fp_kind)
        if(if_svd) then
           ! Use singular value decomposition to solve system of equations.
           ! n.b. scaling is done inside if necessary and jacobian and rhs1_lapack *may*
           ! be modified by these scaling factors.
           ! the returned sol_lapack is the solution of the original unscaled system
           ! of equations.
           call solve_linear_svd(verbosity, njacobian, 1, jacobian_lapack, naux, rhs1_lapack, naux,&
                tolerance_lapack, rcond_lapack, sol_lapack, naux, info_lapack)
        else
           ! Use LU factorization solution from lapack to solve system of equations.
           ! n.b. scaling is done inside if necessary and jacobian and rhs1_lapack *may*
           ! be modified by these scaling factors.
           ! the returned sol_lapack is the solution of the original unscaled system
           ! of equations.
           call dgesvx('E', 'N',&
                njacobian, 1, jacobian_lapack, naux,&
                lu_lapack, naux, ipiv_lapack, equed_lapack,&
                row_lapack, col_lapack, rhs1_lapack, naux,&
                sol_lapack, naux,&
                rcond_lapack, ferr_lapack, berr_lapack,&
                work_lapack, iwork_lapack, info_lapack)
        endif
        ! According to the dgesvx man page,
        ! "If the reciprocal of the condition number is less than machine precision,
        ! INFO = N+1 is returned as a warning, but the routine still goes on
        ! to solve for X and compute error bounds...."
        ! Therefore, the "(.not.if_svd.and.info_lapack.gt.njacobian)" clause
        ! in the if statement below is to ignore this warning if it does occur.

        if(.not.(info_lapack.eq.0.or.(.not.if_svd.and.info_lapack.gt.njacobian))) then
           if(verbosity.ge.1) then
              write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
              write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
              write(stderr,*) 'eps ='
              write(stderr,'(1p5e25.15e4)') eps
              write(stderr,*) 'ioncount, ne, ne_old ='
              write(stderr,'(i5,1p5e25.15e4)') ioncount, ne, ne_old
              write(stderr,'(a,/,i5,1p5e25.15e4)') 'info_lapack, rcond_lapack = ', info_lapack, rcond_lapack
              if(if_svd) then
                 write(stderr ,*) 'free_eos_detailed ERROR: (first info_lapack) no solve_linear_svd solution'
              else
                 write(stderr ,*) 'free_eos_detailed ERROR: (first info_lapack) no dgesvx solution'
              endif
           endif
           info = info_offset_lapack + info_lapack
           return
        endif
        if(verbosity.ge.2.and.rcond_lapack.lt.rcond_lapack_min) then
           write(stderr,'(a,/,i5,1p5e25.15e4)')&
                'free_eos_detailed WARNING: (1) ioncount, rcond_lapack = ', ioncount, rcond_lapack
        elseif(.false..and.verbosity.ge.4) then
           write(stderr,'(a,/,i5,1p5e25.15e4)')&
                'free_eos_detailed: (1) ioncount, rcond_lapack = ', ioncount, rcond_lapack
        endif
        maxfaux = 0._fp_kind
        do index_aux = 1, n_partial_aux
           ! find solution
           faux_nr(index_aux) = real(sol_lapack(index_aux,1),fp_kind)
           iaux = partial_aux(index_aux)
           ! 4th auxiliary variable is already in log form so the if
           ! statement selects all auxiliary variables that in both
           ! their old and new form have been treated as logarithmic.
           if((iaux.eq.rl_aux_index.or.old_new_allow_log(index_aux).or.&
                old_new_allow_log_neg(index_aux)).and.&
                abs(faux_nr(index_aux)).gt.abs(maxfaux)) then
              index_max = index_aux
              maxfaux = faux_nr(index_aux)
           endif
        enddo
        ! Perform unnecessary initialization of pac_nr to suppress
        ! spurious gfortran [-Wmaybe-uninitialized] warning message.
        pac_nr = 0._fp_kind
        if(njacobian.gt.n_partial_aux) then
           ! raw change in fl
           pac_nr = real(sol_lapack(njacobian,1),fp_kind)
           if(abs(pac_nr).gt.abs(maxfaux)) then
              index_max = njacobian
              maxfaux = pac_nr
           endif
        endif
        if(.false..and.verbosity.ge.4) then
           if(njacobian.gt.n_partial_aux) then
              write(stderr,*) 'raw NR solution including fl'
              write(stderr,'(5(1pe15.5e4,l2))')&
                   (faux_nr(index_aux), ifmajor(index_aux), index_aux = 1, n_partial_aux), pac_nr, ifmajor(njacobian)
           else
              write(stderr,*) 'raw NR solution with frozen fl'
              write(stderr,'(5(1pe15.5e4,l2))') (faux_nr(index_aux), ifmajor(index_aux), index_aux = 1, n_partial_aux)
           endif
        endif
        any_ifsimple = .false.
        !temporary put no limit on change by disabling following do loop
        !do index_aux = 1, 0
        do index_aux = 1, n_partial_aux
           ! default is the NR solution.
           faux(index_aux) = faux_nr(index_aux)
           ! to control changes far from the solution it is usually
           ! best to use a lower-order solution, e.g., the simple
           ! iteration solution where auxnew = result delivered by
           ! eos_calc.  We only substitute the lower-order solution if
           ! the NR solution is larger than both faux_limit_diag *and* the
           ! lower-order solution.
           iaux = partial_aux(index_aux)
           if(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)) then
              if(abs(faux_nr(index_aux)).gt.max(faux_limit_diag,abs(log(aux(iaux)) - log(aux_old(iaux))))) then
                 any_ifsimple = .true.
                 simple_lambda(index_aux) = simple_lambda_max
              else
                 simple_lambda(index_aux) = simple_lambda_ratio*simple_lambda(index_aux)
                 if(simple_lambda(index_aux).ge.simple_lambda_min.and.abs(maxfaux).ge.simple_lambda_min) then
                    any_ifsimple = .true.
                 else
                    simple_lambda(index_aux) = 0._fp_kind
                 endif
              endif
           else
              ! no log transformation allowed.
              ! iaux = 4 (already logarithmic) or
              ! aux or aux_old was zero or had opposite signs.
              if(abs(faux_nr(index_aux)).gt.max(faux_limit_diag,&
                   abs(aux(iaux)-aux_old(iaux)))) then
                 any_ifsimple = .true.
                 simple_lambda(index_aux) = simple_lambda_max
              else
                 simple_lambda(index_aux) = simple_lambda_ratio*&
                      simple_lambda(index_aux)
                 if(simple_lambda(index_aux).ge.simple_lambda_min.and.&
                      abs(maxfaux).ge.simple_lambda_min) then
                    any_ifsimple = .true.
                 else
                    simple_lambda(index_aux) = 0._fp_kind
                 endif
              endif
           endif
        enddo
        if(njacobian.gt.n_partial_aux) then
           ! default value.
           pac = pac_nr
           simple_lambda(njacobian) = 0._fp_kind
        endif
        if(any_ifsimple) then
           ! multiply RHS and diagonal by (1+lambda).  For large lambda this
           ! has similar effect to zeroing off-diagonal elements.
           rhs1(1:njacobian) = (1._fp_kind + simple_lambda(1:njacobian))*rhs1_save(1:njacobian)
           jacobian(1:njacobian,1:njacobian) = jacobian_save(1:njacobian,1:njacobian)
           do ijacobian = 1,njacobian
              jacobian(ijacobian,ijacobian) = jacobian(ijacobian,ijacobian) + simple_lambda(ijacobian)
           enddo

           rhs1_lapack(1:njacobian) = real(rhs1(1:njacobian),lapack_fp_kind)
           jacobian_lapack(1:njacobian,1:njacobian) = real(jacobian(1:njacobian,1:njacobian), lapack_fp_kind)
           if(if_svd) then
              ! Use SVD solution to solve system of equations.
              ! n.b. scaling is done inside if necessary and jacobian_lapack and rhs1_lapack
              ! *may* be modified by these scaling factors.
              ! the returned sol_lapack is the solution of the original
              !unscaled system of equations.
              call solve_linear_svd(verbosity, njacobian, 1, jacobian_lapack, naux, rhs1_lapack, naux,&
                   tolerance_lapack, rcond_lapack, sol_lapack, naux, info_lapack)
           else
              ! Use LU factorization solution from lapack to solve system of equations.
              ! n.b. scaling is done inside if necessary and jacobian_lapack and rhs1_lapack
              ! *may* be modified by these scaling factors.
              ! the returned sol_lapack is the solution of the original
              !unscaled system of equations.
              call dgesvx('E', 'N',&
                   njacobian, 1, jacobian_lapack, naux,&
                   lu_lapack, naux, ipiv_lapack, equed_lapack,&
                   row_lapack, col_lapack, rhs1_lapack, naux,&
                   sol_lapack, naux,&
                   rcond_lapack, ferr_lapack, berr_lapack,&
                   work_lapack, iwork_lapack, info_lapack)
           endif
           ! According to the dgesvx man page,
           ! "If the reciprocal of the condition number is less than machine precision,
           ! INFO = N+1 is returned as a warning, but the routine still goes on
           ! to solve for X and compute error bounds...."
           ! Therefore, the "(.not.if_svd.and.info_lapack.gt.njacobian)" clause
           ! in the if statement below is to ignore this warning if it does occur.
           if(.not.(info_lapack.eq.0.or.(.not.if_svd.and.info_lapack.gt.njacobian))) then
              if(verbosity.ge.1) then
                 write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
                 write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
                 write(stderr,*) 'eps ='
                 write(stderr,'(1p5e25.15e4)') eps
                 write(stderr,*) 'ioncount, ne, ne_old ='
                 write(stderr,'(i5,1p5e25.15e4)') ioncount, ne, ne_old
                 write(stderr,'(a,/,i5,1p5e25.15e4)') 'info_lapack, rcond_lapack = ', info_lapack, rcond_lapack
                 if(if_svd) then
                    write(stderr,*) 'free_eos_detailed ERROR: (second info_lapack) no solve_linear_svd solution'
                 else
                    write(stderr,*) 'free_eos_detailed ERROR: (second info_lapack) no dgesvx solution'
                 endif
              endif
              info = info_offset_lapack + info_lapack
              return
           endif
           if(verbosity.ge.2.and.rcond_lapack.lt.rcond_lapack_min) then
              write(stderr,'(a,/,i5,1p5e25.15e4)')&
                   'free_eos_detailed WARNING: (2) ioncount, rcond_lapack = ', ioncount, rcond_lapack
           elseif(.false..and.verbosity.ge.4) then
              write(stderr,'(a,/,i5,1p5e25.15e4)')&
                   'free_eos_detailed: (2) ioncount, rcond_lapack = ', ioncount, rcond_lapack
           endif
           do ijacobian = 1, njacobian
              if(ijacobian.gt.n_partial_aux) then
                 pac = real(sol_lapack(ijacobian,1),fp_kind)
              else
                 faux(ijacobian) = real(sol_lapack(ijacobian,1),fp_kind)
              endif
           enddo
        endif
        if(.false..and.verbosity.ge.4) then
           write(stderr,*) 'simple_lambda'
           write(stderr,'(5(1pe15.5e4,l2))') (simple_lambda(index_aux), simple_lambda(index_aux).gt.0._fp_kind,&
                index_aux = 1, njacobian)
           if(njacobian.gt.n_partial_aux) then
              write(stderr,*) 'Mixed NR and simple-iteration solution including fl'
              write(stderr,'(5(1pe15.5e4,l2))') (faux(index_aux), simple_lambda(index_aux).gt.0._fp_kind,&
                   index_aux = 1, n_partial_aux), pac, simple_lambda(njacobian).gt.0._fp_kind
           else
              write(stderr,*) 'Mixed NR and simple-iteration solution with frozen fl'
              write(stderr,'(5(1pe15.5e4,l2))') (faux(index_aux), simple_lambda(index_aux).gt.0._fp_kind,&
                   index_aux = 1, n_partial_aux)
           endif
        endif
        maxfaux_diag = 0._fp_kind
        do index_aux = 1, n_partial_aux
           ! 4th auxiliary variable is already in log form so the if
           ! statement selects all log variables and only if major
           ! auxiliary variable.
           iaux = partial_aux(index_aux)
           if(ifmajor(index_aux).and.&
                (iaux.eq.rl_aux_index.or.old_new_allow_log(index_aux).or.&
                old_new_allow_log_neg(index_aux)).and.&
                abs(faux(index_aux)).gt.abs(maxfaux_diag)) then
              ! n.b. maxfaux_diag perturbed by diagonal component and
              ! used for scaling, but maxfaux used for iteration control
              maxfaux_diag = faux(index_aux)
           endif
        enddo
        if(njacobian.gt.n_partial_aux.and.&
             ifmajor(njacobian).and.&
             abs(pac).gt.abs(maxfaux_diag))&
             maxfaux_diag = pac

        ! if change in sign, halve faux_limit_small (unless would force
        ! below faux_limit_min)
        if(maxfaux_diag*maxfaux_diag_old.lt.0._fp_kind.and.&
             ioncount.ge.2)&
             faux_limit_small =&
             max(faux_limit_min,0.5_fp_kind*faux_limit_small)

        ! only worry about maximum changes if there would have been
        ! a substantial change in ne proportional to n_e/rho for the
        ! previous *unscaled* solution.  ne is an overall
        ! measure of ionization fractions.  (We also control rho for those
        ! cases [molecular formation important] where it is more independent
        ! of ne.) Also, lambda.  This logic is meant to
        ! take care of the case where there are large ln changes in an
        ! auxiliary variable that have little or no effect on overall
        ! ionization fractions because the auxiliary variable is
        ! approaching zero.
        if(max(abs(rl-rl_old),abs(ne-ne_old)/ne,abs(lambda-lambda_old)/max(1.e-15_fp_kind,lambda))/faux_scale.le.1.e-2_fp_kind) then
           faux_limit = faux_limit_large
           faux_limit_small = faux_limit_small_start
        else
           faux_limit = faux_limit_small
        endif
        faux_scale = min(1._fp_kind,&
             faux_limit/max(1.e-16_fp_kind,abs(maxfaux_diag)))
        if(njacobian.gt.n_partial_aux) then
           ! scale so that maximum change in fl is always <= 0.01
           faux_scale = min(faux_scale,&
                min(faux_limit,0.01_fp_kind)/max(1.e-16_fp_kind,abs(pac)))
        endif
        ! cannot get into much trouble if current fl and NR-predicted fl
        ! are less than -10 so turn off limits to NR change in that case.
        if(njacobian.gt.n_partial_aux) then
           if(max(fl,fl+pac).lt.-10._fp_kind) faux_scale = 1._fp_kind
        else
           ! fixed fl.
           if(fl.lt.-10._fp_kind) faux_scale = 1._fp_kind
        endif
        !temporary (no limit on NR change)
        !faux_scale = 1.d0
        do index_aux = 1, n_partial_aux
           iaux = partial_aux(index_aux)
           ! do nothing when eos_calc returns a zero aux (index.ne.rl_aux_index)
           ! because solution may add some noise.
           if(old_new_allow_log(index_aux).or.&
                old_new_allow_log_neg(index_aux).or.&
                iaux.eq.rl_aux_index) then
              ! these comments pertain to each of branches below:
              ! keep direction of faux vector, but scale it so no
              ! log component exceeds faux_limit.
              faux(index_aux) = faux_scale*faux(index_aux)
              if(old_new_allow_log(index_aux)) then
                 aux(iaux) = log(aux_old(iaux)) + faux(index_aux)
                 if(aux(iaux).gt.ln_aux_underflow) then
                    aux(iaux) = exp(aux(iaux))
                 else
                    aux(iaux) = 0._fp_kind
                 endif
              elseif(old_new_allow_log_neg(index_aux)) then
                 aux(iaux) = log(-aux_old(iaux)) + faux(index_aux)
                 if(aux(iaux).gt.ln_aux_underflow) then
                    aux(iaux) = -exp(aux(iaux))
                 else
                    aux(iaux) = 0._fp_kind
                 endif
              elseif(iaux.eq.rl_aux_index) then
                 aux(iaux) = aux_old(iaux) + faux(index_aux)
                 ! n.b. last two branches are old code that is avoided
                 ! by above outer if statement.
              elseif(sign(1._fp_kind,aux(iaux)).eq.sign(1._fp_kind, aux_old(iaux) + faux(index_aux))) then
                 ! by logic of code this branch occurs only if sign change
                 ! or change from zero to non-zero
                 ! only use Newton-Raphson value if its new sign agrees with
                 ! sign of value delivered by eos_calc.
                 aux(iaux) = aux_old(iaux) + faux(index_aux)
              elseif(aux_old(iaux).ne.0._fp_kind) then
                 ! scale change which involves sign change that disagrees with
                 ! NR sign
                 aux(iaux) = aux_old(iaux) + faux_scale*(aux(iaux)-aux_old(iaux))
              endif
           endif
           ! n.b. above logic drops through (accepts aux delivered by
           ! eos_calc) in case aux_old or aux is zero or
           ! opposite signs (aside from iaux.eq.rl_aux_index which is logarithmic).

           ! if trying to go from non-zero to zero or vice versa
           ! use special limits on aux
           if(iaux.ne.rl_aux_index.and.&
                (aux_old(iaux).ne.aux(iaux).and.&
                (aux_old(iaux).eq.0._fp_kind.or.&
                aux(iaux).eq.0._fp_kind))) then
              ! continue with iteration unless zerolim is zero
              if(zerolim(iaux).gt.0._fp_kind) iflast = 0
              if(aux_old(iaux).ne.0._fp_kind) then
                 ! if trying to go from non-zero to zero...
                 aux(iaux) = aux_old(iaux)/exp(faux_limit)
                 ! define maximum magnitude for subsequent zero to non-zero
                 zerolim(iaux) = abs(aux(iaux))
                 ! mark infinite decrease in log auxiliary variable.
                 maxfaux = -50._fp_kind
                 index_max = index_aux
                 ! if this occurred on previous iteration, increase
                 ! ifzerocount and allow to go to zero if 5 in a row
                 if(ioncount-1.eq.ioncountzero) then
                    ! zero occurred on previous iteration
                    ! n.b. this branch only once per iteration
                    ifzerocount = ifzerocount + 1
                 elseif(ioncount.ne.ioncountzero) then
                    ! zero did not occur for previous iteration.
                    ! n.b. this branch only once per iteration
                    ifzerocount = 1
                 endif
                 ioncountzero = ioncount
                 if(mod(ifzerocount,5).eq.0) aux(iaux) = 0._fp_kind
              else
                 ! if going from zero to non-zero, don't go above
                 ! previous magnitude found on non-zero to zero step
                 ! or if no previous such step use initial zerolim value from above.
                 if(abs(aux(iaux)).gt.zerolim(iaux)) aux(iaux) = sign(zerolim(iaux), aux(iaux))
                 ! mark infinite increase in log auxiliary variable.
                 maxfaux = 50._fp_kind
                 index_max = index_aux
              endif
           endif
        enddo
        if(njacobian.gt.n_partial_aux) then
           ! scale fl change just like raw NR changes were
           ! scaled above.
           pab = faux_scale*pac
           fl = fl + pab
        endif
        ! rho must be updated also.
        rho = exp(rl)
        if(.false..and.verbosity.ge.4) then
           if(njacobian.gt.n_partial_aux) then
              write(stderr,*) 'old value of aux including fl'
              write(stderr,'(5(1pe15.5e4,l2))')&
                   (aux_old(partial_aux(index_aux)), ifmajor(index_aux), index_aux = 1, n_partial_aux), flold, ifmajor(njacobian)
              write(stderr,*) 'new value of aux including fl'
              write(stderr,'(5(1pe15.5e4,l2))')&
                   (aux(partial_aux(index_aux)), ifmajor(index_aux), index_aux = 1, n_partial_aux), fl, ifmajor(njacobian)
           else
              write(stderr,*) 'old value of aux for frozen fl'
              write(stderr,'(5(1pe15.5e4,l2))') (aux_old(partial_aux(index_aux)), ifmajor(index_aux), index_aux = 1, n_partial_aux)
              write(stderr,*) 'new value of aux for frozen fl'
              write(stderr,'(5(1pe15.5e4,l2))') (aux(partial_aux(index_aux)), ifmajor(index_aux), index_aux = 1, n_partial_aux)
           endif
        endif
        if(njacobian.gt.n_partial_aux.and.ieee_support_nan(1._fp_kind)) then
           if(ieee_is_nan(fl)) then
              if(verbosity.ge.1) then
                 pab = faux_scale*pac
                 write(stderr,*) 'match_variable/ln(10), tl/ln(10) ='
                 write(stderr,'(1p5e25.15e4)') match_variable/ln10, tl/ln10
                 write(stderr,*) 'faux_scale, pac, pab = faux_scale*pac, fl = ', faux_scale, pac, pab, fl
                 write(stderr,*) 'free_eos_detailed ERROR: (4) NaN value of fl detected'
              endif
              info = info_offset_free_eos_detailed + 4
              return
           endif
        endif
        if((verbosity.ge.2.and.ioncount.ge.100).or.verbosity.ge.4) then
           write(stderr,*) 'ioncount, index_max, maxfaux, maxfaux_diag, faux_limit, faux_scale, fl'
           write(stderr,*) ioncount, index_max, maxfaux, maxfaux_diag, faux_limit, faux_scale, fl
        endif
     enddo  ! do while(iflast.ne.1.and.ioncount.lt.maxioncount)
     if((ifwarm.and.ioncount.ge.maxioncount).or.ioncount.ge.maxioncount) then
        if(verbosity.ge.1) then
           write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
           write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
           write(stderr,*) 'eps ='
           write(stderr,'(1p5e25.15e4)') eps
           write(stderr,*) 'ioncount, ne, ne_old ='
           write(stderr,'(i5,1p5e25.15e4)') ioncount, ne, ne_old
           write(stderr,*) 'maxfaux ='
           write(stderr,'(1pe25.15e4)') maxfaux
           write(stderr,*) 'free_eos_detailed ERROR: (5) auxiliary variable iteration did not converge'
        endif
        info = info_offset_free_eos_detailed + 5
        return
     endif
     !temporary just to get comparisons, but messes up derivatives.
     if(.false.) then
        if(verbosity.ge.3) write(stderr,*) 'start final eos_bfgs per NR cycle'
        call eos_bfgs(&
             verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
             ifrad, match_variable, kif, fl,&
             aux_old, aux, auxf, auxt, aux_dv,&
             njacobian,&
             rhs1, jacobian, p, pr,&
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
             partial_elements(1:n_partial_elements+2), n_partial_elements, ion_end, mion_end,&
             ifionized, if_pteh, if_mc, ifreducedmass,&
             ifsame_abundances, ifmtrace,&
             iatomic_number, ifpi_local,&
             ifpl, ifmodified, ifh2, ifh2plus,&
             izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax,&
             eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
             r_ion3, r_neutral,&
             ifelement, dvzero, dv, dvf, dvt, nion,&
             ne, nef, net, sion, sionf, siont, uion,&
             pnorad, pnoradf, pnoradt, pnorad_aux(1:n_partial_aux+2),&
             free_fp, free_aux(1:n_partial_aux+2),&
             nextrasum,&
             sumpl1, sumpl1f, sumpl1t, sumpl2,&
             h2, h2f, h2t, h2_dv, info)
        if(info.ne.0) return
        if(verbosity.ge.3) write(stderr,*) 'end final eos_bfgs per NR cycle'
     endif
     ! save aux so that can restore after eos_calc below.
     ! (It is possible that eos_calc could move off of
     ! NR-converged value by propagating numerical
     ! noise.)
     do index_aux = 1,n_partial_aux
        iaux = partial_aux(index_aux)
        aux_old(iaux) = aux(iaux)
     enddo
     ! from partial of RHS wrt match_variable and tl use chain rule
     ! and pre-factored Jacobian to calculate partials of
     ! converged auxiliary variables and fl wrt match_variable and tl.
     do index_aux = 1, n_partial_aux
        iaux = partial_aux(index_aux)
        if(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)) then
           ! transform to log form of faux
           if(njacobian.gt.n_partial_aux) then
              rhs2(index_aux,1) = 0._fp_kind
           else
              rhs2(index_aux,1) = auxf(iaux)/aux_old(iaux)
           endif
           rhs2(index_aux,2) = auxt(iaux)/aux_old(iaux)
        else
           if(njacobian.gt.n_partial_aux) then
              rhs2(index_aux,1) = 0._fp_kind
           else
              rhs2(index_aux,1) = auxf(iaux)
           endif
           rhs2(index_aux,2) = auxt(iaux)
        endif
     enddo
     if(njacobian.gt.n_partial_aux) then
        if(kif.eq.1) then
           if(ifrad.gt.1) then
              ! For kif = 1, RHS = match_variable - log(pnorad)
              ! partial of RHS wrt tl
              rhs2(index_aux,2) = -pnoradt/pnorad
           else
              ! For kif = 1, RHS = match_variable - log(p)
              rhs2(index_aux,2) = -(4._fp_kind*pr + pnoradt)/p
           endif
        elseif(kif.eq.2) then
           ! For kif = 2, RHS = match_variable - rl
           rhs2(index_aux,2) = -rt
        endif
        ! partial of RHS wrt to match_variable.
        rhs2(index_aux,1) = 1._fp_kind
     endif
     rhs2_lapack(1:njacobian,1:2) = real(rhs2(1:njacobian,1:2),lapack_fp_kind)
     if(if_svd) then
        ! SVD solution with jacobian already
        ! scaled (potentially) and decomposed to lu_lapack???
        ! n.b. equed_lapack keeps track of what scaling occurred on
        ! previous call that decomposed the (scaled) jacobian and this
        ! scaling if any is applied to rhs2
        ! the returned sol_lapack is the solution of the original
        ! unscaled system of equations.
        call solve_linear_svd(verbosity, njacobian, 2, jacobian_lapack, naux, rhs2_lapack, naux,&
             tolerance_lapack, rcond_lapack, sol_lapack, naux, info_lapack )
     else
        ! LU factorization solution from lapack with jacobian already
        ! scaled (potentially) and factored to lu_lapack
        ! n.b. equed_lapack keeps track of what scaling occurred on
        ! previous call that factored (scaled) jacobian and this
        ! scaling if any is applied to rhs_lapack
        ! the returned sol_lapack is the solution of the original
        ! unscaled system of equations.
        call dgesvx('F', 'N',&
             njacobian, 2, jacobian_lapack, naux,&
             lu_lapack, naux, ipiv_lapack, equed_lapack,&
             row_lapack, col_lapack, rhs2_lapack, naux,&
             sol_lapack, naux,&
             rcond_lapack, ferr_lapack, berr_lapack,&
             work_lapack, iwork_lapack, info_lapack)
     endif
     ! According to the dgesvx man page,
     ! "If the reciprocal of the condition number is less than machine precision,
     ! INFO = N+1 is returned as a warning, but the routine still goes on
     ! to solve for X and compute error bounds...."
     ! Therefore, the "(.not.if_svd.and.info_lapack.gt.njacobian)" clause
     ! in the if statement below is to ignore this warning if it does occur.
     if(.not.(info_lapack.eq.0.or.(.not.if_svd.and.info_lapack.gt.njacobian))) then
        if(verbosity.ge.1) then
           write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
           write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
           write(stderr,*) 'eps ='
           write(stderr,'(1p5e25.15e4)') eps
           write(stderr,*) 'ioncount, ne, ne_old ='
           write(stderr,'(i5,1p5e25.15e4)') ioncount, ne, ne_old
           write(stderr,'(a,/,i5,1p5e25.15e4)') 'info_lapack, rcond_lapack = ', info_lapack, rcond_lapack
           if(if_svd) then
              write(stderr,*) 'free_eos_detailed ERROR: (third info_lapack) no solve_linear_svd solution'
           else
              write(stderr,*) 'free_eos_detailed ERROR: (third_info_lapack) no dgesvx solution'
           endif
        endif
        info = info_offset_lapack + info_lapack
        return
     endif
     if(verbosity.ge.2.and.rcond_lapack.lt.rcond_lapack_min) then
        write(stderr,'(a,/,i5,1p5e25.15e4)')&
             'free_eos_detailed WARNING: (3) ioncount, rcond_lapack = ', ioncount, rcond_lapack
     elseif(.false..and.verbosity.ge.4) then
        write(stderr,'(a,/,i5,1p5e25.15e4)')&
             'free_eos_detailed: (3) ioncount, rcond_lapack = ', ioncount, rcond_lapack
     endif
     do index_aux = 1, n_partial_aux
        iaux = partial_aux(index_aux)
        if(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)) then
           ! transform from log auxiliary variable derivatives
           auxf(iaux) = real(sol_lapack(index_aux,1),fp_kind)*aux_old(iaux)
           auxt(iaux) = real(sol_lapack(index_aux,2),fp_kind)*aux_old(iaux)
        else
           auxf(iaux) = real(sol_lapack(index_aux,1),fp_kind)
           auxt(iaux) = real(sol_lapack(index_aux,2),fp_kind)
        endif
     enddo
     if(njacobian.gt.n_partial_aux) then
        ! transform independent variable from match_variable to fl.
        ! this could introduce significance loss when we transform
        ! back again below, but this does not appear to be a problem.
        ! sol_lapack(njacobian,1) is the partial of
        ! fl(match_variable,tl) wrt match_variable
        ! sol_lapack(njacobian,2) is the partial of
        ! fl(match_variable,tl) wrt tl.
        ! partial of calculated match_variable(fl,tl) wrt fl
        match_variablef = 1._fp_kind/real(sol_lapack(njacobian,1),fp_kind)
        ! partial of calculated match_variable(fl,tl) wrt tl
        match_variablet = -match_variablef*real(sol_lapack(njacobian,2),fp_kind)
        do index_aux = 1, n_partial_aux
           iaux = partial_aux(index_aux)
           ! auxf and auxt contain partial of aux wrt match_variable and tl.
           ! so transform appropriately to fl and tl derivatives.
           auxt(iaux) = auxf(iaux)*match_variablet + auxt(iaux)
           auxf(iaux) = auxf(iaux)*match_variablef
        enddo
     endif
     if(.false..and.verbosity.ge.4) then
        write(stderr,*) 'auxf = '
        write(stderr,'(5(1pe15.5e4,2x))') (auxf(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
        write(stderr,*) 'auxt = '
        write(stderr,'(5(1pe15.5e4,2x))') (auxt(partial_aux(index_aux)), index_aux = 1, n_partial_aux)
     endif
     !*************************************************************
     ! now that all auxiliary variable derivatives have
     ! been determined from NR procedure, call whole
     ! eos procedure again to determine fl and tl derivatives
     ! of all thermodynamic quantities (tq).
     ! n.b.  eos_tqft returns nu form of many auxiliary variables
     ! since this form is required by the entropy calculation.
     iteration_count = iteration_count + 1
     call eos_tqft(&
          verbosity,&
          dv_aux,&
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
          partial_elements(1:n_partial_elements+2), ion_end, mion_end,&
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
          extrasum(:nextrasum), extrasumf(:nextrasum), extrasumt(:nextrasum), extrasum_dv(:,:nextrasum),&
          sumpl1, sumpl1f, sumpl1t, sumpl2,&
          rl, rf, rt, r_dv,&
          h2, h2f, h2t, h2plus, h2plusf, h2plust, h2plus_dv,&
          xextrasum, xextrasumf, xextrasumt, xextrasum_dv, info)
     if(info.ne.0) return

     ! do not restore aux, since part of it has been put
     ! in nu form by eos_calc.
     !*************************************************************
     ! end of if block to calculate derivatives of all
     ! thermodynamic quantities wrt fl and tl using the chain rule
     ! for NR-converged results.
  endif !if(.not.ifnoaux_iteration)
  ! there is the possibility of pressure ionization and
  ! Planck-Larkin terms.
  ! n.b. For kif = 1, some of the pressure-related stuff (but not all of
  ! it in the ifnoaux_iteration case) is calculated previously.  But just
  ! repeat those calculations here rather than worrying about special logic.
  if(ifpi_local.eq.0) then
     ! there are no p, s, u terms from pressure ionization.
     ppi = 0._fp_kind
     ppif = 0._fp_kind
     ppit = 0._fp_kind
     spi = 0._fp_kind
     spif = 0._fp_kind
     spit = 0._fp_kind
     upi = 0._fp_kind
  elseif(ifpi_local.eq.1) then
     ! According to debug printouts, ifnr can be uninitialized for kif = 0
     ! Fixme, review the kif = 0 logic to verify this debug result.
     ! Must calculate f and t derivatives.
     ifnr = 0
     ifnr03 = .true.
     call pteh_pi_end(&
          ifnr, full_sum1, rho, rf, rt, t, ne, nef, net,&
          ppi, ppif, ppit, spi, spif, spit, upi)
  elseif(ifpi_local.eq.2) then
     ! n.b. call to eos_tqft produces nu = n/(rho*avogadro)
     ! form of h_ion, he_ion, and he_ion2.  Also use nu form
     ! of nx, ny, nz, ne, etc.
     call fjs_pi_end(&
          t, rho, rf, rt,&
          nux, 0._fp_kind, 0._fp_kind, nuy, 0._fp_kind, 0._fp_kind, nuz, 0._fp_kind, 0._fp_kind,&
          h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t, ne, nef, net,&
          ppi, ppif, ppit, spi, spif, spit, upi)
  elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
     ! n.b. call to eos_tqft produces nu = n/(rho*avogadro)
     ! form of extrasum.
     call mdh_pi_end(&
          t, rho, rf, rt,&
          nion,&
          extrasum(:nextrasum), extrasumf(:nextrasum), extrasumt(:nextrasum),&
          ppi, ppif, ppit, spi, spif, spit, upi)
  endif
  if(ifexcited.gt.0) then
     call excitation_pi_end(t, rho, rf, rt,&
          pexcited, pexcitedf, pexcitedt,&
          sexcited, sexcitedf, sexcitedt, uexcited)
  else
     pexcited = 0._fp_kind
     pexcitedf = 0._fp_kind
     pexcitedt = 0._fp_kind
     sexcited = 0._fp_kind
     sexcitedf = 0._fp_kind
     sexcitedt = 0._fp_kind
     uexcited = 0._fp_kind
  endif
  if(ifpl.eq.1) then
     ! no pressure terms for Planck-Larkin occupation probability
     ! but there are entropy and energy terms.
     spi = spi + cr*sumpl1
     upi = upi + cr*t*sumpl2
     if(ifnr03) then
        spif = spif + cr*sumpl1f
        spit = spit + cr*sumpl1t
     endif
  endif
  if(if_pteh.eq.1) then
     ! PTEH (full ionization) approximation to sum0, sum2
     ! reassert this approximation because eos_calc
     ! calls ionize which messes a bit with sum0 and sum2
     sum0ne = full_sum0/full_sum1
     sum0 = sum0ne*n_e
     sum0f = 0._fp_kind
     sum0t = 0._fp_kind
     sum2ne = full_sum2/full_sum1
     sum2 = sum2ne*n_e
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
     ! n.b. call to eos_tqft produces nu = n/(rho*avogadro)
     ! form of h2plus, h_ion, he_ion, he_ion2.
     sum0 = (sum0_mc + h2plus + h_ion*(1._fp_kind+hcon_mc) +&
          he_ion +  he_ion2)*(rho*avogadro)
     sum2 = (sum2_mc + h2plus + h_ion*(1._fp_kind+hcon_mc) +&
          he_ion +  he_ion2*(4._fp_kind+hecon_mc))*(rho*avogadro)
     sum0f = (h2plusf + h_ionf*(1._fp_kind+hcon_mc) +&
          he_ionf + he_ion2f)*(rho*avogadro) + sum0*rf
     sum0t = (h2plust + h_iont*(1._fp_kind+hcon_mc) +&
          he_iont + he_ion2t)*(rho*avogadro) + sum0*rt
     sum2f = (h2plusf + h_ionf*(1._fp_kind+hcon_mc) +&
          he_ionf + he_ion2f*(4._fp_kind+hecon_mc))*(rho*avogadro) + sum2*rf
     sum2t = (h2plust + h_iont*(1._fp_kind+hcon_mc) +&
          he_iont + he_ion2t*(4._fp_kind+hecon_mc))*(rho*avogadro) + sum2*rt
  endif
  call master_coulomb_end(rhostar,&
       sum0, sum0ne, sum0f, sum0t,&
       sum2, sum2ne, sum2f, sum2t,&
       n_e, t, pstar,&
       ifcoulomb_mod, if_dc, if_pteh,&
       dpcoulomb, dpcoulombf, dpcoulombt,&
       dscoulomb, dscoulombf, dscoulombt, ducoulomb)
  call exchange_end(rhostar, pstar,&
       pex, pext, pexf,&
       sex, sexf, sext, uex)
  p0 = cr*rho*t
  ni = full_sum0 - (h2+h2plus)
  pion = ni*p0
  pionf = ni*p0*rf - (h2f+h2plusf)*p0
  pe_cgs = cpe*pe
  pnorad = pe_cgs + pion + ppi + pexcited + dpcoulomb + pex
  p = pnorad + pr
  if(p.le.0._fp_kind) then
     if(verbosity.ge.1) then
        write(stderr,*) 'match_variable/ln(10), fl, tl/ln(10) ='
        write(stderr,'(1p5e25.15e4)') match_variable/ln10, fl, tl/ln10
        write(stderr,*) 'eps ='
        write(stderr,'(1p5e25.15e4)') eps
        write(stderr,*) 'ioncount, ne, ne_old ='
        write(stderr,'(i5,1p5e25.15e4)') ioncount, ne, ne_old
        write(stderr,*) 'free_eos_detailed ERROR: (6) negative p calculated'
     endif
     info = info_offset_free_eos_detailed + 6
     return
  endif
  ! d ln p(f,t)/d ln f
  pnoradf = (pe_cgs*pef + pionf + ppif + pexcitedf +&
       dpcoulombf + pexf)/pnorad
  pf = pnoradf*pnorad/p
  pl = log(p)
  ! pr_ratio, ppie_ratio, pc_ratio, and pex_ratio  calculated only for
  ! diagnostic and teaching purposes.
  pr_ratio = pr/p
  ! Combine pressure-ionization correction to lower states and entire
  ! effect of pressure-ionization corrected Rydberg states since the
  ! two effects tend to offset each other and should be treated together
  ! (or not).
  ppie_ratio = (ppi + pexcited)/p
  pc_ratio = dpcoulomb/p
  pex_ratio = pex/p

  ! pion = ni*p0
  piont = ni*p0*(1._fp_kind + rt) - (h2t+h2plust)*p0
  ! d ln p(f,T)/d ln t
  pnoradt = (pe_cgs*pet + piont + ppit + pexcitedt +&
       dpcoulombt + pext)/pnorad
  pt = (pnorad*pnoradt + 4._fp_kind*pr)/p
  ! rpt = d ln rho(T,P)/ d ln P = 1/chi_rho
  rpt = rf/pf
  ! rtp = - d ln rho(T,P)/ d ln T  n.b. *negative* sign = chi_t/chi_rho
  rtp = rpt*pt - rt
  rlout = rl
  ! Save determined value of fl for next call to this routine.
  flold = fl
  degeneracy(1) = fl
  if(kif.eq.0) then
     degeneracy(2) = 1._fp_kind
     degeneracy(3) = 0._fp_kind
  elseif(kif.eq.1) then
     fm = 1._fp_kind/pf
     ft = -pt/pf
     degeneracy(2) = fm
     degeneracy(3) = ft
     if(ifrad.eq.2) then
        ! kif=1, and ifrad = 2 returns most derivatives wrt ln p, but in the
        ! special case of fm and ft return derivatives wrt ln pnorad
        ! since fm and ft used outside in calling free_eos programme
        ! for Taylor series.
        fm = 1._fp_kind/pnoradf
        ft = -pnoradt/pnoradf
     endif
  elseif(kif.eq.2) then
     fm = 1._fp_kind/rf
     ft = -rt/rf
     degeneracy(2) = fm
     degeneracy(3) = ft
  endif
  pressure(1) = pl
  if(kif.eq.0) then
     pressure(2) = pf
     pressure(3) = pt
     ! special tests
     ! pressure(1) = pe_cgs
     ! pressure(2) = pef*pe_cgs
     ! pressure(3) = pet*pe_cgs
     ! pressure(1) = pion
     ! pressure(2) = pionf
     ! pressure(3) = piont
     ! pressure(1) = pexcited
     ! pressure(2) = pexcitedf
     ! pressure(3) = pexcitedt
     ! index_aux = 12
     ! iaux = partial_aux(index_aux)
     ! pressure(1) = aux(iaux)
     ! pressure(2) = auxf(iaux)
     ! pressure(3) = auxt(iaux)
     ! pressure(1) = dpcoulomb
     ! pressure(2) = dpcoulombf
     ! pressure(3) = dpcoulombt
     ! pressure(1) = pex
     ! pressure(2) = pexf
     ! pressure(3) = pext
  elseif(kif.eq.1) then
     pressure(2) = 1._fp_kind
     pressure(3) = 0._fp_kind
  elseif(kif.eq.2) then
     if(ifrad.gt.1) then
        ! n.b. this makes derivatives of this quantity
        ! inconsistent, but nevertheless some tables have
        ! all quantities calculated with radiation
        ! pressure except for pg = p - pr.
        pressure(1) = log(pnorad)
     else
        pressure(1) = log(p)
     endif
     pressure(2) = 1._fp_kind/rpt
     pressure(3) = rtp/rpt
  endif
  density(1) = rl
  if(kif.eq.0) then
     density(2) = rf
     density(3) = rt
     ! temporary check
     ! density(1) = ne
     ! density(2) = nef
     ! density(3) = net
  elseif(kif.eq.1) then
     density(2) = rpt
     density(3) = -rtp
  elseif(kif.eq.2) then
     density(2) = 1._fp_kind
     density(3) = 0._fp_kind
  endif
  sr = 4._fp_kind*pr/p0
  ! the third term corrects for the dropped ideal terms in
  ! ionize and eos_calc.  (See the notes for those routines).
  s = ne*se + sr + full_sum0*(2.5_fp_kind + 1.5_fp_kind*tl - rl) +  sion
  s = cr*s + spi + sexcited + dscoulomb/rho + sex/rho
  sf = nef*se + ne*se*sef - rf*sr - ni*rf + sionf
  sf = cr*sf + spif + sexcitedf +&
       (dscoulombf + sexf - (dscoulomb+sex)*rf)/rho
  st = net*se + ne*se*set + (3._fp_kind-rt)*sr + ni*(1.5_fp_kind-rt) + siont
  st = cr*st + spit + sexcitedt +&
       (dscoulombt + sext - (dscoulomb+sex)*rt)/rho
  ! by definition of d S(T,P)/d ln T and by chain rule
  cp = st - sf*(pt/pf)
  entropy(1) = s
  if(kif.eq.0) then
     if(iftc.eq.1) then
        ! this derivatives from thermodynamic consistency arguments.
        ! partial s(p,t) wrt ln p = (p/(rho*t)*d ln rho(P,T)/d ln T
        entropy(2) = -pf*(p/(rho*t))*rtp
     else
        ! straight derivative by chain rule
        entropy(2) = sf
     endif
     entropy(3) = st
     ! special tests
     ! entropy(1) = cr*ne*se
     ! entropy(2) = cr*(nef*se + ne*se*sef)
     ! entropy(3) = cr*(net*se + ne*se*set)
     ! ideal part of entropy
     ! entropy(1) = cr*(eps(19)*(2.5d0 + 1.5d0*tl - rl) +  sion)
     ! entropy(2) = cr*(- eps(19)*rf + sionf)
     ! entropy(3) = cr*(eps(19)*(1.5d0-rt) + siont)
     ! entropy(1) = cr*(full_sum0*(2.5d0 + 1.5d0*tl - rl) +  sion)
     ! entropy(2) = cr*(-ni*rf + sionf)
     ! entropy(3) = cr*(ni*(1.5d0-rt) + siont)
     ! entropy(1) = sexcited
     ! entropy(2) = sexcitedf
     ! entropy(3) = sexcitedt
     ! entropy(1) = dv(3) + dvzero(2)
     ! entropy(2) = dvf(3)
     ! entropy(3) = dvt(3)
     ! index_aux = 13
     ! iaux = partial_aux(index_aux)
     ! entropy(1) = aux(iaux)
     ! entropy(2) = auxf(iaux)
     ! entropy(3) = auxt(iaux)
     ! entropy(1) = dscoulomb/rho
     ! entropy(2) = dscoulombf/rho - dscoulomb*rf/rho
     ! entropy(3) = dscoulombt/rho - dscoulomb*rt/rho
     ! entropy(1) = sex/rho
     ! entropy(2) = sexf/rho - sex*rf/rho
     ! entropy(3) = sext/rho - sex*rt/rho
  elseif(kif.eq.1) then
     if(iftc.eq.1) then
        ! this derivatives from thermodynamic consistency arguments.
        ! partial s(p,t) wrt ln p = (p/(rho*t)*d ln rho(P,T)/d ln T
        entropy(2) = -(p/(rho*t))*rtp
     else
        ! straight derivative by chain rule
        entropy(2) = sf/pf
     endif
     ! partial s(p,t) wrt ln t = C_p by definition (see above).
     entropy(3) = cp
  elseif(kif.eq.2) then
     if(iftc.eq.1) then
        ! this derivatives from thermodynamic consistency arguments.
        ! partial s(rho,t) wrt ln rho
        entropy(2) = -(p/(rho*t))*rtp/rpt
     else
        ! straight derivative by chain rule
        entropy(2) = sf/rf
     endif
     ! straight derivatives by chain rule
     ! partial s(rho,t) wrt ln t
     entropy(3) = st - sf*rt/rf
     ! this derivative from thermodynamic consistency arguments.
     ! however, it is straightforward to derive this from straight
     ! derivatives and thermodynamic consistency for entropy(2),
     ! so it is not an independent test.
     ! entropy(3) = cp - rtp*rtp*p/rpt/rho/t
  endif
  ! vdb expression with metal, h2, and h2+, and coulomb contribution added.
  ! note energy of electrons/volume = n_e k T ue, but n_e = ne*rho/H, thus,
  ! energy of electrons/volume = rho * ne R T ue, and
  ! energy/mass = ne R T ue.
  u = cr*t*(ue*ne+1.5_fp_kind*ni+0.75_fp_kind*sr) + (c2*cr)*uion + upi + uexcited + (ducoulomb+uex)/rho
  ! for hydrogen shift energy zero from monatomic to diatomic ground state,
  ! that is add h2diss/2 to monatomics and h2diss to diatomics
  ! n.b. this zero point shift does not affect pressure, entropy, or
  ! chemical potentials, but does keep ideal u positive.
  u = u + 0.5_fp_kind*(c2*cr)*h2diss*eps(1)
  !temporary
  free_rad = -cr*t*sr/4._fp_kind
  free_e = cr*t*(ue*ne-ne*se)
  free_coulomb = (ducoulomb -t*dscoulomb)/rho
  free_ex = (uex -t*sex)/rho
  if(ifpl.eq.1) then
     free_pl = cr*t*(sumpl2 - sumpl1)
  else
     free_pl = 0._fp_kind
  endif
  free_pi = upi -t*spi - free_pl
  free_excited = uexcited -t*sexcited
  free_ion = cr*t*(1.5_fp_kind*ni - full_sum0*&
       (2.5_fp_kind + 1.5_fp_kind*tl - rl) - sion) + (c2*cr)*uion +&
       0.5_fp_kind*(c2*cr)*h2diss*eps(1)
  free =&
       free_rad + free_e + free_coulomb + free_ex +&
       free_pi +&
       free_excited + free_pl + free_ion
  if(.false..and.verbosity.ge.4) then
     write(stderr,'(a,/,(1p2e25.15e4))')&
          'free_eos_detailed: fl, tl, pnorad, (scaled) free = ',&
          fl, tl, pnorad, free/(full_sum0*cr*t)
  endif
  free_rad = free_rad/free
  free_e = free_e/free
  free_coulomb = free_coulomb/free
  free_ex = free_ex/free
  free_pi = free_pi/free
  free_excited = free_excited/free
  free_pl = free_pl/free
  free_ion = free_ion/free
  if(.false..and.verbosity.ge.4) then
     write(stderr,*) 'fractions of rad, e, coulomb, ex, pi, excited, pl, and ion components of free = '
     write(stderr,*) free_rad, free_e, free_coulomb, free_ex, free_pi, free_excited, free_pl, free_ion
     write(stderr,*) 'total rad, e, coulomb, ex, pi, star, pl, and ion components of free = '
     write(stderr,'(1p4e25.15e4)')&
          free*free_rad, free*free_e, free*free_coulomb,&
          free*free_ex, free*free_pi, free*free_excited,&
          free*free_pl, free*free_ion
  endif
  ! adiabiatic gradient d log T/d log P
  grada = -sf/(cp*pf)
  ! do in this way to avoid overflow problems
  vz = (rt/cp)*(sf/pf)
  vx = (rf/cp)*(st/pf)
  gamma1 = 1._fp_kind/(vx-vz)
  if(kif.eq.0) then
     energy(1) = u
     ! these derivatives from thermodynamic consistency arguments.
     energy(2) = pf*(rpt - rtp)*p/rho
     ! it is tedious but straightforward to show this expression has
     ! large significance loss in the radiation dominated case so
     ! don't depend on it in that case.
     energy(3) = (pt*(rpt - rtp)*p/rho - rtp*p/rho) + cp*t
     ! special tests
     ! energy(1) = h2
     ! energy(2) = h2f
     ! energy(3) = h2t
     ! energy(1) = dv(1) + dvzero(1)
     ! energy(2) = dvf(1)
     ! energy(3) = dvt(1)
     ! index_aux = 14
     ! iaux = partial_aux(index_aux)
     ! energy(1) = aux(iaux)
     ! energy(2) = auxf(iaux)
     ! energy(3) = auxt(iaux)
  elseif(kif.eq.1) then
     energy(1) = u
     ! these derivatives from thermodynamic consistency arguments.
     energy(2) = (rpt - rtp)*p/rho
     energy(3) = cp*t - rtp*p/rho
  elseif(kif.eq.2) then
     energy(1) = u
     ! these derivatives from thermodynamic consistency arguments.
     energy(2) = (rpt - rtp)*p/rho/rpt
     energy(3) = cp*t - rtp*rtp*p/rpt/rho
  endif
  if(kif.eq.0) then
     enthalpy(1) = u + p/rho
     ! these derivatives from thermodynamic consistency arguments.
     enthalpy(2) = pf*(1._fp_kind-rtp)*p/rho
     ! it is tedious but straightforward to show this expression has
     ! large significance loss in the radiation dominated case so
     ! don't depend on it in that case.
     enthalpy(3) = pt*(1._fp_kind-rtp)*p/rho + cp*t
     ! special tests
     ! enthalpy(1) = dve + eta
     ! enthalpy(2) = dvef + wf
     ! enthalpy(3) = dvet
     ! enthalpy(1) = dv(2) + dvzero(2)
     ! enthalpy(2) = dvf(2)
     ! enthalpy(3) = dvt(2)
     ! index_aux = 15
     ! iaux = partial_aux(index_aux)
     ! enthalpy(1) = aux(iaux)
     ! enthalpy(2) = auxf(iaux)
     ! enthalpy(3) = auxt(iaux)
     ! enthalpy(1) = h2plus
     ! enthalpy(2) = h2plusf
     ! enthalpy(3) = h2plust
     ! enthalpy(1) = dve_pi
     ! enthalpy(2) = dve_pif
     ! enthalpy(3) = dve_pit
  elseif(kif.eq.1) then
     enthalpy(1) = u + p/rho
     ! these derivatives from thermodynamic consistency arguments.
     enthalpy(2) = (1._fp_kind-rtp)*p/rho
     enthalpy(3) = cp*t
     !! temporary special test of free derivatives.
     ! if(ifnoaux_iteration) then
     !   ! in this case eos_jacobian not called so must combine results
     !   ! for u and s to calculate free energy.
     !   enthalpy(1) = u -t*s
     ! else
     !   ! enthalpy(1) = free_fp
     !   enthalpy(1) = free
     ! endif
     ! enthalpy(2) = (p/rho)*density(2)
     ! enthalpy(3) = -s*t + (p/rho)*density(3)
  elseif(kif.eq.2) then
     enthalpy(1) = u + p/rho
     ! these derivatives from thermodynamic consistency arguments.
     enthalpy(2) = (1._fp_kind-rtp)*p/rho/rpt
     enthalpy(3) = cp*t - (rtp-1._fp_kind)*rtp*p/rpt/rho
     !! temporary special test of free derivatives.
     ! if(ifnoaux_iteration) then
     !   ! in this case eos_jacobian not called so must combine results
     !   ! for u and s to calculate free energy.
     !   enthalpy(1) = u -t*s
     ! else
     !   !enthalpy(1) = free_fp
     !   enthalpy(1) = free
     ! endif
     ! enthalpy(2) = p/(rho)
     ! enthalpy(3) = -s*t
  endif
  if(cp.lt.0._fp_kind.and.ifprintcp.eq.1) then
     ifprintcp = 0
     if(verbosity.ge.2) then
        write(stderr,*) 'free_eos_detailed WARNING: negative C_p occurred at least once.'
        write(stderr,*) 'watch out for convection implications.'
     endif
  endif
  if(grada.lt.0._fp_kind.and.ifprintgrada.eq.1) then
     ifprintgrada = 0
     if(verbosity.ge.2) then
        write(stderr,*) 'free_eos_detailed WARNING: negative grada occurred at least once.'
        write(stderr,*) 'watch out for convection implications.'
     endif
  endif
  if(rtp.lt.0._fp_kind) then
     if(ifprintrtp.eq.1) then
        ifprintrtp = 0
        if(verbosity.ge.2) then
           write(stderr,*) 'free_eos_detailed WARNING: negative rtp occurred at least once and for all such cases cf is set to zero'
           write(stderr,*) 'watch out for convection implications.'
        endif
     endif
     cf = 0._fp_kind
  else
     cf = cp*sqrt(rtp)
  endif
  gamma2 = 1._fp_kind/(1._fp_kind-grada)
  gamma3 = 1._fp_kind+gamma1*grada
  xmu1 = ni+ne
  xmu3 = ne
  i = 2
  if(eps(1).gt.0._fp_kind) i = 1
  vx = dne(i)*ni/ne-dni(i)
  if(eps(1).gt.0._fp_kind) then
     fh2 = h_ion/eps(1)      !fraction of H+
     h2rat = 2._fp_kind*h2/eps(1)    !fraction of H2
     h2plusrat = 2._fp_kind*h2plus/eps(1)  !fraction of H2+
  else
     fh2 = 0._fp_kind
     h2rat = 0._fp_kind
     h2plusrat = 0._fp_kind
  endif
  if(eps(2).gt.0._fp_kind) then
     fhe2 = he_ion/eps(2)    !fraction of He+
     fhe3 = he_ion2/eps(2)    !fraction of He++
  else
     fhe2 = 0._fp_kind
     fhe3 = 0._fp_kind
  endif

  ! Calculate square of sound speed according to CG, p.292, and taking
  ! the total energy per unit volume as
  !
  ! u = rho (energy(1) + c^2),
  !

  ! where energy(1) is the internal energy per unit mass
  ! exclusive of rest mass energy,
  !
  ! the relativistically correct square of the sound speed is
  ! clight**2*gamma1/(1.d0 + rho(energy(1) + clight**2)/p) =
  ! (gamma1*p/rho)/((1.d0 + (energy(1) + p/rho)/clight**2)
  !
  ! n.b. one uncertainty is what zero to use for energy.
  ! I have adopted the ground state of H2 as the H energy
  ! zero, and the ground electronic state of non-H
  ! neutrals for the rest of the elements, but should
  ! check with Landau and Lifshitz (referenced in CG) that
  ! this interpretation is correct.  But that check is a
  ! low priority because *if* a change in energy(1)
  ! zero-point is required it will only have a minute
  ! effect on the above result because energy(1)/clight**2 << 1
  sound2 = (gamma1*p/rho)/(1._fp_kind + (energy(1) + p/rho)/(clight*clight))
  if(sound2.le.0._fp_kind.and.ifprintsound2.eq.1) then
     ifprintsound2 = 0
     if(verbosity.ge.2)&
          write(stderr,*) 'free_eos_detailed WARNING: non-positive sound2 occurred at least once.'
  endif

end subroutine free_eos_detailed
