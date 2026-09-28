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

!> This principal module for the free_eos library provides the public
!> generic free_eos procedure that is normally called by those
!> routines that use that library.  That generic procedure provides
!> interfaces to the private free_eos_legacy API (identical to the
!> FreeEOS-2.2.1 API) and the free_eos_modern API (a more powerful API
!> that is available for FreeEOS-3.0.0 and beyond that we recommend
!> for all modern use of the free_eos library).

module mod_free_eos
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public free_eos, version

  !> Interface for the generic subroutine free_eos
  interface free_eos
     module procedure free_eos_legacy
     module procedure free_eos_modern
  end interface free_eos

contains

  !> \return version string
  function version()
    character(len=:), allocatable :: version
    version = '3.0.0'
  end function version

  !> The legacy variant of the free_eos generic API that has been
  !> available for many different releases of FreeEOS.<br>
  !><br>
  !> Note that fp_kind refers to the type of floating-point arguments,
  !> and by FreeEOS build-system default this corresponds to
  !> double-precision floating point arguments, i.e., a real type that
  !> has a precision of at least 15 decimal places and a dynamic range
  !> of at least 10^{-300} to 10^{300}. The build system allows the
  !> user to specify real types (such as extended precision and
  !> quadruple precision) with higher precision or dynamic range, but
  !> those options can be quite slow so they are normally not used
  !> except to investigate (by comparison of results with double
  !> precision results) if there are severe significance loss or
  !> dynamic range issues for the default double precision case.<br>
  !><br>
  !> Note that the combination of the ifoption, ifmodified, and ifion arguments
  !> determines the free_eos_detailed arguments.  (**See notes about
  !> the free_eos option suites and the corresponding calculated
  !> free_eos_detailed arguments below**.)<br>

  !> \param[in] ifoption
  !>   ifoption = 1 (pteh style)<br>
  !>   ifoption = 2 (fjs style)<br>
  !>   ifoption = 3 (mdh style)<br>
  !>   ifoption = 4 (variant of mdh style to fit SCVH)<br>
  !>   ifoption = 5 (no explicit pteh, fjs, or mdh pressure ionization, just use Planck-Larkin)
  !> \param[in] ifmodified
  !>   ifmodified = 0 (original formulation suitable for stellar interiors calculation in particular style)<br>
  !>   ifmodified = 1 (modified version suitable for stellar interiors calculation<br>
  !>   ifmodified = 2 (modify style for best fit of opal tables + extension)<br>
  !>   ifmodified = many other values for a number of different option
  !>     suites. (**See notes about the free_eos option suites and the
  !>     corresponding calculated free_eos_detailed arguments
  !>     below**.)<br>
  !> \param[in] ifion is a flag that controls the way that ionization is
  !>   done.  In general, the lower ifion, the slower the code, and
  !>   the more ionization details that are calculated<br>
  !>   ifion = -2 sets lw=4 and ifmtrace = 0, which implies all 295
  !>     ionization stages of the 20 elements are treated in
  !>     detail<br>
  !>   ifion = -1 sets lw=3 and ifmtrace = 0, which implies minor
  !>     metals are treated as fully ionized while H, He, and the
  !>     major metals are treated in detail<br>
  !>   ifion = 0 (recommended) sets lw = 3, and ifmtrace = 0 below T =
  !>     1.d6 and ifmtrace = 1 (both major and minor metals treated as
  !>     fully ionized) above T = 1.d6<br>
  !>   ifion = 1 sets lw = 3, and ifmtrace = 1, which implies all
  !>     major and minor metals are always treated as fully
  !>     ionized<br>
  !>   ifion = 2 sets lw = 2, and ifmtrace = 1, which implies all
  !>     elements are treated as fully ionized.<br>
  !>   n.b. major metals controlled by array iftracemetal in
  !>     free_eos_detailed.  currently list includes C, N, O, Ne, Mg,
  !>     Si, S, Fe<br>
  !> \param[in] kif_in = 0, ln f = match_variable and tl are independent variables.<br>
  !>   kif_in = 1, ln p = match_variable and tl are independent variables.<br>
  !>   kif_in = -1 sets kif = 1, and ifrad = 2. (**See notes about the
  !>     free_eos option suites and the corresponding calculated
  !>     free_eos_detailed arguments below**.)<br>
  !>   kif_in = 2, ln rho = match_variable and tl are independent
  !>     variables.<br>
  !>   n.b. match_variable and tl are *always* the independent
  !>     variables, and it is the calling routines responsibility to
  !>     place the correct value consistent with kif_in in
  !>     match_variable for each call.  the fl, pl, and rl values in the
  !>     argument list are used *strictly for output*<br>
  !>   n.b. if kif_in is not equal to 0, then much more computer time
  !>     is required to compute the EOS because an fl iteration must
  !>     be used to to match the match_variable.  However, the initial
  !>     guess for fl is improved by a Taylor series approach in this
  !>     case to reduce these fl iterations to a minimum.<br>
  !> \param[in] eps is an abundance array corresponding to the most abundant
  !>   20 elements in the following order:<br>
  !>   H, He, C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl, A, Ca, Ti, Cr,
  !>   Mn, Fe, and Ni<br>
  !>   The ith element of this array is X_i/A_i where X_i is the
  !>   abundance by weight of the ith element (where X_i is normalized
  !>   so that the sums of those abundances are unity) and A_i is the
  !>   atomic weight of the ith element.<br>
  !>   We use the atomic weight scale where (un-ionized) C(12) has a
  !>   weight of 12.00000000....  All weights are for the un-ionized
  !>   element. The eps value for an element should be the sum of the
  !>   individual isotopic eps values for that element.<br>
  !> \param[in] neps_check = 24 = size of eps<br>
  !> \param[in] match_variable is matched by iterative adjustment of fl when kif_in /= 0.
  !>   for kif_in = 0, match_variable is interpreted as fl
  !> \param[in] tl = ln t<br>
  !> \param[out] fl = EFF degeneracy parameter.  For kif_in = 0 this is merely determined
  !>   from fl = match_variable.  For kif_in /= 0 this is determined
  !>   from a Taylor series approach as an initial guess, then
  !>   refined through iteration inside free_eos_detailed.<br>
  !> \param[out] t = temperature<br>
  !> \param[out] rho = density<br>
  !> \param[out] rl = ln rho<br>
  !> \param[out] p = pressure (total if ifrad = 1, prad subtracted if ifrad
  !>   = 0.  (**See notes about the free_eos option suites and the
  !>   corresponding calculated free_eos_detailed arguments
  !>   below**.)<br>
  !> \param[out] pl = ln p<br>
  !> \param[out] cf = cp (- partial ln rho(T,P)/partial T)^{1/2}<br>
  !> \param[out] cp = specific heat at constant pressure<br>
  !> \param[out] qp
  !> \param[out] qf
  !>   qp and qf are legacy parameters which *were* used to calculate an
  !>   approximate gravitational energy generation rate when
  !>   abundances were changing.  These parameters have now been
  !>   disabled (set to 1.d300) since the exact gravitational energy
  !>   generation rate for the case of changing abundances can be
  !>   calculated using the precepts in Strittmatter, P.A.,
  !>   Faulkner, J.. Robertson, J.W., and Faulkner, D.J. 1970, "A
  !>   Question of Entropy", Ap.J. 161, 369-373.
  !>   <https://ui.adsabs.harvard.edu/link_gateway/1970ApJ...161..369S/ADS_PDF>
  !>   From eq. 10a of that paper the preferred rate form is eps_grav =
  !>   -(eps_nuclear - d L/dm - eps_neutrino) = - du/dt - P d (1/rho)/dt,
  !>   where p, rho, and energy(1) (=u) are variables already returned by
  !>   free_eos_modern.  That paper also gives an alternative form (eq. 10b) of
  !>   the rate equation which is equivalent to Cox and Giuli eq.
  !>   17.75''', but that form is deprecated since it will be invalid for
  !>   discontinuous abundances and is much more complicated than eq.
  !>   10a.<br>
  !>  \param[out] sf = partial entropy(fl, tl)/partial fl<br>
  !>  \param[out] st = partial entropy(fl, tl)/partial tl<br>
  !>  \param[out] grada = the adiabatic temperature gradient<br>
  !>  \param[out] rtp = - partial ln rho(T,P) wrt ln T.  n.b. note negative sign so
  !>     quantity should be positive for normal EOS.<br>
  !>  \param[out] rmue = rho/mu_e, where mu_e is the mean molecular
  !>    weight per electron.<br>
  !>  \param[out] fh2 = n(H+)/n(all H)<br>
  !>  \param[out] fhe2 = n(He+)/n(all He)<br>
  !>  \param[out] fhe3 = n(He++)/n(all He)<br>
  !>  \param[out] xmu1 = 1/mu, where mu is the mean molecular weight per particle.<br>
  !>  \param[out] xmu3 = 1/mu_e, where mu is the mean molecular weight per electron.<br>
  !>  \param[out] eta = degeneracy parameter (Cox and Guili eta), related
  !>    to EFF ln f = fl by<br> wf = sqrt(1.d0 + f) = d eta/d fl<br>
  !>  eta = fl+2.d0*(wf-log(1.d0+wf))<br>
  !> \param[out] gamma1 is Gamma_1 as defined by Cox and Guili, eq. 9.88<br>
  !> \param[out] gamma2 is Gamma_2 as defined by Cox and Guili, eq. 9.88<br>
  !> \param[out] gamma3 is Gamma_3 as defined by Cox and Guili, eq. 9.88<br>
  !> \param[out] h2rat = 2*n(H2)/n(all H)<br>
  !> \param[out] h2plusrat = 2*n(H2+)/n(all H)<br>
  !> \param[out] lambda lambda is the Coulomb interaction parameter.  Note
  !>   lambda has two different definitions depending on ifcoulomb.<br>
  !> \param[out] gamma_e gamma_e is the Coulomb diffraction parameter.<br>
  !> \param[out] degeneracy
  !> \param[out] pressure
  !> \param[out] density
  !> \param[out] energy
  !> \param[out] enthalpy
  !> \param[out] entropy <br>
  !>   degeneracy, pressure, density, energy, enthalpy, and entropy
  !>   are all vectors of size 3, with the second component being the derivative
  !>   of the first component wrt match_variable (except for the case where
  !>   kif_in = -1), and the third component being the derivative of
  !>   the first component wrt tl.<br>
  !>   definitions:<br>
  !>     degeneracy(1) is EFF degeneracy parameter ln f defined above.<br>
  !>     pressure(1) = ln pressure (except for kif_in = 2, ifrad=2.
  !>     (**See notes about the free_eos option suites and the
  !>     corresponding calculated free_eos_detailed arguments
  !>     below**.)<br>
  !>     density(1) = ln density.<br>
  !>     energy(1) = internal energy per unit mass.<br>
  !>     enthalpy(1) = enthalpy per unit mass = energy(1) + p/rho.<br>
  !>     entropy(1) = entropy per unit mass.<br>
  !> \param[out] iteration_count
  !>   if iteration_count > 0 then this calculated value is the total
  !>   number of ionization fraction loops completed to obtain the
  !>   solution of the EOS.<br>
  !>   if iteration_count < 0 then this calculated value is the
  !>   negative of the info status variable returned by
  !>   free_eos_detailed (only if that info > 0, i.e., it
  !>   signals an abormal end).<br>
  subroutine free_eos_legacy(ifoption, ifmodified,&
       ifion, kif_in, eps, neps_check, match_variable, tl, fl,&
       t, rho, rl, p, pl, cf, cp, qf, qp, sf, st, grada, rtp,&
       rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
       gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,&
       degeneracy, pressure, density, energy, enthalpy, entropy,&
       iteration_count)
    use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

    ! Arguments
    integer, intent(in) :: ifoption, ifmodified, ifion, kif_in, neps_check
    integer, intent(out) :: iteration_count

    real(fp_kind), intent(in) :: eps(:), match_variable, tl

    real(fp_kind), intent(out) :: fl,&
         t, rho, rl, p, pl, cf, cp, qf, qp, sf, st, grada, rtp,&
         rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
         gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e,&
         degeneracy(:), pressure(:), density(:), energy(:), enthalpy(:), entropy(:)

    ! Local variables

    ! Minimum required size of intent(in) eps array
    integer, parameter :: neps = 24

    ! These variables used to locally store free_eos_modern results that are not
    ! returned by free_eos_legacy to the calling routine.

    ! Use default verbosity level
    integer, parameter :: verbosity = 3

    integer info

    real(fp_kind) s, qe, qv, sound2

    if(neps.gt.neps_check) error stop 'free_eos_legacy: neps_check argument must be at least 20'
    if(neps.gt.size(eps)) error stop 'free_eos_legacy: actual size of eps array must be at least 20'

    call free_eos_modern(verbosity, ifoption, ifmodified,&
         ifion, kif_in, eps, match_variable, tl, fl,&
         t, rho, rl, p, pl, cf, cp, s, sf, st, grada, rtp,&
         qe, qv, rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
         gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2,&
         iteration_count, info,&
         degeneracy, pressure, density, energy, enthalpy, entropy)

    ! Disable qp and qf just like in the previous non-generic free_eos
    qp = 1.e300_fp_kind
    qf = 1.e300_fp_kind

    ! The logic for free_eos_modern and all routines it calls implies
    ! info == 0 on success, and info > 0 on failure (with the value
    ! describing what the error condition was).
    ! N.B. info == 0 (success) simply drops through the following Boolean logic.
    if(info.lt.0) then
       ! Sanity check
       write(stderr,*) "info returned by free_eos_modern = ", info
       flush(stderr)
       error stop "free_eos_legacy: negative info returned by free_eos_modern is not allowed"
    elseif(info.gt.0) then
       ! Overwrite iteration_count with a negative value on failure which is the
       ! lame way that free_eos_legacy is forced to communicate positive info
       ! because of the limitations of the free_eos_legacy API.
       iteration_count = -info
    endif

  end subroutine free_eos_legacy

  !> The recommended modern variant of the free_eos generic API that
  !> has been available since the release of FreeEOS-3.0.0.<br>
  !><br>
  !> Note that fp_kind refers to the type of floating-point arguments,
  !> and by FreeEOS build-system default this corresponds to
  !> double-precision floating point arguments, i.e., a real type that
  !> has a precision of at least 15 decimal places and a dynamic range
  !> of at least 10^{-300} to 10^{300}. The build system allows the
  !> user to specify real types (such as extended precision and
  !> quadruple precision) with higher precision or dynamic range, but
  !> those options can be quite slow so they are normally not used
  !> except to investigate (by comparison of results with double
  !> precision results) if there are severe significance loss or
  !> dynamic range issues for the default double precision case.<br>
  !><br>
  !> Note that the combination of the ifoption, ifmodified, and ifion arguments
  !> determines the free_eos_detailed arguments.  (**See notes about
  !> the free_eos option suites and the corresponding calculated
  !> free_eos_detailed arguments below**.)<br>

  !> \param[in] verbosity
  !>   verbosity = 0 means that no output from free_eos is allowed
  !>   other than explanatory messages prior to "error stop"
  !>   statements.  Those statements are only associated with internal
  !>   logic errors or a bad call of free_eos (e.g., with incorrect
  !>   array sizes) which rarely happens so verbosity = 0 typically
  !>   means no output.  (In which case, it is recommended the calling
  !>   routine should carefully parse any non-zero info value returned
  !>   to discover the source of errors.)<br>
  !>   verbosity >= 1 means that error messages corresponding
  !>   to non-zero info results are allowed.<br>
  !>   verbosity >= 2 means that warning messages are allowed.<br>
  !>   verbosity >= 3 means that important (e.g., BFGS) debugging messages are allowed.<br>
  !> \param[in] ifoption
  !>   ifoption = 1 (pteh style)<br>
  !>   ifoption = 2 (fjs style)<br>
  !>   ifoption = 3 (mdh style)<br>
  !>   ifoption = 4 (variant of mdh style to fit SCVH)<br>
  !>   ifoption = 5 (no explicit pteh, fjs, or mdh pressure ionization, just use Planck-Larkin)
  !> \param[in] ifmodified
  !>   ifmodified = 0 (original formulation suitable for stellar interiors calculation in particular style)<br>
  !>   ifmodified = 1 (modified version suitable for stellar interiors calculation<br>
  !>   ifmodified = 2 (modify style for best fit of opal tables + extension)<br>
  !>   ifmodified = many other values for a number of different test cases for each ifoption case<br>
  !> \param[in] ifion is a flag that controls the way that ionization is done.  In
  !>   general, the lower ifion, the slower the code, and the more
  !>   ionization details that are calculated<br>
  !>   ifion = -2 sets lw=4 and ifmtrace = 0, which implies all 295 ionization
  !>     stages of the 20 elements are treated in detail<br>
  !>   ifion = -1 sets lw=3 and ifmtrace = 0, which implies minor metals
  !>     are treated as fully ionized while H, He, and the major
  !>     metals are treated in detail<br>
  !>   ifion = 0 (recommended) sets lw = 3, and ifmtrace = 0 below T = 1.d6
  !>     and ifmtrace = 1 (both major and minor metals treated as
  !>     fully ionized) above T = 1.d6<br>
  !>   ifion = 1 sets lw = 3, and ifmtrace = 1, which implies all major
  !>     and minor metals are always treated as fully ionized<br>
  !>   ifion = 2 sets lw = 2, and ifmtrace = 1, which implies all elements
  !>     are treated as fully ionized.<br>
  !>   n.b. major metals controlled by array iftracemetal in free_eos_detailed.
  !>     currently list includes C, N, O, Ne, Mg, Si, S, Fe<br>
  !> \param[in] kif_in = 0, ln f = match_variable and tl are independent variables.<br>
  !>   kif_in = 1, ln p = match_variable and tl are independent variables.<br>
  !>   kif_in = -1 sets kif = 1, and ifrad = 2.  (**See notes about
  !>   the free_eos option suites and the corresponding calculated
  !>   free_eos_detailed arguments below**.)<br>
  !>   kif_in = 2, ln rho = match_variable and tl are independent variables.<br>
  !>   n.b. match_variable and tl are *always* the independent variables,
  !>   and it is the calling routines responsibility to place the correct
  !>   value consistent with kif_in in match_variable for each call.
  !>   the fl, pl, and rl values in the argument list are used *strictly for
  !>   output*<br>
  !>   n.b. if kif_in is not equal to 0, then much more computer time is
  !>   required to compute the EOS because an fl iteration must
  !>   be used to to match the match_variable.  However, the initial guess
  !>   for fl is improved by a Taylor series approach in this case to
  !>   reduce these fl iterations to a minimum.<br>
  !> \param[in] eps is an abundance array corresponding to the most abundant
  !>   20 elements in the following order:<br>
  !>   H, He, C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl, A, Ca, Ti, Cr,
  !>   Mn, Fe, and Ni<br>
  !>   The ith element of this array is X_i/A_i where X_i is the
  !>     abundance by weight of the ith element (where X_i is
  !>     normalized so that the sums of those abundances are unity)
  !>     and A_i is the atomic weight of the ith element.<br>
  !>   We use the atomic weight scale where (un-ionized) C(12)
  !>     has a weight of 12.00000000....  All weights are for the
  !>     un-ionized element. The eps value for an element
  !>     should be the sum of the individual isotopic eps values for
  !>     that element.<br>
  !> \param[in] match_variable is matched by iterative adjustment of fl when kif_in /= 0.
  !>   for kif_in = 0, match_variable is interpreted as fl
  !> \param[in] tl = ln t<br>
  !> \param[out] fl = EFF degeneracy parameter.  For kif_in = 0 this is merely determined
  !>   from fl = match_variable.  For kif_in /= 0 this is determined
  !>   from a Taylor series approach as an initial guess, then
  !>   refined through iteration inside free_eos_detailed.<br>
  !> \param[out] t = temperature<br>
  !> \param[out] rho = density<br>
  !> \param[out] rl = ln rho<br>
  !> \param[out] p = pressure (total if ifrad = 1, prad subtracted if ifrad
  !>   = 0.  (**See notes about the free_eos option suites and the
  !>   corresponding calculated free_eos_detailed arguments
  !>   below**.)<br>
  !> \param[out] pl = ln p<br>
  !> \param[out] cf = cp (- partial ln rho(T,P)/partial T)^{1/2}<br>
  !> \param[out] cp = specific heat at constant pressure<br>
  !> \param[out] s is the entropy per unit mass.<br>
  !> \param[out] sf = partial s(fl, tl)/partial fl<br>
  !> \param[out] st = partial s(fl, tl)/partial tl<br>
  !> \param[out] grada = the adiabatic temperature gradient<br>
  !> \param[out] rtp = - partial ln rho(T,P) wrt ln T.  n.b. note negative sign so
  !>  quantity should be positive for normal EOS.<br>
  !> \param[out] qe = energy per unit mass<br>
  !> \param[out] qv = 1/rho<br>
  !>   qe and qv allow the calling routine to calculate the
  !>   gravitational energy generation rate = -(d qe/dt + p*d qv/dt)
  !>   in a convenient way from first principles (the first law of
  !>   thermodynamics) rather than using deprecated approximations
  !>   (e.g., -T ds/dt) for this rate.  See Strittmatter, P.A.,
  !>   Faulkner, J.. Robertson, J.W., and Faulkner, D.J. 1970, "A
  !>   Question of Entropy", Ap.J. 161, 369-373.
  !>   <https://ui.adsabs.harvard.edu/link_gateway/1970ApJ...161..369S/ADS_PDF>
  !> \param[out] rmue = rho/mu_e, where mu_e is the mean molecular weight
  !>   per electron.<br>
  !> \param[out] fh2 = n(H+)/n(all H)<br>
  !> \param[out] fhe2 = n(He+)/n(all He)<br>
  !> \param[out] fhe3 = n(He++)/n(all He)<br>
  !> \param[out] xmu1 = 1/mu, where mu is the mean molecular weight per particle.<br>
  !> \param[out] xmu3 = 1/mu_e, where mu is the mean molecular weight per electron.<br>
  !> \param[out] eta = degeneracy parameter (Cox and Guili eta), related
  !>    to EFF ln f = fl by<br>
  !>    wf = sqrt(1.d0 + f) = d eta/d fl<br>
  !>    eta = fl+2.d0*(wf-log(1.d0+wf))<br>
  !> \param[out] gamma1 is Gamma_1 as defined by Cox and Guili, eq. 9.88<br>
  !> \param[out] gamma2 is Gamma_2 as defined by Cox and Guili, eq. 9.88<br>
  !> \param[out] gamma3 is Gamma_3 as defined by Cox and Guili, eq. 9.88<br>
  !> \param[out] h2rat = 2*n(H2)/n(all H)<br>
  !> \param[out] h2plusrat = 2*n(H2+)/n(all H)<br>
  !> \param[out] lambda = Coulomb interaction parameter (two definitions depending on ifcoulomb).<br>
  !> \param[out] gamma_e = Coulomb diffraction parameter.<br>
  !> \param[out] sound2 = the square of the relativistically corrected
  !>   sound speed, see CG eq. 14.29 as corrected in a footnote.<br>
  !> \param[out] iteration_count =
  !>   the total number of ionization fraction iterations completed to
  !>   obtain the solution<br>
  !> \param[out] info returns 0 on success and a non-zero value on failure.
  !>   Each unique non-zero value indicates the routine where the
  !>   problem occurred (see the info_offset data in
  !>   src/mod_info_data.f90) and once that offset is subtracted for a
  !>   given routine, the index of the problem that occurred which is
  !>   indicated in each routine in the src directory that returns an
  !>   info value.<br>
  !> \param[out] degeneracy
  !> \param[out] pressure
  !> \param[out] density
  !> \param[out] energy
  !> \param[out] enthalpy
  !> \param[out] entropy <br>
  !>   degeneracy, pressure, density, energy, enthalpy, and entropy
  !>   are all vectors of size 3, with the second component being the derivative
  !>   of the first component wrt match_variable (except for the case where
  !>   kif_in = -1), and the third component being the derivative of
  !>   the first component wrt tl.<br>
  !> definitions:<br>
  !>   degeneracy(1) is EFF degeneracy parameter ln f defined above.<br>
  !>   pressure(1) = ln pressure (except for kif_in = 2, ifrad=2.
  !>     (**See notes about the free_eos option suites and the
  !>     corresponding calculated free_eos_detailed arguments
  !>     below**.)<br>
  !>   density(1) = ln density.<br>
  !>   energy(1) = internal energy per unit mass.<br>
  !>   enthalpy(1) = enthalpy per unit mass = energy(1) + p/rho.<br>
  !>   entropy(1) = entropy per unit mass.<br>
  !>
  !> **Notes about the free_eos option suites and the corresponding
  !>   calculated free_eos_detailed arguments**<br>
  !> The purpose of the ifoption, ifmodified, and ifion arguments of
  !> the free_eos_legacy and free_eos_modern variants of the free_eos
  !> generic routine is to specify the option suite via the calculated
  !> free_eos_detailed arguments. <br>
  !> The free_eos_detailed routine implements many different option suites, but the
  !> noteable ones are EOS1, EOS1a, EOS2, EOS3, and EOS4
  !> which are ordered by increasing speed and decreasing accuracy of
  !> the corresponding free-energy model.<br>
  !> These option suites correspond to the following
  !> values of ifoption, ifmodified, and ifion:<br>
  !> EOS1: ifoption, ifmodified, ifion = 3, 1, -2<br>
  !> EOS1 is our recommended EOS that has been constrained by OPAL and
  !>   SCVH fits.  The ifion=-2 flag means that all 295 ionization
  !>   states of the 20 most abundant elements are treated in
  !>   detail.<br>
  !> EOS1a: ifoption, ifmodified, ifion = 3, 1, -1<br>
  !> EOS1a is identical to EOS1 except that minor metals are
  !>   approximated as fully ionized (because of the different ifion
  !>   flag, see above) which increases the speed of the computation
  !>   by almost a factor of three at the expense of almost negligible
  !>   pressure errors at low temperatures.<br>
  !> EOS2: ifoption, ifmodified, ifion = 2, 1, -1<br>
  !> EOS2 is a fit of EOS1 using a modified form of the Fritz Swenson
  !>   pressure ionization approximation.<br>
  !> EOS3: ifoption, ifmodified, ifion = 1, 1, 0<br>
  !> EOS3 is a fit of EOS1 using a modified PTEH pressure ionization
  !>   and using the PTEH approximation for the Coulomb sums above log
  !>   T = 6.  The ifion flag means the major metals are assumed fully
  !>   ionized above log T = 6.<br>
  !> EOS4: ifoption, ifmodified, ifion = 1, 101, 0<br>
  !> EOS4 is the same as EOS3 except for using the PTEH approximation
  !>   for the Coulomb sums for **all** temperatures. EOS4 produces
  !>   (if a solar calibration is used) good results even for the
  !>   extreme LMS if **just** reliable radii and luminosities are
  !>   required from the interior model.  For more detailed results
  !>   such as vibrational frequencies and the highest quality radii
  !>   and luminosities, EOS1 is recommended instead.<br>
  !>
  !> **Terse notes about calculated free_eos_detailed arguments**.<br>
  !>   See the free_eos_modern source code about the details
  !>   concerning how various combinations of these arguments
  !>   correspond to free_eos_modern (and free_eos_legacy) ifoption,
  !>   ifmodified, ifion, and kif_in arguments.<br>
  !> For additional details about these arguments consult the source
  !>   code for free_eos_detailed and the routines that it calls.<br>
  !> kif = kif_in if kif_in >= 0<br>
  !> kif = 1 and ifrad = 2 if kif_in = -1<br>
  !> ifh2 = 0  no h2<br>
  !> ifh2 = 1  vdb h2<br>
  !> ifh2 = 2  st h2<br>
  !> ifh2 = 3  irwin h2<br>
  !> ifh2plus = 0  no h2plus<br>
  !> ifh2plus = 1  st h2plus<br>
  !> ifh2plus = 2  irwin h2plus<br>
  !> morder = 3, 5, or 8  (or 13, 15, or 18) means 3rd, 5th, or 8th order<br>
  !> fermi_dirac integral approximation following EFF fit (or modified
  !> version of EFF fit which reduces to Cody-Thacher approximation for
  !> low relativistic correction).<br>
  !> morder = -3, -5, or -8 (or -13, -15, or -18) means use above
  !> approximations in non-relativistic limit.<br>
  !> morder = 1 uses Cody-Thacher approximation directly.<br>
  !> morder = 21 calculate Fermi-Dirac integrals with slow, but precise
  !> (~1.d-9 relative errors) numerical integration.<br>
  !> morder = -21 is same as 21 in non-relativistic limit.<br>
  !> morder = 23 means original 3rd order eff result<br>
  !> morder = -23 is same as 23 in non-relativistic limit.<br>
  !> if abs(ifexchange_in) > 100 then use linear transform approximation
  !>   for exchange treatment.  Otherwise, use numerical transform as described in Paper IV.<br>
  !>   ifexchange = mod(ifexchange_in,100):<br>
  !>   One other wrinkle on deciding ifexchange is that large degeneracy
  !>   approximations (mod(ifexchange,10) = 2) only allowed above psi_lim.<br>
  !> ifexchange:<br>
  !> To understand ifexchange, must summarize free-energy model of fex used
  !> in research note.  In general from kapusta relation,<br>
  !> fex is proportional to I - J + 2^1.5 pi^2/3 beta^2 K, where<br>
  !> K and J and known integrals which are functions of psi and beta,<br>
  !> and I = K^2. (see Kovetz et al 1972, ApJ 174, 109, hereafter KLVH)<br>
  !>  0  < ifexchange < 10 --> KLVH treatment (i.e., drop K term).<br>
  !> 10 <= ifexchange < 20 --> Kapusta treatment (i.e., retain K term).<br>
  !> 20 <= ifexchange < 30 --> special test case, use I term alone.<br>
  !> 30 <= ifexchange < 40 --> special test case, use J term alone.<br>
  !> 40 <= ifexchange < 50 --> special test case, use K term alone.<br>
  !> N.B. if ifexchange is negative, then use non-relativistic limit
  !>   of corresponding positive ifexchange option.<br>
  !> ifexchange details:<br>
  !> ifexchange = 1 --> G(psi) + weak relativistic correction from KLVH<br>
  !> ifexchange = 2 --> I-J degenerate expression from KLVH (with corrected
  !>   sign error on a2 and high numerical precision a2 and a3, see
  !>   paper IV)<br>
  !> ifexchange = 11 or 12 K term added (both series from CG).<br>
  !> ifexchange = 21 or 22 I term alone (both series from KLVH)<br>
  !> ifexchange = 31 or 32 -J term alone (both series from KLVH).<br>
  !> ifexchange = 41 or 42 K term alone.<br>
  !> mod(ifexchange,10) = 4 is lowest order fit of J, K<br>
  !> mod(ifexchange,10) = 5 is next higher order fit of J, K<br>
  !> mod(ifexchange,10) = 6 is highest order fit of J, K<br>
  !> ifcoulomb > 9, do diffraction correction<br>
  !> ifcoulomb > 0 means do metal Coulomb contribution to sum0 and sum2 exactly<br>
  !> -10 < ifcoulomb < 0 means do metal Coulomb contribution to sum0 and
  !>   sum2 using the "metal Coulomb" approximation.<br>
  !>   n.b. this approximation is internally replaced by PTEH
  !>   approximation (next line) if either eps(1) or eps(2) are zero.<br>
  !> ifcoulomb < -9 means use PTEH approximation for sum0 and sum2.<br>
  !>   n.b. if this approximation is combined with the PTEH
  !>     approximation for pressure ionization, then the
  !>     routine is much faster because no ionization fraction
  !>     iterations are required.<br>
  !> mod(|ifcoulomb|,10) = 0 ignore Coulomb interaction.<br>
  !> mod(|ifcoulomb|,10) = 1 use Debye-Huckel Coulomb approximation.<br>
  !> mod(|ifcoulomb|,10) = 2 use Debye-Huckel Coulomb approximation with tau(x) correction<br>
  !> mod(|ifcoulomb|,10) = 3 use PTEH Coulomb approximation with their theta_e<br>
  !> mod(|ifcoulomb|,10) = 4 use PTEH Coulomb approximation with fermi-dirac theta_e<br>
  !> mod(|ifcoulomb|,10) = 9 same as 4 with DeWitt definition of lambda
  !>   (using sum0a = sum0 + ne*theta_e).<br>
  !> mod(|ifcoulomb|,10) = 5 use DH smoothly connected to modified OCP
  !>   with DeWitt definition of lambda.<br>
  !> mod(|ifcoulomb|,10) = 6 same as 5 with alternative smooth connection<br>
  !> mod(|ifcoulomb|,10) = 7 DH (Gamma < 1) or OCP using new DeWitt lambda.<br>
  !> mod(|ifcoulomb|,10) = 8 same as 7 with theta_e = 0.<br>
  !> ifpi contains meaning of two flags:<br>
  !> ifpi > 0 means use Planck-Larkin occupation probability, otherwise not.<br>
  !> remaining meaning in absolute value of ifpi<br>
  !> |ifpi| = 0, use no pressure ionization<br>
  !> |ifpi| = 1, use pteh pressure ionization<br>
  !> |ifpi| = 2, use fjs pressure ionization<br>
  !> |ifpi| = 3, use MDH pressure ionization<br>
  !> |ifpi| = 4, use Saumon-like variation of MDH pressure ionization<br>
  !> |ifpi| > 4 same as zero, i.e., use no pressure ionization.<br>
  !> ifrad = 0, no radiation pressure,<br>
  !>   input match_variable is consistent (ln P excluding radiation
  !>   pressure) for kif = 1.<br>
  !> ifrad = 1, radiation pressure included,<br>
  !>   input match_variable is consistent (ln P including radiation
  !>   pressure) for kif = 1.<br>
  !> ifrad = 2, radiation pressure included,<br>
  !>   input match_variable is ln(ptotal-prad) for kif = 1, but
  !>   all output quantities are calculated with radiation included.<br>
  !>   this feature is used to reduce significance loss in regions which
  !>   are dominated by radiation pressure.<br>
  !>   for kif = 2, the input match_variable is ln rho as per normal, but
  !>   the output pressure(1) is ln(ptotal-prad).<br>
  !>   n.b. this latter case is only used for some tables which have
  !>   rho and T as the independent variable and all quantities including
  !>   pressure derivatives calculated with radiation pressure *except for*
  !>   the pressure itself.<br>
  !>   Also note for this latter case
  !>   (ifrad = 2, kif = 2) that pressure(1) = ln (ptotal-prad) will be
  !>   inconsistent with pressure(2) and pressure(3) which will be
  !>   partials of ln ptotal wrt ln rho and ln t.<br>
  !> lw < 3 every element treated as fully ionized.<br>
  !> lw = 3 trace metals treated as fully ionized.<br>
  !> lw > 3 very slow option with all elements treated as partially ionized.<br>
  !> iftc = 1 only used for thermodynamic consistency tests on entropy(2)
  !> (in which case entropy(2) returned via rtp).<br>

  subroutine free_eos_modern(verbosity, ifoption, ifmodified,&
       ifion, kif_in, eps, match_variable, tl, fl,&
       t, rho, rl, p, pl, cf, cp, s, sf, st, grada, rtp,&
       qe, qv, rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
       gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2,&
       iteration_count, info,&
       degeneracy, pressure, density, energy, enthalpy, entropy)
    use mod_free_eos_detailed, only: free_eos_detailed
    use mod_flow_data, only: ln_underflow_limit, ln_overflow_limit
    use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

    ! Arguments
    integer, intent(in) :: verbosity, ifoption, ifmodified, ifion, kif_in
    integer, intent(out) :: iteration_count, info

    real(fp_kind), intent(in) :: eps(:), match_variable, tl

    real(fp_kind), intent(out) :: fl,&
         t, rho, rl, p, pl, cf, cp, s, sf, st, grada, rtp,&
         qe, qv, rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
         gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2

    real(fp_kind), intent(out), optional :: &
         degeneracy(:), pressure(:), density(:), energy(:), enthalpy(:), entropy(:)

    ! Local variables

    ! Parameters

    ! This parameter should be changed to .true. if any of the debug
    ! options are being tried in free_eos_detailed.
    logical, parameter :: debug_any = .false.

    ! Minimum required size of intent(in) eps array
    integer, parameter :: neps = 24

    ! Minimum required size of intent(out) arrays
    integer, parameter :: nderivp1 = 3

    ! make sure change in match_variable and tl is not too large
    ! for first-order Taylor series.  This is my guess of
    ! the larger limit allowed, but if you run into trouble
    ! you can always reduce this limit at the expense of starting
    ! with a safe fl initial value which then will require a
    ! large number of iterations to refine.
    ! the following values are slightly larger than the opal
    ! grid spacing: delta log10 rho = .25d0 (==> delta ln rho < 0.58,
    ! delta log10 t ~ 0.1d0 (==> delta ln t < 0.24).
    ! n.b. maximum delta tl, fl larger than warm start criterion
    ! in free_eos_detailed, but cold start in free_eos_detailed
    ! is reliable and more efficient than cold start for "safe" fl
    ! which then must be iterated a lot more times.
    real(fp_kind), parameter :: dm_nrlim = 0.6_fp_kind, dt_nrlim = 0.25_fp_kind

    logical ifcheck_abrupt
    data ifcheck_abrupt/.true./

    integer &
         ifh2, ifh2plus, morder, ifexchange_in,&
         ifmtrace,&
         ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,&
         ifreducedmass, iftc,&
         kif

    integer ncall, ncall_abrupt, ncall_start
    data ncall, ncall_abrupt, ncall_start /2*0,1/

    real(fp_kind) fl_old, fm, ft,&
         match_variable_old, tl_old, dm, dt, tllim,&
         dfl

    ! must be ridiculous values
    data match_variable_old, tl_old, fl_old/3*1.e30_fp_kind/
    ! must be non-zero
    data fm, ft/2*1._fp_kind/

    real(fp_kind), allocatable ::&
         local_degeneracy(:), local_pressure(:), local_density(:),&
         local_energy(:), local_enthalpy(:), local_entropy(:)
    ! Save everything that is initialized with data statements
    ! N.B. From the logic below (and also from my uninitialized tests)
    ! there doesn't appear to be any other variables associated with
    ! these saved variables that need to be saved.
    save ifcheck_abrupt, ncall, ncall_abrupt, ncall_start,&
         match_variable_old, tl_old, fl_old,&
         fm, ft

    ! Sanity checks
    if(neps.gt.size(eps)) error stop 'free_eos_modern: actual size of eps array must be at least 20'
    allocate(&
         local_degeneracy(nderivp1),&
         local_pressure(nderivp1),&
         local_density(nderivp1),&
         local_energy(nderivp1),&
         local_enthalpy(nderivp1),&
         local_entropy(nderivp1))

    ! takes care of kif_in = -1 case.
    kif = iabs(kif_in)
    ! temperature limit for a number of options
    tllim = log(1.e6_fp_kind)
    ! starting (or final when kif=0) value for fl:
    if(kif.eq.0) then
       fl = match_variable
    elseif(kif.eq.1.or.kif.eq.2) then
       dm = match_variable - match_variable_old
       dt = tl - tl_old
       dfl = fm*dm + ft*dt
       ! N.B. for kif = 1, must be much more conservative than for kif = 2
       ! since non-ideal effects tend to be the same fraction of the gas
       ! pressure regardless of density for intermediate and low densities.
       if((abs(dm).le.dm_nrlim.and.abs(dt).le.dt_nrlim).or.&
            (kif.eq.2.and.max(fl_old, fl_old + dfl).lt.-10._fp_kind.and.abs(dfl).le.-20._fp_kind)) then
          ! cannot get into trouble with Taylor-series approach
          ! for small dm and dt or for kif=2 with small fl_old, and fl.
          fl = fl_old + dfl
       elseif(kif.eq.2.and.abs(dfl).le.20._fp_kind) then
          ! For kif = 2 and not a grotesquely large step, then decrementing
          !   by 20 is a safe option.
          fl = max(-100._fp_kind, fl_old - 20._fp_kind)
       else
          ! this is a safe starting value which is inefficient at high
          ! density but quite efficient (especially for kif=1) at low density.
          fl = -100._fp_kind
       endif
       ! must always specify same starting fl if doing derivative tests
       ! in free_eos_detailed.
       if(debug_any) fl = -100._fp_kind
    else
       error stop 'free_eos_modern: bad kif value'
    endif
    if((2.le.ifoption.and.ifoption.le.4).and.ifcheck_abrupt) then
       ! check for abrupt changes in variables and warn if this
       ! occurs too often when auxiliary variable iterations
       ! are required
       ncall = ncall + 1
       if(abs(match_variable-match_variable_old).gt.dm_nrlim.or.abs(tl-tl_old).gt.dt_nrlim) then
          ncall_abrupt = ncall_abrupt + 1
       elseif(ncall_start.eq.1) then
          ! Start real counting after initial cold start with lots of abrupt changes.
          ncall = 0
          ncall_abrupt = 0
          ncall_start = 0
       endif
       ! sample first 100 real calls after initial cold start to see
       ! what fraction have excessively abrupt changes.
       if(ncall.eq.100) then
          if(ncall_abrupt.ge.10) then
             ifcheck_abrupt = .false.
             if(verbosity.ge.2) then
                write(stderr,*) 'free_eos_modern warning: too many cold starts of the eos'
                write(stderr,*) 'smaller step sizes corresponding to'
                write(stderr,'(a,f3.1,a,f3.1)') '|delta match_variable| < ',dm_nrlim, ' and |delta tl| < ',dt_nrlim
                write(stderr,*) 'will make the eos run *much* faster per call'
             endif
          endif
          ncall = 0
          ncall_abrupt = 0
       endif
    endif
    tl_old = tl
    match_variable_old = match_variable
    ! set ifmtrace and lw according to ifion
    if(ifion.eq.-2) then
       ! all ionization stages of all elements treated in detail.
       lw = 4
       ifmtrace = 0
    elseif(ifion.eq.-1) then
       ! minor metals approximated as fully ionized.
       lw = 3
       ifmtrace = 0
    elseif(ifion.eq.0) then
       lw = 3
       if(tl.lt.tllim) then
          ! only minor metals approximated as fully ionized.
          ifmtrace = 0
       else
          ! minor and major metals approximated as fully ionized.
          ifmtrace = 1
       endif
    elseif(ifion.eq.1) then
       ! minor and major metals approximated as fully ionized.
       lw = 3
       ifmtrace = 1
    elseif(ifion.eq.2) then
       ! all elements approximated as fully ionized.
       lw = 2
       ifmtrace = 1
    else
       error stop 'free_eos_modern: bad ifion flag'
    endif
    !  default behaviour unless otherwise specified below
    ifh2 = 3
    ifh2plus = 2
    morder = 13
    if(ifmodified.le.0) then
       ifexchange_in = 0
    else
       ! Non-linear transform, Kapusta term, lowest-order EFF-style
       ! approximation for J and K integrals.
       ifexchange_in = 14
    endif
    ! default treatment of radiation pressure.
    if(kif_in.eq.-1) then
       ifrad = 2
    else
       ifrad = 1
    endif

    if(lw.le.2.or.(ifmodified.le.0)) then
       ifexcited = 0
    else
       ! excitation approximation works well here
       ifexcited = 2
    endif

    ! Use infinite nmax approximation by default.  In some cases
    ! (e.g., no excitation included in free energy model at all)
    ! nmax has no effect on results at all, but it is used
    ! in free_eos_detailed in comparing old and new values of nmax
    ! regardless of free-energy model so to avoid valgrind complaints
    ! driven by that uninitialized if logic set it here in all cases.
    nmax = 300001

    ! usually iftc = 0 (entropy(2) done with straight derivatives.  there
    ! is some evidence that where entropy derivative is done with
    ! thermodynamic relations, significance loss is more severe (especially
    ! in calculation of Q = rtp), but the difference is usually at the 1.d-13
    ! level or better.
    iftc = 0
    ! in all but one case, use reduced mass of electron in equilibrium
    !   constant equation.
    ifreducedmass = 1
    if(ifoption.eq.1) then
       ! pteh style
       ! default ifcouloumb
       if(tl.ge.tllim) then
          ifcoulomb = -15
       else
          ifcoulomb = 5
       endif
       if(ifmodified.eq.0) then
          ! original version for interior calculations.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
       elseif(ifmodified.eq.1) then
          ! ****EOS3****
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.2) then
          ! modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 1
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 4
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.102) then
          ! modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ! same as 2 except for Debye-Huckel
          ifcoulomb = 1
          ifpi = 1
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 4
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.3) then
          ! modified version to fit old opal table extension
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 1
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 4
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.4) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with no exchange
          ifexchange_in = 0
       elseif(ifmodified.eq.5) then
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.6) then
          ! non-linear, Kovetz et al, no strong degeneracy approximation.
          ! G(eta) + relativistic Kapusta II equiv Kovetz terms
          ifexchange_in = 1
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.7) then
          ! non-linear, Kovetz et al, possible strong degeneracy approximation.
          ifexchange_in = 2
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.8) then
          ! non-linear, Kapusta, no strong degeneracy approximation
          ifexchange_in = 11
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.9) then
          ! linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 112
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.11) then
          ! modified version with enhancements for interior calculations.
          ! but *always* used PTEH Coulomb sum approximation
          ifcoulomb = -15
          ifpi = 1
       elseif(ifmodified.eq.12) then
          ! modified version with enhancements for interior calculations.
          ! but use case b (metal Coulomb sum approximation)
          ifcoulomb = -5
          ifpi = 1
       elseif(ifmodified.eq.13) then
          ! modified version with enhancements for interior calculations.
          ! but use no Coulomb sum approximation (recommended for work
          ! of high precision)
          ifcoulomb = 5
          ifpi = 1
       elseif(ifmodified.eq.15) then
          ! Cody-Thacher
          ! modified version with enhancements for interior calculations.
          morder = 1
          ifpi = 1
       elseif(ifmodified.eq.16) then
          ! non-relativistic
          ! modified version with enhancements for interior calculations.
          morder = -3
          ifpi = 1
       elseif(ifmodified.eq.17) then
          ! non-relativistic
          ! modified version with enhancements for interior calculations.
          morder = -5
          ifpi = 1
       elseif(ifmodified.eq.18) then
          ! non-relativistic
          ! modified version with enhancements for interior calculations.
          morder = -8
          ifpi = 1
       elseif(ifmodified.eq.20) then
          ! drop h2 and h2+
          ! modified version with enhancements for interior calculations.
          ifh2 = 0
          ifh2plus = 0
          ifpi = 1
       elseif(ifmodified.eq.21) then
          ! qh2
          ! modified version with enhancements for interior calculations.
          ifh2 = 4
          ifpi = 1
       elseif(ifmodified.eq.22) then
          ! drop h2+
          ! modified version with enhancements for interior calculations.
          ifh2plus = 0
          ifpi = 1
       elseif(ifmodified.eq.24) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with no Coulomb
          ifcoulomb = 0
       elseif(ifmodified.eq.25) then
          ! no excitation or Planck-Larkin
          ! modified version with enhancements for interior calculations.
          ifpi = -1
          ifexcited = 0
       elseif(ifmodified.eq.26) then
          ! no excitation with Planck-Larkin
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ifexcited = 0
       elseif(ifmodified.eq.27) then
          ! Debye-Hueckel
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 1
          ifpi = 1
       elseif(ifmodified.eq.28) then
          ! Debye-Hueckel with tau correction
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 2
          ifpi = 1
       elseif(ifmodified.eq.29) then
          ! PTEH original Coulomb
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 3
          ifpi = 1
       elseif(ifmodified.eq.30) then
          ! thermodynamic consistency
          iftc = 0
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.31) then
          ! thermodynamic consistency
          iftc = 1
          ! modified version with enhancements for interior calculations.
          ifpi = 1
       elseif(ifmodified.eq.40) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with "exact" excitation calculation for H, He, and molecules
          if(lw.gt.2) then
             ifexcited = 12
             nmax = 10
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.41) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with approximate excitation calculation for H, He only
          ! with default "infinite" nmax
          if(lw.gt.2) then
             ifexcited = 1
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.42) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with approximation excitation calculation for H, He
          !   + molecules + metals and with default "infinite" nmax
          if(lw.gt.2) then
             ifexcited = 3
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.43) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with approximate excitation calculation for H, He, and molecules
          ! with special nmax
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 4
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.44) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with approximate excitation calculation for H, He, and molecules
          ! with special nmax
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 10
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.45) then
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with approximate excitation calculation for H, He, and molecules
          ! with special nmax
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 50
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.101) then
          ! ****EOS4****
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! *always* used PTEH Coulomb sum approximation
          ! (causes substantial errors for LMS, but much quicker)
          ifcoulomb = -15
          ! linearly transformed relativistic exchange
          ! (causes substantial errors for LMS, but quicker)
          ! n.b. found didn't matter that much
          ! linear, Kapusta, possible strong degeneracy approximation.
          ! ifexchange_in = 112
          !  use lower-order approximation for fermi-dirac integrals?
          !  The maximum errors are similar to the higher-order approximation,
          !  (1.d-3 in ln P), but the larger errors are more widespread than
          !  the morder = 5 case.  The fermi-dirac overhead doesn't matter
          !  much for other forms of the EOS, but for the quick form it might
          !  be useful to eliminate this extra source of overhead.
          !  n.b. turned out to matter very little
          ! morder = 3
       elseif(ifmodified.eq.103) then
          ! emulation of Sweigart EOS
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with no Coulomb and no exchange
          ifcoulomb = 0
          ifexchange_in = 0
       elseif(ifmodified.eq.104) then
          ! emulation of Sweigart EOS
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with no Coulomb and no exchange
          ifcoulomb = 0
          ifexchange_in = 0
          !      and no radiation!
          ifrad = 0
       elseif(ifmodified.eq.105) then
          ! emulation of Sweigart EOS
          ! modified version with enhancements for interior calculations.
          ifpi = 1
          ! but with no Coulomb and no exchange
          ifcoulomb = 0
          ifexchange_in = 0
          ! and Cody-Thacher F-D integrals
          morder = 1
       elseif(ifmodified.eq.110) then
          ! EOS3 with no radiation pressure.
          ifpi = 1
          ifrad = 0
       elseif(ifmodified.eq.111) then
          ! EOS4 with no radiation pressure.
          ifpi = 1
          ifrad = 0
          ifcoulomb = -15
       elseif(ifmodified.eq.-1) then
          ! original version for comparison with pteh table.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
          ifrad = 0
       elseif(ifmodified.eq.-101) then
          ! original pteh with radiation pressure added (by default)
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
       elseif(ifmodified.eq.-6) then
          ! modified version with enhancements for interior calculations.
          ! ifmodified = 1 defaults
          ! non-linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 12
          ! but with original (via ifmodified, ifpi) pteh
          ! pressure ionization without planck-larkin occupation probability
          ! or excitation
          ifexcited = 0
          ifpi = -1
       elseif(ifmodified.eq.-7) then
          ! modified version with enhancements for interior calculations.
          ! ifmodified = 1 defaults
          ! non-linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 12
          ! but with original (via ifmodified, ifpi) geff
          ! pressure ionization without planck-larkin occupation probability
          ! or excitation
          ifexcited = 0
          ifpi = -2
       elseif(ifmodified.eq.-8) then
          ! modified version with enhancements for interior calculations.
          ! ifmodified = 1 defaults
          ! non-linear, Kapusta, possible strong degeneracy approximation.
          ifexchange_in = 12
          ! but with original (via ifmodified, ifpi) mhd
          ! pressure ionization without planck-larkin occupation probability
          ! or excitation
          ifexcited = 0
          ifpi = -3
       elseif(ifmodified.eq.-30) then
          ! thermodynamic consistency
          iftc = 0
          ! original version for interior calculations.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
       elseif(ifmodified.eq.-31) then
          ! thermodynamic consistency
          iftc = 1
          ! original version for interior calculations.
          ifreducedmass = 0
          ifh2 = 4
          ifh2plus = 0
          morder = 3
          ifcoulomb = -13
          ifpi = -1
       else
          error stop 'free_eos_modern: bad ifmodified for ifoption=1'
       endif
    elseif(ifoption.eq.2) then
       ! fjs style
       if(ifmodified.eq.0) then
          ! original version for interior calculations
          ifreducedmass = 0
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
       elseif(ifmodified.eq.1) then
          ! EOS2
          ifcoulomb = -5
          ifpi = 2
       elseif(ifmodified.eq.11) then
          ! EOS2 without radiation pressure
          ifcoulomb = -5
          ifpi = 2
          ifrad = 0
       elseif(ifmodified.eq.2) then
          ! modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = -5
          ifpi = 2
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
             ifexcited = 2
             nmax = 4
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.3) then
          ! modified version to fit old opal table extension
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifcoulomb = -5
          ifpi = 2
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(lw.gt.2) then
             ifexcited = 1
             nmax = 4
          else
             ifexcited = 0
          endif
       elseif(ifmodified.eq.4) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
       elseif(ifmodified.eq.26) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
          ! special with no excitation
          ifexcited = 0
       elseif(ifmodified.eq.-1) then
          ! original version for comparison with geff code results
          ifreducedmass = 0
          ! but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.e5_fp_kind) then
             ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
          ! version of geff code I have ignores radiation pressure
          ! so mimic this behaviour with the current code.
          ifrad = 0
       elseif(ifmodified.eq.-101) then
          ! original geff with radiation pressure added
          ifreducedmass = 0
          ! but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.e5_fp_kind) then
             ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
          ! version of geff code I have ignores radiation pressure, but
          ! put it in anyway (by default) for this option of the current code.
       elseif(ifmodified.eq.-2) then
          ! original version for comparison with sireff code results
          ! g(eta) exchange with linear transform
          ifexchange_in = -101
          ifreducedmass = 0
          ! but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.e5_fp_kind) then
             ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
          ! version of sireff code I have ignores radiation pressure
          ! so mimic this behaviour with the current code.
          ifrad = 0
       elseif(ifmodified.eq.-201) then
          ! Same as -2 above except for including radiation pressure.
          ! g(eta) exchange with linear transform
          ifexchange_in = -101
          ifreducedmass = 0
          ! but with different H2 t limit (for strict mimicry)
          if(exp(tl).gt.1.e5_fp_kind) then
             ifh2 = 0
          endif
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
          ! version of sireff code I have ignores radiation pressure
          ! but put it in anyway (by default) for this option of the current code.
       elseif(ifmodified.eq.-30) then
          ! thermodynamic consistency
          iftc = 0
          ! original version for interior calculations
          ifreducedmass = 0
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
       elseif(ifmodified.eq.-31) then
          ! thermodynamic consistency
          iftc = 1
          ! original version for interior calculations
          ifreducedmass = 0
          ifh2plus = 0
          morder = 3
          ifcoulomb = -1
          ifpi = -2
       elseif(ifmodified.eq.30) then
          ! thermodynamic consistency
          iftc = 0
          ! modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
       elseif(ifmodified.eq.31) then
          ! thermodynamic consistency
          iftc = 1
          ! modified version with enhancements for interior calculations.
          ifcoulomb = -5
          ifpi = 2
       else
          error stop 'free_eos_modern: bad ifmodified for ifoption=2'
       endif
    elseif(ifoption.eq.3) then
       ! mdh style
       if(lw.gt.2.and.ifmodified.gt.0) then
          ! n.b. use truncated sum because approximations don't work
          ! very well for mhd case.
          ! n.b.  use excited states for H, He, and molecules but
          ! not metals (for now).
          ifexcited = 12
          nmax = 10
       elseif(lw.gt.2.and.ifmodified.le.0) then
          ! n.b. use truncated sum because approximations don't work very well for mhd case
          ! original mode doesn't have molecular electronic excitation
          ifexcited = 11
          nmax = 10
       else
          ifexcited = 0
       endif
       if(ifmodified.eq.0) then
          ! original version for interior calculations.
          ! we know morder=1 (Cody-Thacher Fermi-Dirac without relativistic
          ! correction) is a bad assumption for red giant cores, but the
          ! point is to demonstrate this for the original mhd mode.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
       elseif(ifmodified.eq.1) then
          ! EOS1
          ifcoulomb = 5
          ifpi = 3
       elseif(ifmodified.eq.11) then
          ! EOS1 without radiation pressure
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
       elseif(ifmodified.eq.2) then
          ! modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          ! follow opal table
          nmax = 4
       elseif(ifmodified.eq.102) then
          ! modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ! same as 2 except use Debye-Huckel
          ifcoulomb = 1
          ifpi = 3
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          ! follow opal table
          nmax = 4
       elseif(ifmodified.eq.103) then
          ! modified version to fit old opal table
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 3
          ! use default ifrad behaviour (only change from 2)
          ! ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          ! follow opal table
          nmax = 4
       elseif(ifmodified.eq.20) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! drop h2 and h2+
          ifh2 = 0
          ifh2plus = 0
       elseif(ifmodified.ge.201.and.ifmodified.le.231) then
          ! Same as ifmodified = 1 except for ifred = 0 to match EOS2005 table
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
          ! only difference for 201 <= ifmodified < 210
          if(ifmodified.eq.201) then
             ! non-linear, Kapusta, possible strong degeneracy approximation.
             ifexchange_in = 12
          elseif(ifmodified.eq.202) then
             ! alternative Coulomb
             ifcoulomb = 6
             ! non-linear, Kapusta, possible strong degeneracy approximation.
             ifexchange_in = 12
          elseif(ifmodified.eq.203) then
             ! no exchange
             ifexchange_in = 0
          elseif(ifmodified.eq.204) then
             ! no Coulomb
             ifcoulomb = 0
             ! non-linear, Kapusta, possible strong degeneracy approximation.
             ifexchange_in = 12
          elseif(ifmodified.eq.205) then
             ! no Coulomb or exchange
             ifcoulomb = 0
             ifexchange_in = 0
          elseif(ifmodified.gt.210) then
             ! Standard Coulomb and exchange treatment
             ! but change Fermi-Dirac morder through all variations.
             if(ifmodified.eq.211) then
                morder = 1
             elseif(ifmodified.eq.212) then
                morder = -3
             elseif(ifmodified.eq.213) then
                morder = -5
             elseif(ifmodified.eq.214) then
                morder = -8
             elseif(ifmodified.eq.215) then
                morder = -21
             elseif(ifmodified.eq.216) then
                morder = -23
             elseif(ifmodified.eq.217) then
                morder = 3
             elseif(ifmodified.eq.218) then
                morder = 5
             elseif(ifmodified.eq.219) then
                morder = 8
             elseif(ifmodified.eq.220) then
                morder = 13
             elseif(ifmodified.eq.221) then
                morder = 15
             elseif(ifmodified.eq.222) then
                morder = 18
             elseif(ifmodified.eq.223) then
                morder = 21
             elseif(ifmodified.eq.224) then
                morder = 23
             endif
          endif
       elseif(301.le.ifmodified.and.ifmodified.le.320) then
          ! This series to be compared with ifmodified.eq.1 (ifcoulomb=5, ifpi=3)
          ifcoulomb = 5
          ifpi = 3
          if(ifmodified.eq.301) then
             ! nocoulomb
             ifcoulomb = 0
          elseif(ifmodified.eq.302) then
             ! dh
             ifcoulomb = 1
          elseif(ifmodified.eq.303) then
             ! dh tau
             ifcoulomb = 2
          elseif(ifmodified.eq.304) then
             ! g_pteh(gamma) with everything else the same as ifcoulomb=5
             ifcoulomb = 9
          elseif(ifmodified.eq.311) then
             ! no exchange
             ifexchange_in = 0
          elseif(ifmodified.eq.312) then
             ! normal exchange (14) with linear inversion approximation
             ifexchange_in = 114
          elseif(ifmodified.eq.313) then
             ! normal exchange (14) with non-relativistic approximation
             ifexchange_in = -14
          elseif(ifmodified.eq.314) then
             ! series approximation(s) for exchange
             ifexchange_in = 12
          elseif(ifmodified.eq.315) then
             ! higher-order approximation for exchange
             ifexchange_in = 15
          elseif(ifmodified.eq.316) then
             ! highest-order approximation for exchange
             ifexchange_in = 16
          endif
       elseif(ifmodified.eq.3) then
          ! modified version to fit old opal table extension
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 3
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          ! follow opal table
          nmax = 4
       elseif(ifmodified.eq.4) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! but with no exchange
          ifexchange_in = 0
       elseif(ifmodified.eq.24) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! but with no Coulomb
          ifcoulomb = 0
       elseif(ifmodified.eq.26) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! but with no excitation
          ifexcited = 0
       elseif(ifmodified.eq.30.or.ifmodified.eq.1030) then
          ! thermodynamic consistency
          iftc = 0
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          if(ifmodified.eq.1030) morder = 23
       elseif(ifmodified.eq.31.or.ifmodified.eq.1031) then
          ! thermodynamic consistency
          iftc = 1
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          if(ifmodified.eq.1031) morder = 23
       elseif(ifmodified.eq.43) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! but with special nmax
          nmax = 4
       elseif(ifmodified.eq.44) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! but with special nmax
          nmax = 25
       elseif(ifmodified.eq.45) then
          ! modified version with enhancements for interior calculations.
          ifcoulomb = 5
          ifpi = 3
          ! but with special nmax
          nmax = 50
       elseif(ifmodified.eq.-1) then
          ! original version for interior calculations.
          ! we know morder=1 (Cody-Thacher Fermi-Dirac without relativistic
          ! correction) is a bad assumption for red giant cores, but the
          ! point is to demonstrate this for the original mhd mode.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
          ! but use no H2+
          ifh2plus = 0
       elseif(ifmodified.eq.-30) then
          ! thermodynamic consistency
          iftc = 0
          ! original version for interior calculations.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
       elseif(ifmodified.eq.-31) then
          ! thermodynamic consistency
          iftc = 1
          ! original version for interior calculations.
          morder = 1
          ifcoulomb = 2
          ifpi = -3
       else
          error stop 'free_eos_modern: bad ifmodified for ifoption=3'
       endif
    elseif(ifoption.eq.4) then
       if(ifmodified.eq.2) then
          ! modified version to fit opal table for Saumon style
          ifh2plus = 0
          morder = 1
          ifcoulomb = 5
          ifpi = 4
          ifrad = 0
          ! non-linear transform of G(eta) exchange
          ifexchange_in = -1
          if(tl.lt.tllim) then
             ifexcited = 12
             nmax = 4
          else
             lw = 2
             ifexcited = 0
          endif
       elseif(ifmodified.eq.3) then
          ! modified version to fit He Saumon tables
          ifh2 = 0
          ifh2plus = 0
          morder = 1
          ifexchange_in = 0
          ! He table calculated with SOCP for log10 rho > 0.5, otherwise DH
          ! theory.  For He keep comparison less than log10 rho = 0.5 and use
          ! DH theory.
          ifcoulomb = 1
          ifpi = -4
          ifrad = 0
          ifexcited = 0
          if(tl.ge.tllim) then
             lw = 2
          endif
       elseif(ifmodified.eq.4 .or.ifmodified.eq.30.or.ifmodified.eq.31) then
          ! modified version to fit H Saumon tables
          ! No H2+
          ifh2plus = 0
          ! Cody-Thacher non-relativistic approximation to Fermi-Dirac
          ! integrals.
          morder = 1
          ! Linear, Kapusta term, lowest-order EFF-style
          ! approximation for J and K integrals.
          ifexchange_in = 114
          ! Best Coulomb treatment.
          ifcoulomb = 5
          ! Fit of Saumon et al pressure ionization with Planck-Larkin dropped.
          ifpi = -4
          ifrad = 0
          if(tl.ge.tllim) then
             ! Full ionization (which also excludes all molecules) forced
             ! for all elements in this case because of deficiencies in the
             ! current pressure-ionization formulation. This constraint
             ! causes a discontinuity (especially at higher densities) at
             ! tllim which makes the SCVH option suite unsuitable for
             ! stellar-interior calculations.  This should all be
             ! reassessed (FIXME) when I implement and calibrate an
             ! improved pressure-ionization formulation which includes
             ! hard-sphere interaction radii for all non-bare positive ions
             ! such as He+, C+, ....
             lw = 2
             ifexcited = 0
          else
             ! Drop excited H2 treatment which requires H2+
             ifexcited = 11
             ! Explicit summation to large principal quantum number
             ! required since the adjustment of the MDH-style pressure
             ! ionization to fit the Saumon et al results does not quench
             ! high principal quantum numbers nearly as much as the other
             ! MDH-style pressure-ionization treatments.
             nmax = 100
          endif
          if(ifmodified.eq.31) iftc = 1
       else
          error stop 'free_eos_modern: bad ifmodified for ifoption=4'
       endif
    elseif(ifoption.eq.5) then
       ! ifpi = 5 (Planck-Larkin pressure ionization, only) options
       if(ifmodified.eq.1) then
          ! same as standard (ifoption.eq.1.and.ifmodified.eq.1) except
          ! use Planck-Larkin only.
          if(tl.ge.tllim) then
             ifcoulomb = -15
             ! all elements approximated as fully ionized.
             lw = 2
             ifmtrace = 1
             ifexcited = 0
             ifpi = 0
          else
             ifcoulomb = 5
             ifpi = 5
          endif
       elseif(ifmodified.eq.2) then
          morder = 1
          ifcoulomb = 5
          ifpi = 5
          ifrad = 0
       else
          error stop 'free_eos_modern: bad ifmodified for ifoption=5'
       endif
    else
       error stop 'free_eos_modern: bad ifoption'
    endif
    if(tl.ge.tllim) then
       ! turn off molecules for ln T > tllim where the limit usually
       ! corresponds to T = 1.d6.
       ! the partition functions are approximated above T = 1.d5
       ! by a Taylor series approach which assures continuity for
       ! second-order thermodynamic quanties (e.g. grada).
       ! This approximation is well behaved (at least), and the
       ! errors in it don't matter very much because n(H2) and n(H2+)
       ! are so small for these temperatures (for the densities
       ! associated with stars).  Eventually the Taylor series would
       ! underflow (negative curvature parabola in ln Q vs ln T in
       ! all cases) so that is why we turn it off above t=10^6.
       ! (also code is more efficient).  We have also tried turning
       ! off above T = 10^5, but leaves tiny but noticable discontinuity.
       ifh2 = 0
       ifh2plus = 0
    endif
    ! N.B. Use array sections for all arrays below to convert input user arrays to
    ! the exact sizes required by free_eos_detailed and the routines that it calls.
    call free_eos_detailed(&
         verbosity, ifh2, ifh2plus, morder, ifexchange_in,&
         ifmtrace, ifcoulomb, ifpi, ifrad, lw, ifexcited, nmax,&
         ifreducedmass, iftc, ifmodified, kif,&
         eps(:neps), match_variable, fl, tl, fm, ft,&
         t, rho, rl, p, pl, cf, cp, sf, st, grada, rtp,&
         rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta,&
         gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2,&
         local_degeneracy(:nderivp1), local_pressure(:nderivp1), local_density(:nderivp1),&
         local_energy(:nderivp1), local_enthalpy(:nderivp1), local_entropy(:nderivp1),&
         iteration_count, info)

    fl_old = fl
    ! special call to test thermodynamic consistency
    if(abs(ifmodified).eq.30.or.abs(ifmodified).eq.31.or.abs(ifmodified).eq.1030.or.abs(ifmodified).eq.1031) then
       if(local_entropy(2).eq.0._fp_kind) then
          pl = 0._fp_kind
       else
          ! this coefficient is normally negative
          pl = log(abs(local_entropy(2)))
       endif
    endif
    qe = local_energy(1)
    ! Sometimes because of optional internal debug testing
    ! local_density(1) contains some value other than ln rho.
    ! Also although overflows and underflows are guarded against for
    ! rho, it is always possible some have slipped through that logic.  Thus, for
    ! both these cases limit the following exp argument limits
    ! to avoid underflows and overflows when calculating qv = 1/rho.
    qv = 1._fp_kind/exp(max(ln_underflow_limit, min(ln_overflow_limit,local_density(1))))
    s = local_entropy(1)

    if(present(degeneracy)) then
       if(nderivp1.gt.size(degeneracy)) error stop 'optional degeneracy array argument has an insufficient size'
       degeneracy = local_degeneracy
    endif
    if(present(pressure)) then
       if(nderivp1.gt.size(pressure)) error stop 'optional pressure array argument has an insufficient size'
       pressure = local_pressure
    endif
    if(present(density)) then
       if(nderivp1.gt.size(density)) error stop 'optional density array argument has an insufficient size'
       density = local_density
    endif
    if(present(energy)) then
       if(nderivp1.gt.size(energy)) error stop 'optional energy array argument has an insufficient size'
       energy = local_energy
    endif
    if(present(enthalpy)) then
       if(nderivp1.gt.size(enthalpy)) error stop 'optional enthalpy array argument has an insufficient size'
       enthalpy = local_enthalpy
    endif
    if(present(entropy)) then
       if(nderivp1.gt.size(entropy)) error stop 'optional entropy array argument has an insufficient size'
       entropy = local_entropy
    endif
  end subroutine free_eos_modern

end module mod_free_eos
