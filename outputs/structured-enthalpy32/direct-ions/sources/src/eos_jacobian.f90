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

!> The purpose of this eos_jacobian subroutine is to calculate the
!> components of the linear system of equations that must be solved in
!> calling routines in order for them to calculate one step in the NR
!> (Newton-Raphson) iterative solution of the EOS (equation of state)
!> that is equivalent to (but typically much faster) than Helmholtz
!> free-energy minimization if the starting solution for the NR
!> iteration is close enough to a minimum of that form of
!> thermodynamic energy.  See <a
!> href="http://freeeos.sourceforge.net/solution.pdf">"The NR Solution
!> Paper"</a> for further details.
!>
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
!> \param[in] old_allow_log old_allow_log(naux)
!>   The first n_partial_aux compact indices of old_allow_log
!>   contain logical flags which show whether it is possible to take
!>   the ln of aux_old(iaux), where iaux is the corresponding
!>   non-compact auxiliary variable index.  That is the corresponding
!>   old auxiliary variable is positive and not the one (iaux = 4)
!>   corresponding to rl = ln (rho) which is already in ln form.<br>
!> \param[in] old_allow_log_neg old_allow_log_neg(naux)
!>   The first n_partial_aux compact indices of old_allow_log_neg
!>   contain logical flags which show whether it is possible to take
!>   the ln of *the negative* of aux_old(iaux), where iaux is the corresponding
!>   non-compact auxiliary variable index.  That is the corresponding
!>   old auxiliary variable is negative and not the one (iaux = 4)
!>   corresponding to rl = ln (rho) which is already in ln form.<br>
!> \param[out] old_new_allow_log old_new_allow_log(naux)
!>   The calculated first n_partial_aux compact indices of
!>   old_new_allow_log contain logical flags which show whether
!>   it is possible to take the ln of *both* aux_old(iaux) and
!>   aux(iaux), where iaux is the corresponding non-compact auxiliary
!>   variable index.  That is both the corresponding old and new
!>   auxiliary variables are positive and not the ones (iaux = 4)
!>   corresponding to rl = ln (rho) which is already in ln form.<br>
!> \param[out] old_new_allow_log_neg old_new_allow_log_neg(naux)
!>   The calculated first n_partial_aux compact indices of
!>   old_new_allow_log_neg contain logical flags which show whether
!>   it is possible to take the ln of *both the negative values* of aux_old(iaux) and
!>   aux(iaux), where iaux is the corresponding non-compact auxiliary
!>   variable index.  That is both the corresponding old and new
!>   auxiliary variables are negative and not the ones (iaux = 4)
!>   corresponding to rl = ln (rho) which is already in ln form.<br>
!> \param[in] partial_aux partial_aux(naux) is a vector whose first
!>   n_partial_aux elements maps from compact to non-compact auxiliary
!>   variable indices for the n_partial_aux auxiliary variables that
!>   are relevant to the adopted free-energy model.<br>
!> \param[in] ifrad
!>   ifrad = 0 means radiation pressure is excluded from the
!>     free-energy model and for kif = 1 the input match_variable is
!>     consistent with this model, i.e., it is ln P excluding
!>     radiation pressure.<br>
!>   ifrad = 1 means radiation pressure is included in the free-energy
!>     model and for kif = 1 the input match_variable is consistent
!>     with this model, i.e., it is ln P including radiation
!>     pressure.<br>
!>   ifrad = 2 is a special case where radiation pressure is included
!>     in the free-energy model (i.e., all output quantities are
!>     calculated with radiation included.), but for kif = 1 (the only
!>     kif value allowed for ifrad = 2) the input match_variable is
!>     *not* consistent with this model, i.e., it is ln(ptotal-prad).
!>     This special case is used to reduce significance loss in
!>     regions which are dominated by radiation pressure.<br>
!> \param[in] match_variable match_variable is one of the two independent variables of the EOS
!>   (with the other one being tl).  The definition of match_variable depends on kif and ifrad,
!>   and for kif > 0 the internal value of the quantity corresponding to match_variable is
!>   iteratively adjusted (by changing fl) to match match_variable.<br>
!>   For kif = 0 match_variable is fl, the EFF degeneracy parameter.<br>
!>   For kif = 1 match_variable is the ln of the total pressure -
!>   radiation pressure (if ifrad = 0 or 2) or the ln of the total
!>   pressure (ifrad = 1).<br>
!>   For kif = 2, match_variable is the ln of the density.<br>
!> \param[in] kif kif (and ifrad) determine the definition of
!>   match_variable.  See the documentation of that argument for
!>   further details.<br>
!> \param[out] fl fl is the calculated EFF degeneracy parameter which for kif = 0 is
!>   determined directly from match_variable and for kif > 0 is
!>   iteratively determined by matching match_variable.  See the
!>   documentation of that argument for further details.<br>
!> \param[in] aux_old
!>   aux_old(naux) contains (using non-compact indices with only the
!>   subset of auxiliary variables defined that are required by the
!>   free-energy model) the input (i.e., old) non-traditional form of
!>   auxiliary variables.  These variables help to determine the
!>   non-ideal equilibrium constants and therefore the dissociation
!>   and ionization balance (i.e., the number densities of all species
!>   included in the free-energy model) and ultimately a new set of
!>   auxiliary variables that depend on those number densities.<br>
!> \param[out] aux
!>   aux(naux) contains the calculated non-traditional form of auxiliary
!>   variables which the call to traditional_to_aux copies from the traditional
!>   form of auxiliary variables supplied by eos_calc.  Note that only
!>   parts of aux are filled in by that routine depending on the
!>   adopted free-energy model, i.e., just the subset of the traditional
!>   auxiliary variables that has been calculated by eos_calc.<br>
!> \param[out] auxf
!>   auxf(naux) contains the calculated non-traditional form of the partial derivatives of auxiliary
!>   variables wrt fl which the call to traditional_to_aux copies from the traditional
!>   form of those partial derivatives supplied by eos_calc.<br>
!> \param[out] auxt
!>   auxt(naux) contains the calculated non-traditional form of the partial derivatives of auxiliary
!>   variables wrt ln t which the call to traditional_to_aux copies from the traditional
!>   form of those partial derivatives supplied by eos_calc.<br>
!> \param[out] aux_dv
!>   aux_dv(nionsp2, naux) contains the calculated non-traditional form of the
!>   partial derivatives of auxiliary variables wrt the elements of
!>   dv which the call to traditional_to_aux copies from the traditional form of
!>   those partial derivatives supplied by eos_calc.  Note the first
!>   index of aux_dv is in compact form with a defined size of max_index while the second index is in
!>   non-compact form.<br>
!> \param[in] njacobian njacobian is the number of linear equations that
!>   must be solved to determine the next Newton-Raphson iterative
!>   step (i.e., the iterative change in aux_old plus (if kif > 0) the
!>   change in fl that would reduce faux to the zero vector if faux
!>   were a linear function of aux_old and fl).  Therefore the
!>   njacobian value indicates the number of faux values (where faux
!>   is the RHS vector corresponding to the linear system of
!>   equations) that must be calculated and the number of rows and
!>   columns of jacobian (the Jacobian matrix corresponding to the
!>   linear system of equations) that must be calculated.<br>
!> \param[out] faux faux(naux) contains the calculated RHS vector of the Newton-Raphson
!>   linear system of equations that must be solved for each NR
!>   iterative step toward the EOS solution.  The number of linear
!>   equations to be solved is njacobion so this routine only
!>   calculates the first njacobian values of faux.<br>
!> \param[out] jacobian jacobian(naux, naux) is the calculated Jacobian matrix of the
!>   Newton-Raphson linear system of equations that must be solved for
!>   each NR iterative step toward the EOS solution.  This matrix
!>   contains the negative of the partial derivatives of the
!>   components of faux wrt the components of aux_old and (if kif > 0)
!>   fl.  The number of linear equations to be solved is njacobion so
!>   this routine only calculates the values of the njacobian by
!>   njacobian submatrix of jacobian.<br>
!> \param[out] p p is the calculated total pressure corresponding to the free-energy
!>   model.<br>
!> \param[in] pr pr is the radiative component of the total pressure
!>   corresponding to the free-energy model.<br>
!> \param[out] dv_aux dv_aux(nionsp2, naux) contains the calculated partial
!>   derivatives of the equilibrium constant vector, dv(nionsp2), wrt
!>   to the auxiliary variable vector aux(naux).  Note both indices of
!>   dv_aux are compact so the defined size of dv_aux is max_index by n_partial_aux.<br>
!> \param[out] lambda lambda is the calculated Coulomb interaction parameter.  Note
!>   lambda has two different definitions depending on ifcoulomb.<br>
!> \param[out] gamma_e gamma_e is the calculated Coulomb diffraction parameter.<br>
!> \param[in] ifnr
!>   ifnr = 0 means calculate partial derivatives of nuvar (communicated via
!>     the mod_nuvar module), hne, sum0, sum2, s, extrasum, h_ion, sh,
!>     he_ion, he_ion2, h_neutral, hequil, and sumpl1 wrt f, t and
!>     sumpl0 wrt f with no other variables fixed, i.e., use input dvf
!>     and dvt and chain rule.  N.B. chain of calling routines depend
!>     on the assertion that ifnr = 0 means no *_dv variables should
!>     be read or written.<br>
!>   ifnr = 1 means calculate partial derivatives of nuvar (communicated via
!>     the mod_nuvar module), hne, sum0, sum2, extrasum, h_ion,
!>     he_ion, and he_ion2 wrt dv with f, t fixed.  n.b. list of
!>     variables is subset of those listed for ifnr = 0 because we
!>     only need dv partial derivatives for variables used to calculate
!>     (output) auxiliary variables.<br>
!>   ifnr = 2 (should not occur).<br>
!>   ifnr = 3 same as combination of ifnr = 0 and ifnr = 1.
!>     Note this is quite different in detail from ifnr = 3
!>     interpretation for many other routines, but the general
!>     motivation is the same for all ifnr = 3 results; calculate
!>     f, t, and other partial derivatives.<br>
!> \param[in] nux nux is the number density (divided by rho*NA) of the
!>   maximum possible ionization electrons for hydrogen.<br>
!> \param[in] nuy nuy is the number density (divided by rho*NA) of the
!>   maximum possible ionization electrons for helium.<br>
!> \param[in] nuz nuz is the sum of the number densities (divided by
!>   rho*NA) of the maximum possible ionization electrons for all
!>   metal (i.e., other than hydrogen or helium) elements.<br>
!> \param[in] n_partial_ions n_partial_ions is the maximum index used for
!>   *monatomic ions" in partial_ions where n_partial_ions <=
!>   max_index <= n_partial_ions+2 depending on whether H2 and H2+ are
!>   included in the free-energy model.<br>
!> \param[in] sum0_mc When if_mc is 1, sum0_mc is defined as the component
!>   of the fully ionized nu form of sum0 for metals which themselves
!>   are approximated as fully ionized and which are therefore not
!>   accounted for by hcon_mc.<br>
!> \param[in] sum2_mc When if_mc is 1, sum2_mc is defined as the component
!>   of the fully ionized nu form of sum2 for metals which themselves
!>   are approximated as fully ionized and which are therefore not
!>   accounted for by hecon_mc.<br>
!> \param[in] hcon_mc When if_mc is 1, hcon_mc is defined as the component
!>   of the fully ionized nu form of sum0 for metals which themselves
!>   are accounted for without approximation (i.e., as partially
!>   ionized) and which are therefore not accounted for by sum0_mc
!>   After this calculation, hcon_mc is further transformed by
!>   dividing it by eps(1).<br>
!> \param[in] hecon_mc When if_mc is 1, hecon_mc is defined as the
!>   component of the fully ionized nu form of sum2 for metals which
!>   themselves are accounted for without approximation (i.e., as
!>   partially ionized) and which are therefore not accounted for by
!>   sum2_mc.  After this calculation, hecon_mc is transformed by
!>   subtracting untransformed hcon_mc and dividing that quantity by
!>   eps(2).<br>
!> \param[in] partial_ions partial_ions(nions+2) is a vector whose first max_index elements map from
!>   a compact to non-compact ion index.<br>
!> \param[in] f f is the EFF degeneracy parameter.<br>
!> \param[in] eta eta is the Cox and Guili degeneracy parameter which
!>   is related to f by<br>
!>   wf = sqrt(1.d0 + f)<br>
!>   eta = fl+2.d0*(wf-log(1.d0+wf)).<br>
!> \param[in] wf wf is d eta/d fl = sqrt(1.d0 + f).<br>
!> \param[in] t t is the temperature on the Kelvin temperature scale.<br>
!> \param[in] n_e n_e is the number density of free electrons.<br>
!> \param[in] pstar
!>   pstar(nstar+6) (where nstar is 3) is an array containing a scaled
!>     fermi-dirac integral and its 8 partial derivatives.<br>
!>   pstar(1) is defined such that p_e = cpe*pstar(1),
!>     where p_e is the partial pressure of free electrons.<br>
!>   pstar(2) through pstar(9) contain the partial derivatives of
!>     ln pstar wrt fl, tcl, fl fl, fl tcl, tcl tcl, fl fl fl, fl fl tcl,
!>     and fl tcl tcl; where fl is ln f and tcl is ln beta = ln kt/mc^2.<br>
!> \param[in] dve_exchange dve_exchange is the *negative* chemical potential of electron/kT from the exchange effect.<br>
!> \param[in] dve_exchangef dve_exchangef is the partial derivative of dve_exchange wrt fl.<br>
!> \param[in] dve_exchanget dve_exchanget is the partial derivative of dve_exchange wrt ln t.<br>
!> \param[in] full_sum0 full_sum0 is the full ionization approximation to sum0/rho*NA,
!>   where sum0 is the sum over positive ion number densities, rho is
!>   the density, and NA is Avogadro's number.<br>
!> \param[in] full_sum1 full_sum1 is the full ionization approximation to sum1/rho*NA,
!>   where sum1 is the sum over positive ion number densities weighted
!>   by the ion charge, rho is the density, and NA is Avogadro's
!>   number.<br>
!> \param[in] full_sum2 full_sum2 is the full ionization approximation to sum2/rho*NA,
!>   where sum2 is the sum over positive ion number densities weighted
!>   by the square of the ion charge, rho is the density, and NA is
!>   Avogadro's number.<br>
!> \param[in] charge charge(nionsp2) contains the charge of the ions in ion order.<br>
!> \param[in] charge2 charge2(nionsp2) contains the square of the charge of the ions in ion order.<br>
!> \param[in] dv_pl dv_pl(nionsp2) contains the Planck-Larkin change in equilibrium constant in species order.<br>
!> \param[in] dv_plt dv_plt(nionsp2) contains the partial derivative of dv_pl wrt ln t.<br>
!> \param[in] ifdv ifdv(nionsp2) only has its last two elements defined and
!>   those elements contain a control flag to help decide whether to
!>   calculate dv for H2 and H2+.
!>   If the hydrogen abundance and ifh2 are both greater than zero
!>   then ifdv(nionsp2-1) is 1 and therefore dv for H2 will be
!>   calculated, but otherwise ifdv(nionsp2-1) is 0 and therefore dv
!>   for H2 will not be calculated.
!>   If the hydrogen abundance and ifh2plus are both greater than zero
!>   then ifdv(nionsp2) is 1 and therefore dv for H2+ will be
!>   calculated, but otherwise ifdv(nionsp2) is 0 and therefore dv
!>   for H2+ will not be calculated.<br>
!> \param[in] ifdvzero ifdvzero = .true. implies dvzero (the quantity
!>   added to all dv in ionize) is zero.  This is the low-temperature
!>   option.  For high temperatures, ifdvzero is .false., and dvzero
!>   is a zero point shift that renders the pressure ionization
!>   contribution to dv of all bare ions zero.  This option gives the
!>   smallest significance loss near full ionization and high
!>   densities.<br>
!> \param[in] ifcoulomb_mod
!>   ifcoulomb_mod = mod(|ifcoulomb|,10) where ifcoulomb is documented in free_eos_detailed.<br>
!>   ifcoulomb_mod = 0 ignore Coulomb interaction.<br>
!>   ifcoulomb_mod = 1 use Debye-Huckel Coulomb approximation.<br>
!>   ifcoulomb_mod = 2 use Debye-Huckel Coulomb approximation with tau(x) correction.<br>
!>   ifcoulomb_mod = 3 use PTEH Coulomb approximation with their theta_e.<br>
!>   ifcoulomb_mod = 4 use PTEH Coulomb approximation with fermi-dirac theta_e.<br>
!>   ifcoulomb_mod = 5 use DH smoothly connected to modified OCP
!>     with DeWitt definition of lambda.<br>
!>   ifcoulomb_mod = 6 same as 5 with alternative smooth connection.<br>
!>   ifcoulomb_mod = 7 DH (Gamma < 1) or OCP using new DeWitt lambda.<br>
!>   ifcoulomb_mod = 8 same as 7 with theta_e = 0.<br>
!>   ifcoulomb_mod = 9 same as 4 with DeWitt definition of lambda
!>     (using sum0a = sum0 + ne*theta_e).<br>
!> \param[in] if_dc
!>   if_dc = 0 means do not use diffraction correction for the Coulomb effect.<br>
!>   if_dc = 1 means do use diffraction correction for the Coulomb effect.<br>
!> \param[in] ifsame_zero_abundances
!>  ifsame_zero_abundances = 1 means this call has the same zero abundance pattern as the last call.<br>
!>  ifsame_zero_abundances /= 1 means this call has a different zero abundance pattern then the last call.<br>
!> \param[in] ifexcited
!>   ifexcited > 0 means use excited states (must have Planck-Larkin or ifpi = 3 or 4).<br>
!>    0 < ifexcited < 10 means use approximation to explicit summation.<br>
!>   10 < ifexcited < 20 means use explicit summation.<br>
!>   mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules) and helium.<br>
!>   mod(ifexcited,10) = 2 same as 1 + H2 and H2+.<br>
!>   mod(ifexcited,10) = 3 same as 2 + partially ionized metals.<br>
!> \param[in] ifsame_under
!>   ifsame_under .eqv. .true. means use same underflow zeroing for each
!>   ion as in the previous call.  This option is useful for removing
!>   small discontinuities caused by variations in the underflow
!>   zeroing which sometimes foil the last stages of
!>   convergence.<br>
!>   ifsame_under .eqv. .false. means calculate underflow zeroing.<br>
!> \param[in] if_free_non_ideal_calc if_free_non_ideal_calc
!>   .eqv. .true. means this eos_jacobian call is made via the
!>   eos_bfgs call chain and therefore free_non_ideal_calc needs to be
!>   called.<br>
!>   Otherwise, free_non_ideal_calc should not be called for
!>   efficiency (and to avoid all ifnr values other than 3).<br>
!> \param[in] inv_aux
!>   inv_aux(naux) is a vector which maps from non-compact to compact
!>   auxiliary variable indices.  Elements of inv_aux that are zero
!>   mark non-compact auxiliary variable indices corresponding to
!>   auxiliary variables that are not relevant to the current
!>   free-energy model.<br>
!> \param[in] inv_ion
!>   inv_ion(nions_inp2) maps ion index to contiguous ion index used
!>     in NR iteration dot products.<br>
!> \param[in] max_index max_index is the actual maximum index used for
!>   partial_ions where n_partial_ions <= max_index <=
!>   n_partial_ions+2.<br>
!> \param[in] partial_elements
!>   partial_elements(n_partial_elements+2) index of elements treated
!>     as partially ionized consistent with ifelement.<br>
!> \param[in] ion_end
!>   ion_end(nelements_in+2) keeps track of largest ion index for each
!>   element and the largest ion index for H2 (ielement = nelements+1
!>   when H2 is included in the EOS) and H2+ (ielement = nelements+2
!>   when H2+ is included in the EOS).<br>
!> \param[in] mion_end mion_end is the maximum index where ion_end is defined.<br>
!> \param[in] ifionized
!>   ifionized = 0 means all elements treated as partially ionized.<br>
!>   ifionized = 1 means trace metals fully ionized.<br>
!>   ifionized = 2 means all elements treated as fully ionized.<br>
!> \param[in] if_pteh
!>   if_pteh = 1 means the PTEH approximation is used for the Coulomb sums and therefore the
!>     actual auxiliary variable Coulomb sums do not have to be calculated.<br>
!>   if_pteh = 0 means the PTEH approximation is not used for the Coulomb sums and therefore the
!>     actual auxiliary variable Coulomb sums have to be calculated if if_mc is also 0.<br>
!> \param[in] if_mc
!>   if_mc = 1 means the metal Coulomb approximation is used to help calculate the Coulomb sums and therefore the
!>     actual auxiliary variable Coulomb sums do not have to be calculated.<br>
!>   if_mc = 0 means the metal Coulomb approximation is not used to help calculate the Coulomb sums and therefore the
!>     actual auxiliary variable Coulomb sums have to be calculated if if_pteh is also 0.<br>
!> \param[in] ifreducedmass
!>   ifreducedmass = 1 means use reduced mass (which is approximately
!>     equal to the electron mass) in equilibrium constant.<br>
!>   ifreducedmass = 0 means use electron mass in equilibrium constant.<br>
!> \param[in] ifsame_abundances
!>   ifsame_abundances = 1 means that this subroutine is called with
!>     the same abundances as the previous call which saves steps in the
!>     ionize calculations.<br>
!>   ifsame_abundances = 0 means that this subroutine *might* be
!>     called with different abundances than the previous call.  This
!>     value makes ionize safe against any abundance changes but
!>     requires extra steps in the ionize calculations.<br>
!> \param[in] ifmtrace
!>   ifmtrace is used to help keep track of changes in ifelement.
!>     That is, if ifionized is different from the previous call or
!>     ifsame_abundances = 0, or ifmtrace different from the previous
!>     call, then local saved variables ionized_hne, ionized_sum0,
!>     ionized_sum2, and ionized_extrasum(2) are recalculated.<br>
!> \param[in] iatomic_number
!>   iatomic_number(nelements_in) is the atomic number of each
!>     element.  The code logic demands there is 1 neutral monatomic
!>     species and iatomic_number ionized monatomic species for each
!>     element.  For example, for elements without molecules there are
!>     iatomic_number(ielement) ionization potentials (for the neutral
!>     species and all ion species other than the bare ion) for the
!>     element indexed by ielement.<br>
!> \param[in] ifpi_local
!>   ifpi_local = 0, 1, 2, 3, 4 corresponds to using no, pteh, geff, mhd, or
!>     Saumon-style pressure ionization.<br>
!> \param[in] ifpl
!>   ifpl = 0, 1 corresponds to using no or planck-larkin occupation probability.<br>
!> \param[in] ifmodified
!>   ifmodified > 0 or not affects ifpi_fit inside excitation_sum.<br>
!> \param[in] ifh2
!>   ifh2 non-zero implies h2 is included.<br>
!> \param[in] ifh2plus
!>   ifh2plus non-zero implies h2plus is included.<br>
!> \param[in] izlo
!>   izlo is lowest core charge (1 for H, He, 2 for He+, etc.) for all
!>     species.<br>
!> \param[in] izhi
!>   izhi is highest core charge (1 for H, He, 2 for He+, etc.) for all
!>     species.<br>
!> \param[in] bmin
!>   bmin(izhi) is the minimum bion for all species with a given core
!>     charge.<br>
!> \param[in] nmin
!>   nmin(izhi) is the minimum excited principal quantum number for
!>     all species with a given core charge.<br>
!> \param[in] nmin_max
!>   nmin_max(izhi) is the largest minimum excited principal quantum
!>     number for all species with a given core charge.<br>
!> \param[in] nmin_species
!>   nmin_species(nions+2) is the minimum excited principal quantum
!>     number organized by species.<br>
!> \param[in] nmax
!>   nmax is maximum principal quantum number included in sum (to be
!>     compatible with opal which used nmax = 4 rather than infinity).
!>     if(nmax > 300000) then treated as infinity in qryd_approx.
!>     otherwise nmax is meant to be used with qryd_calc only (i.e.,
!>     case for mhd approximations not programmed).<br>
!> \param[in] eps
!>   eps(nelements_in) is an array of neps = nelements_in = 24 values
!>     of relative abundance by weight divided by the appropriate
!>     atomic weight. We use the atomic weight scale where
!>     (un-ionized) C(12) has a weight of 12.00000000....  All weights
!>     are for the un-ionized element. The eps value for an element
!>     should be the sum of the individual isotopic eps values for
!>     that element.  The eps array refers to the elements in the
!>     following order:
!>     H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,A,Ca,Ti,Cr,Mn,Fe,Ni.<br>
!> \param[in] tl tl = log(t).<br>
!> \param[in] tc2
!>   tc2 = c2/t.<br>
!> \param[in] bi
!>   bi(nions_inp2) ionization potentials (cm**-1) in ion order from
!>     neutral to last ion before the bare nucleus ion, e.g., H
!>     (excluding the H+ bare nucleus); He, He+ (excluding the He++
!>     bare nucleus); C through C+++++ (excluding the C++++++ bare
!>     nucleus); etc.  The H2 and H2+ ionization potentials are the
!>     last two in this list.<br>
!> \param[in] h2diss h2diss = h2 dissociation energy (cm^{-1}).<br>
!> \param[in] plop
!>   plop(nions_inp2) is the planck-larkin occupation probability in
!>     ion order (starting with neutral in each case and avoiding bare
!>     nuclei as in bi).<br>
!> \param[in] plopt plopt(nions_inp2) is the first derivative of plop wrt ln t.<br>
!> \param[in] plopt2 plopt2(nions_inp2) is the second derivative of plop wrt ln t.<br>
!> \param[in] r_ion3
!>   r_ion3(nions_inp2) is *the cube of the* effective radii (in ion
!>     order but going from neutral to next to bare ion for each
!>     species) corresponding to the MDH interaction between all but
!>     bare nucleii species and ionized species.  Last two are H2 and
!>     H2+ which are not used in ionize.<br>
!> \param[in] r_neutral
!>   r_neutral(nelements_in+2) is the effective radii corresponding to
!>     MDH neutral-neutral interactions.  Last two are H2 and H2+
!>     which are not used in ionize.<br>
!> \param[in] ifelement
!>   ifelement(nelements_in) is an array whose values are 1 if the
!>     corresponding element is treated as partially ionized
!>     (depending on abundance and trace metal treatment) or 0 if the
!>     corresponding element is treated as fully ionized or has no
!>     abundance.<br>
!> \param[out] dvzero
!>   dvzero(nelements_in) contains the calculated zero-point shifts (zero for
!>     hydrogen) that must be added to dv in all uses below.
!>     Sometimes (high density, high ionization) some of these
!>     quantities can be larger than 1.e9_fp_kind so if we avoid using
!>     it in differences of dv quantities we can gain many significant
!>     digits.<br>
!> \param[out] dv
!>   dv(nions_inp2) contains the calculated equilibrium constants (explained above).<br>
!> \param[out] dvf
!>   dvf(nions_inp2) contains the calculated partial derivatives of dv wrt fl.<br>
!> \param[out] dvt
!>   dvt(nions_inp2) contains the calculated partial derivatives of dv wrt ln t.<br>
!> \param[in] nion
!>   nion(nions_inp2) is an array of ion charge in ion order (must be same order
!>     as bi) e.g., for H+, He+, He++, etc.
!> \param[in] rhostar
!>   rhostar(nstar+6) (where nstar is 3) is an array containing a scaled
!>     fermi-dirac integral and its 8 partial derivatives.<br>
!>   rhostar(1) is defined such that n_e/NA = cd*rhostar(1),
!>     where n_e is the number density of free electrons.<br>
!>   rhostar(2) through rhostar(9) contain the partial derivatives of
!>     ln rhostar wrt fl, tcl, fl fl, fl tcl, tcl tcl, fl fl fl, fl fl tcl,<br>
!>     and fl tcl tcl; where fl is ln f and tcl is ln beta = ln kt/mc^2.<br>
!> \param[out] ne ne is the calculated nu form (n/(rho*NA)) of number density of positive ions.<br>
!> \param[out] nef nef is the calculated partial derivative of ne wrt fl.<br>
!> \param[out] net net is the calculated partial derivative of ne wrt ln t.<br>
!> \param[out] sion
!>   R*sion is the calculated entropy per unit mass corresponding to
!>     the ideal component of the free energy.  n.b. for the
!>     completely ionized case sion is completely calculated in the
!>     ionize routine that is called by eos_calc.  For the partially
!>     ionized case without molecular formation, the ionize routine
!>     separates the hydrogen component of sion from the rest of the
!>     elements (see the sh argument of ionize).  For the partially
!>     ionized case with H2 and possibly H2+ formation, sh returned by
!>     ionize is ignored and instead that hydrogen component is
!>     calculated properly in this routine for the case given and
!>     added to sion calculated by ionize for the remaining
!>     elements.<br>
!>   Units check: From the formula below, sion has the same units as
!>     nu = n/(rho*NA), i.e., the inverse of the atomic mass.
!>     Therefore, R sion = k*NA sion has units of energy per degree K
!>     (from Boltzman's constant) per unit cgs mass, i.e., the
!>     expected cgs units for the entropy per unit mass.<br>
!>   (FIXME(2021), the code to implement the entropy below has been
!>     independently confirmed to be consistent with the ideal free
!>     energy component for the non-electron species.  However, the
!>     dated entropy formulas documented below do not appear be
!>     consistent with the code and have terms added to the -5/2 term
!>     that appear not to be dimensionless!  So this documentation
!>     needs to be rewritten starting with the formula for the ideal
!>     free energy component for the non-electron species.<br>
!>     From the ideal component of the free energy sion is the sum over
!>     all non-electron species of<br>
!>       -nu_i*[-5/2 - 3/2 ln T + ln(alpha) - ln Na +<br>
!>       ln (n_i/[A_i^{3/2} Q_i]) - dln Q_i/d ln T],<br>
!>       where nu_i = n_i/(Na rho),<br>
!>       Na is the Avogadro number,<br>
!>       A_i is the atomic weight,<br>
!>       alpha =  (2 pi k/[Na h^2])^(-3/2) Na
!>       (see documentation of alpha^2 in the mod_free_eos_constants
!>       module), and<br>
!>       Q_i is the internal ideal partition function of the
!>       non-Rydberg states.  (Currently, we calculate this partition
!>       function by the statistical weights of the ground states of
!>       monatomic hydrogen and helium and the combined lower states
!>       (roughly approximated) of each of the metal.  This is a
!>       superb approximation for hydrogen, helium is subsequently
!>       corrected for detailed excitation of the non-Rydberg states,
!>       and this crude approximation for the metals is not currently
!>       corrected.  Thus, in all *current* monatomic cases ln Q_i is
!>       a constant, and dln Q_i/d ln T is zero, but this will be
!>       updated for the metals eventually.)<br>
!>       From the equilibrium constant approach we have<br>
!>       -nu_i ln (n_i/[A_i^{3/2} Q_i]) =<br>
!>       -nu_i * [ln (n_neutral/[A_neutral^{3/2} Q_neutral) +<br>
!>       (- chi_i/kT + dv_i)].<br>
!>       If we sum this term over all species of an element we obtain
!>       s_element = - eps * ln (eps*n_neutral/sum(n))<br>
!>       - eps ln (alpha/(A_neutral^{3/2} Q_neutral)) -<br>
!>       - sum over all species of the element of nu_i*(-chi/kT + dv(i))<br>
!>       where we have ignored the term<br>
!>       eps * (-5/2 - 3/2 ln T + ln rho)<br>
!>       (taking into account the first 4 terms above)
!>       until a final correction of the sion zero point
!>       in free_eos_detailed.<br>
!> \param[out] sionf sionf is the calculated partial derivative of sion wrt fl.<br>
!> \param[out] siont siont is the calculated partial derivative of sion wrt ln t.<br>
!> \param[out] uion
!>   c2*cr*uion is the calculated internal energy per unit mass corresponding to
!>     the the ideal component of the free energy.<br>
!>   Units check: uion has the same units as energy (in cm^{-1}) times
!>     nu = n/(rho*NA).  The factor h*c converts energy in cm^{-1} to
!>     energy in ergs.  Therefore, c2*cr*uion = h*c*NA*uion has the
!>     expected cgs units of energy per unit mass.<br>
!> \param[out] pnorad pnorad is the calculated (non-ideal) pressure with the radiation pressure subtracted.<br>
!> \param[out] pnoradf pnoradf is the calculated partial derivative of pnorad wrt fl.<br>
!> \param[out] pnoradt pnoradt is the calculated partial derivative of pnorad wrt ln t.<br>
!> \param[out] pnorad_aux pnorad_aux(n_partial_aux+2) contains the
!> calculated partial derivatives of pnorad wrt aux, fl, and ln t.<br>
!> \param[out] free free is the calculated (non-ideal) free energy.<br>
!> \param[out] free_aux free_aux(n_partial_aux+2) contains the
!> calculated partial derivatives of free wrt aux, fl, and ln t.
!> N.B. That ln t partial derivative is currently unused so it just
!> contains a 0._fp_kind placeholder.<br>
!> \param[in] nextrasum nextrasum is the size of the extrasum and related
!>   arrays that will be locally allocated as traditional auxiliary
!>   variables used by eos_calc and copied back to the appropriate
!>   parts of aux and related non-traditional auxiliary variables by
!>   traditional_to_aux.  nextrasum is either 0 or 9 depending on the
!>   adopted free-energy model, i.e., depending on whether the Boolean
!>   value of (ifpi_local.eq.3.or.ifpi_local.eq.4) which is assigned
!>   to the local variable ifextrasum is either .true. or .false.<br>
!> \param[out] sumpl1
!>   sumpl1 (returned only if ifpl = 1) is the calculated weighted sum over non-H
!>     nu(i) = n(i)/(rho*NA) with weight of w + d ln w/d ln T, where w
!>     is the Planck-Larkin occupation probability.<br>
!> \param[out] sumpl1f
!>  sumpl1f is the calculated partial derivative of sumpl1 wrt fl.<br>
!> \param[out] sumpl1t
!>  sumpl1t is the calculated partial derivative of sumpl1 wrt ln t.<br>
!> \param[out] sumpl2
!>   sumpl2 (returned only if ifpl = 1) is the calculated weighted sum over non-H
!>     nu(i) = n(i)/(rho*NA) with weight of d ln w/d ln T, where w is
!>     the Planck-Larkin occupation probability.<br>
!> \param[out] h2
!>   h2 is the calculated number density (optionally in n form or else
!>     in nu form where nu = n/(rho*avogadro)) of H2.<br>
!> \param[out] h2f h2f is the calculated partial derivative of h2 wrt fl.<br>
!> \param[out] h2t h2t is the calculated partial derivative of h2 wrt ln t.<br>
!> \param[out] h2_dv h2_dv(nions_inp2) is an array of the calculated
!>   partial derivatives of h2 wrt the elements of dv.  This array has a compact
!>   index whose maximum defined size is max_index.<br>
!> \param[out] info
!>   info is the calculated return code for the procedure which is
!>     zero for success and non-zero (with offset to identify which
!>     procedure the error occurred in) if an error occurred within
!>     the procedure or within some procedure it called.<br>
!> \param[in] debug_aux_dv debug_aux_dv controls whether (.true.) or not
!>   (.false.) to debug selected parts of the aux_dv matrix result for
!>   partial derivatives of the aux vector wrt the dv vector by using
!>   numerical differences to approximate aux_dv.  Only one of
!>   debug_aux_dv, debug_dv_aux, or debug_jacobian can be .true. at
!>   the same time.<br>
!> \param[in] debug_dv_aux debug_dv_aux controls whether (.true.) or not
!>   (.false.) to debug selected parts of the dv_aux matrix, dvf
!>   vector, and dvt vector results for partial derivatives of the dv
!>   vector wrt aux, fl, and tl by using numerical differences to
!>   approximate those partial derivatives.  Only one of debug_aux_dv,
!>   debug_dv_aux, or debug_jacobian can be .true. at the same
!>   time.<br>
!> \param[in] debug_jacobian debug_jacobian controls whether (.true.) or
!>   not (.false.) to debug the jacobian matrix (which contains the
!>   negative partial derivative of components of faux wrt the
!>   components of aux_old) or the partial derivatives of components
!>   of faux wrt fl or tl by using numerical differences to
!>   approximate those partial derivatives.  Only one of debug_aux_dv,
!>   debug_dv_aux, or debug_jacobian can be .true. at the same
!>   time.<br>
!> \param[in] itemp1 itemp1 indicates the first independent variable where
!>   derivatives wrt that variable are being debugged.<br>
!>   When debug_aux_dv is .true. the partial derivatives aux(1:n_partial_aux, dv) wrt to the component
!>   of dv corresponding to itemp1 are being debugged.<br>
!>   When debug_dv_aux is .true. we have the following cases:<br>
!>       When itemp1 <= n_partial_aux, the partial derivative dv(1:max_index, aux_old, fl, tl)
!>         wrt the component of aux_old corresponding to itemp1 is being debugged.<br>
!>       When itemp1 == n_partial_aux+1, the partial derivative dv(1:max_index, aux_old, fl, tl)
!>       wrt fl is being debugged.<br>
!>       When itemp1 == n_partial_aux+2, the partial derivative dv(1:max_index, aux_old, fl, tl)
!>       wrt tl is being debugged.<br>
!>   When debug_jacobian is .true. we have the following cases:<br>
!>       When itemp1 <= n_partial_aux, the partial derivative faux(1:njacobian, aux_old, fl, tl)
!>         wrt the component of aux_old corresponding to itemp1 is being debugged.<br>
!>       When itemp1 == n_partial_aux+1, the partial derivative faux(1:njacobian, aux_old, fl, tl)
!>       wrt fl is being debugged.<br>
!>       When itemp1 == n_partial_aux+2, the partial derivative faux(1:njacobian, aux_old, fl, tl)
!>       wrt tl is being debugged.<br>
!> \param[in] itemp2 itemp2 indicates the second independent variable where
!>   derivatives wrt that variable are being debugged.<br>
!>   When debug_aux_dv is .true. the partial derivatives aux(1:n_partial_aux, dv) wrt to the component
!>   of dv corresponding to itemp2 are being debugged.<br>
!>   When debug_dv_aux is .true. we have the following cases:<br>
!>       When itemp2 <= n_partial_aux, the partial derivative dv(1:max_index, aux_old, fl, tl)
!>         wrt the component of aux_old corresponding to itemp2 is being debugged.<br>
!>       When itemp2 == n_partial_aux+1, the partial derivative dv(1:max_index, aux_old, fl, tl)
!>       wrt fl is being debugged.<br>
!>       When itemp2 == n_partial_aux+2, the partial derivative dv(1:max_index, aux_old, fl, tl)
!>       wrt tl is being debugged.<br>
!>   When debug_jacobian is .true. we have the following cases:<br>
!>       When itemp2 <= n_partial_aux, the partial derivative faux(1:njacobian, aux_old, fl, tl)
!>         wrt the component of aux_old corresponding to itemp2 is being debugged.<br>
!>       When itemp2 == n_partial_aux+1, the partial derivative faux(1:njacobian, aux_old, fl, tl)
!>       wrt fl is being debugged.<br>
!>       When itemp2 == n_partial_aux+2, the partial derivative faux(1:njacobian, aux_old, fl, tl)
!>       wrt tl is being debugged.<br>
!> \param[in] deltav1 When debug_aux_dv is .true. deltav1 is the
!>   increment to dv(partial_ions(itemp1)) that is used to to help
!>   debug aux_dv using numerical differences.<br>
!> \param[in] deltav2 When debug_aux_dv is .true. deltav2 is the
!>   increment to dv(partial_ions(itemp1)) that is used to to help
!>   debug aux_dv using numerical differences.<br>
!> \param[out] debug_results debug_results(3, max_index or njacobian+2) contains
!>   calculated debugging information whenever one (and only one!) of
!>   debug_aux_dv, debug_dv_aux, or debug_jacobian is .true.<br>
subroutine eos_jacobian(&
     verbosity, old_allow_log, old_allow_log_neg, old_new_allow_log, old_new_allow_log_neg, partial_aux,&
     ifrad,&
     match_variable, kif, fl,&
     aux_old, aux, auxf, auxt, aux_dv,&
     njacobian,&
     faux, jacobian, p, pr,&
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
     pnorad, pnoradf, pnoradt, pnorad_aux,&
     free, free_aux,&
     nextrasum,&
     sumpl1, sumpl1f, sumpl1t, sumpl2,&
     h2, h2f, h2t, h2_dv, info,&
     debug_aux_dv, debug_dv_aux, debug_jacobian, itemp1, itemp2, deltav1, deltav2, debug_results)

  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  integer, intent(in) ::&
       verbosity, partial_aux(:),&
       ifrad,&
       kif, njacobian,&
       ion_end(:),&
       ifnr,&
       mion_end,&
       n_partial_ions,&
       partial_ions(:),&
       ifdv(:), ifcoulomb_mod, if_dc,&
       ifsame_zero_abundances,&
       ifexcited,&
       inv_aux(:),&
       inv_ion(:), max_index,&
       partial_elements(:),&
       ifionized,&
       if_pteh, if_mc,&
       ifreducedmass,&
       ifsame_abundances,&
       ifmtrace,&
       iatomic_number(:),&
       ifpi_local,&
       ifpl, ifmodified, ifh2, ifh2plus,&
       izlo, izhi,&
       nmin(:), nmin_max(:),&
       nmin_species(:), nmax,&
       ifelement(:),&
       nion(:),&
       nextrasum

  integer, intent(out) :: info

  logical, intent(in) :: ifdvzero, ifsame_under, if_free_non_ideal_calc

  logical, intent(in) :: old_allow_log(:), old_allow_log_neg(:)
  logical, intent(out) :: old_new_allow_log(:), old_new_allow_log_neg(:)

  ! For debugging
  logical, intent(in), optional :: debug_aux_dv, debug_dv_aux, debug_jacobian

  integer, intent(in), optional :: itemp1, itemp2

  real(fp_kind), intent(in), optional :: deltav1, deltav2

  real(fp_kind), intent(out), optional :: debug_results(:,:)

  real(fp_kind), intent(in) ::&
       match_variable, tl,&
       aux_old(:),&
       pr,&
       nux, nuy, nuz,&
       sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
       f, eta, wf, t, n_e,&
       pstar(:),&
       dve_exchange, dve_exchangef, dve_exchanget,&
       full_sum0, full_sum1, full_sum2,&
       charge(:), charge2(:),&
       dv_pl(:), dv_plt(:),&
       bmin(:),&
       eps(:), tc2,&
       bi(:), h2diss,&
       plop(:), plopt(:), plopt2(:),&
       r_ion3(:), r_neutral(:),&
       rhostar(:)

  real(fp_kind), intent(out) ::&
       fl, aux(:), auxf(:), auxt(:), aux_dv(:,:),&
       faux(:), jacobian(:,:), p,&
       dv_aux(:,:),&
       lambda, gamma_e,&
       dvzero(:),&
       dv(:), dvf(:), dvt(:),&
       ne, nef, net,&
       sion, sionf, siont, uion,&
       pnorad, pnoradf, pnoradt, pnorad_aux(:),&
       free, free_aux(:),&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       h2, h2f, h2t, h2_dv(:)

  ! Internal variables
  logical debug_present

  integer, parameter :: rl_aux_index = 4
  integer, parameter :: size_three = 3
  integer, parameter :: nstar = 9
  ! maximum core charge for non-bare ion
  integer, parameter :: maxcorecharge = 28
  integer naux, nelements, nionsp2, n_partial_aux
  !integer, parameter :: maxfjs_aux = 4
  integer, parameter :: nxextrasum = 4
  integer index_aux, index
  integer iaux, jndex_aux, iaux1, iaux2
  integer, parameter :: maxnextrasum_local = 9
  integer jtemp, jaux, jion

  !real(fp_kind)&
       ! dve, dvef, dvet, dve0, dve2,&
       ! dve_pi, dve_pif, dve_pit,&
       ! dve_coulomb, dve_coulombf, dve_coulombt,&
       ! dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,&
  !real(fp_kind)&
       ! rho, rho_new,&
       ! dvh, dvhf, dvht,&
       ! dvhe1, dvhe1f, dvhe1t,&
       ! dvhe2, dvhe2f, dvhe2t,&
       ! dsum0, dsum2,&
       ! nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,&
  real(fp_kind)&
       ! pe_cgs,&
       ! p0, pion, pionf, piont,&
       ! pex, pexf, pext,&
       ! ppi, ppif, ppit, ppir,&
       ! pexcited, pexcitedf, pexcitedt,&
       ! dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2,&
       ! fion, fionf,&
       fion, fionf,&
       sumpl0, sumpl0f,&
       deriv_lnp_factor

  real(fp_kind), allocatable ::&
       fion_dv(:),&
       sumpl0_dv(:),&
       aux_aux(:,:)

  naux = size(partial_aux)
  nelements = size(ion_end) - 2
  nionsp2 = size(partial_ions)
  n_partial_aux = size(pnorad_aux) - 2
  debug_present =&
       present(debug_aux_dv).and.&
       present(debug_dv_aux).and.&
       present(debug_jacobian).and.&
       present(itemp1).and.&
       present(itemp2).and.&
       present(deltav1).and.&
       present(deltav2).and.&
       present(debug_results)

  ! sanity checks.
  if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3)) error stop 'eos_jacobian: ifnr must be 0, 1, or 3'
  if(ifnr.eq.0) error stop 'eos_jacobian: ifnr = 0 is disabled'
  if(ifrad.eq.2.and.kif.ne.1) error stop 'eos_jacobian: for ifrad == 2, kif values other than 1 are disabled.'

  if(nextrasum.gt.maxnextrasum_local) error stop 'eos_jacobian: nextrasum too large'
  if(n_partial_aux.gt.naux) error stop 'eos_jacobian: n_partial_aux too large'
  if(&
       naux.ne.size(inv_aux).or.&
       naux.ne.size(old_allow_log).or.&
       naux.ne.size(old_allow_log_neg).or.&
       naux.ne.size(old_new_allow_log).or.&
       naux.ne.size(old_new_allow_log_neg).or.&
       naux.ne.size(aux_old).or.&
       naux.ne.size(aux).or.&
       naux.ne.size(auxf).or.&
       naux.ne.size(auxt).or.&
       naux.ne.size(aux_dv,2).or.&
       naux.ne.size(faux).or.&
       naux.ne.size(jacobian,1).or.&
       naux.ne.size(jacobian,2).or.&
       naux.ne.size(dv_aux,2))&
       error stop 'eos_jacobian: inconsistent naux array sizes'
  if(&
       nelements.ne.size(iatomic_number).or.&
       nelements.ne.size(ifelement).or.&
       nelements.ne.size(eps).or.&
       nelements+2.ne.size(r_neutral).or.&
       nelements.ne.size(dvzero))&
       error stop 'eos_jacobian: inconsistent nelements array sizes'
  if(&
       nionsp2.ne.size(ifdv).or.&
       nionsp2.ne.size(inv_ion).or.&
       nionsp2.ne.size(nmin_species).or.&
       nionsp2.ne.size(nion).or.&
       nionsp2.ne.size(charge).or.&
       nionsp2.ne.size(charge2).or.&
       nionsp2.ne.size(dv_pl).or.&
       nionsp2.ne.size(dv_plt).or.&
       nionsp2.ne.size(bi).or.&
       nionsp2.ne.size(plop).or.&
       nionsp2.ne.size(plopt).or.&
       nionsp2.ne.size(plopt2).or.&
       nionsp2.ne.size(r_ion3).or.&
       nionsp2.ne.size(aux_dv,1).or.&
       nionsp2.ne.size(dv_aux,1).or.&
       nionsp2.ne.size(dv).or.&
       nionsp2.ne.size(dvf).or.&
       nionsp2.ne.size(dvt).or.&
       nionsp2.ne.size(h2_dv))&
       error stop 'eos_jacobian: inconsistent nionsp2 array sizes'
  if(&
       maxcorecharge.ne.size(nmin).or.&
       maxcorecharge.ne.size(nmin_max).or.&
       maxcorecharge.ne.size(bmin))&
       error stop 'eos_jacobian: inconsistent maxcorecharge array sizes'
  if(&
       nstar.ne.size(pstar).or.&
       nstar.ne.size(rhostar))&
       error stop 'eos_jacobian: inconsistent nstar array sizes'
  if(n_partial_aux+2.ne.size(free_aux))&
       error stop 'eos_jacobian: inconsistent n_partial_aux array sizes'

  if(debug_present) then
     ! First dimension of debug_results must be size_three, and second dimension
     ! must be the compact size of the vector of functions whose partial derivatives
     ! are being calculated.
     if(size_three.ne.size(debug_results,1))&
          error stop 'eos_jacobian: wrong size for first dimension of debug_results'
     if(debug_dv_aux) then
        if(max_index.ne.size(debug_results,2))&
             error stop 'eos_jacobian: debug_dv_aux: wrong size for second dimension of debug_results'
     elseif(debug_aux_dv) then
        if(n_partial_aux.ne.size(debug_results,2))&
             error stop 'eos_jacobian: debug_aux_dv: wrong size for second dimension of debug_results'
     elseif(debug_jacobian) then
        if(njacobian.ne.size(debug_results,2))&
             error stop 'eos_jacobian: debug_jacobian: wrong size for second dimension of debug_results'
     endif
  endif
  allocate(&
       fion_dv(nionsp2),&
       sumpl0_dv(nionsp2),&
       aux_aux(naux, naux)&
       )

  if(if_taint_allocated_real) then
     call taint_allocated_real(fion_dv)
     call taint_allocated_real(sumpl0_dv)
     call taint_allocated_real(aux_aux)
  endif

  call eos_warm_step(&
       verbosity, old_allow_log, old_allow_log_neg, aux_old, partial_aux,&
       dv_aux,&
       lambda, gamma_e,&
       ifnr, nux, nuy, nuz,&
       n_partial_ions,&
       nextrasum, nxextrasum,&
       sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
       partial_ions, f, eta, wf, t, n_e, pstar,&
       dve_exchange, dve_exchangef, dve_exchanget,&
       full_sum0, full_sum1, full_sum2, charge, charge2,&
       dv_pl, dv_plt,&
       ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
       ifsame_zero_abundances,&
       ifexcited, ifsame_under,&
       inv_aux,&
       inv_ion, max_index,&
       partial_elements, ion_end, mion_end,&
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
       pnorad, pnoradf, pnoradt, pnorad_aux,&
       fion, fionf, fion_dv,&
       sumpl0, sumpl0f, sumpl0_dv,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       h2, h2f, h2t, h2_dv,&
       aux, auxf, auxt, aux_dv,&
       info,&
       debug_aux_dv, itemp1, itemp2, deltav1, deltav2)
  if(info.ne.0) return

  if(debug_present) then
     if(debug_dv_aux) then
        if(1.le.itemp1.and.itemp1.le.n_partial_aux) iaux1 = partial_aux(itemp1)
        if(1.le.itemp2.and.itemp2.le.n_partial_aux) iaux2 = partial_aux(itemp2)
        ! N.B. debug_results has a compact second index while dv, dvf, and dvt have corresponding
        ! non-compact indices.
        do jtemp = 1, max_index
           jion = partial_ions(jtemp)
           debug_results(1, jtemp) = dv(jion)
           if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
              debug_results(2, jtemp) = dv_aux(jtemp,itemp1)*aux_old(iaux1)
           elseif(itemp1.eq.n_partial_aux+1) then
              debug_results(2, jtemp) = dvf(jion)
           elseif(itemp1.eq.n_partial_aux+2) then
              debug_results(2, jtemp) = dvt(jion)
           endif
           if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
              debug_results(3, jtemp) = dv_aux(jtemp,itemp2)*aux_old(iaux2)
           elseif(itemp2.eq.n_partial_aux+1) then
              debug_results(3, jtemp) = dvf(jion)
           elseif(itemp2.eq.n_partial_aux+2) then
              debug_results(3, jtemp) = dvt(jion)
           endif
        enddo
        return
     elseif(debug_aux_dv) then
        ! N.B. debug_results has a compact second index while aux and
        ! second index of aux_dv have corresponding non-compact
        ! indices.
        do jtemp = 1, n_partial_aux
           jaux = partial_aux(jtemp)
           debug_results(1, jtemp) = aux(jaux)
           debug_results(2, jtemp) = aux_dv(itemp1,jaux)
           debug_results(3, jtemp) = aux_dv(itemp2,jaux)
        enddo
        return
     endif
  endif

  do jndex_aux = 1,n_partial_aux
     jaux = partial_aux(jndex_aux)
     ! n.b. auxiliary variable corresponding to jaux = rl_aux_index is already in
     ! log form (i.e., rl = log(rho)).
     ! If element of old_new_allow_log is .true. can take log of
     ! corresponding element of both old and new aux.
     old_new_allow_log(jndex_aux) = jaux.ne.rl_aux_index.and.aux_old(jaux).gt.0._fp_kind.and.aux(jaux).gt.0._fp_kind
     ! If element of old_new_allow_log_neg is .true. can take log of
     ! corresponding *negative* element of both old and new aux.
     old_new_allow_log_neg(jndex_aux) = jaux.ne.rl_aux_index.and.aux_old(jaux).lt.0._fp_kind.and.aux(jaux).lt.0._fp_kind
  enddo
  ! form RHS and Jacobian
  do jndex_aux = 1,n_partial_aux
     jaux = partial_aux(jndex_aux)
     ! RHS that will be transformed to change in aux from aux_old.
     if(old_new_allow_log(jndex_aux)) then
        faux(jndex_aux) = log(aux(jaux))
     elseif(old_new_allow_log_neg(jndex_aux)) then
        faux(jndex_aux) = log(-aux(jaux))
     else
        faux(jndex_aux) = aux(jaux)
     endif
     do index_aux = 1,n_partial_aux
        ! negative partial aux(jaux) wrt aux_old(iaux)
        ! n.b. dv_aux(index,index_aux) has compact index for
        ! index and index_aux and
        ! aux_dv(index,jaux)  has compact index for
        ! index and *uncompact*index for jaux
        jacobian(jndex_aux,index_aux) = -dot_product(aux_dv(1:max_index,jaux), dv_aux(1:max_index,index_aux))
        ! partial aux(jaux) wrt aux_old(iaux)
        aux_aux(jndex_aux,index_aux) = -jacobian(jndex_aux,index_aux)
        ! transform to derivative of faux = log(aux) or log(-aux)
        if(old_new_allow_log(jndex_aux).or.old_new_allow_log_neg(jndex_aux)) then
           jacobian(jndex_aux,index_aux) = jacobian(jndex_aux,index_aux)/aux(jaux)
        else
           jacobian(jndex_aux,index_aux) = jacobian(jndex_aux,index_aux)
        endif
        ! transform to derivative with respect to log(auxold) or
        ! log(-auxold)
        iaux = partial_aux(index_aux)
        if(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)) then
           jacobian(jndex_aux,index_aux) = jacobian(jndex_aux,index_aux)*aux_old(iaux)
        endif
     enddo
     if(old_new_allow_log(jndex_aux)) then
        faux(jndex_aux) = faux(jndex_aux) - log(aux_old(jaux))
     elseif(old_new_allow_log_neg(jndex_aux)) then
        faux(jndex_aux) = faux(jndex_aux) - log(-aux_old(jaux))
     else
        faux(jndex_aux) = faux(jndex_aux) - aux_old(jaux)
     endif
     ! take negative derivative wrt to log(auxold), log(-auxold)
     ! or auxold whichever is appropriate
     jacobian(jndex_aux,jndex_aux) = jacobian(jndex_aux,jndex_aux) + 1._fp_kind
     ! negative partial of aux(new) wrt fl with aux(old) fixed.
     if(njacobian.gt.n_partial_aux) then
        jacobian(jndex_aux,njacobian) = -auxf(jaux)
        ! transform to derivative of faux = log(aux) or log(-aux)
        if(old_new_allow_log(jndex_aux).or.old_new_allow_log_neg(jndex_aux)) then
           if((jaux.ne.rl_aux_index.and.(aux(jaux).eq.0._fp_kind.or.aux_old(jaux).eq.0._fp_kind))) then
              ! This branch should not be taken because of above definition
              ! of old_new_allow_log variables referring to both aux_old and aux.
              error stop 'eos_jacobian: bad logic'
           else
              jacobian(jndex_aux,njacobian) = jacobian(jndex_aux,njacobian)/aux(jaux)
           endif
        else
           jacobian(jndex_aux,njacobian) = jacobian(jndex_aux,njacobian)
        endif
     endif
  enddo
  if(njacobian.gt.n_partial_aux) then
     if(kif.eq.1) then
        p = pnorad + pr
        ! Additional RHS function that must be zeroed
        if(ifrad.gt.1) then
           if(pnorad.gt.0._fp_kind) then
              faux(njacobian) = match_variable - log(pnorad)
              deriv_lnp_factor = 1._fp_kind/pnorad
           else
              faux(njacobian) = exp(match_variable) - pnorad
              deriv_lnp_factor = 1._fp_kind
           endif
        else
           if(p.gt.0._fp_kind) then
              faux(njacobian) = match_variable - log(p)
              deriv_lnp_factor = 1._fp_kind/p
           else
              faux(njacobian) = exp(match_variable) - p
              deriv_lnp_factor = 1._fp_kind
           endif
        endif
        ! negative partial derivative of faux wrt fl.
        jacobian(njacobian,njacobian) = deriv_lnp_factor*pnoradf
        ! negative partial derivative of faux wrt old
        ! auxiliary variables.
        do index_aux = 1,n_partial_aux
           ! pnorad_aux already with respect to log(auxold) or log(-auxold)
           ! if appropriate, see eos_warm_start.
           jacobian(njacobian,index_aux) = pnorad_aux(index_aux)
           ! convert to ln(pnorad) or ln(p) derivatives
           jacobian(njacobian,index_aux) = deriv_lnp_factor*jacobian(njacobian,index_aux)
        enddo
     elseif(kif.eq.2) then
        ! Additional RHS function that must be zeroed
        faux(njacobian) = match_variable - aux(rl_aux_index)
        ! negative partial derivative of faux wrt fl.
        jacobian(njacobian,njacobian) = auxf(rl_aux_index)
        ! negative partial derivative of faux wrt old
        ! auxiliary variables.
        do index_aux = 1,n_partial_aux
           jacobian(njacobian,index_aux) = dot_product(aux_dv(1:max_index,rl_aux_index), dv_aux(1:max_index,index_aux))
           ! take derivative with respect to log(auxold) or log(-auxold)
           iaux = partial_aux(index_aux)
           if(old_new_allow_log(index_aux).or.old_new_allow_log_neg(index_aux)) then
              jacobian(njacobian,index_aux) = jacobian(njacobian,index_aux)*aux_old(iaux)
           endif
        enddo
     else
        error stop 'eos_jacobian: bad internal logic'
     endif
  endif
  ! FIXME(2022)? I am not sure whether should use old_new_allow_log or
  ! old_allow_log in the argument list below.  The argument for the
  ! first is consistency in aux derivatives, the argument for the
  ! second is the free_aux derivatives in question are wrt aux_old.
  ! More analysis needed to sort out which is correct.  I have stuck
  ! with the first interpretation for now because that is what the old
  ! version of the code did (before the split of allow_log into
  ! old_allow_log and old_new_allow_log).
  if(if_free_non_ideal_calc) &
       call free_non_ideal_calc(&
       verbosity, ifpi_local, ifnr, ifmodified,&
       ifdvzero, inv_aux, partial_ions, inv_ion,&
       partial_elements,&
       ion_end, mion_end,&
       ifdv, ifrad,&
       r_ion3, nion, r_neutral,&
       tl, f, wf, eta, rhostar, pstar,&
       nux, nuy, nuz,&
       if_pteh, full_sum0, full_sum1, full_sum2,&
       if_mc, hcon_mc, hecon_mc, sum0_mc, sum2_mc,&
       ifcoulomb_mod, if_dc,&
       ifexcited, ifpl, ifsame_zero_abundances, ifh2, ifh2plus,&
       izlo, izhi, bmin, nmin, nmin_max, nmin_species, nmax, bi,&
       plop, plopt, plopt2,&
       sumpl0, sumpl0f, sumpl0_dv,&
       max_index,&
       old_new_allow_log, old_new_allow_log_neg, aux_old, partial_aux,&
       dv_aux,&
       nextrasum, aux, auxf, auxt, aux_dv, aux_aux,&
       fion, fionf, fion_dv,&
       ne, nef,&
       free, free_aux)
  if(debug_present) then
     if(debug_jacobian) then
        do jtemp = 1, njacobian
           jaux = partial_aux(jtemp)
           debug_results(1, jtemp) = faux(jtemp)
           ! For each of the jtemp  values, and if itemp[12] in the appropriate range
           ! transform auxf and auxt from partials of aux_new(dv, fl, tl) to
           ! partials of aux_new(dv(aux_old, fl, tl), fl, tl) =
           ! partials aux_new(dv, fl, tl) + partial aux_new(dv, fl, tl)/partial dv * partials dv(aux_old, fl, tl)
           ! Have not implemented partial of faux(njacobian) wrt fl or tl which implies condition on jtemp.
           if(jtemp.le.n_partial_aux.and.itemp1.gt.n_partial_aux.or.itemp2.gt.n_partial_aux) then
              do index = 1, max_index
                 ! first index of aux_dv is compact, second index is not.
                 ! indices of dvf and dvt are non-compact.
                 auxf(jaux) = auxf(jaux) + aux_dv(index, jaux) * dvf(partial_ions(index))
                 auxt(jaux) = auxt(jaux) + aux_dv(index, jaux) * dvt(partial_ions(index))
              enddo
           endif
           if(1.le.itemp1.and.itemp1.le.n_partial_aux) then
              debug_results(2, jtemp) = -jacobian(jtemp,itemp1)
           elseif(jtemp.le.n_partial_aux.and.itemp1.eq.n_partial_aux+1) then
              if(old_new_allow_log(jtemp).or.old_new_allow_log_neg(jtemp)) then
                 debug_results(2, jtemp) = auxf(jaux)/aux(jaux)
              else
                 debug_results(2, jtemp) = auxf(jaux)
              endif
           elseif(jtemp.le.n_partial_aux.and.itemp1.eq.n_partial_aux+2) then
              if(old_new_allow_log(jtemp).or.old_new_allow_log_neg(jtemp)) then
                 debug_results(2, jtemp) = auxt(jaux)/aux(jaux)
              else
                 debug_results(2, jtemp) = auxt(jaux)
              endif
           endif

           if(1.le.itemp2.and.itemp2.le.n_partial_aux) then
              debug_results(3, jtemp) = -jacobian(jtemp,itemp2)
           elseif(jtemp.le.n_partial_aux.and.itemp2.eq.n_partial_aux+1) then
              if(old_new_allow_log(jtemp).or.old_new_allow_log_neg(jtemp)) then
                 debug_results(3, jtemp) = auxf(jaux)/aux(jaux)
              else
                 debug_results(3, jtemp) = auxf(jaux)
              endif
           elseif(jtemp.le.n_partial_aux.and.itemp2.eq.n_partial_aux+2) then
              if(old_new_allow_log(jtemp).or.old_new_allow_log_neg(jtemp)) then
                 debug_results(3, jtemp) = auxt(jaux)/aux(jaux)
              else
                 debug_results(3, jtemp) = auxt(jaux)
              endif
           endif
        enddo
     endif
     return
  endif
end subroutine eos_jacobian
