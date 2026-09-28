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

!> The purpose of this eos_warm_step subroutine is to calculate a new
!> set of auxiliary variables and their derivatives as a function of
!> an old set of auxiliary variables and temperature.  The calculation
!> proceeds as follows:
!>
!> The old set of auxiliary variables is used to calculate general
!> (i.e., non-ideal) equilibrium constants for each species included
!> in a given free-energy model.  Those equilibrium constants
!> determine the complete set of number densities for those species.
!> Finally, a new set of auxiliary variables (and their derivatives)
!> are calculated from those number densities.
!>
!> Note the only routine that calls this one is eos_jacobian. See the
!> documentation for that routine for a description of how all these
!> calculations are used to support the the NR (Newton-Raphson)
!> iterative solution of the EOS (equation of state) that is
!> equivalent to (but typically much faster) than Helmholtz
!> free-energy minimization if the starting solution for the NR
!> iteration is close enough to a minimum of that form of
!> thermodynamic energy.
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
!> \param[in] old_allow_log
!>   The first n_partial_aux compact indices of old_allow_log(naux)
!>   contain logical flags which show whether it is possible to take
!>   the ln of aux_old(iaux), where iaux is the corresponding
!>   non-compact auxiliary variable index.  That is the corresponding
!>   old auxiliary variable is positive and not the one (iaux = 4)
!>   corresponding to rl = ln (rho) which is already in ln form.<br>
!> \param[in] old_allow_log_neg
!>   The first n_partial_aux compact indices of old_allow_log_neg(naux)
!>   contain logical flags which show whether it is possible to take
!>   the ln of *the negative* of aux_old(iaux), where iaux is the corresponding
!>   non-compact auxiliary variable index.  That is the corresponding
!>   old auxiliary variable is negative and not the one (iaux = 4)
!>   corresponding to rl = ln (rho) which is already in ln form.<br>
!> \param[in] aux_old
!>   aux_old(naux) contains (using non-compact indices with only the
!>   subset of auxiliary variables defined that are required by the
!>   free-energy model) the input (i.e., old) non-traditional form of
!>   auxiliary variables.  These variables help to determine the
!>   non-ideal equilibrium constants and therefore the dissociation
!>   and ionization balance (i.e., the number densities of all species
!>   included in the free-energy model) and ultimately a new set of
!>   auxiliary variables that depend on those number densities.<br>
!> \param[in] partial_aux partial_aux(naux) is a vector whose first
!>   n_partial_aux elements maps from compact to non-compact auxiliary
!>   variable indices for the n_partial_aux auxiliary variables that
!>   are relevant to the adopted free-energy model.<br>
!> \param[out] dv_aux dv_aux(nionsp2, naux) contains the calculated partial
!>   derivatives of the equilibrium constant vector, dv(nionsp2), wrt
!>   to the auxiliary variable vector aux(naux).  Note both indices of
!>   dv_aux are compact.<br>
!> \param[out] lambda lambda is the calculated Coulomb interaction
!>   parameter.  Note lambda has two different definitions depending
!>   on ifcoulomb.<br>
!> \param[out] gamma_e gamma_e is the calculated Coulomb diffraction
!>   parameter.<br>
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
!> \param[in] nextrasum nextrasum is the size of the extrasum and related
!>   arrays that will be locally allocated as traditional auxiliary
!>   variables used by eos_calc and copied back to the appropriate
!>   parts of aux and related non-traditional auxiliary variables by
!>   traditional_to_aux.  nextrasum is either 0 or 9 depending on the
!>   adopted free-energy model, i.e., depending on whether the Boolean
!>   value of (ifpi_local.eq.3.or.ifpi_local.eq.4) which is assigned
!>   to the local variable ifextrasum is either .true. or .false.<br>
!> \param[in] nxextrasum nxextrasum is the size of the xextrasum and related
!>   arrays that will be locally allocated as traditional auxiliary
!>   variables used by eos_calc and copied back to the appropriate
!>   parts of aux and related non-traditional auxiliary variables by
!>   traditional_to_aux.  nxextrasum should always be 4.<br>
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
!> \param[in] dve_exchange dve_exchange is the *negative* chemical
!> potential of electron/kT from the exchange effect.<br>
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
!> \param[in] dv_pl dv_pl(nionsp2) contains the Planck-Larkin change
!> in equilibrium constant in species order.<br>
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
!>  ifsame_zero_abundances /= 1 means this call has a different zero
!>  abundance pattern then the last call.<br>
!> \param[in] ifexcited
!>   ifexcited > 0 means use excited states (must have Planck-Larkin or ifpi = 3 or 4).<br>
!>    0 < ifexcited < 10 means use approximation to explicit summation.<br>
!>   10 < ifexcited < 20 means use explicit summation.<br>
!>   mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules) and helium.<br>
!>   mod(ifexcited,10) = 2 same as 1 + H2 and H2+.<br>
!>   mod(ifexcited,10) = 3 same as 2 + partially ionized metals.<br>
!> \param[in] ifsame_under
!>   ifsame_under .eqv. .true. means use same underflow zeroing for each
!>     ion as in previous call.  This option is useful for removing
!>     small discontinuities caused by variations in the underflow
!>     zeroing which sometimes foil the last stages of
!>     convergence.<br>
!>   ifsame_under .eqv. .false. means calculate underflow zeroing.<br>
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
!>   if_mc = 1 means the metal Coulomb approximation is used to help
!>     calculate the Coulomb sums and therefore the actual auxiliary
!>     variable Coulomb sums do not have to be calculated.<br>
!>   if_mc = 0 means the metal Coulomb approximation is not used to
!>     help calculate the Coulomb sums and therefore the actual
!>     auxiliary variable Coulomb sums have to be calculated if
!>     if_pteh is also 0.<br>
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
!>     That is, if ifionized is different from previous call or
!>     ifsame_abundances = 0, or ifmtrace different from previous
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
!> \param[out] pnorad pnorad is the calculated (non-ideal) pressure
!> with the radiation pressure subtracted.<br>
!> \param[out] pnoradf pnoradf is the calculated partial derivative of pnorad wrt fl.<br>
!> \param[out] pnoradt pnoradt is the calculated partial derivative of pnorad wrt ln t.<br>
!> \param[out] pnorad_aux pnorad_aux(n_partial_aux+2) contains the
!> calculated partial derivatives of pnorad wrt aux.  Note this vector
!> has a compact index.<br>
!> \param[out] fion
!>   fion (returned only if ifpi = 3 or 4) is the calculated ion
!>     free_energy/(rho*kT*NA) (= tc2*u - s at equilibrium).<br>
!> \param[out] fionf fionf is the calculated partial derivative of fion wrt fl.<br>
!> \param[out] fion_dv fion_dv(nions_inp2) is an array of the
!> calculated partial derivatives of fion wrt the elements of dv.<br>
!> \param[out] sumpl0
!>   sumpl0 (returned only if ifpl = 1) is the calculated weighted sum over non-H
!>     nu(i) = n(i)/(rho*NA) with weight of w where w is the Planck-Larkin
!>     occupation probability.<br>
!> \param[out] sumpl0f
!>  sumpl0f is the calculated partial derivative of sumpl0 wrt fl.<br>
!> \param[out] sumpl0_dv
!>  sumpl0_dv(nions_inp2) contains the calculated partial derivatives
!>  of sumpl0 wrt the elements of dv.<br>
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
!>   partial derivatives of h2 wrt the elements of dv.<br>
!> \param[out] aux
!>   aux(naux) contains the non-traditional form of auxiliary
!>   variables which traditional_to_aux copies from the traditional
!>   form of auxiliary variables supplied by eos_calc.  Note that only
!>   parts of aux are filled in by that routine depending on the
!>   adopted free-energy model, i.e., what subset of the traditional
!>   auxiliary variables has been calculated by eos_calc.<br>
!> \param[out] auxf
!>   auxf(naux) contains the non-traditional form of the partial derivatives of auxiliary
!>   variables wrt fl which traditional_to_aux copies from the traditional
!>   form of those partial derivatives supplied by eos_calc.<br>
!> \param[out] auxt
!>   auxt(naux) contains the non-traditional form of the partial derivatives of auxiliary
!>   variables wrt ln t which traditional_to_aux copies from the traditional
!>   form of those partial derivatives supplied by eos_calc.<br>
!> \param[out] aux_dv
!>   aux_dv(nionsp2, naux) contains the non-traditional form of the
!>   partial derivatives of auxiliary variables wrt the elements of
!>   dv which traditional_to_aux copies from the traditional form of
!>   those partial derivatives supplied by eos_calc.  Note the first
!>   index of aux_dv is in compact form while the second index is in
!>   non-compact form.<br>
!> \param[out] info
!>   info is the calculated return code for the procedure which is
!>     zero for success and non-zero (with offset to identify which
!>     procedure the error occurred in) if an error occurred within
!>     the procedure or within some procedure it called.<br>
!> \param[in] debug_aux_dv debug_aux_dv controls whether (.true.) or not
!>   (.false.) to debug the aux_dv vector result by using numerical
!>   differences to approximate aux_dv.<br>
!> \param[in] itemp1 When debug_aux_dv is .true. partial_ions(itemp1) is
!>   the index of the first dv index to be incremented by deltav1 to
!>   help debug aux_dv using numerical differences.<br>
!> \param[in] itemp2 When debug_aux_dv is .true. partial_ions(itemp2) is
!>   the index of the second dv index to be incremented by deltav2 to
!>   help debug aux_dv using numerical differences.<br>
!> \param[in] deltav1 When debug_aux_dv is .true. partial_ions(itemp1) is
!>   the index of the first dv index to be incremented by deltav1 to
!>   help debug aux_dv using numerical differences.<br>
!> \param[in] deltav2 When debug_aux_dv is .true. partial_ions(itemp2) is
!>   the index of the second dv index to be incremented by deltav2 to
!>   help debug aux_dv using numerical differences.<br>
subroutine eos_warm_step(&
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
     ifsame_abundances, ifmtrace, iatomic_number, ifpi_local,&
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

  use mod_free_eos_constants, only: avogadro, boltzmann, c2, cpe, cr

  use mod_eos_calc, only: eos_calc
  use mod_coulomb, only: master_coulomb, master_coulomb_pressure
  use mod_exchange, only: exchange_pressure
  use mod_pi, only:&
       mdh_pi, mdh_pi_pressure_free,&
       fjs_pi, fjs_pi_free,&
       pteh_pi, pteh_pi_pressure
  use mod_flow_data, only: ln_underflow_limit, ln_overflow_limit
  use mod_info_data, only: info_offset_eos_warm_step
  use mod_excitation, only: excitation_pi, excitation_pi_pressure_free
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  logical, intent(in) :: ifsame_under, ifdvzero
  integer, intent(in) :: nion(:), ion_end(:)
  integer, intent(in) :: verbosity, ifreducedmass, ifmodified, ifexcited, nmax
  integer, intent(in) :: nmin(:), nmin_max(:), izlo, izhi,&
       nmin_species(:), iatomic_number(:), ifnr,&
       ifh2, ifh2plus, ifmtrace,&
       if_pteh, ifcoulomb_mod, if_mc, if_dc,&
       ifionized, max_index,&
       ifpi_local, ifpl,&
       ifsame_abundances, ifsame_zero_abundances,&
       ifdv(:), ifelement(:),&
       partial_ions(:), inv_ion(:),&
       partial_elements(:),&
       inv_aux(:)
  integer, intent(in) :: partial_aux(:)
  integer, intent(in) :: mion_end, n_partial_ions, nextrasum, nxextrasum

  integer, intent(out) :: info

  logical, intent(in) :: old_allow_log(:), old_allow_log_neg(:)

  real(fp_kind), intent(in) :: bi(:), h2diss, t, tl, tc2, eps(:), eta
  real(fp_kind), intent(in) :: rhostar(:), pstar(:)
  real(fp_kind), intent(in) ::&
       dve_exchange, dve_exchangef, dve_exchanget,&
       plop(:), plopt(:), plopt2(:),&
       dv_pl(:), dv_plt(:)
  real(fp_kind), intent(in) ::&
       full_sum0, full_sum1, full_sum2,&
       r_ion3(:), r_neutral(:),&
       charge(:), charge2(:)
  real(fp_kind), intent(in) :: f, wf, n_e, bmin(:)
  real(fp_kind), intent(in) :: aux_old(:)
  real(fp_kind), intent(in) ::&
       sum0_mc, sum2_mc, hcon_mc, hecon_mc,&
       nux, nuy, nuz

  real(fp_kind), intent(out) ::&
       ne, nef, net,&
       dv(:), dvf(:), dvt(:), dv_aux(:,:),&
       dvzero(:),&
       sion, sionf, siont, uion,&
       sumpl0, sumpl0f, sumpl0_dv(:),&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       h2, h2f, h2t, h2_dv(:),&
       pnorad, pnoradf, pnoradt, pnorad_aux(:),&
       fion, fionf, fion_dv(:),&
       ! fion_aux(:),&
       lambda, gamma_e

  real(fp_kind), intent(out) :: aux(:), auxf(:), auxt(:), aux_dv(:, :)

  ! Arguments for debugging
  logical, intent(in), optional :: debug_aux_dv
  integer, intent(in), optional :: itemp1, itemp2
  real(fp_kind), intent(in), optional :: deltav1, deltav2

  ! Internal variables
  logical ifnuform
  logical ifnr03, ifnr13
  logical ifcsum, ifextrasum, ifxextrasum
  logical debug_present

  integer, parameter :: maxfjs_aux = 4
  ! iextraoff is the number of non-extrasum auxiliary variables
  integer, parameter :: iextraoff = 7
  integer ion, index_ion, index, ion0, ielement, index_aux,&
       inv_ion_index, ifjs_aux, iaux, jndex_aux
  ! variables associated with pressure and free energy calculation:
  ! Note, this value must be greater than or equal to 4 as well to store
  ! fjs_pi_free and corresponding pressure quantities.
  integer, parameter ::  maxnextrasum_local = 9
  integer nionsp2, nions, naux, nelements, n_partial_elements, n_partial_aux
  ! maximum core charge for non-bare ion
  integer, parameter :: maxcorecharge = 28
  integer, parameter :: nstar = 9

  real(fp_kind)&
       h_ion, h_ionf, h_iont,&
       he_ion, he_ionf, he_iont,&
       he_ion2, he_ion2f, he_ion2t,&
       rl, rt, rf,&
       h2plus, h2plusf, h2plust,&
       sumion0, sumion2,&
       sum0, sum0f, sum0t,&
       sum2, sum2f, sum2t

  real(fp_kind), allocatable ::&
       h_ion_dv(:),&
       he_ion_dv(:),&
       he_ion2_dv(:),&
       r_dv(:),&
       h2plus_dv(:),&
       sum0_dv(:),&
       sum2_dv(:),&
       extrasum(:),&
       extrasumf(:),&
       extrasumt(:),&
       extrasum_dv(:,:),&
       xextrasum(:),&
       xextrasumf(:),&
       xextrasumt(:),&
       xextrasum_dv(:,:)

  real(fp_kind), allocatable ::&
       dve_pi_aux(:),&
       sum0_aux(:),&
       sum2_aux(:),&
       dvh_aux(:),&
       dvhe1_aux(:),&
       dvhe2_aux(:),&
       ! free_pi_dv(:),&
       ! free_pl_dv(:),&
       ! free_dv(:),&
       pexcited_dv(:),&
       pnorad_dv(:),&
       pion_dv(:),&
       ppi_dv(:),&
       free_excited_dv(:),&
       ni_dv(:),&
       ! free_aux(:),&
       ppi2_aux(:),&
       ! free_pi_aux(:),&
       ppi_aux(:),&
       free_pi2_aux(:)

  real(fp_kind)&
       dve, dvef, dvet, dve0, dve2,&
       dve_pi, dve_pif, dve_pit,&
       dve_coulomb, dve_coulombf, dve_coulombt,&
       dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22, sum0ne, sum2ne,&
       dsum0, dsum2,&
       dvh, dvhf, dvht,&
       dvhe1, dvhe1f, dvhe1t,&
       dvhe2, dvhe2f, dvhe2t,&
       nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,&
       rho, rho_new
  real(fp_kind) pe_cgs,&
       p0, pion, pionf, piont,&
       pex, pexf, pext,&
       ppi, ppif, ppit, ppir,&
       pexcited, pexcitedf, pexcitedt,&
       dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2,&
       ! free_rad, free_e, free_ef,&
       free_pi, free_pif,&
       ! free_pit, free_pir,&
       free_excited, free_excitedf,&
       ! free_coulomb, free_coulombf, free_coulomb0, free_coulomb2,&
       ! free_ex, free_exf,&
       ni, nif!,&
       ! free_pl, free_plf,&
       ! free, freef

  nionsp2 = size(nion)
  nions = nionsp2 - 2
  naux = size(partial_aux)
  nelements = size(ion_end) - 2
  n_partial_elements = size(partial_elements) - 2
  n_partial_aux = size(pnorad_aux) - 2

  ! sanity checks.
  if(nextrasum.gt.maxnextrasum_local) error stop 'eos_warm_step: nextrasum.gt.maxnextrasum_local'
  if(.not.(ifnr.eq.0.or.ifnr.eq.1.or.ifnr.eq.3)) error stop 'eos_warm_step: ifnr must be 0, 1, or 3'
  ! Maintenance, 2021.  ifnr.eq.0 was historically implemented, but it
  ! is currently not used (see sanity check just below), and it has
  ! not been maintained.  So we have commented out the ifnr.eq.0
  ! stanzas below.  Note, that one of the problems with those stanzas
  ! (which would have to be fixed if ifnr.eq.0 was enabled again) is
  ! communication issues with rf, rt, h2plusf, h_ionf, he_ionf, and
  ! he_ion2f leading to uninitialized values of these variables.
  if(ifnr.eq.0) error stop 'eos_warm_step: ifnr = 0 is disabled'
  if(n_partial_aux.gt.naux) error stop 'eos_warm_step: n_partial_aux too large'
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
       nionsp2.ne.size(dv_aux,1).or.&
       nionsp2.ne.size(aux_dv,1).or.&
       nionsp2.ne.size(sumpl0_dv).or.&
       nionsp2.ne.size(fion_dv))&
       error stop 'eos_warm_step: inconsistent nionsp2 array sizes'
  if(&
       naux.ne.size(inv_aux).or.&
       naux.ne.size(old_allow_log).or.&
       naux.ne.size(old_allow_log_neg).or.&
       naux.ne.size(aux_old).or.&
       naux.ne.size(aux).or.&
       naux.ne.size(auxf).or.&
       naux.ne.size(auxt).or.&
       naux.ne.size(aux_dv,2).or.&
       naux.ne.size(dv_aux,2))&
       error stop 'eos_warm_step: inconsistent naux array sizes'
  if(&
       nelements.ne.size(iatomic_number).or.&
       nelements.ne.size(ifelement).or.&
       nelements.ne.size(eps).or.&
       nelements.ne.size(r_neutral)-2.or.&
       nelements.ne.size(dvzero))&
       error stop 'eos_warm_step: inconsistent nelements array sizes'
  if(&
       maxcorecharge.ne.size(nmin).or.&
       maxcorecharge.ne.size(nmin_max).or.&
       maxcorecharge.ne.size(bmin))&
       error stop 'eos_warm_step: inconsistent maxcorecharge array sizes'
  if(&
       nstar.ne.size(rhostar).or.&
       nstar.ne.size(pstar))&
       error stop 'eos_warm_step: inconsistent nstar array sizes'
  !if(&
  !     npartial_aux.ne.size(fion_aux)-2.or.&
  !     npartial_aux.ne.size(free_aux)-2)&
  !     error stop 'eos_warm_step: inconsistent npartial_aux array sizes'

  allocate(&
       ! free_pi_dv(nionsp2),&
       ! free_pi_aux(naux),&
       ! free_pl_dv(nionsp2),&
       ! free_dv(nionsp2),&
       ! free_aux(n_partial_aux+2),&
       dve_pi_aux(maxfjs_aux),&
       sum0_aux(maxfjs_aux+1),&
       sum2_aux(maxfjs_aux+1),&
       dvh_aux(maxfjs_aux),&
       dvhe1_aux(maxfjs_aux),&
       dvhe2_aux(maxfjs_aux),&
       pion_dv(nionsp2),&
       ppi_dv(nionsp2),&
       ppi2_aux(maxnextrasum_local),&
       ppi_aux(naux),&
       pexcited_dv(nionsp2),&
       pnorad_dv(nionsp2),&
       free_pi2_aux(maxnextrasum_local),&
       free_excited_dv(nionsp2),&
       ni_dv(nionsp2),&
       h_ion_dv(nionsp2),&
       he_ion_dv(nionsp2),&
       he_ion2_dv(nionsp2),&
       r_dv(nionsp2),&
       h2plus_dv(nionsp2),&
       sum0_dv(nionsp2),&
       sum2_dv(nionsp2),&
       extrasum(nextrasum),&
       extrasumf(nextrasum),&
       extrasumt(nextrasum),&
       extrasum_dv(nionsp2,nextrasum),&
       xextrasum(nxextrasum),&
       xextrasumf(nxextrasum),&
       xextrasumt(nxextrasum),&
       xextrasum_dv(nionsp2, nxextrasum)&
       )

  if(if_taint_allocated_real) then
     ! free_pi_dv = tainted_fp_kind
     ! free_pi_aux = tainted_fp_kind
     ! free_pl_dv = tainted_fp_kind
     ! free_dv = tainted_fp_kind
     ! free_aux = tainted_fp_kind
     call taint_allocated_real(dve_pi_aux)
     call taint_allocated_real(sum0_aux)
     call taint_allocated_real(sum2_aux)
     call taint_allocated_real(dvh_aux)
     call taint_allocated_real(dvhe1_aux)
     call taint_allocated_real(dvhe2_aux)
     call taint_allocated_real(pion_dv)
     call taint_allocated_real(ppi_dv)
     call taint_allocated_real(ppi2_aux)
     call taint_allocated_real(ppi_aux)
     call taint_allocated_real(pexcited_dv)
     call taint_allocated_real(pnorad_dv)
     call taint_allocated_real(free_pi2_aux)
     call taint_allocated_real(free_excited_dv)
     call taint_allocated_real(ni_dv)
     call taint_allocated_real(h_ion_dv)
     call taint_allocated_real(he_ion_dv)
     call taint_allocated_real(he_ion2_dv)
     call taint_allocated_real(r_dv)
     call taint_allocated_real(h2plus_dv)
     call taint_allocated_real(sum0_dv)
     call taint_allocated_real(sum2_dv)
     call taint_allocated_real(extrasum)
     call taint_allocated_real(extrasumf)
     call taint_allocated_real(extrasumt)
     call taint_allocated_real(extrasum_dv)
     call taint_allocated_real(xextrasum)
     call taint_allocated_real(xextrasumf)
     call taint_allocated_real(xextrasumt)
     call taint_allocated_real(xextrasum_dv)
  endif

  ifcsum = .not.(if_mc.eq.1.or.if_pteh.eq.1)
  ifextrasum = ifpi_local.eq.3.or.ifpi_local.eq.4
  ifxextrasum = ifextrasum .and. (ifexcited.gt.0)

  call aux_to_traditional(&
       iextraoff, iextraoff + maxnextrasum_local, inv_aux, aux_old,&
       h_ion, he_ion, he_ion2, rl, h2plus, sum0, sum2, extrasum, xextrasum)

  ifnr03 = ifnr.eq.0.or.ifnr.eq.3
  ifnr13 = ifnr.eq.1.or.ifnr.eq.3
  if(if_mc.eq.1.or.ifpi_local.eq.2) then
     ! Sanity check.
     ! rl only defined for this case.
     if(inv_aux(4).eq.0)&
          error stop "eos_warm_step: inv_aux(4), if_mc, and ifpi_local are not consistent"
     if(rl.lt.ln_underflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_warm_step: calculated mass density too close to underflow limit"
        info = info_offset_eos_warm_step + 1
        return
     elseif(rl.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_warm_step: calculated mass density too close to overflow limit"
        info = info_offset_eos_warm_step + 2
        return
     endif
     ! Calculate rho for this case where rl is available and in the allowed range
     ! and when rho is needed below (for if_mc.eq.1.or.ifpi_local.eq.2)
     rho = exp(rl)
  else
     ! gfortran generates a spurious [-Wmaybe-uninitialized] warning
     ! unless do this unneeded (since rho only used for the
     ! if_mc.eq.1.or.ifpi_local.eq.2 case below) initialization of rho here.
     rho = 0._fp_kind
  endif

  do index = 1, n_partial_elements
     ielement = partial_elements(index)
     dvzero(ielement) = 0._fp_kind
  enddo
  do index_ion = 1,max_index
     ion = partial_ions(index_ion)
     dv(ion) = 0._fp_kind
     if(ifnr03) then
        dvf(ion) = 0._fp_kind
        dvt(ion) = 0._fp_kind
     endif
  enddo
  if(ifnr13) then
     dv_aux(1:max_index, 1:n_partial_aux) = 0._fp_kind
  endif
  ppi_aux(1:n_partial_aux) = 0._fp_kind
  if(ifpi_local.eq.0) then
     dve_pi = 0._fp_kind
     if(ifnr03) then
        dve_pif = 0._fp_kind
        dve_pit = 0._fp_kind
     endif
  elseif(ifpi_local.eq.1) then
     call pteh_pi(ifmodified,&
          rhostar(1), rhostar(2), rhostar(3), t,&
          dve_pi, dve_pif, dve_pit)
  elseif(ifpi_local.eq.2) then
     ! use FJS pressure ionization
     nx = nux*rho*avogadro
     ny = nuy*rho*avogadro
     nz = nuz*rho*avogadro
     !if(ifnr.eq.0) then
        !nxf = nx*rf
        !nxt = nx*rt
        !nyf = ny*rf
        !nyt = ny*rt
        !nzf = nz*rf
        !nzt = nz*rt
     !elseif(ifnr.eq.3) then
     if(ifnr.eq.3) then
        ! For ifnr.eq.3 calculate all f and t derivatives assuming input
        ! auxiliary variables (e.g., rho for ifpi_local.eq.2) are fixed.
        nxf = 0._fp_kind
        nxt = 0._fp_kind
        nyf = 0._fp_kind
        nyt = 0._fp_kind
        nzf = 0._fp_kind
        nzt = 0._fp_kind
     endif
     call fjs_pi(ifnr, ifmodified,&
          t, nx, nxf, nxt, ny, nyf, nyt, nz, nzf, nzt,&
          h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t,&
          n_e, n_e*rhostar(2), n_e*rhostar(3),&
          dvh, dvhf, dvht, dvh_aux,&
          dvhe1, dvhe1f, dvhe1t, dvhe1_aux,&
          dvhe2, dvhe2f, dvhe2t, dvhe2_aux,&
          dve_pi, dve_pif, dve_pit, dve_pi_aux)
     ! this call must occur before call to eos_calc (which
     ! calculates new values of rho (used to calculate nx, ny, and nz),
     ! h_ion, he_ion, and he_ion2.
     ! n.b. ppi and all its derivatives are equivalent to
     ! free_pi and all its derivatives for the fjs form of
     ! pressure ionization.
     call fjs_pi_free(inv_aux, nx, ny, nz, h_ion, he_ion, he_ion2,&
          t, n_e, n_e*rhostar(2), n_e*rhostar(3),&
          ppi, ppif, ppit, ppi_aux)
     ! update dv keeping in mind that calculated dvh etc. have
     ! dve_pi included.  the dve_pi quantity
     ! must be subtracted so that further free_eos_detailed logic
     ! which adds dve_pi works properly.
     if(eps(1).gt.0._fp_kind) then
        dv(1) = dv(1) + dvh - dve_pi
        if(ifnr03) then
           dvf(1) = dvf(1) + dvhf - dve_pif
           dvt(1) = dvt(1) + dvht - dve_pit
        endif
        if(ifnr13) then
           inv_ion_index = inv_ion(1)
           do ifjs_aux = 1,maxfjs_aux
              iaux = inv_aux(ifjs_aux)
              if(iaux.gt.0) then
                 dv_aux(inv_ion_index,iaux) = dv_aux(inv_ion_index,iaux) +&
                      dvh_aux(ifjs_aux) - dve_pi_aux(ifjs_aux)
              endif
           enddo
        endif
     endif
     if(eps(2).gt.0._fp_kind) then
        if(ifdvzero) then
           dv(2) = dv(2) + dvhe1 - dve_pi
           dv(3) = dv(3) + dvhe1 + dvhe2 - 2._fp_kind*dve_pi
        else
           dvzero(2) = dvzero(2) + dvhe1 + dvhe2 - 2._fp_kind*dve_pi
           dv(2) = dv(2) - dvhe2 + dve_pi
        endif
        if(ifnr03) then
           dvf(2) = dvf(2) + dvhe1f - dve_pif
           dvt(2) = dvt(2) + dvhe1t - dve_pit
           dvf(3) = dvf(3) + dvhe1f + dvhe2f - 2._fp_kind*dve_pif
           dvt(3) = dvt(3) + dvhe1t + dvhe2t - 2._fp_kind*dve_pit
        endif
        if(ifnr13) then
           inv_ion_index = inv_ion(2)
           do ifjs_aux = 1,maxfjs_aux
              iaux = inv_aux(ifjs_aux)
              if(iaux.gt.0) then
                 dv_aux(inv_ion_index,iaux) = dv_aux(inv_ion_index,iaux) +&
                      dvhe1_aux(ifjs_aux) - dve_pi_aux(ifjs_aux)
              endif
           enddo
           inv_ion_index = inv_ion(3)
           do ifjs_aux = 1,maxfjs_aux
              iaux = inv_aux(ifjs_aux)
              if(iaux.gt.0) then
                 dv_aux(inv_ion_index,iaux) = dv_aux(inv_ion_index,iaux) +&
                      dvhe1_aux(ifjs_aux) + dvhe2_aux(ifjs_aux) - 2._fp_kind*dve_pi_aux(ifjs_aux)
              endif
           enddo
        endif
     endif
     ! no change in H2 equilibrium constant because fjs_pi is independent
     ! n(H) or n(H2).  Change in H2+ equilibrium constant
     ! via dve_pi which is handled later as part of dve.
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
     ! this call must occur before call to eos_calc (which
     ! updates extrasum to new values.)
     call mdh_pi_pressure_free(t, extrasum,&
          ppi, ppif, ppit, ppi2_aux,&
          free_pi, free_pif, free_pi2_aux)
     ! MDH pressure ionization has no N_e dependence
     dve_pi = 0._fp_kind
     if(ifnr03) then
        dve_pif = 0._fp_kind
        dve_pit = 0._fp_kind
     endif
  endif
  ! calculate change in equilibrium constant due to electron degeneracy
  dve = dve_pi - eta

  ! gfortran generates a spurious [-Wmaybe-uninitialized] warning
  ! unless do this redundant initialization of dvef and dvet here.
  dvef = 0._fp_kind
  dvet = 0._fp_kind

  if(ifnr03) then
     dvef = dve_pif - wf
     dvet = dve_pit
  endif
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

  ! sumion0 and sumion2 only calculated and used below for the if_mc.eq.1 case.
  ! However, gfortran generates a spurious [-Wmaybe-uninitialized] warning
  ! unless do this redundant initialization of sumion0 and sumion2 here.
  sumion0 = 0._fp_kind
  sumion2 = 0._fp_kind
  if(if_mc.eq.1) then
     ! N.B. this code depends on auxiliary variables that are only
     ! defined for ifpi_local = 2 case.  So it only works because
     ! if_mc.eq.1 is correlated with that case.  Check this.
     if(ifpi_local.ne.2) error stop 'eos_warm_step: logical screwup wrt ifpi_local'

     ! For the if_mc.eq.1 case, sum0_mc and sum2_mc should contain the
     ! component of the nu form of sum0 and sum2 metals that is
     ! approximated as fully ionized and which is therefore not
     ! accounted for by hcon_mc and hecon_mc.
     sumion0 = sum0_mc*rho*avogadro
     sumion2 = sum2_mc*rho*avogadro
     if(ifh2plus.gt.0) then
        sum0 = sumion0 + h2plus + h_ion*(1._fp_kind+hcon_mc) + he_ion + he_ion2
        sum2 = sumion2 + h2plus + h_ion*(1._fp_kind+hcon_mc) + he_ion + he_ion2*(4._fp_kind+hecon_mc)
        !if(ifnr.eq.0) then
           !sum0f = sumion0*rf + h2plusf + h_ionf*(1.d0+hcon_mc) + he_ionf + he_ion2f
           !sum0t = sumion0*rt + h2plust + h_iont*(1.d0+hcon_mc) + he_iont + he_ion2t
           !sum2f = sumion2*rf + h2plusf + h_ionf*(1.d0+hcon_mc) + he_ionf + he_ion2f*(4.d0+hecon_mc)
           !sum2t = sumion2*rt + h2plust + h_iont*(1.d0+hcon_mc) + he_iont + he_ion2t*(4.d0+hecon_mc)
        !endif
     else
        ! Exclude uninitialized h2plus and derivatives from the calculation.
        sum0 = sumion0 + h_ion*(1._fp_kind+hcon_mc) + he_ion + he_ion2
        sum2 = sumion2 + h_ion*(1._fp_kind+hcon_mc) + he_ion + he_ion2*(4._fp_kind+hecon_mc)
        !if(ifnr.eq.0) then
           !sum0f = sumion0*rf + h_ionf*(1.d0+hcon_mc) + he_ionf + he_ion2f
           !sum0t = sumion0*rt + h_iont*(1.d0+hcon_mc) + he_iont + he_ion2t
           !sum2f = sumion2*rf + h_ionf*(1.d0+hcon_mc) + he_ionf + he_ion2f*(4.d0+hecon_mc)
           !sum2t = sumion2*rt + h_iont*(1.d0+hcon_mc) + he_iont + he_ion2t*(4.d0+hecon_mc)
        !endif
     endif
     if(ifnr.eq.3) then
        ! For ifnr.eq.3 calculate all f and t derivatives assuming input
        ! auxiliary variables are fixed.
        sum0f = 0._fp_kind
        sum0t = 0._fp_kind
        sum2f = 0._fp_kind
        sum2t = 0._fp_kind
     endif
     if(ifnr13) then
        sum0_aux(1) = (1._fp_kind+hcon_mc)
        sum2_aux(1) = (1._fp_kind+hcon_mc)
        sum0_aux(2) = 1._fp_kind
        sum2_aux(2) = 1._fp_kind
        sum0_aux(3) = 1._fp_kind
        sum2_aux(3) = (4._fp_kind+hecon_mc)
        sum0_aux(4) = sumion0
        sum2_aux(4) = sumion2
        sum0_aux(5) = 1._fp_kind
        sum2_aux(5) = 1._fp_kind
     endif
  endif
  call master_coulomb(ifnr,&
       rhostar, f,&
       sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,&
       n_e, t, cpe*pstar(1), pstar, lambda, gamma_e,&
       ifcoulomb_mod, if_dc, if_pteh,&
       dve_coulomb, dve_coulombf, dve_coulombt, dve0, dve2,&
       dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
  ! N.B. this routine must be called _before_ sum0 and sum2 are changed
  ! in value.
  call master_coulomb_pressure(&
       ifcoulomb_mod, if_dc, if_pteh, sum0, sum2,&
       t, n_e, rhostar, pstar,&
       dpcoulomb, dpcoulombf, dpcoulombt, dpcoulomb0, dpcoulomb2)
  dve = dve + dve_coulomb
  if(ifnr03) then
     dvef = dvef + dve_coulombf
     dvet = dvet + dve_coulombt
  endif
  ! add in pre-calculated exchange effects
  dve = dve + dve_exchange
  if(ifnr03) then
     dvef = dvef + dve_exchangef
     dvet = dvet + dve_exchanget
  endif
  if(if_mc.ne.1) then
     index_aux = inv_aux(iextraoff-1)
     jndex_aux = inv_aux(iextraoff)
     do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if(.not.(ifdvzero.or.ielement.eq.1)) then
           dvzero(ielement) = dvzero(ielement) +&
                dv0 + charge(ion_end(ielement))*(dve + charge(ion_end(ielement))*dv2)
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
        if(ifnr03) then
           dvf(ion0+1:ion_end(ielement)) = dvf(ion0+1:ion_end(ielement)) +&
                dv0f + charge(ion0+1:ion_end(ielement))*(dvef + charge(ion0+1:ion_end(ielement))*dv2f)
           dvt(ion0+1:ion_end(ielement)) = dvt(ion0+1:ion_end(ielement)) +&
                dv0t + charge(ion0+1:ion_end(ielement))*(dvet + charge(ion0+1:ion_end(ielement))*dv2t)
        endif
        do ion = ion0+1,ion_end(ielement)
           if(ifnr13.and.index_aux.ne.0) then
              index_ion = inv_ion(ion)
              dv_aux(index_ion,index_aux) = dv_aux(index_ion,index_aux) +&
                   dv00 + charge(ion)*(dve0 + charge(ion)*dv02)
              dv_aux(index_ion,jndex_aux) = dv_aux(index_ion,jndex_aux) +&
                   dv02 + charge(ion)*(dve2 + charge(ion)*dv22)
           endif
        enddo
     enddo
     ! H2+
     if(ifdv(nions+2).eq.1) then
        dv(nions+2) = dv(nions+2) + dv0 + dve + dv2
        if(ifnr03) then
           dvf(nions+2) = dvf(nions+2) + dv0f + dvef + dv2f
           dvt(nions+2) = dvt(nions+2) + dv0t + dvet + dv2t
        endif
        if(ifnr13.and.index_aux.ne.0) then
           ! max_index is correct pointer when ifdv(nions+2).eq.1
           dv_aux(max_index,index_aux) = dv_aux(max_index,index_aux) + dv00 + dve0 + dv02
           dv_aux(max_index,jndex_aux) = dv_aux(max_index,jndex_aux) + dv02 + dve2 + dv22
        endif
     endif
  else
     ! if_mc *is* 1
     do index = 1, n_partial_elements
        ielement = partial_elements(index)
        if(ifdvzero.or.ielement.eq.1) then
        elseif(ielement.eq.2) then
           dvzero(ielement) = dvzero(ielement) +&
                dv0 + charge(ion_end(ielement))*(dve + charge(ion_end(ielement))*dv2)
        else
           dvzero(ielement) = dvzero(ielement) +&
                charge(ion_end(ielement))*dve
        endif
        if(ielement.gt.1) then
           ion0 = ion_end(ielement-1)
        else
           ion0 = 0
        endif

        ! Sanity check:
        if(ion0.lt.0) error stop 'eos_warm_step: invalid ion0'

        do ion = ion0+1,ion_end(ielement)
           if(ifnr13) then
              ! From above sanity check, ion.ge.1 ==> dsum0 and dsum2 always defined by
              ! the various ion Boolean code blocks below.
              ! However, gfortran generates a spurious [-Wmaybe-uninitialized] warning
              ! unless we redundantly initialize dsum0 and dsum2 here.
              dsum0 = 0._fp_kind
              dsum2 = 0._fp_kind

              if(ion.ge.4) then
                 dsum0 = charge(ion)*dve0
                 dsum2 = charge(ion)*dve2
              elseif(ion.eq.1) then
                 dsum0 = charge(ion)*dve0 + (1._fp_kind + hcon_mc)*(dv00 + charge2(ion)*dv02)
                 dsum2 = charge(ion)*dve2 + (1._fp_kind + hcon_mc)*(dv02 + charge2(ion)*dv22)
              elseif(ion.eq.2) then
                 dsum0 = dv00 + charge(ion)*(dve0 + charge(ion)*dv02)
                 dsum2 = dv02 + charge(ion)*(dve2 + charge(ion)*dv22)
              elseif(ion.eq.3) then
                 dsum0 = charge(ion)*dve0 + dv00 + (hecon_mc + charge2(ion))*dv02
                 dsum2 = charge(ion)*dve2 + dv02 + (hecon_mc + charge2(ion))*dv22
              endif
              index_ion = inv_ion(ion)
              do ifjs_aux = 1, maxfjs_aux+1
                 index_aux = inv_aux(ifjs_aux)
                 if(index_aux.gt.0)&
                      dv_aux(index_ion,index_aux) = dv_aux(index_ion,index_aux) +&
                      dsum0*sum0_aux(ifjs_aux) + dsum2*sum2_aux(ifjs_aux)
              enddo
           endif
           if(ion.ge.4) then
              if(ifdvzero) then
                 dv(ion) = dv(ion) + charge(ion)*dve
              else
                 dv(ion) = dv(ion) +&
                      real((nion(ion)-nion(ion_end(ielement))),fp_kind)*dve
              endif
              if(ifnr03) then
                 dvf(ion) = dvf(ion) + charge(ion)*dvef
                 dvt(ion) = dvt(ion) + charge(ion)*dvet
              endif
           elseif(ion.eq.1) then
              dv(ion) = dv(ion) +&
                   charge(ion)*dve + (1._fp_kind + hcon_mc)*(dv0 + charge2(ion)*dv2)
              if(ifnr03) then
                 dvf(ion) = dvf(ion) + charge(ion)*dvef + (1._fp_kind + hcon_mc)*(dv0f + charge2(ion)*dv2f)
                 dvt(ion) = dvt(ion) + charge(ion)*dvet + (1._fp_kind + hcon_mc)*(dv0t + charge2(ion)*dv2t)
              endif
           elseif(ion.eq.2) then
              if(ifdvzero) then
                 dv(ion) = dv(ion) +&
                      dv0 + charge(ion)*(dve + charge(ion)*dv2)
              else
                 dv(ion) = dv(ion) +&
                      real((nion(ion)-nion(ion_end(ielement))),fp_kind)*&
                      (dve + real((nion(ion)+nion(ion_end(ielement))),fp_kind)*dv2)
              endif
              if(ifnr03) then
                 dvf(ion) = dvf(ion) + dv0f + charge(ion)*(dvef + charge(ion)*dv2f)
                 dvt(ion) = dvt(ion) + dv0t + charge(ion)*(dvet + charge(ion)*dv2t)
              endif
           elseif(ion.eq.3) then
              if(ifdvzero) then
                 dv(ion) = dv(ion) + charge(ion)*dve + dv0 + (hecon_mc + charge2(ion))*dv2
              else
                 dv(ion) = dv(ion) + hecon_mc*dv2
              endif
              if(ifnr03) then
                 dvf(ion) = dvf(ion) + charge(ion)*dvef + dv0f + (hecon_mc + charge2(ion))*dv2f
                 dvt(ion) = dvt(ion) + charge(ion)*dvet + dv0t + (hecon_mc + charge2(ion))*dv2t
              endif
           endif
        enddo
     enddo
     ! H2+
     if(ifdv(nions+2).eq.1) then
        dv(nions+2) = dv(nions+2) + dv0 + dve + dv2
        if(ifnr03) then
           dvf(nions+2) = dvf(nions+2) + dv0f + dvef + dv2f
           dvt(nions+2) = dvt(nions+2) + dv0t + dvet + dv2t
        endif
        if(ifnr13) then
           dsum0 = dv00 + dve0 + dv02
           dsum2 = dv02 + dve2 + dv22
           do ifjs_aux = 1, maxfjs_aux+1
              index_aux = inv_aux(ifjs_aux)
              ! max_index is correct pointer when ifdv(nions+2).eq.1
              if(index_aux.gt.0)&
                   dv_aux(max_index,index_aux) = dv_aux(max_index,index_aux) +&
                   dsum0*sum0_aux(ifjs_aux) + dsum2*sum2_aux(ifjs_aux)
           enddo
        endif
     endif
  endif
  if(ifpi_local.eq.2.and.ifnr13) then
     ! in all monatomic cases added term of charge(ion)*dve to dv
     ! for ifpi_local == 2, there is a pressure ionization
     ! component to dve.
     do index_ion = 1,n_partial_ions
        ion = partial_ions(index_ion)
        do ifjs_aux = 1,maxfjs_aux
           index_aux = inv_aux(ifjs_aux)
           if(index_aux.gt.0) then
              dv_aux(index_ion,index_aux) = dv_aux(index_ion,index_aux) +&
                   charge(ion)*dve_pi_aux(ifjs_aux)
           endif
        enddo
     enddo
     ! for H2+ added term of dve to dv(nions+2)
     if(ifdv(nions+2).eq.1) then
        do ifjs_aux = 1,maxfjs_aux
           index_aux = inv_aux(ifjs_aux)
           if(index_aux.gt.0) then
              ! max_index is correct pointer when ifdv(nions+2).eq.1
              dv_aux(max_index,index_aux) = dv_aux(max_index,index_aux) +&
                   dve_pi_aux(ifjs_aux)
           endif
        enddo
     endif
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

  debug_present = present(debug_aux_dv).and.present(itemp1).and.present(itemp2).and.present(deltav1).and.present(deltav2)

  if(debug_present) then
     if(debug_aux_dv) then
        ! specify change in dv as a function of fl and tl
        dv(partial_ions(itemp1)) = dv(partial_ions(itemp1)) + deltav1
        dv(partial_ions(itemp2)) = dv(partial_ions(itemp2)) + deltav2
     endif
  endif
  ! add in Planck-Larkin occupation probability effect
  if(ifpl.eq.1) then
     do index_ion = 1,max_index
        ion = partial_ions(index_ion)
        dv(ion) = dv(ion) + dv_pl(ion)
        if(ifnr03) then
           dvt(ion) = dvt(ion) + dv_plt(ion)
        endif
     enddo
  endif
  ! eos_warm_step always returns n form rather than nu = n/(rho*avagadro)
  ! form for h2, h2plus, h_ion, he_ion, he_ion2, and extrasum.
  ifnuform = .false.
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
  rho_new = exp(rl)

  ! Calculate pressure and its aux and fl derivatives.
  ! free electrons component
  pe_cgs = cpe*pstar(1)
  ! free ions component
  p0 = cr*rho_new*t
  pion = full_sum0*p0 - (h2+h2plus)*boltzmann*t
  if(ifnr03) then
     pionf = full_sum0*p0*rf - (h2f+h2plusf)*boltzmann*t
     piont = full_sum0*p0*(1._fp_kind+rt) - (h2+h2plus+h2t+h2plust)*boltzmann*t
  endif
  if(ifnr13) then
     pion_dv(1:max_index) = full_sum0*p0*r_dv(1:max_index) - (h2_dv(1:max_index) + h2plus_dv(1:max_index))*boltzmann*t
  endif
  ! exchange component
  call exchange_pressure(rhostar, pstar, pex, pexf, pext)
  ! pressure ionization component
  ppi_dv(1:max_index) = 0._fp_kind
  if(ifpi_local.eq.0) then
     ! there is no free energy term from pressure ionization.
     ppi = 0._fp_kind
     ppif = 0._fp_kind
     ppit = 0._fp_kind
  elseif(ifpi_local.eq.1) then
     call pteh_pi_pressure(&
          ifnr, full_sum1, rho_new, rf, rt, t, ne, nef, net,&
          ppi, ppif, ppit, ppir)
     ppi_dv(1:max_index) = ppir*r_dv(1:max_index)
  elseif(ifpi_local.eq.2) then
     ! ppi_aux already updated inside fjs_pi_free that was previously called.
  elseif(ifpi_local.eq.3.or.ifpi_local.eq.4) then
     ! ppi, ppif, ppit, and ppi2_aux defined as a result of
     ! the call to mdh_pi_pressure_free done earlier in this routine.
     do index = 1,nextrasum
        index_aux = inv_aux(index+iextraoff)
        if(index_aux.gt.0)&
             ppi_aux(index_aux) = ppi2_aux(index)
     enddo
  endif
  if(ifexcited.gt.0) then
     call excitation_pi_pressure_free(&
          t, rho_new, rf, rt, r_dv(1:max_index),&
          pexcited, pexcitedf, pexcitedt, pexcited_dv(1:max_index),&
          free_excited, free_excitedf, free_excited_dv(1:max_index))
  else
     pexcited = 0._fp_kind
     pexcitedf = 0._fp_kind
     pexcitedt = 0._fp_kind
     pexcited_dv(1:max_index) = 0._fp_kind
     free_excited = 0._fp_kind
     free_excitedf = 0._fp_kind
     ! free_excited_dv(1:max_index) = 0.d0
  endif
  ! N.B. there are no pressure terms for Planck-Larkin occupation
  ! probability.
  pnorad = pe_cgs + pion + pex + ppi + pexcited
  pnorad = pnorad + dpcoulomb
  if(ifnr03) then
     pnoradf = pe_cgs*pstar(2) + pionf + pexf + ppif + pexcitedf
     pnoradf = pnoradf + dpcoulombf
     pnoradt = pe_cgs*pstar(3) + piont + pext + ppit + pexcitedt
     pnoradt = pnoradt + dpcoulombt
  endif
  if(ifnr13) then
     pnorad_dv(1:max_index) = pion_dv(1:max_index) + ppi_dv(1:max_index) + pexcited_dv(1:max_index)
     do index_aux = 1,n_partial_aux
        pnorad_aux(index_aux) = ppi_aux(index_aux) + dot_product(pnorad_dv(1:max_index), dv_aux(1:max_index,index_aux))
     enddo
  endif
  if(if_pteh.ne.1) then
     if(if_mc.eq.1) then
        ! sum0 and sum2 calculated according to following formulas:
        ! sum0 = sumion0 + h2plus + h_ion*(1.d0+hcon_mc) + he_ion + he_ion2
        ! sum2 = sumion2 + h2plus + h_ion*(1.d0+hcon_mc) + he_ion + he_ion2*(4.d0+hecon_mc)
        ! where h_ion, he_ion, he_ion2, and h2plus are the first, second, third
        ! and fifth auxiliary variables,
        ! hcon_mc and hecon_mc are constants, and
        ! sumion0 and sumion2 are constants times exp(rl), where rl
        ! is the 4th auxiliary variable.
        ! partial wrt h_ion.
        index_aux = inv_aux(1)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             (1._fp_kind+hcon_mc)*(dpcoulomb0+dpcoulomb2)
        ! partial wrt he_ion.
        index_aux = inv_aux(2)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             (dpcoulomb0+dpcoulomb2)
        ! partial wrt he_ion2.
        index_aux = inv_aux(3)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             dpcoulomb0 + (4._fp_kind+hecon_mc)*dpcoulomb2
        ! partial wrt (old) rl.
        index_aux = inv_aux(4)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             sumion0*dpcoulomb0 + sumion2*dpcoulomb2
        ! partial wrt h2plus.
        index_aux = inv_aux(5)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             (dpcoulomb0+dpcoulomb2)
     else
        ! sum0
        index_aux = inv_aux(6)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             dpcoulomb0
        ! sum2
        index_aux = inv_aux(7)
        if(index_aux.gt.0)&
             pnorad_aux(index_aux) = pnorad_aux(index_aux) +&
             dpcoulomb2
     endif
  endif

  ! Final transformation of fion:
  ! add in (1.d0 + 1.5d0*tl - rl)) term summed over all species,
  ! multiply by RT to obtain free-energy in ergs, and put in
  ! hydrogen species energy zero point shift of h2diss for the
  ! two diatomic species and h2diss/2 for the monatomic species.
  ni = full_sum0 - (h2+h2plus)/(rho_new*avogadro)
  fion = cr*t*(fion - ni*(1._fp_kind + 1.5_fp_kind*tl - rl)) + 0.5_fp_kind*(c2*cr)*h2diss*eps(1)
  if(ifnr03) then
     nif = ((h2+h2plus)*rf - (h2f+h2plusf))/(rho_new*avogadro)
     fionf = cr*t*(fionf - nif*(1._fp_kind + 1.5_fp_kind*tl - rl) + ni*rf)
  endif
  if(ifnr13) then
     ni_dv(1:max_index) = ((h2+h2plus)*r_dv(1:max_index) -&
          (h2_dv(1:max_index)+h2plus_dv(1:max_index)))/(rho_new*avogadro)
     fion_dv(1:max_index) = cr*t*(fion_dv(1:max_index) -&
          ni_dv(1:max_index)*(1._fp_kind + 1.5_fp_kind*tl - rl) + ni*r_dv(1:max_index))
  endif

  do index_aux = 1,n_partial_aux
     ! take derivative with respect to log(auxold) or log(-auxold)
     iaux = partial_aux(index_aux)
     if(old_allow_log(index_aux).or.old_allow_log_neg(index_aux))&
          pnorad_aux(index_aux) = pnorad_aux(index_aux)*aux_old(iaux)
  enddo
  if(ifnr03) then
     pnorad_aux(n_partial_aux+1) = pnoradf
     pnorad_aux(n_partial_aux+2) = pnoradt
  endif
  call traditional_to_aux(&
       if_mc, ifcsum, ifextrasum, ifxextrasum, ifnr03, ifnr13, iextraoff, iextraoff + maxnextrasum_local, max_index,&
       h_ion, h_ionf, h_iont, h_ion_dv,&
       he_ion, he_ionf, he_iont, he_ion_dv,&
       he_ion2, he_ion2f, he_ion2t, he_ion2_dv,&
       rl, rf, rt, r_dv,&
       h2plus, h2plusf, h2plust, h2plus_dv,&
       sum0, sum0f, sum0t, sum0_dv,&
       sum2, sum2f, sum2t, sum2_dv,&
       extrasum, extrasumf, extrasumt, extrasum_dv,&
       xextrasum, xextrasumf, xextrasumt, xextrasum_dv,&
       aux, auxf, auxt, aux_dv)

end subroutine eos_warm_step
