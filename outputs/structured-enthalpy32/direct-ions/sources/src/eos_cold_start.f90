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

!> The purpose of this eos_cold_start subroutine is to calculate the
!> EOS (including all "output" auxiliary variables and their
!> derivatives) for a simplified form of free-energy model that does
!> not depend on "input" auxiliary variables. This routine is called
!> from free_eos_detailed for both the case where the simplified EOS
!> is the final result that is wanted and the case where the "output"
!> auxiliary variables calculated for the simplified EOS are used for
!> a starting approximation for a more realistic EOS.
!>
!> This routine first calculates dv and its partial derivatives wrt fl
!> and tl where dv is the change in the equilibrium constant of an ion
!> relative to the un-ionized reference state and the free electron.
!> dv is defined by<br>
!> dv = (partial F/partial nref - partial F/partial nion - ion*partial
!> F/partial ne)/kT.<br>
!> N.B. with the simplified free-energy model, dv is
!> independent of "input" auxiliary variables.  After the dv, dvf, and
!> dvt values (held in arrays) are determined, eos_calc is called to
!> calculate the ionization/molecular equilibrium and all "output"
!> auxiliary variables and their fl and tl derivatives
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
!> \param[in] ifsame_under
!>   ifsame_under .eqv. .true. means use same underflow zeroing for each
!>     ion as in previous call.  This option is useful for removing
!>     small discontinuities caused by variations in the underflow
!>     zeroing which sometimes foil the last stages of
!>     convergence.<br>
!>   ifsame_under .eqv. .false. means calculate underflow zeroing.<br>
!> \param[out] lambda lambda is the calculated Coulomb interaction parameter.  Note
!>   lambda has two different definitions depending on ifcoulomb.<br>
!> \param[out] gamma_e gamma_e is the calculated Coulomb diffraction parameter.<br>
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
!> \param[in] ifnuform
!>   ifnuform .eqv. .true. means return the nu = n/(rho*NA) form of all number densities and most auxiliary variables.<br>
!>   ifnuform .eqv. .false. means return the n form of all number densities and most auxiliary variables.<br>
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
!>   bi(nionsp2) ionization potentials (cm**-1) in ion order from
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
!>     it in differences of dv quanties we can gain many significant
!>     digits.<br>
!> \param[out] dv
!>   dv(nions_inp2) contains the calculated equilibrium constants (explained above).<br>
!> \param[out] dvf
!>   dvf(nions_inp2) contains the calculated partial derivative of dv wrt fl.<br>
!> \param[out] dvt
!>   dvt(nions_inp2) contains the calculated partial derivative of dv wrt ln t.<br>
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
!> \param[in] nextrasum nextrasum is the size of the extrasum and related
!>   arrays that will be locally allocated as traditional auxiliary
!>   variables used by eos_calc and copied back to the appropriate
!>   parts of aux and related non-traditional auxiliary variables by
!>   traditional_to_aux.  nextrasum is either 0 or 9 depending on the
!>   adopted free-energy model, i.e., depending on whether the Boolean
!>   value of (ifpi_local.eq.3.or.ifpi_local.eq.4) which is assigned
!>   to the local variable ifextrasum is either .true. or .false.<br>
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
!> \param[out] info
!>   info is the calculated return code for the procedure which is
!>     zero for success and non-zero (with offset to identify which
!>     procedure the error occurred in) if an error occurred within
!>     the procedure or within some procedure it called.<br>

subroutine eos_cold_start(&
     verbosity, ifsame_under, lambda, gamma_e,&
     partial_ions, f, eta, wf, t, n_e, pstar,&
     dve_exchange, dve_exchangef, dve_exchanget,&
     full_sum0, full_sum1, full_sum2, charge, dv_pl, dv_plt,&
     ifdv, ifdvzero, ifcoulomb_mod, if_dc,&
     ifsame_zero_abundances, ifexcited, ifnuform,&
     inv_aux,&
     inv_ion, max_index,&
     partial_elements, ion_end,&
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

  use mod_free_eos_constants, only: cpe
  use mod_aux_scale, only: sum0_scale, sum2_scale, extrasum_scale, xextrasum_scale
  use mod_eos_calc, only: eos_calc
  use mod_coulomb, only: master_coulomb
  use mod_pi, only: pteh_pi
  use mod_eos_jacobian, only: traditional_to_aux
  use mod_excitation, only: excitation_pi
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  logical, intent(in) :: ifsame_under, ifdvzero, ifnuform

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
       nextrasum, ifpi_local, ifpl,&
       ifsame_abundances,&
       ifsame_zero_abundances
  integer, intent(in) :: ifdv(:), ifelement(:),&
       partial_ions(:), inv_ion(:),&
       partial_elements(:)
  integer, intent(in) :: inv_aux(:)

  integer, intent(out) :: info

  real(fp_kind), intent(in) :: bi(:), h2diss, t, tl, eps(:), eta
  real(fp_kind), intent(in) :: rhostar(:), pstar(:)
  real(fp_kind), intent(in) ::&
       dve_exchange, dve_exchangef, dve_exchanget,&
       plop(:), plopt(:), plopt2(:),&
       dv_pl(:), dv_plt(:),&
       full_sum0, full_sum1, full_sum2,&
       r_ion3(:), r_neutral(:),&
       charge(:),&
       tc2,&
       f, wf, n_e,&
       bmin(:)

  real(fp_kind), intent(out) ::&
       dvzero(:),&
       dv(:), dvf(:), dvt(:),&
       lambda, gamma_e,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       ne, nef, net,&
       sion, sionf, siont, uion,&
       h2, h2f, h2t,&
       aux(:), auxf(:), auxt(:)

  ! Internal variables
  integer nionsp2, naux, nelements, n_partial_elements
  logical ifcsum, ifextrasum, ifxextrasum, ifnr03, ifnr13

  integer ion, index_ion, ielement, ion0, index, ifnr

  real(fp_kind) dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22,&
       dve, dvef, dvet, dve0, dve2,&
       dve_pi, dve_pif, dve_pit,&
       dve_coulomb, dve_coulombf, dve_coulombt,&
       sum0, sum0ne, sum0f, sum0t,&
       sum2, sum2ne, sum2f, sum2t,&
       fion, fionf,&
       rt, rf,&
       rl,&
       sumpl0, sumpl0f,&
       h_ion, h_ionf, h_iont,&
       he_ion, he_ionf, he_iont,&
       he_ion2, he_ion2f, he_ion2t,&
       h2plus, h2plusf, h2plust

  real(fp_kind), allocatable ::&
       extrasum(:),&
       extrasumf(:),&
       extrasumt(:),&
       extrasum_unused(:),&
       extrasumf_unused(:),&
       extrasumt_unused(:),&
       extrasum_dv_unused(:,:),&
       dv_aux_unused(:, :),&
       sum0_dv_unused(:),&
       sum2_dv_unused(:),&
       fion_dv_unused(:),&
       r_dv_unused(:),&
       sumpl0_dv_unused(:),&
       xextrasum(:),&
       xextrasumf(:),&
       xextrasumt(:),&
       xextrasum_unused(:),&
       xextrasumf_unused(:),&
       xextrasumt_unused(:),&
       xextrasum_dv_unused(:,:),&
       h_ion_dv_unused(:),&
       he_ion_dv_unused(:),&
       he_ion2_dv_unused(:),&
       h2_dv_unused(:),&
       h2plus_dv_unused(:),&
       aux_dv_unused(:, :)

  integer, parameter :: nionsp2_local = 318
  integer, parameter :: naux_local = 21
  ! maximum core charge for non-bare ion
  integer, parameter :: maxcorecharge = 28
  integer, parameter :: nxextrasum_local = 4
  integer, parameter :: maxnextrasum_local = 9
  ! iextraoff is the number of non-extrasum auxiliary variables
  integer, parameter :: iextraoff = 7
  integer, parameter :: nstar = 9

  nionsp2 = size(nion)
  naux = size(inv_aux)
  nelements = size(ion_end) - 2
  n_partial_elements = size(partial_elements) - 2

  ! Sanity checks
  if(nionsp2.ne.nionsp2_local) error stop 'eos_cold_start: nionsp2 must be equal to nionsp2_local'
  if(naux.ne.naux_local) error stop 'eos_cold_start: naux must be equal to naux_local'
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
       nionsp2.ne.size(dv).or.&
       nionsp2.ne.size(dvf).or.&
       nionsp2.ne.size(dvt))&
       error stop 'eos_cold_start: inconsistent nionsp2 array sizes'
  if(&
       naux.ne.size(aux).or.&
       naux.ne.size(auxf).or.&
       naux.ne.size(auxt))&
       error stop 'eos_cold_start: inconsistent naux array sizes'
  if(&
       nelements.ne.size(iatomic_number).or.&
       nelements.ne.size(ifelement).or.&
       nelements.ne.size(eps).or.&
       nelements+2.ne.size(r_neutral).or.&
       nelements.ne.size(dvzero))&
       error stop 'eos_cold_start: inconsistent nelements array sizes'
  if(&
       maxcorecharge.ne.size(nmin).or.&
       maxcorecharge.ne.size(nmin_max).or.&
       maxcorecharge.ne.size(bmin))&
       error stop 'eos_cold_start: inconsistent maxcorecharge array sizes'
  if(nstar.ne.size(rhostar).or.nstar.ne.size(pstar))&
       error stop 'eos_cold_start: inconsistent nstar array sizes'
  if(nextrasum.gt.maxnextrasum_local.or.nextrasum.gt.size(extrasum_scale))&
       error stop 'eos_cold_start: inconsistent nextrasum array sizes'
  if(nxextrasum_local.ne.size(xextrasum_scale))&
       error stop 'eos_cold_start: inconsistent nxextrasum array sizes'

  allocate(&
       extrasum(nextrasum),&
       extrasumf(nextrasum),&
       extrasumt(nextrasum),&
       extrasum_unused(nextrasum),&
       extrasumf_unused(nextrasum),&
       extrasumt_unused(nextrasum),&
       extrasum_dv_unused(nionsp2, nextrasum),&
       dv_aux_unused(nionsp2, naux),&
       sum0_dv_unused(nionsp2),&
       sum2_dv_unused(nionsp2),&
       fion_dv_unused(nionsp2),&
       r_dv_unused(nionsp2),&
       sumpl0_dv_unused(nionsp2),&
       xextrasum(nxextrasum_local),&
       xextrasumf(nxextrasum_local),&
       xextrasumt(nxextrasum_local),&
       xextrasum_unused(nxextrasum_local),&
       xextrasumf_unused(nxextrasum_local),&
       xextrasumt_unused(nxextrasum_local),&
       xextrasum_dv_unused(nionsp2,nxextrasum_local),&
       h_ion_dv_unused(nionsp2),&
       he_ion_dv_unused(nionsp2),&
       he_ion2_dv_unused(nionsp2),&
       h2_dv_unused(nionsp2),&
       h2plus_dv_unused(nionsp2),&
       aux_dv_unused(nionsp2, naux)&
       )

  if(if_taint_allocated_real) then
     call taint_allocated_real(extrasum)
     call taint_allocated_real(extrasumf)
     call taint_allocated_real(extrasumt)
     call taint_allocated_real(extrasum_unused)
     call taint_allocated_real(extrasumf_unused)
     call taint_allocated_real(extrasumt_unused)
     call taint_allocated_real(extrasum_dv_unused)
     call taint_allocated_real(dv_aux_unused)
     call taint_allocated_real(sum0_dv_unused)
     call taint_allocated_real(sum2_dv_unused)
     call taint_allocated_real(fion_dv_unused)
     call taint_allocated_real(r_dv_unused)
     call taint_allocated_real(sumpl0_dv_unused)
     call taint_allocated_real(xextrasum)
     call taint_allocated_real(xextrasumf)
     call taint_allocated_real(xextrasumt)
     call taint_allocated_real(xextrasum_unused)
     call taint_allocated_real(xextrasumf_unused)
     call taint_allocated_real(xextrasumt_unused)
     call taint_allocated_real(xextrasum_dv_unused)
     call taint_allocated_real(h_ion_dv_unused)
     call taint_allocated_real(he_ion_dv_unused)
     call taint_allocated_real(he_ion2_dv_unused)
     call taint_allocated_real(h2_dv_unused)
     call taint_allocated_real(h2plus_dv_unused)
     call taint_allocated_real(aux_dv_unused)
  endif

  ifnr = 0
  ifnr03 = ifnr.eq.0.or.ifnr.eq.3
  ifnr13 = ifnr.eq.1.or.ifnr.eq.3

  ! Always initialize *_scale to unity for eos_cold_start since by definition,
  ! the physical conditions have changed significantly since the last iterative
  ! determination of *_scale.
  sum0_scale = 1._fp_kind
  sum2_scale = 1._fp_kind
  extrasum_scale(1:nextrasum) = 1._fp_kind
  xextrasum_scale(1:nxextrasum_local) = 1._fp_kind

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
  else
     call pteh_pi(ifmodified,&
          rhostar(1), rhostar(2), rhostar(3), t,&
          dve_pi,&
          dve_pif, dve_pit)
  endif
  ! calculate change in equilibrium constant due to electron degeneracy
  dve = dve_pi - eta
  dvef = dve_pif - wf
  dvet = dve_pit
  ! PTEH (full ionization) approximation to sum0, sum2
  sum0ne = full_sum0/full_sum1
  sum0 = sum0ne*n_e
  sum2ne = full_sum2/full_sum1
  sum2 = sum2ne*n_e
  sum0f = 0._fp_kind
  sum0t = 0._fp_kind
  sum2f = 0._fp_kind
  sum2t = 0._fp_kind
  ! n.b. the 1 in the argument list is to use if_pteh = 1 inside
  ! master_coulomb consistent with the above approximation.  This yields a
  ! consistent calculation of dv and derivatives according to this
  ! free-energy model approximation. However, later in the eos_calc
  ! argument list the actual if_pteh value is used so that the output
  ! sum0, sum2, and derivatives are calculated when needed for further
  ! calculations without this cold-start approximation.
  call master_coulomb(ifnr,&
       rhostar, f,&
       sum0, sum0ne, sum0f, sum0t, sum2, sum2ne, sum2f, sum2t,&
       n_e, t, cpe*pstar(1), pstar, lambda, gamma_e,&
       ifcoulomb_mod, if_dc, 1,&
       dve_coulomb, dve_coulombf, dve_coulombt, dve0, dve2,&
       dv0, dv0f, dv0t, dv2, dv2f, dv2t, dv00, dv02, dv22)
  dve = dve + dve_coulomb
  dvef = dvef + dve_coulombf
  dvet = dvet + dve_coulombt
  ! add in pre-calculated exchange effects
  dve = dve + dve_exchange
  dvef = dvef + dve_exchangef
  dvet = dvet + dve_exchanget
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
  ! apply Planck-Larkin excitation for pteh pressure ionization
  ! N.B. note because ifpi = 1 for this call, the extrasum, extrasumf,
  ! extrasumt, xextrasum, xextrasumf, xextrasumt, and dv_aux arguments
  ! are not read or written inside excitation_pi so all these
  ! arguments can be replaced with their "unused" variants.
  if(ifexcited.gt.0.and.ifpl.eq.1)&
       call excitation_pi(verbosity, ifexcited, ifsame_zero_abundances,&
          ifpl, 1, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,&
          partial_elements, ion_end,&
          tl, izlo, bmin(:izhi), nmin(:izhi), nmin_max(:izhi), nmin_species, nmax,&
          bi, plop, plopt, plopt2,&
          r_ion3, nion, r_neutral,&
          inv_aux, iextraoff,&
          extrasum_unused, extrasumf_unused, extrasumt_unused,&
          xextrasum_unused, xextrasumf_unused, xextrasumt_unused,&
          dv, dvf, dvt, dv_aux_unused)
  ! n.b. since ifnr = 0 (see above) all *_dv variables
  ! are not read or written inside eos_calc so all these
  ! arguments can be replaced with their "unused" variants.
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
       ne, nef, net, sion, sionf, siont, uion,&
       h_ion, h_ionf, h_iont, h_ion_dv_unused,&
       he_ion, he_ionf, he_iont, he_ion_dv_unused,&
       he_ion2, he_ion2f, he_ion2t, he_ion2_dv_unused,&
       fion, fionf, fion_dv_unused,&
       sum0, sum0f, sum0t, sum0_dv_unused,&
       sum2, sum2f, sum2t, sum2_dv_unused,&
       extrasum, extrasumf, extrasumt, extrasum_dv_unused,&
       sumpl0, sumpl0f, sumpl0_dv_unused,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       rl, rf, rt, r_dv_unused,&
       h2, h2f, h2t, h2_dv_unused,&
       h2plus, h2plusf, h2plust, h2plus_dv_unused,&
       xextrasum, xextrasumf, xextrasumt, xextrasum_dv_unused, info)

  if(info.ne.0) return
  ifcsum = .not.(if_mc.eq.1.or.if_pteh.eq.1)
  ifextrasum = ifpi_local.eq.3.or.ifpi_local.eq.4
  ! Sanity check:
  if((ifextrasum.and.nextrasum.ne.maxnextrasum_local).or.(.not.ifextrasum.and.nextrasum.ne.0))&
       error stop 'eos_cold_start: nextrasum and ifpi_local are not compatible'
  ifxextrasum = ifextrasum .and. (ifexcited.gt.0)
  ! Since ifnr.eq.0 (see above), ifnr13 is .false.  Therefore all *_dv
  ! variables are not read or written inside traditional_to_aux so all
  ! these arguments can be replaced with their "unused" variants.
  call traditional_to_aux(&
       if_mc, ifcsum, ifextrasum, ifxextrasum, ifnr03, ifnr13, iextraoff, iextraoff + maxnextrasum_local, max_index,&
       h_ion, h_ionf, h_iont, h_ion_dv_unused,&
       he_ion, he_ionf, he_iont, he_ion_dv_unused,&
       he_ion2, he_ion2f, he_ion2t, he_ion2_dv_unused,&
       rl, rf, rt, r_dv_unused,&
       h2plus, h2plusf, h2plust, h2plus_dv_unused,&
       sum0, sum0f, sum0t, sum0_dv_unused,&
       sum2, sum2f, sum2t, sum2_dv_unused,&
       extrasum, extrasumf, extrasumt, extrasum_dv_unused,&
       xextrasum, xextrasumf, xextrasumt, xextrasum_dv_unused,&
       aux, auxf, auxt, aux_dv_unused)

end subroutine eos_cold_start
