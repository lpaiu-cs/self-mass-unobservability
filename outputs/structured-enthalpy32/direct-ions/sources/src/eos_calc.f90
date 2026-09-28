!*******************************************************************************
!    Copyright (C) 1996-2022 Alan W. Irwin
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

!> This eos_calc subroutine calculates number densities (optionally in
!> n form or else in nu form where nu = n/(rho*avogadro)) and
!> auxiliary variables (appropriate weighted sums over those
!> quantities) for all species included in the FreeEOS free-energy
!> model.  These auxiliary variables are ultimately used in an
!> iterative way to calculate chemical potentials according to some
!> free energy model.  The appropriate difference in chemical
!> potentials are then combined to form the dv quantities, the change
!> in the equilibrium constant of the non-reference species (ion or
!> molecule) relative to the un-ionized monatomic reference state and
!> the free electron
!>
!> dv = (partial F/partial nref - partial F/partial nspecies - ion*partial F/partial ne)/kT.
!>
!> These dv values (held in an array) are required input to
!> eos_calc.
!>
!> The free energy model:
!> it is the responsibility of the calling programme to calculate the
!> dv quantities for the various ions and molecules in a consistent manner.  It
!> should be noted that for a cold start of the EOS iteration, the
!> free-energy terms corresponding to pressure ionization and the
!> Coulomb effect are approximated as functions of only ne, so the dv
!> array elements are quite similar to each other.  (The electron
!> exchange term is always just a function of ne.) Later in the cold
!> start iteration or directly if a Taylor series approximation is
!> used to initialize the auxiliary variables more complicated models
!> are typically used for the Coulomb effect and pressure ionization
!> so the dv array elements can become more varied.
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
!> \param[in] ifexcited
!>   ifexcited > 0 means use excited states (must have Planck-Larkin or ifpi = 3 or 4).<br>
!>    0 < ifexcited < 10 means use approximation to explicit summation.<br>
!>   10 < ifexcited < 20 means use explicit summation.<br>
!>   mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules) and helium.<br>
!>   mod(ifexcited,10) = 2 same as 1 + H2 and H2+.<br>
!>   mod(ifexcited,10) = 3 same as 2 + partially ionized metals.<br>
!> \param[in] ifsame_zero_abundances
!>  ifsame_zero_abundances = 1 means this call has the same zero abundance pattern as the last call.<br>
!>  ifsame_zero_abundances /= 1 means this call has a different zero abundance pattern then the last call.<br>
!> \param[in] ifnuform
!>   ifnuform .eqv. .true. means return the nu = n/(rho*NA) form of all number densities and most auxiliary variables.<br>
!>   ifnuform .eqv. .false. means return the n form of all number densities and most auxiliary variables.<br>
!> \param[in] ifsame_under
!>   ifsame_under .eqv. .true. means use same underflow zeroing for each
!>     ion as in previous call.  This option is useful for removing
!>     small discontinuities caused by variations in the underflow
!>     zeroing which sometimes foil the last stages of
!>     convergence.<br>
!>   ifsame_under .eqv. .false. means calculate underflow zeroing.<br>
!> \param[in] ifnr
!>   ifnr = 0 means calculate partial derivatives of nuvar (communicated via
!>     the mod_nuvar module), hne, sum0, sum2, s, extrasum, h_ion, sh,
!>     he_ion, he_ion2, h_neutral, h_ion_equil, and sumpl1 wrt f, t and
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
!> \param[in] ifpi
!>   ifpi = 0, 1, 2, 3, 4 corresponds to using no, pteh, geff, mhd, or
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
!> \param[in] dvzero
!>   dvzero(nelements_in) is an array of zero-point shifts (zero for
!>     hydrogen) that must be added to dv in all uses below.
!>     Sometimes (high density, high ionization) some of these
!>     quantities can be larger than 1.e9_fp_kind so if we avoid using
!>     it in differences of dv quanties we can gain many significant
!>     digits.<br>
!> \param[in] dv
!>   dv(nions_inp2) is an array of equilibrium constants (explained above).<br>
!> \param[in] dvf
!>   dvf(nions_inp2) is an array of the partial derivative of dv wrt fl.<br>
!> \param[in] dvt
!>   dvt(nions_inp2) is an array of the partial derivative of dv wrt ln t.<br>
!> \param[in] nion
!>   nion(nions_inp2) is an array of ion charge in ion order (must be same order
!>     as bi) e.g., for H+, He+, He++, etc.
!> \param[in] re
!>   re is a scaled fermi-dirac integral such that n_e/NA = cd*re,
!>   where n_e is the number density of free electrons.<br>
!> \param[in] ref ref is the partial derivative of ln re (or ln n_e) wrt fl.<br>
!> \param[in] ret ret is the partial derivative of ln re (or ln n_e) wrt ln t.<br>
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
!> \param[out] h_ion
!>   h_ion is the calculated hydrogen ionization fraction (n(H+)/(n(H)+n(H+)))
!>     times eps(1).  If no molecular formation of H2 and H2+, then
!>     n(H)+n(H+) = eps(1)*rho*NA, and h_ion = n(H+)/(rho*NA).<br>
!> \param[out] h_ionf h_ionf is the calculated partial derivative of h_ion wrt fl.<br>
!> \param[out] h_iont h_iont is the calculated partial derivative of h_ion wrt ln t.<br>
!> \param[out] h_ion_dv h_ion_dv(nions_inp2) is an array of the calculated partial derivatives of h_ion wrt the elements of dv.<br>
!> \param[out] he_ion
!>   he_ion is the calculated helium first ionization fraction n(He+)/(rho*NA).<br>
!> \param[out] he_ionf he_ionf is the calculated partial derivative of he_ion wrt fl.<br>
!> \param[out] he_iont he_iont is the calculated partial derivative of he_ion wrt ln t.<br>
!> \param[out] he_ion_dv he_ion_dv(nions_inp2) is an array of the
!> calculated partial derivatives of he_ion wrt the elements of
!> dv.<br>
!> \param[out] he_ion2
!>   he_ion2 is the calculated helium second ionization fraction n(He++)/(rho*NA).<br>
!> \param[out] he_ion2f he_ion2f is the calculated partial derivative of he_ion2 wrt fl.<br>
!> \param[out] he_ion2t he_ion2t is the calculated partial derivative of he_ion2 wrt ln t.<br>
!> \param[out] he_ion2_dv he_ion2_dv(nions_inp2) is an array of the
!> calculated partial derivatives of he_ion2 wrt the elements of
!> dv.<br>
!> \param[out] fion
!>   fion (returned only if ifpi = 3 or 4) is the calculated ion
!>     free_energy/(rho*kT*NA) (= tc2*u - s at equilibrium).<br>
!> \param[out] fionf fionf is the calculated partial derivative of fion wrt fl.<br>
!> \param[out] fion_dv fion_dv(nions_inp2) is an array of the calculated partial derivatives of fion wrt the elements of dv.<br>
!> \param[out] sum0
!>   sum0 (returned only if ifcsums is .true.) is the calculated sum of the
!>     non-hydrogen ion nu values.  The hydrogen component of this
!>     auxiliary variable is added by the eos_calc calling
!>     routine.<br>
!> \param[out] sum0f sum0f is the calculated partial derivative of sum0 wrt fl.<br>
!> \param[out] sum0t sum0t is the calculated partial derivative of sum0 wrt ln t.<br>
!> \param[out] sum0_dv sum0_dv(nions_inp2) is an array of the calculated partial derivatives of sum0 wrt the elements of dv.<br>
!> \param[out] sum2
!>   sum2 (returned only if ifcsums is .true.) is the calculated sum of the
!>     non-hydrogen ion nu values weighted by the square of the number
!>     of positive charges for each of those ions.  The hydrogen
!>     component of this auxiliary variable is added by the eos_calc
!>     calling routine.<br>
!> \param[out] sum2f sum2f is the calculated partial derivative of sum2 wrt fl.<br>
!> \param[out] sum2t sum2t is the calculated partial derivative of sum2 wrt ln t.<br>
!> \param[out] sum2_dv sum2_dv(nions_inp2) is an array of the calculated partial derivatives of sum2 wrt the elements of dv.<br>
!> \param[out] extrasum
!>   extrasum(nextrasum) (returned only if ifpi = 3 or 4) are calculated
!>     auxiliary variables consisting of weighted sums over non-H
!>     nu(i) = n(i)/(rho*NA).  For iextrasum = 1,nextrasum-2, sum is
!>     only over neutral species with weight of
!>     r_neutral^{iextrasum-1}.  For iextrasum = nextrasum-1, sum is
!>     over all ionized species including bare nucleii but excluding
!>     free electrons with weight of Z^1.5. For iextrasum = nextrasum,
!>     sum is over all species excluding bare nucleii and free
!>     electrons with weight of rion^3.<br>
!> \param[out] extrasumf
!>   extrasumf(nextrasum) is the calculated partial derivative of extrasum wrt fl.<br>
!> \param[out] extrasumt
!>   extrasumt(nextrasum) is the calculated partial derivative of extrasum wrt ln t.<br>
!> \param[out] extrasum_dv
!>   extrasum_dv(nions_inp2, nextrasum) is an array of the calculated partial derivatives of
!>     extrasum wrt the elements of dv.<br>
!> \param[out] sumpl0
!>   sumpl0 (returned only if ifpl = 1) is the calculated weighted sum over non-H
!>     nu(i) = n(i)/(rho*NA) with weight of w where w is the Planck-Larkin
!>     occupation probability.<br>
!> \param[out] sumpl0f
!>  sumpl0f is the calculated partial derivative of sumpl0 wrt fl.<br>
!> \param[out] sumpl0_dv
!>  sumpl0_dv(nions_inp2) is an array of the calculated partial derivatives of sumpl0 wrt the elements of dv.<br>
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
!> \param[out] rl rl is the calculated ln mass density.<br>
!> \param[out] rf rf is the calculated partial derivative of rl wrt fl.<br>
!> \param[out] rt rt is the calculated partial derivative of rl wrt ln t.<br>
!> \param[out] r_dv r_dv(nions_inp2) is an array of the calculated partial derivatives of rl wrt the elements of dv.<br>
!> \param[out] h2
!>   h2 is the calculated number density (optionally in n form or else
!>     in nu form where nu = n/(rho*avogadro)) of H2.<br>
!> \param[out] h2f h2f is the calculated partial derivative of h2 wrt fl.<br>
!> \param[out] h2t h2t is the calculated partial derivative of h2 wrt ln t.<br>
!> \param[out] h2_dv h2_dv(nions_inp2) is an array of the calculated
!>   partial derivatives of h2 wrt the elements of dv.<br>
!> \param[out] h2plus
!>   h2plus is the calculated number density (optionally in n form or else
!>     in nu form where nu = n/(rho*avogadro)) of H2+.<br>
!> \param[out] h2plusf h2plusf is the calculated partial derivative of h2plus wrt fl.<br>
!> \param[out] h2plust h2plust is the calculated partial derivative of h2plus wrt ln t.<br>
!> \param[out] h2plus_dv h2plus_dv(nions_inp2) is an array of the
!>   calculated partial derivatives of h2plus wrt the elements of dv.<br>
!> \param[out] xextrasum
!>   xextrasum(nxextrasum = 4) (returned only if ifexcited.gt.0 and ifpi.eq.3.or4)
!>     are calculated auxiliary variables corresponding to the
!>     *negative* sum nuvar/(1 + qratio)* partial qratio/partial
!>     extrasum(k), for k = 1, 2, 3, and nextrasum-1.  These excitation auxiliary
!>     variables are unique in that they depend on other auxiliary variables (i.e., extrasum) so
!>     must be calculated with a call to excitation_sum after extrasum is calculated.
!> \param[out] xextrasumf
!>   xextrasumf(nxextrasum = 4) is the calculated partial derivative of xextrasum wrt fl.<br>
!> \param[out] xextrasumt
!>   xextrasumt(nxextrasum = 4) is the calculated partial derivative of xextrasum wrt ln t.<br>
!> \param[out] xextrasum_dv
!>   xextrasum_dv(nions_inp2,nxextrasum = 4) is an array of the
!>   calculated partial derivatives of xextrasum wrt the elements of dv.<br>
!> \param[out] info
!>   info is the calculated return code for the procedure which is
!>     zero for success and non-zero (with offset to identify which
!>     procedure the error occurred in) if an error occurred within
!>     the procedure or within some procedure it called.<br>
subroutine eos_calc(verbosity, ifexcited, ifsame_zero_abundances,&
     ifnuform, ifsame_under, ifnr,&
     inv_ion, max_index,&
     partial_elements, ion_end,&
     ifionized, if_pteh, if_mc, ifreducedmass,&
     ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,&
     ifmodified, ifh2, ifh2plus,&
     izlo, bmin, nmin, nmin_max, nmin_species, nmax,&
     eps, tl, tc2, bi, h2diss, plop, plopt, plopt2,&
     r_ion3, r_neutral,&
     ifelement, dvzero,&
     dv, dvf, dvt,&
     nion,&
     re, ref, ret,&
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
     h2, h2f, h2t, h2_dv,&
     h2plus, h2plusf, h2plust, h2plus_dv,&
     xextrasum, xextrasumf, xextrasumt, xextrasum_dv,&
     info)

  use mod_free_eos_constants, only: logalpha, avogadro, cd, electron_mass, ln10
  use mod_molecular_hydrogen, only: molecular_hydrogen
  use mod_nuvar, only: nuvar
  use mod_aux_scale, only: aux_scale_limit_factor, aux_underflow,&
       sum0_scale, sum2_scale, extrasum_scale, xextrasum_scale
  use mod_statistical_weight_data, only: iqion, iqneutral
  use mod_excitation_block, only: nx
  use mod_excitation, only: excitation_sum
  use mod_flow_data, only: underflow_limit, ln_underflow_limit, ln_overflow_limit
  use mod_info_data, only: info_offset_eos_calc
  use mod_isotopic_mass_data, only: nelements_iso, isotopic_mass
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  logical, intent(in) :: ifnuform, ifsame_under
  integer, intent(in) :: verbosity, ifexcited, ifsame_zero_abundances, ifnr,&
       partial_elements(:),&
       ifionized, if_pteh, if_mc, ifreducedmass,&
       ifsame_abundances, ifmtrace, ifpi, ifpl,&
       ifmodified, ifh2, ifh2plus,&
       izlo, nmin(:), nmin_max(:),&
       nmin_species(:), nmax,&
       nion(:), inv_ion(:), max_index,&
       iatomic_number(:), ion_end(:), ifelement(:)

  integer, intent(out) :: info

  real(fp_kind), intent(in) :: eps(:), tl, tc2, bi(:), h2diss,&
       plop(:), plopt(:), plopt2(:), bmin(:),&
       dvzero(:), dv(:), dvf(:), dvt(:),&
       re, ref, ret,&
       r_ion3(:), r_neutral(:)
  real(fp_kind), intent(out) ::&
       ne, nef, net, sion, sionf, siont, uion,&
       h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t,&
       h_ion_dv(:), he_ion_dv(:), he_ion2_dv(:),&
       fion, fionf, fion_dv(:),&
       sum0, sum0f, sum0t, sum0_dv(:),&
       sum2, sum2f, sum2t, sum2_dv(:),&
       extrasum(:), extrasumf(:), extrasumt(:), extrasum_dv(:,:),&
       xextrasum(:), xextrasumf(:), xextrasumt(:), xextrasum_dv(:,:),&
       sumpl0, sumpl0f, sumpl0_dv(:),&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       rl, rf, rt, r_dv(:),&
       h2, h2f, h2t, h2plus, h2plusf, h2plust,&
       h2_dv(:), h2plus_dv(:)

  ! Parameters
  integer, parameter :: nxextrasum = 4

  real(fp_kind), parameter :: ln_avogadro = log(avogadro)
  ! Parameter associated with finding the dominant H species.
  ! lnerrcrit is set so that ln(1 + exp(-lnerrcrit)) ~
  ! exp(-lnerrcrit) = 10**(-precision(1._fp_kind)) ==> lnerrcrit = precision(1._fp_kind)*ln(10._fp_kind)
  real(fp_kind), parameter :: lnerrcrit = real(precision(1._fp_kind), kind=fp_kind)*ln10
  real(fp_kind), parameter :: rho_perturb_factor = 1._fp_kind! + 1.e-15_fp_kind

  real(fp_kind), parameter :: isotopic_mass_h = isotopic_mass(1)
  real(fp_kind), parameter :: logmass_h = log(isotopic_mass_h)
  real(fp_kind), parameter :: logmass_h2 = log(2._fp_kind*isotopic_mass_h)
  real(fp_kind), parameter :: h2_log = 1.5_fp_kind*(logmass_h2 - 2._fp_kind*logmass_h)
  ! to be consistent with the slightly lame free-energy model
  ! underlying the entropy and energy calculation which assumes
  ! the mass of H2 and H2+ are equal to each other.
  real(fp_kind), parameter :: logmass_h2plus = logmass_h2
  ! (unused) real(fp_kind), parameter ::  h2plus_log = 1.5_fp_kind*(logmass_h2plus - logmass_h2)
  ! drop partition function from all molecular values since that not
  ! constant with temperature.
  real(fp_kind), parameter :: logqtl_const_h2plus = 1.5_fp_kind*logmass_h2plus - logalpha
  real(fp_kind), parameter :: logqtl_const_h2 = 1.5_fp_kind*logmass_h2 - logalpha
  real(fp_kind), parameter :: logqtl_const_h = log(real((iqneutral(1)),fp_kind)) + 1.5_fp_kind*logmass_h - logalpha
  real(fp_kind), parameter :: logmass_hplus_reduced = log(isotopic_mass_h - electron_mass)
  real(fp_kind), parameter :: logmass_hplus_unreduced = log(isotopic_mass_h)
  real(fp_kind), parameter :: logqtl_const_hplus_reduced =&
       log(real((iqion(1)),fp_kind)) + 1.5_fp_kind*logmass_hplus_reduced - logalpha
  real(fp_kind), parameter :: logqtl_const_hplus_unreduced =&
       log(real((iqion(1)),fp_kind)) + 1.5_fp_kind*logmass_hplus_unreduced - logalpha

  ! Local variables
  logical ifabh, ifabg, ifaah2plus, ifaah2
  logical ifnr03, ifnr13
  logical ifpi34, ifnuform_old, ifcsums

  integer iextrasum,  ixextrasum, jextrasum,iffirst
  data iffirst/1/
  integer index_maxh
  integer index, index2, index3, indexh, indexh2, indexh2plus,&
       nions, nionsp2, nelements, nextrasum

  real(fp_kind) fionh, fionhf
  real(fp_kind) sharg,&
       h_neutral, h_neutralf, h_neutralt,&
       ! nh* variables/arrays for temporary storage of n forms...
       nh2, nh2f, nh2t,&
       nh2plus, nh2plusf, nh2plust,&
       nh_neutral, nh_neutralf, nh_neutralt,&
       nh_ion, nh_ionf, nh_iont,&
       log_h_neutral,&
       sh, shf, sht, hne, hnef, hnet,&
       en, rmue, rho,&
       tlold, qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt,&
       h_ion_equil, h_ion_equilf, h_ion_equilt, h_ion_equil_dv,&
       h2equil, h2equilf, h2equilt, h2equilt0, h2equil_dv,&
       h2plusequil, h2plusequilf, h2plusequilt, h2plusequil_dv,&
       ablog, ablogf, ablogt,&
       abh, abhf, abht, abh_dv,&
       abg, abgf, abgt, abg_dv, abexp,&
       abdlog, abdlogf, abdlogt,&
       abdh, abdhf, abdht, abdh_dv,&
       abdg, abdgf, abdgt, abdg_dv,&
       ac, acf, act,&
       aalog, aalogf, aalogt, aalog_dv,&
       aah2plus, aah2plusf, aah2plust, aah2plus_dv,&
       aah2, aah2f, aah2t, aah2_dv, aaexp,&
       cprime, cprimef, cprimet, cprime_dv,&
       xprimelog, xprimelogf, xprimelogt,&
       rpower, rpowerh2, rpowerh2plus, logqtl_const_hplus, constant_sh

  ! something ridiculous
  data tlold/1.e-300_fp_kind/
  real(fp_kind) maxh, second_maxh

  real(fp_kind), allocatable ::&
       lextrasum(:),&
       lextrasumf(:),&
       lextrasumt(:),&
       lextrasum_dv(:, :),&
       fionh_dv(:),&
       h_neutral_dv(:),&
       ! nh* variables/arrays for temporary storage of n forms...
       nh2_dv(:),&
       nh2plus_dv(:),&
       nh_neutral_dv(:),&
       nh_ion_dv(:),&
       hne_dv(:),&
       ablog_dv(:),&
       abdlog_dv(:),&
       xprimelog_dv(:)

  ! Most/all Fortran compilers specify the save attribute for all variables
  ! intialized by data statements.  But just in case...
  save iffirst, tlold

  ! Save ifnuform_old for obvious reasons.
  save ifnuform_old

  ! Save variables associated with tlold
  save qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt

  ! Assume success by default unless an error condition is discovered below.
  info = 0

  nionsp2 = size(nion)
  nions = nionsp2 - 2
  nelements = size(iatomic_number)
  nextrasum = size(extrasum)

  ! sanity checking:
  if(if_pteh.eq.1.and.if_mc.eq.1) error stop 'eos_calc: if_pteh and if_mc cannot simultaneously be non-zero'
  if(nelements.ne.nelements_iso) error stop 'eos_calc: inconsistent nelements and nelements_iso'
  if(ifnr.lt.0.or.ifnr.eq.2.or.ifnr.gt.3) error stop 'eos_calc: ifnr must be 0, 1, or 3'
  ! ifnuform should be true only on last call to eos_calc which coincides with ifnr = 0.
  if(ifnuform.and.ifnr.ne.0) error stop 'eos_calc: bad combination of ifnuform and ifnr'
  if(ifionized.eq.2.and.ifexcited.gt.0) error stop 'eos_calc: bad combination of ifionized and ifexcited'
  if(iffirst.eq.1.and.ifsame_abundances.eq.1) error stop 'eos_calc: invalid ifsame_abundances on first call'
  if(size(nmin).ne.size(nmin_max).or.size(nmin).ne.size(bmin))&
       error stop 'eos_calc: inconsistent sizes for nmin, nmin_max, or bmin'
  if(&
       nionsp2.ne.size(nmin_species).or.&
       nionsp2.ne.size(inv_ion).or.&
       nionsp2.ne.size(bi).or.&
       nionsp2.ne.size(plop).or.&
       nionsp2.ne.size(plopt).or.&
       nionsp2.ne.size(plopt2).or.&
       nionsp2.ne.size(dv).or.&
       nionsp2.ne.size(dvf).or.&
       nionsp2.ne.size(dvt).or.&
       nionsp2.ne.size(r_ion3).or.&
       nionsp2.ne.size(h_ion_dv).or.&
       nionsp2.ne.size(he_ion_dv).or.&
       nionsp2.ne.size(he_ion2_dv).or.&
       nionsp2.ne.size(fion_dv).or.&
       nionsp2.ne.size(sum0_dv).or.&
       nionsp2.ne.size(sum2_dv).or.&
       nionsp2.ne.size(extrasum_dv,1).or.&
       nionsp2.ne.size(sumpl0_dv).or.&
       nionsp2.ne.size(r_dv).or.&
       nionsp2.ne.size(h2_dv).or.&
       nionsp2.ne.size(h2plus_dv))&
       error stop 'eos_calc: inconsistent sizes for 22 variables'
  if(&
       nelements+2.ne.size(ion_end).or.&
       nelements.ne.size(ifelement).or.&
       nelements.ne.size(eps).or.&
       nelements.ne.size(dvzero).or.&
       nelements+2.ne.size(r_neutral))&
       error stop 'eos_calc: inconsistent sizes for iatomic_number, ion_end, ifelement, eps, dvzero, or r_neutral'
  if(&
       nextrasum.gt.size(extrasum_scale).or.&
       nextrasum.ne.size(extrasumf).or.&
       nextrasum.ne.size(extrasumt).or.&
       nextrasum.ne.size(extrasum_dv,2))&
       error stop 'eos_calc: inconsistent nextrasum sizes'
  if(&
       size(xextrasum_scale).ne.nxextrasum.or.&
       size(xextrasum).ne.nxextrasum.or.&
       size(xextrasumf).ne.nxextrasum.or.&
       size(xextrasumt).ne.nxextrasum.or.&
       size(xextrasum_dv,1).ne.nionsp2.or.&
       size(xextrasum_dv,2).ne.nxextrasum)&
       error stop 'eos_calc: bad sizes for xextrasum, xextrasumf, xextrasumt, or xextrasum_dv'

  if(iffirst.eq.1) then
     iffirst = 0
     ! Always initialize *_scale factors to unity for first call.
     sum0_scale = 1._fp_kind
     sum2_scale = 1._fp_kind
     extrasum_scale(1:nextrasum) = 1._fp_kind
     xextrasum_scale(1:nxextrasum) = 1._fp_kind
     ifnuform_old = ifnuform
  else
     ! Always initialize aux_scale to unity if changed ifnuform
     if(ifnuform.neqv.ifnuform_old) then
        sum0_scale = 1._fp_kind
        sum2_scale = 1._fp_kind
        extrasum_scale(1:nextrasum) = 1._fp_kind
        xextrasum_scale(1:nxextrasum) = 1._fp_kind
        ifnuform_old = ifnuform
     endif
  endif

  ! n.b. the atomic masses and the ground statistical weight for H and H+
  ! used below *must* be identical with the values used in ionize.
  ! Otherwise, you will introduce a discontinuity when switching from
  ! negligible molecular hydrogen to a calculation without molecular
  ! hydrogen.
  if(ifreducedmass.eq.1) then
     logqtl_const_hplus = logqtl_const_hplus_reduced
  else
     logqtl_const_hplus = logqtl_const_hplus_unreduced
  endif
  ! calculate constant component of hydrogen entropy to be used
  ! in molecular case only (see below).
  constant_sh = eps(1)*logqtl_const_h
  ifnr03 = ifnr.eq.0.or.ifnr.eq.3
  ifnr13 = ifnr.eq.1.or.ifnr.eq.3
  ifpi34 = ifpi.eq.3.or.ifpi.eq.4

  allocate(&
       lextrasum(nextrasum),&
       lextrasumf(nextrasum),&
       lextrasumt(nextrasum),&
       lextrasum_dv(nionsp2, nextrasum),&
       fionh_dv(nionsp2),&
       h_neutral_dv(nionsp2),&
       ! nh* variables/arrays for temporary storage of n forms...
       nh2_dv(nionsp2),&
       nh2plus_dv(nionsp2),&
       nh_neutral_dv(nionsp2),&
       nh_ion_dv(nionsp2),&
       hne_dv(nionsp2),&
       ablog_dv(nionsp2),&
       abdlog_dv(nionsp2),&
       xprimelog_dv(nionsp2))

  if(if_taint_allocated_real) then
     call taint_allocated_real(lextrasum)
     call taint_allocated_real(lextrasumf)
     call taint_allocated_real(lextrasumt)
     call taint_allocated_real(lextrasum_dv)
     call taint_allocated_real(fionh_dv)
     call taint_allocated_real(h_neutral_dv)
     ! nh* variables/arrays for temporary storage of n forms...
     call taint_allocated_real(nh2_dv)
     call taint_allocated_real(nh2plus_dv)
     call taint_allocated_real(nh_neutral_dv)
     call taint_allocated_real(nh_ion_dv)
     call taint_allocated_real(hne_dv)
     call taint_allocated_real(ablog_dv)
     call taint_allocated_real(abdlog_dv)
     call taint_allocated_real(xprimelog_dv)
  endif

  ! Analysis of the ionize subroutine shows that if ifionized.ne.2
  ! then the call to ionize defines h_ion_equil, h_neutral, and sh and
  ! (if ifnr03 is additionally .true.)  h_ion_equilf, h_ion_equilft,
  ! h_neutralf, h_neutralt, sh, shf, and sht.  Further analysis of the
  ! current eos_calc code shows those variables are only used in code
  ! blocks after that call where those same conditions are met.
  ! Note that the gfortran -Wmaybe-uninitialized option implied by
  ! gfortran -Wall generates a false warning for this case so to avoid
  ! that could uncomment the following lines, but better yet use the
  ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
  ! code analysis to determine uninitialized variables and instead
  ! rely on run-time analysis (e.g., tainting with NaN's) to discover
  ! actual uninitialized values.

  !h_ion_equil = 0._fp_kind
  !h_ion_equilf = 0._fp_kind
  !h_ion_equilt = 0._fp_kind
  !h_neutral = 0._fp_kind
  !h_neutralf = 0._fp_kind
  !h_neutralt = 0._fp_kind
  !sh = 0._fp_kind
  !shf = 0._fp_kind
  !sht = 0._fp_kind

  ! N.B. sum[02]-related quantities calculated in ionize (and consistently below)
  ! only if ifcsums is .true.
  ifcsums = .not.(if_mc.eq.1.or.if_pteh.eq.1)
  ifcsums = .true.

  call ionize(verbosity, ifexcited, ifsame_under, ifnr, inv_ion,&
       max_index, partial_elements, ion_end,&
       ifionized, ifcsums, ifreducedmass,&
       ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,&
       eps, tc2, bi, plop, plopt, plopt2,&
       r_ion3, r_neutral,&
       ifelement, dvzero,&
       dv, dvf, dvt,&
       nion,&
       hne, hnef, hnet, hne_dv,&
       sion, sionf, siont,&
       uion,&
       h_neutral, h_neutralf, h_neutralt,&
       sh, shf, sht,&
       h_ion_equil, h_ion_equilf, h_ion_equilt,&
       h_ion, h_ionf, h_iont, h_ion_dv,&
       he_ion, he_ionf, he_iont, he_ion_dv,&
       he_ion2, he_ion2f, he_ion2t, he_ion2_dv,&
       fionh, fionhf, fionh_dv,&
       fion, fionf, fion_dv,&
       sum0, sum0f, sum0t, sum0_dv,&
       sum2, sum2f, sum2t, sum2_dv,&
       extrasum, extrasumf, extrasumt, extrasum_dv,&
       sumpl0, sumpl0f, sumpl0_dv,&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       info)
  if(info.ne.0) return

  ! These variables are defined in the "if(ifnr13) then" code block
  ! below, and are only used afterward in identical code blocks.
  ! Note that the gfortran -Wmaybe-uninitialized option implied by
  ! gfortran -Wall generates a false warning for this case so to avoid
  ! that could uncomment the following lines, but better yet use the
  ! gfortran option -Wno-maybe-uninitialized to drop this aspect of
  ! code analysis to determine uninitialized variables and instead
  ! rely on run-time analysis (e.g., tainting with huge values) to
  ! discover actual uninitialized values.

  !indexh = 0
  !indexh2 = 0
  !indexh2plus = 0

  if(ifnr13) then
     ! required by h_ion_equil
     indexh = inv_ion(1)
     ! required by h2equil *and* h2plusequil
     indexh2 = inv_ion(nions+1)
     ! required by h2plusequil
     indexh2plus = inv_ion(nions+2)
     ! As returned from ionize, h_neutral+h_ion = eps(1), hence h_neutral_dv = -h_ion_dv.
     h_neutral_dv(1:max_index) = -h_ion_dv(1:max_index)
  endif

  ! rho/mu_e = n_e H = n_e/N_A = cd*re = rmue
  rmue = cd*re
  if(ifionized.eq.2) then
     ! Start of logic block for case of everything fully ionized.

     ! ignore ifnr (i.e. forget dv derivatives) since no NR iteration
     ! required with full ionization.
     ! n.b. for full ionization, ionize includes hydrogen effects
     ! in hne, sion, uion, and fion.
     ne = hne

     if(ne.lt.underflow_limit) then
        if(verbosity.ge.1) write(stderr,*)&
             "eos_calc: (1) calculated electron number density (in nu form) is too close to underflowing"
        info = info_offset_eos_calc + 1
        return
     endif

     en = 1._fp_kind/ne
     ! rho/mu_e = n_e H = cd*re = rmue
     rl = log(rmue) + log(en)
     if(rl.lt.ln_underflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_calc: (21) calculated mass density too close to underflow limit"
        info = info_offset_eos_calc + 21
        return
     elseif(rl.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_calc: (31) calculated mass density too close to overflow limit"
        info = info_offset_eos_calc + 31
        return
     endif
     rho = exp(rl)*rho_perturb_factor

     if(ifnr03) then
        nef = 0._fp_kind
        net = 0._fp_kind
        rt = ret-en*net
        rf = ref-en*nef
     endif
     ! calculate fully ionized nu forms of h2, h2plus, h_ion, he_ion, he_ion2
     h2 = 0._fp_kind
     h2plus = 0._fp_kind
     h_ion = eps(1)
     he_ion = 0._fp_kind
     he_ion2 = eps(2)
     if(ifnr03) then
        h2f = 0._fp_kind
        h2t = 0._fp_kind
        h2plusf = 0._fp_kind
        h2plust = 0._fp_kind
        h_ionf = 0._fp_kind
        h_iont = 0._fp_kind
        he_ionf = 0._fp_kind
        he_iont = 0._fp_kind
        he_ion2f = 0._fp_kind
        he_ion2t = 0._fp_kind
     endif
     if(ifnr13) then
        h2_dv(1:max_index) = 0._fp_kind
        h2plus_dv(1:max_index) = 0._fp_kind
        h_ion_dv(1:max_index) = 0._fp_kind
        he_ion_dv(1:max_index) = 0._fp_kind
        he_ion2_dv(1:max_index) = 0._fp_kind
     endif
     ! sum0, sum2, extrasum, sumpl0, sumpl1, sumpl2
     ! not needed.
     ! End of logic block for case of everything fully ionized.
  elseif(eps(1).le.0._fp_kind) then
     ! Start of logic block for case of partially ionized but no hydrogen species of any kind.
     ne = hne

     if(ne.lt.underflow_limit) then
        if(verbosity.ge.1) write(stderr,*)&
             "eos_calc: (2) calculated electron number density (in nu form) is too close to underflowing"
        info = info_offset_eos_calc + 2
        return
     endif

     ! sion and derivatives, fion and derivative, and uion
     ! passed back as is from ionize because no h_ion component.
     en = 1._fp_kind/ne
     ! rho/mu_e = n_e H = cd*re = rmue
     rl = log(rmue) + log(en)
     if(rl.lt.ln_underflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_calc: (22) calculated mass density too close to underflow limit"
        info = info_offset_eos_calc + 22
        return
     elseif(rl.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_calc: (32) calculated mass density too close to overflow limit"
        info = info_offset_eos_calc + 32
        return
     endif
     rho = exp(rl)*rho_perturb_factor

     h2 = 0._fp_kind
     h2plus = 0._fp_kind
     if(ifnr03) then
        nef = hnef
        net = hnet
        rf = ref-en*nef
        rt = ret-en*net
        h2f = 0._fp_kind
        h2t = 0._fp_kind
        h2plusf = 0._fp_kind
        h2plust = 0._fp_kind
     endif
     if(ifnr13) then
        r_dv(1:max_index) = -en*hne_dv(1:max_index)
        h2_dv(1:max_index) = 0._fp_kind
        h2plus_dv(1:max_index) = 0._fp_kind
     endif
     ! h_ion and derivatives irrelevant if zero abundance,
     ! but force to be zero because previous
     ! calls may have had non-zero abundance.
     h_ion = 0._fp_kind
     if(ifpi.eq.2.or.if_mc.eq.1) then
        if(ifnr03) then
           h_ionf = 0._fp_kind
           h_iont = 0._fp_kind
        endif
     endif
     if(eps(2).eq.0._fp_kind) then
        he_ion = 0._fp_kind
        he_ion2 = 0._fp_kind
        if(ifpi.eq.2.or.if_mc.eq.1) then
           if(ifnr03) then
              he_ionf = 0._fp_kind
              he_iont = 0._fp_kind
              he_ion2f = 0._fp_kind
              he_ion2t = 0._fp_kind
           endif
        endif
     endif
     ! return sumpl0, sumpl1, derivatives, and sumpl2 as is from ionize
     ! since no h_ion component and keep in nu = n/(rho*avogadro) form.
     ! End of logic block for case of partially ionized but no hydrogen (or molecular hydrogen).
  elseif(ifh2.eq.0) then
     ! Start of logic block for case of partially ionized including monatomic forms of hydrogen,
     ! but excluding H2 and H2+.

     ! keep ne, sion, sion derivatives, uion,
     ! sumpl0, sumpl0f, sumpl1, sumpl1 derivatives,
     ! and sumpl2 in nu = n/(rho*avogadro) form.
     ne = hne + h_ion

     if(ne.lt.underflow_limit) then
        if(verbosity.ge.1) write(stderr,*)&
             "eos_calc: (3) calculated electron number density (in nu form) is too close to underflowing"
        info = info_offset_eos_calc + 3
        return
     endif

     uion = uion + bi(1)*h_ion
     sion = sion + sh
     fion = fion + fionh
     if(ifpl.eq.1) then
        sumpl0 = sumpl0 + h_neutral*plop(1)
        sumpl1 = sumpl1 + h_neutral*(plop(1) + plopt(1))
        sumpl2 = sumpl2 + h_neutral*plopt(1)
        if(ifnr03) then
           sumpl0f = sumpl0f + h_neutralf*plop(1)
           sumpl1f = sumpl1f + h_neutralf*(plop(1) + plopt(1))
           sumpl1t = sumpl1t + h_neutralt*(plop(1) + plopt(1)) + h_neutral*(plopt(1) + plopt2(1))
        endif
        if(ifnr13) then
           sumpl0_dv(1:max_index) = sumpl0_dv(1:max_index) + h_neutral_dv(1:max_index)*plop(1)
        endif
     endif !if(ifpl.eq.1) then
     en = 1._fp_kind/ne
     ! rho/mu_e = n_e H = cd*re = rmue
     rl = log(rmue) + log(en)
     if(rl.lt.ln_underflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_calc: (23) calculated mass density too close to underflow limit"
        info = info_offset_eos_calc + 23
        return
     elseif(rl.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "eos_calc: (33) calculated mass density too close to overflow limit"
        info = info_offset_eos_calc + 33
        return
     endif
     rho = exp(rl)*rho_perturb_factor

     nuvar(1,1) = h_neutral
     nuvar(2,1) = h_ion
     ! store nuh2 and nuh2plus in the 3rd and 4th index
     nuvar(3,1) = 0._fp_kind
     nuvar(4,1) = 0._fp_kind
     h2 = 0._fp_kind
     h2plus = 0._fp_kind
     if(ifnr03) then
        nef = hnef + h_ionf
        net = hnet + h_iont
        sionf = sionf + shf
        siont = siont + sht
        fionf = fionf + fionhf
        rf = ref - en*nef
        rt = ret - en*net
        h2f = 0._fp_kind
        h2t = 0._fp_kind
        h2plusf = 0._fp_kind
        h2plust = 0._fp_kind
     endif
     if(ifnr13) then
        do index = 1, max_index
           if(index.eq.indexh) then
              r_dv(index) = -en*(hne_dv(index) + h_ion_dv(index))
              fion_dv(index) = fion_dv(index) + fionh_dv(index)
           else
              r_dv(index) = -en*hne_dv(index)
           endif
        enddo !do index = 1, max_index
        h2_dv(1:max_index) = 0._fp_kind
        h2plus_dv(1:max_index) = 0._fp_kind
     endif !if(ifnr13) then
     if(eps(2).eq.0._fp_kind) then
        he_ion = 0._fp_kind
        he_ion2 = 0._fp_kind
        if(ifpi.eq.2.or.if_mc.eq.1) then
           if(ifnr03) then
              he_ionf = 0._fp_kind
              he_ion2f = 0._fp_kind
              he_ionf = 0._fp_kind
              he_ion2t = 0._fp_kind
           endif
        endif
     endif
     if(ifcsums) then
        ! Add in effect of monatomic hydrogen ion species to Coulomb sums.
        sum0 = sum0 + h_ion*sum0_scale
        sum2 = sum2 + h_ion*sum2_scale
        if(ifnr03) then
           sum0f = sum0f + h_ionf*sum0_scale
           sum0t = sum0t + h_iont*sum0_scale
           sum2f = sum2f + h_ionf*sum2_scale
           sum2t = sum2t + h_iont*sum2_scale
        endif
        if(ifnr13) then
           sum0_dv(1:max_index) = sum0_dv(1:max_index) + h_ion_dv(1:max_index)*sum0_scale
           sum2_dv(1:max_index) = sum2_dv(1:max_index) + h_ion_dv(1:max_index)*sum2_scale
        endif
     endif !if(ifcsums) then

     if(ifpi34) then
        ! add in effect of neutral monatomic H
        rpower = 1._fp_kind
        ! Cannot replace following do loop with array ops because of rpower.
        do iextrasum = 1,nextrasum-2
           extrasum(iextrasum) = extrasum(iextrasum) + h_neutral*extrasum_scale(iextrasum)*rpower
           if(ifnr03) then
              extrasumf(iextrasum) = extrasumf(iextrasum) + h_neutralf*extrasum_scale(iextrasum)*rpower
              extrasumt(iextrasum) = extrasumt(iextrasum) + h_neutralt*extrasum_scale(iextrasum)*rpower
           endif
           if(ifnr13) then
              extrasum_dv(1:max_index,iextrasum) = extrasum_dv(1:max_index,iextrasum) +&
                   h_neutral_dv(1:max_index)*extrasum_scale(iextrasum)*rpower
           endif
           rpower = rpower*r_neutral(1)
        enddo
        ! sum over ions with weight of Z^1.5 which is unity for H+.
        extrasum(nextrasum-1) = extrasum(nextrasum-1) + h_ion*extrasum_scale(nextrasum-1)
        ! sum over all neutral and ionized states except bare nuclei
        extrasum(nextrasum) = extrasum(nextrasum) + h_neutral*extrasum_scale(nextrasum)*r_ion3(1)
        if(ifnr03) then
           extrasumf(nextrasum-1) = extrasumf(nextrasum-1) + h_ionf*extrasum_scale(nextrasum-1)
           extrasumt(nextrasum-1) = extrasumt(nextrasum-1) + h_iont*extrasum_scale(nextrasum-1)
           extrasumf(nextrasum) = extrasumf(nextrasum) + h_neutralf*extrasum_scale(nextrasum)*r_ion3(1)
           extrasumt(nextrasum) = extrasumt(nextrasum) + h_neutralt*extrasum_scale(nextrasum)*r_ion3(1)
        endif
        if(ifnr13) then
           extrasum_dv(1:max_index,nextrasum-1) = extrasum_dv(1:max_index,nextrasum-1) +&
                h_ion_dv(1:max_index)*extrasum_scale(nextrasum-1)
           extrasum_dv(1:max_index,nextrasum) = extrasum_dv(1:max_index,nextrasum) +&
                h_neutral_dv(1:max_index)*extrasum_scale(nextrasum)*r_ion3(1)
        endif
     endif !if(ifpi34) then
     ! End of logic block for case of partially ionized including monatomic forms of hydrogen,
     ! but excluding H2 and H2+.
  else !if(ifionized.eq.2) then
     ! Start of logic block for (most common) case of partially ionized including hydrogen monatomics and molecules.
     if(tl.ne.tlold) then
        tlold = tl
        call molecular_hydrogen(verbosity, ifh2, ifh2plus, tl, qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt)
     endif

     ! n(H+)/n(H) = exp(h_ion_equil), where h_ion_equil returned from ionize
     ! n(H2)*avogadro/(n(H)*n(H)) = exp(h2equil)
     if(ifnr13) then
        h_ion_equil_dv = 1._fp_kind
     endif
     ! n.b. keep h2equil in ln form for under/over flow reasons.
     h2equil = logalpha - log(4._fp_kind) + qh2 - 1.5_fp_kind*tl + h2diss*tc2 + h2_log + dv(nions+1)
     ! d h2equil/ d ln T (excluding dv effects)
     h2equilt0 = qh2t - 1.5_fp_kind - h2diss*tc2

     ! These variables are defined in the logical code blocks below
     ! and are only used afterward in identical code blocks.
     ! Note that the gfortran -Wmaybe-uninitialized option implied by
     ! gfortran -Wall generates a false warning for this case so to avoid
     ! that could uncomment the following lines, but better yet use the
     ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
     ! code analysis to determine uninitialized variables and instead
     ! rely on run-time analysis (e.g., tainting with NaN's) to discover
     ! actual uninitialized values.

     !h2equilf = 0._fp_kind
     !h2equilt = 0._fp_kind
     !h2equil_dv = 0._fp_kind

     if(ifnr03) then
        ! d h2equil/ d ln f
        h2equilf = dvf(nions+1)
        ! d h2equil/ d ln T
        h2equilt = h2equilt0 + dvt(nions+1)
     endif
     if(ifnr13) then
        h2equil_dv = 1._fp_kind
     endif

     ! These variables are defined in the "if(ifh2plus.gt.0) then" code
     ! block below, and are only used afterward in identical code
     ! blocks.
     ! Note that the gfortran -Wmaybe-uninitialized option implied by
     ! gfortran -Wall generates a false warning for this case so to avoid
     ! that could uncomment the following lines, but better yet use the
     ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
     ! code analysis to determine uninitialized variables and instead
     ! rely on run-time analysis (e.g., tainting with NaN's) to discover
     ! actual uninitialized values.

     !h2plusequil = 0._fp_kind
     !h2plusequilf = 0._fp_kind
     !h2plusequilt = 0._fp_kind

     if(ifh2plus.gt.0) then
        ! n(H2+)*avogadro/(n(H)*n(H)) = exp(h2plusequil)
        h2plusequil = h2equil + qh2plus - qh2 - bi(nions+1)*tc2 + dv(nions+2)
        if(ifnr03) then
           h2plusequilf = h2equilf + dvf(nions+2)
           h2plusequilt = h2equilt + qh2plust - qh2t + bi(nions+1)*tc2 + dvt(nions+2)
        endif
        if(ifnr13) then
           ! this derivative required for *both* indexh2 and indexh2plus
           h2plusequil_dv = 1._fp_kind
        endif
     endif
     ! solve charge equation which is quadratic in x = n(H)/avogadro,
     ! i.e., aa*x^2 + ab*x = ac
     ! this equation derived from
     ! n(H2+)/avogadro + n(H+)/avogadro + hne*rho = n_e/avogadro = cd*re,
     !   where
     !   hne*rho is the sum of non-hydrogen positive charges per
     !   unit volume divided by avogadro.
     !   rho = (2 n(H2) + 2 n(H2+) + n(H) + n(H+))/(eps(H)*avogadro)
     !       = (2 n(H2)/avogadro + 2 n(H2+)/avogadro + n(H)/avogadro + n(H+)/avogadro))/eps(H)
     !   n(H2+)/avagadro = exp(h2plusequil)*(n(H)*n(H))/(avogadro*avogadro) =  exp(h2plusequil)*x^2
     !   n(H2)/avagadro = exp(h2equil)*(n(H)*n(H))/(avogadro*avogadro) = exp(h2equil)*x^2
     !   n(H+)/avogadro = exp(h_ion_equil)*n(H)/avogadro = exp(h_ion_equil)*x
     !   n(H)/avogadro = x
     !
     !   The RHS of this equation is
     !   ac = cd*re
     !
     !   Collecting terms in x and substituting eps(H) = eps(1) implies
     !   ab = exp(h_ion_equil) + (hne/eps(1))*(exp(h_ion_equil) + 1._fp_kind)
     !   Furthermore, if you define
     !   abh = h_ion_equil + log(1._fp_kind + hne/eps(1))
     !   abg = log(hne/eps(1)), and
     !   ablog = log(exp(abh) + exp(abg)), then
     !   ab = exp(ablog)
     !   To avoid significance loss for large h_ion_equil define "d" analogs
     !   of ab, abh, and abg where exp(h_ion_equil) is divided out.
     !   abdh = log(1._fp_kind + hne/eps(1))
     !   abdg = log(hne/eps(1)) - h_ion_equil
     !   abdlog = log(exp(abdh) + exp(abdg))
     !   abd = exp(abdlog) = 1._fp_kind + (hne/eps(1))*(1._fp_kind + 1._fp_kind/exp(h_ion_equil)
     !
     !   Collecting terms in x^2 and substituting eps(H) = eps(1) implies
     !   aa = exp(h2plusequil) + (2._fp_kind*hne/eps(1))*(exp(h2plusequil) + exp(h2equil))
     !   N.B. if no molecular formation, then aa is zero, and the quadratic reduces to
     !   a linear equation in x.
     !

     ! form ln quantities first:
     abh = h_ion_equil + log(1._fp_kind + hne/eps(1))
     abdh = log(1._fp_kind + hne/eps(1))
     if(hne.gt.0._fp_kind) then
        abg = log(hne/eps(1))
        abdg = abg - h_ion_equil
     else
        ! Analysis of the hne.le.0._fp_kind case (where, e.g., ifabh
        ! is always .true. and ifabg is always .false. regardless of
        ! abg value) shows results do not depend on the following
        ! variable for this case.  But must set this variable to an
        ! arbitrary value to avoid uninitialized warning for, e.g.,
        ! out-of-order boolean logic that might occur for some
        ! compilers for Boolean calculation of ifabh or ifabg below.
        abg = 0._fp_kind

        ! Analysis of the hne.le.0._fp_kind case (where, e.g., ifabh
        ! is always .true. and ifabg is always .false. regardless of
        ! abdg value) shows results do not depend on abdg for this
        ! case.
        ! Note that the gfortran -Wmaybe-uninitialized option
        ! implied by gfortran -Wall generates a false warning for this
        ! case so to avoid that could uncomment the following line,
        ! but better yet use the gcc option -Wno-maybe-uninitialized
        ! to drop this aspect of code analysis to determine
        ! uninitialized variables and instead rely on run-time
        ! analysis (e.g., tainting with NaN's) to discover actual
        ! uninitialized values.

        !abdg = 0._fp_kind
     endif
     ! exp(-46) ~ 1.d-20
     ifabh = hne.le.0._fp_kind.or.(hne.gt.0._fp_kind.and.abg.lt.abh-46._fp_kind)
     ifabg = hne.gt.0._fp_kind.and.abh.lt.abg-46._fp_kind

     ! This variable is defined in the "elseif(abh.gt.abg) then" and
     ! "else" code blocks below and is only used afterward in
     ! identical code blocks.
     ! Note that the gfortran -Wmaybe-uninitialized option implied by
     ! gfortran -Wall generates a false warning for this case so to avoid
     ! that could uncomment the following line, but better yet use the
     ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
     ! code analysis to determine uninitialized variables and instead
     ! rely on run-time analysis (e.g., tainting with NaN's) to discover
     ! actual uninitialized values.

     !abexp = 0._fp_kind

     if(ifabh) then
        ablog = abh
        abdlog = abdh
     elseif(ifabg) then
        ablog = abg
        abdlog = abdg
     elseif(abh.gt.abg) then
        abexp = exp(abg-abh)
        ablog = abh + log(1._fp_kind + abexp)
        abdlog = abdh + log(1._fp_kind + abexp)
     else
        abexp = exp(abh-abg)
        ablog = abg + log(1._fp_kind + abexp)
        abdlog = abdg + log(1._fp_kind + abexp)
     endif
     ac = cd*re

     ! These variables are defined in the "if(ifnr03) then" code block
     ! below and are only used afterward in identical code blocks.
     ! Note that the gfortran -Wmaybe-uninitialized option implied by
     ! gfortran -Wall generates a false warning for this case so to avoid
     ! that could uncomment the following lines, but better yet use the
     ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
     ! code analysis to determine uninitialized variables and instead
     ! rely on run-time analysis (e.g., tainting with NaN's) to discover
     ! actual uninitialized values.

     !acf = 0._fp_kind
     !act = 0._fp_kind
     !ablogf = 0._fp_kind
     !ablogt = 0._fp_kind
     !abdlogf = 0._fp_kind
     !abdlogt = 0._fp_kind

     if(ifnr03) then
        ! abh = h_ion_equil + log(1._fp_kind + hne/eps(1))
        abhf = h_ion_equilf + (hnef/eps(1))/(1._fp_kind + hne/eps(1))
        abht = h_ion_equilt + (hnet/eps(1))/(1._fp_kind + hne/eps(1))
        abdhf = (hnef/eps(1))/(1._fp_kind + hne/eps(1))
        abdht = (hnet/eps(1))/(1._fp_kind + hne/eps(1))
        if(hne.gt.0._fp_kind) then
           ! abg = log(hne/eps(1))
           ! abdg = abg - h_ion_equil
           abgf = hnef/hne
           abgt = hnet/hne
           abdgf = abgf - h_ion_equilf
           abdgt = abgt - h_ion_equilt
        else
           ! Analysis of the hne.le.0._fp_kind case (where, e.g., ifabh is
           ! always .true. and ifabg is always .false. regardless of abg
           ! value) shows results do not depend on the following
           ! variables for this case.
           ! Note that the gfortran -Wmaybe-uninitialized option implied by
           ! gfortran -Wall generates a false warning for this case so to avoid
           ! that could uncomment the following lines, but better yet use the
           ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
           ! code analysis to determine uninitialized variables and instead
           ! rely on run-time analysis (e.g., tainting with NaN's) to discover
           ! actual uninitialized values.

           !abgf = 0._fp_kind
           !abgt = 0._fp_kind
           !abdgf = 0._fp_kind
           !abdgt = 0._fp_kind
        endif
        if(ifabh) then
           ! ablog = abh
           ablogf = abhf
           ablogt = abht
           abdlogf = abdhf
           abdlogt = abdht
        elseif(ifabg) then
           ! ablog = abg
           ablogf = abgf
           ablogt = abgt
           abdlogf = abdgf
           abdlogt = abdgt
        elseif(abh.gt.abg) then
           ! abexp = exp(abg-abh)
           ! ablog = abh + log(1._fp_kind + abexp)
           ablogf = abhf + (abgf-abhf)*abexp/(1._fp_kind + abexp)
           ablogt = abht + (abgt-abht)*abexp/(1._fp_kind + abexp)
           abdlogf = abdhf + (abgf-abhf)*abexp/(1._fp_kind + abexp)
           abdlogt = abdht + (abgt-abht)*abexp/(1._fp_kind + abexp)
        else
           ! abexp = exp(abh-abg)
           ! ablog = abg + log(1._fp_kind + abexp)
           ablogf = abgf + (abhf-abgf)*abexp/(1._fp_kind + abexp)
           ablogt = abgt + (abht-abgt)*abexp/(1._fp_kind + abexp)
           abdlogf = abdgf + (abhf-abgf)*abexp/(1._fp_kind + abexp)
           abdlogt = abdgt + (abht-abgt)*abexp/(1._fp_kind + abexp)
        endif
        acf = ac*ref
        act = ac*ret
     endif !if(ifnr03) then
     if(ifnr13) then
        do index = 1,max_index
           ! abh = h_ion_equil + log(1._fp_kind + hne/eps(1))
           abh_dv = (hne_dv(index)/eps(1))/(1._fp_kind + hne/eps(1))
           abdh_dv = abh_dv
           if(index.eq.indexh) abh_dv = abh_dv + h_ion_equil_dv
           if(hne.gt.0._fp_kind) then
              ! abg = log(hne/eps(1))
              abg_dv = hne_dv(index)/hne
              abdg_dv = abg_dv
              if(index.eq.indexh) abdg_dv = abdg_dv - h_ion_equil_dv
           else
              ! Analysis of the hne.le.0._fp_kind case (where, e.g., ifabh is
              ! always .true. and ifabg is always .false. regardless of abg
              ! value) shows results do not depend on the following
              ! variables for this case.
              ! Note that the gfortran -Wmaybe-uninitialized option implied by
              ! gfortran -Wall generates a false warning for this case so to avoid
              ! that could uncomment the following lines, but better yet use the
              ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
              ! code analysis to determine uninitialized variables and instead
              ! rely on run-time analysis (e.g., tainting with NaN's) to discover
              ! actual uninitialized values.

              !abg_dv = 0._fp_kind
              !abdg_dv = 0._fp_kind
           endif
           if(ifabh) then
              ! ablog = abh
              ablog_dv(index) = abh_dv
              abdlog_dv(index) = abdh_dv
           elseif(ifabg) then
              ! ablog = abg
              ablog_dv(index) = abg_dv
              abdlog_dv(index) = abdg_dv
           elseif(abh.gt.abg) then
              ! abexp = exp(abg-abh)
              ! ablog = abh + log(1._fp_kind + abexp)
              ablog_dv(index) = abh_dv + (abg_dv-abh_dv)*abexp/(1._fp_kind + abexp)
              abdlog_dv(index) = abdh_dv + (abg_dv-abh_dv)*abexp/(1._fp_kind + abexp)
           else
              ! abexp = exp(abh-abg)
              ! ablog = abg + log(1._fp_kind + abexp)
              ablog_dv(index) = abg_dv + (abh_dv-abg_dv)*abexp/(1._fp_kind + abexp)
              abdlog_dv(index) = abdg_dv + (abh_dv-abg_dv)*abexp/(1._fp_kind + abexp)
           endif
        enddo ! do index = 1,max_index
     endif !if(ifnr13) then

     ! These variables are defined in the "if(ifnr03) then" code blocks
     ! below, and are only used afterward in identical code blocks
     ! Note that the gfortran -Wmaybe-uninitialized option implied by
     ! gfortran -Wall generates a false warning for this case so to avoid
     ! that could uncomment the following lines, but better yet use the
     ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
     ! code analysis to determine uninitialized variables and instead
     ! rely on run-time analysis (e.g., tainting with NaN's) to discover
     ! actual uninitialized values.

     !xprimelogf = 0._fp_kind
     !xprimelogt = 0._fp_kind

     if(ifh2plus.gt.0.or.hne.gt.0._fp_kind) then
        !   N.B. aa, ab, and ac are all positive.  This suggests transforming
        !   aa*x^2 + ab*x = ac
        !   equation to
        !   c' x'^2 + x' = 1
        !   where
        !   x' = x/xmax <= 1
        !   xmax = ac/ab is the solution if no molecular formation (i.e. aa=0).
        !   c' = (aa/ac)*xmax^2 = aa*ac/(ab*ab)
        ! Reminders:
        ! ab = exp(ablog),
        ! aah2plus = h2plusequil + ln(1._fp_kind + 2._fp_kind*hne/eps(1))
        ! aah2 = h2equil + ln(2._fp_kind*hne/eps(1))
        ! aalog = log(exp(aah2plus) + exp(aah2)), where
        ! aa = exp(aalog),
        ! form ln quantities first:
        if(ifh2plus.gt.0) then
           aah2plus = h2plusequil + log(1._fp_kind + 2._fp_kind*hne/eps(1))
        else
           ! ifh2plus.le.0 corresponds to hne.gt.0._fp_kind (by the
           ! enclosing if).  Thus, this case implies ifaah2plus =
           ! .false. and ifaah2 = .true.  regardless of the values of
           ! aah2plus and aah2, and the subsequent calculations do not
           ! depend on the value of aah2plus.  But must set this variable to
           ! an arbitrary value to avoid uninitialized
           ! warning for, e.g., out-of-order boolean logic
           ! that might occur for some compilers for Boolean calculation
           ! of ifaah2plus or ifaah2 below.

           aah2plus = 0._fp_kind
        endif

        if(hne.gt.0._fp_kind) then
           aah2 = h2equil + log(2._fp_kind*hne/eps(1))
        else
           ! hne.le.0._fp_kind corresponds to ifh2plus.gt.0 (by the
           ! enclosing if).  Thus, this case implies ifaah2plus =
           ! .true. and ifaah2 = .false. regardless of the value of aah2, and
           ! the subsequent calculations do not
           ! depend on the value of aah2.  But must set this variable to
           ! an arbitrary value to avoid uninitialized
           ! warning for, e.g., out-of-order boolean logic
           ! that might occur for some compilers for Boolean calculation
           ! of ifaah2plus or ifaah2 below.

           aah2 = 0._fp_kind
        endif
        ! exp(-46) ~ 1.d-20
        ifaah2plus = hne.le.0._fp_kind.or.(ifh2plus.gt.0.and.hne.gt.0._fp_kind.and.aah2.lt.aah2plus-46._fp_kind)
        ifaah2 = ifh2plus.le.0.or.(ifh2plus.gt.0.and.hne.gt.0._fp_kind.and.aah2plus.lt.aah2-46._fp_kind)

        ! figure out meaning of both ifaah2plus and ifaah being .false. to
        ! figure out what needs to be true to reach third or fourth logical block below.
        ! if(ifaah2plus) then
        !  ...
        ! elseif(ifaah2) then
        !  ...
        ! elseif(aah2plus.gt.aah2) then
        !  ...
        ! else
        !  ...
        ! endif
        !
        ! .not.ifaah2plus.and..not.ifaah =
        ! .not.(ifaah2plus.or.ifaah) =
        ! .not.(hne.le.0._fp_kind.or.ifh2plus.le.0.or.(complex logical expression)) =
        ! hne.gt.0._fp_kind.and.ifh2plus.gt.0.and..not.(complex logical expression)
        ! which can only be true when aah2plus and aah2 are both set to their true
        ! values above rather than the arbitrary 0._fp_kind values above.

        ! This variable is defined in the "elseif(aah2plus.gt.aah2)
        ! then" and "else" code blocks below and is only used
        ! afterward in identical code blocks.
        ! Note that the gfortran -Wmaybe-uninitialized option implied by
        ! gfortran -Wall generates a false warning for this case so to avoid
        ! that could uncomment the following line, but better yet use the
        ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
        ! code analysis to determine uninitialized variables and instead
        ! rely on run-time analysis (e.g., tainting with NaN's) to discover
        ! actual uninitialized values.

        !aaexp = 0._fp_kind

        if(ifaah2plus) then
           aalog = aah2plus
        elseif(ifaah2) then
           aalog = aah2
        elseif(aah2plus.gt.aah2) then
           aaexp = exp(aah2-aah2plus)
           aalog = aah2plus + log(1._fp_kind + aaexp)
        else
           aaexp = exp(aah2plus-aah2)
           aalog = aah2 + log(1._fp_kind + aaexp)
        endif
        ! form cprime in ln form to start
        cprime = aalog + log(ac) - 2._fp_kind*ablog
        if(cprime.le.ln_underflow_limit) then
           ! cprime = 0 implies xprime = 1.
           xprimelog = 0._fp_kind
           if(ifnr03) then
              xprimelogf = 0._fp_kind
              xprimelogt = 0._fp_kind
           endif
           if(ifnr13) then
              xprimelog_dv(1:max_index) = 0._fp_kind
           endif
        elseif(cprime.le.92._fp_kind) then
           ! limit corresponds to about 1.d40
           cprime = exp(cprime)
           ! solve cprime*xprime^2 + xprime = 1._fp_kind
           ! for xprimelog = log(xprime) and derivatives
           ! temporary use which is immediately superseded.
           xprimelog = sqrt(1._fp_kind + 4._fp_kind*cprime)
           ! for now xprimelog actually carries xprime.
           xprimelog = 2._fp_kind/(1._fp_kind + xprimelog)
           if(ifnr03) then
              if(ifh2plus.gt.0) then
                 ! aah2plus = h2plusequil + log(1._fp_kind + 2._fp_kind*hne/eps(1))
                 aah2plusf = h2plusequilf + (2._fp_kind*hnef/eps(1))/(1._fp_kind + 2._fp_kind*hne/eps(1))
                 aah2plust = h2plusequilt + (2._fp_kind*hnet/eps(1))/(1._fp_kind + 2._fp_kind*hne/eps(1))
              else
                 ! ifh2plus.le.0 corresponds to hne.gt.0._fp_kind (by the
                 ! enclosing if).  Thus, this case implies ifaah2plus
                 ! = .false. and ifaah2 = .true. regardless of the
                 ! values of aah2plus and aah2, and the subsequent
                 ! calculations do not depend on the values of the
                 ! following variables.
                 ! Note that the gfortran -Wmaybe-uninitialized option implied by
                 ! gfortran -Wall generates a false warning for this case so to avoid
                 ! that could uncomment the following lines, but better yet use the
                 ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                 ! code analysis to determine uninitialized variables and instead
                 ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                 ! actual uninitialized values.

                 !aah2plusf = 0._fp_kind
                 !aah2plust = 0._fp_kind
              endif
              if(hne.gt.0._fp_kind) then
                 ! aah2 = h2equil + log(2._fp_kind*hne/eps(1))
                 aah2f = h2equilf + hnef/hne
                 aah2t = h2equilt + hnet/hne
              else
                 ! hne.le.0._fp_kind corresponds to ifh2plus.gt.0 (by the
                 ! enclosing if).  Thus, this case implies ifaah2plus
                 ! = .true. and ifaah2 = .false. regardless of the
                 ! values of aah2plus and aah2, and the subsequent
                 ! calculations do not depend on the values of the
                 ! following variables.
                 ! Note that the gfortran -Wmaybe-uninitialized option implied by
                 ! gfortran -Wall generates a false warning for this case so to avoid
                 ! that could uncomment the following lines, but better yet use the
                 ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                 ! code analysis to determine uninitialized variables and instead
                 ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                 ! actual uninitialized values.

                 !aah2f = 0._fp_kind
                 !aah2t = 0._fp_kind
              endif
              if(ifaah2plus) then
                 ! aalog = aah2plus
                 aalogf = aah2plusf
                 aalogt = aah2plust
              elseif(ifaah2) then
                 ! aalog = aah2
                 aalogf = aah2f
                 aalogt = aah2t
              elseif(aah2plus.gt.aah2) then
                 ! aaexp = exp(aah2-aah2plus)
                 ! aalog = aah2plus + log(1._fp_kind + aaexp)
                 aalogf = aah2plusf + (aah2f-aah2plusf)*(aaexp/(1._fp_kind+aaexp))
                 aalogt = aah2plust + (aah2t-aah2plust)*(aaexp/(1._fp_kind+aaexp))
              else
                 ! aaexp = exp(aah2plus-aah2)
                 ! aalog = aah2 + log(1._fp_kind + aaexp)
                 aalogf = aah2f + (aah2plusf-aah2f)*(aaexp/(1._fp_kind+aaexp))
                 aalogt = aah2t + (aah2plust-aah2t)*(aaexp/(1._fp_kind+aaexp))
              endif
              ! cprime = exp(aalog + log(ac) - 2._fp_kind*ablog)
              cprimef = aalogf + acf/ac - 2._fp_kind*ablogf
              cprimet = aalogt + act/ac - 2._fp_kind*ablogt
              cprimef = cprime*cprimef
              cprimet = cprime*cprimet
              ! implicit differentiation of cprime*xprime^2 + xprime = 1._fp_kind
              ! xprimelog currently carries xprime,
              ! but derivatives are of the log variable.
              xprimelogf = -cprimef*xprimelog/(2._fp_kind*cprime*xprimelog + 1._fp_kind)
              xprimelogt = -cprimet*xprimelog/(2._fp_kind*cprime*xprimelog + 1._fp_kind)
           endif !if(ifnr03) then
           if(ifnr13) then
              do index = 1,max_index
                 if(ifh2plus.gt.0) then
                    ! aah2plus = h2plusequil + log(1._fp_kind + 2._fp_kind*hne/eps(1))
                    aah2plus_dv = (2._fp_kind*hne_dv(index)/eps(1))/(1._fp_kind + 2._fp_kind*hne/eps(1))
                    if(index.eq.indexh2.or.index.eq.indexh2plus) aah2plus_dv = aah2plus_dv + h2plusequil_dv
                 else
                    ! ifh2plus.le.0 corresponds to hne.gt.0._fp_kind (by the
                    ! enclosing if).  Thus, this case implies
                    ! ifaah2plus = .false. and ifaah2 = .true.
                    ! regardless of the values of aah2plus and aah2,
                    ! and the subsequent calculations do not depend on
                    ! the value of aah2plus_dv.
                    ! Note that the gfortran -Wmaybe-uninitialized option implied by
                    ! gfortran -Wall generates a false warning for this case so to avoid
                    ! that could uncomment the following line, but better yet use the
                    ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                    ! code analysis to determine uninitialized variables and instead
                    ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                    ! actual uninitialized values.

                    !aah2plus_dv = 0._fp_kind
                 endif
                 if(hne.gt.0._fp_kind) then
                    ! aah2 = h2equil + log(2._fp_kind*hne/eps(1))
                    aah2_dv = hne_dv(index)/hne
                    if(index.eq.indexh2) aah2_dv = aah2_dv + h2equil_dv
                 else
                    ! hne.le.0._fp_kind corresponds to ifh2plus.gt.0 (by the
                    ! enclosing if).  Thus, this case implies
                    ! ifaah2plus = .true. and ifaah2 = .false.
                    ! regardless of the values of aah2plus and aah2,
                    ! and the subsequent calculations do not depend on
                    ! the value of aah2_dv.
                    ! Note that the gfortran -Wmaybe-uninitialized option implied by
                    ! gfortran -Wall generates a false warning for this case so to avoid
                    ! that could uncomment the following line, but better yet use the
                    ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                    ! code analysis to determine uninitialized variables and instead
                    ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                    ! actual uninitialized values.

                    !aah2_dv = 0._fp_kind
                 endif
                 if(ifaah2plus) then
                    ! aalog = aah2plus
                    aalog_dv = aah2plus_dv
                 elseif(ifaah2) then
                    ! aalog = aah2
                    aalog_dv = aah2_dv
                 elseif(aah2plus.gt.aah2) then
                    ! aaexp = exp(aah2-aah2plus)
                    ! aalog = aah2plus + log(1._fp_kind + aaexp)
                    aalog_dv = aah2plus_dv + (aah2_dv-aah2plus_dv)*(aaexp/(1._fp_kind+aaexp))
                 else
                    ! aaexp = exp(aah2plus-aah2)
                    ! aalog = aah2 + log(1._fp_kind + aaexp)
                    aalog_dv = aah2_dv + (aah2plus_dv-aah2_dv)*(aaexp/(1._fp_kind+aaexp))
                 endif
                 ! cprime = exp(aalog + log(ac) - 2._fp_kind*ablog)
                 cprime_dv = aalog_dv - 2._fp_kind*ablog_dv(index)
                 cprime_dv = cprime*cprime_dv
                 ! implicit differentiation of cprime*xprime^2 + xprime = 1._fp_kind
                 ! xprimelog currently carries xprime,
                 ! but derivatives are of the log variable.
                 xprimelog_dv(index) = -cprime_dv*xprimelog/(2._fp_kind*cprime*xprimelog + 1._fp_kind)
              enddo !do index = 1,max_index
           endif
           xprimelog = log(xprimelog)
        else
           ! cprime is greater than about 1.d40.
           ! solve cprime*xprime^2 = 1, but variable cprime
           ! actually carries log of that quantity so
           ! really solve: cprime + 2._fp_kind*xprimelog = 0
           xprimelog = -0.5_fp_kind*cprime
           if(ifnr03) then
              if(ifh2plus.gt.0) then
                 ! aah2plus = h2plusequil + log(1._fp_kind + 2._fp_kind*hne/eps(1))
                 aah2plusf = h2plusequilf + (2._fp_kind*hnef/eps(1))/(1._fp_kind + 2._fp_kind*hne/eps(1))
                 aah2plust = h2plusequilt + (2._fp_kind*hnet/eps(1))/(1._fp_kind + 2._fp_kind*hne/eps(1))
              else
                 ! ifh2plus.le.0 corresponds to hne.gt.0._fp_kind (by the
                 ! enclosing if).  Thus, this case implies ifaah2plus
                 ! = .false. and ifaah2 = .true.  regardless of the
                 ! values of aah2plus and aah2, and the subsequent
                 ! calculations do not depend on the values of the
                 ! following variables.
                 ! Note that the gfortran -Wmaybe-uninitialized option implied by
                 ! gfortran -Wall generates a false warning for this case so to avoid
                 ! that could uncomment the following lines, but better yet use the
                 ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                 ! code analysis to determine uninitialized variables and instead
                 ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                 ! actual uninitialized values.

                 !aah2plusf = 0._fp_kind
                 !aah2plust = 0._fp_kind
              endif
              if(hne.gt.0._fp_kind) then
                 ! aah2 = h2equil + log(2._fp_kind*hne/eps(1))
                 aah2f = h2equilf + hnef/hne
                 aah2t = h2equilt + hnet/hne
              else
                 ! hne.le.0._fp_kind corresponds to ifh2plus.gt.0 (by the
                 ! enclosing if).  Thus, this case implies ifaah2plus
                 ! = .true. and ifaah2 = .false.  regardless of the
                 ! values of aah2plus and aah2, and the subsequent
                 ! calculations do not depend on the values of the
                 ! following variables.
                 ! Note that the gfortran -Wmaybe-uninitialized option implied by
                 ! gfortran -Wall generates a false warning for this case so to avoid
                 ! that could uncomment the following lines, but better yet use the
                 ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                 ! code analysis to determine uninitialized variables and instead
                 ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                 ! actual uninitialized values.

                 !aah2f = 0._fp_kind
                 !aah2t = 0._fp_kind
              endif
              if(ifaah2plus) then
                 ! aalog = aah2plus
                 aalogf = aah2plusf
                 aalogt = aah2plust
              elseif(ifaah2) then
                 ! aalog = aah2
                 aalogf = aah2f
                 aalogt = aah2t
              elseif(aah2plus.gt.aah2) then
                 ! aaexp = exp(aah2-aah2plus)
                 ! aalog = aah2plus + log(1._fp_kind + aaexp)
                 aalogf = aah2plusf + (aah2f-aah2plusf)*(aaexp/(1._fp_kind+aaexp))
                 aalogt = aah2plust + (aah2t-aah2plust)*(aaexp/(1._fp_kind+aaexp))
              else
                 ! aaexp = exp(aah2plus-aah2)
                 ! aalog = aah2 + log(1._fp_kind + aaexp)
                 aalogf = aah2f + (aah2plusf-aah2f)*(aaexp/(1._fp_kind+aaexp))
                 aalogt = aah2t + (aah2plust-aah2t)*(aaexp/(1._fp_kind+aaexp))
              endif
              ! cprime = aalog + log(ac) - 2._fp_kind*ablog
              cprimef = aalogf + acf/ac - 2._fp_kind*ablogf
              cprimet = aalogt + act/ac - 2._fp_kind*ablogt
              xprimelogf = -0.5_fp_kind*cprimef
              xprimelogt = -0.5_fp_kind*cprimet
           endif !if(ifnr03) then
           if(ifnr13) then
              do index = 1,max_index
                 if(ifh2plus.gt.0) then
                    ! aah2plus = h2plusequil + log(1._fp_kind + 2._fp_kind*hne/eps(1))
                    aah2plus_dv = (2._fp_kind*hne_dv(index)/eps(1))/(1._fp_kind + 2._fp_kind*hne/eps(1))
                    if(index.eq.indexh2.or.index.eq.indexh2plus) aah2plus_dv = aah2plus_dv + h2plusequil_dv
                 else
                    ! ifh2plus.le.0 corresponds to hne.gt.0._fp_kind (by the
                    ! enclosing if).  Thus, this case implies
                    ! ifaah2plus = .false. and ifaah2 = .true.
                    ! regardless of the values of aah2plus and aah2,
                    ! and the subsequent calculations do not depend on
                    ! the value of aah2plus_dv.
                    ! Note that the gfortran -Wmaybe-uninitialized option implied by
                    ! gfortran -Wall generates a false warning for this case so to avoid
                    ! that could uncomment the following line, but better yet use the
                    ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                    ! code analysis to determine uninitialized variables and instead
                    ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                    ! actual uninitialized values.

                    !aah2plus_dv = 0._fp_kind
                 endif
                 if(hne.gt.0._fp_kind) then
                    ! aah2 = h2equil + log(2._fp_kind*hne/eps(1))
                    aah2_dv = hne_dv(index)/hne
                    if(index.eq.indexh2) aah2_dv = aah2_dv + h2equil_dv
                 else
                    ! hne.le.0._fp_kind corresponds to ifh2plus.gt.0 (by the
                    ! enclosing if).  Thus, this case implies
                    ! ifaah2plus = .true. and ifaah2 = .false.
                    ! regardless of the values of aah2plus and aah2,
                    ! and the subsequent calculations do not depend on
                    ! the value of aah2_dv.
                    ! Note that the gfortran -Wmaybe-uninitialized option implied by
                    ! gfortran -Wall generates a false warning for this case so to avoid
                    ! that could uncomment the following line, but better yet use the
                    ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
                    ! code analysis to determine uninitialized variables and instead
                    ! rely on run-time analysis (e.g., tainting with NaN's) to discover
                    ! actual uninitialized values.

                    !aah2_dv = 0._fp_kind
                 endif
                 if(ifaah2plus) then
                    ! aalog = aah2plus
                    aalog_dv = aah2plus_dv
                 elseif(ifaah2) then
                    ! aalog = aah2
                    aalog_dv = aah2_dv
                 elseif(aah2plus.gt.aah2) then
                    ! aaexp = exp(aah2-aah2plus)
                    ! aalog = aah2plus + log(1._fp_kind + aaexp)
                    aalog_dv = aah2plus_dv + (aah2_dv-aah2plus_dv)*(aaexp/(1._fp_kind+aaexp))
                 else
                    ! aaexp = exp(aah2plus-aah2)
                    ! aalog = aah2 + log(1._fp_kind + aaexp)
                    aalog_dv = aah2_dv + (aah2plus_dv-aah2_dv)*(aaexp/(1._fp_kind+aaexp))
                 endif
                 ! cprime = aalog + log(ac) - 2._fp_kind*ablog
                 cprime_dv = aalog_dv - 2._fp_kind*ablog_dv(index)
                 xprimelog_dv(index) = -0.5_fp_kind*cprime_dv
              enddo !do index = 1,max_index
           endif
        endif
     else !if(ifh2plus.gt.0.or.hne.gt.0._fp_kind) then
        ! cprime = 0 implies xprime = 1.
        xprimelog = 0._fp_kind
        if(ifnr03) then
           xprimelogf = 0._fp_kind
           xprimelogt = 0._fp_kind
        endif
        if(ifnr13) then
           xprimelog_dv(1:max_index) = 0._fp_kind
        endif
     endif !if(ifh2plus.gt.0.or.hne.gt.0._fp_kind) then

     ! At this point have calculated xprimelog = log(x/xmax) and derivatives.
     ! transform to nh_neutral  = log(xprime*xmax) = log(xprime*ac/ab) where
     ! nh_neutral contains log(n(neutral monatomic hydrogen)/avogadro).
     nh_neutral = xprimelog + log(ac) - ablog

     ! These variables are defined in the "if(ifnr03) then" code
     ! block below, and are only used afterward in identical code
     ! blocks.
     ! Note that the gfortran -Wmaybe-uninitialized option implied by
     ! gfortran -Wall generates a false warning for this case so to avoid
     ! that could uncomment the following lines, but better yet use the
     ! gfortan option -Wno-maybe-uninitialized to drop this aspect of
     ! code analysis to determine uninitialized variables and instead
     ! rely on run-time analysis (e.g., tainting with NaN's) to discover
     ! actual uninitialized values.

     !nh_neutralf = 0._fp_kind
     !nh_neutralt = 0._fp_kind

     if(ifnr03) then
        nh_neutralf = xprimelogf + acf/ac - ablogf
        nh_neutralt = xprimelogt + act/ac - ablogt
     endif
     if(ifnr13) then
        nh_neutral_dv(1:max_index) = xprimelog_dv(1:max_index) - ablog_dv(1:max_index)
     endif

     ! nh2 contains log(n(H2)/avogadro)
     nh2 = 2._fp_kind*nh_neutral + h2equil
     if(ifnr03) then
        nh2f = 2._fp_kind*nh_neutralf + h2equilf
        nh2t = 2._fp_kind*nh_neutralt + h2equilt
     endif
     if(ifnr13) then
        nh2_dv(1:max_index) = 2._fp_kind*nh_neutral_dv(1:max_index)
        if(indexh2.le.max_index) nh2_dv(indexh2) = nh2_dv(indexh2) + h2equil_dv
     endif
     if(ifh2plus.gt.0) then
        ! nh2plus contains log(n(H2+)/avogadro)
        nh2plus = 2._fp_kind*nh_neutral + h2plusequil
        if(ifnr03) then
           nh2plusf = 2._fp_kind*nh_neutralf + h2plusequilf
           nh2plust = 2._fp_kind*nh_neutralt + h2plusequilt
        endif
        if(ifnr13) then
           nh2plus_dv(1:max_index) = 2._fp_kind*nh_neutral_dv(1:max_index)
           if(indexh2.le.max_index) nh2plus_dv(indexh2) = nh2plus_dv(indexh2) + h2plusequil_dv
           if(indexh2plus.le.max_index) nh2plus_dv(indexh2plus) = nh2plus_dv(indexh2plus) + h2plusequil_dv
        endif
     endif
     ! nh_ion contains log(n(H+)/avogadro)
     if(h_ion_equil.gt.10._fp_kind) then
        ! large pressure ionization can force h_ion_equil to large values.
        ! avoid significance loss for this case by using the abd form
        ! of variables.
        ! nh_neutral = xprimelog + log(ac) - ablog = xprimelog + log(ac) - (abdlog + h_ion_equil)
        ! nh_ion = nh_neutral + h_ion_equil
        nh_ion = xprimelog + log(ac) - abdlog
        if(ifnr03) then
           nh_ionf = xprimelogf + acf/ac - abdlogf
           nh_iont = xprimelogt + act/ac - abdlogt
        endif
        if(ifnr13) then
           nh_ion_dv(1:max_index) = xprimelog_dv(1:max_index) - abdlog_dv(1:max_index)
        endif
     else
        nh_ion = nh_neutral + h_ion_equil
        if(ifnr03) then
           nh_ionf = nh_neutralf + h_ion_equilf
           nh_iont = nh_neutralt + h_ion_equilt
        endif
        if(ifnr13) then
           nh_ion_dv(1:max_index) = nh_neutral_dv(1:max_index)
           if(indexh.le.max_index) nh_ion_dv(indexh) = nh_ion_dv(indexh) + h_ion_equil_dv
        endif
     endif
     ! N.B. nh2, nh2plus, nh_neutral, and nh_ion are all in
     ! log(n/avogadro) form where n is the number density of H2, H2+,
     ! neutral monatomic H, and H+.  Therefore, the general
     ! expression for rl = log(rho) is
     ! rl = log(2._fp_kind*exp(nh2) + 2._fp_kind*exp(nh2plus) + exp(nh_ion) + exp(nh_neutral)) - log(eps(1))

     ! Find whether any of nh2, nh2plus, nh_neutral, or nh_ion (all in ln(n/avogadro) form)
     ! dominant the above equation so can calculate rl without fear of
     ! over/underflow where "dominant" means one of nh2, nh2plus,
     ! nh_neutral, or nh_ion is lnerrcrit larger than any of the rest
     ! so that the above log factor can be approximated by the maximum
     ! component with error less than ln(1 + exp(-lnerrcrit)) ~
     ! exp(-lnerrcrit).
     maxh = nh2
     index_maxh = 1
     ! Initialize second_maxh to a small value that insures the index_maxhh = 0 logic
     ! below is NOT executed for the case where nh2 dominates, i.e.,
     ! all of nh2plus, nh_neutral, and nh_ion, are less than this small value.
     second_maxh = maxh - lnerrcrit - 1._fp_kind

     if(ifh2plus.gt.0) then
        if(nh2plus.gt.maxh) then
           index_maxh = 2
           second_maxh = maxh
           maxh = nh2plus
        else
           second_maxh = max(second_maxh, nh2plus)
        endif
     endif
     if(nh_ion.gt.maxh) then
        index_maxh = 3
        second_maxh = maxh
        maxh = nh_ion
     else
        second_maxh = max(second_maxh, nh_ion)
     endif
     if(nh_neutral.gt.maxh) then
        index_maxh = 4
        second_maxh = maxh
        maxh = nh_neutral
     else
        second_maxh = max(second_maxh, nh_neutral)
     endif

     ! index_maxh = 0 implies maxh is not sufficiently larger than
     ! second_maxh, i.e, no single hydrogen species dominates the
     ! hydrogen abundance constraint equation.
     if((maxh - lnerrcrit).lt.second_maxh) index_maxh = 0

     ! *if* a single hydrogen species dominates the hydrogen abundance
     ! equation (i.e., index_maxh.ne.0) then can calculate rl and
     ! derivatives without transforming this dominate term to number
     ! density form so in these cases no fear of overflows/underflows.

     if(index_maxh.eq.1) then
        ! H2 dominates
        rl = nh2 - log(eps(1)/2._fp_kind)
        if(ifnr03) then
           rf = nh2f
           rt = nh2t
        endif
        if(ifnr13) then
           r_dv(1:max_index) = nh2_dv(1:max_index)
        endif
     elseif(index_maxh.eq.2) then
        ! H2+ dominates
        if(ifh2plus.le.0) error stop 'eos_calc: internal ifh2plus error'
        rl = nh2plus - log(eps(1)/2._fp_kind)
        if(ifnr03) then
           rf = nh2plusf
           rt = nh2plust
        endif
        if(ifnr13) then
           r_dv(1:max_index) = nh2plus_dv(1:max_index)
        endif
     elseif(index_maxh.eq.3) then
        ! H+ dominates
        rl = nh_ion - log(eps(1))
        if(ifnr03) then
           rf = nh_ionf
           rt = nh_iont
        endif
        if(ifnr13) then
           r_dv(1:max_index) = nh_ion_dv(1:max_index)
        endif
     elseif(index_maxh.eq.4) then
        ! Neutral monatomic H dominates
        rl = nh_neutral - log(eps(1))
        if(ifnr03) then
           rf = nh_neutralf
           rt = nh_neutralt
        endif
        if(ifnr13) then
           r_dv(1:max_index) = nh_neutral_dv(1:max_index)
        endif
     endif

     ! Transform nh2, nh2plus, nh_neutral, and nh_ion from ln (n/avogadro) form to n form

     if(nh2 + ln_avogadro.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "n(H2) is too close to overflowing"
        info = info_offset_eos_calc + 11
        return
     elseif(nh2 + ln_avogadro.gt.ln_underflow_limit) then
        nh2 = exp(nh2 + ln_avogadro)
        if(ifnr03) then
           nh2f = nh2*nh2f
           nh2t = nh2*nh2t
        endif
        if(ifnr13) then
           nh2_dv(1:max_index) = nh2*nh2_dv(1:max_index)
        endif
     else
        nh2 = 0._fp_kind
        if(ifnr03) then
           nh2f = 0._fp_kind
           nh2t = 0._fp_kind
        endif
        if(ifnr13) then
           nh2_dv(1:max_index) = 0._fp_kind
        endif
     endif

     if(ifh2plus.gt.0) then
        if(nh2plus + ln_avogadro.gt.ln_overflow_limit) then
           if(verbosity.ge.1) write(stderr,*) "n(H2+) is too close to overflowing"
           info = info_offset_eos_calc + 11
           return
        elseif(nh2plus + ln_avogadro.gt.ln_underflow_limit) then
           nh2plus = exp(nh2plus + ln_avogadro)
           if(ifnr03) then
              nh2plusf = nh2plus*nh2plusf
              nh2plust = nh2plus*nh2plust
           endif
           if(ifnr13) then
              nh2plus_dv(1:max_index) = nh2plus*nh2plus_dv(1:max_index)
           endif
        else
           ! nh2plus would underflow.  Zero its value and all relevant derivatives.
           nh2plus = 0._fp_kind
           if(ifnr03) then
              nh2plusf = 0._fp_kind
              nh2plust = 0._fp_kind
           endif
           if(ifnr13) then
              nh2plus_dv(1:max_index) = 0._fp_kind
           endif
        endif
     else
        ! No H2+ included in free-energy model. Zero its value and all relevant derivatives.
        nh2plus = 0._fp_kind
        if(ifnr03) then
           nh2plusf = 0._fp_kind
           nh2plust = 0._fp_kind
        endif
        if(ifnr13) then
           nh2plus_dv(1:max_index) = 0._fp_kind
        endif
     endif

     ! preserve nh_neutral (which currently contains log(n(neutral
     ! monatomic H)/avogadro)) in log_h_neutral variable for future
     ! reference.
     log_h_neutral = nh_neutral

     if(nh_neutral + ln_avogadro.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "n(neutral monatomic H) is too close to overflowing"
        info = info_offset_eos_calc + 11
        return
     elseif(nh_neutral + ln_avogadro.gt.ln_underflow_limit) then
        nh_neutral = exp(nh_neutral + ln_avogadro)
        if(ifnr03) then
           nh_neutralf = nh_neutral*nh_neutralf
           nh_neutralt = nh_neutral*nh_neutralt
        endif
        if(ifnr13) then
           nh_neutral_dv(1:max_index) = nh_neutral*nh_neutral_dv(1:max_index)
        endif
     else
        nh_neutral = 0._fp_kind
        if(ifnr03) then
           nh_neutralf = 0._fp_kind
           nh_neutralt = 0._fp_kind
        endif
        if(ifnr13) then
           nh_neutral_dv(1:max_index) = 0._fp_kind
        endif
     endif

     if(nh_ion + ln_avogadro.gt.ln_overflow_limit) then
        if(verbosity.ge.1) write(stderr,*) "n(H+) is too close to overflowing"
        info = info_offset_eos_calc + 11
        return
     elseif(nh_ion + ln_avogadro.gt.ln_underflow_limit) then
        nh_ion = exp(nh_ion + ln_avogadro)
        if(ifnr03) then
           nh_ionf = nh_ion*nh_ionf
           nh_iont = nh_ion*nh_iont
        endif
        if(ifnr13) then
           nh_ion_dv(1:max_index) = nh_ion*nh_ion_dv(1:max_index)
        endif
     else
        nh_ion = 0._fp_kind
        if(ifnr03) then
           nh_ionf = 0._fp_kind
           nh_iont = 0._fp_kind
        endif
        if(ifnr13) then
           nh_ion_dv(1:max_index) = 0._fp_kind
        endif
     endif

     if(index_maxh.eq.0) then
        ! No one hydrogen species dominates the hydrogen abundance
        ! equation so use general expression for rl and derivatives.

        if(ifh2plus.gt.0) then
           ! Include the H2+ terms
           rl = log(2._fp_kind*nh2 + 2._fp_kind*nh2plus + nh_neutral + nh_ion) - ln_avogadro - log(eps(1))
           if(rl.lt.ln_underflow_limit) then
              if(verbosity.ge.1) write(stderr,*) "eos_calc: (24) calculated mass density too close to underflow limit"
              info = info_offset_eos_calc + 24
              return
           elseif(rl.gt.ln_overflow_limit) then
              if(verbosity.ge.1) write(stderr,*) "eos_calc: (34) calculated mass density too close to overflow limit"
              info = info_offset_eos_calc + 34
              return
           endif
           rho = exp(rl)*rho_perturb_factor
           if(ifnr03) then
              rf = (2._fp_kind*nh2f + 2._fp_kind*nh2plusf + nh_neutralf + nh_ionf)/(rho*avogadro*eps(1))
              rt = (2._fp_kind*nh2t + 2._fp_kind*nh2plust + nh_neutralt + nh_iont)/(rho*avogadro*eps(1))
           endif
           if(ifnr13) then
              r_dv(1:max_index) = (&
                   2._fp_kind*nh2_dv(1:max_index) +&
                   2._fp_kind*nh2plus_dv(1:max_index) +&
                   nh_neutral_dv(1:max_index) +&
                   nh_ion_dv(1:max_index)&
                   )/(rho*avogadro*eps(1))
           endif
        else
           ! Exclude the H2+ terms (and guard against any undefined H2+ derivatives in this case).
           rl = log(2._fp_kind*nh2 + nh_neutral + nh_ion) - ln_avogadro - log(eps(1))
           if(rl.lt.ln_underflow_limit) then
              if(verbosity.ge.1) write(stderr,*) "eos_calc: (24) calculated mass density too close to underflow limit"
              info = info_offset_eos_calc + 24
              return
           elseif(rl.gt.ln_overflow_limit) then
              if(verbosity.ge.1) write(stderr,*) "eos_calc: (34) calculated mass density too close to overflow limit"
              info = info_offset_eos_calc + 34
              return
           endif
           rho = exp(rl)*rho_perturb_factor
           if(ifnr03) then
              rf = (2._fp_kind*nh2f + nh_neutralf + nh_ionf)/(rho*avogadro*eps(1))
              rt = (2._fp_kind*nh2t + nh_neutralt + nh_iont)/(rho*avogadro*eps(1))
           endif
           if(ifnr13) then
              r_dv(1:max_index) = (&
                   2._fp_kind*nh2_dv(1:max_index) +&
                   nh_neutral_dv(1:max_index) +&
                   nh_ion_dv(1:max_index)&
                   )/(rho*avogadro*eps(1))
           endif
        endif
     else
        ! rl calculated above for each case where one of the hydrogen species dominates the
        ! hydrogen abundance equation.
        if(rl.lt.ln_underflow_limit) then
           if(verbosity.ge.1) write(stderr,*) "eos_calc: (25) calculated mass density too close to underflow limit"
           info = info_offset_eos_calc + 25
           return
        elseif(rl.gt.ln_overflow_limit) then
           if(verbosity.ge.1) write(stderr,*) "eos_calc: (35) calculated mass density too close to overflow limit"
           info = info_offset_eos_calc + 35
           return
        endif
        rho = exp(rl)*rho_perturb_factor
     endif

     ! Calculate (without fear of over or underflows) the log nu form
     ! of the neutral monatomic species of H from ln(n/avogadro) form
     ! which was previously preserved in the log_h_neutral variable.
     log_h_neutral = log_h_neutral - rl

     ! Calculate nu form of hydrogen species (saved in h* variables)
     ! from n form (saved in nh* variables).  Assuming the range of (cgs) rho values
     ! of interest for stellar interiors work is roughly between 10^{-10} to 10^{6}, then
     ! rho*avogadro (cgs) is roughly in the range from 10^{14} through 10^{30}.
     ! Thus dividing by this factor may induce underflows but definitely not overflows
     ! (which are already checked in the case of the n form).

     if(nh2.gt.underflow_limit*(rho*avogadro)) then
        h2 = nh2/(rho*avogadro)
        if(ifnr03) then
           h2f = nh2f/(rho*avogadro) - h2*rf
           h2t = nh2t/(rho*avogadro) - h2*rt
        endif
        if(ifnr13) then
           h2_dv(1:max_index) = nh2_dv(1:max_index)/(rho*avogadro) - h2*r_dv(1:max_index)
        endif
     else
        h2 = 0._fp_kind
        if(ifnr03) then
           h2f = 0._fp_kind
           h2t = 0._fp_kind
        endif
        if(ifnr13) then
           h2_dv(1:max_index) = 0._fp_kind
        endif
     endif

     if(ifh2plus.gt.0.and.nh2plus.gt.underflow_limit*(rho*avogadro)) then
        h2plus = nh2plus/(rho*avogadro)
        if(ifnr03) then
           h2plusf = nh2plusf/(rho*avogadro) - h2plus*rf
           h2plust = nh2plust/(rho*avogadro) - h2plus*rt
        endif
        if(ifnr13) then
           h2plus_dv(1:max_index) = nh2plus_dv(1:max_index)/(rho*avogadro) - h2plus*r_dv(1:max_index)
        endif
     else
        h2plus = 0._fp_kind
        if(ifnr03) then
           h2plusf = 0._fp_kind
           h2plust = 0._fp_kind
        endif
        if(ifnr13) then
           h2plus_dv(1:max_index) = 0._fp_kind
        endif
     endif

     if(nh_neutral.gt.underflow_limit*(rho*avogadro)) then
        h_neutral = nh_neutral/(rho*avogadro)
        if(ifnr03) then
           h_neutralf = nh_neutralf/(rho*avogadro) - h_neutral*rf
           h_neutralt = nh_neutralt/(rho*avogadro) - h_neutral*rt
        endif
        if(ifnr13) then
           h_neutral_dv(1:max_index) = nh_neutral_dv(1:max_index)/(rho*avogadro) - h_neutral*r_dv(1:max_index)
        endif
     else
        h_neutral = 0._fp_kind
        if(ifnr03) then
           h_neutralf = 0._fp_kind
           h_neutralt = 0._fp_kind
        endif
        if(ifnr13) then
           h_neutral_dv(1:max_index) = 0._fp_kind
        endif
     endif

     if(nh_ion.gt.underflow_limit*(rho*avogadro)) then
        h_ion = nh_ion/(rho*avogadro)
        if(ifnr03) then
           h_ionf = nh_ionf/(rho*avogadro) - h_ion*rf
           h_iont = nh_iont/(rho*avogadro) - h_ion*rt
        endif
        if(ifnr13) then
           h_ion_dv(1:max_index) = nh_ion_dv(1:max_index)/(rho*avogadro) - h_ion*r_dv(1:max_index)
        endif
     else
        h_ion = 0._fp_kind
        if(ifnr03) then
           h_ionf = 0._fp_kind
           h_iont = 0._fp_kind
        endif
        if(ifnr13) then
           h_ion_dv(1:max_index) = 0._fp_kind
        endif
     endif

     nuvar(1,1) = h_neutral
     nuvar(2,1) = h_ion
     ! store h2 and h2plus in the 3rd and 4th index
     nuvar(3,1) = h2
     nuvar(4,1) = h2plus
     ne = hne + h_ion + h2plus
     if(ne.lt.underflow_limit) then
        if(verbosity.ge.1) write(stderr,*)&
             "eos_calc: (4) calculated electron number density (in nu form) is too close to underflowing"
        info = info_offset_eos_calc + 4
        return
     endif

     if(ifnr03) then
        nef = hnef + h_ionf + h2plusf
        net = hnet + h_iont + h2plust
     endif

     en = 1._fp_kind/ne
     if(eps(2).eq.0._fp_kind) then
        he_ion = 0._fp_kind
        he_ion2 = 0._fp_kind
        if(ifpi.eq.2.or.if_mc.eq.1) then
           if(ifnr03) then
              he_ionf = 0._fp_kind
              he_ion2f = 0._fp_kind
              he_ionf = 0._fp_kind
              he_ion2t = 0._fp_kind
           endif
        endif
     endif
     ! these expressions derived from vdb notes. they use the equations:
     ! 2 nu(H2+) + 2 nu(H2) + nu(H) + nu(H+) = eps(1),
     ! where the fortran variables h2plus, h2, h_neutral, and h_ion
     ! carry n(H2+), n(H2), n(H), n(H+) or the n/(avogadro*rho)
     ! equivalents when ifnuform is true.  It follows that
     ! alpha nu(H+)/(A(H+)^3/2 Q(H+)) =
     !   alpha nu(H)/(A(H)^3/2 Q(H)) *
     !   exp(-tc2*Hion + dv(1))
     ! alpha nu(H2)/(A(H2)^3/2 Q(H2)) =
     !   rho T^{-3/2} [alpha nu(H)/(A(H)^3/2 Q(H))]^2 *
     !   exp(tc2*hdiss + dv(nions+1))
     ! alpha nu(H2+)/(A(H2+)^3/2 Q(H2+)) =
     !   alpha nu(H2)/(A(H2)^3/2 Q(H2)) *
     !   exp(-tc2*H2ion + dv(nions+2))
     !   where alpha = avogadro*(2 pi k/(avogadro*h^2)^{-3/2}
     !   (see commentary in the constants module), and
     !   Q(i) exp(-tc2*E(i)) A(i)^{3/2}/alpha is the
     !   total partition per unit volume divided by avogadro.
     !   Furthermore, we have taken advantage of the fact that
     !   the individual rho, T dependence of n(H2+) and n(H2)
     !   partially cancels the overall rho, T dependence, see
     !   expression for s in free_eos_detailed.f.
     !   note that h2equilt0 = qh2t - 1.5_fp_kind - h2diss*tc2.  Thus, must
     !   subtract 1 from this quantity in sion so that can correct by
     !   full_sum0 in free_eos.f for the 5/2 term.  The equivalent ln T
     !   and ln rho terms are taken care of by the above equilibrium
     !   constant relations.
     sharg = bi(1)*tc2 - dv(1)
     if(ifh2plus.eq.0) then
        ! For this case, qh2plus and derivatives, and the nions+2
        ! indices of dv and derivatives are uninitialized so to avoid
        ! valgrind uninitialized tainting of the following results
        ! (even though h2plus, etc.  are zero) must drop the h2plus
        ! components of these expressions.
        sion = sion + constant_sh  - eps(1)*log_h_neutral + h_ion*sharg + h2*(h2equilt0-dv(nions+1)-1._fp_kind)
        ! correct ideal internal energy cm^-1 per unit mass
        ! per avogadro for hydrogen nu values.
        uion = uion + h_ion*bi(1) + h2*(-h2diss + qh2t/tc2)
        if(ifnr03) then
           sionf = sionf + h_ionf*sharg + h2f*(h2equilt0 - dv(nions+1))
           siont = siont + h_iont*sharg + h2t*(h2equilt0 - dv(nions+1)) + h2*(qh2t + qh2tt)
        endif
     else
        sion = sion + constant_sh  - eps(1)*log_h_neutral + h_ion*sharg + h2*(h2equilt0-dv(nions+1)-1._fp_kind) +&
             h2plus*(h2equilt0 - 1._fp_kind + qh2plust - qh2t + bi(nions+1)*tc2 - dv(nions+1) - dv(nions+2))
        ! correct ideal internal energy cm^-1 per unit mass
        ! per avogadro for hydrogen nu values.
        uion = uion + h_ion*bi(1) + h2*(-h2diss + qh2t/tc2) + h2plus*(bi(nions+1) - h2diss + qh2plust/tc2)
        if(ifnr03) then
           sionf = sionf + h_ionf*sharg + h2f*(h2equilt0 - dv(nions+1)) +&
                h2plusf*(h2equilt0 + qh2plust - qh2t + bi(nions+1)*tc2 - dv(nions+1) - dv(nions+2))
           siont = siont + h_iont*sharg + h2t*(h2equilt0 - dv(nions+1)) + h2*(qh2t + qh2tt) +&
                h2plust*(h2equilt0 + qh2plust - qh2t + bi(nions+1)*tc2 - dv(nions+1) - dv(nions+2)) +&
                h2plus*(qh2plust+qh2plustt)
        endif
     endif
     ! add in ideal free-energy terms due to all hydrogen species
     ! from first principles (MHD II, F1 + F2 term transformed to
     ! free-energy per unit mass and ignoring constant
     ! [1 + 1.5_fp_kind ln T - ln rho] term).
     if(h2plus.gt.0._fp_kind) then
        fion = fion +&
             h2plus*(log(h2plus) + tc2*(bi(nions+1)-h2diss) - qh2plus - logqtl_const_h2plus)
        if(ifnr03)&
             fionf = fionf +&
             h2plusf*(log(h2plus) + tc2*(bi(nions+1)-h2diss) - qh2plus - logqtl_const_h2plus + 1._fp_kind)
        if(ifnr13) then
           fion_dv(1:max_index) = fion_dv(1:max_index) +&
                h2plus_dv(1:max_index)*(log(h2plus) + tc2*(bi(nions+1)-h2diss) - qh2plus - logqtl_const_h2plus + 1._fp_kind)
        endif
     endif
     if(h2.gt.0._fp_kind) then
        fion = fion + h2*(log(h2) + tc2*(-h2diss) - qh2 - logqtl_const_h2)
        if(ifnr03)&
             fionf = fionf + h2f*(log(h2) + tc2*(-h2diss) - qh2 - logqtl_const_h2 + 1._fp_kind)
        if(ifnr13) then
           fion_dv(1:max_index) = fion_dv(1:max_index) +&
                h2_dv(1:max_index)*(log(h2) + tc2*(-h2diss) - qh2 - logqtl_const_h2 + 1._fp_kind)
        endif
     endif
     if(h_neutral.gt.0._fp_kind) then
        fion = fion + h_neutral*(log_h_neutral - logqtl_const_h)
        if(ifnr03)&
             fionf = fionf +&
             h_neutralf*(log_h_neutral - logqtl_const_h + 1._fp_kind)
        if(ifnr13) then
           fion_dv(1:max_index) = fion_dv(1:max_index) +&
                h_neutral_dv(1:max_index)*(log_h_neutral - logqtl_const_h + 1._fp_kind)
        endif
     endif
     if(h_ion.gt.0._fp_kind) then
        fion = fion + h_ion*(log(h_ion) + tc2*bi(1) - logqtl_const_hplus)
        if(ifnr03) fionf = fionf + h_ionf*(log(h_ion) + tc2*bi(1) - logqtl_const_hplus + 1._fp_kind)
        if(ifnr13) then
           fion_dv(1:max_index) = fion_dv(1:max_index) +&
                h_ion_dv(1:max_index)*(log(h_ion) + tc2*bi(1) - logqtl_const_hplus + 1._fp_kind)
        endif
     endif

     if(ifcsums) then
        ! Add in effect of diatomic and monatomic hydrogen ion species to Coulomb sums.
        sum0 = sum0 + sum0_scale*(h_ion + h2plus)
        sum2 = sum2 + sum2_scale*(h_ion + h2plus)
        if(ifnr03) then
           sum0f = sum0f + sum0_scale*(h_ionf + h2plusf)
           sum0t = sum0t + sum0_scale*(h_iont + h2plust)
           sum2f = sum2f + sum2_scale*(h_ionf + h2plusf)
           sum2t = sum2t + sum2_scale*(h_iont + h2plust)
        endif
        if(ifnr13) then
           sum0_dv(1:max_index) = sum0_dv(1:max_index) +&
                sum0_scale*(h_ion_dv(1:max_index) + h2plus_dv(1:max_index))
           sum2_dv(1:max_index) = sum2_dv(1:max_index) +&
                sum2_scale*(h_ion_dv(1:max_index) + h2plus_dv(1:max_index))
        endif
     endif !if(ifcsums)
     if(ifpl.eq.1) then
        ! correct sumpl0, sumpl1 and sumpl2 for hydrogen nu values.
        sumpl0 = sumpl0 +&
             h_neutral*plop(1) +&
             h2*plop(nions+1) +&
             h2plus*plop(nions+2)
        sumpl1 = sumpl1 + h_neutral*(plop(1) + plopt(1)) +&
             h2*(plop(nions+1) + plopt(nions+1)) +&
             h2plus*(plop(nions+2) + plopt(nions+2))
        sumpl2 = sumpl2 + h_neutral*plopt(1) +&
             h2*plopt(nions+1) + h2plus*plopt(nions+2)
        if(ifnr03) then
           sumpl0f = sumpl0f +&
                h_neutralf*plop(1) +&
                h2f*plop(nions+1) +&
                h2plusf*plop(nions+2)
           sumpl1f = sumpl1f + h_neutralf*(plop(1) + plopt(1)) +&
                h2f*(plop(nions+1) + plopt(nions+1)) +&
                h2plusf*(plop(nions+2) + plopt(nions+2))
           sumpl1t = sumpl1t + h_neutralt*(plop(1) + plopt(1)) +&
                h2t*(plop(nions+1) + plopt(nions+1)) +&
                h2plust*(plop(nions+2) + plopt(nions+2)) +&
                h_neutral*(plopt(1) + plopt2(1)) +&
                h2*(plopt(nions+1) + plopt2(nions+1)) +&
                h2plus*(plopt(nions+2) + plopt2(nions+2))
        endif
        if(ifnr13) then
           sumpl0_dv(1:max_index) = sumpl0_dv(1:max_index) +&
                h_neutral_dv(1:max_index)*plop(1) +&
                h2_dv(1:max_index)*plop(nions+1) +&
                h2plus_dv(1:max_index)*plop(nions+2)
        endif
     endif !if(ifpl.eq.1) then
     if(ifpi34) then
        ! correct extrasum for hydrogen nu values.
        rpower = 1._fp_kind
        rpowerh2 = 1._fp_kind
        rpowerh2plus = 1._fp_kind
        ! N.B. Note that extrasum and its derivatives are already scaled above.
        do iextrasum = 1,nextrasum-2
           extrasum(iextrasum) = extrasum(iextrasum) +&
                (h2*rpowerh2 + h2plus*rpowerh2plus + h_neutral*rpower)*extrasum_scale(iextrasum)
           if(ifnr03) then
              extrasumf(iextrasum) = extrasumf(iextrasum) +&
                   (h2f*rpowerh2 + h2plusf*rpowerh2plus + h_neutralf*rpower)*extrasum_scale(iextrasum)
              extrasumt(iextrasum) = extrasumt(iextrasum) +&
                   (h2t*rpowerh2 + h2plust*rpowerh2plus + h_neutralt*rpower)*extrasum_scale(iextrasum)
           endif
           if(ifnr13) then
              extrasum_dv(1:max_index,iextrasum) =&
                   extrasum_dv(1:max_index,iextrasum) +&
                   (h2_dv(1:max_index)*rpowerh2 +&
                   h2plus_dv(1:max_index)*rpowerh2plus +&
                   h_neutral_dv(1:max_index)*rpower)*extrasum_scale(iextrasum)
           endif
           rpower = rpower*r_neutral(1)
           rpowerh2 = rpowerh2*r_neutral(nelements+1)
           rpowerh2plus = rpowerh2plus*r_neutral(nelements+2)
        enddo !do iextrasum = 1,nextrasum-2
        ! sum over ions with weight of Z^1.5 which is unity for H2+ and H+.
        extrasum(nextrasum-1) = extrasum(nextrasum-1) + (h2plus + h_ion)*extrasum_scale(nextrasum-1)
        ! sum over all neutral and ionized states except bare nuclei
        extrasum(nextrasum) = extrasum(nextrasum) +&
             (h2plus*r_ion3(nions+2) + h2*r_ion3(nions+1) + h_neutral*r_ion3(1))*extrasum_scale(nextrasum)
        if(ifnr03) then
           extrasumf(nextrasum-1) = extrasumf(nextrasum-1) + (h2plusf + h_ionf)*extrasum_scale(nextrasum-1)
           extrasumt(nextrasum-1) = extrasumt(nextrasum-1) + (h2plust + h_iont)*extrasum_scale(nextrasum-1)
           extrasumf(nextrasum) = extrasumf(nextrasum) +&
                (h2plusf*r_ion3(nions+2) + h2f*r_ion3(nions+1) + h_neutralf*r_ion3(1))*extrasum_scale(nextrasum)
           extrasumt(nextrasum) = extrasumt(nextrasum) +&
                (h2plust*r_ion3(nions+2) + h2t*r_ion3(nions+1) + h_neutralt*r_ion3(1))*extrasum_scale(nextrasum)
        endif
        if(ifnr13) then
           extrasum_dv(1:max_index,nextrasum-1) = extrasum_dv(1:max_index,nextrasum-1) +&
                (h2plus_dv(1:max_index) + h_ion_dv(1:max_index))*&
                extrasum_scale(nextrasum-1)
           extrasum_dv(1:max_index,nextrasum) = extrasum_dv(1:max_index,nextrasum) +&
                (h2plus_dv(1:max_index)*r_ion3(nions+2) + h2_dv(1:max_index)*r_ion3(nions+1) +&
                h_neutral_dv(1:max_index)*r_ion3(1))*&
                extrasum_scale(nextrasum)
        endif

     endif
     ! End of logic block for (most common) case of partially ionized including hydrogen monatomics and molecules.
  endif

  ! N.B. xextrasum is calculated later so references in the next
  ! paragraph of commentary to auxiliary variables refers to all
  ! auxiliary variables *other than xextrasum*.

  ! The auxiliary variables and their derivatives that are calculated
  ! above are in scaled nu form so unscale to nu form and recalculate
  ! *_scale using
  !
  ! *_scale = aux_scale_limit_factor/abs(*aux)
  !
  ! so that the next time this routine is called
  ! the magnitude of the scaled nuform auxiliary variables will
  ! iteratively approximate aux_scale_limit_factor.

  if(ifcsums) then
     sum0 = sum0/sum0_scale
     if(sum0.lt.aux_underflow) sum0 = 0._fp_kind
     if(ifnr03) then
        if(sum0.gt.0._fp_kind) then
           sum0f = sum0f/sum0_scale
           sum0t = sum0t/sum0_scale
        else
           sum0f = 0._fp_kind
           sum0t = 0._fp_kind
        endif
     endif
     if(ifnr13) then
        if(sum0.gt.0._fp_kind) then
           sum0_dv(1:max_index) = sum0_dv(1:max_index)/sum0_scale
        else
           sum0_dv(1:max_index) = 0._fp_kind
        endif
     endif
     if(sum0.gt.0._fp_kind) sum0_scale = aux_scale_limit_factor/abs(sum0)

     sum2 = sum2/sum2_scale
     if(sum2.lt.aux_underflow) sum2 = 0._fp_kind
     if(ifnr03) then
        if(sum2.gt.0._fp_kind) then
           sum2f = sum2f/sum2_scale
           sum2t = sum2t/sum2_scale
        else
           sum2f = 0._fp_kind
           sum2t = 0._fp_kind
        endif
     endif
     if(ifnr13) then
        if(sum2.gt.0._fp_kind) then
           sum2_dv(1:max_index) = sum2_dv(1:max_index)/sum2_scale
        else
           sum2_dv(1:max_index) = 0._fp_kind
        endif
     endif
     if(sum2.gt.0._fp_kind) sum2_scale = aux_scale_limit_factor/abs(sum2)
  endif !if(ifcsums)
  if(ifpi34) then
     do iextrasum = 1, nextrasum
        extrasum(iextrasum) = extrasum(iextrasum)/extrasum_scale(iextrasum)
        if(extrasum(iextrasum).lt.aux_underflow) extrasum(iextrasum) = 0._fp_kind
        if(ifnr03) then
           if(extrasum(iextrasum).gt.0._fp_kind) then
              extrasumf(iextrasum) = extrasumf(iextrasum)/extrasum_scale(iextrasum)
              extrasumt(iextrasum) = extrasumt(iextrasum)/extrasum_scale(iextrasum)
           else
              extrasumf(iextrasum) = 0._fp_kind
              extrasumt(iextrasum) = 0._fp_kind
           endif
        endif
        if(ifnr13) then
           if(extrasum(iextrasum).gt.0._fp_kind) then
              extrasum_dv(1:max_index,iextrasum) = extrasum_dv(1:max_index,iextrasum)/extrasum_scale(iextrasum)
           else
              extrasum_dv(1:max_index,iextrasum) = 0._fp_kind
           endif
        endif
        if(extrasum(iextrasum).gt.0._fp_kind) extrasum_scale(iextrasum) = aux_scale_limit_factor/abs(extrasum(iextrasum))
     enddo ! do iextrasum = 1, nextrasum
  endif ! if(ifpi34) then

  if(ifexcited.gt.0) then
     ! Calculate xextrasum and derivatives when ifpi34 is .true.

     ! N.B. lextrasum and derivative values only used in excitation_sum call when
     ! ifpi34 is true.
     if(ifpi34) then
        ! Select the subset of the extrasum auxiliary variables and
        ! their derivatives that affect excitation.  Also convert
        ! from the currently calculated nu form of this subset to
        ! the n form which is always required by excitation_sum.
        do jextrasum = 1, nx
           if(jextrasum.lt.nx) then
              iextrasum = jextrasum
           else
              iextrasum = nextrasum-1
           endif
           lextrasum(iextrasum) = extrasum(iextrasum)*(rho*avogadro)
           if(ifnr03) then
              lextrasumf(iextrasum) = extrasumf(iextrasum)*(rho*avogadro) + lextrasum(iextrasum)*rf
              lextrasumt(iextrasum) = extrasumt(iextrasum)*(rho*avogadro) + lextrasum(iextrasum)*rt
           endif
           if(ifnr13) then
              lextrasum_dv(1:max_index,iextrasum) = extrasum_dv(1:max_index,iextrasum)*(rho*avogadro) +&
                   lextrasum(iextrasum)*r_dv(1:max_index)
           endif
        enddo !do jextrasum = 1, nx
     endif ! if(ifpi34) then

     call excitation_sum(verbosity, ifexcited, ifsame_zero_abundances,&
          ifpl, ifpi, ifmodified, ifnr, inv_ion, ifh2, ifh2plus,&
          partial_elements, ion_end,&
          tl, izlo, bmin, nmin, nmin_max, nmin_species, nmax,&
          bi, plop, plopt, plopt2,&
          r_ion3, nion, r_neutral,&
          lextrasum(:nextrasum), lextrasumf(:nextrasum), lextrasumt(:nextrasum), lextrasum_dv(:,:nextrasum),&
          max_index, h_neutral, h_neutralf, h_neutralt, h_neutral_dv,&
          h2, h2f, h2t, h2_dv,&
          h2plus, h2plusf, h2plust, h2plus_dv,&
          xextrasum, xextrasumf, xextrasumt, xextrasum_dv)

     ! xextrasum and its derivatives that are calculated in the
     ! excitation_sum routine are in scaled nu form so unscale to nu
     ! form and recalculate xextrasum_scale using
     !
     ! xextrasum_scale = aux_scale_limit_factor/abs(xextrasum)
     !
     ! so that the next time this routine is called
     ! the magnitude of the scaled nuform xextrasum will
     ! iteratively approximate aux_scale_limit_factor.

     if(ifpi34) then
        do ixextrasum = 1,nxextrasum
           xextrasum(ixextrasum) = xextrasum(ixextrasum)/xextrasum_scale(ixextrasum)
           if(abs(xextrasum(ixextrasum)).lt.aux_underflow) xextrasum(ixextrasum) = 0._fp_kind
           if(ifnr03) then
              if(abs(xextrasum(ixextrasum)).gt.0._fp_kind) then
                 xextrasumf(ixextrasum) = xextrasumf(ixextrasum)/xextrasum_scale(ixextrasum)
                 xextrasumt(ixextrasum) = xextrasumt(ixextrasum)/xextrasum_scale(ixextrasum)
              else
                 xextrasumf(ixextrasum) = 0._fp_kind
                 xextrasumt(ixextrasum) = 0._fp_kind
              endif
           endif
           if(ifnr13) then
              if(abs(xextrasum(ixextrasum)).gt.0._fp_kind) then
                 xextrasum_dv(1:max_index,ixextrasum) = xextrasum_dv(1:max_index,ixextrasum)/xextrasum_scale(ixextrasum)
              else
                 xextrasum_dv(1:max_index,ixextrasum) = 0._fp_kind
              endif
           endif
           if(abs(xextrasum(ixextrasum)).gt.0._fp_kind)&
                xextrasum_scale(ixextrasum) = aux_scale_limit_factor/abs(xextrasum(ixextrasum))
        enddo
     endif
  endif ! if(ifexcited.gt.0) then

  ! If appropriate flags demand it, convert nu form of sum0, sum2, h2,
  ! h2plus, h_ion, he_ion, he_ion2, extrasum, and xextrasum and
  ! derivatives to n form.
  if(ifionized.eq.2) then
     ! Start of logic block for case of everything fully ionized.
     if(.not.ifnuform) then
        ! should never reach this block of code because with
        ! full ionization in free_eos_detailed should have
        ! straight through option which should yield
        ! ifnuform = .true.
        error stop 'eos_calc: bad ifnuform setting for full ionization'

        ! calculate number densities for non-zero components.
        h_ion = eps(1)*rho*avogadro
        h_ionf = eps(1)*rf
        h_iont = eps(1)*rt
        he_ion2 = eps(2)*rho*avogadro
        he_ion2f = eps(2)*rf
        he_ion2t = eps(2)*rt
     endif !if(.not.ifnuform) then
     ! End of logic block for case of everything fully ionized.
  elseif(eps(1).le.0._fp_kind) then
     ! Start of logic block for case of partially ionized but no hydrogen species of any kind.
     ! convert he_ion, he_ion2, and derivatives to number densities.
     if(eps(2).gt.0._fp_kind) then
        ! if ifnuform is true then ifnr is zero and leave he_ion, he_ion2
        ! and their derivatives in the nu form returned by
        ! ionize.
        if(.not.ifnuform) then
           he_ion = he_ion*rho*avogadro
           he_ion2 = he_ion2*rho*avogadro
           if(ifpi.eq.2.or.if_mc.eq.1) then
              if(ifnr03) then
                 he_ionf = he_ionf*rho*avogadro + rf*he_ion
                 he_iont = he_iont*rho*avogadro + rt*he_ion
                 he_ion2f = he_ion2f*rho*avogadro + rf*he_ion2
                 he_ion2t = he_ion2t*rho*avogadro + rt*he_ion2
              endif
              if(ifnr13) then
                 index2 = inv_ion(2)
                 index3 = inv_ion(3)
                 do index = 1, max_index
                    if(index.eq.index2.or.index.eq.index3) then
                       he_ion_dv(index) = he_ion_dv(index)*rho*avogadro + he_ion*r_dv(index)
                       he_ion2_dv(index) = he_ion2_dv(index)*rho*avogadro + he_ion2*r_dv(index)
                    else
                       he_ion_dv(index) = he_ion*r_dv(index)
                       he_ion2_dv(index) = he_ion2*r_dv(index)
                    endif
                 enddo
              endif !if(ifnr13) then
           endif !if(ifpi.eq.2.or.if_mc.eq.1) then
        endif !if(.not.ifnuform) then
     endif !if(eps(2).gt.0._fp_kind) then
     if(ifcsums) then
        ! convert Coulomb sums to number densities.
        sum0 = sum0*(rho*avogadro)
        sum2 = sum2*(rho*avogadro)
        if(ifnr03) then
           sum0f = sum0f*(rho*avogadro) + sum0*rf
           sum0t = sum0t*(rho*avogadro) + sum0*rt
           sum2f = sum2f*(rho*avogadro) + sum2*rf
           sum2t = sum2t*(rho*avogadro) + sum2*rt
        endif
        if(ifnr13) then
           sum0_dv(1:max_index) = sum0_dv(1:max_index)*(rho*avogadro) + sum0*r_dv(1:max_index)
           sum2_dv(1:max_index) = sum2_dv(1:max_index)*(rho*avogadro) + sum2*r_dv(1:max_index)
        endif
     endif !if(ifcsums) then
     if(.not.ifnuform.and.ifpi34) then
        ! if ifnuform is true then leave
        ! extrasum and its derivative in the nu form returned by
        ! ionize.
        ! convert to n form from nu = n/(rho*avogadro) form
        extrasum(1:nextrasum) = extrasum(1:nextrasum)*(rho*avogadro)
        if(ifnr03) then
           extrasumf(1:nextrasum) = extrasumf(1:nextrasum)*(rho*avogadro) + extrasum(1:nextrasum)*rf
           extrasumt(1:nextrasum) = extrasumt(1:nextrasum)*(rho*avogadro) + extrasum(1:nextrasum)*rt
        endif
        if(ifnr13) then
           do iextrasum = 1, nextrasum
              extrasum_dv(1:max_index,iextrasum) = extrasum_dv(1:max_index,iextrasum)*(rho*avogadro) +&
                   extrasum(iextrasum)*r_dv(1:max_index)
           enddo
        endif
     endif !if(.not.ifnuform.and.ifpi34) then
     ! End of logic block for case of partially ionized but no hydrogen (or molecular hydrogen).
  elseif(ifh2.eq.0) then
     ! Start of logic block for case of partially ionized including monatomic forms of hydrogen,
     ! but excluding H2 and H2+.
     if(.not.ifnuform) then
        ! convert h_neutral, h_ion, and derivatives to number densities.
        h_neutral = h_neutral*rho*avogadro
        h_ion = h_ion*rho*avogadro
        if(ifnr03) then
           h_neutralf = h_neutralf*rho*avogadro + rf*h_neutral
           h_neutralt = h_neutralt*rho*avogadro + rt*h_neutral
           h_ionf = h_ionf*rho*avogadro + rf*h_ion
           h_iont = h_iont*rho*avogadro + rt*h_ion
        endif
        if(ifnr13) then
           do index = 1, max_index
              if(index.eq.indexh) then
                 ! before conversion h_neutral+h_ion = eps(1), hence h_neutral_dv = -h_ion_dv.
                 h_neutral_dv(index) = -h_ion_dv(index)*rho*avogadro + h_neutral*r_dv(index)
                 ! the following expression has significance loss so do
                 ! another way.
                 !h_dv(index) = h_ion_dv(index)*rho*avogadro + h_ion*r_dv(index)
                 h_ion_dv(index) = rho*avogadro*en*h_ion_dv(index)*hne - en*h_ion*hne_dv(index)
              else
                 h_neutral_dv(index) = h_neutral*r_dv(index)
                 h_ion_dv(index) = h_ion*r_dv(index)
              endif
           enddo !do index = 1, max_index
        endif !if(ifnr13) then
        if(eps(2).gt.0._fp_kind) then
           he_ion = he_ion*rho*avogadro
           he_ion2 = he_ion2*rho*avogadro
           if(ifpi.eq.2.or.if_mc.eq.1) then
              if(ifnr03) then
                 he_ionf = he_ionf*rho*avogadro + rf*he_ion
                 he_iont = he_iont*rho*avogadro + rt*he_ion
                 he_ion2f = he_ion2f*rho*avogadro + rf*he_ion2
                 he_ion2t = he_ion2t*rho*avogadro + rt*he_ion2
              endif
              if(ifnr13) then
                 index2 = inv_ion(2)
                 index3 = inv_ion(3)
                 do index = 1, max_index
                    if(index.eq.index2.or.index.eq.index3) then
                       he_ion_dv(index) = he_ion_dv(index)*rho*avogadro + he_ion*r_dv(index)
                       he_ion2_dv(index) = he_ion2_dv(index)*rho*avogadro + he_ion2*r_dv(index)
                    else
                       he_ion_dv(index) = he_ion*r_dv(index)
                       he_ion2_dv(index) = he_ion2*r_dv(index)
                    endif
                 enddo
              endif
           endif !if(ifpi.eq.2.or.if_mc.eq.1) then
        endif !if(eps(2).gt.0._fp_kind) then
        if(ifpi34) then
           ! convert to n form from nu = n/(rho*avogadro) form
           extrasum(1:nextrasum) = extrasum(1:nextrasum)*(rho*avogadro)
           if(ifnr03) then
              extrasumf(1:nextrasum) = extrasumf(1:nextrasum)*(rho*avogadro) + extrasum(1:nextrasum)*rf
              extrasumt(1:nextrasum) = extrasumt(1:nextrasum)*(rho*avogadro) + extrasum(1:nextrasum)*rt
           endif
           if(ifnr13) then
              ! Could implement this with array ops and appropriate rank
              ! promotion to form a matrix from product of two
              ! vectors, but do loop probably faster.
              do iextrasum = 1, nextrasum
                 extrasum_dv(1:max_index,iextrasum) =&
                      extrasum_dv(1:max_index,iextrasum)*(rho*avogadro) +&
                      extrasum(iextrasum)*r_dv(1:max_index)
              enddo
           endif
        endif
     endif !if(.not.ifnuform) then
     if(ifcsums) then
        ! convert Coulomb sums to number densities
        sum0 = sum0*(rho*avogadro)
        sum2 = sum2*(rho*avogadro)
        if(ifnr03) then
           sum0f = sum0f*(rho*avogadro) + sum0*rf
           sum0t = sum0t*(rho*avogadro) + sum0*rt
           sum2f = sum2f*(rho*avogadro) + sum2*rf
           sum2t = sum2t*(rho*avogadro) + sum2*rt
        endif
        if(ifnr13) then
           sum0_dv(1:max_index) = sum0_dv(1:max_index)*(rho*avogadro) + sum0*r_dv(1:max_index)
           sum2_dv(1:max_index) = sum2_dv(1:max_index)*(rho*avogadro) + sum2*r_dv(1:max_index)
        endif
     endif !if(ifcsums) then
     ! End of logic block for case of partially ionized including monatomic forms of hydrogen,
     ! but excluding H2 and H2+.
  else !if(ifionized.eq.2) then
     ! Start of logic block for (most common) case of partially ionized including hydrogen monatomics and molecules.
     if(.not.ifnuform) then
        ! replace nu form of h2, h2plus, h_neutral, h_ion, and derivatives with stored n form.
        h2 = nh2
        h2plus = nh2plus
        h_neutral = nh_neutral
        h_ion = nh_ion
        if(ifnr03) then
           h2f = nh2f
           h2t = nh2t
           h2plusf = nh2plusf
           h2plust = nh2plust
           h_neutralf = nh_neutralf
           h_neutralt = nh_neutralt
           h_ionf = nh_ionf
           h_iont = nh_iont
        endif
        if(ifnr13) then
           h2_dv(1:max_index) = nh2_dv(1:max_index)
           h2plus_dv(1:max_index) = nh2plus_dv(1:max_index)
           h_neutral_dv(1:max_index) = nh_neutral_dv(1:max_index)
           h_ion_dv(1:max_index) = nh_ion_dv(1:max_index)
        endif

        ! convert he_ion, he_ion2, and derivatives to number densities.
        if(eps(2).gt.0._fp_kind) then
           he_ion = he_ion*rho*avogadro
           he_ion2 = he_ion2*rho*avogadro
           if(ifpi.eq.2.or.if_mc.eq.1) then
              if(ifnr03) then
                 he_ionf = he_ionf*rho*avogadro + rf*he_ion
                 he_iont = he_iont*rho*avogadro + rt*he_ion
                 he_ion2f = he_ion2f*rho*avogadro + rf*he_ion2
                 he_ion2t = he_ion2t*rho*avogadro + rt*he_ion2
              endif
              if(ifnr13) then
                 index2 = inv_ion(2)
                 index3 = inv_ion(3)
                 do index = 1, max_index
                    if(index.eq.index2.or.index.eq.index3) then
                       he_ion_dv(index) = he_ion_dv(index)*rho*avogadro + he_ion*r_dv(index)
                       he_ion2_dv(index) = he_ion2_dv(index)*rho*avogadro + he_ion2*r_dv(index)
                    else
                       he_ion_dv(index) = he_ion*r_dv(index)
                       he_ion2_dv(index) = he_ion2*r_dv(index)
                    endif
                 enddo !do index = 1,max_index
              endif
           endif
        endif
        if(ifpi34) then
           ! convert to n form from nu = n/(rho*avogadro) form
           extrasum(1:nextrasum) = extrasum(1:nextrasum)*(rho*avogadro)
           if(ifnr03) then
              extrasumf(1:nextrasum) = extrasumf(1:nextrasum)*(rho*avogadro) + extrasum(1:nextrasum)*rf
              extrasumt(1:nextrasum) = extrasumt(1:nextrasum)*(rho*avogadro) + extrasum(1:nextrasum)*rt
           endif
           if(ifnr13) then
              ! Could implement this with array ops and appropriate rank
              ! promotion to form a matrix from product of two
              ! vectors, but do loop probably faster.
              do iextrasum = 1, nextrasum
                 extrasum_dv(1:max_index,iextrasum) =&
                      extrasum_dv(1:max_index,iextrasum)*(rho*avogadro) +&
                      extrasum(iextrasum)*r_dv(1:max_index)
              enddo
           endif
        endif
     endif !if(.not.ifnuform) then
     if(ifcsums) then
        ! convert Coulomb sums to number densities
        sum0 = sum0*(rho*avogadro)
        sum2 = sum2*(rho*avogadro)
        if(ifnr03) then
           sum0f = sum0f*(rho*avogadro) + sum0*rf
           sum0t = sum0t*(rho*avogadro) + sum0*rt
           sum2f = sum2f*(rho*avogadro) + sum2*rf
           sum2t = sum2t*(rho*avogadro) + sum2*rt
        endif
        if(ifnr13) then
           sum0_dv(1:max_index) = sum0_dv(1:max_index)*(rho*avogadro) + sum0*r_dv(1:max_index)
           sum2_dv(1:max_index) = sum2_dv(1:max_index)*(rho*avogadro) + sum2*r_dv(1:max_index)
        endif
     endif !if(ifcsums) then
     ! End of logic block for (most common) case of partially ionized including hydrogen monatomics and molecules.
  endif

  if(ifexcited.gt.0) then
     if(.not.ifnuform.and.ifpi34) then
        xextrasum(1:nxextrasum) = xextrasum(1:nxextrasum)*(rho*avogadro)
        if(ifnr03) then
           xextrasumf(1:nxextrasum) = xextrasumf(1:nxextrasum)*(rho*avogadro) + xextrasum(1:nxextrasum)*rf
           xextrasumt(1:nxextrasum) = xextrasumt(1:nxextrasum)*(rho*avogadro) + xextrasum(1:nxextrasum)*rt
        endif
        if(ifnr13) then
           do ixextrasum = 1,nxextrasum
              xextrasum_dv(1:max_index, ixextrasum) =&
                   xextrasum_dv(1:max_index, ixextrasum)*(rho*avogadro) + xextrasum(ixextrasum)*r_dv(1:max_index)
           enddo
        endif
     endif ! if(.not.ifnuform.and.ifpi34) then
  endif ! if(ifexcited.gt.0) then

  if(.false..and.verbosity.ge.4) then
     write(stderr,*) "h2, h2f, h2t = ", achar(10), h2, h2f, h2t
     write(stderr,*) "h2plus, h2plusf, h2plust = ", achar(10), h2plus, h2plusf, h2plust
     write(stderr,*) "h_neutral, h_neutralf, h_neutralt = ", achar(10), h_neutral, h_neutralf, h_neutralt
     write(stderr,*) "h_ion, h_ionf, h_iont = ", achar(10), h_ion, h_ionf, h_iont
     write(stderr,*) "he_ion, he_ionf, he_iont = ", achar(10), he_ion, he_ionf, he_iont
     write(stderr,*) "he_ion2, he_ion2f, he_ion2t = ", achar(10), he_ion2, he_ion2f, he_ion2t
     write(stderr,*) "h2_dv(1:max_index) = ", achar(10), h2_dv(1:max_index)
     write(stderr,*) "h2plus_dv = ", achar(10), h2plus_dv
     write(stderr,*) "h_neutral_dv = ", achar(10), h_neutral_dv
     write(stderr,*) "h_ion_dv(1:max_index) = ", achar(10), h_ion_dv(1:max_index)
     write(stderr,*) "he_ion_dv(1:max_index) = ", achar(10), he_ion_dv(1:max_index)
     write(stderr,*) "he_ion2_dv(1:max_index) = ", achar(10), he_ion2_dv(1:max_index)
     write(stderr,*) "ne, nef, net = ", achar(10), ne, nef, net
     write(stderr,*) "hne, hnef, hnet = ", achar(10), hne, hnef, hnet
     write(stderr,*) "hne_dv(1:max_index) = ", achar(10), hne_dv(1:max_index)
     write(stderr,*) "sion, sionf, siont, uion = ", achar(10), sion, sionf, siont, uion
     write(stderr,*) "fion, fionf = ", achar(10), fion, fionf
     write(stderr,*) "fion_dv(1:max_index) = ", achar(10), fion_dv(1:max_index)
     write(stderr,*) "sum0, sum0f, sum0t = ", achar(10), sum0, sum0f, sum0t
     write(stderr,*) "sum0_dv(1:max_index) = ", achar(10), sum0_dv(1:max_index)
     write(stderr,*) "sum2, sum2f, sum2t = ", achar(10), sum2, sum2f, sum2t
     write(stderr,*) "sum2_dv(1:max_index) = ", achar(10), sum2_dv(1:max_index)
     write(stderr,*) "extrasum(:) = ", achar(10), extrasum(:)
     write(stderr,*) "extrasumf(:) = ", achar(10), extrasumf(:)
     write(stderr,*) "extrasumt(:) = ", achar(10), extrasumt(:)
     write(stderr,*) "extrasum_dv(1:max_index,:) = ", achar(10), extrasum_dv(1:max_index,:)
     write(stderr,*) "xextrasum(:) = ", achar(10), xextrasum(:)
     write(stderr,*) "xextrasumf(:) = ", achar(10), xextrasumf(:)
     write(stderr,*) "xextrasumt(:) = ", xextrasumt(:)
     write(stderr,*) "xextrasum_dv(1:max_index,:) = ", achar(10), xextrasum_dv(1:max_index,:)
     write(stderr,*) "sumpl0, sumpl0f = ", achar(10), sumpl0, sumpl0f
     write(stderr,*) "sumpl0_dv(1:max_index) = ", achar(10), sumpl0_dv(1:max_index)
     write(stderr,*) "sumpl1, sumpl1f, sumpl1t, sumpl2 = ", achar(10), sumpl1, sumpl1f, sumpl1t, sumpl2
     write(stderr,*) "rl, rf, rt = ", achar(10), rl, rf, rt
     write(stderr,*) "r_dv(1:max_index) = ", achar(10), r_dv(1:max_index)
  endif
end subroutine eos_calc
