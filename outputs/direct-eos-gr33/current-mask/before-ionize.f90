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

!> The purpose of this ionize subroutine is to calculate number
!> densities (in nu form where nu = n/(rho*avogadro)) and
!> non-excitation auxiliary variables (appropriate weighted sums over those
!> quantities) for all monatomic species.  For hydrogen, these results
!> are corrected in the calling (eos_calc) subroutine for the presence
!> of the molecules H2 and H2+ and optionally converted to n form.  These corrected auxiliary variables
!> (and excitation auxiliary variables once those are calculated) are
!> ultimately used in an iterative way to calculate chemical potentials
!> according to some free energy model.  The appropriate difference in
!> chemical potentials are then combined to form the dv quantities, the
!> change in the equilibrium constant of an ion relative to the
!> un-ionized reference state and the free electron
!>
!> dv = (partial F/partial nref - partial F/partial nion - ion*partial F/partial ne)/kT.
!>
!> These dv values (held in an array) are required input to
!> ionize.
!>
!> The free energy model:
!> it is the responsibility of the calling programme to calculate the
!> dv quantities for the various ions in a consistent manner.  It
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
!> \param[in] ifsame_under
!>   ifsame_under = .true. means use same underflow zeroing for each
!>     ion as in previous call.  This option is useful for removing
!>     small discontinuities caused by variations in the underflow
!>     zeroing which sometimes foil the last stages of
!>     convergence.<br>
!>   ifsame_under = .false. means calculate underflow zeroing.<br>
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
!>     element.<br>
!> \param[in] ifionized
!>   ifionized = 0 means all elements treated as partially ionized.<br>
!>   ifionized = 1 means trace metals fully ionized.<br>
!>   ifionized = 2 means all elements treated as fully ionized.<br>
!> \param[in] ifcsums
!>   if ifcsums is .true. then sum0 and sum2 (the Coulomb sum auxiliary variables) should be calculated.<br>
!>   if ifcsums is .false. then sum0 and sum2 (the Coulomb sum auxiliary variables) should *not* be calculated.<br>
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
!> \param[in] tc2
!>   tc2 = c2/t.<br>
!> \param[in] bi
!>   bi(nions_inp2) ionization potentials (cm**-1) in ion order from
!>     neutral to last ion before the bare nucleus ion, e.g., H
!>     (excluding the H+ bare nucleus); He, He+ (excluding the He++
!>     bare nucleus); C through C+++++ (excluding the C++++++ bare
!>     nucleus); etc.  The H2 and H2+ ionization potentials are the
!>     last two in this list.<br>
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
!> \param[out] hne hne is the calculated nu(e) = n_e/(Navogadro*rho) from
!>     the metals and helium (and hydrogen if treated as fully
!>     ionized).<br>
!> \param[out] hnef hnef is the calculated partial derivative of hne wrt lnf.<br>
!> \param[out] hnet hnet is the calculated partial derivative of hne wrt ln t.<br>
!> \param[out] hne_dv hne_dv(nions_inp2) is an array of the calculated partial derivatives of hne wrt the elements of dv.<br>
!> \param[out] s
!>   R*s is the calculated entropy per unit mass corresponding to
!>     the ideal component of the free energy.  n.b. for the
!>     completely ionized case s contains both the hydrogen and non-hydrogen components of the entropy.  For the partially
!>     ionized case without molecular formation, this routine
!>     separates the hydrogen component of s from the rest of the
!>     elements (see the sh argument of this routine).  For the partially
!>     ionized case with H2 and possibly H2+ formation, sh is ignored
!>     by the eos_calc calling routine and instead that hydrogen component is
!>     calculated properly in that routine for the case given and
!>     added to s calculated by ionize for the remaining
!>     elements.<br>
!>   Units check: From the formula below, s has the same units as
!>     nu = n/(rho*NA), i.e., the inverse of the atomic mass.
!>     Therefore, R s = k*NA s has units of energy per degree K
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
!>   From the ideal component of the free energy s is the sum over
!>     all non-electron species of<br>
!>       -nu_i*[-5/2 - 3/2 ln T + ln(alpha) - ln Na +<br>
!>       ln (n_i/[A_i^{3/2} Q_i]) - dln Q_i/d ln T],<br>
!>       where nu_i = n_i/(Na rho),<br>
!>       Na is the Avogadro number,<br>
!>       A_i is the atomic weight,<br>
!>       alpha =  (2 pi k/[Na h^2])^(-3/2) Na
!>       (see documentation of alpha^2 in the mod_free_eos_constants
!>       module), and<br>
!>       Q_i is the internal ideal partition function of
!>       the non-Rydberg states.  (Currently, we calculate this
!>       partition function by the statistical weights of the ground
!>       states of monatomic hydrogen and helium and the combined
!>       lower states (roughly approximated) of each of the metal.
!>       This is a superb approximation for hydrogen, helium is
!>       subsequently corrected for detailed excitation of the
!>       non-Rydberg states, and this crude approximation for the
!>       metals is not currently corrected.  Thus, in all *current*
!>       monatomic cases ln Q_i is a constant, and dln Q_i/d ln T is
!>       zero, but this will be updated for the metals eventually.)<br>
!>       From the equilibrium constant approach and the monatomic species
!>       treated in this subroutine (molecules treated outside) we have<br>
!>       -nu_i ln (n_i/[A_i^{3/2} Q_i]) =<br>
!>       -nu_i * [ln (n_neutral/[A_neutral^{3/2} Q_neutral) +<br>
!>       (- chi_i/kT + dv_i)].<br>
!>       If we sum this term over all species
!>       of an element without molecules (true for all species treated here
!>       since molecular hydrogen treated specially outside when important)
!>       we obtain<br>
!>       s_element = - eps * ln (eps*n_neutral/sum(n))<br>
!>       - eps ln (alpha/(A_neutral^{3/2} Q_neutral)) -<br>
!>       - sum over all species of the element of nu_i*(-chi/kT + dv(i))<br>
!>       where we have ignored the term<br>
!>       eps * (-5/2 - 3/2 ln T + ln rho)<br>
!>       (taking into account the first 4 terms above)
!>       until a final correction of the s zero point
!>       in free_eos_detailed.<br>
!> \param[out] sf sf is the calculated partial derivative of s wrt fl.<br>
!> \param[out] st st is the calculated partial derivative of s wrt ln t.<br>
!> \param[out] u
!>   c2*cr*u is the calculated internal energy per unit mass (excluding the
!>     hydrogen components) corresponding to the the ideal component
!>     of the free energy.<br>
!>   Units check: u has the same units as energy (in cm^{-1}) times nu
!>     = n/(rho*NA).  The factor h*c converts energy in cm^{-1} to
!>     energy in ergs.  Therefore, c2*cr*u = h*c*NA*u has the expected
!>     cgs units of energy per unit mass.<br>
!> \param[out] h_neutral
!>   h_neutral is the calculated hydrogen neutral fraction (n(H)/(n(H)+n(H+)))
!>   times eps(1).  If there is no molecular formation of H2 or H2+,
!>   then<br>
!>   n(H)+n(H+) = eps(1)*rho*NA, and<br>
!>   h_neutral = n(H)/(rho*NA).<br>
!> \param[out] h_neutralf h_neutralf is the calculated partial derivative of h_neutral wrt fl.<br>
!> \param[out] h_neutralt h_neutralt is the calculated partial derivative of h_neutral wrt ln t.<br>
!> \param[out] sh
!>   sh is the calculated hydrogen contribution to the entropy s.  For the full
!>     ionization case, these terms are not calculated but instead are
!>     included in s.  For the case where hydrogen
!>     molecular formation is important, these terms are not valid and
!>     should be ignored.<br>
!> \param[out] shf shf is the calculated partial derivative of sh wrt fl.<br>
!> \param[out] sht sht is the calculated partial derivative of sh wrt ln t.<br>
!> \param[out] h_ion_equil
!>   h_ion_equil is the calculated h_ion equilibrium constant ln [n(H+)/n(H)] (only
!>     used for the case of h2 and h2+ formation).<br>
!> \param[out] h_ion_equilf h_ion_equilf is the calculated partial derivative of h_ion_equil wrt fl.<br>
!> \param[out] h_ion_equilt h_ion_equilt is the calculated partial derivative of h_ion_equil wrt ln t.<br>
!> \param[out] h_ion
!>   h_ion is the calculated hydrogen ionization fraction (n(H+)/(n(H)+n(H+)))
!>     times eps(1).  If no molecular formation of H2 and H2+, then
!>     n(H)+n(H+) = eps(1)*rho*NA, and h_ion = n(H+)/(rho*NA).<br>
!> \param[out] h_ionf h_ionf is the calculated partial derivative of h_ion wrt fl.<br>
!> \param[out] h_iont h_iont is the calculated partial derivative of h_ion wrt ln t.<br>
!> \param[out] h_ion_dv h_ion_dv(nions_inp2) is an array of the calculated partial derivatives of h_ion wrt the elements of dv.<br>
!> \param[out] he_ion
!>   he_ion is the calculated helium first ionization fraction n(He+)/(rho*NA).<br>
!> \param[out] he_ionf he_ionf is the partial derivative of he_ion wrt fl.<br>
!> \param[out] he_iont he_iont is the partial derivative of he_ion wrt ln t.<br>
!> \param[out] he_ion_dv he_ion_dv(nions_inp2) is an array of the partial derivatives of he_ion wrt the elements of dv.<br>
!> \param[out] he_ion2
!>   he_ion2 is the calculated helium second ionization fraction n(He++)/(rho*NA).<br>
!> \param[out] he_ion2f he_ion2f is the partial derivative of he_ion2 wrt fl.<br>
!> \param[out] he_ion2t he_ion2t is the partial derivative of he_ion2 wrt ln t.<br>
!> \param[out] he_ion2_dv he_ion2_dv(nions_inp2) is an array of the partial derivatives of he_ion2 wrt the elements of dv.<br>
!> \param[out] fionh
!>   fionh (returned only if ifpi = 3 or 4) is the calculated hydrogen ion
!>     free_energy/(rho*kT*NA) (= tc2*uh - sh at equilibrium).<br>
!> \param[out] fionhf fionhf is the calculated partial derivative of fionh wrt fl.<br>
!> \param[out] fionh_dv fionh_dv(nions_inp2) is an array of the calculated partial derivatives of fionh wrt the elements of dv.<br>
!> \param[out] fion
!>   fion (returned only if ifpi = 3 or 4) is the calculated ion
!>     free_energy/(rho*kT*NA) (= tc2*u - s at equilibrium).<br>
!> \param[out] fionf fionf is the partial derivative of fion wrt fl.<br>
!> \param[out] fion_dv fion_dv(nions_inp2) is an array of the partial derivatives of fion wrt the elements of dv.<br>
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
!> \param[out] info
!>   info is the calculated return code for the procedure which is
!>     zero for success and non-zero (with offset to identify which
!>     procedure the error occurred in) if an error occurred within
!>     the procedure or within some procedure it called.<br>

subroutine ionize(verbosity, ifexcited, ifsame_under, ifnr, inv_ion,&
     max_index, partial_elements, ion_end,&
     ifionized, ifcsums, ifreducedmass,&
     ifsame_abundances, ifmtrace, iatomic_number, ifpi, ifpl,&
     eps, tc2, bi, plop, plopt, plopt2,&
     r_ion3, r_neutral,&
     ifelement, dvzero,&
     dv, dvf, dvt,&
     nion,&
     hne, hnef, hnet,hne_dv,&
     s, sf, st,&
     u,&
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

  use mod_free_eos_constants, only: logalpha, electron_mass
  use mod_nuvar, only: nuvar, nuvarf, nuvart, nuvar_dv, nuvar_index_element, nuvar_atomic_number, nuvar_nelements
  use mod_aux_scale, only: sum0_scale, sum2_scale, extrasum_scale
  use mod_statistical_weight_data, only: nions_stat, nelements_stat, iqion, iqneutral
  use mod_isotopic_mass_data, only: nelements_iso, isotopic_mass
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real
  use mod_info_data, only: info_offset_ionize
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  logical, intent(in) :: ifsame_under, ifcsums
  integer, intent(in) :: verbosity, ifexcited, ifnr,&
       partial_elements(:),&
       ifionized, ifreducedmass, ifpi, ifpl,&
       ifelement(:),&
       max_index,&
       inv_ion(:),&
       ion_end(:)
  integer, intent(in) :: nion(:)  !ionization state in ion order.
  integer, intent(in) :: ifsame_abundances, ifmtrace, iatomic_number(:)

  integer, intent(out) :: info

  real(fp_kind), intent(in) :: eps(:), tc2, dvzero(:), dv(:), dvf(:), dvt(:), r_ion3(:), r_neutral(:)
  real(fp_kind), intent(out) ::&
       hne, hnef, hnet, hne_dv(:),&
       sum0, sum0f, sum0t, sum0_dv(:),&
       sum2, sum2f, sum2t, sum2_dv(:),&
       s, sf, st, u, h_neutral, h_neutralf, h_neutralt, sh, shf, sht,&
       h_ion_equil, h_ion_equilf, h_ion_equilt,&
       fionh, fionhf, fionh_dv(:),&
       fion, fionf, fion_dv(:),&
       extrasum(:), extrasumf(:), extrasumt(:), extrasum_dv(:,:),&
       sumpl0, sumpl0f, sumpl0_dv(:),&
       sumpl1, sumpl1f, sumpl1t, sumpl2,&
       h_ion, h_ionf, h_iont, he_ion, he_ionf, he_iont, he_ion2, he_ion2f, he_ion2t,&
       h_ion_dv(:), he_ion_dv(:), he_ion2_dv(:)
  ! ionization potentials in cm-1 in ion order for the ions,
  ! H+,
  ! He+, He++,
  ! Li, Be, B missing from list
  ! C+ through C6+,
  ! N+ through N7+,
  ! O+ through O8+,
  ! F missing from list
  ! Ne+ through Ne10+,
  ! Na+ through Na11+,
  ! Mg+ through Mg12+,
  ! Al+ through Al13+,
  ! Si+ through Si14+,
  ! P+ through P15+,
  ! S+ through S16+,
  ! Cl+ through Cl17+,
  ! A+ through A18+,
  ! K missing from list
  ! Ca+ through Ca20+,
  ! Sc missing from list
  ! Ti+ through Ti22+,
  ! V missing from list
  ! Cr+ through Cr24+,
  ! Mn+ through Mn25+,
  ! Fe+ through Fe26+,
  ! Co missing from list
  ! Ni+ through Ni28+,
  real(fp_kind), intent(in) :: bi(:), plop(:), plopt(:), plopt2(:)

  ! Parameters

  ! I believe this parameter chooses which of two different sets of
  ! ionization fraction logic is to be used below.  I assume iffion =
  ! .false. corresponds to a historical version of that logic that has
  ! now been superseded by the iffion = .true.  version of this logic,
  ! and, in fact, iffion = .false. may not even work any more.
  logical, parameter :: iffion = .true.

  ! Number of elements considered.
  integer, parameter :: nelements = 24
  ! Number of ionized monatomic species of those elements to be considered.
  integer, parameter :: nions = 316
  ! Treat up to 28 ions of one element, i.e., Nickel has the maximum
  ! atomic number of all elements treated.
  integer, parameter :: maxionstage = 28

  ! exp(2*arglim) is the maximum ratio of largest to least ionization
  ! component.
  ! n.b. typical overflow limit ~ 1.d305, but we back off
  ! by more than a factor of two because subsequent calculations use
  ! the square of these factors, and we want to leave some room for
  ! bad scaling.
  ! exp(2*2*arglim) = exp(2*2*144) ~ 1.d250 which gives ~50 orders
  ! of magnitude of room for bad scaling.  Thus, 144 is the maximum
  ! arglim.

  ! N.B. 144.d0 is not gross overkill.  The whole EOS is expressed in
  ! terms of ne so *some* ionized species must be present.  Some of
  ! the auxiliary variables depend only on neutrals, only on ions, or
  ! only on non-bare ions.  Thus, near full neutralization or near
  ! full ionization can cause discontinuous behaviour whenever some
  ! species is zeroed because of underflow concerns.  These
  ! discontinuities shouldn't matter for converged EOS results, but
  ! the discontinuities will mess up convergence criteria if the
  ! discontinuity occurs near the converged solution.  Even with a
  ! criterion of 144.d0, we still found this behaviour so ifsame_under
  ! logic was introduced to avoid these discontinuities (see
  ! ifsame_under logic comments above).
  real(fp_kind), parameter :: arglim = 144._fp_kind

  ! Internal variables

  logical ifneutral_zero(nelements), ifion_zero(nions)

  logical ifpi34, ifnr03, ifnr13, ifcsums_old

  real(fp_kind), allocatable ::&
       fract(:),&
       fractf(:),&
       fractt(:),&
       vz_dv(:),&
       rpower_dv(:)

  integer nelements_in
  integer nions_in, nions_inp2
  integer nextrasum
  integer ion, ielement, index_f, jndex_f, ion0, ion_index,&
       inv_ion_index, index_max, index_maxnb, iextrasum,&
       index_element
  integer n_partial_elements

  integer iffirst
  data iffirst/1/

  integer ifreducedmassold
  ! need invalid value to force initialization
  data ifreducedmassold/-1/

  integer ifionized_old, ifpi_old, ifmtrace_old
  ! Need invalid values to force initialization
  data ifionized_old, ifpi_old, ifmtrace_old/3*-10000/
  ! Pick an arbitrary value since initialization is assured in
  ! any case by the above data statements.
  data ifcsums_old /.true./

  real(fp_kind) ionized_hne, ionized_sum0, ionized_sum2, ionized_u, ionized_s
  ! These cannot be allocatable without a memory leak because they are all saved.
  real(fp_kind) bi_ref(nions), ce0(nions), ce(nions),&
       logqtl_const0(nelements), logqtl_const(nelements+nions+2),&
       ionized_extrasum(2)

  real(fp_kind) constant_sh, constant_s
  real(fp_kind) fractneutral, ln_fractmax, fractmax, sum_full, logsum0, logsum0f,&
       arg,&
       sum, sumf, sumt, sumnb,&
       vz, vzf, vzt, sarg, sarg0, epsilon, rpower, rpowerf,&
       rpowert

  ! iffirst initialized data
  save iffirst, bi_ref, ce0

  ! ifreducedmassold initialized data
  save ifreducedmassold, ce, logqtl_const0, logqtl_const

  ! ifsame_abundances data calculated on first but not necessarily all subsequent calls.
  save constant_sh, constant_s

  ! ifionized_old, etc., initialized data
  save ifionized_old, ifpi_old, ifcsums_old, ifmtrace_old,&
       ionized_hne, ionized_u, ionized_s, ionized_sum0, ionized_sum2, ionized_extrasum

  ! ifsame_under saved data
  save ifneutral_zero, ifion_zero

  ! Assume success by default unless an error condition is discovered below.
  info = 0

  n_partial_elements = size(partial_elements) - 2
  nelements_in = size(ifelement)
  nextrasum = size(extrasum)
  nions_inp2 = size(nion)
  nions_in = nions_inp2 - 2

  ! sanity checks:
  if(ifnr.lt.0.or.ifnr.eq.2.or.ifnr.gt.3) error stop 'ionize: ifnr must be 0, 1, or 3'
  if(nions_in.ne.nions.or.nions_in.ne.nions_stat) error stop 'ionize: invalid nions_in'
  if(nelements_in.ne.nelements.or.nelements_in.ne.nelements_stat.or.nelements_in.ne.nelements_iso)&
       error stop 'ionize: inconsistent nelements values'
  if(ifpi.lt.0.or.ifpi.gt.4) error stop 'ionize: invalid ifpi'
  if(.not.((ifpi.ne.3.and.ifpi.ne.4.and.nextrasum.eq.0).or.&
       ((ifpi.eq.3.or.ifpi.eq.4).and.(nextrasum.eq.6.or.nextrasum.eq.9))))&
       error stop 'ionize: invalid nextrasum'
  if(iffirst.eq.1.and.ifsame_abundances.eq.1) error stop 'ionize: invalid ifsame_abundances on first call'
  if(&
       nelements_in+2.ne.size(ion_end).or.&
       nelements_in.ne.size(eps).or.&
       nelements_in.ne.size(dvzero).or.&
       nelements_in+2.ne.size(r_neutral).or.&
       nelements_in.ne.size(iatomic_number))&
       error stop 'ionize: inconsistent nelements sizes'
  if(&
       nextrasum.ne.size(extrasumf).or.&
       nextrasum.ne.size(extrasumt).or.&
       nextrasum.ne.size(extrasum_dv,2))&
       error stop 'ionize: inconsistent nextrasum sizes'
  if(&
       nions_inp2.ne.size(inv_ion).or.&
       nions_inp2.ne.size(r_ion3).or.&
       nions_inp2.ne.size(extrasum_dv,1).or.&
       nions_inp2.ne.size(dv).or.&
       nions_inp2.ne.size(dvf).or.&
       nions_inp2.ne.size(dvt).or.&
       nions_inp2.ne.size(hne_dv).or.&
       nions_inp2.ne.size(sum0_dv).or.&
       nions_inp2.ne.size(sum2_dv).or.&
       nions_inp2.ne.size(fionh_dv).or.&
       nions_inp2.ne.size(fion_dv).or.&
       nions_inp2.ne.size(sumpl0_dv).or.&
       nions_inp2.ne.size(h_ion_dv).or.&
       nions_inp2.ne.size(he_ion_dv).or.&
       nions_inp2.ne.size(he_ion2_dv).or.&
       nions_inp2.ne.size(bi).or.&
       nions_inp2.ne.size(plop).or.&
       nions_inp2.ne.size(plopt).or.&
       nions_inp2.ne.size(plopt2))&
       error stop 'ionize: inconsistent nions_inp2 sizes'

  allocate(&
       fract(maxionstage),&
       fractf(maxionstage),&
       fractt(maxionstage),&
       vz_dv(maxionstage),&
       rpower_dv(maxionstage))

  if(if_taint_allocated_real) then
     call taint_allocated_real(fract)
     call taint_allocated_real(fractf)
     call taint_allocated_real(fractt)
     call taint_allocated_real(vz_dv)
     call taint_allocated_real(rpower_dv)
  endif

  if(iffirst.eq.1) then
     iffirst = 0
     ! refer energies to neutral species.
     do ion = 1,nions
        bi_ref(ion) = bi(ion)
        if(nion(ion).gt.1) bi_ref(ion) = bi_ref(ion) + bi_ref(ion-1)
     enddo
     ielement = 0
     do ion = 1,nions
        if(nion(ion).eq.1) ielement = ielement + 1
        ce0(ion) = log(real((iqion(ion)),fp_kind)/real((iqneutral(ielement)),fp_kind))
     enddo
     if(ifsame_under) error stop 'ionize: ifsame_under must be .false. on first call'
  endif
  ! useful combinations and flags in following tests:
  ifpi34 = ifpi.eq.3.or.ifpi.eq.4
  ifnr03 = ifnr.eq.0.or.ifnr.eq.3
  ifnr13 = ifnr.eq.1.or.ifnr.eq.3
  if(ifreducedmass.ne.ifreducedmassold) then
     ifreducedmassold = ifreducedmass
     ielement = 0
     do ion = 1,nions
        ce(ion) = ce0(ion)
        if(nion(ion).eq.1) then
           ielement = ielement + 1
           ! temperature-independent part of log product of translational and
           ! internal partition functions/avogadro
           ! logqtl_const0 is the part that is constant for each element
           logqtl_const0(ielement) = 1.5_fp_kind*log(isotopic_mass(ielement)) - logalpha
           logqtl_const(ielement) = log(real((iqneutral(ielement)),fp_kind)) + logqtl_const0(ielement)
        endif
        logqtl_const(nelements+ion) = log(real((iqion(ion)),fp_kind)) + logqtl_const0(ielement)
        if(ifreducedmass.eq.1) then
           ce(ion) = ce(ion) +&
                1.5_fp_kind*log((isotopic_mass(ielement) - real((nion(ion)),fp_kind)*electron_mass)/isotopic_mass(ielement))
           logqtl_const(nelements+ion) = logqtl_const(nelements+ion) +&
                1.5_fp_kind*log((isotopic_mass(ielement) - real((nion(ion)),fp_kind)*electron_mass)/isotopic_mass(ielement))
        endif
     enddo !do ion = 1,nions
     ! Sanity check
     if(ielement.ne.nelements)&
          error stop 'ionize: nion array values not consistent with nelements ==> logqtl_const array not calculated correctly'
  endif
  ! calculate constant component of entropy.  Keep hydrogen separate
  ! for now, but add in with constant_s (see below) in full ionization
  ! case.
  if(ifsame_abundances.ne.1) then
     ! Occurs on first call (see above check) but not necessarily all subsequent calls.
     constant_sh = eps(1)*logqtl_const(1)
     constant_s = dot_product(eps(2:nelements), logqtl_const(2:nelements))
  endif
  ! ionized_hne is the contribution of fully ionized elements
  ! (depending on ifelement and ifionized) to hne.
  ! ionized_u is the contribution of fully ionized elements
  ! (depending on ifelement and ifionized) to u.
  ! ionized_s is the contribution of fully ionized elements
  ! (depending on ifelement and ifionized) to s
  ! ionized_sum0 is the contribution of fully ionized elements
  ! (depending on ifelement and ifionized) to sum0.
  ! ionized_sum2 is the contribution of fully ionized elements
  ! (depending on ifelement and ifionized) to sum2.
  ! ionized_extrasum(2) is the contribution of fully ionized elements
  ! (depending on ifelement and ifionized) to the last two elements of
  ! extrasum.
  if(ifionized.ne.ifionized_old.or.ifpi.ne.ifpi_old.or.&
       (ifcsums.neqv.ifcsums_old).or.ifsame_abundances.ne.1.or.&
       ifmtrace.ne.ifmtrace_old) then
     ! if different ifionized flag, or if the pressure ionization
     ! has changed, or if the Coulomb approximation has changed,
     ! or if the abundances have changed,
     ! or ifmtrace (ifelement) has changed,
     ! then recalculate full-ionized sum approximations.
     ifionized_old = ifionized
     ifpi_old = ifpi
     ifcsums_old = ifcsums
     ifmtrace_old = ifmtrace

     ionized_hne = 0._fp_kind
     ionized_u = 0._fp_kind
     ionized_s = 0._fp_kind
     if(ifcsums) then
        ionized_sum0 = 0._fp_kind
        ionized_sum2 = 0._fp_kind
     endif
     if(ifpi34) then
        ionized_extrasum(1) = 0._fp_kind
        ionized_extrasum(2) = 0._fp_kind
     endif
     index_f = 0
     do ielement = 1, nelements
        ! index_f should point to last index of ion for this
        ! ielement.
        index_f = index_f + iatomic_number(ielement)
        ! if element treated as fully ionized (or zero abundance)
        ! ifelement(ielement) = 0.  For ifionized = 2, overides
        ! ifelement and all elements are treated as fully ionized.
        if(eps(ielement).gt.0._fp_kind.and.&
             (ifionized.eq.2.or.ifelement(ielement).eq.0)) then
           ionized_hne = ionized_hne +&
                real((iatomic_number(ielement)),fp_kind)*eps(ielement)
           ionized_u = ionized_u + eps(ielement)*bi_ref(index_f)
           ! this expression is limit of expressions below
           ! when all lower ions may be finite but negligible
           ! relative to bare nucleus. n.b. much cancellation!
           ionized_s = ionized_s + eps(ielement)*(ce(index_f) - log(eps(ielement)))
           if(ifcsums) then
              ionized_sum0 = ionized_sum0 +  eps(ielement)
              ionized_sum2 = ionized_sum2 +&
                   real((iatomic_number(ielement)),fp_kind)*real((iatomic_number(ielement)),fp_kind)*eps(ielement)
           endif
           if(ifpi34) then
              ionized_extrasum(1) = ionized_extrasum(1) +&
                   real((iatomic_number(ielement)),fp_kind)*sqrt(real((iatomic_number(ielement)),fp_kind))*eps(ielement)
              ! n.b. there should be *no* contribution
              ! of bare nuclei to last extrasum
           endif
        endif
     enddo !do ielement = 1, nelements
  endif !if(ifionized.ne.ifionized_old.or.ifpi.ne.ifpi_old.or....
  if(ifnr03) then
     hnef = 0._fp_kind
     hnet = 0._fp_kind
     if(ifcsums) then
        sum0f = 0._fp_kind
        sum0t = 0._fp_kind
        sum2f = 0._fp_kind
        sum2t = 0._fp_kind
     endif
     sf = 0._fp_kind
     st = 0._fp_kind
     fionf = 0._fp_kind
     fionhf = 0._fp_kind
  endif
  if(ifpi34) then
     extrasum(1:nextrasum) = 0._fp_kind
     if(ifnr03) then
        extrasumf(1:nextrasum) = 0._fp_kind
        extrasumt(1:nextrasum) = 0._fp_kind
     endif
  endif
  if(ifpl.eq.1) then
     ! n.b. for fully ionized elements, there is no contribution
     ! to Planck-Larkin sums.
     sumpl0 = 0._fp_kind
     sumpl1 = 0._fp_kind
     sumpl2 = 0._fp_kind
     if(ifnr03) then
        sumpl0f = 0._fp_kind
        sumpl1f = 0._fp_kind
        sumpl1t = 0._fp_kind
     endif
  endif
  if(ifionized.eq.2) then
     ! everything (H, He, metals) fully ionized and ifnr irrelevant
     hne = ionized_hne
     u = ionized_u
     s = ionized_s + constant_s + constant_sh
     ! verified as complete ionization limit for detailed expression below
     ! for fion.
     fion = tc2*u - s
     if(ifcsums) then
        sum0 = ionized_sum0*sum0_scale
        sum2 = ionized_sum2*sum2_scale
     endif
     if(ifpi34) then
        extrasum(nextrasum-1) = ionized_extrasum(1)*extrasum_scale(nextrasum-1)
        extrasum(nextrasum) = ionized_extrasum(2)*extrasum_scale(nextrasum)
     endif
  else
     if(ifionized.eq.1) then
        ! trace metals fully ionized, hydrogen effect split off
        hne = ionized_hne
        u = ionized_u
        s = ionized_s + constant_s
        ! verified as limit for detailed expression below
        ! for fion of completely ionized trace metals.
        fion = tc2*u - s
        fionh = 0._fp_kind
        if(ifcsums) then
           sum0 = ionized_sum0*sum0_scale
           sum2 = ionized_sum2*sum2_scale
        endif
        if(ifpi34) then
           extrasum(nextrasum-1) = ionized_extrasum(1)*extrasum_scale(nextrasum-1)
           extrasum(nextrasum) = ionized_extrasum(2)*extrasum_scale(nextrasum)
        endif
     else
        ! ifionized not 2 or 1 so everything treated as partially
        ! ionized, hydrogen effect split off
        hne = 0._fp_kind
        u = 0._fp_kind
        s = constant_s
        ! fion = tc2*u - s
        fion = 0._fp_kind
        fionh = 0._fp_kind
        if(ifcsums) then
           sum0 = 0._fp_kind
           sum2 = 0._fp_kind
        endif
        ! the following commands are redundant so comment out
        !if(ifpi34) then
        !  extrasum(nextrasum-1) = 0.d0
        !  extrasum(nextrasum) = 0.d0
        !endif
     endif

     if(ifnr13) then
        ! Pre-zero all *_dv arrays that are intent(out).
        ! N.B. This zeroing includes the compact molecular indices (if
        ! any) in the upper part of the range of compact indices,
        ! 1:max_index since every quantity calculated by ionize is by
        ! definition independent of molecular formation.
        hne_dv(1:max_index) = 0._fp_kind
        h_ion_dv(1:max_index) = 0._fp_kind
        he_ion_dv(1:max_index) = 0.d0
        he_ion2_dv(1:max_index) = 0.d0
        fionh_dv(1:max_index) = 0._fp_kind
        fion_dv(1:max_index) = 0._fp_kind
        if(ifcsums) then
           sum0_dv(1:max_index) = 0._fp_kind
           sum2_dv(1:max_index) = 0._fp_kind
        endif
        if(ifpi34) extrasum_dv(1:max_index,1:nextrasum) =  0._fp_kind
        if(ifpl.eq.1) sumpl0_dv(1:max_index) = 0._fp_kind
     endif
     nuvar_nelements = n_partial_elements
     do index_element = 1, n_partial_elements
        ielement = partial_elements(index_element)
        epsilon = eps(ielement)
        if(ielement.gt.1) then
           ion0 = ion_end(ielement-1)
        else
           ion0 = 0
        endif
        ! ln fraction of neutral relative to neutral.  So should be
        ! zero other than the consistent shift of dvzero(ielement).
        ! (Such zero point shifts don't matter so long as they are
        ! constant for all monatomic neutral and ionized species of a
        ! given element).
        ln_fractmax = -dvzero(ielement)
        do index_f = 1, iatomic_number(ielement)
           ion = ion0+index_f
           fract(index_f) = ce(ion) + dv(ion) - bi_ref(ion)*tc2
           ln_fractmax = max(ln_fractmax, fract(index_f))
           if(ifnr03) then
              fractf(index_f) = dvf(ion)
              fractt(index_f) = dvt(ion) + bi_ref(ion)*tc2
           endif
        enddo !do index_f = 1, iatomic_number(ielement)

        ! Convert ionization fractions from ln form to actual form
        ! with a normalization near overflow limit exp(arglim) for
        ! maximum accuracy.
        ln_fractmax = ln_fractmax - arglim
        ! evaluate neutral contribution to sum where
        ! fraction = 1, ln fract = 0.
        arg = -ln_fractmax - dvzero(ielement)
        ! normalized so that maximum fraction is exp(arglim) so
        ! that exp(-arglim) or less can safely be ignored.
        ! n.b. note ifsame_under = .false. check on first call above.
        if(ifsame_under) then
           if(ifneutral_zero(ielement)) then
              fractneutral = 0._fp_kind
           else
              fractneutral = exp(arg)
           endif
        else
           if(arg.le.-arglim) then
              fractneutral = 0._fp_kind
              ifneutral_zero(ielement) = .true.
           else
              fractneutral = exp(arg)
              ifneutral_zero(ielement) = .false.
           endif
        endif

        logsum0 = arg
        ! Detailed analysis of code below indicates sumf and sumt are
        ! always used just for the ifnr03 true case and therefore
        ! their initialization for just that case below is sufficient.
        ! Nevertheless, must use unneeded initialization of these
        ! quantities to quiet spurious gfortran
        ! [-Wmaybe-uninitialized] warning.
        sumf = 0._fp_kind
        sumt = 0._fp_kind

        if(ifnr03) then
           sumf = 0._fp_kind
           sumt = 0._fp_kind
           logsum0f = 0._fp_kind
        endif

        ! Because the ifsame_under = .true. case
        ! can (very rarely) zero what ordinarily would be the maximum
        ! fraction, it is important to calculate index_max based on
        ! the fraction logic here rather than the previous ln fraction logic.
        index_max = 0
        fractmax = fractneutral
        ! n.b., the index_maxnb logic below does not get triggered for the
        ! hydrogen case (iatomic_number(ielement) = 1) or used later on,
        ! but that is an accident waiting to happen (in case the subsequent
        ! code is changed) so always initialize index_maxnb here.
        index_maxnb = 0
        do jndex_f = 1,iatomic_number(ielement)
           ! difference in zero point doesn't matter here
           ! and may reduce significance loss.
           arg = fract(jndex_f) - ln_fractmax
           ! normalized so that maximum fraction is exp(arglim) so
           ! that exp(-arglim) or less can safely be ignored.
           ! n.b. note ifsame_under = .false. check on first call above.
           if(ifsame_under) then
              if(ifion_zero(ion0+jndex_f)) then
                 fract(jndex_f) = 0._fp_kind
              else
                 fract(jndex_f) = exp(arg)
              endif
           else
              if(arg.le.-arglim) then
                 fract(jndex_f) = 0._fp_kind
                 ifion_zero(ion0+jndex_f) = .true.
              else
                 fract(jndex_f) = exp(arg)
                 ifion_zero(ion0+jndex_f) = .false.
              endif
           endif
           if(fract(jndex_f).gt.fractmax) then
              index_max = jndex_f
              fractmax = fract(jndex_f)
           endif
           ! this is the maximum component excluding the
           ! bare nucleus (used only for ifpi.eq.3 or 4)
           if(jndex_f.eq.iatomic_number(ielement)-1) index_maxnb = index_max
        enddo !do jndex_f = 1,iatomic_number(ielement)
        if(fractmax.le.0._fp_kind) then
           if(verbosity.ge.1) then
              write(stderr,"(1x,a,i2)") "ionize: neutral and ionized fractions are all zero for ielement = ", ielement
              write(stderr,"(1x,a)") "which constitutes a complete failure of these fraction calculations"
           endif
           info = info_offset_ionize + 1
           return
        endif
        ! note that sum skips maximum component to solve significance
        ! loss problems.
        if(index_max.ne.0) then
           sum = fractneutral
        else
           sum = 0._fp_kind
        endif
        ! Detailed analysis of code below indicates sumnb is
        ! always used just for the ifpi34 true case and therefore
        ! the initialization for just that case below is sufficient.
        ! Nevertheless, must use unneeded initialization of this
        ! quantity to quiet spurious gfortran
        ! [-Wmaybe-uninitialized] warning.
        sumnb = 0._fp_kind
        if(ifpi34) sumnb = 0._fp_kind

        do jndex_f = 1,iatomic_number(ielement)
           ! sum, sumf, sumt skips maximum component to solve
           ! significance loss problems.
           if(jndex_f.ne.index_max) then
              sum = sum + fract(jndex_f)
              if(ifnr03) then
                 sumf = sumf + fract(jndex_f)*fractf(jndex_f)
                 sumt = sumt + fract(jndex_f)*fractt(jndex_f)
              endif
           endif
           if(ifpi34.and.&
                jndex_f.lt.iatomic_number(ielement))&
                sumnb = sumnb + fract(jndex_f)
        enddo !do jndex_f = 1,iatomic_number(ielement)
        if(index_max.ne.0) then
           sum_full = sum + fract(index_max)
        else
           sum_full = sum + fractneutral
        endif
        ! n.b. -logsum0 = log(sum_full/fractneutral)
        ! = log(fract_max/fractneutral) + log(1+sum/fract_max)
        logsum0 = logsum0 - log(sum_full)
        if(ifnr03) then
           if(index_max.ne.0) then
              logsum0f = logsum0f - (sumf + fract(index_max)*fractf(index_max))/sum_full
           else
              logsum0f = logsum0f - sumf/sum_full
           endif
           ! convert to derivatives of ln(sum_full) with missing maximum
           ! component (except where the maximum is the neutral component
           ! in which case the derivative of that component is zero
           ! in any case.)
           sumf = sumf/sum_full
           sumt = sumt/sum_full
        endif
        rpower = fractneutral*epsilon/sum_full
        if(rpower.gt.0._fp_kind.and.iffion) then
           if(ielement.eq.1) then
              fionh = fionh + rpower*(log(rpower) - logqtl_const(ielement))
           else
              fion = fion + rpower*(log(rpower) - logqtl_const(ielement))
           endif
        endif
        ! in all cases and all indices save n/(rho*avogadro) in nuvar
        ! where nuvar is used by eos_calc and excitation_sum.
        nuvar(1,index_element) = rpower
        if(ielement.eq.1.or.(rpower.gt.0._fp_kind.and.(iffion.or.ifpl.eq.1.or.ifpi34))) then
           ! Detailed analysis of code below indicates rpowerf and
           ! rpowert are always used just for the ifnr03 true case and
           ! therefore their initialization for just that case below
           ! is sufficient.  Nevertheless, must use unneeded
           ! initialization of these quantities to quiet spurious
           ! gfortran [-Wmaybe-uninitialized] warning.
           rpowerf = 0._fp_kind
           rpowert = 0._fp_kind

           if(ifnr03) then
              if(index_max.ne.0) then
                 ! sumf, sumt skip maximum component.
                 rpowerf = -rpower*(sumf + fract(index_max)*fractf(index_max)/sum_full)
                 rpowert = -rpower*(sumt + fract(index_max)*fractt(index_max)/sum_full)
              else
                 ! fractf, fractt = 0, and sumf, sumt complete, for
                 ! index_max = 0.
                 rpowerf = -rpower*sumf
                 rpowert = -rpower*sumt
              endif
           endif
           if(ifnr13.and.(iffion.or.ifpl.eq.1.or.ifpi34)) then
              ! n.b. d *ln* rpower/d dv
              rpower_dv(1:iatomic_number(ielement)) = -fract(1:iatomic_number(ielement))/sum_full
           endif
           if(rpower.gt.0._fp_kind.and.iffion) then
              if(ielement.eq.1) then
                 if(ifnr03) then
                    fionhf = fionhf + rpowerf*(log(rpower) - logqtl_const(ielement) + 1._fp_kind)
                 endif
                 if(ifnr13) then
                    fionh_dv(1:iatomic_number(ielement)) = fionh_dv(1:iatomic_number(ielement)) +&
                         rpower*rpower_dv(1:iatomic_number(ielement))*(&
                         log(rpower) - logqtl_const(ielement) + 1._fp_kind)
                 endif
              else
                 if(ifnr03) then
                    fionf = fionf + rpowerf*(log(rpower) - logqtl_const(ielement) + 1._fp_kind)
                 endif
                 if(ifnr13) then
                    do index_f = 1, iatomic_number(ielement)
                       inv_ion_index = inv_ion(index_f+ion0)
                       fion_dv(inv_ion_index) = fion_dv(inv_ion_index) +&
                            rpower*rpower_dv(index_f)*(log(rpower) - logqtl_const(ielement) + 1._fp_kind)
                    enddo
                 endif
              endif
           endif

           if(ifexcited.gt.0) then
              if(ifnr03) then
                 nuvarf(1,index_element) = rpowerf
                 nuvart(1,index_element) = rpowert
              endif
              if(ifnr13) then
                 nuvar_dv(1:iatomic_number(ielement),1,index_element) = rpower*rpower_dv(1:iatomic_number(ielement))
              endif
           endif
           if(ielement.eq.1) then
              ! save special quantities for hydrogen.
              h_neutral = rpower
              h_ion_equil = ce(1) + dv(1) - bi_ref(1)*tc2
              if(ifnr03) then
                 h_neutralf = rpowerf
                 h_neutralt = rpowert
                 h_ion_equilf = dvf(1)
                 h_ion_equilt = dvt(1) + bi_ref(1)*tc2
              endif
           else
              if(ifpl.eq.1) then
                 ! Planck-Larkin occupation probability sums
                 sumpl0 = sumpl0 + rpower*plop(ion0+1)
                 sumpl1 = sumpl1 + rpower*(plop(ion0+1)+plopt(ion0+1))
                 sumpl2 = sumpl2 + rpower*plopt(ion0+1)
                 if(ifnr03) then
                    sumpl0f = sumpl0f + rpowerf*plop(ion0+1)
                    sumpl1f = sumpl1f + rpowerf*(plop(ion0+1)+plopt(ion0+1))
                    sumpl1t = sumpl1t + rpowert*(plop(ion0+1)+plopt(ion0+1)) +&
                         rpower*(plopt(ion0+1)+plopt2(ion0+1))
                 endif
                 if(ifnr13) then
                    do index_f = 1,iatomic_number(ielement)
                       inv_ion_index = inv_ion(index_f+ion0)
                       sumpl0_dv(inv_ion_index) =&
                            sumpl0_dv(inv_ion_index) +&
                            rpower_dv(index_f)*rpower*plop(ion0+1)
                    enddo
                 endif
              endif
              if(ifpi34) then
                 extrasum(nextrasum) = extrasum(nextrasum) +&
                      rpower*extrasum_scale(nextrasum)*r_ion3(ion0+1)
                 if(ifnr03) then
                    extrasumf(nextrasum) = extrasumf(nextrasum) +&
                         rpowerf*extrasum_scale(nextrasum)*r_ion3(ion0+1)
                    extrasumt(nextrasum) = extrasumt(nextrasum) +&
                         rpowert*extrasum_scale(nextrasum)*r_ion3(ion0+1)
                 endif
                 if(ifnr13) then
                    do index_f = 1,iatomic_number(ielement)
                       inv_ion_index = inv_ion(index_f+ion0)
                       extrasum_dv(inv_ion_index,nextrasum) =&
                            extrasum_dv(inv_ion_index,nextrasum) +&
                            rpower_dv(index_f)*rpower*&
                            extrasum_scale(nextrasum)*r_ion3(ion0+1)
                    enddo
                 endif
                 do iextrasum = 1,nextrasum-2
                    extrasum(iextrasum) = extrasum(iextrasum) +&
                         extrasum_scale(iextrasum)*rpower
                    if(ifnr03) then
                       extrasumf(iextrasum) = extrasumf(iextrasum) +&
                            extrasum_scale(iextrasum)*rpowerf
                       extrasumt(iextrasum) = extrasumt(iextrasum) +&
                            extrasum_scale(iextrasum)*rpowert
                       rpowerf = rpowerf*r_neutral(ielement)
                       rpowert = rpowert*r_neutral(ielement)
                    endif
                    if(ifnr13) then
                       do index_f = 1,iatomic_number(ielement)
                          inv_ion_index = inv_ion(index_f+ion0)
                          extrasum_dv(inv_ion_index,iextrasum) =&
                               extrasum_dv(inv_ion_index,iextrasum) +&
                               rpower_dv(index_f)*&
                               extrasum_scale(iextrasum)*rpower
                       enddo
                    endif
                    rpower = rpower*r_neutral(ielement)
                 enddo !do iextrasum = 1,nextrasum-2
              endif
           endif
        else !if(ielement.eq.1.or.(rpower.gt.0._fp_kind.and.(iffion.or.ifpl.eq.1.or.ifpi34))) then
           ! nuvar(1,index_element) always calculated (see above).
           if(ifexcited.gt.0) then
              if(ifnr03) then
                 nuvarf(1,index_element) = 0._fp_kind
                 nuvart(1,index_element) = 0._fp_kind
              endif
              if(ifnr13) then
                 nuvar_dv(1:iatomic_number(ielement),1,index_element) = 0._fp_kind
              endif
           endif
        endif !if(ielement.eq.1.or.(rpower.gt.0._fp_kind.and.(iffion.or.ifpl.eq.1.or.ifpi34))) then
        ! Detailed analysis of code below indicates sarg0 is
        ! always used just for the index_max.ne.0 case and therefore
        ! the initialization for just that case below is sufficient.
        ! Nevertheless, must use unneeded initialization of this
        ! quantity to quiet spurious gfortran
        ! [-Wmaybe-uninitialized] warning.
        sarg0 = 0._fp_kind
        if(index_max.ne.0) then
           sarg0 = bi_ref(ion0+index_max)*tc2-dv(ion0+index_max)
        endif
        nuvar_atomic_number(index_element) = iatomic_number(ielement)
        nuvar_index_element(index_element) = ielement
        do jndex_f = 1,iatomic_number(ielement)
           if(ielement.eq.1.or.(fract(jndex_f).gt.0._fp_kind)) then
              ! sum_full is the sum of fract values for this element so
              ! fract(jndex_f)/sum_full is the fraction of the number
              ! density of this element = the fraction of
              ! (rho*Navogadro)*epsilon in this particular form.  Thus,
              ! vz = number density/(rho*Navogadro) of this particular form
              ! of the element, and sum of vz over all forms of the element
              ! is equal to epsilon.
              vz = fract(jndex_f)*epsilon/sum_full
              ! Detailed analysis of code below indicates vzf and vzt
              ! are always used just for the ifnr03 true case and
              ! therefore their initialization for just that case
              ! below is sufficient.  Nevertheless, must use unneeded
              ! initialization of these quantities to quiet spurious
              ! gfortran [-Wmaybe-uninitialized] warning.
              vzf = 0._fp_kind
              vzt = 0._fp_kind

              if(ifnr03) then
                 ! sumf and sumt skip the maximum component
                 ! expressions below have negligible significance loss.
                 if(index_max.eq.jndex_f) then
                    vzf = vz*(fractf(jndex_f)*sum/sum_full-sumf)
                    vzt = vz*(fractt(jndex_f)*sum/sum_full-sumt)
                 elseif(index_max.ne.0) then
                    vzf = vz*(fractf(jndex_f) - sumf -&
                         fract(index_max)*fractf(index_max)/sum_full)
                    vzt = vz*(fractt(jndex_f) - sumt -&
                         fract(index_max)*fractt(index_max)/sum_full)
                 else
                    ! for index_max = 0, sumf and sumt are complete.
                    vzf = vz*(fractf(jndex_f)-sumf)
                    vzt = vz*(fractt(jndex_f)-sumt)
                 endif
              endif
              if(ifnr13) then
                 do index_f = 1, iatomic_number(ielement)
                    if(index_f.ne.jndex_f) then
                       vz_dv(index_f) = -vz*fract(index_f)/sum_full
                    endif
                 enddo
                 if(index_max.eq.jndex_f) then
                    vz_dv(jndex_f) = vz*sum/sum_full
                 else
                    vz_dv(jndex_f) = vz*&
                         (1._fp_kind - fract(jndex_f)/sum_full)
                 endif
              endif
              if(iffion.and.vz.gt.0._fp_kind) then
                 if(ielement.eq.1) then
                    fionh = fionh + vz*(&
                         log(vz) - logqtl_const(nelements+ion0+jndex_f) + tc2*bi_ref(ion0+jndex_f))
                    if(ifnr03) then
                       fionhf = fionhf + vzf*(&
                            log(vz) - logqtl_const(nelements+ion0+jndex_f) + tc2*bi_ref(ion0+jndex_f) + 1._fp_kind)
                    endif
                    if(ifnr13) then
                       fionh_dv(1:iatomic_number(ielement)) = fionh_dv(1:iatomic_number(ielement)) +&
                            vz_dv(1:iatomic_number(ielement))*(&
                            log(vz) - logqtl_const(nelements+ion0+jndex_f) + tc2*bi_ref(ion0+jndex_f) + 1._fp_kind)
                    endif
                 else
                    fion = fion + vz*(&
                         log(vz) - logqtl_const(nelements+ion0+jndex_f) + tc2*bi_ref(ion0+jndex_f))
                    if(ifnr03) then
                       fionf = fionf + vzf*(&
                            log(vz) - logqtl_const(nelements+ion0+jndex_f) + tc2*bi_ref(ion0+jndex_f) + 1._fp_kind)
                    endif
                    if(ifnr13) then
                       do index_f = 1, iatomic_number(ielement)
                          inv_ion_index = inv_ion(index_f+ion0)
                          fion_dv(inv_ion_index) = fion_dv(inv_ion_index) +&
                               vz_dv(index_f)*(&
                               log(vz) - logqtl_const(nelements+ion0+jndex_f) + tc2*bi_ref(ion0+jndex_f) + 1._fp_kind)
                       enddo
                    endif
                 endif
              endif
              ! in all cases and all indices save vz in nuvar where
              ! nuvar is used by eos_calc and excitation_sum.
              nuvar(jndex_f+1,index_element) = vz
              ! save derivatives as well for excitation and for all but
              ! bare ion index range.
              if(ifexcited.gt.0.and.&
                   jndex_f.lt.iatomic_number(ielement)) then
                 if(ifnr03) then
                    nuvarf(jndex_f+1,index_element) = vzf
                    nuvart(jndex_f+1,index_element) = vzt
                 endif
                 if(ifnr13) then
                    nuvar_dv(1:iatomic_number(ielement),jndex_f+1,index_element) = vz_dv(1:iatomic_number(ielement))
                 endif
              endif
              if(ielement.eq.1) then
                 ! save special quantities for hydrogen.
                 ! note, jndex_f is 1.
                 sarg = bi_ref(ion0+1)*tc2-dv(ion0+1)
                 if(index_max.eq.0) then
                    ! use regular expression.
                    sh = vz*sarg - epsilon*(log(epsilon) + logsum0)
                 else
                    ! when h_ion is maximum component of sum_full, then
                    ! cancellations occur which cause significance loss
                    ! unless you use this revised expression.
                    ! n.b. -logsum0 = log(sum_full/fractneutral)
                    !               = log(fract_max/fractneutral) + log(1+sum/fract_max)
                    sh = epsilon*(-log(epsilon) + ce(ion0+index_max) +&
                         log(1._fp_kind+sum/fract(index_max)) - (sum/sum_full)*&
                         (bi_ref(ion0+index_max)*tc2-dv(ion0+index_max)))
                 endif
                 sh = sh + constant_sh
                 h_ion = vz
                 if(ifnr03) then
                    shf = vzf*sarg
                    sht = vzt*sarg
                    h_ionf = vzf
                    h_iont = vzt
                 endif
                 if(ifnr13) then
                    inv_ion_index = inv_ion(ion0+1)
                    h_ion_dv(inv_ion_index) = vz_dv(1)
                 endif
              else
                 ! save special quanties for helium.
                 if(ielement.eq.2.and.jndex_f.eq.1) then
                    he_ion = vz
                    if(ifnr03) then
                       he_ionf = vzf
                       he_iont = vzt
                    endif
                    if(ifnr13) then
                       inv_ion_index = inv_ion(2)
                       he_ion_dv(inv_ion_index) = vz_dv(1)
                       inv_ion_index = inv_ion(3)
                       he_ion_dv(inv_ion_index) = vz_dv(2)
                    endif
                 elseif(ielement.eq.2.and.jndex_f.eq.2) then
                    he_ion2 = vz
                    if(ifnr03) then
                       he_ion2f = vzf
                       he_ion2t = vzt
                    endif
                    if(ifnr13) then
                       inv_ion_index = inv_ion(2)
                       he_ion2_dv(inv_ion_index) = vz_dv(1)
                       inv_ion_index = inv_ion(3)
                       he_ion2_dv(inv_ion_index) = vz_dv(2)
                    endif
                 endif
                 ion_index = ion0+jndex_f
                 u = u + bi_ref(ion_index)*vz
                 ! sarg should have dvzero(ielement) subtracted, but
                 ! this split off to avoid significance loss
                 ! for s, sf, and st
                 ! for same reason subtract off sarg of maximum
                 ! component to be added later.
                 if(index_max.eq.0) then
                    sarg = bi_ref(ion_index)*tc2-dv(ion_index)
                 elseif(jndex_f.ne.index_max) then
                    sarg = bi_ref(ion_index)*tc2-dv(ion_index) - sarg0
                 else
                    ! sarg = bi_ref(ion_index)*tc2-dv(ion_index) - sarg0
                    sarg = 0._fp_kind
                 endif
                 ! ignore entropy term due to maximum fraction.  this term
                 ! will be added later in a way which avoids significance
                 ! loss.  the derivatives are unaffected by this significance
                 ! loss so they are done without the complications.
                 ! n.b. we also ignore *complete* sum over vz*(sarg0-dvzero)
                 ! which will be added later.
                 if(jndex_f.ne.index_max) then
                    s = s + vz*sarg
                 endif
                 if(ifcsums) then
                    ! accumulate number density of positive ions/
                    ! (rho/H = rho*Navogadro)
                    ! do the following sum analytically to avoid
                    ! significance loss. also subtract out constant
                    ! component of sum2.
                    ! sum0 = sum0 + vz
                    ! accumulate number density of positive ions *
                    ! charge^2/(rho/H = rho*Navogadro)
                    sum2 = sum2 +&
                         real((jndex_f*jndex_f-index_max*index_max),fp_kind)*&
                         sum2_scale*vz
                    if(ifnr03) then
                       ! see above comment.
                       ! sum0f = sum0f + vzf
                       ! sum0t = sum0t + vzt
                       sum2f = sum2f +&
                            real((jndex_f*jndex_f-index_max*index_max),fp_kind)*&
                            sum2_scale*vzf
                       sum2t = sum2t +&
                            real((jndex_f*jndex_f-index_max*index_max),fp_kind)*&
                            sum2_scale*vzt
                    endif
                 endif
                 ! accumulate number density of positive charges/
                 ! (rho/H = rho*Navogadro)
                 ! add in constant component later except for
                 ! index_max = 0 case.
                 hne = hne + real((jndex_f-index_max),fp_kind)*vz
                 if(ifnr03) then
                    hnef = hnef + real((jndex_f-index_max),fp_kind)*vzf
                    hnet = hnet + real((jndex_f-index_max),fp_kind)*vzt
                 endif
                 if((ifpl.eq.1).and.&
                      jndex_f.lt.iatomic_number(ielement)) then
                    ! Planck-Larkin occupation probability sums
                    ! add into sum if not bare nucleus.
                    sumpl0 = sumpl0 + vz*plop(ion_index+1)
                    sumpl1 = sumpl1 + vz*(plop(ion_index+1) +&
                         plopt(ion_index+1))
                    sumpl2 = sumpl2 + vz*plopt(ion_index+1)
                    if(ifnr03) then
                       sumpl0f = sumpl0f + vzf*plop(ion_index+1)
                       sumpl1f = sumpl1f +&
                            vzf*(plop(ion_index+1)+plopt(ion_index+1))
                       sumpl1t = sumpl1t +&
                            vzt*(plop(ion_index+1)+plopt(ion_index+1)) +&
                            vz*(plopt(ion_index+1)+plopt2(ion_index+1))
                    endif
                    if(ifnr13) then
                       do index_f = 1,iatomic_number(ielement)
                          inv_ion_index = inv_ion(index_f+ion0)
                          sumpl0_dv(inv_ion_index) =&
                               sumpl0_dv(inv_ion_index) +&
                               vz_dv(index_f)*plop(ion_index+1)
                       enddo
                    endif
                 endif
                 if(ifpi34) then
                    rpower = real((jndex_f),fp_kind)*sqrt(real((jndex_f),fp_kind))
                    extrasum(nextrasum-1) = extrasum(nextrasum-1) +&
                         rpower*extrasum_scale(nextrasum-1)*vz
                    if(ifnr03) then
                       extrasumf(nextrasum-1) = extrasumf(nextrasum-1) +&
                            rpower*extrasum_scale(nextrasum-1)*vzf
                       extrasumt(nextrasum-1) = extrasumt(nextrasum-1) +&
                            rpower*extrasum_scale(nextrasum-1)*vzt
                    endif
                    if(ifnr13) then
                       do index_f = 1,iatomic_number(ielement)
                          inv_ion_index = inv_ion(index_f+ion0)
                          extrasum_dv(inv_ion_index,nextrasum-1) =&
                               extrasum_dv(inv_ion_index,nextrasum-1) +&
                               rpower*&
                               extrasum_scale(nextrasum-1)*vz_dv(index_f)
                       enddo
                    endif
                    ! the ion radii indices are offset by one, and the bare
                    ! nucleii should be skipped.
                    ! n.b. the if statement means that hydrogen is always
                    ! skipped, and eos_calc must correct for that.
                    if(jndex_f.lt.iatomic_number(ielement)) then
                       if(index_maxnb.eq.0) then
                          rpower = r_ion3(ion_index+1)
                       elseif(jndex_f.ne.index_maxnb) then
                          rpower = r_ion3(ion_index+1) - r_ion3(ion0+index_maxnb+1)
                       else
                          ! offset value.
                          rpower = r_ion3(ion0+index_maxnb+1)
                       endif
                       if(jndex_f.ne.index_maxnb) then
                          ! n.b. this if block includes
                          ! the case where index_maxnb = 0,
                          ! and no offset value is needed or calculated.
                          extrasum(nextrasum) = extrasum(nextrasum) +&
                               rpower*extrasum_scale(nextrasum)*vz
                          if(ifnr03) then
                             extrasumf(nextrasum) = extrasumf(nextrasum) + rpower*extrasum_scale(nextrasum)*vzf
                             extrasumt(nextrasum) = extrasumt(nextrasum) + rpower*extrasum_scale(nextrasum)*vzt
                          endif
                          if(ifnr13) then
                             do index_f = 1,iatomic_number(ielement)
                                inv_ion_index = inv_ion(index_f+ion0)
                                extrasum_dv(inv_ion_index,nextrasum) = extrasum_dv(inv_ion_index,nextrasum) +&
                                     rpower*extrasum_scale(nextrasum)*vz_dv(index_f)
                             enddo
                          endif
                       else
                          ! offset maximum component is zero.
                          ! maximum component index
                          ! encountered once (if at all) per element.
                          ! for this index add in
                          ! sum over all ions of offset.
                          extrasum(nextrasum) = extrasum(nextrasum) + rpower*epsilon*sumnb*&
                               (extrasum_scale(nextrasum)/sum_full)
                          ! to obtain derivatives take sumnb/sum_full =
                          ! 1-(fractneutral+fract(iatomic_number(ielement)))/
                          ! sum_full
                          if(ifnr03) then
                             ! sumf and sumt skip the maximum component
                             ! expressions below have negligible significance loss.
                             if(index_max.eq.iatomic_number(ielement)) then
                                extrasumf(nextrasum) =&
                                     extrasumf(nextrasum) -&
                                     rpower*epsilon*&
                                     (extrasum_scale(nextrasum)/sum_full)*&
                                     (fract(iatomic_number(ielement))*&
                                     (fractf(iatomic_number(ielement))*&
                                     sum/sum_full-sumf) - fractneutral*&
                                     (sumf + fract(index_max)*&
                                     fractf(index_max)/sum_full))
                                extrasumt(nextrasum) =&
                                     extrasumt(nextrasum) -&
                                     rpower*epsilon*&
                                     (extrasum_scale(nextrasum)/sum_full)*&
                                     (fract(iatomic_number(ielement))*&
                                     (fractt(iatomic_number(ielement))*&
                                     sum/sum_full-sumt) - fractneutral*&
                                     (sumt + fract(index_max)*&
                                     fractt(index_max)/sum_full))
                             else
                                extrasumf(nextrasum) =&
                                     extrasumf(nextrasum) -&
                                     rpower*epsilon*&
                                     (extrasum_scale(nextrasum)/sum_full)*&
                                     (fract(iatomic_number(ielement))*&
                                     fractf(iatomic_number(ielement)) -&
                                     (fractneutral +&
                                     fract(iatomic_number(ielement)))*&
                                     (sumf + fract(index_max)*&
                                     fractf(index_max)/sum_full))
                                extrasumt(nextrasum) =&
                                     extrasumt(nextrasum) -&
                                     rpower*epsilon*&
                                     (extrasum_scale(nextrasum)/sum_full)*&
                                     (fract(iatomic_number(ielement))*&
                                     fractt(iatomic_number(ielement)) -&
                                     (fractneutral +&
                                     fract(iatomic_number(ielement)))*&
                                     (sumt + fract(index_max)*&
                                     fractt(index_max)/sum_full))
                             endif
                          endif
                          if(ifnr13) then
                             do index_f = 1, iatomic_number(ielement)-1
                                inv_ion_index = inv_ion(index_f+ion0)
                                extrasum_dv(inv_ion_index,nextrasum) =&
                                     extrasum_dv(inv_ion_index,nextrasum) +&
                                     rpower*epsilon*&
                                     (extrasum_scale(nextrasum)/sum_full)*&
                                     (fractneutral +&
                                     fract(iatomic_number(ielement)))*&
                                     (fract(index_f)/sum_full)
                             enddo
                             inv_ion_index =&
                                  inv_ion(iatomic_number(ielement)+ion0)
                             extrasum_dv(inv_ion_index,nextrasum) =&
                                  extrasum_dv(inv_ion_index,nextrasum) -&
                                  rpower*epsilon*sumnb*&
                                  (extrasum_scale(nextrasum)/sum_full)*&
                                  (fract(iatomic_number(ielement))/sum_full)
                          endif
                       endif
                    endif
                 endif !if(ifpi34) then
                 ! naively, this expression doesn't seem right, but have
                 ! used constancy of sum of vz over all species of helium and
                 ! metals to derive this expression.
                 if(ifnr03) then
                    ! skip sarg0 - dvzero(ielement) term which will
                    ! be added in later in a way that reduces
                    ! significance loss
                    sf = sf + vzf*sarg
                    st = st + vzt*sarg
                 endif
                 if(ifnr13) then
                    ! to avoid significance loss use special
                    ! pre-summed form for sum0_dv.
                    if(ifcsums) then
                       inv_ion_index = inv_ion(ion0+jndex_f)
                       sum0_dv(inv_ion_index) =&
                            vz*fractneutral*(sum0_scale/sum_full)
                       sum2_dv(inv_ion_index) =&
                            sum2_dv(inv_ion_index) +&
                            real((index_max*index_max),fp_kind)*&
                            vz*fractneutral*(sum2_scale/sum_full)
                    endif
                    do index_f = 1, iatomic_number(ielement)
                       inv_ion_index = inv_ion(ion0+index_f)
                       if(ifcsums) then
                          ! avoid the following commented out expression
                          ! because of bad significance loss.
                          !sum0_dv(inv_ion_index) = sum0_dv(inv_ion_index) + vz_dv(index_f)
                          sum2_dv(inv_ion_index) =&
                               sum2_dv(inv_ion_index) +&
                               real((jndex_f*jndex_f-index_max*index_max),fp_kind)*&
                               sum2_scale*vz_dv(index_f)
                       endif
                       hne_dv(inv_ion_index) = hne_dv(inv_ion_index) +&
                            real((jndex_f-index_max),fp_kind)*vz_dv(index_f)
                    enddo
                    ! use sum rule for vz_dv for constant part to
                    ! avoid significance loss.
                    inv_ion_index = inv_ion(ion0+jndex_f)
                    hne_dv(inv_ion_index) = hne_dv(inv_ion_index) +&
                         real((index_max),fp_kind)*vz*fractneutral/sum_full
                 endif
              endif
           else !if(ielement.eq.1.or.(fract(jndex_f).gt.0._fp_kind)) then
              ! save special quantities for helium.
              if(ielement.eq.2.and.jndex_f.eq.1) then
                 he_ion = 0._fp_kind
                 if(ifnr03) then
                    he_ionf = 0._fp_kind
                    he_iont = 0._fp_kind
                 endif
                 if(ifnr13) then
                    inv_ion_index = inv_ion(2)
                    he_ion_dv(inv_ion_index) = 0._fp_kind
                    inv_ion_index = inv_ion(3)
                    he_ion_dv(inv_ion_index) = 0._fp_kind
                 endif
              elseif(ielement.eq.2.and.jndex_f.eq.2) then
                 he_ion2 = 0._fp_kind
                 if(ifnr03) then
                    he_ion2f = 0._fp_kind
                    he_ion2t = 0._fp_kind
                 endif
                 if(ifnr13) then
                    inv_ion_index = inv_ion(2)
                    he_ion2_dv(inv_ion_index) = 0._fp_kind
                    inv_ion_index = inv_ion(3)
                    he_ion2_dv(inv_ion_index) = 0._fp_kind
                 endif
              endif
              nuvar(jndex_f+1,index_element) = 0._fp_kind
              if(ifexcited.gt.0.and.jndex_f.lt.iatomic_number(ielement)) then
                 if(ifnr03) then
                    nuvarf(jndex_f+1,index_element) = 0._fp_kind
                    nuvart(jndex_f+1,index_element) = 0._fp_kind
                 endif
                 if(ifnr13) then
                    nuvar_dv(1:iatomic_number(ielement),jndex_f+1,index_element) = 0._fp_kind
                 endif
              endif
           endif !if(ielement.eq.1.or.(fract(jndex_f).gt.0._fp_kind)) then
        enddo !do jndex_f = 1,iatomic_number(ielement)
        ! have used constancy of sum of vz over all ionization species
        ! to derive this expression.
        ! n.b. there are no corresponding additions to sf and st, see above.
        if(ielement.gt.1.and.epsilon.gt.0._fp_kind) then
           if(index_max.eq.0) then
              ! maximum term not dropped from entropy sum above so
              ! use regular expression.
              s = s - epsilon*(log(epsilon) + logsum0)
              ! correct for dvzero term that is missing from *complete*
              ! sum over dv.
              s = s - epsilon*(sum/sum_full)*dvzero(ielement)
              if(ifcsums) sum0 = sum0 + epsilon*sum*(sum0_scale/sum_full)
           else
              ! maximum term was dropped from entropy sum above so
              ! use regular expression corrected by this term in
              ! a way that avoids significance loss.
              ! n.b. -logsum0 = log(sum_full/fractneutral) = log(fract_max/fractneutral) + log(1+sum/fract_max)
              s = s + (epsilon*(-log(epsilon) + ce(ion0+index_max) +&
                   log(1._fp_kind+sum/fract(index_max)) - (sum/sum_full)*&
                   (bi_ref(ion0+index_max)*tc2 - dv(ion0+index_max))))
              ! correct for dvzero term that is missing from *complete*
              ! sum over dv.  the -epsilon*dvzero term has already been
              ! analytically subtracted without significance loss
              s = s + epsilon*(fractneutral/sum_full)*dvzero(ielement)
              ! correct for sarg0 term that must be added for all vz
              ! except maximum term
              s = s + epsilon*((sum - fractneutral)/sum_full)*sarg0
              hne = hne +&
                   epsilon*(1._fp_kind - fractneutral/sum_full)*real((index_max),fp_kind)
              if(ifcsums) then
                 sum0 = sum0 + epsilon*sum0_scale*(1._fp_kind - fractneutral/sum_full)
                 sum2 = sum2 + epsilon*sum2_scale*(1._fp_kind - fractneutral/sum_full)*real((index_max*index_max),fp_kind)
              endif
           endif
           if(ifnr03) then
              ! add in sum over vzf and vzt times -dvzero(ielement) term.
              ! these sum derivatives are the negative of the vzneutral
              ! derivatives which I derive from the vzf and vzt expressions
              ! above with vzneutral = fractneutral*epsilon/sum_full
              ! replacing vz.
              ! note from outside logic that dvzero is non-zero only
              ! when fractneutral is small.
              if(index_max.eq.0) then
                 sf = sf + sumf*(fractneutral*epsilon/sum_full)*(-dvzero(ielement))
                 st = st + sumt*(fractneutral*epsilon/sum_full)*(-dvzero(ielement))
                 if(ifcsums) then
                    sum0f = sum0f + sumf*epsilon*fractneutral*(sum0_scale/sum_full)
                    sum0t = sum0t + sumt*epsilon*fractneutral*(sum0_scale/sum_full)
                 endif
              else
                 sf = sf +&
                      (sumf + fract(index_max)*fractf(index_max)/sum_full)*&
                      (fractneutral*epsilon/sum_full)*(sarg0-dvzero(ielement))
                 st = st +&
                      (sumt + fract(index_max)*fractt(index_max)/sum_full)*&
                      (fractneutral*epsilon/sum_full)*(sarg0-dvzero(ielement))
                 hnef = hnef +&
                      (sumf + fract(index_max)*fractf(index_max)/sum_full)*&
                      (fractneutral*epsilon/sum_full)*real((index_max),fp_kind)
                 hnet = hnet +&
                      (sumt + fract(index_max)*fractt(index_max)/sum_full)*&
                      (fractneutral*epsilon/sum_full)*real((index_max),fp_kind)
                 if(ifcsums) then
                    sum0f = sum0f +&
                         (sumf + fract(index_max)*fractf(index_max)/sum_full)*epsilon*&
                         fractneutral*(sum0_scale/sum_full)
                    sum0t = sum0t +&
                         (sumt + fract(index_max)*fractt(index_max)/sum_full)*epsilon*&
                         fractneutral*(sum0_scale/sum_full)
                    sum2f = sum2f +&
                         (sumf + fract(index_max)*fractf(index_max)/sum_full)*epsilon*&
                         fractneutral*(sum2_scale/sum_full)*real((index_max*index_max),fp_kind)
                    sum2t = sum2t +&
                         (sumt + fract(index_max)*fractt(index_max)/sum_full)*epsilon*&
                         fractneutral*(sum2_scale/sum_full)*real((index_max*index_max),fp_kind)
                 endif
              endif
           endif
        endif !if(ielement.gt.1.and.epsilon.gt.0._fp_kind) then
     enddo !do index_element = 1, n_partial_elements
  endif
end subroutine ionize
