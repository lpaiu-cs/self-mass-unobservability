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

! the purpose of this routine is to calculate the change in equilibrium
! constant (relative to ground state monatomic and free electrons)
! from the appropriate combinations of pressure ionization
! free energy/kT.
! For this case, the MDH-like free energy/unit volume model given by:
! f = kt*4 pi/3 *[
!   2(3 alpha1*alpha2 + alpha0*alpha3) +
!   gamma*beta] + f'
! f' is an ad hoc extra term designed (MDH paper II, Appendix B) to
!   force pressure ionization for low temperatures.  It has the
!   right functional form for the next higher order interaction.
!   if one expands out the various powers of r in the MDH expression, then
! f' = quad*kt*(4 pi/3)^2 *[
!   3 alpha0*alpha3*alpha3 +
!   30 alpha1*alpha2*alpha3 +
!   9 alpha2*alpha2*alpha2 +
!   6 alpha0*alpha2*alpha4 +
!   9 alpha1*alpha1*alpha4 +
!   6 alpha0*alpha1*alpha5 +
!   1 alpha0*alpha0*alpha6]
!
! alphak = sum(over all neutrals + H2+) n(i)* r_neutral(i)^k
! beta = sum(over all ions including bare nucleii, but excluding ne)
!   n(i)*nion(i)^1.5.
! gamma = sum (all neutral and ionized species except ne, and bare nucleii) of:
!   n_i * rion(i)^3.
! n_i and rion are in *shifted* ion order, e.g., h, he, he+, ...,
!   where each element starts at neutral and ends with next to
!   bare nucleii.
! n.b. extrasum(1-->nextrasum-2) = alpha0-->alpha{nextrasum-3}.
!   extrasum(nextrasum-1) = beta, and
!   extrasum(nextrasum) = gamma.
! n.b. MDH paper II introduced quad term to force pressure ionization for
!   low T, high densities.  The quad term corresponds to alpha of
!   the MDH paper which was set to 10 for their computations.
! input:
! \param ifdvzero ifdvzero = .true. implies dvzero (the quantity
!   added to all dv in ionize) is zero.  This is the low-temperature
!   option.  For high temperatures, ifdvzero is .false., and dvzero
!   is a zero point shift that renders the pressure ionization
!   contribution to dv of all bare ions zero.  This option gives the
!   smallest significance loss near full ionization and high
!   densities.<br>
! ifnr = 0, calculate fl and tl derivatives of output using chain rule
!   and auxf and auxt
! ifnr = 1, calculate aux derivatives of output with fl, tl fixed
! ifnr = 2, calculate fl and tl derivatives of output with aux fixed.
! ifnr = 3 is combination of ifnr = 1 and ifnr = 2.
! inv_aux(naux) is pointer to actual auxiliary variable index
! iextraoff is index offset between extrasum and aux.
! naux is number of active auxiliary variables
! inv_ion(nions+2) maps ion index to contiguous ion index used in NR
!   iteration dot products.
! partial_elements(n_partial_elements+2) index of elements treated as
!   partially ionized consistent with ifelement.
! ion_end(mion_end) keeps track of largest ion index for each element
!   and largest ion index for H2 (ielement = nelements+1 when
!   H2 is treated in the EOS and H2+ (ielement = nelements+2 when
!   H2+ is treated in the EOS.)
! ifmodified > 0 means quad terms modified.
! r_ion3(nions+2), *the cube of the*
!   effective radii (in ion order but going from neutral
!   to next to bare ion for each species) of MDH interaction between
!   all but bare nucleii species and ionized species.  last 2 are
!   H2 and H2+
! nion(nions), charge on ion in ion order (must be same order as bi)
!   e.g., for H+, He+, He++, etc.
! r_neutral(nelements+2) effective radii for MDH neutral-neutral
!   interactions.  Last two are H2 and H2+ (the only ionic
!   species in the MHD model with a non-zero hard-sphere radius).
! extrasum(nextrasum = 9) weighted sums over n(i)
!   for iextrasum = 1,nextrasum-2,
!   sum is only over neutral species (and H2+) and
!   weight is r_neutral^{iextrasum-1}
!   for iextrasum = nextrasum-1, sum is over all ionized species
!   including bare nucleii, but excluding free electrons, weight is
!   Z^1.5.
!   for iextrasum = nextrasum, sum is over all species excluding
!   bare nucleii and free electrons, the weight is rion^3.
! extrasumf(nextrasum) = partial of extrasum/partial ln f
! extrasumt(nextrasum) = partial of extrasum/partial ln t
! ifdv(nions+2) = 0 means it is safe to exclude this ion from the
!   calculations
! output:
! dvzero(nelements) is the quantity added to dv in ionize.  it
!   is a useful zero point shift to avoid significance loss.
! dv(nions+2) = change in equilibrium constant in ion order with last
!   two H2 and H2+.
! dvf(nions+2) and dvt(nions+2) the ln f and ln t derivatives.

!> This mdh_pi subroutine calculates (for the MDH
!> form of pressure ionization) increments to the
!> equilibrium constant array, dv, and their partial
!> derivatives wrt fl, ln t, and aux_old.
!>
!> \param[in] ifdvzero PARAMETERS NEED DOCUMENTATION
!>
subroutine mdh_pi(ifdvzero, ifnr, inv_aux,&
     iextraoff, inv_ion,&
     partial_elements,&
     ion_end,&
     ifmodified,&
     r_ion3, nion, r_neutral,&
     extrasum, extrasumf, extrasumt,&
     ifdv, dvzero, dv, dvf, dvt, dv_aux)
  use mod_pi_fit, only: quad_original, quad_modified
  use mod_free_eos_constants, only: pi
  use mod_mdh_pi_data, only: mdh_pi_called, quad
  use mod_aux_scale, only: aux_overflow, ln_aux_overflow
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

  ! Arguments
  logical, intent(in) :: ifdvzero
  integer, intent(in) ::&
       ifnr, inv_aux(:), iextraoff,&
       inv_ion(:),&
       partial_elements(:),&
       ion_end(:),&
       ifmodified, nion(:), ifdv(:)

  real(fp_kind), intent(in) :: r_ion3(:), r_neutral(:),&
       extrasum(:), extrasumf(:), extrasumt(:)
  real(fp_kind), intent(inout) ::&
       dvzero(:),&
       dv(:), dvf(:), dvt(:), dv_aux(:,:)

  ! Internal variables
  logical ifnr0, ifnr13, if_neutral_save, if_neutral_save_h, if_neutral_save_h2, if_neutral_save_h2plus
  integer ion,&
       nionsp2, nions, naux, n_partial_elements, mion_end, nelements, nextrasum,&
       ielement,&
       index, max_indexe, ion_start,&
       iextrasum, index_inv_ion, index_aux
  integer, parameter :: maxnextrasum = 9

  real(fp_kind)&
       muneutral, muneutralf, muneutralt, muneutral_aux(maxnextrasum-1),&
       muneutral_h, muneutral_hf, muneutral_ht, muneutral_h_aux(maxnextrasum-1),&
       muneutral_h2, muneutral_h2f, muneutral_h2t, muneutral_h2_aux(maxnextrasum-1),&
       muneutral_h2plus, muneutral_h2plusf, muneutral_h2plust, muneutral_h2plus_aux(maxnextrasum-1),&
       scale_quad, extrasum_s(maxnextrasum-2), muneutralq

  nionsp2 = size(inv_ion)
  nions = nionsp2 - 2
  naux = size(inv_aux)
  n_partial_elements = size(partial_elements) - 2
  mion_end = size(ion_end)
  nelements = size(r_neutral) - 2
  nextrasum = size(extrasum)
  max_indexe = n_partial_elements
  if(ifdv(nions+1).eq.1) max_indexe = max_indexe + 1
  if(ifdv(nions+2).eq.1) max_indexe = max_indexe + 1

  ! sanity checks:
  if(nextrasum.ne.maxnextrasum) error stop 'mdh_pi: nextrasum.ne.maxnextrasum'
  if(&
       nionsp2.ne.size(nion).or.&
       nionsp2.ne.size(ifdv).or.&
       nionsp2.ne.size(r_ion3).or.&
       nionsp2.ne.size(dv).or.&
       nionsp2.ne.size(dvf).or.&
       nionsp2.ne.size(dvt).or.&
       nionsp2.ne.size(dv_aux,1))&
       error stop 'mdh_pi: inconsistent nionsp2 sizes for inv_ion, nion, ifdv, r_ion3, dv, dvf, dvt, or dv_aux'
  if(naux.ne.size(dv_aux,2)) error stop 'mdh_pi: inconsistent naux sizes for inv_aux and dv_aux'
  if(nelements.ne.size(dvzero)) error stop 'mdh_pi: inconsistent nelements sizes for r_neutral and dvzero'
  if(nextrasum.ne.size(extrasumf).or.nextrasum.ne.size(extrasumt))&
       error stop 'mdh_pi: inconsistent nextrasum sizes for extrasum, extrasumf, or extrasumt'
  if(partial_elements(max_indexe).ne.mion_end) then
     write(stderr,*) 'nelements, max_indexe, partial_elements(max_indexe), mion_end = '
     write(stderr,*) nelements, max_indexe, partial_elements(max_indexe), mion_end
     error stop 'mdh_pi: above quantities inconsistent ==> logic error'
  endif

  if(.not.mdh_pi_called) mdh_pi_called = .true.
  ifnr13 = ifnr.eq.1.or.ifnr.eq.3
  ifnr0 = ifnr.eq.0
  if(ifmodified.le.0) then
     ! original mdh value.
     quad = quad_original
  else
     quad = quad_modified
  endif

  if_neutral_save = .false.
  ! Due to if_neutral_save logic below, muneutral, etc., cannot be used uninitialized.
  ! However, gfortran does not recognize that logic so to avoid spurious [-Wmaybe-uninitialized]
  ! warnings, must redundantly initialize here.
  muneutral = 0._fp_kind
  muneutralf = 0._fp_kind
  muneutralt = 0._fp_kind

  if_neutral_save_h = .false.
  ! Due to if_neutral_save_h logic below, muneutral_h, etc., cannot be used uninitialized.
  ! However, gfortran does not recognize that logic so to avoid spurious [-Wmaybe-uninitialized]
  ! warnings, must redundantly initialize here.
  muneutral_h = 0._fp_kind
  muneutral_hf = 0._fp_kind
  muneutral_ht = 0._fp_kind
  muneutral_aux = 0._fp_kind

  if_neutral_save_h2 = .false.
  ! Due to if_neutral_save_h2 logic below, muneutral_h2, etc., cannot be used uninitialized.
  ! However, gfortran does not recognize that logic so to avoid spurious [-Wmaybe-uninitialized]
  ! warnings, must redundantly initialize here.
  muneutral_h2 = 0._fp_kind
  muneutral_h2f = 0._fp_kind
  muneutral_h2t = 0._fp_kind

  if_neutral_save_h2plus = .false.
  ! Due to if_neutral_save_h2plus logic below, muneutral_h2plus, etc., cannot be used uninitialized.
  ! However, gfortran does not recognize that logic so to avoid spurious [-Wmaybe-uninitialized]
  ! warnings, must redundantly initialize here.
  muneutral_h2plus = 0._fp_kind
  muneutral_h2plusf = 0._fp_kind
  muneutral_h2plust = 0._fp_kind

  ! n.b. ifnr=2 component is zero since fixed extrasum.  this routine
  ! therefore should do nothing to dvf and dvt for ifnr = 2.
  do index = 1, max_indexe
     ielement = partial_elements(index)
     if(ielement.gt.1) then
        ion_start = ion_end(ielement-1) + 1
     else
        ion_start = 1
     endif
     do ion = ion_start, ion_end(ielement)
        if(ion.gt.nions.or.nion(min(nions,ion)).eq.1) then
           ! ielement points to index of neutral monatomic or nelement+1
           ! for H2 or nelement+2 for H2+.
           ! Calculate chemical potential of neutral /kt.
           ! neutral-neutral (where "neutral" includes any species with
           ! non-zero radius (all neutrals plus H2+ according to MHD model).
           if_neutral_save = .true.
           muneutral = (4._fp_kind*pi/3._fp_kind)*2._fp_kind*(&
                extrasum(4) +&
                r_neutral(ielement)*(3._fp_kind*extrasum(3) +&
                r_neutral(ielement)*(3._fp_kind*extrasum(2) +&
                r_neutral(ielement)*(extrasum(1)))))
           ! neutral-ion
           muneutral = muneutral + (4._fp_kind*pi/3._fp_kind)*(r_ion3(ion)*extrasum(nextrasum-1))
           if(ifnr0) then
              ! calculate fl, tl derivatives assuming extrasum
              ! is a function of fl and tl with specified
              ! derivatives.
              muneutralf = (4._fp_kind*pi/3._fp_kind)*2._fp_kind*(&
                   extrasumf(4) +&
                   r_neutral(ielement)*(3._fp_kind*extrasumf(3) +&
                   r_neutral(ielement)*(3._fp_kind*extrasumf(2) +&
                   r_neutral(ielement)*(extrasumf(1)))))
              muneutralt = (4._fp_kind*pi/3._fp_kind)*2._fp_kind*(&
                   extrasumt(4) +&
                   r_neutral(ielement)*(3._fp_kind*extrasumt(3) +&
                   r_neutral(ielement)*(3._fp_kind*extrasumt(2) +&
                   r_neutral(ielement)*(extrasumt(1)))))
              muneutralf = muneutralf +&
                   (4._fp_kind*pi/3._fp_kind)*(r_ion3(ion)*extrasumf(nextrasum-1))
              muneutralt = muneutralt +&
                   (4._fp_kind*pi/3._fp_kind)*(r_ion3(ion)*extrasumt(nextrasum-1))
           elseif(ifnr13) then
              muneutral_aux(7) = 0._fp_kind
              muneutral_aux(6) = 0._fp_kind
              muneutral_aux(5) = 0._fp_kind
              muneutral_aux(4) = (4._fp_kind*pi/3._fp_kind)*2._fp_kind
              muneutral_aux(3) = muneutral_aux(4)*3._fp_kind*r_neutral(ielement)
              muneutral_aux(2) = muneutral_aux(3)*r_neutral(ielement)
              muneutral_aux(1) = muneutral_aux(2)*r_neutral(ielement)/3._fp_kind
              muneutral_aux(nextrasum-1) = (4._fp_kind*pi/3._fp_kind)*r_ion3(ion)
           endif
           ! ad hoc neutral-neutral-neutral
           ! Scale factor to help avoid underflow and overflow issues for extreme conditions
           ! or the initial stages of bfgs minimization which can generate extremely wild values
           ! of extrasum.
           scale_quad = maxval(extrasum(1:maxnextrasum-2))
           if(quad.gt.0._fp_kind.and.scale_quad.gt.0._fp_kind) then
              ! Overflow protection required because of wild extrasum values used in bfgs minimization.
              extrasum_s = extrasum(1:maxnextrasum-2)/scale_quad
              muneutralq = 2._fp_kind*log(scale_quad) + log(quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                   3._fp_kind*extrasum_s(4)*extrasum_s(4) +&
                   6._fp_kind*extrasum_s(3)*extrasum_s(5) +&
                   6._fp_kind*extrasum_s(2)*extrasum_s(6) +&
                   2._fp_kind*extrasum_s(1)*extrasum_s(7) +&
                   r_neutral(ielement)*(&
                   30._fp_kind*extrasum_s(3)*extrasum_s(4) +&
                   18._fp_kind*extrasum_s(2)*extrasum_s(5) +&
                   6._fp_kind*extrasum_s(1)*extrasum_s(6) +&
                   r_neutral(ielement)*(&
                   27._fp_kind*extrasum_s(3)*extrasum_s(3) +&
                   30._fp_kind*extrasum_s(2)*extrasum_s(4) +&
                   6._fp_kind*extrasum_s(1)*extrasum_s(5) +&
                   r_neutral(ielement)*(&
                   30._fp_kind*extrasum_s(2)*extrasum_s(3) +&
                   6._fp_kind*extrasum_s(1)*extrasum_s(4) +&
                   r_neutral(ielement)*(&
                   9._fp_kind*extrasum_s(2)*extrasum_s(2) +&
                   6._fp_kind*extrasum_s(1)*extrasum_s(3) +&
                   r_neutral(ielement)*(&
                   6._fp_kind*extrasum_s(1)*extrasum_s(2) +&
                   r_neutral(ielement)*(&
                   1._fp_kind*extrasum_s(1)*extrasum_s(1)))))))))
              if(muneutralq.lt.ln_aux_overflow) then
                 muneutral = muneutral + exp(muneutralq)
                 if(ifnr0) then
                    ! calculate fl, tl derivatives assuming extrasum
                    ! is a function of fl and tl with specified
                    ! derivatives.
                    muneutralf = muneutralf +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         3._fp_kind*(extrasumf(4)*extrasum(4) +&
                         extrasum(4)*extrasumf(4)) +&
                         6._fp_kind*(extrasumf(3)*extrasum(5) +&
                         extrasum(3)*extrasumf(5)) +&
                         6._fp_kind*(extrasumf(2)*extrasum(6) +&
                         extrasum(2)*extrasumf(6)) +&
                         2._fp_kind*(extrasumf(1)*extrasum(7) +&
                         extrasum(1)*extrasumf(7)) +&
                         r_neutral(ielement)*(&
                         30._fp_kind*(extrasumf(3)*extrasum(4) +&
                         extrasum(3)*extrasumf(4)) +&
                         18._fp_kind*(extrasumf(2)*extrasum(5) +&
                         extrasum(2)*extrasumf(5)) +&
                         6._fp_kind*(extrasumf(1)*extrasum(6) +&
                         extrasum(1)*extrasumf(6)) +&
                         r_neutral(ielement)*(&
                         27._fp_kind*(extrasumf(3)*extrasum(3) +&
                         extrasum(3)*extrasumf(3)) +&
                         30._fp_kind*(extrasumf(2)*extrasum(4) +&
                         extrasum(2)*extrasumf(4)) +&
                         6._fp_kind*(extrasumf(1)*extrasum(5) +&
                         extrasum(1)*extrasumf(5)) +&
                         r_neutral(ielement)*(&
                         30._fp_kind*(extrasumf(2)*extrasum(3) +&
                         extrasum(2)*extrasumf(3)) +&
                         6._fp_kind*(extrasumf(1)*extrasum(4) +&
                         extrasum(1)*extrasumf(4)) +&
                         r_neutral(ielement)*(&
                         9._fp_kind*(extrasumf(2)*extrasum(2) +&
                         extrasum(2)*extrasumf(2)) +&
                         6._fp_kind*(extrasumf(1)*extrasum(3) +&
                         extrasum(1)*extrasumf(3)) +&
                         r_neutral(ielement)*(&
                         6._fp_kind*(extrasumf(1)*extrasum(2) +&
                         extrasum(1)*extrasumf(2)) +&
                         r_neutral(ielement)*(&
                         1._fp_kind*(extrasumf(1)*extrasum(1) +&
                         extrasum(1)*extrasumf(1)))))))))
                    muneutralt = muneutralt +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         3._fp_kind*(extrasumt(4)*extrasum(4) +&
                         extrasum(4)*extrasumt(4)) +&
                         6._fp_kind*(extrasumt(3)*extrasum(5) +&
                         extrasum(3)*extrasumt(5)) +&
                         6._fp_kind*(extrasumt(2)*extrasum(6) +&
                         extrasum(2)*extrasumt(6)) +&
                         2._fp_kind*(extrasumt(1)*extrasum(7) +&
                         extrasum(1)*extrasumt(7)) +&
                         r_neutral(ielement)*(&
                         30._fp_kind*(extrasumt(3)*extrasum(4) +&
                         extrasum(3)*extrasumt(4)) +&
                         18._fp_kind*(extrasumt(2)*extrasum(5) +&
                         extrasum(2)*extrasumt(5)) +&
                         6._fp_kind*(extrasumt(1)*extrasum(6) +&
                         extrasum(1)*extrasumt(6)) +&
                         r_neutral(ielement)*(&
                         27._fp_kind*(extrasumt(3)*extrasum(3) +&
                         extrasum(3)*extrasumt(3)) +&
                         30._fp_kind*(extrasumt(2)*extrasum(4) +&
                         extrasum(2)*extrasumt(4)) +&
                         6._fp_kind*(extrasumt(1)*extrasum(5) +&
                         extrasum(1)*extrasumt(5)) +&
                         r_neutral(ielement)*(&
                         30._fp_kind*(extrasumt(2)*extrasum(3) +&
                         extrasum(2)*extrasumt(3)) +&
                         6._fp_kind*(extrasumt(1)*extrasum(4) +&
                         extrasum(1)*extrasumt(4)) +&
                         r_neutral(ielement)*(&
                         9._fp_kind*(extrasumt(2)*extrasum(2) +&
                         extrasum(2)*extrasumt(2)) +&
                         6._fp_kind*(extrasumt(1)*extrasum(3) +&
                         extrasum(1)*extrasumt(3)) +&
                         r_neutral(ielement)*(&
                         6._fp_kind*(extrasumt(1)*extrasum(2) +&
                         extrasum(1)*extrasumt(2)) +&
                         r_neutral(ielement)*(&
                         1._fp_kind*(extrasumt(1)*extrasum(1) +&
                         extrasum(1)*extrasumt(1)))))))))
                 elseif(ifnr13) then
                    muneutral_aux(1) = muneutral_aux(1) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         2._fp_kind*extrasum(7) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(6) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(5) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(4) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(3) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(2) + r_neutral(ielement)*(&
                         2._fp_kind*extrasum(1) )))))))
                    muneutral_aux(2) = muneutral_aux(2) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         6._fp_kind*extrasum(6) + r_neutral(ielement)*(&
                         18._fp_kind*extrasum(5) + r_neutral(ielement)*(&
                         30._fp_kind*extrasum(4) + r_neutral(ielement)*(&
                         30._fp_kind*extrasum(3) + r_neutral(ielement)*(&
                         18._fp_kind*extrasum(2) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(1) ))))))
                    muneutral_aux(3) = muneutral_aux(3) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         6._fp_kind*extrasum(5) + r_neutral(ielement)*(&
                         30._fp_kind*extrasum(4) + r_neutral(ielement)*(&
                         54._fp_kind*extrasum(3) + r_neutral(ielement)*(&
                         30._fp_kind*extrasum(2) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(1) )))))
                    muneutral_aux(4) = muneutral_aux(4) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         6._fp_kind*extrasum(4) + r_neutral(ielement)*(&
                         30._fp_kind*extrasum(3) + r_neutral(ielement)*(&
                         30._fp_kind*extrasum(2) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(1) ))))
                    muneutral_aux(5) = muneutral_aux(5) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         6._fp_kind*extrasum(3) + r_neutral(ielement)*(&
                         18._fp_kind*extrasum(2) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(1) )))
                    muneutral_aux(6) = muneutral_aux(6) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         6._fp_kind*extrasum(2) + r_neutral(ielement)*(&
                         6._fp_kind*extrasum(1) ))
                    muneutral_aux(7) = muneutral_aux(7) +&
                         quad*(16._fp_kind*pi*pi/9._fp_kind)*(&
                         2._fp_kind*extrasum(1) )
                 endif !elseif(ifnr13)....
              else ! if(muneutralq...
                 ! On assumption that aux_overflow dominates everything else contributing
                 ! to muneutral, pin muneutral value to aux_overflow and its derivatives to zero for this case.
                 muneutral = aux_overflow
                 if(ifnr0) then
                    muneutralf = 0._fp_kind
                    muneutralt = 0._fp_kind
                 elseif(ifnr13) then
                    muneutral_aux = 0._fp_kind
                 endif
              endif ! if(muneutralq...
           endif !if(quad.gt.0._fp_kindd....
           if(ion.le.nions) then
              if(ifdvzero.or.ielement.eq.1) then
                 ! significance loss not affected for
                 ! hydrogen with just one neutral and
                 ! one ionized state so don't bother
                 ! with zero point shift.
              else
                 ! dvzero is a zero point shift that
                 ! is added to dv in ionize so must
                 ! be subtracted from dv here.  for this
                 ! case use a zero point shift that
                 ! renders the bare ion pressure ionization
                 ! dv contribution zero.
                 ! Sanity check
                 if(.not.if_neutral_save) error stop 'mdh_pi: (1) muneutral, etc., not initialized ==> logic error'

                 dvzero(ielement) = dvzero(ielement) + muneutral -&
                      (4._fp_kind*pi/3._fp_kind)*(real((nion(ion_end(ielement))),fp_kind)*&
                      sqrt(real((nion(ion_end(ielement))),fp_kind))*extrasum(nextrasum))
              endif
           endif
        endif !if(ion.gt.nions.or.nion.....
        if(ion.le.nions) then
           ! Sanity check
           if(.not.if_neutral_save) error stop 'mdh_pi: (2) muneutral, etc., not initialized ==> logic error'

           if(ifdvzero.or.ielement.eq.1) then
              dv(ion) = dv(ion) + muneutral
              ! add portion of negative of chemical potential of ion species/kt
              dv(ion) = dv(ion) -&
                   (4._fp_kind*pi/3._fp_kind)*(real((nion(ion)),fp_kind)*sqrt(real((nion(ion)),fp_kind))*extrasum(nextrasum))
           else
              ! subtract out dvzero without causing significance loss.
              ! n.b. i^1.5-j^1.5 = (i^3-j^3)/(i^1.5+j^1.5)
              dv(ion) = dv(ion) - (4._fp_kind*pi/3._fp_kind)*extrasum(nextrasum)*&
                   real(nion(ion)*nion(ion)*nion(ion) -&
                   nion(ion_end(ielement))*nion(ion_end(ielement))*nion(ion_end(ielement)),fp_kind)/&
                   (real(nion(ion),fp_kind)*sqrt(real(nion(ion),fp_kind)) +&
                   real(nion(ion_end(ielement)),fp_kind)*sqrt(real(nion(ion_end(ielement)),fp_kind)))
           endif
           if(ifnr0) then
              dvf(ion) = dvf(ion) + muneutralf
              dvt(ion) = dvt(ion) + muneutralt
              dvf(ion) = dvf(ion) -&
                   (4._fp_kind*pi/3._fp_kind)*(real((nion(ion)),fp_kind)*sqrt(real((nion(ion)),fp_kind))*extrasumf(nextrasum))
              dvt(ion) = dvt(ion) -&
                   (4._fp_kind*pi/3._fp_kind)*(real((nion(ion)),fp_kind)*sqrt(real((nion(ion)),fp_kind))*extrasumt(nextrasum))
           elseif(ifnr13) then
              index_inv_ion = inv_ion(ion)
              do iextrasum = 1, nextrasum-1
                 index_aux = inv_aux(iextraoff+iextrasum)
                 dv_aux(index_inv_ion,index_aux) = dv_aux(index_inv_ion,index_aux) + muneutral_aux(iextrasum)
              enddo
              index_aux = inv_aux(iextraoff+nextrasum)
              dv_aux(index_inv_ion,index_aux) = dv_aux(index_inv_ion,index_aux) -&
                   (4._fp_kind*pi/3._fp_kind)*real((nion(ion)),fp_kind)*sqrt(real((nion(ion)),fp_kind))
           endif
           if(ion.lt.nions.and.nion(min(nions,ion+1)).gt.1) then
              ! if not last ion (usually bare nucleus, but when
              ! not bare nucleus this is done consistently) add
              ! in relevant additional component.
              ! add negative of chemical potential of species/kt
              dv(ion) = dv(ion) -&
                   (4._fp_kind*pi/3._fp_kind)*(r_ion3(ion+1)*&
                   extrasum(nextrasum-1))
              if(ifnr0) then
                 dvf(ion) = dvf(ion) -&
                      (4._fp_kind*pi/3._fp_kind)*(r_ion3(ion+1)*&
                      extrasumf(nextrasum-1))
                 dvt(ion) = dvt(ion) -&
                      (4._fp_kind*pi/3._fp_kind)*(r_ion3(ion+1)*&
                      extrasumt(nextrasum-1))
              elseif(ifnr13) then
                 index_aux = inv_aux(iextraoff+nextrasum-1)
                 dv_aux(index_inv_ion,index_aux) =&
                      dv_aux(index_inv_ion,index_aux) -&
                      (4._fp_kind*pi/3._fp_kind)*r_ion3(ion+1)
              endif
           endif
        endif !if(ion.le......
        if(ion.eq.1) then
           ! Sanity check
           if(.not.if_neutral_save) error stop 'mdh_pi: (3) muneutral, etc., not initialized ==> logic error'

           ! save muneutral_h, etc., for later H2, H2+ calculation.
           if_neutral_save_h = .true.
           muneutral_h = muneutral
           if(ifnr0) then
              muneutral_hf = muneutralf
              muneutral_ht = muneutralt
           elseif(ifnr13) then
              muneutral_h_aux(:nextrasum-1) = muneutral_aux(:nextrasum-1)
           endif
        endif
        if(ion.eq.nions+1) then
           ! Sanity check
           if(.not.if_neutral_save) error stop 'mdh_pi: (4) muneutral, etc., not initialized ==> logic error'

           ! save muneutral_h2, etc.,  for later H2, H2+ calculation.
           if_neutral_save_h2 = .true.
           muneutral_h2 = muneutral
           if(ifnr0) then
              muneutral_h2f = muneutralf
              muneutral_h2t = muneutralt
           elseif(ifnr13) then
              muneutral_h2_aux(:nextrasum-1) = muneutral_aux(:nextrasum-1)
           endif
        endif
        if(ion.eq.nions+2) then
           ! Sanity check
           if(.not.if_neutral_save) error stop 'mdh_pi: (5) muneutral, etc., not initialized ==> logic error'

           ! save muneutral_h2plus, etc., for later H2+ calculation.
           if_neutral_save_h2plus = .true.
           muneutral_h2plus = muneutral
           if(ifnr0) then
              muneutral_h2plusf = muneutralf
              muneutral_h2plust = muneutralt
           elseif(ifnr13) then
              muneutral_h2plus_aux(:nextrasum-1) = muneutral_aux(:nextrasum-1)
           endif
        endif
     enddo !do ion = ion_start,....
  enddo !do index = 1,.......
  ! special handling of H2, and H2+.  At this point:
  ! muneutral_h holds partial F partial n(H)/kt,
  ! muneutral_h2 holds partial F partial n(H2)/kt
  ! muneutral_h2plus holds partial F partial n(H2+)/kt
  if(ifdv(nions+1).eq.1) then
     ! Sanity checks
     if(.not.if_neutral_save_h) error stop 'mdh_pi: muneutralh, etc., not initialized ==> logic error'
     if(.not.if_neutral_save_h2) error stop 'mdh_pi: muneutralh2, etc., not initialized ==> logic error'

     ! by definition of H2 equilibrium constant relative to neutral monatomic
     dv(nions+1) = dv(nions+1) + 2._fp_kind*muneutral_h - muneutral_h2
     if(ifnr0) then
        dvf(nions+1) = dvf(nions+1) + 2._fp_kind*muneutral_hf - muneutral_h2f
        dvt(nions+1) = dvt(nions+1) + 2._fp_kind*muneutral_ht - muneutral_h2t
     elseif(ifnr13) then
        index_inv_ion = inv_ion(nions+1)
        do iextrasum = 1, nextrasum-1
           index_aux = inv_aux(iextraoff+iextrasum)
           dv_aux(index_inv_ion,index_aux) = dv_aux(index_inv_ion,index_aux) +&
                2._fp_kind*muneutral_h_aux(iextrasum) - muneutral_h2_aux(iextrasum)
        enddo
     endif
  endif
  if(ifdv(nions+2).eq.1) then
     ! Sanity checks
     if(.not.if_neutral_save_h2) error stop 'mdh_pi: muneutralh2, etc., not initialized ==> logic error'
     if(.not.if_neutral_save_h2plus) error stop 'mdh_pi: muneutralh2plus, etc., not initialized ==> logic error'

     ! by definition of H2+ equilibrium constant relative to H2 and e-.
     ! n.b. the extrasum(nextrasum-1) component is already taken care
     ! of in muneutral_h2plus
     dv(nions+2) = dv(nions+2) + muneutral_h2 - muneutral_h2plus -&
          (4._fp_kind*pi/3._fp_kind)*extrasum(nextrasum)
     if(ifnr0) then
        dvf(nions+2) = dvf(nions+2) + muneutral_h2f - muneutral_h2plusf -&
             (4._fp_kind*pi/3._fp_kind)*extrasumf(nextrasum)
        dvt(nions+2) = dvt(nions+2) + muneutral_h2t - muneutral_h2plust -&
             (4._fp_kind*pi/3._fp_kind)*extrasumt(nextrasum)
     elseif(ifnr13) then
        index_inv_ion = inv_ion(nions+2)
        do iextrasum = 1, nextrasum-1
           index_aux = inv_aux(iextraoff+iextrasum)
           dv_aux(index_inv_ion,index_aux) = dv_aux(index_inv_ion,index_aux) +&
                muneutral_h2_aux(iextrasum) - muneutral_h2plus_aux(iextrasum)
        enddo
        index_aux = inv_aux(iextraoff+nextrasum)
        dv_aux(index_inv_ion,index_aux) = dv_aux(index_inv_ion,index_aux) -&
             (4._fp_kind*pi/3._fp_kind)
     endif
  endif
end subroutine mdh_pi

! Calculate MDH-like free energy/unit volume model given by:
! f = kt*4 pi/3 *[
!   2(3 alpha1*alpha2 + alpha0*alpha3) +
!   gamma*beta] + f'
! f' is an ad hoc extra term designed (MDH paper II, Appendix B) to
!   force pressure ionization for low temperatures.  It has the
!   right functional form for the next higher order interaction.
!   if one expands out the various powers of r in the MDH expression, then
! f' = quad*kt*(4 pi/3)^2 *[
!   3 alpha0*alpha3*alpha3 +
!   30 alpha1*alpha2*alpha3 +
!   9 alpha2*alpha2*alpha2 +
!   6 alpha0*alpha2*alpha4 +
!   9 alpha1*alpha1*alpha4 +
!   6 alpha0*alpha1*alpha5 +
!   1 alpha0*alpha0*alpha6]
! Calculate f and fquad divided by t.
! n.b. extrasum is in number per unit volume form

!> This mdh_pi_pressure_free subroutine calculates (for the MDH
!> form of pressure ionization) the pressure and
!> the Helmholtz free energy per unit volume and the partial
!> derivatives of those quantities wrt fl, ln t, and aux_old.
!>
!> \param[in] t PARAMETERS NEED DOCUMENTATION
!>
subroutine mdh_pi_pressure_free(t, extrasum,&
       ppi, ppif, ppit, ppi2_aux,&
       free_pi, free_pif, free_pi2_aux)
  use mod_free_eos_constants, only: boltzmann, pi
  use mod_mdh_pi_data, only: mdh_pi_called, quad
  use mod_aux_scale, only: aux_overflow, ln_aux_overflow

  ! Arguments
  real(fp_kind), intent(in) :: t, extrasum(:)
  real(fp_kind), intent(out) ::&
       ppi, ppif, ppit, ppi2_aux(:),&
       free_pi, free_pif, free_pi2_aux(:)

  ! Internal variables
  integer, parameter :: maxnextrasum = 9

  integer nextrasum

  real(fp_kind) fdt, fdt_aux(maxnextrasum),&
       fquaddt, fquaddt_aux(maxnextrasum), extrasum_s(maxnextrasum), scale_quad
  ! must zero so that later logic involving these arrays works properly.
  data fdt_aux, fquaddt_aux/18*0._fp_kind/

  ! Most/all Fortran compilers use the save attribute for variables initialized by
  ! a data statement, but just in case ...
  save fdt_aux, fquaddt_aux

  nextrasum = size(extrasum)

  ! sanity checks
  if(.not.mdh_pi_called) error stop 'mdh_pi_pressure_free: mdh_pi must be called first'
  if(nextrasum.ne.maxnextrasum) error stop 'mdh_pi_pressure_free: nextrasum.ne.maxnextrasum'
  if(nextrasum.ne.size(ppi2_aux).or.nextrasum.ne.size(free_pi2_aux))&
       error stop 'mdh_pi_pressure_free: inconsistent nextrasum sizes for extrasum, ppi2_aux, or free_pi2_aux'

  ! Scale factor to help avoid underflow and overflow issues for extreme conditions
  ! or the initial stages of bfgs minimization which can generate extremely wild values
  ! of extrasum.
  scale_quad = maxval(extrasum)
  if(scale_quad.gt.0._fp_kind) then
     ! All calculations below use extrasum_s so after those
     ! calculations must multiply calculated fdt by scale_quad*2, the
     ! calculated fdt derivatives by scale_quad, the calculated
     ! fquaddt by scale_quad**3, and the calculated fquaddt
     ! derivatives by scale_quad**2.
     extrasum_s = extrasum/scale_quad
     fdt = boltzmann*(4._fp_kind*pi/3._fp_kind)*(&
          2._fp_kind*(3._fp_kind*extrasum_s(2)*extrasum_s(3) +&
          extrasum_s(1)*extrasum_s(4)) +&
          extrasum_s(nextrasum-1)*extrasum_s(nextrasum))
     fquaddt = quad*boltzmann*(16._fp_kind*pi*pi/9._fp_kind)*(&
          3._fp_kind*extrasum_s(1)*extrasum_s(4)*extrasum_s(4) +&
          30._fp_kind*extrasum_s(2)*extrasum_s(3)*extrasum_s(4) +&
          9._fp_kind*extrasum_s(3)*extrasum_s(3)*extrasum_s(3) +&
          6._fp_kind*extrasum_s(1)*extrasum_s(3)*extrasum_s(5) +&
          9._fp_kind*extrasum_s(2)*extrasum_s(2)*extrasum_s(5) +&
          6._fp_kind*extrasum_s(1)*extrasum_s(2)*extrasum_s(6) +&
          1._fp_kind*extrasum_s(1)*extrasum_s(1)*extrasum_s(7))
     if(fdt.gt.0._fp_kind.or.fquaddt.gt.0._fp_kind) then
        fdt_aux(1) = boltzmann*(4._fp_kind*pi/3._fp_kind)*(2._fp_kind*extrasum_s(4))
        fdt_aux(2) = boltzmann*(4._fp_kind*pi/3._fp_kind)*(6._fp_kind*extrasum_s(3))
        fdt_aux(3) = boltzmann*(4._fp_kind*pi/3._fp_kind)*(6._fp_kind*extrasum_s(2))
        fdt_aux(4) = boltzmann*(4._fp_kind*pi/3._fp_kind)*(2._fp_kind*extrasum_s(1))
        fdt_aux(nextrasum-1) = boltzmann*(4._fp_kind*pi/3._fp_kind)*extrasum_s(nextrasum)
        fdt_aux(nextrasum) = boltzmann*(4._fp_kind*pi/3._fp_kind)*extrasum_s(nextrasum-1)
        fquaddt_aux(1) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             3._fp_kind*extrasum_s(4)*extrasum_s(4) +&
             6._fp_kind*extrasum_s(3)*extrasum_s(5) +&
             6._fp_kind*extrasum_s(2)*extrasum_s(6) +&
             2._fp_kind*extrasum_s(1)*extrasum_s(7)&
             ))
        fquaddt_aux(2) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             30._fp_kind*extrasum_s(3)*extrasum_s(4) +&
             18._fp_kind*extrasum_s(2)*extrasum_s(5) +&
             6._fp_kind*extrasum_s(1)*extrasum_s(6)&
             ))
        fquaddt_aux(3) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             30._fp_kind*extrasum_s(2)*extrasum_s(4) +&
             27._fp_kind*extrasum_s(3)*extrasum_s(3) +&
             6._fp_kind*extrasum_s(1)*extrasum_s(5)&
             ))
        fquaddt_aux(4) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             6._fp_kind*extrasum_s(1)*extrasum_s(4) +&
             30._fp_kind*extrasum_s(2)*extrasum_s(3) +&
             6._fp_kind*extrasum_s(1)*extrasum_s(5)&
             ))
        fquaddt_aux(5) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             6._fp_kind*extrasum_s(1)*extrasum_s(3) +&
             9._fp_kind*extrasum_s(2)*extrasum_s(2)&
             ))
        fquaddt_aux(6) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             6._fp_kind*extrasum_s(1)*extrasum_s(2)&
             ))
        fquaddt_aux(7) = boltzmann*(4._fp_kind*pi/3._fp_kind)*&
             (&
             quad*(4._fp_kind*pi/3._fp_kind)*(&
             1._fp_kind*extrasum_s(1)*extrasum_s(1)&
             ))
        ! F(T,V,N) = t*V*fdt + t*V*fquaddt.
        ! additional V factors divide to convert extrasum products to N form.
        ! Therefore, first term is proportional to V^{-1} and
        ! second term is proportional to V^{-2}
        ! P = -partial F(T,V,N)/partial V
        ! ppi = scale_quad*scale_quad*t*(fdt + 2.d0*scale_quad*fquaddt)
        ! Overflow protection
        ppi = 2._fp_kind*log(scale_quad) + log(t*(fdt + 2._fp_kind*scale_quad*fquaddt))
        if(ppi.lt.ln_aux_overflow) then
           ppi = exp(ppi)
           ppif = 0._fp_kind
           ppit = ppi
           ! n.b. some parts of this array are initialized to zero by
           ! data statement earlier in routine.
           ppi2_aux(:nextrasum) = scale_quad*t*(fdt_aux(:nextrasum) + 2._fp_kind*scale_quad*fquaddt_aux(:nextrasum))
        else
           ppi = aux_overflow
           ppif = 0._fp_kind
           ppit = 0._fp_kind
           ppi2_aux(:nextrasum) = 0._fp_kind
        endif

        !free_pi = scale_quad*scale_quad*t*(fdt + scale_quad*fquaddt)
        ! Overflow protection
        free_pi = 2._fp_kind*log(scale_quad) + log(t*(fdt + scale_quad*fquaddt))
        if(free_pi.lt.ln_aux_overflow) then
           free_pi = exp(free_pi)
           free_pif = 0._fp_kind
           ! n.b. some parts of this array are initialized to zero by
           ! data statement earlier in routine.
           free_pi2_aux(:nextrasum) = scale_quad*t*(fdt_aux(:nextrasum) + scale_quad*fquaddt_aux(:nextrasum))
        else
           free_pi = aux_overflow
           free_pif = 0._fp_kind
           free_pi2_aux(:nextrasum) = 0._fp_kind
        endif
     else !if(fdt.gt.0._fp_kind.or.fquaddt.gt.0._fp_kind) then
        ppi = 0._fp_kind
        ppif = 0._fp_kind
        ppit = 0._fp_kind
        ppi2_aux(:nextrasum) = 0._fp_kind
        free_pi = 0._fp_kind
        free_pif = 0._fp_kind
        free_pi2_aux(:nextrasum) = 0._fp_kind
     endif
  else
     ppi = 0._fp_kind
     ppif = 0._fp_kind
     ppit = 0._fp_kind
     ppi2_aux(:nextrasum) = 0._fp_kind
     free_pi = 0._fp_kind
     free_pif = 0._fp_kind
     free_pi2_aux(:nextrasum) = 0._fp_kind
  endif
end subroutine mdh_pi_pressure_free

! compute remaining quantities having found ionization balance
! MDH-like free energy/unit volume model given by:
! f = kt*4 pi/3 *[
!   2(3 alpha1*alpha2 + alpha0*alpha3) +
!   gamma*beta] + f'
! f' is an ad hoc extra term designed (MDH paper II, Appendix B) to
!   force pressure ionization for low temperatures.  It has the
!   right functional form for the next higher order interaction.
!   if one expands out the various powers of r in the MDH expression, then
! f' = quad*kt*(4 pi/3)^2 *[
!   3 alpha0*alpha3*alpha3 +
!   30 alpha1*alpha2*alpha3 +
!   9 alpha2*alpha2*alpha2 +
!   6 alpha0*alpha2*alpha4 +
!   9 alpha1*alpha1*alpha4 +
!   6 alpha0*alpha1*alpha5 +
!   1 alpha0*alpha0*alpha6]
! spi = s per unit mass = - partial fV/partial T/m = -f/(rho*T)
! correct for quad part later.
! n.b. extrasum is in nu = n/(rho*avogadro) form for mdh_pi_end call

!> This mdh_pi_end subroutine calculates (for the MDH form of pressure
!> ionization) the pressure, entropy per unit mass, and energy per
!> unit mass as well as the partial derivatives of those first two
!> quantities wrt fl and ln t.
!>
!> \param[in] t PARAMETERS NEED DOCUMENTATION
!>
subroutine mdh_pi_end(t, rho, rf, rt,&
       nion,&
       extrasum, extrasumf, extrasumt,&
       ppi, ppif, ppit, spi, spif, spit, upi)
  use mod_free_eos_constants, only: avogadro, cr, pi
  use mod_mdh_pi_data, only: mdh_pi_called, quad, pi_trace_state, pi_trace_extra

  !Arguments
  integer, intent(in) :: nion(:)
  real(fp_kind), intent(in) :: t, rho, rf, rt,&
       extrasum(:), extrasumf(:), extrasumt(:)
  real(fp_kind), intent(out) :: ppi, ppif, ppit, spi, spif, spit, upi

  ! Internal variables
  integer nextrasum
  real(fp_kind) squad, squadf, squadt

  nextrasum = size(extrasum)

  ! Sanity check
  if(nextrasum.ne.size(extrasumf).or.nextrasum.ne.size(extrasumt))&
       error stop 'mdh_pi_end: inconsistent nextrasum sizes for extrasum, extrasumf, or extrasumf'

  if(.not.mdh_pi_called) error stop 'mdh_pi_end: mdh_pi must be called first'

  ! No overflow issues should need to be addressed below since mdh_pi_end
  ! (unlike mdh_pi_pressure_free above) is only called for converged solutions
  ! (i.e., is not called with wild values of extrasum that can occur for bfgs minimization).

  spi = -cr*(4._fp_kind*pi/3._fp_kind)*(&
       2._fp_kind*(3._fp_kind*extrasum(2)*extrasum(3) + extrasum(1)*extrasum(4)) +&
       extrasum(nextrasum-1)*extrasum(nextrasum))*(rho*avogadro)
  spif = -cr*(4._fp_kind*pi/3._fp_kind)*(&
       2._fp_kind*(3._fp_kind*(extrasumf(2)*extrasum(3) + extrasum(2)*extrasumf(3)) +&
       extrasumf(1)*extrasum(4) + extrasum(1)*extrasumf(4)) +&
       extrasumf(nextrasum-1)*extrasum(nextrasum) +&
       extrasum(nextrasum-1)*extrasumf(nextrasum))*(rho*avogadro) +&
       spi*rf
  spit = -cr*(4._fp_kind*pi/3._fp_kind)*(&
       2._fp_kind*(3._fp_kind*(extrasumt(2)*extrasum(3) + extrasum(2)*extrasumt(3)) +&
       extrasumt(1)*extrasum(4) + extrasum(1)*extrasumt(4)) +&
       extrasumt(nextrasum-1)*extrasum(nextrasum) +&
       extrasum(nextrasum-1)*extrasumt(nextrasum))*(rho*avogadro) +&
       spi*rt
  ! p = -partial (fV)/partial V
  ! fV proportional to V^-1 and ppi = f = -spi*rho*t
  ppi = -spi*rho*t
  ppif = -(spif+spi*rf)*rho*t
  ppit = -(spit+spi*(1._fp_kind+rt))*rho*t
  ! upi = internal energy/unit mass = fpi/rho + t*spi
  !   = -spi*t + t*spi = 0
  upi = 0._fp_kind
  if(quad.le.0._fp_kind) return
  ! s per unit mass = - partial fquad V/partial T/m = -fquad/(rho*t)
  ! n.b. extrasum is in nu = n/(rho*avogadro) form for mdh_pi_end call
  squad = -quad*cr*(16._fp_kind*pi*pi/9._fp_kind)*(&
       3._fp_kind*extrasum(1)*extrasum(4)*extrasum(4) +&
       30._fp_kind*extrasum(2)*extrasum(3)*extrasum(4) +&
       9._fp_kind*extrasum(3)*extrasum(3)*extrasum(3) +&
       6._fp_kind*extrasum(1)*extrasum(3)*extrasum(5) +&
       9._fp_kind*extrasum(2)*extrasum(2)*extrasum(5) +&
       6._fp_kind*extrasum(1)*extrasum(2)*extrasum(6) +&
       1._fp_kind*extrasum(1)*extrasum(1)*extrasum(7))*&
       (rho*avogadro)*(rho*avogadro)
  squadf = -quad*cr*(16._fp_kind*pi*pi/9._fp_kind)*(&
       3._fp_kind*extrasumf(1)*extrasum(4)*extrasum(4) +&
       3._fp_kind*extrasum(1)*extrasumf(4)*extrasum(4) +&
       3._fp_kind*extrasum(1)*extrasum(4)*extrasumf(4) +&
       30._fp_kind*extrasumf(2)*extrasum(3)*extrasum(4) +&
       30._fp_kind*extrasum(2)*extrasumf(3)*extrasum(4) +&
       30._fp_kind*extrasum(2)*extrasum(3)*extrasumf(4) +&
       9._fp_kind*extrasumf(3)*extrasum(3)*extrasum(3) +&
       9._fp_kind*extrasum(3)*extrasumf(3)*extrasum(3) +&
       9._fp_kind*extrasum(3)*extrasum(3)*extrasumf(3) +&
       6._fp_kind*extrasumf(1)*extrasum(3)*extrasum(5) +&
       6._fp_kind*extrasum(1)*extrasumf(3)*extrasum(5) +&
       6._fp_kind*extrasum(1)*extrasum(3)*extrasumf(5) +&
       9._fp_kind*extrasumf(2)*extrasum(2)*extrasum(5) +&
       9._fp_kind*extrasum(2)*extrasumf(2)*extrasum(5) +&
       9._fp_kind*extrasum(2)*extrasum(2)*extrasumf(5) +&
       6._fp_kind*extrasumf(1)*extrasum(2)*extrasum(6) +&
       6._fp_kind*extrasum(1)*extrasumf(2)*extrasum(6) +&
       6._fp_kind*extrasum(1)*extrasum(2)*extrasumf(6) +&
       1._fp_kind*extrasumf(1)*extrasum(1)*extrasum(7) +&
       1._fp_kind*extrasum(1)*extrasumf(1)*extrasum(7) +&
       1._fp_kind*extrasum(1)*extrasum(1)*extrasumf(7))*&
       (rho*avogadro)*(rho*avogadro) + 2._fp_kind*squad*rf
  squadt = -quad*cr*(16._fp_kind*pi*pi/9._fp_kind)*(&
       3._fp_kind*extrasumt(1)*extrasum(4)*extrasum(4) +&
       3._fp_kind*extrasum(1)*extrasumt(4)*extrasum(4) +&
       3._fp_kind*extrasum(1)*extrasum(4)*extrasumt(4) +&
       30._fp_kind*extrasumt(2)*extrasum(3)*extrasum(4) +&
       30._fp_kind*extrasum(2)*extrasumt(3)*extrasum(4) +&
       30._fp_kind*extrasum(2)*extrasum(3)*extrasumt(4) +&
       9._fp_kind*extrasumt(3)*extrasum(3)*extrasum(3) +&
       9._fp_kind*extrasum(3)*extrasumt(3)*extrasum(3) +&
       9._fp_kind*extrasum(3)*extrasum(3)*extrasumt(3) +&
       6._fp_kind*extrasumt(1)*extrasum(3)*extrasum(5) +&
       6._fp_kind*extrasum(1)*extrasumt(3)*extrasum(5) +&
       6._fp_kind*extrasum(1)*extrasum(3)*extrasumt(5) +&
       9._fp_kind*extrasumt(2)*extrasum(2)*extrasum(5) +&
       9._fp_kind*extrasum(2)*extrasumt(2)*extrasum(5) +&
       9._fp_kind*extrasum(2)*extrasum(2)*extrasumt(5) +&
       6._fp_kind*extrasumt(1)*extrasum(2)*extrasum(6) +&
       6._fp_kind*extrasum(1)*extrasumt(2)*extrasum(6) +&
       6._fp_kind*extrasum(1)*extrasum(2)*extrasumt(6) +&
       1._fp_kind*extrasumt(1)*extrasum(1)*extrasum(7) +&
       1._fp_kind*extrasum(1)*extrasumt(1)*extrasum(7) +&
       1._fp_kind*extrasum(1)*extrasum(1)*extrasumt(7))*&
       (rho*avogadro)*(rho*avogadro) + 2._fp_kind*squad*rt
  spi = spi + squad
  spif = spif + squadf
  spit = spit + squadt
  ! n.b. fquadV is proportional to V^-2, so
  !   pquad = 2*fquad = -2*squad*rho*t
  ppi = ppi - 2._fp_kind*squad*rho*t
  ppif = ppif - 2._fp_kind*(squadf+squad*rf)*rho*t
  ppit = ppit - 2._fp_kind*(squadt+squad*(1._fp_kind+rt))*rho*t
  ! upi = internal energy/unit mass = fquad/rho + t*squad
  !   = pquad/(2*rho) + t*squad = -squad*t + t*squad = 0
  ! upi = upi + 0.d0
  pi_trace_extra(:,1)=extrasum
  pi_trace_extra(:,2)=extrasumf
  pi_trace_extra(:,3)=extrasumt
  pi_trace_state=[t,rho,rf,rt,ppi,ppif,ppit,spi,spif,spit,upi,quad]
end subroutine mdh_pi_end
