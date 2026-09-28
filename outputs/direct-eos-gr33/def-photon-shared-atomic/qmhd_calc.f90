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

! calculate mhd partition function sums for nlevel terms
! summand is w_pl*w_mhd*weight*exp(a), where a = c_2*R/neff^2/T
! w_pl is the planck_larkin occupation probability (ifpl true)
! w_mhd is the mhd occupation probability.
! w_pl is given by (1 - exp(-a)*(1+a)).
! The *ln* of the mhd occupation probability is given by
! -4pi/3*[x(1) + x(2)*r + x(3)*r^2 +
! x(4)*r^3+ x(5)*rion^3], where
! r = bohr/2*(3neff^2-ell(ell-1))
! rion^3 = 16*(iz^0.5 e^2 neff^2)^3/(K^0.5 hc Rydberg)
! K(neff) is a quantum correction term

! input quantities:
! ifpi_fit = 2, use best fit to Saumon table
! ifpi_fit = 1, use best fit to opal table + extensions
! ifpi_fit = 0, use best fit to original MDH table.
! ifpl (logical) controls whether to use the planck-larkin occupation probability.
! weight_ion = statistical weight of next higher ion =
!   the factor that multiplies Rydberg part of partition function
! ifapprox controls whether to use approximation or exact sum for rydberg states.
! eps_factor = exp(-c2 bion(ion)/t) is used to terminate rydberg state sum.
! nmin_rydberg is the lowest principal quantum number treated as a rydberg level.
! nmax is maximum principal number cutoff (only an issue for this
!   call to qstar_calc when ifapprox = .false.)  Also used here
!   to cutoff exact non-Rydberg He calculation.
! tl = ln T
! iz = core charge = 1 for neutrals, 2 for first ions, etc.
! rfactor = finite mass factor for Rydberg = 1/(1+m_e/M)
! nlevel = number of levels included in sum before use rydberg formulas
!   (other routines) to finish off to principal quantum number of infinity.
! weight(nlevel) is the statistical weight
! neff(nlevel) is the effective principal quantum number
! ell(nlevel) is the azimuthal quantum number
! x(nx=5) is ordered as extrasum(4), extrasum(3), extrasum(2), extrasum(1)
! extrasum(nextrasum-1)

! output quantities:
! qmhd, qmhdt, qmhdx(nx), qmhdt2, qmhdtx(nx), qmhdx2(nx, nx) are
! the resulting sum plus partial derivatives wrt tl and x.

!> This qmhd_calc subroutine calculates occupation-weighted partition
!> function sums and their first and second (mixed) partial
!> derivatives wrt to x (the subset of auxiliary variables that are
!> relevant) and ln t for excited states.
!>
!> \param[in] ifpi_fit PARAMETERS NEED DOCUMENTATION
!>
subroutine qmhd_calc(ifpi_fit, weight_ion,&
     ifpl, ifapprox, eps_factor, nmin_rydberg, nmax,&
     tl, iz, rfactor, weight, neff, ell, x,&
     qmhd, qmhdt, qmhdx, qmhdt2, qmhdtx, qmhdx2)

  use mod_pi_fit, only: &
       pi_fitx_neutral_original, pi_fitx_neutral, pi_fitx_neutral_saumon,&
       pi_fitx_ion_original, pi_fitx_ion, pi_fitx_ion_saumon
  use mod_free_eos_constants, only: pi, echarge, bohr, boltzmann, ergspercmm1, rydberg
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! arguments
  logical, intent(in) :: ifpl, ifapprox
  integer, intent(in) :: ifpi_fit, iz
  integer, intent(in) :: nmin_rydberg, nmax
  real(fp_kind), intent(in) :: weight_ion, eps_factor, tl, rfactor,&
       weight(:), neff(:), ell(:), x(:)
  real(fp_kind), intent(out) :: qmhd, qmhdt, qmhdx(:), qmhdt2, qmhdtx(:), qmhdx2(:,:)

  ! internal variables

  ! minimum value of principal quantum number where approximations are used
  integer, parameter :: min_nmin_max = 3

  real(fp_kind), parameter :: const = 4._fp_kind*pi/3._fp_kind
  real(fp_kind), parameter :: rionconst3 = const*16._fp_kind*echarge*echarge*echarge*echarge*echarge*echarge

  logical ifneutral

  integer nx, nlevel
  integer ilevel, nmin_max

  real(fp_kind) t, rn2, energy, arg, summand, summandt, summandt2,&
       kn, occupation,&
       lnoccupation,&
       r, rnapprox, ellapprox
  real(fp_kind) neutral_factor, ionized_factor

  real(fp_kind), allocatable ::&
       qstar(:),&
       qstart(:),&
       qstarx(:,:),&
       qstart2(:),&
       qstartx(:,:),&
       qstarx2(:,:,:),&
       dlnoccupation(:),&
       dsummand(:),&
       dsummandt(:),&
       dsummand2(:,:)

  nlevel = size(weight)
  nx = size(x)

  ! Sanity checks:
  if(nlevel.ne.size(neff).or.nlevel.ne.size(ell))&
       error stop 'qmhd_calc: inconsistent dimensions for weight, neff, and ell'
  if(nx.ne.size(qmhdx).or.nx.ne.size(qmhdtx).or.nx.ne.size(qmhdx2,1).or.nx.ne.size(qmhdx2,2))&
       error stop 'qmhd_calc: inconsistent dimensions for x, qmhdx, qmhdtx, or qmhdx2'

  allocate(&
       dlnoccupation(nx),&
       dsummand(nx),&
       dsummandt(nx),&
       dsummand2(nx,nx))

  if(if_taint_allocated_real) then
     call taint_allocated_real(dlnoccupation)
     call taint_allocated_real(dsummand)
     call taint_allocated_real(dsummandt)
     call taint_allocated_real(dsummand2)
  endif

  ifneutral = iz.eq.1
  if(nmax.ge.nmin_rydberg) then
     ! n.b max(...) protects against using approximation
     ! inappropriately, (i.e., for n < 3)
     ! maximum value of mininum principal quantum number
     nmin_max =  max(min_nmin_max,nmin_rydberg)
     allocate(&
          qstar(nmin_max),&
          qstart(nmin_max),&
          qstarx(nx,nmin_max),&
          qstart2(nmin_max),&
          qstartx(nx,nmin_max),&
          qstarx2(nx,nx,nmin_max))

     if(if_taint_allocated_real) then
        call taint_allocated_real(qstar)
        call taint_allocated_real(qstart)
        call taint_allocated_real(qstarx)
        call taint_allocated_real(qstart2)
        call taint_allocated_real(qstartx)
        call taint_allocated_real(qstarx2)
     endif

     call qstar_calc(ifpi_fit, ifpl, .true. , ifneutral, ifapprox,&
          eps_factor, nmin_rydberg, nmax, iz, tl, x,&
          qstar, qstart, qstarx, qstart2, qstartx, qstarx2)
     qmhd = weight_ion*qstar(nmin_rydberg)
     qmhdt = weight_ion*qstart(nmin_rydberg)
     qmhdt2 = weight_ion*qstart2(nmin_rydberg)
     qmhdx(5) = weight_ion*qstarx(5,nmin_rydberg)
     qmhdtx(5) = weight_ion*qstartx(5,nmin_rydberg)
     qmhdx2(5,5) = weight_ion*qstarx2(5,5,nmin_rydberg)
     if(iz.eq.1) then
        qmhdx(1) = weight_ion*qstarx(1,nmin_rydberg)
        qmhdx(2) = weight_ion*qstarx(2,nmin_rydberg)
        qmhdx(3) = weight_ion*qstarx(3,nmin_rydberg)
        qmhdx(4) = weight_ion*qstarx(4,nmin_rydberg)
        qmhdtx(1) = weight_ion*qstartx(1,nmin_rydberg)
        qmhdtx(2) = weight_ion*qstartx(2,nmin_rydberg)
        qmhdtx(3) = weight_ion*qstartx(3,nmin_rydberg)
        qmhdtx(4) = weight_ion*qstartx(4,nmin_rydberg)
        qmhdx2(1,1) = weight_ion*qstarx2(1,1,nmin_rydberg)
        qmhdx2(2,1) = weight_ion*qstarx2(2,1,nmin_rydberg)
        qmhdx2(3,1) = weight_ion*qstarx2(3,1,nmin_rydberg)
        qmhdx2(4,1) = weight_ion*qstarx2(4,1,nmin_rydberg)
        qmhdx2(5,1) = weight_ion*qstarx2(5,1,nmin_rydberg)
        qmhdx2(2,2) = weight_ion*qstarx2(2,2,nmin_rydberg)
        qmhdx2(3,2) = weight_ion*qstarx2(3,2,nmin_rydberg)
        qmhdx2(4,2) = weight_ion*qstarx2(4,2,nmin_rydberg)
        qmhdx2(5,2) = weight_ion*qstarx2(5,2,nmin_rydberg)
        qmhdx2(3,3) = weight_ion*qstarx2(3,3,nmin_rydberg)
        qmhdx2(4,3) = weight_ion*qstarx2(4,3,nmin_rydberg)
        qmhdx2(5,3) = weight_ion*qstarx2(5,3,nmin_rydberg)
        qmhdx2(4,4) = weight_ion*qstarx2(4,4,nmin_rydberg)
        qmhdx2(5,4) = weight_ion*qstarx2(5,4,nmin_rydberg)
     endif
  else
     ! nmax cutoff allows no rydberg part of helium partition function
     qmhd = 0._fp_kind
     qmhdt = 0._fp_kind
     qmhdt2 = 0._fp_kind
     qmhdx(5) = 0._fp_kind
     qmhdtx(5) = 0._fp_kind
     qmhdx2(5,5) = 0._fp_kind
     if(iz.eq.1) then
        qmhdx(1) = 0._fp_kind
        qmhdx(2) = 0._fp_kind
        qmhdx(3) = 0._fp_kind
        qmhdx(4) = 0._fp_kind
        qmhdtx(1) = 0._fp_kind
        qmhdtx(2) = 0._fp_kind
        qmhdtx(3) = 0._fp_kind
        qmhdtx(4) = 0._fp_kind
        qmhdx2(1,1) = 0._fp_kind
        qmhdx2(2,1) = 0._fp_kind
        qmhdx2(3,1) = 0._fp_kind
        qmhdx2(4,1) = 0._fp_kind
        qmhdx2(5,1) = 0._fp_kind
        qmhdx2(2,2) = 0._fp_kind
        qmhdx2(3,2) = 0._fp_kind
        qmhdx2(4,2) = 0._fp_kind
        qmhdx2(5,2) = 0._fp_kind
        qmhdx2(3,3) = 0._fp_kind
        qmhdx2(4,3) = 0._fp_kind
        qmhdx2(5,3) = 0._fp_kind
        qmhdx2(4,4) = 0._fp_kind
        qmhdx2(5,4) = 0._fp_kind
     endif
  endif
  if(ifpi_fit.eq.0) then
     neutral_factor = pi_fitx_neutral_original
     ionized_factor = pi_fitx_ion_original
  elseif(ifpi_fit.eq.1) then
     neutral_factor = pi_fitx_neutral
     ionized_factor = pi_fitx_ion
  elseif(ifpi_fit.eq.2) then
     neutral_factor = pi_fitx_neutral_saumon
     ionized_factor = pi_fitx_ion_saumon
  else
     error stop 'qmhd_calc: bad ifpi_fit value'
  endif
  t = exp(tl)
  do ilevel = 1, nlevel
     ! n.b. checked that for helium energy levels, nint(neff) = n for
     ! *all* excited states. int(neff)+1 doesn't cut it because quantum
     ! defect sometimes slightly negative (level slightly too high).
     if(nint(neff(ilevel)).le.min(nmin_rydberg-1,nmax)) then
        rn2 = neff(ilevel)*neff(ilevel)
        ! energy below ionization continuum
        energy = ergspercmm1*rydberg*rfactor*real((iz*iz),fp_kind)/rn2
        arg = energy/(boltzmann*t)
        if(ifpl) then
           if(arg.gt.0.01_fp_kind) then
              ! lose a maximum of 4 significant digits
              summand = weight(ilevel)*(exp(arg) - (1._fp_kind + arg))
              summandt = -arg*weight(ilevel)*(exp(arg) - 1._fp_kind)
              summandt2 = arg*weight(ilevel)*(exp(arg)*(1._fp_kind+arg) - 1._fp_kind)
           else
              summand = 0.5_fp_kind*weight(ilevel)*arg*arg*(&
                   1._fp_kind + arg/3._fp_kind*(&
                   1._fp_kind + arg/4._fp_kind*(&
                   1._fp_kind + arg/5._fp_kind*(&
                   1._fp_kind + arg/6._fp_kind*(&
                   1._fp_kind + arg/7._fp_kind*(&
                   1._fp_kind + arg/8._fp_kind*(&
                   1._fp_kind + arg/9._fp_kind*(&
                   1._fp_kind + arg/10._fp_kind*(&
                   1._fp_kind + arg/11._fp_kind)))))))))
              summandt = -weight(ilevel)*arg*arg*(&
                   1._fp_kind + arg/2._fp_kind*(&
                   1._fp_kind + arg/3._fp_kind*(&
                   1._fp_kind + arg/4._fp_kind*(&
                   1._fp_kind + arg/5._fp_kind*(&
                   1._fp_kind + arg/6._fp_kind*(&
                   1._fp_kind + arg/7._fp_kind*(&
                   1._fp_kind + arg/8._fp_kind*(&
                   1._fp_kind + arg/9._fp_kind*(&
                   1._fp_kind + arg/10._fp_kind)))))))))
              summandt2 = weight(ilevel)*arg*arg*(&
                   2._fp_kind + arg/2._fp_kind*(&
                   3._fp_kind + arg/3._fp_kind*(&
                   4._fp_kind + arg/4._fp_kind*(&
                   5._fp_kind + arg/5._fp_kind*(&
                   6._fp_kind + arg/6._fp_kind*(&
                   7._fp_kind + arg/7._fp_kind*(&
                   8._fp_kind + arg/8._fp_kind*(&
                   9._fp_kind + arg/9._fp_kind*(&
                   10._fp_kind + arg*11._fp_kind/10._fp_kind)))))))))
           endif
        else
           summand = weight(ilevel)*exp(arg)
           summandt = -weight(ilevel)*arg*exp(arg)
           summandt2 = weight(ilevel)*arg*exp(arg)*(1._fp_kind+arg)
        endif
        ! quantum correction K_n see Hummer and Mihalas eq. 4.24
        if(neff(ilevel).le.3) then
           kn = 1._fp_kind
        else
           kn =&
                (16._fp_kind*rn2*(neff(ilevel) + 7._fp_kind/6._fp_kind))/&
                (3._fp_kind*(neff(ilevel)+1._fp_kind)*(neff(ilevel)+1._fp_kind)*&
                (rn2 + neff(ilevel) + 0.5_fp_kind))
        endif
        !          mhd occupation probability for neutral-ion and ion-ion interactions
        dlnoccupation(5) = -rionconst3*(ionized_factor*sqrt(real((iz),fp_kind)/kn)/energy)**3
        lnoccupation = x(5)*dlnoccupation(5)
        if(iz.eq.1) then
           ! first expression is "exact", but details of this part of
           ! occupation probability don't matter (first part much larger) so
           ! simplify to be consistent with simplified presentation
           ! r = neutral_factor*0.5d0*bohr*&
           !   (3.d0*rn2 - ell(ilevel)*(ell(ilevel)+1.d0))
           rnapprox = real((nint(neff(ilevel))),fp_kind)
           ellapprox = rnapprox - 1._fp_kind
           r = neutral_factor*0.5_fp_kind*bohr*&
                (3._fp_kind*rnapprox*rnapprox - ellapprox*(ellapprox+1._fp_kind))
           dlnoccupation(1) = -const
           dlnoccupation(2) = dlnoccupation(1)*r*3._fp_kind
           dlnoccupation(3) = dlnoccupation(2)*r
           dlnoccupation(4) = dlnoccupation(3)*r/3._fp_kind
           lnoccupation = lnoccupation +&
                dlnoccupation(1)*x(1) +&
                dlnoccupation(2)*x(2) +&
                dlnoccupation(3)*x(3) +&
                dlnoccupation(4)*x(4)
        endif
        occupation = exp(lnoccupation)
        summand = summand*occupation
        summandt = summandt*occupation
        summandt2 = summandt2*occupation
        dsummand(5) = summand*dlnoccupation(5)
        dsummandt(5) = summandt*dlnoccupation(5)
        dsummand2(5,5) = summand*dlnoccupation(5)*dlnoccupation(5)
        if(iz.eq.1) then
           dsummand(1) = summand*dlnoccupation(1)
           dsummand(2) = summand*dlnoccupation(2)
           dsummand(3) = summand*dlnoccupation(3)
           dsummand(4) = summand*dlnoccupation(4)
           dsummandt(1) = summandt*dlnoccupation(1)
           dsummandt(2) = summandt*dlnoccupation(2)
           dsummandt(3) = summandt*dlnoccupation(3)
           dsummandt(4) = summandt*dlnoccupation(4)
           dsummand2(1,1) = summand*dlnoccupation(1)*dlnoccupation(1)
           dsummand2(2,1) = summand*dlnoccupation(2)*dlnoccupation(1)
           dsummand2(3,1) = summand*dlnoccupation(3)*dlnoccupation(1)
           dsummand2(4,1) = summand*dlnoccupation(4)*dlnoccupation(1)
           dsummand2(5,1) = summand*dlnoccupation(5)*dlnoccupation(1)
           dsummand2(2,2) = summand*dlnoccupation(2)*dlnoccupation(2)
           dsummand2(3,2) = summand*dlnoccupation(3)*dlnoccupation(2)
           dsummand2(4,2) = summand*dlnoccupation(4)*dlnoccupation(2)
           dsummand2(5,2) = summand*dlnoccupation(5)*dlnoccupation(2)
           dsummand2(3,3) = summand*dlnoccupation(3)*dlnoccupation(3)
           dsummand2(4,3) = summand*dlnoccupation(4)*dlnoccupation(3)
           dsummand2(5,3) = summand*dlnoccupation(5)*dlnoccupation(3)
           dsummand2(4,4) = summand*dlnoccupation(4)*dlnoccupation(4)
           dsummand2(5,4) = summand*dlnoccupation(5)*dlnoccupation(4)
        endif
        qmhd = qmhd + summand
        qmhdt = qmhdt + summandt
        qmhdt2 = qmhdt2 + summandt2
        qmhdx(5) = qmhdx(5) + dsummand(5)
        qmhdtx(5) = qmhdtx(5) + dsummandt(5)
        qmhdx2(5,5) = qmhdx2(5,5) + dsummand2(5,5)
        if(iz.eq.1) then
           qmhdx(1) = qmhdx(1) + dsummand(1)
           qmhdx(2) = qmhdx(2) + dsummand(2)
           qmhdx(3) = qmhdx(3) + dsummand(3)
           qmhdx(4) = qmhdx(4) + dsummand(4)
           qmhdtx(1) = qmhdtx(1) + dsummandt(1)
           qmhdtx(2) = qmhdtx(2) + dsummandt(2)
           qmhdtx(3) = qmhdtx(3) + dsummandt(3)
           qmhdtx(4) = qmhdtx(4) + dsummandt(4)
           qmhdx2(1,1) = qmhdx2(1,1) + dsummand2(1,1)
           qmhdx2(2,1) = qmhdx2(2,1) + dsummand2(2,1)
           qmhdx2(3,1) = qmhdx2(3,1) + dsummand2(3,1)
           qmhdx2(4,1) = qmhdx2(4,1) + dsummand2(4,1)
           qmhdx2(5,1) = qmhdx2(5,1) + dsummand2(5,1)
           qmhdx2(2,2) = qmhdx2(2,2) + dsummand2(2,2)
           qmhdx2(3,2) = qmhdx2(3,2) + dsummand2(3,2)
           qmhdx2(4,2) = qmhdx2(4,2) + dsummand2(4,2)
           qmhdx2(5,2) = qmhdx2(5,2) + dsummand2(5,2)
           qmhdx2(3,3) = qmhdx2(3,3) + dsummand2(3,3)
           qmhdx2(4,3) = qmhdx2(4,3) + dsummand2(4,3)
           qmhdx2(5,3) = qmhdx2(5,3) + dsummand2(5,3)
           qmhdx2(4,4) = qmhdx2(4,4) + dsummand2(4,4)
           qmhdx2(5,4) = qmhdx2(5,4) + dsummand2(5,4)
        endif
     endif
  enddo
end subroutine qmhd_calc
