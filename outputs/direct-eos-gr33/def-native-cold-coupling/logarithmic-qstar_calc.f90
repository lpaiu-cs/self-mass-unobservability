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

! calculate internal partition function (referred to the ionized
! state) and derivatives using Rydberg energy levels.
!
! The partition function sum is taken for the range
! of minimum principal quantum numbers (given by nmin to nmin_max)
! to infinity of 2 n^2 exp(a/n^2) times
! planck_larkin occupation probability (ifpl true) times
! mhd occupation probability (ifmhd true).
! The planck-larkin occupation probability is given by
!   (1 - exp(-a)*(1+a)).
! The *ln* of the mhd occupation probability is given by
!   -[b(1) + b(2)*n^2*(1+g(n)) + b(3)*n^4*(1+g(n))^2 +
!   b(4)*n^6*(1+g(n))^3 + b(5)*n^(15/2)*(1+h(n))].
!   g(n) = 1/(2n).
!   h(n) = (16/(3*n*K_n))^(3/2) - 1
!        = (16/(3*n))^(3/2) - 1 [for n <= 3]
!        = ((n+1)/n)^3*((n^2 + n + 1/2)/(n^2 + 7/6*n))^(3/2) - 1 [for n >=3]
!   a = c2*R*iz*2/t is calculated internally and
!   the b vector is also calculated internally from the x vector and iz.

! input quantities:
! ifpi_fit = 2, use best fit to Saumon table
! ifpi_fit = 1, use best fit to opal table + extensions
! ifpi_fit = 0, use best fit to original MDH table.
! ifpl (logical) controls whether to use planck-larkin occupation probability.
! ifmhd (logical) controls whether to use mhd occupation probability.
! ifneutral (logical) controls when b(1) through b(4) are employed in the
!   occupation probability calculation.
! ifapprox (logical) controls whether to approximate sum (preferred
!   because much quicker and only 1.d-4 errors) or do actual sum
!   (more exact to test the approximations).
! eps_factor is a factor used to help terminate the principal
!   quantum number sum. eps_factor = exp(-c2*min(chi)/t), where
! min(chi) = the minimum ionization potential for all species
!   with the same a value (i.e. charge).  This definition means
!   that at the limit summand contributes approximately eps (see
!   parameter below) to ln(1+q_excited/q_ground).  With exponential
!   cutoff (MHD), this insures maximum cutoff errors are of order eps.
!   if not MHD (i.e., only Planck-Larkin occupation probability) the
!   remainder of the sum is n*summand.  With nmax of order 300000
!   this means total relative error is of order 1.d-4 which is still
!   fine.
! nmin to nmin_max is the range of minimum principal quantum numbers required.
! nmax is the maximum principal quantum number.  In Planck-Larkin case
!   obtain relative errors of 1.d-5 for nmin of 3 if nmax = 300000.
! iz is the charge number of the *ion* of the species.
! tl is ln T.
! x(nx), where nx = 5 are the auxiliary variables which
!   help to determine b and the mhd occupation probability.
!   n.b. index is *not* in auxiliary variable order and it is the
!   calling routine's responsibility to put x into the following order
!   (and also to reorder the derivatives appropriately):
!   using the notation of the paper on the RHS we have
!   x(i) = sigma^PI_{4-i} for i = 1,4 and
!   x(5) = sigma^PI_Z.

! output quantities:
! qstar(nmin_max), qstart(nmin_max), qstarx(nx,nmin_max)
! qstart2(nmin_max), qstartx(nx,nmin_max), qstarx2(nx,nx,nmin_max) are
!   the partition function plus partial derivatives
!   wrt *tl* and x.

!> This qstar_calc subroutine calculates occupation-weighted partition
!> function sums and their first and second (mixed) partial
!> derivatives wrt to x (the subset of auxiliary variables that are
!> relevant) and ln t for Rydberg excited states.
!>
!> \param[in] ifpi_fit PARAMETERS NEED DOCUMENTATION
!>
subroutine qstar_calc(ifpi_fit,&
     ifpl, ifmhd, ifneutral, ifapprox,&
     eps_factor, nmin, nmax, iz, tl, x,&
     qstar, qstart, qstarx, qstart2, qstartx, qstarx2, logscale)

  use mod_pi_fit, only: &
       pi_fitx_neutral_original, pi_fitx_neutral, pi_fitx_neutral_saumon,&
       pi_fitx_ion_original, pi_fitx_ion, pi_fitx_ion_saumon
  use mod_free_eos_constants, only: bohr, c2, echarge, electron_mass, ergspercmm1, h_mass, pi, rydberg
  use mod_free_eos_debug, only: if_taint_allocated_real, taint_allocated_real

  ! Arguments
  integer, intent(in) :: ifpi_fit, nmin, nmax, iz
  logical, intent(in) :: ifpl, ifmhd, ifneutral, ifapprox
  real(fp_kind), intent(in) :: eps_factor, tl, x(:)
  real(fp_kind), intent(out) ::&
       qstar(:), qstart(:),  qstart2(:),&
       qstarx(:,:), qstartx(:,:), qstarx2(:,:,:)

  real(fp_kind), intent(out), optional :: logscale(:)

  ! Internal variables
  integer, parameter :: nb = 5

  real(fp_kind), parameter :: rionconst = sqrt(3._fp_kind/16._fp_kind)*echarge*echarge/ergspercmm1/rydberg
  real(fp_kind), parameter :: const = 4._fp_kind*pi/3._fp_kind

  integer nmax_reached, n, ib, nx, nmin_max, ntarget

  real(fp_kind) a, rnuconst, neutral_factor, ionized_factor

  real(fp_kind), allocatable :: b(:), bconst(:), sq(:), st(:), sx(:,:), stt(:), stx(:,:), sxx(:,:,:)

  nx = size(x)
  nmin_max = size(qstar)

  ! Sanity checks:
  if(iz.lt.1) error stop 'qstar_calc: bad iz value'
  if(nx.ne.nb) error stop 'qstar_calc: bad x dimension'
  if(&
       nx.ne.size(qstarx,1).or.&
       nx.ne.size(qstartx,1).or.&
       nx.ne.size(qstarx2,1).or.&
       nx.ne.size(qstarx2,2))&
       error stop 'qstar_calc: inconsistent nx dimensions for x, qstarx, qstartx, or qstarx2'
  if(&
       nmin_max.ne.size(qstart).or.&
       nmin_max.ne.size(qstart2).or.&
       nmin_max.ne.size(qstarx,2).or.&
       nmin_max.ne.size(qstartx,2).or.&
       nmin_max.ne.size(qstarx2,3))&
       error stop 'qstar_calc: inconsistent nmin_max dimensions for qstar, qstart, qstart2, qstarx, qstartx, or qstarx2'
  if(nmin.lt.1.or.nmin.gt.nmin_max) error stop 'qstar_calc: bad nmin'

  ! n.b. rnu = rnuconst*n**2*(1+g(n))
  ! At the moment, we only calculate "neutral" radii for neutrals
  ! and H2+, but because of that latter we must include the Z=iz factor
  ! (see p. 117 of Condon and Shortly, "The Theory of Atomic Spectra").

  rnuconst = bohr/real((iz),fp_kind)

  allocate(b(nb), bconst(nb))

  if(if_taint_allocated_real) then
     call taint_allocated_real(b)
     call taint_allocated_real(bconst)
  endif

  ! n.b. rion^3 = 16.d0*(rionconst*iz**(-9/2)*n**7.5*(1+h(n))
  if(ifpi_fit.eq.0) then
     neutral_factor = pi_fitx_neutral_original*rnuconst
     ionized_factor = 16._fp_kind*(pi_fitx_ion_original*rionconst)**3
  elseif(ifpi_fit.eq.1) then
     neutral_factor = pi_fitx_neutral*rnuconst
     ionized_factor = 16._fp_kind*(pi_fitx_ion*rionconst)**3
  elseif(ifpi_fit.eq.2) then
     neutral_factor = pi_fitx_neutral_saumon*rnuconst
     ionized_factor = 16._fp_kind*(pi_fitx_ion_saumon*rionconst)**3
  else
     error stop 'qstar_calc: bad ifpi_fit value'
  endif
  bconst(1) = const
  bconst(2) = 3._fp_kind*const*neutral_factor
  bconst(3) = bconst(2)*neutral_factor
  bconst(4) = bconst(3)*neutral_factor/3._fp_kind
  bconst(5) = const*ionized_factor/(real((iz),fp_kind)**4*sqrt(real((iz),fp_kind)))

  a = c2*rydberg*real((iz*iz),fp_kind)*exp(-tl)
  ! convert rydberg(inf) to rydberg (finite) assuming:
  if(iz.eq.1) then
     ! iz = 1 dominated by hydrogen,
     a = a/(1._fp_kind + electron_mass/h_mass)
     bconst(5) = bconst(5)*(1._fp_kind + electron_mass/h_mass)**3
  elseif(iz.eq.2) then
     ! iz = 2 dominated by He,
     a = a/(1._fp_kind + electron_mass/(4._fp_kind*h_mass))
     bconst(5) = bconst(5)*(1._fp_kind + electron_mass/(4._fp_kind*h_mass))**3
  elseif(iz.gt.2) then
     ! and iz > 2 dominated by N.
     a = a/(1._fp_kind + electron_mass/(14._fp_kind*h_mass))
     bconst(5) = bconst(5)*(1._fp_kind + electron_mass/(14._fp_kind*h_mass))**3
  else
     error stop 'qstar_calc: bad iz'
  endif
  if(ifmhd) then
     if(ifneutral) then
        b = bconst*x
     else
        b(nb) = bconst(nb)*x(nb)
     endif
  endif
  if(present(logscale)) logscale=0._fp_kind
  qryd_logscale=0._fp_kind
  if(present(logscale).and.a/real(nmin*nmin,fp_kind).gt.680._fp_kind) then
     allocate(sq(nmin_max),st(nmin_max),sx(nb,nmin_max),stt(nmin_max),stx(nb,nmin_max),sxx(nb,nb,nmin_max))
     qstar=0._fp_kind;qstart=0._fp_kind;qstarx=0._fp_kind
     qstart2=0._fp_kind;qstartx=0._fp_kind;qstarx2=0._fp_kind
     do ntarget=nmin,nmin_max
        logscale(ntarget)=max(0._fp_kind,a/real(ntarget*ntarget,fp_kind)-400._fp_kind)
        qryd_logscale=logscale(ntarget)
        call qryd_calc(ifpl,ifmhd,ifneutral,eps_factor,ntarget,nmax,a,b,sq,st,sx,stt,stx,sxx,nmax_reached)
        qstar(ntarget)=sq(ntarget);qstart(ntarget)=st(ntarget);qstarx(:,ntarget)=sx(:,ntarget)
        qstart2(ntarget)=stt(ntarget);qstartx(:,ntarget)=stx(:,ntarget);qstarx2(:,:,ntarget)=sxx(:,:,ntarget)
     enddo
     qryd_logscale=0._fp_kind
  elseif(ifapprox) then
     call qryd_approx(ifpl, ifmhd, ifneutral, eps_factor, nmin, nmax, a, b,&
          qstar, qstart, qstarx, qstart2, qstartx, qstarx2, nmax_reached)
  else
     call qryd_calc(ifpl, ifmhd, ifneutral, eps_factor, nmin, nmax, a, b,&
          qstar, qstart, qstarx, qstart2, qstartx, qstarx2, nmax_reached)
  endif
  do n = nmin, nmin_max
     ! convert from b derivative to x derivative.
     if(ifmhd) then
        if(ifneutral) then
           qstarx(:,n) = qstarx(:,n)*bconst
           do ib = 1, nb
              qstarx2(ib:nb,ib,n) = qstarx2(ib:nb,ib,n)*bconst(ib:nb)*bconst(ib)
           enddo
           ! convert from a derivatives to tl derivatives.
           ! a = c2*rydberg*real((iz),fp_kind)*exp(-tl)
           qstartx(:,n) = -a*qstartx(:,n)*bconst
        else
           qstarx(nb,n) = qstarx(nb,n)*bconst(nb)
           qstarx2(nb,nb,n) = qstarx2(nb,nb,n)*bconst(nb)*bconst(nb)
           ! convert from a derivatives to tl derivatives.
           ! a = c2*rydberg*real((iz),fp_kind)*exp(-tl)
           qstartx(nb,n) = -a*qstartx(nb,n)*bconst(nb)
        endif
     endif
     ! convert from a derivatives to tl derivatives.
     ! a = c2*rydberg*real((iz),fp_kind)*exp(-tl)
     qstart2(n) = a*(qstart(n) + a*qstart2(n))
     qstart(n) = -a*qstart(n)
  enddo
end subroutine qstar_calc
