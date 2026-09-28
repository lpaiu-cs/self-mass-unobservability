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
!> This module provides variables that are used to help communicate between
!> the module procedures defined by the mod_exchange module.

module mod_master_exchange_data
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  real(fp_kind), save, public :: last_fl_input=0._fp_kind, last_tl_input=0._fp_kind
  public&
       iforder, ifstart,&
       dpsidf, dpsidf2,&
       dpsiprimedf, dpsiprimedf2,&
       fex, fexf, fext, fexf2, fexft, fext2, n_e, t,&
       fexprime, fexprimef, fexprimet, fexprimef2, fexprimeft, fexprimet2,&
       flprimef, flprimet, flprimef2, flprimeft, flprimet2,&
       muex, muexf, muext,&
       muex2, muex2f, muex2t,&
       p_e, pstarprime, psiprime, psi

  integer iforder, ifstart
  real(fp_kind)&
       dpsidf, dpsidf2,&
       dpsiprimedf, dpsiprimedf2,&
       fex, fexf, fext, fexf2, fexft, fext2, n_e, t,&
       fexprime, fexprimef, fexprimet, fexprimef2, fexprimeft, fexprimet2,&
       flprimef, flprimet, flprimef2, flprimeft, flprimet2,&
       muex, muexf, muext,&
       muex2, muex2f, muex2t,&
       p_e, pstarprime(9), psiprime, psi

  data ifstart/1/

  save&
       iforder, ifstart,&
       dpsidf, dpsidf2,&
       dpsiprimedf, dpsiprimedf2,&
       fex, fexf, fext, fexf2, fexft, fext2, n_e, t,&
       fexprime, fexprimef, fexprimet, fexprimef2, fexprimeft, fexprimet2,&
       flprimef, flprimet, flprimef2, flprimeft, flprimet2,&
       muex, muexf, muext,&
       muex2, muex2f, muex2t,&
       p_e, pstarprime, psiprime, psi
end module mod_master_exchange_data

!> This module provides module procedures to help determine the
!> non-ideal exchange component of the free-energy and all relevant
!> derivatives of same.
module mod_exchange
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public master_exchange, exchange_free, exchange_end, exchange_pressure, exchange_gcpf
contains
#include "master_exchange.f90"
#include "exchange_gcpf.f90"
#include "f_psi.f90"
#include "exct_calc.f90"
#include "fp12_calc.f90"
#include "exchange_coeff.f90"
end module mod_exchange
