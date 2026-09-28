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
!> the module procedures defined by the mod_coulomb module.
module mod_master_coulomb_data
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public&
       dmucoulomb0, dmucoulomb2, dmucoulombf, dmucoulombt, dmucoulomb,&
       fcft0f, fcft0t, fcft2f, fcft2t, fcftxf, fcftxt,&
       fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0,&
       fcoulomb22, fcoulomb2t, fcoulomb2,&
       fcoulombtt, fcoulombtx, fcoulombt, fcoulombx, fcoulomb,&
       ifstart,&
       sum0af, sum0ane, sum0atl,&
       sum2af, sum2ane, sum2atl,&
       theta_eff, theta_eft, theta_ef, thetatf, thetatt,&
       xf, xt, x

  real(fp_kind) ::&
       dmucoulomb0, dmucoulomb2, dmucoulombf, dmucoulombt, dmucoulomb,&
       fcft0f, fcft0t, fcft2f, fcft2t, fcftxf, fcftxt,&
       fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0,&
       fcoulomb22, fcoulomb2t, fcoulomb2,&
       fcoulombtt, fcoulombtx, fcoulombt, fcoulombx, fcoulomb,&
       sum0af, sum0ane, sum0atl,&
       sum2af, sum2ane, sum2atl,&
       theta_eff, theta_eft, theta_ef, thetatf, thetatt,&
       xf, xt, x
  integer ifstart
  data ifstart/1/

  save&
       dmucoulomb0, dmucoulomb2, dmucoulombf, dmucoulombt, dmucoulomb,&
       fcft0f, fcft0t, fcft2f, fcft2t, fcftxf, fcftxt,&
       fcoulomb00, fcoulomb02, fcoulomb0t, fcoulomb0,&
       fcoulomb22, fcoulomb2t, fcoulomb2,&
       fcoulombtt, fcoulombtx, fcoulombt, fcoulombx, fcoulomb,&
       ifstart,&
       sum0af, sum0ane, sum0atl,&
       sum2af, sum2ane, sum2atl,&
       theta_eff, theta_eft, theta_ef, thetatf, thetatt,&
       xf, xt, x

end module mod_master_coulomb_data

!> This module provides module procedures to help determine the
!> non-ideal Coulomb component of the free-energy and all relevant
!> derivatives of same.
module mod_coulomb
  use mod_free_eos_types, only: fp_kind
  implicit none
  private
  public master_coulomb, master_coulomb_free, master_coulomb_end, master_coulomb_pressure, coulomb_adjust
contains
#include "master_coulomb.f90"
#include "coulomb.f90"
#include "coulomb_adjust.f90"
#include "tau_calc.f90"
#include "pteh_theta.f90"
#include "lnx_calc.f90"
end module mod_coulomb
