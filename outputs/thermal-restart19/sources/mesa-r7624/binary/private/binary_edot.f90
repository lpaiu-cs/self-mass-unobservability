! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton and Joris Vos
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

    module binary_edot
    
    use const_def
    use star_lib
    use star_def
    use binary_def
    use binary_tides
    use crlibm_lib

    implicit none

    contains

    real(dp) function get_edot(b) result(edot)

       type (binary_info), pointer :: b
       real(dp) :: edot_tidal, edot_enhance

       if (b% eccentricity == 0d0) then
           b% edot_tidal = 0d0
           b% edot_enhance = 0d0
           b% edot = 0d0
           edot = 0d0
       else
           ! tidal circularisation
           b% edot_tidal = get_edot_tidal(b)
           
           if (b% edot_tidal < -b% max_abs_edot_tidal) then
              b% edot_tidal = -b% max_abs_edot_tidal
           end if
           
           ! eccentricity enhancement
           if (b% use_eccentricity_enhancement) then
              b% edot_enhance = edot_enhancement_Isotropic(b)
           else
              b% edot_enhance = 0d0
           end if
           
           if (b% edot_enhance > b% max_abs_edot_enhance) then
              b% edot_enhance = b% max_abs_edot_enhance
           end if
           
           edot = b% edot_tidal + b% edot_enhance
           
           b% edot = edot
           
       end if

    end function get_edot

    ! ==========================================
    ! edot TIDAL
    ! ==========================================
    real(dp) function get_edot_tidal(b) result(edot_tidal)
       type (binary_info), pointer :: b
       edot_tidal = 0d0
       if (.not. b% do_tidal_circ) return

       if (b% circ_type_1 == "Hut_conv") then
          edot_tidal = edot_tidal_Hut(b, b% s1, .true.)
       else if (b% circ_type_1 == "Hut_rad") then
          edot_tidal = edot_tidal_Hut(b, b% s1, .false.)
       else
          write(*,*) "Unrecognized circ_type_1", b% circ_type_1
       end if
       if (b% evolve_both_stars) then
          if (b% circ_type_2 == "Hut_conv") then
             edot_tidal = edot_tidal + edot_tidal_Hut(b, b% s2, .true.)
          else if (b% circ_type_2 == "Hut_rad") then
             edot_tidal = edot_tidal + edot_tidal_Hut(b, b% s2, .false.)
          else
             write(*,*) "Unrecognized circ_type_2", b% circ_type_2
          end if
       end if
    end function get_edot_tidal

    real(dp) function edot_tidal_Hut(b, s , has_convective_envelope) result(edot_tidal)
        type (binary_info), pointer :: b
        type (star_info), pointer :: s
        logical, intent(in) :: has_convective_envelope
        real(dp) :: m, porb, r_phot, osep, qratio, omega_s, omega_sync

        edot_tidal = 0d0

        porb = b% period
        omega_sync = 2*pi/b% period
        omega_s = s% omega_avg_surf
        osep = b% separation

        qratio = b% m(b% a_i) / b% m(b% d_i)
        if (is_donor(b, s)) then
           m = b% m(b% d_i)
           r_phot = b% r(b% d_i)
        else
           qratio = 1.0/qratio
           m = b% m(b% a_i)
           r_phot = b% r(b% a_i)
        end if

        ! eq. (10) of Hut, P. 1981, A&A, 99, 126
        edot_tidal = -27.0d0*qratio*(1+qratio)*pow8(r_phot/osep) &
            * b% eccentricity*pow_cr(1-b% eccentricity**2,6.5d0)*b% Ftid_1
        ! add multiplication by (k/T), eq. (29) of Hurley et al. 2002
        edot_tidal = edot_tidal*k_div_T(b, s, has_convective_envelope)
        ! add terms dependant on omega
        edot_tidal = edot_tidal*(f3(b% eccentricity) - &
            11d0/18d0 * omega_s / omega_sync * f4(b% eccentricity) * &
            pow_cr(1-b% eccentricity**2,1.5d0))
    
    end function edot_tidal_Hut
    
    ! ==========================================
    ! Edot MASS LOSS
    ! ==========================================
    
    real(dp) function edot_enhancement_Isotropic(b) result(edot_enhance)
        type (binary_info), pointer :: b
        integer :: i
        real(dp) :: de, Mtot, xfer
        real(dp), DIMENSION(b% anomaly_steps):: e1, e2, e3, edot
        
        e1 = (b% eccentricity + cos(b% theta_co)) 
        e2 = 2 * cos(b% theta_co) + b% eccentricity + b% eccentricity * cos(b% theta_co)**2
        e3 = b% eccentricity * sin(b% theta_co)**2
        
!         xfer = min(b% wind_xfer_fraction, b% xfer_fraction)
        Mtot = b% m(1) + b% m(2) ! total mass in gr
        
        edot = - b% mdot_donor_theta / Mtot * e1 !-&
!                b% mdot_donor_theta * xfer / b% m(b% a_i) * (b% m(b% d_i) / Mtot *&
!                ((b% m(b% a_i)**2 / b% m(b% d_i)**2 - 1 ) * e2 - e3 ))
        
        !integrate to get total eccentricity enhancement
         de = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            de = de + 0.5 * (edot(i-1) + edot(i)) * (b% time_co(i) - b% time_co(i-1)) 
         end do
         
         edot_enhance = de
    
    end function edot_enhancement_Isotropic
    

    end module binary_edot

