// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */
#ifdef MPIMODE
    #include "mpi.h"
#endif
#include "Constants.h"
#include <cmath>
#include "Utilities.h"
#include <valarray>
#include <cstdlib>
#include "Orbital_elements.h"
#include "Spline.h"
#include "Delay_brut.h"
#include "IO.h"
#include "Fittriple.h"
#include "Parameters.h"
#include "Diagnostics.h"

using namespace std ;





double Fittriple::Compute_lnposterior(int fractional) {
    /*
     * Compute the logarithm of the posterior probability ( at a factor ) .
     * If fractional is not 0 considers that the phase residuals are the fractional parts of the phases.
//      * If fractional is 0 the  phase residuals are the difference between the phase and the turn numbers in "turns".
    */
    bool nanflag = false;
    value_type pRA = previous_RA; // save values before they are changed by Initialize. See use below .
    value_type pDEC = previous_DEC;

    Initialize() ;


    long int i,j =0;
    value_type errtoe = pow(dix, -14);
    value_type npulsdir = 0.;
    value_type lnposterior = zero ;
    value_type * delay ;
    value_type * delay_geom ;
    value_type ** delay_geom_details = NULL ;
    value_type delay0;
    value_type maxdelay = zero;
    value_type delaymaxtemp = delaymax;//parameters.Pi * 16;
    long int ntisaround = 30 + 2*max(floor( 0.012  / dt_interp ),floor(delaymaxtemp / dt_interp )); // Number of interpolation steps to take around each toa to compute geometric delays. must be even.
    long int firsttis = 0L;

    value_type pospsr[3];
    value_type velpsr[3];

    delay = (value_type*) malloc( sizeof(value_type) * ninterp );
    delay_geom = (value_type*) malloc( sizeof(value_type) * ntisaround );
    if (getdelaysflag > 0)
    {
        delay_geom_details = (value_type**) malloc(sizeof(value_type*) * 6);
        for (i = 0 ; i < 6 ; i++) delay_geom_details[i] = (value_type*) malloc(sizeof(value_type) * ntoas);
    }


  // Recover proper motion and position from tempo
    for (i = 0; i < 3 ; i ++)
        {
            pospsr[i] = static_cast<value_type>(tempo2_psr[0].posPulsar[i]);
            velpsr[i] = static_cast<value_type>(tempo2_psr[0].velPulsar[i]) / static_cast<value_type>(100); // convertion from rad/century to rad/year.
        }



    if ( motion_changed == true)   // Reintegrate equations of motion only if necessary
    {
        Integrate_Allways() ;
            for (i = 0; i < ninterp ; ++i) 
            {
                sp[i][0] = states[i][0] * length ;
                sp[i][1] = states[i][1] * length ;
                sp[i][2] = states[i][2] * length ;
                si[i][0] = states[i][3] * length ;
                si[i][1] = states[i][4] * length ;
                si[i][2] = states[i][5] * length ;
                so[i][0] = states[i][6] * length ;
                so[i][1] = states[i][7] * length ;
                so[i][2] = states[i][8] * length ;
                sp[i][3] = states[i][extrashift + 9] * length / timescale ;
                sp[i][4] = states[i][extrashift + 10] * length / timescale ;
                sp[i][5] = states[i][extrashift + 11] * length / timescale ;              
                si[i][3] = states[i][extrashift + 12] * length / timescale ;
                si[i][4] = states[i][extrashift + 13] * length / timescale ;
                si[i][5] = states[i][extrashift + 14] * length / timescale ;
                so[i][3] = states[i][extrashift + 15] * length / timescale ;
                so[i][4] = states[i][extrashift + 16] * length / timescale ;
                so[i][5] = states[i][extrashift + 17] * length / timescale ;

            }
    }
   
   // test 
//         valarray<value_type> center_of_mass(3);
//         valarray<value_type> impulsion(3);
//         value_type energy;
//         int nbody = parameters.nextra +3;
//         valarray<value_type> sv(nbody*6);
//         for (i = 0 ; i < nbody*3; i++) sv[i] = states[states.size() -1][i] * length;
//         for (i = nbody*3 ; i < nbody*6; i++) sv[i] = states[states.size() -1][i]*length/timescale;
//         valarray<value_type> Ms(3+parameters.nextra);
//         Ms[0] = int_M0;
//         Ms[1] = int_M1;
//         Ms[2] = int_M2;
//         for (i = 0; i < parameters.nextra; i++)
//         {
//             Ms[i+3] = int_M_extra[i];
//         }
//         IntegralePrems3_1PN(sv, Ms,
//                          0.,
//                          int_Gg, int_gammabar, int_betabar,
//                          center_of_mass, impulsion,
//                          energy  );
//         impulsion /= int_M0 + int_M1 + int_M2 + Ms[3];
//         
//         printf("last CDMMMMMMMMMMMMMMM %.5Le %.5Le %.5Le \n\n", center_of_mass[0], center_of_mass[1],center_of_mass[2]);
//         printf("last VCDMMMMMMMMMMMMMM %.5Le %.5Le %.5Le \n\n", impulsion[0], impulsion[1],impulsion[2]);
//         printf("energyyyyyyyyyyyy %.12Le\n\n", energy);
//         Print_table(sv);
//         



// *** Uncomment the following and adjust the variables to test the numerical accuracy by integrating from a different reference time using an initial state vector found using a previous integration
//     printf("\n\n state vec at 100 000 : tinterp [100000] = %.19Le \n", tinterp[100000]);
//     for (i=0 ; i < 18 ; i++) printf("x0[%i] = %.19Le ;\n", i, states[100000][i]);

// end of test



    if (motion_changed == false && (pRA != parameters.RA || pDEC != parameters.DEC))
    {
            for (i = 0; i < ninterp ; ++i) {            // In RA or DEC changed, rotate back and forth the vectors without need for recomputation.
                // This could be otpimezed with only one rotation !
                Rotate_SSB_to_PSB(sp[i], pRA, pDEC);
                Rotate_PSB_to_SSB(sp[i], parameters.RA, parameters.DEC);
                Rotate_SSB_to_PSB(sp[i]+3, pRA, pDEC);
                Rotate_PSB_to_SSB(sp[i]+3, parameters.RA, parameters.DEC);

                Rotate_SSB_to_PSB(si[i], pRA, pDEC);
                Rotate_PSB_to_SSB(si[i], parameters.RA, parameters.DEC);
                Rotate_SSB_to_PSB(si[i]+3, pRA, pDEC);
                Rotate_PSB_to_SSB(si[i]+3, parameters.RA, parameters.DEC);

                Rotate_SSB_to_PSB(so[i], pRA, pDEC);
                Rotate_PSB_to_SSB(so[i], parameters.RA, parameters.DEC);
                Rotate_SSB_to_PSB(so[i]+3, pRA, pDEC);
                Rotate_PSB_to_SSB(so[i]+3, parameters.RA, parameters.DEC);
            }
    }



    // Compute Shapiro and aberration delays interpolation
        Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pospsr,
                 int_M0, int_M1, int_M2, parameters.f,
                 false, parameters.shapiro, parameters.aberration,
                 delay,  parameters.truefreq, freqshift, nanflag,  NULL
                ) ;
                
        if (nanflag == true) // if a delay is NaN, at least reject the location in parameter space
         {
             printf("\n Warning : delay computation returned nan, lnposterior will be -DBL_MAX. Printing parameters below : \n\n");
             parameters.Print();
             int mpirank=0;
#ifdef MPIMODE
             MPI_Comm_rank(MPI_COMM_WORLD , &mpirank);
#endif
             char nanparfile[100];
             sprintf(nanparfile, "/travail/gvoisin/4thbody-20210401/nanparfile-%d", mpirank);
             parameters.Save_parfile(nanparfile);
             return -DBL_MAX;
         }

         Spline delays_interpolation( tinterp, delay, ninterp );



   // Compute Einstein delay interpolation
         Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pospsr,
                 int_M0, int_M1, int_M2, parameters.f,
                 parameters.einstein, false, false,
                 delay,  parameters.truefreq, freqshift, nanflag, NULL
                ) ;
        
         if (nanflag == true) // if a delay is NaN, at least reject the location in parameter space
         {
             printf("\n Warning : delay computation returned nan, lnposterior will be -DBL_MAX. Printing parameters below : \n\n");
             parameters.Print();
             return -DBL_MAX;
         }
                
         Spline delay_einstein_i( tinterp, delay, ninterp);

  // Initialize Geometric delay interpolation
         Spline delay_geom_i (ntisaround);
//          if (getdelaysflag > 0)
//          {
             Spline delay_geom0_i (ntisaround);
             Spline delay_geom1_i (ntisaround);
             Spline delay_geom2_i (ntisaround);
             Spline delay_geom3_i (ntisaround);
             Spline delay_geom4_i (ntisaround);
             Spline delay_geom5_i (ntisaround);
//          }


  // Diagnostics
         if (get_global_interpolated_delays_flag == true)
         {
             for (i=0; i < diag_interpolated_delays_n ; i++)
             {
                 diag_einstein_interp[i] = static_cast<double>( delay_einstein_i(static_cast<value_type>( diag_interpolated_delays_times[i]), 0, ninterp) );
                 diag_shapiro_aberration_interp[i] = static_cast<double>( delays_interpolation(static_cast<value_type>(diag_interpolated_delays_times[i]) , 0, ninterp) );
             }
         }

         // Computation of margins around the interpolation node. This could be optimized (delaymax can become large because of proper motion and einstein delay, see def in Fittriple-init.cpp)
        long int marginmin = static_cast<long int>(max(floor( 0.012  / dt_interp ), floor((delaymaxtemp ) / dt_interp)) ) +2L ;// + 200L ; // We assume that the delays will never be greater than 0.1 days
        const long int marginsup = static_cast<long int>(floor( 0.012  / dt_interp ) ) + 2L; // Take into account the fact tinterp are bats which vary mostly from the roemer delay of the earth (0.012 ~ 10 ligth min)
        value_type ta, Te0, Te1 = zero;
        value_type test = 0.;



               test = 0. ; // ! test ! , supprimer la variable test dans delays à terme


   // Inversion of the timing formula for each toa

        for (i=0 ; i < ntoas ; ++i ){

            ta = toas[i] * parameters.DopplerF;

            firsttis = toa_in_interp[i] - ntisaround / 2 ;

            if ( tinterp[firsttis] > ta or tinterp[firsttis + ntisaround ] < ta ) printf("\n\n ******* Error !! %li %li %.19Le %.19Le \n\n " , firsttis, ntisaround, tinterp[firsttis], ta);


             Delays_Brut_geometric_local(tinterp + firsttis, toas[i], ntisaround,
                 sp + firsttis, si + firsttis, so + firsttis, pospsr, velpsr, parameters.posepoch, obs_ssb_BAT[i],
                 parameters.distance, parameters.distance1,
                 delay_geom, parameters.kopeikin, parameters.shklovskii, test, delay_geom_details
                )  ;

            delay_geom_i.ReSpline( tinterp + firsttis, delay_geom ) ;


            Te0 = ta ;

            delay0 = delays_interpolation( Te0 , (toa_in_interp[i] - marginmin),
                                           toa_in_interp[i] + marginsup) ;
            delay0 += delay_geom_i(Te0, ntisaround/2 - marginmin,
                                            ntisaround/2 + marginmin) ;

            Te1 = ta - delays_interpolation( ta - delay0,
                                             toa_in_interp[i] - marginmin, toa_in_interp[i] + marginsup ) ;
            Te1 -= delay_geom_i( ta - delay0,
                                              ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ;

            j = 0 ;
            while ( abs(Te1 - Te0 ) > errtoe and j < 100) 
            {
                ++j;
                Te0 = Te1 ;

                delay0 = delays_interpolation( Te0 , toa_in_interp[i] - marginmin,
                                            toa_in_interp[i] + marginsup ) ;
                delay0 += delay_geom_i(Te0, ntisaround/2 - marginmin,
                                            ntisaround/2 + marginmin) ;

                Te1 = ta - delays_interpolation( ta - delay0,
                                             toa_in_interp[i] - marginmin, toa_in_interp[i] + marginsup ) ;
                Te1 -= delay_geom_i( ta - delay0,
                                              ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ;

            }

            if (j == 100 ) {
                cout << " Iteration did not converge in _Update_delays_ ! " << endl;
                cout << " ninterp, ta, Te1, delays_interpolation( ta ): " << ninterp << ", " <<
                                                                       ta << ", " <<
                                                                       Te1 << ", " <<
                                                                       delays_interpolation( ta ,  toa_in_interp[i] - marginmin,
                                            toa_in_interp[i] + marginsup ) << endl;
            }


            toes[i] = Te1 - delay_einstein_i( Te1, toa_in_interp[i] - marginmin, toa_in_interp[i] + marginsup ) ;

            if (abs(ta - toes[i] ) > maxdelay ) maxdelay = huit *abs(ta - Te1 ); // reevaluate interpolation margins if necessary.

             if (getdelaysflag >  0) // extract delays for diagnostics
             {
                geomdelay[i] =  static_cast<double>(delay_geom_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
                nogeomdelay[i] = static_cast<double>(delays_interpolation( Te1, toa_in_interp[i] - marginmin,
                                            toa_in_interp[i] + marginsup ) );
                einsteindelay[i] = static_cast<double>(delay_einstein_i( toes[i] , toa_in_interp[i] - marginmin,
                                            toa_in_interp[i] + marginsup ) );
                delay_geom0_i.ReSpline( tinterp + firsttis, delay_geom_details[0] ) ;
                delay_geom1_i.ReSpline( tinterp + firsttis, delay_geom_details[1] ) ;
                delay_geom2_i.ReSpline( tinterp + firsttis, delay_geom_details[2] ) ;
                delay_geom3_i.ReSpline( tinterp + firsttis, delay_geom_details[3] ) ;
                delay_geom4_i.ReSpline( tinterp + firsttis, delay_geom_details[4] ) ;
                delay_geom5_i.ReSpline( tinterp + firsttis, delay_geom_details[5] ) ;
                geom_delay_details[i][0] =  static_cast<double>(delay_geom0_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
                geom_delay_details[i][1] =  static_cast<double>(delay_geom1_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
                geom_delay_details[i][2] =  static_cast<double>(delay_geom2_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
                geom_delay_details[i][3] =  static_cast<double>(delay_geom3_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
                geom_delay_details[i][4] =  static_cast<double>(delay_geom4_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
                geom_delay_details[i][5] =  static_cast<double>(delay_geom5_i( Te1,
                                               ntisaround/2 - marginmin, ntisaround/2 + marginmin ) ) ;
             }
        }

// -------------------------------------------------------------------------------------------
// ----------                      Special cases
// -------------------------------------------------------------------------------------------
//-------------- !!!! Double sinusoid model !!!!---------------------------------------------
        if (strcmp(parameters.specialcase,"RN_PL")==0)
        {
            value_type ds_f = 0.;
            value_type ds_A1 = 0.;
            value_type ds_phi1 = 0.;
            value_type ds_A2 = 0.;
            value_type ds_phi2 = 0.;
            value_type ds_A3 =0.;
            value_type ds_phi3 = 0.;
            value_type ds_A4 =0.;
            value_type ds_phi4 = 0.;
            value_type ds_A5 =0.;
            value_type ds_phi5 = 0.;
            value_type gam=0.;
            if (parameters.quadrupole.size() >= 5)
            {
                ds_f = parameters.quadrupole[0];
                ds_A1 = parameters.quadrupole[1]/(daysec * pow(dix, 6));
                ds_phi1 = parameters.quadrupole[2];
                ds_A2 = parameters.quadrupole[3]/(daysec * pow(dix, 6));
                ds_phi2 = parameters.quadrupole[4];
            }
            if (parameters.quadrupole.size() >= 7)
            {
                ds_A3 = parameters.quadrupole[5]/(daysec * pow(dix, 6));
                ds_phi3 = parameters.quadrupole[6];
            };
            if (parameters.quadrupole.size() >= 8) 
            {
                gam = parameters.quadrupole[7];
                ds_A2 = ds_A1/pow(2.,gam);
                ds_A3 = ds_A1/pow(3.,gam);
                ds_A4 = ds_A1/pow(4.,gam);
                if (parameters.quadrupole.size() >= 9) 
                {
                    ds_phi4=parameters.quadrupole[8];
                    ds_A5= ds_A1/pow(5.,gam);
                    if (parameters.quadrupole.size() >= 10) ds_phi5 = parameters.quadrupole[9];
                }
            }
            for (i = 0; i < ntoas ; ++i ) toes[i] -= (ds_A1 *sin(toas[i] *ds_f + ds_phi1) + ds_A2 *sin(2 * toas[i] *ds_f + ds_phi2) + ds_A3 *sin(3 * toas[i] *ds_f + ds_phi3) + ds_A4 *sin(4 * toas[i] *ds_f + ds_phi4) + ds_A5 *sin(5 * toas[i] *ds_f + ds_phi5));
        }

// -------------------------------------------------------------------------------------------

//-------------- !!!! Special: Circum-Ternary Keplerian orbit with period >> triple system and small mass (such that inversion of timing formula can be neglected  !!!!---------------------------------------------
        if (strcmp(parameters.specialcase,"Circum-ternary_Kepler")==0)
        {
            value_type dR_Kepler;
            value_type orbphase;
            value_type Ean, Eanold;
            value_type xorb, ecc, omp, Tp, orbfreq;
            
            orbfreq = parameters.quadrupole[0]; //  orbital frequency, day^-1 
            xorb = parameters.quadrupole[1]/daysec; // a sin(i) of triple CoM w.r.t extra body, day
            ecc = parameters.quadrupole[2]; // eccentricity
            omp = parameters.quadrupole[3]; // argument of peristron (rad)
            Tp = parameters.quadrupole[4] - parameters.treference - parameters.timeshift; // time of peristron passage (day)
            
            for (i = 0; i < ntoas ; ++i ) 
            {
                orbphase = (toas[i] - Tp ) * orbfreq * deuxpi; // Warning: only ok to use toe here as long as einstein delay not too large (can be cumulative). Otherwise should use coordinate time I.e. before correcting for einstein. 
                // Solving for Eccentric anomaly
                Ean = orbphase;
                Eanold = 0.;
                while (abs(Ean-Eanold) > pow(10., -14))
                {
                    Eanold = Ean;
                    Ean = orbphase + ecc * sin(Ean) ;
                }
                dR_Kepler = xorb * (sin(omp)*(cos(Ean) - 0*ecc) + sqrt(1.-ecc*ecc) * cos(omp)*sin(Ean));
                toes[i] -= dR_Kepler;
            }
        }

// -------------------------------------------------------------------------------------------
// ----------                      End of special cases
// -------------------------------------------------------------------------------------------
  // Compute the log of the posterior probability function
        //value_type spinPeriodmicrosec =  daysec * pow(dix, 6) / spinfreq ; // spin period in microseconds
        value_type mean = zero;
        value_type phase = zero;
        value_type weightsum = zero;
        value_type phase0 = toes[0] - parameters.treference ;
        phase0 = spinfreq /parameters.DopplerF * phase0 + undemi * spinfreq1/parameters.DopplerF*parameters.DopplerF * pow( phase0, 2 ) + parameters.dphase0;
        if (fractional == 0)
        {
            for (i = 0; i < ntoas ; ++i ) {
                phase = toes[i] - parameters.treference ;
                phase = spinfreq/parameters.DopplerF * phase + undemi * spinfreq1 /pow(parameters.DopplerF,2) * pow( phase, 2 ) ; // phase
                residuals[i] =  phase - phase0 - static_cast<value_type>( turns[i] );
                residuals[i] *= spinPeriodmicrosec ; // residuals in microseconds
                mean += residuals[i]*weights[i] ;
            }
        }
        else
        {
           phase0 = phase0 - dix - undemi ; // dix pour éviter que phase - phase0 ne soit négatif, un demi pour avoir la bonne partie entière (pas celle du tour d'avant)
           for (i = 0; i < ntoas ; ++i ) {
                phase = toes[i] - parameters.treference ;
                phase = spinfreq /parameters.DopplerF * phase + undemi * spinfreq1 /pow(parameters.DopplerF,2) * pow( phase, 2 ) ;
                residuals[i] =  (phase - phase0) - floor( phase - phase0 ) - undemi    ;
                residuals[i] *= spinPeriodmicrosec ;
                mean += residuals[i] * weights[i] ;
            }
        }

        // Remove weighted mean from residuals. This is equivalent to marginalizing over 
        // a constant offset. 
	if (remove_mean == true)
        	for (i = 0; i < ntoas ; ++i ) residuals[i] -= mean;

        // Compute the chi2
        for (i = 0; i < ntoas ; ++i ) {
            lnposterior += pow(  ( residuals[i]) / errors[i], 2 );
        }
        // Compute the log of the posterior probability density up to a constant and divided by ntoas
        lnposterior = -( 0.5*lnposterior / (parameters.efac* parameters.efac *ntoas) + log(parameters.efac) );

        // To prevent fits from going to negative quadrupole[0]
//         if (parameters.quadrupole.size() >0 and parameters.quadrupole[0] < 0.) lnposterior -= exp(-parameters.quadrupole[0]*5000.);
//         lnposterior /= ntoas ; //  - Quasi-reduced chi2


        free(delay);
        free(delay_geom);
        if (getdelaysflag > 0) 
        {
            for (i = 0 ; i < 6 ; i++) free(delay_geom_details[i]);
            free(delay_geom_details);
        }

        int p = 0;
        double pshift = 0.;
        if (tracker == true)
        {
            callnumber++;
            printf("Call %i : \n", callnumber);
            for (int i=0 ; i < parameters.fitted_parameters.size() ; i++ ) {
                p = parameters.fitted_parameters[i];
                pshift = abs(parameters.Parameters(p) - parameters.parameters_ini[p]) / parameters.Parameter_scale(p) ;
                printf(" %i:%.4f ", p, pshift);
            }
            printf("\n Log posterior probability  = %.15Le \n\n", lnposterior );
        }

        return double(lnposterior);
}
