// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include "Diagnostics.h"
#include "Fittriple.h"
#include <iostream>
#include "Delay_brut.h"
#include "Spline.h"
#include "Utilities.h"
#include <time.h>


using namespace std;


int Fittriple::Check_tinterp_grid(double tol, double fact)
{
    clock_t tclock;
    int interpstepsini =  parameters.interpsteps_per_period_i ;
    int interpsteps = 0;
    float marginini = 1.;
    float margin = 1.;
    double lnp = zero;
    double oldlnp = 0.;

    tclock = clock();
    lnp = Compute_lnposterior(0);
    tclock = clock() - tclock;
    printf("Computation time with initial grid resolution = %.2f seconds\n",((float)tclock)/CLOCKS_PER_SEC);

    cout  << "\n\n********************* Entering Check_tinterp_grid *********************" << endl << endl ;

    cout << endl << "    *** Checking for grid resolution with " << tol << " relative chi2 accuracy when resolution is multiplied by "<< fact << endl << endl;

    printf(" interpsteps_per_period_i = %i ; Log posterior = %.15f \n" , parameters.interpsteps_per_period_i, lnp) ;

    forcerecomputeinterp = true ;
    do {
        oldlnp = lnp ;
        parameters.interpsteps_per_period_i = floor( parameters.interpsteps_per_period_i * fact ) ;
        tclock = clock();
        lnp = Compute_lnposterior(0);
        tclock = clock() - tclock;
        printf(" interpsteps_per_period_i = %i ; Log posterior = %.15f ; Computation time < %.2fs \n" , parameters.interpsteps_per_period_i, lnp, ((float)tclock)/CLOCKS_PER_SEC) ;
    } while ( abs(lnp - oldlnp) / abs(lnp) > tol ) ;
    interpsteps = parameters.interpsteps_per_period_i ;
    printf("  Final interpsteps_per_period_i = %i -> Accuracy reached for at least interpsteps_per_period_i = %i \n\n", interpsteps , static_cast<int>(floor(interpsteps/fact)));
    printf("  Going back to initial value... \n");
    parameters.interpsteps_per_period_i = interpstepsini ;
    lnp = Compute_lnposterior(0);

     cout << endl << "    *** Checking for grid margin with " << tol << " relative chi2 accuracy when interp_margin is multiplied by "<< fact << endl << endl;
     marginini = parameters.interp_margin;
    do {
        oldlnp = lnp ;
        parameters.interp_margin = parameters.interp_margin * fact ;
        cout << parameters.interp_margin << endl;
        tclock = clock();
        lnp = Compute_lnposterior(0);
        tclock = clock() - tclock;
        printf(" interp_margin = %f ; Log posterior = %.15f ; Computation time < %.2fs\n" , parameters.interp_margin, lnp, ((float)tclock)/CLOCKS_PER_SEC) ;
    } while ( abs(lnp - oldlnp) / abs(lnp) > tol )  ;
    margin = parameters.interp_margin;

     printf("  Final interp_margin = %f -> Accuracy reached for at least interp_margin = %f \n", margin, margin / fact);
     printf("  Going back to initial value... \n");
     parameters.interp_margin = marginini ;
     Compute_lnposterior(0);

     forcerecomputeinterp = false;
     cout << "\n********************* Leaving  Check_tinterp_grid *********************\n\n";

     return interpsteps ;
}


void Fittriple::Compute_integrals_of_motion()
{
    center_of_mass_positions.resize(ninterp);
    center_of_mass_impulsions.resize(ninterp);
    energies.resize(ninterp);

    if (integrator_type == 3 or integrator_type==1)//(integrator_type == 1)
    {
        cout << endl << "Computing 1PN integrals of motion" << endl ;
        valarray<value_type> sv(6* parameters.nbodies_plus_extra);
        valarray<value_type> Ms(parameters.nbodies_plus_extra);
        Ms[0] = int_M0;
        Ms[1] = int_M1;
        Ms[2] = int_M2;
        for (int i = 0; i < parameters.nextra; i++) Ms[i+3] = int_M_extra[i];
        for (long int i = 0 ; i < ninterp ; ++i)
        {
            for (int j = 0 ; j < parameters.nbodies_plus_extra *3 ; j++) sv[j] = states[i][j] * length;
            for (int j = parameters.nbodies_plus_extra * 3; j < parameters.nbodies_plus_extra * 6; j++) sv[j] = states[i][j] * length/timescale;
        
            IntegralePrems3_1PN( sv, Ms, 0.,
                                    int_Gg, int_gammabar, int_betabar,
                                    center_of_mass_positions[i], center_of_mass_impulsions[i],
                                    energies[i]  ) ;
        }
    }
//    {
//         cout << endl << "Computing 1PN integrals of motion bloum" << endl ;
//         for (long int i = 0 ; i < ninterp ; ++i)
//             IntegralePrems3_1PN( sp[i], si[i], so[i],
//                                     int_M0, int_M1, int_M2, 0.,
//                                      int_Gg, int_gammabar, int_betabar,
//                                     center_of_mass_positions[i], center_of_mass_impulsions[i],
//                                     energies[i]  ) ;
//     } 
    else if (integrator_type == 2)
    {
        cout << endl << "Computing 1PN integrals of motion with inner WD quadrupole" << endl ;
        for (long int i = 0 ; i < ninterp ; ++i)
            IntegralePrems3_1PN( sp[i], si[i], so[i],
                                    int_M0, int_M1, int_M2, int_quadrupole[0],
                                     int_Gg, int_gammabar, int_betabar,
                                    center_of_mass_positions[i], center_of_mass_impulsions[i],
                                    energies[i]  ) ;
    }

    else if (integrator_type == 0)
    {
        cout << endl << "Computing Newtonian integrals of motion" << endl ;
        for (long int i = 0 ; i < ninterp ; ++i)
            IntegralePrems3_Newt( sp[i], si[i], so[i],
                                    int_M0, int_M1, int_M2,
                                    center_of_mass_positions[i], center_of_mass_impulsions[i],
                                    energies[i]  ) ;
    }
    else if (integrator_type == 10)
    {
        cout << endl << "Computing 2-body Newtonian with quadrupole term integrals of motion" << endl ;
        for (long int i = 0 ; i < ninterp ; ++i)
        {
            IntegralePrems2_NewtQuad( sp[i], si[i],
                                    int_M0, int_M1,  int_quadrupole[0],
                                    center_of_mass_positions[i], center_of_mass_impulsions[i],
                                    energies[i]  ) ;
        }
    }
    return ;
}


void Fittriple::Test_delays()
{
    long int i = 0L;
    value_type ** statevp;   // Contains the state vertors of the pulsar sp[number of times][6]
    value_type ** statevi;   // idem for inner companion
    value_type ** statevo;   // idem for outer companion
    value_type ** nss;
    value_type ** obs_position;
    value_type spinaxis[3];
    value_type ae = parameters.aB ;
    value_type ai = parameters.ap ;
    value_type Pi = parameters.Pi ;
    value_type Pe = parameters.Po ;
    value_type ra = parameters.RA ;
    value_type rap = 10.L / 365.25 / 1000. / 3600. * pi/180. ; //parameters.RA1 ; overwriting RA1 to make sure it's a big value (10mas/yr converted in rad/day)
    value_type dec = parameters.DEC ;
    value_type decp = 2.L / 365.25 / 1000. / 3600. * pi/180. ; //parameters.DEC1 ;
    value_type d = parameters.distance;
    value_type dprime = 500.L; // parameters.distance1; overwriting distance1 to make sure it's a significant value
    value_type mi = parameters.Mi;
    value_type me = parameters.Mo;
    value_type spinfreq = parameters.f;
    value_type muip = 100.L ;
    value_type muep = 1000.L ;
    value_type t0 = parameters.treference ;

    // Analytical references
    value_type * geometric ;
    value_type * einstein ;
    value_type * shapiro ;
    value_type * aberration ;
    
    // Numerical calculations
    value_type * geo_num ;
    value_type * ein_num ;
    value_type * shap_num ;
    value_type * ab_num ;
    value_type * total_num ;
    bool nanflag= false;

    printf("\n ********************* Diagnostic of validity and accuracy of numerical calculations **************** \n\n");

        statevp = (value_type**) malloc(sizeof(value_type*) * ninterp );
        statevi = (value_type**) malloc( sizeof(value_type*) * ninterp );
        statevo = (value_type**) malloc( sizeof(value_type*) * ninterp );
        nss =  (value_type**) malloc( sizeof(value_type*) * ninterp );
        obs_position = (value_type**) malloc( sizeof(value_type*) * ninterp );

        geometric = (value_type*) malloc( sizeof(value_type) * ninterp );
        einstein = (value_type*) malloc( sizeof(value_type) * ninterp );
        shapiro = (value_type*) malloc( sizeof(value_type) * ninterp );
        aberration = (value_type*) malloc( sizeof(value_type) * ninterp );
        geo_num = (value_type*) malloc( sizeof(value_type) * ninterp );
        ein_num = (value_type*) malloc( sizeof(value_type) * ninterp );
        shap_num = (value_type*) malloc( sizeof(value_type) * ninterp );
        ab_num = (value_type*) malloc( sizeof(value_type) * ninterp );
        total_num = (value_type*) malloc( sizeof(value_type) * ninterp );
//         obs_ssb = (value_type**) realloc(obs_ssb, sizeof(value_type*) * ninterp );
//         roemer_ss = (value_type*) realloc(roemer_ss, ninterp * sizeof(value_type) ) ;

        for (i = 0; i < ninterp ; ++i)
        {
            statevp[i] = (value_type*) malloc(sizeof(value_type) * 6 );
            statevi[i] = (value_type*) malloc(sizeof(value_type) * 6 );
            statevo[i] = (value_type*) malloc(sizeof(value_type) * 6 );
            nss[i] = (value_type*) malloc(sizeof(value_type) * 3 );
            obs_position[i] = (value_type*) malloc(sizeof(value_type) * 3 );
//             obs_ssb[i] = (value_type*) malloc(sizeof(value_type) * 3 );
        }

    Generate_einstein_geometric_state_vectors_diagnostic(tinterp, ninterp,
                                                        Pe, ae, Pi, ai,
                                                        ra, rap, dec, decp,
                                                        statevp, statevi, statevo,
                                                        nss, obs_position, spinaxis,
                                                        muip, muep  ) ;

   Compute_geometric_diagnostic(tinterp , t0 , geometric, ninterp,
                                     Pe, Pi, ae, ai,
                                     ra, rap, dec, decp,
                                     d,  dprime) ;

   Compute_einstein_diagnostic_for_geometric(tinterp, t0 , einstein, ninterp,
                                     Pe, Pi, ae, ai, mi, me,
                                     ra, rap, dec, decp,
                                     muip, muep  );

   Compute_shapiro_diagnostic(tinterp, shapiro, ninterp,
                                     mi, me, Pe, Pi, ae, ai,
                                     ra, rap, dec, decp,
                                     muip, muep);

   Compute_aberration_diagnostic(tinterp, aberration, ninterp,
                                     spinfreq, Pe, Pi, ae, ai,
                                     ra, rap, dec, decp) ;


           Delays_Brut_fullgeometric( tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 statevp, statevi, statevo, nss, obs_position, NULL,
                 mi, me, spinfreq, d, dprime,
                 true, false, false, false,
                 geo_num, nanflag, spinaxis
                ) ;

           Delays_Brut_fullgeometric( tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 statevp, statevi, statevo, nss, obs_position, NULL,
                 mi, me, spinfreq, d, dprime,
                 false, true, false, false,
                 ein_num, nanflag, spinaxis
                ) ;
           Delays_Brut_fullgeometric( tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 statevp, statevi, statevo, nss, obs_position, NULL,
                 mi, me, spinfreq, d, dprime,
                 false, false, true, false,
                 shap_num, nanflag, spinaxis
                ) ;
           Delays_Brut_fullgeometric( tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 statevp, statevi, statevo, nss, obs_position, NULL,
                 mi, me, spinfreq, d, dprime,
                 false, false, false, true,
                 ab_num, nanflag, spinaxis
                ) ;

       value_type dgeo = zero;
       value_type dein = zero;
       value_type dshap = zero;
       value_type dab = zero;

       FILE * myfile ;
       myfile = fopen("delays_diagnostics-results.dat", "w") ;
       for(i=0 ; i < ninterp ; ++i)
       {
            if (abs(geometric[i] - geo_num[i]) > dgeo ) dgeo = abs(geometric[i] - geo_num[i]);
            if (abs(einstein[i] - ein_num[i]) > dein) dein = abs(einstein[i] - ein_num[i] );
            if (abs(shapiro[i] - shap_num[i]) > dshap) dshap = abs(shapiro[i] - shap_num[i]);
            if (abs(aberration[i] - ab_num[i]) > dab ) dab = abs(aberration[i] - ab_num[i]);

            fprintf(myfile, "%.19Le    %.19Le    %.19Le    %.19Le    %.19Le\n", geometric[i] - geo_num[i],
                                                                                einstein[i] - ein_num[i],
                                                                                shapiro[i] - shap_num[i],
                                                                                aberration[i] - ab_num[i],
           geometric[i] - geo_num[i] + einstein[i] - ein_num[i] + shapiro[i] - shap_num[i] + aberration[i] - ab_num[i] );

           total_num[i] = geo_num[i] + ein_num[i] + shap_num[i] + ab_num[i] ;
//            if ( i < 10 )
//                printf(" diag %i : g %.19Le    %.19Le  |e  %.19Le    %.19Le  |s  %.19Le    %.19Le  |a  %.19Le    %.19Le\n",
//                       geometric[i], geo_num[i], einstein[i], ein_num[i], shapiro[i], shap_num[i], aberration[i], ab_num[i]);
       }
       fclose(myfile);

       dgeo *= daysec * pow(dix, 6);
       dein *= daysec * pow(dix, 6);
       dshap *= daysec * pow(dix, 6);
       dab *= daysec * pow(dix, 6);
       printf( "* Max differences (microsec) :\n  Geometric %.19Le\n  Einstein %.19Le\n  Shapiro %.19Le\n  Aberration %.19Le\n",
               dgeo, dein, dshap, dab);

       //***** Checking interpolation accuracy

       long int margin = 0L;
       value_type dinterpmax = zero;
       value_type dinterp ;
       Spline delays_interpolation( tinterp, total_num, ninterp );
       long int ninterpfin = 10 * ninterp;
       value_type * tinterpfin;
       value_type dtinterpfin;
       value_type delais;
       value_type interp;

       tinterpfin = (value_type*) malloc(sizeof(value_type) * ninterpfin );
       geometric = (value_type*) realloc(geometric, sizeof(value_type) * ninterpfin );
       einstein = (value_type*) realloc(einstein, sizeof(value_type) * ninterpfin );
       shapiro = (value_type*) realloc(shapiro, sizeof(value_type) * ninterpfin );
       aberration = (value_type*) realloc(aberration, sizeof(value_type) * ninterpfin );

       dtinterpfin = (tinterpfin[ninterpfin -1] - tinterpfin[0] ) / static_cast<value_type>(ninterpfin) ;
       for (i = 0; i < ninterpfin ; i++)
           tinterpfin[i] = tinterp[0] + i * dtinterpfin;

        Compute_geometric_diagnostic(tinterpfin , t0 , geometric, ninterpfin,
                                     Pe, Pi, ae, ai,
                                     ra, rap, dec, decp,
                                     d,  dprime) ;

        Compute_einstein_diagnostic_for_geometric(tinterpfin, t0 , einstein, ninterpfin,
                                            Pe, Pi, ae, ai, mi, me,
                                            ra, rap, dec, decp,
                                            muip, muep  );

        Compute_shapiro_diagnostic(tinterpfin, shapiro, ninterpfin,
                                            mi, me, Pe, Pi, ae, ai,
                                            ra, rap, dec, decp,
                                            muip, muep);

        Compute_aberration_diagnostic(tinterpfin, aberration, ninterpfin,
                                            spinfreq, Pe, Pi, ae, ai,
                                            ra, rap, dec, decp) ;

        margin = static_cast<long int>(floor( 0.012 / dtinterpfin ) ) + 2L; // Take into account the fact tinterp are sats which vary mostly from the roemer delay of the earth (0.012 ~ 18 ligth min in days)
        for (i = 0; i < ninterpfin ; i++)
        {
            delais = geometric[i] + einstein[i] + shapiro[i] + aberration[i] ;
            interp = delays_interpolation(tinterpfin[i], max(0, floor(i/10.) - margin) , min(ninterp, floor(i/10.) + margin) )  ;
            dinterp = abs( delais - interp) ;
            if (dinterp > dinterpmax )
                dinterpmax = dinterp;
            if (dinterp * daysec * pow(dix,6) > 0.01 and ( i > 10 * parameters.interp_margin or  i < (ninterpfin - 10 * parameters.interp_margin) ))
                printf("Warning : dinterp > 0.1 microsec for time %li. dinterp = %.19Le  %.19Le %.19Le \n", i, dinterp * daysec * pow(dix,6), delais, interp);
//             if (dinterp * daysec * pow(dix,6) > 0.1 and i < (ninterpfin - 10 * parameters.interp_margin) )
//                 printf("Warning : dinterp > 0.1 microsec for time %i. dinterp = %.19Le\n", i);
        }
        printf( "\n* Max difference of total delay after interpolation (microsec) : %.19Le \n", dinterpmax * daysec * pow(dix, 6) );


       free(tinterpfin);
       for (i = 0; i < ninterp ; ++i){
            free(statevp[i]);
            free(statevi[i]);
            free(statevo[i]);
            free(nss[i]);
            free(obs_position[i] ) ;
        }
        free(statevp);
        free(statevi);
        free(statevo);
        free(nss);
        free(obs_position);
        free(geometric);
        free(einstein);
        free(shapiro);
        free(aberration);
        free(geo_num);
        free(ein_num);
        free(shap_num);
        free(ab_num);

        printf("\n ************************************************************************************ \n\n");
}



void Fittriple::Get_global_interpolation_values(long int ntimes, long double * times, double * einstein_interp, double * shapiro_aberration_interp)
{
    //long int i = 0;

    diag_einstein_interp = einstein_interp ; //(double *) malloc( sizeof(double) * ntimes);
    diag_shapiro_aberration_interp = shapiro_aberration_interp; //(double *) malloc( sizeof(double) * ntimes);
    diag_interpolated_delays_n = ntimes;
    diag_interpolated_delays_times = times;//(value_type *) malloc( sizeof(value_type) * ntimes);

    //for (i = 0 ; i < ntimes ; i++) diag_interpolated_delays_times[i] = times[i];

    get_global_interpolated_delays_flag = true ;
    Compute_lnposterior();
    get_global_interpolated_delays_flag = false;

    //for

    //free(diag_einstein_interp);
    //free(diag_shapiro_aberration_interp);
    //free(diag_interpolated_delays_times);
}


void Fittriple::Get_local_interpolation_values(long int ntimes, long double * relativetimes, double * geometrical_inter)
{
}
