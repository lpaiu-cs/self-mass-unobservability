// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include "Constants.h"
#include "Utilities.h"
#include "Fittriple.h"
#include "Orbital_elements.h"
#include "tempo2.h"
#include "IO.h"




void Fittriple::Save_parfile(const char * filename)
{
    parameters.Save_parfile(filename);
    return ;
}


void Fittriple::Save_data(const char * filename)
{
   FILE *fp;

   fp = fopen(filename,"w");

   long int i = 0;

   fprintf(fp, "%s \n", "# toas (MJD)    turn number   uncertainty (microsecs)");
   for (i=0 ; i < ntoas -1 ; ++i){
        fprintf(fp, "%.19Le    %lli    %.3f\n", toas[i] + parameters.timeshift, turns[i], errors[i]);
    }
    fprintf(fp, "%.19Le    %lli    %.3f", toas[i] + parameters.timeshift, turns[i], errors[i]);

    fclose(fp);

    return;
}


void Fittriple::Save_turn_numbers(char * filename)
/*
 * Save the turn numbers used internally : " turns "
 */
{
   FILE *fp;

   fp = fopen(filename,"w");

   long int i = 0;

   fprintf(fp, "%s \n", "#turns");
   for (i=0 ; i < ntoas -1 ; ++i){
        fprintf(fp, "%lli\n",turns[i]);
    }
     fprintf(fp, "%lli",turns[i]);
    fclose(fp);

    return;
}


void Fittriple::Load_turn_numbers(char * filename)
{
        FILE *fp;
    const int sizechar = 50;
    int i;
    turntype tt ;
    char charun[1];
    char chardash[] = "#";
    char line[sizechar];
    char textin[sizechar]; // Should be the same size as line otherwise sscanf fails


    fp = fopen(filename,"r");
    if (fp == NULL)
    {
        printf("\n\n Opening of turn file %s failed ! \n\n", filename);
        return;
    }
     fgets(line, sizechar, fp);
     printf(" \n !!! WARNING Load_turn_numbers : DATA MUST BEGIN ON SECOND LINE !!!! \n\n");
    i = 0;
    while ( fgets(line, sizechar, fp) != NULL && i < ntoas)
    {
        sscanf( line, "%lli", &turns[i]);//
        i++;
    }
    fclose(fp);
    return;
}

void Fittriple::Save_output_timing_data(const char * filename)
/*
 * Save the output data with toas : " toas toes residuals "
 */
{
   FILE *fp;

   fp = fopen(filename,"w");

   long int i = 0;

   fprintf(fp, "%s \n", "# toas (MJD)    toes (MJD)    residuals (microsecs)");
   for (i=0 ; i < ntoas -1 ; ++i){
        fprintf(fp, "%.19Le    %.19Le    %.19Le\n", toas[i] + parameters.timeshift , toes[i] + parameters.timeshift, residuals[i]);
    }
    fprintf(fp, "%.19Le    %.19Le    %.19Le", toas[i] + parameters.timeshift, toes[i] + parameters.timeshift, residuals[i]);

    fclose(fp);

    return;
}



void Fittriple::Save_interp_state_vectors(const char * filename)
/*
 * Save  "tinterp sp si so"
 */
{
   FILE *fp;

   fp = fopen(filename,"w");

   long int i = 0;

   fprintf(fp, "# tinterp(days, timeshift %.19Le )%s \n", parameters.timeshift, "    sp     si     so   with for each state vector (x    y(line of sight)    z) in meters and  (vx   vy   vz) in meters/second )");
   for (i=0 ; i < ninterp -1 ; ++i){
        fprintf(fp, "%.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le\n",
                tinterp[i],
                sp[i][0], sp[i][1], sp[i][2], sp[i][3], sp[i][4], sp[i][5],
                si[i][0], si[i][1], si[i][2], si[i][3], si[i][4], si[i][5],
                so[i][0], so[i][1], so[i][2], so[i][3], so[i][4], so[i][5]      );
    }
    fprintf(fp, "%.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le",
                tinterp[i],
                sp[i][0], sp[i][1], sp[i][2], sp[i][3], sp[i][4], sp[i][5],
                si[i][0], si[i][1], si[i][2], si[i][3], si[i][4], si[i][5],
                so[i][0], so[i][1], so[i][2], so[i][3], so[i][4], so[i][5]      );

    fclose(fp);

    return;
}



void Fittriple::Save_masked_copy_of_tim_file(const char * masked_tim_file, const int * mask)
{
    value_type * SATs;
    double * errs;
    
    SATs = (value_type*) malloc(ntoas * sizeof(value_type));
    errs = (double*) malloc(ntoas * sizeof(double));
    
    for (int i = 0; i < ntoas; i++)
    {
        SATs[i] = tempo2_psr[0].obsn[i].sat;
        errs[i] = tempo2_psr[0].obsn[i].toaErr;
    };
    
    // tempo2_timfile should point to the sorted file already. 
    Write_new_tim_file(tempo2_timfile, masked_tim_file, SATs, errs, mask, 1);
    
    free(SATs);
    free(errs);
}


void Fittriple::Print_initial_state_vector()
/*
 * Convert x0 back to physical units (m and m/s in principal) and print it
 */
{
  value_type x0normal[18] ;

  x0normal[0] = x0[0] * length ;
  x0normal[1] = x0[1] * length ;
  x0normal[2] = x0[2] * length ;
  x0normal[9] = x0[9] * length / timescale ;
  x0normal[10] = x0[10] * length / timescale ;
  x0normal[11] = x0[11] * length / timescale ;

  x0normal[3] = x0[3] * length ;
  x0normal[4] = x0[4] * length ;
  x0normal[5] = x0[5] * length ;
  x0normal[12] = x0[12] * length / timescale ;
  x0normal[13] = x0[13] * length / timescale ;
  x0normal[14] = x0[14] * length / timescale ;

  x0normal[6] = x0[6] * length ;
  x0normal[7] = x0[7] * length ;
  x0normal[8] = x0[8] * length ;
  x0normal[15] = x0[15] * length / timescale ;
  x0normal[16] = x0[16] * length / timescale ;
  x0normal[17] = x0[17] * length / timescale ;

    cout << "*********** Initial state vector at treference in m and m/s***********" << endl ;
    Print_table(x0normal, 18);

    return;
}

void Fittriple::Print()
{
    cout << endl << "********** Par file ***************" << endl ;
    parameters.Print();
//     printf("%s %.19Le \n", "Time of reference, treference + timeshift :", treference + timeshift ) ;
//
//     printf("%s %i \n", "roemer", int(roemer) ) ;
//     printf("%s %i \n", "einstein", int(einstein) ) ;
//     printf("%s %i \n", "shapiro", int(shapiro) ) ;
//     printf("%s %i \n", "aberration", int(aberration) ) ;
//
//     printf("%s %i \n", "integrator", integrator_type ) ;
//     printf("%s %i \n", "interpsteps_per_period_i", interpsteps_per_period_i ) ;
//
//     printf( "%s %.19Le  %.19Le \n", "spinfreq", parameters[0], parameter_shift_scales[0] ) ;
//     printf( "%s %.19Le %.19Le \n", "spinfreq1", parameters[1], parameter_shift_scales[1] ) ;
//     printf( "%s %.19Le %.19Le \n", "eta_p", parameters[2], parameter_shift_scales[2] ) ;
//     printf( "%s %.19Le %.19Le \n", "apsini_i", parameters[3], parameter_shift_scales[3] ) ;
//     printf( "%s %.19Le %.19Le \n", "kappa_p", parameters[4], parameter_shift_scales[4] ) ;
//     printf( "%s %.19Le %.19Le \n", "apcosi_i", parameters[5], parameter_shift_scales[5] ) ;
//     printf( "%s %.19Le %.19Le \n", "tasc_p", parameters[6], parameter_shift_scales[6] ) ;
//     printf( "%s %.19Le %.19Le \n", "period_i", parameters[7], parameter_shift_scales[7] ) ;
//     printf( "%s %.19Le %.19Le \n", "eta_b", parameters[8], parameter_shift_scales[8] ) ;
//     printf( "%s %.19Le %.19Le \n", "absini_o", parameters[9], parameter_shift_scales[9] ) ;
//     printf( "%s %.19Le %.19Le \n", "kappa_b", parameters[10], parameter_shift_scales[10] ) ;
//     printf( "%s %.19Le %.19Le \n", "abcosi_o", parameters[11], parameter_shift_scales[11] ) ;
//     printf( "%s %.19Le %.19Le \n", "tasc_b", parameters[12], parameter_shift_scales[12] ) ;
//     printf( "%s %.19Le %.19Le \n", "period_o", parameters[13], parameter_shift_scales[13] ) ;
//     printf( "%s %.19Le %.19Le \n", "masspar_p", parameters[14], parameter_shift_scales[14] ) ;
//     printf( "%s %.19Le %.19Le \n", "masspar_i", parameters[15], parameter_shift_scales[15] ) ;
//     printf( "%s %.19Le %.19Le \n", "mass_o", parameters[16], parameter_shift_scales[16] ) ;
//     printf( "%s %.19Le %.19Le \n", "deltaoman", parameters[17], parameter_shift_scales[17] ) ;
//
//
//     printf("%s %li \n", "interp_margin", interp_margin ) ;
//     printf("%s %.19Le \n", "Timeshift applied to times internally, timeshift : ", timeshift );
//     printf("%s %.19Le \n", "Tolerance level of integrator, tolint : ", tolint ) ;


    cout << endl << "********** Data ***************" << endl ;
    printf("%s %li \n", "Number of toas, ntoas :", ntoas ) ;
    printf("%s %.19Le \n", "Range of toas, toas[ntoas-1] - toas[0] :", toas[ntoas-1] - toas[0] ) ;
    printf("%s %.19Le \n", "First toa, toas[0] + timeshift : ", toas[0] + parameters.timeshift ) ;
    printf("%s %.19Le \n", "Last toa, toas[ntoas-1] + timeshift : ", toas[ntoas-1] + parameters.timeshift ) ;

    cout << endl << "********** Interpolation ***************" << endl ;
    printf("%s %.19Le \n", "dt_interp", dt_interp ) ;
    printf("%s %li \n", "ninterp", ninterp ) ;
    printf("%s %li \n", "treference_in_interp", treference_in_interp ) ;
    printf("%s %.19Le \n", "First interpolation time :", tinterp[0] ) ;
    printf("%s %.19Le \n", "Last interpolation time : ", tinterp[ninterp-1] ) ;

    cout << endl << "********** Others ***************" << endl ;
    printf("%s %.19Le \n", "spinfreq", spinfreq ) ;
    printf("%s %.19Le \n", "spinfreq1 ", spinfreq1  ) ;
    printf("%s %.19Le \n", "delaymax", delaymax ) ;

    cout << endl << "********** Integrateur ***************" << endl ;
    printf("%s %.19Le \n", "Tolerance level of integrator, tolint : ", tolint ) ;
    printf("%s %.19Le \n", "Mass scale, mass : ", mass ) ;
    printf("%s %.19Le \n", "Length scale, length : ", length ) ;
    printf("%s %.19Le \n", "Time scale (sec), timescale : ", timescale ) ;
    printf("%s %.19Le \n", "Time scale (days), timedays : ", timedays ) ;
    printf("%s %.19Le \n", "Integration starting time, t0 : ", t0 ) ;
    printf("%s %.19Le \n", "Starting step, dt0 ", dt0 ) ;
    printf("%s %li \n", "ntneg ", ntneg ) ;
    printf("%s %.19Le \n", "int_M0 ", int_M0 ) ;
    printf("%s %.19Le \n", "int_M1 ", int_M1 ) ;
    printf("%s %.19Le \n", "int_M2 ", int_M2 ) ;
    printf("%s %li \n", "Number of integration stops, n ", n ) ;
    printf("%s %i \n", "roemer_second_order ", roemer_second_order ) ;
    cout << endl <<  " Initial state vector (at time t0) : " << endl ;
    Print_table(x0);
    cout << endl ;
    cout << endl << "********** Tempo 2 ***************" << endl ;
    Print_tempo();

    return ;
}


void Fittriple::Print_tempo()
{
        int i = 0;
        double maxshklovskii = 0.;
        double test = 1.;

        cout << endl << "*** Infos from tempo2 " << endl;
        printf(" Tempo vector Earth-> pulsar : %.10f %.10f %.10f : \n", tempo2_psr[0].posPulsar[0],  tempo2_psr[0].posPulsar[1],  tempo2_psr[0].posPulsar[2]);
        printf(" Tempo velocity pulsar : %.10f %.10f %.10f : \n", tempo2_psr[0].velPulsar[0],  tempo2_psr[0].velPulsar[1],  tempo2_psr[0].velPulsar[2]);
        printf(" Tempo nobs : %i \n", tempo2_psr[0].nobs);
        for (i = 0 ; i < ntoas ; i++)
        {
            if ( abs(tempo2_psr[0].obsn[i].shklovskii) > maxshklovskii ) maxshklovskii = abs(tempo2_psr[0].obsn[i].shklovskii);
        }


        printf(" Tempo : max shklovskii : %.10f : \n", maxshklovskii);
        printf(" Tempo : ecliptic coordinates : %i \n" , tempo2_psr[0].eclCoord==1) ;
        printf(" Tempo : POSEPOCH : %.19Le\n", tempo2_psr[0].param[param_posepoch].val[0]);
       /* printf(" Tempo : position of Earth wrt SSB at POSEPOCH (added to tempo by G. Voisin) : %.10f %.10f %.10f : \n", tempo2_psr[0].earth_ssb_at_posepoch[0],
                                                                                                                      tempo2_psr[0].earth_ssb_at_posepoch[1],
                                                                                                                      tempo2_psr[0].earth_ssb_at_posepoch[2]);
        printf(" Tempo : velocity of Earth wrt SSB at POSEPOCH (added to tempo by G. Voisin) : %.10f %.10f %.10f : \n", tempo2_psr[0].earth_ssb_at_posepoch[3],
                                                                                                                      tempo2_psr[0].earth_ssb_at_posepoch[4],
                                                                                                                      tempo2_psr[0].earth_ssb_at_posepoch[5]);
        */printf(" Tempo : psr[p].param[param_px].val[0] : %.10Lf \n", tempo2_psr[0].param[param_px].val[0] );
        printf(" Tempo : psr[p].param[param_pmrv].paramSet[0] : %i \n", tempo2_psr[0].param[param_pmrv].paramSet[0] );
        printf(" Tempo : psr[p].t2cMethod - T2C_IAU2000B : %i \n", tempo2_psr[0].t2cMethod - T2C_IAU2000B );
        for (i = 0; i< ntoas; i++)
            test *= tempo2_psr[0].obsn[i].delayCorr;
        printf(" Tempo : Any psr[0].obsn[i].delayCorr = 0 ?  (not good !!) (0 if yes) : %e \n", test);
        test = 1.;
        for (i = 0; i< ntoas; i++)
            test *= strcmp(tempo2_psr[0].obsn[i].telID,"STL_FBAT");
        printf(" Tempo : strcmp(psr[p].obsn[i].telID,'STL_FBAT') = 0 ?  (0 pas bien !) (0 if yes): %e \n", test );
        printf(" Tempo : psr[0].obsn[i].telID : %s\n", tempo2_psr[0].obsn[i].telID);
        printf(" Tempo : psr[p].nTelDX (0 bien) : %i\n", tempo2_psr[0].nTelDX );
        printf(" Tempo : psr[p].nTelDY (0 bien) : %i\n", tempo2_psr[0].nTelDZ );
        printf(" Tempo : psr[p].nTelDZ (0 bien) : %i\n", tempo2_psr[0].nTelDZ );
        printf(" Tempo : JPL ephemeris file : %s \n" , tempo2_psr[0].JPL_EPHEMERIS );

        printf(" Tempo : DM : %.10Le \n", tempo2_psr[0].param[param_dm].val[0] ) ;
        printf(" Tempo : DM1 : %.10Le \n", tempo2_psr[0].param[param_dm].val[1]) ;

        return;
};



void Fittriple::Print_tests()
/*
 * Put here what you want to print for testing purposes.
 */
{
    cout << endl << endl << "************ Print_tests ***************** " << endl << endl ;


//     Print_table(tneg);
//     cout << "------------" << endl;
//     Print_table(tpos);
//     cout << endl <<"  tsssssssss" << endl ;
//     Print_table(ts);
//     cout << endl << " toa  toe" << endl ;
//     Print_table(toas, toes,ntoas);


    cout << endl << " Residuals " << endl ;
    Print_table(residuals, ntoas);


    cout << endl << endl << "************ End Print_tests ************* " << endl << endl ;

}

long int Fittriple::Get_number_of_toas()
{
    return ntoas;
}

double Fittriple::Get_toa_double(long int toa_number)
{

//     number_of_toas = ntoas;
//     times = new double[ntoas];
//     for (long int i = 0L; i< ntoas ; ++i) times[i] = static_cast<double>( toas[i] );
    return static_cast<double>( toas[toa_number] ) ;
}

double Fittriple::Get_residual_double(long int toa_number)
{
// Return residuals computed in Fittriple-compute.cpp. Should be in microseconds.
    return static_cast<double>( residuals[toa_number] );
}


double Fittriple::Get_error_double(long int toa_number)
{
    return static_cast<double>( errors[toa_number] );
}

double Fittriple::Get_parameter_double(int param_number )
{
    return static_cast<double>(parameters[param_number]);
}

double Fittriple::Get_reference_parameter_double(int param_number )
{
    return static_cast<double>(parameters.parameters_ini[param_number]);
}

double Fittriple::Get_parameter_scale_double(int param_number )
{
    return static_cast<double>( parameters.Parameter_scale(param_number) );
}


double Fittriple::Get_mass_double(int bodynumber)
{
    if (bodynumber == 0)
        return static_cast<double>( parameters.Mp );
    else if (bodynumber == 1)
        return static_cast<double>(parameters.Mi);
    else if (bodynumber == 2)
        return static_cast<double>(parameters.Mo );
    else if (bodynumber > 2 and bodynumber < parameters.nbodies_plus_extra)
        return static_cast<double>(parameters.M_extra[bodynumber - 3]);
    else
        cout << "Warning : invalid bodynumber in Get_mass(). bodynumber = " << bodynumber << endl;
        return static_cast<double>( -1. );
}



double Fittriple::Get_delays(int delaynb, long int toa_number)
 // Return delay number delaynb at toa number toa_number in double precision in seconds.
// (Currently disabled )delay < 0 : return some tempo delays 
// delaynb = 0 : return the total delay
// delaynb = 1 : return the geometrical delay
// delaynb = 2 : return the sum of the enabled non-geometrical delays (eintein, shapiro, aberration)
// -- Decomposition of the geometrical delay --
// delaynb = 10 : roemer delay
// delaynb = 11 : schklovski delay
// delaynb = 12 : scklovski correction delay
// delaynb = 13 : kopeikin orbital delay
// delaynb = 14 : kopeikin annual-orbital delay
// delaynb = 15 : kopeikin parallax-orbital delay

    {
        if (delaynb > 0 and getdelaysflag <= 0)
        {
            getdelaysflag = 1;
            geomdelay = (double*) malloc( sizeof(double)*ntoas ) ;
            nogeomdelay = (double*) malloc( sizeof(double)*ntoas ) ;
            einsteindelay = (double*) malloc( sizeof(double)*ntoas ) ;
            geom_delay_details = (double**) malloc(sizeof(double*) * ntoas);
            for (int i = 0 ; i < ntoas ; i++) geom_delay_details[i] = (double*) malloc(sizeof(double) * 6);
            Compute_lnposterior();
        }
        if (delaynb ==0) // return total delays
            return (toas[toa_number] - toes[toa_number] ) * 86400. ;
// THe following block requires tempo2 modifications 
//         else if (delaynb < 0) // Return tempo2 model delays
//         {
//             if (parameters.shapiro == false)
//                 tempo2_psr[0].grrrr_T2shapiro = 0;
//             else
//                 tempo2_psr[0].grrrr_T2shapiro = 1;
//             tempo2_psr[0].grrrr_T2kopeikin = 1;
//             formResiduals(tempo2_psr,1,1);
//             printf("T2 omdot %e \n", tempo2_psr[0].grrrr_T2omdot);
//             if (delaynb == -2 )
//                 return tempo2_psr[0].obsn[toa_number].grrrr_T2dre;
//             else
//                 return  tempo2_psr[0].obsn[toa_number].torb;
//         }
        else if (delaynb ==1)
        {
            return geomdelay[toa_number] * 86400.;
        }
        else if (delaynb ==2)
        {
            return nogeomdelay[toa_number] * 86400.;
        }
        else if (delaynb == 3)
        {
            return einsteindelay[toa_number] * 86400.;
        }
        else if (delaynb == 10) 
            return geom_delay_details[toa_number][0] * 86400.;
        else if (delaynb == 11) 
            return geom_delay_details[toa_number][1] * 86400.;
        else if (delaynb == 12) 
            return geom_delay_details[toa_number][2] * 86400.;
        else if (delaynb == 13) 
            return geom_delay_details[toa_number][3] * 86400.;
        else if (delaynb == 14) 
            return geom_delay_details[toa_number][4] * 86400.;
        else if (delaynb == 15) 
            return geom_delay_details[toa_number][5] * 86400.;
        else
            return 0.;
    };


void Fittriple::Get_fake_bats_and_delays_interp(long int nfakeBats, value_type * fake_BATs, value_type * fake_delay_geom, value_type * fake_delay_ein,
                                                                                   value_type * fake_delay_shap, value_type * fake_delay_aber)
// Use Compute_fake_BATS_and_delays_from_parameters to create fake bats from every interpolation time from 1 to nfakeBats + 1. nfakeBats should be smaller than ninterp - 1.
// Each delay is computed according to the current parameters.
// At actual toas, this should give the correct values, execpt for the geometric delay for which all coupling terms between the earth orbit and the pulsar orbit are neglected.
{

    long int ntisaround =   2* (floor(parameters.interp_margin * parameters.interpsteps_per_period_i) +  floor(delaymax / dt_interp ) ) ; // Number of interpolation steps to take around each toa to compute geometric delays. must be even.
    long int re = (ninterp- 2* ntisaround) % nfakeBats;
    long int step = (ninterp - 2 *ntisaround) / nfakeBats;
    if (step == 0 ) printf("\n Warning : nfakeBats=%li must be lower than (ninterp - 2*ntisaround)=%li \n", nfakeBats, ninterp- 2* ntisaround);
    for (long int i=0 ; i<nfakeBats ; i++)
    {
        fake_BATs[i] = tinterp[ntisaround + re + step * i];
    }
    Compute_fake_BATS_and_delays_from_parameters( nfakeBats,fake_BATs,  fake_delay_geom,  fake_delay_ein,
                                                                                      fake_delay_shap, fake_delay_aber) ;

    return;
}



void Fittriple::Get_analytical_pulsar_statevector_double(long int tinterp_nb, double * statevector)
// Return the analytical, newtonian, state vector computed from the orbital elements at time tinterp[tinterp_nb] with orbel2statevect
{
    value_type sv[6];
// double svdbl[6];
    if (parameters.integrator_type == 1)
        orbel2statevect_1pn(tinterp[tinterp_nb], parameters.ei, parameters.ap, parameters.omp, parameters.anglii, parameters.tperii, parameters.Pi, parameters.delta_oman,
                        parameters.Mp, parameters.Mi,
                            sv,
                            1, 1) ;
    else
        orbel2statevect(tinterp[tinterp_nb], parameters.ei, parameters.ap, parameters.omp, parameters.anglii, parameters.tperii, parameters.Pi, parameters.delta_oman,
                            sv,
                            1, 1) ;
    valarray<value_type> rtest(3), vtest(3);
    rtest[0] = sv[0];rtest[1] = sv[1]; rtest[2] = sv[2];
    vtest[0] = sv[3];vtest[1] = sv[4]; vtest[2] = sv[5];
    Rotate_PSB_to_SSB(rtest, parameters.RA, parameters.DEC) ;
    Rotate_PSB_to_SSB(vtest, parameters.RA, parameters.DEC) ;
    statevector[0] = static_cast<double>(rtest[0]); statevector[1] = static_cast<double>(rtest[1]); statevector[2] = static_cast<double>(rtest[2]);
    statevector[3] = static_cast<double>(vtest[0]); statevector[4] = static_cast<double>(vtest[1]); statevector[5] = static_cast<double>(vtest[2]);
    return;// svdbl[component];
}


void Fittriple::Get_analytical_inner_statevector_double(long int tinterp_nb, double * statevector)
// Return the analytical, newtonian, state vector computed from the orbital elements at time tinterp[tinterp_nb] with orbel2statevect
{
    value_type sv[6];
// double svdbl[6];
    value_type ai = parameters.ap * parameters.Mp / parameters.Mi;
    if (parameters.integrator_type == 1)
        orbel2statevect_1pn(tinterp[tinterp_nb], parameters.ei, ai, parameters.omp + pi, parameters.anglii, parameters.tperii, parameters.Pi, parameters.delta_oman,
                            parameters.Mi, parameters.Mp,
                            sv,
                            1, 1) ;
    else
        orbel2statevect(tinterp[tinterp_nb], parameters.ei, ai, parameters.omp + pi, parameters.anglii, parameters.tperii, parameters.Pi, parameters.delta_oman,
                            sv,
                            1, 1) ;
    valarray<value_type> rtest(3), vtest(3);
    rtest[0] = sv[0];rtest[1] = sv[1]; rtest[2] = sv[2];
    vtest[0] = sv[3];vtest[1] = sv[4]; vtest[2] = sv[5];
    Rotate_PSB_to_SSB(rtest, parameters.RA, parameters.DEC) ;
    Rotate_PSB_to_SSB(vtest, parameters.RA, parameters.DEC) ;
    statevector[0] = static_cast<double>(rtest[0]); statevector[1] = static_cast<double>(rtest[1]); statevector[2] = static_cast<double>(rtest[2]);
    statevector[3] = static_cast<double>(vtest[0]); statevector[4] = static_cast<double>(vtest[1]); statevector[5] = static_cast<double>(vtest[2]);
    return;// svdbl[component];
}



void Fittriple::Compute_initial_state_vectors_double(double ** svs)
// Run parameters.Compute_state_vectors() and return svs such that svs[0] = {rp, vp} svs[1] = {ri, vi} svs[2] = {ro, vo} svs[i>2] = {rextra_i, vextra_i}
{
    int i,j ;
    valarray<value_type> rp(3);
    valarray<value_type> rpt(3);
    valarray<value_type> ri(3);
    valarray<value_type> rit(3);
    valarray<value_type> ro(3);
    valarray<value_type> rot(3);
    valarray<valarray<value_type>> r_extra(parameters.nextra);
    valarray<valarray<value_type>> v_extra(parameters.nextra);
    for (i = 0; i < parameters.nextra; i++)
    {
      r_extra[i] = valarray<value_type>(3);
      v_extra[i] = valarray<value_type>(3);
    }

    parameters.Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
    
// Copy to svs
    for (j = 0 ; j < 3; j++) 
    {
        svs[0][j] = static_cast<double>(rp[j]); 
        svs[0][j+3] = static_cast<double>(rpt[j]); 
        svs[1][j] = static_cast<double>(ri[j]); 
        svs[1][j+3] = static_cast<double>(rit[j]);
        svs[2][j] = static_cast<double>(ro[j]); 
        svs[2][j+3] = static_cast<double>(rot[j]);
    }
    for (i = 3 ; i < parameters.nbodies_plus_extra ; i++) 
    {
        for (j = 0 ; j < 3; j++) 
        {
            svs[i][j] = static_cast<double>(r_extra[i-3][j]); 
            svs[i][j+3] = static_cast<double>(v_extra[i-3][j]);
        }
    }
}


/*
int * Fittriple::Get_absolute_parameter_map(){
        int * parammap;
        parammap = (int *)malloc(sizeof(int) * parameters.n_absparameters);
        for (int i = 0; i < parameters.n_absparameters ; i ++) parammap[i] = parameters.absolute_parameter_map[i];
        return parammap;
    };*/
