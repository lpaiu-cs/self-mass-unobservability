// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include "Constants.h"
#include <cmath>
#include "Utilities.h"
#include <valarray>
#include <cstring>
#include "Orbital_elements.h"
#include "Spline.h"
#include "Delay_brut.h"
#include "IO.h"
#include "Fittriple.h"
#include "tempo2.h"
#include "Diagnostics.h"


Fittriple::Fittriple(char * parfile, char * datafile) {

    tinterp = NULL;
    tempo2_psr = NULL;
    toas = NULL;
    turns = NULL;
    errors = NULL;
    weights=NULL;
    obs_ssb_BAT = NULL;
    roemer_ss_BAT = NULL;
    pulsardirection_BAT = NULL;
    sp = NULL;
    si = NULL;
    so = NULL;

    toa_in_interp = NULL;
    residuals = NULL;
    toes = NULL;

        // Diagnostics
    energies = valarray<value_type> (1) ;
    center_of_mass_positions = valarray<valarray<value_type> >(1);
    center_of_mass_impulsions = valarray<valarray<value_type> >(1);
    center_of_mass_positions[0] = valarray<value_type>(3);
    center_of_mass_impulsions[0] = valarray<value_type>(3);

    // Other initializations
    removed_toas = vector< long int >(0);
    disable_sats = false;
    motion_changed = true;
    remove_mean = true;
//     Fittriple();
    Reconstruct(parfile, datafile);
    return ;
}

Fittriple::Fittriple() {


    tinterp = NULL;
    tempo2_psr = NULL;
    toas = NULL;
    turns = NULL;
    errors = NULL;
    weights=NULL;
    obs_ssb_BAT = NULL;
    roemer_ss_BAT = NULL;
    pulsardirection_BAT = NULL;
    sp = NULL;
    si = NULL;
    so = NULL;

    toa_in_interp = NULL;
    residuals = NULL;
    toes = NULL;

    // Diagnostics
    energies = valarray<value_type> (1) ;
    center_of_mass_positions = valarray<valarray<value_type> >(1);
    center_of_mass_impulsions = valarray<valarray<value_type> >(1);
    center_of_mass_positions[0] = valarray<value_type>(3);
    center_of_mass_impulsions[0] = valarray<value_type>(3);

    // Other initializations
    removed_toas = vector< long int >(0);
    disable_sats = false ;
    motion_changed = true;
    remove_mean = true;

    return ;
}


Fittriple::~Fittriple()
{
    // things to free :
//     toas
// pulsardirection_BAT
// obs_ssb_BAT
// roemer_ss_BAT
// turns
// errors
// sat_to_bat
// sp
// si
// so
// tempo2_psr
// toes
// residuals
// toa_in_interp
// tinterp
// tsat_testmode
//
// sp[i]
// si[i]
// so[i]
// pulsardirection_BAT[i]
// obs_ssb_BAT[i]
// sat_to_bat[i]


     for (int i = 0; i < ninterp ; ++i){
            free(sp[i]);
            free(si[i]);
            free(so[i]);
        }
     free(sp);
     free(si);
     free(so);

    for (int i = 0; i < ntoas ; ++i){
            free(pulsardirection_BAT[i]);
            free(obs_ssb_BAT[i]);
        }
    free(toas);
    free(pulsardirection_BAT);
    free(obs_ssb_BAT);
    free(roemer_ss_BAT);
    free(turns);
    free(errors);
    free(weights);
    free(toes);
    free(residuals);
    free(toa_in_interp);
    free(tinterp);
    //free(tsat_testmode);
    destroyOne(tempo2_psr);
    free(tempo2_psr);

    if (getdelaysflag >0)
    {
        free(geomdelay);
        free(nogeomdelay);
        free(einsteindelay);
        for (int i = 0; i < ntoas; i++) free(geom_delay_details[i]);
        free(geom_delay_details);
    }
}

void Fittriple::Reconstruct_noCompute(char * parfile, char * datafile, const bool sorted_datafile)
/* Can be called after the default constructor Fittriple().
 * Should not be called after a previous Reconstruct or Fittriple(char * parfile, char * datafile) unless memory deaollation is implemented i.e. call a destructor before reconstructing.
 * Does not compute the model. Therefore turn numbers are not calculated !
 */
{
    long int i =0;

    // Diagnostics variables
    get_global_interpolated_delays_flag = false;
    diag_einstein_interp = NULL;
    diag_shapiro_aberration_interp = NULL;
    diag_interpolated_delays_n = 0;
    diag_interpolated_delays_times = NULL;
    getdelaysflag = 0 ;

    motion_changed = true;
    roemer_second_order = -1;
    mass = Msol ;

    // Tracker
    callnumber = 0;
    tracker = false;


    int_quadrupole=vector<value_type>(4,zero);

    parameters.Load_new_parfile(parfile);
    previous_RA = -1. ; // To avoid updating bats with tempo when ra or dec do not change between two calls. Used in "Initialize()"
    previous_DEC = -1. ;
    previous_RA1 = -1. ; // To avoid updating bats with tempo when ra or dec do not change between two calls. Used in "Initialize()"
    previous_DEC1 = -1. ;
    previous_distance1 =-1.;
    previous_DM = -1. ;  // To avoid updating bats with tempo when DM or DM1 do not change between two calls. Used in "Initialize()"
    previous_DM1 = -1. ;
    previous_DMX = vector<value_type>(parameters.DMX.size(), 0.L);
    if (previous_DMX == parameters.DMX and previous_DMX.size() > 0 ) previous_DMX[0] -=1.;
    previous_FD = vector<value_type>(parameters.FD.size(), 0.L);
    if (previous_FD == parameters.FD and previous_FD.size() > 0 ) previous_FD[0] -=1.;

    delaymaxold = 0.;
    forcerecomputeinterp = false ;

    tolint = parameters.tolint ;
    integrator_type = parameters.integrator_type ;



        MAX_PSR = 1;
        MAX_OBSN = 1000000;
        tempo2_npsr = 1;
        tempo2_psr = (pulsar *)realloc(tempo2_psr, sizeof(pulsar)*MAX_PSR);
        if (tempo2_psr == NULL) printf(" \n *** Error : memory allocation failed in Reconstruct for tempo2_psr !\n ");
        tempo2_noWarnings = 1;
        char parFile[tempo2_npsr][MAX_FILELEN];
        char timFile[tempo2_npsr][MAX_FILELEN];
        strcpy(tempo2_timfile, datafile);
        strcpy(tempo2_parfile, parameters.t2parfile);
        if (sorted_datafile == false)
        {
          Sortoutdatafile(datafile, MAX_FILELEN); // Sort out tim file by toa
          strcat(tempo2_timfile,"-sorted"); // use sorted datafile.
        }
        strcpy(parFile[0],tempo2_parfile);
        strcpy(timFile[0],tempo2_timfile);

        // Initialise Tempo2
        initialise(tempo2_psr,tempo2_noWarnings); // SINON SEGFAULT !!!
        readParfile(tempo2_psr,parFile,timFile,tempo2_npsr);
        readTimfile(tempo2_psr,timFile,tempo2_npsr);
        ntoas = static_cast<long int>(tempo2_psr[0].nobs);

        
        // Determine DMX ranges
        value_type DMr1 = zero;
        value_type DMr2 = zero;
        tempo2_psr[0].param[param_pmra].paramSet[0]=1  ;
        tempo2_psr[0].param[param_pmdec].paramSet[0]=1 ;
        tempo2_psr[0].param[param_dshk].paramSet[0] = 0;
        if (parameters.DMX.size() > 0)
        {
            parameters.DMXranges.resize(parameters.DMX.size()+1);
            DMr1 =  tempo2_psr[0].obsn[0].sat - un;
            DMr2 = tempo2_psr[0].obsn[ntoas-1].sat + un; // +/- 1 to ensure that the first and last toa are within range.
            for (i = 0; i < parameters.DMXranges.size() ; i++)
            {
                parameters.DMXranges[i] = DMr1 + i * (DMr2 - DMr1)/parameters.DMX.size();
            }
        }
         tempo2_psr[0].ndmx = parameters.DMX.size();
         for  (i = 0; i < parameters.DMX.size() ; i++)
         {
             tempo2_psr[0].param[param_dmxr1].val[i] = parameters.DMXranges[i];
             tempo2_psr[0].param[param_dmxr2].val[i] = parameters.DMXranges[i+1];
             tempo2_psr[0].param[param_dmxr1].paramSet[i] = 1;
             tempo2_psr[0].param[param_dmxr2].paramSet[i] = 1;
             tempo2_psr[0].param[param_dmx].paramSet[i] = 1;
        };

        tempo2_psr[0].param[param_fd].aSize = parameters.FD.size();
        for (i = 0 ; i < parameters.FD.size() ; i ++) tempo2_psr[0].param[param_fd].paramSet[i] = 1;

         tempo2_psr[0].param[param_raj].val[0]  = static_cast<double>(parameters.RA);
         tempo2_psr[0].param[param_decj].val[0] = static_cast<double>(parameters.DEC);
         tempo2_psr[0].param[param_pmra].val[0] = static_cast<double>(parameters.RA1);    // RA1 in mas of degree/yr
         tempo2_psr[0].param[param_pmdec].val[0] = static_cast<double>(parameters.DEC1) ; // DEC1 in mas of degree /yr
         tempo2_psr[0].param[param_pmrv].val[0] = static_cast<double>(parameters.distance1) ; // distance1 in mas of degree /yr
         tempo2_psr[0].param[param_dm].val[0] = static_cast<double>(parameters.DM) ; // pc/cm^3
         for (i = 0; i < parameters.DMX.size() ; i++) tempo2_psr[0].param[param_dmx].val[i] = parameters.DMX[i]; // pc / cm^3
         for (i = 0 ; i < parameters.FD.size() ; i ++) tempo2_psr[0].param[param_fd].val[i] = parameters.FD[i]; // pc / cm^3
         tempo2_psr[0].param[param_dm].val[1] = static_cast<double>(parameters.DM1) ; // pc / cm^3 / yr
         tempo2_psr[0].param[param_px].val[0] = 3261.563775885328 / static_cast<double>(parameters.distance) ; // light years -> mas

         tempo2_psr[0].param[param_posepoch].val[0] = parameters.posepoch + parameters.timeshift; // POSEPOCH in SSB time.
         tempo2_psr[0].param[param_dm].paramSet[0]=1 ; // Activate DM in tempo2
         tempo2_psr[0].param[param_dm].paramSet[1]=1 ; // Activate DM1 in tempo2

         veryFast = 0; // tempo2 global variable

         
         formBatsAll(tempo2_psr,tempo2_npsr);


        // Table allocations :
        toas = (value_type*) malloc( sizeof(value_type)*ntoas);// realloc(toas, sizeof(value_type)*ntoas);
        pulsardirection_BAT = (value_type**) realloc(pulsardirection_BAT, sizeof(value_type*) * ntoas );
        obs_ssb_BAT = (value_type**) realloc(obs_ssb_BAT, sizeof(value_type*) * ntoas );
        roemer_ss_BAT = (value_type*) realloc(roemer_ss_BAT, sizeof(value_type)*ntoas);
        turns = (turntype*) realloc(turns, sizeof(turntype)*ntoas);
        errors = (errortype*) realloc(errors, sizeof(errortype)*ntoas);
        weights = (errortype*) realloc(weights, sizeof(errortype)*ntoas);
        errortype totweight=0.;
        for (i = 0; i < ntoas ; i++)
        {
            toas[i] =  tempo2_psr[0].obsn[i].bat - parameters.timeshift; // toas Initialized in Initialize
            errors[i] = tempo2_psr[0].obsn[i].toaErr ;
            weights[i] = pow(  1. / errors[i], 2 ) ;
            totweight += weights[i];
            pulsardirection_BAT[i] = (value_type*) malloc( sizeof(value_type) * 3) ;
            obs_ssb_BAT[i] = (value_type*) malloc( sizeof(value_type) * 3) ;
        };
        for (i = 0; i < ntoas ; i++) weights[i] /= totweight;

     toes = (value_type*)realloc(toes, sizeof(value_type)*ntoas);
     residuals = (value_type*) realloc(residuals,sizeof(value_type)*ntoas);
     toa_in_interp = (long int*)realloc(toa_in_interp, sizeof(long int)*ntoas) ;



    ninterp = 0;
    sp = (value_type**) realloc(sp, sizeof(value_type*) * 1 );
    si = (value_type**) realloc(si, sizeof(value_type*) * 1 );
    so = (value_type**) realloc(so, sizeof(value_type*) * 1 );

    for (i = 0; i < 1 ; i++){
        sp[i] = NULL;
        si[i] = NULL;
        so[i] = NULL;
    }
    
    x0.resize(18 + parameters.nextra*6);  // allocate initial state vector
    
    first_interpolation_allocation = true;

    // For test mode
    testmode = false;


    spinPeriodmicrosec =  daysec * pow(dix, 6) / parameters.f ; // spin period in microseconds. doesn't change through fit, allows fixed conversion of phase difference

 // Switch to true spin frequency
 if (parameters.truefreq == 2 ) // assumes parameters are given for truefreq =0 (no einstein delay constant term in particular) and that we want to switch to truefreq=1 (all terms included)
 {
    printf("\n Changing spin frequency to take into account linear component of delays (eintein...) and truefreq changed to 1. \n" );
    parameters.f = parameters.f/(1.-freqshift) ;
    printf("New frequency = %.15Le \n\n", parameters.f);
    parameters.truefreq=1;
 } 
 printf("\n\n Out of Reconstruct_noCompute ! \n");

    return ;
}

void Fittriple::Reconstruct(char * parfile, char * datafile, const bool sorted_datafile)
{
  Reconstruct_noCompute(parfile, datafile, sorted_datafile);

  Compute_lnposterior(1);

  Compute_turns_from_parameters();

  Print();

  return;
}

void Fittriple::Reconstruct(char * parfile, char * datafile, const bool sorted_datafile, turntype * external_turns)
{
  Reconstruct_noCompute(parfile, datafile, sorted_datafile);
  Set_turns(external_turns);
  printf("\n Warning : A first call to Compute_lnposterior is needed to finalize intitialsation !\n");
  return;
}

void Fittriple::Reconstruct(char * parfile, char * datafile)
{
  Reconstruct(parfile, datafile, false);
  return;
}
void Fittriple::Reconstruct(char * parfile, char * datafile, turntype * external_turns)
{
  Reconstruct(parfile, datafile, false, external_turns);
  return;
}

void Fittriple::Initialize()
/*
 * This routine computes all the inner variables that depend on the parameters used .
 * It should be called every time parameters are changed . (it just need all the outer info : par file and data file)
 * It should ensure compatibility with different system of parameters (since the inner variables should be fairly independant on it)
 * Variables set :
 * spinfreq, spinfreq1, delaymax
 * dt_interp, ninterp, toa_in_interp, treference_in_interp
 *
 *  Compute end Initialize all the necessary stuff for integration of the equation of motion       * with the inherited properties from the "Integrateur" class;
 */
{
    long int i = 0L;
    long int k = 0L;


    valarray<value_type> rp(3);
    valarray<value_type> rpt(3);
    valarray<value_type> ri(3);
    valarray<value_type> rit(3);
    valarray<value_type> ro(3);
    valarray<value_type> rot(3);
    valarray<valarray<value_type>> r_extra(parameters.nextra);
    valarray<valarray<value_type>> v_extra(parameters.nextra);
    for (int i = 0; i < parameters.nextra; i++)
    {
      r_extra[i] = valarray<value_type>(3);
      v_extra[i] = valarray<value_type>(3);
    }

    parameters.Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);

    spinfreq = parameters.true_spinfreq ;
    spinfreq1 = parameters.true_spinfreq1 ;

    // Reset posepoch. Useful only to fit posepoch using Minuit_fittriple_interface_fitposepoch.h instead of  Minuit_fittriple_interface.h for tests..
    tempo2_psr[0].param[param_posepoch].val[0] = parameters.posepoch + parameters.timeshift; // POSEPOCH in SSB time.



  // Compute new toas with tempo2, if RA and DEC are defined and if they changed


        if ( (disable_sats == false) and ( previous_RA != parameters.RA or previous_DEC != parameters.DEC or previous_RA1 != parameters.RA1 or previous_DEC1 != parameters.DEC1
          or previous_distance1 != parameters.distance1   or previous_DM != parameters.DM or previous_DM1 != parameters.DM1 or previous_DMX != parameters.DMX or previous_FD != parameters.FD) )
        {

            tempo2_psr[0].param[param_raj].val[0]  = static_cast<double>(parameters.RA);
            tempo2_psr[0].param[param_decj].val[0] = static_cast<double>(parameters.DEC);
            tempo2_psr[0].param[param_pmra].val[0] = static_cast<double>(parameters.RA1); // RA1 in mas of degree/yr
            tempo2_psr[0].param[param_pmdec].val[0] = static_cast<double>(parameters.DEC1) ; // DEC1 in mas of degree /yr
            tempo2_psr[0].param[param_pmrv].val[0] = static_cast<double>(parameters.distance1) ; // distance1 in mas of degree /yr
            tempo2_psr[0].param[param_dm].val[0] = static_cast<double>(parameters.DM) ; // pc/cm^3
            for (i = 0; i < parameters.DMX.size() ; i++) tempo2_psr[0].param[param_dmx].val[i] = parameters.DMX[i]; // pc / cm^3
            for (i = 0 ; i < parameters.FD.size() ; i ++) tempo2_psr[0].param[param_fd].val[i] = parameters.FD[i]; // pc / cm^3
            tempo2_psr[0].param[param_dm].val[1] = static_cast<double>(parameters.DM1) ; // pc / cm^3 / yr
            tempo2_psr[0].param[param_px].val[0] = 3261.563775885328 / static_cast<double>(parameters.distance) ; // light years -> mas

            updateBatsAll(tempo2_psr,tempo2_npsr);  // Ne pas utiliser formBatsAll sinon fuite de mémoire !!! (probablement causé par l'appel à clock_corrections(psr,npsr);   ou ephemeris_routines(psr,npsr); ou extra_delays(psr,npsr);  ...
            for(i = 0; i < ntoas ; i++)
            {
                toas[i] =  tempo2_psr[0].obsn[i].bat - parameters.timeshift  ;
                if (strcmp(tempo2_psr[0].obsn[i].telID,"STL_FBAT")==0){ // Observatory_earth is actually observatory to ssb
                    for (k=0;k<3;k++) obs_ssb_BAT[i][k] = static_cast<value_type>(tempo2_psr[0].obsn[i].observatory_earth[k] ) * clight; // tempo expresses distances in lght-sec
                }
                else {
                    for (k=0;k<3;k++) obs_ssb_BAT[i][k] = static_cast<value_type>( tempo2_psr[0].obsn[i].earth_ssb[k] + tempo2_psr[0].obsn[i].observatory_earth[k] ) * clight; // tempo expresses distances in lght-sec
                }

                for (k=0;k<3;k++) pulsardirection_BAT[i][k] = tempo2_psr[0].obsn[i].psrPos[k];
                
                roemer_ss_BAT[i] = tempo2_psr[0].obsn[i].roemer;
            }



            previous_RA = parameters.RA ;
            previous_DEC = parameters.DEC ;
            previous_RA1 = parameters.RA1 ;
            previous_DEC1 = parameters.DEC1 ;
            previous_distance1 = parameters.distance1;
            previous_DM = parameters.DM;
            previous_DMX = parameters.DMX;
            previous_DM1 = parameters.DM1;
            previous_FD = parameters.FD;
        }


        // Use the following lines to replace toas by model generated toas and check that the residues are numerically compatible with 0.
//         if (first_interpolation_allocation == false) // ! Test !
//             {
//                 printf("\n\n using fake toas ! \n\n");
//                 Compute_fake_BATs_from_parameters(toas);
//                 for(i = 0; i < ntoas ; i++) toas[i] -= parameters.timeshift  ;
//                 parameters.dphase0 = 0.;
//             }
        initialpuslardirection[0] = tempo2_psr[0].posPulsar[0];
        initialpuslardirection[1] = tempo2_psr[0].posPulsar[1];
        initialpuslardirection[2] = tempo2_psr[0].posPulsar[2];


    delaymax = huit * ( parameters.ap + parameters.aB ) / clight / daysec ;//parameters.Pi * 16;//huit * ( parameters.ap + parameters.aB ) / clight / daysec ; // huit = security factor


    if (delaymax > (toas[ntoas - 1] - toas[0] ) * 0.01 )  printf("\n Warning : the delay margin in the interpolation grid is larger than 0.01 * time span . (Time span = %.5Lf, delaymax = %.5Lf) \n\n", toas[ntoas - 1] - toas[0], delaymax);


    //**** Reset interpolation properties if necessary

    long int oldninterp = 0 ;

    if (delaymaxold < unquart * delaymax or delaymaxold > huit * delaymax or forcerecomputeinterp == true)
    {

        delaymaxold = delaymax ;
        dt_interp = parameters.Pi / parameters.interpsteps_per_period_i;
        value_type tm = min( toas[0] - 2*delaymax, parameters.treference ) * parameters.DopplerF  ;   // Beginning of the interpolated interval
        value_type tp = max( toas[ntoas - 1] + 2*delaymax, parameters.treference ) * parameters.DopplerF  ;           // End of the interpolated interval
        tm = dt_interp * floor(tm / dt_interp ) ;
        tp = dt_interp * floor(tp / dt_interp ) ;
        oldninterp = ninterp;
        tm -= floor(parameters.interp_margin * parameters.interpsteps_per_period_i) * dt_interp ;
        tp += floor(parameters.interp_margin * parameters.interpsteps_per_period_i) * dt_interp ;
        ninterp = static_cast<long int>( floor ( ( tp - tm ) / dt_interp ) ) + 1L ;

        for (i = 0; i< ntoas ; ++i )  toa_in_interp[i] = static_cast<long int>( floor( ( toas[i]*parameters.DopplerF - tm ) / dt_interp ) ) ;
        treference_in_interp =  static_cast<long int>( floor( ( parameters.treference - tm ) / dt_interp ) ) ;

      // Affecting both tinterp and Integrateur.ts at the same time (to use less loops)

        tinterp = (value_type*) realloc(tinterp, ninterp * sizeof(value_type) ) ;
        ts.resize(ninterp) ;

        for (i = 0 ; i < ninterp ; ++i ) tinterp[i] = tm + i * dt_interp ;


        if (first_interpolation_allocation == true )
        {
            cout << " Computing and allocating the BAT interpolation grid ... " << endl;
        }
        else
        {
            cout << " Interpolation grid changed, reinitializing it ... " << endl;
            for (i = 0; i < oldninterp ; ++i)
            {
                free(sp[i]);
                free(si[i]);
                free(so[i]);
            }
        }


        sp = (value_type**) realloc(sp, sizeof(value_type*) * ninterp );
        si = (value_type**) realloc(si, sizeof(value_type*) * ninterp );
        so = (value_type**) realloc(so, sizeof(value_type*) * ninterp );


        for (i = 0; i < ninterp ; ++i)
        {
            sp[i] = (value_type*) malloc(sizeof(value_type) * 6 );
            si[i] = (value_type*) malloc(sizeof(value_type) * 6 );
            so[i] = (value_type*) malloc(sizeof(value_type) * 6 );
        }



        first_interpolation_allocation = false ;

    } //****** End of intialization after change of interpolation grid *****


       // Set caracteristic time and length scales
       if (parameters.Mo > 0.) {
        length = ( norm3d(rp) + norm3d(ri) + norm3d(ro) ) / trois ;
        timescale = ( norm3d(rp) / norm3d(rpt) + norm3d(ri) / norm3d(rit) + norm3d(ro) / norm3d(rot) ) / trois ;
       }
       else {
        length = ( norm3d(rp) + norm3d(ri)) / deux ;
        timescale = ( norm3d(rp) / norm3d(rpt) + norm3d(ri) / norm3d(rit) ) / deux ;
       }
        timedays = timescale / daysec;

       // Set other properties inherited from "Integrateur"
        n = ninterp ;
        t0 = parameters.treference / timedays ;
        dt0 = dt_interp / timedays / cinq;

    // ********* Pass on parameters of the equations of motion to the integrator**********************
        int_M0 = parameters.Mp ;
        int_M1 = parameters.Mi ;
        int_M2 = parameters.Mo ;
        int_SEP_D = parameters.SEP_D ;
        int_M_extra = parameters.M_extra;
        int_Gg = parameters.Gg;
        int_gammabar = parameters.gammabar;
        int_betabar = parameters.betabar;
	if (strcmp(parameters.specialcase,"")==0) 
        	for (i =0 ; i < (min(int_quadrupole.size(), parameters.quadrupole.size())) ; i++) int_quadrupole[i] = parameters.quadrupole_kgm2[i];

    // Add parameter here if your parameter is used as a parameter in the equations of motions (like masses)
    // *****************************************************************************

        // Initialize the other time arrays
        while (tinterp[treference_in_interp + 1] < zero   )
        {
            cout << "Initialize : tinterp[treference_in_interp + 1] < zero " << endl ;
            treference_in_interp++ ;
        }
        while (tinterp[treference_in_interp] > zero )
        {
            cout << "Initialize : tinterp[treference_in_interp + 1] > zero  " << endl ;
            treference_in_interp-- ;
        }

        ntneg = treference_in_interp + 1L ; // Gives the position in ts of the first tpos value, tpos[0]. And incidentally the size of tneg
        tneg.resize( ntneg);
        tpos.resize( ninterp - ntneg ) ;

        for (i = 0 ; i < ntneg ; ++i ) {
            ts[i] = tinterp[i] / timedays ;
            tneg[i] = abs( ts[i] ) ;
        }
        for (i = ntneg ; i < ninterp ; ++i ) {
            ts[i] = tinterp[i] / timedays ;
            tpos[i-ntneg] = ts[i] ;
        }


    // Set intial state vector for integration of the equation of motion
        x0[0] = rp[0] / length ;
        x0[1] = rp[1] / length ;
        x0[2] = rp[2] / length ;
        x0[3] = ri[0] / length ;
        x0[4] = ri[1] / length ;
        x0[5] = ri[2] / length ;
        x0[6] = ro[0] / length ;
        x0[7] = ro[1] / length ;
        x0[8] = ro[2] / length ;
        for (i = 0 ; i < parameters.nextra ; i ++) 
        {
            x0[9 + i*3] = r_extra[i][0] / length ;
            x0[10 + i*3] = r_extra[i][1] / length ;
            x0[11 + i*3] = r_extra[i][2] / length ;
        }
        
        extrashift = parameters.nextra * 3;
        
        x0[9 + extrashift] = rpt[0] / length  * timescale ;
        x0[10 + extrashift] = rpt[1] / length * timescale ;
        x0[11 + extrashift] = rpt[2] / length * timescale ;
        x0[12 + extrashift] = rit[0] / length * timescale ;
        x0[13 + extrashift] = rit[1] / length * timescale ;
        x0[14 + extrashift] = rit[2] / length * timescale ;
        x0[15 + extrashift] = rot[0] / length * timescale ;
        x0[16 + extrashift] = rot[1] / length * timescale ;
        x0[17 + extrashift] = rot[2] / length * timescale ;
        for (i = 0 ; i < parameters.nextra ; i ++) 
        {
            x0[18 + extrashift + i*3] = v_extra[i][0] / length * timescale  ;
            x0[19 + extrashift + i*3] = v_extra[i][1] / length * timescale  ;
            x0[20 + extrashift + i*3] = v_extra[i][2] / length * timescale  ;
        }

/* Test        Print();
        printf("\n\n allocation done %.15Le \n", x0[23]);
        
        valarray<value_type> center_of_mass(3);
        valarray<value_type> impulsion(3);
        value_type energy;
        valarray<value_type> sv(x0.size());
        for (i = 0 ; i < x0.size()/2; i++) sv[i] = x0[i]*length;
        for (i = x0.size()/2 ; i < x0.size(); i++) sv[i] = x0[i]*length/timescale;
        valarray<value_type> Ms(3+parameters.nextra);
        Ms[0] = int_M0;
        Ms[1] = int_M1;
        Ms[2] = int_M2;
        for (i = 0; i < parameters.nextra; i++)
        {
            Ms[i+3] = int_M_extra[i];
        }
        IntegralePrems3_1PN(sv, Ms,
                         0.,
                         int_Gg, int_gammabar, int_betabar,
                         center_of_mass, impulsion,
                         energy  );
        impulsion /= int_M0 + int_M1 + int_M2 + Ms[3];
        
        printf("CDMMMMMMMMMMMMMMM %.5Le %.5Le %.5Le \n\n", center_of_mass[0], center_of_mass[1],center_of_mass[2]);
        printf("VCDMMMMMMMMMMMMMM %.5Le %.5Le %.5Le \n\n", impulsion[0], impulsion[1],impulsion[2]);
        printf("energyyyyyyyyyyyy %.12Le\n\n", energy);
        printf("vextra %.5Le %.5Le %.5Le \n\n", v_extra[0][0] / length * timescale , v_extra[0][1] / length * timescale , v_extra[0][2] / length * timescale );
       */ 
// *** Uncomment the following and adjust the variables to test the numerical accuracy by integrating from a different reference time using an initial state vector found using a previous integration

//         value_type treftest= 8.1323044685202865184e+02 ; //
//         value_type nreftest = 100000;
//         ntneg = nreftest + 1L ; // Gives the position in ts of the first tpos value, tpos[0]. And incidentally the size of tneg
//         tneg.resize( ntneg);
//         tpos.resize( ninterp - ntneg ) ;
//
//         for (i = 0 ; i < ntneg ; ++i ) {
//             ts[i] = (tinterp[i]-treftest) / timedays ;
//             tneg[i] = abs( ts[i] ) ;
//         }
//         for (i = ntneg ; i < ninterp ; ++i ) {
//             ts[i] = (tinterp[i]-treftest) / timedays ;
//             tpos[i-ntneg] = ts[i] ;
//         }
//
//
//
// x0[0] = -1.1037274088965150786e-01 ;
// x0[1] = 3.3362342753036770110e-01 ;
// x0[2] = 3.1491364092357522516e-01 ;
// x0[3] = -9.2543389631955620479e-02 ;
// x0[4] = 2.8415546665539705753e-01 ;
// x0[5] = 3.5442495772645276833e-01 ;
// x0[6] = 4.3142209221411534307e-01 ;
// x0[7] = -1.3061890214020829701e+00 ;
// x0[8] = -1.2744436970888303432e+00 ;
// x0[9] = -8.1616934204317992522e-02 ;
// x0[10] = 2.7342922171934991177e-01 ;
// x0[11] = 7.6199490167639091163e-01 ;
// x0[12] = 1.0443327989083603982e+00 ;
// x0[13] = -3.2178732463913938266e+00 ;
// x0[14] = -4.1236212598358207766e+00 ;
// x0[15] = -2.1675807734585752104e-01 ;
// x0[16] = 5.9097096558377754875e-01 ;
// x0[17] = -6.8537984650355310924e-01 ;

//*** end of test



    return ;
}



void Fittriple::Set_parameters_relativeshift(double * relativeshift)
{
    parameters.Set_parameters_relativeshift (relativeshift);

    return;
}

void Fittriple::Initialise_parameters()
{
  valarray<value_type> rp(3);
  valarray<value_type> rpt(3);
  valarray<value_type> ri(3);
  valarray<value_type> rit(3);
  valarray<value_type> ro(3);
  valarray<value_type> rot(3);
  valarray<valarray<value_type>> r_extra(parameters.nextra);
  valarray<valarray<value_type>> v_extra(parameters.nextra);
  for (int i = 0; i < parameters.nextra; i++)
  {
      r_extra[i] = valarray<value_type>(3);
      v_extra[i] = valarray<value_type>(3);
  }
  parameters.Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
  
  return;
}

void Fittriple::Set_parameters_relativeshift(vector<double> relativeshift)
{

   parameters.Set_parameters_relativeshift(relativeshift);
   if (tracker == true)
    {
        printf("\n Values of parameter shifts : \n");
//         printf("%.8e    %.8e    \n", dRA, dDEC);
//         printf("%.8e    %.8e    \n", df, df1);
//         printf("%.8e    %.8e    %.8e    %.8e    %.8e    %.8e    \n", detap, dapsinii , dkappap, dapcosii, dtascp, dPi );
//         printf("%.8e    %.8e    %.8e    %.8e    %.8e    %.8e    \n", detaB, daBsinio , dkappaB, daBcosio, dtascB, dPo );
//         printf("%.8e    %.8e    %.8e    %.8e    %.8e\n",dmassparp, dmasspari, dmo, ddeltaoman, ddphase0);
//         printf("%.8e \n",dSEP_D);
    }
    return;

}


void Fittriple::Set_fitted_parameter_relativeshifts(vector<double> parameter_shifts)  // uses the filter "fitted_parameters"
{
   parameters.Set_fitted_parameter_relativeshifts(parameter_shifts);
   if (tracker == true)
    {
        printf("Call Sfpr %i : ", callnumber + 1);
        for (int i=0 ; i < parameters.fitted_parameters.size() ; i++ ) printf(" %i:%.4f ", parameters.fitted_parameters[i], parameter_shifts[i]);
        printf("\n");
    }
}


void Fittriple::Estimate_parameter_shift_scales(double chi2variationgoal)
{
    int i,j,p=0;
//     const int npars = parameters.nfitparams;
    double relativeshift = 1.;
    double chi2variation ;
    value_type delta = deux;
    value_type deltaini, deltaprev = zero;
    value_type lnpini = Compute_lnposterior(0);
    value_type rescale = deux;

    printf("\n Now estimating scalings of parameters to obtain a variation of %f of the chi2 for a variation of 1 of a parameter. The initial chi2 is %.5Lf .\n", chi2variationgoal, lnpini);


    for(p= 0; p< parameters.fitted_parameters.size()  ; p++) {
        j = 0;
        delta = deux;
        i = parameters.fitted_parameters[p];

        if (parameters.Mo == 0. and (i > 7 and i  < 14 or i == 16 or i == 17 or i == 19) ) continue; // Skip three-body parameters and SEP_D when only 2 bodies present



        if (parameters.Parameter_scale(i) <= 0. )
        {
            printf("Warning in Estimate_parameter_shift_scales : parameters.Parameter_scale(%i) <= 0., it's not gonna work ! Trying to set parameter %i scale to 1.", i, i);
        }

        if ( i == 22 )          // Distance parameter
        {
            relativeshift = -1.;
            chi2variation = chi2variationgoal / 1000. ;
            cout << " Special chi2 variation goal for distance " << chi2variation << endl;
        }
        else
        {
            relativeshift = 1.;
            chi2variation = chi2variationgoal  ;
        }

        parameters.Set_parameter_relativeshift(i, relativeshift);

        deltaini = abs(Compute_lnposterior(0) - lnpini);
        delta = deltaini;

        if (deltaini > chi2variation )
        {
            //parameters.parameter_shift_scales[i] /= deux;
            parameters.Set_parameter_scale(i, parameters.Parameter_scale(i) / rescale );
            while(delta > chi2variation and j <100)
            {
                parameters.Set_parameter_relativeshift(i, relativeshift);
                delta = abs(Compute_lnposterior(0) - lnpini);
           // cout << "delta bas " << delta << " " << i << "   "  ; printf("%.19Le \n", parameters[i] ) ;
                if (delta > chi2variation ) parameters.Set_parameter_scale(i, parameters.Parameter_scale(i) / rescale );
             //   cout << "i " << i << "  " << parameters.Parameter_scale(i) <<endl;
                j++;
            }
            if (j >= 100) cout << "Impossible to converge under 1 in Estimate_parameter_shift_scales for parameter " << i <<endl;
        }
        else
        {
            parameters.Set_parameter_scale(i, parameters.Parameter_scale(i) * rescale );
            while(delta < chi2variation and j <100)
            {
                parameters.Set_parameter_relativeshift(i, relativeshift);
                deltaprev = delta;
                delta = abs(Compute_lnposterior(0) - lnpini);
           //  cout << "delta agr " << delta << " " << i << "   "  ; printf("%.19Le \n", parameters[i] ) ;
                if (delta < chi2variation ) parameters.Set_parameter_scale(i, parameters.Parameter_scale(i) * rescale );
                j++;
            }
            if (j >= 100) cout << "Impossible to converge above 1 in Estimate_parameter_shift_scales for parameter " << i <<endl;
            //parameters.parameter_shift_scales[i] /= deux; // We keep a variation <= 1
            parameters.Set_parameter_scale(i, parameters.Parameter_scale(i) / rescale );
            delta = deltaprev;
        }

        parameters.Set_parameter_relativeshift(i, 0.); // relativeshift[i] = 0.;
        cout << "Parameter " << i << " ,  " << j << " iterations " << endl;
         cout << "     Initial variation : " << deltaini << ", Final variation : "  << delta << endl;
    }

    return;
}


void Fittriple::Create_mask(double threshold)//,  long int * & masked_toas, long int & mask_size)
/*
 * Reset the data by removing all the toas which current residual is larger than threshold (in microsec).
 * Also updates toes and residuals.
 * "masked_toas" of size "mask_size" is an array containing the indexes of the removed toas (so starting at 0). NULL pointer is returned if empty.
 */
{
    long int i,j,k,newntoas = 0L;
    long int mask_size = 0L;
    valarray<bool> mask(ntoas);
    valarray<value_type> oldtoes(toes, ntoas) ;
    valarray<value_type> oldarray(toas, ntoas) ;
    valarray<value_type> oldres(residuals, ntoas) ;
    valarray<turntype> oldturns(turns, ntoas) ;
    valarray<errortype> olderr(errors, ntoas) ;
    valarray<errortype> oldweights(weights, ntoas) ;


    // Create the mask
    newntoas = 0L;
    for(i = 0; i < ntoas ; ++i)
    {
        if (abs(residuals[i]) >= threshold)
            mask[i] = false;
        else
        {
            mask[i] = true;
            ++newntoas;
        }
    }

    // Create the masked arrays
    toas = (value_type*) realloc(toas, newntoas * sizeof(value_type) ) ;
    turns = (turntype*) realloc(turns, newntoas * sizeof(turntype) ) ;
    errors = (errortype*) realloc(errors, newntoas * sizeof(errortype) ) ;
    weights = (errortype*) realloc(weights, newntoas * sizeof(errortype) ) ;
    toes = (value_type*) realloc(toes, newntoas * sizeof(value_type) ) ;
    residuals = (value_type*) realloc(residuals, newntoas * sizeof(value_type) ) ;

    mask_size = ntoas - newntoas ;
//     if ( mask_size > 0L )
//         masked_toas = new long int[mask_size] ;
//     else
//         masked_toas = NULL;
    removed_toas.resize(mask_size);

    j = 0L; k = 0L;
    // Affect the masked arrays
    for (i = 0 ; i < ntoas ; ++i )
    {
        if ( mask[i] == true )
        {
            toas[j] = oldarray[i] ;
            turns[j] = oldturns[i] ;
            errors[j] = olderr[i] ;
            weights[j] = oldweights[i] ;
            toes[j] = oldtoes[j];
            residuals[j] = oldres[j];
            ++j;
        }
        else
        {
            removed_toas[k] = i;
            ++k;
        }
    }

    ntoas = newntoas ;

    Initialize();

    return ;

}



void Fittriple::Print_mask(double threshold)
{
    long int i = 0L;
    long int toremove = 0L;

    printf("\n# Number  SAT  Residual(toas which abs(residual) > %f microsec )\n", threshold);
    for(i = 0; i < ntoas ; ++i)
    {
        if (abs(residuals[i]) > threshold)
        {
            printf("%li %.19Lf %.19Le\n",i+1, tempo2_psr[0].obsn[i].sat, residuals[i] ) ;
            toremove ++;
        }
    }
    printf("Total number of toas to remove : %li out of %li \n", toremove, ntoas);

    return ;

}


void Fittriple::Compute_turns_from_parameters()
/*
 * Replace the turn numbers by those computed from the current times of emission (so Compute_lnposterior must have been run).
 * The first toa correspond to turn 0.
 */
{
    long int i = 0L;
    // Look for the position of the time of reference
   // long int trefintoas = Rank_in_sorted_array(treference, toas, 0L, ntoas-1L);
    value_type phase0 = toes[0] - parameters.treference ;
    value_type phase = zero ;
    phase0 = spinfreq * phase0 + undemi * spinfreq1 * pow( phase0, 2 ) ;
    phase0 = phase0 - dix - undemi ; // dix pour éviter que phase - phase0 ne soit négatif, un demi pour avoir la bonne partie entière (pas celle du tour d'avant)
    for (i = 0; i < ntoas ; ++i ) {
        phase = toes[i] - parameters.treference ;
        phase = spinfreq * phase + undemi * spinfreq1 * pow( phase, 2 ) ;
        turns[i] =  static_cast<turntype>( floor( phase - phase0 )  - dix ) ;
     }


     return ;
}


void Fittriple::Compute_fake_BATs_from_parameters(value_type * fake_BATs)
/*
Should be accurate including for geometric delays, contrary to Compute_fake_BATS_and_delays_from_parameters
But does not return details of delays 
 */
{
    long int i = 0L;
    value_type turn = 0L;
    value_type t0 = 0.;
    value_type t1 = 0.;   
    value_type eps = 1.e-14; // precision to achieve in days in inversion loop
    // Look for the position of the time of reference
   // long int trefintoas = Rank_in_sorted_array(treference, toas, 0L, ntoas-1L);
    long int marginmin = static_cast<long int>(floor(delaymax / dt_interp) )  + 2L ; // We assume that the delays will never be greater than 0.1 days
    long int ntisaround =  2* (floor(parameters.interp_margin * parameters.interpsteps_per_period_i) +  floor(delaymax / dt_interp ) ) ; // Number of interpolation steps to take around each toa to compute geometric delays. must be even.
    value_type * delay ;
    value_type * delay_geom;
    delay = (value_type*) malloc( sizeof(value_type) * ninterp );
    delay_geom = (value_type*) malloc( sizeof(value_type) * ntisaround );
    value_type * oldtoas ;
    oldtoas = (value_type*) malloc( sizeof(value_type) * ntoas );
    value_type pospsr[3];
    value_type velpsr[3];
    value_type maxdiff=zero;
    value_type oldBAT,newBAT;
    bool nanflag= false;

    for (i = 0; i < 3 ; i ++)
     {
        pospsr[i] = static_cast<value_type>(tempo2_psr[0].posPulsar[i]);
        velpsr[i] = static_cast<value_type>(tempo2_psr[0].velPulsar[i])/100.;
     }


    value_type phase0 = toes[0] - parameters.treference ;
    value_type phase = zero ;
    phase0 = spinfreq * phase0 + undemi * spinfreq1 * pow( phase0, 2 ) + parameters.dphase0 ;
    phase0 = phase0 - dix - undemi ; // dix pour éviter que phase - phase0 ne soit négatif, un demi pour avoir la bonne partie entière (pas celle du tour d'avant)
    for (i = 0; i < ntoas ; ++i ) {
        phase = toes[i] - parameters.treference ;
        phase = spinfreq * phase + undemi * spinfreq1 * pow( phase, 2 ) ;
        turn =  ( floor( phase - phase0 ) - dix ) + (phase0 + dix + undemi)  ;
        // Inverse the timing formula. Rmk : the exact formula is not practical because of rounding errors
        t0 = 0.;
        t1 = turn / spinfreq ;
        do
        {
            t0 = t1;
            t1 = turn / spinfreq - spinfreq1 /(2.*spinfreq) * pow(t0,2);
        } while (abs(t1 - t0) / t1 > 1.e-17 ) ;
        t0 = t1;
        t1 = turn / spinfreq - spinfreq1 /(2.*spinfreq) * pow(t0,2);
        // End of inversion
        fake_BATs[i] = t1;// At this stage fake_BATs = fake time of emission
        fake_BATs[i] += parameters.treference ;

        
     }
//** Uncomment to check that the following inversion procedure is ale to give back the original bats
      //  for(i=0; i< ntoas ; i++) fake_BATs[i] = toes[i];

 // Compute Einstein delay interpolation
         Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pospsr,
                 int_M0, int_M1, int_M2, parameters.f,
                 parameters.einstein, false, false,
                 delay,  parameters.truefreq, freqshift, nanflag, NULL
                ) ;

         Spline delay_einstein_i( tinterp, delay, ninterp);

         for(i=0; i< ntoas ; i++)
         {
            t0 = -1.;
            t1 = fake_BATs[i];
            while (abs(t0 - t1) > eps )
            {
                t0 = t1;
                t1 = fake_BATs[i] + delay_einstein_i(t1, 0, ninterp);
            };
            fake_BATs[i] = t1;
         };


 //   Spline delays_interpolation( tinterp, delay, ninterp );

// Compute Shapiro and aberration delays interpolation
    Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pospsr,
                 int_M0, int_M1, int_M2, parameters.f,
                 false, parameters.shapiro, parameters.aberration,
                 delay,  parameters.truefreq, freqshift, nanflag, NULL
                ) ;

    Spline delays_interpolation( tinterp, delay, ninterp );

    long int firsttis = 0L;
    Spline delay_geom_i (ntisaround);
    
//-------------- !!!! Double sinusoid model or RN_PL model!!!!---------------------------------------------
    value_type ds_f = parameters.quadrupole[0];
    value_type ds_A1 = parameters.quadrupole[1]/(daysec * pow(dix, 6));
    value_type ds_phi1 = parameters.quadrupole[2];
    value_type ds_A2 = parameters.quadrupole[3]/(daysec * pow(dix, 6));
    value_type ds_phi2 = parameters.quadrupole[4];
    value_type ds_A3 =0.;
    value_type ds_phi3 = 0.;
    value_type ds_A4 =0.;
    value_type ds_phi4 = 0.;
    value_type ds_A5 =0.;
    value_type ds_phi5 = 0.;
    value_type gam=0.;    
    if (strcmp(parameters.specialcase,"RN_PL")==0)
    {
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
    }
    else if (strcmp(parameters.specialcase,"RN_PL")!=0)
	printf(" Error : special case '%s' not implemented in Compute_fake_BATs_from_parameters", parameters.specialcase);
//-----------------------------------------------------------------------------------------    
    for (i = 0; i < ntoas ; ++i )
    {
        firsttis = toa_in_interp[i] - ntisaround / 2 ;
        oldBAT = toas[i]+1.;
        newBAT = toas[i];
        while (abs(newBAT - oldBAT) > eps)
        {
          Delays_Brut_geometric_local(tinterp + firsttis, newBAT,    ntisaround,
                 sp + firsttis, si + firsttis, so + firsttis, pospsr, velpsr, parameters.posepoch, obs_ssb_BAT[i],
                 parameters.distance, parameters.distance1,
                 delay_geom, parameters.kopeikin, parameters.shklovskii, (i==100), NULL
                )  ;


        delay_geom_i.ReSpline( tinterp + firsttis, delay_geom ) ;
        oldBAT = newBAT;
        newBAT = fake_BATs[i] + delays_interpolation( fake_BATs[i], toa_in_interp[i] - marginmin,
                                            toa_in_interp[i] + marginmin )
                              + delay_geom_i(fake_BATs[i], ntisaround/2 - marginmin,
                                            ntisaround/2 + marginmin) ;
//-------------- !!!! Double sinusoid model or RN_PL model !!!!---------------------------------------
        if (strcmp(parameters.specialcase,"RN_PL")==0)        
            newBAT += ds_A1 *sin(oldBAT *ds_f + ds_phi1) + ds_A2 *sin(2 * oldBAT *ds_f + ds_phi2) + ds_A3 *sin(3 * oldBAT *ds_f + ds_phi3) + ds_A4 *sin(4 * oldBAT *ds_f + ds_phi4) + ds_A5 *sin(5 * oldBAT *ds_f + ds_phi5);
//------------------------------------------------------------------------------------
        // newBAT = fake_BATs[i] + delay_geom_i(fake_BATs[i], ntisaround/2 - marginmin,
        //                                     ntisaround/2 + marginmin) ;
        // if (i==100)
        // {
        //   printf("newbat - oldbat newbat intert %.5Le %.19Le %.19Le\n\n", newBAT-oldBAT, newBAT, fake_BATs[i]);
        //   printf("nogeom %.19Le\n",delays_interpolation( fake_BATs[i], toa_in_interp[i] - marginmin,
        //                                       toa_in_interp[i] + marginmin ));
        //   printf("geom %.19Le\n",  delay_geom_i(fake_BATs[i], ntisaround/2 - marginmin,
        //                                     ntisaround/2 + marginmin));
        // };
       }

        fake_BATs[i] = newBAT;
        fake_BATs[i] += parameters.timeshift ;

    }     
    
    maxdiff=zero;
    for (i = 0; i < ntoas ; ++i ) maxdiff = max(abs(toas[i] + parameters.timeshift - fake_BATs[i]), maxdiff);
    printf("\n+++ maxdiff bats - fakeBats = %.5Le \n\n", maxdiff);

//** Check that fake_BATs are accurate
    for (i = 0; i < ntoas ; ++i )
    {
        oldtoas[i] = toas[i];
        toas[i] = fake_BATs[i] - parameters.timeshift;
    }
    disable_sats = true; // prevent SATs being recomputed
    Compute_lnposterior(1);
    Compute_turns_from_parameters();
    Save_output_timing_data("datafile_Fake_Bats.dat");
    Save_turn_numbers("turns_Fake_Bats.dat");
    maxdiff= zero;
    for (i = 0; i < ntoas ; ++i ) maxdiff = max(abs(residuals[i]), maxdiff);
    printf("\n+++ maxdiff residuals with fakeBats = %.5Le microsec\n\n", maxdiff);
    // Return to the original toas
    for (i = 0; i < ntoas ; ++i ) toas[i] = oldtoas[i];
    Compute_lnposterior(1);
    Compute_turns_from_parameters();
    disable_sats =false; //*/
//** End of check

    free(delay);
    free(delay_geom);
    free(oldtoas);

    return ;
}




void Fittriple::Compute_fake_SATs_from_BATs(value_type * BATs, value_type * SATs)
{
    value_type maxdiff = zero;
    value_type * oldsats ;

    oldsats = (value_type*) malloc( sizeof(value_type) * ntoas );
    for (int i=0; i < ntoas ; i++) oldsats[i]= tempo2_psr[0].obsn[i].sat  ;


    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i]   ;
    formBatsAll(tempo2_psr,tempo2_npsr);

    // This should be a loop but formabatsall bugs when called from within a loop !! (don't know why...)
    // Number of calls should be more than enough though
    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
    formBatsAll(tempo2_psr,tempo2_npsr);
    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
    formBatsAll(tempo2_psr,tempo2_npsr);
    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
    formBatsAll(tempo2_psr,tempo2_npsr);
    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
    formBatsAll(tempo2_psr,tempo2_npsr);
    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
    formBatsAll(tempo2_psr,tempo2_npsr);
    maxdiff = zero;
    for (int i=0; i < ntoas ; i++)
        {
            if (abs(tempo2_psr[0].obsn[i].bat - BATs[i]) > maxdiff ) maxdiff = abs(tempo2_psr[0].obsn[i].bat - BATs[i]) ;
        }
    if (maxdiff > pow(10.,-15) ) printf("\nWarning: in Compute_fake_SATs_from_BATs convergence not reached, consider adding steps to  the loop. Maximum error = %Le days \n\n", maxdiff);

    for (int i=0; i < ntoas ; i++) SATs[i] = tempo2_psr[0].obsn[i].sat;

    maxdiff=0.;
    for (int i=0; i < ntoas ; i++)
        {
            if (abs(SATs[i] - oldsats[i]) > maxdiff ) maxdiff = abs(SATs[i] - oldsats[i]) ;
        }
//      printf("\n--------------- = %Le days \n\n", maxdiff);

    for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = oldsats[i]   ;
    formBatsAll(tempo2_psr,tempo2_npsr);

    maxdiff=0.;
    for (int i=0; i < ntoas ; i++)
        {
            if (abs(SATs[i] - oldsats[i]) > maxdiff ) maxdiff = abs(SATs[i] - oldsats[i]) ;
        }
//     printf("\n--------------- = %Le days \n\n", maxdiff);


    free(oldsats);

    return;

}




void Fittriple::Compute_fake_BATS_and_delays_from_parameters(long int nfakeBats, value_type * fake_BATs, value_type * fake_delay_geom, value_type * fake_delay_ein,
                                                                                     value_type * fake_delay_shap, value_type * fake_delay_aber)
/*
 Based on Compute_fake_BATs_from_parameters.
 All arrays have length nfakeBats
 fake_BATs contains initial guesses for the toes, is rounded to a neighbouring possible toe, and then converted and returned as the corresponding toa assuming
 the current parameters.
 The delays are separated into the different arrays, and correspond to each final fake_BATs bat time.
 fake_BATs  in days
 delays in seconds
 For the geometric delay, the coupling with the position of the earth is neglected, and proper motion is approximated.
 */
{
    long int i = 0L;
    value_type turn = 0L;
    value_type t0 = 0.;
    value_type t1 = 0.;
    // Look for the position of the time of reference
   // long int trefintoas = Rank_in_sorted_array(treference, toas, 0L, ntoas-1L);
    long int marginmin = static_cast<long int>(floor(delaymax / dt_interp) )  + 2L ; // We assume that the delays will never be greater than 0.1 days
    long int ntisaround =  2* (floor(parameters.interp_margin * parameters.interpsteps_per_period_i) +  floor(delaymax / dt_interp ) ) ; // Number of interpolation steps to take around each toa to compute geometric delays. must be even.
    value_type fshift = 0. ; // useless in this routine, just to give an argument to Delays_Brut_nogeometric
    value_type * delay ;
    value_type * delay_geom;
    delay = (value_type*) malloc( sizeof(value_type) * ninterp );
    delay_geom = (value_type*) malloc( sizeof(value_type) * ntisaround );
    value_type pospsr[3];
    value_type velpsr[3];
    bool nanflag = false;

    value_type fake_obs_ssb[3];
    value_type * delay_ein ;
    value_type * delay_shap ;
    value_type * delay_aber ;

    delay_ein = (value_type*) malloc( sizeof(value_type) * ninterp );
    delay_shap = (value_type*) malloc( sizeof(value_type) * ninterp );
    delay_aber = (value_type*) malloc( sizeof(value_type) * ninterp );

    for (i = 0; i < 3 ; i ++)
     {
        pospsr[i] = static_cast<value_type>(tempo2_psr[0].posPulsar[i]);
        velpsr[i] = static_cast<value_type>(tempo2_psr[0].velPulsar[i])/ static_cast<value_type>(100); // convertion from rad/century to rad/year. 
        fake_obs_ssb[i] = 0.;
     }


    value_type phase0 = toes[0] - parameters.treference ;
    value_type phase = zero ;
    phase0 = spinfreq * phase0 + undemi * spinfreq1 * pow( phase0, 2 ) + parameters.dphase0 ;
    phase0 = phase0 - dix - undemi ; // dix pour éviter que phase - phase0 ne soit négatif, un demi pour avoir la bonne partie entière (pas celle du tour d'avant)
    for (i = 0; i < nfakeBats ; ++i ) { // Determine the time of emission corresponding to an integer N the closest to each fake_BATs[i] and put it into fake_BAts[i]
        phase = fake_BATs[i] - parameters.treference ;
        phase = spinfreq * phase + undemi * spinfreq1 * pow( phase, 2 ) ;
        turn =  ( floor( phase - phase0 ) - dix ) + (phase0 + dix + undemi)  ;
        // Inverse the timing formula. Rmk : the exact formula is not practical because of rounding errors
        t0 = 0.;
        t1 = turn / spinfreq ;
        do
        {
            t0 = t1;
            t1 = turn / spinfreq - spinfreq1 /(2.*spinfreq) * pow(t0,2);
        } while (abs(t1 - t0) / t1 > pow(dix, -17) ) ;
        t0 = t1;
        t1 = turn / spinfreq - spinfreq1 /(2.*spinfreq) * pow(t0,2);
        // End of inversion
        fake_BATs[i] = t1;//( - spinfreq + sqrt(pow(spinfreq,2) + 2. * spinfreq1 * turn ) ) / spinfreq1 ;
        fake_BATs[i] += parameters.treference ;
     }

// this one to compute the BATs
    Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
                 int_M0, int_M1, int_M2, parameters.f,
                 parameters.einstein, parameters.shapiro, parameters.aberration,
                 delay,  parameters.truefreq, fshift, nanflag, NULL
                ) ;

// These three to keep track of the delays individually
    Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
                 int_M0, int_M1, int_M2, parameters.f,
                 parameters.einstein, false, false, //parameters.shapiro, parameters.aberration,
                 delay_ein,  parameters.truefreq, fshift, nanflag, NULL
                ) ;

    Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
                 int_M0, int_M1, int_M2, parameters.f,
                 false, parameters.shapiro, false, //parameters.shapiro, parameters.aberration,
                 delay_shap,  parameters.truefreq, fshift, nanflag, NULL
                ) ;

    Delays_Brut_nogeometric(tinterp, parameters.treference,
                 ninterp, treference_in_interp,
                 sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
                 int_M0, int_M1, int_M2, parameters.f,
                 false, false, parameters.aberration,
                 delay_aber,  parameters.truefreq, fshift, nanflag, NULL
                ) ;
        Spline delays_interpolation_ein( tinterp, delay_ein, ninterp );
            Spline delays_interpolation_shap( tinterp, delay_shap, ninterp );
                Spline delays_interpolation_aber( tinterp, delay_aber, ninterp );

    long int firsttis = 0L;
    long int tintinterp = 0L;
    long int maxtis = 0L;
    Spline delays_interpolation( tinterp, delay, ninterp );
    Spline delay_geom_i (ntisaround);

    for (i = 0; i < nfakeBats ; ++i ) // Go from time of emission to time of arrival and record the delays
    {
        while (tinterp[tintinterp] < fake_BATs[i] and tintinterp < ninterp-1) {
            tintinterp += 1;
        }
        if (tintinterp == ninterp -1 or tintinterp == 0 ) {
            printf("Warning : fake_BATs may be out of range of tinterp in Compute_fake_BATS_and_delays_from_parameters !");
        };
        firsttis = max(tintinterp - ntisaround / 2, 0) ;

        maxtis =   min(firsttis + ntisaround, ninterp) - firsttis ;// Need to be defined because the tinterp span may be too narrow.

        if (maxtis != ntisaround ) cout << "Error : maxtis != ntisaround" << maxtis << " " << firsttis << endl;

        Delays_Brut_geometric_local(tinterp + firsttis, fake_BATs[i],  maxtis,
                 sp + firsttis, si + firsttis, so + firsttis, pospsr, velpsr, parameters.posepoch, fake_obs_ssb, 
                 parameters.distance, parameters.distance1,
                 delay_geom,parameters.kopeikin, parameters.shklovskii, 0., NULL
                )  ;


        delay_geom_i.ReSpline( tinterp + firsttis, delay_geom ) ; // !!!!!! Will make problems anyway if maxtis < ntisaround

        fake_delay_geom[i] = delay_geom_i(fake_BATs[i], 0,
                                            ntisaround)  * 86400. ;
        fake_delay_ein[i] = delays_interpolation_ein( fake_BATs[i], max(0,tintinterp - marginmin),
                                            min(tintinterp +marginmin, ninterp) )  * 86400. ;
        fake_delay_shap[i] = delays_interpolation_shap( fake_BATs[i],max(0,tintinterp - marginmin),
                                            min(tintinterp +marginmin, ninterp) )  * 86400. ;
        fake_delay_aber[i] = delays_interpolation_aber( fake_BATs[i], max(0,tintinterp - marginmin),
                                            min(tintinterp +marginmin, ninterp) )   * 86400. ;

        fake_BATs[i] += delays_interpolation( fake_BATs[i],max(0,tintinterp - marginmin),
                                            min(tintinterp +marginmin, ninterp) ) +
                        delay_geom_i(fake_BATs[i], 0,
                                            ntisaround) ;
       // fake_BATs[i] += parameters.timeshift ;


    }
    free(delay);
    free(delay_geom);
    free(delay_ein);
    free(delay_aber);
    free(delay_shap);

    return ;
}



//
//
// void Fittriple::Compute_fake_SATS_and_delays_from_parameters(long int nfakeBats, value_type * fake_BATs, value_type * fake_delay_geom, value_type * fake_delay_ein,
//                                                                                      value_type * fake_delay_shap, value_type * fake_delay_aber)
// /*
//  Based on Compute_fake_BATs_from_parameters.
//  All arrays have length nfakeBats
//  fake_BATs contains initial guesses for the toes, is rounded to a neighbouring possible toe, and then converted and returned as the corresponding toa assuming
//  the current parameters.
//  The delays are separated into the different arrays, and correspond to each final fake_BATs bat time.
//  For the geometric delay, the coupling with the position of the earth is neglected, and proper motion is approximated.
//  */
// {
//     long int i = 0L;
//     value_type turn = 0L;
//     value_type t0 = 0.;
//     value_type t1 = 0.;
//     // Look for the position of the time of reference
//    // long int trefintoas = Rank_in_sorted_array(treference, toas, 0L, ntoas-1L);
//     long int marginmin = static_cast<long int>(floor(delaymax / dt_interp) )  + 2L ; // We assume that the delays will never be greater than 0.1 days
//     long int ntisaround =  2* (floor(parameters.interp_margin * parameters.interpsteps_per_period_i) +  floor(delaymax / dt_interp ) ) ; // Number of interpolation steps to take around each toa to compute geometric delays. must be even.
//     value_type fshift = 0. ; // useless in this routine, just to give an argument to Delays_Brut_nogeometric
//     value_type * delay ;
//     value_type * delay_geom;
//     delay = (value_type*) malloc( sizeof(value_type) * ninterp );
//     delay_geom = (value_type*) malloc( sizeof(value_type) * ntisaround );
//     value_type pospsr[3];
//     value_type velpsr[3];
//
//     value_type fake_obs_ssb[3];
//     value_type fake_roemer_ss =0. ;
//     value_type * delay_ein ;
//     value_type * delay_shap ;
//     value_type * delay_aber ;
//
//     delay_ein = (value_type*) malloc( sizeof(value_type) * ninterp );
//     delay_shap = (value_type*) malloc( sizeof(value_type) * ninterp );
//     delay_aber = (value_type*) malloc( sizeof(value_type) * ninterp );
//
//     for (i = 0; i < 3 ; i ++)
//      {
//         pospsr[i] = static_cast<value_type>(tempo2_psr[0].posPulsar[i]);
//         velpsr[i] = static_cast<value_type>(tempo2_psr[0].velPulsar[i]);
//         fake_obs_ssb[i] = 0.;
//      }
//
//
// // OLD MAKE SATS
//
//
// value_type maxdiff = zero;
// int k = 0;
// value_type * oldsats ;
//
// oldsats = (value_type*) malloc( sizeof(value_type) * ntoas );
// for (int i=0; i < ntoas ; i++) oldsats[i]= tempo2_psr[0].obsn[i].sat  ;
//
//
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i]   ;
// formBatsAll(tempo2_psr,tempo2_npsr);
//
// // This should be a loop but formabatsall bugs when called from within a loop !! (don't know why...)
// // Number of calls should be more than enough though
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
// formBatsAll(tempo2_psr,tempo2_npsr);
// for (i = 0; i < nfakeBats ; ++i ) // Go from time of emission to time of arrival and record the delays
// {
//     while (tinterp[tintinterp] < fake_BATs[i] and tintinterp < ninterp-1) {
//         tintinterp += 1;
//     }
//     if (tintinterp == ninterp -1 or tintinterp == 0 ) {
//         printf("Warning : fake_BATs may be out of range of tinterp in Compute_fake_BATS_and_delays_from_parameters !");
//     };
//     firsttis = max(tintinterp - ntisaround / 2, 0) ;
//     maxtis =   min(firsttis + ntisaround, ninterp) - firsttis ;// Need to be defined because the tinterp span may be too narrow.
//
//     Delays_Brut_geometric_local(tinterp + firsttis, fake_BATs[i],  maxtis,
//              sp + firsttis, si + firsttis, so + firsttis, pospsr, velpsr, parameters.posepoch, fake_obs_ssb, //fake_roemer_ss,
//              parameters.distance, parameters.distance1,
//              delay_geom, parameters.truefreq, parameters.kopeikin, parameters.shklovskii, 0.
//             )  ;
//
//     delay_geom_i.ReSpline( tinterp + firsttis, delay_geom ) ; // !!!!!! Will make problems anyway if maxtis < ntisaround
//
//     newfake_BATs[i] = fakeBats_base[i] + delay_geom_i(fake_BATs[i], 0, ntisaround) ;
// }
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
// formBatsAll(tempo2_psr,tempo2_npsr);
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
// formBatsAll(tempo2_psr,tempo2_npsr);
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
// formBatsAll(tempo2_psr,tempo2_npsr);
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = BATs[i] - (tempo2_psr[0].obsn[i].bat - tempo2_psr[0].obsn[i].sat);
// formBatsAll(tempo2_psr,tempo2_npsr);
// maxdiff = zero;
// for (int i=0; i < ntoas ; i++)
//     {
//         if (abs(tempo2_psr[0].obsn[i].bat - BATs[i]) > maxdiff ) maxdiff = abs(tempo2_psr[0].obsn[i].bat - BATs[i]) ;
//     }
// if (maxdiff > pow(10.,-15) ) printf("\nWarning: in Compute_fake_SATs_from_BATs convergence not reached, consider adding steps to  the loop. Maximum error = %Le days \n\n", maxdiff);
//
// for (int i=0; i < ntoas ; i++) SATs[i] = tempo2_psr[0].obsn[i].sat;
// //      printf("-- %.15Lf %.15Lf %.15Lf %.15Lf\n", SATs[0], SATs[1], SATs[2], SATs[3]);
// //     printf("--- %.15Lf \n", tempo2_psr[0].obsn[0].sat);
//
// maxdiff=0.;
// for (int i=0; i < ntoas ; i++)
//     {
//         if (abs(SATs[i] - oldsats[i]) > maxdiff ) maxdiff = abs(SATs[i] - oldsats[i]) ;
//     }
// //      printf("\n--------------- = %Le days \n\n", maxdiff);
//
// for (int i=0; i < ntoas ; i++) tempo2_psr[0].obsn[i].sat = oldsats[i]   ;
// formBatsAll(tempo2_psr,tempo2_npsr);
//
// maxdiff=0.;
// for (int i=0; i < ntoas ; i++)
//     {
//         if (abs(SATs[i] - oldsats[i]) > maxdiff ) maxdiff = abs(SATs[i] - oldsats[i]) ;
//     }
// //     printf("\n--------------- = %Le days \n\n", maxdiff);
//
//
// free(oldsats);
// // HERE STARTS THE OLD MAKE_BATS
//     value_type phase0 = toes[0] - parameters.treference ;
//     value_type phase = zero ;
//     phase0 = spinfreq * phase0 + undemi * spinfreq1 * pow( phase0, 2 ) + parameters.dphase0 ;
//     phase0 = phase0 - dix - undemi ; // dix pour éviter que phase - phase0 ne soit négatif, un demi pour avoir la bonne partie entière (pas celle du tour d'avant)
//     for (i = 0; i < nfakeBats ; ++i ) { // Determine the time of emission corresponding to an integer N the closest to each fake_BATs[i] and put it into fake_BAts[i]
//         phase = fake_BATs[i] - parameters.treference ;
//         phase = spinfreq * phase + undemi * spinfreq1 * pow( phase, 2 ) ;
//         turn =  ( floor( phase - phase0 ) - dix ) + (phase0 + dix + undemi)  ;
//         // Inverse the timing formula. Rmk : the exact formula is not practical because of rounding errors
//         t0 = 0.;
//         t1 = turn / spinfreq ;
//         do
//         {
//             t0 = t1;
//             t1 = turn / spinfreq - spinfreq1 /(2.*spinfreq) * pow(t0,2);
//         } while (abs(t1 - t0) / t1 > pow(dix, -17) ) ;
//         t0 = t1;
//         t1 = turn / spinfreq - spinfreq1 /(2.*spinfreq) * pow(t0,2);
//         // End of inversion
//         fake_BATs[i] = t1;//( - spinfreq + sqrt(pow(spinfreq,2) + 2. * spinfreq1 * turn ) ) / spinfreq1 ;
//         fake_BATs[i] += parameters.treference ;
//         //printf("fake_bat 1 : %.19Le %.19Le  %.19Le %.5Le %.19Le %.19Le \n", parameters.timeshift, fake_BATs[i], parameters.timeshift + fake_BATs[i], turn,  parameters.timeshift +( - spinfreq + sqrt(pow(spinfreq,2) + 2. * spinfreq1 * turn ) ) / spinfreq1, spinfreq1);
//      }
// //printf("\n\n");
//
// // this one to compute the BATs
//     Delays_Brut_nogeometric(tinterp, parameters.treference,
//                  ninterp, treference_in_interp,
//                  sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
//                  int_M1, int_M2, parameters.f,
//                  parameters.einstein, parameters.shapiro, parameters.aberration,
//                  delay,  parameters.truefreq, fshift, NULL
//                 ) ;
//
// // These three to keep track of the delays individually
//     Delays_Brut_nogeometric(tinterp, parameters.treference,
//                  ninterp, treference_in_interp,
//                  sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
//                  int_M1, int_M2, parameters.f,
//                  parameters.einstein, false, false, //parameters.shapiro, parameters.aberration,
//                  delay_ein,  parameters.truefreq, fshift, NULL
//                 ) ;
//
//     Delays_Brut_nogeometric(tinterp, parameters.treference,
//                  ninterp, treference_in_interp,
//                  sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
//                  int_M1, int_M2, parameters.f,
//                  false, parameters.shapiro, false, //parameters.shapiro, parameters.aberration,
//                  delay_shap,  parameters.truefreq, fshift, NULL
//                 ) ;
//
//     Delays_Brut_nogeometric(tinterp, parameters.treference,
//                  ninterp, treference_in_interp,
//                  sp, si, so, pulsardirection_BAT[static_cast<long int>(floor(0.5*ntoas))],
//                  int_M1, int_M2, parameters.f,
//                  false, false, parameters.aberration,
//                  delay_aber,  parameters.truefreq, fshift, NULL
//                 ) ;
//         Spline delays_interpolation_ein( tinterp, delay_ein, ninterp );
//             Spline delays_interpolation_shap( tinterp, delay_shap, ninterp );
//                 Spline delays_interpolation_aber( tinterp, delay_aber, ninterp );
//
//     long int firsttis = 0L;
//     long int tintinterp = 0L;
//     long int maxtis = 0L;
//     Spline delays_interpolation( tinterp, delay, ninterp );
//     Spline delay_geom_i (ntisaround);
//
//     for (i = 0; i < nfakeBats ; ++i ) // Go from time of emission to time of arrival and record the delays
//     {
//         while (tinterp[tintinterp] < fake_BATs[i] and tintinterp < ninterp-1) {
//             tintinterp += 1;
//         }
//         if (tintinterp == ninterp -1 or tintinterp == 0 ) {
//             printf("Warning : fake_BATs may be out of range of tinterp in Compute_fake_BATS_and_delays_from_parameters !");
//         };
//         firsttis = max(tintinterp - ntisaround / 2, 0) ;
//
//         maxtis =   min(firsttis + ntisaround, ninterp) - firsttis ;// Need to be defined because the tinterp span may be too narrow.
//
//         if (maxtis != ntisaround ) cout << "Error : maxtis != ntisaround" << maxtis << " " << firsttis << endl;
//
//         Delays_Brut_geometric_local(tinterp + firsttis, fake_BATs[i],  maxtis,
//                  sp + firsttis, si + firsttis, so + firsttis, pospsr, velpsr, parameters.posepoch, fake_obs_ssb, //fake_roemer_ss,
//                  parameters.distance, parameters.distance1,
//                  delay_geom, parameters.truefreq, parameters.kopeikin, parameters.shklovskii, 0.
//                 )  ;
//
// //         for (i = maxtis; i < ntisaround ; i++) // delay_geom has to be defined on ntisaround times. If it is not possible because the tinterp span is too narrow, just extend it.
// //         {
// //             delay_geom[i] = delay_geom[maxtis - 1];
// //         }
//
//         delay_geom_i.ReSpline( tinterp + firsttis, delay_geom ) ; // !!!!!! Will make problems anyway if maxtis < ntisaround
//
//         fake_delay_geom[i] = delay_geom_i(fake_BATs[i], 0,
//                                             ntisaround)  * 86400. ;
//         fake_delay_ein[i] = delays_interpolation_ein( fake_BATs[i], max(0,tintinterp - marginmin),
//                                             min(tintinterp +marginmin, ninterp) )  * 86400. ;
//         fake_delay_shap[i] = delays_interpolation_shap( fake_BATs[i],max(0,tintinterp - marginmin),
//                                             min(tintinterp +marginmin, ninterp) )  * 86400. ;
//         fake_delay_aber[i] = delays_interpolation_aber( fake_BATs[i], max(0,tintinterp - marginmin),
//                                             min(tintinterp +marginmin, ninterp) )   * 86400. ;
//
//         fake_BATs[i] += delays_interpolation( fake_BATs[i],max(0,tintinterp - marginmin),
//                                             min(tintinterp +marginmin, ninterp) ) +
//                         delay_geom_i(fake_BATs[i], 0,
//                                             ntisaround) ;
//        // fake_BATs[i] += parameters.timeshift ;
//
//
//      //   printf("fake_bat 2 : %.19Le \n",  fake_BATs[i]);
//     }
//     free(delay);
//     free(delay_geom);
//     free(delay_ein);
//     free(delay_aber);
//     free(delay_shap);
//
//     return ;
// }
