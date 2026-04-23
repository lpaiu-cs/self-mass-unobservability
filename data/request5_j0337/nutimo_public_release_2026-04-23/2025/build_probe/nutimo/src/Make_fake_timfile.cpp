// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * This particular file contains a code that generates fake times of arrival fom a set of parameters and an tim file. 
 * Times of emission (TOEs) are calculated for each SAT in the tim file using the parameters.
 * The closest integer number of turns of the yo each TOE is then used to calculate fake times of arrival i.e. the ideal time of arrival had there not been any noise or pertubation.
 *
 * Written by Guillaume Voisin 2021 , LUTh, Observatoire de Paris, PSL Research University - CNRS (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 *

 */

#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include "Utilities.h"
#include <stdio.h>
#include <string.h>

#include "Fittriple.h"
#include "IO.h"
#include "Constants.h"



// Needs at least gcc 4.8
// To compile without fittriple for monoproc : g++-4.8 -g -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp
// To run without fittriple : ./MCMCtest.exe fakeparfile fakedatfile mcmout.dat nbw_per_proc maxit chain_freq coeur (moy_freq) (previous_chain)
// To compile with fittriple for monoproc : g++-4.8  -g -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp  AllTheories3Bodies.cpp Delay_brut.cpp Fittriple-compute.cpp Fittriple-init.cpp Fittriple-IO.cpp Fittriple-diagnostics.cpp IO.cpp  Orbital_elements.cpp Spline.cpp  Diagnostics.cpp Parameters.cpp -Bstatic -I/usr/local/src/boost_1_55_0/  -I/usr/share/tempo2/include/ -L./libstatictempo2  -ltempo2 -lsofa
// To compile without fittriple and MPI :  mpic++ -D MPIMODE -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp
// To run : mpiexec -n 4 ./MCMCtest.exe parfile timfile outfile.dat nb_of_walkers_per_proc number_of_ensemble_iterations chain_saving_freq



using namespace std;



int main ( int argc, char *argv[] ) {

// Input/output
    char  parfile[500] ;
    char  datafile[500] ;
    value_type noisestddev;
    char faketimfile[500];
    char fakebats[550];
    char turnfile[550];
    char detailfile[550];


    strcpy(fakebats,"BATs-");
    strcpy(turnfile,"turns-");


    if (argc < 5 )
    {
      printf("Make_fake_timfile.exe parfile original_timfile noisestddev faketimfile\n");
      printf("\n noisestddev : std dev of Gaussian noise. If <0, then uses values of original_timfile multiplied by efac\n");
      return 1;
    }
    strcpy(parfile,argv[1]) ;
    strcpy(datafile, argv[2]);
    sscanf(argv[3], "%Le",&noisestddev) ;
    strcpy(faketimfile, argv[4]);
    strcat(fakebats,faketimfile);
    strcat(turnfile,faketimfile);
    strcpy(detailfile, faketimfile);
    strcat(detailfile, "-details");
    
    int i = 0;
// Random generator init
    // obtain a seed from the system clock:
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine
    normal_distribution<double> gaussian_noise(0.,1.);

    Fittriple fonction(parfile,datafile);

    double * noise;
    value_type * fake_BATs;
    value_type * fake_SATs;
    
    int ntoas = fonction.ntoas;
    fake_BATs = (value_type*) malloc(ntoas * sizeof(value_type));
    fake_SATs = (value_type*) malloc(ntoas * sizeof(value_type));
    noise = (double *) malloc(sizeof(double) * fonction.ntoas);
    
    fonction.Compute_fake_BATs_from_parameters(fake_BATs);


    // Uncomment to get fake bats as well
    Savetxt_L(fakebats, fake_BATs,  static_cast<int>(fonction.ntoas)) ;

    // Save turn numbers
    fonction.Save_turn_numbers(turnfile);

    fonction.Compute_fake_SATs_from_BATs(fake_BATs, fake_SATs);
    
    // Compute noise :
    if (noisestddev < 0) 
        for (i = 0; i < fonction.ntoas; i++) noise[i] = static_cast<double>(fonction.errors[i] * fonction.parameters.efac);
    else
        for (i = 0; i < fonction.ntoas; i++) noise[i] = static_cast<double>(noisestddev);
        
    for (i = 0; i < fonction.ntoas; i++) fake_SATs[i] += noise[i]* gaussian_noise(randomgen) * 1.e-6 / daysec ; // add noise in days
    
    
// SAVE RESULTS
    Write_new_tim_file(datafile, faketimfile, fake_SATs, noise);
    
    FILE * myfile ;
    myfile = fopen(detailfile, "w") ;
    fprintf(myfile, "# Original SATs  |  Turn number  |  Fake BATs  |  Fake Sats \n");
    for(i=0 ; i < ntoas ; ++i)
    {
        fprintf(myfile, "%.19Le     %Ld     %.19Le     %.19Le\n", fonction.Get_sat(i), fonction.turns[i], fake_BATs[i], fake_SATs[i]);
    }
    
    fclose(myfile);

    free(fake_BATs);
    free(fake_SATs);


    return 0;
    }
