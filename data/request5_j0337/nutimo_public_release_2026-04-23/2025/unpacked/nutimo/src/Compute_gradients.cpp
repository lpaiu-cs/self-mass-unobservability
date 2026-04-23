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
 * *

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
    char gradientfile[500];
    double min_dlnp, max_dlnp, eps;
    
    if (argc < 7 )
    {
      printf("Compute_gradients.exe parfile timfile gradientfile epsilon min_dlnp max_dlnp\n");
      printf(" Return the gradients of each residual with respect to each fitted parameters normalised to its scale\n");
      printf(" epsilon : fraction of parameter scale to use as initial guess in finite differentiation\n");
      printf(" min_dlnp : minimum (absolute) variation of lnposterior (to small may lead to rounding errors).\n");
      printf(" max_dlnp : maximum aboslution variation of lnposterior\n");
      return 1;
    }
    strcpy(parfile,argv[1]) ;
    strcpy(datafile, argv[2]);
    strcpy(gradientfile, argv[3]);
    sscanf(argv[4], "%le",&eps) ;
    sscanf(argv[5], "%le",&min_dlnp) ;
    sscanf(argv[6], "%le",&max_dlnp) ;
    
    printf("parfile=%s \ntimfile=%s \ngradientfile=%s \nepsilon=%.5e \nmin_dlnp=%.5e \n max_dlnp=%.5e\n\n", parfile, datafile, gradientfile, eps, min_dlnp, max_dlnp);
    
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    
    // Check  variables 
    double dxcheck =0;
    double dtcheck = 0;
    double * res0;
    double dtcheckmax, meancheck, stdcheck;
    int icheckmax;
    
    
    
    Fittriple fonction(parfile,datafile);
    
    fonction.remove_mean=false; // Prevents marginalisation over mean. 

    int npar = fonction.Get_nfitted_params();
    int ntoas = fonction.ntoas;
    double ** dt ;
    vector<double> params(npar, 0.);
//     double params[npar]; // parameter relative shifts
    double scales[npar]; // parameter scales. 
    for (j = 0; j < npar ; j++) 
    {
//         params[j] = 0.;
        scales[j] = fonction.parameters.Parameter_scale(fonction.parameters.fitted_parameters[j]);
    }
    dt = (double**) malloc(sizeof(double*) * ntoas);
    for (i = 0; i < ntoas ; i++) dt[i] = (double*) malloc(sizeof(double) * npar);
    
    
    double dx ;
    double lnp2, lnp1;
    bool accept ;
    
    fonction.Set_fitted_parameter_relativeshifts(params);
    lnp2 = fonction.Compute_lnposterior(0) * ntoas;
    printf(" Initial lnp : %.10e\n", lnp2);
    // For check
    res0 = (double*) malloc(sizeof(double) * ntoas);
    for (i = 0; i < ntoas ; i++) res0[i] = fonction.residuals[i];
     
    for (j = 0; j < npar ; j++)
    {
        accept= false;
        k = 0;
        dx = eps; 
        fonction.motion_changed= true;
        printf("Parameter %d : %s\n", j, fonction.parameters.Get_parameter_name(fonction.parameters.fitted_parameters[j])) ;
        while (accept == false and k < 5)
        {
            k += 1;
            for (l = 0; l < npar ; l++) params[l] = 0.; // Not necessary in principle but failsafe..
            params[j] = 0.5 * dx;
//             lnp1 = fonction(params);
            fonction.Set_fitted_parameter_relativeshifts(params);
            if (fonction.parameters.motion_changed == false ) fonction.motion_changed= false;
            lnp1 = fonction.Compute_lnposterior(0) * ntoas;
            for (i = 0; i < ntoas ; i++) dt[i][j] = fonction.residuals[i];  
            params[j] = -params[j];
//             lnp2 = fonction(params);
            fonction.Set_fitted_parameter_relativeshifts(params);
            lnp2 = fonction.Compute_lnposterior(0) * ntoas;
            for (i = 0; i < ntoas ; i++) dt[i][j] -= fonction.residuals[i];
            if (abs(lnp2 - lnp1) > max_dlnp ) 
                dx /=10;
            else if (abs(lnp2 - lnp1) < min_dlnp  ) 
                dx *= 10;
            else 
                accept = true;
            printf("  %d %.10e\n", k, abs(lnp2 - lnp1));
        }
        if (accept == false ) printf("Warning : Delta chi2 not acceptable. %.10e\n", abs(lnp2 - lnp1));
        params[j] = 0.;
        for (i = 0; i < ntoas ; i++) dt[i][j] /= dx;
        
        
    // Perform checks
        for (l= 0; l < npar ; l++) params[l] = 0.; // Not necessary in principle but failsafe..
        // check short
        dxcheck = 2.;//0.25 * dx;
        params[j] = dxcheck;
        fonction.Set_fitted_parameter_relativeshifts(params);
        lnp2 = fonction.Compute_lnposterior(0) * ntoas;
        dtcheckmax = 0.;
        meancheck = 0.;
        stdcheck = 0.;
        for (i = 0; i < ntoas ; i++) 
        {
            dtcheck = fonction.residuals[i] - res0[i];
            dtcheck -= dt[i][j]*dxcheck;
            if (abs(dtcheck) > abs(dtcheckmax) )
            {
                dtcheckmax = abs(dtcheck);
                icheckmax = i;
            }
            meancheck += abs(dtcheck);
            stdcheck += dtcheck*dtcheck;
        }
        meancheck /= ntoas;
        stdcheck = sqrt(stdcheck/ntoas - meancheck*meancheck);
        printf(" Check short : max = (%.5e, %d) mean = %.5e std = %.5e\n", dtcheckmax, icheckmax, meancheck, stdcheck);
        // Check long
//         dxcheck = - dx;
//         params[j] = dxcheck;
//         fonction.Set_fitted_parameter_relativeshifts(params);
//         lnp2 = fonction.Compute_lnposterior(0) * ntoas;
//         dtcheckmax = 0.;
//         meancheck = 0.;
//         stdcheck = 0.;
//         for (i = 0; i < ntoas ; i++) 
//         {
//             dtcheck = fonction.residuals[i] - res0[i];
//             dtcheck -= dt[i][j]*dxcheck;
//             if (abs(dtcheck) > abs(dtcheckmax) )
//             {
//                 dtcheckmax = dtcheck;
//                 icheckmax = i;
//             }
//             meancheck += abs(dtcheck);
//             stdcheck += dtcheck*dtcheck;
//         }
//         meancheck /= ntoas;
//         stdcheck = sqrt(stdcheck/ntoas - meancheck*meancheck);
//         printf(" Check long : max = (%.5e, %d) mean = %.5e std = %.5e\n", dtcheckmax, icheckmax, meancheck, stdcheck);
    }
    
    Savetxt(gradientfile, dt, ntoas, npar );

    for (i = 0; i < ntoas ; i++) free(dt[i]);
    free(dt);
    free(res0);
    
    return 0;
    }
