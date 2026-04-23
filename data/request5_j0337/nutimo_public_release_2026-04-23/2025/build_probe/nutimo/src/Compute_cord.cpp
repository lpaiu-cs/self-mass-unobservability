// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 * 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 * 
 * * 
 * This program samples a distribution using the algorithm described in Goodman & Weare 2010 (GW10), Communitcations in Applied mathematics and Computational Sciences, Vol5, Issue 1
 * The parallelization scheme is from : Foreman-Mackey et al., 2013, Pulications of the Astronomical Society of the Pacific, Vol 125, Issue 925
 * 
 */

#ifdef MPIMODE
    #include "mpi.h"
#endif
#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include "Utilities.h"
#include "acor.h"
#include <stdio.h>
#include <string.h>

#include "Fittriple.h"



// Needs at least gcc 4.8
// To compile without fittriple for monoproc : g++-4.8 -g -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp
// To run without fittriple : ./MCMCtest.exe fakeparfile fakedatfile mcmout.dat nbw_per_proc maxit chain_freq coeur (moy_freq) (previous_chain)
// To compile with fittriple for monoproc : g++-4.8  -g -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp  AllTheories3Bodies.cpp Delay_brut.cpp Fittriple-compute.cpp Fittriple-init.cpp Fittriple-IO.cpp Fittriple-diagnostics.cpp IO.cpp  Orbital_elements.cpp Spline.cpp  Diagnostics.cpp Parameters.cpp -Bstatic -I/usr/local/src/boost_1_55_0/  -I/usr/share/tempo2/include/ -L./libstatictempo2  -ltempo2 -lsofa
// To compile without fittriple and MPI :  mpic++ -D MPIMODE -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp
// To run : mpiexec -n 4 ./MCMCtest.exe parfile timfile outfile.dat nb_of_walkers_per_proc number_of_ensemble_iterations chain_saving_freq



using namespace std;


class fctgaussienne
{
    int ndim ;
public:
    int verbose=0;
    
    fctgaussienne(int number_of_dim){ndim =number_of_dim;};
    double operator()(double* pars)
    {
        double res = 0.;
        double inter = 0.;
//          res -= pow(pars[0] - 1000.*pars[1],2);
//          res -= pow(pars[0] + pars[1],2);
        for (int i = 0; i < ndim ; i++) 
        {
            inter = pow(pars[i] - 0.2 * static_cast<double>(i),2);
            res -= inter / static_cast<double>(i+1) ;
        }
        return res;
    };
    void Print_MCMC()
    {
        printf("Fonction Gaussienne in %d dimensions! \n", ndim);
    };
};



class fctmultigaussienne
{
    int ndim ;
    int npeaksperdim;
    double stdmainpeak=1.;
    double stdpeaksec;
    double peakstep ;
public:
    int verbose=0;
    
    fctmultigaussienne(int number_of_dim, int half_number_of_peaks_per_dim)
    {
        ndim =number_of_dim; 
        npeaksperdim=half_number_of_peaks_per_dim;
        stdpeaksec = stdmainpeak / (4.* npeaksperdim);
        peakstep = stdmainpeak;
    };
    
    double operator()(double* pars)
    {
        double res = 0.;
        double inter = 0.;
        int k ;
        
//         printf("pars %.5e %.5e \n", pars[0], pars[1]);
        
        for (int i = 0; i < ndim ; i++) 
        {
            for (int j = -npeaksperdim ; j <= npeaksperdim ; j++)
            {
                if (j != 0) 
                {
                    inter = 0.;
                    for (k = 0; k < ndim ; k++) 
                    {
                        if (k == i ) 
                            inter += -pow((pars[i] - peakstep*j)/stdpeaksec,2)*0.5;
                        else 
                            inter += -pow((pars[k])/stdpeaksec,2) * 0.5;
//                         printf("i %d j %d k %d inter %.5e \n", i, j, k, inter);
                    }
                    res += exp(inter);
//                     printf("i %d j  %d exp(inter) %.5e \n", i, j, exp(inter));
                }
            }
            inter = 0.;
            for (k = 0; k < ndim ; k++) inter += -pow((pars[k])/stdmainpeak,2) *0.5;
            //res += exp(inter);
        }
        return log(res);
    };
    
    
    void Print_MCMC()
    {
        printf("Fonction Multi Gaussienne in %d dimensions! \n", ndim);
        printf("npeaksperdim %d stdmainpeak %.5e stdpeaksec %.5e peakstep %.5e \n", npeaksperdim, stdmainpeak, stdpeaksec, peakstep);
    };
};


class fctFittriple : public Fittriple // Wraps Fittriple object with a prior and a different parameter set
{
    const value_type secyr = 31557600.0; // Number of seconds per year
    
    // Position from Ransom et al. 2014
    const value_type posransom_rad = 9.5002822899862595762e-01;
    const value_type posransom_dec = 3.0114118619832732786e-01; 
    const value_type posprior = 2;
    
    // Radial velocity from Kaplan et al. 2014
    const value_type rvkaplan = 29.7e+03; //m/s
    const value_type rverrorkaplan = 0.9e+03; // m/s
    const value_type rvprior = 2;
    
    // Photometric distance from Kaplan et al. 2014
    const value_type dkaplan = 4078.5496710122447 ; // lyr
    const value_type derrorkaplan = 250.98767206229201; // lyr
    const value_type dprior = 2;
    
public :
    
    void Print_MCMC()
    {
        printf("Position prior is +/- %.2Le mas.\n", posprior);
        printf("Radial velocity prior is %.2Le +/- %.2Le km/s.\n", rvkaplan/1000., rverrorkaplan*rvprior/1000.);
        printf("Distance prior is %.3Le +/- %.3Le ly.\n", dkaplan, dprior*derrorkaplan);
        double chi2 = Compute_lnposterior(0);
        printf("Initial reduced chi2 is %.8f, and chi2 = %.8f\n", chi2, chi2*ntoas );
        printf("\n Parameter map : [");
        for (int i=0; i < Get_nfitted_params() ; i++) printf("%i, ", Get_list_of_fitted_parameters(i));
        printf("]\n");
    };
    
    double operator()(double * relativeshift)
    {
        vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
        Set_fitted_parameter_relativeshifts(vectrelativeshift);
        
        if ( (parameters.RA < posransom_rad - posprior *radmasdeg )   or  (parameters.RA > posransom_rad + posprior*radmasdeg) )
        {
          //  printf("prior ra %.15Le %.15Le %.15Le %.15Le\n", parameters.RA, posransom_rad, posransom_rad - posprior *radmasdeg, posransom_rad + posprior*radmasdeg);
            return  -DBL_MAX ;
        }
        else if ( (parameters.DEC < posransom_dec - posprior* radmasdeg) or  (parameters.DEC > posransom_dec + posprior*radmasdeg) )
            {
            //printf("prior rdec\n");
            return  -DBL_MAX ;
            }
        else if ( (parameters.distance < dkaplan - dprior* derrorkaplan) or  (parameters.distance > dkaplan + dprior* derrorkaplan) )
            {
  //          printf("prior d\n");
            return  -DBL_MAX ;
            }
        else if ( (parameters.distance1*parameters.distance*clight < rvkaplan - rvprior* rverrorkaplan) or  (parameters.distance1*parameters.distance*clight > rvkaplan + rvprior* rverrorkaplan) )
        {
//            printf("prior rv %.15Le %.15Le %.15Le %.15Le\n", parameters.distance1*parameters.distance*clight , rvkaplan,  rvkaplan - rvprior* rverrorkaplan,rvkaplan + rvprior* rverrorkaplan);
            return  -DBL_MAX ;
        }
        else
            return Compute_lnposterior(0) * ntoas; // Return the actual log(proba) = -chi2 , not the reduced one
        
    };
};





class fctFittriple_gaussprior : public Fittriple // Wraps Fittriple object with a Gaussian prior 
{
    const value_type secyr = 31557600.0; // Number of seconds per year
    
    // Position from Ransom et al. 2014 (radians)
    const value_type posransom_rad = 9.5002822899862595762e-01;
    const value_type posransom_dec = 3.0114118619832732786e-01; 
    const value_type raerror = 9.45386657846825e-09;
    const value_type decerror = 9.69627362219072e-09;
    const value_type posprior = 2;
    
    // Radial velocity from Kaplan et al. 2014
    const value_type rvkaplan = 29.7e+03; //m/s
    const value_type rverrorkaplan = 0.9e+03; // m/s
    const value_type rvprior = 2;
    
    // Photometric distance from Kaplan et al. 2014
    const value_type dkaplan = 4078.5496710122447 ; // lyr
    const value_type derrorkaplan = 250.98767206229201; // lyr
    const value_type dprior = 2;
    
public :
    int verbose=0;
    double beta =1.; // 1/Temperature
    
    double lastlnlikelyhood;
    
    
    void Print_MCMC()
    {
        printf("Beta = 1/temperature = %.5e \n",beta);
        printf("Gaussian priors at 1sigma :\n");
        printf(" Position prior is RA +/- %.2Le mas.\n", posprior);
        printf(" Radial velocity prior is %.2Le +/- %.2Le km/s.\n", rvkaplan/1000., rverrorkaplan*rvprior/1000.);
        printf(" Distance prior is %.3Le +/- %.3Le ly.\n", dkaplan, dprior*derrorkaplan);
        double chi2 = Compute_lnposterior(0);
        printf("Initial reduced chi2 is %.8f, and chi2 = %.8f\n", chi2, chi2*ntoas );
        printf("\nParameter map : [");
        for (int i=0; i < Get_nfitted_params() ; i++) printf("%i, ", Get_list_of_fitted_parameters(i));
        printf("]\n");
    };
    
    double operator()(double * relativeshift)
    {
        vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
        Set_fitted_parameter_relativeshifts(vectrelativeshift);
        double ragauss = -0.5*pow((parameters.RA - posransom_rad)/(raerror*posprior),2);
        double decgauss = -0.5*pow((parameters.DEC - posransom_dec)/(decerror*posprior),2);
        double rvgauss = -0.5 * pow((parameters.distance1*parameters.distance*clight - rvkaplan)/(rverrorkaplan*rvprior),2);
        double dgauss = -0.5 * pow((parameters.distance - dkaplan)/(derrorkaplan*dprior),2);
        
        lastlnlikelyhood = Compute_lnposterior(0) * ntoas * beta;
        
        if (verbose == 1) printf("lnlikelyhood %.3e ragauss %.3e decgauss %.3e rvgauss %.3e dgauss %.3e\n", lastlnlikelyhood , ragauss , decgauss , rvgauss , dgauss);
        
        return (lastlnlikelyhood + ragauss + decgauss + rvgauss + dgauss) ;
    };
};


class gw10_distribution // the g distribution in gw10
{
    double a;
    double cte;
    double Ng;
    
public:
    
    gw10_distribution(double parameter_a){a = parameter_a; Ng = 2. * (sqrt(a) - sqrt(1./a)); cte = - 2. / (sqrt(a) * Ng );};
    double operator()(double uniform_random_number_01){return 0.25 * pow(Ng *(uniform_random_number_01 - cte), 2 );}; // uniform_random_number_01 is a random number between [0,1[ drawn from a uniform distribution
};


double test_gw10_distribution(int nb_tirages, double a) // test the distribution function defined by gw10_distribution by comparing the ratio of (number events in [a/2, a]) / (number events in [1/a, a/2[) with the theoretical value. there are nb_tirages events.
{
  unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
  mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine
  uniform_real_distribution<double> random01(0.0,1.0);
  gw10_distribution distro(a);
  double s =2.;
  double gmin=0.;
  double gmax = 0.;
  double valth=( sqrt(a) - sqrt(a/s) ) / ( sqrt(a/s) - sqrt(1./a) ) ;
  
  for(int i = 0; i < nb_tirages; i++)
  {
      if ( distro(random01(randomgen)) < a/s) 
          gmin += 1. ;
      else
          gmax += 1. ; 
  }
  printf("\n test_gw10_distribution : gmax / gmin = %f . Différence ac valeur théorique (%f) : %f \n\n", gmax / gmin, valth, gmax / gmin - valth );
  return gmax / gmin - valth;
}



void Save_chain(char * filename, double ** chain, int ndim, int nb_walkers, int chain_size, int chain_freq)
{
    int i , j;
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    fprintf(myfile, "#  %i  %i  %i %i\n",ndim, nb_walkers, chain_size, chain_freq);
    for(i=0 ; i < chain_size ; ++i) {
        for (j = 0; j < ndim+1 ; ++j ) {
             fprintf(myfile, "%.15e    ", chain[i][j] );
         }
         fprintf(myfile, "\n") ; 
     }
     fclose(myfile);
}


void Load_chain(char * filename, double ** chain, int ndim, int nb_walkers, int & chain_size, int chain_freq)
// Load a chain saved with "Save_chain". Allocate chain with chain_size + size of the chain in the file. Update chain_size.
// Check that the characteristics of the loaded chain are the same as those of chain (walkers, frequency, dimension). Issue a Warning or an error if it is not the case. 
{
    int i , j;
    FILE * myfile ;
    myfile = fopen(filename, "r") ;
    int fndim=0; 
    int fchain_freq =0;
    int fnb_walkers =0;
    int fchain_size =0;
    fscanf(myfile, "#  %i  %i  %i %i\n",&fndim, &fnb_walkers, &fchain_size, &fchain_freq);
    if (ndim != fndim)
    {
        printf("\n Error ! Chain in file %s and desired chain have different number of dimensions : %i against %i \n\n", filename, fndim, ndim);
        return;
    }
    chain = (double**) malloc((chain_size + fchain_size) * sizeof(double*) );
    for (i =0  ; i < (chain_size + fchain_size) ; i++)
    {
        chain[i] = (double*) malloc((ndim+1)*sizeof(double));
        for (j= 0 ; j < ndim + 1 ; j++) chain[i][j] = 0.;
    }
    
    for(i=0 ; i < fchain_size ; ++i) {
        for (j = 0; j < fndim ; ++j ) {
             fscanf(myfile, "%le    ", &(chain[i][j]) );
         }
         fscanf(myfile, "%le    \n", &(chain[i][fndim]) ); 
     }
     fclose(myfile);
     
     if (nb_walkers != fnb_walkers)
     {
         printf("\n Warning ! nb_walkers = %i != fnb_walkers = %i \n\n", nb_walkers, fnb_walkers);
     }
     if (chain_freq != fchain_freq)
     {
         printf("\n Warning ! chain_freq = %i != fchain_freq = %i \n\n", chain_freq, fchain_freq);
     }
     chain_size = chain_size + fchain_size;
     
     return;
}

int Init_from_prev_chain(char * prev_chain_file, double * walkers1d, const  int ndim, const int nb_walkers)
{
    
    int i , j;
    FILE * myfile ;
    myfile = fopen(prev_chain_file, "r") ;
    int fndim=0; 
    int fchain_freq =0;
    int fnb_walkers =0;
    int fchain_size =0;
    double ** chain;
    
    fscanf(myfile, "#  %i  %i  %i %i\n",&fndim, &fnb_walkers, &fchain_size, &fchain_freq);
    if (ndim != fndim)
    {
        printf("\n Error ! Chain in file %s and desired chain have different number of dimensions : %i against %i \n\n", prev_chain_file, fndim, ndim);
        return 1 ;
    }
    
    if ( fchain_size < fnb_walkers ) 
    {
        printf("\n Error ! Chain in file %s is too small for initialization : length is %i against %i walkers needed \n\n", prev_chain_file, fchain_size, nb_walkers);
        return 1 ;
    }
    
    chain = (double**) malloc((fchain_size) * sizeof(double*) );
    for (i =0  ; i < (fchain_size) ; i++)
    {
        chain[i] = (double*) malloc((ndim+1)*sizeof(double));
        for (j= 0 ; j < ndim + 1 ; j++) chain[i][j] = 0.;
    }
    
    for(i=0 ; i < fchain_size ; ++i) {
        for (j = 0; j < fndim ; ++j ) {
             fscanf(myfile, "%le    ", &(chain[i][j]) );
         }
         fscanf(myfile, "%le    \n", &(chain[i][fndim]) ); 
     }
     fclose(myfile);
     
     if (nb_walkers != fnb_walkers)
     {
         printf("\n Warning ! nb_walkers = %i != fnb_walkers = %i \n\n", nb_walkers, fnb_walkers);
     }

     
     for (i = 0 ; i < nb_walkers ; i ++)
     {
         for (j = 0 ; j < ndim ; j++) walkers1d[i*ndim + j ] = chain[fchain_size - 1  - i][j] ;
     }
     
     for (i =0  ; i < (fchain_size) ; i++) free(chain[i]);
     free(chain);
     
     return 0 ;
}


int main ( int argc, char *argv[] ) {
 
// MPI variable, by default for 1 proc. (no mpi)
    int code=0;
    int nbproc=1;
    int mpirank =0;
    int masterproc = 0; // process doing stuff alone
    int coeursparnoeud = 10000; // Used to do parallelism per node and not per thread. 

//
    int npts = 0;
    int i,j,k;
    int ndim =0;
    double * chi2;
    double * chi2_local;
    double width = 2.; // in units of sigma, on one one side of params0 ( so total is twice) 

// MPI initialization

#ifdef MPIMODE
    code = MPI_Init(&argc, &argv) ;
    code = MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    code = MPI_Comm_rank(MPI_COMM_WORLD , &mpirank);
    
    printf("mpirank %i \n\n", mpirank);
#endif
    
    
// Input/output
    char  fileout[500] ;        // where the chain is saved
    char  parfile[500] ;        // passed to the function "fonction " object " (only use)
    char  datafile[500] ;       // passed to the function "fonction " object " (only use)
    char  cordfile[500] ; 
    
    strcpy(parfile,argv[1]) ;
    strcpy(datafile, argv[2]);
    strcpy(cordfile, argv[3]);
    sscanf(argv[4], "%le",&width) ;
    sscanf(argv[5], "%i",&npts) ;
    strcpy(fileout, argv[6]);
  
    coeursparnoeud = nbproc;
 
    int npercore = npts/nbproc;
    printf("npercore %i\n",npercore);
    
 fctFittriple_gaussprior fonction = fctFittriple_gaussprior();
 fonction.verbose=0;
 
   for (i = 0 ; i < coeursparnoeud ; i++) // initialization simultaneous on each node but not on each thread because of I/O with files
   {
#ifdef MPIMODE    
       code = MPI_Barrier(MPI_COMM_WORLD);
#endif
       if (i == mpirank%coeursparnoeud) 
       {
           fonction.Reconstruct(parfile, datafile);
           printf("\n ---- > passed %i \n", mpirank);
       }
   }
   fonction.Set_tracker_off();
 
 ndim = fonction.Get_nfitted_params();   // set ndim 
#ifdef MPIMODE    
   code = MPI_Barrier ( MPI_COMM_WORLD); // Not really necessary but synchronizes the prints
#endif 
   
   // Initiaze parameters
   double params[ndim];
   double params0[ndim];
   double step = 2.* width/ npts; // in unit of sigmas 
   double sigmas[ndim]; // std dev of each parameter
   double corddirection[ndim]; // direction of the cord accross the center
   double cordnorm =0.;
   
   if (mpirank==0)
   {
    Loadtxt("meanparams.txt", params0, ndim);
    Loadtxt("sigmaparams.txt", sigmas, ndim);
    Loadtxt(cordfile, corddirection, ndim);
   

   printf("Means:\n");
   Print_table(params0, ndim);
   
    printf("Std dev:\n");
   Print_table(sigmas, ndim);
   
    printf("Cord direction:\n");
   Print_table(corddirection, ndim);

   for (j = 0 ; j < ndim ; j++)  cordnorm += pow(corddirection[j],2);
   cordnorm = sqrt(cordnorm);
//    for (j = 0 ; j < ndim ; j++)  corddirection[j] /= cordnorm;
   for (j = 0 ; j < ndim ; j++)  params0[j] -=  corddirection[j] ;//* sigmas[j] * width;
   
   printf("Cord direction normalized:\n");
   Print_table(corddirection, ndim);
   
   printf("\nWidth %.2e Step %.2e npts %i\n", width, step, npts);
  }
#ifdef MPIMODE 
    code = MPI_Bcast(params0, ndim, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
    code = MPI_Bcast(sigmas, ndim, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
    code = MPI_Bcast(corddirection, ndim, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
#endif
// Compute
    chi2 = (double*) malloc(sizeof(double) * npercore * nbproc);
    chi2_local= (double*) malloc(sizeof(double) * npercore );
    
    for (i=npercore*mpirank ; i < npercore*(mpirank+1) ; i++)
    {
        if (mpirank==0) printf("ne means %i/%i \n", i-npercore*mpirank, npercore);
        for (j =0; j < ndim; j++)
        {
            params[j] = params0[j] + i * step*corddirection[j];//*sigmas[j];
        }
//         printf("means %i \n", i);
//         Print_table(params, ndim);
        chi2_local[i-npercore*mpirank] = fonction(params); 
    }
    
    code = MPI_Allgather(chi2_local , npercore, MPI_DOUBLE, chi2, npercore, MPI_DOUBLE, MPI_COMM_WORLD);
    
    if (mpirank==0)
    {
        Savetxt(fileout, chi2, npercore*nbproc);
    };
    
    free(chi2);
    free(chi2_local);
    
#ifdef MPIMODE
    code = MPI_Finalize();
#endif
    
  return 0;
    }
