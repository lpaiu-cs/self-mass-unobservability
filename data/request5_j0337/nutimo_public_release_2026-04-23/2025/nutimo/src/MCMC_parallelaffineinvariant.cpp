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
#include "Fittriple.h"

#include "MCMC_lnposterior_functions.h"
#include "Utilities_MCMC_affinv.h"
#include "Utilities_MPI.h"

// Needs at least gcc 4.8
// To compile without fittriple for monoproc : g++-4.8 -g -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp
// To run without fittriple : ./MCMCtest.exe fakeparfile fakedatfile mcmout.dat nbw_per_proc maxit chain_freq coeur (moy_freq) (previous_chain)
// To compile with fittriple for monoproc : g++-4.8  -g -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp  AllTheories3Bodies.cpp Delay_brut.cpp Fittriple-compute.cpp Fittriple-init.cpp Fittriple-IO.cpp Fittriple-diagnostics.cpp IO.cpp  Orbital_elements.cpp Spline.cpp  Diagnostics.cpp Parameters.cpp -Bstatic -I/usr/local/src/boost_1_55_0/  -I/usr/share/tempo2/include/ -L./libstatictempo2  -ltempo2 -lsofa
// To compile without fittriple and MPI :  mpic++ -D MPIMODE -std=gnu++11 -o MCMCtest.exe MCMC_parallelaffineinvariant.cpp Utilities.cpp acc.cpp
// To run : mpiexec -n 4 ./MCMCtest.exe parfile timfile outfile.dat nb_of_walkers_per_proc number_of_ensemble_iterations chain_saving_freq



using namespace std;

//
// class fctgaussienne
// {
//     int ndim ;
// public:
//     fctgaussienne(int number_of_dim){ndim =number_of_dim;};
//     double operator()(double* pars)
//     {
//         double res = 0.;
//         double inter = 0.;
//          res -= pow(pars[0] - 1000.*pars[1],2);
//          res -= pow(pars[0] + pars[1],2);
//         for (int i = 2; i < ndim ; i++)
//         {
//             inter = pow(pars[i] - 0.2 * static_cast<double>(i),2);
//             res -= inter / static_cast<double>(i+1) ;
//         }
//         return res;
//     };
//     void Print_MCMC()
//     {
//     };
// };
//
//
// class fctFittriple : public Fittriple // Wraps Fittriple object with a prior and a different parameter set
// {
//     const value_type secyr = 31557600.0; // Number of seconds per year
//
//     // Position from Ransom et al. 2014
//     const value_type posransom_rad = 9.5002822899862595762e-01;
//     const value_type posransom_dec = 3.0114118619832732786e-01;
//     const value_type posprior = 2;
//
//     // Radial velocity from Kaplan et al. 2014
//     const value_type rvkaplan = 29.7e+03; //m/s
//     const value_type rverrorkaplan = 0.9e+03; // m/s
//     const value_type rvprior = 2;
//
//     // Photometric distance from Kaplan et al. 2014
//     const value_type dkaplan = 4078.5496710122447 ; // lyr
//     const value_type derrorkaplan = 250.98767206229201; // lyr
//     const value_type dprior = 2;
//
// public :
//
//     void Print_MCMC()
//     {
//         printf("Position prior is +/- %.2Le mas.\n", posprior);
//         printf("Radial velocity prior is %.2Le +/- %.2Le km/s.\n", rvkaplan/1000., rverrorkaplan*rvprior/1000.);
//         printf("Distance prior is %.3Le +/- %.3Le ly.\n", dkaplan, dprior*derrorkaplan);
//         double chi2 = Compute_lnposterior(0);
//         printf("Initial reduced chi2 is %.8f, and chi2 = %.8f\n", chi2, chi2*ntoas );
//         printf("\n Parameter map : [");
//         for (int i=0; i < Get_nfitted_params() ; i++) printf("%i, ", Get_list_of_fitted_parameters(i));
//         printf("]\n");
//     };
//
//     double operator()(double * relativeshift)
//     {
//         vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
//         Set_fitted_parameter_relativeshifts(vectrelativeshift);
//
//         if ( (parameters.RA < posransom_rad - posprior *radmasdeg )   or  (parameters.RA > posransom_rad + posprior*radmasdeg) )
//         {
//           //  printf("prior ra %.15Le %.15Le %.15Le %.15Le\n", parameters.RA, posransom_rad, posransom_rad - posprior *radmasdeg, posransom_rad + posprior*radmasdeg);
//             return  -DBL_MAX ;
//         }
//         else if ( (parameters.DEC < posransom_dec - posprior* radmasdeg) or  (parameters.DEC > posransom_dec + posprior*radmasdeg) )
//             {
//             //printf("prior rdec\n");
//             return  -DBL_MAX ;
//             }
//         else if ( (parameters.distance < dkaplan - dprior* derrorkaplan) or  (parameters.distance > dkaplan + dprior* derrorkaplan) )
//             {
//   //          printf("prior d\n");
//             return  -DBL_MAX ;
//             }
//         else if ( (parameters.distance1*radmasdeg*parameters.distance*clight < rvkaplan - rvprior* rverrorkaplan) or  (parameters.distance1*radmasdeg*parameters.distance*clight > rvkaplan + rvprior* rverrorkaplan) )
//         {
// //            printf("prior rv %.15Le %.15Le %.15Le %.15Le\n", parameters.distance1*parameters.distance*clight , rvkaplan,  rvkaplan - rvprior* rverrorkaplan,rvkaplan + rvprior* rverrorkaplan);
//             return  -DBL_MAX ;
//         }
//         else
//             return Compute_lnposterior(0) * ntoas; // Return the actual log(proba) = -chi2 , not the reduced one
//
//     };
// };
//
//
//
//
//
// class fctFittriple_gaussprior : public Fittriple // Wraps Fittriple object with a Gaussian prior
// {
//     const value_type secyr = 31557600.0; // Number of seconds per year
//     /*
//     // Position from Ransom et al. 2014 (radians)
//     const value_type posransom_rad = 9.5002822899862595762e-01;
//     const value_type posransom_dec = 3.0114118619832732786e-01;
//     const value_type raerror = 9.45386657846825e-09;
//     const value_type decerror = 9.69627362219072e-09;
//     const value_type posprior = 2;
//
//     // Radial velocity from Kaplan et al. 2014
//     const value_type rvkaplan = 29.7e+03; //m/s
//     const value_type rverrorkaplan = 0.9e+03; // m/s
//     const value_type rvprior = 2;
//
//     // Photometric distance from Kaplan et al. 2014
//     const value_type dkaplan = 4078.5496710122447 ; // lyr
//     const value_type derrorkaplan = 250.98767206229201; // lyr
//     const value_type dprior = 2;*/
//
//     // Position from GAIA DR2 (radians) Attention : POSEPOCH = MJD2015.5
//     const value_type ramean = 9.500283096783616e-01;
//     const value_type decmean = 3.011411347204199-01;
//     const value_type raerror = 8.945972920351486e-10;
//     const value_type decerror = 9.302791001352137e-10;
//     const value_type posprior = 2;
//
//     // Tangential proper motion (mas/yr) from Gaia DR2. Warning : POSEPOCH = MJD2015.5
//     const value_type pmramean = 4.81377140723352;
//     const value_type pmdecmean = -4.42182046193475;
//     const value_type pmraerror = 0.498024673422251;
//     const value_type pmdecerror = 0.42770795176738;
//     const value_type tpmprior = 2;
//
//     // Radial velocity from Kaplan et al. 2014
//     const value_type rvkaplan = 29.7e+03; //m/s
//     const value_type rverrorkaplan = 0.9e+03; // m/s
//     const value_type rvprior = 1;
//
//     // Distance from Gaia DR2. Attention : POSEPOCH = MJD2015.5
//     const value_type dmean = 4604.639595490277 ; // lyr
//     const value_type derror = 1558.9119367898795; // lyr
//     const value_type dprior = 1;
//
// public :
//     int verbose=0;
//     double beta =1.; // 1/Temperature
//
//     double lastlnlikelyhood;
//
//
//     void Print_MCMC()
//     {
//         printf("Beta = 1/temperature = %.5e \n",beta);
//         printf("Gaussian priors at 1sigma :\n");
//         printf(" Position prior is RA +/- %.2Le mas.\n", posprior);
//         printf(" Radial velocity prior is %.2Le +/- %.2Le km/s.\n", rvkaplan/1000., rverrorkaplan*rvprior/1000.);
//         printf(" Distance prior is %.3Le +/- %.3Le ly.\n", dmean, dprior*derror);
//         double chi2 = Compute_lnposterior(0);
//         printf("Initial reduced chi2 is %.8f, and chi2 = %.8f\n", chi2, chi2*ntoas );
//         printf("\nParameter map : [");
//         for (int i=0; i < Get_nfitted_params() ; i++) printf("%i, ", Get_list_of_fitted_parameters(i));
//         printf("]\n");
//     };
//
//     double operator()(double * relativeshift)
//     {
//         vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
//         Set_fitted_parameter_relativeshifts(vectrelativeshift);
//         double ragauss = -0.5*pow((parameters.RA - ramean)/(raerror*posprior),2);
//         double decgauss = -0.5*pow((parameters.DEC - decmean)/(decerror*posprior),2);
//         double rvgauss = -0.5 * pow((parameters.distance1*radmasdeg* parameters.distance*clight - rvkaplan)/(rverrorkaplan*rvprior),2);
//         double dgauss = -0.5 * pow((parameters.distance - dmean)/(derror*dprior),2);
//         double pmragauss = -0.5 * pow((parameters.RA1 - pmramean)/(pmraerror*tpmprior),2);
//         double pmdecgauss = -0.5 * pow((parameters.DEC1 - pmdecmean)/(pmdecerror*tpmprior),2);
//
//         lastlnlikelyhood = Compute_lnposterior(0) * ntoas * beta;
//
//         if (verbose == 1) printf("lnlikelyhood %.3e ragauss %.3e decgauss %.3e rvgauss %.3e dgauss %.3e pmragauss %.3e pmdecgauss %.3e\n", lastlnlikelyhood , ragauss , decgauss , rvgauss , dgauss, pmragauss, pmdecgauss);
//
//         return (lastlnlikelyhood + ragauss + decgauss + rvgauss + dgauss + pmragauss + pmdecgauss) ;
//     };
// };
//
//
// class gw10_distribution // the g distribution in gw10
// {
//     double a;
//     double cte;
//     double Ng;
//
// public:
//
//     gw10_distribution(double parameter_a){a = parameter_a; Ng = 2. * (sqrt(a) - sqrt(1./a)); cte = - 2. / (sqrt(a) * Ng );};
//     double operator()(double uniform_random_number_01){return 0.25 * pow(Ng *(uniform_random_number_01 - cte), 2 );}; // uniform_random_number_01 is a random number between [0,1[ drawn from a uniform distribution
// };
//
//
// double test_gw10_distribution(int nb_tirages, double a) // test the distribution function defined by gw10_distribution by comparing the ratio of (number events in [a/2, a]) / (number events in [1/a, a/2[) with the theoretical value. there are nb_tirages events.
// {
//   unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
//   mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine
//   uniform_real_distribution<double> random01(0.0,1.0);
//   gw10_distribution distro(a);
//   double s =2.;
//   double gmin=0.;
//   double gmax = 0.;
//   double valth=( sqrt(a) - sqrt(a/s) ) / ( sqrt(a/s) - sqrt(1./a) ) ;
//
//   for(int i = 0; i < nb_tirages; i++)
//   {
//       if ( distro(random01(randomgen)) < a/s)
//           gmin += 1. ;
//       else
//           gmax += 1. ;
//   }
//   printf("\n test_gw10_distribution : gmax / gmin = %f . Différence ac valeur théorique (%f) : %f \n\n", gmax / gmin, valth, gmax / gmin - valth );
//   return gmax / gmin - valth;
// }
//
//
//
// void Save_chain(char * filename, double ** chain, int ndim, int nb_walkers, int chain_size, int chain_freq)
// {
//     int i , j;
//     FILE * myfile ;
//     myfile = fopen(filename, "w") ;
//     fprintf(myfile, "#  %i  %i  %i %i\n",ndim, nb_walkers, chain_size, chain_freq);
//     for(i=0 ; i < chain_size ; ++i) {
//         for (j = 0; j < ndim+1 ; ++j ) {
//              fprintf(myfile, "%.15e    ", chain[i][j] );
//          }
//          fprintf(myfile, "\n") ;
//      }
//      fclose(myfile);
// }
//
//
// void Load_chain(char * filename, double ** chain, int ndim, int nb_walkers, int & chain_size, int chain_freq)
// // Load a chain saved with "Save_chain". Allocate chain with chain_size + size of the chain in the file. Update chain_size.
// // Check that the characteristics of the loaded chain are the same as those of chain (walkers, frequency, dimension). Issue a Warning or an error if it is not the case.
// {
//     int i , j;
//     FILE * myfile ;
//     myfile = fopen(filename, "r") ;
//     int fndim=0;
//     int fchain_freq =0;
//     int fnb_walkers =0;
//     int fchain_size =0;
//     fscanf(myfile, "#  %i  %i  %i %i\n",&fndim, &fnb_walkers, &fchain_size, &fchain_freq);
//     if (ndim != fndim)
//     {
//         printf("\n Error ! Chain in file %s and desired chain have different number of dimensions : %i against %i \n\n", filename, fndim, ndim);
//         return;
//     }
//     chain = (double**) malloc((chain_size + fchain_size) * sizeof(double*) );
//     for (i =0  ; i < (chain_size + fchain_size) ; i++)
//     {
//         chain[i] = (double*) malloc((ndim+1)*sizeof(double));
//         for (j= 0 ; j < ndim + 1 ; j++) chain[i][j] = 0.;
//     }
//
//     for(i=0 ; i < fchain_size ; ++i) {
//         for (j = 0; j < fndim ; ++j ) {
//              fscanf(myfile, "%le    ", &(chain[i][j]) );
//          }
//          fscanf(myfile, "%le    \n", &(chain[i][fndim]) );
//      }
//      fclose(myfile);
//
//      if (nb_walkers != fnb_walkers)
//      {
//          printf("\n Warning ! nb_walkers = %i != fnb_walkers = %i \n\n", nb_walkers, fnb_walkers);
//      }
//      if (chain_freq != fchain_freq)
//      {
//          printf("\n Warning ! chain_freq = %i != fchain_freq = %i \n\n", chain_freq, fchain_freq);
//      }
//      chain_size = chain_size + fchain_size;
//
//      return;
// }
//
// int Init_from_prev_chain(char * prev_chain_file, double * walkers1d, const  int ndim, const int nb_walkers)
// {
//
//     int i , j;
//     FILE * myfile ;
//     myfile = fopen(prev_chain_file, "r") ;
//     int fndim=0;
//     int fchain_freq =0;
//     int fnb_walkers =0;
//     int fchain_size =0;
//     double ** chain;
//
//     fscanf(myfile, "#  %i  %i  %i %i\n",&fndim, &fnb_walkers, &fchain_size, &fchain_freq);
//     if (ndim != fndim)
//     {
//         printf("\n Error ! Chain in file %s and desired chain have different number of dimensions : %i against %i \n\n", prev_chain_file, fndim, ndim);
//         return 1 ;
//     }
//
//     if ( fchain_size < fnb_walkers )
//     {
//         printf("\n Error ! Chain in file %s is too small for initialization : length is %i against %i walkers needed \n\n", prev_chain_file, fchain_size, nb_walkers);
//         return 1 ;
//     }
//
//     chain = (double**) malloc((fchain_size) * sizeof(double*) );
//     for (i =0  ; i < (fchain_size) ; i++)
//     {
//         chain[i] = (double*) malloc((ndim+1)*sizeof(double));
//         for (j= 0 ; j < ndim + 1 ; j++) chain[i][j] = 0.;
//     }
//
//     for(i=0 ; i < fchain_size ; ++i) {
//         for (j = 0; j < fndim ; ++j ) {
//              fscanf(myfile, "%le    ", &(chain[i][j]) );
//          }
//          fscanf(myfile, "%le    \n", &(chain[i][fndim]) );
//      }
//      fclose(myfile);
//
//      if (nb_walkers != fnb_walkers)
//      {
//          printf("\n Warning ! nb_walkers = %i != fnb_walkers = %i \n\n", nb_walkers, fnb_walkers);
//      }
//
//
//      for (i = 0 ; i < nb_walkers ; i ++)
//      {
//          for (j = 0 ; j < ndim ; j++) walkers1d[i*ndim + j ] = chain[fchain_size - 1  - i][j] ;
//      }
//
//      for (i =0  ; i < (fchain_size) ; i++) free(chain[i]);
//      free(chain);
//
//      return 0 ;
// }


int main ( int argc, char *argv[] ) {

// MPI variable, by default for 1 proc. (no mpi)
    int code=0;
    int nbproc=1;
    int mpirank =0;
    int masterproc = 0; // process doing stuff alone
    int coeursparnoeud = 10000; // Used to do parallelism per node and not per thread.


// MPI initialization

#ifdef MPIMODE
    code = MPI_Init(&argc, &argv) ;
    code = MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    code = MPI_Comm_rank(MPI_COMM_WORLD , &mpirank);

    printf("mpirank %i \n\n", mpirank);
#endif

/* // This tests the distribution fiunction g
    double test = 0.;
    for (int i = 1; i < 10 ; i++) test += test_gw10_distribution( 1000000*i, 2.);
    printf("test : %f \n", test/10.);*/

// Input/output
    char  fileout[500] ;         // where the chain is saved
    char  parfile[500] ;        // passed to the function "fonction " object " (only use)
    char  datafile[500] ;       // passed to the function "fonction " object " (only use)

// MCMC variables
   // ********* to be set by user
   int ndim= 20;         // number of dimensions of the run, !! reset by Fittriple !!
   int nbw_per_proc = 1; // Number of walkers treated in parallel = Half of the number of walkers per processor
   int maxit = 0;           // max number of moves of the set of walkers
   // *************





   int halfnbw = 0; // Number of walkers divided by two
   int k,i,j = 0  ;      // index variable
   double accept = 0.; // acceptance probability
   double accepted = 0.; // number of accepted moves per core
   double total_accepted=0.; // number of accepted moves summed over all cores
   double acceptance_fraction=0. ; // acceptance fraction
   double acceptance_fraction_local=0. ; // acceptance fraction computed over the last moy_freq calls
   double parama = 2.; // parameter of the g probability function of GW10
  // double Ng = 2. * (sqrt(parama) - sqrt(1./parama)); // to transform uniform random gen to  g random generator
  // double cteg = - 2. / (sqrt(parama) * Ng ) ; // to transform uniform random gen to  g random generator
   gw10_distribution distribg(parama);
   int wj = 0 ;
   double z = 0;
   double lnchi2y = 0.;
   // Other intermediate variables
   double dlnchi2 = 0.;
   double mlnznm1 = 0.;



   double ** walkers;
   double * walkers1d;
   double * walkers1d_local; // for local proc before allgathering
   //double  walkery[ndim];
   double ** chaine;        // contain the chain in a matrix of chaine[maxit][ndim+1] where chaine[k][0:ndim] = walker[k%(2*halfnbw)] and chaine[k][ndim] = lnchi2[k%(2*halfnbw)]
   int chain_freq=1 ; // frequency with which the walkers states are added to the chain
   int chain_size=0; // = maxit * 2 * halfnbw /chaine_freq
   int quo = 0; // for use with remquo
   double rem = 0.; // for use with remquo
   int ensemble_chain_index = 0;

   // Recovery from previous chain
   char * prev_chain_file;
   prev_chain_file = NULL;

   // statistics
   //double variances[ndim];
   double ** moys;
   double moymoy = 0.;
   double varmoy = 0.;
   int moy_size =0;
   int moy_freq = 5000; // number of steps between two statistics (mean, autocorrelation time...) and backup
//     double moyacor[ndim];
//     double sigma[ndim];
//     double acortime[ndim];
    double * ensemble_chain;
    double * ensemble_chain2;
    int ensemble_chain_size = 0;



// Random generator init
    // obtain a seed from the system clock:
  unsigned seed1 = chrono::system_clock::now().time_since_epoch().count() + 1256 *mpirank; // seed different for every mpi process !
  mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine
  uniform_real_distribution<double> random01(0.0,1.0);
  normal_distribution<double> normaldist(0.,0.0001); // used for initialization


// Read arguments from command line and initialise
   if ( argc == 6 )
    {
        strcpy(parfile,argv[1]) ;
        strcpy(datafile, argv[2]);
        strcpy(fileout, argv[3]);
        sscanf(argv[4], "%i",&nbw_per_proc) ;
        sscanf(argv[5], "%i",&maxit) ;
        chain_freq = 1;
    }
    else if (argc == 7)
    {
        strcpy(parfile,argv[1]) ;
        strcpy(datafile, argv[2]);
        strcpy(fileout, argv[3]);
        sscanf(argv[4], "%i",&nbw_per_proc) ;
        sscanf(argv[5], "%i",&maxit) ;
        sscanf(argv[6], "%i",&chain_freq);
    }
    else if (argc == 8)
    {
        strcpy(parfile,argv[1]) ;
        strcpy(datafile, argv[2]);
        strcpy(fileout, argv[3]);
        sscanf(argv[4], "%i",&nbw_per_proc) ;
        sscanf(argv[5], "%i",&maxit) ;
        sscanf(argv[6], "%i",&chain_freq);
        sscanf(argv[7], "%i",&coeursparnoeud);
    }
    else if (argc == 9)
    {
        strcpy(parfile,argv[1]) ;
        strcpy(datafile, argv[2]);
        strcpy(fileout, argv[3]);
        sscanf(argv[4], "%i",&nbw_per_proc) ;
        sscanf(argv[5], "%i",&maxit) ;
        sscanf(argv[6], "%i",&chain_freq);
        sscanf(argv[7], "%i",&coeursparnoeud);
        sscanf(argv[8], "%i",&moy_freq);
    }
    else if (argc == 10)
    {
        strcpy(parfile,argv[1]) ;
        strcpy(datafile, argv[2]);
        strcpy(fileout, argv[3]);
        sscanf(argv[4], "%i",&nbw_per_proc) ;
        sscanf(argv[5], "%i",&maxit) ;
        sscanf(argv[6], "%i",&chain_freq);
        sscanf(argv[7], "%i",&coeursparnoeud);
        sscanf(argv[8], "%i",&moy_freq);
        prev_chain_file = argv[9];
    }
    else
    {
#ifdef MPIMODE
        code = MPI_Finalize();
#endif
        if (mpirank == masterproc) {
            printf("Not enough arguments passed ! Should be : parfile datafile fileout number_of_walkers_per_processor number_of_walker_set_moves\n");
            printf("Or : parfile datafile fileout number_of_walkers_per_processor number_of_walker_set_moves frequency_of_chain_extension\n");
            printf("Or : parfile datafile fileout number_of_walkers_per_processor number_of_walker_set_moves frequency_of_chain_extension number_of_proc_by_node\n");
            printf("Or : parfile datafile fileout number_of_walkers_per_processor number_of_walker_set_moves frequency_of_chain_extension number_of_proc_by_node stat_and_save_frequency\n");
        }
        coeursparnoeud = 10000;
        return 1;
    }


 // Last initializations/declarations depending on command line arguments
    halfnbw = nbw_per_proc * nbproc ;
    // defining what depends on halfnbw
    double lnchi2[2*halfnbw];
    double lnchi2_local[nbw_per_proc*2]; // For computation on one core before allgathering in lnchi2.
    /* lnchi2_local[0:nbw_per_proc] belong to [mpirank * nbw_per_proc : ( mpirank + 1) *nbw_per_proc] and
     * AND *
     * lnchi2_local[nbw_per_proc:2*nbw_per_proc] belong to [mpirank * nbw_per_proc + halfnbw : ( mpirank + 1) *nbw_per_proc + halfnbw]
     */
    uniform_int_distribution<int> randomwalker(0,halfnbw - 1);


   // **********************************************************************************************************************************
   // Declare an object called "fonction" or define a function somewhere that returns the log(chi2)
   // when called with fonction(double * parameters)

  //  ndim = 21; fctgaussienne fonction(ndim);
 fctFittriple_gaussprior fonction = fctFittriple_gaussprior();
 fonction.verbose=1;

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
   //************************************************************************************************************************************

 // Intitializing variables depending on ndim (in case fonction sets ndim)

   // walker varibles depending on ndim
   double  walkery[ndim];
   // statistics depending on ndim
   double variances[ndim];
   double moyacor[ndim];
   double moylocal[ndim];
   double sigma[ndim];
   double acortime[ndim];
   double moyacorhalf[ndim];
   double sigmahalf[ndim];
   double acortimehalf[ndim];

    for (i =0 ; i < ndim ; i++) // mere initialization
    {
        variances[i] = 0.;
        walkery[i] = 0.;
        moyacor[i] = 0.;
        moylocal[i] = 0.;
        sigma[i] = 0.;
        acortime[i] = 0.;
        moyacorhalf[i] = 0.;
        sigmahalf[i] = 0.;
        acortimehalf[i] = 0.;
    }

 // size of the the chain
    rem = fmod(static_cast<double>(maxit), static_cast<double>(chain_freq) );
    chain_size = int ( (static_cast<double>(maxit) - rem) / static_cast<double>(chain_freq) + 1. ) * 2 * halfnbw  ; //chain_size  * 2*halfnbw ; // saves iteration 0 hence the +1

//     if ( prev_chain_file != NULL )
//     {
//         printf("\n **** Initialization from previous chain in file %s \n", prev_chain_file);
//         if (Init_from_prev_chain(prev_chain_file, walkers1d, ndim, 2*halfnbw) == 1)
//         {
//             printf("Error : Initialization failed ! ");
//         }
//         printf("\n");
//     }

// Initializes statistical diagnostics
    if (mpirank == masterproc)
    {
        rem = fmod(static_cast<double>(maxit), static_cast<double>(moy_freq) );
        moy_size = int ( (static_cast<double>(maxit) - rem) / static_cast<double>(moy_freq) )  + 1;
        moys = (double**) malloc(moy_size * sizeof(double*));
        for (k = 0 ; k < moy_size ; k++)
        {
            moys[k] = (double *) malloc(ndim * sizeof(double) );
            for (i = 0 ; i < ndim ; i ++) moys[k][i] = 0.;
        }

        ensemble_chain_size = chain_size / (2 * halfnbw);
        ensemble_chain = (double*) malloc(ensemble_chain_size  * sizeof(double) );
        ensemble_chain2 = (double*) malloc(ensemble_chain_size  * sizeof(double) );
    }

  // Walkers memory allocations
   walkers = (double **) malloc(2*halfnbw * sizeof(double*)) ;
   walkers1d = (double*) malloc(2*halfnbw * ndim * sizeof(double)) ;
   walkers1d_local = (double*) malloc(nbw_per_proc * ndim * sizeof(double) );
   for ( k = 0 ; k < 2*halfnbw ; k++) walkers[k] = &walkers1d[k*ndim];

   // initialization of the walker in a gaussian ball or from previous chain
   if (mpirank == masterproc)
   {
    if ( prev_chain_file != NULL )
        {
            printf("\n **** Initialization from previous chain in file %s \n", prev_chain_file);
            if (Init_from_prev_chain(prev_chain_file, walkers1d, ndim, 2*halfnbw) == 1)
            {
                printf("Error : Initialization failed ! ");
                return 1;
            }
            printf("\n");
        }
        else
        {
            for ( i = 0 ; i < ndim*2*halfnbw ; i++) walkers1d[i] = normaldist(randomgen);
        }
//         // Test // To be removed !
//    for ( k = 0 ; k < 2*halfnbw ; k++) walkers[k][5] = 23.37 + normaldist(randomgen); // force 30km/s for rv
    };
#ifdef MPIMODE
   code = MPI_Bcast(walkers1d, ndim*2*halfnbw, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
#endif



   // initialization of chaine
   if (mpirank==masterproc )
   {
       chaine = (double **) malloc(chain_size *sizeof(double*)) ;
       for (k = 0 ; k < chain_size ; k++) chaine[k] = (double*) malloc((ndim + 1) * sizeof(double));

   }

// end of Memory alocations


   if (mpirank == masterproc) {
   printf("\n*****************************************************************************************************\n");
   printf(" Function parfile= %s and datafile= %s \n", parfile, datafile);
   fonction.Print_MCMC();
  //printf(" After initialization (and before MCMC) the chi2 is : %.5f\n", fonction.Compute_lnposterior(0) );
  printf(" Number of dimensions : %i \n", ndim );
   printf(" There are %i walkers spread on %i processors. \n", 2*halfnbw, nbproc);
   printf(" The whole set of walkers will undergo %i moves.\n", maxit);
   printf(" The state will of the set of walkers will be saved every %i moves, for a total of %i walkers saved.\n", chain_freq, chain_size);
   printf(" The length of the ensemble chain is thus %i \n", ensemble_chain_size);
   printf(" The resulting chain + lnchi2 will be saved in %s every %i iterations. \n", fileout, moy_freq);
   printf("*******************************************************************************************************\n\n");



    printf("\n Initializing the walkers with chi2. \n\n");
   }
    for ( k = mpirank * nbw_per_proc ; k < ( mpirank + 1) *nbw_per_proc ; k++)
    {
        lnchi2_local[k- mpirank * nbw_per_proc] = fonction(walkers[k]);
        lnchi2_local[(k- mpirank * nbw_per_proc) + nbw_per_proc] = fonction(walkers[k + halfnbw]);
        //lnchi2[k + halfnbw] = fonction(walkers[k + halfnbw]);
    }
#ifdef MPIMODE
    code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, lnchi2, nbw_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);
    code = MPI_Allgather(lnchi2_local + nbw_per_proc , nbw_per_proc, MPI_DOUBLE, lnchi2 + halfnbw, nbw_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);
/*
    for ( k = mpirank * nbw_per_proc ; k < ( mpirank + 1) *nbw_per_proc ; k++)
    {
        lnchi2_local[k- mpirank * nbw_per_proc] = fonction(walkers[k + halfnbw]);
        //lnchi2[k + halfnbw] = fonction(walkers[k + halfnbw]);
    }
    code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, lnchi2 + halfnbw, nbw_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);*/
    //code = MPI_Allgather(&lnchi2[mpirank * nbw_per_proc + halfnbw] , nbw_per_proc, MPI_DOUBLE, lnchi2 + halfnbw, nbw_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);
#else
    memcpy(lnchi2, lnchi2_local, 2* sizeof(double)* nbw_per_proc);
#endif
  if (mpirank==masterproc ) // Print initial state of the chain
   {
      //Print_table(walkers1d, ndim*2*halfnbw);

       printf("\n*** Initial chain state : \n");
       for (k = 0 ; k < 2* halfnbw; k++)
        {
            printf("walker%i | ", k);
            for (i = 0; i < ndim ; i++)  printf("%.4e ", walkers[k][i]);
            printf(" | %.5e \n", lnchi2[k]);
        }
       printf("\n");
   }

   fonction.verbose=0; // minimal printing

    if (mpirank == masterproc) printf("\n Starting the actual MCMC... \n\n");


 // **********   Main loop
    for (int it = 0 ; it < maxit ; it++)
    {
#ifdef MPIMODE
        code = MPI_Barrier ( MPI_COMM_WORLD);
#endif
        for ( k = mpirank * nbw_per_proc ; k < ( mpirank + 1) *nbw_per_proc ; k++)
        {
            // Take care of walker "k"
            // Draw "z" for stretch move and newwalker number wj
            wj = randomwalker(randomgen) + halfnbw;
            // G = 2 /sqrt (z)
            z = distribg(random01(randomgen)) ;
            // Determine walkery the possible new position of walker "k"
            for (i=0; i < ndim ; i++) walkery[i] = walkers[wj][i] + z * ( walkers[k][i] - walkers[wj][i] ) ;
            // Compute p(Y) and compute acceptance accept
            lnchi2y = fonction(walkery);

            // Determine if the step is accepted by computing accept = z^(ndim-1)*p(Y) / p(X_k) = pow(z, ndim - 1)*exp(lnchi2y - lnchi2[k])
            dlnchi2 = lnchi2y - lnchi2[k];
            mlnznm1 = -(ndim-1)*log(z);
            if (dlnchi2 > mlnznm1) // necessarily accepted since i this case accept > 1
                accept = 2.;
            else
                accept = exp(-mlnznm1 + dlnchi2);
            // draw r in [0,1] to accept or reject the move
            if (random01(randomgen) <= accept)
            {
                for (i=0; i < ndim ; i++) walkers[k][i] = walkery[i];
                lnchi2[k] = lnchi2y;
                accepted += 1.;
            }
        }
#ifdef MPIMODE
        memcpy(lnchi2_local, &lnchi2[mpirank * nbw_per_proc] , nbw_per_proc*sizeof(double) );
        memcpy(walkers1d_local, &walkers1d[(mpirank * nbw_per_proc )*ndim], nbw_per_proc*ndim * sizeof(double) );
        for(i = 0; i < nbw_per_proc ; i++)
        {
            if (lnchi2_local[i] != lnchi2[mpirank * nbw_per_proc + i]) printf("\n\nCATASTROPHE !! \n\n");
        }
        code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, lnchi2, nbw_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);
        code = MPI_Allgather(walkers1d_local, nbw_per_proc*ndim, MPI_DOUBLE, walkers1d, nbw_per_proc*ndim, MPI_DOUBLE, MPI_COMM_WORLD);
#endif

        for ( k = mpirank * nbw_per_proc + halfnbw ; k < ( mpirank + 1) *nbw_per_proc + halfnbw ; k++)
        {

            // Take care of walker "k"
            // Draw "z" for stretch move and newwalker number wj
            wj = randomwalker(randomgen) ;
            // G = 2 /sqrt (z)
            z = distribg(random01(randomgen)) ;//0.25 * pow(Ng *(random01(randomgen) - cteg), 2 );
            // Determine walkery the possible new position of walker "k"
            for (i=0; i < ndim ; i++) walkery[i] = walkers[wj][i] + z * ( walkers[k][i] - walkers[wj][i] ) ;
            // Compute ln( p(Y) ) and compute acceptance accept
            lnchi2y = fonction(walkery);

            // Determine if the step is accepted by computing accept = z^(ndim-1)*p(Y) / p(X_k) = z^(ndim-1)*exp(lnchi2y - lnchi2[k])
            dlnchi2 = lnchi2y - lnchi2[k];
            mlnznm1 = -(ndim-1)*log(z);
            if (dlnchi2 > mlnznm1) // necessarily accepted since i this case accept > 1
                accept = 2.;
            else
                accept = exp(-mlnznm1 + dlnchi2); //  pow(z, ndim - 1) * exp(lnchi2y - lnchi2[k]) ;
            // draw r in [0,1] to accept or reject the move
            if (random01(randomgen) <= accept)
            {
                for (i=0; i < ndim ; i++) walkers[k][i] = walkery[i];
                lnchi2[k] = lnchi2y;
                accepted += 1.;
            }
        }

#ifdef MPIMODE
        memcpy(lnchi2_local, &lnchi2[mpirank * nbw_per_proc + halfnbw] , nbw_per_proc*sizeof(double) );
        memcpy(walkers1d_local, &walkers1d[(mpirank * nbw_per_proc + halfnbw)*ndim], nbw_per_proc*ndim * sizeof(double) );

        code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, &lnchi2[halfnbw], nbw_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);
        code = MPI_Allgather(walkers1d_local, nbw_per_proc*ndim, MPI_DOUBLE, &walkers1d[halfnbw*ndim], nbw_per_proc*ndim, MPI_DOUBLE, MPI_COMM_WORLD);

#endif



       if (mpirank == masterproc) // Record the state of the chain, do diagnostics..
       {
           rem = fmod(static_cast<double>(it + 1) , static_cast<double>(chain_freq));
           ensemble_chain_index = int( (static_cast<double>(it + 1) - rem ) / static_cast<double>(chain_freq) );

           if (rem == 0. or it == 0)
           {
               if (chain_freq == 1)  ensemble_chain_index -= 1; // ensemble_chain_index = nb of chain_freq done, but index of array starts at 0;
               printf("Saving chain after %i iterations. %i walkers saved.  \n", it, ensemble_chain_index * 2 * halfnbw);
                for(k=0; k < 2*halfnbw ; k++)
                {
                    for (i = 0; i < ndim ; i++) chaine[ensemble_chain_index*2*halfnbw + k][i] = walkers[k][i];
                    chaine[ensemble_chain_index*2*halfnbw + k][ndim] = lnchi2[k];
                }
           }


           rem = fmod(static_cast<double>(it ) , static_cast<double>(moy_freq));
           quo = int( (static_cast<double>(it ) - rem ) / static_cast<double>(moy_freq) );

           for (i = 0 ; i < ndim ; i++)
           {
               for (k = 0 ; k < 2*halfnbw ; k++)
               {
                   moys[quo][i] += walkers[k][i];
                   variances[i] += pow(walkers[k][i],2);
               }
           }
       }


           // Make periodic backup and statistics
           rem = fmod(static_cast<double>(it+1 ) , static_cast<double>(moy_freq));
           //quo = int( (static_cast<double>(it ) - rem ) / static_cast<double>(moy_freq) );
           if (rem ==0.)
           {
#ifdef MPIMODE
               code = MPI_Reduce(&accepted , &total_accepted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // Compute total number of accepted moves
#else
               total_accepted = accepted;
#endif
                if (mpirank == masterproc)
                {
                printf("\n");
                                    printf(" %i over %i and %i \n\n", ensemble_chain_index, ensemble_chain_size, it  );


                printf("* Ensemble chain step : %i \n", ensemble_chain_index);
                for (i = 0 ; i < ndim ; i++)
                {
                    for (k = 0 ; k < (ensemble_chain_index+1) ; k++)
                    {   // Create an ensemble chain with each parameter averaged over the set of walkers
                        ensemble_chain[k] = 0.;
                        for ( j = 0 ; j < 2*halfnbw ; j ++) ensemble_chain[k] += chaine[k* (2*halfnbw) + j][i];
                        ensemble_chain[k] /= 2*halfnbw ;
                    }

                    // Compute the mean on the last moy_freq iterations
                    moylocal[i] = 0.;
//                     for (j = it - moy_freq ; j < it ; j ++)
//                     {
//                         moylocal[i] += chaine[j][i];
//                     }
//                     moylocal[i] /= moy_freq;
                    memcpy(ensemble_chain2, &ensemble_chain[ensemble_chain_index - ensemble_chain_index/2 ],  (ensemble_chain_index/2 +1)*sizeof(double) ); // Necessary to make a copy because acor changes it
                    acor( &moyacorhalf[i], &sigmahalf[i], &acortimehalf[i], ensemble_chain2, (ensemble_chain_index/2 +1));
                    acor( &moyacor[i], &sigma[i], &acortime[i], ensemble_chain, (ensemble_chain_index+1)  ); // Run the GW10 autocorrelation algorithm

            // Run the GW10 autocorrelation algorithm

                    printf("Iteration %i ,  parameter %i : Mean %.5e (last %i iterations %.5e) | Std(Mean) %.5e | AcorTime %f  \n", it, i, moyacor[i], moy_freq, moylocal[i], sigma[i], acortime[i] );
                    printf("                              Second Half of the chain : Mean %.5e | Std(Mean) %.5e | AcorTime %f  \n\n",  moyacorhalf[i], sigmahalf[i], acortimehalf[i] );
                }
                printf("\n");
                printf("Chi2 value at mean of the seconf half of the chain : %.10e\n", fonction(moyacorhalf) );
          //      printf("Chi2 value at mean evaluated on the last %i iterations : %.10e\n", moy_freq, fonction(moyacor) );

                acceptance_fraction_local = ( total_accepted - acceptance_fraction * ((it-moy_freq) * 2 * halfnbw) ) / (moy_freq*2 * halfnbw) ; // acceptance fraction on the last moy_freq iterations
                acceptance_fraction = total_accepted/ (it * 2 * halfnbw);
                printf("\n");
                printf("Acceptance fraction on the last %i iterations is : %.5e\n", moy_freq, acceptance_fraction_local);
                printf("Total acceptance fraction is : %.5e\n", acceptance_fraction);

                printf("\n Saving data to %s \n", fileout);
               // Savetxt(fileout, chaine, (quo+1)*2*halfnbw, ndim+1 );
                Save_chain(fileout, chaine, ndim, 2*halfnbw, ensemble_chain_index * 2 * halfnbw, chain_freq);

                printf("\n\n");


                }
            }



    } // ******** end of main loop

#ifdef MPIMODE
    code = MPI_Reduce(&accepted , &acceptance_fraction, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // Compute total number of accepted moves
#else
    acceptance_fraction = accepted;
#endif

    acceptance_fraction /= (maxit * 2 * halfnbw);

    int acorfail = 0;
    int maxacortimei = 0;
    int burnin = 0;
    double maxacortime = 0.;

    // Save the result "chaine" in "fileout"
    if (mpirank == masterproc)
    {

        printf("\n Acor analysis for final chain : \n");
        acorfail = 0;
        maxacortimei = -1;

               for (i = 0 ; i < ndim ; i++)
               {

                    for (k = 0 ; k < ensemble_chain_size ; k++)
                    {   // Create a ensemble chain with each parameter averaged over the set of walkers
                        ensemble_chain[k] = 0.;
                        for ( j = 0 ; j < 2*halfnbw ; j ++) ensemble_chain[k] += chaine[k * (2*halfnbw) + j][i];
                        ensemble_chain[k] /= 2*halfnbw ;
                    }


                    if (0 == acor( &moyacor[i], &sigma[i], &acortime[i], ensemble_chain, ensemble_chain_size )) // Run the GW10 autocorrelation algorithm
                    {
                        if (acortime[i] > maxacortime )
                        {
                            maxacortimei = i;
                            maxacortime = acortime[i];
                        }
                    }
                    else
                        acorfail = 1;
                    printf("Parameter %i : Mean %f | Std(Mean) %f | AcorTime %f  | r %f | burnin %i \n", i, moyacor[i], sigma[i], acortime[i], acortime[i]/ (ensemble_chain_size - burnin), burnin );
                }
       printf("\n");
       printf("Chi2 value at mean : %.10e\n", fonction(moyacor) );

               if (maxacortimei >= 0)
               {
                   if (20 * acortime[maxacortimei] <  ensemble_chain_size)
                   {
                       if (acorfail == 1)
                           printf("\n acor failed sometimes, retrying after discarding 20 times the longest acor time which is parameter %i, %f... \n\n", maxacortimei, maxacortime);
                       else
                           printf("\n Attempt to estimate burnin time : retrying after discarding 20 times the longest acor time which is parameter %i, %f... \n\n", maxacortimei, maxacortime);
                       burnin = int( maxacortime )*20;
                       for (i = 0 ; i < ndim ; i++)
                       {
                           for (k = 0 ; k < ensemble_chain_size ; k++)
                            {   // Create a ensemble chain with each parameter averaged over the set of walkers
                                ensemble_chain[k] = 0.;
                                for ( j = 0 ; j < 2*halfnbw ; j ++) ensemble_chain[k] += chaine[k * (2*halfnbw) + j][i];
                                ensemble_chain[k] /= 2*halfnbw ;
                            }
                           acor( &moyacor[i], &sigma[i], &acortime[i], &ensemble_chain[burnin], ensemble_chain_size - burnin );
                           printf("Parameter %i : Mean %f | Std(Mean) %f | AcorTime %f  | r %f | burnin %i \n", i, moyacor[i], sigma[i], acortime[i], acortime[i]/ (ensemble_chain_size - burnin), burnin );
                       }
                   } else
                       printf("\nacor failed sometimes but there is not enough data even to remove the longest burning time from valid acor times\n");
               }
               //printf("acor failed sometimes. maxacortime %f %i aacorfail %i \n\n", maxacortime, ensemble_chain_size, acorfail);


        printf("\n");
        printf("Acceptance fraction is : %.5e\n", acceptance_fraction);

        printf("\n\n");
        printf("Saving data to %s \n", fileout);
       // Savetxt(fileout, chaine, chain_size, ndim+1 );
        Save_chain(fileout, chaine, ndim, 2*halfnbw, chain_size, chain_freq);

    }
// // //
    free(walkers);
    free(walkers1d);
    free(walkers1d_local);
    if (mpirank == masterproc)
    {
        for (k = 0 ; k < chain_size ; k++) free(chaine[k]);
        free(chaine);

        for (k = 0 ; k < moy_size ; k++)
        {
            free(moys[k]);
        }
        free(moys);
        free( ensemble_chain);
        free( ensemble_chain2);
    }



#ifdef MPIMODE
    code = MPI_Finalize();
#endif

  return 0;
};
