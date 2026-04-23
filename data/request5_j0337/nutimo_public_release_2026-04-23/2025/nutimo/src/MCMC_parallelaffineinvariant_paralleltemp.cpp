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
#include <stdio.h>
#include <string.h>
#include <time.h>

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

class decovardistance
{
    int ndim;
    double ** canonical2diagonal;
    double ** loadingarray;
    double * pref;
    double * sigs;
public:
    decovardistance(char * filename, int nb_of_dim)
    {
        int i ;

        ndim = nb_of_dim;

        pref = (double*)malloc(sizeof(double)*ndim);
        sigs = (double*)malloc(sizeof(double)*ndim);
        canonical2diagonal = (double**) malloc(sizeof(double*) * ndim);
        loadingarray = (double**) malloc(sizeof(double*) * (ndim+2));
        for (int i = 0 ; i < ndim; i++)
        {
            canonical2diagonal[i] = (double*) malloc(ndim*sizeof(double));
            loadingarray[i] = (double*) malloc(ndim*sizeof(double));
        }
        Loadtxt(filename, loadingarray, ndim+2, ndim ) ;
        memcpy(pref, loadingarray[0], ndim*sizeof(double));
        memcpy(sigs, loadingarray[1], ndim*sizeof(double));
        for (int i = 0 ; i < ndim; i++) memcpy(canonical2diagonal[i], loadingarray[i+2], ndim*sizeof(double));
        for (int i = 0 ; i < ndim+2; i++) free(loadingarray[i]);
        free(loadingarray);
    };

    ~decovardistance()
    {
        for (int i = 0 ; i < ndim; i++)
        {
            free(canonical2diagonal[i]);
        }
        free(canonical2diagonal);
        free(sigs);
        free(pref);
    };

    double operator()(double * pars)
    {
        int i,j;
        double decopars[ndim];

        for ( i = 0 ; i < ndim ; i++)
        {
            decopars[i] = 0.;
            for ( j = 0 ; j < ndim ; j++) decopars[i] += canonical2diagonal[i][j] * (pars[j] - pref[j]);
        }

        double dst = 0.;
        for ( j = 0 ; j < ndim ; j++) dst += decopars[j]* decopars[j] / ( sigs[j] * sigs[j]);

        return sqrt(dst);
    };

};







int main ( int argc, char *argv[] ) {
/* **** Parallelization scheme ****
 *
 *         Temp 0,1..ntemp_per_cluster -1           ntemp_per_cluster - 2*ntemp_per_cluster -1   ...
 *|----------------------------------------------|----------------------------------------------|...
 *        Cluster 0                                                 Cluster 1                    ...Cluster ncluster-1
 *
 * Each cluster has nproc_per_temp processors, and computes ntemp_per_cluster temperatures.
 * Clusters have their own communicator MPI_COMM_tempcluster
 * Each cluster has a master processor characterized by mpirankincluster =0
 *
 * Each temperature has 2 * nbw_per_proc * nproc_per_temp walkers
 *
 * The only time all walkers of all temperatures of all clusters are gathered, it is on the masterproc  for temperature swaps.
 *
 *
 * TODO : Memory optimization: chaine, chaine_refused, Temp_walkers1d, Temp_lnchi2 and others do not need to be fully allocated
 * in every cluster, but only for the number of temperatures dealt with by the cluster.
 *******************************************/


// MPI variable, by default for 1 proc. (no mpi)
    int code=0;
    int nbproc=1;
    int mpirank =0;
    const int masterproc = 0; // process doing stuff alone. NEED TO BE = 0
    int coeursparnoeud = 10000; // Used to do parallelism per node and not per thread.

    // For temperature parallelization
    int nproc_per_temp =-1; // number of processor for each temperature
    int ntemp_per_cluster = -1; // number of temperatures computed by one cluster
    int ncluster = 1;       // number of cluster of processors working on different sets of temperature
    int icluster = 0;       // index of the cluster of the current processor
    int mpirankincluster = 0;   // rank of a processor within the cluster it belongs to
    MPI_Comm MPI_COMM_tempcluster; // Communicator within one cluster


// MPI initialization

#ifdef MPIMODE

    code = MPI_Init(&argc, &argv) ;
    code = MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    code = MPI_Comm_rank(MPI_COMM_WORLD , &mpirank);
#endif

/* // This tests the distribution fiunction g
    double test = 0.;
    for (int i = 1; i < 10 ; i++) test += test_gw10_distribution( 1000000*i, 2.);
    printf("test : %f \n", test/10.);*/

// Input/output
    char  fileout[500] ;        // where the chain is saved
    char  parfile[500] ;        // passed to the function "fonction " object " (only use)
    char  datafile[500] ;       // passed to the function "fonction " object " (only use)
    char * turnfile;         // passed to the function "fonction " object " (only use)
    turnfile = NULL;

// Printing
    int size_strtoprint=2000;
    char * strtoprint;
    char strinter[2000];
    int * printers; // list of process rank that print
    char paramname[100];
    clock_t tclock_start=clock();


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
   int wj = 0 ;
   double z = 0;
   double lnchi2y = 0.;
   // Other intermediate variables
   double dlnchi2 = 0.;
   double mlnznm1 = 0.;

   // Experimental
   // int total_countersmoothing=0;
   // int countersmoothing=0;
   // double smoothinglevel = 0.5;



   double ** walkers;
   double * walkers1d;
   double * walkers1d_local; // for local proc before allgathering
   //double  walkery[ndim];
   double *** chaine;        // contain the chain in a matrix of chaine[ntemp][maxit][ndim+1] where chaine[itemp][k][0:ndim] = walker[k%(2*halfnbw)] and chaine[itemp][k][ndim] = lnchi2[k%(2*halfnbw)] for temperature itemp
   double *** chaine_refused; // optionally contains the refused jumps
   bool record_refused=false; // If true, record refused jumps in chaine_refused
   int irefused=0 ;
   int * irefuseds;
   int chain_freq=1 ; // frequency with which the walkers states are added to the chain
   int chain_size=0; // = maxit * 2 * halfnbw /chaine_freq
   int quo = 0; // for use with remquo
   double rem = 0.; // for use with remquo
   double printstat =0.;
   int ensemble_chain_index = 0;
   int ntry = 0 ; // Trial counter

   // Recovery from previous chain
   char * prev_chain_file;
   prev_chain_file = NULL;

   // statistics
   //double variances[ndim];
   bool calc_mean_chi2 = true;
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
    double * autocorrelation;


// Temperature variables
    int ntemp = 1;  // Number of temperatures
    double tempfactor[100] ; // ratios between two consecutive temperatures
    double * betas; // array of inverses of temperatures
    double targettemp = 1.; // target temperature
    int tempswapfreq =1; // frequency of temperature swap
    double * Temp_accepted; // number of accepted moves per core for each temperature (see ''accepted'')
    double * Temp_acceptance_fraction; // ''acceptance_fraction'' for each temperature
    double ** Temp_lnchi2; // ''lnchi2'' for each temperature
    double ** Temp_walkers1d; // ''walkers1d'' for each temperature
    int itemp=0; // simple indexing variable
    double * walkerbuffer; // Buffer that contains 1 walker at a time;
    char  fileouttemp[500] ;    // buffer for different temperatures (1 file / temperature)
    double * accepted_tempswaps; // average probability of accepting temperature swaps for each temperature
    double * last_accepted_tempswaps; // Same but reset every moy_freq iterations
    int kex=0;
    int kshift=0; // Permutation of swapping walkers
    int * Temp_irefused ;

// Random generator init
    // obtain a seed from the system clock:
  unsigned seed1 = chrono::system_clock::now().time_since_epoch().count() + 1256 *mpirank; // seed different for every mpi process !
  mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine
  uniform_real_distribution<double> random01(0.0,1.0);


// Diagnostics of the MCMC algorithm
  bool dodiagnostics = false;
  double ** diagnosticslog; // see header below for the columns
  char diagnosticsheader[200];
// !!!! THE FOLLOWING LINE CAUSED random01 TO RETURN ABSURD RESULTS ON HYDRUS !!!
//  sprintf(diagnosticsheader, "Iteration   |   Temperature   |   Half-iteration   |   Current Walker  |  Candidate Walker   |   Step factor z   |   Current lnposterior   |   Candidate lnposterior   |   Acceptance probability   |   Acceptance outcome ");
  int ndiagnostics=0;   // number of diagnostics lines (ie 2 / temperature / iteration)
  int idiagnostics = 0;
  const int diagnostics_columns = 10; // nb of diagnostics per line
  char diagnosticsfile[100]; // file name to save diagnostics
  char walkersdiagnosticsfilename[100]; // save walkers to files from each core to check MPi is ok
  FILE * walkersdiagnosticsfile ;

// Sieving of excessively low lnposteriors :
bool dosieving = false;
int maxlnchi2index =0;
double maxlnchi2 = 0.;
int isieved = 0;
int maxsieved = -1;
double maxtension = 0.;
double ** sievinglog;
int ncol_sievinglog=0;
int wk = 0;
int max_sieving_attempts = 15; // TODO : should be user defined.

// Diverse
  double inter = 0.;

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
//     else if (argc == 9)
//     {
//         strcpy(parfile,argv[1]) ;
//         strcpy(datafile, argv[2]);
//         strcpy(fileout, argv[3]);
//         sscanf(argv[4], "%i",&nbw_per_proc) ;
//         sscanf(argv[5], "%i",&maxit) ;
//         sscanf(argv[6], "%i",&chain_freq);
//         sscanf(argv[7], "%i",&coeursparnoeud);
//         sscanf(argv[8], "%i",&moy_freq);
//     }
    else if (argc >= 9)
    {
        strcpy(parfile,argv[1]) ;
        strcpy(datafile, argv[2]);
        strcpy(fileout, argv[3]);
        sscanf(argv[4], "%i",&nbw_per_proc) ;
        sscanf(argv[5], "%i",&maxit) ;
        sscanf(argv[6], "%i",&chain_freq);
        sscanf(argv[7], "%i",&coeursparnoeud);
        sscanf(argv[8], "%i",&moy_freq);
        if (mpirank == masterproc)
        {
            printf("\n *** Reading ordered parameters :\n");
            printf("parfile : %s\n", parfile);
            printf("datafile : %s\n", datafile);
            printf("fileout : %s\n", fileout);
            printf("nbw_per_proc : %i\n", nbw_per_proc);
            printf("maxit : %i\n", maxit);
            printf("chain_freq : %i\n", chain_freq);
            printf("coeursparnoeud : %i\n", coeursparnoeud);
            printf("moy_freq : %i\n", moy_freq);
        }
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
            printf("Optional keyword options :\n");
            printf("  turnfile <filename>\n");
            printf("  previouschain <filename_prefix>\n");
            printf("  turnfile <filename>\n");
            printf("  recordrefused\n");
            printf("  aparam <value>\n");
            printf("  targettemp <float value>\n");
            printf("  tempswapfreq <int value>\n");
            printf("  ntemp <int value>\n");
            printf("  tempfactor <float value> <float value> ...\n");
            printf("  nproc_per_temperature <int value>\n");
            printf("  --sieving <float value> \n");
            printf("  --diagnostics\n");
            printf("  --no_mean_calc\n    : does not calculate chi2 at mean walkers during diagnostics.");
        }
        coeursparnoeud = 10000;
        return 1;
    }


//**  Reading Keyword arguments
        if (mpirank == masterproc) printf("\n *** Reading keyword parameters :\n");

            for (i = 1 ; i < argc ; i++)
            {
//                if (mpirank == masterproc) printf(" argc i argv[i] %i, %i, %s \n" , argc, i, argv[i]); // test
                if (strcmp(argv[i], "turnfile") == 0)
                {
                    if (argc > i+1)
                    {
                        turnfile =  argv[i+1] ;
                        if (mpirank == masterproc) printf("Reading turnfile : %s \n", turnfile);
                    }
                    else
                    {
                        printf("Failed to read turnfile !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "previouschain") == 0)
                {
                    if (argc > i+1)
                    {
                        prev_chain_file =  argv[i+1] ;
                        if (mpirank == masterproc) printf("Reading previouschain : %s \n", prev_chain_file);
                    }
                    else
                    {
                        printf("Failed to read previouschain !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "recordrefused") == 0)
                {
                    record_refused=true;
                    if (mpirank == masterproc) printf("Reading record_refused : %i \n", int(record_refused));
                }
                if (strcmp(argv[i], "aparam") == 0)
                {
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%lf",&parama) ;
                        if (mpirank == masterproc) printf("Reading parama : %.2f \n", parama);
                    }
                    else
                    {
                        printf("Failed to read aparam !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "targettemp") == 0)
                {
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%lf",&targettemp) ;
                        if (mpirank == masterproc) printf("Reading targettemp : %f \n", targettemp);
                    }
                    else
                    {
                        printf("Failed to read targettemp !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "tempswapfreq") == 0)
                {
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%i",&tempswapfreq) ;
                        if (mpirank == masterproc) printf("Reading tempswapfreq : %i \n", tempswapfreq);
                    }
                    else
                    {
                        printf("Failed to read tempswapfreq !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "ntemp") == 0)
                {
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%i",&ntemp) ;
                        if (mpirank == masterproc) printf("Reading ntemp : %i \n", ntemp);
                    }
                    else
                    {
                        printf("Failed to read ntemp !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                    for (j = 0; j < ntemp - 1 ; j++) tempfactor[j] = sqrt(2.);
                    if (ntemp>1)
                    {
                        for (k = 0 ; k < argc ; k++)
                        {
                            if (strcmp(argv[k], "tempfactor") == 0)
                            {
                                j = 0;
                                while ( j < min(ntemp-1, argc - k - 1))
                                {
                                    if (sscanf(argv[k+1+j], "%lf",&tempfactor[j]) == 1)
                                    {
                                        j++ ;
                                        if (mpirank == masterproc) printf("Reading tempfactor %i : %.2f \n", j, tempfactor[j-1]);
                                    }
                                }
                                if (j > 0)
                                {
                                    while(j+1 < ntemp-1)
                                    {
                                        tempfactor[j+1] = tempfactor[j];
                                        j++;
                                    }
                                }
                                else
                                {
                                    printf("Failed to read tempfactor !\n");
    #ifdef MPIMODE
                                    code = MPI_Finalize();
    #endif
                                    return 1;
                                }
                            }
                        }

                    }
                }
                if (strcmp(argv[i], "nproc_per_temperature") == 0)
                {
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%i",&nproc_per_temp) ;
                        if (mpirank == masterproc) printf("Reading nproc_per_temp : %i \n", nproc_per_temp);
                    }
                    else
                    {
                        printf("Failed to read nproc_per_temp !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "--sieving") == 0)
                {
                    dosieving = true;
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%lf",&maxtension) ;
                        // maxtension = 0.5*maxtension; // chi2 -> posterior
                        if (mpirank == masterproc) printf("Reading sieving maxtension : %.5e \n", maxtension); //(DIVIDED/2 cause chi2 -> posterior)
                    }
                    else
                    {
                        printf("Failed to read sieving maxtension !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "--maxsieved") == 0)
                {
                    if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%i",&maxsieved) ;
                        if (mpirank == masterproc) printf("Reading sieving maxsieved : %i \n", maxsieved);
                    }
                    else
                    {
                        printf("Failed to read sieving maxsieved !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
                }
                if (strcmp(argv[i], "--no_mean_calc") == 0)
                {
                    calc_mean_chi2 = false;
                }
                if (strcmp(argv[i], "--diagnostics") == 0)
                {
                    dodiagnostics = true;
                }
            }
        printf("\n");
 // End of reading keywords arguments

// Temperature parallelization variables
    if (nproc_per_temp == -1 ) nproc_per_temp = nbproc; // if nproc_per_temp was not set by user

    icluster = mpirank / nproc_per_temp;
    ncluster = nbproc / nproc_per_temp;
    ntemp_per_cluster = ntemp / ncluster;

    if (ntemp_per_cluster * ncluster < ntemp )
    {
        while ( ntemp%ncluster != 0.) ncluster--; //  find the number of cluster actually used
        ntemp_per_cluster = ntemp/ncluster;
        printf("\n Warning: the number of temperatures is not a multiple of the number of clusters. Only %i/%i processors used ! \n\n", ncluster*nproc_per_temp, nbproc);
    }
    mpirankincluster = mpirank%nproc_per_temp;
    printf("mpirank %i mpirankincluster %i icluster %i ncluster %i ntemp_per_cluster %i \n", mpirank, mpirankincluster, icluster, ncluster, ntemp_per_cluster);
    // Create the communicator for each cluster and define the corresponding printers
#ifdef MPIMODE
    code = MPI_Comm_split(MPI_COMM_WORLD, icluster, mpirank, &MPI_COMM_tempcluster);
    printers = (int*)malloc(sizeof(int) * ncluster);
    for (k = 0 ; k < ncluster; k++)
    {
        printers[k] = k * nproc_per_temp;
    };
#endif

// Last initializations/declarations depending on command line arguments
    halfnbw = nbw_per_proc * nproc_per_temp ;
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

//   ndim = 2; fctgaussienne fonction(ndim);

//      fctsimudata fonction(10000,1., 0.1);
//      ndim = 2; fonction.errorscale =1;
//    ndim=2; fctmultigaussienne fonction(ndim, 3);
//    fonction.errorscale =1;
//    ndim = 3;
//
// double ** xy;
// double xx[2];
// double ys[100];
// xy = (double**)malloc(sizeof(double*) * 100);
// for (i = 0;  i < 100 ; i++) xy[i] = (double*) malloc(sizeof(double) * 2);
// for (i = 0;  i < 100 ; i++)
// {
//     xx[0] = (i -50) * 0.08;
//     xx[1] = 0.;
//     xy[i][0] = (i -50) * 0.04;
//     xy[i][1] = fonction(xx);
// }
//
// Savetxt("test.txt", xy, 100, 2);
//
// for (i = 0;  i < 100 ; i++) free(xy[i]);
// free(xy);
// return 1;

    // Proposal distribution
    gw10_distribution distribg(parama);
    // ------------ Estimate variables ( uncomment lines related to these variables to enable jumps conditionned to estimates
//     fctFittriple_gaussprior fonctionestimate = fctFittriple_gaussprior();
//     fonctionestimate.verbose=0;
//     char datafile_estimate[500];
//     double minestimate =0.;
//     sprintf(datafile_estimate, "%s-estimate", datafile);

// Initial distribution
//     normal_distribution<double> initial_distribution(0.,0.0001); // used for initialization

// Posterior distribution
// fctmultigaussienne fonction = fctmultigaussienne(20,1);

// To run with a dummy model and simulated data
#ifdef SIMUDATA
    fctsimudata fonction =fctsimudata();
    const int simundim = 10;
    double simuparams[simundim];
    int simundata = 1000;
    double simustddev = 0.5;
    double simubasefreq = 2.;
    double simutimespan = 10.;
    double ** simucoeffs;
    double * simudata;

    simudata =(double*) malloc(sizeof(double)*simundata);
    simucoeffs = (double**) malloc(sizeof(double*) * simundata);
    for (i = 0; i < simundata ; i++) simucoeffs[i] = (double*) malloc(sizeof(double)*simundim);

    for (i = 0 ; i < simundim ; i ++) simuparams[i] = 8. - i;
    if (mpirank == 0)
    {
      fonction.generate_data(simundata,  simuparams, simustddev, simundim, simubasefreq,  simutimespan); // int numberofdata,  double * parameters, double stddev, int ndimensions, double basefreq, double timespan)
      for ( i = 0 ; i<simundata ; i ++)
      {
        simudata[i] = fonction.data[i];
        for ( j = 0 ; j<simundim ; j ++)
        {
          simucoeffs[i][j] = fonction.coeffs[i][j];
        }
      }
    }
  #ifdef MPIMODE
      code = MPI_Bcast(simudata, simundata, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for (i = 0 ; i < simundata; i++)
      {
        code = MPI_Bcast(simucoeffs[i], simundim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
      if (mpirank > 0) fonction.external_data(simudata, simucoeffs, simundata, simundim, simustddev);
  #endif


   for (i = 0; i < simundata ; i++) free(simucoeffs[i]);
   free(simucoeffs);
   free(simudata);
#endif
// End of dummy model and simulated data

#ifndef SIMUDATA
 fctFittriple_gaussprior fonction = fctFittriple_gaussprior();
 fonction.verbose=1;
 // turntype * turns;
 // long int nturns=0;

#ifdef MPIMODE
       code = MPI_Barrier(MPI_COMM_WORLD);
#endif
       i = 0;
       if (i == mpirank%coeursparnoeud)
       {
           fonction.Reconstruct(parfile, datafile, true);
           if (turnfile != NULL )
           {
               fonction.Load_turn_numbers(turnfile);
               if (mpirank == masterproc) printf("\n Turns loaded from : %s \n\n", turnfile);
           };
//            fonctionestimate.Reconstruct(parfile, datafile_estimate);
           printf("\n ---- > passed slow intitialisation %i \n", mpirank);
           // nturns = fonction.ntoas;
           // memcpy(turns, fonction.turns, nturns*sizeof(turntype));
       }

       #ifdef MPIMODE
              code = MPI_Barrier(MPI_COMM_WORLD);
       #endif

       for (i = 1 ; i < coeursparnoeud ; i++) // initialization simultaneous on each node but not on each thread because of I/O with files
       {
         if (i == mpirank%coeursparnoeud)
         {
             fonction.Reconstruct_noCompute(parfile, datafile, true);
             printf("\n ---- > passed fast intitialisation %i \n", mpirank);
         }
       }
#ifdef MPIMODE
           code = MPI_Bcast(fonction.turns, fonction.ntoas, MPI_LONG_LONG_INT, masterproc, MPI_COMM_WORLD);
#endif
   fonction.Set_tracker_off();
#endif

 ndim = fonction.Get_ndim();   // sevt ndim
#ifdef MPIMODE
   code = MPI_Barrier ( MPI_COMM_WORLD); // Not really necessary but synchronizes the prints
#endif
   //************************************************************************************************************************************

 // Intitializing variables depending on ndim (in case fonction sets ndim)

   // Printing
   size_strtoprint = max(size_strtoprint * ndim, (ndim+1)*100 * nbw_per_proc * nbproc*ntemp);
   strtoprint = (char*) malloc(sizeof(char)*size_strtoprint);
   strcpy(strtoprint, "");

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

// Initializes temperature variables
    betas = (double*) malloc(ntemp * sizeof(double));
    Temp_accepted = (double*) malloc(ntemp * sizeof(double));
    Temp_acceptance_fraction = (double*) malloc(ntemp * sizeof(double));
    accepted_tempswaps = (double*) malloc((ntemp-1) * sizeof(double));
    last_accepted_tempswaps = (double*) malloc((ntemp-1) * sizeof(double));
    walkerbuffer = (double*) malloc(ndim * sizeof(double));
    Temp_lnchi2 = (double**) malloc(ntemp * sizeof(double*));
    Temp_walkers1d = (double**) malloc(ntemp * sizeof(double*));
    if (record_refused == true)
    {
        Temp_irefused = (int*) malloc(ntemp * sizeof(int));
    }
    for (itemp = 0 ; itemp < ntemp; itemp++)
    {
        Temp_lnchi2[itemp] = (double*) malloc(2* halfnbw * sizeof(double));
        Temp_walkers1d[itemp] = (double*) malloc(2*halfnbw * ndim * sizeof(double));
        Temp_accepted[itemp] = 0.;
        Temp_acceptance_fraction[itemp] = 0.;
        if (record_refused == true)  Temp_irefused[itemp] = 0;
    };
    betas[0] = 1./targettemp;
    for (itemp=1 ; itemp < ntemp ; itemp++)
    {
        betas[itemp] = betas[itemp-1]/tempfactor[itemp-1];
        accepted_tempswaps[itemp-1] = 0.;
        last_accepted_tempswaps[itemp -1] = 0.;
    }



 // size of the chain
    rem = fmod(static_cast<double>(maxit), static_cast<double>(chain_freq) );
    chain_size = int ( (static_cast<double>(maxit) - rem) / static_cast<double>(chain_freq)   + 1) * 2 * halfnbw  ;

// Initializes statistical diagnostics
    if (mpirankincluster == 0)
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
        autocorrelation = (double*) malloc(ensemble_chain_size  * sizeof(double) );
    }

// Initialises MCMC diagnostics
  if (dodiagnostics)
  {
    sprintf(diagnosticsfile, "%s-diagnostics-proc%d.txt", fileout, mpirank);
    sprintf(walkersdiagnosticsfilename, "%s-walkers_diagnostics-proc%d.txt", fileout, mpirank);
    ndiagnostics = maxit * ntemp * 2 * halfnbw * 2;
    diagnosticslog = (double**) malloc( ndiagnostics * sizeof(double*) );
    for (k = 0 ; k < ndiagnostics ; k ++)
    {
      diagnosticslog[k] = (double*) malloc( 10 * sizeof(double) );
    }
  }

// Initialises sieving log
  ncol_sievinglog = ndim + 2;
  if (maxsieved<0) maxsieved = halfnbw;
  sievinglog = (double **) malloc(maxsieved*sizeof(double*));
  for ( k = 0 ; k < maxsieved ; k++) sievinglog[k] = (double * ) malloc(ncol_sievinglog * sizeof(double));

  // Walkers memory allocations
   walkers = (double **) malloc(2*halfnbw * sizeof(double*)) ;
   walkers1d = (double*) malloc(2*halfnbw * ndim * sizeof(double)) ;
   walkers1d_local = (double*) malloc((1+int(dosieving))*nbw_per_proc * ndim * sizeof(double) );
   for ( k = 0 ; k < 2*halfnbw ; k++) walkers[k] = &walkers1d[k*ndim];


// Initialization of the walkers
   for (itemp = 0 ; itemp < ntemp ; itemp++)
   {
    if (mpirank == masterproc)
    {
        if ( prev_chain_file != NULL ) //  from a previous chain
            {
                sprintf(fileouttemp, "%s--temp=%.2f.res", prev_chain_file, 1./betas[itemp]);
                printf("\n **** Initialization from previous chain in file %s \n", fileouttemp);
                printf("\n Temperature %.5e \n", 1./betas[itemp]);
                if (Init_from_prev_chain(fileouttemp, walkers1d, Temp_lnchi2[itemp], ndim, 2*halfnbw) == 1)
                {
                    printf("Error : Initialization failed ! ");
                    return 1;
                }
        /* ! Test !
                printf("\n\n test temperature %.5e\n", 1./betas[itemp]);
                for (k = 0 ; k < 2* halfnbw; k++)

             printf("ffwalker%i | ", k);
             for (i = 0; i < ndim ; i++)  printf("%.4e ", walkers[k][i]);
             printf(" | %.5e \n", lnchi2[k]);
        */
       printf("\n");
            }
            else // in a random gaussian ball
            {
                for ( i = 0 ; i < ndim*2*halfnbw ; i++) walkers1d[i] = fonction.initial_distribution(randomgen, i%(ndim));
            }
    };
    #ifdef MPIMODE
    code = MPI_Bcast(walkers1d, ndim*2*halfnbw, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
    code = MPI_Bcast(Temp_lnchi2[itemp], 2*halfnbw, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
    #endif
    memcpy(Temp_walkers1d[itemp], walkers1d, ndim*2*halfnbw*sizeof(double));
    //printf("\n walkers 1d 1\n");
    //Print_table(walkers,2*halfnbw, ndim);
   }

   // initialization of chaine
   if (mpirankincluster==0 )
   {
       chaine = (double ***) malloc(ntemp *sizeof(double**)) ;
       for (itemp = 0 ; itemp < ntemp ; itemp++)
       {
            chaine[itemp] = (double **) malloc(chain_size *sizeof(double*)) ;
            for (k = 0 ; k < chain_size ; k++) chaine[itemp][k] = (double*) malloc((ndim + 1) * sizeof(double));
       };
   }

   if (record_refused == true )
   {
        irefuseds = (int*) malloc(sizeof(int) * nbproc);
        chaine_refused = (double ***) malloc(ntemp *sizeof(double**)) ;
        for (itemp = 0 ; itemp < ntemp ; itemp++)
        {
            chaine_refused[itemp] = (double **) malloc(chain_size *sizeof(double*)) ;
            for (k = 0 ; k < chain_size ; k++) chaine_refused[itemp][k] = (double*) malloc((ndim + 1) * sizeof(double));
        };
    }

// end of Memory alocations


   if (mpirank == masterproc) {
   // printf("\n!!!!!!!!!!! EXPERIMENTAL SMOOTHING at level %.5e \n", smoothinglevel); // Experimental test !!
   printf("\n*****************************************************************************************************\n");
   printf(" Function parfile= %s and datafile= %s \n", parfile, datafile);
   fonction.Print_MCMC();
//    printf(" Function estimate parfile= %s and datafile= %s \n", parfile, datafile_estimate);
//    fonctionestimate.Print_MCMC();
   //printf(" After initialization (and before MCMC) the chi2 is : %.5f\n", fonction.Compute_lnposterior(0) );
   printf(" Parameter of the proposal distribution a=%.5e \n",parama);
   printf(" Temperatures : ");
   for (itemp=0; itemp < ntemp ; itemp++) printf("%.5e,",1./betas[itemp]);
   printf("\n");
   printf(" Swap temperatures every %i iterations.\n ",tempswapfreq);
   printf(" Number of dimensions : %i \n", ndim );
   printf(" There are %i walkers per temperature spread on %i processors per cluster. \n", 2*halfnbw, nproc_per_temp);
   printf(" There are %i clusters each treating %i temperatures. \n", ncluster, ntemp_per_cluster);
   printf(" The whole set of walkers will undergo %i moves.\n", maxit);
   printf(" The state will of the set of walkers will be saved every %i moves, for a total of %i walkers saved.\n", chain_freq, chain_size);
   printf(" The length of the ensemble chain is thus %i \n", ensemble_chain_size);
   printf(" The resulting chain + lnchi2 will be saved with prefix %s every %i iterations. \n", fileout, moy_freq);
   printf("*******************************************************************************************************\n\n");



    printf("\n Initializing the walkers with chi2. \n\n");
   }

    for (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++) //(itemp=0; itemp < ntemp ; itemp++)
    {
        memcpy(walkers1d, Temp_walkers1d[itemp], 2* sizeof(double) * halfnbw * ndim) ;
        //printf("\n walkers 1d 2\n");
        //Print_table(walkers,2*halfnbw, ndim);

        for ( k = mpirankincluster * nbw_per_proc ; k < ( mpirankincluster + 1) *nbw_per_proc ; k++)
        {
            lnchi2_local[k- mpirankincluster * nbw_per_proc] = fonction(walkers[k], betas[itemp]);
            lnchi2_local[(k- mpirankincluster * nbw_per_proc) + nbw_per_proc] = fonction(walkers[k + halfnbw], betas[itemp]);

            if ( prev_chain_file != NULL )  // Check that the previous chain had the same lnposteriors as the freshly calculated ones
            {
              if (abs(lnchi2_local[k- mpirankincluster * nbw_per_proc] -  Temp_lnchi2[itemp][k]) > 1.e-3)
                printf("\nWarning !(1) lnposterior from previous chain file (=%.15e) does not match the calculated one (=%.15e) for walker %i in temperature %i on process %i\n", Temp_lnchi2[itemp][k], lnchi2_local[k- mpirankincluster * nbw_per_proc], k, itemp, mpirank);
              if (abs(lnchi2_local[k- mpirankincluster * nbw_per_proc + nbw_per_proc] -  Temp_lnchi2[itemp][k + halfnbw]) > 1.e-3)
                printf("\nWarning !(2) lnposterior from previous chain file (=%.15e) does not match the calculated one (=%.15e) for walker %i in temperature %i on process %i\n", Temp_lnchi2[itemp][k + halfnbw], lnchi2_local[k- mpirankincluster * nbw_per_proc + nbw_per_proc], k + halfnbw, itemp, mpirank);
            }
        }
    #ifdef MPIMODE
        code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, Temp_lnchi2[itemp], nbw_per_proc, MPI_DOUBLE, MPI_COMM_tempcluster);
        code = MPI_Allgather(lnchi2_local + nbw_per_proc , nbw_per_proc, MPI_DOUBLE, Temp_lnchi2[itemp] + halfnbw, nbw_per_proc, MPI_DOUBLE, MPI_COMM_tempcluster);
    #else
        memcpy(Temp_lnchi2[itemp], lnchi2_local, 2* sizeof(double)* nbw_per_proc);
    #endif
    }

  if (mpirankincluster==0 ) // Print initial state of the chain
   {
       strcpy(strtoprint,"");
       sprintf(strinter, "\n *** mpirank  %i %i %i*** \n\n", mpirank, icluster, ntemp_per_cluster);
       strcat(strtoprint, strinter);
     for  (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++)//(itemp=0; itemp < ntemp ; itemp++)
     {
       sprintf(strinter, "\n*** Initial chain state, temperature =%.5e : \n", 1./betas[itemp]);
       strcat(strtoprint, strinter);

       memcpy(walkers1d, Temp_walkers1d[itemp], 2* sizeof(double) * halfnbw * ndim) ;
       memcpy(lnchi2, Temp_lnchi2[itemp], 2* halfnbw * sizeof(double));

       for (k = 0 ; k < 2* halfnbw; k++)
        {
            sprintf(strinter, "walker%i | ", k);
            strcat(strtoprint, strinter);
            for (i = 0; i < ndim ; i++)
            {
                sprintf(strinter, "%.4e ", walkers[k][i]);
                strcat(strtoprint, strinter);
            }
            sprintf(strinter," | %.5e \n", lnchi2[k]);
            strcat(strtoprint, strinter);
        }
       sprintf(strinter,"\n");
       strcat(strtoprint, strinter);
     }

   }

#ifdef MPIMODE
   MPI_Collect_and_print(strtoprint, size_strtoprint, masterproc, printers, ncluster, MPI_COMM_WORLD);
#else
   printf("%s", strtoprint);
#endif

   fonction.verbose=0; // minimal printing
//    fonctionestimate.verbose=0;

    if (mpirank == masterproc)
    {
      printf("\n Running time: %.0f seconds.\n", ((float)(clock()  - tclock_start)/CLOCKS_PER_SEC));
      printf("###########################################################################################");
      printf("\n                      Starting the actual MCMC... \n\n");
      printf("###########################################################################################");
    }

//############################################################################################
 // **********   Main loop
for (int it = 0 ; it < maxit ; it++)
{
    printstat = fmod(static_cast<double>(it+1 ) , static_cast<double>(moy_freq)); // if = 0 stats will be printed during this iteration
    strcpy(strtoprint,"");

// Send all the results of each temperature cluster to masterproc so that it has them all before swapping
 #ifdef MPIMODE
    if (mpirankincluster == 0 and mpirank != masterproc)
    {
//             printf("\n *** mpirank  %i %i %i*** \n\n", mpirank, icluster, ntemp_per_cluster);
        for (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++)
        {
//               printf("\n *** Sending itemp %i on proc %i \n\n", itemp, mpirank);
             code = MPI_Send(Temp_walkers1d[itemp], 2* halfnbw * ndim, MPI_DOUBLE, masterproc, itemp, MPI_COMM_WORLD);
             code = MPI_Send(Temp_lnchi2[itemp], 2* halfnbw , MPI_DOUBLE, masterproc, itemp, MPI_COMM_WORLD);
        }
    }
    if (mpirank == masterproc)
    {
        for (k = 0; k < ncluster ; k++)
        {   if (nproc_per_temp* k != masterproc)
            {
            for (itemp = k * ntemp_per_cluster; itemp < min((k+1) * ntemp_per_cluster, ntemp) ; itemp++)
            {

//               printf("\n *** Receiving itemp %i on proc %i from proc %i \n\n", itemp, mpirank, nproc_per_temp * k);
                 code = MPI_Recv(Temp_walkers1d[itemp], 2* halfnbw * ndim, MPI_DOUBLE, nproc_per_temp * k, itemp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 code = MPI_Recv(Temp_lnchi2[itemp], 2* halfnbw , MPI_DOUBLE, nproc_per_temp * k, itemp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            }
        }
    }
#endif

// Do temperature swaps
    if (mpirank == masterproc)
    {
       // Swap every tempswapfreq iterations
        rem = fmod(static_cast<double>(it+1 ) , static_cast<double>(tempswapfreq));
       if (rem == 0.)
        {
        if (ntemp > 1) printf("\n *** Swapping temperatures *** \n\n");
        kshift += 1;
        for (itemp = 0 ; itemp < ntemp -1; itemp++)
        {
            for  ( k= 0 ; k < 2*halfnbw ; k++)
            {
                kex =(k + kshift) % (2 * halfnbw);

                accept = fonction.Get_chi2(Temp_lnchi2[itemp][k], &Temp_walkers1d[itemp][k*ndim], betas[itemp] );
                accept -= fonction.Get_chi2(Temp_lnchi2[itemp+1][kex], &Temp_walkers1d[itemp + 1][kex*ndim], betas[itemp+1] );
                accept *= (betas[itemp] - betas[itemp+1]);

//                 accept = (Temp_lnchi2[itemp+1][kex]/betas[itemp+1] - Temp_lnchi2[itemp][k]/betas[itemp]) * (betas[itemp] - betas[itemp+1]); // ln of the acceptance probability. divide by betas[itemp(+1)] because temperature included in lnchi2

                if (accept < 0.)
                {
                    accept = exp(accept);
                    accepted_tempswaps[itemp] += accept;
                    last_accepted_tempswaps[itemp] += accept;
                    if (random01(randomgen) <= accept)
                    {
                        accept = 1.;
                    }
                }
                else
                {
                    accepted_tempswaps[itemp] += 1.;
                    last_accepted_tempswaps[itemp] += 1.;
                    accept = 1.;
                };
                if (accept == 1.) // swap the walkers
                {
                    fonction.Swap_temperatures(&Temp_walkers1d[itemp][k * ndim], Temp_lnchi2[itemp][k], betas[itemp], &Temp_walkers1d[itemp+1][kex * ndim], Temp_lnchi2[itemp+1][kex], betas[itemp+1] );

                    /*lnchi2y = Temp_lnchi2[itemp][k];
                    Temp_lnchi2[itemp][k] = Temp_lnchi2[itemp+1][kex] * betas[itemp] / betas[itemp+1];
                    Temp_lnchi2[itemp+1][kex] = lnchi2y * betas[itemp+1] / betas[itemp];;

                    memcpy(walkerbuffer, &Temp_walkers1d[itemp][k * ndim], ndim*sizeof(double));
                    memcpy(&Temp_walkers1d[itemp][k*ndim], &Temp_walkers1d[itemp+1][kex*ndim], ndim*sizeof(double));
                    memcpy(&Temp_walkers1d[itemp+1][kex*ndim], walkerbuffer, ndim*sizeof(double));
                */};
            };
        };
        };
    };
// Broadcast back the result of the swap. Here one could also do a send/receive to each submaster (mpirankincluster=0).
#ifdef MPIMODE
        for (itemp = 0 ; itemp < ntemp -1; itemp++)
        {
            code = MPI_Bcast(Temp_lnchi2[itemp], 2*halfnbw, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);
            code = MPI_Bcast(Temp_walkers1d[itemp], ndim*2*halfnbw, MPI_DOUBLE, masterproc, MPI_COMM_WORLD);

        };
#endif
    // End of temperature swaps

    for(int itemp = icluster * ntemp_per_cluster ; itemp < min((icluster+1) * ntemp_per_cluster, ntemp) ; itemp++)
    { // Loop over temperatures

        // Select the temperature chain to advance
         memcpy(walkers1d, Temp_walkers1d[itemp], 2* sizeof(double) * halfnbw * ndim) ;
         memcpy(lnchi2, Temp_lnchi2[itemp], 2* halfnbw * sizeof(double));
         accepted = Temp_accepted[itemp];
         acceptance_fraction = Temp_acceptance_fraction[itemp];
         if (record_refused == true) irefused = Temp_irefused[itemp];

#ifdef MPIMODE
        code = MPI_Barrier ( MPI_COMM_tempcluster);
#endif

        for ( k = mpirankincluster * nbw_per_proc ; k < ( mpirankincluster + 1) *nbw_per_proc ; k++)
        {
            // Take care of walker "k"
            ntry = 0;
            do
            {
                // Draw "z" for stretch move and newwalker number wj
                wj = randomwalker(randomgen) + halfnbw;
                // G = 2 /sqrt (z)
                z = distribg(random01(randomgen)) ;
                // Determine walkery the possible new position of walker "k"
                for (i=0; i < ndim ; i++) walkery[i] = walkers[wj][i] + z * ( walkers[k][i] - walkers[wj][i] ) ;
                
                // Check that this walker is compatible with priors
                ntry += 1;
            } while (fonction.is_valid(walkery) == false and ntry < 100);
            if (ntry >= 100) printf("\n Warning : No valid walker could be found in 100 trials \n\n");

            // Make an estimate of the new log likelyhood of the new walker
//             lnchi2y = fonctionestimate(walkery) * betas[itemp];
//             } while (lnchi2y < minestimate);
            // Compute p(Y) and compute acceptance accept

            // Print_walker(walkers[k], lnchi2[k], ndim, strinter);
            // printf("proc%i1 wcurrent | %s || z=%.2f\n",mpirank, strinter,z);

            // Print_walker(walkery, 0., ndim, strinter);
            // printf("<<proc%i1 w | %s\n",mpirank, strinter);

            lnchi2y = fonction(walkery, betas[itemp]);
            // Print_walker(walkery, lnchi2y, ndim, strinter);
            // printf(">>proc%i1 w | %s\n",mpirank, strinter);
            // Determine if the step is accepted by computing accept = z^(ndim-1)*p(Y) / p(X_k) = pow(z, ndim - 1)*exp(lnchi2y - lnchi2[k])
            dlnchi2 = lnchi2y - lnchi2[k];
            // Experimental ! Test !
            // if (abs(dlnchi2) < smoothinglevel)
            // {
            //   dlnchi2 = 0.;
            //   countersmoothing +=1;
            // }
            // End Experimental
            mlnznm1 = -(ndim-1)*log(z);
            if (dlnchi2 > mlnznm1) // necessarily accepted since i this case accept > 1
                accept = 2.;
            else
                accept = exp(-mlnznm1 + dlnchi2);

            // Record MCMC algorithm diagnostics
            if (dodiagnostics)
            {
              diagnosticslog[idiagnostics][0] = double(it);
              diagnosticslog[idiagnostics][1] = double(itemp);
              diagnosticslog[idiagnostics][2] = 0.;
              diagnosticslog[idiagnostics][3] = double(k);
              diagnosticslog[idiagnostics][4] = double(wj);
              diagnosticslog[idiagnostics][5] = z;
              diagnosticslog[idiagnostics][6] = lnchi2[k];
              diagnosticslog[idiagnostics][7] = lnchi2y;
              diagnosticslog[idiagnostics][8] = accept;
              diagnosticslog[idiagnostics][9] = 0.;          // see below
              idiagnostics +=1;
            }

            // draw r in [0,1] to accept or reject the move
            if (random01(randomgen) <= accept)
            {
                for (i=0; i < ndim ; i++) walkers[k][i] = walkery[i];
                lnchi2[k] = lnchi2y;
                accepted += 1.;
                if (dodiagnostics) diagnosticslog[idiagnostics-1][9] = 1.;
            } else
            {
                if (record_refused == true)
                {
                    for (i = 0; i < ndim ; i++) chaine_refused[itemp][irefused][i] = walkery[i];
                    chaine_refused[itemp][irefused][ndim] = lnchi2y;
                    irefused++;
                }
            }
        }
#ifdef MPIMODE
        memcpy(lnchi2_local, &lnchi2[mpirankincluster * nbw_per_proc] , nbw_per_proc*sizeof(double) );
        memcpy(walkers1d_local, &walkers1d[(mpirankincluster * nbw_per_proc )*ndim], nbw_per_proc*ndim * sizeof(double) );
        for(i = 0; i < nbw_per_proc ; i++)
        {
            if (lnchi2_local[i] != lnchi2[mpirankincluster * nbw_per_proc + i]) printf("\n\nCATASTROPHE !! \n\n");
        }
        code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, lnchi2, nbw_per_proc, MPI_DOUBLE, MPI_COMM_tempcluster);
        code = MPI_Allgather(walkers1d_local, nbw_per_proc*ndim, MPI_DOUBLE, walkers1d, nbw_per_proc*ndim, MPI_DOUBLE, MPI_COMM_tempcluster);
#endif

        for ( k = mpirankincluster * nbw_per_proc + halfnbw ; k < ( mpirankincluster + 1) *nbw_per_proc + halfnbw ; k++)
        {
            // Take care of walker "k"
            ntry = 0;
            do {
                // Draw "z" for stretch move and newwalker number wj
                wj = randomwalker(randomgen) ;
                // G = 2 /sqrt (z)
                z = distribg(random01(randomgen)) ;//0.25 * pow(Ng *(random01(randomgen) - cteg), 2 );
                // Determine walkery the possible new position of walker "k"
                for (i=0; i < ndim ; i++) walkery[i] = walkers[wj][i] + z * ( walkers[k][i] - walkers[wj][i] ) ;
            
                // Check that this walker is compatible with priors
                ntry += 1;
            } while (fonction.is_valid(walkery) == false and ntry < 100);
            if (ntry >= 100) printf("\n Warning : No valid walker could be found in 100 trials \n\n");
            // Make an estimate of the new log likelyhood of the new walker
//             lnchi2y = fonctionestimate(walkery) * betas[itemp];
//             } while (lnchi2y < minestimate);
            // Compute ln( p(Y) ) and compute acceptance accept
            // Print_walker(walkers[k], lnchi2[k], ndim, strinter);
            // printf("proc%i2 wcurrent | %s || z=%.2f\n",mpirank, strinter,z);

            // Print_walker(walkery, 0., ndim, strinter);
            // printf("<<proc%i2 w | %s\n",mpirank,strinter);
            lnchi2y = fonction(walkery, betas[itemp]);
            // Print_walker(walkery, lnchi2y, ndim, strinter);
            // printf(">>proc%i2 w | %s\n",mpirank,strinter);
            // Determine if the step is accepted by computing accept = z^(ndim-1)*p(Y) / p(X_k) = z^(ndim-1)*exp(lnchi2y - lnchi2[k])
            dlnchi2 = lnchi2y - lnchi2[k];
            // Experimental ! Test !
            // if (abs(dlnchi2) < smoothinglevel)
            // {
            //   dlnchi2 = 0.;
            //   countersmoothing +=1;
            // }
            // End Experimental
            mlnznm1 = -(ndim-1)*log(z);
            if (dlnchi2 > mlnznm1) // necessarily accepted since i this case accept > 1
                accept = 2.;
            else
                accept = exp(-mlnznm1 + dlnchi2); //  pow(z, ndim - 1) * exp(lnchi2y - lnchi2[k]) ;

          // Record MCMC algorithm diagnostics
            if (dodiagnostics)
            {
              diagnosticslog[idiagnostics][0] = double(it);
              diagnosticslog[idiagnostics][1] = double(itemp);
              diagnosticslog[idiagnostics][2] = 1.;
              diagnosticslog[idiagnostics][3] = double(k);
              diagnosticslog[idiagnostics][4] = double(wj);
              diagnosticslog[idiagnostics][5] = z;
              diagnosticslog[idiagnostics][6] = lnchi2[k];
              diagnosticslog[idiagnostics][7] = lnchi2y;
              diagnosticslog[idiagnostics][8] = accept;
              diagnosticslog[idiagnostics][9] = 0.;          // see below
              idiagnostics +=1;
            }

            // draw r in [0,1] to accept or reject the move
            if (random01(randomgen) <= accept)
            {
                for (i=0; i < ndim ; i++) walkers[k][i] = walkery[i];
                lnchi2[k] = lnchi2y;
                accepted += 1.;
                if (dodiagnostics) diagnosticslog[idiagnostics - 1][9] = 1.;
            }  else
            {
                if (record_refused == true)
                {
                    for (i = 0; i < ndim ; i++) chaine_refused[itemp][irefused][i] = walkery[i];
                    chaine_refused[itemp][irefused][ndim] = lnchi2y;
                    irefused++;
                }
           }
        }

#ifdef MPIMODE
        memcpy(lnchi2_local, &lnchi2[mpirankincluster * nbw_per_proc + halfnbw] , nbw_per_proc*sizeof(double) );
        memcpy(walkers1d_local, &walkers1d[(mpirankincluster * nbw_per_proc + halfnbw)*ndim], nbw_per_proc*ndim * sizeof(double) );

        code = MPI_Allgather(lnchi2_local , nbw_per_proc, MPI_DOUBLE, &lnchi2[halfnbw], nbw_per_proc, MPI_DOUBLE, MPI_COMM_tempcluster);
        code = MPI_Allgather(walkers1d_local, nbw_per_proc*ndim, MPI_DOUBLE, &walkers1d[halfnbw*ndim], nbw_per_proc*ndim, MPI_DOUBLE, MPI_COMM_tempcluster);
#endif



// Record the state of the chain, do diagnostics..
       if (mpirankincluster == 0)
       {
           rem = fmod(static_cast<double>(it + 1), static_cast<double>(chain_freq));
           ensemble_chain_index = int( (static_cast<double>(it + 1) - rem ) / static_cast<double>(chain_freq) );

           if (rem == 0. or it == 0)
           {
                if (chain_freq == 1)  ensemble_chain_index -= 1; // ensemble_chain_index = nb of chain_freq done, but index of array starts at 0;
//                printf("Saving chain after %i iterations. %i walkers saved.  \n", it, ensemble_chain_index * 2 * halfnbw);
                for(k=0; k < 2*halfnbw ; k++)
                {
                    for (i = 0; i < ndim ; i++) chaine[itemp][ensemble_chain_index*2*halfnbw + k][i] = walkers[k][i];
                    chaine[itemp][ensemble_chain_index*2*halfnbw + k][ndim] = lnchi2[k];
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
           if (printstat ==0.)
           {

      // Save MCMC diagnostics
             if (dodiagnostics)
             {
               Savetxt(diagnosticsfile, diagnosticslog, idiagnostics, diagnostics_columns);
               //Savetxt(walkersdiagnosticsfile, walkers, 2*halfnbw, ndim);
               walkersdiagnosticsfile = fopen(walkersdiagnosticsfilename, "w") ;
               for(k=0 ; k < 2*halfnbw ; k++)
               {
                 for (i = 0; i < ndim ; i++)  fprintf(walkersdiagnosticsfile, "%.15e    ", walkers[k][i] );
                 fprintf(walkersdiagnosticsfile, "%.15e\n", lnchi2[k]);
               }
               fclose(walkersdiagnosticsfile);
             }
      // End of MCMC diagnostics


      // Do the sieving of the excessively low lnposteriors
      // TODO : This action should probably have its own periodicity and not be attached to diagnostics
             if (dosieving)
             {
               maxlnchi2index = Maxofarray(lnchi2, 2*halfnbw);
               maxlnchi2 = lnchi2[maxlnchi2index];
               isieved = 0;

               k = 0;
               while(k < 2 * halfnbw and isieved < maxsieved)
               {
                 if (maxlnchi2 - lnchi2[k] > maxtension) // Save up to the sieved walker log.
                 {
                   memcpy(sievinglog[isieved], walkers[k], sizeof(double)*ndim);
                   sievinglog[isieved][ndim] = lnchi2[k];
                   sievinglog[isieved][ndim+1] = maxlnchi2 - lnchi2[k];
                   isieved ++;
                 }
                 k++;
               }
               if (mpirankincluster == 0) printf("\n* Now sieving... (maxlnposterior = lnchi2[%i]= %.5e, maxtension = %.5e)\n", maxlnchi2index, maxlnchi2, maxtension);
               if (isieved >= maxsieved)
               {
                 if (mpirankincluster == 0) printf(" More than maxsieved=%i walkers to sieve. Aborting...\n", maxsieved);
               }
               else
               {
                 for (k = mpirankincluster * nbw_per_proc * 2; k < (mpirankincluster+1) * nbw_per_proc * 2 ; k ++)
                 {
                   j =0;
                   lnchi2y = lnchi2[k];
                   while (maxlnchi2 - lnchi2y > maxtension and j < max_sieving_attempts)
                   {
                     j++;
                     wj = randomwalker(randomgen) ;
                     while (maxlnchi2 - lnchi2[wj] > maxtension ) wj = randomwalker(randomgen);
                     wk = halfnbw + randomwalker(randomgen) ;
                     while (maxlnchi2 - lnchi2[wk] > maxtension or wk == wj ) wk = randomwalker(randomgen);
                     // Propose replacement walker. The swap and .2 factor are empirical ways of ensuring the proposal will be successful.
                     ntry=0;
                     do {
                        z = distribg(random01(randomgen)) ;//0.25 * pow(Ng *(random01(randomgen) - cteg), 2 );
                        if (lnchi2[wj] > lnchi2[wk]) Swap(wk,wj); // ie wk = wj and wj = wk
                        for (i=0; i < ndim ; i++) walkery[i] = walkers[wk][i] + 0.2*z * ( walkers[wk][i] - walkers[wj][i] ) ;
                        ntry += 1;
                     } while (fonction.is_valid(walkery) == false and ntry < 100);
                     if (ntry >= 100) printf("\n Warning : No valid walker could be found in 100 trials \n\n");
                     lnchi2y = fonction(walkery, betas[itemp]);
                     printf("   Rep candidate for lnchi2[%i]=%.5e on proc%i: lnchi2y=%.5e from lnchi2[%i]=%.5e and lnchi2[%i]=%.5e\n", k, lnchi2[k], mpirank,lnchi2y, wk, lnchi2[wk], wj, lnchi2[wj]);

                   }
                   if (j>0)
                   {
                     if (j == max_sieving_attempts)
                     {
                        printf("\n WARNING : no satisfying replacement was found for sieving in %i attempts. Is maxtension set too low ? mpirank %i\n\n",max_sieving_attempts, mpirank);
                     }
                     else
                     {
                        lnchi2[k] = lnchi2y;
                        memcpy(walkers[k], walkery,ndim*sizeof(double));
                     }
                  }
                 }
                 if (mpirankincluster == 0)
                 {
                   printf("Sieving done : %i walkers were sieved.\n", isieved);
                   sprintf(fileouttemp, "%s-Sieving_log--temp=%.2f.res", fileout, 1./betas[itemp]);
                   Appendtxt(fileouttemp, sievinglog, isieved, ncol_sievinglog);
                 }
                 #ifdef MPIMODE
                  memcpy(lnchi2_local, &lnchi2[mpirankincluster * 2 * nbw_per_proc] , 2*nbw_per_proc*sizeof(double) );
                  memcpy(walkers1d_local, &walkers1d[(mpirankincluster * 2 * nbw_per_proc)*ndim], 2*nbw_per_proc*ndim * sizeof(double) );

                  code = MPI_Allgather(lnchi2_local , 2*nbw_per_proc, MPI_DOUBLE, lnchi2, 2*nbw_per_proc, MPI_DOUBLE, MPI_COMM_tempcluster);
                  code = MPI_Allgather(walkers1d_local, 2*nbw_per_proc*ndim, MPI_DOUBLE, walkers1d, 2*nbw_per_proc*ndim, MPI_DOUBLE, MPI_COMM_tempcluster);

                 #endif
                 // MPI check
              //   printf("\n\nBRAAAAAAAAAAAAAAAAAAAATTTTT\n");
                 // for (k = mpirankincluster * nbw_per_proc * 2; k < (mpirankincluster+1) * nbw_per_proc * 2 ; k ++)
                 // {
                 //   lnchi2y = fonction(walkers[k], betas[itemp]);
                 //   if (lnchi2y != lnchi2[k]) printf("\nMPI PRONLEM 1 %i %i %.6e != %.6e!!! \n\n ", k, mpirank, lnchi2y, lnchi2[k]); //(lnchi2_local[k-mpirankincluster * nbw_per_proc * 2] != lnchi2[k]) printf("\nMPI PRONLEM 1 !!! \n\n ");
                 //   for (i = 0; i < ndim ; i++)
                 //   {
                 //     if (walkers1d_local[ndim*(k-mpirankincluster * nbw_per_proc * 2) + i] != walkers1d[k*ndim + i]) printf("\nMPI PRONLEM 2 !!! \n\n ");
                 //   }
                 // }
               }
             }
      // End of sieving


#ifdef MPIMODE
               code = MPI_Reduce(&accepted , &total_accepted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_tempcluster); // Compute total number of accepted moves
 // Experimental test !!
               // code = MPI_Reduce(&countersmoothing , &total_countersmoothing, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_tempcluster); // Compute total number of smoothed jumps
#else
               total_accepted = accepted;
// Experimental test !!
               // total_countersmoothing =countersmoothing;
#endif
                if (mpirankincluster == 0)//(mpirank == masterproc)
                {

                sprintf(strinter,"* Temperature : %.5e |  Ensemble chain step : %i \n", 1./betas[itemp], ensemble_chain_index);
                strcat(strtoprint, strinter);
                for (i = 0 ; i < ndim ; i++)
                {
                    for (k = 0 ; k < (ensemble_chain_index+1) ; k++)
                    {   // Create an ensemble chain with each parameter averaged over the set of walkers
                        ensemble_chain[k] = 0.;
                        for ( j = 0 ; j < 2*halfnbw ; j ++) ensemble_chain[k] += chaine[itemp][k* (2*halfnbw) + j][i];
                        ensemble_chain[k] /= 2*halfnbw ;
                    }

                    // Compute the mean on the last moy_freq iterations
                     moylocal[i] = 0.;
//                     for (j = it - moy_freq ; j < it ; j ++)
//                     {
//                         moylocal[i] += chaine[itemp][j][i];
//                     }
//                     moylocal[i] /= moy_freq;

            // Calculate and save the autocorrelation function of the ensemble chain
                    Correlate(ensemble_chain, ensemble_chain, autocorrelation, ensemble_chain_index+1, (ensemble_chain_index+1)/2);
                    sprintf(fileouttemp, "autocorrelation_ensemblechain%i--temp=%.2f.res", i, 1./betas[itemp]);
                    Savetxt(fileouttemp, autocorrelation, (ensemble_chain_index+1)/2);


                    memcpy(ensemble_chain2, &ensemble_chain[ensemble_chain_index - ensemble_chain_index/2 ],  (ensemble_chain_index/2 +1)*sizeof(double) ); // Necessary to make a copy because acor changes it
                    acor( &moyacorhalf[i], &sigmahalf[i], &acortimehalf[i], ensemble_chain2, (ensemble_chain_index/2 +1));
                    acor( &moyacor[i], &sigma[i], &acortime[i], ensemble_chain, (ensemble_chain_index+1)  ); // Run the GW10 autocorrelation algorithm


            // Run the GW10 autocorrelation algorithm
                    fonction.Get_parameter_name(i, paramname);
                    sprintf(strinter,"Iteration %i ,  parameter %i '%s' : Mean %.5e (last %i iterations %.5e) | Std(Mean) %.5e | AcorTime %f  \n", it, i, paramname, moyacor[i], moy_freq, moylocal[i], sigma[i], acortime[i] );
                    strcat(strtoprint, strinter);
                    sprintf(strinter,"                              Second Half of the chain : Mean %.5e | Std(Mean) %.5e | AcorTime %f  \n\n",  moyacorhalf[i], sigmahalf[i], acortimehalf[i] );
                    strcat(strtoprint, strinter);
                    sprintf(strinter,"Autocorrelation half=%.2e last=%.2e %i \n", autocorrelation[(ensemble_chain_index+1)/2/2], autocorrelation[(ensemble_chain_index+1)/2 -1], (ensemble_chain_index+1)/2);
                    strcat(strtoprint, strinter);
                }

                if (calc_mean_chi2 == true)
                {
                  lnchi2y = fonction(moyacorhalf, betas[itemp]);
                  strcat(strtoprint, "\n");
                  sprintf(strinter,"Lnposterior (chi2) value at mean of the seconf half of the chain : %.10e (%.10e)\n",  lnchi2y, fonction.Get_chi2(lnchi2y, moyacorhalf, betas[itemp]) );
                  strcat(strtoprint, strinter);
                }
                maxlnchi2index = Maxofarray(lnchi2, 2*halfnbw);
                maxlnchi2 = lnchi2[maxlnchi2index];
                strcat(strtoprint, "\n");
                sprintf(strinter,"Maximum lnposterior in sample = %.10e\n",  maxlnchi2);
                strcat(strtoprint, strinter);


                acceptance_fraction_local = ( total_accepted - acceptance_fraction * ((it-moy_freq) * 2 * halfnbw) ) / (moy_freq*2 * halfnbw) ; // acceptance fraction on the last moy_freq iterations
                acceptance_fraction = total_accepted/ (it * 2 * halfnbw);
                strcat(strtoprint, "\n");
                sprintf(strinter,"Acceptance fraction on the last %i iterations is : %.5e\n", moy_freq, acceptance_fraction_local);
                strcat(strtoprint, strinter);
                sprintf(strinter,"Total acceptance fraction is : %.5e\n", acceptance_fraction);
                strcat(strtoprint, strinter);

// Experimental test !!
                // strcat(strtoprint, "\n");
                // sprintf(strinter, "Smoothed jumps: %i/%i\n", total_countersmoothing, (it * 2 * halfnbw));
                // strcat(strtoprint, strinter);

                sprintf(fileouttemp, "%s--temp=%.2f.res", fileout, 1./betas[itemp]);
                sprintf(strinter,"\n Saving data to %s \n", fileouttemp);
                strcat(strtoprint, strinter);

                Save_chain(fileouttemp, chaine[itemp], ndim, 2*halfnbw, ensemble_chain_index * 2 * halfnbw, chain_freq);

//                 if (record_refused == true)
//                 {
//                     sprintf(fileouttemp, "%s--temp=%.2f-refused.res", fileout, 1./betas[itemp]);
//                     Save_chain(fileouttemp, chaine_refused[itemp], ndim, 2*halfnbw, irefused, chain_freq);
//                 }
//
                strcat(strtoprint, "\n\n");

                }

            }

    // Copy back the modified temperature into the all-temperature array
            memcpy(Temp_walkers1d[itemp], walkers1d,  2* sizeof(double) * halfnbw * ndim) ;
            memcpy(Temp_lnchi2[itemp], lnchi2, 2* halfnbw * sizeof(double));
            Temp_accepted[itemp] = accepted;
            Temp_acceptance_fraction[itemp] = acceptance_fraction;
            if (record_refused == true) Temp_irefused[itemp]  = irefused;

        } // ############################ End of temperature loop #############################################################################################

        if (printstat ==0.)
        {

            if (mpirank == masterproc)
            {
                printf("\n###########################################################################");
                printf(" %i over %i and %i \n\n", ensemble_chain_index, ensemble_chain_size, it  );
                printf(" Time running : %.0f seconds\n\n", ((float)(clock()  - tclock_start)/CLOCKS_PER_SEC));
            }
#ifdef MPIMODE
            MPI_Collect_and_print(strtoprint, size_strtoprint, masterproc, printers, ncluster, MPI_COMM_WORLD);
#else
            printf("%s", strtoprint);
#endif

            // Print temperature swap statistics
            if (mpirank == masterproc)
            {
                printf("* Swap statistics:\n");
                inter = 0.;
                for (i = 0 ; i < ntemp -1 ; i++) inter += last_accepted_tempswaps[i];
                printf("Total Temperature swap average acceptance probability on the last %i iterations : %.5e \n", moy_freq, inter/(moy_freq*(ntemp-1)*2 * halfnbw/tempswapfreq));
                inter = 0.;
                for (i = 0 ; i < ntemp -1 ; i++) inter += accepted_tempswaps[i];
                printf("Total Temperature swap average acceptance probability : %.5e \n", inter/(it*(ntemp-1)*2 * halfnbw/tempswapfreq));

                for (i = 0 ; i < ntemp -1 ; i++)
                {
                    printf("Temperature swap average acceptance probability on the last %i iterations  temperature %i : %.5e \n", moy_freq, i, last_accepted_tempswaps[i]/(moy_freq*2 * halfnbw/tempswapfreq));
                    last_accepted_tempswaps[i] = 0.;
                    printf("Temperature swap average acceptance probability temperature %i: %.5e \n",  i, accepted_tempswaps[i]/(it*2 * halfnbw/tempswapfreq));

                }
            }
        } // Endif of printing

    } // ******** end of main loop

#ifdef MPIMODE
    code = MPI_Reduce(&accepted , &acceptance_fraction, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_tempcluster); // Compute total number of accepted moves
#else
    acceptance_fraction = accepted;
#endif

    acceptance_fraction /= (maxit * 2 * halfnbw);

    int acorfail = 0;
    int maxacortimei = 0;
    int burnin = 0;
    double maxacortime = 0.;


// Gather the refused chain
    if (record_refused == true)
    {
        for (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++)//(itemp = 0 ; itemp < ntemp ; itemp++)
        {
                code = MPI_Allgather( &Temp_irefused[itemp] , 1, MPI_INTEGER, irefuseds, 1, MPI_INTEGER, MPI_COMM_tempcluster);

                irefused=0;
                for (i = 0 ;i < nproc_per_temp; i++) irefused += irefuseds[i];
                if (irefused >=chain_size)
                {
                    if (mpirankincluster==0) printf("\nError : Too many refused jumps for temperature %i !\n\n", itemp);
                } else
                {
                    k=0;
                    for (i=0 ; i < chain_size; i++)
                    {
                    for(j = 1; j < nproc_per_temp ; j++)
                    {
                        if (i < irefuseds[j])
                        {
                                if (mpirankincluster == j)
                                {
                                    code = MPI_Send(chaine_refused[itemp][i] , (ndim+1), MPI_DOUBLE, 0, j,   MPI_COMM_tempcluster);
                                }
                                if (mpirankincluster ==0)
                                {
                                    code = MPI_Recv(chaine_refused[itemp][irefuseds[0]+k] , (ndim+1), MPI_DOUBLE, j,  j, MPI_COMM_tempcluster, MPI_STATUS_IGNORE);
                                    k+=1;
                                }
                        }

                        }
                    }
                    Temp_irefused[itemp] = k + irefuseds[0];
                    }
                    /*
                k=0;
                for (i=1; i< nbproc; i++)
                {
                    k += irefuseds[i-1];
                    for (j=0; j < irefuseds[i] ; j++) memcpy(chaine_refused[itemp][k+j], chaine_refused[itemp][i*(ndim+1)*chain_size/4+j], (ndim+1)*sizeof(double));
                }
                Temp_irefused[itemp] = k + irefuseds[nbproc-1];*/
        }
    }

    // Save the result "chaine" in "fileout"

    if (mpirankincluster==0)//(mpirank == masterproc)
    {
        strcpy(strtoprint,"");
        for (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++)//(itemp = 0 ; itemp < ntemp ; itemp++)
        {


            sprintf(strinter,"\n Acor analysis for final chain, Temperature= %.5e : \n", 1./betas[itemp]);
            strcat(strtoprint, strinter);
            acorfail = 0;
            maxacortimei = -1;

                for (i = 0 ; i < ndim ; i++)
                {

                        for (k = 0 ; k < ensemble_chain_size ; k++)
                        {   // Create a ensemble chain with each parameter averaged over the set of walkers
                            ensemble_chain[k] = 0.;
                            for ( j = 0 ; j < 2*halfnbw ; j ++) ensemble_chain[k] += chaine[itemp][k * (2*halfnbw) + j][i];
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
                        fonction.Get_parameter_name(i, paramname);
                        sprintf(strinter,"Parameter %i '%s' : Mean %f | Std(Mean) %f | AcorTime %f  | r %f | burnin %i \n", i, paramname, moyacor[i], sigma[i], acortime[i], acortime[i]/ (ensemble_chain_size - burnin), burnin );
                        strcat(strtoprint, strinter);
                    }

        lnchi2y = fonction(moyacorhalf, betas[itemp]);
        strcat(strtoprint, "\n");
        sprintf(strinter,"Lnposterior (chi2) value at mean : %.10e(%.10e)\n", lnchi2y, fonction.Get_chi2(lnchi2y, moyacorhalf, betas[itemp]) );
        strcat(strtoprint, strinter);

                if (maxacortimei >= 0)
                {
                    if (20 * acortime[maxacortimei] <  ensemble_chain_size)
                    {
                        if (acorfail == 1)
                        {
                            sprintf(strinter,"\n acor failed sometimes, retrying after discarding 20 times the longest acor time which is parameter %i, %f... \n\n", maxacortimei, maxacortime);
                            strcat(strtoprint, strinter);
                        }
                        else
                        {
                            sprintf(strinter,"\n Attempt to estimate burnin time : retrying after discarding 20 times the longest acor time which is parameter %i, %f... \n\n", maxacortimei, maxacortime);
                            strcat(strtoprint, strinter);
                        }
                        burnin = int( maxacortime )*20;
                        for (i = 0 ; i < ndim ; i++)
                        {
                            for (k = 0 ; k < ensemble_chain_size ; k++)
                                {   // Create a ensemble chain with each parameter averaged over the set of walkers
                                    ensemble_chain[k] = 0.;
                                    for ( j = 0 ; j < 2*halfnbw ; j ++) ensemble_chain[k] += chaine[itemp][k * (2*halfnbw) + j][i];
                                    ensemble_chain[k] /= 2*halfnbw ;
                                }
                            acor( &moyacor[i], &sigma[i], &acortime[i], &ensemble_chain[burnin], ensemble_chain_size - burnin );
                            sprintf(strinter,"Parameter %i : Mean %f | Std(Mean) %f | AcorTime %f  | r %f | burnin %i \n", i, moyacor[i], sigma[i], acortime[i], acortime[i]/ (ensemble_chain_size - burnin), burnin );
                            strcat(strtoprint, strinter);
                        }
                    } else
                        sprintf(strinter,"\nacor failed sometimes but there is not enough data even to remove the longest burning time from valid acor times\n");
                        strcat(strtoprint, strinter);
                }



            sprintf(strinter,"Acceptance fraction is : %.5e\n", acceptance_fraction);
            strcat(strtoprint, strinter);

            strcat(strtoprint, "\n\n");
            sprintf(fileouttemp, "%s--temp=%.2f.res", fileout, 1./betas[itemp]);
            sprintf(strinter,"Saving data to %s \n", fileouttemp);
            strcat(strtoprint, strinter);
            Save_chain(fileouttemp, chaine[itemp], ndim, 2*halfnbw, chain_size, chain_freq);

            if (record_refused == true)
            {
                sprintf(fileouttemp, "%s--temp=%.2f-refused.res", fileout, 1./betas[itemp]);
                Save_chain(fileouttemp, chaine_refused[itemp], ndim, 2*halfnbw, Temp_irefused[itemp], chain_freq);
            }

        } // end of temperature loop

    }

// Now actually print
    if (mpirank == masterproc)
    {
        printf("\n###########################################################################");
        printf("\n            Final Report \n");
        printf("\n###########################################################################");
    }
#ifdef MPIMODE
    MPI_Collect_and_print(strtoprint, size_strtoprint, masterproc, printers, ncluster, MPI_COMM_WORLD);

#else
    printf("%s",strtoprint);
#endif

    if (mpirank == masterproc)
    {
         printf( "\n ");
         printf("* Swap statistics:\n");
                inter = 0.;
                for (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++) inter += last_accepted_tempswaps[itemp];
                printf("Total Temperature swap average acceptance probability on the last %i iterations : %.5e \n", moy_freq, inter/(moy_freq*(ntemp-1)*2 * halfnbw/tempswapfreq));
                inter = 0.;
                for (itemp = icluster * ntemp_per_cluster; itemp <  min((icluster+1) * ntemp_per_cluster, ntemp); itemp++) inter += accepted_tempswaps[itemp];
                printf("Total Temperature swap average acceptance probability : %.5e \n", inter/(maxit*(ntemp-1)*2 * halfnbw/tempswapfreq));

                for (itemp = 0; itemp <   ntemp-1; itemp++)
                {
                    printf("Temperature swap average acceptance probability on the last %i iterations  temperature %i : %.5e \n", moy_freq, itemp, last_accepted_tempswaps[itemp]/(moy_freq*2 * halfnbw/tempswapfreq));
                    last_accepted_tempswaps[itemp] = 0.;
                    printf("Temperature swap average acceptance probability temperature %i: %.5e \n",  itemp, accepted_tempswaps[itemp]/(maxit*2 * halfnbw/tempswapfreq));

                }
    }

    if (record_refused == true)
    {
        if (mpirank == masterproc)  printf("\n* Number of refused walkers for each temperature: \n");

        if (mpirankincluster == 0) sPrint_table(strtoprint,Temp_irefused + ntemp_per_cluster*icluster, min((icluster+1) * ntemp_per_cluster, ntemp) - ntemp_per_cluster*icluster);
#ifdef MPIMODE
        MPI_Collect_and_print(strtoprint, size_strtoprint, masterproc, printers, ncluster, MPI_COMM_WORLD);
#else
        printf("%s", strtoprint);
#endif
    }
// End of printings

#ifdef MPIMODE
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    free(walkers);
    free(walkers1d);
    free(walkers1d_local);
    if (record_refused == true)
    {
        free(chaine_refused);
        free(irefuseds);
        free(Temp_irefused);
    }
    if (mpirankincluster == 0)
    {
        for (itemp = 0 ; itemp < ntemp ; itemp++)
        {
            for (k = 0 ; k < chain_size ; k++) free(chaine[itemp][k]);
            free(chaine[itemp]);
        }
        free(chaine);

        for (k = 0 ; k < moy_size ; k++)
        {
            free(moys[k]);
        }
        free(moys);
        free( ensemble_chain);
        free( ensemble_chain2);
        free(autocorrelation);
    }

//     free(tempfactor);
    free(betas);
    free(Temp_accepted);
    free(Temp_acceptance_fraction);
    free(last_accepted_tempswaps);
    free(accepted_tempswaps);
    free(walkerbuffer);
    for (itemp=0; itemp < ntemp ; itemp++)
    {
        free(Temp_lnchi2[itemp]);
        free(Temp_walkers1d[itemp]);
    }

    free(strtoprint);

    if (dodiagnostics)
    {
      for (k = 0 ; k < ndiagnostics ; k++)
      {
        free(diagnosticslog[k]);
      }
      free(diagnosticslog);
    }

    if (dosieving)
    {
      for (k = 0 ; k < maxsieved ; k ++) free(sievinglog[k]);
      free(sievinglog);
    }

#ifdef MPIMODE
    free(printers);
    code = MPI_Comm_free(&MPI_COMM_tempcluster);
    code = MPI_Finalize();
#endif

  return 0;
    }
