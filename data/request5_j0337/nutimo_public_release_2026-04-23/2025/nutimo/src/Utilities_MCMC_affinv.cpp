// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
*/
#include "Utilities_MCMC_affinv.h"


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


// !!! THERE IS A SEGFAULT PROBLEM WITH THIS WHICH I CANNOT FIND OUT. IT WORKS ID THE CODE IS COPY PASTED DIRECTLY INSTEAD OF CALLING THE ROUTINE 
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
        for (j = 0; j < fndim ; ++j )
        {
           fscanf(myfile, "%le    ", &(chain[i][j]) );
        }
        fscanf(myfile, "%le    \n", &(chain[i][fndim]) );
        for (int j = 0; j < ndim+1; j++) printf("%.10e  ", chain[i][j]);
        printf("\n");
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

int Init_from_prev_chain(char * prev_chain_file, double * walkers1d, double * lnposteriors, const  int ndim, const int nb_walkers)
{

    int i , j;
    FILE * myfile ;
    myfile = fopen(prev_chain_file, "r") ;
    if (myfile == NULL) return 1;
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
         lnposteriors[i] = chain[fchain_size - 1  - i][ndim] ;
         // printf("reading previous lnpost %.15e %.5e\n", lnposteriors[i], chain[fchain_size - 1  - i][ndim]);
     }
    // printf("\nPrinting read chain !\n");
     //Print_table(&chain[fchain_size - nb_walkers], nb_walkers, ndim+1);
     for (i =0  ; i < (fchain_size) ; i++) free(chain[i]);
     free(chain);

     return 0 ;
}


void Print_walker(double * walker, double lnposterior, int ndim, char * printedwalker)
{
  char strinter[100];
  sprintf(printedwalker,"");
  sprintf(strinter,"");
  for (int i = 0; i < ndim ; i++)
  {
      sprintf(strinter, "%.4e ", walker[i]);
      strcat(printedwalker,strinter);
  }
  sprintf(strinter,"| %.5e \n", lnposterior);
  strcat(printedwalker,strinter);
}
