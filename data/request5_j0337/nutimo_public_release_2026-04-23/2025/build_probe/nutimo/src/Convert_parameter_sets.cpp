// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
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

#include "Parameters.h"
#include "Utilities_MCMC_affinv.h"
#include "IO.h"


using namespace std;

// TO COMPILE :
// g++ -std=c++11 Convert_parameter_sets.cpp Parameters.cpp Utilities_MCMC_affinv.cpp Utilities.cpp Orbital_elements.cpp -o Convert_parameter_sets.exe


int main ( int argc, char *argv[] ) {

  if ( argc < 3 )
  {
    printf("Syntax : Convert_parameter_sets.exe parfile target_parameter_set [MCMC_chain_to_convert] \n");
    return 1;
  }

  char parfile[200];
  strcpy(parfile, argv[1]);
  char parfile_converted[200];
  strcpy(parfile_converted, parfile);
  strcat(parfile_converted, "-converted");

  printf("%s to %s \n",parfile, parfile_converted);
  Parametres param(parfile);

  int pset_ini = param.parameter_set;
  int pset_target ;
  sscanf(argv[2], "%i", &pset_target);

  // Create a new paramter file with the new parameter set
  param.Convert_in_place(pset_target);
  param.Save_parfile(parfile_converted);

  double ** chain;
  // Convert the MCMC file if present
  if (argc >= 4 )
  {
    char MCMC_file[200];
    char MCMC_file_converted[200];
    strcpy(MCMC_file, argv[3]);
    strcpy(MCMC_file_converted, MCMC_file);
    strcat(MCMC_file_converted, "-converted");
    printf("%s to %s \n",MCMC_file, MCMC_file_converted);

    Parametres parammc(parfile);
    int ndim = parammc.fitted_parameters.size();
    int nb_walkers = 0;
    int chain_size = 0;
    int chain_freq = 0;
    double converted[ndim];

int i,j;
//###### Load chain ###########################################################################
FILE * myfile ;
myfile = fopen(MCMC_file, "r") ;
int fndim=0;
int fchain_freq =0;
int fnb_walkers =0;
fscanf(myfile, "#  %i  %i  %i %i\n",&fndim, &fnb_walkers, &chain_size, &fchain_freq);
if (ndim != fndim)
{
    printf("\n Error ! Chain in file %s and desired chain have different number of dimensions : %i against %i \n\n", MCMC_file, fndim, ndim);
    return 1;
}
chain = (double**) malloc((chain_size) * sizeof(double*) );
for (i =0  ; i < (chain_size) ; i++)
{
    chain[i] = (double*) malloc((ndim+1)*sizeof(double));
    for (j= 0 ; j < ndim + 1 ; j++) chain[i][j] = 0.;
}

for(i=0 ; i < chain_size ; ++i) {
    for (j = 0; j < fndim ; ++j )
    {
       fscanf(myfile, "%le    ", &(chain[i][j]) );
    }
    fscanf(myfile, "%le    \n", &(chain[i][fndim]) );
 }
 fclose(myfile);
//####################################################################################


    printf("chain size %i %.10e \n", chain_size, chain[0][0]);

    for (i = 0; i < chain_size; i++)
    {
        parammc.Convert_relative_shifts(pset_target, chain[i], converted);
      for (j = 0; j < ndim; j++) chain[i][j] = converted[j];
    }

    Save_chain(MCMC_file_converted, chain, ndim, nb_walkers, chain_size, chain_freq);

    for (i = 0; i < chain_size; i++) free(chain[i]);
    free(chain);

  }

}
