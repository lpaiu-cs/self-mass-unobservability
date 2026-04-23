// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

// Written by Guillaume Voisin, 2021
// email : guillaume.voisin@obspm.fr ;  astro.guillaume.voisin@gmail.com
//
// Converts a MCMC chain file to refer to a different parfile.
#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include <stdio.h>
#include <string.h>
#include "MCMC_lnposterior_functions.h"
#include "Utilities_MCMC_affinv.h"
#include "Constants.h"

int main ( int argc, char *argv[] )
{
    if (argc != 6)
    {
      printf("Syntax is : 'Convert_MCMC_chain.exe current_parfile timfile MCMC_result_file new_parfile new_MCMC_result_file'\n");
      printf("Note: a timfile is needed although it does not matter which one\n");
      return 1;
    }
    fctFittriple_gaussprior fonctionnew = fctFittriple_gaussprior();
    fonctionnew.Reconstruct(argv[4], argv[2]); // Standard version should be fonction.Load(char *, char *) 
    fctFittriple_gaussprior fonction = fctFittriple_gaussprior();
    fonction.Reconstruct(argv[1], argv[2]); // Standard version should be fonction.Load(char *, char *) 
   
   
    printf("\n Loaded parfile '%s' and timfile '%s'\n\n", argv[1], argv[2]);
    int ndim = fonction.Get_ndim();
    double params[ndim];

    //------------------ Load chain --- (This should be Load_chain drom Utilities_MCMC_affinv.h but there is a segfault problem)
    double ** chain ;
    int i , j;
    FILE * myfile ;
    myfile = fopen(argv[3], "r") ;
    int fndim=0;
    int fchain_freq =0;
    int fnb_walkers =0;
    int fchain_size =0;
    if (myfile == NULL )
    {
      printf("Could not read MCMC result file.");
      return 1;
    }
    fscanf(myfile, "#  %i  %i  %i %i\n",&fndim, &fnb_walkers, &fchain_size, &fchain_freq);
    if (ndim != fndim)
    {
        printf("\n Error ! Chain in file %s and desired chain have different number of dimensions : %i against %i \n\n", argv[3], fndim, ndim);
        return 1 ;
    }
    chain = (double**) malloc(( fchain_size) * sizeof(double*) );
    for (i =0  ; i < (fchain_size) ; i++)
    {
        chain[i] = (double*) malloc((ndim+1)*sizeof(double));
        for (j= 0 ; j < ndim + 1 ; j++) chain[i][j] = 0.;
    }

    for(i=0 ; i < fchain_size ; ++i)
    {
        if ( fscanf(myfile, "%le    ", &(chain[i][0]) ) == 1)
        {
          for (j = 1; j < fndim ; ++j ) {
               fscanf(myfile, "%le    ", &(chain[i][j]) );
           }
           fscanf(myfile, "%le    \n", &(chain[i][fndim]) );
       } else // if read fails line is skipped
       {
         i -= 1;
         fnextline(myfile);
       }
    }
    if (i < fchain_size)
    {
      printf("\nWarning : only %i/%i lines could be read from file %s\n", i, fchain_size, argv[3]);
      fchain_size = i;
    }
    // ---- end of loading chain

    // ----- Now creating new chain 
    int newndim = fonctionnew.Get_ndim();
    double ** newchain;
    
    newchain = (double**) malloc(( fchain_size) * sizeof(double*) );
    for (i =0  ; i < (fchain_size) ; i++)
    {
        newchain[i] = (double*) malloc((newndim+1)*sizeof(double));
        for (j= 0 ; j < newndim + 1 ; j++) newchain[i][j] = 0.;
    }
    
    vector<int> fpmap_current = fonction.parameters.Get_fitted_parameter_map();
    vector<int> fpmap = fonctionnew.parameters.Get_fitted_parameter_map();
    printf("\n Fitted parameter maps : \n");
    printf(" Current:");
    for (j=0; j < ndim; j++) printf("  %d", fpmap_current[j]);
    printf("\n     New:");
    for (j=0; j < newndim; j++) printf("  %d", fpmap[j]);
    printf("\n");
    
    int abspar, index, fittedpar, p, pcur, jcur;
    value_type shift, scale;
    // Random generator init
    // obtain a seed from the system clock:
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count() ;
    mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine

    for (j = 0; j < newndim ; j++)
    {
        printf("\n* New fitted parameter index %d (Fitted map index %d) : '%s' \n", j, fpmap[j], fonctionnew.parameters.fitparams_names[fpmap[j]]);
        fonctionnew.parameters.ParamList2AbsParam(fpmap[j], abspar, index); // Get absolute parameter number and index (for array parameters) of fonctionnew
        jcur = fonction.parameters.Get_fitted_parameter_index(abspar); // Get position in the fitted parameter list of fonction (of the first element if array-like parameter)
        
        p = fpmap[j];
        pcur = fpmap_current[jcur] + index;
        
        printf("  Current fitted parameter index %d (Fitted map index %d) \n", jcur + index, pcur);
        
        if (jcur > -1 and fonction.parameters.absolute_parameter_map[abspar+1] - fonction.parameters.absolute_parameter_map[abspar] > index) // Check parameter exists in "fonction"
        {
            scale = fonction.parameters.Parameter_scale(pcur)/fonctionnew.parameters.Parameter_scale(p);
            shift = (fonction.parameters.parameters_ini[pcur] - fonctionnew.parameters.parameters_ini[p])/fonctionnew.parameters.Parameter_scale(p);
            printf("  Init values : Current =%.5Le  New = %.5Le\n",  fonction.parameters.parameters_ini[pcur], fonctionnew.parameters.parameters_ini[p]);
            printf("                Scale = %.5Le  Shift = %.5Le\n", scale, shift);
            for (i = 0; i < fchain_size ; i++) newchain[i][j] = static_cast<double>(chain[i][jcur]*scale + shift);
        }
        else // If not attribute a random value to it
        {
            pcur = fonction.parameters.absolute_parameter_map[abspar] + index; 
            shift = (fonction.parameters.parameters_ini[pcur] - fonctionnew.parameters.parameters_ini[p])/fonctionnew.parameters.Parameter_scale(p);
            printf("  Init values : Current =%.5Le  New = %.5Le\n",  fonction.parameters.parameters_ini[pcur], fonctionnew.parameters.parameters_ini[p]);
            printf("                Scale = NA  Shift = %.5Le\n", shift);
            for (i = 0; i < fchain_size ; i++) newchain[i][j] = fonctionnew.initial_distribution(randomgen, j) + shift;
        }
    };
        
    // --------- Saving new chain
    
    Save_chain(argv[5], newchain, newndim, fnb_walkers, fchain_size, fchain_freq);
    
    free(chain);
    free(newchain);
     return 0;
}
