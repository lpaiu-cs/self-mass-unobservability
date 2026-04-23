// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

// Written by Guillaume Voisin, 2021
// email : guillaume.voisin@obspm.fr;  astro.guillaume.voisin@gmail.com
// Create a parfile from a MCMC chain according to different rules (mean, max..)
#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include <stdio.h>
#include <string.h>
#include "MCMC_lnposterior_functions.h"

int main ( int argc, char *argv[] )
{
    if (argc != 6)
    {
      printf("Syntax is : 'New_parfile_from_MCMC.exe parfile timfile MCMC_result_file new_parfile key'\n");
      printf("Create a new parfile from a MCMC chain. The 'key' decides what statistics to apply on the chain to derive a single set. Parameter scales are set equal to chain standard deviations. \n");
      printf("    parfile : parfile associated with the MCMC chain. \n");
      printf("    timfile : only needed for initialisation. Does not matter which one.\n");
      printf("    MCMC_result_file : chain as obtained using 'MCMC_pai_pt.exe'.\n");
      printf("    new_parfile : filename of the new parfile.\n");
      printf("    key : what parameters to extract from the MCMC chain result. Can be 'mean' or 'max' (max lnposterior)\n");
      return 1;
    }
    fctFittriple_gaussprior fonction = fctFittriple_gaussprior();
    fonction.Reconstruct(argv[1], argv[2]); // Standard version should be fonction.Load(char *, char *) 
    
    printf("\n Loaded parfile '%s' and timfile '%s' with chi2 = %.5e \n\n", argv[1], argv[2], fonction());
    int ndim = fonction.Get_ndim();
    double params[ndim];
    double means[ndim];
    double stds[ndim];
    
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

    
// Compute (relative) means
     for (j = 0 ; j < ndim ; j ++) means[j] = 0.;
       for (i = 0; i < fchain_size; i ++)
       {
         for (j = 0 ; j < ndim ; j ++)
         {
           means[j] += chain[i][j];
         }
       }
       for (j = 0 ; j < ndim ; j ++) means[j] /= fchain_size;
       
// Compute (absolute) standard deviations 
    for (j = 0 ; j < ndim ; j ++) stds[j] = 0.;
    for (i = 0; i < fchain_size; i ++)
    {
        for (j = 0 ; j < ndim ; j ++)
        {
            stds[j] += chain[i][j] * chain[i][j];
        }
    }
    for (j = 0 ; j < ndim ; j ++) stds[j] = sqrt(stds[j]/fchain_size - means[j]*means[j]) * fonction.parameters.Parameter_scale(fonction.parameters.fitted_parameters[j]);
    
// Choose parameter key
     if (strcmp(argv[5],"mean") == 0) // Calculate mean of the chain
     {
       for (j = 0 ; j < ndim ; j ++) params[j] = means[j];
     }
     else if (strcmp(argv[5],"max") == 0) // Return the parameter set with max lnposterior
     {
        int argmaxlnpost = 0;
        for (i = 0; i < fchain_size; i ++)
        {
            if (chain[argmaxlnpost][fndim] < chain[i][fndim]) argmaxlnpost = i;
        }
        for (j = 0 ; j < ndim ; j ++)
        {
           params[j] = chain[argmaxlnpost][j];
        }
     }
     else
     {
       printf("You need to choose what parameters to create a parfile from ! \n");
       return 1;
     }
    
// Apply params
     printf("Writing parameters '%s' to parfile at %s\n", argv[5], argv[4]);
     printf("Chi2 = %.3e\n", fonction(params));
    
// Set new scales from standard deviations
    for (j = 0 ; j < ndim ; j ++) fonction.parameters.Set_parameter_scale(fonction.parameters.fitted_parameters[j], static_cast<value_type>(stds[j]));
    
// Save to the new parfile
    fonction.Save_parfile(argv[4]);
     
     free(chain);
     return 0;
}
