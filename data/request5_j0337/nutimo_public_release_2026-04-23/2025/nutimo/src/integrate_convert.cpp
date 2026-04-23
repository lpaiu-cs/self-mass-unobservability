// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code converts position and velocity vectors such as output by integrate.exe to other quantities
 * in particular orbital elements.
 * Written by Guillaume Voisin 2021 , LUTh, Observatoire de Paris, PSL Research University, CNRS (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifdef MPIMODE
    #include "mpi.h"
#endif
#include "Constants.h"
#include <cmath>
#include "Utilities.h"
#include <valarray>
#include <cstdlib>
#include "Parameters.h"
#include "AllTheories3Bodies.h"
#include "Diagnostics.h"
#include <chrono>
#include <string.h>
#include "Orbital_elements.h"
#include "Utilities_integrate.h"


int main ( int argc, char *argv[] ) {

    char parfile[500];
    char filein[500] ;
    char fileout[500] ;
    char fileout_reverse[500] ;
    char fileout_intofmotion[500] ;
    char file_masses[500] ;
    bool flag_imotion =false;
    int nflag_meanlong = 10 ;            // max number of meanlong flags, should be >= number of bodies
    bool flag_meanlong[nflag_meanlong];
    bool flag_masses = false;
    bool flag_reverse = false;
    int nt = -1;
    int i,j,k;
    
    for (i=0; i<nflag_meanlong ; i++) flag_meanlong[i] = false; 
    
    
// Read command line 
    
    if (argc >= 2)
    {
        strcpy(parfile,argv[1]) ;
        strcpy(filein,argv[2]) ;
        strcpy(fileout, filein);
        strcat(fileout, "-orbels");
        strcpy(fileout_reverse, filein);
        strcat(fileout_reverse, "-reverse");
        
        for (i = 3 ; i < argc ; i++)
        {
            if (strcmp(argv[i], "-o") == 0) 
            {
                if (argc > i+1)
                    {
                        strcpy(fileout,  argv[i+1]) ;
                        strcpy(fileout_reverse,  argv[i+1]) ;
                        strcat(fileout_reverse, "-reverse");
//                         if (mpirank == masterproc) printf("Reading turnfile : %s \n", turnfile);
                    }
                    else
                    {
                        printf("Failed to read output suffix !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
            }
            if (strcmp(argv[i], "-m") == 0) 
            {
                if (argc > i+1)
                    {
                        flag_masses = true;
                        strcpy(file_masses,  argv[i+1]) ;
//                         if (mpirank == masterproc) printf("Reading turnfile : %s \n", turnfile);
                    }
                    else
                    {
                        printf("Failed to read mass-file name after -m!\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
            }
            if (strcmp(argv[i], "-i") == 0) 
            {
                flag_imotion = true;
            }
            if (strcmp(argv[i], "-r") == 0) 
            {
                flag_reverse = true;
            }
            if (strcmp(argv[i], "-l") == 0) 
            {
                if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%i",&j) ;
                        if (j >= nflag_meanlong) 
                        {
                            printf("\n Error : integer following '-l' larger than maximum allowed (%d)\n\n", nflag_meanlong);
                        }
                        else
                        {
                            flag_meanlong[j] = true; 
                        }
                    }
                    else
                    {
                        printf("Failed to read index of mean longitude after argument -l !\n");
#ifdef MPIMODE
                        code = MPI_Finalize();
#endif
                        return 1;
                    }
            }
        }
    }
    else 
    {
        printf("Syntax: integrate_convert.exe parfile trajfile [-o outputprefix][-i][-l integer]\n");
        printf("Convert state vectors in trajfile (same format as output by integrate.exe) into other quantities, in particular orbital elements. By default, orbital elements are calculated w.r.t the invariant plane of the system.\n");
        printf("    parfile : nutimo type parfile \n");
        printf("    trajfile : file of the form 't(days) r_1 v_1 ... r_n v_n' where 'n' is the number of bodies and 'r' and 'v' are position and velocities. These coordinates will be converted to orbital elements and other derived quantities according to options. If -m is given then masses are taken from a separate file, otherwise a single set of masses is used thoughout and taken from the parfile. \n");
        printf("    -o outputprefix : prefix of output file(s). default=trajectories\n");
        printf("    -i : compute and save integrals of motion in a separate file\n");
        printf("    -l integer : compute mean longitude of the corresponding subsystem in hierarchical order and save resuld to a separate file with two columns (time in days | mean long in radians). Mean longitude = argument of asc. node (Oman) +  arg. of periastron (omperi) + mean anomaly (m). The option can be called several times creating one separate file each time. (0<= integer < nbody -1)\n");
        printf("    -m massfilename : each line of this file contain the set of masses (in solar mass) to use with the corresponding line of trajfile\n");
        printf("    -r : reverse action. Trajfile contains orbital elements that are converted back to position and velocities.\n");
        return 0;
    }
    
    strcpy(fileout_intofmotion, fileout);
    strcat(fileout_intofmotion, "-intofmotion");
    
    
    
    Parametres parameters(parfile);
    
// Load initial conditions
// Define these variables only to give output arguments to parameters.Compute_state_vectors which is needed to compute masses 
    valarray<value_type> rp(3);
    valarray<value_type> rpt(3);
    valarray<value_type> ri(3);
    valarray<value_type> rit(3);
    valarray<value_type> ro(3);
    valarray<value_type> rot(3);
    valarray<valarray<value_type>> r_extra(parameters.nextra);
    valarray<valarray<value_type>> v_extra(parameters.nextra);
    int nbody = 3 + parameters.nextra;

// Conserved quantities arrays 
    valarray<valarray<value_type>> center_of_mass_positions;
    valarray<valarray<value_type>> center_of_mass_impulsions;
    valarray<value_type>energies;
    valarray<valarray<value_type>> angular_momentum;
    
// Initialise parameter object    
    for (i = 0; i < parameters.nextra; i++)
    {
      r_extra[i] = valarray<value_type>(3);
      v_extra[i] = valarray<value_type>(3);
    }
    parameters.Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
// Create Mass array 
    valarray<value_type>Ms(nbody); // In units of Msun
    Ms[0] = parameters.Mp;
    Ms[1] = parameters.Mi;
    Ms[2] = parameters.Mo;
    for (i = 0; i < parameters.nextra; i++) Ms[3+i] = parameters.M_extra[i];
    
    int ncomments, nempty;
    Count_lines_in_file(filein, nt, ncomments, nempty);
    
    valarray<long double>  ts_day;
    valarray<valarray<valarray<long double>>>  svs;
    value_type *** state;
    value_type t;
    
    if (flag_reverse ==true)
    {
        orbel_t ** orbelsr;
        orbelsr = (orbel_t **) malloc(sizeof(orbel_t*)*nt);
        for (i=0; i < nt; i++) orbelsr[i] = (orbel_t *) malloc(sizeof(orbel_t)*(nbody-1));

        // Calculate dir_obs using parfile info (not necessarily accurate if orbels have been changed)
        state = (value_type***) malloc(sizeof(value_type**)*nbody);
        for (j=0; j < nbody; j++) 
        {
            state[j] = (value_type**) malloc(sizeof(value_type*));
            state[j][0] = (value_type*) malloc(sizeof(value_type)*6);
        }
        for (k=0; k < 3; k++) 
        {
            state[0][0][k] = rp[k];
            state[0][0][3+k] = rpt[k];
            state[1][0][k] = ri[k];
            state[1][0][3+k] = rit[k];
            state[2][0][k] = ro[k];
            state[2][0][3+k] = rot[k];
            for (j =3; j < nbody; j++)
            {
                state[j][0][k] = r_extra[j-3][k];
                state[j][0][3+k] = v_extra[j-3][k];
            }
        }
        Compute_integrals_of_motion(state, 1, Ms, parameters.Gg, parameters.gammabar, parameters.betabar, parameters.integrator_type, center_of_mass_positions, center_of_mass_impulsions,energies, angular_momentum);
        valarray<value_type> dir_obs = angular_momentum[0] / norm3d(angular_momentum[0]);
        for (j = 0 ; j < nbody; j++)
        {
            free(state[j][0]);
            free(state[j]);
        }
        free(state);
        // End of calculating dir_obs

        svs.resize(nt);
        Read_orbels(filein, nbody, ts_day, orbelsr, nt);
        for (i=0; i < nt; i++) 
        {
            orbel2statevects_bis_nbody(nbody, orbelsr[i], Ms, dir_obs, true, svs[i], t);
        }
        Save_traj(fileout_reverse, ts_day, svs);
        for (j = 0 ; j < nbody-1; j++) free(orbelsr[j]);
        free(orbelsr);
    }
    else
    {
        Read_trajectory(filein, nbody,  ts_day, svs);
    }
    
// Uncomment to print svs
//     for (i = 0 ; i < nt; i++)
//     {
//         for (j = 0 ; j< nbody ; j ++) 
//         {
//             for (k = 0 ; k< 6 ; k ++) printf("%.10Le  ", svs[i][j][k] );
//             printf("  |  ");
//         };
//         printf("\n");
//     }
    

   
//---- Compute orbital elements and optionnally mean longitudes

        printf("\nComputing orbital elements...\n");
        char fname[100];
        orbel_t ** orbels;
        orbels = (orbel_t **) malloc(sizeof(orbel_t*)*nt);
        for (i=0; i < nt; i++) orbels[i] = (orbel_t *) malloc(sizeof(orbel_t)*(nbody-1));
        valarray<valarray<value_type>> l(nbody-1);
        state = (value_type***) malloc(sizeof(value_type**)*nbody);
        for (j=0; j < nbody; j++) 
        {
            state[j] = (value_type**) malloc(sizeof(value_type*));
            state[j][0] = (value_type*) malloc(sizeof(value_type)*6);
            for (k=0; k < 6; k++) state[j][0][k] = svs[0][j][k];
        }
        
        Compute_integrals_of_motion(state, 1, Ms, parameters.Gg, parameters.gammabar, parameters.betabar, parameters.integrator_type, center_of_mass_positions, center_of_mass_impulsions,energies, angular_momentum);
        valarray<value_type> dir_obs = angular_momentum[0] / norm3d(angular_momentum[0]);
        valarray<value_type> masses(nbody);
        valarray<valarray<value_type>> Ms_file;
        if (flag_masses == true)
        {
            Loadtxt(file_masses, Ms_file, nt, nbody ) ;
        }
        for (j=0; j < nbody -1; j++) l[j] = valarray<value_type>(nt);
        for (i=0 ; i < nt ; i ++)
        {
            if (i%10000 ==0 ) printf("Converted %d/%d\n", i, nt);
            if (flag_masses == true) 
            {
                masses = Ms_file[i];
                // if masses change one must recompute dir_obs
                Compute_integrals_of_motion(state, 1, masses, parameters.Gg, parameters.gammabar, parameters.betabar, parameters.integrator_type, center_of_mass_positions, center_of_mass_impulsions,energies, angular_momentum);
                dir_obs = angular_momentum[0] / norm3d(angular_momentum[0]);
                
            }
            else
            {
                masses = Ms;
            }
            statevect2orbel_nbody(nbody, svs[i], masses, ts_day[i]*daysec, dir_obs, orbels[i]);
            for (j=0; j<nbody -1; j++) 
            {
                l[j][i] = orbels[i][j].Oman + orbels[i][j].omperi + orbels[i][j].m;
            }
        }
// Saving to file
        Save_orbels(fileout, ts_day, orbels, nbody);
        for (j=0; j < nbody -1; j++)
        {
            if (flag_meanlong[j] == true) 
            {
                sprintf(fname, "%s-l%d.dat", fileout, j);
                printf("Saving mean longitude %d to file: %s\n", j, fname);
                Save_meanlongs(fname, ts_day, l[j]);
            }
        }
        
        
    //---- Option : Compute integrals of motion
    if (flag_imotion == true)
    {
        printf("\nComputing integrals of motion: OPTION NOT IMPLEMENTED !\n");
//         printf("\nComputing integrals of motion...\n");
//         Compute_integrals_of_motion(states, nt, Ms, parameters.Gg, parameters.gammabar, parameters.betabar, parameters.integrator_type, center_of_mass_positions, center_of_mass_impulsions,energies, angular_momentum);
//         // Save to file 
//         printf("Saving integrals of motion to file: %s\n", fileout_intofmotion);
//         fout = fopen(fileout_intofmotion, "w");
//         fprintf(fout, "# t    E_Tot    r_cof   v_cof   L_Tot(Newt!)\n");
//         for (i=0 ; i < nt ; i ++)
//         {
//             fprintf(fout, "%.19Le    %.19Le    ", ts_day[i], energies[i]);
//             for (k = 0; k < 3 ; ++k ) fprintf(fout, "%.19Le    ", center_of_mass_positions[i][k] );
//             for (k = 0; k < 3 ; ++k ) fprintf(fout, "%.19Le    ", center_of_mass_impulsions[i][k] );
//             for (k = 0; k < 3 ; ++k ) fprintf(fout, "%.19Le    ", angular_momentum[i][k] );
//             fprintf(fout, "\n") ;
//         }
//         fclose(fout);
    }
    
    
// Cleaning 
    for (j = 0 ; j < nbody; j++)
    {
        free(state[j][0]);
        free(state[j]);
    }
    free(state);
    
    for (j = 0 ; j < nbody-1; j++) free(orbels[j]);
    free(orbels);
    

}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

