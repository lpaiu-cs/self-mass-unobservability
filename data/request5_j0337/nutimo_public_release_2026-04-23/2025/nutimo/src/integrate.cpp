// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code integrates the equations of motion for the triple system J0337+1715 (Ransom et al. 2013). The inputs are Nutimo parameters
 *
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
    char fileout[500] ;
    char fileout_intofmotion[500] ;
    char fileout_resume[500] ;
    char resumefile[500] ;
    bool flag_imotion =false;
    bool flag_savetraj =true;
    int nflag_meanlong = 10 ;            // max number of meanlong flags, should be >= number of bodies
    bool flag_meanlong[nflag_meanlong];
    bool flag_anymeanlong = false; 
    bool flag_resume = false;
    value_type tstart ;
    value_type tend; 
    int npts = -1;
    int i,j,k;
    
    for (i=0; i<nflag_meanlong ; i++) flag_meanlong[i] = false; 
    
    strcpy(fileout, "trajectories");
    strcpy(fileout_intofmotion, fileout);
    strcat(fileout_intofmotion, "-intofmotion");
    strcpy(fileout_resume, fileout);
    strcat(fileout_resume, "-resume");
    
    
// Read command line 
    
    if (argc >= 2)
    {
        strcpy(parfile,argv[1]) ;
        sscanf(argv[2], "%Lf",&tstart) ;
        sscanf(argv[3], "%Lf",&tend) ;
        
        for (i = 2 ; i < argc ; i++)
        {
            if (strcmp(argv[i], "-o") == 0) 
            {
                if (argc > i+1)
                    {
                        strcpy(fileout,  argv[i+1]) ;
                        strcpy(fileout_intofmotion, fileout);
                        strcat(fileout_intofmotion, "-intofmotion");
                        strcpy(fileout_resume, fileout);
                        strcat(fileout_resume, "-resume");
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
            if (strcmp(argv[i], "-n") == 0) 
            {
                if (argc > i+1)
                    {
                        sscanf(argv[i+1], "%i",&npts) ;
//                         if (mpirank == masterproc) printf("Reading turnfile : %s \n", turnfile);
                    }
                    else
                    {
                        printf("Failed to read number of points argument -n !\n");
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
            if (strcmp(argv[i], "-no") == 0) 
            {
                flag_savetraj = false;
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
                            flag_anymeanlong = true;                       
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
            if (strcmp(argv[i], "-r") == 0) 
            {
                if (argc > i+1)
                    {
                        strcpy(resumefile,  argv[i+1]) ;
                        flag_resume =true;
                    }
                    else
                    {
                        printf("Failed to read resume file after argument -r !\n");
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
        printf("Syntax: integrate.exe parfile start end [-o outputprefix][-n integer][-i][-l integer][-r resumefile]\n");
        printf("Integrate equations of motion using parameters from 'parfile' and save results to a file.\n");
        printf("The output format is r0_x r0_y r0_z v0_x v0_y v0_z ... rn_x .. vn_z for n+1 bodies.");
        printf("    parfile : nutimo type parfile \n");
        printf("    start : start date in days <=0 relative to the reference time in parfile \n");
        printf("    end : end date in days >=0 relative to the reference time in parfile \n");
        printf("    -o outputprefix : prefix of output file(s). default=trajectories\n");
        printf("    -no : No Output of trajectories to file. Does prevent output from -i, -l or other options. \n");
        printf("    -n integer : Number of times to record on an evenly sampled grid. default=interpolation parameters from parfile\n");
        printf("    -i : compute and save integrals of motion in a separate file\n");
        printf("    -l integer : compute mean longitude of the corresponding subsystem in hierarchical order and save resuld to a separate file with two columns (time in days | mean long in radians). Mean longitude = argument of asc. node (Oman) +  arg. of periastron (omperi) + mean anomaly (m). The option can be called several times creating one separate file each time. (0<= integer < nbody -1)\n");
        printf("    -r resumefile : initial state vector is overwritten by the one found in 'resumefile', e.g. a file containing the last state vector of a previous trajetory file.\n");
        return 0;
    }
    
    
    
    Parametres parameters(parfile);
    Integrateur intgr;
    
// Load initial conditions
     
    intgr.x0.resize(18 + parameters.nextra*6);  // allocate initial state vector
    intgr.tolint = parameters.tolint ;
    intgr.integrator_type = parameters.integrator_type ;

    valarray<value_type> rp(3);
    valarray<value_type> rpt(3);
    valarray<value_type> ri(3);
    valarray<value_type> rit(3);
    valarray<value_type> ro(3);
    valarray<value_type> rot(3);
    valarray<valarray<value_type>> r_extra(parameters.nextra);
    valarray<valarray<value_type>> v_extra(parameters.nextra);
    value_type timedays;
    for (i = 0; i < parameters.nextra; i++)
    {
      r_extra[i] = valarray<value_type>(3);
      v_extra[i] = valarray<value_type>(3);
    }
    parameters.Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
    if (flag_resume == true) 
    {
        j =  Read_resume_file(resumefile, rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
        if (intgr.x0.size() != j)
        {
            printf("\nFAILED TO RESUME : only %d out of %d components read successfully.\n", j, intgr.x0.size());
            return 1;
        }
    }

// // Uncomment to print initial state vectors 
//      printf("INITIAL STATE VECTOR\n");
//      printf("%.10Le  %.10Le  %.10Le  %.10Le  %.10Le  %.10Le |  ", rp[0], rp[1], rp[2], rpt[0], rpt[1], rpt[2]);
//     printf("%.10Le  %.10Le  %.10Le  %.10Le  %.10Le  %.10Le |  ", ri[0], ri[1], ri[2], rit[0], rit[1], rit[2]);
//     printf("%.10Le  %.10Le  %.10Le  %.10Le  %.10Le  %.10Le |  ", ro[0], ro[1], ro[2], rot[0], rot[1], rot[2]);
//     for (i = 0; i < parameters.nextra; i++)
//     {
//         printf("%.10Le  %.10Le  %.10Le  %.10Le  %.10Le  %.10Le |  ", r_extra[i][0], r_extra[i][1], r_extra[i][2], v_extra[i][0], v_extra[i][1], v_extra[i][2]);
//     }
    
    
    // Set caracteristic time and length scales
    intgr.mass = Msol ;
    if (parameters.Mo > 0.) {
    intgr.length = ( norm3d(rp) + norm3d(ri) + norm3d(ro) ) / trois ;
    intgr.timescale = ( norm3d(rp) / norm3d(rpt) + norm3d(ri) / norm3d(rit) + norm3d(ro) / norm3d(rot) ) / trois ;
    }
    else {
    intgr.length = ( norm3d(rp) + norm3d(ri)) / deux ;
    intgr.timescale = ( norm3d(rp) / norm3d(rpt) + norm3d(ri) / norm3d(rit) ) / deux ;
    }
    timedays = intgr.timescale / daysec;

    // Set other properties inherited from "Integrateur"
    value_type dt_interp =  parameters.Pi / parameters.interpsteps_per_period_i;
    if (npts <= 0)
        intgr.n = static_cast<long int>( abs(tstart - tend) / dt_interp )  ; 
    else 
        intgr.n = npts;
    value_type dt = abs(tend - tstart)/intgr.n;
    intgr.t0 = parameters.treference / timedays ;
    intgr.dt0 = dt_interp / timedays / cinq;

    // ********* Pass on parameters of the equations of motion to the integrator**********************
    intgr.int_M0 = parameters.Mp ;
    intgr.int_M1 = parameters.Mi ;
    intgr.int_M2 = parameters.Mo ;
    intgr.int_SEP_D = parameters.SEP_D ;
    intgr.int_M_extra = parameters.M_extra;
    intgr.int_Gg = parameters.Gg;
    intgr.int_gammabar = parameters.gammabar;
    intgr.int_betabar = parameters.betabar;
    for (i =0 ; i < (min(intgr.int_quadrupole.size(), parameters.quadrupole.size())) ; i++) intgr.int_quadrupole[i] = parameters.quadrupole_kgm2[i];

    // Add parameter here if your parameter is used as a parameter in the equations of motions (like masses)
    // *****************************************************************************
        // Initialize the other time arrays
              
        valarray<value_type> ts_day(intgr.n); // integration steps in days. 
        intgr.ts.resize(intgr.n);
        if (tend == 0) 
            intgr.ntneg = intgr.n - 1;
        else if (tstart == 0)
            intgr.ntneg = 0;
        else 
            intgr.ntneg = intgr.n - static_cast<int>(floor(tend / dt)) -1 ;
        intgr.tneg.resize( intgr.ntneg);
        intgr.tpos.resize( intgr.n - intgr.ntneg ) ;
        for (i = 1 ; i <= intgr.ntneg ; ++i ) {
            ts_day[intgr.ntneg -i] = -dt * i ; 
            intgr.ts[intgr.ntneg -i] = ts_day[intgr.ntneg -i] / timedays ;
            intgr.tneg[intgr.ntneg-i] = abs( intgr.ts[intgr.ntneg-i] ) ;
        }
        for (i = intgr.ntneg ; i < intgr.n ; ++i ) {
            ts_day[i] = dt * (i-intgr.ntneg) ;
            intgr.ts[i] = ts_day[i] / timedays ;
            intgr.tpos[i-intgr.ntneg] = intgr.ts[i] ;
        }
        printf("\nReference time at index %ld in [1,%ld]\n\n", intgr.ntneg, intgr.n);


    // Set intial state vector for integration of the equation of motion
        
        intgr.x0[0] = rp[0] / intgr.length ;
        intgr.x0[1] = rp[1] / intgr.length ;
        intgr.x0[2] = rp[2] / intgr.length ;
        intgr.x0[3] = ri[0] / intgr.length ;
        intgr.x0[4] = ri[1] / intgr.length ;
        intgr.x0[5] = ri[2] / intgr.length ;
        intgr.x0[6] = ro[0] / intgr.length ;
        intgr.x0[7] = ro[1] / intgr.length ;
        intgr.x0[8] = ro[2] / intgr.length ;
        for (i = 0 ; i < parameters.nextra ; i ++) 
        {
            intgr.x0[9 + i*3] = r_extra[i][0] / intgr.length ;
            intgr.x0[10 + i*3] = r_extra[i][1] / intgr.length ;
            intgr.x0[11 + i*3] = r_extra[i][2] / intgr.length ;
        }
        
        int extrashift = parameters.nextra * 3;
        
        intgr.x0[9 + extrashift] = rpt[0] / intgr.length  * intgr.timescale ;
        intgr.x0[10 + extrashift] = rpt[1] / intgr.length * intgr.timescale ;
        intgr.x0[11 + extrashift] = rpt[2] / intgr.length * intgr.timescale ;
        intgr.x0[12 + extrashift] = rit[0] / intgr.length * intgr.timescale ;
        intgr.x0[13 + extrashift] = rit[1] / intgr.length * intgr.timescale ;
        intgr.x0[14 + extrashift] = rit[2] / intgr.length * intgr.timescale ;
        intgr.x0[15 + extrashift] = rot[0] / intgr.length * intgr.timescale ;
        intgr.x0[16 + extrashift] = rot[1] / intgr.length * intgr.timescale ;
        intgr.x0[17 + extrashift] = rot[2] / intgr.length * intgr.timescale ;
        for (i = 0 ; i < parameters.nextra ; i ++) 
        {
            intgr.x0[18 + extrashift + i*3] = v_extra[i][0] / intgr.length * intgr.timescale  ;
            intgr.x0[19 + extrashift + i*3] = v_extra[i][1] / intgr.length * intgr.timescale  ;
            intgr.x0[20 + extrashift + i*3] = v_extra[i][2] / intgr.length * intgr.timescale  ;
        }
        
        Print_table(intgr.x0);
        printf("t0 = %Le \n", intgr.t0);
        printf("Integrator type = %d \n", intgr.integrator_type);
        
        
// Integrate equations of motion
        
        int nbody = 3 + parameters.nextra;
        
        value_type *** states;
        states = (value_type***) malloc(sizeof(value_type**) * nbody); // states[body][time][x y z vx vy vz]
        for (j=0; j< nbody ; j ++) 
        {
            states[j] = (value_type**) malloc(sizeof(value_type*) * intgr.n);
            for (i=0; i< intgr.n ; i ++)
            {
                states[j][i] = (value_type*) malloc(sizeof(value_type) * 6);
            }
        }

        printf("\nRunning integration...\n");
        clock_t tclock_start=clock();
        
        intgr.Integrate_Allways() ;

        for (i = 0; i < intgr.n ; ++i) 
        {
            for (j = 0 ; j < nbody ; j ++) 
            {
                states[j][i][0] = intgr.states[i][j*3] * intgr.length ;
                states[j][i][1] = intgr.states[i][1 + j*3] * intgr.length ;
                states[j][i][2] = intgr.states[i][2 + j*3] * intgr.length ;
                states[j][i][3] = intgr.states[i][3*nbody + j*3] * intgr.length / intgr.timescale;
                states[j][i][4] = intgr.states[i][3*nbody + 1 + j*3] * intgr.length / intgr.timescale;
                states[j][i][5] = intgr.states[i][3*nbody + 2 + j*3] * intgr.length / intgr.timescale;
            }
        }

    printf("Integration took %f seconds\n", ((float)(clock()  - tclock_start)/CLOCKS_PER_SEC));
    
    
// Save trajectories to file   
    FILE *fout;
    if (flag_savetraj == true)
    {
        printf("\nSaving travectories to file: %s\n", fileout);
        fout = fopen(fileout, "w");
        for (i=0 ; i < intgr.n ; i ++)
        {
            fprintf(fout, "%.19Le            ", ts_day[i]);
            for (j = 0; j < nbody ; ++j ) 
            {
                for (k = 0; k < 6 ; ++k ) fprintf(fout, "%.19Le    ", states[j][i][k] );
                fprintf(fout, "            ");
            }
            fprintf(fout, "\n") ;
        }
        fclose(fout);
    }
    
// Create resume file (i.e. last line of trajectory file)
    fout = fopen(fileout_resume, "w");
    for (j = 0; j < nbody ; ++j ) 
    {
        for (k = 0; k < 6 ; ++k ) fprintf(fout, "%.19Le    ", states[j][intgr.n-1][k] );
        fprintf(fout, "            ");
    }
    fclose(fout);

// ------------------------ Options ------------------------
// Create Mass array 
    valarray<value_type>Ms(nbody); // In units of Msun
    Ms[0] = parameters.Mp;
    Ms[1] = parameters.Mi;
    Ms[2] = parameters.Mo;
    for (i = 0; i < parameters.nextra; i++) Ms[3+i] = parameters.M_extra[i];
// Conserved quantities arrays 
    valarray<valarray<value_type>> center_of_mass_positions;
    valarray<valarray<value_type>> center_of_mass_impulsions;
    valarray<value_type>energies;
    valarray<valarray<value_type>> angular_momentum;
        
//---- Option : Compute integrals of motion
    if (flag_imotion == true)
    {
        printf("\nComputing integrals of motion...\n");
        Compute_integrals_of_motion(states, intgr.n, Ms, intgr.int_Gg, intgr.int_gammabar, intgr.int_betabar, intgr.integrator_type, center_of_mass_positions, center_of_mass_impulsions,energies, angular_momentum);
        printf("Saving integrals of motion to file: %s\n", fileout_intofmotion);
        fout = fopen(fileout_intofmotion, "w");
        fprintf(fout, "# t    E_Tot    r_cof   v_cof   L_Tot(Newt!)\n");
        for (i=0 ; i < intgr.n ; i ++)
        {
            fprintf(fout, "%.19Le    %.19Le    ", ts_day[i], energies[i]);
            for (k = 0; k < 3 ; ++k ) fprintf(fout, "%.19Le    ", center_of_mass_positions[i][k] );
            for (k = 0; k < 3 ; ++k ) fprintf(fout, "%.19Le    ", center_of_mass_impulsions[i][k] );
            for (k = 0; k < 3 ; ++k ) fprintf(fout, "%.19Le    ", angular_momentum[i][k] );
            fprintf(fout, "\n") ;
        }
        fclose(fout);
    }
//---- Option : Compute mean longitude
    if (flag_anymeanlong == true)
    {
        printf("\nComputing mean longitudes...\n");
        char fname[100];
        orbel_t orbels[nbody-1];
        valarray<valarray<value_type>> l(nbody-1);
        valarray<valarray<value_type>> svs(nbody);
        value_type *** state;
        state = (value_type***) malloc(sizeof(value_type**)*nbody);
        for (j=0; j < nbody; j++) 
        {
            state[j] = (value_type**) malloc(sizeof(value_type*));
            state[j][0] = (value_type*) malloc(sizeof(value_type)*6);
            for (k=0; k < 6; k++) state[j][0][k] = states[j][intgr.ntneg][k];
        }
        
        Compute_integrals_of_motion(state, 1, Ms, intgr.int_Gg, intgr.int_gammabar, intgr.int_betabar, intgr.integrator_type, center_of_mass_positions, center_of_mass_impulsions,energies, angular_momentum);
        valarray<value_type> dir_obs = angular_momentum[0] / norm3d(angular_momentum[0]);
        for (j=0; j<nbody -1; j++) l[j] = valarray<value_type>(intgr.n);
        for (i=0 ; i < intgr.n ; i ++)
        {
            for (j=0; j < nbody; j++) svs[j] = valarray<value_type>(states[j][i],6);
            statevect2orbel_nbody(nbody, svs, Ms, ts_day[i]*daysec, dir_obs, orbels);
            for (j=0; j<nbody -1; j++) 
            {
                l[j][i] = orbels[j].Oman + orbels[j].omperi + orbels[j].m;
            }
        }
        for (j=0; j < nbody -1; j++)
        {
            if (flag_meanlong[j] == true) 
            {
                sprintf(fname, "%s-l%d.dat", fileout, j);
                printf("Saving mean longitude %d to file: %s\n", j, fname);
                Save_meanlongs(fname, ts_day, l[j]);
            }
        }
        
        for (j = 0 ; j < nbody; j++)
        {
            free(state[j][0]);
            free(state[j]);
        }
        free(state); 
    }
    
// Cleaning 

    for (j = 0 ; j < nbody; j++)
    {
        for (i = 0 ; i< intgr.n ; i ++) free(states[j][i]);
        free(states[j]);
    }
    free(states);
    
//
    printf("Total time: %f seconds\n", ((float)(clock()  - tclock_start)/CLOCKS_PER_SEC));

}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

