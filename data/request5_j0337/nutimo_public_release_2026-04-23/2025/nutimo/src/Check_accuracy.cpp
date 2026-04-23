// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

// This file is inspired by Quad4FMain.cxx coming with the Minuit distribution
// Written by Guillaume Voisin, LUTh, Observatoire de Paris, March 2016


// To run correctly first do :  export OMP_NUM_THREADS=1
#include "Fittriple.h"
#include "Utilities.h"
#include "Constants.h"
#include "IO.h"
// MODIF LG
#include "tempo2.h"


int main ( int argc, char *argv[] ) {


  char * parfile;
  char * datafile ;
  double tol = 0.01;
  double factor =1.5;
//   char * fileout ;

  long int i = 0;

  valarray<value_type> maxcof(3);
  valarray<value_type> mincof(3);
  valarray<value_type> maxmom(3);
  valarray<value_type> minmom(3);
  valarray<value_type> inicomposition(3);
  valarray<value_type> inicomvelocity(3);
  valarray<value_type> meancomvelocity(3);
  value_type Mtot;
  value_type minnrj=zero;
  value_type maxnrj=zero;



    if ( argc == 5 )
    {
        parfile = argv[1] ;
        datafile = argv[2];
        sscanf(argv[3], "%lf", &factor);
        sscanf(argv[4], "%lf", &tol);
//         fileout = argv[3];
    }
    else
    {
        printf("The syntax is: \n Check_accuracy parfile timfile factor tolerance\n");
        return 1;
    }
    Fittriple fit(parfile, datafile);

    printf("\n*** Doing the equation of motion conserved quantities check ***\n\n");
    fit.Compute_integrals_of_motion();

    for (i =0; i <3; i++)
    {
        maxcof[i] = fit.center_of_mass_positions[0][i];
        maxmom[i] =  fit.center_of_mass_impulsions[0][i];
        mincof[i] = fit.center_of_mass_positions[0][i];
        minmom[i] =  fit.center_of_mass_impulsions[0][i];
        inicomposition[i] = fit.center_of_mass_positions[fit.treference_in_interp][i];
        inicomvelocity[i] = fit.center_of_mass_impulsions[fit.treference_in_interp][i];
    }

    maxnrj = fit.energies[0];
    minnrj = fit.energies[0];

    Mtot = (fit.parameters.Mp + fit.parameters.Mi + fit.parameters.Mo) + fit.energies[fit.treference_in_interp] / (clight*clight);
    inicomvelocity /= Mtot;

    printf("\n Initial position of the centre of mass : %.5Le %.5Le %.5Le m\n", inicomposition[0], inicomposition[1], inicomposition[2]);
    printf(" Initial velocity of the centre of mass : %.5Le %.5Le %.5Le m/s\n", inicomvelocity[0], inicomvelocity[1], inicomvelocity[2]);

    meancomvelocity = 0;

    for (i = 0; i < fit.ninterp; i++)
    {
        if (maxcof[0] < fit.center_of_mass_positions[i][0]) maxcof[0] = fit.center_of_mass_positions[i][0];
        if (maxcof[1] < fit.center_of_mass_positions[i][1]) maxcof[1] = fit.center_of_mass_positions[i][1];
        if (maxcof[2] < fit.center_of_mass_positions[i][2]) maxcof[2] = fit.center_of_mass_positions[i][2];

        if (mincof[0] > fit.center_of_mass_positions[i][0]) mincof[0] = fit.center_of_mass_positions[i][0];
        if (mincof[1] > fit.center_of_mass_positions[i][1]) mincof[1] = fit.center_of_mass_positions[i][1];
        if (mincof[2] > fit.center_of_mass_positions[i][2]) mincof[2] = fit.center_of_mass_positions[i][2];

        if (maxmom[0] < fit.center_of_mass_impulsions[i][0]) maxmom[0] = fit.center_of_mass_impulsions[i][0];
        if (maxmom[1] < fit.center_of_mass_impulsions[i][1]) maxmom[1] = fit.center_of_mass_impulsions[i][1];
        if (maxmom[2] < fit.center_of_mass_impulsions[i][2]) maxmom[2] = fit.center_of_mass_impulsions[i][2];

        if (minmom[0] > fit.center_of_mass_impulsions[i][0]) minmom[0] = fit.center_of_mass_impulsions[i][0];
        if (minmom[1] > fit.center_of_mass_impulsions[i][1]) minmom[1] = fit.center_of_mass_impulsions[i][1];
        if (minmom[2] > fit.center_of_mass_impulsions[i][2]) minmom[2] = fit.center_of_mass_impulsions[i][2];

        if (maxnrj < fit.energies[i]) maxnrj = fit.energies[i];
        if (minnrj > fit.energies[i]) minnrj = fit.energies[i];

        meancomvelocity[0] += fit.center_of_mass_impulsions[i][0];
        meancomvelocity[1] += fit.center_of_mass_impulsions[i][1];
        meancomvelocity[2] += fit.center_of_mass_impulsions[i][2];

    }
    meancomvelocity /= fit.ninterp * Mtot;
    printf(" Mean velocity of the centre of mass : %.5Le %.5Le %.5Le m/s\n\n", meancomvelocity[0], meancomvelocity[1], meancomvelocity[2]);

    printf(" Initial value of the energy (Joules/MSol): %.3Le \n",fit.energies[fit.treference_in_interp]);
    printf(" Initial values of the center-of-mass momentum (Msol.m/s): %.3Le  %.3Le  %.3Le \n\n",fit.center_of_mass_impulsions[fit.treference_in_interp][0],fit.center_of_mass_impulsions[fit.treference_in_interp][1],fit.center_of_mass_impulsions[fit.treference_in_interp][2]);

    printf("The maximum energy variation is %.15Le Joules/Msol representing %.15Le percent of the initial value \n", maxnrj - minnrj, (maxnrj - minnrj)/fit.energies[fit.treference_in_interp]);
    for (i =0; i <3; i++)

    {
        printf("The maximum center-of-mass position variation component %Li is maximum %.3Le meters. \n",i, maxcof[i] - mincof[i]);
        printf("                                     With mean velocity correction: %.3Le meters. \n",i, fit.center_of_mass_positions[fit.ninterp-1][i] - fit.center_of_mass_positions[0][i] - meancomvelocity[i] * (fit.tinterp[fit.ninterp-1] - fit.tinterp[0]) * daysec );
        printf("The maximum center-of-mass velocity component %Li is %.3Le m/s \n", i, (maxmom[i] - minmom[i])/Mtot);
    }

    printf("\n*** Doing the interpolation grid check ***\n\n");
    fit.Check_tinterp_grid(tol, factor);
}
