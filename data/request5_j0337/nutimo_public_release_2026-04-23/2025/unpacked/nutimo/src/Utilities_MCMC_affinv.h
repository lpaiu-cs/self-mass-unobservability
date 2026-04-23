/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Utilities_MCMC_affinv_h
# define Utilities_MCMC_affinv_h


#ifdef MPIMODE
    #include "mpi.h"
#endif
#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include "Utilities.h"
#include <stdio.h>
#include <string.h>

class gw10_distribution // the g distribution in gw10
{
    double a;
    double cte;
    double Ng;

public:

    gw10_distribution(double parameter_a){a = parameter_a; Ng = 2. * (sqrt(a) - sqrt(1./a)); cte = - 2. / (sqrt(a) * Ng );};
    double operator()(double uniform_random_number_01){return 0.25 * pow(Ng *(uniform_random_number_01 - cte), 2 );}; // uniform_random_number_01 is a random number between [0,1[ drawn from a uniform distribution
};

double test_gw10_distribution(int nb_tirages, double a); // test the distribution function defined by gw10_distribution by comparing the ratio of (number events in [a/2, a]) / (number events in [1/a, a/2[) with the theoretical value. there are nb_tirages events.

void Save_chain(char * filename, double ** chain, int ndim, int nb_walkers, int chain_size, int chain_freq);

// !!! THERE IS A SEGFAULT PROBLEM WITH THIS WHICH I CANNOT FIND OUT. IT WORKS IF THE CODE IS COPY PASTED DIRECTLY INSTEAD OF CALLING THE ROUTINE 
void Load_chain(char * filename, double ** chain, int ndim, int nb_walkers, int & chain_size, int chain_freq);

int Init_from_prev_chain(char * prev_chain_file, double * walkers1d, double * lnposteriors, const  int ndim, const int nb_walkers);


void Print_walker(double * walker, double lnposterior, int ndim, char * printedwalker);

#endif
