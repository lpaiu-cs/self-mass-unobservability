/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>

 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 * 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifndef Orbital_elements_h
# define Orbital_elements_h

#include <cmath>
#include "Constants.h"

struct orbel_t {
    value_type h[3];    // Angular Momentum
    value_type qlap[3]; // Laplace vector
    value_type easc[3]; // unit vector of the line of ascending nodes
    value_type Oman;
    value_type i;   // inclination
    value_type ecc; // Eccentricity
    value_type a;   // Separation
    value_type omperi ; // argument of periastron
    value_type Porb ;   // orbital period
    value_type norb ;   // Mean motion
    value_type v;       // True anomaly
    value_type E;       // Eccentric anomaly
    value_type m;       // Mean anomaly
    value_type tperi;   // Time of periastron passage 
    value_type tasc_approx; // Approximate time of ascending node (to zeroth order in eccentricity)
}; 

void printorbels(orbel_t orbel);


long double Solve_Kepler(long double dt, long double P, long double e, long double err = pow(10.L,-19)) ;


void orbel2statevect(value_type t, value_type e, value_type a1, value_type om1, value_type angli, value_type tperi1, value_type OrbPeriod, value_type om_an1 ,
                     value_type statevector1[6],
                     int jours=1, int solarframe=1) ;
                     
void orbel2statevect_1pn(value_type t, value_type et, value_type ar, value_type om1, value_type angli, value_type tperi1, value_type OrbPeriod, value_type om_an1 , 
                        value_type m1, value_type mc,
                     value_type statevector1[6],
                     int jours=1, int solarframe=1) ;
                     
// ##########################################################################################################
// Deals with arbitrary observer line of sight and hierarchical systems 
// These 4 routines are consistent one with another

void statevect2orbel(valarray<value_type> statevector, value_type mu, value_type t, valarray<value_type> dir_obs, orbel_t & orbel);

void statevect2orbel_nbody(int nbody, valarray<valarray<value_type>> statevectors, valarray<value_type> Ms, value_type t, valarray<value_type> dir_obs, orbel_t * orbels); // Jacobi coordinates

void orbel2statevects_bis(orbel_t orbel, valarray<value_type> dir_obs, valarray<value_type> & statevector, value_type & t); 

void orbel2statevects_bis_nbody(int nbody, orbel_t * orbels, valarray<value_type> Ms, valarray<value_type> dir_obs,  bool com_1PN, valarray<valarray<value_type>> & statevectors, value_type &  t); // Jacobi coordinates

// ##########################################################################################################
#endif
