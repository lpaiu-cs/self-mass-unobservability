/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 *
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>


// The numerical data type which is used within the stepper
typedef long double value_type;
typedef long long int turntype; // if this is changed need to change MPI commands refering to long long int
typedef float errortype;

const long double zero = 0.L;
const long double un = 1.L;
const long double deux = 2.L;
const long double trois = 3.L;
const long double quatre = 4.L;
const long double cinq = 5.L;
const long double six = 6.L;
const long double sept = 7.L;
const long double huit = 8.L;
const long double neuf = 9.L;
const long double dix = 10.L;
const long double douze = 12.L;
const long double undemi = 0.5L;
const long double deuxtiers = 2.L/3.L;
const long double troisdemis = 1.5L;
const long double cinqdemis = 2.5L;
const long double septdemis = 3.5L;
const long double neufdemis = 4.5L;
const long double unquart = 0.25L;

const long double deuxpi = 2.L * acos(-1.L) ;
const long double pi = acos(-1.L) ;


const long double clight = 299792458.L; // m/s
const long double clight2 = pow(clight,2); // m^2/s^2
const long double Ggrav = 6.67408e-11L;      // m^3 s^-2 kg^-1 (CODATA2014) ; previously used value //6.67384e-11L
const long double Msol = 1.9884754153381438e30L;   // Masse du Soleil (kg) (IAU 2015 Resolution B 3 + CODATA 2014) // previously used value 1.988435e30L;
const long double GMsol = Ggrav * Msol ;
const long double daysec = 86400.L; // number of sec by day
const long double yrsec = 86400.L * 365.25L; // number of sec by year
const long double yrday =  365.25L; // number of days per year
const long double raddeg = pi/180.L; // radians per degree
const long double radmasdeg = raddeg / (3600.L * 1000.L); // radian per mas of degree
const long double masdegrad = un/radmasdeg;
const long double radmash = deuxpi / (86400.L*1000.L); // radian per mas of hour
const long double lightyrmeter = clight * yrsec; // light year in meters

#endif
