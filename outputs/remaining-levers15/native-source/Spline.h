/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

# ifndef Spline_h
# define Spline_h

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 * 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */


#include <cmath>
#include <iostream>
#include "Constants.h"

   using namespace std;

   // Déclaration des types








class Spline

{    
    value_type * xs;
    value_type * ys;
    
    long int n;         // number of interpolation points = _xs.size
    
    value_type origin ; // Contains the shift internally applied to the abscissa _xs to avoid numerical problems

    public :

    value_type * y2; // computed values of the second derivative at interpolation points

    Spline( long int ) ; // Just allocate memory 
    Spline (value_type *, value_type *, long int ) ; // Allocate and spline 
    ~Spline(){delete[] xs ; delete[] ys ; delete[] y2; };
    
    void ReSpline(value_type [], value_type []) ; // Perform another spline with the same number of nodes as provided initially with the constructor. (identical to constructor but without memory alloc, "reconstruct")
    
    value_type operator () (const value_type&, const long int&, const long int& );
    
    value_type Integrate(const value_type& , const value_type& , 
                         const long int& , const long int& , const long int&, const long int& );
};


#endif
