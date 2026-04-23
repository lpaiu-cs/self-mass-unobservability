// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

// To compile : g++  -o Compare_parameters.exe Compare_parameters.cpp Parameters.cpp Utilities.cpp Orbital_elements.cpp

#include "Parameters.h"
#include "Utilities.h"




int main ( int argc, char *argv[] ) {


  char * parfile1;
  char * parfile2 ;
 // char * fileout ;

  int i = 0 ;
  
    if ( argc == 3 ) 
    {
        parfile1 = argv[1] ;
        parfile2 = argv[2];
       // fileout = argv[3];
    }
//     else
//     { 
//         cout << " Wrong numbers of arguments ! " << endl;
//          cout << "********************** Syntax : 'program parfile1 parfile1 '" << "************************" << endl;
//         //return 1;
//     }
    
    Parametres param1(parfile1);
    Parametres param2(parfile2);
    cout << "Variation of parameters in proportion of the initial parameter (p2 - p1)/abs(p1)" << endl;
    for (i = 0 ; i < param1.nfitparams ; i++) 
    {
        printf("%i %s %.19Le \n", i, param1.fitparams_names[i], (param2[i] - param1[i]) / abs(param1[i]) );
    }
    
    cout << endl << "Variation of parameters in proportion of the scale (p2 - p1)/abs(scale2)" << endl;
    for (i = 0 ; i < param1.nfitparams ; i++) 
    {
        printf("%i %s %.19Le \n", i, param1.fitparams_names[i], (param2[i] - param1[i]) / abs(param2.Parameter_scale(i)) );
    }
    
        cout << endl << "Variation of parameters in proportion of the scale (p2 - p1)/abs(scale1)" << endl;
    for (i = 0 ; i < param1.nfitparams ; i++) 
    {
        printf("%i %s %.19Le \n", i, param1.fitparams_names[i], (param2[i] - param1[i]) / abs(param1.Parameter_scale(i)) );
    }
    
    return 0;
   }

