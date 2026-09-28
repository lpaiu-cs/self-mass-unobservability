// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 *
 */

# include "Parameters.h"
#include "Diagnostics.h"
# include <cctype>
#include "Constants.h"
#include <cmath>
#include "Utilities.h"
#include <valarray>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "Orbital_elements.h"
#include <iostream>


using namespace std;

Parametres::Parametres(const char * filename)
{
    allocation_flag = false;
    Load_new_parfile( filename)   ;
    motion_changed = true;
    return ;
}

void Parametres::Set_fitted_parameter_map(const vector<int> list_of_fitted_parameters)
{
    fitted_parameters = list_of_fitted_parameters;
    return;
}


vector<int> Parametres::Get_fitted_parameter_map()
{
    return fitted_parameters;
}

int Parametres::Get_fitted_parameter_index(int absolute_index)
// Return the index of parameter absolute_index  in the fitted parameter list fitted_parameters.
// If this parameter is not fitted, returns -1
{
    int index = absolute_parameter_map[absolute_index];
    int i = 0;
    while (i < fitted_parameters.size() and fitted_parameters[i] != index) i++;
    if (i >= fitted_parameters.size())
        return -1;
    else
        return i;
}


void Parametres::ParamList2AbsParam(int paramnumber, int & absparamnumber, int & index)
/* Map: paramnumber in Parameter List -> (absolute parameter, index)
 * Basically inverts  absolute_parameter_map
 */
{
    int i = 0;
    while(i < absolute_parameter_map.size() and absolute_parameter_map[i] <= paramnumber)  i++;
    
    absparamnumber = i - 1 ;
    index = paramnumber - absolute_parameter_map[absparamnumber];
};


int Parametres::Is_in_parameter_set(int pset)
// Returns -1 if all fitted parameter are in parameter_sets[pset]. Otherwise returns the absolute index of the first fitted parameter not in the list.
{
  int absparam;
  int subindex;
  const int * psetarray = parameter_sets[pset];
  for (int i = 0; i < fitted_parameters.size(); i++)
      {
        ParamList2AbsParam(fitted_parameters[i], absparam, subindex);
        if (SearchArray<int>(absparam, psetarray, n_absparameters) < 0 )
        {
          return absparam;
        }
      }
    return -1;
}


value_type Parametres::Parameters(int paramnumber)
{
    if (paramnumber >= nfitparams or paramnumber < 0)
    {
        cout << "Error : Wrong parnumber in Parametres::Parameters(int paramnumber) " << endl ;
        return zero;
    };

    int absparamnumber=0;
    int index = 0;

    ParamList2AbsParam(paramnumber, absparamnumber, index);


    if (absparamnumber == 0) return f;
    if (absparamnumber == 1) return f1;

    if (absparamnumber == 2) return etap;
    if (absparamnumber == 3) return apsinii;
    if (absparamnumber == 4) return kappap;
    if (absparamnumber == 5) return apcosii;
    if (absparamnumber == 6) return tascp;
    if (absparamnumber == 7) return Pi;

    if (absparamnumber == 8) return etaB;
    if (absparamnumber == 9) return aBsinio;
    if (absparamnumber == 10) return kappaB;
    if (absparamnumber == 11) return aBcosio;
    if (absparamnumber == 12) return tascB;
    if (absparamnumber == 13) return Po;

    if (absparamnumber == 14) return mui;
    if (absparamnumber == 15) return muio;
    if (absparamnumber == 16) return Mo;

    if (absparamnumber == 17) return oman;
    if (absparamnumber == 18) return delta_oman;

    if (absparamnumber == 19) return dphase0;

    if (absparamnumber == 20) return SEP_D;
    if (absparamnumber == 21) return SEP_gamma ;
    if (absparamnumber == 22) return SEP_beta_0pp  ;
    if (absparamnumber == 23) return SEP_beta_p00 ;
    if (absparamnumber == 24) return SEP_beta_0p0 ;

    if (absparamnumber == 25) return RA;
    if (absparamnumber == 26) return DEC;
    if (absparamnumber == 27) return distance;
    if (absparamnumber == 28) return RA1;
    if (absparamnumber == 29) return DEC1;
    if (absparamnumber == 30) return distance1;

    if (absparamnumber == 31) return DM;
    if (absparamnumber == 32) return DM1;
    if (absparamnumber == 33) return DMX[index];
    if (absparamnumber == 34) return FD[index];

    if (absparamnumber == 35) return quadrupole[index];

    if (absparamnumber == 36) return efac ;

    if (absparamnumber == 37) return delta_i;
    
    if (absparamnumber == 38) return eta_extra[index];
    if (absparamnumber == 39) return kappa_extra[index];
    if (absparamnumber == 40) return asini_extra[index];
    if (absparamnumber == 41) return acosi_extra[index];
    if (absparamnumber == 42) return tasc_extra[index];
    if (absparamnumber == 43) return P_extra[index];
    if (absparamnumber == 44) return oman_extra[index];
    
    // Add parameter here


    return zero;

}



void Parametres::Set_parameters(int paramnumber, const value_type newvalue)
{
    if (paramnumber >= nfitparams or paramnumber < 0)
    {
        cout << "Error : Wrong paramnumber in Parametres::Set_parameters(int paramnumber, const value_type newvalue) " << endl ;
        return;
    };

    int absparamnumber=0;
    int index = 0;

    ParamList2AbsParam(paramnumber, absparamnumber, index);

    if (absparamnumber == 0) f = newvalue;
    if (absparamnumber == 1) f1 = newvalue ;

    if (absparamnumber == 2) etap = newvalue ;
    if (absparamnumber == 3) apsinii = newvalue ;
    if (absparamnumber == 4) kappap = newvalue ;
    if (absparamnumber == 5) apcosii = newvalue ;
    if (absparamnumber == 6) tascp = newvalue ;
    if (absparamnumber == 7) Pi = newvalue ;

    if (absparamnumber == 8) etaB = newvalue ;
    if (absparamnumber == 9) aBsinio = newvalue ;
    if (absparamnumber == 10) kappaB = newvalue ;
    if (absparamnumber == 11) aBcosio = newvalue ;
    if (absparamnumber == 12) tascB = newvalue ;
    if (absparamnumber == 13) Po = newvalue ;

    if (absparamnumber == 14) mui = newvalue ;
    if (absparamnumber == 15) muio = newvalue ;
    if (absparamnumber == 16) Mo = newvalue ;

    if (absparamnumber == 17) oman = newvalue ;
    if (absparamnumber == 18) delta_oman = newvalue ;

    if (absparamnumber == 19) dphase0 = newvalue ;

    if (absparamnumber == 20) SEP_D = newvalue ;
    if (absparamnumber == 21) SEP_gamma = newvalue ;
    if (absparamnumber == 22) SEP_beta_0pp = newvalue ;
    if (absparamnumber == 23) SEP_beta_p00 = newvalue ;
    if (absparamnumber == 24) SEP_beta_0p0 = newvalue ;

    if (absparamnumber == 25) RA = newvalue ;
    if (absparamnumber == 26) DEC = newvalue ;
    if (absparamnumber == 27) distance = newvalue;
    if (absparamnumber == 28) RA1 = newvalue ;
    if (absparamnumber == 29) DEC1 = newvalue ;
    if (absparamnumber == 30) distance1 = newvalue;

    if (absparamnumber == 31) DM = newvalue;
    if (absparamnumber == 32) DM1 = newvalue;
    if (absparamnumber == 33) DMX[index] = newvalue;
    if (absparamnumber == 34) FD[index] = newvalue;

    if (absparamnumber == 35) quadrupole[index] = newvalue;

    if (absparamnumber == 36) efac = newvalue ;

    if (absparamnumber == 37) delta_i = newvalue ;

    if (absparamnumber == 38) eta_extra[index] = newvalue;
    if (absparamnumber == 39) kappa_extra[index] = newvalue;
    if (absparamnumber == 40) asini_extra[index] = newvalue;
    if (absparamnumber == 41) acosi_extra[index] = newvalue;
    if (absparamnumber == 42) tasc_extra[index] = newvalue;
    if (absparamnumber == 43) P_extra[index] = newvalue;
    if (absparamnumber == 44) oman_extra[index] = newvalue;
    
    // Add parameter here
    return;
}



value_type Parametres::Parameter_scale(int paramnumber)
{// Return scale of parameter "paramnumber"
    if (paramnumber >= nfitparams or paramnumber < 0)
    {
        cout << "Error : Wrong parnumber in Parametres::Parameter_scale(int paramnumber) " << endl ;
        return zero;
    }

    int absparamnumber=0;
    int index = 0;

    ParamList2AbsParam(paramnumber, absparamnumber, index);

    if (absparamnumber == 0) return scale_f;
    if (absparamnumber == 1) return scale_f1;

    if (absparamnumber == 2) return scale_etap;
    if (absparamnumber == 3) return scale_apsinii;
    if (absparamnumber == 4) return scale_kappap;
    if (absparamnumber == 5) return scale_apcosii;
    if (absparamnumber == 6) return scale_tascp;
    if (absparamnumber == 7) return scale_Pi;

    if (absparamnumber == 8) return scale_etaB;
    if (absparamnumber == 9) return scale_aBsinio;
    if (absparamnumber == 10) return scale_kappaB;
    if (absparamnumber == 11) return scale_aBcosio;
    if (absparamnumber == 12) return scale_tascB;
    if (absparamnumber == 13) return scale_Po;

    if (absparamnumber == 14) return scale_mui;
    if (absparamnumber == 15) return scale_muio;
    if (absparamnumber == 16) return scale_Mo;

    if (absparamnumber == 17) return scale_oman;
    if (absparamnumber == 18) return scale_delta_oman;

    if (absparamnumber == 19) return scale_dphase0;

    if (absparamnumber == 20) return scale_SEP_D;
    if (absparamnumber == 21) return scale_SEP_gamma ;
    if (absparamnumber == 22) return scale_SEP_beta_0pp  ;
    if (absparamnumber == 23) return scale_SEP_beta_p00  ;
    if (absparamnumber == 24) return scale_SEP_beta_0p0  ;


    if (absparamnumber == 25) return scale_RA;
    if (absparamnumber == 26) return scale_DEC;
    if (absparamnumber == 27) return scale_distance;
    if (absparamnumber == 28) return scale_RA1;
    if (absparamnumber == 29) return scale_DEC1;
    if (absparamnumber == 30) return scale_distance1;

    if (absparamnumber == 31) return scale_DM;
    if (absparamnumber == 32) return scale_DM1;
    if (absparamnumber == 33) return scale_DMX[index];
    if (absparamnumber == 34) return scale_FD[index];

    if (absparamnumber == 35) return scale_quadrupole[index];

    if (absparamnumber == 36) return scale_efac  ;

    if (absparamnumber == 37) return scale_delta_i;

    if (absparamnumber == 38) return scale_eta_extra[index];
    if (absparamnumber == 39) return scale_kappa_extra[index];
    if (absparamnumber == 40) return scale_asini_extra[index];
    if (absparamnumber == 41) return scale_acosi_extra[index];
    if (absparamnumber == 42) return scale_tasc_extra[index];
    if (absparamnumber == 43) return scale_P_extra[index];
    if (absparamnumber == 44) return scale_oman_extra[index];
    
    // Add parameter here
}



void Parametres::Set_parameter_scale(int paramnumber, const value_type newvalue)
{
    if (paramnumber >= nfitparams or paramnumber < 0)
    {
        cout << "Error : Wrong parnumber in Parametres::Set_parameter_scale(int paramnumber, const value_type newvalue) " << endl ;
        return;
    }

    int absparamnumber=0;
    int index = 0;

    ParamList2AbsParam(paramnumber, absparamnumber, index);


    if (absparamnumber == 0) scale_f = newvalue;
    if (absparamnumber == 1) scale_f1 = newvalue ;

    if (absparamnumber == 2) scale_etap = newvalue ;
    if (absparamnumber == 3) scale_apsinii = newvalue ;
    if (absparamnumber == 4) scale_kappap = newvalue ;
    if (absparamnumber == 5) scale_apcosii = newvalue ;
    if (absparamnumber == 6) scale_tascp = newvalue ;
    if (absparamnumber == 7) scale_Pi = newvalue ;

    if (absparamnumber == 8) scale_etaB = newvalue ;
    if (absparamnumber == 9) scale_aBsinio = newvalue ;
    if (absparamnumber == 10) scale_kappaB = newvalue ;
    if (absparamnumber == 11) scale_aBcosio = newvalue ;
    if (absparamnumber == 12) scale_tascB = newvalue ;
    if (absparamnumber == 13) scale_Po = newvalue ;

    if (absparamnumber == 14) scale_mui = newvalue ;
    if (absparamnumber == 15) scale_muio = newvalue ;
    if (absparamnumber == 16) scale_Mo = newvalue ;

    if (absparamnumber == 17) scale_oman = newvalue ;
    if (absparamnumber == 18) scale_delta_oman = newvalue ;

    if (absparamnumber == 19) scale_dphase0 = newvalue ;

    if (absparamnumber == 20) scale_SEP_D = newvalue ;
    if (absparamnumber == 21) scale_SEP_gamma = newvalue ;
    if (absparamnumber == 22) scale_SEP_beta_0pp = newvalue ;
    if (absparamnumber == 23) scale_SEP_beta_p00 = newvalue ;
    if (absparamnumber == 24) scale_SEP_beta_0p0 = newvalue ;

    if (absparamnumber == 25) scale_RA = newvalue ;
    if (absparamnumber == 26) scale_DEC = newvalue ;
    if (absparamnumber == 27) scale_distance = newvalue;
    if (absparamnumber == 28) scale_RA1 = newvalue ;
    if (absparamnumber == 29) scale_DEC1 = newvalue ;
    if (absparamnumber == 30) scale_distance1 = newvalue;

    if (absparamnumber == 31) scale_DM = newvalue;
    if (absparamnumber == 32) scale_DM1 = newvalue;
    if (absparamnumber == 33) scale_DMX[index] = newvalue;
    if (absparamnumber == 34) scale_FD[index] = newvalue;

    if (absparamnumber == 35) scale_quadrupole[index] = newvalue;

    if (absparamnumber == 36) scale_efac = newvalue ;

    if (absparamnumber == 37) scale_delta_i = newvalue;
    
    if (absparamnumber == 38) scale_eta_extra[index];
    if (absparamnumber == 39) scale_kappa_extra[index];
    if (absparamnumber == 40) scale_asini_extra[index];
    if (absparamnumber == 41) scale_acosi_extra[index];
    if (absparamnumber == 42) scale_tasc_extra[index];
    if (absparamnumber == 43) scale_P_extra[index];
    if (absparamnumber == 44) scale_oman_extra[index];
    
    // Add parameter here
    return;
}


void Parametres::Load_new_parfile(const char * filename)
{
  int i = 0;
  int j,k;
  char snumber[30];
  char paramname_with_number[50];

  //****** Default values before reading parfile

    strcpy(t2parfile,"");
    strcpy(specialcase,"");

    // *** default values of fittable parameters
        // Fittable parameters          Index in array
    f = zero ;            //        0
    f1 = zero ;           //        1

    etap = zero ;         //        2
    apsinii = zero ;      //        3
    kappap = zero ;       //        4
    apcosii = zero ;      //        5
    tascp = zero ;        //        6
    Pi = zero ;           //        7

    etaB = zero ;         //        8
    aBsinio = zero ;      //        9
    kappaB = zero ;       //        10
    aBcosio = zero ;      //        11
    tascB = zero ;        //        12
    Po = zero ;           //        13

    mui = zero ;          //        14
    muio = zero ;         //        15
    Mo = zero ;           //        16

    oman = zero;          //        17
    delta_oman = zero ;   //        18

    dphase0 = zero ;      //        19

    SEP_D = zero ;        //        20
    SEP_gamma = zero;       //      21
    SEP_beta_0pp = zero;    //      22
    SEP_beta_p00 = zero;    //      23
    SEP_beta_0p0 = zero;    //      24

    RA = zero;            //        25
    DEC = zero;           //        26
    distance = pow(dix, 12) ;     //27 // by default set distance to 10^12 light years
    RA1 = 0.;            //         28
    DEC1 = 0.;           //         29
    distance1 = zero ;     //       30

    DM = zero;            //        31
    DM1 = zero;           //        32
    DMX.resize(0);        //        33
    DMXranges.resize(0);
    FD.resize(0);         //        34

    quadrupole.resize(0);   //      35

    efac = un;              //      36

    delta_i  = zero;        //      37
    
    eta_extra.resize(0);    //      38
    kappa_extra.resize(0);  //      39
    asini_extra.resize(0);  //      40
    acosi_extra.resize(0);  //      41
    tasc_extra.resize(0);   //      42
    P_extra.resize(0);      //      43
    oman_extra.resize(0);      //     44
    
    // Add parameter here

// For relative modification of parameters : newpar = oldpar + scale_par * dpar

        scale_f = un ;            //        0
        scale_f1 = un ;           //        1

        scale_etap = un ;         //        2
        scale_apsinii = un ;      //        3
        scale_kappap = un ;       //        4
        scale_apcosii = un ;      //        5
        scale_tascp = un ;        //        6
        scale_Pi = un ;           //        7

        scale_etaB = un ;         //        8
        scale_aBsinio = un ;      //        9
        scale_kappaB = un ;       //        10
        scale_aBcosio = un ;      //        11
        scale_tascB = un ;        //        12
        scale_Po = un ;           //        13

        scale_mui = un ;          //        14
        scale_muio = un ;         //        15
        scale_Mo = un ;           //        16

        scale_oman = un ;         //        17
        scale_delta_oman = un ;   //        18

        scale_dphase0 = un ;      //        19

        scale_SEP_D = un ;        //        20
        scale_SEP_gamma = un;       //    21
        scale_SEP_beta_0pp = un;    //    22
        scale_SEP_beta_p00 = un;    //    23
        scale_SEP_beta_0p0 = un;    //    24

        scale_RA = un;            //        25
        scale_DEC = un;           //        26
        scale_distance = 100 ;     //        27  // Big uncertainty on the distance (in light years)
        scale_RA1 = un;            //        28
        scale_DEC1 = un;           //        29
        scale_distance1 = un ;     //        30

        scale_DM = un ;           //        31
        scale_DM1 = un ;          //        32
        scale_DMX.resize(0);      //        33
        scale_FD.resize(0);      //         34

        scale_quadrupole.resize(0); //      35

        scale_efac = un;            //      36

        scale_delta_i = zero;       //      37
        
        scale_eta_extra.resize(0);    //      38
        scale_kappa_extra.resize(0);  //      39
        scale_asini_extra.resize(0);  //      40
        scale_acosi_extra.resize(0);  //      41
        scale_tasc_extra.resize(0);   //      42
        scale_P_extra.resize(0);      //      43
        scale_oman_extra.resize(0);      //     44
        
        // Add parameter here


   // *********** End of default values of fittable parameters

  // *************** Initialise absolute parameter map ***************
  // Add parameter here
    absolute_parameter_map.resize(n_absparameters+1); // "+1" such that absolute_parameter_map[n_absparameters] = absolute_parameter_map[n_absparameters-1] + size of parameter (n_absparameters-1)
    for (i=0; i < 33; i++)
    {
        absolute_parameter_map[i] = i;
    }
    //a parameter stored in a vector has the same fitting index as the next parameter if its size is zero
    for (i=33; i < 36; i++) absolute_parameter_map[i] = 33;
    absolute_parameter_map[36] = 33; // efac
    absolute_parameter_map[37] = 34; // delta_i
    for (i=38; i < 45; i++) absolute_parameter_map[i] = 35; // extra_... parameters 
    absolute_parameter_map[n_absparameters] = 35;    // Accounts for size of last parameters
    nfitparams = 35; // Sum of the "sizes" of each parameter (see also below)
  // *************** End absolute parameter map ***************


    // Delays parameters
    treference = zero ;
    posepoch = zero ;
    geometric = true ;
    kopeikin=true;
    shklovskii = false;
    roemer = false ;
    einstein  = true ;
    shapiro = true ;
    aberration = false ;

    truefreq = 1;
    DopplerF = un ;

    // Integrator and computation parameters
    interpsteps_per_period_i = 100 ;
    interp_margin = 1.;
    integrator_type = 0;
    tolint = pow(10.L, -16) ;

    // Deduced parameters
     Mp = 0. ;
     Mi = 0. ;
     ei = 0. ;
     eo = 0. ;
     omp = 0. ;
     omB  = 0. ;
     tperii = 0. ;
     tperio = 0. ;
     ap  = 0. ;
     aB = 0. ;
     anglii = 0. ;
     anglio = 0. ;
//      for (i = 0 ; i < 3 ; i++)
//      {
//         for (j = 0 ; j < 3 ; j++)
//         {
//             Gg[i][j] = un;;
//             gammabar[i][j] = zero;
//             for (k = 0 ; k < 3 ; k++)
//             {
//                 betabar[i][j][k] = zero;
//             }
//         }
//      }
    quadrupole_kgm2.resize(0);

   // Parameter set
    parameter_set = -1;

    Read_parfile(filename) ;

    // Initialises those parameters that depend on the number of extra bodies
    nbodies_plus_extra = 3+nextra;
    Gg.resize(nbodies_plus_extra);
    gammabar.resize(nbodies_plus_extra);
    betabar.resize(nbodies_plus_extra);
    for (i= 0 ; i < nbodies_plus_extra ; i++)
    {
        Gg[i].resize(nbodies_plus_extra, un);
        gammabar[i].resize(nbodies_plus_extra, zero);
        betabar[i].resize(nbodies_plus_extra);
        for (j= 0 ; j < nbodies_plus_extra ; j++)  betabar[i][j].resize(nbodies_plus_extra, zero);
    }

  // Initialize and check parameter set
  int inpset ;
    if (parameter_set > -1)
    {
      inpset = Is_in_parameter_set(parameter_set);
      if (inpset >= 0)
      {
        printf("\n WARNING : parameter with absolute index %i is not in the current parameter set (%d) !\n\n", inpset, parameter_set);
      }
    }

    // Apply timeshift to all time references

    timeshift = treference ;

    treference -= timeshift ;
    posepoch -= timeshift ;
    tascp -= timeshift ;
    tascB -= timeshift ;
    for (i = 0; i < nextra ; i++) tasc_extra[i] -= timeshift;

    // Max number of fittable parameters
    nfitparams = nfitparams + DMX.size() + FD.size() + quadrupole.size();
    nfitparams += eta_extra.size() + kappa_extra.size() + acosi_extra.size() + asini_extra.size() + P_extra.size() + tasc_extra.size() + oman_extra.size();
    // Add parameter here : if it is a vector parameter, its size must be part of the sum above

    if ( allocation_flag == false )
    {
            // Initialize parameter names
        fitparams_names = (char**) malloc(nfitparams * sizeof( char *) );
        for (i = 0 ; i< nfitparams ; i++) fitparams_names[i] = (char*) malloc(50 * sizeof( char ) );
        //for (int i = 0 ; i < nfitparams ; i++) fitparams_names[i] = (char*) malloc(100* sizeof( char ) );
        strcpy(fitparams_names[absolute_parameter_map[0]], "spinfreq");
        strcpy(fitparams_names[absolute_parameter_map[1]], "spinfreq1");
        strcpy(fitparams_names[absolute_parameter_map[2]], "eta_p");
        strcpy(fitparams_names[absolute_parameter_map[3]] , "apsini_i");
        strcpy(fitparams_names[absolute_parameter_map[4]] , "kappa_p");
        strcpy(fitparams_names[absolute_parameter_map[5]] , "apcosi_i");
        strcpy(fitparams_names[absolute_parameter_map[6]] , "tasc_p");
        strcpy(fitparams_names[absolute_parameter_map[7]] , "period_i");
        strcpy(fitparams_names[absolute_parameter_map[8]] , "eta_b");
        strcpy(fitparams_names[absolute_parameter_map[9]] , "absini_o");
        strcpy(fitparams_names[absolute_parameter_map[10]] , "kappa_b");
        strcpy(fitparams_names[absolute_parameter_map[11]] , "abcosi_o");
        strcpy(fitparams_names[absolute_parameter_map[12]] , "tasc_b");
        strcpy(fitparams_names[absolute_parameter_map[13]] , "period_o");
        strcpy(fitparams_names[absolute_parameter_map[14]] , "masspar_p");
        strcpy(fitparams_names[absolute_parameter_map[15]] , "masspar_i");
        strcpy(fitparams_names[absolute_parameter_map[16]] , "mass_o");
        strcpy(fitparams_names[absolute_parameter_map[17]] , "oman");
        strcpy(fitparams_names[absolute_parameter_map[18]] , "deltaoman");
        strcpy(fitparams_names[absolute_parameter_map[19]] , "dphase0");
        strcpy(fitparams_names[absolute_parameter_map[20]] , "SEP_D");
        strcpy(fitparams_names[absolute_parameter_map[21]] , "SEP_gamma");
        strcpy(fitparams_names[absolute_parameter_map[22]] , "SEP_beta_0pp");
        strcpy(fitparams_names[absolute_parameter_map[23]] , "SEP_beta_p00");
        strcpy(fitparams_names[absolute_parameter_map[24]] , "SEP_beta_0p0");
        strcpy(fitparams_names[absolute_parameter_map[25]] , "right_ascension");
        strcpy(fitparams_names[absolute_parameter_map[26]] , "declination");
        strcpy(fitparams_names[absolute_parameter_map[27]] , "distance");
        strcpy(fitparams_names[absolute_parameter_map[28]] , "right_ascension1");
        strcpy(fitparams_names[absolute_parameter_map[29]] , "declination1");
        strcpy(fitparams_names[absolute_parameter_map[30]] , "distance1");
        strcpy(fitparams_names[absolute_parameter_map[31]] , "dispmeasure");
        strcpy(fitparams_names[absolute_parameter_map[32]] , "dispmeasure1");
        for (i = 0; i  <  DMX.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "dispmeasure_X") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[33]+i] , paramname_with_number);
            };
        for (i = 0; i  <  FD.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "dispmeasure_FD") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[34]+i] , paramname_with_number);
            };
        for (i = 0; i  <  quadrupole.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "quadrupole") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[35]+i] , paramname_with_number);
            };
        strcpy(fitparams_names[absolute_parameter_map[36]] , "efac");
        strcpy(fitparams_names[absolute_parameter_map[37]] , "delta_i");
        for (i = 0; i  <  eta_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "eta_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[38]+i] , paramname_with_number);
            };
        for (i = 0; i  <  kappa_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "kappa_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[39]+i] , paramname_with_number);
            };
        for (i = 0; i  <  asini_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "asini_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[40]+i] , paramname_with_number);
            };
        for (i = 0; i  <  acosi_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "acosi_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[41]+i] , paramname_with_number);
            };
        for (i = 0; i  <  tasc_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "tasc_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[42]+i] , paramname_with_number);
            };
        for (i = 0; i  <  P_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "P_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[43]+i] , paramname_with_number);
            };
        for (i = 0; i  <  oman_extra.size() ; i++)
            {
                sprintf(snumber, "%d", i+1);
                strcpy(paramname_with_number, "oman_extra") ;
                strcat(paramname_with_number, snumber);
                strcpy(fitparams_names[absolute_parameter_map[44]+i] , paramname_with_number);
            };
        // Add parameter here
        };


        parameters_ini  = (value_type*) malloc( nfitparams * sizeof(value_type ) );

    // Adjusts size of derived-parameter  for extra bodies
        a_extra.resize(nextra);
        e_extra.resize(nextra);
        om_extra.resize(nextra);
        angli_extra.resize(nextra);
        tperi_extra.resize(nextra);
        M_extra.resize(nextra);

                
        allocation_flag = true;        
    Print();
    Set_reference_to_current_parameters();

    return;
}


void Parametres::Set_reference_to_current_parameters()
{
    for (int i = 0 ; i < nfitparams ; ++i ) parameters_ini[i] = Parameters(i);
    return;
}

void Parametres::Reset_to_initial_parameters()
{
    for (int i = 0 ; i < nfitparams ; ++i ) Set_parameters(i, parameters_ini[i]);
    return;
}


void Parametres::Insert_in_parameter_maps(int absolute_parameter_map_index, int index )
// Update the paramater maps ( fittable, and fitted) when a new parameter is inserted at absolute parameter map index  absolute_parameter_map_index, and subindex index.
// No need to change this routine to add a parameter
{
    int shift_absparam=0;
    int min_shifted_index=0;
    int j=0;

    if (absolute_parameter_map_index < absolute_parameter_map.size() )
    {
        shift_absparam = index - (absolute_parameter_map[absolute_parameter_map_index+1] - absolute_parameter_map[absolute_parameter_map_index]);
        min_shifted_index = absolute_parameter_map[absolute_parameter_map_index+1];
    }
    else
    {
        shift_absparam = 0;
        min_shifted_index = nfitparams;
    }

    // Update the fitted-parameter map
    for (j=0; j < fitted_parameters.size() ; j++)
    {
        if (fitted_parameters[j] >= min_shifted_index ) fitted_parameters[j] += shift_absparam;
    };

    // Update the absolute parameter map
    for (j=absolute_parameter_map_index + 1; j < absolute_parameter_map.size() ; j++)
    {
        absolute_parameter_map[j] += shift_absparam;
    };
};


void Parametres::Read_parfile(const char * filename)
{
    /*To add a parameter, add a block like the following :
     *
     * else if (str.find("paramname") != string::npos){
            sscanf( line, "%s %Le %Le ", textin, &param, &scale_param, &active );
            param *= convert_coeff // convert to the internal system of units if necessary
            if ( active == 1 ) fitted_parameters.push_back(param number);
        }
     * Where paramname should be as in fitparams_names, param number as defined by convention in the header.
     * Be careful with name redundances : if the "param name " pattern appears in another name it might be confused (this would need improvement), like for "spinfreq" and "spinfreq1".
     */

    FILE *fp;
    const int sizechar = 100;
    //int returnvalue = 0;
    //value_type val, dval = zero;
    int ival ;
    int i, j;
    char line[sizechar];
    char car[1];
    char textin[sizechar]; // Should be the same size as line otherwise sscanf fails
    string str ;

    int active = 1;
    int shift_absparam =0;
    int min_shifted_index = 0;
    int index = 0;
    value_type xread = zero;
    value_type xread_scale = zero;
    
    fitted_parameters.resize(0);

    fp = fopen(filename,"r");
    if (fp == NULL)
    {
        printf("\n\n Opening of parfile %s failed ! \n\n", filename);
        return;
    }

    while ( fgets(line, sizechar, fp) != NULL ){
        active = 1;
        

        if (sscanf( line, "%s ", textin) == EOF )  continue ;
        str = textin ;

        if (str.find("tempo2parfile") != string::npos) {
            sscanf( line, "%s %s", textin, t2parfile );
        }
        
        if (str.find("specialcase") != string::npos) {
            sscanf( line, "%s %s", textin, specialcase );
        }
        
         if (sscanf( line, "spinfreq1 %Le %Le %i", &f1, &scale_f1, &active  ) == 3) {
            f1 *= pow(daysec,2);
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[1]);
        }
        else  if (sscanf( line, "spinfreq %Le %Le %i", &f, &scale_f, &active  ) == 3){ // ! Be careful with the order !
                f *= daysec;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[0]);
        }
        if (sscanf( line, "eta_p %Le %Le %i", &etap, &scale_etap, &active  ) == 3){ //str.find("eta_p") != string::npos){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[2]);
        }
        if (str.find("apsini_i") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &apsinii, &scale_apsinii , &active );
                apsinii *= clight;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[3]);
        }
        if (str.find("kappa_p") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &kappap, &scale_kappap , &active );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[4]);
        }
        if (str.find("apcosi_i") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &apcosii, &scale_apcosii , &active );
                apcosii *= clight;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[5]);
        }
        if (str.find("tasc_p") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &tascp, &scale_tascp , &active );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[6]);
        }
        if (str.find("period_i") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &Pi, &scale_Pi, &active  );
                if ( active == 1 )
                {
                    fitted_parameters.push_back(absolute_parameter_map[7]);
                }
        }
        if (str.find("eta_b") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &etaB, &scale_etaB, &active  );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[8]);
        }
        if (str.find("absini_o") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &aBsinio, &scale_aBsinio, &active  );
                aBsinio *= clight;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[9]);
        }
        if (str.find("kappa_b") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &kappaB, &scale_kappaB, &active  );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[10]);
        }
        if (str.find("abcosi_o") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &aBcosio, &scale_aBcosio, &active  );
                aBcosio *= clight;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[11]);
        }
        if (str.find("tasc_b") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &tascB, &scale_tascB, &active  );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[12]);
        }
        if (str.find("period_o") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &Po, &scale_Po, &active  );
                if ( active == 1 )
                {
                    fitted_parameters.push_back(absolute_parameter_map[13]);
                }
        }
        if (str.find("masspar_p") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &mui, &scale_mui, &active  );
                if ( active == 1 )
                {
                    fitted_parameters.push_back(absolute_parameter_map[14]);
                }
        }
        if (str.find("masspar_i") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &muio, &scale_muio, &active  );
                if ( active == 1 )
                {
                    fitted_parameters.push_back(absolute_parameter_map[15]);
                }
        }
        if (str.find("mass_o") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &Mo, &scale_Mo, &active  );
                if ( active == 1 )
                {
                    fitted_parameters.push_back(absolute_parameter_map[16]);
                }
        }
        if (str.find("deltaoman") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &delta_oman, &scale_delta_oman, &active  );
                delta_oman *= raddeg;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[18]);
        }
        else if (str.find("oman") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &xread, &xread_scale, &active  );
                if (strcmp(textin, "oman") == 0 ) // check in particular that it is not an "oman_extra" line
                {
                    oman = xread * raddeg;
                    scale_oman = xread_scale;
                    if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[17]);
                }
        }
        if (str.find("dphase0") != string::npos){
                sscanf( line, "%s %Le %Le %i", textin, &dphase0, &scale_dphase0, &active  );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[19]);
        }
        if (str.find("SEP_D") != string::npos){
                sscanf( line, "%s %Le %Le  %i", textin, &SEP_D, &scale_SEP_D, &active  );
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[20]);
        }
        if (sscanf( line, "SEP_gamma %Le %Le %i", &SEP_gamma, &scale_SEP_gamma, &active  ) == 3){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[21]);
        }
        if (sscanf( line, "SEP_beta_0pp %Le %Le %i", &SEP_beta_0pp, &scale_SEP_beta_0pp, &active  ) == 3){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[22]);
        }
         if (sscanf( line, "SEP_beta_p00 %Le %Le %i", &SEP_beta_p00, &scale_SEP_beta_p00, &active  ) == 3){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[23]);
        }
        if (sscanf( line, "SEP_beta_0p0 %Le %Le %i", &SEP_beta_0p0, &scale_SEP_beta_0p0, &active  ) == 3){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[24]);
        }
        if (str.find("right_ascension1") != string::npos){
            sscanf( line, "%s %Le %Le %i", textin, &RA1, &scale_RA1, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[28]);
        }
        else if (str.find("right_ascension") != string::npos){
            sscanf( line, "%s %Le %Le %i", textin, &RA, &scale_RA, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[25]);
        }
        if (str.find("declination1") != string::npos){
            sscanf( line, "%s %Le %Le %i ", textin, &DEC1, &scale_DEC1, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[29]);
        }
        else if (str.find("declination") != string::npos){
            sscanf( line, "%s %Le %Le %i", textin, &DEC, &scale_DEC, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[26]);
        }
        if (str.find("distance1") != string::npos){
            sscanf( line, "%s %Le %Le %i", textin, &distance1, &scale_distance1, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[30]);
        }
        else if (str.find("distance") != string::npos){
            sscanf( line, "%s %Le %Le %i", textin, &distance, &scale_distance, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[27]);
        }
        if (str.find("dispmeasure1") != string::npos){
            sscanf( line, "%s %Le %Le %i", textin, &DM1, &scale_DM1, &active );
            if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[32]);
        }
        else if (str.find("dispmeasure") != string::npos){
            if (str.find("dispmeasure_X") != string::npos){
                sscanf(line, "dispmeasure_X%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > DMX.size())
                {
                    DMX.resize(index, zero );
                    scale_DMX.resize(index, un);
                    DMXranges.resize(index+1,0);
                    Insert_in_parameter_maps(33, index);
                };
                DMX[index-1] = xread;
                scale_DMX[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[33] + index - 1 );
            }
            else if (str.find("dispmeasure_FD") != string::npos){
                sscanf(line, "dispmeasure_FD%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > FD.size())
                {
                    FD.resize(index, zero );
                    scale_FD.resize(index, un);
                    Insert_in_parameter_maps(34, index );
                    printf("fd size %i\n", FD.size());
                    printf("fd size %i\n", FD.size());
                };
                FD[index-1] = xread;
                scale_FD[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[34] + index - 1 );

            }
            else
            {
                sscanf( line, "%s %Le %Le %i", textin, &xread, &xread_scale, &active );
                DM = xread;
                scale_DM = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[31]);
            }
        }
        if (str.find("quadrupole") != string::npos){
                sscanf(line, "quadrupole%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > 0)
                {
                    if (index > quadrupole.size())
                    {
                        quadrupole.resize(index, zero );
                        scale_quadrupole.resize(index, un);
                        Insert_in_parameter_maps(35, index );
                    };
                    quadrupole[index-1] = xread;
                    scale_quadrupole[index-1] = xread_scale;
                    if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[35] + index - 1 );
                } else
                {
                    printf("\nError ! quadrupole index invalid : %i\n\n", index);
                }
        }
        if (sscanf( line, "efac %Le %Le %i", &efac, &scale_efac, &active  ) == 3){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[36]);
        }
        if (sscanf( line, "delta_i %Le %Le %i", &delta_i, &scale_delta_i, &active  ) == 3){
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[37]);
        }
        if (str.find("eta_extra") != string::npos){
                sscanf(line, "eta_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > eta_extra.size())
                {
                    eta_extra.resize(index, zero );
                    scale_eta_extra.resize(index, un);
                    Insert_in_parameter_maps(38, index);
                };
                eta_extra[index-1] = xread;
                scale_eta_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[38] + index - 1 );
        }
        if (str.find("kappa_extra") != string::npos){
                sscanf(line, "kappa_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > kappa_extra.size())
                {
                    kappa_extra.resize(index, zero );
                    scale_kappa_extra.resize(index, un);
                    Insert_in_parameter_maps(39, index);
                };
                kappa_extra[index-1] = xread;
                scale_kappa_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[39] + index - 1 );
        }
        
        if (str.find("asini_extra") != string::npos){
                sscanf(line, "asini_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > asini_extra.size())
                {
                    asini_extra.resize(index, zero );
                    scale_asini_extra.resize(index, un);
                    Insert_in_parameter_maps(40, index);
                };
                asini_extra[index-1] = xread * clight;
                scale_asini_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[40] + index - 1 );
        }
        
        if (str.find("acosi_extra") != string::npos){
                sscanf(line, "acosi_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > acosi_extra.size())
                {
                    acosi_extra.resize(index, zero );
                    scale_acosi_extra.resize(index, un);
                    Insert_in_parameter_maps(41, index);
                };
                acosi_extra[index-1] = xread * clight;
                scale_acosi_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[41] + index - 1 );
        }
        if (str.find("tasc_extra") != string::npos){
                sscanf(line, "tasc_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > tasc_extra.size())
                {
                    tasc_extra.resize(index, zero );
                    scale_tasc_extra.resize(index, un);
                    Insert_in_parameter_maps(42, index);
                };
                tasc_extra[index-1] = xread;
                scale_tasc_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[42] + index - 1 );
        }
        if (str.find("P_extra") != string::npos){
                sscanf(line, "P_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > P_extra.size())
                {
                    P_extra.resize(index, zero );
                    scale_P_extra.resize(index, un);
                    Insert_in_parameter_maps(43, index);
                };
                P_extra[index-1] = xread;
                scale_P_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[43] + index - 1 );
        }
        if (str.find("oman_extra") != string::npos){
                sscanf(line, "oman_extra%i %Le %Le %i", &index, &xread, &xread_scale, &active );
                if (index > oman_extra.size())
                {
                    oman_extra.resize(index, zero );
                    scale_oman_extra.resize(index, un);
                    Insert_in_parameter_maps(44, index);
                };
                oman_extra[index-1] = xread * raddeg;
                scale_oman_extra[index-1] = xread_scale;
                if ( active == 1 ) fitted_parameters.push_back(absolute_parameter_map[44] + index - 1 );
        }
        // Add parameter here : "else if ... " beware if a parameter name is included in another parameter name you need to take care of the reading order. See the example of "dispmeasure" or "spinfreq".
        if (str.find("treference") != string::npos){
            sscanf( line, "%s %Le ", textin, &treference );
        }
        if (str.find("geometric") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            geometric = ival;
        }
        if (str.find("shklovskii") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            shklovskii = ival;
        }
        if (str.find("kopeikin") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            kopeikin = ival;
        }
        if (str.find("roemer") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            roemer = ival;
        }
        if (str.find("einstein") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            einstein = ival;
        }
        if (str.find("shapiro") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            shapiro = ival;
        }
        if (str.find("aberration") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            aberration = ival;
        }
        if (str.find("truefreq") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            truefreq = ival;
        }
        if (str.find("interpsteps") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            interpsteps_per_period_i = ival;
        }
        if (str.find("interpolation_margin") != string::npos){
            sscanf( line, "%s %f ", textin, &interp_margin );
        }
        if (str.find("integrator") != string::npos){
            sscanf( line, "%s %i ", textin, &ival );
            integrator_type = ival;
        }
        if (str.find("integration_tolerance") != string::npos){
            sscanf( line, "%s %Le ", textin, &tolint );
        }
        if (str.find("posepoch") != string::npos){
            sscanf(line, "%s %Le ", textin, &posepoch);
        }
        if (str.find("parameter_set") != string::npos){
            sscanf(line, "%s %i ", textin, &parameter_set);
        }

    }
    fclose(fp);

    
    
    // Check consistency of the "extra" parameters
    nextra = acosi_extra.size();
    min(nextra = min(min(min(min(min(nextra, static_cast<int>(kappa_extra.size())), static_cast<int>(eta_extra.size())), static_cast<int>(tasc_extra.size())), static_cast<int>(P_extra.size())), static_cast<int>(oman_extra.size())), static_cast<int>(asini_extra.size()));

    if (nextra != kappa_extra.size() || nextra != asini_extra.size() || nextra != eta_extra.size() || nextra != tasc_extra.size() || nextra != P_extra.size() || nextra != oman_extra.size())
    {
        printf("\n  Warning : Some 'extra_' parameters are missing. Only %d extra bodies will be taken into account\n\n", nextra);
    }
    // End of consistency check of "extra" parameters
    

    return;
}


void Parametres::Save_parfile(const char * filename)
{
    FILE *fp;
    value_type inter = 0.;

    fp = fopen(filename,"w");

    vector<int> fitstatus(nfitparams, 0);
    for (int i = 0 ; i < fitted_parameters.size() ; i++)
    {
        fitstatus[ fitted_parameters[i] ] = 1;
    }

    fprintf( fp, "# Tempo2 parfile\n");
    fprintf( fp, "%s %s \n", "tempo2parfile", t2parfile ) ;

    fprintf( fp, "\n# Time of reference\n");
    fprintf( fp, "%s %.19Le \n", "treference", treference + timeshift ) ;
    fprintf( fp, "%s %.19Le \n", "posepoch", posepoch + timeshift ) ;

    fprintf( fp, "\n# Enabled delays\n");
    fprintf( fp, "%s %i \n", "geometric", int(geometric) ) ;
    fprintf( fp, "%s %i \n", "kopeikin", int(kopeikin) ) ;
    fprintf( fp, "%s %i \n", "shklovskii", int(shklovskii) ) ;
    fprintf( fp, "%s %i \n", "roemer",  int(roemer) ) ;
    fprintf( fp, "%s %i \n", "einstein", int(einstein) ) ;
    fprintf( fp, "%s %i \n", "shapiro", int(shapiro) ) ;
    fprintf( fp, "%s %i \n", "aberration", int(aberration) ) ;
    fprintf( fp, "%s %i \n", "truefreq", int(truefreq) ) ;

    fprintf( fp, "\n# Integrator parameters \n");
    fprintf( fp, "%s %i \n", "integrator", integrator_type ) ;
    fprintf( fp, "%s %.19Le \n", "integration_tolerance", tolint ) ;

    fprintf( fp, "\n# Interpolation parameters \n");
    fprintf( fp, "%s %i \n", "interpsteps", interpsteps_per_period_i ) ;
    fprintf( fp, "%s %f \n", "interpolation_margin", interp_margin ) ;

    fprintf( fp, "\n# Parameter set \n");
    fprintf( fp, "%s %i \n", "parameter_set", parameter_set ) ;
    
    fprintf( fp, "\n# Special case \n");
    fprintf( fp, "%s %s \n", "specialcase", specialcase ) ;

    fprintf( fp, "\n# ***** Fittable parameters ****\n");
    fprintf( fp, "\n# Position in the sky\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "right_ascension", RA, scale_RA, fitstatus[absolute_parameter_map[25]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "declination", DEC, scale_DEC , fitstatus[absolute_parameter_map[26]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "distance", distance, scale_distance , fitstatus[absolute_parameter_map[27]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "right_ascension1", RA1, scale_RA1, fitstatus[absolute_parameter_map[28]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "declination1", DEC1, scale_DEC1 , fitstatus[absolute_parameter_map[29]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "distance1", distance1, scale_distance1 , fitstatus[absolute_parameter_map[30]]) ;

    fprintf( fp, "\n# Dispersion Measure\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "dispmeasure",  DM, scale_DM, fitstatus[absolute_parameter_map[31]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "dispmeasure1", DM1, scale_DM1 , fitstatus[absolute_parameter_map[32]]) ;
    if (DMX.size() >= 1)
    {
        for (int i = 0 ; i < DMX.size() ; i++)
        {
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "dispmeasure_X",i+1,  DMX[i], scale_DMX[i], fitstatus[absolute_parameter_map[33] + i]) ;
        };
    };

    if (FD.size() >= 1)
    {
        for (int i = 0 ; i < FD.size() ; i++)
        {
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "dispmeasure_FD",i+1,  FD[i], scale_FD[i], fitstatus[absolute_parameter_map[34] + i]) ;
        };
    };


    fprintf( fp, "\n#Pulsar system parameters \n");
    fprintf( fp, "%s %.19Le  %.19Le %i \n", "spinfreq", Parameters(absolute_parameter_map[0])/daysec, Parameter_scale(absolute_parameter_map[0]) , fitstatus[absolute_parameter_map[0]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "spinfreq1", Parameters(absolute_parameter_map[1])/pow(daysec,2), Parameter_scale(absolute_parameter_map[1]) , fitstatus[absolute_parameter_map[1]]) ;
    fprintf( fp, "\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "eta_p", Parameters(absolute_parameter_map[2]), Parameter_scale(absolute_parameter_map[2]) , fitstatus[absolute_parameter_map[2]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "apsini_i", Parameters(absolute_parameter_map[3])/clight, Parameter_scale(absolute_parameter_map[3]) , fitstatus[absolute_parameter_map[3]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "kappa_p", Parameters(absolute_parameter_map[4]), Parameter_scale(absolute_parameter_map[4]) , fitstatus[absolute_parameter_map[4]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "apcosi_i", Parameters(absolute_parameter_map[5])/clight, Parameter_scale(absolute_parameter_map[5]) , fitstatus[absolute_parameter_map[5]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "delta_i", Parameters(absolute_parameter_map[37]), Parameter_scale(absolute_parameter_map[37]) , fitstatus[absolute_parameter_map[37]]) ;

    inter = Parameters(absolute_parameter_map[6]) + timeshift ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "tasc_p", inter, Parameter_scale(absolute_parameter_map[6]) , fitstatus[absolute_parameter_map[6]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "period_i", Parameters(absolute_parameter_map[7]), Parameter_scale(absolute_parameter_map[7]) , fitstatus[absolute_parameter_map[7]]) ;
    fprintf( fp, "\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "eta_b", Parameters(absolute_parameter_map[8]), Parameter_scale(absolute_parameter_map[8]) , fitstatus[absolute_parameter_map[8]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "absini_o", Parameters(absolute_parameter_map[9])/clight, Parameter_scale(absolute_parameter_map[9]) , fitstatus[absolute_parameter_map[9]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "kappa_b", Parameters(absolute_parameter_map[10]), Parameter_scale(absolute_parameter_map[10]) , fitstatus[absolute_parameter_map[10]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "abcosi_o", Parameters(absolute_parameter_map[11])/clight, Parameter_scale(absolute_parameter_map[11]) , fitstatus[absolute_parameter_map[11]]) ;
    inter = Parameters(absolute_parameter_map[12]) + timeshift ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "tasc_b", inter, Parameter_scale(absolute_parameter_map[12]) , fitstatus[absolute_parameter_map[12]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "period_o", Parameters(absolute_parameter_map[13]), Parameter_scale(absolute_parameter_map[13]) , fitstatus[absolute_parameter_map[13]]) ;
    fprintf( fp, "\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "masspar_p", Parameters(absolute_parameter_map[14]), Parameter_scale(absolute_parameter_map[14]) , fitstatus[absolute_parameter_map[14]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "masspar_i", Parameters(absolute_parameter_map[15]), Parameter_scale(absolute_parameter_map[15]), fitstatus[absolute_parameter_map[15]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "mass_o", Parameters(absolute_parameter_map[16]), Parameter_scale(absolute_parameter_map[16]), fitstatus[absolute_parameter_map[16]]) ;
    fprintf( fp, "\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "oman", Parameters(absolute_parameter_map[17])/raddeg, Parameter_scale(absolute_parameter_map[17]) , fitstatus[absolute_parameter_map[17]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "deltaoman", Parameters(absolute_parameter_map[18])/raddeg, Parameter_scale(absolute_parameter_map[18]) , fitstatus[absolute_parameter_map[18]]) ;
    fprintf( fp, "\n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "dphase0", Parameters(absolute_parameter_map[19]), Parameter_scale(absolute_parameter_map[19]) , fitstatus[absolute_parameter_map[19]]) ;
    fprintf( fp, "\n");
    fprintf( fp, "\n#SEP violation parameters \n");
    fprintf( fp, "%s %.19Le %.19Le %i \n", "SEP_D", Parameters(absolute_parameter_map[20]), Parameter_scale(absolute_parameter_map[20]) , fitstatus[absolute_parameter_map[20]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "SEP_gamma", Parameters(absolute_parameter_map[21]), Parameter_scale(absolute_parameter_map[21]) , fitstatus[absolute_parameter_map[21]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "SEP_beta_0pp", Parameters(absolute_parameter_map[22]), Parameter_scale(absolute_parameter_map[22]) , fitstatus[absolute_parameter_map[22]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "SEP_beta_p00", Parameters(absolute_parameter_map[23]), Parameter_scale(absolute_parameter_map[23]) , fitstatus[absolute_parameter_map[23]]) ;
    fprintf( fp, "%s %.19Le %.19Le %i \n", "SEP_beta_0p0", Parameters(absolute_parameter_map[24]), Parameter_scale(absolute_parameter_map[24]) , fitstatus[absolute_parameter_map[24]]) ;
    fprintf( fp, "\n");
    if (quadrupole.size() >= 1)
    {
        fprintf( fp, "\n#Quadrupole parameters \n");
        for (int i = 0 ; i < quadrupole.size() ; i++)
        {
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "quadrupole",i+1,  quadrupole[i], scale_quadrupole[i], fitstatus[absolute_parameter_map[35] + i]) ;
        };
    };
    fprintf( fp, "\n");
    if (asini_extra.size() >= 1)
    {
        fprintf( fp, "\n#Extra bodies \n");
        for (int i = 0 ; i < asini_extra.size() ; i++)
        {
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "eta_extra",i+1,  eta_extra[i], scale_eta_extra[i], fitstatus[absolute_parameter_map[38] + i]) ;
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "kappa_extra",i+1,  kappa_extra[i], scale_kappa_extra[i], fitstatus[absolute_parameter_map[39] + i]) ;
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "asini_extra",i+1,  asini_extra[i]/clight, scale_asini_extra[i], fitstatus[absolute_parameter_map[40] + i]) ;
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "acosi_extra",i+1,  acosi_extra[i]/clight, scale_acosi_extra[i], fitstatus[absolute_parameter_map[41] + i]) ;
            inter =  tasc_extra[i] + timeshift ;
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "tasc_extra",i+1, inter, scale_tasc_extra[i], fitstatus[absolute_parameter_map[42] + i]) ;
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "P_extra",i+1,  P_extra[i], scale_P_extra[i], fitstatus[absolute_parameter_map[43] + i]) ;
            fprintf( fp, "%s%i %.19Le %.19Le %i \n", "oman_extra",i+1,  oman_extra[i]/raddeg, scale_oman_extra[i], fitstatus[absolute_parameter_map[44] + i]) ;
            fprintf( fp, "\n");
        };
    }
    fprintf( fp, "\n#Others \n");
    fprintf( fp, "%s %.19Le %.19Le %i", "efac", Parameters(absolute_parameter_map[36]), Parameter_scale(absolute_parameter_map[36]) , fitstatus[absolute_parameter_map[36]]) ;
    
    // Add parameter here

    fclose(fp);

    return ;
}




void Parametres::Print()
{
    int i = 0;
    vector<int> fitstatus(nfitparams, 0);
    for (i = 0 ; i < fitted_parameters.size() ; i++)
    {
        fitstatus[ fitted_parameters[i] ] = 1;
    }
    printf("%s %s \n", "Tempo 2 parfile :", t2parfile) ;
    printf("%s %.19Le \n", "Time of reference, treference + timeshift :", treference + timeshift ) ;
    printf("%s %.19Le \n", "Time of reference for position in the sky, posepoch + timeshift :", posepoch + timeshift ) ;

    printf("%s %i \n", "geometric", int(geometric) ) ;
    printf("%s %i \n", "kopeikin", int(kopeikin) ) ;
    printf("%s %i \n", "shklovskii", int(shklovskii) ) ;
    printf("%s %i \n", "roemer", int(roemer) ) ;
    printf("%s %i \n", "einstein", int(einstein) ) ;
    printf("%s %i \n", "shapiro", int(shapiro) ) ;
    printf("%s %i \n", "aberration", int(aberration) ) ;
    printf("%s %i \n", "truefreq", int(truefreq) ) ;
    
    printf("%s '%s' \n", "Special case", specialcase ) ;
    
    printf("%s %i \n", "integrator", integrator_type ) ;
    printf("%s %i \n", "interpsteps_per_period_i", interpsteps_per_period_i ) ;
    printf("%s %f \n", "interp_margin", interp_margin ) ;
    printf("%s %.19Le \n", "Tolerance level of integrator, tolint : ", tolint ) ;

    printf("%s %.19Le \n", "Timeshift applied to times internally, timeshift : ", timeshift );

    printf("%s %i \n", "Parameter set (0 : Pi and Po fitted ; 1 : masspar_p and mo fitted) :", parameter_set ) ;

    printf("\n*** Fittable parameters ***\n");
    for (i = 0 ; i < nfitparams ; i++)
    {
        if (fitstatus[i] == 1)
            printf( "%i : %s %.19Le %.19Le fitted \n", i, fitparams_names[i], Parameters(i), Parameter_scale(i) ) ;
        else
            printf( "%i : %s %.19Le %.19Le \n", i, fitparams_names[i], Parameters(i), Parameter_scale(i) ) ;
    }
    printf("\nNumber of fitted parameters : %i \n", static_cast<int>(fitted_parameters.size()));

    printf("\n   DMX ranges : %d \n", DMXranges.size());
    for (i = 0 ; i < DMX.size() ; i++)
    {
         printf("DMX %i %.0Lf %.0Lf\n", i+1, DMXranges[i], DMXranges[i+1]);
    }

    printf("\n*** Derived parameters ***\n");

    printf("Mp %.19Le \n", Mp);
    printf("Mi %.19Le \n", Mi);
    printf("ei %.19Le \n", ei);
    printf("eo %.19Le \n", eo);
    printf("omp %.19Le \n", omp);
    printf("omB %.19Le \n", omB);
    printf("tperii %.19Le \n", tperii);
    printf("tperio %.19Le \n", tperio);
    printf("ap %.19Le \n", ap);
    printf("aB %.19Le \n", aB);
    printf("anglii %.19Le \n", anglii);
    printf("anglio %.19Le \n", anglio);
    printf("DopplerF %.19Le \n", DopplerF);
    for (i=0 ; i < quadrupole_kgm2.size(); i++) printf("quadrupole%i (kg.m^2) %.19Le \n", i, quadrupole_kgm2[i]);
    
    printf("\n*** Derived parameters for %d extra bodies ***\n", nextra);
    for (i = 0; i < nextra ; i++)
    {
        printf("a_extra%d %.19Le \n", i+1, a_extra[i]);
        printf("angli_extra%d %.19Le \n", i+1, angli_extra[i]);
        printf("e_extra%d %.19Le \n", i+1, e_extra[i]);
        printf("om_extra%d %.19Le \n", i+1, om_extra[i]);
        printf("tperi_extra%d %.19Le \n", i+1, tperi_extra[i]);
        printf("M_extra%d %.19Le \n", i+1, M_extra[i]);
        printf("\n");
    }

    return;

}



void Rotate_PSB_to_SSB(valarray<value_type> & vector3, const value_type & alpha, const value_type & delta)
// Rotate the vectors expressed in the PSB frame (as expressed from orbel2statevect) in the SSB frame. No Lorentz boost just rotation.
{
    value_type sa = sin(alpha);
    value_type ca = cos(alpha);
    value_type sd = sin(delta);
    value_type cd = cos(delta);
    value_type vi[3] ;

    vi[0] = vector3[0] ;
    vi[1] = vector3[1] ;
    vi[2] = vector3[2] ;


    vector3[0] = vi[0] * (-sa)    +  vi[1] * (  -ca * sd ) + vi[2] *  ca*cd    ;
    vector3[1] = vi[0] * ca       +  vi[1] * (  -sa * sd ) + vi[2] *  sa*cd    ;
    vector3[2] =                     vi[1] * cd            + vi[2] *  sd       ;

    return;
}

void Rotate_PSB_to_SSB(value_type * vector3, const value_type & alpha, const value_type & delta)
// Rotate the vectors expressed in the PSB frame (as expressed from orbel2statevect) in the SSB frame. No Lorentz boost just rotation.
{
    value_type sa = sin(alpha);
    value_type ca = cos(alpha);
    value_type sd = sin(delta);
    value_type cd = cos(delta);
    value_type vi[3] ;

    vi[0] = vector3[0] ;
    vi[1] = vector3[1] ;
    vi[2] = vector3[2] ;

    vector3[0] = vi[0] * (-sa)    +  vi[1] * (  -ca * sd ) + vi[2] *  ca*cd    ;
    vector3[1] = vi[0] * ca       +  vi[1] * (  -sa * sd ) + vi[2] *  sa*cd    ;
    vector3[2] =                     vi[1] * cd            + vi[2] *  sd       ;

    return;
}

void Rotate_SSB_to_PSB(valarray<value_type> & vector3, const value_type & alpha, const value_type & delta)
// Inverse of Rotate_PSB_to_SSB
{
    value_type sa = sin(alpha);
    value_type ca = cos(alpha);
    value_type sd = sin(delta);
    value_type cd = cos(delta);
    value_type vi[3] ;

    vi[0] = vector3[0] ;
    vi[1] = vector3[1] ;
    vi[2] = vector3[2] ;

    vector3[0] = vi[0] * (-sa)          +  vi[1] * ca     ;
    vector3[1] = vi[0] * (-ca*sd)       +  vi[1] * (  -sa * sd ) + vi[2] *  cd    ;
    vector3[2] = vi[0] * ca*cd          +  vi[1] * sa*cd         + vi[2] *  sd       ;

    return;
}

void Rotate_SSB_to_PSB(value_type * vector3, const value_type & alpha, const value_type & delta)
// Inverse of Rotate_PSB_to_SSB
{
    value_type sa = sin(alpha);
    value_type ca = cos(alpha);
    value_type sd = sin(delta);
    value_type cd = cos(delta);
    value_type vi[3] ;

    vi[0] = vector3[0] ;
    vi[1] = vector3[1] ;
    vi[2] = vector3[2] ;

    vector3[0] = vi[0] * (-sa)          +  vi[1] * ca     ;
    vector3[1] = vi[0] * (-ca*sd)       +  vi[1] * (  -sa * sd ) + vi[2] *  cd    ;
    vector3[2] = vi[0] * ca*cd          +  vi[1] * sa*cd         + vi[2] *  sd       ;

    return;
}





void Parametres::Compute_state_vectors(valarray<value_type> & rp, valarray<value_type> & rpt,
                                       valarray<value_type> & ri, valarray<value_type> & rit,
                                       valarray<value_type> & ro, valarray<value_type> & rot,
                                       valarray<valarray<value_type>> & r_extra, valarray<valarray<value_type>> & v_extra)
/*
 * This routine computes all the inner variables that depend on the parameters used .
 * It should be called every time parameters are changed . (it just need all the outer info : par file and data file)
 * It should ensure compatibility with different system of parameters (since the inner variables should be fairly independant on it)
 * Variables set :
 * spinfreq, spinfreq1, delaymax
 * dt_interp, ninterp, toa_in_interp, treference_in_interp
 *
 *  Compute end Initialize all the necessary stuff for integration of the equation of motion       * with the inherited properties from the "Integrateur" class;
 *
 *  Add parameter here: if your parameter is used to initialize some seocndary quantities.
 *
 * TODO : Ideally the part that initialises derived quatities should be separate from the actual computation of state vectors.
 */
{
        int i,j;
        
        value_type Mpi = zero ;
        valarray<value_type> rBi(3);
        valarray<value_type> rBit(3);
        valarray<value_type> nip(3), nio(3), npo(3);
        valarray<value_type> Ppi(3), MXip(3);
        
        value_type rip, rpo, rio = zero ;
        value_type nrpt2, nrit2, nrot2 = zero ;
        value_type M1PN = zero ;
        value_type fmasse = zero;
        value_type Moold = zero;
        value_type Mb = zero;
        value_type nu = zero;
        value_type norbi = zero;
        value_type sepi = zero;

        valarray<valarray<value_type>>  rcdmextra(nextra);
        valarray<valarray<value_type>>  vcdmextra(nextra);
        for (i =0; i < nextra ; i++)
        {
            rcdmextra = valarray<value_type>(3);
            vcdmextra = valarray<value_type>(3);
        }
        
        
        DopplerF = un;

    // SEP parameters, assuming index 0 = neutron star and index>0 weak field body (eg white dwarf):

        for (i= 1 ; i < nbodies_plus_extra ; i++)
        {
            Gg[0][i] = un + SEP_D;
            Gg[i][0] = un + SEP_D;
            gammabar[0][i] = SEP_gamma;
            gammabar[i][0] = SEP_gamma;
            betabar[i][0][0] = SEP_beta_0pp;
            for (j= 1 ; j < nbodies_plus_extra ; j++) 
            {
                betabar[0][i][j] = SEP_beta_p00;
                betabar[i][0][j] = SEP_beta_0p0;
            }
        }

      // Orbital parameters
        // Outer orbit
             if (Mo > zero) {// Three-body case
                 eo = sqrt( pow( etaB, 2) + pow(kappaB, 2) ) ;
                 omB = inversetrigo( kappaB / eo, etaB / eo ) ;
//                  tperio = ( tascB + Po * omB / deuxpi ) * DopplerF;
                 aB = sqrt( pow( aBsinio, 2 ) + pow(aBcosio, 2 ) ) ;
                 anglio = acos( aBcosio / aB ) * sgn( aBsinio ) ;
             }
        // Inner orbit
             ei = sqrt( pow(etap,2) + pow(kappap,2) ) ;
             omp = inversetrigo( kappap / ei , etap / ei ) ;
             // tperii moved below in case Pi determined from masses
             if (parameter_set==3 or parameter_set==4 or parameter_set==5 or parameter_set == 6)
             {
               anglii = anglio + delta_i * raddeg;
               ap =  apsinii / sin(anglii);
               apcosii = ap * cos(anglii);
             }
             else
             {
               ap = sqrt( pow(apsinii, 2) + pow( apcosii , 2 ) ) ;
               anglii = acos( apcosii / ap ) * sgn( apsinii ) ;
             }
        // Checks on eccentricities
             if (ei > un )
             {
                 cout << endl << "Warning ! ei > 1" << endl;
                 printf("   etap kappap ei omp %.19Le %.19Le %.19Le %.19Le \n ", etap, kappap, ei, omp);
                 cout << "setting ei = 0.99999999999" << endl << endl;
                 ei = 0.99999999999;
             }
             if (eo > un )
             {
                 cout << endl << "Warning ! eo > 1" << endl;
                 printf("   etaB kappaB eo omB %.19Le %.19Le %.19Le %.19Le \n ", etaB, kappaB, eo, omB);
                 cout << "setting eo = 0.99999999999" << endl << endl;
                 eo = 0.99999999999;
             }




        // Deduce masses from parameters in the same way as in parameter system 1 and 10
            if (Mo > zero) {// Three-body case
                if (parameter_set == 0)
                {
                // Uncomment the following lines to use kepler's third law to deduce the masses
                    Mi = muio ;
                    fmasse = quatre*pi*pi *pow(ap,3) / (Pi*Pi * daysec * daysec * GMsol);
                    Mp = sqrt( pow(Mi,3) / fmasse ) - Mi ;
                    Mb = Mp + Mi ;
                    fmasse = quatre*pi*pi *pow(aB,3) / (Po*Po * daysec * daysec * GMsol);
                    // Solve kepler's third law for third body
                    Mo = 0.4; // initial guess
                    do
                    {
                        Moold = Mo;
                        Mo = pow(fmasse * (Mo + Mb) * (Mo+Mb) , un/trois );
                    }
                    while (abs(Mo - Moold) > pow(dix, -15));
                } else if (parameter_set == 1)
                {
                    Mp = muio + mui; // mui = (Mp + Mi) * 0.5  et muio = (Mp - Mi) * 0.5
                    Mi = mui - muio;
                    // uncomment to deduce periods from masses
                    Pi = sqrt(quatre*pi*pi *pow(ap,3) / ( daysec * daysec * GMsol) / pow(Mi,3)) * (Mp + Mi);
                    Mb = Mp + Mi ;
                    Po = sqrt(quatre*pi*pi *pow(aB,3) / ( daysec * daysec * GMsol) / pow(Mo,3)) * (Mb + Mo);
                }
               else if (parameter_set == 2 or parameter_set == 3 or parameter_set==4 or parameter_set==5  or parameter_set == 6) // mp, mi, mo given by : mui = Mi/Mp, Pi, Po
                {
                  sepi = ap * (1. + 1./mui); // separation of the inner binary
                  norbi = deuxpi/(Pi*daysec);
                  Mb = pow(sepi,3) * pow(norbi,2) / (Ggrav*Msol);
                  Mp = Mb / (1. + mui);
                  Mi = mui * Mp;

                  fmasse = quatre*pi*pi *pow(aB,3) / (Po*Po * daysec * daysec * GMsol);
                  // Solve kepler's third law for third body
                  Mo = 0.4; // initial guess
                  do
                  {
                      Moold = Mo;
                      Mo = pow(fmasse * (Mo + Mb) * (Mo+Mb) , un/trois );
                  }
                  while (abs(Mo - Moold) > pow(dix, -15));
                }
            // version with masses fitted
//                 Mp = undemi * ( muio/ Mo + sqrt( pow( muio/Mo , 2 ) - quatre * mui ) ) ;
//                 Mi = undemi * ( muio/ Mo - sqrt( pow( muio/Mo , 2 ) - quatre * mui ) ) ;
             // version with masses fitted  and a different parametrization

           }
            else { // 2-body case
                Mp = mui ;
                Mi = muio ;
                if (integrator_type == 1 ) // 1PN case
                {
                    fmasse = quatre*pi*pi *pow(ap,3) / (Pi*Pi * daysec * daysec * GMsol);
                    Mp = 1.4;
                    do
                    {
                        nu = Mp * Mi / pow(Mp + Mi, 2);
                        Moold = Mp ;
                        Mp = sqrt( pow(Mi,3) / fmasse )* (un  + GMsol * Mi / (deux*ap * clight2) * ( -neuf + nu) )  - Mi ;
                    } while (abs( Moold - Mp) > pow(dix, -15));
                }
                else //Newtonian case
                {
                    if ( parameter_set ==1 ) 
                    {
                        // Deduce pulsar mass from companion's mass
                        fmasse = quatre*pi*pi *pow(ap,3) / (Pi*Pi * daysec * daysec * GMsol);
                        Mp = sqrt( pow(Mi,3) / fmasse ) - Mi ;
                    }
                    else
                    {
                        fmasse = quatre*pi*pi *pow(ap,3) / (Pi*Pi * daysec * daysec * GMsol);
                        // Solve kepler's equation for companion
                        Mi = 0.4; // initial guess
                        do
                        {
                            Moold = Mi;
                            Mi = pow(fmasse * (Mi + Mp) * (Mi+Mp) , un/trois );
                        }
                        while (abs(Mi - Moold) > pow(dix, -15));
                    }
                }

            }

            value_type vp2, vb2;
            value_type Uo, Ui;
            value_type DT;
            value_type muperp, mupara;

            if (parameter_set == 4 or parameter_set==5) // pset 4 is equivalent to pset3 with Einstein delay correction to spin frequency
            {
              vp2 = deuxpi * ap / (Pi *daysec) / clight;
              vp2 *= vp2;
              vb2 = deuxpi * aB / (Po*daysec)  / clight;
              vb2 *= vb2;
              Uo = GMsol * Mo / (clight2 * aB *(1 + Mb / Mo));
              Ui = GMsol * Mi / (clight2 * ap *(1 + Mp / Mi));
              true_spinfreq = f/(un - 0.5*(vp2 + vb2) - Ui - Uo);
            }
            else
            {
              true_spinfreq = f;
            }

            if (parameter_set == 5 ) // pset 5 is equivalent to pset 4 with Schklovskii correction to spin frequency derivative
            {
              value_type muperp = sqrt(RA1 * RA1 + DEC1 * DEC1) * radmasdeg / yrsec;
              value_type mupara = distance1 * radmasdeg / yrsec;
              true_spinfreq1 = f1 + f* distance * lightyrmeter * muperp * muperp *(1 - muperp * mupara)/clight * daysec; //  there's a mistake here. See pset 6 below for the correct expression
            }
            else if (parameter_set == 6) // pset 6  is equivalent to pset 5 but with a corrrected expression for the spin freq derivative and shklovskii correction to the spin frequency too
            {
              vp2 = deuxpi * ap / (Pi *daysec) / clight;
              vp2 *= vp2;
              vb2 = deuxpi * aB / (Po*daysec)  / clight;
              vb2 *= vb2;
              Uo = GMsol * Mo / (clight2 * aB *(1 + Mb / Mo));
              Ui = GMsol * Mi / (clight2 * ap *(1 + Mp / Mi));
              muperp = sqrt(RA1 * RA1 + DEC1 * DEC1) * radmasdeg / yrsec;
              mupara = distance1 * radmasdeg / yrsec;
              DT = (treference - posepoch) * daysec;

              true_spinfreq = f/(un - 0.5*(vp2 + vb2) - Ui - Uo - distance * lightyrmeter * muperp * muperp * DT * (1 - 1.5*DT * mupara)/clight);

              true_spinfreq1 = f1 +  true_spinfreq * distance * lightyrmeter * muperp * muperp *(1 - 3*DT * mupara)/clight * daysec; // day^-2
            }
            else
            {
              true_spinfreq1 = f1;
            }

            // Derive internal quadrupole parameters
            quadrupole_kgm2.assign(quadrupole.begin(), quadrupole.end());
            if (parameter_set == 3 or parameter_set == 4 or parameter_set==5 or parameter_set==6) // quadrupole[0] = J2 = 3/2 Q/(M1 a^2) where a is the separation of the inner binary
            {
              for (i = 0 ; i < quadrupole.size(); i++)
              {
                quadrupole_kgm2[i] *=  deuxtiers * Mi*Msol * pow(ap *(1+ Mp/Mi),2); // in kg.m^2
              }
              if (quadrupole.size() > 0  and integrator_type != 2 and msg_count < 10)
              {
                msg_count++;
                printf("\n WARNING : Bad combination of quadrupole and integrator type !\n");
              }
            }
            else
            {
              printf("\n WARNING : with the current parameter set quadrupole_kgm2 might be in a different unit than kg.m^2. Implementation to do ?\n");
            }


            // Determine tperi here in case Pi and Po are determined from masses (above)
            tperii = (tascp   + Pi * omp / deuxpi) * DopplerF ; // Temps de passage au périastre par rapport au temps de passge au noeud ascendant
            if (Mo > zero)  {
                tperio = ( tascB + Po * omB / deuxpi ) * DopplerF;
            }
            Mpi = Mp + Mi ;

            value_type sv[6], sv1[6], sv2[6] ; // temporary state vector

            // Initial orbital elements to state vector for inner binary
            if (integrator_type == 1 or integrator_type==2 or integrator_type == 3)
            {
                 orbel2statevect_1pn(treference, ei, DopplerF*ap, omp, anglii, tperii, DopplerF*Pi, oman +delta_oman, Mp, Mi,
                            sv,
                            1, 1) ;}
            else
                orbel2statevect(treference, ei, DopplerF*ap, omp, anglii, tperii, DopplerF*Pi, oman +delta_oman,
                            sv,
                            1, 1) ;


            rp[0] = sv[0] ; rp[1] = sv[1] ; rp[2] = sv[2] ;
            rpt[0] = sv[3] ; rpt[1] = sv[4] ; rpt[2] = sv[5] ;

            // Initial orbital elements to state vector for outer binary
            if (Mo > zero) {
                if (integrator_type == 1 or integrator_type==2  or integrator_type == 3)
                    orbel2statevect_1pn(treference , eo,DopplerF* aB, omB, anglio, tperio,DopplerF* Po, oman, Mpi, Mo,
                                sv,
                                1, 1) ;
                 else
                    orbel2statevect(treference , eo,DopplerF* aB, omB, anglio, tperio,DopplerF* Po, oman,
                                sv,
                                1, 1) ;



                rBi[0] = sv[0] ; rBi[1] = sv[1] ; rBi[2] = sv[2] ;
                rBit[0] = sv[3] ; rBit[1] = sv[4] ; rBit[2] = sv[5] ;
                
            // Computing derived quantities and masses for extra bodies 
            value_type Mcdm = Mp + Mi + Mo;
            for (i=0 ; i < nextra ; i++)
            {
                e_extra[i] = sqrt( pow(eta_extra[i],2) + pow(kappa_extra[i],2) ) ;
                if (e_extra[i] > 0)
                    om_extra[i] = atan2(eta_extra[i]/e_extra[i], kappa_extra[i]/e_extra[i]);
                else
                    om_extra[i] = 0.;
                a_extra[i] = sqrt( pow(asini_extra[i],2) + pow(acosi_extra[i],2) ) ;
                if (a_extra[i] > 0.)
                    angli_extra[i] = atan2(asini_extra[i]/a_extra[i], acosi_extra[i]/a_extra[i]);
                else 
                    angli_extra[i] = 0.;
                tperi_extra[i] = tasc_extra[i] + om_extra[i]/deuxpi * P_extra[i] ;
                if (P_extra[i] == 0. or a_extra[i] == 0)
                    M_extra[i] == 0.;
                else 
                {
                    fmasse = quatre*pi*pi *pow(a_extra[i],3) / (P_extra[i]* P_extra[i] * daysec * daysec * GMsol);
                    // Solve kepler's third law for third body
                    M_extra[i] = Mcdm * pow(fmasse / Mcdm, un/trois); // initial guess, accurate if M_extra[i] << Mcdm
                    do
                    {
                      Moold = M_extra[i];
                      M_extra[i] = pow(fmasse * (M_extra[i] + Mcdm) * (M_extra[i] + Mcdm) , un/trois );
                    }
                    while (abs(M_extra[i] - Moold) > pow(dix, -15));
                    Mcdm += M_extra[i];
//                     M_extra[i] = sqrt(GMsol * pow(M3/a_extra[i],3) * pow(P_extra[i]*daysec/deuxpi,2)) - M3; 
                }
            }
                    
            // Computing initial state vectors of extra bodies
            for (i = 0; i < nextra; i++)
            {
                if (M_extra[i] > 0) {
                    orbel2statevect(treference , e_extra[i],DopplerF* a_extra[i], om_extra[i], angli_extra[i], tperi_extra[i],DopplerF* P_extra[i], oman_extra[i],
                                    sv,
                                    1, 1) ;
                    for (j = 0; j <3 ; j++) 
                    {
                        rcdmextra[i][j] =  sv[j];
                        vcdmextra[i][j] =  sv[j+3];
                    }
                }
                else 
                {
                    r_extra[i] = 0.; 
                    v_extra[i] = 0.;
                }
            }

                // Changement d'origine de rp vers le centre de masse du système interne
                if (integrator_type == 1  or integrator_type==2 or integrator_type==3) {    // Respect the 1PN conserved center of mass and momentum set to zero
                    ri = - Mp / Mi * rp ;
                    rit = - Mp / Mi * rpt ;
                    rip = norm3d<value_type>(ri - rp) ;
                    nrit2 = pow( rit[0], 2) + pow( rit[1], 2) + pow( rit[2], 2) ;
                    nrpt2 = pow( rpt[0], 2) + pow( rpt[1], 2) + pow( rpt[2], 2) ;
                    nip = (rp - ri ) / rip ;

                    rit = rit - un/(deux * pow( clight, 2 ) * Mi) * ( Mp * rpt * nrpt2 + Mi * rit * nrit2 -
                                                    GMsol * Gg[0][1] * Mi * Mp / rip * ( rit + rpt + dotprod3d<value_type>(rit + rpt, nip) * nip) ) ;
                    // Go to the observer's frame
                    rpt = rpt + rBit - un/(deux* pow( clight, 2) ) * ( deux * dotprod3d<value_type>(rpt, rBit) - GMsol * Gg[0][1] * Mi / rip * dotprod3d<value_type>( rBit, nip ) * nip ) ;
                    rit = rit + rBit - un/(deux* pow( clight, 2) ) * ( deux * dotprod3d<value_type>(rit, rBit) - GMsol * Gg[0][1] * Mp / rip * dotprod3d<value_type>( rBit, nip ) * nip ) ;

                    ri = ri + rBi ;
                    rp = rp + rBi ;
                    
//                                 printf(" \nb rp0 %.10Le\n", rp[0]); // ! Test !

                                
                    // Adjust the state vector of the third body to respect the nullity of the total momentum and center of mass
                    ro = -rBi * Mpi/Mo ;
                    rot = -rBit * Mpi/Mo ;
                    rpo = norm3d<value_type>(ro - rp) ;
                    rio = norm3d<value_type>(ro - ri) ;
                    nrit2 = pow( rit[0], 2) + pow( rit[1], 2) + pow( rit[2], 2) ;
                    nrpt2 = pow( rpt[0], 2) + pow( rpt[1], 2) + pow( rpt[2], 2) ;
                    nrot2 = pow( rot[0], 2) + pow( rot[1], 2) + pow( rot[2], 2) ;
                    nip = (rp - ri ) / rip ;
                    nio = (ro - ri ) / rio ;
                    npo = (ro - rp ) / rpo ;

                    M1PN = ( Mpi ) * ( un +
                                            ( Mp * nrpt2 + Mi * nrit2 - GMsol * Gg[0][1] * Mi * Mp / rip ) / (deux * Mpi * pow(clight, 2 ) ) ) ;

                    ri = un / Mi * (rBi * M1PN - Mp * rp - un / (deux * clight2) * ( Mp * nrpt2 * rp + Mi * nrit2 * ri - GMsol* Gg[0][1] * Mp * Mi / rip * ( rp + ri ) ) ) ;

                    // Part of the impulsion depending only on the inner binary
                    Ppi = rpt * Mp + rit * Mi + un/(deux* clight2 ) * ( Mp * rpt * nrpt2 + Mi * rit * nrit2 -
                                                    GMsol * Gg[0][1] * Mi * Mp / rip * ( rit + rpt + dotprod3d<value_type>(rit + rpt, nip) * nip) ) ;

                    rot = - un/Mo * ( Ppi + un/(deux * clight2 ) * ( Mo * nrot2 * rot -
                                                            Mp * GMsol * Gg[0][2] * Mo / rpo * ( rpt + rot ) - Mi * Gg[1][2] * GMsol * Mo / rio * ( rit + rot ) -
                                                            GMsol * Gg[0][2] * Mp * Mo / rpo * dotprod3d<value_type>(rpt + rot, npo) * npo -
                                                            GMsol * Gg[1][2] * Mi * Mo / rio * dotprod3d<value_type>(rit + rot, nio) * nio ) ) ;

                    // Part of the " Mass * center of mass " depending only on the inner binary
                    MXip = Mp * rp + Mi * ri + un / (deux * clight2) * (Mp * nrpt2 * rp + Mi * nrit2 * rit - GMsol * Gg[0][1] * Mi * Mp / rip * ( ri + rp ) ) ;

                    ro = - un / Mo * ( MXip + un/ (deux * clight2) * ( Mo * nrot2 * ro  - GMsol * Gg[0][2] * Mp * Mo / rpo * ( rp + ro ) - GMsol * Gg[1][2] * Mi * Mo / rio * ( ri + ro ) )) ;

            // Accounting for extra bodies ----------------------
                    if (nextra > 0)
                    {
                        Mcdm = Mp + Mi + Mo; // Mass of virtual centre of mass (cdm) object
                        
                        for (i = 0 ; i < nextra ; i++) 
                        {
                            rp += rcdmextra[i] ;
                            rpt += vcdmextra[i] ;
                            ri += rcdmextra[i] ;
                            rit += vcdmextra[i] ;
                            ro += rcdmextra[i] ;
                            rot += vcdmextra[i] ;
                            r_extra[i] = - rcdmextra[i] * Mcdm / M_extra[i];
                            v_extra[i] = - vcdmextra[i] * Mcdm / M_extra[i];
                            Mcdm += M_extra[i];
                        }
                        
                        // Fine-tune the global centre of mass to 0 at 1PN
                        valarray<value_type> center_of_mass(3);
                        valarray<value_type> impulsion(3);
                        value_type energy;
                        int nbody = 3+nextra;
                        valarray<value_type> svfull(nbody*6);
                        valarray<value_type> Ms(nbody);
                        Ms[0] = Mp;
                        Ms[1] = Mi;
                        Ms[2] = Mo;
                        for (j = 0; j < nextra ; j++) Ms[3 + j] = M_extra[j];
                        
                        int count = 0;
                        
                        do // iterate until the centre of mass is 0 enough
                        {
                            for (i = 0 ; i < 3; i++)  // just for conversion to a valarray 
                            {
                                svfull[i] = rp[i];
                                svfull[i + 3*nbody] = rpt[i];
                                
                                svfull[i + 3] = ri[i];
                                svfull[i + 3 + 3*nbody] = rit[i];
                                
                                svfull[i + 6] = ro[i];
                                svfull[i + 6 + 3*nbody] = rot[i];
                                
                                for (j = 0; j < nextra ; j++)
                                {
                                    svfull[i + 9 + 3 * j] = r_extra[j][i];
                                    svfull[i + 9 + 3 * j + 3*nbody] = v_extra[j][i];
                                }
                            }
                            IntegralePrems3_1PN(svfull, Ms,
                            0.,
                            Gg, gammabar, betabar,
                            center_of_mass, impulsion,
                            energy  );
                            
                            impulsion /= Mcdm; // approximate center of mass velocity
                            
                            rp -= center_of_mass;
                            rpt -= impulsion ;
                            ri -= center_of_mass ;
                            rit -= impulsion ;
                            ro -= center_of_mass ;
                            rot -= impulsion ;
                            for (j=0; j < nextra ; j++) 
                                {
                                    r_extra[j] -= center_of_mass;
                                    v_extra[j] -= impulsion;
                                }
                            count +=1;
                        }
                        while ((norm3d(center_of_mass) > 1.e-1 or norm3d(impulsion) > 1e-14) and count < 10);
/*                     
 *
 *   valarray<value_type> Re(3, zero), Ve(3,zero);
                        valarray<value_type> R3(3, zero), P3(3,zero);
                        value_type E3 = 0.;
                        value_type Me = 0.;
                        
                        for (i = 0 ; i < nextra ; i++) // Compute centre of mass position and velocity of extra bodies
                        {
                            Re += M_extra[i] * r_extra[i]  ;
                            Ve += M_extra[i] * v_extra[i]  ;
                            Me += M_extra[i];
                        }
                        if (Me > 0.) 
                        {
                            Re /= Me ;
                            Ve /= Me ;
                            Compute centre of mass position and velocity of triple system
                            for (j = 0; j < 3 ; j++) {sv[j] = rp[j]; sv[j+3] = rpt[j];};
                            for (j = 0; j < 3 ; j++) {sv1[j] = ri[j]; sv1[j+3] = rit[j];};
                            for (j = 0; j < 3 ; j++) {sv2[j] = ro[j]; sv2[j+3] = rot[j];};
                            IntegralePrems3_1PN(sv, sv1, sv2, Mp, Mi, Mo,  0., Gg, gammabar, betabar, R3, P3, E3 );
                            Update position and velocity of triple system so that the total centre of mass is zero at Newtonian order
                            rp += - Re * Me / M3 - R3/M3;
                            rpt += - Ve * Me / M3 - P3/M3;
                            ri += - Re * Me / M3 - R3/M3;
                            rit += - Ve * Me / M3 - P3/M3;
                            ro += - Re * Me / M3 - R3/M3;
                            rot += - Ve * Me / M3 - P3/M3;
                            Recompute centre of mass position and velocity of triple system
                            for (j = 0; j < 3 ; j++) {sv[j] = rp[j]; sv[j+3] = rpt[j];};
                            for (j = 0; j < 3 ; j++) {sv1[j] = ri[j]; sv1[j+3] = rit[j];};
                            for (j = 0; j < 3 ; j++) {sv2[j] = ro[j]; sv2[j+3] = rot[j];};
                            IntegralePrems3_1PN(sv, sv1, sv2, Mp, Mi, Mo,  0., Gg, gammabar, betabar, R3, P3, E3 );
                            Update extra bodies position and velocities to compensate for the shift introduced by 1PN terms in the last calculation
                            for (i = 0 ; i < nextra ; i++)
                            {
                                r_extra[i] -= Re + R3 * M3 / Me ;
                                v_extra[i] -= Ve + P3 / Me ;
                            }
                        }
                        else 
                        {
                            printf("\n Warning : total mass of extra bodies <= 0 (Me = %.19Le)\n\n", Me);
                        }*/
                    }
                }
                else {
                    rp = rp + rBi ;
                    rpt = rpt + rBit ;
                    ro = -rBi * Mpi/Mo ;
                    rot = -rBit * Mpi/Mo ;
                    ri = -(rp - rBi) * Mp / Mi + rBi ;
                    rit = -(rpt - rBit) * Mp / Mi + rBit ;
                    if (nextra > 0)
                    {
                        Mcdm = Mp + Mi + Mo; // Mass of virtual centre of mass (cdm) object
                        
                        for (i = 0 ; i < nextra ; i++) 
                        {
                            rp += rcdmextra[i] ;
                            rpt += vcdmextra[i] ;
                            ri += rcdmextra[i] ;
                            rit += vcdmextra[i] ;
                            ro += rcdmextra[i] ;
                            rot += vcdmextra[i] ;
                            r_extra[i] = - rcdmextra[i] * Mcdm / M_extra[i];
                            v_extra[i] = - vcdmextra[i] * Mcdm / M_extra[i];
                            Mcdm += M_extra[i];
                        }
                    }
                }
            }
            else { // 2-body case
                if (integrator_type == 1 ) {
                    // In the inner binary frame
                    rBi = zero;
                    rBit = zero;
                    ro = zero;
                    rot=zero;

                    ri = - Mp / Mi * rp ;
                    rit = - Mp / Mi * rpt ;
                    rip = norm3d<value_type>(ri - rp) ;
                    nrit2 = pow( rit[0], 2) + pow( rit[1], 2) + pow( rit[2], 2) ;
                    nrpt2 = pow( rpt[0], 2) + pow( rpt[1], 2) + pow( rpt[2], 2) ; //np.sum(rpt**2)
                    nip = (rp - ri ) / rip ;



                    rit = rit - un/(deux * pow( clight, 2 ) * Mi) * ( Mp * rpt * nrpt2 + Mi * rit * nrit2 -
                                                    GMsol * Gg[0][1] * Mi * Mp / rip * ( rit + rpt + dotprod3d<value_type>(rit + rpt, nip) * nip) ) ;

                    // Adjust the state vector of the third body to respect the nullity of the total momentum and center of mass
                    nrit2 = pow( rit[0], 2) + pow( rit[1], 2) + pow( rit[2], 2) ;
                    nrpt2 = pow( rpt[0], 2) + pow( rpt[1], 2) + pow( rpt[2], 2) ;
                    nrot2 = pow( rot[0], 2) + pow( rot[1], 2) + pow( rot[2], 2) ;
                    nip = (rp - ri ) / rip ;


                    M1PN = ( Mpi ) * ( un +
                                            ( Mp * nrpt2 + Mi * nrit2 - GMsol * Gg[0][1] * Mi * Mp / rip ) / (deux * Mpi * pow(clight, 2 ) ) ) ;

                    ri = un / Mi * ( - Mp * rp - un / (deux * clight2) * ( Mp * nrpt2 * rp + Mi * nrit2 * ri - GMsol * Gg[0][1] * Mp * Mi / rip * ( rp + ri ) ) ) ;



////         Uncomment the following lines to check that momentum and center of mass are zero
// //                     cout << "\n Test pos and momentum \n" ; // ! Test !
// //                     Print_table(MXip/ Mpi);
// //                     Print_table(Ppi / Mpi);
                } else {
                    rBi = zero ;
                    rBit = zero ;
                    // Changement d'origine de rp vers le centre de masse du système interne
                    rp = rp + rBi ;
                    rpt = rpt + rBit ;
                    ro = rBi ;
                    rot = rBit ;
                    ri = -(rp - rBi) * Mp / Mi + rBi    ;
                    rit = -(rpt - rBit) * Mp / Mi + rBit ;
                }
            }

    // Rotation to put everything in the same frame as the solar system in tempo i.e. with the pulsar oriented in direction gven by (RA,DEC)
    if (geometric == true )
    {
        Rotate_PSB_to_SSB(rp, RA, DEC) ;
        Rotate_PSB_to_SSB(rpt, RA, DEC) ;
        Rotate_PSB_to_SSB(ri, RA, DEC) ;
        Rotate_PSB_to_SSB(rit, RA, DEC) ;
        Rotate_PSB_to_SSB(ro, RA, DEC) ;
        Rotate_PSB_to_SSB(rot, RA, DEC) ;
        for (i = 0 ; i < nextra ; i++)
        {
            Rotate_PSB_to_SSB(r_extra[i], RA, DEC) ;
            Rotate_PSB_to_SSB(v_extra[i], RA, DEC) ;
        }
    }

  // Set some inner variables
            if (Mo <= zero ) aB = zero;


    return ;
}





void Parametres::Set_parameter_relativeshift(int paramnumber, const double newvalue)
{
    Set_parameters(paramnumber, parameters_ini[paramnumber] + Parameter_scale(paramnumber) * static_cast<value_type>( newvalue ) );
    return;
}


void Parametres::Set_fitted_parameter_relativeshifts(vector<double> parameter_shifts)
{   // parameter_shifts contains as many parameters as fitted_parameters.size()
    int p =0;
    value_type param = 0.;
    motion_changed = false;
    for (int i = 0 ; i < fitted_parameters.size() ;  i++)
    {
        p = fitted_parameters[i];
        // Add parameter here : in the below if statement, all the parameters which change implies recomputing the orbits must be in the "if", otherwise in the "else".
        if  ( (p >= absolute_parameter_map[2] && p <= absolute_parameter_map[18]) || p >= absolute_parameter_map[20] && p <= absolute_parameter_map[24] ||  p == absolute_parameter_map[37] || (p >= absolute_parameter_map[38] && p <= absolute_parameter_map[44]))  { // i == absolute_parameter_map[35]  ||  removed on this branch ! 
            param = Parameters(p);
            Set_parameters(p, parameters_ini[p] + Parameter_scale(p) * static_cast<value_type>(parameter_shifts[i]) );
            if (Parameters(p) != param ) 
            {
                motion_changed = true ;
            }
            
        }
        else
        {
            Set_parameters(p, parameters_ini[p] + Parameter_scale(p) * static_cast<value_type>(parameter_shifts[i]) );
        }
    }
    return;
}

void Parametres::Set_fitted_parameter_relativeshifts(double * parameter_shifts)
{
  vector<double> paramshifts(parameter_shifts, parameter_shifts + fitted_parameters.size());
  Set_fitted_parameter_relativeshifts(paramshifts);
}


void Parametres::Set_parameters_relativeshift(double * parameter_shifts)
{
    value_type param = 0.;
    motion_changed = false;
    for (int i = 0 ; i < nfitparams ; i++)
    {
        // Add parameter here : in the below if statement, all the parameters which change implies recomputing the orbits must be in the "if", otherwise in the "else".
        if  (  (i >= absolute_parameter_map[2] && i <= absolute_parameter_map[18]) || i >= absolute_parameter_map[20] && i <= absolute_parameter_map[24] ||  i == absolute_parameter_map[37] || (i >= absolute_parameter_map[38] && i <= absolute_parameter_map[44]))  {// i == absolute_parameter_map[35]  ||  removed on this branch ! 
            param = Parameters(i);
            Set_parameter_relativeshift(i, parameter_shifts[i]);
            if (Parameters(i) != param ) motion_changed = true ;
        }
        else {
            Set_parameter_relativeshift(i, parameter_shifts[i]);
        }
    }
    return;
}

void Parametres::Set_parameters_relativeshift(vector<double>  parameter_shifts)
{
    value_type param = 0.;
    motion_changed = false;
    for (int i = 0 ; i < nfitparams ; i++)
    {
        // Add parameter here : in the below if statement, all the parameters which change implies recomputing the orbits must be in the "if", otherwise in the "else".
        if  ( (i >= absolute_parameter_map[2] && i <= absolute_parameter_map[18]) || i >= absolute_parameter_map[20] && i <= absolute_parameter_map[24] ||  i == absolute_parameter_map[37] || (i >= absolute_parameter_map[38] && i <= absolute_parameter_map[44]))  { // i == absolute_parameter_map[35]  ||  removed on this branch ! 
            param = Parameters(i);
            Set_parameter_relativeshift(i, parameter_shifts[i]);
            if (Parameters(i) != param ) motion_changed = true ;
        }
        else {
            Set_parameter_relativeshift(i, parameter_shifts[i]);
        };
    }
    return;
}



void Parametres::Convert(int target_pset, int param_number,  value_type & converted_value)
{
  // Run Compute_state_vectors before use !
    int abs_param_number, index;
    if (target_pset == 5 and parameter_set == 4)
    {
      ParamList2AbsParam(param_number, abs_param_number, index);
      if (abs_param_number == 1)
      {
        value_type muperp = sqrt(RA1 * RA1 + DEC1 * DEC1) * radmasdeg / yrsec;
        value_type mupara = distance1 * radmasdeg / yrsec;
        converted_value = f1 -   f* distance * lightyrmeter * muperp * muperp *(1 - muperp * mupara)/clight * daysec; //  there's a mistake here. See pset 6 below for the correct expression
      }
      else
      {
        converted_value = Parameters(param_number);
      }
    }
    else if (target_pset == 6 and parameter_set == 5)
    {
      value_type vp2, vb2;
      value_type Uo, Ui;
      ParamList2AbsParam(param_number, abs_param_number, index);
      if (abs_param_number == 0)
      {
        vp2 = deuxpi * ap / (Pi *daysec) / clight;
        vp2 *= vp2;
        vb2 = deuxpi * aB / (Po*daysec)  / clight;
        vb2 *= vb2;
        Uo = GMsol * Mo / (clight2 * aB *(1 + (Mp + Mi) / Mo));
        Ui = GMsol * Mi / (clight2 * ap *(1 + Mp / Mi));
        value_type muperp = sqrt(RA1 * RA1 + DEC1 * DEC1) * radmasdeg / yrsec; // in rad/sec
        value_type mupara = distance1 * radmasdeg / yrsec;
        value_type DT = (treference - posepoch) * daysec; // in sec
        converted_value = true_spinfreq * ( 1 - ( 0.5*(vp2 + vb2) + Ui + Uo) - distance * lightyrmeter * muperp * muperp * DT * (1 - 1.5*DT * mupara)/clight ); // day^-1
      }
      else if (abs_param_number == 1)
      {
        value_type muperp = sqrt(RA1 * RA1 + DEC1 * DEC1) * radmasdeg / yrsec; // in rad/sec
        value_type mupara = distance1 * radmasdeg / yrsec;
        value_type DT = (treference - posepoch) * daysec; // in sec
        converted_value = true_spinfreq1 -   true_spinfreq * distance * lightyrmeter * muperp * muperp *(1 - 3*DT * mupara)/clight * daysec; // day^-2
      }
      else
      {
        converted_value = Parameters(param_number);
      }
    }
    else
    {
      printf("\n Error : Parameter conversion from %i to %i not implemented !\n", parameter_set,target_pset);
    }
}

void Parametres::Convert_in_place(int target_pset)
{
  // This bit is really not optimal
   valarray<value_type>  rp(3);
   valarray<value_type>  ri(3);
   valarray<value_type>  ro(3);
   valarray<value_type>  rpt(3);
   valarray<value_type>  rit(3);
   valarray<value_type>  rot(3);
   valarray<valarray<value_type>> r_extra(nextra);
   valarray<valarray<value_type>> v_extra(nextra);
   for (int i = 0; i < nextra; i++)
   {
       r_extra[i] = valarray<value_type>(3);
       v_extra[i] = valarray<value_type>(3);
   }

   Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
   // end of non optimal bit
   int abs_param_number, index;
  value_type converted_value;
  if (target_pset == 5 and parameter_set == 4 or target_pset == 6 and parameter_set == 5)
  {
    for (int i = 0 ; i < nfitparams ; i ++)
    {
      Convert(target_pset, i, converted_value);
      ParamList2AbsParam(i, abs_param_number, index);
      Set_parameters(i, converted_value);


    }
    parameter_set = target_pset;
    Set_reference_to_current_parameters();
    Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
  }
  else
  {
    printf("\n Error : Parameter conversion from %i to %i not implemented !\n", parameter_set,target_pset);
  }
}

void Parametres::Convert_relative_shifts(int target_pset, double * parameter_shifts, double * converted_parameter_shifts)
{
  // This bit is really not optimal
   valarray<value_type>  rp(3);
   valarray<value_type>  ri(3);
   valarray<value_type>  ro(3);
   valarray<value_type>  rpt(3);
   valarray<value_type>  rit(3);
   valarray<value_type>  rot(3);
   valarray<valarray<value_type>> r_extra(nextra);
   valarray<valarray<value_type>> v_extra(nextra);
   for (int i = 0; i < nextra; i++)
   {
      r_extra[i] = valarray<value_type>(3);
      v_extra[i] = valarray<value_type>(3);
   }
  
   Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
   // end of non optimal bit

  value_type converted_value;
  if (target_pset == 5 and parameter_set == 4)
  {
    for (int i = 0 ; i < fitted_parameters.size() ; i++) converted_parameter_shifts[i] = parameter_shifts[i];
    int fittedindex = Get_fitted_parameter_index(1);
    int paramnumber = absolute_parameter_map[1];
    value_type refvalue;

    if (fittedindex >= 0)
    {
      Reset_to_initial_parameters();
      Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
      Convert(target_pset, paramnumber,  refvalue);

      Set_fitted_parameter_relativeshifts(parameter_shifts);
      Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
      Convert(target_pset, paramnumber,  converted_value);
      converted_parameter_shifts[fittedindex] = (converted_value - refvalue)/Parameter_scale(paramnumber);
   }
  }

  if (target_pset == 6 and parameter_set == 5)
  {
    for (int i = 0 ; i < fitted_parameters.size() ; i++) converted_parameter_shifts[i] = parameter_shifts[i];
// ------------------------------------ For f
    int fittedindex = Get_fitted_parameter_index(0); // f
    int paramnumber = absolute_parameter_map[0];     // f
    value_type refvalue;

    if (fittedindex >= 0)
    {
      Reset_to_initial_parameters();
      Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
      Convert(target_pset, paramnumber,  refvalue);

      Set_fitted_parameter_relativeshifts(parameter_shifts);
      Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
      Convert(target_pset, paramnumber,  converted_value);
      converted_parameter_shifts[fittedindex] = (converted_value - refvalue)/Parameter_scale(paramnumber);
   }
// ------------------For f1
    fittedindex = Get_fitted_parameter_index(1); // f1
    paramnumber = absolute_parameter_map[1];     // f1

    if (fittedindex >= 0)
    {
     Reset_to_initial_parameters();
     Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
     Convert(target_pset, paramnumber,  refvalue);

     Set_fitted_parameter_relativeshifts(parameter_shifts);
     Compute_state_vectors(rp, rpt, ri, rit, ro, rot, r_extra, v_extra);
     Convert(target_pset, paramnumber,  converted_value);
     converted_parameter_shifts[fittedindex] = (converted_value - refvalue)/Parameter_scale(paramnumber);
    }
  }
}
