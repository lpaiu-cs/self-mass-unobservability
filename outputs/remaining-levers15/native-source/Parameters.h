/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifndef Parameters_h
# define Parameters_h

#include <cstdlib>
#include "Constants.h"
#include <valarray>
#include<vector>
#include<string>

using namespace std;

void Rotate_PSB_to_SSB(valarray<value_type> & vector3, const value_type & alpha, const value_type & delta) ;
void Rotate_PSB_to_SSB(value_type * vector3, const value_type & alpha, const value_type & delta) ;

void Rotate_SSB_to_PSB(valarray<value_type> & vector3, const value_type & alpha, const value_type & delta) ;
void Rotate_SSB_to_PSB(value_type * vector3, const value_type & alpha, const value_type & delta) ;


class Parametres
{
    bool allocation_flag;
    int msg_count;
    void Read_parfile(const char * filename);


public :
    
/* Fittable parameters can be identified in 3 different ways : 
 *  - Absolute : "n" and "i", n being the "absolute" index (e.g. "0" for "f") and "i" the array index if this parameter is in fact an array (such as DMX) containing several values. 
 *  - List: the 2-dimensional absolute representation is converted into a 1D list described by a single index k = L(n,i). If an array parameter has size 0 then it is not included, therefore List describes only a subset of parameters. The size of arrays depends on whether parameters are provided in the input parfile. If 
 *  - Fitted : this is the subset of list which is fitted. It can be described by an index f such that the fth fitted parameter correspond to the list index k = F(f).
 * 
 *  
 *  Name        Absolute (N elements)   List (N elements)       Fitted  (M < N elements) 
 *  ----------------------------------------------------------------------------------
 *  p0            0(,0)                 0 = L(0,0)                   f0
 *  ...           ...                    ...                         ...
 *  pn[i]         n,i                   k = L(n,i)                    fM
 * 
 * Maps : s
 * L(n,i=0)    : k = absolute_parameter_map[n]  for i = 0 (if non-array or array not empty)
 * L(n,i)      : k = absolute_parameter_map[n] + i (if array not empty)
 * L^-1 : n,i = ParamList2AbsParam(k) (this is however a complete map)
 * F    : k = fitted_paramers(f) (also returned by Get_fitted_parameter_map)
 * F^-1 : not defined
 */

    // Corresponding Tempo2 parfile
    char t2parfile[500];
    
  // Add parameter here
 // Fittable parameters          Absolute index
    value_type f;           //        0
    value_type f1;          //        1

    value_type etap;        //        2
    value_type apsinii;     //        3
    value_type kappap;      //        4
    value_type apcosii;     //        5
    value_type tascp;       //        6
    value_type Pi;          //        7

    value_type etaB;        //        8
    value_type aBsinio;     //        9
    value_type kappaB;      //        10
    value_type aBcosio;     //        11
    value_type tascB;       //        12
    value_type Po;          //        13

    value_type mui;         //        14
    value_type muio;        //        15
    value_type Mo;          //        16

    value_type oman;        //        17
    value_type delta_oman;  //        18

    value_type dphase0;     //        19

    value_type SEP_D;       //        20
    value_type SEP_gamma;   //        21    // gammabar_p0 of Damour Scigrav2007
    value_type SEP_beta_0pp;//        22    // betabar^0_pp of Damour Scigrav2007
    value_type SEP_beta_p00;//        23    // betabar^p_00 of Damour Scigrav2007
    value_type SEP_beta_0p0;//        24    // betabar^0_p0=betabar^0_0p of Damour Scigrav2007

    // Position in the sky, interaction with tempo2
    value_type RA;          //        25
    value_type DEC;         //        26
    value_type distance;    //        27        // distance between the pulsar system barycenter and the solar system barycenter in light years
    value_type RA1;         //        28        // first derivative with time of RA in mas/yr
    value_type DEC1;        //        29        // idem of DEc in mas/yr
    value_type distance1;   //        30        // mas/yr ( ie speed / distance converted from rad to mas/yr)

    value_type  DM ;        //        31      // Dispersion measure pc.cm^-3.
    value_type DM1 ;        //        32        // First derivative of the dispersion measure DM / year
    vector<value_type>  DMX ;   //    33      //cTime dependent Corrections to DM pc.cm^-3.
    vector<value_type> FD;      //    34      // Frequency dependent DM corrections pc.cm^-3

    vector<value_type> quadrupole; // 35        // Parameters for quadrupolar moment variations

    value_type efac;        //        36        // Global prefactor for the errors. (like in tempo2)

    value_type delta_i;     //        37
    
    // Parameters of extra bodies
    vector<value_type> eta_extra;   // 38      
    vector<value_type> kappa_extra; // 39
    vector<value_type> asini_extra; // 40      
    vector<value_type> acosi_extra; // 41      
    vector<value_type> tasc_extra;  // 42      
    vector<value_type> P_extra;     // 43  
    vector<value_type> oman_extra;     // 44  
    

    // Bookkeeping variables
    vector<value_type> DMXranges; // Contains DM.size +1 MJD delimiting the time ranges where each DMX applies.
    value_type * parameters_ini; // Stores initial values of fittable parameters
    int nfitparams ;             // Number of fittable parameters


    // Add parameter here
    // Scales for relative modification of parameters : newpar = oldpar + scale_par * dpar

        value_type scale_f;           //        0
        value_type scale_f1;          //        1

        value_type scale_etap;        //        2
        value_type scale_apsinii;     //        3
        value_type scale_kappap;      //        4
        value_type scale_apcosii;     //        5
        value_type scale_tascp;       //        6
        value_type scale_Pi;          //        7

        value_type scale_etaB;        //        8
        value_type scale_aBsinio;     //        9
        value_type scale_kappaB;      //        10
        value_type scale_aBcosio;     //        11
        value_type scale_tascB;       //        12
        value_type scale_Po;          //        13

        value_type scale_mui;         //        14
        value_type scale_muio;        //        15
        value_type scale_Mo;          //        16

        value_type scale_oman;        //        17
        value_type scale_delta_oman;  //        18

        value_type scale_dphase0;     //        19

        value_type scale_SEP_D;       //        20
        value_type scale_SEP_gamma;   //        21    // gammabar_p0 of Damour Scigrav2007
        value_type scale_SEP_beta_0pp;//        22    // betabar^0_pp of Damour Scigrav2007
        value_type scale_SEP_beta_p00;//        23    // betabar^p_00 of Damour Scigrav2007
        value_type scale_SEP_beta_0p0;//        24    // betabar^0_p0=betabar^0_0p of Damour Scigrav2007

        value_type scale_RA;          //        25
        value_type scale_DEC;         //        26
        value_type scale_distance;    //        27
        value_type scale_RA1;         //        28        // first derivative with time of RA
        value_type scale_DEC1;        //        29        // idem of DEc
        value_type scale_distance1;   //        30        // iem of distance

        value_type scale_DM ;         //        31        // Dispersion measure pc.cm^-3 .
        value_type scale_DM1 ;        //        32        // First derivative of the dispersion measure DM/year.
        vector<value_type>  scale_DMX ;   //    33      // Time dependent corrections to DM pc.cm^-3.
        vector<value_type> scale_FD;      //    34      // Frequency dependent DM corrections pc.cm^-3

        vector<value_type> scale_quadrupole; // 35       // Parameters for quadrupolar moment variations
        value_type scale_efac;               //  36      // Global prefactor for the errors. (like in tempo2)

        value_type scale_delta_i;           // 37
        
        // Parameters of extra bodies
        vector<value_type> scale_eta_extra;   // 38      
        vector<value_type> scale_kappa_extra; // 39
        vector<value_type> scale_asini_extra; // 40      
        vector<value_type> scale_acosi_extra; // 41      
        vector<value_type> scale_tasc_extra;  // 42      
        vector<value_type> scale_P_extra;     // 43
        vector<value_type> scale_oman_extra;    // 44 
        
    // Add parameter here
    static const int n_absparameters =45;
    // Add parameter here
    vector<int> absolute_parameter_map; // absolute_parameter_map[absolute_index] gives the position in a 1D array of the first element of each parameter ( see absolute indexes in comments of the definitions above). If it is a simple parameter this is its index in the 1D array. If it is a vector parameter (like DMX), then it is the index of the first element of the vector and the subsequent vector.size()-1 indexes correspond to the other elements of the vector.
    vector<int> fitted_parameters ;   // List of parameters that will be fitted (ie [0, 1, 3, 20, 21], where the indexes are the 1D indexes defined from absolute_parameter_map ! SHOULD NOT BE CHANGED WITHOUT ALSO UPDATING parameter_set !

    
    int nextra; // Number of extra bodies. This is deduced from the parfile in Load_new_parfile
    int nbodies_plus_extra ; // Convenience variable equal to total number of bodies dealt with = 3 + nextra (min 3, but one can have 0 mass..)
    
  // Parameter sets : messy stuff
    int parameter_set;                // Define which parameters can be fitted, and which are derived.. Defined  in "Read_parfile".
    static const int nparameter_sets =7; // Number of parameter sets
    const int parameter_sets [nparameter_sets][n_absparameters] = {
                               {0,1,2,3,4,5,6,7,8,9,10,11,12,13,-1,15,-1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,-1,38,39,40,41,42,43,44},  // muio = mi
                               {0,1,2,3,4,5,6,-1,8,9,10,11,12,-1,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,-1,38,39,40,41,42,43,44}, // mui = (Mp + Mi) * 0.5  et muio = (Mp - Mi) * 0.5
                               {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,-1,-1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,-1,38,39,40,41,42,43,44},  // mui = mi/mp, delta_i instead of abcosi
                               {0,1,2,3,4,-1,6,7,8,9,10,11,12,13,14,-1,-1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44},  // mui = mi/mp, delta_i instead of apcosi,
                               {0,1,2,3,4,-1,6,7,8,9,10,11,12,13,14,-1,-1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44},  // mui = mi/mp, delta_i instead of apcosi, f without einstein contribution
                               {0,1,2,3,4,-1,6,7,8,9,10,11,12,13,14,-1,-1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44},  // mui = mi/mp, delta_i instead of apcosi, f without einstein contribution and f' Shcklovskii
                               {0,1,2,3,4,-1,6,7,8,9,10,11,12,13,14,-1,-1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44}  // mui = mi/mp, delta_i instead of apcosi, f without einstein contribution and Schklovskii and f' Shcklovskii with correction with respect to 5
                             };


    char ** fitparams_names; // array of the names of the fittable parameters. Set in "Load_new_parfile", have to change it when adding or removing a parameter from above.

    // Delays parameters
    value_type treference; // reference time for delays -> t0 of Integrateur is set to that
    value_type posepoch;   // Reference time (in SSB time ) for the time at which the position in the sky was determined. (used by tempo2 to determine position) .
    bool geometric; // if geometric is true roemer is ignored
    bool kopeikin;  // add the kopeikin terms to the geometric delay if activated
    bool shklovskii; // add the shklovskii part of the geometric delay up to the first k_para term if activated
    bool roemer ;
    bool einstein  ;
    bool shapiro ;
    bool aberration ;

    int truefreq ; // If 1 takes into account the constant component of the einstein delay and the shklovskii delay (giving a"truer" freq), if "0" remove linear components in some delays, if "2" assumes user provided parameters for 0 and wants to switch to 1

    // Integrator parameters
    int integrator_type;
    value_type tolint ;
    int interpsteps_per_period_i;
    float interp_margin;           // interp_margin * interpsteps_per_period_i points of margin on each side of interpolations
    
    // Keyword for special cases
    char specialcase[500];
    
    // End of parameters

  // Derived parameters
    value_type true_spinfreq;
    value_type true_spinfreq1;
    value_type SSB_RA1;
    value_type SSB_DEC1;
    value_type Mp;
    value_type Mi;
    value_type ei;
    value_type eo  ;
    value_type omp ;
    value_type omB  ;
    value_type tperii  ;
    value_type tperio;
    value_type ap  ;
    value_type aB  ;
    value_type anglii  ;
    value_type anglio ;
    value_type DopplerF ; // Doppler factor
    valarray<valarray<value_type>> gammabar;
    valarray<valarray<valarray<value_type>>> betabar;
    valarray<valarray<value_type>> Gg; 
//     value_type gammabar[3][3];
//     value_type betabar[3][3][3];
//     value_type Gg[3][3];
    vector<value_type> quadrupole_kgm2;
    value_type timeshift; // Time shift applied to all toas and interpolation to center the time span around zero
    valarray<value_type> a_extra;
    valarray<value_type> angli_extra;
    valarray<value_type> e_extra;
    valarray<value_type> om_extra;
    valarray<value_type> tperi_extra;
    valarray<value_type> M_extra;

    bool motion_changed; // True if parameters affecting orbits are changed when using a Set_parameter_relativeshift  routine


    Parametres()
    {
        allocation_flag = false;
        msg_count = 0;
    };
    Parametres(const char * parfile);
    ~Parametres()
    {
        if ( allocation_flag == true )
        {
            free(parameters_ini);
            for (int i = 0; i < nfitparams ; i++) free(fitparams_names[i]);
            free(fitparams_names);
        }
    } ;


    int Get_fitted_parameter_index(int absolute_index); // Return the index of parameter absolute_index in the fitted parameter list fitted_parameters, -1 if not fitted.
    void ParamList2AbsParam(int paramnumber, int & absparamnumber, int & index); // Inverts absolute_parameter_map
    void Insert_in_parameter_maps(int absolute_parameter_map_index, int index );
    int Is_in_parameter_set(int pset) ; // Returns -1 if all fitted parameter are in parameter_sets[pset]. Otherwise returns the absolute index of the first fitted parameter not in the list.


    value_type Parameters(int paramnumber);
    void Set_parameters(int paramnumber, const value_type newvalue);
    void Set_parameter_relativeshift(int paramnumber, const double newvalue); // Change parameter "paramnumber" (as indicated above) by a relativeshift "newvalue"
    void Set_parameters_relativeshift(double * parameter_shifts) ;
    void Set_parameters_relativeshift(vector<double>  parameter_shifts);

    void Convert(int target_pset, int param_number,  value_type & converted_value); // Run Compute_state_vectors before use !
    void Convert_in_place(int target_pset);
    void Convert_relative_shifts(int target_pset, double * parameter_shifts, double * converted_parameter_shifts);valarray<valarray<valarray<value_type> >> ns_extra;
valarray<valarray<value_type> >  rs_extra;
valarray<value_type> Ms_extra; //  Masses of extra bodies


    value_type Parameter_scale(int paramnumber);                        // Return the current scale_parameter of parameter number "paramnumber"
    void Set_parameter_scale(int paramnumber, const value_type newvalue);

    value_type operator[](const int i ) { return Parameters(i) ; };

    void Load_new_parfile(const char * filename);

    void Save_parfile(const char * filename);
    void Print();

    void Set_reference_to_current_parameters(); // Set the reference for the relative shift to the current set of parameters
    void Reset_to_initial_parameters();

    // With parameter map
    vector<int> Get_fitted_parameter_map();

    void Set_fitted_parameter_map(const vector<int> list_of_fitted_parameters); // Set a parameter map such that changed parameters when calling Set_parameters_relativeshift(vector paramshifts)
                                                                             // will be parameter[ list_of_fitted_parameters[i] ] = paramshifts[i]

    void Set_fitted_parameter_relativeshifts(vector<double> parameter_shifts) ; // uses the filter "fitted_parameters"
    void Set_fitted_parameter_relativeshifts(double * parameter_shifts);

    value_type Get_fit_parameter(int paramnumber){ return Parameters( fitted_parameters[paramnumber] ); }; // paramnumber goes from 0 to fitted_parameters.size()
    void Set_fit_parameter(int paramnumber, const value_type newvalue){ Set_parameters( fitted_parameters[paramnumber], newvalue ); }; // paramnumber goes from 0 to fitted_parameters.size()
    char* Get_fit_parameter_name(int paramnumber) { return fitparams_names[ fitted_parameters[paramnumber] ] ;};
    char * Get_parameter_name (int paramnumber) { return fitparams_names[paramnumber ] ;}; // return the name of parameter paramnumber
    // End "with parameter map"



    void Compute_state_vectors(        valarray<value_type> & rp, valarray<value_type> & rpt,
                                       valarray<value_type> & ri, valarray<value_type> & rit,
                                       valarray<value_type> & ro, valarray<value_type> & rot,
                                       valarray<valarray<value_type>> & r_extra, valarray<valarray<value_type>> & v_extra
                              );

};


#endif
