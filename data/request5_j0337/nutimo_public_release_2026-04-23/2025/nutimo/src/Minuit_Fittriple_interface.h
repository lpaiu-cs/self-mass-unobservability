/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */
// This file is inspired by Quad4.h coming with the Minuit distribution
// Written by Guillaume Voisin, LUTh, Observatoire de Paris, March 2016


#ifndef Minuit_Fittriple_interface_h
# define Minuit_Fittriple_interface_h


#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FCNGradientBase.h"

#include "Fittriple.h"

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
  #define omp_get_num_threads() 1
#endif

namespace ROOT {

   namespace Minuit2 {




// # FIRSTPRIVATE(fitwithminuit)

class Fittriple_with_Minuit : public FCNBase {

public:
  int * callnb;
  char * outparfile;

  bool prior_on = true;
  // Position from GAIA DR2 (radians) Warning : POSEPOCH = MJD2015.5
  value_type ramean = 9.500283096783616e-01;
  value_type decmean = 3.011411347204199e-01;
  value_type raerror = 8.945972920351486e-10; // =  0.184523937110807 mas
  value_type decerror = 9.302791001352137e-10; // = 0.191883838345113 mas
  value_type posprior = 2;

  // Tangential proper motion (mas/yr) from Gaia DR2. Warning : POSEPOCH = MJD2015.5
  value_type pmramean = 4.81377140723352;
  value_type pmdecmean = -4.42182046193475;
  value_type pmraerror = 0.498024673422251;
  value_type pmdecerror = 0.42770795176738;
  value_type tpmprior = 2;

  // Radial velocity from Kaplan et al. 2014
  value_type rvkaplan = 29.7e+03; //m/s
  value_type rverrorkaplan = 0.9e+03; // m/s
  value_type rvprior = 1;

  // Prior centered on Kaplan 2014 with ~2error kaplan.
  value_type dmean = 4350.;//4604.639595490277 ; // lyr
  value_type derror = 250.; //1558.9119367898795; // lyr
  value_type dprior = 2;


    // -------------------- Basic Minuit interface  ------------------------------------------------

  Fittriple_with_Minuit(char * parfile, char * datafile, char * parfile_postfit)
    {

        number_of_ompthreads  = omp_get_max_threads();
        fitwithminuit= new Fittriple[number_of_ompthreads];
        upvalue=1.;
        cout << std::endl << "******** Max number of threads (OpenMP) threads expected in Fittriple/Minuit interface : " << number_of_ompthreads << " **********" << std::endl << std::endl;
             //   omp_set_num_threads(1); // Make sure the number of threads is one  because tempo2 does not support more

        for (int i = 0 ; i < number_of_ompthreads; i++)
        {
          fitwithminuit[i].Reconstruct(parfile, datafile);
        };

        callnb = new int[1];
        *callnb = 0;

        outparfile = parfile_postfit;

        return;
    }

  ~Fittriple_with_Minuit() {delete callnb; delete[] fitwithminuit ;}; //delete fitwithminuit ; }  // Ce truc est responsable d'un segmentation fault je ne comprends pas trop pourquoi, mais on s'en fout un peu aussi...



  double operator()(const std::vector<double>& par) const {

    /* Voici comment j'imagine que l'intéraction avec tempo pourrait se faire ici :
     * Set_tempo_with_par(par) // Transmet à tempo les nouveaux paramètres initiaux, en particulier position
     * Run tempo
     * Get_results_from_tempo // Récupère les bats, ainsi que la position de la Terre / SSB pour calcul des termes de retard de Kopeikin
     * Attention !! si openMP a plus de 1 thread, ils pourraient tous modifier tempo en même temps et faire n'importe quoi ...(export OMP_NUM_THREADS=1 avant d'exécuter pour pas être embêté )
    */

     double chi2 = 0;
     //double chi2test = 0.;
     bool thread = false;

     *callnb = *callnb +1 ;

      if ( ( omp_get_thread_num() == 0 ) and ( (*callnb)% static_cast<int>(ceil(20./number_of_ompthreads)) == 0 ) )
      {
          cout << "** Saving current state to '"<< outparfile << "'" << endl;
          fitwithminuit[0].Save_parfile(outparfile) ;
          cout << " *** saved params = \n " ;
           fitwithminuit[0].Print() ;
          cout << std::endl;
           fitwithminuit[0].Print_tempo() ;
           cout << "\n";
           fitwithminuit[0].Print_initial_state_vector() ;
      }

    for (int i = 0 ; i < number_of_ompthreads; i++)
        {
            if (omp_get_thread_num() == i )
            {
                thread = true ;
                fitwithminuit[i].Set_fitted_parameter_relativeshifts(par) ;
                if (fitwithminuit[i].parameters.motion_changed == false && *callnb > 1)
                    {fitwithminuit[i].motion_changed = false;
                    //chi2test = fitwithminuit[i].Compute_lnposterior(0) ;
                   // fitwithminuit[i].motion_changed = true;
                 //   chi2 = fitwithminuit[i].Compute_lnposterior(0) ;
               //     printf("\n motion changed false : %.15e %.15e %.15e ", (chi2 -chi2test)/chi2, chi2, chi2test);
                    }
                else
                    {fitwithminuit[i].motion_changed = true ;}// chi2 = fitwithminuit[i].Compute_lnposterior(0) ;}
                chi2 = fitwithminuit[i].Compute_lnposterior(0) ;
                if (prior_on ==true)
                {
                    if ( ( omp_get_thread_num() == 0 ) and ( (*callnb)% static_cast<int>(ceil(20./number_of_ompthreads)) == 0 ) )
                        chi2 += Compute_prior(i, true);
                    else
                        chi2 += Compute_prior(i, false);
//                   if ( ( omp_get_thread_num() == 0 ) and ( (*callnb)% static_cast<int>(ceil(20./number_of_ompthreads)) == 0 ) ) Print_prior();
                }
                printf("Call nb : %i | Thread nb : %i | chi2 = %.15f | ", *callnb, i, chi2) ;
                cout << " params = " ;
                for (int j = 0 ; j < fitwithminuit[i].parameters.fitted_parameters.size() ; j++) cout << par[j] << "  " ;
                cout << std::endl;
            }
        }

    if (thread == false) cout << "************************************* !!! Too many threads in OpenMP for Fittriple_with_Minuit !!! *********************" << std::endl ;

    return - chi2 ;

  };



    void Set_parameters_from_minuit(const MnUserParameters & mnparameters)
       {

        /* Voici comment j'imagine que l'intéraction avec tempo pourrait se faire ici :
         * Set_tempo_with_par(par) // Transmet à tempo les nouveaux paramètres initiaux, en particulier position
         * Run tempo
         * Get_results_from_tempo // Récupère les bats, ainsi que la position de la Terre / SSB pour calcul des termes de retard de Kopeikin
        */
        std::vector<double> params;
        for (unsigned int i = 0; i < fitwithminuit[0].parameters.fitted_parameters.size() ; i++) params.push_back( mnparameters.Value(i) ) ;
        fitwithminuit[0].Set_fitted_parameter_relativeshifts(params);

        return ;
      };

  double Up() const {return upvalue;}

  // ----------------------------------------------------------------------------------------------------

    double operator()() const  // Compute the chi2 of the current state of the Fittriple object (on thread 0 if relevant).
    {
     double chi2 = 0;
    fitwithminuit[0].motion_changed = true;
     chi2 = fitwithminuit[0].Compute_lnposterior(0);
     if (prior_on ==true) chi2 += Compute_prior(0);
     return - chi2 ;
    };


  void Print() {      fitwithminuit[0].Print(); Compute_prior(0,true); return; };

  void Initialize_minuit_parameters(MnUserParameters & mnparameters)
  {
      for (int i = 0; i < fitwithminuit[0].parameters.fitted_parameters.size() ; i++)
      {
          mnparameters.Add(fitwithminuit[0].parameters.Get_fit_parameter_name(i) , 0. , 1. );
      };

      return;
  };

  void Save_parfile(const char * parfile ) // Save dans un parfile les paramètres du thread principal
  {
      fitwithminuit[0].Save_parfile(parfile) ; return;
};


  void Set_upvalue(double value) // Ajuste l'objet fittriple pour "upvalue". Propage ça sur tous les threads si nécessaire.
  {
      upvalue = value;
      fitwithminuit[0].Estimate_parameter_shift_scales(upvalue);
      for (int i = 0 ; i < number_of_ompthreads ; i++)
      {
          for (int j = 0 ; j < fitwithminuit[0].parameters.nfitparams ; j++) fitwithminuit[i].parameters.Set_parameter_scale(j, fitwithminuit[0].parameters.Parameter_scale(j) ) ;
      }
      return;
  };

  double Get_upvalue() {return upvalue;};

  Fittriple * fitwithminuit;//= new Fittriple[number_of_ompthreads];


  void Set_parameter_map(vector<int> parameter_map)
  {
      for (int i = 0 ; i < number_of_ompthreads; i++)
        {
          fitwithminuit[i].parameters.Set_fitted_parameter_map(parameter_map);
        };
      return;
  };

  double Compute_prior(const int omp_i, bool print=false) const 
  {
    double ragauss = -0.5*pow((fitwithminuit[omp_i].parameters.RA - ramean)/(raerror*posprior),2);
    double decgauss = -0.5*pow((fitwithminuit[omp_i].parameters.DEC - decmean)/(decerror*posprior),2);
    double rvgauss = -0.5 * pow((fitwithminuit[omp_i].parameters.distance1* radmasdeg * fitwithminuit[omp_i].parameters.distance *clight - rvkaplan)/(rverrorkaplan*rvprior),2);
    double dgauss = -0.5 * pow((fitwithminuit[omp_i].parameters.distance - dmean)/(derror*dprior),2);
    double pmragauss = -0.5 * pow((fitwithminuit[omp_i].parameters.RA1 - pmramean)/(pmraerror*tpmprior),2);
    double pmdecgauss = -0.5 * pow((fitwithminuit[omp_i].parameters.DEC1 - pmdecmean)/(pmdecerror*tpmprior),2);
    double totalprior = (ragauss + decgauss + rvgauss + dgauss + pmragauss + pmdecgauss)/ fitwithminuit[omp_i].ntoas;
    if (print==true)
    {
      printf("\nPrior :  RA   DEC  RV   D   PMRA PMDEC\n");
      printf(" %.2e    %.2f %.2f %.2f %.2f %.2f %.2f\n\n", -totalprior, 2*sqrt(-ragauss),2*sqrt(-decgauss), 2*sqrt(-rvgauss), 2*sqrt(-dgauss), 2*sqrt(-pmragauss), 2*sqrt(-pmdecgauss)); 
    }
    
    return totalprior;
  };

private:

    double upvalue;
    int number_of_ompthreads  ;//= omp_get_max_threads();

};


   }
}

#endif
