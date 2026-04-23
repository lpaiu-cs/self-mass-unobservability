// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

// This file is inspired by Quad4FMain.cxx coming with the Minuit distribution
// Written by Guillaume Voisin, LUTh, Observatoire de Paris, March 2016


// To run correctly first do :  export OMP_NUM_THREADS=1

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnStrategy.h"

#include "Minuit_Fittriple_interface.h"
#include "Utilities.h"
#include "IO.h"
// MODIF LG
#include "tempo2.h"

using namespace ROOT::Minuit2;


int main ( int argc, char *argv[] ) {



  unsigned int maxfcn = 10000 ; // Nombre maximal d'appel à la fonction avant que minuit s'arrête, convergence ou non.
  double maxedm = 0.1;  // critère d'arrêt : minuit stop quand il pense que la distance au minimum (edm) < maxedm * 0.0001 * upvalue. j'ai mis 10 pour que les tests soient rapides, mais 0.1 en pratique c'est mieux.
  double upvalue = 1.; // Les erreurs passées à minuit correspondent aux déplacements des paramètres tels que la fonction change de "upvalue"
                       // Ici les erreurs (passée par fcn.Initialize_minuit_parameters(parameters) ) sont toujours 1., et les objets "Fittriple" ajustent leur système de paramètre interne en fonction de "upvalue" avec "fcn.Set_upvalue(upvalue)".

  //char * parfile = "parfile-Post12fitJ0337_11741toas-Newt_dphase0-scaled";
  //char * parfile = "post13fit-isma10_fin-rescaled"; // Ce sont "mes" parfile, pas ceux de tempo
  //char * datafile = "isma10-0337.20160212-TRACK-ncyobs-1400MHz-acheval";//TNS_masked_isma6.tns"; // Donnees déjà barycentrées, le fichier comporte trois colonnes "bat turns uncertainties"
 // char * datafile = "isma10-0337.20160212-TRACK-ncyobs-1400MHz-fin";
//  char * datafile = "isma10-0337.20160212-TRACK-ncyobs-1400MHz";

  char * parfile;
  char * datafile ;
  char * turnfile = NULL;
  char * fileout ;
  char postfitdat[100] = "postfit-toa_toe_residu-" ;
  char postfitpar[100] = "parfile-";
  char fakebats[150] = "Fake_BATS-";
  char fakesats[150] = "Fake_SATS-";

//   cout << argc << "   argv  " << argv[1] << endl;

    if ( argc == 4 )
    {
        parfile = argv[1] ;
        datafile = argv[2];
        fileout = argv[3];
    }
    else if (argc == 5)
    {
        parfile = argv[1] ;
        datafile = argv[2];
        turnfile = argv[3];
        fileout = argv[4];
    }
    else
    {
        printf("Syntax is: \n Minuit_fit.exe parfile datafile output_file_suffix \n or \n Minuit_fit.exe parfile datafile turnfile output_file_suffix");
        return 1;
    }

    // Names of the output files
    strcat(postfitdat, fileout);
    strcat(postfitdat,".dat") ;
    strcat(postfitpar, fileout) ;
    strcat(fakebats, parfile) ;
    strcat(fakebats, "-") ;
    strcat(fakebats, datafile) ;
    strcat(fakesats, parfile) ;
    strcat(fakesats, "-") ;
    strcat(fakesats, datafile) ;

    cout << "********************** Syntax : 'program parfile datafile outfilesid'" << "************************" << endl;


  // fcn is the function to pass to Minuit (it is defined by a class inheriting Minuit's FCNBase

  Fittriple_with_Minuit fcn(parfile, datafile, postfitpar) ;

 // ***** To enable / disable parameters during the fit ******
//   This avoid fitting on parameter 19 i.e. SEP_D
  vector<int> fittedpar;
//   for (int i = 0; i < 19 ; i ++) fittedpar.push_back(i);
//   fittedpar.push_back(20);
//   fittedpar.push_back(21);
// //  fittedpar.push_back(22); // distance
//   fittedpar.push_back(23);  // RA1
//   fittedpar.push_back(24); // DEC1
//   //   // fittedpar.push_back(25); // distance1
//   fittedpar.push_back(26); // DM
//   fittedpar.push_back(27); // DM1
//     fcn.Set_parameter_map(fittedpar);


//***** Uncomment below lines to switch from Newtonian parfile to a 1PN  one.
//(The turn numbers are calculated with Newtonian if parameters of the originial parfile are correct in this theory and a 1PN  fit is peformed)
//  cout << "Chi2 value before settheory " << fcn() << std::endl ;
//  fcn.fitwithminuit[0].Set_theory(1);
//  cout << "Chi2 value after settheory " << fcn() << std::endl ;
//*****------------------------------

  if ( turnfile != NULL )
  {
      printf("\nLoading turns from %s\n\n", turnfile);
      fcn.fitwithminuit[0].Load_turn_numbers(turnfile);
  }


  //***  Uncomment to set the value "up" as defined by Minuit and run "Fittriple.Estimate_parameter_shift_scales(upvalue)" to adjust the internal parameter system of fittriple to this value.
  // fcn.Set_upvalue(upvalue);




  //**** This recovers a previous crash. SHOULD BE COMMENTED MOST OF THE TIME !
//     double arrayparamsreco[] = {-0.908294 , 0.534511 , 0.0639505 , -0.0560163,  0.0302361 , 0.194546 , 0.220091 , -0.0767039 , -4.1409e-09 , -3.23335 , 0.00832329 , -1.60192 , -0.92846 , -0.155919 , -0.257937 , 0.0443228  ,0.219591 , -0.9851 , -0.68663 , -5.16677,  8.97969};
//     vector<double> paramsreco(arrayparamsreco, arrayparamsreco + sizeof(arrayparamsreco) / sizeof(double) );
//     fcn.fitwithminuit[0].Set_parameters_relativeshift( paramsreco );
//     cout << "Chi2 value before fit " << fcn() << std::endl ;
//   fcn.Save_parfile("parfile_avant_minuit-reco2");
//    fcn.Set_upvalue(upvalue);
//
//    double arrayparamsreco2[] = {   1.01563 , -0.72606 , 0.00285762 , 0.0462593 , 0.0323343 , 0.0317392  ,-0.151921,  -0.230423 , -1.17674e-09 , -0.676175,  -0.00275716,
//        -0.213696,  -0.912556  ,-0.254665 , -0.205754,  0.0969866,  -0.0674727,  -0.771582,  0.459456 , -6.36218 , 4.9577};
//     vector<double> paramsreco2(arrayparamsreco2, arrayparamsreco2 + sizeof(arrayparamsreco2) / sizeof(double) );
//     fcn.fitwithminuit[0].Set_parameters_relativeshift( paramsreco2 );
//               cout << " lksqfhlfkh " << endl;
//         cout << "Chi2 value before fit " << fcn() << std::endl ;
//    fcn.Save_parfile("parfile_avant_minuit-reco2");
//
//    fcn.Set_upvalue(upvalue);
//
//
  // End of recovery


  fcn.fitwithminuit[0].parameters.Set_reference_to_current_parameters();
  fcn.Save_parfile("parfile_avant_minuit");


  //Definition of parameters
  MnUserParameters parameters;
  fcn.Initialize_minuit_parameters(parameters);

  // Stragégie de Minuit
  unsigned int strategie = 1 ;
  MnStrategy mnstrategie(strategie) ; // 1 est la stratégie par défaut, 2 est la stratégie lente mais plus sûre, à envisager...



  fcn.Print(); // Print the internal state of the fittriple object

  fcn.fitwithminuit[0].Print_initial_state_vector();

  cout << "********* Récapitulatif*********" << std::endl;
  cout << "parfile = " << parfile << std::endl;
  cout << "datafile = " << datafile << std::endl;
  cout << "maxfcn = " << maxfcn << std::endl;
  cout << "maxedm = " << maxedm << std::endl;
  cout << "upvalue = " << upvalue << std::endl;
  cout << "strategie = " << strategie << std::endl;
  cout << "********************************" << std::endl;


  printf("Logposterior value before fit = %.15e \n", fcn() ) ;
  if (fcn.prior_on == true ) printf(" > Logposterior includes priors.\n");
  printf("\n*** Parameter map : ");
  for (int i=0; i < fcn.fitwithminuit[0].Get_nfitted_params() ; i++) printf("%i ", fcn.fitwithminuit[0].Get_list_of_fitted_parameters(i));
  printf("\n");

  fcn.fitwithminuit[0].Save_output_timing_data("datafile_avant_minuit.dat");
  fcn.fitwithminuit[0].Save_turn_numbers("turns.dat");

//********* Uncomment to obtain ideal BATs and SATs (no kopeikin) generated from the model
  // long double * toas;
  // long double * sats_toas;
  // toas = (long double *) malloc(sizeof(long double) * fcn.fitwithminuit[0].ntoas);
  // sats_toas = (long double *) malloc(sizeof(long double) * fcn.fitwithminuit[0].ntoas);
  // fcn.fitwithminuit[0].Compute_fake_BATs_from_parameters(toas);
  // Savetxt_L(fakebats, toas,  static_cast<int>(fcn.fitwithminuit[0].ntoas)) ;
  //
  //
  // fcn.fitwithminuit[0].Compute_fake_SATs_from_BATs(toas, sats_toas);
  // Write_new_tim_file(datafile, fakesats, sats_toas);
  //
  // free(sats_toas);
  // free(toas);
  //


//**********---------------------------------

  // Fit !
   MnMigrad migrad(fcn, parameters, mnstrategie); // voir fcn::operator() pour plus d'infos sur une possibilité d'intéraction avec tempo
   FunctionMinimum min = migrad(maxfcn, maxedm);

// Récupère le résultat et enregispartre le nouveau parfile :
  fcn.Set_parameters_from_minuit(min.UserParameters());
  fcn.Save_parfile(postfitpar) ;

// Sauvegarde les résidus postfit
  fcn.fitwithminuit[0].Save_output_timing_data(postfitdat);

// Affiche le résultat de minuit
  std::cout<<"Minimum: "<<min<<std::endl;

//      try to run hesse
//      MnHesse hesse;
//      hesse( fcn, min);
//      std::cout<<"minimum after hesse: "<<min<<std::endl;
   }
