/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */


#ifndef Fittriple_h
# define Fittriple_h

# include "AllTheories3Bodies.h"
#include <valarray>
#include "Parameters.h"
#include "tempo2.h"

#ifndef MAX_FILELEN
# define MAX_FILELEN  500   /*!< Maximum filename length (should be already defined in tempo2.h     */
#endif

/*
double tempo2_getCorrectionTT(observation *obs);
void tempo2_formBats(pulsar *psr,int npsr);*/



class Fittriple : private Integrateur
{
    // Inner variables derived from "parameters" by "Parameters_to_innervariables"
    value_type spinfreq ; // Spin frequency of the pulsar in days^-1
    value_type spinfreq1 ; // Derivative of the spin frequency days^-2
    //value_type period_i ; // period of the inner binary
    value_type delaymax ; // estimate of the maximum amplitude of the delays (in days)
    value_type delaymaxold ; // to check if delaymax varies too much
    bool forcerecomputeinterp ; // To force recomputing the interpolation grid for sat -> bat in Initialize
    value_type freqshift ;  // if needed contains the frequency shift due to delays with linear time components (e.g. einstein delay). No unit, effective freq = (1-freqshift)*freq

    value_type spinPeriodmicrosec; // Conversion factor between phase and time fixed during fit

    int extrashift; // contains the shift in indices of velocity components when extra bodies are present
    
    // Interaction with tempo2
     pulsar *tempo2_psr;

    int tempo2_npsr;
    int tempo2_noWarnings;
    char  tempo2_timfile[MAX_FILELEN]; // Tim file sent to tempo2. According to Fittriple::Reconstruct it can be first sorted
    char  tempo2_parfile[MAX_FILELEN];
    value_type previous_RA; // To avoid updating bats with tempo when ra or dec do not change between two calls. Used in "Initialize()"
    value_type previous_DEC;
    value_type previous_RA1; // To avoid updating bats with tempo when ra or dec do not change between two calls. Used in "Initialize()"
    value_type previous_DEC1;
    value_type previous_distance1;
    value_type previous_DM; // To avoid updating bats with tempo when DM or DM1 do not change between two calls. Used in "Initialize()"
    value_type previous_DM1;
    vector<value_type> previous_DMX;
    vector<value_type> previous_FD;
    value_type initialpuslardirection[3];
    value_type ** pulsardirection_BAT; // same as above but for bats of the real times of arrival (and not interpolation)
    value_type ** obs_ssb_BAT; // same as above but for bats of the real times of arrival (and not interpolation)
    value_type * roemer_ss_BAT ; // same as above but for bats of the real times of arrival (and not interpolation)
    int ndebarpsr;
    value_type * tsat_testmode ; // testmode :
    value_type ** sat_to_bat ; // testmode


    void Initialize() ;

    bool first_interpolation_allocation; // True if the tinterp was never defined and allocated yet. See Fittriple-init.cpp




    // For Diagnostics
    short int getdelaysflag ;
    double * geomdelay;
    double * nogeomdelay;
    double * einsteindelay;
    double ** geom_delay_details;

    // Used by Get_global_interpolation_values
    double * diag_einstein_interp;
    double * diag_shapiro_aberration_interp;
    long int diag_interpolated_delays_n;
    long double * diag_interpolated_delays_times;
    bool get_global_interpolated_delays_flag;


public :
     bool motion_changed;
     bool remove_mean;


    // Test mode
    bool testmode ; // if true enables some comments, variables, tests...



    // Diagnostics variables
    valarray<value_type> energies;
    valarray<valarray<value_type> > center_of_mass_positions;
    valarray<valarray<value_type> > center_of_mass_impulsions;

    // for memory
    vector< long int > removed_toas ; // Contains the indices of the toas removes by "Create_mask(value_type threshold);"

    value_type timedays ; // =time from Integrateur / daysec -> time is in sec and timedays in days


    long int ninterp; // Number of interpolation times
    value_type dt_interp ; // interpolation time step in days
    value_type * tinterp ; // times of interpolation in real tiume (days in principle)

    long int ntoas ; // Number of time of arrival
    long int treference_in_interp ; // Position of treference in the interpolation times


    value_type * toas ; // Times of arrival (should be in days MJD)
    value_type * toes ; // Times of emission (should be in days MJD)
    value_type * residuals ; // Timings residuals (microseconds )
    long int * toa_in_interp ; // Position of treference in the interpolation times "ts" (inherited from Integrateur)
    turntype * turns ;
    errortype * errors ;
    errortype * weights ;

    bool disable_sats ; // if true bats are not recomputed from sats using tempo2

    value_type ** sp;   // Contains the state vertors of the pulsar sp[number of times][6]
    value_type ** si;   // idem for inner companion
    value_type ** so;   // idem for outer companion

    Parametres parameters;


    bool tracker;
    int callnumber;

    Fittriple();  // Must be followed by a call to reconstruct
    Fittriple(char * parfile, char * datafile); // Call Reconstruct
    ~Fittriple(); // a développer

    void Reconstruct_noCompute(char * parfile, char * datafile, bool sorted_datafile); // should not be called after a previous Reconstruct or Fittriple(char * parfile, char * datafile) unless memory deaollation is implemented i.e. call a destructor before reconstructing.
    void Reconstruct(char * parfile, char * datafile); // Same as above + call Compute_lnposterior(1) and Compute_turns_from_parameters() to initialize turn numbers.
    void Reconstruct(char * parfile, char * datafile, turntype * external_turns); // same but DOES NOT call Compute_lnposterior. Load external turns instead.
    void Reconstruct(char * parfile, char * datafile, bool sorted_datafile); // Same as above. if sorted_datafile == true, does not create a sorted timfile
    void Reconstruct(char * parfile, char * datafile, bool sorted_datafile, turntype * external_turns); //  Same as above. if sorted_datafile == true, does not create a sorted timfile

    double Compute_lnposterior(int fractional = 0) ; // Fittriple-compute.cpp : Compute the logarithm of the posterior probability that the parameters in parfile are correct

//------------------------------------- Parameters --------------------------------------------------
    void Set_parameters(value_type params[18]); // Fittriple-init.cpp :
    void Initialise_parameters(); // Derive all the internal parameters. To be run after Set_parameters or equivalent. In effect, interface to paramters.Compute_state_vectors.  Fittriple-init.cpp

// ----- Overloads of Set_parameters_relativeshift ---- // Fittriple-init.cpp
  // These are just quick interfaces to parameters.Set_parameters_relativeshift


    void Set_parameters_relativeshift(double * relativeshift); // Fittriple-init.cpp :
    double operator()(double * relativeshift)
    {
        vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
        Set_fitted_parameter_relativeshifts(vectrelativeshift);
        return Compute_lnposterior(0) * ntoas; // Return the actual log(proba) = -chi2/2 , not the reduced one

    }
    

    void Set_parameters_relativeshift(vector<double> relativeshift); // Fittriple-init.cpp : // The only one that really takes into account parameters.fitted_parameters (list of unblocked parameters)
                                                                                             // Interface pour minuit cpp
    void Set_fitted_parameter_relativeshifts(vector<double> parameter_shifts);  // uses the map "parameters.fitted_parameters"

   // int Get_number_of_fitted_parameters(){ return parameters.fitted_parameters.size() ; } ; // return size of the map "parameters.fitted_parameters". Max number of fittable parameters : see Get_nfitparams()
// -------------------------------------------------------------------------------------------------------

    void Estimate_parameter_shift_scales(double chi2variation=0.01);// Fittriple-init.cpp :

    void Set_aberration_delay_on(){parameters.aberration = 1;};
    void Set_aberration_delay_off(){parameters.aberration = 0;};
    void Set_geometric_delay_on(){parameters.geometric = 1;};
    void Set_geometric_delay_off(){parameters.geometric = 0;};
    void Set_einstein_delay_on(){parameters.einstein = 1;};
    void Set_einstein_delay_off(){parameters.einstein = 0;};
    void Set_shapiro_delay_on(){parameters.shapiro = 1;};
    void Set_shapiro_delay_off(){parameters.shapiro = 0;};

    void Create_mask(double threshold); //  Fittriple-init.cpp : Reset the data by removing all the toas which current residual is larger than threshold.
    void Print_mask(double threshold); //  Fittriple-init.cpp : Idem but just display the toas that should be removed and their numbers, starting at one
    void Compute_turns_from_parameters();// Fittriple-init.cpp : Replace the turn numbers by those computed from the current set of parameters. The first toa correspond to turn 0.
    void Compute_fake_BATs_from_parameters(value_type * fake_BATs);  // Fittriple-init.cpp : create fake BATs
    void Compute_fake_BATS_and_delays_from_parameters(long int nfakeBats, value_type * fake_BATs, value_type * fake_delay_geom, value_type * fake_delay_ein,
                                                                                     value_type * fake_delay_shap, value_type * fake_delay_aber) ;
    void Compute_fake_SATs_from_BATs(value_type * BATs, value_type * SATs) ;


// ------Tracker status
    void Set_tracker_on(){tracker = true; callnumber = 0;};
    void Set_tracker_off(){tracker = false;};
    int Get_tracker_status(){ return int(tracker); } ;

// ------- Theory of gravity in use
    void Set_theory(int theorynumber){parameters.integrator_type = theorynumber; integrator_type = parameters.integrator_type;};
    int Get_theory(){return parameters.integrator_type;};
// ---------

// ------- Other parameters
    void Get_t2parfile(char * t2parfile){strcpy(t2parfile, parameters.t2parfile);};  
    void Get_specialcase(char * t2parfile){strcpy(t2parfile, parameters.specialcase);};
    
// ------------ Fittriple-IO.cpp
    long int Get_number_of_toas() ;
    double Get_toa_double(long int toa_number);
    double Get_timeshift_double(){ return static_cast<double>(parameters.timeshift);};
    double Get_residual_double(long int toa_number) ;
    double Get_error_double(long int toa_number) ;
    value_type Get_sat(long int toa_number){return tempo2_psr[0].obsn[toa_number].sat; }; // return the site arrival time
    double Get_sat_double(long int toa_number){return static_cast<double>(tempo2_psr[0].obsn[toa_number].sat); }; // return the site arrival time
    double Get_frequency(long int toa_number){return tempo2_psr[0].obsn[toa_number].freq; }; // return the frequency of the observation
    double Get_delays(int delaynb, long int toa_number) ;
    void Get_fake_bats_and_delays_interp(long int nfakeBats, value_type * fake_BATs, value_type * fake_delay_geom, value_type * fake_delay_ein,
                                                                                     value_type * fake_delay_shap, value_type * fake_delay_aber)  ;
    int Get_parameter_set(){return parameters.parameter_set;};
    double Get_treference_double(){ return static_cast<double>(parameters.treference + parameters.timeshift);};
    double Get_posepoch_double(){ return static_cast<double>(parameters.posepoch + parameters.timeshift);};
    double Get_DMX_range_double(int DMXrange_number){return static_cast<double>(parameters.DMXranges[DMXrange_number]);};
    int Get_number_of_DMX(){return parameters.DMXranges.size() -1;};

    long int Get_turn(long int toa_number) {return static_cast<long int>(turns[toa_number]);};
    void Set_turns(turntype * external_turns) {memcpy(turns, external_turns, ntoas * sizeof(turntype));};

    int Get_n_absparameters(){ return parameters.n_absparameters;};
    int Get_absolute_parameter_map(int param_number){return parameters.absolute_parameter_map[param_number];};

    double Get_nfitparams() { return parameters.nfitparams;}; // Return the current number of fittable parameters
    double Get_parameter_double(int param_number );
    double Get_reference_parameter_double(int param_number);
    double Get_parameter_scale_double(int param_number ) ;
    int Get_nfitted_params() {return parameters.fitted_parameters.size();} // Return the number of fitted parameters
    int Get_list_of_fitted_parameters(int param_number) {return parameters.fitted_parameters[param_number];}; // Return the indexes of fitted parameters. paramènumber < Get_nfitted_params()
    void Get_parameter_name (int paramnumber, char * name){ strcpy(name, parameters.Get_parameter_name(paramnumber)) ;};
    int Get_nextra(){return parameters.nextra;};

    double Get_mass_double(int bodynumber); // return mass in solar masses, if bodynumber = 0 : pulsar ; 1 : inner companion ; 2 outer companion
    void Get_spinfreq_double(double & spinfreq, double & spinfreq1) // Return the true spinfreq and derivative
    {spinfreq=parameters.true_spinfreq; spinfreq1 = parameters.true_spinfreq1;};
    double Get_RA_double(){return  static_cast<double>(parameters.RA)  ;} ;
    double Get_DEC_double(){return  static_cast<double>(parameters.DEC)  ;} ;

    
    long int Get_number_of_tinterp(){ return ninterp; } ;
    double Get_tinterp_double(long int tinterp_number ) { return static_cast<double>( tinterp[tinterp_number] ) ;} ;


    long int Get_number_of_removed_toas() { return removed_toas.size() ;} ;
    long int Get_removed_toa(long int toa_number) { return removed_toas[toa_number] ;};

    double Get_pulsar_statevector_double(long int tinterp_nb, int component){return static_cast<double>( sp[tinterp_nb][component] ) ;};
    double Get_inner_statevector_double(long int tinterp_nb, int component){return static_cast<double>( si[tinterp_nb][component] ) ;};
    double Get_outer_statevector_double(long int tinterp_nb, int component){return static_cast<double>( so[tinterp_nb][component] ) ;};
    double Get_extra_statevector_double(int extra, long int tinterp_nb, int component)
    {
        if (component < 3) 
            return static_cast<double>( states[tinterp_nb][9 + extra*3 + component] * length) ;
        else if (component >= 3) 
            return static_cast<double>( states[tinterp_nb][18 + extra*3 + extrashift + component-3] * length / timescale ) ;
        else
        {
            printf("\n Error in Get_extra_statevector_double : invalid component (%d)\n", component);
            return 0.;
        }
    };
    void Compute_initial_state_vectors_double(double ** svs); // Run parameters.Compute_state_vectors() and return svs such that svs[0] = {rp, vp} svs[1] = {ri, vi} svs[2] = {ro, vo} svs[i>2] = {rextra_i, vextra_i}
    void Get_analytical_pulsar_statevector_double(long int tinterp_nb, double * statevector); // Return the analytical, newtonian, state vector computed from the orbital elements at time tinterp[tinterp_nb] with orbel2statevect
    void Get_analytical_inner_statevector_double(long int tinterp_nb, double * statevector);

// -----------------------------------------------
    //void Sortoutdatafile(const char * datafile) ; //  Fittriple-IO.cpp : Sort out a tim file by toas and record it in datafile+"-sorted".

    void Save_parfile(const char * filename) ; // Fittriple-IO.cpp : save in use parameters in a file
    void Save_data(const char * filename) ; // Fittriple-IO.cpp : save in use "tns" data in a file : " toa turn error "
    void Save_output_timing_data(const char * filename)  ;// Fittriple-IO.cpp : output "toas toes residuals"
    void Save_interp_state_vectors(const char * filename)  ;// Fittriple-IO.cpp : output "tinterp sp si so"
    void Save_masked_copy_of_tim_file(const char * masked_tim_file, const int * mask); // Fittriple-IO.cpp : Uses the sorted version of the original tim file and adds # where specified
    
    void Load_turn_numbers(char * filename) ;// Fittriple-IO.cpp : loas turn numbers from file and put them into "turns"
    void Save_turn_numbers(char * filename) ; // Fittriple-IO.cpp :save turn numbers


    void Print_initial_state_vector(); //  Fittriple-IO.cpp : Convert x0 back to physical units (m and m/s in principal) and print it
    void Print() ; // Fittriple-IO.cpp : print to stdio all the inner variables except data arrays
    void Print_tempo() ; // Print information from the tempo2 library

    void Print_tests() ; // Fittriple-IO.cpp : contains what you want to print for test purposes



    // Diagnostics
    void Compute_integrals_of_motion(); // Compute the integrals of motion at each interpolation step energies, center_of_mass_positions, center_of_mass_impulsions
    double Get_energy_double(long int tinterp_number){return static_cast<double>( energies[tinterp_number] ); }; ;
    double Get_center_of_mass_position_double(long int tinterp_number, int component){return static_cast<double>(center_of_mass_positions[tinterp_number][component] ); } ;
    double Get_center_of_mass_impulsion_double(long int tinterp_number, int component){return static_cast<double>(center_of_mass_impulsions[tinterp_number][component] ); } ;
    int Check_tinterp_grid(double tol, double fact) ; // Check accuracy of the tinterp grid
    void Test_delays();
    void Get_global_interpolation_values(long int ntimes, long double * times, double * einstein_interp, double * shapiro_aberration_interp);
    void Get_local_interpolation_values(long int ntimes, long double * relativetimes, double * geometrical_inter);
};

# endif
