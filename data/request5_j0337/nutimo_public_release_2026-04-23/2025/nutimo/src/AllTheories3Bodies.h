/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

# ifndef AllTheories3Bodies_h
# define AllTheories3Bodies_h
/*
******* !!!!! **********
Erreur avec gcc 4.8 pour la génération du boost::python !
Marche avec gcc 4.6, voir http://itk-users.7.n7.nabble.com/ITK-Build-Failing-td32665.html
pour la marche suivie pour faire cohabiter les deux versions
Erreur renvoyée par py++ :
gccxml-0.9/GCC/4.8/xmmintrin.h:1215: error: '__builtin_ia32_pause' was not declared in this scope


Ce code permet d'effectuer l'intégration des équations du mouvement de trois corps en intéraction gravitationnelle au premier ordre post-newtonien de la relativité générale.
cf Will, "Theory and experiment in gravitational physics".
Il s'agit d'une modification de Newtonian3Bodies.cpp et Newtonian3Bodies.hpp où la routine " Integrateur::rhs" est modifiée.
*/
//#include <valarray>
#include <boost/numeric/odeint.hpp>
#include "Constants.h"
#include <valarray>

   using namespace std;
    using namespace boost::numeric::odeint;
// Déclaration des types


// Type of the state vector
typedef vector<value_type> state_type;
/* components      object property          ( 0:3 means from 0 included to 3 excluded )
    0:3          pulsar position
    3:6          inner companion position
    6:9          outer companion position
    9:12         pulsar velocity
    12:15        inner companion velocity
    15:18        outer companion velocity
*/


// Type of the time-derivative of the state vector
typedef vector<value_type> deriv_type;


// A tag type characterizing the category of the stepper. This type must be convertible to dense_output_stepper_tag.
typedef dense_output_stepper_tag stepper_category ;

// Type of the stepper
typedef bulirsch_stoer_dense_out<
               state_type, value_type,
               deriv_type, value_type,
               range_algebra,   // default : Algebra = range_algebra
               default_operations, // default : Operation = default_operation
               initially_resizer>  // default : Resizer = initially_resizer
                                    Stepper;

/*bulirsch_stoer_dense_out(value_type eps_abs = 1E-6, value_type eps_rel = 1E-6,
                         value_type factor_x = 1.0,
                         value_type factor_dxdt = 1.0,
                         time_type max_dt = static_cast< time_type >(0),
                         bool control_interpolation = false);
     http://www.boost.org/doc/libs/1_60_0/libs/numeric/odeint/doc/html/boost/numeric/odeint/bulirsch_stoer_dense_out.html
*/
// typedef runge_kutta_dopri5< state_type, value_type,
// deriv_type, value_type,
// range_algebra,   // default : Algebra = range_algebra
// default_operations, // default : Operation = default_operation
// initially_resizer>  // default : Resizer = initially_resizer
// errstepper;

typedef pair <value_type, value_type> intervals_type;

//namespace IntegratorScope{




value_type trimfloat(value_type x, int mantissa);


/* Here is the syntax of the bulirsch stoer stepper constructor
    bulirsch_stoer_dense_out(
        value_type eps_abs = 1E-6 , value_type eps_rel = 1E-6 ,
        value_type factor_x = 1.0 , value_type factor_dxdt = 1.0 ,
        bool control_interpolation = false )
*/

// Class of the GR integrator
class Integrateur
{

// protected :
public:
    long int n;                      //Number of times
    int roemer_second_order;    // Number of the indice of the state vector to record the derivatice at each time, allows do get the second derivative for computation of roemer delays
    //Stepper stepper;            // Stepper from boost
    state_type x0;              // Initial State vector at initial time
    value_type dt0;             // Initial time step for the stepper
    value_type t0;              // Reference time.
    vector<value_type> ts;       // Times when the state vectors must be recorded
    vector<value_type> tneg;       // Times below t0
    long int ntneg ;         // Size of tneg. Gives the position the time preceding t0 in the table ts or the position of t0 if t0 belongs to ts. Affected by the children class Fittriple.
    vector<value_type> tpos;       // Times above or equal to t0
    vector<state_type> states;  // State vectors of the differential system at each time ts
    vector<value_type>  derivatives;    // Derivatives recorded along integration if roemer_second_order > -1
    vector<value_type> retro_derivatives; // used as intermediate for derivatives when integrating backward

    // Dimensional coefficients of the differential equation system
    value_type mass; // should be in kg most of the time = 1 solar mass
    value_type length;  // should be in meters most of the time
    value_type timescale;     // should be in seconds most of the time

    static const int n3body =3; // Just a constant pointing to the number of the 3 main bodies i.e. "3"
// Body masses (in solar masses)
    value_type int_M0, int_M1, int_M2;
    
// Masses of extra bodies 
    valarray<value_type> int_M_extra;

// SEP violation parameters
    value_type int_SEP_D; // In newtonian : M0_gravitational = (1 + SEP_D) M0_inertial
//     value_type int_gammabar[n3body][n3body];
//     value_type int_betabar[n3body][n3body][n3body];
//     value_type int_Gg[n3body][n3body];
    valarray<valarray<value_type>> int_gammabar;
    valarray<valarray<valarray<value_type>>> int_betabar;
    valarray<valarray<value_type>> int_Gg;
    vector<value_type> int_quadrupole;

    // Tolerance of the boost integrator
    value_type tolint;

    int integrator_type;

//optimize
    // void SecondDerivative(state_type& xvect); //
    // void Retro_SecondDerivative(state_type& xvect);


//     public :

    static void retro_rhs(const state_type& , deriv_type& , const value_type& );
    static void rhs(const state_type& , deriv_type& , const value_type& ); //rhsOH(const state_type& , deriv_type& , const value_type& );
    
    static void retro_rhs_GR(const state_type& , deriv_type& , const value_type& );
    static void rhs_GR(const state_type& , deriv_type& , const value_type& );
    
    static void retro_rhs_GR_nbody(const state_type& , deriv_type& , const value_type& );
    static void rhs_GR_nbody(const state_type& , deriv_type& , const value_type& );
    
    static void retro_rhs_0PN_nbody(const state_type& , deriv_type& , const value_type& );
    static void rhs_0PN_nbody(const state_type& , deriv_type& , const value_type& );
    
    static void retro_rhs_GR_quadrupole(const state_type& , deriv_type& , const value_type& );
    static void rhs_GR_quadrupole(const state_type& , deriv_type& , const value_type& );
    
    static void retro_rhs_Newt(const state_type& , deriv_type& , const value_type& );
    static void rhs_Newt(const state_type& , deriv_type& , const value_type& );
    
    static void retro_rhs_NewtQuad(const state_type& , deriv_type& , const value_type& );
    static void rhs_NewtQuad(const state_type& , deriv_type& , const value_type& );
    
    static void rhs_biKeplerian(const state_type& , deriv_type& , const value_type& );
    static void retro_rhs_biKeplerian(const state_type& , deriv_type& , const value_type& );
    
    void Retro_Integre();

    void Integrate_Allways();

    void Initialise();

    Integrateur() :
        tolint(10E-14),
        n (0),
        t0(0.),
        dt0 (0.),
        mass (1.),
        length (1.),
        timescale (1.),
        int_M0 (0.),
        int_M1 (0.),
        int_M2 (0.),
        int_SEP_D(zero),
        integrator_type (0),
        roemer_second_order (-1)
        {
            Initialise();
            x0.resize(18); states.push_back(x0) ;
            int_quadrupole=vector<value_type>(4,zero);
            int_SEP_D =zero;
            int_Gg.resize(n3body);
            int_gammabar.resize(n3body);
            int_betabar.resize(n3body);
            for (int i= 0 ; i < n3body ; i++)
            {
                int_Gg[i].resize(n3body);
                int_gammabar[i].resize(n3body);
                int_betabar[i].resize(n3body);
                for (int j= 0 ; j < n3body ; j++) int_betabar[i][j].resize(n3body);
            }
            for (int i = 0 ; i < n3body ; i++)
            {
                for (int j = 0 ; j < n3body ; j++)
                {
                    int_Gg[i][j] = un;
                    int_gammabar[i][j] = zero;
                    for (int k = 0 ; k < n3body ; k++)
                    {
                        int_betabar[i][j][k] = zero;
                    }
                }
            }
        };

    Integrateur ( const long int NbTimes, value_type tini, value_type dtini, value_type AdimM, value_type AdimL, value_type AdimT,
                                                    value_type Mass0, value_type Mass1, value_type Mass2, int roemer2order,
                                                    value_type tolerance=10E-10, int integrator=0) :
        tolint(tolerance),
        //stepper(tolerance, tolerance, 1. ,1., true),
        n (NbTimes),
        t0(tini),
        dt0 (dtini),
        mass (AdimM),
        length (AdimL),
        timescale (AdimT),
        int_M0 (Mass0),
        int_M1 (Mass1),
        int_M2 (Mass2),
        int_SEP_D(zero),
        integrator_type (integrator),
        roemer_second_order (roemer2order)
     {
        Initialise();
        int_quadrupole=vector<value_type>(4,zero);
         int_SEP_D =zero;
//             for (int i = 0 ; i < n3body ; i++)
//             {
//                 for (int j = 0 ; j < n3body ; j++)
//                 {
//                     int_Gg[i][j] = un;
//                     int_gammabar[i][j] = zero;
//                     for (int k = 0 ; k < n3body ; k++)
//                     {
//                         int_betabar[i][j][k] = zero;
//                     }
//                 }
//             }
        ts.reserve(n);
        tpos.reserve(n); //tneg should in principle have a small number of component
        states.reserve(n);
        if (roemer2order > -1) derivatives.reserve(n);
    };

//    void Integre();

////optimize    void Append_initial_component(const value_type& state_component);// Do not use with standard constructor because it preallocates space for x0
////optimize    void Append_time(const value_type& t);
    void Clear();
    long double Get_time(long int time_nb){ return ts[time_nb];};
    long double Get_state(long int time_nb, long int component_nb) {return states[time_nb][component_nb];};
    long double Get_derivatives(long int time_nb) {return derivatives[time_nb];};
};


#endif
