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


#ifndef Diagnostics_h
# define Diagnostics_h

#include<valarray>

#include "Constants.h"

using namespace std ;


void  Generate_einstein_roemer_state_vectors_diagnostic(const value_type t[], int nt,
                                                        value_type ** sp, value_type ** si, value_type ** so ,
                                                        value_type Pext, value_type aext, value_type Pin) ;

void Compute_roemer_diagnostic(const value_type t[], value_type roemer[], int nt,
                               value_type Pext, value_type aext, value_type Pin) ;

void Compute_einstein_diagnostic(const value_type ta[], value_type einstein[], int nt,
                               value_type Pext, value_type aext, value_type Pin,
                               value_type Mi, value_type Mo,
                               value_type t0) ;


void IntegralePrems3_1PN( const value_type s0[6], const value_type s1[6], const value_type s2[6],
                         value_type m0, value_type m1, value_type m2, value_type quadrupole1,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion, // <value_type>
                         value_type &  energy  );

// -----------------  n-body cases (TODO rename IntegralePrems3 into IntegralePrems -----------------------------------------------------------------------------------------------------------
// Newtonian order, including angular momentum
void IntegralePrems_0PN(const valarray<value_type> sv, const valarray<value_type> Ms, 
                         value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy,  valarray<value_type> & angular_momentum  );

// 1PN order 
void IntegralePrems3_1PN(const valarray<value_type> sv, const valarray<value_type> Ms, 
                         value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg, const valarray<valarray<value_type>>  gammabar, const valarray<valarray<valarray<value_type>>>  betabar,

                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  );
// Version that includes the SEP violation parameters Gg, gammabar, betabar

void IntegralePrems3_1PN(const valarray<value_type> sv, const valarray<value_type> Ms, 
                         value_type quadrupole1,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  );
// Same as IntegralePrems3_1PN, just assumes plain GR.

// -----------------  end of n-body cases -----------------------------------------------------------------------------------------------------------


void IntegralePrems3_1PN(const value_type sv0[6], const value_type sv1[6], const value_type sv2[6],
                         value_type m0, value_type m1, value_type m2,  value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg, const valarray<valarray<value_type>>  gammabar, const valarray<valarray<valarray<value_type>>>  betabar,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  );


void IntegralePrems3_1PN_double ( const double s0[6], const double s1[6], const double s2[6],
                         double m0, double m1, double m2, double quadrupole1,
                         double  center_of_mass[3], double  impulsion[3],
                         double &  energy  );
// Double version of IntegralePrems3_1PN, using only basic arrays (good for interfacing with python)




void IntegralePrems3_1PN_extra(const valarray<value_type> sv, const valarray<value_type> Ms,
                         value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg, const valarray<valarray<value_type>>  gammabar, const valarray<valarray<valarray<value_type>>>  betabar,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  );

void IntegralePrems3_Newt(const value_type sv0[6], const value_type sv1[6], const value_type sv2[6],
                         value_type m0, value_type m1, value_type m2,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  );

void IntegralePrems2_NewtQuad(const value_type sv0[6], const value_type sv1[6],
                         value_type m0, value_type m1, value_type q1,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  );


// To be called together to test Delays_Brut_fullgeometric
    void  Generate_einstein_geometric_state_vectors_diagnostic(const value_type t[], int nt,
                                                        value_type Pext, value_type aext, value_type Pin, value_type ai,
                                                        value_type ra, value_type rap, value_type dec, value_type decp,
                                                        value_type ** sp, value_type ** si, value_type ** so ,
                                                        value_type * nss[3], value_type * earth_ssb[3], value_type spinaxis[3],
                                                        value_type muip = 100.L, value_type muep = 1000.L  ) ;

   void Compute_einstein_diagnostic_for_geometric(const value_type ts[], const value_type t0 , value_type einstein[], int nt,
                                     value_type Pext, value_type Pi, value_type ae, value_type ai, value_type mi, value_type me,
                                     value_type ra, value_type rap, value_type dec, value_type decp, value_type muip = 100.L, value_type muep = 1000.L ) ;

   void Compute_geometric_diagnostic(const value_type ts[], const value_type t0 , value_type geometric[], int nt,
                                     value_type Pext, value_type Pi, value_type ae, value_type ai,
                                     value_type ra, value_type rap, value_type dec, value_type decp,
                                     value_type d, value_type dprime) ;

   void Compute_shapiro_diagnostic(const value_type ts[], value_type shapiro[], int nt,
                                     value_type mi, value_type me, value_type Pext, value_type Pi, value_type ae, value_type ai,
                                     value_type ra, value_type rap, value_type dec, value_type decp,
                                     value_type muip = 100.L, value_type muep = 1000.L);

   void Compute_aberration_diagnostic(const value_type ts[], value_type aberration[], int nt,
                                     value_type spinfreq, value_type Pext, value_type Pi, value_type ae, value_type ai,
                                     value_type ra, value_type rap, value_type dec, value_type decp) ;

#endif
