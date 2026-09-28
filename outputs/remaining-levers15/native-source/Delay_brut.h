/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifndef Delays_Brut_h
# define Delays_Brut_h

#include "Constants.h"


// void Delays_Brut(value_type toas[], value_type tis[], const value_type& t0,
//                  const int& ntoas, const int& ntis, const int& nt0,
//                  value_type sp[][6], value_type si[][6], value_type so[][6],
//                  const value_type& Mi, const value_type& Mo, const value_type& freq,
//                  bool roemer, bool einstein, bool shapiro, bool aberration ) ;

void Delays_Brut(value_type tis[], const value_type& t0,
                 const long int& ntis, const long int& nt0,
                 value_type ** sp, value_type **  si, value_type ** so,        // value_type sp[][6], value_type si[][6], value_type so[][6],
                 const value_type& Mi, const value_type& Mo, const value_type& freq,
                 bool roemer, bool einstein, bool shapiro, bool aberration,
                 value_type delay[]
                ) ;

void Delays_Brut_fullgeometric(value_type tis[], const value_type& t0,
                 const long int& ntis, const long int& nt0,
                 value_type ** sp, value_type **  si, value_type ** so, value_type ** SSB_to_PSB, value_type ** r_obs,  value_type * roemer_ss,
                 const value_type& Mi, const value_type& Mo, const value_type& freq, const value_type distance, const value_type distance_derivative,
                 bool geometric, bool einstein, bool shapiro, bool aberration,
                 value_type delay[], bool & nanflag, value_type * spinaxis = NULL
                ) ;


void Delays_Brut_geometric_local(value_type * tis, const value_type& currentBAT,  const long int& ntis,
                 value_type ** sp, value_type **  si, value_type ** so, value_type * SSB_to_PSB, value_type * SSB_to_PSB_t, value_type posepoch,
                 value_type * r_obs,
                 const value_type distance, const value_type distance_derivative,
                 value_type * delay, const bool kopeikin, const bool shklovskii, value_type test, value_type ** delay_details=NULL
                ) ;

void Delays_Brut_nogeometric(value_type * tis, const value_type& t0,
                 const long int& ntis, const long  int& nt0,
                 value_type ** sp, value_type **  si, value_type ** so, value_type * SSB_to_PSB,
                 const value_type& Mp, const value_type& Mi, const value_type& Mo, const value_type& freq,
                 bool einstein, bool shapiro, bool aberration,
                 value_type * delay , const int truefreq, value_type & moydeindt, bool & nanflag, 
                 value_type * spinaxis = NULL
                ) ;

#endif
