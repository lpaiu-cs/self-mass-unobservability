/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifndef IO_h
# define IO_h

#include "Constants.h"

void Read_parfile(const char * filename, value_type parameters[18], value_type parameter_shift_scales[18],
                  value_type& treference, int& integrator,int& interpsteps_per_period_i,
                  bool& roemer, bool& einstein, bool& shapiro, bool& aberration);

void Read_TNS_file(char * filename, value_type * &toas, turntype * &turns, errortype * &errors, long int& ntoas );

void Write_new_tim_file(const char * datafile_in, const char * datafile_out, const value_type * SATs);
void Write_new_tim_file(const char * datafile_in, const char * datafile_out, const value_type * SATs, const double * stddev);
void Write_new_tim_file(const char * datafile_in, const char * datafile_out, const value_type * SATs, const double * stddev, const int * mask, const int warning);

void Sortoutdatafile(const char * datafile, const int maxlinesize); // Sort out a tim file by toas and record it in datafile+"-sorted".

#endif
