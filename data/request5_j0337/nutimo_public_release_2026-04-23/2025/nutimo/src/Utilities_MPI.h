/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifdef MPIMODE
#ifndef Utilities_MPI_h
#define Utilities_MPI_h

#ifdef MPIMODE
    #include "mpi.h"
#endif
#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include "Utilities.h"
#include "acor.h"
#include <stdio.h>
#include <string.h>

void MPI_Collect_and_print(const char * strtoprint, const int size_strtoprint, const int printingproc, const int * printers, const int nprinters, const MPI_Comm MPI_COMM_printers);


#endif
#endif
