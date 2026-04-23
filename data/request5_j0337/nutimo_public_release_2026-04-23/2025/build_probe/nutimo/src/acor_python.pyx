# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#coding:utf8

# distutils: language = c++

#/* This code wraps the acor routine by J. Goodman 2009
 #* 
 #* Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 #* 
 #*/
 
# To create libacorcpp.so : 
#   g++ -c -fpic acc.cpp
#   g++ -shared -o libacorcpp.so acc.o
# Then :
#   python setup_acor_python.py build_ext --inplace 
 
 
import numpy as np
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free

cdef extern from "acor.h" :
    int acor(         # The version of acor that takes the data as an argument. 
                  # This is the one that does the actual work.
                  
   double *mean,    #   A return value -- the mean of X, or possibly a slightly shorter sequence.
   double *sigma,   #   A return value -- an estimate of the standard devation of the sample mean.
   double *tau,     #   A return value -- an estimate of the autocorrelation time.
   double X[],      #   The time series to be analized.
   int L)         #   The length of the array X[].
    
    
def  acorpython(chain):
    cdef double mean
    cdef double sigma
    cdef double tau
    cdef int L 
    cdef double * X
    
    L = len(chain)
    X = <double *> malloc(sizeof(double) * L)
    
    for i in range(L):
        X[i] = chain[i]
    code = acor(&mean, &sigma, &tau, X, L)
    
    free(X)
    
    return tau, mean, sigma
    
