/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#ifndef Utilities_h
# define Utilities_h

#include <iostream>
#include <cstdio>
#include<valarray>
#include<vector>

using namespace std ;

// template<typename T>
// class myvalarray<T> : valarray<T>
// {
//
// }


int max(int m, int n) ;

int min(int m, int n) ;


void rot2D(long double angle, long double vector[],
           const int index1, const int index2 ) ;

// template <typename T> valarray<T> rodrigues(valarray<T> vector, valarray<T> axis, T angle_rad) ; // Rotate vector right-handedly around axis with angle angle_rad.


long double inversetrigo(long double cosv , long double sinv) ;

void Savetxt(const char * filename, double ** table, int nlines, int nrows );
void Savetxt(char * filename, double * table, int nlines ) ;
void Savetxt_L(char * filename, long double * table, int nlines) ;
void Savetxt_L(char * filename, valarray<long double> table, int nlines ) ;
void Savetxt_L(char * filename, long double ** table, int nlines, int nrows ) ;


void Appendtxt(const char * filename, double ** table, int nlines, int nrows );


void Count_lines_cols_in_file(char * fname, int & nlines, int & ncols, int & ncomments, int & nempty, int max_line_size=3000);
// Total number of lines in file "fname" = nlines + ncomments + nempty
// Number of columns = ncols (counted only on the first line non-commented and non-empty line)
// Comment character is "#"
// Separating character is " " or "\t" or "\n"
// max_line_size : maximum number of characters expected on a line  

void Loadtxt(char * filename, double * table, int nlines ) ;
void Loadtxt(char * filename, double ** table, int nlines, int ncols ) ;
void Loadtxt(char * filename, valarray<valarray<long double>> &table, int nlines=-1, int ncols=-1, int max_line_size =3000 ) ;

void Arange(long double table[], int ntable ,
            long double firstelement =0.L, long double spacing = 1.L ) ;


void Print_table(long double table[], long int ntable);
void Print_table(double table[], long int ntable);
void Print_table(long double table1[], long double table2[], long int ntable);
void Print_table(float table[], long int ntable);
void Print_table(long long int table[], long int ntable);
void Print_table(long int table[], long int ntable);
void Print_table(const int table[], const int ntable);
void Print_table(vector<long double> table);
void Print_table(vector<long int> table);
void Print_table(vector<int> table);
void Print_table(valarray<long double> table);
void Print_table(double ** table, long int ntable, int ncols);
void sPrint_table(char * str, int table[], int ntable);

// void Print_table(long double table[], long int ntable){
//     long int i = 0;
//     for (i=0 ; i < ntable ; ++i){
//         printf("%.19Le\n", table[i]);
//     }
//     return;
// }
//
//
// void Print_table(float table[], long int ntable){
//     long int i = 0;
//     for (i=0 ; i < ntable ; ++i){
//         printf("%f\n", table[i]);
//     }
//     return;
// }
//
//
// void Print_table(long long int table[], long int ntable){
//     long int i = 0;
//     for (i=0 ; i < ntable ; ++i){
//         printf("%lli\n", table[i]);
//     }
//     return;
// }


// void Print_table(vector<long double> table){
//     long int i = 0;
//     for (i=0 ; i < table.size() ; ++i){
//         printf("%.19Le\n", table[i]);
//     }
//     return;
// }




template < typename T > T norm3d( valarray<T> vector3 ) {
    return sqrt ( pow( vector3[0] , 2 ) + pow( vector3[1] , 2 ) +  pow( vector3[2] , 2 ) ) ;
}

template < typename T > T norm3d( T vector3[3] ) {
    return sqrt ( pow( vector3[0] , 2 ) + pow( vector3[1] , 2 ) +  pow( vector3[2] , 2 ) ) ;
}

template < typename T > T sumsquares3d( valarray<T> vector3 ) {
    return pow( vector3[0] , 2 ) + pow( vector3[1] , 2 ) +  pow( vector3[2] , 2 ) ;
}

template < typename T > T sumsquares3d( T * vector3 ) {
    return pow( vector3[0] , 2 ) + pow( vector3[1] , 2 ) +  pow( vector3[2] , 2 ) ;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T > T dotprod3d (T x[3], T y[3] )
{
    return x[0] * y[0] + x[1] *  y[1] +  x[2] *  y[2] ;
}

template <typename T> T dotprod3d (valarray< T > x, valarray< T > y )
{
    return x[0] * y[0] + x[1] *  y[1] +  x[2] *  y[2] ;
}

template <typename T> T dotprod (valarray< T > x, valarray< T > y )
{
    T result = static_cast<T>(0) ;

    for (int i = 0; i < x.size() ; i++) { result += x[i] * y[i] ; }

    return result ;
}

template <typename T> valarray<T> crossprod (valarray< T > x, valarray< T > y )
{
    valarray<T> z(3) ;

    z[0] = x[1] * y[2] - x[2] * y[1] ;
    z[1] = x[2] * y[0] - x[0] * y[2] ;
    z[2] = x[0] * y[1] - x[1] * y[0] ;

    return z;
}


template <typename T> void crossprod ( T * x,  T * y, T * resultat )
{
    resultat[0] = x[1] * y[2] - x[2] * y[1] ;
    resultat[1] = x[2] * y[0] - x[0] * y[2] ;
    resultat[2] = x[0] * y[1] - x[1] * y[0] ;
}

template <typename T> valarray<T> rodrigues(valarray<T> vector, valarray<T> axis, T angle_rad) 
// Rotate vector right-handedly around axis with angle angle_rad.
{
    T cost = cos(angle_rad);
    valarray<T> vrot = vector * cost;
    vrot += crossprod(axis, vector) * sin(angle_rad);
    vrot += axis * (dotprod3d(axis, vector) * (1-cost));
    return vrot;
}

template<typename T > long int Rank_in_sorted_array(T x , T * sortedarray, const long int& minrank , const long int& maxrank )
/*
 * Return the rank in the array at which to insert x such that the new array be still sorted (in increasing order).
 * Use a bisection method with starting points minrank and maxrank.
 */
{
    long int khi = maxrank ;
    long int klo = minrank ;
    long int k = 0;
    while (khi - klo > 1 ) {
        k = (khi + klo ) / 2 ;
        if  ( sortedarray[k] > x )
                khi = k ;
        else
            klo = k ;
    }
    return khi ;
}

template<typename T> int SearchArray(const T value, const T * array, const int array_size)
// Return the rank at which "value" is found in "array". Otherwise returns -1.
{
  for (int i =0; i < array_size; i ++)
  {
    if (array[i] == value) return i;
  }
  return -1;
}


void Correlate(double * var1, double * var2, double* correlation, int size, int size_correlation);

int Maxofarray(double * array, int size);
//
// void Maxofarray(double * array, int size, double& maxvalue, int& maxindex);

void Swap(int & val1, int & val2); // Exchange the values of val1 and val2

int fnextline(FILE * afile) ; // Read a file until the end of the current line. Return first charcter of next line or EOF.

// void Eucldiv(int numer, int denom, int quot, int rem)
// {
//     quot = static_cast<int>(trunc(numer/denom));
//     rem = 
// }

#endif
