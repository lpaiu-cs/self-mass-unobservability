// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "Constants.h"
#include<vector>
#include<valarray>


using namespace std ;





int min(int m, int n){
    if (m > n)
        return n ;
    else
        return m ;
}

int max(int m, int n){
    if (m > n)
        return m ;
    else
        return n ;
}


void Print_table(long double table[], long int ntable){
    long int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%.19Le\n", table[i]);
    }
    return;
}


void Print_table(double table[], long int ntable){
    long int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%.15e\n", table[i]);
    }
    return;
}


void Print_table(long double table1[], long double table2[], long int ntable){
    long int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%.19Le    %.19Le\n", table1[i], table2[i]);
    }
    return;
}


void Print_table(float table[], long int ntable){
    long int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%f\n", table[i]);
    }
    return;
}


void Print_table(const int table[], const int ntable){
    int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%i\n", table[i]);
    }
    return;
}

void Print_table(long int table[], long int ntable){
    long int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%li\n", table[i]);
    }
    return;
}



void Print_table(long long int table[], long int ntable){
    long int i = 0;
    for (i=0 ; i < ntable ; ++i){
        printf("%lli\n", table[i]);
    }
    return;
}


void Print_table(vector<long double> table){
    long int i = 0;
    for (i=0 ; i < table.size() ; ++i){
        printf("%.19Le\n", table[i]);
    }
    return;
}

void Print_table(vector<long int> table){
    long int i = 0;
    for (i=0 ; i < table.size() ; ++i){
        printf("%li\n", table[i]);
    }
    return;
}

void Print_table(vector<int> table){
    long int i = 0;
    for (i=0 ; i < table.size() ; ++i){
        printf("%i\n", table[i]);
    }
    return;
}


void Print_table(valarray<long double> table){
    long int i = 0;
    for (i=0 ; i < table.size() ; ++i){
        printf("%.19Le\n", table[i]);
    }
    return;
}


void Print_table(double ** table, long int ntable, int ncols)
{
  long int i = 0;
  long int j = 0;
  for (i=0 ; i < ntable ; ++i)
  {
    for (j = 0; j < ncols-1 ; j++) printf("%.15e    ", table[i][j]);
    printf("%.15e\n", table[i][ncols-1]);
  }
  return;
}


void sPrint_table(char * str, int table[], int ntable)
{
    int i = 0;
    strcpy(str,"");
    for (i=0 ; i < ntable ; ++i){
        sprintf(str, "%s%i\n", str,table[i]);
    }
    return;
}

long double inversetrigo(long double cosv , long double sinv){

//     cdef np.ndarray cosv, sinv
//     cdef np.ndarray v
//     cdef int i = 0

//     if (type(cosx) != np.ndarray):
//         cosv=np.array([0.], dtype=DTYPE)
//         sinv= np.array([0.], dtype=DTYPE)     // to be sure there is no random tail in the numbers
//         cosv=np.array([cosx], dtype=DTYPE)
//         sinv= np.array([sinx], dtype=DTYPE)
//     else:
//         cosv=np.zeros([len(cosx)], dtype=DTYPE)
//         sinv= np.zeros([len(sinx)], dtype=DTYPE)
//         cosv = cosx
//         sinv = sinx

    long double v = zero ;

//     for i in range(len(cosv)) :
        if ( abs(sinv ) > un or abs(cosv) > un ) {
            cout << "inversetrigo : Sinus ou cosinus supérieur à 1 !!!sin = " << sinv << " et cos = " << cosv << endl ;
            sinv = fmin(sinv, un);
            sinv = fmax(sinv, -un);
            cosv = fmin(cosv, un);
            cosv = fmax(cosv, -un);
            cout << "inversetrigo : utilise sin = " << sinv << " et cos = " << cosv << endl ;
        }
        if (sinv >= zero and cosv >= zero )
            v = asin( sinv ) ;
        else if (sinv >= zero and cosv < zero )
            v = acos( cosv ) ;
        else if (sinv  < zero and cosv >= zero )
            v = asin( sinv ) ;
        else if (sinv < zero and cosv < zero )
            v = pi + atan( sinv / cosv ) ;

        return v ;
}


void rot2D(long double angle, long double vector[],
           const int index1, const int index2){
    /* fait une rotation de angle sur le vecteur vecteur2D sur les composantes données par indexes
        rot2D[indexes[0]] =    | cos(angle)   -sin(angle)   0    | |x|   avec |x| = |vector[indexes[0]] |
        rot2D[indexes[1]]      | sin(angle)    cos(angle)   0    | |y|        |y|   |vector[indexes[1]] |
        Autre composantes      |   0            0           1    | |z|        |z|   |autres composantes le cas échéant |
    */
    int n = index1 ;
    int m = index2 ;
    int vectorlength = max(n,m) + 1;
    long double vecteur[vectorlength];
    long double rotmat[2][2] = { { cos(angle), -sin( angle ) } , { sin( angle ) , cos( angle ) } } ;

    memcpy(vecteur, vector, sizeof(long double) * vectorlength ) ;

    vector[n] = rotmat[0][0]*vecteur[n] + rotmat[0][1]*vecteur[m] ;
    vector[m] = rotmat[1][0]*vecteur[n] + rotmat[1][1]*vecteur[m] ;

    return ;
}




void Count_lines_cols_in_file(char * fname, int & nlines, int & ncols, int & ncomments, int & nempty, int max_line_size)
// Total number of lines in file "fname" = nlines + ncomments + nempty
// Number of columns = ncols (counted only on the first line non-commented and non-empty line)
// Comment character is "#"
// Separating character is " " or "\t" or "\n"
// max_line_size : maximum number of characters expected on a line  
{
    int l,i;
    FILE * fp = fopen(fname, "r");
    char line[max_line_size];
    char * pch;
    char com='#';
    
    nlines=0;
    nempty=0;
    ncomments=0;
    ncols=0;
    
    if (fp == NULL)
    {
        printf("\n\n Opening of trajectory file %s failed ! \n\n", fname);
        return;
    }
    l=0;
    while (fgets(line, max_line_size, fp) != NULL)
    {
        pch = strtok(line, " \t\n");
        if (pch == NULL)
            nempty++;
        else if (pch[0] == com) 
            ncomments++;
        else 
        {
            nlines++;
            if (l ==0)
            {
                i=0;
                while (pch !=NULL and pch[0] != com)
                {
                    i++;
                    pch = strtok(NULL, " \t\n");
                }
                ncols = i;
            }
        }
        l++;
    }
    fclose(fp);
};


void Savetxt_L(char * filename, valarray<long double> table, int nlines ) {
    int i ;
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    for(i=0 ; i < nlines ; ++i)
    {
        fprintf(myfile, "%.19Le\n", table[i] );
    }
    fclose(myfile);
}


void Savetxt_L(char * filename, long double * table, int nlines ) {
    int i ;
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    for(i=0 ; i < nlines ; ++i)
    {
        fprintf(myfile, "%.19Le\n", table[i] );
    }
    fclose(myfile);
}

void Savetxt_L(char * filename, long double ** table, int nlines, int nrows ) {
    int i , j;
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    for(i=0 ; i < nlines ; ++i) {
        for (j = 0; j < nrows ; ++j ) {
             fprintf(myfile, "%.19Le    ", table[i][j] );
         }
         fprintf(myfile, "\n") ;
     }
     fclose(myfile);
}

void Savetxt(const char * filename, double ** table, int nlines, int nrows ) {
    int i ,j;
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    for(i=0 ; i < nlines ; ++i) {
        for (j = 0; j < nrows ; ++j ) {
             fprintf(myfile, "%.15e    ", table[i][j] );
         }
         fprintf(myfile, "\n") ;
     }
     fclose(myfile);
}

void Savetxt(char * filename, double * table, int nlines ) {
    int i ;
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    for(i=0 ; i < nlines ; ++i)
    {
        fprintf(myfile, "%.15e\n", table[i] );
    }
    fclose(myfile);
}

void Appendtxt(const char * filename, double ** table, int nlines, int nrows ) {
    int i ,j;
    FILE * myfile ;
    myfile = fopen(filename, "a") ;
    for(i=0 ; i < nlines ; ++i) {
        for (j = 0; j < nrows ; ++j ) {
             fprintf(myfile, "%.15e    ", table[i][j] );
         }
         fprintf(myfile, "\n") ;
     }
     fclose(myfile);
}

void Loadtxt(char * filename, double * table, int nlines ) {
    int i ;
    const int sizechar = 50;
    char line[sizechar];
    FILE * myfile ;
    myfile = fopen(filename, "r") ;
//     for(i=0 ; i < nlines ; ++i)
    i=0;
    while (fgets(line, sizechar, myfile) != NULL && i < nlines)
    {
        sscanf(line, "%le", &table[i] );
        i++;
    }
    fclose(myfile);
}

void Loadtxt(char * filename, double ** table, int nlines, int ncols ) {
    int i,j ;
    FILE * myfile ;
    myfile = fopen(filename, "r") ;
    for (i=0; i < nlines ; i++)
    {
        for (j = 0 ; j < ncols-1; j++) fscanf(myfile, "%le ", &table[i][j] );
        fscanf(myfile, "%le\n", &table[i][j] );
    }
    fclose(myfile);
}

void Loadtxt(char * filename, valarray<valarray<long double>> & table, int nlines, int ncols, int max_line_size) 
{
    int i,j,l;
    char line[max_line_size];
    char * pch;
    char com='#';

// Initialise arrays
    if (nlines < 0 or ncols < 0) // Count number of lines if not given 
    {
        int ncomments, nempty, nl, nc;
        Count_lines_cols_in_file(filename, nl, nc, ncomments, nempty, max_line_size);
        if (ncols < 0) ncols = nc;
        if (nlines < 0) nlines = nl;
    }
    
    value_type components[ncols];
    
    table.resize(nlines);
    for (l=0; l < nlines; l++) table[l].resize(ncols);
    
    if (nlines ==0 ) 
    {
        printf("\nFile %s empty !\n\n", filename);
        return;
    }
    
// Open file 
    FILE * fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("\n\n Opening of file %s failed ! \n\n", filename);
        return ;
    };

// Read file 
    l=0;
    while (fgets(line, max_line_size, fp) != NULL and l < nlines)
    {
        i=0;
        pch = strtok(line, " \t\n");
        
        while (i < ncols and pch !=NULL and pch[0] != com)
        {
            sscanf(pch, "%Le", &components[i]);
            pch = strtok(NULL, " \t\n");
            i++;
        }
        if (i < ncols)
        {
            if (pch[0] == com)
                printf("Line %d : skipping due to comment sign in '%s' \n", l, filename);
            else
                printf("\n Line %d : Warning : only %d components read out of %d in '%s' !\n\n", l, i, ncols, filename);
        }
        else
        {
            for (j = 0; j < ncols; j++)
            {
                table[l][j] = components[j];
            };
            printf("lll %d\n", l);
            l++;
        }
    }
    fclose(fp);
    
    if (l < nlines) printf("Warning: only %d lines read out of %d expected in file '%s'.\n", l, nlines, filename);
};
//     //--------------------------------------------------------------------------------------------
//     int i,j ;
//     FILE * myfile ;
//     myfile = fopen(filename, "r") ;
//     long double inter;
//     table.resize(nlines);
//     for (i=0; i < nlines ; i++) table[i].resize(ncols);
//     for (i=0; i < nlines ; i++)
//     {
//         for (j = 0 ; j < ncols-1; j++) 
//         {
//             if (fscanf(myfile, "%Le ", &inter ) <1) 
//             {
//                 printf("\nError while reading file '%s' at row %d column %d\n\n", filename, i, j);
//                 return ;
//             }
//             table[i][j] = inter;
//         }
//         if (fscanf(myfile, "%Le\n", &inter ) <1) 
//         {
//             printf("\nError while reading file '%s' at row %d column %d\n\n", filename, i, ncols -1);
//             return ;
//         }
//         table[i][j] = inter;
//     }
//     fclose(myfile);
// }

void Arange(long double table[], int ntable ,
            long double firstelement =0.L, long double spacing = 1.L ){
    int i ;
    for (i = 0 ; i < ntable ; ++i) {
       table[i] = i * spacing + firstelement ;
    }
}


void Correlate(double * var1, double * var2, double* correlation, int size, int size_correlation)
{
    int k = 0;
    double mean1 =0.;
    double mean2 = 0.;

    for (k=0; k < size ; k++)
    {
        mean1 += var1[k];
        mean2 += var2[k];
    }
    mean1 /= size;
    mean2 /= size;

    for (int i = 0; i < size_correlation ; i++)
    {
        correlation[i] = 0;
        for (k = 0; k + i < size ; k ++)
        {
            correlation[i] += (var1[k] - mean1) * (var2[k+i] - mean2);
        }
        correlation[i] /= (size -i);
        correlation[i] /= correlation[0];
    }
}



// void rad2dms(long double rad, long double & d, long double & m, long double & s)
// {
//     long double angle = abs(rad) ;
//     long double minute = raddeg / 60.L ;
//     long double seconde = minute / 60.L ;
//     d = floor(angle / raddeg) ;
//     m = floor( (angle - d * raddeg) / minute) ;
//     s = ( angle - d*degre - m * minute ) / seconde ;
// }



int Maxofarray(double * array, int size)
{
  double maxvalue = array[0];
  int maxindex = 0;
  for (int k = 1; k < size ; k++)
  {
    if (array[k] > maxvalue)
      {
        maxindex = k;
        maxvalue = array[k];
      }

  }
  return maxindex;
}

// void Maxofarray(double * array, int size, double& maxvalue, int& maxindex)
// {
//   maxindex = 0;
//   for (int k = 1; k < size ; k++)
//   {
//     if (array[k] > maxvalue) maxindex = k;
//   }
//   maxvalue = array[maxindex];
// }


void Swap(int & val1, int & val2)
{
  int inter = val1;
  val1 = val2;
  val2 = inter;
}


int fnextline(FILE * afile)
// Read a file until the end of the current line. The stream indicator points on the first character of the next line.
// Return first character of next line, or EOF.
{
  int buffer = 0;
  while (buffer != '\n' and buffer != EOF) buffer = fgetc(afile);
  return buffer;
}

