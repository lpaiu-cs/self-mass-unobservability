// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 * 
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */
            
#include "Spline.h"
#include <cstdio>

    using namespace std;

    
Spline::Spline( long int npoints ) 
{
    n = npoints ;
    y2 = new value_type [n] ;
    xs = new value_type[n] ;
    ys = new value_type[n] ;
}

    
Spline::Spline(value_type *  x, value_type * y, long int npoints ) {
    //value_type dy0 = first_dy
    //value_type dyn = last_dy
    
    long int i = 0 ;
    long int k = 0 ;

    
    value_type p = zero ;
    value_type qn = zero ;
    value_type sig = zero ;
    value_type lastu = zero ;

    n = npoints ;   
    value_type * u ;
    u = new value_type[n];
//    value_type x[n]
    //value_type d2ys[n]
    y2 = new value_type [n] ;

    origin = undemi * ( x[0] + x[n-1] ) ;
    
    xs = new value_type[n] ;
    ys = new value_type[n] ;
    for (i = 0 ; i < n ; ++i ) {
        xs[i] = x[i] - origin ;
        ys[i] = y[i];
    }
    
                
    // Set first derivative of left edge
//     if natural == 1 :
//             d2ys[0] = 0.
//             u[0] = 0.
//         else :
//             d2ys[0] = -0.5
//             u[0] =  3. / ( x [1] - x [0] ) * ( ( ys[1] - ys[0] ) / ( x [1] - x [0] ) - dy0 ) 
//    
    // natural spline for this version
    y2[0] = zero ;
    u[0] = zero ;

    // Tridiagonal algorithm for inversion of the system
    for (i=1 ; i < n-1 ; ++i){
        sig = ( x [i] - x [i-1] ) / ( x [i + 1] - x [i-1] ) ;
        p = sig *  y2[i-1] + deux ;
        y2[i] = (sig - 1.) / p ;
        u[i] =  ( ( six * ( ( ys[i + 1] - ys[i] ) / (xs[i+1] - xs[i] ) - ( ys[i] - ys[i-1] )  /
                    ( xs[i] - xs[i-1] ) ) / ( xs[i+1] - xs[i-1] ) - sig * u[i-1] ) / p ) ;
    } 
    // Set first derivative of the right edge
//         if natural == 1 :
//             qn = 0.
//             un = 0.
//         else :
//             qn = -0.5
//             un =  3. / ( x [n-1] - x [n-2] ) * ( dyn - ( ys[n-1] - ys[n-2] ) / ( x [n-1] - x [n-2] ) )
    // Natural spline for this 
    qn = zero ;
    lastu = zero ;
        
    y2[n-1] = ( lastu - qn * u[n-2] ) / ( qn  * y2[n-2] + un ) ;
            
    // Backsubstituion loop for the trigdiagonal algorithm
    for (k = n-2 ; k > -1 ; --k) {
        y2[k] = y2[k] * y2[k + 1] + u[k] ;
    }
    
    delete[] u;
    
};



void Spline::ReSpline(value_type  x[], value_type y[]) {
    // Perform another spline with the same number of nodes.
    
    //value_type dy0 = first_dy
    //value_type dyn = last_dy
    
    long int i = 0 ;
    long int k = 0 ;

    value_type p = zero ;
    value_type qn = zero ;
    value_type sig = zero ;
    value_type lastu = zero ;   
    
    value_type * u ;
    u = new value_type[n];

    
    origin = undemi * ( x[0] + x[n-1] ) ;
    
    
    for (i = 0 ; i < n ; ++i ) {
        xs[i] = x[i] - origin ;
        ys[i] = y[i];
    }
    
                
  
    // natural spline for this version
    y2[0] = zero ;
    u[0] = zero ;
    
    // Tridiagonal algorithm for inversion of the system
    for (i=1 ; i < n-1 ; ++i){
        sig = ( x [i] - x [i-1] ) / ( x [i + 1] - x [i-1] ) ;
        p = sig *  y2[i-1] + deux ;
        y2[i] = (sig - 1.) / p ;
        u[i] =  ( ( six * ( ( ys[i + 1] - ys[i] ) / (xs[i+1] - xs[i] ) - ( ys[i] - ys[i-1] )  /
                    ( xs[i] - xs[i-1] ) ) / ( xs[i+1] - xs[i-1] ) - sig * u[i-1] ) / p ) ;
    } 
  
  
    // Natural spline for this 
    qn = zero ;
    lastu = zero ;
        
    y2[n-1] = ( lastu - qn * u[n-2] ) / ( qn  * y2[n-2] + un ) ;
            
    // Backsubstituion loop for the trigdiagonal algorithm
    for (k = n-2 ; k > -1 ; --k) {
        y2[k] = y2[k] * y2[k + 1] + u[k] ;
    }
    delete[] u;
            
};



value_type Spline::operator()(const value_type& x, const long int& kmin, const long int& kmax) {
     //   cdef np.ndarray[DTYPE_t,ndim =1] xa = self._xs
   //     cdef np.ndarray[DTYPE_t,ndim =1] ya = self._ys
  //      cdef np.ndarray[DTYPE_t,ndim =1] y2a = self._y2
    value_type xr = x - origin ; // Shift to the shifted abscissa of the spline
        //cdef int n = self._n
    value_type y = zero ;
    long int k = 0 ;
    long int khi = 0 ;
    long int klo = 0 ;
    value_type a = zero ;
    value_type b = zero ;
    value_type h = zero ;
        
    klo = kmin ;
    if (kmax < 0) 
        khi = n-1 ;
    else 
        khi = kmax ;
    
    // Look for the right place in the table by bisection
    while (khi - klo > 1 ) {
        k = (khi + klo ) / 2 ;
        if  ( xs[k] > xr ) 
            khi = k ;
        else 
            klo = k ;
    }
    
    if (xs[klo] > xr ) cout << "Error in Spline() : klo pas bon ! (klo, x, xs[klo]  " << klo << "  " << xr << "   " << xs[klo] << endl;
    if (xs[khi] < xr ) cout << "Error in Spline() : khi pas bon ! (khi, x, xs[khi]  " << khi << "  " << xr << "   " << xs[khi] << endl;
    if ( (xs[klo] > xr) or (xs[khi] < xr ) )
    {
        printf("Restarting the determination of neighbor interpolation points with klo =0 and khi = max = %li \n", n-1);
        // Look for the right place in the table by bisection
        khi = n - 1;
        klo = 0;
        while (khi - klo > 1 ) {
            k = (khi + klo ) / 2 ;
            if  ( xs[k] > xr ) 
                khi = k ;
            else 
                klo = k ;
        }
        printf(" klo -kmin %li %li \n" , klo - kmin, khi - kmax);
     }
// 1.894567926601465
//1.894555941069689 avec erreur

    h = xs[khi] - xs[klo] ;
    if (h == zero) { cout << "Bad xa input in splint, the xas must be distinct" << endl ;}
    a = (xs[khi] - xr ) / h ;
    b = ( xr - xs[klo] ) / h ;
    y = a * ys[klo] + b * ys[khi] + ( ( pow(a,trois) - a) * y2[klo]  + (pow(b,trois) - b) * y2[khi] ) * pow(h,deux) / six ;
    
    return y ;
    
}






value_type Spline::Integrate(const value_type& xstart, const value_type& xend, 
                             const long int& kminstart , const long int& kmaxstart , const long int& kminend , const long int& kmaxend ){
        //cdef np.ndarray[DTYPE_t,ndim =1] xa = self._xs
        //cdef np.ndarray[DTYPE_t,ndim =1] ya = self._ys
        //cdef np.ndarray[DTYPE_t,ndim =1] y2a = self._y2
        //cdef int n = self._n
        long int k = 0 ;
        long int khi = 0 ;
        long int klo = 0 ;
        long int khi1, klo1, klo2, khi2, k1, k2 ;
        long int i = 0 ;
        //cdef int khi1, khi2, klo1, klo2
        value_type ia = zero ;
        value_type ib = zero ;
        value_type ia3 = zero ;
        value_type ib3 = zero ;
        value_type h = zero ;
        value_type x1 = zero ;
        value_type x2 = zero ;
        value_type iy = zero ;
        value_type x = zero ;

        value_type dx = zero ;
        value_type dlo = zero ;
        value_type dhi = zero ;
        
        if ( xstart >= xend ) {
            x1 = xend - origin ;
            x2 = xstart - origin ;
            klo1 = kminend ;
            khi1 = kmaxend ;
            klo2 = kminstart ;
            khi2 = kmaxstart ;
        }
        else {
            x1 = xstart - origin ;
            x2 = xend - origin ;
            klo1 = kminstart ;
            khi1 = kmaxstart ;
            klo2 = kminend ;
            khi2 = kmaxend ;
        }
            
        klo = klo1 ;
        if ( khi1 <= 0 ) 
            khi = n - 1 ;
        else 
            khi = khi1 ;
        
        // Look for the right place in the table by bisection
        while (khi - klo > 1 ) {
            k = (khi + klo ) / 2 ;
            if  ( xs[k] > x1 ) 
                khi = k ;
            else 
                klo = k ;
        }
            
        k1 = klo ;
        
        klo = klo2 ;
        if ( khi2 <= 0 ) 
            khi = n - 1 ;
        else 
            khi = khi2 ;
            
        while (khi - klo > 1 ) {
            k = (khi + klo ) / 2 ;
            if  ( xs[k] > x2 ) {
                khi = k ;
            }
            else {
                klo = k ;
            }
        }
        k2 = khi ;
        
        if (xs[k2] < x2 or xs[k2-1] > x2 ) cout << "Error in Spline.Integrate : k2 pas bon ! (k2, x2, xs[k2]  " << k2 << "  " << x2 << "   " << xs[k2] << endl;
        if (xs[k1] > x1 or xs[k1+1] < x1 ) cout << "Error in Spline.Integrate : k1 pas bon ! (k1, x1, xs[k1]  " << k1 << "  " << x1 << "   " << xs[k1] << endl;
       
        
        
        if ( (k2 - k1) > 1 ) { // if the boundaries of the integral are not between the same two interpolation points
            // Integrate the first partial interval from x1 to x[khi]
            klo = k1 ;
            khi = k1 + 1 ;
            h = xs[khi] - xs[klo] ;
            if (h == zero ) { cout << "Bad xs input in integrate, the xas must be distinct" << endl ;} 
            x = x1 ;
            dx = (xs[khi] - x) ;
            ia = undemi * ( dx / h ) * dx ;
            ib =  dx * ( un - undemi * ( dx / h ) ) ;
            ia3 = unquart * pow( dx / h , trois) * dx ;
            ib3 = dx * ( un - troisdemis * ( dx / h ) + pow( dx / h , deux ) - unquart * pow( dx / h , trois ) );
            iy += ia * ys[klo] + ib * ys[khi] + ( ( ia3 - ia) * y2[klo]  + (ib3 - ib) * y2[khi] ) * pow( h, deux ) / six ;
        
                       
            // Integrate all the full interpolation intervals
            for (i = k1 + 1; i < k2 - 1; ++i) {
                klo = i ;
                khi = i + 1 ;
                h = xs[khi] - xs[klo] ;
           //     if ( h == zero ) { cout << "Bad xs input in integrate, the xas must be distinct" << endl ; }
                iy += undemi * h * ( ( ys[khi] + ys[klo] ) - pow( h, deux ) / douze * ( y2[khi] + y2[klo] ) ) ;
            }
            
            // Integrate the last partial interval between x[klo] and x2
            klo = k2 - 1 ;
            khi = k2 ;
            h = xs[khi] - xs[klo] ;
            if ( h == zero ) { cout << "Bad xa input in integrate, the xas must be distinct" << endl ;}
            x = x2 ;
            dx = x - xs[klo] ;
            ia = dx  * ( un - undemi * ( dx / h ) ) ;
            ib = undemi * ( dx / h ) * dx ;
            ia3 = dx * ( un - troisdemis * ( dx / h ) + pow( dx / h , deux ) - unquart * pow( dx / h , trois ) ) ;
            ib3 = unquart * dx * pow( dx / h , trois ) ;
            
            iy +=  ia * ys[klo] + ib * ys[khi] + ( ( ia3 - ia) * y2[klo]  + (ib3 - ib) * y2[khi] ) * pow( h, deux ) / six ;
        }    
        else { // if the boundaries are between the same two interpolation points 
            
            klo = k1 ;
            khi = k2 ;
            h = xs[khi] - xs[klo] ;
            dx = ( x2 - x1 ) / h ;
            dlo = ( x1 - xs[klo] ) / h ;
            dhi = ( xs[khi] - x1 ) / h ;
            
            ia = h * dx * ( dhi - undemi * dx ) ;
            ib = h * dx * ( dlo + undemi * dx ) ;
            ia3 = h * dx * ( pow( dhi, trois ) - troisdemis * pow( dhi, deux ) * dx + dhi * pow( dx, deux ) - unquart * pow(dx, trois ) );
            ib3 = h * dx * ( pow( dlo, trois ) + troisdemis * pow( dlo, deux) * dx + dlo * pow( dx, deux ) + unquart * pow( dx, trois ) );
            
            iy +=  ia * ys[klo] + ib * ys[khi] + ( ( ia3 - ia) * y2[klo]  + (ib3 - ib) * y2[khi] ) * pow( h, deux ) / six ;
        }
        
                
        if ( xstart > xend ) 
            return -iy ;
        else 
            return iy ;
     
};
