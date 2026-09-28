// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 *
 */

#include "Constants.h"
# include <cmath>
#include <iostream>
#include <string>
#include <valarray>
#include "Utilities.h"

using namespace std;

void  Generate_einstein_roemer_state_vectors_diagnostic(const value_type t[], int nt,
                                                        value_type ** sp, value_type ** si, value_type ** so ,
                                                        value_type Pext, value_type aext, value_type Pin
                                                        ) {
 //       cdef np.ndarray[DTYPE_t, ndim = 2] statevects = np.ndarray((t.size, 18), dtype = DTYPE)
  //      cdef np.ndarray[DTYPE_t, ndim = 1] d2roeXc = np.ndarray((t.size), dtype = DTYPE)
        //cdef np.ndarray[DTYPE_t, ndim = 1] ni = np.array([1.,2.,3.], dtype = DTYPE ) / np.sqrt(14.)
        //cdef np.ndarray[DTYPE_t, ndim = 1] no = np.array([3.,2., 1.], dtype = DTYPE ) / np.sqrt(14.)
        //cdef DTYPE_t Pext = self._checking_parameters[0]
        //cdef DTYPE_t aext = self._checking_parameters[2]
        value_type ni[] = { un / sqrt(14.L), deux / sqrt(14.L), trois  / sqrt(14.L) } ;
        value_type no[] = { trois / sqrt(14.L), deux / sqrt(14.L), un  / sqrt(14.L) } ;
        value_type omegaext = deuxpi / ( Pext  );
        value_type omegain = deuxpi / ( Pin ) ;    // Pin and Pext are in days
        value_type const0 = aext / pow( Pext, deux / trois ) ;
        value_type ain = const0 * pow( Pin , deux / trois ) ;
        value_type muip = 100.L ; //un ; //10. / 0.01
        value_type muop = 1000.L ; //un ; //12. / 0.001
        value_type inter = 0.L ;
        int i = 0 ;


        for (i = 0 ; i < nt ; ++i ) {
        // statevects[:, 1] = ( aext * cos( omegaext * t ) + ain * sin( omegain * t) ) / self._rscale
            //statevects[:, 10] =  ( - omegaext * aext * sin( omegaext * t ) + ain * omegain * cos(omegain * t) ) / daysec * (self._pscale * daysec)  / self._rscale
            sp[ i ][ 1 ] = ( aext * cos( omegaext * t[i] ) + ain * sin( omegain * t[i]) ) ;
            sp[ i ][ 4 ] = ( - omegaext * aext * sin( omegaext * t[i] ) + ain * omegain * cos(omegain * t[i] ) );
            //d2roeXc = ( -  pow( omegaext , 2 )  * aext * cos( omegaext * t ) - ain *  pow( omegain , 2 )  * sin(omegain * t) ) / daysec**2  * (self._pscale * daysec)**2  / self._rscale

            // Just to have something but should not play a role in einstein or roemer
            //statevects[:, 0] = aext / 2. * cos(omegaext * t) / self._rscale
            //statevects[:, 2] = aext / 2. * sin(omegaext * t) / self._rscale
            //statevects[:, 9] = - aext / 2. * omegaext / daysec * sin( omegaext * t ) * (self._pscale * daysec)  / self._rscale
            //statevects[:, 11] =  aext / 2. * omegaext / daysec * cos( omegaext * t ) * (self._pscale * daysec)  / self._rscale
            sp[ i ][ 0 ] = aext / deux * cos(omegaext * t[i] )    ;
            sp[ i ][ 2 ] = aext / deux * sin(omegaext * t[i] )  ;
            sp[ i ][ 3 ] = - aext / deux * omegaext * sin( omegaext * t[i] ) ;
            sp[ i ][ 5 ] =  aext / deux * omegaext * cos( omegaext * t[i] ) ;

        // statevects[i, 3:6] = statevects[i, 0:3] + muip * ain / ( ain / aext * cos( omegaext * t[i] ) + sin( omegain * t[i]  ) )**2 * ni / self._rscale
        // statevects[i, 6:9] = statevects[i, 0:3] + muop * aext / (  cos( omegaext * t[i] ) + ain / aext * sin( omegain * t[i]  ) )**2 * no / self._rscale
            // Just to have something but should play a role in einstein or roemer
        // statevects[i, 12:15] = ( statevects[i, 3:6] -2. * muip * ain * ( - ain / aext * omegaext * sin( omegaext * t[i] ) + omegain * cos( omegain * t[i]  ) ) /
            //                                                            ( ain / aexpowt * cos( omegaext * t[i] ) + sin( omegain * t[i]  ) )**4 / daysec * ni *
            //                                                          (self._pscale * daysec)  / self._rscale )
            //statevects[i, 15:18] = (statevects[i, 3:6] - 2. * muop * aext * ( - omegaext * sin( omegaext * t[i] ) + ain / aext * omegain * cos( omegain * t[i]  ) ) /
            //                           (  cos( omegaext * t[i] ) + ain / aext * sin( omegain * t[i]  ) )**4 / daysec * no *
            //                          (self._pscale * daysec)  / self._rscale)
            inter = muip * ain  / pow( ain / aext * cos( omegaext * t[i] ) + sin( omegain * t[i]  ) , 2)  ;
            si[i ][ 0 ] = sp[ i ][ 0 ] + inter * ni[0] ;
            si[i ][ 1] = sp[ i ][ 1 ] + inter * ni[1] ;
            si[i ][ 2] = sp[ i ][ 2 ] + inter * ni[2] ;

            inter = muop * aext / pow (  cos( omegaext * t[i] ) + ain / aext * sin( omegain * t[i]  ) , 2 ) ;
            so[i ][ 0] = sp[ i ][ 0 ] + inter * no[0] ;
            so[i ][ 1] = sp[ i ][ 1 ] + inter * no[1] ;
            so[i ][ 2] = sp[ i ][ 2 ] + inter * no[2] ;

         // Just to have something but should play a role in einstein or roemer
            inter = - deux * muip * ain * ( - ain / aext * omegaext * sin( omegaext * t[i] ) + omegain * cos( omegain * t[i]  ) ) /
                                                                        pow( ain / aext * cos( omegaext * t[i] ) + sin( omegain * t[i]  ) , 4 ) ;
            si[i ][ 3] = sp[ i ][ 3 ] + inter * ni[0] ;
            si[i ][ 4] = sp[ i ][ 4 ] + inter * ni[1] ;
            si[i ][ 5] = sp[ i ][ 5 ] + inter * ni[2] ;

            inter =- deux * muop * aext * ( - omegaext * sin( omegaext * t[i] ) + ain / aext * omegain * cos( omegain * t[i]  ) ) /
                                        pow(  cos( omegaext * t[i] ) + ain / aext * sin( omegain * t[i]  ) , 4 ) ;
            so[i ][ 3] = sp[ i ][ 3 ] + inter * no[0] ;
            so[i ][ 4] = sp[ i ][ 4 ] + inter * no[1] ;
            so[i ][ 5] = sp[ i ][ 5 ] + inter * no[2] ;

         // Putting units meter/sec in velocities
            sp[ i ][ 3 ] /= daysec ;
            sp[ i ][ 4 ] /= daysec ;
            sp[ i ][ 5 ] /= daysec ;
            si[ i ][ 3 ] /= daysec ;
            si[ i ][ 4 ] /= daysec ;
            si[ i ][ 5 ] /= daysec ;
            so[ i ][ 3 ] /= daysec ;
            so[ i ][ 4 ] /= daysec ;
            so[ i ][ 5 ] /= daysec ;

        }

        return ;
}



void  Generate_einstein_geometric_state_vectors_diagnostic(const value_type t[], int nt,
                                                        value_type Pext, value_type aext, value_type Pin, value_type ain,
                                                        value_type ra, value_type rap, value_type dec, value_type decp,
                                                        value_type ** sp, value_type ** si, value_type ** so ,
                                                        value_type * nss[3], value_type * earth_ssb[3], value_type spinaxis[3],
                                                        value_type muip = 100.L, value_type muep = 1000.L  )
{

        value_type ni[] = { un / sqrt(14.L), deux / sqrt(14.L), trois  / sqrt(14.L) } ;
        value_type no[] = { trois / sqrt(14.L), deux / sqrt(14.L), un  / sqrt(14.L) } ;
        value_type omegaext = deuxpi / ( Pext  );
        value_type omegain = deuxpi / ( Pin ) ;    // Pin and Pext are in days
        value_type omt = deuxpi / 365.25 ;
        //value_type const0 = aext / pow( Pext, deux / trois ) ;
        //value_type ain = const0 * pow( Pin , deux / trois ) ;
        value_type at = 150.e9 ; // about one astronomical unit (earth orbit)
        //value_type muip = 100.L ; //un ; //10. / 0.01
        value_type muop = muep ; //un ; //12. / 0.001
        value_type inter = 0.L ;
        int i = 0 ;


        for (i = 0 ; i < nt ; ++i ) {

            sp[ i ][ 1 ] = ( aext * cos( omegaext * t[i] ) + ain * sin( omegain * t[i]) ) ;
            sp[ i ][ 4 ] = ( - omegaext * aext * sin( omegaext * t[i] ) + ain * omegain * cos(omegain * t[i] ) );


            // Just to have something but should not play a role in einstein or roemer

            sp[ i ][ 0 ] = aext / deux * cos(omegaext * t[i] ) + ain / deux * sin( omegain * t[i])   ;
            sp[ i ][ 2 ] = aext / deux * sin(omegaext * t[i] ) +  + ain / deux * cos( omegain * t[i]) ;
            sp[ i ][ 3 ] = undemi*ain*omegain*cos(omegain*t[i]) - undemi*aext*omegaext*sin(omegaext*t[i]);//- aext / deux * omegaext * sin( omegaext * t[i] ) ;
            sp[ i ][ 5 ] = undemi*aext*omegaext*cos(omegaext*t[i]) - undemi*ain*omegain*sin(omegain*t[i]) ;// aext / deux * omegaext * cos( omegaext * t[i] ) ;

            inter = muip * ain  / pow( ain / aext * cos( omegaext * t[i] ) + sin( omegain * t[i]  ) , 2)  ;
            si[i ][ 0 ] = sp[ i ][ 0 ] + inter * ni[0] ;
            si[i ][ 1] = sp[ i ][ 1 ] + inter * ni[1] ;
            si[i ][ 2] = sp[ i ][ 2 ] + inter * ni[2] ;

            inter = muop * aext / pow (  cos( omegaext * t[i] ) + ain / aext * sin( omegain * t[i]  ) , 2 ) ;
            so[i ][ 0] = sp[ i ][ 0 ] + inter * no[0] ;
            so[i ][ 1] = sp[ i ][ 1 ] + inter * no[1] ;
            so[i ][ 2] = sp[ i ][ 2 ] + inter * no[2] ;

         // Just to have something but should play a role in einstein or roemer
            inter = - deux * muip * ain * ( - ain / aext * omegaext * sin( omegaext * t[i] ) + omegain * cos( omegain * t[i]  ) ) /
                                                                        pow( ain / aext * cos( omegaext * t[i] ) + sin( omegain * t[i]  ) , 4 ) ;
            si[i ][ 3] = sp[ i ][ 3 ] + inter * ni[0] ;
            si[i ][ 4] = sp[ i ][ 4 ] + inter * ni[1] ;
            si[i ][ 5] = sp[ i ][ 5 ] + inter * ni[2] ;

            inter =- deux * muop * aext * ( - omegaext * sin( omegaext * t[i] ) + ain / aext * omegain * cos( omegain * t[i]  ) ) /
                                        pow(  cos( omegaext * t[i] ) + ain / aext * sin( omegain * t[i]  ) , 4 ) ;
            so[i ][ 3] = sp[ i ][ 3 ] + inter * no[0] ;
            so[i ][ 4] = sp[ i ][ 4 ] + inter * no[1] ;
            so[i ][ 5] = sp[ i ][ 5 ] + inter * no[2] ;

         // Putting units meter/sec in velocities
            sp[ i ][ 3 ] /= daysec ;
            sp[ i ][ 4 ] /= daysec ;
            sp[ i ][ 5 ] /= daysec ;
            si[ i ][ 3 ] /= daysec ;
            si[ i ][ 4 ] /= daysec ;
            si[ i ][ 5 ] /= daysec ;
            so[ i ][ 3 ] /= daysec ;
            so[ i ][ 4 ] /= daysec ;
            so[ i ][ 5 ] /= daysec ;

         // Computing the nss and earth_ssb
            nss[i][0] = cos(decp*t[i] + dec)*cos(rap*t[i] + ra) ;
            nss[i][1] = cos(decp*t[i] + dec)*sin(rap*t[i] + ra) ;
            nss[i][2] = sin(decp*t[i] + dec) ;

            earth_ssb[i][0] = at*cos(omt*t[i] ) ;
            earth_ssb[i][1] = zero ;
            earth_ssb[i][2] = at*sin( omt*t[i] ) ;

        }

        spinaxis[0] = zero;
        spinaxis[1] = un;
        spinaxis[2] = zero;

        return ;
}

void Compute_geometric_diagnostic(const value_type ts[], const value_type t0 , value_type geometric[], int nt,
                                     value_type Pext, value_type Pi, value_type ae, value_type ai,
                                     value_type ra, value_type rap, value_type dec, value_type decp,
                                     value_type distance, value_type distance1)
// TO use with Generate_einstein_geometric_state_vectors_diagnostic
/* ts, t0 in julian days
 * Pext, Pi in julian daysfcn.fitwithminuit[0].
 * ae, ai in meters
 * ra , dec in radians
 * rap, decp in radians per days
 * distance in light years
 * distance1 in km/s
 * geometric is returned in days
 */
{
    value_type t = 0.;
    value_type omt = deuxpi / 365.25;
    value_type at = 150.e9 ; // about one astronomical unit (earth orbit)
    value_type omi = deuxpi / Pi;
    value_type ome = deuxpi/ Pext;
    value_type d = distance * ( clight * daysec * 365.25 ) ; // ltyr -> meters
    value_type dprime = 1000.*distance1 * daysec ; // km/s -> meters / day
//     value_type clightm_d = clight * daysec ; // speed of light in meters/ day
    value_type geo0, geo1, geo2, geo3;

    for (long int i = 0 ; i < nt ; i++)
    {
        t = ts[i];
        // this formula is also ok i think, but  unpractical to compare term by term
//         geometric[i] = pow(dprime,deux)*pow(t,deux)/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) -
//                         2.*pow(dprime,deux)*t*t0/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) +
//                         pow(dprime,deux)*pow(t0,deux)/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) +
//                         undemi*( (ae*cos(ome*t) - 2.*at*cos(omt*t) + ai*sin(omi*t))*cos(decp*t + dec)*cos(rap*t + ra) + 2*(ae*cos(ome*t) + ai*sin(omi*t))*cos(decp*t + dec)*sin(rap*t + ra) +
//                         (ai*cos(omi*t) + ae*sin(ome*t) - 2.*at*sin(omt*t))*sin(decp*t + dec) ) * dprime*t/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) +
//                         d*dprime*t/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) -
//                         undemi*( (ae*cos(ome*t) - 2.*at*cos(omt*t) + ai*sin(omi*t))*cos(decp*t + dec)*cos(rap*t + ra) + 2*(ae*cos(ome*t) +
//                         ai*sin(omi*t) ) * cos(decp*t + dec)*sin(rap*t + ra) + (ai*cos(omi*t) + ae*sin(ome*t) -
//                         2.*at*sin(omt*t))*sin(decp*t + dec) ) *dprime*t0/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) -
//                         d*dprime*t0/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) -
//                         1./8.*pow((ae*cos(ome*t) - 2*at*cos(omt*t) + ai*sin(omi*t))*cos(decp*t + dec)*cos(rap*t + ra) + 2*(ae*cos(ome*t) +
//                         ai*sin(omi*t))*cos(decp*t + dec)*sin(rap*t + ra) + (ai*cos(omi*t) + ae*sin(ome*t) -
//                         2.*at*sin(omt*t) ) * sin(decp*t + dec),deux) / (clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) +
//                         undemi*((ae*cos(ome*t) - 2.*at*cos(omt*t) + ai*sin(omi*t))*cos(decp*t + dec)*cos(rap*t + ra) + 2.*(ae*cos(ome*t) +
//                         ai*sin(omi*t))*cos(decp*t + dec)*sin(rap*t + ra) +
//                         (ai*cos(omi*t) + ae*sin(ome*t) - 2.*at*sin(omt*t) ) *sin(decp*t + dec))*d/(clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d) +
//                         1./8.*( pow( ae*cos(ome*t) - 2*at*cos(omt*t) + ai*sin(omi*t),deux) + 4.*pow(ae*cos(ome*t) + ai*sin(omi*t),deux) + pow( ai*cos(omi*t) + ae*sin(ome*t) - 2.*at*sin(omt*t),deux) ) /
//                         (clightm_d*dprime*t - clightm_d*dprime*t0 + clightm_d*d)   ;

       geo0 =  undemi *((ae*cos(ome*t) - 2.*at*cos(omt*t) + ai*sin(omi*t))*cos(decp*t + dec)*cos(rap*t + ra) +
                2.*(ae*cos(ome*t) + ai*sin(omi*t))*cos(decp*t + dec)*sin(rap*t + ra) + (ai*cos(omi*t) + ae*sin(ome*t) -
                2.*at*sin(omt*t))*sin(decp*t + dec)) / clight ;

       geo1 = 0.1250*(pow(ae*cos(ome*t) - 2.*at*cos(omt*t) + ai*sin(omi*t),2) + 4.*pow(ae*cos(ome*t) + ai*sin(omi*t),2) +
                pow(ai*cos(omi*t) + ae*sin(ome*t) - 2.*at*sin(omt*t),2))/((dprime*(t - t0) + d)*clight) ;

       geo2 = -0.125*pow((ae*cos(ome*t) - 2.*at*cos(omt*t) + ai*sin(omi*t))*cos(decp*t + dec)*cos(rap*t + ra) + 2.*(ae*cos(ome*t) +
                ai*sin(omi*t))*cos(decp*t + dec)*sin(rap*t + ra) + (ai*cos(omi*t) + ae*sin(ome*t) -
                2.*at*sin(omt*t))*sin(decp*t + dec),2) / ((dprime*(t - t0) + d)*clight) ;

       geo3 = dprime*(t - t0)/clight ;

       //if (i < 10) printf(" geo %i  %.19Le %.19Le  %.19Le  %.19Le \n", i, geo0/daysec, geo1/daysec, geo2/daysec, geo3/daysec);
       geometric[i] = ( geo0 + geo1 + geo2 + geo3 ) / daysec ;

    }
    return;
}



void Compute_einstein_diagnostic_for_geometric(const value_type ts[], const value_type t0 , value_type einstein[], int nt,
                                     value_type Pext, value_type Pi, value_type ae, value_type ai, value_type mi, value_type me,
                                     value_type ra, value_type rap, value_type dec, value_type decp,
                                              value_type muip = 100.L, value_type muep = 1000.L  )
// TO use with Generate_einstein_geometric_state_vectors_diagnostic
/* ts, t0 in julian days
 * Pext, Pi in julian daysaberration
 * mi and me are in solar masses
 * ae, ai in meters
 * ra , dec in radians
 * rap, decp in radians per dayss
 * einstein is returned in days
 */
{
    long int i = 0L;
    value_type t = 0.;
    value_type omi = deuxpi / Pi;
    value_type ome = deuxpi/ Pext;
    value_type einui, einue, einvp;
    value_type GMsol = Ggrav * Msol;
    value_type ae2 = pow(ae,2);
    value_type ai2 = pow(ai,2);
    value_type ae3 = pow(ae,3);

    for (i = 0 ; i < nt ; i ++)
        {
        t = ts[i] ;
        einui = -unquart*(4.*GMsol*ae*ai*mi*(cos((ome + omi)*t)/(ome + omi) - cos(-(ome - omi)*t)/(ome - omi)) -
        (2.*ome*t + sin(2*ome*t))*GMsol*ai2*mi/ome - (2.*omi*t - sin(2.*omi*t))*GMsol*ae2*mi/omi) / (ae2*ai* pow(clight,2) * muip ) ;
        t = t0 ;
        einui -=  -unquart*(4.*GMsol*ae*ai*mi*(cos((ome + omi)*t)/(ome + omi) - cos(-(ome - omi)*t)/(ome - omi)) -
        (2.*ome*t + sin(2*ome*t))*GMsol*ai2*mi/ome - (2.*omi*t - sin(2.*omi*t))*GMsol*ae2*mi/omi) / (ae2*ai* pow(clight,2) * muip ) ;

        t = ts[i] ;
        einue = -unquart*(4.*GMsol*ae*ai*me*(cos((ome + omi)*t)/(ome + omi) - cos(-(ome - omi)*t)/(ome - omi)) -
                (2.*ome*t + sin(2.*ome*t))*GMsol*ae2*me/ome - (2.*omi*t - sin(2.*omi*t))*GMsol*ai2*me/omi)/(ae3*pow(clight,2) * muep) ;
        t = t0 ;
        einue -= -unquart*(4.*GMsol*ae*ai*me*(cos((ome + omi)*t)/(ome + omi) - cos(-(ome - omi)*t)/(ome - omi)) -
                (2.*ome*t + sin(2.*ome*t))*GMsol*ae2*me/ome - (2.*omi*t - sin(2.*omi*t))*GMsol*ai2*me/omi)/(ae3*pow(clight,2) * muep) ;

        t = ts[i] ;
        einvp = 1./16.L*(20.L*ae*ai*ome*omi*(cos((ome + omi)*t)/(ome + omi) + cos((ome - omi)*t)/(ome - omi)) +
                4.*ae*ai*ome*omi*(cos((ome + omi)*t)/(ome + omi) - cos(-(ome - omi)*t)/(ome - omi)) +
                (2.*ome*t + sin(2.*ome*t))*ae2*ome + 5.*(2.*ome*t - sin(2.*ome*t))*ae2*ome +
                5.*(2.*omi*t + sin(2.L*omi*t))*ai2*omi + (2.L*omi*t - sin(2.L*omi*t))*ai2*omi)/pow(clight*daysec,2) ;

        t = t0 ;
        einvp -= 1./16.L*(20.L*ae*ai*ome*omi*(cos((ome + omi)*t)/(ome + omi) + cos((ome - omi)*t)/(ome - omi)) +
                4.*ae*ai*ome*omi*(cos((ome + omi)*t)/(ome + omi) - cos(-(ome - omi)*t)/(ome - omi)) +
                (2.*ome*t + sin(2.*ome*t))*ae2*ome + 5.*(2.*ome*t - sin(2.*ome*t))*ae2*ome +
                5.*(2.*omi*t + sin(2.L*omi*t))*ai2*omi + (2.L*omi*t - sin(2.L*omi*t))*ai2*omi)/pow(clight*daysec,2) ;

        einstein[i]  = einui + einue + undemi * einvp ;
    }
    return;

}


void Compute_aberration_diagnostic(const value_type ts[], value_type aberration[], int nt,
                                     value_type spinfreq, value_type Pext, value_type Pi, value_type ae, value_type ai,
                                     value_type ra, value_type rap, value_type dec, value_type decp)
// TO use with Generate_einstein_geometric_state_vectors_diagnostic
/* ts, t0 in julian days
 * spinfreq in days^-1
 * Pext, Pi in julian days
 * ae, ai in meters
 * ra , dec in radians
 * rap, decp in radians per dayss
 * aberration is returned in days
 */
{
    value_type t = 0.;
    value_type om = deuxpi * spinfreq;
    value_type omi = deuxpi / Pi;
    value_type ome = deuxpi/ Pext;
    value_type clightm_d = clight * daysec ; // speed of light in meters/ day

    for (long int i = 0 ; i < nt ; i++)
    {
        t = ts[i];
        aberration[i] = -undemi* ( (ae*ome*cos(ome*t) - ai*omi*sin(omi*t))*om*cos(decp*t + dec)*cos(rap*t + ra) -
                        (ai*omi*cos(omi*t) - ae*ome*sin(ome*t))*om*sin(decp*t + dec)) /
                        ((pow(om*cos(decp*t + dec)*cos(rap*t + ra),deux) + pow(om*sin(decp*t + dec),deux))*clightm_d);
        //aberration[i] /= daysec; spinfreq en jours^-1
    }
    return;
}



void Compute_shapiro_diagnostic(const value_type ts[], value_type shapiro[], int nt,
                                     value_type mi, value_type me, value_type Pext, value_type Pi, value_type ae, value_type ai,
                                     value_type ra, value_type rap, value_type dec, value_type decp,
                                     value_type muip = 100.L, value_type muep = 1000.L)
// TO use with Generate_einstein_geometric_state_vectors_diagnostic
/* ts, t0 in julian days
 * Pext, Pi in julian days
 * ae, ai in meters
 * mi and me in solar masses
 * ra , dec in radians
 * rap, decp in radians per days
 * shapiro is returned in days
 */
{
    value_type t = 0.;
    value_type omi = deuxpi / Pi;
    value_type ome = deuxpi/ Pext;
    value_type inter = zero;

    for (long int i = 0 ; i < nt ; i++)
    {
        t = ts[i];
        inter  = 1./14.*sqrt(14.)*ai*muip*cos(decp*t + dec)*cos(rap*t + ra)/pow(ai*cos(ome*t)/ae +
                 sin(omi*t),deux) + 1/7.*sqrt(14.)*ai*muip*cos(decp*t + dec)*sin(rap*t + ra)/pow(ai*cos(ome*t)/ae + sin(omi*t),deux) +
                 3./14.*sqrt(14.)*ai*muip*sin(decp*t + dec)/pow(ai*cos(ome*t)/ae + sin(omi*t),deux) +
                 ai*muip/pow(ai*cos(ome*t)/ae + sin(omi*t),deux);
        shapiro[i] = - deux * mi * Ggrav * Msol / pow(clight,3) * log(inter / clight) ;
        inter = 3./14.*sqrt(14.)*ae*muep*cos(decp*t + dec)*cos(rap*t + ra)/pow(ai*sin(omi*t)/ae + cos(ome*t),deux) +
                1./7.*sqrt(14.)*ae*muep*cos(decp*t + dec)*sin(rap*t + ra)/pow(ai*sin(omi*t)/ae + cos(ome*t),deux) +
                1./14.*sqrt(14)*ae*muep*sin(decp*t + dec)/pow(ai*sin(omi*t)/ae + cos(ome*t),deux) +
                ae*muep/pow(ai*sin(omi*t)/ae + cos(ome*t),deux)  ;
        shapiro[i] += - deux * me * Ggrav * Msol / pow(clight,3) * log(inter / clight) ;
        shapiro[i] /= daysec ;
    }
    return;
}



void Compute_roemer_diagnostic(const value_type t[], value_type roemer[], int nt,
                               value_type Pext, value_type aext, value_type Pin) {
/*
             Identical to _Generate_check_roemer_state_vectors_ but return only the roemer delay (in days) at (not retarded) time t.
 */
        int i ;
        value_type omegaext = deuxpi / ( Pext );
        value_type omegain = deuxpi / ( Pin ) ;    // Pin and Pext are in days
        value_type const0 = aext / pow( Pext, deux / trois ) ;
        value_type ain = const0 * pow( Pin , deux / trois ) ;

        for(i = 0 ; i < nt ; ++i) {
            roemer[i] = aext * cos(omegaext * t[i])    +  ain * sin( omegain * t[i]) ;
            roemer[i] /= (clight * daysec ) ;
        }

        return ;
}


void Compute_einstein_diagnostic(const value_type ta[], value_type einstein[], int nt,
                               value_type Pext, value_type aext, value_type Pin,
                               value_type Mi, value_type Mo,
                               value_type t0)
{
        int i ;
        value_type omegaext = deuxpi / ( Pext *  daysec);
        value_type omegain = deuxpi / ( Pin * daysec) ;    // Pin and Pext are in days
        value_type const0 = aext / pow( Pext, deux / trois ) ;
        value_type ain = const0 * pow( Pin , deux / trois ) ;
        value_type muip = 100.L ;   // un ; //10. / 0.01
        value_type muop = 1000.L ;  // un ; //12. / 0.001
        value_type GMsol = Ggrav * Msol ;
        value_type DeltaEi, DeltaEo, DeltaEv ;
        value_type tref = t0 * daysec ;
//         value_type t[nt] ;
        value_type tm ;

        for (i = 0 ; i < nt ; ++i ) {
            tm = ta[i] * daysec;

            DeltaEi = (GMsol*Mi*(2*tm + (2* pow( ain,2)*tm)/ pow( aext,2) - 2*tref -
       (2* pow( ain,2)*tref)/ pow( aext,2) +
       (8*ain*omegain* cos( omegaext*tm)* cos( omegain*tm))/
        (aext* pow( omegaext,2) - aext* pow( omegain,2)) +
       (8*ain*omegain* cos( omegaext*tref)* cos( omegain*tref))/
        (-(aext* pow( omegaext,2)) + aext* pow( omegain,2)) +
       ( pow( ain,2)* sin( 2*omegaext*tm))/( pow( aext,2)*omegaext) +
       (8*ain*omegaext* sin( omegaext*tm)* sin( omegain*tm))/
        (aext* pow( omegaext,2) - aext* pow( omegain,2)) -
        sin( 2*omegain*tm)/omegain -
       ( pow( ain,2)* sin( 2*omegaext*tref))/( pow( aext,2)*omegaext) -
       (8*ain*omegaext* sin( omegaext*tref)* sin( omegain*tref))/
        (aext* pow( omegaext,2) - aext* pow( omegain,2)) +
        sin( 2*omegain*tref)/omegain))/(4.*ain*muip) ;

//             DeltaEi = ( (GMsol*Mi*(2*t[i] + (2*pow( ain, 2 ) * t[i] ) / pow( aext, 2 ) - 2*tref -
//             (2* pow( ain, 2) *tref)/ pow( aext, 2 ) +
//             (8*ain*omegain* cos( omegaext*t[i])* cos( omegain*t[i]))/
//                 (aext*pow( omegaext, 2 ) - aext * pow( omegain , 2 )) +
//             (8*ain*omegain* cos( omegaext*tref)* cos( omegain*tref))/
//                 (-(aext* pow( omegaext ,2 ) ) + aext* pow( omegain , 2 ) ) +
//             (pow( ain , 2 )* sin( 2*omegaext*t[i]))/(pow( aext , 2 ) * omegaext) +
//             (8*ain*omegaext* sin( omegaext*t[i])* sin( omegain*t[i]))/
//                 (aext*pow( omegaext , 2 ) - aext* pow( omegain , 2 ) ) -
//                 sin( 2*omegain*t[i])/omegain -
//             ( pow( ain , 2 ) * sin( 2*omegaext*tref))/( pow( aext , 2 ) *omegaext) -
//             (8*ain*omegaext* sin( omegaext*tref)* sin( omegain*tref))/
//                 (aext* pow( omegaext , 2 )  - aext* pow( omegain , 2 ) ) +
//                 sin( 2*omegain*tref)/omegain))/(quatre*ain*muip) ) ;

       DeltaEo = (GMsol*Mo*(2* pow( aext,2)*tm + 2* pow( ain,2)*tm - 2* pow( aext,2)*tref -
       2* pow( ain,2)*tref + (8*aext*ain*omegain* cos( omegaext*tm)*
           cos( omegain*tm))/( pow( omegaext,2) -  pow( omegain,2)) +
       (8*aext*ain*omegain* cos( omegaext*tref)* cos( omegain*tref))/
        (- pow( omegaext,2) +  pow( omegain,2)) +
       ( pow( aext,2)* sin( 2*omegaext*tm))/omegaext +
       (8*aext*ain*omegaext* sin( omegaext*tm)* sin( omegain*tm))/
        ( pow( omegaext,2) -  pow( omegain,2)) -
       ( pow( ain,2)* sin( 2*omegain*tm))/omegain -
       ( pow( aext,2)* sin( 2*omegaext*tref))/omegaext -
       (8*aext*ain*omegaext* sin( omegaext*tref)* sin( omegain*tref))/
        ( pow( omegaext,2) -  pow( omegain,2)) +
       ( pow( ain,2)* sin( 2*omegain*tref))/omegain))/(4.* pow( aext,3)*muop) ;

//             DeltaEo = ( (GMsol*Mo*(2* pow( aext , 2 ) *t[i] + 2* pow( ain , 2 ) *t[i] - 2* pow( aext , 2 ) *tref -
//             2* pow( ain , 2 ) *tref + (8*aext*ain*omegain* cos( omegaext*t[i])*
//                 cos( omegain*t[i]))/( pow( omegaext , 2 )  -  pow( omegain , 2 ) ) +
//             (8*aext*ain*omegain* cos( omegaext*tref)* cos( omegain*tref))/
//                 (- pow( omegaext , 2 )  +  pow( omegain , 2 ) ) +
//             ( pow( aext , 2 ) * sin( 2*omegaext*t[i]))/omegaext +
//             (8*aext*ain*omegaext* sin( omegaext*t[i])* sin( omegain*t[i]))/
//                 ( pow( omegaext , 2 )  -  pow( omegain , 2 ) ) -
//             ( pow( ain , 2 ) * sin( 2*omegain*t[i]))/omegain -
//             ( pow( aext , 2 ) * sin( 2*omegaext*tref))/omegaext -
//             (8*aext*ain*omegaext* sin( omegaext*tref)* sin( omegain*tref))/
//                 ( pow( omegaext , 2 )  -  pow( omegain , 2 ) ) +
//             ( pow( ain , 2 ) * sin( 2*omegain*tref))/omegain))/(quatre* pow( aext , 3 ) *muop) ) ;
//
//            DeltaEv = 0.L;

       DeltaEv =0.03125* pow( aext,2)*omegaext*
    (2*omegaext*(tm - tref) +  sin( 2*omegaext*tm) -  sin( 2*omegaext*tref)) +
   0.15625* pow( aext,2)*omegaext*
    (2*omegaext*(tm - tref) -  sin( 2*omegaext*tm) +  sin( 2*omegaext*tref)) +
   (1.*aext*ain*omegaext*omegain*
      (omegaext* cos( omegaext*tm)* cos( omegain*tm) -
        1.*omegaext* cos( omegaext*tref)* cos( omegain*tref) +
        omegain* sin( omegaext*tm)* sin( omegain*tm) -
        1.*omegain* sin( omegaext*tref)* sin( omegain*tref)))/
    ( pow( omegaext,2) - 1.* pow( omegain,2)) +
   0.125* pow( ain,2)*omegain*(2*omegain*(tm - tref) +  sin( 2*omegain*tm) -
       sin( 2*omegain*tref)) ;

//             DeltaEv = ( 0.03125L* pow( aext , 2 )  * omegaext *
//             (2 * omegaext * (t[i] - tref) +  sin( 2*omegaext*t[i]) -  sin( 2*omegaext*tref))
//             + 0.15625L *  pow( aext , 2 )  * omegaext *
//             (2*omegaext*(t[i] - tref) -  sin( 2*omegaext*t[i]) +  sin( 2*omegaext*tref))
//             + (un  * aext*ain*omegaext*omegain*
//             (omegaext* cos( omegaext*t[i])* cos( omegain*t[i]) -
//                 un *omegaext* cos( omegaext*tref)* cos( omegain*tref) +
//                 omegain* sin( omegaext*t[i])* sin( omegain*t[i]) -
//                 un * omegain* sin( omegaext*tref)* sin( omegain*tref)))/
//             ( pow( omegaext , 2 )  - un * pow( omegain , 2 ) ) +
//         0.125L *  pow( ain , 2 )  * omegain*(2*omegain*(t[i] - tref) +  sin( 2 * omegain * t[i]) -
//             sin( 2*omegain*tref)) ) ;

           // einstein[i] = ( DeltaEv ) / pow( clight , 2 ) / daysec ;

          einstein[i] = ( DeltaEi + DeltaEo + DeltaEv ) / pow( clight , 2 ) / daysec ; //pow( clight * daysec , 2 ) ;
        }

        return ;
}







void IntegralePrems3_Newt(const value_type sv0[6], const value_type sv1[6], const value_type sv2[6],
                         value_type m0, value_type m1, value_type m2,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
{
            valarray<value_type> s0(sv0, 6);
        valarray<value_type> s1(sv1, 6);
        valarray<value_type> s2(sv2, 6);

        valarray<valarray<value_type> > xs(3) ;

            xs[0] = valarray<value_type>( s0[slice(0,3,1)] );
            xs[1] = valarray<value_type>( s1[slice(0,3,1)] );
            xs[2] = valarray<value_type>( s2[slice(0,3,1)] );


        valarray<valarray<value_type> > v(3) ;
            v[0] = valarray<value_type>( s0[slice(3,3,1)] );
            v[1] = valarray<value_type>( s1[slice(3,3,1)] );
            v[2] = valarray<value_type>( s2[slice(3,3,1)] );

        valarray<value_type> vv(3) ;
            vv[0] = sumsquares3d<value_type>( v[0] );
            vv[1] = sumsquares3d<value_type>( v[1] );
            vv[2] = sumsquares3d<value_type>( v[2] );

        valarray<value_type> ms(3) ;
            ms[0] = m0 ;
            ms[1] = m1 ;
            ms[2] = m2 ;

        valarray<valarray<value_type> > r(3) ;
            r[0] = valarray<value_type>( 3 ) ;
                r[0][0] = un ;
                r[0][1] = norm3d<value_type>( xs[1] - xs[0] ) ;
                r[0][2] = norm3d<value_type>( xs[2] - xs[0] ) ;
            r[1] = valarray<value_type>(3);
                r[1][0] = norm3d<value_type>(xs[0] - xs[1]) ;
                r[1][1] = un ;
                r[1][2] = norm3d<value_type>(xs[2] - xs[1] ) ;
            r[2] = valarray<value_type>(3);
                r[2][0] = norm3d<value_type>(xs[0] - xs[2]) ;
                r[2][1] = norm3d<value_type>(xs[1] - xs[2]) ;
                r[2][2] = un ;

        value_type Mt = m0 + m1 + m2 ;
        value_type GMsol = Ggrav * Msol ;

        int n1 = 0;
        int n2 = 1;
        int n3 = 2;

        impulsion =  ms[0]*v[0] + ms[1]*v[1] + ms[2]*v[2] ;

        energy = (ms[n1]*vv[n1])/ deux  + (ms[n2]*vv[n2])/ deux  + (ms[n3]*vv[n3])/ deux -
                    ( GMsol*ms[n1]*ms[n2] / r[n1][n2] + GMsol*ms[n2]*ms[n3] / r[n2][n3] + GMsol*ms[n1]*ms[n3] / r[n1][n3] ) ;

        center_of_mass =(ms[n1]*xs[n1] + ms[n2]*xs[n2] + ms[n3]*xs[n3] ) / Mt;

        return;


}


void IntegralePrems2_NewtQuad(const value_type sv0[6], const value_type sv1[6],
                         value_type m0, value_type m1, value_type q1,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
// q1 is the quadrupole moment along the 01 axis in units of Msol.ls^2
{
        valarray<value_type> s0(sv0, 6);
        valarray<value_type> s1(sv1, 6);

        valarray<valarray<value_type> > xs(2) ;

            xs[0] = valarray<value_type>( s0[slice(0,3,1)] );
            xs[1] = valarray<value_type>( s1[slice(0,3,1)] );
//             xs[2] = valarray<value_type>( s2[slice(0,3,1)] );


        valarray<valarray<value_type> > v(2) ;
            v[0] = valarray<value_type>( s0[slice(3,3,1)] );
            v[1] = valarray<value_type>( s1[slice(3,3,1)] );
//             v[2] = valarray<value_type>( s2[slice(3,3,1)] );

        valarray<value_type> vv(2) ;
            vv[0] = sumsquares3d<value_type>( v[0] );
            vv[1] = sumsquares3d<value_type>( v[1] );
//             vv[2] = sumsquares3d<value_type>( v[2] );

        valarray<value_type> ms(2) ;
            ms[0] = m0 ;
            ms[1] = m1 ;
//             ms[2] = m2 ;

        valarray<valarray<value_type> > r(2) ;
            r[0] = valarray<value_type>( 2 ) ;
                r[0][0] = un ;
                r[0][1] = norm3d<value_type>( xs[1] - xs[0] ) ;
//                 r[0][2] = norm3d<value_type>( xs[2] - xs[0] ) ;
            r[1] = valarray<value_type>(2);
                r[1][0] = norm3d<value_type>(xs[0] - xs[1]) ;
                r[1][1] = un ;
//                 r[1][2] = norm3d<value_type>(xs[2] - xs[1] ) ;
//             r[2] = valarray<value_type>(3);
//                 r[2][0] = norm3d<value_type>(xs[0] - xs[2]) ;
//                 r[2][1] = norm3d<value_type>(xs[1] - xs[2]) ;
//                 r[2][2] = un ;

        value_type Mt = m0 + m1;// + m2 ;
        value_type GMsol = Ggrav * Msol ;

        int n1 = 0;
        int n2 = 1;
//         int n3 = 2;

        impulsion =  ms[0]*v[0] + ms[1]*v[1] ;//+ ms[2]*v[2] ;
        energy = (ms[n1]*vv[n1])/ deux  + (ms[n2]*vv[n2])/ deux - ( GMsol*ms[n1]*ms[n2] / r[n1][n2] + troisdemis * GMsol * m0 *(q1*clight*clight) / (r[n1][n2] * r[n1][n2] * r[n1][n2]) ) ;

        center_of_mass =(ms[n1]*xs[n1] + ms[n2]*xs[n2] ) / Mt;

        return;


}



void IntegralePrems3_1PN(const value_type sv0[6], const value_type sv1[6], const value_type sv2[6],
                         value_type m0, value_type m1, value_type m2, value_type quadrupole1,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
{
        /*
            R, L, P, E = IntegralePrems3(sv0,sv1,sv2,m0,m1,m2) \n
            sv0, sv1, sv2 are the state vectors (postion and velocity) of bodies of masses m0, m1, m2 \n
            R is the center of mass, L the angular momentum, P the impulsion and E the energy of the system.
            To compute the first integrals for a serie of n times, sv0 must have the shape (n times, 6 components), then the first integrals will have the same shape.
            If sv0.shape = (6), then the vectorial first integrals will also have this shape.
        */

//         R = (m0*s0[:,0:3] + m1*s1[:,0:3] + m2*s2[:,0:3])/(m0+m1+m2)
//             #if self._Get_theory_() == 1 :
//         R_GR = (m0 * np.linalg.norm(s0[:,3:6])**2/c**2 * s0[:,0:3] +        # the grav potential term is missing but should also be negligible
//                 m1 * np.linalg.norm(s1[:,3:6])**2/c**2 * s1[:,0:3] +
//                 m2 * np.linalg.norm(s2[:,3:6])**2/c**2 * s2[:,0:3] ) / (m0+m1+m2)
//         print "rgr : ", np.max(R_GR)
//         GR_cor = np.linalg.norm(R_GR, axis = 1)  / np.linalg.norm(R, axis=1)
//         print "grcor ", np.max(GR_cor)
//         normGR  =  (m0+m1+m2) / ( (m0+m1+m2)  + m0 * np.linalg.norm(s0[:,3:6])**2/c**2  +       # the grav potential term is missing but should also be negligible
//                     m1 * np.linalg.norm(s1[:,3:6])**2/c**2 +
//                     m2 * np.linalg.norm(s2[:,3:6])**2/c**2 )
//         #R = (R + R_GR) * normGR
//         print 'GR not activated in IntegralePrems3 !! '
//         print "Average and maximum of GR relative correction on barycenter : ", np.mean(GR_cor), np.max(GR_cor)
//         L = m0*np.cross(s0[:,0:3],s0[:,3:6]) + m1*np.cross(s1[:,0:3],s1[:,3:6]) + m2*np.cross(s2[:,0:3],s2[:,3:6])
//         P = m0*s0[:,3:6] + m1*s1[:,3:6] + m2*s2[:,3:6]
//         E = np.ones(len(s0[:,0]))

        valarray<value_type> s0(sv0, 6);
        valarray<value_type> s1(sv1, 6);
        valarray<value_type> s2(sv2, 6);

        valarray<valarray<value_type> > xs(3) ;

            xs[0] = valarray<value_type>( s0[slice(0,3,1)] );
            xs[1] = valarray<value_type>( s1[slice(0,3,1)] );
            xs[2] = valarray<value_type>( s2[slice(0,3,1)] );


        valarray<valarray<value_type> > v(3) ;
            v[0] = valarray<value_type>( s0[slice(3,3,1)] );
            v[1] = valarray<value_type>( s1[slice(3,3,1)] );
            v[2] = valarray<value_type>( s2[slice(3,3,1)] );

        valarray<value_type> vv(3) ;
            vv[0] = sumsquares3d<value_type>( v[0] );
            vv[1] = sumsquares3d<value_type>( v[1] );
            vv[2] = sumsquares3d<value_type>( v[2] );

        valarray<value_type> ms(3) ;
            ms[0] = m0 ;
            ms[1] = m1 ;
            ms[2] = m2 ;

        value_type Mt = m0 + m1 + m2 ;
        value_type GMsol = Ggrav * Msol ;

        valarray<value_type> VV0( ( ms[0]*v[0] + ms[1]*v[1] + ms[2]*v[2] ) / Mt );

        /////V = np.ndarray((len(s0[:,0]), 3), dtype = np.float)

//                vv = np.array([ np.linalg.norm(s0[i, 3:6])**2 , np.linalg.norm(s1[i, 3:6])**2 ,  np.linalg.norm(s2[i, 3:6])**2 ])
//                 v = np.array([s0[i, 3:6],s1[i, 3:6], s2[i, 3:6] ])
                //xs = np.array([s0[i, 3:6] , s1[i, 3:6] ,  s2[i, 3:6] ])
//                 #np.array([[np.zeros(3)*1., s1[i, 0:3] - s0[i, 0:3], s2[i, 0:3] - s0[i, 0:3]],
//                               #[s0[i, 0:3] - s1[i, 0:3], np.zeros(3)*1. , s2[i, 0:3] - s1[i, 0:3]],
//                               #[s0[i, 0:3] - s2[i, 0:3], s1[i, 0:3] - s2[i, 0:3], np.zeros(3)*1.]])

        valarray<valarray<value_type> > r(3) ;
            r[0] = valarray<value_type>( 3 ) ;
                r[0][0] = un ;
                r[0][1] = norm3d<value_type>( xs[1] - xs[0] ) ;
                r[0][2] = norm3d<value_type>( xs[2] - xs[0] ) ;
            r[1] = valarray<value_type>(3);
                r[1][0] = norm3d<value_type>(xs[0] - xs[1]) ;
                r[1][1] = un ;
                r[1][2] = norm3d<value_type>(xs[2] - xs[1] ) ;
            r[2] = valarray<value_type>(3);
                r[2][0] = norm3d<value_type>(xs[0] - xs[2]) ;
                r[2][1] = norm3d<value_type>(xs[1] - xs[2]) ;
                r[2][2] = un ;

        valarray<valarray<valarray<value_type> > > ns(3) ;
            ns[0] = valarray<valarray<value_type> >(3);
                ns[0][0] = valarray<value_type>(3);
                ns[0][1] = valarray<value_type>( xs[1] - xs[0] ) ;
                ns[0][2] = valarray<value_type>( xs[2] - xs[0] ) ;
            ns[1] = valarray<valarray<value_type> >(3);
                ns[1][0] = valarray<value_type>( xs[0] - xs[1] ) ;
                ns[1][1] = valarray<value_type>(3);
                ns[1][2] = valarray<value_type>( xs[2] - xs[1] ) ;
            ns[2] = valarray<valarray<value_type> >(3);
                ns[2][0] = valarray<value_type>( xs[0] - xs[2] ) ;//valarray<value_type>( xs[2] - xs[0] ) ; // !!!!!!!!!!!!!!!!!!
                ns[2][1] = valarray<value_type>( xs[1] - xs[2] );
                ns[2][2] = valarray<value_type>(3);
            for (int j = 0 ; j<3 ; ++j)
            {
                for (int k = 0 ; k<3 ; ++k) ns[j][k] /= r[j][k];
            }

        int n1 = 0 ;
        int n2 = 1 ;
        int n3 = 2 ;
        const value_type cl = clight ;
        const value_type clm2 = pow(cl,-2) ;



// This expression is derived from the Lagrangian in Soffel's book "Relativity in Astrometry, Celestial Mechanics and Geodesy" , also in Damour And Taylor 1992, Eq_motion_1PN-Will_TEGP.nb.
// Agrees with the formula given by Will, 2014, "Incorporating post-Newtonian effects in N-body dynamics", formula 3.2
// !!! Don't forget to remove the mass energy from it !!!
//
        energy = -(pow(cl,-2)*(-4* dotprod<value_type>( v[n3],v[n3])*ms[n3]*vv[n3] - 8* dotprod<value_type>( v[n3],v[n3])*ms[n3]*pow(cl,2) + 4*ms[n1]*vv[n1]*pow(cl,2) + 4*ms[n2]*vv[n2]*pow(cl,2) +
        4*ms[n3]*vv[n3]*pow(cl,2) + ms[n1]*pow(vv[n1],2) + ms[n2]*pow(vv[n2],2) + ms[n3]*pow(vv[n3],2) - 4*ms[n1]*pow(GMsol,2)*pow(ms[n2],2)*pow(r[n1][n2],-2) +
        14*GMsol* dotprod<value_type>( v[n2],v[n1])*ms[n1]*ms[n2]*pow(r[n1][n2],-1) +
        2*GMsol* dotprod<value_type>( v[n1],ns[n1][n2])* dotprod<value_type>( v[n2],ns[n1][n2])*ms[n1]*ms[n2]*pow(r[n1][n2],-1) + 6*GMsol*ms[n1]*ms[n2]*vv[n1]*pow(r[n1][n2],-1) +
        6*GMsol*ms[n1]*ms[n2]*vv[n2]*pow(r[n1][n2],-1) + 4*GMsol*ms[n1]*ms[n2]*pow(cl,2)*pow(r[n1][n2],-1) -
        4*ms[n1]*pow(GMsol,2)*pow(ms[n3],2)*pow(r[n1][n3],-2) + 14*GMsol* dotprod<value_type>( v[n3],v[n1])*ms[n1]*ms[n3]*pow(r[n1][n3],-1) -
        12*GMsol* dotprod<value_type>( v[n3],v[n3])*ms[n1]*ms[n3]*pow(r[n1][n3],-1) +
        2*GMsol* dotprod<value_type>( v[n1],ns[n1][n3])* dotprod<value_type>( v[n3],ns[n1][n3])*ms[n1]*ms[n3]*pow(r[n1][n3],-1) + 6*GMsol*ms[n1]*ms[n3]*vv[n1]*pow(r[n1][n3],-1) +
        6*GMsol*ms[n1]*ms[n3]*vv[n3]*pow(r[n1][n3],-1) + 4*GMsol*ms[n1]*ms[n3]*pow(cl,2)*pow(r[n1][n3],-1) -
        8*ms[n1]*ms[n2]*ms[n3]*pow(GMsol,2)*pow(r[n1][n2],-1)*pow(r[n1][n3],-1) - 4*ms[n2]*pow(GMsol,2)*pow(ms[n1],2)*pow(r[n2][n1],-2) +
        14*GMsol* dotprod<value_type>( v[n1],v[n2])*ms[n1]*ms[n2]*pow(r[n2][n1],-1) +
        2*GMsol* dotprod<value_type>( v[n1],ns[n2][n1])* dotprod<value_type>( v[n2],ns[n2][n1])*ms[n1]*ms[n2]*pow(r[n2][n1],-1) + 6*GMsol*ms[n1]*ms[n2]*vv[n1]*pow(r[n2][n1],-1) +
        6*GMsol*ms[n1]*ms[n2]*vv[n2]*pow(r[n2][n1],-1) + 4*GMsol*ms[n1]*ms[n2]*pow(cl,2)*pow(r[n2][n1],-1) -
        4*ms[n2]*pow(GMsol,2)*pow(ms[n3],2)*pow(r[n2][n3],-2) + 14*GMsol* dotprod<value_type>( v[n3],v[n2])*ms[n2]*ms[n3]*pow(r[n2][n3],-1) -
        12*GMsol* dotprod<value_type>( v[n3],v[n3])*ms[n2]*ms[n3]*pow(r[n2][n3],-1) +
        2*GMsol* dotprod<value_type>( v[n2],ns[n2][n3])* dotprod<value_type>( v[n3],ns[n2][n3])*ms[n2]*ms[n3]*pow(r[n2][n3],-1) + 6*GMsol*ms[n2]*ms[n3]*vv[n2]*pow(r[n2][n3],-1) +
        6*GMsol*ms[n2]*ms[n3]*vv[n3]*pow(r[n2][n3],-1) + 4*GMsol*ms[n2]*ms[n3]*pow(cl,2)*pow(r[n2][n3],-1) -
        8*ms[n1]*ms[n2]*ms[n3]*pow(GMsol,2)*pow(r[n2][n1],-1)*pow(r[n2][n3],-1) - 4*ms[n3]*pow(GMsol,2)*pow(ms[n1],2)*pow(r[n3][n1],-2) +
        14*GMsol* dotprod<value_type>( v[n1],v[n3])*ms[n1]*ms[n3]*pow(r[n3][n1],-1) - 12*GMsol* dotprod<value_type>( v[n3],v[n3])*ms[n1]*ms[n3]*pow(r[n3][n1],-1) +
        2*GMsol* dotprod<value_type>( v[n1],ns[n3][n1])* dotprod<value_type>( v[n3],ns[n3][n1])*ms[n1]*ms[n3]*pow(r[n3][n1],-1) + 6*GMsol*ms[n1]*ms[n3]*vv[n1]*pow(r[n3][n1],-1) +
        6*GMsol*ms[n1]*ms[n3]*vv[n3]*pow(r[n3][n1],-1) + 4*GMsol*ms[n1]*ms[n3]*pow(cl,2)*pow(r[n3][n1],-1) -
        4* dotprod<value_type>( v[n1],v[n1])*ms[n1]*(vv[n1] + 2*pow(cl,2) + 3*GMsol*ms[n2]*(pow(r[n1][n2],-1) + pow(r[n2][n1],-1)) +
           3*GMsol*ms[n3]*(pow(r[n1][n3],-1) + pow(r[n3][n1],-1))) - 4*ms[n3]*pow(GMsol,2)*pow(ms[n2],2)*pow(r[n3][n2],-2) +
        2*GMsol*ms[n2]*ms[n3]*(-4*GMsol*ms[n1] + r[n3][n1]*(7* dotprod<value_type>( v[n2],v[n3]) - 6* dotprod<value_type>( v[n3],v[n3]) +
               dotprod<value_type>( v[n2],ns[n3][n2])* dotprod<value_type>( v[n3],ns[n3][n2]) + 3*vv[n2] + 3*vv[n3] + 2*pow(cl,2)))*pow(r[n3][n1],-1)*pow(r[n3][n2],-1) -
        4* dotprod<value_type>( v[n2],v[n2])*ms[n2]*(vv[n2] + 2*pow(cl,2) + 3*GMsol*ms[n1]*(pow(r[n1][n2],-1) + pow(r[n2][n1],-1)) +
           3*GMsol*ms[n3]*(pow(r[n2][n3],-1) + pow(r[n3][n2],-1)))))/8. ;

        energy -=  troisdemis * quadrupole1 * ms[2] * GMsol / pow(r[1][2],3);
        energy -=  troisdemis * quadrupole1 * ms[0] * GMsol / pow(r[1][0],3);

        /* Old version (made with alternative 2 in the mathematica file Eq_motion_1PN-Will_TEGP.nb) : numerically equal to the one above.
                energy = ((ms[n1]*vv[n1])/ deux  + (ms[n2]*vv[n2])/ deux  + (ms[n3]*vv[n3])/ deux   + ( trois *ms[n1]* clm2 *pow(vv[n1],2))/ huit  +
   ( trois *ms[n2]* clm2 *pow(vv[n2],2))/ huit  + ( trois *ms[n3]* clm2 *pow(vv[n3],2))/ huit  +
   (ms[n1]* clm2 *pow(GMsol,2)*pow(ms[n2],2)*pow(r[n1][n2],-2) -
      GMsol*ms[n1]*ms[n2]*( un  + (deux + troisdemis)* dotprod3d<value_type>( v[n1],v[n2])* clm2  + ( dotprod3d<value_type>( v[n1],ns[n1][n2])* dotprod3d<value_type>( v[n2],ns[n1][n2])* clm2 )/ deux  -
          trois *vv[n1]* clm2 )*pow(r[n1][n2],-1) + ms[n1]* clm2 *pow(GMsol,2)*pow(ms[n3],2)*pow(r[n1][n3],-2) -
      GMsol*ms[n1]*ms[n3]*( un  + (deux + troisdemis)* dotprod3d<value_type>( v[n1],v[n3])* clm2  + ( dotprod3d<value_type>( v[n1],ns[n1][n3])* dotprod3d<value_type>( v[n3],ns[n1][n3])* clm2 )/ deux  -
          trois *vv[n1]* clm2 )*pow(r[n1][n3],-1) + 2*ms[n1]*ms[n2]*ms[n3]* clm2 *pow(GMsol,2)*pow(r[n1][n2],-1)*pow(r[n1][n3],-1))/ deux  +
   (ms[n2]* clm2 *pow(GMsol,2)*pow(ms[n1],2)*pow(r[n2][n1],-2) -
      GMsol*ms[n1]*ms[n2]*( un  + (deux + troisdemis)* dotprod3d<value_type>( v[n2],v[n1])* clm2  + ( dotprod3d<value_type>( v[n1],ns[n2][n1])* dotprod3d<value_type>( v[n2],ns[n2][n1])* clm2 )/ deux  -
          trois *vv[n2]* clm2 )*pow(r[n2][n1],-1) + ms[n2]* clm2 *pow(GMsol,2)*pow(ms[n3],2)*pow(r[n2][n3],-2) -
      GMsol*ms[n2]*ms[n3]*( un  + (deux + troisdemis)* dotprod3d<value_type>( v[n2],v[n3])* clm2  + ( dotprod3d<value_type>( v[n2],ns[n2][n3])* dotprod3d<value_type>( v[n3],ns[n2][n3])* clm2 )/ deux  -
          trois *vv[n2]* clm2 )*pow(r[n2][n3],-1) + 2*ms[n1]*ms[n2]*ms[n3]* clm2 *pow(GMsol,2)*pow(r[n2][n1],-1)*pow(r[n2][n3],-1))/ deux  +
   (ms[n3]* clm2 *pow(GMsol,2)*pow(ms[n1],2)*pow(r[n3][n1],-2) -
      GMsol*ms[n1]*ms[n3]*( un  + (deux + troisdemis)* dotprod3d<value_type>( v[n3],v[n1])* clm2  + ( dotprod3d<value_type>( v[n1],ns[n3][n1])* dotprod3d<value_type>( v[n3],ns[n3][n1])* clm2 )/ deux  -
          trois *vv[n3]* clm2 )*pow(r[n3][n1],-1) + ms[n3]* clm2 *pow(GMsol,2)*pow(ms[n2],2)*pow(r[n3][n2],-2) -
      GMsol*ms[n2]*ms[n3]*( un  + (deux + troisdemis)* dotprod3d<value_type>( v[n3],v[n2])* clm2  + ( dotprod3d<value_type>( v[n2],ns[n3][n2])* dotprod3d<value_type>( v[n3],ns[n3][n2])* clm2 )/ deux  -
          trois *vv[n3]* clm2 )*pow(r[n3][n2],-1) +  deux *ms[n1]*ms[n2]*ms[n3]* clm2 *pow(GMsol,2)*pow(r[n3][n1],-1)*pow(r[n3][n2],-1))/ deux ) ;*/

          impulsion   = ( ms[n1] * v[n1] + ms[n2] * v[n2] + ms[n3] * v[n3] + (ms[n1] * v[n1] * vv[n1] *  clm2 ) / deux  + (ms[n2]*v[n2]*vv[n2]* clm2 )/ deux  + (ms[n3]*v[n3]*vv[n3]* clm2 )/ deux  +
   (-(GMsol*ms[n1]*ms[n2]*(v[n1] *  clm2  +  dotprod3d<value_type>( v[n1],ns[n1][n2])*ns[n1][n2]* clm2 )*pow(r[n1][n2],-1)) -
      GMsol*ms[n1]*ms[n3]*(v[n1]* clm2  +  dotprod3d<value_type>( v[n1],ns[n1][n3])*ns[n1][n3]* clm2 )*pow(r[n1][n3],-1))/ deux  +
   (-(GMsol*ms[n1]*ms[n2]*(v[n2]* clm2  +  dotprod3d<value_type>( v[n2],ns[n2][n1])*ns[n2][n1]* clm2 )*pow(r[n2][n1],-1)) -
      GMsol*ms[n2]*ms[n3]*(v[n2]* clm2  +  dotprod3d<value_type>( v[n2],ns[n2][n3])*ns[n2][n3]* clm2 )*pow(r[n2][n3],-1))/ deux  +
   (-(GMsol*ms[n1]*ms[n3]*(v[n3]* clm2  +  dotprod3d<value_type>( v[n3],ns[n3][n1])*ns[n3][n1]* clm2 )*pow(r[n3][n1],-1)) -
      GMsol*ms[n2]*ms[n3]*(v[n3]* clm2  +  dotprod3d<value_type>( v[n3],ns[n3][n2])*ns[n3][n2]* clm2 )*pow(r[n3][n2],-1))/ deux  ) ;

                //V[i,:] = P[i,:] / Mt * ( un  - 0.5*E[i]/(Mt*cl**2) - 0.5 * VV0/cl**2 )

                center_of_mass =  un /Mt * ( un  - ( energy * clm2 *pow(Mt,-1))/ deux  - ( clm2 *VV0)/ deux )*(ms[n1]*xs[n1] + ms[n2]*xs[n2] + ms[n3]*xs[n3] + (ms[n1]*vv[n1]*xs[n1]* clm2 )/ deux  +
     (ms[n2]*vv[n2]*xs[n2]* clm2 )/ deux  + (ms[n3]*vv[n3]*xs[n3]* clm2 )/ deux  -
     (GMsol* clm2 *(ms[n1]*ms[n2]*xs[n1]*pow(r[n1][n2],-1) + ms[n1]*ms[n3]*xs[n1]*pow(r[n1][n3],-1)))/ deux  -
     (GMsol* clm2 *(ms[n1]*ms[n2]*xs[n2]*pow(r[n2][n1],-1) + ms[n2]*ms[n3]*xs[n2]*pow(r[n2][n3],-1)))/ deux  -
     (GMsol* clm2 *(ms[n1]*ms[n3]*xs[n3]*pow(r[n3][n1],-1) + ms[n2]*ms[n3]*xs[n3]*pow(r[n3][n2],-1)))/ deux ) ;

      return;
}

/*
    cpdef _Compute_emission_times_check_(self) :
        cdef DTYPE_t Te0 = 0.
        cdef int i, j = 0
        cdef np.ndarray[DTYPE_t, ndim = 1] ta
        cdef np.ndarray[DTYPE_t, ndim = 1] Te1

        if self._interpolation :
            ta = self._interp_times
        else :
            ta = self._toas

        Te1 = np.zeros(ta.size, dtype = DTYPE)

        for i in range(ta.size) :
            Te0 = ta[i]
            Te1[i] = ta[i] - self._Compute_total_delay_check_( ta[i] - self._Compute_total_delay_check_(Te0) )
            j = 0
            while ( abs(Te1[i] - Te0 ) > 10.**(-14) and j < 100) :
                j += 1
                Te0 = Te1[i]
                Te1[i] = ta[i] - self._Compute_total_delay_check_( ta[i] - self._Compute_total_delay_check_(Te0) )

            if j == 100 : print " Iteration did not converge in Retarded_times ! Time number : ", i

        return Te1*/




void IntegralePrems3_1PN_double ( const double s0[6], const double s1[6], const double s2[6],
                         double m0, double m1, double m2, double quadrupole1,
                         double  center_of_mass[3], double  impulsion[3],
                         double &  energy  )
// Double version of IntegralePrems3_1PN, using only basic arrays (good for interfacing with python)
{
    valarray<value_type>  cm;
    valarray<value_type>  imp;
    value_type nrj;
    value_type ss0[6], ss1[6], ss2[6] ;

    for(int i = 0; i < 6 ; i++)
    {
        ss0[i] = static_cast<value_type>(s0[i]);
        ss1[i] = static_cast<value_type>(s1[i]);
        ss2[i] = static_cast<value_type>(s2[i]);
    }

    IntegralePrems3_1PN( ss0, ss1, ss2,
                         static_cast<value_type>(m0), static_cast<value_type>(m1), static_cast<value_type>(m2), static_cast<value_type>(quadrupole1),
                         cm, imp,
                         nrj  );

    energy = static_cast<double>(nrj);

    for(int i = 0; i < 3 ; i++)
    {
        center_of_mass[i] = static_cast<double>(cm[i]);
        impulsion[i] = static_cast<double>(imp[i]);
    };


}



void IntegralePrems3_1PN(const value_type sv0[6], const value_type sv1[6], const value_type sv2[6],
                         value_type m0, value_type m1, value_type m2, value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg, const valarray<valarray<value_type>>  gammabar, const valarray<valarray<valarray<value_type>>>  betabar,
//                          const value_type Gg[3][3], const value_type gammabar[3][3], const value_type betabar[3][3][3],
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
{
        /*
            R, L, P, E = IntegralePrems3(sv0,sv1,sv2,m0,m1,m2) \n
            sv0, sv1, sv2 are the state vectors (postion and velocity) of bodies of masses m0, m1, m2 \n
            Body 1 has a spin-induced quadrupole moment "quadrupole1" with spin axis perpendicular to the orbital plane of the inner binary. Exact if 3 the 3 bodies have coplanar orbits.
            R is the center of mass, L the angular momentum, P the impulsion and E the energy of the system.
            To compute the first integrals for a serie of n times, sv0 must have the shape (n times, 6 components), then the first integrals will have the same shape.
            If sv0.shape = (6), then the vectorial first integrals will also have this shape.

            Energy is returned in units of Msol (m/s)^2
        */

        static const int nbody = 3; // change this to generalize to any number of bodies

        int n1;
        int n2;
        int n3;
        const value_type cl = clight ;
        const value_type clm2 = pow(cl,-2) ;

        static const value_type troishuitiemes = trois/huit;
        value_type sinter = zero;
        value_type sinter1 = zero;
        value_type sinter2 = zero;
        value_type sinter21 = zero;
        valarray<value_type> inter(3);
        valarray<value_type> inter1(3);
        valarray<value_type> inter2(3);


        valarray<value_type> s0(sv0, 6);
        valarray<value_type> s1(sv1, 6);
        valarray<value_type> s2(sv2, 6);

        valarray<valarray<value_type> > xs(3) ;

            xs[0] = valarray<value_type>( s0[slice(0,3,1)] );
            xs[1] = valarray<value_type>( s1[slice(0,3,1)] );
            xs[2] = valarray<value_type>( s2[slice(0,3,1)] );


        valarray<valarray<value_type> > v(3) ;
            v[0] = valarray<value_type>( s0[slice(3,3,1)] );
            v[1] = valarray<value_type>( s1[slice(3,3,1)] );
            v[2] = valarray<value_type>( s2[slice(3,3,1)] );

        valarray<value_type> vv(3) ;
            vv[0] = sumsquares3d<value_type>( v[0] );
            vv[1] = sumsquares3d<value_type>( v[1] );
            vv[2] = sumsquares3d<value_type>( v[2] );

        valarray<value_type> ms(3) ;
            ms[0] = m0 ;
            ms[1] = m1 ;
            ms[2] = m2 ;

        value_type Mt = m0 + m1 + m2 ;
        value_type GMsol = Ggrav * Msol ;

        valarray<value_type> VV0( ( ms[0]*v[0] + ms[1]*v[1] + ms[2]*v[2] ) / Mt );

        valarray<valarray<value_type> > r(3) ;
            r[0] = valarray<value_type>( 3 ) ;
                r[0][0] = un ;
                r[0][1] = norm3d<value_type>( xs[1] - xs[0] ) ;
                r[0][2] = norm3d<value_type>( xs[2] - xs[0] ) ;
            r[1] = valarray<value_type>(3);
                r[1][0] = norm3d<value_type>(xs[0] - xs[1]) ;
                r[1][1] = un ;
                r[1][2] = norm3d<value_type>(xs[2] - xs[1] ) ;
            r[2] = valarray<value_type>(3);
                r[2][0] = norm3d<value_type>(xs[0] - xs[2]) ;
                r[2][1] = norm3d<value_type>(xs[1] - xs[2]) ;
                r[2][2] = un ;

        valarray<valarray<valarray<value_type> > > ns(3) ;
            ns[0] = valarray<valarray<value_type> >(3);
                ns[0][0] = valarray<value_type>(3);
                ns[0][1] = valarray<value_type>( xs[1] - xs[0] ) ;
                ns[0][2] = valarray<value_type>( xs[2] - xs[0] ) ;
            ns[1] = valarray<valarray<value_type> >(3);
                ns[1][0] = valarray<value_type>( xs[0] - xs[1] ) ;
                ns[1][1] = valarray<value_type>(3);
                ns[1][2] = valarray<value_type>( xs[2] - xs[1] ) ;
            ns[2] = valarray<valarray<value_type> >(3);
                ns[2][0] = valarray<value_type>( xs[0] - xs[2] ) ;
                ns[2][1] = valarray<value_type>( xs[1] - xs[2] );
                ns[2][2] = valarray<value_type>(3);
            for (int j = 0 ; j<3 ; ++j)
            {
                for (int k = 0 ; k<3 ; ++k) ns[j][k] /= r[j][k];
            }



    // Compute the Hamiltonian
     energy = zero;
    for (n1 = 0 ; n1 < nbody ; n1++)
    {
        sinter1 = undemi  * vv[n1] + troishuitiemes* vv[n1] * vv[n1]   * clm2;

        sinter2 = zero;
        for (n2 = 0 ; n2 < nbody ; n2++)
        {
          if (n2 != n1)
          {
              sinter = -dotprod3d<value_type>(ns[n1][n2], v[n1])* dotprod3d<value_type>(ns[n1][n2], v[n2]) * undemi;
              sinter -=   dotprod3d<value_type>(v[n1],v[n2])*(septdemis + deux * gammabar[n1][n2]);
              sinter += (vv[n1]) *  (trois + deux*gammabar[n1][n2]) ; // could be put in the parent loop...
              sinter *= clm2;
              sinter -= un;

              sinter21 = zero;
              for (n3 = 0 ; n3 < nbody ; n3++)
              {
                  if (n3 != n1)
                  {
                      sinter21 +=Gg[n1][n3] * ms[n3] / r[n1][n3]* (un + deux * betabar[n1][n2][n3] ) ;
                  }
              }
              sinter21 *= clm2 * GMsol;

              sinter2 +=  undemi*(sinter + sinter21 ) * (GMsol * Gg[n1][n2]*ms[n2] / r[n1][n2]);

//               printf("jjj %.5Le %.5Le %.5Le %.5Le \n", undemi  * vv[n1], troishuitiemes* vv[n1] * vv[n1]   * clm2, sinter2, sinter21* (GMsol * Gg[n1][n2]*ms[n2] / r[n1][n2])*0.5);
          }
        }

        energy += ms[n1] * (sinter1 + sinter2);
    }

    // Add quadrupole component due to object 1
    energy -=  troisdemis * quadrupole1 * ms[2] * Ggrav / pow(r[1][2],3);
    energy -=  troisdemis * quadrupole1 * ms[0] * Ggrav / pow(r[1][0],3);

// This expression is derived from the Lagrangian in Soffel's book "Relativity in Astrometry, Celestial Mechanics and Geodesy" , also in Damour And Taylor 1992, Eq_motion_1PN-Will_TEGP.nb.
// Agrees with the formula given by Will, 2014, "Incorporating post-Newtonian effects in N-body dynamics", formula 3.2
// !!! Don't forget to remove the mass energy from it !!!

        impulsion = valarray<value_type>(zero,3);
        for (n1 = 0 ; n1 < nbody ; n1++)
        {
            sinter = un + undemi*vv[n1]   * clm2;
            inter1 = v[n1] * sinter;
            inter2 = zero;
            for (n2 = 0 ; n2 < nbody ; n2++)
            {
                if (n2 != n1)
                {
                    inter = - dotprod3d<value_type>(v[n2], ns[n1][n2]) * ns[n1][n2];
                    inter -= v[n2];
                    inter2 += (Gg[n1][n2] * GMsol * ms[n2] *clm2 * undemi/ r[n1][n2]) * inter;
                }
            }
            impulsion += (inter1 + inter2) * ms[n1];
        }

        center_of_mass = valarray<value_type>(zero,3);
        for (n1 = 0 ; n1 < nbody ; n1++)
        {
            sinter1 = un + undemi*vv[n1]   * clm2;
            sinter2 = zero;
            for (n2 = 0 ; n2 < nbody ; n2++)
            {
                if (n2 != n1)
                {
                    sinter2 -= (Gg[n1][n2] * ms[n2] / r[n1][n2]);
                }
            }
            sinter2 *= GMsol  *clm2 * undemi;
            center_of_mass += (sinter1 + sinter2) * ms[n1] * xs[n1];
        }
        center_of_mass /= (energy*clm2 + Mt);

      return;
}



void IntegralePrems_0PN(const valarray<value_type> sv, const valarray<value_type> Ms, 
                         value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy,  valarray<value_type> & angular_momentum  )
{
        /*  sv is a state vector such that : sv = [x0,y0,z0,x1,y1,z1...,vx0,vy0,vz0,vx1,vy1,vz1,...] where (x0,y0,z0) is the position of body 1 and (vx0,vy0,vz0) its velocity.
         * There can be as many bodies as desired. 
         * Ms lists the corresponding masses (in units of Msun). 
            
            Body 1 has a spin-induced quadrupole moment "quadrupole1" with spin axis perpendicular to the orbital plane of the inner binary. Exact if 3 the 3 bodies have coplanar orbits.
            R is the center of mass, L the angular momentum, P the impulsion and E the energy of the system.
            To compute the first integrals for a serie of n times, sv0 must have the shape (n times, 6 components), then the first integrals will have the same shape.
            If sv0.shape = (6), then the vectorial first integrals will also have this shape.

            Energy is returned in units of Msol (m/s)^2
        */

        int nbody = sv.size() / 6; 
        int n1;
        int n2;

        value_type sinter1 = zero;
        value_type sinter2 = zero;

        valarray<valarray<value_type> > xs(nbody) ;
        valarray<valarray<value_type> > v(nbody) ;
        valarray<value_type> vv(nbody) ;
        
        for (n1 = 0; n1 < nbody ; n1++)
        {
            xs[n1] = sv[slice(n1*3,3,1)];
            v[n1] =  sv[slice(nbody*3 + n1*3,3,1)];
            vv[n1] = sumsquares3d<value_type>(v[n1]);
        }
    
        value_type Mt = Ms.sum() ;
        value_type GMsol = Ggrav * Msol ;

        valarray<valarray<value_type> > r(nbody) ;
        
        for (n1 = 0; n1 < nbody ; n1++)
        {
            r[n1].resize(nbody);
            for (n2 = 0; n2 < nbody ; n2++)
            {
                r[n1][n2] = norm3d<value_type>( xs[n1] - xs[n2] ) ;
            }
            r[n1][n1] = un;
            
        }
           
    // Compute the Hamiltonian
     energy = zero;
    for (n1 = 0 ; n1 < nbody ; n1++)
    {
        sinter1 = undemi  * vv[n1];

        sinter2 = zero;
        for (n2 = 0 ; n2 < nbody ; n2++)
        {
          if (n2 != n1)
          {
             sinter2 +=  undemi*(-un ) * (GMsol * Gg[n1][n2]*Ms[n2] / r[n1][n2]);
          }
        }

        energy += Ms[n1] * (sinter1 + sinter2);
    }

    // Add quadrupole component due to object 1
    energy -=  troisdemis * quadrupole1 * Ms[2] * Ggrav / pow(r[1][2],3);
    energy -=  troisdemis * quadrupole1 * Ms[0] * Ggrav / pow(r[1][0],3);

        impulsion = valarray<value_type>(zero,3);
        for (n1 = 0 ; n1 < nbody ; n1++)
        {
            impulsion += v[n1] * Ms[n1];
        }

        center_of_mass = valarray<value_type>(zero,3);
        for (n1 = 0 ; n1 < nbody ; n1++)
        {
            center_of_mass += Ms[n1] * xs[n1];
        }
        center_of_mass /= Mt;

    // Computing angular momentum to Newtonian level
       angular_momentum = valarray<value_type>(zero,3);
       for (n1 = 0; n1 < nbody; n1++)
        {
            angular_momentum += Ms[n1] * crossprod(xs[n1],v[n1]);               
        }
      return;
}


void IntegralePrems3_1PN(const valarray<value_type> sv, const valarray<value_type> Ms, 
                         value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg, const valarray<valarray<value_type>>  gammabar, const valarray<valarray<valarray<value_type>>>  betabar,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
{
        /*  sv is a state vector such that : sv = [x0,y0,z0,x1,y1,z1...,vx0,vy0,vz0,vx1,vy1,vz1,...] where (x0,y0,z0) is the position of body 1 and (vx0,vy0,vz0) its velocity.
         * There can be as many bodies as desired. 
         * Ms lists the corresponding masses (in units of Msun). 
            
            Body 1 has a spin-induced quadrupole moment "quadrupole1" with spin axis perpendicular to the orbital plane of the inner binary. Exact if 3 the 3 bodies have coplanar orbits.
            R is the center of mass, L the angular momentum, P the impulsion and E the energy of the system.
            To compute the first integrals for a serie of n times, sv0 must have the shape (n times, 6 components), then the first integrals will have the same shape.
            If sv0.shape = (6), then the vectorial first integrals will also have this shape.

            Energy is returned in units of Msol (m/s)^2
        */

        int nbody = sv.size() / 6; 
        int n1;
        int n2;
        int n3;
        const value_type cl = clight ;
        const value_type clm2 = pow(cl,-2) ;

        static const value_type troishuitiemes = trois/huit;
        value_type sinter = zero;
        value_type sinter1 = zero;
        value_type sinter2 = zero;
        value_type sinter21 = zero;
        valarray<value_type> inter(3);
        valarray<value_type> inter1(3);
        valarray<value_type> inter2(3);

        valarray<valarray<value_type> > xs(nbody) ;
        valarray<valarray<value_type> > v(nbody) ;
        valarray<value_type> vv(nbody) ;
        
        for (n1 = 0; n1 < nbody ; n1++)
        {
            xs[n1] = sv[slice(n1*3,3,1)];
            v[n1] =  sv[slice(nbody*3 + n1*3,3,1)];
            vv[n1] = sumsquares3d<value_type>(v[n1]);
        }
    
        value_type Mt = Ms.sum() ;
        value_type GMsol = Ggrav * Msol ;

        valarray<value_type> VV0(3,0.);
        valarray<valarray<value_type> > r(nbody) ;
        
        for (n1 = 0; n1 < nbody ; n1++)
        {
            VV0 += Ms[n1] * v[n1];
            r[n1].resize(nbody);
            for (n2 = 0; n2 < nbody ; n2++)
            {
                r[n1][n2] = norm3d<value_type>( xs[n1] - xs[n2] ) ;
            }
            r[n1][n1] = un;
            
        }
        VV0 /= Mt;
           

        valarray<valarray<valarray<value_type> > > ns(nbody) ;
        for (n1 = 0; n1 < nbody ; n1++)
        {
            ns[n1] = valarray<valarray<value_type> >(nbody);
            for (n2 = 0; n2 < nbody ; n2++)
            {
                ns[n1][n2] = xs[n2] - xs[n1];
                ns[n1][n2] /= r[n1][n2];
            }
            ns[n1][n1].resize(3,zero);
        }

    // Compute the Hamiltonian
     energy = zero;
    for (n1 = 0 ; n1 < nbody ; n1++)
    {
        sinter1 = undemi  * vv[n1] + troishuitiemes* vv[n1] * vv[n1]   * clm2;

        sinter2 = zero;
        for (n2 = 0 ; n2 < nbody ; n2++)
        {
          if (n2 != n1)
          {
              sinter = -dotprod3d<value_type>(ns[n1][n2], v[n1])* dotprod3d<value_type>(ns[n1][n2], v[n2]) * undemi;
              sinter -=   dotprod3d<value_type>(v[n1],v[n2])*(septdemis + deux * gammabar[n1][n2]);
              sinter += (vv[n1]) *  (trois + deux*gammabar[n1][n2]) ; 
              sinter *= clm2;
              sinter -= un;

              sinter21 = zero;
              for (n3 = 0 ; n3 < nbody ; n3++)
              {
                  if (n3 != n1)
                  {
                      sinter21 +=Gg[n1][n3] * Ms[n3] / r[n1][n3]* (un + deux * betabar[n1][n2][n3] ) ;
                  }
              }
              sinter21 *= clm2 * GMsol;

              sinter2 +=  undemi*(sinter + sinter21 ) * (GMsol * Gg[n1][n2]*Ms[n2] / r[n1][n2]);
          }
        }

        energy += Ms[n1] * (sinter1 + sinter2);
    }

    // Add quadrupole component due to object 1
    energy -=  troisdemis * quadrupole1 * Ms[2] * Ggrav / pow(r[1][2],3);
    energy -=  troisdemis * quadrupole1 * Ms[0] * Ggrav / pow(r[1][0],3);

// This expression is derived from the Lagrangian in Soffel's book "Relativity in Astrometry, Celestial Mechanics and Geodesy" , also in Damour And Taylor 1992, Eq_motion_1PN-Will_TEGP.nb.
// Agrees with the formula given by Will, 2014, "Incorporating post-Newtonian effects in N-body dynamics", formula 3.2
// !!! Don't forget to remove the mass energy from it !!!

        impulsion = valarray<value_type>(zero,3);
        for (n1 = 0 ; n1 < nbody ; n1++)
        {
            sinter = un + undemi*vv[n1]   * clm2;
            inter1 = v[n1] * sinter;
            inter2 = zero;
            for (n2 = 0 ; n2 < nbody ; n2++)
            {
                if (n2 != n1)
                {
                    inter = - dotprod3d<value_type>(v[n2], ns[n1][n2]) * ns[n1][n2];
                    inter -= v[n2];
                    inter2 += (Gg[n1][n2] * GMsol * Ms[n2] *clm2 * undemi/ r[n1][n2]) * inter;
                }
            }
            impulsion += (inter1 + inter2) * Ms[n1];
        }

        center_of_mass = valarray<value_type>(zero,3);
        for (n1 = 0 ; n1 < nbody ; n1++)
        {
            sinter1 = un + undemi*vv[n1]   * clm2;
            sinter2 = zero;
            for (n2 = 0 ; n2 < nbody ; n2++)
            {
                if (n2 != n1)
                {
                    sinter2 -= (Gg[n1][n2] * Ms[n2] / r[n1][n2]);
                }
            }
            sinter2 *= GMsol  *clm2 * undemi;
            center_of_mass += (sinter1 + sinter2) * Ms[n1] * xs[n1];
        }
        center_of_mass /= (energy*clm2 + Mt);

      return;
}

void IntegralePrems3_1PN(const valarray<value_type> sv, const valarray<value_type> Ms, 
                         value_type quadrupole1,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
// Same as IntegralePrems3_1PN, just assume plain GR.
{
    int i,j,nbody;
    nbody = Ms.size();
    valarray<valarray<value_type>> Gg(nbody);
    valarray<valarray<value_type>>  gammabar(nbody);
    valarray<valarray<valarray<value_type>>>  betabar(nbody);
    
    for (i = 0; i < nbody; i++)
    {
        gammabar[i] = valarray<value_type>(0.,nbody);
        Gg[i] = valarray<value_type>(0.,nbody);
        betabar[i].resize(nbody);
        for (j = 0; j < nbody; j++)
        {
            Gg[i][j] = 1.;
            betabar[i][j] = valarray<value_type>(0.,nbody);
        }
    }
    
    IntegralePrems3_1PN(sv, Ms, quadrupole1, Gg, gammabar, betabar, center_of_mass, impulsion, energy);
}


void IntegralePrems3_1PN_extra(const valarray<value_type> sv, const valarray<value_type> Ms,
                         value_type quadrupole1,
                         const valarray<valarray<value_type>> Gg, const valarray<valarray<value_type>>  gammabar, const valarray<valarray<valarray<value_type>>>  betabar,
                         valarray<value_type> & center_of_mass, valarray<value_type> & impulsion,
                         value_type & energy  )
// Uses IntegralePrems3_1PN for the first 3 bodies and consider any extra body as Newtonian
// sv = [position body 1, position body 2,...., velocity body 1, velocity body 2...]
// It is assumed that there are at least three bodies
{
    value_type sv0[6], sv1[6], sv2[6];
    int i,j;
    int nbody = Ms.size();
    value_type m0 = Ms[0];
    value_type m1 = Ms[1];
    value_type m2 = Ms[2];
    valarray<value_type> velo(3);
    valarray<value_type> dr(3);
    const value_type cl = clight ;
    const value_type clm2 = pow(cl,-2) ;
    
    for (i = 0 ; i < 3; i++) 
    {
        sv0[i] = sv[i];
        sv0[i+3] = sv[i + 3 * nbody];
        
        sv1[i] = sv[i+3];
        sv1[i+3] = sv[i + 3 + 3 * nbody];
        
        sv2[i] = sv[i+6];
        sv2[i+3] = sv[i + 6 + 3 * nbody];
    }
    
    IntegralePrems3_1PN(sv0, sv1, sv2, m0,m1,m2, quadrupole1, Gg, gammabar, betabar, 
                        center_of_mass, impulsion, energy);
    
    value_type nrj_1pn = energy; 
    center_of_mass *= nrj_1pn* clm2 + m0 + m1 + m2;
    
    for (i = 3; i < nbody ; i++)
    {
        center_of_mass += Ms[i] * sv[std::slice(i *3, 3, 1)] ;
        velo = sv[std::slice(3*nbody + i *3, 3, 1)];
        impulsion += Ms[i] * velo;
        energy += 0.5*Ms[i]*sumsquares3d<value_type>(velo);
        for (j = 0 ; j < i ; j++)
        {
            dr = sv[std::slice(i *3, 3, 1)];
            dr -= sv[std::slice(j *3, 3, 1)];
            energy -= GMsol * Ms[i] * Ms[j] / norm3d(dr); 
        }
    }
    center_of_mass /= Ms.sum() + nrj_1pn*clm2;
}
