// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 *
 */

# include <cmath>
# include "Spline.h"
# include "AllTheories3Bodies.h"
# include "Utilities.h"
# include <valarray>

/* Routines Delays_Brut and Delays_Brut_fullgeometric are legacies from previous versions of the code.
 * The current version uses  "Delays_Brut_geometric_local" and "Delays_Brut_nogeometric".
 * G. Voisin, October 2017
 */

void Delays_Brut(value_type tis[], const value_type& t0,
                 const long int& ntis, const long int& nt0,
                 value_type ** sp, value_type **  si, value_type ** so,
                 const value_type& Mi, const value_type& Mo, const value_type& freq,
                 bool roemer, bool einstein, bool shapiro, bool aberration,
                 value_type delay[]
                ) {

    long int i = 0;
    //value_type delay[ntis] ;
    value_type * deindt ;

    value_type ein ;
    value_type nvp2 ;
    value_type nrpi ;
    value_type nrpo ;


    // Compute Roemer
  if ( roemer == true ) {
    for (i=0; i < ntis ; ++i) {
        delay[i] = sp[i][1] / ( clight  * daysec ) ;
    }
  }
  else {
    for (i=0; i < ntis ; ++i) {
        delay[i] = zero ;
    }
  }

    // Compute Einstein delay
   if ( einstein == true ) {
       deindt = new value_type[ntis];

   //cdef np.ndarray[DTYPE_t, ndim =1] t = np.zeros(self._ninterp, DTYPE)
      //  cdef np.ndarray[DTYPE_t, ndim =1] tis = np.zeros(self._ninterp, DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] nrpi = np.zeros(self._ninterp, DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] nrpo = np.zeros(self._ninterp, DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] gravpot = np.zeros(self._ninterp, DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] nvp2 = np.zeros(self._ninterp, DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] knots = np.zeros(self._ninterp, DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] ttau = np.zeros(self._ntoas,DTYPE)
//         cdef np.ndarray[DTYPE_t, ndim =1] toas = self._toas
//         cdef np.ndarray[np.int_t, ndim =1] ks
//         cdef np.ndarray[np.float64_t, ndim =1] test = np.zeros(self._ninterp,np.float64) #! test !
//
//         cdef ntoas = self._ntoas
//         cdef DTYPE_t t0, toa = 0.
//         cdef int i = 0


//        ttau = np.zeros(ntoas,DTYPE)
  //      if times is None :
    //        toas = self._interp_times
//             ks = np.arange(self._ninterp, dtype=np.int)
//             ks[self._ninterp-1] = ks[self._ninterp-2] # Deal with the limit condition
//         else :    Arange(tis, ntis, a, dtis) ;

//             ks = - np.ones(times.size, dtype=np.int)
      //      toas = times

        //tis = self._interp_times

        for (i = 0 ; i < ntis ; ++i) {
            nrpi = pow( sp[i][0] - si[i][0], deux) + pow( sp[i][1] - si[i][1], deux)  + pow( sp[i][2] - si[i][2], deux) ;
            nrpi = sqrt( nrpi ) ;
            nrpo = pow( sp[i][0] - so[i][0], deux) + pow( sp[i][1] - so[i][1], deux)  + pow( sp[i][2] - so[i][2], deux) ;
            nrpo = sqrt( nrpo ) ;
            nvp2 = pow( sp[i][3], deux ) + pow( sp[i][4], deux ) + pow( sp[i][5], deux ) ;
            deindt[i] = Mi / nrpi ;
            deindt[i] += Mo / nrpo ;
            deindt[i] *= Msol * Ggrav / pow( clight, 2 ) ;
            deindt[i] += undemi * nvp2 / pow( clight, 2 ) ;
        };
       // gravpot =  self._mi / nrpi
        //gravpot +=  self._mo / nrpo
        //gravpot *= Msol * Ggrav / c**2
        //nvp2 /= c**2
        //knots = 0.5 * nvp2 + gravpot

        Spline splineein( tis, deindt, ntis ) ;

        //dttau_dtau = cMySpline( tis, knots ) # UnivariateSpline(self._interp_times, knots, k=2, s=0.) # d(t - tau)/ dt where tau is the proper time of an observer at the pulsar center.

        //self._knots_test = knots #! Test !
        //self._splineein_y2_test = dttau_dtau.Get_y2s()

        //t0 = self._treference
        ein = splineein.Integrate( t0, tis[0],
                                    nt0, nt0 + 1,
                                    0, 1 ) ;
        delay[0] += ein ;
        //ttau[0] = dttau_dtau.Integrate( t0, toas[0] , kminstart = self._treference_in_interp, kmaxstart = self._treference_in_interp + 1, kminend = 0 , kmaxend =  1  )#dttau_dtau.integral( t0, toas[i] )   # could be opitmized summing the differences

        for (i = 1 ; i < ntis - 1 ; ++i ) {// in range(1,ntoas) :
            //deindt[i] = deindt[ i - 1 ] ;
            ein += splineein.Integrate( tis[ i - 1 ] , tis[ i ],
                                        i - 1, i,
                                        i, i + 1 ) ;
            delay[i] += ein ;
            //ttau[i] = ttau[ i - 1 ] + dttau_dtau.Integrate( toas[i-1], toas[i],
              //                                             kminstart = max(ks[i-1], 0), kmaxstart = ks[i-1] + 1, kminend = max( ks[i], 0 ), kmaxend = ks[i] + 1  )#dttau_dtau.integral( toas[i-1], toas[i] )
        }

        // Deal with the particular case (for indices) of the last interpolation point
        //deindt[ ntis - 1 ] = deindt[ ntis - 2 ] ;
        ein += splineein.Integrate( tis[ ntis - 2 ] , tis[ ntis - 1 ],
                                    ntis - 2, ntis - 1,
                                    ntis - 2, ntis - 1 ) ;
        delay[ ntis - 1 ] += ein ;

        //delete &splineein;
        delete[] deindt;
                  printf( "ein : %.19Le \n", ein ) ;
   }


    // Compute Shapiro delay         (    For more ref : Backers and Hellings , Annual review of astronomy and astrophysics, 1986 )

    if (shapiro == true ) {
        //cdef np.ndarray[DTYPE_t, ndim = 2] pulsar = self._state_vectors[:,0:3]/c  # Position of the pulsar
        //cdef np.ndarray[DTYPE_t, ndim = 2] incomp = self._state_vectors[:,3:6]/c  # Position of the inner companion
        //cdef np.ndarray[DTYPE_t, ndim = 2] oucomp = self._state_vectors[:,6:9]/c  # Position of the outer companion
        //cdef np.ndarray[DTYPE_t, ndim = 1] Shap = np.zeros(self._ntoas, dtype=DTYPE)
        //cdef np.ndarray[DTYPE_t, ndim = 2] pi = np.zeros((self._ntoas, 3), dtype=DTYPE)
        //cdef np.ndarray[DTYPE_t, ndim = 2] po = np.zeros((self._ntoas, 3), dtype=DTYPE)
        value_type pi_dot_nss, po_dot_nss ;
        value_type Shap ;
        //cdef DTYPE_t npi, npo
        const value_type cst = 4.92521372097374e-06L ;            // = Ggrav * Msol / c**3 (unit = second)
        //cdef np.ndarray[DTYPE_t, ndim = 1] nss = np.array([0., 1., 0.], dtype=DTYPE)
        value_type PNParameter = 1.L ; // post newtonian parameter gamma, temporary


        //pi = pulsar - incomp
        //po = pulsar - oucomp
        for (i=0 ; i < ntis ; ++i ) {// in range(len(pulsar[:,1])):
            pi_dot_nss = sp[i][1] - si[i][1] ; //  np.dot(pi[i,:], nss)
            po_dot_nss = sp[i][1] - so[i][1] ; //  np.dot(po[i,:], nss)
            nrpo = sqrt( pow( sp[i][0] - so[i][0], 2) + pow( po_dot_nss, 2) + pow( sp[i][2] - so[i][2], 2) ) ;
            nrpi = sqrt( pow( sp[i][0] - si[i][0], 2) + pow( pi_dot_nss, 2) + pow( sp[i][2] - si[i][2], 2) ) ;
            Shap = Mi * log(  ( - pi_dot_nss + nrpi ) / clight ) ; // clight brings in a scaling parameter of one light second
            Shap += Mo * log( (nrpo - po_dot_nss) / clight ) ;
            Shap *= -(1. + PNParameter) * cst ;
            delay[i] += Shap / daysec ;
            //npo = np.linalg.norm(po[i,:])
            //npi = np.linalg.norm(pi[i,:])
        }
            //Shap[i] = self._mi * np.log( ( - pi_dot_nss + npi ) / scale )              # Computing Shapiro delay at arrival time (zeroth order)
            //Shap[i] += self._mo * np.log( (npo - po_dot_nss) / scale )
            //Shap[i] *= -(1. + PNParameter) * cst #! Test !
 printf( "shap : %.19Le \n", Shap ) ;
    }


//     // Compute aberration delay ( Aberration delay as in Smarr and Taylor 1976 or Damour & Deruelle 1986  )
    if ( aberration == true )    {

  //      pulsar = self._interp_state_vectors[:,0:3]/c  # Position of the pulsar in light seconds
    //    pulsarv = self._interp_state_vectors[:,9:12]/c  # Velocity of the pulsar in units of c

//        nbts = self._ninterp
  //      if times is None :
    //        ts = self._interp_times
      //  else :
        //    ts = times

        value_type rotaxis[3] ;//np.cross(pulsar[0,:], pulsarv[0,:] ) # Assume the spinning axis of the pulsar is aligned with its angular momentum.
        value_type nrotaxis ;
        value_type transverse0, transverse2 ;

        rotaxis[0] = ( sp[0][1] * sp[0][5] - sp[0][2]*sp[0][4] ) / pow( clight, 2 );
        rotaxis[1] = ( sp[0][2] * sp[0][3] - sp[0][0]*sp[0][5] ) / pow( clight, 2 ) ;
        rotaxis[2] = ( sp[0][0] * sp[0][4] - sp[0][1]*sp[0][3] ) / pow( clight, 2 ) ;

        nrotaxis = sqrt( pow( rotaxis[0], 2) + pow( rotaxis[1], 2) + pow( rotaxis[2], 2) ) ; //np.linalg.norm(rotaxis)
        nrotaxis = nrotaxis / ( pow( rotaxis[2] , 2 ) + pow ( rotaxis[0], 2 ) ) ;
        transverse0 = - rotaxis[2] * nrotaxis ; //np.cross(rotaxis, nss)
        transverse2 = rotaxis[0] * nrotaxis ;
        //transverse /= np.sum(transverse**2)

        //aberration = np.ndarray(nbts, dtype=DTYPE)
        //knots = np.ndarray(self._ninterp, dtype=DTYPE)

        for (i= 0 ; i < ntis ; ++i ) { // in range(self._ninterp) :
            //knots[i] = np.dot( transverse, pulsarv[i] )
            delay[i] += ( transverse0 * sp[i][3] + transverse2 * sp[i][5] ) / clight / freq ;
        }

        //spline = cMySpline(self._interp_times, knots ) # UnivariateSpline(self._interp_times, knots, k=2, s=0.)

        //for i in range(nbts) :
          //  aberration[i] = spline( ts[i] ) / self._fitted_params[0]
    }

        return ;
}







void Rotate_PSB_to_SSBh(valarray<value_type> & vector3, value_type alpha, value_type delta)
{
    value_type sa = sin(alpha);
    value_type ca = cos(alpha);
    value_type sd = sin(delta);
    value_type cd = cos(delta);
    value_type vi[3] ;

    vi[0] = vector3[0] ;
    vi[1] = vector3[1] ;
    vi[2] = vector3[2] ;

    vector3[0] = vi[0] * sa    +  vi[1] *  ca*cd  +  vi[2] * ( - ca * sd );
    vector3[1] = vi[0] * (-ca) +  vi[1] *  sa*cd  +  vi[2] * ( - sa * sd );
    vector3[2] =                  vi[1] *  sd     +  vi[2] * cd;

    return;
}







void Delays_Brut_fullgeometric(value_type * tis, const value_type& t0,
                 const long int& ntis, const long  int& nt0,
                 value_type ** sp, value_type **  si, value_type ** so, value_type ** SSB_to_PSB, value_type ** r_obs, value_type * roemer_ss,
                 const value_type& Mi, const value_type& Mo, const value_type& freq, const value_type distance, const value_type distance_derivative,
                 bool geometric, bool einstein, bool shapiro, bool aberration,
                 value_type delay[], bool & nanflag, value_type * spinaxis = NULL
                )
{
    /* SSB_to_PSB[i] = unit vector between the solar system barycenter (SSB) and the pulsar system barycenter (PSB) at time i
     * r_obs[i] = position of the earth geocenter at time i with respect to the SSB.
     * roemer_ss : array containing the component of roemer delays corresponding to the solar system as calculated by tempo (or anything else) to make the bats.
     * Of course sp, si, so , SSB_to_PSB and r_obs must be expressed in the same frame.
     * distance = distance to the pulsar at t0 (in meters)
     * distance_derivative = first derivative of "distance" , at t0, in m/s
     * spinaxis : unitary vector, rotation axis of the pulsar. If "NULL" (default), then it is taken parallel to the angular momentum of the pulsar.
     * nanflag : is set to True if at least one of the computed delays is NaN
     */

    long int i = 0;
    long int j = 0;
    //value_type delay[ntis] ;
    value_type * deindt ;

    value_type ein ;
    value_type nvp2 ;
    value_type nrpi ;
    value_type nrpo ;


    valarray<value_type> pe(3);
    valarray<value_type> nss(3);
    value_type pe_dot_nss = zero;
    value_type nsquare_pe = zero;
    value_type dist = zero;

    value_type distance_meters = distance * ( daysec * 365.25 *clight ) ; // convert light years into meters
    value_type distance_derivative_ms = distance_derivative * 1000.L  ; // convert from km/s to m / s


    // Compute the geometric delay
  if ( geometric == true ) {

    for (i=0; i < ntis ; ++i) {


        dist = distance_meters + (tis[i] - t0 ) * daysec * distance_derivative_ms; //norm3d( SSB_to_PSB[i] );
        for (j = 0; j < 3 ; j++)
        {
            pe[j] = sp[i][j] - r_obs[i][j]; // pulsar position - earth position with respect to their own barycenters
            nss[j] = SSB_to_PSB[i][j]; // / dist ;
        }

        pe_dot_nss = dotprod3d<value_type>(pe, nss) ;
        nsquare_pe = sumsquares3d<value_type>(pe); // square of the norm of pe
        delay[i] = pe_dot_nss  ;
        delay[i] += undemi *  nsquare_pe  / dist  ;
        delay[i] += - undemi * pow(pe_dot_nss,2)  / dist ;
        delay[i] /= daysec * clight ;
        delay[i] += (tis[i] - t0 ) * ( distance_derivative_ms / clight ) ;

        if (roemer_ss != NULL ) delay[i] += roemer_ss[i] / daysec;

    }
  }
  else {
    for (i=0; i < ntis ; ++i) {
        delay[i] = zero ;
    }
  }

      for (i = 0; i < ntis ; i ++)
          {
              if ( isnan(delay[i]) )
              {
                  cout << " isnan geo " << i << endl ;
                  nanflag = true;
                  break;
              }
          };

    // Compute Einstein delay
   if ( einstein == true ) {
       deindt = new value_type[ntis];

        for (i = 0 ; i < ntis ; ++i) {
            nrpi = pow( sp[i][0] - si[i][0], deux) + pow( sp[i][1] - si[i][1], deux)  + pow( sp[i][2] - si[i][2], deux) ;
            nrpi = sqrt( nrpi ) ;
            nrpo = pow( sp[i][0] - so[i][0], deux) + pow( sp[i][1] - so[i][1], deux)  + pow( sp[i][2] - so[i][2], deux) ;
            nrpo = sqrt( nrpo ) ;
            nvp2 = pow( sp[i][3], deux ) + pow( sp[i][4], deux ) + pow( sp[i][5], deux ) ;
            deindt[i] = Mi / nrpi ;
            deindt[i] += Mo / nrpo ;
            deindt[i] *= Msol * Ggrav / pow( clight, 2 ) ;
            deindt[i] += undemi * nvp2 / pow( clight, 2 ) ;
        };
 
        Spline splineein( tis, deindt, ntis ) ;

        ein = splineein.Integrate( t0, tis[0],
                                    nt0, nt0 +1,
                                    0, 1 ) ;
        delay[0] += ein ;

        for (i = 1 ; i < ntis - 1 ; ++i ) {// in range(1,ntoas) :
            ein += splineein.Integrate( tis[ i - 1 ] , tis[ i ],
                                        i - 1, i,
                                        i, i + 1 ) ;
            delay[i] += ein ;
        }

        // Deal with the particular case (for indices) of the last interpolation point
        ein += splineein.Integrate( tis[ ntis - 2 ] , tis[ ntis - 1 ],
                                    ntis - 2, ntis - 1,
                                    ntis - 2, ntis - 1 ) ;
        delay[ ntis - 1 ] += ein ;



          for (i = 0; i < ntis ; i ++)
          {
              if ( isnan(delay[i]) )
              {
                  cout << " isnan ein " << i << endl ;
                  nrpi = pow( sp[i][0] - si[i][0], deux) + pow( sp[i][1] - si[i][1], deux)  + pow( sp[i][2] - si[i][2], deux) ;
                nrpi = sqrt( nrpi ) ;
                nrpo = pow( sp[i][0] - so[i][0], deux) + pow( sp[i][1] - so[i][1], deux)  + pow( sp[i][2] - so[i][2], deux) ;
                nrpo = sqrt( nrpo ) ;
                nvp2 = pow( sp[i][3], deux ) + pow( sp[i][4], deux ) + pow( sp[i][5], deux ) ;
                deindt[i] = Mi / nrpi ;
                deindt[i] += Mo / nrpo ;
                deindt[i] *= Msol * Ggrav / pow( clight, 2 ) ;
                deindt[i] += undemi * nvp2 / pow( clight, 2 ) ;
                printf( " --> Rapport isnan ein : %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le   \n", nrpi, nrpo, nvp2, deindt[i],
                        splineein.Integrate( tis[ i - 1 ] , tis[ i ],
                                        i - 1, i,
                                        i, i + 1 ) ,
                        tis[i-1], tis[i] );
                  nanflag = true;
                  break;
              }
          };

        delete[] deindt;
   }

    // Compute Shapiro delay         (    For more ref : Backers and Hellings , Annual review of astronomy and astrophysics, 1986 )

    if (shapiro == true ) {
        value_type pi_dot_nss, po_dot_nss ;
        value_type Shap ;
        const value_type cst = 4.92521372097374e-06L ;            // = Ggrav * Msol / c**3 (unit = second)
        value_type PNParameter = 1.L ; // post newtonian parameter gamma, temporary


        for (i=0 ; i < ntis ; ++i ) {
            pi_dot_nss = ( sp[i][0] - si[i][0] ) * SSB_to_PSB[i][0] +
                         ( sp[i][1] - si[i][1] ) * SSB_to_PSB[i][1] +
                         ( sp[i][2] - si[i][2] ) * SSB_to_PSB[i][2] ;
            po_dot_nss = ( sp[i][0] - so[i][0] ) * SSB_to_PSB[i][0] +
                         ( sp[i][1] - so[i][1] ) * SSB_to_PSB[i][1] +
                         ( sp[i][2] - so[i][2] ) * SSB_to_PSB[i][2] ;
            nrpo = sqrt( pow( sp[i][0] - so[i][0], 2) + pow( sp[i][1] - so[i][1], 2) + pow( sp[i][2] - so[i][2], 2) ) ;
            nrpi = sqrt( pow( sp[i][0] - si[i][0], 2) + pow( sp[i][1] - si[i][1], 2) + pow( sp[i][2] - si[i][2], 2) ) ;
            Shap = Mi * log(  ( - pi_dot_nss + nrpi ) / clight ) ; // clight brings in a scaling parameter of one light second
            Shap += Mo * log( (nrpo - po_dot_nss) / clight ) ;
            Shap *= -(1. + PNParameter) * cst ;
            delay[i] += Shap / daysec ;
            if ( isnan(delay[i]) )
            {
                printf(" --> Rapport : %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le  %.19Le  \n", pi_dot_nss, po_dot_nss, nrpo, nrpi,
                       Mi * log(  ( - pi_dot_nss + nrpi ) / clight ),
                       Mo * log( (nrpo - po_dot_nss) / clight ) ,
                       Shap, delay[i] );
                nanflag = true;
                break;
            }
        }


          for (i = 0; i < ntis ; i ++)
          {
              if ( isnan(delay[i]) )
              {
                  cout << " isnan shap " << i << endl ;
                  pi_dot_nss = ( sp[i][0] - si[i][0] ) * SSB_to_PSB[i][0] +
                         ( sp[i][1] - si[i][1] ) * SSB_to_PSB[i][1] +
                         ( sp[i][2] - si[i][2] ) * SSB_to_PSB[i][2] ;
                  po_dot_nss = ( sp[i][0] - so[i][0] ) * SSB_to_PSB[i][0] +
                         ( sp[i][1] - so[i][1] ) * SSB_to_PSB[i][1] +
                         ( sp[i][2] - so[i][2] ) * SSB_to_PSB[i][2] ;
                  nrpo = sqrt( pow( sp[i][0] - so[i][0], 2) + pow( sp[i][1] - so[i][1], 2) + pow( sp[i][2] - so[i][2], 2) ) ;
                  nrpi = sqrt( pow( sp[i][0] - si[i][0], 2) + pow( sp[i][1] - si[i][1], 2) + pow( sp[i][2] - si[i][2], 2) ) ;
                  printf(" --> Rapport isnan shap : %.19Le %.19Le %.19Le %.19Le \n", pi_dot_nss, po_dot_nss, nrpo, nrpi);
                  nanflag = true;
                  break;
              }
          };

    }

//     // Compute aberration delay ( Aberration delay as in Smarr and Taylor 1976 or Damour & Deruelle 1986  )
    if ( aberration == true )    {

        value_type rotaxis[3] ; // Assume the spinning axis of the pulsar is aligned with its angular momentum.
        value_type nrotaxis ;
        value_type spin_cross_nss[3] ;
        value_type vp_dot_spin_cross_nss = 0.;

        if ( spinaxis == NULL )
        {
            rotaxis[0] = ( sp[0][1] * sp[0][5] - sp[0][2]*sp[0][4] ) ; // / pow( clight, 2 );
            rotaxis[1] = ( sp[0][2] * sp[0][3] - sp[0][0]*sp[0][5] ) ; // / pow( clight, 2 ) ;
            rotaxis[2] = ( sp[0][0] * sp[0][4] - sp[0][1]*sp[0][3] ) ; // / pow( clight, 2 ) ;
        }
        else
        {
            for (i=0 ; i<3 ; i++) rotaxis[i] = spinaxis[i];
        }
        nrotaxis = sqrt( pow( rotaxis[0], 2) + pow( rotaxis[1], 2) + pow( rotaxis[2], 2) ) ; 

        rotaxis[0] /= nrotaxis ;
        rotaxis[1] /= nrotaxis ;
        rotaxis[2] /= nrotaxis ;


        for (i= 0 ; i < ntis ; ++i ) {
            crossprod( rotaxis, SSB_to_PSB[i], spin_cross_nss ) ;
            vp_dot_spin_cross_nss = spin_cross_nss[0] * sp[i][3] +
                                   spin_cross_nss[1] * sp[i][4] +
                                   spin_cross_nss[2] * sp[i][5] ;
            delay[i] += vp_dot_spin_cross_nss / sumsquares3d(spin_cross_nss) / (deuxpi * freq ) / clight ; 

            if ( isnan(delay[i]) )
              {
                  cout << " isnan aberration " << i << endl ;
                  printf(" --> Rapport isnan aberration : %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le \n", spin_cross_nss[0], vp_dot_spin_cross_nss, rotaxis[0], rotaxis[1], rotaxis[2], nrotaxis);
                  nanflag = true;
                  break;
              }

        }


        //spline = cMySpline(self._interp_times, knots ) # UnivariateSpline(self._interp_times, knots, k=2, s=0.)

        //for i in range(nbts) :
          //  aberration[i] = spline( ts[i] ) / self._fitted_params[0]
    }

        return ;
}










void Delays_Brut_nogeometric(value_type * tis, const value_type& t0,
                 const long int& ntis, const long  int& nt0,
                 value_type ** sp, value_type **  si, value_type ** so, value_type * SSB_to_PSB,
                 const value_type& Mp, const value_type& Mi, const value_type& Mo, const value_type& freq,
                 bool einstein, bool shapiro, bool aberration,
                 value_type * delay,  const int truefreq, value_type & moydeindt, bool & nanflag,
                 value_type * spinaxis = NULL
                )
{
    /* t0 = time of reference from which the einstein delay is counted. ( this is a ssb time, for what it matters.. )
     * SSB_to_PSB = unit vector between the solar system barycenter (SSB) and the pulsar system barycenter (PSB)
     * Of course sp, si, so , SSB_to_PSB and r_obs must be expressed in the same frame.
     * spinaxis : unitary vector, rotation axis of the pulsar. If "NULL" (default), then it is taken parallel to the angular momentum of the pulsar.
     * nanflag : is set to True if at least one of the computed delays is NaN
     */

    long int i = 0;
    long int j = 0;
    value_type * deindt ;

    value_type ein ;
    value_type nvp2 ;
    value_type nrpi ;
    value_type nrpo ;


    valarray<value_type> pe(3);
    valarray<value_type> nss(3);
    value_type pe_dot_nss = zero;
    value_type nsquare_pe = zero;
    value_type dist = zero;

    nanflag = false; // set to true only if a delay is nan in the tests below

    // Compute Einstein delay
   if ( einstein == true ) {
       deindt = new value_type[ntis];
       
        for (i = 0 ; i < ntis ; ++i) {
            nrpi = pow( sp[i][0] - si[i][0], deux) + pow( sp[i][1] - si[i][1], deux)  + pow( sp[i][2] - si[i][2], deux) ;
            nrpi = sqrt( nrpi ) ;
            nrpo = pow( sp[i][0] - so[i][0], deux) + pow( sp[i][1] - so[i][1], deux)  + pow( sp[i][2] - so[i][2], deux) ;
            nrpo = sqrt( nrpo ) ;
            nvp2 = pow( sp[i][3], deux ) + pow( sp[i][4], deux ) + pow( sp[i][5], deux ) ;
            deindt[i] = Mi / nrpi ;
            deindt[i] += Mo / nrpo ;
            deindt[i] *= Msol * Ggrav / pow( clight, 2 ) ;
            deindt[i] += undemi * nvp2 / pow( clight, 2 ) ;
        };

        Spline splineein( tis, deindt, ntis ) ;
        ein = splineein.Integrate( t0, tis[0],
                                    nt0, nt0 +1,
                                    0, 1 ) ;
        delay[0] = ein ;

        for (i = 1 ; i < ntis - 1 ; ++i ) {
            ein += splineein.Integrate( tis[ i - 1 ] , tis[ i ],
                                        i - 1, i,
                                        i, i + 1 ) ;
            delay[i] = ein ;
        }

        // Deal with the particular case (for indices) of the last interpolation point
        ein += splineein.Integrate( tis[ ntis - 2 ] , tis[ ntis - 1 ],
                                    ntis - 2, ntis - 1,
                                    ntis - 2, ntis - 1 ) ;
        delay[ ntis - 1 ] = ein ;

          for (i = 0; i < ntis ; i ++)
          {
              if ( isnan(delay[i]) )
              {
                  cout << " isnan ein " << i << endl ;
                  nrpi = pow( sp[i][0] - si[i][0], deux) + pow( sp[i][1] - si[i][1], deux)  + pow( sp[i][2] - si[i][2], deux) ;
                nrpi = sqrt( nrpi ) ;
                nrpo = pow( sp[i][0] - so[i][0], deux) + pow( sp[i][1] - so[i][1], deux)  + pow( sp[i][2] - so[i][2], deux) ;
                nrpo = sqrt( nrpo ) ;
                nvp2 = pow( sp[i][3], deux ) + pow( sp[i][4], deux ) + pow( sp[i][5], deux ) ;
                deindt[i] = Mi / nrpi ;
                deindt[i] += Mo / nrpo ;
                deindt[i] *= Msol * Ggrav / pow( clight, 2 ) ;
                deindt[i] += undemi * nvp2 / pow( clight, 2 ) ;
                printf( " --> Rapport isnan ein : %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le   \n", nrpi, nrpo, nvp2, deindt[i],
                        splineein.Integrate( tis[ i - 1 ] , tis[ i ],
                                        i - 1, i,
                                        i, i + 1 ) ,
                        tis[i-1], tis[i] );
                nanflag = true;
                return;
              }
          };

        moydeindt = 0.;
        if (truefreq == 0 or truefreq == 2) { // Remove the constant component of the einstein delay
            for (i = 0; i < ntis ; i++) moydeindt += deindt[i];
            moydeindt /= ntis;
            for (i = 0; i < ntis ; i++) delay[i] -= moydeindt * (tis[i] - t0) ;
        }

        delete[] deindt;
   }
   else {
    for (i=0; i < ntis ; ++i) {
        delay[i] = zero ;
    }
  }

    // Compute Shapiro delay         (    For more ref : Backers and Hellings , Annual review of astronomy and astrophysics, 1986 )

    if (shapiro == true ) {

        value_type pi_dot_nss, po_dot_nss ;
        value_type Shap ;
        value_type PNParameter = 1.L ; // post newtonian parameter gamma, temporary
        const value_type cst = -(un + PNParameter) * 4.92521372097374e-06L / daysec ;            // = Ggrav * Msol / c**3 (unit = second)


        for (i=0 ; i < ntis ; ++i ) {
            pi_dot_nss = ( sp[i][0] - si[i][0] ) * SSB_to_PSB[0] +
                         ( sp[i][1] - si[i][1] ) * SSB_to_PSB[1] +
                         ( sp[i][2] - si[i][2] ) * SSB_to_PSB[2] ;
            po_dot_nss = ( sp[i][0] - so[i][0] ) * SSB_to_PSB[0] +
                         ( sp[i][1] - so[i][1] ) * SSB_to_PSB[1] +
                         ( sp[i][2] - so[i][2] ) * SSB_to_PSB[2] ;
            nrpo = sqrt( pow( sp[i][0] - so[i][0], 2) + pow( sp[i][1] - so[i][1], 2) + pow( sp[i][2] - so[i][2], 2) ) ;
            nrpi = sqrt( pow( sp[i][0] - si[i][0], 2) + pow( sp[i][1] - si[i][1], 2) + pow( sp[i][2] - si[i][2], 2) ) ;
            Shap = Mi * log(  ( - pi_dot_nss + nrpi ) / clight ) ; // clight brings in a scaling factor of one light second
            Shap += Mo * log( (nrpo - po_dot_nss) / clight ) ;
            Shap *=  cst ;
            delay[i] += Shap ;
            if ( isnan(delay[i]) )
            {
                printf(" --> Rapport isnan shap: %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le  %.19Le  \n", pi_dot_nss, po_dot_nss, nrpo, nrpi,
                       Mi * log(  ( - pi_dot_nss + nrpi ) / clight ),
                       Mo * log( (nrpo - po_dot_nss) / clight ) ,
                       Shap, delay[i] );
                nanflag = true;
                return;
            }

        }

    }

//     // Compute aberration delay ( Aberration delay as in Smarr and Taylor 1976 or Damour & Deruelle 1986  )
    if ( aberration == true )    {

        value_type rotaxis[3] ;
        value_type nrotaxis ;
        value_type spin_cross_nss[3] ;
        value_type vp_dot_spin_cross_nss = 0.;

        if ( spinaxis == NULL ) // Assume the spinning axis of the pulsar is aligned with the total orbital angular momentum.
        {
            rotaxis[0] = Mp*( sp[0][1] * sp[0][5] - sp[0][2]*sp[0][4] ) ;
            rotaxis[1] = Mp*( sp[0][2] * sp[0][3] - sp[0][0]*sp[0][5] ) ;
            rotaxis[2] = Mp*( sp[0][0] * sp[0][4] - sp[0][1]*sp[0][3] ) ;

            rotaxis[0] += Mi*( si[0][1] * si[0][5] - si[0][2]*si[0][4] ) ;
            rotaxis[1] += Mi*( si[0][2] * si[0][3] - si[0][0]*si[0][5] ) ;
            rotaxis[2] += Mi*( si[0][0] * si[0][4] - si[0][1]*si[0][3] ) ;

            rotaxis[0] += Mo*( so[0][1] * so[0][5] - so[0][2]*so[0][4] ) ;
            rotaxis[1] += Mo*( so[0][2] * so[0][3] - so[0][0]*so[0][5] ) ;
            rotaxis[2] += Mo*( so[0][0] * so[0][4] - so[0][1]*so[0][3] ) ;

            nrotaxis = sqrt( pow( rotaxis[0], 2) + pow( rotaxis[1], 2) + pow( rotaxis[2], 2) ) ;
            rotaxis[0] /= nrotaxis ;
            rotaxis[1] /= nrotaxis ;
            rotaxis[2] /= nrotaxis ;
        }
        else
        {
            for (i=0 ; i<3 ; i++) rotaxis[i] = spinaxis[i];
        }

        crossprod<value_type>( rotaxis, SSB_to_PSB, spin_cross_nss ) ;

        for (i= 0 ; i < ntis ; ++i ) {
            vp_dot_spin_cross_nss = spin_cross_nss[0] * sp[i][3] +
                                   spin_cross_nss[1] * sp[i][4] +
                                   spin_cross_nss[2] * sp[i][5] ;

            delay[i] += vp_dot_spin_cross_nss / sumsquares3d(spin_cross_nss) / (deuxpi * freq ) / clight ;

            if ( isnan(delay[i]) )
              {
                  cout << " isnan aberration " << i << endl ;
                  printf(" --> Rapport isnan aberration : %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le \n", spin_cross_nss[0], vp_dot_spin_cross_nss, rotaxis[0], rotaxis[1], rotaxis[2]);
                  nanflag = true;
                  return;
              }
        }
    }

        return ;
}





void Delays_Brut_geometric_local(value_type * tis, const value_type& currentBAT,  const long int& ntis,
                value_type ** sp, value_type **  si, value_type ** so, value_type * SSB_to_PSB, value_type * proper_motion,  value_type posepoch,
                 value_type * r_obs,
                 const value_type distance, const value_type distance_derivative,
                 value_type * delay, const bool kopeikin, const bool shklovskii, value_type test, value_type ** delay_details=NULL
                )
{
    /* tis[] = interpolation times for the pulsar motion (times of emission)
     * currentBAT = BAT (tempo2) for which the position of the observatory (r_obs + SSB_to_PSB ) are taken.
     * t0 = epoch of reference (in tempo2 ssb time in principle) for the position in the sky.
     * ntis = number of times in tis
     * sp, si , so = position and velocities (in m and m/s) of the pulsar, inner white dwarf and outer white dwarf at each time in "tis"
     * SSB_to_PSB = unit vector between the solar system barycenter (SSB) and the pulsar system barycenter (PSB) at a given BAT time (from tempo2)
     * proper_motion = proper motion in the plane of the sky  (unit = rad/year)
     * distance_derivative: proper motion along the radial axis (unit = mas/year). distance_derivative = (radial velocity/ distance)
     * r_obs = position of the observer with respect to the SSB  at a given time SAT (given by tempo2) corresponding to the above BAT.
     * Of course sp, si, so , SSB_to_PSB and r_obs must be expressed in the same frame.
     * distance = distance to the pulsar at BAT t0 (in ligth years)
     * shklovskii : if true include the shklovskii terms to the delay
     * test : if =1, prints a bunch of diagnostics
     * delay_details : if null, nothing happens. Otherwise returns the decomposition of "delays" into its components delay_details[:][i] = [roemer, shlovski, shklovski correction, kopeikin_orbital, kopeikin_annual_orbital, kopeikin_parallax_orbital]
     */

    long int i = 0;

    value_type dt = zero;

    value_type clightdays = clight * daysec;

    valarray<value_type> pm(proper_motion,3);
    pm = pm * distance; // unit : rad/ yr * (distance_meters / (clight * yrsec) = v/c
    valarray<value_type> nss(SSB_to_PSB,3);

    value_type distance_meters = distance * ( yrsec *clight ) ; // convert light years into meters

    // Initialize roemer
    value_type roemer = zero;


    // Initialize Proper-motion variables

    dt = (currentBAT - posepoch ) ; // in days

    value_type beta_para = (distance_derivative*radmasdeg) * distance ; // lyr / yr = d / (clight*yrsec) * (yrsec*d'/d) = d'/c with d' and d in m/s and m
    value_type k_para = beta_para * dt;// diplacement parallel to the direction of the pulsar at tpos divided by clight. In days
    valarray<value_type> beta_perp = pm ;
    value_type nbeta_perp = norm3d(beta_perp);// v_perp/c
    valarray<value_type> k_perp = beta_perp * dt;


    // Compute Shkloskii terms
    value_type shklo = dotprod3d<value_type>(k_perp,k_perp) / (deux * distance_meters)*clightdays;
    value_type shklo_cor = - k_para *clightdays / distance_meters * shklo; // corrective term to shklovskii

    // Initialize Kopeikin terms
    valarray<value_type> rp_perp(3);
    //value_type ro_perp = 0.L;
    valarray<value_type> ro_perp =  valarray<value_type>(r_obs,3) - dotprod3d<value_type>(r_obs, SSB_to_PSB) * nss ;
    value_type kopeikin_1 = zero;
    value_type kopeikin_2 = zero;
    value_type kopeikin_3 = zero;

    valarray<value_type> rp(3);

    for (i=0; i < ntis ; ++i) {
        rp = valarray<value_type>(sp[i],3);

    // Compute Roemer pulsar
        roemer = dotprod3d<value_type>(sp[i], SSB_to_PSB) / clightdays; // in days

        delay[i] = roemer;

    // Compute Kopeikin terms
        if (kopeikin == true)
        {
            rp_perp = rp - dotprod3d<value_type>(sp[i], SSB_to_PSB) * nss;
            kopeikin_1 = dotprod3d<value_type>(rp_perp,rp_perp) / (deux * distance_meters * clightdays);
            kopeikin_2 = - dotprod3d<value_type>(rp_perp , ro_perp) / (distance_meters * clightdays);
            kopeikin_3 = dotprod3d<value_type>(k_perp , rp_perp) / (distance_meters);

            delay[i] += kopeikin_1 + kopeikin_2 + kopeikin_3;
        }

    // Add Shklovskii terms
        if (shklovskii == true) delay[i] += shklo + shklo_cor;
        
    // Write delay components if requested :
        if (delay_details != NULL) 
        {
            delay_details[0][i] = roemer;
            delay_details[1][i] = shklo;
            delay_details[2][i] = shklo_cor;
            delay_details[3][i] = kopeikin_1;
            delay_details[4][i] = kopeikin_2;
            delay_details[5][i] = kopeikin_3;
        }
    }
     if (test == un) {
         printf("test geo pm : %.5Le  %.5Le  %.5Le ls %.5Le ls \n", beta_para, nbeta_perp,   k_para * daysec  , norm3d(k_perp)* daysec);
         printf("test geo roe : %.5Le s \n", roemer * 86400);
         printf("test kopeikin : %.5Le %.5Le %.5Le microsec \n ", kopeikin_1* 86400 *pow(10.,6), kopeikin_2* 86400 *pow(10.,6), kopeikin_3* 86400 *pow(10.,6));
         printf("test shklo : %.5Le %.5Le microsec \n", shklo* 86400 *pow(10.,6), shklo_cor* 86400 *pow(10.,6));
         printf("test dtpm : %.5Le %.5Le\n", dotprod3d<value_type>(k_perp  , ro_perp) / distance_meters * daysec * pow(10.,6), norm3d<value_type>(ro_perp));
         printf("dt ta tpos : %.5Le %.5Le %.5Le \n", dt, currentBAT, posepoch);
     }


        return ;
}
