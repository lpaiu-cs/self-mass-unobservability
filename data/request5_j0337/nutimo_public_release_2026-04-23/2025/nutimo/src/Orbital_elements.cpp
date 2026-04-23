// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
// SPDX-FileCopyrightText: 2024 Guillaume VOISIN, LUTH, Observatoire de Paris, PSL, CNRS (guillaume.voisin@obspm.fr astro.guillaume.voisin@gmail.com)
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */

#include "Constants.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include "Utilities.h"
#include "Diagnostics.h"
#include "Orbital_elements.h"

using namespace std ;

long double Solve_Kepler(long double dt, long double P, long double e, long double err/* = pow(10.L,-19)*/){
    /*
        Fixed point method to solve Kepler's equation : see Beutler's book
        dt : time elapsed since last pariastron passage
        P  : Period
        e  : eccentricity
        err : rounding error to achieve
    */
    long double  U = deuxpi * dt / P ;
    long double  dU = dix ;
    long double  sigma = zero ;
    int i = 0 ;

    while ( abs(dU) > err or i > 100 ) {
        sigma = 2.*pi*dt/P ;
        dU = (sigma - (U - e * sin(U) ) ) / ( un - e * cos(U) ) ;
        U += dU ;
        i += 1 ;
    }
    if (i == 101) cout << "Warning in 'Solve_Kepler' : maximum number of iterations reached. Last residual was : " << dU << endl ;

    return U ;
}


void orbel2statevect(value_type t, value_type e, value_type a1, value_type om1, value_type angli, value_type tperi1, value_type OrbPeriod, value_type om_an1 ,
                     value_type statevector1[6],
                     int jours/*=1*/, int solarframe/*=1*/) {
    /*
    Renvoie les positions et vitesses de l'objet 1 membre d'un système gravitationnel à deux corps de masses m1 et m0

        t       : temps (en jours, ou en s si jours=False)
        e       : excentricité
        a1      : demi grand axe de l'objet 1
        om1     : argument du périastre de l'objet 1, par rapport à la ligne de noeuds ascendants (rad)
        angli   : angle 0 >= (h,nss) >= -pi  avec h le moment angulaire du système, nss le vecteur pointant du barycentre du système solaire vers celui du système
        tperi1  : époque d'un passage au périastre de l'objet 1 (en jours, ou en s si jours=False)
        jours   : True si les temps sont spécifiés en jours, False s'ils le sont en secondes


    Renvoie la position et la vitesse à chaque temps en m et m/s, dans le référentiel (longitude ref x, y, nss) par défaut et dans (na,-n3,h/nh) si solarFrame=False.

        r1 = np.ndarray(len(t),3)   : position de 1 en m
        r1_t = np.ndarray(len(t),3) : vitesse de 1 en m/s

    Ressource : Beutler, Celestial Mechanics 1, $4.2.2
    */


    value_type et = e ;
    value_type a1t = a1 ;
    value_type om1t = om1 ;
    value_type anglit = angli ;
    value_type tperi1t = tperi1 ;

    value_type errKepler = pow(10.L, -18) ;

    value_type P;
    value_type dt ;
    value_type U = zero ;
    value_type sinv = zero ;
    value_type cosv = zero ;
    value_type v = zero ;
    value_type thet = zero ;
    value_type nr1 = zero ;
    value_type U_t = zero ;
    

    P = OrbPeriod ;  //   2.*pi*a**2*np.sqrt(1.-et**2) / h  # Période


    // Résout Kepler

    dt = fmod( abs( t-tperi1t ),  P ) * sgn( t-tperi1t ) ;
    U = Solve_Kepler(dt, P, et, errKepler)  ;

    sinv = sqrt( un - pow(et, 2) ) * sin(U) / ( un - et * cos(U) );
    cosv = ( cos(U) - et ) / ( un - et * cos(U) ) ;
    v = inversetrigo(cosv , sinv) ;
    // ! Test !
    //b = sqrt((un + et)/(un -et) ) ;
   // v = deux * atan(b * tan(U*undemi) ) ;
  //  cout << " test " <<  abs(v - deux * atan(b * tan(U*undemi) )) << endl;
    //
    thet = v + om1t ;
    nr1 = a1t * (un - et * cos( U ) ) ;

    value_type r1[3] = { cos( thet ) , sin( thet ) , zero } ;

    U_t = deuxpi * a1t / ( P * nr1 ) ;

    value_type r1_t[3] = { - sin( U ) , sqrt(un  - pow(et, 2) ) * cos( U ) , zero } ;

// // ! Test !
//     value_type nr1_t =  et *sin(U) ;//* U_t;
//
//     //cout << "test " << ar << "  " << U_t << "  " << er << "  " << et << "  " << etheta << endl ;
//     //printf("test %.19Le %.19Le %.19Le \n", J2 , K , eRR);
//     value_type v_t =  nr1 / a1t * b* deux / (deux*pow(cos(undemi*U),2) * (pow(b*tan(undemi*U),2) + un) ) ;//* U_t;
//     value_type r1_tbis[3] = { -v_t * sin( v ) + cos(v) *  nr1_t ,
//                             v_t * cos(v) + sin(v) * nr1_t,
//                             zero } ;
//     // printf(" test r1t , %.19Le, %.19Le, %.19Le \n",  r1_t[0], r1_t[1],   r1_t[2]) ;
//    //  printf(" test r1t , %.19Le, %.19Le, %.19Le \n", (r1_tbis[0] - r1_t[0]) /r1_t[0],(r1_tbis[1] - r1_t[1]) /r1_t[1],   (r1_tbis[2] - r1_t[2]) );///r1_t[2]) ;
// // fin test

     r1_t[0] *= U_t * a1t ;     r1_t[1] *= U_t * a1t ;     r1_t[2] *= U_t * a1t ;
     rot2D( om1t,r1_t , 0 , 1 ) ;       // Passe dans le repère (na,-n3,h/nh)

     if (solarframe == 1) {
        rot2D(anglit, r1_t, 1 , 2 ) ; // Passe dans le repère (na,n3',nss)
        rot2D(anglit,r1 , 1 , 2 ) ;        // Passe dans le repère (na,n3', nss)
        rot2D(om_an1, r1, 0 , 1 ) ;      // Passe dans le repère (x, y, nss) (ou x et y sont arbitraires)
        rot2D(om_an1, r1_t, 0 , 1) ;   // idem
        r1[0] *= nr1 ; r1[1] *= nr1 ; r1[2] *= nr1 ;
     }


    if (jours == 1 ) {
        r1_t[0] /= daysec ; r1_t[1] /= daysec ; r1_t[2] /= daysec ;
    }

    statevector1[0] = r1[0] ; statevector1[1] = r1[1] ; statevector1[2] = r1[2] ;
    statevector1[3] = r1_t[0] ; statevector1[4] = r1_t[1] ; statevector1[5] = r1_t[2] ;

    return ;

}



void orbel2statevect_1pn(value_type t, value_type et, value_type ar, value_type om1, value_type angli, value_type tperi1, value_type OrbPeriod, value_type om_an1 ,
                        value_type m1, value_type mc,
                     value_type statevector1[6],
                     int jours/*=1*/, int solarframe/*=1*/) {
    /*
     * Based on Damour & Deruelle 1985
     *
    Renvoie les positions et vitesses de l'objet 1 membre d'un système gravitationnel à deux corps de masses m1 et m0

        t       : temps (en jours, ou en s si jours=False)
        et      : excentricité "temporelle
        ar      : demi grand axe de l'objet 1 par rapport au barycentre 1PN
        om1     : longitude du périastre de l'objet 1, par rapport à la ligne de noeuds ascendants (rad)
        angli   : angle 0 >= (h,nss) >= -pi  avec h le moment angulaire du système, nss le vecteur pointant du barycentre du système solaire vers celui du système
        tperi1  : époque d'un passage au périastre de l'objet 1 (en jours, ou en s si jours=False)
        jours   : True si les temps sont spécifiés en jours, False s'ils le sont en secondes
        m1 : masse de l'objet 1 (en masse solaire)
        mc : masse du compagnon (masse 0)

    Renvoie la position et la vitesse à chaque temps en m et m/s, dans le référentiel  (longitude ref x, y, nss) par défaut et dans (na,-n3,h/nh)si solarFrame=False:
        nss est np.array([0,1,0])

        r1 = np.ndarray(len(t),3)   : position de 1 en m
        r1_t = np.ndarray(len(t),3) : vitesse de 1 en m/s

    Ressource : Damour & Deruelle 1985, Beutler, Celestial Mechanics 1, $4.2.2
    */

    const value_type GMsol = Ggrav  * Msol ;

    value_type Mtot = m1 + mc;  // after formula DD85 2.4b
    value_type nu = m1 * mc / pow(Mtot,2);
    value_type aRR = Mtot / mc * ar ; // DD85 6.3a
    value_type eRR = et * (un + GMsol * Mtot / (aRR*clight2) * (quatre -troisdemis * nu) ) ; // DD85 3.8b
    value_type er = eRR * ( un - GMsol * m1 * (m1 - mc) / (deux * Mtot * aRR * clight2) ) ; // DD85 6.3b
    value_type etheta = eRR * (un + GMsol * m1 * mc / (deux * Mtot * aRR * clight2 ) ) ; // DD85 4.13


    value_type anglit = angli ;

    value_type errKepler = pow(10.L, -18) ;

    value_type P;
    value_type dt ;
    value_type U = zero ;
    value_type v = zero ;
    value_type v_t = zero;
    value_type thet = zero ;
    value_type nr1 = zero ;
    value_type nr1_t = zero ;
    value_type U_t = zero ;
//     value_type b = zero;
    value_type J2 =  zero;
    value_type K = zero;



    //Vérifie que angli est dans la bonne plage ou tente une correction

//     if (anglit > zero and fmod(anglit, deuxpi) < pi)   anglit = -(fmod(anglit,deuxpi)) ;
//     if (anglit > zero and fmod(anglit, deuxpi) >= pi )  anglit = fmod(anglit,deuxpi) - deuxpi ;
//     if (anglit <= zero and fmod((-anglit),deuxpi) < pi ) anglit = -fmod((-anglit),deuxpi) ;
//     if (anglit <= zero and fmod((-anglit),deuxpi) >= pi ) anglit = -(-fmod((-anglit),deuxpi) + deuxpi) ;
//     if ( anglit != angli )  cout << "orbel2statevect : 0 >= angli >= -pi not satisfied. Now assuming angli =  " << anglit << ". Old angli : " << angli << endl;
//

    P = OrbPeriod ;  // Période


    // Résout Kepler

    dt = fmod( abs( t-tperi1 ),  P ) * sgn( t-tperi1 ) ;
    //dt = t - tperi1;
    U = Solve_Kepler(dt, P, et, errKepler)  ;

    // Compute Ae with formula DD86 17d
    value_type Ae = U + deuxpi * ( t-tperi1  - dt) / P ;
    value_type Aeold = zero;
    value_type Ae_inter = deux;
    value_type Aefactor = etheta / (un + sqrt(un - pow(etheta,2)) );
//    printf("\n test aez deuxpi %.19Le\n\n", ( t-tperi1  - dt) / P);
    int i = 1;
    do
    {
        Aeold = Ae;
        Ae_inter *= Aefactor;
        Ae += sin(i*U) *Ae_inter / i;
        i+=1;
    } while ( abs((Ae - Aeold)/ Ae) > errKepler) ;

//     //! Test !
//     value_type nn = 86400.*sqrt(GMsol*Mtot/pow(aRR,3)) * ( un + (m1*mc/pow(Mtot,2) - neuf) * GMsol *Mtot / (deux*aRR*clight2) );
//     printf("\n n test %.19Le, %.19Le, %.19Le \n\n", P,deuxpi / nn, P - deuxpi/nn );
//     // fin test

    //sinv = sqrt( un - pow(et, 2) ) * sin(U) / ( un - et * cos(U) );
    //cosv = ( cos(U) - et ) / ( un - et * cos(U) ) ;
    J2 = aRR * (un - pow(eRR,2) ) * ( GMsol * Mtot ) ; // DD85 4.15
    K = sqrt(J2 / ( J2 - six* pow(GMsol * Mtot / clight,2) ) ) ; //DD85 4.14
  //  K *= deux; // factor 2 in the function Ae(u)
//     b = sqrt((un + etheta)/(un -etheta) ) ;
    v = K * Ae;//* atan(b * tan(U*undemi) ) ;
    thet = v + om1 ;
    nr1 = ar * (un - er * cos( U ) ) ;
//             printf("\ntest ae %.19Le %.19Le %.19Le %i \n\n", tan(0.5*Ae), b * tan(U*undemi) , (tan(0.5*Ae) - b * tan(U*undemi))/(b * tan(U*undemi)), i );

    value_type r1[3] = { cos( thet ) , sin( thet ) , zero } ;

    U_t = deuxpi / ( P *  (un - et * cos( U ) )) ;
    nr1_t = ar * U_t * er *sin(U);

    //cout << "test " << ar << "  " << U_t << "  " << er << "  " << et << "  " << etheta << endl ;
    //printf("test %.19Le %.19Le %.19Le \n", J2 , K , eRR);

    // Compute derivative of Ae with respect to U (here variable Ae = dAe/dU) using DD86 17d
    i = 1;
    Ae_inter = deux;
    Ae = un;
    do
    {
        Aeold = Ae;
        Ae_inter *= Aefactor;
        Ae += cos(i*U) *Ae_inter ;
        i+=1;
    } while ( abs((Ae - Aeold)/ Ae) > errKepler) ;

    v_t = K* U_t * Ae * nr1 ; // derivative of v times nr1

    //v_t = U_t * nr1 * b*K / (deux*pow(cos(undemi*U),2) * (pow(b*tan(undemi*U),2) + un) ) ;
    value_type r1_t[3] = { -v_t * sin( v ) + cos(v) *  nr1_t ,
                            v_t * cos(v) + sin(v) * nr1_t,
                            zero } ;


  //   r1_t[0] *= U_t * a1t ;     r1_t[1] *= U_t * a1t ;     r1_t[2] *= U_t * a1t ;
     rot2D( om1,r1_t , 0 , 1 ) ;       // Passe dans le repère (na,-n3,h/nh)

 if (solarframe == 1) {
        rot2D(anglit, r1_t, 1 , 2 ) ; // Passe dans le repère (na,n3',nss)
        rot2D(anglit,r1 , 1 , 2 ) ;        // Passe dans le repère (na,n3', nss)
        rot2D(om_an1, r1, 0 , 1 ) ;      // Passe dans le repère (x, y, nss) (ou x et z sont arbitraires)
        rot2D(om_an1, r1_t, 0 , 1) ;   // idem
        r1[0] *= nr1 ; r1[1] *= nr1 ; r1[2] *= nr1 ;
     }


    if (jours == 1 ) {
        r1_t[0] /= daysec ; r1_t[1] /= daysec ; r1_t[2] /= daysec ;
    }

    statevector1[0] = r1[0] ; statevector1[1] = r1[1] ; statevector1[2] = r1[2] ;
    statevector1[3] = r1_t[0] ; statevector1[4] = r1_t[1] ; statevector1[5] = r1_t[2] ;

    return ;

}

// ############################################################################################################
void printorbels(orbel_t orbel)
{
    printf("a %.19Le\n", orbel.a);
    printf("ecc %.19Le\n", orbel.ecc);
    printf("Porb %.19Le\n", orbel.Porb);
    printf("norb %.19Le\n", orbel.norb);
    printf("i %.19Le\n", orbel.i);
    printf("Oman %.19Le\n", orbel.Oman);
    printf("omperi %.19Le\n", orbel.omperi);
    printf("tperi %.19Le\n", orbel.tperi);
    printf("tasc_approx %.19Le\n", orbel.tasc_approx);
    printf("m %.19Le\n", orbel.m);
    printf("v %.19Le\n", orbel.v);
    printf("E %.19Le\n", orbel.E);
    printf("h %.19Le %.19Le %.19Le\n", orbel.h[0], orbel.h[1], orbel.h[2]);
    printf("qlap %.19Le %.19Le %.19Le\n", orbel.qlap[0],orbel.qlap[1],orbel.qlap[2]);
    printf("easc %.19Le %.19Le %.19Le\n", orbel.easc[0],orbel.easc[1],orbel.easc[2]);
}
    
    
void statevect2orbel(valarray<value_type> statevector, value_type mu, value_type t, valarray<value_type> dir_obs, orbel_t & orbel)
// SI units assumed for all variables except otherwise mention
// statevector : [x, y, z, vx, vy, vz] (m and m/s)
// mu : G * Mtot
// t : time corresponding to the state vector (s)
// dir_obs : vector giving the direction from the observer to the system.
//           dir_obs defines the z' axis. Oman is defined wrt the x' and y' axis which are the projections 
//           of x and y along z'. 
// orbel : output of orbital elements and associated quantities 
//        - tperi is defined as the first periastron passage after time 0. tasc_approx defined wrt tperi
//        - Angles v, E, m are defined to be continuously varying from there value at t=0 (no phase jump)
//    
{
    valarray<value_type> kobs(dir_obs);
    kobs /= norm3d(kobs); // vector orthogonal to plane of sky
    valarray<value_type> r(statevector[slice(0,3,1)]);
    value_type nr = norm3d(r);
    valarray<value_type> rn = r / nr ; 
    valarray<value_type> v(statevector[slice(3,3,1)]) ; 
    valarray<value_type> h = crossprod(r,v) ; // Angular Momentum
    value_type nh = norm3d(h) ;
    valarray<value_type> rdot_x_h = crossprod(v,h) ;
    int Nper=0;
    
    valarray<value_type> qlap = rdot_x_h - mu * rn ;// Laplace vector
    value_type nqlap = norm3d(qlap);
    valarray<value_type> easc = crossprod(kobs, h); // unit vector of the line of ascending nodes
    easc /= norm3d(easc) ;

    // Defining the frame of the observer (xobs, yobs, kobs) by projection of x,y,z on the observer's plane
    // defined by dir_obs
    valarray<value_type> xobs(0., 3); // projection x base vector on the plane of the observer (or plane of reference)
    valarray<value_type> yobs(0., 3); // projection y base vector on the plane of the observer (or plane of reference)
    valarray<value_type> axis(0.,3); // Axis of rotation between coord frame (in which dir_obs is expressed) and observer frame (where zobs = dir_obs)
    valarray<value_type> zcoord(0.,3);
    zcoord[2] = 1.;
    xobs[0] = 1.;
    yobs[1] = 1.;
    axis = crossprod(zcoord, kobs);
    value_type sinthet = norm3d(axis);
    if (sinthet != 0.)
    {
        value_type costhet = dotprod3d(zcoord, kobs);
        value_type thet = atan2(sinthet, costhet); 
        axis /= sinthet; 
        xobs = rodrigues<value_type>(xobs, axis, thet); // Rotate to frame of obs
        yobs = rodrigues(yobs, axis, thet);
    }
//     zcoord = rodrigues(zcoord, axis, thet);
//     printf("check zcoord : %.5Le %.5Le %.5Le\n", zcoord[0] - kobs[0], zcoord[1] - kobs[1], zcoord[2] - kobs[2]);
    /*
    xobs[0] = 1.;
    xobs -= kobs * dotprod3d(kobs,xobs);
    xobs /= norm3d(xobs);
//     yobs[1] = 1.;
//     yobs -= kobs * dotprod3d(kobs,yobs);
//     yobs /= norm3d(yobs);
    yobs = crossprod(kobs, xobs);
    */
//     printf("sv2oe xobs %.5Le  %.5Le  %.5Le\n", xobs[0], xobs[1], xobs[2]);
//     printf("sv2oe yobs %.5Le  %.5Le  %.5Le\n", yobs[0], yobs[1], yobs[2]);
//     printf("sv2oe zobs %.5Le  %.5Le  %.5Le\n\n", kobs[0], kobs[1], kobs[2]);
    
    // Frame of the orbit (qlap, yorb, h)
    valarray<value_type> xorb = qlap / nqlap ;
    valarray<value_type> yorb = crossprod(h, qlap) ;
    yorb /= norm3d(yorb) ;
    
    // Longitude of ascending node
    orbel.Oman = atan2(dotprod3d(easc,yobs), dotprod3d(easc,xobs)); 
//     orbel.Oman = copysign(orbel.Oman, kobs[2]);  // if kobs[2]>0 then orientation of the frame (xobs, yobs, kobs) is the same as that of (x,y,z)

    value_type slr = nh*nh / mu ; // semi-latus rectum
    orbel.i = acos(dotprod3d(h,kobs) / nh) ; // inclination
    orbel.ecc = nqlap / mu ;                 // eccentricity
    orbel.a = slr/ (1. - pow(orbel.ecc,2)) ; // separation
    orbel.omperi = -atan2(dotprod3d(easc, yorb), dotprod3d(easc, xorb)) ; // argument of periastron, minus sign is because angle is oriented from easc to xorb
    orbel.Porb = deuxpi * sqrt(pow(orbel.a,3)/mu) ;  // Orbital period
    orbel.norb = deuxpi / orbel.Porb ;              // Mean motion
    
    orbel.v = atan2(dotprod3d(r, yorb), dotprod3d(r, xorb)); // True anomaly
    value_type sinE = dotprod3d(r, yorb) / (orbel.a * sqrt(1 - pow(orbel.ecc,2))) ;
    value_type cosE = dotprod3d(r, xorb)/ orbel.a + orbel.ecc ;
    orbel.E = atan2(sinE, cosE); // eccentric anomaly
    orbel.m = orbel.E - sinE * orbel.ecc; // Mean anomaly
    orbel.tperi = t - orbel.m / orbel.norb; // Time of passage at periastron
    Nper = floor(orbel.tperi / orbel.Porb) ;
    orbel.tperi -= Nper * orbel.Porb ; // Set tperi to be the first periastron passage after time 0 (assuming periodic orbit)
    orbel.E += Nper * deuxpi;
    orbel.m += Nper * deuxpi; 
    orbel.tasc_approx = orbel.tperi - orbel.omperi / orbel.norb; // Time of passage at ascending node to zeroth order in eccentricity
    
    for (int i =0 ; i < 3 ; i++)
    {
        orbel.h[i] = h[i];
        orbel.easc[i] = easc[i];
        orbel.qlap[i] = qlap[i];
    };
}


void statevect2orbel_nbody(int nbody, valarray<valarray<value_type>> statevectors, valarray<value_type> Ms, value_type t, valarray<value_type> dir_obs, orbel_t * orbels)
// SI units assumed except otherwise mention
// Ms : Masses in units of Msun
{
    value_type mu = GMsol * Ms[0];
    value_type Mtotp = 0.;
    value_type Mtot = 0.;
    valarray<value_type> s1(0., 6);  
    valarray<value_type> ds(0., 6);  
//     valarray<value_type> sv(0., 6);
    
    for (int nb= 1 ; nb < nbody ; nb++)
    {
        mu += GMsol * Ms[nb];
        Mtot += Ms[nb-1];
        s1 = (s1 * Mtotp + statevectors[nb-1] * Ms[nb-1])/Mtot;
        Mtotp = Mtot;
        ds = statevectors[nb] - s1;

        statevect2orbel(ds, mu, t, dir_obs, orbels[nb-1]);
//         orbel2statevects_bis(orbels[nb-1], dir_obs, sv, t);
        
//         printf("\n##### sv - ds %d : \n", nb) ;
//         printf("%.8Le  %.8Le  %.8Le  %.8Le %.8Le %.8Le \n", statevectors[nb][0], statevectors[nb][1],  statevectors[nb][2], statevectors[nb][3], statevectors[nb][4], statevectors[nb][5]);
//         printf("%.8Le  %.8Le  %.8Le  %.8Le %.8Le %.8Le \n", s1[0], s1[1],  s1[2], s1[3], s1[4], s1[5]);
//         Print_table((sv -ds)/sv);
    }
}


void orbel2statevects_bis(orbel_t orbel, valarray<value_type> dir_obs, valarray<value_type> & statevector, value_type & t)
// SI units assumed 
// Return the statevector corresponding to orbel in the coordinate system in which dir_obs is expressed. 
// The state vector corresponds is the effective one-body. 
{
    valarray<value_type> kobs(dir_obs);
    kobs /= norm3d(kobs); // vector orthogonal to plane of sky
    // Defining the frame of the observer (xobs, yobs, kobs) by projection of x,y,z on the observer's plane
    // defined by dir_obs
    valarray<value_type> xobs(0., 3); // projection x base vector on the plane of the observer (or plane of reference)
    valarray<value_type> yobs(0., 3); // projection y base vector on the plane of the observer (or plane of reference)
    valarray<value_type> axis(0.,3); // Axis of rotation between coord frame (in which dir_obs is expressed) and observer frame (where zobs = dir_obs)
    valarray<value_type> zcoord(0.,3);
    zcoord[2] = 1.;
    xobs[0] = 1.;
    yobs[1] = 1.;
    axis = crossprod(zcoord, kobs);
    value_type sinthet = norm3d(axis);
    if (sinthet != 0.)
    {
        value_type costhet = dotprod3d(zcoord, kobs);
        value_type thet = atan2(sinthet, costhet); 
        axis /= sinthet; 
        xobs = rodrigues<value_type>(xobs, axis, thet);
        yobs = rodrigues(yobs, axis, thet);
    }
//     zcoord = rodrigues(zcoord, axis, thet);
//     printf("check zcoord : %.5Le %.5Le %.5Le\n", zcoord[0] - kobs[0], zcoord[1] - kobs[1], zcoord[2] - kobs[2]);
    
    /*
    xobs[0] = 1.;
    xobs -= kobs * dotprod3d(kobs,xobs);
    xobs /= norm3d(xobs);
    yobs = crossprod(kobs, xobs);*/
//     yobs[1] = 1.;
//     yobs -= kobs * dotprod3d(kobs,yobs);
//     yobs /= norm3d(yobs);
//     
    // Matrix converting from observer's frame coordinates to coordinates in which dir_obs_is given
    // It is unity if dir_obs = [0,0,1]
//     for (j = 0 ; j < 3 ; j++)
//     {
//         R[j][0] = xobs[j]; 
//         R[j][1] = yobs[j]; 
//         R[j][2] = kobs[j]; 
//     }
    
    t = orbel.m / orbel.norb + orbel.tperi;
    
    //-------------------------------
    value_type et = orbel.ecc ;
    value_type a1t = orbel.a ;
    value_type om1t = orbel.omperi ;
    value_type anglit = orbel.i ;
    value_type tperi1t = orbel.tperi ;
    value_type om_an1 = orbel.Oman ;

    value_type errKepler = pow(10.L, -18) ;

    value_type P = orbel.Porb;
    value_type dt ;
    value_type U = zero ;
    value_type cosU, sinU;
    value_type nr1 = zero ;
    value_type U_t = zero ;


    // Résout Kepler

    dt = fmod( abs( t-tperi1t ),  P ) * sgn( t-tperi1t ) ;
    U = Solve_Kepler(dt, P, et, errKepler)  ;
    cosU = cos(U);
    sinU = sin(U);
//     printf("kepler : %.10Le  %.10Le  %.10Le\n", U - et*sinU - orbel.m, U, dt*deuxpi/P);
//     printf("kepler : %.10Le\n", U - et*sinU - dt*deuxpi/P);
    nr1 = a1t * (un - et * cosU ) ;
//     valarray<value_type> r1(0.,3);
//     valarray<value_type> v1(0.,3);
//     value_type se2 = sqrt(un  - et*et );
//     r1[0] =a1t*(cosU -et); 
//     r1[1] = sinU * se2;
//     v1[0] = -sinU * U_t *a1t;
//     v1[1] = cosU*se2 * U_t *a1t;
    
    value_type r1[] = { cosU - et , sqrt(un  - pow(et, 2) ) * sinU , zero } ;

    U_t = deuxpi * a1t / ( P * nr1 ) ;

    value_type r1_t[] = { - sinU , sqrt(un  - pow(et, 2) ) * cosU , zero } ;

    r1_t[0] *= U_t * a1t ;     r1_t[1] *= U_t * a1t ;     r1_t[2] *= U_t * a1t ;
    //----- TEST
//     valarray<value_type> r(r1,3);
//     valarray<value_type> rn = r / norm3d(r) ; 
//     r *= a1t;
//     valarray<value_type> v(r1_t,3) ; 
//     valarray<value_type> h = crossprod(r,v) ; // Angular Momentum
//     value_type nh = norm3d(h) ;
//     valarray<value_type> rdot_x_h = crossprod(v,h) ;
//     value_type mu = pow(a1t,3)*pow(deuxpi/P,2);
//     valarray<value_type> qlap = rdot_x_h - mu * rn ;// Laplace vector
//     value_type nqlap = norm3d(qlap);
//     printf("qlap : ");
//     for (int i =0; i < qlap.size(); i++) printf("%.5Le  ", rdot_x_h[i]);
//     printf("\n");
//     for (int i =0; i < qlap.size(); i++) printf("%.5Le  ", mu*rn[i]);
//     printf("\n");
//     for (int i =0; i < qlap.size(); i++) printf("%.5Le  ", qlap[i]);
//     printf("e = %.10Le ( vs %.10Le)\n", nqlap/mu, et);
    //---------------------------------------------------------------------------------
    rot2D( om1t,r1_t , 0 , 1 ) ;       // Passe dans le repère (na, n2,h/nh) where na is ascending node and h x na = n2
    rot2D( om1t,r1 , 0 , 1 ) ;
    rot2D(anglit, r1_t, 1 , 2 ) ; // Passe dans le repère (na,n3',nss)
    rot2D(anglit, r1 , 1 , 2 ) ;        // Passe dans le repère (na,n3', nss)
    rot2D(om_an1, r1, 0 , 1 ) ;      // Passe dans le repère (x, y, nss) (ou x et y sont arbitraires)
    rot2D(om_an1, r1_t, 0 , 1) ;   // idem
    r1[0] *= a1t ; r1[1] *= a1t ; r1[2] *= a1t ;
    //-------------------------------------------------------------------------------------
//     printf(" U dt %.10Le %.10Le\n", U, dt);
//     printf("\nsvv %.5Le  %.5Le  %.5Le %.5Le  %.5Le  %.5Le\n", sv[0], sv[1], sv[2], sv[3], sv[4], sv[5]);
//     printf("xobs %.5Le  %.5Le  %.5Le\n", xobs[0], xobs[1], xobs[2]);
//     printf("yobs %.5Le  %.5Le  %.5Le\n", yobs[0], yobs[1], yobs[2]);
//     printf("zobs %.5Le  %.5Le  %.5Le\n\n", kobs[0], kobs[1], kobs[2]);

//     printf("kobs2 %.5Le %.5Le\n", kobs[2], dotprod3d(crossprod(xobs, yobs), kobs));

    statevector.resize(6);
    statevector[slice(0,3,1)] = xobs * r1[0];//sv[0];
    statevector[slice(0,3,1)] += yobs * r1[1];// sv[1];
    statevector[slice(0,3,1)] += kobs * r1[2]; //sv[2];    
    statevector[slice(3,3,1)] = xobs * r1_t[0]; //sv[3];
    statevector[slice(3,3,1)] += yobs * r1_t[1]; //sv[4];
    statevector[slice(3,3,1)] += kobs * r1_t[2]; //sv[5];    
}


void orbel2statevects_bis_nbody(int nbody, orbel_t * orbels, valarray<value_type> Ms, valarray<value_type> dir_obs, bool com_1PN,valarray<valarray<value_type>> & statevectors, value_type &  t)
// SI units assumed except otherwise mention
// Ms : Masses in units of Msun
// com_1PN : True sets the 1PN center of mass to zero (instead of Newtonian by default). 
{
    int nb=0;
    valarray<value_type> Msum(0., nbody);
    valarray<value_type> sv(0.,6);
    
    Msum[0] = Ms[0]; 
    for (nb=1; nb < nbody ; nb++) 
    {
        Msum[nb] =  Msum[nb-1] + Ms[nb];
//         printf("Msum%d = %.5Le\n", nb,Msum[nb]);
    }
    
    statevectors.resize(nbody);
    statevectors[nbody -1] = valarray<value_type>(0.,6);
    
    for (nb= nbody-2 ; nb >= 0 ; nb--)
    {
        orbel2statevects_bis(orbels[nb], dir_obs, sv, t);
        
        statevectors[nb] = statevectors[nb+1] - sv * Ms[nb+1] / Msum[nb+1];
        statevectors[nb+1] = statevectors[nb+1] + sv * Msum[nb] / Msum[nb+1];
        
//         printf("%.5Le\n", t);
//         printf("%.5Le  %.5Le  %.5Le  %.5Le  %.5Le  %.5Le\n", sv[0], 
//                sv[1], sv[2], sv[3], 
//                sv[4], sv[5]);
//         printf("%.5Le  %.5Le  %.5Le  %.5Le  %.5Le  %.5Le\n", statevectors[nb][0], 
//                statevectors[nb+1][1], statevectors[nb+1][2], statevectors[nb+1][3], 
//                statevectors[nb+1][4], statevectors[nb+1][5]);
//         printf("%.5Le  %.5Le  %.5Le  %.5Le  %.5Le  %.5Le\n\n", statevectors[nb+1][0], 
//                statevectors[nb][1], statevectors[nb][2], statevectors[nb][3], 
//                statevectors[nb][4], statevectors[nb][5]);
        
        
    }
    // Compute CoM and Energy
    if (com_1PN == true )
    {
        valarray<value_type> fullsv(0., 6*nbody);
        valarray<value_type> com_pos(0.,3);
        valarray<value_type> com_vel(0.,3);
        value_type nrj;
        for (nb=0; nb < nbody; nb++)
        {
            fullsv[slice(3*nb, 3,1)] = statevectors[nb][slice(0,3,1)];
            fullsv[slice(3*nb + 3*nbody, 3,1)] = statevectors[nb][slice(3,3,1)];
        }
        
        IntegralePrems3_1PN( fullsv, Ms, 0., com_pos, com_vel, nrj) ;
        
        sv[slice(0,3,1)] = com_pos;
        sv[slice(3,3,1)] = com_vel / Msum[nbody-1];
        
        for (nb=0; nb < nbody; nb++) statevectors[nb] -= sv;
    }
    // Test
//     for (nb=0; nb < nbody; nb++)
//     {
//         fullsv[slice(3*nb, 3,1)] = statevectors[nb][slice(0,3,1)];
//         fullsv[slice(3*nb + 3*nbody, 3,1)] = statevectors[nb][slice(3,3,1)];
//     }
//     IntegralePrems3_1PN( fullsv, Ms, 0., com_pos, com_vel, nrj) ;
//     printf("com_pos1 %.10Le   %.10Le   %.10Le \n ", com_pos[0], com_pos[1], com_pos[2]);
//     printf("com_vel1 %.10Le   %.10Le   %.10Le \n ", com_vel[0], com_vel[1], com_vel[2]);
    
}
