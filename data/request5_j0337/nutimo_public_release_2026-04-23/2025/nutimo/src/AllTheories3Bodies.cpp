// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
#include "AllTheories3Bodies.h"
# include "Constants.h"
#include "Utilities.h"

    using namespace std;
    using namespace boost::numeric::odeint;
  //  using namespace IntegratorScope;

// test 
    int printrhs = 0;
    
// ***** Implementation Integrateur (newtonian-3-body integrator)
const int n3body =3; // Convenience variable, the code is written such that this is n3body=3

valarray<valarray<valarray<value_type> > > dx(3);
valarray<valarray<valarray<value_type> > > ns(3);
valarray<valarray<value_type> >  rs(3);
valarray<valarray<value_type> > vs(3);
valarray<value_type>  vvs(3);
valarray<value_type> inter(3);
valarray<value_type> onebodyinter(3);
valarray<value_type> inter1(3);
valarray<value_type> inter2(3);
valarray<value_type> inter3(3);
value_type sinter = zero;
value_type sinter1 = zero;
value_type interflt = 0.;
value_type rmb = 0.;
int n1,n2,n3;
valarray<value_type> x (18);
valarray<value_type> dxdt(0., 18);
value_type Q1Ad = 1.;

// Variables for extra bodies
int nbody_plus_extra ; // 3 bodies + number of extra bodies  

// Dimensional coefficients of the differential equation system
    value_type k0 =1. ;
    value_type k1= 1. ;
    value_type k2 = 1. ;
// Body masses (in solar masses)
    value_type M0, M1, M2;
    valarray<value_type> Ms(3); // = {M0, M1, M2}; // Masses
    

// General SEP violation parameters
//     value_type gammabar[3][3];
//     value_type betabar[3][3][3];
//     value_type Gg[3][3];
    valarray<valarray<value_type>> gammabar;
    valarray<valarray<valarray<value_type>>> betabar;
    valarray<valarray<value_type>> Gg;
    
// SEP violation parameters
    value_type SEP_D; // In newtonian : M0_gravitational = (1 + SEP_D) M0_inertial

// Quadrupole parameters
    vector<value_type> quadrupole(4,zero);
    bool quadrupole_enabled;

// // Physical constants
//     const value_type Ggrav = 6.67384e-11 ;      // m^3 s^-2 kg^-1
//     const value_type clight = 299792458. ;     // m/s

// Adimensioned constants (prior to integration thanks to parameters given by the user)
    value_type GgravAd = Ggrav ;
    value_type G2oC2 = pow( Ggrav/clight , 2) ;
    value_type GoC2 = Ggrav / pow(clight , 2) ;
    value_type clightAd = clight;
    value_type clightAd2 = clight * clight;

//  Rough estimate of the quadrupolar momentum of the 0.2Msol white dwarf, from K. Boshkayev, H. Quevedo and B. Zhami 2017
    const value_type Qwd = 7.395682251336515e+48 ; // kg.m^2 / (2pi/Spin Period)^2
    value_type QwdAd = Qwd; // QwdAd = dimensionless version of Qwd, done further in the code


// Natural coordinates base vectors :
    valarray<value_type> e_xnat(3);
    valarray<value_type> e_ynat(3);
    valarray<value_type> e_znat(3);

// Type of integrator to be used (value equal to Integrateur::integrator_type )
    int integrator_choice;

// RHS counter : number of times a rhs has been called
    int rhs_counter =0;

// ------------- Small utilities
    int mantissa = 15;
value_type trimfloat(value_type x, int mantissa)
{
    value_type test;
    int puissance;

    puissance = -50;
    test = abs(x)/pow(dix,puissance);
    while((test > 10. or test <= 1.) and test != 0.)
    {
        test = abs(x)/pow(dix,puissance);
        puissance++;
    }
    puissance--;
    test = x * pow(dix,mantissa - 1 - puissance);
    if (test < 0.) test = test + floor(abs(test));
    if (test > 0.) test = test - floor(test);
    test = x - test/ pow(dix,mantissa - 1 - puissance);
    return test;
}



void Integrateur::Initialise()
{
  // Initialise computing variaibles (used in rhs)
  for (int i = 0 ; i < n3body ; i++)
  {
  dx[i] =  valarray<valarray<value_type> >(3) ;
  ns[i] =  valarray<valarray<value_type> >(3) ;
  rs[i] = valarray<value_type> (3) ;
  vs[i] = valarray<value_type> (3) ;
      for (int j = 0 ; j < n3body ; j++)
      {
      dx[i][j] = valarray<value_type>(3) ;    // Differences
      ns[i][j] = valarray<value_type>(3) ;    // unit vectors
      }
  }
}



// Main function : integrate the and put the results in "states".
// Integrate both forwards and backwards from the inital condition depending depending on the times recquired by the user.
void Integrateur::Integrate_Allways()
{
    // Declarations
    intervals_type interv;
    const state_type xini = x0;
    const value_type tini = 0. ;
    const value_type dtini = dt0;
    state_type current_state = x0;
    int compteur=0;
    int groscompt = 0;
    const int nt = tpos.size();
//     value_type pasmin=1000. ; // ! Test !
    int i,j,k;
    int tn=0;

    valarray<value_type> xx(3);
    valarray<value_type> vv(3);


    printrhs = 0; // ! test !
    
    // Initialize state vectors
    states.resize( ts.size() );


    // Initialize constants used in the right-hand-side routine ("rhs" and "retro_rhs")
    M0 = int_M0;
    M1 = int_M1;
    M2 = int_M2;

    SEP_D = int_SEP_D;
    Gg = int_Gg;
    gammabar = int_gammabar;
    betabar = int_betabar;


// Defining natural coordinates :
    xx[0] = x0[0]; xx[1] = x0[1]; xx[2] = x0[2];
    vv[0] = x0[9]; vv[1] = x0[10]; vv[2] = x0[11];
    e_znat = crossprod(xx, vv) * M0;
    xx[0] = x0[3]; xx[1] = x0[4]; xx[2] = x0[5];
    vv[0] = x0[12]; vv[1] = x0[13]; vv[2] = x0[14];
    e_znat += crossprod(xx, vv) * M1;
    xx[0] = x0[6]; xx[1] = x0[7]; xx[2] = x0[8];
    vv[0] = x0[15]; vv[1] = x0[16]; vv[2] = x0[17];
    e_znat += crossprod(xx, vv) * M2;
    e_znat /= norm3d<value_type>(e_znat);

    e_ynat[0] = -e_znat[1];
    e_ynat[1] = e_znat[0] - e_znat[2];
    e_ynat[1] = e_znat[1];
    e_ynat /= norm3d<value_type>(e_ynat); // defined arbitrarily

    e_xnat = crossprod(e_ynat, e_znat); // defined arbitrarily, could be Laplace vector

// Uncomment to see what the natural basis is
    /*
    printf("\n Natural base unit vectors : \ne_x = %.4Le %.4Le %.4Le\ne_y = %.4Le %.4Le %.4Le\ne_z = %.4Le %.4Le %.4Le\n\n", e_xnat[0], e_xnat[1], e_xnat[2], e_ynat[0], e_ynat[1], e_ynat[2], e_znat[0], e_znat[1], e_znat[2]);
    */
    k0 = pow(length, 3) / pow(timescale, 2) / mass / Ggrav ;  // Proportionality constant in units of "length", "timescale" and "mass"
    G2oC2 = pow(Ggrav / clight, 2) * pow(timescale * mass, 2) /  pow( length, 4) ;
    GoC2 = Ggrav / pow(clight, 2) * mass / length ;
    GgravAd = Ggrav * pow(timescale, 2) * mass / pow(length, 3) ;
    clightAd = clight * timescale / length ;
    clightAd2 = clightAd * clightAd;
    integrator_choice = integrator_type;

    Ms[0] = M0 * GgravAd / clightAd2; Ms[1] = M1 * GgravAd / clightAd2; Ms[2] = M2 * GgravAd / clightAd2; // Takes into account that c = G = 1 in the expressions below
    
    
// Integrator specific initialsarions

  if (integrator_type == 10)
  {
      quadrupole = int_quadrupole;
      quadrupole[2] *= daysec/timescale;
      quadrupole[3] *= daysec/timescale;
      quadrupole[0] *= clight * clight / length / length ; // convert Msol.ls^2 -> Msol.length^2
  //     printf("\nquadrupole %.5Le %.5Le %.5Le %.5Le \n\n",quadrupole[0], quadrupole[1], quadrupole[2], quadrupole[3] );
  }
  if (integrator_type == 1)
  {
//       if (int_quadrupole.size() > 0)
//           quadrupole = int_quadrupole;
//       else
//           quadrupole[0] = zero;
//       QwdAd = Qwd / (mass * length * length * timescale * timescale);
//       quadrupole[0] *= timescale * timescale;
// 
//       if ( abs(quadrupole[0]) > zero )
//           quadrupole_enabled = true;
//       else
//           quadrupole_enabled = false;
      quadrupole[0] = zero;
      quadrupole_enabled = false;
      
  }

  if (integrator_type == 2)
  {
    quadrupole[0] = int_quadrupole[0] / (mass * length * length);
    Q1Ad =  GgravAd *  neufdemis  * quadrupole[0] /  Ms[1] ;
  }
  
  if (integrator_type ==3 or integrator_type ==30) // Core triple system at 0PN or 1PN + extra bodies also at 1PN
  {
      nbody_plus_extra = n3body + int_M_extra.size();
      Ms.resize(nbody_plus_extra);
      rs.resize(nbody_plus_extra);
      vvs.resize(nbody_plus_extra);
      vs.resize(nbody_plus_extra);
      dx.resize(nbody_plus_extra);
      ns.resize(nbody_plus_extra);
      for (i = 0; i < nbody_plus_extra ; i++) 
      {
          rs[i].resize(nbody_plus_extra);
          ns[i].resize(nbody_plus_extra);
          dx[i].resize(nbody_plus_extra);
          vs[i].resize(3);
          for (j = 0; j < nbody_plus_extra ; j++) 
          {
              ns[i][j].resize(3);
              dx[i][j].resize(3);
          }
      }
      Ms[0] = M0 * GgravAd / clightAd2; Ms[1] = M1 * GgravAd / clightAd2; Ms[2] = M2 * GgravAd / clightAd2; // Takes into account that c = G = 1 in the expressions below
      for (i = n3body; i < nbody_plus_extra ; i++)  
      {
        Ms[i] = int_M_extra[i - n3body] * GgravAd/ clightAd2;
      }
      x.resize(nbody_plus_extra * 6);
      dxdt.resize(nbody_plus_extra*6, 0.);
      
      for (i = 0; i < ts.size(); i++) states[i].resize((nbody_plus_extra)*6);
           
  }


// Do the backward integration
    Retro_Integre();

// Now integrate forward

    // Initialize the Boost stepper.
     Stepper stepper(tolint, tolint, 1. ,1., true);

    stepper.initialize(xini , tini, dtini );

    if (tini == tpos[0]) // & tn == 0
    {
        states[ntneg] = x0 ;
        tn=1; // tn++
    }

    while (tn < nt) //for (tn = 0; tn < nt; tn=tn)
    {

        if (compteur > 10000)
        {
            cout << "Warning in Integrateur::Integre_Allways() : more than 10000 iteration to reach the next time step ! " ;
            cout << " Time number is : " << tn << " in " << nt << endl;
            compteur = 0;
            groscompt++;
            if (groscompt > 10)
            {
                cout << "Error in Integrateur::Integre_Allways() : more than 100 000 iteration to reach the next time step ! Breaking the loop.." << endl;
                cout << " Time number is : " << tn << " in " << nt << " and time = " << tpos[tn] << endl;
                cout << "Initial time step was : " << dt0 << endl;
                cout << "Last interpolation interval is : " << interv.first << "  " << interv.second << endl;
                break;
            }
        }

        // Do the actual step forward
        switch (integrator_type) {
            case 1:
              interv = stepper.do_step(rhs_GR) ;
              break;
            case 2:
              interv = stepper.do_step(rhs_GR_quadrupole);
              break;
            case 3:
              interv = stepper.do_step(rhs_GR_nbody); 
              break;
            case 30:
              interv = stepper.do_step(rhs_0PN_nbody); 
              break;
            case -1:
              interv = stepper.do_step( rhs_biKeplerian ) ;
              break;
            case 0:
              interv = stepper.do_step(rhs_Newt) ;
              break;
            case 10:
              interv = stepper.do_step(rhs_NewtQuad) ;
              break;
            default :
              cout << "Integrator type does not match anything known ! Using Newtonian." << "\n";
        }

        // Uncomment to print the width of each step
//          cout << "--> dif : " << interv.second - interv.first << endl;
//         if (interv.second - interv.first < pasmin) //  ! Test !
//         {
//             pasmin = interv.second - interv.first ;
//         }

        while ( interv.second > tpos[tn] and tn < nt) // interv.first <= tpos[tn] and
        {
            stepper.calc_state( tpos[tn] , current_state) ;
            states[ ntneg + tn ] = current_state ;
            compteur=0;
            tn++;
        }
        compteur++;
    }
    // Uncomment to print the minimum width between two steps
    //cout << "--------> pasmin : " << pasmin << endl;
    return;
}



//  Do the backward integration and store the result in the "ntneg" first terms of "states". Called by "Integrateur::Integrate_Allways()"
void Integrateur::Retro_Integre()
{
    // Declarations
    intervals_type interv;
    const state_type xini = x0;
    const value_type tini = 0. ;
    const value_type dtini = dt0;
    state_type current_state = x0;
    int compteur=0;
    int groscompt = 0;
    const int nt = tneg.size();


    // Initialisation du stepper et de la récurrence
    Stepper stepper(tolint, tolint, 1. ,1., true);
    stepper.initialize(xini , tini, dtini );
    int tn = nt -1;
    while (tn >= 0) //for (int tn = nt - 1; tn >= 0; tn=tn)
    {
        if (compteur > 10000)
        {
            cout << "Warning in Integrateur::Retro_Integre() : more than 10000 iteration to reach the next time step ! " ;
            cout << " Time number is : " << tn << " in " << nt << endl;
            cout << "Current initial time: " << tini *timescale << "| Current time target: " << (tini+dtini)*timescale << " timescale "<< timescale << endl;
            compteur = 0;
            groscompt++;
            if (groscompt > 10)
            {
                cout << "Error in Integrateur::Retro_Integre : more than 100 000 iteration to reach the next time step ! Breaking the loop.." << endl;
                cout << " Time number is : " << tn << " in " << nt << " and time = " << tneg[tn] << endl;
                cout << "Initial time step was : " << dt0 << endl;
                cout << "Last interpolation interval is : " << interv.first << "  " << interv.second << endl;
                break;
            }
        }
        // Do the actual step
        switch (integrator_type) {
            case 1:
              interv = stepper.do_step(retro_rhs_GR) ;
              break;
            case 2:
              interv = stepper.do_step(retro_rhs_GR_quadrupole);
              break;
            case 3:
              interv = stepper.do_step(retro_rhs_GR_nbody); 
              break;
            case 30:
              interv = stepper.do_step(retro_rhs_0PN_nbody); 
              break;
            case -1:
              interv = stepper.do_step( retro_rhs_biKeplerian ) ;
              break;
            case 0:
              interv = stepper.do_step(retro_rhs_Newt) ;
              break;
            case 10:
              interv = stepper.do_step(retro_rhs_NewtQuad) ;
              break;
            default :
              cout << "Integrator type does not match anything known ! Using Newtonian." << "\n";
        }

        while ( interv.second > tneg[tn] and tn >= 0) // interv.first <= tneg[tn]  and
        {
            stepper.calc_state( tneg[tn] , current_state) ;
            states[tn ] = current_state;
            compteur=0;
            tn--;
        }
        compteur++;
    }

    return;


}

// Right-hand side for forward integration
// This routine is an interface that chooses which right-hand-side to use depending on the theory defined by the user.
void Integrateur::rhs(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    rhs_GR(xvect, dxdtvect, t) ;
}

// Right-hand side for backward integration
// This routine is an interface that chooses which right-hand-side to use depending on the theory defined by the user.
void Integrateur::retro_rhs(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    rhs(xvect, dxdtvect, t);
    for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
    {
        dxdtvect[n1] = - dxdtvect[n1]; // Only change from rhs here !!
    }
}



// RHS for 3-body integration at first post-newtonian order
void Integrateur::rhs_GR(const state_type& xvect, state_type& dxdtvect, const value_type& t) // Version with equations from equation from Soffel, Relativity in astrometry, Celestial mechanics and Geodesy 1989
{

    for ( n1=0 ; n1 < xvect.size() ; n1++)
    {
        x[n1] = xvect[n1];
    }
/* Initialize matrixes for the differences between the locations of the bodies, two by two ;
   the distances between two bodies; and unit vectors. */
    for ( n1 = 0 ; n1 < 3 ; n1++)
    {
    vs[n1] = x[slice(9 + n1 * 3, 3, 1)] ;
    vs[n1] /= clightAd;    // Takes into account that c = G = 1 in the expressions below
    inter = pow(vs[n1], deux);
    vvs[n1] = inter.sum() ;
        for ( n2 = 0 ; n2 < 3 ; n2++)
        {
        dx[n1][n2] = x[slice(n2 * 3,3,1)];     // !! Different convention from Will's !
        dx[n1][n2] -= x[slice(n1 * 3,3,1)];

        inter = pow(dx[n1][n2], deux);
        rs[n1][n2] = sqrt ( inter.sum() );    // Distances
        ns[n1][n2] = dx[n1][n2];
        ns[n1][n2] /= rs[n1][n2];
        }
    }


    dxdt[slice(0,9,1)] = x[slice(9,9,1)]; // Set first derivatives

// Computation of the acceleration :

    for (n1 = 0 ; n1 < n3body ; n1++)
     {
        onebodyinter =zero;

        for (n2 = 0 ; n2 < n3body ; n2++)
        {
          inter1 = zero;
          if (n2 != n1)
          {
            // 2 body terms in the direction of ns[n1][n2]
              sinter = un;
              sinter1 = dotprod3d<value_type>(vs[n1], vs[n2]);
              sinter -= (quatre + deux *  gammabar[n1][n2])*sinter1 ;
              sinter += vvs[n1] * ( un  + gammabar[n1][n2]);
              sinter += (deux + gammabar[n1][n2]) * vvs[n2];
              sinter1 = dotprod3d<value_type>(vs[n2], ns[n1][n2]);
              sinter -= troisdemis * sinter1 * sinter1;
              inter1 += ns[n1][n2] * sinter;

            // 2-body terms in the direction of (vs[n2] - vs[n1])
              inter = (vs[n2] - vs[n1]);
              inter *= - deux * gammabar[n1][n2];
              inter +=  quatre * vs[n1] - trois * vs[n2];
              sinter = dotprod3d<value_type>(ns[n1][n2], inter);//quatre * vs[n1] - trois * vs[n2] - deux * gammabar[n1][n2] * (vs[n2] - vs[n1]));
              inter1 += (vs[n2] - vs[n1]) * sinter ;
            // Prefactor 2-body
              inter1 *= Gg[n1][n2] * Ms[n2] / (rs[n1][n2]*rs[n1][n2]) ;

            // 3-body terms
            inter2 = zero;
            inter3 = zero;
            for (n3 = 0 ; n3 < n3body ; n3++)
            {
            // 3-body terms with n3 != n2
                if (n3 != n2)
                {
                    // terms in the direction of ns[n1][n2]
                    sinter = undemi * dotprod3d<value_type>(ns[n1][n2], ns[n2][n3]) / rs[n2][n3] ;
                    sinter -= (un + betabar[n2][n3][n1] + betabar[n2][n1][n3]) / rs[n1][n2];
                    inter = sinter * ns[n1][n2];

                    // terms in the direction of ns[n2][n3]
                    sinter = (septdemis + deux * gammabar[n1][n2]) / rs[n2][n3] ;
                    inter += sinter * ns[n2][n3];

                    // prefactor
                    inter *=  Gg[n1][n2] * Gg[n2][n3] * Ms[n2] * Ms[n3] / (rs[n1][n2]*rs[n2][n3]);

                    inter2 += inter;
                }
             // 3-body terms with n3 != n1
                if (n3 != n1)
                {
                    sinter = (quatre + deux * gammabar[n1][n3] + betabar[n1][n2][n3] + betabar[n1][n3][n2]) ;
                    sinter *= Gg[n1][n2] * Gg[n1][n3] * Ms[n2] * Ms[n3] / (rs[n1][n2]*rs[n1][n2]*rs[n1][n3]);
                    inter3 -= ns[n1][n2] * sinter ;
                }
            }
            onebodyinter += inter1 + inter2 + inter3;
          }
        }
        onebodyinter *= clightAd2;
        dxdt[std::slice(9 + n1 * 3,3,1)] = onebodyinter;
    }

    for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
    {
        dxdtvect[n1] = dxdt[n1];
    }
 };

 void Integrateur::retro_rhs_GR(const state_type& xvect, state_type& dxdtvect, const value_type& t)
 {
     rhs_GR(xvect, dxdtvect, t);
     for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
     {
         dxdtvect[n1] = - dxdtvect[n1]; // Only change from rhs here !!
     }
 }


 

// RHS for n-body integration at first post-newtonian order 
void Integrateur::rhs_GR_nbody(const state_type& xvect, state_type& dxdtvect, const value_type& t) // Version with equations from equation from Soffel, Relativity in astrometry, Celestial mechanics and Geodesy 1989
{
    for ( n1=0 ; n1 < xvect.size() ; n1++)
    {
        x[n1] = xvect[n1];
    }
/* Initialize matrixes for the differences between the locations of the bodies, two by two ;
   the distances between two bodies; and unit vectors. */
    for ( n1 = 0 ; n1 < nbody_plus_extra ; n1++)
    {
    vs[n1] = x[slice(3*nbody_plus_extra + n1 * 3, 3, 1)] ;
    vs[n1] /= clightAd;    // Takes into account that c = G = 1 in the expressions below
    inter = pow(vs[n1], deux);
    vvs[n1] = inter.sum() ;
        for ( n2 = 0 ; n2 < n1 ; n2++)
        {
            dx[n1][n2] = x[slice(n2 * 3,3,1)];     // !! Different convention from Will's !
            dx[n1][n2] -= x[slice(n1 * 3,3,1)];
            
            inter = pow(dx[n1][n2], deux);
            rs[n1][n2] = sqrt ( inter.sum() );    // Distances
            ns[n1][n2] = dx[n1][n2];
            ns[n1][n2] /= rs[n1][n2];
            
            // Symetric components
            dx[n2][n1] = -dx[n1][n2];
            rs[n2][n1] = rs[n1][n2];
            ns[n2][n1] = -ns[n1][n2];
        }
    }
    

    dxdt[slice(0,3*(nbody_plus_extra),1)] = x[slice(3*(nbody_plus_extra),3*(nbody_plus_extra),1)]; // Set first derivatives

// Computation of the acceleration :

    for (n1 = 0 ; n1 < nbody_plus_extra ; n1++)
     {
        onebodyinter =zero;

        for (n2 = 0 ; n2 < nbody_plus_extra ; n2++)
        {
          inter1 = zero;
          if (n2 != n1)
          {
            // 2 body terms in the direction of ns[n1][n2]
              sinter = un;
              sinter1 = dotprod3d<value_type>(vs[n1], vs[n2]);
              sinter -= (quatre + deux *  gammabar[n1][n2])*sinter1 ;
              sinter += vvs[n1] * ( un  + gammabar[n1][n2]);
              sinter += (deux + gammabar[n1][n2]) * vvs[n2];
              sinter1 = dotprod3d<value_type>(vs[n2], ns[n1][n2]);
              sinter -= troisdemis * sinter1 * sinter1;
              inter1 += ns[n1][n2] * sinter;

            // 2-body terms in the direction of (vs[n2] - vs[n1])
              inter = (vs[n2] - vs[n1]);
              inter *= - deux * gammabar[n1][n2];
              inter +=  quatre * vs[n1] - trois * vs[n2];
              sinter = dotprod3d<value_type>(ns[n1][n2], inter);//quatre * vs[n1] - trois * vs[n2] - deux * gammabar[n1][n2] * (vs[n2] - vs[n1]));
              inter1 += (vs[n2] - vs[n1]) * sinter ;
            // Prefactor 2-body
              inter1 *= Gg[n1][n2] * Ms[n2] / (rs[n1][n2]*rs[n1][n2]) ;

            // 3-body terms
            inter2 = zero;
            inter3 = zero;
            for (n3 = 0 ; n3 < nbody_plus_extra ; n3++)
            {
            // 3-body terms with n3 != n2
                if (n3 != n2)
                {
                    // terms in the direction of ns[n1][n2]
                    sinter = undemi * dotprod3d<value_type>(ns[n1][n2], ns[n2][n3]) / rs[n2][n3] ;
                    sinter -= (un + betabar[n2][n3][n1] + betabar[n2][n1][n3]) / rs[n1][n2];
                    inter = sinter * ns[n1][n2];

                    // terms in the direction of ns[n2][n3]dxdt
                    sinter = (septdemis + deux * gammabar[n1][n2]) / rs[n2][n3] ;
                    inter += sinter * ns[n2][n3];

                    // prefactor
                    inter *=  Gg[n1][n2] * Gg[n2][n3] * Ms[n2] * Ms[n3] / (rs[n1][n2]*rs[n2][n3]);

                    inter2 += inter;
                }
             // 3-body terms with n3 != n1
                if (n3 != n1)
                {
                    sinter = (quatre + deux * gammabar[n1][n3] + betabar[n1][n2][n3] + betabar[n1][n3][n2]) ;
                    sinter *= Gg[n1][n2] * Gg[n1][n3] * Ms[n2] * Ms[n3] / (rs[n1][n2]*rs[n1][n2]*rs[n1][n3]);
                    inter3 -= ns[n1][n2] * sinter ;
                }
            }
            onebodyinter += inter1 + inter2 + inter3;
//              if (printrhs < 10)
//             {
//                 printf("n1=%d n2=%d Ms %.2Le rs %.2Le ns %.2Le %.2Le %.2Le \n", n1, n2, Ms[n2], rs[n1][n2],ns[n1][n2][0], ns[n1][n2][1], ns[n1][n2][2]);
//             }
          }
        }
        onebodyinter *= clightAd2;
        dxdt[std::slice(3*(nbody_plus_extra) + n1 * 3,3,1)] = onebodyinter;
    }


    for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
    {
        dxdtvect[n1] = dxdt[n1];
    }

 };

 void Integrateur::retro_rhs_GR_nbody(const state_type& xvect, state_type& dxdtvect, const value_type& t)
 {
     rhs_GR_nbody(xvect, dxdtvect, t);
     for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
     {
         dxdtvect[n1] = - dxdtvect[n1]; // Only change from rhs here !!
     }
 }


 
// RHS for n-body integration at Newtonian order 
void Integrateur::rhs_0PN_nbody(const state_type& xvect, state_type& dxdtvect, const value_type& t) 
{
    for ( n1=0 ; n1 < xvect.size() ; n1++)
    {
        x[n1] = xvect[n1];
    }
/* Initialize matrixes for the differences between the locations of the bodies, two by two ;
   the distances between two bodies; and unit vectors. */
    for ( n1 = 0 ; n1 < nbody_plus_extra ; n1++)
    {
    vs[n1] = x[slice(3*nbody_plus_extra + n1 * 3, 3, 1)] ;
    vs[n1] /= clightAd;    // Takes into account that c = G = 1 in the expressions below
    inter = pow(vs[n1], deux);
    vvs[n1] = inter.sum() ;
        for ( n2 = 0 ; n2 < n1 ; n2++)
        {
            dx[n1][n2] = x[slice(n2 * 3,3,1)];     // !! Different convention from Will's !
            dx[n1][n2] -= x[slice(n1 * 3,3,1)];
            
            inter = pow(dx[n1][n2], deux);
            rs[n1][n2] = sqrt ( inter.sum() );    // Distances
            ns[n1][n2] = dx[n1][n2];
            ns[n1][n2] /= rs[n1][n2];
            
            // Symetric components
            dx[n2][n1] = -dx[n1][n2];
            rs[n2][n1] = rs[n1][n2];
            ns[n2][n1] = -ns[n1][n2];
        }
    }
    

    dxdt[slice(0,3*(nbody_plus_extra),1)] = x[slice(3*(nbody_plus_extra),3*(nbody_plus_extra),1)]; // Set first derivatives

// Computation of the acceleration :

    for (n1 = 0 ; n1 < nbody_plus_extra ; n1++)
     {
        onebodyinter =zero;

        for (n2 = 0 ; n2 < nbody_plus_extra ; n2++)
        {
          inter1 = zero;
          if (n2 != n1)
          {
            // 2 body terms in the direction of ns[n1][n2]
            inter1 += ns[n1][n2];

            // Prefactor 2-body
            inter1 *= Gg[n1][n2] * Ms[n2] / (rs[n1][n2]*rs[n1][n2]) ;

            onebodyinter += inter1 ;
          }
        }
        onebodyinter *= clightAd2;
        dxdt[std::slice(3*(nbody_plus_extra) + n1 * 3,3,1)] = onebodyinter;
    }


    for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
    {
        dxdtvect[n1] = dxdt[n1];
    }

 };

 void Integrateur::retro_rhs_0PN_nbody(const state_type& xvect, state_type& dxdtvect, const value_type& t)
 {
     rhs_0PN_nbody(xvect, dxdtvect, t);
     for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
     {
         dxdtvect[n1] = - dxdtvect[n1]; // Only change from rhs here !!
     }
 }
 

 // RHS for 3-body integration at first post-newtonian order
 void Integrateur::rhs_GR_quadrupole(const state_type& xvect, state_type& dxdtvect, const value_type& t) // Version with equations from equation from Soffel, Relativity in astrometry, Celestial mechanics and Geodesy 1989
 {

 //    // Natural coordinates ( ie z = total angular momentum, defined in Integrate_Allways)

     for ( n1=0 ; n1 < xvect.size() ; n1++)
     {
         x[n1] = xvect[n1];
     }
     /* Initialize matrixes for the differences between the locations of the bodies, two by two ;
    the distances between two bodies; and unit vectors. */
     for ( n1 = 0 ; n1 < 3 ; n1++)
     {
         vs[n1] = x[slice(9 + n1 * 3, 3, 1)] ;
         vs[n1] /= clightAd;    // Takes into account that c = G = 1 in the expressions below
         inter = pow(vs[n1], deux);
         vvs[n1] = inter.sum() ;
         for ( n2 = 0 ; n2 < 3 ; n2++)
         {
         dx[n1][n2] = x[slice(n2 * 3,3,1)];     // !! Different convention from Will's !
         dx[n1][n2] -= x[slice(n1 * 3,3,1)];

         inter = pow(dx[n1][n2], deux);
         rs[n1][n2] = sqrt ( inter.sum() );    // Distances

         ns[n1][n2] = dx[n1][n2];
             ns[n1][n2] /= rs[n1][n2];
         }
     }


     dxdt[slice(0,9,1)] = x[slice(9,9,1)]; // Set first derivatives

 // Computation of the acceleration :
 //optimization  if (Ms[2] > 0.) {   // Three-body case


     for (n1 = 0 ; n1 < n3body ; n1++)
      {
         onebodyinter =zero;

         for (n2 = 0 ; n2 < n3body ; n2++)
         {
           inter1 = zero;
           if (n2 != n1)
           {
             // 2 body terms in the direction of ns[n1][n2]
               sinter = un;
               sinter1 = dotprod3d<value_type>(vs[n1], vs[n2]);
               sinter -= (quatre + deux *  gammabar[n1][n2])*sinter1 ;
               sinter += vvs[n1] * ( un  + gammabar[n1][n2]);
               sinter += (deux + gammabar[n1][n2]) * vvs[n2];
               sinter1 = dotprod3d<value_type>(vs[n2], ns[n1][n2]);
               sinter -= troisdemis * sinter1 * sinter1;
               inter1 += ns[n1][n2] * sinter;

             // 2-body terms in the direction of (vs[n2] - vs[n1])
               inter = (vs[n2] - vs[n1]);
               inter *= - deux * gammabar[n1][n2];
               inter +=  quatre * vs[n1] - trois * vs[n2];
               sinter = dotprod3d<value_type>(ns[n1][n2], inter);//quatre * vs[n1] - trois * vs[n2] - deux * gammabar[n1][n2] * (vs[n2] - vs[n1]));
               inter1 += (vs[n2] - vs[n1]) * sinter ;
             // Prefactor 2-body
               inter1 *= Gg[n1][n2] * Ms[n2] / (rs[n1][n2]*rs[n1][n2]) ;

             // 3-body terms
             inter2 = zero;
             inter3 = zero;
             for (n3 = 0 ; n3 < n3body ; n3++)
             {
             // 3-body terms with n3 != n2
                 if (n3 != n2)
                 {
                     // terms in the direction of ns[n1][n2]
                     sinter = undemi * dotprod3d<value_type>(ns[n1][n2], ns[n2][n3]) / rs[n2][n3] ;
                     sinter -= (un + betabar[n2][n3][n1] + betabar[n2][n1][n3]) / rs[n1][n2];
                     inter = sinter * ns[n1][n2];

                     // terms in the direction of ns[n2][n3]
                     sinter = (septdemis + deux * gammabar[n1][n2]) / rs[n2][n3] ;
                     inter += sinter * ns[n2][n3];

                     // prefactor
                     inter *=  Gg[n1][n2] * Gg[n2][n3] * Ms[n2] * Ms[n3] / (rs[n1][n2]*rs[n2][n3]);

                     inter2 += inter;
                 }
              // 3-body terms with n3 != n1
                 if (n3 != n1)
                 {
                     sinter = (quatre + deux * gammabar[n1][n3] + betabar[n1][n2][n3] + betabar[n1][n3][n2]) ;
                     sinter *= Gg[n1][n2] * Gg[n1][n3] * Ms[n2] * Ms[n3] / (rs[n1][n2]*rs[n1][n2]*rs[n1][n3]);
                     inter3 -= ns[n1][n2] * sinter ;
                 }
             }
             onebodyinter += inter1 + inter2 + inter3;
           }
         }
         onebodyinter *= clightAd2;
         dxdt[std::slice(9 + n1 * 3,3,1)] = onebodyinter;
     }


 // Interaction quadrupolar momentum of inner WD with NS.
    inter =  ns[0][1];
    inter *=   (Q1Ad / pow(rs[0][1], 4)) ;
    
     dxdt[std::slice(9 + 0 * 3,3,1)] += inter*Ms[1];
     dxdt[std::slice(9 + 1 * 3,3,1)] -= inter*Ms[0];

 // Interaction quadrupolar momentum of inner WD with outer WD, shoud be negligible in J0337+1715
     inter =  ns[1][2] ;
    inter *= ( Q1Ad / pow(rs[1][2], 4) );
     dxdt[std::slice(9 + 1 * 3,3,1)] += inter*Ms[2];
     dxdt[std::slice(9 + 2 * 3,3,1)] -= inter*Ms[1];



     for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
     {
         dxdtvect[n1] = dxdt[n1];
     }

  };

void Integrateur::retro_rhs_GR_quadrupole(const state_type& xvect, state_type& dxdtvect, const value_type& t)
  {
      rhs_GR_quadrupole(xvect, dxdtvect, t);
      for ( n1=0 ; n1 < dxdtvect.size() ; n1++)
      {
          dxdtvect[n1] = - dxdtvect[n1]; // Only change from rhs here !!
      }
  }



// // RHS for 3-body integration at newtonian order
// Includes two-body newtonian integration.
void Integrateur::rhs_Newt(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    valarray<value_type> f01(3);
    valarray<value_type> f02(3);
    valarray<value_type> f12(3);
    valarray<value_type> inter(3);

    valarray<value_type> x (18);
    valarray<value_type> dxdt(18);
    value_type dist01;
    value_type dist02;
    value_type dist12;

    for (int i=0 ; i < xvect.size() ; i++)
    {
        x[i] = xvect[i];
    }


    dxdt[slice(0,9,1)] = x[slice(9,9,1)]; // Set first derivatives
    //cout << "dxdt : " << dxdt[0] << "  " << x[9] << endl;

    f01 = x[slice(0,3,1)];
    f01 -= x[slice(3,3,1)];
    f02 = x[slice(0,3,1)];
    f02 -= x[slice(6,3,1)];
    f12 = x[slice(3,3,1)];
    f12 -= x[slice(6,3,1)];

    inter = pow(f01, deux); dist01 = sqrt( inter.sum() );
    inter = pow(f02, deux); dist02 = sqrt( inter.sum() );
    inter = pow(f12, deux); dist12 = sqrt( inter.sum() );

    f01 /= pow(dist01, trois);  // force of 0 on 1 (without masses)
    f02 /= pow(dist02, trois);  //
    f12 /= pow(dist12, trois);  //

    dxdt[std::slice(9,3,1)] = un/k0 * ( - M1*f01 - M2*f02 ) *  (un + SEP_D) ;
    dxdt[std::slice(12,3,1)] = un/k0 * ( (un + SEP_D)* M0 * f01 - M2*f12 );
    if (M2 > zero)
        dxdt[std::slice(15,3,1)] = un/k0 * ( (un + SEP_D) * M0*f02 + M1*f12 ) ;
    else
        dxdt[std::slice(15,3,1)] = zero ;
    //cout << "dxdt 9: " << dxdt[10] << "  " << f01[1] << " " << k0 << " " << M1 << " f02 :" << f02[1] << endl;
    if (dxdtvect.size() < xvect.size()) dxdtvect.resize(18);
    for (int i=0 ; i < dxdtvect.size() ; i++)
    {
        dxdtvect[i] = dxdt[i];
    }

};
void Integrateur::retro_rhs_Newt(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    rhs_Newt(xvect, dxdtvect, t);
    for (int i=0 ; i < dxdtvect.size() ; i++)
    {
        dxdtvect[i] = - dxdtvect[i]; // Only change from rhs here !!
    }
}




// RHS bikeplerian : Newtonian without crossed terms !!!!!!!!! UNCOMPLETE
void Integrateur::rhs_biKeplerian(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    valarray<value_type> f01(3);
    valarray<value_type> f02(3);
    valarray<value_type> inter(3);

    valarray<value_type> x (18);
    valarray<value_type> dxdt(18);
    value_type dist01;
    value_type dist02;
    value_type dist12;


    for (int i=0 ; i < xvect.size() ; i++)
    {
        x[i] = xvect[i];
    }


    dxdt[slice(0,9,1)] = x[slice(9,9,1)]; // Set first derivatives
    //cout << "dxdt : " << dxdt[0] << "  " << x[9] << endl;

    f01 = x[slice(0,3,1)];
    f01 -= x[slice(3,3,1)];
    f02 = x[slice(0,3,1)];
    f02 -= x[slice(6,3,1)];

    inter = pow(f01, deux); dist01 = sqrt( inter.sum() );
    inter = pow(f02, deux); dist02 = sqrt( inter.sum() );

    f01 /= pow(dist01, trois);  // force of 0 on 1 (without masses)
    f02 /= pow(dist02, trois);  //

    dxdt[std::slice(9,3,1)] = un/k0 * ( - M1*f01 - M2*f02 );
    dxdt[std::slice(12,3,1)] = un/k0 * ( M0*f01 );
    if (M2 > zero)
        dxdt[std::slice(15,3,1)] = un/k0 * ( M0*f02 ) ;
    else
        dxdt[std::slice(15,3,1)] = zero ;
    //cout << "dxdt 9: " << dxdt[10] << "  " << f01[1] << " " << k0 << " " << M1 << " f02 :" << f02[1] << endl;
    if (dxdtvect.size() < xvect.size()) dxdtvect.resize(18);
    for (int i=0 ; i < dxdtvect.size() ; i++)
    {
        dxdtvect[i] = dxdt[i];
    }

};

void Integrateur::retro_rhs_biKeplerian(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    rhs_biKeplerian(xvect, dxdtvect, t);
    for (int i=0 ; i < dxdtvect.size() ; i++)
    {
        dxdtvect[i] = - dxdtvect[i]; // Only change from rhs here !!
    }
}


void Integrateur::Clear()
{
    ts.clear();
    x0.clear();
    //ts = 0.;
    //x0 = 0.;
};

// ****** Fin implémentation GR



// ******************** RHS 2-body Newtonian with quadrupole moment ***********

// // RHS 2-body integration with quadrupolar interaction at Newtonian order
// Includes two-body newtonian integration.
void Integrateur::rhs_NewtQuad(const state_type& xvect, state_type& dxdtvect, const value_type& t)
// All 3-body lines are commented
{
    valarray<value_type> f01(3);
//     valarray<value_type> f02(3);
//     valarray<value_type> f12(3);
    valarray<value_type> fq01(3);
    valarray<value_type> x (18);
    valarray<value_type> dxdt(18);
    value_type dist01;
//     value_type dist02;
//     value_type dist12;

    value_type qi=zero;
    value_type dtq = zero;

    for (int i=0 ; i < xvect.size() ; i++)
    {
        x[i] = xvect[i];
    }


    dxdt[slice(0,9,1)] = x[slice(9,9,1)]; // Set first derivatives


    f01 = x[slice(0,3,1)];
    f01 -= x[slice(3,3,1)];
//     f02 = x[slice(0,3,1)];
//     f02 -= x[slice(6,3,1)];
//     f12 = x[slice(3,3,1)];
//     f12 -= x[slice(6,3,1)];

    dist01 = sqrt(dotprod3d<value_type>(f01,f01));
//     inter = pow(f02, 2); dist02 = sqrt( inter.sum() );
//     inter = pow(f12, 2); dist12 = sqrt( inter.sum() );

    // Quadrupole force between 0 and 1
    fq01 = f01;
    fq01 *=  (neufdemis * quadrupole[0] / (pow(dist01, 5) *M1));

    // Quadrupole time dependence
    dtq = (quadrupole[3] - quadrupole[2])/trois;
    if (t > quadrupole[2] && t <= quadrupole[2] + dtq) fq01 *= un - quadrupole[1] * ( t - quadrupole[2])/( quadrupole[3] - quadrupole[2])  ;
    qi = quadrupole[1] * ( dtq )/( quadrupole[3] - quadrupole[2])  ;
    if (t > quadrupole[2] + dtq  && t <= quadrupole[3] ) fq01 *= un - qi + quadrupole[1] * ( t - quadrupole[2] - dtq )/( quadrupole[3] - quadrupole[2])  ;
    qi = quadrupole[1] * ( quadrupole[3] - quadrupole[2] - dtq )/( quadrupole[3] - quadrupole[2]) - qi;
    if (t > quadrupole[3] ) fq01 *= un + qi;
//     else
//         printf(" not it quad3\n");

    // Normal Newtonian forces
    f01 /= pow(dist01, 3);  // force of 0 on 1 (without masses)
//     f02 /= pow(dist02, 3);  //
//     f12 /= pow(dist12, 3);  //


    dxdt[std::slice(9,3,1)] = GgravAd * ( - M1*f01    - fq01 * M1)  ;
    dxdt[std::slice(12,3,1)] = GgravAd * (  M0 * f01   + fq01*M0) ;
    dxdt[std::slice(15,3,1)] = zero ;
//     dxdt[std::slice(9,3,1)] = un/k0 * ( - M1*f01 - M2*f02 ) *  (un + SEP_D) - fq01*M1  ;
//     dxdt[std::slice(12,3,1)] = un/k0 * ( (un + SEP_D)* M0 * f01 - M2*f12 ) + fq01*M0 ;
//     if (M2 > zero)
//         dxdt[std::slice(15,3,1)] = un/k0 * ( (un + SEP_D) * M0*f02 + M1*f12 ) ;
//     else
//         dxdt[std::slice(15,3,1)] = zero ;
    if (dxdtvect.size() < xvect.size()) dxdtvect.resize(18);
    for (int i=0 ; i < dxdtvect.size() ; i++)
    {
        dxdtvect[i] = dxdt[i];
    }

};

void Integrateur::retro_rhs_NewtQuad(const state_type& xvect, state_type& dxdtvect, const value_type& t)
{
    rhs_NewtQuad(xvect, dxdtvect, t);
    for (int i=0 ; i < dxdtvect.size() ; i++)
    {
        dxdtvect[i] = - dxdtvect[i]; // Only change from rhs here !!
    }
}

// ******************** End RHS 2-body Newtonian with quadrupole moment ***********
