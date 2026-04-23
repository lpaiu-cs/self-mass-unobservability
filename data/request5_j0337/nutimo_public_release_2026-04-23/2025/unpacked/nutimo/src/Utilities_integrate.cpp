// SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

/* This code integrates the equations of motion for the triple system J0337+1715 (Ransom et al. 2013). The inputs are Nutimo parameters
 *
 * Written by Guillaume Voisin 2021 , LUTh, Observatoire de Paris, PSL Research University, CNRS (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 *
 */

#ifdef MPIMODE
    #include "mpi.h"
#endif
#include "Utilities_integrate.h"



int Read_resume_file(char * resumefile, valarray<value_type> & rp, valarray<value_type> &rpt, 
                      valarray<value_type> &ri, valarray<value_type> &rit, 
                      valarray<value_type>& ro, valarray<value_type> &rot, 
                      valarray<valarray<value_type>>& r_extra, 
                      valarray<valarray<value_type>>& v_extra)
// resumefile : filename of a file from which the first line is read 
// rp, rpt : position and velocity vector of the pulsar 
// ri, rit : position and velocity vector of the inner WD
// ro, rot : position and velocity vector of the outer WD
// r_extra[0], v_extra[0] : position and vleocity vector of the first extra body. Must be pre-allocated
// The expected format is rp[0] rp[1] rp[2] rpt[0] ... rot[2] r_extra[0][0] ... v_extra[nextra-1][2]
{
    FILE * fp; 
    int ncomp = (6* 3 + 6*r_extra.size());
    int sizechar = ncomp*30;
    char line[sizechar];
    char * pch;
    int i, k, rval ; 
    value_type components[ncomp];
    
    for (i=0; i< ncomp; i++) components[i]=0.;
    
    fp = fopen(resumefile,"r");
    if (fp == NULL)
    {
        printf("\n\n Opening of resumefile %s failed ! \n\n", resumefile);
        return 0;
    }
    if (fgets(line, sizechar, fp) != NULL)
    {
        i=0;
        pch = strtok(line, " \t");
        while (i < ncomp and pch !=NULL)
        {
            sscanf(pch, "%Le", &components[i]);
//             printf("%s  %.15Le\n", pch, components[i]);
            pch = strtok(NULL, " \t\n");
            i++;
        }
        if (i < ncomp) printf("\n Warning : only %d components read out of %d in Read_resume_file !\n\n", i, ncomp);
    }
    fclose(fp);
    rval = i;
    
    for (i=0; i<3; i++) 
    {
        rp[i] = components[i];
        rpt[i] = components[i+3];
        ri[i] = components[i+6];
        rit[i] = components[i+9];
        ro[i] = components[i+12];
        rot[i] = components[i+15];
        for(k=0; k < r_extra.size(); k++)
        {
            r_extra[k][i] = components[i+18+k*6];
            v_extra[k][i] = components[i+21+k*6];
        }
    }
    return rval;
}


void Compute_integrals_of_motion(value_type *** states, const int n, 
                                 const valarray<value_type>  Ms, 
                                 const valarray<valarray<value_type>>  int_Gg, 
                                 const valarray<valarray<value_type>>  int_gammabar, 
                                 const valarray<valarray<valarray<value_type>>> int_betabar, 
                                 const int integrator_type, 
                                 valarray<valarray<value_type>> & center_of_mass_positions,
                                 valarray<valarray<value_type>> & center_of_mass_impulsions,
                                 valarray<value_type>  & energies, valarray<valarray<value_type>> & angular_momentum)
// !! Angular momentum is only computed at Newtonian level 
{
    value_type nrj;
    valarray<value_type> x(0.,3);
    valarray<value_type> v(0.,3);
    valarray<value_type> inter(0.,3);
    value_type mbin, mtri, eta2, eta3, r, vv;
    
    mbin = Ms[0] + Ms[1];
    mtri = mbin + Ms[2];
    eta2 = Ms[0] *Ms[1]/pow(mbin,2);
    eta3 = mbin *Ms[2]/pow(mtri,2);
        
    center_of_mass_positions.resize(n);
    center_of_mass_impulsions.resize(n);
    energies.resize(n);
    angular_momentum.resize(n);
    
    if (integrator_type == 3 or integrator_type==1)
    {
        valarray<value_type> sv(6* Ms.size());
        int ncomp = 3*Ms.size(); 
        for (long int i = 0 ; i < n ; ++i)
        {
            angular_momentum[i] = valarray<value_type>(0.,3);
            for (int j = 0; j < Ms.size(); j++) 
            {
                for (int k = 0; k < 3; k++) sv[3*j+k] = states[j][i][k];
                for (int k = 0; k < 3; k++) sv[3*j+ ncomp + k] = states[j][i][k+3];
            }
            // Computing angular momentum to Newtonian level
            for (int j = 0; j < Ms.size(); j++)
            {
                x = sv[slice(3*j,3,1)];
                v = sv[slice(3*j+ncomp,3,1)];
                angular_momentum[i] += Ms[j] * crossprod(x,v);               
            }
            // 1PN contribution of the inner binary to the angular momentum:
            x = sv[slice(0,3,1)];
            x -= sv[slice(3,3,1)]; // x1 -x0
            v = sv[slice(ncomp, 3,1)];
            v -= sv[slice(ncomp+ 3, 3,1)]; // v1 - v0
            vv = sumsquares3d(v);
            r = norm3d(x);
            angular_momentum[i] += mbin*eta2*(0.5*(1-3*eta2)*vv + (3+eta2) *GMsol*mbin/r)/clight2 *crossprod(x,v);
            
            // 1PN contribution of the binary between the inner binary and the outer WD to the angular momentum:
            mtri = mbin + Ms[2];
            x = sv[slice(0,3,1)];
            x *= Ms[0]/mbin;
            inter = sv[slice(3,3,1)];
            inter *= Ms[1]/mbin;
            x += inter; // xbin
            x -= sv[slice(6,3,1)]; // xbin -x2
            v = sv[slice(ncomp,3,1)];
            v *= Ms[0]/mbin;
            inter = sv[slice(ncomp+3,3,1)];
            inter *= Ms[1]/mbin;
            v += inter; // vbin
            v -= sv[slice(ncomp+6,3,1)]; // vbin -v2
            vv = sumsquares3d(v);
            r = norm3d(x);
            angular_momentum[i] += mtri*eta3*(0.5*(1-3*eta3)*vv + (3+eta3) *GMsol*mtri/r)/clight2 *crossprod(x,v);
            

            // Compute COF and Energy
            IntegralePrems3_1PN( sv, Ms, 0.,
                                    int_Gg, int_gammabar, int_betabar,
                                    center_of_mass_positions[i], center_of_mass_impulsions[i],
                                    nrj  ) ;
            
            energies[i] = nrj;
        }
    }
    else if (integrator_type == 30)
    {
        valarray<value_type> sv(6* Ms.size());
        int ncomp = 3*Ms.size(); 
        for (long int i = 0 ; i < n ; ++i)
        {
            for (int j = 0; j < Ms.size(); j++) 
            {
                for (int k = 0; k < 3; k++) sv[3*j+k] = states[j][i][k];
                for (int k = 0; k < 3; k++) sv[3*j+ ncomp + k] = states[j][i][k+3];
            }
            
            // Compute COF and Energy and ang momentum
            IntegralePrems_0PN( sv, Ms, 0.,
                                    int_Gg, 
                                    center_of_mass_positions[i], center_of_mass_impulsions[i],
                                    nrj, angular_momentum[i]  ) ;
            
            energies[i] = nrj;
        }
    }
    else
    {
        printf("\n integratortype not implemented in Compute_integrals_of_motion\n\n");
    }
//     else if (integrator_type == 2)
//     {
//         cout << endl << "Computing 1PN integrals of motion with inner WD quadrupole" << endl ;
//         for (long int i = 0 ; i < ninterp ; ++i)
//             IntegralePrems3_1PN( sp[i], si[i], so[i],
//                                     int_M0, int_M1, int_M2, int_quadrupole[0],
//                                      int_Gg, int_gammabar, int_betabar,
//                                     center_of_mass_positions[i], center_of_mass_impulsions[i],
//                                     energies[i]  ) ;
//     }
// 
//     else if (integrator_type == 0)
//     {
//         cout << endl << "Computing Newtonian integrals of motion" << endl ;
//         for (long int i = 0 ; i < ninterp ; ++i)
//             IntegralePrems3_Newt( sp[i], si[i], so[i],
//                                     int_M0, int_M1, int_M2,
//                                     center_of_mass_positions[i], center_of_mass_impulsions[i],
//                                     energies[i]  ) ;
//     }
//     else if (integrator_type == 10)
//     {
//         cout << endl << "Computing 2-body Newtonian with quadrupole term integrals of motion" << endl ;
//         for (long int i = 0 ; i < ninterp ; ++i)
//         {
//             IntegralePrems2_NewtQuad( sp[i], si[i],
//                                     int_M0, int_M1,  int_quadrupole[0],
//                                     center_of_mass_positions[i], center_of_mass_impulsions[i],
//                                     energies[i]  ) ;
//         }
//     }
    return ;
}



void Save_meanlongs(char * fname, valarray<value_type> ts_day, valarray<value_type>  l)
// Used to save to file mean longitudes 
{
     FILE * fout = fopen(fname, "w");
     int i=0;
     for (i=0 ; i < l.size() - 1; i ++) fprintf(fout, "%.19Le    %.19Le\n", ts_day[i], l[i]);
     fprintf(fout, "%.19Le    %.19Le", ts_day[i], l[i]);
     fclose(fout);
}


void Save_orbels(char * fname, valarray<value_type> & ts_day, orbel_t ** orbels, int nbody)
// Used to save to file nbody orbital elements
{
     FILE * fout = fopen(fname, "w");
     int i=0;
     int j;
     fprintf(fout, "#t   (a Porb ecc i Oman omperi tperi)_{1..nbody-1}");
     for (i=0 ; i < ts_day.size() ; i ++) 
     {
         fprintf(fout, "\n%.19Le    ", ts_day[i]);
         for (j = 0; j < nbody -1 ; j ++) 
         {
             fprintf(fout, "    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    %.19Le    ", orbels[i][j].a, orbels[i][j].Porb, orbels[i][j].ecc, orbels[i][j].i, orbels[i][j].Oman, orbels[i][j].omperi, orbels[i][j].tperi);
         }
     }
     fclose(fout);
}


void Count_lines_in_file(char * fname, int & nlines, int & ncomments, int & nempty, int max_line_size)
// Total number of lines in file "fname" = nlines + ncomments + nempty
// Comment character is "#"
// Separating character is " " or "\t" or "\n"
{
    int l;
    FILE * fp = fopen(fname, "r");
    char line[max_line_size];
    char * pch;
    char com='#';
    
    nlines=0;
    nempty=0;
    ncomments=0;
    
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
            nlines++;
        l++;
    }
    fclose(fp);
};

void Read_trajectory(char * fname, int nbody,  valarray<long double> & ts_day, valarray<valarray<valarray<long double>>> &  svs)
// Used to read trajectory files generated by "integrate.exe" 
{
    int i,j,l;
    int ncomp = nbody*6 + 1;
    int sizechar = ncomp*100;
    char line[sizechar];
    char * pch;
    value_type components[ncomp];
    char com='#';

// Initialise arrays
    int nlines, ncomments, nempty;
    Count_lines_in_file(fname, nlines, ncomments, nempty, sizechar);
    ts_day.resize(nlines);
    svs.resize(nlines);
    for (l=0; l < nlines; l++) 
    {
        svs[l].resize(nbody);
        for (i=0; i < nbody; i++) svs[l][i].resize(6);
    };
    
    if (nlines ==0 ) 
    {
        printf("\nFile %s empty !\n\n", fname);
        return;
    }
    
// Open file 
    FILE * fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("\n\n Opening of trajectory file %s failed ! \n\n", fname);
        return ;
    };

// Read file 
    l=0;
    while (fgets(line, sizechar, fp) != NULL)
    {
        i=0;
        pch = strtok(line, " \t\n");
        
        while (i < ncomp and pch !=NULL and pch[0] != com)
        {
            sscanf(pch, "%Le", &components[i]);
            pch = strtok(NULL, " \t\n");
            i++;
        }
        if (i < ncomp)
        {
            if (pch[0] != com)
                printf("Line %d : skipping due to comment sign in '%s' \n", l, fname);
            else
                printf("\n Line %d : Warning : only %d components read out of %d in '%s' !\n\n", l, i, ncomp, fname);
        }
        else
        {
            ts_day[l] = components[0];
            for (j = 1; j < ncomp; j++)
            {
                svs[l][(j-1)/6][(j-1)%6] = components[j];
            };
            l++;
        }
    }
    fclose(fp);
};


void Save_traj(char * filename, valarray<value_type> ts_day, valarray<valarray<valarray<value_type>>> trajs) 
{
    int i, j, k, nsys, nlines;
    nsys= trajs[0].size();
    nlines = trajs.size();
    FILE * myfile ;
    myfile = fopen(filename, "w") ;
    for(i=0 ; i < nlines ; ++i) 
    {
        fprintf(myfile, "%.19Le    ", ts_day[i] );
        for (j = 0; j < nsys ; j++ ) 
        {
             for (k=0; k < 5 ; k ++) fprintf(myfile, "%.19Le    ", trajs[i][j][k] );
             if (j == nsys -1)
                 fprintf(myfile, "%.19Le\n", trajs[i][j][5] );
             else 
                 fprintf(myfile, "%.19Le        ", trajs[i][j][5] );
        }
     }
     fclose(myfile);
}

void Read_orbels(char * fname, int nbody,  valarray<long double> & ts_day, orbel_t ** orbels, int nlines  )
// Used to read trajectory-orbels files generated by "integrate_convert.exe" 
{
    int i,j,k;
    int ncols = 7*(nbody-1)+1;
    int max_line_size = ncols *100;
    valarray<valarray<value_type>> table;
    Loadtxt(fname, table, nlines, ncols, max_line_size);
    
    ts_day.resize(nlines); 
    for (i = 0; i < nlines; i++)
    {
        ts_day[i] = table[i][0];
        for (j = 0 ; j < nbody -1; j++)
        { //(a Porb ecc i Oman omperi tperi)_{1..nbody-1}
            orbels[i][j].a = table[i][1 + j*7 + 0];  
            orbels[i][j].Porb = table[i][1 + j*7 + 1];
            orbels[i][j].ecc = table[i][1 + j*7 + 2];
            orbels[i][j].i = table[i][1 + j*7 + 3];
            orbels[i][j].Oman = table[i][1 + j*7 + 4];
            orbels[i][j].omperi = table[i][1 + j*7 + 5];
            orbels[i][j].tperi = table[i][1 + j*7 + 6];
            orbels[i][j].norb = deuxpi/orbels[i][j].Porb;
            orbels[i][j].m = (ts_day[i]* daysec - orbels[i][j].tperi) * orbels[i][j].norb;
            orbels[i][j].v = 0.;
            orbels[i][j].E = 0.;
            orbels[i][j].tasc_approx = 0.;
            for (k = 0; k<3; k++)
            {
                orbels[i][j].h[k] = .0;
                orbels[i][j].qlap[k] = .0;
                orbels[i][j].easc[k] = .0;
            }
        }
    }
};
