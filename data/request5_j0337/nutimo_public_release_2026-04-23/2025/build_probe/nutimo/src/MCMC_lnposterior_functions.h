/*
 * SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 *
 * Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 */


#ifndef MCMC_lnposterior_functions_h
# define MCMC_lnposterior_functions_h

#ifdef MPIMODE
    #include "mpi.h"
#endif
#include <iostream>
#include <chrono>
#include <random>
#include <cfloat>
#include "Utilities.h"
#include <stdio.h>
#include <string.h>

#include "Fittriple.h"


class fctgaussienne
{
    int ndim ;
public:
    int verbose=0;

    fctgaussienne(int number_of_dim){ndim =number_of_dim;};
    double operator()(double* pars)
    {
        double res = 0.;
        double inter = 0.;
//          res -= pow(pars[0] - 1000.*pars[1],2);
//          res -= pow(pars[0] + pars[1],2);
        for (int i = 0; i < ndim ; i++)
        {
            inter = pow(pars[i] - 0.2 * static_cast<double>(i),2);
            res -= inter / static_cast<double>(i+1) ;
        }
        return res;
    };
    void Print_MCMC()
    {
        printf("Fonction Gaussienne in %d dimensions! \n", ndim);
    };
};



class fctsimudata
{
    double std ;
    int ndim=1;
    int ndata;
    double efac = 1.;
    double * ts ;
    double ** posterior_cormat;
public:
    double * data;
    double ** coeffs;
    int verbose=0;
    int errorscale=0; // activate the hyperparameter error scale, as the las parameter of pars

    fctsimudata()
    {
    };

    void external_data(double * simudata, double ** coefficients, int numberofdata, int ndimensions, double stddev)
    {
      int i,j;
      ndata = numberofdata;
      data = (double*) malloc(sizeof(double) * ndata);
      std = stddev;
      ndim = ndimensions;
      coeffs = (double **)malloc(sizeof(double*) * ndata);
      for (i =0 ; i < ndata ; i ++) coeffs [i] = (double*) malloc(sizeof(double) * ndim);
      for (i =0 ; i < ndata ; i ++)
      {
        data[i] = simudata[i];
        for (j =0 ; j < ndim ; j ++) coeffs[i][j] =  coefficients[i][j];
      }
    }

    void generate_data( int numberofdata,  double * parameters, double stddev, int ndimensions, double basefreq, double timespan)
    {
        int i,j;
        double mean ;
        ndata = numberofdata;
        data = (double*) malloc(sizeof(double) * ndata);
        std = stddev;
        ndim = ndimensions;
        ts = (double *) malloc(sizeof(double) * ndata);
        coeffs = (double **)malloc(sizeof(double*) * ndata);
        for (i =0 ; i < ndata ; i ++) coeffs [i] = (double*) malloc(sizeof(double) * ndim);
        for (i =0 ; i < ndata ; i ++)
        {
          ts[i] = timespan / ndata * i;
          for (j =0 ; j < ndim ; j ++) coeffs[i][j] =  cos(j*basefreq * ts[i]);
        }
            // Random generator init
        // obtain a seed from the system clock:
        unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
        mt19937 randomgen(seed1);  // mt19937 is a standard mersenne_twister_engine
        normal_distribution<double> normaldist(0.,std); // used for initialization

        for ( i = 0; i < ndata ; i++)
        {
          mean = 0.;
          for (j = 0; j < ndim  ; j ++ ) mean += parameters[j] *coeffs[i][j];
          data[i] = mean + normaldist(randomgen);
        }
        double var = 0;
        for (i = 0; i < ndata ; i++) var += (data[i] - mean)*(data[i] - mean);
        var /= ndata;
        var = sqrt(var);
        printf("--------var %f %f %f\n", var, var/ std, std );
      //  std /=1.2;
        printf("std = %f\n\n", std);

// Compute correlation matrix of the posterior distribution
        int k,l;
        posterior_cormat = (double **)malloc(sizeof(double*) * ndata);
        for (i =0 ; i < ndim ; i ++) posterior_cormat[i] = (double*) malloc(sizeof(double) * ndim);
        for (k =0 ; k < ndim ; k ++)
        {
          for (l =0 ; l < ndim ; l ++)
          {
            posterior_cormat[k][l] = 0;
            for (i = 0; i < ndata; i ++) posterior_cormat[k][l] += coeffs[i][k] * coeffs[i][l] ;
            posterior_cormat[k][l] /=   std * std;
          };
        };

// Save Intial conditions and posterior correlation
        int mpirank = 0;
        char  filename[50];
        sprintf(filename,"simudata.txt");
#ifdef MPIMODE
        MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
#endif
        if (mpirank == 0)
        {
          sprintf(filename,"simu-data.txt");
          Savetxt(filename, data, ndata);
          sprintf(filename,"simu-times.txt");
          Savetxt(filename, ts, ndata);
          sprintf(filename,"simu-model_coeffs.txt");
          Savetxt(filename, coeffs, ndata, ndim);
          sprintf(filename,"simu-posterior_correlation_matrix.txt");
          Savetxt(filename, posterior_cormat, ndim, ndim);
        }
        free(ts);
        for (i =0 ; i < ndim ; i ++) free(posterior_cormat [i]);
        free(posterior_cormat);
    };

    ~fctsimudata()
    {
      int i = 0;
      free(data);
      for (i =0 ; i < ndata ; i ++) free(coeffs [i]);
      free(coeffs);
    };

    double operator()(double* pars, double beta)
    {
        double parfinal[ndim+1];
        memcpy(parfinal, pars, sizeof(double) * ndim);
        if (errorscale >0)
        {
            if (beta != 1.)
            {
                parfinal[ndim]=0.;
                pars[ndim] =efac -1.;
            }
            else
            {
                parfinal[ndim] = (pars[ndim]);
            }
        }

        return lnposterior(parfinal) * beta;
    }


    double operator()(double* pars)
    {
        return lnposterior(pars);
    }

    double lnposterior(double * pars)
    {

        double res = 0.;
        double inter = 0.;
        double currentefac = 1.;
        double mean ;
        int j ;

        if (errorscale >0) currentefac = 1. + pars[ndim];

        res = 0.;
        for (int i = 0; i < ndata ; i++)
        {
            mean = 0.;
            for (j = 0; j < ndim  ; j ++ ) mean += pars[j] *coeffs[i][j];
            res -= 0.5*(data[i] - mean)*(data[i] - mean)/( std*std * currentefac *currentefac);
        }
        res -= ndata * log(currentefac);

//         printf("std efac %f %f \n", std, efac);

        return res;
    };


    double Swap_temperatures(double * params_lowtemp, double & lnposterior_lowtemp, const double beta_lowtemp, double * params_hightemp, double & lnposterior_hightemp, const double beta_hightemp )
    {

        int locndim = ndim;
        double walkerbuffer[locndim];
        double lnposteriory = lnposterior_lowtemp;
        double efac_temp1 = 0.;
        double efac_lowtemp = 1.;
        double efac_hightemp = 1.;
        double logefac =0.;


        int iefac =-1;
        if (errorscale >0)
        {
            iefac = 1;
            locndim +=1;
        }

        if (iefac >= 0)
        {
            efac_temp1 = 1 + params_lowtemp[iefac];
            efac_lowtemp = 1 + params_lowtemp[iefac];
            efac_hightemp = 1 + params_hightemp[iefac];
            logefac = ndata* log(efac_lowtemp / efac_hightemp);
        }

        lnposterior_lowtemp = lnposterior_hightemp * beta_lowtemp / beta_hightemp;
        lnposterior_hightemp = lnposteriory * beta_hightemp / beta_lowtemp;

        memcpy(walkerbuffer, params_lowtemp, locndim*sizeof(double));
        memcpy(params_lowtemp, params_hightemp, locndim*sizeof(double));
        memcpy(params_hightemp, walkerbuffer, locndim*sizeof(double));

        if (beta_lowtemp == 1. and iefac >= 0)
        {
            lnposterior_lowtemp += logefac;
            lnposterior_hightemp -= beta_hightemp * logefac;
            params_hightemp[iefac] = efac- 1.;
            params_lowtemp[iefac] = efac_temp1 -1.;
        }
    };


    double Get_chi2(const double lnposterior, const double * params, const double beta)
    // Return the chi2 based on the log of the posterior probability density, therefore removing any efac contribution.
    // Does not recompute the lnposterior : fast
    // Remove effect of temperature. ( i.e. gives chi2 with temperature = 1)
    {
        int iefac = 1;
        double currentefac =1.;
        if (errorscale >0 ) currentefac = 1. + params[iefac];

        return -(lnposterior/beta - ndata*log(currentefac));
    };



    void Print_MCMC()
    {
        printf("Fonction Simu Data in %d dimensions! \n", ndim);

    };

    void Get_parameter_name(int paramnumber, char * name)
    {
        sprintf(name, "p%i", paramnumber);
    }

    double initial_distribution(mt19937& random_generator, int param_number)
    {
        normal_distribution<double> initial_distribution(0.,0.0001);

        return initial_distribution(random_generator);
    }


    int Get_ndim() // return number of dimensions
    {
        return ndim;
    };

};


class fctmultigaussienne
{
    int ndim ;
    int npeaksperdim;
    double stdmainpeak=1.;
    double stdpeaksec;
    double peakstep ;
public:
    int verbose=0;
    int errorscale=0; // activate the hyperparameter errror scale, as the las parameter of pars

    fctmultigaussienne(int number_of_dim, int half_number_of_peaks_per_dim)
    {
        ndim =number_of_dim;
        npeaksperdim=half_number_of_peaks_per_dim;
        stdpeaksec = stdmainpeak / (4.* npeaksperdim);
        peakstep = stdmainpeak;
    };

    double operator()(double* pars, double beta)
    {
        double parfinal[ndim+1];
        memcpy(parfinal, pars, sizeof(double) * ndim);
        if (errorscale >0)
        {
            if (beta != 1.)
                parfinal[ndim]=1.;
            else
                parfinal[ndim] = (1.+ pars[ndim]);
        }

        return lnposterior(parfinal) * beta;
    }


    double operator()(double* pars)
    {
        return lnposterior(pars);
    }

    double lnposterior(double * pars)
    {
        double res = 0.;
        double inter = 0.;
        int k ;

//         printf("pars %.5e %.5e \n", pars[0], pars[1]);

        for (int i = 0; i < ndim ; i++)
        {
            for (int j = -npeaksperdim ; j <= npeaksperdim ; j++)
            {
                if (j != 0)
                {
                    inter = 0.;
                    for (k = 0; k < ndim ; k++)
                    {
                        if (k == i )
                            inter += -pow((pars[i] - peakstep*j)/stdpeaksec,2)*0.5;
                        else
                            inter += -pow((pars[k])/stdpeaksec,2) * 0.5;
//                         printf("i %d j %d k %d inter %.5e \n", i, j, k, inter);
                    }
                    res += exp(inter);
//                     printf("i %d j  %d exp(inter) %.5e \n", i, j, exp(inter));
                }
            }
            inter = 0.;
            for (k = 0; k < ndim ; k++) inter += -pow((pars[k])/stdmainpeak,2) *0.5;
            //res += exp(inter);
        }
        res = log(res);
        if (errorscale >0)
        {
                res /= pars[ndim] * pars[ndim];
                res += log(pars[ndim]);
        }
        return res;
    };

    double Get_chi2(const double lnposterior, const double * params, const double beta)
    {
      return lnposterior;
    };

    void Print_MCMC()
    {
        printf("Fonction Multi Gaussienne in %d dimensions! \n", ndim);
        printf("npeaksperdim %d stdmainpeak %.5e stdpeaksec %.5e peakstep %.5e \n", npeaksperdim, stdmainpeak, stdpeaksec, peakstep);
    };

    double initial_distribution(mt19937& random_generator, int param_number)
    {
        normal_distribution<double> initial_distribution(0.,0.0001);

        return initial_distribution(random_generator);
    }

    void Get_parameter_name(int paramnumber, char * name)
    {

        sprintf(name, "%i", paramnumber);
    }

    int Get_ndim() // return number of dimensions
    {
        return ndim;
    };



    double Swap_temperatures(double * params_lowtemp, double & lnposterior_lowtemp, const double beta_lowtemp, double * params_hightemp, double & lnposterior_hightemp, const double beta_hightemp )
    {
      printf("Swap_temperatures NOT IMPLEMENTED !");
    };


};


class fctFittriple : public Fittriple // Wraps Fittriple object with a prior and a different parameter set
{
    const value_type secyr = 31557600.0; // Number of seconds per year

    // Position from Ransom et al. 2014
    const value_type posransom_rad = 9.5002822899862595762e-01;
    const value_type posransom_dec = 3.0114118619832732786e-01;
    const value_type posprior = 2;

    // Radial velocity from Kaplan et al. 2014
    const value_type rvkaplan = 29.7e+03; //m/s
    const value_type rverrorkaplan = 0.9e+03; // m/s
    const value_type rvprior = 2;

    // Photometric distance from Kaplan et al. 2014
    const value_type dkaplan = 4078.5496710122447 ; // lyr
    const value_type derrorkaplan = 250.98767206229201; // lyr
    const value_type dprior = 2;

public :

    void Print_MCMC()
    {
        printf("Position prior is +/- %.2Le mas.\n", posprior);
        printf("Radial velocity prior is %.2Le +/- %.2Le km/s.\n", rvkaplan/1000., rverrorkaplan*rvprior/1000.);
        printf("Distance prior is %.3Le +/- %.3Le ly.\n", dkaplan, dprior*derrorkaplan);
        double chi2 = Compute_lnposterior(0);
        printf("Initial reduced chi2 is %.8f, and chi2 = %.8f\n", chi2, chi2*ntoas );
        printf("\n Parameter map : [");
        for (int i=0; i < Get_nfitted_params() ; i++) printf("%i, ", Get_list_of_fitted_parameters(i));
        printf("]\n");
    };

    double operator()(double * relativeshift)
    {
        vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
        Set_fitted_parameter_relativeshifts(vectrelativeshift);

        if ( (parameters.RA < posransom_rad - posprior *radmasdeg )   or  (parameters.RA > posransom_rad + posprior*radmasdeg) )
        {
          //  printf("prior ra %.15Le %.15Le %.15Le %.15Le\n", parameters.RA, posransom_rad, posransom_rad - posprior *radmasdeg, posransom_rad + posprior*radmasdeg);
            return  -DBL_MAX ;
        }
        else if ( (parameters.DEC < posransom_dec - posprior* radmasdeg) or  (parameters.DEC > posransom_dec + posprior*radmasdeg) )
            {
            //printf("prior rdec\n");
            return  -DBL_MAX ;
            }
        else if ( (parameters.distance < dkaplan - dprior* derrorkaplan) or  (parameters.distance > dkaplan + dprior* derrorkaplan) )
            {
  //          printf("prior d\n");
            return  -DBL_MAX ;
            }
        else if ( (parameters.distance1*radmasdeg*parameters.distance*clight < rvkaplan - rvprior* rverrorkaplan) or  (parameters.distance1*parameters.distance*clight > rvkaplan + rvprior* rverrorkaplan) )
        {
//            printf("prior rv %.15Le %.15Le %.15Le %.15Le\n", parameters.distance1*parameters.distance*clight , rvkaplan,  rvkaplan - rvprior* rverrorkaplan,rvkaplan + rvprior* rverrorkaplan);
            return  -DBL_MAX ;
        }
        else
            return Compute_lnposterior(0) * ntoas; // Return the actual log(proba) = -chi2 , not the reduced one

    };



};





class fctFittriple_gaussprior : public Fittriple // Wraps Fittriple object with a Gaussian prior
{
    const value_type secyr = 31557600.0; // Number of seconds per year

//     // Position from Ransom et al. 2014 (radians)
//     const value_type posransom_rad = 9.5002822899862595762e-01;
//     const value_type posransom_dec = 3.0114118619832732786e-01;
//     const value_type raerror = 9.45386657846825e-09;
//     const value_type decerror = 9.69627362219072e-09;
//     const value_type posprior = 2;
//
//     // Radial velocity from Kaplan et al. 2014
//     const value_type rvkaplan = 29.7e+03; //m/s
//     const value_type rverrorkaplan = 0.9e+03; // m/s
//     const value_type rvprior = 2;
//
//     // Photometric distance from Kaplan et al. 2014
//     const value_type dkaplan = 4078.5496710122447 ; // lyr
//     const value_type derrorkaplan = 250.98767206229201; // lyr
//     const value_type dprior = 2;

//     // Distance from Gaia DR2. Warning : POSEPOCH = MJD2015.5
//     const value_type dmean = 4604.639595490277 ; // lyr
//     const value_type derror = 1558.9119367898795; // lyr
//     const value_type dprior = 1;

    // Position from GAIA DR2 (radians) Warning : POSEPOCH = MJD2015.5
    const value_type ramean = 9.500283096783616e-01;
    const value_type decmean = 3.011411347204199e-01;
    const value_type raerror = 8.945972920351486e-10; // =  0.184523937110807 mas
    const value_type decerror = 9.302791001352137e-10; // = 0.191883838345113 mas
    const value_type posprior = 2;

    // Tangential proper motion (mas/yr) from Gaia DR2. Warning : POSEPOCH = MJD2015.5
    const value_type pmramean = 4.81377140723352;
    const value_type pmdecmean = -4.42182046193475;
    const value_type pmraerror = 0.498024673422251;
    const value_type pmdecerror = 0.42770795176738;
    const value_type tpmprior = 2;

    // Radial velocity from Kaplan et al. 2014
    const value_type rvkaplan = 29.7e+03; //m/s
    const value_type rverrorkaplan = 0.9e+03; // m/s
    const value_type rvprior = 1;

    // Prior centered on Kaplan 2014 with ~2error kaplan.
    const value_type dmean = 4350.;//4604.639595490277 ; // lyr
    const value_type derror = 250.; //1558.9119367898795; // lyr
    const value_type dprior = 2;

    // Normalisation constants for the (acos(i), asin(i)) priors
    const value_type nctei = 911596921.5018525; // = ai / sinii at mean for end of DataMarch2018/Run20180509
    const value_type ncteo = 55954106240.898094; // = ao / sinio at mean for end of DataMarch2018/Run20180509

    // Constant EFAC for temperatures > 1
    const value_type ref_efac = 1.57878505e+00 - 0.25; //for quad on coma 04/06/18 // mean for end of DataMarch2018/Run20180509

public :
    int verbose=0;
    double beta =1.; // 1/Temperature

    const int name_size=100; // size of parameter names

    double lastlnlikelyhood;
    double lastprior;
    
    int counter =0; // counts number of calls to lnposterior

    int Get_ndim() // return number of dimensions
    {
        return Get_nfitted_params();
    };


    void Print_MCMC()
    {
        double chi2 =zero;
        double params[Get_ndim()] ;

        for (int i = 0; i < Get_ndim() ; i++) params[i] = 0.;

        double totlnposterior = lnposterior(params);

        printf("Beta = 1/temperature = %.5e \n",beta);
        printf("Gaussian priors at 1sigma :\n");
        printf(" Position prior is RA +/- %.2Le mas, DEC +/- %.2Le mas.\n", posprior*raerror *masdegrad, posprior*decerror *masdegrad);
        printf(" Transverse proper motion prior is pmra +/- %.2Le mas/yr, pmdec +/- %.2Le mas/yr.\n", posprior*pmraerror, posprior*pmdecerror);
        printf(" Radial velocity prior is %.2Le +/- %.2Le km/s.\n", rvkaplan/1000., rverrorkaplan*rvprior/1000.);
        printf(" Distance prior is %.3Le +/- %.3Le ly.\n", dmean, dprior*derror);
        printf(" (acos(i), asin(i)) prior is sin(i)/a (from  Prior(i) propto sin(i) and  Jacobian[(a,i)->(acos(i), asin(i))] )\n" );
        printf("\nOther priors :\n");
        printf(" Quadrupole prior : None\n");
        printf("\n");
        printf("Initial lnposterior/ntoas = %.8f, lnposterior = %.8f\n", totlnposterior/ntoas, totlnposterior );
        printf("Initial (lnposterior - prior)/ntoas = %.8f, (lnposterior-prior) = %.8f\n", lastlnlikelyhood/ntoas, lastlnlikelyhood );
        printf("Initial prior = %.8f\n", lastprior);
        printf("EFAC fitted parameter index is %i, absolute index %i . Reference EFAC value = %.5Lf \n", parameters.Get_fitted_parameter_index(36),parameters.absolute_parameter_map[36], ref_efac);
        chi2 = Get_chi2(lastlnlikelyhood, params, 1.);
        printf("Initial chi2 /ntoas is %.8f, and chi2 = %.8f\n", chi2/ntoas, chi2);

        printf("\nParameter map : [");
        for (int i=0; i < Get_nfitted_params() ; i++) printf("%i, ", Get_list_of_fitted_parameters(i));
        printf("]\n");
    };

    void Get_parameter_name(int paramnumber, char * name)
    {

        strcpy(name, parameters.fitparams_names[parameters.fitted_parameters[paramnumber]]);
    }

    double operator()(double * relativeshift)
    {
        return lnposterior(relativeshift);
    };

    double operator()(double * relativeshift, const double betatemp)
    {

        int iefac_abs = parameters.absolute_parameter_map[36];
        int iefac_fitted = parameters.Get_fitted_parameter_index(36);

        if (betatemp != 1. and iefac_fitted > -1 ) relativeshift[iefac_fitted] = (ref_efac - parameters.parameters_ini[iefac_abs])/parameters.scale_efac; // make sure efac remains constant for temperatures other than 1
//         printf("\nreeeeeeeeeeeeeeeeeee betatemp %.2f iefac %i %i ref_efac %.5Lf ini %.5Lf scale %.5Lf rshift %.5f \n", betatemp, iefac_abs, iefac_fitted, ref_efac, parameters.parameters_ini[iefac_abs], parameters.scale_efac,  relativeshift[iefac_fitted]);
        return lnposterior(relativeshift) * betatemp;
    };
    
    double operator()() // Return lnposterior at initial parameter values
    {
        int ndim = Get_ndim();
        double shifts[ndim];
        for (int i = 0; i < ndim; i++) shifts[i] = 0.;
        return operator()(shifts);

    }
    
    
    bool is_valid(double * relativeshift)
    // Return true if the proposed parameter set is valid, false otherwise (e.g. infinite prior)
    {
        vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
        Set_fitted_parameter_relativeshifts(vectrelativeshift);
        if (parameters.Mo < 0)
        {
            return false;
        }
        else 
        {
            return true;
        }
    };

    double lnposterior(double * relativeshift)
    {
        counter +=1;
        
        vector<double> vectrelativeshift(relativeshift, relativeshift + parameters.fitted_parameters.size() );
        Set_fitted_parameter_relativeshifts(vectrelativeshift);
        double ragauss = -0.5*pow((parameters.RA - ramean)/(raerror*posprior),2);
        double decgauss = -0.5*pow((parameters.DEC - decmean)/(decerror*posprior),2);
        double rvgauss = -0.5 * pow((parameters.distance1* radmasdeg * parameters.distance *clight - rvkaplan)/(rverrorkaplan*rvprior),2);
        double dgauss = -0.5 * pow((parameters.distance - dmean)/(derror*dprior),2);
        double pmragauss = -0.5 * pow((parameters.RA1 - pmramean)/(pmraerror*tpmprior),2);
        double pmdecgauss = -0.5 * pow((parameters.DEC1 - pmdecmean)/(pmdecerror*tpmprior),2);
        

       // double a = parameters.apsinii*parameters.apsinii + parameters.apcosii*parameters.apcosii;
        //double sini = parameters.apsinii /a;
        double ai_i = 0.; //log( sini * nctei / a ); // nctei just to provide normalisation
        //a = parameters.aBsinio*parameters.aBsinio + parameters.aBcosio*parameters.aBcosio;
        //sini = parameters.aBsinio /a;
        double ai_o = 0.; //log( sini * ncteo / a ); // ncteo just to provide normalisation

        // Prior for quadrupole parameter
         double quadrupolepositive = 0.;
//         if (parameters.Get_fitted_parameter_index(35) > -1 ) // check parameter is fitted
//         {
//             quadrupolepositive += parameters.quadrupole[0] * 1.e9; // enable to favor exploration of larger moments
//             //             if (parameters.quadrupole[0] <0.)
// //                 quadrupolepositive += -DBL_MAX;
// //             else
// //                 quadrupolepositive += 0.;// 6.907 - log(0.0000000001+parameters.quadrupole[0]); // Modified Jeffrey's prior . Rmq : log(0.001) ~ -6.907
//          }
//
         // THIS IS NOT WORKING FOR SOME REASON motion_changed is always "false" 
//         if (parameters.motion_changed == false and counter > 2)
//         {   
//             motion_changed = false;
//             printf("\n MOTION CHANGED FALSE %d ", counter);
//         }
//         else
//         {
//             motion_changed = true ;
//             
//         }
                    
        lastlnlikelyhood = Compute_lnposterior(0) * ntoas * beta;
        lastprior = ragauss + decgauss + rvgauss + dgauss + pmragauss + pmdecgauss + ai_i + ai_o + quadrupolepositive;


        if (verbose == 1) printf("lnlikelyhood %.3e ragauss %.3e decgauss %.3e rvgauss %.3e dgauss %.3e pmragauss %.3e pmdecgauss %.3e ai_i %.3e ai_o %.3e quadrupolepositive %.3e\n", lastlnlikelyhood , ragauss , decgauss , rvgauss , dgauss, pmragauss, pmdecgauss, ai_i, ai_o, quadrupolepositive);

        return (lastlnlikelyhood + lastprior) ;
    };

    double Swap_temperatures(double * params_lowtemp, double & lnposterior_lowtemp, const double beta_lowtemp, double * params_hightemp, double & lnposterior_hightemp, const double beta_hightemp )
    {
        const int ndim = Get_ndim();
        double walkerbuffer[ndim];
        double lnposteriory = lnposterior_lowtemp;
        double efac_temp1 = 0.;
        double efac_lowtemp = 1.;
        double efac_hightemp = 1.;
        double logefac =0.;


        int iefac_abs = parameters.absolute_parameter_map[36];
        int iefac_fitted = parameters.Get_fitted_parameter_index(36);

        if (iefac_fitted >= 0)
        {
            efac_temp1 = params_lowtemp[iefac_fitted];
            efac_lowtemp = (params_lowtemp[iefac_fitted] * parameters.scale_efac)  + parameters.parameters_ini[iefac_abs];
            efac_hightemp = (params_hightemp[iefac_fitted] * parameters.scale_efac)  + parameters.parameters_ini[iefac_abs];
            logefac = ntoas* log(efac_lowtemp / efac_hightemp);
        }

        lnposterior_lowtemp = lnposterior_hightemp * beta_lowtemp / beta_hightemp;
        lnposterior_hightemp = lnposteriory * beta_hightemp / beta_lowtemp;

        memcpy(walkerbuffer, params_lowtemp, ndim*sizeof(double));
        memcpy(params_lowtemp, params_hightemp, ndim*sizeof(double));
        memcpy(params_hightemp, walkerbuffer, ndim*sizeof(double));

        if (beta_lowtemp == 1. and iefac_fitted >= 0)
        {
            lnposterior_lowtemp += logefac;
            lnposterior_hightemp -= beta_hightemp * logefac;
            params_hightemp[iefac_fitted] = (ref_efac - parameters.parameters_ini[iefac_abs])/parameters.scale_efac;
            params_lowtemp[iefac_fitted] = efac_temp1;
        }
    };

    double Get_chi2(const double lnposterior, const double * params, const double beta)
    // Return the chi2 based on the log of the posterior probability density, therefore removing the Ntoas*ln(efac) contribution.
    // Does not recompute the lnposterior : fast
    // Remove effect of temperature. ( i.e. gives chi2 with temperature = 1)
    {
        int iefac_abs = parameters.absolute_parameter_map[36];
        int iefac_fitted = parameters.Get_fitted_parameter_index(36);
        double currentefac = ref_efac;

        if (iefac_fitted >=0 and beta ==1. ) currentefac = (params[iefac_fitted] * parameters.scale_efac)  + parameters.parameters_ini[iefac_abs];

        return -2.*(lnposterior/beta + ntoas*log(currentefac));
    };


    double initial_distribution(mt19937& random_generator, int param_number)
    {
        normal_distribution<double> initial_distribution(0.,0.0001);
        int iquadrupole1 = parameters.Get_fitted_parameter_index(35);
        if (param_number == iquadrupole1)
            return abs(initial_distribution(random_generator));
        else
            return initial_distribution(random_generator);
    }
    

};


#endif
