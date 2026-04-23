# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#coding:utf8
#/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 #*
 #* Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 #*/

import corner as cn
import numpy as np
from pylab import *
from acor_python import acorpython
from outils import Print_numbers_with_errors_latex, Print_numbers
from scipy.stats import ks_2samp
ioff()

names = ["Spin frequency",
                 "Spin frequency derivative",
                 "Laplace-Lagrange sin",        # inner system
                 "Semimajor axis line of sight",
                 "Laplace-Lagrange cos",
                 "Semimajor axis plane of sky",
                 "Time of ascending node",
                 "Orbital period",
                 "Laplace-Lagrange sin",        # outer system
                 "Semimajor axis line of sight",
                 "Laplace-Lagrange cos",
                 "Semimajor axis plane of sky",
                 "Time of ascending node",
                 "Orbital period",
                 "Pulsar mass",
                 "Inner companion mass",
                 "Outer companion mass",
                 "Longitude of outer ascending node",
                 "Difference in long. of asc. nodes",
                 "dphase0",
                 r"SEP $\Delta$",
                 "Right ascension",
                 "Declination",
                 "Distance",
                 "Right-ascension proper motion",
                 "Declination  proper motion",
                 "Distance  proper motion",
                 "Dispersion measure",
                 "Dispersion measure variation"
                 ]

#
# symbols = [r"$f$",                      #0
#                    r"$f'$",
#                    r"$e_I\sin \omega_I$",
#                    r"$a_p\sin i_I$",            #3
#                    r"$e_I\cos\omega_I$",
#                    r"$a_p\cos i_I$",            #5
#                    r"${t_{\mathrm{asc}}}_I$",
#                    r"$P_I$",
#                    r"$e_O\sin \omega_O$",       #8
#                    r"$a_b\sin i_O$",
#                    r"$e_O\cos\omega_O$",
#                    r"$a_b\cos i_O$",            #11
#                    r"${t_{\mathrm{asc}}}_O$",
#                    r"$P_O$",
#                    r"$(m_p+m_i)/2$",                    #14
#                    r"$(m_p-m_i)/2$",#r"$(m_p-m_i)/2$",
#                    r"$m_o$",
#                    r"$\Omega_o$",                #17
#                    r"$\delta\Omega$",           #18
#                    r"$\delta\phi_0$",
#                    r"$\Delta$",
#                    r"$\alpha$",                 #21
#                    r"$\delta$",
#                    r"$d$",
#                    r"$\mu_\alpha$",             #24
#                    r"$\mu_\delta$",
#                    r"$\mu_d$",
#                    r"$\mathrm{DM}$",            #27
#                    r"$\mathrm{DM}'$"
#                    ]


unitstr= [r"\mathrm{s}^{-1}", r"\mathrm{s}^{-2}",
                  r"0",r"\mathrm{ls}",r"0",r"\mathrm{ls}",r"\mathrm{d}",r"\mathrm{d}",
                  r"0",r"\mathrm{ls}",r"0",r"\mathrm{ls}",r"\mathrm{d}",r"\mathrm{d}",
                  r"M_{\odot}",r"M_{\odot}",r"M_{\odot}",r"\degree", r"\degree",r"0", r"0",r"0",r"0",r"\mathrm{ly}",
                  r"\mathrm{mas}/\mathrm{yr}",r"\mathrm{mas}/\mathrm{yr}", r"\mathrm{mas}/\mathrm{yr}",
                  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}",  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}\mathrm{yr}^{-1}",
                  ]

names = ["Spin frequency",
                 "Spin frequency derivative",
                 "Laplace-Lagrange sin",        # inner system
                 "Semimajor axis line of sight",
                 "Laplace-Lagrange cos",
                 "Semimajor axis plane of sky",
                 "Time of ascending node",
                 "Orbital period",
                 "Laplace-Lagrange sin",        # outer system
                 "Semimajor axis line of sight",
                 "Laplace-Lagrange cos",
                 "Semimajor axis plane of sky",
                 "Time of ascending node",
                 "Orbital period",
                 "Pulsar mass",
                 "Inner companion mass",
                 "Outer companion mass",
                 "Longitude of outer ascending node",
                 "Difference in long. of asc. nodes",
                 "dphase0",
                 r"SEP $\Delta$",
                 r"SEP $\bar{\gamma}_{0p}$",
                 r"SEP $\bar{\beta}_{0pp}$",
                 r"SEP $\bar{\beta}_{p00}$",
                 r"SEP $\bar{\beta}_{0p0}$",
                 "Right ascension",
                 "Declination",
                 "Distance",
                 "Right-ascension proper motion",
                 "Declination  proper motion",
                 "Distance  proper motion",
                 "Dispersion measure",
                 "Dispersion measure variation",
                 #"DMX", # (33)
                 #"FD",  # (34)
                 "Quadrupole moment", # (35)
                 "EFAC",
                 r"$\Delta_i$" # 37
                 ]

symbols = [r"$f$",                      #0
                   r"$f'$",
                   r"$e_I\sin \omega_I$",
                   r"$a_p\sin i_I$",            #3
                   r"$e_I\cos\omega_I$",
                   r"$a_p\cos i_I$",            #5
                   r"${t_{\mathrm{asc}}}_I$",
                   r"$P_I$",
                   r"$e_O\sin \omega_O$",       #8
                   r"$a_b\sin i_O$",
                   r"$e_O\cos\omega_O$",
                   r"$a_b\cos i_O$",            #11
                   r"${t_{\mathrm{asc}}}_O$",
                   r"$P_O$",
                   r"$(m_p+m_i)/2$",                    #14
                   r"$(m_p-m_i)/2$",#r"$(m_p-m_i)/2$",
                   r"$m_o$",
                   r"$\Omega_o$",                #17
                   r"$\delta\Omega$",           #18
                   r"$\delta\phi_0$",
                   r"$\Delta$",
                   r"$\bar{\gamma}_{0p}$",
                   r"$\bar{\beta}_{0pp}$",
                   r"$\bar{\beta}_{p00}$",
                   r"$\bar{\beta}_{0p0}$",
                   r"$\alpha$",                 #25
                   r"$\delta$",
                   r"$d$",
                   r"$\mu_\alpha$",             #28
                   r"$\mu_\delta$",
                   r"$\mu_d$",
                   r"$\mathrm{DM}$",            #31
                   r"$\mathrm{DM}'$",           #32
                   #r"$\mathrm{DMX}$",           #(33)
                   #r"$\mathrm{FD}$",            #(34)
                   r"$Q$",                      #33(35)
                   "EFAC",                      #34 (36)
                   r"$\Delta i$"                #   (37)
                   ]

symbols2 = [r"$\bar{f}$",                      #0
                   r"$\bar{f'}$",
                   r"$e_I\sin \omega_I$",
                   r"$a_p\sin i_I$",            #3
                   r"$e_I\cos\omega_I$",
                   r"$a_p\cos i_I$",            #5
                   r"${t_{\mathrm{asc}}}_I$",
                   r"$P_I$",
                   r"$e_O\sin \omega_O$",       #8
                   r"$a_b\sin i_O$",
                   r"$e_O\cos\omega_O$",
                   r"$a_b\cos i_O$",            #11
                   r"${t_{\mathrm{asc}}}_O$",
                   r"$P_O$",
                   r"$m_i/m_p$",                    #14
                   r"$(m_p-m_i)/2$",#r"$(m_p-m_i)/2$",
                   r"$m_o$",
                   r"$\Omega_o$",                #17
                   r"$\delta\Omega$",           #18
                   r"$\delta\phi_0$",
                   r"$\Delta$",
                   r"$\bar{\gamma}_{0p}$",
                   r"$\bar{\beta}_{0pp}$",
                   r"$\bar{\beta}_{p00}$",
                   r"$\bar{\beta}_{0p0}$",
                   r"$\alpha$",                 #25
                   r"$\delta$",
                   r"$d$",
                   r"$\mu_\alpha$",             #28
                   r"$\mu_\delta$",
                   r"$\mu_d$",
                   r"$\mathrm{DM}$",            #31
                   r"$\mathrm{DM}'$",           #32
                   #r"$\mathrm{DMX}$",           #(33)
                   #r"$\mathrm{FD}$",            #(34)
                   r"$Q$",                      #33(35)
                   "EFAC",                      #34 (36)
                   r"$\delta i$",               #   (37)
                   r"$e_x \sin \omega_x$",            # 38      
                   r"$e_x \cos \omega_x$",          # 39
                   r"$a_x \sin i_x$",          # 40      
                   r"$a_x \cos i_x$",          # 41      
                   r"${t_{\mathrm{asc}}}_x$",           # 42      
                   r"$P_x$",              # 43  
                   r"$\Omega_x$"            # 44  
                   ]

unitstr= [r"\mathrm{s}^{-1}", r"\mathrm{s}^{-2}",
                  r"0",r"\mathrm{ls}",r"0",r"\mathrm{ls}",r"\mathrm{d}",r"\mathrm{d}",
                  r"0",r"\mathrm{ls}",r"0",r"\mathrm{ls}",r"\mathrm{d}",r"\mathrm{d}",
                  r"M_{\odot}",r"M_{\odot}",r"M_{\odot}",r"\degree", r"\degree",r"0",
                  r"0", r"0",r"0", r"0",r"0", # SEP
                  r"0",r"0",r"\mathrm{ly}",   # Position
                  r"\mathrm{mas}/\mathrm{yr}",r"\mathrm{mas}/\mathrm{yr}", r"\mathrm{mas}/\mathrm{yr}",
                  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}",  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}\mathrm{yr}^{-1}",
                  r"0", #(33)
                  r"0", #(34)
                  r"0", #(35)
                  r"0",
                  r"\degree", # (37)
                  r"0", # 38
                  r"0", # 39
                  r"ls", # 40
                  r"ls", # 41
                  r"d", # 42
                  r"d", # 43
                  r"\degree", # 44
                  ]


def funcprior(params):
    raddeg = np.arccos(-1.)/180. # radians per degree
    radmasdeg = raddeg / (3600. * 1000.) # radian per mas of degree
    # Position from GAIA DR2 (radians) Warning : POSEPOCH = MJD2015.5
    ramean = 9.500283096783616e-01
    decmean = 3.011411347204199e-01
    raerror = 8.945972920351486e-10 # =  0.184523937110807 mas
    decerror = 9.302791001352137e-10 # = 0.191883838345113 mas
    posprior = 2

    # Tangential proper motion (mas/yr) from Gaia DR2. Warning : POSEPOCH = MJD2015.5
    pmramean = 4.81377140723352
    pmdecmean = -4.42182046193475
    pmraerror = 0.498024673422251
    pmdecerror = 0.42770795176738
    tpmprior = 2

    # Radial velocity from Kaplan et al. 2014
    rvkaplan = 29.7e+03 #m/s
    rverrorkaplan = 0.9e+03 # m/s
    rvprior = 1

   # Arbitrary prior centered between Gaia DR2 and Kaplan 2014 with ~2error kaplan. Warning : POSEPOCH = MJD2015.5
    dmean = 4350.#4604.639595490277  # lyr
    derror = 500. #1558.9119367898795 # lyr
    dprior = 1


    RA = params[25]
    DEC = params[26]
    distance = params[27]
    RA1 = params[28]
    DEC1 = params[29]
    distance1 = params[30]

    ao2 = params[9]**2 + params[11]**2
    ai2 = params[3]**2 + params[5]**2
    sinii = params[3] /ai2
    sinio = params[9] /ao2
    priors=np.zeros(params.size)
    priors[0] = np.log(sinii /ai2 *  911596921.5018525)
    priors[1] = np.log(sinio /ao2 * 55954106240.898094)
    ragauss = -0.5*((RA - ramean)/(raerror*posprior))**2
    decgauss = -0.5*((DEC - decmean)/(decerror*posprior))**2
    rvgauss = -0.5 * ((distance1* 1.4534348513185598 * distance - rvkaplan)/(rverrorkaplan*rvprior))**2 # clight * radmasdeg = 1.4534348513185598
    dgauss = -0.5 * ((distance - dmean)/(derror*dprior))**2
    pmragauss = -0.5 * ((RA1 - pmramean)/(pmraerror*tpmprior))**2
    pmdecgauss = -0.5 * ((DEC1 - pmdecmean)/(pmdecerror*tpmprior))**2

    priors[2] = ragauss
    priors[3] = decgauss
    priors[4] = dgauss
    priors[5] = pmragauss
    priors[6] = pmdecgauss

    return priors


paramnames = [symbols[i] + " ($" + unitstr[i] + "$) " for i in range(len(symbols))]
### mcmc-testnewmodel0_mcmc-288x200000d50.res
#parameter_map = [21, 22, 24, 25, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 20]

### mcmc-testnewmodel0posransom_mcmc-192x1000000d50
#parameter_map = [24, 25, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 20]

from astropy.modeling import models,fitting
from scipy.special import erf

def Fit_and_plot_distribution(sample, astropy_law, astropy_law_args={}, bins=50, xlabel='Value', ylabel='Number of occurences of value', write_legend=None, write_legend_args=(), loc_legend=0, fontsize='x-large', figsize=(8,6)):
    # Create the Cumulative Distribution Function to fit on
    sampleordered = np.array(sample)
    sampleordered.sort()
    cdf = np.zeros(sample.size)
    x = np.zeros(sample.size)

    j=0
    x[j] = sampleordered[0]
    cdf[j] = 1
    for i in range(1, sample.size):
        if (sampleordered[i] == sampleordered[i-1]):
            cdf[j] += 1.
        else :
            j += 1
            x[j] = sampleordered[i]
            cdf[j] +=1 #= cdf[j-1] + 1.
    x = x[:j+1]
    cdf = cdf[:j+1]
    for i in range(1, cdf.size):
        cdf[i] = cdf[i] + cdf[i-1]
    cdf /= cdf[-1]
    # Initialise the parameters of the fit ( ultimately should provide an external routine to do that)

    # Do the Fit
    fitc = fitting.LevMarLSQFitter()
    model = astropy_law(**astropy_law_args)
    fitres = fitc(model, x, cdf)

    # Plot the fit
    fig = figure(figsize=figsize)
    plt = fig.add_subplot(111)
    plt.plot(x, cdf, '.')
    plt.plot(x, fitres(x))
    fig.show()

    # Plot the fitted histogram of the Distribution.
    fig = figure(figsize=figsize)
    plt = fig.add_subplot(111)
    counts, bins, patches = plt.hist(sample, bins=bins)
    if write_legend is not None:
        legend_fit= write_legend(fitres, *write_legend_args)
    else :
        legend_fit = ''
        for i in range(fitres.parameters.size):
            legend_fit += fitres.param_names[i] + '=' + '{:.2e} '.format(fitres.parameters[i])
    plt.plot(bins[1:] - 0.5*(bins[1] - bins[0]), (fitres(bins[1:]) - fitres(bins[:bins.size -1]))*counts.sum(), label=legend_fit)

    #plt.ticklabel_format(style='sci', scilimits=(1,1), axis='y')
    formatter = matplotlib.ticker.ScalarFormatter(useMathText='True')
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,1))
    plt.xaxis.set_major_formatter(formatter)
    plt.yaxis.set_major_formatter(formatter)

    plt.set_xlabel(xlabel)
    plt.set_ylabel(ylabel)
    plt.xaxis.label.set_fontsize(fontsize)
    plt.yaxis.label.set_fontsize(fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.legend(loc=loc_legend, fontsize=fontsize)

    fig.set_tight_layout(True)

    fig.show()

    return fig, fitres, cdf


def erf_distribution(x,mean=0., std=1.):
    return 0.5* ( erf((x-mean)/(std*np.sqrt(2.))) + 1.)
Gaussian1D_cdf = models.custom_model(erf_distribution)

def doubleerf_distribution(x,mean1=0., std1=1., mean2=0., std2=1., ratio=0.5):
    return ratio *0.5* ( erf((x-mean1)/(std1*np.sqrt(2.))) + 1.) + (1.-ratio)*0.5* ( erf((x-mean2)/(std2*np.sqrt(2.))) + 1.)
doubleGaussian1D_cdf = models.custom_model(doubleerf_distribution)

def write_Gaussian1D_legend(model, quantity_name='X'):
    '''
        model : any astropy model with parameters "mean" and "std"
    '''
    return r'$\left<' + quantity_name + r'\right>' +'=' + Print_numbers(model.mean.value,chiffre_significatif=2, limite_notasci=3) + r'; \sigma_{' + quantity_name + r'} =' + Print_numbers(model.std.value,chiffre_significatif=2, limite_notasci=3) + '$'

class mcmcresult :
    def __init__(self,chains=None, chi2=None, resfile=None, parameter_object = None, parameter_map=None, parameter_scales=None, parameter_initials=None, symbols=None, derived_parameters=None, print_parameters=None, res_dtype=np.float):
        '''
            parameter_object is an object with methods "Get_fitted_parameter_map", "Get_parameter_scales", "Get_fitted_parameter_names", "derive_parameters", "Get_parameters".
            If parameter_map, parameter_scales, symbols, derived_parameters, parameter_initials are given the values obtained from parameter_object are overwritten.

            If parameter_scales is provided, the chains are rescaled.

            If a parameter map is provided, the chain is NOT shuffled according to it.

        '''
        self.symbols = None
        if resfile is not None :
            res = np.loadtxt(resfile, comments='#', dtype = res_dtype)
            self._ndim = res.shape[1] - 1
            self.res = np.ndarray((res.shape[0], self._ndim))
            self.chi2 = res[:,-1]
            self.res = res[:,:self._ndim]

            # Now look for parameter names
            file = open(resfile,"r")
            line =file.readline()
            while (line != '') :
                sline = line.split()
                if (sline[0] =='#' and sline[-1] == 'lnposterior') : break
                line = file.readline()
            file.close()
            symcandidate=[sline[1]]
            if (len(sline) >= self._ndim + 2) :
                j = 0
                for i in range(2,len(sline)-1):
                    if sline[i][0] == '(' :
                        symcandidate[j] += ' ' + sline[i]
                    else :
                        symcandidate.append(sline[i])
                        j += 1
                if len(symcandidate) == self._ndim:
                    print("Parameter names found in result file.")
                    self.symbols = symcandidate
                else :
                    print('Warning : It seems that the result file contains parameter names but they could not be extracted')
        else :
            if type(chains) is dict : 
                self.__from_dict(chains, chi2)
            else:
                self.res = chains.copy()
                self.chi2 = chi2.copy()
                self._ndim = self.res.shape[1]

        if parameter_object is not None :
            self.print_parameters = parameter_object.Print_parameters_latex
            self.derive_parameters = parameter_object.Derive_parameters
            self.parammap = np.array(parameter_object.Get_fitted_parameter_map())
            self.paramscales = np.array(parameter_object.Get_parameter_scales())
            self.symbols = parameter_object.Get_parameter_names()
            self.paraminitials = np.array(parameter_object.Get_parameters())
            print(self.parammap)
            print(self.paramscales)
            print(self.symbols)

        if derived_parameters is not None :
            self.derive_parameters = derived_parameters
        elif parameter_object is None:
            self.derive_parameters = None

        if print_parameters is not None :
            self.print_parameters = print_parameters
        elif parameter_object is None:
            self.print_parameters = None

        if parameter_map is not None :
            self.parammap = np.array(parameter_map)
            if len(parameter_map) < self._ndim :
                print("Error : parameter map too short ! %i < %i"%( len(parameter_map), self._ndim))
        elif parameter_object is None :
            self.parammap = np.arange(self._ndim)

        if symbols is not None :
            self.symbols = symbols[0: max(self.parammap)+1]
        elif (parameter_object is None) and (self.symbols is None) :
            self.symbols = ["{:d}".format(i) for i in range(max(self.parammap)+1)]

        #if names is not None :
            #self.names = names[0: max(self.parammap)+1]
        #elif parameter_object is None :
            #self.names = ["{:i}".format(i) for i in range(max(self.parammap)+1)]

        if parameter_scales is not None:
            self.paramscales = np.array(parameter_scales)
        elif parameter_object is None :
            self.paramscales = np.ones(max(self.parammap)+1)

        if parameter_initials is not None:
            self.paraminitials = np.array(parameter_initials)
        elif parameter_object is None :
            self.paraminitials = np.zeros(max(self.parammap)+1)

        for i in range(self._ndim) :
            self.res[:, i] *= self.paramscales[self.parammap[i]]
            self.res[:, i] += self.paraminitials[self.parammap[i]]


        # Create inverse parameter map
        self.inversparammap = np.ones(self.parammap.max()+1,dtype=int) * (self.parammap.size ) # the multiplication factor ensures that any call to a parameter that is not in the chain will provoke an error.
        self.inversparammap[self.parammap] = np.arange(self.parammap.size, dtype=int)

        #if parameter_map is not None :
           #for i in range(self._ndim) :
               #self.res[:, parameter_map[i]] = res[:, i ]
        #else :
            #self.res = res[:, :ndim]


        return
    
    def __from_dict(self, chains_dict, chi2 = None):
        '''
        Initialise from a chain past as dictionary.
        '''
        self.symbols = list(chains_dict.keys())
        self._ndim = len(self.symbols)
        nsamp = len(chains_dict[self.symbols[0]])
        self.res = np.zeros((nsamp,self._ndim))
        for i,k in enumerate(self.symbols):
            self.res[:,i] = chains_dict[k]
        if chi2 is not None : 
            print("bla ", chi2)
            self.chi2 = chi2.copy()
        else: 
            self.chi2 = np.zeros(nsamp)
        return 
        

    def Get_keys(self):
        '''
            Return keys taken from self.symbols given the current self.parammap, 
            that is to say corresponding to the order of the samples in self.res.
            These keys are used in Get_dictionary and are the default labels in plots.
        '''
        return [self.symbols[self.parammap[i]] for i in range(self._ndim) ]
        
    def Get_dictionary(self, start=0, end =None, return_keys=False):
        '''
            return dico
            
            Return the sample as a dictionary the key of which are taken from self.Get_keys.
            The sample is restricted to the range start:end (negative values allowed for end)
        '''
        names = self.Get_keys()
        dico = {}
        for i in range(self._ndim):
            dico[names[i]] = self.res[start:end, i]
        return dico
            
    def Mean_parameters(self, start=0, end=None):
        '''
            Give an array of the chain averaged between the iterations "start" and "end".
        '''
        mean = np.zeros(self._ndim)
        for i in range(self._ndim):
            mean[i] = self.res[start:end,i].mean()
        return mean


    def Errors_derived_parameters(self, start=0, end=-1, levelpercentage=95.):
        '''
            Compute the derived chain using the derive_parameters method and compute there errors as in the Errors method.
        '''
        if self.derive_parameters is None :
            return np.zeros(self._ndim), np.zeros(self._ndim)

        if end < 0 :
            npts = self.res.shape[0] + 1 + end - start
        else :
            npts = end - start

        testsize = self.derive_parameters(self.res[start, :], external_param_map= self.parammap)
        derivedchain = np.ndarray((npts, testsize.size))

        for i in range(npts ):
            derivedchain[i] = self.derive_parameters(self.res[start+i, :self._ndim], external_param_map= self.parammap)

        means = np.zeros(derivedchain.shape[1])
        for i in range(derivedchain.shape[1]):
            means[i] = derivedchain[:,i].mean()
            #print(means[i])

        errorplus = np.ndarray(testsize.size)
        errorminus = np.ndarray(testsize.size)
        for  i in range(testsize.size):
            print(i, testsize.size, derivedchain[:,i])
            halfdistro = np.compress(derivedchain[:,i] >  means[i],  derivedchain[:,i] - means[i])
            halfdistro.sort()
            if halfdistro.size ==0:
                print('Warning : derived parameter {:} has no dispersion.'.format(i))
                errorplus[i]=0.
            else :
                errorplus[i]  = halfdistro[np.int(np.floor(levelpercentage/100. * halfdistro.size -1))]


            halfdistro = np.compress(derivedchain[:,i] < means[i],  - derivedchain[:,i] + means[i])
            #print(halfdistro.size)
            halfdistro.sort()
            if halfdistro.size ==0:
                print('Warning : derived parameter {:} has no dispersion.'.format(i))
                errorminus[i]=0.
            else :
                errorminus[i]  = halfdistro[np.int(np.floor(levelpercentage/100. * halfdistro.size-1))]


        return errorplus, errorminus


    def Errors(self, start=0, end=-1, levelpercentage=95.):
        '''
            return errorplus, errorminus

            Estimate the value errorplus/errorminus of each quantity in the chain for which levelpercentage/100 of the sample resp above/below its mean is resp below/above.
        '''
        means = self.Mean_parameters(start=start, end=end)
        errorplus = np.ndarray(self._ndim)
        errorminus = np.ndarray(self._ndim)
        for  i in range(self._ndim):
            halfdistro = np.compress(self.res[start:end,i] > means[i],  self.res[start:end,i] - means[i])
            halfdistro.sort()
            errorplus[i]  = halfdistro[np.int(np.floor(levelpercentage/100. * halfdistro.size -1))]

            halfdistro = - np.compress(self.res[start:end,i] < means[i],  self.res[start:end,i] - means[i])
            halfdistro.sort()
            errorminus[i]  = halfdistro[np.int(np.floor(levelpercentage/100. * halfdistro.size-1))]

        return errorplus, errorminus


    def Standard_deviations_half(self, start=0, end=-1):
        '''
            return sigmasplus, sigmasminus

        '''
        means = self.Mean_parameters(start=start, end=end)
        sigmasplus = np.ndarray(self._ndim)
        sigmasminus = np.ndarray(self._ndim)

        for  i in range(self._ndim):
            halfdistro = np.compress(self.res[start:end,i] > means[i],  self.res[start:end,i] - means[i])
            sample = np.ndarray(halfdistro.size*2)
            sample[:halfdistro.size] = halfdistro
            sample[halfdistro.size:] = -halfdistro
            sigmasplus[i] = sample.std()

            halfdistro = - np.compress(self.res[start:end,i] < means[i],  self.res[start:end,i] - means[i])
            sample = np.ndarray(halfdistro.size*2)
            sample[:halfdistro.size] = halfdistro
            sample[halfdistro.size:] = -halfdistro
            sigmasminus[i] = sample.std()

        return sigmasplus, sigmasminus

    def Standard_deviations(self, start=0, end=-1):
        '''
            return sigmasplus, sigmasminus

        '''
        sigs = np.zeros(self._ndim)
        for  i in range(self._ndim):
            sigs[i]=self.res[start:end,i].std()

        return sigs


    def Covariance_matrix(self, start=0, end=-1, normalized=False):
        '''
        return np.cov(samples), samples
        where "samples" is np.transpose(self.res ) with the mean of each parameter removed.

        normalized: if True, each parameter is normalized by its standard deviation so that the diagonal elements
                    of the covariance matrix are 1.
        '''
        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        samples = (self.res[start:maxit]).copy()
        samples = np.transpose(samples)

        for i in range(self._ndim):
            samples[i] -= samples[i].mean()
            if normalized is True:
                sig = samples[i].std()
            else:
                sig =1.
            samples[i] /= sig

        return np.cov(samples), samples


    def Plot_Covariance_matrix(self):

        return

    def Decorrelate_chain(self, start=0, end=-1, normalized=False):
        '''
        return decorrelated_chain
        '''
        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end
        decorchain = np.zeros((maxit - start, self.res.shape[1]))
        cov, samp = self.Covariance_matrix(start=start, end=end, normalized=normalized)
        samp = samp.transpose()
        eigvals, eigvect = np.linalg.eig(cov)
        if normalized :
            for i in range(eigvals.size):
                eigvect[:,i] /= np.sqrt(eigvals[i])
        decormat = np.linalg.inv(cov)
        for i in range(start, maxit):
            decorchain[i-start] = np.dot(decormat, samp[i-start])
        return decorchain.copy()

    def Autocorrelation_time(self, start=0, end=-1) :
        '''
            autocorrelation times, standard deviations = Autocorrelation_time(start=0, end=-1)
            start and end give the range of steps to take into account in the computation of autocorrelation times
        '''
        autocor = np.ndarray(self._ndim)
        sigmas = np.ndarray(self._ndim)
        for i in range(self._ndim) :
            autocor[i],mean, sigmas[i]  = acorpython(self.res[start:end, i])
        return autocor, sigmas

    def Corner_plot(self, start=0, end=None, show_means=True, figsize=None,
                    tight_layout=False, show_tick_labels=True, rotationx=0, rotationy=0,
                    subplotadjust=(0.05,0.05, 0.99,0.99,0.,0.), labels=None, submap=None):
        '''
            Create a corner plot using the chain between iterations "start" and "end"
            
            end : if None take the whole chain
            submap : list of keys to use in self.Get_dictionary() to plot a subcorner
            labels : overwrites default labels
        '''
        if show_means:
            truths = self.Mean_parameters(start=start, end = end)
        else:
            truths=None
        
        if end is None : 
            end = self.chi2.size
        
        if submap is not None:
            allkeys = self.Get_keys()
            dico = self.Get_dictionary()
            npar = len(submap)
            npts = self.chi2.size
            npts = (end < 0)*(npts+end) + (end>=0)*end - start
            sample = np.zeros((npts, npar))
            submap_idx =[] # map of indices
            for i in range(npar):
                sample[:,i] = dico[submap[i]][start:end]
                j=0
                while (submap[i] != allkeys[j]):
                    j+=1
                submap_idx.append(j)
            submap_idx = np.array(submap_idx)
            if truths is not None:
                truths = truths[submap_idx]
        else:
            sample = self.res[start:end]
            
        if labels is None : 
            if submap is None : 
                labels = self.Get_keys()
            else:
                labels = submap
            
        fig = cn.corner(sample, labels = labels, label_kwargs={'fontsize':'x-large', },truths = truths)
        if figsize is not None:
            fig.set_size_inches(figsize)
        fig.set_tight_layout(tight_layout)
        plts = fig.get_axes()
        formatter = matplotlib.ticker.ScalarFormatter(useMathText='True')
        formatter.set_scientific(True)
        formatter.set_powerlimits((0,0))
        xmin = 10.
        ymin = 10.
        for ax in plts :
            if ax.get_position().xmin < xmin : xmin = ax.get_position().xmin
            if ax.get_position().ymin < ymin : ymin = ax.get_position().ymin

        for ax in plts :
            if ax.get_position().xmin == xmin : # plots in left column
                formatter = matplotlib.ticker.ScalarFormatter(useMathText='True')
                formatter.set_scientific(True)
                formatter.set_powerlimits((0,0))
                ax.yaxis.set_major_formatter(formatter)
                ax.set_ylabel(ax.get_ylabel(), rotation=90 - rotationy)
                ax.xaxis.label.set_in_layout(False) # for tight_layout
                if not show_tick_labels :
                    ax.set_yticklabels([])
            if ax.get_position().ymin == ymin : # plot in bottom row
                formatter = matplotlib.ticker.ScalarFormatter(useMathText='True')
                formatter.set_scientific(True)
                formatter.set_powerlimits((0,0))
                ax.xaxis.set_major_formatter(formatter)
                ax.set_xlabel(ax.get_xlabel(), rotation=rotationx)
                ax.yaxis.label.set_in_layout(False)
                if not show_tick_labels :
                    ax.set_xticklabels([])
            else:
                ax.xaxis.label.set_in_layout(False)
                ax.yaxis.label.set_in_layout(False)


        fig.set_tight_layout(tight_layout)
        fig.subplots_adjust(*subplotadjust)
            #ax.yaxis.set_major_formatter(formatter)
            #ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
        return fig



    def Walker_traces(self, walker_list, walker_number, start = 0, end =-1):
        '''
        return walkertraces, walkerindexes, walkerslimits,  walkerjumps, walkeracceptfrac

        walker_list : list of indexes of the walkers to follow
        walker_number : Total number of walkers that were present in the run
        '''
        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        # Computing walker traces
        nwalkerstep = self.res.shape[0] // walker_number
        walkertraces = np.ndarray((len(walker_list), nwalkerstep, self._ndim+1))
        walkerindexes = np.ndarray((len(walker_list), nwalkerstep), dtype=int)
        walkerslimits = np.ndarray((len(walker_list),2), dtype=int)
        walkerjumps = np.zeros((len(walker_list), nwalkerstep), dtype=int)
        walkeracceptfrac = np.zeros(len(walker_list))
        for w in range(len(walker_list)):
            ist = np.int((start - walker_list[w]) / walker_number +  ((start - walker_list[w]) / walker_number >0.) )
            ie = np.int((maxit - walker_list[w]) / walker_number)
            walkerindexes[w] = walker_list[w] + np.arange(nwalkerstep)* walker_number
            walkerslimits[w]= [ist,ie+1]
            for i in range(nwalkerstep):
                walkertraces[w, i, :self._ndim] = self.res[walkerindexes[w,i]]
                walkertraces[w, i, self._ndim] = self.chi2[walkerindexes[w,i]]
            for i in range(1,nwalkerstep):
                walkerjumps[w,i] = np.any(walkertraces[w, i, :self._ndim] != walkertraces[w, i-1, :self._ndim])*10
            walkeracceptfrac[w] = walkerjumps[w,walkerslimits[w,0]:walkerslimits[w,1]].sum()/(walkerslimits[w,1]-walkerslimits[w,0]) /10

        return walkertraces, walkerindexes, walkerslimits,  walkerjumps, walkeracceptfrac




    def Chain_plot(self, start = 0, end =-1, average_window=1, parameter_list = None, walker_list=None, walker_number=None, show_initial=True, show_mean=True,
                    show_std=True, figsize=(16,40), alpha = 0.5):
        '''
            fig = Chain_plot(start = 0, end =-1, average_window=1, parameter_list = None,walker_list=None, walker_number=None, show_initial=True, figsize=(16,40)))

            start, end : Plot the chain between step "start" and step "end".
            show_initial: If show_initial = True : Plot in blue horizontal lines the value of parameter in mcmcresult.parameter_initials
            average_window: lot in red the value of the chain averaged over a sliding window of size "average_window". If average_window <=1, do not plot
            figsize: tuple (width, height) of the size of the figure
            alpha: Transparency coefficient to use for walkers
            walker_list : list of indexes of the walkers to follow
            walker_number : Total number of walkers that were present in the run

        '''
        cols = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        if parameter_list is not None:
            paramlist = self.inversparammap[parameter_list] # convert parameter_list to chain indexes
        else :
            paramlist = np.arange(self._ndim)


        if walker_list is not None and walker_number is not None :
            walkertraces, walkerindexes, walkerslimits,  walkerjumps, walkeracceptfrac = self.Walker_traces(walker_list, walker_number, start=start, end=end)


        labels = [self.symbols[self.parammap[paramlist[i]]] for i in range(paramlist.size) ]
        fig = figure(figsize=figsize)
        naverage = np.int((maxit -start) / average_window)
        print(naverage)
        iternb = start + np.arange(maxit - start)
        aviter = start + average_window//2 + np.arange(naverage)*average_window
        avchaine = np.ndarray((paramlist.size+1, naverage))
        avchaine_std = np.ndarray((paramlist.size+1, naverage))
        means = self.Mean_parameters(start=start, end = end)
        stds = self.Standard_deviations(start=start, end = end)
        for i in range(paramlist.size):
            plt = fig.add_subplot(paramlist.size+1,1,i+1)
            plt.plot(iternb, self.res[start:maxit,paramlist[i]], ',g', alpha =alpha)
            if show_initial :
                plt.axhline(self.paraminitials[self.parammap[paramlist[i]]], linestyle = ':', color= 'b')#, start, maxit,  'b')
            if show_mean:
                plt.axhline(means[i], linestyle = '-', color= 'b')# , start, maxit,  'b')
            if show_std:
                plt.axhline((means+stds)[i], linestyle = '--', color= 'b')# , start, maxit,  'b')
                plt.axhline((means-stds)[i], linestyle = '--', color= 'b')# , start, maxit,  'b')
            if average_window > 1 :
                for j in range(naverage):
                    avchaine[i,j] = (self.res[start + average_window * j : start + average_window*(j+1) -1, paramlist[i]]).mean()
                    avchaine_std[i,j] = (self.res[start + average_window * j : start + average_window*(j+1) -1, paramlist[i]]).std()
                plt.plot(aviter, avchaine[i], '-r', label=r'$\hat{r}_{\hat{m}}'+' = {:.2f}$'.format(avchaine[i].std()/stds[i]))
                plt.plot(aviter, avchaine[i] + avchaine_std[i], '--r', label=r'$\hat{r}_{\hat{\sigma}}'+' = {:.2f}$'.format(avchaine_std[i].std()/stds[i]))
                plt.plot(aviter, avchaine[i] - avchaine_std[i], '--r')
            if walker_list is not None :
                for w in range(len(walker_list)):
                    plt.plot(walkerindexes[w, walkerslimits[w,0]:walkerslimits[w,1]], walkertraces[w, walkerslimits[w,0]:walkerslimits[w,1], paramlist[i]], '-', color=cols[w%7])
                    plt.scatter(walkerindexes[w, walkerslimits[w,0]:walkerslimits[w,1]], walkertraces[w, walkerslimits[w,0]:walkerslimits[w,1], paramlist[i]], s= walkerjumps[w, walkerslimits[w,0]:walkerslimits[w,1]], color='black')
            plt.set_ylabel(labels[i])
        #    plt.yaxis.label.set_fontsize('x-large')
            plt.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            plt.legend(loc='upper left')

        # Chi2
        plt = fig.add_subplot(paramlist.size+1,1,paramlist.size+1)
        i = paramlist.size
        plt.plot(iternb, self.chi2[start :maxit], ',g', alpha=alpha)
        chi2mean = self.chi2[start :maxit].mean()
        chi2std = self.chi2[start :maxit].std()
        if show_mean:
            plt.axhline(chi2mean, linestyle = '-', color= 'b')
        if show_std:
            plt.axhline(chi2mean+chi2std, linestyle = '--', color= 'b')# , start, maxit,  'b')
            plt.axhline(chi2mean-chi2std, linestyle = '--', color= 'b')# , start, maxit,  'b')
        if average_window > 1 :
                for j in range(naverage):
                    avchaine[i,j] = (self.chi2[start + average_window * j : start + average_window*(j+1) -1]).mean()
                    avchaine_std[i,j] = (self.chi2[start + average_window * j : start + average_window*(j+1) -1]).std()
                plt.plot(aviter, avchaine[i], '-r', label=r'$\hat{r}_{\hat{m}}'+' = {:.2f}$'.format(avchaine[i].std()/chi2std))
                plt.plot(aviter, avchaine[i] + avchaine_std[i], '--r', label=r'$\hat{r}_{\hat{\sigma}}'+' = {:.2f}$'.format(avchaine_std[i].std()/chi2std))
                plt.plot(aviter, avchaine[i] - avchaine_std[i], '--r')
        if walker_list is not None :
            for w in range(len(walker_list)):
                    plt.plot(walkerindexes[w, walkerslimits[w,0]:walkerslimits[w,1]], walkertraces[w, walkerslimits[w,0]:walkerslimits[w,1],self._ndim], '-', color=cols[w%7])
                    plt.scatter(walkerindexes[w, walkerslimits[w,0]:walkerslimits[w,1]], walkertraces[w, walkerslimits[w,0]:walkerslimits[w,1], self._ndim], s= walkerjumps[w, walkerslimits[w,0]:walkerslimits[w,1]], c='black', linestyle='solid', linewidths=3)
        plt.set_ylabel("$-\ln(P)$")
        plt.legend(loc='upper left')

        #plt.yaxis.label.set_fontsize('x-large')
        plt.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

        #fig.tight_layout(rect=(0,0,1,0.90), pad =1.11)
        fig.set_tight_layout(True)

        return fig#,walkerslimits, walkerindexes# , avchaine, aviter



    def Plot_chi2_vs_parameter(self, parameternb, chi2_chain=None, start=0, end=-1, figsize=(12,6), alpha = 0.05):
        maxpa = self.parammap.max()
        i=0
        while (i <= maxpa and self.parammap[i] != parameternb): # inverse parameter map
                i+=1
        if i > maxpa :
            print("Error : parameter number does not exist")
        paramnb = i
        print("Chain parameter number:", i, self.symbols[parameternb])

        fig = figure(figsize = figsize)
        plt = fig.add_subplot(111)
        if chi2_chain is None:
            chi2plot = self.chi2[start:end]
        else:
            chi2plot = chi2_chain[start:end]
        plt.plot(self.res[start:end, paramnb],chi2plot, '.', color = 'black', alpha=alpha)

        plt.set_xlabel(self.symbols[parameternb])
        plt.set_ylabel(r'$-\chi^2/2$')

        plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='x')
        plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='y')

        fig.tight_layout(rect=(0,0,1,0.90), pad =1.11)
        fig.set_tight_layout(True)
        return fig

    def Plot_parameter_histogram(self, parameter_number, start=0, end=-1, bins=50,
                                figsize=(12,6), show_gaussian=True, show_errorbars=True,
                                plt = None, errorbars=95., fontsize='large', legendfontsize='large'):
        '''
         plt : axes object that you can pass to draw on
        '''
        if (end < 0):
            endd = self.chi2.size + end
        else:
            endd = end
        if plt is None:
            fig = figure(figsize = figsize)
            plt = fig.add_subplot(111)
        else:
            fig=None
        distro = self.res[start:endd,self.inversparammap[parameter_number]]
        cnts, bins, patches = plt.hist(distro, bins=bins)
        plt.set_xlabel(self.symbols[parameter_number], fontsize=fontsize)
        plt.set_ylabel('Counts', fontsize=fontsize)

        if show_gaussian :
            xs = 0.5*(bins[1:] + bins[:bins.size-1])
            mean = distro.mean()
            std = distro.std()
            gausslaw = cnts.max() * np.exp(-0.5*((xs - mean)/ std)**2)
            plt.plot(xs, gausslaw, color = 'orange', label='Normal')
            plt.legend(loc=0, fontsize=legendfontsize)

        if show_errorbars:
            errplus, errminus = self.Errors(start = start, end = end, levelpercentage=errorbars)
            mean = distro.mean()
            vals = [mean - errminus[self.inversparammap[parameter_number]], mean, mean + errplus[self.inversparammap[parameter_number]]]
            plt.axvline(vals[1], linestyle=':', color='red', linewidth=1.5)
            plt.axvline(vals[2], linestyle='--', color='red')
            plt.axvline(vals[0], linestyle='--', color='red')
            arrowy = 0.1*cnts.max()
            plt.annotate("", xy=(vals[0], arrowy), xytext=(vals[2], arrowy), arrowprops=dict(arrowstyle="<->"))
            plt.text(vals[1]+ 0.05 * (vals[2] - mean), 1.05*arrowy, "{:.2g}%".format(errorbars), fontsize=fontsize)
            plttop = plt.twiny()
            plttop.set_xlim(plt.get_xlim())
            plttop.set_xticks(vals)
            plttop.set_xticklabels(['{:.2g}'.format(val) for val in vals], fontsize=fontsize)

        plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='x')
        plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='y')
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        plt.tick_params(axis='both', which='minor', labelsize=fontsize)

        if fig is not None:
            fig.tight_layout(rect=(0,0,1,0.90), pad =1.11)
            fig.set_tight_layout(True)
            self.plt = plt
            self.fig = fig
            fig.show()

        return #fig #, xs, gausslaw

    def Plot_chi2_histogram(self, start=0, end=-1, bins=50, figsize=(12,6), show_chi2_distro=True, temperature=1., return_bins=False):
        '''
        show_chi2_distro : if True, show an estimate of the corresponding chi2 distribution with the same mean

        '''
        if (end < 0):
            endd = self.chi2.size + end
        else:
            endd = end
        fig = figure(figsize = figsize)
        plt = fig.add_subplot(111)
        # // lnposteriors are provided in the chain, not chi2 !!
        chi2counts, chi2bins, patches = plt.hist(abs(self.chi2[start:endd])*2, bins=bins)

        if show_chi2_distro:
            dof = self.res.shape[1] #ndatapoints - self.res.shape[1] # number of degrees of freedom
            chi2bins = 0.5*(chi2bins[:chi2bins.size-1] + chi2bins[1:])
            chi2bins_mean = np.average(chi2bins, weights=chi2counts)
            chi2bins_mean = min(chi2bins_mean, chi2bins[0]+dof-0.1)
            chi2s = chi2bins - chi2bins_mean + dof #*temperature
            chi2distro = np.log(chi2s)*(dof/2. -1.) + (-chi2s/2.)
            chi2distro = exp(chi2distro)#exp(chi2distro- chi2distro.max())
            chi2distro *= chi2counts.max() / chi2distro.max()#[chi2counts.argmax()] # rough scaling
            plt.plot(chi2bins, chi2distro, label=r'Estimate of $\chi^2$ distribution {:d} dof'.format(dof))
            fig.legend()

        plt.set_xlabel(r'$\chi^2$')
        plt.set_ylabel('Counts')

        plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='x')
        plt.ticklabel_format(style='sci', scilimits=(-1,1), axis='y')

        fig.tight_layout(rect=(0,0,1,0.90), pad =1.11)
        fig.set_tight_layout(True)

        if return_bins:
            return fig, chi2bins, chi2counts, chi2distro
        else:
            return fig

    #def Invert_parammap(self):
        #paramnb1 = np.arange(self.parammap.size, dtype=int)
        #paramnb = np.ndarray(paramnb1.shape, dtype=int)

        #maxpa = self.parammap.max()
        #for idim in range(paramnb.size): # Invert param_nb using param_map
            #i=0
            #while (i <= maxpa and self.parammap[i] != paramnb1[idim]):
                #i+=1
            #if i > maxpa :
                #print("Error : parameter number does not exist")
            #paramnb[idim] = i

    def _Compute_2D_maps(self,p1,p2, gridres=100):
        chainindex1 = self.inversparammap[p1]
        chainindex2 = self.inversparammap[p2]
        avpixels = np.zeros((gridres, gridres))
        maxpixels = np.ones((gridres, gridres))* self.chi2.min()
        nbs = np.zeros((gridres, gridres))
        extent = [self.res[:,chainindex1].min(), self.res[:,chainindex1].max(), self.res[:,chainindex2].min(), self.res[:,chainindex2].max()]
        # range1 = np.linspace(self.res[:,chainindex1].min(), self.res[:,chainindex1].max(), num=gridres+1,endpoint=True)
        # range2 = np.linspace(self.res[:,chainindex2].min(), self.res[:,chainindex2].max(), num=gridres+1,endpoint=True)
        step1 = (extent[1] - extent[0])/gridres #(range1[-1] - range1[0]) / gridres
        step2 = (extent[3] - extent[2])/gridres # (range2[-1] - range2[0]) / gridres
        # print(range1)
        # print(range2)
        for i in range(self.chi2.size):
            i1 = min(gridres-1,np.int(np.floor((self.res[i,chainindex1] - extent[0])/step1)))
            i2 = min(gridres-1, np.int(np.floor((self.res[i,chainindex2] - extent[2])/step2)))
            avpixels[i1,i2] += self.chi2[i]
            maxpixels[i1,i2] = max(maxpixels[i1,i2], self.chi2[i])
            nbs[i1,i2] += 1
        for i1 in range(gridres):
            for i2 in range(gridres):
                if (nbs[i1,i2]>0):
                    avpixels[i1,i2] /= nbs[i1,i2]
        minpix= avpixels.min()
        for i1 in range(gridres):
            for i2 in range(gridres):
                if (nbs[i1,i2]==0):
                    avpixels[i1,i2] = minpix
                else:
                    nbs[i1,i2] /= self.chi2.size
        # np.where(nbs>0, pixels/nbs, 0)
        # np.where(pixels==0, pixels.min(), pixels)
        return avpixels, maxpixels, nbs, extent #range1, range2

    def Plot_2D_maps(self,p1,p2,gridres=100):
        fig = figure()

        avmap, maxmap, densmap, extent = self._Compute_2D_maps(p1,p2, gridres=gridres)
        pltav = fig.add_subplot(221)
        pltav.imshow(avmap)#, extent=extent)
        fig.colorbar(pltav.get_images()[0])
        pltav.set_xlabel(self.symbols[p1])
        pltav.set_ylabel(self.symbols[p2])
        pltav.set_title(r"Average $-\chi^2/2$")
        # maxmap = = self._Compute_2D_map(p1,p2, gridres=gridres, maptype='chi2max')
        pltmax = fig.add_subplot(222)
        pltmax.imshow(maxmap)#, extent=extent)
        fig.colorbar(pltmax.get_images()[0])
        pltmax.set_title(r"Max $-\chi^2/2$")

        # densmap = self._Compute_2D_mapw(p1,p2, gridres=gridres, maptype='density')
        pltdens = fig.add_subplot(223)
        pltdens.imshow(densmap)#, extent=extent)
        fig.colorbar(pltdens.get_images()[0])
        pltdens.set_title(r"Walker density")

        fig.set_tight_layout(True)

        return fig




    def Unpriorize(self, funcprior):
        '''
        return chi2 - funcprior(priorparameters).sum(), funcprior(priorparameters)
        Returns the chain of chi2 with the effect of the prior given by funcprior removed
        funcprior(priorparameters) returns the ln(prior probability) for each parameter in priorparameters (and in the same order).
        priorparameters are given as an array containing of size mcmc_results.parammap.max() +1. The chain parameters are converted to priorparameters using mcmc_results.parammap.
        '''
        chainnoprior = self.chi2.copy()
        params = np.zeros(self.parammap.max()+1)
        priors = np.zeros((self.chi2.size, params.size))
        for i in range(self.chi2.size):
            for j in range(self._ndim):
                params[self.parammap[j]] = self.res[i,j]
            priors[i] = funcprior(params)
            chainnoprior[i] -= funcprior(params).sum()
        print(params)
        return chainnoprior, priors


    def Parameter_filter(self, param_nb, minvalue, maxvalue ):
        '''
            Return a mcmcresult object with only the steps of the chain for which min value < parameters[param_nb] < maxvalue

            param_nb : a list of absolute parameter index, and mcmcresult.param_map is used to find the corresponding index in the chain.

            minvalue,maxvalue : lists, these values include the scaling and initial value shifts

            If param_nb, minvalue, and maxvalue are a list or an array, then do it for the set of parameters in paramnb
            TODO : scalar case nor implemented
        '''
        paramnb1 = np.array(param_nb, dtype=int)
        paramnb = np.ndarray(paramnb1.shape, dtype=int)
        print(paramnb1)
        print(paramnb)
        maxpa = self.parammap.max()
        for idim in range(paramnb.size): # Invert param_nb using param_map
            i=0
            while (i <= maxpa and self.parammap[i] != paramnb1[idim]):
                i+=1
            if i > maxpa :
                print("Error : parameter number does not exist")
            paramnb[idim] = i

        mins = np.array(minvalue)
        maxs = np.array(maxvalue)

        for idim in range(paramnb.size) :
            im = 0
            chi2i = np.ndarray(self.chi2.shape)
            resi = np.ndarray(self.res.shape)
            for i in range(self.res.shape[0]):
                if (self.res[i, paramnb[idim]] < maxs[idim]) and  (self.res[i, paramnb[idim]] > mins[idim]):
                    resi[im] = self.res[i]
                    chi2i[im] = self.chi2[i]
                    im += 1
            chi2i = chi2i[:im]
            resi = resi[:im]


        for i in range(self._ndim) :
            resi[:, i] -= self.paraminitials[self.parammap[i]]

        return mcmcresult(chains=resi, chi2=chi2i, resfile=None,  parameter_map=self.parammap.copy(),
                          parameter_scales=np.ones(self.parammap.max()+1), parameter_initials=self.paraminitials.copy(), symbols=self.symbols, derived_parameters=self.derive_parameters, print_parameters=self.print_parameters)



    def Print_parameters_latex(self, levelpercentage=95, start=0, end=-1, print_param_map = None):
        '''
            Print the list parameters with error bars
        '''
        pars = self.Mean_parameters(start=start,  end=end)
        errsplus, errsminus = self.Errors(start=start,  end=end, levelpercentage=levelpercentage)
        errsplusderived, errsminusderived = self.Errors_derived_parameters(start=start,  end=end, levelpercentage=levelpercentage)


        if self.print_parameters is not None :
            return self.print_parameters(external_parameters=pars, external_errorsplus=errsplus, external_errorsminus=errsminus, external_param_map=self.parammap, significant_digits=2, print_param_map=print_param_map, derived_errorplus = errsplusderived, derived_errorminus=errsminusderived)
        else :
            print("Not implemenented yet !!!")

        return


    def Convergence_KS(self,windowsize, reference_subsample=1., start=0, end=-1, parameter_list = None,figsize=(10,10)):
        '''
        figure, ks = Convergence_KS(windowsize, start=0, end=-1, parameter_list = None,figsize=(10,10))

        Check the stationarity of the marginalized distribution of each parameter. Each 1-parameter chain is divided in subsamples of size ''windowsize'', which is compared to a reference subsample ( by default the last one) though a Kolmogorov-Smirnof test.

        ** Input
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * reference_subsample : between 0 and 1, position of the subsample taken as reference
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1
            * parameter_list : if None, all parameters of the chain are plotted, otherwise parameter list contains the user indices (inversed to chain indices using self.parammap)
            * figsize : size of the figure

        ** Return
            * figure : Plot of ks
            * ks : ks.shape = (number of parameters, number of subsamples of size windowsize)
            ks[parameter, subsample] = p-value of Kolmogorov-Smirnov test between the subsample and the reference subsample
        '''
        if parameter_list is not None:
            paramlist = self.inversparammap[parameter_list] # convert parameter_list to chain indexes
        else :
            paramlist = np.arange(self._ndim)

        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        nsamples = (maxit -start) // windowsize
        ks = np.zeros((len(paramlist), nsamples))
        labels = [self.symbols[self.parammap[paramlist[i]]] for i in range(paramlist.size) ]
        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)

        for ip in range(paramlist.size):
            testchain = self.res[start:maxit,paramlist[ip]]
            i=0
            refsample = testchain[(nsamples-i-1)*windowsize: (nsamples-i) * windowsize]
            for i in range(nsamples-1):
                    ks_stat, ks[ip,nsamples - i-1] = ks_2samp(refsample,testchain[(nsamples-i-1)*windowsize: (nsamples-i) * windowsize])

            plt.plot(ks[ip],'.-', label=labels[ip])

        plt.set_xlabel('Subsample index')
        plt.set_ylabel('KS p-value')
        plt.legend(loc=0)
        fig.set_tight_layout(True)

        return fig, ks


    def Plot_convergence_std(self,windowsize, reference_subsample=1., start=0, end=-1, parameter_list = None,figsize=(10,10)):
        '''
        figure, ks = Convergence_KS(windowsize, start=0, end=-1, parameter_list = None,figsize=(10,10))

        Check the stationarity of the marginalized distribution of each parameter. Each 1-parameter chain is divided in subsamples of size ''windowsize'', which is compared to a reference subsample ( by default the last one) though a Kolmogorov-Smirnof test.

        ** Input
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * reference_subsample : between 0 and 1, position of the subsample taken as reference
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1
            * parameter_list : if None, all parameters of the chain are plotted, otherwise parameter list contains the user indices (inversed to chain indices using self.parammap)
            * figsize : size of the figure

        ** Return
            * figure : Plot of ks
            * ks : ks.shape = (number of parameters, number of subsamples of size windowsize)
            ks[parameter, subsample] = p-value of Kolmogorov-Smirnov test between the subsample and the reference subsample
        '''
        if parameter_list is not None:
            paramlist = self.inversparammap[parameter_list] # convert parameter_list to chain indexes
        else :
            paramlist = np.arange(self._ndim)

        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        nsamples = (maxit -start) // windowsize
        stds = np.zeros((len(paramlist), nsamples))
        labels = [self.symbols[self.parammap[paramlist[i]]] for i in range(paramlist.size) ]
        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)

        for ip in range(paramlist.size):
            testchain = self.res[start:maxit,paramlist[ip]]
            i=np.int(np.floor(reference_subsample* nsamples))
            refsample = testchain.std() #testchain[max((i-1),0)*windowsize: (max(i-1,0) +1) * windowsize].std()
            for i in range(nsamples-1):
                stds[ip,nsamples-i-1] = testchain[(nsamples-i-1)*windowsize: min((nsamples-i) * windowsize, testchain.shape[0])].std()

            plt.plot((stds[ip]-refsample)/refsample,'.-', label=labels[ip])

        plt.set_xlabel('Subsample index')
        plt.set_ylabel('Standard deviation relative whole sample standard deviation')
        plt.set_ylim(-10., 10.)
        plt.legend(loc=0)
        fig.set_tight_layout(True)

        return fig, stds


    def Plot_convergence_mean(self,windowsize, start=0, end=-1, parameter_list = None,figsize=(10,10)):
        '''
        figure, ks = Convergence_KS(windowsize, start=0, end=-1, parameter_list = None,figsize=(10,10))

        Check the stationarity of the marginalized distribution of each parameter. Each 1-parameter chain is divided in subsamples of size ''windowsize'', which is compared to a reference subsample ( by default the last one) though a Kolmogorov-Smirnof test.

        ** Input
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * reference_subsample : between 0 and 1, position of the subsample taken as reference
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1
            * parameter_list : if None, all parameters of the chain are plotted, otherwise parameter list contains the user indices (inversed to chain indices using self.parammap)
            * figsize : size of the figure

        ** Return
            * figure : Plot of ks
            * ks : ks.shape = (number of parameters, number of subsamples of size windowsize)
            ks[parameter, subsample] = p-value of Kolmogorov-Smirnov test between the subsample and the reference subsample
        '''
        if parameter_list is not None:
            paramlist = self.inversparammap[parameter_list] # convert parameter_list to chain indexes
        else :
            paramlist = np.arange(self._ndim)

        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        nsamples = (maxit -start) // windowsize
        means = np.zeros((len(paramlist), nsamples))
        labels = [self.symbols[self.parammap[paramlist[i]]] for i in range(paramlist.size) ]
        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)

        for ip in range(paramlist.size):
            testchain = self.res[start:maxit,paramlist[ip]]
            refsample = testchain.mean()
            refsample_std = testchain.std()
            if refsample_std == 0.:
                printf("Error ! Std of whole sample = 0  for parameter {:d}".format(paramlist[ip]))
                refsample_std=1.
            for i in range(nsamples-1):
                means[ip,nsamples-i-1] = testchain[(nsamples-i-1)*windowsize: min((nsamples-i) * windowsize, testchain.shape[0])].mean()

            plt.plot((means[ip]-refsample)/refsample_std,'.-', label=labels[ip])

        plt.set_xlabel('Subsample index')
        plt.set_ylabel('Mean variation relative to whole sample standard deviation')
        plt.set_ylim(-10., 10.)
        plt.legend(loc=0)
        fig.set_tight_layout(True)

        return fig, means


    def Plot_convergence_autocorrelation(self,windowsize=1, start=0, end=-1, parameter_list = None,figsize=(10,10)):
        '''
        figure, ks = Convergence_KS(windowsize, start=0, end=-1, parameter_list = None,figsize=(10,10))

        Check the stationarity of the marginalized distribution of each parameter. Each 1-parameter chain is divided in subsamples of size ''windowsize'', which is compared to a reference subsample ( by default the last one) though a Kolmogorov-Smirnof test.

        ** Input
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * reference_subsample : between 0 and 1, position of the subsample taken as reference
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1
            * parameter_list : if None, all parameters of the chain are plotted, otherwise parameter list contains the user indices (inversed to chain indices using self.parammap)
            * figsize : size of the figure

        ** Return
            * figure : Plot of ks
            * ks : ks.shape = (number of parameters, number of subsamples of size windowsize)
            ks[parameter, subsample] = p-value of Kolmogorov-Smirnov test between the subsample and the reference subsample
        '''
        if parameter_list is not None:
            paramlist = self.inversparammap[parameter_list] # convert parameter_list to chain indexes
        else :
            paramlist = np.arange(self._ndim)

        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        nsamples = (maxit -start) // windowsize
        autocor = np.zeros((len(paramlist), nsamples))
        labels = [self.symbols[self.parammap[paramlist[i]]] for i in range(paramlist.size) ]
        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)

        for ip in range(paramlist.size):
            testchain = np.zeros(nsamples)
            for i in range(nsamples):
                testchain[i] = self.res[start  + i*windowsize : start + (i+1)*windowsize, paramlist[ip]].mean()
                #testchain -= testchain.mean()
            autocor[ip][:testchain.size//2] = myautocor(testchain, testchain)#np.correlate(testchain, testchain, 'same') / np.correlate(np.ones(nsamples), np.ones(nsamples))
            #autocor /= autocor.max()

            plt.plot(autocor[ip,autocor.argmax():], '.-', label=labels[ip])

        plt.set_xlabel('Subsample index lag')
        plt.set_ylabel('Autocorrelation')
        plt.legend(loc=0)
        fig.set_tight_layout(True)

        return fig, autocor




    def Convergence_stds(self, param_nb, windowsize, start=0, end=-1):
        '''
        stds = Convergence_std(windowsize, start=0, end=-1)

        Check the stationarity of the marginalized distribution of each parameter, by computing the windowed standard deviation of the chain.
        ** Input
            * param_nb : parameter index in the reference parameter space (the chain-space parameter index is self.inversparammap[param_nb])
                         If param_nb <= 0 : take the self.chi2 chain.
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1

        ** Return
            * stds : stds.shape = (number of subsamples of size windowsize)
            stds[subsample] = relative difference betzeen the standard deviation and the reference standard deviation for the subsample number "subsample"
        '''

        if param_nb >= 0 :
            paramnb = self.inversparammap[param_nb]
        else :
            paramnb = param_nb

        if end < 0 :
            maxit = self.res.shape[0]
        else :
            maxit = end

        nsamples = (maxit -start) // windowsize
        stds = np.zeros(nsamples-1)

        if paramnb >= 0 :
            testchain = self.res[start:maxit,paramnb]
        else :
            testchain = self.chi2[start:maxit]
        refsample = testchain.std() #testchain[max((i-1),0)*windowsize: (max(i-1,0) +1) * windowsize].std()
        for i in range(nsamples-1):
            stds[i] = testchain[i*windowsize: (i+1) * windowsize].std()
            stds[i] = (stds[i]-refsample)/refsample

        return stds


    def Convergence_means(self, param_nb, windowsize, start=0, end=-1):
            '''
            means = Convergence_means(param_nb, windowsize, start=0, end=-1)

            Check the stationarity of the marginalized distribution of each parameter param_nb, by Computing the windowed mean of the chain.
            ** Input
                * param_nb : parameter index in the reference parameter space (the chain-space parameter index is self.inversparammap[param_nb])
                             If param_nb <= 0 : take the self.chi2 chain.
                * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
                * start : index of the first chain element to take into account
                * end : index of the last chain element to take into account plus 1

            ** Return
                * means : means.shape = (number of subsamples of size windowsize)
                means[subsample] = relative difference betzeen the standard deviation and the reference standard deviation for the subsample number "subsample"
            '''

            if param_nb >= 0 :
                paramnb = self.inversparammap[param_nb]
            else :
                paramnb = param_nb

            if end < 0 :
                maxit = self.res.shape[0]
            else :
                maxit = end

            nsamples = (maxit -start) // windowsize
            means = np.zeros(nsamples-1)

            if paramnb >= 0 :
                testchain = self.res[start:maxit,paramnb]
            else :
                testchain = self.chi2[start:maxit]
            scalesample = testchain.std() #testchain[max((i-1),0)*windowsize: (max(i-1,0) +1) * windowsize].std()
            refsample = testchain.mean()
            for i in range(nsamples-1):
                means[i] = testchain[i*windowsize: (i+1) * windowsize].mean()
                means[i] = (means[i]-refsample)/scalesample

            return means

    def Convergence_mean_autocorrelations(self, param_nb, windowsize, start=0, end=-1):
        '''
        autocors = Convergence_autocors(param_nb, windowsize, start=0, end=-1)

        Check the stationarity of the marginalized distribution of each parameter param_nb, by computing the autocorrelation function of the windowed mean of the chain.


        ** Input
            * param_nb : parameter index in the reference parameter space (the chain-space parameter index is self.inversparammap[param_nb])
                         If param_nb <= 0 : take the self.chi2 chain.
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1

        ** Return
            * autocors : autocors.shape = (number of subsamples of size windowsize // 2)
        '''
        means = self.Convergence_means(param_nb, windowsize, start=start, end=end)

        nsamples = means.size
        autocors = np.zeros(nsamples//2)

        autocors = myautocor(means, means)

        return autocors


    def Convergence_std_autocorrelations(self, param_nb, windowsize, start=0, end=-1):
        '''
        stdautocors = Convergence_std_autocorrelations(param_nb, windowsize, start=0, end=-1)

        Check the stationarity of the marginalized distribution of each parameter param_nb, by computing the autocorrelation function of the windowed stamdard deviation of the chain.


        ** Input
            * param_nb : parameter index in the reference parameter space (the chain-space parameter index is self.inversparammap[param_nb])
                         If param_nb <= 0 : take the self.chi2 chain.
            * windowsize : lenght of the subsample. One may take at least the number of walkers in an affine invariant MCMC
            * start : index of the first chain element to take into account
            * end : index of the last chain element to take into account plus 1

        ** Return
            * autocors : autocors.shape = (number of subsamples of size windowsize // 2)
        '''
        stds = self.Convergence_stds(param_nb, windowsize, start=start, end=end)

        nsamples = stds.size
        autocors = np.zeros(nsamples//2)

        autocors = myautocor(stds, stds)

        return autocors


    def Plot_convergence(self, windowsize, start=0, end=-1, parameter_list=None, figsize=(18,14)):
        '''
            Plot for each parameter in parameter_list a figure with 3 panels : the mean and the standard deviation in unit of the whole-chain standard deviation, and the autocorrelation function of the mean.
            Means and standard deviations are shown as the variations with respect to the whole-sample quantities.

            ** Input
                * windowsize : window length to use to calculate the local mean and standard deviations
                * start : index of the first chain element to take into account
                * end : index of the last chain element to take into account plus 1
                * parameter_list : list of parameters to plot with indexes given in the reference space (the chain-space parameter index is self.inversparammap[param_list[i]]).
                                    If "None", plot all parameters.
                                    If parameter_list[i] < 0 : take the self.chi2 chain.
                * figsize : size of each figure given as a tuple (width, height)
        '''
        if parameter_list is None :
            paramlist = np.zeros(self.parammap.size + 1, dtype = int)
            paramlist[0] = -1
            paramlist[1:] = self.parammap.copy()
            paramlist[:].sort()
        else :
            paramlist = []
            for i in range(len(parameter_list)):
                if parameter_list[i] >= 0:
                    paramlist.append(self.inversparammap[parameter_list[i]])
                else :
                    paramlist.append(parameter_list[i])

        for i in range(len(paramlist)):
            fig = figure(figsize=figsize)

            means = self.Convergence_means(paramlist[i], windowsize, start=start, end=end)
            stds = self.Convergence_stds(paramlist[i], windowsize, start=start, end=end)
            autocors = self.Convergence_mean_autocorrelations(paramlist[i], windowsize, start=start, end=end)
            stdautocors = self.Convergence_std_autocorrelations(paramlist[i], windowsize, start=start, end=end)

            if paramlist[i] >= 0:
                symbol = self.symbols[paramlist[i]]
            else :
                symbol = "$\chi^2$"

            pltm = fig.add_subplot(221)
            plts = fig.add_subplot(223)
            pltsa = fig.add_subplot(224)
            plta = fig.add_subplot(222)

            pltm.axhline(0., color="red", alpha=0.5)
            plts.axhline(0., color="red", alpha=0.5)
            plta.axhline(0., color="red", alpha=0.5)

            pltm.set_xlabel('Window index')
            pltm.set_ylabel('Win mean / std ')
            pltm.legend(loc=0)
            pltm.plot(means,label=symbol)
            pltm.legend(loc=0)

            plts.set_xlabel('Window index')
            plts.set_ylabel('Win std / std ')
            plts.legend(loc=0)
            plts.plot(stds,label=symbol)
            plts.legend(loc=0)

            plta.set_xlabel('Window index')
            plta.set_ylabel('Win mean autocor / std ')
            plta.legend(loc=0)
            plta.plot(autocors[:autocors.size //2],label=symbol)
            ax2 = plta.twinx()
            ax2.set_ylim(plta.get_ylim())
            ax2.set_yticks(plta.get_yticks())
            plta.legend(loc=0)

            pltsa.set_xlabel('Window index')
            pltsa.set_ylabel('Win std autocor / std ')
            pltsa.legend(loc=0)
            pltsa.plot(stdautocors[:stdautocors.size //2],label=symbol)
            ax2 = plta.twinx()
            ax2.set_ylim(plta.get_ylim())
            ax2.set_yticks(plta.get_yticks())
            pltsa.legend(loc=0)

            fig.set_tight_layout(True)
            fig.show()

        return



class mcmcdiagnostics:
    def __init__(self, diagnosticslog=None, diagnostics_file=None, first=None, last = None, filename_suffix= '.txt'):
        '''
            Analyse diagnostics logs  as produced by MCMC_pai_pt.exe.
            Logs are arrays of dimension (nlogs, ndiagnostics_per_log).
            Each column contains:
                 0                  1                   2                       3                       4                           5               6                           7                           8                               9
            Iteration nb  |   Temperature index  |   Half-iteration   |   Current Walker index |  Mover Walker  index |   Step factor z   |   Current lnposterior   |   Candidate lnposterior   |   Acceptance probability   |   Acceptance outcome

            - diagnosticslog : array of diagnostics
            - diagnostics_file : if no  diagnosticslog provided, load from this filename
            - first, last, filename_suffix: if provided, merge diagnostics logs from all files named
                            diagnostics_file + i + filename_suffix, where first <= i <= last

            Mover = walker randomly picked out to perform the next move such
            that the candidate is determined by :
            candidate = mover + z * (current - mover )
            and z is a random number.
        '''
        if diagnosticslog is not None:
            self.diag = diagnosticslog.copy()
            return
        if first is None or last is None :
            self.diag = np.loadtxt(diagnostics_file)
        else:
            nfiles = last - first +1
            print('Loading {:d} files from index {:d} to {:d} included.'.format(nfiles, first, last))
            interdiag = np.loadtxt(diagnostics_file+"{:d}".format(first)+filename_suffix)
            ndiag = interdiag.shape[0]
            ndiagcol = interdiag.shape[1]
            self.diag=np.zeros((nfiles * ndiag, ndiagcol))
            self.diag[:ndiag, :] = interdiag
            for i in range(first+1,last+1):
                interdiag = np.loadtxt(diagnostics_file+"{:d}".format(i)+filename_suffix)
                if interdiag.shape != (ndiag, ndiagcol) :
                    print("Warning ! File {:d} has diffferent dimensions : ({:d},{:d}) instead of ({:d},{:d})".format(i, interdiag.shape[0], interdiag.shape[1], ndiag, ndiagcol))
                self.diag[i* ndiag :ndiag * (i+1), :] = interdiag
        return


    def Plot_mover_histogram(self, bins=20):
        '''
            Mover = walker randomly picked out to perform the next move such
            that the candidate is determined by :
            candidate = mover + z * (current - mover )
            and z is a random number.
        '''
        fig = figure()
        plt = fig.add_subplot(111)
        plt.hist(self.diag[:,4], bins = bins)
        plt.set_xlabel('Mover index')
        return fig

    def Plot_mover_vs_current(self, alpha=0.1):
        '''
            Mover = walker randomly picked out to perform the next move such
            that the candidate is determined by :
            candidate = mover + z * (current - mover )
            and z is a random number.
        '''
        fig = figure()
        plt = fig.add_subplot(111)
        plt.scatter(self.diag[:,3], self.diag[:,4], alpha=alpha)
        plt.set_xlabel('Current index')
        plt.set_ylabel('Mover index')
        return fig

    def Plot_proposition_distribution(self, bins = 20, aparam =2.):
        fig = figure()
        plt = fig.add_subplot(111)
        vals, xbins, patches = plt.hist(self.diag[:,5], bins = bins)
        distroth = aparam / np.sqrt(xbins[1:])
        distroth *=  vals[0] / distroth[0]
        plt.plot(xbins[1:], distroth)
        plt.set_xlabel('Current index')
        plt.set_ylabel('Mover index')
        return fig

    def Plot_current_lnposterior_histogram(self, bins=50):
        '''
        This plots the histogram of current lnposteriors.
        '''
        fig = figure()
        plt = fig.add_subplot(111)
        plt.hist(self.diag[:,6], bins = bins)
        plt.set_xlabel('Current lnposterior')
        return fig

    def Plot_candidate_lnposterior_histogram(self, bins=50):
        '''
        This plots the histogram of candidate lnposteriors.
        '''
        fig = figure()
        plt = fig.add_subplot(111)
        plt.hist(self.diag[:,7], bins = bins)
        plt.set_xlabel('Candidate lnposterior')
        return fig


    def Plot_candidate_vs_current_lnposterior(self, alpha=0.02):
        fig = figure()
        plt = fig.add_subplot(111)
        plt.scatter(self.diag[:,6], self.diag[:,7], alpha=alpha)
        plt.set_xlabel('Current lnposterior')
        plt.set_ylabel('Candidate lnposterior')
        return fig

    def Plot_acceptance_probability_vs_dlnposterior(self, alpha=0.02):
        fig = figure()
        plt = fig.add_subplot(111)
        acceptproba = self.diag[:,8]
        for i in range(len(acceptproba)):
            acceptproba[i] = min(1., acceptproba[i])
        plt.scatter(self.diag[:,7] - self.diag[:,6], acceptproba, alpha=alpha)
        plt.set_xlabel('Candidate - Current lnposterior')
        plt.set_ylabel('Acceptance probability')
        return fig

    def Plot_acceptance_probabilities(self, bins=50):
        fig = figure()
        pltacc = fig.add_subplot(131)
        pltrej = fig.add_subplot(132)
        plttot = fig.add_subplot(133)
        acc_accproba = []
        rej_accproba = []
        for i in range(self.diag.shape[0]):
            if (self.diag[i,9] == 1.):
                acc_accproba.append(min(1.,self.diag[i,8]))
            else:
                rej_accproba.append(self.diag[i,8])
        pltacc.hist(acc_accproba, bins = bins)
        pltrej.hist(rej_accproba, bins = bins)
        plttot.hist(self.diag[i,8], bins = bins)
        pltacc.set_xlabel('Accepted acceptance probability')
        pltrej.set_xlabel('Rejected acceptance probability')
        plttot.set_xlabel('All acceptance probability')
        return fig




def myautocor(var1, var2):
    '''
     return autocor

     autocor.size = var1.size // 2

     Compute autocorrelation function of var1 - var1.mean() and var2 - var2.mean() normalized to 1
     var1 and var2 are assumed to be the same length

     NOTE : numpy.correlate does weird things
    '''
    autocor = np.zeros(var1.size)
    mean1 = var1.mean()
    mean2 = var2.mean()
    for ilag in range(var1.size//2):
        for i in range(var2.size - ilag):
             autocor[ilag] += (var1[i]- mean1)* (var2[i+ilag]-mean2)
        autocor[ilag]/= (var1.size - ilag)
    return (autocor[:var1.size//2]/autocor.max()).copy()


def Remapchain(filename, newfilename, parammap, newparammap):
    '''
        return chain_with_newparammap

        Load a chain from filename to which is associated a parameter map parammap : i_chain -> i_user.
        Remap it to a chain with map newparammap and save it to newfilename

        newparammap must me included in parammap
    '''
    dat = np.loadtxt(filename)
    ndim = dat.shape[1] - 1
    nit = dat.shape[0]
    newdat = np.zeros((nit, len(newparammap) +1))

    intermap = np.zeros(len(newparammap), dtype=int)
    for i in range(len(newparammap)):
        j = 0
        while (parammap[j] != newparammap[i]) :
            j +=1
        intermap[i] = j

    for i in range(nit):
        for j in range(len(newparammap)):
            newdat[i,j] = dat[i,intermap[j]]
        newdat[i,len(newparammap)] = dat[i, ndim]

    np.savetxt(newfilename, newdat)

    print("You need to change the header by hand !")

    return newdat
