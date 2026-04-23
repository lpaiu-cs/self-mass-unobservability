# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#coding:utf8
#/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 #*
 #* Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 #*/
from python_Fittriple_interface import PyFittriple
import numpy as np
from scipy import interpolate
from scipy.signal import lombscargle
from pylab import *
from Constantes_physiques import rad_to_hms, rad_to_dms, arcsec
from outils import Print_numbers_with_errors_latex, Print_numbers
from astropy.modeling import models,fitting
from astropy import units as apu
import astropy.time as apt
import astropy.coordinates as apco
from astropy.timeseries import LombScargle

ioff() # stop matplotlib.pyplot interactive mode
#from outils import Print_numbers

####
####
# To convert a tempo parfile to a parfile suited for this program see Donnees/mypar_from_tempopar.py
####
####

clight = 299792458. #    ! m/s Par déf
daysec = 86400.
Ggrav = 6.67408e-11      # (31) m^3 s^-2 kg^-1 (CODATA 2014)
Msol = 1.988435e30         # Masse du Soleil (kg) (CODATA 2014 / IAU)
GMsol = Ggrav * Msol



def Print_value_with_error_bar(value,error) :
    '''
        Prints in Latex format value as x.xxx(yy)\cdot 10^(z) where yy are two digits
        corresponding to error and x and z are the base and exponent of value.
    '''
    if value == 0.:
        expo = 0
    else :
        expo = np.int(np.floor(np.log10(np.abs(value))))
    nf = np.int(np.floor(np.log10(error)))
    base = np.sign(value) * value/ 10**expo
    prec = max(-(nf - expo) -1, 0)
    chne = "${:.{prec}f}".format(base, prec= prec)
    if prec > 0 :
        chne+="({:d})".format(np.int(error *10 /10**nf))
    if expo == 0:
        chne += r"$"
    else :
        chne+= r"\cdot 10^{" + "{:d}".format(expo) + r"}$"
    return chne



def covariance_minuit(covar_file, nbdigits=2):
    cov = np.loadtxt(covar_file)
    covnorm = cov.copy()
    sig = np.ndarray(cov.shape[0])
    for i in range(sig.size) :
       sig[i] = cov[i,i]

    chaine = ''
    for i in range(sig.size):
        for j in range(sig.size):
            covnorm[i,j] /= np.sqrt(sig[i]*sig[j])
            chaine += '{cov:.{nbdigits:d}}   '.format(cov = covnorm[i,j], nbdigits = nbdigits)
        chaine += '\n'

    print(chaine)

    return cov, covnorm, sig




class Fittriple(PyFittriple):
    
    fontsizes={'xx-large':'xx-large', 'x-large':'x-large', 'large':'large'} # for plots 
    
    
    def Print_parameters_latex(self, external_parameters= None, external_errorsplus= None, external_errorsminus=None, print_param_map=None, external_param_map=None, significant_digits=2, derived_errorplus=None, derived_errorminus=None ):
        '''
            Print the list of fitted, fittable and deduced parameters with error bars as a Latex table
            Can be used as an independent routine by specifying the parameters through ''external_parameters''
            external_parameters : array of parameters
            external_param_map : map iext in external_parameters -> external_parameter_map[iext] = iin in internal representation ; If not provided, assumed to be identity.
            external_errorsplus : error above mean on each parameters transformed through external_parameter_map. If not provided, the internal scale value is used.
            external_errorsminus : error above mean on each parameters transformed through external_parameter_map. If not provided, assumed to be equal to external_errorsplus.
            print_param_map : map iprint in printing order -> external_parameter_map[iprint] = iin in internal representation
            significant_digits : Number of digits to print for the error
        '''
        parammap1 = np.array([0,1,
                              3,5,2,4,6,7,
                              9, 11, 8, 10, 12, 13,
                              17, 18,
                              14,15,16,
                              21,22,23,24,25,26,27,28,
                              20])

        print("\n WARNING : ONLY WORKS IF Q IS IN THE PARFILE (EVEN IF NOT FITTED) BUT NO DMX OR FD \n\n")

        if self.Get_parameter_set()==1 :
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
                 "Mass parameter 1",
                 "Mass parameter 2",
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

            symbols = [r"$\bar{f}$",                     #0
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
                   r"$\frac{m_p+m_i}{2}$",      #14
                   r"$\frac{m_p - m_i}{2}$",
                   r"$m_o$",
                   r"$\Omega_o$",                #17
                   r"$\delta\Omega$",           #18
                   r"$\delta\phi_0$",
                   r"$\Delta$",
                   r"$\alpha$",                 #21
                   r"$\delta$",
                   r"$d$",
                   r"$\mu_\alpha$",             #24
                   r"$\mu_\delta$",
                   r"$\mu_d$",
                   r"$\mathrm{DM}$",            #27
                   r"$\mathrm{DM}'$"
                   ]

            unitstr= [r"\si{s^{-1}}", r"\si{s^{-2}}",
                  r"",r"\mathrm{ls}",r"\si{}",r"\mathrm{ls}",r"\mathrm{days}",r"\mathrm{days}",
                  r"",r"\mathrm{ls}",r"\si{}",r"\mathrm{ls}",r"\mathrm{days}",r"\mathrm{days}",
                  r"\Msol",r"\Msol",r"\Msol",r"^{\circ}", r"^{\circ}",r"", r"",r"",r"",r"\mathrm{ly}",
                  r"\mathrm{mas}/\mathrm{yr}",r"\mathrm{mas}/\mathrm{yr}", r"\mathrm{mas}/\mathrm{yr}",
                  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}",  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}\mathrm{yr}^{-1}"
                  ]

        else:
            names = ["Rescaled spin frequency",
                 "Rescaled spin frequency derivative",
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
                 r"SEP $\gamma$",
                 r"SEP $\beta_{pp}^0$",
                 r"SEP $\beta_{00}^p$",
                 r"SEP $\beta_{0p}^0$",
                 "Right ascension",
                 "Declination",
                 "Distance",
                 "Right-ascension proper motion",
                 "Declination  proper motion",
                 "Distance  proper motion",
                 "Dispersion measure",
                 "Dispersion measure variation",
                 #"Dispersion measure piecewise",
                 #"Dispersion measure frequency",
                 "Inner WD quadrupole moment",
                 "TOA uncertainty rescaling",
                 "Difference in inclinations"
                 ]

            symbols = [r"$f$",                     #0
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
                   r"$m_p$",                    #14
                   r"$m_i$",
                   r"$m_o$",
                   r"$\Omega_O$",                #17
                   r"$\delta\Omega = \Omega_I - \Omega_O$",           #18
                   r"$\delta\phi_0$",
                   r"$\Delta$",                 #20
                   r"$\gamma$",
                   r"$\beta_{pp}^0$",
                   r"$\beta_{00}^p$",
                   r"$\beta_{0p}^0$",
                   r"$\alpha$",                 #25
                   r"$\delta$",
                   r"$d$",
                   r"$\mu_\alpha$",             #28
                   r"$\mu_\delta$",
                   r"$\mu_d$",
                   r"$\mathrm{DM}$",            #31
                   r"$\mathrm{DM}'$",
                   #r"$\mathrm{DMX}'$",
                   #r"$\mathrm{FD}'$",
                   r"$\mathrm{Q}$",
                   r"$\mathrm{EFAC}$",          # 36
                   r"$\delta i = i_I - i_O$",
                   ]

            unitstr= [r"\si{s^{-1}}", r"\si{s^{-2}}",
                  r"",r"\mathrm{ls}",r"\si{}",r"\mathrm{ls}",r"\mathrm{days}",r"\mathrm{days}",
                  r"",r"\mathrm{ls}",r"\si{}",r"\mathrm{ls}",r"\mathrm{days}",r"\mathrm{days}",
                  r"\Msol",r"\Msol",r"\Msol",r"^{\circ}", r"^{\circ}",r"", r"", r"", r"", r"", r"", # beta_0p^0
                  r"",r"",r"\mathrm{kpc}", # alpha, delta, d
                  r"\mathrm{mas}/\mathrm{yr}",r"\mathrm{mas}/\mathrm{yr}", r"\mathrm{mas}/\mathrm{yr}",
                  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}",  r"\mathrm{pc}\cdot\mathrm{cm}^{-3}\mathrm{yr}^{-1}",
                  #r"\mathrm{pc}\cdot\mathrm{cm}^{-3}", r"\mathrm{pc}\cdot\mathrm{cm}^{-3}/\mathrm{Hz}",
                  r"check unit", r"", r"^{\circ}"
                  ]

        # For derived parameters
        names_derived = ['Pulsar semimajor axis',
                         'Inner-orbit inclination',
                         'Inner-orbit eccentricity',
                         'Inner-orbit longitude of periastron',
                         'Inner-orbit long. of asc. node',
                         'Inner orbit time of periastron passage',
                         'Inner-orbit period',
                         'Inner-binary semimajor axis',
                         'Outer-orbit inclination',
                         'Outer-orbit eccentricity',
                         'Outer-orbit longitude of periastron',
                         'Outer orbit time of periastron passage',
                         'Outer-orbit period',
                         'Pulsar mass',
                         'Inner-companion mass',
                         'Outer-companion mass',
                         'Parallel proper motion',
                         'Plane-of-sky proper motion',
                         'Spin frequency',
                         'Spin frequency derivative']
        symbols_derived = [r'$a_p$',
                         r'$i_I$',
                         r'$e_I$',
                         r'$\omega_I$',
                         r'$\Omega_I$',
                         r'${t_p}_I$',
                         r'$P_I$',
                         r'$a_b$',
                         r'$i_O$',
                         r'$e_O$',
                         r'$\omega_O$',
                         r'${t_p}_O$',
                         r'$P_O$',
                         r'$m_p$',
                         r'$m_i$',
                         r'$m_o$',
                         r'$V_{\parallel}$',
                         r'$V_{\bot}$',
                         r'$f$',
                         r"$f'$"]
        unitstr_derived = [r"\mathrm{ls}",
                           r"^{\circ}",
                           r"",
                           r"^{\circ}",
                           r"^{\circ}",
                           r"\mathrm{MJD}",
                           r"\mathrm{days}",
                           r"\mathrm{ls}",
                           r"^{\circ}",
                           r"",
                           r"^{\circ}",
                           r"\mathrm{MJD}",
                           r"\mathrm{days}",
                           r"\Msol",
                           r"\Msol",
                           r"\Msol",
                           r"\mathrm{km}/\mathrm{s}",
                           r"\mathrm{km}/\mathrm{s}",
                           r"\si{s^{-1}}",
                           r"\si{s^{-2}}"]

        if self.Get_parameter_set() == 6:
            symbols[14] = r"$m_i/m_p$"
            names[14] = "Inner mass ratio"
        else:
            print("Beware this hasnt been thoroughly tested with parameter_set != 6")

        if external_param_map is not None:
            external_parammap = external_param_map.copy()
        else :
            external_parammap = None

        pars = self.Get_parameters()
        if external_parameters is not None:
            if external_parammap is None :
                external_parammap = np.arange(external_parameters.size, dtype=int)
            for i in range(external_parameters.size):
                pars[external_parammap[i]] = external_parameters[i]

        errplus = self.Get_parameter_scales()
        errminus = errplus.copy()

        if external_errorsplus is not None:
            if external_parammap is None : external_parammap = np.arange(external_parameters.size, dtype=int)
            for i in range(external_errorsplus.size):
                errplus[external_parammap[i]] = external_errorsplus[i]

        if external_errorsminus is not None :
            for i in range(external_errorsplus.size):
                errminus[external_parammap[i]] = external_errorsminus[i]
        else :
            errminus = errplus.copy()

        if derived_errorplus is not None :
            errplusderived = derived_errorplus.copy()
        else :
            errplusderived = np.zeros(14)

        if derived_errorminus is not None :
            errminusderived = derived_errorminus.copy()
        else :
            errminusderived = errplusderived.copy()



        fittedlist = self.Get_fitted_parameters_list()
        fitted = np.zeros(pars.size)
        for i in range(fittedlist.size):
            fitted[fittedlist[i]] = 1

        if print_param_map is None :
            parammap = np.arange(pars.size)
        else :
            parammap = print_param_map

# Derived parameters need to be calculated before units are changed
        paramsderived = self.Derive_parameters(pars)

# Changing the units
        unitcoeffs = np.ones(pars.size)
        unitcoeffs[0] /= daysec
        unitcoeffs[1] /= daysec**2
        unitcoeffs[3] /= clight
        unitcoeffs[5] /= clight
        unitcoeffs[9] /= clight
        unitcoeffs[11] /= clight
        #unitcoeffs[17] /= pi/180.
        unitcoeffs[27] *= 0.3066013937879527/1000 # lt-y / kpc

        pars = pars * unitcoeffs
        errminus *= unitcoeffs
        errplus *= unitcoeffs
        print("\n pars 17 \n", pars[17])
        pars[17] = np.sign(pars[17]) * (np.abs(pars[17])%(360)) # Omega_o modulo 360deg
        print("\n pars 17 \n", pars[17])

        pars[6] += self.Get_timeshift() # put back to MJD the times of ascending node
        pars[12] += self.Get_timeshift()

        paramsderived[0] /= clight
        paramsderived[7] /= clight
        paramsderived[5] += self.Get_timeshift() # put back to MJD the times of periastron passage
        paramsderived[11] += self.Get_timeshift()
        errminusderived[0] /= clight # a_p in ls
        errminusderived[7] /= clight # a_b in ls
        errplusderived[0] /= clight # a_p in ls
        errplusderived[7] /= clight # a_b in ls


# Writting of the numbers for fitted and fixed quantities
        values = np.ndarray(pars.size, dtype='U100')
        for i in range(pars.size):
            if fitted[i] == 1 :
                values[i] = "$" + Print_numbers_with_errors_latex(pars[i],
                    errplus[i], errorminus=errminus[i], chiffre_significatif_error=significant_digits,
                    not_scientific=3,epsint = 10.**(-14)) + "$"
            else:
                values[i] = "$" + Print_numbers(pars[i], limite_notasci=3, chiffre_significatif = 3)+"$"
        h,m,s = rad_to_hms(pars[25])
        s_errminus = rad_to_hms(errminus[25])[2]
        s_errplus = rad_to_hms(errplus[25])[2]
        #values[21] = "${:2.0f}h{:2.0f}m{:2.4f}s$".format(h,m,s)
        values[25] = '${:2.0f}'.format(h) + r'^\mathrm{h}' + '{:2.0f}'.format(m) + r'^\mathrm{m}'
        values[25] += Print_numbers_with_errors_latex(s,
                    s_errplus, errorminus=s_errminus, chiffre_significatif_error=significant_digits,
                    not_scientific =3,epsint = 10.**(-14) ) + r'^\mathrm{s}$'
        d,m,s = rad_to_dms(pars[26])
        s_errminus = rad_to_dms(errminus[26])[2]
        s_errplus = rad_to_dms(errplus[26])[2]
        #values[22] = "${:2.0f}".format(d) + r"^{\circ} " + "{:2.0f}'{:2.4f}''$".format(m,s)
        values[26] = "${:2.0f}".format(d) + r"^\circ" + "{:2.0f}'".format(m)
        values[26] += Print_numbers_with_errors_latex(s,
                    s_errplus, errorminus=s_errminus, chiffre_significatif_error=significant_digits,
                    not_scientific =3,epsint = 10.**(-14) ) + r"{}'{}'$"

# Writting of table lines for fixed and fitted quantities
        strfitted = ""
        strfixed = ""
        strderived=""

        strfixed += 'Reference epoch' + ' & ' + r'$T_{\mathrm{ref}}$' + ' & ' + '{:.0f}'.format(self.Get_treference()) + r" $" + r'\mathrm{MJD}' +  r"$ \\" + "\n"
        strfixed += 'Position epoch' + ' & ' + r'$T_{\mathrm{pos}}$' + ' & ' + '{:.0f}'.format(self.Get_posepoch()) + r" $" + r'\mathrm{MJD}' +  r"$ \\" + "\n"

        for i in parammap :
            if fitted[i] == 1 :
                if unitstr[i] == "":
                    ustr = ""
                else:
                    ustr =   r" $" + unitstr[i] +  r"$"
                strfitted += names[i] + " & " + symbols[i] + " & " + values[i] + ustr +  r" \\" + "\n"
            else :
                if self.Get_parameter_set() ==1 :
                    if i != 7 and i != 13 : # for periods
                        if unitstr[i] == "":
                            ustr = ""
                        else:
                            ustr =   r" $" + unitstr[i] +  r"$"
                        strfixed += names[i] + " & " + symbols[i] + " & " + values[i] + ustr +  r" \\" + "\n"

        for i in range(pars.size):
            if fitted[i] == 1 :
                values[i] = "$" + Print_numbers_with_errors_latex(pars[i],
                    errplus[i], errorminus=errminus[i], chiffre_significatif_error=significant_digits,
                    not_scientific=3,epsint = 10.**(-14)) + "$"

# Writting of derived parameters:

        for i in range(paramsderived.size):
            value = "$" + Print_numbers_with_errors_latex(paramsderived[i],
                    errplusderived[i], errorminus=errminusderived[i], chiffre_significatif_error=significant_digits,
                    not_scientific=3,epsint = 10.**(-14)) + "$"
            if unitstr_derived[i] == "":
                ustr = ""
            else:
                ustr =   r" $" + unitstr_derived[i] +  r"$"
            strderived += names_derived[i] + " & " + symbols_derived[i] + " & " + value + ustr +  r" \\" + "\n"

# Do the general writting

        str = r"""
\usepackage{gensymb} % \degree
\usepackage{siunitx} % units

\begin{tabular}{lcl}
Parameter & Symbol & Value \\
\hline
\multicolumn{3}{c}{Fixed values} \\
\hline
                """

        str = str + strfixed
        str += r"""\hline
                \multicolumn{3}{c}{Fitted values} \\
                \hline
                """
        str += strfitted
        str +=  r"""\hline
                \multicolumn{3}{c}{Derived values} \\
                \hline
                """
        str += strderived
        str += r"""\end{tabular}"""

        print(str)

        return



    def Plot_all_conserved_quantities(self):
        tinterp = self.Get_tinterp()
        fig = figure(figsize=(10,12))
        plte = fig.add_subplot(311)
        plti = fig.add_subplot(312)
        pltx = fig.add_subplot(313)

        x = self.Get_energies()
        x = ( x - x[0] ) / x[0]
        plte.plot(tinterp, x, label = '(E - Eini) / Eini')
        #plte.set_xlabel('Times of interpolation (MJD)')
        plte.set_ylabel('Error ')
        plte.legend(fontsize=self.fontsizes['large'], loc = 1)
        plte.xaxis.label.set_fontsize(self.fontsizes['large'])
        plte.yaxis.label.set_fontsize(self.fontsizes['large'])
        plte.tick_params(axis='both', which='major', labelsize=self.fontsizes['x-large'])
        plte.tick_params(axis='both', which='minor', labelsize=self.fontsizes['x-large'])
        plte.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        x = self.Get_center_of_mass_impulsions()
        ms = self.Get_masses()
        x  = x / (ms.sum() * clight)
        plti.plot(tinterp, x[:,0], label = 'Px/Mc')
        plti.plot(tinterp, x[:,1], label = 'Py/Mc (to Earth)')
        plti.plot(tinterp, x[:,2], label = 'Pz/Mc')
        #plti.set_xlabel('Times of interpolation (MJD)')
        plti.set_ylabel('Momentum/Mc')
        plti.legend(fontsize=self.fontsizes['large'], loc = 2)
        plti.xaxis.label.set_fontsize(self.fontsizes['large'])
        plti.yaxis.label.set_fontsize(self.fontsizes['large'])
        plti.tick_params(axis='both', which='major', labelsize=self.fontsizes['x-large'])
        plti.tick_params(axis='both', which='minor', labelsize=self.fontsizes['x-large'])
        plti.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        x = self.Get_center_of_mass_positions()
        x = (x - x[0])/clight * 10.**6
        pltx.plot(tinterp, x[:,0], label = 'x')
        pltx.plot(tinterp, x[:,1], label = 'y (to Earth)')
        pltx.plot(tinterp, x[:,2], label = 'z')
        pltx.set_xlabel('Times of interpolation (in days from reference day)')
        pltx.set_ylabel('Position (light-$\mu$s)')
        pltx.legend(fontsize=self.fontsizes['large'], loc =2)
        pltx.xaxis.label.set_fontsize(self.fontsizes['large'])
        pltx.yaxis.label.set_fontsize(self.fontsizes['large'])
        pltx.tick_params(axis='both', which='major', labelsize=self.fontsizes['x-large'])
        pltx.tick_params(axis='both', which='minor', labelsize=self.fontsizes['x-large'])
        pltx.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        return fig


    def Derive_parameters(self,external_parameters=None, external_param_map=None):
        '''
            return derived_quantities
            Compute derived quantities as an array in order :
            a_p i_I e_I omega_I  tperi_I P_I a_b i_O e_O omega_O  tperi_O P_O M_p M_I v_parallelpm v_totpm

            angles are in degrees.

            if external_parameters is None, use the current internal parameters.

            external_param_map : map iext in external_parameters -> external_parameter_map[iext] = iin in internal representation ; If not provided, assumed to be identity.
        '''
        parammap = self.Get_fitted_parameter_map()
        inversparammap = np.ones(parammap.max()+1,dtype=int) * (parammap.size ) # the multiplication factor ensures that any call to a parameter that is not in the chain will provoke an error.
        inversparammap[parammap] = np.arange(parammap.size, dtype=int)
        parscale = self.Get_parameter_scales()

        parref = self.Get_reference_parameters()

        parini = self.Get_parameters()
        pars = parini.copy()

        relshift = np.zeros(len(self.Get_fitted_parameters_list()))

        if external_parameters is not None:
            if external_param_map is not None:
                for i in range(external_parameters.size):
                    pars[external_param_map[i]] = external_parameters[i]
            else :
                pars[:external_parameters.size] = external_parameters.copy()

        # print(parref)
        # print(parini)
        # print(pars)
        # print(external_param_map)
        for i in range(parammap.size):
            relshift[i] = (pars[parammap[i]] - parref[parammap[i]]) / parscale[parammap[i]]

        self.Set_fitted_parameter_relativeshifts(relshift)
        self.Initialise_parameters()

        paramsderived = np.ndarray(20) # truespinfreq truespinfreq1 a_p i_I e_I omega_I  tperi_I P_I a_b i_O e_O omega_O  tperi_O P_O M_p M_I vparallel vtot

        paramsderived[0] = np.sqrt(pars[3]**2 + pars[5]**2)                 # a_p
        paramsderived[1] = np.arcsin(pars[3]/paramsderived[0]) * 180./ pi   # i_I
        paramsderived[2] = np.sqrt(pars[2]**2 + pars[4]**2)                 # e_I
        paramsderived[3] = np.arctan2(pars[2]/paramsderived[2], pars[4]/paramsderived[2]) * 180./(pi) # omega_I
        paramsderived[4] = pars[17] + pars[18] # Omega_I
        paramsderived[5] = pars[6] + paramsderived[3] / (360.) * pars[7]  # tperi_I

        paramsderived[7] = np.sqrt(pars[9]**2 + pars[11]**2)                 # a_b
        paramsderived[8] = np.arcsin(pars[9]/paramsderived[7]) * 180./ pi   # i_O
        paramsderived[9] = np.sqrt(pars[8]**2 + pars[10]**2)                 # e_O
        paramsderived[10] = np.arctan2(pars[8]/paramsderived[9], pars[10]/paramsderived[9]) * 180./(pi) # omega_O
        paramsderived[11] = pars[12] + paramsderived[10] / (360.) * pars[13]   # tperi_O

        mp, mi, mo = self.Get_masses()

        if self.Get_parameter_set() <=4:
            paramsderived[13] = mi                              # m_i
            paramsderived[14] = mo                              # m_o
        else :
            paramsderived[13] = mp                              # m_p
            paramsderived[14] = mi                              # m_i
            paramsderived[15] = mo                              # m_o


        paramsderived[6] = np.sqrt(4.*pi**2 *paramsderived[0]**3 / ( daysec**2 * GMsol) / paramsderived[14]**3) * (mp+mi)         # P_I
        paramsderived[12] = np.sqrt(4.*pi**2 *paramsderived[7]**3 / ( daysec**2 * GMsol) / pars[16]**3) * (mp+mi + pars[16])     #P_O

        paramsderived[16] = pars[30] * pars[27] * clight * arcsec /1.e6  # proper motion vparallel km/s
        paramsderived[17] = np.sqrt((pars[28])**2 + pars[29]**2)*arcsec * pars[27] * clight /1.e6  # proper motion plane of sky in km/s
        #np.sqrt((pars[28])**2 + pars[29]**2)/(1000.*3600.) * (pi/180.) * pars[27] * clight/1000.   # proper motion plane of sky in km/s

        truespinfreq, truespinfreq1 = self.Get_true_spinfreq()
        paramsderived[18] = truespinfreq / daysec # true_spinfreq
        paramsderived[19] = truespinfreq1 /daysec**2 # truespinfreq1

        #self.Set_parameters(parini)

        return paramsderived


    def Plot_all_conserved_quantities_NandPN(self):
        '''
            Identical to Plot_all_conserved_quantities but with two columns, one for Newtonian and one for PostNewtonian
        '''
        if self.Get_theory() == 1 :
            self.Set_theory(0)
            self.Compute_lnposterior()

        self.Compute_integrals_of_motion()

        tinterp = self.Get_tinterp()
        fig = figure(figsize=(15,12))
        plte = fig.add_subplot(321)
        plti = fig.add_subplot(323)
        pltx = fig.add_subplot(325)
        pltepn = fig.add_subplot(322)
        pltipn = fig.add_subplot(324)
        pltxpn = fig.add_subplot(326)

        fig.text(0.20, 0.97, 'Newtonian', fontsize=self.fontsizes['x-large'])
        fig.text(0.65, 0.97, 'Post-Newtonian (1PN)', fontsize=self.fontsizes['x-large'])
        fig.tight_layout(rect=(0,0,1,0.90), pad =1.11)
        fig.set_tight_layout(True)

        x = self.Get_energies()
        x = ( x - x[0] ) / x[0]
        plte.plot(tinterp, x, label = '(E - Eini) / Eini')
        #plte.set_xlabel('Times of interpolation (MJD)')
        plte.set_ylabel('Error ')
        plte.legend(fontsize=self.fontsizes['large'], loc = 1)
        plte.xaxis.label.set_fontsize(self.fontsizes['large'])
        plte.yaxis.label.set_fontsize(self.fontsizes['large'])
        plte.tick_params(axis='both', which='major', labelsize=self.fontsizes['large'])
        plte.tick_params(axis='both', which='minor', labelsize=self.fontsizes['large'])
        plte.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        x = self.Get_center_of_mass_impulsions()
        ms = self.Get_masses()
        x  = x / (ms.sum() * clight)
        plti.plot(tinterp, x[:,0], label = 'Px/Mc')
        plti.plot(tinterp, x[:,1], label = 'Py/Mc (to Earth)')
        plti.plot(tinterp, x[:,2], label = 'Pz/Mc')
        #plti.set_xlabel('Times of interpolation (MJD)')
        plti.set_ylabel('Momentum/Mc')
        plti.legend(fontsize=self.fontsizes['large'], loc = 2)
        plti.xaxis.label.set_fontsize(self.fontsizes['large'])
        plti.yaxis.label.set_fontsize(self.fontsizes['large'])
        plti.tick_params(axis='both', which='major', labelsize=self.fontsizes['large'])
        plti.tick_params(axis='both', which='minor', labelsize=self.fontsizes['large'])
        plti.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        x = self.Get_center_of_mass_positions()
        x = (x - x[0])/clight * 10.**6
        pltx.plot(tinterp, x[:,0], label = 'x')
        pltx.plot(tinterp, x[:,1], label = 'y (to Earth)')
        pltx.plot(tinterp, x[:,2], label = 'z')
        pltx.set_xlabel('Times of interpolation (in days from reference day)')
        pltx.set_ylabel('Position (light-$\mu$s)')
        pltx.legend(fontsize=self.fontsizes['large'], loc =2)
        pltx.xaxis.label.set_fontsize(self.fontsizes['large'])
        pltx.yaxis.label.set_fontsize(self.fontsizes['large'])
        pltx.tick_params(axis='both', which='major', labelsize=self.fontsizes['large'])
        pltx.tick_params(axis='both', which='minor', labelsize=self.fontsizes['large'])
        pltx.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis


        self.Set_theory(1)
        self.Compute_lnposterior()
        self.Compute_integrals_of_motion()

        x = self.Get_energies()
        x = ( x - x[0] ) / x[0]
        pltepn.plot(tinterp, x, label = '(E - Eini) / Eini')
        #pltepn.set_xlabel('Times of interpolation (MJD)')
        pltepn.set_ylabel('Error ')
        pltepn.legend(fontsize=self.fontsizes['large'], loc = 1)
        pltepn.xaxis.label.set_fontsize(self.fontsizes['large'])
        pltepn.yaxis.label.set_fontsize(self.fontsizes['large'])
        pltepn.tick_params(axis='both', which='major', labelsize=self.fontsizes['large'])
        pltepn.tick_params(axis='both', which='minor', labelsize=self.fontsizes['large'])
        pltepn.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        x = self.Get_center_of_mass_impulsions()
        ms = self.Get_masses()
        x  = x / (ms.sum() * clight)
        pltipn.plot(tinterp, x[:,0], label = 'Px/Mc')
        pltipn.plot(tinterp, x[:,1], label = 'Py/Mc (to Earth)')
        pltipn.plot(tinterp, x[:,2], label = 'Pz/Mc')
        #pltipn.set_xlabel('Times of interpolation (MJD)')
        pltipn.set_ylabel('Momentum/Mc')
        pltipn.legend(fontsize=self.fontsizes['large'], loc = 2)
        pltipn.xaxis.label.set_fontsize(self.fontsizes['large'])
        pltipn.yaxis.label.set_fontsize(self.fontsizes['large'])
        pltipn.tick_params(axis='both', which='major', labelsize=self.fontsizes['large'])
        pltipn.tick_params(axis='both', which='minor', labelsize=self.fontsizes['large'])
        pltipn.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis

        x = self.Get_center_of_mass_positions()
        x = (x - x[0])/clight * 10.**6
        pltxpn.plot(tinterp, x[:,0], label = 'x')
        pltxpn.plot(tinterp, x[:,1], label = 'y (to Earth)')
        pltxpn.plot(tinterp, x[:,2], label = 'z')
        pltxpn.set_xlabel('Times of interpolation (in days from reference day)')
        pltxpn.set_ylabel('Position (light-$\mu$s)')
        pltxpn.legend(fontsize=self.fontsizes['large'], loc =2)
        pltxpn.xaxis.label.set_fontsize(self.fontsizes['large'])
        pltxpn.yaxis.label.set_fontsize(self.fontsizes['large'])
        pltxpn.tick_params(axis='both', which='major', labelsize=self.fontsizes['large'])
        pltxpn.tick_params(axis='both', which='minor', labelsize=self.fontsizes['large'])
        pltxpn.ticklabel_format(style='sci', scilimits=(0,0), axis='y') # forces scientification notation on y axis



        return fig


    def Plot_center_of_mass_impulsions(self):
        tinterp = self.Get_tinterp()
        x = self.Get_center_of_mass_impulsions()
        fig = figure(figsize=(16,8))
        plt = fig.add_subplot(111)
        plt.plot(tinterp, x[:,0], label = 'Px')
        plt.plot(tinterp, x[:,1], label = 'Py (to Earth)')
        plt.plot(tinterp, x[:,2], label = 'Pz')
        plt.set_xlabel('Times of interpolation (MJD)')
        plt.set_ylabel('Impulsion (kg.m/s)')
        plt.legend()
        return fig


    def Plot_center_of_mass_positions(self, fit_quadratic_out=False, figsize=(12,6), fontsize = 'large'):
        tinterp = self.Get_tinterp()
        x = self.Get_center_of_mass_positions()
        print()
        print("Initial center of mass position (m): ", x[0])

        labels = ['x', 'y (to Earth)', 'z']
        if fit_quadratic_out :
            altlabels = ['$X', '$Y', '$Z']
            iref = tinterp.searchsorted(0)
            x -= x[iref]
            dt = tinterp[-1] - tinterp[0]
            fitc = fitting.LinearLSQFitter()
            model = models.Polynomial1D(2, domain=[tinterp[0], tinterp[-1]], window=[tinterp[0], tinterp[-1]])
            for k in range(3):
                model.c0.value = x[iref, k]
                model.c1.value = (x[-1, k] - x[0, k])/ dt
                model.c2.value = 0
                fitres = fitc(model, tinterp, x[:,k])
                x[:,k] -=  fitres(tinterp)
                strc = ['{:.2g}'.format(fitres.c0.value), '{:.2g}'.format(fitres.c1.value), '{:.2g}'.format(fitres.c2.value)]
                for j in range(3):
                    if strc[j].rfind('e') > 0 :
                        strc[j] = strc[j].replace('e', r'\cdot 10^{') + r'}'

                labels[k] = altlabels[k] + ' - \\left( {:}  + {:} t + {:} t^2 \\right)$'.format(strc[0], strc[1], strc[2])
                print("Fit parameters of component {:} : ".format(k), fitres)

        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)
        plt.plot(tinterp, x[:,0], label = labels[0])
        plt.plot(tinterp, x[:,1], label = labels[1])
        plt.plot(tinterp, x[:,2], label = labels[2])
        plt.set_xlabel('Times of interpolation (day)', fontsize = fontsize)
        plt.set_ylabel('Position shift (m)', fontsize = fontsize)
        plt.xaxis.label.set_fontsize(fontsize)
        plt.yaxis.label.set_fontsize(fontsize)
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        plt.tick_params(axis='both', which='minor', labelsize=fontsize)
        fig.set_tight_layout(True)
        plt.legend(loc=0, fontsize = fontsize)
        return fig

    def Plot_energies(self, figsize=(12,6), fontsize='large'):
        tinterp = self.Get_tinterp()
        x = self.Get_energies()
        x = ( x - x[0] ) / x[0]

        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)
        plt.plot(tinterp, x, label = '$(E - E_\mathrm{ini}) / E_\mathrm{ini}$')
        plt.set_xlabel('Times of interpolation (days)')
        plt.set_ylabel('Relative error ')
        plt.legend(loc = 0, fontsize = fontsize)
        plt.xaxis.label.set_fontsize(fontsize)
        plt.yaxis.label.set_fontsize(fontsize)
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        plt.tick_params(axis='both', which='minor', labelsize=fontsize)
        fig.set_tight_layout(True)
        return fig


    def Compute_WRMS(self, residuals=None , errors = None) :
        '''
        wrms, wrmstempo = Compute_WRMS()
        Compute the weighted root mean square of the residuals with my formula and that of tempo2 in microsec.
        
        If "residuals" and "errors" are given use these instead of the internal values for time residuals and errors. 
        '''
        if residuals is None:
            res  = self.Get_time_residuals()
            sig  = self.Get_errors()
            toas = self.Get_toas()
        else :
            res = residuals
            sig=errors

        wrms = np.sqrt( np.sum(res**2 / sig) / np.sum(1./sig) - ( sum(res / sig) / np.sum(1./sig) )**2 )
        wrmstempo = np.sqrt( (np.sum(res**2 / sig**2)  - ( sum(res / sig**2))**2 / np.sum(1./sig**2)) / np.sum(1./sig**2) )

        return wrms, wrmstempo
    
    def Compute_weighted_average(self, residuals=None, errors=None):
        '''
        Return weighted average in microsec
        
        Compute the weighted average assuming a Gaussian uncertainty on Toas (weight = 1/sigma**2, where sigma is the Gaussian RMS). 
        
        If "residuals" and "errors" are given use these instead of the internal values for time residuals and errors. 
        '''
        if residuals is None:
            res  = self.Get_time_residuals()
            sig  = self.Get_errors()
            toas = self.Get_toas()
        else :
            res = residuals
            sig=errors

        av = np.average(res, weights=1/sig**2)
                        
        return av
    

    def Plot_time_residuals(self, show_both_wrms=True, markersize = 4, figsize=(16,8), mjd=False, show_line=False, show_error_bars=True, show_mean=False, plt =None):
        '''
        fig = Plot_time_residuals(self, show_both_wrms=True, markersize=4)
        Return a matplotlib figure and plot object showing the time residuals for the current fitted parameters.
        If show_both_wrms is True, put the two computations of wrms (see Compute_WRMS) in the label, otherwiseput only Tempo's.

        figsize : tuple (width, height) as taken from the matplotlib figure class
        mjd : if True, the horizontal axis is in MJD, otherwise in days from the first toas
        markersize : size of the dots
        show_line: plot a line between toas
        show_error_bars : if True, add the error bars read from the timfile to the plot for each TOA.
        show_mean : Plot a horizontal red line at the weighted average of residuals (see Compute_weighted_average). If show_mean is a positive integer then this value is used as line width. 
        plt : axis object to draw upon. If None, return a figure instead.
        '''
        if plt is None:
            noplt = True
            
        res  = self.Get_time_residuals()
        if show_error_bars:
            sig  = self.Get_errors()
        else:
            sig = np.zeros(res.size)
        toas = self.Get_toas()
        if mjd is True :
            toas += self.Get_timeshift()

        marginx = (toas[-1] - toas[0]) * 0.05

        #wrms = np.sqrt( np.sum(res**2 / sig) / np.sum(1./sig) - ( sum(res / sig) / np.sum(1./sig) )**2 )
        wrms, wrmstempo =  self.Compute_WRMS()
        if show_both_wrms == True :
            label = 'Timing residuals (WRMS = $%.3f\mathrm{\mu s}$, WRMStempo = $%.3f\mathrm{\mu s}$)'%(wrms,wrmstempo)
        else :
            label = 'Timing residuals (WRMS = $%.3f\mathrm{\mu s}$)'%(wrmstempo)
        if noplt:
            fig = figure(figsize=figsize)
            plt = fig.add_subplot(111)
        if show_line is True:
            fmt = '.-'
        else:
            fmt = '.'
        plt.errorbar(toas, res, yerr=sig, fmt = fmt, markersize= markersize,  label = label) #'Timing residuals (WRMS = $%.2f\mathrm{\mu s}$)'%(wrmstempo))#
        print(" Should we take sig or sig /2 for the error bars ? Seems tempo2 is taking sig")
        if mjd is True :
            plt.set_xlabel('Time of arrival (Modified Julian Day)')
        else :
            plt.set_xlabel('Time since first TOA (in days)')
        plt.set_ylabel('Residuals ($\mu$s)')
        
        if show_mean:
            plt.axhline(y=self.Compute_weighted_average(), color='red', linewidth=np.int(show_mean))
            
        # To make text large
        plt.legend(fontsize=self.fontsizes['xx-large'])
        plt.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        plt.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        plt.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        plt.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
        plt.set_xlim(toas[0] - marginx, toas[-1] + marginx)

        if noplt:
            fig.set_tight_layout(True)
            return fig
        else:
            return plt



    def Plot_time_residuals_vs_orbitphases_moy(self, ni =1, no = 1, nt=1, alarmlowpts=500):
      '''
            return plot of residuals versus inner phase, outer phase and earth phase after averaging over a window of size Ti / ni , To/no and Tt/nt where t are the orbital periods.
            if the number of points in the window is lower than alarmlowpts, a warning is printed.
            matplotlib figure = Plot_time_residuals_vs_orbitphases(self,ni =1, no = 1, nt=1, alarmlowpts=500)
            This routine can be greedy, need to be otpimized

      '''
      toas = self.Get_toas()
      res = self.Get_time_residuals()
      params = self.Get_parameters()

      Ti = params[7]
      omi = params[4]
      To = params[13]
      omo = params[10]
      Tt = 365.25
      tperii = params[6]
      tperio = params[12]
      phasei = np.zeros(len(res))
      phaseo = np.zeros(len(res))
      phaset = np.zeros(len(res))
      if tperii > toas[0] : tperii = tperii - Ti * ( floor( (tperii - toas[0]) / Ti ) + 1 )
      if tperio > toas[0] : tperio = tperio - To * ( floor( (tperio - toas[0]) / To ) + 1 )
      #phasei = (toas - tperii) % Ti
      #phaseo = (toas - tperio) % To
      phasei = (toas - toas[0]) % Ti
      phaseo = (toas - toas[0]) % To
      phaset = (toas - toas[0] ) % Tt
      #Tiss = ( Ti + (3.*pi/2. - omi) * Ti / (2. * pi ) )  % Ti
      #Toss = (To + (3.*pi/2. - omo) * To / (2.*pi) ) % To
      fig = figure(figsize=(16,8))
      plti = fig.add_subplot(311)
      plto = fig.add_subplot(312)
      pltt = fig.add_subplot(313)
      pltt.set_xlabel('Earth orbital phase (MJD)')
      pltt.set_ylabel('Residuals ($\mu$s)')
      plti.set_xlabel('Inner orbital phase (MJD)')
      plti.set_ylabel('Residuals ($\mu$s)')
      #plti.plot([Tiss, Tiss],[res.min(),res.max()], 'r-', label='Solar system direction')
      #plto.plot([Toss, Toss],[res.min(),res.max()], 'r-', label='Solar system direction')
      plto.set_xlabel('Outer orbital phase (MJD)')
      plto.set_ylabel('Residuals ($\mu$s)')

      pmoyi = np.zeros(ni)
      pmoyo = np.zeros(no)
      pmoyt = np.zeros(nt)
      moyi = np.zeros(ni)
      moyo = np.zeros(no)
      moyt = np.zeros(nt)

      args = phasei.argsort()
      resi = res[args]
      phasei = phasei[args]
      args = phaseo.argsort()
      reso = res[args]
      phaseo = phaseo[args]
      args = phaset.argsort()
      rest = res[args]
      phaset = phaset[args]

      for i in range(ni): pmoyi[i] = (Ti/ni * i + Ti/ni/2.)
      for i in range(no): pmoyo[i] = (To/no * i + To/no/2.)
      for i in range(nt): pmoyt[i] = (Tt/nt * i + Tt/nt/2.)
      for i in range(ni):
          pos1 = phasei.searchsorted(i*Ti/ni)
          pos2 = phasei.searchsorted((i+1)*Ti/ni)
          if (pos2 - pos1) < alarmlowpts : print("Too few points in phasei, range number %i : %i"%(i, pos2 - pos1))
          moyi[i] = resi[pos1 : pos2].mean()
      for i in range(no):
          pos1 = phaseo.searchsorted(i*To/no)
          pos2 = phaseo.searchsorted((i+1)*To/no)
          if (pos2 - pos1) < alarmlowpts : print("Too few points in phaseo, range number %i : %i"%(i, pos2 - pos1))
          moyo[i] = reso[pos1 : pos2].mean()
      for i in range(nt):
          pos1 = phaset.searchsorted(i*Tt/nt)
          pos2 = phaset.searchsorted((i+1)*Tt/nt)
          if (pos2 - pos1) < alarmlowpts : print("Too few points in phaset, range number %i : %i"%(i, pos2 - pos1))
          moyt[i] = rest[pos1 : pos2].mean()

      plti.plot(pmoyi,moyi, '.-')
      plto.plot(pmoyo,moyo, '.-')
      pltt.plot(pmoyt,moyt,'.-')
      plti.legend()
      plto.legend()
      pltt.legend()
      return fig

    def Get_distance_from_sun(self):
        '''
            Return the distance between the sun and the pulsar for each TOA in degrees. 
        '''
        toas = self.Get_toas()
        toas += self.Get_timeshift()
        toas = apt.Time(toas, format="mjd")
        coords = self.Get_parameters()[25:27]
        psr = apco.SkyCoord(coords[0]*apu.radian, coords[1]*apu.radian, frame="icrs")
        dists = apco.get_sun(toas).separation(psr) # WARNING : psr.separation(apco.get_sun(toas)) returns crap !
        dists = np.array(dists.value)
        return dists

    def Create_masked_tim_file(self, filename=None, max_residual=np.inf, min_sun_distance=0.):
        '''
            Creates a mask of the data according to the angular distance between the pulsar and the sun and the amplitude of the residual at a given date 
            If residual > max_residual or distance_to_sun < min_sun_distance the mask is set to 1, 0 otherwise.
            If a filename is given a new timfile is created where all the lines with mask=1 are commented. 
            
            Parameters : 
                - filename : path to timfile
                - max_residual : maximum value before a toa is considered an outlier in microsec
                - min_sun_distance : minimum angular distance under which the sun is considered too close to the line of sight (in degrees)
            
            Returns :
                mask : array of 0 and 1. 
        '''
        residuals = self.Get_time_residuals()
        mask = np.zeros(residuals.size, dtype=int)
        sun_dist = self.Get_distance_from_sun()
        
        mask = np.where(np.abs(residuals)<max_residual, 0, 1) 
        mask = np.where(sun_dist > min_sun_distance, mask, 1) 
        
        if filename is not None :
            self.Save_masked_copy_of_tim_file(filename, mask)
            
        return mask

    def Plot_time_residuals_vs_sun_distance(self, window=None, show_error_bars=True, markersize=4, figsize=(16,8)):
        '''
            return fig

            window : if not None then add to the plot the mean obtained by a sliding window (in degrees)
            figsize : tuple (width, height) as taken from the matplotlib figure class
            markersize : size of the dots
            show_error_bars : if True, add the error bars read from the timfile to the plot for each TOA.
        '''
        res  = self.Get_time_residuals()
        res -= res.mean()
        if show_error_bars:
            sig  = self.Get_errors()
        else:
            sig = np.zeros(res.size)
        toas = self.Get_toas()
        toas += self.Get_timeshift()

        toas = apt.Time(toas, format="mjd")
        coords = self.Get_parameters()[25:27]
        psr = apco.SkyCoord(coords[0]*apu.radian, coords[1]*apu.radian, frame="icrs")
        dists = apco.get_sun(toas).separation(psr) # WARNING : psr.separation(apco.get_sun(toas)) returns crap !
        dists = np.array(dists.value)
        #marginx = (toas[-1] - toas[0]) * 0.05

        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)
        fmt = '.'
        plt.errorbar(dists, res, yerr=sig, fmt = fmt, markersize= markersize) #'Timing residuals (WRMS = $%.2f\mathrm{\mu s}$)'%(wrmstempo))#

        if window is not None :
            distargsort = dists.argsort()
            distsorted = dists[distargsort]
            ressorted = res[distargsort]
            avpts = dists.min() + np.arange(np.int(np.floor((dists.max()-dists.min())/window)))* window + window
            winres = np.zeros(avpts.size)
            oldarg=0
            for i in range(avpts.size):
                newarg = distsorted.searchsorted(avpts[i])
                winres[i] = ressorted[oldarg:newarg].mean()
                oldarg = newarg
            plt.plot(avpts - 0.5*window, winres, '.-', color='red', linewidth=3, label = 'Window {:.1g}deg'.format(window))
            plt.legend(fontsize=self.fontsizes['xx-large'])

        plt.set_xlabel('Angular separation between the Sun and the pulsar (in degrees)')
        plt.set_ylabel('Residuals ($\mu$s)')
        # To make text large
        plt.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        plt.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        plt.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        plt.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])

        fig.set_tight_layout(True)

        return fig

    def Plot_time_residuals_vs_orbitphases(self, idphase=False, innerphaseranges= None, outerphaseranges = None, markersize=2, external_residuals=None):
      '''
            matplotlib figure = Plot_time_residuals_vs_orbitphases(self, idphase=False, innerphaseranges= None, outerphaseranges = None)
            This routine can be greedy, need to be otpimized
            If idphase = True, then display each revolution with a different color going from green (first) to red (last)
            If innerphasesranges or outerphaseranges are specified as array then two possibilities :
                - phaseranges is 1D, then it contains all the revolution numbers you want to display
                - phaseranges.shape=(n,2) then it contains n doublets of [first revolution, last revolution], all the toas within that range being displayed with the same color.

             markersize : size of the dots

            external_residuals : if not None, then must be an array of residuals of same size as self.Get_time_residuals() that replaces internally calculated residuals.
      '''
      toas = self.Get_toas()
      if external_residuals is None:
          res = self.Get_time_residuals()
      else:
          res = external_residuals.copy()
      params = self.Get_parameters()

      Ti = params[7]
      omi = params[4]
      To = params[13]
      omo = params[10]
      Tt = 365.25
      tperii = params[6]
      tperio = params[12]
      phasei = np.zeros(len(res))
      phaseo = np.zeros(len(res))
      phaset = np.zeros(len(res))
      if tperii > toas[0] : tperii = tperii - Ti * ( floor( (tperii - toas[0]) / Ti ) + 1 )
      if tperio > toas[0] : tperio = tperio - To * ( floor( (tperio - toas[0]) / To ) + 1 )
      phasei = (toas - tperii) % Ti
      phaseo = (toas - tperio) % To
      phaset = (toas - toas[0] ) % Tt
      print("tperio", tperio, toas[0], To)
      if idphase :
          orbitnbi = np.ndarray(len(toas), dtype=np.int)
          orbitnbo = np.ndarray(len(toas), dtype=np.int)
          for i in range(len(toas)) :
            orbitnbi[i] = np.int((toas[i] - tperii - phasei[i]) / Ti )
            orbitnbo[i] = np.int( (toas[i] - tperio  - phaseo[i] ) / To)
          orbitnbi -= orbitnbi[0]
          orbitnbo -= orbitnbo[0]
          nbi = np.float(orbitnbi[-1])
          nbo = np.float(orbitnbo[-1])
      Tiss = ( Ti + (3.*pi/2. - omi) * Ti / (2. * pi ) )  % Ti
      Toss = (To + (3.*pi/2. - omo) * To / (2.*pi) ) % To
      print('tiss ', (3.*pi/2. - omi) * Ti / (2. * pi ), Tiss)
      print(omi, Ti, tperii)
      fig = figure(figsize=(16,11))
      plti = fig.add_subplot(311)
      plto = fig.add_subplot(312)
      pltt = fig.add_subplot(313)
      pltt.set_xlabel('Earth orbital phase (days)')
      pltt.set_ylabel('Residuals ($\mu$s)')
      pltt.plot(phaset,res,'.', markersize =markersize)
      plti.set_xlabel('Inner orbital phase (days)')
      plti.set_ylabel('Residuals ($\mu$s)')
      plti.plot([Tiss, Tiss],[res.min(),res.max()], 'r-', label='Solar system direction')
      if not idphase :
          plti.plot(phasei,res, '.',markersize =markersize)
          plto.plot(phaseo,res, '.',markersize =markersize)
      else :
          print("Number of inner revolutions : ", nbi + 1)
          print("Number of outer revolutions : ", nbo + 1)
          if innerphaseranges is not None :

              if type(innerphaseranges) is int :
                  nbranges = innerphaseranges
                  innerphaseranges = np.ndarray((nbranges + 1, 2), dtype=int)
                  for i in range(nbranges) :
                      innerphaseranges[i,0] = floor( nbi / nbranges ) * i
                      innerphaseranges[i,1] = floor( nbi / nbranges ) * (i + 1)
                  innerphaseranges[nbranges,0] = floor( nbi / nbranges ) * nbranges
                  innerphaseranges[nbranges,1] = nbi
                  print("Generated inner phase ranges (in number of inner revolutions achieved) : ")
                  print(innerphaseranges)
                  print()

              if len(innerphaseranges.shape) >= 2 :
                  for rangenb in range(innerphaseranges.shape[0]) :
                    innerphasenb = np.arange( innerphaseranges[rangenb, 0], innerphaseranges[rangenb, -1] )
                    deb = orbitnbi.searchsorted(innerphaseranges[rangenb, 0], side='left')
                    fin = orbitnbi.searchsorted(innerphaseranges[rangenb, -1], side='right')
                    plti.plot(phasei[deb:fin], res[deb:fin], '.', markersize =markersize,color = (rangenb/max(1.,np.float(innerphaseranges.shape[0] -1)), 1. - rangenb/max(1.,np.float(innerphaseranges.shape[0] - 1)), 0.) )
              else :
                innerphasenb = innerphaseranges
                for ip in range(len(innerphasenb)) :
                    deb = orbitnbi.searchsorted(innerphasenb[ip], side='left')
                    if innerphasenb[ip] + 1 == orbitnbi[-1] :
                        fin = len(toas)
                    else :
                        fin = orbitnbi.searchsorted(innerphasenb[ip] + 1, side='right')
                        for i in range(deb, fin):
                            plti.plot(phasei[i], res[i], '.', markersize =markersize,color = (ip/np.float(innerphaseranges.shape[0] - 1), 1. - ip/np.float(innerphaseranges.shape[0] - 1), 0.) )
          else :
            for i in range(len(toas)) :
                plti.plot(phasei[i], res[i], '.',markersize =markersize, color = (orbitnbi[i]/nbi, 1. - orbitnbi[i]/nbi, 0.) )


          if outerphaseranges is not None :

              if type(outerphaseranges) is int :
                  nbranges = outerphaseranges
                  outerphaseranges = np.ndarray((nbranges + 1, 2), dtype=int)
                  for i in range(nbranges) :
                      outerphaseranges[i,0] = floor( nbo / nbranges ) * i
                      outerphaseranges[i,1] = floor( nbo / nbranges ) * (i + 1)
                  outerphaseranges[nbranges,0] = floor( nbo / nbranges ) * nbranges
                  outerphaseranges[nbranges,1] = nbo
                  print("Generated outer phase ranges (in number of outer revolutions achieved) : ")
                  print(outerphaseranges)
                  print()

              if len(outerphaseranges.shape) >= 2 :
                  for rangenb in range(outerphaseranges.shape[0]) :
                    outerphasenb = np.arange( outerphaseranges[rangenb, 0], outerphaseranges[rangenb, -1] )
                    deb = orbitnbo.searchsorted( outerphaseranges[rangenb, 0], side='left')
                    fin = orbitnbo.searchsorted( outerphaseranges[rangenb, -1] + 1, side='right')
                    plto.plot(phaseo[deb:fin], res[deb:fin], '.', markersize =markersize,color = (rangenb/max(1.,np.float(innerphaseranges.shape[0] -1)), 1. - rangenb/max(1.,np.float(innerphaseranges.shape[0] -1)), 0.) )
              else :
                outerphasenb = outerphaseranges
                for op in range(len(outerphasenb)) :
                    deb = orbitnbo.searchsorted(outerphasenb[op], side='left')
                    if outerphasenb[op] + 1 == orbitnbo[-1] :
                        fin = len(toas)
                    else :
                        fin = orbitnbo.searchsorted(outerphasenb[op] + 1, side='right')
                    plto.plot(phaseo[deb:fin], res[deb:fin], '.', markersize =markersize,color = (op/np.float(outerphaseranges.shape[0]), 1. - op/np.float(outerphaseranges.shape[0]), 0.) )
          else :
              for op in range( orbitnbo[-1] + 1 ) :
                    deb = orbitnbo.searchsorted(op, side='left')
                    if op + 1 == orbitnbo[-1] :
                        fin = len(toas)
                    else :
                        fin = orbitnbo.searchsorted(op + 1, side='right')
                    plto.plot(phaseo[deb:fin], res[deb:fin], '.', markersize =markersize, color = (op/(nbo), 1. - op/nbo, 0.) )
            #for i in range(len(toas)) :
                #plto.plot(phaseo[i], res[i], '.', color = (orbitnbo[i]/nbo, 1. - orbitnbo[i]/nbo, 0. ) )

      plto.plot([Toss, Toss],[res.min(),res.max()], 'r-', label='Solar system direction')
      plto.set_xlabel('Outer orbital phase (days)')
      plto.set_ylabel('Residuals ($\mu$s)')

      # Find the point of closes approach to the sun
      psr = self.Get_SkyCoord()
      ts = np.linspace(0, 365.25, num = 1000)
      apts = apt.Time(ts+ self.Get_timeshift() + toas.min(), format="mjd")
      dists = np.array(apco.get_sun(apts).separation(psr).value)
      closest_approach = ts[dists.argmin()]
      pltt.plot([closest_approach,closest_approach],[res.min(),res.max()], 'r-', label='PSR closest to the Sun')

      plti.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
      plti.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
      plti.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
      plti.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
      plto.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
      plto.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
      plto.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
      plto.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
      pltt.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
      pltt.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
      pltt.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
      pltt.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])

      plti.legend(loc=0, fontsize=self.fontsizes['large'])
      plto.legend(loc=0, fontsize=self.fontsizes['large'])
      pltt.legend(loc=0, fontsize=self.fontsizes['large'])

      fig.set_tight_layout(True)

      return fig


    def Get_SkyCoord(self):
        '''
        return an astropy.coordinates.SkyCoord object corresponding to the position of the pulsar at posepoch
        '''
        coords = self.Get_parameters()[25:27]
        psr = apco.SkyCoord(coords[0]*apu.radian, coords[1]*apu.radian, frame="icrs")
        return psr

    def Plot_Lombscargle_of_residuals(self, sampling_factor=1, loglog = True, maxfreq = None, minfreq=None,
                                        figsize=(15,10), linewidth=1.5, normalize = 'standard', external_residuals=None,
                                        max_low_freq=1/365.25, show_low_peak=5, show_low_plus_earth=5, extra_frequencies=[], 
                                        mean_noise=None, min_noise=None, max_noise=None, noise_legends=['Mean noise', 'Min noise', 'Max noise'], 
                                        legend_kwargs={'loc':'best', 'fontsize':"large"}, implementation='astropy'
                                        ):
        '''
        external_residuals : if not None, then must be an array of residuals of same size as self.Get_time_residuals() that replaces internally calculated residuals.
        max_low_freq : gives the maximum frequency below which to search for a low frequency peak (i.e. red noise fondamental)
        show_low_peak : if >0 show the lowest frequency peak and its harmonics show_low_peak -1 harmonics
        show_low_plus_earth : idem with earth frequency + (show_low_plus_earth -1 ) lowest frequency.
        extra_frequencies : list of frequencies that will be shown as vertical dashed brown lines
        mean_noise : Curve corresponding to the theoretical mean (the corresponding frequencies are in self.fs, after Plot_Lombscargle_of_residuals has been run)
        max_noise : curve for which only 2.5% of realisations are above (the corresponding frequencies are in self.fs, after Plot_Lombscargle_of_residuals has been run)
        max_noise : curve for which only 2.5% of realisations are below (the corresponding frequencies are in self.fs, after Plot_Lombscargle_of_residuals has been run)
        noise_legends : legends of the noise curves
        implementation: either 'scipy' for compatibility with version until 13/05/2022 or 'astropy'. 
        normalize : if implementation 'scipy' then can be True or False. True mean normalisation akin to Fourier transform. If implementation "astropy", then can be any of astropy's normalizations 'standard', 'psd', 'model', 'log'. 
        '''
        ts = self.Get_toas()
        if external_residuals is None:
            res = self.Get_time_residuals()
        else:
            res = external_residuals.copy()
        #res -= res.mean()
        errs = self.Get_errors()
        
        fi = 1./ self.Get_parameters()[7]
        fo = 1./ self.Get_parameters()[13]
        ft = 1./ 365.25

        deltat = ts.max() - ts.min()
        dt = np.abs((ts[1:] - ts[:ts.size-1]).min())
        df = 1. / (deltat * sampling_factor)

        if maxfreq is None :
            maxf = 10*fi
        else : 
            maxf = maxfreq
        if minfreq is None :
            minf = 1./deltat
        else : 
            minf = minfreq
        
        # Compute LombScargle
        self.fs = minf + np.arange(np.int(maxf/df)) * df
        if implementation == 'scipy':
            self.periodogram  = lombscargle(ts, res, 2*np.pi*self.fs)
            labely = r"Power (scipy normalization)"
            if normalize:
                self.periodogram /= 0.25*ts.size
                labely = r"Power ($\mu$s${}^2$)"
        else : 
            self.periodogram  = LombScargle(ts, res, errs, normalization=normalize).power(self.fs)
            labely = r"Power (astropy {:} normalization)".format(normalize)
        
        # Plot figure
        self.fig = figure(figsize=figsize)
        self.plt = self.fig.add_subplot(111)
        if loglog :
            self.plt.loglog(self.fs, self.periodogram, '-', linewidth=linewidth)
        else :
            self.plt.plot(self.fs, self.periodogram, '-', linewidth=linewidth)

        if deltat > 365.25 :
            flowfpeak = self.fs[self.periodogram[:self.fs.searchsorted(max_low_freq)].argmax()]
            if show_low_peak > 0:
                self.plt.axvline(flowfpeak, linestyle='--', label='$f_{{\mathrm{{RN}}}}={:.3e}$'.format(flowfpeak), color='black', linewidth=linewidth)
                if show_low_peak >1:
                    self.plt.axvline(2*flowfpeak, linestyle=':', label='$f_{\mathrm{RN}}$ harmonics', color='black', linewidth=linewidth)
                for i in range(3, show_low_peak):
                    self.plt.axvline(i*flowfpeak, linestyle=':', color='black', linewidth=linewidth)
                if show_low_plus_earth > 0:
                    self.plt.axvline(ft + -flowfpeak, linestyle='-.', color='gray', linewidth=linewidth)
                    self.plt.axvline(ft + flowfpeak, linestyle='-.', label='$f_E + k f_{\mathrm{RN}}$', color='gray', linewidth=linewidth)
                    for i in range(3, show_low_plus_earth):
                        self.plt.axvline(ft + i*flowfpeak, linestyle='-.', color='gray', linewidth=linewidth)

        for freq in extra_frequencies:
            self.plt.axvline(freq, linestyle='--', label='{:.3g}'.format(freq), color='brown', linewidth=lindtewidth)
            
        if mean_noise is not None :
            if loglog:
                self.plt.loglog(self.fs, mean_noise, '--', linewidth=linewidth, color='grey', label=noise_legends[0])
            else:
                self.plt.plot(self.fs, mean_noise, '--', linewidth=linewidth, color='grey', label=noise_legends[0])
        if min_noise is not None : 
            if loglog:
                self.plt.loglog(self.fs, min_noise, ':', linewidth=linewidth, color='grey', label=noise_legends[1])   
            else:
                 self.plt.plot(self.fs, min_noise, ':', linewidth=linewidth, color='grey', label=noise_legends[1])            
        if max_noise is not None : 
            if loglog:
                self.plt.loglog(self.fs, max_noise, ':', linewidth=linewidth, color='grey', label=noise_legends[2])          
            else : 
                self.plt.plot(self.fs, max_noise, ':', linewidth=linewidth, color='grey', label=noise_legends[2])        



        self.plt.axvline(fi, linestyle='--', label='$f_I$', color='red', linewidth=linewidth)
        self.plt.axvline(fo, linestyle='--', label='$f_O$', color='orange', linewidth=linewidth)
        self.plt.axvline(ft, linestyle='--', label='$f_E$', color='green', linewidth=linewidth)
        #self.plt.axvline(2*ft, linestyle=':', label='$2f_E$', color='green')
        self.plt.axvline(2*fi - fo, linestyle='--', label='$2f_I - f_O$', color='purple', linewidth=linewidth)


        #
        # fmarks = np.array([fi, fo, ft])
        # fmarks.sort()
        # plttop = self.plt.twiny()
        # plttop.set_xticks(fmarks)
        # plttop.set_xticklabels([{:}fmark])
        self.plt.set_xlabel(r"Frequency (day${}^{-1}$)", fontsize="xx-large")
        self.plt.set_ylabel(labely, fontsize="xx-large")
        self.plt.tick_params(axis='both', which='major', labelsize=self.fontsizes['x-large'])
        self.plt.tick_params(axis='both', which='minor', labelsize=self.fontsizes['x-large'])

        self.fig.set_tight_layout(True)

        self.plt.legend(**legend_kwargs)

        self.fig.show()
        return
    

    def Plot_delays(self, nfakes=10000, edgecut = 5., zoom_aberration=False, markersize=15, legendloc = [0,0,0,0]):
        '''
            Compute delays at evenly-sampled fake TOAs using the current solution, BUT under the approximation that coupling between solar system and pulsar system motion are neglected (i.e. no Kopeikin).
            
            If you want a fully accurate decomposition of delays at quasi-evenly sampled TOAs, first generate a fake tim file using "Make_fake_timfile.exe", then use Get_delays(). 
            
            NOTE : This is different from Get_delays() which computes the delays for the actual TOAs only. 
        '''
        res  = self.Get_time_residuals()
        sig  = self.Get_errors()
        toas = self.Get_toas()
        wrmstempo = np.sqrt( (np.sum(res**2 / sig**2)  - ( sum(res / sig**2))**2 / np.sum(1./sig**2)) / np.sum(1./sig**2) )
        faketoas, geo, ein, shap, aber = self.Get_fake_bats_and_delays_interp(nfakes)

        marginx = 0.05*(faketoas[-1] - faketoas[0] )


        # remove the first and last toas to avoid interpolation issues on the edges. (alternative : have a larger faketoas span)
        print(" Removing the %f first and last days of toas. "%edgecut)
        beg = toas.searchsorted(toas[0] + edgecut)
        end = toas.searchsorted(toas[-1] - edgecut)
        toas = toas[ beg : end ]
        res = res[ beg : end]

        # ---------Geometric
        res /= 10.**6 # in seconds

        s = interpolate.InterpolatedUnivariateSpline(faketoas, geo)

        figgeo = figure(figsize=(16,8))
        pltgeo = figgeo.add_subplot(111)

        pltgeo.plot(faketoas, geo, '-', label = 'Geometric delay model')
        pltgeo.scatter(toas, res + s(toas), s = markersize, c='r', marker = 'o',   label = 'Timing residuals')

        pltgeo.set_xlabel('Time since first TOA (in days)')
        pltgeo.set_ylabel('Residuals (s)')
        # To make text large
        pltgeo.legend(loc=legendloc[0], fontsize=self.fontsizes['xx-large'])
        pltgeo.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltgeo.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltgeo.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        pltgeo.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
        pltgeo.set_xlim(faketoas[0] - marginx, faketoas[-1] + marginx)
        marginy = 0.05 * ( (res + s(toas)).max() - (res + s(toas)).min() )
        pltgeo.set_ylim((res + s(toas)).min() - marginy , (res + s(toas)).max() + marginy)

        figgeo.set_tight_layout(True)

        # ---------Einstein
        # put everything in microsec now
        ein = 10.**6 * ein
        shap = 10.**6 *shap
        aber = 10.**6 *aber
        res = 10**6 *res

        # remove the mean
        ein -= (ein[-1] - ein[0])/(faketoas[-1] - faketoas[0]) * faketoas

        s = interpolate.InterpolatedUnivariateSpline(faketoas, ein)

        figein = figure(figsize=(16,8))
        pltein = figein.add_subplot(111)

        pltein.plot(faketoas, ein, '-', label = 'Einstein delay model')
        pltein.scatter(toas, res + s(toas), s = markersize, c='r', marker = 'o',  label = 'Timing residuals')

        pltein.set_xlabel('Time since first TOA (in days)')
        pltein.set_ylabel('Residuals ($\mu$s)')
        # To make text large
        pltein.legend(loc=legendloc[1], fontsize=self.fontsizes['xx-large'])
        pltein.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltein.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltein.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        pltein.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
        pltein.set_xlim(faketoas[0] - marginx, faketoas[-1] + marginx)
        marginy = 0.05 * ( (res + s(toas)).max() - (res + s(toas)).min() )
        pltein.set_ylim((res + s(toas)).min() - marginy , (res + s(toas)).max() + marginy)

        figein.set_tight_layout(True)

        # ---------- Shapiro
        s = interpolate.InterpolatedUnivariateSpline(faketoas, shap)

        figshap = figure(figsize=(16,8))
        pltshap = figshap.add_subplot(111)

        pltshap.plot(faketoas, shap, '-', label = 'Shapiro delay model')
        pltshap.scatter(toas, res + s(toas), s = markersize, c='r', marker = 'o',  label = 'Timing residuals')

        pltshap.set_xlabel('Time since first TOA (in days)')
        pltshap.set_ylabel('Residuals ($\mu$s)')
        # To make text large
        pltshap.legend(loc=legendloc[2], fontsize=self.fontsizes['xx-large'])
        pltshap.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltshap.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltshap.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        pltshap.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
        pltshap.set_xlim(faketoas[0] - marginx, faketoas[-1] + marginx)
        marginy = 0.05 * ( (res + s(toas)).max() - (res + s(toas)).min() )
        pltshap.set_ylim((res + s(toas)).min() - marginy , (res + s(toas)).max() + marginy)

        figshap.set_tight_layout(True)

        # ---------- Aberration
        s = interpolate.InterpolatedUnivariateSpline(faketoas, aber)

        if zoom_aberration == True :
            figaber = figure(figsize=(16,12))
            pltaber = figaber.add_subplot(212)
            pltaberzoom = figaber.add_subplot(211)
        else :
            figaber = figure(figsize=(16,8))
            pltaber = figaber.add_subplot(111)

        pltaber.plot(faketoas, aber, '-', label = 'Aberration delay model')
        pltaber.scatter(toas, res + s(toas), s = markersize, c='r', marker = 'o',  label = 'Timing residuals')

        pltaber.set_xlabel('Time since first TOA (in days)')
        pltaber.set_ylabel('Residuals ($\mu$s)')
        # To make text large
        pltaber.legend(loc=legendloc[3], fontsize=self.fontsizes['xx-large'])
        pltaber.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltaber.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        pltaber.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        pltaber.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
        pltaber.set_xlim(faketoas[0] - marginx, faketoas[-1] + marginx)
        marginy = 0.05 * ( (res + s(toas)).max() - (res + s(toas)).min() )
        pltaber.set_ylim((res + s(toas)).min() - marginy , (res + s(toas)).max() + marginy)

        if zoom_aberration == True :
            pltaberzoom.plot(faketoas, aber, '-', label = 'Aberration delay model')
            pltaberzoom.scatter(toas, res + s(toas), s = markersize, c='r', marker = 'o',  label = 'Timing residuals')
            pltaberzoom.set_ylabel('Residuals ($\mu$s)')
            # To make text large
          #  pltaberzoom.legend(loc=legendloc[3], fontsize=self.fontsizes['xx-large'])
            pltaberzoom.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
            pltaberzoom.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
            pltaberzoom.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
            pltaberzoom.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])
            pltaberzoom.set_xlim(faketoas[0] - marginx, faketoas[-1] + marginx)
            marginy = 0.05 *  ( aber.max() - aber.min() )
            pltaberzoom.set_ylim(aber.min() - marginy , aber.max() + marginy)


        figaber.set_tight_layout(True)

        figgeo.show()
        figein.show()
        figshap.show()
        figaber.show()

        return figgeo, figein, figshap, figaber


    def Plot_orbits_only(self, show=True, rotate_to_psb=True, xy=[0,1]):
        '''
            figure = Plot_orbits()
            
            rotate_to_psb : if True, rotate state vectors to the Pulsar Binary frame (PSB) such that the line of sight is along the z axis
            xy : the two axis to show x axis= xy[0], y axis = xy[1]. If rotate_to_psb=True and xy=[0,1] then the plane of the plot is the plane of the sky

            Returns plots of the orbits.
        '''
        Msol = 1.988435e30         # Masse du Soleil (kg)
        c = 299792458.
        self.Compute_integrals_of_motion()
        svs = self.Get_interp_state_vectors(PSB = rotate_to_psb)
        sp = svs[0]
        sc1 = svs[1]
        sc2 = svs[2]
        cdm = self.Get_center_of_mass_positions()
        cdm /= c
        Mp, Mi, Mo = self.Get_masses()        
        
        ax = xy[0]
        ay = xy[1]
        
        dessin = figure(figsize=(12,12))
        orbs = dessin.add_subplot(111)
        #if (Mo >0.) :
            #xt = [sc2[:,0].min()/c, sc2[:,0].max()/c]
        #else :
            #xt = [sc1[:,0].min()/c, sc1[:,0].max()/c]
        #yt = [sc1[:,2].min()/c, sc1[:,2].max()/c ]
        #xt[0] += np.abs(xt[1] - xt[0])*0.1
        #xt[1] -= np.abs(xt[1] - xt[0])*0.1
        #yt[0] += np.abs(yt[1] - yt[0])*0.1
        #yt[1] -= np.abs(yt[1] - yt[0])*0.1
        orbs.plot(sp[:,ax]/c,sp[:,ay]/c,'b-')
        orbs.plot([cdm[0,ax],sp[0,ax]/c],[cdm[0,ay],sp[0,ay]/c],'b--')
        orbs.plot(sc1[:,ax]/c,sc1[:,ay]/c,'g-')
        orbs.plot([cdm[0,ax],sc1[0,ax]/c],[cdm[0,ay],sc1[0,ay]/c],'g--')
        if (Mo >0.) :
            orbs.plot(sc2[:,ax]/c,sc2[:,ay]/c,'r-')
            orbs.plot([cdm[0,ax],sc2[0,ax]/c],[cdm[0,ay],sc2[0,ay]/c],'r--')
            orbs.plot([cdm[-1,ax],sc2[-1,ax]/c],[cdm[-1,ay],sc2[-1,ay]/c],'r-.')
        if len(svs) > 3 :
            for i in range(0, svs[3].shape[0]):
                orbs.plot(svs[3][i,:,ax]/c,svs[3][i,:,ay]/c,'k-')
                orbs.plot([cdm[0,ax],svs[3][i,0,ax]/c],[cdm[0,ay],svs[3][i,0,ay]/c],'k--')
                orbs.plot([cdm[-1,ax],svs[3][i,-1,ax]/c],[cdm[-1,ay],svs[3][i,-1,ay]/c],'k-.')
        orbs.plot(cdm[:,ax],cdm[:,ay],'-+')
        orbs.set_xlabel('x={:} (ls)'.format(ax))
        orbs.set_ylabel('y={:} (ls)'.format(ay))
        #orbs.set_xticks(xt)
        #orbs.set_yticks(yt)
        orbs.set_title('Orbits')

        if show: 
            dessin.show()
        
        return dessin



    def Plot_orbits(self):
        '''
            figure = Plot_orbits()

            Returns plots of the orbits as well as plots of the errors on some of the integrals of motion.
        '''
        Msol = 1.988435e30         # Masse du Soleil (kg)
        c = 299792458.
        self.Compute_integrals_of_motion()
        svs = self.Get_interp_state_vectors(PSB = rotate_to_psb)
        sp = svs[0]
        sc1 = svs[1]
        sc2 = svs[2]
        cdm = self.Get_center_of_mass_positions()
        Mp, Mi, Mo = self.Get_masses()
        t = self.Get_tinterp()
        pE = self.Get_energies()
        pR = self.Get_center_of_mass_positions()
        pL = np.zeros( (pE.size, 3), dtype = np.float )
        pP = self.Get_center_of_mass_impulsions()#np.zeros( (pE.size, 3), dtype = np.float )

        dessin = figure(figsize=(16,8))
        orbs = dessin.add_subplot(331)
        if (Mo >0.) :
            xt = [sc2[:,0].min()/c, sc2[:,0].max()/c]
        else :
            xt = [sc1[:,0].min()/c, sc1[:,0].max()/c]
        yt = [sc1[:,2].min()/c, sc1[:,2].max()/c ]
        xt[0] += np.abs(xt[1] - xt[0])*0.1
        xt[1] -= np.abs(xt[1] - xt[0])*0.1
        yt[0] += np.abs(yt[1] - yt[0])*0.1
        yt[1] -= np.abs(yt[1] - yt[0])*0.1
        orbs.plot(sp[:,0]/c,sp[:,2]/c,'b-')
        orbs.plot([cdm[0,0],sp[0,0]/c],[cdm[0,2],sp[0,2]/c],'b--')
        orbs.plot(sc1[:,0]/c,sc1[:,2]/c,'g-')
        orbs.plot([cdm[0,0],sc1[0,0]/c],[cdm[0,2],sc1[0,2]/c],'g--')
        if (Mo >0.) :
            orbs.plot(sc2[:,0]/c,sc2[:,2]/c,'r-')
            orbs.plot([cdm[0,0],sc2[0,0]/c],[cdm[0,2],sc2[0,2]/c],'r--')
            orbs.plot([cdm[-1,0],sc2[-1,0]/c],[cdm[-1,2],sc2[-1,2]/c],'r-.')
        orbs.plot(cdm[:,0]/c,cdm[:,2]/c,'-+')
        orbs.set_xlabel('ls')
        orbs.set_ylabel('ls')
        orbs.set_xticks(xt)
        orbs.set_yticks(yt)
        orbs.set_title('Orbits from Earth')

        orbsp = dessin.add_subplot(332)
        orbsp.plot(sp[:,0]/c,sp[:,2]/c,'b-')
        orbsp.plot(sc1[:,0]/c,sc1[:,2]/c,'g-')
        #orbsp.plot(sc2[:,0]/c,sc2[:,2]/c,'-')
        orbsp.plot(cdm[:,0]/c,cdm[:,2]/c,'p-')
        xt = [ sc1[:,0].min()/c, sc1[:,0].max()/c ]
        yt = [ sc1[:,2].min()/c, sc1[:,2].max()/c ]
        xt[0] += np.abs(xt[1] - xt[0])*0.1
        xt[1] -= np.abs(xt[1] - xt[0])*0.1
        yt[0] += np.abs(yt[1] - yt[0])*0.1
        yt[1] -= np.abs(yt[1] - yt[0])*0.1
        orbsp.set_xticks(xt )
        orbsp.set_yticks(yt )
        orbsp.set_xlabel('ls')
        orbsp.set_ylabel('ls')
        orbsp.set_title('Inner orbits')

        nrj = dessin.add_subplot(333)
        plE = (pE-pE[0])/pE[0]
        xt = [t.min(), t.max()]
        yt = [plE.min(), plE.max()]
        nrj.plot(t, plE, '-')
        nrj.set_xticks(xt)
        nrj.set_yticks(yt)
        nrj.set_xlabel('Time since first toa (days)')
        nrj.set_ylabel('Relative error')
        nrj.set_title('Energy')

        desP = dessin.add_subplot(334)
        desP.plot(pP[:,0]/Msol,pP[:,2]/Msol,'-')
        xt = [pP[:,0].min()/Msol, pP[:,0].max()/Msol]
        yt = [pP[:,2].min()/Msol, pP[:,2].max()/Msol]
        xt[0] += np.abs(xt[1] - xt[0])*0.1
        xt[1] -= np.abs(xt[1] - xt[0])*0.1
        yt[0] += np.abs(yt[1] - yt[0])*0.1
        yt[1] -= np.abs(yt[1] - yt[0])*0.1
        desP.set_xlabel('m/s')
        desP.set_ylabel('m/s')
        desP.set_xticks(xt)
        desP.set_yticks(yt)
        desP.set_title('P/$M_\odot$')

        #desL = dessin.add_subplot(335)
        #plL = (pL[:,1]-pL[0,1])/max(1.,abs(pL[0,1]))
        #xt = [ t.min(), t.max() ]
        #yt = [plL.min(), plL.max()]
        #xt[0] += np.abs(xt[1] - xt[0])*0.1
        #xt[1] -= np.abs(xt[1] - xt[0])*0.1
        #yt[0] += np.abs(yt[1] - yt[0])*0.1
        #yt[1] -= np.abs(yt[1] - yt[0])*0.1
        #desL.plot(t,plL,'-')
        #desL.set_xticks(xt)
        #desL.set_yticks(yt)
        #desL.set_xlabel('Time since first toa (days)')
        #desL.set_ylabel('Relative error')
        #desL.set_title('L')

        desR = dessin.add_subplot(336)
        desR.plot((pR[:,0]-pR[0,0]),(pR[:,2]-pR[0,2]),'x')
        xt = [ (pR[:,0]-pR[0,0]).min(), (pR[:,0] - pR[0,0]).max() ]
        yt = [ (pR[:,2]-pR[0,2]).min(), (pR[:,2]-pR[0,2]).max() ]
        xt[0] += np.abs(xt[1] - xt[0])*0.1
        xt[1] -= np.abs(xt[1] - xt[0])*0.1
        yt[0] += np.abs(yt[1] - yt[0])*0.1
        yt[1] -= np.abs(yt[1] - yt[0])*0.1
        desR.set_xticks( xt )
        desR.set_yticks( yt)
        desR.set_xlabel('m')
        desR.set_ylabel('m')
        desR.set_title('Barycenter')

        dessin.subplots_adjust(hspace = 0.5, wspace = 0.5, bottom = 0., top = 0.95, left = 0.1, right = 0.95)


        return dessin


    def Get_DMX(self):
        '''
            Return DMX, DMXranges
        '''
        absmap = self.Get_absolute_parameter_map()
        params = self.Get_parameters()
        DMXranges = self.Get_DMX_ranges()
        maxpar = absmap[29] + DMXranges.size -1
        DMX = params[absmap[29]:maxpar]
        return DMX, DMXranges


    def Plot_DMX(self, start_sun =56415., end_sun=56445, figsize=(12,6)):
        # 15 deg fron the sun = 3rd of may to 2nd of june ie 2456415.500000 to 2456445.500000

        Pearth = 365.25

        DMX, DMXranges = self.Get_DMX()

        fig = figure(figsize=figsize)
        plt = fig.add_subplot(111)
        plt.plot(DMXranges[1:], DMX)

        # Fitting a line
        fitc = fitting.LinearLSQFitter()
        model = models.Linear1D(slope=0., intercept=0.)
        fitres = fitc(model, DMXranges[1:], DMX)
        plt.plot(DMXranges[1:], fitres(DMXranges[1:]), 'x-', color='green', label='{:.1e} + {:.1e}MJD'.format(fitres.intercept.value, fitres.slope.value))


        ni = (DMXranges[0] - start_sun) / Pearth
        if ni == np.floor(ni) :
            ni = np.int(np.floor(ni))
        else :
            ni =  np.int(np.floor(ni))+1

        nf = (DMXranges[-1] - end_sun) / Pearth
        nf =  np.int(np.floor(nf))
        for i in range(ni, nf+1):
            if i == ni:
                plt.axvline(x=start_sun + Pearth*i, linestyle='dashed', label='Sun proximity', color='red')
            else :
                plt.axvline(x=start_sun + Pearth*i, linestyle='dashed', color='red')
            plt.axvline(x=end_sun + Pearth*i, linestyle='dashed', color='red')


        plt.axhline(y=DMX.mean(), label='DMX mean={:.1e}'.format(DMX.mean()))
        plt.set_xlabel('End of DMX range (MJD)')
        plt.set_ylabel('DMX correction (pc/cm3)')
        plt.xaxis.label.set_fontsize(self.fontsizes['xx-large'])
        plt.yaxis.label.set_fontsize(self.fontsizes['xx-large'])
        plt.legend(loc=0, fontsize=self.fontsizes['xx-large'])
        plt.tick_params(axis='both', which='major', labelsize=self.fontsizes['xx-large'])
        plt.tick_params(axis='both', which='minor', labelsize=self.fontsizes['xx-large'])

        return fig


    def Get_outliers(self,max_residual, min_residual=None, first_toa=None, last_toa = None, mjd=False):
        '''
            return the list of indexes of outliers

            outliers are defined as toas with residual (as per Get_time_residuals()) min_residual > residual  or resdual > max_residual

            max_residual : in microsec (as per Get_time_residuals)
            min_residual : if None = - abs(max_residual)
            first_toa : start of  (barycentric) toa range where outliers are searched for. In days. If mjd = True, in mjds. If None start from first.
            last_toa : end of  (barycentric) toa range where outliers are searched for. In days. If mjd = True, in mjds. If None corresponds to the last toa.
            mjd : determines if first and last toa are given in Modified Julian Days (mjd=True) or in time since Get_timeshift() (mjd=False)

        '''
        toas = self.Get_toas()
        res  = self.Get_time_residuals()
        maxres = abs(max_residual)
        if min_residual is None :
            minres = - maxres
        else :
            minres = -abs(min_residual)
        # res -= res.mean()
        outindexes =[]
        timeshift = self.Get_timeshift()
        if first_toa is not None:
            firsttoa = toas.searchsorted(first_toa - mjd*timeshift)
        else:
            firsttoa = 0

        if last_toa is not None:
            lasttoa = toas.searchsorted(last_toa - mjd*timeshift)
        else:
            lasttoa = toas.size

        for i in range(firsttoa, lasttoa):
            if res[i] > maxres or res[i]< minres :
                outindexes.append(i)

        return outindexes
