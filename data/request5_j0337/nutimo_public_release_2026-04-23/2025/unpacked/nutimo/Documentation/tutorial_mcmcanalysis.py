#coding:utf8
#
# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later
################################################
# Small tutorial written by Guillaume Voisin, LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr ; astro.guillaume.voisin@gmail.com)
# October 2017
# 
# Load the necessary modules
import mcmanalysis as mc
import numpy as np
import Fittriple as ft

# Create a "Fittriple" object (containing the timing model) with the parameters used for the MCMC 
fit = ft.Fittriple("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/parfile-reference-isma11_newt_posransom_forMCMC", "0337-FDM.20160521-trimmed")

# Create a "mcmcresult" object to analyse the MCMC chain. Assume all initials parameters are "zero" so that one can see the deviation. 
# "parameter_object=fit" is not compulsory. Only without the parameter rescaling has to be done manually using the scales in the parameter file anf keyword "parameter_scale".
# Similarly the parameter map, giving the mapping between the mcmc parameters and the parameter indexes in the Fittriple object is not correctly set. 
mcres = mc.mcmcresult(resfile="mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25.res", 
                      parameter_object=fit, parameter_initials=np.zeros(29), symbols=mc.symbols)

# Alternative version without "parameter_object=fit", but with correct parameter map and scales
paramscales = np.array([  8.19200000e-07,   2.44140625e-09,   1.60000000e-06,
         1.63840000e-01,   3.20000000e-06,   1.02400000e-01,
         1.95312500e-07,   4.76837158e-10,   3.90625000e-08,
         4.19430400e+01,   3.12500000e-08,   5.24288000e+01,
         1.56250000e-06,   3.12500000e-07,   4.76837158e-11,
         2.98023224e-11,   3.81469727e-10,   1.60000000e-01,
         3.20000000e-05,   4.88281250e-04,   1.00000000e-07,
         4.00000000e-09,   1.60000000e-08,   2.00000000e+02,
         2.50000000e-01,   1.00000000e+00,   1.00000000e-09,
         1.25000000e-01,   9.76562500e-04])
parammap = np.array([21, 22, 23, 24, 25, 26,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
       11, 12, 13, 15, 17, 18, 20])

mcres = mc.mcmcresult(resfile="mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25.res",
                      parameter_map=parammap,parameter_scales=paramscales, parameter_initials=np.zeros(29), symbols=mc.symbols)

# Note : mc.symbols contains the latex symbols of the parameters in the order of the "Fittriple" objects. 

# The Log(posterior proba) is in the variable chi2. To see the number of steps/iterations just look at :
mcres.chi2.size

# When the chain has converged the autocorrelation time is small. If it is not enough then an error is printed
mcres.Autocorrelation_time(start=300000)
# So one tries to start the chain at a step further :
mcres.Autocorrelation_time(start=600000)
# This one is correct, it gives small autocorrelation times with respect to the number of steps still taken into account. 

# Show the mean parameters ( with a shift that puts all the initials parameters to zero) computed on the end of the chain
mcres.Mean_parameters(start=600000)

# Make a corner/triangle plot
fig = mcres.Corner_plot(start=600000)
# The plot is too large to be shown directly in python, it is better to save it and look with a thrid-party viewer
fig.savefig("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/cornerplot-start600000.png")

# Plot the end of the chain and save it
fig = mcres.Chain_plot(start=600000)
fig.savefig("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/chainplot-start600000.png")


# Create a new mcmcresult object by selectring only the iterations of the chain for which parameter 20 (SEP delta) is between -1e-6 and 1e-6
mcres_sep0 = mcres.Parameter_filter(20, -1.e-6, 1.e-6)
mcres_sep0.chi2.size
# Look for the possibily converged area
mcres_sep0.Autocorrelation_time(50000)
mcres_sep0.Autocorrelation_time(74000)
# In this case there are only ~25000 steps possibily converged. It is small given the high number of dimensions. But we continue and make the plots :
fig = mcres_sep0.Corner_plot(start=74000)
fig.savefig("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/cornerplot-sep0-start74000.png")
fig = mcres_sep0.Chain_plot(start=74000)
fig.savefig("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/chainplot-sep0-start74000.png")

# Do the same for the second peak of the SEP Delta parameter distribution, around 3.5e-6
mcres_sep1 = mcres.Parameter_filter(20, 1.e-6, 1.e-5)
mcres_sep1.chi2.size
mcres_sep1.Autocorrelation_time(500000)
fig = mcres_sep1.Corner_plot(start=500000)
fig.savefig("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/cornerplot-sep1-start500000.png")
fig = mcres_sep1.Chain_plot(start=500000)
fig.savefig("mcmc-results/mcmc-reference-isma11_newt_posransom_forMCMC-432x250000d25/chainplot-sep1-start500000.png")
