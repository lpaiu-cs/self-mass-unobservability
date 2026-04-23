<!--
SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>

SPDX-License-Identifier: GPL-3.0-or-later
-->

# What ?
MCMC_pai_pt.exe provides a parallelized, affine-invariant MCMC for Nutimo with the possibility of using parallel tempering.

This program samples a distribution using the algorithm described in Goodman & Weare 2010 (GW10), Communitcations in Applied mathematics and Computational Sciences, Vol5, Issue 1
The parallelization scheme is from : Foreman-Mackey et al., 2013, Pulications of the Astronomical Society of the Pacific, Vol 125, Issue 925

Used in :

- Voisin et al. 2020
- Voisin et al. 2024

# Who ? 
Written by Guillaume Voisin 2017-2019 ,
JBCA, The University of Manchester, UK (guillaume.voisin@manchester.ac.uk)
LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)


# How to use it ? 

If called without arguments "MCMC_pai_pt.exe" will prompt a brief summary of the syntax and arguments.  

The syntax should be : 

`[mpirun -n number__cores] MCMC_pai_pt.exe parfile datafile fileout half_number_of_number_of_walkers_per_processor number_of_walker_set_moves`

Or : 

`[mpirun -n number__cores] MCMC_pai_pt.exe parfile datafile fileout half_number_of_number_of_walkers_per_processor number_of_walker_set_moves frequency_of_chain_extension`

Or : 

`[mpirun -n number__cores] MCMC_pai_pt.exe parfile datafile fileout half_number_of_number_of_walkers_per_processor number_of_walker_set_moves frequency_of_chain_extension number_of_proc_by_node`

Or : 

`[mpirun -n number__cores] MCMC_pai_pt.exe parfile datafile fileout half_number_of_number_of_walkers_per_processor number_of_walker_set_moves frequency_of_chain_extension number_of_proc_by_node stat_and_save_frequency [Optional keyword options]`

Note that only the last syntax accepts optional keyword options. 


## Compulsory arguments :

  parfile : Nutimo parfile

  datafile : Tempo2 timfile

  fileout : prefix of the name of the file containing the chain. The full filename is of the form fileout--temp=X.XX.res where X.XX is the temperature. If parallel tempering is used, then one file is created for each temperature.

  half_number_of_walkers_per_processor: each computing core computes 2 *  half_number_of_walkers_per_processor walkers

  number_of_walker_set_moves: total number of iterations

  frequency_of_chain_extension : number of iteration between two extension of the chain (i.e. recording of the state of chain) with the current walkers

  number_of_proc_per_node : number of cores on one node of a computer cluster. (if using a single node/compute/laptop just write the number of cores you want to use)

  stat_and_save_frequency : number of iterations between computing and displaying statistics, and saving the chain to file.


## Note on total number of walkers :
    The number of walkers of a run is 2 *  half_number_of_walkers_per_processor * number_of_cores where number_of_cores is 1 if run in monoprocessor mode and specified with the mpirun command otherwise. If running on several nodes, specifying number_of_proc_per_node allows a (much) faster initialization, but may lead to bugs if incorrect.

## Optional keyword options :

  turnfile <filename> : "filename" contains the integer number of turns (from an arbitrary reference) each ToA in "datafile" corresponds to. Essential if initial    parameters are very incaccurate. Useless otherwise. 

  previouschain <filename_prefix> : to restart from a previous run which file is named according to <filename_prefix>--temp=X.XX.res

  recordrefused : enables the recording of candidate walkers to which a jump has been refused. Useful for diagnostics

  aparam <value> : value of the parameter "a" of Goodman and Weare (2010) that characterises the proposal distribution. Default is 2.

  targettemp <float value> : Temperature of the run, lowest temperature if parallel tempering is enabled. Default is 1.

  tempswapfreq <int value> : Number of iterations between two attempts to swap walkers between temperatures

  ntemp <int value> : Number of parallel temperatures to run. Default is 1.

  tempfactor <float value> <float value> ... : list of ratios of Tn+1/Tn where Tn are the parallel temperatures. There should be ntemp-1. If there is only one then the same ratio is applied everywhere, otherwise there needs to be ntemp - 1 values given. Default ratio is sqrt(2) (see e.g. Earl and Deem 2005, DOI 10.1039/b509983h).

  nproc_per_temperature <int value> : For parallelisation of temperature. Gives the number of processors dedicated to one temperature. 

  --sieving <float value> : remove the walkers for which Best_log_posterior - log(posterior) > value. Removal happens every stat_and_save_frequency iterations. Help removing stuck walkers. To use with care. 

  --diagnostics : save some extra files for diagnostics

  --no_mean_calc: does not calculate chi2 at mean during diagnostics (i.e. every stat_and_save_frequency steps). This can be used to save computing time.

