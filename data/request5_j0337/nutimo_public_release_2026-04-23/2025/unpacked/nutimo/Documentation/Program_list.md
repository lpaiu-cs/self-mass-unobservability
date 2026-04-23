<!--
SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>

SPDX-License-Identifier: GPL-3.0-or-later
-->

# List of all the programs and scripts coming with Nutimo.

 For all programs, a simple call `Program.exe` without arguments normally shows a quick help. 


## `Make_fake_timfile.exe parfile original_timfile noisestddev faketimfile`
Make a mock/fake timfile, with SATs based on original timfile but shifted to exactly correspond to the model with parameters from parfile. A Gaussian noise is added on top of it.
- noisestddev is the standard deviation of the added Gaussian noise
Output : 
    - tempo2 timfile saved in faketimfile
    - a text file with detailed intermediate quantities in faketimfile-details

## `mpirun Cord-2points.exe parfile timfile npoints outfile`
compute a cord between two points with relative shifts in cord_point1.txt and cord_point2.txt.
The cord has npoints
The lnposterior of each point is stored in outfile.


## `Check_accuracy.exe parfile timfile factor tolerance`
Perform various accuracy checks.

## `Minuit_fit.exe. original_parfile timfile postfit`
Fit timfile data starting from parameters in original_parfile and saving it in "parfile-postfit"


## `mcmc2parfile.exe parfile timfile mcmcresultfile outputparfile`
mcmcresultfile is a file containing results from a mcmc simulation as done with MCMC_pai_pt.exe using parfile and timfile.
The last element of the mcmc chain is converted into a parfile saved to oututparfile


## `Convert_parameter_sets.exe parfile pset [MCMC_chain_to_convert]`
Convert "parfile" to the parameter set given by pset. If present, also convert all the lines of "MCMC_chain_to_convert" which must be a result file produced by MCMC_pai_pt.exe.
The output files bear identical names with "-converted" suffixes


## `MCMC_pai_pt.exe`
See the dedicated "MCMC" file in this folder. Does the MCMC inference for Nutimo. 


## `New_parfile_from_MCMC.exe parfile timfile MCMC_result_file new_parfile key`
Create a new parfile from a MCMC chain. The 'key' decides what statistics to apply on the chain to derive a single set. Parameter scales are set equal to chain standard deviations. 
    parfile : parfile associated with the MCMC chain. 
    timfile : needed for only initialisation. Does not matter which one.
    MCMC_result_file : chain as obtained using MCMC_pai_pt.exe.
    new_parfile : filename of the new parfile.
    key : what parameters to extract from the MCMC chain result. Can be 'mean' or 'max' (max lnposterior)

## `Convert_MCMC_chain.exe current_parfile timfile MCMC_result_file new_parfile new_MCMC_result_file`
Create a new MCMC chain file relative to the parameters in "new_parfile". 
If a parameter is fitted in "new_parfile" but not in "current_parfile" then this parameters is added with randomly generated values in "new_MCMC_result_file" according to the "initial_distribution" method of the lnposterior object (see MCMC_lnposterior_functions.h). If the central value or the scale change between the two parfiles, the values in the MCMC file is converted accordingly. 
Note: a timfile is needed although it does not matter which one


## `integrate.exe parfile start end [-o outputprefix][-n integer][-i][-l integer][-r resumefile]`
Integrate equations of motion using parameters from "parfile" and save results to a file.
The output format is r0_x r0_y r0_z v0_x v0_y v0_z ... rn_x .. vn_z for n+1 bodies. 
    parfile : nutimo type parfile 
    start : start date in days <=0 relative to the reference time in parfile 
    end : end date in days >=0 relative to the reference time in parfile
    -o outputprefix : prefix of output file(s). default=trajectories
    -no : No Output of trajectories to file. Does prevent output from -i, -l or other options.
    -n integer : Number of times to record on an evenly sampled grid. default=interpolation parameters from parfile
    -i : compute and save integrals of motion in a separate file
    -l integer : compute mean longitude of the corresponding subsystem in hierarchical order and save resuld to a separate file with two columns (time in days | mean long in radians). Mean longitude = argument of asc. node (Oman) +  arg. of periastron (omperi) + mean anomaly (m). The option can be called several times creating one separate file each time. (0<= integer < nbody -1)
    -r resumefile : initial state vector is overwritten by the one found in 'resumefile', e.g. a file containing the last state vector of a previous trajetory file.

## `integrate_convert.exe parfile trajfile [-o outputprefix][-i][-l integer]`
Convert state vectors in trajfile (same format as output by integrate.exe) into other quantities, in particular orbital elements. By default, orbital elements are calculated w.r.t the invariant plane of the system.
    parfile : nutimo type parfile 
    trajfile : file of the form 't(days) r_1 v_1 ... r_n v_n' where 'n' is the number of bodies and 'r' and 'v' are position and velocities. These coordinates will be converted to orbital elements and other derived quantities according to options. If -m is given then masses are taken from a separate file, otherwise a single set of masses is used thoughout and taken from the parfile. 
    -o outputprefix : prefix of output file(s). default=trajectories
    -i : compute and save integrals of motion in a separate file
    -l integer : compute mean longitude of the corresponding subsystem in hierarchical order and save resuld to a separate file with two columns (time in days | mean long in radians). Mean longitude = argument of asc. node (Oman) +  arg. of periastron (omperi) + mean anomaly (m). The option can be called several times creating one separate file each time. (0<= integer < nbody -1)
    -m massfilename : each line of this file contain the set of masses (in solar mass) to use with the corresponding line of trajfile


## `integrate.py`
Load and plot results from integrate.exe

## `Fittriple.py`
Python interface to Nutimo, designed with the triple system J0337+1715 in mind. See also "Small_tutorial_Fittriple_python"


## `mcmcanalysis.py`
Python module to analyse the results of MCMC runs, such as those produced by "MCMC_pai_pt.exe".
