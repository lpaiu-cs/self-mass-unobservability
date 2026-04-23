<!--
SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>

SPDX-License-Identifier: GPL-3.0-or-later
-->

## To fit with Minuit :
* First set the number of threads of open mp to 1 : `export OMP_NUM_THREADS=1`
* The syntax is one of : 
    1. `./Minuit_fit.exe parfile timfile turnnumbers outparfile_suffix`
    2. `./Minuit_fit.exe parfile timfile`

In case 1) the number of turns of the pulsar corresponding to each TOA is given by the file "turnnumbers". This kind of file is automatically generated at every run (see below). 
     
In case 2) the number of turns of the pulsar corresponding to each TOA is computed assuming the timing model on the parfile is accurate to less than a spin period. 
        
The output parfile filename is the concatenation "parfile-" + "outparfile_suffix"


## The output files of `./Minuit_fit.exe` are :
- "turns.dat" : the number of turns of the pulsar corresponding to each TOA is computed assuming the timing model on the parfile is accurate to less than a spin period. 
- "parfile_avant_minuit" : parameter file before running Minuit. Contains the same parameters as "parfile". Useful to make a backup if you overwrite your initial parameter file "parfile". Also, this writes the most complete parfile possible, which avoids leaving implicit parameters implicitely determined by the program.
- "datafile_avant_minuit.dat" : contains three columns "toas", "toes" (Times of emission) and the "residuals" before Minuit is run.
- "timfile" + "-sorted" : timfile sorted by increasing time of arrival. Internally, "Fittriple" needs sorted times of arrival, which is not necessarily the case in tempo2 .tim files, in particular when the several frequency bands are present. 
- "Fake_Bats" + "timfile" : Contains the times of arrival predicted by the model with the parfile parameters before Minuit is run. 
- "parfile-" + "outparfile_suffix": the parfile at the end of fit with Minuit. 


