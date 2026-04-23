<!--
SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>

SPDX-License-Identifier: GPL-3.0-or-later
-->

#Directory Content: 

This folder contains the Nutimo's (Numerical Timing Model) source code. This software can analyse times of arrivals from pulsars in systems with multiple companions. The core is dedicated to deal with a triple system the orbits of which are numerically integrated at using first-post-newtonian-order equations of motion. It was initially developed to deal with the case of PSR J0337+1715 which is in a hierarchical triple stellar system with two white-dwarf companions. 

- `src/` : source code of the Nutimo software
- `resources/` : third-party files such as clock files or earth ephemeredis 
- `Documentation/` : some help
- `third_party/` : contains third party software if you decide to put it there (see `install_script.sh`)
- `install_script.sh` : hopefully should allow to compile and install Nutimo

# Contact 
Guillaume Voisin, LUX, Observatoire de Paris, PSL Research University  CNRS (guillaume.voisin@obspm.fr or  astro.guillaume.voisin@google.com )

# Licence
Copyright 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) under GPL-3.0-or-later

See  `Licenses/`

# References and how to cite 
If using this code please cite at least [Voisin+2024][4] and others below if relevant:

- Chapter 5 of [Voisin17][1] : detailed description of main implemented physical models. 
- [Voisin+2020][2] : further descprition of the physics, first application to PSR J0337+1715. Code release in [Zenodo2020][3].
- [Voisin+2024][4] : n-body implementation and application to a 4th object in the PSR J0337+1715 system.Code and data release in [Zenodo2024][5]


[1]: https://hal.archives-ouvertes.fr/tel-01677325 "G. Voisin's PhD thesis"
[2]: https://doi.org/10.1051/0004-6361/202038104 "Voisin et al. 2020"
[3]: https://zenodo.org/records/3778978 "Nutimo 2020"
[4]: https://doi.org/10.1051/0004-6361/202452100
[5]: https://doi.org/10.5281/zenodo.13899771

# Installation
The script "install_script.sh" should allow the installation of Nutimo.
The user should manually edit the relevant paths at the beginning of the script. 
The script is made to install the dependencies below as well. Tempo2 in particular requires some specific compilation options. Of course one can also skip these sections of the script provided as long as linking directories are appropriate. 

# Dependencies:
Tested on Linux Mint 18.3 64-bit and Linux Mint 20.1. with `gcc` version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.2), and `python` version 3.8.5.

- Boost library: www.boost.org. Version 1.55.0 was used in all papers. 
- Tempo2: https://bitbucket.org/psrsoft/tempo2/src/master/
    - Voisin+20 version: 2018.02.01
    - Voisin+24 version: downloaded on Git repo on 2021.04
- Minuit: http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/
Version for all papers: 5.34.14 
- `Acor`, Jonathan Goodman, sping 2009. It can be downloaded from J. Goodman's website https://math.nyu.edu/~goodman/software/acor/index.html. Note that `acor` is only used for diagnostics and Nutimo can run without it simply by commenting lines where `acor` is used. 


# Troubleshooting :
Tempo2's dependencies are of course required as well. Nutimo uses Tempo2 as an external library, and it is possible that libraries required by Tempo2 are not properly linked in Nutimo's makefile. If so, identify the missing library in error messages and append relevant linking information to the `LDFLAGS` and `LDFLAGSfpic`variable in `src/makefile-original`. For instance, the presence of `-llapack` in these variables is uniquely to accomodate tempo2's requirements. 
