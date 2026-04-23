#! /bin/bash
# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later
rootpath=$(pwd)/ # DO NOT CHANGE THIS

#####################################################################################################
# THIS SECTION IS THE ONLY ONE THAT A USER SHOULD HAVE TO CHANGE ! (otherwise, sorry for the bug!) 

# Write here the paths to the directories containing the source files of the third party softwares:
minuitsrc=$rootpath/third_party/Minuit2-5.34.14/ # Minuit sources
tempo2src=$rootpath/third_party/tempo2/          # Tempo2 sources.
boostsrc=$rootpath/third_party/boost_1_55_0/     # Boost sources.  
acorsrc=$rootpath/third_party/acor/              # Should contain acc.cpp and acor.h. Nothing else needed
resourcespath=$rootpath/ressources/

# Write here the path where you want to install executables and python scripts:
installdir=/path/to/installation/directory/ 

#####################################################################################################

export nutimo=$installdir

sourcepath=$rootpath/src


minuitinstall=$installdir/third_party/Minuit2-5.34.14/
tempo2install=$installdir/third_party/tempo2/
tempo2staticlib=$installdir/third_party/libstatictempo2/
tempo2staticfpiclib=$tempo2staticlib/with_fpic/
boostinstall=$installdir/third_party/boost_1_55_0/

mkdir $installdir
mkdir $installdir/third_party


### Create link to tempo2 sources. This is a quick and dirty trick to make sure Nutimo "make python" finds tempo2's headers. Another way to fix it is to change paths by hand in src/setup_cppFittriple.py.
mkdir $rootpath/third_party
ln -s tempo2src $rootpath/third_party/

# ### Minuit
echo "#######################################################################################################################"
echo
echo
echo                      INSTALLING MINUIT
echo
echo
echo "#######################################################################################################################"

cd $minuitsrc
./configure --prefix=$minuitinstall
autoreconf -i # # to comply with the installed version of autoconf
make
make install
make clean
#

# ### Tempo2
echo "#######################################################################################################################"
echo
echo
echo                      INSTALLING TEMPO2
echo
echo
echo "#######################################################################################################################"

mkdir $tempo2staticlib
mkdir $tempo2staticfpiclib
cd $tempo2src

autoreconf -f -i # to comply with the installed version of autoconf
./bootstrap
./configure --prefix=$tempo2install

sed -i 's/-DHAVE_OPENMP/-UHAVE_OPENMP/g' Makefile # disable openmp in tempo2

cp Makefile Makefile-without_fpic

# add the fpic option to the relevant lines in the Makefiles
sed -i 's/CFLAGS = -g -O2 -Wno-error/CFLAGS = -g -O2 -Wno-error -fpic/g' Makefile
sed -i 's/CPPFLAGS = -I/CPPFLAGS = -I -fpic/g' Makefile
sed -i 's/CXXFLAGS = -g -O2/CXXFLAGS = -g -O2 -fpic/g' Makefile

make
make install

cp $tempo2install/lib/*.a $tempo2staticfpiclib/

make clean

# now the same but without fpic

cp Makefile Makefile-with_fpic
cp Makefile-without_fpic Makefile

make
make install

cp $tempo2install/lib/*.a $tempo2staticlib

make clean

# Copy the runtime data of tempo2 to the installdir
cp -r $rootpath/third_party/tempo2/T2runtime/* $tempo2install/

# Update tempo2 with earth motion ephemeris and clock corrections located in ressources

cp -r resourcespathearth/* $tempo2install/earth/
cp -r resourcespathclock/* $tempo2install/clock/
cd $tempo2install/earth/
./update_eop.sh

# Link the tempo libraries in the src source folder to the src forlder for cython


#### Boost
echo "#######################################################################################################################"
echo
echo
echo                      INSTALLING BOOST
echo
echo
echo "#######################################################################################################################"

mkdir $boostinstall/
cp -r $boostsrc* $boostinstall/


### Nutimo itself
echo "#######################################################################################################################"
echo
echo
echo                      INSTALLING NUTIMO
echo
echo
echo "#######################################################################################################################"

cd $rootpath/src/

# Link acor source files 
ln -s $acorsrc/acor.h ./
ln -s $acorsrc/acc.cpp ./



# Backup the original makefile before modification
cp makefile-original makefile


# put the right paths into the makefile
sedargument="s|^TEMPO2 ="
sedargument=$sedargument'|TEMPO2 = '
sedargument=$sedargument$tempo2install
sedargument=$sedargument"/|"
echo $sedargument
sed -e "${sedargument}" makefile > makefile-inter

sedargument="s|^MINUIT ="
sedargument=$sedargument'|MINUIT = '
sedargument=$sedargument$minuitinstall
sedargument=$sedargument"/|"
sed -e "${sedargument}" makefile-inter >makefile

sedargument="s|^BOOST ="
sedargument=$sedargument'|BOOST = '
sedargument=$sedargument$boostinstall
sedargument=$sedargument"/|"
sed -e "${sedargument}" makefile > makefile-inter

sedargument="s|^LIBSTATICTEMPO2 ="
sedargument=$sedargument'|LIBSTATICTEMPO2 = '
sedargument=$sedargument$tempo2staticlib
sedargument=$sedargument"/|"
sed -e "${sedargument}" makefile-inter > makefile

sedargument="s|^LIBSTATICTEMPO2FPIC ="
sedargument=$sedargument'|LIBSTATICTEMPO2FPIC = '
sedargument=$sedargument$tempo2staticfpiclib
sedargument=$sedargument"/|"
sed -e "${sedargument}" makefile > makefile-inter

sedargument="s|^INSTALLDIR ="
sedargument=$sedargument'|INSTALLDIR = '
sedargument=$sedargument$installdir
sedargument=$sedargument"/|"
sed -e "${sedargument}" makefile-inter > makefile



#mv makefile-inter makefile

echo
echo
echo " ********Compile Minuit fitting program ******** "
echo
echo

make

echo
echo
echo " ********Compile MCMC_pai.exe for mcmc simulations ******** "
echo
echo

make mcmc

echo
echo
echo " ********Compile Python extensions ******** "
echo
echo

export PYTHONPATH=$PYTHONPATH:"./"
make python

echo " ******** Make install ******** "

make install



echo
echo
echo
echo
echo
echo " To run without linking problems you will probablyy need the followings :"
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$installdir" as otherwise python_Fittriple_interface.so does not see libFittriplecpp.so"
echo "export PYTHONPATH=$PYTHONPATH:$installdir"
echo "export TEMPO2="$installdir"/third_party/tempo2"
echo "ln -s $installdir/*.exe ~/.local/bin/"
