# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Guillaume VOISIN , LUTh, Observatoire de Paris - PSL (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)  
#coding:utf8
from distutils.core import setup
#from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.extension import Extension
import os
#
# Run  python setup_cppFittriple.py build_ext --inplace après avoir généré les librairie .a en faisant un "make lib" 
# Attention  ! make lib must be run with a -fpic version of the tempo2 libraries
#

installpath = os.getenv('nutimo')

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension(name="python_Fittriple_interface", 
                  sources=["python_Fittriple_interface.pyx",],
                  include_dirs = ["./", installpath+ "/third_party/boost_1_55_0/", installpath+"/third_party/tempo2/include/"],
                  library_dirs = ["./"], #"libstatictempo2/with_fpic/"],
                  libraries =["Fittriplecpp"],#, "tempo2", "sofa"],   # refers to "libexternlib.so"
                  language="c++",                   # remove this if C and not C++
                  extra_compile_args=["-O2"],#["-fopenmp", "-O2"],#,"-I."]
                  extra_link_args=[]#"-DSOME_DEFINE_OPT", 
                   #                "-L./some/extra/dependency/dir/"]
             )
        ]
)     

