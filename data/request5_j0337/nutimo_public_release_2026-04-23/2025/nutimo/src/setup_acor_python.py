# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.extension import Extension
#
# Run  python setup_cppFittriple.py build_ext --inplace après avoir généré les librairie .a en faisant un "make lib" 
# Attention  ! make lib must be run with a -fpic version of the tempo2 libraries
#
# To create libacorcpp.so : 
#   g++ -c -fpic acc.cpp
#   g++ -shared -o libacorcpp.so acc.o
# Then :
#   python setup_acor_python.py build_ext --inplace 

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension(name="acor_python", 
                  sources=["acor_python.pyx",],
                  include_dirs = ["./"],
                  library_dirs = ["./"],
                  libraries =["acorcpp"],#, "tempo2", "sofa"],   # refers to "libexternlib.so"
                  language="c++",                   # remove this if C and not C++
                  extra_compile_args=["-O2"],#["-fopenmp", "-O2"],#,"-I."]
                  extra_link_args=[]
             )
        ]
)     

