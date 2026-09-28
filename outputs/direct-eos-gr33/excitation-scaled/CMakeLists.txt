# src/CMakeLists.txt
# Copyright (C) 2006-2022 Alan W. Irwin
#
#
#This file is part of FreeEOS.
#
#This file is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.
#
#This file is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Library General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with the file; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

# Build free_eos library.

set(supported_fp_precision 15 CACHE STRING "This floating-point decimal precision will be supported by the kind=fp_kind value that is selected in the Fortran source code")
if(supported_fp_precision LESS 15)
  message(FATAL_ERROR "supported_fp_precision = ${supported_fp_precision} must be >= 15")
endif(supported_fp_precision LESS 15)

set(supported_fp_exponent 300 CACHE STRING "This floating-point positive decimal exponent will be supported by the kind=fp_kind value that is selected in the Fortran source code")

if(supported_fp_exponent LESS 300)
  message(FATAL_ERROR "supported_fp_exponent = ${supported_fp_exponent} must be >= 300")
endif(supported_fp_exponent LESS 300)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/included_fortran_types.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/included_fortran_types.f90
  @ONLY
)

set(if_taint_allocated_real ".false." CACHE STRING "For gfortran only taint real arrays that are either allocated or equivalenced to work around a gfortran -finit limitation")
if(NOT (if_taint_allocated_real STREQUAL ".false." OR if_taint_allocated_real STREQUAL ".true."))
  message(FATAL_ERROR "if_taint_allocated_real = ${if_taint_allocated_real} must be either \".true.\" or \".false.\"")
endif(NOT (if_taint_allocated_real STREQUAL ".false." OR if_taint_allocated_real STREQUAL ".true."))

set(if_taint_allocated_integer ".false." CACHE STRING "For gfortran only taint integer arrays that are either allocated or equivalenced to work around a gfortran -finit limitation")
if(NOT (if_taint_allocated_integer STREQUAL ".false." OR if_taint_allocated_integer STREQUAL ".true."))
  message(FATAL_ERROR "if_taint_allocated_integer = ${if_taint_allocated_integer} must be either \".true.\" or \".false.\"")
endif(NOT (if_taint_allocated_integer STREQUAL ".false." OR if_taint_allocated_integer STREQUAL ".true."))

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/included_tainting_control.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/included_tainting_control.f90
  @ONLY
)

set(module_basenames
  # Primary module
  mod_free_eos
  # Other modules that are used by the primary module.
  mod_free_eos_constants
  mod_pi_fit
  mod_ionization_data
  mod_nuvar
  mod_aux_scale
  mod_diagnostics
  mod_excitation_block
  mod_helium1_data
  mod_statistical_weight_data
  mod_isotopic_mass_data
  mod_fermi_dirac_integrator
  mod_free_eos_types
  mod_bfgs
  mod_poly_sum
  mod_molecular_hydrogen
  mod_excitation
  mod_eos_calc
  mod_fermi_dirac
  mod_lapack
  mod_coulomb
  mod_exchange
  mod_pi
  mod_pl
  mod_free_eos_detailed
  mod_eos_bfgs
  mod_eos_jacobian
  mod_info_data
  mod_flow_data
  mod_free_eos_debug
  )

set(source_files)
set(module_files)
foreach(module_basename ${module_basenames})
  list(APPEND source_files ${module_basename}.f90)
  list(APPEND module_files ${CMAKE_CURRENT_BINARY_DIR}/${module_basename}.mod)
endforeach(module_basename ${module_basenames})

set(tll_list ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
#message(STATUS "DEBUG: tll_list = ${tll_list}")

set(LIB_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}" "${LAPACK_RPATH}")

configure_library_build(free_eos "" "${source_files}" "${tll_list}" "${LIB_INSTALL_RPATH}")

target_include_directories(${WRITEABLE_TARGET}free_eos
  # Everything that links to this library will need access to the module files.
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODDIR}>
  PRIVATE
  # The library itself depends on included source files from these directories.
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

set_target_properties(${WRITEABLE_TARGET}free_eos PROPERTIES Fortran_PREPROCESS ON)

# While user programmes should only use mod_free_eos
# some compilers (e.g. Intel Fortran, but not gfortran) require access
# to all the modules that mod_free_eos depends upon.
install(FILES ${module_files} DESTINATION ${CMAKE_INSTALL_MODDIR})

install(
FILES README.developers
DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/free_eos
)

set_target_properties(${WRITEABLE_TARGET}free_eos PROPERTIES OUTPUT_NAME free_eos_direct24_excitation_scaled)
