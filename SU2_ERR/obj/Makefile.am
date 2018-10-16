################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD
# \author M. Colonno, T. Economon, F. Palacios
# \version 6.1.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS =

if BUILD_NORMAL
bin_PROGRAMS += ../bin/SU2_ERR
endif

if BUILD_DIRECTDIFF
bin_PROGRAMS += ../bin/SU2_ERR_DIRECTDIFF
endif

if BUILD_REVERSE
bin_PROGRAMS += ../bin/SU2_ERR_AD
endif

noinst_LIBRARIES= 

libSU2Core_sources = \
  ../../SU2_CFD/include/definition_structure.hpp \
  ../../SU2_CFD/include/fluid_model.hpp \
  ../../SU2_CFD/include/fluid_model.inl \
  ../../SU2_CFD/include/integration_structure.hpp \
  ../../SU2_CFD/include/integration_structure.inl \
  ../../SU2_CFD/include/driver_structure.hpp \
  ../../SU2_CFD/include/iteration_structure.hpp \
  ../../SU2_CFD/include/numerics_structure.hpp \
  ../../SU2_CFD/include/numerics_structure.inl \
  ../../SU2_CFD/include/output_structure.hpp \
  ../../SU2_CFD/include/output_structure.inl \
  ../../SU2_CFD/include/solver_structure.hpp \
  ../../SU2_CFD/include/solver_structure.inl \
  ../../SU2_CFD/include/solver_boom.hpp \
  ../../SU2_CFD/include/transport_model.hpp \
  ../../SU2_CFD/include/transport_model.inl \
  ../../SU2_CFD/include/variable_structure.hpp \
  ../../SU2_CFD/include/variable_structure.inl \
  ../../SU2_CFD/include/transfer_structure.hpp \
  ../../SU2_CFD/include/transfer_structure.inl \
  ../../SU2_CFD/src/definition_structure.cpp \
  ../../SU2_CFD/src/fluid_model.cpp \
  ../../SU2_CFD/src/fluid_model_pig.cpp \
  ../../SU2_CFD/src/fluid_model_pvdw.cpp \
  ../../SU2_CFD/src/fluid_model_ppr.cpp \
  ../../SU2_CFD/src/fluid_model_inc.cpp \
  ../../SU2_CFD/src/integration_structure.cpp \
  ../../SU2_CFD/src/integration_time.cpp \
  ../../SU2_CFD/src/driver_structure.cpp \
  ../../SU2_CFD/src/iteration_structure.cpp \
  ../../SU2_CFD/src/numerics_adjoint_mean.cpp \
  ../../SU2_CFD/src/numerics_adjoint_turbulent.cpp \
  ../../SU2_CFD/src/numerics_direct_poisson.cpp \
  ../../SU2_CFD/src/numerics_direct_heat.cpp \
  ../../SU2_CFD/src/numerics_direct_mean.cpp \
  ../../SU2_CFD/src/numerics_direct_mean_inc.cpp \
  ../../SU2_CFD/src/numerics_direct_transition.cpp \
  ../../SU2_CFD/src/numerics_direct_turbulent.cpp \
  ../../SU2_CFD/src/numerics_direct_wave.cpp \
  ../../SU2_CFD/src/numerics_direct_elasticity_nonlinear.cpp \
  ../../SU2_CFD/src/numerics_direct_elasticity_linear.cpp \
  ../../SU2_CFD/src/numerics_direct_elasticity.cpp \
  ../../SU2_CFD/src/numerics_structure.cpp \
  ../../SU2_CFD/src/numerics_template.cpp \
  ../../SU2_CFD/src/output_cgns.cpp \
  ../../SU2_CFD/src/output_structure.cpp \
  ../../SU2_CFD/src/output_physics.cpp \
  ../../SU2_CFD/src/output_tecplot.cpp \
  ../../SU2_CFD/src/output_fieldview.cpp \
  ../../SU2_CFD/src/output_su2.cpp \
  ../../SU2_CFD/src/output_paraview.cpp \
  ../../SU2_CFD/src/python_wrapper_structure.cpp \
  ../../SU2_CFD/src/solver_adjoint_mean.cpp \
  ../../SU2_CFD/src/solver_adjoint_turbulent.cpp \
  ../../SU2_CFD/src/solver_adjoint_discrete.cpp \
  ../../SU2_CFD/src/solver_adjoint_elasticity.cpp \
  ../../SU2_CFD/src/solver_direct_poisson.cpp \
  ../../SU2_CFD/src/solver_direct_heat.cpp \
  ../../SU2_CFD/src/solver_direct_mean.cpp \
  ../../SU2_CFD/src/solver_direct_mean_inc.cpp \
  ../../SU2_CFD/src/solver_direct_transition.cpp \
  ../../SU2_CFD/src/solver_direct_turbulent.cpp \
  ../../SU2_CFD/src/solver_direct_wave.cpp \
  ../../SU2_CFD/src/solver_direct_elasticity.cpp \
  ../../SU2_CFD/src/solver_structure.cpp \
  ../../SU2_CFD/src/solver_template.cpp \
  ../../SU2_CFD/src/solver_boom.cpp \
  ../../SU2_CFD/src/transfer_physics.cpp \
  ../../SU2_CFD/src/transfer_structure.cpp \
  ../../SU2_CFD/src/transport_model.cpp \
  ../../SU2_CFD/src/variable_adjoint_mean.cpp \
  ../../SU2_CFD/src/variable_adjoint_turbulent.cpp \
  ../../SU2_CFD/src/variable_adjoint_discrete.cpp \
  ../../SU2_CFD/src/variable_adjoint_elasticity.cpp \
  ../../SU2_CFD/src/variable_direct_poisson.cpp \
  ../../SU2_CFD/src/variable_direct_heat.cpp \
  ../../SU2_CFD/src/variable_direct_mean.cpp \
  ../../SU2_CFD/src/variable_direct_mean_inc.cpp \
  ../../SU2_CFD/src/variable_direct_transition.cpp \
  ../../SU2_CFD/src/variable_direct_turbulent.cpp \
  ../../SU2_CFD/src/variable_direct_wave.cpp \
  ../../SU2_CFD/src/variable_direct_elasticity.cpp \
  ../../SU2_CFD/src/variable_structure.cpp \
  ../../SU2_CFD/src/variable_template.cpp

su2_err_sources = \
  ../include/SU2_ERR.hpp \
  ../src/SU2_ERR.cpp

libSU2Core_cxx_flags = -fPIC
libSU2Core_libadd = 

su2_err_cxx_flags =
su2_err_ldadd =

# always link to built dependencies from ./externals
su2_err_cxx_flags += @su2_externals_INCLUDES@
su2_err_ldadd += @su2_externals_LIBS@
libSU2Core_cxx_flags += @su2_externals_INCLUDES@

# if BUILD_MUTATIONPP
su2_err_cxx_flags += @MUTATIONPP_CXX@
su2_err_ldadd += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
su2_err_cxx_flags += @JSONCPP_CXX@
su2_err_ldadd += @JSONCPP_LD@
# endif

# if BUILD_HDF5
su2_err_cxx_flags += @HDF5_CXX@
su2_err_ldadd += @HDF5_LD@
# endif

# if BUILD_SZIP
su2_err_cxx_flags += @SZIP_CXX@
su2_err_ldadd += @SZIP_LD@
# endif

# if BUILD_ZLIB
su2_err_cxx_flags += @ZLIB_CXX@
su2_err_ldadd += @ZLIB_LD@
# endif

if BUILD_NORMAL
___bin_SU2_ERR_SOURCES = ${su2_err_sources}
___bin_SU2_ERR_CXXFLAGS = ${su2_err_cxx_flags}
___bin_SU2_ERR_LDADD =  ../../SU2_CFD/obj/libSU2Core.a ../../Common/lib/libSU2.a ${su2_err_ldadd}
endif

if BUILD_REVERSE
___bin_SU2_ERR_AD_SOURCES = $(su2_err_sources) ${libSU2Core_sources}
___bin_SU2_ERR_AD_CXXFLAGS = @REVERSE_CXX@ ${su2_err_cxx_flags}
___bin_SU2_ERR_AD_LDADD = @REVERSE_LIBS@ ../../Common/lib/libSU2_AD.a ${su2_err_ldadd}
endif
