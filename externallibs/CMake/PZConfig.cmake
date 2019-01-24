# - Config file for the PZ package
# It defines the following variables
#  PZ_INCLUDE_DIRS - include directories for using PZ
#  PZ_LIBRARIES    - PZ library to link against

# In the future, if project executables are to be installed, their names should also be defined here:
#  <exec01_name>_EXECUTABLE   
#  <exec02_name>_EXECUTABLE  
#  ...
#
# so that they can be used below
 
# Compute paths
get_filename_component(PZ_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(PZ_INCLUDE_DIRS "${PZ_CMAKE_DIR}/../pzlib/include;/../pzlib/include/Python;${PZ_CMAKE_DIR}/../pzlib/include/Util;${PZ_CMAKE_DIR}/../pzlib/include/PerfUtil;${PZ_CMAKE_DIR}/../pzlib/include/Common;${PZ_CMAKE_DIR}/../pzlib/include/Save;${PZ_CMAKE_DIR}/../pzlib/include/Matrix;${PZ_CMAKE_DIR}/../pzlib/include/Topology;${PZ_CMAKE_DIR}/../pzlib/include/Geom;${PZ_CMAKE_DIR}/../pzlib/include/SpecialMaps;${PZ_CMAKE_DIR}/../pzlib/include/Refine;${PZ_CMAKE_DIR}/../pzlib/include/Shape;${PZ_CMAKE_DIR}/../pzlib/include/Material;${PZ_CMAKE_DIR}/../pzlib/include/Material/REAL;${PZ_CMAKE_DIR}/../pzlib/include/Material/REAL/Plasticity;${PZ_CMAKE_DIR}/../pzlib/include/Material/Complex;${PZ_CMAKE_DIR}/../pzlib/include/Multigrid;${PZ_CMAKE_DIR}/../pzlib/include/Mesh;${PZ_CMAKE_DIR}/../pzlib/include/StrMatrix;${PZ_CMAKE_DIR}/../pzlib/include/Integral;${PZ_CMAKE_DIR}/../pzlib/include/Frontal;${PZ_CMAKE_DIR}/../pzlib/include/Pre;${PZ_CMAKE_DIR}/../pzlib/include/Post;${PZ_CMAKE_DIR}/../pzlib/include/Random;${PZ_CMAKE_DIR}/../pzlib/include/Optimization;${PZ_CMAKE_DIR}/../pzlib/include/Analysis;${PZ_CMAKE_DIR}/../pzlib/include/SubStruct;${PZ_CMAKE_DIR}/../pzlib/include/LinearSolvers;${PZ_CMAKE_DIR}/../pzlib/include/External;${PZ_CMAKE_DIR}/../pzlib/include/External/sloan;${PZ_CMAKE_DIR}/../pzlib/include/Publications;${PZ_CMAKE_DIR}/../pzlib/include/Exception;${PZ_CMAKE_DIR}/../pzlib/include/External/FAD;${PZ_CMAKE_DIR}/../pzlib/include/External/FAD/Fad;${PZ_CMAKE_DIR}/../pzlib/include/External/FAD/TinyFad;${PZ_CMAKE_DIR}/../pzlib/include/External/FAD/TinyFadET")
set(PZ_BRANCH "master")
set(PZ_REVISION "0f240136e")
set(PZ_REVISION_DATE "Thu Jan 24 09:42:47 2019")
 
# Our library dependencies (contains definitions for IMPORTED targets)
if(NOT TARGET pz AND NOT PZ_BINARY_DIR)
  include("${PZ_CMAKE_DIR}/PZTargets.cmake")
endif()
 
# These are IMPORTED targets created by PZTargets.cmake
set(PZ_LIBRARIES pz)

set(PZ_BUILD_PLASTICITY_MATERIALS OFF)
set(PZ_REAL_TYPE double)
set(PZ_STATE_TYPE double)
set(PZ_USING_BOOST OFF)
set(PZ_USING_FAD ON)
set(PZ_USING_LOG4CXX OFF)
set(PZ_USING_METIS OFF)
set(PZ_USING_OPENSSL OFF)
set(PZ_USING_TBB OFF)
set(PZ_USING_OPENMP OFF)
set(PZ_USING_LIKWID OFF)
set(PZ_USING_LIBNUMA OFF)
set(PZ_USING_MATLAB_ENGINE OFF)
set(PZ_USING_LAPACK OFF)
set(PZ_USING_BLAS OFF)
set(PZ_USING_PAPI OFF)
set(PZ_USING_HWLOC OFF)
set(PZ_BUILD_PYTHON_BINDING OFF)
set(PZ_PTHREAD_LIB "/usr/lib/libpthread.dylib")
set(PZ_PTHREAD_INCLUDE "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk/usr/include")
set(PZ_USING_MKL OFF)

if( EXISTS "${PZ_CMAKE_DIR}/PZConfig_Debug.cmake" )
	include("${PZ_CMAKE_DIR}/PZConfig_Debug.cmake")
endif()

#(for future decision!)
#set(<exec01_name>_EXECUTABLE <name> )   

