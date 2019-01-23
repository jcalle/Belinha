# @file CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

project (BELINHA)

set(BELSOURCEDIR ${PROJECT_SOURCE_DIR})

#include directories to neopz
include_directories(. ./externallibs/pzlib/include/Mesh ./externallibs/pzlib/include/Analysis ./externallibs/pzlib/include/Save ./externallibs/pzlib/include/Common ./externallibs/pzlib/include/Frontal ./externallibs/pzlib/include/Geom ./externallibs/pzlib/include/Integral ./externallibs/pzlib/include/LinearSolvers ./externallibs/pzlib/include/Material ./externallibs/pzlib/include/Matrix ./externallibs/pzlib/include/Multigrid ./externallibs/pzlib/include/Post ./externallibs/pzlib/include/Pre ./externallibs/pzlib/include/Refine ./externallibs/pzlib/include/Refine/RefPatterns ./externallibs/pzlib/include/Shape ./externallibs/pzlib/include/StrMatrix ./externallibs/pzlib/include/SubStruct ./externallibs/pzlib/include/Topology ./externallibs/pzlib/include/Util ./externallibs/pzlib/include/External ./externallibs/pzlib/include/)

link_directories( ./externallibs/pzlib/lib)

#To include Adaptive library
	add_subdirectory(AdaptivityLib)

#To include Projects implemented in Juca Costa Laboratory
	add_subdirectory(Projects_LabJC)


