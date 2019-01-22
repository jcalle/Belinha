# @file CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

project (BELINHA)

#To include NeoPZ Repository
	add_subdirectory(neopz)

#To include Juca Costa Lab Projects
#if(NOT BUILD_TUTORIAL AND NOT BUILD_PROJECTS AND NOT BUILD_PLASTICITY_MATERIALS AND NOT BUILD_TESTING)
#	add_subdirectory(Projects_LabJC Projects_JC)
#endif()

#To include Adaptive library
	add_subdirectory(AdaptivityLib)

#To include Projects implemented in Juca Costa Laboratory
	add_subdirectory(Projects_LabJC)

