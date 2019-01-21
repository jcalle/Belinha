# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

#To include Error Estimation Projects
if(NOT BUILD_TUTORIAL)
	add_subdirectory(ErrorEstimation)
	add_subdirectory(ErrorEstimationHDiv)
endif()

