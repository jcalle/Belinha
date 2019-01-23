# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

#To include Juca Costa Lab Projects
if(NOT BUILD_TUTORIAL AND NOT BUILD_PROJECTS)
	add_subdirectory(2017Arequipa)
	add_subdirectory(ErrorEstimation)
	add_subdirectory(Projects_LabJC)
endif()

#Choi Sung-Boing