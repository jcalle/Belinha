# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

#To include Juca Costa Lab Projects
if(NOT BUILD_TUTORIAL)
	add_subdirectory(Laplacian2D_GID)
endif()

