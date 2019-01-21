# @file CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

project (BELINHA)

#specifying the path to neopz source code
set(BELSOURCEDIR ${PROJECT_SOURCE_DIR})

#To include NeoPZ Repository
	add_subdirectory(neopz)

