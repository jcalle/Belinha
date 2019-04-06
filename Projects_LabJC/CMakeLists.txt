# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

#To include Juca Costa Lab Projects
if(NOT BUILD_TUTORIAL AND NOT BUILD_PROJECTS)

	#Curso em Arequipa 2017
	#source_group(Aqp_2017 FILES ${aqp2017})
	add_subdirectory(2017Arequipa)

	#Jorge 2017 - Implementing error analysis
	#add_subdirectory(ABelinha)
	add_subdirectory(ErrorEstimation)

	add_subdirectory(Projects_LabJC)

	endif()


#Choi Sung-Boing