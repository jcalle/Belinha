# @file Belinha/Adaptivity/CMakeList.txt  -- File to read for CMake

#disable excessive warning messages by Microsoft VC compiler
set (CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_LIST_DIR}/CompilerOptionsMSVC.cmake")

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE  INTERNAL FORCE )

IF (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
ELSE (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
	if ( NOT APPLE AND NOT CMAKE_BUILD_TYPE )
		set (CMAKE_BUILD_TYPE Release CACHE STRING "Options: Release Debug None MinSizeRel RelWithDebInfo" FORCE) 
	endif()
ENDIF()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## inserting main definitions

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include usr/lib /usr/lib64 /usr/local/lib ./../externallibs/pzlib ./../externallibs/pzlib/include ./../externallibs/pzlib/include/Analysis ./../externallibs/pzlib/include/Geom ./../externallibs/pzlib/include/Mesh ./../externallibs/pzlib/lib)


# neopz library
find_path(NEOPZ_INCLUDE NAMES pz_config.h PATHS ${SEARCH_DIRS} )
find_library(NEOPZ_LIB NAMES libpz.a PATHS ${SEARCH_DIRS})
if(NEOPZ_INCLUDE-NOTFOUND)
	set (NEOPZ_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(NEOPZ_LIB-NOTFOUND)
	set (NEOPZ_LIB "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${NEOPZ_INCLUDE})


#target_link_libraries(pz ${SEARCH_DIRS})
#target_include_directories(pz PUBLIC ${SEARCH_DIRS})

add_subdirectory(ErrorEstimators)
add_subdirectory(StrategyRef)
add_subdirectory(WaveletCoef)

#Creating groups for organized display of files into PZ modules (subfolders)
source_group(ErrorEstimators FILES ${estimators_STAT_SRCS})
source_group(StrategyRefinements FILES ${refinement_STAT_SRCS})
source_group(WaveletCoefficient FILES ${wavelet_STAT_SRCS})


