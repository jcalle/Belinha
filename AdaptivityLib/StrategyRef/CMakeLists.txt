# @file CMakeList.txt  -- First file to read for CMake

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(refinements_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
