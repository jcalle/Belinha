# @file CMakeList.txt  -- First file to read for CMake

file(GLOB headers *.h)
file(GLOB sources *.cpp)

install(FILES ${headers} DESTINATION ${LABJC_INCLUDE}/Analysis)
set(wavelet_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
