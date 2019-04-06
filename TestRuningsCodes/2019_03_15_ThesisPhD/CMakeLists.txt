# @file neopz/Projects/CursoPZ/ThesisPhD/CMakeLists.txt  -- Simple project to test bugs in three dimensions

#   add_definitions(-DTPZEQNARRAY)
#   add_definitions(-DPZENVIRONMENT)
   add_definitions(-D__VC__)

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(ThesisPhD ${main_SRCS})

target_link_libraries(ThesisPhD  pz)

