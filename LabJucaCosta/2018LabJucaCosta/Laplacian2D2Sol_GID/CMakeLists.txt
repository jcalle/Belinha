file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Laplacian2D2Sol ${main_SRCS})

target_link_libraries(Laplacian2D2Sol pz)
