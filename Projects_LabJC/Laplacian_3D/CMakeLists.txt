file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Laplacian3D ${main_SRCS})

target_link_libraries(Laplacian3D pz)
