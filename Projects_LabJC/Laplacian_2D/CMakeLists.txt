file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Laplacian2D_GID ${main_SRCS})

target_link_libraries(Laplacian2D_GID pz)
