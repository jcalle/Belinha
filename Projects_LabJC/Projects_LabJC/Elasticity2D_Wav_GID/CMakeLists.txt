file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Elasticity2D_Wav_GID ${main_SRCS})

target_link_libraries(Elasticity2D_Wav_GID pz)
