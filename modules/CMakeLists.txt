#voro++
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
add_subdirectory(voro-plus-plus)

# testing
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()

# FFTW
#include("cmake/modules/fftw3.cmake")
include("cmake/modules/gcc.cmake")

# AlgoPacking implicitly includes compiler options
add_subdirectory(AlgoPacking)
add_subdirectory(merope_core)


