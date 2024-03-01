find_package(OpenMP)

add_executable(Optimization_test Main.cxx Tests.cxx)
install(TARGETS Optimization_test DESTINATION bin)

