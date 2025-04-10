find_package(OpenMP)

add_library(RSAalgolib AlgoPacking.cxx Path.cxx MultiDArrayObject.cxx Ellipse.cxx BasicGeometricOperations.cxx Area.cxx)
target_include_directories(RSAalgolib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(RSAalgolib PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(RSAalgolib PRIVATE -fPIC)
install(TARGETS RSAalgolib DESTINATION lib)

add_executable(RSAalgo Main_pour_test.cxx)
target_include_directories(RSAalgo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(RSAalgo PUBLIC OpenMP::OpenMP_CXX RSAalgolib)
install(TARGETS RSAalgo DESTINATION bin)

