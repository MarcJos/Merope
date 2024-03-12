add_subdirectory(pybind11)

# clean name
pybind11_add_module(merope Interface.cxx ../../AlgoPacking/Interface_python/Interface_AlgoPacking.cxx)
target_link_libraries(merope PUBLIC GeomVER Grid VTKinout VoroInterface VoroToMesh GeneOrientations OptiLaguerre)
if(DEFINED MEROPE_USE_OPENMP_FOR_FFT)
    target_compile_definitions(merope PUBLIC USE_OPENMP_FOR_FFT=1)
    MESSAGE(STATUS "-> bad interaction with python lambdas in the gaussian function")
endif(DEFINED MEROPE_USE_OPENMP_FOR_FFT)


target_compile_definitions(merope PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_compile_options(merope PRIVATE -fPIC)
### for MKL
install(TARGETS merope DESTINATION lib)

# for retrocompatibility
pybind11_add_module(Merope Interface.cxx ../../AlgoPacking/Interface_python/Interface_AlgoPacking.cxx)
target_link_libraries(Merope PUBLIC GeomVER Grid VTKinout VoroInterface VoroToMesh GeneOrientations OptiLaguerre)
target_compile_definitions(Merope PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
if(DEFINED MEROPE_USE_OPENMP_FOR_FFT)
    target_compile_definitions(Merope PUBLIC USE_OPENMP_FOR_FFT=1)
    MESSAGE(STATUS "-> bad interaction with python lambdas in the gaussian function")
endif(DEFINED MEROPE_USE_OPENMP_FOR_FFT)

target_compile_options(Merope PRIVATE -fPIC)
### for MKL
install(TARGETS Merope DESTINATION lib)

